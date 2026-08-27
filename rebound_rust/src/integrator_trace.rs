//! integrator_trace.rs — the TRACE hybrid (almost) time-reversible
//! integrator (from integrator_trace.c/h; Lu, Hernandez & Rein 2024,
//! after Hernandez & Dehnen 2023). WHFast for the long-term evolution,
//! switching to BS or IAS15 for close encounters, with a reversible
//! pre/post timestep switching check.
//!
//! Ownership notes (same conventions as integrator_mercurius.rs):
//! - `r->map` aliases `trace->encounter_map` during the encounter BS
//!   loop; here the Vec is moved into `r.map` and moved back.
//! - While the BS sub-integration runs, the trace state is stored in
//!   `r.integrator` so the custom gravity routine, the TRACE nbody
//!   derivatives and the particle hooks can reach it.
//! - The C switching functions read `r->integrator.state`; the Rust fn
//!   pointer types receive the state explicitly as a parameter.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Tiger Lu, Hanno Rein and contributors. See crate root.

use crate::collision::reb_collision_search;
use crate::integrator_bs::{
    reb_integrator_bs_nbody_derivatives, reb_integrator_bs_state, reb_integrator_bs_step_odes,
    reb_integrator_bs_update_particles, reb_ode, reb_ode_create, reb_ode_free,
};
use crate::integrator_ias15::{reb_integrator_ias15_state, reb_integrator_ias15_step_state};
use crate::integrator_whfast::reb_integrator_whfast_kepler_solver;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;

/// integrator_trace.c `#define MAX(a, b)` — the exact ternary.
fn MAX(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// Switching function for close encounters between non-central bodies
/// (C: `int (*S)(r, i, j)`; the state is explicit here).
pub type reb_trace_S_fn =
    fn(r: &reb_simulation, trace: &reb_integrator_trace_state, i: usize, j: usize) -> i32;
/// Switching function for close encounters with the central body
/// (C: `int (*S_peri)(r, j)`; the state is explicit here).
pub type reb_trace_S_peri_fn =
    fn(r: &reb_simulation, trace: &reb_integrator_trace_state, j: usize) -> i32;

/// integrator_trace.h `REB_INTEGRATOR_TRACE_PERIMODE` enum.
pub const REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS: i32 = 0;
pub const REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS: i32 = 1;
pub const REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15: i32 = 2;

/// integrator_trace.h internal mode enum.
pub const REB_INTEGRATOR_TRACE_MODE_INTERACTION: u32 = 0;
pub const REB_INTEGRATOR_TRACE_MODE_KEPLER: u32 = 1;
pub const REB_INTEGRATOR_TRACE_MODE_FULL: u32 = 3;

/// integrator_trace.h `struct reb_integrator_trace_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_trace_state {
    /// Switching function (None: use the default).
    pub S: Option<reb_trace_S_fn>,
    /// Pericenter switching function (None: use the default).
    pub S_peri: Option<reb_trace_S_peri_fn>,
    /// How TRACE integrates close approaches with the central star.
    pub peri_mode: i32,
    /// Critical switchover distance in units of the modified Hill radius.
    pub r_crit_hill: f64,
    /// Pericenter approach criterion (Pham, Rein & Spiegel 2024).
    pub peri_crit_eta: f64,
    // Internal use
    pub mode: u32,
    /// Number of particles currently having an encounter.
    pub encounter_N: usize,
    /// Number of active particles currently having an encounter.
    pub encounter_N_active: usize,
    /// C `N_allocated`.
    pub N_allocated: usize,
    /// 0 if any encounters are between two massive bodies. 1 if
    /// encounters only involve test particles.
    pub tponly_encounter: u32,
    /// Coordinates before the entire step.
    pub particles_backup: Vec<reb_particle>,
    /// Coordinates before the Kepler step.
    pub particles_backup_kepler: Vec<reb_particle>,
    /// Backup around additional_forces evaluation.
    pub particles_backup_additional_forces: Vec<reb_particle>,
    /// Map to represent which particles are integrated with BS. Moved
    /// into `r.map` during the encounter BS loop.
    pub encounter_map: Vec<usize>,
    /// Encounter map from after the pre-timestep check.
    pub encounter_map_backup: Vec<usize>,
    /// Centre of mass during the timestep.
    pub com_pos: reb_vec3d,
    pub com_vel: reb_vec3d,
    /// Tracking K_ij for the entire timestep (N*N row-major).
    pub current_Ks: Vec<i32>,
    /// Tracking C for the entire timestep.
    pub current_C: u32,
    /// Force accept for irreversible steps: collisions and adding particles.
    pub force_accept: u32,
}

impl Default for reb_integrator_trace_state {
    /// integrator_trace.c `reb_integrator_trace_create`.
    fn default() -> Self {
        reb_integrator_trace_state {
            S: None,
            S_peri: None,
            peri_mode: REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS,
            r_crit_hill: 3.,
            peri_crit_eta: 1.0,
            mode: REB_INTEGRATOR_TRACE_MODE_INTERACTION,
            encounter_N: 0,
            encounter_N_active: 0,
            N_allocated: 0,
            tponly_encounter: 0,
            particles_backup: Vec::new(),
            particles_backup_kepler: Vec::new(),
            particles_backup_additional_forces: Vec::new(),
            encounter_map: Vec::new(),
            encounter_map_backup: Vec::new(),
            com_pos: reb_vec3d::default(),
            com_vel: reb_vec3d::default(),
            current_Ks: Vec::new(),
            current_C: 0,
            force_accept: 0,
        }
    }
}

/// integrator_trace.c `reb_integrator_trace_switch_default` — returns 1
/// for a close encounter between i and j, 0 otherwise.
pub fn reb_integrator_trace_switch_default(
    r: &reb_simulation,
    trace: &reb_integrator_trace_state,
    i: usize,
    j: usize,
) -> i32 {
    let h2 = r.dt / 2.;

    let dxi = r.particles[i].x;
    let dyi = r.particles[i].y;
    let dzi = r.particles[i].z;

    let dxj = r.particles[j].x;
    let dyj = r.particles[j].y;
    let dzj = r.particles[j].z;

    let dx = dxi - dxj;
    let dy = dyi - dyj;
    let dz = dzi - dzj;
    let rp = dx * dx + dy * dy + dz * dz;

    let mut dcriti6 = 0.0;
    let mut dcritj6 = 0.0;

    let m0 = r.particles[0].m;

    // Check central body for physical radius ONLY
    if i == 0 && r.particles[i].r != 0. {
        let rs = r.particles[0].r;
        dcriti6 = rs * rs * rs * rs * rs * rs;
    } else if r.particles[i].m != 0. {
        let di2 = dxi * dxi + dyi * dyi + dzi * dzi;
        let mr = r.particles[i].m / (3. * m0);
        dcriti6 = di2 * di2 * di2 * mr * mr;
    }

    if r.particles[j].m != 0. {
        let dj2 = dxj * dxj + dyj * dyj + dzj * dzj;
        let mr = r.particles[j].m / (3. * m0);
        dcritj6 = dj2 * dj2 * dj2 * mr * mr;
    }

    let r_crit_hill2 = trace.r_crit_hill * trace.r_crit_hill;
    let dcritmax6 = r_crit_hill2 * r_crit_hill2 * r_crit_hill2 * MAX(dcriti6, dcritj6);

    if rp * rp * rp < dcritmax6 {
        return 1;
    }

    let dvx = r.particles[i].vx - r.particles[j].vx;
    let dvy = r.particles[i].vy - r.particles[j].vy;
    let dvz = r.particles[i].vz - r.particles[j].vz;
    let v2 = dvx * dvx + dvy * dvy + dvz * dvz;

    let qv = dx * dvx + dy * dvy + dz * dvz;
    let d: i32;

    if qv == 0.0 {
        // Small
        // minimum is at present, which is already checked for
        return 0;
    } else if qv < 0. {
        d = 1;
    } else {
        d = -1;
    }

    let dmin2: f64;
    let tmin = -(d as f64) * qv / v2;
    if tmin < h2 {
        // minimum is in the window
        dmin2 = rp - qv * qv / v2;
    } else {
        dmin2 = rp + 2. * (d as f64) * qv * h2 + v2 * h2 * h2;
    }

    (dmin2 * dmin2 * dmin2 < dcritmax6) as i32
}

/// integrator_trace.c `reb_integrator_trace_switch_peri_default` —
/// following Pham et al (2024).
pub fn reb_integrator_trace_switch_peri_default(
    r: &reb_simulation,
    trace: &reb_integrator_trace_state,
    j: usize,
) -> i32 {
    let GM = r.G * r.particles[0].m; // Not sure if this is the right mass to use.

    let x = r.particles[j].x;
    let y = r.particles[j].y;
    let z = r.particles[j].z;
    let d2 = x * x + y * y + z * z;
    let d = d2.sqrt();

    // first derivative
    let dx = r.particles[j].vx;
    let dy = r.particles[j].vy;
    let dz = r.particles[j].vz;

    // second derivative
    let prefact2 = -GM / (d2 * d);
    let ddx = prefact2 * x;
    let ddy = prefact2 * y;
    let ddz = prefact2 * z;
    // need sqrt for this one...
    let dd = (ddx * ddx + ddy * ddy + ddz * ddz).sqrt();

    // third derivative
    let prefact3 = GM / (d2 * d2 * d);
    let dddx = prefact3 * (-dx * (y * y + z * z) + 2. * x * x * dx + 3. * x * (y * dy + z * dz));
    let dddy = prefact3 * (-dy * (x * x + z * z) + 2. * y * y * dy + 3. * y * (x * dx + z * dz));
    let dddz = prefact3 * (-dz * (x * x + y * y) + 2. * z * z * dz + 3. * z * (x * dx + y * dy));

    let ddd2 = dddx * dddx + dddy * dddy + dddz * dddz;

    // fourth derivative
    let prefact4 = GM / (d2 * d2 * d2 * d);
    let ddddx = prefact4
        * (d2 * (-ddx * (y * y + z * z)
            + 2. * x * x * ddx
            + dx * (y * dy + z * dz)
            + x * (4. * dx * dx + 3. * (y * ddy + dy * dy + z * ddz + dz * dz)))
            - 5. * (x * dx + y * dy + z * dz)
                * (-dx * (y * y + z * z) + 2. * x * x * dx + 3. * x * (y * dy + z * dz)));
    let ddddy = prefact4
        * (d2 * (-ddy * (x * x + z * z)
            + 2. * y * y * ddy
            + dy * (x * dx + z * dz)
            + y * (4. * dy * dy + 3. * (x * ddx + dx * dx + z * ddz + dz * dz)))
            - 5. * (y * dy + x * dx + z * dz)
                * (-dy * (x * x + z * z) + 2. * y * y * dy + 3. * y * (x * dx + z * dz)));
    let ddddz = prefact4
        * (d2 * (-ddz * (y * y + x * x)
            + 2. * z * z * ddz
            + dz * (y * dy + x * dx)
            + z * (4. * dz * dz + 3. * (y * ddy + dy * dy + x * ddx + dx * dx)))
            - 5. * (z * dz + y * dy + x * dx)
                * (-dz * (y * y + x * x) + 2. * z * z * dz + 3. * z * (y * dy + x * dx)));
    let dddd = (ddddx * ddddx + ddddy * ddddy + ddddz * ddddz).sqrt();

    let tau_prs2 = 2. * dd * dd / (ddd2 + dd * dddd); // Eq 16
    let dt_prs2 = trace.peri_crit_eta * trace.peri_crit_eta * tau_prs2;

    if r.dt * r.dt > dt_prs2 {
        1
    } else {
        0
    }
}

/// integrator_trace.c `reb_integrator_trace_switch_peri_none` — no
/// pericenter flags.
pub fn reb_integrator_trace_switch_peri_none(
    _r: &reb_simulation,
    _trace: &reb_integrator_trace_state,
    _j: usize,
) -> i32 {
    0
}

/// integrator_trace.c `reb_integrator_trace_inertial_to_dh`
/// (state-explicit).
pub fn reb_integrator_trace_inertial_to_dh(
    r: &mut reb_simulation,
    trace: &mut reb_integrator_trace_state,
) {
    let particles = &mut r.particles;
    let mut com_pos = reb_vec3d::default();
    let mut com_vel = reb_vec3d::default();
    let mut mtot = 0.;
    let N_active = if r.N_active == usize::MAX || r.testparticle_type == 1 {
        r.N
    } else {
        r.N_active
    };
    let N = r.N;
    for i in 0..N_active {
        let m = particles[i].m;
        com_pos.x += m * particles[i].x;
        com_pos.y += m * particles[i].y;
        com_pos.z += m * particles[i].z;
        com_vel.x += m * particles[i].vx;
        com_vel.y += m * particles[i].vy;
        com_vel.z += m * particles[i].vz;
        mtot += m;
    }
    com_pos.x /= mtot;
    com_pos.y /= mtot;
    com_pos.z /= mtot;
    com_vel.x /= mtot;
    com_vel.y /= mtot;
    com_vel.z /= mtot;
    // Particle 0 is also changed to allow for easy collision detection
    let p0 = particles[0];
    for i in 0..N {
        particles[i].x -= p0.x;
        particles[i].y -= p0.y;
        particles[i].z -= p0.z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    trace.com_pos = com_pos;
    trace.com_vel = com_vel;
}

/// integrator_trace.c `reb_integrator_trace_dh_to_inertial`
/// (state-explicit; only reads the state).
pub fn reb_integrator_trace_dh_to_inertial(
    r: &mut reb_simulation,
    trace: &reb_integrator_trace_state,
) {
    let particles = &mut r.particles;
    let mut temp = reb_particle::default();
    let N = r.N;
    let N_active = if r.N_active == usize::MAX || r.testparticle_type == 1 {
        r.N
    } else {
        r.N_active
    };
    for i in 1..N_active {
        let m = particles[i].m;
        temp.x += m * particles[i].x;
        temp.y += m * particles[i].y;
        temp.z += m * particles[i].z;
        temp.vx += m * particles[i].vx;
        temp.vy += m * particles[i].vy;
        temp.vz += m * particles[i].vz;
        temp.m += m;
    }
    temp.m += particles[0].m;
    temp.x /= temp.m;
    temp.y /= temp.m;
    temp.z /= temp.m;
    temp.vx /= particles[0].m;
    temp.vy /= particles[0].m;
    temp.vz /= particles[0].m;
    // Use com to calculate central object's position.
    // This ignores previous values stored in particles[0].
    // Should not matter unless collisions occurred.
    particles[0].x = trace.com_pos.x - temp.x;
    particles[0].y = trace.com_pos.y - temp.y;
    particles[0].z = trace.com_pos.z - temp.z;

    for i in 1..N {
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += trace.com_vel.x;
        particles[i].vy += trace.com_vel.y;
        particles[i].vz += trace.com_vel.z;
    }
    particles[0].vx = trace.com_vel.x - temp.vx;
    particles[0].vy = trace.com_vel.y - temp.vy;
    particles[0].vz = trace.com_vel.z - temp.vz;
}

/// integrator_trace.c static
/// `reb_integrator_trace_calculate_acceleration_mode_interaction`
/// (state-explicit).
fn reb_integrator_trace_calculate_acceleration_mode_interaction(
    r: &mut reb_simulation,
    trace: &mut reb_integrator_trace_state,
) {
    let N = r.N;
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let _testparticle_type = r.testparticle_type;
    // (OPENMP branch and the `reb_sigint` signal check of the C omitted:
    // the reference Windows build has neither.)
    for i in 0..N {
        r.particles[i].ax = 0.;
        r.particles[i].ay = 0.;
        r.particles[i].az = 0.;
    }
    for i in 2..N_active {
        for j in 1..i {
            if trace.current_Ks[j * N + i] != 0 {
                continue;
            }
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let prefact = G / (_r * _r * _r);
            let prefactj = -prefact * r.particles[j].m;
            let prefacti = prefact * r.particles[i].m;
            r.particles[i].ax += prefactj * dx;
            r.particles[i].ay += prefactj * dy;
            r.particles[i].az += prefactj * dz;
            r.particles[j].ax += prefacti * dx;
            r.particles[j].ay += prefacti * dy;
            r.particles[j].az += prefacti * dz;
        }
    }
    let startitestp = std::cmp::max(N_active, 2);
    for i in startitestp..N {
        for j in 1..N_active {
            if trace.current_Ks[j * N + i] != 0 {
                continue;
            }
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let prefact = G / (_r * _r * _r);
            let prefactj = -prefact * r.particles[j].m;
            r.particles[i].ax += prefactj * dx;
            r.particles[i].ay += prefactj * dy;
            r.particles[i].az += prefactj * dz;
            if _testparticle_type != 0 {
                let prefacti = prefact * r.particles[i].m;
                r.particles[j].ax += prefacti * dx;
                r.particles[j].ay += prefacti * dy;
                r.particles[j].az += prefacti * dz;
            }
        }
    }

    // Handle Additional forces
    if r.additional_forces.is_some() {
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        // We should NOT do this in FULL mode, already in inertial frame
        if r.N > trace.particles_backup_additional_forces.len() {
            trace
                .particles_backup_additional_forces
                .resize(r.N, reb_particle::default());
        }
        trace.particles_backup_additional_forces[..r.N].copy_from_slice(&r.particles[..r.N]);
        reb_integrator_trace_dh_to_inertial(r, trace);
        (r.additional_forces.unwrap())(r);
        let particles = &mut r.particles;
        let backup = &trace.particles_backup_additional_forces;
        for i in 0..r.N {
            particles[i].x = backup[i].x;
            particles[i].y = backup[i].y;
            particles[i].z = backup[i].z;
            particles[i].vx = backup[i].vx;
            particles[i].vy = backup[i].vy;
            particles[i].vz = backup[i].vz;
        }
    }
}

/// Fetch entry `i` of the encounter map (`r.map` while the BS loop runs,
/// otherwise the state's own Vec) — see integrator_mercurius.rs.
fn encounter_map_get(
    r: &reb_simulation,
    trace: &reb_integrator_trace_state,
    i: usize,
) -> usize {
    match &r.map {
        Some(m) => m[i],
        None => trace.encounter_map[i],
    }
}

/// integrator_trace.c static
/// `reb_integrator_trace_calculate_acceleration_mode_kepler`
/// (state-explicit).
fn reb_integrator_trace_calculate_acceleration_mode_kepler_state(
    r: &mut reb_simulation,
    trace: &reb_integrator_trace_state,
) {
    let N = r.N;
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let _testparticle_type = r.testparticle_type;
    let m0 = r.particles[0].m;
    let encounter_N = trace.encounter_N;
    let encounter_N_active = trace.encounter_N_active;
    // (OPENMP branch of the C omitted.)
    r.particles[0].ax = 0.; // map[0] is always 0
    r.particles[0].ay = 0.;
    r.particles[0].az = 0.;

    // Acceleration due to star
    for i in 1..encounter_N {
        let mi = encounter_map_get(r, trace, i);
        let x = r.particles[mi].x;
        let y = r.particles[mi].y;
        let z = r.particles[mi].z;
        let _r = (x * x + y * y + z * z + softening2).sqrt();
        let prefact = -G * m0 / (_r * _r * _r);
        r.particles[mi].ax = prefact * x;
        r.particles[mi].ay = prefact * y;
        r.particles[mi].az = prefact * z;
    }

    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Interactions between active-active
    if encounter_N_active > 2 {
        // if two or less, no active-active planets
        for i in 2..encounter_N_active {
            let mi = encounter_map_get(r, trace, i);
            for j in 1..i {
                let mj = encounter_map_get(r, trace, j);
                if trace.current_Ks[mj * N + mi] == 0 {
                    continue;
                }
                let dx = r.particles[mi].x - r.particles[mj].x;
                let dy = r.particles[mi].y - r.particles[mj].y;
                let dz = r.particles[mi].z - r.particles[mj].z;
                let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
                let prefact = G / (_r * _r * _r);
                let prefactj = -prefact * r.particles[mj].m;
                let prefacti = prefact * r.particles[mi].m;

                r.particles[mi].ax += prefactj * dx;
                r.particles[mi].ay += prefactj * dy;
                r.particles[mi].az += prefactj * dz;
                r.particles[mj].ax += prefacti * dx;
                r.particles[mj].ay += prefacti * dy;
                r.particles[mj].az += prefacti * dz;
            }
        }
    }

    // Interactions between active-testparticle
    let startitestp = std::cmp::max(encounter_N_active, 2);
    for i in startitestp..encounter_N {
        let mi = encounter_map_get(r, trace, i);
        for j in 1..encounter_N_active {
            let mj = encounter_map_get(r, trace, j);
            if trace.current_Ks[mj * N + mi] == 0 {
                continue;
            }
            let dx = r.particles[mi].x - r.particles[mj].x;
            let dy = r.particles[mi].y - r.particles[mj].y;
            let dz = r.particles[mi].z - r.particles[mj].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let prefact = G / (_r * _r * _r);
            let prefactj = -prefact * r.particles[mj].m;
            r.particles[mi].ax += prefactj * dx;
            r.particles[mi].ay += prefactj * dy;
            r.particles[mi].az += prefactj * dz;

            if _testparticle_type != 0 {
                let prefacti = prefact * r.particles[mi].m;
                r.particles[mj].ax += prefacti * dx;
                r.particles[mj].ay += prefacti * dy;
                r.particles[mj].az += prefacti * dz;
            }
        }
    }
}

/// Fetch-version installed as `r.gravity_custom` (the C reads
/// `r->integrator.state` inside the same function).
pub fn reb_integrator_trace_calculate_acceleration_mode_kepler(r: &mut reb_simulation) {
    let trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_trace_calculate_acceleration_mode_kepler_state(r, &trace);
    r.integrator = reb_integrator_state::trace(trace);
}

/// integrator_trace.c `reb_integrator_trace_interaction_step`
/// (state-explicit).
pub fn reb_integrator_trace_interaction_step(
    r: &mut reb_simulation,
    trace: &mut reb_integrator_trace_state,
    dt: f64,
) {
    let N = r.N;
    trace.mode = REB_INTEGRATOR_TRACE_MODE_INTERACTION;
    reb_integrator_trace_calculate_acceleration_mode_interaction(r, trace);
    let particles = &mut r.particles;
    for i in 1..N {
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
    }
}

/// integrator_trace.c `reb_integrator_trace_jump_step` (state-explicit).
pub fn reb_integrator_trace_jump_step(
    r: &mut reb_simulation,
    trace: &reb_integrator_trace_state,
    dt: f64,
) {
    let current_C = trace.current_C;
    if current_C != 0 {
        return; // No jump step for pericenter approaches
    }

    let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };

    // If TP type 1, use r->N. Else, use N_active.
    let N = if r.testparticle_type == 0 { N_active } else { r.N };

    let mut px = 0.;
    let mut py = 0.;
    let mut pz = 0.;
    for i in 1..N {
        px += r.particles[i].vx * r.particles[i].m; // in dh
        py += r.particles[i].vy * r.particles[i].m;
        pz += r.particles[i].vz * r.particles[i].m;
    }
    px *= dt / r.particles[0].m;
    py *= dt / r.particles[0].m;
    pz *= dt / r.particles[0].m;

    let N_all = r.N;
    let particles = &mut r.particles;
    for i in 1..N_all {
        particles[i].x += px;
        particles[i].y += py;
        particles[i].z += pz;
    }
}

/// integrator_trace.c `reb_integrator_trace_com_step` (state-explicit).
pub fn reb_integrator_trace_com_step(trace: &mut reb_integrator_trace_state, dt: f64) {
    trace.com_pos.x += dt * trace.com_vel.x;
    trace.com_pos.y += dt * trace.com_vel.y;
    trace.com_pos.z += dt * trace.com_vel.z;
}

/// integrator_trace.c `reb_integrator_trace_whfast_step` — Kepler
/// evolution of all particles around the central mass (the C passes
/// r=NULL to the solver: no warnings, no variational particles).
pub fn reb_integrator_trace_whfast_step(r: &mut reb_simulation, dt: f64) {
    let N = r.N;
    let mu = r.G * r.particles[0].m;
    let mut particles = std::mem::take(&mut r.particles);
    let mut no_var: [reb_particle; 0] = [];
    for i in 1..N {
        reb_integrator_whfast_kepler_solver(None, &mut particles, &mut no_var, i, mu, dt);
    }
    r.particles = particles;
}

/// integrator_trace.c `reb_integrator_trace_update_particles`
/// (state-explicit).
fn reb_integrator_trace_update_particles_state(
    r: &mut reb_simulation,
    trace: &reb_integrator_trace_state,
    y: &[f64],
) {
    let N = trace.encounter_N;
    for i in 0..N {
        let mi = encounter_map_get(r, trace, i);
        let p = &mut r.particles[mi];
        p.x = y[i * 6];
        p.y = y[i * 6 + 1];
        p.z = y[i * 6 + 2];
        p.vx = y[i * 6 + 3];
        p.vy = y[i * 6 + 4];
        p.vz = y[i * 6 + 5];
    }
}

/// integrator_trace.c `reb_integrator_trace_update_particles` — fetches
/// the state from `r.integrator` like the C does.
pub fn reb_integrator_trace_update_particles(r: &mut reb_simulation, y: &[f64]) {
    let trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_trace_update_particles_state(r, &trace, y);
    r.integrator = reb_integrator_state::trace(trace);
}

/// integrator_trace.c `reb_integrator_trace_nbody_derivatives` — the
/// derivatives function of the temporary encounter ODE. Fetches the
/// trace state from `r.integrator` (the C uses `ode->r` +
/// `r->integrator.state`).
pub fn reb_integrator_trace_nbody_derivatives(
    r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    y: &[f64],
    _t: f64,
) {
    let trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    // TRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_trace_update_particles_state(r, &trace, y);
    reb_integrator_trace_calculate_acceleration_mode_kepler_state(r, &trace);

    let mut px = 0.;
    let mut py = 0.;
    let mut pz = 0.;
    let N = trace.encounter_N;

    if r.map.is_none() && trace.encounter_map.is_empty() {
        reb_simulation_error(r, "Cannot access TRACE map from BS.");
        r.integrator = reb_integrator_state::trace(trace);
        return;
    }

    // Kepler Step
    // This is only for pericenter approach
    if trace.current_C != 0 {
        for i in 1..r.N {
            // all particles
            px += r.particles[i].vx * r.particles[i].m; // in dh
            py += r.particles[i].vy * r.particles[i].m;
            pz += r.particles[i].vz * r.particles[i].m;
        }
        px /= r.particles[0].m;
        py /= r.particles[0].m;
        pz /= r.particles[0].m;
    }
    yDot[0] = 0.0;
    yDot[1] = 0.0;
    yDot[2] = 0.0;
    yDot[3] = 0.0;
    yDot[4] = 0.0;
    yDot[5] = 0.0;

    for i in 1..N {
        let mi = encounter_map_get(r, &trace, i);
        let p = r.particles[mi];
        yDot[i * 6] = p.vx + px; // Already checked for current_L
        yDot[i * 6 + 1] = p.vy + py;
        yDot[i * 6 + 2] = p.vz + pz;
        yDot[i * 6 + 3] = p.ax;
        yDot[i * 6 + 4] = p.ay;
        yDot[i * 6 + 5] = p.az;
    }
    r.integrator = reb_integrator_state::trace(trace);
}

/// Peek helpers for values needed while the state is stashed in
/// `r.integrator` during the BS loop.
fn peek_encounter_N(r: &reb_simulation) -> usize {
    match &r.integrator {
        reb_integrator_state::trace(t) => t.encounter_N,
        _ => 0,
    }
}

/// integrator_trace.c `reb_integrator_trace_bs_step` — BS integration
/// of the encounter subsystem. Takes/returns the state by value; while
/// the BS loop runs the state lives in `r.integrator`.
pub fn reb_integrator_trace_bs_step(
    r: &mut reb_simulation,
    mut trace: reb_integrator_trace_state,
    mut dt: f64,
) -> reb_integrator_trace_state {
    if trace.encounter_N < 2 {
        // No close encounters, skip
        return trace;
    }

    let mut i_enc: usize = 0;
    let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };
    trace.encounter_N_active = 0;
    for i in 0..r.N {
        if trace.encounter_map[i] != 0 {
            let tmp = r.particles[i]; // Copy for potential use for tponly_encounter
            r.particles[i] = trace.particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            trace.encounter_map[i_enc] = i;
            i_enc += 1;
            if i < N_active {
                trace.encounter_N_active += 1;
                if trace.tponly_encounter != 0 {
                    trace.particles_backup_kepler[i] = tmp; // Make copy of particles after the kepler step.
                                                            // used to restore the massive objects' states in the case
                                                            // of only massless test-particle encounters
                }
            }
        }
    }

    trace.mode = REB_INTEGRATOR_TRACE_MODE_KEPLER;
    r.map = Some(std::mem::take(&mut trace.encounter_map)); // for collision search
    r.N_map = trace.encounter_N;
    r.gravity = REB_GRAVITY::CUSTOM;
    r.gravity_custom = Some(reb_integrator_trace_calculate_acceleration_mode_kepler);

    let peri_mode = trace.peri_mode;
    let current_C = trace.current_C;

    // Store the state where the derivatives, gravity routine and
    // particle hooks find it.
    r.integrator = reb_integrator_state::trace(trace);

    // Only Partial BS uses this step
    if peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS || current_C == 0 {
        // run
        let old_dt = r.dt;
        let old_t = r.t;
        let t_needed = r.t + dt;
        let mut bs = reb_integrator_bs_state::default(); // reb_integrator_bs.create()

        // Temporarily remove all odes for BS step
        let odes_backup = std::mem::take(&mut r.odes);

        // Temporarily add new nbody ode for BS step
        let mut nbody_ode: Option<usize> = None;

        // TODO: Support backwards integrations
        while r.t < t_needed && (dt / old_dt).abs() > 1e-14 && r.status <= 0 {
            let encounter_N = peek_encounter_N(r);
            let stale = match nbody_ode {
                Some(id) => match r.odes.iter().find(|o| o.id == id) {
                    Some(o) => o.length != encounter_N * 3 * 2,
                    None => true,
                },
                None => true,
            };
            if stale {
                // (re)create the ODE
                if let Some(id) = nbody_ode {
                    reb_ode_free(r, id);
                }
                let id = reb_ode_create(r, encounter_N * 3 * 2);
                if let Some(ode) = r.odes.iter_mut().find(|o| o.id == id) {
                    ode.derivatives = Some(reb_integrator_trace_nbody_derivatives);
                    ode.needs_nbody = 0;
                }
                nbody_ode = Some(id);
                bs.first_or_last_step = 1;
            }

            // In case of overshoot
            if r.t + dt > t_needed {
                dt = t_needed - r.t;
            }

            let mut star = r.particles[0]; // backup velocity
            r.particles[0].vx = 0.; // star does not move in dh
            r.particles[0].vy = 0.;
            r.particles[0].vz = 0.;

            if let Some(id) = nbody_ode {
                if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
                    let mut y = std::mem::take(&mut r.odes[pos].y);
                    for i in 0..encounter_N {
                        let mi = match &r.map {
                            Some(m) => m[i],
                            None => 0,
                        };
                        let p = r.particles[mi];
                        y[i * 6] = p.x;
                        y[i * 6 + 1] = p.y;
                        y[i * 6 + 2] = p.z;
                        y[i * 6 + 3] = p.vx;
                        y[i * 6 + 4] = p.vy;
                        y[i * 6 + 5] = p.vz;
                    }
                    r.odes[pos].y = y;
                }
            }

            let success = reb_integrator_bs_step_odes(r, &mut bs, dt);
            if success != 0 {
                r.t += dt;
            }
            dt = bs.dt_proposed;
            if let Some(id) = nbody_ode {
                if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
                    let y = std::mem::take(&mut r.odes[pos].y);
                    reb_integrator_trace_update_particles(r, &y);
                    r.odes[pos].y = y;
                }
            }

            r.particles[0].vx = star.vx; // restore every timestep for collisions
            r.particles[0].vy = star.vy;
            r.particles[0].vz = star.vz;

            if success != 0 {
                // Only do a collision search for accepted steps.
                reb_collision_search(r);
                if r.N_collisions != 0 {
                    if let reb_integrator_state::trace(ref mut t) = r.integrator {
                        t.force_accept = 1;
                    }
                }
            }

            let p0 = r.particles[0];
            star.vx = p0.vx; // keep track of changed star velocity for later collisions
            star.vy = p0.vy;
            star.vz = p0.vz;
            let _ = star;

            if r.particles[0].x != 0. || r.particles[0].y != 0. || r.particles[0].z != 0. {
                // Collision with star occurred
                // Shift all particles back to heliocentric coordinates
                // Ignore stars velocity:
                //   - will not be used after this
                //   - com velocity is unchanged. this velocity will be used
                //     to reconstruct star's velocity later.
                for i in 0..r.N {
                    r.particles[i].x -= p0.x;
                    r.particles[i].y -= p0.y;
                    r.particles[i].z -= p0.z;
                }
            }
        }

        // Take the state back (hooks may have modified it).
        let trace_back = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
            reb_integrator_state::trace(s) => s,
            other => {
                r.integrator = other;
                reb_integrator_trace_state::default()
            }
        };

        // if only test particles encountered massive bodies, reset the
        // massive body coordinates to their post Kepler step state
        if trace_back.tponly_encounter != 0 {
            for i in 1..trace_back.encounter_N_active {
                let mi = encounter_map_get(r, &trace_back, i);
                r.particles[mi] = trace_back.particles_backup_kepler[mi];
            }
        }

        // Restore odes
        if let Some(id) = nbody_ode {
            reb_ode_free(r, id);
        }
        r.odes = odes_backup;

        r.t = old_t;

        // Resetting BS here reduces binary file size.
        drop(bs); // reb_integrator_bs.free(bs)

        r.integrator = reb_integrator_state::trace(trace_back);
    }

    let mut trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            reb_integrator_trace_state::default()
        }
    };
    trace.encounter_map = r.map.take().unwrap_or_default(); // r->map = NULL
    r.N_map = 0;
    trace
}

/// integrator_trace.c `reb_integrator_trace_kepler_step` (by-value like
/// the BS step it wraps).
pub fn reb_integrator_trace_kepler_step(
    r: &mut reb_simulation,
    mut trace: reb_integrator_trace_state,
    _dt: f64,
) -> reb_integrator_trace_state {
    trace.particles_backup_kepler[..r.N].copy_from_slice(&r.particles[..r.N]);
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, trace, _dt)
}

/// integrator_trace.c `reb_integrator_trace_pre_ts_check`.
pub fn reb_integrator_trace_pre_ts_check(
    r: &mut reb_simulation,
    trace: &mut reb_integrator_trace_state,
) {
    let N = r.N;
    let Nactive = if r.N_active == usize::MAX { r.N } else { r.N_active };
    let _switch: reb_trace_S_fn = match trace.S {
        Some(f) => f,
        None => reb_integrator_trace_switch_default,
    };
    let _switch_peri: reb_trace_S_peri_fn = match trace.S_peri {
        Some(f) => f,
        None => reb_integrator_trace_switch_peri_default,
    };

    // Clear encounter map
    for i in 1..r.N {
        trace.encounter_map[i] = 0;
    }
    trace.encounter_map[0] = 1;
    trace.encounter_N = 1;

    // Reset encounter triggers.
    trace.current_C = 0;

    for i in 0..N {
        for j in (i + 1)..N {
            trace.current_Ks[i * N + j] = 0;
        }
    }

    if r.testparticle_type == 1 {
        trace.tponly_encounter = 0; // testparticles affect massive particles
    } else {
        trace.tponly_encounter = 1;
    }

    // Check for pericenter CE
    for j in 1..Nactive {
        if _switch_peri(r, trace, j) != 0 {
            trace.current_C = 1;
            if trace.peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS
                || trace.peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15
            {
                // Everything will be integrated with BS/IAS15. No need to check any further.
                return;
            }
            if j < Nactive {
                // Two massive particles have a close encounter
                trace.tponly_encounter = 0;
                break; // No need to check other particles
            }
        }
    }

    if trace.current_C != 0 {
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        trace.encounter_N = N;
        for i in 1..N {
            trace.encounter_map[i] = 1; // trigger encounter
        }
    }

    // Body-body
    // there cannot be TP-TP CEs
    for i in 0..Nactive {
        // Check central body, for collisions
        for j in (i + 1)..N {
            if _switch(r, trace, i, j) != 0 {
                trace.current_Ks[i * N + j] = 1;
                if trace.encounter_map[i] == 0 {
                    trace.encounter_map[i] = 1; // trigger encounter
                    trace.encounter_N += 1;
                }
                if trace.encounter_map[j] == 0 {
                    trace.encounter_map[j] = 1; // trigger encounter
                    trace.encounter_N += 1;
                }

                if j < Nactive {
                    // Two massive particles have a close encounter
                    trace.tponly_encounter = 0;
                }
            }
        }
    }
    let (backup, map) = (&mut trace.encounter_map_backup, &trace.encounter_map);
    backup[..N].copy_from_slice(&map[..N]);
}

/// integrator_trace.c `reb_integrator_trace_post_ts_check` — returns 1
/// if any new encounters occurred.
pub fn reb_integrator_trace_post_ts_check(
    r: &mut reb_simulation,
    trace: &mut reb_integrator_trace_state,
) -> usize {
    let N = r.N;
    let Nactive = if r.N_active == usize::MAX { r.N } else { r.N_active };
    let _switch: reb_trace_S_fn = match trace.S {
        Some(f) => f,
        None => reb_integrator_trace_switch_default,
    };
    let _switch_peri: reb_trace_S_peri_fn = match trace.S_peri {
        Some(f) => f,
        None => reb_integrator_trace_switch_peri_default,
    };
    let mut new_close_encounter: usize = 0; // New CEs

    // Set this from pre-ts encounter map. I don't think we need to reset encounter_N here.
    {
        let (map, backup) = (&mut trace.encounter_map, &trace.encounter_map_backup);
        map[..N].copy_from_slice(&backup[..N]);
    }

    if trace.current_C == 0 {
        // Check for pericenter CE if not already triggered from pre-timestep.
        for j in 1..Nactive {
            if _switch_peri(r, trace, j) != 0 {
                trace.current_C = 1;
                new_close_encounter = 1;
                if trace.peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS
                    || trace.peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15
                {
                    // Everything will be integrated with BS/IAS15. No need to check any further.
                    return new_close_encounter;
                }

                if j < Nactive {
                    // Two massive particles have a close encounter
                    trace.tponly_encounter = 0;
                    break; // No need to check other particles
                }
            }
        }
    }
    if trace.current_C != 0 {
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        trace.encounter_N = N;
        for i in 0..N {
            trace.encounter_map[i] = 1; // trigger encounter
        }
    }

    // Body-body
    // there cannot be TP-TP CEs
    for i in 0..Nactive {
        // Do not check for central body anymore
        for j in (i + 1)..N {
            if _switch(r, trace, i, j) != 0 {
                if trace.current_Ks[i * N + j] == 0 {
                    new_close_encounter = 1;
                }
                trace.current_Ks[i * N + j] = 1;
                if trace.encounter_map[i] == 0 {
                    trace.encounter_map[i] = 1; // trigger encounter
                    trace.encounter_N += 1;
                }
                if trace.encounter_map[j] == 0 {
                    trace.encounter_map[j] = 1; // trigger encounter
                    trace.encounter_N += 1;
                }

                if j < Nactive {
                    // Two massive particles have a close encounter
                    trace.tponly_encounter = 0;
                }
            }
        }
    }

    new_close_encounter
}

/// integrator_trace.c static `reb_integrator_trace_step_try` (by-value
/// like the Kepler/BS steps it wraps).
fn reb_integrator_trace_step_try(
    r: &mut reb_simulation,
    mut trace: reb_integrator_trace_state,
) -> reb_integrator_trace_state {
    if trace.current_C == 0 || trace.peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS {
        reb_integrator_trace_interaction_step(r, &mut trace, r.dt / 2.);
        reb_integrator_trace_jump_step(r, &trace, r.dt / 2.);
        let dt = r.dt;
        trace = reb_integrator_trace_kepler_step(r, trace, dt);
        reb_integrator_trace_com_step(&mut trace, r.dt);
        reb_integrator_trace_jump_step(r, &trace, r.dt / 2.);
        reb_integrator_trace_interaction_step(r, &mut trace, r.dt / 2.);
    } else {
        // Pericenter approach with one of the FULL prescriptions
        let t_needed = r.t + r.dt;
        let old_dt = r.dt;
        let old_t = r.t;
        r.gravity = REB_GRAVITY::BASIC;
        trace.mode = REB_INTEGRATOR_TRACE_MODE_FULL;
        reb_integrator_trace_dh_to_inertial(r, &trace);
        match trace.peri_mode {
            REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15 => {
                let mut ias15 = reb_integrator_ias15_state::default();
                while r.t < t_needed && (r.dt / old_dt).abs() > 1e-14 && r.status <= 0 {
                    reb_integrator_ias15_step_state(r, &mut ias15);
                    if r.t + r.dt > t_needed {
                        r.dt = t_needed - r.t;
                    }
                    reb_collision_search(r);
                    if r.N_collisions != 0 {
                        trace.force_accept = 1;
                    }
                }
                drop(ias15); // reb_integrator_ias15.free(ias15)
            }
            REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS => {
                let mut bs = reb_integrator_bs_state::default();
                let mut nbody_ode: Option<usize> = None;

                while r.t < t_needed && (r.dt / old_dt).abs() > 1e-14 && r.status <= 0 {
                    let stale = match nbody_ode {
                        Some(id) => match r.odes.iter().find(|o| o.id == id) {
                            Some(o) => o.length != 6 * r.N,
                            None => true,
                        },
                        None => true,
                    };
                    if stale {
                        // (re)create the ODE
                        if let Some(id) = nbody_ode {
                            reb_ode_free(r, id);
                        }
                        let id = reb_ode_create(r, 6 * r.N);
                        if let Some(ode) = r.odes.iter_mut().find(|o| o.id == id) {
                            ode.derivatives = Some(reb_integrator_bs_nbody_derivatives);
                            ode.needs_nbody = 0;
                        }
                        nbody_ode = Some(id);
                        bs.first_or_last_step = 1;
                    }

                    if let Some(id) = nbody_ode {
                        if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
                            let mut y = std::mem::take(&mut r.odes[pos].y);
                            for i in 0..r.N {
                                let p = r.particles[i];
                                y[i * 6] = p.x;
                                y[i * 6 + 1] = p.y;
                                y[i * 6 + 2] = p.z;
                                y[i * 6 + 3] = p.vx;
                                y[i * 6 + 4] = p.vy;
                                y[i * 6 + 5] = p.vz;
                            }
                            r.odes[pos].y = y;
                        }
                    }

                    let dt_step = r.dt;
                    let success = reb_integrator_bs_step_odes(r, &mut bs, dt_step);
                    if success != 0 {
                        r.t += r.dt;
                    }
                    r.dt = bs.dt_proposed;
                    if r.t + r.dt > t_needed {
                        r.dt = t_needed - r.t;
                    }

                    if let Some(id) = nbody_ode {
                        if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
                            let y = std::mem::take(&mut r.odes[pos].y);
                            reb_integrator_bs_update_particles(r, &y);
                            r.odes[pos].y = y;
                        }
                    }

                    if success != 0 {
                        // Only do a collision search for accepted steps.
                        reb_collision_search(r);
                        if r.N_collisions != 0 {
                            trace.force_accept = 1;
                        }
                    }
                }
                if let Some(id) = nbody_ode {
                    reb_ode_free(r, id);
                }
                drop(bs); // reb_integrator_bs.free(bs)
            }
            _ => {
                reb_simulation_error(r, "Unsupported peri_mode encountered\n");
            }
        }
        r.t = old_t; // final time will be set later
        r.dt = old_dt;
        reb_integrator_trace_inertial_to_dh(r, &mut trace);
    }
    trace
}

/// integrator_trace.c `reb_integrator_trace_did_add_particle` (hook).
pub fn reb_integrator_trace_did_add_particle(r: &mut reb_simulation) {
    // TRACE can add particles mid-timestep now
    let mut trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    if trace.mode == REB_INTEGRATOR_TRACE_MODE_KEPLER {
        let old_N = r.N - 1;
        if trace.N_allocated < r.N {
            trace.current_Ks.resize(r.N * r.N, 0);
            trace.particles_backup.resize(r.N, reb_particle::default());
            trace
                .particles_backup_kepler
                .resize(r.N, reb_particle::default());
            match &mut r.map {
                Some(m) => m.resize(r.N, 0), // r->map = trace->encounter_map (alias)
                None => trace.encounter_map.resize(r.N, 0),
            }
            trace.encounter_map_backup.resize(r.N, 0);
            trace.N_allocated = r.N;
        }

        // First reshuffle existing Ks
        let mut i = old_N;
        while i > 0 {
            i -= 1;
            let mut j = old_N;
            while j > 0 {
                j -= 1;
                trace.current_Ks[i * old_N + j + i] = trace.current_Ks[i * old_N + j];
            }
        }

        // add in new particle, we want it to interact with all currently interacting particles
        // exclude star
        for i in 1..trace.encounter_N {
            let mi = match &r.map {
                Some(m) => m[i],
                None => trace.encounter_map[i],
            };
            trace.current_Ks[mi * r.N + old_N] = 1;
        }

        let encounter_N = trace.encounter_N;
        match &mut r.map {
            Some(m) => m[encounter_N] = old_N,
            None => trace.encounter_map[encounter_N] = old_N,
        }
        trace.encounter_N += 1;
        r.N_map += 1;

        if r.N_active == usize::MAX {
            // If global N_active is not set, then all particles are active, so the new one as well.
            // Otherwise, assume we're adding non active particle.
            trace.encounter_N_active += 1;
        }
    }
    r.integrator = reb_integrator_state::trace(trace);
}

/// integrator_trace.c `reb_integrator_trace_will_remove_particle` (hook).
pub fn reb_integrator_trace_will_remove_particle(r: &mut reb_simulation, index: usize) {
    let mut trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    if trace.mode == REB_INTEGRATOR_TRACE_MODE_KEPLER {
        // Only removed mid-timestep if collision - BS Step!
        let mut after_to_be_removed_particle = 0;
        let mut encounter_index = usize::MAX;
        {
            let map: &mut Vec<usize> = match r.map.as_mut() {
                Some(m) => m,
                None => &mut trace.encounter_map,
            };
            for i in 0..trace.encounter_N {
                if after_to_be_removed_particle == 1 {
                    map[i - 1] = map[i] - 1;
                }
                if map[i] == index {
                    encounter_index = i;
                    after_to_be_removed_particle = 1;
                }
            }
        }
        if encounter_index == usize::MAX {
            reb_simulation_error(r, "Cannot find particle in encounter map. Did not remove particle.");
            r.integrator = reb_integrator_state::trace(trace);
            return;
        }

        // reshuffle current_Ks
        let mut counter: usize = 0;
        let new_N = r.N - 1;
        for i in 0..new_N {
            if i == index {
                counter += r.N;
            }
            for j in 0..new_N {
                if j == index {
                    counter += 1;
                }
                trace.current_Ks[i * new_N + j] = trace.current_Ks[i * new_N + j + counter];
            }
        }
        if encounter_index < trace.encounter_N_active {
            trace.encounter_N_active -= 1;
        }
        trace.encounter_N -= 1;
        r.N_map -= 1;
    }
    r.integrator = reb_integrator_state::trace(trace);
}

/// integrator_trace.c `reb_integrator_trace_step` (state-explicit).
pub fn reb_integrator_trace_step_state(
    r: &mut reb_simulation,
    mut trace: reb_integrator_trace_state,
) -> reb_integrator_trace_state {
    // Do memory management and consistency checks
    let N = r.N;

    if r.N_var != 0 {
        reb_simulation_warning(r, "TRACE does not work with variational equations.");
    }

    if trace.N_allocated < N {
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        trace.particles_backup.resize(N, reb_particle::default());
        trace
            .particles_backup_kepler
            .resize(N, reb_particle::default());
        trace.current_Ks.resize(N * N, 0);
        trace.encounter_map.resize(N, 0);
        trace.encounter_map_backup.resize(N, 0);
        trace.N_allocated = N;
    }

    // Calculate collisions only with DIRECT or LINE method
    if r.collision != REB_COLLISION::NONE
        && (r.collision != REB_COLLISION::DIRECT && r.collision != REB_COLLISION::LINE)
    {
        reb_simulation_warning(r, "TRACE only works with a direct or line collision search.");
    }
    r.N_targets = usize::MAX; // Search for collisions between all particles in encounter step or full steps.

    // Calculate gravity with special function
    if r.gravity != REB_GRAVITY::BASIC && r.gravity != REB_GRAVITY::CUSTOM {
        reb_simulation_warning(
            r,
            "TRACE has its own gravity routine. Gravity routine set by the user will be ignored.",
        );
    }

    reb_integrator_trace_inertial_to_dh(r, &mut trace);

    // Create copy of all particle to allow for the step to be rejected.
    trace.particles_backup[..N].copy_from_slice(&r.particles[..N]);

    // This will be set to 1 if a collision occurred.
    trace.force_accept = 0;

    // Check if there are any close encounters
    reb_integrator_trace_pre_ts_check(r, &mut trace);

    // Attempt one step.
    let mut trace = reb_integrator_trace_step_try(r, trace);

    // We always accept the step if a collision occurred as it is impossible to undo the collision.
    if trace.force_accept == 0 {
        // We check again for close encounters to ensure time reversibility.
        if reb_integrator_trace_post_ts_check(r, &mut trace) != 0 {
            // New encounters were found. Will reject the step.
            // Revert particles to the beginning of the step.
            r.particles[..N].copy_from_slice(&trace.particles_backup[..N]);

            // Do step again
            trace = reb_integrator_trace_step_try(r, trace);
        }
    }
    reb_integrator_trace_dh_to_inertial(r, &trace);

    r.t += r.dt;
    r.dt_last_done = r.dt;
    r.N_targets = 1; // Only search for collisions with star after complete timestep.
    trace
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_trace_step(r: &mut reb_simulation) {
    let trace = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::trace(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    let trace = reb_integrator_trace_step_state(r, trace);
    r.integrator = reb_integrator_state::trace(trace);
}
