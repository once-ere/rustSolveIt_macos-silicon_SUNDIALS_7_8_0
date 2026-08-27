//! integrator_mercurius.rs — the hybrid symplectic MERCURIUS integrator
//! (from integrator_mercurius.c/h; Rein et al. 2019, after Chambers
//! 1999). WHFast democratic-heliocentric splitting for the long-term
//! evolution, switching smoothly to IAS15 for close encounters.
//!
//! Ownership notes (all Rust mechanics, no arithmetic changes):
//! - The C keeps `r->map` as an alias of `mercurius->encounter_map`
//!   while the encounter sub-integration runs. Here the Vec is *moved*
//!   into `r.map` for the duration of the encounter loop and moved back
//!   afterwards, so the collision/IAS15 code sees exactly the array the
//!   C sees.
//! - The encounter gravity routine is installed as `r.gravity_custom`
//!   (a plain fn pointer, as in C). It reaches the mercurius state by
//!   temporarily taking it out of `r.integrator`, which holds the state
//!   during the encounter loop.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Dan Tamayo and contributors. See crate root.

use crate::collision::reb_collision_search;
use crate::integrator_ias15::{reb_integrator_ias15_state, reb_integrator_ias15_step_state};
use crate::integrator_whfast::reb_integrator_whfast_kepler_solver;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;

/// integrator_mercurius.c `#define MIN(a, b)` — the exact ternary.
fn MIN(a: f64, b: f64) -> f64 {
    if a > b {
        b
    } else {
        a
    }
}

/// integrator_mercurius.c `#define MAX(a, b)` — the exact ternary.
fn MAX(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// Switching-function pointer type (C: `double (*L)(const struct
/// reb_simulation* const, double d, double dcrit)`).
pub type reb_integrator_mercurius_Lfunc = fn(&reb_simulation, f64, f64) -> f64;

/// integrator_mercurius.h anonymous mode enum.
pub const REB_INTEGRATOR_MERCURIUS_MODE_WH: u32 = 0;
pub const REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER: u32 = 1;

/// integrator_mercurius.h `struct reb_integrator_mercurius_state`.
/// `dcrit.len()` plays the role of `N_allocated_dcrit` and
/// `particles_backup_additional_forces.len()` the role of
/// `N_allocated_additional_forces` (realloc only ever grows them, so
/// the Vec length carries the same information).
#[derive(Clone, Debug)]
pub struct reb_integrator_mercurius_state {
    /// Switching function (default same as Mercury).
    pub L: reb_integrator_mercurius_Lfunc,
    /// Critical switching distance in units of Hill radii.
    pub r_crit_hill: f64,
    /// Combine kick steps at beginning and end of timestep.
    pub safe_mode: u32,
    // Internal use
    pub mode: u32,
    /// Number of particles currently having an encounter.
    pub encounter_N: usize,
    /// Number of active particles currently having an encounter.
    pub encounter_N_active: usize,
    /// 0 if any encounters are between two massive bodies. 1 if
    /// encounters only involve test particles.
    pub tponly_encounter: u32,
    /// C `N_allocated` (size of particles_backup / encounter_map).
    pub N_allocated: usize,
    /// Precalculated switching radii for particles.
    pub dcrit: Vec<f64>,
    /// Coordinates before Kepler step for encounter prediction.
    pub particles_backup: Vec<reb_particle>,
    /// Coordinates backup used around additional_forces evaluation.
    pub particles_backup_additional_forces: Vec<reb_particle>,
    /// Map to represent which particles are integrated with ias15.
    /// Moved into `r.map` while the encounter sub-integration runs.
    pub encounter_map: Vec<usize>,
    /// Used to keep track of the center of mass during the timestep.
    pub com_pos: reb_vec3d,
    pub com_vel: reb_vec3d,
}

impl Default for reb_integrator_mercurius_state {
    /// integrator_mercurius.c `reb_integrator_mercurius_create`.
    fn default() -> Self {
        reb_integrator_mercurius_state {
            L: reb_integrator_mercurius_L_mercury,
            r_crit_hill: 3.,
            safe_mode: 1,
            mode: REB_INTEGRATOR_MERCURIUS_MODE_WH,
            encounter_N: 0,
            encounter_N_active: 0,
            tponly_encounter: 0,
            N_allocated: 0,
            dcrit: Vec::new(),
            particles_backup: Vec::new(),
            particles_backup_additional_forces: Vec::new(),
            encounter_map: Vec::new(),
            com_pos: reb_vec3d::default(),
            com_vel: reb_vec3d::default(),
        }
    }
}

/// integrator_mercurius.c `reb_integrator_mercurius_L_mercury` — the
/// changeover function used by the Mercury integrator.
pub fn reb_integrator_mercurius_L_mercury(_r: &reb_simulation, d: f64, dcrit: f64) -> f64 {
    let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
    if y < 0. {
        0.
    } else if y > 1. {
        1.
    } else {
        10. * (y * y * y) - 15. * (y * y * y * y) + 6. * (y * y * y * y * y)
    }
}

/// integrator_mercurius.c `reb_integrator_mercurius_L_C4` — the
/// changeover function C4 proposed by Hernandez (2019).
pub fn reb_integrator_mercurius_L_C4(_r: &reb_simulation, d: f64, dcrit: f64) -> f64 {
    let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
    if y < 0. {
        0.
    } else if y > 1. {
        1.
    } else {
        (70. * y * y * y * y - 315. * y * y * y + 540. * y * y - 420. * y + 126.) * y * y * y * y * y
    }
}

/// integrator_mercurius.c `reb_integrator_mercurius_L_C5` — the
/// changeover function C5 proposed by Hernandez (2019).
pub fn reb_integrator_mercurius_L_C5(_r: &reb_simulation, d: f64, dcrit: f64) -> f64 {
    let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
    if y < 0. {
        0.
    } else if y > 1. {
        1.
    } else {
        (-252. * y * y * y * y * y + 1386. * y * y * y * y - 3080. * y * y * y + 3465. * y * y
            - 1980. * y
            + 462.)
            * y * y * y * y * y * y
    }
}

/// integrator_mercurius.c static `f`.
fn f(x: f64) -> f64 {
    if x < 0. {
        return 0.;
    }
    (-1. / x).exp()
}

/// integrator_mercurius.c `reb_integrator_mercurius_L_infinity` —
/// infinitely differentiable switching function.
pub fn reb_integrator_mercurius_L_infinity(_r: &reb_simulation, d: f64, dcrit: f64) -> f64 {
    let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
    if y < 0. {
        0.
    } else if y > 1. {
        1.
    } else {
        f(y) / (f(y) + f(1. - y))
    }
}

/// integrator_mercurius.c `reb_integrator_mercurius_inertial_to_dh`.
pub fn reb_integrator_mercurius_inertial_to_dh(
    r: &mut reb_simulation,
    mercurius: &mut reb_integrator_mercurius_state,
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
    mercurius.com_pos = com_pos;
    mercurius.com_vel = com_vel;
}

/// integrator_mercurius.c `reb_integrator_mercurius_dh_to_inertial`.
pub fn reb_integrator_mercurius_dh_to_inertial(
    r: &mut reb_simulation,
    mercurius: &mut reb_integrator_mercurius_state,
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
    particles[0].x = mercurius.com_pos.x - temp.x;
    particles[0].y = mercurius.com_pos.y - temp.y;
    particles[0].z = mercurius.com_pos.z - temp.z;

    for i in 1..N {
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += mercurius.com_vel.x;
        particles[i].vy += mercurius.com_vel.y;
        particles[i].vz += mercurius.com_vel.z;
    }
    particles[0].vx = mercurius.com_vel.x - temp.vx;
    particles[0].vy = mercurius.com_vel.y - temp.vy;
    particles[0].vz = mercurius.com_vel.z - temp.vz;
}

/// integrator_mercurius.c static `reb_mercurius_encounter_predict` —
/// predicts close encounters during the timestep from the old and new
/// positions/velocities around the Kepler step.
fn reb_mercurius_encounter_predict(
    r: &mut reb_simulation,
    mercurius: &mut reb_integrator_mercurius_state,
) {
    let particles = &r.particles;
    let particles_backup = &mercurius.particles_backup;
    let dcrit = &mercurius.dcrit;
    let N = r.N;
    let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };
    let dt = r.dt;
    let mut encounter_N: usize = 1;
    mercurius.encounter_map[0] = 1;
    let mut tponly_encounter: u32 = if r.testparticle_type == 1 {
        0 // testparticles affect massive particles
    } else {
        1
    };
    for i in 1..N {
        mercurius.encounter_map[i] = 0;
    }
    for i in 0..N_active {
        for j in (i + 1)..N {
            let dxn = particles[i].x - particles[j].x;
            let dyn_ = particles[i].y - particles[j].y;
            let dzn = particles[i].z - particles[j].z;
            let dvxn = particles[i].vx - particles[j].vx;
            let dvyn = particles[i].vy - particles[j].vy;
            let dvzn = particles[i].vz - particles[j].vz;
            let rn = dxn * dxn + dyn_ * dyn_ + dzn * dzn;
            let dxo = particles_backup[i].x - particles_backup[j].x;
            let dyo = particles_backup[i].y - particles_backup[j].y;
            let dzo = particles_backup[i].z - particles_backup[j].z;
            let dvxo = particles_backup[i].vx - particles_backup[j].vx;
            let dvyo = particles_backup[i].vy - particles_backup[j].vy;
            let dvzo = particles_backup[i].vz - particles_backup[j].vz;
            let ro = dxo * dxo + dyo * dyo + dzo * dzo;

            let drndt = (dxn * dvxn + dyn_ * dvyn + dzn * dvzn) * 2.;
            let drodt = (dxo * dvxo + dyo * dvyo + dzo * dvzo) * 2.;

            let a = 6. * (ro - rn) + 3. * dt * (drodt + drndt);
            let b = 6. * (rn - ro) - 2. * dt * (2. * drodt + drndt);
            let c = dt * drodt;

            let mut rmin = MIN(rn, ro);

            let s = b * b - 4. * a * c;
            let sr = MAX(0., s).sqrt();
            let tmin1 = (-b + sr) / (2. * a);
            let tmin2 = (-b - sr) / (2. * a);
            if tmin1 > 0. && tmin1 < 1. {
                let rmin1 = (1. - tmin1) * (1. - tmin1) * (1. + 2. * tmin1) * ro
                    + tmin1 * tmin1 * (3. - 2. * tmin1) * rn
                    + tmin1 * (1. - tmin1) * (1. - tmin1) * dt * drodt
                    - tmin1 * tmin1 * (1. - tmin1) * dt * drndt;
                rmin = MIN(MAX(rmin1, 0.), rmin);
            }
            if tmin2 > 0. && tmin2 < 1. {
                let rmin2 = (1. - tmin2) * (1. - tmin2) * (1. + 2. * tmin2) * ro
                    + tmin2 * tmin2 * (3. - 2. * tmin2) * rn
                    + tmin2 * (1. - tmin2) * (1. - tmin2) * dt * drodt
                    - tmin2 * tmin2 * (1. - tmin2) * dt * drndt;
                rmin = MIN(MAX(rmin2, 0.), rmin);
            }

            let mut dcritmax2 = MAX(dcrit[i], dcrit[j]);
            dcritmax2 *= 1.21 * dcritmax2;
            if rmin < dcritmax2 {
                if mercurius.encounter_map[i] == 0 {
                    mercurius.encounter_map[i] = i;
                    encounter_N += 1;
                }
                if mercurius.encounter_map[j] == 0 {
                    mercurius.encounter_map[j] = j;
                    encounter_N += 1;
                }
                if j < N_active {
                    // Two massive particles have a close encounter
                    tponly_encounter = 0;
                }
            }
        }
    }
    mercurius.encounter_N = encounter_N;
    mercurius.tponly_encounter = tponly_encounter;
}

/// integrator_mercurius.c `reb_integrator_mercurius_interaction_step`.
pub fn reb_integrator_mercurius_interaction_step(r: &mut reb_simulation, dt: f64) {
    let particles = &mut r.particles;
    let N = r.N;
    for i in 1..N {
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
    }
}

/// integrator_mercurius.c `reb_integrator_mercurius_jump_step`.
pub fn reb_integrator_mercurius_jump_step(r: &mut reb_simulation, dt: f64) {
    let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };
    let N = if r.testparticle_type == 0 { N_active } else { r.N };
    let mut px = 0.;
    let mut py = 0.;
    let mut pz = 0.;
    for i in 1..N {
        px += r.particles[i].vx * r.particles[i].m; // in dh
        py += r.particles[i].vy * r.particles[i].m;
        pz += r.particles[i].vz * r.particles[i].m;
    }
    px /= r.particles[0].m;
    py /= r.particles[0].m;
    pz /= r.particles[0].m;
    let N_all = r.N;
    let particles = &mut r.particles;
    for i in 1..N_all {
        particles[i].x += dt * px;
        particles[i].y += dt * py;
        particles[i].z += dt * pz;
    }
}

/// integrator_mercurius.c static `reb_integrator_mercurius_kepler_step`.
/// The C calls the WHFast Kepler solver directly on `r->particles`; the
/// Rust solver takes the particle buffer explicitly, so the Vec is
/// temporarily moved out of the simulation (pure ownership mechanics).
fn reb_integrator_mercurius_kepler_step(r: &mut reb_simulation, dt: f64) {
    let N = r.N;
    let mu = r.G * r.particles[0].m;
    let mut particles = std::mem::take(&mut r.particles);
    let mut no_var: [reb_particle; 0] = [];
    for i in 1..N {
        reb_integrator_whfast_kepler_solver(Some(&mut *r), &mut particles, &mut no_var, i, mu, dt); // in dh
    }
    r.particles = particles;
}

/// Fetch entry `i` of the encounter map. During the encounter loop the
/// map Vec lives in `r.map` (the C aliases `r->map` to
/// `mercurius->encounter_map`); otherwise it is still in the state.
fn encounter_map_get(
    r: &reb_simulation,
    mercurius: &reb_integrator_mercurius_state,
    i: usize,
) -> usize {
    match &r.map {
        Some(m) => m[i],
        None => mercurius.encounter_map[i],
    }
}

/// integrator_mercurius.c static `reb_mercurius_encounter_step` — only
/// particles having a close encounter are integrated by IAS15. Takes
/// and returns the state by value: while the IAS15 sub-loop runs, the
/// state is stored in `r.integrator` so that the custom gravity routine
/// and the add/remove particle hooks can reach it (the C shares it via
/// `r->integrator.state`).
fn reb_mercurius_encounter_step(
    r: &mut reb_simulation,
    mut mercurius: reb_integrator_mercurius_state,
    _dt: f64,
) -> reb_integrator_mercurius_state {
    if mercurius.encounter_N < 2 {
        return mercurius; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }
    let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };
    let mut i_enc: usize = 0;
    mercurius.encounter_N_active = 0;
    for i in 0..r.N {
        if mercurius.encounter_map[i] != 0 {
            let tmp = r.particles[i]; // Copy for potential use for tponly_encounter
            r.particles[i] = mercurius.particles_backup[i]; // Use coordinates before whfast step
            mercurius.encounter_map[i_enc] = i;
            i_enc += 1;
            if i < N_active {
                mercurius.encounter_N_active += 1;
                if mercurius.tponly_encounter != 0 {
                    mercurius.particles_backup[i] = tmp; // Make copy of particles after the kepler step.
                                                         // used to restore the massive objects' states in the case
                                                         // of only massless test-particle encounters
                }
            }
        }
    }

    mercurius.mode = REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER;

    // run
    let old_dt = r.dt;
    let dtsign = if old_dt >= 0. { 1. } else { -1. };
    let old_t = r.t;
    let t_needed = r.t + _dt;

    let mut ias15 = reb_integrator_ias15_state::default(); // reb_integrator_ias15.create()
    r.map = Some(std::mem::take(&mut mercurius.encounter_map)); // r->map = mercurius->encounter_map
    r.N_map = mercurius.encounter_N;
    r.N_targets = usize::MAX; // Search for any possible collisions between N_map particles.

    r.dt = 0.0001 * _dt; // start with a small timestep.

    // No additional forces during encounter
    r.gravity = REB_GRAVITY::CUSTOM;
    r.gravity_custom = Some(reb_integrator_mercurius_calculate_acceleration_mode_encounter);

    // Store the state where the gravity routine and particle hooks find it.
    r.integrator = reb_integrator_state::mercurius(mercurius);

    while dtsign * r.t < dtsign * t_needed && (r.dt / old_dt).abs() > 1e-14 && r.status <= 0 {
        let mut star = r.particles[0]; // backup velocity
        r.particles[0].vx = 0.; // star does not move in dh
        r.particles[0].vy = 0.;
        r.particles[0].vz = 0.;
        reb_integrator_ias15_step_state(r, &mut ias15);
        r.particles[0].vx = star.vx; // restore every timestep for collisions
        r.particles[0].vy = star.vy;
        r.particles[0].vz = star.vz;

        if dtsign * (r.t + r.dt) > dtsign * t_needed {
            r.dt = t_needed - r.t;
        }

        // Search and resolve collisions
        reb_collision_search(r);

        let p0 = r.particles[0];
        star.vx = p0.vx; // keep track of changed star velocity for later collisions
        star.vy = p0.vy;
        star.vz = p0.vz;
        let _ = star;
        if p0.x != 0. || p0.y != 0. || p0.z != 0. {
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

    // Take the state back (the hooks may have modified it).
    let mut mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return reb_integrator_mercurius_state::default();
        }
    };

    // if only test particles encountered massive bodies, reset the
    // massive body coordinates to their post Kepler step state
    if mercurius.tponly_encounter != 0 {
        for i in 1..mercurius.encounter_N_active {
            let mi = encounter_map_get(r, &mercurius, i);
            r.particles[mi] = mercurius.particles_backup[mi];
        }
    }

    drop(ias15); // reb_integrator_ias15.free(ias15)

    r.t = old_t;
    r.dt = old_dt;
    mercurius.mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;
    mercurius.encounter_map = r.map.take().unwrap_or_default(); // r->map = NULL
    r.N_map = 0;
    mercurius
}

/// integrator_mercurius.c
/// `reb_integrator_mercurius_calculate_dcrit_for_particle`.
pub fn reb_integrator_mercurius_calculate_dcrit_for_particle(
    r: &reb_simulation,
    mercurius: &reb_integrator_mercurius_state,
    i: usize,
) -> f64 {
    let m0 = r.particles[0].m;
    let dx = r.particles[i].x; // in dh
    let dy = r.particles[i].y;
    let dz = r.particles[i].z;
    let dvx = r.particles[i].vx - r.particles[0].vx;
    let dvy = r.particles[i].vy - r.particles[0].vy;
    let dvz = r.particles[i].vz - r.particles[0].vz;
    let _r = (dx * dx + dy * dy + dz * dz).sqrt();
    let v2 = dvx * dvx + dvy * dvy + dvz * dvz;

    let GM = r.G * (m0 + r.particles[i].m);
    let a = GM * _r / (2. * GM - _r * v2);
    let vc = (GM / a.abs()).sqrt();
    let mut dcrit = 0.;
    // Criteria 1: average velocity
    dcrit = MAX(dcrit, vc * 0.4 * r.dt);
    // Criteria 2: current velocity
    dcrit = MAX(dcrit, v2.sqrt() * 0.4 * r.dt);
    // Criteria 3: Hill radius
    dcrit = MAX(
        dcrit,
        mercurius.r_crit_hill * a * (r.particles[i].m / (3. * r.particles[0].m)).cbrt(),
    );
    // Criteria 4: physical radius
    dcrit = MAX(dcrit, 2. * r.particles[i].r);
    dcrit
}

/// integrator_mercurius.c
/// `reb_integrator_mercurius_calculate_acceleration_mode_encounter`.
/// Installed as `r.gravity_custom` during the encounter sub-integration;
/// fetches the mercurius state from `r.integrator` (C reads
/// `r->integrator.state`).
pub fn reb_integrator_mercurius_calculate_acceleration_mode_encounter(r: &mut reb_simulation) {
    let mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let _testparticle_type = r.testparticle_type;
    let _L = mercurius.L;
    let m0 = r.particles[0].m;
    let encounter_N = mercurius.encounter_N;
    let encounter_N_active = mercurius.encounter_N_active;
    // (OPENMP branch of the C omitted: the reference Windows build is
    // compiled without OpenMP.)
    r.particles[0].ax = 0.; // map[0] is always 0
    r.particles[0].ay = 0.;
    r.particles[0].az = 0.;
    // Acceleration due to star
    for i in 1..encounter_N {
        let mi = encounter_map_get(r, &mercurius, i);
        let x = r.particles[mi].x;
        let y = r.particles[mi].y;
        let z = r.particles[mi].z;
        let _r = (x * x + y * y + z * z + softening2).sqrt();
        let prefact = -G / (_r * _r * _r) * m0;
        r.particles[mi].ax = prefact * x;
        r.particles[mi].ay = prefact * y;
        r.particles[mi].az = prefact * z;
    }
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Interactions between active-active
    for i in 2..encounter_N_active {
        let mi = encounter_map_get(r, &mercurius, i);
        for j in 1..i {
            let mj = encounter_map_get(r, &mercurius, j);
            let dx = r.particles[mi].x - r.particles[mj].x;
            let dy = r.particles[mi].y - r.particles[mj].y;
            let dz = r.particles[mi].z - r.particles[mj].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let dcritmax = MAX(mercurius.dcrit[mi], mercurius.dcrit[mj]);
            let L = _L(r, _r, dcritmax);
            let prefact = G * (1. - L) / (_r * _r * _r);
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
    // Interactions between active-testparticle
    let startitestp = std::cmp::max(encounter_N_active, 2);
    for i in startitestp..encounter_N {
        let mi = encounter_map_get(r, &mercurius, i);
        for j in 1..encounter_N_active {
            let mj = encounter_map_get(r, &mercurius, j);
            let dx = r.particles[mi].x - r.particles[mj].x;
            let dy = r.particles[mi].y - r.particles[mj].y;
            let dz = r.particles[mi].z - r.particles[mj].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let dcritmax = MAX(mercurius.dcrit[mi], mercurius.dcrit[mj]);
            let L = _L(r, _r, dcritmax);
            let prefact = G * (1. - L) / (_r * _r * _r);
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
    r.integrator = reb_integrator_state::mercurius(mercurius);
}

/// integrator_mercurius.c static
/// `reb_integrator_mercurius_calculate_acceleration_mode_wh`.
fn reb_integrator_mercurius_calculate_acceleration_mode_wh(
    r: &mut reb_simulation,
    mercurius: &mut reb_integrator_mercurius_state,
) {
    let N = r.N;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let _testparticle_type = r.testparticle_type;
    let _L = mercurius.L;
    // (OPENMP branch and the `reb_sigint` signal check of the C omitted:
    // the reference Windows build has neither.)
    for i in 0..N {
        r.particles[i].ax = 0.;
        r.particles[i].ay = 0.;
        r.particles[i].az = 0.;
    }
    for i in 2..N_active {
        for j in 1..i {
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let dcritmax = MAX(mercurius.dcrit[i], mercurius.dcrit[j]);
            let L = _L(r, _r, dcritmax);
            let prefact = G * L / (_r * _r * _r);
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
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
            let dcritmax = MAX(mercurius.dcrit[i], mercurius.dcrit[j]);
            let L = _L(r, _r, dcritmax);
            let prefact = G * L / (_r * _r * _r);
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

    if r.additional_forces.is_some() {
        // Additional forces are only calculated in the kick step, not during close encounter
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        if r.N > mercurius.particles_backup_additional_forces.len() {
            mercurius
                .particles_backup_additional_forces
                .resize(r.N, reb_particle::default());
        }
        mercurius.particles_backup_additional_forces[..r.N].copy_from_slice(&r.particles[..r.N]);
        reb_integrator_mercurius_dh_to_inertial(r, mercurius);

        (r.additional_forces.unwrap())(r);

        let particles = &mut r.particles;
        let backup = &mercurius.particles_backup_additional_forces;
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

/// integrator_mercurius.c `reb_integrator_mercurius_did_add_particle`
/// (hook, dispatched from `reb_integrator_did_add_particle`).
pub fn reb_integrator_mercurius_did_add_particle(r: &mut reb_simulation) {
    let mut mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    match mercurius.mode {
        REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER => {
            if mercurius.dcrit.len() < r.N {
                mercurius.dcrit.resize(r.N, 0.);
            }
            let dcrit_new =
                reb_integrator_mercurius_calculate_dcrit_for_particle(r, &mercurius, r.N - 1);
            mercurius.dcrit[r.N - 1] = dcrit_new;
            if mercurius.N_allocated < r.N {
                mercurius.particles_backup.resize(r.N, reb_particle::default());
                match &mut r.map {
                    Some(m) => m.resize(r.N, 0), // r->map = mercurius->encounter_map (alias)
                    None => mercurius.encounter_map.resize(r.N, 0),
                }
                mercurius.N_allocated = r.N;
            }
            let encounter_N = mercurius.encounter_N;
            match &mut r.map {
                Some(m) => m[encounter_N] = r.N - 1,
                None => mercurius.encounter_map[encounter_N] = r.N - 1,
            }
            mercurius.encounter_N += 1;
            r.N_map += 1;
            if r.N_active == usize::MAX {
                // If global N_active is not set, then all particles are active, so the new one as well.
                // Otherwise, assume we're adding non active particle.
                mercurius.encounter_N_active += 1;
            }
        }
        _ => {
            // REB_INTEGRATOR_MERCURIUS_MODE_WH:
            // Nothing to do here. r->did_modify_particles will get set automatically
        }
    }
    r.integrator = reb_integrator_state::mercurius(mercurius);
}

/// integrator_mercurius.c `reb_integrator_mercurius_will_remove_particle`
/// (hook, dispatched from `reb_integrator_will_remove_particle`).
pub fn reb_integrator_mercurius_will_remove_particle(r: &mut reb_simulation, index: usize) {
    let mut mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    if !mercurius.dcrit.is_empty() && index < mercurius.dcrit.len() {
        for i in 0..(r.N - 1) {
            if i >= index {
                mercurius.dcrit[i] = mercurius.dcrit[i + 1];
            }
        }
    }
    if mercurius.mode == REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER {
        let mut after_to_be_removed_particle = 0;
        let mut encounter_index = usize::MAX;
        {
            let map: &mut Vec<usize> = match r.map.as_mut() {
                Some(m) => m,
                None => &mut mercurius.encounter_map,
            };
            for i in 0..mercurius.encounter_N {
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
            reb_simulation_error(r, "Cannot find particle in encounter map.");
            r.integrator = reb_integrator_state::mercurius(mercurius);
            return;
        }
        if encounter_index < mercurius.encounter_N_active {
            mercurius.encounter_N_active -= 1;
        }
        mercurius.encounter_N -= 1;
        r.N_map -= 1;
    }
    r.integrator = reb_integrator_state::mercurius(mercurius);
}

/// integrator_mercurius.c `reb_integrator_mercurius_step` (state-explicit).
pub fn reb_integrator_mercurius_step_state(
    r: &mut reb_simulation,
    mut mercurius: reb_integrator_mercurius_state,
) -> reb_integrator_mercurius_state {
    if !r.var_config.is_empty() {
        reb_simulation_warning(r, "Mercurius does not work with variational equations.");
    }

    let N = r.N;

    if mercurius.dcrit.len() < N {
        // Need to safe these arrays in Simulationarchive
        mercurius.dcrit.resize(N, 0.);
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        r.did_modify_particles = 1;
    }
    if mercurius.N_allocated < N {
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility
        mercurius.particles_backup.resize(N, reb_particle::default());
        mercurius.encounter_map.resize(N, 0);
        mercurius.N_allocated = N;
    }
    if mercurius.safe_mode != 0 || r.did_modify_particles != 0 {
        if r.is_synchronized == 0 {
            reb_integrator_mercurius_synchronize_state(r, &mut mercurius);
            reb_simulation_warning(
                r,
                "Particles were modified while simulation was not synchronized.",
            );
        }
        reb_integrator_mercurius_inertial_to_dh(r, &mut mercurius);
    }

    if r.did_modify_particles != 0 {
        if r.is_synchronized == 0 {
            reb_integrator_mercurius_synchronize_state(r, &mut mercurius);
            reb_integrator_mercurius_inertial_to_dh(r, &mut mercurius);
            reb_simulation_warning(
                r,
                "MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.",
            );
        }
        mercurius.dcrit[0] = 2. * r.particles[0].r; // central object only uses physical radius
        for i in 1..N {
            let dcrit_i = reb_integrator_mercurius_calculate_dcrit_for_particle(r, &mercurius, i);
            mercurius.dcrit[i] = dcrit_i;
        }
    }

    // Calculate collisions only with DIRECT method
    if r.collision != REB_COLLISION::NONE && r.collision != REB_COLLISION::DIRECT {
        reb_simulation_warning(r, "Mercurius only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if r.gravity != REB_GRAVITY::BASIC && r.gravity != REB_GRAVITY::CUSTOM {
        reb_simulation_warning(
            r,
            "Mercurius has its own gravity routine. Gravity routine set by the user will be ignored.",
        );
    }
    mercurius.mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;

    reb_integrator_mercurius_calculate_acceleration_mode_wh(r, &mut mercurius);

    if r.is_synchronized != 0 {
        reb_integrator_mercurius_interaction_step(r, r.dt / 2.);
    } else {
        reb_integrator_mercurius_interaction_step(r, r.dt);
    }
    reb_integrator_mercurius_jump_step(r, r.dt / 2.);

    // COM step
    mercurius.com_pos.x += r.dt * mercurius.com_vel.x;
    mercurius.com_pos.y += r.dt * mercurius.com_vel.y;
    mercurius.com_pos.z += r.dt * mercurius.com_vel.z;

    // Make copy of particles before the kepler step.
    // Then evolve all particles in kepler step.
    // Result will be used in encounter prediction.
    // Particles having a close encounter will be overwritten
    // later by encounter step.
    mercurius.particles_backup[..N].copy_from_slice(&r.particles[..N]);
    reb_integrator_mercurius_kepler_step(r, r.dt);

    reb_mercurius_encounter_predict(r, &mut mercurius);

    let dt = r.dt;
    let mut mercurius = reb_mercurius_encounter_step(r, mercurius, dt);

    reb_integrator_mercurius_jump_step(r, r.dt / 2.);

    r.is_synchronized = 0;
    if mercurius.safe_mode != 0 {
        reb_integrator_mercurius_synchronize_state(r, &mut mercurius);
    }

    r.t += r.dt;
    r.dt_last_done = r.dt;
    r.N_targets = 1; // Only search for collisions with star in-between timesteps.
    mercurius
}

/// integrator_mercurius.c `reb_integrator_mercurius_synchronize`
/// (state-explicit).
pub fn reb_integrator_mercurius_synchronize_state(
    r: &mut reb_simulation,
    mercurius: &mut reb_integrator_mercurius_state,
) {
    if r.is_synchronized == 0 {
        mercurius.mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;
        reb_integrator_mercurius_calculate_acceleration_mode_wh(r, mercurius);
        reb_integrator_mercurius_interaction_step(r, r.dt / 2.);

        reb_integrator_mercurius_dh_to_inertial(r, mercurius);

        r.is_synchronized = 1;
    }
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_mercurius_step(r: &mut reb_simulation) {
    let mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    let mercurius = reb_integrator_mercurius_step_state(r, mercurius);
    r.integrator = reb_integrator_state::mercurius(mercurius);
}

/// Synchronize entry point for the dispatcher.
pub fn reb_integrator_mercurius_synchronize(r: &mut reb_simulation) {
    let mut mercurius = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::mercurius(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_mercurius_synchronize_state(r, &mut mercurius);
    r.integrator = reb_integrator_state::mercurius(mercurius);
}
