//! modify_orbits_forces.rs — translation of REBOUNDx modify_orbits_forces.c
//! Update orbital elements with prescribed timescales using forces.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Header (modify_orbits_forces.c)
//!
//! ```text
//! @file    modify_orbits_forces.c
//! @brief   Update orbital elements with prescribed timescales using forces.
//! @author  Dan Tamayo <tamayo.daniel@gmail.com>
//!
//! $Orbit Modifications$
//!
//! ======================= ===============================================
//! Authors                 D. Tamayo, H. Rein
//! Implementation Paper    Kostov et al., 2016
//! Based on                Papaloizou & Larwood 2000
//! C Example               :ref:`c_example_modify_orbits`
//! Python Example          Migration.ipynb
//!                         EccAndIncDamping.ipynb
//! ======================= ===============================================
//! ```
//!
//! This applies physical forces that orbit-average to give exponential
//! growth/decay of the semimajor axis, eccentricity and inclination.
//! The eccentricity damping keeps the angular momentum constant
//! (corresponding to `p=1` in modify_orbits_direct), which means that
//! eccentricity damping will induce some semimajor axis evolution.
//! Additionally, eccentricity/inclination damping will induce
//! pericenter/nodal precession. Both these effects are physical, and the
//! method is more robust for strongly perturbed systems.
//!
//! **Effect Parameters**
//!
//! If coordinates not, defaults to using Jacobi coordinates.
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
//!                                          See the examples for usage.
//! ============================ =========== ==================================================================
//! ```
//!
//! **Particle Parameters**
//!
//! One can pick and choose which particles have which parameters set.
//! For each particle, any unset parameter is ignored.
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! tau_a (double)               No          Semimajor axis exponential growth/damping timescale
//! tau_e (double)               No          Eccentricity exponential growth/damping timescale
//! tau_inc (double)             No          Inclination axis exponential growth/damping timescale
//! ============================ =========== ==================================================================
//! ```
//!
//! The force additionally reads the two inner-disk-edge ("planet trap")
//! parameters off its own parameter list, exactly as the C does:
//!
//! ```text
//! ============================ =========== ==================================================================
//! ide_position (double)        No          Position of the inner disk edge (planet trap). Only used if
//!                                          ide_width is also set.
//! ide_width (double)           No          Width of the inner disk edge. Only used if ide_position is
//!                                          also set.
//! ============================ =========== ==================================================================
//! ```

use crate::core::rebx_get_param_double;
use crate::core::rebx_get_param_int;
use crate::inner_disk_edge::rebx_calculate_planet_trap;
use crate::rebxtools::{rebx_com_force, REBX_COORDINATES, REBX_COORDINATES_JACOBI};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_orbit_from_particle_err, reb_particle, reb_simulation, reb_vec3d};

/// modify_orbits_forces.c `rebx_calculate_modify_orbits_forces`
/// (C: `static`).
///
/// The C takes `struct reb_particle* p` and `struct reb_particle* source`
/// and reads `p->ap` for the per-particle timescales; here the particle's
/// parameter list is named by `p_index` (C: `&p->ap`) and the force's own
/// list by `force_idx` (C: `force->ap`). Both particles are passed by
/// value because the C only reads them.
fn rebx_calculate_modify_orbits_forces(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    p_index: usize,
    p: reb_particle,
    source: reb_particle,
) -> reb_vec3d {
    let mut invtau_a = 0.0;
    let mut tau_e = f64::INFINITY;
    let mut tau_inc = f64::INFINITY;

    let tau_a_ptr = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_a");
    let tau_e_ptr = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_e");
    let tau_inc_ptr = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_inc");

    //Implement the planet trap
    let dedge = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ide_position");
    let hedge = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ide_width");

    let dvx = p.vx - source.vx;
    let dvy = p.vy - source.vy;
    let dvz = p.vz - source.vz;
    let dx = p.x - source.x;
    let dy = p.y - source.y;
    let dz = p.z - source.z;
    let r2 = dx * dx + dy * dy + dz * dz;

    if let Some(tau_a_ptr) = tau_a_ptr {
        invtau_a = 1.0 / (tau_a_ptr);
        // C: if ((dedge!=NULL)&(hedge!=NULL)) — a bitwise & of the two
        // NULL tests, i.e. both must be present (no short circuit, but
        // neither side has a side effect).
        if let (Some(dedge), Some(hedge)) = (dedge, hedge) {
            let mut err = 0;
            let o = reb_orbit_from_particle_err(sim.G, p, source, &mut err);
            let a0 = o.a;
            invtau_a *= rebx_calculate_planet_trap(a0, dedge, hedge);
        }
    }
    if let Some(tau_e_ptr) = tau_e_ptr {
        tau_e = tau_e_ptr;
    }
    if let Some(tau_inc_ptr) = tau_inc_ptr {
        tau_inc = tau_inc_ptr;
    }

    let mut a = reb_vec3d { x: 0., y: 0., z: 0. }; // C: struct reb_vec3d a = {0};

    a.x = dvx * invtau_a / (2.);
    a.y = dvy * invtau_a / (2.);
    a.z = dvz * invtau_a / (2.);

    if tau_e < f64::INFINITY || tau_inc < f64::INFINITY {
        let vdotr = dx * dvx + dy * dvy + dz * dvz;
        let prefac = 2. * vdotr / r2 / tau_e;
        a.x += prefac * dx;
        a.y += prefac * dy;
        a.z += prefac * dz + 2. * dvz / tau_inc;
    }
    a
}

/// modify_orbits_forces.c `rebx_modify_orbits_forces` — the force's
/// `update_accelerations` entry point.
pub fn rebx_modify_orbits_forces(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let ptr = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "coordinates");
    let mut coordinates: REBX_COORDINATES = REBX_COORDINATES_JACOBI; // Default
    if let Some(ptr) = ptr {
        coordinates = ptr;
    }
    let back_reactions_inclusive = 1;
    let reference_name = "primary";
    rebx_com_force(
        sim,
        rebx,
        force_idx,
        coordinates,
        back_reactions_inclusive,
        reference_name,
        rebx_calculate_modify_orbits_forces,
        N,
    );
}
