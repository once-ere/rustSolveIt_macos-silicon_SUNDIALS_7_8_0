//! exponential_migration.rs — translation of REBOUNDx exponential_migration.c
//! Continuous velocity kicks leading to exponential change in the object's semimajor axis.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! ```text
//! @file    exponential_migration.c
//! @brief   Continuous velocity kicks leading to exponential change in the object's semimajor axis.
//! @author  Mohamad Ali-Dib <mma9132@nyu.edu>
//!
//! @section     LICENSE
//! Copyright (c) 2021 Mohamad Ali-Dib
//! ```
//!
//! $Orbit Modifications$
//!
//! ```text
//! ======================= ===============================================
//! Author                  Mohamad Ali-Dib
//! Implementation Paper    `Ali-Dib et al., 2021 AJ <https://arxiv.org/abs/2104.04271>`_.
//! Based on                `Hahn & Malhotra 2005 <https://ui.adsabs.harvard.edu/abs/2005AJ....130.2392H/abstract>`_.
//! C Example               :ref:`c_example_exponential_migration`
//! Python Example          `ExponentialMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ExponentialMigration.ipynb>`_.
//! ======================= ===============================================
//! ```
//!
//! Continuous velocity kicks leading to exponential change in the object's semimajor axis.
//! One of the standard prescriptions often used in Neptune migration & Kuiper Belt formation models.
//! Does not directly affect the eccentricity or inclination of the object.
//!
//! **Particle Parameters**
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! em_tau_a (double)              Yes          Semimajor axis exponential growth/damping timescale
//! em_aini (double)               Yes          Object's initial semimajor axis
//! em_afin (double)               Yes          Object's final semimajor axis
//! ============================ =========== ==================================================================
//! ```
//!
//! # Deviations from the C, all mechanical
//!
//! * The C's `struct rebx_force*` is an index into
//!   `rebx_extras::allocated_forces`, and the REBOUNDx state travels as an
//!   explicit `&mut rebx_extras` rather than through `sim->extras`.
//! * The per-particle callback reaches `p->ap` through the particle index
//!   `p_index` that [`rebx_com_force`] hands it (see rebxtools).

use crate::core::{rebx_get_param_double, rebx_get_param_int};
use crate::rebxtools::{rebx_com_force, REBX_COORDINATES, REBX_COORDINATES_JACOBI};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_orbit_from_particle, reb_particle, reb_simulation, reb_vec3d};

/// exponential_migration.c `rebx_calculate_modify_orbits_forces_new`
/// (C: `static`). Per-particle acceleration callback handed to
/// [`rebx_com_force`].
fn rebx_calculate_modify_orbits_forces_new(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    p_index: usize,
    p: reb_particle,
    source: reb_particle,
) -> reb_vec3d {
    let o = reb_orbit_from_particle(sim.G, p, source);

    let mut em_tau_a = f64::INFINITY;
    let mut em_aini = 24.;
    let mut em_afin = 30.;

    let em_tau_a_ptr = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "em_tau_a");
    let em_ainipoint = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "em_aini");
    let em_afinpoint = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "em_afin");

    let dvx = p.vx - source.vx;
    let dvy = p.vy - source.vy;
    let dvz = p.vz - source.vz;
    //const double dx = p->x-source->x;
    //const double dy = p->y-source->y;
    //const double dz = p->z-source->z;
    //const double r2 = dx*dx + dy*dy + dz*dz;

    if let Some(v) = em_tau_a_ptr {
        em_tau_a = v;
    }
    if let Some(v) = em_ainipoint {
        em_aini = v;
    }
    if let Some(v) = em_afinpoint {
        em_afin = v;
    }

    let mut a = reb_vec3d { x: 0., y: 0., z: 0. }; // C: struct reb_vec3d a = {0};

    a.x = (dvx / (2. * em_tau_a)) * ((em_afin - em_aini) / (o.a)) * (-(sim.t) / em_tau_a).exp();
    a.y = (dvy / (2. * em_tau_a)) * ((em_afin - em_aini) / (o.a)) * (-(sim.t) / em_tau_a).exp();
    a.z = (dvz / (2. * em_tau_a)) * ((em_afin - em_aini) / (o.a)) * (-(sim.t) / em_tau_a).exp();

    a
}

/// exponential_migration.c `rebx_exponential_migration`.
///
/// Reads the effect's `coordinates` parameter (default
/// `REBX_COORDINATES_JACOBI`) and applies the kicks through
/// [`rebx_com_force`] with inclusive back reactions relative to the
/// `"primary"` reference particle.
pub fn rebx_exponential_migration(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let ptr = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "coordinates");
    let mut coordinates: REBX_COORDINATES = REBX_COORDINATES_JACOBI; // Default
    if let Some(v) = ptr {
        coordinates = v;
    }
    let back_reactions_inclusive: i32 = 1;
    let reference_name = "primary";
    rebx_com_force(
        sim,
        rebx,
        force_idx,
        coordinates,
        back_reactions_inclusive,
        reference_name,
        rebx_calculate_modify_orbits_forces_new,
        N,
    );
}
