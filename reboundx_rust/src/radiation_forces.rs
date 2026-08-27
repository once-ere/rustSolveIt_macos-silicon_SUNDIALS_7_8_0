//! radiation_forces.rs — translation of REBOUNDx radiation_forces.c
//! Adds radiation forces: both radiation pressure and Poynting-Robertson drag.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Radiation Forces
//!
//! ======================= ===============================================
//! Authors                 H. Rein, D. Tamayo
//! Implementation Paper    Tamayo, Rein, Shi and Hernandez, 2019
//!                         <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>
//! Based on                Burns et al. 1979
//!                         <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>
//! C Example               `c_example_rad_forces_debris_disk`,
//!                         `c_example_rad_forces_circumplanetary`.
//! Python Example          `Radiation_Forces_Debris_Disk.ipynb`,
//!                         `Radiation_Forces_Circumplanetary_Dust.ipynb`.
//! ======================= ===============================================
//!
//! This applies radiation forces to particles in the simulation.
//! It incorporates both radiation pressure and Poynting-Robertson drag.
//! Only particles whose `beta` parameter is set will feel the radiation.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! c (double)                   Yes         Speed of light in the units used for the simulation.
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! If no particles have radiation_source set, effect will assume the particle
//! at index 0 in the particles array is the source.
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! radiation_source (int)       No          Flag identifying the particle as the source of radiation.
//! beta (float)                 Yes         Ratio of radiation pressure force to gravitational force.
//!                                          Particles without beta set feel no radiation forces.
//! ============================ =========== ==================================================================

use crate::core::{rebx_get_param_double, rebx_get_param_int};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_simulation, reb_simulation_error};

/// radiation_forces.c `rebx_calculate_radiation_forces`.
///
/// The C takes `particles` as a separate argument; here it is
/// `sim.particles`, mutated in place (see the crate's translation notes).
fn rebx_calculate_radiation_forces(
    rebx: &rebx_extras,
    sim: &mut reb_simulation,
    c: f64,
    source_index: usize,
    N: usize,
) {
    let source = sim.particles[source_index];
    let mu = sim.G * source.m;

    for i in 0..N {
        if i == source_index {
            continue;
        }

        // only particles with beta set feel radiation forces
        let beta = match rebx_get_param_double(rebx, rebx_ap::particle(i), "beta") {
            Some(beta) => beta,
            None => continue,
        };

        let p = sim.particles[i];
        let dx = p.x - source.x;
        let dy = p.y - source.y;
        let dz = p.z - source.z;
        let dr = (dx * dx + dy * dy + dz * dz).sqrt(); // distance to star

        let dvx = p.vx - source.vx;
        let dvy = p.vy - source.vy;
        let dvz = p.vz - source.vz;
        let rdot = (dx * dvx + dy * dvy + dz * dvz) / dr; // radial velocity
        let a_rad = beta * mu / (dr * dr);

        // Equation (5) of Burns, Lamy & Soter (1979)

        sim.particles[i].ax += a_rad * ((1. - rdot / c) * dx / dr - dvx / c);
        sim.particles[i].ay += a_rad * ((1. - rdot / c) * dy / dr - dvy / c);
        sim.particles[i].az += a_rad * ((1. - rdot / c) * dz / dr - dvz / c);
    }
}

/// radiation_forces.c `rebx_radiation_forces` — the force's
/// `update_accelerations` entry point.
pub fn rebx_radiation_forces(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let c = match rebx_get_param_double(rebx, rebx_ap::force(force_idx), "c") {
        Some(c) => c,
        None => {
            reb_simulation_error(
                sim,
                "Need to set speed of light in radiation_forces effect.  See examples in documentation.\n",
            );
            return;
        }
    };

    let mut source_found = 0;
    for i in 0..N {
        if rebx_get_param_int(rebx, rebx_ap::particle(i), "radiation_source").is_some() {
            source_found = 1;
            rebx_calculate_radiation_forces(rebx, sim, c, i, N);
        }
    }
    if source_found == 0 {
        // default source to index 0 if "radiation_source" not found on any particle
        rebx_calculate_radiation_forces(rebx, sim, c, 0, N);
    }
}

/// radiation_forces.c `rebx_rad_calc_beta`. Returns the ratio of the
/// radiation pressure force to the gravitational force for a grain of the
/// given radius, density and radiation pressure coefficient `Q_pr`.
pub fn rebx_rad_calc_beta(
    G: f64,
    c: f64,
    source_mass: f64,
    source_luminosity: f64,
    radius: f64,
    density: f64,
    Q_pr: f64,
) -> f64 {
    3. * source_luminosity * Q_pr
        / (16. * std::f64::consts::PI * G * source_mass * c * density * radius)
}

/// radiation_forces.c `rebx_rad_calc_particle_radius`. Inverse of
/// `rebx_rad_calc_beta`: the grain radius that yields the given `beta`.
pub fn rebx_rad_calc_particle_radius(
    G: f64,
    c: f64,
    source_mass: f64,
    source_luminosity: f64,
    beta: f64,
    density: f64,
    Q_pr: f64,
) -> f64 {
    3. * source_luminosity * Q_pr
        / (16. * std::f64::consts::PI * G * source_mass * c * density * beta)
}
