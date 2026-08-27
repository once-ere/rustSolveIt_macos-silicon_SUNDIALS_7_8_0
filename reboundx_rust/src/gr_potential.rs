//! gr_potential.rs — translation of REBOUNDx gr_potential.c
//! Post-newtonian general relativity corrections using a simple potential that
//! gets the pericenter precession right.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Authors: Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo.
//! Copyright (c) 2015 Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo.
//!
//! # General Relativity
//!
//! ======================= ===============================================
//! Authors                 H. Rein, D. Tamayo
//! Implementation Paper    Tamayo, Rein, Shi and Hernandez, 2019
//!                         <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>
//! Based on                Nobili and Roxburgh 1986
//!                         <http://labs.adsabs.harvard.edu/adsabs/abs/1986IAUS..114..105N/>
//! C Example               `c_example_gr`
//! Python Example          `GeneralRelativity.ipynb`
//! ======================= ===============================================
//!
//! This is the simplest potential you can use for general relativity.
//! It assumes that the masses are dominated by a single central body.
//! It gets the precession right, but gets the mean motion wrong by
//! O(GM/ac^2). It's the fastest option, and because it's not
//! velocity-dependent, it automatically keeps WHFast symplectic.
//! Nice if you have a single-star system, don't need to get GR exactly right,
//! and want speed.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! c (double)                   Yes         Speed of light, needs to be specified in the units used for the simulation.
//! ============================ =========== ==================================================================

use rebound_rs::{reb_particle, reb_simulation, reb_simulation_error};

use crate::core::{rebx_error, rebx_get_param_double};
use crate::types::{rebx_ap, rebx_extras};

/// C: `static void rebx_calculate_gr_potential(struct reb_particle* const
/// particles, const int N, const double C2, const double G)`.
fn rebx_calculate_gr_potential(particles: &mut [reb_particle], N: usize, C2: f64, G: f64) {
    let source = particles[0];
    let prefac1 = 6. * (G * source.m) * (G * source.m) / C2;
    for i in 1..N {
        let p = particles[i];
        let dx = p.x - source.x;
        let dy = p.y - source.y;
        let dz = p.z - source.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        let prefac = prefac1 / (r2 * r2);

        particles[i].ax -= prefac * dx;
        particles[i].ay -= prefac * dy;
        particles[i].az -= prefac * dz;
        particles[0].ax += p.m / source.m * prefac * dx;
        particles[0].ay += p.m / source.m * prefac * dy;
        particles[0].az += p.m / source.m * prefac * dz;
    }
}

/// C: `void rebx_gr_potential(struct reb_simulation* const sim,
/// struct rebx_force* const gr_potential, struct reb_particle* const
/// particles, const int N)`.
pub fn rebx_gr_potential(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "c");
    match c {
        None => {
            reb_simulation_error(sim, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
        }
        Some(c) => {
            let C2 = c * c;
            let G = sim.G;
            rebx_calculate_gr_potential(&mut sim.particles, N, C2, G);
        }
    }
}

/// C: `static double rebx_calculate_gr_potential_potential(
/// struct reb_simulation* const sim, const double C2)`.
fn rebx_calculate_gr_potential_potential(sim: &reb_simulation, C2: f64) -> f64 {
    let particles = &sim.particles;
    let N = sim.N;
    let G = sim.G;
    let source = particles[0];
    let mu = G * source.m;
    let prefac = 3. * mu * mu / C2;
    let mut H = 0.;

    for i in 1..N {
        let pi = particles[i];
        let dx = pi.x - source.x;
        let dy = pi.y - source.y;
        let dz = pi.z - source.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        H -= prefac * pi.m / r2;
    }

    H
}

/// C: `double rebx_gr_potential_potential(struct rebx_extras* const rebx,
/// const struct rebx_force* const gr_potential)`.
///
/// The C reaches the simulation through `rebx->sim`; here `sim` is passed
/// explicitly, since `rebx_extras` holds no back-pointer.
pub fn rebx_gr_potential_potential(
    sim: &reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
) -> f64 {
    let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "c");
    let c = match c {
        None => {
            rebx_error(rebx, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
            return 0.;
        }
        Some(c) => c,
    };
    let C2 = c * c;
    rebx_calculate_gr_potential_potential(sim, C2)
}
