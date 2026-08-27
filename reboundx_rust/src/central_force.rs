//! central_force.rs — translation of REBOUNDx central_force.c
//! A general central force: an acceleration a = Acentral*r^gammacentral,
//! outward along the direction from a central particle to the body.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Central Force
//!
//! | | |
//! | --- | --- |
//! | Authors | D. Tamayo |
//! | Implementation Paper | Tamayo, Rein, Shi and Hernandez, 2019 |
//! | Based on | None |
//! | C Example | `c_example_central_force` |
//! | Python Example | `CentralForce.ipynb` |
//!
//! Adds a general central acceleration of the form
//! a = Acentral*r^gammacentral, outward along the direction from a
//! central particle to the body. Effect is turned on by adding Acentral
//! and gammacentral parameters to a particle, which will act as the
//! central body for the effect, and will act on all other particles.
//!
//! ## Effect Parameters
//!
//! None
//!
//! ## Particle Parameters
//!
//! | Field (C type) | Required | Description |
//! | --- | --- | --- |
//! | Acentral (double) | Yes | Normalization for central acceleration. |
//! | gammacentral (double) | Yes | Power index for central acceleration. |

use rebound_rs::{reb_orbit_from_particle, reb_particle, reb_simulation, reb_simulation_error};

use crate::core::rebx_get_param_double;
use crate::types::{rebx_ap, rebx_extras};

/// central_force.c `rebx_calculate_central_force`.
///
/// Deviation: the C takes `struct reb_simulation* const sim` as its
/// first argument but never reads it (the particle array is passed
/// separately); it is dropped here because `particles` is borrowed out
/// of the simulation and the two borrows cannot coexist in safe Rust.
fn rebx_calculate_central_force(
    particles: &mut [reb_particle],
    N: usize,
    A: f64,
    gamma: f64,
    source_index: usize,
) {
    // The C copies the source particle *before* the loop; the copy is
    // therefore stale with respect to the accelerations accumulated
    // into particles[source_index] below. Keep it that way.
    let source = particles[source_index];
    for i in 0..N {
        if i == source_index {
            continue;
        }
        let p = particles[i];
        let dx = p.x - source.x;
        let dy = p.y - source.y;
        let dz = p.z - source.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        let prefac = A * r2.powf((gamma - 1.) / 2.);

        particles[i].ax += prefac * dx;
        particles[i].ay += prefac * dy;
        particles[i].az += prefac * dz;
        particles[source_index].ax -= p.m / source.m * prefac * dx;
        particles[source_index].ay -= p.m / source.m * prefac * dy;
        particles[source_index].az -= p.m / source.m * prefac * dz;
    }
}

/// central_force.c `rebx_central_force` — the force's
/// `update_accelerations` callback.
pub fn rebx_central_force(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    N: usize,
) {
    for i in 0..N {
        let Acentral = rebx_get_param_double(rebx, rebx_ap::particle(i), "Acentral");
        if let Some(Acentral) = Acentral {
            let gammacentral = rebx_get_param_double(rebx, rebx_ap::particle(i), "gammacentral");
            if let Some(gammacentral) = gammacentral {
                // only calculates force if a particle has both Acentral
                // and gammacentral parameters set.
                rebx_calculate_central_force(&mut sim.particles, N, Acentral, gammacentral, i);
            }
        }
    }
}

/// central_force.c `rebx_calculate_central_force_potential`.
fn rebx_calculate_central_force_potential(
    sim: &reb_simulation,
    A: f64,
    gamma: f64,
    source_index: usize,
) -> f64 {
    let particles = &sim.particles;
    let N = sim.N;
    let source = particles[source_index];
    let mut H = 0.;
    for i in 0..N {
        if i == source_index {
            continue;
        }
        let p = particles[i];
        let dx = p.x - source.x;
        let dy = p.y - source.y;
        let dz = p.z - source.z;
        let r2 = dx * dx + dy * dy + dz * dz;

        if (gamma + 1.).abs() < f64::EPSILON {
            // F propto 1/r
            H -= p.m * A * r2.sqrt().ln();
        } else {
            H -= p.m * A * r2.powf((gamma + 1.) / 2.) / (gamma + 1.);
        }
    }
    H
}

/// central_force.c `rebx_central_force_potential`.
///
/// Deviation: the C takes only `struct rebx_extras* rebx` and reaches
/// the simulation through `rebx->sim`. `rebx_extras` holds no
/// back-pointer here, so the simulation is passed explicitly and the C's
/// `if (rebx->sim == NULL) { rebx_error(rebx, ""); return 0; }` guard
/// has no reachable counterpart — a caller that can name both arguments
/// necessarily has them both attached.
pub fn rebx_central_force_potential(sim: &reb_simulation, rebx: &rebx_extras) -> f64 {
    // The C also caches `struct reb_particle* const particles =
    // sim->particles` here, only to reach `particles[i].ap`; the
    // parameter lists live in `rebx` in this translation, so the cache
    // has nothing left to do.
    let N = sim.N;
    let mut Htot = 0.;
    for i in 0..N {
        let Acentral = rebx_get_param_double(rebx, rebx_ap::particle(i), "Acentral");
        if let Some(Acentral) = Acentral {
            let gammacentral = rebx_get_param_double(rebx, rebx_ap::particle(i), "gammacentral");
            if let Some(gammacentral) = gammacentral {
                Htot += rebx_calculate_central_force_potential(sim, Acentral, gammacentral, i);
            }
        }
    }
    Htot
}

/// central_force.c `rebx_central_force_Acentral` — initialize Acentral
/// from a desired pericenter precession rate `pomegadot`.
///
/// Deviation: the C reads the simulation through the particle's `sim`
/// back-pointer (`p.sim`), which `reb_particle` does not carry in this
/// translation; the simulation is passed explicitly instead. It is
/// `&mut` because the gamma = -2 branch calls `reb_simulation_error`.
pub fn rebx_central_force_Acentral(
    sim: &mut reb_simulation,
    p: reb_particle,
    primary: reb_particle,
    pomegadot: f64,
    gamma: f64,
) -> f64 {
    let G = sim.G;
    let o = reb_orbit_from_particle(G, p, primary);
    if (gamma + 2.).abs() < f64::EPSILON {
        // precession goes to 0 at r^-2, so A diverges for gamma=-2
        reb_simulation_error(sim, "Precession vanishes for force law varying as r^-2, so can't initialize Acentral from a precession rate for gamma=-2)\n");
        return 0.;
    }
    G * primary.m * pomegadot / (1. + gamma / 2.) / o.d.powf(gamma + 2.) / o.n
}
