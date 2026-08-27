//! stochastic_forces.rs — translation of REBOUNDx stochastic_forces.c
//! Add stochastic forces to particles in the simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Author: Hanno Rein <hanno@hanno-rein.de>.
//! Copyright (c) 2022 Hanno Rein, Dan Tamayo.
//!
//! # Stochastic Forces
//!
//! ======================= ===============================================
//! Authors                 H. Rein
//! Based on                Rein and Papaloizou 2009
//!                         <https://ui.adsabs.harvard.edu/abs/2009A%26A...497..595R/abstract>
//! Implementation Paper    Rein and Choksi 2022
//!                         <https://iopscience.iop.org/article/10.3847/2515-5172/ac6e41>
//! C Example               `c_example_stochastic_forces`
//! Python Example          `StochasticForces.ipynb`, `StochasticForcesCartesian.ipynb`
//! ======================= ===============================================
//!
//! This applies stochastic forces to particles in the simulation.
//!
//! **Effect Parameters**
//!
//! None
//!
//! **Particle Parameters**
//!
//! All particles which have the field kappa set, will experience stochastic forces.
//! The particle with index 0 cannot experience stochastic forces.
//!
//! ============================ =========== ==================================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================================
//! kappa (double)               Yes         Strength of stochastic forces relative to gravity from central object
//! tau_kappa (double)           No          Auto-correlation time of stochastic forces. Defaults to orbital period if not set.
//!                                          The units are relative to the current orbital period.
//! ============================ =========== ==================================================================================
//!
//! The Cartesian variant additionally uses `kappa_x`/`kappa_y`/`kappa_z`
//! (each optional, but each requiring the matching `tau_kappa_x`,
//! `tau_kappa_y` or `tau_kappa_z` to be set). Unlike `kappa`, these apply
//! to every particle, including index 0.
//!
//! ## RNG note
//!
//! The C helper `rebx_random_normal2` normalises `rand_r` by the C
//! library's `RAND_MAX`. REBOUND vendors glibc's 31-bit `rand_r` on
//! Windows but `<stdlib.h>`'s `RAND_MAX` there is 32767, which makes the
//! rejection loop in the C never terminate in practice. This translation
//! uses REBOUND's own `REB_RAND_MAX` (2147483647, the value `RAND_MAX`
//! takes on Linux and macOS, and the value REBOUND's `reb_random_normal`
//! uses on every platform), so the random stream matches the C exactly
//! wherever the C is usable.

use rebound_rs::{
    rand_r, reb_orbit_from_particle_err, reb_particle_com_of_pair, reb_simulation,
    reb_simulation_error, REB_RAND_MAX,
};

use crate::core::{rebx_get_param_double, rebx_set_param_double};
use crate::types::{rebx_ap, rebx_extras};

/// C: `static void rebx_random_normal2(struct reb_simulation* r,
/// double* n0, double* n1)`.
///
/// Marsaglia polar method, keeping *both* variates (unlike REBOUND's
/// `reb_random_normal`, which discards the second). Returns `(n0, n1)`.
fn rebx_random_normal2(r: &mut reb_simulation) -> (f64, f64) {
    let mut v1 = 0.;
    let mut v2 = 0.;
    let mut rsq = 1.;
    while rsq >= 1. || rsq < 1.0e-12 {
        v1 = 2. * (rand_r(&mut r.rand_seed) as f64) / (REB_RAND_MAX as f64) - 1.0;
        v2 = 2. * (rand_r(&mut r.rand_seed) as f64) / (REB_RAND_MAX as f64) - 1.0;
        rsq = v1 * v1 + v2 * v2;
    }
    let n0 = v1 * (-2. * rsq.ln() / rsq).sqrt();
    let n1 = v2 * (-2. * rsq.ln() / rsq).sqrt();
    (n0, n1)
}

/// C: `void rebx_stochastic_forces(struct reb_simulation* const sim,
/// struct rebx_force* const radiation_forces,
/// struct reb_particle* const particles, const int N)`.
pub fn rebx_stochastic_forces(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    N: usize,
) {
    let mut com = sim.particles[0];

    for i in 0..N {
        let kappa = rebx_get_param_double(rebx, rebx_ap::particle(i), "kappa");
        if i > 0 && kappa.is_some() {
            let kappa = kappa.unwrap();

            // First run? The C sets the parameter to 0. and re-reads it.
            let mut stochastic_force_r =
                match rebx_get_param_double(rebx, rebx_ap::particle(i), "stochastic_force_r") {
                    Some(v) => v,
                    None => {
                        rebx_set_param_double(rebx, rebx_ap::particle(i), "stochastic_force_r", 0.);
                        0.
                    }
                };
            let mut stochastic_force_phi =
                match rebx_get_param_double(rebx, rebx_ap::particle(i), "stochastic_force_phi") {
                    Some(v) => v,
                    None => {
                        rebx_set_param_double(
                            rebx,
                            rebx_ap::particle(i),
                            "stochastic_force_phi",
                            0.,
                        );
                        0.
                    }
                };

            let p = sim.particles[i];

            // Get auto-correlation time
            let mut err: i32 = 0;
            let o = reb_orbit_from_particle_err(sim.G, sim.particles[i], com, &mut err);
            if err != 0 {
                reb_simulation_error(
                    sim,
                    "An error occured during the orbit calculation in rebx_stochastic_forces.\n",
                );
                return;
            }
            let mut tau = o.P; // Default is current orbital period.

            let tau_kappa = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau_kappa");
            if let Some(tau_kappa) = tau_kappa {
                tau *= tau_kappa;
            }

            let dt = sim.dt_last_done;

            let prefac = (-dt / tau).exp();

            // Decay
            stochastic_force_r = stochastic_force_r * prefac;
            stochastic_force_phi = stochastic_force_phi * prefac;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_r",
                stochastic_force_r,
            );
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_phi",
                stochastic_force_phi,
            );

            let variance = 1. - prefac * prefac;
            if variance < 0. {
                reb_simulation_error(
                    sim,
                    "Timestep is larger than the correlation time for stochastic forces.\n",
                );
                return;
            }
            let std = variance.sqrt();

            let (n0, n1) = rebx_random_normal2(sim);

            // Excitation
            stochastic_force_r = stochastic_force_r + n0 * std;
            stochastic_force_phi = stochastic_force_phi + n1 * std;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_r",
                stochastic_force_r,
            );
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_phi",
                stochastic_force_phi,
            );

            let dx = p.x - com.x;
            let dy = p.y - com.y;
            let dz = p.z - com.z;
            let dr = (dx * dx + dy * dy + dz * dz).sqrt();

            let dvx = p.vx - com.vx;
            let dvy = p.vy - com.vy;
            let dvz = p.vz - com.vz;
            let dv = (dvx * dvx + dvy * dvy + dvz * dvz).sqrt();

            let force_prefac = kappa * sim.G / (dr * dr) * com.m;
            sim.particles[i].ax +=
                force_prefac * (stochastic_force_r * dx / dr + stochastic_force_phi * dvx / dv);
            sim.particles[i].ay +=
                force_prefac * (stochastic_force_r * dy / dr + stochastic_force_phi * dvy / dv);
            sim.particles[i].az +=
                force_prefac * (stochastic_force_r * dz / dr + stochastic_force_phi * dvz / dv);

            com = reb_particle_com_of_pair(com, p);
        }

        let kappa_x = rebx_get_param_double(rebx, rebx_ap::particle(i), "kappa_x");
        if let Some(kappa_x) = kappa_x {
            let mut stochastic_force_x =
                match rebx_get_param_double(rebx, rebx_ap::particle(i), "stochastic_force_x") {
                    Some(v) => v,
                    None => {
                        // First run?
                        rebx_set_param_double(rebx, rebx_ap::particle(i), "stochastic_force_x", 0.);
                        0.
                    }
                };

            let tau_kappa_x = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau_kappa_x");
            let tau_kappa_x = match tau_kappa_x {
                Some(v) => v,
                None => {
                    reb_simulation_error(
                        sim,
                        "Need to set tau_kappa_x to enable stochastic forces.\n",
                    );
                    return;
                }
            };

            let dt = sim.dt_last_done;
            let prefac = (-dt / tau_kappa_x).exp();

            // Decay
            stochastic_force_x = stochastic_force_x * prefac;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_x",
                stochastic_force_x,
            );

            // Excitation
            let variance = 1. - prefac * prefac;
            if variance < 0. {
                reb_simulation_error(
                    sim,
                    "Timestep is larger than the correlation time for stochastic forces.\n",
                );
                return;
            }
            let std = kappa_x * variance.sqrt();
            let (n0, _n1) = rebx_random_normal2(sim);
            stochastic_force_x = stochastic_force_x + n0 * std;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_x",
                stochastic_force_x,
            );

            sim.particles[i].ax += stochastic_force_x;
        }

        let kappa_y = rebx_get_param_double(rebx, rebx_ap::particle(i), "kappa_y");
        if let Some(kappa_y) = kappa_y {
            let mut stochastic_force_y =
                match rebx_get_param_double(rebx, rebx_ap::particle(i), "stochastic_force_y") {
                    Some(v) => v,
                    None => {
                        // First run?
                        rebx_set_param_double(rebx, rebx_ap::particle(i), "stochastic_force_y", 0.);
                        0.
                    }
                };

            let tau_kappa_y = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau_kappa_y");
            let tau_kappa_y = match tau_kappa_y {
                Some(v) => v,
                None => {
                    reb_simulation_error(
                        sim,
                        "Need to set tau_kappa_y to enable stochastic forces.\n",
                    );
                    return;
                }
            };

            let dt = sim.dt_last_done;
            let prefac = (-dt / tau_kappa_y).exp();

            // Decay
            stochastic_force_y = stochastic_force_y * prefac;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_y",
                stochastic_force_y,
            );

            // Excitation
            let variance = 1. - prefac * prefac;
            if variance < 0. {
                reb_simulation_error(
                    sim,
                    "Timestep is larger than the correlation time for stochastic forces.\n",
                );
                return;
            }
            let std = kappa_y * variance.sqrt();
            let (n0, _n1) = rebx_random_normal2(sim);
            stochastic_force_y = stochastic_force_y + n0 * std;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_y",
                stochastic_force_y,
            );

            sim.particles[i].ay += stochastic_force_y;
        }

        let kappa_z = rebx_get_param_double(rebx, rebx_ap::particle(i), "kappa_z");
        if let Some(kappa_z) = kappa_z {
            let mut stochastic_force_z =
                match rebx_get_param_double(rebx, rebx_ap::particle(i), "stochastic_force_z") {
                    Some(v) => v,
                    None => {
                        // First run?
                        rebx_set_param_double(rebx, rebx_ap::particle(i), "stochastic_force_z", 0.);
                        0.
                    }
                };

            let tau_kappa_z = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau_kappa_z");
            let tau_kappa_z = match tau_kappa_z {
                Some(v) => v,
                None => {
                    reb_simulation_error(
                        sim,
                        "Need to set tau_kappa_z to enable stochastic forces.\n",
                    );
                    return;
                }
            };

            let dt = sim.dt_last_done;
            let prefac = (-dt / tau_kappa_z).exp();

            // Decay
            stochastic_force_z = stochastic_force_z * prefac;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_z",
                stochastic_force_z,
            );

            // Excitation
            let variance = 1. - prefac * prefac;
            if variance < 0. {
                reb_simulation_error(
                    sim,
                    "Timestep is larger than the correlation time for stochastic forces.\n",
                );
                return;
            }
            let std = kappa_z * variance.sqrt();
            let (n0, _n1) = rebx_random_normal2(sim);
            stochastic_force_z = stochastic_force_z + n0 * std;
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(i),
                "stochastic_force_z",
                stochastic_force_z,
            );

            sim.particles[i].az += stochastic_force_z;
        }
    }
}
