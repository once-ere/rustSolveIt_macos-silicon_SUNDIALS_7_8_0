//! lense_thirring.rs — translation of REBOUNDx lense_thirring.c
//! Adds the Lense-Thirring effect due to a rotating central body in the
//! simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Author: Arya Akmal <akmala@gmail.com>.
//! Copyright (c) 2023 Arya Akmal.
//!
//! # General Relativity
//!
//! ======================= ===============================================
//! Authors                 A. Akmal
//! Implementation Paper    None
//! Based on                `Park et al.
//!                         <https://iopscience.iop.org/article/10.3847/1538-3881/abd414/>`_.
//! C Example               `c_example_lense_thirring`
//! Python Example          `LenseThirring.ipynb`
//! ======================= ===============================================
//!
//! Adds Lense-Thirring effect due to rotating central body in the simulation.
//! Assumes the source body is particles[0]
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! lt_c (double)                Yes         Speed of light in the units used for the simulation.
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! I (double)                   Yes         Moment of Inertia of source body
//! Omega (reb_vec3d)            Yes         Angular rotation frequency (Omega_x, Omega_y, Omega_z)
//! ============================ =========== ==================================================================

use rebound_rs::{reb_particle, reb_simulation, reb_simulation_error, reb_vec3d};

use crate::core::{rebx_get_param_double, rebx_get_param_vec3d};
use crate::types::{rebx_ap, rebx_extras};

/// C: `static void rebx_calculate_LT_force(struct reb_simulation* const sim,
/// struct reb_particle* const particles, const int N,
/// const struct reb_vec3d Omega, const double I, const double C2)`.
///
/// The C reads `sim->G` as its first statement and never touches `sim`
/// again; `G` is therefore passed in directly, exactly as the C's
/// `const double G = sim->G;` would have bound it.
fn rebx_calculate_LT_force(
    G: f64,
    particles: &mut [reb_particle],
    N: usize,
    Omega: reb_vec3d,
    I: f64,
    C2: f64,
) {
    let gamma = 1.000021; // hard-coded Eddington-Robertson-Shiff parameter for now
    let source = particles[0]; // hard-code particles[0] as source particle
    for i in 1..N {
        let p = particles[i];
        let dx = p.x - source.x;
        let dy = p.y - source.y;
        let dz = p.z - source.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        let r = r2.sqrt();
        let r3 = r2 * r;
        let dvx = p.vx - source.vx;
        let dvy = p.vy - source.vy;
        let dvz = p.vz - source.vz;
        let Jx = I * Omega.x; //C_fac*source.m * R_eq*R_eq*omega*p_hat_x ;
        let Jy = I * Omega.y; //C_fac*source.m * R_eq*R_eq*omega*p_hat_y ;
        let Jz = I * Omega.z; //C_fac*source.m * R_eq*R_eq*omega*p_hat_z ;
        let ms = source.m;
        let mt = p.m;
        let mtot = ms + mt;
        let mratio = mt / ms;
        let Omega_fac = ms / mtot * (1. + gamma) * G / 2. / C2 / r3;
        let Jdotr = Jx * dx + Jy * dy + Jz * dz;
        let Omega_x = Omega_fac * (-Jx + 3. * Jdotr * dx / r2);
        let Omega_y = Omega_fac * (-Jy + 3. * Jdotr * dy / r2);
        let Omega_z = Omega_fac * (-Jz + 3. * Jdotr * dz / r2);

        particles[i].ax += 2. * (Omega_y * dvz - Omega_z * dvy);
        particles[i].ay += 2. * (Omega_z * dvx - Omega_x * dvz);
        particles[i].az += 2. * (Omega_x * dvy - Omega_y * dvx);
        particles[0].ax -= mratio * 2. * (Omega_y * dvz - Omega_z * dvy);
        particles[0].ay -= mratio * 2. * (Omega_z * dvx - Omega_x * dvz);
        particles[0].az -= mratio * 2. * (Omega_x * dvy - Omega_y * dvx);
    }
}

/// C: `void rebx_lense_thirring(struct reb_simulation* const sim,
/// struct rebx_force* const force, struct reb_particle* const particles,
/// const int N)`.
///
/// Deviation: after reporting the missing `lt_c` the C falls through and
/// dereferences the NULL pointer (`(*c)*(*c)`), i.e. it crashes. Here the
/// error is reported and the force contributes nothing this step; there is
/// no defined C behavior to reproduce.
pub fn rebx_lense_thirring(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "lt_c");
    let c = match c {
        None => {
            reb_simulation_error(sim, "REBOUNDx Error: Need to set speed of light in LT effect.  See examples in documentation.\n");
            return;
        }
        Some(c) => c,
    };
    let C2 = c * c;

    let I = rebx_get_param_double(rebx, rebx_ap::particle(0), "I");
    if let Some(I) = I {
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(0), "Omega");
        if let Some(Omega) = Omega {
            let G = sim.G;
            rebx_calculate_LT_force(G, &mut sim.particles, N, Omega, I, C2);
        }
    }
}
