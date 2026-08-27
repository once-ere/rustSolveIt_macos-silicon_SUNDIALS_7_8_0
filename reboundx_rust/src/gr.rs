//! gr.rs — translation of REBOUNDx gr.c
//! Post-newtonian general relativity corrections arising from a single massive body.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! @author Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein (tamayo.daniel@gmail.com)
//!
//! # General Relativity
//!
//! ======================= ===============================================
//! Authors                 P. Shi, D. Tamayo, H. Rein
//! Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
//! Based on                `Anderson et al. 1975 <http://labs.adsabs.harvard.edu/adsabs/abs/1975ApJ...200..221A/>`_.
//! C Example               :ref:`c_example_gr`
//! Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
//! ======================= ===============================================
//!
//! This assumes that the masses are dominated by a single central body, and
//! should be good enough for most applications with planets orbiting single
//! stars. It ignores terms that are smaller by of order the mass ratio with
//! the central body. It gets both the mean motion and precession correct, and
//! will be significantly faster than `gr_full`, particularly with several
//! bodies. Adding this effect to several bodies is NOT equivalent to using
//! gr_full.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! c (double)                   Yes         Speed of light, needs to be specified in the units used for the simulation.
//! ============================ =========== ==================================================================
//!
//! # Deviations from the C
//!
//! * The C's `malloc`'d scratch arrays `ps` / `ps_j` are `Vec`s. `ps_j` is
//!   value-initialised rather than left uninitialised; every field this code
//!   reads is written first by the Jacobi transformation, so the values used
//!   are identical.
//! * `reb_transformations_jacobi_to_inertial_acc(ps, ps_j, ps, ...)` aliases
//!   `ps` as both the (mutable) output array and the (read-only) mass array.
//!   Safe Rust forbids that, so a snapshot clone of `ps` is passed as the mass
//!   array. The routine only ever *reads* `.m` from it and only ever *writes*
//!   `ax`/`ay`/`az` to the output, and `gr` never changes any mass, so the
//!   snapshot is bit-identical to the aliased read.
//! * `rebx_gr_hamiltonian` takes the simulation explicitly; the C reaches it
//!   through `rebx->sim`, a back-pointer this translation does not carry.

use crate::core::{rebx_error, rebx_get_param_double, rebx_get_param_int};
use crate::rebxtools::rebx_calculate_jacobi_masses;
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{
    reb_particle, reb_simulation, reb_simulation_error, reb_simulation_warning,
    reb_transformations_inertial_to_jacobi_posvel,
    reb_transformations_inertial_to_jacobi_posvelacc, reb_transformations_jacobi_to_inertial_acc,
    reb_vec3d,
};

fn rebx_calculate_gr(
    sim: &mut reb_simulation,
    N: usize,
    C2: f64,
    G: f64,
    max_iterations: i32,
) {
    let mut ps: Vec<reb_particle> = sim.particles[0..N].to_vec();
    let mut ps_j: Vec<reb_particle> = vec![reb_particle::default(); N];

    let N_active = if sim.N_active == usize::MAX {
        N
    } else {
        sim.N_active
    };
    // Calculate Newtonian accelerations
    for i in 0..N {
        ps[i].ax = 0.;
        ps[i].ay = 0.;
        ps[i].az = 0.;
    }

    for i in 0..N_active {
        let pi = ps[i];
        for j in (i + 1)..N {
            let pj = ps[j];
            let dx = pi.x - pj.x;
            let dy = pi.y - pj.y;
            let dz = pi.z - pj.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            let prefac = G / (r2 * r);
            ps[i].ax -= prefac * pj.m * dx;
            ps[i].ay -= prefac * pj.m * dy;
            ps[i].az -= prefac * pj.m * dz;
            ps[j].ax += prefac * pi.m * dx;
            ps[j].ay += prefac * pi.m * dy;
            ps[j].az += prefac * pi.m * dz;
        }
    }

    // Transform to Jacobi coordinates
    let source = ps[0];
    let mu = G * source.m;
    reb_transformations_inertial_to_jacobi_posvelacc(&ps, &mut ps_j, &ps, N, N_active);

    for i in 1..N {
        let p = ps_j[i];
        let mut vi = reb_vec3d::default();
        vi.x = p.vx;
        vi.y = p.vy;
        vi.z = p.vz;
        let mut vi2 = vi.x * vi.x + vi.y * vi.y + vi.z * vi.z;
        let ri = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        let mut q: i32 = 0;
        let mut A = (0.5 * vi2 + 3. * mu / ri) / C2;
        let mut old_v = reb_vec3d::default();
        while q < max_iterations {
            old_v.x = vi.x;
            old_v.y = vi.y;
            old_v.z = vi.z;
            vi.x = p.vx / (1. - A);
            vi.y = p.vy / (1. - A);
            vi.z = p.vz / (1. - A);
            vi2 = vi.x * vi.x + vi.y * vi.y + vi.z * vi.z;
            A = (0.5 * vi2 + 3. * mu / ri) / C2;
            let dvx = vi.x - old_v.x;
            let dvy = vi.y - old_v.y;
            let dvz = vi.z - old_v.z;
            if (dvx * dvx + dvy * dvy + dvz * dvz) / vi2 < f64::EPSILON * f64::EPSILON {
                break;
            }
            q += 1;
        }
        let default_max_iterations: i32 = 10;
        if q == default_max_iterations {
            reb_simulation_warning(sim, "REBOUNDx Warning: 10 iterations in gr.c failed to converge. This is typically because the perturbation is too strong for the current implementation.");
        }

        let B = (mu / ri - 1.5 * vi2) * mu / (ri * ri * ri) / C2;
        let rdotrdot = p.x * p.vx + p.y * p.vy + p.z * p.vz;

        let mut vidot = reb_vec3d::default();
        vidot.x = p.ax + B * p.x;
        vidot.y = p.ay + B * p.y;
        vidot.z = p.az + B * p.z;

        let vdotvdot = vi.x * vidot.x + vi.y * vidot.y + vi.z * vidot.z;
        let D = (vdotvdot - 3. * mu / (ri * ri * ri) * rdotrdot) / C2;

        ps_j[i].ax = B * (1. - A) * p.x - A * p.ax - D * vi.x;
        ps_j[i].ay = B * (1. - A) * p.y - A * p.ay - D * vi.y;
        ps_j[i].az = B * (1. - A) * p.z - A * p.az - D * vi.z;
    }

    ps_j[0].ax = 0.;
    ps_j[0].ay = 0.;
    ps_j[0].az = 0.;

    // The C aliases `ps` as both the output and the mass array here; see the
    // module-level deviation note. Masses are read-only in the callee.
    let ps_mass = ps.clone();
    reb_transformations_jacobi_to_inertial_acc(&mut ps, &ps_j, &ps_mass, N, N_active);
    for i in 0..N {
        sim.particles[i].ax += ps[i].ax;
        sim.particles[i].ay += ps[i].ay;
        sim.particles[i].az += ps[i].az;
    }
}

/// gr.c `rebx_gr` — the force callback registered under the name `"gr"`.
pub fn rebx_gr(sim: &mut reb_simulation, rebx: &mut rebx_extras, force_idx: usize, N: usize) {
    let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "c");
    let c = match c {
        None => {
            reb_simulation_error(
                sim,
                "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n",
            );
            return;
        }
        Some(c) => c,
    };
    let C2 = c * c;
    let G = sim.G;
    let max_iterations = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "max_iterations");
    match max_iterations {
        Some(max_iterations) => {
            rebx_calculate_gr(sim, N, C2, G, max_iterations);
        }
        None => {
            let default_max_iterations: i32 = 10;
            rebx_calculate_gr(sim, N, C2, G, default_max_iterations);
        }
    }
}

fn rebx_calculate_gr_hamiltonian(
    _rebx: &rebx_extras,
    sim: &reb_simulation,
    C2: f64,
) -> f64 {
    let N = sim.N;
    let G = sim.G;

    let mut ps_j: Vec<reb_particle> = vec![reb_particle::default(); N];
    let ps: &Vec<reb_particle> = &sim.particles;
    // Calculate Newtonian potentials

    let mut V_newt = 0.;
    for i in 0..N {
        let pi = ps[i];
        for j in (i + 1)..N {
            let pj = ps[j];
            let dx = pi.x - pj.x;
            let dy = pi.y - pj.y;
            let dz = pi.z - pj.z;
            let softening2 = sim.softening * sim.softening;
            let r2 = dx * dx + dy * dy + dz * dz + softening2;
            V_newt -= G * pi.m * pj.m / r2.sqrt();
        }
    }

    // Transform to Jacobi coordinates
    let source = ps[0];
    let mu = G * source.m;
    let mut m_j: Vec<f64> = vec![0.; N];
    rebx_calculate_jacobi_masses(ps, &mut m_j, N);
    reb_transformations_inertial_to_jacobi_posvel(ps, &mut ps_j, ps, N, N);

    let mut T =
        0.5 * m_j[0] * (ps_j[0].vx * ps_j[0].vx + ps_j[0].vy * ps_j[0].vy + ps_j[0].vz * ps_j[0].vz);
    let mut V_PN = 0.;
    for i in 1..N {
        let p = ps_j[i];
        let rdoti2 = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
        let mut vtildei2 = rdoti2;
        let ri = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        let vscale2 = mu / ri; // characteristic v^2
        for _q in 0..10 {
            let old_vtildei2 = vtildei2;
            let A = (0.5 * vtildei2 + 3. * vscale2) / C2;
            vtildei2 = rdoti2 / ((1. - A) * (1. - A));
            if (vtildei2 - old_vtildei2) / vtildei2 < f64::EPSILON {
                break;
            }
        }

        V_PN += m_j[i]
            * (0.5 * mu * mu / (ri * ri) - 0.125 * vtildei2 * vtildei2 - 1.5 * mu * vtildei2 / ri);
        T += 0.5 * m_j[i] * vtildei2;
    }
    V_PN /= C2;

    T + V_newt + V_PN
}

/// gr.c `rebx_gr_hamiltonian`. Returns the value of the GR Hamiltonian
/// (kinetic + Newtonian potential + post-Newtonian potential).
///
/// `gr` is an index into `rebx_extras::allocated_forces` (C: the
/// `struct rebx_force*`); `sim` replaces the C's `rebx->sim` back-pointer.
pub fn rebx_gr_hamiltonian(rebx: &mut rebx_extras, sim: &reb_simulation, gr: usize) -> f64 {
    let c = rebx_get_param_double(rebx, rebx_ap::force(gr), "c");
    let c = match c {
        None => {
            rebx_error(
                rebx,
                "Need to set speed of light in gr effect.  See examples in documentation.\n",
            );
            return 0.;
        }
        Some(c) => c,
    };
    let C2 = c * c;
    rebx_calculate_gr_hamiltonian(rebx, sim, C2)
}
