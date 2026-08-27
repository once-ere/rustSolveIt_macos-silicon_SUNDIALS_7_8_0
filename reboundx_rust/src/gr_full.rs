//! gr_full.rs — translation of REBOUNDx gr_full.c
//! Post-newtonian general relativity corrections for all bodies in the simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # General Relativity
//!
//! ======================= ===============================================
//! Authors                 P. Shi, H. Rein, D. Tamayo
//! Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
//! Based on                `Newhall et al. 1983 <http://labs.adsabs.harvard.edu/adsabs/abs/1983A%26A...125..150N/>`_.
//! C Example               :ref:`c_example_gr`
//! Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
//! ======================= ===============================================
//!
//! This algorithm incorporates the first-order post-newtonian effects from all
//! bodies in the system, and is necessary for multiple massive bodies like
//! stellar binaries.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! c (double)                   Yes         Speed of light, needs to be specified in the units used for the simulation.
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! *None*
//!
//! (The C also reads an undocumented effect parameter `max_iterations` (int)
//! and passes it to the worker; see [`rebx_gr_full`]. As in the C, the worker
//! ignores it and always performs at most 10 substitution passes.)

use crate::core::{rebx_error, rebx_get_param_double, rebx_get_param_int};
use crate::types::{rebx_ap, rebx_extras};

use rebound_rs::{
    reb_particle, reb_particle_isub, reb_simulation, reb_simulation_com, reb_simulation_error,
    reb_simulation_warning, reb_vec3d, REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1,
};

/// gr_full.c `rebx_calculate_gr_full` (static).
///
/// The C takes `particles` as a separate pointer; it always aliases
/// `sim->particles`, so here it is reached through `sim` (`sim.particles`).
/// `max_iterations` and `gravity_ignore_10` are accepted and ignored exactly
/// as in the C body.
fn rebx_calculate_gr_full(
    sim: &mut reb_simulation,
    N: usize,
    C2: f64,
    G: f64,
    max_iterations: i32,
    gravity_ignore_10: u32,
) {
    // Both are unused in the C body as well; bound here so the signature can
    // stay identical to the C's without tripping the unused-variable lint.
    let _ = max_iterations;
    let _ = gravity_ignore_10;

    // array that stores the value of the constant term
    let mut a_const: Vec<[f64; 3]> = vec![[0.; 3]; N];
    let mut ps_b: Vec<reb_particle> = sim.particles[0..N].to_vec();

    // Calculate Newtonian accelerations
    for i in 0..N {
        ps_b[i].ax = 0.;
        ps_b[i].ay = 0.;
        ps_b[i].az = 0.;
    }

    for i in 0..N {
        let pi = ps_b[i];
        for j in (i + 1)..N {
            let pj = ps_b[j];
            let dx = pi.x - pj.x;
            let dy = pi.y - pj.y;
            let dz = pi.z - pj.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            let prefac = G / (r2 * r);
            ps_b[i].ax -= prefac * pj.m * dx;
            ps_b[i].ay -= prefac * pj.m * dy;
            ps_b[i].az -= prefac * pj.m * dz;
            ps_b[j].ax += prefac * pi.m * dx;
            ps_b[j].ay += prefac * pi.m * dy;
            ps_b[j].az += prefac * pi.m * dz;
        }
    }

    // Transform to barycentric coordinates
    let com = reb_simulation_com(sim);
    for i in 0..N {
        reb_particle_isub(&mut ps_b[i], &com);
    }
    for i in 0..N {
        // then compute the constant terms:
        let mut a_constx = 0.;
        let mut a_consty = 0.;
        let mut a_constz = 0.;
        // 1st constant part
        for j in 0..N {
            if j != i {
                let dxij = ps_b[i].x - ps_b[j].x;
                let dyij = ps_b[i].y - ps_b[j].y;
                let dzij = ps_b[i].z - ps_b[j].z;
                let rij2 = dxij * dxij + dyij * dyij + dzij * dzij;
                let rij = rij2.sqrt();
                let rij3 = rij2 * rij;

                let mut a1 = 0.;
                for k in 0..N {
                    if k != i {
                        let dxik = ps_b[i].x - ps_b[k].x;
                        let dyik = ps_b[i].y - ps_b[k].y;
                        let dzik = ps_b[i].z - ps_b[k].z;
                        let rik = (dxik * dxik + dyik * dyik + dzik * dzik).sqrt();
                        a1 += (4. / (C2)) * G * sim.particles[k].m / rik;
                    }
                }

                let mut a2 = 0.;
                for l in 0..N {
                    if l != j {
                        let dxlj = ps_b[l].x - ps_b[j].x;
                        let dylj = ps_b[l].y - ps_b[j].y;
                        let dzlj = ps_b[l].z - ps_b[j].z;
                        let rlj = (dxlj * dxlj + dylj * dylj + dzlj * dzlj).sqrt();
                        a2 += (1. / (C2)) * G * sim.particles[l].m / rlj;
                    }
                }

                let vi2 =
                    ps_b[i].vx * ps_b[i].vx + ps_b[i].vy * ps_b[i].vy + ps_b[i].vz * ps_b[i].vz;
                let a3 = -vi2 / (C2);

                let vj2 =
                    ps_b[j].vx * ps_b[j].vx + ps_b[j].vy * ps_b[j].vy + ps_b[j].vz * ps_b[j].vz;
                let a4 = -2. * vj2 / (C2);

                let a5 = (4. / (C2))
                    * (ps_b[i].vx * ps_b[j].vx
                        + ps_b[i].vy * ps_b[j].vy
                        + ps_b[i].vz * ps_b[j].vz);

                let a6_0 = dxij * ps_b[j].vx + dyij * ps_b[j].vy + dzij * ps_b[j].vz;
                let a6 = (3. / (2. * C2)) * a6_0 * a6_0 / rij2;

                // Newtonian piece of first ddot(r) piece
                let a7 = (dxij * ps_b[j].ax + dyij * ps_b[j].ay + dzij * ps_b[j].az) / (2. * C2);

                let factor1 = a1 + a2 + a3 + a4 + a5 + a6 + a7;

                a_constx += G * sim.particles[j].m * dxij * factor1 / rij3;
                a_consty += G * sim.particles[j].m * dyij * factor1 / rij3;
                a_constz += G * sim.particles[j].m * dzij * factor1 / rij3;

                // 2nd constant part

                let dvxij = ps_b[i].vx - ps_b[j].vx;
                let dvyij = ps_b[i].vy - ps_b[j].vy;
                let dvzij = ps_b[i].vz - ps_b[j].vz;

                let factor2 = dxij * (4. * ps_b[i].vx - 3. * ps_b[j].vx)
                    + dyij * (4. * ps_b[i].vy - 3. * ps_b[j].vy)
                    + dzij * (4. * ps_b[i].vz - 3. * ps_b[j].vz);

                a_constx += G * sim.particles[j].m / C2
                    * (factor2 * dvxij / rij3 + 7. / 2. * ps_b[j].ax / rij);
                a_consty += G * sim.particles[j].m / C2
                    * (factor2 * dvyij / rij3 + 7. / 2. * ps_b[j].ay / rij);
                a_constz += G * sim.particles[j].m / C2
                    * (factor2 * dvzij / rij3 + 7. / 2. * ps_b[j].az / rij);
            }
        }

        a_const[i][0] = a_constx;
        a_const[i][1] = a_consty;
        a_const[i][2] = a_constz;
    }
    for i in 0..N {
        ps_b[i].ax = a_const[i][0];
        ps_b[i].ay = a_const[i][1];
        ps_b[i].az = a_const[i][2];
    }

    // Now running the substitution again and again through the loop below
    for k in 0..10 {
        // you can set k as how many substitution you want to make
        // initialize an arry that stores the information of previousu calculated accleration
        let mut a_old: Vec<[f64; 3]> = vec![[0.; 3]; N];
        for i in 0..N {
            a_old[i][0] = ps_b[i].ax;
            a_old[i][1] = ps_b[i].ay;
            a_old[i][2] = ps_b[i].az;
        }
        // now add on the non-constant term
        for i in 0..N {
            // a_j is used to update a_i and vice versa
            let mut non_constx = 0.;
            let mut non_consty = 0.;
            let mut non_constz = 0.;
            for j in 0..N {
                if j != i {
                    let dxij = ps_b[i].x - ps_b[j].x;
                    let dyij = ps_b[i].y - ps_b[j].y;
                    let dzij = ps_b[i].z - ps_b[j].z;
                    let rij = (dxij * dxij + dyij * dyij + dzij * dzij).sqrt();
                    let rij3 = rij * rij * rij;
                    let dotproduct = dxij * ps_b[j].ax + dyij * ps_b[j].ay + dzij * ps_b[j].az;

                    non_constx += (G * sim.particles[j].m * dxij / rij3) * dotproduct / (2. * C2)
                        + (7. / (2. * C2)) * G * sim.particles[j].m * ps_b[j].ax / rij;
                    non_consty += (G * sim.particles[j].m * dyij / rij3) * dotproduct / (2. * C2)
                        + (7. / (2. * C2)) * G * sim.particles[j].m * ps_b[j].ay / rij;
                    non_constz += (G * sim.particles[j].m * dzij / rij3) * dotproduct / (2. * C2)
                        + (7. / (2. * C2)) * G * sim.particles[j].m * ps_b[j].az / rij;
                }
            }
            ps_b[i].ax = a_const[i][0] + non_constx;
            ps_b[i].ay = a_const[i][1] + non_consty;
            ps_b[i].az = a_const[i][2] + non_constz;
        }

        // break out loop if ps_b is converging
        let mut maxdev = 0.;
        // The C declares dx/dy/dz once outside this loop; each is fully
        // assigned before it is read on every iteration, so scoping them to
        // the loop body is equivalent.
        for i in 0..N {
            let dx = if ps_b[i].ax.abs() < f64::EPSILON {
                0.
            } else {
                ((ps_b[i].ax - a_old[i][0]) / ps_b[i].ax).abs()
            };
            let dy = if ps_b[i].ay.abs() < f64::EPSILON {
                0.
            } else {
                ((ps_b[i].ay - a_old[i][1]) / ps_b[i].ay).abs()
            };
            let dz = if ps_b[i].az.abs() < f64::EPSILON {
                0.
            } else {
                ((ps_b[i].az - a_old[i][2]) / ps_b[i].az).abs()
            };

            if dx > maxdev {
                maxdev = dx;
            }
            if dy > maxdev {
                maxdev = dy;
            }
            if dz > maxdev {
                maxdev = dz;
            }
        }

        if maxdev < f64::EPSILON {
            break;
        }
        if k == 9 {
            reb_simulation_warning(sim, "10 loops in rebx_gr_full did not converge.\n");
            eprintln!("Fractional Error: {:e}", maxdev);
        }
    }

    for i in 0..N {
        sim.particles[i].ax += ps_b[i].ax;
        sim.particles[i].ay += ps_b[i].ay;
        sim.particles[i].az += ps_b[i].az;
    }

    // (C frees ps_b here; the Vec is dropped at end of scope.)
}

/// gr_full.c `rebx_gr_full`.
pub fn rebx_gr_full(sim: &mut reb_simulation, rebx: &mut rebx_extras, force_idx: usize, N: usize) {
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
    let gravity_ignore_10: u32 =
        (sim.gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1) as u32;
    let max_iterations = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "max_iterations");
    match max_iterations {
        Some(max_iterations) => {
            let G = sim.G;
            rebx_calculate_gr_full(sim, N, C2, G, max_iterations, gravity_ignore_10);
        }
        None => {
            let default_max_iterations: i32 = 10;
            let G = sim.G;
            rebx_calculate_gr_full(sim, N, C2, G, default_max_iterations, gravity_ignore_10);
        }
    }
}

/// gr_full.c `rebx_gr_full_hamiltonian`.
///
/// The C reaches the simulation through `rebx->sim` and returns 0 (after
/// `rebx_error`) when that back-pointer is NULL. `rebx_extras` carries no
/// back-pointer here, so the simulation is passed explicitly and that
/// unreachable branch is dropped.
pub fn rebx_gr_full_hamiltonian(
    sim: &reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
) -> f64 {
    let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "c");
    let c = match c {
        None => {
            rebx_error(
                rebx,
                "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n",
            );
            return 0.;
        }
        Some(c) => c,
    };
    let C2 = c * c;
    let N = sim.N;
    let G = sim.G;
    let particles = &sim.particles;

    let mut e_kin = 0.;
    let mut e_pot = 0.;
    let mut e_pn = 0.;

    let mut vtilde: Vec<reb_vec3d> = vec![reb_vec3d::default(); N];
    // The C also allocates `vtilde_old` here; it is never read or written.
    for j in 0..N {
        vtilde[j].x = particles[j].vx;
        vtilde[j].y = particles[j].vy;
        vtilde[j].z = particles[j].vz;
    }

    for _q in 0..10 {
        for i in 0..N {
            let pi = particles[i];

            let vtildei2 =
                vtilde[i].x * vtilde[i].x + vtilde[i].y * vtilde[i].y + vtilde[i].z * vtilde[i].z;
            let A = 1. - 0.5 * vtildei2 / C2;

            // The result of the following calculation is never used.
            //double sumk = 0.;
            //for (int k=0;k<N;k++){
            //    if (k!=i){
            //        struct reb_particle pk = particles[k];
            //        double xik = pk.x - pi.x;
            //        double yik = pk.y - pi.y;
            //        double zik = pk.z - pi.z;
            //        double rik = sqrt(xik*xik + yik*yik + zik*zik);
            //        sumk -= 2.*G*pk.m/rik;
            //    }
            //}

            let mut dv_pn = reb_vec3d {
                x: 0.,
                y: 0.,
                z: 0.,
            };
            for j in 0..N {
                if j != i {
                    let pj = particles[j];
                    let xij = pj.x - pi.x;
                    let yij = pj.y - pi.y;
                    let zij = pj.z - pi.z;
                    let rij2 = xij * xij + yij * yij + zij * zij;
                    let rij = rij2.sqrt();
                    let rijdotvj = vtilde[j].x * xij + vtilde[j].y * yij + vtilde[j].z * zij;
                    let pfac = pj.m / rij;

                    dv_pn.x += pfac * (6. * vtilde[i].x - 7. * vtilde[j].x - rijdotvj * xij / rij2);
                    dv_pn.y += pfac * (6. * vtilde[i].y - 7. * vtilde[j].y - rijdotvj * yij / rij2);
                    dv_pn.z += pfac * (6. * vtilde[i].z - 7. * vtilde[j].z - rijdotvj * zij / rij2);
                }
            }

            dv_pn.x *= G / (2. * C2);
            dv_pn.y *= G / (2. * C2);
            dv_pn.z *= G / (2. * C2);

            vtilde[i].x = (pi.vx + dv_pn.x) / A;
            vtilde[i].y = (pi.vy + dv_pn.y) / A;
            vtilde[i].z = (pi.vz + dv_pn.z) / A;
        }
    }

    for i in 0..N {
        let p = particles[i];
        let vtildei2 =
            vtilde[i].x * vtilde[i].x + vtilde[i].y * vtilde[i].y + vtilde[i].z * vtilde[i].z;
        e_kin += 0.5 * p.m * vtildei2;
    }

    for i in 0..N {
        let pi = particles[i];
        let mut sumk = 0.;
        for k in 0..N {
            if k != i {
                let pk = particles[k];
                let xik = pk.x - pi.x;
                let yik = pk.y - pi.y;
                let zik = pk.z - pi.z;
                let rik = (xik * xik + yik * yik + zik * zik).sqrt();
                sumk -= 2. * G * pk.m / rik;
            }
        }

        let vtildei2 =
            vtilde[i].x * vtilde[i].x + vtilde[i].y * vtilde[i].y + vtilde[i].z * vtilde[i].z;

        for j in 0..N {
            if j != i {
                let pj = particles[j];
                let xij = pj.x - pi.x;
                let yij = pj.y - pi.y;
                let zij = pj.z - pi.z;
                let rij2 = xij * xij + yij * yij + zij * zij;
                let rij = rij2.sqrt();
                let rijdotvj = vtilde[j].x * xij + vtilde[j].y * yij + vtilde[j].z * zij;
                let rijdotvi = vtilde[i].x * xij + vtilde[i].y * yij + vtilde[i].z * zij;
                let vidotvj = vtilde[i].x * vtilde[j].x
                    + vtilde[i].y * vtilde[j].y
                    + vtilde[i].z * vtilde[j].z;

                e_pn -= G / (4. * C2) * pi.m * pj.m / rij
                    * (6. * vtildei2 - 7. * vidotvj - rijdotvi * rijdotvj / rij2 + sumk);
            }
        }

        e_pn -= pi.m / (8. * C2) * vtildei2 * vtildei2;

        for j in (i + 1)..N {
            // classic full
            let pj = particles[j];
            let dx = pi.x - pj.x;
            let dy = pi.y - pj.y;
            let dz = pi.z - pj.z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            e_pot -= G * pi.m * pj.m / r;
        }
    }

    e_kin + e_pot + e_pn
}
