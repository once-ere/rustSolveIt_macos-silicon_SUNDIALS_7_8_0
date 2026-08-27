//! gravity.rs — gravity modules (from gravity.c, serial/non-OpenMP
//! paths; the OpenMP branches change summation order and are excluded
//! by the Windows C build as well).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::boundary::{reb_boundary_check, reb_boundary_get_ghostbox};
use crate::tools::reb_simulation_error;
use crate::tree::*;
use crate::types::*;

/// gravity.c `reb_gravity_tree_calculate_acceleration`.
pub fn reb_gravity_tree_calculate_acceleration(r: &mut reb_simulation) {
    // Check if particles are in box
    reb_boundary_check(r);

    reb_tree_construct(r);
    if r.tree_root.is_empty() {
        reb_simulation_error(r, "Tree does not exist. Cannot calculate accelerations.");
        return;
    }

    // Update center of mass in tree in preparation of force calculation.
    reb_tree_calculate_gravity_data(r);

    let N = r.N;
    for i in 0..N {
        r.particles[i].ax = 0.;
        r.particles[i].ay = 0.;
        r.particles[i].az = 0.;
    }
    // Summing over all Ghost Boxes
    for gbx in -r.N_ghost_x..=r.N_ghost_x {
        for gby in -r.N_ghost_y..=r.N_ghost_y {
            for gbz in -r.N_ghost_z..=r.N_ghost_z {
                // Summing over all particle pairs
                for i in 0..N {
                    let mut gb = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                    // Precalculated shifted position
                    gb.x += r.particles[i].x;
                    gb.y += r.particles[i].y;
                    gb.z += r.particles[i].z;
                    reb_tree_calculate_acceleration_for_particle(r, i, &gb);
                }
            }
        }
    }
    // Delete tree (if it exists)
    reb_tree_delete(r);
}

/// gravity.c `reb_gravity_jacobi_calculate_acceleration`.
pub fn reb_gravity_jacobi_calculate_acceleration(r: &mut reb_simulation) {
    let N = r.N;
    let G = r.G;
    let mut Rjx = 0.;
    let mut Rjy = 0.;
    let mut Rjz = 0.;
    let mut Mj = 0.;
    for j in 0..N {
        r.particles[j].ax = 0.;
        r.particles[j].ay = 0.;
        r.particles[j].az = 0.;
        for i in 0..(j + 1) {
            if j > 1 {
                // Jacobi term (j==1 terms cancel and are skipped, as in C)
                let Qjx = r.particles[j].x - Rjx / Mj;
                let Qjy = r.particles[j].y - Rjy / Mj;
                let Qjz = r.particles[j].z - Rjz / Mj;
                let dr = (Qjx * Qjx + Qjy * Qjy + Qjz * Qjz).sqrt();
                let mut dQjdri = Mj;
                if i < j {
                    dQjdri = -r.particles[j].m;
                }
                let prefact = G * dQjdri / (dr * dr * dr);
                r.particles[i].ax += prefact * Qjx;
                r.particles[i].ay += prefact * Qjy;
                r.particles[i].az += prefact * Qjz;
            }
            if i != j && (i != 0 || j != 1) {
                // Direct term
                let dx = r.particles[i].x - r.particles[j].x;
                let dy = r.particles[i].y - r.particles[j].y;
                let dz = r.particles[i].z - r.particles[j].z;
                let dr = (dx * dx + dy * dy + dz * dz).sqrt();
                let prefact = G / (dr * dr * dr);
                let prefacti = prefact * r.particles[i].m;
                let prefactj = prefact * r.particles[j].m;

                r.particles[i].ax -= prefactj * dx;
                r.particles[i].ay -= prefactj * dy;
                r.particles[i].az -= prefactj * dz;
                r.particles[j].ax += prefacti * dx;
                r.particles[j].ay += prefacti * dy;
                r.particles[j].az += prefacti * dz;
            }
        }
        Rjx += r.particles[j].m * r.particles[j].x;
        Rjy += r.particles[j].m * r.particles[j].y;
        Rjz += r.particles[j].m * r.particles[j].z;
        Mj += r.particles[j].m;
    }
}

/// gravity.c `reb_gravity_basic_calculate_acceleration` (serial).
pub fn reb_gravity_basic_calculate_acceleration(r: &mut reb_simulation) {
    let N = r.N;
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let gravity_ignore_terms = r.gravity_ignore_terms;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let _testparticle_type = r.testparticle_type;
    let N_ghost_x = r.N_ghost_x;
    let N_ghost_y = r.N_ghost_y;
    let N_ghost_z = r.N_ghost_z;
    let starti: usize = if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_NONE { 1 } else { 2 };
    let startj: usize = if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 { 1 } else { 0 };
    for i in 0..N {
        r.particles[i].ax = 0.;
        r.particles[i].ay = 0.;
        r.particles[i].az = 0.;
    }
    // Summing over all Ghost Boxes
    for gbx in -N_ghost_x..=N_ghost_x {
        for gby in -N_ghost_y..=N_ghost_y {
            for gbz in -N_ghost_z..=N_ghost_z {
                let gb = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                // All active particle pairs, O(N^2/2)
                for i in starti..N_active {
                    for j in startj..i {
                        let dx = (gb.x + r.particles[i].x) - r.particles[j].x;
                        let dy = (gb.y + r.particles[i].y) - r.particles[j].y;
                        let dz = (gb.z + r.particles[i].z) - r.particles[j].z;
                        let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
                        let prefact = G / (_r * _r * _r);
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
                // Interactions of test particles with active particles
                let startitestp = N_active.max(starti);
                for i in startitestp..N {
                    for j in startj..N_active {
                        let dx = (gb.x + r.particles[i].x) - r.particles[j].x;
                        let dy = (gb.y + r.particles[i].y) - r.particles[j].y;
                        let dz = (gb.z + r.particles[i].z) - r.particles[j].z;
                        let _r = (dx * dx + dy * dy + dz * dz + softening2).sqrt();
                        let prefact = G / (_r * _r * _r);
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
            }
        }
    }
}

/// One Kahan-compensated accumulation into (ax,ay,az)/(cs) — the inner
/// block reused throughout gravity.c's compensated routine.
fn compensated_add(p: &mut reb_particle, cs: &mut reb_vec3d, ix: f64, iy: f64, iz: f64) {
    let yx = ix - cs.x;
    let tx = p.ax + yx;
    cs.x = (tx - p.ax) - yx;
    p.ax = tx;

    let yy = iy - cs.y;
    let ty = p.ay + yy;
    cs.y = (ty - p.ay) - yy;
    p.ay = ty;

    let yz = iz - cs.z;
    let tz = p.az + yz;
    cs.z = (tz - p.az) - yz;
    p.az = tz;
}

/// gravity.c `reb_gravity_compensated_calculate_acceleration` (serial).
pub fn reb_gravity_compensated_calculate_acceleration(r: &mut reb_simulation) {
    let N = r.N;
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let gravity_ignore_terms = r.gravity_ignore_terms;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let _testparticle_type = r.testparticle_type;
    if r.gravity_cs.len() < N {
        r.gravity_cs.resize(N, reb_vec3d::default());
    }
    for i in 0..N {
        r.particles[i].ax = 0.;
        r.particles[i].ay = 0.;
        r.particles[i].az = 0.;
        r.gravity_cs[i] = reb_vec3d::default();
    }
    // Summing over all massive particle pairs
    for i in 0..N_active {
        for j in (i + 1)..N_active {
            if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1
                && ((j == 1 && i == 0) || (i == 1 && j == 0))
            {
                continue;
            }
            if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && (j == 0 || i == 0) {
                continue;
            }
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let r2 = dx * dx + dy * dy + dz * dz + softening2;
            let _r = r2.sqrt();
            let prefact = G / (r2 * _r);
            let prefacti = prefact * r.particles[i].m;
            let prefactj = -prefact * r.particles[j].m;

            {
                let (particles, cs) = (&mut r.particles, &mut r.gravity_cs);
                compensated_add(&mut particles[i], &mut cs[i], prefactj * dx, prefactj * dy, prefactj * dz);
                compensated_add(&mut particles[j], &mut cs[j], prefacti * dx, prefacti * dy, prefacti * dz);
            }
        }
    }
    // Testparticles
    for i in N_active..N {
        for j in 0..N_active {
            if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1
                && ((j == 1 && i == 0) || (i == 1 && j == 0))
            {
                continue;
            }
            if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && (j == 0 || i == 0) {
                continue;
            }
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;
            let r2 = dx * dx + dy * dy + dz * dz + softening2;
            let _r = r2.sqrt();
            let prefact = G / (r2 * _r);
            let prefactj = -prefact * r.particles[j].m;

            {
                let (particles, cs) = (&mut r.particles, &mut r.gravity_cs);
                compensated_add(&mut particles[i], &mut cs[i], prefactj * dx, prefactj * dy, prefactj * dz);
            }
            if _testparticle_type != 0 {
                let prefacti = prefact * r.particles[i].m;
                let (particles, cs) = (&mut r.particles, &mut r.gravity_cs);
                compensated_add(&mut particles[j], &mut cs[j], prefacti * dx, prefacti * dy, prefacti * dz);
            }
        }
    }
}

/// gravity.c `reb_gravity_basic_calculate_acceleration_var`
/// (first- and second-order variational accelerations; serial paths).
pub fn reb_gravity_basic_calculate_acceleration_var(r: &mut reb_simulation) {
    let G = r.G;
    let gravity_ignore_terms = r.gravity_ignore_terms;
    let _testparticle_type = r.testparticle_type;
    let N = r.N;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let starti: usize = if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_NONE { 1 } else { 2 };
    let startj: usize = if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 { 1 } else { 0 };
    for v in 0..r.var_config.len() {
        let vc = r.var_config[v];
        if vc.order == 1 {
            let base = vc.index;
            if vc.testparticle < 0 {
                for i in 0..N {
                    r.particles_var[base + i].ax = 0.;
                    r.particles_var[base + i].ay = 0.;
                    r.particles_var[base + i].az = 0.;
                }
                for i in starti..N_active {
                    for j in startj..i {
                        var1_pair(r, G, base, i, j, true);
                    }
                }
                for i in N_active..N {
                    for j in startj..N_active {
                        var1_pair(r, G, base, i, j, _testparticle_type != 0);
                    }
                }
            } else {
                // testparticle
                let i = vc.testparticle as usize;
                r.particles_var[base].ax = 0.;
                r.particles_var[base].ay = 0.;
                r.particles_var[base].az = 0.;
                for j in 0..N {
                    if i == j {
                        continue;
                    }
                    if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1
                        && ((j == 1 && i == 0) || (i == 1 && j == 0))
                    {
                        continue;
                    }
                    if gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0
                        && (j == 0 || i == 0)
                    {
                        continue;
                    }
                    let dx = r.particles[i].x - r.particles[j].x;
                    let dy = r.particles[i].y - r.particles[j].y;
                    let dz = r.particles[i].z - r.particles[j].z;
                    let r2 = dx * dx + dy * dy + dz * dz;
                    let _r = r2.sqrt();
                    let r3inv = 1. / (r2 * _r);
                    let r5inv = 3. * r3inv / r2;
                    let ddx = r.particles_var[base].x;
                    let ddy = r.particles_var[base].y;
                    let ddz = r.particles_var[base].z;
                    let Gmj = G * r.particles[j].m;

                    let dxdx = dx * dx * r5inv - r3inv;
                    let dydy = dy * dy * r5inv - r3inv;
                    let dzdz = dz * dz * r5inv - r3inv;
                    let dxdy = dx * dy * r5inv;
                    let dxdz = dx * dz * r5inv;
                    let dydz = dy * dz * r5inv;
                    let dax = ddx * dxdx + ddy * dxdy + ddz * dxdz;
                    let day = ddx * dxdy + ddy * dydy + ddz * dydz;
                    let daz = ddx * dxdz + ddy * dydz + ddz * dzdz;

                    r.particles_var[base].ax += Gmj * dax;
                    r.particles_var[base].ay += Gmj * day;
                    r.particles_var[base].az += Gmj * daz;
                }
            }
        } else if vc.order == 2 {
            if _testparticle_type != 0 {
                reb_simulation_error(
                    r,
                    "testparticletype=1 not implemented for second order variational equations.",
                );
            }
            if gravity_ignore_terms != REB_GRAVITY_IGNORE_TERMS_NONE {
                reb_simulation_error(
                    r,
                    "Second order variational equations do not support gravity_ignore_terms",
                );
            }
            let base2 = vc.index;
            let base1a = vc.index_1st_order_a;
            let base1b = vc.index_1st_order_b;
            if vc.testparticle < 0 {
                for i in 0..N {
                    r.particles_var[base2 + i].ax = 0.;
                    r.particles_var[base2 + i].ay = 0.;
                    r.particles_var[base2 + i].az = 0.;
                }
                for i in 0..N {
                    for j in (i + 1)..N {
                        var2_pair(r, G, base2, base1a, base1b, i, j);
                    }
                }
            } else {
                // testparticle
                let i = vc.testparticle as usize;
                r.particles_var[base2].ax = 0.;
                r.particles_var[base2].ay = 0.;
                r.particles_var[base2].az = 0.;
                for j in 0..N {
                    if i == j {
                        continue;
                    }
                    let dx = r.particles[i].x - r.particles[j].x;
                    let dy = r.particles[i].y - r.particles[j].y;
                    let dz = r.particles[i].z - r.particles[j].z;
                    let r2 = dx * dx + dy * dy + dz * dz;
                    let _r = r2.sqrt();
                    let r3inv = 1. / (r2 * _r);
                    let r5inv = r3inv / r2;
                    let r7inv = r5inv / r2;
                    let ddx = r.particles_var[base2].x;
                    let ddy = r.particles_var[base2].y;
                    let ddz = r.particles_var[base2].z;
                    let Gmj = G * r.particles[j].m;

                    let mut dax = ddx * (3. * dx * dx * r5inv - r3inv)
                        + ddy * (3. * dx * dy * r5inv)
                        + ddz * (3. * dx * dz * r5inv);
                    let mut day = ddx * (3. * dy * dx * r5inv)
                        + ddy * (3. * dy * dy * r5inv - r3inv)
                        + ddz * (3. * dy * dz * r5inv);
                    let mut daz = ddx * (3. * dz * dx * r5inv)
                        + ddy * (3. * dz * dy * r5inv)
                        + ddz * (3. * dz * dz * r5inv - r3inv);

                    let dk1dx = r.particles_var[base1a].x;
                    let dk1dy = r.particles_var[base1a].y;
                    let dk1dz = r.particles_var[base1a].z;
                    let dk2dx = r.particles_var[base1b].x;
                    let dk2dy = r.particles_var[base1b].y;
                    let dk2dz = r.particles_var[base1b].z;

                    let rdk1 = dx * dk1dx + dy * dk1dy + dz * dk1dz;
                    let rdk2 = dx * dk2dx + dy * dk2dy + dz * dk2dz;
                    let dk1dk2 = dk1dx * dk2dx + dk1dy * dk2dy + dk1dz * dk2dz;
                    dax += 3. * r5inv * dk2dx * rdk1
                        + 3. * r5inv * dk1dx * rdk2
                        + 3. * r5inv * dx * dk1dk2
                        - 15. * dx * r7inv * rdk1 * rdk2;
                    day += 3. * r5inv * dk2dy * rdk1
                        + 3. * r5inv * dk1dy * rdk2
                        + 3. * r5inv * dy * dk1dk2
                        - 15. * dy * r7inv * rdk1 * rdk2;
                    daz += 3. * r5inv * dk2dz * rdk1
                        + 3. * r5inv * dk1dz * rdk2
                        + 3. * r5inv * dz * dk1dk2
                        - 15. * dz * r7inv * rdk1 * rdk2;

                    r.particles_var[base2].ax += Gmj * dax;
                    r.particles_var[base2].ay += Gmj * day;
                    r.particles_var[base2].az += Gmj * daz;
                }
            }
        }
    }
}

/// The 1st-order pair interaction of
/// `reb_gravity_basic_calculate_acceleration_var` (the two C loops
/// share this body; `apply_j` is the `_testparticle_type` gate of the
/// second loop and constant `true` in the first).
fn var1_pair(r: &mut reb_simulation, G: f64, base: usize, i: usize, j: usize, apply_j: bool) {
    let dx = r.particles[i].x - r.particles[j].x;
    let dy = r.particles[i].y - r.particles[j].y;
    let dz = r.particles[i].z - r.particles[j].z;
    let r2 = dx * dx + dy * dy + dz * dz;
    let _r = r2.sqrt();
    let r3inv = 1. / (r2 * _r);
    let r5inv = 3. * r3inv / r2;
    let ddx = r.particles_var[base + i].x - r.particles_var[base + j].x;
    let ddy = r.particles_var[base + i].y - r.particles_var[base + j].y;
    let ddz = r.particles_var[base + i].z - r.particles_var[base + j].z;
    let Gmi = G * r.particles[i].m;
    let Gmj = G * r.particles[j].m;

    let dxdx = dx * dx * r5inv - r3inv;
    let dydy = dy * dy * r5inv - r3inv;
    let dzdz = dz * dz * r5inv - r3inv;
    let dxdy = dx * dy * r5inv;
    let dxdz = dx * dz * r5inv;
    let dydz = dy * dz * r5inv;
    let dax = ddx * dxdx + ddy * dxdy + ddz * dxdz;
    let day = ddx * dxdy + ddy * dydy + ddz * dydz;
    let daz = ddx * dxdz + ddy * dydz + ddz * dzdz;

    let dGmi = G * r.particles_var[base + i].m;
    let dGmj = G * r.particles_var[base + j].m;

    r.particles_var[base + i].ax += Gmj * dax - dGmj * r3inv * dx;
    r.particles_var[base + i].ay += Gmj * day - dGmj * r3inv * dy;
    r.particles_var[base + i].az += Gmj * daz - dGmj * r3inv * dz;

    if apply_j {
        r.particles_var[base + j].ax -= Gmi * dax - dGmi * r3inv * dx;
        r.particles_var[base + j].ay -= Gmi * day - dGmi * r3inv * dy;
        r.particles_var[base + j].az -= Gmi * daz - dGmi * r3inv * dz;
    }
}

/// The 2nd-order pair interaction of
/// `reb_gravity_basic_calculate_acceleration_var`.
fn var2_pair(r: &mut reb_simulation, G: f64, base2: usize, base1a: usize, base1b: usize, i: usize, j: usize) {
    let dx = r.particles[i].x - r.particles[j].x;
    let dy = r.particles[i].y - r.particles[j].y;
    let dz = r.particles[i].z - r.particles[j].z;
    let r2 = dx * dx + dy * dy + dz * dz;
    let _r = r2.sqrt();
    let r3inv = 1. / (r2 * _r);
    let r5inv = r3inv / r2;
    let r7inv = r5inv / r2;
    let ddx = r.particles_var[base2 + i].x - r.particles_var[base2 + j].x;
    let ddy = r.particles_var[base2 + i].y - r.particles_var[base2 + j].y;
    let ddz = r.particles_var[base2 + i].z - r.particles_var[base2 + j].z;
    let Gmi = G * r.particles[i].m;
    let Gmj = G * r.particles[j].m;
    let ddGmi = G * r.particles_var[base2 + i].m;
    let ddGmj = G * r.particles_var[base2 + j].m;

    let mut dax = ddx * (3. * dx * dx * r5inv - r3inv)
        + ddy * (3. * dx * dy * r5inv)
        + ddz * (3. * dx * dz * r5inv);
    let mut day = ddx * (3. * dy * dx * r5inv)
        + ddy * (3. * dy * dy * r5inv - r3inv)
        + ddz * (3. * dy * dz * r5inv);
    let mut daz = ddx * (3. * dz * dx * r5inv)
        + ddy * (3. * dz * dy * r5inv)
        + ddz * (3. * dz * dz * r5inv - r3inv);

    let dk1dx = r.particles_var[base1a + i].x - r.particles_var[base1a + j].x;
    let dk1dy = r.particles_var[base1a + i].y - r.particles_var[base1a + j].y;
    let dk1dz = r.particles_var[base1a + i].z - r.particles_var[base1a + j].z;
    let dk2dx = r.particles_var[base1b + i].x - r.particles_var[base1b + j].x;
    let dk2dy = r.particles_var[base1b + i].y - r.particles_var[base1b + j].y;
    let dk2dz = r.particles_var[base1b + i].z - r.particles_var[base1b + j].z;

    let rdk1 = dx * dk1dx + dy * dk1dy + dz * dk1dz;
    let rdk2 = dx * dk2dx + dy * dk2dy + dz * dk2dz;
    let dk1dk2 = dk1dx * dk2dx + dk1dy * dk2dy + dk1dz * dk2dz;
    dax += 3. * r5inv * dk2dx * rdk1
        + 3. * r5inv * dk1dx * rdk2
        + 3. * r5inv * dx * dk1dk2
        - 15. * dx * r7inv * rdk1 * rdk2;
    day += 3. * r5inv * dk2dy * rdk1
        + 3. * r5inv * dk1dy * rdk2
        + 3. * r5inv * dy * dk1dk2
        - 15. * dy * r7inv * rdk1 * rdk2;
    daz += 3. * r5inv * dk2dz * rdk1
        + 3. * r5inv * dk1dz * rdk2
        + 3. * r5inv * dz * dk1dk2
        - 15. * dz * r7inv * rdk1 * rdk2;

    let dk1Gmi = G * r.particles_var[base1a + i].m;
    let dk1Gmj = G * r.particles_var[base1a + j].m;
    let dk2Gmi = G * r.particles_var[base1b + i].m;
    let dk2Gmj = G * r.particles_var[base1b + j].m;

    r.particles_var[base2 + i].ax += Gmj * dax
        - ddGmj * r3inv * dx
        - dk2Gmj * r3inv * dk1dx + 3. * dk2Gmj * r5inv * dx * rdk1
        - dk1Gmj * r3inv * dk2dx + 3. * dk1Gmj * r5inv * dx * rdk2;
    r.particles_var[base2 + i].ay += Gmj * day
        - ddGmj * r3inv * dy
        - dk2Gmj * r3inv * dk1dy + 3. * dk2Gmj * r5inv * dy * rdk1
        - dk1Gmj * r3inv * dk2dy + 3. * dk1Gmj * r5inv * dy * rdk2;
    r.particles_var[base2 + i].az += Gmj * daz
        - ddGmj * r3inv * dz
        - dk2Gmj * r3inv * dk1dz + 3. * dk2Gmj * r5inv * dz * rdk1
        - dk1Gmj * r3inv * dk2dz + 3. * dk1Gmj * r5inv * dz * rdk2;

    r.particles_var[base2 + j].ax -= Gmi * dax
        - ddGmi * r3inv * dx
        - dk2Gmi * r3inv * dk1dx + 3. * dk2Gmi * r5inv * dx * rdk1
        - dk1Gmi * r3inv * dk2dx + 3. * dk1Gmi * r5inv * dx * rdk2;
    r.particles_var[base2 + j].ay -= Gmi * day
        - ddGmi * r3inv * dy
        - dk2Gmi * r3inv * dk1dy + 3. * dk2Gmi * r5inv * dy * rdk1
        - dk1Gmi * r3inv * dk2dy + 3. * dk1Gmi * r5inv * dy * rdk2;
    r.particles_var[base2 + j].az -= Gmi * daz
        - ddGmi * r3inv * dz
        - dk2Gmi * r3inv * dk1dz + 3. * dk2Gmi * r5inv * dz * rdk1
        - dk1Gmi * r3inv * dk2dz + 3. * dk1Gmi * r5inv * dz * rdk2;
}

/// gravity.c `reb_gravity_basic_calculate_and_apply_jerk` (serial).
pub fn reb_gravity_basic_calculate_and_apply_jerk(r: &mut reb_simulation, v: f64) {
    let N = r.N;
    let G = r.G;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let _testparticle_type = r.testparticle_type;
    let starti: usize = if r.gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_NONE { 1 } else { 2 };
    let startj: usize = if r.gravity_ignore_terms == REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 { 1 } else { 0 };
    // All interactions between active particles
    for i in starti..N_active {
        for j in startj..i {
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;

            let dax = r.particles[i].ax - r.particles[j].ax;
            let day = r.particles[i].ay - r.particles[j].ay;
            let daz = r.particles[i].az - r.particles[j].az;

            let dr = (dx * dx + dy * dy + dz * dz).sqrt();
            let alphasum = dax * dx + day * dy + daz * dz;
            let prefact2 = 2. * v * G / (dr * dr * dr);
            let prefact2i = prefact2 * r.particles[j].m;
            let prefact2j = prefact2 * r.particles[i].m;
            let prefact1 = alphasum * prefact2 / dr * 3. / dr;
            let prefact1i = prefact1 * r.particles[j].m;
            let prefact1j = prefact1 * r.particles[i].m;
            r.particles[i].vx += dx * prefact1i - dax * prefact2i;
            r.particles[i].vy += dy * prefact1i - day * prefact2i;
            r.particles[i].vz += dz * prefact1i - daz * prefact2i;
            r.particles[j].vx += dax * prefact2j - dx * prefact1j;
            r.particles[j].vy += day * prefact2j - dy * prefact1j;
            r.particles[j].vz += daz * prefact2j - dz * prefact1j;
        }
    }
    // Interactions between active particles and test particles
    for i in N_active..N {
        for j in startj..i {
            let dx = r.particles[i].x - r.particles[j].x;
            let dy = r.particles[i].y - r.particles[j].y;
            let dz = r.particles[i].z - r.particles[j].z;

            let dax = r.particles[i].ax - r.particles[j].ax;
            let day = r.particles[i].ay - r.particles[j].ay;
            let daz = r.particles[i].az - r.particles[j].az;

            let dr = (dx * dx + dy * dy + dz * dz).sqrt();
            let alphasum = dax * dx + day * dy + daz * dz;
            let prefact2 = 2. * v * G / (dr * dr * dr);
            let prefact1 = alphasum * prefact2 / dr * 3. / dr;
            let prefact1i = prefact1 * r.particles[j].m;
            let prefact2i = prefact2 * r.particles[j].m;
            r.particles[i].vx += dx * prefact1i - dax * prefact2i;
            r.particles[i].vy += dy * prefact1i - day * prefact2i;
            r.particles[i].vz += dz * prefact1i - daz * prefact2i;
            if _testparticle_type != 0 {
                let prefact1j = prefact1 * r.particles[i].m;
                let prefact2j = prefact2 * r.particles[i].m;
                r.particles[j].vx += dax * prefact2j - dx * prefact1j;
                r.particles[j].vy += day * prefact2j - dy * prefact1j;
                r.particles[j].vz += daz * prefact2j - dz * prefact1j;
            }
        }
    }
}
