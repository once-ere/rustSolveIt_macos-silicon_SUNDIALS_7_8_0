//! collision.rs — collision search and resolution (from collision.c;
//! serial paths, MPI branches excluded like the Windows C build).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::boundary::reb_boundary_get_ghostbox;
use crate::particle::{reb_simulation_remove_particle, reb_simulation_two_largest_particles};
use crate::tools::rand_r;
use crate::tree::{reb_tree_construct, reb_tree_delete};
use crate::types::*;

/// collision.c `reb_collision_search` — finds collisions with the
/// selected module, shuffles them with `rand_r` (consuming the
/// simulation's RNG stream exactly like the C), then resolves each.
pub fn reb_collision_search(r: &mut reb_simulation) {
    r.N_collisions = 0;
    r.collisions.clear();
    let N_root = r.N_root_x * r.N_root_y * r.N_root_z;

    let N_projectiles = if r.map.is_some() { r.N_map } else { r.N };
    let N_targets = if r.N_targets != usize::MAX { r.N_targets } else { N_projectiles };

    match r.collision {
        REB_COLLISION::NONE => {}
        REB_COLLISION::DIRECT => {
            // Loop over ghost boxes, but only the inner most ring.
            let N_ghost_xcol = if r.N_ghost_x > 1 { 1 } else { r.N_ghost_x };
            let N_ghost_ycol = if r.N_ghost_y > 1 { 1 } else { r.N_ghost_y };
            let N_ghost_zcol = if r.N_ghost_z > 1 { 1 } else { r.N_ghost_z };
            for gbx in -N_ghost_xcol..=N_ghost_xcol {
                for gby in -N_ghost_ycol..=N_ghost_ycol {
                    for gbz in -N_ghost_zcol..=N_ghost_zcol {
                        // Loop over all projectiles
                        for i in 0..N_projectiles {
                            let ip = match &r.map {
                                Some(m) => m[i],
                                None => i,
                            };
                            let p1 = r.particles[ip];
                            let gborig = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                            let mut gb = gborig;
                            // Precalculate shifted position
                            gb.x += p1.x;
                            gb.y += p1.y;
                            gb.z += p1.z;
                            gb.vx += p1.vx;
                            gb.vy += p1.vy;
                            gb.vz += p1.vz;
                            // Loop over all targets
                            for j in 0..N_targets {
                                // Do not collide particle with itself.
                                if i == j {
                                    continue;
                                }
                                let jp = match &r.map {
                                    Some(m) => m[j],
                                    None => j,
                                };
                                let p2 = r.particles[jp];
                                let dx = gb.x - p2.x;
                                let dy = gb.y - p2.y;
                                let dz = gb.z - p2.z;
                                let sr = p1.r + p2.r;
                                let r2 = dx * dx + dy * dy + dz * dz;
                                // Check if particles are overlapping
                                if r2 > sr * sr {
                                    continue;
                                }
                                let dvx = gb.vx - p2.vx;
                                let dvy = gb.vy - p2.vy;
                                let dvz = gb.vz - p2.vz;
                                // Check if particles are approaching each other
                                if dvx * dx + dvy * dy + dvz * dz > 0. {
                                    continue;
                                }
                                r.collisions.push(reb_collision {
                                    p1: ip,
                                    p2: jp,
                                    gb: gborig,
                                    ri: 0,
                                });
                                r.N_collisions += 1;
                            }
                        }
                    }
                }
            }
        }
        REB_COLLISION::LINE => {
            let dt_last_done = r.dt_last_done;
            let N_ghost_xcol = if r.N_ghost_x > 1 { 1 } else { r.N_ghost_x };
            let N_ghost_ycol = if r.N_ghost_y > 1 { 1 } else { r.N_ghost_y };
            let N_ghost_zcol = if r.N_ghost_z > 1 { 1 } else { r.N_ghost_z };
            for gbx in -N_ghost_xcol..=N_ghost_xcol {
                for gby in -N_ghost_ycol..=N_ghost_ycol {
                    for gbz in -N_ghost_zcol..=N_ghost_zcol {
                        for i in 0..N_projectiles {
                            let ip = match &r.map {
                                Some(m) => m[i],
                                None => i,
                            };
                            let p1 = r.particles[ip];
                            let gborig = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                            let mut gb = gborig;
                            gb.x += p1.x;
                            gb.y += p1.y;
                            gb.z += p1.z;
                            gb.vx += p1.vx;
                            gb.vy += p1.vy;
                            gb.vz += p1.vz;
                            for j in (i + 1)..N_projectiles {
                                let jp = match &r.map {
                                    Some(m) => m[j],
                                    None => j,
                                };
                                let p2 = r.particles[jp];
                                let dx1 = gb.x - p2.x; // distance at end
                                let dy1 = gb.y - p2.y;
                                let dz1 = gb.z - p2.z;
                                let r1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
                                let dvx1 = gb.vx - p2.vx;
                                let dvy1 = gb.vy - p2.vy;
                                let dvz1 = gb.vz - p2.vz;
                                let dx2 = dx1 - dt_last_done * dvx1; // distance at beginning
                                let dy2 = dy1 - dt_last_done * dvy1;
                                let dz2 = dz1 - dt_last_done * dvz1;
                                let r2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
                                let t_closest = (dx1 * dvx1 + dy1 * dvy1 + dz1 * dvz1)
                                    / (dvx1 * dvx1 + dvy1 * dvy1 + dvz1 * dvz1);

                                let mut rmin2_ab = r1.min(r2);
                                if t_closest / dt_last_done >= 0. && t_closest / dt_last_done <= 1. {
                                    let dx3 = dx1 - t_closest * dvx1; // closest approach
                                    let dy3 = dy1 - t_closest * dvy1;
                                    let dz3 = dz1 - t_closest * dvz1;
                                    let r3 = dx3 * dx3 + dy3 * dy3 + dz3 * dz3;
                                    rmin2_ab = rmin2_ab.min(r3);
                                }
                                let rsum = p1.r + p2.r;
                                if rmin2_ab > rsum * rsum {
                                    continue;
                                }
                                r.collisions.push(reb_collision {
                                    p1: ip,
                                    p2: jp,
                                    gb: gborig,
                                    ri: 0,
                                });
                                r.N_collisions += 1;
                            }
                        }
                    }
                }
            }
        }
        REB_COLLISION::TREE => {
            // Construct tree
            reb_tree_construct(r);

            let N_ghost_xcol = if r.N_ghost_x > 1 { 1 } else { r.N_ghost_x };
            let N_ghost_ycol = if r.N_ghost_y > 1 { 1 } else { r.N_ghost_y };
            let N_ghost_zcol = if r.N_ghost_z > 1 { 1 } else { r.N_ghost_z };
            // Find second largest radius
            let mut l1 = usize::MAX;
            let mut l2 = usize::MAX;
            reb_simulation_two_largest_particles(r, &mut l1, &mut l2);
            let mut second_largest_radius = 0.;
            if l2 != usize::MAX {
                second_largest_radius = r.particles[l2].r;
            }

            // Loop over all particles
            for i in 0..N_projectiles {
                let p1 = r.particles[i];
                let mut collision_nearest = reb_collision {
                    p1: i,
                    p2: usize::MAX,
                    gb: reb_vec6d::default(),
                    ri: 0,
                };
                let p1_r = p1.r;
                // Loop over ghost boxes.
                for gbx in -N_ghost_xcol..=N_ghost_xcol {
                    for gby in -N_ghost_ycol..=N_ghost_ycol {
                        for gbz in -N_ghost_zcol..=N_ghost_zcol {
                            // Calculated shifted position (for speedup).
                            let mut gb = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                            let gbunmod = gb;
                            gb.x += p1.x;
                            gb.y += p1.y;
                            gb.z += p1.z;
                            gb.vx += p1.vx;
                            gb.vy += p1.vy;
                            gb.vz += p1.vz;
                            // Loop over all root boxes.
                            for ri in 0..N_root {
                                let rootcell = r.tree_root[ri];
                                if rootcell != REB_TREECELL_NONE {
                                    get_nearest_neighbour_in_cell(
                                        r,
                                        &gb,
                                        &gbunmod,
                                        ri,
                                        p1_r,
                                        second_largest_radius,
                                        &mut collision_nearest,
                                        rootcell,
                                    );
                                }
                            }
                        }
                    }
                }
                // Continue if no collision was found
                if collision_nearest.p2 == usize::MAX {
                    continue;
                }
            }
            reb_tree_delete(r);
        }
        REB_COLLISION::LINETREE => {
            // Calculate max drift
            let mut vmax2: f64 = 0.;
            for i in 0..N_projectiles {
                let p1 = r.particles[i];
                vmax2 = vmax2.max(p1.vx * p1.vx + p1.vy * p1.vy + p1.vz * p1.vz);
            }
            let maxdrift = r.dt_last_done * vmax2.sqrt();
            // Construct tree
            reb_tree_construct(r);

            let N_ghost_xcol = if r.N_ghost_x > 1 { 1 } else { r.N_ghost_x };
            let N_ghost_ycol = if r.N_ghost_y > 1 { 1 } else { r.N_ghost_y };
            let N_ghost_zcol = if r.N_ghost_z > 1 { 1 } else { r.N_ghost_z };
            for i in 0..N_projectiles {
                let p1 = r.particles[i];
                let mut collision_nearest = reb_collision {
                    p1: i,
                    p2: usize::MAX,
                    gb: reb_vec6d::default(),
                    ri: 0,
                };
                let p1_r = p1.r;
                // Add drift during last timestep
                let p1_r_plus_dtv = p1_r
                    + r.dt_last_done * (p1.vx * p1.vx + p1.vy * p1.vy + p1.vz * p1.vz).sqrt();
                for gbx in -N_ghost_xcol..=N_ghost_xcol {
                    for gby in -N_ghost_ycol..=N_ghost_ycol {
                        for gbz in -N_ghost_zcol..=N_ghost_zcol {
                            let mut gb = reb_boundary_get_ghostbox(r, gbx, gby, gbz);
                            let gbunmod = gb;
                            gb.x += p1.x;
                            gb.y += p1.y;
                            gb.z += p1.z;
                            gb.vx += p1.vx;
                            gb.vy += p1.vy;
                            gb.vz += p1.vz;
                            for ri in 0..N_root {
                                let rootcell = r.tree_root[ri];
                                if rootcell != REB_TREECELL_NONE {
                                    check_for_overlapping_trajectories_in_cell(
                                        r,
                                        &gb,
                                        &gbunmod,
                                        ri,
                                        p1_r,
                                        p1_r_plus_dtv,
                                        &mut collision_nearest,
                                        rootcell,
                                        maxdrift,
                                    );
                                }
                            }
                        }
                    }
                }
                if collision_nearest.p2 == usize::MAX {
                    continue;
                }
            }
            reb_tree_delete(r);
        }
    }

    // randomize
    for i in 0..r.N_collisions {
        let new = (rand_r(&mut r.rand_seed) as usize) % r.N_collisions;
        let c1 = r.collisions[i];
        r.collisions[i] = r.collisions[new];
        r.collisions[new] = c1;
    }
    // Loop over all collisions previously found in reb_collision_search().

    let resolve = match r.collision_resolve {
        Some(f) => f,
        None => reb_collision_resolve_halt, // Default is to throw an exception
    };

    for i in 0..r.N_collisions {
        let mut c = r.collisions[i];
        if c.p1 != usize::MAX && c.p2 != usize::MAX {
            // Resolve collision
            let outcome = resolve(r, c);

            // Remove particles
            if outcome & REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1 != 0 {
                // Remove p1
                let removedp1 = reb_simulation_remove_particle(r, c.p1) == 0;
                if removedp1 {
                    if c.p2 > c.p1 && c.p2 != usize::MAX {
                        c.p2 -= 1;
                    }
                    for j in (i + 1)..r.N_collisions {
                        // Update other collisions
                        let cp = &mut r.collisions[j];
                        // Skip collisions which involve the removed particle
                        if cp.p1 == c.p1 || cp.p2 == c.p1 {
                            cp.p1 = usize::MAX;
                            cp.p2 = usize::MAX;
                        }
                        // Adjust collisions
                        if cp.p1 > c.p1 && cp.p1 != usize::MAX {
                            cp.p1 -= 1;
                        }
                        if cp.p2 > c.p1 && cp.p2 != usize::MAX {
                            cp.p2 -= 1;
                        }
                    }
                }
            }
            if outcome & REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2 != 0 {
                // Remove p2
                let removedp2 = reb_simulation_remove_particle(r, c.p2) == 0;
                if removedp2 {
                    for j in (i + 1)..r.N_collisions {
                        let cp = &mut r.collisions[j];
                        if cp.p1 == c.p2 || cp.p2 == c.p2 {
                            cp.p1 = usize::MAX;
                            cp.p2 = usize::MAX;
                        }
                        if cp.p1 > c.p2 && cp.p1 != usize::MAX {
                            cp.p1 -= 1;
                        }
                        if cp.p2 > c.p2 && cp.p2 != usize::MAX {
                            cp.p2 -= 1;
                        }
                    }
                }
            }
        }
    }
}

/// collision.c `reb_tree_get_nearest_neighbour_in_cell` (recursive).
fn get_nearest_neighbour_in_cell(
    r: &mut reb_simulation,
    gb: &reb_vec6d,
    gbunmod: &reb_vec6d,
    ri: usize,
    p1_r: f64,
    second_largest_radius: f64,
    collision_nearest: &mut reb_collision,
    c: usize,
) {
    let cell = r.tree_cells[c];
    if cell.pt >= 0 {
        // c is a leaf node
        let condition = cell.pt as usize != collision_nearest.p1;
        if condition {
            let p2 = r.particles[cell.pt as usize];

            let dx = gb.x - p2.x;
            let dy = gb.y - p2.y;
            let dz = gb.z - p2.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            let rp = p1_r + p2.r;
            // reb_particles are not overlapping
            if r2 > rp * rp {
                return;
            }
            let dvx = gb.vx - p2.vx;
            let dvy = gb.vy - p2.vy;
            let dvz = gb.vz - p2.vz;
            // reb_particles are not approaching each other
            if dvx * dx + dvy * dy + dvz * dz > 0. {
                return;
            }
            // Found a new nearest neighbour. Save it for later.
            collision_nearest.ri = ri;
            collision_nearest.p2 = cell.pt as usize;
            collision_nearest.gb = *gbunmod;
            // Save collision in collisions array.
            r.collisions.push(*collision_nearest);
            r.N_collisions += 1;
        }
    } else {
        // c is not a leaf node
        let dx = gb.x - cell.x;
        let dy = gb.y - cell.y;
        let dz = gb.z - cell.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        let rp = p1_r + second_largest_radius + 0.86602540378443 * cell.w;
        // Check if we need to descend into daughter cells
        if r2 < rp * rp {
            for o in 0..8 {
                let d = cell.oct[o];
                if d != REB_TREECELL_NONE {
                    get_nearest_neighbour_in_cell(
                        r,
                        gb,
                        gbunmod,
                        ri,
                        p1_r,
                        second_largest_radius,
                        collision_nearest,
                        d,
                    );
                }
            }
        }
    }
}

/// collision.c `reb_tree_check_for_overlapping_trajectories_in_cell`.
fn check_for_overlapping_trajectories_in_cell(
    r: &mut reb_simulation,
    gb: &reb_vec6d,
    gbunmod: &reb_vec6d,
    ri: usize,
    p1_r: f64,
    p1_r_plus_dtv: f64,
    collision_nearest: &mut reb_collision,
    c: usize,
    maxdrift: f64,
) {
    let cell = r.tree_cells[c];
    if cell.pt >= 0 {
        // c is a leaf node
        if cell.pt as usize != collision_nearest.p1 {
            let p2 = r.particles[cell.pt as usize];
            let dt_done_last = r.dt_last_done;
            let dx1 = gb.x - p2.x; // distance at beginning
            let dy1 = gb.y - p2.y;
            let dz1 = gb.z - p2.z;
            let r1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
            let dvx1 = gb.vx - p2.vx;
            let dvy1 = gb.vy - p2.vy;
            let dvz1 = gb.vz - p2.vz;
            let dx2 = dx1 - dt_done_last * dvx1; // distance at end
            let dy2 = dy1 - dt_done_last * dvy1;
            let dz2 = dz1 - dt_done_last * dvz1;
            let r2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
            let t_closest = (dx1 * dvx1 + dy1 * dvy1 + dz1 * dvz1)
                / (dvx1 * dvx1 + dvy1 * dvy1 + dvz1 * dvz1);

            let mut rmin2_ab = r1.min(r2);
            if t_closest / dt_done_last >= 0. && t_closest / dt_done_last <= 1. {
                let dx3 = dx1 - t_closest * dvx1; // closest approach
                let dy3 = dy1 - t_closest * dvy1;
                let dz3 = dz1 - t_closest * dvz1;
                let r3 = dx3 * dx3 + dy3 * dy3 + dz3 * dz3;
                rmin2_ab = rmin2_ab.min(r3);
            }
            let rsum = p1_r + p2.r;
            if rmin2_ab > rsum * rsum {
                return;
            }
            collision_nearest.ri = ri;
            collision_nearest.p2 = cell.pt as usize;
            collision_nearest.gb = *gbunmod;
            r.collisions.push(*collision_nearest);
            r.N_collisions += 1;
        }
    } else {
        // c is not a leaf node
        let dx = gb.x - cell.x;
        let dy = gb.y - cell.y;
        let dz = gb.z - cell.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        let rp = p1_r_plus_dtv + maxdrift + 0.86602540378443 * cell.w;
        if r2 < rp * rp {
            for o in 0..8 {
                let d = cell.oct[o];
                if d != REB_TREECELL_NONE {
                    check_for_overlapping_trajectories_in_cell(
                        r,
                        gb,
                        gbunmod,
                        ri,
                        p1_r,
                        p1_r_plus_dtv,
                        collision_nearest,
                        d,
                        maxdrift,
                    );
                }
            }
        }
    }
}

/// collision.c `reb_collision_resolve_hardsphere`.
pub fn reb_collision_resolve_hardsphere(
    r: &mut reb_simulation,
    c: reb_collision,
) -> REB_COLLISION_RESOLVE_OUTCOME {
    let p1 = r.particles[c.p1];
    let p2 = r.particles[c.p2];
    let gb = c.gb;
    let x21 = p1.x + gb.x - p2.x;
    let y21 = p1.y + gb.y - p2.y;
    let z21 = p1.z + gb.z - p2.z;
    let rp = p1.r + p2.r;
    let oldvyouter;
    if x21 > 0. {
        oldvyouter = p1.vy;
    } else {
        oldvyouter = p2.vy;
    }
    if rp * rp < x21 * x21 + y21 * y21 + z21 * z21 {
        return 0;
    }
    let vx21 = p1.vx + gb.vx - p2.vx;
    let vy21 = p1.vy + gb.vy - p2.vy;
    let vz21 = p1.vz + gb.vz - p2.vz;
    if vx21 * x21 + vy21 * y21 + vz21 * z21 > 0. {
        return 0; // not approaching
    }
    // Bring the two balls in the xy plane.
    let theta = z21.atan2(y21);
    let stheta = theta.sin();
    let ctheta = theta.cos();
    let vy21n = ctheta * vy21 + stheta * vz21;
    let y21n = ctheta * y21 + stheta * z21;

    // Bring the two balls onto the positive x axis.
    let phi = y21n.atan2(x21);
    let cphi = phi.cos();
    let sphi = phi.sin();
    let vx21nn = cphi * vx21 + sphi * vy21n;

    // Coefficient of restitution
    let mut eps = 1.; // perfect bouncing by default
    if let Some(cor) = r.coefficient_of_restitution {
        eps = cor(r, vx21nn);
    }
    let mut dvx2 = -(1.0 + eps) * vx21nn;
    let minr = if p1.r > p2.r { p2.r } else { p1.r };
    let maxr = if p1.r < p2.r { p2.r } else { p1.r };
    let mut mindv = minr * r.minimum_collision_velocity;
    let _r = (x21 * x21 + y21 * y21 + z21 * z21).sqrt();
    mindv *= 1. - (_r - maxr) / minr;
    if mindv > maxr * r.minimum_collision_velocity {
        mindv = maxr * r.minimum_collision_velocity;
    }
    if dvx2 < mindv {
        dvx2 = mindv;
    }
    // Now we are rotating backwards
    let dvx2n = cphi * dvx2;
    let dvy2n = sphi * dvx2;
    let dvy2nn = ctheta * dvy2n;
    let dvz2nn = stheta * dvy2n;

    // Applying the changes to the particles.
    let p2pf = p1.m / (p1.m + p2.m);
    r.particles[c.p2].vx -= p2pf * dvx2n;
    r.particles[c.p2].vy -= p2pf * dvy2nn;
    r.particles[c.p2].vz -= p2pf * dvz2nn;
    let p1pf = p2.m / (p1.m + p2.m);
    r.particles[c.p1].vx += p1pf * dvx2n;
    r.particles[c.p1].vy += p1pf * dvy2nn;
    r.particles[c.p1].vz += p1pf * dvz2nn;

    // Return y-momentum change
    if x21 > 0. {
        r.collisions_plog += -x21.abs() * (oldvyouter - r.particles[c.p1].vy) * p1.m;
        r.collisions_log_n += 1;
    } else {
        r.collisions_plog += -x21.abs() * (oldvyouter - r.particles[c.p2].vy) * p2.m;
        r.collisions_log_n += 1;
    }
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE
}

/// collision.c `reb_collision_resolve_halt`.
pub fn reb_collision_resolve_halt(
    r: &mut reb_simulation,
    _c: reb_collision,
) -> REB_COLLISION_RESOLVE_OUTCOME {
    r.status = REB_STATUS_COLLISION;
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE
}

/// collision.c `reb_collision_resolve_merge` (serial path).
pub fn reb_collision_resolve_merge(
    r: &mut reb_simulation,
    c: reb_collision,
) -> REB_COLLISION_RESOLVE_OUTCOME {
    // Always remove particle with larger index; merge into lower index.
    let mut swap = false;
    let mut i = c.p1;
    let mut j = c.p2;
    if j < i {
        swap = true;
        i = c.p2;
        j = c.p1;
    }

    let pi = r.particles[i];
    let pj = r.particles[j];

    let invmass = 1.0 / (pi.m + pj.m);

    // Scale out energy from collision - initial energy
    let mut Ei = 0.;
    let mut Ef = 0.;
    if r.track_energy_offset != 0 {
        let vix = pi.vx;
        let viy = pi.vy;
        let viz = pi.vz;
        Ei += 0.5 * pi.m * (vix * vix + viy * viy + viz * viz);
        let vjx = pj.vx;
        let vjy = pj.vy;
        let vjz = pj.vz;
        Ei += 0.5 * pj.m * (vjx * vjx + vjy * vjy + vjz * vjz);
        let N_active = if r.N_active == usize::MAX { r.N } else { r.N_active };
        // No potential energy between test particles
        if i < N_active || j < N_active {
            let x = pi.x - pj.x;
            let y = pi.y - pj.y;
            let z = pi.z - pj.z;
            let _r = (x * x + y * y + z * z).sqrt();

            Ei += -r.G * pi.m * pj.m / _r;
        }
    }

    // Merge by conserving mass, volume and momentum
    {
        let p = &mut r.particles[i];
        p.vx = (pi.vx * pi.m + pj.vx * pj.m) * invmass;
        p.vy = (pi.vy * pi.m + pj.vy * pj.m) * invmass;
        p.vz = (pi.vz * pi.m + pj.vz * pj.m) * invmass;
        p.x = (pi.x * pi.m + pj.x * pj.m) * invmass;
        p.y = (pi.y * pi.m + pj.y * pj.m) * invmass;
        p.z = (pi.z * pi.m + pj.z * pj.m) * invmass;
        p.m = pi.m + pj.m;
        p.r = (pi.r * pi.r * pi.r + pj.r * pj.r * pj.r).cbrt();
    }

    // Keeping track of energy offset
    if r.track_energy_offset != 0 {
        let p = r.particles[i];
        Ef += 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
        r.energy_offset += Ei - Ef;
    }

    if swap {
        REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1
    } else {
        REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2
    }
}
