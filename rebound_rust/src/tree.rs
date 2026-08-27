//! tree.rs — octree construction and tree-walk force evaluation
//! (from tree.c). The C allocates each `reb_treecell` with calloc and
//! links them with pointers; here the cells live in an index arena
//! (`r.tree_cells`) that is cleared when the tree is deleted — the
//! insertion order, cell geometry and traversal order are identical.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::particle::reb_get_rootbox_for_particle;
use crate::tools::reb_simulation_error;
use crate::types::*;

/// tree.c `reb_reb_tree_get_octant_for_particle_in_cell`.
fn get_octant_for_particle_in_cell(p: &reb_particle, node: &reb_treecell) -> usize {
    let mut octant = 0;
    if p.x < node.x {
        octant += 1;
    }
    if p.y < node.y {
        octant += 2;
    }
    if p.z < node.z {
        octant += 4;
    }
    octant
}

fn new_cell(r: &mut reb_simulation, cell: reb_treecell) -> usize {
    r.tree_cells.push(cell);
    r.tree_cells.len() - 1
}

/// tree.c `reb_tree_add_particle_to_cell` (recursive; `node` and
/// `parent` are arena indices, REB_TREECELL_NONE = C NULL).
fn add_particle_to_cell(
    r: &mut reb_simulation,
    node: usize,
    pt: usize,
    parent: usize,
    o: usize,
) -> usize {
    // Initialize a new node
    if node == REB_TREECELL_NONE {
        let p = r.particles[pt];
        let mut cell = reb_treecell {
            x: 0.,
            y: 0.,
            z: 0.,
            w: 0.,
            m: 0.,
            mx: 0.,
            my: 0.,
            mz: 0.,
            oct: [REB_TREECELL_NONE; 8],
            pt: 0,
            remote: 0,
        };
        if parent == REB_TREECELL_NONE {
            // The new node is a root
            cell.w = r.root_size;
            let boxsize = reb_vec3d {
                x: r.root_size * (r.N_root_x as f64),
                y: r.root_size * (r.N_root_y as f64),
                z: r.root_size * (r.N_root_z as f64),
            };
            let i = (((p.x + boxsize.x / 2.) / r.root_size).floor() as i64)
                .rem_euclid(r.N_root_x as i64) as f64;
            let j = (((p.y + boxsize.y / 2.) / r.root_size).floor() as i64)
                .rem_euclid(r.N_root_y as i64) as f64;
            let k = (((p.z + boxsize.z / 2.) / r.root_size).floor() as i64)
                .rem_euclid(r.N_root_z as i64) as f64;
            cell.x = -boxsize.x / 2. + r.root_size * (0.5 + i);
            cell.y = -boxsize.y / 2. + r.root_size * (0.5 + j);
            cell.z = -boxsize.z / 2. + r.root_size * (0.5 + k);
        } else {
            // The new node is a normal node
            let pw = r.tree_cells[parent].w;
            let px = r.tree_cells[parent].x;
            let py = r.tree_cells[parent].y;
            let pz = r.tree_cells[parent].z;
            cell.w = pw / 2.;
            cell.x = px + cell.w / 2. * if (o >> 0) % 2 == 0 { 1. } else { -1. };
            cell.y = py + cell.w / 2. * if (o >> 1) % 2 == 0 { 1. } else { -1. };
            cell.z = pz + cell.w / 2. * if (o >> 2) % 2 == 0 { 1. } else { -1. };
        }
        if cell.w <= 0.0 {
            reb_simulation_error(r, "Tree cell has size zero.");
            return REB_TREECELL_NONE;
        }
        cell.pt = pt as i32;
        return new_cell(r, cell);
    }
    // In an existing node
    if r.tree_cells[node].pt >= 0 {
        // It's a leaf node
        let node_pt = r.tree_cells[node].pt as usize;
        let o1 = get_octant_for_particle_in_cell(&r.particles[node_pt], &r.tree_cells[node]);
        let o2 = get_octant_for_particle_in_cell(&r.particles[pt], &r.tree_cells[node]);
        if o1 == o2 {
            // Same octant: check same coordinates to avoid infinite recursion
            if r.particles[pt].x == r.particles[node_pt].x
                && r.particles[pt].y == r.particles[node_pt].y
                && r.particles[pt].z == r.particles[node_pt].z
            {
                reb_simulation_error(
                    r,
                    "Cannot add two particles with the same coordinates to the tree.",
                );
                return node;
            }
        }
        let c1 = add_particle_to_cell(r, r.tree_cells[node].oct[o1], node_pt, node, o1);
        r.tree_cells[node].oct[o1] = c1;
        let c2 = add_particle_to_cell(r, r.tree_cells[node].oct[o2], pt, node, o2);
        r.tree_cells[node].oct[o2] = c2;
        r.tree_cells[node].pt = -2;
    } else {
        // It's not a leaf
        r.tree_cells[node].pt -= 1;
        let o = get_octant_for_particle_in_cell(&r.particles[pt], &r.tree_cells[node]);
        let c = add_particle_to_cell(r, r.tree_cells[node].oct[o], pt, node, o);
        r.tree_cells[node].oct[o] = c;
    }
    node
}

/// tree.c `reb_tree_add_particle_to_tree`.
pub fn reb_tree_add_particle_to_tree(r: &mut reb_simulation, pt: usize) {
    let p = r.particles[pt];
    if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
        reb_simulation_error(r, "Particle has non-finite coordinates. Cannot add to tree.");
        return;
    }
    let rootbox = reb_get_rootbox_for_particle(r, p) as usize;
    let root = add_particle_to_cell(r, r.tree_root[rootbox], pt, REB_TREECELL_NONE, 0);
    r.tree_root[rootbox] = root;
}

/// tree.c `reb_tree_calculate_gravity_data_in_cell` (recursive).
fn calculate_gravity_data_in_cell(r: &mut reb_simulation, node: usize) {
    if r.tree_cells[node].pt < 0 {
        // Non-leaf nodes
        r.tree_cells[node].m = 0.;
        r.tree_cells[node].mx = 0.;
        r.tree_cells[node].my = 0.;
        r.tree_cells[node].mz = 0.;
        for o in 0..8 {
            let d = r.tree_cells[node].oct[o];
            if d != REB_TREECELL_NONE {
                calculate_gravity_data_in_cell(r, d);
                // Calculate the total mass and the center of mass
                let d_m = r.tree_cells[d].m;
                let dmx = r.tree_cells[d].mx;
                let dmy = r.tree_cells[d].my;
                let dmz = r.tree_cells[d].mz;
                r.tree_cells[node].mx += dmx * d_m;
                r.tree_cells[node].my += dmy * d_m;
                r.tree_cells[node].mz += dmz * d_m;
                r.tree_cells[node].m += d_m;
            }
        }
        let m_tot = r.tree_cells[node].m;
        if m_tot > 0. {
            r.tree_cells[node].mx /= m_tot;
            r.tree_cells[node].my /= m_tot;
            r.tree_cells[node].mz /= m_tot;
        }
    } else {
        // Leaf nodes
        let p = r.particles[r.tree_cells[node].pt as usize];
        r.tree_cells[node].m = p.m;
        r.tree_cells[node].mx = p.x;
        r.tree_cells[node].my = p.y;
        r.tree_cells[node].mz = p.z;
    }
}

/// tree.c `reb_tree_calculate_gravity_data`.
pub fn reb_tree_calculate_gravity_data(r: &mut reb_simulation) {
    let N_root = r.N_root_x * r.N_root_y * r.N_root_z;
    for i in 0..N_root {
        if !r.tree_root.is_empty() && r.tree_root[i] != REB_TREECELL_NONE {
            let root = r.tree_root[i];
            calculate_gravity_data_in_cell(r, root);
        }
    }
}

/// tree.c `reb_tree_delete` (frees all cells; the arena is cleared).
pub fn reb_tree_delete(r: &mut reb_simulation) {
    if !r.tree_root.is_empty() {
        let N_root = r.N_root_x * r.N_root_y * r.N_root_z;
        for i in 0..N_root {
            r.tree_root[i] = REB_TREECELL_NONE;
        }
    }
    r.tree_cells.clear();
}

/// tree.c `reb_tree_construct`.
pub fn reb_tree_construct(r: &mut reb_simulation) {
    if r.root_size <= 0.0 {
        reb_simulation_error(
            r,
            "Set root_size to a finite value to use a tree based gravity or collision solver.",
        );
        return;
    }
    if r.tree_root.is_empty() {
        let N_root = r.N_root_x * r.N_root_y * r.N_root_z;
        r.tree_root = vec![REB_TREECELL_NONE; N_root];
    }
    for i in 0..r.N {
        let p = r.particles[i];
        if p.x.abs() > r.root_size * (r.N_root_x as f64) / 2.
            || p.y.abs() > r.root_size * (r.N_root_y as f64) / 2.
            || p.z.abs() > r.root_size * (r.N_root_z as f64) / 2.
        {
            reb_simulation_error(r, "Particle is outside of simulation box. Cannot add to tree.");
            return;
        }
        reb_tree_add_particle_to_tree(r, i);
    }
}

/// tree.c `reb_tree_calculate_acceleration_for_particle_from_cell`
/// (recursive tree walk with the opening-angle criterion).
fn calculate_acceleration_for_particle_from_cell(
    r: &mut reb_simulation,
    pt: usize,
    node: usize,
    gb: &reb_vec6d,
) {
    let G = r.G;
    let softening2 = r.softening * r.softening;
    let cell = r.tree_cells[node];
    let dx = gb.x - cell.mx;
    let dy = gb.y - cell.my;
    let dz = gb.z - cell.mz;
    let r2 = dx * dx + dy * dy + dz * dz;
    if cell.pt < 0 {
        // Not a leaf
        if cell.w * cell.w > r.opening_angle2 * r2 {
            for o in 0..8 {
                if cell.oct[o] != REB_TREECELL_NONE {
                    calculate_acceleration_for_particle_from_cell(r, pt, cell.oct[o], gb);
                }
            }
        } else {
            let _r = (r2 + softening2).sqrt();
            let prefact = -G / (_r * _r * _r) * cell.m;
            r.particles[pt].ax += prefact * dx;
            r.particles[pt].ay += prefact * dy;
            r.particles[pt].az += prefact * dz;
        }
    } else {
        // It's a leaf node
        if cell.remote == 0 && cell.pt as usize == pt {
            return;
        }
        let _r = (r2 + softening2).sqrt();
        let prefact = -G / (_r * _r * _r) * cell.m;
        r.particles[pt].ax += prefact * dx;
        r.particles[pt].ay += prefact * dy;
        r.particles[pt].az += prefact * dz;
    }
}

/// tree.c `reb_tree_calculate_acceleration_for_particle`.
pub fn reb_tree_calculate_acceleration_for_particle(
    r: &mut reb_simulation,
    pt: usize,
    gb: &reb_vec6d,
) {
    let N_root = r.N_root_x * r.N_root_y * r.N_root_z;
    for i in 0..N_root {
        let node = r.tree_root[i];
        if node != REB_TREECELL_NONE {
            calculate_acceleration_for_particle_from_cell(r, pt, node, gb);
        }
    }
}
