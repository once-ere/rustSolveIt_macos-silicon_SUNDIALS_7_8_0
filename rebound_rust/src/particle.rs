//! particle.rs — particle add/remove and lookup (from particle.c).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::boundary::reb_boundary_particle_is_in_box;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;

/// particle.c `reb_simulation_add`. The C reallocs the particle array
/// (growth factor 2, start 8) — Vec does the equivalent; the stored
/// values are identical.
pub fn reb_simulation_add(r: &mut reb_simulation, pt: reb_particle) {
    if !reb_boundary_particle_is_in_box(r, pt) {
        reb_simulation_error(r, "Particle outside of box boundaries. Did not add particle.");
        return;
    }
    r.particles.push(pt);
    r.N += 1;

    // Check if any integrators need to do extra work
    crate::simulation::reb_integrator_did_add_particle(r);

    r.did_modify_particles = 1;
}

/// particle.c `reb_particle_cmp` (pointer `name` comparison becomes
/// index comparison). Returns true if the particles differ.
pub fn reb_particle_cmp(p1: reb_particle, p2: reb_particle) -> bool {
    let mut differ = false;
    differ = differ || (p1.x != p2.x);
    differ = differ || (p1.y != p2.y);
    differ = differ || (p1.z != p2.z);
    differ = differ || (p1.vx != p2.vx);
    differ = differ || (p1.vy != p2.vy);
    differ = differ || (p1.vz != p2.vz);
    differ = differ || (p1.ax != p2.ax);
    differ = differ || (p1.ay != p2.ay);
    differ = differ || (p1.az != p2.az);
    differ = differ || (p1.m != p2.m);
    differ = differ || (p1.r != p2.r);
    differ = differ || (p1.name != p2.name);
    differ
}

/// particle.c `reb_particle_check_testparticles`.
pub fn reb_particle_check_testparticles(r: &reb_simulation) -> bool {
    if r.N_active == r.N || r.N_active == usize::MAX {
        return false;
    }
    if r.testparticle_type == 0 {
        let mut found_issue = false;
        for i in r.N_active..r.N {
            if r.particles[i].m != 0. {
                found_issue = true;
            }
        }
        if found_issue {
            return true;
        }
    }
    false
}

/// particle.c `reb_get_rootbox_for_particle`. C integer semantics:
/// `((int)floor(...) + N_root) % N_root` with truncating `%` — for the
/// values reachable in-box the operand is non-negative, so plain `%`
/// on i64 is exact.
pub fn reb_get_rootbox_for_particle(r: &reb_simulation, pt: reb_particle) -> i32 {
    if r.root_size == -1. {
        return 0;
    }
    let nx = r.N_root_x as i64;
    let ny = r.N_root_y as i64;
    let nz = r.N_root_z as i64;
    let i = ((((pt.x + r.root_size * (nx as f64) / 2.) / r.root_size).floor() as i64) + nx) % nx;
    let j = ((((pt.y + r.root_size * (ny as f64) / 2.) / r.root_size).floor() as i64) + ny) % ny;
    let k = ((((pt.z + r.root_size * (nz as f64) / 2.) / r.root_size).floor() as i64) + nz) % nz;
    let index = (k * ny + j) * nx + i;
    index as i32
}

/// particle.c `reb_simulation_register_name` (interning into the
/// simulation's name list; returns the index used by
/// `reb_particle::name`).
pub fn reb_simulation_register_name(r: &mut reb_simulation, name: &str) -> usize {
    for (i, n) in r.name_list.iter().enumerate() {
        if n == name {
            return i;
        }
    }
    r.name_list.push(name.to_string());
    r.name_list.len() - 1
}

/// particle.c `reb_particle_set_name` — takes (simulation, particle
/// index) instead of the C particle pointer with `sim` back-pointer.
pub fn reb_particle_set_name(r: &mut reb_simulation, index: usize, name: Option<&str>) {
    match name {
        None => r.particles[index].name = None,
        Some(n) => {
            let id = reb_simulation_register_name(r, n);
            r.particles[index].name = Some(id);
        }
    }
}

/// particle.c `reb_simulation_get_particle_by_name` — returns the
/// particle's index (the C returns a pointer).
pub fn reb_simulation_get_particle_by_name(r: &reb_simulation, name: &str) -> Option<usize> {
    for i in 0..r.N {
        if let Some(id) = r.particles[i].name {
            if r.name_list[id] == name {
                return Some(i);
            }
        }
    }
    None
}

/// particle.c `reb_simulation_remove_all_particles`.
pub fn reb_simulation_remove_all_particles(r: &mut reb_simulation) {
    r.N = 0;
    r.N_active = usize::MAX;
    r.N_var = 0;
    r.particles.clear();
    r.particles_var.clear();
}

/// particle.c `reb_simulation_remove_particle`. Returns 0 on success,
/// 1 on failure (C convention).
pub fn reb_simulation_remove_particle(r: &mut reb_simulation, index: usize) -> i32 {
    if r.N_var != 0 {
        reb_simulation_error(
            r,
            "Removing particles not supported when variational particles are in use. Did not remove particle.",
        );
        return 1;
    }
    crate::simulation::reb_integrator_will_remove_particle(r, index);

    if r.N == 1 {
        r.N = 0;
        r.particles.clear();
        r.did_modify_particles = 1;
        reb_simulation_warning(r, "Last particle removed.");
        return 0;
    }
    if index >= r.N {
        let warning = format!(
            "Index {} passed to particles_remove was out of range (N={}).  Did not remove particle.",
            index, r.N
        );
        reb_simulation_error(r, &warning);
        return 1;
    }
    r.N -= 1;
    if index < r.N_active && r.N_active != usize::MAX {
        r.N_active -= 1;
    }
    for j in index..r.N {
        r.particles[j] = r.particles[j + 1];
    }
    r.particles.truncate(r.N);

    r.did_modify_particles = 1;
    0
}

/// particle.c `reb_simulation_remove_particle_by_name`.
pub fn reb_simulation_remove_particle_by_name(r: &mut reb_simulation, name: &str) -> i32 {
    match reb_simulation_get_particle_by_name(r, name) {
        None => {
            reb_simulation_error(r, "Particle not found.");
            1
        }
        Some(index) => reb_simulation_remove_particle(r, index),
    }
}

/// particle.c `reb_particle_isub`.
pub fn reb_particle_isub(p1: &mut reb_particle, p2: &reb_particle) {
    p1.x -= p2.x;
    p1.y -= p2.y;
    p1.z -= p2.z;
    p1.vx -= p2.vx;
    p1.vy -= p2.vy;
    p1.vz -= p2.vz;
    p1.m -= p2.m;
}

/// particle.c `reb_particle_iadd`.
pub fn reb_particle_iadd(p1: &mut reb_particle, p2: &reb_particle) {
    p1.x += p2.x;
    p1.y += p2.y;
    p1.z += p2.z;
    p1.vx += p2.vx;
    p1.vy += p2.vy;
    p1.vz += p2.vz;
    p1.m += p2.m;
}

/// particle.c `reb_particle_imul`.
pub fn reb_particle_imul(p1: &mut reb_particle, value: f64) {
    p1.x *= value;
    p1.y *= value;
    p1.z *= value;
    p1.vx *= value;
    p1.vy *= value;
    p1.vz *= value;
    p1.m *= value;
}

/// particle.c `reb_particle_distance`.
pub fn reb_particle_distance(p1: &reb_particle, p2: &reb_particle) -> f64 {
    let dx = p1.x - p2.x;
    let dy = p1.y - p2.y;
    let dz = p1.z - p2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// simulation.c `reb_simulation_two_largest_particles` (serial branch).
pub fn reb_simulation_two_largest_particles(r: &reb_simulation, p1: &mut usize, p2: &mut usize) {
    *p1 = usize::MAX;
    *p2 = usize::MAX;
    let mut largest1 = -1.0;
    let mut largest2 = -1.0;
    for i in 0..r.N {
        if r.particles[i].r > largest1 {
            largest2 = largest1;
            *p2 = *p1;
            largest1 = r.particles[i].r;
            *p1 = i;
        } else if r.particles[i].r > largest2 {
            largest2 = r.particles[i].r;
            *p2 = i;
        }
    }
}
