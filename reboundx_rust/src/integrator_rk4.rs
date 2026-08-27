//! integrator_rk4.rs — translation of REBOUNDx integrator_rk4.c
//! 4th order Runge Kutta method, used by `integrate_force` (and by any
//! effect that asks for `REBX_INTEGRATOR_RK4`) to advance the velocities
//! under a single REBOUNDx force across one timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Parameters
//!
//! This integrator takes no user parameters. It allocates two scratch
//! `struct reb_particle` buffers on its first call and keeps them on the
//! force's parameter list (C: `force->ap`) between calls:
//!
//! | Name      | Type      | Required | Description                                        |
//! |-----------|-----------|----------|----------------------------------------------------|
//! | `rk4_k2`  | pointer   | internal | Scratch particles for RK4 stages 2 and 4           |
//! | `rk4_k3`  | pointer   | internal | Scratch particles for RK4 stage 3, then k2+k3 sums  |
//!
//! Both names are registered as `REBX_TYPE_POINTER` by
//! [`crate::core::rebx_attach`]; here their payload is a
//! `rebx_param_value::particles` (an owned `Vec<reb_particle>`).
//!
//! # Deviations from the C
//!
//! 1. The C calls `force->update_accelerations(sim, force, ps, N)` with
//!    `ps` pointing either at `sim->particles` or at one of the scratch
//!    buffers. The translated [`crate::types::rebx_force_fn`] has no
//!    particles argument — every force works on `sim.particles` — so a
//!    stage that must be evaluated on a scratch buffer swaps that buffer
//!    into `sim.particles` for the duration of the call and swaps it
//!    straight back. The force therefore sees exactly the array the C
//!    hands it, element for element, and per-particle parameters still
//!    resolve because `rebx_ap::particle(i)` indexes the same `i` (the
//!    C's `particles[i].ap` is likewise memcpy'd from particle `i`).
//! 2. `free()` has no counterpart: the buffers are owned by the
//!    parameter list and are dropped with the extras. The C's
//!    `force->free_memory` hook does not exist on `rebx_force`, so
//!    [`rebx_rk4_free_memory`] is provided for the same effect (it
//!    releases both buffers) but nothing installs it.
//! 3. The C dereferences `force->update_accelerations` unconditionally;
//!    here it is an `Option` (`rebx_add_force` refuses a force without
//!    one), and a `None` simply performs no stage evaluation rather than
//!    crashing.
//! 4. The scratch buffers are sized to `sim.N` on every call. The C
//!    mallocs `N*sizeof(struct reb_particle)` once and would run off the
//!    end of that allocation if particles were added afterwards; the
//!    resize is a no-op whenever `N` has not changed, and the buffers'
//!    first `N` entries are overwritten from `sim.particles` immediately
//!    afterwards, so no value the C computes is affected.

use rebound_rs::{reb_particle, reb_simulation};

use crate::core::{rebx_get_param_particles, rebx_reset_accelerations, rebx_set_param_particles};
use crate::types::{rebx_ap, rebx_extras, rebx_param_value};

/// integrator_rk4.c `rebx_rk4_free_memory`.
///
/// The C looks both scratch buffers up on the force's parameter list and
/// `free()`s them; it is installed as `force->free_memory` and called
/// from `rebx_free`. Here the buffers are owned `Vec`s, so releasing
/// them means putting the parameters back into the "registered but never
/// set" state (C: `param->value == NULL`) and letting the `Vec`s drop.
/// Nothing calls this in the Rust port — dropping the extras frees the
/// buffers — but it is kept so the C function has a counterpart.
pub fn rebx_rk4_free_memory(rebx: &mut rebx_extras, force: usize) {
    rebx_rk4_release_buffer(rebx, force, "rk4_k2");
    rebx_rk4_release_buffer(rebx, force, "rk4_k3");
}

/// Drop the `Vec` held by one scratch parameter, leaving the node in
/// place with no value (C: `free(ptr)`).
fn rebx_rk4_release_buffer(rebx: &mut rebx_extras, force: usize, name: &str) {
    let ap = match rebx.ap_mut(rebx_ap::force(force)) {
        Some(ap) => ap,
        None => return,
    };
    for param in ap.iter_mut() {
        if param.name == name {
            param.value = rebx_param_value::none;
            return;
        }
    }
}

/// Move one scratch buffer out of the force's parameter list so that it
/// can be handed to the force functions alongside `rebx` (C: reading the
/// `struct reb_particle*` out of the parameter — the C keeps using the
/// list's copy of the pointer, which safe Rust cannot alias).
///
/// The traversal is `rebx_get_param_struct`'s: from index 0, the C list
/// head, taking the first name match. The node keeps its place in the
/// list, so [`rebx_rk4_store_buffer`] puts the buffer back without
/// changing the list order.
fn rebx_rk4_take_buffer(rebx: &mut rebx_extras, force: usize, name: &str) -> Vec<reb_particle> {
    let ap = match rebx.ap_mut(rebx_ap::force(force)) {
        Some(ap) => ap,
        None => return Vec::new(),
    };
    for param in ap.iter_mut() {
        if param.name == name {
            if let rebx_param_value::particles(buf) = &mut param.value {
                return std::mem::take(buf);
            }
            return Vec::new();
        }
    }
    Vec::new()
}

/// Put back what [`rebx_rk4_take_buffer`] removed. `rebx_set_param_*`
/// reuses the existing node when the name is already on the list, so the
/// list — and therefore every later traversal order — is unchanged.
fn rebx_rk4_store_buffer(
    rebx: &mut rebx_extras,
    force: usize,
    name: &str,
    buf: Vec<reb_particle>,
) {
    rebx_set_param_particles(rebx, rebx_ap::force(force), name, buf);
}

/// C: `force->update_accelerations(sim, force, sim->particles, N)`.
fn rebx_rk4_update_accelerations(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force: usize,
    N: usize,
) {
    // Read the function pointer out in its own statement: it must not
    // still be borrowed from `rebx` when `rebx` is handed to it.
    let update_accelerations = rebx.allocated_forces[force].update_accelerations;
    if let Some(update_accelerations) = update_accelerations {
        update_accelerations(sim, rebx, force, N);
    }
}

/// C: `force->update_accelerations(sim, force, ps, N)` where `ps` is one
/// of the scratch buffers rather than `sim->particles`.
///
/// See deviation 1 in the module docs: the buffer takes the place of
/// `sim.particles` for exactly the duration of the call.
fn rebx_rk4_update_accelerations_on(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force: usize,
    N: usize,
    ps: &mut Vec<reb_particle>,
) {
    std::mem::swap(&mut sim.particles, ps);
    rebx_rk4_update_accelerations(sim, rebx, force, N);
    std::mem::swap(&mut sim.particles, ps);
}

/// integrator_rk4.c `rebx_integrator_rk4_integrate`.
///
/// Classic RK4 on the velocities only: the force is evaluated four
/// times (on `sim.particles`, then on `k2`, `k3` and — reusing `k2` —
/// `k4`), and the four accelerations are combined with weights
/// `dt/6 * (k1 + k4 + 2*(k2 + k3))`. As in the C, `k3` accumulates
/// `k2 + k3` in place so that the last stage can reuse the `k2` buffer
/// without a memcpy.
///
/// C signature: `(struct reb_simulation* sim, const double dt,
/// struct rebx_force* force)`; `rebx` is explicit here and `force` is an
/// index into `rebx.allocated_forces`.
pub fn rebx_integrator_rk4_integrate(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    dt: f64,
    force: usize,
) {
    let N = sim.N;
    rebx_reset_accelerations(&mut sim.particles, N);

    // C: k2 = rebx_get_param(rebx, force->ap, "rk4_k2");
    //    if (k2 == NULL){ malloc both, set both, install free_memory }
    // The two setters run in the C's order, so the parameters land on
    // the force's list in the C's order too.
    if rebx_get_param_particles(rebx, rebx_ap::force(force), "rk4_k2").is_none() {
        let k2 = vec![reb_particle::default(); N];
        let k3 = vec![reb_particle::default(); N];
        rebx_set_param_particles(rebx, rebx_ap::force(force), "rk4_k2", k2);
        rebx_set_param_particles(rebx, rebx_ap::force(force), "rk4_k3", k3);
        // C: force->free_memory = rebx_rk4_free_memory;
        // No `free_memory` slot exists here; see deviation 2.
    }

    let mut k2 = rebx_rk4_take_buffer(rebx, force, "rk4_k2");
    let mut k3 = rebx_rk4_take_buffer(rebx, force, "rk4_k3");
    // Deviation 4: a no-op unless sim.N changed since the buffers were
    // allocated. Every one of the N entries is overwritten just below.
    k2.resize(N, reb_particle::default());
    k3.resize(N, reb_particle::default());

    // memcpy(k2, sim->particles, N*sizeof(*k2));
    k2[..N].copy_from_slice(&sim.particles[..N]);
    // memcpy(k3, sim->particles, N*sizeof(*k3));
    k3[..N].copy_from_slice(&sim.particles[..N]);

    let dt2 = dt / 2.;
    rebx_rk4_update_accelerations(sim, rebx, force, N); // k1 = sim.particles.a

    for i in 0..N {
        k2[i].vx = sim.particles[i].vx + dt2 * sim.particles[i].ax;
        k2[i].vy = sim.particles[i].vy + dt2 * sim.particles[i].ay;
        k2[i].vz = sim.particles[i].vz + dt2 * sim.particles[i].az;
    }
    rebx_rk4_update_accelerations_on(sim, rebx, force, N, &mut k2);

    for i in 0..N {
        k3[i].vx = sim.particles[i].vx + dt2 * k2[i].ax;
        k3[i].vy = sim.particles[i].vy + dt2 * k2[i].ay;
        k3[i].vz = sim.particles[i].vz + dt2 * k2[i].az;
    }
    rebx_rk4_update_accelerations_on(sim, rebx, force, N, &mut k3);

    for i in 0..N {
        // store k2+k3 in k3 and reuse k2 for k4 to avoid a memcpy
        k2[i].vx = sim.particles[i].vx + dt * k3[i].ax;
        k2[i].vy = sim.particles[i].vy + dt * k3[i].ay;
        k2[i].vz = sim.particles[i].vz + dt * k3[i].az;
        k3[i].ax += k2[i].ax;
        k3[i].ay += k2[i].ay;
        k3[i].az += k2[i].az;
    }
    rebx_reset_accelerations(&mut k2, N);
    rebx_rk4_update_accelerations_on(sim, rebx, force, N, &mut k2);

    let dt6 = dt / 6.;
    for i in 0..N {
        sim.particles[i].vx += dt6 * (sim.particles[i].ax + k2[i].ax + 2. * k3[i].ax);
        sim.particles[i].vy += dt6 * (sim.particles[i].ay + k2[i].ay + 2. * k3[i].ay);
        sim.particles[i].vz += dt6 * (sim.particles[i].az + k2[i].az + 2. * k3[i].az);
    }

    // The C never gives its pointers back because it never took them
    // away; here the buffers return to the force's parameter list so the
    // next call finds them (and `rebx_get_param_particles` above sees
    // them as already allocated).
    rebx_rk4_store_buffer(rebx, force, "rk4_k2", k2);
    rebx_rk4_store_buffer(rebx, force, "rk4_k3", k3);
}
