//! integrator_rk2.rs — translation of REBOUNDx integrator_rk2.c
//! 2nd order Runge Kutta method (Ralston's method), used to integrate a
//! single REBOUNDx force across a timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Parameters
//!
//! This file implements an integrator, not an effect, so there are no
//! user-facing parameters. It keeps one *internal* parameter on the
//! force it is integrating:
//!
//! Name (key)   | Type            | Description
//! ------------ | --------------- | -----------------------------------
//! rk2_k2       | REBX_TYPE_POINTER | Scratch buffer of `sim.N` particles holding the second Runge-Kutta stage. Allocated on first use and reused afterwards; not meant to be set or read by users.
//!
//! Ralston's method advances the velocities with the two stage
//! evaluations of the force at weights `b1 = dt/4` and `b2 = 3 dt/4`,
//! the second stage being evaluated at `a21 = 2 dt/3`.
//!
//! # Deviations from the C
//!
//! * The C keeps `k2` as a `malloc`'d `struct reb_particle*` stored in a
//!   `REBX_TYPE_POINTER` parameter, and registers
//!   `rebx_rk2_free_memory` in `force->free_memory` so that `rebx_free`
//!   releases it. Here the buffer is a `Vec<reb_particle>` owned by the
//!   parameter (`rebx_param_value::particles`), so it is released when
//!   the REBOUNDx state is dropped; `rebx_force` accordingly has no
//!   `free_memory` slot. `rebx_rk2_free_memory` is still provided below
//!   and does exactly what the C's does — release the buffer.
//! * The C hands the scratch buffer to the force as the `particles`
//!   argument (`force->update_accelerations(sim, force, k2, N)`). The
//!   Rust force signature has no `particles` argument — a force always
//!   acts on `sim.particles` — so the buffer is swapped into
//!   `sim.particles` for the duration of that call and swapped back
//!   afterwards. The force therefore sees exactly the same particle
//!   state, at the same indices, as the C does.
//! * The buffer is moved out of the parameter list for the duration of
//!   the step and stored back at the end (the crate's take/use/put
//!   pattern), because safe Rust cannot hold `&mut sim` and a `&mut`
//!   into `rebx` at the same time. On the first call this makes the
//!   parameter appear at the end of the step rather than at its start;
//!   names in an `ap` list are unique, so nothing observes the
//!   difference.
//! * If `sim.N` has grown since the buffer was allocated, it is grown to
//!   match. The C would `memcpy` past the end of its `malloc`.

use rebound_rs::{reb_particle, reb_simulation};

use crate::core::rebx_set_param_particles;
use crate::types::{rebx_ap, rebx_extras, rebx_param_value};

/// Move the `"rk2_k2"` buffer out of the force's parameter list,
/// leaving the parameter behind with no value.
///
/// Returns `None` in exactly the cases where the C's
/// `rebx_get_param(rebx, force->ap, "rk2_k2")` returns `NULL`: the
/// parameter is absent, or it exists but was never given a value.
fn rebx_rk2_take_k2(rebx: &mut rebx_extras, force_idx: usize) -> Option<Vec<reb_particle>> {
    let ap = rebx.ap_mut(rebx_ap::force(force_idx))?;
    let param = ap.iter_mut().find(|param| param.name == "rk2_k2")?;
    match std::mem::replace(&mut param.value, rebx_param_value::none) {
        rebx_param_value::particles(k2) => Some(k2),
        // Some other payload under this name: put it back untouched.
        other => {
            param.value = other;
            None
        }
    }
}

/// integrator_rk2.c `rebx_rk2_free_memory`.
///
/// The C installs this in `force->free_memory` and it runs `free(k2)`.
/// There is no `free_memory` slot here (see the module docs): the buffer
/// dies with the REBOUNDx state. Calling this releases it early, exactly
/// as the C's `free` does; the next integration step reallocates it.
pub fn rebx_rk2_free_memory(rebx: &mut rebx_extras, force_idx: usize) {
    let k2 = rebx_rk2_take_k2(rebx, force_idx);
    drop(k2); // C: free(k2);
}

/// integrator_rk2.c `rebx_integrator_rk2_integrate`.
///
/// C: `rebx_integrator_rk2_integrate(sim, dt, force)`. `force_idx`
/// indexes `rebx.allocated_forces`, and the REBOUNDx state the C reaches
/// through `sim->extras` is passed explicitly.
pub fn rebx_integrator_rk2_integrate(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    dt: f64,
    force_idx: usize,
) {
    let N = sim.N;
    let mut k2 = match rebx_rk2_take_k2(rebx, force_idx) {
        Some(k2) => k2,
        None => {
            // C: k2 = malloc(N*sizeof(*k2)); the contents are
            // uninitialized there and overwritten by the memcpy below,
            // so zeroed particles are equivalent.
            // C also sets force->free_memory = rebx_rk2_free_memory here.
            vec![reb_particle::default(); N]
        }
    };
    if k2.len() < N {
        // Not in the C: it would memcpy past the end of the buffer it
        // allocated when sim->N was smaller.
        k2.resize(N, reb_particle::default());
    }
    // C: memcpy(k2, sim->particles, N*sizeof(*k2));
    k2[..N].copy_from_slice(&sim.particles[..N]);

    // C: force->update_accelerations(sim, force, sim->particles, N);
    // Read the function pointer out first: it must not still be borrowed
    // from `rebx` when `rebx` is handed to it.
    let update_accelerations = rebx.allocated_forces[force_idx].update_accelerations;
    if let Some(update_accelerations) = update_accelerations {
        update_accelerations(sim, rebx, force_idx, N);
    }
    let a21 = 2. * dt / 3.;
    for i in 0..N {
        k2[i].vx = sim.particles[i].vx + a21 * sim.particles[i].ax;
        k2[i].vy = sim.particles[i].vy + a21 * sim.particles[i].ay;
        k2[i].vz = sim.particles[i].vz + a21 * sim.particles[i].az;
    }

    // C: force->update_accelerations(sim, force, k2, N);
    // The force acts on sim.particles here, so k2 takes that place for
    // the call and the real particles are put back straight after. Note
    // that, as in the C, the k2 accelerations are NOT reset first: they
    // are the values copied from sim->particles by the memcpy above
    // (zero, since rebx_integrate_force calls rebx_reset_accelerations
    // before invoking the integrator).
    std::mem::swap(&mut sim.particles, &mut k2);
    let update_accelerations = rebx.allocated_forces[force_idx].update_accelerations;
    if let Some(update_accelerations) = update_accelerations {
        update_accelerations(sim, rebx, force_idx, N);
    }
    std::mem::swap(&mut sim.particles, &mut k2);

    let b1 = dt / 4.;
    let b2 = 3. * dt / 4.;
    for i in 0..N {
        let p = &mut sim.particles[i];
        p.vx += b1 * p.ax + b2 * k2[i].ax;
        p.vy += b1 * p.ay + b2 * k2[i].ay;
        p.vz += b1 * p.az + b2 * k2[i].az;
    }

    // Put the scratch buffer back where the C left it (in the C it never
    // left the parameter list; see the module docs).
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "rk2_k2", k2);
}
