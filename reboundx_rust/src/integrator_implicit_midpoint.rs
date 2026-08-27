//! integrator_implicit_midpoint.rs — translation of REBOUNDx
//! integrator_implicit_midpoint.c
//! Symplectic numerical integration scheme: it steps a single REBOUNDx
//! force with the implicit midpoint method, iterating the midpoint
//! velocities to convergence.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Parameters
//!
//! This file is one of REBOUNDx's own integrators rather than an effect,
//! so it exposes no user-facing parameters. It is selected by setting
//! the `integrator` parameter (`REBX_TYPE_INT`) of an `integrate_force`
//! operator to `REBX_INTEGRATOR_IMPLICIT_MIDPOINT`, and it keeps three
//! internal scratch buffers on the force it integrates:
//!
//! | Name          | Type    | Description                                            |
//! |---------------|---------|--------------------------------------------------------|
//! | `im_ps_final` | Pointer | Internal. Particle buffer holding the iterate's result. |
//! | `im_ps_prev`  | Pointer | Internal. Previous iterate, for the convergence test.   |
//! | `im_ps_avg`   | Pointer | Internal. Midpoint (average) state the force acts on.   |
//!
//! # Deviations from the C
//!
//! * The C hands the force its scratch buffer as the `particles`
//!   argument of `update_accelerations`. The Rust force signature has no
//!   such argument — a force always acts on `sim.particles` — so the
//!   buffer is swapped into `sim.particles` for the duration of the call
//!   and swapped straight back out. The force therefore sees exactly the
//!   array the C passes it, element for element.
//! * `force->free_memory = rebx_im_free_memory` has no counterpart:
//!   `rebx_force` owns no function pointer for it, and the buffers are
//!   `Vec`s that are dropped with the force. [`rebx_im_free_memory`] is
//!   ported for callers that want to release them early.
//! * The C `malloc`s the buffers once at `N` particles and would read
//!   and write past the end if `sim->N` later grew. Here they are grown
//!   to `N` instead.

use crate::core::{rebx_get_param_particles, rebx_set_param_particles};
use crate::types::{rebx_ap, rebx_extras, rebx_param_value};
use rebound_rs::{reb_particle, reb_simulation, reb_simulation_warning};

/// C: `static void avg_particles(...)`. Midpoint state: the average of
/// the state at the beginning of the step and the current iterate.
fn avg_particles(
    ps_avg: &mut [reb_particle],
    ps1: &[reb_particle],
    ps2: &[reb_particle],
    N: usize,
) {
    for i in 0..N {
        ps_avg[i].x = 0.5 * (ps1[i].x + ps2[i].x);
        ps_avg[i].y = 0.5 * (ps1[i].y + ps2[i].y);
        ps_avg[i].z = 0.5 * (ps1[i].z + ps2[i].z);
        ps_avg[i].vx = 0.5 * (ps1[i].vx + ps2[i].vx);
        ps_avg[i].vy = 0.5 * (ps1[i].vy + ps2[i].vy);
        ps_avg[i].vz = 0.5 * (ps1[i].vz + ps2[i].vz);
        ps_avg[i].ax = 0.;
        ps_avg[i].ay = 0.;
        ps_avg[i].az = 0.;
        ps_avg[i].m = 0.5 * (ps1[i].m + ps2[i].m);
    }
}

/// C: `static int compare(...)`. Returns 1 once the fractional change in
/// the total squared velocity between successive iterates has dropped
/// below `DBL_EPSILON*DBL_EPSILON`.
fn compare(ps1: &[reb_particle], ps2: &[reb_particle], N: usize) -> i32 {
    let mut tot2 = 0.;
    let mut deltatot2 = 0.;
    for i in 0..N {
        let dvx = ps1[i].vx - ps2[i].vx;
        let dvy = ps1[i].vy - ps2[i].vy;
        let dvz = ps1[i].vz - ps2[i].vz;
        deltatot2 += dvx * dvx + dvy * dvy + dvz * dvz;
        tot2 += ps1[i].vx * ps1[i].vx + ps1[i].vy * ps1[i].vy + ps1[i].vz * ps1[i].vz;
    }
    if deltatot2 / tot2 < f64::EPSILON * f64::EPSILON {
        1
    } else {
        0
    }
}

/// C: `static void rebx_im_free_memory(...)`, installed as the force's
/// `free_memory` callback.
///
/// The C `free`s the three malloc'd buffers. Here they are owned `Vec`s
/// that are released when the force is dropped, so this only has to drop
/// the parameters that hold them; it is public because `rebx_force` has
/// no `free_memory` slot to install it in.
pub fn rebx_im_free_memory(rebx: &mut rebx_extras, force_idx: usize) {
    if let Some(ap) = rebx.ap_mut(rebx_ap::force(force_idx)) {
        ap.retain(|param| {
            param.name != "im_ps_final" && param.name != "im_ps_prev" && param.name != "im_ps_avg"
        });
    }
}

/// C: `static struct reb_particle* setup(...)`. Allocates the three
/// scratch buffers and stores them on the force.
///
/// The C returns `ps_final`; here the caller takes the buffers back out
/// of the parameter list itself, so there is nothing to hand back. The C
/// also sets `force->free_memory` here — see the module docs.
fn setup(rebx: &mut rebx_extras, force_idx: usize, N: usize) {
    // C: malloc(N*sizeof(*ps_final)) — left uninitialized there, and
    // fully overwritten by the memcpys below before it is ever read.
    let ps_final: Vec<reb_particle> = vec![reb_particle::default(); N];
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_final", ps_final);
    let ps_prev: Vec<reb_particle> = vec![reb_particle::default(); N];
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_prev", ps_prev);
    let ps_avg: Vec<reb_particle> = vec![reb_particle::default(); N];
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_avg", ps_avg);
}

/// Move one scratch buffer out of the force's parameter list so it can
/// be written while `rebx` is handed to `update_accelerations`. The
/// parameter itself stays in the list (holding an empty buffer) exactly
/// as the C's pointer parameter stays non-NULL, so the `is_none()` check
/// in the integrator still fires only on the very first step.
fn rebx_im_take_buffer(
    rebx: &mut rebx_extras,
    force_idx: usize,
    param_name: &str,
    N: usize,
) -> Vec<reb_particle> {
    let mut buf: Vec<reb_particle> = Vec::new();
    if let Some(ap) = rebx.ap_mut(rebx_ap::force(force_idx)) {
        for param in ap.iter_mut() {
            if param.name == param_name {
                if let rebx_param_value::particles(v) = &mut param.value {
                    std::mem::swap(&mut buf, v);
                }
                break;
            }
        }
    }
    if buf.len() < N {
        buf.resize(N, reb_particle::default());
    }
    buf
}

/// C: `void rebx_integrator_implicit_midpoint_integrate(struct
/// reb_simulation* const sim, const double dt, struct rebx_force* const
/// force)`.
pub fn rebx_integrator_implicit_midpoint_integrate(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    dt: f64,
    force_idx: usize,
) {
    let N = sim.N;
    // C: ps_final = rebx_get_param(rebx, force->ap, "im_ps_final");
    //    if (ps_final == NULL){ ps_final = setup(rebx, force, N); }
    if rebx_get_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_final").is_none() {
        setup(rebx, force_idx, N);
    }
    let mut ps_final = rebx_im_take_buffer(rebx, force_idx, "im_ps_final", N);
    // These should not fail since we check above and setup if not there
    let mut ps_prev = rebx_im_take_buffer(rebx, force_idx, "im_ps_prev", N);
    let mut ps_avg = rebx_im_take_buffer(rebx, force_idx, "im_ps_avg", N);
    // C: ps_orig = sim->particles — an alias, not a copy. sim.particles
    // is not written until the very last loop, so it is read directly.
    ps_final[..N].copy_from_slice(&sim.particles[..N]);
    ps_avg[..N].copy_from_slice(&sim.particles[..N]);
    // C: `int n, converged;` — `n` outlives the loop (it decides the
    // warning below), `converged` does not.
    let mut n: i32 = 0;
    while n < 10 {
        ps_prev[..N].copy_from_slice(&ps_final[..N]);
        // C: force->update_accelerations(sim, force, ps_avg, N);
        // The force acts on sim.particles here, so ps_avg takes its
        // place for the duration of the call.
        let update_accelerations = rebx.allocated_forces[force_idx].update_accelerations;
        std::mem::swap(&mut sim.particles, &mut ps_avg);
        if let Some(update_accelerations) = update_accelerations {
            update_accelerations(sim, rebx, force_idx, N);
        }
        std::mem::swap(&mut sim.particles, &mut ps_avg);
        for i in 0..N {
            ps_final[i].vx = sim.particles[i].vx + dt * ps_avg[i].ax;
            ps_final[i].vy = sim.particles[i].vy + dt * ps_avg[i].ay;
            ps_final[i].vz = sim.particles[i].vz + dt * ps_avg[i].az;
        }
        let converged = compare(&ps_final, &ps_prev, N);
        if converged != 0 {
            break;
        }
        avg_particles(&mut ps_avg, &sim.particles, &ps_final, N);
        n += 1;
    }
    let default_max_iterations: i32 = 10;
    if n == default_max_iterations {
        reb_simulation_warning(sim, "REBOUNDx: 10 iterations in integrator_implicit_midpoint.c failed to converge. This is typically because the perturbation is too strong for the current implementation.");
    }
    for i in 0..N {
        sim.particles[i].vx = ps_final[i].vx;
        sim.particles[i].vy = ps_final[i].vy;
        sim.particles[i].vz = ps_final[i].vz;
    }

    // Hand the scratch buffers back to the force's parameter list (the C
    // never moved them; they stayed behind its pointers all along).
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_final", ps_final);
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_prev", ps_prev);
    rebx_set_param_particles(rebx, rebx_ap::force(force_idx), "im_ps_avg", ps_avg);
}
