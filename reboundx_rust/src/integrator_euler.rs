//! integrator_euler.rs — translation of REBOUNDx integrator_euler.c
//! Euler's method: one force evaluation, one first-order velocity update.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! The C file (`euler.c`, "Euler's method", Dan Tamayo & Hanno Rein) carries
//! no parameter table: this integrator takes no parameters of its own. It is
//! selected by setting the `integrator` parameter of the `integrate_force`
//! operator to `REBX_INTEGRATOR_EULER`, which is also the default.

use crate::types::rebx_extras;
use rebound_rs::reb_simulation;

/// core.h `rebx_integrator_euler_integrate`.
///
/// C signature:
/// `void rebx_integrator_euler_integrate(struct reb_simulation* const sim,
/// const double dt, struct rebx_force* const force)`.
///
/// `force_idx` names the force by its index into
/// `rebx_extras::allocated_forces` (the C passes the pointer). The C calls
/// `force->update_accelerations(sim, force, sim->particles, N)`; here the
/// particles are `sim.particles`, mutated in place, so the callee sees the
/// same array the C does.
pub fn rebx_integrator_euler_integrate(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    dt: f64,
    force_idx: usize,
) {
    let N = sim.N;
    // Read the function pointer out in its own statement: it must not still
    // be borrowed from `rebx` when `rebx` is handed to it.
    let update_accelerations = rebx.allocated_forces[force_idx].update_accelerations;
    if let Some(update_accelerations) = update_accelerations {
        update_accelerations(sim, rebx, force_idx, N);
    }
    for i in 0..N {
        sim.particles[i].vx += dt * sim.particles[i].ax;
        sim.particles[i].vy += dt * sim.particles[i].ay;
        sim.particles[i].vz += dt * sim.particles[i].az;
    }
}
