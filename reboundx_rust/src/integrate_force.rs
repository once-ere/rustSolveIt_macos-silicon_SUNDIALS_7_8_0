//! integrate_force.rs — translation of REBOUNDx integrate_force.c
//! Generic operator for integrating a force across a timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! `@file integrate_force.c`, `@brief Generic operator for integrating a
//! force across a timestep`, `@author Dan Tamayo
//! <tamayo.daniel@gmail.com>, Hanno Rein`.
//!
//! # Parameters
//!
//! The C file carries no parameter table of its own; these are the two
//! parameters the operator reads off its own `ap` list (both are
//! registered by `rebx_register_default_params` in core.c).
//!
//! | Name         | Type              | Required | Description                                                                                                            |
//! |--------------|-------------------|----------|------------------------------------------------------------------------------------------------------------------------|
//! | `force`      | `REBX_TYPE_FORCE` | Yes      | The force this operator integrates across the timestep.                                                                  |
//! | `integrator` | `REBX_TYPE_INT`   | No       | Which REBOUNDx integrator to use (a `rebx_integrator` value). Defaults to `REBX_INTEGRATOR_EULER` when it is not set.     |
//!
//! # Deviations from the C
//!
//! * When the `force` parameter is not set, the C reports the error with
//!   `reb_simulation_error` and then falls straight through into the
//!   switch, handing the chosen integrator a `NULL` force pointer, which
//!   it dereferences. There is no safe counterpart to that, so this port
//!   reports the identical error message and then performs no
//!   integration. Everything the C does before the dereference — the
//!   `integrator` lookup and `rebx_reset_accelerations` — still happens,
//!   in the C's order.
//! * The REBOUNDx integrators take `(sim, rebx, dt, force_idx)` here:
//!   `rebx` is the state the C reaches through `sim->extras`, and
//!   `force_idx` indexes `rebx.allocated_forces` in place of the C's
//!   `struct rebx_force*`.

use rebound_rs::{reb_simulation, reb_simulation_error};

use crate::core::{rebx_get_param_force, rebx_get_param_int, rebx_reset_accelerations};
use crate::integrator_euler::rebx_integrator_euler_integrate;
use crate::integrator_implicit_midpoint::rebx_integrator_implicit_midpoint_integrate;
use crate::integrator_rk2::rebx_integrator_rk2_integrate;
use crate::integrator_rk4::rebx_integrator_rk4_integrate;
use crate::types::{
    rebx_ap, rebx_extras, rebx_integrator, REBX_INTEGRATOR_EULER,
    REBX_INTEGRATOR_IMPLICIT_MIDPOINT, REBX_INTEGRATOR_NONE, REBX_INTEGRATOR_RK2,
    REBX_INTEGRATOR_RK4,
};

/// integrate_force.c `rebx_integrate_force`.
///
/// C: `void rebx_integrate_force(struct reb_simulation* const sim,
/// struct rebx_operator* const operator, const double dt)`.
///
/// The C's `struct rebx_extras* rebx = sim->extras;` is the `rebx`
/// argument, and `operator` is `operator_idx`, an index into
/// `rebx.allocated_operators` (so C's `operator->ap` is
/// `rebx_ap::operator_(operator_idx)`).
pub fn rebx_integrate_force(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    operator_idx: usize,
    dt: f64,
) {
    // C: struct rebx_force* force = rebx_get_param(rebx, operator->ap, "force");
    let force = rebx_get_param_force(rebx, rebx_ap::operator_(operator_idx), "force");
    if force.is_none() {
        reb_simulation_error(sim, "REBOUNDx Error: Force parameter not set in rebx_integrate operator. See examples for how to add as a parameter.\n");
    }
    let mut integrator: rebx_integrator = REBX_INTEGRATOR_EULER; // default
    let integratorparam = rebx_get_param_int(rebx, rebx_ap::operator_(operator_idx), "integrator");
    if let Some(integratorparam) = integratorparam {
        integrator = integratorparam;
    }

    // C: rebx_reset_accelerations(sim->particles, sim->N);
    let N = sim.N;
    rebx_reset_accelerations(&mut sim.particles, N);

    // See the module docs: the C would dereference the NULL force here.
    let force = match force {
        Some(force) => force,
        None => return,
    };

    match integrator {
        REBX_INTEGRATOR_IMPLICIT_MIDPOINT => {
            rebx_integrator_implicit_midpoint_integrate(sim, rebx, dt, force);
        }
        REBX_INTEGRATOR_RK2 => {
            rebx_integrator_rk2_integrate(sim, rebx, dt, force);
        }
        REBX_INTEGRATOR_RK4 => {
            rebx_integrator_rk4_integrate(sim, rebx, dt, force);
        }
        REBX_INTEGRATOR_EULER => {
            rebx_integrator_euler_integrate(sim, rebx, dt, force);
        }
        REBX_INTEGRATOR_NONE => {}
        _ => {}
    }
}
