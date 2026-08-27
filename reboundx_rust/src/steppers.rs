//! steppers.rs — translation of REBOUNDx steppers.c
//! Wrapper operators that take a step with one of REBOUND's integrators,
//! so custom operator-splitting schemes can be assembled from the pieces.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Misc
//!
//! ======================= ===============================================
//! Authors                 D. Tamayo, H. Rein
//! Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
//! Based on                `Rein and Liu, 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.128R/abstract>`_.
//! C Example               None
//! Python Example          `CustomSplittingIntegrationSchemes.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CustomSplittingIntegrationSchemes.ipynb>`_.
//! ======================= ===============================================
//!
//! These are wrapper functions to taking steps with several of REBOUND's
//! integrators in order to build custom splitting schemes.
//!
//! **Effect Parameters**
//!
//! None
//!
//! **Particle Parameters**
//!
//! None

use rebound_rs::integrator_ias15::{reb_integrator_ias15_state, reb_integrator_ias15_step_state};
use rebound_rs::integrator_whfast::{
    reb_integrator_whfast_com_step, reb_integrator_whfast_from_inertial,
    reb_integrator_whfast_init, reb_integrator_whfast_interaction_step,
    reb_integrator_whfast_jump_step, reb_integrator_whfast_kepler_step,
    reb_integrator_whfast_state, reb_integrator_whfast_to_inertial,
};
use rebound_rs::{reb_simulation, reb_simulation_update_acceleration, REB_GRAVITY_IGNORE_TERMS_NONE};

use crate::types::rebx_extras;

// will do IAS with gravity + any additional_forces

/// steppers.c `rebx_ias15_step`.
pub fn rebx_ias15_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    let old_t = sim.t;
    let t_needed = old_t + dt;
    let old_dt = sim.dt;
    sim.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;

    // C: reb_integrator_ias15.create()
    let mut ias15 = reb_integrator_ias15_state::default();

    sim.dt = 0.0001 * dt; // start with a small timestep.

    while sim.t < t_needed && (sim.dt / old_dt).abs() > 1e-14 {
        reb_integrator_ias15_step_state(sim, &mut ias15);
        if sim.t + sim.dt > t_needed {
            sim.dt = t_needed - sim.t;
        }
    }
    // C: reb_integrator_ias15.free(ias15) — the owned state drops here.
    sim.t = old_t;
    sim.dt = old_dt; // reset in case this is part of a chain of steps
}

/// steppers.c `rebx_kepler_step`.
pub fn rebx_kepler_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    // C: reb_integrator_whfast.create()
    let mut whfast = reb_integrator_whfast_state::default();
    reb_integrator_whfast_init(sim, &mut whfast);
    let coordinates = whfast.coordinates;
    reb_integrator_whfast_from_inertial(
        sim,
        &mut whfast.p_jh,
        &mut whfast.p_jh_var,
        coordinates,
    );
    reb_integrator_whfast_kepler_step(
        sim,
        &mut whfast.p_jh,
        &mut whfast.p_jh_var,
        coordinates,
        dt,
    );
    reb_integrator_whfast_com_step(sim, &mut whfast.p_jh, &mut whfast.p_jh_var, dt);
    reb_integrator_whfast_to_inertial(sim, &whfast.p_jh, &whfast.p_jh_var, coordinates);
    // C: reb_integrator_whfast.free(whfast) — the owned state drops here.
}

/// steppers.c `rebx_jump_step`.
pub fn rebx_jump_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    // TODO: This will never do anything because by default whfast coordinates are Jacobi which do not have a jump step.
    let mut whfast = reb_integrator_whfast_state::default();
    reb_integrator_whfast_init(sim, &mut whfast);
    let coordinates = whfast.coordinates;
    reb_integrator_whfast_from_inertial(
        sim,
        &mut whfast.p_jh,
        &mut whfast.p_jh_var,
        coordinates,
    );
    reb_integrator_whfast_jump_step(sim, &mut whfast.p_jh, coordinates, dt);
    reb_integrator_whfast_to_inertial(sim, &whfast.p_jh, &whfast.p_jh_var, coordinates);
    // C: reb_integrator_whfast.free(whfast) — the owned state drops here.
}

/// steppers.c `rebx_interaction_step`.
pub fn rebx_interaction_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    let mut whfast = reb_integrator_whfast_state::default();
    reb_integrator_whfast_init(sim, &mut whfast);
    let coordinates = whfast.coordinates;
    reb_integrator_whfast_from_inertial(
        sim,
        &mut whfast.p_jh,
        &mut whfast.p_jh_var,
        coordinates,
    );
    reb_simulation_update_acceleration(sim);
    reb_integrator_whfast_interaction_step(
        sim,
        &mut whfast.p_jh,
        &mut whfast.p_jh_var,
        coordinates,
        dt,
    );
    reb_integrator_whfast_to_inertial(sim, &whfast.p_jh, &whfast.p_jh_var, coordinates);
    // C: reb_integrator_whfast.free(whfast) — the owned state drops here.
}

/// steppers.c `rebx_drift_step`.
pub fn rebx_drift_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    let N = sim.N;
    let particles = &mut sim.particles;
    for i in 0..N {
        particles[i].x += dt * particles[i].vx;
        particles[i].y += dt * particles[i].vy;
        particles[i].z += dt * particles[i].vz;
    }
}

/// steppers.c `rebx_kick_step`.
pub fn rebx_kick_step(
    sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    reb_simulation_update_acceleration(sim);
    let N = sim.N;
    let particles = &mut sim.particles;
    for i in 0..N {
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
    }
}
