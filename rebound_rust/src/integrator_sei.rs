//! integrator_sei.rs — the Symplectic Epicycle Integrator
//! (from integrator_sei.c; Rein & Tremaine 2011).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::simulation::reb_simulation_update_acceleration;
use crate::types::*;

/// integrator_sei.c `struct reb_integrator_sei_state` — cached
/// trigonometric values for the current dt.
#[derive(Clone, Copy, Debug, Default)]
pub struct reb_integrator_sei_state {
    pub lastdt: f64,
    pub sindt: f64,
    pub tandt: f64,
    pub sindtz: f64,
    pub tandtz: f64,
}

/// integrator_sei.c `operator_H012` — the exact epicyclic propagator,
/// its rotation implemented as three shears to avoid round-off.
fn operator_H012(dt: f64, ri_sei: &reb_integrator_sei_state, p: &mut reb_particle, OMEGA: f64, OMEGAZ: f64) {
    // Integrate vertical motion
    let zx = p.z * OMEGAZ;
    let zy = p.vz;

    let zt1 = zx - ri_sei.tandtz * zy;
    let zyt = ri_sei.sindtz * zt1 + zy;
    let zxt = zt1 - ri_sei.tandtz * zyt;
    p.z = zxt / OMEGAZ;
    p.vz = zyt;

    // Integrate motion in xy directions
    let aO = 2. * p.vy + 4. * p.x * OMEGA; // Center of epicyclic motion
    let bO = p.y * OMEGA - 2. * p.vx;

    let ys = (p.y * OMEGA - bO) / 2.; // Epicycle vector
    let xs = p.x * OMEGA - aO;

    let xst1 = xs - ri_sei.tandt * ys;
    let yst = ri_sei.sindt * xst1 + ys;
    let xst = xst1 - ri_sei.tandt * yst;

    p.x = (xst + aO) / OMEGA;
    p.y = (yst * 2. + bO) / OMEGA - 3. / 4. * aO * dt;
    p.vx = yst;
    p.vy = -xst * 2. - 3. / 2. * aO;
}

/// integrator_sei.c `operator_phi1` — kick.
fn operator_phi1(dt: f64, p: &mut reb_particle) {
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.vz += p.az * dt;
}

/// integrator_sei.c `reb_integrator_sei_step`.
pub fn reb_integrator_sei_step(r: &mut reb_simulation) {
    r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    let N = r.N;
    let mut sei = match &r.integrator {
        reb_integrator_state::sei(s) => *s,
        _ => reb_integrator_sei_state::default(),
    };
    if sei.lastdt != r.dt {
        // Pre-calculates sin() and tan() needed for SEI.
        if r.OMEGAZ == -1. {
            r.OMEGAZ = r.OMEGA;
        }
        sei.sindt = (r.OMEGA * (-r.dt / 2.)).sin();
        sei.tandt = (r.OMEGA * (-r.dt / 4.)).tan();
        sei.sindtz = (r.OMEGAZ * (-r.dt / 2.)).sin();
        sei.tandtz = (r.OMEGAZ * (-r.dt / 4.)).tan();
        sei.lastdt = r.dt;
    }
    for i in 0..N {
        operator_H012(r.dt, &sei, &mut r.particles[i], r.OMEGA, r.OMEGAZ);
    }
    r.t += r.dt / 2.;

    reb_simulation_update_acceleration(r);

    for i in 0..N {
        operator_phi1(r.dt, &mut r.particles[i]);
        operator_H012(r.dt, &sei, &mut r.particles[i], r.OMEGA, r.OMEGAZ);
    }
    r.t += r.dt / 2.;
    r.dt_last_done = r.dt;
    r.integrator = reb_integrator_state::sei(sei);
}
