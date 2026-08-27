//! integrator_leapfrog.rs — the standard leapfrog integrator and its
//! 4th/6th/8th-order generalizations (from integrator_leapfrog.c;
//! Yoshida 4th order, Blanes & Casas 2016 p91 for 6/8).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::reb_simulation_error;
use crate::types::*;

/// integrator_leapfrog.c `struct reb_integrator_leapfrog_state`.
#[derive(Clone, Copy, Debug)]
pub struct reb_integrator_leapfrog_state {
    /// Order of the integrator. Default is 2. Other allowed values are
    /// 4, 6 and 8.
    pub order: u32,
}

impl Default for reb_integrator_leapfrog_state {
    fn default() -> Self {
        reb_integrator_leapfrog_state { order: 2 }
    }
}

pub const reb_integrator_leapfrog_lf4_a: f64 = 0.675603595979828817023843904485;
pub const reb_integrator_leapfrog_lf6_a: [f64; 5] = [
    0.1867,
    0.5554970237124784,
    0.1294669489134754,
    -0.843265623387734,
    0.9432033015235604,
];
pub const reb_integrator_leapfrog_lf8_a: [f64; 9] = [
    0.128865979381443,
    0.581514087105251,
    -0.410175371469850,
    0.1851469357165877,
    -0.4095523434208514,
    0.1444059410800120,
    0.2783355003936797,
    0.3149566839162949,
    -0.6269948254051343979,
];

fn drift(r: &mut reb_simulation, dt: f64) {
    let N = r.N;
    for i in 0..N {
        r.particles[i].x += dt * r.particles[i].vx;
        r.particles[i].y += dt * r.particles[i].vy;
        r.particles[i].z += dt * r.particles[i].vz;
    }
    r.t += dt; // advance time so that force evaluations are correct
}

fn kick(r: &mut reb_simulation, dt: f64) {
    let N = r.N;
    for i in 0..N {
        r.particles[i].vx += dt * r.particles[i].ax;
        r.particles[i].vy += dt * r.particles[i].ay;
        r.particles[i].vz += dt * r.particles[i].az;
    }
}

/// integrator_leapfrog.c `reb_integrator_leapfrog_step`
/// (Drift-Kick-Drift, non-rotating frame).
pub fn reb_integrator_leapfrog_step(r: &mut reb_simulation) {
    r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    let dt = r.dt;
    let order = match &r.integrator {
        reb_integrator_state::leapfrog(s) => s.order,
        _ => 2,
    };
    match order {
        2 => {
            drift(r, dt * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt);
            drift(r, dt * 0.5);
        }
        4 => {
            let a = reb_integrator_leapfrog_lf4_a;
            drift(r, dt * a);
            reb_simulation_update_acceleration(r);
            kick(r, dt * 2. * a);
            drift(r, dt * (0.5 - a));
            reb_simulation_update_acceleration(r);
            kick(r, dt * (1. - 4. * a));
            drift(r, dt * (0.5 - a));
            reb_simulation_update_acceleration(r);
            kick(r, dt * 2. * a);
            drift(r, dt * a);
        }
        6 => {
            let a = &reb_integrator_leapfrog_lf6_a;
            drift(r, dt * a[0] * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[0]);
            drift(r, dt * (a[0] + a[1]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[1]);
            drift(r, dt * (a[1] + a[2]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[2]);
            drift(r, dt * (a[2] + a[3]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[3]);
            drift(r, dt * (a[3] + a[4]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[4]);
            drift(r, dt * (a[3] + a[4]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[3]);
            drift(r, dt * (a[2] + a[3]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[2]);
            drift(r, dt * (a[1] + a[2]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[1]);
            drift(r, dt * (a[0] + a[1]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[0]);
            drift(r, dt * a[0] * 0.5);
        }
        8 => {
            let a = &reb_integrator_leapfrog_lf8_a;
            drift(r, dt * a[0] * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[0]);
            drift(r, dt * (a[0] + a[1]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[1]);
            drift(r, dt * (a[1] + a[2]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[2]);
            drift(r, dt * (a[2] + a[3]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[3]);
            drift(r, dt * (a[3] + a[4]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[4]);
            drift(r, dt * (a[4] + a[5]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[5]);
            drift(r, dt * (a[5] + a[6]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[6]);
            drift(r, dt * (a[6] + a[7]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[7]);
            drift(r, dt * (a[7] + a[8]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[8]);
            drift(r, dt * (a[7] + a[8]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[7]);
            drift(r, dt * (a[6] + a[7]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[6]);
            drift(r, dt * (a[5] + a[6]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[5]);
            drift(r, dt * (a[4] + a[5]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[4]);
            drift(r, dt * (a[3] + a[4]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[3]);
            drift(r, dt * (a[2] + a[3]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[2]);
            drift(r, dt * (a[1] + a[2]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[1]);
            drift(r, dt * (a[0] + a[1]) * 0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt * a[0]);
            drift(r, dt * a[0] * 0.5);
        }
        _ => {
            reb_simulation_error(r, "Leapfrog order not supported.");
            return;
        }
    }
    r.dt_last_done = dt;
}
