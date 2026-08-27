//! integrator_janus.rs — the bit-wise time-reversible JANUS integrator
//! (from integrator_janus.c/h; Rein & Tamayo 2018). Positions and
//! velocities live on an int64 grid; drift and kick round to the grid
//! with C's truncating double→int64 conversion, which Rust's `as i64`
//! reproduces exactly for all in-range values.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Daniel Tamayo and contributors. See crate root.

use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::reb_simulation_error;
use crate::types::*;

/// integrator_janus.h `struct reb_particle_int`.
#[derive(Clone, Copy, Debug, Default)]
pub struct reb_particle_int {
    pub x: i64,
    pub y: i64,
    pub z: i64,
    pub vx: i64,
    pub vy: i64,
    pub vz: i64,
}

/// integrator_janus.h `struct reb_integrator_janus_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_janus_state {
    /// Scale of position grid. Default 1e-16.
    pub scale_pos: f64,
    /// Scale of velocity grid. Default 1e-16.
    pub scale_vel: f64,
    /// Order: 2, 4, 6 (default), 8, 10.
    pub order: u32,
    /// Set to 1 if particles have been modified.
    pub recalculate_integer_coordinates_this_timestep: u32,
    // Internal use
    pub p_int: Vec<reb_particle_int>,
}

impl Default for reb_integrator_janus_state {
    /// integrator_janus.c `reb_integrator_janus_create`.
    fn default() -> Self {
        reb_integrator_janus_state {
            scale_pos: 1e-16,
            scale_vel: 1e-16,
            order: 6,
            recalculate_integer_coordinates_this_timestep: 0,
            p_int: Vec::new(),
        }
    }
}

/// One specific JANUS scheme (order, stages, gamma coefficients).
struct reb_janus_scheme {
    stages: u32,
    gamma: [f64; 17],
}

const s1odr2: reb_janus_scheme = reb_janus_scheme {
    stages: 1,
    gamma: [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
};

const s5odr4: reb_janus_scheme = reb_janus_scheme {
    stages: 5,
    gamma: [
        0.41449077179437573714,
        0.41449077179437573714,
        -0.65796308717750294857,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    ],
};

const s9odr6a: reb_janus_scheme = reb_janus_scheme {
    stages: 9,
    gamma: [
        0.39216144400731413928,
        0.33259913678935943860,
        -0.70624617255763935981,
        0.082213596293550800230,
        0.79854399093482996340,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    ],
};

const s15odr8: reb_janus_scheme = reb_janus_scheme {
    stages: 15,
    gamma: [
        0.74167036435061295345,
        -0.40910082580003159400,
        0.19075471029623837995,
        -0.57386247111608226666,
        0.29906418130365592384,
        0.33462491824529818378,
        0.31529309239676659663,
        -0.79688793935291635402,
        0., 0., 0., 0., 0., 0., 0., 0., 0.,
    ],
};

const s33odr10c: reb_janus_scheme = reb_janus_scheme {
    stages: 33,
    gamma: [
        0.12313526870982994083,
        0.77644981696937310520,
        0.14905490079567045613,
        -0.17250761219393744420,
        -0.54871240818800177942,
        0.14289765421841842100,
        -0.31419193263986861997,
        0.12670943739561041022,
        0.17444734584181312998,
        0.44318544665428572929,
        -0.81948900568299084419,
        0.13382545738489583020,
        0.64509023524410605020,
        -0.71936337169922060719,
        0.20951381813463649682,
        -0.26828113140636051966,
        0.83647216092348048955,
    ],
};

/// integrator_janus.c `gg` — symmetric coefficient lookup.
fn gg(s: &reb_janus_scheme, stage: u32) -> f64 {
    if stage < (s.stages + 1) / 2 {
        s.gamma[stage as usize]
    } else {
        s.gamma[((s.stages - 1 - stage) % 17) as usize]
    }
}

/// integrator_janus.c `to_int`.
fn to_int(psi: &mut [reb_particle_int], ps: &[reb_particle], N: usize, scale_pos: f64, scale_vel: f64) {
    for i in 0..N {
        psi[i].x = (ps[i].x / scale_pos) as i64;
        psi[i].y = (ps[i].y / scale_pos) as i64;
        psi[i].z = (ps[i].z / scale_pos) as i64;
        psi[i].vx = (ps[i].vx / scale_vel) as i64;
        psi[i].vy = (ps[i].vy / scale_vel) as i64;
        psi[i].vz = (ps[i].vz / scale_vel) as i64;
    }
}

/// integrator_janus.c `to_double`.
fn to_double(ps: &mut [reb_particle], psi: &[reb_particle_int], N: usize, scale_pos: f64, scale_vel: f64) {
    for i in 0..N {
        ps[i].x = (psi[i].x as f64) * scale_pos;
        ps[i].y = (psi[i].y as f64) * scale_pos;
        ps[i].z = (psi[i].z as f64) * scale_pos;
        ps[i].vx = (psi[i].vx as f64) * scale_vel;
        ps[i].vy = (psi[i].vy as f64) * scale_vel;
        ps[i].vz = (psi[i].vz as f64) * scale_vel;
    }
}

/// integrator_janus.c `drift`. C's signed-integer addition; wrapping_add
/// matches the two's-complement machine behavior.
fn drift(janus: &mut reb_integrator_janus_state, N: usize, dt: f64, scale_pos: f64, scale_vel: f64) {
    for i in 0..N {
        let p = &mut janus.p_int[i];
        p.x = p.x.wrapping_add((dt * (p.vx as f64) * scale_vel / scale_pos) as i64);
        p.y = p.y.wrapping_add((dt * (p.vy as f64) * scale_vel / scale_pos) as i64);
        p.z = p.z.wrapping_add((dt * (p.vz as f64) * scale_vel / scale_pos) as i64);
    }
}

/// integrator_janus.c `kick`.
fn kick(r: &reb_simulation, janus: &mut reb_integrator_janus_state, dt: f64, scale_vel: f64) {
    let N = r.N;
    for i in 0..N {
        let p = &mut janus.p_int[i];
        p.vx = p.vx.wrapping_add((dt * r.particles[i].ax / scale_vel) as i64);
        p.vy = p.vy.wrapping_add((dt * r.particles[i].ay / scale_vel) as i64);
        p.vz = p.vz.wrapping_add((dt * r.particles[i].az / scale_vel) as i64);
    }
}

/// integrator_janus.c `reb_integrator_janus_step` (state-explicit).
pub fn reb_integrator_janus_step_state(
    r: &mut reb_simulation,
    janus: &mut reb_integrator_janus_state,
) {
    r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    let N = r.N;
    let dt = r.dt;
    let scale_vel = janus.scale_vel;
    let scale_pos = janus.scale_pos;
    if janus.p_int.len() != N {
        janus.p_int.resize(N, reb_particle_int::default());
        janus.recalculate_integer_coordinates_this_timestep = 1;
    }

    if janus.recalculate_integer_coordinates_this_timestep == 1 {
        to_int(&mut janus.p_int, &r.particles, N, scale_pos, scale_vel);
        janus.recalculate_integer_coordinates_this_timestep = 0;
    }

    let s: &reb_janus_scheme = match janus.order {
        2 => &s1odr2,
        4 => &s5odr4,
        6 => &s9odr6a,
        8 => &s15odr8,
        10 => &s33odr10c,
        _ => {
            reb_simulation_error(r, "Order not supported in JANUS.");
            &s1odr2
        }
    };

    drift(janus, N, gg(s, 0) * dt / 2., scale_pos, scale_vel);
    to_double(&mut r.particles, &janus.p_int, N, scale_pos, scale_vel);

    reb_simulation_update_acceleration(r);

    kick(r, janus, gg(s, 0) * dt, scale_vel);
    for i in 1..s.stages {
        drift(janus, N, (gg(s, i - 1) + gg(s, i)) * dt / 2., scale_pos, scale_vel);
        to_double(&mut r.particles, &janus.p_int, N, scale_pos, scale_vel);
        reb_simulation_update_acceleration(r);
        kick(r, janus, gg(s, i) * dt, scale_vel);
    }
    drift(janus, N, gg(s, s.stages - 1) * dt / 2., scale_pos, scale_vel);

    // Always get positions and velocities in floating point at the end.
    reb_integrator_janus_synchronize_state(r, janus);

    r.t += r.dt;
}

/// integrator_janus.c `reb_integrator_janus_synchronize`.
pub fn reb_integrator_janus_synchronize_state(
    r: &mut reb_simulation,
    janus: &mut reb_integrator_janus_state,
) {
    if janus.p_int.len() == r.N {
        to_double(&mut r.particles, &janus.p_int, r.N, janus.scale_pos, janus.scale_vel);
    }
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_janus_step(r: &mut reb_simulation) {
    let mut janus = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::janus(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_janus_step_state(r, &mut janus);
    r.integrator = reb_integrator_state::janus(janus);
}

/// Synchronize entry point for the dispatcher.
pub fn reb_integrator_janus_synchronize(r: &mut reb_simulation) {
    let mut janus = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::janus(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_janus_synchronize_state(r, &mut janus);
    r.integrator = reb_integrator_state::janus(janus);
}
