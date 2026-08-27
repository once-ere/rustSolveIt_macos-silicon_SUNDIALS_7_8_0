//! integrator_eos.rs — the Embedded Operator Splitting (EOS) methods
//! (from integrator_eos.c/h; Rein 2019). Two nested splitting schemes
//! phi0 (outer) and phi1 (inner, n sub-steps), each selectable from
//! LF, LF4, LF6, LF8, LF4_2, LF8_6_4, PLF7_6_4, PMLF4, PMLF6.
//!
//! The C passes drift/interaction function pointers into the pre- and
//! postprocessors; here the shell (0 or 1) selects them, with identical
//! call sequences.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein and contributors. See crate root.

use crate::gravity::reb_gravity_basic_calculate_and_apply_jerk;
use crate::integrator_leapfrog::{
    reb_integrator_leapfrog_lf4_a, reb_integrator_leapfrog_lf6_a, reb_integrator_leapfrog_lf8_a,
};
use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::reb_simulation_warning;
use crate::types::*;

pub const REB_INTEGRATOR_EOS_TYPE_LF: i32 = 0;
pub const REB_INTEGRATOR_EOS_TYPE_LF4: i32 = 1;
pub const REB_INTEGRATOR_EOS_TYPE_LF6: i32 = 2;
pub const REB_INTEGRATOR_EOS_TYPE_LF8: i32 = 3;
pub const REB_INTEGRATOR_EOS_TYPE_LF4_2: i32 = 4;
pub const REB_INTEGRATOR_EOS_TYPE_LF8_6_4: i32 = 5;
pub const REB_INTEGRATOR_EOS_TYPE_PLF7_6_4: i32 = 6;
pub const REB_INTEGRATOR_EOS_TYPE_PMLF4: i32 = 7;
pub const REB_INTEGRATOR_EOS_TYPE_PMLF6: i32 = 8;

/// integrator_eos.h `struct reb_integrator_eos_state`.
#[derive(Clone, Copy, Debug)]
pub struct reb_integrator_eos_state {
    /// Outer operator splitting method.
    pub phi0: i32,
    /// Inner operator splitting method.
    pub phi1: i32,
    /// Number of inner splittings per outer splitting. Default: 2.
    pub n: u32,
    /// Combine kick steps at beginning and end of timestep.
    pub safe_mode: u32,
}

impl Default for reb_integrator_eos_state {
    /// integrator_eos.c `reb_integrator_eos_create`.
    fn default() -> Self {
        reb_integrator_eos_state {
            phi0: REB_INTEGRATOR_EOS_TYPE_LF,
            phi1: REB_INTEGRATOR_EOS_TYPE_LF,
            n: 2,
            safe_mode: 1,
        }
    }
}

const lf4_2_a: f64 = 0.211324865405187117745425609749;

const lf8_6_4_a: [f64; 4] = [
    0.0711334264982231177779387300061549964174,
    0.241153427956640098736487795326289649618,
    0.521411761772814789212136078067994229991,
    -0.333698616227678005726562603400438876027,
];
const lf8_6_4_b: [f64; 4] = [
    0.183083687472197221961703757166430291072,
    0.310782859898574869507522291054262796375,
    -0.0265646185119588006972121379164987592663,
    0.0653961422823734184559721793911134363710,
];

const pmlf6_a: [f64; 2] = [-0.0682610383918630, 0.568261038391863038121699];
const pmlf6_b: [f64; 2] = [0.2621129352517028, 0.475774129496594366806050];
const pmlf6_c: [f64; 2] = [0., 0.0164011128160783];
const pmlf6_z: [f64; 6] = [
    0.07943288242455420,
    0.02974829169467665,
    -0.7057074964815896,
    0.3190423451260838,
    -0.2869147334299646,
    0.564398710666239478150885,
];
const pmlf6_y: [f64; 6] = [
    1.3599424487455264,
    -0.6505973747535132,
    -0.033542814598338416,
    -0.040129915275115030,
    0.044579729809902803,
    -0.680252073928462652752103,
];
const pmlf6_v: [f64; 6] = [
    -0.034841228074994859,
    0.031675672097525204,
    -0.005661054677711889,
    0.004262222269023640,
    0.005,
    -0.005,
];

const pmlf4_y: [f64; 3] = [0.1859353996846055, 0.0731969797858114, -0.1576624269298081];
const pmlf4_z: [f64; 3] = [0.8749306155955435, -0.237106680151022, -0.5363539829039128];

const plf7_6_4_a: [f64; 2] = [0.5600879810924619, -0.060087981092461900000];
const plf7_6_4_b: [f64; 2] = [1.5171479707207228, -2.0342959414414456000];
const plf7_6_4_z: [f64; 6] = [
    -0.3346222298730800,
    1.0975679907321640,
    -1.0380887460967830,
    0.6234776317921379,
    -1.1027532063031910,
    -0.0141183222088869,
];
const plf7_6_4_y: [f64; 6] = [
    -1.6218101180868010,
    0.0061709468110142,
    0.8348493592472594,
    -0.0511253369989315,
    0.5633782670698199,
    -0.5,
];

/// integrator_eos.c `reb_integrator_eos_interaction_shell0`.
fn interaction_shell0(r: &mut reb_simulation, _eos: &reb_integrator_eos_state, y: f64, v: f64) {
    // Calculate gravity using standard gravity routine
    r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_INVOLVING_0;
    r.gravity = REB_GRAVITY::BASIC;
    reb_simulation_update_acceleration(r);
    if v != 0. {
        reb_gravity_basic_calculate_and_apply_jerk(r, v);
    }
    // Apply acceleration (jerk already applied)
    let N = r.N;
    for i in 0..N {
        r.particles[i].vx += y * r.particles[i].ax;
        r.particles[i].vy += y * r.particles[i].ay;
        r.particles[i].vz += y * r.particles[i].az;
    }
}

/// integrator_eos.c `reb_integrator_eos_interaction_shell1`
/// (Kepler-shell interactions with the central object only).
fn interaction_shell1(r: &mut reb_simulation, _eos: &reb_integrator_eos_state, y: f64, v: f64) {
    let N = r.N;
    let N_active = if r.N_active == usize::MAX { N } else { r.N_active };
    let testparticle_type = r.testparticle_type;
    let G = r.G;

    if v != 0. {
        // Normal force calculation
        r.particles[0].ax = 0.;
        r.particles[0].ay = 0.;
        r.particles[0].az = 0.;
        // Interactions between central object and all other active particles
        for j in 1..N_active {
            let dx = r.particles[0].x - r.particles[j].x;
            let dy = r.particles[0].y - r.particles[j].y;
            let dz = r.particles[0].z - r.particles[j].z;
            let dr = (dx * dx + dy * dy + dz * dz).sqrt();

            let prefact = G / (dr * dr * dr);
            let prefactj = -prefact * r.particles[j].m;
            r.particles[0].ax += prefactj * dx;
            r.particles[0].ay += prefactj * dy;
            r.particles[0].az += prefactj * dz;
            let prefacti = prefact * r.particles[0].m;
            r.particles[j].ax = prefacti * dx;
            r.particles[j].ay = prefacti * dy;
            r.particles[j].az = prefacti * dz;
        }
        // Interactions between central object and all test particles
        for j in N_active..N {
            let dx = r.particles[0].x - r.particles[j].x;
            let dy = r.particles[0].y - r.particles[j].y;
            let dz = r.particles[0].z - r.particles[j].z;
            let dr = (dx * dx + dy * dy + dz * dz).sqrt();

            let prefact = G / (dr * dr * dr);
            let prefacti = prefact * r.particles[0].m;
            r.particles[j].ax = prefacti * dx;
            r.particles[j].ay = prefacti * dy;
            r.particles[j].az = prefacti * dz;
            if testparticle_type != 0 {
                let prefactj = -prefact * r.particles[j].m;
                r.particles[0].ax += prefactj * dx;
                r.particles[0].ay += prefactj * dy;
                r.particles[0].az += prefactj * dz;
            }
        }
        // Jerk calculation
        // Interactions between central object and all other active particles
        for i in 1..N_active {
            let dx = r.particles[0].x - r.particles[i].x;
            let dy = r.particles[0].y - r.particles[i].y;
            let dz = r.particles[0].z - r.particles[i].z;

            let dax = r.particles[0].ax - r.particles[i].ax;
            let day = r.particles[0].ay - r.particles[i].ay;
            let daz = r.particles[0].az - r.particles[i].az;

            let dr = (dx * dx + dy * dy + dz * dz).sqrt();
            let alphasum = dax * dx + day * dy + daz * dz;
            let prefact2 = 2. * v * G / (dr * dr * dr);
            let prefact2i = prefact2 * r.particles[i].m;
            let prefact2j = prefact2 * r.particles[0].m;
            let prefact1 = alphasum * prefact2 / dr * 3. / dr;
            let prefact1i = prefact1 * r.particles[i].m;
            let prefact1j = prefact1 * r.particles[0].m;
            r.particles[0].vx += -dax * prefact2i + dx * prefact1i;
            r.particles[0].vy += -day * prefact2i + dy * prefact1i;
            r.particles[0].vz += -daz * prefact2i + dz * prefact1i;
            r.particles[i].vx += y * r.particles[i].ax + dax * prefact2j - dx * prefact1j;
            r.particles[i].vy += y * r.particles[i].ay + day * prefact2j - dy * prefact1j;
            r.particles[i].vz += y * r.particles[i].az + daz * prefact2j - dz * prefact1j;
        }
        // Interactions between central object and all test particles
        for i in N_active..N {
            let dx = r.particles[0].x - r.particles[i].x;
            let dy = r.particles[0].y - r.particles[i].y;
            let dz = r.particles[0].z - r.particles[i].z;

            let dax = r.particles[0].ax - r.particles[i].ax;
            let day = r.particles[0].ay - r.particles[i].ay;
            let daz = r.particles[0].az - r.particles[i].az;

            let dr = (dx * dx + dy * dy + dz * dz).sqrt();
            let alphasum = dax * dx + day * dy + daz * dz;
            let prefact2 = 2. * v * G / (dr * dr * dr);
            let prefact2j = prefact2 * r.particles[0].m;
            let prefact1 = alphasum * prefact2 / dr * 3. / dr;
            let prefact1j = prefact1 * r.particles[0].m;
            if testparticle_type != 0 {
                let prefact2i = prefact2 * r.particles[i].m;
                let prefact1i = prefact1 * r.particles[i].m;
                r.particles[0].vx += -dax * prefact2i + dx * prefact1i;
                r.particles[0].vy += -day * prefact2i + dy * prefact1i;
                r.particles[0].vz += -daz * prefact2i + dz * prefact1i;
            }
            r.particles[i].vx += y * r.particles[i].ax + dax * prefact2j - dx * prefact1j;
            r.particles[i].vy += y * r.particles[i].ay + day * prefact2j - dy * prefact1j;
            r.particles[i].vz += y * r.particles[i].az + daz * prefact2j - dz * prefact1j;
        }
        r.particles[0].vx += y * r.particles[0].ax;
        r.particles[0].vy += y * r.particles[0].ay;
        r.particles[0].vz += y * r.particles[0].az;
    } else {
        // Normal force calculation
        for j in 1..N_active {
            let dx = r.particles[0].x - r.particles[j].x;
            let dy = r.particles[0].y - r.particles[j].y;
            let dz = r.particles[0].z - r.particles[j].z;
            let dr = (dx * dx + dy * dy + dz * dz).sqrt();

            let prefact = y * G / (dr * dr * dr);
            let prefactj = -prefact * r.particles[j].m;
            r.particles[0].vx += prefactj * dx;
            r.particles[0].vy += prefactj * dy;
            r.particles[0].vz += prefactj * dz;
            let prefacti = prefact * r.particles[0].m;
            r.particles[j].vx += prefacti * dx;
            r.particles[j].vy += prefacti * dy;
            r.particles[j].vz += prefacti * dz;
        }
        for j in N_active..N {
            let dx = r.particles[0].x - r.particles[j].x;
            let dy = r.particles[0].y - r.particles[j].y;
            let dz = r.particles[0].z - r.particles[j].z;
            let dr = (dx * dx + dy * dy + dz * dz).sqrt();

            let prefact = y * G / (dr * dr * dr);
            let prefacti = prefact * r.particles[0].m;
            r.particles[j].vx += prefacti * dx;
            r.particles[j].vy += prefacti * dy;
            r.particles[j].vz += prefacti * dz;
            if testparticle_type != 0 {
                let prefactj = -prefact * r.particles[j].m;
                r.particles[0].vx += prefactj * dx;
                r.particles[0].vy += prefactj * dy;
                r.particles[0].vz += prefactj * dz;
            }
        }
    }
}

/// integrator_eos.c `reb_integrator_eos_drift_shell1`.
fn drift_shell1(r: &mut reb_simulation, _eos: &reb_integrator_eos_state, dt: f64) {
    let N = r.N;
    for i in 0..N {
        r.particles[i].x += dt * r.particles[i].vx;
        r.particles[i].y += dt * r.particles[i].vy;
        r.particles[i].z += dt * r.particles[i].vz;
    }
}

/// The shell selects which drift/interaction pair the pre- and
/// postprocessors call (C: function pointers).
#[derive(Clone, Copy, PartialEq)]
enum Shell {
    Zero,
    One,
}

fn drift_step(r: &mut reb_simulation, eos: &reb_integrator_eos_state, shell: Shell, a: f64) {
    match shell {
        Shell::Zero => drift_shell0(r, eos, a),
        Shell::One => drift_shell1(r, eos, a),
    }
}

fn interaction_step(r: &mut reb_simulation, eos: &reb_integrator_eos_state, shell: Shell, y: f64, v: f64) {
    match shell {
        Shell::Zero => interaction_shell0(r, eos, y, v),
        Shell::One => interaction_shell1(r, eos, y, v),
    }
}

/// integrator_eos.c `reb_integrator_eos_preprocessor`.
fn preprocessor(r: &mut reb_simulation, eos: &reb_integrator_eos_state, dt: f64, type_: i32, shell: Shell) {
    match type_ {
        REB_INTEGRATOR_EOS_TYPE_PMLF6 => {
            for i in 0..6 {
                drift_step(r, eos, shell, dt * pmlf6_z[i]);
                interaction_step(r, eos, shell, dt * pmlf6_y[i], dt * dt * dt * pmlf6_v[i]);
            }
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF4 => {
            for i in 0..3 {
                interaction_step(r, eos, shell, dt * pmlf4_y[i], 0.);
                drift_step(r, eos, shell, dt * pmlf4_z[i]);
            }
        }
        REB_INTEGRATOR_EOS_TYPE_PLF7_6_4 => {
            for i in 0..6 {
                drift_step(r, eos, shell, dt * plf7_6_4_z[i]);
                interaction_step(r, eos, shell, dt * plf7_6_4_y[i], 0.);
            }
        }
        _ => {}
    }
}

/// integrator_eos.c `reb_integrator_eos_postprocessor`.
fn postprocessor(r: &mut reb_simulation, eos: &reb_integrator_eos_state, dt: f64, type_: i32, shell: Shell) {
    match type_ {
        REB_INTEGRATOR_EOS_TYPE_PMLF6 => {
            for i in (0..6).rev() {
                interaction_step(r, eos, shell, -dt * pmlf6_y[i], -dt * dt * dt * pmlf6_v[i]);
                drift_step(r, eos, shell, -dt * pmlf6_z[i]);
            }
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF4 => {
            for i in (0..3).rev() {
                drift_step(r, eos, shell, -dt * pmlf4_z[i]);
                interaction_step(r, eos, shell, -dt * pmlf4_y[i], 0.);
            }
        }
        REB_INTEGRATOR_EOS_TYPE_PLF7_6_4 => {
            for i in (0..6).rev() {
                interaction_step(r, eos, shell, -dt * plf7_6_4_y[i], 0.);
                drift_step(r, eos, shell, -dt * plf7_6_4_z[i]);
            }
        }
        _ => {}
    }
}

/// integrator_eos.c `reb_integrator_eos_drift_shell0` — one outer
/// drift = a full inner phi1 integration over n sub-steps.
fn drift_shell0(r: &mut reb_simulation, eos: &reb_integrator_eos_state, _dt: f64) {
    let n = eos.n as usize;
    let dt = _dt / (n as f64);
    preprocessor(r, eos, dt, eos.phi1, Shell::One);
    match eos.phi1 {
        REB_INTEGRATOR_EOS_TYPE_LF => {
            drift_shell1(r, eos, dt * 0.5);
            for i in 0..n {
                interaction_shell1(r, eos, dt, 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, dt);
                }
            }
            drift_shell1(r, eos, dt * 0.5);
        }
        REB_INTEGRATOR_EOS_TYPE_LF4 => {
            let a = reb_integrator_leapfrog_lf4_a;
            drift_shell1(r, eos, dt * a);
            for i in 0..n {
                interaction_shell1(r, eos, dt * 2. * a, 0.);
                drift_shell1(r, eos, dt * (0.5 - a));
                interaction_shell1(r, eos, dt * (1. - 4. * a), 0.);
                drift_shell1(r, eos, dt * (0.5 - a));
                interaction_shell1(r, eos, dt * 2. * a, 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, dt * 2. * a);
                }
            }
            drift_shell1(r, eos, dt * a);
        }
        REB_INTEGRATOR_EOS_TYPE_LF6 => {
            let a = &reb_integrator_leapfrog_lf6_a;
            drift_shell1(r, eos, dt * a[0] * 0.5);
            for i in 0..n {
                interaction_shell1(r, eos, dt * a[0], 0.);
                drift_shell1(r, eos, dt * (a[0] + a[1]) * 0.5);
                interaction_shell1(r, eos, dt * a[1], 0.);
                drift_shell1(r, eos, dt * (a[1] + a[2]) * 0.5);
                interaction_shell1(r, eos, dt * a[2], 0.);
                drift_shell1(r, eos, dt * (a[2] + a[3]) * 0.5);
                interaction_shell1(r, eos, dt * a[3], 0.);
                drift_shell1(r, eos, dt * (a[3] + a[4]) * 0.5);
                interaction_shell1(r, eos, dt * a[4], 0.);
                drift_shell1(r, eos, dt * (a[3] + a[4]) * 0.5);
                interaction_shell1(r, eos, dt * a[3], 0.);
                drift_shell1(r, eos, dt * (a[2] + a[3]) * 0.5);
                interaction_shell1(r, eos, dt * a[2], 0.);
                drift_shell1(r, eos, dt * (a[1] + a[2]) * 0.5);
                interaction_shell1(r, eos, dt * a[1], 0.);
                drift_shell1(r, eos, dt * (a[0] + a[1]) * 0.5);
                interaction_shell1(r, eos, dt * a[0], 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, dt * a[0]);
                }
            }
            drift_shell1(r, eos, dt * a[0] * 0.5);
        }
        REB_INTEGRATOR_EOS_TYPE_LF8 => {
            let a = &reb_integrator_leapfrog_lf8_a;
            drift_shell1(r, eos, dt * a[0] * 0.5);
            for i in 0..n {
                interaction_shell1(r, eos, dt * a[0], 0.);
                drift_shell1(r, eos, dt * (a[0] + a[1]) * 0.5);
                interaction_shell1(r, eos, dt * a[1], 0.);
                drift_shell1(r, eos, dt * (a[1] + a[2]) * 0.5);
                interaction_shell1(r, eos, dt * a[2], 0.);
                drift_shell1(r, eos, dt * (a[2] + a[3]) * 0.5);
                interaction_shell1(r, eos, dt * a[3], 0.);
                drift_shell1(r, eos, dt * (a[3] + a[4]) * 0.5);
                interaction_shell1(r, eos, dt * a[4], 0.);
                drift_shell1(r, eos, dt * (a[4] + a[5]) * 0.5);
                interaction_shell1(r, eos, dt * a[5], 0.);
                drift_shell1(r, eos, dt * (a[5] + a[6]) * 0.5);
                interaction_shell1(r, eos, dt * a[6], 0.);
                drift_shell1(r, eos, dt * (a[6] + a[7]) * 0.5);
                interaction_shell1(r, eos, dt * a[7], 0.);
                drift_shell1(r, eos, dt * (a[7] + a[8]) * 0.5);
                interaction_shell1(r, eos, dt * a[8], 0.);
                drift_shell1(r, eos, dt * (a[7] + a[8]) * 0.5);
                interaction_shell1(r, eos, dt * a[7], 0.);
                drift_shell1(r, eos, dt * (a[6] + a[7]) * 0.5);
                interaction_shell1(r, eos, dt * a[6], 0.);
                drift_shell1(r, eos, dt * (a[5] + a[6]) * 0.5);
                interaction_shell1(r, eos, dt * a[5], 0.);
                drift_shell1(r, eos, dt * (a[4] + a[5]) * 0.5);
                interaction_shell1(r, eos, dt * a[4], 0.);
                drift_shell1(r, eos, dt * (a[3] + a[4]) * 0.5);
                interaction_shell1(r, eos, dt * a[3], 0.);
                drift_shell1(r, eos, dt * (a[2] + a[3]) * 0.5);
                interaction_shell1(r, eos, dt * a[2], 0.);
                drift_shell1(r, eos, dt * (a[1] + a[2]) * 0.5);
                interaction_shell1(r, eos, dt * a[1], 0.);
                drift_shell1(r, eos, dt * (a[0] + a[1]) * 0.5);
                interaction_shell1(r, eos, dt * a[0], 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, dt * a[0]);
                }
            }
            drift_shell1(r, eos, dt * a[0] * 0.5);
        }
        REB_INTEGRATOR_EOS_TYPE_LF4_2 => {
            drift_shell1(r, eos, dt * lf4_2_a);
            for i in 0..n {
                interaction_shell1(r, eos, dt * 0.5, 0.);
                drift_shell1(r, eos, dt * (1. - 2. * lf4_2_a));
                interaction_shell1(r, eos, dt * 0.5, 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, 2. * dt * lf4_2_a);
                }
            }
            drift_shell1(r, eos, dt * lf4_2_a);
        }
        REB_INTEGRATOR_EOS_TYPE_LF8_6_4 => {
            drift_shell1(r, eos, dt * lf8_6_4_a[0]);
            for i in 0..n {
                interaction_shell1(r, eos, lf8_6_4_b[0] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[1] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[1] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[2] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[2] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[3] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[3] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[3] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[2] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[2] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[1] * dt, 0.);
                drift_shell1(r, eos, lf8_6_4_a[1] * dt);
                interaction_shell1(r, eos, lf8_6_4_b[0] * dt, 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, 2. * dt * lf8_6_4_a[0]);
                }
            }
            drift_shell1(r, eos, dt * lf8_6_4_a[0]);
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF4 => {
            drift_shell1(r, eos, dt * 0.5);
            for i in 0..n {
                interaction_shell1(r, eos, dt, dt * dt * dt / 24.);
                if i < n - 1 {
                    drift_shell1(r, eos, dt);
                }
            }
            drift_shell1(r, eos, dt * 0.5);
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF6 => {
            drift_shell1(r, eos, dt * pmlf6_a[0]);
            for i in 0..n {
                interaction_shell1(r, eos, dt * pmlf6_b[0], dt * dt * dt * pmlf6_c[0]);
                drift_shell1(r, eos, dt * pmlf6_a[1]);
                interaction_shell1(r, eos, dt * pmlf6_b[1], dt * dt * dt * pmlf6_c[1]);
                drift_shell1(r, eos, dt * pmlf6_a[1]);
                interaction_shell1(r, eos, dt * pmlf6_b[0], dt * dt * dt * pmlf6_c[0]);
                if i < n - 1 {
                    drift_shell1(r, eos, 2. * dt * pmlf6_a[0]);
                }
            }
            drift_shell1(r, eos, dt * pmlf6_a[0]);
        }
        REB_INTEGRATOR_EOS_TYPE_PLF7_6_4 => {
            drift_shell1(r, eos, dt * plf7_6_4_a[0]);
            for i in 0..n {
                interaction_shell1(r, eos, plf7_6_4_b[0] * dt, 0.);
                drift_shell1(r, eos, plf7_6_4_a[1] * dt);
                interaction_shell1(r, eos, plf7_6_4_b[1] * dt, 0.);
                drift_shell1(r, eos, plf7_6_4_a[1] * dt);
                interaction_shell1(r, eos, plf7_6_4_b[0] * dt, 0.);
                if i < n - 1 {
                    drift_shell1(r, eos, 2. * dt * plf7_6_4_a[0]);
                }
            }
            drift_shell1(r, eos, dt * plf7_6_4_a[0]);
        }
        _ => {}
    }
    postprocessor(r, eos, dt, eos.phi1, Shell::One);
}

/// integrator_eos.c `reb_integrator_eos_step` (state-explicit).
pub fn reb_integrator_eos_step_state(r: &mut reb_simulation, eos: &mut reb_integrator_eos_state) {
    if r.gravity != REB_GRAVITY::BASIC {
        reb_simulation_warning(r, "EOS only supports the BASIC gravity routine.");
    }
    if r.N_var != 0 {
        reb_simulation_warning(r, "Variational particles/MEGNO in EOS no longer supported since REBOUND version 5.");
    }
    r.gravity = REB_GRAVITY::NONE;

    let dt = r.dt;
    let eosc = *eos;

    let mut dtfac = 1.;
    if r.is_synchronized != 0 {
        preprocessor(r, &eosc, r.dt, eosc.phi0, Shell::Zero);
    } else {
        dtfac = 2.;
    }
    match eosc.phi0 {
        REB_INTEGRATOR_EOS_TYPE_LF => {
            drift_shell0(r, &eosc, dt * 0.5 * dtfac);
            interaction_shell0(r, &eosc, dt, 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_LF4 => {
            let a = reb_integrator_leapfrog_lf4_a;
            drift_shell0(r, &eosc, dt * a * dtfac);
            interaction_shell0(r, &eosc, dt * 2. * a, 0.);
            drift_shell0(r, &eosc, dt * (0.5 - a));
            interaction_shell0(r, &eosc, dt * (1. - 4. * a), 0.);
            drift_shell0(r, &eosc, dt * (0.5 - a));
            interaction_shell0(r, &eosc, dt * 2. * a, 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_LF6 => {
            let a = &reb_integrator_leapfrog_lf6_a;
            drift_shell0(r, &eosc, dt * a[0] * 0.5 * dtfac);
            interaction_shell0(r, &eosc, dt * a[0], 0.);
            drift_shell0(r, &eosc, dt * (a[0] + a[1]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[1], 0.);
            drift_shell0(r, &eosc, dt * (a[1] + a[2]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[2], 0.);
            drift_shell0(r, &eosc, dt * (a[2] + a[3]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[3], 0.);
            drift_shell0(r, &eosc, dt * (a[3] + a[4]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[4], 0.);
            drift_shell0(r, &eosc, dt * (a[3] + a[4]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[3], 0.);
            drift_shell0(r, &eosc, dt * (a[2] + a[3]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[2], 0.);
            drift_shell0(r, &eosc, dt * (a[1] + a[2]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[1], 0.);
            drift_shell0(r, &eosc, dt * (a[0] + a[1]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[0], 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_LF8 => {
            let a = &reb_integrator_leapfrog_lf8_a;
            drift_shell0(r, &eosc, dt * a[0] * 0.5 * dtfac);
            interaction_shell0(r, &eosc, dt * a[0], 0.);
            drift_shell0(r, &eosc, dt * (a[0] + a[1]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[1], 0.);
            drift_shell0(r, &eosc, dt * (a[1] + a[2]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[2], 0.);
            drift_shell0(r, &eosc, dt * (a[2] + a[3]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[3], 0.);
            drift_shell0(r, &eosc, dt * (a[3] + a[4]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[4], 0.);
            drift_shell0(r, &eosc, dt * (a[4] + a[5]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[5], 0.);
            drift_shell0(r, &eosc, dt * (a[5] + a[6]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[6], 0.);
            drift_shell0(r, &eosc, dt * (a[6] + a[7]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[7], 0.);
            drift_shell0(r, &eosc, dt * (a[7] + a[8]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[8], 0.);
            drift_shell0(r, &eosc, dt * (a[7] + a[8]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[7], 0.);
            drift_shell0(r, &eosc, dt * (a[6] + a[7]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[6], 0.);
            drift_shell0(r, &eosc, dt * (a[5] + a[6]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[5], 0.);
            drift_shell0(r, &eosc, dt * (a[4] + a[5]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[4], 0.);
            drift_shell0(r, &eosc, dt * (a[3] + a[4]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[3], 0.);
            drift_shell0(r, &eosc, dt * (a[2] + a[3]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[2], 0.);
            drift_shell0(r, &eosc, dt * (a[1] + a[2]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[1], 0.);
            drift_shell0(r, &eosc, dt * (a[0] + a[1]) * 0.5);
            interaction_shell0(r, &eosc, dt * a[0], 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_LF4_2 => {
            drift_shell0(r, &eosc, dt * lf4_2_a * dtfac);
            interaction_shell0(r, &eosc, dt * 0.5, 0.);
            drift_shell0(r, &eosc, dt * (1. - 2. * lf4_2_a));
            interaction_shell0(r, &eosc, dt * 0.5, 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_LF8_6_4 => {
            drift_shell0(r, &eosc, dt * lf8_6_4_a[0] * dtfac);
            interaction_shell0(r, &eosc, lf8_6_4_b[0] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[1] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[1] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[2] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[2] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[3] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[3] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[3] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[2] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[2] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[1] * dt, 0.);
            drift_shell0(r, &eosc, lf8_6_4_a[1] * dt);
            interaction_shell0(r, &eosc, lf8_6_4_b[0] * dt, 0.);
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF4 => {
            drift_shell0(r, &eosc, dt * 0.5 * dtfac);
            interaction_shell0(r, &eosc, dt, dt * dt * dt / 24.);
        }
        REB_INTEGRATOR_EOS_TYPE_PMLF6 => {
            drift_shell0(r, &eosc, dt * pmlf6_a[0] * dtfac);
            interaction_shell0(r, &eosc, dt * pmlf6_b[0], dt * dt * dt * pmlf6_c[0]);
            drift_shell0(r, &eosc, dt * pmlf6_a[1]);
            interaction_shell0(r, &eosc, dt * pmlf6_b[1], dt * dt * dt * pmlf6_c[1]);
            drift_shell0(r, &eosc, dt * pmlf6_a[1]);
            interaction_shell0(r, &eosc, dt * pmlf6_b[0], dt * dt * dt * pmlf6_c[0]);
        }
        REB_INTEGRATOR_EOS_TYPE_PLF7_6_4 => {
            drift_shell0(r, &eosc, dt * plf7_6_4_a[0] * dtfac);
            interaction_shell0(r, &eosc, plf7_6_4_b[0] * dt, 0.);
            drift_shell0(r, &eosc, plf7_6_4_a[1] * dt);
            interaction_shell0(r, &eosc, plf7_6_4_b[1] * dt, 0.);
            drift_shell0(r, &eosc, plf7_6_4_a[1] * dt);
            interaction_shell0(r, &eosc, plf7_6_4_b[0] * dt, 0.);
        }
        _ => {}
    }

    r.is_synchronized = 0;
    if eosc.safe_mode != 0 {
        reb_integrator_eos_synchronize_state(r, eos);
    }

    r.t += r.dt;
    r.dt_last_done = r.dt;
}

/// integrator_eos.c `reb_integrator_eos_synchronize` (state-explicit).
pub fn reb_integrator_eos_synchronize_state(r: &mut reb_simulation, eos: &mut reb_integrator_eos_state) {
    let dt = r.dt;
    let eosc = *eos;
    if r.is_synchronized == 0 {
        match eosc.phi0 {
            REB_INTEGRATOR_EOS_TYPE_PMLF4 | REB_INTEGRATOR_EOS_TYPE_LF => {
                drift_shell0(r, &eosc, dt * 0.5);
            }
            REB_INTEGRATOR_EOS_TYPE_PMLF6 => {
                drift_shell0(r, &eosc, dt * pmlf6_a[0]);
            }
            REB_INTEGRATOR_EOS_TYPE_LF4 => {
                drift_shell0(r, &eosc, dt * reb_integrator_leapfrog_lf4_a);
            }
            REB_INTEGRATOR_EOS_TYPE_LF4_2 => {
                drift_shell0(r, &eosc, dt * lf4_2_a);
            }
            REB_INTEGRATOR_EOS_TYPE_PLF7_6_4 => {
                drift_shell0(r, &eosc, dt * plf7_6_4_a[0]);
            }
            REB_INTEGRATOR_EOS_TYPE_LF8_6_4 => {
                drift_shell0(r, &eosc, dt * lf8_6_4_a[0]);
            }
            REB_INTEGRATOR_EOS_TYPE_LF6 => {
                drift_shell0(r, &eosc, dt * reb_integrator_leapfrog_lf6_a[0] * 0.5);
            }
            REB_INTEGRATOR_EOS_TYPE_LF8 => {
                drift_shell0(r, &eosc, dt * reb_integrator_leapfrog_lf8_a[0] * 0.5);
            }
            _ => {}
        }
        postprocessor(r, &eosc, r.dt, eosc.phi0, Shell::Zero);
        r.is_synchronized = 1;
    }
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_eos_step(r: &mut reb_simulation) {
    let mut eos = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::eos(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_eos_step_state(r, &mut eos);
    r.integrator = reb_integrator_state::eos(eos);
}

/// Synchronize entry point for the dispatcher.
pub fn reb_integrator_eos_synchronize(r: &mut reb_simulation) {
    let mut eos = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::eos(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_eos_synchronize_state(r, &mut eos);
    r.integrator = reb_integrator_state::eos(eos);
}
