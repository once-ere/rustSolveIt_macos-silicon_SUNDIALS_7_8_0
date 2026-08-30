#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

//! Independent cross-check (verification gate 2.8): the Mercury tidal despin
//! rate from the workspace's pure-Rust REBOUNDx port (`tides_spin`, the
//! Eggleton 1998 / Lu 2023 constant-time-lag model integrated as a real
//! N-body force) must agree with the Hut (1981) orbit-averaged secular
//! formula that mercury_rs integrates. Two entirely different codes, two
//! entirely different formulations, one physical answer.
//!
//! GPL-3.0-or-later (links the GPL rebound/reboundx ports). Standalone —
//! shares no code with mercury_rs; the Hut formula is re-derived inline.

use std::f64::consts::PI;

use rebound_rs::*;
use reboundx_rs::*;

const G_SI: f64 = 6.67430e-11;
const M_SUN: f64 = 1.98847e30;
const M_MERC: f64 = 3.3011e23;
const R_MERC: f64 = 2.4397e6;
const A0: f64 = 5.790905e10;
const E0: f64 = 0.20563;
const C_FACTOR: f64 = 0.34;
const K2: f64 = 0.12;
/// The movie (S = 1000) HUT time lag [s].
const TAU_HUT: f64 = 1.0e5;
/// REBOUNDx tides_spin defines its lag through sigma = 4*tau*G/(3*R^5*k2)
/// (Eggleton 1998 convention). Deriving the circular-orbit spin torque by
/// hand and matching it to Hut's 3*G*Msun^2*R^5*k2*tau/a^6*(Omega - n) shows
/// the two "tau"s differ by EXACTLY a factor of two:
///     tau_reboundx = tau_Hut / 2.
/// (First run of this cross-check measured precisely 2.005x before the
/// mapping — the 0.5% being within-orbit sampling — which is exactly the
/// kind of convention trap a cross-check exists to catch.)
const TAU_REBX: f64 = TAU_HUT / 2.0;
const SPIN_RATIO: f64 = 100.0; // fast-spin era, far from every resonance
const ORBITS: f64 = 200.0;

fn hut_f1(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let om = 1.0 - e2;
    let om2 = om * om;
    (1.0 + 3.0 * e2 + 0.375 * e4) / (om2 * om2 * om.sqrt())
}

fn hut_f2(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let e6 = e4 * e2;
    let om = 1.0 - e2;
    let om2 = om * om;
    (1.0 + 7.5 * e2 + 5.625 * e4 + 0.3125 * e6) / (om2 * om2 * om2)
}

fn main() {
    println!("mercury_crosscheck: REBOUNDx tides_spin vs the Hut (1981) secular despin rate");
    let mut sim = reb_simulation_create();
    sim.G = G_SI;
    reb_simulation_add_fmt(&mut sim, "m", &[reb_fmt_arg::d(M_SUN)]);
    reb_simulation_add_fmt(
        &mut sim,
        "m a e r",
        &[
            reb_fmt_arg::d(M_MERC),
            reb_fmt_arg::d(A0),
            reb_fmt_arg::d(E0),
            reb_fmt_arg::d(R_MERC),
        ],
    );
    sim.N_active = 2;
    reb_simulation_set_integrator(&mut sim, "ias15");
    let n0 = (G_SI * (M_SUN + M_MERC) / (A0 * A0 * A0)).sqrt();
    let p_orb = 2.0 * PI / n0;
    sim.dt = p_orb / 100.0;

    rebx_attach(&mut sim);
    let effect = match rebx_load_force(&mut sim, "tides_spin") {
        Some(e) => e,
        None => {
            println!("FAIL - crosscheck.setup: the tides_spin force is missing from reboundx_rs");
            println!("FAILURE");
            std::process::exit(1);
        }
    };
    rebx_add_force(&mut sim, effect);
    let omega0 = SPIN_RATIO * n0;
    let inertia = C_FACTOR * M_MERC * (R_MERC * R_MERC);
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "k2", K2);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau", TAU_REBX);
        rebx_set_param_vec3d(
            rebx,
            rebx_ap::particle(1),
            "Omega",
            reb_vec3d {
                x: 0.0,
                y: 0.0,
                z: omega0,
            },
        );
        rebx_set_param_double(rebx, rebx_ap::particle(1), "I", inertia);
    }
    reb_simulation_move_to_com(&mut sim);
    rebx_with(&mut sim, |sim, rebx| {
        rebx_spin_initialize_ode(sim, rebx, effect);
    });

    let t_end = ORBITS * p_orb;
    reb_simulation_integrate(&mut sim, t_end);
    if sim.t < 0.99 * t_end {
        println!(
            "FAIL - crosscheck.integration_incomplete: reached t = {:.3e} of {:.3e} s",
            sim.t, t_end
        );
        println!("FAILURE");
        std::process::exit(1);
    }

    let om_end = match rebx_extras_ref(&sim)
        .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(1), "Omega"))
    {
        Some(v) => v,
        None => {
            println!("FAIL - crosscheck.readback: Mercury's Omega parameter could not be read back");
            println!("FAILURE");
            std::process::exit(1);
        }
    };
    let measured = (om_end.z - omega0) / sim.t;

    // Hut (1981) secular rate at the initial state (Omega changes < 0.02%
    // over the window, so the initial-state rate is the right comparator):
    let a2 = A0 * A0;
    let a6 = a2 * a2 * a2;
    let r2 = R_MERC * R_MERC;
    let r5 = r2 * r2 * R_MERC;
    let k_brake = 3.0 * G_SI * (M_SUN * M_SUN) * r5 * (K2 * TAU_HUT) / a6;
    let predicted = -k_brake * (omega0 * hut_f1(E0) - n0 * hut_f2(E0)) / inertia;

    let rel = (measured - predicted).abs() / predicted.abs();
    println!("  spin ratio {SPIN_RATIO} x n, e = {E0}, {ORBITS} orbits with IAS15 + tides_spin");
    println!("  measured  dOmega/dt = {measured:+.6e} rad/s^2  (REBOUNDx tides_spin)");
    println!("  predicted dOmega/dt = {predicted:+.6e} rad/s^2  (Hut 1981 secular formula)");
    println!("  relative difference = {rel:.3e}");
    let ok = rel < 0.02;
    println!(
        "{} - crosscheck.despin_rate: two independent formulations agree within 2%",
        if ok { "PASS" } else { "FAIL" }
    );
    if ok {
        println!("SUCCESS");
    } else {
        println!("FAILURE");
        std::process::exit(1);
    }
}
