#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

//! Analytic-expectation tests (Gate 1 of the verification plan). Every
//! expectation is a number derivable with pencil and paper — never a golden
//! snapshot of previous program output.

use std::f64::consts::PI;

use mercury_rs::driver::{
    integrate_segment, wrap_pi, k2tau_movie, ObserverCmd, Sample, SegmentSpec, State, Stats,
};
use mercury_rs::hut;
use mercury_rs::kepler;
use mercury_rs::output;
use mercury_rs::params;

fn noop() -> impl FnMut(&Sample) -> ObserverCmd {
    |_s: &Sample| ObserverCmd::Continue
}

/// 1.1 — n(a0) from Kepler's third law and the resulting orbital period.
#[test]
fn mean_motion_matches_keplers_third_law() {
    let n = params::mean_motion(params::A0);
    assert!(
        (n - 8.2669e-7).abs() / 8.2669e-7 < 1.0e-4,
        "n(a0) = {n}, expected 8.2669e-7 rad/s"
    );
    let p_orb_days = 2.0 * PI / n / 86400.0;
    assert!(
        (p_orb_days - 87.969).abs() < 0.01,
        "P_orb = {p_orb_days} d, expected ~87.969 d"
    );
}

/// 1.2 — solving M = E - e sinE back recovers E; f and r match closed forms.
#[test]
fn kepler_solver_inverts_keplers_equation() {
    for &e in &[0.0, 0.05, 0.20563, 0.285, 0.6] {
        for k in 0..24 {
            let ecc_true = -PI + (2.0 * PI) * (k as f64 + 0.5) / 24.0;
            let m = ecc_true - e * ecc_true.sin();
            let sol = kepler::solve(m, e, params::A0).expect("kepler::solve failed");
            assert!(
                (sol.ecc_anom - ecc_true).abs() < 1.0e-12,
                "e = {e}, E_true = {ecc_true}: solver E = {}",
                sol.ecc_anom
            );
            let r_from_E = params::A0 * (1.0 - e * sol.ecc_anom.cos());
            assert!(
                (sol.radius - r_from_E).abs() / r_from_E < 1.0e-12,
                "radius mismatch at e = {e}, E = {ecc_true}"
            );
            let tan_half_f = ((1.0 + e) / (1.0 - e)).sqrt() * (ecc_true / 2.0).tan();
            let f_closed = 2.0 * tan_half_f.atan();
            assert!(
                wrap_pi(sol.true_anom - f_closed).abs() < 1.0e-10,
                "true anomaly mismatch at e = {e}, E = {ecc_true}"
            );
        }
    }
}

/// 1.3 — at e = 0: E = M, f = M, r = a exactly.
#[test]
fn kepler_solver_circular_orbit_is_identity() {
    for &m in &[0.0, 0.7, -1.3, 3.0] {
        let sol = kepler::solve(m, 0.0, params::A0).expect("kepler::solve failed");
        assert!((sol.ecc_anom - m).abs() < 1.0e-14);
        assert!(wrap_pi(sol.true_anom - m).abs() < 1.0e-13);
        assert!((sol.radius - params::A0).abs() / params::A0 < 1.0e-15);
    }
}

/// 1.4 — f1(0) = f2(0) = f3(0) = f4(0) = f5(0) = 1 exactly.
#[test]
fn hut_polynomials_at_zero_eccentricity_are_one() {
    assert_eq!(hut::f1(0.0), 1.0);
    assert_eq!(hut::f2(0.0), 1.0);
    assert_eq!(hut::f3(0.0), 1.0);
    assert_eq!(hut::f4(0.0), 1.0);
    assert_eq!(hut::f5(0.0), 1.0);
}

/// 1.5 — the pseudo-synchronous ratio at Mercury's eccentricity is 1.256.
#[test]
fn pseudo_synchronous_ratio_at_mercury_e() {
    let r = hut::pseudo_synchronous_ratio(params::E0);
    assert!(
        (r - 1.2560).abs() < 0.0005,
        "f2/f1 at e = 0.20563 is {r}, expected 1.2560 +/- 0.0005"
    );
}

/// 1.6 — the pseudo-synchronous ratio is exactly 3/2 at e = 0.285.
#[test]
fn pseudo_synchronous_equals_three_halves_at_e285() {
    let r = hut::pseudo_synchronous_ratio(0.285);
    assert!(
        (r - 1.5000).abs() < 0.0005,
        "f2/f1 at e = 0.285 is {r}, expected 1.5000 +/- 0.0005"
    );
}

/// 1.7 — the brake slows a fast spin and spins up a slow one
/// (Omega_eq = 1.256 n lies above n).
#[test]
fn tidal_torque_sign_brakes_fast_spin_and_spins_up_slow() {
    let n = params::mean_motion(params::A0);
    let k = params::tidal_k(params::A0, k2tau_movie());
    assert!(hut::tidal_torque(k, 2.0 * n, n, params::E0) < 0.0);
    assert!(hut::tidal_torque(k, 1.0 * n, n, params::E0) > 0.0);
}

/// 1.8 — a short secular CVODE run despins at the hand-computed rate.
#[test]
fn despin_rate_matches_linearized_prediction() {
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, params::E0, 0.0, 0.0, 100.0 * n0],
    };
    let t_end = 1000.0 * params::YEAR;
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let spec = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: false,
        t_end,
        cadence: 50.0 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let mut obs = noop();
    integrate_segment(&ic, &spec, &mut obs, &mut samples, &mut stats).expect("segment failed");
    let last = samples.last().expect("no samples");
    let measured = (last.omega - ic.y[4]) / last.t;
    let mid = &samples[samples.len() / 2];
    let k = params::tidal_k(mid.a, k2tau_movie());
    let predicted =
        -k * (mid.omega * hut::f1(mid.e) - mid.n * hut::f2(mid.e)) / params::moment_of_inertia();
    let rel = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel < 1.0e-3,
        "despin rate: measured {measured}, predicted {predicted}, rel diff {rel}"
    );
}

/// 1.9 — the corrected Hut trio conserves C*Omega + L_orb EXACTLY at the
/// derivative level (the uncorrected source-spec da/dt fails wildly).
#[test]
fn angular_momentum_ledger_closes_in_secular_stage() {
    let states: [(f64, f64, f64); 4] = [
        (params::A0, params::E0, 181.4),
        (params::A0, 0.10, 2.0),
        (1.1 * params::A0, 0.30, 1.7),
        (params::A0, params::E0, 0.5),
    ];
    for &(a, e, ratio) in &states {
        let n = params::mean_motion(a);
        let omega = ratio * n;
        let k = params::tidal_k(a, k2tau_movie());
        let m = params::M_MERCURY;
        let da = (2.0 * k / (m * n * a)) * (omega * hut::f2(e) - n * hut::f3(e));
        let de = (9.0 * k * e / (m * n * (a * a)))
            * ((11.0 / 18.0) * omega * hut::f4(e) - n * hut::f5(e));
        let dom = -k * (omega * hut::f1(e) - n * hut::f2(e)) / params::moment_of_inertia();
        // dL_orb/dt with L_orb = m n a^2 sqrt(1-e^2) and dn/da = -(3/2) n/a:
        let sq = (1.0 - e * e).sqrt();
        let dl_orb = m * (0.5 * n * a * sq * da - n * (a * a) * e * de / sq);
        let dl_spin = params::moment_of_inertia() * dom;
        let closure = (dl_spin + dl_orb).abs() / dl_spin.abs();
        assert!(
            closure < 1.0e-12,
            "ledger closure at (a={a}, e={e}, ratio={ratio}): rel residual {closure}"
        );
        // The SOURCE spec's uncorrected da/dt (wrong prefactor AND sign) breaks it:
        let da_wrong = -(2.0 * a * k / (params::G * params::M_SUN * m))
            * (omega * hut::f2(e) - n * hut::f3(e));
        let dl_orb_wrong = m * (0.5 * n * a * sq * da_wrong - n * (a * a) * e * de / sq);
        let closure_wrong = (dl_spin + dl_orb_wrong).abs() / dl_spin.abs();
        assert!(
            closure_wrong > 1.0e-3,
            "the uncorrected spec equation unexpectedly closed the ledger ({closure_wrong})"
        );
    }
}

/// 1.10 — re-anchoring leaves gamma = 2 theta - 3 M unchanged.
#[test]
fn reanchoring_preserves_gamma_exactly() {
    let st = State {
        t: 1.0e12,
        y: [
            params::A0,
            params::E0,
            4321.75,
            6482.6255,
            1.3 * params::mean_motion(params::A0),
        ],
    };
    let before = wrap_pi(2.0 * st.y[3] - 3.0 * st.y[2]);
    let re = st.reanchored();
    let after = wrap_pi(2.0 * re.y[3] - 3.0 * re.y[2]);
    assert!(re.y[2] >= 0.0 && re.y[2] < 2.0 * PI, "M not reduced: {}", re.y[2]);
    assert!(
        wrap_pi(after - before).abs() < 1.0e-9,
        "gamma changed by re-anchoring: {before} -> {after}"
    );
}

/// 1.11 — the handle torque vanishes when the long axis points at the Sun and
/// peaks at 45 degrees off.
#[test]
fn triaxial_torque_zero_when_aligned() {
    let r = params::A0;
    assert_eq!(hut::triaxial_torque(0.7, 0.7, r), 0.0);
    let peak = hut::triaxial_torque(PI / 4.0, 0.0, r);
    let amplitude = 1.5 * params::G * params::M_SUN * params::b_minus_a() / (r * r * r);
    assert!(
        (peak.abs() - amplitude).abs() / amplitude < 1.0e-14,
        "peak torque {peak}, expected magnitude {amplitude}"
    );
    assert!(peak < 0.0, "restoring sign: at theta - f = +45 deg the torque must be negative");
}

/// 1.12 — restart rows (fmt_e 17) round-trip bit-exactly; sample rows
/// (fmt_e 12) round-trip to 1e-12 relative.
#[test]
fn csv_row_roundtrips_through_fmt_e() {
    use sundials_core::sundials_utils::fmt_e;
    let x = 5.790905e10_f64 * (1.0 + 1.0e-13) / 3.0;
    let x17: f64 = fmt_e(x, 17).trim().parse().expect("parse");
    assert_eq!(x17.to_bits(), x.to_bits(), "%.17e must round-trip bit-exactly");
    let x12: f64 = fmt_e(x, 12).trim().parse().expect("parse");
    assert!(((x12 - x) / x).abs() < 1.0e-12);
}

/// 1.13 — restart.csv reload is bit-exact and continuation stays consistent.
#[test]
fn restart_reproduces_bit_identical_state() {
    let dir = std::env::temp_dir().join("mercury_rs_test_restart");
    std::fs::create_dir_all(&dir).expect("temp dir");
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, params::E0, 1.234, 2.345, 50.0 * n0],
    };
    let spec = |t_end: f64| SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: false,
        t_end,
        cadence: 10.0 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: false,
    };
    let mut v: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let mut obs = noop();
    let half = integrate_segment(&ic, &spec(100.0 * params::YEAR), &mut obs, &mut v, &mut stats)
        .expect("first half failed");
    output::write_restart(&dir, &half.state).expect("write restart");
    let reloaded = output::read_restart(&dir).expect("read restart");
    assert_eq!(reloaded.t.to_bits(), half.state.t.to_bits());
    for i in 0..5 {
        assert_eq!(
            reloaded.y[i].to_bits(),
            half.state.y[i].to_bits(),
            "restart component {i} not bit-exact"
        );
    }
    let mut obs2 = noop();
    let cont =
        integrate_segment(&reloaded, &spec(200.0 * params::YEAR), &mut obs2, &mut v, &mut stats)
            .expect("continuation failed");
    let mut obs3 = noop();
    let straight =
        integrate_segment(&ic, &spec(200.0 * params::YEAR), &mut obs3, &mut v, &mut stats)
            .expect("straight run failed");
    for i in 0..5 {
        let scale = straight.state.y[i].abs().max(1.0e-30);
        assert!(
            ((cont.state.y[i] - straight.state.y[i]) / scale).abs() < 1.0e-8,
            "continuation drifted in component {i}: {} vs {}",
            cont.state.y[i],
            straight.state.y[i]
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// 1.14 — with tides off, gamma librates at the Goldreich-Peale frequency
/// omega_lib = n sqrt(3 (B-A)/C H(e)), H = (7/2) e - (123/16) e^3.
#[test]
fn libration_frequency_matches_goldreich_peale() {
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, params::E0, 0.0, 0.05, 1.5 * n0],
    };
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let spec = SegmentSpec {
        k2tau: 0.0, // clean pendulum: handle torque only
        triaxial_on: true,
        t_end: 60.0 * params::YEAR,
        cadence: 0.05 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'R',
        record: true,
    };
    let mut obs = noop();
    integrate_segment(&ic, &spec, &mut obs, &mut samples, &mut stats).expect("segment failed");
    let meang: f64 = samples.iter().map(|s| s.gamma).sum::<f64>() / (samples.len() as f64);
    let mut crossings: Vec<f64> = Vec::new();
    for pair in samples.windows(2) {
        if (pair[0].gamma - meang) <= 0.0 && (pair[1].gamma - meang) > 0.0 {
            crossings.push(pair[1].t);
        }
    }
    assert!(crossings.len() >= 3, "too few libration cycles observed: {}", crossings.len());
    let measured =
        (crossings[crossings.len() - 1] - crossings[0]) / ((crossings.len() - 1) as f64);
    let h_e = 3.5 * params::E0 - 7.6875 * (params::E0 * params::E0 * params::E0);
    let predicted = 2.0 * PI / (n0 * (3.0 * params::B_MINUS_A_OVER_C * h_e).sqrt());
    let rel = (measured - predicted).abs() / predicted;
    assert!(
        rel < 0.05,
        "libration period: measured {} yr vs predicted {} yr (rel {rel})",
        measured / params::YEAR,
        predicted / params::YEAR
    );
}

/// 1.15 — the CVODE root function stops exactly on the requested ratio.
#[test]
fn root_function_stops_at_the_crossing() {
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, params::E0, 0.0, 0.0, 1.601 * n0],
    };
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let spec = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: false,
        t_end: 1.0e5 * params::YEAR,
        cadence: 100.0 * params::YEAR,
        root_ratio: 1.6,
        stop_on_root: true,
        reanchor: false,
        stage_tag: 'S',
        record: false,
    };
    let mut obs = noop();
    let end = integrate_segment(&ic, &spec, &mut obs, &mut samples, &mut stats)
        .expect("segment failed");
    let root = end.root_state.expect("root never fired");
    let ratio = root.y[4] / params::mean_motion(root.y[0]);
    assert!(
        (ratio - 1.6).abs() < 1.0e-9,
        "root ratio {ratio}, expected 1.6 to 1e-9"
    );
}
