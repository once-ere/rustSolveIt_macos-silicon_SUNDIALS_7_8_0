//! Integration tests for the integrators_adaptive module group of rebound_rs.
//! Part of rebound_rs, GPL-3.0-or-later.
#![allow(non_snake_case)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. See rebound_rust.md section 17.
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::manual_swap)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::misrefactored_assign_op)]
#![allow(clippy::neg_multiply)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_late_init)]
#![allow(clippy::while_let_loop)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::seek_from_current)]
#![allow(clippy::drop_non_drop)]
#![allow(clippy::redundant_closure)]
#![allow(clippy::approx_constant)]
#![allow(clippy::useless_vec)]
#![allow(clippy::type_complexity)]
use rebound_rs::*;

use rebound_rs::integrator_bs::{
    reb_integrator_bs_state, reb_integrator_bs_step_odes, reb_integrator_bs_update_particles,
    reb_ode, reb_ode_create,
};
use rebound_rs::integrator_ias15::{
    reb_integrator_ias15_timescale, REB_IAS15_ADAPTIVEMODE_AARSETH85,
    REB_IAS15_ADAPTIVEMODE_GLOBAL, REB_IAS15_ADAPTIVEMODE_INDIVIDUAL, REB_IAS15_ADAPTIVEMODE_PRS23,
};

const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// shared helpers
// ---------------------------------------------------------------------------

/// Star of mass 1 at the origin plus one body on the requested orbit.
/// `move_to_com` is applied so the barycentre is at rest.
fn two_body(m_planet: f64, a: f64, e: f64, inc: f64) -> reb_simulation {
    let mut r = reb_simulation_create();
    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(&mut r, star);
    let primary = r.particles[0];
    let p = reb_particle_from_orbit(r.G, primary, m_planet, a, e, inc, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, p);
    reb_simulation_move_to_com(&mut r);
    r
}

/// Star of mass `M` at the origin at rest plus a *massless* body on a
/// circular orbit of radius `a`. The star therefore feels no force and
/// the whole configuration is exactly the analytic Kepler problem.
fn massless_circular(M: f64, a: f64) -> reb_simulation {
    let mut r = reb_simulation_create();
    let mut star = reb_particle::default();
    star.m = M;
    reb_simulation_add(&mut r, star);
    let primary = r.particles[0];
    let p = reb_particle_from_orbit(r.G, primary, 0.0, a, 0.0, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, p);
    r
}

fn ias15_set<F: FnOnce(&mut rebound_rs::integrator_ias15::reb_integrator_ias15_state)>(
    r: &mut reb_simulation,
    f: F,
) {
    if let reb_integrator_state::ias15(ref mut i15) = r.integrator {
        f(i15);
    } else {
        panic!("simulation is not using the IAS15 integrator");
    }
}

fn ias15_iterations_max_exceeded(r: &reb_simulation) -> u64 {
    match r.integrator {
        reb_integrator_state::ias15(ref i15) => i15.iterations_max_exceeded,
        _ => panic!("simulation is not using the IAS15 integrator"),
    }
}

/// Every f64 that defines the dynamical state, as raw bits.
fn state_bits(r: &reb_simulation) -> Vec<u64> {
    let mut v = Vec::new();
    for p in &r.particles {
        for f in [
            p.x, p.y, p.z, p.vx, p.vy, p.vz, p.ax, p.ay, p.az, p.m, p.r,
        ] {
            v.push(f.to_bits());
        }
    }
    v.push(r.t.to_bits());
    v.push(r.dt.to_bits());
    v
}

fn rel_energy_drift(r: &reb_simulation, e0: f64) -> f64 {
    ((reb_simulation_energy(r) - e0) / e0).abs()
}

// ---------------------------------------------------------------------------
// IAS15 — Gauss-Radau spacings / quadrature exactness
// ---------------------------------------------------------------------------

/// Coefficients of a degree-7 polynomial acceleration; all are exact
/// binary fractions so the reference series below is unambiguous.
const POLY_A: [f64; 8] = [
    0.5,
    -0.25,
    0.125,
    -0.0625,
    0.03125,
    -0.015625,
    0.0078125,
    -0.00390625,
];

fn poly_force_deg7(r: &mut reb_simulation) {
    let t = r.t;
    // Horner from the top coefficient down.
    let mut a = 0.0;
    let mut k = POLY_A.len();
    while k > 0 {
        k -= 1;
        a = a * t + POLY_A[k];
    }
    r.particles[0].ax = a;
}

fn force_t13(r: &mut reb_simulation) {
    let t = r.t;
    r.particles[0].ax = t.powi(13);
}

fn force_t15(r: &mut reb_simulation) {
    let t = r.t;
    r.particles[0].ax = t.powi(15);
}

/// One particle at rest at the origin, no gravity, driven only by a
/// time-dependent `additional_forces` callback, fixed timestep.
fn driven_particle(force: fn(&mut reb_simulation), dt: f64) -> reb_simulation {
    let mut r = reb_simulation_create();
    r.gravity = REB_GRAVITY::NONE;
    r.additional_forces = Some(force);
    reb_simulation_add(&mut r, reb_particle::default());
    ias15_set(&mut r, |i15| i15.epsilon = 0.0);
    r.dt = dt;
    r
}

#[test]
fn ias15_integrates_degree7_polynomial_acceleration_exactly() {
    // a(t) = sum_k A_k t^k with deg 7. IAS15 fits the acceleration with
    // a0 + b0 h + ... + b6 h^7 through the eight Gauss-Radau nodes, so a
    // degree-7 acceleration is represented exactly and the closed-form
    // double integral must come back.
    let mut r = driven_particle(poly_force_deg7, 0.25);
    reb_simulation_steps(&mut r, 4);

    // 0.25 and 4 are exact, so the naive time accumulation is exact.
    assert_eq!(r.t.to_bits(), 1.0f64.to_bits(), "t after 4 fixed steps of 0.25");

    // v(t) = sum A_k t^{k+1}/(k+1); x(t) = sum A_k t^{k+2}/((k+1)(k+2))
    let t = 1.0f64;
    let mut v_ref = 0.0;
    let mut x_ref = 0.0;
    for k in 0..POLY_A.len() {
        let kf = k as f64;
        v_ref += POLY_A[k] * t.powi(k as i32 + 1) / (kf + 1.0);
        x_ref += POLY_A[k] * t.powi(k as i32 + 2) / ((kf + 1.0) * (kf + 2.0));
    }

    let x = r.particles[0].x;
    let v = r.particles[0].vx;
    assert!(
        (x - x_ref).abs() < 1e-14,
        "IAS15 position for degree-7 acceleration: got {:.17e}, exact {:.17e}, diff {:.3e}",
        x,
        x_ref,
        x - x_ref
    );
    assert!(
        (v - v_ref).abs() < 1e-14,
        "IAS15 velocity for degree-7 acceleration: got {:.17e}, exact {:.17e}, diff {:.3e}",
        v,
        v_ref,
        v - v_ref
    );
    // No transverse motion may appear out of nothing.
    assert_eq!(r.particles[0].y.to_bits(), 0.0f64.to_bits(), "y stays 0");
    assert_eq!(r.particles[0].z.to_bits(), 0.0f64.to_bits(), "z stays 0");
}

#[test]
fn ias15_radau_quadrature_is_exact_to_degree_13() {
    // The eight Gauss-Radau nodes (h[0] = 0 fixed) give a quadrature that
    // is exact for degree 2*8-2 = 14. The end-of-step position uses
    // int_0^1 (1-u) a(u) du, so a(u) = u^13 makes the integrand
    // (1-u)(u^13 - p7(u)) of degree 14 and it vanishes at every node:
    // the result must be exact. Reference: int_0^1 (1-u) u^13 du
    // = 1/14 - 1/15 = 1/210, and int_0^1 u^13 du = 1/14.
    let mut r = driven_particle(force_t13, 1.0);
    reb_simulation_steps(&mut r, 1);

    let x_ref = 1.0 / 14.0 - 1.0 / 15.0;
    let v_ref = 1.0 / 14.0;
    let x = r.particles[0].x;
    let v = r.particles[0].vx;
    // Only round-off separates the two: the g/b recursion divides by the
    // small Gauss-Radau differences rr[], which costs a couple of digits.
    assert!(
        ((x - x_ref) / x_ref).abs() < 1e-12,
        "Radau position quadrature of u^13: got {:.17e}, exact {:.17e}, rel diff {:.3e}",
        x,
        x_ref,
        (x - x_ref) / x_ref
    );
    assert!(
        ((v - v_ref) / v_ref).abs() < 1e-12,
        "Radau velocity quadrature of u^13: got {:.17e}, exact {:.17e}, rel diff {:.3e}",
        v,
        v_ref,
        (v - v_ref) / v_ref
    );
}

#[test]
fn ias15_radau_quadrature_saturates_above_degree_14() {
    // Same construction with a(u) = u^15: now (1-u)(u^15 - p7(u)) has
    // degree 16 > 14, so the Radau quadrature is no longer exact and the
    // position must miss the closed form by far more than round-off.
    // int_0^1 (1-u) u^15 du = 1/16 - 1/17 = 1/272.
    let mut r = driven_particle(force_t15, 1.0);
    reb_simulation_steps(&mut r, 1);

    let x_ref = 1.0 / 16.0 - 1.0 / 17.0;
    let err = (r.particles[0].x - x_ref).abs();
    // The degree-13 companion test above lands within 3e-16 of its own
    // closed form; here the truncation error is many orders larger.
    assert!(
        err > 1e-10,
        "degree-15 acceleration must NOT be integrated exactly, but the error was only {:.3e}",
        err
    );
    // ... yet the method is still 15th order, so the miss stays small.
    assert!(
        err < 1e-5,
        "degree-15 quadrature error unexpectedly large: {:.3e}",
        err
    );
}

// ---------------------------------------------------------------------------
// IAS15 — compensated summation
// ---------------------------------------------------------------------------

/// A vanishingly small but *normal* acceleration. It keeps the
/// predictor-corrector convergence test (|b6|/|a|) well defined while
/// being far below the round-off of the position accumulator, so the
/// motion is uniform to the last bit.
fn tiny_force(r: &mut reb_simulation) {
    r.particles[0].ax = 1e-300;
}

#[test]
fn ias15_compensated_summation_beats_naive_accumulation() {
    // Uniform motion, fixed step: every step adds exactly fl(0.1) to x.
    // add_cs() is Kahan compensated summation, so after n steps x must
    // agree with the correctly rounded n*fl(0.1) -- which a single
    // correctly rounded multiply reproduces -- to within an ulp, while a
    // naive running sum of the same increments drifts away.
    const N_STEPS: usize = 100_000;
    let dt = 0.1f64;

    let mut r = driven_particle(tiny_force, dt);
    r.particles[0].vx = 1.0;
    reb_simulation_steps(&mut r, N_STEPS);

    let exact = (N_STEPS as f64) * dt; // correctly rounded sum of n copies of fl(0.1)
    let x = r.particles[0].x;

    let mut naive = 0.0f64;
    for _ in 0..N_STEPS {
        naive += dt;
    }

    let ulp = (exact.abs() * f64::EPSILON).max(f64::MIN_POSITIVE);
    let err_cs = (x - exact).abs();
    let err_naive = (naive - exact).abs();

    assert!(
        err_cs <= 2.0 * ulp,
        "IAS15 compensated position after {} steps: got {:.17e}, exact {:.17e}, err {:.3e} (> 2 ulp = {:.3e})",
        N_STEPS,
        x,
        exact,
        err_cs,
        2.0 * ulp
    );
    assert!(
        err_naive > 10.0 * err_cs.max(ulp),
        "the test is not discriminating: naive summation err {:.3e} is not clearly worse than compensated err {:.3e}",
        err_naive,
        err_cs
    );
    // The velocity is untouched by the 1e-300 kick at double precision.
    assert_eq!(
        r.particles[0].vx.to_bits(),
        1.0f64.to_bits(),
        "vx must stay exactly 1.0"
    );
}

// ---------------------------------------------------------------------------
// IAS15 — PRS23 timescale
// ---------------------------------------------------------------------------

#[test]
fn ias15_timescale_prs23_is_exact_for_circular_orbits() {
    // For a massless body on a circular orbit of radius R about a mass M
    // the pairwise sums collapse to |y2| = (n^2 R)^2, |y3| = (n^3 R)^2,
    // |y4| = (n^4 R)^2 with n = sqrt(GM/R^3), hence
    //   timescale^2 = 2 y2 / (y3 + sqrt(y4 y2)) = 1/n^2.
    // Chosen so every intermediate is a power of two: the answer is bit
    // exact.
    let cases: [(f64, f64, f64); 3] = [(1.0, 1.0, 1.0), (1.0, 4.0, 8.0), (4.0, 1.0, 0.5)];
    for &(M, a, expect) in &cases {
        let mut r = massless_circular(M, a);
        let ts = reb_integrator_ias15_timescale(&mut r);
        // n = sqrt(G M / a^3), timescale = 1/n = sqrt(a^3/(G M))
        let n = (r.G * M / (a * a * a)).sqrt();
        assert_eq!(
            expect.to_bits(),
            (1.0 / n).to_bits(),
            "reference 1/n for M={} a={} is not the power of two claimed",
            M,
            a
        );
        assert_eq!(
            ts.to_bits(),
            expect.to_bits(),
            "reb_integrator_ias15_timescale for M={} a={}: got {:.17e}, expected {:.17e}",
            M,
            a,
            ts,
            expect
        );
    }
}

#[test]
fn ias15_timescale_is_infinite_without_interactions() {
    // N = 1: no pair, the single particle has no acceleration, so the
    // "skip particles with non-normal y2" branch skips everything and
    // min_timescale2 stays INFINITY.
    let mut r = reb_simulation_create();
    reb_simulation_add(&mut r, reb_particle::default());
    let ts1 = reb_integrator_ias15_timescale(&mut r);
    assert!(
        ts1.is_infinite() && ts1 > 0.0,
        "timescale for N=1 should be +inf, got {}",
        ts1
    );

    // N = 0: the loop body never runs.
    let mut r0 = reb_simulation_create();
    let ts0 = reb_integrator_ias15_timescale(&mut r0);
    assert!(
        ts0.is_infinite() && ts0 > 0.0,
        "timescale for N=0 should be +inf, got {}",
        ts0
    );
}

// ---------------------------------------------------------------------------
// IAS15 — step size control
// ---------------------------------------------------------------------------

#[test]
fn ias15_equilibrium_timestep_matches_prs23_formula() {
    // In PRS23 mode the accepted step obeys
    //   dt_new = sqrt(min_timescale2) * dt_done * (5040 eps)^(1/7)
    // and sqrt(min_timescale2)*dt_done is the physical timescale, which
    // for the circular orbit above is exactly 1/n. So the timestep must
    // relax onto T * (5040 eps)^(1/7) independently of where it started.
    for &eps in &[1e-9f64, 1e-12] {
        let mut r = massless_circular(1.0, 1.0); // n = 1, timescale = 1
        ias15_set(&mut r, |i15| {
            i15.epsilon = eps;
            i15.adaptive_mode = REB_IAS15_ADAPTIVEMODE_PRS23;
        });
        r.dt = 0.01;
        reb_simulation_steps(&mut r, 200);

        let predicted = (5040.0 * eps).powf(1.0 / 7.0);
        let got = r.dt;
        assert!(
            (got - predicted).abs() < 0.02 * predicted,
            "equilibrium IAS15 timestep for eps={:e}: got {:.6e}, PRS23 formula gives {:.6e}",
            eps,
            got,
            predicted
        );
    }

    // Halving-in-log: eps down by 1e3 shrinks dt by 1000^(1/7).
    let mut ra = massless_circular(1.0, 1.0);
    ias15_set(&mut ra, |i15| i15.epsilon = 1e-9);
    ra.dt = 0.01;
    reb_simulation_steps(&mut ra, 200);
    let mut rb = massless_circular(1.0, 1.0);
    ias15_set(&mut rb, |i15| i15.epsilon = 1e-12);
    rb.dt = 0.01;
    reb_simulation_steps(&mut rb, 200);
    let ratio = ra.dt / rb.dt;
    let ratio_ref = 1000.0f64.powf(1.0 / 7.0);
    assert!(
        (ratio - ratio_ref).abs() < 0.05 * ratio_ref,
        "dt(eps=1e-9)/dt(eps=1e-12) = {:.5}, expected 1000^(1/7) = {:.5}",
        ratio,
        ratio_ref
    );
}

#[test]
fn ias15_min_dt_clamps_the_timestep() {
    // The natural step for this problem is ~0.175 (see the test above).
    // With min_dt = 0.5 the controller must clamp to exactly 0.5 and,
    // since dt_new/dt_done == 1, neither reject nor grow the step.
    let mut r = massless_circular(1.0, 1.0);
    ias15_set(&mut r, |i15| {
        i15.epsilon = 1e-9;
        i15.min_dt = 0.5;
    });
    r.dt = 0.5;
    reb_simulation_steps(&mut r, 5);
    assert_eq!(
        r.dt.to_bits(),
        0.5f64.to_bits(),
        "min_dt must pin dt to 0.5, got {:.17e}",
        r.dt
    );
    assert_eq!(
        r.dt_last_done.to_bits(),
        0.5f64.to_bits(),
        "the completed step must also be 0.5"
    );
    assert_eq!(r.t.to_bits(), 2.5f64.to_bits(), "5 steps of 0.5 gives t=2.5");
}

#[test]
fn ias15_epsilon_zero_gives_a_fixed_timestep() {
    // epsilon <= 0 skips the whole step size controller.
    let mut r = two_body(1e-3, 1.0, 0.2, 0.0);
    ias15_set(&mut r, |i15| i15.epsilon = 0.0);
    r.dt = 0.125; // exact binary
    reb_simulation_steps(&mut r, 32);
    assert_eq!(
        r.dt.to_bits(),
        0.125f64.to_bits(),
        "dt must not change when epsilon == 0"
    );
    assert_eq!(
        r.t.to_bits(),
        4.0f64.to_bits(),
        "t must be exactly 32*0.125 = 4"
    );
}

#[test]
fn ias15_timestep_shrinks_towards_pericentre() {
    // An e = 0.9 orbit needs a far smaller step at pericentre than at
    // apocentre. The body starts at pericentre (f = 0); sample dt over a
    // whole orbit and compare the extremes with the ratio of the
    // Keplerian timescale 1/n_local ~ r^{3/2}: (1+e)/(1-e) = 19 in
    // radius, so 19^{3/2} ~ 82 in timescale.
    let mut r = two_body(0.0, 1.0, 0.9, 0.0);
    ias15_set(&mut r, |i15| i15.epsilon = 1e-9);
    r.dt = 1e-3;
    let period = 2.0 * PI;

    let mut dt_min = f64::INFINITY;
    let mut dt_max = 0.0f64;
    // one full orbit, sampling the accepted step
    while r.t < period {
        reb_simulation_steps(&mut r, 1);
        let d = r.dt_last_done.abs();
        if d < dt_min {
            dt_min = d;
        }
        if d > dt_max {
            dt_max = d;
        }
    }
    assert!(
        dt_max / dt_min > 10.0,
        "adaptive step should vary strongly over an e=0.9 orbit: dt_max/dt_min = {:.3}",
        dt_max / dt_min
    );
    assert!(
        dt_max < period,
        "no accepted step may span a whole orbit: dt_max = {:.4}, P = {:.4}",
        dt_max,
        period
    );
}

#[test]
fn ias15_recovers_from_an_absurdly_large_initial_timestep() {
    // Start ten orbits per step. The controller must reject repeatedly
    // and settle near the PRS23 equilibrium 1*(5040*eps)^(1/7).
    let mut r = massless_circular(1.0, 1.0);
    r.save_messages = 1; // capture the non-convergence warnings
    ias15_set(&mut r, |i15| i15.epsilon = 1e-9);
    r.dt = 10.0 * 2.0 * PI;
    reb_simulation_steps(&mut r, 40);

    let predicted = (5040.0 * 1e-9f64).powf(1.0 / 7.0);
    assert!(
        (r.dt - predicted).abs() < 0.05 * predicted,
        "after recovering from a huge dt the step should be ~{:.5e}, got {:.5e}",
        predicted,
        r.dt
    );
    // The completed step must be finite and forward.
    assert!(
        r.dt_last_done > 0.0 && r.dt_last_done.is_finite(),
        "dt_last_done = {}",
        r.dt_last_done
    );
}

#[test]
fn ias15_predictor_corrector_converges_on_a_smooth_orbit() {
    // A well resolved circular orbit must never hit the 12-iteration cap.
    let eps = 1e-9f64;
    let n_orbits = 100.0;
    let period = 2.0 * PI; // a = 1, G = 1, M = 1 => n = 1
    let mut r = massless_circular(1.0, 1.0);
    ias15_set(&mut r, |i15| i15.epsilon = eps);
    reb_simulation_integrate(&mut r, n_orbits * period);
    assert_eq!(
        ias15_iterations_max_exceeded(&r),
        0,
        "IAS15 predictor-corrector should converge for every step of a circular orbit"
    );

    // Because the timescale of a circular orbit is constant, every step
    // is accepted at the PRS23 equilibrium size T*(5040 eps)^(1/7) = 1 *
    // (5040 eps)^(1/7); the step count therefore follows from the
    // integration span alone. steps_done also counts rejected steps, so
    // any thrashing in the controller would show up as a surplus.
    let dt_eq = (5040.0 * eps).powf(1.0 / 7.0);
    let steps_ref = n_orbits * period / dt_eq;
    let got = r.steps_done as f64;
    assert!(
        (got - steps_ref).abs() < 0.03 * steps_ref,
        "IAS15 took {} steps over {} circular orbits; the equilibrium step {:.6e} predicts {:.1}",
        r.steps_done,
        n_orbits,
        dt_eq,
        steps_ref
    );
}

// ---------------------------------------------------------------------------
// IAS15 — physical invariants
// ---------------------------------------------------------------------------

#[test]
fn ias15_conserves_energy_and_angular_momentum() {
    let mut r = two_body(1e-3, 1.0, 0.0, 0.0);
    let e0 = reb_simulation_energy(&r);
    let l0 = reb_simulation_angular_momentum(&r);
    reb_simulation_integrate(&mut r, 200.0 * 2.0 * PI);
    let drift = rel_energy_drift(&r, e0);
    assert!(
        drift < 1e-12,
        "IAS15 relative energy drift over 200 circular orbits: {:.3e}",
        drift
    );
    let l1 = reb_simulation_angular_momentum(&r);
    let ldrift = ((l1.z - l0.z) / l0.z).abs();
    assert!(
        ldrift < 1e-13,
        "IAS15 relative Lz drift over 200 circular orbits: {:.3e}",
        ldrift
    );
    assert!(
        l1.x.abs() < 1e-14 && l1.y.abs() < 1e-14,
        "a coplanar orbit must keep Lx = Ly = 0, got ({:.3e}, {:.3e})",
        l1.x,
        l1.y
    );
}

#[test]
fn ias15_conserves_energy_at_high_eccentricity() {
    let mut r = two_body(1e-3, 1.0, 0.99, 0.0);
    let e0 = reb_simulation_energy(&r);
    reb_simulation_integrate(&mut r, 50.0 * 2.0 * PI);
    let drift = rel_energy_drift(&r, e0);
    assert!(
        drift < 1e-10,
        "IAS15 relative energy drift over 50 orbits at e=0.99: {:.3e}",
        drift
    );
}

#[test]
fn ias15_conserves_the_orbital_elements_of_a_hyperbolic_flyby() {
    // a < 0, e > 1. The two-body elements are integrals of the motion,
    // so recomputing them from the final state must reproduce the input.
    let a0 = -2.0;
    let e0 = 1.5;
    let mut r = two_body(0.0, a0, e0, 0.0);
    // Start well before pericentre so the flyby is actually integrated.
    let start = reb_particle_from_orbit(r.G, r.particles[0], 0.0, a0, e0, 0.0, 0.0, 0.0, -1.0);
    r.particles[1] = start;

    let orb_in = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);
    assert!(orb_in.e > 1.0, "setup must be hyperbolic, e = {}", orb_in.e);
    let d_in = reb_particle_distance(&r.particles[0], &r.particles[1]);
    // f = -1 rad is before pericentre, so the body is falling inwards.
    let rv_in = {
        let p = r.particles[1];
        let s = r.particles[0];
        (p.x - s.x) * (p.vx - s.vx) + (p.y - s.y) * (p.vy - s.vy) + (p.z - s.z) * (p.vz - s.vz)
    };
    assert!(rv_in < 0.0, "the flyby must start inbound, r.v = {:.6}", rv_in);

    // Track the closest approach; it must match the analytic pericentre
    // distance q = a(1-e) = -2*(1-1.5) = 1.
    let q = a0 * (1.0 - e0);
    let mut d_min = f64::INFINITY;
    while r.t < 20.0 {
        reb_simulation_steps(&mut r, 1);
        let d = reb_particle_distance(&r.particles[0], &r.particles[1]);
        if d < d_min {
            d_min = d;
        }
    }
    let orb_out = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);

    assert!(
        ((orb_out.a - orb_in.a) / orb_in.a).abs() < 1e-12,
        "hyperbolic a not conserved: {:.17e} -> {:.17e}",
        orb_in.a,
        orb_out.a
    );
    assert!(
        ((orb_out.e - orb_in.e) / orb_in.e).abs() < 1e-12,
        "hyperbolic e not conserved: {:.17e} -> {:.17e}",
        orb_in.e,
        orb_out.e
    );
    // The sampled minimum can only ever sit at or above the true
    // pericentre, and with the adaptive step it must sit close to it.
    assert!(
        d_min >= q - 1e-12,
        "closest approach {:.12} fell below the analytic pericentre {:.12}",
        d_min,
        q
    );
    assert!(
        d_min - q < 1e-2,
        "closest approach {:.12} missed the analytic pericentre {:.12} by {:.3e}",
        d_min,
        q,
        d_min - q
    );
    // ... and the body must leave again.
    let rv_out = {
        let p = r.particles[1];
        let s = r.particles[0];
        (p.x - s.x) * (p.vx - s.vx) + (p.y - s.y) * (p.vy - s.vy) + (p.z - s.z) * (p.vz - s.vz)
    };
    assert!(
        rv_out > 0.0,
        "the flyby must end outbound, r.v = {:.6}",
        rv_out
    );
    assert!(
        orb_out.d > d_in,
        "an unbound orbit must end further out: {:.6} -> {:.6}",
        d_in,
        orb_out.d
    );
}

#[test]
fn ias15_preserves_a_retrograde_inclination() {
    // inc = pi is retrograde and coplanar: Lz must stay negative and the
    // recovered inclination must stay pi.
    let inc = PI;
    let mut r = two_body(1e-4, 1.0, 0.1, inc);
    let l0 = reb_simulation_angular_momentum(&r);
    assert!(l0.z < 0.0, "retrograde orbit must have Lz < 0, got {}", l0.z);
    reb_simulation_integrate(&mut r, 30.0 * 2.0 * PI);
    let orb = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);
    assert!(
        (orb.inc - inc).abs() < 1e-10,
        "retrograde inclination drifted: {:.17e} vs {:.17e}",
        orb.inc,
        inc
    );
    let l1 = reb_simulation_angular_momentum(&r);
    assert!(
        ((l1.z - l0.z) / l0.z).abs() < 1e-12,
        "retrograde Lz drift: {:.3e}",
        ((l1.z - l0.z) / l0.z).abs()
    );
}

#[test]
fn ias15_zero_and_one_particle_edge_cases() {
    // N = 0: the integrator just advances time.
    let mut r = reb_simulation_create();
    r.dt = 0.25;
    reb_simulation_steps(&mut r, 4);
    assert_eq!(r.t.to_bits(), 1.0f64.to_bits(), "N=0 must advance t by dt");
    assert_eq!(
        r.dt_last_done.to_bits(),
        0.25f64.to_bits(),
        "N=0 must record dt_last_done"
    );

    // N = 1 with no forces: exact straight line motion.
    let mut r1 = reb_simulation_create();
    r1.gravity = REB_GRAVITY::NONE;
    let mut p = reb_particle::default();
    p.x = 0.5;
    p.vx = 2.0;
    p.vy = -0.25;
    reb_simulation_add(&mut r1, p);
    ias15_set(&mut r1, |i15| i15.epsilon = 0.0);
    r1.dt = 0.0625;
    r1.save_messages = 1; // a force-free particle exhausts the convergence test
    reb_simulation_steps(&mut r1, 16);
    assert_eq!(r1.t.to_bits(), 1.0f64.to_bits(), "N=1 t after 16*0.0625");
    assert!(
        (r1.particles[0].x - 2.5).abs() < 1e-15,
        "N=1 x should be 0.5 + 2*1 = 2.5, got {:.17e}",
        r1.particles[0].x
    );
    assert!(
        (r1.particles[0].y + 0.25).abs() < 1e-15,
        "N=1 y should be -0.25, got {:.17e}",
        r1.particles[0].y
    );
    assert_eq!(
        r1.particles[0].vx.to_bits(),
        2.0f64.to_bits(),
        "force free vx unchanged"
    );
}

#[test]
fn ias15_integrates_backwards_in_time() {
    let mut r = two_body(1e-3, 1.0, 0.3, 0.0);
    let e0 = reb_simulation_energy(&r);
    r.dt = -0.01;
    reb_simulation_integrate(&mut r, -20.0);
    assert!(
        (r.t + 20.0).abs() < 1e-12,
        "backward integration must land on t = -20, got {:.17e}",
        r.t
    );
    assert!(r.dt < 0.0, "dt must stay negative, got {}", r.dt);
    let drift = rel_energy_drift(&r, e0);
    assert!(
        drift < 1e-12,
        "backward IAS15 relative energy drift: {:.3e}",
        drift
    );
}

#[test]
fn ias15_fixed_step_is_nearly_time_reversible() {
    // With epsilon = 0 the same sequence of steps run backwards must
    // return to the starting state to well below the truncation error.
    let mut r = two_body(1e-3, 1.0, 0.2, 0.0);
    ias15_set(&mut r, |i15| i15.epsilon = 0.0);
    r.dt = 0.03125;
    let start: Vec<reb_particle> = r.particles.clone();

    reb_simulation_steps(&mut r, 320); // t = 10
    r.dt = -0.03125;
    reb_simulation_steps(&mut r, 320); // back to t = 0

    assert!(r.t.abs() < 1e-13, "t should return to 0, got {:.3e}", r.t);
    for i in 0..r.N {
        let dx = r.particles[i].x - start[i].x;
        let dy = r.particles[i].y - start[i].y;
        let dvx = r.particles[i].vx - start[i].vx;
        let dvy = r.particles[i].vy - start[i].vy;
        assert!(
            dx.abs() < 1e-11 && dy.abs() < 1e-11,
            "particle {} position not recovered after forward+backward: ({:.3e}, {:.3e})",
            i,
            dx,
            dy
        );
        assert!(
            dvx.abs() < 1e-11 && dvy.abs() < 1e-11,
            "particle {} velocity not recovered after forward+backward: ({:.3e}, {:.3e})",
            i,
            dvx,
            dvy
        );
    }
}

#[test]
fn ias15_all_adaptive_modes_conserve_energy() {
    for &mode in &[
        REB_IAS15_ADAPTIVEMODE_INDIVIDUAL,
        REB_IAS15_ADAPTIVEMODE_GLOBAL,
        REB_IAS15_ADAPTIVEMODE_PRS23,
        REB_IAS15_ADAPTIVEMODE_AARSETH85,
    ] {
        let mut r = two_body(1e-3, 1.0, 0.3, 0.0);
        ias15_set(&mut r, |i15| {
            i15.epsilon = 1e-9;
            i15.adaptive_mode = mode;
        });
        let e0 = reb_simulation_energy(&r);
        reb_simulation_integrate(&mut r, 50.0 * 2.0 * PI);
        let drift = rel_energy_drift(&r, e0);
        assert!(
            drift < 1e-11,
            "adaptive_mode {} relative energy drift over 50 orbits: {:.3e}",
            mode,
            drift
        );
        assert!(
            r.dt > 0.0 && r.dt < 2.0 * PI,
            "adaptive_mode {} produced an implausible step {:.3e}",
            mode,
            r.dt
        );
    }
}

#[test]
fn ias15_is_bit_deterministic() {
    let run = || {
        let mut r = two_body(1e-3, 1.0, 0.4, 0.3);
        ias15_set(&mut r, |i15| i15.epsilon = 1e-9);
        reb_simulation_integrate(&mut r, 37.0);
        r
    };
    let a = run();
    let b = run();
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical IAS15 runs must be bit identical"
    );
}

#[test]
fn ias15_agrees_with_whfast_on_a_two_body_problem() {
    // WHFast is exact for the Kepler part; its error is second order in
    // dt and first order in the planet mass ratio. With dt = P/2000 and
    // m/M = 1e-3 the two integrators must agree far better than 1e-6.
    let tmax = 3.0 * 2.0 * PI;

    let mut a = two_body(1e-3, 1.0, 0.15, 0.0);
    reb_simulation_integrate(&mut a, tmax);

    let mut b = two_body(1e-3, 1.0, 0.15, 0.0);
    reb_simulation_set_integrator(&mut b, "whfast");
    b.dt = 2.0 * PI / 2000.0;
    reb_simulation_integrate(&mut b, tmax);

    let dx = a.particles[1].x - b.particles[1].x;
    let dy = a.particles[1].y - b.particles[1].y;
    let d = (dx * dx + dy * dy).sqrt();
    assert!(
        d < 1e-6,
        "IAS15 vs WHFast separation after 3 orbits: {:.3e} (positions {:?} / {:?})",
        d,
        (a.particles[1].x, a.particles[1].y),
        (b.particles[1].x, b.particles[1].y)
    );

    // Semi-major axes must agree much more tightly than the phase.
    let oa = reb_orbit_from_particle(a.G, a.particles[1], a.particles[0]);
    let ob = reb_orbit_from_particle(b.G, b.particles[1], b.particles[0]);
    assert!(
        ((oa.a - ob.a) / oa.a).abs() < 1e-9,
        "IAS15 vs WHFast semi-major axis: {:.17e} vs {:.17e}",
        oa.a,
        ob.a
    );
}

// ---------------------------------------------------------------------------
// BS — ODE plumbing
// ---------------------------------------------------------------------------

fn ho_derivatives(
    _r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    y: &[f64],
    _t: f64,
) {
    yDot[0] = y[1];
    yDot[1] = -y[0];
}

fn decay_derivatives(
    _r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    y: &[f64],
    _t: f64,
) {
    yDot[0] = -y[0];
}

fn cos_derivatives(
    _r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    _y: &[f64],
    t: f64,
) {
    yDot[0] = t.cos();
}

fn const_derivatives(
    _r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    _y: &[f64],
    _t: f64,
) {
    yDot[0] = 0.75;
}

/// A simulation carrying a single standalone ODE (no N-body coupling),
/// plus a freshly created BS state to drive it with.
fn ode_sim(
    length: usize,
    f: rebound_rs::integrator_bs::reb_ode_derivatives_fn,
    y0: &[f64],
) -> (reb_simulation, reb_integrator_bs_state, usize) {
    let mut r = reb_simulation_create();
    reb_simulation_set_integrator(&mut r, "bs");
    let id = reb_ode_create(&mut r, length);
    {
        let ode = r.odes.iter_mut().find(|o| o.id == id).unwrap();
        ode.derivatives = Some(f);
        ode.needs_nbody = 0;
        ode.y[..length].copy_from_slice(&y0[..length]);
    }
    let bs = reb_integrator_bs_state::default();
    (r, bs, id)
}

fn ode_y(r: &reb_simulation, id: usize) -> Vec<f64> {
    r.odes.iter().find(|o| o.id == id).unwrap().y.clone()
}

/// Advance a standalone ODE to `t_end` with BS' own step size control.
fn bs_run(
    r: &mut reb_simulation,
    bs: &mut reb_integrator_bs_state,
    _id: usize,
    t_end: f64,
    dt0: f64,
) -> usize {
    let mut dt = dt0;
    let mut steps = 0usize;
    while r.t < t_end && steps < 200_000 {
        steps += 1;
        let mut d = dt;
        if r.t + d > t_end {
            d = t_end - r.t;
            bs.first_or_last_step = 1;
        }
        if reb_integrator_bs_step_odes(r, bs, d) != 0 {
            r.t += d;
        }
        dt = bs.dt_proposed;
        assert!(dt > 0.0 && dt.is_finite(), "BS proposed a bad dt: {}", dt);
    }
    steps
}

/// The modified midpoint rule exactly as integrator_bs.c `tryStep`
/// writes it, for the harmonic oscillator y' = (y1, -y0).
fn ho_modified_midpoint(y0: [f64; 2], step: f64, n: i32) -> [f64; 2] {
    let sub = step / (n as f64);
    let y0dot = [y0[1], -y0[0]];
    let mut y1 = [y0[0] + sub * y0dot[0], y0[1] + sub * y0dot[1]];
    let mut ydot = [y1[1], -y1[0]];
    let mut ytmp = y0;
    for _j in 1..n {
        let middle = y1;
        y1 = [
            ytmp[0] + 2. * sub * ydot[0],
            ytmp[1] + 2. * sub * ydot[1],
        ];
        ytmp = middle;
        ydot = [y1[1], -y1[0]];
    }
    [
        0.5 * (ytmp[0] + y1[0] + sub * ydot[0]),
        0.5 * (ytmp[1] + y1[1] + sub * ydot[1]),
    ]
}

#[test]
fn bs_sequence_and_cost_tables_are_as_documented() {
    let bs = reb_integrator_bs_state::default();
    // step size sequence 2, 6, 10, ... = 4k+2
    let seq_ref = [2i32, 6, 10, 14, 18, 22, 26, 30, 34];
    assert_eq!(bs.sequence, seq_ref, "BS substep sequence");
    // cost[0] = seq[0]+1, cost[k] = cost[k-1] + seq[k]
    let cost_ref = [3i32, 9, 19, 33, 51, 73, 99, 129, 163];
    assert_eq!(bs.cost_per_step, cost_ref, "BS cost per step");
    // coeff[k] = (1/seq[k])^2 -- the squared substep size, which is what
    // the Richardson extrapolation eliminates.
    assert_eq!(
        bs.coeff[0].to_bits(),
        0.25f64.to_bits(),
        "coeff[0] must be exactly (1/2)^2"
    );
    for k in 0..bs.sequence.len() {
        let n = bs.sequence[k] as f64;
        assert!(
            (bs.coeff[k] * n * n - 1.0).abs() < 1e-15,
            "coeff[{}] = {:.17e} is not 1/{}^2",
            k,
            bs.coeff[k],
            bs.sequence[k]
        );
    }
    assert_eq!(bs.first_or_last_step, 1, "a fresh BS state starts flagged");
    assert_eq!(bs.target_iter, 0, "order is selected on the first step");
}

#[test]
fn bs_modified_midpoint_is_exact_for_a_constant_derivative() {
    // y' = 0.75 => y(t) = y0 + 0.75 t. Every substep of the modified
    // midpoint reproduces this and the extrapolation of identical values
    // is the identity, so a single step must land on the closed form.
    let (mut r, mut bs, id) = ode_sim(1, const_derivatives, &[0.0]);
    bs.eps_abs = 1e-10;
    bs.eps_rel = 1e-10;
    let dt = 0.5;
    assert_eq!(
        reb_integrator_bs_step_odes(&mut r, &mut bs, dt),
        1,
        "a constant derivative must never be rejected"
    );
    r.t += dt;
    let y = ode_y(&r, id)[0];
    let expect = 0.75 * 0.5;
    assert!(
        (y - expect).abs() < 1e-16,
        "BS on y'=0.75 after dt=0.5: got {:.17e}, exact {:.17e}",
        y,
        expect
    );
}

#[test]
fn bs_modified_midpoint_matches_the_hand_computed_two_substep_rule() {
    // Hand evaluation of tryStep for y=(1,0), step=1, n=2 (subStep=0.5):
    //   y1 = (1, -0.5); yDot = (-0.5, -1); yTmp = (1, 0)
    //   j=1: y1 = (1,0) + 1*(-0.5,-1) = (0.5, -1); yTmp = (1, -0.5)
    //        yDot = (-1, -0.5)
    //   corr: 0.5*((1,-0.5) + (0.5,-1) + 0.5*(-1,-0.5)) = (0.5, -0.875)
    // Every value here is an exact binary fraction.
    let m = ho_modified_midpoint([1.0, 0.0], 1.0, 2);
    assert_eq!(m[0].to_bits(), 0.5f64.to_bits(), "midpoint q, n=2, step=1");
    assert_eq!(
        m[1].to_bits(),
        (-0.875f64).to_bits(),
        "midpoint p, n=2, step=1"
    );
    // ... and it really is a poor approximation to (cos 1, -sin 1).
    let err = ((m[0] - 1.0f64.cos()).powi(2) + (m[1] + 1.0f64.sin()).powi(2)).sqrt();
    assert!(
        err > 1e-2,
        "the two-substep midpoint should be crude here, err = {:.3e}",
        err
    );
}

#[test]
fn bs_richardson_extrapolation_combines_the_two_lowest_substeps() {
    // With eps_rel = 1e-2 the initial order selection gives
    //   target_iter = clamp(floor(0.5 - 0.6*log10(1e-2)), 1, 7)
    //               = clamp(floor(1.7), 1, 7) = 1,
    // so a small easy step converges exactly at k = 1. There the state
    // is the Richardson combination of the n=2 and n=6 modified midpoint
    // results with the weights built from coeff[] :
    //   facD = coeff[1]/(coeff[0]-coeff[1]),  y = facD*(M6-M2) + M6
    // which is the classic (9*M6 - M2)/8 for h^2 error elimination.
    let y0 = [1.0f64, 0.0];
    let step = 0.05f64;
    let (mut r, mut bs, id) = ode_sim(2, ho_derivatives, &y0);
    bs.eps_abs = 1e-2;
    bs.eps_rel = 1e-2;
    assert_eq!(
        reb_integrator_bs_step_odes(&mut r, &mut bs, step),
        1,
        "the easy step must be accepted"
    );
    // Converging at k == 1 with no prior rejection sets optimalIter = 2
    // and, since previous_rejected == 0, target_iter = optimalIter. Any
    // other exit from the loop would leave a different value here.
    assert_eq!(
        bs.target_iter, 2,
        "the k == 1 exit must leave target_iter at 2"
    );
    assert_eq!(
        bs.optimal_step[2].to_bits(),
        0.0f64.to_bits(),
        "order 2 must never have been attempted"
    );

    let m2 = ho_modified_midpoint(y0, step, 2);
    let m6 = ho_modified_midpoint(y0, step, 6);
    let c0 = {
        let q = 1.0 / 2.0f64;
        q * q
    };
    let c1 = {
        let q = 1.0 / 6.0f64;
        q * q
    };
    let facd = c1 / (c0 - c1);
    // sanity: the analytic weights are 9/8 and -1/8
    assert!(
        (facd - 0.125).abs() < 1e-15,
        "Richardson weight facD should be 1/8, got {:.17e}",
        facd
    );

    let y = ode_y(&r, id);
    for i in 0..2 {
        let expect = facd * (m6[i] - m2[i]) + m6[i];
        assert!(
            (y[i] - expect).abs() <= 1e-15 * expect.abs().max(1e-3),
            "BS component {}: got {:.17e}, Richardson(M2,M6) = {:.17e}, diff {:.3e}",
            i,
            y[i],
            expect,
            y[i] - expect
        );
    }
    // The extrapolation must also be far better than either input.
    let exact = [step.cos(), -step.sin()];
    let e_bs = ((y[0] - exact[0]).powi(2) + (y[1] - exact[1]).powi(2)).sqrt();
    let e_m2 = ((m2[0] - exact[0]).powi(2) + (m2[1] - exact[1]).powi(2)).sqrt();
    assert!(
        e_bs < e_m2 / 100.0,
        "extrapolated error {:.3e} should be far below the n=2 midpoint error {:.3e}",
        e_bs,
        e_m2
    );
}

#[test]
fn bs_harmonic_oscillator_accuracy_follows_the_tolerance() {
    let y0 = [1.0f64, 0.0];
    let t_end = 2.0 * PI;
    let mut errs = Vec::new();
    for &tol in &[1e-6f64, 1e-9, 1e-12] {
        let (mut r, mut bs, id) = ode_sim(2, ho_derivatives, &y0);
        bs.eps_abs = tol;
        bs.eps_rel = tol;
        bs_run(&mut r, &mut bs, id, t_end, 0.05);
        let y = ode_y(&r, id);
        // exact solution at 2 pi is (1, 0)
        let err = ((y[0] - 1.0).powi(2) + y[1].powi(2)).sqrt();
        assert!(
            err < 1e4 * tol,
            "BS harmonic oscillator error {:.3e} at tol {:.0e} exceeds 1e4*tol",
            err,
            tol
        );
        errs.push(err);
    }
    assert!(
        errs[0] > errs[1] && errs[1] > errs[2],
        "tightening the BS tolerance must reduce the error: {:?}",
        errs
    );
    assert!(
        errs[2] < 1e-9,
        "BS at tol 1e-12 should reach at least 1e-9 over one period, got {:.3e}",
        errs[2]
    );
}

#[test]
fn bs_integrates_a_time_dependent_right_hand_side() {
    // y' = cos t, y(0) = 0 => y(t) = sin t. This exercises the t
    // bookkeeping inside the modified midpoint (t += subStep).
    let (mut r, mut bs, id) = ode_sim(1, cos_derivatives, &[0.0]);
    bs.eps_abs = 1e-12;
    bs.eps_rel = 1e-12;
    let t_end = 3.0;
    bs_run(&mut r, &mut bs, id, t_end, 0.05);
    let y = ode_y(&r, id)[0];
    let expect = t_end.sin();
    assert!(
        (y - expect).abs() < 1e-10,
        "BS on y'=cos t at t=3: got {:.17e}, exact sin(3) = {:.17e}",
        y,
        expect
    );
}

#[test]
fn bs_integrates_exponential_decay() {
    let (mut r, mut bs, id) = ode_sim(1, decay_derivatives, &[1.0]);
    bs.eps_abs = 1e-13;
    bs.eps_rel = 1e-13;
    let t_end = 5.0;
    bs_run(&mut r, &mut bs, id, t_end, 0.05);
    let y = ode_y(&r, id)[0];
    let expect = (-t_end).exp();
    assert!(
        ((y - expect) / expect).abs() < 1e-9,
        "BS on y'=-y at t=5: got {:.17e}, exact exp(-5) = {:.17e}",
        y,
        expect
    );
}

#[test]
fn bs_rejects_an_oversized_step_and_reduces_dt() {
    // Two and a half periods in a single step at a tight tolerance can
    // not converge within the extrapolation table.
    let (mut r, mut bs, _id) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
    r.save_messages = 1;
    bs.eps_abs = 1e-10;
    bs.eps_rel = 1e-10;
    let dt = 5.0 * PI;
    let ok = reb_integrator_bs_step_odes(&mut r, &mut bs, dt);
    assert_eq!(ok, 0, "a {:.3}-long step must be rejected", dt);
    assert_eq!(bs.previous_rejected, 1, "rejection must be recorded");
    assert!(
        bs.dt_proposed > 0.0 && bs.dt_proposed < dt,
        "the proposed step {:.3e} must be positive and smaller than {:.3e}",
        bs.dt_proposed,
        dt
    );
    assert_eq!(
        bs.first_or_last_step, 1,
        "a rejected step must not clear first_or_last_step"
    );
}

#[test]
fn bs_min_dt_and_max_dt_clamp_the_proposed_step() {
    // max_dt: the natural proposal for this easy step is well above
    // 1e-4, so the clamp must bite exactly.
    let (mut r, mut bs, _id) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
    r.save_messages = 1;
    bs.eps_abs = 1e-8;
    bs.eps_rel = 1e-8;
    bs.max_dt = 1e-4;
    assert_eq!(reb_integrator_bs_step_odes(&mut r, &mut bs, 1e-4), 1);
    assert_eq!(
        bs.dt_proposed.to_bits(),
        1e-4f64.to_bits(),
        "max_dt must pin dt_proposed, got {:.17e}",
        bs.dt_proposed
    );

    // min_dt: the natural proposal is well below 10, so the floor bites.
    let (mut r2, mut bs2, _id2) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
    r2.save_messages = 1;
    bs2.eps_abs = 1e-8;
    bs2.eps_rel = 1e-8;
    bs2.min_dt = 10.0;
    assert_eq!(reb_integrator_bs_step_odes(&mut r2, &mut bs2, 0.05), 1);
    assert_eq!(
        bs2.dt_proposed.to_bits(),
        10.0f64.to_bits(),
        "min_dt must pin dt_proposed, got {:.17e}",
        bs2.dt_proposed
    );

    // A negative input step must keep the sign of the proposal.
    let (mut r3, mut bs3, _id3) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
    bs3.eps_abs = 1e-8;
    bs3.eps_rel = 1e-8;
    assert_eq!(reb_integrator_bs_step_odes(&mut r3, &mut bs3, -0.05), 1);
    assert!(
        bs3.dt_proposed < 0.0,
        "a backward step must propose a backward dt, got {:.17e}",
        bs3.dt_proposed
    );
}

#[test]
fn bs_order_control_bookkeeping_is_consistent() {
    // cost_per_time_unit[k] is defined as cost_per_step[k]/optimal_step[k];
    // check the identity for every order the step actually visited, and
    // that the retained order stays inside [1, sequence_length-2].
    let (mut r, mut bs, _id) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
    bs.eps_abs = 1e-10;
    bs.eps_rel = 1e-10;
    assert_eq!(reb_integrator_bs_step_odes(&mut r, &mut bs, 0.2), 1);

    let mut visited = 0;
    for k in 0..bs.sequence.len() {
        if bs.optimal_step[k] != 0.0 {
            visited += 1;
            let expect = (bs.cost_per_step[k] as f64) / bs.optimal_step[k];
            assert_eq!(
                bs.cost_per_time_unit[k].to_bits(),
                expect.to_bits(),
                "cost_per_time_unit[{}] = {:.17e} != cost_per_step/optimal_step = {:.17e}",
                k,
                bs.cost_per_time_unit[k],
                expect
            );
            assert!(
                bs.optimal_step[k] > 0.0,
                "optimal_step[{}] must be positive",
                k
            );
        }
    }
    assert!(visited >= 1, "at least one extrapolation order must be used");
    assert!(
        bs.target_iter >= 1 && bs.target_iter as usize <= bs.sequence.len() - 2,
        "target_iter {} outside [1, {}]",
        bs.target_iter,
        bs.sequence.len() - 2
    );
    assert_eq!(
        bs.first_or_last_step, 0,
        "an accepted step clears first_or_last_step"
    );
    assert_eq!(bs.previous_rejected, 0, "no rejection was expected");
}

#[test]
fn bs_order_control_raises_the_order_as_the_tolerance_tightens() {
    // BS trades extrapolation order against step size, so a tighter
    // tolerance does not simply mean a smaller step -- it means a longer
    // extrapolation table. Drive the same oscillator over the same span
    // and watch the retained order climb.
    let mut orders = Vec::new();
    for &tol in &[1e-4f64, 1e-6, 1e-8, 1e-10] {
        let (mut r, mut bs, id) = ode_sim(2, ho_derivatives, &[1.0, 0.0]);
        r.save_messages = 1;
        bs.eps_abs = tol;
        bs.eps_rel = tol;
        bs_run(&mut r, &mut bs, id, 10.0, 0.05);
        assert!(
            bs.target_iter >= 1 && bs.target_iter as usize <= bs.sequence.len() - 2,
            "target_iter {} out of range at tol {:e}",
            bs.target_iter,
            tol
        );
        orders.push(bs.target_iter);
    }
    for i in 1..orders.len() {
        assert!(
            orders[i] > orders[i - 1],
            "BS order should increase with tolerance, got {:?}",
            orders
        );
    }
}

#[test]
fn bs_update_particles_round_trips_the_state_vector() {
    let mut r = two_body(1e-3, 1.0, 0.3, 0.2);
    let before: Vec<reb_particle> = r.particles.clone();
    let mut y = vec![0.0f64; r.N * 6];
    for i in 0..r.N {
        let p = r.particles[i];
        y[i * 6] = p.x;
        y[i * 6 + 1] = p.y;
        y[i * 6 + 2] = p.z;
        y[i * 6 + 3] = p.vx;
        y[i * 6 + 4] = p.vy;
        y[i * 6 + 5] = p.vz;
    }
    // Scramble, then restore from y.
    for i in 0..r.N {
        r.particles[i].x = f64::NAN;
        r.particles[i].vz = f64::NAN;
    }
    reb_integrator_bs_update_particles(&mut r, &y);
    for i in 0..r.N {
        for (got, want, name) in [
            (r.particles[i].x, before[i].x, "x"),
            (r.particles[i].y, before[i].y, "y"),
            (r.particles[i].z, before[i].z, "z"),
            (r.particles[i].vx, before[i].vx, "vx"),
            (r.particles[i].vy, before[i].vy, "vy"),
            (r.particles[i].vz, before[i].vz, "vz"),
        ] {
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "particle {} {} must round trip through the BS state vector",
                i,
                name
            );
        }
    }
}

// ---------------------------------------------------------------------------
// BS — N-body
// ---------------------------------------------------------------------------

fn bs_two_body(tol: f64, e: f64) -> reb_simulation {
    let mut r = two_body(1e-3, 1.0, e, 0.0);
    reb_simulation_set_integrator(&mut r, "bs");
    if let reb_integrator_state::bs(ref mut bs) = r.integrator {
        bs.eps_abs = tol;
        bs.eps_rel = tol;
    } else {
        panic!("integrator is not BS");
    }
    r.dt = 0.01;
    r
}

#[test]
fn bs_conserves_energy_on_a_two_body_problem() {
    let mut r = bs_two_body(1e-13, 0.1);
    let e0 = reb_simulation_energy(&r);
    reb_simulation_integrate(&mut r, 20.0 * 2.0 * PI);
    let drift = rel_energy_drift(&r, e0);
    assert!(
        drift < 1e-9,
        "BS relative energy drift over 20 orbits at tol 1e-13: {:.3e}",
        drift
    );
    assert!(
        (r.t - 20.0 * 2.0 * PI).abs() < 1e-10,
        "BS must finish exactly on tmax, got {:.17e}",
        r.t
    );
}

#[test]
fn bs_energy_drift_shrinks_with_tolerance() {
    let mut drifts = Vec::new();
    for &tol in &[1e-7f64, 1e-11] {
        let mut r = bs_two_body(tol, 0.1);
        let e0 = reb_simulation_energy(&r);
        reb_simulation_integrate(&mut r, 5.0 * 2.0 * PI);
        drifts.push(rel_energy_drift(&r, e0));
    }
    assert!(
        drifts[1] < drifts[0],
        "tightening the BS tolerance must reduce the energy drift: {:?}",
        drifts
    );
    assert!(
        drifts[1] < 1e-9,
        "BS energy drift at tol 1e-11 over 5 orbits: {:.3e}",
        drifts[1]
    );
}

#[test]
fn bs_agrees_with_ias15_on_a_two_body_problem() {
    let tmax = 5.0 * 2.0 * PI;

    let mut a = two_body(1e-3, 1.0, 0.2, 0.0);
    reb_simulation_integrate(&mut a, tmax);

    let mut b = bs_two_body(1e-13, 0.2);
    reb_simulation_integrate(&mut b, tmax);

    let dx = a.particles[1].x - b.particles[1].x;
    let dy = a.particles[1].y - b.particles[1].y;
    let d = (dx * dx + dy * dy).sqrt();
    assert!(
        d < 1e-7,
        "IAS15 vs BS separation after 5 orbits: {:.3e}",
        d
    );

    let oa = reb_orbit_from_particle(a.G, a.particles[1], a.particles[0]);
    let ob = reb_orbit_from_particle(b.G, b.particles[1], b.particles[0]);
    assert!(
        ((oa.a - ob.a) / oa.a).abs() < 1e-10,
        "IAS15 vs BS semi-major axis: {:.17e} vs {:.17e}",
        oa.a,
        ob.a
    );
    assert!(
        (oa.e - ob.e).abs() < 1e-10,
        "IAS15 vs BS eccentricity: {:.17e} vs {:.17e}",
        oa.e,
        ob.e
    );
}

#[test]
fn bs_is_bit_deterministic() {
    let run = || {
        let mut r = bs_two_body(1e-11, 0.35);
        reb_simulation_integrate(&mut r, 13.0);
        r
    };
    let a = run();
    let b = run();
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical BS runs must be bit identical"
    );
}

#[test]
fn bs_handles_a_high_eccentricity_orbit() {
    let mut r = bs_two_body(1e-12, 0.9);
    let e0 = reb_simulation_energy(&r);
    let o0 = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);
    reb_simulation_integrate(&mut r, 5.0 * 2.0 * PI);
    let drift = rel_energy_drift(&r, e0);
    assert!(
        drift < 1e-7,
        "BS relative energy drift over 5 orbits at e=0.9: {:.3e}",
        drift
    );
    let o1 = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);
    assert!(
        ((o1.a - o0.a) / o0.a).abs() < 1e-7,
        "BS semi-major axis at e=0.9: {:.17e} -> {:.17e}",
        o0.a,
        o1.a
    );
}
