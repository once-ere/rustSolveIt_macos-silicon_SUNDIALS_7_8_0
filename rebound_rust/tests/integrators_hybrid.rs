//! Integration tests for the integrators_hybrid module group of rebound_rs.
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

use std::f64::consts::PI;

use integrator_eos as eos;
use integrator_janus as janus;
use integrator_mercurius as merc;
use integrator_trace as trace;

// ---------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------

/// A simulation with messages captured instead of printed to stderr.
fn quiet_sim() -> reb_simulation {
    let mut r = reb_simulation_create();
    r.save_messages = 1;
    r
}

/// Build a particle from raw state (mass, radius, position, velocity).
fn mkp(m: f64, rad: f64, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> reb_particle {
    let mut p = reb_particle::default();
    p.m = m;
    p.r = rad;
    p.x = x;
    p.y = y;
    p.z = z;
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    p
}

/// Raw bit pattern of every particle coordinate plus the time, for
/// bit-identity (determinism / reversibility) comparisons.
fn state_bits(r: &reb_simulation) -> Vec<u64> {
    let mut v = Vec::with_capacity(r.N * 8 + 1);
    for p in &r.particles {
        for q in [p.x, p.y, p.z, p.vx, p.vy, p.vz, p.m, p.r] {
            v.push(q.to_bits());
        }
    }
    v.push(r.t.to_bits());
    v
}

/// A mercurius state whose internal arrays are already sized for `n`
/// particles (the sizing `reb_integrator_mercurius_step_state` performs
/// on its first call).
fn mercurius_state_for(n: usize) -> merc::reb_integrator_mercurius_state {
    let mut m = merc::reb_integrator_mercurius_state::default();
    m.dcrit = vec![0.; n];
    m.particles_backup = vec![reb_particle::default(); n];
    m.encounter_map = vec![0; n];
    m.N_allocated = n;
    m
}

/// A trace state whose internal arrays are already sized for `n`
/// particles (as `reb_integrator_trace_step_state` sizes them).
fn trace_state_for(n: usize) -> trace::reb_integrator_trace_state {
    let mut t = trace::reb_integrator_trace_state::default();
    t.particles_backup = vec![reb_particle::default(); n];
    t.particles_backup_kepler = vec![reb_particle::default(); n];
    t.current_Ks = vec![0; n * n];
    t.encounter_map = vec![0; n];
    t.encounter_map_backup = vec![0; n];
    t.N_allocated = n;
    t
}

/// Star + two planets, heliocentric orbital elements, moved to the COM.
fn two_planet_sim(integrator: &str, dt: f64, m_pl: f64) -> reb_simulation {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_set_integrator(&mut r, integrator);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0., 0., 0., 0., 0., 0.));
    let star = r.particles[0];
    let p1 = reb_particle_from_orbit(r.G, star, m_pl, 1.0, 0.05, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, p1);
    let p2 = reb_particle_from_orbit(r.G, star, m_pl, 2.3, 0.02, 0.01, 0.0, 0.3, 1.7);
    reb_simulation_add(&mut r, p2);
    reb_simulation_move_to_com(&mut r);
    r
}

/// Relative energy drift of `r` after integrating to `tmax`.
fn energy_drift(r: &mut reb_simulation, tmax: f64) -> f64 {
    let e0 = reb_simulation_energy(r);
    reb_simulation_integrate(r, tmax);
    let e1 = reb_simulation_energy(r);
    ((e1 - e0) / e0).abs()
}

// =====================================================================
// MERCURIUS: changeover (switching) functions
// =====================================================================

// The four MERCURIUS changeover functions all map d -> L in [0,1] via
// y = (d - 0.1*dcrit)/(0.9*dcrit).  With dcrit = 10 both 0.1*dcrit and
// 0.9*dcrit are exact in binary (1.0 and 9.0), so y is exact too:
//   d = 1.0  -> y = 0     d = 5.5 -> y = 0.5     d = 10.0 -> y = 1
const DCRIT_EXACT: f64 = 10.0;

#[test]
fn test_mercurius_L_endpoints_are_exactly_zero_and_one() {
    let r = quiet_sim();
    let fns: [(&str, merc::reb_integrator_mercurius_Lfunc); 4] = [
        ("mercury", merc::reb_integrator_mercurius_L_mercury),
        ("C4", merc::reb_integrator_mercurius_L_C4),
        ("C5", merc::reb_integrator_mercurius_L_C5),
        ("infinity", merc::reb_integrator_mercurius_L_infinity),
    ];
    for (name, L) in fns {
        // Below the switch-on radius 0.1*dcrit: fully handled by the
        // close-encounter integrator, L == 0.
        assert_eq!(L(&r, 0.5, DCRIT_EXACT), 0.0, "L_{} below 0.1*dcrit", name);
        // y == 0 exactly: every polynomial has a zero of high order at 0.
        assert_eq!(L(&r, 1.0, DCRIT_EXACT), 0.0, "L_{} at y=0", name);
        // y == 1 exactly: the polynomials are normalised to 1 there.
        assert_eq!(L(&r, 10.0, DCRIT_EXACT), 1.0, "L_{} at y=1", name);
        // Beyond dcrit: fully handled by the WH kick, L == 1.
        assert_eq!(L(&r, 20.0, DCRIT_EXACT), 1.0, "L_{} above dcrit", name);
    }
}

#[test]
fn test_mercurius_L_midpoint_is_exactly_one_half() {
    // y = 1/2 is a dyadic rational, so each polynomial evaluates in
    // exact binary arithmetic:
    //   mercury: 10/8 - 15/16 + 6/32                      = 1/2
    //   C4:      (4.375-39.375+135-210+126)/2^5 = 16/32    = 1/2
    //   C5:      (-7.875+86.625-385+866.25-990+462)/2^6    = 1/2
    //   infty:   f(1/2)/(f(1/2)+f(1/2))                    = 1/2
    let r = quiet_sim();
    let d = 5.5; // (5.5 - 1.0)/9.0 == 0.5 exactly
    assert_eq!(
        merc::reb_integrator_mercurius_L_mercury(&r, d, DCRIT_EXACT),
        0.5,
        "L_mercury at y=1/2"
    );
    assert_eq!(
        merc::reb_integrator_mercurius_L_C4(&r, d, DCRIT_EXACT),
        0.5,
        "L_C4 at y=1/2"
    );
    assert_eq!(
        merc::reb_integrator_mercurius_L_C5(&r, d, DCRIT_EXACT),
        0.5,
        "L_C5 at y=1/2"
    );
    assert_eq!(
        merc::reb_integrator_mercurius_L_infinity(&r, d, DCRIT_EXACT),
        0.5,
        "L_infinity at y=1/2"
    );
}

#[test]
fn test_mercurius_L_are_symmetric_and_monotone() {
    // Every changeover here is a symmetric smoothstep: L(y) + L(1-y) = 1
    // and L is non-decreasing.  Both follow from the polynomial
    // coefficients (regularised incomplete beta functions) and, for
    // L_infinity, from f(y)/(f(y)+f(1-y)) + f(1-y)/(f(1-y)+f(y)) = 1.
    //
    // In double precision the identity is limited by cancellation in the
    // alternating C5 polynomial, whose largest coefficient is 3465: the
    // absolute error near y = 1 is ~3465*2^-53 ~ 4e-13.  TOL is set an
    // order of magnitude above that; SLACK covers the same round-off in
    // the monotonicity check, where the true increments vanish at the
    // endpoints.
    const TOL: f64 = 1e-11;
    const SLACK: f64 = 1e-12;
    let r = quiet_sim();
    let fns: [(&str, merc::reb_integrator_mercurius_Lfunc); 4] = [
        ("mercury", merc::reb_integrator_mercurius_L_mercury),
        ("C4", merc::reb_integrator_mercurius_L_C4),
        ("C5", merc::reb_integrator_mercurius_L_C5),
        ("infinity", merc::reb_integrator_mercurius_L_infinity),
    ];
    for (name, L) in fns {
        let mut prev = -1.0f64;
        for k in 0..=200 {
            let y = (k as f64) / 200.0;
            let d = 1.0 + 9.0 * y; // y = (d-1)/9
            let d_mirror = 1.0 + 9.0 * (1.0 - y);
            let a = L(&r, d, DCRIT_EXACT);
            let b = L(&r, d_mirror, DCRIT_EXACT);
            assert!(
                (a + b - 1.0).abs() < TOL,
                "L_{} symmetry at y={}: {} + {} != 1",
                name,
                y,
                a,
                b
            );
            assert!(
                (0.0..=1.0).contains(&a),
                "L_{} out of range at y={}: {}",
                name,
                y,
                a
            );
            assert!(
                a >= prev - SLACK,
                "L_{} not monotone at y={}: {} follows {}",
                name,
                y,
                a,
                prev
            );
            prev = a;
        }
    }
}

#[test]
fn test_mercurius_L_infinity_matches_closed_form() {
    // L_infinity(y) = e^{-1/y} / (e^{-1/y} + e^{-1/(1-y)}), evaluated
    // here by an independent expression.
    let r = quiet_sim();
    for k in 1..100 {
        let y = (k as f64) / 100.0;
        let d = 1.0 + 9.0 * y;
        let got = merc::reb_integrator_mercurius_L_infinity(&r, d, DCRIT_EXACT);
        let fy = (-1.0 / y).exp();
        let f1y = (-1.0 / (1.0 - y)).exp();
        let want = fy / (fy + f1y);
        assert!(
            (got - want).abs() < 1e-14,
            "L_infinity closed form at y={}: {} vs {}",
            y,
            got,
            want
        );
        // Strictly inside (0,1) only where both exponentials are
        // representable relative to each other: e^{-1/y}/e^{-1/(1-y)}
        // exceeds 1/eps once y > ~0.85, and the quotient then rounds to
        // exactly 1 (and to exactly 0 for y < ~0.15).
        if (0.2..=0.8).contains(&y) {
            assert!(
                got > 0.0 && got < 1.0,
                "L_infinity strictly inside (0,1) at y={}: {}",
                y,
                got
            );
        }
    }
}

#[test]
fn test_mercurius_L_C4_C5_switch_on_later_than_mercury() {
    // C4 and C5 have higher-order zeros at y=0 (y^5 and y^6 versus y^3),
    // so on the lower half of the changeover they are strictly smaller
    // than the Mercury changeover, and by symmetry strictly larger on
    // the upper half.
    let r = quiet_sim();
    for k in 1..100 {
        let y = (k as f64) / 200.0; // y in (0, 0.5)
        let d = 1.0 + 9.0 * y;
        let lm = merc::reb_integrator_mercurius_L_mercury(&r, d, DCRIT_EXACT);
        let l4 = merc::reb_integrator_mercurius_L_C4(&r, d, DCRIT_EXACT);
        let l5 = merc::reb_integrator_mercurius_L_C5(&r, d, DCRIT_EXACT);
        assert!(l4 < lm, "L_C4 < L_mercury on lower half at y={}", y);
        assert!(l5 < l4, "L_C5 < L_C4 on lower half at y={}", y);
    }
}

// =====================================================================
// MERCURIUS: dcrit criteria
// =====================================================================

#[test]
fn test_mercurius_dcrit_criterion1_average_velocity() {
    // Circular massless orbit: r = 1, v^2 = 1, GM = 1 =>
    // a = GM*r/(2GM - r v^2) = 1, vc = sqrt(GM/|a|) = 1 exactly.
    // Criteria 3 and 4 vanish (m = 0, radius = 0), criterion 2 equals
    // criterion 1, so dcrit = 1.0 * 0.4 * dt.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.03125; // exact binary
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.0, 0.0, 1., 0., 0., 0., 1., 0.));
    let m = mercurius_state_for(2);
    let dcrit = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m, 1);
    let want = 1.0f64 * 0.4 * r.dt;
    assert_eq!(dcrit, want, "dcrit from circular velocity criterion");
}

#[test]
fn test_mercurius_dcrit_criterion2_current_velocity() {
    // Hyperbolic: r = 1, v = 2 => v^2 = 4, GM = 1, a = 1/(2-4) = -0.5,
    // vc = sqrt(1/0.5) = sqrt(2) < 2.  Criterion 2 (current speed)
    // therefore wins: dcrit = 2 * 0.4 * dt.  Criterion 3 is negative
    // (a < 0, m = 0) and criterion 4 is zero.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.0625;
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.0, 0.0, 1., 0., 0., 0., 2., 0.));
    let m = mercurius_state_for(2);
    let dcrit = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m, 1);
    let want = 2.0f64 * 0.4 * r.dt;
    assert_eq!(dcrit, want, "dcrit from current-velocity criterion");
    assert!(want > (2.0f64).sqrt() * 0.4 * r.dt, "criterion 2 dominates 1");
}

#[test]
fn test_mercurius_dcrit_criterion3_hill_radius() {
    // Circular orbit with GM = G(m0+m): v^2 = GM/r gives a = r = 1, so
    // criterion 3 = r_crit_hill * a * cbrt(m/(3 m0)).  With m = 1e-3 and
    // m0 = 1 that is 3 * cbrt(1e-3/3) ~ 0.208, which dwarfs the velocity
    // criteria (~0.4*dt = 4e-3) and the (zero) radius criterion.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.01;
    let m_pl = 1e-3;
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0., 0., 0., 0., 0., 0.));
    let v = (1.0f64 + m_pl).sqrt();
    reb_simulation_add(&mut r, mkp(m_pl, 0.0, 1., 0., 0., 0., v, 0.));
    let m = mercurius_state_for(2);
    let dcrit = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m, 1);
    let want = 3.0 * 1.0 * (m_pl / 3.0).cbrt();
    assert!(
        (dcrit - want).abs() < 1e-12 * want,
        "dcrit Hill criterion: {} vs {}",
        dcrit,
        want
    );
    // The Hill criterion scales linearly with r_crit_hill.
    let mut m6 = mercurius_state_for(2);
    m6.r_crit_hill = 6.0;
    let dcrit6 = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m6, 1);
    assert!(
        (dcrit6 - 2.0 * dcrit).abs() < 1e-12 * dcrit,
        "dcrit doubles when r_crit_hill doubles"
    );
}

#[test]
fn test_mercurius_dcrit_criterion4_physical_radius() {
    // A large physical radius sets the floor: dcrit = 2*r.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 1e-6; // makes the velocity criteria negligible
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.0, 0.5, 1., 0., 0., 0., 1., 0.));
    let m = mercurius_state_for(2);
    let dcrit = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m, 1);
    assert_eq!(dcrit, 1.0, "dcrit from physical radius (2*0.5)");
}

#[test]
fn test_mercurius_dcrit_is_the_max_of_all_four_criteria() {
    // Independently evaluate the four criteria of
    // reb_integrator_mercurius_calculate_dcrit_for_particle from the
    // vis-viva relation and check the routine returns their maximum.
    let mut r = quiet_sim();
    r.G = 1.3;
    r.dt = 0.017;
    reb_simulation_add(&mut r, mkp(2.0, 0.03, 0., 0., 0., 0.1, -0.2, 0.05));
    reb_simulation_add(&mut r, mkp(4e-3, 0.011, 0.7, -0.3, 0.2, 0.4, 1.1, -0.15));
    let mut m = mercurius_state_for(2);
    m.r_crit_hill = 4.0;
    let got = merc::reb_integrator_mercurius_calculate_dcrit_for_particle(&r, &m, 1);

    let p = r.particles[1];
    let s = r.particles[0];
    let d = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    let dvx = p.vx - s.vx;
    let dvy = p.vy - s.vy;
    let dvz = p.vz - s.vz;
    let v2 = dvx * dvx + dvy * dvy + dvz * dvz;
    let GM = r.G * (s.m + p.m);
    // vis-viva: 1/a = 2/d - v^2/GM
    let a = 1.0 / (2.0 / d - v2 / GM);
    let c1 = (GM / a.abs()).sqrt() * 0.4 * r.dt;
    let c2 = v2.sqrt() * 0.4 * r.dt;
    let c3 = m.r_crit_hill * a * (p.m / (3.0 * s.m)).cbrt();
    let c4 = 2.0 * p.r;
    let want = c1.max(c2).max(c3).max(c4);
    assert!(
        (got - want).abs() < 1e-12 * want,
        "dcrit = max(criteria): {} vs {} (c1={},c2={},c3={},c4={})",
        got,
        want,
        c1,
        c2,
        c3,
        c4
    );
}

// =====================================================================
// MERCURIUS: democratic-heliocentric transform
// =====================================================================

#[test]
fn test_mercurius_inertial_to_dh_and_back_roundtrip() {
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0.3, -0.2, 0.1, 0.01, -0.02, 0.03));
    reb_simulation_add(&mut r, mkp(1e-3, 0., 1.4, 0.5, -0.2, -0.3, 0.9, 0.05));
    reb_simulation_add(&mut r, mkp(2e-4, 0., -2.1, 1.1, 0.4, -0.5, -0.6, -0.02));
    let before: Vec<reb_particle> = r.particles.clone();

    let mut m = mercurius_state_for(3);
    merc::reb_integrator_mercurius_inertial_to_dh(&mut r, &mut m);

    // Particle 0 sits at the origin in democratic-heliocentric
    // coordinates and the barycentric velocity has been removed.
    assert_eq!(r.particles[0].x, 0.0, "dh puts particle 0 at x=0");
    assert_eq!(r.particles[0].y, 0.0, "dh puts particle 0 at y=0");
    assert_eq!(r.particles[0].z, 0.0, "dh puts particle 0 at z=0");
    for i in 1..3 {
        assert!(
            (r.particles[i].x - (before[i].x - before[0].x)).abs() < 1e-15,
            "dh position of particle {} is heliocentric",
            i
        );
    }
    // Total momentum in dh velocities is zero.
    let mut px = 0.0;
    for i in 0..3 {
        px += r.particles[i].m * r.particles[i].vx;
    }
    assert!(px.abs() < 1e-15, "dh velocities carry no net momentum");

    merc::reb_integrator_mercurius_dh_to_inertial(&mut r, &mut m);
    for i in 0..3 {
        let a = r.particles[i];
        let b = before[i];
        assert!((a.x - b.x).abs() < 1e-14, "roundtrip x of particle {}", i);
        assert!((a.y - b.y).abs() < 1e-14, "roundtrip y of particle {}", i);
        assert!((a.z - b.z).abs() < 1e-14, "roundtrip z of particle {}", i);
        assert!((a.vx - b.vx).abs() < 1e-14, "roundtrip vx of particle {}", i);
        assert!((a.vy - b.vy).abs() < 1e-14, "roundtrip vy of particle {}", i);
        assert!((a.vz - b.vz).abs() < 1e-14, "roundtrip vz of particle {}", i);
    }
}

#[test]
fn test_mercurius_dh_com_matches_simulation_com() {
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0.3, -0.2, 0.1, 0.01, -0.02, 0.03));
    reb_simulation_add(&mut r, mkp(1e-2, 0., 1.4, 0.5, -0.2, -0.3, 0.9, 0.05));
    reb_simulation_add(&mut r, mkp(5e-3, 0., -2.1, 1.1, 0.4, -0.5, -0.6, -0.02));
    let com = reb_simulation_com(&r);
    let mut m = mercurius_state_for(3);
    merc::reb_integrator_mercurius_inertial_to_dh(&mut r, &mut m);
    assert!(
        (m.com_pos.x - com.x).abs() < 1e-15
            && (m.com_pos.y - com.y).abs() < 1e-15
            && (m.com_pos.z - com.z).abs() < 1e-15,
        "mercurius com_pos equals reb_simulation_com position"
    );
    assert!(
        (m.com_vel.x - com.vx).abs() < 1e-15
            && (m.com_vel.y - com.vy).abs() < 1e-15
            && (m.com_vel.z - com.vz).abs() < 1e-15,
        "mercurius com_vel equals reb_simulation_com velocity"
    );
}

// =====================================================================
// MERCURIUS: individual sub-steps
// =====================================================================

#[test]
fn test_mercurius_jump_step_exact_value() {
    // Jump step drifts every non-central particle by dt * P/m0, where P
    // is the total dh momentum of the (massive) particles.
    // m1 = 0.5, vx = 2 => P = 1, m0 = 1 => shift = dt*1 = 0.25.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.25;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.5, 0., 1., 2., 4., 2., 0., 0.));
    let dt = r.dt;
    merc::reb_integrator_mercurius_jump_step(&mut r, dt);
    assert_eq!(r.particles[1].x, 1.25, "jump step x");
    assert_eq!(r.particles[1].y, 2.0, "jump step leaves y alone (vy=0)");
    assert_eq!(r.particles[1].z, 4.0, "jump step leaves z alone (vz=0)");
    assert_eq!(r.particles[0].x, 0.0, "jump step never moves the star");
    assert_eq!(r.particles[1].vx, 2.0, "jump step never changes velocity");
}

#[test]
fn test_mercurius_interaction_step_exact_value() {
    // v_i += dt * a_i for i >= 1; particle 0 is untouched.
    let mut r = quiet_sim();
    r.dt = 0.5;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0., 1., 0., 0., 0., 0., 0.));
    r.particles[0].ax = 8.0;
    r.particles[1].ax = 4.0;
    r.particles[1].ay = -2.0;
    r.particles[1].az = 1.0;
    merc::reb_integrator_mercurius_interaction_step(&mut r, 0.5);
    assert_eq!(r.particles[1].vx, 2.0, "interaction step vx");
    assert_eq!(r.particles[1].vy, -1.0, "interaction step vy");
    assert_eq!(r.particles[1].vz, 0.5, "interaction step vz");
    assert_eq!(r.particles[0].vx, 0.0, "interaction step skips particle 0");
}

#[test]
fn test_mercurius_encounter_acceleration_full_pair_force_when_L_is_zero() {
    // In encounter mode the pair force enters with weight (1-L).  With
    // dcrit = 100 and separation 1 the changeover argument y is negative
    // so L = 0 and the full pair force is applied on top of the stellar
    // term.  Star at origin, m0 = m1 = m2 = 1, G = 1:
    //   a1 = -G m0 x1/|x1|^3 + G m2 (x2-x1)/|x2-x1|^3 = -1 + 1 = 0
    //   a2 = -G m0 x2/|x2|^3 - G m1 (x2-x1)/|x2-x1|^3 = -0.25 - 1
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0., 1., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0., 2., 0., 0., 0., 0., 0.));
    let mut m = mercurius_state_for(3);
    m.dcrit = vec![100.0, 100.0, 100.0];
    m.encounter_map = vec![0, 1, 2];
    m.encounter_N = 3;
    m.encounter_N_active = 3;
    r.integrator = reb_integrator_state::mercurius(m);
    merc::reb_integrator_mercurius_calculate_acceleration_mode_encounter(&mut r);
    assert_eq!(r.particles[0].ax, 0.0, "star feels no acceleration in dh");
    assert_eq!(r.particles[1].ax, 0.0, "a1 with L=0");
    assert_eq!(r.particles[2].ax, -1.25, "a2 with L=0");
}

#[test]
fn test_mercurius_encounter_acceleration_drops_pair_force_when_L_is_one() {
    // With dcrit = 0.1 and separation 1 the changeover argument is
    // y = (1 - 0.01)/0.09 = 11 > 1, so L = 1 and (1-L) = 0: only the
    // stellar term survives.
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0., 1., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0., 2., 0., 0., 0., 0., 0.));
    let mut m = mercurius_state_for(3);
    m.dcrit = vec![0.1, 0.1, 0.1];
    m.encounter_map = vec![0, 1, 2];
    m.encounter_N = 3;
    m.encounter_N_active = 3;
    r.integrator = reb_integrator_state::mercurius(m);
    merc::reb_integrator_mercurius_calculate_acceleration_mode_encounter(&mut r);
    assert_eq!(r.particles[1].ax, -1.0, "a1 = stellar term only");
    assert_eq!(r.particles[2].ax, -0.25, "a2 = stellar term only");
}

// =====================================================================
// MERCURIUS: full integrations
// =====================================================================

#[test]
fn test_mercurius_single_particle_drifts_at_constant_velocity() {
    // N = 1 edge case: the whole step reduces to the centre-of-mass
    // drift, so a lone star must move ballistically.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.1;
    reb_simulation_set_integrator(&mut r, "mercurius");
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0.3, -0.4, 0.0));
    reb_simulation_integrate(&mut r, 1.0);
    assert!((r.t - 1.0).abs() < 1e-14, "reached tmax");
    assert!(
        (r.particles[0].x - 0.3 * r.t).abs() < 1e-14,
        "lone star x: {} vs {}",
        r.particles[0].x,
        0.3 * r.t
    );
    assert!(
        (r.particles[0].y + 0.4 * r.t).abs() < 1e-14,
        "lone star y: {}",
        r.particles[0].y
    );
    assert_eq!(r.particles[0].vx, 0.3, "lone star velocity is unchanged");
    assert_eq!(r.particles[0].vy, -0.4, "lone star velocity is unchanged");
}

#[test]
fn test_mercurius_conserves_energy_for_two_planets() {
    // No close encounters: MERCURIUS is a democratic-heliocentric
    // WH map, whose relative energy error is O((m/M)(dt/P)^2)
    // ~ 1e-5 * (1/200)^2 ~ 2.5e-10 and does not drift secularly.
    let mut r = two_planet_sim("mercurius", 2.0 * PI / 200.0, 1e-5);
    let drift = energy_drift(&mut r, 100.0);
    assert!(
        drift < 1e-8,
        "MERCURIUS relative energy drift {} exceeds 1e-8",
        drift
    );
}

#[test]
fn test_mercurius_conserves_angular_momentum_for_two_planets() {
    // Angular momentum is conserved exactly by every operator of the
    // splitting (drift, Kepler, kick), so only round-off accumulates.
    let mut r = two_planet_sim("mercurius", 2.0 * PI / 200.0, 1e-5);
    let l0 = reb_simulation_angular_momentum(&r);
    reb_simulation_integrate(&mut r, 100.0);
    let l1 = reb_simulation_angular_momentum(&r);
    let n0 = (l0.x * l0.x + l0.y * l0.y + l0.z * l0.z).sqrt();
    let dl = ((l1.x - l0.x).powi(2) + (l1.y - l0.y).powi(2) + (l1.z - l0.z).powi(2)).sqrt();
    assert!(
        dl / n0 < 1e-12,
        "MERCURIUS relative angular momentum drift {}",
        dl / n0
    );
}

#[test]
fn test_mercurius_forced_encounter_branch_still_conserves_energy() {
    // A huge r_crit_hill forces every pair into the IAS15 encounter
    // sub-step, exercising encounter prediction, the encounter map and
    // the (1-L) encounter gravity every timestep.
    let mut r = two_planet_sim("mercurius", 2.0 * PI / 200.0, 1e-5);
    if let reb_integrator_state::mercurius(ref mut m) = r.integrator {
        m.r_crit_hill = 5000.0;
    }
    let e0 = reb_simulation_energy(&r);
    reb_simulation_integrate(&mut r, 20.0);
    let e1 = reb_simulation_energy(&r);
    let drift = ((e1 - e0) / e0).abs();
    // The encounter machinery must actually have engaged.
    if let reb_integrator_state::mercurius(ref m) = r.integrator {
        assert!(
            m.encounter_N >= 3,
            "expected all 3 particles flagged as encountering, got {}",
            m.encounter_N
        );
        assert_eq!(m.tponly_encounter, 0, "massive-massive encounter flagged");
    } else {
        panic!("integrator is no longer mercurius");
    }
    assert!(
        drift < 1e-9,
        "MERCURIUS forced-encounter energy drift {}",
        drift
    );
}

#[test]
fn test_mercurius_agrees_with_ias15_on_a_two_body_problem() {
    // Same Kepler problem, two very different integrators.
    let mut a = quiet_sim();
    a.G = 1.0;
    a.dt = 2.0 * PI / 2000.0;
    reb_simulation_set_integrator(&mut a, "mercurius");
    reb_simulation_add(&mut a, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    let star = a.particles[0];
    let p = reb_particle_from_orbit(1.0, star, 1e-6, 1.0, 0.3, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut a, p);
    reb_simulation_move_to_com(&mut a);

    let mut b = quiet_sim();
    b.G = 1.0;
    b.dt = 2.0 * PI / 2000.0;
    reb_simulation_set_integrator(&mut b, "ias15");
    reb_simulation_add(&mut b, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut b, p);
    reb_simulation_move_to_com(&mut b);

    reb_simulation_integrate(&mut a, 4.0 * PI);
    reb_simulation_integrate(&mut b, 4.0 * PI);
    let dx = a.particles[1].x - b.particles[1].x;
    let dy = a.particles[1].y - b.particles[1].y;
    let dz = a.particles[1].z - b.particles[1].z;
    let d = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        d < 1e-6,
        "MERCURIUS vs IAS15 separation after 2 orbits: {}",
        d
    );
}

#[test]
fn test_mercurius_is_deterministic() {
    let mut a = two_planet_sim("mercurius", 0.03, 1e-4);
    let mut b = two_planet_sim("mercurius", 0.03, 1e-4);
    reb_simulation_integrate(&mut a, 12.0);
    reb_simulation_integrate(&mut b, 12.0);
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical MERCURIUS runs must agree bit for bit"
    );
}

#[test]
fn test_mercurius_handles_high_eccentricity_orbit() {
    // e = 0.9: the Kepler solver has to cope with a sharp pericentre
    // passage.  Energy must still be conserved by the WH map.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 2.0 * PI / 4000.0;
    reb_simulation_set_integrator(&mut r, "mercurius");
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    let star = r.particles[0];
    let p = reb_particle_from_orbit(1.0, star, 1e-6, 1.0, 0.9, 0.0, 0.0, 0.0, PI);
    reb_simulation_add(&mut r, p);
    reb_simulation_move_to_com(&mut r);
    let drift = energy_drift(&mut r, 6.0 * PI);
    assert!(drift < 1e-7, "e=0.9 MERCURIUS energy drift {}", drift);
    let o = reb_orbit_from_particle(r.G, r.particles[1], r.particles[0]);
    assert!((o.a - 1.0).abs() < 1e-6, "semi-major axis preserved: {}", o.a);
    assert!((o.e - 0.9).abs() < 1e-6, "eccentricity preserved: {}", o.e);
}

#[test]
fn test_mercurius_handles_retrograde_orbit() {
    // inc = pi: the orbit normal points along -z, so Lz must be negative
    // and stay so.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 2.0 * PI / 500.0;
    reb_simulation_set_integrator(&mut r, "mercurius");
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    let star = r.particles[0];
    let p = reb_particle_from_orbit(1.0, star, 1e-5, 1.0, 0.1, PI, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, p);
    reb_simulation_move_to_com(&mut r);
    let l0 = reb_simulation_angular_momentum(&r);
    assert!(l0.z < 0.0, "retrograde orbit has Lz < 0: {}", l0.z);
    let drift = energy_drift(&mut r, 20.0);
    assert!(drift < 1e-9, "retrograde MERCURIUS energy drift {}", drift);
    let l1 = reb_simulation_angular_momentum(&r);
    assert!(
        (l1.z - l0.z).abs() < 1e-14 * l0.z.abs(),
        "Lz preserved for a retrograde orbit"
    );
}

// =====================================================================
// TRACE: body-body switching function
// =====================================================================

/// Star (m=1, radius rs) + two bodies; returns the simulation and a
/// sized trace state.
fn trace_pair(
    rs: f64,
    m1: f64,
    p1: [f64; 6],
    m2: f64,
    p2: [f64; 6],
    dt: f64,
) -> (reb_simulation, trace::reb_integrator_trace_state) {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_add(&mut r, mkp(1.0, rs, 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(m1, 0., p1[0], p1[1], p1[2], p1[3], p1[4], p1[5]));
    reb_simulation_add(&mut r, mkp(m2, 0., p2[0], p2[1], p2[2], p2[3], p2[4], p2[5]));
    let t = trace_state_for(3);
    (r, t)
}

#[test]
fn test_trace_switch_default_triggers_inside_the_hill_sphere() {
    // Body 1 has m = 3e-3 at distance 1 from a unit-mass star, so its
    // modified Hill radius is 1 * (3e-3/3)^(1/3) = 0.1 and, with
    // r_crit_hill = 3, dcritmax = 0.3.  A companion 0.2 away is inside;
    // one 0.5 away (with zero relative velocity, so no look-ahead) is
    // not.
    let (r, t) = trace_pair(
        0.0,
        3e-3,
        [1., 0., 0., 0., 0., 0.],
        0.0,
        [1.2, 0., 0., 0., 0., 0.],
        0.01,
    );
    assert_eq!(
        trace::reb_integrator_trace_switch_default(&r, &t, 1, 2),
        1,
        "separation 0.2 < 0.3 must switch on"
    );

    let (r, t) = trace_pair(
        0.0,
        3e-3,
        [1., 0., 0., 0., 0., 0.],
        0.0,
        [1.5, 0., 0., 0., 0., 0.],
        0.01,
    );
    assert_eq!(
        trace::reb_integrator_trace_switch_default(&r, &t, 1, 2),
        0,
        "separation 0.5 > 0.3 with zero relative velocity must not switch"
    );
}

#[test]
fn test_trace_switch_default_lookahead_over_half_a_timestep() {
    // Body 1 (m=3e-3, distance 1) has dcritmax = 0.3, dcritmax^6 = 7.29e-4.
    // The companion sits 0.5 away and closes at unit speed, so at time
    // tau the squared separation is (0.5 - tau)^2.  The routine samples
    // it at tau = min(t_closest, dt/2) = min(0.5, dt/2):
    //   dt = 2.0  -> tau = 0.5  -> dmin^2 = 0      -> switch on
    //   dt = 0.5  -> tau = 0.25 -> dmin^2 = 0.0625 -> 2.44e-4 < 7.29e-4 -> on
    //   dt = 0.2  -> tau = 0.1  -> dmin^2 = 0.16   -> 4.10e-3 > 7.29e-4 -> off
    for (dt, want) in [(2.0, 1), (0.5, 1), (0.2, 0)] {
        let (r, t) = trace_pair(
            0.0,
            3e-3,
            [1., 0., 0., 0., 0., 0.],
            0.0,
            [1.5, 0., 0., -1., 0., 0.],
            dt,
        );
        assert_eq!(
            trace::reb_integrator_trace_switch_default(&r, &t, 1, 2),
            want,
            "look-ahead with dt = {}",
            dt
        );
    }
}

#[test]
fn test_trace_switch_default_is_symmetric_in_its_two_bodies() {
    // For i, j >= 1 the formula only involves |dx|^2, |dv|^2 and dx.dv,
    // all invariant under swapping i and j, and MAX(dcriti6, dcritj6) is
    // symmetric too.
    for k in 0..40 {
        let s = 0.05 + 0.02 * (k as f64);
        let (r, t) = trace_pair(
            0.0,
            3e-3,
            [1., 0., 0., 0., 0.4, 0.1],
            2e-3,
            [1. + s, 0.1, -0.05, -0.7, 0.2, -0.3],
            0.3,
        );
        let a = trace::reb_integrator_trace_switch_default(&r, &t, 1, 2);
        let b = trace::reb_integrator_trace_switch_default(&r, &t, 2, 1);
        assert_eq!(a, b, "switch(1,2) == switch(2,1) at separation {}", s);
    }
}

#[test]
fn test_trace_switch_default_is_invariant_under_velocity_reversal() {
    // The direction flag d = -sign(dx.dv) enters the look-ahead only
    // through the products -d*qv and 2*d*qv; both are invariant when all
    // velocities flip sign.  This is exactly what makes the pre/post
    // timestep check reversible.
    for k in 0..60 {
        let s = 0.02 + 0.02 * (k as f64);
        let fwd = trace_pair(
            0.0,
            3e-3,
            [1., 0., 0., 0.05, 0.4, 0.1],
            2e-3,
            [1. + s, 0.1, -0.05, -0.7, 0.2, -0.3],
            0.4,
        );
        let rev = trace_pair(
            0.0,
            3e-3,
            [1., 0., 0., -0.05, -0.4, -0.1],
            2e-3,
            [1. + s, 0.1, -0.05, 0.7, -0.2, 0.3],
            0.4,
        );
        assert_eq!(
            trace::reb_integrator_trace_switch_default(&fwd.0, &fwd.1, 1, 2),
            trace::reb_integrator_trace_switch_default(&rev.0, &rev.1, 1, 2),
            "velocity reversal changed the switch at separation {}",
            s
        );
    }
}

#[test]
fn test_trace_switch_default_uses_only_the_physical_radius_for_the_star() {
    // For i == 0 the criterion is the star's physical radius (times
    // r_crit_hill), never a Hill radius: rs = 0.2 gives dcritmax = 0.6.
    let (r, t) = trace_pair(
        0.2,
        0.0,
        [0.5, 0., 0., 0., 0., 0.],
        0.0,
        [9., 0., 0., 0., 0., 0.],
        0.01,
    );
    assert_eq!(
        trace::reb_integrator_trace_switch_default(&r, &t, 0, 1),
        1,
        "body at 0.5 is inside 3*rs = 0.6"
    );

    let (r, t) = trace_pair(
        0.2,
        0.0,
        [0.7, 0., 0., 0., 0., 0.],
        0.0,
        [9., 0., 0., 0., 0., 0.],
        0.01,
    );
    assert_eq!(
        trace::reb_integrator_trace_switch_default(&r, &t, 0, 1),
        0,
        "body at 0.7 is outside 3*rs = 0.6 and is not approaching"
    );

    // With a point-like star (rs = 0) even a very close massless body
    // does not trigger the star-body switch.
    let (r, t) = trace_pair(
        0.0,
        0.0,
        [1e-6, 0., 0., 0., 0., 0.],
        0.0,
        [9., 0., 0., 0., 0., 0.],
        0.01,
    );
    assert_eq!(
        trace::reb_integrator_trace_switch_default(&r, &t, 0, 1),
        0,
        "a zero-radius star never triggers the body-body switch"
    );
}

#[test]
fn test_trace_switch_default_scales_with_r_crit_hill() {
    // dcritmax = r_crit_hill * Hill radius (= 0.1 here), so a companion
    // 0.5 away switches on only once r_crit_hill exceeds 5.
    for (rch, want) in [(3.0, 0), (4.9, 0), (5.1, 1), (20.0, 1)] {
        let (r, mut t) = trace_pair(
            0.0,
            3e-3,
            [1., 0., 0., 0., 0., 0.],
            0.0,
            [1.5, 0., 0., 0., 0., 0.],
            0.01,
        );
        t.r_crit_hill = rch;
        assert_eq!(
            trace::reb_integrator_trace_switch_default(&r, &t, 1, 2),
            want,
            "r_crit_hill = {}",
            rch
        );
    }
}

// =====================================================================
// TRACE: pericentre switching function
// =====================================================================

/// The pericentre switch evaluated for a circular orbit of radius `d`
/// about a star of mass `GM` (with G = 1).
fn peri_switch(d: f64, GM: f64, eta: f64, dt: f64) -> i32 {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_add(&mut r, mkp(GM, 0., 0., 0., 0., 0., 0., 0.));
    let v = (GM / d).sqrt();
    reb_simulation_add(&mut r, mkp(0.0, 0., d, 0., 0., 0., v, 0.));
    let mut t = trace_state_for(2);
    t.peri_crit_eta = eta;
    trace::reb_integrator_trace_switch_peri_default(&r, &t, 1)
}

#[test]
fn test_trace_peri_switch_threshold_is_the_inverse_mean_motion() {
    // For a circular orbit the Pham/Rein/Spiegel time scale of Eq. 16
    // reduces analytically to tau^2 = 2*|r''|^2/(|r'''|^2 + |r''||r''''|)
    // with |r''| = GM/d^2, |r'''| = GM v/d^3 and |r''''| = (GM)^2/d^5,
    // which for v^2 = GM/d gives tau^2 = d^3/(GM) = 1/n^2.  The switch
    // fires when dt > eta/n, i.e. dt > eta * P/(2 pi).
    let cases: [(f64, f64, f64); 5] = [
        (1.0, 1.0, 1.0),
        (4.0, 1.0, 1.0),
        (1.0, 1.0, 0.5),
        (2.5, 3.0, 1.7),
        (0.3, 0.7, 0.9),
    ];
    for (d, GM, eta) in cases {
        let want = eta * (d * d * d / GM).sqrt();
        assert_eq!(peri_switch(d, GM, eta, 0.0), 0, "dt = 0 never fires");
        // Bisect for the dt at which the flag flips.
        let mut lo = 0.0f64;
        let mut hi = want * 1e3;
        assert_eq!(peri_switch(d, GM, eta, hi), 1, "a huge dt must fire");
        for _ in 0..90 {
            let mid = 0.5 * (lo + hi);
            if peri_switch(d, GM, eta, mid) != 0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        assert!(
            (hi - want).abs() < 1e-10 * want,
            "peri switch threshold for d={} GM={} eta={}: {} vs 1/n = {}",
            d,
            GM,
            eta,
            hi,
            want
        );
    }
}

#[test]
fn test_trace_peri_switch_none_never_fires() {
    // The "none" prescription must not flag anything, whatever dt is.
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.0, 0., 1., 0., 0., 0., 1., 0.));
    let t = trace_state_for(2);
    for dt in [1e-6, 1.0, 1e3, 1e9] {
        r.dt = dt;
        assert_eq!(
            trace::reb_integrator_trace_switch_peri_none(&r, &t, 1),
            0,
            "peri_none with dt = {}",
            dt
        );
        // ... while the default prescription does fire for the large dt.
        let expected_default = if dt > 1.0 { 1 } else { 0 };
        assert_eq!(
            trace::reb_integrator_trace_switch_peri_default(&r, &t, 1),
            expected_default,
            "peri_default with dt = {}",
            dt
        );
    }
}

// =====================================================================
// TRACE: reversible pre-/post-timestep checks
// =====================================================================

/// Star + two planets on well-separated circular orbits, one on the
/// x axis and one on the y axis (so no accidental qv == 0 shortcuts).
fn trace_three_body(dt: f64, sep: f64) -> (reb_simulation, trace::reb_integrator_trace_state) {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1e-3, 0., 1., 0., 0., 0., 1., 0.));
    let v2 = (1.0f64 / 5.0).sqrt();
    if sep > 0.0 {
        // Second planet close to the first one.
        reb_simulation_add(&mut r, mkp(1e-3, 0., 1.0 + sep, 0., 0., 0., 1., 0.));
    } else {
        reb_simulation_add(&mut r, mkp(1e-3, 0., 0., 5., 0., -v2, 0., 0.));
    }
    let t = trace_state_for(3);
    (r, t)
}

#[test]
fn test_trace_pre_ts_check_finds_nothing_when_well_separated() {
    let (mut r, mut t) = trace_three_body(0.01, 0.0);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    assert_eq!(t.encounter_N, 1, "only the star is in the encounter map");
    assert_eq!(t.encounter_map[0], 1, "map slot 0 is always flagged");
    assert_eq!(t.encounter_map[1], 0, "planet 1 not flagged");
    assert_eq!(t.encounter_map[2], 0, "planet 2 not flagged");
    assert_eq!(t.current_C, 0, "no pericentre encounter");
    assert_eq!(t.tponly_encounter, 1, "no massive-massive encounter");
    for k in &t.current_Ks {
        assert_eq!(*k, 0, "no K_ij set");
    }
    assert_eq!(
        t.encounter_map_backup[..3],
        t.encounter_map[..3],
        "backup mirrors the encounter map"
    );
}

#[test]
fn test_trace_pre_ts_check_flags_a_close_pair() {
    // The two planets have Hill radius 1*(1e-3/3)^(1/3) = 0.0693, so
    // dcritmax = 3*0.0693 = 0.208; a separation of 0.1 is inside it.
    let (mut r, mut t) = trace_three_body(0.01, 0.1);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    let n = r.N;
    assert_eq!(t.encounter_N, 3, "star and both planets are flagged");
    assert_eq!(t.encounter_map[1], 1, "planet 1 flagged");
    assert_eq!(t.encounter_map[2], 1, "planet 2 flagged");
    assert_eq!(t.current_Ks[1 * n + 2], 1, "K_12 set");
    assert_eq!(t.current_Ks[0 * n + 1], 0, "K_01 not set (point-like star)");
    assert_eq!(t.current_Ks[0 * n + 2], 0, "K_02 not set (point-like star)");
    assert_eq!(t.tponly_encounter, 0, "two massive bodies encountered");
    assert_eq!(t.current_C, 0, "no pericentre encounter at this dt");
}

#[test]
fn test_trace_post_ts_check_reports_no_new_encounter_when_state_is_unchanged() {
    // Running the post check straight after the pre check on unchanged
    // coordinates must find exactly the same encounter set.
    for sep in [0.0, 0.1] {
        let (mut r, mut t) = trace_three_body(0.01, sep);
        trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
        let ks_before = t.current_Ks.clone();
        let n = trace::reb_integrator_trace_post_ts_check(&mut r, &mut t);
        assert_eq!(n, 0, "no new encounters for sep = {}", sep);
        assert_eq!(t.current_Ks, ks_before, "K_ij unchanged for sep = {}", sep);
    }
}

#[test]
fn test_trace_post_ts_check_detects_a_newly_close_pair() {
    // Pre-check with the planets far apart, then move them together:
    // the post check must reject the step by returning a nonzero count.
    let (mut r, mut t) = trace_three_body(0.01, 0.0);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    assert_eq!(t.encounter_N, 1, "nothing flagged before the step");
    // Move planet 2 next to planet 1.
    r.particles[2] = mkp(1e-3, 0., 1.08, 0.0, 0., 0., 1., 0.);
    let new_ce = trace::reb_integrator_trace_post_ts_check(&mut r, &mut t);
    assert_eq!(new_ce, 1, "post check must report the new encounter");
    let n = r.N;
    assert_eq!(t.current_Ks[1 * n + 2], 1, "K_12 set by the post check");
    assert_eq!(t.encounter_N, 3, "both planets flagged after the post check");
}

#[test]
fn test_trace_post_ts_check_restores_the_map_from_the_pre_check_backup() {
    // The post check starts from the pre-check encounter map, not from
    // whatever the encounter map happened to contain.
    let (mut r, mut t) = trace_three_body(0.01, 0.1);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    let backup = t.encounter_map_backup.clone();
    // Scribble over the live map the way reb_integrator_trace_bs_step does.
    t.encounter_map = vec![0, 1, 2];
    trace::reb_integrator_trace_post_ts_check(&mut r, &mut t);
    assert_eq!(
        t.encounter_map[..3],
        backup[..3],
        "post check reloads the map from encounter_map_backup"
    );
}

#[test]
fn test_trace_pre_ts_check_pericentre_full_bs_returns_immediately() {
    // Planet 1 is on a unit circular orbit, so the pericentre switch
    // fires for dt > 1.  With the FULL_BS prescription the routine sets
    // current_C and returns before touching the encounter map.
    let (mut r, mut t) = trace_three_body(2.0, 0.0);
    t.peri_mode = trace::REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS;
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    assert_eq!(t.current_C, 1, "pericentre encounter flagged");
    assert_eq!(t.encounter_N, 1, "early return leaves encounter_N at 1");
    for k in &t.current_Ks {
        assert_eq!(*k, 0, "early return leaves K_ij untouched");
    }
}

#[test]
fn test_trace_pre_ts_check_pericentre_partial_bs_flags_everything() {
    // With PARTIAL_BS a pericentre encounter puts the whole simulation
    // into the BS sub-step: encounter_N = N and every map slot set.
    let (mut r, mut t) = trace_three_body(2.0, 0.0);
    t.peri_mode = trace::REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS;
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    assert_eq!(t.current_C, 1, "pericentre encounter flagged");
    assert_eq!(t.encounter_N, r.N, "all particles in the encounter map");
    for i in 0..r.N {
        assert_eq!(t.encounter_map[i], 1, "map slot {} flagged", i);
    }
    assert_eq!(t.tponly_encounter, 0, "massive pericentre encounter");
}

#[test]
fn test_trace_pre_ts_check_with_peri_none_never_flags_a_pericentre() {
    let (mut r, mut t) = trace_three_body(2.0, 0.0);
    t.S_peri = Some(trace::reb_integrator_trace_switch_peri_none);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    assert_eq!(t.current_C, 0, "peri_none suppresses the pericentre flag");
}

#[test]
fn test_trace_pre_ts_check_honours_a_custom_switching_function() {
    fn always_on(
        _r: &reb_simulation,
        _t: &trace::reb_integrator_trace_state,
        _i: usize,
        _j: usize,
    ) -> i32 {
        1
    }
    let (mut r, mut t) = trace_three_body(0.01, 0.0);
    t.S = Some(always_on);
    t.S_peri = Some(trace::reb_integrator_trace_switch_peri_none);
    trace::reb_integrator_trace_pre_ts_check(&mut r, &mut t);
    let n = r.N;
    assert_eq!(t.encounter_N, 3, "custom switch flags everything");
    for i in 0..n {
        for j in (i + 1)..n {
            assert_eq!(t.current_Ks[i * n + j], 1, "K_{}{} set", i, j);
        }
    }
}

// =====================================================================
// TRACE: individual sub-steps
// =====================================================================

#[test]
fn test_trace_com_step_exact_value() {
    let mut t = trace_state_for(1);
    t.com_pos = reb_vec3d { x: 1.0, y: 1.0, z: 1.0 };
    t.com_vel = reb_vec3d { x: 2.0, y: 4.0, z: 8.0 };
    trace::reb_integrator_trace_com_step(&mut t, 0.5);
    assert_eq!(t.com_pos.x, 2.0, "com_step x");
    assert_eq!(t.com_pos.y, 3.0, "com_step y");
    assert_eq!(t.com_pos.z, 5.0, "com_step z");
    assert_eq!(t.com_vel.x, 2.0, "com_step leaves com_vel alone");
}

#[test]
fn test_trace_jump_step_exact_value_and_pericentre_suppression() {
    // m1 = 0.5 with vx = 2 gives P = 1; the shift is P*dt/m0 = 0.25.
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.25;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.5, 0., 1., 2., 4., 2., 0., 0.));
    let mut t = trace_state_for(2);
    let dt = r.dt;
    trace::reb_integrator_trace_jump_step(&mut r, &t, dt);
    assert_eq!(r.particles[1].x, 1.25, "trace jump step x");
    assert_eq!(r.particles[0].x, 0.0, "trace jump step never moves the star");

    // During a pericentre approach the jump step is skipped entirely.
    t.current_C = 1;
    trace::reb_integrator_trace_jump_step(&mut r, &t, dt);
    assert_eq!(
        r.particles[1].x, 1.25,
        "jump step is a no-op while current_C is set"
    );
}

#[test]
fn test_trace_interaction_step_skips_the_central_body_and_flagged_pairs() {
    // Only non-central pairs contribute.  With m1 = m2 = 1 unit apart
    // and G = 1 the mutual acceleration has magnitude 1, so after
    // dt = 0.5 we expect vx = +0.5 and -0.5.  Setting K_12 removes the
    // pair from the interaction Hamiltonian (BS handles it instead).
    for (k12, want) in [(0i32, 0.5f64), (1i32, 0.0f64)] {
        let mut r = quiet_sim();
        r.G = 1.0;
        reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
        reb_simulation_add(&mut r, mkp(1.0, 0., 1., 0., 0., 0., 0., 0.));
        reb_simulation_add(&mut r, mkp(1.0, 0., 2., 0., 0., 0., 0., 0.));
        let mut t = trace_state_for(3);
        t.current_Ks[1 * 3 + 2] = k12;
        trace::reb_integrator_trace_interaction_step(&mut r, &mut t, 0.5);
        assert_eq!(r.particles[1].vx, want, "v1 with K_12 = {}", k12);
        assert_eq!(r.particles[2].vx, -want, "v2 with K_12 = {}", k12);
        assert_eq!(r.particles[0].vx, 0.0, "the star is never kicked");
    }
}

#[test]
fn test_trace_kepler_acceleration_uses_only_flagged_pairs() {
    // Mirror of the interaction step: in the Kepler shell the star term
    // is always present and the pair term is present only when K_ij = 1.
    for (k12, extra) in [(1i32, -1.0f64), (0i32, 0.0f64)] {
        let mut r = quiet_sim();
        r.G = 1.0;
        reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
        reb_simulation_add(&mut r, mkp(1.0, 0., 1., 0., 0., 0., 0., 0.));
        reb_simulation_add(&mut r, mkp(1.0, 0., 2., 0., 0., 0., 0., 0.));
        let mut t = trace_state_for(3);
        t.current_Ks[1 * 3 + 2] = k12;
        t.encounter_map = vec![0, 1, 2];
        t.encounter_N = 3;
        t.encounter_N_active = 3;
        r.integrator = reb_integrator_state::trace(t);
        trace::reb_integrator_trace_calculate_acceleration_mode_kepler(&mut r);
        assert_eq!(r.particles[0].ax, 0.0, "star feels nothing in dh");
        assert_eq!(r.particles[1].ax, -1.0 - extra, "a1 with K_12 = {}", k12);
        assert_eq!(r.particles[2].ax, -0.25 + extra, "a2 with K_12 = {}", k12);
    }
}

#[test]
fn test_trace_whfast_step_advances_a_circular_orbit_by_one_period() {
    // a = 1, GM = 1 => P = 2 pi.  A Kepler step of one period is the
    // identity up to the solver's round-off.
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(0.0, 0., 1., 0., 0., 0., 1., 0.));
    trace::reb_integrator_trace_whfast_step(&mut r, 2.0 * PI);
    assert!(
        (r.particles[1].x - 1.0).abs() < 1e-12,
        "x after one period: {}",
        r.particles[1].x
    );
    assert!(
        r.particles[1].y.abs() < 1e-12,
        "y after one period: {}",
        r.particles[1].y
    );
    assert!(
        (r.particles[1].vy - 1.0).abs() < 1e-12,
        "vy after one period: {}",
        r.particles[1].vy
    );
    assert_eq!(r.particles[0].x, 0.0, "the Kepler step never moves the star");
}

#[test]
fn test_trace_whfast_step_is_reversible() {
    // Kepler evolution by +dt then -dt returns the original state.
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1e-3, 0., 1.0, 0., 0., 0., 0.7, 0.2));
    reb_simulation_add(&mut r, mkp(1e-4, 0., -0.4, 1.9, -0.3, -0.6, -0.1, 0.05));
    let before = r.particles.clone();
    trace::reb_integrator_trace_whfast_step(&mut r, 0.37);
    let moved = (r.particles[1].x - before[1].x).abs();
    assert!(moved > 1e-3, "the Kepler step actually moved something");
    trace::reb_integrator_trace_whfast_step(&mut r, -0.37);
    for i in 0..r.N {
        assert!(
            (r.particles[i].x - before[i].x).abs() < 1e-13,
            "reversed Kepler step, x of particle {}",
            i
        );
        assert!(
            (r.particles[i].vx - before[i].vx).abs() < 1e-13,
            "reversed Kepler step, vx of particle {}",
            i
        );
    }
}

#[test]
fn test_trace_inertial_to_dh_and_back_roundtrip() {
    let mut r = quiet_sim();
    r.G = 1.0;
    reb_simulation_add(&mut r, mkp(1.0, 0., -0.1, 0.4, 0.2, 0.02, -0.01, 0.005));
    reb_simulation_add(&mut r, mkp(3e-3, 0., 1.1, -0.6, 0.3, 0.4, 0.8, -0.05));
    reb_simulation_add(&mut r, mkp(7e-4, 0., -1.8, 0.9, -0.7, -0.3, -0.55, 0.02));
    let before = r.particles.clone();
    let mut t = trace_state_for(3);
    trace::reb_integrator_trace_inertial_to_dh(&mut r, &mut t);
    assert_eq!(r.particles[0].x, 0.0, "dh puts particle 0 at the origin");
    trace::reb_integrator_trace_dh_to_inertial(&mut r, &t);
    for i in 0..3 {
        for (got, want, name) in [
            (r.particles[i].x, before[i].x, "x"),
            (r.particles[i].y, before[i].y, "y"),
            (r.particles[i].z, before[i].z, "z"),
            (r.particles[i].vx, before[i].vx, "vx"),
            (r.particles[i].vy, before[i].vy, "vy"),
            (r.particles[i].vz, before[i].vz, "vz"),
        ] {
            assert!(
                (got - want).abs() < 1e-14,
                "trace dh roundtrip {} of particle {}: {} vs {}",
                name,
                i,
                got,
                want
            );
        }
    }
}

// =====================================================================
// TRACE: full integrations
// =====================================================================

#[test]
fn test_trace_single_particle_drifts_at_constant_velocity() {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = 0.1;
    reb_simulation_set_integrator(&mut r, "trace");
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., -0.25, 0.5, 0.0));
    reb_simulation_integrate(&mut r, 1.0);
    assert!((r.t - 1.0).abs() < 1e-14, "reached tmax");
    assert!(
        (r.particles[0].x + 0.25 * r.t).abs() < 1e-14,
        "lone star x: {}",
        r.particles[0].x
    );
    assert_eq!(r.particles[0].vx, -0.25, "lone star velocity unchanged");
}

#[test]
fn test_trace_conserves_energy_for_two_planets() {
    let mut r = two_planet_sim("trace", 2.0 * PI / 200.0, 1e-5);
    let drift = energy_drift(&mut r, 100.0);
    assert!(drift < 1e-8, "TRACE relative energy drift {}", drift);
}

#[test]
fn test_trace_conserves_angular_momentum_for_two_planets() {
    let mut r = two_planet_sim("trace", 2.0 * PI / 200.0, 1e-5);
    let l0 = reb_simulation_angular_momentum(&r);
    reb_simulation_integrate(&mut r, 100.0);
    let l1 = reb_simulation_angular_momentum(&r);
    let n0 = (l0.x * l0.x + l0.y * l0.y + l0.z * l0.z).sqrt();
    let dl = ((l1.x - l0.x).powi(2) + (l1.y - l0.y).powi(2) + (l1.z - l0.z).powi(2)).sqrt();
    assert!(
        dl / n0 < 1e-11,
        "TRACE relative angular momentum drift {}",
        dl / n0
    );
}

#[test]
fn test_trace_forced_encounter_branch_conserves_energy() {
    // A huge r_crit_hill drives every pair through the BS encounter
    // sub-step, exercising the encounter map, the K_ij bookkeeping and
    // the pre/post reversibility check on every timestep.
    let mut r = two_planet_sim("trace", 2.0 * PI / 200.0, 1e-5);
    if let reb_integrator_state::trace(ref mut t) = r.integrator {
        t.r_crit_hill = 5000.0;
    }
    let e0 = reb_simulation_energy(&r);
    reb_simulation_integrate(&mut r, 20.0);
    let e1 = reb_simulation_energy(&r);
    let drift = ((e1 - e0) / e0).abs();
    if let reb_integrator_state::trace(ref t) = r.integrator {
        assert_eq!(t.encounter_N, 3, "all particles routed through BS");
        assert_eq!(t.current_Ks[1 * 3 + 2], 1, "K_12 set every step");
    } else {
        panic!("integrator is no longer trace");
    }
    assert!(drift < 1e-9, "TRACE forced-encounter energy drift {}", drift);
}

#[test]
fn test_trace_agrees_with_ias15_on_a_two_body_problem() {
    let mut a = quiet_sim();
    a.G = 1.0;
    a.dt = 2.0 * PI / 2000.0;
    reb_simulation_set_integrator(&mut a, "trace");
    reb_simulation_add(&mut a, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    let star = a.particles[0];
    let p = reb_particle_from_orbit(1.0, star, 1e-6, 1.0, 0.3, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut a, p);
    reb_simulation_move_to_com(&mut a);

    let mut b = quiet_sim();
    b.G = 1.0;
    b.dt = 2.0 * PI / 2000.0;
    reb_simulation_set_integrator(&mut b, "ias15");
    reb_simulation_add(&mut b, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut b, p);
    reb_simulation_move_to_com(&mut b);

    reb_simulation_integrate(&mut a, 4.0 * PI);
    reb_simulation_integrate(&mut b, 4.0 * PI);
    let dx = a.particles[1].x - b.particles[1].x;
    let dy = a.particles[1].y - b.particles[1].y;
    let dz = a.particles[1].z - b.particles[1].z;
    let d = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(d < 1e-6, "TRACE vs IAS15 separation after 2 orbits: {}", d);
}

#[test]
fn test_trace_is_deterministic() {
    let mut a = two_planet_sim("trace", 0.03, 1e-4);
    let mut b = two_planet_sim("trace", 0.03, 1e-4);
    reb_simulation_integrate(&mut a, 12.0);
    reb_simulation_integrate(&mut b, 12.0);
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical TRACE runs must agree bit for bit"
    );
}

#[test]
fn test_trace_full_ias15_peri_mode_matches_full_bs_peri_mode() {
    // Both FULL prescriptions integrate the whole pericentre passage
    // with an adaptive high-order scheme, so on a smooth eccentric orbit
    // they must agree far more closely than the WH splitting error.
    let make = |mode: i32| {
        let mut r = quiet_sim();
        r.G = 1.0;
        r.dt = 2.0 * PI / 30.0; // large enough to trigger the pericentre switch
        reb_simulation_set_integrator(&mut r, "trace");
        if let reb_integrator_state::trace(ref mut t) = r.integrator {
            t.peri_mode = mode;
        }
        reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
        let star = r.particles[0];
        let p = reb_particle_from_orbit(1.0, star, 1e-6, 1.0, 0.7, 0.0, 0.0, 0.0, 0.0);
        reb_simulation_add(&mut r, p);
        reb_simulation_move_to_com(&mut r);
        r
    };
    let mut a = make(trace::REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS);
    let mut b = make(trace::REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15);
    reb_simulation_integrate(&mut a, 4.0 * PI);
    reb_simulation_integrate(&mut b, 4.0 * PI);
    let dx = a.particles[1].x - b.particles[1].x;
    let dy = a.particles[1].y - b.particles[1].y;
    let d = (dx * dx + dy * dy).sqrt();
    assert!(d < 1e-5, "FULL_BS vs FULL_IAS15 separation: {}", d);
}

// =====================================================================
// JANUS: integer grid arithmetic
// =====================================================================

fn janus_sim(order: u32, dt: f64) -> reb_simulation {
    let mut r = quiet_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_set_integrator(&mut r, "janus");
    if let reb_integrator_state::janus(ref mut j) = r.integrator {
        j.order = order;
    }
    reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1e-3, 0., 1.0, 0., 0., 0., 1.0, 0.));
    reb_simulation_add(&mut r, mkp(1e-3, 0., 0., 2.0, 0., -0.7, 0., 0.02));
    reb_simulation_move_to_com(&mut r);
    r
}

fn janus_ints(r: &reb_simulation) -> Vec<janus::reb_particle_int> {
    match &r.integrator {
        reb_integrator_state::janus(j) => j.p_int.clone(),
        _ => panic!("not a janus simulation"),
    }
}

/// integrator_janus.c `to_int`, reproduced here so the expected initial
/// grid coordinates are derived independently of the library.
fn to_int_ref(p: &reb_particle, scale_pos: f64, scale_vel: f64) -> janus::reb_particle_int {
    janus::reb_particle_int {
        x: (p.x / scale_pos) as i64,
        y: (p.y / scale_pos) as i64,
        z: (p.z / scale_pos) as i64,
        vx: (p.vx / scale_vel) as i64,
        vy: (p.vy / scale_vel) as i64,
        vz: (p.vz / scale_vel) as i64,
    }
}

#[test]
fn test_janus_is_bitwise_time_reversible() {
    // The JANUS drift and kick both add trunc(x) to an integer, and C's
    // double->int64 conversion truncates towards zero, so trunc(-x) is
    // exactly -trunc(x).  Together with the palindromic gamma sequence
    // this makes a step with -dt the exact inverse of a step with +dt.
    for order in [2u32, 4, 6, 8, 10] {
        let mut r = janus_sim(order, 0.01);
        let expected: Vec<janus::reb_particle_int> = {
            let (sp, sv) = match &r.integrator {
                reb_integrator_state::janus(j) => (j.scale_pos, j.scale_vel),
                _ => unreachable!(),
            };
            r.particles.iter().map(|p| to_int_ref(p, sp, sv)).collect()
        };
        reb_simulation_steps(&mut r, 40);
        let mid = janus_ints(&r);
        assert!(
            mid[1].x != expected[1].x,
            "order {}: forward integration must move the grid coordinates",
            order
        );
        r.dt = -r.dt;
        reb_simulation_steps(&mut r, 40);
        let back = janus_ints(&r);
        for i in 0..r.N {
            assert_eq!(back[i].x, expected[i].x, "order {} reversed x[{}]", order, i);
            assert_eq!(back[i].y, expected[i].y, "order {} reversed y[{}]", order, i);
            assert_eq!(back[i].z, expected[i].z, "order {} reversed z[{}]", order, i);
            assert_eq!(
                back[i].vx, expected[i].vx,
                "order {} reversed vx[{}]",
                order, i
            );
            assert_eq!(
                back[i].vy, expected[i].vy,
                "order {} reversed vy[{}]",
                order, i
            );
            assert_eq!(
                back[i].vz, expected[i].vz,
                "order {} reversed vz[{}]",
                order, i
            );
        }
    }
}

#[test]
fn test_janus_particles_live_exactly_on_the_integer_grid() {
    // After synchronisation the floating point state is exactly
    // p_int * scale, so every coordinate is an exact grid point.
    let mut r = janus_sim(6, 0.01);
    reb_simulation_steps(&mut r, 17);
    let (sp, sv) = match &r.integrator {
        reb_integrator_state::janus(j) => (j.scale_pos, j.scale_vel),
        _ => unreachable!(),
    };
    let ints = janus_ints(&r);
    for i in 0..r.N {
        assert_eq!(
            r.particles[i].x.to_bits(),
            ((ints[i].x as f64) * sp).to_bits(),
            "particle {} x is on the position grid",
            i
        );
        assert_eq!(
            r.particles[i].vy.to_bits(),
            ((ints[i].vy as f64) * sv).to_bits(),
            "particle {} vy is on the velocity grid",
            i
        );
    }
}

#[test]
fn test_janus_to_int_of_to_double_loses_at_most_a_few_grid_units() {
    // n -> n*scale -> n/scale is NOT the identity: fl(n*s) carries a
    // relative error up to u = 2^-53, fl(d/s) another, and the C cast
    // truncates towards zero, so the reconstruction can miss by up to
    // 2*u*|n| + 1 grid units.  With |n| ~ 1e16 that is about 5.  This is
    // precisely why JANUS carries p_int across timesteps instead of
    // re-deriving it, and why the recalculate flag has to be explicit.
    let mut r = janus_sim(6, 0.01);
    reb_simulation_steps(&mut r, 5);
    let (sp, sv) = match &r.integrator {
        reb_integrator_state::janus(j) => (j.scale_pos, j.scale_vel),
        _ => unreachable!(),
    };
    let ints = janus_ints(&r);
    for i in 0..r.N {
        let again = to_int_ref(&r.particles[i], sp, sv);
        for (a, b, name) in [
            (again.x, ints[i].x, "x"),
            (again.y, ints[i].y, "y"),
            (again.z, ints[i].z, "z"),
            (again.vx, ints[i].vx, "vx"),
            (again.vy, ints[i].vy, "vy"),
            (again.vz, ints[i].vz, "vz"),
        ] {
            let bound = 2 + (2.0f64 * (b as f64).abs() * f64::EPSILON) as i64;
            assert!(
                (a - b).abs() <= bound,
                "to_int(to_double({})) for {}[{}] drifted by {} (bound {})",
                b,
                name,
                i,
                a - b,
                bound
            );
        }
    }
}

#[test]
fn test_janus_conserves_energy() {
    // JANUS with order 6 on a well separated system: the discretisation
    // error dominates the integer round-off at these scales.
    let mut r = janus_sim(6, 0.01);
    let drift = energy_drift(&mut r, 20.0);
    assert!(drift < 1e-9, "JANUS order-6 energy drift {}", drift);
}

#[test]
fn test_janus_higher_order_is_more_accurate() {
    // Compare against IAS15 on the same initial conditions.
    let reference = {
        let mut r = janus_sim(6, 0.05);
        reb_simulation_set_integrator(&mut r, "ias15");
        reb_simulation_integrate(&mut r, 20.0);
        (r.particles[1].x, r.particles[1].y)
    };
    let err = |order: u32| {
        let mut r = janus_sim(order, 0.05);
        reb_simulation_integrate(&mut r, 20.0);
        let dx = r.particles[1].x - reference.0;
        let dy = r.particles[1].y - reference.1;
        (dx * dx + dy * dy).sqrt()
    };
    let e2 = err(2);
    let e4 = err(4);
    let e6 = err(6);
    assert!(e4 < e2 * 1e-2, "JANUS order 4 ({}) beats order 2 ({})", e4, e2);
    assert!(e6 < e4, "JANUS order 6 ({}) beats order 4 ({})", e6, e4);
}

#[test]
fn test_janus_recalculate_flag_resets_the_grid_from_the_floats() {
    // Modifying the particles is only picked up when the recalculate
    // flag is set; that is the documented JANUS contract.
    let mut r = janus_sim(6, 0.01);
    reb_simulation_steps(&mut r, 3);
    r.particles[1].x = 1.5;
    if let reb_integrator_state::janus(ref mut j) = r.integrator {
        j.recalculate_integer_coordinates_this_timestep = 1;
    }
    let (sp, sv) = match &r.integrator {
        reb_integrator_state::janus(j) => (j.scale_pos, j.scale_vel),
        _ => unreachable!(),
    };
    let want = to_int_ref(&r.particles[1], sp, sv);
    // Take a step with dt = 0 so the grid is re-seeded but not evolved.
    r.dt = 0.0;
    reb_simulation_steps(&mut r, 1);
    let ints = janus_ints(&r);
    assert_eq!(ints[1].x, want.x, "grid re-seeded from the modified float");
    if let reb_integrator_state::janus(ref j) = r.integrator {
        assert_eq!(
            j.recalculate_integer_coordinates_this_timestep, 0,
            "the recalculate flag is cleared after use"
        );
    }
}

#[test]
fn test_janus_is_deterministic() {
    let mut a = janus_sim(6, 0.01);
    let mut b = janus_sim(6, 0.01);
    reb_simulation_steps(&mut a, 200);
    reb_simulation_steps(&mut b, 200);
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical JANUS runs must agree bit for bit"
    );
    assert_eq!(
        janus_ints(&a)
            .iter()
            .map(|p| p.x)
            .collect::<Vec<i64>>(),
        janus_ints(&b)
            .iter()
            .map(|p| p.x)
            .collect::<Vec<i64>>(),
        "grid coordinates must agree too"
    );
}

// =====================================================================
// EOS: embedded operator splitting schemes
// =====================================================================

const EOS_TYPES: [(&str, i32); 9] = [
    ("LF", eos::REB_INTEGRATOR_EOS_TYPE_LF),
    ("LF4", eos::REB_INTEGRATOR_EOS_TYPE_LF4),
    ("LF6", eos::REB_INTEGRATOR_EOS_TYPE_LF6),
    ("LF8", eos::REB_INTEGRATOR_EOS_TYPE_LF8),
    ("LF4_2", eos::REB_INTEGRATOR_EOS_TYPE_LF4_2),
    ("LF8_6_4", eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4),
    ("PLF7_6_4", eos::REB_INTEGRATOR_EOS_TYPE_PLF7_6_4),
    ("PMLF4", eos::REB_INTEGRATOR_EOS_TYPE_PMLF4),
    ("PMLF6", eos::REB_INTEGRATOR_EOS_TYPE_PMLF6),
];

fn eos_sim(phi0: i32, phi1: i32, n: u32, dt: f64, safe_mode: u32) -> reb_simulation {
    let mut r = two_planet_sim("eos", dt, 1e-4);
    if let reb_integrator_state::eos(ref mut e) = r.integrator {
        e.phi0 = phi0;
        e.phi1 = phi1;
        e.n = n;
        e.safe_mode = safe_mode;
    }
    r
}

#[test]
fn test_eos_every_splitting_scheme_conserves_energy() {
    // All nine schemes are symplectic (or symplectic with pre/post
    // processing), so none of them may drift secularly.  At dt = P1/100
    // with n = 2 inner sub-steps the bounded energy error is set by each
    // scheme's own order: ~(dt/(n P))^p, i.e. ~1e-5 for the second-order
    // shells (LF, LF4_2), ~1e-8..1e-9 for the fourth-order ones and the
    // double-precision floor (~1e-13) for sixth order and above.  Each
    // bound below is that estimate rounded up; a scheme that silently
    // lost an order of accuracy would break through it.
    let bounds: [(&str, f64); 9] = [
        ("LF", 1e-5),
        ("LF4", 1e-7),
        ("LF6", 1e-12),
        ("LF8", 1e-12),
        ("LF4_2", 1e-5),
        ("LF8_6_4", 1e-10),
        ("PLF7_6_4", 1e-8),
        ("PMLF4", 1e-8),
        ("PMLF6", 1e-12),
    ];
    for (i, (name, ty)) in EOS_TYPES.into_iter().enumerate() {
        assert_eq!(bounds[i].0, name, "bound table lines up with EOS_TYPES");
        let mut r = eos_sim(ty, ty, 2, 2.0 * PI / 100.0, 1);
        let drift = energy_drift(&mut r, 40.0);
        assert!(
            drift < bounds[i].1,
            "EOS scheme {} relative energy drift {} exceeds {}",
            name,
            drift,
            bounds[i].1
        );
    }
}

#[test]
fn test_eos_higher_order_schemes_are_more_accurate() {
    // Compare the final position of the inner planet with an IAS15
    // reference.  LF is second order, LF4 fourth, LF6 sixth, so at a
    // fixed dt the errors must fall steeply.
    let dt = 2.0 * PI / 60.0;
    let tmax = 20.0;
    let reference = {
        let mut r = two_planet_sim("ias15", dt, 1e-4);
        reb_simulation_integrate(&mut r, tmax);
        (r.particles[1].x, r.particles[1].y, r.particles[1].z)
    };
    let err = |ty: i32| {
        let mut r = eos_sim(ty, ty, 2, dt, 1);
        reb_simulation_integrate(&mut r, tmax);
        let dx = r.particles[1].x - reference.0;
        let dy = r.particles[1].y - reference.1;
        let dz = r.particles[1].z - reference.2;
        (dx * dx + dy * dy + dz * dz).sqrt()
    };
    let e_lf = err(eos::REB_INTEGRATOR_EOS_TYPE_LF);
    let e_lf4 = err(eos::REB_INTEGRATOR_EOS_TYPE_LF4);
    let e_lf6 = err(eos::REB_INTEGRATOR_EOS_TYPE_LF6);
    assert!(e_lf > 0.0, "the LF error must be nonzero");
    assert!(
        e_lf4 < 0.1 * e_lf,
        "EOS LF4 error {} not clearly below LF error {}",
        e_lf4,
        e_lf
    );
    assert!(
        e_lf6 < e_lf4,
        "EOS LF6 error {} not below LF4 error {}",
        e_lf6,
        e_lf4
    );
}

#[test]
fn test_eos_lf_converges_at_second_order() {
    // Halving dt must reduce the error of the second-order LF splitting
    // by roughly a factor of four.
    let tmax = 8.0;
    let dt = 2.0 * PI / 40.0;
    let reference = {
        let mut r = two_planet_sim("ias15", dt, 1e-4);
        reb_simulation_integrate(&mut r, tmax);
        (r.particles[1].x, r.particles[1].y, r.particles[1].z)
    };
    let err = |h: f64| {
        let mut r = eos_sim(
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            1,
            h,
            1,
        );
        reb_simulation_integrate(&mut r, tmax);
        let dx = r.particles[1].x - reference.0;
        let dy = r.particles[1].y - reference.1;
        let dz = r.particles[1].z - reference.2;
        (dx * dx + dy * dy + dz * dz).sqrt()
    };
    let e1 = err(dt);
    let e2 = err(dt / 2.0);
    let ratio = e1 / e2;
    assert!(
        ratio > 3.2 && ratio < 5.0,
        "EOS LF convergence ratio {} is not ~4 (errors {} and {})",
        ratio,
        e1,
        e2
    );
}

#[test]
fn test_eos_more_inner_substeps_improve_accuracy() {
    // The inner shell resolves the Kepler part with n sub-steps; raising
    // n must reduce the error of the inner splitting.
    let dt = 2.0 * PI / 20.0;
    let tmax = 12.0;
    let reference = {
        let mut r = two_planet_sim("ias15", dt, 1e-4);
        reb_simulation_integrate(&mut r, tmax);
        (r.particles[1].x, r.particles[1].y)
    };
    let err = |n: u32| {
        let mut r = eos_sim(
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            n,
            dt,
            1,
        );
        reb_simulation_integrate(&mut r, tmax);
        let dx = r.particles[1].x - reference.0;
        let dy = r.particles[1].y - reference.1;
        (dx * dx + dy * dy).sqrt()
    };
    let e1 = err(1);
    let e8 = err(8);
    assert!(e1 > 1e-12, "the n=1 run has a resolvable error: {}", e1);
    assert!(
        e8 < e1,
        "EOS with n=8 ({}) must beat n=1 ({})",
        e8,
        e1
    );
}

#[test]
fn test_eos_safe_mode_one_leaves_the_simulation_synchronized() {
    let mut r = eos_sim(
        eos::REB_INTEGRATOR_EOS_TYPE_LF,
        eos::REB_INTEGRATOR_EOS_TYPE_LF,
        2,
        0.01,
        1,
    );
    reb_simulation_step(&mut r);
    assert_eq!(r.is_synchronized, 1, "safe_mode = 1 synchronises each step");

    let mut r0 = eos_sim(
        eos::REB_INTEGRATOR_EOS_TYPE_LF,
        eos::REB_INTEGRATOR_EOS_TYPE_LF,
        2,
        0.01,
        0,
    );
    reb_simulation_step(&mut r0);
    assert_eq!(r0.is_synchronized, 0, "safe_mode = 0 defers the drift");
    reb_simulation_synchronize(&mut r0);
    assert_eq!(r0.is_synchronized, 1, "explicit synchronise catches up");
}

#[test]
fn test_eos_synchronize_is_idempotent() {
    // A second synchronise on an already synchronised state must be a
    // bit-for-bit no-op.
    let mut r = eos_sim(
        eos::REB_INTEGRATOR_EOS_TYPE_PMLF6,
        eos::REB_INTEGRATOR_EOS_TYPE_LF4,
        2,
        0.01,
        0,
    );
    reb_simulation_steps(&mut r, 7);
    reb_simulation_synchronize(&mut r);
    let a = state_bits(&r);
    reb_simulation_synchronize(&mut r);
    let b = state_bits(&r);
    assert_eq!(a, b, "repeated EOS synchronise must not move anything");
}

#[test]
fn test_eos_unsafe_mode_stays_within_the_schemes_own_truncation_error() {
    // safe_mode = 0 merges the drift at the end of one step with the
    // drift at the start of the next (dtfac = 2).  That is a different
    // composition, not a round-off-level change.  drift_shell0(h) is
    // itself an inner leapfrog with n sub-steps, so its error over a
    // drifted time h is ~C*h*(h/n)^2; replacing two calls with h = dt/2
    // (total C*dt^3/(4n^2)) by one call with h = dt (C*dt^3/n^2) inflates
    // exactly that piece by a factor of four.  So the expected penalty
    // for dropping safe mode is ~4x, and FACTOR is set to twice that.
    const FACTOR: f64 = 8.0;
    let tmax = 10.0;
    let dt = 2.0 * PI / 100.0;
    let reference = {
        let mut r = two_planet_sim("ias15", dt, 1e-4);
        reb_simulation_integrate(&mut r, tmax);
        (r.particles[1].x, r.particles[1].y)
    };
    let run = |safe_mode: u32| {
        let mut r = eos_sim(
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            eos::REB_INTEGRATOR_EOS_TYPE_LF,
            2,
            dt,
            safe_mode,
        );
        let e0 = reb_simulation_energy(&r);
        reb_simulation_integrate(&mut r, tmax);
        let e1 = reb_simulation_energy(&r);
        (r.particles[1].x, r.particles[1].y, ((e1 - e0) / e0).abs(), e1, e0)
    };
    let a = run(1);
    let b = run(0);
    assert_eq!(a.4, b.4, "both runs start from the same energy");

    let err = |p: &(f64, f64, f64, f64, f64)| {
        ((p.0 - reference.0).powi(2) + (p.1 - reference.1).powi(2)).sqrt()
    };
    let ea = err(&a);
    let eb = err(&b);
    assert!(ea > 0.0, "the safe-mode LF run has a nonzero truncation error");
    let d = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
    assert!(
        d <= FACTOR * ea.max(eb),
        "safe_mode 0 vs 1 differ by {}, more than the truncation errors {} / {}",
        d,
        ea,
        eb
    );
    assert!(
        eb <= FACTOR * ea,
        "safe_mode 0 error {} is far worse than safe_mode 1 error {}",
        eb,
        ea
    );
    // Same argument for the energy: dropping safe mode must not degrade
    // the (bounded) energy error, and the two runs' final energies must
    // agree to within that error.
    assert!(a.2 > 0.0, "the safe-mode LF run has a nonzero energy error");
    assert!(
        b.2 <= FACTOR * a.2,
        "safe_mode 0 energy error {} is far worse than safe_mode 1 error {}",
        b.2,
        a.2
    );
    let de = ((a.3 - b.3) / a.4).abs();
    assert!(
        de <= FACTOR * a.2.max(b.2),
        "EOS safe_mode 0 vs 1 energy mismatch {} exceeds their own errors {} / {}",
        de,
        a.2,
        b.2
    );
}

#[test]
fn test_eos_is_deterministic() {
    let mut a = eos_sim(
        eos::REB_INTEGRATOR_EOS_TYPE_PLF7_6_4,
        eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4,
        3,
        0.02,
        1,
    );
    let mut b = eos_sim(
        eos::REB_INTEGRATOR_EOS_TYPE_PLF7_6_4,
        eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4,
        3,
        0.02,
        1,
    );
    reb_simulation_integrate(&mut a, 9.0);
    reb_simulation_integrate(&mut b, 9.0);
    assert_eq!(
        state_bits(&a),
        state_bits(&b),
        "two identical EOS runs must agree bit for bit"
    );
}

#[test]
fn test_eos_conserves_angular_momentum() {
    // Every EOS operator is a drift or a kick derived from a central
    // potential, so the total angular momentum is preserved up to
    // round-off for all nine schemes.
    for (name, ty) in EOS_TYPES {
        let mut r = eos_sim(ty, ty, 2, 2.0 * PI / 100.0, 1);
        let l0 = reb_simulation_angular_momentum(&r);
        reb_simulation_integrate(&mut r, 20.0);
        let l1 = reb_simulation_angular_momentum(&r);
        let n0 = (l0.x * l0.x + l0.y * l0.y + l0.z * l0.z).sqrt();
        let dl = ((l1.x - l0.x).powi(2) + (l1.y - l0.y).powi(2) + (l1.z - l0.z).powi(2)).sqrt();
        assert!(
            dl / n0 < 1e-11,
            "EOS scheme {} angular momentum drift {}",
            name,
            dl / n0
        );
    }
}

// =====================================================================
// Cross-integrator agreement across the whole hybrid group
// =====================================================================

#[test]
fn test_all_hybrid_integrators_agree_on_the_same_kepler_orbit() {
    // One eccentric two-body problem, four integrators, one reference.
    let build = |name: &str, dt: f64| {
        let mut r = quiet_sim();
        r.G = 1.0;
        r.dt = dt;
        reb_simulation_set_integrator(&mut r, name);
        reb_simulation_add(&mut r, mkp(1.0, 0., 0., 0., 0., 0., 0., 0.));
        let star = r.particles[0];
        let p = reb_particle_from_orbit(1.0, star, 1e-6, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0);
        reb_simulation_add(&mut r, p);
        reb_simulation_move_to_com(&mut r);
        r
    };
    let tmax = 6.0 * PI;
    let dt = 2.0 * PI / 3000.0;
    let mut refsim = build("ias15", dt);
    reb_simulation_integrate(&mut refsim, tmax);
    let rx = refsim.particles[1].x;
    let ry = refsim.particles[1].y;

    for name in ["mercurius", "trace", "janus", "eos"] {
        let mut r = build(name, dt);
        reb_simulation_integrate(&mut r, tmax);
        let dx = r.particles[1].x - rx;
        let dy = r.particles[1].y - ry;
        let d = (dx * dx + dy * dy).sqrt();
        assert!(
            d < 1e-5,
            "{} disagrees with IAS15 by {} after three orbits",
            name,
            d
        );
    }
}
