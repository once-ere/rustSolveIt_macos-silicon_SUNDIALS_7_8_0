//! Integration tests for the forces_dissipative group of reboundx_rs.
//! Part of reboundx_rs, GPL-3.0-or-later.
#![allow(non_snake_case)]
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::redundant_closure)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::too_many_arguments)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. Each waiver below carries its own reason; the same
// list and the rationale are in README.md under "Building and testing".
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::assign_op_pattern)]
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
#![allow(clippy::ptr_arg)]
#![allow(clippy::seek_from_current)]
#![allow(clippy::drop_non_drop)]
#![allow(clippy::approx_constant)]
#![allow(clippy::useless_vec)]
#![allow(clippy::type_complexity)]
use rebound_rs::*;
use reboundx_rs::*;

use std::f64::consts::PI;

const TWO_PI: f64 = 2.0 * PI;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Star at the origin plus one body on the given osculating orbit.
///
/// When `m_p == 0` the body is made a test particle (`N_active = 1`), which
/// makes the REBOUNDx back-reaction mass ratio identically zero and leaves the
/// star exactly at rest at the origin — so `reb_orbit_from_particle(G, p,
/// particles[0])` is the same two-body element set the effects themselves use
/// (their Jacobi reference is the centre of mass with the body removed, i.e.
/// the star).
fn star_and_body(
    m_star: f64,
    r_star: f64,
    m_p: f64,
    r_p: f64,
    a: f64,
    e: f64,
    inc: f64,
    Omega: f64,
    omega: f64,
    f: f64,
    integrator: &str,
    dt: f64,
) -> reb_simulation {
    let mut sim = reb_simulation_create();
    reb_simulation_set_integrator(&mut sim, integrator);
    sim.dt = dt;

    let mut star = reb_particle::default();
    star.m = m_star;
    star.r = r_star;
    reb_simulation_add(&mut sim, star);

    let mut p = reb_particle_from_orbit(sim.G, sim.particles[0], m_p, a, e, inc, Omega, omega, f);
    p.r = r_p;
    reb_simulation_add(&mut sim, p);

    if m_p == 0.0 {
        sim.N_active = 1;
    }
    sim
}

/// Osculating orbit of particle `i` relative to particle 0.
fn orb(sim: &reb_simulation, i: usize) -> reb_orbit {
    reb_orbit_from_particle(sim.G, sim.particles[i], sim.particles[0])
}

/// Keplerian period of particle `i` about particle 0.
fn period(sim: &reb_simulation, i: usize) -> f64 {
    let o = orb(sim, i);
    TWO_PI * (o.a * o.a * o.a / (sim.G * (sim.particles[0].m + sim.particles[i].m))).sqrt()
}

/// Distance of particle `i` from particle 0.
fn dist(sim: &reb_simulation, i: usize) -> f64 {
    let d = sim.particles[i];
    let s = sim.particles[0];
    ((d.x - s.x).powi(2) + (d.y - s.y).powi(2) + (d.z - s.z).powi(2)).sqrt()
}

/// Magnitude of the specific relative angular momentum |r x v| of particle `i`
/// about particle 0. A purely radial perturbation exerts no torque, so this is
/// an exactly conserved quantity for those effects.
fn hmag(sim: &reb_simulation, i: usize) -> f64 {
    let p = sim.particles[i];
    let s = sim.particles[0];
    let (x, y, z) = (p.x - s.x, p.y - s.y, p.z - s.z);
    let (vx, vy, vz) = (p.vx - s.vx, p.vy - s.vy, p.vz - s.vz);
    let hx = y * vz - z * vy;
    let hy = z * vx - x * vz;
    let hz = x * vy - y * vx;
    (hx * hx + hy * hy + hz * hz).sqrt()
}

/// Integrate to `t_end` in `n` equal chunks, sampling `f` after each chunk.
fn sample<F>(sim: &mut reb_simulation, t_end: f64, n: usize, mut f: F) -> Vec<f64>
where
    F: FnMut(&reb_simulation) -> f64,
{
    let t0 = sim.t;
    let mut out = Vec::with_capacity(n);
    for k in 1..=n {
        reb_simulation_integrate(sim, t0 + (t_end - t0) * (k as f64) / (n as f64));
        out.push(f(sim));
    }
    out
}

/// Mean of `f` over one orbital period, sampled at `nsub` equally spaced
/// phases. Averaging over a full period removes the periodic wobble that a
/// conservative perturbation puts into the osculating elements, leaving only
/// the secular drift.
fn mean_over_period<F>(sim: &mut reb_simulation, nsub: usize, mut f: F) -> f64
where
    F: FnMut(&reb_simulation) -> f64,
{
    let p = period(sim, 1);
    let t0 = sim.t;
    let mut s = 0.0;
    for k in 1..=nsub {
        reb_simulation_integrate(sim, t0 + p * (k as f64) / (nsub as f64));
        s += f(sim);
    }
    s / (nsub as f64)
}

fn assert_monotone(label: &str, vals: &[f64], want_increase: bool) {
    assert!(
        vals.len() >= 2,
        "{}: need at least two samples, got {}",
        label,
        vals.len()
    );
    for k in 1..vals.len() {
        let ok = if want_increase {
            vals[k] > vals[k - 1]
        } else {
            vals[k] < vals[k - 1]
        };
        assert!(
            ok,
            "{}: sample[{}] = {:.17e} is not {} sample[{}] = {:.17e} (series = {:?})",
            label,
            k,
            vals[k],
            if want_increase {
                "greater than"
            } else {
                "less than"
            },
            k - 1,
            vals[k - 1],
            vals
        );
    }
}

/// Every particle coordinate as raw IEEE-754 bits, plus the simulation time.
fn state_bits(sim: &reb_simulation) -> Vec<u64> {
    let mut v = Vec::with_capacity(6 * sim.N + 1);
    for i in 0..sim.N {
        let p = sim.particles[i];
        for x in [p.x, p.y, p.z, p.vx, p.vy, p.vz] {
            v.push(x.to_bits());
        }
    }
    v.push(sim.t.to_bits());
    v
}

fn set_p(sim: &mut reb_simulation, sel: rebx_ap, name: &str, v: f64) {
    let rebx = rebx_extras_mut(sim).expect("REBOUNDx attached");
    rebx_set_param_double(rebx, sel, name, v);
}

// ---------------------------------------------------------------------------
// modify_orbits_forces
// ---------------------------------------------------------------------------

/// `modify_orbits_forces` with only `tau_a` set applies a = dv * (1/tau_a) / 2.
/// For a body on a circular orbit around a fixed primary the specific energy
/// obeys dE/dt = v.a = v^2/(2 tau_a); with v^2 = mu/a and E = -mu/(2a) this is
/// exactly da/dt = a/tau_a, so a(t) = a0 exp(t/tau_a) with no orbit-averaging
/// approximation at all. tau_a < 0 therefore shrinks a, monotonically.
#[test]
fn mof_negative_tau_a_shrinks_a_as_a_pure_exponential() {
    let tau_a = -1.0e4;
    let a0_nominal = 1.0;
    let mut sim = star_and_body(
        1.0, 0.0, 0.0, 0.0, a0_nominal, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2,
    );
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(modify_orbits_forces) must report success"
    );
    set_p(&mut sim, rebx_ap::particle(1), "tau_a", tau_a);

    let a0 = orb(&sim, 1).a;
    let t_end = 200.0 * TWO_PI;
    let series = sample(&mut sim, t_end, 10, |s| orb(s, 1).a);

    assert_monotone("modify_orbits_forces tau_a<0: a", &series, false);

    let a_end = *series.last().unwrap();
    let predicted = a0 * (sim.t / tau_a).exp();
    let rel = (a_end - predicted).abs() / predicted;
    assert!(
        rel < 1.0e-6,
        "modify_orbits_forces tau_a<0: a(t) = {:.17e} but a0*exp(t/tau_a) = {:.17e} \
         (a0 = {:.17e}, t = {:.17e}, tau_a = {:.17e}, rel err = {:.3e})",
        a_end,
        predicted,
        a0,
        sim.t,
        tau_a,
        rel
    );
    assert!(
        a_end < a0,
        "modify_orbits_forces tau_a<0: a must shrink, a0 = {:.17e} -> a = {:.17e}",
        a0,
        a_end
    );
}

/// The same effect with +tau_a must grow a at exactly the reciprocal rate.
/// exp(+t/T) * exp(-t/T) = 1, so a_grow * a_shrink = a0^2 independently of the
/// value of T or of the integration accuracy model.
#[test]
fn mof_tau_a_sign_flip_is_reciprocal() {
    let tau = 5.0e3;
    let t_end = 100.0 * TWO_PI;

    let run = |tau_a: f64| -> (f64, f64) {
        let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2);
        rebx_attach(&mut sim);
        let force = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
        rebx_add_force(&mut sim, force);
        set_p(&mut sim, rebx_ap::particle(1), "tau_a", tau_a);
        let a0 = orb(&sim, 1).a;
        reb_simulation_integrate(&mut sim, t_end);
        (a0, orb(&sim, 1).a)
    };

    let (a0_p, a_plus) = run(tau);
    let (a0_m, a_minus) = run(-tau);
    assert!(
        (a0_p - a0_m).abs() < 1.0e-15,
        "the two runs must start from the same a: {:.17e} vs {:.17e}",
        a0_p,
        a0_m
    );
    assert!(
        a_plus > a0_p,
        "modify_orbits_forces tau_a>0 must grow a: {:.17e} -> {:.17e}",
        a0_p,
        a_plus
    );
    assert!(
        a_minus < a0_m,
        "modify_orbits_forces tau_a<0 must shrink a: {:.17e} -> {:.17e}",
        a0_m,
        a_minus
    );

    let product = a_plus * a_minus;
    let rel = (product - a0_p * a0_p).abs() / (a0_p * a0_p);
    assert!(
        rel < 1.0e-6,
        "modify_orbits_forces: a(+tau)*a(-tau) = {:.17e} should equal a0^2 = {:.17e} \
         (rel err = {:.3e}; a_plus = {:.17e}, a_minus = {:.17e})",
        product,
        a0_p * a0_p,
        rel,
        a_plus,
        a_minus
    );
}

/// tau_a = INFINITY makes invtau_a = 1/inf = 0, and with tau_e and tau_inc unset
/// the whole eccentricity/inclination branch is skipped, so the effect adds
/// exactly 0.0 to every acceleration component. The trajectory must then be the
/// bit-for-bit identical plain-gravity trajectory.
#[test]
fn mof_infinite_tau_a_is_a_bit_exact_no_op() {
    // Deliberately generic angles so that no coordinate or velocity component
    // is exactly zero (adding +0.0 to a -0.0 would flip a sign bit).
    let mk = || {
        star_and_body(
            1.0, 0.0, 1.0e-4, 0.0, 1.0, 0.17, 0.31, 0.63, 0.41, 0.77, "ias15", 1.0e-2,
        )
    };
    let t_end = 40.0 * TWO_PI;

    let mut plain = mk();
    reb_simulation_integrate(&mut plain, t_end);

    let mut damped = mk();
    rebx_attach(&mut damped);
    let force = rebx_load_force(&mut damped, "modify_orbits_forces").expect("force loads");
    rebx_add_force(&mut damped, force);
    set_p(&mut damped, rebx_ap::particle(1), "tau_a", f64::INFINITY);
    reb_simulation_integrate(&mut damped, t_end);

    let a = state_bits(&plain);
    let b = state_bits(&damped);
    assert_eq!(
        a.len(),
        b.len(),
        "state vectors must have the same length: {} vs {}",
        a.len(),
        b.len()
    );
    for k in 0..a.len() {
        assert_eq!(
            a[k], b[k],
            "modify_orbits_forces with tau_a = INFINITY must reproduce plain gravity \
             bit-for-bit; word {} differs: plain {:016x} vs damped {:016x}",
            k, a[k], b[k]
        );
    }
}

/// With only tau_e set the acceleration is prefac*r with
/// prefac = 2 (v.r)/(r^2 tau_e): a purely RADIAL force. Radial forces exert no
/// torque, so |r x v| is exactly conserved. Writing the radial excursion as a
/// harmonic oscillator of frequency n, the force is 2 rdot/tau_e along rhat, so
/// the epicyclic energy (proportional to e^2) obeys dE/dt = 2E/tau_e and hence
/// e(t) = e0 exp(t/tau_e). Negative tau_e therefore damps e.
#[test]
fn mof_negative_tau_e_damps_e_at_constant_angular_momentum() {
    let tau_e = -3.0e3;
    let e0_nominal = 0.05;
    let mut sim = star_and_body(
        1.0, 0.0, 0.0, 0.0, 1.0, e0_nominal, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2,
    );
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
    rebx_add_force(&mut sim, force);
    set_p(&mut sim, rebx_ap::particle(1), "tau_e", tau_e);

    let e0 = orb(&sim, 1).e;
    let a0 = orb(&sim, 1).a;
    let h0 = hmag(&sim, 1);

    let t_end = 3.0e3; // exactly one e-folding
    let series = sample(&mut sim, t_end, 10, |s| orb(s, 1).e);
    assert_monotone("modify_orbits_forces tau_e<0: e", &series, false);

    // Exact invariant: no torque from a radial force.
    let h1 = hmag(&sim, 1);
    let hrel = (h1 - h0).abs() / h0;
    assert!(
        hrel < 1.0e-10,
        "modify_orbits_forces tau_e: |r x v| must be exactly conserved by the purely \
         radial e-damping force; {:.17e} -> {:.17e} (rel change {:.3e})",
        h0,
        h1,
        hrel
    );

    let e_end = *series.last().unwrap();
    let predicted = e0 * (sim.t / tau_e).exp();
    let rel = (e_end - predicted).abs() / predicted;
    assert!(
        rel < 0.03,
        "modify_orbits_forces tau_e<0: e(t) = {:.17e} but e0*exp(t/tau_e) = {:.17e} \
         (e0 = {:.17e}, t = {:.17e}, tau_e = {:.17e}, rel err = {:.3e})",
        e_end,
        predicted,
        e0,
        sim.t,
        tau_e,
        rel
    );

    // h = sqrt(mu a (1-e^2)) conserved with e falling forces a to fall too.
    let a1 = orb(&sim, 1).a;
    assert!(
        a1 < a0,
        "modify_orbits_forces tau_e<0: constant-h eccentricity damping must also shrink a; \
         a {:.17e} -> {:.17e}",
        a0,
        a1
    );
    let a_from_h = a0 * (1.0 - e0 * e0) / (1.0 - e_end * e_end);
    let arel = (a1 - a_from_h).abs() / a_from_h;
    assert!(
        arel < 1.0e-6,
        "modify_orbits_forces tau_e<0: a must follow a0(1-e0^2)/(1-e^2) = {:.17e}, got {:.17e} \
         (rel err {:.3e})",
        a_from_h,
        a1,
        arel
    );
}

/// With only tau_inc set the acceleration is (0, 0, 2 vz / tau_inc). The
/// vertical motion is a harmonic oscillator of frequency n with energy
/// proportional to inc^2, and dEz/dt = vz * (2 vz/tau_inc) = 2 Ez/tau_inc gives
/// inc(t) = inc0 exp(t/tau_inc). Negative tau_inc damps the inclination.
#[test]
fn mof_negative_tau_inc_damps_inclination() {
    let tau_inc = -4.0e3;
    let mut sim = star_and_body(
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.05, 0.0, 0.0, 0.0, "ias15", 1.0e-2,
    );
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
    rebx_add_force(&mut sim, force);
    set_p(&mut sim, rebx_ap::particle(1), "tau_inc", tau_inc);

    let inc0 = orb(&sim, 1).inc;
    let t_end = 4.0e3; // one e-folding
    let series = sample(&mut sim, t_end, 10, |s| orb(s, 1).inc);
    assert_monotone("modify_orbits_forces tau_inc<0: inc", &series, false);

    let inc_end = *series.last().unwrap();
    let predicted = inc0 * (sim.t / tau_inc).exp();
    let rel = (inc_end - predicted).abs() / predicted;
    assert!(
        rel < 0.03,
        "modify_orbits_forces tau_inc<0: inc(t) = {:.17e} but inc0*exp(t/tau_inc) = {:.17e} \
         (inc0 = {:.17e}, t = {:.17e}, tau_inc = {:.17e}, rel err = {:.3e})",
        inc_end,
        predicted,
        inc0,
        sim.t,
        tau_inc,
        rel
    );
}

/// inner_disk_edge: `rebx_calculate_planet_trap` returns -10 for a body well
/// inside the trap (r < dedge*(1-hedge)), and modify_orbits_forces multiplies
/// invtau_a by that factor. So a body at a = 0.5 with the edge at 1.0 and an
/// inward tau_a must migrate OUTWARD at ten times the nominal rate.
#[test]
fn mof_planet_trap_reverses_and_amplifies_migration() {
    let tau_a = -1.0e4;
    let dedge = 1.0;
    let hedge = 0.1;
    assert!(
        0.5 < dedge * (1.0 - hedge),
        "the test body at a = 0.5 must be inside the trap edge {:.17e}",
        dedge * (1.0 - hedge)
    );
    assert!(
        (rebx_calculate_planet_trap(0.5, dedge, hedge) + 10.0).abs() < 1.0e-15,
        "planet trap factor inside the edge must be exactly -10, got {:.17e}",
        rebx_calculate_planet_trap(0.5, dedge, hedge)
    );

    let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 5.0e-3);
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
    rebx_add_force(&mut sim, force);
    set_p(&mut sim, rebx_ap::particle(1), "tau_a", tau_a);
    set_p(&mut sim, rebx_ap::force(force), "ide_position", dedge);
    set_p(&mut sim, rebx_ap::force(force), "ide_width", hedge);

    let a0 = orb(&sim, 1).a;
    let t_end = 500.0;
    let series = sample(&mut sim, t_end, 8, |s| orb(s, 1).a);
    assert_monotone("planet trap: a", &series, true);

    let a_end = *series.last().unwrap();
    // invtau_a = trap/tau_a = -10/tau_a, so a(t) = a0 exp(-10 t/tau_a).
    let predicted = a0 * (-10.0 * sim.t / tau_a).exp();
    let rel = (a_end - predicted).abs() / predicted;
    assert!(
        a_end > a0,
        "planet trap: an inward tau_a inside the trap must push the body OUT; \
         a {:.17e} -> {:.17e}",
        a0,
        a_end
    );
    assert!(
        rel < 1.0e-5,
        "planet trap: a(t) = {:.17e} but a0*exp(-10 t/tau_a) = {:.17e} (rel err {:.3e})",
        a_end,
        predicted,
        rel
    );

    // Far outside the trap the factor is exactly 1, i.e. no amplification.
    assert!(
        (rebx_calculate_planet_trap(5.0, dedge, hedge) - 1.0).abs() < 1.0e-15,
        "planet trap factor far outside the edge must be exactly 1, got {:.17e}",
        rebx_calculate_planet_trap(5.0, dedge, hedge)
    );
}

// ---------------------------------------------------------------------------
// modify_orbits_direct
// ---------------------------------------------------------------------------

/// The operator does `o.a += a0*dt*invtau_a`, i.e. a <- a*(1 + dt/tau_a), once
/// before and once after every WHFast step with dt = sim.dt/2 each. WHFast is
/// an exact Kepler map for a single test particle, so after n steps the
/// semimajor axis must be exactly a0*(1 + (dt/2)/tau_a)^(2n) — a closed form
/// that can be evaluated without touching the library.
#[test]
fn mod_direct_tau_a_follows_the_exact_geometric_progression() {
    let tau_a = -2.0e3;
    let dt = 1.0e-2;
    let nsteps = 20_000usize;

    let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, 1.0, 0.02, 0.0, 0.0, 0.0, 0.9, "whfast", dt);
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_orbits_direct").expect("operator loads");
    assert_eq!(
        rebx_add_operator(&mut sim, op),
        1,
        "rebx_add_operator(modify_orbits_direct) must report success with WHFast"
    );
    set_p(&mut sim, rebx_ap::particle(1), "tau_a", tau_a);

    let a0 = orb(&sim, 1).a;
    reb_simulation_steps(&mut sim, nsteps);
    let a_end = orb(&sim, 1).a;

    let factor = 1.0 + (dt / 2.0) / tau_a;
    let predicted = a0 * factor.powi(2 * nsteps as i32);
    let rel = (a_end - predicted).abs() / predicted;
    assert!(
        rel < 1.0e-9,
        "modify_orbits_direct tau_a: a = {:.17e} but a0*(1+dt/(2 tau_a))^(2n) = {:.17e} \
         (a0 = {:.17e}, dt = {:.17e}, tau_a = {:.17e}, n = {}, rel err = {:.3e})",
        a_end,
        predicted,
        a0,
        dt,
        tau_a,
        nsteps,
        rel
    );
    assert!(
        a_end < a0,
        "modify_orbits_direct tau_a<0 must shrink a: {:.17e} -> {:.17e}",
        a0,
        a_end
    );
}

/// With only tau_omega set the operator adds 2*pi*dt/tau_omega to the argument
/// of pericentre and rebuilds the particle from the (otherwise unchanged)
/// elements. Kepler motion does not change omega, so after a time t the total
/// advance must be exactly 2*pi*t/tau_omega while a, e and inc are untouched —
/// an exact element round-trip test.
#[test]
fn mod_direct_tau_omega_is_pure_apsidal_precession() {
    let dt = 1.0e-2;
    let nsteps = 20_000usize;
    let t_total = dt * nsteps as f64;
    // Quarter turn over the run, so no 2*pi wrapping to undo.
    let tau_omega = 4.0 * t_total;

    let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.1, 0.3, 0.4, 1.1, "whfast", dt);
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_orbits_direct").expect("operator loads");
    rebx_add_operator(&mut sim, op);
    set_p(&mut sim, rebx_ap::particle(1), "tau_omega", tau_omega);

    let o0 = orb(&sim, 1);
    reb_simulation_steps(&mut sim, nsteps);
    let o1 = orb(&sim, 1);

    let expected = TWO_PI * t_total / tau_omega;
    let got = o1.omega - o0.omega;
    assert!(
        (got - expected).abs() < 1.0e-9,
        "modify_orbits_direct tau_omega: omega advanced by {:.17e}, expected 2*pi*t/tau_omega \
         = {:.17e} (t = {:.17e}, tau_omega = {:.17e})",
        got,
        expected,
        t_total,
        tau_omega
    );
    assert!(
        got > 0.0,
        "modify_orbits_direct with tau_omega > 0 must precess forwards, got {:.17e}",
        got
    );

    for (name, x0, x1, tol) in [
        ("a", o0.a, o1.a, 1.0e-9),
        ("e", o0.e, o1.e, 1.0e-9),
        ("inc", o0.inc, o1.inc, 1.0e-9),
    ] {
        let rel = (x1 - x0).abs() / x0.abs();
        assert!(
            rel < tol,
            "modify_orbits_direct tau_omega: {} must be untouched, {:.17e} -> {:.17e} \
             (rel change {:.3e})",
            name,
            x0,
            x1,
            rel
        );
    }
}

/// modify_orbits_direct sets `o.e += e0*dt/tau_e` each half step, so with
/// nothing else set e follows the exact geometric progression
/// e0*(1 + (dt/2)/tau_e)^(2n) and tau_e < 0 damps it. Because the operator
/// rebuilds the particle from the elements, a and inc must be unchanged (the
/// `p` coupling parameter is not set, so the e-a coupling term is skipped).
#[test]
fn mod_direct_negative_tau_e_damps_e_and_leaves_a_alone() {
    let tau_e = -1.0e3;
    let dt = 1.0e-2;
    let nsteps = 10_000usize;

    let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.1, 0.3, 0.4, 1.1, "whfast", dt);
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_orbits_direct").expect("operator loads");
    rebx_add_operator(&mut sim, op);
    set_p(&mut sim, rebx_ap::particle(1), "tau_e", tau_e);

    let o0 = orb(&sim, 1);
    reb_simulation_steps(&mut sim, nsteps);
    let o1 = orb(&sim, 1);

    let factor = 1.0 + (dt / 2.0) / tau_e;
    let predicted = o0.e * factor.powi(2 * nsteps as i32);
    let rel = (o1.e - predicted).abs() / predicted;
    assert!(
        o1.e < o0.e,
        "modify_orbits_direct tau_e<0 must damp e: {:.17e} -> {:.17e}",
        o0.e,
        o1.e
    );
    assert!(
        rel < 1.0e-8,
        "modify_orbits_direct tau_e: e = {:.17e} but e0*(1+dt/(2 tau_e))^(2n) = {:.17e} \
         (e0 = {:.17e}, rel err = {:.3e})",
        o1.e,
        predicted,
        o0.e,
        rel
    );
    let arel = (o1.a - o0.a).abs() / o0.a;
    assert!(
        arel < 1.0e-9,
        "modify_orbits_direct tau_e without the `p` coupling must leave a alone: \
         {:.17e} -> {:.17e} (rel change {:.3e})",
        o0.a,
        o1.a,
        arel
    );
}

// ---------------------------------------------------------------------------
// exponential_migration
// ---------------------------------------------------------------------------

/// exponential_migration applies a = dv/(2 tau) * (afin-aini)/a * exp(-t/tau).
/// By the same energy argument as modify_orbits_forces this is
/// da/dt = a * [(afin-aini)/a] exp(-t/tau)/tau = (afin-aini) exp(-t/tau)/tau,
/// which integrates to the Hahn & Malhotra form
/// a(t) = a(0) + (afin-aini)(1 - exp(-t/tau)).
#[test]
fn exponential_migration_matches_the_closed_form_track() {
    let tau = 2.0e3;
    let aini = 1.0;
    let afin = 1.5;

    let mut sim = star_and_body(1.0, 0.0, 0.0, 0.0, aini, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2);
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "exponential_migration").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(exponential_migration) must report success"
    );
    set_p(&mut sim, rebx_ap::particle(1), "em_tau_a", tau);
    set_p(&mut sim, rebx_ap::particle(1), "em_aini", aini);
    set_p(&mut sim, rebx_ap::particle(1), "em_afin", afin);

    let a_start = orb(&sim, 1).a;
    let mut series = Vec::new();
    for k in 1..=8 {
        let t = 4.0 * tau * (k as f64) / 8.0;
        reb_simulation_integrate(&mut sim, t);
        let a = orb(&sim, 1).a;
        let predicted = a_start + (afin - aini) * (1.0 - (-t / tau).exp());
        let rel = (a - predicted).abs() / predicted;
        assert!(
            rel < 3.0e-3,
            "exponential_migration: at t = {:.17e} a = {:.17e} but the closed form \
             a0+(afin-aini)(1-exp(-t/tau)) = {:.17e} (rel err {:.3e})",
            t,
            a,
            predicted,
            rel
        );
        series.push(a);
    }
    assert_monotone("exponential_migration outward: a", &series, true);

    // After 4 e-foldings the body must have essentially arrived at afin.
    let a_end = *series.last().unwrap();
    assert!(
        (a_end - afin).abs() / afin < 0.02,
        "exponential_migration: after 4 tau the body should sit at afin = {:.17e}, got {:.17e}",
        afin,
        a_end
    );
}

/// The same effect run with afin < aini must migrate inward, and the two runs
/// must be mirror images: (a_out - a0) = -(a_in - a0) when |afin-aini| matches,
/// because da/dt is linear in (afin - aini).
#[test]
fn exponential_migration_direction_follows_afin_minus_aini() {
    let tau = 2.0e3;
    let aini = 1.0;
    let delta = 0.2;
    let t_end = tau;

    let run = |afin: f64| -> (f64, f64) {
        let mut sim =
            star_and_body(1.0, 0.0, 0.0, 0.0, aini, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2);
        rebx_attach(&mut sim);
        let force = rebx_load_force(&mut sim, "exponential_migration").expect("force loads");
        rebx_add_force(&mut sim, force);
        set_p(&mut sim, rebx_ap::particle(1), "em_tau_a", tau);
        set_p(&mut sim, rebx_ap::particle(1), "em_aini", aini);
        set_p(&mut sim, rebx_ap::particle(1), "em_afin", afin);
        let a0 = orb(&sim, 1).a;
        reb_simulation_integrate(&mut sim, t_end);
        (a0, orb(&sim, 1).a)
    };

    let (a0_out, a_out) = run(aini + delta);
    let (a0_in, a_in) = run(aini - delta);
    let d_out = a_out - a0_out;
    let d_in = a_in - a0_in;

    assert!(
        d_out > 0.0,
        "exponential_migration with afin > aini must move outward, delta a = {:.17e}",
        d_out
    );
    assert!(
        d_in < 0.0,
        "exponential_migration with afin < aini must move inward, delta a = {:.17e}",
        d_in
    );
    let asym = (d_out + d_in).abs() / d_out.abs();
    assert!(
        asym < 0.05,
        "exponential_migration is linear in (afin-aini), so the inward and outward drifts \
         must mirror: +{:.17e} vs {:.17e} (asymmetry {:.3e})",
        d_out,
        d_in,
        asym
    );
}

// ---------------------------------------------------------------------------
// type_I_migration
// ---------------------------------------------------------------------------

fn type_I_sim(sd0: f64, mp: f64, e: f64, inc: f64) -> (reb_simulation, usize) {
    let mut sim = star_and_body(1.0, 0.0, mp, 0.0, 1.0, e, inc, 0.0, 0.0, 0.0, "ias15", 1.0e-2);
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "type_I_migration").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(type_I_migration) must report success"
    );
    set_p(&mut sim, rebx_ap::force(force), "tIm_surface_density_1", sd0);
    set_p(&mut sim, rebx_ap::force(force), "tIm_scale_height_1", 0.05);
    set_p(&mut sim, rebx_ap::force(force), "tIm_surface_density_exponent", 1.0);
    set_p(&mut sim, rebx_ap::force(force), "tIm_flaring_index", 0.25);
    (sim, force)
}

/// type_I_migration is a pure sink: the migration term is -dv*invtau_mig (a
/// drag, with invtau_mig > 0 for a nearly circular orbit because t_mig reduces
/// to 2*wave/((2.7+1.1 s) h^2) > 0), and the eccentricity/inclination terms
/// carry the opposite sign to modify_orbits_forces, so a, e and inc must all
/// decrease monotonically.
#[test]
fn type_I_migration_damps_a_e_and_inc() {
    let (mut sim, _f) = type_I_sim(1.0e-3, 1.0e-5, 0.02, 0.02);
    let o0 = orb(&sim, 1);

    let t_end = 660.0;
    let n = 10;
    let mut a_series = Vec::new();
    let mut e_series = Vec::new();
    let mut i_series = Vec::new();
    for k in 1..=n {
        reb_simulation_integrate(&mut sim, t_end * (k as f64) / (n as f64));
        let o = orb(&sim, 1);
        a_series.push(o.a);
        e_series.push(o.e);
        i_series.push(o.inc);
    }

    assert_monotone("type_I_migration: a", &a_series, false);
    assert_monotone("type_I_migration: e", &e_series, false);
    assert_monotone("type_I_migration: inc", &i_series, false);

    assert!(
        *a_series.last().unwrap() < o0.a,
        "type_I_migration must migrate inward: a {:.17e} -> {:.17e}",
        o0.a,
        a_series.last().unwrap()
    );
    // Eccentricity damps on t_e ~ wave/0.78 while a migrates on t_mig, which is
    // ~1/h^2 = 400 times longer; so over the run e must fall much further than a.
    let a_frac = 1.0 - a_series.last().unwrap() / o0.a;
    let e_frac = 1.0 - e_series.last().unwrap() / o0.e;
    assert!(
        e_frac > 10.0 * a_frac,
        "type_I_migration: e damps on t_e ~ h^2 t_mig so it must decay far faster than a; \
         fractional decay e {:.3e} vs a {:.3e}",
        e_frac,
        a_frac
    );
}

/// The damping timescale wave = sqrt(ms^3) h^4 / (mp sd sqrt(G a)) is inversely
/// proportional to the disk surface density, and every timescale built from it
/// (t_mig included) inherits that. So the migration RATE must be exactly
/// proportional to tIm_surface_density_1: doubling the density doubles da/dt.
#[test]
fn type_I_migration_rate_scales_linearly_with_surface_density() {
    let t_end = 660.0;
    let run = |sd0: f64| -> f64 {
        let (mut sim, _f) = type_I_sim(sd0, 1.0e-5, 0.0, 0.0);
        let a0 = orb(&sim, 1).a;
        reb_simulation_integrate(&mut sim, t_end);
        a0 - orb(&sim, 1).a
    };

    let d1 = run(1.0e-3);
    let d2 = run(2.0e-3);
    assert!(
        d1 > 0.0 && d2 > 0.0,
        "type_I_migration must shrink a in both runs, got drops {:.17e} and {:.17e}",
        d1,
        d2
    );
    let ratio = d2 / d1;
    assert!(
        (ratio - 2.0).abs() < 0.05,
        "type_I_migration: doubling tIm_surface_density_1 must double the migration rate; \
         delta a = {:.17e} at sd0 and {:.17e} at 2*sd0 (ratio {:.6})",
        d1,
        d2,
        ratio
    );
}

/// wave is also inversely proportional to the planet mass, so a planet of twice
/// the mass migrates at twice the rate.
#[test]
fn type_I_migration_rate_scales_linearly_with_planet_mass() {
    let t_end = 660.0;
    let run = |mp: f64| -> f64 {
        let (mut sim, _f) = type_I_sim(1.0e-3, mp, 0.0, 0.0);
        let a0 = orb(&sim, 1).a;
        reb_simulation_integrate(&mut sim, t_end);
        a0 - orb(&sim, 1).a
    };

    let d1 = run(1.0e-6);
    let d2 = run(2.0e-6);
    let ratio = d2 / d1;
    assert!(
        (ratio - 2.0).abs() < 0.05,
        "type_I_migration: doubling the planet mass must double the migration rate; \
         delta a = {:.17e} at mp and {:.17e} at 2*mp (ratio {:.6})",
        d1,
        d2,
        ratio
    );
}

// ---------------------------------------------------------------------------
// tides_constant_time_lag
// ---------------------------------------------------------------------------

fn tides_sim(k2: Option<f64>, tau: Option<f64>, e: f64, r_p: f64) -> reb_simulation {
    let mut sim = star_and_body(
        1.0, 0.005, 1.0e-3, r_p, 0.02, e, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-4,
    );
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "tides_constant_time_lag").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(tides_constant_time_lag) must report success"
    );
    if let Some(k2) = k2 {
        set_p(&mut sim, rebx_ap::particle(1), "tctl_k2", k2);
    }
    if let Some(tau) = tau {
        set_p(&mut sim, rebx_ap::particle(1), "tctl_tau", tau);
        set_p(&mut sim, rebx_ap::particle(1), "OmegaMag", 0.0);
    }
    sim
}

/// With tctl_tau unset the effect reduces to the conservative piece of the
/// tidal potential. Its acceleration on the target is
/// rfac*ms*dx with rfac = -3G fac/dr^8, and the library's own
/// `rebx_tides_constant_time_lag_potential` returns
/// U = -G ms mt fac/(2 dr^6); differentiating U really does give that force, so
/// gravitational energy PLUS U must be conserved even though the gravitational
/// energy alone is not.
#[test]
fn tides_ctl_conservative_piece_conserves_total_energy() {
    let mut sim = tides_sim(Some(0.5), None, 0.1, 0.002);

    let total = |s: &reb_simulation| -> f64 {
        let rebx = rebx_extras_ref(s).expect("REBOUNDx attached");
        reb_simulation_energy(s) + rebx_tides_constant_time_lag_potential(s, rebx)
    };

    let h0 = total(&sim);
    let e0 = reb_simulation_energy(&sim);
    assert!(
        h0 < 0.0,
        "the test binary must be bound, H0 = {:.17e}",
        h0
    );

    let p = period(&sim, 1);
    let mut worst_h = 0.0f64;
    let mut worst_grav = 0.0f64;
    for k in 1..=200 {
        reb_simulation_integrate(&mut sim, p * (k as f64));
        worst_h = worst_h.max(((total(&sim) - h0) / h0).abs());
        worst_grav = worst_grav.max(((reb_simulation_energy(&sim) - e0) / e0).abs());
    }

    assert!(
        worst_h < 1.0e-9,
        "tides_constant_time_lag with tau unset is conservative: gravity + tidal potential \
         must be conserved, worst relative drift = {:.3e} (H0 = {:.17e})",
        worst_h,
        h0
    );
    assert!(
        worst_grav > 1.0e-5,
        "the tidal potential must actually matter here: the gravitational energy alone should \
         swing by far more than the {:.3e} seen (tidal term too weak to make the test \
         meaningful)",
        worst_grav
    );
    assert!(
        worst_grav > 1.0e3 * worst_h,
        "adding the tidal potential must be what conserves the total: gravity-only drift {:.3e} \
         vs total drift {:.3e}",
        worst_grav,
        worst_h
    );
}

/// With tctl_tau > 0 and the body not rotating (OmegaMag = 0) the lag term is a
/// pure drag: the acceleration contains -thetafac*ms*(thetadot x r), which for a
/// circular orbit is exactly antiparallel to the relative velocity. Orbital
/// energy can then only decrease, so a must fall monotonically.
#[test]
fn tides_ctl_time_lag_drag_shrinks_a_monotonically() {
    let tau = 1.0e-5;
    let mut sim = tides_sim(Some(0.5), Some(tau), 0.0, 0.002);
    let a0 = orb(&sim, 1).a;

    // The decay is violently self-accelerating -- the drag goes as 1/r^8, so
    // da/dt ~ -K a^-7 and the body reaches the star in finite time. Stay well
    // inside the linear-ish regime (collapse is ~150 orbits away here).
    let p = period(&sim, 1);
    let series = sample(&mut sim, 40.0 * p, 8, |s| orb(s, 1).a);
    assert_monotone("tides_constant_time_lag drag: a", &series, false);

    let a_end = *series.last().unwrap();
    assert!(
        a_end < a0,
        "tides_constant_time_lag with tau > 0 and Omega = 0 must shrink a: {:.17e} -> {:.17e}",
        a0,
        a_end
    );
    // Sanity: the decay has to be resolvable but still perturbative.
    let frac = (a0 - a_end) / a0;
    assert!(
        frac > 1.0e-3 && frac < 0.5,
        "tides_constant_time_lag drag: fractional decay {:.3e} is outside the range the test \
         was designed for (want 1e-3 .. 0.5)",
        frac
    );
}

/// For a circular orbit the radial correction 1 + 3 tau (r.v)/r^2 vanishes
/// identically (r.v = 0), so the ONLY tau-dependent term left is the drag
/// thetafac = -prefac*tau, which is exactly linear in tau. Doubling tctl_tau
/// must therefore double the orbital decay.
#[test]
fn tides_ctl_decay_rate_is_linear_in_tau() {
    let run = |tau: f64| -> f64 {
        let mut sim = tides_sim(Some(0.5), Some(tau), 0.0, 0.002);
        let a0 = orb(&sim, 1).a;
        let p = period(&sim, 1);
        // Short arc: da/dt ~ a^-7 makes the decay self-accelerating, so the
        // linear-in-tau statement only holds while delta a / a stays small.
        reb_simulation_integrate(&mut sim, 4.0 * p);
        a0 - orb(&sim, 1).a
    };

    let d1 = run(1.0e-5);
    let d2 = run(2.0e-5);
    assert!(
        d1 > 0.0 && d2 > 0.0,
        "both tides runs must decay, got {:.17e} and {:.17e}",
        d1,
        d2
    );
    let ratio = d2 / d1;
    assert!(
        (ratio - 2.0).abs() < 0.05,
        "tides_constant_time_lag: for a circular orbit the dissipative term is exactly linear \
         in tau, so doubling tau must double the decay; {:.17e} vs {:.17e} (ratio {:.6})",
        d1,
        d2,
        ratio
    );
}

/// Eccentricity damping: with Omega = 0 and tau > 0 the lag always removes
/// energy, and pseudo-synchronisation drives e down. The osculating e wobbles
/// once per orbit under the conservative tidal term, so compare the mean over a
/// full period rather than instantaneous values.
#[test]
fn tides_ctl_time_lag_damps_eccentricity() {
    let mut sim = tides_sim(Some(0.5), Some(2.0e-6), 0.1, 0.002);
    let p = period(&sim, 1);
    let a_start = orb(&sim, 1).a;

    let mut means = Vec::new();
    for _ in 0..6 {
        means.push(mean_over_period(&mut sim, 24, |s| orb(s, 1).e));
        let next = sim.t + 15.0 * p;
        reb_simulation_integrate(&mut sim, next);
    }
    // Guard: the run must stay in the perturbative regime (see the drag test --
    // the 1/r^8 force otherwise drives a finite-time plunge).
    let a_end = orb(&sim, 1).a;
    assert!(
        (a_start - a_end) / a_start < 0.2,
        "tides eccentricity test drifted out of the perturbative regime: a {:.17e} -> {:.17e}",
        a_start,
        a_end
    );
    assert_monotone("tides_constant_time_lag: <e> per orbit", &means, false);
    let drop = (means[0] - means[means.len() - 1]) / means[0];
    assert!(
        drop > 0.01,
        "tides_constant_time_lag with tau > 0 must damp e by a resolvable amount: \
         <e> {:.17e} -> {:.17e} (fractional drop {:.3e}); a smaller drop would not be \
         distinguishable from the per-orbit wobble",
        means[0],
        means[means.len() - 1],
        drop
    );
}

/// A body with no tctl_k2 set is skipped entirely by the effect (the C returns
/// before touching any acceleration), so loading and adding the force must
/// reproduce plain gravity bit-for-bit.
#[test]
fn tides_ctl_without_k2_is_a_bit_exact_no_op() {
    let t_end = 20.0;

    let mut plain = star_and_body(
        1.0, 0.005, 1.0e-3, 0.002, 0.02, 0.13, 0.29, 0.55, 0.37, 0.81, "ias15", 1.0e-4,
    );
    reb_simulation_integrate(&mut plain, t_end);

    let mut tidal = star_and_body(
        1.0, 0.005, 1.0e-3, 0.002, 0.02, 0.13, 0.29, 0.55, 0.37, 0.81, "ias15", 1.0e-4,
    );
    rebx_attach(&mut tidal);
    let force = rebx_load_force(&mut tidal, "tides_constant_time_lag").expect("force loads");
    rebx_add_force(&mut tidal, force);
    reb_simulation_integrate(&mut tidal, t_end);

    let a = state_bits(&plain);
    let b = state_bits(&tidal);
    for k in 0..a.len() {
        assert_eq!(
            a[k], b[k],
            "tides_constant_time_lag with no tctl_k2 set must reproduce plain gravity \
             bit-for-bit; word {} differs: plain {:016x} vs tidal {:016x}",
            k, a[k], b[k]
        );
    }
    // And the potential of an effect nobody opted into must be exactly zero.
    let rebx = rebx_extras_ref(&tidal).expect("REBOUNDx attached");
    let u = rebx_tides_constant_time_lag_potential(&tidal, rebx);
    assert_eq!(
        u, 0.0,
        "tides_constant_time_lag potential with no tctl_k2 anywhere must be exactly 0, got {:.17e}",
        u
    );
}

// ---------------------------------------------------------------------------
// radiation_forces
// ---------------------------------------------------------------------------

fn radiation_sim(beta: f64, c: f64, x: f64, vy: f64) -> reb_simulation {
    let mut sim = reb_simulation_create();
    reb_simulation_set_integrator(&mut sim, "ias15");
    sim.dt = 1.0e-3;

    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(&mut sim, star);

    let mut dust = reb_particle::default();
    dust.x = x;
    dust.vy = vy;
    reb_simulation_add(&mut sim, dust);
    sim.N_active = 1;

    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "radiation_forces").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(radiation_forces) must report success"
    );
    set_p(&mut sim, rebx_ap::force(force), "c", c);
    set_p(&mut sim, rebx_ap::particle(1), "beta", beta);
    sim
}

/// The radial part of the radiation force is beta*mu/r^2 outward, so beta = 1
/// cancels the star's gravity exactly and leaves only the O(v/c)
/// Poynting-Robertson terms. With c large the grain must travel in a straight
/// line at constant velocity.
#[test]
fn radiation_beta_one_cancels_gravity_for_a_test_particle() {
    let c = 1.0e9;
    let vy = 0.5;
    let mut sim = radiation_sim(1.0, c, 1.0, vy);
    let t_end = 3.0;
    reb_simulation_integrate(&mut sim, t_end);

    let p = sim.particles[1];
    let dx = (p.x - 1.0).abs();
    let dy = (p.y - vy * t_end).abs();
    let dz = p.z.abs();
    // Residual acceleration is of order (mu/r^2)(|rdot|+|v|)/c ~ 1/c; over time
    // t_end that displaces the grain by at most ~ t_end^2/(2c) ~ 5e-9.
    let tol = 1.0e-6;
    assert!(
        dx < tol && dy < tol && dz < tol,
        "radiation_forces with beta = 1 must cancel gravity: expected the grain at \
         ({:.17e}, {:.17e}, 0) but found ({:.17e}, {:.17e}, {:.17e})",
        1.0,
        vy * t_end,
        p.x,
        p.y,
        p.z
    );

    // Control: without the radiation force the same grain is strongly deflected.
    let mut plain = reb_simulation_create();
    reb_simulation_set_integrator(&mut plain, "ias15");
    plain.dt = 1.0e-3;
    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(&mut plain, star);
    let mut dust = reb_particle::default();
    dust.x = 1.0;
    dust.vy = vy;
    reb_simulation_add(&mut plain, dust);
    plain.N_active = 1;
    reb_simulation_integrate(&mut plain, t_end);
    let q = plain.particles[1];
    let deflection = ((q.x - 1.0).powi(2) + (q.y - vy * t_end).powi(2)).sqrt();
    assert!(
        deflection > 0.5,
        "control run: plain gravity should deflect the grain a long way from the straight \
         line, but the offset was only {:.17e}",
        deflection
    );
}

/// beta < 1 leaves an effective central mass mu(1-beta). A grain launched at
/// the circular speed for that reduced mass, sqrt(mu(1-beta)/r), must stay on a
/// circle of radius r — a number derived from the physics, not from the code.
#[test]
fn radiation_beta_reduces_the_effective_central_mass() {
    let beta = 0.5;
    let c = 1.0e9;
    let r0 = 1.0;
    let mu = 1.0f64;
    let v_circ = (mu * (1.0 - beta) / r0).sqrt();
    let mut sim = radiation_sim(beta, c, r0, v_circ);

    let t_orbit = TWO_PI * (r0 * r0 * r0 / (1.0 * (1.0 - beta))).sqrt();
    let mut worst = 0.0f64;
    for k in 1..=100 {
        reb_simulation_integrate(&mut sim, t_orbit * (k as f64) / 10.0);
        worst = worst.max((dist(&sim, 1) - r0).abs());
    }
    assert!(
        worst < 1.0e-6,
        "radiation_forces with beta = {:.3}: the grain launched at sqrt(mu(1-beta)/r) = {:.17e} \
         must stay on a circle of radius {:.17e}; worst deviation {:.3e}",
        beta,
        v_circ,
        r0,
        worst
    );
}

/// beta > 1 makes the net central force repulsive, so a grain launched at the
/// unmodified circular speed must fly away, monotonically.
#[test]
fn radiation_beta_above_one_unbinds_the_grain() {
    let mut sim = radiation_sim(2.0, 1.0e9, 1.0, 1.0);
    let series = sample(&mut sim, 20.0, 10, |s| dist(s, 1));
    assert_monotone("radiation_forces beta = 2: r", &series, true);
    assert!(
        *series.last().unwrap() > 5.0,
        "radiation_forces with beta = 2 gives a net outward force; the grain should have \
         escaped well past r = 5, but ended at r = {:.17e}",
        series.last().unwrap()
    );
}

/// `rebx_rad_calc_beta` and `rebx_rad_calc_particle_radius` are exact inverses
/// of one another: both are 3 L Q / (16 pi G M c rho X) with X the quantity
/// being solved for. Round-tripping must return the input to within rounding.
#[test]
fn radiation_beta_and_radius_helpers_round_trip() {
    let (G, c, m, l, rho, q) = (1.0, 1.0e4, 1.0, 0.5, 3.0e6, 1.0);
    for &radius in &[1.0e-10, 1.0e-8, 1.0e-6] {
        let beta = rebx_rad_calc_beta(G, c, m, l, radius, rho, q);
        let back = rebx_rad_calc_particle_radius(G, c, m, l, beta, rho, q);
        let rel = (back - radius).abs() / radius;
        assert!(
            rel < 1.0e-14,
            "rebx_rad_calc_particle_radius(rebx_rad_calc_beta(r)) must return r: \
             r = {:.17e} -> beta = {:.17e} -> r = {:.17e} (rel err {:.3e})",
            radius,
            beta,
            back,
            rel
        );
        assert!(
            beta > 0.0,
            "beta must be positive for positive luminosity, got {:.17e}",
            beta
        );
    }
    // beta scales as 1/radius: halving the grain doubles beta.
    let b1 = rebx_rad_calc_beta(G, c, m, l, 1.0e-8, rho, q);
    let b2 = rebx_rad_calc_beta(G, c, m, l, 0.5e-8, rho, q);
    let ratio = b2 / b1;
    assert!(
        (ratio - 2.0).abs() < 1.0e-12,
        "beta must scale as 1/radius: {:.17e} vs {:.17e} (ratio {:.17e})",
        b1,
        b2,
        ratio
    );
}

// ---------------------------------------------------------------------------
// yarkovsky_effect
// ---------------------------------------------------------------------------

/// Build a Sun + asteroid system with the simple (`ye_flag = +-1`) Yarkovsky
/// version, tuned so that the acceleration magnitude at r = 1 is `amp`.
fn yarkovsky_sim(flag: i32, amp: f64, c: f64, radius: f64, density: f64) -> (reb_simulation, f64) {
    // magnitude = 3 q L / (64 pi R rho c d^2) with q = 1 - albedo = 1.
    let lstar = amp * (64.0 * PI * radius * density * c) / 3.0;

    let mut sim = star_and_body(1.0, 0.0, 0.0, radius, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, "ias15", 1.0e-2);
    rebx_attach(&mut sim);
    let force = rebx_load_force(&mut sim, "yarkovsky_effect").expect("force loads");
    assert_eq!(
        rebx_add_force(&mut sim, force),
        1,
        "rebx_add_force(yarkovsky_effect) must report success"
    );
    set_p(&mut sim, rebx_ap::force(force), "ye_lstar", lstar);
    set_p(&mut sim, rebx_ap::force(force), "ye_c", c);
    set_p(&mut sim, rebx_ap::particle(1), "ye_body_density", density);
    set_p(&mut sim, rebx_ap::particle(1), "ye_albedo", 0.0);
    {
        let rebx = rebx_extras_mut(&mut sim).expect("REBOUNDx attached");
        rebx_set_param_int(rebx, rebx_ap::particle(1), "ye_flag", flag);
    }
    (sim, lstar)
}

/// The simple version puts a single 1 into the rotation matrix. For
/// `ye_flag = 1` the acceleration is A*(x/r) yhat, which for a prograde
/// circular orbit is always along +v; for `ye_flag = -1` it is A*(y/r) xhat,
/// always against +v. So +1 must push the body out and -1 must pull it in, and
/// for a circular orbit the two orbit-averaged rates are equal and opposite
/// (<cos^2> = <sin^2> = 1/2).
#[test]
fn yarkovsky_flag_sign_sets_the_migration_direction() {
    let amp = 1.0e-6;
    let t_end = 200.0 * TWO_PI;

    let run = |flag: i32| -> (f64, Vec<f64>) {
        let (mut sim, _l) = yarkovsky_sim(flag, amp, 1.0e4, 1.0e-8, 1.0e8);
        let a0 = orb(&sim, 1).a;
        let series = sample(&mut sim, t_end, 10, |s| orb(s, 1).a);
        (a0, series)
    };

    let (a0_out, out) = run(1);
    let (a0_in, inw) = run(-1);
    assert_monotone("yarkovsky ye_flag = +1: a", &out, true);
    assert_monotone("yarkovsky ye_flag = -1: a", &inw, false);

    let d_out = out.last().unwrap() - a0_out;
    let d_in = inw.last().unwrap() - a0_in;
    assert!(
        d_out > 0.0,
        "yarkovsky ye_flag = +1 must migrate outward, delta a = {:.17e}",
        d_out
    );
    assert!(
        d_in < 0.0,
        "yarkovsky ye_flag = -1 must migrate inward, delta a = {:.17e}",
        d_in
    );
    let asym = (d_out + d_in).abs() / d_out;
    assert!(
        asym < 0.02,
        "yarkovsky: on a circular orbit the +1 and -1 rates differ only by <cos^2> vs <sin^2>, \
         both 1/2, so the drifts must mirror: +{:.17e} vs {:.17e} (asymmetry {:.3e})",
        d_out,
        d_in,
        asym
    );
}

/// Independent check of the size of the drift. For a circular orbit the
/// tangential acceleration averages to A/2 along v, so
/// dE/dt = A v / 2 with E = -mu/(2a), giving da/dt = a^(3/2) A / sqrt(mu).
/// A here is fixed by construction of the test simulation, so the predicted
/// drift is computed without reference to the library's own arithmetic.
#[test]
fn yarkovsky_drift_matches_the_energy_balance_prediction() {
    let amp = 1.0e-6;
    let (mut sim, _l) = yarkovsky_sim(1, amp, 1.0e4, 1.0e-8, 1.0e8);
    let mu = sim.G * sim.particles[0].m;
    let a0 = orb(&sim, 1).a;

    let t_end = 200.0 * TWO_PI;
    reb_simulation_integrate(&mut sim, t_end);
    let a1 = orb(&sim, 1).a;

    let predicted = a0.powf(1.5) * amp / mu.sqrt() * t_end;
    let got = a1 - a0;
    let rel = (got - predicted).abs() / predicted;
    assert!(
        rel < 0.03,
        "yarkovsky: delta a = {:.17e} but the energy balance da/dt = a^1.5 A/sqrt(mu) predicts \
         {:.17e} over t = {:.17e} (A = {:.17e}, a0 = {:.17e}, rel err {:.3e})",
        got,
        predicted,
        t_end,
        amp,
        a0,
        rel
    );
}

/// The magnitude is proportional to the luminosity and inversely proportional
/// to the grain radius and density, all through the single factor
/// 3 q L / (64 pi R rho c d^2). Doubling the luminosity while holding
/// everything else fixed must double the drift.
#[test]
fn yarkovsky_drift_scales_linearly_with_luminosity() {
    let t_end = 100.0 * TWO_PI;
    let run = |amp: f64| -> f64 {
        let (mut sim, _l) = yarkovsky_sim(1, amp, 1.0e4, 1.0e-8, 1.0e8);
        let a0 = orb(&sim, 1).a;
        reb_simulation_integrate(&mut sim, t_end);
        orb(&sim, 1).a - a0
    };
    let d1 = run(1.0e-6);
    let d2 = run(2.0e-6);
    let ratio = d2 / d1;
    assert!(
        (ratio - 2.0).abs() < 0.01,
        "yarkovsky: the acceleration is linear in ye_lstar, so doubling it must double the \
         drift; {:.17e} vs {:.17e} (ratio {:.6})",
        d1,
        d2,
        ratio
    );
}

/// A body with no `ye_flag` (or no radius) fails the guard in
/// `rebx_yarkovsky_effect` and is skipped, so the run must be plain gravity
/// bit-for-bit.
#[test]
fn yarkovsky_without_flag_is_a_bit_exact_no_op() {
    let t_end = 20.0 * TWO_PI;
    let mk = || {
        star_and_body(
            1.0, 0.0, 0.0, 1.0e-8, 1.0, 0.19, 0.23, 0.47, 0.53, 0.71, "ias15", 1.0e-2,
        )
    };

    let mut plain = mk();
    reb_simulation_integrate(&mut plain, t_end);

    let mut yark = mk();
    rebx_attach(&mut yark);
    let force = rebx_load_force(&mut yark, "yarkovsky_effect").expect("force loads");
    rebx_add_force(&mut yark, force);
    set_p(&mut yark, rebx_ap::force(force), "ye_lstar", 1.0);
    set_p(&mut yark, rebx_ap::force(force), "ye_c", 1.0e4);
    set_p(&mut yark, rebx_ap::particle(1), "ye_body_density", 1.0e8);
    set_p(&mut yark, rebx_ap::particle(1), "ye_albedo", 0.0);
    // ye_flag deliberately left unset -> the guard rejects the body.
    reb_simulation_integrate(&mut yark, t_end);

    let a = state_bits(&plain);
    let b = state_bits(&yark);
    for k in 0..a.len() {
        assert_eq!(
            a[k], b[k],
            "yarkovsky_effect without ye_flag must reproduce plain gravity bit-for-bit; \
             word {} differs: plain {:016x} vs yarkovsky {:016x}",
            k, a[k], b[k]
        );
    }
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// Two identical set-ups combining three dissipative effects must produce
/// bit-identical state. Nothing in these code paths may depend on address
/// order, hashing or uninitialised memory.
#[test]
fn stacked_dissipative_effects_are_bit_reproducible() {
    let build_and_run = || -> Vec<u64> {
        let mut sim = star_and_body(
            1.0, 0.0, 1.0e-5, 1.0e-8, 1.0, 0.11, 0.07, 0.3, 0.6, 1.3, "ias15", 1.0e-2,
        );
        rebx_attach(&mut sim);

        let mo = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force loads");
        rebx_add_force(&mut sim, mo);
        let em = rebx_load_force(&mut sim, "exponential_migration").expect("force loads");
        rebx_add_force(&mut sim, em);
        let ye = rebx_load_force(&mut sim, "yarkovsky_effect").expect("force loads");
        rebx_add_force(&mut sim, ye);

        set_p(&mut sim, rebx_ap::particle(1), "tau_a", -5.0e3);
        set_p(&mut sim, rebx_ap::particle(1), "tau_e", -2.0e3);
        set_p(&mut sim, rebx_ap::particle(1), "em_tau_a", 4.0e3);
        set_p(&mut sim, rebx_ap::particle(1), "em_aini", 1.0);
        set_p(&mut sim, rebx_ap::particle(1), "em_afin", 1.2);
        set_p(&mut sim, rebx_ap::force(ye), "ye_lstar", 0.7);
        set_p(&mut sim, rebx_ap::force(ye), "ye_c", 1.0e4);
        set_p(&mut sim, rebx_ap::particle(1), "ye_body_density", 1.0e8);
        set_p(&mut sim, rebx_ap::particle(1), "ye_albedo", 0.0);
        {
            let rebx = rebx_extras_mut(&mut sim).expect("REBOUNDx attached");
            rebx_set_param_int(rebx, rebx_ap::particle(1), "ye_flag", 1);
        }

        reb_simulation_integrate(&mut sim, 60.0 * TWO_PI);
        state_bits(&sim)
    };

    let a = build_and_run();
    let b = build_and_run();
    assert_eq!(
        a.len(),
        b.len(),
        "the two runs must produce the same number of state words: {} vs {}",
        a.len(),
        b.len()
    );
    for k in 0..a.len() {
        assert_eq!(
            a[k], b[k],
            "stacked dissipative effects must be deterministic; state word {} differs between \
             identical runs: {:016x} vs {:016x}",
            k, a[k], b[k]
        );
    }
    // And the stacked run must actually have done something.
    assert!(
        a != vec![0u64; a.len()],
        "the determinism test must not be comparing an all-zero state"
    );
}
