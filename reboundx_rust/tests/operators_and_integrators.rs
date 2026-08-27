//! Integration tests for the operators_and_integrators group of reboundx_rs.
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

// ---------------------------------------------------------------------
// Custom linear forces used to give REBOUNDx's own integrators an ODE
// with a closed-form solution.
//
// The REBOUNDx integrators advance ONLY the velocities across a step;
// positions are left alone (integrate_force.c is an operator that kicks).
// So for a force that depends only on v, one call of `rebx_integrate_force`
// solves  dv/dt = A v  over [0, dt] with A constant, and every one of the
// four integrators is exactly its own linear stability function R(dt*A)
// applied to v.  Those stability functions follow from the Butcher
// tableaux alone, so they are an INDEPENDENT prediction, not a transcript
// of the library.
// ---------------------------------------------------------------------

/// Linear drag with unit rate: `a += -v`, i.e. dv/dt = -v.
/// Exact solution v(t) = v(0) * exp(-t).
fn drag_force(sim: &mut reb_simulation, _rebx: &mut rebx_extras, _force_idx: usize, N: usize) {
    for i in 0..N {
        sim.particles[i].ax -= sim.particles[i].vx;
        sim.particles[i].ay -= sim.particles[i].vy;
        sim.particles[i].az -= sim.particles[i].vz;
    }
}

/// Unit rotation about z: `a += z_hat x v`, i.e. dv/dt = i*v in the
/// complex plane v = vx + i*vy. Exact solution rotates v with |v| fixed.
fn rotation_force(sim: &mut reb_simulation, _rebx: &mut rebx_extras, _force_idx: usize, N: usize) {
    for i in 0..N {
        let vx = sim.particles[i].vx;
        let vy = sim.particles[i].vy;
        sim.particles[i].ax -= vy;
        sim.particles[i].ay += vx;
    }
}

/// A force that depends on POSITION only: `a += -r`. Positions do not
/// move during a REBOUNDx integrator step, so every stage sees exactly
/// the same acceleration.
fn harmonic_force(sim: &mut reb_simulation, _rebx: &mut rebx_extras, _force_idx: usize, N: usize) {
    for i in 0..N {
        sim.particles[i].ax -= sim.particles[i].x;
        sim.particles[i].ay -= sim.particles[i].y;
        sim.particles[i].az -= sim.particles[i].z;
    }
}

/// A force that adds nothing: the "effect switched off" control.
fn null_force(_sim: &mut reb_simulation, _rebx: &mut rebx_extras, _force_idx: usize, _N: usize) {}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

const NAME_TEST_FORCE: &str = "test_linear_force";

/// One particle, REBOUNDx attached, an `integrate_force` operator whose
/// `force` parameter points at a freshly created custom force.
/// `integrator` = `None` leaves the `integrator` parameter unset, which
/// integrate_force.c documents as defaulting to Euler.
fn linear_ode_sim(
    force_fn: rebx_force_fn,
    force_type: rebx_force_type,
    p0: reb_particle,
    integrator: Option<i32>,
) -> (reb_simulation, usize) {
    let mut sim = reb_simulation_create();
    sim.save_messages = 1; // keep the test output clean; messages still recorded
    reb_simulation_add(&mut sim, p0);
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "integrate_force").expect("integrate_force operator");
    rebx_with(&mut sim, |_s, rebx| {
        let f = rebx_create_force(rebx, NAME_TEST_FORCE);
        rebx.allocated_forces[f].update_accelerations = Some(force_fn);
        rebx.allocated_forces[f].force_type = force_type;
        rebx_set_param_force(rebx, rebx_ap::operator_(op), "force", f);
        if let Some(id) = integrator {
            rebx_set_param_int(rebx, rebx_ap::operator_(op), "integrator", id);
        }
    })
    .expect("extras attached");
    (sim, op)
}

fn moving_particle(vx: f64, vy: f64, vz: f64) -> reb_particle {
    let mut p = reb_particle::default();
    p.m = 1.0;
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    p
}

fn step_force(sim: &mut reb_simulation, op: usize, dt: f64) {
    rebx_with(sim, |s, rebx| rebx_integrate_force(s, rebx, op, dt)).expect("extras attached");
}

fn vbits(p: &reb_particle) -> (u64, u64, u64) {
    (p.vx.to_bits(), p.vy.to_bits(), p.vz.to_bits())
}

fn rel_err(got: f64, want: f64) -> f64 {
    ((got - want) / want).abs()
}

/// All four REBOUNDx integrators, with a printable name.
const INTEGRATORS: [(&str, i32); 4] = [
    ("euler", REBX_INTEGRATOR_EULER),
    ("rk2", REBX_INTEGRATOR_RK2),
    ("implicit_midpoint", REBX_INTEGRATOR_IMPLICIT_MIDPOINT),
    ("rk4", REBX_INTEGRATOR_RK4),
];

/// Star + eccentric planet, ready to step. `dt` is fixed; WHFast with a
/// two-body problem has no interaction term, so the trajectory is the
/// exact Kepler ellipse to roundoff.
fn eccentric_sim(a: f64, e: f64, f0: f64) -> reb_simulation {
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;
    sim.G = 1.0;
    reb_simulation_set_integrator(&mut sim, "whfast");
    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(&mut sim, star);
    reb_particle_set_name(&mut sim, 0, Some("star"));
    let planet = reb_particle_from_orbit(sim.G, sim.particles[0], 1e-4, a, e, 0., 0., 0., f0);
    reb_simulation_add(&mut sim, planet);
    reb_simulation_move_to_com(&mut sim);
    sim
}

/// The separation exactly as track_min_distance.c computes it.
fn separation(sim: &reb_simulation, i: usize, j: usize) -> f64 {
    let p = sim.particles[i];
    let s = sim.particles[j];
    let dx = p.x - s.x;
    let dy = p.y - s.y;
    let dz = p.z - s.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn get_double(sim: &reb_simulation, sel: rebx_ap, name: &str) -> Option<f64> {
    rebx_extras_ref(sim).and_then(|rebx| rebx_get_param_double(rebx, sel, name))
}

// =====================================================================
// modify_mass
// =====================================================================

/// Run n steps of `modify_mass` on particle 1 with the given tau and dt,
/// as a single POST step of the full timestep (dt_fraction = 1), using
/// the `none` integrator so nothing but the operator touches the state.
/// Returns (m_star, m_planet) at the end.
fn run_modify_mass(tau: f64, dt: f64, n: usize) -> (f64, f64) {
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;
    reb_simulation_set_integrator(&mut sim, "none");
    sim.dt = dt;
    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(&mut sim, star);
    let mut planet = reb_particle::default();
    planet.m = 1e-3;
    planet.x = 1.0;
    planet.vy = 1.0;
    reb_simulation_add(&mut sim, planet);

    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("modify_mass operator");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_mass", tau);
    }
    let ok = rebx_add_operator_step(&mut sim, op, 1.0, rebx_timing::REBX_TIMING_POST);
    assert_eq!(ok, 1, "rebx_add_operator_step(modify_mass) returned {}", ok);

    let _ = reb_simulation_steps(&mut sim, n);
    (sim.particles[0].m, sim.particles[1].m)
}

#[test]
fn modify_mass_reproduces_the_first_order_recurrence_exactly() {
    // modify_mass.c does exactly  p->m += p->m*dt/tau_mass  once per
    // applied step. Reproduce that recurrence here, in the same
    // arithmetic order, and demand bit-for-bit agreement. dt and tau are
    // powers of two so nothing in the setup is itself rounded.
    let tau = -32.0;
    let dt = 0.5;
    let n = 20;
    let m0 = 1e-3;

    let (m_star, m_planet) = run_modify_mass(tau, dt, n);

    let mut expected = m0;
    for _ in 0..n {
        expected = expected + expected * dt / tau;
    }

    assert_eq!(
        m_planet.to_bits(),
        expected.to_bits(),
        "modify_mass planet mass after {} steps: got {:.17e} (bits {:016x}), \
         hand-iterated recurrence gives {:.17e} (bits {:016x})",
        n,
        m_planet,
        m_planet.to_bits(),
        expected,
        expected.to_bits()
    );
    // The star has no tau_mass parameter, so the C's `if (tau_mass ==
    // NULL) continue` must leave it untouched, bit for bit.
    assert_eq!(
        m_star.to_bits(),
        1.0f64.to_bits(),
        "star has no tau_mass and must be untouched: got {:.17e}",
        m_star
    );
}

#[test]
fn modify_mass_ratio_matches_exp_minus_t_over_tau() {
    // Negative tau_mass = mass loss. Over total time T the exact
    // continuous solution is m(T)/m(0) = exp(T/tau) = exp(-T/|tau|).
    let tau = -32.0;
    let dt = 0.5;
    let n = 20usize;
    let T = dt * (n as f64); // 10.0
    let m0 = 1e-3;

    let (_, m_planet) = run_modify_mass(tau, dt, n);
    let ratio = m_planet / m0;
    let exact = (T / tau).exp();

    assert!(
        ratio < 1.0,
        "tau_mass < 0 must LOSE mass: ratio = {:.17e}",
        ratio
    );
    // First-order (explicit Euler) scheme: the relative deviation from
    // the exponential is ~ (T/tau)^2 / (2n) = 0.3125^2/40 = 2.4e-3.
    let dev = rel_err(ratio, exact);
    assert!(
        dev < 5e-3,
        "m(T)/m(0) = {:.17e} vs exp(T/tau) = {:.17e}; relative deviation {:.3e} \
         exceeds the 5e-3 allowed for a first-order scheme with T/tau = {:.4} over {} steps",
        ratio,
        exact,
        dev,
        T / tau,
        n
    );
    // ... and it really is the expected size, not accidentally zero:
    // a correct first-order update cannot match the exponential exactly.
    assert!(
        dev > 1e-3,
        "deviation from the exponential is {:.3e}; the C's explicit first-order \
         update should deviate by about (T/tau)^2/(2n) = {:.3e}",
        dev,
        (T / tau) * (T / tau) / (2.0 * n as f64)
    );
}

#[test]
fn modify_mass_deviation_from_exponential_is_first_order_in_dt() {
    // Same total time T = 10, half the step, twice the steps. An
    // explicit first-order scheme must halve its deviation.
    let tau = -32.0;
    let m0 = 1e-3;
    let T = 10.0f64;
    let exact = (T / tau).exp();

    let (_, m_coarse) = run_modify_mass(tau, 0.5, 20);
    let (_, m_fine) = run_modify_mass(tau, 0.25, 40);

    let dev_coarse = rel_err(m_coarse / m0, exact);
    let dev_fine = rel_err(m_fine / m0, exact);
    let ratio = dev_coarse / dev_fine;

    assert!(
        (ratio - 2.0).abs() < 0.2,
        "halving dt must halve the first-order deviation: dev(dt=0.5) = {:.6e}, \
         dev(dt=0.25) = {:.6e}, ratio = {:.4} (expected ~2)",
        dev_coarse,
        dev_fine,
        ratio
    );
}

#[test]
fn modify_mass_positive_tau_grows_mass_by_the_mirror_factor() {
    // tau > 0 is growth, tau < 0 is loss: the per-step factors are
    // (1 + dt/|tau|) and (1 - dt/|tau|). Since ln(1+u) < u for every
    // u != 0, BOTH runs end below their own exponential.
    let dt = 0.5;
    let n = 20usize;
    let m0 = 1e-3;
    let T = dt * (n as f64);

    let (_, m_grow) = run_modify_mass(32.0, dt, n);
    let (_, m_shrink) = run_modify_mass(-32.0, dt, n);

    assert!(
        m_grow > m0,
        "tau_mass > 0 must grow the mass: {:.17e} vs m0 = {:.17e}",
        m_grow,
        m0
    );
    assert!(
        m_shrink < m0,
        "tau_mass < 0 must shrink the mass: {:.17e} vs m0 = {:.17e}",
        m_shrink,
        m0
    );
    // (1+x)^n and (1-x)^n both sit below their exponentials, since
    // ln(1+x) < x and ln(1-x) < -x.
    let up = (T / 32.0).exp();
    let down = (-T / 32.0).exp();
    assert!(
        m_grow / m0 < up,
        "(1+dt/tau)^n = {:.17e} must be below exp(T/tau) = {:.17e}",
        m_grow / m0,
        up
    );
    assert!(
        m_shrink / m0 < down,
        "(1-dt/|tau|)^n = {:.17e} must be below exp(-T/|tau|) = {:.17e}",
        m_shrink / m0,
        down
    );
    // Both are within the first-order error band of their exponentials.
    for (label, got, want) in [("growth", m_grow / m0, up), ("loss", m_shrink / m0, down)] {
        let dev = rel_err(got, want);
        assert!(
            dev < 5e-3,
            "{}: ratio {:.17e} vs exp {:.17e}, relative deviation {:.3e} > 5e-3",
            label,
            got,
            want,
            dev
        );
    }
}

// =====================================================================
// track_min_distance
// =====================================================================

/// Build the eccentric two-body sim with a track_min_distance recorder
/// watching particle 1. Returns (sim, nsteps, pericentre distance); the
/// simulation's `dt` is one 2048th of the orbital period.
fn min_distance_sim(from: Option<&str>, with_orbit: bool) -> (reb_simulation, usize, f64) {
    let a = 1.0;
    let e = 0.9;
    let mut sim = eccentric_sim(a, e, 1.0);
    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    let nsteps = 2048usize;
    sim.dt = orb.P / (nsteps as f64);
    let q = orb.a * (1.0 - orb.e);

    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "track_min_distance").expect("track_min_distance");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        // The C only tracks particles whose min_distance is already set;
        // it is the running minimum, so it must start above every sample.
        rebx_set_param_double(rebx, rebx_ap::particle(1), "min_distance", 100.0);
        if let Some(name) = from {
            rebx_set_param_string(rebx, rebx_ap::particle(1), "min_distance_from", name);
        }
        if with_orbit {
            rebx_set_param_orbit(
                rebx,
                rebx_ap::particle(1),
                "min_distance_orbit",
                reb_orbit::default(),
            );
        }
    }
    let ok = rebx_add_operator(&mut sim, op);
    assert_eq!(ok, 1, "rebx_add_operator(track_min_distance) returned {}", ok);
    (sim, nsteps, q)
}

#[test]
fn track_min_distance_equals_the_minimum_over_post_step_samples() {
    // track_min_distance is a RECORDER run once after every timestep, so
    // the value it leaves behind must be exactly the smallest of the
    // post-step separations. Compute those separations independently in
    // a REBOUNDx-free twin and demand bit-for-bit equality — the sqrt is
    // the same expression in both places.
    let (mut sim, nsteps, _q) = min_distance_sim(None, false);
    let dt = sim.dt;
    let _ = reb_simulation_steps(&mut sim, nsteps);
    let recorded = get_double(&sim, rebx_ap::particle(1), "min_distance").expect("min_distance");

    let mut twin = eccentric_sim(1.0, 0.9, 1.0);
    twin.dt = dt;
    let mut best = 100.0f64;
    for _ in 0..nsteps {
        reb_simulation_step(&mut twin);
        let d = separation(&twin, 1, 0);
        if d < best {
            best = d;
        }
    }

    assert_eq!(
        recorded.to_bits(),
        best.to_bits(),
        "recorded min_distance {:.17e} (bits {:016x}) != independently sampled minimum \
         {:.17e} (bits {:016x}) over {} steps",
        recorded,
        recorded.to_bits(),
        best,
        best.to_bits(),
        nsteps
    );
}

#[test]
fn track_min_distance_brackets_the_pericentre_distance() {
    // Physics: no point of a Kepler ellipse lies inside q = a(1-e), so
    // the recorded minimum can never fall below q. It can only exceed q
    // by the sampling miss: the nearest sample is at most dt/2 from the
    // pericentre passage, and near pericentre
    //     d2r/dt2 = h^2/r^3 - GM/r^2 = GM*e/q^2,
    // so the overshoot is at most (1/2)(GM e/q^2)(dt/2)^2.
    let (mut sim, nsteps, q) = min_distance_sim(None, false);
    let dt = sim.dt;
    let GM = sim.G * (sim.particles[0].m + sim.particles[1].m);
    let e = 0.9;
    let _ = reb_simulation_steps(&mut sim, nsteps);
    let recorded = get_double(&sim, rebx_ap::particle(1), "min_distance").expect("min_distance");

    assert!(
        recorded >= q * (1.0 - 1e-9),
        "recorded minimum {:.17e} lies inside the pericentre distance q = {:.17e}",
        recorded,
        q
    );
    let rddot = GM * e / (q * q);
    let bound = 0.5 * rddot * (dt / 2.0) * (dt / 2.0);
    // Allow 3x the leading-order bound to cover the higher-order terms
    // of the expansion about pericentre.
    assert!(
        recorded - q < 3.0 * bound,
        "recorded minimum {:.17e} overshoots q = {:.17e} by {:.3e}; the sampling bound \
         (1/2)*GM*e/q^2*(dt/2)^2 with dt = {:.3e} is {:.3e}",
        recorded,
        q,
        recorded - q,
        dt,
        bound
    );
    // Sanity on the other side: it must actually be near pericentre, not
    // stuck at the 100.0 seed or at apocentre a(1+e) = 1.9.
    assert!(
        recorded < 0.5 * (q + 1.9),
        "recorded minimum {:.17e} is not near pericentre q = {:.17e} at all",
        recorded,
        q
    );
}

#[test]
fn track_min_distance_recorder_does_not_perturb_the_trajectory() {
    // REBX_OPERATOR_RECORDER promises the operator leaves the state
    // alone. Compare every particle coordinate against a run with no
    // REBOUNDx attached at all, bit for bit.
    let (mut sim, nsteps, _q) = min_distance_sim(None, false);
    let dt = sim.dt;
    let _ = reb_simulation_steps(&mut sim, nsteps);

    let mut twin = eccentric_sim(1.0, 0.9, 1.0);
    twin.dt = dt;
    let _ = reb_simulation_steps(&mut twin, nsteps);

    assert_eq!(
        sim.t.to_bits(),
        twin.t.to_bits(),
        "time diverged: {:.17e} vs {:.17e}",
        sim.t,
        twin.t
    );
    for i in 0..sim.N {
        let a = sim.particles[i];
        let b = twin.particles[i];
        for (name, x, y) in [
            ("x", a.x, b.x),
            ("y", a.y, b.y),
            ("z", a.z, b.z),
            ("vx", a.vx, b.vx),
            ("vy", a.vy, b.vy),
            ("vz", a.vz, b.vz),
        ] {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "particle {} {}: recorder run {:.17e} != plain run {:.17e}",
                i,
                name,
                x,
                y
            );
        }
    }
}

#[test]
fn track_min_distance_from_named_particle_matches_the_default_source() {
    // With min_distance_from unset the C measures from particles[0].
    // Naming particles[0] explicitly must therefore change nothing.
    let (mut a, nsteps, _q) = min_distance_sim(None, false);
    let (mut b, nsteps_b, _) = min_distance_sim(Some("star"), false);
    assert_eq!(nsteps, nsteps_b);
    let _ = reb_simulation_steps(&mut a, nsteps);
    let _ = reb_simulation_steps(&mut b, nsteps_b);

    let da = get_double(&a, rebx_ap::particle(1), "min_distance").expect("min_distance a");
    let db = get_double(&b, rebx_ap::particle(1), "min_distance").expect("min_distance b");
    assert_eq!(
        da.to_bits(),
        db.to_bits(),
        "min_distance_from=\"star\" gave {:.17e} but the default source gave {:.17e}",
        db,
        da
    );
}

#[test]
fn track_min_distance_unknown_source_leaves_the_value_untouched_and_warns() {
    // reb_simulation_get_particle_by_name returns NULL: the C warns and
    // skips the update, so the seed value must survive unchanged.
    let (mut sim, nsteps, _q) = min_distance_sim(Some("no_such_particle"), false);
    let _ = reb_simulation_steps(&mut sim, nsteps);
    let d = get_double(&sim, rebx_ap::particle(1), "min_distance").expect("min_distance");
    assert_eq!(
        d.to_bits(),
        100.0f64.to_bits(),
        "unknown min_distance_from must leave the seed 100.0 alone, got {:.17e}",
        d
    );
    let warned = sim
        .messages
        .iter()
        .any(|(_, m)| m.contains("min_distance_from cannot find particle"));
    assert!(
        warned,
        "expected a 'min_distance_from cannot find particle' warning; messages were {:?}",
        sim.messages
    );
}

#[test]
fn track_min_distance_orbit_records_the_state_at_the_minimum() {
    // When min_distance_orbit has been set, the C stores
    // reb_orbit_from_particle(sim->G, *p, *source) at the same instant it
    // stores min_distance. reb_orbit_from_particle computes
    // o.d = sqrt(dx*dx+dy*dy+dz*dz) with the identical expression, so
    // o.d must equal min_distance bit for bit.
    let (mut sim, nsteps, q) = min_distance_sim(None, true);
    let _ = reb_simulation_steps(&mut sim, nsteps);

    let d = get_double(&sim, rebx_ap::particle(1), "min_distance").expect("min_distance");
    let orb = rebx_extras_ref(&sim)
        .and_then(|rebx| rebx_get_param_orbit(rebx, rebx_ap::particle(1), "min_distance_orbit"))
        .expect("min_distance_orbit");

    assert_eq!(
        orb.d.to_bits(),
        d.to_bits(),
        "stored orbit.d = {:.17e} != stored min_distance = {:.17e}",
        orb.d,
        d
    );
    // The stored orbit is a snapshot of the same Kepler ellipse.
    assert!(
        rel_err(orb.a, 1.0) < 1e-9,
        "orbit recorded at the minimum has a = {:.17e}, expected the setup value 1.0",
        orb.a
    );
    assert!(
        rel_err(orb.e, 0.9) < 1e-9,
        "orbit recorded at the minimum has e = {:.17e}, expected the setup value 0.9",
        orb.e
    );
    // At (nearly) pericentre the true anomaly is near 0 mod 2*pi, so
    // cos f must be very close to +1. Derive it from d itself:
    // d = a(1-e^2)/(1+e cos f).
    let cosf = (orb.a * (1.0 - orb.e * orb.e) / d - 1.0) / orb.e;
    assert!(
        cosf > 0.999,
        "recorded orbit is not at pericentre: cos f = {:.6} (d = {:.17e}, q = {:.17e})",
        cosf,
        d,
        q
    );
}

#[test]
fn track_min_distance_is_bit_deterministic() {
    let mut out = Vec::new();
    for _ in 0..2 {
        let (mut sim, nsteps, _q) = min_distance_sim(None, false);
        let _ = reb_simulation_steps(&mut sim, nsteps);
        let d = get_double(&sim, rebx_ap::particle(1), "min_distance").expect("min_distance");
        out.push((d.to_bits(), sim.particles[1].x.to_bits(), sim.t.to_bits()));
    }
    assert_eq!(
        out[0], out[1],
        "two identical track_min_distance runs disagreed: {:?} vs {:?}",
        out[0], out[1]
    );
}

// =====================================================================
// integrate_force and REBOUNDx's own integrators
// =====================================================================

#[test]
fn integrate_force_defaults_to_euler_when_the_integrator_param_is_unset() {
    // integrate_force.c: `int integrator = REBX_INTEGRATOR_EULER;` then
    // overwrite only if the parameter is present.
    let dt = 0.015625;
    let (mut a, opa) = linear_ode_sim(
        drag_force,
        rebx_force_type::REBX_FORCE_VEL,
        moving_particle(1.0, 0.25, -0.5),
        None,
    );
    let (mut b, opb) = linear_ode_sim(
        drag_force,
        rebx_force_type::REBX_FORCE_VEL,
        moving_particle(1.0, 0.25, -0.5),
        Some(REBX_INTEGRATOR_EULER),
    );
    for _ in 0..8 {
        step_force(&mut a, opa, dt);
        step_force(&mut b, opb, dt);
    }
    assert_eq!(
        vbits(&a.particles[0]),
        vbits(&b.particles[0]),
        "unset `integrator` gave v = ({:.17e}, {:.17e}, {:.17e}) but explicit EULER gave \
         ({:.17e}, {:.17e}, {:.17e})",
        a.particles[0].vx,
        a.particles[0].vy,
        a.particles[0].vz,
        b.particles[0].vx,
        b.particles[0].vy,
        b.particles[0].vz
    );
}

#[test]
fn integrate_force_resets_accelerations_before_integrating() {
    // integrate_force.c calls rebx_reset_accelerations(sim->particles,
    // sim->N) before handing the force to an integrator, so whatever was
    // sitting in a is discarded and cannot leak into the kick.
    let dt = 0.015625;
    for (name, id) in INTEGRATORS {
        let (mut sim, op) =
            linear_ode_sim(drag_force, rebx_force_type::REBX_FORCE_VEL, moving_particle(1.0, 0.0, 0.0), Some(id));
        sim.particles[0].ax = 1e3;
        sim.particles[0].ay = -7.0;
        sim.particles[0].az = 42.0;
        step_force(&mut sim, op, dt);

        let (mut clean, opc) =
            linear_ode_sim(drag_force, rebx_force_type::REBX_FORCE_VEL, moving_particle(1.0, 0.0, 0.0), Some(id));
        step_force(&mut clean, opc, dt);

        assert_eq!(
            vbits(&sim.particles[0]),
            vbits(&clean.particles[0]),
            "{}: a stale acceleration leaked into the step; vx {:.17e} vs {:.17e}",
            name,
            sim.particles[0].vx,
            clean.particles[0].vx
        );
        // Only the drag acted, so vy and vz (which started at 0) stay 0.
        assert_eq!(
            sim.particles[0].vy, 0.0,
            "{}: vy became {:.17e} from a stale ay",
            name, sim.particles[0].vy
        );
        assert_eq!(
            sim.particles[0].vz, 0.0,
            "{}: vz became {:.17e} from a stale az",
            name, sim.particles[0].vz
        );
    }
}

#[test]
fn integrate_force_without_a_force_parameter_only_resets_accelerations() {
    // The port documents this deviation: the C would dereference NULL, so
    // here nothing is integrated. Everything the C does BEFORE that
    // dereference — including rebx_reset_accelerations — still happens.
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;
    reb_simulation_add(&mut sim, moving_particle(3.0, -2.0, 1.0));
    sim.particles[0].ax = 5.0;
    sim.particles[0].ay = 6.0;
    sim.particles[0].az = 7.0;
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "integrate_force").expect("integrate_force operator");
    let before = vbits(&sim.particles[0]);

    step_force(&mut sim, op, 0.5);

    assert_eq!(
        vbits(&sim.particles[0]),
        before,
        "no force parameter: velocities changed to ({:.17e}, {:.17e}, {:.17e})",
        sim.particles[0].vx,
        sim.particles[0].vy,
        sim.particles[0].vz
    );
    assert_eq!(
        (
            sim.particles[0].ax,
            sim.particles[0].ay,
            sim.particles[0].az
        ),
        (0.0, 0.0, 0.0),
        "rebx_reset_accelerations must still have run; a = ({:.17e}, {:.17e}, {:.17e})",
        sim.particles[0].ax,
        sim.particles[0].ay,
        sim.particles[0].az
    );
    let errored = sim
        .messages
        .iter()
        .any(|(_, m)| m.contains("Force parameter not set in rebx_integrate operator"));
    assert!(
        errored,
        "expected the 'Force parameter not set' error; messages were {:?}",
        sim.messages
    );
}

#[test]
fn integrate_force_with_integrator_none_leaves_velocities_alone() {
    // REBX_INTEGRATOR_NONE hits the empty switch arm.
    let (mut sim, op) = linear_ode_sim(
        drag_force,
        rebx_force_type::REBX_FORCE_VEL,
        moving_particle(1.0, 2.0, 3.0),
        Some(REBX_INTEGRATOR_NONE),
    );
    let before = vbits(&sim.particles[0]);
    sim.particles[0].ax = 9.0;
    step_force(&mut sim, op, 0.5);
    assert_eq!(
        vbits(&sim.particles[0]),
        before,
        "REBX_INTEGRATOR_NONE changed the velocity to ({:.17e}, {:.17e}, {:.17e})",
        sim.particles[0].vx,
        sim.particles[0].vy,
        sim.particles[0].vz
    );
    assert_eq!(
        sim.particles[0].ax, 0.0,
        "accelerations are still reset before the (empty) NONE branch: ax = {:.17e}",
        sim.particles[0].ax
    );
}

#[test]
fn a_force_that_adds_nothing_leaves_the_state_bit_identical() {
    // The "switch the effect off" control: with a force whose
    // update_accelerations contributes zero, every integrator must be the
    // identity on the velocities.
    for (name, id) in INTEGRATORS {
        let p = moving_particle(0.3, -1.25, 0.75);
        let (mut sim, op) =
            linear_ode_sim(null_force, rebx_force_type::REBX_FORCE_VEL, p, Some(id));
        let before = vbits(&sim.particles[0]);
        for _ in 0..4 {
            step_force(&mut sim, op, 0.125);
        }
        assert_eq!(
            vbits(&sim.particles[0]),
            before,
            "{}: a null force changed v to ({:.17e}, {:.17e}, {:.17e}) from ({:.17e}, {:.17e}, {:.17e})",
            name,
            sim.particles[0].vx,
            sim.particles[0].vy,
            sim.particles[0].vz,
            p.vx,
            p.vy,
            p.vz
        );
        assert_eq!(
            sim.particles[0].x.to_bits(),
            0.0f64.to_bits(),
            "{}: the REBOUNDx integrators must not move positions; x = {:.17e}",
            name,
            sim.particles[0].x
        );
    }
}

#[test]
fn one_step_amplification_factors_match_the_butcher_tableaux() {
    // dv/dt = -v with v(0) = 1: after one step of size x each integrator
    // applies its own linear stability function. These come from the
    // tableaux (Ralston RK2: b = (1/4, 3/4), c2 = 2/3; classic RK4;
    // implicit midpoint = the (1,1) Pade approximant), NOT from the code.
    let x = 0.015625_f64; // 2^-6
    let exact = (-x).exp();
    let cases: [(&str, i32, f64); 4] = [
        ("euler", REBX_INTEGRATOR_EULER, 1.0 - x),
        ("rk2", REBX_INTEGRATOR_RK2, 1.0 - x + x * x / 2.0),
        (
            "implicit_midpoint",
            REBX_INTEGRATOR_IMPLICIT_MIDPOINT,
            (1.0 - x / 2.0) / (1.0 + x / 2.0),
        ),
        (
            "rk4",
            REBX_INTEGRATOR_RK4,
            1.0 - x + x * x / 2.0 - x * x * x / 6.0 + x * x * x * x / 24.0,
        ),
    ];
    for (name, id, predicted) in cases {
        let (mut sim, op) = linear_ode_sim(
            drag_force,
            rebx_force_type::REBX_FORCE_VEL,
            moving_particle(1.0, 0.0, 0.0),
            Some(id),
        );
        step_force(&mut sim, op, x);
        let got = sim.particles[0].vx;
        let err = rel_err(got, predicted);
        assert!(
            err < 1e-14,
            "{}: one step of dt = {} on dv/dt = -v gave v = {:.17e}, the stability \
             function predicts {:.17e} (relative difference {:.3e}); exp(-dt) = {:.17e}",
            name,
            x,
            got,
            predicted,
            err,
            exact
        );
        // Each stability function is a distinct truncation, so the four
        // answers are genuinely different from one another.
        assert!(
            (got - exact).abs() > 0.0,
            "{}: a finite-order method cannot land exactly on exp(-dt)",
            name
        );
    }
}

/// Integrate dv/dt = -v from v = 1 over t in [0, 1] with `n` equal steps
/// and return the signed relative error against exp(-1).
fn decay_relative_error(integrator: i32, n: usize) -> f64 {
    let (mut sim, op) = linear_ode_sim(
        drag_force,
        rebx_force_type::REBX_FORCE_VEL,
        moving_particle(1.0, 0.0, 0.0),
        Some(integrator),
    );
    let dt = 1.0 / (n as f64); // n is a power of two, so dt and n*dt are exact
    for _ in 0..n {
        step_force(&mut sim, op, dt);
    }
    let exact = (-1.0f64).exp();
    (sim.particles[0].vx - exact) / exact
}

#[test]
fn integrators_show_their_expected_order_of_convergence() {
    // Global order p: halving the step must divide the error by 2^p.
    // Euler is order 1; Ralston RK2 and the implicit midpoint rule are
    // order 2; classic RK4 is order 4.
    let cases: [(&str, i32, i32); 4] = [
        ("euler", REBX_INTEGRATOR_EULER, 1),
        ("rk2", REBX_INTEGRATOR_RK2, 2),
        (
            "implicit_midpoint",
            REBX_INTEGRATOR_IMPLICIT_MIDPOINT,
            2,
        ),
        ("rk4", REBX_INTEGRATOR_RK4, 4),
    ];
    for (name, id, order) in cases {
        let e_coarse = decay_relative_error(id, 64).abs();
        let e_fine = decay_relative_error(id, 128).abs();
        assert!(
            e_fine > 0.0,
            "{}: error at n = 128 underflowed to zero, the test can say nothing",
            name
        );
        let ratio = e_coarse / e_fine;
        let expected = 2f64.powi(order);
        assert!(
            rel_err(ratio, expected) < 0.15,
            "{} (order {}): |err| went {:.6e} -> {:.6e} when the step was halved, \
             ratio {:.4}, expected about {:.1}",
            name,
            order,
            e_coarse,
            e_fine,
            ratio,
            expected
        );
    }
}

#[test]
fn higher_order_integrators_are_more_accurate_on_exponential_decay() {
    // Same problem, same step: the error must fall as the order rises.
    let n = 64;
    let euler = decay_relative_error(REBX_INTEGRATOR_EULER, n).abs();
    let rk2 = decay_relative_error(REBX_INTEGRATOR_RK2, n).abs();
    let im = decay_relative_error(REBX_INTEGRATOR_IMPLICIT_MIDPOINT, n).abs();
    let rk4 = decay_relative_error(REBX_INTEGRATOR_RK4, n).abs();

    assert!(
        euler > rk2,
        "euler error {:.6e} should exceed rk2 error {:.6e}",
        euler,
        rk2
    );
    assert!(
        euler > im,
        "euler error {:.6e} should exceed implicit_midpoint error {:.6e}",
        euler,
        im
    );
    assert!(
        rk2 > rk4,
        "rk2 error {:.6e} should exceed rk4 error {:.6e}",
        rk2,
        rk4
    );
    assert!(
        im > rk4,
        "implicit_midpoint error {:.6e} should exceed rk4 error {:.6e}",
        im,
        rk4
    );
    // The gaps are the order gaps, not noise: with n = 64 steps of
    // x = 1/64 the leading terms are n*x^2/2 ~ 7.8e-3 (euler),
    // n*x^3/6 ~ 4.1e-5 (rk2), n*x^3/12 ~ 2.0e-5 (midpoint) and
    // n*x^5/120 ~ 5.0e-10 (rk4).
    assert!(
        euler / rk2 > 50.0,
        "euler/rk2 error ratio is only {:.2} ({:.6e} vs {:.6e}); an order gap of one \
         at x = 1/64 should be ~190",
        euler / rk2,
        euler,
        rk2
    );
    assert!(
        rk2 / rk4 > 1e4,
        "rk2/rk4 error ratio is only {:.3e} ({:.6e} vs {:.6e}); an order gap of two \
         at x = 1/64 should be ~8e4",
        rk2 / rk4,
        rk2,
        rk4
    );
}

/// Run n steps of the unit-rotation force with step `theta` starting from
/// v = (1, 0, 0) and return (vx, vy).
fn rotate_n(integrator: i32, theta: f64, n: usize) -> (f64, f64) {
    let (mut sim, op) = linear_ode_sim(
        rotation_force,
        rebx_force_type::REBX_FORCE_VEL,
        moving_particle(1.0, 0.0, 0.0),
        Some(integrator),
    );
    for _ in 0..n {
        step_force(&mut sim, op, theta);
    }
    (sim.particles[0].vx, sim.particles[0].vy)
}

#[test]
fn rotation_step_matches_each_methods_complex_amplification_factor() {
    // dv/dt = i*v in the complex plane v = vx + i*vy. One step multiplies
    // v by R(i*theta), so after n steps v = R^n with v(0) = 1. R comes
    // straight from each tableau.
    let theta = 0.015625_f64; // 2^-6
    let n = 64usize;
    let t2 = theta * theta;
    let t3 = t2 * theta;
    let t4 = t2 * t2;
    let cases: [(&str, i32, f64, f64); 4] = [
        ("euler", REBX_INTEGRATOR_EULER, 1.0, theta),
        ("rk2", REBX_INTEGRATOR_RK2, 1.0 - t2 / 2.0, theta),
        (
            "rk4",
            REBX_INTEGRATOR_RK4,
            1.0 - t2 / 2.0 + t4 / 24.0,
            theta - t3 / 6.0,
        ),
        // (1 + i*theta/2)/(1 - i*theta/2), rationalised.
        (
            "implicit_midpoint",
            REBX_INTEGRATOR_IMPLICIT_MIDPOINT,
            (1.0 - t2 / 4.0) / (1.0 + t2 / 4.0),
            theta / (1.0 + t2 / 4.0),
        ),
    ];
    for (name, id, re, im) in cases {
        let rho = (re * re + im * im).sqrt();
        let phi = im.atan2(re);
        let want_r = rho.powi(n as i32);
        let want_x = want_r * ((n as f64) * phi).cos();
        let want_y = want_r * ((n as f64) * phi).sin();

        let (vx, vy) = rotate_n(id, theta, n);
        let got_r = (vx * vx + vy * vy).sqrt();

        assert!(
            rel_err(got_r, want_r) < 1e-11,
            "{}: |v| after {} steps is {:.17e}, the amplification factor predicts \
             rho^n = {:.17e} (rho = {:.17e})",
            name,
            n,
            got_r,
            want_r,
            rho
        );
        assert!(
            (vx - want_x).abs() < 1e-11 && (vy - want_y).abs() < 1e-11,
            "{}: v after {} steps is ({:.17e}, {:.17e}), R^n predicts ({:.17e}, {:.17e})",
            name,
            n,
            vx,
            vy,
            want_x,
            want_y
        );
    }
}

#[test]
fn implicit_midpoint_conserves_speed_while_the_explicit_methods_do_not() {
    // dv/dt = omega x v is norm preserving. The implicit midpoint rule
    // inherits that exactly: (v1 - v0).(v1 + v0) = dt*(omega x vavg).(2
    // vavg) = 0, so |v1| = |v0| identically. The explicit methods do not:
    // their amplification factors have |R| = 1 + O(theta^(2p)).
    let theta = 0.015625_f64;
    let n = 64usize;
    let t2 = theta * theta;

    let speed = |id: i32| {
        let (vx, vy) = rotate_n(id, theta, n);
        (vx * vx + vy * vy).sqrt()
    };

    let s_im = speed(REBX_INTEGRATOR_IMPLICIT_MIDPOINT);
    let s_eu = speed(REBX_INTEGRATOR_EULER);
    let s_rk2 = speed(REBX_INTEGRATOR_RK2);
    let s_rk4 = speed(REBX_INTEGRATOR_RK4);

    assert!(
        (s_im - 1.0).abs() < 1e-13,
        "implicit midpoint must preserve |v| = 1 to roundoff; got {:.17e} \
         (deviation {:.3e}) after {} steps",
        s_im,
        s_im - 1.0,
        n
    );
    // Euler: |R|^2 = 1 + theta^2 exactly, so |v| = (1+theta^2)^(n/2).
    let want_eu = (1.0 + t2).powf((n as f64) / 2.0);
    assert!(
        rel_err(s_eu, want_eu) < 1e-12,
        "euler |v| = {:.17e}, closed form (1+theta^2)^(n/2) = {:.17e}",
        s_eu,
        want_eu
    );
    assert!(
        s_eu > 1.0,
        "euler must GAIN energy on a rotation: |v| = {:.17e}",
        s_eu
    );
    // Ralston RK2: |R|^2 = (1-theta^2/2)^2 + theta^2 = 1 + theta^4/4.
    let want_rk2 = (1.0 + t2 * t2 / 4.0).powf((n as f64) / 2.0);
    assert!(
        rel_err(s_rk2, want_rk2) < 1e-12,
        "rk2 |v| = {:.17e}, closed form (1+theta^4/4)^(n/2) = {:.17e}",
        s_rk2,
        want_rk2
    );
    // The deviations line up with the order of each method.
    let d_eu = (s_eu - 1.0).abs();
    let d_rk2 = (s_rk2 - 1.0).abs();
    let d_rk4 = (s_rk4 - 1.0).abs();
    assert!(
        d_eu > d_rk2 && d_rk2 > d_rk4,
        "speed drift must shrink with order: euler {:.3e}, rk2 {:.3e}, rk4 {:.3e}",
        d_eu,
        d_rk2,
        d_rk4
    );
    assert!(
        d_rk4 < (s_im - 1.0).abs().max(1e-9),
        "rk4 speed drift {:.3e} should be tiny compared with rk2's {:.3e}",
        d_rk4,
        d_rk2
    );
}

#[test]
fn a_position_only_force_gives_every_integrator_the_same_kick() {
    // Positions are frozen for the duration of a REBOUNDx integrator step,
    // so a REBX_FORCE_POS force is a CONSTANT acceleration over the step
    // and every consistent method must return exactly dt*a. With a = -r
    // and r = (1, 0, 0) that is dv = (-dt, 0, 0).
    let dt = 0.015625_f64;
    let mut results = Vec::new();
    for (name, id) in INTEGRATORS {
        let mut p = reb_particle::default();
        p.m = 1.0;
        p.x = 1.0;
        let (mut sim, op) =
            linear_ode_sim(harmonic_force, rebx_force_type::REBX_FORCE_POS, p, Some(id));
        step_force(&mut sim, op, dt);
        let vx = sim.particles[0].vx;
        assert!(
            rel_err(vx, -dt) < 1e-15,
            "{}: a constant acceleration -1 over dt = {} must give vx = {:.17e}, got {:.17e}",
            name,
            dt,
            -dt,
            vx
        );
        assert_eq!(
            sim.particles[0].x.to_bits(),
            1.0f64.to_bits(),
            "{}: the position must not move during a force integration; x = {:.17e}",
            name,
            sim.particles[0].x
        );
        results.push((name, vx));
    }
    // Euler and the implicit midpoint rule reach dt*a by the very same
    // floating-point route (the midpoint iterate converges after one
    // refinement because the acceleration never changes), so they agree
    // bit for bit.
    let euler = results
        .iter()
        .find(|(n, _)| *n == "euler")
        .expect("euler result")
        .1;
    let im = results
        .iter()
        .find(|(n, _)| *n == "implicit_midpoint")
        .expect("implicit_midpoint result")
        .1;
    assert_eq!(
        euler.to_bits(),
        im.to_bits(),
        "for a position-only force euler ({:.17e}) and implicit midpoint ({:.17e}) \
         must agree bit for bit",
        euler,
        im
    );
}

#[test]
fn integrator_scratch_buffers_are_pure_scratch() {
    // rk2/rk4/implicit_midpoint cache reb_particle buffers on the force's
    // parameter list between steps. They are memcpy'd from sim->particles
    // at the top of every step, so releasing them between steps (what the
    // C's free_memory hooks do) cannot change a single bit of the answer.
    let dt = 0.015625_f64;
    let nsteps = 6;
    let cases: [(&str, i32); 3] = [
        ("rk2", REBX_INTEGRATOR_RK2),
        ("rk4", REBX_INTEGRATOR_RK4),
        ("implicit_midpoint", REBX_INTEGRATOR_IMPLICIT_MIDPOINT),
    ];
    for (name, id) in cases {
        let p = moving_particle(1.0, -0.5, 0.25);
        let (mut kept, op1) = linear_ode_sim(drag_force, rebx_force_type::REBX_FORCE_VEL, p, Some(id));
        for _ in 0..nsteps {
            step_force(&mut kept, op1, dt);
        }

        let (mut freed, op2) = linear_ode_sim(drag_force, rebx_force_type::REBX_FORCE_VEL, p, Some(id));
        for _ in 0..nsteps {
            step_force(&mut freed, op2, dt);
            rebx_with(&mut freed, |_s, rebx| {
                let f = rebx_get_force(rebx, NAME_TEST_FORCE).expect("force by name");
                match id {
                    REBX_INTEGRATOR_RK2 => rebx_rk2_free_memory(rebx, f),
                    REBX_INTEGRATOR_RK4 => rebx_rk4_free_memory(rebx, f),
                    _ => rebx_im_free_memory(rebx, f),
                }
            })
            .expect("extras attached");
        }

        assert_eq!(
            vbits(&kept.particles[0]),
            vbits(&freed.particles[0]),
            "{}: releasing the scratch buffers between steps changed the answer; \
             v = ({:.17e}, {:.17e}, {:.17e}) vs ({:.17e}, {:.17e}, {:.17e})",
            name,
            kept.particles[0].vx,
            kept.particles[0].vy,
            kept.particles[0].vz,
            freed.particles[0].vx,
            freed.particles[0].vy,
            freed.particles[0].vz
        );
    }
}

#[test]
fn the_scratch_buffers_land_on_the_force_parameter_list() {
    // Each REBOUNDx integrator stores its buffers under the names the C
    // registers, on the force's own ap list (not the operator's).
    let dt = 0.015625_f64;
    let cases: [(i32, &[&str]); 3] = [
        (REBX_INTEGRATOR_RK2, &["rk2_k2"]),
        (REBX_INTEGRATOR_RK4, &["rk4_k2", "rk4_k3"]),
        (
            REBX_INTEGRATOR_IMPLICIT_MIDPOINT,
            &["im_ps_final", "im_ps_prev", "im_ps_avg"],
        ),
    ];
    for (id, names) in cases {
        let (mut sim, op) = linear_ode_sim(
            drag_force,
            rebx_force_type::REBX_FORCE_VEL,
            moving_particle(1.0, 0.0, 0.0),
            Some(id),
        );
        step_force(&mut sim, op, dt);
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        let f = rebx_get_force(rebx, NAME_TEST_FORCE).expect("force by name");
        for name in names {
            let buf = rebx_get_param_particles(rebx, rebx_ap::force(f), name);
            let len = buf.map(|b| b.len()).unwrap_or(0);
            assert_eq!(
                len,
                sim.N,
                "integrator {}: scratch buffer '{}' should hold {} particles, holds {}",
                id,
                name,
                sim.N,
                len
            );
        }
    }
}

#[test]
fn integrate_force_is_bit_deterministic() {
    // Same setup twice must give identical bits for every integrator.
    let dt = 0.015625_f64;
    for (name, id) in INTEGRATORS {
        let mut runs = Vec::new();
        for _ in 0..2 {
            let (mut sim, op) = linear_ode_sim(
                rotation_force,
                rebx_force_type::REBX_FORCE_VEL,
                moving_particle(1.0, 0.0, 0.0),
                Some(id),
            );
            for _ in 0..37 {
                step_force(&mut sim, op, dt);
            }
            runs.push(vbits(&sim.particles[0]));
        }
        assert_eq!(
            runs[0], runs[1],
            "{}: two identical runs gave different bits, {:?} vs {:?}",
            name, runs[0], runs[1]
        );
    }
}

#[test]
fn integrate_force_runs_as_a_post_timestep_operator_in_a_real_simulation() {
    // Wire integrate_force into the simulation the way a user does, with
    // a POST step of the whole timestep, and check the drag it applies is
    // the one the operator applies by hand. With the `none` integrator
    // nothing else touches the velocities, so after n steps the velocity
    // must be exactly the Euler amplification (1 - dt)^n applied to v0
    // -- computed here by the same repeated scalar update.
    let dt = 0.015625_f64;
    let n = 16usize;
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;
    reb_simulation_set_integrator(&mut sim, "none");
    sim.dt = dt;
    reb_simulation_add(&mut sim, moving_particle(1.0, 0.0, 0.0));
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "integrate_force").expect("integrate_force operator");
    rebx_with(&mut sim, |_s, rebx| {
        let f = rebx_create_force(rebx, NAME_TEST_FORCE);
        rebx.allocated_forces[f].update_accelerations = Some(drag_force);
        rebx.allocated_forces[f].force_type = rebx_force_type::REBX_FORCE_VEL;
        rebx_set_param_force(rebx, rebx_ap::operator_(op), "force", f);
        rebx_set_param_int(
            rebx,
            rebx_ap::operator_(op),
            "integrator",
            REBX_INTEGRATOR_EULER,
        );
    })
    .expect("extras attached");
    let ok = rebx_add_operator_step(&mut sim, op, 1.0, rebx_timing::REBX_TIMING_POST);
    assert_eq!(ok, 1, "rebx_add_operator_step(integrate_force) gave {}", ok);

    let _ = reb_simulation_steps(&mut sim, n);

    // Euler on dv/dt = -v: v <- v + dt*(-v), in that order.
    let mut expected = 1.0f64;
    for _ in 0..n {
        expected = expected + dt * (-expected);
    }
    assert_eq!(
        sim.particles[0].vx.to_bits(),
        expected.to_bits(),
        "operator-driven drag gave vx = {:.17e} (bits {:016x}); the hand-iterated Euler \
         update gives {:.17e} (bits {:016x})",
        sim.particles[0].vx,
        sim.particles[0].vx.to_bits(),
        expected,
        expected.to_bits()
    );
    // And it is close to the continuous answer exp(-T).
    let T = dt * (n as f64);
    assert!(
        rel_err(sim.particles[0].vx, (-T).exp()) < 2e-3,
        "vx = {:.17e} after T = {} should be near exp(-T) = {:.17e}",
        sim.particles[0].vx,
        T,
        (-T).exp()
    );
}
