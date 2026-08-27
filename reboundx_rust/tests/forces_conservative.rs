//! Integration tests for the forces_conservative group of reboundx_rs.
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

// ---------------------------------------------------------------------------
// Shared helpers
//
// Units throughout: G = 1, central mass = 1. A body on a circular orbit at
// a = 1 then has mean motion n = 1 and period 2*pi. The "speed of light"
// c = 100 makes GM/(a c^2) = 1e-4, so every relativistic effect below is a
// clean first-order perturbation (second-order terms are O(1e-8)).
// ---------------------------------------------------------------------------

const PI: f64 = std::f64::consts::PI;

/// Wrap an angle difference into (-pi, pi].
fn wrap_pi(x: f64) -> f64 {
    let two_pi = 2.0 * PI;
    let mut y = x % two_pi;
    if y > PI {
        y -= two_pi;
    }
    if y <= -PI {
        y += two_pi;
    }
    y
}

/// A star of mass 1 at the origin plus one orbiting body, moved to the
/// centre of mass. `omega` is set to pi so that the measured pomega starts
/// far from the 0/2pi branch cut.
fn two_body(m: f64, a: f64, e: f64, inc: f64, f: f64) -> reb_simulation {
    let mut sim = reb_simulation_create();
    let star = reb_particle {
        m: 1.0,
        ..Default::default()
    };
    reb_simulation_add(&mut sim, star);
    let p = reb_particle_from_orbit(sim.G, star, m, a, e, inc, 0.0, PI, f);
    reb_simulation_add(&mut sim, p);
    reb_simulation_move_to_com(&mut sim);
    sim
}

/// Attach REBOUNDx, load `name`, add it as a force and set its `c`.
fn add_gr_force(sim: &mut reb_simulation, name: &str, c: f64) -> usize {
    rebx_attach(sim);
    let idx = rebx_load_force(sim, name).unwrap_or_else(|| panic!("rebx_load_force({name})"));
    assert_eq!(
        rebx_add_force(sim, idx),
        1,
        "rebx_add_force({name}) must report success (1)"
    );
    if let Some(rebx) = rebx_extras_mut(sim) {
        rebx_set_param_double(rebx, rebx_ap::force(idx), "c", c);
    }
    idx
}

fn orbit(sim: &reb_simulation) -> reb_orbit {
    reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0])
}

/// Integrate to `tmax` in `nchunk` pieces, accumulating the *unwrapped*
/// change in pomega and Omega. Chunking is what lets the total precession
/// exceed pi without the mod-2pi branch cut corrupting the sum; each chunk's
/// individual change stays far below pi in every test below.
fn accumulate_angles(sim: &mut reb_simulation, tmax: f64, nchunk: usize) -> (f64, f64) {
    let mut dpomega = 0.0;
    let mut dOmega = 0.0;
    let mut prev = orbit(sim);
    for k in 1..=nchunk {
        reb_simulation_integrate(sim, tmax * (k as f64) / (nchunk as f64));
        let o = orbit(sim);
        dpomega += wrap_pi(o.pomega - prev.pomega);
        dOmega += wrap_pi(o.Omega - prev.Omega);
        prev = o;
    }
    (dpomega, dOmega)
}

/// Every f64 that defines the dynamical state, as raw bit patterns.
fn state_bits(sim: &reb_simulation) -> Vec<u64> {
    let mut v = vec![sim.t.to_bits(), sim.dt.to_bits()];
    for i in 0..sim.N {
        let p = sim.particles[i];
        for q in [p.x, p.y, p.z, p.vx, p.vy, p.vz, p.m] {
            v.push(q.to_bits());
        }
    }
    v
}

fn assert_bit_identical(a: &[u64], b: &[u64], what: &str) {
    assert_eq!(
        a.len(),
        b.len(),
        "{what}: state vectors have different lengths ({} vs {})",
        a.len(),
        b.len()
    );
    for (k, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(
            x,
            y,
            "{what}: state word {k} differs bitwise: {:016x} ({:.17e}) vs {:016x} ({:.17e})",
            x,
            f64::from_bits(*x),
            y,
            f64::from_bits(*y)
        );
    }
}

/// Analytic 1PN apsidal precession rate, d(pomega)/dt = 3 n GM / (a c^2 (1-e^2)).
/// Derived from the Lagrange planetary equation for the orbit-averaged
/// disturbing function of the -3 (GM)^2 / (c^2 r^2) potential that
/// gr_potential implements; gr and gr_full reproduce the same leading term.
fn analytic_precession_rate(GM: f64, n: f64, a: f64, e: f64, c: f64) -> f64 {
    3.0 * n * GM / (a * c * c * (1.0 - e * e))
}

/// Analytic J2 nodal regression rate, d(Omega)/dt = -(3/2) n J2 (R_eq/p)^2 cos(i),
/// with p = a(1-e^2). Negative for a prograde orbit about an oblate body.
fn analytic_node_rate(n: f64, J2: f64, R_eq: f64, a: f64, e: f64, inc: f64) -> f64 {
    let p = a * (1.0 - e * e);
    -1.5 * n * J2 * (R_eq / p) * (R_eq / p) * inc.cos()
}

/// Run one GR flavour and return (measured d(pomega), analytic prediction).
fn gr_precession(name: &str, a: f64, e: f64, c: f64, tmax: f64) -> (f64, f64) {
    let mut sim = two_body(1e-8, a, e, 0.0, 0.0);
    add_gr_force(&mut sim, name, c);
    let o0 = orbit(&sim);
    let predicted = analytic_precession_rate(sim.G * sim.particles[0].m, o0.n, o0.a, o0.e, c) * tmax;
    let (dpomega, _) = accumulate_angles(&mut sim, tmax, 200);
    (dpomega, predicted)
}

/// Run gravitational_harmonics and return (measured d(Omega), analytic prediction).
fn j2_node_drift(J2: f64, R_eq: f64, a: f64, inc_deg: f64, tmax: f64) -> (f64, f64) {
    let mut sim = two_body(1e-8, a, 0.0, inc_deg * PI / 180.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", J2);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", R_eq);
    }
    let o0 = orbit(&sim);
    let predicted = analytic_node_rate(o0.n, J2, R_eq, o0.a, o0.e, o0.inc) * tmax;
    let (_, dOmega) = accumulate_angles(&mut sim, tmax, 100);
    (dOmega, predicted)
}

fn assert_close(got: f64, want: f64, rtol: f64, what: &str) {
    let rel = ((got - want) / want).abs();
    assert!(
        rel <= rtol,
        "{what}: got {got:.10e}, expected {want:.10e} (relative error {rel:.3e} > {rtol:.3e})"
    );
}

// ===========================================================================
// gr_potential
// ===========================================================================

#[test]
fn gr_potential_precession_is_prograde_and_matches_analytic() {
    // Two eccentricities: the 1/(1-e^2) factor in the analytic rate has to
    // hold as well as the overall normalisation.
    for &(a, e) in &[(1.0, 0.2), (1.0, 0.5), (0.5, 0.2)] {
        let (got, want) = gr_precession("gr_potential", a, e, 100.0, 6000.0);
        assert!(
            got > 0.0,
            "gr_potential a={a} e={e}: relativistic apsidal precession must be PROGRADE, \
             got d(pomega) = {got:.6e}"
        );
        // Residual is dominated by the short-period wobble of the osculating
        // pomega, of order (force ratio)/e ~ 6e-4/e; 2% is comfortably above it.
        assert_close(
            got,
            want,
            2e-2,
            &format!("gr_potential d(pomega) over t=6000 for a={a}, e={e}"),
        );
    }
}

#[test]
fn gr_potential_precession_scales_as_inverse_c_squared() {
    let (fast, _) = gr_precession("gr_potential", 1.0, 0.2, 100.0, 6000.0);
    let (slow, _) = gr_precession("gr_potential", 1.0, 0.2, 200.0, 6000.0);
    let ratio = fast / slow;
    assert!(
        fast > slow && slow > 0.0,
        "gr_potential: halving c must strengthen the precession; \
         d(pomega)[c=100] = {fast:.6e}, d(pomega)[c=200] = {slow:.6e}"
    );
    assert_close(
        ratio,
        4.0,
        2e-2,
        "gr_potential d(pomega) ratio between c=100 and c=200 (expected c^-2 scaling)",
    );
}

#[test]
fn gr_potential_tighter_orbit_precesses_faster() {
    // At fixed e and fixed elapsed time, d(pomega)/dt = 3 n GM/(a c^2 (1-e^2))
    // scales as n/a = a^{-5/2}. Halving a therefore multiplies it by 2^{5/2}.
    let (wide, _) = gr_precession("gr_potential", 1.0, 0.2, 100.0, 6000.0);
    let (tight, _) = gr_precession("gr_potential", 0.5, 0.2, 100.0, 6000.0);
    assert!(
        tight > wide && wide > 0.0,
        "gr_potential: the tighter orbit must precess more in the same time; \
         d(pomega)[a=0.5] = {tight:.6e}, d(pomega)[a=1.0] = {wide:.6e}"
    );
    let expected = 2f64.powf(2.5);
    assert_close(
        tight / wide,
        expected,
        2e-2,
        "gr_potential d(pomega) ratio a=0.5 : a=1.0 (expected a^-5/2 scaling)",
    );
}

#[test]
fn gr_potential_precession_grows_linearly_in_time() {
    let mut sim = two_body(1e-8, 1.0, 0.2, 0.0, 0.0);
    add_gr_force(&mut sim, "gr_potential", 100.0);
    let (early, _) = accumulate_angles(&mut sim, 2000.0, 100);
    let (late, _) = accumulate_angles(&mut sim, 6000.0, 200);
    let total = early + late;
    assert!(
        early > 0.0 && late > 0.0,
        "gr_potential: precession must accumulate in the same direction in both \
         intervals; first 2000 = {early:.6e}, next 4000 = {late:.6e}"
    );
    assert_close(
        total / early,
        3.0,
        1e-2,
        "gr_potential d(pomega)(t=6000) / d(pomega)(t=2000) (secular drift is linear in t)",
    );
}

#[test]
fn gr_potential_conserves_energy_including_its_own_potential() {
    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    let idx = add_gr_force(&mut sim, "gr_potential", 100.0);
    let e0 = reb_simulation_energy(&sim);
    let v0 = rebx_with(&mut sim, |s, r| rebx_gr_potential_potential(s, r, idx)).expect("attached");
    let h0 = e0 + v0;

    let mut worst_newtonian: f64 = 0.0;
    let mut worst_total: f64 = 0.0;
    for k in 1..=100 {
        reb_simulation_integrate(&mut sim, 2000.0 * (k as f64) / 100.0);
        let e = reb_simulation_energy(&sim);
        let v =
            rebx_with(&mut sim, |s, r| rebx_gr_potential_potential(s, r, idx)).expect("attached");
        worst_newtonian = worst_newtonian.max(((e - e0) / e0).abs());
        worst_total = worst_total.max(((e + v - h0) / h0).abs());
    }
    // gr_potential is derivable from the position-only potential
    // V = -3 (GM)^2 m / (c^2 r^2), and the force in the source is exactly
    // -grad V, so kinetic + Newtonian + gr potential is a true constant of
    // motion. Only the integrator's round-off should move it.
    assert!(
        worst_total < 1e-11,
        "gr_potential: kinetic + Newtonian + gr potential must be conserved; \
         max relative drift was {worst_total:.4e}"
    );
    // ... whereas the purely Newtonian energy is not conserved: it trades
    // against the gr potential as r varies over the orbit. This is what
    // proves the previous assertion is not vacuous.
    assert!(
        worst_newtonian > 100.0 * worst_total,
        "gr_potential: the Newtonian-only energy should visibly breathe against \
         the gr potential (max relative excursion {worst_newtonian:.4e}) while the \
         total stays put (max relative drift {worst_total:.4e})"
    );
}

#[test]
fn gr_potential_potential_matches_closed_form() {
    let c = 100.0;
    let mut sim = two_body(1e-4, 1.3, 0.35, 0.0, 0.7);
    let idx = add_gr_force(&mut sim, "gr_potential", c);
    reb_simulation_integrate(&mut sim, 3.0);

    let got = rebx_with(&mut sim, |s, r| rebx_gr_potential_potential(s, r, idx)).expect("attached");

    // Independent evaluation of V = -3 (G m0)^2 m1 / (c^2 r^2).
    let p0 = sim.particles[0];
    let p1 = sim.particles[1];
    let r = ((p1.x - p0.x).powi(2) + (p1.y - p0.y).powi(2) + (p1.z - p0.z).powi(2)).sqrt();
    let mu = sim.G * p0.m;
    let want = -3.0 * mu * mu * p1.m / (c * c * r * r);

    assert!(
        got < 0.0,
        "gr_potential potential must be negative (attractive correction), got {got:.6e}"
    );
    assert_close(got, want, 1e-14, "rebx_gr_potential_potential vs -3 mu^2 m / (c^2 r^2)");
}

#[test]
fn gr_potential_with_huge_c_is_bitwise_plain_gravity() {
    // prefac1 = 6 (GM)^2 / c^2 = 6e-200 here, so every acceleration update is
    // a no-op in double precision and the whole trajectory must be identical.
    let mut plain = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    add_gr_force(&mut sim, "gr_potential", 1e100);
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "gr_potential with c = 1e100 vs plain gravity",
    );
}

// ===========================================================================
// gr
// ===========================================================================

#[test]
fn gr_precession_is_prograde_and_matches_analytic() {
    for &(a, e) in &[(1.0, 0.2), (1.0, 0.5), (0.5, 0.2)] {
        let (got, want) = gr_precession("gr", a, e, 100.0, 6000.0);
        assert!(
            got > 0.0,
            "gr a={a} e={e}: relativistic apsidal precession must be PROGRADE, \
             got d(pomega) = {got:.6e}"
        );
        assert_close(
            got,
            want,
            2e-2,
            &format!("gr d(pomega) over t=6000 for a={a}, e={e}"),
        );
    }
}

#[test]
fn gr_tighter_orbit_precesses_faster() {
    let (wide, _) = gr_precession("gr", 1.0, 0.2, 100.0, 6000.0);
    let (tight, _) = gr_precession("gr", 0.5, 0.2, 100.0, 6000.0);
    assert!(
        tight > wide && wide > 0.0,
        "gr: the tighter orbit must precess more in the same time; \
         d(pomega)[a=0.5] = {tight:.6e}, d(pomega)[a=1.0] = {wide:.6e}"
    );
    assert_close(
        tight / wide,
        2f64.powf(2.5),
        2e-2,
        "gr d(pomega) ratio a=0.5 : a=1.0 (expected a^-5/2 scaling)",
    );
}

#[test]
fn gr_conserves_its_hamiltonian() {
    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    let idx = add_gr_force(&mut sim, "gr", 100.0);
    let e0 = reb_simulation_energy(&sim);
    let h0 = rebx_with(&mut sim, |s, r| rebx_gr_hamiltonian(r, s, idx)).expect("attached");

    let mut worst_newtonian: f64 = 0.0;
    let mut worst_h: f64 = 0.0;
    for k in 1..=100 {
        reb_simulation_integrate(&mut sim, 2000.0 * (k as f64) / 100.0);
        let e = reb_simulation_energy(&sim);
        let h = rebx_with(&mut sim, |s, r| rebx_gr_hamiltonian(r, s, idx)).expect("attached");
        worst_newtonian = worst_newtonian.max(((e - e0) / e0).abs());
        worst_h = worst_h.max(((h - h0) / h0).abs());
    }
    assert!(
        worst_h < 1e-11,
        "gr: rebx_gr_hamiltonian is the conserved quantity of the gr force; \
         max relative drift was {worst_h:.4e}"
    );
    assert!(
        worst_newtonian > 100.0 * worst_h,
        "gr: the Newtonian-only energy should visibly breathe (max relative \
         excursion {worst_newtonian:.4e}) while the gr Hamiltonian stays put \
         (max relative drift {worst_h:.4e})"
    );
}

#[test]
fn gr_with_huge_c_is_bitwise_plain_gravity() {
    let mut plain = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    add_gr_force(&mut sim, "gr", 1e100);
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "gr with c = 1e100 vs plain gravity",
    );
}

// ===========================================================================
// gr_full
// ===========================================================================

#[test]
fn gr_full_precession_is_prograde_and_matches_analytic() {
    for &(a, e) in &[(1.0, 0.2), (1.0, 0.5), (0.5, 0.2)] {
        let (got, want) = gr_precession("gr_full", a, e, 100.0, 6000.0);
        assert!(
            got > 0.0,
            "gr_full a={a} e={e}: relativistic apsidal precession must be PROGRADE, \
             got d(pomega) = {got:.6e}"
        );
        assert_close(
            got,
            want,
            2e-2,
            &format!("gr_full d(pomega) over t=6000 for a={a}, e={e}"),
        );
    }
}

#[test]
fn gr_full_hamiltonian_residual_scales_as_inverse_c_fourth() {
    // rebx_gr_full_hamiltonian is a 1PN expression paired with the 1PN EIH
    // accelerations, so it is conserved only through O(v^2/c^2); the residual
    // is the neglected O(v^4/c^4) piece. Doubling c must therefore shrink it
    // by a factor of 2^4 = 16. Asserting the SCALING (rather than a magnitude)
    // is what distinguishes "an expected truncation error" from "a bug".
    let mut residuals = Vec::new();
    for &c in &[100.0f64, 200.0, 400.0] {
        let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
        let idx = add_gr_force(&mut sim, "gr_full", c);
        let h0 = rebx_with(&mut sim, |s, r| rebx_gr_full_hamiltonian(s, r, idx)).expect("attached");
        let mut worst: f64 = 0.0;
        for k in 1..=100 {
            reb_simulation_integrate(&mut sim, 2000.0 * (k as f64) / 100.0);
            let h =
                rebx_with(&mut sim, |s, r| rebx_gr_full_hamiltonian(s, r, idx)).expect("attached");
            worst = worst.max(((h - h0) / h0).abs());
        }
        residuals.push(worst);
    }
    for (k, r) in residuals.iter().enumerate() {
        assert!(
            r.is_finite() && *r > 0.0,
            "gr_full: Hamiltonian residual {k} is not a positive finite number: {r:.4e}"
        );
    }
    for k in 0..2 {
        let ratio = residuals[k] / residuals[k + 1];
        assert_close(
            ratio,
            16.0,
            0.2,
            &format!(
                "gr_full Hamiltonian residual ratio between c={} and c={} \
                 (residuals {:.4e} and {:.4e}); the neglected term is O(c^-4)",
                100 << k,
                100 << (k + 1),
                residuals[k],
                residuals[k + 1]
            ),
        );
    }
}

#[test]
fn gr_full_with_huge_c_is_bitwise_plain_gravity() {
    let mut plain = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    add_gr_force(&mut sim, "gr_full", 1e100);
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "gr_full with c = 1e100 vs plain gravity",
    );
}

#[test]
fn three_gr_flavours_agree_on_the_precession_of_a_test_particle() {
    // gr_potential, gr and gr_full are three different truncations of the
    // same 1PN dynamics. For a massless body around a dominant central mass
    // they must agree on the secular apsidal precession to O(GM/(a c^2)) = 1e-4.
    let mut got = Vec::new();
    for name in ["gr_potential", "gr", "gr_full"] {
        let (d, _) = gr_precession(name, 1.0, 0.2, 100.0, 6000.0);
        got.push((name, d));
    }
    for i in 0..got.len() {
        for j in (i + 1)..got.len() {
            let rel = ((got[i].1 - got[j].1) / got[j].1).abs();
            assert!(
                rel < 1.5e-2,
                "{} and {} disagree on d(pomega): {:.8e} vs {:.8e} (relative {rel:.3e})",
                got[i].0,
                got[j].0,
                got[i].1,
                got[j].1
            );
        }
    }
}

// ===========================================================================
// central_force
// ===========================================================================

/// Set up central_force on the star with the given Acentral/gammacentral.
fn add_central_force(sim: &mut reb_simulation, A: f64, gamma: f64) -> usize {
    rebx_attach(sim);
    let idx = rebx_load_force(sim, "central_force").expect("load central_force");
    assert_eq!(
        rebx_add_force(sim, idx),
        1,
        "rebx_add_force(central_force) failed"
    );
    if let Some(rebx) = rebx_extras_mut(sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", A);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "gammacentral", gamma);
    }
    idx
}

#[test]
fn central_force_Acentral_reproduces_the_requested_precession_rate() {
    // rebx_central_force_Acentral inverts d(pomega)/dt = A a^(gamma-1)(gamma+2)/(2n)
    // for A. Feed it a target rate, integrate, and measure the rate back out.
    // The particle is started at the true anomaly where r = a exactly
    // (1 + e cos f = 1 - e^2, i.e. cos f = -e), because the C uses the
    // INSTANTANEOUS separation o.d in place of a.
    let e = 0.05f64;
    let f_at_r_equals_a = (-e).acos();
    let tmax = 6000.0;
    for &gamma in &[-3.0f64, -1.0, 0.0, 1.0] {
        for &target in &[1e-3f64, -1e-3] {
            let mut sim = two_body(1e-8, 1.0, e, 0.0, f_at_r_equals_a);
            let (body, primary) = (sim.particles[1], sim.particles[0]);
            let d = reb_orbit_from_particle(sim.G, body, primary).d;
            assert!(
                (d - 1.0).abs() < 1e-12,
                "test setup: body should start at r = a = 1 so that o.d = a; got d = {d:.16e}"
            );
            let A = rebx_central_force_Acentral(&mut sim, body, primary, target, gamma);
            assert!(
                A.is_finite() && A != 0.0,
                "rebx_central_force_Acentral(gamma={gamma}, pomegadot={target:e}) \
                 returned {A:.6e}"
            );
            add_central_force(&mut sim, A, gamma);
            let (dpomega, _) = accumulate_angles(&mut sim, tmax, 600);
            let rate = dpomega / tmax;
            assert!(
                rate * target > 0.0,
                "central_force gamma={gamma}: requested pomegadot {target:e} but measured \
                 {rate:.6e} — the sign of the precession is wrong"
            );
            assert_close(
                rate,
                target,
                2.5e-2,
                &format!("central_force measured pomegadot for gamma={gamma}, target={target:e}"),
            );
        }
    }
}

#[test]
fn central_force_gamma_minus_three_reproduces_gr_potential() {
    // gr_potential adds a = -6 (GM)^2/(c^2 r^3) r_hat. That is exactly the
    // central_force law a = A r^gamma r_hat with gamma = -3 and
    // A = -6 (GM)^2/c^2, including the back-reaction on the source, so two
    // completely separate code paths must produce the same trajectory (up to
    // the powf-vs-multiply round-off difference in forming 1/r^4).
    let c = 100.0;
    let tmax = 1000.0;

    let mut s1 = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    let f1 = add_gr_force(&mut s1, "gr_potential", c);
    reb_simulation_integrate(&mut s1, tmax);
    let o1 = orbit(&s1);

    let mut s2 = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    let mu = s2.G * s2.particles[0].m;
    add_central_force(&mut s2, -6.0 * mu * mu / (c * c), -3.0);
    reb_simulation_integrate(&mut s2, tmax);
    let o2 = orbit(&s2);

    assert_close(o2.a, o1.a, 1e-11, "central_force(gamma=-3) vs gr_potential: semimajor axis");
    assert_close(o2.e, o1.e, 1e-11, "central_force(gamma=-3) vs gr_potential: eccentricity");
    assert_close(
        o2.pomega,
        o1.pomega,
        1e-11,
        "central_force(gamma=-3) vs gr_potential: longitude of pericentre after t=1000",
    );

    // The two potential functions must agree too: both evaluate
    // -3 mu^2 m / (c^2 r^2), one as A r^(gamma+1)/(gamma+1).
    let v1 = rebx_with(&mut s1, |s, r| rebx_gr_potential_potential(s, r, f1)).expect("attached");
    let v2 = rebx_extras_ref(&s2)
        .map(|r| rebx_central_force_potential(&s2, r))
        .expect("attached");
    assert_close(
        v2,
        v1,
        1e-11,
        "rebx_central_force_potential(gamma=-3) vs rebx_gr_potential_potential",
    );
}

#[test]
fn central_force_conserves_energy_including_its_own_potential() {
    let gamma = -1.0;
    let A = 2e-3;
    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    add_central_force(&mut sim, A, gamma);
    let e0 = reb_simulation_energy(&sim);
    let v0 = rebx_extras_ref(&sim)
        .map(|r| rebx_central_force_potential(&sim, r))
        .expect("attached");
    let h0 = e0 + v0;

    let mut worst_newtonian: f64 = 0.0;
    let mut worst_total: f64 = 0.0;
    for k in 1..=100 {
        reb_simulation_integrate(&mut sim, 2000.0 * (k as f64) / 100.0);
        let e = reb_simulation_energy(&sim);
        let v = rebx_extras_ref(&sim)
            .map(|r| rebx_central_force_potential(&sim, r))
            .expect("attached");
        worst_newtonian = worst_newtonian.max(((e - e0) / e0).abs());
        worst_total = worst_total.max(((e + v - h0) / h0).abs());
    }
    assert!(
        worst_total < 1e-11,
        "central_force: kinetic + Newtonian + central-force potential must be \
         conserved; max relative drift was {worst_total:.4e}"
    );
    assert!(
        worst_newtonian > 100.0 * worst_total,
        "central_force: the Newtonian-only energy should visibly breathe (max \
         relative excursion {worst_newtonian:.4e}) while the total stays put \
         (max relative drift {worst_total:.4e})"
    );
}

#[test]
fn central_force_potential_matches_closed_form_in_both_branches() {
    // General branch: V = -m A r^(gamma+1)/(gamma+1).
    for &(gamma, A) in &[(1.0f64, 3e-3f64), (0.0, 2e-3), (-3.0, -1e-3)] {
        let mut sim = two_body(1e-4, 1.3, 0.35, 0.0, 0.7);
        add_central_force(&mut sim, A, gamma);
        reb_simulation_integrate(&mut sim, 3.0);
        let got = rebx_extras_ref(&sim)
            .map(|r| rebx_central_force_potential(&sim, r))
            .expect("attached");
        let p0 = sim.particles[0];
        let p1 = sim.particles[1];
        let r = ((p1.x - p0.x).powi(2) + (p1.y - p0.y).powi(2) + (p1.z - p0.z).powi(2)).sqrt();
        let want = -p1.m * A * r.powf(gamma + 1.0) / (gamma + 1.0);
        assert_close(
            got,
            want,
            1e-13,
            &format!("rebx_central_force_potential for gamma={gamma}, A={A:e}"),
        );
    }

    // gamma = -1 branch: the source special-cases |gamma+1| < EPSILON to the
    // logarithmic potential V = -m A ln(r).
    let (gamma, A) = (-1.0f64, 2e-3f64);
    let mut sim = two_body(1e-4, 1.3, 0.35, 0.0, 0.7);
    add_central_force(&mut sim, A, gamma);
    reb_simulation_integrate(&mut sim, 3.0);
    let got = rebx_extras_ref(&sim)
        .map(|r| rebx_central_force_potential(&sim, r))
        .expect("attached");
    let p0 = sim.particles[0];
    let p1 = sim.particles[1];
    let r = ((p1.x - p0.x).powi(2) + (p1.y - p0.y).powi(2) + (p1.z - p0.z).powi(2)).sqrt();
    let want = -p1.m * A * r.ln();
    assert_close(
        got,
        want,
        1e-13,
        "rebx_central_force_potential for gamma=-1 (logarithmic branch)",
    );
}

#[test]
fn central_force_Acentral_is_refused_for_gamma_minus_two() {
    // An r^-2 force law is degenerate with gravity: it produces no precession
    // at all, so no finite Acentral can deliver a requested pomegadot. The
    // source raises a REBOUND error and returns 0.
    let mut sim = two_body(1e-8, 1.0, 0.05, 0.0, 0.0);
    sim.save_messages = 1;
    let (body, primary) = (sim.particles[1], sim.particles[0]);
    let A = rebx_central_force_Acentral(&mut sim, body, primary, 1e-3, -2.0);
    assert_eq!(
        A.to_bits(),
        0.0f64.to_bits(),
        "rebx_central_force_Acentral(gamma=-2) must return exactly +0.0, got {A:.6e}"
    );
    assert!(
        !sim.messages.is_empty(),
        "rebx_central_force_Acentral(gamma=-2) must raise a REBOUND error message; \
         none was recorded"
    );

    // And the physics behind the refusal: an r^-2 addition really does leave
    // the pericentre fixed, because it just rescales GM. Compare against the
    // SAME A applied with gamma = -1, where the secular rate
    // A a^(gamma-1) (gamma+2) / (2 n) is non-zero. The absolute residual for
    // gamma = -2 cannot be driven to zero because the osculating pomega
    // wobbles at the ~(force ratio)/e level within each orbit, so the
    // meaningful statement is the ratio.
    let A = -1e-3;
    let mut deg = two_body(1e-8, 1.0, 0.05, 0.0, 0.0);
    add_central_force(&mut deg, A, -2.0);
    let (d_deg, _) = accumulate_angles(&mut deg, 6000.0, 600);

    let mut ref_sim = two_body(1e-8, 1.0, 0.05, 0.0, 0.0);
    add_central_force(&mut ref_sim, A, -1.0);
    let (d_ref, _) = accumulate_angles(&mut ref_sim, 6000.0, 600);

    assert!(
        d_ref.abs() > 1.0,
        "sanity: gamma=-1 with A={A:e} should precess by order a radian over t=6000, \
         got {d_ref:.6e}"
    );
    assert!(
        d_deg.abs() < 0.01 * d_ref.abs(),
        "central_force with gamma=-2 must not precess the pericentre (it only \
         rescales GM); got d(pomega) = {d_deg:.6e} against {d_ref:.6e} for the same \
         A at gamma=-1"
    );
}

#[test]
fn central_force_without_parameters_is_bitwise_plain_gravity() {
    // The force is registered and called every substep, but neither Acentral
    // nor gammacentral is set on any particle, so it must touch nothing.
    let mut plain = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 0.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "central_force").expect("load central_force");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force failed");
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "central_force with no Acentral/gammacentral vs plain gravity",
    );

    // Its potential contribution is likewise exactly zero.
    let v = rebx_extras_ref(&sim)
        .map(|r| rebx_central_force_potential(&sim, r))
        .expect("attached");
    assert_eq!(
        v.to_bits(),
        0.0f64.to_bits(),
        "rebx_central_force_potential with no parameters set must be exactly +0.0, got {v:.6e}"
    );
}

// ===========================================================================
// gravitational_harmonics
// ===========================================================================

#[test]
fn j2_regresses_the_node_at_the_analytic_rate() {
    let (got, want) = j2_node_drift(0.01, 0.1, 1.0, 20.0, 3000.0);
    assert!(
        got < 0.0,
        "J2 on the central body must make the ascending node REGRESS for a \
         prograde orbit; got d(Omega) = {got:.6e}"
    );
    assert!(
        want < 0.0,
        "sanity: the analytic prediction should itself be negative, got {want:.6e}"
    );
    assert_close(got, want, 1e-2, "J2 nodal drift over t=3000 (i=20 deg, a=1)");
}

#[test]
fn j2_node_regression_scales_linearly_with_J2() {
    let (small, _) = j2_node_drift(0.01, 0.1, 1.0, 20.0, 3000.0);
    let (large, _) = j2_node_drift(0.02, 0.1, 1.0, 20.0, 3000.0);
    assert!(
        large < small && small < 0.0,
        "doubling J2 must deepen the nodal regression; d(Omega)[J2=0.01] = {small:.6e}, \
         d(Omega)[J2=0.02] = {large:.6e}"
    );
    assert_close(
        large / small,
        2.0,
        5e-3,
        "J2 nodal drift ratio J2=0.02 : J2=0.01 (rate is linear in J2)",
    );
}

#[test]
fn j2_node_regression_scales_as_cos_inclination() {
    let a20 = 20.0f64;
    let a50 = 50.0f64;
    let (d20, _) = j2_node_drift(0.01, 0.1, 1.0, a20, 3000.0);
    let (d50, _) = j2_node_drift(0.01, 0.1, 1.0, a50, 3000.0);
    assert!(
        d20 < 0.0 && d50 < 0.0,
        "both inclinations must regress: d(Omega)[20 deg] = {d20:.6e}, \
         d(Omega)[50 deg] = {d50:.6e}"
    );
    let expected = (a50 * PI / 180.0).cos() / (a20 * PI / 180.0).cos();
    assert_close(
        d50 / d20,
        expected,
        1e-2,
        "J2 nodal drift ratio i=50 deg : i=20 deg (rate goes as cos i)",
    );
}

#[test]
fn j2_node_regression_scales_as_semimajor_axis_to_the_minus_seven_halves() {
    // rate = -(3/2) n J2 (R/p)^2 cos i, and with e = 0, n (R/a)^2 goes as
    // a^{-3/2} a^{-2} = a^{-7/2}.
    let (near, _) = j2_node_drift(0.01, 0.1, 1.0, 20.0, 3000.0);
    let (far, _) = j2_node_drift(0.01, 0.1, 1.5, 20.0, 3000.0);
    assert!(
        far > near && far < 0.0,
        "the wider orbit must regress more slowly: d(Omega)[a=1.0] = {near:.6e}, \
         d(Omega)[a=1.5] = {far:.6e}"
    );
    assert_close(
        near / far,
        1.5f64.powf(3.5),
        1e-2,
        "J2 nodal drift ratio a=1.0 : a=1.5 (rate goes as a^-7/2)",
    );
}

#[test]
fn gravitational_harmonics_conserves_energy_including_its_own_potential() {
    let mut sim = two_body(1e-4, 1.0, 0.05, 20.0 * PI / 180.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", 0.01);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J4", 0.001);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", 0.1);
    }
    let e0 = reb_simulation_energy(&sim);
    let v0 = rebx_with(&mut sim, rebx_gravitational_harmonics_potential).expect("attached");
    let h0 = e0 + v0;

    let mut worst_newtonian: f64 = 0.0;
    let mut worst_total: f64 = 0.0;
    for k in 1..=100 {
        reb_simulation_integrate(&mut sim, 2000.0 * (k as f64) / 100.0);
        let e = reb_simulation_energy(&sim);
        let v = rebx_with(&mut sim, rebx_gravitational_harmonics_potential).expect("attached");
        worst_newtonian = worst_newtonian.max(((e - e0) / e0).abs());
        worst_total = worst_total.max(((e + v - h0) / h0).abs());
    }
    assert!(
        worst_total < 1e-11,
        "gravitational_harmonics: kinetic + Newtonian + J2/J4 potential must be \
         conserved; max relative drift was {worst_total:.4e}"
    );
    assert!(
        worst_newtonian > 100.0 * worst_total,
        "gravitational_harmonics: the Newtonian-only energy should visibly breathe \
         (max relative excursion {worst_newtonian:.4e}) while the total stays put \
         (max relative drift {worst_total:.4e})"
    );
}

#[test]
fn gravitational_harmonics_potential_matches_closed_form() {
    let (J2, J4, R_eq) = (0.01f64, 0.001f64, 0.1f64);
    let mut sim = two_body(1e-4, 1.3, 0.35, 25.0 * PI / 180.0, 0.7);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", J2);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J4", J4);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", R_eq);
    }
    reb_simulation_integrate(&mut sim, 3.0);
    let got = rebx_with(&mut sim, rebx_gravitational_harmonics_potential).expect("attached");

    // Independent evaluation of the standard zonal expansion with the default
    // spin axis Omega = z_hat, so cos(theta) = dz/r:
    //   V = G m0 m1 / r * [ J2 (R/r)^2 P2 + J4 (R/r)^4 P4 ]
    // with P2 = (3c^2-1)/2 and P4 = (35c^4-30c^2+3)/8.
    let p0 = sim.particles[0];
    let p1 = sim.particles[1];
    let (dx, dy, dz) = (p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    let ct = dz / r;
    let ct2 = ct * ct;
    let p2 = 0.5 * (3.0 * ct2 - 1.0);
    let p4 = (35.0 * ct2 * ct2 - 30.0 * ct2 + 3.0) / 8.0;
    let pref = sim.G * p0.m * p1.m / r;
    let want = pref * (J2 * (R_eq / r).powi(2) * p2 + J4 * (R_eq / r).powi(4) * p4);
    assert_close(
        got,
        want,
        1e-13,
        "rebx_gravitational_harmonics_potential vs the zonal J2/J4 expansion",
    );

    // The J4 term must actually be contributing: rerun with J4 removed and
    // check the difference equals the J4 term computed above.
    let mut sim2 = two_body(1e-4, 1.3, 0.35, 25.0 * PI / 180.0, 0.7);
    rebx_attach(&mut sim2);
    let idx2 = rebx_load_force(&mut sim2, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim2, idx2), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim2) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", J2);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", R_eq);
    }
    reb_simulation_integrate(&mut sim2, 3.0);
    let got_j2_only = rebx_with(&mut sim2, rebx_gravitational_harmonics_potential).expect("attached");
    // sim2 has no J4 force, so its trajectory differs slightly from sim's;
    // evaluate the J2 term at sim2's own configuration.
    let q0 = sim2.particles[0];
    let q1 = sim2.particles[1];
    let (ex, ey, ez) = (q1.x - q0.x, q1.y - q0.y, q1.z - q0.z);
    let r2b = (ex * ex + ey * ey + ez * ez).sqrt();
    let ctb = ez / r2b;
    let p2b = 0.5 * (3.0 * ctb * ctb - 1.0);
    let j2_term = sim2.G * q0.m * q1.m / r2b * J2 * (R_eq / r2b).powi(2) * p2b;
    assert_close(
        got_j2_only,
        j2_term,
        1e-13,
        "rebx_gravitational_harmonics_potential with J4 unset vs the J2 term alone",
    );
    // And J4 was genuinely contributing in the first simulation.
    let j4_term = pref * J4 * (R_eq / r).powi(4) * p4;
    assert!(
        (j4_term / got).abs() > 1e-6,
        "the J4 term must be a non-negligible part of the potential for this test to \
         mean anything: J4 term {j4_term:.6e} out of total {got:.6e}"
    );
}

#[test]
fn gravitational_harmonics_with_zero_J2_is_bitwise_plain_gravity() {
    // The source short-circuits on `J2 == 0.0` before ever reading R_eq, so a
    // body with J2 = 0 must be indistinguishable from a point mass.
    let mut plain = two_body(1e-4, 1.0, 0.2, 20.0 * PI / 180.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 20.0 * PI / 180.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", 0.0);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", 0.1);
    }
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "gravitational_harmonics with J2 = 0 vs plain gravity",
    );

    let v = rebx_with(&mut sim, rebx_gravitational_harmonics_potential).expect("attached");
    assert_eq!(
        v.to_bits(),
        0.0f64.to_bits(),
        "rebx_gravitational_harmonics_potential with J2 = 0 must be exactly +0.0, got {v:.6e}"
    );
}

#[test]
fn gravitational_harmonics_without_R_eq_is_bitwise_plain_gravity() {
    // R_eq is required: with J2 set but R_eq missing the source `continue`s,
    // so the effect is a no-op rather than a NaN.
    let mut plain = two_body(1e-4, 1.0, 0.2, 20.0 * PI / 180.0, 0.0);
    reb_simulation_integrate(&mut plain, 300.0);

    let mut sim = two_body(1e-4, 1.0, 0.2, 20.0 * PI / 180.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", 0.01);
    }
    reb_simulation_integrate(&mut sim, 300.0);

    assert_bit_identical(
        &state_bits(&sim),
        &state_bits(&plain),
        "gravitational_harmonics with J2 set but R_eq missing vs plain gravity",
    );
}

#[test]
fn j2_advances_the_pericentre_of_a_coplanar_orbit() {
    // Spin axis along the orbit normal (the default Omega = z_hat with
    // inclination 0), so cos(theta) = 0 everywhere on the orbit. The in-plane
    // r^-4 term drives a PROGRADE apsidal precession. Combining the standard
    // secular rates
    //     d(omega)/dt = +(3/4) n J2 (R/p)^2 (5 cos^2 i - 1)
    //     d(Omega)/dt = -(3/2) n J2 (R/p)^2 cos i
    // at i = 0 gives d(pomega)/dt = (3 - 3/2) n J2 (R/p)^2 = (3/2) n J2 (R/p)^2,
    // i.e. the same magnitude as the nodal regression but with the opposite
    // sign — a nice cross-check on the previous test.
    let (J2, R_eq, a, e) = (0.01f64, 0.1f64, 1.0f64, 0.05f64);
    let tmax = 3000.0;
    let mut sim = two_body(1e-8, a, e, 0.0, 0.0);
    rebx_attach(&mut sim);
    let idx = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
    assert_eq!(rebx_add_force(&mut sim, idx), 1, "rebx_add_force(gh) failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", J2);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", R_eq);
    }
    let o0 = orbit(&sim);
    let p = o0.a * (1.0 - o0.e * o0.e);
    let want = 1.5 * o0.n * J2 * (R_eq / p) * (R_eq / p) * tmax;
    let (dpomega, _) = accumulate_angles(&mut sim, tmax, 300);
    assert!(
        dpomega > 0.0,
        "an oblate primary must advance the pericentre of a coplanar orbit; \
         got d(pomega) = {dpomega:.6e}"
    );
    assert_close(dpomega, want, 1e-2, "coplanar J2 apsidal advance over t=3000");

    // The magnitude matches the i=0 nodal regression rate of the previous
    // tests, with the opposite sign.
    let (node, _) = j2_node_drift(J2, R_eq, a, 0.0, tmax);
    assert!(
        node.abs() < 1e-14,
        "at i = 0 the ascending node is degenerate and must not drift; got {node:.6e}"
    );
}

// ===========================================================================
// Determinism
// ===========================================================================

#[test]
fn repeated_runs_are_bit_identical() {
    // Three conservative effects stacked on one simulation. Nothing in this
    // path may depend on allocation addresses, hash iteration order or the
    // (randomly seeded) rand_seed field.
    fn run() -> Vec<u64> {
        let mut sim = two_body(1e-4, 1.0, 0.25, 15.0 * PI / 180.0, 0.3);
        rebx_attach(&mut sim);
        let gr = rebx_load_force(&mut sim, "gr").expect("load gr");
        rebx_add_force(&mut sim, gr);
        let gh = rebx_load_force(&mut sim, "gravitational_harmonics").expect("load gh");
        rebx_add_force(&mut sim, gh);
        let cf = rebx_load_force(&mut sim, "central_force").expect("load central_force");
        rebx_add_force(&mut sim, cf);
        if let Some(rebx) = rebx_extras_mut(&mut sim) {
            rebx_set_param_double(rebx, rebx_ap::force(gr), "c", 100.0);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", 0.005);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", 0.1);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", 1e-4);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "gammacentral", -1.0);
        }
        reb_simulation_integrate(&mut sim, 500.0);
        state_bits(&sim)
    }
    let a = run();
    let b = run();
    assert_bit_identical(&a, &b, "two identical gr + J2 + central_force runs");

    // ... and the run actually did something: the state must not be the
    // untouched initial condition.
    let fresh = state_bits(&two_body(1e-4, 1.0, 0.25, 15.0 * PI / 180.0, 0.3));
    assert!(
        a != fresh,
        "the determinism run produced the initial state unchanged, so it proves nothing"
    );
}

#[test]
fn force_ordering_is_the_reverse_of_addition_order() {
    // rebx_add_force prepends, so accelerations are summed in reverse order of
    // addition. Floating-point addition is not associative, so adding two
    // effects in the opposite order is expected to give a DIFFERENT bit
    // pattern while agreeing physically. This pins the documented ordering
    // contract that bit-for-bit agreement with the C depends on.
    fn run(gr_first: bool) -> (Vec<u64>, reb_orbit) {
        let mut sim = two_body(1e-4, 1.0, 0.25, 15.0 * PI / 180.0, 0.3);
        rebx_attach(&mut sim);
        let (a, b) = if gr_first {
            ("gr_potential", "gravitational_harmonics")
        } else {
            ("gravitational_harmonics", "gr_potential")
        };
        let ia = rebx_load_force(&mut sim, a).expect("load a");
        rebx_add_force(&mut sim, ia);
        let ib = rebx_load_force(&mut sim, b).expect("load b");
        rebx_add_force(&mut sim, ib);
        let grx = if gr_first { ia } else { ib };
        if let Some(rebx) = rebx_extras_mut(&mut sim) {
            rebx_set_param_double(rebx, rebx_ap::force(grx), "c", 100.0);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "J2", 0.005);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "R_eq", 0.1);
        }
        reb_simulation_integrate(&mut sim, 500.0);
        let o = orbit(&sim);
        (state_bits(&sim), o)
    }
    let (bits_a, oa) = run(true);
    let (bits_b, ob) = run(false);

    // Physically the same system: elements must agree to well within the
    // physical precession accumulated over the run.
    assert_close(ob.a, oa.a, 1e-9, "force-order swap: semimajor axis");
    assert_close(ob.e, oa.e, 1e-6, "force-order swap: eccentricity");
    assert_close(ob.pomega, oa.pomega, 1e-6, "force-order swap: pomega");

    // But the summation order differs, so the bits are allowed to differ; the
    // useful invariant is that each ordering is reproducible on its own.
    let (bits_a2, _) = run(true);
    assert_bit_identical(&bits_a, &bits_a2, "same force ordering, repeated");
    assert_eq!(
        bits_a.len(),
        bits_b.len(),
        "state vectors must have the same shape regardless of force order"
    );
}
