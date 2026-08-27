//! Integration tests for the integrators_symplectic module group of rebound_rs.
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

// =====================================================================
// helpers
// =====================================================================

/// A zeroed particle, used as the base of struct-update literals.
const P0: reb_particle = reb_particle {
    x: 0.,
    y: 0.,
    z: 0.,
    vx: 0.,
    vy: 0.,
    vz: 0.,
    ax: 0.,
    ay: 0.,
    az: 0.,
    m: 0.,
    r: 0.,
    name: None,
};

fn particle(m: f64, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> reb_particle {
    reb_particle { m, x, y, z, vx, vy, vz, ..P0 }
}

fn norm(p: &reb_particle) -> f64 {
    (p.x * p.x + p.y * p.y + p.z * p.z).sqrt()
}

fn vnorm(p: &reb_particle) -> f64 {
    (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz).sqrt()
}

fn dist(a: &reb_particle, b: &reb_particle) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn vdist(a: &reb_particle, b: &reb_particle) -> f64 {
    let dx = a.vx - b.vx;
    let dy = a.vy - b.vy;
    let dz = a.vz - b.vz;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Difference of two particles, i.e. the relative coordinate.
fn relative(a: &reb_particle, b: &reb_particle) -> reb_particle {
    particle(
        a.m,
        a.x - b.x,
        a.y - b.y,
        a.z - b.z,
        a.vx - b.vx,
        a.vy - b.vy,
        a.vz - b.vz,
    )
}

/// Advance one particle along a Keplerian orbit of gravitational
/// parameter `mu` for a time `dt`, driving WHFast's universal-variable
/// (Stumpff / Stiefel) solver directly. The focus sits at the origin.
fn kepler_advance(p: reb_particle, mu: f64, dt: f64) -> reb_particle {
    let mut p_jh = [P0, p];
    let mut p_var: [reb_particle; 0] = [];
    integrator_whfast::reb_integrator_whfast_kepler_solver(None, &mut p_jh, &mut p_var, 1, mu, dt);
    p_jh[1]
}

/// Two-body specific orbital energy, v^2/2 - mu/r.
fn kepler_energy(p: &reb_particle, mu: f64) -> f64 {
    let v2 = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
    0.5 * v2 - mu / norm(p)
}

/// Two-body specific angular momentum r x v.
fn kepler_h(p: &reb_particle) -> reb_vec3d {
    reb_vec3d {
        x: p.y * p.vz - p.z * p.vy,
        y: p.z * p.vx - p.x * p.vz,
        z: p.x * p.vy - p.y * p.vx,
    }
}

fn vec3d_norm(v: &reb_vec3d) -> f64 {
    (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

fn vec3d_dist(a: &reb_vec3d, b: &reb_vec3d) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Independent Keplerian propagation, used as the reference route in
/// the Kepler-solver tests. It never touches the universal-variable
/// machinery: it converts the true anomaly to a mean anomaly in closed
/// form, advances the mean anomaly analytically, and rebuilds the state
/// vector via `reb_particle_from_orbit` (which solves Kepler's equation
/// with Newton's method in `reb_M_to_E`).
fn reference_advance(
    G: f64,
    primary: reb_particle,
    m: f64,
    a: f64,
    e: f64,
    inc: f64,
    Omega: f64,
    omega: f64,
    f0: f64,
    dt: f64,
) -> reb_particle {
    let mu = G * (m + primary.m);
    let (M0, n) = if e < 1. {
        let E0 = 2. * (((1. - e) / (1. + e)).sqrt() * (0.5 * f0).tan()).atan();
        (E0 - e * E0.sin(), (mu / (a * a * a)).sqrt())
    } else {
        let H0 = 2. * (((e - 1.) / (e + 1.)).sqrt() * (0.5 * f0).tan()).atanh();
        let q = -a;
        (e * H0.sinh() - H0, (mu / (q * q * q)).sqrt())
    };
    let M1 = M0 + n * dt;
    let f1 = reb_M_to_f(e, M1);
    reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f1)
}

/// Bit-level comparison of the dynamical state of two particles.
fn same_bits(a: &reb_particle, b: &reb_particle) -> bool {
    a.x.to_bits() == b.x.to_bits()
        && a.y.to_bits() == b.y.to_bits()
        && a.z.to_bits() == b.z.to_bits()
        && a.vx.to_bits() == b.vx.to_bits()
        && a.vy.to_bits() == b.vy.to_bits()
        && a.vz.to_bits() == b.vz.to_bits()
        && a.m.to_bits() == b.m.to_bits()
}

fn states_same_bits(a: &[reb_particle], b: &[reb_particle]) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| same_bits(x, y))
}

/// A deterministic, mildly hierarchical five-body configuration used by
/// the coordinate-transformation tests. Nothing here is symmetric, so a
/// transposed index or a dropped term shows up immediately.
fn five_body() -> Vec<reb_particle> {
    vec![
        particle(1.0, 0.13, -0.07, 0.02, -0.011, 0.023, -0.004),
        particle(1e-3, 1.31, 0.44, -0.09, -0.31, 0.87, 0.05),
        particle(4e-4, -2.05, 1.17, 0.33, -0.44, -0.62, 0.017),
        particle(7e-5, 3.41, -2.62, -0.55, 0.39, 0.41, -0.11),
        particle(2e-6, -5.13, -0.88, 1.04, 0.13, -0.36, 0.07),
    ]
}

/// Center of mass of `p[first..last]`, computed straightforwardly.
fn com_of(p: &[reb_particle]) -> reb_particle {
    let mut m = 0.;
    let mut c = P0;
    for q in p {
        c.x += q.m * q.x;
        c.y += q.m * q.y;
        c.z += q.m * q.z;
        c.vx += q.m * q.vx;
        c.vy += q.m * q.vy;
        c.vz += q.m * q.vz;
        m += q.m;
    }
    c.x /= m;
    c.y /= m;
    c.z /= m;
    c.vx /= m;
    c.vy /= m;
    c.vz /= m;
    c.m = m;
    c
}

/// Scramble positions and velocities so that a transform that forgets
/// to write an entry cannot silently pass a round-trip test.
fn scrambled(src: &[reb_particle]) -> Vec<reb_particle> {
    src.iter()
        .enumerate()
        .map(|(i, p)| reb_particle {
            x: 1e3 + i as f64,
            y: -2e3 - i as f64,
            z: 3e3 + i as f64,
            vx: -4e3 - i as f64,
            vy: 5e3 + i as f64,
            vz: -6e3 - i as f64,
            ..*p
        })
        .collect()
}

/// Star + two planets, built from orbital elements. `integrator` is set
/// before the particles are added so that the integrator's
/// `did_add_particle` hook (if any) runs the same way in every test.
fn three_body_sim(integrator: &str, dt: f64) -> reb_simulation {
    let mut r = reb_simulation_create();
    r.G = 1.;
    r.save_messages = 1;
    reb_simulation_set_integrator(&mut r, integrator);
    let star = reb_particle { m: 1., ..P0 };
    reb_simulation_add(&mut r, star);
    let p1 = reb_particle_from_orbit(r.G, star, 1e-3, 1.0, 0.05, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, p1);
    let com = reb_simulation_com(&r);
    let p2 = reb_particle_from_orbit(r.G, com, 2e-3, 2.35, 0.12, 0.03, 0.4, 0.9, 1.3);
    reb_simulation_add(&mut r, p2);
    reb_simulation_move_to_com(&mut r);
    r.dt = dt;
    r
}

/// Star + one planet on an orbit of semi-major axis `a`, eccentricity
/// `e`, started at periapsis in the x-y plane.
fn two_body_sim(integrator: &str, m1: f64, a: f64, e: f64, dt: f64) -> reb_simulation {
    let mut r = reb_simulation_create();
    r.G = 1.;
    r.save_messages = 1;
    reb_simulation_set_integrator(&mut r, integrator);
    let star = reb_particle { m: 1., ..P0 };
    reb_simulation_add(&mut r, star);
    let p = reb_particle_from_orbit(r.G, star, m1, a, e, 0., 0., 0., 0.);
    reb_simulation_add(&mut r, p);
    reb_simulation_move_to_com(&mut r);
    r.dt = dt;
    r
}

/// Relative energy drift of `r` over `n_steps` fixed timesteps.
fn energy_drift(r: &mut reb_simulation, n_steps: usize) -> f64 {
    let e0 = reb_simulation_energy(r);
    reb_simulation_steps(r, n_steps);
    let e1 = reb_simulation_energy(r);
    ((e1 - e0) / e0).abs()
}

fn has_message(r: &reb_simulation, needle: &str) -> bool {
    r.messages.iter().any(|(_, m)| m.contains(needle))
}

fn has_error(r: &reb_simulation) -> bool {
    r.messages
        .iter()
        .any(|(t, _)| *t == REB_MESSAGE_TYPE::ERROR)
}

// =====================================================================
// Stumpff / Stiefel series and the Kepler solver
// =====================================================================

/// For a circular orbit of unit radius about mu = 1 the universal
/// anomaly is exactly X = dt (eta0 = x.v = 0 and zeta0 = mu - beta*r0 =
/// 0 make the very first Newton correction a fixed point), so the
/// solver's f and g reduce to
///     1 + f = 1 - G2 = cos X,   g = dt - G3 = sin X,
///     fd    = -G1    = -sin X,  1 + gd = 1 - G2 = cos X.
/// The Kepler step therefore has to be a plane rotation by dt, which
/// pins the Stumpff series c0..c3 against the standard trig functions.
/// dt values above 0.316 push |z| = dt^2 past 0.1 and so also exercise
/// the argument-halving/doubling loop in `stumpff_cs3`.
#[test]
fn stumpff_series_makes_circular_kepler_step_a_rotation() {
    let mu = 1.0;
    let start = particle(0., 1., 0., 0., 0., 1., 0.);
    for &dt in &[
        -6.0_f64, -3.0, -1.0, -0.3, -0.05, 0.05, 0.3, 1.0, 2.0, 3.0, 6.0,
    ] {
        let p = kepler_advance(start, mu, dt);
        let (c, s) = (dt.cos(), dt.sin());
        // g = dt - G3 subtracts two quantities of size |dt|, so the
        // achievable absolute accuracy grows linearly with |dt|.
        let tol = 2e-14 * (1. + dt.abs());
        assert!(
            (p.x - c).abs() < tol,
            "circular Kepler step x != cos(dt) for dt={}: got {}, want {}",
            dt,
            p.x,
            c
        );
        assert!(
            (p.y - s).abs() < tol,
            "circular Kepler step y != sin(dt) for dt={}: got {}, want {}",
            dt,
            p.y,
            s
        );
        assert!(
            p.z == 0.,
            "circular Kepler step left the orbital plane for dt={}: z={}",
            dt,
            p.z
        );
        assert!(
            (p.vx + s).abs() < tol,
            "circular Kepler step vx != -sin(dt) for dt={}: got {}, want {}",
            dt,
            p.vx,
            -s
        );
        assert!(
            (p.vy - c).abs() < tol,
            "circular Kepler step vy != cos(dt) for dt={}: got {}, want {}",
            dt,
            p.vy,
            c
        );
    }
    // Small steps stay within a couple of ulps of the trig functions.
    for &dt in &[-1.0_f64, -0.25, 0.25, 1.0] {
        let p = kepler_advance(start, mu, dt);
        assert!(
            (p.x - dt.cos()).abs() < 4e-16 && (p.y - dt.sin()).abs() < 4e-16,
            "circular Kepler step is not accurate to a few ulps for dt={}: ({}, {}) vs ({}, {})",
            dt,
            p.x,
            p.y,
            dt.cos(),
            dt.sin()
        );
    }
}

/// The reference route used throughout this file has to be validated on
/// its own terms: `reb_M_to_E` must return a root of Kepler's equation.
#[test]
fn M_to_E_solves_keplers_equation() {
    for &e in &[0.0, 0.01, 0.3, 0.7, 0.95, 0.999] {
        for k in 0..24 {
            let M = -3.0 * PI + 7.0 * PI * (k as f64) / 23.0;
            let E = reb_M_to_E(e, M);
            let residual = E - e * E.sin() - reb_mod2pi(M);
            // The root may be reported one revolution away from mod2pi(M).
            let residual = residual - (residual / (2. * PI)).round() * 2. * PI;
            assert!(
                residual.abs() < 1e-9,
                "reb_M_to_E(e={}, M={}) = {} does not satisfy Kepler's equation (residual {})",
                e,
                M,
                E,
                residual
            );
        }
    }
    // Hyperbolic branch: e sinh(H) - H = M.
    for &e in &[1.2, 2.0, 5.0] {
        for k in 1..12 {
            let M = 0.25 * (k as f64);
            for &sgn in &[-1.0_f64, 1.0] {
                let H = reb_M_to_E(e, sgn * M);
                let residual = e * H.sinh() - H - sgn * M;
                assert!(
                    residual.abs() < 1e-9,
                    "reb_M_to_E(e={}, M={}) = {} does not satisfy the hyperbolic Kepler equation (residual {})",
                    e,
                    sgn * M,
                    H,
                    residual
                );
            }
        }
    }
}

/// `reb_E_to_f` must satisfy the conic relation cos f = (cos E - e) /
/// (1 - e cos E), which is an entirely different formula from the
/// half-angle tangent it actually evaluates.
#[test]
fn E_to_f_matches_the_conic_relation() {
    for &e in &[0.0, 0.2, 0.6, 0.9] {
        for k in 0..17 {
            let E = -2.9 + 5.8 * (k as f64) / 16.0;
            let f = reb_E_to_f(e, E);
            let want = (E.cos() - e) / (1. - e * E.cos());
            assert!(
                (f.cos() - want).abs() < 1e-12,
                "reb_E_to_f(e={}, E={}) = {}: cos f = {} but conic relation gives {}",
                e,
                E,
                f,
                f.cos(),
                want
            );
        }
    }
}

/// `reb_mod2pi` maps into [0, 2pi) and is the identity modulo 2pi.
#[test]
fn mod2pi_maps_into_the_canonical_interval() {
    assert!(
        reb_mod2pi(0.0) == 0.0,
        "reb_mod2pi(0) should be exactly 0, got {}",
        reb_mod2pi(0.0)
    );
    for &f in &[
        0.0,
        1.0,
        PI,
        2. * PI,
        -2. * PI,
        4. * PI,
        -1.0,
        -7.5,
        123.456,
        -123.456,
    ] {
        let m = reb_mod2pi(f);
        assert!(
            (0.0..2. * PI).contains(&m),
            "reb_mod2pi({}) = {} is outside [0, 2pi)",
            f,
            m
        );
        let k = ((f - m) / (2. * PI)).round();
        assert!(
            (f - m - k * 2. * PI).abs() < 1e-12 * (1. + f.abs()),
            "reb_mod2pi({}) = {} differs from the input by a non-multiple of 2pi",
            f,
            m
        );
    }
    // Exact multiples of 2pi collapse to zero.
    for k in -4..=4 {
        let m = reb_mod2pi(2. * PI * (k as f64));
        assert!(
            m == 0.0,
            "reb_mod2pi({}*2pi) should be exactly 0, got {}",
            k,
            m
        );
    }
}

/// The universal-variable solver has to agree with the classical
/// element route over a wide range of eccentricity and step size. The
/// small steps exercise Newton's method; the steps that are a sizable
/// fraction of the period make the second-order initial guess bad
/// enough that the solver switches to the quartic (Laguerre) branch.
#[test]
fn kepler_solver_matches_element_propagation_elliptic() {
    let G = 1.0;
    let mstar = 1.0;
    let m = 1e-3;
    let primary = reb_particle { m: mstar, ..P0 };
    let mu = G * (mstar + m);
    let a = 1.7;
    let period = 2. * PI * (a * a * a / mu).sqrt();
    let (inc, Omega, omega, f0) = (0.31, 0.77, 1.23, 0.62);
    for &e in &[0.0, 0.1, 0.5, 0.9] {
        let p0 = reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f0);
        // Sanity gate on the reference route itself: at dt = 0 it must
        // reproduce the initial state.
        let r0 = reference_advance(G, primary, m, a, e, inc, Omega, omega, f0, 0.0);
        assert!(
            dist(&p0, &r0) < 1e-12 * a,
            "reference route is inconsistent at dt=0 for e={}: offset {}",
            e,
            dist(&p0, &r0)
        );
        for &frac in &[0.005, 0.05, 0.25, 0.4, 0.6, 0.9, -0.25, -0.6] {
            let dt = frac * period;
            let got = kepler_advance(p0, mu, dt);
            let want = reference_advance(G, primary, m, a, e, inc, Omega, omega, f0, dt);
            let scale = norm(&want).max(a);
            assert!(
                dist(&got, &want) < 1e-11 * scale,
                "Kepler solver position disagrees with element propagation (e={}, dt={} = {} P): |dr| = {}",
                e,
                dt,
                frac,
                dist(&got, &want)
            );
            let vscale = vnorm(&want).max((mu / a).sqrt());
            assert!(
                vdist(&got, &want) < 1e-11 * vscale,
                "Kepler solver velocity disagrees with element propagation (e={}, dt={} = {} P): |dv| = {}",
                e,
                dt,
                frac,
                vdist(&got, &want)
            );
        }
    }
}

/// Same cross-check on hyperbolic orbits (beta < 0), which take the
/// hyperbolic initial guess X = 0 and never see `X_per_period`.
#[test]
fn kepler_solver_matches_element_propagation_hyperbolic() {
    let G = 1.0;
    let mstar = 1.0;
    let m = 0.0;
    let primary = reb_particle { m: mstar, ..P0 };
    let mu = G * (mstar + m);
    let a = -2.0;
    let (inc, Omega, omega, f0) = (0.44, 1.1, 2.2, 0.0);
    for &e in &[1.2, 2.0, 5.0] {
        let p0 = reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f0);
        assert!(
            p0.x.is_finite() && p0.vx.is_finite(),
            "hyperbolic initial condition is not finite for e={}",
            e
        );
        for &dt in &[-3.0, -1.0, -0.2, 0.2, 1.0, 3.0] {
            let got = kepler_advance(p0, mu, dt);
            let want = reference_advance(G, primary, m, a, e, inc, Omega, omega, f0, dt);
            let scale = norm(&want).max(1.0);
            assert!(
                dist(&got, &want) < 1e-10 * scale,
                "hyperbolic Kepler solver disagrees with element propagation (e={}, dt={}): |dr| = {}",
                e,
                dt,
                dist(&got, &want)
            );
            assert!(
                vdist(&got, &want) < 1e-10 * vnorm(&want).max(1.0),
                "hyperbolic Kepler solver velocity disagrees with element propagation (e={}, dt={}): |dv| = {}",
                e,
                dt,
                vdist(&got, &want)
            );
        }
    }
}

/// A dt of exactly zero drives the solver down its slowest path: the
/// first Newton correction returns X = 0, which is not a normal float,
/// so the Newton loop breaks without setting `converged` and the
/// bisection fallback runs. Bisection has to terminate and to leave the
/// particle untouched, bit for bit.
#[test]
fn kepler_solver_zero_timestep_is_bitwise_identity() {
    let G = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    // Elliptic (beta > 0) and hyperbolic (beta < 0) fallback branches.
    let cases = [
        reb_particle_from_orbit(G, primary, 1e-3, 1.4, 0.35, 0.2, 0.5, 0.9, 1.1),
        reb_particle_from_orbit(G, primary, 1e-3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        reb_particle_from_orbit(G, primary, 0.0, -2.0, 1.6, 0.3, 0.1, 0.2, 0.4),
    ];
    for (i, p0) in cases.iter().enumerate() {
        let mu = G * (1. + p0.m);
        let p1 = kepler_advance(*p0, mu, 0.0);
        assert!(
            same_bits(p0, &p1),
            "Kepler step with dt=0 changed case {}: r ({}, {}, {}) -> ({}, {}, {})",
            i,
            p0.x,
            p0.y,
            p0.z,
            p1.x,
            p1.y,
            p1.z
        );
    }
}

/// The Keplerian flow is a one-parameter group: phi_dt must equal
/// phi_{dt/n} applied n times. This is an internal consistency check
/// that needs no external reference and stays sharp at eccentricities
/// where solving Kepler's equation by Newton iteration is delicate.
#[test]
fn kepler_solver_flow_obeys_the_group_property() {
    let G = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    let mu = G;
    let orbits = [
        (1.0, 0.0),
        (1.3, 0.5),
        (1.3, 0.99),
        (1.3, 0.999),
        (-2.0, 1.4),
        (-2.0, 4.0),
    ];
    for &(a, e) in &orbits {
        let f0 = if e < 1. { 0.9 } else { 0.0 };
        let p0 = reb_particle_from_orbit(G, primary, 0.0, a, e, 0.21, 0.32, 0.43, f0);
        let period = if a > 0. {
            2. * PI * (a * a * a / mu).sqrt()
        } else {
            1.0
        };
        for &dt in &[0.13 * period, -0.37 * period] {
            let one = kepler_advance(p0, mu, dt);
            let mut many = p0;
            for _ in 0..8 {
                many = kepler_advance(many, mu, dt / 8.);
            }
            let scale = norm(&one).max(norm(&p0));
            assert!(
                dist(&one, &many) < 1e-10 * scale,
                "Kepler flow is not a group for (a={}, e={}, dt={}): one step vs eight gives |dr| = {}",
                a,
                e,
                dt,
                dist(&one, &many)
            );
            assert!(
                vdist(&one, &many) < 1e-10 * vnorm(&one).max(1e-3),
                "Kepler flow velocities disagree for (a={}, e={}, dt={}): |dv| = {}",
                a,
                e,
                dt,
                vdist(&one, &many)
            );
        }
    }
}

/// Stepping forward by dt and back by -dt has to return the original
/// state (the exact flow is invertible; only the iterative solve is
/// inexact).
#[test]
fn kepler_solver_is_time_reversible() {
    let G = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    let mu = G;
    for &(a, e) in &[(1.0, 0.0), (1.5, 0.7), (1.5, 0.98), (-3.0, 2.5)] {
        let f0 = if e < 1. { 1.4 } else { 0.0 };
        let p0 = reb_particle_from_orbit(G, primary, 0.0, a, e, 0.1, 0.2, 0.3, f0);
        for &dt in &[0.05, 0.7, 2.9] {
            let fwd = kepler_advance(p0, mu, dt);
            let back = kepler_advance(fwd, mu, -dt);
            let scale = norm(&p0);
            assert!(
                dist(&p0, &back) < 1e-11 * scale,
                "Kepler step is not reversible for (a={}, e={}, dt={}): |dr| = {}",
                a,
                e,
                dt,
                dist(&p0, &back)
            );
            assert!(
                vdist(&p0, &back) < 1e-11 * vnorm(&p0),
                "Kepler velocities are not reversible for (a={}, e={}, dt={}): |dv| = {}",
                a,
                e,
                dt,
                vdist(&p0, &back)
            );
        }
    }
}

/// A step of exactly one orbital period is the identity map on a bound
/// orbit.
#[test]
fn kepler_solver_full_period_returns_to_start() {
    let G: f64 = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    for &e in &[0.0_f64, 0.3, 0.85] {
        for &a in &[0.4_f64, 1.0, 3.7] {
            let mu = G;
            let period = 2. * PI * (a * a * a / mu).sqrt();
            let p0 = reb_particle_from_orbit(G, primary, 0.0, a, e, 0.5, 1.0, 1.5, 0.8);
            let p1 = kepler_advance(p0, mu, period);
            assert!(
                dist(&p0, &p1) < 1e-11 * a,
                "one full period is not the identity (a={}, e={}): |dr| = {}",
                a,
                e,
                dist(&p0, &p1)
            );
        }
    }
}

/// The Kepler step is the exact two-body flow, so the specific orbital
/// energy and the angular momentum vector are invariants of it.
#[test]
fn kepler_solver_conserves_energy_and_angular_momentum() {
    let G = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    let mu = G;
    for &(a, e) in &[
        (1.0, 0.0),
        (1.0, 0.3),
        (1.0, 0.9),
        (1.0, 0.999),
        (-1.5, 1.3),
        (-1.5, 3.0),
    ] {
        let f0 = if e < 1. { 0.35 } else { 0.0 };
        let mut p = reb_particle_from_orbit(G, primary, 0.0, a, e, 0.6, 0.2, 1.9, f0);
        let e0 = kepler_energy(&p, mu);
        let h0 = kepler_h(&p);
        for k in 0..60 {
            let dt = 0.017 + 0.003 * (k as f64);
            p = kepler_advance(p, mu, dt);
            let e1 = kepler_energy(&p, mu);
            assert!(
                ((e1 - e0) / e0).abs() < 1e-11,
                "Kepler step does not conserve energy (a={}, e={}, step {}): E {} -> {}",
                a,
                e,
                k,
                e0,
                e1
            );
            let h1 = kepler_h(&p);
            assert!(
                vec3d_dist(&h0, &h1) < 1e-11 * vec3d_norm(&h0),
                "Kepler step does not conserve angular momentum (a={}, e={}, step {}): |dh| = {}",
                a,
                e,
                k,
                vec3d_dist(&h0, &h1)
            );
        }
    }
}

/// Exactly rectilinear (zero angular momentum) hyperbolic motion. This
/// is the degenerate case the solver guards with its `ri.is_nan()`
/// branch: h = 0 makes the periapsis distance used by the bisection
/// fallback zero, so `X_min` becomes NaN and `X_max` becomes +inf.
/// For steps the Newton iteration can handle, the answer is still the
/// exact rectilinear Kepler solution.
#[test]
fn kepler_solver_handles_exactly_rectilinear_motion() {
    let mu = 1.0;
    // v^2 = 9 > 2 mu / r = 2 (unbound), and r x v = 0 identically.
    let start = particle(0., 1., 0., 0., 3., 0., 0.);
    let e0 = kepler_energy(&start, mu);
    for &dt in &[0.05_f64, 0.1, 0.25, 0.5, 1.0, -0.1] {
        let p = kepler_advance(start, mu, dt);
        assert!(
            p.y.to_bits() == 0u64
                && p.z.to_bits() == 0u64
                && p.vy.to_bits() == 0u64
                && p.vz.to_bits() == 0u64,
            "rectilinear motion left the x axis for dt={}: y={}, z={}, vy={}, vz={}",
            dt,
            p.y,
            p.z,
            p.vy,
            p.vz
        );
        let e1 = kepler_energy(&p, mu);
        assert!(
            ((e1 - e0) / e0).abs() < 1e-14,
            "rectilinear Kepler step does not conserve energy for dt={}: {} -> {}",
            dt,
            e0,
            e1
        );
        // Outbound motion decelerates but never reverses on a hyperbola.
        if dt > 0. {
            assert!(
                p.x > start.x && p.vx < start.vx && p.vx > 0.,
                "outbound rectilinear motion is wrong for dt={}: x {} -> {}, vx {} -> {}",
                dt,
                start.x,
                p.x,
                start.vx,
                p.vx
            );
        }
    }
    // The rectilinear flow is a group just like a curved one.
    let one = kepler_advance(start, mu, 1.0);
    let mut many = start;
    for _ in 0..4 {
        many = kepler_advance(many, mu, 0.25);
    }
    assert!(
        dist(&one, &many) < 1e-12 * norm(&one),
        "rectilinear Kepler flow is not a group: |dr| = {}",
        dist(&one, &many)
    );
}

/// When the angular momentum is exactly zero and the step is too large
/// for the Newton iteration, the solver runs its bisection fallback on
/// NaN bounds, exits after one pass, and lands in the C's documented
/// "Exception for (almost) straight line motion in hyperbolic case".
/// That exception zeroes ri and G1..G3, which leaves f = fd = gd = 0
/// and g = dt: the step degrades to exact free streaming. This pins
/// that escape hatch, which is unreachable for any h != 0.
#[test]
fn kepler_solver_straight_line_exception_is_exact_free_streaming() {
    let mu = 1.0;
    let start = particle(0., 1., 0., 0., 3., 0., 0.);
    let mut triggered = 0;
    for &dt in &[2.0_f64, 3.0, 5.0, -0.5, -1.0, -2.0] {
        let p = kepler_advance(start, mu, dt);
        let is_free = p.x.to_bits() == (start.x + dt * start.vx).to_bits()
            && p.vx.to_bits() == start.vx.to_bits();
        if !is_free {
            // The Newton iteration coped; then it must be a real Kepler
            // step, i.e. it must conserve energy.
            let rel = ((kepler_energy(&p, mu) - kepler_energy(&start, mu))
                / kepler_energy(&start, mu))
            .abs();
            assert!(
                rel < 1e-14,
                "dt={} took neither the straight-line exception nor an energy-conserving Kepler step (dE/E = {})",
                dt,
                rel
            );
            continue;
        }
        triggered += 1;
        assert!(
            p.y.to_bits() == 0u64 && p.z.to_bits() == 0u64,
            "the straight-line exception moved the particle off the x axis for dt={}",
            dt
        );
        assert!(
            p.vy.to_bits() == 0u64 && p.vz.to_bits() == 0u64,
            "the straight-line exception changed the transverse velocity for dt={}",
            dt
        );
    }
    assert!(
        triggered > 0,
        "none of the probed timesteps reached the straight-line exception, so that branch is untested"
    );
}

/// Near-rectilinear hyperbolic motion: the transverse velocity is only
/// 1.7% of the radial one, so |h| is far below |r||v| and the orbit is
/// an extremely eccentric hyperbola. The transverse velocity has to
/// stay above sqrt(eps)*|v| or `h2 = r0^2 v^2 - eta0^2` cancels to
/// exactly zero and the orbit is rectilinear after all.
#[test]
fn kepler_solver_handles_near_rectilinear_hyperbolic_motion() {
    let mu = 1.0;
    // v^2 > 2 mu / r = 2, so beta < 0; h = 0.05 is small but well above
    // the cancellation threshold.
    let p0 = particle(0., 1., 0., 0., 3., 0.05, 0.);
    let e0 = kepler_energy(&p0, mu);
    let h0 = kepler_h(&p0);
    let mut p = p0;
    for k in 0..20 {
        p = kepler_advance(p, mu, 0.1);
        assert!(
            p.z == 0. && p.vz == 0.,
            "planar motion left the z = 0 plane at step {}: z={}, vz={}",
            k,
            p.z,
            p.vz
        );
        assert!(
            p.x > 1.0,
            "outbound near-radial motion should move away from the focus, x = {} at step {}",
            p.x,
            k
        );
        let e1 = kepler_energy(&p, mu);
        assert!(
            ((e1 - e0) / e0).abs() < 1e-11,
            "near-rectilinear Kepler step does not conserve energy at step {}: {} -> {}",
            k,
            e0,
            e1
        );
        let h1 = kepler_h(&p);
        assert!(
            vec3d_dist(&h0, &h1) < 1e-9 * vec3d_norm(&h0),
            "near-rectilinear Kepler step does not conserve angular momentum at step {}: |dh| = {}",
            k,
            vec3d_dist(&h0, &h1)
        );
    }
    // And it composes with itself the same way a curved orbit does.
    let one = kepler_advance(p0, mu, 2.0);
    let mut many = p0;
    for _ in 0..4 {
        many = kepler_advance(many, mu, 0.5);
    }
    assert!(
        dist(&one, &many) < 1e-10 * norm(&one),
        "near-rectilinear Kepler flow is not a group: |dr| = {}",
        dist(&one, &many)
    );
}

/// The solver warns at most once about a timestep exceeding one orbital
/// period, and records the fact in `messages_timestep_warning`.
#[test]
fn kepler_solver_warns_once_about_oversized_timestep() {
    let mut r = reb_simulation_create();
    r.G = 1.;
    r.save_messages = 1;
    let primary = reb_particle { m: 1., ..P0 };
    let p0 = reb_particle_from_orbit(1., primary, 0., 1.0, 0.2, 0., 0., 0., 0.);
    let mut p_jh = [P0, p0];
    let mut p_var: [reb_particle; 0] = [];
    let period = 2. * PI;

    // A step well inside one period must not warn.
    integrator_whfast::reb_integrator_whfast_kepler_solver(
        Some(&mut r),
        &mut p_jh,
        &mut p_var,
        1,
        1.0,
        0.1 * period,
    );
    assert!(
        r.messages.is_empty(),
        "Kepler solver warned about a timestep of 0.1 P: {:?}",
        r.messages
    );
    assert!(
        r.messages_timestep_warning == 0,
        "messages_timestep_warning was set by a short timestep"
    );

    // A step longer than one period warns exactly once.
    integrator_whfast::reb_integrator_whfast_kepler_solver(
        Some(&mut r),
        &mut p_jh,
        &mut p_var,
        1,
        1.0,
        1.7 * period,
    );
    assert!(
        has_message(&r, "larger than one orbital period"),
        "Kepler solver did not warn about a timestep of 1.7 P: {:?}",
        r.messages
    );
    assert!(
        r.messages.len() == 1,
        "expected exactly one warning, got {:?}",
        r.messages
    );
    assert!(
        (r.messages_timestep_warning & 1) == 1,
        "messages_timestep_warning bit 0 was not set, value {}",
        r.messages_timestep_warning
    );

    integrator_whfast::reb_integrator_whfast_kepler_solver(
        Some(&mut r),
        &mut p_jh,
        &mut p_var,
        1,
        1.0,
        1.7 * period,
    );
    assert!(
        r.messages.len() == 1,
        "the oversized-timestep warning repeated: {:?}",
        r.messages
    );
}

/// A retrograde orbit (inc = pi) has a negative z angular momentum and
/// must still be propagated correctly.
#[test]
fn kepler_solver_handles_retrograde_orbits() {
    let G = 1.0;
    let primary = reb_particle { m: 1., ..P0 };
    let mu = G;
    let (a, e) = (1.6, 0.4);
    let p0 = reb_particle_from_orbit(G, primary, 0., a, e, PI, 0., 0., 0.9);
    let h = kepler_h(&p0);
    assert!(
        h.z < 0.,
        "an inc = pi orbit must have h_z < 0, got {}",
        h.z
    );
    let period = 2. * PI * (a * a * a / mu).sqrt();
    for &frac in &[0.03, 0.3, 0.75] {
        let dt = frac * period;
        let got = kepler_advance(p0, mu, dt);
        let want = reference_advance(G, primary, 0., a, e, PI, 0., 0., 0.9, dt);
        assert!(
            dist(&got, &want) < 1e-10 * a,
            "retrograde Kepler step disagrees with element propagation at dt = {} P: |dr| = {}",
            frac,
            dist(&got, &want)
        );
        assert!(
            kepler_h(&got).z < 0.,
            "retrograde orbit flipped orientation at dt = {} P",
            frac
        );
    }
}

// =====================================================================
// coordinate systems
// =====================================================================

#[test]
fn jacobi_transform_round_trips_to_the_original_state() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_j = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_posvel(&inertial, &mut p_j, &inertial, n, n);
    let mut back = scrambled(&inertial);
    reb_transformations_jacobi_to_inertial_posvel(&mut back, &p_j, &inertial, n, n);
    for i in 0..n {
        assert!(
            dist(&back[i], &inertial[i]) < 1e-14,
            "inertial -> Jacobi -> inertial changed the position of particle {}: |dr| = {}",
            i,
            dist(&back[i], &inertial[i])
        );
        assert!(
            vdist(&back[i], &inertial[i]) < 1e-14,
            "inertial -> Jacobi -> inertial changed the velocity of particle {}: |dv| = {}",
            i,
            vdist(&back[i], &inertial[i])
        );
    }
}

#[test]
fn jacobi_zeroth_coordinate_is_the_center_of_mass() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_j = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_posvel(&inertial, &mut p_j, &inertial, n, n);
    let com = com_of(&inertial);
    assert!(
        (p_j[0].m - com.m).abs() < 1e-15 * com.m,
        "Jacobi coordinate 0 should carry the total mass {}, got {}",
        com.m,
        p_j[0].m
    );
    assert!(
        dist(&p_j[0], &com) < 1e-15,
        "Jacobi coordinate 0 is not the center of mass: |dr| = {}",
        dist(&p_j[0], &com)
    );
    assert!(
        vdist(&p_j[0], &com) < 1e-15,
        "Jacobi coordinate 0 velocity is not the center-of-mass velocity: |dv| = {}",
        vdist(&p_j[0], &com)
    );
}

/// Jacobi coordinate i is the offset of particle i from the center of
/// mass of particles 0..i-1. That definition is computed here directly
/// and compared against the recurrence the transform actually uses.
#[test]
fn jacobi_coordinates_are_offsets_from_the_inner_center_of_mass() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_j = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_posvel(&inertial, &mut p_j, &inertial, n, n);
    for i in 1..n {
        let inner = com_of(&inertial[0..i]);
        let want = relative(&inertial[i], &inner);
        assert!(
            dist(&p_j[i], &want) < 1e-14,
            "Jacobi coordinate {} is not r_{} - com(0..{}): |dr| = {}",
            i,
            i,
            i - 1,
            dist(&p_j[i], &want)
        );
        assert!(
            vdist(&p_j[i], &want) < 1e-14,
            "Jacobi velocity {} is not v_{} - v_com(0..{}): |dv| = {}",
            i,
            i,
            i - 1,
            vdist(&p_j[i], &want)
        );
        assert!(
            p_j[i].m == inertial[i].m,
            "Jacobi coordinate {} must carry the particle's own mass",
            i
        );
    }
}

/// `inertial_to_jacobi_acc` performs exactly the same linear map as the
/// position half of `inertial_to_jacobi_posvel`. Feeding accelerations
/// that equal the positions must therefore produce Jacobi accelerations
/// that equal the Jacobi positions, bit for bit.
#[test]
fn jacobi_acceleration_transform_matches_the_position_transform() {
    let mut inertial = five_body();
    for p in inertial.iter_mut() {
        p.ax = p.x;
        p.ay = p.y;
        p.az = p.z;
    }
    let n = inertial.len();
    let mut p_pos = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_posvel(&inertial, &mut p_pos, &inertial, n, n);
    let mut p_acc = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_acc(&inertial, &mut p_acc, &inertial, n, n);
    for i in 0..n {
        assert!(
            p_acc[i].ax.to_bits() == p_pos[i].x.to_bits()
                && p_acc[i].ay.to_bits() == p_pos[i].y.to_bits()
                && p_acc[i].az.to_bits() == p_pos[i].z.to_bits(),
            "Jacobi acceleration transform differs from the position transform at {}: ({}, {}, {}) vs ({}, {}, {})",
            i,
            p_acc[i].ax,
            p_acc[i].ay,
            p_acc[i].az,
            p_pos[i].x,
            p_pos[i].y,
            p_pos[i].z
        );
    }
    // ... and the inverse map, too. Only the posvel transform writes
    // the Jacobi masses, and the inverse maps read p_j[0].m, so the
    // acceleration array is given the same masses it would carry inside
    // WHFast's p_jh buffer.
    for i in 0..n {
        p_acc[i].m = p_pos[i].m;
    }
    let mut back_pos = scrambled(&inertial);
    reb_transformations_jacobi_to_inertial_pos(&mut back_pos, &p_pos, &inertial, n, n);
    let mut back_acc = scrambled(&inertial);
    reb_transformations_jacobi_to_inertial_acc(&mut back_acc, &p_acc, &inertial, n, n);
    for i in 0..n {
        assert!(
            back_acc[i].ax.to_bits() == back_pos[i].x.to_bits()
                && back_acc[i].ay.to_bits() == back_pos[i].y.to_bits()
                && back_acc[i].az.to_bits() == back_pos[i].z.to_bits(),
            "Jacobi-to-inertial acceleration transform differs from the position transform at {}",
            i
        );
    }
}

/// Test particles (index >= N_active) get the "outside" branch of the
/// Jacobi recurrence, which must still invert exactly.
#[test]
fn jacobi_transform_round_trips_with_test_particles() {
    let mut inertial = five_body();
    let n = inertial.len();
    let n_active = 3;
    for p in inertial.iter_mut().skip(n_active) {
        p.m = 0.;
    }
    let mut p_j = scrambled(&inertial);
    reb_transformations_inertial_to_jacobi_posvel(&inertial, &mut p_j, &inertial, n, n_active);
    let mut back = scrambled(&inertial);
    reb_transformations_jacobi_to_inertial_posvel(&mut back, &p_j, &inertial, n, n_active);
    for i in 0..n {
        assert!(
            dist(&back[i], &inertial[i]) < 1e-14,
            "Jacobi round trip with test particles moved particle {}: |dr| = {}",
            i,
            dist(&back[i], &inertial[i])
        );
        assert!(
            vdist(&back[i], &inertial[i]) < 1e-14,
            "Jacobi round trip with test particles changed v of particle {}: |dv| = {}",
            i,
            vdist(&back[i], &inertial[i])
        );
    }
    // Test particles are referred to the total center of mass of the
    // active bodies, so their Jacobi offset is r_i - com(active).
    let inner = com_of(&inertial[0..n_active]);
    for i in n_active..n {
        let want = relative(&inertial[i], &inner);
        assert!(
            dist(&p_j[i], &want) < 1e-14,
            "test-particle Jacobi coordinate {} is not measured from the active center of mass: |dr| = {}",
            i,
            dist(&p_j[i], &want)
        );
    }
}

/// Democratic-heliocentric positions are literal heliocentric
/// differences, so the transform must reproduce them bit for bit.
#[test]
fn democraticheliocentric_positions_are_exactly_heliocentric() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_h = scrambled(&inertial);
    reb_transformations_inertial_to_democraticheliocentric_posvel(&inertial, &mut p_h, n, n);
    for i in 1..n {
        assert!(
            p_h[i].x.to_bits() == (inertial[i].x - inertial[0].x).to_bits()
                && p_h[i].y.to_bits() == (inertial[i].y - inertial[0].y).to_bits()
                && p_h[i].z.to_bits() == (inertial[i].z - inertial[0].z).to_bits(),
            "democratic-heliocentric position {} is not exactly r_{} - r_0",
            i,
            i
        );
    }
    let com = com_of(&inertial);
    assert!(
        dist(&p_h[0], &com) < 1e-15 && vdist(&p_h[0], &com) < 1e-15,
        "democratic-heliocentric coordinate 0 is not the barycenter: |dr| = {}, |dv| = {}",
        dist(&p_h[0], &com),
        vdist(&p_h[0], &com)
    );
    // The heliocentric momenta measured against the barycentric
    // velocity must sum, with the star, to zero total momentum offset.
    // The democratic-heliocentric momenta are barycentric momenta, so
    // the planets' must balance the star's: summing m_i (v_i - v_bary)
    // over every body has to give zero.
    let scale: f64 = inertial.iter().map(|p| p.m * vnorm(p)).sum();
    let mut px = inertial[0].m * (inertial[0].vx - p_h[0].vx);
    let mut py = inertial[0].m * (inertial[0].vy - p_h[0].vy);
    let mut pz = inertial[0].m * (inertial[0].vz - p_h[0].vz);
    for i in 1..n {
        px += inertial[i].m * p_h[i].vx;
        py += inertial[i].m * p_h[i].vy;
        pz += inertial[i].m * p_h[i].vz;
    }
    let residual = (px * px + py * py + pz * pz).sqrt();
    assert!(
        residual < 1e-15 * scale,
        "democratic-heliocentric momenta do not sum to zero in the barycentric frame: |p| = {} (scale {})",
        residual,
        scale
    );
}

#[test]
fn democraticheliocentric_transform_round_trips() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_h = scrambled(&inertial);
    reb_transformations_inertial_to_democraticheliocentric_posvel(&inertial, &mut p_h, n, n);
    let mut back = scrambled(&inertial);
    // The inverse transform reads the masses out of the target array.
    for i in 0..n {
        back[i].m = inertial[i].m;
    }
    reb_transformations_democraticheliocentric_to_inertial_posvel(&mut back, &p_h, n, n);
    for i in 0..n {
        assert!(
            dist(&back[i], &inertial[i]) < 1e-14,
            "democratic-heliocentric round trip moved particle {}: |dr| = {}",
            i,
            dist(&back[i], &inertial[i])
        );
        assert!(
            vdist(&back[i], &inertial[i]) < 1e-14,
            "democratic-heliocentric round trip changed v of particle {}: |dv| = {}",
            i,
            vdist(&back[i], &inertial[i])
        );
    }
}

#[test]
fn whds_transform_round_trips() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_h = scrambled(&inertial);
    reb_transformations_inertial_to_whds_posvel(&inertial, &mut p_h, n, n);
    let mut back = scrambled(&inertial);
    for i in 0..n {
        back[i].m = inertial[i].m;
    }
    reb_transformations_whds_to_inertial_posvel(&mut back, &p_h, n, n);
    for i in 0..n {
        assert!(
            dist(&back[i], &inertial[i]) < 1e-14,
            "WHDS round trip moved particle {}: |dr| = {}",
            i,
            dist(&back[i], &inertial[i])
        );
        assert!(
            vdist(&back[i], &inertial[i]) < 1e-14,
            "WHDS round trip changed v of particle {}: |dv| = {}",
            i,
            vdist(&back[i], &inertial[i])
        );
    }
    // WHDS shares its positions with the democratic-heliocentric set.
    let mut p_dh = scrambled(&inertial);
    reb_transformations_inertial_to_democraticheliocentric_posvel(&inertial, &mut p_dh, n, n);
    for i in 1..n {
        assert!(
            p_h[i].x.to_bits() == p_dh[i].x.to_bits()
                && p_h[i].y.to_bits() == p_dh[i].y.to_bits()
                && p_h[i].z.to_bits() == p_dh[i].z.to_bits(),
            "WHDS position {} differs from the democratic-heliocentric one",
            i
        );
    }
}

#[test]
fn barycentric_transform_round_trips() {
    let inertial = five_body();
    let n = inertial.len();
    let mut p_b = scrambled(&inertial);
    reb_transformations_inertial_to_barycentric_posvel(&inertial, &mut p_b, n, n);
    let com = com_of(&inertial);
    assert!(
        dist(&p_b[0], &com) < 1e-15 && vdist(&p_b[0], &com) < 1e-15,
        "barycentric coordinate 0 is not the barycenter: |dr| = {}, |dv| = {}",
        dist(&p_b[0], &com),
        vdist(&p_b[0], &com)
    );
    for i in 1..n {
        let want = relative(&inertial[i], &com);
        assert!(
            dist(&p_b[i], &want) < 1e-14,
            "barycentric coordinate {} is not r_{} - com: |dr| = {}",
            i,
            i,
            dist(&p_b[i], &want)
        );
    }
    let mut back = scrambled(&inertial);
    for i in 0..n {
        back[i].m = inertial[i].m;
    }
    reb_transformations_barycentric_to_inertial_posvel(&mut back, &p_b, n, n);
    for i in 0..n {
        assert!(
            dist(&back[i], &inertial[i]) < 1e-13,
            "barycentric round trip moved particle {}: |dr| = {}",
            i,
            dist(&back[i], &inertial[i])
        );
        assert!(
            vdist(&back[i], &inertial[i]) < 1e-13,
            "barycentric round trip changed v of particle {}: |dv| = {}",
            i,
            vdist(&back[i], &inertial[i])
        );
    }
}

// =====================================================================
// WHFast
// =====================================================================

/// For N = 2 in Jacobi coordinates WHFast's interaction step vanishes
/// identically (the 0-1 gravity term is switched off and there is no
/// third body), so the integrator reduces to the exact Kepler flow of
/// the relative coordinate. dt is a power of two so that r.t is
/// accumulated without rounding.
#[test]
fn whfast_jacobi_two_body_reproduces_the_analytic_orbit() {
    let (a, e, m1) = (1.0, 0.3, 1e-3);
    let dt = 1.0 / 64.0;
    let n_steps = 2000usize;
    let mut r = two_body_sim("whfast", m1, a, e, dt);
    let star = reb_particle { m: 1., ..P0 };
    let mu = r.G * (1. + m1);

    reb_simulation_steps(&mut r, n_steps);
    let t = r.t;
    assert!(
        t == dt * (n_steps as f64),
        "expected t = {}, got {}",
        dt * (n_steps as f64),
        t
    );
    let got = relative(&r.particles[1], &r.particles[0]);
    let want = reference_advance(1.0, star, m1, a, e, 0., 0., 0., 0., t);
    assert!(
        dist(&got, &want) < 1e-10 * a,
        "WHFast/Jacobi two-body separation drifted from the analytic orbit after {} steps: |dr| = {}",
        n_steps,
        dist(&got, &want)
    );
    assert!(
        vdist(&got, &want) < 1e-10 * (mu / a).sqrt(),
        "WHFast/Jacobi two-body velocity drifted from the analytic orbit: |dv| = {}",
        vdist(&got, &want)
    );
}

/// WHDS also makes the isolated two-body problem exact: its jump step
/// cancels identically for N = 2 and its Kepler step already uses
/// mu = G(m0 + m1).
#[test]
fn whfast_whds_two_body_reproduces_the_analytic_orbit() {
    let (a, e, m1) = (1.0, 0.45, 1e-3);
    let dt = 1.0 / 64.0;
    let n_steps = 1000usize;
    let mut r = two_body_sim("whfast", m1, a, e, dt);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.coordinates = integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_WHDS;
    }
    let star = reb_particle { m: 1., ..P0 };
    reb_simulation_steps(&mut r, n_steps);
    let got = relative(&r.particles[1], &r.particles[0]);
    let want = reference_advance(1.0, star, m1, a, e, 0., 0., 0., 0., r.t);
    assert!(
        dist(&got, &want) < 1e-10 * a,
        "WHFast/WHDS two-body separation drifted from the analytic orbit: |dr| = {}",
        dist(&got, &want)
    );
}

#[test]
fn whfast_conserves_energy_and_angular_momentum_on_a_three_body_system() {
    let mut r = three_body_sim("whfast", 2. * PI / 200.);
    let l0 = reb_simulation_angular_momentum(&r);
    let drift = energy_drift(&mut r, 4000);
    assert!(
        drift < 1e-7,
        "WHFast relative energy drift over 4000 steps is {}, expected < 1e-7",
        drift
    );
    // Angular momentum is not merely bounded but conserved to roundoff:
    // every Kepler step and every kick preserves it exactly.
    let l1 = reb_simulation_angular_momentum(&r);
    assert!(
        vec3d_dist(&l0, &l1) < 1e-12 * vec3d_norm(&l0),
        "WHFast did not conserve angular momentum: |dL|/|L| = {}",
        vec3d_dist(&l0, &l1) / vec3d_norm(&l0)
    );
    assert!(
        !has_error(&r),
        "WHFast reported an error: {:?}",
        r.messages
    );
}

/// WHFast without correctors is a second-order scheme, so halving the
/// timestep must divide the energy error by four. Anything else means
/// the operator splitting has lost its symmetry.
#[test]
fn whfast_energy_error_scales_as_the_square_of_the_timestep() {
    let mut drifts = Vec::new();
    for &n in &[100usize, 200, 400] {
        // Same physical end time in every run: 4 orbits of the inner
        // planet.
        let mut r = three_body_sim("whfast", 2. * PI / (n as f64));
        drifts.push(energy_drift(&mut r, 4 * n));
    }
    for k in 0..2 {
        let ratio = drifts[k] / drifts[k + 1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "WHFast energy error ratio between {} and {} steps per orbit is {}, expected ~4 for a second-order scheme (drifts {:e} and {:e})",
            100 << k,
            200 << k,
            ratio,
            drifts[k],
            drifts[k + 1]
        );
    }
}

/// Every legal first-corrector order must run and conserve energy, and
/// a corrector has to earn its cost: switching one on must shrink the
/// energy error by orders of magnitude relative to the uncorrected run.
#[test]
fn whfast_all_corrector_orders_conserve_energy() {
    let mut drifts = Vec::new();
    for &order in &[0u32, 3, 5, 7, 11, 17] {
        let mut r = three_body_sim("whfast", 2. * PI / 200.);
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            wh.corrector = order;
        }
        let drift = energy_drift(&mut r, 1000);
        assert!(
            !has_error(&r),
            "WHFast corrector order {} reported an error: {:?}",
            order,
            r.messages
        );
        assert!(
            drift < 1e-7,
            "WHFast corrector order {}: relative energy drift {} exceeds 1e-7",
            order,
            drift
        );
        drifts.push((order, drift));
    }
    let uncorrected = drifts[0].1;
    for &(order, drift) in drifts.iter().skip(1) {
        assert!(
            drift < uncorrected / 100.,
            "corrector order {} gives energy drift {:e}, not at least 100x better than the uncorrected {:e}",
            order,
            drift,
            uncorrected
        );
    }
}

/// The second symplectic corrector is a separate code path (operators
/// C, Y and U) and must also run and conserve energy.
#[test]
fn whfast_second_corrector_conserves_energy() {
    let mut r = three_body_sim("whfast", 2. * PI / 200.);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.corrector = 11;
        wh.corrector2 = 1;
    }
    let drift = energy_drift(&mut r, 500);
    assert!(
        !has_error(&r),
        "WHFast with corrector2 reported an error: {:?}",
        r.messages
    );
    assert!(
        drift < 1e-9,
        "WHFast with the second corrector: relative energy drift {} exceeds 1e-9",
        drift
    );
    // The second corrector is applied on top of an 11th-order first
    // corrector, so the result must be far better than no corrector.
    let mut plain = three_body_sim("whfast", 2. * PI / 200.);
    let plain_drift = energy_drift(&mut plain, 500);
    assert!(
        drift < plain_drift / 100.,
        "corrector2 gives energy drift {:e}, not at least 100x better than the uncorrected {:e}",
        drift,
        plain_drift
    );
}

/// A corrector is a change of variables, so it has to move the state:
/// if the corrector code were a no-op this test would fail. The
/// corrected run must still track the uncorrected one closely, because
/// the two differ only by a near-identity transformation.
#[test]
fn whfast_correctors_change_the_trajectory_without_changing_the_orbit() {
    let dt = 2. * PI / 60.;
    let n = 300usize;
    let mut plain = three_body_sim("whfast", dt);
    reb_simulation_steps(&mut plain, n);

    let mut corrected = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = corrected.integrator {
        wh.corrector = 17;
    }
    reb_simulation_steps(&mut corrected, n);

    let d = dist(&plain.particles[1], &corrected.particles[1]);
    assert!(
        d > 1e-13,
        "corrector 17 produced a bit-identical trajectory (|dr| = {}), so the corrector code did not run",
        d
    );
    assert!(
        d < 1e-2,
        "corrector 17 moved the planet by {}, far more than a near-identity change of variables should",
        d
    );
}

/// The `init` guard rejects unsupported corrector orders, kernels and
/// coordinate/kernel combinations.
#[test]
fn whfast_rejects_unsupported_configurations() {
    // Corrector order 4 is not one of 0, 3, 5, 7, 11, 17.
    let mut r = three_body_sim("whfast", 0.01);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.corrector = 4;
    }
    reb_simulation_steps(&mut r, 1);
    assert!(
        has_message(&r, "First symplectic correctors are only available"),
        "corrector order 4 was accepted; messages: {:?}",
        r.messages
    );

    // Kernel 4 does not exist.
    let mut r = three_body_sim("whfast", 0.01);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.kernel = 4;
    }
    reb_simulation_steps(&mut r, 1);
    assert!(
        has_message(&r, "Kernel method must be"),
        "kernel 4 was accepted; messages: {:?}",
        r.messages
    );

    // Non-default kernels need Jacobi coordinates.
    let mut r = three_body_sim("whfast", 0.01);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.kernel = integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION;
        wh.coordinates = integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    }
    reb_simulation_steps(&mut r, 1);
    assert!(
        has_message(&r, "Non-standard kernel requires Jacobi coordinates"),
        "a composition kernel in democratic-heliocentric coordinates was accepted; messages: {:?}",
        r.messages
    );

    // Correctors need Jacobi or barycentric coordinates.
    let mut r = three_body_sim("whfast", 0.01);
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.corrector = 5;
        wh.coordinates = integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_WHDS;
    }
    reb_simulation_steps(&mut r, 1);
    assert!(
        has_message(&r, "Symplectic correctors are only compatible"),
        "a corrector in WHDS coordinates was accepted; messages: {:?}",
        r.messages
    );
}

/// All four kernels must run and conserve energy. Kernels 1 and 3
/// switch the gravity routine to the Jacobi one, kernel 2 runs the
/// composition scheme with its extra Kepler sub-steps.
#[test]
fn whfast_all_kernels_conserve_energy() {
    for &kernel in &[
        integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT,
        integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK,
        integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION,
        integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_LAZY,
    ] {
        let mut r = three_body_sim("whfast", 2. * PI / 200.);
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            wh.kernel = kernel;
        }
        let drift = energy_drift(&mut r, 800);
        assert!(
            !has_error(&r),
            "WHFast kernel {} reported an error: {:?}",
            kernel,
            r.messages
        );
        assert!(
            drift < 1e-7,
            "WHFast kernel {}: relative energy drift {} exceeds 1e-7",
            kernel,
            drift
        );
    }
}

/// All four coordinate systems must run and conserve energy.
#[test]
fn whfast_all_coordinate_systems_conserve_energy() {
    for &coord in &[
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_WHDS,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC,
    ] {
        let mut r = three_body_sim("whfast", 2. * PI / 400.);
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            wh.coordinates = coord;
        }
        let drift = energy_drift(&mut r, 800);
        assert!(
            !has_error(&r),
            "WHFast coordinates {} reported an error: {:?}",
            coord,
            r.messages
        );
        assert!(
            drift < 1e-6,
            "WHFast coordinates {}: relative energy drift {} exceeds 1e-6",
            coord,
            drift
        );
    }
}

/// The four coordinate systems are four splittings of the same
/// Hamiltonian, so with a small timestep they have to agree with a
/// high-accuracy IAS15 run of the same initial condition.
#[test]
fn whfast_coordinate_systems_agree_with_ias15() {
    let dt = 2. * PI / 800.;
    let n = 1600usize;
    let mut reference = three_body_sim("ias15", dt);
    let tmax = dt * (n as f64);
    reb_simulation_integrate(&mut reference, tmax);

    for &coord in &[
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_WHDS,
        integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC,
    ] {
        let mut r = three_body_sim("whfast", dt);
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            wh.coordinates = coord;
        }
        reb_simulation_integrate(&mut r, tmax);
        for i in 0..r.N {
            assert!(
                dist(&r.particles[i], &reference.particles[i]) < 1e-6,
                "WHFast coordinates {} disagree with IAS15 on particle {} after t = {}: |dr| = {}",
                coord,
                i,
                tmax,
                dist(&r.particles[i], &reference.particles[i])
            );
        }
    }
}

/// One WHFast step is a palindromic operator sequence, so reversing the
/// sign of dt reverses the map.
#[test]
fn whfast_is_time_reversible() {
    let dt = 2. * PI / 100.;
    // 2.5 orbits of the inner planet, so it ends up on the far side of
    // the star rather than back where it started.
    let n = 250usize;
    let mut r = three_body_sim("whfast", dt);
    let start: Vec<reb_particle> = r.particles.clone();
    reb_simulation_steps(&mut r, n);
    let midway = dist(&r.particles[1], &start[1]);
    assert!(
        midway > 1.0,
        "the forward integration barely moved (|dr| = {}), the reversibility test would be vacuous",
        midway
    );
    r.dt = -dt;
    reb_simulation_steps(&mut r, n);
    for i in 0..r.N {
        assert!(
            dist(&r.particles[i], &start[i]) < 1e-11,
            "WHFast is not time reversible: particle {} is off by {} after {} steps out and back",
            i,
            dist(&r.particles[i], &start[i]),
            n
        );
        assert!(
            vdist(&r.particles[i], &start[i]) < 1e-11,
            "WHFast is not time reversible in velocity: particle {} is off by {}",
            i,
            vdist(&r.particles[i], &start[i])
        );
    }
    assert!(
        r.t.abs() < 1e-12,
        "the reversed integration did not return to t = 0, got {}",
        r.t
    );
}

#[test]
fn whfast_is_bitwise_deterministic() {
    let dt = 2. * PI / 137.0;
    let mut a = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = a.integrator {
        wh.corrector = 11;
        wh.kernel = integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_LAZY;
    }
    reb_simulation_steps(&mut a, 500);

    let mut b = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = b.integrator {
        wh.corrector = 11;
        wh.kernel = integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_LAZY;
    }
    reb_simulation_steps(&mut b, 500);

    assert!(
        states_same_bits(&a.particles, &b.particles),
        "two identical WHFast runs did not produce bit-identical particles"
    );
    assert!(
        a.t.to_bits() == b.t.to_bits(),
        "two identical WHFast runs disagree on the final time: {} vs {}",
        a.t,
        b.t
    );
    assert!(
        reb_simulation_energy(&a).to_bits() == reb_simulation_energy(&b).to_bits(),
        "two identical WHFast runs disagree on the final energy"
    );
}

/// safe_mode = 0 merges the trailing and leading half drift steps of
/// consecutive timesteps. Once the run is synchronized the result is
/// the same map, so it may differ only by rounding.
#[test]
fn whfast_safe_mode_zero_agrees_with_safe_mode_one() {
    let dt = 2. * PI / 200.;
    let n = 1000usize;
    let mut safe = three_body_sim("whfast", dt);
    reb_simulation_steps(&mut safe, n);

    let mut fast = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = fast.integrator {
        wh.safe_mode = 0;
    }
    reb_simulation_steps(&mut fast, n);

    assert!(
        fast.is_synchronized == 1,
        "reb_simulation_steps should leave the simulation synchronized"
    );
    for i in 0..safe.N {
        assert!(
            dist(&safe.particles[i], &fast.particles[i]) < 1e-11,
            "safe_mode 0 and 1 disagree on particle {}: |dr| = {}",
            i,
            dist(&safe.particles[i], &fast.particles[i])
        );
    }
}

/// With keep_unsynchronized = 1 a synchronization is a pure read-out:
/// the internal Jacobi state is restored afterwards, so continuing the
/// integration must be bit-identical to never having synchronized.
#[test]
fn whfast_keep_unsynchronized_continues_identically() {
    let dt = 2. * PI / 150.;
    let mut split = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = split.integrator {
        wh.safe_mode = 0;
        wh.keep_unsynchronized = 1;
    }
    reb_simulation_steps(&mut split, 250);
    let readout: Vec<reb_particle> = split.particles.clone();
    assert!(
        split.is_synchronized == 0,
        "keep_unsynchronized should leave is_synchronized at 0, got {}",
        split.is_synchronized
    );
    reb_simulation_steps(&mut split, 250);

    let mut straight = three_body_sim("whfast", dt);
    if let reb_integrator_state::whfast(ref mut wh) = straight.integrator {
        wh.safe_mode = 0;
        wh.keep_unsynchronized = 1;
    }
    reb_simulation_steps(&mut straight, 500);

    assert!(
        states_same_bits(&split.particles, &straight.particles),
        "keep_unsynchronized changed the trajectory: 250+250 steps differ from 500 steps"
    );
    // The intermediate read-out is a real synchronized state, not the
    // internal half-drifted one.
    assert!(
        dist(&readout[1], &split.particles[1]) > 1e-6,
        "the simulation did not advance between the two halves of the run"
    );
}

/// A single particle has no interactions: WHFast must move it on a
/// straight line and leave its velocity untouched, bit for bit.
#[test]
fn whfast_single_particle_drifts_linearly() {
    let mut r = reb_simulation_create();
    r.G = 1.;
    r.save_messages = 1;
    reb_simulation_set_integrator(&mut r, "whfast");
    // A power-of-two mass keeps the inertial <-> Jacobi round trip of a
    // single particle exact (it is a multiply by m followed by a divide
    // by m), so the velocity can be checked bit for bit.
    let p = particle(2.0, 0.25, -0.5, 0.75, 0.3, -0.2, 0.1);
    reb_simulation_add(&mut r, p);
    r.dt = 1.0 / 32.0;
    let n = 512usize;
    reb_simulation_steps(&mut r, n);
    let t = r.dt * (n as f64);
    assert!(
        r.t == t,
        "single-particle run: expected t = {}, got {}",
        t,
        r.t
    );
    let q = r.particles[0];
    assert!(
        q.vx.to_bits() == p.vx.to_bits()
            && q.vy.to_bits() == p.vy.to_bits()
            && q.vz.to_bits() == p.vz.to_bits(),
        "a lone particle's velocity changed: ({}, {}, {}) -> ({}, {}, {})",
        p.vx,
        p.vy,
        p.vz,
        q.vx,
        q.vy,
        q.vz
    );
    let want = particle(p.m, p.x + t * p.vx, p.y + t * p.vy, p.z + t * p.vz, 0., 0., 0.);
    assert!(
        dist(&q, &want) < 1e-11,
        "a lone particle did not drift linearly: |dr| = {}",
        dist(&q, &want)
    );
    assert!(
        !has_error(&r),
        "the single-particle run reported an error: {:?}",
        r.messages
    );
}

/// An empty simulation exits with REB_STATUS_NO_PARTICLES instead of
/// stepping.
#[test]
fn whfast_with_no_particles_exits_cleanly() {
    let mut r = reb_simulation_create();
    r.G = 1.;
    r.save_messages = 1;
    reb_simulation_set_integrator(&mut r, "whfast");
    r.dt = 0.01;
    let status = reb_simulation_integrate(&mut r, 1.0);
    assert!(
        status == REB_STATUS_NO_PARTICLES,
        "an empty simulation returned status {}, expected {}",
        status,
        REB_STATUS_NO_PARTICLES
    );
    assert!(r.t == 0., "an empty simulation advanced time to {}", r.t);
    assert!(
        has_message(&r, "No particles in simulation"),
        "no warning about the empty simulation: {:?}",
        r.messages
    );
}

/// WHFast tracks a high-eccentricity orbit as long as the timestep
/// resolves periapsis; energy must stay bounded and the semi-major axis
/// must not wander.
#[test]
fn whfast_conserves_energy_for_an_eccentric_two_body_orbit() {
    let (a, e, m1) = (1.0, 0.9, 1e-3);
    let mut r = two_body_sim("whfast", m1, a, e, 2. * PI / 4000.);
    let star = reb_particle { m: 1., ..P0 };
    let o0 = reb_orbit_from_particle(r.G, relative(&r.particles[1], &r.particles[0]), star);
    let drift = energy_drift(&mut r, 8000);
    // The isolated pair is integrated by the exact Kepler flow, so only
    // rounding accumulates over the 8000 steps.
    assert!(
        drift < 1e-12,
        "WHFast on an e = 0.9 two-body orbit drifted in energy by {}, expected < 1e-12",
        drift
    );
    let o1 = reb_orbit_from_particle(r.G, relative(&r.particles[1], &r.particles[0]), star);
    assert!(
        ((o1.a - o0.a) / o0.a).abs() < 1e-11,
        "semi-major axis drifted from {} to {}",
        o0.a,
        o1.a
    );
    assert!(
        (o1.e - o0.e).abs() < 1e-11,
        "eccentricity drifted from {} to {}",
        o0.e,
        o1.e
    );
}

// =====================================================================
// SABA
// =====================================================================

const SABA_TYPES: [i32; 18] = [
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_1,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_2,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_3,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CM_1,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CM_2,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CM_3,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CM_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_1,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_2,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_3,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_10_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_8_6_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_10_6_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_H_8_4_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_H_8_6_4,
    integrator_saba::REB_INTEGRATOR_SABA_TYPE_H_10_6_4,
];

#[test]
fn saba_all_types_conserve_energy() {
    for &t in SABA_TYPES.iter() {
        let mut r = three_body_sim("saba", 2. * PI / 200.);
        if let reb_integrator_state::saba(ref mut s) = r.integrator {
            s.type_ = t;
        }
        let drift = energy_drift(&mut r, 800);
        assert!(
            !has_error(&r),
            "SABA type 0x{:x} reported an error: {:?}",
            t,
            r.messages
        );
        assert!(
            drift < 1e-7,
            "SABA type 0x{:x}: relative energy drift {} exceeds 1e-7",
            t,
            drift
        );
    }
}

/// Every SABA variant except the three one-stage schemes (SABA1 and its
/// two corrector flavours, which are all second order) must beat plain
/// SABA1 by at least four orders of magnitude at the same timestep.
#[test]
fn saba_multistage_variants_beat_saba1() {
    let dt = 2. * PI / 200.;
    let mut base = three_body_sim("saba", dt);
    if let reb_integrator_state::saba(ref mut s) = base.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_1;
    }
    let base_drift = energy_drift(&mut base, 800);
    let second_order = [
        integrator_saba::REB_INTEGRATOR_SABA_TYPE_1,
        integrator_saba::REB_INTEGRATOR_SABA_TYPE_CM_1,
        integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_1,
    ];
    for &t in SABA_TYPES.iter() {
        if second_order.contains(&t) {
            continue;
        }
        let mut r = three_body_sim("saba", dt);
        if let reb_integrator_state::saba(ref mut s) = r.integrator {
            s.type_ = t;
        }
        let drift = energy_drift(&mut r, 800);
        assert!(
            drift < base_drift / 1e3,
            "SABA type 0x{:x} gives energy drift {:e}, not at least 1000 times better than SABA1's {:e}",
            t,
            drift,
            base_drift
        );
    }
}

/// Like WHFast, SABA without correctors reduces to the exact Kepler
/// flow for an isolated pair.
#[test]
fn saba_two_body_reproduces_the_analytic_orbit() {
    let (a, e, m1) = (1.0, 0.25, 1e-3);
    let dt = 1.0 / 64.0;
    let n_steps = 1000usize;
    let mut r = two_body_sim("saba", m1, a, e, dt);
    if let reb_integrator_state::saba(ref mut s) = r.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_4;
    }
    let star = reb_particle { m: 1., ..P0 };
    reb_simulation_steps(&mut r, n_steps);
    let got = relative(&r.particles[1], &r.particles[0]);
    let want = reference_advance(1.0, star, m1, a, e, 0., 0., 0., 0., r.t);
    assert!(
        dist(&got, &want) < 1e-10 * a,
        "SABA4 two-body separation drifted from the analytic orbit: |dr| = {}",
        dist(&got, &want)
    );
}

/// SABA4 is a fourth-order scheme, SABA1 is second order, so at the
/// same timestep SABA4 must conserve energy far better.
#[test]
fn saba_high_order_beats_low_order() {
    let dt = 2. * PI / 100.;
    let mut low = three_body_sim("saba", dt);
    if let reb_integrator_state::saba(ref mut s) = low.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_1;
    }
    let d_low = energy_drift(&mut low, 400);

    let mut high = three_body_sim("saba", dt);
    if let reb_integrator_state::saba(ref mut s) = high.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_4;
    }
    let d_high = energy_drift(&mut high, 400);

    assert!(
        d_high < 0.1 * d_low,
        "SABA4 energy drift ({}) is not clearly smaller than SABA1's ({})",
        d_high,
        d_low
    );
}

#[test]
fn saba_rejects_an_invalid_type() {
    let mut r = three_body_sim("saba", 0.01);
    if let reb_integrator_state::saba(ref mut s) = r.integrator {
        s.type_ = 0x42;
    }
    let before: Vec<reb_particle> = r.particles.clone();
    reb_simulation_steps(&mut r, 1);
    assert!(
        has_message(&r, "Invalid SABA integrator type"),
        "SABA accepted type 0x42; messages: {:?}",
        r.messages
    );
    assert!(
        states_same_bits(&before, &r.particles),
        "SABA moved the particles despite rejecting the integrator type"
    );
}

#[test]
fn saba_is_bitwise_deterministic() {
    let dt = 2. * PI / 111.0;
    let mut a = three_body_sim("saba", dt);
    if let reb_integrator_state::saba(ref mut s) = a.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_4;
    }
    reb_simulation_steps(&mut a, 400);

    let mut b = three_body_sim("saba", dt);
    if let reb_integrator_state::saba(ref mut s) = b.integrator {
        s.type_ = integrator_saba::REB_INTEGRATOR_SABA_TYPE_CL_4;
    }
    reb_simulation_steps(&mut b, 400);

    assert!(
        states_same_bits(&a.particles, &b.particles),
        "two identical SABA runs did not produce bit-identical particles"
    );
}

// =====================================================================
// EOS (embedded operator splitting)
// =====================================================================

const EOS_TYPES: [i32; 9] = [
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF4,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF6,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF8,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF4_2,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_PLF7_6_4,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_PMLF4,
    integrator_eos::REB_INTEGRATOR_EOS_TYPE_PMLF6,
];

#[test]
fn eos_all_outer_splittings_conserve_energy() {
    for &phi0 in EOS_TYPES.iter() {
        let mut r = three_body_sim("eos", 2. * PI / 2000.);
        if let reb_integrator_state::eos(ref mut e) = r.integrator {
            e.phi0 = phi0;
        }
        let drift = energy_drift(&mut r, 2000);
        assert!(
            !has_error(&r),
            "EOS phi0 = {} reported an error: {:?}",
            phi0,
            r.messages
        );
        assert!(
            drift < 1e-7,
            "EOS phi0 = {}: relative energy drift {} exceeds 1e-7",
            phi0,
            drift
        );
    }
}

#[test]
fn eos_inner_splittings_conserve_energy() {
    for &phi1 in EOS_TYPES.iter() {
        let mut r = three_body_sim("eos", 2. * PI / 2000.);
        if let reb_integrator_state::eos(ref mut e) = r.integrator {
            e.phi1 = phi1;
            e.n = 2;
        }
        let drift = energy_drift(&mut r, 1000);
        assert!(
            !has_error(&r),
            "EOS phi1 = {} reported an error: {:?}",
            phi1,
            r.messages
        );
        assert!(
            drift < 1e-7,
            "EOS phi1 = {}: relative energy drift {} exceeds 1e-7",
            phi1,
            drift
        );
    }
}

#[test]
fn eos_is_bitwise_deterministic() {
    let dt = 2. * PI / 2000.;
    let mut a = three_body_sim("eos", dt);
    if let reb_integrator_state::eos(ref mut e) = a.integrator {
        e.phi0 = integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4;
        e.phi1 = integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF4;
        e.n = 3;
    }
    reb_simulation_steps(&mut a, 500);

    let mut b = three_body_sim("eos", dt);
    if let reb_integrator_state::eos(ref mut e) = b.integrator {
        e.phi0 = integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF8_6_4;
        e.phi1 = integrator_eos::REB_INTEGRATOR_EOS_TYPE_LF4;
        e.n = 3;
    }
    reb_simulation_steps(&mut b, 500);

    assert!(
        states_same_bits(&a.particles, &b.particles),
        "two identical EOS runs did not produce bit-identical particles"
    );
}

// =====================================================================
// leapfrog
// =====================================================================

#[test]
fn leapfrog_conserves_energy_and_momentum() {
    let mut r = three_body_sim("leapfrog", 2. * PI / 4000.);
    let l0 = reb_simulation_angular_momentum(&r);
    let drift = energy_drift(&mut r, 4000);
    assert!(
        drift < 1e-7,
        "leapfrog relative energy drift {} exceeds 1e-7",
        drift
    );
    let l1 = reb_simulation_angular_momentum(&r);
    assert!(
        vec3d_dist(&l0, &l1) < 1e-12 * vec3d_norm(&l0),
        "leapfrog did not conserve angular momentum: |dL| = {}",
        vec3d_dist(&l0, &l1)
    );
}

// =====================================================================
// cross-integrator agreement
// =====================================================================

/// The same two-body problem integrated by IAS15 and by WHFast has to
/// end up in the same place. WHFast/Jacobi is exact here, so the
/// residual is IAS15's own error budget.
#[test]
fn whfast_and_ias15_agree_on_a_two_body_problem() {
    let (a, e, m1) = (1.0, 0.2, 1e-3);
    let tmax = 20. * 2. * PI;
    let mut wh = two_body_sim("whfast", m1, a, e, 2. * PI / 200.);
    reb_simulation_integrate(&mut wh, tmax);
    let mut ias = two_body_sim("ias15", m1, a, e, 2. * PI / 200.);
    reb_simulation_integrate(&mut ias, tmax);
    for i in 0..2 {
        assert!(
            dist(&wh.particles[i], &ias.particles[i]) < 1e-8,
            "WHFast and IAS15 disagree on particle {} after 20 orbits: |dr| = {}",
            i,
            dist(&wh.particles[i], &ias.particles[i])
        );
    }
}

/// WHFast and SABA are different splittings of the same Hamiltonian and
/// must converge to the same trajectory as the timestep shrinks; both
/// are compared against IAS15.
#[test]
fn whfast_and_saba_agree_with_ias15_on_a_three_body_problem() {
    let dt = 2. * PI / 1000.;
    let n = 2000usize;
    let tmax = dt * (n as f64);

    let mut ias = three_body_sim("ias15", dt);
    reb_simulation_integrate(&mut ias, tmax);

    let mut wh = three_body_sim("whfast", dt);
    reb_simulation_integrate(&mut wh, tmax);

    let mut saba = three_body_sim("saba", dt);
    reb_simulation_integrate(&mut saba, tmax);

    for i in 0..ias.N {
        let dwh = dist(&wh.particles[i], &ias.particles[i]);
        let dsaba = dist(&saba.particles[i], &ias.particles[i]);
        assert!(
            dwh < 1e-7,
            "WHFast disagrees with IAS15 on particle {}: |dr| = {}",
            i,
            dwh
        );
        assert!(
            dsaba < 1e-7,
            "SABA disagrees with IAS15 on particle {}: |dr| = {}",
            i,
            dsaba
        );
        assert!(
            dist(&wh.particles[i], &saba.particles[i]) < 1e-7,
            "WHFast and SABA disagree on particle {}: |dr| = {}",
            i,
            dist(&wh.particles[i], &saba.particles[i])
        );
    }
}

