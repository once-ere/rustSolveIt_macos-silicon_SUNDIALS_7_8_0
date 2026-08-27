//! Integration tests for the interpolation_and_spin group of reboundx_rs.
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

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/// A bare simulation with REBOUNDx attached, messages captured instead of
/// printed (several code paths under test warn on every timestep).
fn attached_sim() -> reb_simulation {
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;
    rebx_attach(&mut sim);
    sim
}

/// `n` knots of `f` spread uniformly over `[lo, hi]` (`times`, `values`).
fn sample<F: Fn(f64) -> f64>(lo: f64, hi: f64, n: usize, f: F) -> (Vec<f64>, Vec<f64>) {
    let mut times = Vec::with_capacity(n);
    let mut values = Vec::with_capacity(n);
    for i in 0..n {
        let t = lo + (hi - lo) * (i as f64) / ((n - 1) as f64);
        times.push(t);
        values.push(f(t));
    }
    (times, values)
}

fn vlen(v: reb_vec3d) -> f64 {
    (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

/// Total angular momentum (orbital + spin) of a simulation with REBOUNDx.
fn total_L(sim: &mut reb_simulation) -> reb_vec3d {
    let L_orb = reb_simulation_angular_momentum(sim);
    let L_spin = rebx_with(sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx))
        .expect("extras attached");
    reb_vec3d_add(L_orb, L_spin)
}

/// Star + hot Jupiter with full tides_spin parameters on both bodies,
/// rotated into the invariable plane exactly as the REBOUNDx examples do.
/// `with_k2` selects whether the bodies are given a Love number at all: with
/// no `k2` anywhere the tides_spin force is a no-op by construction (see
/// `rebx_tides_spin`, which `continue`s unless `k2` is set).
fn hot_jupiter(with_k2: bool, evolve_spins: bool) -> reb_simulation {
    hot_jupiter_dt(with_k2, evolve_spins, 1e-3)
}

fn hot_jupiter_dt(with_k2: bool, evolve_spins: bool, dt: f64) -> reb_simulation {
    let mut sim = reb_simulation_create();
    sim.save_messages = 1;

    let solar_mass = 1.;
    let solar_rad = 0.00465;
    reb_simulation_add_fmt(
        &mut sim,
        "m r",
        &[reb_fmt_arg::d(solar_mass), reb_fmt_arg::d(solar_rad)],
    );

    let p1_mass = 9.55e-4;
    let p1_rad = 4.676e-4;
    reb_simulation_add_fmt(
        &mut sim,
        "m a e inc r",
        &[
            reb_fmt_arg::d(p1_mass),
            reb_fmt_arg::d(0.04072),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(p1_rad),
        ],
    );

    sim.N_active = 2;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = dt;

    rebx_attach(&mut sim);
    let effect = rebx_load_force(&mut sim, "tides_spin").expect("tides_spin force");
    rebx_add_force(&mut sim, effect);

    let solar_spin_period = 27. * 2. * PI / 365.;
    let solar_spin = (2. * PI) / solar_spin_period;
    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    let solar_tau = 1. / (2. * 1e6 * orb.n);

    let spin_period_1 = 0.5 * 2. * PI / 365.;
    let spin_1 = (2. * PI) / spin_period_1;
    let Omega_1 = reb_tools_spherical_to_xyz(spin_1, 30. * (PI / 180.), 0.);

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        if with_k2 {
            rebx_set_param_double(rebx, rebx_ap::particle(0), "k2", 0.07);
            rebx_set_param_double(rebx, rebx_ap::particle(1), "k2", 0.3);
            rebx_set_param_double(rebx, rebx_ap::particle(0), "tau", solar_tau);
            rebx_set_param_double(rebx, rebx_ap::particle(1), "tau", 1. / (2. * 10000. * orb.n));
        }
        if evolve_spins {
            rebx_set_param_vec3d(
                rebx,
                rebx_ap::particle(0),
                "Omega",
                reb_vec3d { x: 0., y: 0., z: solar_spin },
            );
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(0),
                "I",
                0.07 * solar_mass * solar_rad * solar_rad,
            );
            rebx_set_param_vec3d(rebx, rebx_ap::particle(1), "Omega", Omega_1);
            rebx_set_param_double(
                rebx,
                rebx_ap::particle(1),
                "I",
                0.25 * p1_mass * p1_rad * p1_rad,
            );
        }
    }

    reb_simulation_move_to_com(&mut sim);

    if evolve_spins {
        let L_orb = reb_simulation_angular_momentum(&sim);
        let L_spin = rebx_with(&mut sim, |sim, rebx| {
            rebx_tools_spin_angular_momentum(sim, rebx)
        })
        .expect("extras attached");
        let newz = reb_vec3d_add(L_orb, L_spin);
        let newx = reb_vec3d_cross(reb_vec3d { x: 0., y: 0., z: 1. }, newz);
        let rot = reb_rotation_init_to_new_axes(newz, newx);
        rebx_with(&mut sim, |sim, rebx| {
            rebx_simulation_irotate(sim, rebx, rot);
        });
        rebx_with(&mut sim, |sim, rebx| {
            rebx_spin_initialize_ode(sim, rebx, effect);
        });
    }

    sim
}

// ===========================================================================
// rebx_interpolator
// ===========================================================================

/// The cubic-spline evaluator (`rebx_splint`) reduces, at `x == xa[klo]`, to
/// `1*ya[klo] + 0*ya[klo+1] + (0*y2[klo] + 0*y2[klo+1])*h*h/6`, so a spline
/// must return the tabulated value at a knot with no error at all — not
/// "close to", but the identical double. Checked bit-for-bit, sweeping the
/// knots forwards (the `klo` cache walks forward one interval per call).
#[test]
fn interpolator_is_bit_exact_at_its_knots_sweeping_forwards() {
    let (times, values) = sample(0., 2. * PI, 65, |t| (t).sin() + 0.5 * (2.3 * t).cos());
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");
    let mut interp = rebx_create_interpolator(
        rebx,
        times.len() as i32,
        &times,
        &values,
        rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
    );

    for k in 0..times.len() {
        let got = rebx_interpolate(rebx, &mut interp, times[k]);
        assert_eq!(
            got.to_bits(),
            values[k].to_bits(),
            "spline at knot {} (t = {}): got {:?} (bits {:016x}), tabulated {:?} (bits {:016x})",
            k,
            times[k],
            got,
            got.to_bits(),
            values[k],
            values[k].to_bits()
        );
    }
    rebx_free_interpolator(interp);
}

/// Same exactness claim, sweeping the knots backwards. This drives the
/// `xa[klo] > x` branch of `rebx_splint`, whose two-stage decrement must land
/// on `klo == k` when `x == xa[k]`; if it landed one interval off, `a` would
/// be 0 and the returned value would be the *neighbouring* knot's.
#[test]
fn interpolator_is_bit_exact_at_its_knots_sweeping_backwards() {
    let (times, values) = sample(-3., 7., 41, |t| (0.4 * t).exp() * (0.7 * t).sin());
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");
    let mut interp = rebx_create_interpolator(
        rebx,
        times.len() as i32,
        &times,
        &values,
        rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
    );

    // Walk the cache up to the far end first, then come back down.
    let _ = rebx_interpolate(rebx, &mut interp, times[times.len() - 1]);
    for k in (0..times.len()).rev() {
        let got = rebx_interpolate(rebx, &mut interp, times[k]);
        assert_eq!(
            got.to_bits(),
            values[k].to_bits(),
            "spline at knot {} (t = {}) on the backward sweep: got {:?}, tabulated {:?}",
            k,
            times[k],
            got,
            values[k]
        );
    }
    rebx_free_interpolator(interp);
}

/// Between the knots the natural cubic spline must reproduce a smooth
/// function. `sin` on `[0, 2*pi]` is the fair test case: its second
/// derivative vanishes at both endpoints, which is exactly the "natural"
/// boundary condition `rebx_spline` imposes, so the classical interior error
/// bound `(5/384) h^4 max|f''''|` holds over the whole interval. With 201
/// knots, h = 2*pi/200 = 0.0314, that bound is 1.3e-8; we assert 1e-7.
/// We also compare against linear interpolation on the very same knots
/// (error h^2/8 = 1.2e-4): a cubic spline must beat it by orders of
/// magnitude, which no constant-return or off-by-one-interval implementation
/// could.
#[test]
fn interpolator_reproduces_a_smooth_function_between_knots() {
    let n = 201;
    let (times, values) = sample(0., 2. * PI, n, |t| t.sin());
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");
    let mut interp = rebx_create_interpolator(
        rebx,
        n as i32,
        &times,
        &values,
        rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
    );

    let nq = 1997;
    let mut max_spline_err = 0.0f64;
    let mut max_linear_err = 0.0f64;
    let mut worst_t = 0.0f64;
    for q in 0..nq {
        let t = 2. * PI * ((q as f64) + 0.5) / (nq as f64);
        let got = rebx_interpolate(rebx, &mut interp, t);
        let err = (got - t.sin()).abs();
        if err > max_spline_err {
            max_spline_err = err;
            worst_t = t;
        }
        // independent linear interpolation on the same table
        let h = times[1] - times[0];
        let k = ((t - times[0]) / h).floor() as usize;
        let k = k.min(n - 2);
        let w = (t - times[k]) / (times[k + 1] - times[k]);
        let lin = (1. - w) * values[k] + w * values[k + 1];
        max_linear_err = max_linear_err.max((lin - t.sin()).abs());
    }

    assert!(
        max_spline_err < 1e-7,
        "max |spline - sin| over (0, 2pi) with {} knots = {:e} (worst at t = {}), \
         expected below the natural-spline bound 1e-7",
        n,
        max_spline_err,
        worst_t
    );
    assert!(
        max_spline_err * 100. < max_linear_err,
        "spline error {:e} is not at least 100x smaller than linear-interpolation \
         error {:e} on the identical knots",
        max_spline_err,
        max_linear_err
    );
    rebx_free_interpolator(interp);
}

/// Fourth-order convergence: halving the knot spacing must cut the maximum
/// interpolation error by about 2^4 = 16. This is the defining property of a
/// cubic spline and cannot be faked by a lower-order rule (linear would give
/// 4, quadratic 8). Same function and boundary-condition argument as above.
#[test]
fn interpolator_error_scales_as_the_fourth_power_of_knot_spacing() {
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");

    let mut max_err_for = |n: usize| -> f64 {
        let (times, values) = sample(0., 2. * PI, n, |t| t.sin());
        let mut interp = rebx_create_interpolator(
            rebx,
            n as i32,
            &times,
            &values,
            rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
        );
        let nq = 1997;
        let mut m = 0.0f64;
        for q in 0..nq {
            let t = 2. * PI * ((q as f64) + 0.5) / (nq as f64);
            let got = rebx_interpolate(rebx, &mut interp, t);
            m = m.max((got - t.sin()).abs());
        }
        rebx_free_interpolator(interp);
        m
    };

    let coarse = max_err_for(51); // h = 2pi/50
    let fine = max_err_for(101); // h = 2pi/100
    let ratio = coarse / fine;
    assert!(
        (12.0..20.0).contains(&ratio),
        "error ratio on halving h is {} (coarse {:e}, fine {:e}); a cubic spline must \
         give ~16, linear would give ~4",
        ratio,
        coarse,
        fine
    );
}

/// A natural cubic spline through samples of a straight line reproduces that
/// line exactly: the tridiagonal right-hand side `u` is identically zero
/// (equal successive slopes), so every `y2` is zero and `rebx_splint`
/// degenerates to `a*y[klo] + b*y[klo+1]` with `a + b == 1`. With integer
/// knots and unit spacing this is exact in binary floating point.
#[test]
fn interpolator_reproduces_a_straight_line() {
    let n = 11;
    let times: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let values: Vec<f64> = times.iter().map(|t| 3. * t - 2.).collect();
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");
    let mut interp = rebx_create_interpolator(
        rebx,
        n as i32,
        &times,
        &values,
        rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
    );

    for y2 in interp.y2.iter() {
        assert!(
            y2.abs() < 1e-300,
            "second derivative of a spline through a straight line should be 0, got {:e}",
            y2
        );
    }
    for q in 0..=40 {
        let t = (q as f64) * 0.25;
        let got = rebx_interpolate(rebx, &mut interp, t);
        let want = 3. * t - 2.;
        assert!(
            (got - want).abs() < 1e-14,
            "spline through the line y = 3t - 2 at t = {}: got {}, want {} (diff {:e})",
            t,
            got,
            want,
            got - want
        );
    }
    rebx_free_interpolator(interp);
}

/// `REBX_INTERPOLATION_NONE` is the C's `return 0.;` branch — it must return
/// a true zero regardless of the table it was built from, and must not touch
/// the `klo` cache or allocate `y2`.
#[test]
fn interpolator_none_returns_exact_zero_and_builds_no_second_derivatives() {
    let (times, values) = sample(0., 10., 13, |t| 5. + t * t);
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");
    let mut interp = rebx_create_interpolator(
        rebx,
        times.len() as i32,
        &times,
        &values,
        rebx_interpolation_type::REBX_INTERPOLATION_NONE,
    );
    assert!(
        interp.y2.is_empty(),
        "REBX_INTERPOLATION_NONE must leave y2 unallocated (C: y2 == NULL), got {} entries",
        interp.y2.len()
    );
    for q in 0..20 {
        let t = 0.5 * (q as f64);
        let got = rebx_interpolate(rebx, &mut interp, t);
        assert_eq!(
            got.to_bits(),
            0.0f64.to_bits(),
            "REBX_INTERPOLATION_NONE at t = {} returned {} (bits {:016x}), expected +0.0",
            t,
            got,
            got.to_bits()
        );
    }
    assert_eq!(
        interp.klo, 0,
        "REBX_INTERPOLATION_NONE must not advance the klo cache, got klo = {}",
        interp.klo
    );
    rebx_free_interpolator(interp);
}

/// Two interpolators built from identical tables, queried in the identical
/// order, must agree bit-for-bit; and a forward sweep and a backward sweep
/// over the same query points must also agree bit-for-bit, because the `klo`
/// cache only selects an interval and the arithmetic in that interval is
/// independent of how the cache got there.
#[test]
fn interpolator_is_bit_deterministic_and_sweep_direction_independent() {
    let (times, values) = sample(0., 12., 33, |t| (0.6 * t).sin() * (1. + 0.1 * t));
    let mut sim = attached_sim();
    let rebx = rebx_extras_mut(&mut sim).expect("extras");

    let build = |rebx: &mut rebx_extras| {
        rebx_create_interpolator(
            rebx,
            times.len() as i32,
            &times,
            &values,
            rebx_interpolation_type::REBX_INTERPOLATION_SPLINE,
        )
    };
    let mut a = build(rebx);
    let mut b = build(rebx);
    let mut c = build(rebx);

    let queries: Vec<f64> = (0..500).map(|q| 12. * (q as f64) / 499.).collect();
    let fwd_a: Vec<f64> = queries.iter().map(|&t| rebx_interpolate(rebx, &mut a, t)).collect();
    let fwd_b: Vec<f64> = queries.iter().map(|&t| rebx_interpolate(rebx, &mut b, t)).collect();
    for (i, (&x, &y)) in fwd_a.iter().zip(fwd_b.iter()).enumerate() {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "two identically-built interpolators disagree at query {} (t = {}): {:016x} vs {:016x}",
            i,
            queries[i],
            x.to_bits(),
            y.to_bits()
        );
    }

    // walk c to the top, then sweep down through the same points
    let _ = rebx_interpolate(rebx, &mut c, queries[queries.len() - 1]);
    let mut bwd: Vec<f64> = queries.iter().rev().map(|&t| rebx_interpolate(rebx, &mut c, t)).collect();
    bwd.reverse();
    for (i, (&x, &y)) in fwd_a.iter().zip(bwd.iter()).enumerate() {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "forward and backward sweeps disagree at query {} (t = {}): {:016x} vs {:016x}",
            i,
            queries[i],
            x.to_bits(),
            y.to_bits()
        );
    }

    rebx_free_interpolator(a);
    rebx_free_interpolator(b);
    rebx_free_interpolator(c);
}

// ===========================================================================
// tides_spin helpers: spin angular momentum
// ===========================================================================

/// `rebx_tools_spin_angular_momentum` is by definition `sum_i I_i * Omega_i`
/// over the bodies that carry BOTH parameters. Verified against an
/// independent sum, and against the linearity that definition implies:
/// scaling every moment of inertia by two scales L by exactly two (a factor
/// of two is exact in binary floating point).
#[test]
fn spin_angular_momentum_is_the_sum_of_I_times_Omega() {
    let mut sim = attached_sim();
    for _ in 0..3 {
        reb_simulation_add_fmt(&mut sim, "m", &[reb_fmt_arg::d(1e-3)]);
    }

    let I = [2.5e-7, 7.25e-8, 1.5e-6];
    let Om = [
        reb_vec3d { x: 11., y: -3.5, z: 240. },
        reb_vec3d { x: -60., y: 7., z: 19.5 },
        reb_vec3d { x: 0.75, y: 130., z: -8. },
    ];

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        for i in 0..3 {
            rebx_set_param_double(rebx, rebx_ap::particle(i), "I", I[i]);
            rebx_set_param_vec3d(rebx, rebx_ap::particle(i), "Omega", Om[i]);
        }
    }

    let L = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx)).unwrap();
    let want = reb_vec3d {
        x: I[0] * Om[0].x + I[1] * Om[1].x + I[2] * Om[2].x,
        y: I[0] * Om[0].y + I[1] * Om[1].y + I[2] * Om[2].y,
        z: I[0] * Om[0].z + I[1] * Om[1].z + I[2] * Om[2].z,
    };
    for (name, got, exp) in [("x", L.x, want.x), ("y", L.y, want.y), ("z", L.z, want.z)] {
        assert!(
            (got - exp).abs() <= 1e-15 * exp.abs(),
            "spin L{} = {:e}, independent sum I*Omega = {:e} (rel diff {:e})",
            name,
            got,
            exp,
            (got - exp).abs() / exp.abs()
        );
    }

    // linearity in I: doubling every I doubles L exactly
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        for i in 0..3 {
            rebx_set_param_double(rebx, rebx_ap::particle(i), "I", 2. * I[i]);
        }
    }
    let L2 = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx)).unwrap();
    for (name, got, exp) in [
        ("x", L2.x, 2. * L.x),
        ("y", L2.y, 2. * L.y),
        ("z", L2.z, 2. * L.z),
    ] {
        assert_eq!(
            got.to_bits(),
            exp.to_bits(),
            "doubling every I must double spin L{} exactly: got {:e}, expected {:e}",
            name,
            got,
            exp
        );
    }
}

/// The C guards the accumulation with `if (Omega != NULL && I != NULL)`: a
/// body carrying only one of the two contributes nothing. Checked by adding
/// the half-specified bodies to a simulation and confirming L is unchanged.
#[test]
fn spin_angular_momentum_skips_bodies_missing_I_or_Omega() {
    let mut sim = attached_sim();
    for _ in 0..4 {
        reb_simulation_add_fmt(&mut sim, "m", &[reb_fmt_arg::d(1e-3)]);
    }
    let Om0 = reb_vec3d { x: 3., y: -4., z: 12. };
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        // fully specified
        rebx_set_param_double(rebx, rebx_ap::particle(0), "I", 5e-7);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(0), "Omega", Om0);
        // Omega but no I
        rebx_set_param_vec3d(
            rebx,
            rebx_ap::particle(1),
            "Omega",
            reb_vec3d { x: 1e3, y: 1e3, z: 1e3 },
        );
        // I but no Omega
        rebx_set_param_double(rebx, rebx_ap::particle(2), "I", 9e-3);
        // neither: particle 3
    }

    let L = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx)).unwrap();
    let want = reb_vec3d { x: 5e-7 * Om0.x, y: 5e-7 * Om0.y, z: 5e-7 * Om0.z };
    assert_eq!(
        (L.x.to_bits(), L.y.to_bits(), L.z.to_bits()),
        (want.x.to_bits(), want.y.to_bits(), want.z.to_bits()),
        "only the fully-specified body may contribute: L = ({:e}, {:e}, {:e}), \
         expected ({:e}, {:e}, {:e})",
        L.x,
        L.y,
        L.z,
        want.x,
        want.y,
        want.z
    );

    // sanity on the magnitude of the half-specified body: were it counted,
    // L would be swamped by its Omega = (1e3, 1e3, 1e3).
    assert!(
        vlen(L) < 1e-4,
        "|L| = {:e} is far larger than 5e-7*13 = 6.5e-6, so a half-specified body \
         was counted",
        vlen(L)
    );
}

// ===========================================================================
// tides_spin helpers: rebx_simulation_irotate
// ===========================================================================

/// A rigid rotation of the whole system cannot change the magnitude of the
/// total (orbital + spin) angular momentum. `rebx_simulation_irotate` rotates
/// orbits through REBOUND and then the `Omega` vectors itself, so if either
/// half were skipped or rotated by a different quaternion, the magnitude of
/// the *sum* would move even though each half's magnitude stayed put. The
/// spin contribution here is a few percent of the orbital one, so the test
/// really does bite.
#[test]
fn irotate_leaves_the_total_angular_momentum_magnitude_invariant() {
    let mut sim = hot_jupiter(true, true);
    let L0 = total_L(&mut sim);
    let L_spin0 = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx))
        .unwrap();
    assert!(
        vlen(L_spin0) > 1e-3 * vlen(L0),
        "test is vacuous: spin |L| = {:e} is negligible next to total |L| = {:e}",
        vlen(L_spin0),
        vlen(L0)
    );

    let q = reb_rotation_init_angle_axis(0.7391, reb_vec3d { x: 0.3, y: -0.5, z: 0.8 });
    rebx_with(&mut sim, |sim, rebx| rebx_simulation_irotate(sim, rebx, q));
    let L1 = total_L(&mut sim);

    let rel = (vlen(L1) - vlen(L0)).abs() / vlen(L0);
    assert!(
        rel < 1e-13,
        "|L_tot| changed under rotation: before {:e}, after {:e} (rel {:e})",
        vlen(L0),
        vlen(L1),
        rel
    );
}

/// Stronger than the magnitude test: because a proper rotation R satisfies
/// R(a x b) = (Ra) x (Rb), the total angular momentum must transform as a
/// vector, L -> R L, component by component. Compared against R applied
/// directly to the pre-rotation L.
#[test]
fn irotate_transforms_the_total_angular_momentum_as_a_vector() {
    let mut sim = hot_jupiter(true, true);
    let L0 = total_L(&mut sim);

    let q = reb_rotation_init_angle_axis(-1.234, reb_vec3d { x: -0.2, y: 0.9, z: 0.35 });
    let want = reb_vec3d_rotate(L0, q);

    rebx_with(&mut sim, |sim, rebx| rebx_simulation_irotate(sim, rebx, q));
    let L1 = total_L(&mut sim);

    let scale = vlen(L0);
    for (name, got, exp) in [("x", L1.x, want.x), ("y", L1.y, want.y), ("z", L1.z, want.z)] {
        assert!(
            (got - exp).abs() < 1e-13 * scale,
            "rotated total L{} = {:e}, but R applied to the original L gives {:e} \
             (diff {:e}, |L| = {:e})",
            name,
            got,
            exp,
            got - exp,
            scale
        );
    }
}

/// Rotating by `q` and then by `q^-1` is the identity, so every position,
/// velocity and spin vector must come back to where it started. This is the
/// round trip that catches a spin vector rotated by the wrong quaternion (or
/// not at all): an un-rotated Omega would return unchanged by luck, so the
/// test also checks that the intermediate state really did move.
#[test]
fn irotate_round_trip_restores_positions_velocities_and_spins() {
    let mut sim = hot_jupiter(true, true);

    let p_before: Vec<reb_particle> = sim.particles[..sim.N].to_vec();
    let om_before: Vec<reb_vec3d> = (0..sim.N)
        .map(|i| {
            rebx_extras_ref(&sim)
                .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega"))
                .expect("Omega set")
        })
        .collect();

    let q = reb_rotation_init_angle_axis(2.0, reb_vec3d { x: 0.1, y: 0.2, z: -0.97 });
    rebx_with(&mut sim, |sim, rebx| rebx_simulation_irotate(sim, rebx, q));

    // the rotation must actually have done something to the spins
    let om_mid: Vec<reb_vec3d> = (0..sim.N)
        .map(|i| {
            rebx_extras_ref(&sim)
                .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega"))
                .unwrap()
        })
        .collect();
    let moved = om_mid
        .iter()
        .zip(om_before.iter())
        .map(|(a, b)| vlen(reb_vec3d { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }))
        .fold(0.0f64, f64::max);
    let om_scale = om_before.iter().map(|o| vlen(*o)).fold(0.0f64, f64::max);
    assert!(
        moved > 0.1 * om_scale,
        "a 2 rad rotation barely moved the spin vectors (max displacement {:e}, \
         |Omega| scale {:e}) — were they rotated at all?",
        moved,
        om_scale
    );

    rebx_with(&mut sim, |sim, rebx| {
        rebx_simulation_irotate(sim, rebx, reb_rotation_inverse(q))
    });

    for i in 0..sim.N {
        let a = sim.particles[i];
        let b = p_before[i];
        let rscale = (b.x * b.x + b.y * b.y + b.z * b.z).sqrt().max(1e-30);
        let vscale = (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz).sqrt().max(1e-30);
        for (name, got, exp, sc) in [
            ("x", a.x, b.x, rscale),
            ("y", a.y, b.y, rscale),
            ("z", a.z, b.z, rscale),
            ("vx", a.vx, b.vx, vscale),
            ("vy", a.vy, b.vy, vscale),
            ("vz", a.vz, b.vz, vscale),
        ] {
            assert!(
                (got - exp).abs() < 1e-13 * sc,
                "particle {} {} after q then q^-1: got {:e}, started at {:e} (diff {:e})",
                i,
                name,
                got,
                exp,
                got - exp
            );
        }

        let om = rebx_extras_ref(&sim)
            .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega"))
            .unwrap();
        let sc = vlen(om_before[i]).max(1e-30);
        for (name, got, exp) in [
            ("Omega.x", om.x, om_before[i].x),
            ("Omega.y", om.y, om_before[i].y),
            ("Omega.z", om.z, om_before[i].z),
        ] {
            assert!(
                (got - exp).abs() < 1e-13 * sc,
                "particle {} {} after q then q^-1: got {:e}, started at {:e} (diff {:e})",
                i,
                name,
                got,
                exp,
                got - exp
            );
        }
    }
}

// ===========================================================================
// tides_spin: rebx_spin_initialize_ode
// ===========================================================================

/// The spin ODE holds three state variables per tracked body, and a body is
/// tracked exactly when it has both `I` and `Omega`. The state vector must
/// also be seeded with those spin vectors in particle order (the C leaves
/// that to `rebx_spin_sync_pre`; this port does it at creation — see the
/// module notes — so the values integrated from are the same either way).
/// Re-initializing must replace the old ODE rather than accumulate ODEs, and
/// with nothing to track no ODE may be created at all.
#[test]
fn spin_initialize_ode_has_three_state_variables_per_tracked_body() {
    let mut sim = attached_sim();
    reb_simulation_set_integrator(&mut sim, "ias15");
    for _ in 0..3 {
        reb_simulation_add_fmt(&mut sim, "m r", &[reb_fmt_arg::d(1e-3), reb_fmt_arg::d(1e-4)]);
    }
    let effect = rebx_load_force(&mut sim, "tides_spin").expect("tides_spin force");
    rebx_add_force(&mut sim, effect);

    let om = [
        reb_vec3d { x: 1.5, y: -2.5, z: 30.25 },
        reb_vec3d { x: -7., y: 0.125, z: 11. },
        reb_vec3d { x: 4., y: 4.5, z: -60. },
    ];

    // Nothing tracked yet -> no ODE at all.
    rebx_with(&mut sim, |sim, rebx| rebx_spin_initialize_ode(sim, rebx, effect));
    assert!(
        sim.odes.is_empty(),
        "no body has I and Omega, so no spin ODE may be created; found {} ODE(s)",
        sim.odes.len()
    );

    // Two tracked bodies (particle 2 gets I only) -> length 6.
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "I", 1e-9);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(0), "Omega", om[0]);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "I", 2e-9);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(1), "Omega", om[1]);
        rebx_set_param_double(rebx, rebx_ap::particle(2), "I", 3e-9);
    }
    rebx_with(&mut sim, |sim, rebx| rebx_spin_initialize_ode(sim, rebx, effect));
    assert_eq!(
        sim.odes.len(),
        1,
        "expected exactly one spin ODE, found {}",
        sim.odes.len()
    );
    assert_eq!(
        sim.odes[0].length, 6,
        "2 tracked bodies x 3 spin components = 6 state variables, got {}",
        sim.odes[0].length
    );
    let want: Vec<f64> = vec![om[0].x, om[0].y, om[0].z, om[1].x, om[1].y, om[1].z];
    for (k, (&got, &exp)) in sim.odes[0].y.iter().zip(want.iter()).enumerate() {
        assert_eq!(
            got.to_bits(),
            exp.to_bits(),
            "spin ODE state[{}] = {:e}, expected the Omega component {:e}",
            k,
            got,
            exp
        );
    }

    // Re-initializing replaces, never accumulates; a third tracked body
    // lengthens the state vector to 9.
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_vec3d(rebx, rebx_ap::particle(2), "Omega", om[2]);
    }
    rebx_with(&mut sim, |sim, rebx| rebx_spin_initialize_ode(sim, rebx, effect));
    assert_eq!(
        sim.odes.len(),
        1,
        "re-initializing must free the previous spin ODE, found {} ODEs",
        sim.odes.len()
    );
    assert_eq!(
        sim.odes[0].length, 9,
        "3 tracked bodies x 3 spin components = 9 state variables, got {}",
        sim.odes[0].length
    );
    assert_eq!(
        sim.odes[0].y.len(),
        9,
        "state vector length {} does not match the declared ODE length 9",
        sim.odes[0].y.len()
    );
    for (k, &exp) in [om[2].x, om[2].y, om[2].z].iter().enumerate() {
        let got = sim.odes[0].y[6 + k];
        assert_eq!(
            got.to_bits(),
            exp.to_bits(),
            "spin ODE state[{}] = {:e}, expected the third body's Omega component {:e}",
            6 + k,
            got,
            exp
        );
    }
}

// ===========================================================================
// tides_spin: end-to-end behaviour
// ===========================================================================

/// Tidal forces and torques are internal: the acceleration applied to the
/// target is `-(ms/mtot) F` and to the source `+(mt/mtot) F`, and the spin
/// ODE removes from the spin exactly the angular momentum that the same `F`
/// deposits into the orbit. Total angular momentum (orbital + spin) is
/// therefore an exact constant of the motion, whatever the dissipation does
/// to the energy — the only thing that spoils it numerically is the
/// first-order operator splitting between the WHFast orbit step and the BS
/// spin step.
///
/// So the assertion is not just "the drift is small" (a number that would
/// depend on the timestep chosen here) but "the drift is proportional to the
/// timestep": halving `dt` must halve it. A missing back reaction, a wrong
/// sign, or a torque that is not the moment of the applied force would leave
/// a residue that does NOT vanish with `dt`.
#[test]
fn tides_spin_conserves_total_angular_momentum_in_the_small_timestep_limit() {
    let tmax = 1.0;
    let drift_at = |dt: f64| -> (f64, f64) {
        let mut sim = hot_jupiter_dt(true, true, dt);
        let L0 = total_L(&mut sim);
        let s0 = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx))
            .unwrap();
        reb_simulation_integrate(&mut sim, tmax);
        let L1 = total_L(&mut sim);
        let s1 = rebx_with(&mut sim, |sim, rebx| rebx_tools_spin_angular_momentum(sim, rebx))
            .unwrap();
        let dL = vlen(reb_vec3d { x: L1.x - L0.x, y: L1.y - L0.y, z: L1.z - L0.z }) / vlen(L0);
        let dS = vlen(reb_vec3d { x: s1.x - s0.x, y: s1.y - s0.y, z: s1.z - s0.z }) / vlen(L0);
        (dL, dS)
    };

    let (d_coarse, s_coarse) = drift_at(1e-3);
    let (d_fine, _) = drift_at(5e-4);

    assert!(
        d_coarse < 1e-6,
        "total (orbit + spin) angular momentum drifted by a relative {:e} over t = {} \
         at dt = 1e-3, which is far more than operator splitting can explain",
        d_coarse,
        tmax
    );
    let ratio = d_coarse / d_fine;
    assert!(
        (1.7..2.3).contains(&ratio),
        "halving dt changed the angular-momentum drift by a factor {} (dt = 1e-3: {:e}, \
         dt = 5e-4: {:e}); first-order splitting of an exactly conserved quantity \
         must give ~2",
        ratio,
        d_coarse,
        d_fine
    );

    // The exchange must be real, not a frozen spin: the spin angular momentum
    // has to move by far more than the total-L bookkeeping error, otherwise
    // conservation would be trivially satisfied.
    assert!(
        s_coarse > 100. * d_coarse,
        "spin angular momentum barely moved (relative change {:e}) compared with the \
         total-L drift {:e}; the spin ODE does not appear to be evolving",
        s_coarse,
        d_coarse
    );
}

/// `rebx_tides_spin` skips any body without a `k2` parameter, so a simulation
/// carrying the force but no Love numbers must integrate as pure Newtonian
/// gravity. Compared bit-for-bit against the same initial conditions run
/// without REBOUNDx at all: if the force ever wrote a non-zero acceleration
/// the trajectories would separate.
#[test]
fn tides_spin_without_k2_reproduces_plain_gravity_bit_for_bit() {
    let tmax = 0.5;

    let mut plain = reb_simulation_create();
    plain.save_messages = 1;
    reb_simulation_add_fmt(&mut plain, "m r", &[reb_fmt_arg::d(1.), reb_fmt_arg::d(0.00465)]);
    reb_simulation_add_fmt(
        &mut plain,
        "m a e inc r",
        &[
            reb_fmt_arg::d(9.55e-4),
            reb_fmt_arg::d(0.04072),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(4.676e-4),
        ],
    );
    plain.N_active = 2;
    reb_simulation_set_integrator(&mut plain, "ias15");
    reb_simulation_move_to_com(&mut plain);
    reb_simulation_integrate(&mut plain, tmax);

    let mut withx = reb_simulation_create();
    withx.save_messages = 1;
    reb_simulation_add_fmt(&mut withx, "m r", &[reb_fmt_arg::d(1.), reb_fmt_arg::d(0.00465)]);
    reb_simulation_add_fmt(
        &mut withx,
        "m a e inc r",
        &[
            reb_fmt_arg::d(9.55e-4),
            reb_fmt_arg::d(0.04072),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(4.676e-4),
        ],
    );
    withx.N_active = 2;
    reb_simulation_set_integrator(&mut withx, "ias15");
    rebx_attach(&mut withx);
    let effect = rebx_load_force(&mut withx, "tides_spin").expect("tides_spin force");
    rebx_add_force(&mut withx, effect);
    reb_simulation_move_to_com(&mut withx);
    reb_simulation_integrate(&mut withx, tmax);

    assert_eq!(
        withx.t.to_bits(),
        plain.t.to_bits(),
        "final times differ: {} vs {}",
        withx.t,
        plain.t
    );
    for i in 0..plain.N {
        let a = withx.particles[i];
        let b = plain.particles[i];
        for (name, got, exp) in [
            ("x", a.x, b.x),
            ("y", a.y, b.y),
            ("z", a.z, b.z),
            ("vx", a.vx, b.vx),
            ("vy", a.vy, b.vy),
            ("vz", a.vz, b.vz),
        ] {
            assert_eq!(
                got.to_bits(),
                exp.to_bits(),
                "particle {} {} with a k2-less tides_spin force = {:e} ({:016x}), \
                 plain gravity = {:e} ({:016x})",
                i,
                name,
                got,
                got.to_bits(),
                exp,
                exp.to_bits()
            );
        }
    }
}

/// Same setup, same code path, twice: the spin-tides integration must be
/// bit-reproducible, spins included. Any dependence on hash iteration order
/// in the parameter lookups, or on uninitialized state, would show up here.
#[test]
fn tides_spin_integration_is_bit_reproducible() {
    let run = || {
        let mut sim = hot_jupiter(true, true);
        reb_simulation_integrate(&mut sim, 0.3);
        let mut out: Vec<u64> = vec![sim.t.to_bits()];
        for i in 0..sim.N {
            let p = sim.particles[i];
            for v in [p.x, p.y, p.z, p.vx, p.vy, p.vz] {
                out.push(v.to_bits());
            }
            let om = rebx_extras_ref(&sim)
                .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega"))
                .expect("Omega set");
            for v in [om.x, om.y, om.z] {
                out.push(v.to_bits());
            }
        }
        out
    };

    let a = run();
    let b = run();
    assert_eq!(
        a.len(),
        b.len(),
        "the two runs produced different numbers of state variables: {} vs {}",
        a.len(),
        b.len()
    );
    for (k, (&x, &y)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(
            x, y,
            "state variable {} differs between two identical runs: {:016x} ({:e}) vs \
             {:016x} ({:e})",
            k,
            x,
            f64::from_bits(x),
            y,
            f64::from_bits(y)
        );
    }
}

/// A super-synchronous, dissipative hot Jupiter must spin DOWN: its spin
/// period (0.5 d) is far shorter than its orbital period (~3 d), so the tidal
/// lag torque removes rotational angular momentum until pseudo-synchronism.
/// The sign of that change is a physical prediction independent of any number
/// in the library, and it flips if the lag enters with the wrong sign.
#[test]
fn dissipative_tides_spin_down_a_super_synchronous_planet() {
    let mut sim = hot_jupiter(true, true);

    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    let P_orb = 2. * PI / orb.n;
    let om0 = rebx_extras_ref(&sim)
        .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(1), "Omega"))
        .unwrap();
    let P_spin = 2. * PI / vlen(om0);
    assert!(
        P_spin < P_orb,
        "setup is not super-synchronous: spin period {:e} vs orbital period {:e}",
        P_spin,
        P_orb
    );

    reb_simulation_integrate(&mut sim, 2.0);

    let om1 = rebx_extras_ref(&sim)
        .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(1), "Omega"))
        .unwrap();
    assert!(
        vlen(om1) < vlen(om0),
        "a super-synchronous dissipative planet must spin down: |Omega| went from \
         {:e} to {:e} (change {:e})",
        vlen(om0),
        vlen(om1),
        vlen(om1) - vlen(om0)
    );
}
