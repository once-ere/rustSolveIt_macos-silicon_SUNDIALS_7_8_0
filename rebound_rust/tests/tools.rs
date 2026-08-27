//! Integration tests for the tools module group of rebound_rs.
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

// ---------------------------------------------------------------------------
// small local helpers (nothing here is part of the library under test)
// ---------------------------------------------------------------------------

/// Relative difference |a-b|/max(|b|,tiny).
fn relerr(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// Shorthand for a numeric `reb_simulation_add_fmt` vararg.
fn d(v: f64) -> reb_fmt_arg {
    reb_fmt_arg::d(v)
}

/// The 11 f64 fields of a particle, as raw bits, for exact comparisons.
fn pbits(p: &reb_particle) -> [u64; 11] {
    [
        p.x.to_bits(),
        p.y.to_bits(),
        p.z.to_bits(),
        p.vx.to_bits(),
        p.vy.to_bits(),
        p.vz.to_bits(),
        p.ax.to_bits(),
        p.ay.to_bits(),
        p.az.to_bits(),
        p.m.to_bits(),
        p.r.to_bits(),
    ]
}

/// An INDEPENDENT model of glibc's `rand_r`: the same LCG written with
/// 64-bit arithmetic and an explicit mask instead of Rust's wrapping u32
/// operations, so agreement is evidence about the algorithm and not just
/// a restatement of the implementation.
fn rand_r_model(seed: &mut u64) -> i64 {
    const M: u64 = 0xFFFF_FFFF;
    let mut next = *seed & M;
    next = (next * 1103515245 + 12345) & M;
    let mut result: i64 = ((next / 65536) % 2048) as i64;
    next = (next * 1103515245 + 12345) & M;
    result = (result * 1024) ^ (((next / 65536) % 1024) as i64);
    next = (next * 1103515245 + 12345) & M;
    result = (result * 1024) ^ (((next / 65536) % 1024) as i64);
    *seed = next;
    result
}

/// Angular difference of two angles, taken to the shorter way round.
fn ang_diff(a: f64, b: f64) -> f64 {
    let mut dd = (a - b) % (2. * PI);
    if dd > PI {
        dd -= 2. * PI;
    }
    if dd < -PI {
        dd += 2. * PI;
    }
    dd.abs()
}

/// Total energy computed a second, completely independent way (an O(N^2)
/// double sum written from the physics, not from tools.rs).
fn energy_by_hand(r: &reb_simulation) -> f64 {
    let mut ek = 0.;
    let mut ep = 0.;
    for i in 0..r.N {
        let p = r.particles[i];
        ek += 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
        for j in 0..r.N {
            if j == i {
                continue;
            }
            let q = r.particles[j];
            let dx = p.x - q.x;
            let dy = p.y - q.y;
            let dz = p.z - q.z;
            ep -= 0.5 * r.G * p.m * q.m / (dx * dx + dy * dy + dz * dz).sqrt();
        }
    }
    ek + ep
}

/// A Sun + two planets test system (G = 1, masses in solar units).
fn three_body() -> reb_simulation {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
    reb_simulation_add_fmt(
        &mut r,
        "m a e inc Omega omega f",
        &[d(1e-3), d(1.0), d(0.05), d(0.02), d(0.3), d(0.7), d(1.1)],
    );
    reb_simulation_add_fmt(
        &mut r,
        "m a e inc Omega omega f",
        &[d(3e-4), d(2.4), d(0.09), d(0.04), d(1.3), d(2.7), d(4.1)],
    );
    reb_simulation_move_to_com(&mut r);
    r
}

// ===========================================================================
// 1. rand_r / the random number generators
// ===========================================================================

#[test]
fn rand_r_reproduces_the_glibc_lcg_exactly() {
    // Anchors computed from s <- (s*1103515245 + 12345) mod 2^32 in exact
    // integer arithmetic, three rounds per draw, result composed as
    // ((r1<<10)^r2)<<10)^r3 with r1 = (s/65536)%2048 and r2,r3 = (s/65536)%1024.
    let expected_from_42: [i32; 8] = [
        681191333, 928546885, 1457394273, 941445650, 2129613237, 1661015563, 2071432601, 222443696,
    ];
    let expected_seed_after_42: [u32; 8] = [
        3148160401, 2219150180, 1314989459, 1108520142, 1974836613, 2852860648, 4120479207,
        3970988082,
    ];
    let mut s: u32 = 42;
    for i in 0..8 {
        let v = rand_r(&mut s);
        assert_eq!(
            v, expected_from_42[i],
            "rand_r draw {} from seed 42: got {}, expected {}",
            i, v, expected_from_42[i]
        );
        assert_eq!(
            s, expected_seed_after_42[i],
            "rand_r seed state after draw {}",
            i
        );
    }

    // The first draw's top 11 bits are the first LCG round's (s/65536)%2048.
    // 42*1103515245 + 12345 = 3397979675 (mod 2^32); 3397979675/65536 = 51849;
    // 51849 % 2048 = 649; and 681191333 >> 20 == 649.
    assert_eq!(
        expected_from_42[0] >> 20,
        649,
        "top 11 bits of the first draw must be the first LCG round"
    );
}

#[test]
fn rand_r_matches_an_independent_64bit_model_over_many_seeds() {
    for &start in &[0u32, 1, 42, 12345, 2147483647, 4294967295, 3141592653] {
        let mut a: u32 = start;
        let mut b: u64 = start as u64;
        for i in 0..64 {
            let va = rand_r(&mut a) as i64;
            let vb = rand_r_model(&mut b);
            assert_eq!(
                va, vb,
                "rand_r vs independent 64-bit LCG model, seed {} draw {}",
                start, i
            );
            assert_eq!(a as u64, b, "seed state, seed {} draw {}", start, i);
            assert!(
                va >= 0 && va <= REB_RAND_MAX as i64,
                "rand_r output {} out of [0, REB_RAND_MAX]",
                va
            );
        }
    }
}

#[test]
fn rand_r_is_a_pure_function_of_its_seed() {
    // Same seed in, same value and same successor state out.
    let mut s1: u32 = 987654321;
    let mut s2: u32 = 987654321;
    let v1 = rand_r(&mut s1);
    let v2 = rand_r(&mut s2);
    assert_eq!(v1, v2, "rand_r must be deterministic in its seed");
    assert_eq!(s1, s2, "rand_r successor state must be deterministic");
    // and the seed really advances (the generator is not stuck)
    assert_ne!(s1, 987654321, "rand_r must advance the seed");
}

#[test]
fn random_uniform_is_the_exact_documented_expression() {
    // reb_random_uniform == rand_r(seed)/REB_RAND_MAX*(max-min)+min, bit for bit.
    let mut r = reb_simulation_create();
    r.rand_seed = 20240101;
    let mut s: u32 = 20240101;
    for i in 0..200 {
        let (min, max) = (-3.5_f64, 11.25_f64);
        let got = reb_random_uniform(Some(&mut r), min, max);
        let want = (rand_r(&mut s) as f64) / (REB_RAND_MAX as f64) * (max - min) + min;
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "reb_random_uniform draw {}: {} vs {}",
            i,
            got,
            want
        );
        assert!(
            got >= min && got <= max,
            "uniform draw {} = {} outside [{}, {}]",
            i,
            got,
            min,
            max
        );
        assert_eq!(r.rand_seed, s, "simulation rand_seed must track rand_r");
    }
}

#[test]
fn random_uniform_covers_its_range_evenly() {
    let mut r = reb_simulation_create();
    r.rand_seed = 5;
    let n = 200_000;
    let (min, max) = (2.0_f64, 6.0_f64);
    let mut sum = 0.0;
    let mut bins = [0usize; 8];
    for _ in 0..n {
        let v = reb_random_uniform(Some(&mut r), min, max);
        sum += v;
        let b = (((v - min) / (max - min)) * 8.0) as usize;
        bins[b.min(7)] += 1;
    }
    let mean = sum / (n as f64);
    // sigma of the mean = (max-min)/sqrt(12 n) = 4/sqrt(2.4e6) = 2.6e-3
    assert!(
        (mean - 4.0).abs() < 0.02,
        "uniform(2,6) sample mean {} should be near 4",
        mean
    );
    for (i, &c) in bins.iter().enumerate() {
        let frac = (c as f64) / (n as f64);
        assert!(
            (frac - 0.125).abs() < 0.01,
            "uniform octile {} holds fraction {}, expected ~0.125",
            i,
            frac
        );
    }
}

#[test]
fn random_powerlaw_slope_zero_reduces_to_uniform() {
    // ((max^1 - min^1) y + min^1)^(1/1) is exactly the uniform map.
    let mut r1 = reb_simulation_create();
    let mut r2 = reb_simulation_create();
    r1.rand_seed = 777;
    r2.rand_seed = 777;
    for i in 0..100 {
        let got = reb_random_powerlaw(Some(&mut r1), 1.0, 5.0, 0.0);
        let y = reb_random_uniform(Some(&mut r2), 0.0, 1.0);
        let want = (5.0 - 1.0) * y + 1.0;
        assert!(
            relerr(got, want) < 1e-15,
            "powerlaw(slope=0) draw {}: {} vs uniform image {}",
            i,
            got,
            want
        );
        assert!(
            got >= 1.0 && got <= 5.0,
            "powerlaw draw {} = {} outside [1,5]",
            i,
            got
        );
    }
}

#[test]
fn random_powerlaw_slope_minus_one_is_log_uniform() {
    let mut r = reb_simulation_create();
    r.rand_seed = 31337;
    let (min, max) = (1.0_f64, 1000.0_f64);
    let n = 100_000;
    let mut logsum = 0.0;
    for _ in 0..n {
        let v = reb_random_powerlaw(Some(&mut r), min, max, -1.0);
        assert!(
            v >= min * (1.0 - 1e-12) && v <= max * (1.0 + 1e-12),
            "powerlaw(-1) draw {} outside [{}, {}]",
            v,
            min,
            max
        );
        logsum += v.ln();
    }
    // ln(v) must be uniform on [ln min, ln max] = [0, ln 1000]
    let mean_log = logsum / (n as f64);
    let want = 0.5 * max.ln();
    assert!(
        (mean_log - want).abs() < 0.03,
        "mean ln(v) = {} should be {} for a log-uniform draw",
        mean_log,
        want
    );
}

#[test]
fn random_powerlaw_slope_two_has_the_right_cdf() {
    // For dN/dx ~ x^slope on [min,max], the inverse-CDF map used here gives
    // <x> = (s+1)/(s+2) * (max^(s+2)-min^(s+2))/(max^(s+1)-min^(s+1)).
    let slope = 2.0_f64;
    let (min, max) = (1.0_f64, 4.0_f64);
    let mut r = reb_simulation_create();
    r.rand_seed = 24680;
    let n = 200_000;
    let mut sum = 0.0;
    for _ in 0..n {
        let v = reb_random_powerlaw(Some(&mut r), min, max, slope);
        assert!(
            v >= min - 1e-12 && v <= max + 1e-12,
            "powerlaw(2) draw {} outside [{}, {}]",
            v,
            min,
            max
        );
        sum += v;
    }
    let mean = sum / (n as f64);
    let want = (slope + 1.) / (slope + 2.) * (max.powf(slope + 2.) - min.powf(slope + 2.))
        / (max.powf(slope + 1.) - min.powf(slope + 1.));
    assert!(
        relerr(mean, want) < 0.01,
        "powerlaw(slope=2) sample mean {} vs analytic {}",
        mean,
        want
    );
}

#[test]
fn random_normal_has_zero_mean_and_the_requested_variance() {
    let variance = 2.5_f64;
    let mut r = reb_simulation_create();
    r.rand_seed = 987;
    let n = 200_000;
    let mut sum = 0.0;
    let mut sum2 = 0.0;
    let mut sum4 = 0.0;
    for _ in 0..n {
        let v = reb_random_normal(Some(&mut r), variance);
        assert!(v.is_finite(), "reb_random_normal produced {}", v);
        sum += v;
        sum2 += v * v;
        sum4 += v * v * v * v;
    }
    let mean = sum / (n as f64);
    let var = sum2 / (n as f64) - mean * mean;
    // sigma(mean) = sqrt(variance/n) = 3.5e-3 -> 0.03 is ~8 sigma
    assert!(
        mean.abs() < 0.03,
        "reb_random_normal sample mean {} should be ~0",
        mean
    );
    assert!(
        relerr(var, variance) < 0.03,
        "reb_random_normal sample variance {} vs requested {}",
        var,
        variance
    );
    // Gaussian kurtosis: <x^4> = 3 sigma^4
    let kurt = (sum4 / (n as f64)) / (variance * variance);
    assert!(
        (kurt - 3.0).abs() < 0.15,
        "reb_random_normal kurtosis {} should be ~3 for a Gaussian",
        kurt
    );
}

#[test]
fn random_normal_scales_as_sqrt_variance() {
    // The same rand_r stream with variance v and variance 4v must differ by
    // exactly a factor 2 in the scale, because only the sqrt() factor changes.
    let mut r1 = reb_simulation_create();
    let mut r2 = reb_simulation_create();
    r1.rand_seed = 4242;
    r2.rand_seed = 4242;
    for i in 0..50 {
        let a = reb_random_normal(Some(&mut r1), 1.0);
        let b = reb_random_normal(Some(&mut r2), 4.0);
        assert!(
            relerr(b, 2.0 * a) < 1e-14,
            "normal draw {}: variance 4 gave {}, twice variance-1 draw is {}",
            i,
            b,
            2.0 * a
        );
    }
}

#[test]
fn random_rayleigh_is_the_exact_documented_expression_and_has_the_right_mean() {
    let sigma = 0.3_f64;
    let mut r = reb_simulation_create();
    r.rand_seed = 555;
    let mut s: u32 = 555;
    // exact composition against the underlying rand_r stream
    for i in 0..64 {
        let got = reb_random_rayleigh(Some(&mut r), sigma);
        let y = (rand_r(&mut s) as f64) / (REB_RAND_MAX as f64) * (1.0 - 0.0) + 0.0;
        let want = sigma * (-2.0 * y.ln()).sqrt();
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "reb_random_rayleigh draw {}: {} vs {}",
            i,
            got,
            want
        );
    }
    // statistics: <x> = sigma sqrt(pi/2), <x^2> = 2 sigma^2
    let n = 200_000;
    let mut sum = 0.0;
    let mut sum2 = 0.0;
    let mut cnt = 0.0;
    for _ in 0..n {
        let v = reb_random_rayleigh(Some(&mut r), sigma);
        if !v.is_finite() {
            continue; // y == 0 exactly (probability 2^-31) gives +inf, as in C
        }
        assert!(v >= 0.0, "Rayleigh draw {} must be non-negative", v);
        sum += v;
        sum2 += v * v;
        cnt += 1.0;
    }
    let mean = sum / cnt;
    let m2 = sum2 / cnt;
    assert!(
        relerr(mean, sigma * (PI / 2.0).sqrt()) < 0.01,
        "Rayleigh mean {} vs sigma*sqrt(pi/2) = {}",
        mean,
        sigma * (PI / 2.0).sqrt()
    );
    assert!(
        relerr(m2, 2.0 * sigma * sigma) < 0.02,
        "Rayleigh second moment {} vs 2 sigma^2 = {}",
        m2,
        2.0 * sigma * sigma
    );
}

// ===========================================================================
// 2. hashing, string compare, fp contract
// ===========================================================================

#[test]
fn reb_hash_is_djb2_with_known_values() {
    // djb2: h <- h*33 + c starting from 5381.
    assert_eq!(reb_hash(""), 5381, "djb2 of the empty string is its seed");
    assert_eq!(
        reb_hash("a"),
        5381u32.wrapping_mul(33).wrapping_add(b'a' as u32),
        "djb2('a')"
    );
    assert_eq!(reb_hash("a"), 177670, "djb2('a') = 5381*33 + 97");
    assert_eq!(reb_hash("Jupiter"), 1562838504, "djb2('Jupiter')");
    // independent model, in u64 with an explicit mod 2^32, over a long string
    let long: String = (0..500).map(|i| (b'A' + (i % 26) as u8) as char).collect();
    let mut h: u64 = 5381;
    for c in long.bytes() {
        h = (h * 33 + c as u64) % 4294967296;
    }
    assert_eq!(
        reb_hash(&long) as u64,
        h,
        "djb2 must wrap mod 2^32 like the C unsigned int"
    );
    // different strings, different hashes (djb2 has no collision here)
    assert_ne!(reb_hash("Earth"), reb_hash("Mars"));
}

#[test]
fn strcmp_ignore_whitespace_matches_the_c_contract() {
    assert_eq!(
        reb_strcmp_ignore_whitespace("solar system", "solarsystem"),
        0,
        "whitespace must be ignored on both sides"
    );
    assert_eq!(reb_strcmp_ignore_whitespace(" \t a \n b ", "ab"), 0);
    assert_eq!(reb_strcmp_ignore_whitespace("", "   "), 0);
    // first differing byte decides, as a byte difference
    assert_eq!(
        reb_strcmp_ignore_whitespace("abc", "abd"),
        (b'c' as i32) - (b'd' as i32)
    );
    // otherwise the length difference decides
    assert_eq!(reb_strcmp_ignore_whitespace("abc", "ab"), 1);
    assert_eq!(reb_strcmp_ignore_whitespace("ab", "abc"), -1);
    // the two built-in dataset names must not alias each other
    assert_ne!(
        reb_strcmp_ignore_whitespace("outer solar system", "solar system"),
        0
    );
}

#[test]
fn floating_point_contraction_is_off() {
    // Bit-exactness with the C build requires that a*b+c is NOT contracted
    // into a single fused multiply-add.
    assert_eq!(
        reb_check_fp_contract(),
        0,
        "the compiler contracted a*b+c into an FMA; results cannot be bit-exact"
    );
}

// ===========================================================================
// 3. reb_mod2pi and the Kepler solvers
// ===========================================================================

#[test]
fn mod2pi_exact_values_and_range() {
    let two_pi = 2.0 * PI;
    assert_eq!(reb_mod2pi(0.0).to_bits(), 0.0f64.to_bits(), "mod2pi(0) == 0");
    assert_eq!(
        reb_mod2pi(two_pi).to_bits(),
        0.0f64.to_bits(),
        "mod2pi(2pi) == 0"
    );
    assert_eq!(
        reb_mod2pi(-two_pi).to_bits(),
        0.0f64.to_bits(),
        "mod2pi(-2pi) == 0"
    );
    assert_eq!(
        reb_mod2pi(2.0 * two_pi).to_bits(),
        0.0f64.to_bits(),
        "mod2pi(4pi) == 0"
    );
    // (2pi + f) never needs a reduction here, so the result is exactly 2pi-0.5
    assert_eq!(
        reb_mod2pi(-0.5).to_bits(),
        (two_pi - 0.5).to_bits(),
        "mod2pi(-0.5) == 2pi - 0.5 exactly"
    );

    // range and periodicity
    let mut x = -6.2;
    while x < 12.5 {
        let m = reb_mod2pi(x);
        assert!(
            m >= 0.0 && m < two_pi,
            "mod2pi({}) = {} outside [0, 2pi)",
            x,
            m
        );
        let k = ((x - m) / two_pi).round();
        assert!(
            (x - m - k * two_pi).abs() < 1e-13,
            "mod2pi({}) = {} differs from x by {} which is not a multiple of 2pi",
            x,
            m,
            x - m
        );
        x += 0.37;
    }
}

#[test]
fn M_to_E_solves_keplers_equation_elliptic() {
    for &e in &[0.0, 1e-9, 0.1, 0.5, 0.9, 0.99] {
        let mut M = -9.0;
        while M < 9.0 {
            let E = reb_M_to_E(e, M);
            assert!(
                E >= 0.0 && E < 2.0 * PI,
                "reb_M_to_E(e={}, M={}) = {} must be reduced to [0,2pi)",
                e,
                M,
                E
            );
            // Kepler's equation, both sides reduced to [0,2pi)
            let resid = ang_diff(E - e * E.sin(), reb_mod2pi(M));
            assert!(
                resid < 1e-11,
                "Kepler residual {} for e={}, M={} (E={})",
                resid,
                e,
                M,
                E
            );
            M += 0.31;
        }
    }
    // e == 0 collapses exactly onto M
    for &M in &[0.0_f64, 1.0, 3.0, 6.0] {
        assert_eq!(
            reb_M_to_E(0.0, M).to_bits(),
            reb_mod2pi(M).to_bits(),
            "for e=0 the eccentric anomaly is the mean anomaly"
        );
    }
}

#[test]
fn M_to_E_solves_the_hyperbolic_kepler_equation() {
    for &e in &[1.05, 1.5, 3.0, 10.0] {
        for &M in &[-40.0_f64, -7.0, -0.5, 0.5, 7.0, 40.0] {
            let H = reb_M_to_E(e, M);
            // C solves F = H - e sinh H + M -> 0, i.e. e sinh H - H = M
            let lhs = e * H.sinh() - H;
            let scale = 1.0 + lhs.abs();
            assert!(
                (lhs - M).abs() < 1e-9 * scale,
                "hyperbolic Kepler residual {} for e={}, M={} (H={})",
                (lhs - M).abs(),
                e,
                M,
                H
            );
            assert!(
                H * M > 0.0,
                "H and M must share a sign (e={}, M={}, H={})",
                e,
                M,
                H
            );
        }
    }
}

#[test]
fn E_to_f_matches_the_true_anomaly_identity() {
    // Elliptic: cos f = (cos E - e)/(1 - e cos E), sin f = sqrt(1-e^2) sin E/(1-e cos E)
    for &e in &[0.0, 0.2, 0.6, 0.95] {
        let mut E = 0.05;
        while E < 2.0 * PI {
            let f = reb_E_to_f(e, E);
            let den = 1.0 - e * E.cos();
            let cf = (E.cos() - e) / den;
            let sf = (1.0 - e * e).sqrt() * E.sin() / den;
            assert!(
                (f.cos() - cf).abs() < 1e-11 && (f.sin() - sf).abs() < 1e-11,
                "E_to_f(e={}, E={}) = {}: cos/sin {} {} vs identity {} {}",
                e,
                E,
                f,
                f.cos(),
                f.sin(),
                cf,
                sf
            );
            E += 0.21;
        }
    }
    // Hyperbolic: cos f = (cosh H - e)/(1 - e cosh H)
    for &e in &[1.2, 2.0, 5.0] {
        let mut H = -2.0;
        while H < 2.0 {
            let f = reb_E_to_f(e, H);
            let cf = (H.cosh() - e) / (1.0 - e * H.cosh());
            assert!(
                (f.cos() - cf).abs() < 1e-11,
                "E_to_f(e={}, H={}) = {}: cos {} vs identity {}",
                e,
                H,
                f,
                f.cos(),
                cf
            );
            H += 0.23;
        }
    }
}

#[test]
fn M_to_f_is_M_to_E_then_E_to_f() {
    for &e in &[0.0, 0.3, 0.8] {
        let mut M = -4.0;
        while M < 4.0 {
            let f = reb_M_to_f(e, M);
            let f2 = reb_E_to_f(e, reb_M_to_E(e, M));
            assert_eq!(
                f.to_bits(),
                f2.to_bits(),
                "reb_M_to_f must compose exactly (e={}, M={})",
                e,
                M
            );
            // and the true anomaly must satisfy Kepler through the r(f) relation
            let E = reb_M_to_E(e, M);
            let resid = ang_diff(E - e * E.sin(), reb_mod2pi(M));
            assert!(resid < 1e-11, "residual {} for e={}, M={}", resid, e, M);
            M += 0.29;
        }
    }
}

// ===========================================================================
// 4. orbit <-> particle conversions
// ===========================================================================

#[test]
fn orbit_to_particle_satisfies_the_conic_and_vis_viva() {
    let G = 1.3;
    let primary = reb_particle {
        m: 2.0,
        x: 0.4,
        y: -0.2,
        z: 0.1,
        vx: 0.05,
        vy: -0.03,
        vz: 0.02,
        ..reb_particle::default()
    };
    let m = 0.01;
    let cases: [(f64, f64, f64, f64, f64, f64); 6] = [
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (2.5, 0.3, 0.6, 1.2, 0.4, 2.0),
        (0.7, 0.9, 2.5, 5.0, 1.1, 3.7),   // retrograde, high e
        (5.0, 0.999, 0.001, 0.2, 0.9, 0.5), // e near 1
        (-3.0, 1.4, 0.8, 2.2, 0.6, 0.4),  // hyperbolic
        (1.0, 0.2, PI / 2.0, 1.0, 2.0, 1.0), // polar
    ];
    for (i, &(a, e, inc, Omega, omega, f)) in cases.iter().enumerate() {
        let p = reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f);
        let o = reb_orbit_from_particle(G, p, primary);
        let mu = G * (m + primary.m);

        // the conic equation r = a(1-e^2)/(1+e cos f)
        let rconic = a * (1. - e * e) / (1. + e * f.cos());
        assert!(
            relerr(o.d, rconic) < 1e-12,
            "case {}: separation {} vs conic {}",
            i,
            o.d,
            rconic
        );
        // vis-viva v^2 = mu(2/r - 1/a)
        let v2 = mu * (2.0 / o.d - 1.0 / a);
        assert!(
            relerr(o.v * o.v, v2) < 1e-11,
            "case {}: v^2 {} vs vis-viva {}",
            i,
            o.v * o.v,
            v2
        );
        // specific angular momentum h^2 = mu a (1-e^2)
        let h2 = mu * a * (1. - e * e);
        assert!(
            relerr(o.h * o.h, h2) < 1e-11,
            "case {}: h^2 {} vs mu a (1-e^2) = {}",
            i,
            o.h * o.h,
            h2
        );
        // the elements themselves come back
        assert!(relerr(o.a, a) < 1e-11, "case {}: a {} vs {}", i, o.a, a);
        assert!(
            (o.e - e).abs() < 1e-11 * (1.0 + e),
            "case {}: e {} vs {}",
            i,
            o.e,
            e
        );
        assert!(
            (o.inc - inc).abs() < 1e-11,
            "case {}: inc {} vs {}",
            i,
            o.inc,
            inc
        );
    }
}

#[test]
fn particle_orbit_particle_round_trip() {
    let G = 0.9;
    let primary = reb_particle {
        m: 1.5,
        x: -0.3,
        y: 0.7,
        z: -0.11,
        vx: 0.02,
        vy: 0.04,
        vz: -0.01,
        ..reb_particle::default()
    };
    let m = 5e-4;
    let cases: [(f64, f64, f64, f64, f64, f64); 8] = [
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),          // circular, planar (degenerate angles)
        (1.0, 0.0, 0.6, 1.1, 0.0, 2.3),          // circular, inclined
        (2.0, 1e-10, 0.0, 0.0, 0.0, 1.0),        // e -> 0 planar
        (2.0, 0.4, 1e-10, 0.0, 0.9, 1.0),        // inc -> 0
        (3.0, 0.4, 2.9, 0.8, 0.9, 1.0),          // retrograde, nearly planar
        (3.0, 0.4, 2.0, 0.8, 0.9, 1.0),          // retrograde, inclined
        (4.0, 0.995, 0.5, 2.0, 1.5, 0.3),        // e near 1
        (-2.0, 2.5, 1.0, 0.5, 0.2, 0.1),         // hyperbolic
    ];
    for (i, &(a, e, inc, Omega, omega, f)) in cases.iter().enumerate() {
        let p = reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f);
        let o = reb_orbit_from_particle(G, p, primary);
        let p2 = reb_particle_from_orbit(G, primary, m, o.a, o.e, o.inc, o.Omega, o.omega, o.f);
        let scale = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt().max(1e-12);
        let vscale = (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz).sqrt().max(1e-12);
        for (name, gv, wv, s) in [
            ("x", p2.x, p.x, scale),
            ("y", p2.y, p.y, scale),
            ("z", p2.z, p.z, scale),
            ("vx", p2.vx, p.vx, vscale),
            ("vy", p2.vy, p.vy, vscale),
            ("vz", p2.vz, p.vz, vscale),
        ] {
            assert!(
                (gv - wv).abs() < 1e-10 * s,
                "case {}: round-tripped {} = {} vs original {} (scale {})",
                i,
                name,
                gv,
                wv,
                s
            );
        }
    }
}

#[test]
fn orbit_angular_momentum_direction_follows_inc_and_Omega() {
    // For the Murray & Dermott convention used by reb_particle_from_orbit the
    // unit angular momentum is (sin i sin Om, -sin i cos Om, cos i).
    let G = 1.0;
    let primary = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    for &(inc, Omega) in &[
        (0.0_f64, 0.0_f64),
        (0.4, 0.0),
        (0.4, 1.7),
        (1.2, 4.0),
        (2.6, 2.2),
    ] {
        let p = reb_particle_from_orbit(G, primary, 1e-6, 1.0, 0.25, inc, Omega, 0.6, 1.4);
        let o = reb_orbit_from_particle(G, p, primary);
        let hx = o.hvec.x / o.h;
        let hy = o.hvec.y / o.h;
        let hz = o.hvec.z / o.h;
        assert!(
            (hx - inc.sin() * Omega.sin()).abs() < 1e-12
                && (hy + inc.sin() * Omega.cos()).abs() < 1e-12
                && (hz - inc.cos()).abs() < 1e-12,
            "inc={}, Omega={}: hhat = ({}, {}, {}) vs expected ({}, {}, {})",
            inc,
            Omega,
            hx,
            hy,
            hz,
            inc.sin() * Omega.sin(),
            -inc.sin() * Omega.cos(),
            inc.cos()
        );
        // and the eccentricity vector must have length e and point at pericentre
        let ev = (o.evec.x * o.evec.x + o.evec.y * o.evec.y + o.evec.z * o.evec.z).sqrt();
        assert!(
            relerr(ev, o.e) < 1e-12,
            "|evec| = {} vs o.e = {}",
            ev,
            o.e
        );
        // e and h are perpendicular
        let dot = o.evec.x * o.hvec.x + o.evec.y * o.hvec.y + o.evec.z * o.hvec.z;
        assert!(
            dot.abs() < 1e-12 * o.h * o.e.max(1e-12),
            "evec.hvec = {} should vanish",
            dot
        );
    }
}

#[test]
fn orbit_period_and_mean_motion_are_consistent() {
    let G = 4.0 * PI * PI; // AU, yr, solar masses
    let primary = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    for &a in &[0.4_f64, 1.0, 5.2, 30.0] {
        let p = reb_particle_from_orbit(G, primary, 0.0, a, 0.1, 0.0, 0.0, 0.0, 0.0);
        let o = reb_orbit_from_particle(G, p, primary);
        // Kepler's third law with m = 0: P = a^{3/2} years
        assert!(
            relerr(o.P, a.powf(1.5)) < 1e-12,
            "a = {}: period {} vs Kepler's third law {}",
            a,
            o.P,
            a.powf(1.5)
        );
        assert!(
            relerr(o.n * o.P, 2.0 * PI) < 1e-14,
            "n*P = {} should be 2pi",
            o.n * o.P
        );
    }
}

#[test]
fn orbit_planar_orbit_has_exactly_zero_z() {
    // sin(inc) with inc = 0.0 is exactly 0, so z and vz must be exactly the
    // primary's.
    let primary = reb_particle {
        m: 1.0,
        z: 0.25,
        vz: -0.125,
        ..reb_particle::default()
    };
    let p = reb_particle_from_orbit(1.0, primary, 1e-3, 2.0, 0.4, 0.0, 1.1, 0.7, 2.2);
    assert_eq!(
        p.z.to_bits(),
        0.25f64.to_bits(),
        "inc = 0 must leave z exactly at the primary's z"
    );
    assert_eq!(
        p.vz.to_bits(),
        (-0.125f64).to_bits(),
        "inc = 0 must leave vz exactly at the primary's vz"
    );
    assert_eq!(p.m.to_bits(), 1e-3f64.to_bits());
    assert_eq!(p.ax, 0.0);
    assert_eq!(p.ay, 0.0);
    assert_eq!(p.az, 0.0);
}

#[test]
fn particle_from_orbit_error_codes() {
    let good = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    let massless = reb_particle::default();
    let cases: [(reb_particle, f64, f64, f64, i32); 6] = [
        (good, 1.0, 1.0, 0.0, 1),   // e exactly 1
        (good, 1.0, -0.1, 0.0, 2),  // e < 0
        (good, 1.0, 1.5, 0.0, 3),   // bound a with e > 1
        (good, -1.0, 0.5, 0.0, 4),  // unbound a with e < 1
        (good, -1.0, 2.0, PI, 5),   // f beyond the asymptote (e cos f < -1)
        (massless, 1.0, 0.5, 0.0, 6), // primary has no mass
    ];
    for (i, &(primary, a, e, f, want)) in cases.iter().enumerate() {
        let mut err = -99;
        let p = reb_particle_from_orbit_err(1.0, primary, 1e-3, a, e, 0.0, 0.0, 0.0, f, &mut err);
        assert_eq!(
            err, want,
            "case {} (a={}, e={}, f={}): error code {} expected {}",
            i, a, e, f, err, want
        );
        assert!(
            p.x.is_nan() && p.vx.is_nan() && p.m.is_nan(),
            "case {}: a rejected orbit must give the NaN particle",
            i
        );
    }
    // and the good case leaves err untouched at 0
    let mut err = 0;
    let p = reb_particle_from_orbit_err(1.0, good, 1e-3, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, &mut err);
    assert_eq!(err, 0, "a valid orbit must not set an error code");
    assert!(p.x.is_finite());
}

#[test]
fn orbit_from_particle_error_codes() {
    let massless = reb_particle::default();
    let p = reb_particle {
        m: 1e-3,
        x: 1.0,
        vy: 1.0,
        ..reb_particle::default()
    };
    let mut err = 0;
    let o = reb_orbit_from_particle_err(1.0, p, massless, &mut err);
    assert_eq!(err, 1, "a massless primary must be error 1");
    assert!(o.a.is_nan() && o.e.is_nan(), "error 1 must give a NaN orbit");

    // particle exactly on top of the primary -> separation below TINY
    let primary = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    let coincident = reb_particle {
        m: 1e-3,
        vy: 1.0,
        ..reb_particle::default()
    };
    let mut err2 = 0;
    let o2 = reb_orbit_from_particle_err(1.0, coincident, primary, &mut err2);
    assert_eq!(err2, 2, "a zero separation must be error 2");
    assert!(o2.d.is_nan(), "error 2 must give a NaN orbit");
}

#[test]
fn orbit_rhill_and_T_are_the_documented_expressions() {
    let G = 1.0;
    let primary = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    let m = 3e-6;
    let a = 1.0;
    let e = 0.2;
    let p = reb_particle_from_orbit(G, primary, m, a, e, 0.0, 0.0, 0.0, 0.9);
    let o = reb_orbit_from_particle(G, p, primary);
    let want = o.a * (m / (3.0 * primary.m)).cbrt();
    assert!(
        relerr(o.rhill, want) < 1e-15,
        "rhill {} vs a (m/3M)^(1/3) = {}",
        o.rhill,
        want
    );
    // With t0 = 0 the time of pericentre passage is -M/n.
    assert!(
        relerr(o.T, -o.M / o.n.abs()) < 1e-15,
        "T {} vs -M/|n| = {}",
        o.T,
        -o.M / o.n.abs()
    );
    // pericentre distance
    assert!(
        relerr(o.a * (1.0 - o.e), a * (1.0 - e)) < 1e-11,
        "pericentre {} vs {}",
        o.a * (1.0 - o.e),
        a * (1.0 - e)
    );
}

// ===========================================================================
// 5. Pal (2009) coordinates
// ===========================================================================

#[test]
fn solve_kepler_pal_satisfies_its_defining_equations_in_both_branches() {
    // Pal's (p,q) solve
    //   q cos p + p sin p = k cos L + h sin L
    //  -q sin p + p cos p = k sin L - h cos L
    // (equivalent to Kepler's equation with p = e sin E, q = e cos E).
    // e^2 < 0.09 takes the Newton branch, e^2 >= 0.09 the M_to_E branch.
    // Note on tolerances: the low-e branch is a hand-rolled Newton iteration
    // whose inverse-Jacobian entries are applied in the layout the C uses
    // (fd00*f0 + fd10*f1 for q, fd01*f0 + fd11*f1 for p). That is not the
    // exact Newton step, so the iteration converges only linearly and stops
    // after 50 rounds; the residual bottoms out near 1e-8 rather than the
    // 1e-15 the loop asks for. The high-e branch goes through reb_M_to_E and
    // reaches full double precision. Both bounds below are measured, and the
    // Rust reproduces the C expression for expression.
    let low_e: [(f64, f64); 3] = [(0.0, 0.0), (0.02, -0.05), (0.15, 0.2)];
    let high_e: [(f64, f64); 3] = [(0.3, 0.4), (-0.5, 0.6), (0.0, 0.85)];
    for (cases, tol, e2tol, label) in [
        (low_e, 1e-7, 1e-7, "low-e Newton branch"),
        (high_e, 1e-13, 1e-13, "high-e reb_M_to_E branch"),
    ] {
        for &(h, k) in cases.iter() {
            let e2 = h * h + k * k;
            assert_eq!(
                e2 < 0.09,
                label.starts_with("low"),
                "case (h={}, k={}) must exercise the {}",
                h,
                k,
                label
            );
            let mut lambda = -3.0;
            while lambda < 3.0 {
                let mut p = 0.0;
                let mut q = 0.0;
                reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
                let f0 = q * p.cos() + p * p.sin() - (k * lambda.cos() + h * lambda.sin());
                let f1 = -q * p.sin() + p * p.cos() - (k * lambda.sin() - h * lambda.cos());
                assert!(
                    f0.abs() < tol && f1.abs() < tol,
                    "{}: solve_kepler_pal(h={}, k={}, lambda={}) -> (p={}, q={}) residuals {} {}",
                    label,
                    h,
                    k,
                    lambda,
                    p,
                    q,
                    f0,
                    f1
                );
                // p = e sin E and q = e cos E, so p^2 + q^2 = e^2
                assert!(
                    (p * p + q * q - e2).abs() < e2tol,
                    "{}: p^2+q^2 = {} should be e^2 = {}",
                    label,
                    p * p + q * q,
                    e2
                );
                lambda += 0.41;
            }
        }
    }
}

#[test]
fn pal_round_trip_through_a_particle() {
    let G = 1.1;
    let primary = reb_particle {
        m: 1.0,
        x: 0.2,
        y: -0.4,
        z: 0.05,
        vx: -0.01,
        vy: 0.03,
        vz: 0.02,
        ..reb_particle::default()
    };
    let m = 1e-3;
    // (a, lambda, k, h, ix, iy); the last two rows sit in the high-e branch.
    let cases: [(f64, f64, f64, f64, f64, f64); 6] = [
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.7, 0.9, 0.12, -0.05, 0.3, -0.2),
        (2.4, -1.3, -0.08, 0.02, 0.0, 0.0),
        (0.8, 2.7, 0.05, 0.05, 0.9, 0.4),
        (3.1, 1.4, 0.4, 0.3, 0.25, 0.1),
        (1.2, -2.2, -0.5, 0.35, 0.0, 0.6),
    ];
    for (i, &(a, lambda, k, h, ix, iy)) in cases.iter().enumerate() {
        let p = reb_particle_from_pal(G, primary, m, a, lambda, k, h, ix, iy);
        let (mut a2, mut l2, mut k2, mut h2, mut ix2, mut iy2) = (0., 0., 0., 0., 0., 0.);
        reb_tools_particle_to_pal(
            G, p, primary, &mut a2, &mut l2, &mut k2, &mut h2, &mut ix2, &mut iy2,
        );
        assert!(relerr(a2, a) < 1e-11, "case {}: a {} vs {}", i, a2, a);
        assert!(
            ang_diff(l2, lambda) < 1e-10,
            "case {}: lambda {} vs {}",
            i,
            l2,
            lambda
        );
        assert!((k2 - k).abs() < 1e-11, "case {}: k {} vs {}", i, k2, k);
        assert!((h2 - h).abs() < 1e-11, "case {}: h {} vs {}", i, h2, h);
        assert!((ix2 - ix).abs() < 1e-11, "case {}: ix {} vs {}", i, ix2, ix);
        assert!((iy2 - iy).abs() < 1e-11, "case {}: iy {} vs {}", i, iy2, iy);
    }
}

#[test]
fn pal_elements_cross_check_against_classical_elements() {
    // Independent identities:
    //   e^2 = h^2 + k^2 and ix^2 + iy^2 = 2(1 - cos inc).
    let G = 1.0;
    let primary = reb_particle {
        m: 1.0,
        ..reb_particle::default()
    };
    let m = 1e-4;
    for &(a, e, inc, Omega, omega, f) in &[
        (1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64),
        (1.5, 0.2, 0.4, 1.1, 0.6, 2.0),
        (2.0, 0.45, 1.2, 4.0, 2.2, 5.0),
        (0.6, 0.7, 2.4, 0.3, 1.0, 1.0),
    ] {
        let p = reb_particle_from_orbit(G, primary, m, a, e, inc, Omega, omega, f);
        let (mut pa, mut pl, mut pk, mut ph, mut pix, mut piy) = (0., 0., 0., 0., 0., 0.);
        reb_tools_particle_to_pal(
            G, p, primary, &mut pa, &mut pl, &mut pk, &mut ph, &mut pix, &mut piy,
        );
        let o = reb_orbit_from_particle(G, p, primary);
        assert!(
            relerr(pa, o.a) < 1e-11,
            "Pal a {} vs classical a {}",
            pa,
            o.a
        );
        let e_pal = (ph * ph + pk * pk).sqrt();
        assert!(
            (e_pal - o.e).abs() < 1e-11,
            "Pal sqrt(h^2+k^2) = {} vs classical e = {}",
            e_pal,
            o.e
        );
        let inc_pal = (1.0 - (pix * pix + piy * piy) / 2.0).acos();
        assert!(
            (inc_pal - o.inc).abs() < 1e-10,
            "Pal acos(1-(ix^2+iy^2)/2) = {} vs classical inc = {}",
            inc_pal,
            o.inc
        );
        // reb_orbit carries the same Pal quantities, computed by the same
        // expressions, so they must agree bit for bit.
        assert_eq!(o.pal_h.to_bits(), ph.to_bits(), "orbit pal_h vs to_pal h");
        assert_eq!(o.pal_k.to_bits(), pk.to_bits(), "orbit pal_k vs to_pal k");
        assert_eq!(
            o.pal_ix.to_bits(),
            pix.to_bits(),
            "orbit pal_ix vs to_pal ix"
        );
        assert_eq!(
            o.pal_iy.to_bits(),
            piy.to_bits(),
            "orbit pal_iy vs to_pal iy"
        );
    }
}

#[test]
fn pal_zero_inclination_vector_gives_an_exactly_planar_orbit() {
    let primary = reb_particle {
        m: 2.0,
        z: 0.75,
        vz: -0.5,
        ..reb_particle::default()
    };
    let p = reb_particle_from_pal(1.0, primary, 1e-3, 1.6, 0.8, 0.1, -0.2, 0.0, 0.0);
    // W = eta*ix - xi*iy vanishes exactly when ix = iy = 0
    assert_eq!(
        p.z.to_bits(),
        0.75f64.to_bits(),
        "ix = iy = 0 must leave z exactly at the primary's z"
    );
    assert_eq!(
        p.vz.to_bits(),
        (-0.5f64).to_bits(),
        "ix = iy = 0 must leave vz exactly at the primary's vz"
    );
    assert_eq!(p.m.to_bits(), 1e-3f64.to_bits());
}

// ===========================================================================
// 6. energy and angular momentum
// ===========================================================================

#[test]
fn two_body_energy_is_minus_G_m1_m2_over_2a() {
    for &(G, m1, m2, a, e) in &[
        (1.0_f64, 1.0_f64, 1e-3_f64, 1.0_f64, 0.0_f64),
        (1.0, 1.0, 1e-3, 1.0, 0.7),
        (39.4784176043574, 1.0, 9.54e-4, 5.2, 0.048),
        (1.0, 3.0, 2.0, 0.5, 0.9),
    ] {
        let mut r = reb_simulation_create();
        r.G = G;
        r.save_messages = 1;
        reb_simulation_add_fmt(&mut r, "m", &[d(m1)]);
        reb_simulation_add_fmt(&mut r, "m a e f", &[d(m2), d(a), d(e), d(1.234)]);
        reb_simulation_move_to_com(&mut r);
        let want = -G * m1 * m2 / (2.0 * a);
        let got = reb_simulation_energy(&r);
        assert!(
            relerr(got, want) < 1e-12,
            "two-body energy {} vs -G m1 m2/(2a) = {} (G={}, a={}, e={})",
            got,
            want,
            G,
            a,
            e
        );
        // and it agrees with an independently written double sum
        assert!(
            relerr(got, energy_by_hand(&r)) < 1e-14,
            "reb_simulation_energy {} vs hand-rolled sum {}",
            got,
            energy_by_hand(&r)
        );
    }
}

#[test]
fn two_body_angular_momentum_is_the_reduced_mass_times_h() {
    let G = 1.0;
    let (m1, m2) = (1.0_f64, 2e-3_f64);
    let (a, e, inc, Omega) = (2.0_f64, 0.35_f64, 0.9_f64, 1.7_f64);
    let mut r = reb_simulation_create();
    r.G = G;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(m1)]);
    reb_simulation_add_fmt(
        &mut r,
        "m a e inc Omega omega f",
        &[d(m2), d(a), d(e), d(inc), d(Omega), d(0.4), d(2.1)],
    );
    reb_simulation_move_to_com(&mut r);
    let L = reb_simulation_angular_momentum(&r);
    let mag = (L.x * L.x + L.y * L.y + L.z * L.z).sqrt();
    let mu_red = m1 * m2 / (m1 + m2);
    let want = mu_red * (G * (m1 + m2) * a * (1. - e * e)).sqrt();
    assert!(
        relerr(mag, want) < 1e-11,
        "|L| = {} vs mu_red*sqrt(G M a (1-e^2)) = {}",
        mag,
        want
    );
    // direction: (sin i sin Om, -sin i cos Om, cos i)
    assert!(
        (L.x / mag - inc.sin() * Omega.sin()).abs() < 1e-11
            && (L.y / mag + inc.sin() * Omega.cos()).abs() < 1e-11
            && (L.z / mag - inc.cos()).abs() < 1e-11,
        "Lhat = ({}, {}, {}) vs expected ({}, {}, {})",
        L.x / mag,
        L.y / mag,
        L.z / mag,
        inc.sin() * Omega.sin(),
        -inc.sin() * Omega.cos(),
        inc.cos()
    );
}

#[test]
fn angular_momentum_of_a_single_particle_is_exact() {
    let mut r = reb_simulation_create();
    r.save_messages = 1;
    reb_simulation_add_fmt(
        &mut r,
        "m x y z vx vy vz",
        &[d(2.0), d(1.0), d(0.0), d(0.0), d(0.0), d(1.0), d(0.0)],
    );
    let L = reb_simulation_angular_momentum(&r);
    assert_eq!(L.x.to_bits(), 0.0f64.to_bits(), "Lx must be exactly 0");
    assert_eq!(L.y.to_bits(), 0.0f64.to_bits(), "Ly must be exactly 0");
    assert_eq!(L.z, 2.0, "Lz = m (x vy - y vx) = 2*1*1");
}

#[test]
fn energy_and_angular_momentum_degenerate_particle_counts() {
    // N = 0
    let mut r = reb_simulation_create();
    assert_eq!(
        reb_simulation_energy(&r),
        0.0,
        "an empty simulation has zero energy"
    );
    let L = reb_simulation_angular_momentum(&r);
    assert_eq!((L.x, L.y, L.z), (0.0, 0.0, 0.0));
    // energy_offset passes straight through when there are no particles
    r.energy_offset = 3.25;
    assert_eq!(
        reb_simulation_energy(&r).to_bits(),
        3.25f64.to_bits(),
        "energy_offset must be added verbatim"
    );

    // N = 1: kinetic energy only, no self interaction
    let mut r1 = reb_simulation_create();
    r1.save_messages = 1;
    reb_simulation_add_fmt(
        &mut r1,
        "m x vx vy",
        &[d(4.0), d(1.0), d(0.5), d(-1.5)],
    );
    let want = 0.5 * 4.0 * (0.25 + 2.25);
    assert!(
        relerr(reb_simulation_energy(&r1), want) < 1e-15,
        "single-particle energy {} vs 0.5 m v^2 = {}",
        reb_simulation_energy(&r1),
        want
    );
    let L1 = reb_simulation_angular_momentum(&r1);
    assert_eq!(L1.z, 4.0 * (1.0 * -1.5 - 0.0 * 0.5), "Lz for one particle");
}

#[test]
fn energy_respects_N_active_and_testparticle_type() {
    // Star + planet + a small third body that is flagged as a test particle.
    let build = || {
        let mut r = reb_simulation_create();
        r.G = 1.0;
        r.save_messages = 1;
        reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
        reb_simulation_add_fmt(&mut r, "m a e f", &[d(1e-3), d(1.0), d(0.1), d(0.5)]);
        reb_simulation_add_fmt(&mut r, "m a e f", &[d(1e-6), d(2.0), d(0.2), d(1.5)]);
        r
    };

    // Reference: the same two massive bodies with no third particle at all.
    let mut two = reb_simulation_create();
    two.G = 1.0;
    two.save_messages = 1;
    reb_simulation_add_fmt(&mut two, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut two, "m a e f", &[d(1e-3), d(1.0), d(0.1), d(0.5)]);
    let e_two = reb_simulation_energy(&two);

    // testparticle_type = 0: the third body is invisible to the energy.
    let mut r0 = build();
    r0.N_active = 2;
    r0.testparticle_type = 0;
    assert_eq!(
        reb_simulation_energy(&r0).to_bits(),
        e_two.to_bits(),
        "a type-0 test particle must not enter the energy at all"
    );

    // testparticle_type = 1: its kinetic energy and its potential against the
    // N_active massive bodies do enter.
    let mut r1 = build();
    r1.N_active = 2;
    r1.testparticle_type = 1;
    let e1 = reb_simulation_energy(&r1);
    let p = r1.particles[2];
    let mut extra = 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    for i in 0..2 {
        let q = r1.particles[i];
        extra -= r1.G * p.m * q.m
            / ((q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y) + (q.z - p.z) * (q.z - p.z))
                .sqrt();
    }
    assert!(
        relerr(e1, e_two + extra) < 1e-13,
        "type-1 test particle energy {} vs two-body {} plus its own {}",
        e1,
        e_two,
        extra
    );
    assert!(
        (e1 - e_two).abs() > 1e-12,
        "the type-1 contribution must be measurable, got {}",
        e1 - e_two
    );
}

#[test]
fn energy_and_angular_momentum_are_conserved_by_ias15() {
    let mut r = three_body();
    let e0 = reb_simulation_energy(&r);
    let l0 = reb_simulation_angular_momentum(&r);
    let l0m = (l0.x * l0.x + l0.y * l0.y + l0.z * l0.z).sqrt();
    r.dt = 0.01;
    reb_simulation_integrate(&mut r, 200.0);
    let e1 = reb_simulation_energy(&r);
    let l1 = reb_simulation_angular_momentum(&r);
    let l1m = (l1.x * l1.x + l1.y * l1.y + l1.z * l1.z).sqrt();
    assert!(
        relerr(e1, e0) < 1e-14,
        "IAS15 relative energy drift {} over ~30 orbits",
        relerr(e1, e0)
    );
    assert!(
        relerr(l1m, l0m) < 1e-14,
        "IAS15 relative |L| drift {}",
        relerr(l1m, l0m)
    );
    assert!(
        (l1.z - l0.z).abs() < 1e-14 * l0m,
        "Lz drift {} (|L| = {})",
        l1.z - l0.z,
        l0m
    );
}

#[test]
fn energy_is_conserved_by_whfast() {
    let mut r = three_body();
    reb_simulation_set_integrator(&mut r, "whfast");
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.corrector = 17;
        wh.safe_mode = 1;
    }
    r.dt = 1.0 / 200.0;
    let e0 = reb_simulation_energy(&r);
    reb_simulation_integrate(&mut r, 100.0);
    let e1 = reb_simulation_energy(&r);
    assert!(
        relerr(e1, e0) < 1e-12,
        "WHFast relative energy drift {} over ~16 orbits at dt = P/200",
        relerr(e1, e0)
    );
}

#[test]
fn ias15_and_whfast_agree_on_the_same_two_body_problem() {
    // A single Kepler orbit is handled exactly by WHFast's Kepler solver and
    // to machine precision by IAS15, so the two must land on the same state.
    let build = || {
        let mut r = reb_simulation_create();
        r.G = 1.0;
        r.save_messages = 1;
        reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
        reb_simulation_add_fmt(&mut r, "m a e f", &[d(1e-3), d(1.0), d(0.3), d(0.4)]);
        reb_simulation_move_to_com(&mut r);
        r
    };
    let tmax = 20.0 * 2.0 * PI;

    let mut a = build();
    reb_simulation_integrate(&mut a, tmax);

    let mut w = build();
    reb_simulation_set_integrator(&mut w, "whfast");
    w.dt = 2.0 * PI / 1000.0;
    reb_simulation_integrate(&mut w, tmax);

    let sep = reb_particle_distance(&a.particles[1], &w.particles[1]);
    assert!(
        sep < 1e-9,
        "IAS15 and WHFast disagree by {} in position after 20 orbits",
        sep
    );
    let oa = reb_orbit_from_particle(a.G, a.particles[1], a.particles[0]);
    let ow = reb_orbit_from_particle(w.G, w.particles[1], w.particles[0]);
    assert!(
        relerr(ow.a, oa.a) < 1e-11,
        "semi-major axis: IAS15 {} vs WHFast {}",
        oa.a,
        ow.a
    );
    assert!(
        (ow.e - oa.e).abs() < 1e-11,
        "eccentricity: IAS15 {} vs WHFast {}",
        oa.e,
        ow.e
    );
    assert!(
        ang_diff(ow.f, oa.f) < 1e-9,
        "true anomaly: IAS15 {} vs WHFast {}",
        oa.f,
        ow.f
    );
}

// ===========================================================================
// 7. the centre-of-mass family
// ===========================================================================

#[test]
fn com_of_pair_is_the_mass_weighted_average() {
    let p1 = reb_particle {
        m: 1.0,
        x: 0.0,
        vy: 4.0,
        ax: 8.0,
        ..reb_particle::default()
    };
    let p2 = reb_particle {
        m: 3.0,
        x: 4.0,
        vy: 0.0,
        ax: 0.0,
        ..reb_particle::default()
    };
    let c = reb_particle_com_of_pair(p1, p2);
    assert_eq!(c.m, 4.0, "total mass");
    assert_eq!(c.x, 3.0, "(0*1 + 4*3)/4 = 3 exactly");
    assert_eq!(c.vy, 1.0, "(4*1 + 0*3)/4 = 1 exactly");
    assert_eq!(c.ax, 2.0, "accelerations are averaged too");
    // order does not matter for this pair (all the arithmetic is exact here)
    let c2 = reb_particle_com_of_pair(p2, p1);
    assert_eq!(pbits(&c), pbits(&c2), "com_of_pair must be symmetric here");
}

#[test]
fn com_of_two_massless_particles_stays_at_the_origin() {
    // The C guards the division with `if (p1.m > 0.)`; with zero total mass
    // the un-normalised (mass weighted) sums are returned, which are all zero.
    let p1 = reb_particle {
        x: 5.0,
        vy: -2.0,
        ..reb_particle::default()
    };
    let p2 = reb_particle {
        x: 7.0,
        vy: 9.0,
        ..reb_particle::default()
    };
    let c = reb_particle_com_of_pair(p1, p2);
    assert_eq!(c.m, 0.0, "massless pair has zero mass");
    assert_eq!(c.x, 0.0, "x = 5*0 + 7*0 = 0, and no division happens");
    assert_eq!(c.vy, 0.0);
}

#[test]
fn simulation_com_and_com_range_agree_with_a_direct_sum() {
    let mut r = reb_simulation_create();
    r.save_messages = 1;
    let data = [
        (1.0, 1.0, 2.0, 3.0, 0.5, -0.25, 0.125),
        (2.0, -1.0, 0.0, 1.0, 0.25, 0.5, -0.5),
        (0.5, 4.0, -2.0, 0.0, -1.0, 0.0, 0.25),
        (4.0, 0.0, 0.5, -1.5, 0.125, 0.25, 0.0),
    ];
    for row in data.iter() {
        reb_simulation_add_fmt(
            &mut r,
            "m x y z vx vy vz",
            &[
                d(row.0),
                d(row.1),
                d(row.2),
                d(row.3),
                d(row.4),
                d(row.5),
                d(row.6),
            ],
        );
    }
    let com = reb_simulation_com(&r);
    let mtot: f64 = data.iter().map(|q| q.0).sum();
    let mut sx = 0.0;
    let mut svz = 0.0;
    for q in data.iter() {
        sx += q.0 * q.1;
        svz += q.0 * q.6;
    }
    assert_eq!(com.m, mtot, "total mass");
    assert!(
        relerr(com.x, sx / mtot) < 1e-15,
        "com.x {} vs sum(m x)/sum(m) {}",
        com.x,
        sx / mtot
    );
    assert!(
        relerr(com.vz, svz / mtot) < 1e-15,
        "com.vz {} vs sum(m vz)/sum(m) {}",
        com.vz,
        svz / mtot
    );

    // com_range over the whole array is reb_simulation_com
    assert_eq!(
        pbits(&reb_simulation_com_range(&r, 0, r.N)),
        pbits(&com),
        "com_range(0,N) must be com"
    );
    // an empty range gives the zero particle
    assert_eq!(
        pbits(&reb_simulation_com_range(&r, 2, 2)),
        pbits(&reb_particle::default()),
        "an empty com range is the default particle"
    );
    // a single-element range gives that particle's mass and position back
    let one = reb_simulation_com_range(&r, 1, 2);
    assert_eq!(one.m, data[1].0);
    assert!(relerr(one.x, data[1].1) < 1e-15);

    // splitting the range: com(0..4) is com_of_pair(com(0..2), com(2..4))
    let left = reb_simulation_com_range(&r, 0, 2);
    let right = reb_simulation_com_range(&r, 2, 4);
    let joined = reb_particle_com_of_pair(left, right);
    assert!(
        relerr(joined.x, com.x) < 1e-14 && relerr(joined.m, com.m) < 1e-15,
        "com must be associative over a split: {} vs {}",
        joined.x,
        com.x
    );
}

#[test]
fn jacobi_com_is_the_com_of_the_preceding_particles() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(0.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(2.0), d(3.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(-4.0)]);
    // index 0: nothing precedes it
    assert_eq!(
        pbits(&reb_simulation_jacobi_com(&r, 0)),
        pbits(&reb_particle::default()),
        "jacobi_com of particle 0 is the zero particle"
    );
    let j1 = reb_simulation_jacobi_com(&r, 1);
    assert_eq!(j1.m, 1.0);
    assert_eq!(j1.x, 0.0);
    let j2 = reb_simulation_jacobi_com(&r, 2);
    assert_eq!(j2.m, 3.0, "mass of particles 0 and 1");
    assert_eq!(j2.x, 2.0, "(0*1 + 3*2)/3 = 2 exactly");
    // out of range gives the NaN particle
    let bad = reb_simulation_jacobi_com(&r, 3);
    assert!(
        bad.x.is_nan() && bad.m.is_nan(),
        "an out-of-range index must give the NaN particle"
    );
}

#[test]
fn move_to_com_zeroes_the_centre_of_mass_and_conserves_the_relative_state() {
    let mut r = three_body();
    // three_body() already moved to the COM; check it took.
    let com = reb_simulation_com(&r);
    let scale = r.particles.iter().map(|p| p.x.abs()).fold(0.0, f64::max);
    assert!(
        com.x.abs() < 1e-15 * scale
            && com.y.abs() < 1e-15 * scale
            && com.z.abs() < 1e-15 * scale,
        "COM position ({}, {}, {}) should vanish",
        com.x,
        com.y,
        com.z
    );
    assert!(
        com.vx.abs() < 1e-15 && com.vy.abs() < 1e-15 && com.vz.abs() < 1e-15,
        "COM velocity ({}, {}, {}) should vanish",
        com.vx,
        com.vy,
        com.vz
    );
    // total linear momentum vanishes too
    let px: f64 = (0..r.N).map(|i| r.particles[i].m * r.particles[i].vx).sum();
    assert!(px.abs() < 1e-17, "total x-momentum {} should vanish", px);

    // moving to the COM a second time is a no-op at the 1e-16 level, and the
    // separations are untouched
    let before: Vec<reb_particle> = r.particles.clone();
    let sep_before = reb_particle_distance(&before[1], &before[2]);
    reb_simulation_move_to_com(&mut r);
    for i in 0..r.N {
        assert!(
            (r.particles[i].x - before[i].x).abs() < 1e-15 * scale,
            "second move_to_com shifted particle {} by {}",
            i,
            r.particles[i].x - before[i].x
        );
    }
    let sep_after = reb_particle_distance(&r.particles[1], &r.particles[2]);
    assert!(
        relerr(sep_after, sep_before) < 1e-14,
        "move_to_com must not change separations: {} vs {}",
        sep_after,
        sep_before
    );
    // and the energy is a COM-frame invariant only up to the bulk kinetic
    // term, which is already zero here
    assert!(
        relerr(reb_simulation_energy(&r), energy_by_hand(&r)) < 1e-13,
        "energy after move_to_com"
    );
}

#[test]
fn move_to_com_handles_second_order_variational_particles() {
    // This replaces an earlier test that asserted move_to_com bailed out
    // when a 2nd-order variational configuration was present. That was
    // asserting a PORT DEFECT: tools.c has a full 2nd-order COM-shift
    // block, and the Rust had been returning early instead -- which also
    // skipped the ordinary particle shift and the boundary check. The
    // block is now translated, so the correct expectation is the C's:
    // every particle IS shifted and the centre of mass lands on the
    // origin. Verified bit-for-bit against the C build by
    // porttest/movetocom_var_c.c.
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(1.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(-1.0)]);
    // shift the whole system so a COM move is visible
    r.particles[0].x += 10.0;
    r.particles[1].x += 10.0;
    let before: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();

    // One 1st-order set (index 0) and a 2nd-order set differentiating
    // twice with respect to it. Note we deliberately use a SINGLE
    // 1st-order configuration: with two of them the second gets
    // index > 0, and tools.c then reads particles[i+index] past the
    // end of the particle array (undefined behaviour upstream, which
    // this port reports rather than imitates).
    reb_simulation_add_variation_1st_order(&mut r, -1);
    let a = r.var_config[0].index;
    reb_simulation_add_variation_2nd_order(&mut r, -1, a, a);
    r.messages.clear();

    reb_simulation_move_to_com(&mut r);

    let after: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_ne!(
        before, after,
        "move_to_com must shift the particles even with a 2nd-order          variational configuration present (the C does)"
    );
    let com = reb_simulation_com(&r);
    assert!(
        com.x.abs() < 1e-15 && com.y.abs() < 1e-15 && com.z.abs() < 1e-15,
        "centre of mass must be at the origin after move_to_com, got ({}, {}, {})",
        com.x,
        com.y,
        com.z
    );
    assert!(
        com.vx.abs() < 1e-15 && com.vy.abs() < 1e-15 && com.vz.abs() < 1e-15,
        "centre-of-mass velocity must be zero after move_to_com"
    );
    assert!(
        !r.messages
            .iter()
            .any(|(t, _)| *t == REB_MESSAGE_TYPE::ERROR),
        "no error should be reported now that the 2nd-order block is          translated, got {:?}",
        r.messages
    );
}

#[test]
fn move_to_hel_puts_particle_zero_exactly_at_the_origin() {
    let mut r = three_body();
    let before: Vec<reb_particle> = r.particles.clone();
    reb_simulation_move_to_hel(&mut r);
    let p0 = r.particles[0];
    assert_eq!(
        pbits(&p0)[0..6],
        [0u64; 6],
        "after move_to_hel particle 0 must be exactly at rest at the origin"
    );
    for i in 1..r.N {
        assert_eq!(
            r.particles[i].x.to_bits(),
            (before[i].x - before[0].x).to_bits(),
            "particle {} x must be exactly the original difference",
            i
        );
        assert_eq!(
            r.particles[i].vz.to_bits(),
            (before[i].vz - before[0].vz).to_bits(),
            "particle {} vz must be exactly the original difference",
            i
        );
        // masses are untouched
        assert_eq!(r.particles[i].m.to_bits(), before[i].m.to_bits());
    }
}

// ===========================================================================
// 8. imul / iadd / isub
// ===========================================================================

#[test]
fn imul_by_powers_of_two_round_trips_bit_exactly() {
    let mut r = three_body();
    let before: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    reb_simulation_imul(&mut r, 2.0, 0.5);
    // the scaling really happened
    for i in 0..r.N {
        assert_eq!(
            r.particles[i].x.to_bits(),
            (f64::from_bits(before[i][0]) * 2.0).to_bits(),
            "particle {} x must be doubled",
            i
        );
        assert_eq!(
            r.particles[i].vx.to_bits(),
            (f64::from_bits(before[i][3]) * 0.5).to_bits(),
            "particle {} vx must be halved",
            i
        );
        assert_eq!(
            r.particles[i].m.to_bits(),
            before[i][9],
            "imul must not touch the masses"
        );
    }
    reb_simulation_imul(&mut r, 0.5, 2.0);
    let after: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_eq!(
        before, after,
        "scaling by 2 then 1/2 must be bit-exactly the identity"
    );
}

#[test]
fn isub_of_a_copy_zeroes_positions_and_velocities_exactly() {
    let a = three_body();
    let mut b = three_body();
    assert_eq!(reb_simulation_isub(&mut b, &a), 0, "isub must report success");
    for i in 0..b.N {
        let p = b.particles[i];
        assert_eq!(
            (p.x, p.y, p.z, p.vx, p.vy, p.vz),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "particle {} must be exactly zeroed by subtracting itself",
            i
        );
        assert_eq!(
            p.m.to_bits(),
            a.particles[i].m.to_bits(),
            "isub must not touch the masses"
        );
    }
}

#[test]
fn iadd_then_isub_with_exact_binary_values_is_the_identity() {
    let mut r = reb_simulation_create();
    let mut r2 = reb_simulation_create();
    r.save_messages = 1;
    r2.save_messages = 1;
    for k in 0..3 {
        let k = k as f64;
        reb_simulation_add_fmt(
            &mut r,
            "m x y z vx vy vz",
            &[d(1.0), d(1.0 + k), d(2.0), d(-0.5), d(0.25), d(0.5), d(-0.125)],
        );
        reb_simulation_add_fmt(
            &mut r2,
            "m x y z vx vy vz",
            &[d(1.0), d(0.5), d(0.25), d(0.125), d(0.0625), d(0.5), d(0.25)],
        );
    }
    let before: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_eq!(reb_simulation_iadd(&mut r, &r2), 0);
    assert_eq!(
        r.particles[0].x, 1.5,
        "1.0 + 0.5 must be exactly 1.5"
    );
    assert_eq!(reb_simulation_isub(&mut r, &r2), 0);
    let after: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_eq!(
        before, after,
        "adding then subtracting exactly representable values is the identity"
    );
}

#[test]
fn iadd_and_isub_reject_mismatched_particle_counts() {
    let mut r = three_body();
    let mut small = reb_simulation_create();
    small.save_messages = 1;
    reb_simulation_add_fmt(&mut small, "m x", &[d(1.0), d(1.0)]);
    let before: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_eq!(
        reb_simulation_iadd(&mut r, &small),
        -1,
        "iadd must fail on an N mismatch"
    );
    assert_eq!(
        reb_simulation_isub(&mut r, &small),
        -1,
        "isub must fail on an N mismatch"
    );
    let after: Vec<[u64; 11]> = r.particles.iter().map(pbits).collect();
    assert_eq!(before, after, "a rejected iadd/isub must change nothing");
}

// ===========================================================================
// 9. spherical <-> cartesian
// ===========================================================================

#[test]
fn spherical_to_xyz_and_back() {
    // theta = 0 is exact: sin(0) = 0 and cos(0) = 1.
    let v = reb_tools_spherical_to_xyz(2.0, 0.0, 1.234);
    assert_eq!(v.x, 0.0, "theta = 0 gives x exactly 0");
    assert_eq!(v.y, 0.0, "theta = 0 gives y exactly 0");
    assert_eq!(v.z, 2.0, "theta = 0 gives z exactly the magnitude");
    let (mut mag, mut th, mut ph) = (0., 0., 0.);
    reb_tools_xyz_to_spherical(v, &mut mag, &mut th, &mut ph);
    assert_eq!(mag, 2.0, "magnitude of (0,0,2)");
    assert_eq!(th, 0.0, "acos2(2,2,1) short-circuits to 0");
    assert_eq!(ph, 0.0, "atan2(0,0) = 0");

    // general round trip
    for &m in &[0.5_f64, 3.0, 100.0] {
        let mut theta = 0.1;
        while theta < PI {
            let mut phi = -3.0;
            while phi < 3.0 {
                let xyz = reb_tools_spherical_to_xyz(m, theta, phi);
                let (mut m2, mut t2, mut p2) = (0., 0., 0.);
                reb_tools_xyz_to_spherical(xyz, &mut m2, &mut t2, &mut p2);
                assert!(
                    relerr(m2, m) < 1e-14,
                    "magnitude {} vs {} (theta={}, phi={})",
                    m2,
                    m,
                    theta,
                    phi
                );
                assert!(
                    (t2 - theta).abs() < 1e-12,
                    "theta {} vs {} (phi={})",
                    t2,
                    theta,
                    phi
                );
                assert!(
                    ang_diff(p2, phi) < 1e-12,
                    "phi {} vs {} (theta={})",
                    p2,
                    phi,
                    theta
                );
                // and the cartesian vector really has the requested length
                let len = (xyz.x * xyz.x + xyz.y * xyz.y + xyz.z * xyz.z).sqrt();
                assert!(relerr(len, m) < 1e-14, "|xyz| = {} vs {}", len, m);
                phi += 0.7;
            }
            theta += 0.4;
        }
    }
}

// ===========================================================================
// 10. reb_simulation_add_fmt
// ===========================================================================

#[test]
fn add_fmt_cartesian_stores_values_verbatim() {
    let mut r = reb_simulation_create();
    r.save_messages = 1;
    reb_simulation_add_fmt(
        &mut r,
        "m r x y z vx vy vz",
        &[
            d(2.0),
            d(0.125),
            d(3.0),
            d(-4.0),
            d(5.0),
            d(0.25),
            d(-0.5),
            d(0.0625),
        ],
    );
    assert_eq!(r.N, 1, "one particle added");
    let p = r.particles[0];
    assert_eq!(
        (p.m, p.r, p.x, p.y, p.z, p.vx, p.vy, p.vz),
        (2.0, 0.125, 3.0, -4.0, 5.0, 0.25, -0.5, 0.0625),
        "add_fmt must store cartesian varargs verbatim"
    );
    assert_eq!((p.ax, p.ay, p.az), (0.0, 0.0, 0.0));
    assert_eq!(p.name, None);

    // unmentioned cartesian coordinates default to zero, not NaN
    let mut r2 = reb_simulation_create();
    r2.save_messages = 1;
    reb_simulation_add_fmt(&mut r2, "m x", &[d(1.0), d(7.0)]);
    let q = r2.particles[0];
    assert_eq!(
        (q.x, q.y, q.z, q.vx, q.vy, q.vz),
        (7.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        "missing coordinates must be zero"
    );

    // an empty format string gives the zero particle
    let mut r3 = reb_simulation_create();
    r3.save_messages = 1;
    reb_simulation_add_fmt(&mut r3, "", &[]);
    assert_eq!(r3.N, 1);
    assert_eq!(
        pbits(&r3.particles[0]),
        pbits(&reb_particle::default()),
        "an empty format gives the default particle"
    );

    // separators: spaces, tabs, commas and semicolons are all delimiters
    let mut r4 = reb_simulation_create();
    r4.save_messages = 1;
    reb_simulation_add_fmt(&mut r4, "m,x;y\tz", &[d(1.0), d(2.0), d(3.0), d(4.0)]);
    let s = r4.particles[0];
    assert_eq!(
        (s.m, s.x, s.y, s.z),
        (1.0, 2.0, 3.0, 4.0),
        "',', ';' and tab must delimit format tokens"
    );
}

#[test]
fn add_fmt_orbital_elements_match_particle_from_orbit_bit_for_bit() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
    // with no `primary` token the primary is the COM of what is already there
    let primary = reb_simulation_com(&r);
    reb_simulation_add_fmt(
        &mut r,
        "m a e inc Omega omega f",
        &[d(1e-3), d(2.5), d(0.4), d(0.3), d(1.1), d(0.7), d(2.0)],
    );
    let want = reb_particle_from_orbit(1.0, primary, 1e-3, 2.5, 0.4, 0.3, 1.1, 0.7, 2.0);
    assert_eq!(
        pbits(&r.particles[1]),
        pbits(&want),
        "add_fmt must reduce exactly to reb_particle_from_orbit"
    );

    // an explicit `primary` token overrides the COM
    let other = reb_particle {
        m: 5.0,
        x: 1.0,
        ..reb_particle::default()
    };
    reb_simulation_add_fmt(
        &mut r,
        "m a e primary",
        &[d(1e-5), d(1.0), d(0.2), reb_fmt_arg::primary(other)],
    );
    let want2 = reb_particle_from_orbit(1.0, other, 1e-5, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(
        pbits(&r.particles[2]),
        pbits(&want2),
        "the `primary` token must select the primary"
    );

    // omitted elements default to zero
    let mut r2 = reb_simulation_create();
    r2.G = 1.0;
    r2.save_messages = 1;
    reb_simulation_add_fmt(&mut r2, "m", &[d(1.0)]);
    let prim2 = reb_simulation_com(&r2);
    reb_simulation_add_fmt(&mut r2, "m a", &[d(1e-4), d(3.0)]);
    let want3 = reb_particle_from_orbit(1.0, prim2, 1e-4, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(
        pbits(&r2.particles[1]),
        pbits(&want3),
        "e, inc, Omega, omega and f default to zero"
    );
}

#[test]
fn add_fmt_period_reproduces_the_requested_period() {
    let mut r = reb_simulation_create();
    r.G = 4.0 * PI * PI;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
    for &P in &[0.5_f64, 1.0, 11.86] {
        let n0 = r.N;
        // the primary of an add_fmt without a `primary` token is the COM of
        // everything already in the simulation, so it must be re-read here
        let primary = reb_simulation_com(&r);
        reb_simulation_add_fmt(&mut r, "m P e", &[d(1e-3), d(P), d(0.1)]);
        assert_eq!(r.N, n0 + 1, "particle for P = {} must be added", P);
        let o = reb_orbit_from_particle(r.G, r.particles[n0], primary);
        assert!(
            relerr(o.P, P) < 1e-12,
            "requested period {} came back as {}",
            P,
            o.P
        );
        // a = (P^2 G M / 4 pi^2)^(1/3)
        let a_want = (P * P * r.G * (primary.m + 1e-3) / (4.0 * PI * PI)).cbrt();
        assert!(
            relerr(o.a, a_want) < 1e-12,
            "semi-major axis {} vs Kepler's third law {}",
            o.a,
            a_want
        );
    }
}

#[test]
fn add_fmt_angle_aliases_reduce_to_omega_and_f() {
    let mut base = reb_simulation_create();
    base.G = 1.0;
    base.save_messages = 1;
    reb_simulation_add_fmt(&mut base, "m", &[d(1.0)]);

    // pomega with inc = 0 (cos inc > 0) is omega + Omega, and Omega defaults to 0
    let mut r1 = reb_simulation_create();
    r1.G = 1.0;
    r1.save_messages = 1;
    reb_simulation_add_fmt(&mut r1, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r1, "m a e pomega", &[d(1e-3), d(1.0), d(0.3), d(0.9)]);
    let mut r2 = reb_simulation_create();
    r2.G = 1.0;
    r2.save_messages = 1;
    reb_simulation_add_fmt(&mut r2, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r2, "m a e omega", &[d(1e-3), d(1.0), d(0.3), d(0.9)]);
    assert_eq!(
        pbits(&r1.particles[1]),
        pbits(&r2.particles[1]),
        "with Omega = 0 and inc = 0, pomega and omega must agree exactly"
    );

    // retrograde: omega = Omega - pomega
    let mut r3 = reb_simulation_create();
    r3.G = 1.0;
    r3.save_messages = 1;
    reb_simulation_add_fmt(&mut r3, "m", &[d(1.0)]);
    reb_simulation_add_fmt(
        &mut r3,
        "m a e inc Omega pomega",
        &[d(1e-3), d(1.0), d(0.3), d(3.0), d(0.5), d(0.9)],
    );
    let mut r4 = reb_simulation_create();
    r4.G = 1.0;
    r4.save_messages = 1;
    reb_simulation_add_fmt(&mut r4, "m", &[d(1.0)]);
    reb_simulation_add_fmt(
        &mut r4,
        "m a e inc Omega omega",
        &[d(1e-3), d(1.0), d(0.3), d(3.0), d(0.5), d(0.5 - 0.9)],
    );
    assert_eq!(
        pbits(&r3.particles[1]),
        pbits(&r4.particles[1]),
        "for a retrograde orbit omega must be Omega - pomega"
    );

    // theta (true longitude) with Omega = omega = 0 is just f
    let mut r5 = reb_simulation_create();
    r5.G = 1.0;
    r5.save_messages = 1;
    reb_simulation_add_fmt(&mut r5, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r5, "m a e theta", &[d(1e-3), d(1.0), d(0.3), d(1.7)]);
    let mut r6 = reb_simulation_create();
    r6.G = 1.0;
    r6.save_messages = 1;
    reb_simulation_add_fmt(&mut r6, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r6, "m a e f", &[d(1e-3), d(1.0), d(0.3), d(1.7)]);
    assert_eq!(
        pbits(&r5.particles[1]),
        pbits(&r6.particles[1]),
        "with Omega = omega = 0, theta must reduce to f"
    );

    // E (eccentric anomaly) goes through reb_E_to_f
    let mut r7 = reb_simulation_create();
    r7.G = 1.0;
    r7.save_messages = 1;
    reb_simulation_add_fmt(&mut r7, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r7, "m a e E", &[d(1e-3), d(1.0), d(0.3), d(1.2)]);
    let mut r8 = reb_simulation_create();
    r8.G = 1.0;
    r8.save_messages = 1;
    reb_simulation_add_fmt(&mut r8, "m", &[d(1.0)]);
    reb_simulation_add_fmt(
        &mut r8,
        "m a e f",
        &[d(1e-3), d(1.0), d(0.3), d(reb_E_to_f(0.3, 1.2))],
    );
    assert_eq!(
        pbits(&r7.particles[1]),
        pbits(&r8.particles[1]),
        "E must be converted with reb_E_to_f"
    );

    // M goes through reb_M_to_f, and the resulting orbit reports it back
    let mut r9 = reb_simulation_create();
    r9.G = 1.0;
    r9.save_messages = 1;
    reb_simulation_add_fmt(&mut r9, "m", &[d(1.0)]);
    let prim = reb_simulation_com(&r9);
    reb_simulation_add_fmt(&mut r9, "m a e M", &[d(1e-3), d(1.0), d(0.3), d(2.0)]);
    let o = reb_orbit_from_particle(1.0, r9.particles[1], prim);
    assert!(
        ang_diff(o.M, 2.0) < 1e-10,
        "mean anomaly came back as {} instead of 2.0",
        o.M
    );

    // T (pericentre time): with r.t == T the particle sits at pericentre
    let mut rt = reb_simulation_create();
    rt.G = 1.0;
    rt.save_messages = 1;
    rt.t = 1.75;
    reb_simulation_add_fmt(&mut rt, "m", &[d(1.0)]);
    let primt = reb_simulation_com(&rt);
    reb_simulation_add_fmt(&mut rt, "m a e T", &[d(1e-3), d(2.0), d(0.4), d(1.75)]);
    let ot = reb_orbit_from_particle(1.0, rt.particles[1], primt);
    assert!(
        relerr(ot.d, 2.0 * (1.0 - 0.4)) < 1e-11,
        "with t == T the particle must be at pericentre a(1-e) = {}, got d = {}",
        2.0 * 0.6,
        ot.d
    );
}

#[test]
fn add_fmt_pal_coordinates() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
    let primary = reb_simulation_com(&r);
    reb_simulation_add_fmt(
        &mut r,
        "m a l k h ix iy",
        &[d(1e-3), d(1.6), d(0.9), d(0.12), d(-0.05), d(0.3), d(-0.2)],
    );
    let want = reb_particle_from_pal(1.0, primary, 1e-3, 1.6, 0.9, 0.12, -0.05, 0.3, -0.2);
    assert_eq!(
        pbits(&r.particles[1]),
        pbits(&want),
        "the Pal path of add_fmt must reduce to reb_particle_from_pal"
    );
    // and the Pal elements read back
    let (mut a2, mut l2, mut k2, mut h2, mut ix2, mut iy2) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(
        1.0,
        r.particles[1],
        primary,
        &mut a2,
        &mut l2,
        &mut k2,
        &mut h2,
        &mut ix2,
        &mut iy2,
    );
    assert!(relerr(a2, 1.6) < 1e-11, "Pal a {} vs 1.6", a2);
    assert!(ang_diff(l2, 0.9) < 1e-10, "Pal lambda {} vs 0.9", l2);
    assert!((k2 - 0.12).abs() < 1e-11, "Pal k {} vs 0.12", k2);
    assert!((h2 + 0.05).abs() < 1e-11, "Pal h {} vs -0.05", h2);
    assert!((ix2 - 0.3).abs() < 1e-11, "Pal ix {} vs 0.3", ix2);
    assert!((iy2 + 0.2).abs() < 1e-11, "Pal iy {} vs -0.2", iy2);
}

#[test]
fn add_fmt_registers_and_finds_names() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m name", &[d(1.0), reb_fmt_arg::name("Sun".into())]);
    reb_simulation_add_fmt(
        &mut r,
        "m a name",
        &[d(1e-3), d(5.2), reb_fmt_arg::name("Jupiter".into())],
    );
    assert_eq!(r.name_list, vec!["Sun".to_string(), "Jupiter".to_string()]);
    assert_eq!(r.particles[0].name, Some(0));
    assert_eq!(r.particles[1].name, Some(1));
    assert_eq!(
        reb_simulation_get_particle_by_name(&r, "Jupiter"),
        Some(1),
        "the named particle must be findable"
    );
    assert_eq!(reb_simulation_get_particle_by_name(&r, "Pluto"), None);
}

#[test]
fn add_fmt_builtin_datasets() {
    let mut outer = reb_simulation_create();
    outer.G = 1.0;
    outer.save_messages = 1;
    reb_simulation_add_fmt(&mut outer, "outer solar system", &[]);
    assert_eq!(outer.N, 5, "Sun + Jupiter..Neptune");
    let want_outer = [0usize, 5, 6, 7, 8];
    for (i, &src) in want_outer.iter().enumerate() {
        assert_eq!(
            pbits(&outer.particles[i]),
            pbits(&reb_particle_solarsystem[src]),
            "outer solar system particle {} must be reb_particle_solarsystem[{}]",
            i,
            src
        );
    }
    // whitespace in the dataset name is ignored
    let mut outer2 = reb_simulation_create();
    outer2.G = 1.0;
    outer2.save_messages = 1;
    reb_simulation_add_fmt(&mut outer2, "outersolarsystem", &[]);
    assert_eq!(outer2.N, 5, "the dataset name ignores whitespace");

    let mut full = reb_simulation_create();
    full.G = 1.0;
    full.save_messages = 1;
    reb_simulation_add_fmt(&mut full, "solar system", &[]);
    assert_eq!(full.N, 9, "Sun + eight planets");
    for i in 0..9 {
        assert_eq!(
            pbits(&full.particles[i]),
            pbits(&reb_particle_solarsystem[i]),
            "solar system particle {}",
            i
        );
    }
    // the Sun really dominates, and the whole thing is bound
    reb_simulation_move_to_com(&mut full);
    assert!(
        reb_simulation_energy(&full) < 0.0,
        "the Solar System must be bound"
    );
    // Jupiter's semi-major axis is ~5.2 AU
    let sun = full.particles[0];
    let o = reb_orbit_from_particle(full.G, full.particles[5], sun);
    assert!(
        (o.a - 5.2).abs() < 0.1,
        "Jupiter's semi-major axis came out as {} AU",
        o.a
    );
    // G != 1 raises a warning
    let mut warn = reb_simulation_create();
    warn.G = 2.0;
    warn.save_messages = 1;
    reb_simulation_add_fmt(&mut warn, "solar system", &[]);
    assert!(
        warn.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::WARNING && m.contains("G should be 1.0")),
        "a non-unit G must warn, got {:?}",
        warn.messages
    );
}

#[test]
fn add_fmt_error_paths_report_and_add_nothing() {
    // (format, args, expected substring of the error message)
    struct Case {
        fmt: &'static str,
        args: Vec<reb_fmt_arg>,
        msg: &'static str,
    }
    let cases = vec![
        Case {
            fmt: "m a e h",
            args: vec![d(1e-3), d(1.0), d(0.1), d(0.05)],
            msg: "Cannot mix Pal coordinates",
        },
        Case {
            fmt: "m x a",
            args: vec![d(1e-3), d(1.0), d(2.0)],
            msg: "Cannot pass cartesian coordinates and orbital elements",
        },
        Case {
            fmt: "m e",
            args: vec![d(1e-3), d(0.2)],
            msg: "Need to pass either semi-major axis or orbital period to initialize",
        },
        Case {
            fmt: "m a P",
            args: vec![d(1e-3), d(1.0), d(1.0)],
            msg: "but not both",
        },
        Case {
            fmt: "m a ix iy",
            args: vec![d(1e-3), d(1.0), d(2.0), d(2.0)],
            msg: "Squared sum exceeds 4",
        },
        Case {
            fmt: "m a omega pomega",
            args: vec![d(1e-3), d(1.0), d(0.1), d(0.2)],
            msg: "Cannot pass both (omega, pomega)",
        },
        Case {
            fmt: "m a f M",
            args: vec![d(1e-3), d(1.0), d(0.1), d(0.2)],
            msg: "Can only pass one longitude/anomaly",
        },
        Case {
            fmt: "m a e",
            args: vec![d(1e-3), d(1.0), d(1.0)],
            msg: "Cannot set e exactly to 1",
        },
        Case {
            fmt: "m a e",
            args: vec![d(1e-3), d(1.0), d(1.5)],
            msg: "Bound orbit (a > 0) must have e < 1",
        },
        Case {
            fmt: "m a e",
            args: vec![d(1e-3), d(-1.0), d(0.5)],
            msg: "Unbound orbit (a < 0) must have e > 1",
        },
    ];
    for (i, c) in cases.iter().enumerate() {
        let mut r = reb_simulation_create();
        r.G = 1.0;
        r.save_messages = 1;
        reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]); // a primary with mass
        r.messages.clear();
        reb_simulation_add_fmt(&mut r, c.fmt, &c.args);
        assert_eq!(
            r.N, 1,
            "case {} ('{}'): a rejected particle must not be added",
            i, c.fmt
        );
        assert!(
            r.messages
                .iter()
                .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m.contains(c.msg)),
            "case {} ('{}'): expected an error containing {:?}, got {:?}",
            i,
            c.fmt,
            c.msg,
            r.messages
        );
    }

    // A massless primary (an empty simulation) is error 6.
    let mut empty = reb_simulation_create();
    empty.G = 1.0;
    empty.save_messages = 1;
    reb_simulation_add_fmt(&mut empty, "m a e", &[d(1e-3), d(1.0), d(0.1)]);
    assert_eq!(empty.N, 0, "nothing may be added without a primary");
    assert!(
        empty
            .messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m.contains("Primary has no mass")),
        "expected a 'Primary has no mass' error, got {:?}",
        empty.messages
    );
}

// ===========================================================================
// 11. variational particles and MEGNO
// ===========================================================================

#[test]
fn add_variation_allocates_the_right_number_of_particles() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    for _ in 0..3 {
        reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(1.0)]);
    }
    assert_eq!(r.N, 3);
    let i1 = reb_simulation_add_variation_1st_order(&mut r, -1);
    assert_eq!(i1, 0, "the first variational set starts at index 0");
    assert_eq!(r.N_var, 3, "a full first-order set is N particles");
    assert_eq!(r.particles_var.len(), 3);
    assert_eq!(r.var_config.len(), 1);
    assert_eq!(r.var_config[0].order, 1);
    assert_eq!(r.var_config[0].testparticle, -1);
    assert_eq!(r.var_config[0].lrescale, 0.0);

    // a test-particle variation only needs one entry
    let i2 = reb_simulation_add_variation_1st_order(&mut r, 2);
    assert_eq!(i2, 3, "the second set starts where the first ended");
    assert_eq!(r.N_var, 4);
    assert_eq!(r.var_config[1].testparticle, 2);

    // second order
    let i3 = reb_simulation_add_variation_2nd_order(&mut r, -1, i1, i2);
    assert_eq!(i3, 4);
    assert_eq!(r.N_var, 7);
    assert_eq!(r.var_config[2].order, 2);
    assert_eq!(r.var_config[2].index_1st_order_a, i1);
    assert_eq!(r.var_config[2].index_1st_order_b, i2);
    // freshly allocated variational particles are all zero
    for i in 0..r.N_var {
        assert_eq!(
            pbits(&r.particles_var[i]),
            pbits(&reb_particle::default()),
            "variational particle {} must start at zero",
            i
        );
    }
}

#[test]
fn init_megno_seed_is_deterministic_and_normalises_the_deviation_vectors() {
    let build = || {
        let mut r = reb_simulation_create();
        r.G = 1.0;
        r.save_messages = 1;
        reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
        reb_simulation_add_fmt(&mut r, "m a e", &[d(1e-3), d(1.0), d(0.1)]);
        reb_simulation_add_fmt(&mut r, "m a e", &[d(1e-3), d(2.0), d(0.2)]);
        r
    };
    let mut a = build();
    let mut b = build();
    reb_simulation_init_megno_seed(&mut a, 42);
    reb_simulation_init_megno_seed(&mut b, 42);
    assert_eq!(a.N_var, 3, "MEGNO adds one first-order set");
    assert_eq!(a.calculate_megno, 1);
    assert_eq!(a.megno_initial_t, a.t);
    for i in 0..a.N_var {
        assert_eq!(
            pbits(&a.particles_var[i]),
            pbits(&b.particles_var[i]),
            "the same seed must give bit-identical deviation vectors (particle {})",
            i
        );
        let p = a.particles_var[i];
        assert_eq!(p.m, 0.0, "deviation vectors carry no mass");
        let n =
            (p.x * p.x + p.y * p.y + p.z * p.z + p.vx * p.vx + p.vy * p.vy + p.vz * p.vz).sqrt();
        assert!(
            (n - 1.0).abs() < 1e-14,
            "deviation vector {} has 6-norm {}, must be normalised to 1",
            i,
            n
        );
    }
    // a different seed gives a different vector
    let mut c = build();
    reb_simulation_init_megno_seed(&mut c, 43);
    assert_ne!(
        pbits(&c.particles_var[0]),
        pbits(&a.particles_var[0]),
        "a different seed must give a different deviation vector"
    );
    // MEGNO is exactly zero before any time has passed
    assert_eq!(
        reb_simulation_megno(&a),
        0.0,
        "MEGNO must be 0 at t == megno_initial_t"
    );
    assert_eq!(
        reb_simulation_lyapunov(&a),
        0.0,
        "the Lyapunov estimate must be 0 while megno_var_t == 0"
    );
}

#[test]
fn megno_of_a_regular_two_body_orbit_tends_to_two() {
    // A Kepler orbit is integrable, so MEGNO -> 2 and the Lyapunov estimate
    // -> 0 (Cincotta & Simo 2000).
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m", &[d(1.0)]);
    reb_simulation_add_fmt(&mut r, "m a e", &[d(1e-3), d(1.0), d(0.1)]);
    reb_simulation_move_to_com(&mut r);
    reb_simulation_set_integrator(&mut r, "whfast");
    let period = 2.0 * PI;
    r.dt = period / 100.0;
    reb_simulation_init_megno_seed(&mut r, 4242);
    reb_simulation_integrate(&mut r, 400.0 * period);
    let megno = reb_simulation_megno(&r);
    assert!(
        (megno - 2.0).abs() < 0.05,
        "MEGNO of a regular Kepler orbit came out as {}, expected ~2",
        megno
    );
    let lyap = reb_simulation_lyapunov(&r);
    assert!(
        lyap.abs() < 5e-4,
        "the Lyapunov estimate of a regular orbit came out as {}, expected ~0",
        lyap
    );
    assert!(r.megno_n > 100, "MEGNO must have been updated many times");
}

#[test]
fn megno_update_keeps_a_correct_running_mean() {
    // reb_tools_megno_update advances the Welford mean of (t - t0); after n
    // updates megno_mean_t must be the plain arithmetic mean of the samples.
    let mut r = reb_simulation_create();
    r.megno_initial_t = 0.5;
    r.t = 0.5;
    let mut samples: Vec<f64> = Vec::new();
    let mut ysum = 0.0;
    for i in 1..=40 {
        r.t = 0.5 + 0.25 * (i as f64);
        samples.push(r.t - r.megno_initial_t);
        ysum += 0.01;
        reb_tools_megno_update(&mut r, 0.01, 0.25);
        assert_eq!(
            r.megno_n, i as i64,
            "megno_n must count the updates"
        );
        // megno_Ys is the running sum of dY
        assert!(
            relerr(r.megno_Ys, ysum) < 1e-13,
            "megno_Ys {} vs sum of dY {}",
            r.megno_Ys,
            ysum
        );
    }
    let mean: f64 = samples.iter().sum::<f64>() / (samples.len() as f64);
    assert!(
        relerr(r.megno_mean_t, mean) < 1e-13,
        "megno_mean_t {} vs the arithmetic mean {}",
        r.megno_mean_t,
        mean
    );
    // reb_simulation_megno is Yss/(t - t0) by definition
    assert!(
        relerr(reb_simulation_megno(&r), r.megno_Yss / (r.t - r.megno_initial_t)) < 1e-15,
        "reb_simulation_megno must be Yss/(t-t0)"
    );
    // and the Lyapunov estimate is cov/var
    assert!(
        r.megno_var_t > 0.0,
        "megno_var_t must have accumulated, got {}",
        r.megno_var_t
    );
    assert!(
        relerr(reb_simulation_lyapunov(&r), r.megno_cov_Yt / r.megno_var_t) < 1e-15,
        "reb_simulation_lyapunov must be cov_Yt/var_t"
    );
}

#[test]
fn megno_deltad_delta_is_the_ratio_of_the_two_inner_products() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(0.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(1.0)]);
    reb_simulation_add_variation_1st_order(&mut r, -1);
    // Choose exact binary values so the ratio is exactly computable.
    r.particles_var[0] = reb_particle {
        x: 1.0,
        y: 0.0,
        z: 0.0,
        vx: 2.0,
        vy: 0.0,
        vz: 0.0,
        ax: 0.5,
        ay: 0.0,
        az: 0.0,
        ..reb_particle::default()
    };
    r.particles_var[1] = reb_particle {
        x: 0.0,
        y: 4.0,
        z: 0.0,
        vx: 0.0,
        vy: 0.25,
        vz: 0.0,
        ay: 8.0,
        ..reb_particle::default()
    };
    // deltad = sum(v.x + a.v) = (2*1 + 0.5*2) + (0.25*4 + 8*0.25) = 3 + 3 = 6
    // delta2 = sum(x^2 + v^2)  = (1 + 4)      + (16 + 0.0625)     = 5 + 16.0625
    let want = 6.0 / 21.0625;
    let got = reb_tools_megno_deltad_delta(&r);
    assert!(
        relerr(got, want) < 1e-15,
        "megno_deltad_delta {} vs the hand-computed ratio {}",
        got,
        want
    );
}

#[test]
fn rescale_var_normalises_a_diverged_deviation_vector() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(0.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(1.0)]);
    reb_simulation_add_variation_1st_order(&mut r, -1);
    r.is_synchronized = 1;

    // below the 1e100 threshold nothing happens at all
    r.particles_var[0].x = 1e50;
    r.particles_var[1].vy = 2e40;
    let before: Vec<[u64; 11]> = r.particles_var.iter().map(pbits).collect();
    r.did_modify_particles = 0;
    reb_simulation_rescale_var(&mut r);
    let after: Vec<[u64; 11]> = r.particles_var.iter().map(pbits).collect();
    assert_eq!(before, after, "no rescaling below the 1e100 threshold");
    assert_eq!(r.var_config[0].lrescale, 0.0);
    assert_eq!(r.did_modify_particles, 0);

    // above it, everything is divided by the largest component
    let scale = 4e105_f64;
    r.particles_var[0].x = scale;
    r.particles_var[0].vy = scale / 4.0;
    r.particles_var[1].z = -scale / 2.0;
    reb_simulation_rescale_var(&mut r);
    assert_eq!(
        r.particles_var[0].x, 1.0,
        "the largest component must become exactly 1"
    );
    assert_eq!(r.particles_var[0].vy, 0.25, "scale/4 / scale = 0.25 exactly");
    assert_eq!(r.particles_var[1].z, -0.5, "-scale/2 / scale = -0.5 exactly");
    assert_eq!(
        r.var_config[0].lrescale.to_bits(),
        scale.ln().to_bits(),
        "lrescale must accumulate ln(scale)"
    );
    assert_eq!(
        r.did_modify_particles, 1,
        "rescaling must flag the particles as modified"
    );
}

#[test]
fn rescale_var_refuses_when_the_integrator_is_unsynchronized() {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.save_messages = 1;
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(0.0)]);
    reb_simulation_add_fmt(&mut r, "m x", &[d(1.0), d(1.0)]);
    reb_simulation_add_variation_1st_order(&mut r, -1);
    r.particles_var[0].x = 1e110;
    r.is_synchronized = 0;
    r.messages.clear();
    let before: Vec<[u64; 11]> = r.particles_var.iter().map(pbits).collect();
    reb_simulation_rescale_var(&mut r);
    let after: Vec<[u64; 11]> = r.particles_var.iter().map(pbits).collect();
    assert_eq!(
        before, after,
        "an unsynchronized integrator must block the rescaling"
    );
    assert_eq!(r.var_config[0].lrescale, 0.0);
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::WARNING && m.contains("not synchronized")),
        "a warning must be recorded, got {:?}",
        r.messages
    );
}

// ===========================================================================
// 12. Plummer sphere and determinism
// ===========================================================================

#[test]
fn add_plummer_is_deterministic_and_has_the_requested_mass() {
    let build = |seed: u32| {
        let mut r = reb_simulation_create();
        r.G = 1.0;
        r.save_messages = 1;
        r.rand_seed = seed;
        reb_simulation_add_plummer(&mut r, 60, 1.0, 1.0);
        r
    };
    let a = build(7);
    let b = build(7);
    assert_eq!(a.N, 60, "60 stars must be added");
    for i in 0..a.N {
        assert_eq!(
            pbits(&a.particles[i]),
            pbits(&b.particles[i]),
            "the same rand_seed must give a bit-identical Plummer sphere (star {})",
            i
        );
        assert!(
            a.particles[i].x.is_finite()
                && a.particles[i].y.is_finite()
                && a.particles[i].z.is_finite()
                && a.particles[i].vx.is_finite()
                && a.particles[i].vy.is_finite()
                && a.particles[i].vz.is_finite(),
            "star {} has a non-finite coordinate",
            i
        );
        assert_eq!(
            a.particles[i].m.to_bits(),
            (1.0f64 / 60.0).to_bits(),
            "every star carries M/N"
        );
    }
    let c = build(8);
    assert_ne!(
        pbits(&c.particles[0]),
        pbits(&a.particles[0]),
        "a different seed must give a different sphere"
    );
    // total mass and the virial-scale sanity check
    let mtot: f64 = a.particles.iter().map(|p| p.m).sum();
    assert!(relerr(mtot, 1.0) < 1e-14, "total mass {} vs 1.0", mtot);
    let com = reb_simulation_com(&a);
    assert!(
        com.x.abs() < 0.5 && com.y.abs() < 0.5 && com.z.abs() < 0.5,
        "the Plummer sphere's COM ({}, {}, {}) drifted far from the origin",
        com.x,
        com.y,
        com.z
    );
    // The model is bound, and stays bound after a shift into the COM frame
    // (which can only lower the kinetic energy).
    let e_lab = reb_simulation_energy(&a);
    assert!(
        e_lab < 0.0,
        "a Plummer sphere must have negative total energy, got {}",
        e_lab
    );
    let mut shifted = build(7);
    reb_simulation_move_to_com(&mut shifted);
    let e_com = reb_simulation_energy(&shifted);
    assert!(
        e_com < 0.0 && e_com <= e_lab * (1.0 - 1e-12),
        "COM-frame energy {} must be negative and no larger than the lab-frame {}",
        e_com,
        e_lab
    );
    // the removed piece is exactly the bulk kinetic energy 1/2 M_tot v_com^2
    let bulk = 0.5 * mtot * (com.vx * com.vx + com.vy * com.vy + com.vz * com.vz);
    assert!(
        relerr(e_lab - e_com, bulk) < 1e-9,
        "the lab/COM energy difference {} must be the bulk kinetic energy {}",
        e_lab - e_com,
        bulk
    );
}

#[test]
fn identical_simulations_integrate_to_bit_identical_states() {
    let run = || {
        let mut r = three_body();
        reb_simulation_set_integrator(&mut r, "whfast");
        r.dt = 1.0 / 137.0;
        reb_simulation_integrate(&mut r, 37.0);
        r
    };
    let a = run();
    let b = run();
    assert_eq!(a.N, b.N);
    for i in 0..a.N {
        assert_eq!(
            pbits(&a.particles[i]),
            pbits(&b.particles[i]),
            "two identical runs must agree bit for bit (particle {})",
            i
        );
    }
    assert_eq!(
        reb_simulation_energy(&a).to_bits(),
        reb_simulation_energy(&b).to_bits(),
        "the energy of two identical runs must agree bit for bit"
    );
    assert_eq!(a.t.to_bits(), b.t.to_bits());
    assert_eq!(a.steps_done, b.steps_done);
}

// ===========================================================================
// 13. miscellaneous tools
// ===========================================================================

#[test]
fn particle_nan_is_nan_everywhere() {
    let p = reb_particle_nan();
    assert!(
        p.x.is_nan()
            && p.y.is_nan()
            && p.z.is_nan()
            && p.vx.is_nan()
            && p.vy.is_nan()
            && p.vz.is_nan()
            && p.ax.is_nan()
            && p.ay.is_nan()
            && p.az.is_nan()
            && p.m.is_nan()
            && p.r.is_nan(),
        "every field of the NaN particle must be NaN: {:?}",
        p
    );
    assert_eq!(p.name, None);
}

#[test]
fn messages_are_stored_and_capped_when_save_messages_is_set() {
    let mut r = reb_simulation_create();
    r.save_messages = 1;
    for i in 0..(reb_messages_max_N + 5) {
        reb_simulation_warning(&mut r, &format!("warning {}", i));
    }
    assert_eq!(
        r.messages.len(),
        reb_messages_max_N,
        "the message ring must be capped at reb_messages_max_N"
    );
    // the oldest were dropped: the buffer ends with the newest message
    assert_eq!(
        r.messages.last().unwrap().1,
        format!("warning {}", reb_messages_max_N + 4)
    );
    assert_eq!(r.messages.last().unwrap().0, REB_MESSAGE_TYPE::WARNING);
    reb_simulation_error(&mut r, "boom");
    assert_eq!(r.messages.last().unwrap().0, REB_MESSAGE_TYPE::ERROR);
    reb_simulation_info(&mut r, "hello");
    assert_eq!(r.messages.last().unwrap().0, REB_MESSAGE_TYPE::INFO);
    assert_eq!(r.messages.last().unwrap().1, "hello");
}
