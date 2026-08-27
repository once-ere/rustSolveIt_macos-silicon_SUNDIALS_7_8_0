//! kepler_rectilinear.rs — Rust twin of
//! `porttest/kepler_rectilinear_c.c`.
//!
//! Calls the Kepler solver with (near-)rectilinear hyperbolic motion —
//! zero or almost-zero angular momentum, the regime where the
//! quartic/Newton iteration and its bisection fallback are most
//! strained — and prints markers before and after, so a
//! non-terminating loop is unambiguous, plus the final state as raw
//! bits for exact comparison with the C build.
//!
//!   argv[1] (optional) = vy, the transverse velocity.
//!     vy = 0     -> h = 0 exactly (purely radial)
//!     vy = 1e-12 -> h tiny but non-zero
//!
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

use rebound_rs::integrator_whfast::reb_integrator_whfast_kepler_solver;
use rebound_rs::*;
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let vy: f64 = if args.len() > 1 {
        args[1].parse().unwrap_or(0.0)
    } else {
        0.0
    };

    // r = (1,0,0), v = (3,vy,0), mu = 1.
    // h = r x v = (0,0,vy).  v^2 = 9 > 2*mu/r = 2, so hyperbolic.
    let mut p = reb_particle::default();
    p.x = 1.;
    p.y = 0.;
    p.z = 0.;
    p.vx = 3.;
    p.vy = vy;
    p.vz = 0.;

    println!("BEFORE: 20 kepler steps, vy={:.17e}", vy);
    std::io::stdout().flush().unwrap();

    let mut buf = [reb_particle::default(), p];
    let mut no_var: [reb_particle; 0] = [];
    let nsteps: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20);
    let dt: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.1);
    for _ in 0..nsteps {
        reb_integrator_whfast_kepler_solver(None, &mut buf, &mut no_var, 1, 1.0, dt);
    }
    let p = buf[1];

    println!(
        "AFTER: x={:016x} y={:016x} vx={:016x} vy={:016x}",
        bits(p.x),
        bits(p.y),
        bits(p.vx),
        bits(p.vy)
    );
}
