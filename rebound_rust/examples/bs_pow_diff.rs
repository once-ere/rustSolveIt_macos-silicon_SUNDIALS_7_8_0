//! bs_pow_diff.rs — Rust twin of porttest/bs_pow_diff.c.
//!
//! Evaluates exactly the `pow` calls the BS step-size controller makes
//! and dumps raw bit patterns, so the C and Rust results can be
//! compared bit-for-bit. See rebound_rust.md, "Why the BS integrator's
//! proposed timestep can differ by one ULP".
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

use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let stepControl1 = 0.65_f64;
    let _stepControl2 = 0.94_f64;
    let stepControl3 = 0.02_f64;
    let mut f = std::io::BufWriter::new(std::fs::File::create("bs_pow_rust.txt").unwrap());
    let mut n: i64 = 0;
    for k in 1..=8i32 {
        let e = 1.0 / ((2 * k + 1) as f64);
        writeln!(f, "P {} {:016x}", k, bits(stepControl3.powf(e))).unwrap();
        for i in 0..25000i64 {
            let error = 1e-8 * 10.0_f64.powf(8.0 * (i as f64) / 25000.0);
            let v = (error / stepControl1).powf(e);
            writeln!(f, "E {} {} {:016x}", k, i, bits(v)).unwrap();
            n += 1;
        }
    }
    println!("bs_pow_rust: {} samples", n);
}
