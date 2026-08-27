//! frequency_test.rs — Rust twin of porttest/frequency_test.c.
//! Part of rebound_rs, GPL-3.0-or-later.
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

use rebound_rs::frequency_analysis::*;
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let ndata: usize = 256;
    let nfreq: usize = 3;
    let mut input = vec![0.0_f64; 2 * ndata];
    // Quasi-periodic signal with three frequencies (rad per sample).
    let (f1, a1, p1) = (0.30_f64, 1.00_f64, 0.40_f64);
    let (f2, a2, p2) = (0.55_f64, 0.35_f64, 1.90_f64);
    let (f3, a3, p3) = (0.11_f64, 0.10_f64, 5.10_f64);
    for i in 0..ndata {
        let t = i as f64;
        input[2 * i] = a1 * (f1 * t + p1).cos() + a2 * (f2 * t + p2).cos() + a3 * (f3 * t + p3).cos();
        input[2 * i + 1] =
            a1 * (f1 * t + p1).sin() + a2 * (f2 * t + p2).sin() + a3 * (f3 * t + p3).sin();
    }

    let mut f = std::fs::File::create("frequency_rust.txt").unwrap();
    let types = [
        REB_FREQUENCY_ANALYSIS_MFT,
        REB_FREQUENCY_ANALYSIS_FMFT,
        REB_FREQUENCY_ANALYSIS_FMFT2,
    ];
    let names = ["MFT", "FMFT", "FMFT2"];
    for ti in 0..3 {
        let mut output = [0.0_f64; 9];
        let ret = reb_frequency_analysis(&mut output, nfreq, 0.05, 1.0, types[ti], &input, ndata);
        writeln!(f, "{} ret {}", names[ti], ret).unwrap();
        for k in 0..3 * nfreq {
            writeln!(f, "{} {:016x}", k, bits(output[k])).unwrap();
        }
    }
    println!("frequency_test done");
}
