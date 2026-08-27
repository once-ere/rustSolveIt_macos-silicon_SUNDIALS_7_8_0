//! addfmt_test.rs — Rust twin of porttest/addfmt_test.c.
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
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let mut sim = reb_simulation_create();
    let r = &mut sim;
    r.G = 1.0;
    reb_simulation_add_fmt(r, "solar system", &[]);
    reb_simulation_add_fmt(
        r,
        "m a e inc Omega omega f",
        &[
            reb_fmt_arg::d(1e-9),
            reb_fmt_arg::d(12.5),
            reb_fmt_arg::d(0.3),
            reb_fmt_arg::d(0.2),
            reb_fmt_arg::d(0.6),
            reb_fmt_arg::d(1.1),
            reb_fmt_arg::d(2.5),
        ],
    );
    reb_simulation_add_fmt(
        r,
        "m a l h k ix iy",
        &[
            reb_fmt_arg::d(2e-9),
            reb_fmt_arg::d(15.5),
            reb_fmt_arg::d(0.7),
            reb_fmt_arg::d(0.05),
            reb_fmt_arg::d(-0.03),
            reb_fmt_arg::d(0.01),
            reb_fmt_arg::d(0.02),
        ],
    );
    reb_simulation_add_fmt(
        r,
        "m P e M",
        &[
            reb_fmt_arg::d(3e-9),
            reb_fmt_arg::d(100.0),
            reb_fmt_arg::d(0.1),
            reb_fmt_arg::d(0.5),
        ],
    );

    let mut f = std::fs::File::create("addfmt_rust.txt").unwrap();
    for i in 0..r.N {
        let p = r.particles[i];
        writeln!(
            f,
            "{} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            i,
            bits(p.m),
            bits(p.x),
            bits(p.y),
            bits(p.z),
            bits(p.vx),
            bits(p.vy),
            bits(p.vz)
        )
        .unwrap();
    }
    println!("addfmt done N={}", r.N);
}
