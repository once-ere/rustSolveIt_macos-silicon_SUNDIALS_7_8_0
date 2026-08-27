//! libm_diff.rs — Rust twin of porttest/libm_diff.c. Same xorshift
//! corpus, same functions, bit-pattern output for diffing.
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

struct Xs(u64);
impl Xs {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn main() {
    let mut s = Xs(88172645463325252u64);
    let mut f = std::io::BufWriter::new(std::fs::File::create("libm_rust.txt").unwrap());
    for _ in 0..200000 {
        let x = ((s.next() % 2000000000u64) as f64) / 1e6 - 1000.0;
        let y = ((s.next() % 2000000000u64) as f64) / 1e6 - 1000.0;
        let xp = x.abs() + 1e-9;
        writeln!(
            f,
            "{:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            bits(x.sin()),
            bits(x.cos()),
            bits(x.tan()),
            bits(y.atan2(x)),
            bits(xp.powf(-0.234)),
            bits(xp.sqrt()),
            bits(y % 3.7),
            bits((x / 100.).exp()),
            bits(xp.ln())
        )
        .unwrap();
    }
    println!("done");
}
