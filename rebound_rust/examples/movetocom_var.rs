//! movetocom_var.rs — Rust twin of `porttest/movetocom_var_c.c`.
//!
//! Regression probe for the variational centre-of-mass shift. The port
//! audit found that the first-order `dm` accumulator summed
//! `particles_var` where the C sums `particles`, which dropped a whole
//! term of the shift and changed every MEGNO / Lyapunov result. This
//! program dumps the post-`move_to_com` state as raw IEEE-754 bits so it
//! can be compared byte-for-byte with the C reference.
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

use rebound_rs::*;
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let mut sim = reb_simulation_create();
    let r = &mut sim;
    r.G = 1.0;
    r.dt = 0.01;
    reb_simulation_set_integrator(r, "ias15");

    let mut sun = reb_particle::default();
    sun.m = 1.0;
    reb_simulation_add(r, sun);
    reb_simulation_add_fmt(
        r,
        "m a e",
        &[
            reb_fmt_arg::d(9.54579e-4),
            reb_fmt_arg::d(5.2),
            reb_fmt_arg::d(0.0489),
        ],
    );

    reb_simulation_init_megno_seed(r, 12345);

    let com = reb_simulation_com(r);
    println!("com.m {:016x} com.x {:016x}", bits(com.m), bits(com.x));

    reb_simulation_move_to_com(r);

    let mut f = std::fs::File::create("movetocom_var_rust.txt").unwrap();
    writeln!(f, "N {} N_var {}", r.N, r.N_var).unwrap();
    for i in 0..r.N {
        let p = r.particles[i];
        writeln!(
            f,
            "p {} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            i,
            bits(p.x),
            bits(p.y),
            bits(p.z),
            bits(p.vx),
            bits(p.vy),
            bits(p.vz)
        )
        .unwrap();
    }
    for i in 0..r.N_var {
        let p = r.particles_var[i];
        writeln!(
            f,
            "v {} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            i,
            bits(p.x),
            bits(p.y),
            bits(p.z),
            bits(p.vx),
            bits(p.vy),
            bits(p.vz)
        )
        .unwrap();
    }
    println!("movetocom_var_rust done");
}
