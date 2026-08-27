//! Rust twin of porttest/movetocom_var_test.c — probes the 1st-order
//! variational block of reb_simulation_move_to_com.
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

fn main() {
    let mut r = reb_simulation_create();
    r.G = 1.;
    let mut sun = reb_particle::default();
    sun.m = 1.;
    reb_simulation_add(&mut r, sun);
    let mut jup = reb_particle::default();
    jup.m = 0.000954588;
    jup.x = 5.2;
    jup.vy = 0.4396;
    reb_simulation_add(&mut r, jup);

    reb_simulation_init_megno_seed(&mut r, 12345);

    let com = reb_simulation_com(&r);
    println!("com.m = {:.17e} (bits {:016x})", com.m, com.m.to_bits());
    println!("com.x = {:.17e}", com.x);
    println!(
        "N={} N_var={} N_var_config={} index={}",
        r.N,
        r.N_var,
        r.var_config.len(),
        r.var_config[0].index
    );
    for i in 0..r.N_var {
        println!(
            "BEFORE var[{}] m={:.17e} x={:.17e} vx={:.17e}",
            i, r.particles_var[i].m, r.particles_var[i].x, r.particles_var[i].vx
        );
    }
    let mut dm_real = 0.;
    let mut dm_var = 0.;
    for i in 0..r.N {
        dm_real += r.particles[i].m;
        dm_var += r.particles_var[i].m;
    }
    println!("dm(particles)={:.17e}  dm(particles_var)={:.17e}", dm_real, dm_var);

    reb_simulation_move_to_com(&mut r);

    for i in 0..r.N_var {
        println!(
            "AFTER  var[{}] x={:.17e} (bits {:016x}) vx={:.17e} (bits {:016x})",
            i,
            r.particles_var[i].x,
            r.particles_var[i].x.to_bits(),
            r.particles_var[i].vx,
            r.particles_var[i].vx.to_bits()
        );
    }
}
