//! shearing_sheet port test — Rust side. Mirrors porttest/problem_test.c
//! exactly: seed 42, no server, no heartbeat, 400 steps, raw-bit dump.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein and contributors.
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
use std::f64::consts::PI as M_PI;
use std::io::Write;

/// Same Bridges law as the C harness twin: the single `pow` is written
/// via exp/log on BOTH sides because `pow` is the one libm function
/// where Rust's runtime does not defer to the UCRT (see
/// shearing_sheet_port_test.md). The stock example keeps `powf`.
fn coefficient_of_restitution_bridges(_r: &reb_simulation, v: f64) -> f64 {
    let mut eps = 0.32 * (-0.234 * (v.abs() * 100.).ln()).exp();
    if eps > 1. {
        eps = 1.;
    }
    if eps < 0. {
        eps = 0.;
    }
    eps
}

fn bits(x: f64) -> u64 {
    x.to_bits()
}

/// C's `%.17e` (17 fractional digits, sign-carrying two-digit exponent).
fn fmt_e17(x: f64) -> String {
    let s = format!("{:.17e}", x);
    if let Some(pos) = s.rfind('e') {
        let (mantissa, exp) = s.split_at(pos);
        let exp = &exp[1..];
        let (sign, digits) = if let Some(stripped) = exp.strip_prefix('-') {
            ("-", stripped)
        } else {
            ("+", exp)
        };
        if digits.len() < 2 {
            format!("{}e{}0{}", mantissa, sign, digits)
        } else {
            format!("{}e{}{}", mantissa, sign, digits)
        }
    } else {
        s
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nsteps: usize = if args.len() > 1 {
        args[1].parse().unwrap_or(400)
    } else {
        400
    };
    let mut sim = reb_simulation_create();
    let r = &mut sim;
    r.rand_seed = 42; // CONTROLLED SEED
    r.opening_angle2 = 0.5;
    reb_simulation_set_integrator(r, "sei");
    r.boundary = REB_BOUNDARY::SHEAR;
    r.gravity = REB_GRAVITY::TREE;
    r.collision = REB_COLLISION::TREE;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.OMEGA = 0.00013143527; // 1/s
    r.G = 6.67428e-11; // N / (1e-5 kg)^2 m^2
    r.softening = 0.1; // m
    r.dt = 1e-3 * 2. * M_PI / r.OMEGA; // s
    let surfacedensity = 400.; // kg/m^2
    let particle_density = 400.; // kg/m^3
    let particle_radius_min = 1.; // m
    let particle_radius_max = 4.; // m
    let particle_radius_slope = -3.;
    let root_size = 100.; // m
    r.root_size = root_size;
    r.N_root_x = 2;
    r.N_root_y = 2;
    r.N_ghost_x = 2;
    r.N_ghost_y = 2;
    r.N_ghost_z = 0;
    let boxsize = reb_vec3d {
        x: r.root_size * (r.N_root_x as f64),
        y: r.root_size * (r.N_root_y as f64),
        z: r.root_size * (r.N_root_z as f64),
    };

    println!(
        "Toomre wavelength: {:.6}",
        4. * M_PI * M_PI * surfacedensity / r.OMEGA / r.OMEGA * r.G
    );
    r.coefficient_of_restitution = Some(coefficient_of_restitution_bridges);
    r.minimum_collision_velocity = particle_radius_min * r.OMEGA * 0.001;

    // Add all ring particles
    let total_mass = surfacedensity * boxsize.x * boxsize.y;
    let mut mass = 0.;
    while mass < total_mass {
        let mut pt = reb_particle::default();
        pt.x = reb_random_uniform(Some(r), -boxsize.x / 2., boxsize.x / 2.);
        pt.y = reb_random_uniform(Some(r), -boxsize.y / 2., boxsize.y / 2.);
        pt.z = reb_random_normal(Some(r), 1.); // m
        pt.vx = 0.;
        pt.vy = -1.5 * pt.x * r.OMEGA;
        pt.vz = 0.;
        pt.ax = 0.;
        pt.ay = 0.;
        pt.az = 0.;
        let radius = reb_random_powerlaw(
            Some(r),
            particle_radius_min,
            particle_radius_max,
            particle_radius_slope,
        );
        pt.r = radius; // m
        let particle_mass = particle_density * 4. / 3. * M_PI * radius * radius * radius;
        pt.m = particle_mass; // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
    }

    println!("N after init: {}", r.N);

    // Dump the initial conditions too (verifies the RNG stream match).
    let mut f0 = std::fs::File::create("state_rust_init.txt").unwrap();
    writeln!(f0, "N {}", r.N).unwrap();
    for i in 0..r.N {
        let p = r.particles[i];
        writeln!(
            f0,
            "{} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            i,
            bits(p.x),
            bits(p.y),
            bits(p.z),
            bits(p.vx),
            bits(p.vy),
            bits(p.vz),
            bits(p.m),
            bits(p.r)
        )
        .unwrap();
    }
    drop(f0);

    reb_simulation_steps(r, nsteps);

    let mut f = std::fs::File::create("state_rust_final.txt").unwrap();
    writeln!(f, "N {}", r.N).unwrap();
    writeln!(f, "t {:016x} {}", bits(r.t), fmt_e17(r.t)).unwrap();
    writeln!(f, "steps_done {}", r.steps_done).unwrap();
    writeln!(f, "collisions_log_n {}", r.collisions_log_n).unwrap();
    writeln!(
        f,
        "collisions_plog {:016x} {}",
        bits(r.collisions_plog),
        fmt_e17(r.collisions_plog)
    )
    .unwrap();
    writeln!(f, "rand_seed {}", r.rand_seed).unwrap();
    for i in 0..r.N {
        let p = r.particles[i];
        writeln!(
            f,
            "{} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
            i,
            bits(p.x),
            bits(p.y),
            bits(p.z),
            bits(p.vx),
            bits(p.vy),
            bits(p.vz),
            bits(p.m),
            bits(p.r)
        )
        .unwrap();
    }
    drop(f);
    println!(
        "final: t={} steps={} collisions={}",
        fmt_e17(r.t),
        r.steps_done,
        r.collisions_log_n
    );
}
