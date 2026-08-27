//! Shearing sheet (Hill's approximation) — the Rust port of
//! rebound/examples/shearing_sheet/problem.c.
//!
//! Simulates a small patch of Saturn's rings in shearing sheet
//! coordinates, exactly as the C example: SEI integrator, shear
//! boundary, tree gravity, tree collisions, Bridges et al. velocity-
//! dependent coefficient of restitution.
//!
//! The C example also starts the visualization web server
//! (`reb_simulation_start_server(r, 1234)`); the server subsystem is
//! not part of this phase of the port, so this example integrates
//! headlessly — the physics is unchanged.
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

// This example is using a custom velocity dependent coefficient of restitution
fn coefficient_of_restitution_bridges(_r: &reb_simulation, v: f64) -> f64 {
    // assumes v in units of [m/s]
    let mut eps = 0.32 * (v.abs() * 100.).powf(-0.234);
    if eps > 1. {
        eps = 1.;
    }
    if eps < 0. {
        eps = 0.;
    }
    eps
}

fn heartbeat(r: &mut reb_simulation) {
    if reb_simulation_output_check(r, 1e-1 * 2. * M_PI / r.OMEGA) {
        reb_simulation_output_timing(r, 0.);
    }
    if reb_simulation_output_check(r, 2. * M_PI / r.OMEGA) {
        //reb_simulation_output_ascii(r, "position.txt");
    }
}

fn main() {
    let mut sim = reb_simulation_create();
    let r = &mut sim;

    // Setup constants
    r.opening_angle2 = 0.5; // Precision of the tree code gravity calculation.
    reb_simulation_set_integrator(r, "sei");
    r.boundary = REB_BOUNDARY::SHEAR;
    r.gravity = REB_GRAVITY::TREE;
    r.collision = REB_COLLISION::TREE;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.OMEGA = 0.00013143527; // 1/s
    r.G = 6.67428e-11; // N / (1e-5 kg)^2 m^2
    r.softening = 0.1; // m
    r.dt = 1e-3 * 2. * M_PI / r.OMEGA; // s
    r.heartbeat = Some(heartbeat);
    // This example uses two root boxes in the x and y direction.
    let surfacedensity = 400.; // kg/m^2
    let particle_density = 400.; // kg/m^3
    let particle_radius_min = 1.; // m
    let particle_radius_max = 4.; // m
    let particle_radius_slope = -3.;
    let mut root_size = 100.; // m
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        // Try to read root_size from command line
        root_size = args[1].parse().unwrap_or(0.0);
    }
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
    let _ = boxsize.z;

    // Initial conditions
    println!(
        "Toomre wavelength: {}",
        crate_fmt_f(4. * M_PI * M_PI * surfacedensity / r.OMEGA / r.OMEGA * r.G)
    );
    // Use Bridges et al coefficient of restitution.
    r.coefficient_of_restitution = Some(coefficient_of_restitution_bridges);
    // Prevent particles from sinking into each other.
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
    reb_simulation_integrate(r, f64::INFINITY);
}

/// C's `printf("%f", x)` (the Toomre wavelength line).
fn crate_fmt_f(x: f64) -> String {
    format!("{:.6}", x)
}
