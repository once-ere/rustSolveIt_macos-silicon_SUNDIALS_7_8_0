//! archive_test.rs — Rust side of the Simulationarchive cross-language
//! round-trip verification (twin of porttest/archive_test.c, with the
//! file roles swapped: writes archive_rust_<integrator>.bin, continues
//! from archive_c_<integrator>.bin, dumps to archive_state_rust.txt).
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

fn setup_particles(r: &mut reb_simulation) {
    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(r, star);
    let mut planet = reb_particle::default();
    planet.m = 1e-3;
    planet.x = 1.6;
    planet.vy = 0.5;
    reb_simulation_add(r, planet);
    let mut moon = reb_particle::default();
    moon.m = 1e-7;
    moon.x = 1.7;
    moon.vy = 0.6;
    moon.z = 0.01;
    moon.vz = 0.001;
    reb_simulation_add(r, moon);
}

fn configure(r: &mut reb_simulation, integrator: &str) {
    if integrator == "whfast-usafe" {
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            wh.safe_mode = 0;
        }
    }
}

fn dump_state(r: &reb_simulation, integrator: &str) {
    let mut f = std::fs::File::create("archive_state_rust.txt").unwrap();
    writeln!(f, "integrator {}", integrator).unwrap();
    writeln!(f, "t {:016x}", bits(r.t)).unwrap();
    writeln!(f, "dt {:016x}", bits(r.dt)).unwrap();
    for i in 0..r.N {
        let p = r.particles[i];
        writeln!(
            f,
            "{} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
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
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let integrator = if args.len() > 1 { args[1].clone() } else { "whfast-usafe".to_string() };
    let mode = if args.len() > 2 { args[2].clone() } else { "write".to_string() };

    let real_integrator = if integrator.starts_with("whfast") { "whfast" } else { integrator.as_str() };

    if mode == "write" {
        let mut sim = reb_simulation_create();
        let r = &mut sim;
        reb_simulation_set_integrator(r, real_integrator);
        configure(r, &integrator);
        r.G = 1.0;
        r.dt = 0.01;
        setup_particles(r);
        let fname = format!("archive_rust_{}.bin", integrator);
        let _ = std::fs::remove_file(&fname);
        for _ in 0..3 {
            reb_simulation_steps(r, 100);
            reb_simulation_save_to_file(r, Some(&fname));
        }
        dump_state(r, &integrator);
        println!("write done: t={:e}", r.t);
    } else {
        let fname = format!("archive_c_{}.bin", integrator);
        let mut sim = match reb_simulation_create_from_file(&fname, 1) {
            Some(s) => s,
            None => {
                println!("Failed to load {}", fname);
                std::process::exit(1);
            }
        };
        let r = &mut sim;
        reb_simulation_steps(r, 100);
        dump_state(r, &integrator);
        println!("continue done: t={:e}", r.t);
    }
}
