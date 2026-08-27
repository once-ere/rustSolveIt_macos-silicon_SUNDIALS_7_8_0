//! tides_spin_migration.rs — Rust twin of
//! `porttest/tides_spin_migration_c.c`, itself
//! `reboundx/examples/tides_spin_migration_driven_obliquity_tides/problem.c`.
//!
//! Two Earth-like planets migrating in a disk while tides evolve their
//! spins (Millholland & Laughlin 2019 parameters). Exercises TWO
//! simultaneous REBOUNDx forces (`tides_spin` + `modify_orbits_forces`)
//! and mid-run parameter mutation: migration is switched off halfway by
//! setting `tau_a` to infinity.
//!
//! Writes `state_migration_rust.txt` as raw IEEE-754 bit patterns.
//!
//! Part of reboundx_rs, GPL-3.0-or-later. Based on REBOUNDx
//! (c) Dan Tamayo, Hanno Rein et al.
#![allow(non_snake_case)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. Each waiver below carries its own reason; the same
// list and the rationale are in README.md under "Building and testing".
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
use reboundx_rs::*;
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let tmax: f64 = if args.len() > 1 {
        args[1].parse().unwrap_or(100. * 2. * std::f64::consts::PI)
    } else { 100. * 2. * std::f64::consts::PI };

    let mut sim = reb_simulation_create();

    let solar_mass = 1.;
    let solar_rad = 0.00465;
    reb_simulation_add_fmt(&mut sim, "m r",
        &[reb_fmt_arg::d(solar_mass), reb_fmt_arg::d(solar_rad)]);

    let p1_mass = 5. * 3.0e-6;
    let p1_rad = 2.5 * 4.26e-5;
    reb_simulation_add_fmt(&mut sim, "m a e r inc Omega pomega M", &[
        reb_fmt_arg::d(p1_mass), reb_fmt_arg::d(0.17308688), reb_fmt_arg::d(0.01),
        reb_fmt_arg::d(p1_rad), reb_fmt_arg::d(0.5 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.))]);

    let p2_mass = 5. * 3.0e-6;
    let p2_rad = 2.5 * 4.26e-5;
    reb_simulation_add_fmt(&mut sim, "m a e r inc Omega pomega M", &[
        reb_fmt_arg::d(p2_mass), reb_fmt_arg::d(0.23290608), reb_fmt_arg::d(0.01),
        reb_fmt_arg::d(p2_rad), reb_fmt_arg::d(-0.431 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.)),
        reb_fmt_arg::d(0.0 * (std::f64::consts::PI / 180.))]);

    sim.N_active = 3;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = 1e-3;

    rebx_attach(&mut sim);
    let effect = rebx_load_force(&mut sim, "tides_spin").expect("tides_spin");
    rebx_add_force(&mut sim, effect);

    let solar_spin_period = 20. * 2. * std::f64::consts::PI / 365.;
    let solar_spin = (2. * std::f64::consts::PI) / solar_spin_period;
    let solar_Q = 1000000.;
    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    let orb2 = reb_orbit_from_particle(sim.G, sim.particles[2], sim.particles[0]);

    let spin_period_1 = 5. * 2. * std::f64::consts::PI / 365.;
    let spin_1 = (2. * std::f64::consts::PI) / spin_period_1;
    let planet_Q = 10000.;
    let spin_period_2 = 3. * 2. * std::f64::consts::PI / 365.;
    let spin_2 = (2. * std::f64::consts::PI) / spin_period_2;

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "k2", 0.07);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "I",
                              0.07 * solar_mass * solar_rad * solar_rad);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(0), "Omega",
                             reb_vec3d { x: 0., y: 0., z: solar_spin });
        rebx_set_param_double(rebx, rebx_ap::particle(0), "tau",
                              1. / (2. * orb.n * solar_Q));

        rebx_set_param_double(rebx, rebx_ap::particle(1), "k2", 0.4);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "I",
                              0.25 * p1_mass * p1_rad * p1_rad);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(1), "Omega",
                             reb_vec3d { x: 0., y: spin_1 * -0.0261769, z: spin_1 * 0.99965732 });
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau",
                              1. / (2. * orb.n * planet_Q));

        rebx_set_param_double(rebx, rebx_ap::particle(2), "k2", 0.4);
        rebx_set_param_double(rebx, rebx_ap::particle(2), "I",
                              0.25 * p2_mass * p2_rad * p2_rad);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(2), "Omega",
                             reb_vec3d { x: 0., y: spin_2 * 0.0249736, z: spin_2 * 0.99968811 });
        rebx_set_param_double(rebx, rebx_ap::particle(2), "tau",
                              1. / (2. * orb2.n * planet_Q));
    }

    let mo = rebx_load_force(&mut sim, "modify_orbits_forces").expect("modify_orbits_forces");
    rebx_add_force(&mut sim, mo);

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_a",
                              -5e6 * 2. * std::f64::consts::PI);
        rebx_set_param_double(rebx, rebx_ap::particle(2), "tau_a",
                              (-5e6 * 2. * std::f64::consts::PI) / 1.1);
    }

    reb_simulation_move_to_com(&mut sim);

    let L_orb = reb_simulation_angular_momentum(&sim);
    let L_spin = rebx_with(&mut sim, |sim, rebx| {
        rebx_tools_spin_angular_momentum(sim, rebx)
    }).expect("extras attached");
    let newz = reb_vec3d_add(L_orb, L_spin);
    let newx = reb_vec3d_cross(reb_vec3d { x: 0., y: 0., z: 1. }, newz);
    let rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_with(&mut sim, |sim, rebx| { rebx_simulation_irotate(sim, rebx, rot); });
    rebx_with(&mut sim, |sim, rebx| { rebx_spin_initialize_ode(sim, rebx, effect); });

    reb_simulation_integrate(&mut sim, tmax / 2.);

    // Migration switching off
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_a", f64::INFINITY);
        rebx_set_param_double(rebx, rebx_ap::particle(2), "tau_a", f64::INFINITY);
    }

    reb_simulation_integrate(&mut sim, tmax);

    let mut f = std::fs::File::create("state_migration_rust.txt").unwrap();
    writeln!(f, "example tides_spin_migration_driven_obliquity_tides tmax {:016x}", bits(tmax)).unwrap();
    writeln!(f, "t {:016x}", bits(sim.t)).unwrap();
    writeln!(f, "dt {:016x}", bits(sim.dt)).unwrap();
    writeln!(f, "N {}", sim.N).unwrap();
    for i in 0..sim.N {
        let p = sim.particles[i];
        writeln!(f, "p {} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
                 i, bits(p.x), bits(p.y), bits(p.z),
                 bits(p.vx), bits(p.vy), bits(p.vz), bits(p.m)).unwrap();
        let om = rebx_extras_ref(&sim)
            .and_then(|rebx| rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega"));
        match om {
            Some(o) => writeln!(f, "Omega {} {:016x} {:016x} {:016x}",
                                i, bits(o.x), bits(o.y), bits(o.z)).unwrap(),
            None => writeln!(f, "Omega {} none", i).unwrap(),
        }
    }
    println!("migration done t={:.17e}", sim.t);
}
