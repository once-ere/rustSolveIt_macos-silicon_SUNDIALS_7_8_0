//! tides_spin_kozai.rs — Rust twin of `porttest/tides_spin_kozai_c.c`,
//! itself `reboundx/examples/tides_spin_kozai/problem.c`.
//!
//! A Lidov-Kozai cycle: a planet perturbed by a distant stellar
//! companion, with tides and spin. Uses the ADAPTIVE IAS15 integrator,
//! so this also checks that the REBOUNDx force feeds IAS15's error
//! estimator identically — the whole adaptive timestep sequence has to
//! match bit-for-bit. Also exercises a second force (`gr_potential`).
//!
//! Writes `state_kozai_rust.txt` as raw IEEE-754 bit patterns.
//!
//! Part of reboundx_rs, GPL-3.0-or-later. Based on REBOUNDx
//! (c) Dan Tamayo, Hanno Rein et al.
#![allow(non_snake_case)]
// Mirrors the C's `struct reb_particle star = {0}; star.m = 1.;`.
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::field_reassign_with_default)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. Each waiver below carries its own reason; the same
// list and the rationale are in README.md under "Building and testing".
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::assign_op_pattern)]
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
    let tmax: f64 = if args.len() > 1 { args[1].parse().unwrap_or(1e5) } else { 1e5 };

    let mut sim = reb_simulation_create();
    sim.dt = std::f64::consts::PI * 1e-1;
    reb_simulation_set_integrator(&mut sim, "ias15");

    let mut star = reb_particle::default();
    star.m = 1.;
    star.r = 0.00465;
    reb_simulation_add(&mut sim, star);

    let planet_m = 0.054 * 9.55e-4;
    let planet_r = 0.3 * 4.676e-4;
    let planet_a = 2.;
    let planet_e = 0.001;
    reb_simulation_add_fmt(&mut sim, "m r a e", &[
        reb_fmt_arg::d(planet_m), reb_fmt_arg::d(planet_r),
        reb_fmt_arg::d(planet_a), reb_fmt_arg::d(planet_e)]);

    let perturber_m = 1.;
    let perturber_a = 50.;
    let perturber_e = 0.7 * std::f64::consts::PI / 180.;
    let perturber_inc = 80. * std::f64::consts::PI / 180.;
    reb_simulation_add_fmt(&mut sim, "m a e inc", &[
        reb_fmt_arg::d(perturber_m), reb_fmt_arg::d(perturber_a),
        reb_fmt_arg::d(perturber_e), reb_fmt_arg::d(perturber_inc)]);

    rebx_attach(&mut sim);
    let effect = rebx_load_force(&mut sim, "tides_spin").expect("tides_spin");
    rebx_add_force(&mut sim, effect);

    let solar_k2 = 0.07;
    let solar_spin_period = 4.6 * 2. * std::f64::consts::PI / 365.;
    let solar_spin = (2. * std::f64::consts::PI) / solar_spin_period;
    let solar_Q = 1e6;
    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    let solar_tau = 1. / (2. * solar_Q * orb.n);

    let planet_k2 = 0.4;
    let spin_period_p = 1. * 2. * std::f64::consts::PI / 365.;
    let spin_p = (2. * std::f64::consts::PI) / spin_period_p;
    let theta_p = 0. * std::f64::consts::PI / 180.;
    let phi_p = 0. * std::f64::consts::PI / 180.;
    let Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    let planet_Q = 3e5;

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "k2", solar_k2);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "I", 0.07 * star.m * star.r * star.r);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(0), "Omega",
                             reb_vec3d { x: 0., y: 0., z: solar_spin });
        rebx_set_param_double(rebx, rebx_ap::particle(0), "tau", solar_tau);

        rebx_set_param_double(rebx, rebx_ap::particle(1), "k2", planet_k2);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "I",
                              0.25 * planet_m * planet_r * planet_r);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(1), "Omega", Omega_sv);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau",
                              1. / (2. * planet_Q * orb.n));
    }

    let gr = rebx_load_force(&mut sim, "gr_potential").expect("gr_potential");
    rebx_add_force(&mut sim, gr);
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::force(gr), "c", 10065.32);
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

    reb_simulation_integrate(&mut sim, tmax);

    let mut f = std::fs::File::create("state_kozai_rust.txt").unwrap();
    writeln!(f, "example tides_spin_kozai tmax {:016x}", bits(tmax)).unwrap();
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
    println!("kozai done t={:.17e}", sim.t);
}
