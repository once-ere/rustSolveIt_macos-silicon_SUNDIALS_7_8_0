//! integrators_test.rs — Rust twin of porttest/integrators_test.c.
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

use rebound_rs::*;
use std::io::Write;

fn bits(x: f64) -> u64 {
    x.to_bits()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let integrator = if args.len() > 1 { args[1].clone() } else { "ias15".to_string() };
    let order: u32 = if args.len() > 2 { args[2].parse().unwrap_or(2) } else { 2 };
    let nsteps: usize = if args.len() > 3 { args[3].parse().unwrap_or(1000) } else { 1000 };

    let mut sim = reb_simulation_create();
    let r = &mut sim;
    // whfast configurations are encoded as pseudo names; see the C twin.
    let real_integrator = if integrator.starts_with("whfast") {
        "whfast"
    } else if integrator.starts_with("saba") {
        "saba"
    } else if integrator.starts_with("janus") {
        "janus"
    } else if integrator.starts_with("eos") {
        "eos"
    } else if integrator.starts_with("mercurius") {
        "mercurius"
    } else if integrator.starts_with("bs") {
        "bs"
    } else if integrator.starts_with("trace") {
        "trace"
    } else {
        integrator.as_str()
    };
    reb_simulation_set_integrator(r, real_integrator);
    if integrator == "leapfrog" {
        if let reb_integrator_state::leapfrog(ref mut lf) = r.integrator {
            lf.order = order;
        }
    }
    if integrator.starts_with("whfast") {
        if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
            match integrator.as_str() {
                "whfast-c11" => wh.corrector = 11,
                "whfast-c17" => {
                    wh.corrector = 17;
                    wh.corrector2 = 1;
                }
                "whfast-dh" => {
                    wh.coordinates =
                        rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC
                }
                "whfast-whds" => {
                    wh.coordinates =
                        rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_WHDS
                }
                "whfast-bary" => {
                    wh.coordinates =
                        rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC
                }
                "whfast-mk" => {
                    wh.kernel =
                        rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK
                }
                "whfast-comp" => {
                    wh.kernel =
                        rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION
                }
                "whfast-lazy" => {
                    wh.kernel = rebound_rs::integrator_whfast::REB_INTEGRATOR_WHFAST_KERNEL_LAZY
                }
                "whfast-usafe" => wh.safe_mode = 0,
                _ => {}
            }
        }
    }
    if integrator.starts_with("saba") {
        if let reb_integrator_state::saba(ref mut sb) = r.integrator {
            use rebound_rs::integrator_saba::*;
            match integrator.as_str() {
                "saba-1" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_1,
                "saba-2" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_2,
                "saba-3" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_3,
                "saba-4" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_4,
                "saba-cm2" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_CM_2,
                "saba-cl2" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_CL_2,
                "saba-104" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_10_4,
                "saba-864" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_8_6_4,
                "saba-h844" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_H_8_4_4,
                "saba-h864" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_H_8_6_4,
                "saba-h1064" => sb.type_ = REB_INTEGRATOR_SABA_TYPE_H_10_6_4,
                "saba-usafe" => sb.safe_mode = 0,
                _ => {}
            }
        }
    }
    if integrator.starts_with("janus") {
        if let reb_integrator_state::janus(ref mut jn) = r.integrator {
            match integrator.as_str() {
                "janus-2" => jn.order = 2,
                "janus-4" => jn.order = 4,
                "janus-8" => jn.order = 8,
                "janus-10" => jn.order = 10,
                _ => {}
            }
        }
    }
    if integrator.starts_with("eos") {
        if let reb_integrator_state::eos(ref mut es) = r.integrator {
            // eos-<phi0>-<phi1> with numeric type ids 0-8, and eos-usafe
            if integrator == "eos-usafe" {
                es.safe_mode = 0;
            } else if integrator.len() == 7 {
                let b = integrator.as_bytes();
                if b[3] == b'-' && b[5] == b'-' {
                    es.phi0 = (b[4] - b'0') as i32;
                    es.phi1 = (b[6] - b'0') as i32;
                }
            }
        }
    }
    if integrator.starts_with("mercurius") {
        if let reb_integrator_state::mercurius(ref mut mc) = r.integrator {
            use rebound_rs::integrator_mercurius::*;
            // see the C twin for the pseudo-name catalogue
            match integrator.as_str() {
                "mercurius-usafe" => mc.safe_mode = 0,
                "mercurius-c4" => mc.L = reb_integrator_mercurius_L_C4,
                "mercurius-c5" => mc.L = reb_integrator_mercurius_L_C5,
                "mercurius-inf" => mc.L = reb_integrator_mercurius_L_infinity,
                "mercurius-hill01" => mc.r_crit_hill = 0.1,
                _ => {}
            }
        }
    }
    if integrator.starts_with("bs") {
        if let reb_integrator_state::bs(ref mut b) = r.integrator {
            // see the C twin for the pseudo-name catalogue
            match integrator.as_str() {
                "bs-tight" => {
                    b.eps_abs = 1e-11;
                    b.eps_rel = 1e-11;
                }
                "bs-loose" => {
                    b.eps_abs = 1e-6;
                    b.eps_rel = 1e-6;
                }
                "bs-maxdt" => b.max_dt = 0.02,
                _ => {}
            }
        }
    }
    if integrator.starts_with("trace") {
        if let reb_integrator_state::trace(ref mut tr) = r.integrator {
            use rebound_rs::integrator_trace::*;
            // see the C twin for the pseudo-name catalogue
            match integrator.as_str() {
                "trace-pbs" => tr.peri_mode = REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS,
                "trace-ias15" => tr.peri_mode = REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15,
                "trace-hill1" => tr.r_crit_hill = 1.,
                "trace-perinone" => tr.S_peri = Some(reb_integrator_trace_switch_peri_none),
                "trace-eta001" => tr.peri_crit_eta = 0.01,
                _ => {}
            }
        }
    }
    r.G = 1.0;
    r.dt = 0.01;

    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(r, star);

    let mut planet = reb_particle::default();
    planet.m = 1e-3;
    planet.x = 1.6; // apocenter of a=1, e=0.6 orbit
    planet.vy = 0.5; // roughly the apocenter speed
    reb_simulation_add(r, planet);

    let mut moon = reb_particle::default();
    moon.m = 1e-7;
    moon.x = 1.7;
    moon.vy = 0.6;
    moon.z = 0.01;
    moon.vz = 0.001;
    reb_simulation_add(r, moon);

    reb_simulation_steps(r, nsteps);

    let mut f = std::fs::File::create("state_rust_final.txt").unwrap();
    writeln!(f, "integrator {} order {} steps {}", integrator, order, nsteps).unwrap();
    writeln!(f, "t {:016x}", bits(r.t)).unwrap();
    writeln!(f, "dt {:016x}", bits(r.dt)).unwrap();
    writeln!(f, "steps_done {}", r.steps_done).unwrap();
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
    drop(f);
    println!("{} done: t={:e} steps={}", integrator, r.t, r.steps_done);
}
