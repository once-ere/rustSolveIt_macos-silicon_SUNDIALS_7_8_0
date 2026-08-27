//! rebx_binary_roundtrip.rs — self-check for the REBOUNDx binary
//! format (src/output.rs and src/input.rs).
//!
//! Builds a simulation with three particles, attaches REBOUNDx, adds
//! two forces and two operators with parameters of several types, and
//! writes the state with `rebx_output_binary`. It then reads the file
//! back into a *fresh* simulation with `rebx_create_extras_from_binary`
//! and compares everything that came back against the original —
//! floating-point values by their IEEE-754 bit patterns, so a value
//! that survived the round trip only approximately would still fail.
//!
//! Prints PASS/FAIL per item and exits non-zero if anything failed.
//!
//! Given a filename argument, it reads *that* file back instead of the
//! one it just wrote, which is how the reader is checked against a
//! binary produced by the C `libreboundx` (build a C program that makes
//! the same state, call the C's `rebx_output_binary`, then pass its
//! file here). The file this example writes is still produced, and is
//! still what the final re-serialization check compares against.
//!
//! Note on types: REBOUNDx 5.1.0's `rebx_sizeof` (core.c) has no case
//! for `REBX_TYPE_STRING`, `REBX_TYPE_ORBIT` or `REBX_TYPE_UINT32`, so
//! the C writes a zero-length payload for parameters of those types and
//! they cannot round trip. This translation reproduces that byte
//! format exactly, so those types are deliberately not exercised here.
//! `REBX_TYPE_POINTER` and `REBX_TYPE_ODE` are skipped by the C writer
//! by design (a pointer means nothing in a later process).
//!
//! Part of reboundx_rs, GPL-3.0-or-later. Based on REBOUNDx
//! (c) Dan Tamayo, Hanno Rein et al.
#![allow(non_snake_case)]
// Mirrors the C's `struct reb_particle p = {0}; p.m = ..;`.
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

/// Three particles: a star and two planets. Only the parameter lists
/// matter to the binary format, but the particle count has to match for
/// the PARTICLE_INDEX fields to land on real particles.
fn build_sim() -> reb_simulation {
    let mut sim = reb_simulation_create();

    let mut star = reb_particle::default();
    star.m = 1.;
    star.r = 0.005;
    reb_simulation_add(&mut sim, star);

    let mut p1 = reb_particle::default();
    p1.m = 1e-3;
    p1.x = 1.;
    p1.vy = 1.;
    reb_simulation_add(&mut sim, p1);

    let mut p2 = reb_particle::default();
    p2.m = 5e-4;
    p2.x = 2.;
    p2.vy = 0.65;
    reb_simulation_add(&mut sim, p2);

    sim
}

struct Report {
    failures: usize,
    checks: usize,
}

impl Report {
    fn check(&mut self, label: &str, ok: bool) {
        self.checks += 1;
        if ok {
            println!("PASS  {}", label);
        } else {
            self.failures += 1;
            println!("FAIL  {}", label);
        }
    }
}

/// Bit-exact comparison of two optional doubles. A value that is
/// missing on either side fails.
fn eq_f64(a: Option<f64>, b: Option<f64>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.to_bits() == b.to_bits(),
        _ => false,
    }
}

fn eq_vec3d(a: Option<reb_vec3d>, b: Option<reb_vec3d>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => {
            a.x.to_bits() == b.x.to_bits()
                && a.y.to_bits() == b.y.to_bits()
                && a.z.to_bits() == b.z.to_bits()
        }
        _ => false,
    }
}

fn eq_i32(a: Option<i32>, b: Option<i32>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a == b,
        _ => false,
    }
}

fn main() {
    let filename = "rebx_binary_roundtrip.bin";
    // Optional: read back a binary written by something else (the C
    // libreboundx) rather than the one written here.
    let args: Vec<String> = std::env::args().collect();
    let read_from: String = if args.len() > 1 {
        args[1].clone()
    } else {
        filename.to_string()
    };

    /*****************************************************
     Build the original state
     *****************************************************/

    let mut sim = build_sim();
    rebx_attach(&mut sim);

    // Two forces, so that the order of allocated_forces and of
    // additional_forces is actually observable.
    let gr = rebx_load_force(&mut sim, "gr_potential").expect("gr_potential");
    rebx_add_force(&mut sim, gr);
    let cf = rebx_load_force(&mut sim, "central_force").expect("central_force");
    rebx_add_force(&mut sim, cf);

    // Two operators, with a pre-timestep step, a post-timestep step and
    // a second post-timestep step, so both step lists are non-trivial.
    let mm = rebx_load_operator(&mut sim, "modify_mass").expect("modify_mass");
    let dr = rebx_load_operator(&mut sim, "drift").expect("drift");
    rebx_add_operator_step(&mut sim, mm, 0.5, rebx_timing::REBX_TIMING_PRE);
    rebx_add_operator_step(&mut sim, mm, 0.5, rebx_timing::REBX_TIMING_POST);
    rebx_add_operator_step(&mut sim, dr, 1.0, rebx_timing::REBX_TIMING_POST);

    // Parameter values, chosen so that every one has a distinctive bit
    // pattern (nothing that a zeroed buffer could accidentally match).
    let c_gr = 10065.320005560323;
    let src_gr = 7_i32;
    let Acentral = 1.2345678901234567e-8;
    let c_op = -3.14159265358979e12;
    let tau_mass = 1.7976931348623157e30;
    let primary = -12345_i32;
    let Omega = reb_vec3d {
        x: 1.5e-7,
        y: -2.5e13,
        z: 0.30000000000000004,
    };
    let beta = 0.1 + 0.2; // 0.30000000000000004, not 0.3
    let coordinates = 2_i32;

    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        // Force parameters (two on one force, to check list order).
        rebx_set_param_double(rebx, rebx_ap::force(gr), "c", c_gr);
        rebx_set_param_int(rebx, rebx_ap::force(gr), "gr_source", src_gr);
        rebx_set_param_double(rebx, rebx_ap::force(cf), "Acentral", Acentral);

        // Operator parameter.
        rebx_set_param_double(rebx, rebx_ap::operator_(mm), "c", c_op);

        // Per-particle parameters of several types.
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_mass", tau_mass);
        rebx_set_param_int(rebx, rebx_ap::particle(1), "primary", primary);
        rebx_set_param_vec3d(rebx, rebx_ap::particle(1), "Omega", Omega);
        // REBX_TYPE_FORCE: stored on disk as the force's name.
        rebx_set_param_force(rebx, rebx_ap::particle(1), "force", gr);

        rebx_set_param_double(rebx, rebx_ap::particle(2), "beta", beta);
        rebx_set_param_int(rebx, rebx_ap::particle(2), "coordinates", coordinates);
    }

    // Snapshot of the original lists for the order comparisons below.
    let (reg0, forces0, ops0, add0, pre0, post0) = {
        let rebx = rebx_extras_ref(&sim).expect("extras");
        let reg: Vec<(String, rebx_param_type)> = rebx
            .registered_params
            .iter()
            .map(|p| (p.name.clone(), p.type_))
            .collect();
        let forces: Vec<String> = rebx.allocated_forces.iter().map(|f| f.name.clone()).collect();
        let ops: Vec<String> = rebx
            .allocated_operators
            .iter()
            .map(|o| o.name.clone())
            .collect();
        let add: Vec<String> = rebx
            .additional_forces
            .iter()
            .map(|i| rebx.allocated_forces[*i].name.clone())
            .collect();
        let pre: Vec<(String, u64)> = rebx
            .pre_timestep_modifications
            .iter()
            .map(|s| {
                (
                    rebx.allocated_operators[s.operator_].name.clone(),
                    s.dt_fraction.to_bits(),
                )
            })
            .collect();
        let post: Vec<(String, u64)> = rebx
            .post_timestep_modifications
            .iter()
            .map(|s| {
                (
                    rebx.allocated_operators[s.operator_].name.clone(),
                    s.dt_fraction.to_bits(),
                )
            })
            .collect();
        (reg, forces, ops, add, pre, post)
    };

    /*****************************************************
     Write, then read back into a fresh simulation
     *****************************************************/

    rebx_output_binary(&mut sim, filename);
    let written = match std::fs::metadata(filename) {
        Ok(m) => m.len(),
        Err(e) => {
            println!("FAIL  rebx_output_binary did not produce a file: {}", e);
            std::process::exit(1);
        }
    };
    println!("wrote {} ({} bytes)", filename, written);
    if read_from != filename {
        println!("reading back {} instead (foreign binary)", read_from);
    }

    let mut sim2 = build_sim();
    rebx_create_extras_from_binary(&mut sim2, &read_from);

    /*****************************************************
     Compare
     *****************************************************/

    let mut r = Report {
        failures: 0,
        checks: 0,
    };

    let rebx2 = match rebx_extras_ref(&sim2) {
        Some(rebx2) => rebx2,
        None => {
            println!("FAIL  no REBOUNDx state attached after rebx_create_extras_from_binary");
            std::process::exit(1);
        }
    };

    // --- structure -------------------------------------------------
    let reg1: Vec<(String, rebx_param_type)> = rebx2
        .registered_params
        .iter()
        .map(|p| (p.name.clone(), p.type_))
        .collect();
    r.check(
        &format!("registered_params ({} entries, names+types+order)", reg0.len()),
        reg0 == reg1,
    );

    let forces1: Vec<String> = rebx2.allocated_forces.iter().map(|f| f.name.clone()).collect();
    r.check(
        &format!("allocated_forces {:?}", forces0),
        forces0 == forces1,
    );

    let ops1: Vec<String> = rebx2
        .allocated_operators
        .iter()
        .map(|o| o.name.clone())
        .collect();
    r.check(&format!("allocated_operators {:?}", ops0), ops0 == ops1);

    let add1: Vec<String> = rebx2
        .additional_forces
        .iter()
        .map(|i| rebx2.allocated_forces[*i].name.clone())
        .collect();
    r.check(
        &format!("additional_forces {:?} (order is load-bearing)", add0),
        add0 == add1,
    );

    let pre1: Vec<(String, u64)> = rebx2
        .pre_timestep_modifications
        .iter()
        .map(|s| {
            (
                rebx2.allocated_operators[s.operator_].name.clone(),
                s.dt_fraction.to_bits(),
            )
        })
        .collect();
    r.check(
        &format!("pre_timestep_modifications {:?}", pre0),
        pre0 == pre1,
    );

    let post1: Vec<(String, u64)> = rebx2
        .post_timestep_modifications
        .iter()
        .map(|s| {
            (
                rebx2.allocated_operators[s.operator_].name.clone(),
                s.dt_fraction.to_bits(),
            )
        })
        .collect();
    r.check(
        &format!("post_timestep_modifications {:?}", post0),
        post0 == post1,
    );

    // The re-created forces/operators are at the same indices here, but
    // look them up by name the way a user would.
    let gr2 = rebx_get_force(rebx2, "gr_potential");
    let cf2 = rebx_get_force(rebx2, "central_force");
    let mm2 = rebx_get_operator(rebx2, "modify_mass");
    r.check("force 'gr_potential' present", gr2.is_some());
    r.check("force 'central_force' present", cf2.is_some());
    r.check("operator 'modify_mass' present", mm2.is_some());
    let gr2 = gr2.unwrap_or(0);
    let cf2 = cf2.unwrap_or(0);
    let mm2 = mm2.unwrap_or(0);

    r.check(
        "force 'gr_potential' update_accelerations restored",
        rebx2.allocated_forces[gr2].update_accelerations.is_some(),
    );
    r.check(
        "operator 'modify_mass' step_function restored",
        rebx2.allocated_operators[mm2].step_function.is_some(),
    );

    // --- force / operator parameters -------------------------------
    r.check(
        &format!("force gr_potential 'c' double bits {:016x}", c_gr.to_bits()),
        eq_f64(
            Some(c_gr),
            rebx_get_param_double(rebx2, rebx_ap::force(gr2), "c"),
        ),
    );
    r.check(
        &format!("force gr_potential 'gr_source' int {}", src_gr),
        eq_i32(
            Some(src_gr),
            rebx_get_param_int(rebx2, rebx_ap::force(gr2), "gr_source"),
        ),
    );
    r.check(
        &format!(
            "force central_force 'Acentral' double bits {:016x}",
            Acentral.to_bits()
        ),
        eq_f64(
            Some(Acentral),
            rebx_get_param_double(rebx2, rebx_ap::force(cf2), "Acentral"),
        ),
    );
    r.check(
        &format!("operator modify_mass 'c' double bits {:016x}", c_op.to_bits()),
        eq_f64(
            Some(c_op),
            rebx_get_param_double(rebx2, rebx_ap::operator_(mm2), "c"),
        ),
    );

    // The parameter list order must survive too (the C writes tail to
    // head and the reader prepends).
    let gr_names0: Vec<String> = rebx_extras_ref(&sim)
        .expect("extras")
        .ap(rebx_ap::force(gr))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    let gr_names1: Vec<String> = rebx2
        .ap(rebx_ap::force(gr2))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    r.check(
        &format!("force gr_potential param order {:?}", gr_names0),
        gr_names0 == gr_names1,
    );

    // --- particle parameters ---------------------------------------
    r.check(
        &format!(
            "particle 1 'tau_mass' double bits {:016x}",
            tau_mass.to_bits()
        ),
        eq_f64(
            Some(tau_mass),
            rebx_get_param_double(rebx2, rebx_ap::particle(1), "tau_mass"),
        ),
    );
    r.check(
        &format!("particle 1 'primary' int {}", primary),
        eq_i32(
            Some(primary),
            rebx_get_param_int(rebx2, rebx_ap::particle(1), "primary"),
        ),
    );
    r.check(
        &format!(
            "particle 1 'Omega' vec3d bits {:016x} {:016x} {:016x}",
            Omega.x.to_bits(),
            Omega.y.to_bits(),
            Omega.z.to_bits()
        ),
        eq_vec3d(
            Some(Omega),
            rebx_get_param_vec3d(rebx2, rebx_ap::particle(1), "Omega"),
        ),
    );
    r.check(
        "particle 1 'force' REBX_TYPE_FORCE relinked to gr_potential",
        match rebx_get_param_force(rebx2, rebx_ap::particle(1), "force") {
            Some(idx) => rebx2.allocated_forces[idx].name == "gr_potential",
            None => false,
        },
    );
    r.check(
        &format!("particle 2 'beta' double bits {:016x}", beta.to_bits()),
        eq_f64(
            Some(beta),
            rebx_get_param_double(rebx2, rebx_ap::particle(2), "beta"),
        ),
    );
    r.check(
        &format!("particle 2 'coordinates' int {}", coordinates),
        eq_i32(
            Some(coordinates),
            rebx_get_param_int(rebx2, rebx_ap::particle(2), "coordinates"),
        ),
    );
    r.check(
        "particle 0 has no parameters (as written)",
        rebx2.ap(rebx_ap::particle(0)).is_empty(),
    );

    let p1_names0: Vec<String> = rebx_extras_ref(&sim)
        .expect("extras")
        .ap(rebx_ap::particle(1))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    let p1_names1: Vec<String> = rebx2
        .ap(rebx_ap::particle(1))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    r.check(
        &format!("particle 1 param order {:?}", p1_names0),
        p1_names0 == p1_names1,
    );

    // --- writing the reloaded state must reproduce the same bytes ---
    let bytes0 = std::fs::read(filename).unwrap_or_default();
    let filename2 = "rebx_binary_roundtrip2.bin";
    rebx_output_binary(&mut sim2, filename2);
    let bytes1 = std::fs::read(filename2).unwrap_or_default();
    r.check(
        &format!("re-serializing the reloaded state is byte-identical ({} bytes)", bytes0.len()),
        !bytes0.is_empty() && bytes0 == bytes1,
    );

    println!(
        "\n{} / {} checks passed",
        r.checks - r.failures,
        r.checks
    );
    if r.failures > 0 {
        println!("ROUND TRIP FAILED");
        std::process::exit(1);
    }
    println!("ROUND TRIP PASSED");
}
