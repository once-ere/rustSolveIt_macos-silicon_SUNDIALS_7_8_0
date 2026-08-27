//! derivatives_test.rs — Rust twin of porttest/derivatives_test.c.
//! Calls all 65 reb_particle_derivative_* functions for two fixed
//! configurations and dumps the raw IEEE-754 bit patterns of the
//! resulting particle (x,y,z,vx,vy,vz,m) to derivatives_rust.txt.
//! Part of the rebound_rs port verification. GPL-3.0-or-later.

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

type DerivFn = fn(f64, reb_particle, reb_particle) -> reb_particle;

fn main() {
    let fns: [(&str, DerivFn); 65] = [
        ("lambda", reb_particle_derivative_lambda),
        ("h", reb_particle_derivative_h),
        ("k", reb_particle_derivative_k),
        ("k_k", reb_particle_derivative_k_k),
        ("h_h", reb_particle_derivative_h_h),
        ("lambda_lambda", reb_particle_derivative_lambda_lambda),
        ("k_lambda", reb_particle_derivative_k_lambda),
        ("h_lambda", reb_particle_derivative_h_lambda),
        ("k_h", reb_particle_derivative_k_h),
        ("a", reb_particle_derivative_a),
        ("a_a", reb_particle_derivative_a_a),
        ("ix", reb_particle_derivative_ix),
        ("ix_ix", reb_particle_derivative_ix_ix),
        ("iy", reb_particle_derivative_iy),
        ("iy_iy", reb_particle_derivative_iy_iy),
        ("k_ix", reb_particle_derivative_k_ix),
        ("h_ix", reb_particle_derivative_h_ix),
        ("lambda_ix", reb_particle_derivative_lambda_ix),
        ("lambda_iy", reb_particle_derivative_lambda_iy),
        ("h_iy", reb_particle_derivative_h_iy),
        ("k_iy", reb_particle_derivative_k_iy),
        ("ix_iy", reb_particle_derivative_ix_iy),
        ("a_ix", reb_particle_derivative_a_ix),
        ("a_iy", reb_particle_derivative_a_iy),
        ("a_lambda", reb_particle_derivative_a_lambda),
        ("a_h", reb_particle_derivative_a_h),
        ("a_k", reb_particle_derivative_a_k),
        ("m", reb_particle_derivative_m),
        ("m_a", reb_particle_derivative_m_a),
        ("m_lambda", reb_particle_derivative_m_lambda),
        ("m_h", reb_particle_derivative_m_h),
        ("m_k", reb_particle_derivative_m_k),
        ("m_ix", reb_particle_derivative_m_ix),
        ("m_iy", reb_particle_derivative_m_iy),
        ("m_m", reb_particle_derivative_m_m),
        ("e", reb_particle_derivative_e),
        ("e_e", reb_particle_derivative_e_e),
        ("inc", reb_particle_derivative_inc),
        ("inc_inc", reb_particle_derivative_inc_inc),
        ("Omega", reb_particle_derivative_Omega),
        ("Omega_Omega", reb_particle_derivative_Omega_Omega),
        ("omega", reb_particle_derivative_omega),
        ("omega_omega", reb_particle_derivative_omega_omega),
        ("f", reb_particle_derivative_f),
        ("f_f", reb_particle_derivative_f_f),
        ("a_e", reb_particle_derivative_a_e),
        ("a_inc", reb_particle_derivative_a_inc),
        ("a_Omega", reb_particle_derivative_a_Omega),
        ("a_omega", reb_particle_derivative_a_omega),
        ("a_f", reb_particle_derivative_a_f),
        ("e_inc", reb_particle_derivative_e_inc),
        ("e_Omega", reb_particle_derivative_e_Omega),
        ("e_omega", reb_particle_derivative_e_omega),
        ("e_f", reb_particle_derivative_e_f),
        ("m_e", reb_particle_derivative_m_e),
        ("inc_Omega", reb_particle_derivative_inc_Omega),
        ("inc_omega", reb_particle_derivative_inc_omega),
        ("inc_f", reb_particle_derivative_inc_f),
        ("m_inc", reb_particle_derivative_m_inc),
        ("omega_Omega", reb_particle_derivative_omega_Omega),
        ("Omega_f", reb_particle_derivative_Omega_f),
        ("m_Omega", reb_particle_derivative_m_Omega),
        ("omega_f", reb_particle_derivative_omega_f),
        ("m_omega", reb_particle_derivative_m_omega),
        ("m_f", reb_particle_derivative_m_f),
    ];

    let G = 1.0;

    let mut primary1 = reb_particle::default();
    primary1.m = 1.0;
    primary1.x = 0.1;
    primary1.y = -0.2;
    primary1.z = 0.05;
    primary1.vx = 0.01;
    primary1.vy = -0.03;
    primary1.vz = 0.002;

    let mut po1 = reb_particle::default();
    po1.m = 1e-3;
    po1.x = 1.3;
    po1.y = 0.4;
    po1.z = 0.1;
    po1.vx = -0.2;
    po1.vy = 0.9;
    po1.vz = 0.03;

    let mut primary2 = reb_particle::default();
    primary2.m = 2.3;

    let mut po2 = reb_particle::default();
    po2.m = 1e-5;
    po2.x = 0.7;
    po2.y = -0.5;
    po2.z = 0.2;
    po2.vx = 0.4;
    po2.vy = 1.1;
    po2.vz = -0.05;

    let mut out = std::fs::File::create("derivatives_rust.txt").expect("create derivatives_rust.txt");
    for cfg in 1..=2 {
        let (primary, po) = if cfg == 1 { (primary1, po1) } else { (primary2, po2) };
        for (name, f) in fns.iter() {
            let np = f(G, primary, po);
            writeln!(
                out,
                "{} cfg{} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x} {:016x}",
                name,
                cfg,
                np.x.to_bits(),
                np.y.to_bits(),
                np.z.to_bits(),
                np.vx.to_bits(),
                np.vy.to_bits(),
                np.vz.to_bits(),
                np.m.to_bits()
            )
            .expect("write");
        }
    }
    println!("derivatives_rust.txt written ({} functions x 2 configs)", fns.len());
}
