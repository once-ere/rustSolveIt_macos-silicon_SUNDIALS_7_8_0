//! Integration tests for the transforms_rotations_derivatives module group of rebound_rs.
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

// ===========================================================================
// Shared helpers
// ===========================================================================

const PI: f64 = std::f64::consts::PI;

fn mkp(m: f64, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> reb_particle {
    let mut p = reb_particle::default();
    p.m = m;
    p.x = x;
    p.y = y;
    p.z = z;
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    p
}

fn close(name: &str, got: f64, want: f64, tol: f64) {
    let d = (got - want).abs();
    assert!(
        d <= tol,
        "{}: got {:.17e}, want {:.17e}, |diff| = {:.3e} > tol {:.3e}",
        name,
        got,
        want,
        d,
        tol
    );
}

fn bits(name: &str, got: f64, want: f64) {
    assert!(
        got.to_bits() == want.to_bits(),
        "{}: got {:.17e} (bits {:#018x}), want {:.17e} (bits {:#018x})",
        name,
        got,
        got.to_bits(),
        want,
        want.to_bits()
    );
}

fn close_vec(tag: &str, got: reb_vec3d, want: reb_vec3d, tol: f64) {
    close(&format!("{}.x", tag), got.x, want.x, tol);
    close(&format!("{}.y", tag), got.y, want.y, tol);
    close(&format!("{}.z", tag), got.z, want.z, tol);
}

fn pos_of(p: reb_particle) -> reb_vec3d {
    reb_vec3d { x: p.x, y: p.y, z: p.z }
}

fn vel_of(p: reb_particle) -> reb_vec3d {
    reb_vec3d { x: p.vx, y: p.vy, z: p.vz }
}

fn sub(a: reb_vec3d, b: reb_vec3d) -> reb_vec3d {
    reb_vec3d { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }
}

/// Positions/velocities agree to `tol` (absolute); mass is not compared
/// because several of the C transforms deliberately leave it untouched.
fn close_state(tag: &str, got: reb_particle, want: reb_particle, tol: f64) {
    close(&format!("{}.x", tag), got.x, want.x, tol);
    close(&format!("{}.y", tag), got.y, want.y, tol);
    close(&format!("{}.z", tag), got.z, want.z, tol);
    close(&format!("{}.vx", tag), got.vx, want.vx, tol);
    close(&format!("{}.vy", tag), got.vy, want.vy, tol);
    close(&format!("{}.vz", tag), got.vz, want.vz, tol);
}

/// Four massive bodies (the first with mass exactly 1.0, which makes some
/// of the transform steps exactly invertible in binary floating point)
/// followed by two massless test particles. Accelerations are filled in so
/// that the `_acc` variants of the transforms have something to chew on.
fn sample_system() -> Vec<reb_particle> {
    let mut v = vec![
        mkp(1.0, 0.10, -0.20, 0.05, 0.0100, -0.0200, 0.0030),
        mkp(3.00e-3, 1.30, 0.40, -0.15, -0.2300, 0.8100, 0.0410),
        mkp(9.50e-4, -2.10, 1.70, 0.33, 0.5100, -0.6200, -0.0170),
        mkp(2.75e-5, 0.40, -3.30, -0.90, 0.3300, 0.1400, 0.0750),
        mkp(0.0, 2.60, 1.10, 0.62, -0.1900, 0.4400, -0.0250),
        mkp(0.0, -1.40, -2.20, 0.08, 0.2700, -0.3100, 0.0620),
    ];
    for i in 0..v.len() {
        let f = (i as f64) + 1.0;
        v[i].ax = 0.031 * f;
        v[i].ay = -0.017 * f;
        v[i].az = 0.0091 * f;
    }
    v
}

const N_ACTIVE: usize = 4;

/// Same geometry as `sample_system`, but the two trailing particles carry
/// mass. They still sit beyond `N_ACTIVE`, so every transform must ignore
/// their mass entirely — which is what makes the N_active cut observable.
fn heavy_tail_system() -> Vec<reb_particle> {
    let mut v = sample_system();
    v[4].m = 7.0e-3;
    v[5].m = 4.0e-3;
    v
}

// ===========================================================================
// Jacobi coordinates
// ===========================================================================

#[test]
fn jacobi_round_trip_recovers_inertial_state() {
    let particles = sample_system();
    let N = particles.len();
    for &N_active in &[1usize, 2, N_ACTIVE, N] {
        let mut p_j = vec![reb_particle::default(); N];
        reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, N, N_active);
        let mut back = vec![reb_particle::default(); N];
        reb_transformations_jacobi_to_inertial_posvel(&mut back, &p_j, &particles, N, N_active);
        for i in 0..N {
            close_state(
                &format!("jacobi round trip N_active={} particle {}", N_active, i),
                back[i],
                particles[i],
                1.0e-13,
            );
        }
    }
}

#[test]
fn jacobi_first_coordinate_is_exact_relative_position() {
    // eta starts at p_mass[0].m == 1.0 exactly, so s_x*ei reduces to
    // particles[0].x bit-for-bit and p_j[1] is exactly the separation.
    let particles = sample_system();
    let N = particles.len();
    let mut p_j = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, N, N_ACTIVE);
    bits("p_j[1].x", p_j[1].x, particles[1].x - particles[0].x);
    bits("p_j[1].y", p_j[1].y, particles[1].y - particles[0].y);
    bits("p_j[1].z", p_j[1].z, particles[1].z - particles[0].z);
    bits("p_j[1].vx", p_j[1].vx, particles[1].vx - particles[0].vx);
    bits("p_j[1].vy", p_j[1].vy, particles[1].vy - particles[0].vy);
    bits("p_j[1].vz", p_j[1].vz, particles[1].vz - particles[0].vz);
    bits("p_j[1].m", p_j[1].m, particles[1].m);
}

#[test]
fn jacobi_zeroth_coordinate_is_center_of_mass() {
    // p_j[0] must be the barycentre of the whole system; derive it here by
    // a plain (independent) mass-weighted sum accumulated in reverse order.
    let particles = sample_system();
    let N = particles.len();
    let mut p_j = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, N, N);

    let (mut mx, mut my, mut mz, mut mvx, mut mvy, mut mvz, mut mt) = (0., 0., 0., 0., 0., 0., 0.);
    for i in (0..N).rev() {
        let p = particles[i];
        mx += p.m * p.x;
        my += p.m * p.y;
        mz += p.m * p.z;
        mvx += p.m * p.vx;
        mvy += p.m * p.vy;
        mvz += p.m * p.vz;
        mt += p.m;
    }
    close("jacobi com m", p_j[0].m, mt, 1.0e-15);
    close("jacobi com x", p_j[0].x, mx / mt, 1.0e-15);
    close("jacobi com y", p_j[0].y, my / mt, 1.0e-15);
    close("jacobi com z", p_j[0].z, mz / mt, 1.0e-15);
    close("jacobi com vx", p_j[0].vx, mvx / mt, 1.0e-15);
    close("jacobi com vy", p_j[0].vy, mvy / mt, 1.0e-15);
    close("jacobi com vz", p_j[0].vz, mvz / mt, 1.0e-15);
}

#[test]
fn jacobi_with_single_active_particle_is_heliocentric() {
    // With N_active == 1 the running sum never advances, so every Jacobi
    // coordinate is measured against particle 0. Mass 0 is exactly 1.0,
    // hence 1/eta == 1.0 and the subtraction is exact.
    let particles = sample_system();
    let N = particles.len();
    let mut p_j = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, N, 1);
    bits("p_j[0].x", p_j[0].x, particles[0].x);
    bits("p_j[0].vz", p_j[0].vz, particles[0].vz);
    bits("p_j[0].m", p_j[0].m, particles[0].m);
    for i in 1..N {
        bits(&format!("p_j[{}].x", i), p_j[i].x, particles[i].x - particles[0].x);
        bits(&format!("p_j[{}].y", i), p_j[i].y, particles[i].y - particles[0].y);
        bits(&format!("p_j[{}].z", i), p_j[i].z, particles[i].z - particles[0].z);
        bits(&format!("p_j[{}].vx", i), p_j[i].vx, particles[i].vx - particles[0].vx);
    }
}

#[test]
fn jacobi_single_particle_round_trip_is_bit_exact() {
    // N == 1 edge case: both loops are empty and the transform degenerates
    // to a multiply by 1.0 / divide by 1.0 pair.
    let particles = vec![sample_system()[0]];
    let mut p_j = vec![reb_particle::default(); 1];
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, 1, 1);
    let mut back = vec![reb_particle::default(); 1];
    reb_transformations_jacobi_to_inertial_posvel(&mut back, &p_j, &particles, 1, 1);
    bits("N=1 x", back[0].x, particles[0].x);
    bits("N=1 y", back[0].y, particles[0].y);
    bits("N=1 z", back[0].z, particles[0].z);
    bits("N=1 vx", back[0].vx, particles[0].vx);
    bits("N=1 vy", back[0].vy, particles[0].vy);
    bits("N=1 vz", back[0].vz, particles[0].vz);
}

#[test]
fn jacobi_acc_only_matches_posvelacc_bitwise() {
    // The acceleration-only kernel repeats the arithmetic of the combined
    // kernel expression for expression, so the results must be identical
    // bit patterns, not merely close.
    let particles = sample_system();
    let N = particles.len();
    let mut full = vec![reb_particle::default(); N];
    let mut only = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvelacc(&particles, &mut full, &particles, N, N_ACTIVE);
    reb_transformations_inertial_to_jacobi_acc(&particles, &mut only, &particles, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("acc kernel [{}].ax", i), only[i].ax, full[i].ax);
        bits(&format!("acc kernel [{}].ay", i), only[i].ay, full[i].ay);
        bits(&format!("acc kernel [{}].az", i), only[i].az, full[i].az);
    }
}

#[test]
fn jacobi_posvel_matches_posvelacc_bitwise() {
    let particles = sample_system();
    let N = particles.len();
    let mut full = vec![reb_particle::default(); N];
    let mut only = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvelacc(&particles, &mut full, &particles, N, N_ACTIVE);
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut only, &particles, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("posvel [{}].x", i), only[i].x, full[i].x);
        bits(&format!("posvel [{}].y", i), only[i].y, full[i].y);
        bits(&format!("posvel [{}].z", i), only[i].z, full[i].z);
        bits(&format!("posvel [{}].vx", i), only[i].vx, full[i].vx);
        bits(&format!("posvel [{}].vy", i), only[i].vy, full[i].vy);
        bits(&format!("posvel [{}].vz", i), only[i].vz, full[i].vz);
    }
}

#[test]
fn jacobi_to_inertial_pos_matches_posvel_bitwise() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_j = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_posvel(&particles, &mut p_j, &particles, N, N_ACTIVE);
    let mut a = vec![reb_particle::default(); N];
    let mut b = vec![reb_particle::default(); N];
    reb_transformations_jacobi_to_inertial_posvel(&mut a, &p_j, &particles, N, N_ACTIVE);
    reb_transformations_jacobi_to_inertial_pos(&mut b, &p_j, &particles, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("to_inertial_pos [{}].x", i), b[i].x, a[i].x);
        bits(&format!("to_inertial_pos [{}].y", i), b[i].y, a[i].y);
        bits(&format!("to_inertial_pos [{}].z", i), b[i].z, a[i].z);
    }
}

#[test]
fn jacobi_acceleration_round_trip() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_j = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_jacobi_acc(&particles, &mut p_j, &particles, N, N_ACTIVE);
    // jacobi_to_inertial_acc reads p_j[0].m as the total mass; the
    // acceleration-only forward kernel never writes it.
    let mut eta = 0.0;
    for i in 0..N_ACTIVE {
        eta += particles[i].m;
    }
    p_j[0].m = eta;
    let mut back = vec![reb_particle::default(); N];
    reb_transformations_jacobi_to_inertial_acc(&mut back, &p_j, &particles, N, N_ACTIVE);
    for i in 0..N {
        close(&format!("jacobi acc round trip [{}].ax", i), back[i].ax, particles[i].ax, 1.0e-14);
        close(&format!("jacobi acc round trip [{}].ay", i), back[i].ay, particles[i].ay, 1.0e-14);
        close(&format!("jacobi acc round trip [{}].az", i), back[i].az, particles[i].az, 1.0e-14);
    }
}

// ===========================================================================
// Democratic heliocentric coordinates
// ===========================================================================

#[test]
fn democraticheliocentric_relative_positions_are_exact() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_h = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_democraticheliocentric_posvel(
        &particles, &mut p_h, N, N_ACTIVE,
    );
    for i in 1..N {
        bits(&format!("dh [{}].x", i), p_h[i].x, particles[i].x - particles[0].x);
        bits(&format!("dh [{}].y", i), p_h[i].y, particles[i].y - particles[0].y);
        bits(&format!("dh [{}].z", i), p_h[i].z, particles[i].z - particles[0].z);
        bits(&format!("dh [{}].m", i), p_h[i].m, particles[i].m);
        // Velocities are measured against the barycentre of the active set.
        bits(&format!("dh [{}].vx", i), p_h[i].vx, particles[i].vx - p_h[0].vx);
    }
}

#[test]
fn democraticheliocentric_center_is_active_barycenter() {
    // Use the variant whose trailing (non-active) particles carry mass, so
    // that a transform which wrongly summed over all N would be caught.
    let particles = heavy_tail_system();
    let N = particles.len();
    let mut p_h = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_democraticheliocentric_posvel(
        &particles, &mut p_h, N, N_ACTIVE,
    );
    // Independent reverse-order accumulation of the same barycentre.
    let (mut mx, mut mvz, mut mt) = (0., 0., 0.);
    for i in (0..N_ACTIVE).rev() {
        mx += particles[i].m * particles[i].x;
        mvz += particles[i].m * particles[i].vz;
        mt += particles[i].m;
    }
    close("dh centre m", p_h[0].m, mt, 1.0e-16);
    close("dh centre x", p_h[0].x, mx / mt, 1.0e-16);
    close("dh centre vz", p_h[0].vz, mvz / mt, 1.0e-17);
    // The two non-active particles must not have contributed their mass.
    let all: f64 = particles.iter().map(|p| p.m).sum();
    assert!(
        all - p_h[0].m > 1.0e-2,
        "dh centre mass {} picked up the {} of mass sitting beyond N_active (total {})",
        p_h[0].m,
        all - p_h[0].m,
        all
    );
}

#[test]
fn democraticheliocentric_round_trip() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_h = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_democraticheliocentric_posvel(
        &particles, &mut p_h, N, N_ACTIVE,
    );
    // The reverse transform reads particles[0].m (the central mass), which
    // the C keeps in the destination array.
    let mut back = particles.clone();
    reb_transformations_democraticheliocentric_to_inertial_posvel(&mut back, &p_h, N, N_ACTIVE);
    for i in 0..N {
        close_state(&format!("dh round trip [{}]", i), back[i], particles[i], 1.0e-14);
    }
}

// ===========================================================================
// WHDS coordinates
// ===========================================================================

#[test]
fn whds_positions_match_democraticheliocentric_bitwise() {
    let particles = sample_system();
    let N = particles.len();
    let mut dh = vec![reb_particle::default(); N];
    let mut wh = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_democraticheliocentric_posvel(&particles, &mut dh, N, N_ACTIVE);
    reb_transformations_inertial_to_whds_posvel(&particles, &mut wh, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("whds/dh [{}].x", i), wh[i].x, dh[i].x);
        bits(&format!("whds/dh [{}].y", i), wh[i].y, dh[i].y);
        bits(&format!("whds/dh [{}].z", i), wh[i].z, dh[i].z);
    }
    bits("whds/dh centre vx", wh[0].vx, dh[0].vx);
    bits("whds/dh centre m", wh[0].m, dh[0].m);
}

#[test]
fn whds_velocities_are_mass_scaled_democraticheliocentric() {
    let particles = sample_system();
    let N = particles.len();
    let mut dh = vec![reb_particle::default(); N];
    let mut wh = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_democraticheliocentric_posvel(&particles, &mut dh, N, N_ACTIVE);
    reb_transformations_inertial_to_whds_posvel(&particles, &mut wh, N, N_ACTIVE);
    let m0 = particles[0].m;
    for i in 1..N_ACTIVE {
        let mf = (m0 + particles[i].m) / m0;
        bits(&format!("whds [{}].vx", i), wh[i].vx, mf * dh[i].vx);
        bits(&format!("whds [{}].vy", i), wh[i].vy, mf * dh[i].vy);
        bits(&format!("whds [{}].vz", i), wh[i].vz, mf * dh[i].vz);
        assert!(mf > 1.0, "mass factor for particle {} should exceed 1, got {}", i, mf);
    }
    for i in N_ACTIVE..N {
        // Test particles carry no mass factor.
        bits(&format!("whds test [{}].vx", i), wh[i].vx, dh[i].vx);
        bits(&format!("whds test [{}].vy", i), wh[i].vy, dh[i].vy);
        bits(&format!("whds test [{}].vz", i), wh[i].vz, dh[i].vz);
    }
}

#[test]
fn whds_round_trip() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_h = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_whds_posvel(&particles, &mut p_h, N, N_ACTIVE);
    let mut back = particles.clone();
    reb_transformations_whds_to_inertial_posvel(&mut back, &p_h, N, N_ACTIVE);
    for i in 0..N {
        close_state(&format!("whds round trip [{}]", i), back[i], particles[i], 1.0e-14);
    }
}

#[test]
fn whds_to_inertial_pos_matches_democraticheliocentric_pos() {
    // The C implements whds_to_inertial_pos by calling the democratic
    // heliocentric routine; the two must therefore agree bit for bit.
    let particles = sample_system();
    let N = particles.len();
    let mut p_h = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_whds_posvel(&particles, &mut p_h, N, N_ACTIVE);
    let mut a = particles.clone();
    let mut b = particles.clone();
    reb_transformations_whds_to_inertial_pos(&mut a, &p_h, N, N_ACTIVE);
    reb_transformations_democraticheliocentric_to_inertial_pos(&mut b, &p_h, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("whds pos [{}].x", i), a[i].x, b[i].x);
        bits(&format!("whds pos [{}].y", i), a[i].y, b[i].y);
        bits(&format!("whds pos [{}].z", i), a[i].z, b[i].z);
    }
}

// ===========================================================================
// Barycentric coordinates
// ===========================================================================

#[test]
fn barycentric_offsets_are_exact() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_b = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_barycentric_posvel(&particles, &mut p_b, N, N_ACTIVE);
    for i in 1..N {
        bits(&format!("bary [{}].x", i), p_b[i].x, particles[i].x - p_b[0].x);
        bits(&format!("bary [{}].y", i), p_b[i].y, particles[i].y - p_b[0].y);
        bits(&format!("bary [{}].z", i), p_b[i].z, particles[i].z - p_b[0].z);
        bits(&format!("bary [{}].vx", i), p_b[i].vx, particles[i].vx - p_b[0].vx);
    }
}

#[test]
fn barycentric_center_matches_independent_com() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_b = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_barycentric_posvel(&particles, &mut p_b, N, N_ACTIVE);

    let (mut mx, mut my, mut mz, mut mt) = (0., 0., 0., 0.);
    for i in (0..N_ACTIVE).rev() {
        mx += particles[i].m * particles[i].x;
        my += particles[i].m * particles[i].y;
        mz += particles[i].m * particles[i].z;
        mt += particles[i].m;
    }
    close("bary centre m", p_b[0].m, mt, 1.0e-16);
    close("bary centre x", p_b[0].x, mx / mt, 1.0e-16);
    close("bary centre y", p_b[0].y, my / mt, 1.0e-16);
    close("bary centre z", p_b[0].z, mz / mt, 1.0e-16);

    // Cross-check against the library's own centre-of-mass routine, which
    // uses a completely different (pairwise) accumulation.
    let mut sim = reb_simulation_create();
    for i in 0..N_ACTIVE {
        reb_simulation_add(&mut sim, particles[i]);
    }
    let com = reb_simulation_com(&sim);
    close("bary centre vs com x", p_b[0].x, com.x, 1.0e-15);
    close("bary centre vs com y", p_b[0].y, com.y, 1.0e-15);
    close("bary centre vs com z", p_b[0].z, com.z, 1.0e-15);
    close("bary centre vs com m", p_b[0].m, com.m, 1.0e-15);
}

#[test]
fn barycentric_round_trip() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_b = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_barycentric_posvel(&particles, &mut p_b, N, N_ACTIVE);
    let mut back = vec![reb_particle::default(); N];
    reb_transformations_barycentric_to_inertial_posvel(&mut back, &p_b, N, N_ACTIVE);
    for i in 0..N {
        close_state(&format!("bary round trip [{}]", i), back[i], particles[i], 1.0e-14);
    }
    for i in 0..N_ACTIVE {
        close(&format!("bary round trip [{}].m", i), back[i].m, particles[i].m, 1.0e-16);
    }
}

#[test]
fn barycentric_to_inertial_pos_matches_posvel_bitwise() {
    let particles = sample_system();
    let N = particles.len();
    let mut p_b = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_barycentric_posvel(&particles, &mut p_b, N, N_ACTIVE);
    let mut a = vec![reb_particle::default(); N];
    let mut b = vec![reb_particle::default(); N];
    reb_transformations_barycentric_to_inertial_posvel(&mut a, &p_b, N, N_ACTIVE);
    reb_transformations_barycentric_to_inertial_pos(&mut b, &p_b, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("bary pos [{}].x", i), b[i].x, a[i].x);
        bits(&format!("bary pos [{}].y", i), b[i].y, a[i].y);
        bits(&format!("bary pos [{}].z", i), b[i].z, a[i].z);
        bits(&format!("bary pos [{}].m", i), b[i].m, a[i].m);
    }
}

#[test]
fn barycentric_to_inertial_acc_mirrors_pos_kernel() {
    // barycentric_to_inertial_acc is the position kernel with x,y,z
    // replaced by ax,ay,az. Feeding the acceleration field in through the
    // position slots must therefore reproduce it bit for bit.
    let particles = sample_system();
    let N = particles.len();
    let mut p_b = vec![reb_particle::default(); N];
    reb_transformations_inertial_to_barycentric_posvel(&particles, &mut p_b, N, N_ACTIVE);
    for i in 0..N {
        p_b[i].ax = 0.7 + 0.13 * (i as f64);
        p_b[i].ay = -0.4 + 0.29 * (i as f64);
        p_b[i].az = 0.05 - 0.11 * (i as f64);
    }
    let mut shuffled = p_b.clone();
    for i in 0..N {
        shuffled[i].x = p_b[i].ax;
        shuffled[i].y = p_b[i].ay;
        shuffled[i].z = p_b[i].az;
    }
    let mut viaacc = vec![reb_particle::default(); N];
    let mut viapos = vec![reb_particle::default(); N];
    reb_transformations_barycentric_to_inertial_acc(&mut viaacc, &p_b, N, N_ACTIVE);
    reb_transformations_barycentric_to_inertial_pos(&mut viapos, &shuffled, N, N_ACTIVE);
    for i in 0..N {
        bits(&format!("bary acc [{}].ax", i), viaacc[i].ax, viapos[i].x);
        bits(&format!("bary acc [{}].ay", i), viaacc[i].ay, viapos[i].y);
        bits(&format!("bary acc [{}].az", i), viaacc[i].az, viapos[i].z);
    }
}

// ===========================================================================
// reb_vec3d algebra
// ===========================================================================

#[test]
fn vec3d_basis_cross_and_dot_are_exact() {
    let x = reb_vec3d { x: 1., y: 0., z: 0. };
    let y = reb_vec3d { x: 0., y: 1., z: 0. };
    let z = reb_vec3d { x: 0., y: 0., z: 1. };
    let c = reb_vec3d_cross(x, y);
    bits("x cross y .x", c.x, 0.0);
    bits("x cross y .y", c.y, 0.0);
    bits("x cross y .z", c.z, 1.0);
    let c = reb_vec3d_cross(y, z);
    bits("y cross z .x", c.x, 1.0);
    let c = reb_vec3d_cross(z, x);
    bits("z cross x .y", c.y, 1.0);
    bits("x dot x", reb_vec3d_dot(x, x), 1.0);
    bits("x dot y", reb_vec3d_dot(x, y), 0.0);
    bits("|z|^2", reb_vec3d_length_squared(z), 1.0);
}

#[test]
fn vec3d_cross_is_antisymmetric_and_orthogonal() {
    let a = reb_vec3d { x: 1.5, y: -2.25, z: 0.75 };
    let b = reb_vec3d { x: -0.5, y: 3.125, z: 2.0 };
    let ab = reb_vec3d_cross(a, b);
    let ba = reb_vec3d_cross(b, a);
    // Every product in the cross-product is exactly representable here, so
    // the antisymmetry holds bit for bit.
    bits("cross antisym x", ab.x, -ba.x);
    bits("cross antisym y", ab.y, -ba.y);
    bits("cross antisym z", ab.z, -ba.z);
    close("a . (a x b)", reb_vec3d_dot(a, ab), 0.0, 1.0e-15);
    close("b . (a x b)", reb_vec3d_dot(b, ab), 0.0, 1.0e-15);
    // Lagrange identity |a x b|^2 = |a|^2|b|^2 - (a.b)^2
    let lhs = reb_vec3d_length_squared(ab);
    let rhs = reb_vec3d_length_squared(a) * reb_vec3d_length_squared(b)
        - reb_vec3d_dot(a, b) * reb_vec3d_dot(a, b);
    close("Lagrange identity", lhs, rhs, 1.0e-13);
}

#[test]
fn vec3d_normalize_produces_unit_vector() {
    let v = reb_vec3d { x: 3.0, y: -4.0, z: 12.0 }; // length exactly 13
    let n = reb_vec3d_normalize(v);
    // 1/13 is not exact in binary, so compare with a tight tolerance.
    close("normalize x", n.x, 3.0 / 13.0, 1.0e-16);
    close("normalize y", n.y, -4.0 / 13.0, 1.0e-16);
    close("normalize z", n.z, 12.0 / 13.0, 1.0e-16);
    close("normalized length", reb_vec3d_length_squared(n), 1.0, 1.0e-15);
    // Scaling by a power of two leaves the direction bit-identical.
    let n2 = reb_vec3d_normalize(reb_vec3d_mul(v, 4.0));
    bits("normalize scale invariance x", n2.x, n.x);
    bits("normalize scale invariance y", n2.y, n.y);
    bits("normalize scale invariance z", n2.z, n.z);
    let s = reb_vec3d_add(v, reb_vec3d_mul(v, -1.0));
    bits("v + (-1)v", s.x, 0.0);
}

// ===========================================================================
// Quaternion algebra
// ===========================================================================

fn unit_quat() -> reb_rotation {
    // ix^2+iy^2+iz^2+r^2 = 4 * 0.25 = 1.0 exactly in binary.
    reb_rotation { ix: 0.5, iy: 0.5, iz: 0.5, r: 0.5 }
}

fn quat_close(tag: &str, got: reb_rotation, want: reb_rotation, tol: f64) {
    close(&format!("{}.ix", tag), got.ix, want.ix, tol);
    close(&format!("{}.iy", tag), got.iy, want.iy, tol);
    close(&format!("{}.iz", tag), got.iz, want.iz, tol);
    close(&format!("{}.r", tag), got.r, want.r, tol);
}

#[test]
fn rotation_identity_leaves_vector_bit_identical() {
    let q = reb_rotation_identity();
    bits("identity r", q.r, 1.0);
    bits("identity ix", q.ix, 0.0);
    let v = reb_vec3d { x: 1.25, y: -3.5, z: 0.125 };
    let w = reb_vec3d_rotate(v, q);
    bits("identity rotate x", w.x, v.x);
    bits("identity rotate y", w.y, v.y);
    bits("identity rotate z", w.z, v.z);
}

#[test]
fn rotation_conjugate_and_unit_inverse_are_exact() {
    let q = unit_quat();
    let c = reb_rotation_conjugate(q);
    bits("conj ix", c.ix, -q.ix);
    bits("conj r", c.r, q.r);
    let cc = reb_rotation_conjugate(c);
    bits("conj conj ix", cc.ix, q.ix);
    bits("conj conj iy", cc.iy, q.iy);
    bits("conj conj iz", cc.iz, q.iz);
    bits("conj conj r", cc.r, q.r);
    // |q|^2 == 1.0 exactly, so inverse == conjugate bit for bit and
    // normalize is the identity.
    let inv = reb_rotation_inverse(q);
    bits("inverse ix", inv.ix, c.ix);
    bits("inverse iy", inv.iy, c.iy);
    bits("inverse iz", inv.iz, c.iz);
    bits("inverse r", inv.r, c.r);
    let n = reb_rotation_normalize(q);
    bits("normalize ix", n.ix, q.ix);
    bits("normalize r", n.r, q.r);
    // q * q^-1 must be the identity quaternion, exactly for this q.
    let prod = reb_rotation_mul(q, inv);
    bits("q*qinv r", prod.r, 1.0);
    close("q*qinv ix", prod.ix, 0.0, 1.0e-17);
    close("q*qinv iy", prod.iy, 0.0, 1.0e-17);
    close("q*qinv iz", prod.iz, 0.0, 1.0e-17);
}

#[test]
fn rotation_normalize_rescales_to_unit_length() {
    let q = reb_rotation { ix: 1.0, iy: 2.0, iz: -2.0, r: 4.0 }; // |q| = 5
    let n = reb_rotation_normalize(q);
    close("normalized ix", n.ix, 0.2, 1.0e-16);
    close("normalized iy", n.iy, 0.4, 1.0e-16);
    close("normalized iz", n.iz, -0.4, 1.0e-16);
    close("normalized r", n.r, 0.8, 1.0e-16);
    let l2 = n.r * n.r + n.ix * n.ix + n.iy * n.iy + n.iz * n.iz;
    close("normalized length^2", l2, 1.0, 1.0e-15);
    // A non-unit quaternion still rotates like its normalized version once
    // the inverse (which divides by |q|^2) is used.
    let v = reb_vec3d { x: 0.3, y: -1.1, z: 2.4 };
    let rotated = reb_vec3d_rotate(v, n);
    close("length preserved", reb_vec3d_length_squared(rotated), reb_vec3d_length_squared(v), 1.0e-14);
}

#[test]
fn rotation_inverse_undoes_rotation() {
    let q = reb_rotation_init_angle_axis(0.937, reb_vec3d { x: 1.0, y: -2.0, z: 0.5 });
    let qi = reb_rotation_inverse(q);
    let v = reb_vec3d { x: 0.3, y: -1.1, z: 2.4 };
    let back = reb_vec3d_rotate(reb_vec3d_rotate(v, q), qi);
    close_vec("rotate then inverse", back, v, 1.0e-15);
    // reb_vec3d_rotate must not touch its argument; irotate must.
    let mut m = v;
    reb_vec3d_irotate(&mut m, q);
    let r = reb_vec3d_rotate(v, q);
    bits("rotate vs irotate x", m.x, r.x);
    bits("rotate vs irotate y", m.y, r.y);
    bits("rotate vs irotate z", m.z, r.z);
    assert!(m != v, "irotate left the vector unchanged for a 0.937 rad rotation");
}

#[test]
fn rotation_mul_composes_rotations() {
    let q1 = reb_rotation_init_angle_axis(0.4, reb_vec3d { x: 0.0, y: 0.0, z: 1.0 });
    let q2 = reb_rotation_init_angle_axis(1.3, reb_vec3d { x: 1.0, y: 0.0, z: 0.0 });
    let v = reb_vec3d { x: 0.7, y: 2.1, z: -0.9 };
    let step = reb_vec3d_rotate(reb_vec3d_rotate(v, q1), q2);
    let once = reb_vec3d_rotate(v, reb_rotation_mul(q2, q1));
    close_vec("mul composes", once, step, 1.0e-15);
    // Two rotations about the same axis add their angles.
    let a = reb_rotation_init_angle_axis(0.4, reb_vec3d { x: 0.0, y: 0.0, z: 1.0 });
    let b = reb_rotation_init_angle_axis(0.9, reb_vec3d { x: 0.0, y: 0.0, z: 1.0 });
    let sum = reb_rotation_init_angle_axis(1.3, reb_vec3d { x: 0.0, y: 0.0, z: 1.0 });
    quat_close("same-axis composition", reb_rotation_mul(b, a), sum, 1.0e-15);
}

#[test]
fn rotation_angle_axis_quarter_turn() {
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let q = reb_rotation_init_angle_axis(PI / 2.0, z);
    let x = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    let r = reb_vec3d_rotate(x, q);
    close_vec("quarter turn of x", r, reb_vec3d { x: 0.0, y: 1.0, z: 0.0 }, 1.0e-15);
    let r2 = reb_vec3d_rotate(reb_vec3d { x: 0.0, y: 1.0, z: 0.0 }, q);
    close_vec("quarter turn of y", r2, reb_vec3d { x: -1.0, y: 0.0, z: 0.0 }, 1.0e-15);
    // A rotation about z leaves z alone, exactly.
    let rz = reb_vec3d_rotate(z, q);
    close_vec("quarter turn of z", rz, z, 1.0e-16);
    // Four quarter turns are the identity.
    let mut w = reb_vec3d { x: 1.3, y: -0.4, z: 2.2 };
    for _ in 0..4 {
        reb_vec3d_irotate(&mut w, q);
    }
    close_vec("four quarter turns", w, reb_vec3d { x: 1.3, y: -0.4, z: 2.2 }, 1.0e-15);
}

#[test]
fn rotation_init_from_to_maps_vectors() {
    let from = reb_vec3d { x: 2.0, y: 1.0, z: -0.5 };
    let to = reb_vec3d { x: -1.0, y: 3.0, z: 2.0 };
    let q = reb_rotation_init_from_to(from, to);
    let mapped = reb_vec3d_rotate(reb_vec3d_normalize(from), q);
    close_vec("from->to", mapped, reb_vec3d_normalize(to), 1.0e-15);
    // The same must work for a pair more than 90 degrees apart, which
    // takes the two-stage branch.
    let to2 = reb_vec3d { x: -2.0, y: -1.5, z: 0.4 };
    assert!(
        reb_vec3d_dot(reb_vec3d_normalize(from), reb_vec3d_normalize(to2)) < 0.0,
        "test vectors were meant to be more than 90 degrees apart"
    );
    let q2 = reb_rotation_init_from_to(from, to2);
    let mapped2 = reb_vec3d_rotate(reb_vec3d_normalize(from), q2);
    close_vec("from->to obtuse", mapped2, reb_vec3d_normalize(to2), 1.0e-14);
}

#[test]
fn rotation_init_from_to_antipodal_is_exact() {
    // from = +x, to = -x: the half vector is (0,0,0), normalize yields NaN
    // and the C falls through to the "pick an orthogonal axis" branch,
    // which returns the pi rotation about +z. Every step is exact.
    let x = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    let mx = reb_vec3d { x: -1.0, y: 0.0, z: 0.0 };
    let q = reb_rotation_init_from_to(x, mx);
    bits("antipodal ix", q.ix, 0.0);
    bits("antipodal iy", q.iy, 0.0);
    bits("antipodal iz", q.iz, 1.0);
    bits("antipodal r", q.r, 0.0);
    let mapped = reb_vec3d_rotate(x, q);
    bits("antipodal maps x to -x (x)", mapped.x, -1.0);
    bits("antipodal maps x to -x (y)", mapped.y, 0.0);
    bits("antipodal maps x to -x (z)", mapped.z, 0.0);
    // z is on the rotation axis and must survive untouched.
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let zz = reb_vec3d_rotate(z, q);
    bits("antipodal fixes z", zz.z, 1.0);
}

#[test]
fn rotation_init_to_new_axes_builds_frame() {
    // newz must already be a unit vector: the C computes the
    // orthogonalisation dot product before normalizing it.
    let newz = reb_vec3d_normalize(reb_vec3d { x: 0.0, y: 0.6, z: 0.8 });
    let newx = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    let q = reb_rotation_init_to_new_axes(newz, newx);
    close_vec("newz -> z", reb_vec3d_rotate(newz, q), reb_vec3d { x: 0., y: 0., z: 1. }, 1.0e-14);
    close_vec("newx -> x", reb_vec3d_rotate(newx, q), reb_vec3d { x: 1., y: 0., z: 0. }, 1.0e-14);
    // Rotations preserve cross products, so the implied y axis follows.
    let newy = reb_vec3d_cross(newz, newx);
    close_vec("newy -> y", reb_vec3d_rotate(newy, q), reb_vec3d { x: 0., y: 1., z: 0. }, 1.0e-14);
}

#[test]
fn rotation_init_orbit_matches_particle_from_orbit() {
    // reb_rotation_init_orbit(Omega, inc, omega) is exactly the rotation
    // that carries the in-plane (Omega=inc=omega=0) orbit into the
    // three-dimensional one built by reb_particle_from_orbit.
    let G = 1.0;
    let primary = mkp(1.0, 0., 0., 0., 0., 0., 0.);
    for &(a, e, inc, Om, om, f) in &[
        (1.0, 0.15, 0.40, 0.70, 1.10, 2.00),
        (2.5, 0.62, 2.30, 5.10, 0.35, 4.10),
        (1.7, 0.00, 1.05, 3.10, 2.20, 0.50),
    ] {
        let flat = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, 0.0, 0.0, 0.0, f);
        let full = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let q = reb_rotation_init_orbit(Om, inc, om);
        let mut rotated = flat;
        reb_particle_irotate(&mut rotated, q);
        close_state(&format!("init_orbit inc={}", inc), rotated, full, 1.0e-14);
    }
}

#[test]
fn rotation_to_orbital_round_trips_generic_angles() {
    let Om = 0.7;
    let inc = 0.4;
    let om = 1.1;
    let q = reb_rotation_init_orbit(Om, inc, om);
    let (mut gO, mut gi, mut go) = (0., 0., 0.);
    reb_rotation_to_orbital(q, &mut gO, &mut gi, &mut go);
    close("to_orbital Omega", gO, Om, 1.0e-14);
    close("to_orbital inc", gi, inc, 1.0e-14);
    close("to_orbital omega", go, om, 1.0e-14);
    // Rebuilding from the recovered angles reproduces the quaternion.
    quat_close("orbital round trip", reb_rotation_init_orbit(gO, gi, go), q, 1.0e-14);
}

#[test]
fn rotation_to_orbital_degenerate_branches_are_exact() {
    // Identity: inclination zero and both nodes at zero.
    let (mut gO, mut gi, mut go) = (1., 1., 1.);
    reb_rotation_to_orbital(reb_rotation_identity(), &mut gO, &mut gi, &mut go);
    bits("identity Omega", gO, 0.0);
    bits("identity inc", gi, 0.0);
    bits("identity omega", go, 0.0);

    // pi rotation about z: inc = 0 (the |inc| <= MIN_INC branch), and only
    // the sum Omega+omega is determined, which the C reports as omega.
    let qz = reb_rotation { ix: 0.0, iy: 0.0, iz: 1.0, r: 0.0 };
    reb_rotation_to_orbital(qz, &mut gO, &mut gi, &mut go);
    bits("z-flip Omega", gO, 0.0);
    bits("z-flip inc", gi, 0.0);
    close("z-flip omega", go, PI, 1.0e-15);

    // pi rotation about x: inc = pi (the |inc - pi| <= MIN_INC branch).
    let qx = reb_rotation { ix: 1.0, iy: 0.0, iz: 0.0, r: 0.0 };
    reb_rotation_to_orbital(qx, &mut gO, &mut gi, &mut go);
    bits("x-flip Omega", gO, 0.0);
    close("x-flip inc", gi, PI, 1.0e-15);
    bits("x-flip omega", go, 0.0);
    // Sanity: that quaternion really is the retrograde flip.
    close_vec(
        "x-flip maps z to -z",
        reb_vec3d_rotate(reb_vec3d { x: 0., y: 0., z: 1. }, qx),
        reb_vec3d { x: 0., y: 0., z: -1. },
        1.0e-16,
    );
}

#[test]
fn rotation_slerp_endpoints_and_midpoint() {
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let q1 = reb_rotation_identity();
    let q2 = reb_rotation_init_angle_axis(PI / 2.0, z);
    quat_close("slerp t=0", reb_rotation_slerp(q1, q2, 0.0), q1, 1.0e-15);
    quat_close("slerp t=1", reb_rotation_slerp(q1, q2, 1.0), q2, 1.0e-15);
    // Halfway between the identity and a quarter turn is an eighth turn.
    let mid = reb_rotation_slerp(q1, q2, 0.5);
    quat_close("slerp midpoint", mid, reb_rotation_init_angle_axis(PI / 4.0, z), 1.0e-15);
    // Identical inputs short-circuit and return q1 unchanged, bit for bit.
    let u = unit_quat();
    let same = reb_rotation_slerp(u, u, 0.37);
    bits("slerp q==q ix", same.ix, u.ix);
    bits("slerp q==q r", same.r, u.r);
}

#[test]
fn rotation_slerp_small_angle_branch() {
    // sin(halfTheta) below QUATERNION_EPS (1e-4) takes the linear-blend
    // branch; for a 1e-5 rad separation the blend is accurate to O(angle^3).
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let q1 = reb_rotation_identity();
    let q2 = reb_rotation_init_angle_axis(1.0e-5, z);
    let cosHalfTheta = q1.r * q2.r + q1.ix * q2.ix + q1.iy * q2.iy + q1.iz * q2.iz;
    assert!(cosHalfTheta.abs() < 1.0, "test setup no longer reaches the slerp blend branch");
    assert!(
        (1.0 - cosHalfTheta * cosHalfTheta).sqrt() < 1.0e-4,
        "test setup no longer reaches the small-angle slerp branch"
    );
    let mid = reb_rotation_slerp(q1, q2, 0.5);
    let want = reb_rotation_init_angle_axis(0.5e-5, z);
    // The branch returns the *unnormalised* midpoint 0.5*(q1+q2), whose
    // length is cos(theta/4) = 1 - theta^2/32 rather than 1. Its direction
    // is the true half rotation, so normalising recovers it.
    let l2 = mid.r * mid.r + mid.ix * mid.ix + mid.iy * mid.iy + mid.iz * mid.iz;
    let deficit = 1.0 - l2.sqrt();
    assert!(
        deficit > 0.0 && deficit < 1.0e-10,
        "linear-blend slerp length deficit {:e} is not the expected O(theta^2/32) = 3.1e-12",
        deficit
    );
    quat_close("small-angle slerp (normalised)", reb_rotation_normalize(mid), want, 1.0e-15);
    let v = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    close_vec(
        "small-angle slerp acting on x",
        reb_vec3d_rotate(v, reb_rotation_normalize(mid)),
        reb_vec3d_rotate(v, want),
        1.0e-15,
    );
}

#[test]
fn particle_irotate_rotates_position_and_velocity() {
    let q = reb_rotation_init_angle_axis(1.234, reb_vec3d { x: 0.3, y: -1.0, z: 2.0 });
    let p = mkp(1.5e-3, 0.7, -2.1, 0.55, 0.11, 0.42, -0.07);
    let mut rp = p;
    reb_particle_irotate(&mut rp, q);
    close_vec("particle pos rotated", pos_of(rp), reb_vec3d_rotate(pos_of(p), q), 1.0e-16);
    close_vec("particle vel rotated", vel_of(rp), reb_vec3d_rotate(vel_of(p), q), 1.0e-16);
    bits("particle mass untouched", rp.m, p.m);
    // Lengths, the dot product and the angular momentum magnitude are all
    // rotation invariants.
    close(
        "|r| invariant",
        reb_vec3d_length_squared(pos_of(rp)),
        reb_vec3d_length_squared(pos_of(p)),
        1.0e-15,
    );
    close(
        "|v| invariant",
        reb_vec3d_length_squared(vel_of(rp)),
        reb_vec3d_length_squared(vel_of(p)),
        1.0e-15,
    );
    close(
        "r.v invariant",
        reb_vec3d_dot(pos_of(rp), vel_of(rp)),
        reb_vec3d_dot(pos_of(p), vel_of(p)),
        1.0e-15,
    );
    close(
        "|r x v| invariant",
        reb_vec3d_length_squared(reb_vec3d_cross(pos_of(rp), vel_of(rp))),
        reb_vec3d_length_squared(reb_vec3d_cross(pos_of(p), vel_of(p))),
        1.0e-15,
    );
}

fn build_sim() -> reb_simulation {
    let mut sim = reb_simulation_create();
    sim.G = 1.0;
    let sys = sample_system();
    for i in 0..N_ACTIVE {
        reb_simulation_add(&mut sim, sys[i]);
    }
    sim
}

#[test]
fn simulation_irotate_preserves_energy_and_angular_momentum() {
    let mut sim = build_sim();
    let e0 = reb_simulation_energy(&sim);
    let l0 = reb_simulation_angular_momentum(&sim);
    let l0mag = reb_vec3d_length_squared(l0).sqrt();
    assert!(l0mag > 1.0e-6, "test system has a degenerate angular momentum ({})", l0mag);
    let q = reb_rotation_init_angle_axis(0.77, reb_vec3d { x: 1.0, y: 2.0, z: -0.5 });
    reb_simulation_irotate(&mut sim, q);
    let e1 = reb_simulation_energy(&sim);
    let l1 = reb_simulation_angular_momentum(&sim);
    close("energy under rotation", e1, e0, 1.0e-14 * e0.abs().max(1.0e-6));
    close(
        "|L| under rotation",
        reb_vec3d_length_squared(l1).sqrt(),
        l0mag,
        1.0e-14 * l0mag,
    );
    // The angular momentum vector itself must follow the rotation.
    close_vec("L rotates", l1, reb_vec3d_rotate(l0, q), 1.0e-14 * l0mag);
}

#[test]
fn simulation_irotate_inverse_round_trip() {
    let mut sim = build_sim();
    let before: Vec<reb_particle> = sim.particles.clone();
    let q = reb_rotation_init_angle_axis(2.13, reb_vec3d { x: -0.4, y: 0.9, z: 1.7 });
    reb_simulation_irotate(&mut sim, q);
    let mut moved = false;
    for i in 0..sim.N {
        if sim.particles[i].x != before[i].x {
            moved = true;
        }
    }
    assert!(moved, "reb_simulation_irotate left every particle in place");
    reb_simulation_irotate(&mut sim, reb_rotation_inverse(q));
    for i in 0..sim.N {
        close_state(&format!("sim irotate round trip [{}]", i), sim.particles[i], before[i], 1.0e-14);
    }
}

// ===========================================================================
// reb_mat4df helpers (single precision, used by the visualisation code)
// ===========================================================================

fn sample_mat() -> reb_mat4df {
    reb_mat4df {
        m: [
            1.0, 2.0, 3.0, 4.0, //
            5.0, 6.0, 7.0, 8.0, //
            9.0, 10.0, 11.0, 12.0, //
            13.0, 14.0, 15.0, 16.0,
        ],
    }
}

#[test]
fn mat4df_identity_is_multiplicative_unit() {
    let id = reb_mat4df_identity();
    let a = sample_mat();
    assert!(reb_mat4df_eq(reb_mat4df_multiply(a, id), a), "A * I != A");
    assert!(reb_mat4df_eq(reb_mat4df_multiply(id, a), a), "I * A != A");
    assert!(reb_mat4df_eq(id, reb_mat4df_identity()), "identity is not self-equal");
    let mut b = a;
    b.m[7] = 99.0;
    assert!(!reb_mat4df_eq(a, b), "reb_mat4df_eq did not notice a changed entry");
}

#[test]
fn mat4df_scale_and_get_scale() {
    let s = reb_mat4df_scale(reb_mat4df_identity(), 2.0, 3.0, 4.0);
    // Column-major storage: the diagonal picks up the scale factors.
    assert!(s.m[0] == 2.0, "scale m[0] = {}", s.m[0]);
    assert!(s.m[5] == 3.0, "scale m[5] = {}", s.m[5]);
    assert!(s.m[10] == 4.0, "scale m[10] = {}", s.m[10]);
    assert!(s.m[15] == 1.0, "scale m[15] = {}", s.m[15]);
    for &j in &[1usize, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14] {
        assert!(s.m[j] == 0.0, "scale m[{}] = {} should be zero", j, s.m[j]);
    }
    // 2, 3 and 4 all have exact squares and square roots in binary32.
    let g = reb_mat4df_get_scale(s);
    assert!(g.x == 2.0 && g.y == 3.0 && g.z == 4.0, "get_scale returned {:?}", g);
    // Scaling twice multiplies the factors.
    let s2 = reb_mat4df_scale(s, 0.5, 2.0, 0.25);
    let g2 = reb_mat4df_get_scale(s2);
    assert!(g2.x == 1.0 && g2.y == 6.0 && g2.z == 1.0, "double scale returned {:?}", g2);
}

#[test]
fn mat4df_translate_writes_the_translation_column() {
    let t = reb_mat4df_translate(reb_mat4df_identity(), 1.5, -2.25, 0.75);
    assert!(t.m[3] == 1.5, "translate m[3] = {}", t.m[3]);
    assert!(t.m[7] == -2.25, "translate m[7] = {}", t.m[7]);
    assert!(t.m[11] == 0.75, "translate m[11] = {}", t.m[11]);
    for &j in &[0usize, 5, 10, 15] {
        assert!(t.m[j] == 1.0, "translate diagonal m[{}] = {}", j, t.m[j]);
    }
    // Translating twice adds the offsets.
    let t2 = reb_mat4df_translate(t, 0.5, 0.25, -0.75);
    assert!(t2.m[3] == 2.0, "double translate m[3] = {}", t2.m[3]);
    assert!(t2.m[7] == -2.0, "double translate m[7] = {}", t2.m[7]);
    assert!(t2.m[11] == 0.0, "double translate m[11] = {}", t2.m[11]);
}

#[test]
fn mat4df_ortho_maps_the_box_onto_the_clip_cube() {
    // l,r,b,t,n,f all chosen so every quotient is exact in binary32.
    let m = reb_mat4df_ortho(0.0, 4.0, 0.0, 2.0, 1.0, 3.0);
    assert!(m.m[0] == 0.5, "ortho m[0] = {}", m.m[0]);
    assert!(m.m[3] == -1.0, "ortho m[3] = {}", m.m[3]);
    assert!(m.m[5] == 1.0, "ortho m[5] = {}", m.m[5]);
    assert!(m.m[7] == -1.0, "ortho m[7] = {}", m.m[7]);
    assert!(m.m[10] == -1.0, "ortho m[10] = {}", m.m[10]);
    assert!(m.m[11] == -2.0, "ortho m[11] = {}", m.m[11]);
    assert!(m.m[15] == 1.0, "ortho m[15] = {}", m.m[15]);
    // Row i of the matrix maps a point to clip space: x' = m0 x + m3, etc.
    let clip_x = |x: f32| m.m[0] * x + m.m[3];
    assert!(clip_x(0.0) == -1.0, "left edge maps to {}", clip_x(0.0));
    assert!(clip_x(4.0) == 1.0, "right edge maps to {}", clip_x(4.0));
    let clip_y = |y: f32| m.m[5] * y + m.m[7];
    assert!(clip_y(0.0) == -1.0, "bottom edge maps to {}", clip_y(0.0));
    assert!(clip_y(2.0) == 1.0, "top edge maps to {}", clip_y(2.0));
    let clip_z = |z: f32| m.m[10] * z + m.m[11];
    assert!(clip_z(1.0) == -3.0, "near plane maps to {}", clip_z(1.0));
    assert!(clip_z(3.0) == -5.0, "far plane maps to {}", clip_z(3.0));
}

#[test]
fn rotation_to_mat4df_agrees_with_quaternion_rotation() {
    // Identity quaternion -> identity matrix, exactly.
    assert!(
        reb_mat4df_eq(reb_rotation_to_mat4df(reb_rotation_identity()), reb_mat4df_identity()),
        "identity quaternion did not produce the identity matrix"
    );
    // A general rotation: the 3x3 block must reproduce reb_vec3d_rotate to
    // single precision.
    let q = reb_rotation_init_angle_axis(0.83, reb_vec3d { x: 0.3, y: -1.1, z: 2.0 });
    let m = reb_rotation_to_mat4df(q);
    for v in [
        reb_vec3d { x: 1.0, y: 0.0, z: 0.0 },
        reb_vec3d { x: 0.0, y: 1.0, z: 0.0 },
        reb_vec3d { x: 0.0, y: 0.0, z: 1.0 },
        reb_vec3d { x: 0.4, y: -0.7, z: 1.3 },
    ] {
        let w = reb_vec3d_rotate(v, q);
        let mx = m.m[0] as f64 * v.x + m.m[1] as f64 * v.y + m.m[2] as f64 * v.z;
        let my = m.m[4] as f64 * v.x + m.m[5] as f64 * v.y + m.m[6] as f64 * v.z;
        let mz = m.m[8] as f64 * v.x + m.m[9] as f64 * v.y + m.m[10] as f64 * v.z;
        close("mat4df rotate x", mx, w.x, 1.0e-6);
        close("mat4df rotate y", my, w.y, 1.0e-6);
        close("mat4df rotate z", mz, w.z, 1.0e-6);
    }
    // The homogeneous row/column must be untouched.
    assert!(m.m[3] == 0.0 && m.m[7] == 0.0 && m.m[11] == 0.0, "translation column is not zero");
    assert!(
        m.m[12] == 0.0 && m.m[13] == 0.0 && m.m[14] == 0.0 && m.m[15] == 1.0,
        "bottom row is not (0,0,0,1)"
    );
}

// ===========================================================================
// Orbital derivatives: finite-difference machinery
// ===========================================================================

type DerivFn = fn(f64, reb_particle, reb_particle) -> reb_particle;
type Eval = fn(Ctx, [f64; 7]) -> [f64; 6];

#[derive(Clone, Copy)]
struct Ctx {
    G: f64,
    primary: reb_particle,
}

// Pal (2009) parameter slots.
const P_A: usize = 0;
const P_LAMBDA: usize = 1;
const P_K: usize = 2;
const P_H: usize = 3;
const P_IX: usize = 4;
const P_IY: usize = 5;
const P_M: usize = 6;

// Keplerian parameter slots.
const K_A: usize = 0;
const K_E: usize = 1;
const K_INC: usize = 2;
const K_OMEGA_BIG: usize = 3;
const K_OMEGA: usize = 4;
const K_F: usize = 5;
const K_M: usize = 6;

fn pal_eval(c: Ctx, v: [f64; 7]) -> [f64; 6] {
    let p = reb_particle_from_pal(c.G, c.primary, v[6], v[0], v[1], v[2], v[3], v[4], v[5]);
    [p.x, p.y, p.z, p.vx, p.vy, p.vz]
}

fn kep_eval(c: Ctx, v: [f64; 7]) -> [f64; 6] {
    let p = reb_particle_from_orbit(c.G, c.primary, v[6], v[0], v[1], v[2], v[3], v[4], v[5]);
    [p.x, p.y, p.z, p.vx, p.vy, p.vz]
}

/// The Pal element vector the derivative functions themselves recover
/// from `po` (so the analytic and numerical derivatives are taken at
/// exactly the same point).
fn pal_base(G: f64, primary: reb_particle, po: reb_particle) -> [f64; 7] {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);
    [a, lambda, k, h, ix, iy, po.m]
}

fn kep_base(G: f64, primary: reb_particle, po: reb_particle) -> [f64; 7] {
    let o = reb_orbit_from_particle(G, po, primary);
    [o.a, o.e, o.inc, o.Omega, o.omega, o.f, po.m]
}

fn fd_first(ev: Eval, c: Ctx, v: [f64; 7], i: usize, d: f64) -> [f64; 6] {
    let mut vp = v;
    vp[i] += d;
    let mut vm = v;
    vm[i] -= d;
    let a = ev(c, vp);
    let b = ev(c, vm);
    let mut o = [0.0; 6];
    for t in 0..6 {
        o[t] = (a[t] - b[t]) / (2.0 * d);
    }
    o
}

fn fd_second(ev: Eval, c: Ctx, v: [f64; 7], i: usize, j: usize, d: f64) -> [f64; 6] {
    let mut o = [0.0; 6];
    if i == j {
        let mut vp = v;
        vp[i] += d;
        let mut vm = v;
        vm[i] -= d;
        let a = ev(c, vp);
        let b = ev(c, v);
        let e = ev(c, vm);
        for t in 0..6 {
            o[t] = ((a[t] - b[t]) - (b[t] - e[t])) / (d * d);
        }
    } else {
        let mut vpp = v;
        vpp[i] += d;
        vpp[j] += d;
        let mut vpm = v;
        vpm[i] += d;
        vpm[j] -= d;
        let mut vmp = v;
        vmp[i] -= d;
        vmp[j] += d;
        let mut vmm = v;
        vmm[i] -= d;
        vmm[j] -= d;
        let app = ev(c, vpp);
        let apm = ev(c, vpm);
        let amp = ev(c, vmp);
        let amm = ev(c, vmm);
        for t in 0..6 {
            o[t] = ((app[t] - apm[t]) - (amp[t] - amm[t])) / (4.0 * d * d);
        }
    }
    o
}

fn compare_deriv(tag: &str, analytic: reb_particle, fd: [f64; 6], rtol: f64) {
    let an = [analytic.x, analytic.y, analytic.z, analytic.vx, analytic.vy, analytic.vz];
    let comp = ["x", "y", "z", "vx", "vy", "vz"];
    let mut scale = 0.0f64;
    for t in 0..6 {
        scale = scale.max(an[t].abs()).max(fd[t].abs());
    }
    assert!(
        scale > 0.0 && scale.is_finite(),
        "{}: analytic and finite-difference derivatives are both identically zero (or not finite)",
        tag
    );
    let tol = rtol * scale;
    for t in 0..6 {
        let d = (an[t] - fd[t]).abs();
        assert!(
            d <= tol,
            "{} d/{}: analytic {:.12e} vs finite difference {:.12e} (|diff| {:.3e} > tol {:.3e}, scale {:.3e})",
            tag,
            comp[t],
            an[t],
            fd[t],
            d,
            tol,
            scale
        );
    }
}

// Central-difference step sizes and the relative tolerances they support.
// For the first derivatives the measured analytic-vs-numeric disagreement
// sits below 2e-10 of the derivative's own scale; for the second
// derivatives, where the O(eps/d^2) roundoff and the O(d^2) truncation
// terms meet, it sits near 5e-8. The tolerances below keep a wide margin
// over those floors while still being far tighter than any real error in
// the analytic expressions could hide under.
const FD1_STEP: f64 = 1.0e-6;
const FD2_STEP: f64 = 1.0e-4;
const FD1_RTOL: f64 = 1.0e-8;
const FD2_RTOL: f64 = 2.0e-6;

fn deriv_primary() -> reb_particle {
    mkp(1.0, 0.30, -0.20, 0.15, 0.020, -0.010, 0.030)
}

fn pal_first_table() -> Vec<(&'static str, usize, DerivFn)> {
    vec![
        ("a", P_A, reb_particle_derivative_a as DerivFn),
        ("lambda", P_LAMBDA, reb_particle_derivative_lambda as DerivFn),
        ("k", P_K, reb_particle_derivative_k as DerivFn),
        ("h", P_H, reb_particle_derivative_h as DerivFn),
        ("ix", P_IX, reb_particle_derivative_ix as DerivFn),
        ("iy", P_IY, reb_particle_derivative_iy as DerivFn),
        ("m", P_M, reb_particle_derivative_m as DerivFn),
    ]
}

fn pal_second_table() -> Vec<(&'static str, usize, usize, DerivFn)> {
    vec![
        ("a_a", P_A, P_A, reb_particle_derivative_a_a as DerivFn),
        ("lambda_lambda", P_LAMBDA, P_LAMBDA, reb_particle_derivative_lambda_lambda as DerivFn),
        ("k_k", P_K, P_K, reb_particle_derivative_k_k as DerivFn),
        ("h_h", P_H, P_H, reb_particle_derivative_h_h as DerivFn),
        ("ix_ix", P_IX, P_IX, reb_particle_derivative_ix_ix as DerivFn),
        ("iy_iy", P_IY, P_IY, reb_particle_derivative_iy_iy as DerivFn),
        ("m_m", P_M, P_M, reb_particle_derivative_m_m as DerivFn),
        ("k_lambda", P_K, P_LAMBDA, reb_particle_derivative_k_lambda as DerivFn),
        ("h_lambda", P_H, P_LAMBDA, reb_particle_derivative_h_lambda as DerivFn),
        ("k_h", P_K, P_H, reb_particle_derivative_k_h as DerivFn),
        ("k_ix", P_K, P_IX, reb_particle_derivative_k_ix as DerivFn),
        ("h_ix", P_H, P_IX, reb_particle_derivative_h_ix as DerivFn),
        ("lambda_ix", P_LAMBDA, P_IX, reb_particle_derivative_lambda_ix as DerivFn),
        ("lambda_iy", P_LAMBDA, P_IY, reb_particle_derivative_lambda_iy as DerivFn),
        ("h_iy", P_H, P_IY, reb_particle_derivative_h_iy as DerivFn),
        ("k_iy", P_K, P_IY, reb_particle_derivative_k_iy as DerivFn),
        ("ix_iy", P_IX, P_IY, reb_particle_derivative_ix_iy as DerivFn),
        ("a_ix", P_A, P_IX, reb_particle_derivative_a_ix as DerivFn),
        ("a_iy", P_A, P_IY, reb_particle_derivative_a_iy as DerivFn),
        ("a_lambda", P_A, P_LAMBDA, reb_particle_derivative_a_lambda as DerivFn),
        ("a_h", P_A, P_H, reb_particle_derivative_a_h as DerivFn),
        ("a_k", P_A, P_K, reb_particle_derivative_a_k as DerivFn),
        ("m_a", P_M, P_A, reb_particle_derivative_m_a as DerivFn),
        ("m_lambda", P_M, P_LAMBDA, reb_particle_derivative_m_lambda as DerivFn),
        ("m_h", P_M, P_H, reb_particle_derivative_m_h as DerivFn),
        ("m_k", P_M, P_K, reb_particle_derivative_m_k as DerivFn),
        ("m_ix", P_M, P_IX, reb_particle_derivative_m_ix as DerivFn),
        ("m_iy", P_M, P_IY, reb_particle_derivative_m_iy as DerivFn),
    ]
}

fn kep_first_table() -> Vec<(&'static str, usize, DerivFn)> {
    vec![
        ("e", K_E, reb_particle_derivative_e as DerivFn),
        ("inc", K_INC, reb_particle_derivative_inc as DerivFn),
        ("Omega", K_OMEGA_BIG, reb_particle_derivative_Omega as DerivFn),
        ("omega", K_OMEGA, reb_particle_derivative_omega as DerivFn),
        ("f", K_F, reb_particle_derivative_f as DerivFn),
    ]
}

fn kep_second_table() -> Vec<(&'static str, usize, usize, DerivFn)> {
    vec![
        ("e_e", K_E, K_E, reb_particle_derivative_e_e as DerivFn),
        ("inc_inc", K_INC, K_INC, reb_particle_derivative_inc_inc as DerivFn),
        ("Omega_Omega", K_OMEGA_BIG, K_OMEGA_BIG, reb_particle_derivative_Omega_Omega as DerivFn),
        ("omega_omega", K_OMEGA, K_OMEGA, reb_particle_derivative_omega_omega as DerivFn),
        ("f_f", K_F, K_F, reb_particle_derivative_f_f as DerivFn),
        ("a_e", K_A, K_E, reb_particle_derivative_a_e as DerivFn),
        ("a_inc", K_A, K_INC, reb_particle_derivative_a_inc as DerivFn),
        ("a_Omega", K_A, K_OMEGA_BIG, reb_particle_derivative_a_Omega as DerivFn),
        ("a_omega", K_A, K_OMEGA, reb_particle_derivative_a_omega as DerivFn),
        ("a_f", K_A, K_F, reb_particle_derivative_a_f as DerivFn),
        ("e_inc", K_E, K_INC, reb_particle_derivative_e_inc as DerivFn),
        ("e_Omega", K_E, K_OMEGA_BIG, reb_particle_derivative_e_Omega as DerivFn),
        ("e_omega", K_E, K_OMEGA, reb_particle_derivative_e_omega as DerivFn),
        ("e_f", K_E, K_F, reb_particle_derivative_e_f as DerivFn),
        ("m_e", K_M, K_E, reb_particle_derivative_m_e as DerivFn),
        ("inc_Omega", K_INC, K_OMEGA_BIG, reb_particle_derivative_inc_Omega as DerivFn),
        ("inc_omega", K_INC, K_OMEGA, reb_particle_derivative_inc_omega as DerivFn),
        ("inc_f", K_INC, K_F, reb_particle_derivative_inc_f as DerivFn),
        ("m_inc", K_M, K_INC, reb_particle_derivative_m_inc as DerivFn),
        ("omega_Omega", K_OMEGA, K_OMEGA_BIG, reb_particle_derivative_omega_Omega as DerivFn),
        ("Omega_f", K_OMEGA_BIG, K_F, reb_particle_derivative_Omega_f as DerivFn),
        ("m_Omega", K_M, K_OMEGA_BIG, reb_particle_derivative_m_Omega as DerivFn),
        ("omega_f", K_OMEGA, K_F, reb_particle_derivative_omega_f as DerivFn),
        ("m_omega", K_M, K_OMEGA, reb_particle_derivative_m_omega as DerivFn),
        ("m_f", K_M, K_F, reb_particle_derivative_m_f as DerivFn),
    ]
}

/// Orbits used for the Pal-parameterised derivatives. The first stays
/// inside the |e| < 0.3 Newton branch of reb_tools_solve_kepler_pal, the
/// second forces the reb_M_to_E branch, and the third is exactly circular.
const PAL_CASES: [(f64, f64, f64, f64, f64, f64); 3] = [
    (1.0, 0.15, 0.40, 0.70, 1.10, 2.00),
    (2.5, 0.45, 1.05, 5.30, 0.35, 4.10),
    (0.8, 0.00, 0.60, 1.90, 2.70, 3.30),
];

/// Orbits used for the Keplerian derivatives (one prograde, one
/// retrograde). Eccentricities stay well away from 0 and 1 so that the
/// perturbed evaluations remain valid ellipses.
const KEP_CASES: [(f64, f64, f64, f64, f64, f64); 2] = [
    (1.0, 0.15, 0.40, 0.70, 1.10, 2.00),
    (2.5, 0.45, 2.30, 5.30, 0.35, 4.10),
];

#[test]
fn pal_first_derivatives_match_finite_differences() {
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let table = pal_first_table();
    assert!(table.len() == 7, "pal first-derivative table lost entries");
    for &(a, e, inc, Om, om, f) in PAL_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let base = pal_base(G, primary, po);
        for &(name, idx, func) in table.iter() {
            let analytic = func(G, primary, po);
            let fd = fd_first(pal_eval, c, base, idx, FD1_STEP);
            compare_deriv(&format!("pal {} (a={}, e={})", name, a, e), analytic, fd, FD1_RTOL);
        }
    }
}

#[test]
fn pal_second_derivatives_match_finite_differences() {
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let table = pal_second_table();
    assert!(table.len() == 28, "pal second-derivative table lost entries");
    for &(a, e, inc, Om, om, f) in PAL_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let base = pal_base(G, primary, po);
        for &(name, i, j, func) in table.iter() {
            let analytic = func(G, primary, po);
            let fd = fd_second(pal_eval, c, base, i, j, FD2_STEP);
            compare_deriv(&format!("pal {} (a={}, e={})", name, a, e), analytic, fd, FD2_RTOL);
        }
    }
}

#[test]
fn kepler_first_derivatives_match_finite_differences() {
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let table = kep_first_table();
    assert!(table.len() == 5, "kepler first-derivative table lost entries");
    for &(a, e, inc, Om, om, f) in KEP_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let base = kep_base(G, primary, po);
        for &(name, idx, func) in table.iter() {
            let analytic = func(G, primary, po);
            let fd = fd_first(kep_eval, c, base, idx, FD1_STEP);
            compare_deriv(&format!("kepler {} (a={}, e={})", name, a, e), analytic, fd, FD1_RTOL);
        }
    }
}

#[test]
fn kepler_second_derivatives_match_finite_differences() {
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let table = kep_second_table();
    assert!(table.len() == 25, "kepler second-derivative table lost entries");
    for &(a, e, inc, Om, om, f) in KEP_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let base = kep_base(G, primary, po);
        for &(name, i, j, func) in table.iter() {
            let analytic = func(G, primary, po);
            let fd = fd_second(kep_eval, c, base, i, j, FD2_STEP);
            compare_deriv(&format!("kepler {} (a={}, e={})", name, a, e), analytic, fd, FD2_RTOL);
        }
    }
}

fn all_finite(p: reb_particle) -> bool {
    p.x.is_finite()
        && p.y.is_finite()
        && p.z.is_finite()
        && p.vx.is_finite()
        && p.vy.is_finite()
        && p.vz.is_finite()
}

#[test]
fn hyperbolic_orbit_derivatives_match_finite_differences() {
    // e > 1 with a < 0: reb_orbit_from_particle takes the acosh branch and
    // reb_particle_from_orbit still produces a well-defined conic, so the
    // finite differences are meaningful.
    //
    // The C derivatives.c is written for bound orbits: several routines
    // evaluate sqrt(G*M/a), sqrt(1-e^2) or sqrt(a^3) on their own rather
    // than as one non-negative quotient, so they return NaN once a < 0 and
    // e > 1. This test checks every routine that stays real, and pins down
    // exactly which routines drop out — a regression that moved that
    // boundary in either direction would be caught here.
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let po = reb_particle_from_orbit(G, primary, 1.0e-3, -2.0, 1.5, 0.55, 1.30, 0.85, 0.50);
    let o = reb_orbit_from_particle(G, po, primary);
    assert!(o.e > 1.0, "expected a hyperbolic orbit, got e = {}", o.e);
    assert!(o.a < 0.0, "expected a negative semi-major axis, got a = {}", o.a);
    let base = kep_base(G, primary, po);

    let mut skipped: Vec<&'static str> = Vec::new();
    let mut checked = 0usize;
    for &(name, idx, func) in kep_first_table().iter() {
        let analytic = func(G, primary, po);
        if !all_finite(analytic) {
            skipped.push(name);
            continue;
        }
        let fd = fd_first(kep_eval, c, base, idx, FD1_STEP);
        compare_deriv(&format!("hyperbolic {}", name), analytic, fd, FD1_RTOL);
        checked += 1;
    }
    for &(name, i, j, func) in kep_second_table().iter() {
        let analytic = func(G, primary, po);
        if !all_finite(analytic) {
            skipped.push(name);
            continue;
        }
        let fd = fd_second(kep_eval, c, base, i, j, FD2_STEP);
        compare_deriv(&format!("hyperbolic {}", name), analytic, fd, FD2_RTOL);
        checked += 1;
    }
    skipped.sort_unstable();
    let mut expected = vec![
        "e", "a_inc", "a_Omega", "a_omega", "a_f", "e_inc", "e_Omega", "e_omega", "e_f", "m_e",
    ];
    expected.sort_unstable();
    assert!(
        skipped == expected,
        "the set of routines that go non-real on a hyperbolic orbit changed: got {:?}, expected {:?}",
        skipped,
        expected
    );
    assert!(checked == 20, "expected 20 real hyperbolic derivatives, checked {}", checked);
}

// ===========================================================================
// Orbital derivatives: closed-form identities
// ===========================================================================

#[test]
fn pal_round_trip_reproduces_particle() {
    let G = 1.0;
    let primary = deriv_primary();
    for &(a, e, inc, Om, om, f) in PAL_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
        let b = pal_base(G, primary, po);
        close(&format!("pal a (e={})", e), b[P_A], a, 1.0e-13);
        let rec = reb_particle_from_pal(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        close_state(&format!("pal round trip (e={})", e), rec, po, 1.0e-12);
        // k and h are the Pal eccentricity vector components.
        close(&format!("pal |e| (e={})", e), (b[P_K] * b[P_K] + b[P_H] * b[P_H]).sqrt(), e, 1.0e-13);
        // ix^2 + iy^2 = 4 sin^2(inc/2)
        let s = 2.0 * (inc / 2.0).sin();
        close(
            &format!("pal |i| (inc={})", inc),
            (b[P_IX] * b[P_IX] + b[P_IY] * b[P_IY]).sqrt(),
            s,
            1.0e-13,
        );
    }
}

#[test]
fn derivative_a_and_m_have_closed_forms() {
    // With the remaining Pal elements held fixed, position scales like a
    // and velocity like a^-1/2; mass enters only through a^-1/2 * sqrt(mu).
    let G = 1.0;
    let primary = deriv_primary();
    for &(a0, e, inc, Om, om, f) in PAL_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a0, e, inc, Om, om, f);
        let b = pal_base(G, primary, po);
        let a = b[P_A];
        let rec = reb_particle_from_pal(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        let dr = sub(pos_of(rec), pos_of(primary));
        let dv = sub(vel_of(rec), vel_of(primary));
        let mu = G * (po.m + primary.m);
        let rs = reb_vec3d_length_squared(dr).sqrt();
        let vs = reb_vec3d_length_squared(dv).sqrt();

        let da = reb_particle_derivative_a(G, primary, po);
        close_vec(&format!("d/da position (a={})", a0), pos_of(da), reb_vec3d_mul(dr, 1.0 / a), 1.0e-12 * rs / a);
        close_vec(&format!("d/da velocity (a={})", a0), vel_of(da), reb_vec3d_mul(dv, -0.5 / a), 1.0e-12 * vs / a);

        let daa = reb_particle_derivative_a_a(G, primary, po);
        close_vec(&format!("d2/da2 position (a={})", a0), pos_of(daa), reb_vec3d::default(), 1.0e-14 * rs);
        close_vec(
            &format!("d2/da2 velocity (a={})", a0),
            vel_of(daa),
            reb_vec3d_mul(dv, 0.75 / (a * a)),
            1.0e-12 * vs / (a * a),
        );

        let dm = reb_particle_derivative_m(G, primary, po);
        bits("d/dm mass", dm.m, 1.0);
        close_vec(&format!("d/dm position (a={})", a0), pos_of(dm), reb_vec3d::default(), 1.0e-16);
        close_vec(
            &format!("d/dm velocity (a={})", a0),
            vel_of(dm),
            reb_vec3d_mul(dv, 0.5 / (mu / G)),
            1.0e-12 * vs,
        );

        let dmm = reb_particle_derivative_m_m(G, primary, po);
        close_vec(
            &format!("d2/dm2 velocity (a={})", a0),
            vel_of(dmm),
            reb_vec3d_mul(dv, -0.25 / ((mu / G) * (mu / G))),
            1.0e-12 * vs,
        );

        let dma = reb_particle_derivative_m_a(G, primary, po);
        close_vec(
            &format!("d2/dmda velocity (a={})", a0),
            vel_of(dma),
            reb_vec3d_mul(dv, -0.25 / (a * (mu / G))),
            1.0e-12 * vs / a,
        );
    }
}

#[test]
fn derivative_lambda_advances_the_orbit_in_time() {
    // Increasing the mean longitude at fixed (a, k, h, ix, iy) is exactly
    // a translation in time, so d(position)/dlambda = v/n and
    // d(velocity)/dlambda = acceleration/n with n the mean motion.
    let G = 1.0;
    let primary = deriv_primary();
    for &(a0, e, inc, Om, om, f) in PAL_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a0, e, inc, Om, om, f);
        let b = pal_base(G, primary, po);
        let a = b[P_A];
        let rec = reb_particle_from_pal(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        let dr = sub(pos_of(rec), pos_of(primary));
        let dv = sub(vel_of(rec), vel_of(primary));
        let mu = G * (po.m + primary.m);
        let n = (mu / (a * a * a)).sqrt();
        let d = reb_particle_derivative_lambda(G, primary, po);

        let vs = reb_vec3d_length_squared(dv).sqrt();
        close_vec(
            &format!("d/dlambda position (a={})", a0),
            pos_of(d),
            reb_vec3d_mul(dv, 1.0 / n),
            1.0e-12 * vs / n,
        );
        let rmag = reb_vec3d_length_squared(dr).sqrt();
        let acc = reb_vec3d_mul(dr, -mu / (rmag * rmag * rmag));
        let as_ = reb_vec3d_length_squared(acc).sqrt();
        close_vec(
            &format!("d/dlambda velocity (a={})", a0),
            vel_of(d),
            reb_vec3d_mul(acc, 1.0 / n),
            1.0e-9 * as_ / n,
        );
    }
}

#[test]
fn derivative_f_scales_with_the_angular_momentum() {
    // df/dt = h/r^2, hence d(position)/df = v r^2/h and
    // d(velocity)/df = -mu r_vec/(r h).
    let G = 1.0;
    let primary = deriv_primary();
    for &(a0, e0, inc, Om, om, f0) in KEP_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a0, e0, inc, Om, om, f0);
        let b = kep_base(G, primary, po);
        let (a, e, f) = (b[K_A], b[K_E], b[K_F]);
        let rec = reb_particle_from_orbit(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        let dr = sub(pos_of(rec), pos_of(primary));
        let dv = sub(vel_of(rec), vel_of(primary));
        let mu = G * (po.m + primary.m);
        let r = a * (1. - e * e) / (1. + e * f.cos());
        let h = (mu * a * (1. - e * e)).sqrt();
        let d = reb_particle_derivative_f(G, primary, po);

        let vs = reb_vec3d_length_squared(dv).sqrt();
        close_vec(
            &format!("d/df position (a={})", a0),
            pos_of(d),
            reb_vec3d_mul(dv, r * r / h),
            1.0e-12 * vs * r * r / h,
        );
        let want_v = reb_vec3d_mul(dr, -mu / (r * h));
        let ws = reb_vec3d_length_squared(want_v).sqrt();
        close_vec(&format!("d/df velocity (a={})", a0), vel_of(d), want_v, 1.0e-12 * ws);
    }
}

#[test]
fn derivative_Omega_is_a_rotation_about_z() {
    // Omega enters reb_particle_from_orbit only through a rotation about
    // the z axis, so the first derivative is z x state and the second is
    // z x (z x state).
    let G = 1.0;
    let primary = deriv_primary();
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    for &(a0, e0, inc, Om, om, f0) in KEP_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a0, e0, inc, Om, om, f0);
        let b = kep_base(G, primary, po);
        let rec = reb_particle_from_orbit(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        let dr = sub(pos_of(rec), pos_of(primary));
        let dv = sub(vel_of(rec), vel_of(primary));
        let rs = reb_vec3d_length_squared(dr).sqrt();
        let vs = reb_vec3d_length_squared(dv).sqrt();

        let d = reb_particle_derivative_Omega(G, primary, po);
        close_vec(&format!("d/dOmega position (a={})", a0), pos_of(d), reb_vec3d_cross(z, dr), 1.0e-13 * rs);
        close_vec(&format!("d/dOmega velocity (a={})", a0), vel_of(d), reb_vec3d_cross(z, dv), 1.0e-13 * vs);

        let dd = reb_particle_derivative_Omega_Omega(G, primary, po);
        close_vec(
            &format!("d2/dOmega2 position (a={})", a0),
            pos_of(dd),
            reb_vec3d_cross(z, reb_vec3d_cross(z, dr)),
            1.0e-13 * rs,
        );
        close_vec(
            &format!("d2/dOmega2 velocity (a={})", a0),
            vel_of(dd),
            reb_vec3d_cross(z, reb_vec3d_cross(z, dv)),
            1.0e-13 * vs,
        );
    }
}

#[test]
fn derivative_inc_and_omega_are_rotations_about_node_and_normal() {
    // inc rotates the orbit about the line of nodes (cos Omega, sin Omega, 0);
    // omega rotates it about the orbit normal (sin i sin Omega, -sin i cos Omega, cos i).
    let G = 1.0;
    let primary = deriv_primary();
    for &(a0, e0, inc0, Om0, om0, f0) in KEP_CASES.iter() {
        let po = reb_particle_from_orbit(G, primary, 1.0e-3, a0, e0, inc0, Om0, om0, f0);
        let b = kep_base(G, primary, po);
        let (inc, Om) = (b[K_INC], b[K_OMEGA_BIG]);
        let rec = reb_particle_from_orbit(G, primary, po.m, b[0], b[1], b[2], b[3], b[4], b[5]);
        let dr = sub(pos_of(rec), pos_of(primary));
        let dv = sub(vel_of(rec), vel_of(primary));
        let rs = reb_vec3d_length_squared(dr).sqrt();
        let vs = reb_vec3d_length_squared(dv).sqrt();

        let node = reb_vec3d { x: Om.cos(), y: Om.sin(), z: 0.0 };
        let di = reb_particle_derivative_inc(G, primary, po);
        close_vec(&format!("d/dinc position (a={})", a0), pos_of(di), reb_vec3d_cross(node, dr), 1.0e-13 * rs);
        close_vec(&format!("d/dinc velocity (a={})", a0), vel_of(di), reb_vec3d_cross(node, dv), 1.0e-13 * vs);

        let normal = reb_vec3d {
            x: inc.sin() * Om.sin(),
            y: -inc.sin() * Om.cos(),
            z: inc.cos(),
        };
        // Cross-check the normal against r x v of the reconstructed orbit.
        let hvec = reb_vec3d_cross(dr, dv);
        let hhat = reb_vec3d_normalize(hvec);
        close_vec(&format!("orbit normal (a={})", a0), hhat, normal, 1.0e-12);

        let dom = reb_particle_derivative_omega(G, primary, po);
        close_vec(
            &format!("d/domega position (a={})", a0),
            pos_of(dom),
            reb_vec3d_cross(normal, dr),
            1.0e-13 * rs,
        );
        close_vec(
            &format!("d/domega velocity (a={})", a0),
            vel_of(dom),
            reb_vec3d_cross(normal, dv),
            1.0e-13 * vs,
        );
    }
}

#[test]
fn second_derivatives_are_symmetric_in_their_arguments() {
    // Mixed partials must not depend on the differentiation order; compare
    // each analytic mixed derivative against the finite difference taken
    // with the two parameters swapped.
    let G = 1.0;
    let primary = deriv_primary();
    let c = Ctx { G, primary };
    let (a, e, inc, Om, om, f) = PAL_CASES[0];
    let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
    let pbase = pal_base(G, primary, po);
    for &(name, i, j, func) in pal_second_table().iter() {
        if i == j {
            continue;
        }
        let analytic = func(G, primary, po);
        let fd = fd_second(pal_eval, c, pbase, j, i, FD2_STEP);
        compare_deriv(&format!("pal {} (swapped)", name), analytic, fd, FD2_RTOL);
    }
    let (a, e, inc, Om, om, f) = KEP_CASES[0];
    let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
    let kbase = kep_base(G, primary, po);
    for &(name, i, j, func) in kep_second_table().iter() {
        if i == j {
            continue;
        }
        let analytic = func(G, primary, po);
        let fd = fd_second(kep_eval, c, kbase, j, i, FD2_STEP);
        compare_deriv(&format!("kepler {} (swapped)", name), analytic, fd, FD2_RTOL);
    }
}

#[test]
fn derivatives_are_deterministic_bit_for_bit() {
    // The derivative routines are pure; two evaluations from identical
    // inputs must agree in every bit.
    let G = 1.0;
    let primary = deriv_primary();
    let (a, e, inc, Om, om, f) = PAL_CASES[0];
    let po = reb_particle_from_orbit(G, primary, 1.0e-3, a, e, inc, Om, om, f);
    for &(name, _, func) in pal_first_table().iter() {
        let p1 = func(G, primary, po);
        let p2 = func(G, primary, po);
        bits(&format!("{} determinism x", name), p1.x, p2.x);
        bits(&format!("{} determinism vz", name), p1.vz, p2.vz);
    }
    for &(name, _, _, func) in kep_second_table().iter() {
        let p1 = func(G, primary, po);
        let p2 = func(G, primary, po);
        bits(&format!("{} determinism y", name), p1.y, p2.y);
        bits(&format!("{} determinism vx", name), p1.vx, p2.vx);
    }
}
