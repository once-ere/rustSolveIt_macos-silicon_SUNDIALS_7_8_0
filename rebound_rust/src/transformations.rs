//! transformations.rs — transformations between inertial, Jacobi,
//! democratic-heliocentric (Duncan, Levison & Lee 1998), WHDS
//! (Hernandez) and barycentric coordinates (from transformations.c;
//! serial paths).
//!
//! The C passes separate `particles`, `p_j` and `p_mass` arrays; here
//! they are slices, with the same element order and arithmetic.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Dan Tamayo and contributors. See crate root.

use crate::types::reb_particle;

// ---- Jacobi --------------------------------------------------------------

pub fn reb_transformations_inertial_to_jacobi_posvel(
    particles: &[reb_particle],
    p_j: &mut [reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_mass[0].m;
    let mut s_x = eta * particles[0].x;
    let mut s_y = eta * particles[0].y;
    let mut s_z = eta * particles[0].z;
    let mut s_vx = eta * particles[0].vx;
    let mut s_vy = eta * particles[0].vy;
    let mut s_vz = eta * particles[0].vz;
    for i in 1..N_active {
        let ei = 1. / eta;
        let pi = particles[i];
        eta += p_mass[i].m;
        let pme = eta * ei;
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x * ei;
        p_j[i].y = pi.y - s_y * ei;
        p_j[i].z = pi.z - s_z * ei;
        p_j[i].vx = pi.vx - s_vx * ei;
        p_j[i].vy = pi.vy - s_vy * ei;
        p_j[i].vz = pi.vz - s_vz * ei;
        s_x = s_x * pme + p_mass[i].m * p_j[i].x;
        s_y = s_y * pme + p_mass[i].m * p_j[i].y;
        s_z = s_z * pme + p_mass[i].m * p_j[i].z;
        s_vx = s_vx * pme + p_mass[i].m * p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m * p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m * p_j[i].vz;
    }
    let ei = 1. / eta;
    for i in N_active..N {
        let pi = particles[i];
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x * ei;
        p_j[i].y = pi.y - s_y * ei;
        p_j[i].z = pi.z - s_z * ei;
        p_j[i].vx = pi.vx - s_vx * ei;
        p_j[i].vy = pi.vy - s_vy * ei;
        p_j[i].vz = pi.vz - s_vz * ei;
    }
    let Mtotal = eta;
    let Mtotali = 1. / Mtotal;
    p_j[0].m = Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
}

pub fn reb_transformations_inertial_to_jacobi_posvelacc(
    particles: &[reb_particle],
    p_j: &mut [reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_mass[0].m;
    let mut s_x = eta * particles[0].x;
    let mut s_y = eta * particles[0].y;
    let mut s_z = eta * particles[0].z;
    let mut s_vx = eta * particles[0].vx;
    let mut s_vy = eta * particles[0].vy;
    let mut s_vz = eta * particles[0].vz;
    let mut s_ax = eta * particles[0].ax;
    let mut s_ay = eta * particles[0].ay;
    let mut s_az = eta * particles[0].az;
    for i in 1..N_active {
        let ei = 1. / eta;
        let pi = particles[i];
        eta += p_mass[i].m;
        let pme = eta * ei;
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x * ei;
        p_j[i].y = pi.y - s_y * ei;
        p_j[i].z = pi.z - s_z * ei;
        p_j[i].vx = pi.vx - s_vx * ei;
        p_j[i].vy = pi.vy - s_vy * ei;
        p_j[i].vz = pi.vz - s_vz * ei;
        p_j[i].ax = pi.ax - s_ax * ei;
        p_j[i].ay = pi.ay - s_ay * ei;
        p_j[i].az = pi.az - s_az * ei;
        s_x = s_x * pme + p_mass[i].m * p_j[i].x;
        s_y = s_y * pme + p_mass[i].m * p_j[i].y;
        s_z = s_z * pme + p_mass[i].m * p_j[i].z;
        s_vx = s_vx * pme + p_mass[i].m * p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m * p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m * p_j[i].vz;
        s_ax = s_ax * pme + p_mass[i].m * p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m * p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m * p_j[i].az;
    }
    let ei = 1. / eta;
    for i in N_active..N {
        let pi = particles[i];
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x * ei;
        p_j[i].y = pi.y - s_y * ei;
        p_j[i].z = pi.z - s_z * ei;
        p_j[i].vx = pi.vx - s_vx * ei;
        p_j[i].vy = pi.vy - s_vy * ei;
        p_j[i].vz = pi.vz - s_vz * ei;
        p_j[i].ax = pi.ax - s_ax * ei;
        p_j[i].ay = pi.ay - s_ay * ei;
        p_j[i].az = pi.az - s_az * ei;
    }
    let Mtotal = eta;
    let Mtotali = 1. / Mtotal;
    p_j[0].m = Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

pub fn reb_transformations_inertial_to_jacobi_acc(
    particles: &[reb_particle],
    p_j: &mut [reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_mass[0].m;
    let mut s_ax = eta * particles[0].ax;
    let mut s_ay = eta * particles[0].ay;
    let mut s_az = eta * particles[0].az;
    for i in 1..N_active {
        let ei = 1. / eta;
        let pi = particles[i];
        eta += p_mass[i].m;
        let pme = eta * ei;
        p_j[i].ax = pi.ax - s_ax * ei;
        p_j[i].ay = pi.ay - s_ay * ei;
        p_j[i].az = pi.az - s_az * ei;
        s_ax = s_ax * pme + p_mass[i].m * p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m * p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m * p_j[i].az;
    }
    let ei = 1. / eta;
    for i in N_active..N {
        let pi = particles[i];
        p_j[i].ax = pi.ax - s_ax * ei;
        p_j[i].ay = pi.ay - s_ay * ei;
        p_j[i].az = pi.az - s_az * ei;
    }
    let Mtotal = eta;
    let Mtotali = 1. / Mtotal;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

pub fn reb_transformations_jacobi_to_inertial_posvel(
    particles: &mut [reb_particle],
    p_j: &[reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_j[0].m;
    let mut s_x = p_j[0].x * eta;
    let mut s_y = p_j[0].y * eta;
    let mut s_z = p_j[0].z * eta;
    let mut s_vx = p_j[0].vx * eta;
    let mut s_vy = p_j[0].vy * eta;
    let mut s_vz = p_j[0].vz * eta;
    let mut i = N;
    while i > N_active {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        particles[i].x = pji.x + s_x * ei;
        particles[i].y = pji.y + s_y * ei;
        particles[i].z = pji.z + s_z * ei;
        particles[i].vx = pji.vx + s_vx * ei;
        particles[i].vy = pji.vy + s_vy * ei;
        particles[i].vz = pji.vz + s_vz * ei;
    }
    let mut i = N_active;
    while i > 1 {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        s_x = (s_x - p_mass[i].m * pji.x) * ei;
        s_y = (s_y - p_mass[i].m * pji.y) * ei;
        s_z = (s_z - p_mass[i].m * pji.z) * ei;
        s_vx = (s_vx - p_mass[i].m * pji.vx) * ei;
        s_vy = (s_vy - p_mass[i].m * pji.vy) * ei;
        s_vz = (s_vz - p_mass[i].m * pji.vz) * ei;
        particles[i].x = pji.x + s_x;
        particles[i].y = pji.y + s_y;
        particles[i].z = pji.z + s_z;
        particles[i].vx = pji.vx + s_vx;
        particles[i].vy = pji.vy + s_vy;
        particles[i].vz = pji.vz + s_vz;
        eta -= p_mass[i].m;
        s_x *= eta;
        s_y *= eta;
        s_z *= eta;
        s_vx *= eta;
        s_vy *= eta;
        s_vz *= eta;
    }
    let mi = 1. / eta;
    particles[0].x = s_x * mi;
    particles[0].y = s_y * mi;
    particles[0].z = s_z * mi;
    particles[0].vx = s_vx * mi;
    particles[0].vy = s_vy * mi;
    particles[0].vz = s_vz * mi;
}

pub fn reb_transformations_jacobi_to_inertial_pos(
    particles: &mut [reb_particle],
    p_j: &[reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_j[0].m;
    let mut s_x = p_j[0].x * eta;
    let mut s_y = p_j[0].y * eta;
    let mut s_z = p_j[0].z * eta;
    let mut i = N;
    while i > N_active {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        particles[i].x = pji.x + s_x * ei;
        particles[i].y = pji.y + s_y * ei;
        particles[i].z = pji.z + s_z * ei;
    }
    let mut i = N_active;
    while i > 1 {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        s_x = (s_x - p_mass[i].m * pji.x) * ei;
        s_y = (s_y - p_mass[i].m * pji.y) * ei;
        s_z = (s_z - p_mass[i].m * pji.z) * ei;
        particles[i].x = pji.x + s_x;
        particles[i].y = pji.y + s_y;
        particles[i].z = pji.z + s_z;
        eta -= p_mass[i].m;
        s_x *= eta;
        s_y *= eta;
        s_z *= eta;
    }
    let mi = 1. / eta;
    particles[0].x = s_x * mi;
    particles[0].y = s_y * mi;
    particles[0].z = s_z * mi;
}

pub fn reb_transformations_jacobi_to_inertial_acc(
    particles: &mut [reb_particle],
    p_j: &[reb_particle],
    p_mass: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut eta = p_j[0].m;
    let mut s_ax = p_j[0].ax * eta;
    let mut s_ay = p_j[0].ay * eta;
    let mut s_az = p_j[0].az * eta;
    let mut i = N;
    while i > N_active {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        particles[i].ax = pji.ax + s_ax * ei;
        particles[i].ay = pji.ay + s_ay * ei;
        particles[i].az = pji.az + s_az * ei;
    }
    let mut i = N_active;
    while i > 1 {
        i -= 1;
        let pji = p_j[i];
        let ei = 1. / eta;
        s_ax = (s_ax - p_mass[i].m * pji.ax) * ei;
        s_ay = (s_ay - p_mass[i].m * pji.ay) * ei;
        s_az = (s_az - p_mass[i].m * pji.az) * ei;
        particles[i].ax = pji.ax + s_ax;
        particles[i].ay = pji.ay + s_ay;
        particles[i].az = pji.az + s_az;
        eta -= p_mass[i].m;
        s_ax *= eta;
        s_ay *= eta;
        s_az *= eta;
    }
    let mi = 1. / eta;
    particles[0].ax = s_ax * mi;
    particles[0].ay = s_ay * mi;
    particles[0].az = s_az * mi;
}

// ---- WHDS (Hernandez) ----------------------------------------------------

pub fn reb_transformations_inertial_to_whds_posvel(
    particles: &[reb_particle],
    p_h: &mut [reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut x0 = 0.;
    let mut y0 = 0.;
    let mut z0 = 0.;
    let mut vx0 = 0.;
    let mut vy0 = 0.;
    let mut vz0 = 0.;
    let mut m0 = 0.;
    for i in 0..N_active {
        let m = particles[i].m;
        x0 += particles[i].x * m;
        y0 += particles[i].y * m;
        z0 += particles[i].z * m;
        vx0 += particles[i].vx * m;
        vy0 += particles[i].vy * m;
        vz0 += particles[i].vz * m;
        m0 += m;
    }
    p_h[0].x = x0 / m0;
    p_h[0].y = y0 / m0;
    p_h[0].z = z0 / m0;
    p_h[0].vx = vx0 / m0;
    p_h[0].vy = vy0 / m0;
    p_h[0].vz = vz0 / m0;
    p_h[0].m = m0;

    let m0 = particles[0].m;
    for i in 1..N_active {
        p_h[i].x = particles[i].x - particles[0].x;
        p_h[i].y = particles[i].y - particles[0].y;
        p_h[i].z = particles[i].z - particles[0].z;
        let mi = particles[i].m;
        let mf = (m0 + mi) / m0;
        p_h[i].vx = mf * (particles[i].vx - p_h[0].vx);
        p_h[i].vy = mf * (particles[i].vy - p_h[0].vy);
        p_h[i].vz = mf * (particles[i].vz - p_h[0].vz);
        p_h[i].m = mi;
    }
    for i in N_active..N {
        p_h[i].x = particles[i].x - particles[0].x;
        p_h[i].y = particles[i].y - particles[0].y;
        p_h[i].z = particles[i].z - particles[0].z;
        p_h[i].vx = particles[i].vx - p_h[0].vx;
        p_h[i].vy = particles[i].vy - p_h[0].vy;
        p_h[i].vz = particles[i].vz - p_h[0].vz;
        p_h[i].m = particles[i].m;
    }
}

pub fn reb_transformations_whds_to_inertial_pos(
    particles: &mut [reb_particle],
    p_h: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    // Same as in heliocentric case.
    reb_transformations_democraticheliocentric_to_inertial_pos(particles, p_h, N, N_active);
}

pub fn reb_transformations_whds_to_inertial_posvel(
    particles: &mut [reb_particle],
    p_h: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    reb_transformations_whds_to_inertial_pos(particles, p_h, N, N_active);
    let m0 = particles[0].m;
    for i in 1..N_active {
        let mi = particles[i].m;
        let mf = (m0 + mi) / m0;
        particles[i].vx = p_h[i].vx / mf + p_h[0].vx;
        particles[i].vy = p_h[i].vy / mf + p_h[0].vy;
        particles[i].vz = p_h[i].vz / mf + p_h[0].vz;
    }
    for i in N_active..N {
        particles[i].vx = p_h[i].vx + p_h[0].vx;
        particles[i].vy = p_h[i].vy + p_h[0].vy;
        particles[i].vz = p_h[i].vz + p_h[0].vz;
    }
    let mut vx0 = 0.;
    let mut vy0 = 0.;
    let mut vz0 = 0.;
    for i in 1..N_active {
        let m = particles[i].m;
        vx0 += p_h[i].vx * m / (m0 + m);
        vy0 += p_h[i].vy * m / (m0 + m);
        vz0 += p_h[i].vz * m / (m0 + m);
    }
    particles[0].vx = p_h[0].vx - vx0;
    particles[0].vy = p_h[0].vy - vy0;
    particles[0].vz = p_h[0].vz - vz0;
}

// ---- Democratic heliocentric (Duncan, Levison & Lee 1998) ---------------

pub fn reb_transformations_inertial_to_democraticheliocentric_posvel(
    particles: &[reb_particle],
    p_h: &mut [reb_particle],
    N: usize,
    N_active: usize,
) {
    let mut x0 = 0.;
    let mut y0 = 0.;
    let mut z0 = 0.;
    let mut vx0 = 0.;
    let mut vy0 = 0.;
    let mut vz0 = 0.;
    let mut m0 = 0.;
    for i in 0..N_active {
        let m = particles[i].m;
        x0 += particles[i].x * m;
        y0 += particles[i].y * m;
        z0 += particles[i].z * m;
        vx0 += particles[i].vx * m;
        vy0 += particles[i].vy * m;
        vz0 += particles[i].vz * m;
        m0 += m;
    }
    p_h[0].x = x0 / m0;
    p_h[0].y = y0 / m0;
    p_h[0].z = z0 / m0;
    p_h[0].vx = vx0 / m0;
    p_h[0].vy = vy0 / m0;
    p_h[0].vz = vz0 / m0;
    p_h[0].m = m0;

    for i in 1..N {
        p_h[i].x = particles[i].x - particles[0].x;
        p_h[i].y = particles[i].y - particles[0].y;
        p_h[i].z = particles[i].z - particles[0].z;
        p_h[i].vx = particles[i].vx - p_h[0].vx;
        p_h[i].vy = particles[i].vy - p_h[0].vy;
        p_h[i].vz = particles[i].vz - p_h[0].vz;
        p_h[i].m = particles[i].m;
    }
}

pub fn reb_transformations_democraticheliocentric_to_inertial_pos(
    particles: &mut [reb_particle],
    p_h: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    let mtot = p_h[0].m;
    let mut x0 = 0.;
    let mut y0 = 0.;
    let mut z0 = 0.;
    for i in 1..N_active {
        let m = p_h[i].m;
        x0 += p_h[i].x * m / mtot;
        y0 += p_h[i].y * m / mtot;
        z0 += p_h[i].z * m / mtot;
        particles[i].m = m; // in case of merger/mass change
    }
    particles[0].x = p_h[0].x - x0;
    particles[0].y = p_h[0].y - y0;
    particles[0].z = p_h[0].z - z0;
    for i in 1..N {
        particles[i].x = p_h[i].x + particles[0].x;
        particles[i].y = p_h[i].y + particles[0].y;
        particles[i].z = p_h[i].z + particles[0].z;
    }
}

pub fn reb_transformations_democraticheliocentric_to_inertial_posvel(
    particles: &mut [reb_particle],
    p_h: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    reb_transformations_democraticheliocentric_to_inertial_pos(particles, p_h, N, N_active);
    let m0 = particles[0].m;
    for i in 1..N {
        particles[i].vx = p_h[i].vx + p_h[0].vx;
        particles[i].vy = p_h[i].vy + p_h[0].vy;
        particles[i].vz = p_h[i].vz + p_h[0].vz;
    }
    let mut vx0 = 0.;
    let mut vy0 = 0.;
    let mut vz0 = 0.;
    for i in 1..N_active {
        let m = particles[i].m;
        vx0 += p_h[i].vx * m / m0;
        vy0 += p_h[i].vy * m / m0;
        vz0 += p_h[i].vz * m / m0;
    }
    particles[0].vx = p_h[0].vx - vx0;
    particles[0].vy = p_h[0].vy - vy0;
    particles[0].vz = p_h[0].vz - vz0;
}

// ---- Barycentric ---------------------------------------------------------

pub fn reb_transformations_barycentric_to_inertial_pos(
    particles: &mut [reb_particle],
    p_b: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    particles[0].x = p_b[0].m * p_b[0].x;
    particles[0].y = p_b[0].m * p_b[0].y;
    particles[0].z = p_b[0].m * p_b[0].z;
    particles[0].m = p_b[0].m;

    let mut s_x = 0.0;
    let mut s_y = 0.0;
    let mut s_z = 0.0;
    let mut s_m = 0.0;
    for i in 1..N {
        particles[i].x = p_b[i].x + p_b[0].x;
        particles[i].y = p_b[i].y + p_b[0].y;
        particles[i].z = p_b[i].z + p_b[0].z;
        if i < N_active {
            let m = p_b[i].m;
            particles[i].m = m;
            s_x += particles[i].x * m;
            s_y += particles[i].y * m;
            s_z += particles[i].z * m;
            s_m += m;
        }
    }
    particles[0].x -= s_x;
    particles[0].y -= s_y;
    particles[0].z -= s_z;
    particles[0].m -= s_m;
    let m0i = 1. / particles[0].m;
    particles[0].x *= m0i;
    particles[0].y *= m0i;
    particles[0].z *= m0i;
}

pub fn reb_transformations_barycentric_to_inertial_acc(
    particles: &mut [reb_particle],
    p_b: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    particles[0].ax = p_b[0].m * p_b[0].ax;
    particles[0].ay = p_b[0].m * p_b[0].ay;
    particles[0].az = p_b[0].m * p_b[0].az;
    particles[0].m = p_b[0].m;

    let mut s_ax = 0.0;
    let mut s_ay = 0.0;
    let mut s_az = 0.0;
    let mut s_m = 0.0;
    for i in 1..N {
        particles[i].ax = p_b[i].ax + p_b[0].ax;
        particles[i].ay = p_b[i].ay + p_b[0].ay;
        particles[i].az = p_b[i].az + p_b[0].az;
        if i < N_active {
            let m = p_b[i].m;
            particles[i].m = m;
            s_ax += particles[i].ax * m;
            s_ay += particles[i].ay * m;
            s_az += particles[i].az * m;
            s_m += m;
        }
    }
    particles[0].ax -= s_ax;
    particles[0].ay -= s_ay;
    particles[0].az -= s_az;
    particles[0].m -= s_m;
    let m0i = 1. / particles[0].m;
    particles[0].ax *= m0i;
    particles[0].ay *= m0i;
    particles[0].az *= m0i;
}

pub fn reb_transformations_barycentric_to_inertial_posvel(
    particles: &mut [reb_particle],
    p_b: &[reb_particle],
    N: usize,
    N_active: usize,
) {
    particles[0].x = p_b[0].m * p_b[0].x;
    particles[0].y = p_b[0].m * p_b[0].y;
    particles[0].z = p_b[0].m * p_b[0].z;
    particles[0].vx = p_b[0].m * p_b[0].vx;
    particles[0].vy = p_b[0].m * p_b[0].vy;
    particles[0].vz = p_b[0].m * p_b[0].vz;
    particles[0].m = p_b[0].m;

    let mut s_x = 0.0;
    let mut s_y = 0.0;
    let mut s_z = 0.0;
    let mut s_vx = 0.0;
    let mut s_vy = 0.0;
    let mut s_vz = 0.0;
    let mut s_m = 0.0;
    for i in 1..N {
        particles[i].x = p_b[i].x + p_b[0].x;
        particles[i].y = p_b[i].y + p_b[0].y;
        particles[i].z = p_b[i].z + p_b[0].z;
        particles[i].vx = p_b[i].vx + p_b[0].vx;
        particles[i].vy = p_b[i].vy + p_b[0].vy;
        particles[i].vz = p_b[i].vz + p_b[0].vz;
        if i < N_active {
            let m = p_b[i].m;
            particles[i].m = m;
            s_x += particles[i].x * m;
            s_y += particles[i].y * m;
            s_z += particles[i].z * m;
            s_vx += particles[i].vx * m;
            s_vy += particles[i].vy * m;
            s_vz += particles[i].vz * m;
            s_m += m;
        }
    }
    particles[0].x -= s_x;
    particles[0].y -= s_y;
    particles[0].z -= s_z;
    particles[0].vx -= s_vx;
    particles[0].vy -= s_vy;
    particles[0].vz -= s_vz;
    particles[0].m -= s_m;
    let m0i = 1. / particles[0].m;
    particles[0].x *= m0i;
    particles[0].y *= m0i;
    particles[0].z *= m0i;
    particles[0].vx *= m0i;
    particles[0].vy *= m0i;
    particles[0].vz *= m0i;
}

pub fn reb_transformations_inertial_to_barycentric_posvel(
    particles: &[reb_particle],
    p_b: &mut [reb_particle],
    N: usize,
    N_active: usize,
) {
    p_b[0].x = particles[0].m * particles[0].x;
    p_b[0].y = particles[0].m * particles[0].y;
    p_b[0].z = particles[0].m * particles[0].z;
    p_b[0].vx = particles[0].m * particles[0].vx;
    p_b[0].vy = particles[0].m * particles[0].vy;
    p_b[0].vz = particles[0].m * particles[0].vz;
    p_b[0].m = particles[0].m;

    let mut s_x = 0.0;
    let mut s_y = 0.0;
    let mut s_z = 0.0;
    let mut s_vx = 0.0;
    let mut s_vy = 0.0;
    let mut s_vz = 0.0;
    let mut s_m = 0.0;
    for i in 1..N_active {
        let m = particles[i].m;
        p_b[i].m = m;
        s_x += particles[i].x * m;
        s_y += particles[i].y * m;
        s_z += particles[i].z * m;
        s_vx += particles[i].vx * m;
        s_vy += particles[i].vy * m;
        s_vz += particles[i].vz * m;
        s_m += m;
    }
    p_b[0].x += s_x;
    p_b[0].y += s_y;
    p_b[0].z += s_z;
    p_b[0].vx += s_vx;
    p_b[0].vy += s_vy;
    p_b[0].vz += s_vz;
    p_b[0].m += s_m;
    let mi = 1.0 / p_b[0].m;
    p_b[0].x *= mi;
    p_b[0].y *= mi;
    p_b[0].z *= mi;
    p_b[0].vx *= mi;
    p_b[0].vy *= mi;
    p_b[0].vz *= mi;
    for i in 1..N {
        p_b[i].x = particles[i].x - p_b[0].x;
        p_b[i].y = particles[i].y - p_b[0].y;
        p_b[i].z = particles[i].z - p_b[0].z;
        p_b[i].vx = particles[i].vx - p_b[0].vx;
        p_b[i].vy = particles[i].vy - p_b[0].vy;
        p_b[i].vz = particles[i].vz - p_b[0].vz;
    }
}
