//! rotations.rs — vector manipulation, quaternion rotations and the
//! single-precision matrix helpers (from rotations.c; conventions of
//! the Apple SIMD quaternion framework).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Dan Tamayo and contributors. See crate root.

use crate::types::*;
use std::f64::consts::PI as M_PI;

// ---- reb_vec3d manipulation functions ------------------------------------

pub fn reb_vec3d_mul(v: reb_vec3d, s: f64) -> reb_vec3d {
    reb_vec3d { x: s * v.x, y: s * v.y, z: s * v.z }
}

pub fn reb_vec3d_add(v: reb_vec3d, w: reb_vec3d) -> reb_vec3d {
    reb_vec3d { x: v.x + w.x, y: v.y + w.y, z: v.z + w.z }
}

pub fn reb_vec3d_cross(a: reb_vec3d, b: reb_vec3d) -> reb_vec3d {
    reb_vec3d {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}

pub fn reb_vec3d_dot(a: reb_vec3d, b: reb_vec3d) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

pub fn reb_vec3d_length_squared(v: reb_vec3d) -> f64 {
    reb_vec3d_dot(v, v)
}

pub fn reb_vec3d_normalize(v: reb_vec3d) -> reb_vec3d {
    reb_vec3d_mul(v, 1. / reb_vec3d_length_squared(v).sqrt())
}

// ---- reb_rotation manipulation functions ---------------------------------

fn reb_rotation_imag(q: reb_rotation) -> reb_vec3d {
    reb_vec3d { x: q.ix, y: q.iy, z: q.iz }
}

pub fn reb_rotation_mul(p: reb_rotation, q: reb_rotation) -> reb_rotation {
    // v_rot = p * ( q * v)
    reb_rotation {
        r: p.r * q.r - p.ix * q.ix - p.iy * q.iy - p.iz * q.iz,
        ix: p.r * q.ix + p.ix * q.r + p.iy * q.iz - p.iz * q.iy,
        iy: p.r * q.iy - p.ix * q.iz + p.iy * q.r + p.iz * q.ix,
        iz: p.r * q.iz + p.ix * q.iy - p.iy * q.ix + p.iz * q.r,
    }
}

fn reb_rotation_length_squared(q: reb_rotation) -> f64 {
    q.r * q.r + q.ix * q.ix + q.iy * q.iy + q.iz * q.iz
}

pub fn reb_rotation_conjugate(q: reb_rotation) -> reb_rotation {
    reb_rotation { ix: -q.ix, iy: -q.iy, iz: -q.iz, r: q.r }
}

pub fn reb_rotation_normalize(q: reb_rotation) -> reb_rotation {
    let l = 1. / reb_rotation_length_squared(q).sqrt();
    reb_rotation { ix: q.ix * l, iy: q.iy * l, iz: q.iz * l, r: q.r * l }
}

pub fn reb_rotation_inverse(q: reb_rotation) -> reb_rotation {
    let mut c = reb_rotation_conjugate(q);
    let rl2 = 1. / reb_rotation_length_squared(q);
    c.r *= rl2;
    c.ix *= rl2;
    c.iy *= rl2;
    c.iz *= rl2;
    c
}

// ---- Object rotation functions -------------------------------------------

pub fn reb_vec3d_rotate(v: reb_vec3d, q: reb_rotation) -> reb_vec3d {
    let mut r = v;
    reb_vec3d_irotate(&mut r, q);
    r
}

pub fn reb_vec3d_irotate(v: &mut reb_vec3d, q: reb_rotation) {
    let imag = reb_rotation_imag(q);
    let t = reb_vec3d_mul(reb_vec3d_cross(imag, *v), 2.);
    let res = reb_vec3d_add(*v, reb_vec3d_add(reb_vec3d_mul(t, q.r), reb_vec3d_cross(imag, t)));
    v.x = res.x;
    v.y = res.y;
    v.z = res.z;
}

pub fn reb_particle_irotate(p: &mut reb_particle, q: reb_rotation) {
    let mut pos = reb_vec3d { x: p.x, y: p.y, z: p.z };
    reb_vec3d_irotate(&mut pos, q);
    p.x = pos.x;
    p.y = pos.y;
    p.z = pos.z;
    let mut vel = reb_vec3d { x: p.vx, y: p.vy, z: p.vz };
    reb_vec3d_irotate(&mut vel, q);
    p.vx = vel.x;
    p.vy = vel.y;
    p.vz = vel.z;
}

pub fn reb_simulation_irotate(sim: &mut reb_simulation, q: reb_rotation) {
    let N = sim.N;
    for i in 0..N {
        reb_particle_irotate(&mut sim.particles[i], q);
    }
}

// ---- Alternate ways of initializing a rotation ---------------------------

pub fn reb_rotation_identity() -> reb_rotation {
    reb_rotation { ix: 0.0, iy: 0.0, iz: 0.0, r: 1.0 }
}

#[inline]
fn reb_rotation_init_from_to_reduced(from: reb_vec3d, to: reb_vec3d) -> reb_rotation {
    // Internal use only
    let mut half = reb_vec3d { x: from.x + to.x, y: from.y + to.y, z: from.z + to.z };
    half = reb_vec3d_normalize(half);
    let cross = reb_vec3d_cross(from, half);
    let dot = reb_vec3d_dot(from, half);
    reb_rotation { ix: cross.x, iy: cross.y, iz: cross.z, r: dot }
}

pub fn reb_rotation_init_from_to(from: reb_vec3d, to: reb_vec3d) -> reb_rotation {
    let from = reb_vec3d_normalize(from);
    let to = reb_vec3d_normalize(to);

    if reb_vec3d_dot(from, to) >= 0. {
        // small angle
        return reb_rotation_init_from_to_reduced(from, to);
    }

    // More than 90 degrees apart: rotate in two stages (from->half), (half->to)
    let mut half = reb_vec3d { x: from.x + to.x, y: from.y + to.y, z: from.z + to.z };
    half = reb_vec3d_normalize(half);

    if !reb_vec3d_length_squared(half).is_normal() {
        // from and to point in nearly opposite directions; rotation is
        // numerically underspecified. Pick an orthogonal axis, angle pi.
        let abs_from = reb_vec3d { x: from.x.abs(), y: from.y.abs(), z: from.z.abs() };
        if abs_from.x <= abs_from.y && abs_from.x <= abs_from.z {
            let axis = reb_vec3d_cross(from, reb_vec3d { x: 1., y: 0., z: 0. });
            return reb_rotation { ix: axis.x, iy: axis.y, iz: axis.z, r: 0.0 };
        }
        if abs_from.y <= abs_from.z {
            let axis = reb_vec3d_cross(from, reb_vec3d { x: 0., y: 1., z: 0. });
            return reb_rotation { ix: axis.x, iy: axis.y, iz: axis.z, r: 0.0 };
        }
        let axis = reb_vec3d_cross(from, reb_vec3d { x: 0., y: 0., z: 1. });
        return reb_rotation { ix: axis.x, iy: axis.y, iz: axis.z, r: 0.0 };
    }

    reb_rotation_mul(
        reb_rotation_init_from_to_reduced(from, half),
        reb_rotation_init_from_to_reduced(half, to),
    )
}

pub fn reb_rotation_init_angle_axis(angle: f64, axis: reb_vec3d) -> reb_rotation {
    let axis = reb_vec3d_normalize(axis);
    let cos2 = (angle / 2.0).cos();
    let sin2 = (angle / 2.0).sin();
    let imag = reb_vec3d_mul(axis, sin2);
    reb_rotation { ix: imag.x, iy: imag.y, iz: imag.z, r: cos2 }
}

pub fn reb_rotation_init_to_new_axes(newz: reb_vec3d, newx: reb_vec3d) -> reb_rotation {
    let dotprod = reb_vec3d_dot(newz, newx);
    let newz = reb_vec3d_normalize(newz);
    // orthogonalize: newx = newx - (newx dot newz) newzhat
    let mut newx = reb_vec3d_add(newx, reb_vec3d_mul(newz, -dotprod));
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let q1 = reb_rotation_init_from_to(newz, z);
    let x = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    reb_vec3d_irotate(&mut newx, q1); // rotate newx to post-first-rotation frame
    let q2 = reb_rotation_init_from_to(newx, x);
    reb_rotation_mul(q2, q1)
}

pub fn reb_rotation_init_orbit(Omega: f64, inc: f64, omega: f64) -> reb_rotation {
    // Murray and Dermott Eq. 2.121 (left hand side)
    let x = reb_vec3d { x: 1.0, y: 0.0, z: 0.0 };
    let z = reb_vec3d { x: 0.0, y: 0.0, z: 1.0 };
    let P1 = reb_rotation_init_angle_axis(omega, z);
    let P2 = reb_rotation_init_angle_axis(inc, x);
    let P3 = reb_rotation_init_angle_axis(Omega, z);
    reb_rotation_mul(P3, reb_rotation_mul(P2, P1))
}

const MIN_INC: f64 = 1.0e-8;

pub fn reb_rotation_to_orbital(q: reb_rotation, Omega: &mut f64, inc: &mut f64, omega: &mut f64) {
    // Bernardes & Viollet (2022); angles may land outside the usual quadrant.
    let ap = q.r;
    let bp = q.iz;
    let cp = q.ix;
    let dp = q.iy;
    *inc = (2.0 * (ap * ap + bp * bp) - 1.0).acos();
    let safe1 = inc.abs() > MIN_INC;
    let safe2 = (*inc - M_PI).abs() > MIN_INC;

    if safe1 && safe2 {
        let half_sum = bp.atan2(ap);
        let half_diff = dp.atan2(cp);
        *omega = half_sum - half_diff;
        *Omega = half_sum + half_diff;
    } else {
        *Omega = 0.;
        if !safe1 {
            let half_sum = bp.atan2(ap);
            *omega = 2.0 * half_sum;
        } else {
            let half_diff = dp.atan2(cp);
            *omega = 2.0 * half_diff;
        }
    }
    if *omega < 0. {
        *omega += M_PI * 2.0;
    }
    if *Omega < 0. {
        *Omega += M_PI * 2.0;
    }
}

const QUATERNION_EPS: f64 = 1e-4; // enough for visualizations

pub fn reb_rotation_slerp(q1: reb_rotation, q2: reb_rotation, t: f64) -> reb_rotation {
    let cosHalfTheta = q1.r * q2.r + q1.ix * q2.ix + q1.iy * q2.iy + q1.iz * q2.iz;

    // if q1=q2 or q1=-q2 then theta = 0 and we can return q1
    if cosHalfTheta.abs() >= 1.0 {
        return q1;
    }

    let halfTheta = cosHalfTheta.acos();
    let sinHalfTheta = (1.0 - cosHalfTheta * cosHalfTheta).sqrt();
    let mut result = reb_rotation::default();
    // theta = 180 degrees: result not fully defined
    if sinHalfTheta.abs() < QUATERNION_EPS {
        result.r = q1.r * 0.5 + q2.r * 0.5;
        result.ix = q1.ix * 0.5 + q2.ix * 0.5;
        result.iy = q1.iy * 0.5 + q2.iy * 0.5;
        result.iz = q1.iz * 0.5 + q2.iz * 0.5;
    } else {
        let ratioA = ((1. - t) * halfTheta).sin() / sinHalfTheta;
        let ratioB = (t * halfTheta).sin() / sinHalfTheta;
        result.r = q1.r * ratioA + q2.r * ratioB;
        result.ix = q1.ix * ratioA + q2.ix * ratioB;
        result.iy = q1.iy * ratioA + q2.iy * ratioB;
        result.iz = q1.iz * ratioA + q2.iz * ratioB;
    }
    result
}

// ---- Matrix methods (single precision, used for visualization) -----------

/// rebound.h `struct reb_vec3df`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct reb_vec3df {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

/// rebound.h `struct reb_mat4df`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct reb_mat4df {
    pub m: [f32; 16],
}

pub fn reb_mat4df_identity() -> reb_mat4df {
    reb_mat4df {
        m: [
            1., 0., 0., 0., //
            0., 1., 0., 0., //
            0., 0., 1., 0., //
            0., 0., 0., 1.,
        ],
    }
}

pub fn reb_mat4df_scale(m: reb_mat4df, x: f32, y: f32, z: f32) -> reb_mat4df {
    let mut nm = m;
    for j in 0..4 {
        nm.m[0 + j * 4] *= x;
        nm.m[1 + j * 4] *= y;
        nm.m[2 + j * 4] *= z;
    }
    nm
}

pub fn reb_mat4df_eq(A: reb_mat4df, B: reb_mat4df) -> bool {
    for j in 0..16 {
        if A.m[j] != B.m[j] {
            return false;
        }
    }
    true
}

pub fn reb_mat4df_get_scale(m: reb_mat4df) -> reb_vec3df {
    reb_vec3df {
        x: (m.m[0] * m.m[0] + m.m[4] * m.m[4] + m.m[8] * m.m[8]).sqrt(),
        y: (m.m[1] * m.m[1] + m.m[5] * m.m[5] + m.m[9] * m.m[9]).sqrt(),
        z: (m.m[2] * m.m[2] + m.m[6] * m.m[6] + m.m[10] * m.m[10]).sqrt(),
    }
}

pub fn reb_mat4df_translate(m: reb_mat4df, x: f32, y: f32, z: f32) -> reb_mat4df {
    let mut nm = m;
    nm.m[3 + 4 * 0] += x * m.m[0 + 4 * 0] + y * m.m[1 + 4 * 0] + z * m.m[2 + 4 * 0];
    nm.m[3 + 4 * 1] += x * m.m[0 + 4 * 1] + y * m.m[1 + 4 * 1] + z * m.m[2 + 4 * 1];
    nm.m[3 + 4 * 2] += x * m.m[0 + 4 * 2] + y * m.m[1 + 4 * 2] + z * m.m[2 + 4 * 2];
    nm
}

pub fn reb_mat4df_multiply(A: reb_mat4df, B: reb_mat4df) -> reb_mat4df {
    let mut C = reb_mat4df { m: [0.0; 16] };
    for i in 0..4 {
        for j in 0..4 {
            C.m[i + 4 * j] = 0.;
            for k in 0..4 {
                C.m[i + 4 * j] += A.m[k + 4 * j] * B.m[i + 4 * k];
            }
        }
    }
    C
}

pub fn reb_rotation_to_mat4df(A: reb_rotation) -> reb_mat4df {
    let mut m = reb_mat4df { m: [0.0; 16] };
    let xx = (A.ix * A.ix) as f32;
    let xy = (A.ix * A.iy) as f32;
    let xz = (A.ix * A.iz) as f32;
    let xw = (A.ix * A.r) as f32;
    let yy = (A.iy * A.iy) as f32;
    let yz = (A.iy * A.iz) as f32;
    let yw = (A.iy * A.r) as f32;
    let zz = (A.iz * A.iz) as f32;
    let zw = (A.iz * A.r) as f32;
    m.m[0] = 1. - 2. * (yy + zz);
    m.m[1] = 2. * (xy - zw);
    m.m[2] = 2. * (xz + yw);
    m.m[4] = 2. * (xy + zw);
    m.m[5] = 1. - 2. * (xx + zz);
    m.m[6] = 2. * (yz - xw);
    m.m[8] = 2. * (xz - yw);
    m.m[9] = 2. * (yz + xw);
    m.m[10] = 1. - 2. * (xx + yy);
    m.m[3] = 0.;
    m.m[7] = 0.;
    m.m[11] = 0.;
    m.m[12] = 0.;
    m.m[13] = 0.;
    m.m[14] = 0.;
    m.m[15] = 1.;
    m
}

pub fn reb_mat4df_ortho(l: f32, r: f32, b: f32, t: f32, n: f32, f: f32) -> reb_mat4df {
    let mut m = reb_mat4df { m: [0.0; 16] };
    m.m[0] = 2. / (r - l);
    m.m[1] = 0.;
    m.m[2] = 0.;
    m.m[3] = -(r + l) / (r - l);
    m.m[4] = 0.;
    m.m[5] = 2. / (t - b);
    m.m[6] = 0.;
    m.m[7] = -(t + b) / (t - b);
    m.m[8] = 0.;
    m.m[9] = 0.;
    m.m[10] = -2. / (f - n);
    m.m[11] = -(f + n) / (f - n);
    m.m[12] = 0.;
    m.m[13] = 0.;
    m.m[14] = 0.;
    m.m[15] = 1.;
    m
}
