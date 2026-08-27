//! gravitational_harmonics.rs — translation of REBOUNDx gravitational_harmonics.c
//! Adds azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Gravity Fields
//!
//! ```text
//! ======================= ===============================================
//! Authors                 M. Broz
//! Implementation Paper    Tamayo, Rein, Shi and Hernandez, 2019
//! Based on                None
//! C Example               :ref:`c_example_j2`
//! Python Example          J2.ipynb
//! ======================= ===============================================
//! ```
//!
//! Allows the user to add azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation.
//! These interact with all other bodies in the simulation (treated as point masses).
//! The implementation allows the user to specify an arbitrary spin axis orientation for each oblate body, which defines the axis of symmetry.
//! This is specified through the angular rotation rate vector Omega.
//! The rotation rate Omega is not currently used other than to specify the spin axis orientation.
//! In particular, the current implementation applies the appropriate torque from the body's oblateness to the orbits of all the other planets, but does not account for the equal and opposite torque on the body's spin angular momentum.
//! The bodies spins therefore remain constant in the current implementation.
//! This is a good approximation in the limit where the bodies' spin angular momenta are much greater than the orbital angular momenta involved.
//!
//! **Effect Parameters**
//!
//! None
//!
//! **Particle Parameters**
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! J2 (double)                  No          J2 coefficient
//! J4 (double)                  No          J4 coefficient
//! R_eq (double)                No          Equatorial radius of nonspherical body used for calculating Jn harmonics
//! Omega (reb_vec3d)            No          Angular rotation frequency (Omega_x, Omega_y, Omega_z)
//! ============================ =========== ==================================================================
//! ```

use crate::core::{rebx_get_param_double, rebx_get_param_vec3d};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_particle, reb_simulation, reb_vec3d};

/// C: `#define DEFAULTOMEGA {0.0, 0.0, 1.0}`
const DEFAULTOMEGA: reb_vec3d = reb_vec3d {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};

/// gravitational_harmonics.c `j2_func`.
///
/// `J2` is `Option` because the C tests its pointer against NULL and
/// returns without touching the accelerations. `R_eq` arrives already
/// dereferenced: every call site has verified it is non-NULL (the C
/// `continue`s otherwise), so the C never reaches this function with a
/// NULL `R_eq`.
fn j2_func(
    G: f64,
    m: f64,
    J2: Option<f64>,
    R_eq: f64,
    r: f64,
    r2: f64,
    costheta2: f64,
    du: f64,
    dv: f64,
    dw: f64,
    au: &mut f64,
    av: &mut f64,
    aw: &mut f64,
) {
    let J2 = match J2 {
        None => {
            return;
        }
        Some(J2) => J2,
    };
    if J2 == 0.0 {
        return;
    }

    let f1 = 3.0 / 2.0 * G * m * (J2) * (R_eq) * (R_eq) / r2 / r2 / r;
    let f2 = 5.0 * costheta2 - 1.0;
    let f3 = f2 - 2.0;

    *au += f1 * f2 * du;
    *av += f1 * f2 * dv;
    *aw += f1 * f3 * dw;
}

/// gravitational_harmonics.c `j4_func`. See [`j2_func`] for the
/// `Option`/`f64` argument convention.
fn j4_func(
    G: f64,
    m: f64,
    J4: Option<f64>,
    R_eq: f64,
    r: f64,
    r2: f64,
    costheta2: f64,
    du: f64,
    dv: f64,
    dw: f64,
    au: &mut f64,
    av: &mut f64,
    aw: &mut f64,
) {
    let J4 = match J4 {
        None => {
            return;
        }
        Some(J4) => J4,
    };
    if J4 == 0.0 {
        return;
    }

    let f1 = 5.0 / 8.0 * G * m * (J4) * (R_eq) * (R_eq) * (R_eq) * (R_eq) / r2 / r2 / r2 / r;
    let f2 = 63.0 * costheta2 * costheta2 - 42.0 * costheta2 + 3.0;
    let f3 = f2 - 28.0 * costheta2 + 12.0;

    *au += f1 * f2 * du;
    *av += f1 * f2 * dv;
    *aw += f1 * f3 * dw;
}

/// gravitational_harmonics.c `uvw`.
///
/// Builds the body-fixed orthonormal basis (hatu, hatv, hatw) whose
/// hatw axis is the spin axis Omega/|Omega|.
fn uvw(Omega: reb_vec3d, hatu: &mut reb_vec3d, hatv: &mut reb_vec3d, hatw: &mut reb_vec3d) {
    let omega2 = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;
    let omega = omega2.sqrt();

    let mut s = reb_vec3d {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    s.x = Omega.x / omega;
    s.y = Omega.y / omega;
    s.z = Omega.z / omega;

    hatw.x = s.x;
    hatw.y = s.y;
    hatw.z = s.z;

    let fac = (s.x * s.x + s.y * s.y).sqrt();
    if fac != 0.0 {
        hatu.x = -s.y / fac;
        hatu.y = s.x / fac;
        hatu.z = 0.0;
    } else {
        hatu.x = 1.0;
        hatu.y = 0.0;
        hatu.z = 0.0;
    }

    hatv.x = -(hatu.y * hatw.z - hatu.z * hatw.y);
    hatv.y = -(hatu.z * hatw.x - hatu.x * hatw.z);
    hatv.z = -(hatu.x * hatw.y - hatu.y * hatw.x);
}

/// gravitational_harmonics.c `rebx_gravitational_harmonics`.
///
/// C: `void rebx_gravitational_harmonics(struct reb_simulation* const sim,
/// struct rebx_force* const gh, struct reb_particle* const particles, const int N)`.
/// `particles` is `sim.particles` here and `gh` is the index of this
/// force (unused by the effect: it reads no effect parameters).
pub fn rebx_gravitational_harmonics(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    N: usize,
) {
    let G = sim.G;

    for i in 0..N {
        let J2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "J2");
        let J2 = match J2 {
            None => {
                continue;
            }
            Some(J2) => J2,
        };
        if J2 == 0.0 {
            continue;
        }
        let J4 = rebx_get_param_double(rebx, rebx_ap::particle(i), "J4");
        let R_eq = rebx_get_param_double(rebx, rebx_ap::particle(i), "R_eq");
        let R_eq = match R_eq {
            None => {
                continue;
            }
            Some(R_eq) => R_eq,
        };
        let mut Omega = DEFAULTOMEGA;
        let Omegaptr = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        if let Some(Omegaptr) = Omegaptr {
            Omega.x = Omegaptr.x;
            Omega.y = Omegaptr.y;
            Omega.z = Omegaptr.z;
        }
        // The C caches particles[i] here, BEFORE the inner loop writes
        // back into particles[i].a{x,y,z}; pi therefore stays stale.
        let pi: reb_particle = sim.particles[i];

        /* new coordinate basis (body-fixed) */
        let mut hatu = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hatv = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hatw = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        uvw(Omega, &mut hatu, &mut hatv, &mut hatw);

        /* old basis (') in new coordinates */
        let mut hatx_ = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut haty_ = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hatz_ = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        hatx_.x = hatu.x;
        hatx_.y = hatv.x;
        hatx_.z = hatw.x;
        haty_.x = hatu.y;
        haty_.y = hatv.y;
        haty_.z = hatw.y;
        hatz_.x = hatu.z;
        hatz_.y = hatv.z;
        hatz_.z = hatw.z;

        for j in 0..N {
            if j == i {
                continue;
            }
            let pj: reb_particle = sim.particles[j];

            let dx = pj.x - pi.x;
            let dy = pj.y - pi.y;
            let dz = pj.z - pi.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();

            /* new coordinates */
            let du = hatu.x * dx + hatu.y * dy + hatu.z * dz;
            let dv = hatv.x * dx + hatv.y * dy + hatv.z * dz;
            let dw = hatw.x * dx + hatw.y * dy + hatw.z * dz;
            let costheta = dw / r;
            let costheta2 = costheta * costheta;

            let mut au = 0.0;
            let mut av = 0.0;
            let mut aw = 0.0;

            j2_func(
                G, pi.m, Some(J2), R_eq, r, r2, costheta2, du, dv, dw, &mut au, &mut av, &mut aw,
            );
            j4_func(
                G, pi.m, J4, R_eq, r, r2, costheta2, du, dv, dw, &mut au, &mut av, &mut aw,
            );

            /* old coordinates */
            let ax = hatx_.x * au + hatx_.y * av + hatx_.z * aw;
            let ay = haty_.x * au + haty_.y * av + haty_.z * aw;
            let az = hatz_.x * au + hatz_.y * av + hatz_.z * aw;

            sim.particles[j].ax += ax;
            sim.particles[j].ay += ay;
            sim.particles[j].az += az;

            let fac = pj.m / pi.m;

            sim.particles[i].ax -= fac * ax;
            sim.particles[i].ay -= fac * ay;
            sim.particles[i].az -= fac * az;
        }
    }
}

/// gravitational_harmonics.c `j2_potential_func`. See [`j2_func`] for
/// the `Option`/`f64` argument convention.
fn j2_potential_func(
    G: f64,
    mi: f64,
    mj: f64,
    J2: Option<f64>,
    R_eq: f64,
    r: f64,
    r2: f64,
    costheta2: f64,
    H: &mut f64,
) {
    let J2 = match J2 {
        None => {
            return;
        }
        Some(J2) => J2,
    };
    if J2 == 0.0 {
        return;
    }

    let f1 = G * mi * mj * (J2) * (R_eq) * (R_eq) / r2 / r;
    let f2 = 1.0 / 2.0 * (3.0 * costheta2 - 1.0);

    *H += f1 * f2;
}

/// gravitational_harmonics.c `j4_potential_func`. See [`j2_func`] for
/// the `Option`/`f64` argument convention.
fn j4_potential_func(
    G: f64,
    mi: f64,
    mj: f64,
    J4: Option<f64>,
    R_eq: f64,
    r: f64,
    r2: f64,
    costheta2: f64,
    H: &mut f64,
) {
    let J4 = match J4 {
        None => {
            return;
        }
        Some(J4) => J4,
    };
    if J4 == 0.0 {
        return;
    }

    let f1 = G * mi * mj * (J4) * (R_eq) * (R_eq) * (R_eq) * (R_eq) / r2 / r2 / r;
    let f2 = 1.0 / 8.0 * (35.0 * costheta2 * costheta2 - 30.0 * costheta2 + 3.0);

    *H += f1 * f2;
}

/// gravitational_harmonics.c `rebx_gravitational_harmonics_potential`.
///
/// Calculates the potential for all particles with additional gravity
/// field harmonics beyond the monopole (i.e., J2, J4).
///
/// C: `double rebx_gravitational_harmonics_potential(struct rebx_extras* const rebx)`.
/// The C reaches the simulation through `rebx->sim` and bails out with
/// `rebx_error(rebx, "")` returning 0 when that back-pointer is NULL.
/// There is no such back-pointer here, so the simulation is passed
/// explicitly — **sim first, then rebx**, like the `rebx_tools_*`
/// siblings — and the detached branch cannot arise.
pub fn rebx_gravitational_harmonics_potential(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
) -> f64 {
    let G = sim.G;
    let N = sim.N;
    let mut H = 0.0;

    for i in 0..N {
        let J2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "J2");
        let J2 = match J2 {
            None => {
                continue;
            }
            Some(J2) => J2,
        };
        if J2 == 0.0 {
            continue;
        }
        let J4 = rebx_get_param_double(rebx, rebx_ap::particle(i), "J4");
        let R_eq = rebx_get_param_double(rebx, rebx_ap::particle(i), "R_eq");
        let R_eq = match R_eq {
            None => {
                continue;
            }
            Some(R_eq) => R_eq,
        };
        let mut Omega = DEFAULTOMEGA;
        let Omegaptr = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        if let Some(Omegaptr) = Omegaptr {
            Omega.x = Omegaptr.x;
            Omega.y = Omegaptr.y;
            Omega.z = Omegaptr.z;
        }
        let pi: reb_particle = sim.particles[i];

        /* new coordinate basis (body-fixed) */
        let mut hatu = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hatv = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hatw = reb_vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        uvw(Omega, &mut hatu, &mut hatv, &mut hatw);

        for j in 0..N {
            if j == i {
                continue;
            }
            let pj: reb_particle = sim.particles[j];

            let dx = pj.x - pi.x;
            let dy = pj.y - pi.y;
            let dz = pj.z - pi.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();

            /* new coordinates */
            let dw = hatw.x * dx + hatw.y * dy + hatw.z * dz;
            let costheta = dw / r;
            let costheta2 = costheta * costheta;

            j2_potential_func(G, pi.m, pj.m, Some(J2), R_eq, r, r2, costheta2, &mut H);
            j4_potential_func(G, pi.m, pj.m, J4, R_eq, r, r2, costheta2, &mut H);
        }
    }
    H
}
