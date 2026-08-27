//! tides_constant_time_lag.rs — translation of REBOUNDx tides_constant_time_lag.c
//! Adds constant time lag tides raised on the primary, on the orbiting bodies, or both.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Authors: Stanley A. Baronett <stanley.a.baronett@gmail.com>,
//! Dan Tamayo <tamayo.daniel@gmail.com>, Noah Ferich.
//!
//! # $Tides$
//!
//! ======================= ===============================================
//! Authors                 Stanley A. Baronett, D. Tamayo, Noah Ferich
//! Implementation Paper    `Baronett et al., 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.6001B/abstract>`_.
//! Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_, `Bolmont et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...583A.116B/abstract>`_.
//! C Example               :ref:`c_example_tides_constant_time_lag`.
//! Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
//! ======================= ===============================================
//!
//! This adds constant time lag tidal interactions between orbiting bodies in the
//! simulation and the primary, both from tides raised on the primary and on the
//! other bodies. In all cases, we need to set masses for all the particles that
//! will feel these tidal forces. After that, we can choose to include tides
//! raised on the primary, on the "planets", or both, by setting the respective
//! bodies' physical radius `particles[i].r`, `k2` (potential Love number of
//! degree 2), constant time lag `tau`, and rotation rate `Omega`. See Baronett
//! et al. (2021), Hut (1981), and Bolmont et al. 2015 above.
//!
//! If tau is not set, it will default to zero and yield the conservative piece
//! of the tidal potential.
//!
//! **Effect Parameters**
//!
//! None
//!
//! **Particle Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
//! tctl_k2 (float)              Yes         Potential Love number of degree 2.
//! tctl_tau (float)             No          Constant time lag. If not set will default to 0 and give conservative tidal potential.
//! OmegaMag (float)             No          Angular rotation frequency. If not set will default to 0.
//! ============================ =========== ==================================================================

use rebound_rs::{reb_particle, reb_simulation};

use crate::core::rebx_get_param_double;
use crate::types::{rebx_ap, rebx_extras};

/// C: `static void rebx_calculate_tides(struct reb_particle* source,
/// struct reb_particle* target, ...)`.
///
/// The two C pointers become indices into the one particle array, so
/// that both bodies can be written through safe Rust. `source` and
/// `target` are never the same particle at any call site (one of them is
/// always `particles[0]` and the other has index `i >= 1`).
fn rebx_calculate_tides(
    particles: &mut [reb_particle],
    source: usize,
    target: usize,
    G: f64,
    k2: f64,
    tau: f64,
    Omega: f64,
) {
    let ms = particles[source].m;
    let mt = particles[target].m;
    let Rt = particles[target].r;

    let mratio = ms / mt; // have already checked for 0 and inf
    let fac = mratio * k2 * Rt * Rt * Rt * Rt * Rt;

    let dx = particles[target].x - particles[source].x;
    let dy = particles[target].y - particles[source].y;
    let dz = particles[target].z - particles[source].z;
    let dr2 = dx * dx + dy * dy + dz * dz;
    let prefac = -3. * G / (dr2 * dr2 * dr2 * dr2) * fac;
    let mut rfac = prefac;

    if tau != 0. {
        let dvx = particles[target].vx - particles[source].vx;
        let dvy = particles[target].vy - particles[source].vy;
        let dvz = particles[target].vz - particles[source].vz;

        rfac *= 1. + 3. * tau / dr2 * (dx * dvx + dy * dvy + dz * dvz);
        let thetafac = -prefac * tau;

        let hx = dy * dvz - dz * dvy;
        let hy = dz * dvx - dx * dvz;
        let hz = dx * dvy - dy * dvx;

        let thetadotcrossrx = (hy * dz - hz * dy) / dr2; // vec(thetadot) cross vec(r) Bolmont Eq. 7
        let thetadotcrossry = (hz * dx - hx * dz) / dr2;
        let thetadotcrossrz = (hx * dy - hy * dx) / dr2;

        // Assumes all spins vec(Omega) = Omega zhat, i.e., spin fixed along z axis
        let Omegacrossrx = -Omega * dy;
        let Omegacrossry = Omega * dx;
        let Omegacrossrz = 0.;

        particles[target].ax += thetafac * ms * (Omegacrossrx - thetadotcrossrx);
        particles[target].ay += thetafac * ms * (Omegacrossry - thetadotcrossry);
        particles[target].az += thetafac * ms * (Omegacrossrz - thetadotcrossrz);
        particles[source].ax -= thetafac * mt * (Omegacrossrx - thetadotcrossrx);
        particles[source].ay -= thetafac * mt * (Omegacrossry - thetadotcrossry);
        particles[source].az -= thetafac * mt * (Omegacrossrz - thetadotcrossrz);
    }

    particles[target].ax += rfac * ms * dx;
    particles[target].ay += rfac * ms * dy;
    particles[target].az += rfac * ms * dz;
    particles[source].ax -= rfac * mt * dx;
    particles[source].ay -= rfac * mt * dy;
    particles[source].az -= rfac * mt * dz;
}

/// C: `void rebx_tides_constant_time_lag(struct reb_simulation* const sim,
/// struct rebx_force* const tides, struct reb_particle* const particles,
/// const int N)`.
///
/// The C `tides` force pointer is unused in the body (this effect has no
/// effect parameters), so `_force_idx` is likewise unused here.
pub fn rebx_tides_constant_time_lag(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    N: usize,
) {
    let G = sim.G;

    // Calculate tides raised on star
    // assumes nearly Keplerian motion around a single primary (particles[0])
    if sim.particles[0].m == 0. {
        // nothing makes sense if primary has no mass
        return;
    }
    let k2 = rebx_get_param_double(rebx, rebx_ap::particle(0), "tctl_k2");
    // tides on star only nonzero if k2 and finite size are set
    if let Some(k2) = k2 {
        if sim.particles[0].r != 0. {
            // We don't require time lag tau to be set. Might just want conservative piece of tidal potential
            let mut tau = 0.;
            let mut Omega = 0.;
            let tauptr = rebx_get_param_double(rebx, rebx_ap::particle(0), "tctl_tau");
            if let Some(tauptr) = tauptr {
                tau = tauptr;
                let Omegaptr = rebx_get_param_double(rebx, rebx_ap::particle(0), "OmegaMag");
                if let Some(Omegaptr) = Omegaptr {
                    Omega = Omegaptr;
                }
            }
            for i in 1..N {
                // particles[i] is the planet raising the tides on the star
                if sim.particles[i].m == 0. {
                    continue;
                }
                rebx_calculate_tides(&mut sim.particles, i, 0, G, k2, tau, Omega);
            }
        }
    }

    // Calculate tides raised on the planets
    // Source is always the star, particles[0] (no planet-planet tides)
    for i in 1..N {
        let k2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "tctl_k2");
        let k2 = match k2 {
            None => continue,
            Some(k2) => k2,
        };
        if sim.particles[i].r == 0. || sim.particles[i].m == 0. {
            continue;
        }
        let mut tau = 0.;
        let mut Omega = 0.;
        let tauptr = rebx_get_param_double(rebx, rebx_ap::particle(i), "tctl_tau");
        if let Some(tauptr) = tauptr {
            tau = tauptr;
            let Omegaptr = rebx_get_param_double(rebx, rebx_ap::particle(i), "OmegaMag");
            if let Some(Omegaptr) = Omegaptr {
                Omega = Omegaptr;
            }
        }
        rebx_calculate_tides(&mut sim.particles, 0, i, G, k2, tau, Omega);
    }
}

/// Calculate potential of conservative piece of tidal interaction.
///
/// C: `static double rebx_calculate_tides_potential(struct reb_particle* source,
/// struct reb_particle* target, const double G, const double k2)`.
fn rebx_calculate_tides_potential(
    particles: &[reb_particle],
    source: usize,
    target: usize,
    G: f64,
    k2: f64,
) -> f64 {
    let ms = particles[source].m;
    let mt = particles[target].m;
    let Rt = particles[target].r;

    let mratio = ms / mt; // have already checked for 0 and inf
    let fac = mratio * k2 * Rt * Rt * Rt * Rt * Rt;

    let dx = particles[target].x - particles[source].x;
    let dy = particles[target].y - particles[source].y;
    let dz = particles[target].z - particles[source].z;
    let dr2 = dx * dx + dy * dy + dz * dz;

    -1. / 2. * G * ms * mt / (dr2 * dr2 * dr2) * fac
}

/// C: `double rebx_tides_constant_time_lag_potential(struct rebx_extras* const rebx)`.
///
/// The C reaches the simulation through the `rebx->sim` back-pointer and
/// bails out with `rebx_error` when it is NULL. `rebx_extras` carries no
/// back-pointer here, so the simulation is passed explicitly and that
/// NULL branch has no counterpart.
pub fn rebx_tides_constant_time_lag_potential(sim: &reb_simulation, rebx: &rebx_extras) -> f64 {
    let N = sim.N;
    let particles = &sim.particles;
    let G = sim.G;
    let mut H = 0.;

    // Calculate tides raised on star
    // assumes nearly Keplerian motion around a single primary (particles[0])
    if particles[0].m == 0. {
        // No potential with massless primary
        return 0.;
    }
    let k2 = rebx_get_param_double(rebx, rebx_ap::particle(0), "tctl_k2");
    // tides on star only nonzero if k2 and finite size are set
    if let Some(k2) = k2 {
        if particles[0].r != 0. {
            for i in 1..N {
                // particles[i] is the planet raising the tides on the star
                if particles[i].m == 0. {
                    continue;
                }
                H += rebx_calculate_tides_potential(particles, i, 0, G, k2);
            }
        }
    }

    // Calculate tides raised on the planets
    // Source is always the star, particles[0] (no planet-planet tides)
    for i in 1..N {
        let k2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "tctl_k2");
        let k2 = match k2 {
            None => continue,
            Some(k2) => k2,
        };
        if particles[i].r == 0. || particles[i].m == 0. {
            continue;
        }
        H += rebx_calculate_tides_potential(particles, 0, i, G, k2);
    }

    H
}
