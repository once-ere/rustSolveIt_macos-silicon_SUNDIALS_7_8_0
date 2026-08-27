//! gas_damping_timescale.rs — translation of REBOUNDx gas_damping_timescale.c
//! Update orbits with prescribed timescales by directly changing orbital
//! elements after each timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Header (gas_damping_timescale.c)
//!
//! ```text
//! @file    gas_damping_timescale.c
//! @brief   Update orbits with prescribed timescales by directly changing orbital elements after each timestep
//! @author  Phoebe Sandhaus <pjs5535@psu.edu>
//! ```
//!
//! ```text
//! $Gas Effects$       // Effect category
//!
//! ======================= ===============================================
//! Authors                 Phoebe Sandhaus
//! Implementation Paper    `Sandhaus et al. 2025 <https://ui.adsabs.harvard.edu/abs/2025ApJ...990...61S/abstract>`_
//! Based on                `Dawson et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...822...54D/abstract>`_, `Kominami & Ida 2002 <https://ui.adsabs.harvard.edu/abs/2002Icar..157...43K/abstract>`_
//! C Example               :ref:`c_example_gas_damping_timescale`
//! Python Example          `GasDampingTimescale.ipynb <https://github.com/PhoebeSandhaus/reboundx_gas_damping/tree/main/ipython_examples/GasDampingTimescale.ipynb>`_
//! ======================= ===============================================
//! ```
//!
//! This updates particles' positions and velocities between timesteps by
//! first calculating a damping timescale for each individual particle, and
//! then applying the timescale to damp both the eccentricity and
//! inclination of the particle. Note: The timescale of damping should be
//! much greater than a particle's orbital period. The damping force should
//! also be small as compared to the gravitational forces on the particle.
//!
//! **Effect Parameters**
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! None                         -           -
//! ============================ =========== ==================================================================
//! ```
//!
//! **Particle Parameters**
//!
//! ```text
//! ============================ =========== ==============================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==============================================================================
//! d_factor (double)            Yes         Depletion factor d in Equation 16 from Dawson et al. 2016; d=1 corresponds
//!                                          roughly to the full minimum mass solar nebula with Sigma_gas (surface gas
//!                                          density) = 1700 g cm^-2 at 1 AU [for d=1]; d>1 corresponds to a more depleted
//!                                          nebula
//! cs_coeff (double)            Yes         Sound speed coefficient; Changing the value will change assumed units;
//!                                          Example: If you are using the units: AU, M_sun, yr -->
//!                                          cs_coeff = 0.272  # AU^(3/4) yr^-1
//! tau_coeff (double)           Yes         Timescale coefficient; Changing the value will change assumed units;
//!                                          Example: If you are using the units: AU, M_sun, yr -->
//!                                          tau_coeff = 0.003  #  yr AU^-2
//! ============================ =========== ==============================================================================
//! ```
//!
//! Note that although the C documentation table above lists `d_factor`,
//! `cs_coeff` and `tau_coeff` together under "Particle Parameters", the C
//! code reads `d_factor` off the *planet* (`planet->ap`) and both
//! `cs_coeff` and `tau_coeff` off the *force* (`force->ap`). The
//! translation reproduces that split exactly.

use crate::core::{rebx_error, rebx_get_param_double, rebx_get_param_int};
use crate::rebxtools::{rebx_com_force, REBX_COORDINATES, REBX_COORDINATES_JACOBI};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_orbit_from_particle, reb_particle, reb_simulation, reb_vec3d};

/// gas_damping_timescale.c `rebx_calculate_gas_damping_timescale`
/// (C: `static`).
///
/// C signature: `(struct reb_simulation* const sim, struct rebx_force*
/// const force, struct reb_particle* planet, struct reb_particle* star)`.
/// Here the force is the index `force_idx`, and `p_index` names the
/// planet's parameter list (C: `planet->ap`).
fn rebx_calculate_gas_damping_timescale(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    p_index: usize,
    planet: reb_particle,
    star: reb_particle,
) -> reb_vec3d {
    let o = reb_orbit_from_particle(sim.G, planet, star);

    let d_factor = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "d_factor");
    let cs_coeff = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "cs_coeff");
    let tau_coeff = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "tau_coeff");

    let mut a = reb_vec3d { x: 0., y: 0., z: 0. }; // C: struct reb_vec3d a = {0};

    let (d_factor, cs_coeff, tau_coeff) = match (d_factor, cs_coeff, tau_coeff) {
        (Some(d_factor), Some(cs_coeff), Some(tau_coeff)) => (d_factor, cs_coeff, tau_coeff),
        _ => {
            rebx_error(
                rebx,
                "Need to set d_factor, cs_coeff, tau_coeff parameters.  See examples in documentation.\n",
            );
            return a;
        }
    };

    // initialize positions and velocities
    let dvx = planet.vx - star.vx;
    let dvy = planet.vy - star.vy;
    let dvz = planet.vz - star.vz;
    let dx = planet.x - star.x;
    let dy = planet.y - star.y;
    let dz = planet.z - star.z;
    let r2 = dx * dx + dy * dy + dz * dz;

    // initial semimajor axis, eccentricity, and inclination
    let a0 = o.a;
    let e0 = o.e;
    let inc0 = o.inc;
    let starMass = star.m;
    let planetMass = planet.m;

    // eccentricity and inclination timescales from Dawson+16 Eqn 16
    let coeff;

    let vk = (sim.G * starMass / a0).sqrt();
    let v = (e0 * e0 + inc0 * inc0).sqrt() * vk;
    let cs = cs_coeff / a0.sqrt().sqrt();
    let v_over_cs = v / cs;

    if v <= cs {
        coeff = 1.;
    } else {
        if inc0 < cs / vk {
            coeff = v_over_cs * v_over_cs * v_over_cs;
        } else {
            coeff = v_over_cs * v_over_cs * v_over_cs * v_over_cs;
        }
    }

    let tau_e = -tau_coeff * (d_factor) * a0 * a0 * (starMass / planetMass) * coeff;
    let tau_inc = 2. * tau_e; // from Kominami & Ida 2002 [Eqs. 2.9 and 2.10]

    if tau_e < f64::INFINITY || tau_inc < f64::INFINITY {
        let vdotr = dx * dvx + dy * dvy + dz * dvz;
        let prefac = 2. * vdotr / r2 / tau_e;
        a.x += prefac * dx;
        a.y += prefac * dy;
        a.z += prefac * dz + 2. * dvz / tau_inc;
    }
    a
}

/// gas_damping_timescale.c `rebx_gas_damping_timescale` — the force
/// entry point registered by `rebx_load_force("gas_damping_timescale")`.
pub fn rebx_gas_damping_timescale(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let ptr = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "coordinates");
    let mut coordinates: REBX_COORDINATES = REBX_COORDINATES_JACOBI; // Default
    if let Some(ptr) = ptr {
        coordinates = ptr;
    }
    let back_reactions_inclusive: i32 = 1;
    let reference_name = "primary";
    rebx_com_force(
        sim,
        rebx,
        force_idx,
        coordinates,
        back_reactions_inclusive,
        reference_name,
        rebx_calculate_gas_damping_timescale,
        N,
    );
}
