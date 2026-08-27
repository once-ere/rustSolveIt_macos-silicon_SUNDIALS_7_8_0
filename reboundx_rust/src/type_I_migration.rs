//! type_I_migration.rs — translation of REBOUNDx type_I_migration.c
//! Applies Type I migration, damping eccentricity, angular momentum and inclination.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Authors: Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>,
//! Gabriele Pichierri <gabrielepichierri@gmail.com>
//!
//! # $Orbit Modifications$
//!
//! ======================= ===============================================
//! Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
//! Implementation Paper    `Kajtazi et al 2022 <https://ui.adsabs.harvard.edu/abs/2022arXiv221106181K/abstract>`_.
//! Based on                `Cresswell & Nelson 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..677C/abstract>`_, and `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>`_.
//! C example               :ref:`c_example_type_I_migration`
//! Python example          `TypeIMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TypeIMigration.ipynb>`_.
//! ======================= ===============================================
//!
//! This applies Type I migration, damping eccentricity, angular momentum and inclination.
//! The base of the code is the same as the modified orbital forces one written by D. Tamayo, H. Rein.
//! It also allows for parameters describing an inner disc edge, modeled using the implementation in inner_disk_edge.c.
//! Note that this code is not machine independent since power laws were not possible to avoid all together.
//!
//! # Effect Parameters
//!
//! ===================================== =========== ==================================================================================================================
//! Field (C type)                        Required    Description
//! ===================================== =========== ==================================================================================================================
//! ide_position (double)                 No          The position of the inner disk edge in code units
//! ide_width (double)                    No          The disk edge width (planet will stop within ide_width of ide_position)
//! tIm_surface_density_1 (double)        Yes         Disk surface density at one code unit from the star; used to find the surface density at any distance from the star
//! tIm_scale_height_1 (double)           Yes         The scale height at one code unit from the star; used to find the aspect ratio at any distance from the star
//! tIm_surface_density_exponent (double) Yes         Exponent of disk surface density, indicative of the surface density profile of the disk
//! tIm_flaring_index (double)            Yes         The flaring index; 1 means disk is irradiated by only the stellar flux
//! ===================================== =========== ==================================================================================================================

use crate::core::{rebx_get_param_double, rebx_get_param_int};
use crate::inner_disk_edge::rebx_calculate_planet_trap;
use crate::rebxtools::{rebx_com_force, REBX_COORDINATES, REBX_COORDINATES_JACOBI};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_orbit_from_particle_err_t, reb_particle, reb_simulation, reb_vec3d};

/* Calculating the t_wave: damping timescale or orbital evolution timescale, from Tanaka & Ward 2004.
h = aspect ratio, h2 = aspect ratio squared, sma = semi-major axis, sd = disk surface denisty to be calculated at every r, ms = stellar mass, mp = planet mass */
/// type_I_migration.c `rebx_calculate_damping_timescale`.
pub fn rebx_calculate_damping_timescale(
    G: f64,
    sd0: f64,
    r: f64,
    s: f64,
    ms: f64,
    mp: f64,
    sma: f64,
    h2: f64,
) -> f64 {
    let sd: f64;
    let t_wave: f64;

    sd = sd0 * r.powf(-s);
    t_wave = ((ms * ms * ms).sqrt() * h2 * h2) / (mp * sd * (sma * G).sqrt());

    t_wave
}

/* Calculating the eccentricity damping timescale t_e = -e/(de/dt), from Cresswell & Nelson 2008.
eh=e/h, ih = i/h, wave is a funcion variable name which will be the t_wave function*/
/// type_I_migration.c `rebx_calculate_eccentricity_damping_timescale`.
pub fn rebx_calculate_eccentricity_damping_timescale(wave: f64, eh: f64, ih: f64) -> f64 {
    let t_e: f64;
    t_e = (wave / 0.780) * (1. - (0.14 * eh * eh) + (0.06 * eh * eh * eh) + (0.18 * eh * ih * ih));

    t_e
}

/* Calculating the migration timescale t_mig = - angmom/torque, from Cresswell & Nelson 2008*/
/// type_I_migration.c `rebx_calculate_migration_timescale`.
pub fn rebx_calculate_migration_timescale(wave: f64, eh: f64, ih: f64, h2: f64, s: f64) -> f64 {
    let Pe: f64;
    let t_mig: f64;
    let term: f64;
    let term2: f64;
    let term3: f64;
    term = eh / 2.25;
    term2 = (eh / 2.84) * (eh / 2.84);
    term3 = (eh / 2.02) * (eh / 2.02);
    Pe = (1. + term.powf(1.2) + term2 * term2 * term2) / (1. - term3 * term3);
    t_mig = ((2. * wave) / (2.7 + 1.1 * s))
        * (1. / h2)
        * (Pe + (Pe / Pe.abs()) * ((0.070 * ih) + (0.085 * ih * ih * ih * ih) - (0.080 * eh * ih * ih)));

    t_mig
}

/* Calculating the inclination damping timescale t_i = -i/(di/dt), from Cresswell & Nelson 2008*/
/// type_I_migration.c `rebx_calculate_inclination_damping_timescale`.
pub fn rebx_calculate_inclination_damping_timescale(wave: f64, eh: f64, ih: f64) -> f64 {
    let t_i: f64;
    t_i = (wave / 0.544) * (1. - (0.30 * ih * ih) + (0.24 * ih * ih * ih) + (0.14 * eh * eh * ih));

    t_i
}

/// type_I_migration.c `rebx_calculate_modify_orbits_with_type_I_migration`
/// (C: `static`). Matches [`crate::rebxtools::rebx_calculate_force_fn`];
/// the C's `struct reb_particle* p` and `* source` are passed by value
/// here and the C's `force` pointer as `force_idx`.
fn rebx_calculate_modify_orbits_with_type_I_migration(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    _p_index: usize,
    p: reb_particle,
    source: reb_particle,
) -> reb_vec3d {
    let invtau_mig: f64;
    let tau_e: f64;
    let tau_inc: f64;

    /* Default values for the parameters in case the user forgets to define them when using this code */
    let mut beta = 0.0;
    let mut h0 = 0.01;
    let mut sd0 = 0.0;
    let mut s = 0.0;
    let mut dedge = 0.0;
    let mut hedge = 0.0;

    /* Parameters that should be changed/set in Python notebook or in C outside of this */
    let dedge_ptr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ide_position");
    let hedge_ptr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ide_width");
    let beta_ptr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "tIm_flaring_index");
    let s_ptr = rebx_get_param_double(
        rebx,
        rebx_ap::force(force_idx),
        "tIm_surface_density_exponent",
    );
    let sd0_ptr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "tIm_surface_density_1");
    let h0_ptr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "tIm_scale_height_1");

    /* Accessing the calculated semi-major axis, eccentricity and inclination for each integration step, via modify_orbits_direct where they are calculated and returned*/
    let mut err: i32 = 0;
    let o = reb_orbit_from_particle_err_t(sim.G, p, source, sim.t, &mut err);

    let a0 = o.a;
    let e0 = o.e;
    let inc0 = o.inc;
    let mp = p.m;
    let ms = source.m;

    let dvx = p.vx - source.vx;
    let dvy = p.vy - source.vy;
    let dvz = p.vz - source.vz;
    let dx = p.x - source.x;
    let dy = p.y - source.y;
    let dz = p.z - source.z;
    let r2 = dx * dx + dy * dy + dz * dz;

    if let Some(v) = beta_ptr {
        beta = v;
    }
    if let Some(v) = s_ptr {
        s = v;
    }
    if let Some(v) = sd0_ptr {
        sd0 = v;
    }
    if let Some(v) = h0_ptr {
        h0 = v;
    }
    if let Some(v) = dedge_ptr {
        dedge = v;
    }
    if let Some(v) = hedge_ptr {
        hedge = v;
    }

    /* Calculating the aspect ratio evaluated at the position of the planet, r and defining other variables */

    let h = (h0) * r2.powf(beta / 2.);
    let h2 = h * h;

    let eh = e0 / h;
    let ih = inc0 / h;

    let G = sim.G;
    let wave = rebx_calculate_damping_timescale(G, sd0, r2.sqrt(), s, ms, mp, a0, h2);
    invtau_mig = rebx_calculate_planet_trap(a0, dedge, hedge)
        / (rebx_calculate_migration_timescale(wave, eh, ih, h2, s));
    tau_e = rebx_calculate_eccentricity_damping_timescale(wave, eh, ih);
    tau_inc = rebx_calculate_inclination_damping_timescale(wave, eh, ih);

    let mut a = reb_vec3d::default(); // C: struct reb_vec3d a = {0};

    if invtau_mig != 0.0 {
        a.x = -dvx * (invtau_mig);
        a.y = -dvy * (invtau_mig);
        a.z = -dvz * (invtau_mig);
    }

    if tau_e < f64::INFINITY || tau_inc < f64::INFINITY {
        let vdotr = dx * dvx + dy * dvy + dz * dvz;
        let prefac = -2. * vdotr / r2 / tau_e;
        a.x += prefac * dx;
        a.y += prefac * dy;
        a.z += prefac * dz - 2. * dvz / tau_inc;
    }
    a
}

/// type_I_migration.c `rebx_modify_orbits_with_type_I_migration` — the
/// force entry point registered under the name `"type_I_migration"`.
pub fn rebx_modify_orbits_with_type_I_migration(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let ptr = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "coordinates");
    let mut coordinates: REBX_COORDINATES = REBX_COORDINATES_JACOBI; // Default
    if let Some(v) = ptr {
        coordinates = v;
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
        rebx_calculate_modify_orbits_with_type_I_migration,
        N,
    );
}
