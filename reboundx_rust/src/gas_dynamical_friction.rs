//! gas_dynamical_friction.rs — translation of REBOUNDx gas_dynamical_friction.c
//! Gas drag from a thin disk with a power-law density profile.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! `@author` Aleksey Generozov
//!
//! # Gas Effects
//!
//! ```text
//! ======================= ===============================================
//! Authors                 A. Generozov, H. Perets
//! Implementation Paper    Generozov and Perets 2022 <https://arxiv.org/abs/2212.11301>
//! Based on                Ostriker 1999 (with simplifications)
//!                         <https://ui.adsabs.harvard.edu/abs/1999ApJ...513..252O/abstract>,
//!                         Just et al 2012
//!                         <https://ui.adsabs.harvard.edu/abs/2012ApJ...758...51J/abstract>
//! C Example               c_example_gas_dynamical_friction
//! Python Example          GasDynamicalFriction.ipynb
//! ======================= ===============================================
//! ```
//!
//! **Effect Parameters**
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! rhog (double)                Yes         Normalization of density. Density in the disk midplane is rhog*r^alpha_rhog
//! alpha_rhog (double)          Yes         Power-law slope of the power-law density profile.
//! cs (double)                  Yes         Normalization of the sound speed. Sound speed has profile cs*r^alpha_cs
//! alpha_cs (double)            Yes         Power-law slope of the sound speed
//! xmin (double)                Yes         Dimensionless parameter that determines the Coulomb logarithm (ln(L) =log (1/xmin))
//! hr (double)                  Yes         Aspect ratio of the disk
//! Qd (double)                  Yes         Prefactor for geometric drag
//! ============================ =========== ==================================================================
//! ```
//!
//! The parameters are looked up on the force under the names
//! `gas_df_rhog`, `gas_df_alpha_rhog`, `gas_df_cs`, `gas_df_alpha_cs`,
//! `gas_df_xmin`, `gas_df_hr` and `gas_df_Qd`.
//!
//! **Particle Parameters**
//!
//! None.

use std::f64::consts::PI;

use rebound_rs::{reb_particle, reb_particle_isub, reb_simulation, reb_simulation_error};

use crate::core::rebx_get_param_double;
use crate::types::{rebx_ap, rebx_extras};

/// gas_dynamical_friction.c `mach_piece_sub`.
fn mach_piece_sub(mach: f64) -> f64 {
    //Using powerlaw expansion at low mach numbers for numerical reasons
    if mach < 0.02 {
        return mach * mach * mach / 3. + mach * mach * mach * mach * mach / 5.;
    }
    //Subsonic expression from Ostriker...
    0.5 * ((1.0 + mach) / (1.0 - mach)).ln() - mach
}

/// gas_dynamical_friction.c `calculate_pre_factor`.
///
/// `vel` and `t` are part of the C signature but unused in its body;
/// they are kept so the correspondence to the C stays visible.
#[allow(unused_variables)]
fn calculate_pre_factor(mach: f64, vel: f64, t: f64, xmin: f64) -> f64 {
    //Simplified version of the Ostriker dynamical friction formula...
    let coul = (1.0 / xmin).ln();
    if mach >= 1.0 {
        coul
    } else {
        coul.min(mach_piece_sub(mach))
    }
}

/// gas_dynamical_friction.c `get_vrel_disk`.
fn get_vrel_disk(p: reb_particle, GMBH: f64, vrel: &mut [f64; 3], hr: f64) {
    let rcyl = (p.x * p.x + p.y * p.y).sqrt();
    let sin_phi = p.y / rcyl;
    let cos_phi = p.x / rcyl;

    let vk = (GMBH / rcyl).sqrt() * (1.0 - hr * hr);
    vrel[0] = p.vx + vk * sin_phi;
    vrel[1] = p.vy - vk * cos_phi;
    vrel[2] = p.vz;
}

/// gas_dynamical_friction.c `rebx_calculate_gas_dynamical_friction`.
///
/// The C takes `particles` separately; here it is `sim.particles`, which
/// is the same array at every call site of this effect. `N` is part of
/// the C signature but unused in its body (the loop bound is `sim->N`).
#[allow(unused_variables)]
fn rebx_calculate_gas_dynamical_friction(
    sim: &mut reb_simulation,
    N: usize,
    rhog: f64,
    alpha_rhog: f64,
    cs: f64,
    alpha_cs: f64,
    xmin: f64,
    hr: f64,
    Qd: f64,
) {
    let N_real = sim.N;
    let bh = sim.particles[0];
    let G = sim.G;
    for i in 1..N_real {
        let p = sim.particles[i];
        let mut diff = p;
        let bh2 = bh;
        reb_particle_isub(&mut diff, &bh2);
        let rcyl = (diff.x * diff.x + diff.y * diff.y).sqrt();
        let mut vrel: [f64; 3] = [0.; 3];
        get_vrel_disk(diff, G * bh.m, &mut vrel, hr);
        let vrel_norm = (vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]).sqrt();
        let mach = vrel_norm / (cs * rcyl.powf(alpha_cs));
        let t = sim.t;

        let integ = calculate_pre_factor(mach, vrel_norm, t, xmin);
        //Accounting for vertical dependence of the density with a Gaussian function
        //scale height is defined by user-defined aspect ratio. Truncate the disc vertically
        //at 10 scale heights...
        let h = hr * rcyl;
        let vert = if diff.z.abs() < (10. * h) {
            (-diff.z * diff.z / (2.0 * h * h)).exp()
        } else {
            0.
        };
        let rhog_loc = rhog * rcyl.powf(alpha_rhog) * vert;
        let mp = p.m;
        let rstar = p.r;
        let fc = 4. * PI * (G * G) * mp * (rhog_loc) / (vrel_norm * vrel_norm * vrel_norm) * integ
            + PI * rhog_loc * rstar * rstar * vrel_norm * Qd / mp;

        sim.particles[i].ax -= fc * vrel[0];
        sim.particles[i].ay -= fc * vrel[1];
        sim.particles[i].az -= fc * vrel[2];
    }
}

/// gas_dynamical_friction.c `rebx_gas_dynamical_friction`.
///
/// Deviation from the C: after reporting a missing required parameter
/// the C goes on to dereference the NULL pointer it just tested. Here
/// every missing parameter is reported in the same order and with the
/// same message, and the function then returns without applying the
/// force.
pub fn rebx_gas_dynamical_friction(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let rhog = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_rhog");
    if rhog.is_none() {
        reb_simulation_error(sim, "Need to specify a gas density\n");
    }
    let alpha_rhog = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_alpha_rhog");
    if alpha_rhog.is_none() {
        reb_simulation_error(sim, "Need to specify a profile for gas density\n");
    }
    let cs = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_cs");
    if cs.is_none() {
        reb_simulation_error(sim, "Need to set a sound speed.\n");
    }
    let alpha_cs = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_alpha_cs");
    if alpha_cs.is_none() {
        reb_simulation_error(sim, "Need to specify a profile for the sound speed\n");
    }
    let xmin = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_xmin");
    if xmin.is_none() {
        reb_simulation_error(sim, "Need to set a cutoff.\n");
    }
    let hr = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_hr");
    if hr.is_none() {
        reb_simulation_error(sim, "Need an aspect ratio.\n");
    }
    let Qd = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "gas_df_Qd");
    if Qd.is_none() {
        reb_simulation_error(sim, "Need to specify Qd");
    }

    let (rhog, alpha_rhog, cs, alpha_cs, xmin, hr, Qd) =
        match (rhog, alpha_rhog, cs, alpha_cs, xmin, hr, Qd) {
            (
                Some(rhog),
                Some(alpha_rhog),
                Some(cs),
                Some(alpha_cs),
                Some(xmin),
                Some(hr),
                Some(Qd),
            ) => (rhog, alpha_rhog, cs, alpha_cs, xmin, hr, Qd),
            _ => return,
        };

    rebx_calculate_gas_dynamical_friction(
        sim, N, rhog, alpha_rhog, cs, alpha_cs, xmin, hr, Qd,
    );
}
