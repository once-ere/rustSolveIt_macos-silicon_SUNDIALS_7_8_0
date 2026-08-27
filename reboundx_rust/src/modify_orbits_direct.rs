//! modify_orbits_direct.rs — translation of REBOUNDx modify_orbits_direct.c
//! Update orbital elements with prescribed timescales by directly changing orbital elements after each timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! @author  Dan Tamayo <tamayo.daniel@gmail.com>
//!
//! # $Orbit Modifications$
//!
//! ======================= ===============================================
//! Authors                 D. Tamayo
//! Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
//! Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_.
//! C Example               :ref:`c_example_modify_orbits`
//! Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_,
//!                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
//! ======================= ===============================================
//!
//! This updates particles' positions and velocities between timesteps to achieve the desired changes to the
//! osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for
//! Omega/omega.
//! This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting
//! dynamics.
//! One can also adjust the coupling parameter `p` between eccentricity and semimajor axis evolution, as well as
//! whether the damping is done on Jacobi, barycentric or heliocentric elements.
//! Since this method changes osculating (i.e., two-body) elements, it can give unphysical results in highly
//! perturbed systems.
//!
//! **Effect Parameters**
//!
//! If p is not set, it defaults to 0.  If coordinates not set, defaults to using Jacobi coordinates.
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! p (double)                   No          Coupling parameter between eccentricity and semimajor axis evolution
//!                                          (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
//!                                          eccentricity evolution at constant angular momentum.
//! coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
//!                                          See the examples for usage.
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! One can pick and choose which particles have which parameters set.
//! For each particle, any unset parameter is ignored.
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! tau_a (double)               No          Semimajor axis exponential growth/damping timescale
//! tau_e (double)               No          Eccentricity exponential growth/damping timescale
//! tau_inc (double)             No          Inclination axis exponential growth/damping timescale
//! tau_Omega (double)           No          Period of linear nodal precession/regression
//! tau_omega (double)           No          Period of linear apsidal precession/regression
//! ============================ =========== ==================================================================

use std::f64::consts::PI;

use crate::core::{rebx_get_param_double, rebx_get_param_int};
use crate::inner_disk_edge::rebx_calculate_planet_trap;
use crate::rebxtools::{rebx_tools_com_ptm, REBX_COORDINATES, REBX_COORDINATES_JACOBI};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{
    reb_orbit_from_particle_err, reb_particle, reb_particle_from_orbit, reb_simulation,
};

/// modify_orbits_direct.c `rebx_calculate_modify_orbits_direct` (C: `static`).
///
/// Per-particle callback handed to [`rebx_tools_com_ptm`]. The C reads the
/// particle's parameter list through `p->ap`; here `p_index` names that list
/// (`rebx_ap::particle(p_index)`), and the operator's own parameters are read
/// through `rebx_ap::operator_(operator_idx)` (C: `operator->ap`).
fn rebx_calculate_modify_orbits_direct(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    operator_idx: usize,
    p_index: usize,
    p: reb_particle,
    primary: reb_particle,
    dt: f64,
) -> reb_particle {
    let mut err: i32 = 0;
    let mut o = reb_orbit_from_particle_err(sim.G, p, primary, &mut err);
    if err != 0 {
        // mass of primary was 0 or p = primary.  Return same particle without doing anything.
        return p;
    }

    let tau_a_ptr = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_a");
    let tau_e = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_e");
    let tau_inc = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_inc");
    let tau_omega = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_omega");
    let tau_Omega = rebx_get_param_double(rebx, rebx_ap::particle(p_index), "tau_Omega");

    //Implement the planet trap
    // C: `double invtau_a = 0.0;` here. The 0.0 is never read — invtau_a is
    // overwritten before its only use and is local to the tau_a branch below —
    // so the declaration moves inside that branch to keep the build
    // warning-free. No arithmetic changes.
    let dedge = rebx_get_param_double(rebx, rebx_ap::operator_(operator_idx), "ide_position");
    let hedge = rebx_get_param_double(rebx, rebx_ap::operator_(operator_idx), "ide_width");

    let a0 = o.a;
    let e0 = o.e;
    let inc0 = o.inc;

    if let Some(tau_a_ptr) = tau_a_ptr {
        let mut invtau_a = 1.0 / tau_a_ptr;
        if let (Some(dedge), Some(hedge)) = (dedge, hedge) {
            invtau_a *= rebx_calculate_planet_trap(a0, dedge, hedge);
        }
        o.a += a0 * dt * invtau_a;
    }
    if let Some(tau_e) = tau_e {
        o.e += e0 * dt / tau_e;
    }
    if let Some(tau_inc) = tau_inc {
        o.inc += inc0 * dt / tau_inc;
    }
    if let Some(tau_omega) = tau_omega {
        o.omega += 2. * PI * dt / tau_omega;
    }
    if let Some(tau_Omega) = tau_Omega {
        o.Omega += 2. * PI * dt / tau_Omega;
    }

    if let Some(tau_e) = tau_e {
        let p_param = rebx_get_param_double(rebx, rebx_ap::operator_(operator_idx), "p");
        if let Some(p_param) = p_param {
            o.a += 2. * a0 * e0 * e0 * p_param * dt / tau_e; // Coupling term between e and a
        }
    }
    reb_particle_from_orbit(
        sim.G, primary, p.m, o.a, o.e, o.inc, o.Omega, o.omega, o.f,
    )
}

/// modify_orbits_direct.c `rebx_modify_orbits_direct` — the operator step
/// function registered under the name `"modify_orbits_direct"`.
pub fn rebx_modify_orbits_direct(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    operator_idx: usize,
    dt: f64,
) {
    let ptr = rebx_get_param_int(rebx, rebx_ap::operator_(operator_idx), "coordinates");
    let mut coordinates: REBX_COORDINATES = REBX_COORDINATES_JACOBI;
    if let Some(ptr) = ptr {
        coordinates = ptr;
    }
    let back_reactions_inclusive = 1;
    let reference_name = "primary";
    rebx_tools_com_ptm(
        sim,
        rebx,
        operator_idx,
        coordinates,
        back_reactions_inclusive,
        reference_name,
        rebx_calculate_modify_orbits_direct,
        dt,
    );
}
