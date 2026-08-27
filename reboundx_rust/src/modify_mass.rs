//! modify_mass.rs — translation of REBOUNDx modify_mass.c
//! Add exponential mass loss/growth between timesteps in the simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Mass Modifications
//!
//! ```text
//! ======================= ===============================================
//! Authors                 D. Tamayo
//! Implementation Paper    Kostov et al., 2016
//!                         https://ui.adsabs.harvard.edu/abs/2016ApJ...832..183K/abstract
//! Based on                None
//! C Example               c_example_modify_mass
//! Python Example          ModifyMass.ipynb
//! ======================= ===============================================
//! ```
//!
//! This adds exponential mass growth/loss to individual particles every
//! timestep. Set particles' `tau_mass` parameter to a negative value for
//! mass loss, positive for mass growth.
//!
//! **Effect Parameters**
//!
//! *None*
//!
//! **Particle Parameters**
//!
//! Only particles with their `tau_mass` parameter set will have their
//! masses affected.
//!
//! ```text
//! ============================ =========== =======================================================
//! Name (C type)                Required    Description
//! ============================ =========== =======================================================
//! tau_mass (double)            Yes         e-folding mass loss (<0) or growth (>0) timescale
//! ============================ =========== =======================================================
//! ```

use rebound_rs::{reb_simulation, reb_simulation_move_to_com};

use crate::core::rebx_get_param_double;
use crate::types::{rebx_ap, rebx_extras};

/// modify_mass.c `rebx_modify_mass`.
///
/// C signature: `void rebx_modify_mass(struct reb_simulation* const sim,
/// struct rebx_operator* const operator, const double dt)`. The operator
/// carries no parameters of its own, so `_operator_idx` is unused, exactly
/// as `operator` is unused in the C.
pub fn rebx_modify_mass(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _operator_idx: usize,
    dt: f64,
) {
    let N = sim.N;
    for i in 0..N {
        // C: struct reb_particle* const p = &sim->particles[i];
        let tau_mass = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau_mass");
        if let Some(tau_mass) = tau_mass {
            let p = &mut sim.particles[i];
            p.m += p.m * dt / tau_mass;
        }
    }
    reb_simulation_move_to_com(sim);
}
