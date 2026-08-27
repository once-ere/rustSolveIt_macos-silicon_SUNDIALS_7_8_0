//! track_min_distance.rs — translation of REBOUNDx track_min_distance.c
//! Track minimum distance of secondaries from primary each timestep and log results.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Misc
//!
//! ```text
//! ======================= ===============================================
//! Authors                 D. Tamayo
//! Implementation Paper    Tamayo, Rein, Shi and Hernandez, 2019
//!                         https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract
//! Based on                None
//! C Example               c_example_track_min_distance
//! Python Example          TrackMinDistance.ipynb
//! ======================= ===============================================
//! ```
//!
//! For a given particle, this keeps track of that particle's minimum distance
//! from another body in the simulation. User should add parameters to the
//! particular particle whose distance should be tracked.
//!
//! **Effect Parameters**
//!
//! *None*
//!
//! **Particle Parameters**
//!
//! Only particles with their `min_distance` parameter set initially will track
//! their minimum distance. The effect will update this parameter when the
//! particle gets closer than the value of `min_distance`, so the user has to
//! set it initially. By default, distance is measured from `sim.particles[0]`,
//! but you can specify a different particle by setting the
//! `min_distance_from` parameter to the name of the target particle.
//!
//! ```text
//! ================================ =========== =======================================================
//! Name (C type)                    Required    Description
//! ================================ =========== =======================================================
//! min_distance (double)            Yes         Particle's mininimum distance.
//! min_distance_from (char*)        No          Name for particle from which to measure distance
//! min_distance_orbit (reb_orbit)   No          Parameter to store orbital elements at moment corresponding to min_distance (heliocentric)
//! ================================ =========== =======================================================
//! ```

use rebound_rs::{
    reb_orbit_from_particle, reb_simulation, reb_simulation_get_particle_by_name,
    reb_simulation_warning,
};

use crate::core::{
    rebx_get_param_double, rebx_get_param_orbit, rebx_get_param_string, rebx_set_param_double,
    rebx_set_param_orbit,
};
use crate::types::{rebx_ap, rebx_extras};

/// track_min_distance.c `rebx_track_min_distance`.
///
/// C signature: `void rebx_track_min_distance(struct reb_simulation* const sim,
/// struct rebx_operator* const operator, const double dt)`. The operator
/// carries no parameters of its own and `dt` is unused, exactly as in the C.
///
/// The C holds `double* min_distance` and `struct reb_orbit* orbit` pointers
/// into the parameter storage and writes through them; here the write-back is
/// the matching `rebx_set_param_*` call, which stores into the same parameter.
pub fn rebx_track_min_distance(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _operator_idx: usize,
    _dt: f64,
) {
    // C: struct rebx_extras* const rebx = sim->extras;  (passed in here)
    let N = sim.N;
    for i in 0..N {
        // C: struct reb_particle* const p = &sim->particles[i];
        let min_distance = rebx_get_param_double(rebx, rebx_ap::particle(i), "min_distance");
        if let Some(min_distance) = min_distance {
            let target = rebx_get_param_string(rebx, rebx_ap::particle(i), "min_distance_from");
            // C: struct reb_particle* source;
            let source = match &target {
                None => Some(0usize),
                Some(target) => reb_simulation_get_particle_by_name(sim, target),
            };
            if let Some(source) = source {
                let p = sim.particles[i];
                let source = sim.particles[source];
                let dx = p.x - source.x;
                let dy = p.y - source.y;
                let dz = p.z - source.z;
                let r2 = dx * dx + dy * dy + dz * dz;
                // C: if (r2 < *min_distance*(*min_distance)) — the parentheses
                // there group the dereference, not the multiplication.
                if r2 < min_distance * min_distance {
                    rebx_set_param_double(rebx, rebx_ap::particle(i), "min_distance", r2.sqrt());
                    let orbit = rebx_get_param_orbit(rebx, rebx_ap::particle(i), "min_distance_orbit");
                    if orbit.is_some() {
                        let orbit = reb_orbit_from_particle(sim.G, p, source);
                        rebx_set_param_orbit(
                            rebx,
                            rebx_ap::particle(i),
                            "min_distance_orbit",
                            orbit,
                        );
                    }
                }
            } else {
                reb_simulation_warning(sim, "min_distance_from cannot find particle");
            }
        }
    }
}
