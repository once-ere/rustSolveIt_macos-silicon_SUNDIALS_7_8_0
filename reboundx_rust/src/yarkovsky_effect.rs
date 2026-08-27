//! yarkovsky_effect.rs — translation of REBOUNDx yarkovsky_effect.c
//! Adds the perturbations from the Yarkovsky effect to one or more of the orbiting bodies.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Authors: Noah Ferich <nofe4108@colorado.edu>, Dan Tamayo <tamayo.daniel@gmail.com>
//!
//! # $Radiation Forces$
//!
//! ======================= ===============================================
//! Authors                 Noah Ferich, D. Tamayo
//! Implementation Paper    Ferich et al., in prep.
//! Based on                `Veras et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2814V/abstract>`_, `Veras et al., 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..708V/abstract>`_.
//! C Example               :ref:`c_example_yarkovsky_effect`.
//! Python Example          `YarkovskyEffect.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/YarkovskyEffect.ipynb>`_.
//! ======================= ===============================================
//!
//! Adds the accelerations and orbital perturbations created by the Yarkovsky
//! effect onto one or more bodies in the simulation. There are two distinct
//! versions of this effect that can be used: the 'full version' and the
//! 'simple version'. The full version uses the full equations found in Veras
//! et al. (2015) to accurately calculate the Yarkovsky effect on a particle.
//! However, this version slows down simulations and requies a large amount of
//! parameters. For these reasons, the simple version of the effect (based on
//! Veras et al. (2019)) is available. While the magnitude of the acceleration
//! created by the effect will be the same, this version places constant values
//! in a crucial rotation matrix to simplify the push from the Yarkovsky effect
//! on a body. This version is faster and requires less parameters and can be
//! used to get an upper bound on how much the Yarkovsky effect can push an
//! object's orbit inwards or outwards. The lists below describes which
//! parameters are needed for one or both versions of this effect. For more
//! information, please visit the papers and examples linked above.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! ye_lstar (float)             Yes         Luminosity of sim's star (Required for both versions).
//! ye_c (float)                 Yes         Speed of light (Required for both versions).
//! ye_stef_boltz (float)        No          Stefan-Boltzmann constant (Required for full version).
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! particles[i].r (float)       Yes         Physical radius of a body (Required for both versions).
//! ye_flag (int)                Yes         0 sets full version of effect. 1 uses simple version with outward migration. -1 uses the simple version with inward migration (see examples and paper).
//! ye_body_density (float)      Yes         Density of an object (Required for both versions)
//! ye_rotation_period (float)   No          Rotation period of a spinning object (Required for full version)
//! ye_albedo (float)            Yes         Albedo of an object (Reuired for both versions)
//! ye_emissivity (float)        No          Emissivity of an object (Required for full version)
//! ye_thermal_inertia (float)   No          Thermal inertia of an object (Required for full version)
//! ye_k (float)                 No          A constant that gets a value between 0 and 1/4 based on the object's rotation - see Veras et al. (2015) for more information on it (Required for full version)
//! ye_spin_axis_x (float)       No          The x value for the spin axis vector of an object (Required for full version)
//! ye_spin_axis_y (float)       No          The y value for the spin axis vector of an object (Required for full version)
//! ye_spin_axis_z (float)       No          The z value for the spin axis vector of an object (Required for full version)
//! ============================ =========== ==================================================================

use crate::core::{rebx_get_param_double, rebx_get_param_int};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{reb_orbit_from_particle, reb_particle, reb_simulation, reb_simulation_error};

/// C's `M_PI`.
const M_PI: f64 = std::f64::consts::PI;

/// yarkovsky_effect.c `rebx_calculate_yarkovsky_effect` (static).
///
/// The C takes `double*`/`int*` parameter pointers straight out of
/// `rebx_get_param` and tests them against NULL; here they are the `Option`s
/// returned by the `rebx_get_param_*` getters, and `None` is exactly C's NULL.
///
/// `target` is the C's `struct reb_particle* target`: a copy of
/// `particles[i]` that the caller writes back after the call (only
/// `ax`/`ay`/`az` are touched, so this is equivalent to the C's in-place
/// update). `star` is `particles[0]`, read but never written.
#[allow(clippy::too_many_arguments)]
fn rebx_calculate_yarkovsky_effect(
    sim: &mut reb_simulation,
    target: &mut reb_particle,
    star: &reb_particle,
    G: f64,
    density: Option<f64>,
    lstar: Option<f64>,
    rotation_period: Option<f64>,
    Gamma: Option<f64>,
    albedo: Option<f64>,
    emissivity: Option<f64>,
    k: Option<f64>,
    c: Option<f64>,
    stef_boltz: Option<f64>,
    yark_flag: Option<i32>,
    sx: Option<f64>,
    sy: Option<f64>,
    sz: Option<f64>,
) {
    // The C dereferences density, lstar, albedo, c and yark_flag
    // unconditionally: the single caller, rebx_yarkovsky_effect, has already
    // checked all five against NULL before calling. Unpacking them here
    // reproduces that guarantee without an unchecked dereference.
    let (density, lstar, albedo, c, yark_flag) = match (density, lstar, albedo, c, yark_flag) {
        (Some(density), Some(lstar), Some(albedo), Some(c), Some(yark_flag)) => {
            (density, lstar, albedo, c, yark_flag)
        }
        _ => return,
    };

    let unit_matrix: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let radius = target.r;

    let q_yar = 1.0 - albedo;

    let dx = target.x - star.x;
    let dy = target.y - star.y;
    let dz = target.z - star.z;

    let dvx = target.vx - star.vx;
    let dvy = target.vy - star.vy;
    let dvz = target.vz - star.vz;

    // distance of asteroid from the star
    let distance = ((dx * dx) + (dy * dy) + (dz * dz)).sqrt();

    // dot product of position and velocity vectors- the term in the
    // denominator is needed when calculating the i-vector
    let rdotv = ((dx * dvx) + (dy * dvy) + (dz * dvz)) / (c * distance);

    let mut i_vector: [[f64; 1]; 3] = [[0.0], [0.0], [0.0]];

    i_vector[0][0] = ((1.0 - rdotv) * (dx / distance)) - (dvx / c);
    i_vector[1][0] = ((1.0 - rdotv) * (dy / distance)) - (dvy / c);
    i_vector[2][0] = ((1.0 - rdotv) * (dz / distance)) - (dvz / c);

    // magnitude of force created by the effect.
    // The C leaves this uninitialized; it is written by every branch that can
    // reach the code below (yark_flag of 1, -1 or 0 — and the 0 branch returns
    // early if it cannot). A yark_flag outside {-1, 0, 1} reads an
    // indeterminate value in C; here it deterministically reads 0.0.
    let mut yarkovsky_magnitude: f64 = 0.0;

    let mut yark_matrix: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];

    if yark_flag == 1 {
        yark_matrix[1][0] = 1.0; // maximizes the effect pushing outwards

        yarkovsky_magnitude =
            (3.0 * q_yar * lstar) / (64.0 * M_PI * radius * density * c * distance * distance);
    }

    if yark_flag == -1 {
        yark_matrix[0][1] = 1.0; //maximizes the effect pushing inwards

        yarkovsky_magnitude =
            (3.0 * q_yar * lstar) / (64.0 * M_PI * radius * density * c * distance * distance);
    }

    //will run through full equations to create the yark_matrix
    if yark_flag == 0 {
        //makes sure all necessary parameters have been entered
        // (the C also tests `albedo == NULL` here; the caller has already
        // guaranteed it is non-NULL, so that disjunct can never fire.)
        if stef_boltz.is_none()
            || rotation_period.is_none()
            || Gamma.is_none()
            || emissivity.is_none()
            || k.is_none()
            || sx.is_none()
            || sy.is_none()
            || sz.is_none()
        {
            reb_simulation_error(sim, "REBOUNDx Error: One or more parameters missing for this version of the Yarkovsky effect in Rebx. Please make sure you've given values to all variables for this version before running simulations. See documentation and YarkovskyEffect.ipynb. If you'd rather use the simplified version of this effect (requires fewer parameters), then please set 'yark_flag' to -1 or 1.\n\n");
            return;
        }

        let stef_boltz = match stef_boltz {
            Some(v) => v,
            None => return,
        };
        let rotation_period = match rotation_period {
            Some(v) => v,
            None => return,
        };
        let Gamma = match Gamma {
            Some(v) => v,
            None => return,
        };
        let emissivity = match emissivity {
            Some(v) => v,
            None => return,
        };
        let k = match k {
            Some(v) => v,
            None => return,
        };
        let sx = match sx {
            Some(v) => v,
            None => return,
        };
        let sy = match sy {
            Some(v) => v,
            None => return,
        };
        let sz = match sz {
            Some(v) => v,
            None => return,
        };

        let o = reb_orbit_from_particle(G, *target, *star);

        yarkovsky_magnitude =
            (3.0 * k * q_yar * lstar) / (16.0 * M_PI * radius * density * c * distance * distance);

        let Smag = ((sx * sx) + sy * sy + sz * sz).sqrt();

        let hx = (dy * dvz) - (dz * dvy);
        let hy = (dz * dvx) - (dx * dvz);
        let hz = (dx * dvy) - (dy * dvx);
        let Hmag = ((hx * hx) + hy * hy + hz * hz).sqrt();

        let inv_smag = 1.0 / Smag;
        let inv_mag_sqrd = 1.0 / (Smag * Smag);
        let inv_hmag = 1.0 / Hmag;
        let inv_hmag_sqrd = 1.0 / (Hmag * Hmag);

        let R1s: [[f64; 3]; 3] = [
            [0.0, -sz * inv_smag, sy * inv_smag],
            [sz * inv_smag, 0.0, -sx * inv_smag],
            [-sy * inv_smag, sx * inv_smag, 0.0],
        ];

        let R2s: [[f64; 3]; 3] = [
            [
                sx * sx * inv_mag_sqrd,
                sx * sy * inv_mag_sqrd,
                sx * sz * inv_mag_sqrd,
            ],
            [
                sx * sy * inv_mag_sqrd,
                sy * sy * inv_mag_sqrd,
                sy * sz * inv_mag_sqrd,
            ],
            [
                sx * sz * inv_mag_sqrd,
                sy * sz * inv_mag_sqrd,
                sz * sz * inv_mag_sqrd,
            ],
        ];

        let R1h: [[f64; 3]; 3] = [
            [0.0, -hz * inv_hmag, hy * inv_hmag],
            [hz * inv_hmag, 0.0, -hx * inv_hmag],
            [-hy * inv_hmag, hx * inv_hmag, 0.0],
        ];

        let R2h: [[f64; 3]; 3] = [
            [
                hx * hx * inv_hmag_sqrd,
                hx * hy * inv_hmag_sqrd,
                hx * hz * inv_hmag_sqrd,
            ],
            [
                hx * hy * inv_hmag_sqrd,
                hy * hy * inv_hmag_sqrd,
                hy * hz * inv_hmag_sqrd,
            ],
            [
                hx * hz * inv_hmag_sqrd,
                hy * hz * inv_hmag_sqrd,
                hz * hz * inv_hmag_sqrd,
            ],
        ];

        let tanPhi = 1.0
            / (1.0
                + (0.5
                    * ((stef_boltz * emissivity) / (M_PI * M_PI * M_PI * M_PI * M_PI)).powf(0.25))
                    * (rotation_period / (Gamma * Gamma)).sqrt()
                    * ((lstar * q_yar) / (distance * distance)).powf(0.75));

        let tanEpsilon = 1.0
            / (1.0
                + (0.5
                    * ((stef_boltz * emissivity) / (M_PI * M_PI * M_PI * M_PI * M_PI)).powf(0.25))
                    * ((o.P) / (Gamma * Gamma)).sqrt()
                    * ((lstar * q_yar) / (distance * distance)).powf(0.75));

        let Phi = tanPhi.atan();
        let Epsilon = tanEpsilon.atan();

        let cos_phi = Phi.cos();
        let sin_phi = Phi.sin();
        let cos_epsilon = Epsilon.cos();
        let sin_epsilon = Epsilon.sin();

        let mut Rys: [[f64; 3]; 3] = [[0.0; 3]; 3]; //diurnal conntribution for effect

        for i in 0..3 {
            for j in 0..3 {
                Rys[i][j] = (cos_phi * unit_matrix[i][j])
                    + (sin_phi * R1s[i][j])
                    + ((1.0 - cos_phi) * R2s[i][j]);
            }
        }

        let mut Ryh: [[f64; 3]; 3] = [[0.0; 3]; 3];

        for i in 0..3 {
            //seasonal contribution for effect
            for j in 0..3 {
                Ryh[i][j] = (cos_epsilon * unit_matrix[i][j]) - (sin_epsilon * R1h[i][j])
                    + ((1.0 - cos_epsilon) * R2h[i][j]);
            }
        }

        for i in 0..3 {
            for j in 0..3 {
                yark_matrix[i][j] =
                    (Rys[i][0] * Ryh[0][j]) + (Rys[i][1] * Ryh[1][j]) + (Rys[i][2] * Ryh[2][j]);
            }
        }
    }

    let mut direction_matrix: [[f64; 1]; 3] = [[0.0], [0.0], [0.0]];

    //calcuates a vector which gives the direction of the acceleration created by the effect
    for i in 0..3 {
        direction_matrix[i][0] = (yark_matrix[i][0] * i_vector[0][0])
            + (yark_matrix[i][1] * i_vector[1][0])
            + (yark_matrix[i][2] * i_vector[2][0]);
    }

    let mut yarkovsky_acceleration: [[f64; 1]; 3] = [[0.0], [0.0], [0.0]];

    for i in 0..3 {
        //final result for particle's change in acceleration due to the effect
        yarkovsky_acceleration[i][0] = yarkovsky_magnitude * direction_matrix[i][0];
    }

    //adds Yarkovsky aceleration to the asteroid's acceleration in the sim
    target.ax += yarkovsky_acceleration[0][0];
    target.ay += yarkovsky_acceleration[1][0];
    target.az += yarkovsky_acceleration[2][0];
}

/// yarkovsky_effect.c `rebx_yarkovsky_effect`.
///
/// C signature:
/// `void rebx_yarkovsky_effect(struct reb_simulation* const sim,
/// struct rebx_force* const force, struct reb_particle* const particles,
/// const int N)`.
pub fn rebx_yarkovsky_effect(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    N: usize,
) {
    let G = sim.G;

    for i in 1..N {
        // C: struct reb_particle* target = &particles[i];
        //    struct reb_particle* star   = &particles[0];
        let mut target = sim.particles[i];
        let star = sim.particles[0];

        let density = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_body_density");
        let lstar = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ye_lstar");
        let rotation_period =
            rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_rotation_period");
        let Gamma = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_thermal_inertia");
        let albedo = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_albedo");
        let emissivity = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_emissivity");
        let k = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_k");
        let c = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ye_c");
        let stef_boltz = rebx_get_param_double(rebx, rebx_ap::force(force_idx), "ye_stef_boltz");
        let yark_flag = rebx_get_param_int(rebx, rebx_ap::particle(i), "ye_flag");
        let sx = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_spin_axis_x");
        let sy = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_spin_axis_y");
        let sz = rebx_get_param_double(rebx, rebx_ap::particle(i), "ye_spin_axis_z");

        //if these necessary conditions are met the Yarkovsky effect will be calculated for a particle in the sim
        if density.is_some()
            && target.r != 0.0
            && albedo.is_some()
            && lstar.is_some()
            && c.is_some()
            && yark_flag.is_some()
        {
            rebx_calculate_yarkovsky_effect(
                sim,
                &mut target,
                &star,
                G,
                density,
                lstar,
                rotation_period,
                Gamma,
                albedo,
                emissivity,
                k,
                c,
                stef_boltz,
                yark_flag,
                sx,
                sy,
                sz,
            );
            // The C mutates particles[i] through `target`; the only fields it
            // touches are ax/ay/az, written back here.
            sim.particles[i].ax = target.ax;
            sim.particles[i].ay = target.ay;
            sim.particles[i].az = target.az;
        }
    }
}
