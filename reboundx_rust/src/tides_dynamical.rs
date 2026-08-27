//! tides_dynamical.rs — translation of REBOUNDx tides_dynamical.c
//! Update body's orbital and modal evolution due to the presence of dynamical tides.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Author: Donald J. Liveoak <dliveoak@umich.edu>
//!
//! # Tides
//!
//! ======================= ===============================================
//! Authors                 D. Liveoak, S. Millholland, M. Vick, D. Tamayo
//! Implementation Paper    Liveoak et al., 2025
//! Based on                Vick et al. 2019
//! C Example               `c_example_tides_dynamical`
//! Python Example          `TidesDynamical.ipynb`
//! ======================= ===============================================
//!
//! This updates body's orbital and modal evolution due to the presence of dynamical tides.
//! Particles are modeled by a gamma=2 polytrope, and the f-mode is evolved at each pericentre passage.
//! The dissipation of orbital energy due to dynamical tides is modeled as an angular
//! momentum-conserving kick at periapse.
//! When mode energy grows to exceed `td_E_max`, it is non-linearly dissipated in one orbital
//! period to `td_E_resid`.
//! To isolate the effects of chaotic model evolution, one can set `dP_hat_crit` to disable
//! dynamical tides whenever chaos is unlikely (see Vick et al. (2019)).
//! Implementation is only applied to particles[1] in the simulation.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! td_disruption_flag (int)     No          Raise error if a planet becomes tidally disrupted (default:0)
//! ============================ =========== ==================================================================
//!
//! **Particle Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! particles[1].m (float)       Yes         Mass
//! particles[1].r (float)       Yes         Physical radius
//! td_E_max (float)             No          Threshold mode energy for non-linear dissipation (default: 0.1 * E_bind)
//! td_E_resid (float)           No          Residual mode energy after non-linear dissipation (default: 0.001 * E_bind)
//! td_c_real (float)            No          Real component of mode (default: 0)
//! td_c_imag (float)            No          Imaginary component of mode (default: 0)
//! td_dP_crit (float)           No          Critical change in mode phase to enable dynamical tides (default: 0)
//! ============================ =========== ==================================================================

use std::f64::consts::PI as M_PI;

use rebound_rs::{
    reb_orbit_from_particle, reb_particle, reb_simulation, reb_simulation_error,
    reb_simulation_warning,
};

use crate::core::{
    rebx_get_param_double, rebx_get_param_int, rebx_set_param_double, rebx_set_param_int,
};
use crate::types::{rebx_ap, rebx_extras};

/// reboundx.h `struct rebx_tides_dynamical_params`.
#[derive(Clone, Copy, Debug)]
pub struct rebx_tides_dynamical_params {
    pub dP: f64,
    pub dE_alpha: f64,
    pub sigma: f64,
}

/// reboundx.h `struct rebx_tides_dynamical_mode`.
///
/// The C struct carries a third member, `char mode`, which the library
/// never writes and never reads. It is kept here for structural fidelity
/// and initialized to 0 (the C leaves it indeterminate).
#[derive(Clone, Copy, Debug)]
pub struct rebx_tides_dynamical_mode {
    pub real: f64,
    pub imag: f64,
    pub mode: u8,
}

/// tides_dynamical.c `rebx_calculate_tides_dynamical_params` (static).
///
/// Deviation from the C signature: the C reads `sim->extras` and
/// `p->ap` off the pointers it is handed. Here `rebx` and the selector
/// `p_ap` naming particle `p`'s parameter list are passed explicitly,
/// and the particles are passed by value (they are only read).
fn rebx_calculate_tides_dynamical_params(
    sim: &mut reb_simulation,
    rebx: &rebx_extras,
    p: reb_particle,
    p_ap: rebx_ap,
    primary: reb_particle,
    raise: i32,
) -> rebx_tides_dynamical_params {
    // The C hard-codes this truncated value; it is NOT f64::consts::E and
    // must not be replaced by it, or the results stop matching bit for bit.
    #[allow(clippy::approx_constant)]
    let EulerConstant: f64 = 2.718281828459;

    // Calculate orbital elements
    let o = reb_orbit_from_particle(sim.G, p, primary);
    let e = o.e;
    let a = o.a;
    let n = o.n;
    let P = o.P;

    // Calculate some useful distances
    let R = p.r; // radius of planet
    let R_tide = R * (primary.m / p.m).powf(0.33333); // tidal radius
    let R_p = a * (1.0 - e); // pericenter distance
    let eta = R_p / R_tide; // pericenter distance in units of tidal radius

    if eta <= 2.5
    // Tidal disruption occurred
    {
        if raise != 0 {
            reb_simulation_error(
                sim,
                "REBOUNDx Error: Planet was disrupted in tides_dynamical.\n",
            );
        }
        reb_simulation_warning(
            sim,
            "Planet was tidally disrupted. No further evolution modeled.\n",
        );
        // C: `struct rebx_tides_dynamical_params toReturn;` with dP,
        // dE_alpha and sigma all set to 0, then returned.
        return rebx_tides_dynamical_params {
            dP: 0.0,
            dE_alpha: 0.0,
            sigma: 0.0,
        };
    }

    // Timescales/frequencies
    let Omega_peri = (sim.G * (p.m + primary.m) / (R_p * R_p * R_p)).powf(0.5); // pericenter frequency
    let time_unit = (sim.G * p.m / (R * R * R)).powf(0.5); // default units for mode parameters

    // Calculate pseudo-synchronous orbital frequency
    // NOTE: (15/2), (45/8), (5/16) and (3/8) are INTEGER divisions in the
    // C source, i.e. 7, 5, 0 and 0. Reproduced digit for digit below.
    let f2 = 1.0 + 7.0 * e.powf(2.0) + 5.0 * e.powf(4.0) + 0.0 * e.powf(6.0);
    let f5 = 1.0 + 3.0 * e.powf(2.0) + 0.0 * e.powf(4.0);
    let Omega_s = n * f2 / ((1.0 - e * e).powf(1.5) * f5);

    // Calculate f-mode parameters, gamma=2 polytrope (see Vick et al. (2019))
    // let omega = (1.22 - Omega_s / time_unit) * time_unit;
    let sigma = (1.22 + Omega_s / time_unit) * time_unit;
    let epsilon = 1.22 * time_unit;
    let Q = 0.56; // overlap integral

    // Calculate K_22 and T
    let z = (2.0_f64).powf(0.5) * sigma / Omega_peri;
    let K_22 = 2.0
        * z.powf(1.5)
        * eta.powf(1.5)
        * EulerConstant.powf(-2.0 * z / 3.0)
        * (1.0 - M_PI.powf(0.5) / (4.0 * z.powf(0.5)))
        / ((15.0_f64).powf(0.5));
    let T = 2.0 * M_PI * M_PI * Q * Q * K_22 * K_22 * sigma / epsilon;

    // Calculate change in mode energy, assuming 0 mode amplitude
    let dE_alpha = sim.G * primary.m * primary.m * R.powf(5.0) * T / R_p.powf(6.0);

    // Calculate dP
    // The C dereferences these three pointers unconditionally; the caller
    // (rebx_tides_dynamical) guarantees all three parameters are set.
    let EB0 = rebx_get_param_double(rebx, p_ap, "td_EB0").unwrap();
    let c_real = rebx_get_param_double(rebx, p_ap, "td_c_real").unwrap();
    let c_imag = rebx_get_param_double(rebx, p_ap, "td_c_imag").unwrap();
    let maxE = dE_alpha
        + 2.0 * (-dE_alpha * (c_real.powf(2.0) + c_imag.powf(2.0)) * EB0).powf(0.5);
    let EBk = -sim.G * p.m * primary.m / (2.0 * a);
    let dP = 1.5 * sigma * P * maxE / (-EBk);

    rebx_tides_dynamical_params {
        dP,
        dE_alpha,
        sigma,
    }
}

/// tides_dynamical.c `rebx_calculate_tides_dynamical_mode_evolution`.
pub fn rebx_calculate_tides_dynamical_mode_evolution(
    old_real: f64,
    old_imag: f64,
    dc_tilde: f64,
    P: f64,
    sigma: f64,
) -> rebx_tides_dynamical_mode {
    let new_real = (old_real + dc_tilde) * (sigma * P).cos() + old_imag * (sigma * P).sin();
    let new_imag = -(old_real + dc_tilde) * (sigma * P).sin() + old_imag * (sigma * P).cos();

    rebx_tides_dynamical_mode {
        real: new_real,
        imag: new_imag,
        mode: 0,
    }
}

/// tides_dynamical.c `rebx_calculate_tides_dynamical_drag_integral`
/// (static): drag integral, eq. 6 of Samsing et al. (2018).
fn rebx_calculate_tides_dynamical_drag_integral(sim: &mut reb_simulation, e: f64, n: f64) -> f64 {
    if n == 10.0 {
        return M_PI
            * (128.0
                + 2944.0 * e * e
                + 10528.0 * e.powf(4.0)
                + 8960.0 * e.powf(6.0)
                + 1715.0 * e.powf(8.0)
                + 35.0 * e.powf(10.0))
            / 128.0;
    }
    if n == 3.0 {
        return M_PI * (1.0 + 2.0 * e * e);
    }
    reb_simulation_error(sim, "REBOUNDx Error: unsupported value for n encountered in rebx_calculate_tides_dynamical_drag_integral().\n");
    0.0
}

/// tides_dynamical.c `rebx_tides_dynamical`.
pub fn rebx_tides_dynamical(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    _N: usize,
) {
    // compute orbit
    // C: `struct reb_particle* const source = &sim->particles[0];` and
    // `p = &sim->particles[1];`. Neither particle's state is modified
    // before the final acceleration writes, so copies read the same
    // values the C's pointers would.
    let source = sim.particles[0];
    let p = sim.particles[1];
    let p_ap = rebx_ap::particle(1);
    let o = reb_orbit_from_particle(sim.G, p, source);

    if p.m == 0.0 || p.r == 0.0 {
        reb_simulation_error(
            sim,
            "REBOUNDx Error: mass and radius must be set for particles[1] in tides_dynamical.\n",
        );
    }

    // Set default parameter values

    let mut raise = 0;
    let raiseptr = rebx_get_param_int(rebx, rebx_ap::force(force_idx), "td_disruption_flag");
    if let Some(raiseptr) = raiseptr {
        raise = raiseptr;
    }

    if rebx_get_param_double(rebx, p_ap, "td_EB0").is_none() {
        let EB0 = -sim.G * p.m * source.m / (2.0 * o.a);
        rebx_set_param_double(rebx, p_ap, "td_EB0", EB0);
    }
    if rebx_get_param_int(rebx, p_ap, "td_num_apoapsis").is_none() {
        rebx_set_param_int(rebx, p_ap, "td_num_apoapsis", 0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_c_real").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_c_real", 0.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_c_imag").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_c_imag", 0.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_dP_crit").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_dP_crit", 0.01);
    }
    if rebx_get_param_double(rebx, p_ap, "td_E_max").is_none() {
        let E_bind = sim.G * p.m * p.m / p.r;
        rebx_set_param_double(rebx, p_ap, "td_E_max", E_bind / 10.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_E_resid").is_none() {
        let E_bind = sim.G * p.m * p.m / p.r;
        rebx_set_param_double(rebx, p_ap, "td_E_resid", E_bind / 1000.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_dP_hat").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_dP_hat", 0.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_drag_coef").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_drag_coef", 0.0);
    }
    if rebx_get_param_double(rebx, p_ap, "td_last_apoapsis").is_none() {
        rebx_set_param_double(rebx, p_ap, "td_last_apoapsis", 0.0);
    }

    let n = 10.0;

    if rebx_get_param_double(rebx, p_ap, "td_M_last").is_some() {
        let M_last = rebx_get_param_double(rebx, p_ap, "td_M_last").unwrap();
        let mut drag = 0.0;
        let last_apoapsis_time = rebx_get_param_double(rebx, p_ap, "td_last_apoapsis").unwrap();
        if (o.M >= M_PI && M_last < M_PI)
            && o.M - M_PI <= 1.0
            && sim.t - last_apoapsis_time >= sim.dt
        {
            // Count apoapsis passages
            let num_apoapsis = rebx_get_param_int(rebx, p_ap, "td_num_apoapsis").unwrap();
            rebx_set_param_int(rebx, p_ap, "td_num_apoapsis", num_apoapsis + 1);

            let dP_crit = rebx_get_param_double(rebx, p_ap, "td_dP_crit").unwrap();
            let dynamical_params =
                rebx_calculate_tides_dynamical_params(sim, rebx, p, p_ap, source, raise);
            let dP = dynamical_params.dP;
            let dE_alpha = dynamical_params.dE_alpha;

            rebx_set_param_double(rebx, p_ap, "td_dP_hat", dP);
            rebx_set_param_double(rebx, p_ap, "td_dE_last", dE_alpha);

            // If system is in chaotic regime, evolve dynamical tides
            if dP >= dP_crit {
                // Calculate map parameters
                let EB0 = rebx_get_param_double(rebx, p_ap, "td_EB0").unwrap();
                // let EBk = -sim.G * p.m * source.m / (2.0 * o.a);
                let dc_tilde = (dE_alpha / -EB0).powf(0.5);
                // let dE_alpha_tilde = dE_alpha / -EB0;
                let c_real = rebx_get_param_double(rebx, p_ap, "td_c_real").unwrap();
                let c_imag = rebx_get_param_double(rebx, p_ap, "td_c_imag").unwrap();

                // Evolve modes
                let sigma = dynamical_params.sigma;
                let mut new_modes = rebx_calculate_tides_dynamical_mode_evolution(
                    c_real, c_imag, dc_tilde, o.P, sigma,
                );
                let dEb = (-EB0)
                    * (new_modes.real * new_modes.real + new_modes.imag * new_modes.imag
                        - c_real * c_real
                        - c_imag * c_imag);

                // If mode energy is too high, non-linear dissipation
                let E_max = rebx_get_param_double(rebx, p_ap, "td_E_max").unwrap();
                let E_resid = rebx_get_param_double(rebx, p_ap, "td_E_resid").unwrap();
                if -(new_modes.real.powf(2.0) + new_modes.imag.powf(2.0)) * EB0 >= E_max {
                    // re-scale modes so that E_mode = E_resid
                    let E_dis_ratio = -E_resid / EB0;
                    new_modes.real = (E_dis_ratio
                        / (1.0 + new_modes.imag.powf(2.0) / new_modes.real.powf(2.0)))
                    .powf(0.5);
                    // NOTE: as in the C, this line uses the *already
                    // overwritten* new_modes.real from the line above.
                    new_modes.imag = (E_dis_ratio
                        / (1.0 + new_modes.real.powf(2.0) / new_modes.imag.powf(2.0)))
                    .powf(0.5);
                }

                rebx_set_param_double(rebx, p_ap, "td_c_real", new_modes.real);
                rebx_set_param_double(rebx, p_ap, "td_c_imag", new_modes.imag);
                rebx_set_param_double(rebx, p_ap, "td_last_apoapsis", sim.t);

                // Compute drag parameter
                let I = rebx_calculate_tides_dynamical_drag_integral(sim, o.e, n);
                drag = dEb * ((o.a) * (1.0 - o.e * o.e)).powf(n - 0.5)
                    / (2.0 * (sim.G * (p.m + source.m)).powf(0.5) * I);
            }
            rebx_set_param_double(rebx, p_ap, "td_drag_coef", drag);
        }
    }

    rebx_set_param_double(rebx, p_ap, "td_M_last", o.M);

    // The C dereferences this unconditionally; the default block above
    // guarantees "td_drag_coef" is set.
    let drag_coef = rebx_get_param_double(rebx, p_ap, "td_drag_coef").unwrap();

    // Compute CoM
    let mut comx = 0.0;
    let mut comy = 0.0;
    let mut comz = 0.0;
    let mut comvx = 0.0;
    let mut comvy = 0.0;
    let mut comvz = 0.0;
    let mut total_m = 0.0;

    for i in 0..2 {
        comx += sim.particles[i].m * sim.particles[i].x;
        comy += sim.particles[i].m * sim.particles[i].y;
        comz += sim.particles[i].m * sim.particles[i].z;
        comvx += sim.particles[i].m * sim.particles[i].vx;
        comvy += sim.particles[i].m * sim.particles[i].vy;
        comvz += sim.particles[i].m * sim.particles[i].vz;
        total_m += sim.particles[i].m;
    }

    let x = sim.particles[1].x - comx / total_m;
    let y = sim.particles[1].y - comy / total_m;
    let z = sim.particles[1].z - comz / total_m;
    let vx = sim.particles[1].vx - comvx / total_m;
    let vy = sim.particles[1].vy - comvy / total_m;
    let vz = sim.particles[1].vz - comvz / total_m;

    let r = (x * x + y * y + z * z).powf(0.5);
    let Fx = -drag_coef * vx / r.powf(n);
    let Fy = -drag_coef * vy / r.powf(n);
    let Fz = -drag_coef * vz / r.powf(n);

    // Apply drag to particle
    let m1 = sim.particles[1].m;
    sim.particles[1].ax += Fx / m1;
    sim.particles[1].ay += Fy / m1;
    sim.particles[1].az += Fz / m1;

    // Apply equal and opposite force to primary
    let m0 = sim.particles[0].m;
    sim.particles[0].ax -= Fx / m0;
    sim.particles[0].ay -= Fy / m0;
    sim.particles[0].az -= Fz / m0;
}
