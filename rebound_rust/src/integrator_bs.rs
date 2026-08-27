//! integrator_bs.rs — the Gragg-Bulirsch-Stoer integration scheme and
//! the general ODE framework (from integrator_bs.c/h; reimplementation
//! of the Hairer & Wanner fortran code via the Hipparchus JAVA
//! implementation; (c) 2021 Hanno Rein, (c) 2004 Ernst Hairer).
//!
//! The C stores `struct reb_ode*` pointers in `r->odes` with weak back
//! references `ode->r`; callbacks receive the ode pointer and reach the
//! simulation through it. In Rust the odes are owned by value in
//! `r.odes` (identified by a per-simulation `id`), and the callback fn
//! pointers receive `&mut reb_simulation` explicitly alongside the ode.
//! During the BS step the odes Vec is temporarily moved out of the
//! simulation so derivatives can mutate particles — pure ownership
//! mechanics, the arithmetic is untouched.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;

/// integrator_bs.c `#define MAX(a, b)` — the exact ternary.
fn MAX(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// integrator_bs.c `#define MIN(a, b)` — the exact ternary.
fn MIN(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

// Default configuration parameters (hard coded in the C as well).
const sequence_length: usize = 9; // = maxOrder / 2
const stepControl1: f64 = 0.65;
const stepControl2: f64 = 0.94;
const stepControl3: f64 = 0.02;
const stepControl4: f64 = 4.0;
const orderControl1: f64 = 0.8;
const orderControl2: f64 = 0.9;
const stabilityReduction: f64 = 0.5;
const maxIter: i32 = 2; // maximal number of iterations for which checks are performed
const maxChecks: i32 = 1; // maximal number of checks for each iteration

/// rebound.h `derivatives` member of `struct reb_ode` (the C signature
/// is `(ode, yDot, y, t)`; `r` is explicit here instead of `ode->r`).
pub type reb_ode_derivatives_fn =
    fn(r: &mut reb_simulation, ode: &mut reb_ode, yDot: &mut [f64], y: &[f64], t: f64);
/// rebound.h `getscale` member (optional; sets `ode.scale`).
pub type reb_ode_getscale_fn =
    fn(r: &mut reb_simulation, ode: &mut reb_ode, y0: &[f64], y1: &[f64]);
/// rebound.h `pre_timestep`/`post_timestep` members (optional).
pub type reb_ode_prepost_fn = fn(r: &mut reb_simulation, ode: &mut reb_ode, y0: &[f64]);

/// rebound.h `struct reb_ode` — one Ordinary Differential Equation set
/// integrated with BS. The C's `ref` (user data pointer) and `r` (weak
/// simulation reference) members have no Rust counterpart; callbacks
/// receive the simulation explicitly.
#[derive(Clone, Debug)]
pub struct reb_ode {
    /// Number of components / dimension.
    pub length: usize,
    /// Current state.
    pub y: Vec<f64>,
    /// 1: ODE needs N-body particles to calculate RHS.
    pub needs_nbody: u32,
    /// Right hand side of the ODE.
    pub derivatives: Option<reb_ode_derivatives_fn>,
    /// Sets scales for components (optional).
    pub getscale: Option<reb_ode_getscale_fn>,
    /// Called just before the ODE integration (optional).
    pub pre_timestep: Option<reb_ode_prepost_fn>,
    /// Called just after the ODE integration (optional).
    pub post_timestep: Option<reb_ode_prepost_fn>,
    // Internal use
    pub scale: Vec<f64>,
    /// Temporary internal array (extrapolation).
    pub C: Vec<f64>,
    /// Temporary internal array (extrapolation), sequence_length rows.
    pub D: Vec<Vec<f64>>,
    /// Temporary internal array (state during the step).
    pub y1: Vec<f64>,
    /// Temporary internal array (derivatives at beginning of step).
    pub y0Dot: Vec<f64>,
    /// Temporary internal array (derivatives).
    pub yDot: Vec<f64>,
    /// Temporary internal array (midpoint method).
    pub yTmp: Vec<f64>,
    /// Rust-side identity (the C identifies odes by pointer).
    pub id: usize,
}

/// integrator_bs.h `struct reb_integrator_bs_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_bs_state {
    /// Allowed absolute scalar error.
    pub eps_abs: f64,
    /// Allowed relative scalar error.
    pub eps_rel: f64,
    /// Minimum timestep.
    pub min_dt: f64,
    /// Maximum timestep.
    pub max_dt: f64,
    // Internal use
    /// id (in `r.odes`) of the ODE corresponding to the N-body system.
    pub nbody_ode: Option<usize>,
    /// Stepsize sequence.
    pub sequence: [i32; sequence_length],
    /// Overall cost of applying step reduction up to iteration k + 1,
    /// in number of calls.
    pub cost_per_step: [i32; sequence_length],
    /// Cost per unit step.
    pub cost_per_time_unit: [f64; sequence_length],
    /// Optimal steps for each order.
    pub optimal_step: [f64; sequence_length],
    /// Extrapolation coefficients.
    pub coeff: [f64; sequence_length],
    pub dt_proposed: f64,
    pub first_or_last_step: i32,
    pub previous_rejected: i32,
    pub target_iter: i32,
    /// Do not set manually. Use needs_nbody in reb_ode instead.
    pub user_ode_needs_nbody: i32,
}

impl Default for reb_integrator_bs_state {
    /// integrator_bs.c `reb_integrator_bs_create`.
    fn default() -> Self {
        let mut bs = reb_integrator_bs_state {
            eps_abs: 1e-8,
            eps_rel: 1e-8,
            max_dt: 0.,
            min_dt: 0.,
            nbody_ode: None,
            sequence: [0; sequence_length],
            cost_per_step: [0; sequence_length],
            cost_per_time_unit: [0.; sequence_length],
            optimal_step: [0.; sequence_length],
            coeff: [0.; sequence_length],
            dt_proposed: 0.,
            first_or_last_step: 1,
            previous_rejected: 0,
            target_iter: 0,
            user_ode_needs_nbody: 0,
        };
        for k in 0..sequence_length {
            // step size sequence: 2, 6, 10, 14, ...
            bs.sequence[k] = 4 * (k as i32) + 2;
            // initialize the extrapolation tables
            let r = 1. / (bs.sequence[k] as f64);
            bs.coeff[k] = r * r;
        }
        // initialize the order selection cost array
        // (number of function calls for each column of the extrapolation table)
        bs.cost_per_step[0] = bs.sequence[0] + 1;
        for k in 1..sequence_length {
            bs.cost_per_step[k] = bs.cost_per_step[k - 1] + bs.sequence[k];
        }
        bs.cost_per_time_unit[0] = 0.;
        bs
    }
}

/// integrator_bs.c `reb_integrator_bs_update_particles`.
pub fn reb_integrator_bs_update_particles(r: &mut reb_simulation, y: &[f64]) {
    for i in 0..r.N {
        let p = &mut r.particles[i];
        p.x = y[i * 6];
        p.y = y[i * 6 + 1];
        p.z = y[i * 6 + 2];
        p.vx = y[i * 6 + 3];
        p.vy = y[i * 6 + 4];
        p.vz = y[i * 6 + 5];
    }
    for i in 0..r.N_var {
        let io = (i + r.N) * 6;
        let p = &mut r.particles_var[i];
        p.x = y[io];
        p.y = y[io + 1];
        p.z = y[io + 2];
        p.vx = y[io + 3];
        p.vy = y[io + 4];
        p.vz = y[io + 5];
    }
}

/// Call one ode's `derivatives` callback with the given target/source
/// buffers taken out of the ode for the duration of the call (the C
/// passes raw pointers into the same struct).
fn call_derivatives(
    r: &mut reb_simulation,
    odes: &mut [reb_ode],
    s: usize,
    target_y0dot: bool, // true: y0Dot <- f(y); false: yDot <- f(y1)
    t: f64,
) {
    let df = match odes[s].derivatives {
        Some(df) => df,
        None => return,
    };
    if target_y0dot {
        let mut yDot = std::mem::take(&mut odes[s].y0Dot);
        let y = std::mem::take(&mut odes[s].y);
        df(r, &mut odes[s], &mut yDot, &y, t);
        odes[s].y0Dot = yDot;
        odes[s].y = y;
    } else {
        let mut yDot = std::mem::take(&mut odes[s].yDot);
        let y = std::mem::take(&mut odes[s].y1);
        df(r, &mut odes[s], &mut yDot, &y, t);
        odes[s].yDot = yDot;
        odes[s].y1 = y;
    }
}

/// integrator_bs.c static `tryStep` — modified midpoint method.
/// Returns 1 on success, 0 if the stability check failed.
fn tryStep(
    r: &mut reb_simulation,
    bs: &mut reb_integrator_bs_state,
    odes: &mut [reb_ode],
    nbody_index: Option<usize>,
    Ns: usize,
    k: i32,
    n: i32,
    t0: f64,
    step: f64,
) -> i32 {
    let subStep = step / (n as f64);
    let mut t = t0;
    let needs_nbody = bs.user_ode_needs_nbody;

    // Modified Midpoint method
    // first substep
    t += subStep;
    for s in 0..Ns {
        let ode = &mut odes[s];
        let length = ode.length;
        for i in 0..length {
            ode.y1[i] = ode.y[i] + subStep * ode.y0Dot[i];
        }
    }

    // other substeps
    if needs_nbody != 0 {
        // Particle do not get updated if user provided ODE does not need N-body data.
        // Also not getting updates if integrator for N-body ODEs is not BS.
        if let Some(ni) = nbody_index {
            let y1 = std::mem::take(&mut odes[ni].y1);
            reb_integrator_bs_update_particles(r, &y1);
            odes[ni].y1 = y1;
        }
    }
    for s in 0..Ns {
        call_derivatives(r, odes, s, false, t);
    }
    for s in 0..Ns {
        let ode = &mut odes[s];
        let length = ode.length;
        for i in 0..length {
            ode.yTmp[i] = ode.y[i];
        }
    }

    for j in 1..n {
        // Note: iterating n substeps, not 2n substeps as in Eq. (9.13)
        t += subStep;
        for s in 0..Ns {
            let ode = &mut odes[s];
            let length = ode.length;
            for i in 0..length {
                let middle = ode.y1[i];
                ode.y1[i] = ode.yTmp[i] + 2. * subStep * ode.yDot[i];
                ode.yTmp[i] = middle;
            }
        }

        if needs_nbody != 0 {
            if let Some(ni) = nbody_index {
                let y1 = std::mem::take(&mut odes[ni].y1);
                reb_integrator_bs_update_particles(r, &y1);
                odes[ni].y1 = y1;
            }
        }
        for s in 0..Ns {
            call_derivatives(r, odes, s, false, t);
        }

        // stability check
        if j <= maxChecks && k < maxIter {
            let mut initialNorm = 0.0;
            let mut deltaNorm = 0.0;
            for s in 0..Ns {
                let ode = &odes[s];
                let length = ode.length;
                for l in 0..length {
                    let ratio1 = ode.y0Dot[l] / ode.scale[l];
                    initialNorm += ratio1 * ratio1;
                    let ratio2 = (ode.yDot[l] - ode.y0Dot[l]) / ode.scale[l];
                    deltaNorm += ratio2 * ratio2;
                }
            }
            if deltaNorm > 4. * MAX(1.0e-15, initialNorm) {
                return 0;
            }
        }
    }

    // correction of the last substep (at t0 + step)
    for s in 0..Ns {
        let ode = &mut odes[s];
        let length = ode.length;
        for i in 0..length {
            ode.y1[i] = 0.5 * (ode.yTmp[i] + ode.y1[i] + subStep * ode.yDot[i]);
        }
    }

    1
}

/// integrator_bs.c static `extrapolate`.
fn extrapolate(ode: &mut reb_ode, coeff: &[f64; sequence_length], k: usize) {
    let length = ode.length;
    for j in 0..k {
        let xi = coeff[k - j - 1];
        let xim1 = coeff[k];
        let facC = xi / (xi - xim1);
        let facD = xim1 / (xi - xim1);
        for i in 0..length {
            let CD = ode.C[i] - ode.D[k - j - 1][i];
            ode.C[i] = facC * CD;
            ode.D[k - j - 1][i] = facD * CD;
        }
    }
    for i in 0..length {
        ode.y1[i] = ode.D[0][i];
    }
    for j in 1..=k {
        for i in 0..length {
            ode.y1[i] += ode.D[j][i];
        }
    }
}

/// integrator_bs.c `reb_integrator_bs_nbody_derivatives`.
pub fn reb_integrator_bs_nbody_derivatives(
    r: &mut reb_simulation,
    _ode: &mut reb_ode,
    yDot: &mut [f64],
    y: &[f64],
    t: f64,
) {
    let t_backup = r.t;
    r.t = t; // Set correct time for time dependent additional forces
    reb_integrator_bs_update_particles(r, y);
    reb_simulation_update_acceleration(r);
    r.t = t_backup;

    for i in 0..r.N {
        let p = r.particles[i];
        yDot[i * 6] = p.vx;
        yDot[i * 6 + 1] = p.vy;
        yDot[i * 6 + 2] = p.vz;
        yDot[i * 6 + 3] = p.ax;
        yDot[i * 6 + 4] = p.ay;
        yDot[i * 6 + 5] = p.az;
    }
    for i in 0..r.N_var {
        let p = r.particles_var[i];
        let io = (i + r.N) * 6;
        yDot[io] = p.vx;
        yDot[io + 1] = p.vy;
        yDot[io + 2] = p.vz;
        yDot[io + 3] = p.ax;
        yDot[io + 4] = p.ay;
        yDot[io + 5] = p.az;
    }
}

/// integrator_bs.c static `reb_integrator_bs_default_scale`. The C
/// passes two state pointers `y1`, `y2` (every current call site passes
/// `ode->y` for both); here the buffers are borrowed field-wise.
fn reb_integrator_bs_default_scale(ode: &mut reb_ode, relTol: f64, absTol: f64) {
    let length = ode.length;
    for i in 0..length {
        // y1 == y2 == ode.y at every call site in the C.
        ode.scale[i] = absTol + relTol * MAX(ode.y[i].abs(), ode.y[i].abs());
    }
}

/// integrator_bs.c `reb_integrator_bs_step_odes` — performs one step on
/// all ODEs in `r.odes`. Gravity not included automatically; it is
/// added in `reb_integrator_bs_step`. Returns 1 if the step was
/// successful, 0 if rejected.
pub fn reb_integrator_bs_step_odes(
    r: &mut reb_simulation,
    bs: &mut reb_integrator_bs_state,
    dt: f64,
) -> i32 {
    let mut odes = std::mem::take(&mut r.odes);
    let ret = reb_integrator_bs_step_odes_inner(r, bs, &mut odes, dt);
    r.odes = odes;
    ret
}

fn reb_integrator_bs_step_odes_inner(
    r: &mut reb_simulation,
    bs: &mut reb_integrator_bs_state,
    odes: &mut Vec<reb_ode>,
    mut dt: f64,
) -> i32 {
    let t = r.t;
    bs.dt_proposed = dt; // In case of early fail

    // initial order selection
    if bs.target_iter == 0 {
        let tol = bs.eps_rel;
        let log10R = MAX(1.0e-10, tol).log10();
        bs.target_iter = std::cmp::max(
            1,
            std::cmp::min(
                (sequence_length as i32) - 2,
                (0.5 - 0.6 * log10R).floor() as i32,
            ),
        );
    }

    let Ns = odes.len(); // Number of ode sets
    let nbody_index: Option<usize> = match bs.nbody_ode {
        Some(id) => odes.iter().position(|o| o.id == id),
        None => None,
    };
    #[allow(unused_assignments)]
    let mut error: f64 = 0.;
    let mut reject = 0;

    // Check if ODEs have been set up correctly
    for s in 0..Ns {
        if odes[s].derivatives.is_none() {
            reb_simulation_error(
                r,
                "A user-specified set of ODEs has not been provided with a derivatives function.",
            );
            r.status = REB_STATUS_GENERIC_ERROR;
            return 0;
        }
    }

    for s in 0..Ns {
        // Check if ODEs need pre timestep setup
        if let Some(pre) = odes[s].pre_timestep {
            let y = std::mem::take(&mut odes[s].y);
            pre(r, &mut odes[s], &y);
            odes[s].y = y;
        }
        // Scaling
        if let Some(gs) = odes[s].getscale {
            // initial scaling: getscale(ode, y, y)
            let y = std::mem::take(&mut odes[s].y);
            gs(r, &mut odes[s], &y, &y);
            odes[s].y = y;
        } else {
            let (eps_rel, eps_abs) = (bs.eps_rel, bs.eps_abs);
            reb_integrator_bs_default_scale(&mut odes[s], eps_rel, eps_abs);
        }
    }

    // first evaluation, at the beginning of the step
    for s in 0..Ns {
        call_derivatives(r, odes, s, true, t);
    }

    let forward = dt >= 0.;

    // iterate over several substep sizes
    let mut k: i32 = -1;
    let mut looping = true;
    while looping {
        k += 1;

        // modified midpoint integration with the current substep
        if tryStep(
            r,
            bs,
            odes,
            nbody_index,
            Ns,
            k,
            bs.sequence[k as usize],
            t,
            dt,
        ) == 0
        {
            // the stability check failed, we reduce the global step
            dt = (dt * stabilityReduction).abs();
            reject = 1;
            looping = false;
        } else {
            for s in 0..Ns {
                let ode = &mut odes[s];
                let length = ode.length;
                for i in 0..length {
                    let CD = ode.y1[i];
                    ode.C[i] = CD;
                    ode.D[k as usize][i] = CD;
                }
            }

            // the substep was computed successfully
            if k > 0 {
                // extrapolate the state at the end of the step
                // using last iteration data
                for s in 0..Ns {
                    let coeff = bs.coeff;
                    extrapolate(&mut odes[s], &coeff, k as usize);
                    if let Some(gs) = odes[s].getscale {
                        let y = std::mem::take(&mut odes[s].y);
                        let y1 = std::mem::take(&mut odes[s].y1);
                        gs(r, &mut odes[s], &y, &y1);
                        odes[s].y = y;
                        odes[s].y1 = y1;
                    } else {
                        let (eps_rel, eps_abs) = (bs.eps_rel, bs.eps_abs);
                        reb_integrator_bs_default_scale(&mut odes[s], eps_rel, eps_abs);
                    }
                }

                // estimate the error at the end of the step.
                error = 0.;
                for s in 0..Ns {
                    let ode = &odes[s];
                    let length = ode.length;
                    for j in 0..length {
                        let e = ode.C[j] / ode.scale[j];
                        error = MAX(error, e * e);
                    }
                }
                // Note: Used to be: error = sqrt(error / combined_length). But for N-body applications it might be more consistent to use:
                error = error.sqrt();
                if error.is_nan() {
                    reb_simulation_error(r, "NaN appearing during ODE integration.");
                    r.status = REB_STATUS_GENERIC_ERROR;
                    return 0;
                }

                if error > 1.0e25 {
                    // error is too big, we reduce the global step
                    dt = (dt * stabilityReduction).abs();
                    reject = 1;
                    looping = false;
                } else {
                    // compute optimal stepsize for this order
                    let exp = 1.0 / ((2 * k + 1) as f64);
                    let mut fac = stepControl2 / (error / stepControl1).powf(exp);
                    let power = stepControl3.powf(exp);
                    fac = MAX(power / stepControl4, MIN(1. / power, fac));
                    bs.optimal_step[k as usize] = (dt * fac).abs();
                    bs.cost_per_time_unit[k as usize] =
                        (bs.cost_per_step[k as usize] as f64) / bs.optimal_step[k as usize];

                    // check convergence
                    match k - bs.target_iter {
                        -1 => {
                            // one before target
                            if bs.target_iter > 1 && bs.previous_rejected == 0 {
                                // check if we can stop iterations now
                                if error <= 1.0 {
                                    // convergence have been reached just before target_iter
                                    looping = false;
                                } else {
                                    // estimate if there is a chance convergence will
                                    // be reached on next iteration, using the
                                    // asymptotic evolution of error
                                    let ratio = ((bs.sequence[bs.target_iter as usize]
                                        * bs.sequence[(bs.target_iter + 1) as usize])
                                        as f64)
                                        / ((bs.sequence[0] * bs.sequence[0]) as f64);
                                    if error > ratio * ratio {
                                        // we don't expect to converge on next iteration
                                        // we reject the step immediately and reduce order
                                        reject = 1;
                                        looping = false;
                                        bs.target_iter = k;
                                        if bs.target_iter > 1
                                            && bs.cost_per_time_unit[(bs.target_iter - 1) as usize]
                                                < orderControl1
                                                    * bs.cost_per_time_unit
                                                        [bs.target_iter as usize]
                                        {
                                            bs.target_iter -= 1;
                                        }
                                        dt = bs.optimal_step[bs.target_iter as usize];
                                    }
                                }
                            }
                        }
                        0 => {
                            // exactly on target
                            if error <= 1.0 {
                                // convergence has been reached exactly at target_iter
                                looping = false;
                            } else {
                                // estimate if there is a chance convergence will
                                // be reached on next iteration, using the
                                // asymptotic evolution of error
                                let ratio = (bs.sequence[(k + 1) as usize] as f64)
                                    / (bs.sequence[0] as f64);
                                if error > ratio * ratio {
                                    // we don't expect to converge on next iteration
                                    // we reject the step immediately
                                    reject = 1;
                                    looping = false;
                                    if bs.target_iter > 1
                                        && bs.cost_per_time_unit[(bs.target_iter - 1) as usize]
                                            < orderControl1
                                                * bs.cost_per_time_unit[bs.target_iter as usize]
                                    {
                                        bs.target_iter -= 1;
                                    }
                                    dt = bs.optimal_step[bs.target_iter as usize];
                                }
                            }
                        }
                        1 => {
                            // one past target
                            if error > 1.0 {
                                reject = 1;
                                if bs.target_iter > 1
                                    && bs.cost_per_time_unit[(bs.target_iter - 1) as usize]
                                        < orderControl1
                                            * bs.cost_per_time_unit[bs.target_iter as usize]
                                {
                                    bs.target_iter -= 1;
                                }
                                dt = bs.optimal_step[bs.target_iter as usize];
                            }
                            looping = false;
                        }
                        _ => {
                            if bs.first_or_last_step != 0 && error <= 1.0 {
                                looping = false;
                            }
                        }
                    }
                }
            }
        }
    }

    if reject == 0 {
        // Swap arrays
        for s in 0..Ns {
            let ode = &mut odes[s];
            std::mem::swap(&mut ode.y, &mut ode.y1);
            // Check if ODEs need post timestep call
            if let Some(post) = odes[s].post_timestep {
                let y = std::mem::take(&mut odes[s].y);
                post(r, &mut odes[s], &y);
                odes[s].y = y;
            }
        }

        let optimalIter: i32;
        if k == 1 {
            optimalIter = if bs.previous_rejected != 0 { 1 } else { 2 };
        } else if k <= bs.target_iter {
            // Converged before or on target
            let mut oi = k;
            if bs.cost_per_time_unit[(k - 1) as usize]
                < orderControl1 * bs.cost_per_time_unit[k as usize]
            {
                oi = k - 1;
            } else if bs.cost_per_time_unit[k as usize]
                < orderControl2 * bs.cost_per_time_unit[(k - 1) as usize]
            {
                oi = std::cmp::min(k + 1, (sequence_length as i32) - 2);
            }
            optimalIter = oi;
        } else {
            // converged after target
            let mut oi = k - 1;
            if k > 2
                && bs.cost_per_time_unit[(k - 2) as usize]
                    < orderControl1 * bs.cost_per_time_unit[(k - 1) as usize]
            {
                oi = k - 2;
            }
            if bs.cost_per_time_unit[k as usize] < orderControl2 * bs.cost_per_time_unit[oi as usize]
            {
                oi = std::cmp::min(k, (sequence_length as i32) - 2);
            }
            optimalIter = oi;
        }

        if bs.previous_rejected != 0 {
            // after a rejected step neither order nor stepsize
            // should increase
            bs.target_iter = std::cmp::min(optimalIter, k);
            dt = MIN(dt.abs(), bs.optimal_step[bs.target_iter as usize]);
        } else {
            // stepsize control
            if optimalIter <= k {
                dt = bs.optimal_step[optimalIter as usize];
            } else {
                if k < bs.target_iter
                    && bs.cost_per_time_unit[k as usize]
                        < orderControl2 * bs.cost_per_time_unit[(k - 1) as usize]
                {
                    dt = bs.optimal_step[k as usize]
                        * (bs.cost_per_step[(optimalIter + 1) as usize] as f64)
                        / (bs.cost_per_step[k as usize] as f64);
                } else {
                    dt = bs.optimal_step[k as usize]
                        * (bs.cost_per_step[optimalIter as usize] as f64)
                        / (bs.cost_per_step[k as usize] as f64);
                }
            }

            bs.target_iter = optimalIter;
        }
    }

    dt = dt.abs();

    if bs.min_dt != 0.0 && dt < bs.min_dt {
        dt = bs.min_dt;
        reb_simulation_warning(r, "Minimal stepsize reached during ODE integration.");
    }

    if bs.max_dt != 0.0 && dt > bs.max_dt {
        dt = bs.max_dt;
        reb_simulation_warning(r, "Maximum stepsize reached during ODE integration.");
    }

    if !forward {
        dt = -dt;
    }
    bs.dt_proposed = dt;

    if reject != 0 {
        bs.previous_rejected = 1;
        0
    } else {
        bs.previous_rejected = 0;
        bs.first_or_last_step = 0;
        1
    }
}

/// integrator_bs.c `reb_ode_create` — allocates an ODE set, attaches it
/// to the simulation and returns its Rust-side id (the C returns the
/// pointer).
pub fn reb_ode_create(r: &mut reb_simulation, length: usize) -> usize {
    let id = r.ode_id_next;
    r.ode_id_next += 1;
    let ode = reb_ode {
        length,
        y: vec![0.; length],
        needs_nbody: 1,
        derivatives: None,
        getscale: None,
        pre_timestep: None,
        post_timestep: None,
        scale: vec![0.; length],
        C: vec![0.; length],
        D: vec![vec![0.; length]; sequence_length],
        y1: vec![0.; length],
        y0Dot: vec![0.; length],
        yDot: vec![0.; length],
        yTmp: vec![0.; length],
        id,
    };
    r.odes.push(ode);

    if let reb_integrator_state::bs(ref mut bs) = r.integrator {
        bs.first_or_last_step = 1;
    }
    id
}

/// integrator_bs.c `reb_ode_free` — detaches the ODE with the given id
/// from the simulation and drops it.
pub fn reb_ode_free(r: &mut reb_simulation, id: usize) {
    if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
        r.odes.remove(pos);
    }
}

/// integrator_bs.c `reb_integrator_bs_step` (state-explicit).
pub fn reb_integrator_bs_step_state(r: &mut reb_simulation, bs: &mut reb_integrator_bs_state) {
    if r.calculate_megno != 0 {
        reb_simulation_error(r, "The BS integrator does currently not support MEGNO.");
    }

    let Ns = r.odes.len();
    for s in 0..Ns {
        let ode = &mut r.odes[s];
        let length = ode.length;
        for i in 0..length {
            ode.y1[i] = ode.y[i];
        }
    }

    let nbody_length = (r.N + r.N_var) * 3 * 2;
    // Check if particle numbers changed, if so delete and recreate ode.
    if let Some(id) = bs.nbody_ode {
        let stale = match r.odes.iter().find(|o| o.id == id) {
            Some(o) => o.length != nbody_length,
            None => true,
        };
        if stale {
            reb_ode_free(r, id);
            bs.nbody_ode = None;
        }
    }
    if bs.nbody_ode.is_none() {
        let id = reb_ode_create(r, nbody_length);
        bs.nbody_ode = Some(id);
        if let Some(ode) = r.odes.iter_mut().find(|o| o.id == id) {
            ode.derivatives = Some(reb_integrator_bs_nbody_derivatives);
            ode.needs_nbody = 0; // No need to update unless there's another ode
        }
        bs.first_or_last_step = 1;
    }

    if let Some(id) = bs.nbody_ode {
        if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
            let mut y = std::mem::take(&mut r.odes[pos].y);
            for i in 0..r.N {
                let p = r.particles[i];
                y[i * 6] = p.x;
                y[i * 6 + 1] = p.y;
                y[i * 6 + 2] = p.z;
                y[i * 6 + 3] = p.vx;
                y[i * 6 + 4] = p.vy;
                y[i * 6 + 5] = p.vz;
            }
            for i in 0..r.N_var {
                let p = r.particles_var[i];
                let io = (i + r.N) * 6;
                y[io] = p.x;
                y[io + 1] = p.y;
                y[io + 2] = p.z;
                y[io + 3] = p.vx;
                y[io + 4] = p.vy;
                y[io + 5] = p.vz;
            }
            r.odes[pos].y = y;
        }
    }

    bs.user_ode_needs_nbody = 0;
    for s in 0..r.odes.len() {
        if r.odes[s].needs_nbody != 0 {
            bs.user_ode_needs_nbody = 1;
        }
    }

    let success = reb_integrator_bs_step_odes(r, bs, r.dt);
    if success != 0 {
        r.t += r.dt;
        r.dt_last_done = r.dt;
    }
    r.dt = bs.dt_proposed;

    if let Some(id) = bs.nbody_ode {
        if let Some(pos) = r.odes.iter().position(|o| o.id == id) {
            let y = std::mem::take(&mut r.odes[pos].y);
            reb_integrator_bs_update_particles(r, &y);
            r.odes[pos].y = y;
        }
    }
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_bs_step(r: &mut reb_simulation) {
    let mut bs = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::bs(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_bs_step_state(r, &mut bs);
    r.integrator = reb_integrator_state::bs(bs);
}
