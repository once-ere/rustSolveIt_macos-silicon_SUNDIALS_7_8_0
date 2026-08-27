//! integrator_whfast.rs — the Wisdom-Holman integrator WHFast (from
//! integrator_whfast.c/h; Rein & Tamayo 2015, correctors of Wisdom et
//! al. 1996, kernels of Rein, Tamayo & Brown 2019). Supports Jacobi,
//! democratic-heliocentric, WHDS and barycentric coordinates, first and
//! second symplectic correctors, all four kernels, safe_mode /
//! keep_unsynchronized, first-order variational equations and MEGNO.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Daniel Tamayo and contributors. See crate root.

use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::{
    reb_simulation_error, reb_simulation_warning, reb_tools_megno_deltad_delta,
    reb_tools_megno_update,
};
use crate::transformations::*;
use crate::types::*;

use std::f64::consts::PI as M_PI;

pub const REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT: u32 = 0;
pub const REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK: u32 = 1;
pub const REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION: u32 = 2;
pub const REB_INTEGRATOR_WHFAST_KERNEL_LAZY: u32 = 3;

pub const REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI: u32 = 0;
pub const REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC: u32 = 1;
pub const REB_INTEGRATOR_WHFAST_COORDINATES_WHDS: u32 = 2;
pub const REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC: u32 = 3;

/// integrator_whfast.h `struct reb_integrator_whfast_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_whfast_state {
    /// Order of first symplectic corrector: 0 (default), 3, 5, 7, 11, 17.
    pub corrector: u32,
    /// 0: no second corrector, 1: use second corrector.
    pub corrector2: u32,
    /// Kernel type (Rein, Tamayo & Brown 2019).
    pub kernel: u32,
    /// Coordinate system used in the Hamiltonian splitting.
    pub coordinates: u32,
    /// 0: DKD scheme, 1 (default): combine first and last sub-step.
    pub safe_mode: u32,
    /// 1: continue from unsynchronized state after synchronization.
    pub keep_unsynchronized: u32,
    // Internal use
    pub p_jh: Vec<reb_particle>,
    pub p_jh_var: Vec<reb_particle>,
    pub recalculate_coordinates_but_not_synchronized_warning: u32,
}

impl Default for reb_integrator_whfast_state {
    /// integrator_whfast.c `reb_integrator_whfast_create`.
    fn default() -> Self {
        reb_integrator_whfast_state {
            corrector: 0,
            corrector2: 0,
            kernel: REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT,
            coordinates: REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI,
            safe_mode: 1,
            keep_unsynchronized: 0,
            p_jh: Vec::new(),
            p_jh_var: Vec::new(),
            recalculate_coordinates_but_not_synchronized_warning: 0,
        }
    }
}

// Corrector coefficients
const reb_whfast_corrector_a_1: f64 = 0.41833001326703777398908601289259374469640768464934;
const reb_whfast_corrector_a_2: f64 = 0.83666002653407554797817202578518748939281536929867;
const reb_whfast_corrector_a_3: f64 = 1.2549900398011133219672580386777812340892230539480;
const reb_whfast_corrector_a_4: f64 = 1.6733200530681510959563440515703749787856307385973;
const reb_whfast_corrector_a_5: f64 = 2.0916500663351888699454300644629687234820384232467;
const reb_whfast_corrector_a_6: f64 = 2.5099800796022266439345160773555624681784461078960;
const reb_whfast_corrector_a_7: f64 = 2.9283100928692644179236020902481562128748537925454;
const reb_whfast_corrector_a_8: f64 = 3.3466401061363021919126881031407499575712614771947;
const reb_whfast_corrector_b_31: f64 = -0.024900596027799867499350357910273437184309981229127;
const reb_whfast_corrector_b_51: f64 = -0.0083001986759332891664501193034244790614366604097090;
const reb_whfast_corrector_b_52: f64 = 0.041500993379666445832250596517122395307183302048545;
const reb_whfast_corrector_b_71: f64 = 0.0024926811426922105779030593952776964450539008582219;
const reb_whfast_corrector_b_72: f64 = -0.018270923246702131478062356884535264841652263842597;
const reb_whfast_corrector_b_73: f64 = 0.053964399093127498721765893493510877532452806339655;
const reb_whfast_corrector_b_111: f64 = 0.00020361579647854651301632818774633716473696537436847;
const reb_whfast_corrector_b_112: f64 = -0.0023487215292295354188307328851055489876255097419754;
const reb_whfast_corrector_b_113: f64 = 0.012309078592019946317544564763237909911330686448336;
const reb_whfast_corrector_b_114: f64 = -0.038121613681288650508647613260247372125243616270670;
const reb_whfast_corrector_b_115: f64 = 0.072593394748842738674253180742744961827622366521517;
const reb_whfast_corrector_b_178: f64 = 0.093056103771425958591541059067553547100903397724386;
const reb_whfast_corrector_b_177: f64 = -0.065192863576377893658290760803725762027864651086787;
const reb_whfast_corrector_b_176: f64 = 0.032422198864713580293681523029577130832258806467604;
const reb_whfast_corrector_b_175: f64 = -0.012071760822342291062449751726959664253913904872527;
const reb_whfast_corrector_b_174: f64 = 0.0033132577069380655655490196833451994080066801611459;
const reb_whfast_corrector_b_173: f64 = -0.00063599983075817658983166881625078545864140848560259;
const reb_whfast_corrector_b_172: f64 = 0.000076436355227935738363241846979413475106795392377415;
const reb_whfast_corrector_b_171: f64 = -0.0000043347415473373580190650223498124944896789841432241;
const reb_whfast_corrector2_b: f64 = 0.03486083443891981449909050107438281205803;

/// Fast inverse factorial lookup table (1/n! for n = 0..34).
const invfactorial: [f64; 35] = [
    1.,
    1.,
    1. / 2.,
    1. / 6.,
    1. / 24.,
    1. / 120.,
    1. / 720.,
    1. / 5040.,
    1. / 40320.,
    1. / 362880.,
    1. / 3628800.,
    1. / 39916800.,
    1. / 479001600.,
    1. / 6227020800.,
    1. / 87178291200.,
    1. / 1307674368000.,
    1. / 20922789888000.,
    1. / 355687428096000.,
    1. / 6402373705728000.,
    1. / 121645100408832000.,
    1. / 2432902008176640000.,
    1. / 51090942171709440000.,
    1. / 1124000727777607680000.,
    1. / 25852016738884976640000.,
    1. / 620448401733239439360000.,
    1. / 15511210043330985984000000.,
    1. / 403291461126605635584000000.,
    1. / 10888869450418352160768000000.,
    1. / 304888344611713860501504000000.,
    1. / 8841761993739701954543616000000.,
    1. / 265252859812191058636308480000000.,
    1. / 8222838654177922817725562880000000.,
    1. / 263130836933693530167218012160000000.,
    1. / 8683317618811886495518194401280000000.,
    1. / 295232799039604140847618609643520000000.,
];

#[inline]
fn fastabs(x: f64) -> f64 {
    if x > 0. {
        x
    } else {
        -x
    }
}

/// integrator_whfast.c `stumpff_cs`.
fn stumpff_cs(cs: &mut [f64; 6], mut z: f64) {
    let mut n: u32 = 0;
    while fastabs(z) > 0.1 {
        z = z / 4.;
        n += 1;
    }
    let nmax = 15usize;
    let mut c_odd = invfactorial[nmax];
    let mut c_even = invfactorial[nmax - 1];
    let mut np = nmax as i32 - 2;
    while np >= 5 {
        c_odd = invfactorial[np as usize] - z * c_odd;
        c_even = invfactorial[np as usize - 1] - z * c_even;
        np -= 2;
    }
    cs[5] = c_odd;
    cs[4] = c_even;
    cs[3] = invfactorial[3] - z * cs[5];
    cs[2] = invfactorial[2] - z * cs[4];
    cs[1] = invfactorial[1] - z * cs[3];
    while n > 0 {
        z = z * 4.;
        cs[5] = (cs[5] + cs[4] + cs[3] * cs[2]) * 0.0625;
        cs[4] = (1. + cs[1]) * cs[3] * 0.125;
        cs[3] = 1. / 6. - z * cs[5];
        cs[2] = 0.5 - z * cs[4];
        cs[1] = 1. - z * cs[3];
        n -= 1;
    }
    cs[0] = invfactorial[0] - z * cs[2];
}

/// integrator_whfast.c `stumpff_cs3`.
fn stumpff_cs3(cs: &mut [f64; 6], mut z: f64) {
    let mut n: u32 = 0;
    while z.abs() > 0.1 {
        z = z / 4.;
        n += 1;
    }
    let nmax = 13usize;
    let mut c_odd = invfactorial[nmax];
    let mut c_even = invfactorial[nmax - 1];
    let mut np = nmax as i32 - 2;
    while np >= 3 {
        c_odd = invfactorial[np as usize] - z * c_odd;
        c_even = invfactorial[np as usize - 1] - z * c_even;
        np -= 2;
    }
    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1] - z * c_odd;
    cs[0] = invfactorial[0] - z * c_even;
    while n > 0 {
        cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
        cs[2] = cs[1] * cs[1] * 0.5;
        cs[1] = cs[0] * cs[1];
        cs[0] = 2. * cs[0] * cs[0] - 1.;
        n -= 1;
    }
}

/// integrator_whfast.c `stiefel_Gs`.
fn stiefel_Gs(Gs: &mut [f64; 6], beta: f64, X: f64) {
    let X2 = X * X;
    stumpff_cs(Gs, beta * X2);
    Gs[1] *= X;
    Gs[2] *= X2;
    let mut _pow = X2 * X;
    Gs[3] *= _pow;
    _pow *= X;
    Gs[4] *= _pow;
    _pow *= X;
    Gs[5] *= _pow;
}

/// integrator_whfast.c `stiefel_Gs3`.
fn stiefel_Gs3(Gs: &mut [f64; 6], beta: f64, X: f64) {
    let X2 = X * X;
    stumpff_cs3(Gs, beta * X2);
    Gs[1] *= X;
    Gs[2] *= X2;
    Gs[3] *= X2 * X;
}

const WHFAST_NMAX_QUART: usize = 64;
const WHFAST_NMAX_NEWT: usize = 32;

#[inline]
fn effective_N_active(r: &reb_simulation) -> usize {
    if r.N_active == usize::MAX || r.testparticle_type == 1 {
        r.N
    } else {
        r.N_active
    }
}

/// integrator_whfast.c `reb_integrator_whfast_kepler_solver` — advances
/// one particle of `whfast.p_jh` (index `pindex`) along a Keplerian
/// orbit for time dt, plus the attached variational particles. The C
/// takes a particle pointer and finds `pindex` by pointer arithmetic;
/// warnings are emitted through `r` (the C casts away const for this).
/// `r` may be `None` (the C passes NULL from TRACE) — then no warning
/// is emitted and no variational particles are advanced.
pub fn reb_integrator_whfast_kepler_solver(
    mut r: Option<&mut reb_simulation>,
    p_jh: &mut [reb_particle],
    p_jh_var: &mut [reb_particle],
    pindex: usize,
    mu: f64,
    dt: f64,
) {
    let p1 = p_jh[pindex]; // Copy of particle

    let r0 = (p1.x * p1.x + p1.y * p1.y + p1.z * p1.z).sqrt();
    let r0i = 1. / r0;
    let v2 = p1.vx * p1.vx + p1.vy * p1.vy + p1.vz * p1.vz;
    let beta = 2. * mu * r0i - v2;
    let eta0 = p1.x * p1.vx + p1.y * p1.vy + p1.z * p1.vz;
    let zeta0 = mu - beta * r0;
    let mut X: f64;
    let mut Gs = [0.0_f64; 6];
    let mut invperiod = 0.; // only used for beta>0
    let mut X_per_period = f64::NAN; // nan triggers Newton's method for beta<0

    if beta > 0. {
        // Elliptic orbit
        let sqrt_beta = beta.sqrt();
        invperiod = sqrt_beta * beta / (2. * M_PI * mu);
        X_per_period = 2. * M_PI / sqrt_beta;
        if dt.abs() * invperiod > 1. {
            if let Some(rr) = r.as_deref_mut() {
                if (rr.messages_timestep_warning & 1) == 0 {
                    rr.messages_timestep_warning |= 1;
                    reb_simulation_warning(
                        rr,
                        "Possible convergence issue. Timestep in Kepler solver is larger than one orbital period.",
                    );
                }
            }
        }
        let dtr0i = dt * r0i;
        X = dtr0i * (1. - dtr0i * eta0 * 0.5 * r0i); // second order guess
    } else {
        // Hyperbolic orbit
        X = 0.; // Initial guess
    }

    let mut converged = false;
    let mut oldX = X;

    // Do one Newton step
    stiefel_Gs3(&mut Gs, beta, X);
    let eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
    let mut ri = 1. / (r0 + eta0Gs1zeta0Gs2);
    X = ri * (X * eta0Gs1zeta0Gs2 - eta0 * Gs[2] - zeta0 * Gs[3] + dt);

    // Choose solver depending on estimated step size
    if fastabs(X - oldX) > 0.01 * X_per_period {
        // Quartic solver, linear initial guess
        X = beta * dt / mu;
        let mut prevX = [0.0_f64; WHFAST_NMAX_QUART + 1];
        let mut n_lag = 1usize;
        while n_lag < WHFAST_NMAX_QUART {
            stiefel_Gs3(&mut Gs, beta, X);
            let f = r0 * X + eta0 * Gs[2] + zeta0 * Gs[3] - dt;
            let fp = r0 + eta0 * Gs[1] + zeta0 * Gs[2];
            let fpp = eta0 * Gs[0] + zeta0 * Gs[1];
            let denom = fp + (16. * fp * fp - 20. * f * fpp).abs().sqrt();
            X = (X * denom - 5. * f) / denom;
            if !X.is_normal() {
                break;
            }
            let mut hit = false;
            for i in 1..n_lag {
                if X == prevX[i] {
                    // Converged. Exit.
                    hit = true;
                    converged = true;
                    break;
                }
            }
            if hit {
                break;
            }
            prevX[n_lag] = X;
            n_lag += 1;
        }
        let eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
        ri = 1. / (r0 + eta0Gs1zeta0Gs2);
    } else {
        // Newton's method
        let mut oldX2;
        let mut n_hg = 1usize;
        while n_hg < WHFAST_NMAX_NEWT {
            oldX2 = oldX;
            oldX = X;
            stiefel_Gs3(&mut Gs, beta, X);
            let eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
            ri = 1. / (r0 + eta0Gs1zeta0Gs2);
            X = ri * (X * eta0Gs1zeta0Gs2 - eta0 * Gs[2] - zeta0 * Gs[3] + dt);
            if !X.is_normal() {
                break;
            }
            if X == oldX || X == oldX2 {
                // Converged. Exit.
                converged = true;
                break;
            }
            n_hg += 1;
        }
    }

    // If solver did not work, fallback to bisection
    if !converged {
        let mut X_min: f64;
        let mut X_max: f64;
        if beta > 0. {
            // Elliptic
            X_min = X_per_period * (dt * invperiod).floor();
            X_max = X_min + X_per_period;
        } else {
            // Hyperbolic
            let h2 = r0 * r0 * v2 - eta0 * eta0;
            let q = h2 / mu / (1. + (1. - h2 * beta / (mu * mu)).sqrt());
            let vq = (h2.sqrt() / q).copysign(dt);
            X_min = dt / (fastabs(vq * dt) + r0);
            X_max = dt / q;
            if dt < 0. {
                std::mem::swap(&mut X_min, &mut X_max);
            }
        }
        X = (X_max + X_min) / 2.;
        loop {
            stiefel_Gs3(&mut Gs, beta, X);
            let s = r0 * X + eta0 * Gs[2] + zeta0 * Gs[3] - dt;
            if s >= 0. {
                X_max = X;
            } else {
                X_min = X;
            }
            X = (X_max + X_min) / 2.;
            // C: `} while (fastabs(X_max-X_min) > fastabs((X_max+X_min)*1e-15));`
            //
            // The negation MUST be written as `!(a > b)`, not as `a <= b`.
            // They differ when either side is NaN, which really happens here:
            // for (near-)rectilinear hyperbolic motion the pericentre q is
            // ~0, so X_max = dt/q is +inf and X_min = dt/(|vq*dt|+r0) is NaN.
            // With NaN, `a > b` is false and the C loop exits after one pass,
            // whereas `a <= b` is ALSO false, so an `if a <= b { break }`
            // never breaks and the solver hangs.
            if !(fastabs(X_max - X_min) > fastabs((X_max + X_min) * 1e-15)) {
                break;
            }
        }
        let eta0Gs1zeta0Gs2 = eta0 * Gs[1] + zeta0 * Gs[2];
        ri = 1. / (r0 + eta0Gs1zeta0Gs2);
    }
    if ri.is_nan() {
        // Exception for (almost) straight line motion in hyperbolic case
        ri = 0.;
        Gs[1] = 0.;
        Gs[2] = 0.;
        Gs[3] = 0.;
    }

    // Note: These are not the traditional f and g functions.
    let f = -mu * Gs[2] * r0i;
    let g = dt - mu * Gs[3];
    let fd = -mu * Gs[1] * r0i * ri;
    let gd = -mu * Gs[2] * ri;

    {
        let p = &mut p_jh[pindex];
        p.x += f * p1.x + g * p1.vx;
        p.y += f * p1.y + g * p1.vy;
        p.z += f * p1.z + g * p1.vz;

        p.vx += fd * p1.x + gd * p1.vx;
        p.vy += fd * p1.y + gd * p1.vy;
        p.vz += fd * p1.z + gd * p1.vz;
    }

    // Variations
    let n_var_config = match r.as_deref() {
        Some(rr) => rr.var_config.len(),
        None => 0,
    };
    for v in 0..n_var_config {
        let vc = match r.as_deref() {
            Some(rr) => rr.var_config[v],
            None => continue,
        };
        let dp1 = p_jh_var[pindex + vc.index];
        stiefel_Gs(&mut Gs, beta, X); // Recalculate (to get Gs[4] and Gs[5])
        let dr0 = (dp1.x * p1.x + dp1.y * p1.y + dp1.z * p1.z) * r0i;
        let dbeta = -2. * mu * dr0 * r0i * r0i
            - 2. * (dp1.vx * p1.vx + dp1.vy * p1.vy + dp1.vz * p1.vz);
        let deta0 = dp1.x * p1.vx + dp1.y * p1.vy + dp1.z * p1.vz
            + p1.x * dp1.vx + p1.y * dp1.vy + p1.z * dp1.vz;
        let dzeta0 = -beta * dr0 - r0 * dbeta;
        let G3beta = 0.5 * (3. * Gs[5] - X * Gs[4]);
        let G2beta = 0.5 * (2. * Gs[4] - X * Gs[3]);
        let G1beta = 0.5 * (Gs[3] - X * Gs[2]);
        let tbeta = eta0 * G2beta + zeta0 * G3beta;
        let dX = -1. * ri * (X * dr0 + Gs[2] * deta0 + Gs[3] * dzeta0 + tbeta * dbeta);
        let dG1 = Gs[0] * dX + G1beta * dbeta;
        let dG2 = Gs[1] * dX + G2beta * dbeta;
        let dG3 = Gs[2] * dX + G3beta * dbeta;
        let dr = dr0 + Gs[1] * deta0 + Gs[2] * dzeta0 + eta0 * dG1 + zeta0 * dG2;
        let df = mu * Gs[2] * dr0 * r0i * r0i - mu * dG2 * r0i;
        let dg = -mu * dG3;
        let dfd = -mu * dG1 * r0i * ri + mu * Gs[1] * (dr0 * r0i + dr * ri) * r0i * ri;
        let dgd = -mu * dG2 * ri + mu * Gs[2] * dr * ri * ri;

        let dp1p = &mut p_jh_var[pindex + vc.index];
        dp1p.x += f * dp1.x + g * dp1.vx + df * p1.x + dg * p1.vx;
        dp1p.y += f * dp1.y + g * dp1.vy + df * p1.y + dg * p1.vy;
        dp1p.z += f * dp1.z + g * dp1.vz + df * p1.z + dg * p1.vz;

        dp1p.vx += fd * dp1.x + gd * dp1.vx + dfd * p1.x + dgd * p1.vx;
        dp1p.vy += fd * dp1.y + gd * dp1.vy + dfd * p1.y + dgd * p1.vy;
        dp1p.vz += fd * dp1.z + gd * dp1.vz + dfd * p1.z + dgd * p1.vz;
    }
}

/// integrator_whfast.c `reb_integrator_whfast_interaction_step`.
pub fn reb_integrator_whfast_interaction_step(
    r: &mut reb_simulation,
    p_jh: &mut [reb_particle],
    p_jh_var: &mut [reb_particle],
    coordinates: u32,
    _dt: f64,
) {
    let N = r.N;
    let N_active = effective_N_active(r);
    let G = r.G;
    let m0 = r.particles[0].m;
    match coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            reb_transformations_inertial_to_jacobi_acc(
                &r.particles,
                p_jh,
                &r.particles,
                N,
                N_active,
            );
            for v in 0..r.var_config.len() {
                let vc = r.var_config[v];
                let (pv, jv) = (&r.particles_var[vc.index..], &mut p_jh_var[vc.index..]);
                reb_transformations_inertial_to_jacobi_acc(pv, jv, &r.particles, N, N_active);
            }
            let mut eta = m0;
            for i in 1..N {
                // Eq 132
                let pji = p_jh[i];
                if i < N_active {
                    eta += pji.m;
                }
                p_jh[i].vx += _dt * pji.ax;
                p_jh[i].vy += _dt * pji.ay;
                p_jh[i].vz += _dt * pji.az;
                if r.gravity != REB_GRAVITY::JACOBI {
                    // Jacobi terms not added in update_acceleration: add here
                    if i > 1 {
                        let rj2i = 1. / (pji.x * pji.x + pji.y * pji.y + pji.z * pji.z);
                        let rji = rj2i.sqrt();
                        let rj3iM = rji * rj2i * G * eta;
                        let prefac1 = _dt * rj3iM;
                        p_jh[i].vx += prefac1 * pji.x;
                        p_jh[i].vy += prefac1 * pji.y;
                        p_jh[i].vz += prefac1 * pji.z;
                        for v in 0..r.var_config.len() {
                            let vc = r.var_config[v];
                            let index = vc.index;
                            let rj5M = rj3iM * rj2i;
                            let pv = p_jh_var[i + index];
                            let rdr = pv.x * pji.x + pv.y * pji.y + pv.z * pji.z;
                            let prefac2 = -_dt * 3. * rdr * rj5M;
                            p_jh_var[i + index].vx += prefac1 * pv.x + prefac2 * pji.x;
                            p_jh_var[i + index].vy += prefac1 * pv.y + prefac2 * pji.y;
                            p_jh_var[i + index].vz += prefac1 * pv.z + prefac2 * pji.z;
                        }
                    }
                    for v in 0..r.var_config.len() {
                        let vc = r.var_config[v];
                        let index = vc.index;
                        let pv = p_jh_var[i + index];
                        p_jh_var[i + index].vx += _dt * pv.ax;
                        p_jh_var[i + index].vy += _dt * pv.ay;
                        p_jh_var[i + index].vz += _dt * pv.az;
                    }
                }
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
            for i in 1..N {
                p_jh[i].vx += _dt * r.particles[i].ax;
                p_jh[i].vy += _dt * r.particles[i].ay;
                p_jh[i].vz += _dt * r.particles[i].az;
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
            for i in 1..N_active {
                let mi = r.particles[i].m;
                p_jh[i].vx += _dt * (m0 + mi) * r.particles[i].ax / m0;
                p_jh[i].vy += _dt * (m0 + mi) * r.particles[i].ay / m0;
                p_jh[i].vz += _dt * (m0 + mi) * r.particles[i].az / m0;
            }
            for i in N_active..N {
                p_jh[i].vx += _dt * r.particles[i].ax;
                p_jh[i].vy += _dt * r.particles[i].ay;
                p_jh[i].vz += _dt * r.particles[i].az;
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            for i in 1..N {
                let pji = p_jh[i];
                let dr = (pji.x * pji.x + pji.y * pji.y + pji.z * pji.z).sqrt();
                let prefac = G * p_jh[0].m / (dr * dr * dr);
                p_jh[i].vx += _dt * (prefac * pji.x + r.particles[i].ax);
                p_jh[i].vy += _dt * (prefac * pji.y + r.particles[i].ay);
                p_jh[i].vz += _dt * (prefac * pji.z + r.particles[i].az);
            }
        }
        _ => {}
    }
}

/// integrator_whfast.c `reb_integrator_whfast_jump_step`.
pub fn reb_integrator_whfast_jump_step(
    r: &reb_simulation,
    p_jh: &mut [reb_particle],
    coordinates: u32,
    _dt: f64,
) {
    let N = r.N;
    let N_active = effective_N_active(r);
    let m0 = r.particles[0].m;
    match coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            // Nothing to be done.
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
            let mut px = 0.;
            let mut py = 0.;
            let mut pz = 0.;
            for i in 1..N_active {
                let m = r.particles[i].m;
                px += m * p_jh[i].vx;
                py += m * p_jh[i].vy;
                pz += m * p_jh[i].vz;
            }
            for i in 1..N {
                p_jh[i].x += _dt * (px / m0);
                p_jh[i].y += _dt * (py / m0);
                p_jh[i].z += _dt * (pz / m0);
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
            let mut px = 0.;
            let mut py = 0.;
            let mut pz = 0.;
            for i in 1..N_active {
                let m = r.particles[i].m;
                px += m * p_jh[i].vx / (m0 + m);
                py += m * p_jh[i].vy / (m0 + m);
                pz += m * p_jh[i].vz / (m0 + m);
            }
            for i in 1..N_active {
                let m = r.particles[i].m;
                let pv = p_jh[i];
                p_jh[i].x += _dt * (px - (m * pv.vx / (m0 + m)));
                p_jh[i].y += _dt * (py - (m * pv.vy / (m0 + m)));
                p_jh[i].z += _dt * (pz - (m * pv.vz / (m0 + m)));
            }
            for i in N_active..N {
                p_jh[i].x += _dt * px;
                p_jh[i].y += _dt * py;
                p_jh[i].z += _dt * pz;
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            // Nothing to be done.
        }
        _ => {}
    }
}

/// integrator_whfast.c `reb_integrator_whfast_kepler_step` (serial).
pub fn reb_integrator_whfast_kepler_step(
    r: &mut reb_simulation,
    p_jh: &mut [reb_particle],
    p_jh_var: &mut [reb_particle],
    coordinates: u32,
    _dt: f64,
) {
    let m0 = r.particles[0].m;
    let G = r.G;
    let N = r.N;
    let N_active = effective_N_active(r);
    let mut eta = m0;
    match coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            for i in 1..N {
                if i < N_active {
                    eta += p_jh[i].m;
                }
                reb_integrator_whfast_kepler_solver(Some(&mut *r), p_jh, p_jh_var, i, eta * G, _dt);
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
            for i in 1..N {
                reb_integrator_whfast_kepler_solver(Some(&mut *r), p_jh, p_jh_var, i, eta * G, _dt); // eta = m0
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
            for i in 1..N {
                if i < N_active {
                    eta = m0 + p_jh[i].m;
                } else {
                    eta = m0;
                }
                reb_integrator_whfast_kepler_solver(Some(&mut *r), p_jh, p_jh_var, i, eta * G, _dt);
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            eta = p_jh[0].m;
            for i in 1..N {
                reb_integrator_whfast_kepler_solver(Some(&mut *r), p_jh, p_jh_var, i, eta * G, _dt);
            }
        }
        _ => {}
    }
}

/// integrator_whfast.c `reb_integrator_whfast_com_step`.
pub fn reb_integrator_whfast_com_step(
    r: &reb_simulation,
    p_jh: &mut [reb_particle],
    p_jh_var: &mut [reb_particle],
    _dt: f64,
) {
    p_jh[0].x += _dt * p_jh[0].vx;
    p_jh[0].y += _dt * p_jh[0].vy;
    p_jh[0].z += _dt * p_jh[0].vz;
    // Only WHFast supports variational equations
    for v in 0..r.var_config.len() {
        let vc = r.var_config[v];
        let pv = p_jh_var[vc.index];
        p_jh_var[vc.index].x += _dt * pv.vx;
        p_jh_var[vc.index].y += _dt * pv.vy;
        p_jh_var[vc.index].z += _dt * pv.vz;
    }
}

/// integrator_whfast.c `reb_whfast_corrector_Z`.
fn reb_whfast_corrector_Z(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    a: f64,
    b: f64,
) {
    let N = r.N;
    let N_active = effective_N_active(r);
    match whfast.coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
            reb_transformations_jacobi_to_inertial_pos_sim(r, whfast, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -b);
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -2. * a);
            reb_transformations_jacobi_to_inertial_pos_sim(r, whfast, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, b);
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
            reb_transformations_barycentric_to_inertial_pos(&mut r.particles, &whfast.p_jh, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -b);
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -2. * a);
            reb_transformations_barycentric_to_inertial_pos(&mut r.particles, &whfast.p_jh, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, b);
            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
        }
        _ => {
            reb_simulation_error(r, "Coordinate system not supported.");
        }
    }
}

/// Jacobi-to-inertial position transform including variational sets
/// (the C repeats this three-line pattern inline).
fn reb_transformations_jacobi_to_inertial_pos_sim(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    N: usize,
    N_active: usize,
) {
    // The C uses `particles` both as the target and as the mass array;
    // the masses are not modified by the transform, so a temporary
    // clone of the particle array serves as the p_mass argument.
    let masses = r.particles.clone();
    reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &whfast.p_jh, &masses, N, N_active);
    for v in 0..r.var_config.len() {
        let vc = r.var_config[v];
        let (pv, jv) = (
            &mut r.particles_var[vc.index..],
            &whfast.p_jh_var[vc.index..],
        );
        reb_transformations_jacobi_to_inertial_pos(pv, jv, &masses, N, N_active);
    }
}

/// integrator_whfast.c `reb_whfast_apply_corrector`.
fn reb_whfast_apply_corrector(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    inv: f64,
    order: u32,
) {
    let dt = r.dt;
    if order == 3 {
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_1 * dt, -inv * reb_whfast_corrector_b_31 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_1 * dt, inv * reb_whfast_corrector_b_31 * dt);
    }
    if order == 5 {
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_2 * dt, -inv * reb_whfast_corrector_b_51 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_1 * dt, -inv * reb_whfast_corrector_b_52 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_1 * dt, inv * reb_whfast_corrector_b_52 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_2 * dt, inv * reb_whfast_corrector_b_51 * dt);
    }
    if order == 7 {
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_3 * dt, -inv * reb_whfast_corrector_b_71 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_2 * dt, -inv * reb_whfast_corrector_b_72 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_1 * dt, -inv * reb_whfast_corrector_b_73 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_1 * dt, inv * reb_whfast_corrector_b_73 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_2 * dt, inv * reb_whfast_corrector_b_72 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_3 * dt, inv * reb_whfast_corrector_b_71 * dt);
    }
    if order == 11 {
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_5 * dt, -inv * reb_whfast_corrector_b_111 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_4 * dt, -inv * reb_whfast_corrector_b_112 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_3 * dt, -inv * reb_whfast_corrector_b_113 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_2 * dt, -inv * reb_whfast_corrector_b_114 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_1 * dt, -inv * reb_whfast_corrector_b_115 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_1 * dt, inv * reb_whfast_corrector_b_115 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_2 * dt, inv * reb_whfast_corrector_b_114 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_3 * dt, inv * reb_whfast_corrector_b_113 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_4 * dt, inv * reb_whfast_corrector_b_112 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_5 * dt, inv * reb_whfast_corrector_b_111 * dt);
    }
    if order == 17 {
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_8 * dt, -inv * reb_whfast_corrector_b_171 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_7 * dt, -inv * reb_whfast_corrector_b_172 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_6 * dt, -inv * reb_whfast_corrector_b_173 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_5 * dt, -inv * reb_whfast_corrector_b_174 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_4 * dt, -inv * reb_whfast_corrector_b_175 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_3 * dt, -inv * reb_whfast_corrector_b_176 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_2 * dt, -inv * reb_whfast_corrector_b_177 * dt);
        reb_whfast_corrector_Z(r, whfast, -reb_whfast_corrector_a_1 * dt, -inv * reb_whfast_corrector_b_178 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_1 * dt, inv * reb_whfast_corrector_b_178 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_2 * dt, inv * reb_whfast_corrector_b_177 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_3 * dt, inv * reb_whfast_corrector_b_176 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_4 * dt, inv * reb_whfast_corrector_b_175 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_5 * dt, inv * reb_whfast_corrector_b_174 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_6 * dt, inv * reb_whfast_corrector_b_173 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_7 * dt, inv * reb_whfast_corrector_b_172 * dt);
        reb_whfast_corrector_Z(r, whfast, reb_whfast_corrector_a_8 * dt, inv * reb_whfast_corrector_b_171 * dt);
    }
}

/// integrator_whfast.c `reb_whfast_operator_C`.
fn reb_whfast_operator_C(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    a: f64,
    b: f64,
) {
    reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
    let N = r.N;
    let N_active = effective_N_active(r);
    {
        let masses = r.particles.clone();
        reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &whfast.p_jh, &masses, N, N_active);
    }
    // Note: variational particles not implemented (as in C).
    reb_simulation_update_acceleration(r);
    reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, b);
    reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -a);
}

/// integrator_whfast.c `reb_whfast_operator_Y`.
fn reb_whfast_operator_Y(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    a: f64,
    b: f64,
) {
    reb_whfast_operator_C(r, whfast, a, b);
    reb_whfast_operator_C(r, whfast, -a, -b);
}

/// integrator_whfast.c `reb_whfast_operator_U`.
fn reb_whfast_operator_U(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    a: f64,
    b: f64,
) {
    reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, a);
    reb_whfast_operator_Y(r, whfast, a, b);
    reb_whfast_operator_Y(r, whfast, a, -b);
    reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -a);
}

/// integrator_whfast.c `reb_whfast_apply_corrector2`.
fn reb_whfast_apply_corrector2(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    inv: f64,
) {
    let a = 0.5 * inv * r.dt;
    let b = reb_whfast_corrector2_b * inv * r.dt;
    reb_whfast_operator_U(r, whfast, a, b);
    reb_whfast_operator_U(r, whfast, -a, b);
}

/// integrator_whfast.c `reb_integrator_whfast_calculate_jerk` — writes
/// the "jerk" into the ax/ay/az fields of `jerk`.
pub fn reb_integrator_whfast_calculate_jerk(r: &reb_simulation, jerk: &mut [reb_particle]) {
    // Assume particles.a calculated.
    let N = r.N;
    let G = r.G;
    let mut Rjx = 0.; // com
    let mut Rjy = 0.;
    let mut Rjz = 0.;
    let mut Mj = 0.;
    let mut Ajx = 0.; // sort of Jacobi acceleration
    let mut Ajy = 0.;
    let mut Ajz = 0.;
    for j in 0..N {
        jerk[j].ax = 0.;
        jerk[j].ay = 0.;
        jerk[j].az = 0.;
        for i in 0..(j + 1) {
            // Jacobi Term (j==1 terms cancel and are skipped, as in C)
            if j > 1 {
                let mut dQkrj = Mj;
                if i < j {
                    dQkrj = -r.particles[j].m;
                }
                let Qkx = r.particles[j].x - Rjx / Mj;
                let Qky = r.particles[j].y - Rjy / Mj;
                let Qkz = r.particles[j].z - Rjz / Mj;
                let dax = r.particles[j].ax - Ajx / Mj;
                let day = r.particles[j].ay - Ajy / Mj;
                let daz = r.particles[j].az - Ajz / Mj;

                let dr = (Qkx * Qkx + Qky * Qky + Qkz * Qkz).sqrt();

                let prefact2 = G * dQkrj / (dr * dr * dr);
                jerk[i].ax += prefact2 * dax;
                jerk[i].ay += prefact2 * day;
                jerk[i].az += prefact2 * daz;

                let alphasum = dax * Qkx + day * Qky + daz * Qkz;
                let prefact1 = 3. * alphasum * prefact2 / (dr * dr);
                jerk[i].ax -= prefact1 * Qkx;
                jerk[i].ay -= prefact1 * Qky;
                jerk[i].az -= prefact1 * Qkz;
            }
            // Direct Term (i==0 && j==1 skipped, cancels)
            if j != i && (i != 0 || j != 1) {
                let dx = r.particles[j].x - r.particles[i].x;
                let dy = r.particles[j].y - r.particles[i].y;
                let dz = r.particles[j].z - r.particles[i].z;

                let dax = r.particles[j].ax - r.particles[i].ax;
                let day = r.particles[j].ay - r.particles[i].ay;
                let daz = r.particles[j].az - r.particles[i].az;

                let dr = (dx * dx + dy * dy + dz * dz).sqrt();
                let alphasum = dax * dx + day * dy + daz * dz;
                let prefact2 = G / (dr * dr * dr);
                let prefact2i = prefact2 * r.particles[i].m;
                let prefact2j = prefact2 * r.particles[j].m;
                jerk[j].ax -= dax * prefact2i;
                jerk[j].ay -= day * prefact2i;
                jerk[j].az -= daz * prefact2i;
                jerk[i].ax += dax * prefact2j;
                jerk[i].ay += day * prefact2j;
                jerk[i].az += daz * prefact2j;
                let prefact1 = 3. * alphasum * prefact2 / (dr * dr);
                let prefact1i = prefact1 * r.particles[i].m;
                let prefact1j = prefact1 * r.particles[j].m;
                jerk[j].ax += dx * prefact1i;
                jerk[j].ay += dy * prefact1i;
                jerk[j].az += dz * prefact1i;
                jerk[i].ax -= dx * prefact1j;
                jerk[i].ay -= dy * prefact1j;
                jerk[i].az -= dz * prefact1j;
            }
        }
        Ajx += r.particles[j].ax * r.particles[j].m;
        Ajy += r.particles[j].ay * r.particles[j].m;
        Ajz += r.particles[j].az * r.particles[j].m;
        Rjx += r.particles[j].x * r.particles[j].m;
        Rjy += r.particles[j].y * r.particles[j].m;
        Rjz += r.particles[j].z * r.particles[j].m;
        Mj += r.particles[j].m;
    }
}

/// integrator_whfast.c `reb_integrator_whfast_init`. Returns 1 on a
/// non-recoverable error (C convention).
pub fn reb_integrator_whfast_init(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
) -> i32 {
    if r.N_var != 0 {
        if whfast.coordinates != REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI {
            reb_simulation_error(r, "Variational particles are only compatible with Jacobi coordinates.");
            return 1;
        }
        if whfast.kernel != REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT {
            reb_simulation_error(r, "Variational particles are only compatible with the standard kernel.");
            return 1;
        }
        if whfast.corrector2 != 0 {
            reb_simulation_error(r, "Variational particles not compatible with 2nd corrector.");
            return 1;
        }
        if whfast.safe_mode == 0 && r.calculate_megno != 0 {
            reb_simulation_error(r, "MEGNO is not compatible with WHFast's safe_mode=0.");
            return 1;
        }
        if whfast.p_jh_var.len() != r.N_var {
            whfast.p_jh_var.resize(r.N_var, reb_particle::default());
        }
        for v in 0..r.var_config.len() {
            let vc = r.var_config[v];
            if vc.order != 1 {
                reb_simulation_error(r, "WHFast only supports first order variational equations.");
                return 1;
            }
            if vc.testparticle >= 0 {
                reb_simulation_error(r, "Test particle variations not supported with WHFast. Use IAS15.");
                return 1;
            }
        }
    }
    if whfast.kernel != REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT
        && whfast.coordinates != REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI
    {
        reb_simulation_error(r, "Non-standard kernel requires Jacobi coordinates.");
        return 1;
    }
    if whfast.kernel > 3 {
        reb_simulation_error(r, "Kernel method must be 0 (default), 1 (exact modified kick), 2 (composition kernel), or 3 (lazy implementer's modified kick). ");
        return 1;
    }
    if whfast.corrector != 0
        && whfast.coordinates != REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI
        && whfast.coordinates != REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC
    {
        reb_simulation_error(r, "Symplectic correctors are only compatible with Jacobi and Barycentric coordinates.");
        return 1;
    }
    if whfast.corrector != 0
        && whfast.corrector != 3
        && whfast.corrector != 5
        && whfast.corrector != 7
        && whfast.corrector != 11
        && whfast.corrector != 17
    {
        reb_simulation_error(r, "First symplectic correctors are only available in the following orders: 0, 3, 5, 7, 11, 17.");
        return 1;
    }
    if whfast.keep_unsynchronized == 1 && whfast.safe_mode == 1 {
        reb_simulation_error(r, "whfast->keep_unsynchronized == 1 is not compatible with safe_mode. Must set whfast->safe_mode = 0.");
    }
    if whfast.kernel == REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK
        || whfast.kernel == REB_INTEGRATOR_WHFAST_KERNEL_LAZY
    {
        r.gravity = REB_GRAVITY::JACOBI;
    } else if whfast.coordinates == REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI {
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1;
    } else if whfast.coordinates == REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC {
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    } else {
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_INVOLVING_0;
    }
    let N = r.N;
    if whfast.p_jh.len() != N {
        whfast.p_jh.resize(N, reb_particle::default());
        r.did_modify_particles = 1;
    }
    0
}

/// integrator_whfast.c `reb_integrator_whfast_from_inertial`.
pub fn reb_integrator_whfast_from_inertial(
    r: &mut reb_simulation,
    p_jh: &mut [reb_particle],
    p_jh_var: &mut [reb_particle],
    coordinates: u32,
) {
    let N = r.N;
    let N_active = effective_N_active(r);
    match coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            let masses = r.particles.clone();
            reb_transformations_inertial_to_jacobi_posvel(&r.particles.clone(), p_jh, &masses, N, N_active);
            for v in 0..r.var_config.len() {
                let vc = r.var_config[v];
                let pv = r.particles_var[vc.index..].to_vec();
                reb_transformations_inertial_to_jacobi_posvel(&pv, &mut p_jh_var[vc.index..], &masses, N, N_active);
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
            reb_transformations_inertial_to_democraticheliocentric_posvel(&r.particles, p_jh, N, N_active);
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
            reb_transformations_inertial_to_whds_posvel(&r.particles, p_jh, N, N_active);
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            reb_transformations_inertial_to_barycentric_posvel(&r.particles, p_jh, N, N_active);
        }
        _ => {}
    }
}

/// integrator_whfast.c `reb_integrator_whfast_to_inertial`. The C's
/// velocity-dependent and independent branches differ only in the
/// variational Jacobi transform (posvel vs pos) — carried exactly.
pub fn reb_integrator_whfast_to_inertial(
    r: &mut reb_simulation,
    p_jh: &[reb_particle],
    p_jh_var: &[reb_particle],
    coordinates: u32,
) {
    let N = r.N;
    let N_active = effective_N_active(r);
    let velocity_dependent = r.force_is_velocity_dependent != 0;
    match coordinates {
        REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
            let masses = r.particles.clone();
            reb_transformations_jacobi_to_inertial_posvel(&mut r.particles, p_jh, &masses, N, N_active);
            for v in 0..r.var_config.len() {
                let vc = r.var_config[v];
                let (pv, jv) = (&mut r.particles_var[vc.index..], &p_jh_var[vc.index..]);
                if velocity_dependent {
                    reb_transformations_jacobi_to_inertial_posvel(pv, jv, &masses, N, N_active);
                } else {
                    reb_transformations_jacobi_to_inertial_pos(pv, jv, &masses, N, N_active);
                }
            }
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
            reb_transformations_democraticheliocentric_to_inertial_posvel(&mut r.particles, p_jh, N, N_active);
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
            reb_transformations_whds_to_inertial_posvel(&mut r.particles, p_jh, N, N_active);
        }
        REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
            reb_transformations_barycentric_to_inertial_posvel(&mut r.particles, p_jh, N, N_active);
        }
        _ => {}
    }
}

/// integrator_whfast.c `reb_integrator_whfast_synchronize` (internal
/// form taking the state explicitly).
pub fn reb_integrator_whfast_synchronize_state(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
) {
    if reb_integrator_whfast_init(r, whfast) != 0 {
        return;
    }
    if r.is_synchronized == 0 {
        let N = r.N;
        let N_active = effective_N_active(r);
        let mut sync_pj: Option<Vec<reb_particle>> = None;
        if whfast.keep_unsynchronized != 0 {
            sync_pj = Some(whfast.p_jh.clone());
        }
        match whfast.kernel {
            REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT
            | REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK
            | REB_INTEGRATOR_WHFAST_KERNEL_LAZY => {
                let half = r.dt / 2.;
                reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, half);
                reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, half);
            }
            REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION => {
                let dt38 = 3. * r.dt / 8.;
                reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt38);
                reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, dt38);
            }
            _ => {
                reb_simulation_error(r, "WHFast kernel not implemented.");
                return;
            }
        }
        if whfast.corrector2 != 0 {
            reb_whfast_apply_corrector2(r, whfast, -1.);
        }
        if whfast.corrector != 0 {
            reb_whfast_apply_corrector(r, whfast, -1., whfast.corrector);
        }
        match whfast.coordinates {
            REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI => {
                let masses = r.particles.clone();
                reb_transformations_jacobi_to_inertial_posvel(&mut r.particles, &whfast.p_jh, &masses, N, N_active);
                for v in 0..r.var_config.len() {
                    let vc = r.var_config[v];
                    let (pv, jv) = (&mut r.particles_var[vc.index..], &whfast.p_jh_var[vc.index..]);
                    reb_transformations_jacobi_to_inertial_posvel(pv, jv, &masses, N, N_active);
                }
            }
            REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC => {
                reb_transformations_democraticheliocentric_to_inertial_posvel(&mut r.particles, &whfast.p_jh, N, N_active);
            }
            REB_INTEGRATOR_WHFAST_COORDINATES_WHDS => {
                reb_transformations_whds_to_inertial_posvel(&mut r.particles, &whfast.p_jh, N, N_active);
            }
            REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC => {
                reb_transformations_barycentric_to_inertial_posvel(&mut r.particles, &whfast.p_jh, N, N_active);
            }
            _ => {}
        }
        if let Some(saved) = sync_pj {
            whfast.p_jh = saved;
        } else {
            r.is_synchronized = 1;
        }
    }
}

/// integrator_whfast.c `reb_integrator_whfast_step` (state-explicit).
pub fn reb_integrator_whfast_step_state(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
) {
    let dt = r.dt;
    let N = r.N;
    let N_active = effective_N_active(r);
    if reb_integrator_whfast_init(r, whfast) != 0 {
        return;
    }

    // Only recalculate Jacobi coordinates if needed
    if whfast.safe_mode != 0 || r.did_modify_particles != 0 {
        if r.is_synchronized == 0 {
            reb_integrator_whfast_synchronize_state(r, whfast);
            if whfast.recalculate_coordinates_but_not_synchronized_warning == 0 {
                reb_simulation_warning(r, "Particles were modified while simulation was not synchronized.");
                whfast.recalculate_coordinates_but_not_synchronized_warning += 1;
            }
        }
        reb_integrator_whfast_from_inertial(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates);
    }
    if r.is_synchronized != 0 {
        // First half DRIFT step
        if whfast.corrector != 0 {
            reb_whfast_apply_corrector(r, whfast, 1., whfast.corrector);
        }
        if whfast.corrector2 != 0 {
            reb_whfast_apply_corrector2(r, whfast, 1.);
        }
        match whfast.kernel {
            REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT
            | REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK
            | REB_INTEGRATOR_WHFAST_KERNEL_LAZY => {
                let half = r.dt / 2.;
                reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, half);
                reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, half);
            }
            REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION => {
                let dt58 = 5. * r.dt / 8.;
                reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt58);
                reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, dt58);
            }
            _ => {
                reb_simulation_error(r, "WHFast kernel not implemented.");
                return;
            }
        }
    } else {
        // Combined DRIFT step
        let fulldt = r.dt;
        reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, fulldt); // full timestep
        reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, fulldt);
    }
    let halfdt = r.dt / 2.;
    reb_integrator_whfast_jump_step(r, &mut whfast.p_jh, whfast.coordinates, halfdt);

    reb_integrator_whfast_to_inertial(r, &whfast.p_jh, &whfast.p_jh_var, whfast.coordinates);

    r.t += dt / 2.;

    reb_simulation_update_acceleration(r);

    match whfast.kernel {
        REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT => {
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt);
            reb_integrator_whfast_jump_step(r, &mut whfast.p_jh, whfast.coordinates, dt / 2.);
        }
        REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK => {
            // p_jh used as a temporary buffer for "jerk"
            {
                let mut jerk = std::mem::take(&mut whfast.p_jh);
                reb_integrator_whfast_calculate_jerk(r, &mut jerk);
                whfast.p_jh = jerk;
            }
            for i in 0..N {
                let prefact = dt * dt / 12.;
                r.particles[i].ax += prefact * whfast.p_jh[i].ax;
                r.particles[i].ay += prefact * whfast.p_jh[i].ay;
                r.particles[i].az += prefact * whfast.p_jh[i].az;
            }
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt);
        }
        REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION => {
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -dt / 6.);

            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -dt / 4.);
            reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, -dt / 4.);

            composition_transform_and_forces(r, whfast, N, N_active);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt / 6.);

            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt / 8.);
            reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, dt / 8.);

            composition_transform_and_forces(r, whfast, N, N_active);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt);

            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -dt / 8.);
            reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, -dt / 8.);

            composition_transform_and_forces(r, whfast, N, N_active);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, -dt / 6.);

            reb_integrator_whfast_kepler_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt / 4.);
            reb_integrator_whfast_com_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, dt / 4.);

            composition_transform_and_forces(r, whfast, N, N_active);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt / 6.);
        }
        REB_INTEGRATOR_WHFAST_KERNEL_LAZY => {
            // Accelerations already calculated. WHT Eq 10.6
            for i in 1..N {
                let prefac1 = dt * dt / 12.;
                r.particles[i].x += prefac1 * r.particles[i].ax;
                r.particles[i].y += prefac1 * r.particles[i].ay;
                r.particles[i].z += prefac1 * r.particles[i].az;
            }
            // Position will be overwritten in next jacobi_to_inertial transformation.

            // recalculate kick
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, &mut whfast.p_jh, &mut whfast.p_jh_var, whfast.coordinates, dt);
        }
        _ => {
            return;
        }
    }

    r.is_synchronized = 0;
    if whfast.safe_mode != 0 {
        reb_integrator_whfast_synchronize_state(r, whfast);
    }

    r.t += r.dt / 2.;
    r.dt_last_done = r.dt;

    if r.calculate_megno != 0 {
        // Need x,v,a synchronized to calculate ddot/d for MEGNO.
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE; // Need all terms.
        crate::gravity::reb_gravity_basic_calculate_acceleration_var(r);
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_INVOLVING_0;

        let dY = r.dt * 2. * (r.t - r.megno_initial_t) * reb_tools_megno_deltad_delta(r);
        reb_tools_megno_update(r, dY, dt);
    }
}

fn composition_transform_and_forces(
    r: &mut reb_simulation,
    whfast: &mut reb_integrator_whfast_state,
    N: usize,
    N_active: usize,
) {
    let masses = r.particles.clone();
    reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &whfast.p_jh, &masses, N, N_active);
    reb_simulation_update_acceleration(r);
}

/// Entry point used by the step dispatcher: takes the state out of the
/// integrator enum for the duration of the step.
pub fn reb_integrator_whfast_step(r: &mut reb_simulation) {
    let mut whfast = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::whfast(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_whfast_step_state(r, &mut whfast);
    r.integrator = reb_integrator_state::whfast(whfast);
}

/// Synchronize entry point for the dispatcher.
pub fn reb_integrator_whfast_synchronize(r: &mut reb_simulation) {
    let mut whfast = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::whfast(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_whfast_synchronize_state(r, &mut whfast);
    r.integrator = reb_integrator_state::whfast(whfast);
}
