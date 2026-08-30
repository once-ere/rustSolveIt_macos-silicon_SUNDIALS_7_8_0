//! The five-equation right-hand side handed to CVODE, and the root function
//! used to stop exactly at a chosen spin/orbit ratio.
//!
//! State vector y = [a, e, M, theta, Omega]:
//!   da/dt     = (2K/(m n a))    [ Omega f2 - n f3 ]
//!   de/dt     = (9Ke/(m n a^2)) [ (11/18) Omega f4 - n f5 ]
//!   dM/dt     = n
//!   dtheta/dt = Omega
//!   dOmega/dt = ( T_tri + <T_tidal> ) / C
//! with n = sqrt(G(M_sun+m)/a^3), K = 3 G M_sun^2 R^5 k2 tau / a^6,
//! C = 0.34 m R^2, and T_tri active only when the parameter block says so
//! (stage R). These are the CORRECTED Hut (1981) forms (source-spec errata
//! E1-E3 fixed); an exact identity d(C*Omega)/dt + dL_orb/dt = 0 holds for
//! the secular part and is unit-tested.

use std::any::Any;

use cvode_rs::prelude::*;

use crate::hut;
use crate::kepler;
use crate::params;
use crate::params::RhsParams;

/// The CVRhsFn given to CVodeInit. Returns 0 on success, +1 (recoverable —
/// CVODE retries with a smaller step) on a bad intermediate state, -1
/// (unrecoverable) if the parameter block is missing.
pub fn rhs(
    _t: f64,
    y: &N_Vector,
    ydot: &N_Vector,
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let p = match user_data.as_mut().and_then(|b| b.downcast_mut::<RhsParams>()) {
        Some(p) => p,
        None => return -1,
    };
    let s = {
        let d = match N_VGetArrayPointer(y) {
            Some(d) => d,
            None => return -1,
        };
        [d[0], d[1], d[2], d[3], d[4]]
    };
    let (a, e, m_anom, theta, omega) = (s[0], s[1], s[2], s[3], s[4]);
    if !(a > 0.0) || !(0.0..1.0).contains(&e) || !omega.is_finite() {
        return 1; // recoverable: let CVODE shrink the step
    }

    let n = params::mean_motion(a);
    let k = params::tidal_k(a, p.k2tau);
    let f2v = hut::f2(e);
    let f3v = hut::f3(e);
    let f4v = hut::f4(e);
    let f5v = hut::f5(e);

    let da = (2.0 * k / (params::M_MERCURY * n * a)) * (omega * f2v - n * f3v);
    let de = (9.0 * k * e / (params::M_MERCURY * n * (a * a)))
        * ((11.0 / 18.0) * omega * f4v - n * f5v);
    let dm = n;
    let dth = omega;

    let mut torque = hut::tidal_torque(k, omega, n, e);
    if p.triaxial_on {
        match kepler::solve(m_anom, e, a) {
            Ok(sol) => torque += hut::triaxial_torque(theta, sol.true_anom, sol.radius),
            Err(_) => return 1, // recoverable
        }
    }
    let dom = torque / params::moment_of_inertia();

    let mut o = match N_VGetArrayPointer(ydot) {
        Some(d) => d,
        None => return -1,
    };
    o[0] = da;
    o[1] = de;
    o[2] = dm;
    o[3] = dth;
    o[4] = dom;
    0
}

/// The CVRootFn given to CVodeRootInit: g(t, y) = Omega - root_ratio * n(a).
/// CVODE reports exactly where the spin/orbit ratio crosses the target.
pub fn ratio_root(
    _t: f64,
    y: &N_Vector,
    gout: &mut [f64],
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let p = match user_data.as_mut().and_then(|b| b.downcast_mut::<RhsParams>()) {
        Some(p) => p,
        None => return -1,
    };
    let (a, omega) = {
        let d = match N_VGetArrayPointer(y) {
            Some(d) => d,
            None => return -1,
        };
        (d[0], d[4])
    };
    if !(a > 0.0) {
        return -1;
    }
    gout[0] = omega - p.root_ratio * params::mean_motion(a);
    0
}
