//! The five Hut (1981) eccentricity polynomials f1..f5 (the CORRECTED forms —
//! source-spec errata E1/E3 fixed against Hut 1981 itself) and the two torque
//! laws built from them:
//!
//!   tidal brake  <T_tidal> = -K [ Omega f1(e) - n f2(e) ]
//!   handle torque  T_tri   = -(3/2) G M_sun (B-A) / r^3 * sin(2 (theta - f))
//!
//! Every coefficient below is a dyadic rational (exact in f64). The evaluation
//! order is fixed and must not be "simplified" — floating point is not
//! associative.

use crate::params;

/// f1(e) = (1 + 3 e^2 + (3/8) e^4) / (1 - e^2)^(9/2)
pub fn f1(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let num = 1.0 + 3.0 * e2 + 0.375 * e4;
    let om = 1.0 - e2;
    let om2 = om * om;
    num / (om2 * om2 * om.sqrt()) // (1-e^2)^4 * (1-e^2)^(1/2) = (1-e^2)^(9/2)
}

/// f2(e) = (1 + (15/2) e^2 + (45/8) e^4 + (5/16) e^6) / (1 - e^2)^6
pub fn f2(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let e6 = e4 * e2;
    let num = 1.0 + 7.5 * e2 + 5.625 * e4 + 0.3125 * e6;
    let om = 1.0 - e2;
    let om2 = om * om;
    num / (om2 * om2 * om2)
}

/// f3(e) = (1 + (31/2) e^2 + (255/8) e^4 + (185/16) e^6 + (25/64) e^8)
///         / (1 - e^2)^(15/2)
pub fn f3(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let e6 = e4 * e2;
    let e8 = e4 * e4;
    let num = 1.0 + 15.5 * e2 + 31.875 * e4 + 11.5625 * e6 + 0.390625 * e8;
    let om = 1.0 - e2;
    let om2 = om * om;
    let om4 = om2 * om2;
    num / (om4 * om2 * om * om.sqrt()) // (1-e^2)^7 * (1-e^2)^(1/2)
}

/// f4(e) = (1 + (3/2) e^2 + (1/8) e^4) / (1 - e^2)^5
pub fn f4(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let num = 1.0 + 1.5 * e2 + 0.125 * e4;
    let om = 1.0 - e2;
    let om2 = om * om;
    num / (om2 * om2 * om)
}

/// f5(e) = (1 + (15/4) e^2 + (15/8) e^4 + (5/64) e^6) / (1 - e^2)^(13/2)
pub fn f5(e: f64) -> f64 {
    let e2 = e * e;
    let e4 = e2 * e2;
    let e6 = e4 * e2;
    let num = 1.0 + 3.75 * e2 + 1.875 * e4 + 0.078125 * e6;
    let om = 1.0 - e2;
    let om2 = om * om;
    let om4 = om2 * om2;
    num / (om4 * om2 * om.sqrt()) // (1-e^2)^6 * (1-e^2)^(1/2)
}

/// Orbit-averaged tidal torque on the spin [N m]:
/// <T_tidal> = -K (Omega f1 - n f2), K = tidal_k(a, k2tau).
pub fn tidal_torque(k: f64, omega: f64, n: f64, e: f64) -> f64 {
    -k * (omega * f1(e) - n * f2(e))
}

/// Instantaneous triaxial "handle" torque [N m]:
/// T_tri = -(3/2) G M_sun (B-A) / r^3 * sin(2 (theta - f)).
pub fn triaxial_torque(theta: f64, true_anom: f64, radius: f64) -> f64 {
    let r3 = radius * radius * radius;
    -1.5 * params::G * params::M_SUN * params::b_minus_a() / r3
        * (2.0 * (theta - true_anom)).sin()
}

/// The pseudo-synchronous spin ratio f2(e)/f1(e) — where the tidal brake
/// alone would park the spin, in units of the mean motion.
pub fn pseudo_synchronous_ratio(e: f64) -> f64 {
    f2(e) / f1(e)
}
