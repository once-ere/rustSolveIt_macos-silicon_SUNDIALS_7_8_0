//! Physical constants, initial conditions, and solver settings — verbatim from
//! the source specification ("Mercury 3:2 Spin-Orbit Resonant Capture
//! Provenance Specification", consolidated 2026-08-25), plus the documented
//! plan parameters (compression factor, stage thresholds, cadences).
//!
//! Units are SI throughout: meters, kilograms, seconds, radians.

/// Gravitational constant [m^3 kg^-1 s^-2].
pub const G: f64 = 6.67430e-11;
/// Sun mass [kg].
pub const M_SUN: f64 = 1.98847e30;
/// Mercury mass [kg].
pub const M_MERCURY: f64 = 3.3011e23;
/// Mercury mean radius [m].
pub const R_MERCURY: f64 = 2.4397e6;
/// Moment of inertia factor C / (m R^2).
pub const C_FACTOR: f64 = 0.34;
/// Triaxial asymmetry ratio (B - A) / C.
pub const B_MINUS_A_OVER_C: f64 = 1.0e-4;
/// Secular Love number of degree 2 (tidal "squishiness").
pub const K2_LOVE: f64 = 0.12;
/// Tidal constant time lag [s] — the SPEC-LITERAL value.
pub const TAU_SPEC: f64 = 100.0;
/// Documented time-compression factor for the movie runs (plan decision D2).
pub const COMPRESSION_MOVIE: f64 = 1000.0;

/// Initial semi-major axis [m] (~0.387098 AU, today's value).
pub const A0: f64 = 5.790905e10;
/// Initial orbital eccentricity (today's value).
pub const E0: f64 = 0.20563;
/// Initial mean anomaly [rad].
pub const M0: f64 = 0.0;
/// Initial spin angle [rad] (the phase-sweep knob adds offsets to this).
pub const THETA0: f64 = 0.0;
/// Initial fast rotation rate [rad/s] (period ~11.6 hours, ~181x orbital rate).
pub const OMEGA0: f64 = 1.5e-4;

/// CVODE relative tolerance.
pub const REL_TOL: f64 = 1.0e-12;
/// CVODE per-component absolute tolerances for [a, e, M, theta, Omega].
pub const ABS_TOL: [f64; 5] = [1.0e-3, 1.0e-6, 1.0e-10, 1.0e-10, 1.0e-14];
/// CVODE maximum internal step [s] = 10 days (spec).
pub const MAX_STEP: f64 = 864000.0;
/// CVODE maximum internal steps per CVode() call (spec's global budget; each
/// output interval uses far fewer).
pub const MAX_STEPS_PER_CALL: i64 = 500_000_000;
/// Simulation window end [s] = 10 million years (spec).
pub const T_FINAL: f64 = 3.15576e14;
/// One year [s] as the spec defines it (T_FINAL / 1e7).
pub const YEAR: f64 = 3.15576e7;

/// Stage handover: the triaxial torque turns on when Omega/n falls to this.
pub const STAGE_HANDOVER_RATIO: f64 = 2.2;
/// Run B saves the sweep restart state when Omega/n first falls to this.
pub const RESTART_RATIO: f64 = 1.6;
/// Sweep size (plan decision D4).
pub const SWEEP_BRANCHES: usize = 64;

/// Polar moment of inertia C = C_FACTOR * m * R^2 [kg m^2].
pub fn moment_of_inertia() -> f64 {
    C_FACTOR * M_MERCURY * (R_MERCURY * R_MERCURY)
}

/// Permanent equatorial asymmetry (B - A) [kg m^2].
pub fn b_minus_a() -> f64 {
    B_MINUS_A_OVER_C * moment_of_inertia()
}

/// Mean motion n = sqrt(G (M_sun + m) / a^3) [rad/s].
pub fn mean_motion(a: f64) -> f64 {
    (G * (M_SUN + M_MERCURY) / (a * a * a)).sqrt()
}

/// Tidal-brake strength K = 3 G M_sun^2 R^5 k2*tau / a^6 [kg m^2 / s].
/// `k2tau` is the product k2 * tau actually in force for the run.
pub fn tidal_k(a: f64, k2tau: f64) -> f64 {
    let a2 = a * a;
    let a6 = a2 * a2 * a2;
    let r2 = R_MERCURY * R_MERCURY;
    let r5 = r2 * r2 * R_MERCURY;
    3.0 * G * (M_SUN * M_SUN) * r5 * k2tau / a6
}

/// The parameter block handed to CVODE as user_data (Option<Box<dyn Any>>).
#[derive(Clone, Debug)]
pub struct RhsParams {
    /// k2 * tau in force (12.0 spec-literal; 12000.0 for the S=1000 movie).
    pub k2tau: f64,
    /// Whether the triaxial "handle" torque is active (stage R) or averaged
    /// away (stage S).
    pub triaxial_on: bool,
    /// Root function target: CVODE flags where Omega - root_ratio*n = 0.
    /// A value <= 0.0 means the root function is unused.
    pub root_ratio: f64,
}
