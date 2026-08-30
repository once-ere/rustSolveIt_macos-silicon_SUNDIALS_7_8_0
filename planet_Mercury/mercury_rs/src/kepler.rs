//! Kepler's equation M = E - e sin E, solved by Newton's method, plus the
//! standard conversions to true anomaly f and Sun-Mercury distance r:
//!
//!   tan(f/2) = sqrt((1+e)/(1-e)) tan(E/2),   r = a (1 - e^2) / (1 + e cos f).
//!
//! The mean anomaly may be any real number; it is reduced to [-pi, pi]
//! internally (every consumer of f uses only 2*pi-periodic expressions).

use std::f64::consts::PI;

/// The solved orbital geometry at one instant.
#[derive(Clone, Copy, Debug)]
pub struct KeplerSolution {
    /// Eccentric anomaly E [rad], on the same branch as the reduced M.
    pub ecc_anom: f64,
    /// True anomaly f [rad] in (-pi, pi].
    pub true_anom: f64,
    /// Sun-Mercury distance r [m].
    pub radius: f64,
}

/// Absolute Newton tolerance on E [rad].
pub const KEPLER_TOL: f64 = 1.0e-14;
/// Iteration cap; non-convergence becomes a named error, never a spin.
pub const KEPLER_MAX_ITER: usize = 50;

/// Solve Kepler's equation for eccentric anomaly, true anomaly, and radius.
pub fn solve(mean_anom: f64, e: f64, a: f64) -> Result<KeplerSolution, String> {
    if !(0.0..1.0).contains(&e) {
        return Err(format!("kepler::solve needs 0 <= e < 1, got e = {e}"));
    }
    let two_pi = 2.0 * PI;
    // Reduce M to [-pi, pi]; sin/cos below then keep full precision.
    let m_red = mean_anom - two_pi * (mean_anom / two_pi).round();

    // Newton iteration, standard starter E0 = M + e sin M.
    let mut ecc_anom = m_red + e * m_red.sin();
    for _ in 0..KEPLER_MAX_ITER {
        let residual = ecc_anom - e * ecc_anom.sin() - m_red;
        let slope = 1.0 - e * ecc_anom.cos();
        let delta = residual / slope;
        ecc_anom -= delta;
        if delta.abs() < KEPLER_TOL {
            let half = 0.5 * ecc_anom;
            // Quadrant-safe form of tan(f/2) = sqrt((1+e)/(1-e)) tan(E/2).
            let true_anom = 2.0
                * ((1.0 + e).sqrt() * half.sin()).atan2((1.0 - e).sqrt() * half.cos());
            let radius = a * (1.0 - e * e) / (1.0 + e * true_anom.cos());
            return Ok(KeplerSolution {
                ecc_anom,
                true_anom,
                radius,
            });
        }
    }
    Err(format!(
        "kepler::solve did not converge in {KEPLER_MAX_ITER} iterations (M = {mean_anom}, e = {e})"
    ))
}
