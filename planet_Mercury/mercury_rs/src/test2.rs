//! TEST 2 — Jupiter + Einstein. The six-state extension of the Mercury model:
//!
//!   y = [ a, e, M, pomega, theta, Omega ]
//!
//! adding the perihelion longitude pomega, which precesses under
//! (1) Einstein's general-relativistic advance
//!         dpomega/dt|GR = 3 n G M_sun / (c^2 a (1 - e^2))    (~43"/century)
//! and (2) Jupiter's Laplace-Lagrange secular terms
//!         dpomega/dt|LL = A11 + A12 (e_J/e) cos(pomega - pomega_J)
//!         de/dt|LL      = A12 e_J sin(pomega - pomega_J)
//! with A11 = +n (1/4)(m_J/M_sun) alpha^2 b32_1(alpha),
//!      A12 = -n (1/4)(m_J/M_sun) alpha^2 b32_2(alpha),
//! alpha = a/a_J, and b32_j the Laplace coefficients b_{3/2}^{(j)}.
//! Jupiter's own orbit is held fixed (its back-reaction from Mercury is
//! negligible), so Mercury's eccentricity oscillates by +/- |A12/A11| e_J
//! with period 2 pi / A11 — the clean two-frequency version of the
//! Correia-Laskar mechanism. The triaxial torque argument becomes
//! 2 (theta - f - pomega) and the resonance angle gamma2 = 2 theta - 3 M -
//! 2 pomega, so the locked mean spin ratio sits at 1.5 + pomega_dot/n: the
//! lock follows the precessing ellipse, not the stars.
//!
//! Test 1's five-state code path is deliberately untouched (its results are
//! published and bit-frozen); this module carries its own six-state RHS,
//! segment integrator, detector, and writers. Angular momentum is NOT
//! ledger-checked here: the Laplace-Lagrange terms exchange orbital angular
//! momentum with Jupiter, which this model does not track (documented).
//!
//! NOTE on the movie clock: the tidal strength keeps test 1's documented
//! 1000x compression, while the GR and Jupiter precession rates are REAL —
//! so Mercury experiences ~12 eccentricity cycles during the compressed
//! braking instead of the ~5750 of the uncompressed history. The
//! qualitative physics (eccentricity varying across the resonance crossing)
//! is present; its statistics are not to scale, and every document says so.

use std::any::Any;
use std::collections::VecDeque;
use std::f64::consts::PI;

use cvode_rs::prelude::*;
use sundials_core::sundials_utils::fmt_e;

use crate::driver::{wrap_pi, Event};
use crate::hut;
use crate::kepler;
use crate::params;

// --- Jupiter and relativity constants (SI) --------------------------------
/// Jupiter mass [kg].
pub const M_JUP: f64 = 1.89813e27;
/// Jupiter semi-major axis [m] (5.2038 AU).
pub const A_JUP: f64 = 7.78479e11;
/// Jupiter orbital eccentricity.
pub const E_JUP: f64 = 0.0489;
/// Jupiter perihelion longitude [rad] — the fixed reference direction.
pub const POMEGA_JUP: f64 = 0.0;
/// Speed of light [m/s].
pub const C_LIGHT: f64 = 2.99792458e8;
/// Radians/second -> arcseconds/century.
pub const RAD_S_TO_ARCSEC_CY: f64 = 206264.80624709636 * 3.15576e9;

/// Laplace coefficient b_{3/2}^{(j)}(alpha) by 4096-point trapezoid over the
/// periodic integrand (spectrally accurate, fully deterministic):
/// b = (1/pi) * INT_0^{2pi} cos(j psi) / (1 - 2 a cos psi + a^2)^{3/2} dpsi.
pub fn laplace_b32(j: i32, alpha: f64) -> f64 {
    let n = 4096usize;
    let mut sum = 0.0f64;
    for i in 0..n {
        let psi = 2.0 * PI * (i as f64) / (n as f64);
        let den = 1.0 - 2.0 * alpha * psi.cos() + alpha * alpha;
        sum += ((j as f64) * psi).cos() / (den * den * den).sqrt();
    }
    sum * (2.0 * PI / (n as f64)) / PI
}

/// The two Laplace-Lagrange secular rates [rad/s], computed at a0.
#[derive(Clone, Copy, Debug)]
pub struct Ll {
    pub a11: f64,
    pub a12: f64,
}

pub fn ll_rates() -> Ll {
    let n0 = params::mean_motion(params::A0);
    let alpha = params::A0 / A_JUP;
    let pref = n0 * 0.25 * (M_JUP / params::M_SUN) * (alpha * alpha);
    Ll {
        a11: pref * laplace_b32(1, alpha),
        a12: -pref * laplace_b32(2, alpha),
    }
}

/// Einstein's apsidal advance [rad/s] at (a, e).
pub fn gr_pomega_dot(a: f64, e: f64) -> f64 {
    let n = params::mean_motion(a);
    3.0 * n * params::G * params::M_SUN / (C_LIGHT * C_LIGHT * a * (1.0 - e * e))
}

// --- the six-state system --------------------------------------------------

#[derive(Clone, Debug)]
pub struct Rhs2Params {
    pub k2tau: f64,
    pub triaxial_on: bool,
    pub gr_on: bool,
    pub jupiter_on: bool,
    pub a11: f64,
    pub a12: f64,
    pub root_ratio: f64,
}

pub fn rhs6(
    _t: f64,
    y: &N_Vector,
    ydot: &N_Vector,
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let p = match user_data.as_mut().and_then(|b| b.downcast_mut::<Rhs2Params>()) {
        Some(p) => p,
        None => return -1,
    };
    let s = {
        let d = match N_VGetArrayPointer(y) {
            Some(d) => d,
            None => return -1,
        };
        [d[0], d[1], d[2], d[3], d[4], d[5]]
    };
    let (a, e, m_anom, pw, theta, omega) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    if !(a > 0.0) || !(0.0..1.0).contains(&e) || e < 1.0e-6 || !omega.is_finite() {
        return 1;
    }

    let n = params::mean_motion(a);
    let k = params::tidal_k(a, p.k2tau);
    let f2v = hut::f2(e);
    let f3v = hut::f3(e);
    let f4v = hut::f4(e);
    let f5v = hut::f5(e);

    let da = (2.0 * k / (params::M_MERCURY * n * a)) * (omega * f2v - n * f3v);
    let mut de = (9.0 * k * e / (params::M_MERCURY * n * (a * a)))
        * ((11.0 / 18.0) * omega * f4v - n * f5v);
    let mut dpw = 0.0;
    if p.gr_on {
        dpw += gr_pomega_dot(a, e);
    }
    if p.jupiter_on {
        de += p.a12 * E_JUP * (pw - POMEGA_JUP).sin();
        dpw += p.a11 + p.a12 * (E_JUP / e) * (pw - POMEGA_JUP).cos();
    }
    let dm = n;
    let dth = omega;

    let mut torque = hut::tidal_torque(k, omega, n, e);
    if p.triaxial_on {
        match kepler::solve(m_anom, e, a) {
            Ok(sol) => {
                torque += hut::triaxial_torque(theta, sol.true_anom + pw, sol.radius);
            }
            Err(_) => return 1,
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
    o[3] = dpw;
    o[4] = dth;
    o[5] = dom;
    0
}

pub fn ratio_root6(
    _t: f64,
    y: &N_Vector,
    gout: &mut [f64],
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let p = match user_data.as_mut().and_then(|b| b.downcast_mut::<Rhs2Params>()) {
        Some(p) => p,
        None => return -1,
    };
    let (a, omega) = {
        let d = match N_VGetArrayPointer(y) {
            Some(d) => d,
            None => return -1,
        };
        (d[0], d[5])
    };
    if !(a > 0.0) {
        return -1;
    }
    gout[0] = omega - p.root_ratio * params::mean_motion(a);
    0
}

/// One instant of the six-variable system.
#[derive(Clone, Copy, Debug)]
pub struct State2 {
    pub t: f64,
    /// y = [a, e, M, pomega, theta, Omega]
    pub y: [f64; 6],
}

impl State2 {
    /// Same re-anchoring as test 1 (M -= 2 pi j, theta -= 3 pi j; pomega is
    /// small and untouched) — gamma2 = 2 theta - 3 M - 2 pomega changes by
    /// 2(-3 pi j) - 3(-2 pi j) = 0 exactly.
    pub fn reanchored(&self) -> State2 {
        let two_pi = 2.0 * PI;
        let j = (self.y[2] / two_pi).floor();
        let mut y = self.y;
        y[2] -= two_pi * j;
        y[4] -= 3.0 * PI * j;
        State2 { t: self.t, y }
    }
}

/// One output sample of the six-state system.
#[derive(Clone, Copy, Debug)]
pub struct Sample2 {
    pub t: f64,
    pub a: f64,
    pub e: f64,
    pub m_anom: f64,
    pub pomega: f64,
    pub theta: f64,
    pub omega: f64,
    pub n: f64,
    pub ratio: f64,
    pub gamma2: f64,
    pub p_orb: f64,
    pub p_rot: f64,
    pub l_spin: f64,
    pub l_orb: f64,
    pub l_tot: f64,
    pub e_spin: f64,
    pub e_orb: f64,
    pub stage: char,
}

impl Sample2 {
    pub fn from_state(t: f64, y: &[f64; 6], stage: char) -> Sample2 {
        let (a, e, m_anom, pomega, theta, omega) = (y[0], y[1], y[2], y[3], y[4], y[5]);
        let n = params::mean_motion(a);
        let c = params::moment_of_inertia();
        let l_spin = c * omega;
        let l_orb = params::M_MERCURY * n * (a * a) * (1.0 - e * e).sqrt();
        Sample2 {
            t,
            a,
            e,
            m_anom,
            pomega,
            theta,
            omega,
            n,
            ratio: omega / n,
            gamma2: wrap_pi(2.0 * theta - 3.0 * m_anom - 2.0 * pomega),
            p_orb: 2.0 * PI / n,
            p_rot: 2.0 * PI / omega,
            l_spin,
            l_orb,
            l_tot: l_spin + l_orb,
            e_spin: 0.5 * c * (omega * omega),
            e_orb: -params::G * params::M_SUN * params::M_MERCURY / (2.0 * a),
            stage,
        }
    }
}

pub enum Cmd2 {
    Continue,
    SetCadence(f64),
    Stop,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Stats2 {
    pub n_steps: i64,
    pub n_rhs: i64,
    pub n_reanchor: i64,
}

pub struct Segment2 {
    pub p: Rhs2Params,
    pub t_end: f64,
    pub cadence: f64,
    pub stop_on_root: bool,
    pub reanchor: bool,
    pub stage_tag: char,
    pub record: bool,
}

pub struct Segment2End {
    pub state: State2,
    pub root_state: Option<State2>,
}

fn with_data_mut<R>(v: &N_Vector, f: impl FnOnce(&mut [f64]) -> R) -> Option<R> {
    let mut d = N_VGetArrayPointer(v)?;
    Some(f(&mut d))
}

fn read_state6(v: &N_Vector) -> Result<[f64; 6], String> {
    let d = N_VGetArrayPointer(v)
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    Ok([d[0], d[1], d[2], d[3], d[4], d[5]])
}

fn harvest2(cv: &CVodeMem, stats: &mut Stats2) -> Result<(), String> {
    let mut s: i64 = 0;
    let mut r: i64 = 0;
    let f = CVodeGetNumSteps(cv, &mut s);
    if f != CV_SUCCESS {
        return Err(format!("CVodeGetNumSteps failed: {f}"));
    }
    let f = CVodeGetNumRhsEvals(cv, &mut r);
    if f != CV_SUCCESS {
        return Err(format!("CVodeGetNumRhsEvals failed: {f}"));
    }
    stats.n_steps += s;
    stats.n_rhs += r;
    Ok(())
}

/// Six-state CVODE segment integrator: BDF + Newton + dense at the spec
/// tolerances (pomega tolerated like the other angles at 1e-10).
pub fn integrate_segment6(
    start: &State2,
    spec: &Segment2,
    observer: &mut dyn FnMut(&Sample2) -> Cmd2,
    samples: &mut Vec<Sample2>,
    stats: &mut Stats2,
) -> Result<Segment2End, String> {
    if spec.t_end <= start.t {
        return Err(format!(
            "integrate_segment6: t_end ({}) must exceed start.t ({})",
            spec.t_end, start.t
        ));
    }
    let mut ctx_out: Option<SUNContext> = None;
    let rc = SUNContext_Create(SUN_COMM_NULL, &mut ctx_out);
    if rc != 0 {
        return Err(format!("SUNContext_Create failed: {rc}"));
    }
    let ctx = ctx_out.ok_or_else(|| "SUNContext_Create returned no context".to_string())?;
    let y = N_VNew_Serial(6, &ctx).ok_or_else(|| "N_VNew_Serial(y) returned None".to_string())?;
    with_data_mut(&y, |d| d.copy_from_slice(&start.y))
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    let abstol =
        N_VNew_Serial(6, &ctx).ok_or_else(|| "N_VNew_Serial(abstol) returned None".to_string())?;
    with_data_mut(&abstol, |d| {
        d.copy_from_slice(&[1.0e-3, 1.0e-6, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-14])
    })
    .ok_or_else(|| "N_VGetArrayPointer returned None for abstol".to_string())?;
    let cv = CVodeCreate(CV_BDF, &ctx)
        .ok_or_else(|| "CVodeCreate(CV_BDF) returned None".to_string())?;
    let mut f = CVodeInit(&cv, rhs6, start.t, &y);
    if f != CV_SUCCESS {
        return Err(format!("CVodeInit failed: {f}"));
    }
    f = CVodeSVtolerances(&cv, params::REL_TOL, &abstol);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSVtolerances failed: {f}"));
    }
    let a_mat =
        SUNDenseMatrix(6, 6, &ctx).ok_or_else(|| "SUNDenseMatrix returned None".to_string())?;
    let ls = SUNLinSol_Dense(&y, &a_mat, &ctx)
        .ok_or_else(|| "SUNLinSol_Dense returned None".to_string())?;
    f = CVodeSetLinearSolver(&cv, &ls, Some(&a_mat));
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetLinearSolver failed: {f}"));
    }
    let user: Box<dyn Any> = Box::new(spec.p.clone());
    f = CVodeSetUserData(&cv, Some(user));
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetUserData failed: {f}"));
    }
    f = CVodeSetMaxNumSteps(&cv, params::MAX_STEPS_PER_CALL);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetMaxNumSteps failed: {f}"));
    }
    f = CVodeSetMaxStep(&cv, params::MAX_STEP);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetMaxStep failed: {f}"));
    }
    f = CVodeSetStopTime(&cv, spec.t_end);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetStopTime failed: {f}"));
    }
    let mut roots_armed = false;
    if spec.p.root_ratio > 0.0 {
        f = CVodeRootInit(&cv, 1, Some(ratio_root6));
        if f != CV_SUCCESS {
            return Err(format!("CVodeRootInit failed: {f}"));
        }
        roots_armed = true;
    }

    let mut cadence = spec.cadence;
    let mut t = start.t;
    let mut root_state: Option<State2> = None;
    let mut tout = t + cadence;
    loop {
        if tout > spec.t_end {
            tout = spec.t_end;
        }
        let cflag = CVode(&cv, tout, &y, &mut t, CV_NORMAL);
        if cflag < 0 {
            return Err(format!("CVode failed with flag {cflag} at t = {t} s"));
        }
        let ys = read_state6(&y)?;
        let smp = Sample2::from_state(t, &ys, spec.stage_tag);
        let at_root = cflag == CV_ROOT_RETURN;
        if spec.record && !at_root {
            samples.push(smp);
        }
        let cmd = observer(&smp);

        let mut done = false;
        if at_root {
            if root_state.is_none() {
                root_state = Some(State2 { t, y: ys });
            }
            if spec.stop_on_root {
                done = true;
            } else if roots_armed {
                let rf = CVodeRootInit(&cv, 0, None);
                if rf != CV_SUCCESS {
                    return Err(format!("CVodeRootInit(disarm) failed: {rf}"));
                }
                roots_armed = false;
            }
        }
        if matches!(cmd, Cmd2::Stop) {
            done = true;
        }
        if cflag == CV_TSTOP_RETURN || t >= spec.t_end - 0.5 {
            done = true;
        }
        if done {
            harvest2(&cv, stats)?;
            let end = Segment2End {
                state: State2 { t, y: ys },
                root_state,
            };
            let mut cv_opt = Some(cv);
            CVodeFree(&mut cv_opt);
            let _ = SUNLinSolFree(Some(ls));
            SUNMatDestroy(a_mat);
            N_VDestroy(y);
            N_VDestroy(abstol);
            let mut ctx_opt = Some(ctx);
            let _ = SUNContext_Free(&mut ctx_opt);
            return Ok(end);
        }
        if let Cmd2::SetCadence(c) = cmd {
            cadence = c;
        }
        if at_root {
            if (tout - t) < 1.0 {
                tout = t + cadence;
            }
            continue;
        }
        if spec.reanchor {
            let two_pi = 2.0 * PI;
            let j = (ys[2] / two_pi).floor();
            if j >= 1.0 {
                harvest2(&cv, stats)?;
                let ynew = [
                    ys[0],
                    ys[1],
                    ys[2] - two_pi * j,
                    ys[3],
                    ys[4] - 3.0 * PI * j,
                    ys[5],
                ];
                with_data_mut(&y, |d| d.copy_from_slice(&ynew))
                    .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
                let rf = CVodeReInit(&cv, t, &y);
                if rf != CV_SUCCESS {
                    return Err(format!("CVodeReInit failed: {rf}"));
                }
                stats.n_reanchor += 1;
            }
        }
        tout = t + cadence;
    }
}

// --- capture detection on gamma2 -------------------------------------------

pub struct Detector2 {
    dense_trigger_ratio: f64,
    dense_cadence: f64,
    passed_ratio: f64,
    window: f64,
    post_capture_dense: f64,
    base_cadence: f64,
    stop_on_decision: bool,
    dense: bool,
    prev_gamma: Option<f64>,
    gu: f64,
    ring: VecDeque<(f64, f64, f64)>,
    sum_ratio: f64,
    scan: usize,
    pub decision: Option<crate::driver::Decision>,
}

impl Detector2 {
    /// The test-1 configuration (path-identity numbers), on gamma2.
    pub fn standard(stop_on_decision: bool) -> Detector2 {
        Detector2 {
            dense_trigger_ratio: 1.53,
            dense_cadence: 2.0 * params::YEAR,
            passed_ratio: 1.45,
            window: 50_000.0 * params::YEAR,
            post_capture_dense: 50_000.0 * params::YEAR,
            base_cadence: 100.0 * params::YEAR,
            stop_on_decision,
            dense: false,
            prev_gamma: None,
            gu: 0.0,
            ring: VecDeque::new(),
            sum_ratio: 0.0,
            scan: 0,
            decision: None,
        }
    }

    pub fn observe(&mut self, s: &Sample2) -> Cmd2 {
        if let Some(dec) = self.decision {
            if let crate::driver::Decision::Captured { t_capture } = dec {
                if self.dense && s.t >= t_capture + self.post_capture_dense {
                    self.dense = false;
                    return Cmd2::SetCadence(self.base_cadence);
                }
            }
            return Cmd2::Continue;
        }
        if !self.dense {
            if s.ratio <= self.dense_trigger_ratio {
                self.dense = true;
                self.prev_gamma = Some(s.gamma2);
                self.gu = 0.0;
                self.ring.clear();
                self.sum_ratio = 0.0;
                return Cmd2::SetCadence(self.dense_cadence);
            }
            return Cmd2::Continue;
        }
        let prev = self.prev_gamma.unwrap_or(s.gamma2);
        self.gu += wrap_pi(s.gamma2 - prev);
        self.prev_gamma = Some(s.gamma2);
        self.ring.push_back((s.t, self.gu, s.ratio));
        self.sum_ratio += s.ratio;
        while let Some(&(t0, _, r0)) = self.ring.front() {
            if t0 < s.t - self.window {
                self.ring.pop_front();
                self.sum_ratio -= r0;
            } else {
                break;
            }
        }
        if s.ratio < self.passed_ratio {
            self.decision = Some(crate::driver::Decision::Passed {
                t: s.t,
                ratio: s.ratio,
            });
            return if self.stop_on_decision { Cmd2::Stop } else { Cmd2::Continue };
        }
        self.scan += 1;
        if self.scan % 200 == 0 && self.ring.len() >= 2 {
            let span = s.t - self.ring.front().map(|&(t0, _, _)| t0).unwrap_or(s.t);
            if span >= 0.999 * self.window {
                let mut gmin = f64::INFINITY;
                let mut gmax = f64::NEG_INFINITY;
                for &(_, g, _) in &self.ring {
                    if g < gmin {
                        gmin = g;
                    }
                    if g > gmax {
                        gmax = g;
                    }
                }
                let mean_ratio = self.sum_ratio / (self.ring.len() as f64);
                if (gmax - gmin) < 4.0 * PI && (1.49..=1.51).contains(&mean_ratio) {
                    let t_capture = self.ring.front().map(|&(t0, _, _)| t0).unwrap_or(s.t);
                    self.decision = Some(crate::driver::Decision::Captured { t_capture });
                    return if self.stop_on_decision { Cmd2::Stop } else { Cmd2::Continue };
                }
            }
        }
        Cmd2::Continue
    }
}

// --- output ----------------------------------------------------------------

pub const SAMPLES2_HEADER: &str = "t_s,a_m,e,M_rad,theta_rad,Omega_rad_s,n_rad_s,ratio,gamma_rad,P_orb_s,P_rot_s,L_spin_kgm2s,L_orb_kgm2s,L_tot_kgm2s,E_spin_j,E_orb_j,stage,pomega_rad";

/// Write samples.csv rows: the test-1 header order with pomega_rad APPENDED
/// (so the database's ALTER TABLE ADD COLUMN pomega_rad lines up).
pub fn write_samples2(dir: &std::path::Path, samples: &[Sample2]) -> Result<usize, String> {
    use std::io::Write as IoWrite;
    let path = dir.join("samples.csv");
    let f = std::fs::File::create(&path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut w = std::io::BufWriter::new(f);
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "{SAMPLES2_HEADER}").map_err(werr)?;
    for s in samples {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            fmt_e(s.t, 12),
            fmt_e(s.a, 12),
            fmt_e(s.e, 12),
            fmt_e(s.m_anom, 12),
            fmt_e(s.theta, 12),
            fmt_e(s.omega, 12),
            fmt_e(s.n, 12),
            fmt_e(s.ratio, 12),
            fmt_e(s.gamma2, 12),
            fmt_e(s.p_orb, 12),
            fmt_e(s.p_rot, 12),
            fmt_e(s.l_spin, 12),
            fmt_e(s.l_orb, 12),
            fmt_e(s.l_tot, 12),
            fmt_e(s.e_spin, 12),
            fmt_e(s.e_orb, 12),
            s.stage,
            fmt_e(s.pomega, 12)
        )
        .map_err(werr)?;
    }
    w.flush().map_err(werr)?;
    Ok(samples.len())
}

pub fn write_restart6(dir: &std::path::Path, state: &State2) -> Result<(), String> {
    use std::io::Write as IoWrite;
    let path = dir.join("restart.csv");
    let f = std::fs::File::create(&path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut w = std::io::BufWriter::new(f);
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "t_s,a_m,e,M_rad,pomega_rad,theta_rad,Omega_rad_s").map_err(werr)?;
    writeln!(
        w,
        "{},{},{},{},{},{},{}",
        fmt_e(state.t, 17),
        fmt_e(state.y[0], 17),
        fmt_e(state.y[1], 17),
        fmt_e(state.y[2], 17),
        fmt_e(state.y[3], 17),
        fmt_e(state.y[4], 17),
        fmt_e(state.y[5], 17)
    )
    .map_err(werr)?;
    w.flush().map_err(werr)?;
    Ok(())
}

pub fn read_restart6(dir: &std::path::Path) -> Result<State2, String> {
    let path = dir.join("restart.csv");
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let line = text
        .lines()
        .nth(1)
        .ok_or_else(|| format!("{} has no data row", path.display()))?;
    let vals: Vec<f64> = line
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<f64>()
                .map_err(|e| format!("bad number {s:?} in {}: {e}", path.display()))
        })
        .collect::<Result<Vec<f64>, String>>()?;
    if vals.len() != 7 {
        return Err(format!(
            "{} data row has {} fields, expected 7",
            path.display(),
            vals.len()
        ));
    }
    Ok(State2 {
        t: vals[0],
        y: [vals[1], vals[2], vals[3], vals[4], vals[5], vals[6]],
    })
}

/// Shared parameter block for the full test-2 physics (movie tides).
pub fn full_params(triaxial_on: bool, root_ratio: f64) -> Rhs2Params {
    let ll = ll_rates();
    Rhs2Params {
        k2tau: params::K2_LOVE * params::TAU_SPEC * params::COMPRESSION_MOVIE,
        triaxial_on,
        gr_on: true,
        jupiter_on: true,
        a11: ll.a11,
        a12: ll.a12,
        root_ratio,
    }
}

/// Test-2 initial state: the spec initial conditions plus pomega = 0.
pub fn initial_state2() -> State2 {
    State2 {
        t: 0.0,
        y: [
            params::A0,
            params::E0,
            params::M0,
            0.0,
            params::THETA0,
            params::OMEGA0,
        ],
    }
}

/// Events reuse test 1's vocabulary.
pub type Event2 = Event;

/// Predicted mean perihelion drift at lock [rad/s] (GR at (a0, e0) plus
/// Jupiter's A11; the A12 term averages to ~0 over the circulating pomega).
pub fn predicted_pomega_dot() -> f64 {
    gr_pomega_dot(params::A0, params::E0) + ll_rates().a11
}

/// Format a rate in arcseconds per century.
pub fn arcsec_cy(rate_rad_s: f64) -> f64 {
    rate_rad_s * RAD_S_TO_ARCSEC_CY
}
