//! The CVODE driver: one generic segment integrator (BDF + Newton + dense,
//! spec tolerances), the angle re-anchoring rule, the capture detector, and
//! shared pieces of the six science runs (A, B, C sweep, B-final, D, E).
//!
//! Everything numerical goes through the vendored pure-Rust SUNDIALS 7.8.0;
//! this file never steps an ODE by hand.

use std::any::Any;
use std::collections::VecDeque;
use std::f64::consts::PI;

use cvode_rs::prelude::*;

use crate::params;
use crate::params::RhsParams;
use crate::rhs;

/// One instant of the five-variable system.
#[derive(Clone, Copy, Debug)]
pub struct State {
    pub t: f64,
    /// y = [a, e, M, theta, Omega]
    pub y: [f64; 5],
}

impl State {
    /// The re-anchoring rule: with j = floor(M / 2pi), subtract 2*pi*j from M
    /// and 3*pi*j from theta. This leaves gamma = 2*theta - 3*M exactly
    /// unchanged (2*(-3pi j) - 3*(-2pi j) = 0) while keeping angles small.
    pub fn reanchored(&self) -> State {
        let two_pi = 2.0 * PI;
        let j = (self.y[2] / two_pi).floor();
        let mut y = self.y;
        y[2] -= two_pi * j;
        y[3] -= 3.0 * PI * j;
        State { t: self.t, y }
    }
}

/// One output sample with every derived quantity the database stores.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub t: f64,
    pub a: f64,
    pub e: f64,
    pub m_anom: f64,
    pub theta: f64,
    pub omega: f64,
    pub n: f64,
    pub ratio: f64,
    pub gamma: f64,
    pub p_orb: f64,
    pub p_rot: f64,
    pub l_spin: f64,
    pub l_orb: f64,
    pub l_tot: f64,
    pub e_spin: f64,
    pub e_orb: f64,
    pub stage: char,
}

/// Reduce an angle to [-pi, pi].
pub fn wrap_pi(x: f64) -> f64 {
    let two_pi = 2.0 * PI;
    x - two_pi * (x / two_pi).round()
}

impl Sample {
    pub fn from_state(t: f64, y: &[f64; 5], stage: char) -> Sample {
        let (a, e, m_anom, theta, omega) = (y[0], y[1], y[2], y[3], y[4]);
        let n = params::mean_motion(a);
        let c = params::moment_of_inertia();
        let l_spin = c * omega;
        let l_orb = params::M_MERCURY * n * (a * a) * (1.0 - e * e).sqrt();
        Sample {
            t,
            a,
            e,
            m_anom,
            theta,
            omega,
            n,
            ratio: omega / n,
            gamma: wrap_pi(2.0 * theta - 3.0 * m_anom),
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

/// A notable moment, written to events.csv.
#[derive(Clone, Debug)]
pub struct Event {
    pub t: f64,
    pub name: String,
    pub value: f64,
}

/// What the per-sample observer may ask the segment loop to do.
pub enum ObserverCmd {
    Continue,
    SetCadence(f64),
    Stop,
}

/// Accumulated CVODE statistics (summed across ReInits and segments).
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub n_steps: i64,
    pub n_rhs: i64,
    pub n_reanchor: i64,
}

/// Configuration of one CVODE integration segment.
pub struct SegmentSpec {
    pub k2tau: f64,
    pub triaxial_on: bool,
    pub t_end: f64,
    pub cadence: f64,
    /// > 0.0 arms a CVODE root function on (Omega - root_ratio * n).
    pub root_ratio: f64,
    /// Stop the segment at the first root (true) or record it and continue,
    /// disabling further root returns (false).
    pub stop_on_root: bool,
    pub reanchor: bool,
    pub stage_tag: char,
    pub record: bool,
}

/// How a segment ended.
pub struct SegmentEnd {
    pub state: State,
    pub root_state: Option<State>,
    pub stopped_by_observer: bool,
}

fn with_data_mut<R>(v: &N_Vector, f: impl FnOnce(&mut [f64]) -> R) -> Option<R> {
    let mut d = N_VGetArrayPointer(v)?;
    Some(f(&mut d))
}

fn read_state(v: &N_Vector) -> Result<[f64; 5], String> {
    let d = N_VGetArrayPointer(v)
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    Ok([d[0], d[1], d[2], d[3], d[4]])
}

fn harvest(cv: &CVodeMem, stats: &mut Stats) -> Result<(), String> {
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

/// Integrate one segment with CVODE (BDF + Newton + dense, spec tolerances),
/// calling `observer` on every output sample. On error paths the solver
/// objects are reclaimed by Rust's ownership (the port is pure Rust — no
/// foreign memory can leak); the success path tears down in the C order.
pub fn integrate_segment(
    start: &State,
    spec: &SegmentSpec,
    observer: &mut dyn FnMut(&Sample) -> ObserverCmd,
    samples: &mut Vec<Sample>,
    stats: &mut Stats,
) -> Result<SegmentEnd, String> {
    if spec.t_end <= start.t {
        return Err(format!(
            "integrate_segment: t_end ({}) must exceed start.t ({})",
            spec.t_end, start.t
        ));
    }
    // --- construction (each None becomes a named error, never a panic) ---
    let mut ctx_out: Option<SUNContext> = None;
    let rc = SUNContext_Create(SUN_COMM_NULL, &mut ctx_out);
    if rc != 0 {
        return Err(format!("SUNContext_Create failed: {rc}"));
    }
    let ctx = ctx_out.ok_or_else(|| "SUNContext_Create returned no context".to_string())?;
    let y = N_VNew_Serial(5, &ctx).ok_or_else(|| "N_VNew_Serial(y) returned None".to_string())?;
    with_data_mut(&y, |d| d.copy_from_slice(&start.y))
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    let abstol =
        N_VNew_Serial(5, &ctx).ok_or_else(|| "N_VNew_Serial(abstol) returned None".to_string())?;
    with_data_mut(&abstol, |d| d.copy_from_slice(&params::ABS_TOL))
        .ok_or_else(|| "N_VGetArrayPointer returned None for abstol".to_string())?;
    let cv = CVodeCreate(CV_BDF, &ctx)
        .ok_or_else(|| "CVodeCreate(CV_BDF) returned None".to_string())?;
    let mut f = CVodeInit(&cv, rhs::rhs, start.t, &y);
    if f != CV_SUCCESS {
        return Err(format!("CVodeInit failed: {f}"));
    }
    f = CVodeSVtolerances(&cv, params::REL_TOL, &abstol);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSVtolerances failed: {f}"));
    }
    let a_mat =
        SUNDenseMatrix(5, 5, &ctx).ok_or_else(|| "SUNDenseMatrix returned None".to_string())?;
    let ls = SUNLinSol_Dense(&y, &a_mat, &ctx)
        .ok_or_else(|| "SUNLinSol_Dense returned None".to_string())?;
    f = CVodeSetLinearSolver(&cv, &ls, Some(&a_mat));
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetLinearSolver failed: {f}"));
    }
    let user: Box<dyn Any> = Box::new(RhsParams {
        k2tau: spec.k2tau,
        triaxial_on: spec.triaxial_on,
        root_ratio: spec.root_ratio,
    });
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
    if spec.root_ratio > 0.0 {
        f = CVodeRootInit(&cv, 1, Some(rhs::ratio_root));
        if f != CV_SUCCESS {
            return Err(format!("CVodeRootInit failed: {f}"));
        }
        roots_armed = true;
    }

    // --- the output loop ---
    let mut cadence = spec.cadence;
    let mut t = start.t;
    let mut root_state: Option<State> = None;
    let mut tout = t + cadence;
    loop {
        if tout > spec.t_end {
            tout = spec.t_end;
        }
        let cflag = CVode(&cv, tout, &y, &mut t, CV_NORMAL);
        if cflag < 0 {
            return Err(format!("CVode failed with flag {cflag} at t = {t} s"));
        }
        let ys = read_state(&y)?;
        let smp = Sample::from_state(t, &ys, spec.stage_tag);
        let at_root = cflag == CV_ROOT_RETURN;
        if spec.record && !at_root {
            samples.push(smp);
        }
        let cmd = observer(&smp);

        let mut done = false;
        let mut stopped_by_observer = false;
        if at_root {
            if root_state.is_none() {
                root_state = Some(State { t, y: ys });
            }
            if spec.stop_on_root {
                done = true;
            } else if roots_armed {
                // First root recorded; disarm so libration-era re-crossings
                // do not keep interrupting the march to each output time.
                // A failed disarm is a named error like every other flag —
                // silently staying armed would corrupt the output grid.
                let rf = CVodeRootInit(&cv, 0, None);
                if rf != CV_SUCCESS {
                    return Err(format!("CVodeRootInit(disarm) failed: {rf}"));
                }
                roots_armed = false;
            }
        }
        if matches!(cmd, ObserverCmd::Stop) {
            stopped_by_observer = true;
            done = true;
        }
        if cflag == CV_TSTOP_RETURN || t >= spec.t_end - 0.5 {
            done = true;
        }
        if done {
            harvest(&cv, stats)?;
            let end = SegmentEnd {
                state: State { t, y: ys },
                root_state,
                stopped_by_observer,
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
        if let ObserverCmd::SetCadence(c) = cmd {
            cadence = c;
        }
        if at_root {
            // Continue toward the pending output time; guard the degenerate
            // case of a root landing essentially on it.
            if (tout - t) < 1.0 {
                tout = t + cadence;
            }
            continue;
        }
        if spec.reanchor {
            let two_pi = 2.0 * PI;
            let j = (ys[2] / two_pi).floor();
            if j >= 1.0 {
                harvest(&cv, stats)?; // CVodeReInit zeroes the counters
                let ynew = [
                    ys[0],
                    ys[1],
                    ys[2] - two_pi * j,
                    ys[3] - 3.0 * PI * j,
                    ys[4],
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

// ---------------------------------------------------------------------------
// Capture detection
// ---------------------------------------------------------------------------

/// The decided fate of a resonance crossing.
#[derive(Clone, Copy, Debug)]
pub enum Decision {
    /// Locked into 3:2: the resonance angle librated (bounded swing) with the
    /// window-mean ratio at 1.5. `t_capture` = start of the qualifying window.
    Captured { t_capture: f64 },
    /// Sailed through: the ratio fell below the passed threshold, heading for
    /// the pseudo-synchronous rate.
    Passed { t: f64, ratio: f64 },
}

/// NOTE on the unwrap margin (review finding, documented not changed): the
/// per-sample gamma step at the dense trigger is 2*n*(trigger - 1.5)*cadence
/// = 3.13 rad for trigger 1.53 at 2-yr cadence — inside but close to the pi
/// aliasing limit (and briefly beyond it under run D's 1.535 trigger). This
/// is provably harmless as configured: aliased early samples age out of the
/// sliding window long before the mean-ratio criterion can be met, and the
/// "passed" decision never uses gamma. Any increase of the trigger, the
/// dense cadence, or the eccentricity must re-derive this margin — and any
/// such change alters the bit-level path, so the sweep must be re-run.
pub struct CaptureConfig {
    pub dense_trigger_ratio: f64,
    pub dense_cadence: f64,
    pub base_cadence: f64,
    pub passed_ratio: f64,
    /// Decision window [s]: gamma must stay bounded (< 4 pi swing) with mean
    /// ratio in [1.49, 1.51] across this span to declare capture.
    pub window: f64,
    /// After capture, keep the dense cadence this long, then revert.
    pub post_capture_dense: f64,
    /// Stop as soon as a decision is made (sweep branches).
    pub stop_on_decision: bool,
    /// If > 0, stop this long after capture (run D); 0 = run to segment end.
    pub stop_after_capture_span: f64,
}

pub struct CaptureDetector {
    cfg: CaptureConfig,
    dense: bool,
    prev_gamma: Option<f64>,
    gu: f64,
    ring: VecDeque<(f64, f64, f64)>, // (t, unwrapped gamma, ratio)
    sum_ratio: f64,
    scan_counter: usize,
    pub decision: Option<Decision>,
}

impl CaptureDetector {
    pub fn new(cfg: CaptureConfig) -> CaptureDetector {
        CaptureDetector {
            cfg,
            dense: false,
            prev_gamma: None,
            gu: 0.0,
            ring: VecDeque::new(),
            sum_ratio: 0.0,
            scan_counter: 0,
            decision: None,
        }
    }

    pub fn observe(&mut self, s: &Sample) -> ObserverCmd {
        if let Some(dec) = self.decision {
            match dec {
                Decision::Captured { t_capture } => {
                    if self.cfg.stop_after_capture_span > 0.0
                        && s.t >= t_capture + self.cfg.stop_after_capture_span
                    {
                        return ObserverCmd::Stop;
                    }
                    if self.dense && s.t >= t_capture + self.cfg.post_capture_dense {
                        self.dense = false;
                        return ObserverCmd::SetCadence(self.cfg.base_cadence);
                    }
                }
                Decision::Passed { .. } => {}
            }
            return ObserverCmd::Continue;
        }
        if !self.dense {
            if s.ratio <= self.cfg.dense_trigger_ratio {
                self.dense = true;
                self.prev_gamma = Some(s.gamma);
                self.gu = 0.0;
                self.ring.clear();
                self.sum_ratio = 0.0;
                return ObserverCmd::SetCadence(self.cfg.dense_cadence);
            }
            return ObserverCmd::Continue;
        }
        // Dense and undecided: unwrap gamma, maintain the sliding window.
        let prev = self.prev_gamma.unwrap_or(s.gamma);
        self.gu += wrap_pi(s.gamma - prev);
        self.prev_gamma = Some(s.gamma);
        self.ring.push_back((s.t, self.gu, s.ratio));
        self.sum_ratio += s.ratio;
        while let Some(&(t0, _, r0)) = self.ring.front() {
            if t0 < s.t - self.cfg.window {
                self.ring.pop_front();
                self.sum_ratio -= r0;
            } else {
                break;
            }
        }
        if s.ratio < self.cfg.passed_ratio {
            self.decision = Some(Decision::Passed {
                t: s.t,
                ratio: s.ratio,
            });
            return if self.cfg.stop_on_decision {
                ObserverCmd::Stop
            } else {
                ObserverCmd::Continue
            };
        }
        self.scan_counter += 1;
        if self.scan_counter % 200 == 0 && self.ring.len() >= 2 {
            let span = s.t - self.ring.front().map(|&(t0, _, _)| t0).unwrap_or(s.t);
            if span >= 0.999 * self.cfg.window {
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
                    self.decision = Some(Decision::Captured { t_capture });
                    return if self.cfg.stop_on_decision {
                        ObserverCmd::Stop
                    } else {
                        ObserverCmd::Continue
                    };
                }
            }
        }
        ObserverCmd::Continue
    }
}

// ---------------------------------------------------------------------------
// Shared run pieces
// ---------------------------------------------------------------------------

/// The spec-literal k2*tau product [s].
pub fn k2tau_spec() -> f64 {
    params::K2_LOVE * params::TAU_SPEC
}

/// The movie (S = 1000) k2*tau product [s].
pub fn k2tau_movie() -> f64 {
    params::K2_LOVE * params::TAU_SPEC * params::COMPRESSION_MOVIE
}

/// The spec initial state.
pub fn initial_state() -> State {
    State {
        t: 0.0,
        y: [
            params::A0,
            params::E0,
            params::M0,
            params::THETA0,
            params::OMEGA0,
        ],
    }
}

/// One sweep branch's outcome (a branches.csv row).
#[derive(Clone, Copy, Debug)]
pub struct BranchOutcome {
    pub branch_id: usize,
    pub theta_offset: f64,
    pub captured: bool,
    pub t_outcome: f64,
    pub final_ratio: f64,
}
