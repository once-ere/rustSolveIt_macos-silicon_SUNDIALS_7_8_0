#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

//! mercury_rs command-line entry: the six science runs of the planet_Mercury
//! project. Every run ends with explicit PASS/FAIL self-check lines and a
//! final SUCCESS or FAILURE verdict (nonzero exit on failure). All numeric
//! output goes through the engine's C-style formatters for byte-identical
//! re-runs.

use std::f64::consts::PI;

use sundials_core::sundials_utils::{fmt_e, fmt_f};

use mercury_rs::driver::{
    initial_state, integrate_segment, k2tau_movie, k2tau_spec, BranchOutcome, CaptureConfig,
    CaptureDetector, Decision, Event, ObserverCmd, Sample, SegmentSpec, State, Stats,
};
use mercury_rs::output;
use mercury_rs::params;

const OBSERVED_P_ORB_D: f64 = 87.969;
const OBSERVED_P_ROT_D: f64 = 58.646;
const OBSERVED_P_SOLAR_D: f64 = 175.938;
const DAY: f64 = 86400.0;
const L_BUDGET: f64 = 1.0e-9;

fn myr(t: f64) -> f64 {
    t / (1.0e6 * params::YEAR)
}

fn check(ok: bool, name: &str, detail: &str) -> bool {
    if ok {
        println!("PASS - {name}: {detail}");
    } else {
        println!("FAIL - {name}: {detail}");
    }
    ok
}

fn print_config() {
    println!("mercury_rs configuration (all SI units)");
    println!("  G                = {} m^3 kg^-1 s^-2", fmt_e(params::G, 5));
    println!("  M_sun            = {} kg", fmt_e(params::M_SUN, 5));
    println!("  m_mercury        = {} kg", fmt_e(params::M_MERCURY, 5));
    println!("  R_mercury        = {} m", fmt_e(params::R_MERCURY, 5));
    println!("  C_factor         = {}", fmt_e(params::C_FACTOR, 5));
    println!("  C (moment)       = {} kg m^2", fmt_e(params::moment_of_inertia(), 5));
    println!("  (B-A)/C          = {}", fmt_e(params::B_MINUS_A_OVER_C, 5));
    println!("  k2               = {}", fmt_e(params::K2_LOVE, 5));
    println!("  tau_lag (spec)   = {} s", fmt_e(params::TAU_SPEC, 5));
    println!("  compression S    = {}", fmt_e(params::COMPRESSION_MOVIE, 5));
    println!("  tau_lag (movie)  = {} s", fmt_e(params::TAU_SPEC * params::COMPRESSION_MOVIE, 5));
    println!("  a0               = {} m", fmt_e(params::A0, 6));
    println!("  e0               = {}", fmt_e(params::E0, 5));
    println!("  M0               = {} rad", fmt_e(params::M0, 5));
    println!("  theta0           = {} rad", fmt_e(params::THETA0, 5));
    println!("  Omega0           = {} rad/s", fmt_e(params::OMEGA0, 5));
    println!("  n(a0)            = {} rad/s", fmt_e(params::mean_motion(params::A0), 6));
    println!("  Omega0/n(a0)     = {}", fmt_f(params::OMEGA0 / params::mean_motion(params::A0), 4));
    println!("  rel_tol          = {}", fmt_e(params::REL_TOL, 2));
    println!(
        "  abs_tol[a,e,M,theta,Omega] = [{}, {}, {}, {}, {}]",
        fmt_e(params::ABS_TOL[0], 2),
        fmt_e(params::ABS_TOL[1], 2),
        fmt_e(params::ABS_TOL[2], 2),
        fmt_e(params::ABS_TOL[3], 2),
        fmt_e(params::ABS_TOL[4], 2)
    );
    println!("  max_step         = {} s (10 days)", fmt_e(params::MAX_STEP, 5));
    println!("  t_final          = {} s (10 Myr)", fmt_e(params::T_FINAL, 5));
    println!("  solver           = CVODE BDF + NEWTON + DENSE (sundials_rs 7.8.0, pure Rust)");
    println!("  stage handover   = Omega/n <= {}", fmt_f(params::STAGE_HANDOVER_RATIO, 2));
    println!("  restart save     = Omega/n = {}", fmt_f(params::RESTART_RATIO, 2));
    println!("SUCCESS");
}

fn base_manifest(run_id: &str, description: &str, k2tau: f64, ic: &State) -> output::Manifest {
    output::Manifest {
        run_id: run_id.to_string(),
        description: description.to_string(),
        k2: params::K2_LOVE,
        tau_lag_s: k2tau / params::K2_LOVE,
        compression: (k2tau / params::K2_LOVE) / params::TAU_SPEC,
        a0_m: ic.y[0],
        e0: ic.y[1],
        M0_rad: ic.y[2],
        theta0_rad: ic.y[3],
        Omega0_rad_s: ic.y[4],
        t_final_s: 0.0,
        n_steps: 0,
        n_rhs_evals: 0,
        n_reanchor: 0,
        verdict: "FAILURE".to_string(),
        extras: Vec::new(),
    }
}

fn finish_manifest(
    dir: &std::path::Path,
    mut m: output::Manifest,
    stats: &Stats,
    t_end: f64,
    ok: bool,
) -> Result<bool, String> {
    m.t_final_s = t_end;
    m.n_steps = stats.n_steps;
    m.n_rhs_evals = stats.n_rhs;
    m.n_reanchor = stats.n_reanchor;
    m.verdict = if ok { "SUCCESS" } else { "FAILURE" }.to_string();
    output::write_manifest(dir, &m)?;
    Ok(ok)
}

// ---------------------------------------------------------------------------
// Run A — spec-literal honesty run
// ---------------------------------------------------------------------------

fn run_a() -> Result<bool, String> {
    println!("run A (A_spec_literal): tau = 100 s, 10 Myr — the honest spec-literal result");
    let dir = output::fresh_run_dir("A_spec_literal")?;
    let ic = initial_state();
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let mut count = 0usize;
    let mut obs = |s: &Sample| {
        count += 1;
        if count % 2000 == 0 {
            println!("  t = {} Myr  Omega/n = {}", fmt_f(myr(s.t), 3), fmt_f(s.ratio, 4));
        }
        ObserverCmd::Continue
    };
    let spec1 = SegmentSpec {
        k2tau: k2tau_spec(),
        triaxial_on: false,
        t_end: params::T_FINAL,
        cadence: 1000.0 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let end1 = integrate_segment(&ic, &spec1, &mut obs, &mut samples, &mut stats)?;
    let final_ratio = samples.last().map(|s| s.ratio).unwrap_or(0.0);

    // Short full-model proof segment (triaxial torque on, spec tau, 50 years).
    let n_before = samples.len();
    let spec2 = SegmentSpec {
        k2tau: k2tau_spec(),
        triaxial_on: true,
        t_end: 50.0 * params::YEAR,
        cadence: 0.1 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'R',
        record: true,
    };
    let mut noop = |_s: &Sample| ObserverCmd::Continue;
    let _end2 = integrate_segment(&ic, &spec2, &mut noop, &mut samples, &mut stats)?;
    let proof_rows = samples.len() - n_before;
    let proof_ratio = samples.last().map(|s| s.ratio).unwrap_or(0.0);

    let mut ok = true;
    ok &= check(
        (170.0..=181.5).contains(&final_ratio),
        "A.final_ratio_barely_moved",
        &format!(
            "final Omega/n = {} after 10 Myr (Finding F1: spec-strength tides need ~4.7 Gyr)",
            fmt_f(final_ratio, 4)
        ),
    );
    ok &= check(
        final_ratio > 100.0,
        "A.no_capture",
        "spin nowhere near any resonance — no capture event",
    );
    ok &= check(
        proof_rows > 100 && (proof_ratio - params::OMEGA0 / params::mean_motion(params::A0)).abs() < 0.01,
        "A.full_model_proof_segment",
        &format!(
            "{proof_rows} full-model rows over 50 yr, Omega/n = {}",
            fmt_f(proof_ratio, 6)
        ),
    );
    let rows = output::write_samples(&dir, &samples)?;
    output::write_events(&dir, &[])?;
    println!("  wrote {rows} samples, 0 events");
    let mut m = base_manifest(
        "A_spec_literal",
        "Spec-literal run: tau = 100 s over 10 Myr (secular stage) + 50 yr full-model proof segment",
        k2tau_spec(),
        &ic,
    );
    m.extras.push(("cadence_secular_s".into(), fmt_e(1000.0 * params::YEAR, 6)));
    m.extras.push(("cadence_proof_s".into(), fmt_e(0.1 * params::YEAR, 6)));
    finish_manifest(&dir, m, &stats, end1.state.t, ok)
}

// ---------------------------------------------------------------------------
// Run B — the movie braking run, to the restart save at Omega/n = 1.6
// ---------------------------------------------------------------------------

fn run_b() -> Result<bool, String> {
    println!("run B (B_movie): S = 1000 movie braking to the restart save at Omega/n = 1.6");
    let dir = output::fresh_run_dir("B_movie")?;
    let ic = initial_state();
    let mut samples: Vec<Sample> = Vec::new();
    let mut events: Vec<Event> = Vec::new();
    let mut stats = Stats::default();

    // Stage S: secular braking until the handover root at Omega/n = 2.2.
    let mut count = 0usize;
    let mut obs1 = |s: &Sample| {
        count += 1;
        if count % 500 == 0 {
            println!("  [S] t = {} Myr  Omega/n = {}", fmt_f(myr(s.t), 3), fmt_f(s.ratio, 4));
        }
        ObserverCmd::Continue
    };
    let spec1 = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: false,
        t_end: 8.0e6 * params::YEAR,
        cadence: 1000.0 * params::YEAR,
        root_ratio: params::STAGE_HANDOVER_RATIO,
        stop_on_root: true,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let end1 = integrate_segment(&ic, &spec1, &mut obs1, &mut samples, &mut stats)?;
    let mut seen_52 = false;
    for s in &samples {
        if !seen_52 && s.ratio < 2.5 {
            seen_52 = true;
            events.push(Event { t: s.t, name: "cross_5:2".into(), value: s.ratio });
        }
    }
    let handover = end1
        .root_state
        .ok_or_else(|| "run B: the Omega/n = 2.2 handover root never fired within 8 Myr".to_string())?;
    events.push(Event {
        t: handover.t,
        name: "stage_handover".into(),
        value: handover.y[4] / params::mean_motion(handover.y[0]),
    });
    println!("  stage handover at t = {} Myr", fmt_f(myr(handover.t), 4));

    // Stage R: full model until the restart root at Omega/n = 1.6.
    let start2 = handover.reanchored();
    let mut seen_21_at: Option<(f64, f64)> = None;
    let mut count2 = 0usize;
    let mut obs2 = |s: &Sample| {
        count2 += 1;
        if seen_21_at.is_none() && s.ratio < 2.0 {
            seen_21_at = Some((s.t, s.ratio));
        }
        if count2 % 2000 == 0 {
            println!("  [R] t = {} Myr  Omega/n = {}", fmt_f(myr(s.t), 3), fmt_f(s.ratio, 4));
        }
        ObserverCmd::Continue
    };
    let spec2 = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: true,
        t_end: handover.t + 3.0e6 * params::YEAR,
        cadence: 100.0 * params::YEAR,
        root_ratio: params::RESTART_RATIO,
        stop_on_root: true,
        reanchor: true,
        stage_tag: 'R',
        record: true,
    };
    let end2 = integrate_segment(&start2, &spec2, &mut obs2, &mut samples, &mut stats)?;
    if let Some((t, v)) = seen_21_at {
        events.push(Event { t, name: "cross_2:1".into(), value: v });
    }
    let restart_raw = match end2.root_state {
        Some(st) => st,
        None => {
            let r = end2.state.y[4] / params::mean_motion(end2.state.y[0]);
            if (r - 2.0).abs() < 0.05 {
                return Err(
                    "run B: spin appears CAPTURED INTO THE 2:1 RESONANCE (a ~1% chance per the plan) — \
                     documented contingency: record this in the provenance and rerun with a slightly \
                     different theta0"
                        .to_string(),
                );
            }
            return Err(format!(
                "run B: the Omega/n = 1.6 restart root never fired (final ratio {})",
                fmt_f(r, 4)
            ));
        }
    };
    let restart = restart_raw.reanchored();
    output::write_restart(&dir, &restart)?;
    let restart_ratio = restart.y[4] / params::mean_motion(restart.y[0]);
    events.push(Event { t: restart.t, name: "restart_saved".into(), value: restart_ratio });

    let mut ok = true;
    ok &= check(
        (restart_ratio - params::RESTART_RATIO).abs() < 1.0e-6,
        "B.restart_ratio_exact",
        &format!("root-found Omega/n = {} (target 1.6)", fmt_e(restart_ratio, 12)),
    );
    ok &= check(
        (2.5e6..=6.5e6).contains(&(restart.t / params::YEAR)),
        "B.restart_time_plausible",
        &format!("restart at t = {} Myr (predicted ~4.5 Myr +/- 20%)", fmt_f(myr(restart.t), 4)),
    );
    ok &= check(
        seen_52 && seen_21_at.is_some(),
        "B.crossings_recorded",
        "5:2 and 2:1 crossings in the event log",
    );
    ok &= check(!samples.is_empty(), "B.samples_recorded", &format!("{} rows", samples.len()));
    let rows = output::write_samples(&dir, &samples)?;
    output::write_events(&dir, &events)?;
    println!("  wrote {rows} samples, {} events", events.len());
    let mut m = base_manifest(
        "B_movie",
        "Movie braking run (S = 1000): secular stage to Omega/n = 2.2, full model to the restart save at 1.6",
        k2tau_movie(),
        &ic,
    );
    m.extras.push(("cadence_secular_s".into(), fmt_e(1000.0 * params::YEAR, 6)));
    m.extras.push(("cadence_resonant_s".into(), fmt_e(100.0 * params::YEAR, 6)));
    m.extras.push(("handover_ratio".into(), fmt_e(params::STAGE_HANDOVER_RATIO, 6)));
    m.extras.push(("restart_ratio".into(), fmt_e(params::RESTART_RATIO, 6)));
    finish_manifest(&dir, m, &stats, restart.t, ok)
}

// ---------------------------------------------------------------------------
// Run C — the phase sweep
// ---------------------------------------------------------------------------

/// The one capture-detection configuration shared by every sweep branch AND
/// by run B-final. Capture outcomes are bit-level path-sensitive: the
/// detector's cadence switches (and the armed 1.5-root) change CVODE's step
/// sequence, so the sweep's promise only transfers to the continuation if
/// BOTH use the numerically identical path up to the decision moment.
/// (`stop_on_decision` is the only permitted difference — it acts AT the
/// decision, after the coin has landed.)
fn branch_config(stop_on_decision: bool) -> CaptureConfig {
    CaptureConfig {
        dense_trigger_ratio: 1.53,
        dense_cadence: 2.0 * params::YEAR,
        base_cadence: 100.0 * params::YEAR,
        passed_ratio: 1.45,
        window: 50_000.0 * params::YEAR,
        post_capture_dense: 50_000.0 * params::YEAR,
        stop_on_decision,
        stop_after_capture_span: 0.0,
    }
}

fn run_one_branch(
    restart: State,
    branch_id: usize,
    n_branches: usize,
) -> Result<(BranchOutcome, Stats), String> {
    let offset = (branch_id as f64) * (PI / (n_branches as f64));
    let mut st = restart;
    st.y[3] += offset;
    let mut det = CaptureDetector::new(branch_config(true));
    let mut stats = Stats::default();
    let mut samples: Vec<Sample> = Vec::new();
    let spec = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: true,
        t_end: restart.t + 1.5e6 * params::YEAR,
        cadence: 100.0 * params::YEAR,
        // The 1.5-root is armed here for path identity with run B-final —
        // an armed root changes CVODE's step sequence, and the branch's
        // capture coin flip must reproduce bit-for-bit in the continuation.
        root_ratio: 1.5,
        stop_on_root: false,
        reanchor: true,
        stage_tag: 'R',
        record: false,
    };
    let mut obs = |s: &Sample| det.observe(s);
    let end = integrate_segment(&st, &spec, &mut obs, &mut samples, &mut stats)?;
    let final_ratio = end.state.y[4] / params::mean_motion(end.state.y[0]);
    match det.decision {
        Some(Decision::Captured { t_capture }) => Ok((
            BranchOutcome {
                branch_id,
                theta_offset: offset,
                captured: true,
                t_outcome: t_capture,
                final_ratio,
            },
            stats,
        )),
        Some(Decision::Passed { t, ratio }) => Ok((
            BranchOutcome {
                branch_id,
                theta_offset: offset,
                captured: false,
                t_outcome: t,
                final_ratio: ratio,
            },
            stats,
        )),
        None => Err(format!(
            "sweep branch {branch_id}: undecided after 1.5 Myr (final ratio {})",
            fmt_f(final_ratio, 4)
        )),
    }
}

fn run_sweep(n_branches: usize) -> Result<bool, String> {
    println!("run C (C_sweep): {n_branches}-branch spin-phase sweep from the saved restart");
    let b_dir = output::run_dir("B_movie")?;
    let restart = output::read_restart(&b_dir)?;
    let dir = output::fresh_run_dir("C_sweep")?;

    // Worker count is machine-dependent and deliberately NOT printed (it
    // would make transcripts differ between hosts); branch results and the
    // integer stat sums are provably order- and worker-count-independent.
    let workers = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(8);
    println!("  restart: t = {} Myr", fmt_f(myr(restart.t), 4));
    let (tx, rx) = std::sync::mpsc::channel::<(usize, Result<(BranchOutcome, Stats), String>)>();
    let mut spawned = 0usize;
    for w in 0..workers {
        let tx = tx.clone();
        let ids: Vec<usize> = (0..n_branches).filter(|k| k % workers == w).collect();
        spawned += ids.len();
        std::thread::spawn(move || {
            for k in ids {
                let res = run_one_branch(restart, k, n_branches);
                if tx.send((k, res)).is_err() {
                    return;
                }
            }
        });
    }
    drop(tx);
    // Collect EVERYTHING first, then handle errors in branch order, so even
    // a failing sweep produces a deterministic transcript.
    let mut collected: Vec<(usize, Result<(BranchOutcome, Stats), String>)> =
        rx.iter().take(spawned).collect();
    collected.sort_by_key(|(k, _)| *k);
    let mut outcomes: Vec<BranchOutcome> = Vec::with_capacity(n_branches);
    let mut stats = Stats::default();
    for (k, res) in collected {
        let (b, s) = res.map_err(|e| format!("branch {k}: {e}"))?;
        outcomes.push(b);
        stats.n_steps += s.n_steps;
        stats.n_rhs += s.n_rhs;
        stats.n_reanchor += s.n_reanchor;
    }
    outcomes.sort_by_key(|b| b.branch_id);
    let captured: Vec<usize> = outcomes.iter().filter(|b| b.captured).map(|b| b.branch_id).collect();
    let canonical = captured.first().copied();
    for b in &outcomes {
        println!(
            "  branch {:2}  offset = {} rad  {}  t = {} Myr  final ratio = {}",
            b.branch_id,
            fmt_e(b.theta_offset, 6),
            if b.captured { "CAPTURED" } else { "passed  " },
            fmt_f(myr(b.t_outcome), 4),
            fmt_f(b.final_ratio, 4)
        );
    }
    let frac = (captured.len() as f64) / (outcomes.len().max(1) as f64);
    println!(
        "  capture fraction: {}/{} = {} (Goldreich-Peale's simple estimate is ~0.070; the velocity-dependent tidal term makes true constant-time-lag odds somewhat higher)",
        captured.len(),
        outcomes.len(),
        fmt_f(frac, 4)
    );
    let mut ok = true;
    ok &= check(
        outcomes.len() == n_branches,
        "C.all_branches_decided",
        &format!("{} of {n_branches}", outcomes.len()),
    );
    let any = !captured.is_empty();
    ok &= check(
        any,
        "C.at_least_one_capture",
        if any {
            "canonical branch selected"
        } else {
            "ZERO captures (~1% chance) — run the documented finer-grid contingency re-sweep"
        },
    );
    output::write_branches(&dir, &outcomes, canonical)?;
    output::write_events(&dir, &[])?;
    output::write_samples(&dir, &[])?;
    let mut m = base_manifest(
        "C_sweep",
        "64-branch spin-phase sweep at the 3:2 crossing (offsets k*pi/64 on theta at the restart state)",
        k2tau_movie(),
        &restart,
    );
    m.extras.push(("branches".into(), format!("{n_branches}")));
    m.extras.push(("offset_step_rad".into(), fmt_e(PI / (n_branches as f64), 17)));
    m.extras.push(("captured_count".into(), format!("{}", captured.len())));
    m.extras.push((
        "canonical_branch".into(),
        canonical.map(|k| format!("{k}")).unwrap_or_else(|| "-1".to_string()),
    ));
    m.extras.push(("n_steps_note".into(), "\"summed over all branch integrations\"".into()));
    let ok = finish_manifest(&dir, m, &stats, restart.t + 1.5e6 * params::YEAR, ok)?;
    if !any {
        // Special exit signal for the orchestrator's contingency path.
        println!("FAILURE");
        std::process::exit(3);
    }
    Ok(ok)
}

// ---------------------------------------------------------------------------
// Run B-final — the canonical captured branch, continued to 10 Myr
// ---------------------------------------------------------------------------

fn run_b_final(branch: usize, n_branches: usize) -> Result<bool, String> {
    if branch >= n_branches {
        return Err(format!(
            "run-b-final: --branch {branch} was never swept (valid: 0..{})",
            n_branches - 1
        ));
    }
    println!("run B-final (B_final): branch {branch} continued from the restart to 10 Myr");
    let b_dir = output::run_dir("B_movie")?;
    let restart = output::read_restart(&b_dir)?;
    let dir = output::fresh_run_dir("B_final")?;
    let offset = (branch as f64) * (PI / (n_branches as f64));
    let mut st = restart;
    st.y[3] += offset;

    // Identical detector configuration to the sweep branches (path identity;
    // see branch_config) — only stop_on_decision differs, which acts after
    // the capture coin has already landed.
    let mut det = CaptureDetector::new(branch_config(false));
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let mut count = 0usize;
    let mut obs = |s: &Sample| {
        count += 1;
        if count % 5000 == 0 {
            println!("  t = {} Myr  Omega/n = {}", fmt_f(myr(s.t), 3), fmt_f(s.ratio, 6));
        }
        det.observe(s)
    };
    let spec = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: true,
        t_end: params::T_FINAL,
        cadence: 100.0 * params::YEAR,
        root_ratio: 1.5,
        stop_on_root: false,
        reanchor: true,
        stage_tag: 'R',
        record: true,
    };
    let end = integrate_segment(&st, &spec, &mut obs, &mut samples, &mut stats)?;

    let mut events: Vec<Event> = Vec::new();
    if let Some(cross) = end.root_state {
        events.push(Event {
            t: cross.t,
            name: "cross_3:2".into(),
            value: cross.y[4] / params::mean_motion(cross.y[0]),
        });
    }
    let decision = det.decision;
    let t_capture = match decision {
        Some(Decision::Captured { t_capture }) => {
            events.push(Event { t: t_capture, name: "capture_detected".into(), value: 1.5 });
            Some(t_capture)
        }
        _ => None,
    };
    events.push(Event { t: end.state.t, name: "reanchor".into(), value: stats.n_reanchor as f64 });

    let mut ok = true;
    ok &= check(
        t_capture.is_some(),
        "Bf.capture_detected",
        &t_capture
            .map(|t| format!("captured at t = {} Myr", fmt_f(myr(t), 4)))
            .unwrap_or_else(|| "branch did not capture — reselect the canonical branch".to_string()),
    );
    if let Some(tc) = t_capture {
        // Locked-mean check over the well-settled era. A locked planet still
        // LIBRATES (rocks) forever with tiny amplitude, so instantaneous
        // samples wiggle by ~1e-4 in the ratio; the lock statement — spin
        // ratio exactly 3/2 — is about the TIME AVERAGE (deviation DEV-4).
        let settled: Vec<&Sample> =
            samples.iter().filter(|s| s.t > tc + 1.0e6 * params::YEAR).collect();
        let mut mean_ratio = 0.0;
        if settled.is_empty() {
            ok = check(false, "Bf.locked_mean_ratio", "no settled-era samples");
        } else {
            mean_ratio =
                settled.iter().map(|s| s.ratio).sum::<f64>() / (settled.len() as f64);
            ok &= check(
                (mean_ratio - 1.5).abs() <= 1.0e-4,
                "Bf.locked_mean_ratio",
                &format!("settled-era mean Omega/n = {}", fmt_e(mean_ratio, 10)),
            );
        }
        // Final periods vs the observed Mercury.
        let last = samples.last().ok_or_else(|| "run B-final produced no samples".to_string())?;
        let p_orb_d = last.p_orb / DAY;
        let p_rot_d = last.p_rot / DAY;
        let p_solar_d = 1.0 / (1.0 / p_rot_d - 1.0 / p_orb_d).abs();
        ok &= check(
            (p_orb_d - OBSERVED_P_ORB_D).abs() / OBSERVED_P_ORB_D < 1.0e-3,
            "Bf.P_orb",
            &format!("final P_orb = {} d (observed 87.969)", fmt_f(p_orb_d, 4)),
        );
        ok &= check(
            (p_rot_d - OBSERVED_P_ROT_D).abs() / OBSERVED_P_ROT_D < 1.0e-3,
            "Bf.P_rot",
            &format!("final P_rot = {} d (observed 58.646)", fmt_f(p_rot_d, 4)),
        );
        ok &= check(
            (p_solar_d - OBSERVED_P_SOLAR_D).abs() / OBSERVED_P_SOLAR_D < 1.0e-3,
            "Bf.P_solar",
            &format!(
                "final solar day = {} d (observed 175.938, = 2 Mercury years)",
                fmt_f(p_solar_d, 4)
            ),
        );
        // The 2/3 lock to 5 significant figures is a statement about the
        // time-averaged spin (instantaneous samples carry the residual
        // libration wiggle of ~2e-4; deviation DEV-4). Measured on its OWN
        // statistic — the settled-era mean of P_rot/P_orb — distinct from
        // the mean-of-Omega/n check above (review hardening).
        let settled2: Vec<&Sample> =
            samples.iter().filter(|s| s.t > tc + 1.0e6 * params::YEAR).collect();
        let mean_pr = if settled2.is_empty() {
            0.0
        } else {
            settled2.iter().map(|s| s.p_rot / s.p_orb).sum::<f64>() / (settled2.len() as f64)
        };
        let ratio_23 = (mean_pr - 2.0 / 3.0).abs() / (2.0 / 3.0);
        ok &= check(
            ratio_23 < 1.0e-5,
            "Bf.two_thirds_lock",
            &format!(
                "settled-era mean P_rot/P_orb - 2/3 relative = {} (instantaneous final sample wiggles at {} from residual libration)",
                fmt_e(ratio_23, 4),
                fmt_e((last.p_rot / last.p_orb - 2.0 / 3.0).abs() / (2.0 / 3.0), 4)
            ),
        );
        let _ = mean_ratio; // retained by the locked-mean check above
        // Angular-momentum ledger — two eras with different physics (review
        // finding, DEV-7). BEFORE capture the handle torque averages to ~zero
        // over the circulating resonance angle, so the model conserves
        // L_tot = C*Omega + L_orb exactly (secular part): budget 1e-9.
        // AFTER capture the lock is held by a NONZERO mean handle torque
        // <T_tri> = -<T_tidal> = K n (1.5 f1 - f2), and because the source
        // spec deliberately gives the handle torque no orbital back-reaction,
        // the model itself leaks L_tot at that secular rate; the locked era
        // is therefore checked against the PREDICTED leak bound, not against
        // zero. (Numerically the recorded leak is even smaller: at lock the
        // per-step change in a falls below f64 resolution of a ~ 5.8e10 m,
        // freezing the recorded orbit — also documented in the provenance.)
        let l0 = samples.first().map(|s| s.l_tot).unwrap_or(1.0);
        let pre_drift = samples
            .iter()
            .filter(|s| s.t <= tc)
            .map(|s| (s.l_tot - l0).abs() / l0)
            .fold(0.0f64, f64::max);
        ok &= check(
            pre_drift <= L_BUDGET,
            "Bf.ledger_precapture",
            &format!(
                "max |L_tot - L_tot(start)|/L_tot = {} before capture (budget {})",
                fmt_e(pre_drift, 4),
                fmt_e(L_BUDGET, 2)
            ),
        );
        let l_cap = samples
            .iter()
            .find(|s| s.t >= tc)
            .map(|s| s.l_tot)
            .unwrap_or(l0);
        let locked_drift = samples
            .iter()
            .filter(|s| s.t > tc)
            .map(|s| (s.l_tot - l_cap).abs() / l0)
            .fold(0.0f64, f64::max);
        let e_cap = samples.iter().find(|s| s.t >= tc).map(|s| s.e).unwrap_or(params::E0);
        let k_movie = params::tidal_k(params::A0, k2tau_movie());
        let n0 = params::mean_motion(params::A0);
        let t_tri_mean = k_movie * n0 * (1.5 * mercury_rs::hut::f1(e_cap) - mercury_rs::hut::f2(e_cap));
        let predicted_leak = t_tri_mean * (end.state.t - tc) / l0;
        ok &= check(
            locked_drift <= 1.5 * predicted_leak + 2.0e-10,
            "Bf.ledger_locked_era",
            &format!(
                "locked-era |dL_tot|/L_tot = {} vs the model's own predicted secular leak {} (the spec's handle torque has no orbital back-reaction; recorded drift is further suppressed by f64 orbit quantization)",
                fmt_e(locked_drift, 4),
                fmt_e(predicted_leak, 4)
            ),
        );
        // Libration amplitude decays: compare gamma swing in early vs late
        // 100-kyr bins after capture.
        let bin = 100_000.0 * params::YEAR;
        let swing = |t0: f64| -> f64 {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for s in samples.iter().filter(|s| s.t >= t0 && s.t < t0 + bin) {
                if s.gamma < lo {
                    lo = s.gamma;
                }
                if s.gamma > hi {
                    hi = s.gamma;
                }
            }
            if hi >= lo {
                hi - lo
            } else {
                0.0
            }
        };
        let early = swing(tc + bin);
        let late = swing(end.state.t - 2.0 * bin);
        ok &= check(
            late < early && early > 0.0,
            "Bf.libration_decays",
            &format!(
                "gamma swing {} rad (early bin) -> {} rad (late bin)",
                fmt_f(early, 4),
                fmt_f(late, 4)
            ),
        );
        println!(
            "  capture at t = {} Myr of movie time; at the real (uncompressed) tidal strength the same story takes ~1000x longer",
            fmt_f(myr(tc), 4)
        );
    } else {
        ok = false;
    }

    let rows = output::write_samples(&dir, &samples)?;
    output::write_events(&dir, &events)?;
    println!("  wrote {rows} samples, {} events", events.len());
    let mut m = base_manifest(
        "B_final",
        "Canonical captured branch continued from the restart to 10 Myr: capture, libration, locked present day",
        k2tau_movie(),
        &st,
    );
    m.extras.push(("branch".into(), format!("{branch}")));
    m.extras.push(("theta_offset_rad".into(), fmt_e(offset, 17)));
    m.extras.push(("restart_t_s".into(), fmt_e(restart.t, 17)));
    m.extras.push(("cadence_base_s".into(), fmt_e(100.0 * params::YEAR, 6)));
    m.extras.push(("cadence_dense_s".into(), fmt_e(2.0 * params::YEAR, 6)));
    finish_manifest(&dir, m, &stats, end.state.t, ok)
}

// ---------------------------------------------------------------------------
// Run D — guaranteed capture at e = 0.285
// ---------------------------------------------------------------------------

fn run_d() -> Result<bool, String> {
    println!("run D (D_high_e): e = 0.285 — the pseudo-synchronous rate IS 1.5 n, capture guaranteed");
    let dir = output::fresh_run_dir("D_high_e")?;
    let e_hi = 0.285;
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, e_hi, 0.0, 0.0, 1.55 * n0],
    };
    let mut det = CaptureDetector::new(CaptureConfig {
        dense_trigger_ratio: 1.535,
        dense_cadence: 2.0 * params::YEAR,
        base_cadence: 50.0 * params::YEAR,
        passed_ratio: 1.40,
        window: 30_000.0 * params::YEAR,
        post_capture_dense: 1.0e18,
        stop_on_decision: false,
        stop_after_capture_span: 60_000.0 * params::YEAR,
    });
    let mut samples: Vec<Sample> = Vec::new();
    let mut stats = Stats::default();
    let mut count = 0usize;
    let mut obs = |s: &Sample| {
        count += 1;
        if count % 5000 == 0 {
            println!(
                "  t = {} kyr  Omega/n = {}",
                fmt_f(s.t / (1000.0 * params::YEAR), 1),
                fmt_f(s.ratio, 5)
            );
        }
        det.observe(s)
    };
    let spec = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: true,
        t_end: 2.0e6 * params::YEAR,
        cadence: 50.0 * params::YEAR,
        root_ratio: 1.5,
        stop_on_root: false,
        reanchor: true,
        stage_tag: 'R',
        record: true,
    };
    let end = integrate_segment(&ic, &spec, &mut obs, &mut samples, &mut stats)?;
    let mut events: Vec<Event> = Vec::new();
    if let Some(cross) = end.root_state {
        events.push(Event {
            t: cross.t,
            name: "cross_3:2".into(),
            value: cross.y[4] / params::mean_motion(cross.y[0]),
        });
    }
    let t_capture = match det.decision {
        Some(Decision::Captured { t_capture }) => {
            events.push(Event { t: t_capture, name: "capture_detected".into(), value: 1.5 });
            Some(t_capture)
        }
        _ => None,
    };

    let mut ok = true;
    ok &= check(
        t_capture.is_some(),
        "D.guaranteed_capture",
        &t_capture
            .map(|t| format!("captured at t = {} kyr", fmt_f(t / (1000.0 * params::YEAR), 1)))
            .unwrap_or_else(|| "no capture — should be impossible at e = 0.285".to_string()),
    );
    if let Some(tc) = t_capture {
        // Libration period vs the Goldreich-Peale pendulum. A fresh capture
        // librates at LARGE amplitude (near the separatrix), so the exact
        // pendulum period T = (4/omega_lib) * K(sin(gamma_max/2)) is used,
        // with K the complete elliptic integral of the first kind (evaluated
        // by the arithmetic-geometric mean) and the amplitude measured from
        // the same data window (documented deviation DEV-2: the plan's
        // small-amplitude formula applies only to settled libration, and is
        // separately unit-tested at small amplitude).
        let h_e = 3.5 * e_hi - 7.6875 * (e_hi * e_hi * e_hi);
        let omega_lib = n0 * (3.0 * params::B_MINUS_A_OVER_C * h_e).sqrt();
        fn elliptic_k(k: f64) -> f64 {
            let mut a = 1.0f64;
            let mut b = (1.0 - k * k).sqrt();
            for _ in 0..60 {
                let an = 0.5 * (a + b);
                b = (a * b).sqrt();
                a = an;
            }
            PI / (2.0 * a)
        }
        let win: Vec<&Sample> = samples
            .iter()
            .filter(|s| s.t >= tc + 30_000.0 * params::YEAR && s.t <= tc + 50_000.0 * params::YEAR)
            .collect();
        let mut measured = 0.0;
        let mut p_pred = 0.0;
        let mut amp = 0.0;
        if win.len() > 100 {
            let meang: f64 = win.iter().map(|s| s.gamma).sum::<f64>() / (win.len() as f64);
            amp = win
                .iter()
                .map(|s| (s.gamma - meang).abs())
                .fold(0.0f64, f64::max);
            p_pred = 4.0 * elliptic_k((0.5 * amp).sin()) / omega_lib;
            let mut crossings: Vec<f64> = Vec::new();
            for pair in win.windows(2) {
                let (p, q) = (pair[0], pair[1]);
                if (p.gamma - meang) <= 0.0 && (q.gamma - meang) > 0.0 {
                    crossings.push(q.t);
                }
            }
            if crossings.len() >= 2 {
                measured = (crossings[crossings.len() - 1] - crossings[0])
                    / ((crossings.len() - 1) as f64);
            }
        }
        let rel = if measured > 0.0 && p_pred > 0.0 {
            (measured - p_pred).abs() / p_pred
        } else {
            1.0
        };
        ok &= check(
            rel < 0.15,
            "D.libration_period_vs_goldreich_peale",
            &format!(
                "measured {} yr vs pendulum-exact {} yr at amplitude {} rad (small-amplitude formula: {} yr; rel diff {})",
                fmt_f(measured / params::YEAR, 2),
                fmt_f(p_pred / params::YEAR, 2),
                fmt_f(amp, 3),
                fmt_f(2.0 * PI / omega_lib / params::YEAR, 2),
                fmt_e(rel, 3)
            ),
        );
    }
    let rows = output::write_samples(&dir, &samples)?;
    output::write_events(&dir, &events)?;
    println!("  wrote {rows} samples, {} events", events.len());
    let mut m = base_manifest(
        "D_high_e",
        "Guaranteed-capture encore at e = 0.285 (pseudo-synchronous rate = 1.5 n) — why eccentricity is the secret",
        k2tau_movie(),
        &ic,
    );
    m.extras.push(("cadence_base_s".into(), fmt_e(50.0 * params::YEAR, 6)));
    m.extras.push(("cadence_dense_s".into(), fmt_e(2.0 * params::YEAR, 6)));
    finish_manifest(&dir, m, &stats, end.state.t, ok)
}

// ---------------------------------------------------------------------------
// Run E — the staging-seam validation at Omega/n = 3
// ---------------------------------------------------------------------------

fn decile_rate(samples: &[Sample], pick: fn(&Sample) -> f64) -> f64 {
    let n = samples.len();
    if n < 2 {
        return f64::NAN; // caller checks sample counts first (run_e)
    }
    let k = (n / 10).max(1);
    let mean = |sl: &[Sample]| -> (f64, f64) {
        let m = sl.iter().map(pick).sum::<f64>() / (sl.len() as f64);
        let t = sl.iter().map(|s| s.t).sum::<f64>() / (sl.len() as f64);
        (m, t)
    };
    let (v1, t1) = mean(&samples[..k]);
    let (v2, t2) = mean(&samples[n - k..]);
    (v2 - v1) / (t2 - t1)
}

fn run_e() -> Result<bool, String> {
    // NOTE (documented deviation DEV-1): the plan first placed this test at
    // Omega/n = 3 — which is itself the 3:1 spin-orbit resonance (resonances
    // sit at every half-integer ratio), so the full model locked there and
    // the drift comparison was meaningless. The seam is validated at the
    // non-resonant ratio 2.7 instead.
    println!("run E (E_seam): secular vs full model over the same 5000-yr window at Omega/n = 2.7");
    let dir = output::fresh_run_dir("E_seam")?;
    let n0 = params::mean_motion(params::A0);
    let ic = State {
        t: 0.0,
        y: [params::A0, params::E0, 0.0, 0.0, 2.7 * n0],
    };
    let t_end = 5000.0 * params::YEAR;
    let mut stats = Stats::default();
    let mut noop = |_s: &Sample| ObserverCmd::Continue;

    let mut sec: Vec<Sample> = Vec::new();
    let spec_s = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: false,
        t_end,
        cadence: 10.0 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    integrate_segment(&ic, &spec_s, &mut noop, &mut sec, &mut stats)?;

    let mut full: Vec<Sample> = Vec::new();
    let spec_r = SegmentSpec {
        k2tau: k2tau_movie(),
        triaxial_on: true,
        t_end,
        cadence: 10.0 * params::YEAR,
        root_ratio: -1.0,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'R',
        record: true,
    };
    integrate_segment(&ic, &spec_r, &mut noop, &mut full, &mut stats)?;
    if sec.len() < 20 || full.len() < 20 {
        return Err(format!(
            "run E: too few samples for a drift comparison (secular {}, full {})",
            sec.len(),
            full.len()
        ));
    }

    let mut ok = true;
    // The staged models differ ONLY in the dOmega/dt equation (the handle
    // torque); da/dt and de/dt are the identical code in both. So the seam
    // is validated on the despin rate, plus end-state sanity on a and e.
    let rs = decile_rate(&sec, |s| s.omega);
    let rr = decile_rate(&full, |s| s.omega);
    let rel = (rs - rr).abs() / rs.abs().max(1.0e-300);
    ok &= check(
        rel < 0.01,
        "E.seam_dOmega_dt",
        &format!(
            "secular {} vs full {} rad/s^2 (rel diff {})",
            fmt_e(rs, 6),
            fmt_e(rr, 6),
            fmt_e(rel, 3)
        ),
    );
    let (a_s, a_r) = (sec.last().map(|s| s.a).unwrap_or(0.0), full.last().map(|s| s.a).unwrap_or(0.0));
    let (e_s, e_r) = (sec.last().map(|s| s.e).unwrap_or(0.0), full.last().map(|s| s.e).unwrap_or(0.0));
    ok &= check(
        (a_s - a_r).abs() < 5.0,
        "E.seam_a_end",
        &format!("end-state |a_S - a_R| = {} m (bound 5 m over 5000 yr)", fmt_e((a_s - a_r).abs(), 3)),
    );
    ok &= check(
        (e_s - e_r).abs() < 1.0e-10,
        "E.seam_e_end",
        &format!("end-state |e_S - e_R| = {} (bound 1e-10)", fmt_e((e_s - e_r).abs(), 3)),
    );
    let mut samples = sec;
    samples.extend(full);
    let rows = output::write_samples(&dir, &samples)?;
    output::write_events(&dir, &[])?;
    println!("  wrote {rows} samples (both stages), 0 events");
    let mut m = base_manifest(
        "E_seam",
        "Staging-seam validation: secular vs full model despin rates over 5000 yr at the non-resonant Omega/n = 2.7 (DEV-1)",
        k2tau_movie(),
        &ic,
    );
    m.extras.push(("window_s".into(), fmt_e(t_end, 6)));
    m.extras.push(("cadence_s".into(), fmt_e(10.0 * params::YEAR, 6)));
    finish_manifest(&dir, m, &stats, t_end, ok)
}

// ---------------------------------------------------------------------------
// TEST 2 — Jupiter + Einstein (six-state runs; see src/test2.rs)
// ---------------------------------------------------------------------------

use mercury_rs::test2;
use mercury_rs::test2::{Cmd2, Detector2, Sample2, Segment2, State2, Stats2};

fn manifest2(run_id: &str, description: &str, k2tau: f64, ic: &State2) -> output::Manifest {
    let mut m = base_manifest(
        run_id,
        description,
        k2tau.max(1.0e-300),
        &State {
            t: ic.t,
            y: [ic.y[0], ic.y[1], ic.y[2], ic.y[4], ic.y[5]],
        },
    );
    if k2tau <= 0.0 {
        m.tau_lag_s = 0.0;
        m.compression = 0.0;
    }
    m.extras.push(("pomega0_rad".into(), fmt_e(ic.y[3], 17)));
    m.extras.push(("m_jupiter_kg".into(), fmt_e(test2::M_JUP, 6)));
    m.extras.push(("a_jupiter_m".into(), fmt_e(test2::A_JUP, 6)));
    m.extras.push(("e_jupiter".into(), fmt_e(test2::E_JUP, 6)));
    let ll = test2::ll_rates();
    m.extras.push(("ll_A11_rad_s".into(), fmt_e(ll.a11, 8)));
    m.extras.push(("ll_A12_rad_s".into(), fmt_e(ll.a12, 8)));
    m.extras.push((
        "gr_rate_rad_s".into(),
        fmt_e(test2::gr_pomega_dot(params::A0, params::E0), 8),
    ));
    m
}

fn finish2(
    dir: &std::path::Path,
    m: output::Manifest,
    stats: &Stats2,
    t_end: f64,
    ok: bool,
) -> Result<bool, String> {
    let s = Stats {
        n_steps: stats.n_steps,
        n_rhs: stats.n_rhs,
        n_reanchor: stats.n_reanchor,
    };
    finish_manifest(dir, m, &s, t_end, ok)
}

/// TEST 2 gate A — Einstein alone must reproduce the famous 43"/century.
fn t2_gr_check() -> Result<bool, String> {
    println!("test 2 gate A (T2_gr_check): GR perihelion advance vs the famous 43\"/century");
    let dir = output::fresh_run_dir("T2_gr_check")?;
    let ic = test2::initial_state2();
    let p = test2::Rhs2Params {
        k2tau: 0.0,
        triaxial_on: false,
        gr_on: true,
        jupiter_on: false,
        a11: 0.0,
        a12: 0.0,
        root_ratio: -1.0,
    };
    let mut samples: Vec<Sample2> = Vec::new();
    let mut stats = Stats2::default();
    let mut noop = |_s: &Sample2| Cmd2::Continue;
    let spec = Segment2 {
        p,
        t_end: 1000.0 * params::YEAR,
        cadence: 10.0 * params::YEAR,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let end = test2::integrate_segment6(&ic, &spec, &mut noop, &mut samples, &mut stats)?;
    let measured = end.state.y[3] / end.state.t;
    let predicted = test2::gr_pomega_dot(params::A0, params::E0);
    let rel = (measured - predicted).abs() / predicted;
    let mut ok = true;
    ok &= check(
        rel < 1.0e-3,
        "T2a.gr_precession_rate",
        &format!(
            "measured {} \"/century vs analytic {} \"/century (rel diff {})",
            fmt_f(test2::arcsec_cy(measured), 3),
            fmt_f(test2::arcsec_cy(predicted), 3),
            fmt_e(rel, 3)
        ),
    );
    let rows = test2::write_samples2(&dir, &samples)?;
    output::write_events(&dir, &[])?;
    println!("  wrote {rows} samples");
    let m = manifest2(
        "T2_gr_check",
        "TEST 2 gate A: GR apsidal precession alone (tides, Jupiter, triaxial all off)",
        0.0,
        &ic,
    );
    finish2(&dir, m, &stats, end.state.t, ok)
}

/// TEST 2 gate B — Jupiter alone must reproduce the Laplace-Lagrange forced
/// eccentricity oscillation (amplitude |A12/A11| e_J, period 2 pi / A11).
fn t2_jupiter_check() -> Result<bool, String> {
    println!("test 2 gate B (T2_jupiter): Jupiter's secular forcing vs Laplace-Lagrange theory");
    let dir = output::fresh_run_dir("T2_jupiter")?;
    let ic = test2::initial_state2();
    let ll = test2::ll_rates();
    let p = test2::Rhs2Params {
        k2tau: 0.0,
        triaxial_on: false,
        gr_on: false,
        jupiter_on: true,
        a11: ll.a11,
        a12: ll.a12,
        root_ratio: -1.0,
    };
    let mut samples: Vec<Sample2> = Vec::new();
    let mut stats = Stats2::default();
    let mut noop = |_s: &Sample2| Cmd2::Continue;
    let spec = Segment2 {
        p,
        t_end: 2.0e6 * params::YEAR,
        cadence: 200.0 * params::YEAR,
        stop_on_root: false,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let end = test2::integrate_segment6(&ic, &spec, &mut noop, &mut samples, &mut stats)?;
    let e_min = samples.iter().map(|s| s.e).fold(f64::INFINITY, f64::min);
    let e_max = samples.iter().map(|s| s.e).fold(f64::NEG_INFINITY, f64::max);
    let amp_measured = 0.5 * (e_max - e_min);
    let amp_predicted = (ll.a12 / ll.a11).abs() * test2::E_JUP;
    let mut maxima: Vec<f64> = Vec::new();
    for w in samples.windows(3) {
        if w[1].e > w[0].e && w[1].e > w[2].e {
            maxima.push(w[1].t);
        }
    }
    let period_measured = if maxima.len() >= 2 {
        (maxima[maxima.len() - 1] - maxima[0]) / ((maxima.len() - 1) as f64)
    } else {
        0.0
    };
    let period_predicted = 2.0 * PI / ll.a11;
    let mut ok = true;
    ok &= check(
        (amp_measured - amp_predicted).abs() / amp_predicted < 0.05,
        "T2b.forced_eccentricity_amplitude",
        &format!(
            "e oscillates {} .. {} -> amplitude {} vs LL prediction {}",
            fmt_f(e_min, 5),
            fmt_f(e_max, 5),
            fmt_e(amp_measured, 4),
            fmt_e(amp_predicted, 4)
        ),
    );
    ok &= check(
        period_measured > 0.0
            && (period_measured - period_predicted).abs() / period_predicted < 0.05,
        "T2b.secular_period",
        &format!(
            "measured {} kyr vs LL prediction {} kyr (Jupiter-driven perihelion circulation)",
            fmt_f(period_measured / (1000.0 * params::YEAR), 1),
            fmt_f(period_predicted / (1000.0 * params::YEAR), 1)
        ),
    );
    let rows = test2::write_samples2(&dir, &samples)?;
    output::write_events(&dir, &[])?;
    println!("  wrote {rows} samples");
    let m = manifest2(
        "T2_jupiter",
        "TEST 2 gate B: Jupiter's Laplace-Lagrange secular forcing alone (tides, GR, triaxial off)",
        0.0,
        &ic,
    );
    finish2(&dir, m, &stats, end.state.t, ok)
}

/// TEST 2 braking movie: full physics, to the restart save at Omega/n = 1.6.
fn t2_movie() -> Result<bool, String> {
    println!("test 2 (T2_movie): braking with tides x1000 + GR + Jupiter, to the restart at 1.6");
    let dir = output::fresh_run_dir("T2_movie")?;
    let ic = test2::initial_state2();
    let mut samples: Vec<Sample2> = Vec::new();
    let mut events: Vec<Event> = Vec::new();
    let mut stats = Stats2::default();

    let mut count = 0usize;
    let mut obs1 = |s: &Sample2| {
        count += 1;
        if count % 500 == 0 {
            println!(
                "  [S] t = {} Myr  Omega/n = {}  e = {}",
                fmt_f(myr(s.t), 3),
                fmt_f(s.ratio, 4),
                fmt_f(s.e, 5)
            );
        }
        Cmd2::Continue
    };
    let spec1 = Segment2 {
        p: test2::full_params(false, params::STAGE_HANDOVER_RATIO),
        t_end: 8.0e6 * params::YEAR,
        cadence: 1000.0 * params::YEAR,
        stop_on_root: true,
        reanchor: false,
        stage_tag: 'S',
        record: true,
    };
    let end1 = test2::integrate_segment6(&ic, &spec1, &mut obs1, &mut samples, &mut stats)?;
    let mut seen_52 = false;
    for s in &samples {
        if !seen_52 && s.ratio < 2.5 {
            seen_52 = true;
            events.push(Event { t: s.t, name: "cross_5:2".into(), value: s.ratio });
        }
    }
    let handover = end1
        .root_state
        .ok_or_else(|| "T2_movie: the 2.2 handover root never fired within 8 Myr".to_string())?;
    events.push(Event {
        t: handover.t,
        name: "stage_handover".into(),
        value: handover.y[5] / params::mean_motion(handover.y[0]),
    });
    println!("  stage handover at t = {} Myr", fmt_f(myr(handover.t), 4));

    let start2 = handover.reanchored();
    let mut seen_21_at: Option<(f64, f64)> = None;
    let mut count2 = 0usize;
    let mut obs2 = |s: &Sample2| {
        count2 += 1;
        if seen_21_at.is_none() && s.ratio < 2.0 {
            seen_21_at = Some((s.t, s.ratio));
        }
        if count2 % 2000 == 0 {
            println!(
                "  [R] t = {} Myr  Omega/n = {}  e = {}",
                fmt_f(myr(s.t), 3),
                fmt_f(s.ratio, 4),
                fmt_f(s.e, 5)
            );
        }
        Cmd2::Continue
    };
    let spec2 = Segment2 {
        p: test2::full_params(true, params::RESTART_RATIO),
        t_end: handover.t + 3.0e6 * params::YEAR,
        cadence: 100.0 * params::YEAR,
        stop_on_root: true,
        reanchor: true,
        stage_tag: 'R',
        record: true,
    };
    let end2 = test2::integrate_segment6(&start2, &spec2, &mut obs2, &mut samples, &mut stats)?;
    if let Some((t, v)) = seen_21_at {
        events.push(Event { t, name: "cross_2:1".into(), value: v });
    }
    let restart_raw = end2.root_state.ok_or_else(|| {
        format!(
            "T2_movie: the 1.6 restart root never fired (final ratio {})",
            fmt_f(end2.state.y[5] / params::mean_motion(end2.state.y[0]), 4)
        )
    })?;
    let restart = restart_raw.reanchored();
    test2::write_restart6(&dir, &restart)?;
    let restart_ratio = restart.y[5] / params::mean_motion(restart.y[0]);
    events.push(Event { t: restart.t, name: "restart_saved".into(), value: restart_ratio });

    let mut ok = true;
    ok &= check(
        (restart_ratio - params::RESTART_RATIO).abs() < 1.0e-6,
        "T2m.restart_ratio_exact",
        &format!("root-found Omega/n = {} (target 1.6)", fmt_e(restart_ratio, 12)),
    );
    ok &= check(
        (2.0e6..=7.5e6).contains(&(restart.t / params::YEAR)),
        "T2m.restart_time_plausible",
        &format!("restart at t = {} Myr", fmt_f(myr(restart.t), 4)),
    );
    ok &= check(
        seen_52 && seen_21_at.is_some(),
        "T2m.crossings_recorded",
        "5:2 and 2:1 crossings in the event log",
    );
    let rows = test2::write_samples2(&dir, &samples)?;
    output::write_events(&dir, &events)?;
    println!("  wrote {rows} samples, {} events", events.len());
    let m = manifest2(
        "T2_movie",
        "TEST 2 braking movie: tides x1000 + GR + Jupiter secular terms, to the restart at Omega/n = 1.6",
        k2tau_movie(),
        &ic,
    );
    finish2(&dir, m, &stats, restart.t, ok)
}

fn t2_one_branch(
    restart: State2,
    branch_id: usize,
    n_branches: usize,
) -> Result<(BranchOutcome, Stats2), String> {
    let offset = (branch_id as f64) * (PI / (n_branches as f64));
    let mut st = restart;
    st.y[4] += offset;
    let mut det = Detector2::standard(true);
    let mut stats = Stats2::default();
    let mut samples: Vec<Sample2> = Vec::new();
    let spec = Segment2 {
        p: test2::full_params(true, 1.5),
        t_end: restart.t + 1.5e6 * params::YEAR,
        cadence: 100.0 * params::YEAR,
        stop_on_root: false,
        reanchor: true,
        stage_tag: 'R',
        record: false,
    };
    let mut obs = |s: &Sample2| det.observe(s);
    let end = test2::integrate_segment6(&st, &spec, &mut obs, &mut samples, &mut stats)?;
    let final_ratio = end.state.y[5] / params::mean_motion(end.state.y[0]);
    match det.decision {
        Some(Decision::Captured { t_capture }) => Ok((
            BranchOutcome {
                branch_id,
                theta_offset: offset,
                captured: true,
                t_outcome: t_capture,
                final_ratio,
            },
            stats,
        )),
        Some(Decision::Passed { t, ratio }) => Ok((
            BranchOutcome {
                branch_id,
                theta_offset: offset,
                captured: false,
                t_outcome: t,
                final_ratio: ratio,
            },
            stats,
        )),
        None => Err(format!(
            "T2 branch {branch_id}: undecided after 1.5 Myr (final ratio {})",
            fmt_f(final_ratio, 4)
        )),
    }
}

fn t2_sweep(n_branches: usize) -> Result<bool, String> {
    println!("test 2 (T2_sweep): {n_branches}-branch phase sweep with Jupiter + GR active");
    let b_dir = output::run_dir("T2_movie")?;
    let restart = test2::read_restart6(&b_dir)?;
    let dir = output::fresh_run_dir("T2_sweep")?;
    let workers = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(8);
    println!("  restart: t = {} Myr", fmt_f(myr(restart.t), 4));
    let (tx, rx) =
        std::sync::mpsc::channel::<(usize, Result<(BranchOutcome, Stats2), String>)>();
    let mut spawned = 0usize;
    for w in 0..workers {
        let tx = tx.clone();
        let ids: Vec<usize> = (0..n_branches).filter(|k| k % workers == w).collect();
        spawned += ids.len();
        std::thread::spawn(move || {
            for k in ids {
                let res = t2_one_branch(restart, k, n_branches);
                if tx.send((k, res)).is_err() {
                    return;
                }
            }
        });
    }
    drop(tx);
    let mut collected: Vec<(usize, Result<(BranchOutcome, Stats2), String>)> =
        rx.iter().take(spawned).collect();
    collected.sort_by_key(|(k, _)| *k);
    let mut outcomes: Vec<BranchOutcome> = Vec::with_capacity(n_branches);
    let mut stats = Stats2::default();
    for (k, res) in collected {
        let (b, s) = res.map_err(|e| format!("branch {k}: {e}"))?;
        outcomes.push(b);
        stats.n_steps += s.n_steps;
        stats.n_rhs += s.n_rhs;
        stats.n_reanchor += s.n_reanchor;
    }
    outcomes.sort_by_key(|b| b.branch_id);
    let captured: Vec<usize> =
        outcomes.iter().filter(|b| b.captured).map(|b| b.branch_id).collect();
    let canonical = captured.first().copied();
    for b in &outcomes {
        println!(
            "  branch {:2}  offset = {} rad  {}  t = {} Myr  final ratio = {}",
            b.branch_id,
            fmt_e(b.theta_offset, 6),
            if b.captured { "CAPTURED" } else { "passed  " },
            fmt_f(myr(b.t_outcome), 4),
            fmt_f(b.final_ratio, 4)
        );
    }
    println!(
        "  capture fraction with Jupiter's eccentricity cycle in play: {}/{}",
        captured.len(),
        outcomes.len()
    );
    let mut ok = true;
    ok &= check(
        outcomes.len() == n_branches,
        "T2s.all_branches_decided",
        &format!("{} of {n_branches}", outcomes.len()),
    );
    let any = !captured.is_empty();
    ok &= check(
        any,
        "T2s.at_least_one_capture",
        if any {
            "canonical branch selected"
        } else {
            "ZERO captures — run the documented finer-grid contingency re-sweep"
        },
    );
    output::write_branches(&dir, &outcomes, canonical)?;
    output::write_events(&dir, &[])?;
    let mut m = manifest2(
        "T2_sweep",
        "TEST 2 phase sweep at the 3:2 crossing with GR + Jupiter active",
        k2tau_movie(),
        &restart,
    );
    m.extras.push(("branches".into(), format!("{n_branches}")));
    m.extras.push((
        "canonical_branch".into(),
        canonical.map(|k| format!("{k}")).unwrap_or_else(|| "-1".to_string()),
    ));
    let ok = finish2(&dir, m, &stats, restart.t + 1.5e6 * params::YEAR, ok)?;
    if !any {
        println!("FAILURE");
        std::process::exit(3);
    }
    Ok(ok)
}

/// TEST 2 canonical continuation: capture and lock WITH the precessing
/// ellipse — the settled mean spin ratio must be 1.5 + pomega_dot/n.
fn t2_final(branch: usize, n_branches: usize) -> Result<bool, String> {
    if branch >= n_branches {
        return Err(format!(
            "t2-final: --branch {branch} was never swept (valid: 0..{})",
            n_branches - 1
        ));
    }
    println!("test 2 (T2_final): branch {branch} continued to 10 Myr with Jupiter + GR");
    let b_dir = output::run_dir("T2_movie")?;
    let restart = test2::read_restart6(&b_dir)?;
    let dir = output::fresh_run_dir("T2_final")?;
    let offset = (branch as f64) * (PI / (n_branches as f64));
    let mut st = restart;
    st.y[4] += offset;

    let mut det = Detector2::standard(false);
    let mut samples: Vec<Sample2> = Vec::new();
    let mut stats = Stats2::default();
    let mut count = 0usize;
    let mut obs = |s: &Sample2| {
        count += 1;
        if count % 5000 == 0 {
            println!(
                "  t = {} Myr  Omega/n = {}  e = {}",
                fmt_f(myr(s.t), 3),
                fmt_f(s.ratio, 6),
                fmt_f(s.e, 5)
            );
        }
        det.observe(s)
    };
    let spec = Segment2 {
        p: test2::full_params(true, 1.5),
        t_end: params::T_FINAL,
        cadence: 100.0 * params::YEAR,
        stop_on_root: false,
        reanchor: true,
        stage_tag: 'R',
        record: true,
    };
    let end = test2::integrate_segment6(&st, &spec, &mut obs, &mut samples, &mut stats)?;

    let mut events: Vec<Event> = Vec::new();
    if let Some(cross) = end.root_state {
        events.push(Event {
            t: cross.t,
            name: "cross_3:2".into(),
            value: cross.y[5] / params::mean_motion(cross.y[0]),
        });
    }
    let t_capture = match det.decision {
        Some(Decision::Captured { t_capture }) => {
            events.push(Event { t: t_capture, name: "capture_detected".into(), value: 1.5 });
            Some(t_capture)
        }
        _ => None,
    };
    events.push(Event { t: end.state.t, name: "reanchor".into(), value: stats.n_reanchor as f64 });

    let n0 = params::mean_motion(params::A0);
    let pw_dot_pred = test2::predicted_pomega_dot();
    let expected_ratio = 1.5 + pw_dot_pred / n0;

    let mut ok = true;
    ok &= check(
        t_capture.is_some(),
        "T2f.capture_detected",
        &t_capture
            .map(|t| format!("captured at t = {} Myr with Jupiter + GR active", fmt_f(myr(t), 4)))
            .unwrap_or_else(|| "branch did not capture — reselect the canonical branch".to_string()),
    );
    if let Some(tc) = t_capture {
        let settled: Vec<&Sample2> =
            samples.iter().filter(|s| s.t > tc + 1.0e6 * params::YEAR).collect();
        if settled.is_empty() {
            ok = check(false, "T2f.locked_mean_ratio", "no settled-era samples");
        } else {
            let mean_ratio: f64 =
                settled.iter().map(|s| s.ratio).sum::<f64>() / (settled.len() as f64);
            ok &= check(
                (mean_ratio - expected_ratio).abs() <= 1.5e-7,
                "T2f.lock_follows_the_precessing_ellipse",
                &format!(
                    "settled mean Omega/n = {} vs predicted 1.5 + pomega_dot/n = {} (GR {} + Jupiter {} \"/century)",
                    fmt_e(mean_ratio, 10),
                    fmt_e(expected_ratio, 10),
                    fmt_f(test2::arcsec_cy(test2::gr_pomega_dot(params::A0, params::E0)), 2),
                    fmt_f(test2::arcsec_cy(test2::ll_rates().a11), 2)
                ),
            );
            ok &= check(
                (mean_ratio - 1.5) >= 2.0e-7,
                "T2f.shift_off_exact_three_halves",
                &format!(
                    "mean ratio exceeds exactly 1.5 by {} — the lock tracks the precessing perihelion, not the stars",
                    fmt_e(mean_ratio - 1.5, 4)
                ),
            );
        }
        let e_min = samples.iter().map(|s| s.e).fold(f64::INFINITY, f64::min);
        let e_max = samples.iter().map(|s| s.e).fold(f64::NEG_INFINITY, f64::max);
        let ll = test2::ll_rates();
        let amp_pred = (ll.a12 / ll.a11).abs() * test2::E_JUP;
        ok &= check(
            ((0.5 * (e_max - e_min)) - amp_pred).abs() / amp_pred < 0.3,
            "T2f.jupiter_eccentricity_cycle",
            &format!(
                "e oscillates {} .. {} through braking, capture, and lock (LL amplitude {})",
                fmt_f(e_min, 5),
                fmt_f(e_max, 5),
                fmt_e(amp_pred, 4)
            ),
        );
        let pw_rate = (end.state.y[3] - st.y[3]) / (end.state.t - st.t);
        ok &= check(
            (pw_rate - pw_dot_pred).abs() / pw_dot_pred < 0.05,
            "T2f.perihelion_advance",
            &format!(
                "perihelion advanced {} rad over the run = {} \"/century (predicted {}; 43.0 of it is Einstein's)",
                fmt_f(end.state.y[3] - st.y[3], 2),
                fmt_f(test2::arcsec_cy(pw_rate), 2),
                fmt_f(test2::arcsec_cy(pw_dot_pred), 2)
            ),
        );
        let last = samples.last().ok_or_else(|| "T2_final produced no samples".to_string())?;
        ok &= check(
            (last.p_orb / DAY - OBSERVED_P_ORB_D).abs() / OBSERVED_P_ORB_D < 1.0e-3,
            "T2f.P_orb",
            &format!("final P_orb = {} d (observed 87.969)", fmt_f(last.p_orb / DAY, 4)),
        );
        println!("  NOTE: angular momentum is deliberately NOT ledger-checked in test 2 — the");
        println!("  Laplace-Lagrange terms exchange orbital angular momentum with Jupiter,");
        println!("  which this model does not track.");
    } else {
        ok = false;
    }

    let rows = test2::write_samples2(&dir, &samples)?;
    output::write_events(&dir, &events)?;
    println!("  wrote {rows} samples, {} events", events.len());
    let mut m = manifest2(
        "T2_final",
        "TEST 2 canonical continuation: capture and lock with GR + Jupiter — the lock follows the precessing ellipse",
        k2tau_movie(),
        &st,
    );
    m.extras.push(("branch".into(), format!("{branch}")));
    m.extras.push(("theta_offset_rad".into(), fmt_e(offset, 17)));
    finish2(&dir, m, &stats, end.state.t, ok)
}

// ---------------------------------------------------------------------------

fn getopt(args: &[String], name: &str, default: usize) -> usize {
    for i in 0..args.len() {
        if args[i] == name {
            match args.get(i + 1).map(|s| s.parse::<usize>()) {
                Some(Ok(v)) => return v,
                _ => {
                    eprintln!(
                        "error: {name} needs a whole-number value, got {:?}",
                        args.get(i + 1)
                    );
                    std::process::exit(2);
                }
            }
        }
    }
    default
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let sub = args.get(1).map(|s| s.as_str()).unwrap_or("");
    let result: Result<bool, String> = match sub {
        "print-config" => {
            print_config();
            return;
        }
        "run-a" => run_a(),
        "run-b" => run_b(),
        "sweep" => run_sweep(getopt(&args, "--branches", params::SWEEP_BRANCHES)),
        "run-b-final" => {
            let k = getopt(&args, "--branch", usize::MAX);
            if k == usize::MAX {
                Err("run-b-final needs --branch <k> (the captured branch from the sweep)".to_string())
            } else {
                run_b_final(k, getopt(&args, "--branches", params::SWEEP_BRANCHES))
            }
        }
        "run-d" => run_d(),
        "run-e" => run_e(),
        "t2-gr-check" => t2_gr_check(),
        "t2-jupiter-check" => t2_jupiter_check(),
        "t2-movie" => t2_movie(),
        "t2-sweep" => t2_sweep(getopt(&args, "--branches", 16)),
        "t2-final" => {
            let k = getopt(&args, "--branch", usize::MAX);
            if k == usize::MAX {
                Err("t2-final needs --branch <k> (the captured branch from t2-sweep)".to_string())
            } else {
                t2_final(k, getopt(&args, "--branches", 16))
            }
        }
        _ => {
            eprintln!(
                "usage: mercury_rs <print-config|run-a|run-b|sweep [--branches N]|run-b-final --branch K|run-d|run-e|t2-gr-check|t2-jupiter-check|t2-movie|t2-sweep [--branches N]|t2-final --branch K>"
            );
            std::process::exit(2);
        }
    };
    match result {
        Ok(true) => {
            println!("SUCCESS");
        }
        Ok(false) => {
            println!("FAILURE");
            std::process::exit(1);
        }
        Err(e) => {
            eprintln!("error: {e}");
            println!("FAILURE");
            std::process::exit(1);
        }
    }
}
