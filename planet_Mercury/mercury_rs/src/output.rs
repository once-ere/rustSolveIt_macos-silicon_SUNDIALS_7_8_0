//! CSV and manifest writers/readers. All floating-point text goes through the
//! engine's C-style formatters (fmt_e), never Rust's `{:e}` — re-runs must be
//! byte-identical.
//!
//! Files per run directory data/runs/<run_id>/:
//!   samples.csv  — the time history (header matches the database exactly)
//!   events.csv   — notable moments (t_s,event,value)
//!   branches.csv — sweep outcomes (run C only)
//!   restart.csv  — full-precision restart state (run B only)
//!   manifest.json— configuration echo + solver statistics + verdict

use std::fs;
use std::io::Write as IoWrite;
use std::path::{Path, PathBuf};

use sundials_core::sundials_utils::fmt_e;

use crate::driver::{BranchOutcome, Event, Sample, State};
use crate::params;

/// The samples.csv header — one for one with the database `sample` table.
pub const SAMPLES_HEADER: &str = "t_s,a_m,e,M_rad,theta_rad,Omega_rad_s,n_rad_s,ratio,gamma_rad,P_orb_s,P_rot_s,L_spin_kgm2s,L_orb_kgm2s,L_tot_kgm2s,E_spin_j,E_orb_j,stage";

/// Resolve the data directory: $MERCURY_DATA_DIR, else ../data/runs relative
/// to the current directory (works from mercury_rs/ and from notebook/).
pub fn data_dir() -> PathBuf {
    match std::env::var("MERCURY_DATA_DIR") {
        Ok(d) if !d.is_empty() => PathBuf::from(d),
        _ => PathBuf::from("../data/runs"),
    }
}

/// Create (if needed) and return the directory for one run — READ access:
/// never touches existing contents (a run reading a SIBLING's directory,
/// like the sweep reading B_movie's restart, must not disturb it).
///
/// Data-shape contract for samples.csv (review documentation): rows sit on
/// the requested output cadence; the instants where a CVODE root fires
/// (stage handover, restart save, first 3:2 crossing) are recorded in
/// events.csv / restart.csv, NOT as samples.csv rows — and when a root lands
/// within 1 s of a scheduled output time, the subsequent grid re-anchors to
/// root-time + k*cadence. Consumers must key on t_s, not on row index
/// arithmetic.
pub fn run_dir(run_id: &str) -> Result<PathBuf, String> {
    let dir = data_dir().join(run_id);
    fs::create_dir_all(&dir)
        .map_err(|e| format!("cannot create run directory {}: {e}", dir.display()))?;
    Ok(dir)
}

/// The OWNING run's entry point: like run_dir, but also removes any manifest
/// from a previous run so a later failure can never leave a stale SUCCESS
/// verdict beside fresh outputs (review hardening). Call this ONLY for the
/// directory the current run itself writes — the first build of this
/// hardening cleared manifests on mere read access and thereby deleted
/// B_movie's manifest during the sweep (caught by the end-to-end pipeline;
/// recorded as deviation DEV-9).
pub fn fresh_run_dir(run_id: &str) -> Result<PathBuf, String> {
    let dir = run_dir(run_id)?;
    let stale = dir.join("manifest.json");
    if stale.exists() {
        fs::remove_file(&stale)
            .map_err(|e| format!("cannot remove stale {}: {e}", stale.display()))?;
    }
    Ok(dir)
}

fn open_writer(path: &Path) -> Result<std::io::BufWriter<fs::File>, String> {
    let f = fs::File::create(path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    Ok(std::io::BufWriter::new(f))
}

/// Write samples.csv; returns the number of rows written.
pub fn write_samples(dir: &Path, samples: &[Sample]) -> Result<usize, String> {
    let path = dir.join("samples.csv");
    let mut w = open_writer(&path)?;
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "{SAMPLES_HEADER}").map_err(werr)?;
    for s in samples {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            fmt_e(s.t, 12),
            fmt_e(s.a, 12),
            fmt_e(s.e, 12),
            fmt_e(s.m_anom, 12),
            fmt_e(s.theta, 12),
            fmt_e(s.omega, 12),
            fmt_e(s.n, 12),
            fmt_e(s.ratio, 12),
            fmt_e(s.gamma, 12),
            fmt_e(s.p_orb, 12),
            fmt_e(s.p_rot, 12),
            fmt_e(s.l_spin, 12),
            fmt_e(s.l_orb, 12),
            fmt_e(s.l_tot, 12),
            fmt_e(s.e_spin, 12),
            fmt_e(s.e_orb, 12),
            s.stage
        )
        .map_err(werr)?;
    }
    w.flush().map_err(werr)?;
    Ok(samples.len())
}

/// Write events.csv (always written, even when empty, so ingest is uniform).
pub fn write_events(dir: &Path, events: &[Event]) -> Result<usize, String> {
    let path = dir.join("events.csv");
    let mut w = open_writer(&path)?;
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "t_s,event,value").map_err(werr)?;
    for ev in events {
        writeln!(w, "{},{},{}", fmt_e(ev.t, 12), ev.name, fmt_e(ev.value, 12)).map_err(werr)?;
    }
    w.flush().map_err(werr)?;
    Ok(events.len())
}

/// Write branches.csv (run C). `canonical` marks the branch continued as
/// B_final (the first captured one).
pub fn write_branches(
    dir: &Path,
    branches: &[BranchOutcome],
    canonical: Option<usize>,
) -> Result<usize, String> {
    let path = dir.join("branches.csv");
    let mut w = open_writer(&path)?;
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "branch_id,theta_offset_rad,captured,t_outcome_s,final_ratio,canonical")
        .map_err(werr)?;
    for b in branches {
        writeln!(
            w,
            "{},{},{},{},{},{}",
            b.branch_id,
            fmt_e(b.theta_offset, 12),
            i32::from(b.captured),
            fmt_e(b.t_outcome, 12),
            fmt_e(b.final_ratio, 12),
            i32::from(canonical == Some(b.branch_id))
        )
        .map_err(werr)?;
    }
    w.flush().map_err(werr)?;
    Ok(branches.len())
}

/// Write restart.csv with full round-trip precision (17 decimals of %.*e).
pub fn write_restart(dir: &Path, state: &State) -> Result<(), String> {
    let path = dir.join("restart.csv");
    let mut w = open_writer(&path)?;
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    writeln!(w, "t_s,a_m,e,M_rad,theta_rad,Omega_rad_s").map_err(werr)?;
    writeln!(
        w,
        "{},{},{},{},{},{}",
        fmt_e(state.t, 17),
        fmt_e(state.y[0], 17),
        fmt_e(state.y[1], 17),
        fmt_e(state.y[2], 17),
        fmt_e(state.y[3], 17),
        fmt_e(state.y[4], 17)
    )
    .map_err(werr)?;
    w.flush().map_err(werr)?;
    Ok(())
}

/// Read restart.csv back, bit-exact (%.17e round-trips every f64).
pub fn read_restart(dir: &Path) -> Result<State, String> {
    let path = dir.join("restart.csv");
    let text = fs::read_to_string(&path)
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
    if vals.len() != 6 {
        return Err(format!(
            "{} data row has {} fields, expected 6",
            path.display(),
            vals.len()
        ));
    }
    Ok(State {
        t: vals[0],
        y: [vals[1], vals[2], vals[3], vals[4], vals[5]],
    })
}

/// Everything the manifest records about one run.
pub struct Manifest {
    pub run_id: String,
    pub description: String,
    pub k2: f64,
    pub tau_lag_s: f64,
    pub compression: f64,
    pub a0_m: f64,
    pub e0: f64,
    pub M0_rad: f64,
    pub theta0_rad: f64,
    pub Omega0_rad_s: f64,
    pub t_final_s: f64,
    pub n_steps: i64,
    pub n_rhs_evals: i64,
    pub n_reanchor: i64,
    pub verdict: String,
    /// Free-form extras (cadences, thresholds, offsets) as (key, value-text).
    pub extras: Vec<(String, String)>,
}

/// Write manifest.json by hand (fixed field order, fmt_e floats — the crate
/// has no JSON dependency and needs none).
pub fn write_manifest(dir: &Path, m: &Manifest) -> Result<(), String> {
    let path = dir.join("manifest.json");
    let mut w = open_writer(&path)?;
    let werr = |e: std::io::Error| format!("write error on {}: {e}", path.display());
    let num = |x: f64| fmt_e(x, 17);
    writeln!(w, "{{").map_err(werr)?;
    writeln!(w, "  \"run_id\": \"{}\",", m.run_id).map_err(werr)?;
    writeln!(w, "  \"description\": \"{}\",", m.description).map_err(werr)?;
    writeln!(w, "  \"k2\": {},", num(m.k2)).map_err(werr)?;
    writeln!(w, "  \"tau_lag_s\": {},", num(m.tau_lag_s)).map_err(werr)?;
    writeln!(w, "  \"compression\": {},", num(m.compression)).map_err(werr)?;
    writeln!(w, "  \"a0_m\": {},", num(m.a0_m)).map_err(werr)?;
    writeln!(w, "  \"e0\": {},", num(m.e0)).map_err(werr)?;
    writeln!(w, "  \"M0_rad\": {},", num(m.M0_rad)).map_err(werr)?;
    writeln!(w, "  \"theta0_rad\": {},", num(m.theta0_rad)).map_err(werr)?;
    writeln!(w, "  \"Omega0_rad_s\": {},", num(m.Omega0_rad_s)).map_err(werr)?;
    writeln!(w, "  \"t_final_s\": {},", num(m.t_final_s)).map_err(werr)?;
    writeln!(w, "  \"rel_tol\": {},", num(params::REL_TOL)).map_err(werr)?;
    writeln!(
        w,
        "  \"abs_tol\": [{}, {}, {}, {}, {}],",
        num(params::ABS_TOL[0]),
        num(params::ABS_TOL[1]),
        num(params::ABS_TOL[2]),
        num(params::ABS_TOL[3]),
        num(params::ABS_TOL[4])
    )
    .map_err(werr)?;
    writeln!(w, "  \"max_step_s\": {},", num(params::MAX_STEP)).map_err(werr)?;
    writeln!(w, "  \"solver\": \"CVODE_BDF_NEWTON_DENSE\",").map_err(werr)?;
    writeln!(w, "  \"n_steps\": {},", m.n_steps).map_err(werr)?;
    writeln!(w, "  \"n_rhs_evals\": {},", m.n_rhs_evals).map_err(werr)?;
    writeln!(w, "  \"n_reanchor\": {},", m.n_reanchor).map_err(werr)?;
    writeln!(w, "  \"verdict\": \"{}\",", m.verdict).map_err(werr)?;
    writeln!(
        w,
        "  \"engine\": \"sundials_rs 7.8.0 (pure Rust, macOS arm64)\","
    )
    .map_err(werr)?;
    writeln!(w, "  \"extras\": {{").map_err(werr)?;
    for (i, (k, v)) in m.extras.iter().enumerate() {
        let comma = if i + 1 < m.extras.len() { "," } else { "" };
        writeln!(w, "    \"{k}\": {v}{comma}").map_err(werr)?;
    }
    writeln!(w, "  }}").map_err(werr)?;
    writeln!(w, "}}").map_err(werr)?;
    w.flush().map_err(werr)?;
    Ok(())
}
