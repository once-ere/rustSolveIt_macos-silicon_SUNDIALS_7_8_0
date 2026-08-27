//! Integration tests for the core_io module group of rebound_rs.
//! Part of rebound_rs, GPL-3.0-or-later.
#![allow(non_snake_case)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. See rebound_rust.md section 17.
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::manual_swap)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::misrefactored_assign_op)]
#![allow(clippy::neg_multiply)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_late_init)]
#![allow(clippy::while_let_loop)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::seek_from_current)]
#![allow(clippy::drop_non_drop)]
#![allow(clippy::redundant_closure)]
#![allow(clippy::approx_constant)]
#![allow(clippy::useless_vec)]
#![allow(clippy::type_complexity)]
use rebound_rs::*;

use std::io::Cursor;
use std::sync::atomic::{AtomicUsize, Ordering};

// =====================================================================
// Shared helpers
// =====================================================================

/// A scratch file path unique to this process and the given tag. Any
/// leftover from a previous run is removed so that
/// `reb_simulation_save_to_file` takes its "file does not exist" branch.
fn temp_path(tag: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!(
        "rebound_rs_core_io_{}_{}.bin",
        std::process::id(),
        tag
    ));
    let _ = std::fs::remove_file(&p);
    p
}

fn path_str(p: &std::path::Path) -> String {
    p.to_str().expect("temp path is valid UTF-8").to_string()
}

/// A simulation whose every serialized field is deterministic (the
/// default `rand_seed` is time+pid based).
fn deterministic_sim() -> reb_simulation {
    let mut r = reb_simulation_create();
    r.rand_seed = 42;
    r.save_messages = 1;
    r
}

/// Star + one planet on a mildly eccentric orbit, moved to the centre of
/// mass frame.
fn two_body(integrator: &str, dt: f64, e: f64) -> reb_simulation {
    let mut r = deterministic_sim();
    r.G = 1.0;
    r.dt = dt;
    reb_simulation_set_integrator(&mut r, integrator);
    reb_simulation_add(
        &mut r,
        reb_particle {
            m: 1.0,
            ..reb_particle::default()
        },
    );
    let planet = reb_particle_from_orbit(r.G, r.particles[0], 1e-3, 1.0, e, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(&mut r, planet);
    reb_simulation_move_to_com(&mut r);
    r
}

/// The 11 doubles of a particle as raw bits, for bit-exactness checks.
fn particle_bits(p: &reb_particle) -> [u64; 11] {
    [
        p.x.to_bits(),
        p.y.to_bits(),
        p.z.to_bits(),
        p.vx.to_bits(),
        p.vy.to_bits(),
        p.vz.to_bits(),
        p.ax.to_bits(),
        p.ay.to_bits(),
        p.az.to_bits(),
        p.m.to_bits(),
        p.r.to_bits(),
    ]
}

fn all_particle_bits(r: &reb_simulation) -> Vec<[u64; 11]> {
    (0..r.N).map(|i| particle_bits(&r.particles[i])).collect()
}

// --- binary buffer walking -------------------------------------------

#[derive(Debug)]
struct Field {
    name: String,
    size_name: u64,
    size_data: u64,
    data_pos: usize,
}

/// Walk a serialized simulation (past the 64-byte header) and return
/// every field in file order, up to and including "end".
fn parse_fields(buf: &[u8]) -> Vec<Field> {
    let mut out = Vec::new();
    let mut pos = 64usize;
    while pos + 16 <= buf.len() {
        let size_name = u64::from_le_bytes(buf[pos..pos + 8].try_into().unwrap());
        let size_data = u64::from_le_bytes(buf[pos + 8..pos + 16].try_into().unwrap());
        pos += 16;
        if pos + (size_name as usize) > buf.len() {
            break;
        }
        let nb = &buf[pos..pos + size_name as usize];
        let nul = nb.iter().position(|&b| b == 0).unwrap_or(nb.len());
        let name = String::from_utf8_lossy(&nb[..nul]).to_string();
        pos += size_name as usize;
        let data_pos = pos;
        pos += size_data as usize;
        let is_end = name == "end";
        out.push(Field {
            name,
            size_name,
            size_data,
            data_pos,
        });
        if is_end {
            break;
        }
    }
    out
}

fn field_names(buf: &[u8]) -> Vec<String> {
    parse_fields(buf).into_iter().map(|f| f.name).collect()
}

fn find_field<'a>(fields: &'a [Field], name: &str) -> &'a Field {
    fields
        .iter()
        .find(|f| f.name == name)
        .unwrap_or_else(|| panic!("field {} missing from binary buffer", name))
}

fn le_f64(buf: &[u8], off: usize) -> f64 {
    f64::from_le_bytes(buf[off..off + 8].try_into().unwrap())
}
fn le_u64(buf: &[u8], off: usize) -> u64 {
    u64::from_le_bytes(buf[off..off + 8].try_into().unwrap())
}
fn le_i32(buf: &[u8], off: usize) -> i32 {
    i32::from_le_bytes(buf[off..off + 4].try_into().unwrap())
}

/// The three int32 of a `reb_simulationarchive_blob` at `off`.
fn blob_at(buf: &[u8], off: usize) -> (i32, i32, i32) {
    (le_i32(buf, off), le_i32(buf, off + 4), le_i32(buf, off + 8))
}

// --- frequency analysis ----------------------------------------------

/// Quasi-periodic complex signal sum_k A_k exp(i (f_k j + phi_k)),
/// sampled at j = 0..ndata-1 and interleaved as x0,y0,x1,y1,...
fn tone_signal(ndata: usize, tones: &[(f64, f64, f64)]) -> Vec<f64> {
    let mut v = vec![0.0f64; 2 * ndata];
    for j in 0..ndata {
        let mut xr = 0.0f64;
        let mut yi = 0.0f64;
        for &(f, a, phi) in tones {
            let arg = f * (j as f64) + phi;
            xr += a * arg.cos();
            yi += a * arg.sin();
        }
        v[2 * j] = xr;
        v[2 * j + 1] = yi;
    }
    v
}

/// Wrap an angle difference into (-pi, pi] so phases can be compared
/// across the 0/2pi branch cut used by `reb_frequency_analysis`.
fn angle_diff(a: f64, b: f64) -> f64 {
    let mut d = a - b;
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d <= -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

// =====================================================================
// Frequency analysis: FFT / MFT / FMFT
// =====================================================================

#[test]
fn freq_analysis_argument_errors() {
    let input = tone_signal(64, &[(0.3, 1.0, 0.0)]);
    let mut out = vec![0.0f64; 3];

    // minfreq must be smaller than maxfreq  -> -1
    let e1 = reb_frequency_analysis(
        &mut out,
        1,
        1.0,
        1.0,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        64,
    );
    assert_eq!(e1, -1, "minfreq == maxfreq must return -1");

    // nfreq must be > 0 -> -2
    let e2 = reb_frequency_analysis(
        &mut out,
        0,
        -1.0,
        1.0,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        64,
    );
    assert_eq!(e2, -2, "nfreq == 0 must return -2");

    // ndata must be a power of two -> -3
    let input63 = tone_signal(63, &[(0.3, 1.0, 0.0)]);
    let e3 = reb_frequency_analysis(
        &mut out,
        1,
        -1.0,
        1.0,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input63,
        63,
    );
    assert_eq!(e3, -3, "non power-of-two ndata must return -3");
}

#[test]
fn mft_recovers_single_tone() {
    // A single pure tone: the MFT must return exactly the frequency,
    // amplitude and phase that built the signal.
    let ndata = 512usize;
    let f0 = 0.30;
    let a0 = 1.5;
    let p0 = 0.7;
    let input = tone_signal(ndata, &[(f0, a0, p0)]);
    let mut out = vec![0.0f64; 3];
    let err = reb_frequency_analysis(
        &mut out,
        1,
        0.05,
        1.0,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        ndata,
    );
    assert_eq!(err, 0, "MFT of a single tone must succeed");

    assert!(
        (out[0] - f0).abs() < 1e-9,
        "MFT frequency {} differs from the constructed {} by {}",
        out[0],
        f0,
        (out[0] - f0).abs()
    );
    assert!(
        (out[1] - a0).abs() < 1e-9,
        "MFT amplitude {} differs from the constructed {}",
        out[1],
        a0
    );
    // The phase is read off at sample 0 but the frequency error is
    // amplified by the record length, so it tracks ndata/2 * df.
    assert!(
        angle_diff(out[2], p0).abs() < 1e-6,
        "MFT phase {} differs from the constructed {}",
        out[2],
        p0
    );
}

#[test]
fn mft_single_tone_reconstruction_roundtrip() {
    // Independent route: rebuild the signal from the returned
    // (f, A, psi) and compare it against the original samples.
    let ndata = 512usize;
    let f0 = 0.8123;
    let a0 = 0.37;
    let p0 = -1.2;
    let input = tone_signal(ndata, &[(f0, a0, p0)]);
    let mut out = vec![0.0f64; 3];
    let err = reb_frequency_analysis(
        &mut out,
        1,
        0.1,
        1.5,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        ndata,
    );
    assert_eq!(err, 0, "MFT must succeed");

    let recon = tone_signal(ndata, &[(out[0], out[1], out[2])]);
    let mut worst = 0.0f64;
    for i in 0..2 * ndata {
        let d = (recon[i] - input[i]).abs();
        if d > worst {
            worst = d;
        }
    }
    assert!(
        worst < 1e-7,
        "reconstruction from the MFT output differs from the input by {}",
        worst
    );
}

#[test]
fn mft_recovers_negative_frequency() {
    // `bracket` has two branches (peak in the lower / upper half of the
    // DFT). A negative frequency exercises the other one.
    let ndata = 512usize;
    let f0 = -0.44;
    let a0 = 2.0;
    let p0 = 2.5;
    let input = tone_signal(ndata, &[(f0, a0, p0)]);
    let mut out = vec![0.0f64; 3];
    let err = reb_frequency_analysis(
        &mut out,
        1,
        -1.0,
        -0.05,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        ndata,
    );
    assert_eq!(err, 0, "MFT of a negative-frequency tone must succeed");
    assert!(
        (out[0] - f0).abs() < 1e-9,
        "MFT frequency {} differs from the constructed {}",
        out[0],
        f0
    );
    assert!(
        (out[1] - a0).abs() < 1e-9,
        "MFT amplitude {} differs from the constructed {}",
        out[1],
        a0
    );
    assert!(
        angle_diff(out[2], p0).abs() < 1e-6,
        "MFT phase {} differs from the constructed {}",
        out[2],
        p0
    );
}

#[test]
fn mft_two_tones_sorted_by_decreasing_amplitude() {
    // sort3 must return the strongest component first, regardless of the
    // order the frequencies were found in.
    let ndata = 1024usize;
    let tones = [(0.21, 0.4, 0.3), (0.63, 1.0, -0.8)];
    let input = tone_signal(ndata, &tones);
    let mut out = vec![0.0f64; 6];
    let err = reb_frequency_analysis(
        &mut out,
        2,
        0.05,
        1.2,
        REB_FREQUENCY_ANALYSIS_MFT,
        &input,
        ndata,
    );
    assert_eq!(err, 0, "MFT of a two-tone signal must succeed");

    assert!(
        out[2] >= out[3],
        "amplitudes must come back in decreasing order, got {} then {}",
        out[2],
        out[3]
    );
    // Strongest tone is the one with amplitude 1.0 at f = 0.63.
    assert!(
        (out[0] - 0.63).abs() < 1e-4,
        "leading frequency {} should be the strong tone at 0.63",
        out[0]
    );
    assert!(
        (out[1] - 0.21).abs() < 1e-4,
        "second frequency {} should be the weak tone at 0.21",
        out[1]
    );
    assert!(
        (out[2] - 1.0).abs() < 1e-3,
        "leading amplitude {} should be near 1.0",
        out[2]
    );
    assert!(
        (out[3] - 0.4).abs() < 1e-3,
        "second amplitude {} should be near 0.4",
        out[3]
    );
}

#[test]
fn fmft_is_more_accurate_than_mft_on_two_tones() {
    // The whole point of the FMFT is that one MFT pass leaves a residual
    // frequency error caused by the interference of the components, and
    // the extra pass extrapolates it away.
    let ndata = 1024usize;
    let tones = [(0.35, 1.0, 0.2), (0.42, 0.6, 1.1)];
    let input = tone_signal(ndata, &tones);

    let mut mft = vec![0.0f64; 6];
    assert_eq!(
        reb_frequency_analysis(
            &mut mft,
            2,
            0.05,
            1.0,
            REB_FREQUENCY_ANALYSIS_MFT,
            &input,
            ndata
        ),
        0,
        "MFT must succeed"
    );
    let mut fmft = vec![0.0f64; 6];
    assert_eq!(
        reb_frequency_analysis(
            &mut fmft,
            2,
            0.05,
            1.0,
            REB_FREQUENCY_ANALYSIS_FMFT,
            &input,
            ndata
        ),
        0,
        "FMFT must succeed"
    );

    // Strongest tone first in both cases: f = 0.35, A = 1.0.
    let mft_err = (mft[0] - 0.35).abs();
    let fmft_err = (fmft[0] - 0.35).abs();
    assert!(
        mft_err < 1e-4,
        "MFT leading frequency error {} is unexpectedly large",
        mft_err
    );
    assert!(
        fmft_err < mft_err,
        "FMFT frequency error {} should improve on the MFT error {}",
        fmft_err,
        mft_err
    );
    assert!(
        fmft_err < 1e-8,
        "FMFT leading frequency error {} is unexpectedly large",
        fmft_err
    );
}

#[test]
fn fmft2_recovers_two_tones() {
    let ndata = 1024usize;
    let tones = [(0.35, 1.0, 0.2), (0.42, 0.6, 1.1)];
    let input = tone_signal(ndata, &tones);
    let mut out = vec![0.0f64; 6];
    assert_eq!(
        reb_frequency_analysis(
            &mut out,
            2,
            0.05,
            1.0,
            REB_FREQUENCY_ANALYSIS_FMFT2,
            &input,
            ndata
        ),
        0,
        "FMFT2 must succeed"
    );
    assert!(
        (out[0] - 0.35).abs() < 1e-8,
        "FMFT2 leading frequency {} differs from 0.35",
        out[0]
    );
    assert!(
        (out[1] - 0.42).abs() < 1e-8,
        "FMFT2 second frequency {} differs from 0.42",
        out[1]
    );
    assert!(
        (out[2] - 1.0).abs() < 1e-6,
        "FMFT2 leading amplitude {} differs from 1.0",
        out[2]
    );
    assert!(
        (out[3] - 0.6).abs() < 1e-6,
        "FMFT2 second amplitude {} differs from 0.6",
        out[3]
    );
}

#[test]
fn freq_analysis_is_bitwise_deterministic() {
    let ndata = 512usize;
    let input = tone_signal(ndata, &[(0.31, 1.0, 0.1), (0.77, 0.5, -0.4)]);
    let mut a = vec![0.0f64; 6];
    let mut b = vec![0.0f64; 6];
    assert_eq!(
        reb_frequency_analysis(&mut a, 2, 0.05, 1.5, REB_FREQUENCY_ANALYSIS_FMFT, &input, ndata),
        0
    );
    assert_eq!(
        reb_frequency_analysis(&mut b, 2, 0.05, 1.5, REB_FREQUENCY_ANALYSIS_FMFT, &input, ndata),
        0
    );
    for k in 0..6 {
        assert_eq!(
            a[k].to_bits(),
            b[k].to_bits(),
            "output[{}] is not bit-identical between two identical FMFT calls",
            k
        );
    }
}

#[test]
fn mft_frequency_invariant_under_exact_amplitude_doubling() {
    // Scaling the signal by 2 scales every windowed sum, every power
    // spectral density value and every phisqr by an exact power of two,
    // so the golden-section search must take the identical branches and
    // return the identical frequency bits; the amplitude must scale by
    // exactly the same factor.
    let ndata = 512usize;
    let base = tone_signal(ndata, &[(0.29, 1.0, 0.55)]);
    let doubled: Vec<f64> = base.iter().map(|v| v * 2.0).collect();

    let mut o1 = vec![0.0f64; 3];
    let mut o2 = vec![0.0f64; 3];
    assert_eq!(
        reb_frequency_analysis(&mut o1, 1, 0.05, 1.0, REB_FREQUENCY_ANALYSIS_MFT, &base, ndata),
        0
    );
    assert_eq!(
        reb_frequency_analysis(
            &mut o2,
            1,
            0.05,
            1.0,
            REB_FREQUENCY_ANALYSIS_MFT,
            &doubled,
            ndata
        ),
        0
    );
    assert_eq!(
        o1[0].to_bits(),
        o2[0].to_bits(),
        "frequency must be bit-identical under an exact factor-of-two rescaling ({} vs {})",
        o1[0],
        o2[0]
    );
    assert!(
        (o2[1] - 2.0 * o1[1]).abs() <= 1e-14 * o2[1].abs(),
        "amplitude {} should be exactly twice {}",
        o2[1],
        o1[1]
    );
}

#[test]
fn mft_phase_is_normalized_into_zero_two_pi() {
    // The C normalizes the returned phases into [0, 2pi) before sorting.
    let ndata = 512usize;
    let two_pi = 2.0 * std::f64::consts::PI;
    for &p0 in &[-3.0f64, -1.0, 0.25, 2.0, 3.0] {
        let input = tone_signal(ndata, &[(0.37, 1.0, p0)]);
        let mut out = vec![0.0f64; 3];
        assert_eq!(
            reb_frequency_analysis(
                &mut out,
                1,
                0.05,
                1.0,
                REB_FREQUENCY_ANALYSIS_MFT,
                &input,
                ndata
            ),
            0
        );
        assert!(
            out[2] >= 0.0 && out[2] < two_pi,
            "phase {} (input phase {}) is outside [0, 2pi)",
            out[2],
            p0
        );
        assert!(
            angle_diff(out[2], p0).abs() < 1e-6,
            "normalized phase {} does not match input phase {} modulo 2pi",
            out[2],
            p0
        );
    }
}

// =====================================================================
// Binary format: header, field order, byte layout
// =====================================================================

#[test]
fn binary_header_layout() {
    let mut r = deterministic_sim();
    let buf = reb_binarydata_simulation_to_stream(&mut r);
    assert!(buf.len() > 64, "buffer must contain a 64-byte header");

    // The C writes sprintf(header, "REBOUND Binary File. Version: %s"),
    // then snprintf(header+cwritten+1, 64-cwritten-1, "%s", githash).
    let vstr = format!("REBOUND Binary File. Version: {}", reb_version_str);
    let cwritten = vstr.len();
    assert_eq!(
        &buf[..cwritten],
        vstr.as_bytes(),
        "header version string mismatch"
    );
    assert_eq!(buf[cwritten], 0, "version string must be NUL terminated");

    let gstart = cwritten + 1;
    let gcap = 64 - gstart - 1; // snprintf reserves the trailing NUL
    let g = reb_githash_str.as_bytes();
    let glen = std::cmp::min(g.len(), gcap);
    assert_eq!(
        &buf[gstart..gstart + glen],
        &g[..glen],
        "githash must follow the NUL of the version string"
    );
    assert_eq!(buf[63], 0, "the 64-byte header must end with a NUL");

    // The first eight header bytes spell "REBOUND " and are what
    // `reb_binarydata_header` is: the reader uses them to tell a header
    // apart from a field descriptor.
    assert_eq!(
        le_u64(&buf, 0),
        reb_binarydata_header,
        "first 8 header bytes must equal reb_binarydata_header"
    );
    assert_eq!(&buf[..8], b"REBOUND ", "reb_binarydata_header spells 'REBOUND '");
}

#[test]
fn binary_field_order_matches_c_descriptor_list() {
    // Order taken from reb_binarydata_field_descriptor_list in the C
    // binarydata.c, followed by the whfast descriptor list from
    // integrator_whfast.c, then "functionpointers" and "end".
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "whfast");
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(
        &mut r,
        reb_particle {
            m: 1e-3,
            x: 1.0,
            vy: 1.0,
            ..reb_particle::default()
        },
    );
    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let got = field_names(&buf);

    let expected: Vec<&str> = vec![
        "t",
        "G",
        "softening",
        "dt",
        "N",
        "N_var",
        "N_active",
        "testparticle_type",
        "opening_angle2",
        "status",
        "exact_finish_time",
        "force_is_velocity_dependent",
        "gravity_ignore_terms",
        "output_timing_last",
        "save_messages",
        "exit_max_distance",
        "exit_min_distance",
        "usleep",
        "track_energy_offset",
        "energy_offset",
        "root_size",
        "N_root_x",
        "N_root_y",
        "N_root_z",
        "N_ghost_x",
        "N_ghost_y",
        "N_ghost_z",
        "minimum_collision_velocity",
        "collisions_plog",
        "collisions_log_n",
        "calculate_megno",
        "megno_Ys",
        "megno_Yss",
        "megno_cov_Yt",
        "megno_var_t",
        "megno_mean_t",
        "megno_mean_Y",
        "megno_initial_t",
        "megno_n",
        "simulationarchive_auto_interval",
        "simulationarchive_auto_walltime",
        "simulationarchive_next",
        "collision",
        "integrator.name",
        "boundary",
        "gravity",
        "OMEGA",
        "OMEGAZ",
        "is_synchronized",
        "did_modify_particles",
        // particles_var / var_config / display_settings / name_list are
        // REB_POINTER fields with size_data == 0 here and are skipped.
        "particles",
        "simulationarchive_version",
        "walltime",
        "walltime_last_steps",
        "python_unit_l",
        "python_unit_m",
        "python_unit_t",
        "simulationarchive_auto_step",
        "simulationarchive_next_step",
        "steps_done",
        "dt_last_done",
        "rand_seed",
        "testparticle_hidewarnings",
        // whfast descriptor list
        "integrator.whfast.corrector",
        "integrator.whfast.safe_mode",
        "integrator.whfast.coordinates",
        "integrator.whfast.corrector2",
        "integrator.whfast.kernel",
        "integrator.whfast.keep_unsynchronized",
        "functionpointers",
        "end",
    ];
    assert_eq!(
        got.len(),
        expected.len(),
        "field count mismatch: got {:?}",
        got
    );
    for (i, e) in expected.iter().enumerate() {
        assert_eq!(
            &got[i], e,
            "field {} should be {} but is {}",
            i, e, got[i]
        );
    }
}

#[test]
fn binary_field_header_sizes_are_exact() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 2.0, x: 0.5, ..reb_particle::default() });
    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let fields = parse_fields(&buf);

    for f in &fields {
        assert_eq!(
            f.size_name,
            (f.name.len() + 1) as u64,
            "size_name of {} must be strlen+1",
            f.name
        );
    }

    // Widths follow the C descriptor types.
    for name in ["t", "G", "softening", "dt", "OMEGA", "OMEGAZ", "energy_offset"] {
        assert_eq!(
            find_field(&fields, name).size_data,
            8,
            "REB_DOUBLE field {} must be 8 bytes",
            name
        );
    }
    for name in ["N", "N_var", "N_active", "N_root_x", "N_root_y", "N_root_z"] {
        assert_eq!(
            find_field(&fields, name).size_data,
            8,
            "REB_SIZE_T field {} must be 8 bytes",
            name
        );
    }
    for name in [
        "status",
        "exact_finish_time",
        "save_messages",
        "N_ghost_x",
        "collision",
        "boundary",
        "gravity",
        "calculate_megno",
        "simulationarchive_version",
        "testparticle_type",
        "testparticle_hidewarnings",
        "functionpointers",
    ] {
        assert_eq!(
            find_field(&fields, name).size_data,
            4,
            "REB_INT field {} must be 4 bytes",
            name
        );
    }
    for name in ["collisions_log_n", "megno_n"] {
        assert_eq!(
            find_field(&fields, name).size_data,
            8,
            "REB_INT64 field {} must be 8 bytes",
            name
        );
    }
    for name in ["steps_done", "simulationarchive_auto_step", "simulationarchive_next_step"] {
        assert_eq!(
            find_field(&fields, name).size_data,
            8,
            "REB_UINT64 field {} must be 8 bytes",
            name
        );
    }
    for name in ["is_synchronized", "did_modify_particles", "rand_seed", "python_unit_l"] {
        assert_eq!(
            find_field(&fields, name).size_data,
            4,
            "REB_UINT field {} must be 4 bytes",
            name
        );
    }
    // "integrator.name" is a REB_STRING: strlen+1 of the integrator name.
    let iname = find_field(&fields, "integrator.name");
    assert_eq!(
        iname.size_data,
        (r.integrator.name().len() + 1) as u64,
        "integrator.name payload must be the name plus NUL"
    );
    assert_eq!(
        &buf[iname.data_pos..iname.data_pos + iname.size_data as usize - 1],
        r.integrator.name().as_bytes()
    );
    assert_eq!(buf[iname.data_pos + iname.size_data as usize - 1], 0);

    // Scalar values must be the simulation's own bits.
    assert_eq!(le_f64(&buf, find_field(&fields, "G").data_pos).to_bits(), r.G.to_bits());
    assert_eq!(le_f64(&buf, find_field(&fields, "dt").data_pos).to_bits(), r.dt.to_bits());
    assert_eq!(le_u64(&buf, find_field(&fields, "N").data_pos), r.N as u64);
    assert_eq!(
        le_u64(&buf, find_field(&fields, "N_active").data_pos),
        usize::MAX as u64,
        "the default N_active sentinel is SIZE_MAX"
    );
}

#[test]
fn binary_particles_payload_layout() {
    // sizeof(struct reb_particle) on x86-64 is 112: eleven doubles then
    // the name / ap / sim pointers.
    let mut r = deterministic_sim();
    let p0 = reb_particle {
        x: 1.5,
        y: -2.25,
        z: 0.125,
        vx: 3.0,
        vy: -4.5,
        vz: 0.0625,
        ax: 7.0,
        ay: 8.0,
        az: 9.0,
        m: 1.0,
        r: 0.5,
        name: None,
    };
    reb_simulation_add(&mut r, p0);
    reb_particle_set_name(&mut r, 0, Some("star"));
    reb_simulation_add(&mut r, reb_particle { m: 2.0, ..reb_particle::default() });

    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let fields = parse_fields(&buf);
    let pf = find_field(&fields, "particles");
    assert_eq!(
        pf.size_data,
        (2 * REB_PARTICLE_RAW_SIZE) as u64,
        "particles payload must be N * 112 bytes"
    );

    let base = pf.data_pos;
    let want = [
        r.particles[0].x,
        r.particles[0].y,
        r.particles[0].z,
        r.particles[0].vx,
        r.particles[0].vy,
        r.particles[0].vz,
        r.particles[0].ax,
        r.particles[0].ay,
        r.particles[0].az,
        r.particles[0].m,
        r.particles[0].r,
    ];
    for (i, w) in want.iter().enumerate() {
        assert_eq!(
            le_f64(&buf, base + 8 * i).to_bits(),
            w.to_bits(),
            "particle double at offset {} is wrong",
            8 * i
        );
    }
    // name pointer slot (offset 88): non-zero for a named particle, and
    // it must match the pointer recorded in the name_list blob.
    let name_ptr = le_u64(&buf, base + 88);
    assert_ne!(name_ptr, 0, "named particle must store a non-zero name slot");
    assert_eq!(le_u64(&buf, base + 96), 0, "ap pointer slot must be zero");
    assert_eq!(le_u64(&buf, base + 104), 0, "sim pointer slot must be zero");

    // Second particle is unnamed -> zero name slot.
    let base2 = base + REB_PARTICLE_RAW_SIZE;
    assert_eq!(
        le_u64(&buf, base2 + 88),
        0,
        "unnamed particle must store a zero name slot"
    );

    // The name_list blob is "star\0" followed by the same pointer value.
    let nl = find_field(&fields, "name_list");
    assert_eq!(
        nl.size_data,
        (4 + 1 + 8) as u64,
        "name_list entry is strlen+1+sizeof(char*)"
    );
    assert_eq!(&buf[nl.data_pos..nl.data_pos + 4], b"star");
    assert_eq!(buf[nl.data_pos + 4], 0);
    assert_eq!(
        le_u64(&buf, nl.data_pos + 5),
        name_ptr,
        "name_list pointer must equal the particle's name slot"
    );
}

#[test]
fn binary_trailer_is_end_field_plus_zero_blob() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    let buf = reb_binarydata_simulation_to_stream(&mut r);

    // ... [16-byte field header][ "end\0" ][12-byte zeroed blob]
    let n = buf.len();
    assert!(n > 16 + 4 + REB_SA_BLOB_SIZE);
    assert_eq!(
        &buf[n - REB_SA_BLOB_SIZE..],
        &[0u8; REB_SA_BLOB_SIZE],
        "a fresh binary ends with an all-zero simulationarchive blob"
    );
    let end_name = n - REB_SA_BLOB_SIZE - 4;
    assert_eq!(&buf[end_name..end_name + 4], b"end\0");
    let end_hdr = end_name - REB_BINARYDATA_FIELD_SIZE;
    assert_eq!(le_u64(&buf, end_hdr), 4, "size_name of 'end' is 4");
    assert_eq!(le_u64(&buf, end_hdr + 8), 0, "size_data of 'end' is 0");
}

#[test]
fn binary_roundtrip_is_bit_exact() {
    let mut r = two_body("whfast", 0.01, 0.2);
    reb_simulation_steps(&mut r, 37);
    let before = all_particle_bits(&r);
    let t_before = r.t.to_bits();
    let dt_before = r.dt.to_bits();
    let steps_before = r.steps_done;

    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let mut loaded = reb_simulation_create();
    loaded.save_messages = 1;
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let mut cur = Cursor::new(buf);
    reb_binarydata_input_fields(&mut loaded, &mut cur, &mut warnings);
    assert_eq!(
        warnings & REB_BINARYDATA_WARNING_CORRUPTFILE,
        0,
        "round-tripped buffer must not look corrupt"
    );
    assert_eq!(
        warnings & REB_BINARYDATA_WARNING_VERSION,
        0,
        "round-tripped buffer must not report a version mismatch"
    );
    assert_eq!(
        warnings & REB_BINARYDATA_WARNING_FIELD_UNKNOWN,
        0,
        "every written field must be understood by the reader"
    );

    assert_eq!(loaded.N, r.N, "N must survive the round trip");
    assert_eq!(loaded.t.to_bits(), t_before, "t must survive bit-exactly");
    assert_eq!(loaded.dt.to_bits(), dt_before, "dt must survive bit-exactly");
    assert_eq!(loaded.steps_done, steps_before);
    assert_eq!(loaded.G.to_bits(), r.G.to_bits());
    assert_eq!(loaded.integrator.name(), "whfast");
    let after = all_particle_bits(&loaded);
    assert_eq!(after, before, "particle bits must survive the round trip");
}

#[test]
fn binary_roundtrip_preserves_particle_names() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 2.0, x: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 3.0, x: 2.0, ..reb_particle::default() });
    reb_particle_set_name(&mut r, 0, Some("sun"));
    reb_particle_set_name(&mut r, 2, Some("outer"));

    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let mut loaded = reb_simulation_create();
    loaded.save_messages = 1;
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let mut cur = Cursor::new(buf);
    reb_binarydata_input_fields(&mut loaded, &mut cur, &mut warnings);

    assert_eq!(
        reb_simulation_get_particle_by_name(&loaded, "sun"),
        Some(0),
        "'sun' must resolve back to particle 0"
    );
    assert_eq!(
        reb_simulation_get_particle_by_name(&loaded, "outer"),
        Some(2),
        "'outer' must resolve back to particle 2"
    );
    assert_eq!(
        loaded.particles[1].name, None,
        "the unnamed particle must stay unnamed"
    );
    assert_eq!(
        reb_simulation_get_particle_by_name(&loaded, "nobody"),
        None,
        "an unregistered name must not resolve"
    );
}

#[test]
fn binary_roundtrip_preserves_integrator_state() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "whfast");
    if let reb_integrator_state::whfast(ref mut wh) = r.integrator {
        wh.corrector = 17;
        wh.corrector2 = 1;
        wh.kernel = 3;
        wh.coordinates = 1;
        wh.safe_mode = 0;
        wh.keep_unsynchronized = 1;
    }
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let mut loaded = reb_simulation_create();
    loaded.save_messages = 1;
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let mut cur = Cursor::new(buf);
    reb_binarydata_input_fields(&mut loaded, &mut cur, &mut warnings);

    match loaded.integrator {
        reb_integrator_state::whfast(ref wh) => {
            assert_eq!(wh.corrector, 17, "corrector must round trip");
            assert_eq!(wh.corrector2, 1, "corrector2 must round trip");
            assert_eq!(wh.kernel, 3, "kernel must round trip");
            assert_eq!(wh.coordinates, 1, "coordinates must round trip");
            assert_eq!(wh.safe_mode, 0, "safe_mode must round trip");
            assert_eq!(wh.keep_unsynchronized, 1, "keep_unsynchronized must round trip");
        }
        ref other => panic!("integrator must be whfast after load, got {}", other.name()),
    }
}

#[test]
fn binary_diff_of_identical_buffers_is_empty() {
    let mut r = two_body("ias15", 0.01, 0.1);
    let b1 = reb_binarydata_simulation_to_stream(&mut r);
    let b2 = reb_binarydata_simulation_to_stream(&mut r);
    assert_eq!(b1, b2, "serializing the same simulation twice must be byte-identical");

    let (differ, stream) = reb_binarydata_diff(&b1, &b2, REB_BINARYDATA_OUTPUT_STREAM);
    assert_eq!(differ, 0, "identical buffers must not be reported different");
    assert!(stream.is_empty(), "identical buffers must produce an empty diff stream");
}

#[test]
fn binary_diff_single_scalar_has_exact_stream_bytes() {
    // Changing only `t` must produce a diff stream of exactly one field:
    // 8 bytes size_name (2) + 8 bytes size_data (8) + "t\0" + the double.
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    let b1 = reb_binarydata_simulation_to_stream(&mut r);
    r.t = 5.0;
    let b2 = reb_binarydata_simulation_to_stream(&mut r);

    let (differ, stream) = reb_binarydata_diff(&b1, &b2, REB_BINARYDATA_OUTPUT_STREAM);
    assert_eq!(differ, 1, "a changed t must be reported as different");
    assert_eq!(stream.len(), 26, "diff stream must be 16 + 2 + 8 bytes, got {:?}", stream);
    assert_eq!(le_u64(&stream, 0), 2, "size_name of 't' is 2");
    assert_eq!(le_u64(&stream, 8), 8, "size_data of 't' is 8");
    assert_eq!(&stream[16..18], b"t\0");
    assert_eq!(le_f64(&stream, 18).to_bits(), 5.0f64.to_bits());
}

#[test]
fn binary_diff_ignores_walltime_for_the_return_value() {
    // The C deliberately excludes every field whose name starts with
    // "walltime" from `are_different`, but still writes it to the stream.
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    let b1 = reb_binarydata_simulation_to_stream(&mut r);
    r.walltime = 1.25;
    let b2 = reb_binarydata_simulation_to_stream(&mut r);
    assert_ne!(b1, b2, "changing walltime must change the serialized bytes");

    let (differ, stream) = reb_binarydata_diff(&b1, &b2, REB_BINARYDATA_OUTPUT_STREAM);
    assert_eq!(differ, 0, "a walltime-only change must not count as different");
    assert_eq!(
        stream.len(),
        16 + 9 + 8,
        "the stream must still carry the walltime field"
    );
    assert_eq!(&stream[16..25], b"walltime\0");
    assert_eq!(le_f64(&stream, 25).to_bits(), 1.25f64.to_bits());
}

#[test]
fn binary_diff_detects_particle_changes() {
    let mut r = two_body("leapfrog", 0.01, 0.0);
    let b1 = reb_binarydata_simulation_to_stream(&mut r);
    r.particles[1].x += 1e-12;
    let b2 = reb_binarydata_simulation_to_stream(&mut r);
    let (differ, _stream) = reb_binarydata_diff(&b1, &b2, REB_BINARYDATA_OUTPUT_NONE);
    assert_eq!(differ, 1, "a moved particle must be reported as different");

    // The C compares particles with reb_particle_cmp, which ignores the
    // ap and sim pointer slots but not the 11 doubles.
    let mut r2 = two_body("leapfrog", 0.01, 0.0);
    let b3 = reb_binarydata_simulation_to_stream(&mut r2);
    let (same, _s) = reb_binarydata_diff(&b1, &b3, REB_BINARYDATA_OUTPUT_NONE);
    assert_eq!(same, 0, "two identically constructed simulations must not differ");
}

#[test]
fn binary_input_flags_version_mismatch() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    let mut buf = reb_binarydata_simulation_to_stream(&mut r);
    // "REBOUND Binary File. Version: " is 30 bytes; corrupt the first
    // character of the version itself.
    let vpos = "REBOUND Binary File. Version: ".len();
    assert_eq!(buf[vpos], reb_version_str.as_bytes()[0]);
    buf[vpos] = b'0';

    let mut loaded = reb_simulation_create();
    loaded.save_messages = 1;
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let mut cur = Cursor::new(buf);
    reb_binarydata_input_fields(&mut loaded, &mut cur, &mut warnings);
    assert_ne!(
        warnings & REB_BINARYDATA_WARNING_VERSION,
        0,
        "a different version string must raise REB_BINARYDATA_WARNING_VERSION"
    );
    // The remaining fields are still read.
    assert_eq!(loaded.N, 1, "fields after the header must still be read");

    assert_eq!(
        reb_binarydata_process_warnings(&mut loaded, warnings),
        0,
        "a version warning is not fatal"
    );
    assert!(
        loaded.messages.iter().any(|(t, m)| *t == REB_MESSAGE_TYPE::WARNING
            && m.contains("different version of REBOUND")),
        "process_warnings must record the version warning, got {:?}",
        loaded.messages
    );
}

#[test]
fn binary_input_flags_unknown_field() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    let good = reb_binarydata_simulation_to_stream(&mut r);

    // Splice an unknown field in right after the header.
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(&good[..64]);
    let name = b"not_a_real_field\0";
    buf.extend_from_slice(&(name.len() as u64).to_le_bytes());
    buf.extend_from_slice(&8u64.to_le_bytes());
    buf.extend_from_slice(name);
    buf.extend_from_slice(&0u64.to_le_bytes());
    buf.extend_from_slice(&good[64..]);

    let mut loaded = reb_simulation_create();
    loaded.save_messages = 1;
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let mut cur = Cursor::new(buf);
    reb_binarydata_input_fields(&mut loaded, &mut cur, &mut warnings);
    assert_ne!(
        warnings & REB_BINARYDATA_WARNING_FIELD_UNKNOWN,
        0,
        "an unrecognized field name must raise REB_BINARYDATA_WARNING_FIELD_UNKNOWN"
    );
    assert_eq!(
        reb_binarydata_process_warnings(&mut loaded, warnings),
        0,
        "an unknown field is a warning, not a fatal error"
    );
}

#[test]
fn binary_process_warnings_reports_fatal_codes() {
    let mut r = deterministic_sim();
    for code in [
        REB_BINARYDATA_ERROR_NOFILE,
        REB_BINARYDATA_ERROR_FILENOTOPEN,
        REB_BINARYDATA_ERROR_OUTOFRANGE,
        REB_BINARYDATA_ERROR_SEEK,
        REB_BINARYDATA_ERROR_OLD,
    ] {
        assert_eq!(
            reb_binarydata_process_warnings(&mut r, code),
            -1,
            "error code {} must be fatal",
            code
        );
    }
    for code in [
        REB_BINARYDATA_WARNING_NONE,
        REB_BINARYDATA_WARNING_VERSION,
        REB_BINARYDATA_WARNING_POINTERS,
        REB_BINARYDATA_WARNING_PARTICLES,
        REB_BINARYDATA_WARNING_FIELD_UNKNOWN,
        REB_BINARYDATA_WARNING_CORRUPTFILE,
        REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR,
    ] {
        assert_eq!(
            reb_binarydata_process_warnings(&mut r, code),
            0,
            "warning code {} must not be fatal",
            code
        );
    }
}

// =====================================================================
// Simulationarchive: blob chain and snapshot index
// =====================================================================

#[test]
fn archive_first_save_writes_a_full_binary() {
    let p = temp_path("first_save");
    let f = path_str(&p);
    let mut r = two_body("whfast", 0.01, 0.1);
    reb_simulation_save_to_file(&mut r, Some(&f));

    let on_disk = std::fs::read(&p).expect("archive file must exist after the first save");
    let expected = reb_binarydata_simulation_to_stream(&mut r);
    assert_eq!(
        on_disk, expected,
        "the first save must be exactly the full binary serialization"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_blob_chain_offsets_are_consistent() {
    let p = temp_path("blob_chain");
    let f = path_str(&p);
    let mut r = two_body("whfast", 0.01, 0.1);
    reb_simulation_save_to_file(&mut r, Some(&f));
    let len0 = std::fs::metadata(&p).unwrap().len() as usize;

    reb_simulation_steps(&mut r, 5);
    reb_simulation_save_to_file(&mut r, Some(&f));
    let bytes = std::fs::read(&p).unwrap();
    let len1 = bytes.len();
    assert!(len1 > len0, "appending a snapshot must grow the file");

    // The first snapshot's trailing blob sits at len0-12 and now points
    // at the appended snapshot; the new trailing blob closes the chain.
    let (i0, prev0, next0) = blob_at(&bytes, len0 - REB_SA_BLOB_SIZE);
    let (i1, prev1, next1) = blob_at(&bytes, len1 - REB_SA_BLOB_SIZE);
    assert_eq!(i0, 0, "the initial blob has index 0");
    assert_eq!(prev0, 0, "the initial blob has no predecessor");
    assert_eq!(
        next0 as usize,
        len1 - len0 - REB_SA_BLOB_SIZE,
        "offset_next must be the appended snapshot's size without its own blob"
    );
    assert_eq!(i1, 1, "the appended blob has index 1");
    assert_eq!(prev1, next0, "offset_prev of blob 1 must mirror offset_next of blob 0");
    assert_eq!(next1, 0, "the last blob terminates the chain with offset_next == 0");

    // The reader's checksum: offset_prev + sizeof(blob) == distance from
    // the start of the snapshot to just past its blob.
    let mut warnings: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let sa = reb_simulationarchive_create_from_file_with_messages(&f, &mut warnings);
    assert_eq!(sa.nblobs, 2, "the archive must index two snapshots");
    assert_eq!(sa.offset[0], 0, "snapshot 0 starts at the beginning of the file");
    assert_eq!(
        sa.offset[1] as usize, len0,
        "snapshot 1 starts right after the first snapshot's blob"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_index_records_snapshot_times() {
    let p = temp_path("index_times");
    let f = path_str(&p);
    let mut r = deterministic_sim();
    r.dt = 1.0;
    reb_simulation_set_integrator(&mut r, "none");
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let mut expected_t: Vec<f64> = Vec::new();
    for _ in 0..4 {
        reb_simulation_save_to_file(&mut r, Some(&f));
        expected_t.push(r.t);
        reb_simulation_steps(&mut r, 3);
    }

    let sa = reb_simulationarchive_create_from_file(&f).expect("archive must open");
    assert_eq!(sa.nblobs, 4, "four saves must produce four snapshots");
    assert_eq!(sa.version, 5, "the archive version field must be 5");
    // reb_version_major is checked separately: see
    // archive_header_reports_the_major_rebound_version.
    assert_eq!(
        sa.reb_version_minor, 1,
        "the minor version is parsed from the header string"
    );
    assert_eq!(
        sa.reb_version_patch, 1,
        "the patch version is parsed from the header string"
    );
    for i in 0..4 {
        assert_eq!(
            sa.t[i].to_bits(),
            expected_t[i].to_bits(),
            "snapshot {} time must be the simulation time at save",
            i
        );
    }
    // "none" advances t by exactly dt, so the times are exact multiples.
    assert_eq!(sa.t, vec![0.0, 3.0, 6.0, 9.0]);
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_header_reports_the_major_rebound_version() {
    // simulationarchive.c splits the header at ':' and the two '.'
    // separators, so cmajor is the substring " 5" (note the leading
    // space that follows the colon in "Version: 5.1.1"). C's atoi()
    // skips leading whitespace and yields 5. The Rust translation's
    // atoi closure collects a leading run of ASCII digits instead, hits
    // the space first and yields 0, so reb_version_major is always 0
    // while reb_version_minor/patch (which have no leading space) are
    // parsed correctly. Fix belongs in
    // simulationarchive.rs::reb_simulationarchive_read_from_stream_with_messages.
    let p = temp_path("version_major");
    let f = path_str(&p);
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_save_to_file(&mut r, Some(&f));

    let sa = reb_simulationarchive_create_from_file(&f).expect("archive must open");
    let major: i32 = reb_version_str
        .split('.')
        .next()
        .unwrap()
        .parse()
        .expect("the crate version string starts with the major number");
    assert_eq!(
        sa.reb_version_major, major,
        "the archive must report the major version from its own header"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_snapshots_reload_bit_exactly() {
    let p = temp_path("reload");
    let f = path_str(&p);
    let mut r = two_body("whfast", 0.005, 0.15);

    reb_simulation_save_to_file(&mut r, Some(&f));
    let bits0 = all_particle_bits(&r);
    let t0 = r.t.to_bits();

    reb_simulation_steps(&mut r, 40);
    reb_simulation_save_to_file(&mut r, Some(&f));
    let bits1 = all_particle_bits(&r);
    let t1 = r.t.to_bits();
    assert_ne!(bits0, bits1, "the simulation must actually have moved");

    let s0 = reb_simulation_create_from_file(&f, 0).expect("snapshot 0 must load");
    assert_eq!(s0.t.to_bits(), t0, "snapshot 0 time must be bit-exact");
    assert_eq!(all_particle_bits(&s0), bits0, "snapshot 0 particles must be bit-exact");

    let s1 = reb_simulation_create_from_file(&f, 1).expect("snapshot 1 must load");
    assert_eq!(s1.t.to_bits(), t1, "snapshot 1 time must be bit-exact");
    assert_eq!(all_particle_bits(&s1), bits1, "snapshot 1 particles must be bit-exact");
    assert_eq!(s1.steps_done, r.steps_done, "steps_done must be carried in the diff");

    // Negative indices count from the end, like the C.
    let sm1 = reb_simulation_create_from_file(&f, -1).expect("snapshot -1 must load");
    assert_eq!(all_particle_bits(&sm1), bits1, "snapshot -1 is the last snapshot");
    let sm2 = reb_simulation_create_from_file(&f, -2).expect("snapshot -2 must load");
    assert_eq!(all_particle_bits(&sm2), bits0, "snapshot -2 is the first snapshot");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_rejects_out_of_range_snapshots() {
    let p = temp_path("out_of_range");
    let f = path_str(&p);
    let mut r = two_body("leapfrog", 0.01, 0.0);
    reb_simulation_save_to_file(&mut r, Some(&f));

    assert!(
        reb_simulation_create_from_file(&f, 1).is_none(),
        "snapshot 1 does not exist in a one-snapshot archive"
    );
    assert!(
        reb_simulation_create_from_file(&f, -2).is_none(),
        "snapshot -2 does not exist in a one-snapshot archive"
    );
    assert!(
        reb_simulation_create_from_file(&f, 0).is_some(),
        "snapshot 0 does exist"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_missing_file_returns_none() {
    let p = temp_path("does_not_exist");
    let f = path_str(&p);
    assert!(
        reb_simulation_create_from_file(&f, 0).is_none(),
        "a missing archive must yield None, not a default simulation"
    );
    assert!(
        reb_simulationarchive_create_from_file(&f).is_none(),
        "opening a missing archive must yield None"
    );
}

#[test]
fn archive_auto_interval_snapshots_at_expected_times() {
    // dt = 1 and the "none" integrator make every time exact, so the
    // snapshot times are fully determined: the heartbeat runs before
    // each step for t = 0..19 and fires whenever next <= t.
    let p = temp_path("auto_interval");
    let f = path_str(&p);
    let mut r = deterministic_sim();
    r.dt = 1.0;
    reb_simulation_set_integrator(&mut r, "none");
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_save_to_file_interval(&mut r, &f, 5.0);
    assert_eq!(r.simulationarchive_next.to_bits(), 0.0f64.to_bits());

    let status = reb_simulation_integrate(&mut r, 20.0);
    assert_eq!(status, REB_STATUS_SUCCESS);
    assert_eq!(r.t.to_bits(), 20.0f64.to_bits(), "t must land exactly on tmax");
    assert_eq!(r.steps_done, 20, "20 steps of dt = 1 reach t = 20");

    let sa = reb_simulationarchive_create_from_file(&f).expect("archive must open");
    assert_eq!(sa.nblobs, 4, "snapshots are taken at t = 0, 5, 10 and 15");
    assert_eq!(sa.t, vec![0.0, 5.0, 10.0, 15.0]);
    assert_eq!(
        sa.auto_interval.to_bits(),
        5.0f64.to_bits(),
        "the archive records the interval it was created with"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_auto_step_snapshots_at_expected_steps() {
    let p = temp_path("auto_step");
    let f = path_str(&p);
    let mut r = deterministic_sim();
    r.dt = 1.0;
    reb_simulation_set_integrator(&mut r, "none");
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_save_to_file_step(&mut r, &f, 5);
    assert_eq!(r.simulationarchive_next_step, 0);

    let status = reb_simulation_integrate(&mut r, 20.0);
    assert_eq!(status, REB_STATUS_SUCCESS);

    let sa = reb_simulationarchive_create_from_file(&f).expect("archive must open");
    assert_eq!(sa.nblobs, 4, "snapshots at steps_done = 0, 5, 10 and 15");
    assert_eq!(sa.t, vec![0.0, 5.0, 10.0, 15.0]);
    assert_eq!(sa.auto_step, 5, "the archive records the step interval");
    assert_eq!(
        r.simulationarchive_next_step, 20,
        "next_step advances past the last snapshot"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn archive_from_buffer_matches_archive_from_file() {
    let p = temp_path("from_buffer");
    let f = path_str(&p);
    let mut r = two_body("whfast", 0.01, 0.1);
    reb_simulation_save_to_file(&mut r, Some(&f));
    reb_simulation_steps(&mut r, 7);
    reb_simulation_save_to_file(&mut r, Some(&f));

    let bytes = std::fs::read(&p).unwrap();
    let mut w1: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let sa_file = reb_simulationarchive_create_from_file_with_messages(&f, &mut w1);
    let mut w2: REB_BINARYDATA_ERROR_CODE = REB_BINARYDATA_WARNING_NONE;
    let sa_buf = reb_simulationarchive_create_from_buffer_with_messages(&bytes, &mut w2);

    assert_eq!(w1, w2, "file and buffer readers must agree on warnings");
    assert_eq!(sa_file.nblobs, sa_buf.nblobs, "both must find the same snapshots");
    assert_eq!(sa_file.offset, sa_buf.offset, "both must find the same offsets");
    assert_eq!(sa_file.t, sa_buf.t, "both must find the same times");
    assert_eq!(sa_file.nblobs, 2);
    let _ = std::fs::remove_file(&p);
}

// =====================================================================
// Particle add / remove
// =====================================================================

#[test]
fn particle_add_appends_and_flags_modification() {
    let mut r = deterministic_sim();
    assert_eq!(r.N, 0, "a fresh simulation has no particles");
    assert_eq!(r.did_modify_particles, 0);

    for i in 0..5 {
        reb_simulation_add(
            &mut r,
            reb_particle {
                m: (i + 1) as f64,
                x: i as f64,
                ..reb_particle::default()
            },
        );
        assert_eq!(r.N, i + 1, "N must track the number of added particles");
        assert_eq!(r.particles.len(), r.N, "the array length must equal N");
    }
    assert_eq!(r.did_modify_particles, 1, "adding must taint the particle array");
    for i in 0..5 {
        assert_eq!(r.particles[i].m.to_bits(), ((i + 1) as f64).to_bits());
        assert_eq!(r.particles[i].x.to_bits(), (i as f64).to_bits());
    }
}

#[test]
fn particle_add_outside_the_box_is_refused() {
    let mut r = deterministic_sim();
    r.boundary = REB_BOUNDARY::PERIODIC;
    r.root_size = 2.0; // box spans [-1, 1] in each direction
    r.N_root_x = 1;
    r.N_root_y = 1;
    r.N_root_z = 1;

    reb_simulation_add(&mut r, reb_particle { m: 1.0, x: 0.5, ..reb_particle::default() });
    assert_eq!(r.N, 1, "a particle inside the box is accepted");

    reb_simulation_add(&mut r, reb_particle { m: 1.0, x: 1.5, ..reb_particle::default() });
    assert_eq!(r.N, 1, "a particle outside the box must be refused");
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m.contains("outside of box boundaries")),
        "refusing a particle must record an error, got {:?}",
        r.messages
    );

    // Exactly on the face is inside (the C uses strict > / <).
    reb_simulation_add(&mut r, reb_particle { m: 1.0, y: 1.0, ..reb_particle::default() });
    assert_eq!(r.N, 2, "a particle exactly on the box face is accepted");
}

#[test]
fn particle_remove_shifts_the_tail_bit_exactly() {
    let mut r = deterministic_sim();
    for i in 0..5 {
        reb_simulation_add(
            &mut r,
            reb_particle {
                m: (i + 1) as f64,
                x: 0.125 * ((i + 1) as f64),
                vy: -0.25 * ((i + 1) as f64),
                ..reb_particle::default()
            },
        );
    }
    let before = all_particle_bits(&r);

    assert_eq!(reb_simulation_remove_particle(&mut r, 1), 0, "removing index 1 must succeed");
    assert_eq!(r.N, 4, "N must drop by one");
    assert_eq!(r.particles.len(), 4, "the array must be truncated to N");
    let after = all_particle_bits(&r);
    // Particles 0, 2, 3, 4 shift into slots 0, 1, 2, 3 unchanged.
    assert_eq!(after[0], before[0], "particle 0 must be untouched");
    assert_eq!(after[1], before[2], "particle 2 must move into slot 1");
    assert_eq!(after[2], before[3], "particle 3 must move into slot 2");
    assert_eq!(after[3], before[4], "particle 4 must move into slot 3");
    assert_eq!(r.did_modify_particles, 1);
}

#[test]
fn particle_remove_out_of_range_is_refused() {
    let mut r = deterministic_sim();
    for i in 0..3 {
        reb_simulation_add(&mut r, reb_particle { m: (i + 1) as f64, ..reb_particle::default() });
    }
    let before = all_particle_bits(&r);
    assert_eq!(
        reb_simulation_remove_particle(&mut r, 7),
        1,
        "an out-of-range index must return the C failure code 1"
    );
    assert_eq!(r.N, 3, "a refused removal must leave N alone");
    assert_eq!(all_particle_bits(&r), before, "a refused removal must not touch particles");
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m.contains("out of range")),
        "an out-of-range removal must record an error, got {:?}",
        r.messages
    );
}

#[test]
fn particle_remove_updates_N_active() {
    let mut r = deterministic_sim();
    for i in 0..4 {
        reb_simulation_add(&mut r, reb_particle { m: (i + 1) as f64, ..reb_particle::default() });
    }
    r.N_active = 2;

    // Removing a massive (index < N_active) particle shrinks N_active.
    assert_eq!(reb_simulation_remove_particle(&mut r, 0), 0);
    assert_eq!(r.N, 3);
    assert_eq!(r.N_active, 1, "removing an active particle must decrement N_active");

    // Removing a test particle (index >= N_active) leaves it alone.
    assert_eq!(reb_simulation_remove_particle(&mut r, 2), 0);
    assert_eq!(r.N, 2);
    assert_eq!(r.N_active, 1, "removing a test particle must not change N_active");
}

#[test]
fn particle_remove_last_particle_empties_the_simulation() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    // The C checks N == 1 before the range check, so any index works.
    assert_eq!(reb_simulation_remove_particle(&mut r, 0), 0);
    assert_eq!(r.N, 0, "removing the only particle empties the simulation");
    assert_eq!(r.particles.len(), 0);
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::WARNING && m == "Last particle removed."),
        "removing the last particle must warn, got {:?}",
        r.messages
    );
}

#[test]
fn particle_remove_all_resets_the_counters() {
    let mut r = deterministic_sim();
    for i in 0..6 {
        reb_simulation_add(&mut r, reb_particle { m: (i + 1) as f64, ..reb_particle::default() });
    }
    r.N_active = 3;
    reb_simulation_remove_all_particles(&mut r);
    assert_eq!(r.N, 0, "N must be zero");
    assert_eq!(r.particles.len(), 0, "the particle array must be cleared");
    assert_eq!(r.N_var, 0, "N_var must be zero");
    assert_eq!(
        r.N_active,
        usize::MAX,
        "N_active must go back to the SIZE_MAX sentinel"
    );
}

#[test]
fn particle_remove_is_blocked_by_variational_particles() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 1e-3, x: 1.0, vy: 1.0, ..reb_particle::default() });
    let idx = reb_simulation_add_variation_1st_order(&mut r, -1);
    assert_eq!(idx, 0, "the first variational set starts at index 0");
    assert_eq!(r.N_var, 2, "a full-order variation adds one partner per particle");

    let before = all_particle_bits(&r);
    assert_eq!(
        reb_simulation_remove_particle(&mut r, 1),
        1,
        "removal must fail while variational particles exist"
    );
    assert_eq!(r.N, 2, "the refused removal must leave N alone");
    assert_eq!(all_particle_bits(&r), before);
    assert!(
        r.messages.iter().any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR
            && m.contains("Removing particles not supported")),
        "the refusal must be reported, got {:?}",
        r.messages
    );
}

#[test]
fn particle_names_are_interned_and_looked_up() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 2.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 3.0, ..reb_particle::default() });

    let a = reb_simulation_register_name(&mut r, "alpha");
    let b = reb_simulation_register_name(&mut r, "beta");
    let a2 = reb_simulation_register_name(&mut r, "alpha");
    assert_eq!(a, 0, "the first registered name gets index 0");
    assert_eq!(b, 1, "the second registered name gets index 1");
    assert_eq!(a2, a, "registering the same name twice must intern to one index");
    assert_eq!(r.name_list.len(), 2, "the name list must hold two entries");

    reb_particle_set_name(&mut r, 0, Some("alpha"));
    reb_particle_set_name(&mut r, 2, Some("gamma"));
    assert_eq!(r.particles[0].name, Some(0));
    assert_eq!(r.particles[1].name, None);
    assert_eq!(r.particles[2].name, Some(2), "'gamma' is the third interned name");

    assert_eq!(reb_simulation_get_particle_by_name(&r, "alpha"), Some(0));
    assert_eq!(reb_simulation_get_particle_by_name(&r, "gamma"), Some(2));
    assert_eq!(
        reb_simulation_get_particle_by_name(&r, "beta"),
        None,
        "a registered but unassigned name matches no particle"
    );

    reb_particle_set_name(&mut r, 0, None);
    assert_eq!(r.particles[0].name, None, "clearing a name must unset it");
    assert_eq!(reb_simulation_get_particle_by_name(&r, "alpha"), None);
}

#[test]
fn particle_remove_by_name() {
    let mut r = deterministic_sim();
    for i in 0..3 {
        reb_simulation_add(
            &mut r,
            reb_particle { m: (i + 1) as f64, x: i as f64, ..reb_particle::default() },
        );
    }
    reb_particle_set_name(&mut r, 1, Some("middle"));
    let before = all_particle_bits(&r);

    assert_eq!(
        reb_simulation_remove_particle_by_name(&mut r, "nope"),
        1,
        "removing an unknown name must fail"
    );
    assert_eq!(r.N, 3, "a failed name lookup must not remove anything");

    assert_eq!(reb_simulation_remove_particle_by_name(&mut r, "middle"), 0);
    assert_eq!(r.N, 2, "removing by name must drop exactly one particle");
    assert_eq!(all_particle_bits(&r)[0], before[0]);
    assert_eq!(all_particle_bits(&r)[1], before[2], "particle 2 shifts into slot 1");
}

#[test]
fn particle_cmp_detects_every_compared_member() {
    let base = reb_particle {
        x: 1.0, y: 2.0, z: 3.0,
        vx: 4.0, vy: 5.0, vz: 6.0,
        ax: 7.0, ay: 8.0, az: 9.0,
        m: 10.0, r: 11.0,
        name: None,
    };
    assert!(!reb_particle_cmp(base, base), "identical particles must compare equal");

    let mutators: [(&str, fn(&mut reb_particle)); 12] = [
        ("x", |p| p.x += 1.0),
        ("y", |p| p.y += 1.0),
        ("z", |p| p.z += 1.0),
        ("vx", |p| p.vx += 1.0),
        ("vy", |p| p.vy += 1.0),
        ("vz", |p| p.vz += 1.0),
        ("ax", |p| p.ax += 1.0),
        ("ay", |p| p.ay += 1.0),
        ("az", |p| p.az += 1.0),
        ("m", |p| p.m += 1.0),
        ("r", |p| p.r += 1.0),
        ("name", |p| p.name = Some(3)),
    ];
    for (field, f) in mutators {
        let mut other = base;
        f(&mut other);
        assert!(
            reb_particle_cmp(base, other),
            "reb_particle_cmp must notice a change in {}",
            field
        );
    }
}

#[test]
fn particle_arithmetic_helpers_are_exact() {
    // Values chosen so that every operation is exact in binary64.
    let a = reb_particle {
        x: 1.0, y: 2.0, z: 4.0,
        vx: 0.5, vy: 0.25, vz: 0.125,
        ax: 3.0, ay: 5.0, az: 7.0,
        m: 8.0, r: 0.75,
        name: None,
    };
    let b = reb_particle {
        x: 0.5, y: -1.0, z: 0.25,
        vx: 0.25, vy: 0.5, vz: 0.0625,
        ax: 100.0, ay: 200.0, az: 300.0,
        m: 2.0, r: 99.0,
        name: None,
    };

    let mut s = a;
    reb_particle_iadd(&mut s, &b);
    assert_eq!(s.x.to_bits(), (a.x + b.x).to_bits());
    assert_eq!(s.y.to_bits(), (a.y + b.y).to_bits());
    assert_eq!(s.z.to_bits(), (a.z + b.z).to_bits());
    assert_eq!(s.vx.to_bits(), (a.vx + b.vx).to_bits());
    assert_eq!(s.vy.to_bits(), (a.vy + b.vy).to_bits());
    assert_eq!(s.vz.to_bits(), (a.vz + b.vz).to_bits());
    assert_eq!(s.m.to_bits(), (a.m + b.m).to_bits());
    assert_eq!(s.ax.to_bits(), a.ax.to_bits(), "iadd must not touch accelerations");
    assert_eq!(s.r.to_bits(), a.r.to_bits(), "iadd must not touch the radius");

    // iadd then isub of the same particle returns the original bits for
    // these exactly representable values.
    reb_particle_isub(&mut s, &b);
    assert_eq!(particle_bits(&s), particle_bits(&a), "iadd then isub must be exact here");

    let mut m = a;
    reb_particle_imul(&mut m, 2.0);
    reb_particle_imul(&mut m, 0.5);
    assert_eq!(
        particle_bits(&m),
        particle_bits(&a),
        "scaling by 2 then 1/2 is exact for binary64"
    );

    // reb_particle_distance against an independently computed norm.
    let d = reb_particle_distance(&a, &b);
    let (dx, dy, dz) = (a.x - b.x, a.y - b.y, a.z - b.z);
    assert_eq!(d.to_bits(), (dx * dx + dy * dy + dz * dz).sqrt().to_bits());
}

#[test]
fn two_largest_particles_by_radius() {
    let mut r = deterministic_sim();
    for radius in [0.3f64, 1.0, 0.7, 2.5, 0.1] {
        reb_simulation_add(&mut r, reb_particle { m: 1.0, r: radius, ..reb_particle::default() });
    }
    let mut p1 = 0usize;
    let mut p2 = 0usize;
    reb_simulation_two_largest_particles(&r, &mut p1, &mut p2);
    assert_eq!(p1, 3, "the largest radius (2.5) is at index 3");
    assert_eq!(p2, 1, "the second largest radius (1.0) is at index 1");

    // With no particles the C sentinels are returned.
    let empty = deterministic_sim();
    let mut q1 = 0usize;
    let mut q2 = 0usize;
    reb_simulation_two_largest_particles(&empty, &mut q1, &mut q2);
    assert_eq!(q1, usize::MAX, "an empty simulation returns the SIZE_MAX sentinel");
    assert_eq!(q2, usize::MAX);

    // With exactly one particle only the first slot is filled.
    let mut one = deterministic_sim();
    reb_simulation_add(&mut one, reb_particle { m: 1.0, r: 4.0, ..reb_particle::default() });
    let mut s1 = 0usize;
    let mut s2 = 0usize;
    reb_simulation_two_largest_particles(&one, &mut s1, &mut s2);
    assert_eq!(s1, 0);
    assert_eq!(s2, usize::MAX);
}

#[test]
fn rootbox_index_for_particle() {
    let mut r = deterministic_sim();
    // No box configured: the C returns 0 immediately.
    assert_eq!(r.root_size.to_bits(), (-1.0f64).to_bits());
    assert_eq!(
        reb_get_rootbox_for_particle(&r, reb_particle { x: 99.0, ..reb_particle::default() }),
        0,
        "without a root box every particle lives in box 0"
    );

    // 2x2x2 root boxes of size 1 spanning [-1, 1] in each direction.
    r.root_size = 1.0;
    r.N_root_x = 2;
    r.N_root_y = 2;
    r.N_root_z = 2;
    // index = (k * N_root_y + j) * N_root_x + i with
    // i = floor((x + root_size*N/2)/root_size) mod N.
    let cases: [(f64, f64, f64, i32); 4] = [
        (-0.25, -0.25, -0.25, 0), // i=j=k=0
        (0.25, -0.25, -0.25, 1),  // i=1
        (-0.25, 0.25, -0.25, 2),  // j=1
        (0.25, 0.25, 0.25, 7),    // i=j=k=1
    ];
    for (x, y, z, expected) in cases {
        let p = reb_particle { x, y, z, ..reb_particle::default() };
        assert_eq!(
            reb_get_rootbox_for_particle(&r, p),
            expected,
            "particle at ({}, {}, {}) belongs to root box {}",
            x,
            y,
            z,
            expected
        );
    }
}

// =====================================================================
// Step dispatch and exit checks
// =====================================================================

#[test]
fn set_integrator_dispatches_every_name() {
    let names = [
        "none", "sei", "leapfrog", "ias15", "whfast", "saba", "janus", "eos", "mercurius", "bs",
        "trace", "whfast512",
    ];
    for n in names {
        let mut r = deterministic_sim();
        reb_simulation_set_integrator(&mut r, n);
        assert_eq!(r.integrator.name(), n, "set_integrator({}) must select it", n);
        assert!(
            r.messages.is_empty(),
            "selecting {} must not raise a message, got {:?}",
            n,
            r.messages
        );
    }

    // An unknown name is an error and leaves the current integrator.
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "whfast");
    reb_simulation_set_integrator(&mut r, "not_an_integrator");
    assert_eq!(
        r.integrator.name(),
        "whfast",
        "an unknown integrator name must not change the selection"
    );
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m == "Integrator not found."),
        "an unknown integrator name must record 'Integrator not found.', got {:?}",
        r.messages
    );

    // The default integrator of a fresh simulation is IAS15.
    let fresh = reb_simulation_create();
    assert_eq!(fresh.integrator.name(), "ias15", "IAS15 is the C default");
}

#[test]
fn step_with_integrator_none_advances_time_exactly() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.25; // exactly representable
    reb_simulation_add(
        &mut r,
        reb_particle { m: 1.0, x: 1.0, vy: 2.0, ..reb_particle::default() },
    );
    let before = all_particle_bits(&r);

    for k in 1..=8u64 {
        reb_simulation_step(&mut r);
        assert_eq!(
            r.t.to_bits(),
            (0.25 * (k as f64)).to_bits(),
            "after {} steps t must be exactly {}",
            k,
            0.25 * (k as f64)
        );
        assert_eq!(r.steps_done, k, "steps_done must count every step");
        assert_eq!(r.dt_last_done.to_bits(), 0.25f64.to_bits());
    }
    assert_eq!(
        all_particle_bits(&r),
        before,
        "the 'none' integrator must not move any particle"
    );
}

static HEARTBEAT_CALLS: AtomicUsize = AtomicUsize::new(0);
fn count_heartbeat(_r: &mut reb_simulation) {
    HEARTBEAT_CALLS.fetch_add(1, Ordering::SeqCst);
}

#[test]
fn steps_runs_the_heartbeat_once_per_step_plus_once() {
    // simulation.c: run_heartbeat() before the loop and after each step.
    HEARTBEAT_CALLS.store(0, Ordering::SeqCst);
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.5;
    r.heartbeat = Some(count_heartbeat);
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let status = reb_simulation_steps(&mut r, 6);
    assert_eq!(status, REB_STATUS_SUCCESS, "a completed step run reports success");
    assert_eq!(r.steps_done, 6, "exactly six steps must run");
    assert_eq!(
        HEARTBEAT_CALLS.load(Ordering::SeqCst),
        7,
        "the heartbeat runs once before the loop and once after each of the 6 steps"
    );
    assert_eq!(r.t.to_bits(), 3.0f64.to_bits(), "6 steps of dt = 0.5 reach t = 3");
}

static PRE_CALLS: AtomicUsize = AtomicUsize::new(0);
static POST_CALLS: AtomicUsize = AtomicUsize::new(0);
fn count_pre(_r: &mut reb_simulation) {
    PRE_CALLS.fetch_add(1, Ordering::SeqCst);
}
fn count_post(_r: &mut reb_simulation) {
    POST_CALLS.fetch_add(1, Ordering::SeqCst);
}

#[test]
fn pre_and_post_timestep_modifications_run_once_per_step() {
    PRE_CALLS.store(0, Ordering::SeqCst);
    POST_CALLS.store(0, Ordering::SeqCst);
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.5;
    r.pre_timestep_modifications = Some(count_pre);
    r.post_timestep_modifications = Some(count_post);
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    reb_simulation_steps(&mut r, 4);
    assert_eq!(PRE_CALLS.load(Ordering::SeqCst), 4, "pre hook runs once per step");
    assert_eq!(POST_CALLS.load(Ordering::SeqCst), 4, "post hook runs once per step");

    // The presence of these callbacks is recorded in the binary file.
    let buf = reb_binarydata_simulation_to_stream(&mut r);
    let fields = parse_fields(&buf);
    let fp = find_field(&fields, "functionpointers");
    assert_eq!(
        le_i32(&buf, fp.data_pos),
        1,
        "post_timestep_modifications must set the functionpointers flag"
    );
}

#[test]
fn step_clears_the_did_modify_particles_flag() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    assert_eq!(r.did_modify_particles, 1, "adding a particle taints the array");
    reb_simulation_step(&mut r);
    assert_eq!(r.did_modify_particles, 0, "a timestep clears the taint flag");
}

#[test]
fn integrate_with_exact_finish_time_lands_on_tmax_and_restores_dt() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.1;
    r.exact_finish_time = 1;
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let tmax = 1.0f64;
    let status = reb_simulation_integrate(&mut r, tmax);
    assert_eq!(status, REB_STATUS_SUCCESS, "the integration must succeed");
    // reb_check_exit accepts |t - tmax| < 1e-12 * |tmax| as "arrived".
    let tscale = 1e-12 * tmax.abs();
    assert!(
        (r.t - tmax).abs() < tscale,
        "exact_finish_time must land within {} of tmax, got t = {}",
        tscale,
        r.t
    );
    assert_eq!(
        r.dt.to_bits(),
        0.1f64.to_bits(),
        "the shrunk last timestep must be replaced by the last full dt"
    );
}

#[test]
fn integrate_without_exact_finish_time_overshoots_by_less_than_one_step() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.1;
    r.exact_finish_time = 0;
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let tmax = 1.0f64;
    let status = reb_simulation_integrate(&mut r, tmax);
    assert_eq!(status, REB_STATUS_SUCCESS);
    assert!(
        r.t >= tmax,
        "without exact_finish_time the loop stops at or past tmax, got {}",
        r.t
    );
    assert!(
        r.t - r.dt_last_done < tmax,
        "the step before the last one must still have been below tmax (t = {}, dt = {})",
        r.t,
        r.dt_last_done
    );
    assert_eq!(
        r.dt.to_bits(),
        0.1f64.to_bits(),
        "dt is never shrunk when exact_finish_time is 0"
    );
}

#[test]
fn integrate_backwards_flips_the_timestep_sign() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.25; // positive; integrating to a negative tmax must flip it
    r.exact_finish_time = 1;
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    let tmax = -2.0f64;
    let status = reb_simulation_integrate(&mut r, tmax);
    assert_eq!(status, REB_STATUS_SUCCESS);
    assert!(r.dt < 0.0, "dt must be negative after a backwards integration, got {}", r.dt);
    assert_eq!(r.dt.to_bits(), (-0.25f64).to_bits(), "|dt| must be unchanged");
    assert_eq!(
        r.t.to_bits(),
        tmax.to_bits(),
        "with dt = -0.25 the time steps land exactly on -2.0"
    );
    assert_eq!(r.steps_done, 8, "8 steps of 0.25 cover 2.0");
}

#[test]
fn exit_max_distance_reports_escape() {
    // Massless particles under G = 1 feel no force, so the trajectory is
    // a straight line and the escape time is known analytically.
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "leapfrog");
    r.dt = 0.1;
    r.exit_max_distance = 3.0;
    reb_simulation_add(&mut r, reb_particle::default());
    reb_simulation_add(
        &mut r,
        reb_particle { x: 1.0, vx: 1.0, ..reb_particle::default() },
    );

    let status = reb_simulation_integrate(&mut r, 100.0);
    assert_eq!(status, REB_STATUS_ESCAPE, "a particle beyond exit_max_distance escapes");
    let p = r.particles[1];
    let d = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    assert!(d > 3.0, "the escaping particle must be past exit_max_distance, d = {}", d);
    assert!(
        d < 3.0 + 0.1 + 1e-9,
        "the run must stop at the first step past the threshold, d = {}",
        d
    );
    assert!(r.t < 100.0, "the integration must stop early, t = {}", r.t);
}

#[test]
fn exit_min_distance_reports_encounter() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "leapfrog");
    r.dt = 0.1;
    r.exit_min_distance = 0.5;
    reb_simulation_add(&mut r, reb_particle::default());
    reb_simulation_add(
        &mut r,
        reb_particle { x: 2.0, vx: -1.0, ..reb_particle::default() },
    );

    let status = reb_simulation_integrate(&mut r, 100.0);
    assert_eq!(status, REB_STATUS_ENCOUNTER, "a close pair reports an encounter");
    let d = reb_particle_distance(&r.particles[0], &r.particles[1]);
    assert!(d < 0.5, "the pair must be closer than exit_min_distance, d = {}", d);
    assert!(
        d > 0.5 - 0.1 - 1e-9,
        "the run must stop at the first step inside the threshold, d = {}",
        d
    );
}

#[test]
fn integrate_without_particles_exits_immediately() {
    let mut r = deterministic_sim();
    assert_eq!(r.N, 0);
    let status = reb_simulation_integrate(&mut r, 1.0);
    assert_eq!(
        status,
        REB_STATUS_NO_PARTICLES,
        "an empty simulation must report REB_STATUS_NO_PARTICLES"
    );
    assert_eq!(r.steps_done, 0, "no step may run without particles");
    assert_eq!(r.t.to_bits(), 0.0f64.to_bits(), "time must not advance");
    assert!(
        r.messages.iter().any(|(t, m)| *t == REB_MESSAGE_TYPE::WARNING
            && m == "No particles in simulation. Will exit."),
        "the C warns before exiting, got {:?}",
        r.messages
    );
}

fn stop_at_half(r: &mut reb_simulation) {
    if r.t >= 0.5 {
        reb_simulation_stop(r);
    }
}

#[test]
fn heartbeat_stop_ends_the_integration_with_user_status() {
    let mut r = deterministic_sim();
    reb_simulation_set_integrator(&mut r, "none");
    r.dt = 0.25; // exact, so t = 0, 0.25, 0.5, ... exactly
    r.heartbeat = Some(stop_at_half);
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });

    // reb_simulation_stop only sets status; reb_simulation_integrate
    // resets status to RUNNING on entry, so the stop has to come from
    // the heartbeat, which runs after every step.
    let status = reb_simulation_integrate(&mut r, 10.0);
    assert_eq!(status, REB_STATUS_USER, "a heartbeat stop reports REB_STATUS_USER");
    assert_eq!(r.t.to_bits(), 0.5f64.to_bits(), "the run stops at t = 0.5 exactly");
    assert_eq!(r.steps_done, 2, "two steps of dt = 0.25 reach t = 0.5");
}

#[test]
fn update_acceleration_honours_the_gravity_selection() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    reb_simulation_add(&mut r, reb_particle { m: 1.0, x: 1.0, ..reb_particle::default() });
    r.particles[0].ax = 123.0;
    r.particles[1].ay = -456.0;

    r.gravity = REB_GRAVITY::NONE;
    reb_simulation_update_acceleration(&mut r);
    for i in 0..r.N {
        assert_eq!(r.particles[i].ax.to_bits(), 0.0f64.to_bits());
        assert_eq!(r.particles[i].ay.to_bits(), 0.0f64.to_bits());
        assert_eq!(r.particles[i].az.to_bits(), 0.0f64.to_bits());
    }

    // Two unit masses one unit apart, G = 1: |a| = 1 towards each other.
    r.gravity = REB_GRAVITY::BASIC;
    reb_simulation_update_acceleration(&mut r);
    assert_eq!(
        r.particles[0].ax.to_bits(),
        1.0f64.to_bits(),
        "particle 0 is pulled towards +x with |a| = G m / d^2 = 1"
    );
    assert_eq!(r.particles[1].ax.to_bits(), (-1.0f64).to_bits());

    // CUSTOM without a callback is an error and leaves accelerations alone.
    r.gravity = REB_GRAVITY::CUSTOM;
    r.messages.clear();
    reb_simulation_update_acceleration(&mut r);
    assert!(
        r.messages
            .iter()
            .any(|(t, m)| *t == REB_MESSAGE_TYPE::ERROR && m.contains("gravity_custom")),
        "CUSTOM gravity without a callback must be an error, got {:?}",
        r.messages
    );
}

#[test]
fn repeated_runs_are_bit_identical() {
    let run = || {
        let mut r = two_body("whfast", 0.0123, 0.3);
        reb_simulation_steps(&mut r, 500);
        (r.t.to_bits(), all_particle_bits(&r))
    };
    let (t1, b1) = run();
    let (t2, b2) = run();
    assert_eq!(t1, t2, "the simulation time must be bit-identical between runs");
    assert_eq!(b1, b2, "particle state must be bit-identical between runs");
}

#[test]
fn whfast_conserves_energy_and_angular_momentum() {
    // Symplectic integrator on a two-body problem: the energy error must
    // stay bounded and the z angular momentum is conserved to roundoff.
    let period = 2.0 * std::f64::consts::PI / (1.001f64).sqrt();
    let mut r = two_body("whfast", period / 200.0, 0.1);
    let e0 = reb_simulation_energy(&r);
    let l0 = reb_simulation_angular_momentum(&r);
    assert!(e0 < 0.0, "a bound two-body system has negative energy, got {}", e0);

    reb_simulation_steps(&mut r, 20_000); // 100 orbits
    let e1 = reb_simulation_energy(&r);
    let l1 = reb_simulation_angular_momentum(&r);

    let de = ((e1 - e0) / e0).abs();
    assert!(
        de < 1e-9,
        "relative energy drift after 100 orbits is {} (bound 1e-9)",
        de
    );
    let dl = ((l1.z - l0.z) / l0.z).abs();
    assert!(
        dl < 1e-12,
        "relative Lz drift after 100 orbits is {} (bound 1e-12)",
        dl
    );
    assert!(
        l1.x.abs() < 1e-14 && l1.y.abs() < 1e-14,
        "a planar orbit must keep Lx = Ly = 0, got ({}, {})",
        l1.x,
        l1.y
    );
}

// =====================================================================
// ASCII output helpers
// =====================================================================

#[test]
fn output_check_phase_known_values() {
    let mut r = deterministic_sim();
    reb_simulation_add(&mut r, reb_particle { m: 1.0, ..reb_particle::default() });
    r.dt = 0.1;

    r.t = 0.0;
    assert!(
        reb_simulation_output_check(&r, 1.0),
        "the C always outputs at t == 0"
    );

    // floor(t/interval) == floor((t-dt)/interval) -> no output.
    r.t = 0.5;
    assert!(
        !reb_simulation_output_check(&r, 1.0),
        "t = 0.5, dt = 0.1 stays inside the same interval bin"
    );

    // Crossing an interval boundary -> output.
    r.t = 1.05;
    assert!(
        reb_simulation_output_check(&r, 1.0),
        "t = 1.05, dt = 0.1 crosses the boundary at 1.0"
    );

    // A phase of 0.5 shifts the boundary by half an interval: shift =
    // t + interval*phase, so the check fires when t crosses 0.5.
    r.t = 0.55;
    assert!(
        reb_simulation_output_check_phase(&r, 1.0, 0.5),
        "with phase 0.5 the boundary sits at t = 0.5 and t = 0.55 just crossed it"
    );
    r.t = 0.3;
    assert!(
        !reb_simulation_output_check_phase(&r, 1.0, 0.5),
        "t = 0.3 with phase 0.5 stays inside a bin"
    );
    assert!(
        !reb_simulation_output_check(&r, 1.0),
        "t = 0.3 with phase 0 also stays inside a bin"
    );
}

#[test]
fn fmt_e_matches_c_printf() {
    // C's "%e": six digits of precision and at least two exponent digits.
    assert_eq!(fmt_e(1.0), "1.000000e+00");
    assert_eq!(fmt_e(0.0), "0.000000e+00");
    assert_eq!(fmt_e(-1.5e-5), "-1.500000e-05");
    assert_eq!(fmt_e(1234567.0), "1.234567e+06");
    assert_eq!(fmt_e(1e100), "1.000000e+100");
    assert_eq!(fmt_e(-2.5), "-2.500000e+00");
}

#[test]
fn fmt_fixed_width_matches_c_printf() {
    // C's "%- 9f": width 9, left justified, space for a positive sign.
    assert_eq!(fmt_f_space_left9(1.5), " 1.500000");
    assert_eq!(fmt_f_space_left9(-1.5), "-1.500000");
    assert_eq!(fmt_f_space_left9(0.0), " 0.000000");
    assert_eq!(fmt_f_space_left9(-0.0), "-0.000000");
    for v in [1.5f64, -1.5, 0.0, -0.0, 12345.6789] {
        assert!(
            fmt_f_space_left9(v).len() >= 9,
            "the field is at least 9 characters wide for {}",
            v
        );
    }
    // C's "%- 9d".
    assert_eq!(fmt_d_space_left9(42), " 42      ");
    assert_eq!(fmt_d_space_left9(-7), "-7       ");
    assert_eq!(fmt_d_space_left9(0), " 0       ");
    assert_eq!(fmt_d_space_left9(123456789), " 123456789");
}
