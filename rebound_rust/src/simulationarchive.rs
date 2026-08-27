//! simulationarchive.rs — tools for creating and reading
//! Simulationarchive binary files (from simulationarchive.c/h;
//! (c) 2016 Hanno Rein).
//!
//! An archive is the initial full binary (binarydata.rs format)
//! followed by incremental diff blobs, separated by 12-byte
//! `reb_simulationarchive_blob` markers. The formats are byte
//! compatible with the C build (up to serialized pointer values, see
//! binarydata.rs).
//!
//! The C keeps a `FILE*` open in `struct reb_simulationarchive`; here
//! the stream is a boxed `Read + Seek` (a `File`, or a `Cursor` over a
//! byte buffer replacing the C's non-portable fmemopen).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use crate::binarydata::*;
use crate::simulation::reb_simulation_create;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;
use crate::reb_version_str;
use std::io::{Cursor, Read, Seek, SeekFrom, Write};

/// Stream abstraction (C: `FILE*` from fopen or fmemopen).
pub trait ReadSeek: Read + Seek {}
impl<T: Read + Seek> ReadSeek for T {}

/// simulationarchive.h `struct reb_simulationarchive_blob` — used in
/// the binary file to identify data blobs (snapshots).
#[derive(Clone, Copy, Debug, Default)]
pub struct reb_simulationarchive_blob {
    /// Index of previous blob (binary file is 0, first blob is 1).
    pub index: i32,
    /// Offset to beginning of previous blob (size of previous blob).
    pub offset_prev: i32,
    /// Offset to end of following blob (size of following blob).
    pub offset_next: i32,
}

fn blob_to_bytes(b: &reb_simulationarchive_blob) -> [u8; REB_SA_BLOB_SIZE] {
    let mut out = [0u8; REB_SA_BLOB_SIZE];
    out[0..4].copy_from_slice(&b.index.to_le_bytes());
    out[4..8].copy_from_slice(&b.offset_prev.to_le_bytes());
    out[8..12].copy_from_slice(&b.offset_next.to_le_bytes());
    out
}

fn blob_from_bytes(d: &[u8]) -> reb_simulationarchive_blob {
    reb_simulationarchive_blob {
        index: i32::from_le_bytes(d[0..4].try_into().unwrap_or([0; 4])),
        offset_prev: i32::from_le_bytes(d[4..8].try_into().unwrap_or([0; 4])),
        offset_next: i32::from_le_bytes(d[8..12].try_into().unwrap_or([0; 4])),
    }
}

/// simulationarchive.h `struct reb_simulationarchive`.
pub struct reb_simulationarchive {
    /// Open stream (C: `FILE* inf`).
    pub inf: Option<Box<dyn ReadSeek>>,
    /// Filename of open file (None for memory-backed archives).
    pub filename: Option<String>,
    /// Simulationarchive version.
    pub version: i32,
    /// REBOUND version used to save the SA.
    pub reb_version_major: i32,
    pub reb_version_minor: i32,
    pub reb_version_patch: i32,
    /// Interval setting used to create the SA (if used).
    pub auto_interval: f64,
    /// Walltime setting used to create the SA (if used).
    pub auto_walltime: f64,
    /// Steps in-between SA snapshots (if used).
    pub auto_step: u64,
    /// Total number of snapshots (including initial binary).
    pub nblobs: i64,
    /// Index of offsets in file (length nblobs).
    pub offset: Vec<u64>,
    /// Index of simulation times in file (length nblobs).
    pub t: Vec<f64>,
}

impl Default for reb_simulationarchive {
    fn default() -> Self {
        reb_simulationarchive {
            inf: None,
            filename: None,
            version: 0,
            reb_version_major: 0,
            reb_version_minor: 0,
            reb_version_patch: 0,
            auto_interval: 0.,
            auto_walltime: 0.,
            auto_step: 0,
            nblobs: 0,
            offset: Vec::new(),
            t: Vec::new(),
        }
    }
}

fn read_exact_or_eof<R: Read + ?Sized>(inf: &mut R, out: &mut [u8]) -> bool {
    let mut got = 0;
    while got < out.len() {
        match inf.read(&mut out[got..]) {
            Ok(0) => return false,
            Ok(n) => got += n,
            Err(_) => return false,
        }
    }
    true
}

/// simulationarchive.c
/// `reb_simulationarchive_read_from_stream_with_messages` — builds the
/// snapshot index. Assumes `sa.inf` is set to an open stream.
pub fn reb_simulationarchive_read_from_stream_with_messages(
    sa: &mut reb_simulationarchive,
    warnings: &mut REB_BINARYDATA_ERROR_CODE,
) {
    let inf = match sa.inf.as_mut() {
        Some(f) => f,
        None => {
            *warnings |= REB_BINARYDATA_ERROR_NOFILE;
            return;
        }
    };

    // Get version
    let _ = inf.seek(SeekFrom::Start(0));
    sa.version = 0;
    let mut t0 = 0.0_f64;
    let _ = t0;
    sa.reb_version_major = 0;
    sa.reb_version_minor = 0;
    sa.reb_version_patch = 0;

    loop {
        let mut fh = [0u8; REB_BINARYDATA_FIELD_SIZE];
        if !read_exact_or_eof(inf.as_mut(), &mut fh) {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break;
        }
        let size_name = u64::from_le_bytes(fh[0..8].try_into().unwrap_or([0; 8]));
        let size_data = u64::from_le_bytes(fh[8..16].try_into().unwrap_or([0; 8]));

        if size_name == reb_binarydata_header {
            // Input header.
            let bufsize = 64 - REB_BINARYDATA_FIELD_SIZE;
            let mut readbuf = vec![0u8; bufsize];
            let header = "REBOUND Binary File. Version: ";
            let curv = format!("{}{}", &header[REB_BINARYDATA_FIELD_SIZE..], reb_version_str);
            let objects = read_exact_or_eof(inf.as_mut(), &mut readbuf);
            // Finding version_major/minor/patch
            let (mut c1, mut c2, mut c3) = (0usize, 0usize, 0usize);
            for c in 0..bufsize {
                if c2 != 0 && c3 == 0 && readbuf[c] == b'.' {
                    c3 = c;
                }
                if c1 != 0 && c2 == 0 && readbuf[c] == b'.' {
                    c2 = c;
                }
                if c1 == 0 && readbuf[c] == b':' {
                    c1 = c;
                }
            }
            if c1 == 0 || c2 == 0 || c3 == 0 {
                *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            } else {
                // C's `atoi`, faithfully: skip leading whitespace, take an
                // optional sign, then the run of digits, and yield 0 if
                // there are none.
                //
                // The whitespace skip is load-bearing. The header reads
                // "REBOUND Binary File. Version: 5.1.1" and the C splits it
                // at the ':' and the two '.'s, so the major-version substring
                // is " 5" WITH the leading space. C's atoi(" 5") is 5; a
                // digits-only scan stops at the space and yields 0, which
                // silently reported major version 0 for every archive.
                let atoi = |bytes: &[u8]| -> i32 {
                    let mut it = bytes.iter().peekable();
                    while let Some(b) = it.peek() {
                        // C: isspace()
                        if matches!(**b, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r') {
                            it.next();
                        } else {
                            break;
                        }
                    }
                    let mut sign = 1i32;
                    if let Some(b) = it.peek() {
                        if **b == b'-' {
                            sign = -1;
                            it.next();
                        } else if **b == b'+' {
                            it.next();
                        }
                    }
                    let s: String = it
                        .take_while(|b| b.is_ascii_digit())
                        .map(|&b| b as char)
                        .collect();
                    sign * s.parse::<i32>().unwrap_or(0)
                };
                sa.reb_version_patch = atoi(&readbuf[c3 + 1..std::cmp::min(c3 + 4, bufsize)]);
                sa.reb_version_minor = atoi(&readbuf[c2 + 1..c3]);
                sa.reb_version_major = atoi(&readbuf[c1 + 1..c2]);
            }
            if !objects {
                *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            } else {
                let mut matches = true;
                for (i, b) in curv.as_bytes().iter().enumerate() {
                    if i >= bufsize || readbuf[i] != *b {
                        matches = false;
                        break;
                    }
                }
                if matches && curv.len() < bufsize && readbuf[curv.len()] != 0 {
                    matches = false;
                }
                if !matches {
                    *warnings |= REB_BINARYDATA_WARNING_VERSION;
                }
            }
            continue;
        }

        if size_name == 0 || size_name as usize > REB_STRING_SIZE_MAX {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break;
        }

        let mut name_bytes = vec![0u8; size_name as usize];
        if !read_exact_or_eof(inf.as_mut(), &mut name_bytes) {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break;
        }
        let nul = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
        let name = String::from_utf8_lossy(&name_bytes[..nul]).to_string();

        let mut data = vec![0u8; size_data as usize];
        match name.as_str() {
            "t" => {
                if read_exact_or_eof(inf.as_mut(), &mut data) && data.len() >= 8 {
                    t0 = f64::from_le_bytes(data[0..8].try_into().unwrap_or([0; 8]));
                    let _ = t0;
                }
            }
            "simulationarchive_version" => {
                if read_exact_or_eof(inf.as_mut(), &mut data) && data.len() >= 4 {
                    sa.version = i32::from_le_bytes(data[0..4].try_into().unwrap_or([0; 4]));
                }
            }
            "simulationarchive_auto_walltime" => {
                if read_exact_or_eof(inf.as_mut(), &mut data) && data.len() >= 8 {
                    sa.auto_walltime = f64::from_le_bytes(data[0..8].try_into().unwrap_or([0; 8]));
                }
            }
            "simulationarchive_auto_interval" => {
                if read_exact_or_eof(inf.as_mut(), &mut data) && data.len() >= 8 {
                    sa.auto_interval = f64::from_le_bytes(data[0..8].try_into().unwrap_or([0; 8]));
                }
            }
            "simulationarchive_auto_step" => {
                if read_exact_or_eof(inf.as_mut(), &mut data) && data.len() >= 8 {
                    sa.auto_step = u64::from_le_bytes(data[0..8].try_into().unwrap_or([0; 8]));
                }
            }
            _ => {
                let _ = inf.seek(SeekFrom::Current(size_data as i64));
            }
        }
        if name == "end" {
            break;
        }
    }

    // Make index
    if sa.version < 5 {
        // No longer supported
        sa.filename = None;
        sa.inf = None;
        *warnings |= REB_BINARYDATA_ERROR_OLD;
        return;
    }
    // New version
    let inf = sa.inf.as_mut().unwrap();
    let _ = inf.seek(SeekFrom::Start(0));
    sa.nblobs = 0;
    sa.t.clear();
    sa.offset.clear();
    let mut read_error = false;
    loop {
        let i = sa.t.len();
        let here = inf.seek(SeekFrom::Current(0)).unwrap_or(0);
        sa.offset.push(here);
        sa.t.push(0.);
        let mut blob_finished = false;
        while !blob_finished && !read_error {
            let mut fh = [0u8; REB_BINARYDATA_FIELD_SIZE];
            let r1 = read_exact_or_eof(inf.as_mut(), &mut fh);
            let size_name = u64::from_le_bytes(fh[0..8].try_into().unwrap_or([0; 8]));
            let size_data = u64::from_le_bytes(fh[8..16].try_into().unwrap_or([0; 8]));
            if r1 && size_name == reb_binarydata_header {
                if inf
                    .seek(SeekFrom::Current((64 - REB_BINARYDATA_FIELD_SIZE) as i64))
                    .is_err()
                {
                    read_error = true;
                }
                continue;
            }
            if !r1 || size_name == 0 || size_name as usize > REB_STRING_SIZE_MAX {
                read_error = true;
                break;
            }
            let mut name_bytes = vec![0u8; size_name as usize];
            if !read_exact_or_eof(inf.as_mut(), &mut name_bytes) {
                read_error = true;
                break;
            }
            let nul = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
            let fname = String::from_utf8_lossy(&name_bytes[..nul]).to_string();
            if fname == "t" {
                let mut d = vec![0u8; size_data as usize];
                if !read_exact_or_eof(inf.as_mut(), &mut d) || d.len() < 8 {
                    read_error = true;
                } else {
                    sa.t[i] = f64::from_le_bytes(d[0..8].try_into().unwrap_or([0; 8]));
                }
            } else if fname == "end" {
                blob_finished = true;
            } else if inf.seek(SeekFrom::Current(size_data as i64)).is_err() {
                read_error = true;
            }
        }
        if read_error {
            // Error during reading. Current snapshot is corrupt.
            sa.offset.pop();
            sa.t.pop();
            break;
        }
        // Attempt to read next blob marker.
        let mut bb = [0u8; REB_SA_BLOB_SIZE];
        let r3 = read_exact_or_eof(inf.as_mut(), &mut bb);
        let blob = blob_from_bytes(&bb);
        let mut next_blob_is_corrupted = false;
        if !r3 {
            next_blob_is_corrupted = true;
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
        }
        if i > 0 {
            // Checking the offsets. Acts like a checksum.
            let here = inf.seek(SeekFrom::Current(0)).unwrap_or(0);
            if (blob.offset_prev as i64) + (REB_SA_BLOB_SIZE as i64)
                != (here as i64) - (sa.offset[i] as i64)
            {
                // Offsets don't work. Next snapshot is corrupted; assume current one as well.
                sa.offset.pop();
                sa.t.pop();
                read_error = true;
                break;
            }
        }
        // All tests passed. Accept current snapshot.
        sa.nblobs = (i + 1) as i64;
        if blob.offset_next == 0 || next_blob_is_corrupted {
            // Last blob.
            break;
        }
    }
    if read_error {
        if sa.nblobs > 0 {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
        } else {
            sa.inf = None;
            sa.filename = None;
            sa.t.clear();
            sa.offset.clear();
            *warnings |= REB_BINARYDATA_ERROR_SEEK;
        }
    }
}

/// simulationarchive.c `reb_simulationarchive_create_from_file_with_messages`.
pub fn reb_simulationarchive_create_from_file_with_messages(
    filename: &str,
    warnings: &mut REB_BINARYDATA_ERROR_CODE,
) -> reb_simulationarchive {
    let mut sa = reb_simulationarchive::default();
    match std::fs::File::open(filename) {
        Ok(f) => sa.inf = Some(Box::new(f)),
        Err(_) => sa.inf = None,
    }
    sa.filename = Some(filename.to_string());
    if sa.inf.is_none() {
        *warnings |= REB_BINARYDATA_ERROR_NOFILE;
        return sa;
    }
    reb_simulationarchive_read_from_stream_with_messages(&mut sa, warnings);
    sa
}

/// simulationarchive.c `reb_simulationarchive_create_from_file` —
/// returns None if the file does not exist (matching the C's NULL).
pub fn reb_simulationarchive_create_from_file(filename: &str) -> Option<reb_simulationarchive> {
    let mut warnings = REB_BINARYDATA_WARNING_NONE;
    let sa = reb_simulationarchive_create_from_file_with_messages(filename, &mut warnings);
    if warnings & REB_BINARYDATA_ERROR_NOFILE != 0 {
        return None;
    }
    let mut tmp = crate::simulation::reb_simulation_create(); // C passes NULL; messages go nowhere
    let _ = reb_binarydata_process_warnings(&mut tmp, warnings);
    Some(sa)
}

/// simulationarchive.c `reb_simulationarchive_create_from_buffer_with_messages`
/// (the C uses fmemopen; Rust uses a Cursor over a copy of the buffer).
pub fn reb_simulationarchive_create_from_buffer_with_messages(
    buffer: &[u8],
    warnings: &mut REB_BINARYDATA_ERROR_CODE,
) -> reb_simulationarchive {
    let mut sa = reb_simulationarchive::default();
    sa.inf = Some(Box::new(Cursor::new(buffer.to_vec())));
    sa.filename = None;
    reb_simulationarchive_read_from_stream_with_messages(&mut sa, warnings);
    sa
}

/// simulationarchive.c `reb_simulationarchive_free` (Drop handles the
/// resources; provided for name parity with the C API).
pub fn reb_simulationarchive_free(sa: reb_simulationarchive) {
    drop(sa);
}

/// simulationarchive.c
/// `reb_simulation_init_from_simulationarchive_with_messages`.
pub fn reb_simulation_init_from_simulationarchive_with_messages(
    r: &mut reb_simulation,
    sa: &mut reb_simulationarchive,
    snapshot: i64,
    warnings: &mut REB_BINARYDATA_ERROR_CODE,
) {
    if sa.inf.is_none() {
        *warnings |= REB_BINARYDATA_ERROR_FILENOTOPEN;
        return;
    }
    let mut snapshot = snapshot;
    if snapshot < 0 {
        snapshot += sa.nblobs;
    }
    if snapshot >= sa.nblobs || snapshot < 0 {
        *warnings |= REB_BINARYDATA_ERROR_OUTOFRANGE;
        return;
    }

    // load original binary file
    r.simulationarchive_filename = None;
    {
        let inf = sa.inf.as_mut().unwrap();
        let _ = inf.seek(SeekFrom::Start(0));
        reb_binarydata_input_fields(r, inf, warnings);
    }

    // Done?
    if snapshot == 0 {
        return;
    }

    // Read SA snapshot
    let inf = sa.inf.as_mut().unwrap();
    if inf.seek(SeekFrom::Start(sa.offset[snapshot as usize])).is_err() {
        *warnings |= REB_BINARYDATA_ERROR_SEEK;
        return;
    }
    if r.simulationarchive_version < 5 {
        *warnings |= REB_BINARYDATA_ERROR_OLD;
        return;
    }
    // Version 5 or higher
    reb_binarydata_input_fields(r, inf, warnings);
}

/// simulationarchive.c `reb_simulation_create_from_simulationarchive`.
/// Returns None if an error occurred (matching the C's NULL).
pub fn reb_simulation_create_from_simulationarchive(
    sa: &mut reb_simulationarchive,
    snapshot: i64,
) -> Option<reb_simulation> {
    let mut warnings = REB_BINARYDATA_WARNING_NONE;
    let mut r = reb_simulation_create();
    reb_simulation_init_from_simulationarchive_with_messages(&mut r, sa, snapshot, &mut warnings);
    if reb_binarydata_process_warnings(&mut r, warnings) < 0 {
        return None;
    }
    Some(r)
}

/// simulationarchive.c `reb_simulation_create_from_file`.
pub fn reb_simulation_create_from_file(filename: &str, snapshot: i64) -> Option<reb_simulation> {
    let mut warnings = REB_BINARYDATA_WARNING_NONE;
    let mut r = reb_simulation_create();

    let mut sa = reb_simulationarchive_create_from_file_with_messages(filename, &mut warnings);
    if warnings & REB_BINARYDATA_ERROR_NOFILE != 0 {
        // Don't output an error if file does not exist, just return None.
        return None;
    }
    let _ = reb_binarydata_process_warnings(&mut r, warnings);
    reb_simulation_init_from_simulationarchive_with_messages(&mut r, &mut sa, snapshot, &mut warnings);
    reb_simulationarchive_free(sa);
    if reb_binarydata_process_warnings(&mut r, warnings) < 0 {
        return None;
    }
    Some(r)
}

/// simulationarchive.c `reb_simulationarchive_heartbeat` — internal
/// function to handle outputs for the Simulationarchive.
pub fn reb_simulationarchive_heartbeat(r: &mut reb_simulation) {
    if r.simulationarchive_filename.is_some() {
        let mut modes = 0;
        if r.simulationarchive_auto_interval != 0. {
            modes += 1;
        }
        if r.simulationarchive_auto_walltime != 0. {
            modes += 1;
        }
        if r.simulationarchive_auto_step != 0 {
            modes += 1;
        }
        if modes > 1 {
            reb_simulation_error(r, "Only use one of simulationarchive_auto_interval, simulationarchive_auto_walltime, or simulationarchive_auto_step");
        }
        if r.simulationarchive_auto_interval != 0. {
            let sign = if r.dt > 0. { 1. } else { -1. };
            if sign * r.simulationarchive_next <= sign * r.t {
                r.simulationarchive_next += sign * r.simulationarchive_auto_interval;
                // Snap
                reb_simulation_save_to_file(r, None);
            }
        }
        if r.simulationarchive_auto_step != 0 {
            if r.simulationarchive_next_step <= r.steps_done {
                r.simulationarchive_next_step += r.simulationarchive_auto_step;
                // Snap
                reb_simulation_save_to_file(r, None);
            }
        }
        if r.simulationarchive_auto_walltime != 0. {
            if r.simulationarchive_next <= r.walltime {
                r.simulationarchive_next += r.simulationarchive_auto_walltime;
                // Snap
                reb_simulation_save_to_file(r, None);
            }
        }
    }
}

/// simulationarchive.c `reb_simulation_save_to_file` — writes a full
/// binary if the file does not exist, otherwise appends a diff
/// snapshot. `filename: None` uses `r.simulationarchive_filename`.
pub fn reb_simulation_save_to_file(r: &mut reb_simulation, filename: Option<&str>) {
    if r.simulationarchive_version < 5 {
        reb_simulation_error(
            r,
            "Writing Simulationarchives with a version < 5 is no longer supported.\n",
        );
        return;
    }
    let filename: String = match filename {
        Some(f) => f.to_string(),
        None => match &r.simulationarchive_filename {
            Some(f) => f.clone(),
            None => {
                reb_simulation_error(r, "Can not open file.");
                return;
            }
        },
    };
    if !std::path::Path::new(&filename).exists() {
        // File does not exist. Output binary.
        let buf = reb_binarydata_simulation_to_stream(r);
        match std::fs::File::create(&filename) {
            Ok(mut of) => {
                let _ = of.write_all(&buf);
            }
            Err(_) => {
                reb_simulation_error(r, "Can not open file.");
            }
        }
        return;
    }
    // File exists, append snapshot.
    let mut existing = match std::fs::read(&filename) {
        Ok(b) => b,
        Err(_) => {
            reb_simulation_error(r, "Can not open file.");
            return;
        }
    };
    // Find the end of the initial binary (scan fields from offset 64).
    let mut pos = 64usize;
    let size_old;
    loop {
        if pos + REB_BINARYDATA_FIELD_SIZE > existing.len() {
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }
        let size_name = u64::from_le_bytes(existing[pos..pos + 8].try_into().unwrap_or([0; 8]));
        let size_data =
            u64::from_le_bytes(existing[pos + 8..pos + 16].try_into().unwrap_or([0; 8]));
        if size_name != reb_binarydata_header
            && (size_name == 0 || size_name as usize > REB_STRING_SIZE_MAX)
        {
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }
        pos += REB_BINARYDATA_FIELD_SIZE;
        if size_name == reb_binarydata_header {
            // (Cannot occur mid-binary in well-formed files; the C
            // would misparse here as well. Continue past.)
            pos += 64 - REB_BINARYDATA_FIELD_SIZE;
            continue;
        }
        let name_end = pos + size_name as usize;
        if name_end > existing.len() {
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }
        let name_bytes = &existing[pos..name_end];
        let nul = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
        let is_end = &name_bytes[..nul] == b"end";
        pos = name_end + size_data as usize;
        if pos > existing.len() {
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }
        if is_end {
            size_old = pos;
            break;
        }
    }

    if size_old + REB_SA_BLOB_SIZE > existing.len() {
        reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
        return;
    }

    let buf_old = &existing[..size_old];

    // Create buffer containing current binary file.
    let buf_new = reb_binarydata_simulation_to_stream(r);

    // Create buffer containing diff.
    let (_differ, buf_diff) =
        reb_binarydata_diff(buf_old, &buf_new, REB_BINARYDATA_OUTPUT_STREAM);

    // Validity check of the tail blob (simplified recovery: the C walks
    // the blob chain to find the last valid snapshot; here a corrupt
    // tail aborts the append with the same warning).
    if existing.len() < REB_SA_BLOB_SIZE {
        reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
        return;
    }
    let tail_start = existing.len() - REB_SA_BLOB_SIZE;
    let mut tail_blob = blob_from_bytes(&existing[tail_start..]);
    if tail_blob.offset_next != 0 {
        reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
        return;
    }

    // Update blob info and write diff to binary file.
    let end_len = "end".len() + 1;
    tail_blob.offset_next =
        (buf_diff.len() + REB_BINARYDATA_FIELD_SIZE + end_len) as i32;
    existing[tail_start..].copy_from_slice(&blob_to_bytes(&tail_blob));
    existing.extend_from_slice(&buf_diff);
    // end field
    existing.extend_from_slice(&(end_len as u64).to_le_bytes());
    existing.extend_from_slice(&0u64.to_le_bytes());
    existing.extend_from_slice(b"end\0");
    // trailing blob
    let new_blob = reb_simulationarchive_blob {
        index: tail_blob.index + 1,
        offset_prev: tail_blob.offset_next,
        offset_next: 0,
    };
    existing.extend_from_slice(&blob_to_bytes(&new_blob));

    if std::fs::write(&filename, &existing).is_err() {
        reb_simulation_error(r, "Can not open file.");
    }
}

/// simulationarchive.c static `check_and_set_simulationarchive_filename`.
fn check_and_set_simulationarchive_filename(r: &mut reb_simulation, filename: &str) -> i32 {
    if filename.is_empty() {
        reb_simulation_error(r, "Filename missing.");
        return -1;
    }
    if std::path::Path::new(filename).exists() {
        reb_simulation_warning(
            r,
            "File in use for Simulationarchive already exists. Snapshots will be appended.",
        );
    }
    r.simulationarchive_filename = Some(filename.to_string());
    0
}

/// simulationarchive.c `reb_simulation_save_to_file_interval`.
pub fn reb_simulation_save_to_file_interval(
    r: &mut reb_simulation,
    filename: &str,
    interval: f64,
) {
    if check_and_set_simulationarchive_filename(r, filename) < 0 {
        return;
    }
    if r.simulationarchive_auto_interval != interval {
        // Only update simulationarchive_next if interval changed.
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r.simulationarchive_auto_interval = interval;
        r.simulationarchive_next = r.t;
    }
}

/// simulationarchive.c `reb_simulation_save_to_file_walltime`.
pub fn reb_simulation_save_to_file_walltime(
    r: &mut reb_simulation,
    filename: &str,
    walltime: f64,
) {
    if check_and_set_simulationarchive_filename(r, filename) < 0 {
        return;
    }
    // Note that this will create two snapshots if restarted.
    r.simulationarchive_auto_walltime = walltime;
    r.simulationarchive_next = r.walltime;
}

/// simulationarchive.c `reb_simulation_save_to_file_step`.
pub fn reb_simulation_save_to_file_step(r: &mut reb_simulation, filename: &str, step: u64) {
    if check_and_set_simulationarchive_filename(r, filename) < 0 {
        return;
    }
    if r.simulationarchive_auto_step != step {
        // Only update simulationarchive_next if interval changed.
        r.simulationarchive_auto_step = step;
        r.simulationarchive_next_step = r.steps_done;
    }
}
