//! input.rs — translation of REBOUNDx input.c
//! Reading a REBOUNDx binary file back into a `rebx_extras`.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx
//! 5.1.0 (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! The byte format, and the reasoning behind every layout decision, is
//! documented in the `output` module; this module is its exact inverse.
//!
//! # How the C's arguments are carried
//!
//! * The C's `FILE* inf` becomes [`rebx_binary_input`], which wraps any
//!   `Read + Seek` stream (a file, a `Cursor` over a buffer, a
//!   Simulationarchive slice) together with the size of the on-disk
//!   `struct rebx_binary_field` this particular file was written with.
//!   That size is 8 or 16 bytes depending on whether the C that wrote
//!   the file had a 4- or 8-byte `long`; the C can only ever read files
//!   written by its own build, while this reader detects which dialect
//!   a file uses and reads either (see
//!   [`rebx_binary_input::rebx_detect_field_size`]).
//! * The C's `struct rebx_node** ap` argument of `rebx_load_list` names
//!   the list that the objects being read must be added to. Here that
//!   is [`rebx_load_target`], which names the same lists by index (the
//!   `rebx_ap` substitution described in the `types` module docs) and
//!   keeps the C's ability to tell `&rebx->pre_timestep_modifications`
//!   from `&rebx->post_timestep_modifications` by identity.
//! * The C reaches the simulation through `rebx->sim`; here the extras
//!   live inside `sim.extras`, so the simulation is the argument and
//!   the extras are borrowed out of it for short stretches
//!   (`rebx_extras_ref` / `rebx_extras_mut`) — the loaders also have to
//!   call `rebx_load_force`, `rebx_add_force` and friends, which need
//!   the simulation itself.
//!
//! # Deliberate departures from the C, all on corrupt input
//!
//! The C reads a field's payload with `fread(malloc(field.size), ...)`
//! and casts it to the type the parameter claims, without checking that
//! the recorded size matches that type; a truncated or hand-edited file
//! therefore reads uninitialised or out-of-bounds memory. Here a
//! payload whose length does not match its type is refused with
//! `REBX_INPUT_BINARY_ERROR_CORRUPT` and the parameter is dropped.
//! Likewise a `PARTICLE_INDEX` outside the simulation's particle array
//! is refused rather than indexed. Well-formed files are unaffected.

use crate::core::{
    rebx_add_force, rebx_add_operator_step, rebx_add_param, rebx_create_param, rebx_error,
    rebx_extras_mut, rebx_extras_ref, rebx_get_force, rebx_get_operator, rebx_initialize,
    rebx_load_force, rebx_load_operator, rebx_with,
};
use crate::output::*;
use crate::rebx_version_str;
use crate::types::rebx_param_type::*;
use crate::types::*;
use rebound_rs::{
    reb_orbit, reb_simulation, reb_simulation_error, reb_simulation_warning, reb_vec3d,
};
use std::io::{Read, Seek, SeekFrom};

/*****************************
 reboundx.h `enum rebx_input_binary_messages`
 ****************************/

/// reboundx.h `enum rebx_input_binary_messages` — a bit field, so it is
/// carried as the underlying integer rather than a Rust enum (the same
/// treatment `rebound_rs` gives `REB_BINARYDATA_ERROR_CODE`).
pub type rebx_input_binary_messages = u32;

pub const REBX_INPUT_BINARY_WARNING_NONE: rebx_input_binary_messages = 0;
pub const REBX_INPUT_BINARY_ERROR_NOFILE: rebx_input_binary_messages = 1;
pub const REBX_INPUT_BINARY_ERROR_CORRUPT: rebx_input_binary_messages = 2;
pub const REBX_INPUT_BINARY_ERROR_NO_MEMORY: rebx_input_binary_messages = 4;
pub const REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED: rebx_input_binary_messages = 8;
pub const REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED: rebx_input_binary_messages = 16;
pub const REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED: rebx_input_binary_messages = 32;
pub const REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED: rebx_input_binary_messages = 64;
pub const REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED: rebx_input_binary_messages = 128;
pub const REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED: rebx_input_binary_messages = 256;
pub const REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED: rebx_input_binary_messages = 512;
pub const REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED: rebx_input_binary_messages = 1024;
pub const REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN: rebx_input_binary_messages = 2048;
pub const REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN: rebx_input_binary_messages = 4096;
pub const REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL: rebx_input_binary_messages = 8192;
pub const REBX_INPUT_BINARY_WARNING_VERSION: rebx_input_binary_messages = 16384;
pub const REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED: rebx_input_binary_messages = 32768;

/// The C's `struct rebx_node** ap` argument: the list that whatever is
/// being read should be added to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_load_target {
    /// C: `NULL`. Forces, operators, additional forces and particles
    /// all know their own destination.
    none,
    /// C: `&force->ap`, `&operator->ap`, `&p->ap`.
    ap(rebx_ap),
    /// C: `&rebx->registered_params`.
    registered_params,
    /// C: `&rebx->pre_timestep_modifications`.
    pre_timestep_modifications,
    /// C: `&rebx->post_timestep_modifications`.
    post_timestep_modifications,
}

/*****************************
 The input stream (C: `FILE* inf`)
 ****************************/

/// The C's `FILE* inf`, plus the on-disk size of one
/// `struct rebx_binary_field` in the file being read.
///
/// See the `output` module docs: that size is 8 bytes where the writing
/// C build had a 4-byte `long` and 16 where it had an 8-byte one. It
/// starts as [`REBX_BINARY_FIELD_SIZE`] — what a C compiled for this
/// same target would use — and is corrected by
/// [`rebx_binary_input::rebx_detect_field_size`] when the header is read.
pub struct rebx_binary_input<R: Read + Seek> {
    /// The stream, positioned where the C's file pointer would be.
    pub inf: R,
    /// `sizeof(struct rebx_binary_field)` in this file.
    pub field_size: usize,
}

impl<R: Read + Seek> rebx_binary_input<R> {
    /// Wrap a stream, assuming for now the layout this target's C would
    /// write.
    pub fn new(inf: R) -> Self {
        rebx_binary_input {
            inf,
            field_size: REBX_BINARY_FIELD_SIZE,
        }
    }

    /// Work out which of the two `struct rebx_binary_field` layouts the
    /// file uses, without moving the stream position.
    ///
    /// Call with the stream just past the 64-byte header, where the
    /// opening `SNAPSHOT` record sits. In both layouts the first four
    /// bytes are the type. The next four are either the whole 4-byte
    /// `long` size (LLP64) or the struct's padding (LP64) — and a
    /// `SNAPSHOT` always has content, so its size is never 0. Zero
    /// there therefore means padding, i.e. the 16-byte layout.
    ///
    /// A file that does not start with a `SNAPSHOT` is left on this
    /// target's layout and will be reported corrupt by the loaders,
    /// exactly as the C reports it.
    pub fn rebx_detect_field_size(&mut self) {
        let pos = match self.inf.stream_position() {
            Ok(pos) => pos,
            Err(_) => return,
        };
        let mut b = [0u8; 8];
        let ok = self.inf.read_exact(&mut b).is_ok();
        if self.inf.seek(SeekFrom::Start(pos)).is_err() {
            return;
        }
        if !ok {
            return;
        }
        let type_ = rebx_read_i32_le(&b, 0);
        let tail = rebx_read_i32_le(&b, 4);
        if type_ == REBX_BINARY_FIELD_TYPE_SNAPSHOT {
            self.field_size = if tail == 0 {
                REBX_BINARY_FIELD_SIZE_LP64
            } else {
                REBX_BINARY_FIELD_SIZE_LLP64
            };
        }
    }

    /// C: `fread(&field, sizeof(field), 1, inf)`. `None` when the
    /// record could not be read in full, which is the C's return of 0.
    fn rebx_fread_binary_field(&mut self) -> Option<rebx_binary_field> {
        let mut b = [0u8; REBX_BINARY_FIELD_SIZE_LP64];
        let n = self.field_size;
        if self.inf.read_exact(&mut b[..n]).is_err() {
            return None;
        }
        let type_ = rebx_read_i32_le(&b, 0);
        // offsetof(struct rebx_binary_field, size) is n/2 in both
        // layouts; the size itself is a 4- or 8-byte `long`.
        let size = if n == REBX_BINARY_FIELD_SIZE_LLP64 {
            rebx_read_i32_le(&b, 4) as i64
        } else {
            let mut s = [0u8; 8];
            s.copy_from_slice(&b[8..16]);
            i64::from_le_bytes(s)
        };
        Some(rebx_binary_field { type_, size })
    }

    /// C: `fread(buf, size, 1, inf)`. Returns `None` exactly where the
    /// C's `fread` returns 0 — including for `size == 0`, since `fread`
    /// counts complete objects read and a zero-sized object is never
    /// "read".
    ///
    /// The bytes are pulled through a `take` so that a corrupt
    /// (absurdly large) size cannot make this allocate more than the
    /// stream actually holds.
    fn rebx_fread_bytes(&mut self, size: i64) -> Option<Vec<u8>> {
        if size <= 0 {
            return None;
        }
        let mut v: Vec<u8> = Vec::new();
        if self
            .inf
            .by_ref()
            .take(size as u64)
            .read_to_end(&mut v)
            .is_err()
        {
            return None;
        }
        if v.len() as i64 != size {
            return None;
        }
        Some(v)
    }

    /// C: `fseek(inf, field_size, SEEK_CUR)`.
    fn rebx_skip(&mut self, field_size: i64) {
        let _ = self.inf.seek(SeekFrom::Current(field_size));
    }
}

/*****************************
 Low-level byte helpers
 ****************************/

/// The C string in `d` as a `String` (up to the first NUL). Bytes that
/// are not valid UTF-8 are replaced, as REBOUNDx names are ASCII.
fn rebx_string_from_bytes(d: &[u8]) -> String {
    let end = match d.iter().position(|b| *b == 0) {
        Some(i) => i,
        None => d.len(),
    };
    String::from_utf8_lossy(&d[..end]).into_owned()
}

fn rebx_read_f64_le(d: &[u8], off: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[off..off + 8]);
    f64::from_le_bytes(b)
}

fn rebx_read_i32_le(d: &[u8], off: usize) -> i32 {
    let mut b = [0u8; 4];
    b.copy_from_slice(&d[off..off + 4]);
    i32::from_le_bytes(b)
}

fn rebx_read_u32_le(d: &[u8], off: usize) -> u32 {
    let mut b = [0u8; 4];
    b.copy_from_slice(&d[off..off + 4]);
    u32::from_le_bytes(b)
}

fn rebx_vec3d_from_bytes(d: &[u8], off: usize) -> reb_vec3d {
    reb_vec3d {
        x: rebx_read_f64_le(d, off),
        y: rebx_read_f64_le(d, off + 8),
        z: rebx_read_f64_le(d, off + 16),
    }
}

/// The inverse of `output::rebx_orbit_bytes`.
fn rebx_orbit_from_bytes(d: &[u8]) -> reb_orbit {
    let f = |i: usize| rebx_read_f64_le(d, 8 * i);
    reb_orbit {
        d: f(0),
        v: f(1),
        h: f(2),
        P: f(3),
        n: f(4),
        a: f(5),
        e: f(6),
        inc: f(7),
        Omega: f(8),
        omega: f(9),
        pomega: f(10),
        f: f(11),
        M: f(12),
        l: f(13),
        theta: f(14),
        T: f(15),
        rhill: f(16),
        pal_h: f(17),
        pal_k: f(18),
        pal_ix: f(19),
        pal_iy: f(20),
        hvec: rebx_vec3d_from_bytes(d, 8 * 21),
        evec: rebx_vec3d_from_bytes(d, 8 * 24),
    }
}

/// The stored `enum rebx_param_type` as the Rust enum. An integer this
/// version does not know becomes `REBX_TYPE_NONE`, which is what the
/// C's own "type is not a type we handle" checks then reject.
fn rebx_param_type_from_i32(v: i32) -> rebx_param_type {
    match v {
        x if x == REBX_TYPE_DOUBLE as i32 => REBX_TYPE_DOUBLE,
        x if x == REBX_TYPE_INT as i32 => REBX_TYPE_INT,
        x if x == REBX_TYPE_POINTER as i32 => REBX_TYPE_POINTER,
        x if x == REBX_TYPE_FORCE as i32 => REBX_TYPE_FORCE,
        x if x == REBX_TYPE_UINT32 as i32 => REBX_TYPE_UINT32,
        x if x == REBX_TYPE_ORBIT as i32 => REBX_TYPE_ORBIT,
        x if x == REBX_TYPE_ODE as i32 => REBX_TYPE_ODE,
        x if x == REBX_TYPE_VEC3D as i32 => REBX_TYPE_VEC3D,
        x if x == REBX_TYPE_STRING as i32 => REBX_TYPE_STRING,
        _ => REBX_TYPE_NONE,
    }
}

/// The stored `PARAM_VALUE` payload as a typed value — the C's cast of
/// the malloc'd buffer to the type the parameter claims.
///
/// `None` when the payload's length does not match the type (the C does
/// not check; see the module docs). `REBX_TYPE_FORCE` is handled by the
/// caller, which resolves the stored force name against
/// `allocated_forces`.
fn rebx_param_value_from_bytes(type_: rebx_param_type, d: &[u8]) -> Option<rebx_param_value> {
    match type_ {
        REBX_TYPE_DOUBLE => {
            if d.len() != 8 {
                return None;
            }
            Some(rebx_param_value::double(rebx_read_f64_le(d, 0)))
        }
        REBX_TYPE_INT => {
            if d.len() != 4 {
                return None;
            }
            Some(rebx_param_value::int(rebx_read_i32_le(d, 0)))
        }
        REBX_TYPE_UINT32 => {
            if d.len() != 4 {
                return None;
            }
            Some(rebx_param_value::uint32(rebx_read_u32_le(d, 0)))
        }
        REBX_TYPE_VEC3D => {
            if d.len() != 24 {
                return None;
            }
            Some(rebx_param_value::vec3d(rebx_vec3d_from_bytes(d, 0)))
        }
        REBX_TYPE_ORBIT => {
            if d.len() != REBX_ORBIT_RAW_SIZE {
                return None;
            }
            Some(rebx_param_value::orbit(rebx_orbit_from_bytes(d)))
        }
        REBX_TYPE_STRING => Some(rebx_param_value::string(rebx_string_from_bytes(d))),
        // REBX_TYPE_FORCE is resolved by the caller; POINTER and ODE are
        // never written (output.c `rebx_write_param` refuses them) and
        // could not be reconstructed from bytes anyway.
        _ => None,
    }
}

/*****************************
 Reading objects
 ****************************/

/// What `rebx_read_param` recovered from a `PARAM` object: the C's
/// `struct rebx_param` before its `void* value` has been cast.
struct rebx_read_param_out {
    type_: rebx_param_type,
    name: String,
    /// C: `param->value`, still the raw malloc'd bytes. `None` is the
    /// C's `param->value == NULL`.
    value: Option<Vec<u8>>,
}

/// input.c `rebx_read_param`.
fn rebx_read_param<R: Read + Seek>(
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> Option<rebx_read_param_out> {
    let mut type_ = REBX_TYPE_NONE;
    let mut name: Option<String> = None;
    let mut value: Option<Vec<u8>> = None;

    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                // means we didn't reach an END field. Corrupt
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                break;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_PARAM_TYPE => match inf.rebx_fread_bytes(field.size) {
                Some(d) if d.len() == 4 => type_ = rebx_param_type_from_i32(rebx_read_i32_le(&d, 0)),
                _ => *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT,
            },
            REBX_BINARY_FIELD_TYPE_NAME => match inf.rebx_fread_bytes(field.size) {
                Some(d) => name = Some(rebx_string_from_bytes(&d)),
                None => *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT,
            },
            REBX_BINARY_FIELD_TYPE_PARAM_VALUE => match inf.rebx_fread_bytes(field.size) {
                Some(d) => value = Some(d),
                None => *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT,
            },
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            // Might have added new fields, saved with a new version and
            // loaded with an old version. Note the C does *not* skip the
            // unknown field's payload here (unlike its other loops), so
            // neither do we.
            _ => *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN,
        }
    }

    // Check type and name after param has been loaded. Check value later
    // (registered params should have value=NULL)
    let name = match name {
        Some(name) => name,
        None => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return None;
        }
    };
    if type_ == REBX_TYPE_NONE {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return None;
    }
    Some(rebx_read_param_out { type_, name, value })
}

/// input.c `rebx_load_param`.
fn rebx_load_param<R: Read + Seek>(
    sim: &mut reb_simulation,
    target: rebx_load_target,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let read = match rebx_read_param(inf, warnings) {
        Some(read) => read,
        None => return 0,
    };

    let data = match read.value {
        Some(data) => data,
        None => {
            *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL;
            return 0;
        }
    };

    let value = if read.type_ == REBX_TYPE_FORCE {
        // The force's name was stored in param->value; look up the force
        // that ALLOCATED_FORCES has already re-created.
        let force_name = rebx_string_from_bytes(&data);
        let force = match rebx_extras_ref(sim) {
            Some(rebx) => rebx_get_force(rebx, &force_name),
            None => None,
        };
        match force {
            Some(force) => rebx_param_value::force(force),
            None => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED;
                return 0;
            }
        }
    } else {
        match rebx_param_value_from_bytes(read.type_, &data) {
            Some(value) => value,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                return 0;
            }
        }
    };

    let mut param = rebx_create_param(&read.name, read.type_);
    param.value = value;

    let rebx = match rebx_extras_mut(sim) {
        Some(rebx) => rebx,
        None => return 0,
    };
    match target {
        rebx_load_target::ap(sel) => match rebx.ap_mut(sel) {
            Some(ap) => {
                rebx_add_param(ap, param);
                1
            }
            None => 0,
        },
        rebx_load_target::registered_params => {
            rebx_add_param(&mut rebx.registered_params, param);
            1
        }
        // C: rebx_add_param with a NULL apptr, which it refuses.
        _ => 0,
    }
}

/// input.c `rebx_load_registered_param`.
fn rebx_load_registered_param<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let read = match rebx_read_param(inf, warnings) {
        Some(read) => read,
        None => return 0,
    };

    // C: the param is added exactly as read, value included — which for
    // a registered parameter is NULL, since output.c writes only its
    // type and name.
    let mut param = rebx_create_param(&read.name, read.type_);
    if let Some(data) = &read.value {
        if let Some(value) = rebx_param_value_from_bytes(read.type_, data) {
            param.value = value;
        }
    }

    let rebx = match rebx_extras_mut(sim) {
        Some(rebx) => rebx,
        None => return 0,
    };
    rebx_add_param(&mut rebx.registered_params, param);
    1
}

/// input.c `rebx_load_name`.
fn rebx_load_name<R: Read + Seek>(
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> Option<String> {
    let field = match inf.rebx_fread_binary_field() {
        Some(field) => field,
        None => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return None;
        }
    };
    if field.type_ != REBX_BINARY_FIELD_TYPE_NAME {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return None;
    }
    match inf.rebx_fread_bytes(field.size) {
        Some(d) => Some(rebx_string_from_bytes(&d)),
        None => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            None
        }
    }
}

/// input.c `rebx_load_force_field`.
fn rebx_load_force_field<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    // Name of force always comes first so that we can load it
    let name = match rebx_load_name(inf, warnings) {
        Some(name) => name,
        None => return 0,
    };
    let force = match rebx_load_force(sim, &name) {
        Some(force) => force,
        None => {
            *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED;
            return 0;
        }
    };

    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                return 0;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_PARAM_LIST => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_PARAM,
                    rebx_load_target::ap(rebx_ap::force(force)),
                    inf,
                    warnings,
                ) == 0
                {
                    return 0;
                }
            }
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }

    1
}

/// input.c `rebx_load_additional_force_field`. The force is already
/// loaded in `allocated_forces`; get it from that list and add it to the
/// simulation.
fn rebx_load_additional_force_field<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let name = match rebx_load_name(inf, warnings) {
        Some(name) => name,
        None => return 0,
    };
    let force = match rebx_extras_ref(sim) {
        Some(rebx) => rebx_get_force(rebx, &name),
        None => None,
    };
    let force = match force {
        Some(force) => force,
        None => return 0,
    };

    // Just catches END for now. This makes it flexible to addition of fields
    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                return 0;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }

    rebx_add_force(sim, force) // add to additional_forces
}

/// input.c `rebx_load_operator_field`.
fn rebx_load_operator_field<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    // Name of operator always comes first so that we can load it
    let name = match rebx_load_name(inf, warnings) {
        Some(name) => name,
        None => return 0,
    };
    let operator_ = match rebx_load_operator(sim, &name) {
        Some(operator_) => operator_,
        None => {
            *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
            return 0;
        }
    };

    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                break;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_PARAM_LIST => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_PARAM,
                    rebx_load_target::ap(rebx_ap::operator_(operator_)),
                    inf,
                    warnings,
                ) == 0
                {
                    return 0;
                }
            }
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }

    1
}

/// input.c `rebx_load_step_field`.
fn rebx_load_step_field<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
    target: rebx_load_target,
) -> i32 {
    let name = match rebx_load_name(inf, warnings) {
        Some(name) => name,
        None => return 0,
    };
    let operator_ = match rebx_extras_ref(sim) {
        Some(rebx) => rebx_get_operator(rebx, &name),
        None => None,
    };
    let operator_ = match operator_ {
        Some(operator_) => operator_,
        None => {
            *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
            return 0;
        }
    };

    let mut dt_fraction = 0.;
    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                break;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_STEP_DT_FRACTION => match inf.rebx_fread_bytes(field.size) {
                Some(d) if d.len() == 8 => dt_fraction = rebx_read_f64_le(&d, 0),
                _ => *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT,
            },
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            // As in rebx_read_param, the C does not skip here.
            _ => *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN,
        }
    }

    if dt_fraction == 0. {
        return 0;
    }

    // The C compares the `ap` pointer it was handed against
    // `&rebx->pre_timestep_modifications` and
    // `&rebx->post_timestep_modifications`; rebx_load_target names the
    // same two lists.
    match target {
        rebx_load_target::pre_timestep_modifications => {
            rebx_add_operator_step(sim, operator_, dt_fraction, rebx_timing::REBX_TIMING_PRE)
        }
        rebx_load_target::post_timestep_modifications => {
            rebx_add_operator_step(sim, operator_, dt_fraction, rebx_timing::REBX_TIMING_POST)
        }
        _ => 0,
    }
}

/// input.c `rebx_load_particle`.
fn rebx_load_particle<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let field = match inf.rebx_fread_binary_field() {
        Some(field) => field,
        None => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
    };
    if field.type_ != REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    let index = match inf.rebx_fread_bytes(field.size) {
        Some(d) if d.len() == 4 => rebx_read_i32_le(&d, 0),
        _ => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
    };

    // C: p = &rebx->sim->particles[index], with no range check. Here the
    // parameter list is `rebx.particle_params[index]`, and attaching
    // parameters to a particle the simulation does not have would be
    // silently useless, so the index is checked.
    if index < 0 || index as usize >= sim.N {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }
    let index = index as usize;

    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                return 0;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_PARAM_LIST => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_PARAM,
                    rebx_load_target::ap(rebx_ap::particle(index)),
                    inf,
                    warnings,
                ) == 0
                {
                    return 0;
                }
            }
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }
    1
}

/// input.c `rebx_load_rebx`.
fn rebx_load_rebx<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                break;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            REBX_BINARY_FIELD_TYPE_REGISTERED_PARAMETERS => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM,
                    rebx_load_target::registered_params,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_ALLOCATED_FORCES => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_FORCE,
                    rebx_load_target::none,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_ALLOCATED_OPERATORS => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_OPERATOR,
                    rebx_load_target::none,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCES => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE,
                    rebx_load_target::none,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_PRE_TIMESTEP_MODIFICATIONS => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_STEP,
                    rebx_load_target::pre_timestep_modifications,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_POST_TIMESTEP_MODIFICATIONS => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_STEP,
                    rebx_load_target::post_timestep_modifications,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }
    1
}

/// input.c `rebx_load_snapshot`.
fn rebx_load_snapshot<R: Read + Seek>(
    sim: &mut reb_simulation,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    let field = match inf.rebx_fread_binary_field() {
        Some(field) => field,
        None => {
            *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
            return 0;
        }
    };
    if field.type_ != REBX_BINARY_FIELD_TYPE_SNAPSHOT {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return 0;
    }

    let mut reading_fields = true;
    while reading_fields {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => {
                *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                break;
            }
        };
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE => {
                if rebx_load_rebx(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_PARTICLES => {
                if rebx_load_list(
                    sim,
                    REBX_BINARY_FIELD_TYPE_PARTICLE,
                    rebx_load_target::none,
                    inf,
                    warnings,
                ) == 0
                {
                    *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_END => reading_fields = false,
            _ => {
                *warnings |= REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN;
                inf.rebx_skip(field.size);
            }
        }
    }

    1
}

/// input.c `rebx_load_list`. Only fails (returns 0) if the binary is in
/// the wrong format.
fn rebx_load_list<R: Read + Seek>(
    sim: &mut reb_simulation,
    expected_type: rebx_binary_field_type,
    target: rebx_load_target,
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) -> i32 {
    loop {
        let field = match inf.rebx_fread_binary_field() {
            Some(field) => field,
            None => return 0,
        };

        // Check whether we've reached end before checking for expected type
        if field.type_ == REBX_BINARY_FIELD_TYPE_END {
            break;
        }

        if field.type_ != expected_type {
            return 0;
        }

        // Only will have fields of expected_type, check function to call
        match field.type_ {
            REBX_BINARY_FIELD_TYPE_PARAM => {
                if rebx_load_param(sim, target, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM => {
                if rebx_load_registered_param(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_FORCE => {
                if rebx_load_force_field(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE => {
                if rebx_load_additional_force_field(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_OPERATOR => {
                if rebx_load_operator_field(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_STEP => {
                if rebx_load_step_field(sim, inf, warnings, target) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            REBX_BINARY_FIELD_TYPE_PARTICLE => {
                if rebx_load_particle(sim, inf, warnings) == 0 {
                    *warnings |= REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED;
                    inf.rebx_skip(field.size);
                }
            }
            _ => {
                if let Some(rebx) = rebx_extras_mut(sim) {
                    rebx_error(rebx, "REBOUNDx Error. Reached default in rebx_load_list reading binary. Should never reach this case. Means we added a list to rebx and didn't add new case to load_list. Please report bug as Github issue.\n");
                }
                return 0;
            }
        }
    }
    1
}

/// input.c `rebx_input_read_header`.
///
/// The C builds the header it would have written and `strcmp`s the two,
/// which stops at the NUL that follows the version string: the githash
/// is never compared. This additionally works out which
/// `struct rebx_binary_field` layout the rest of the file uses.
fn rebx_input_read_header<R: Read + Seek>(
    inf: &mut rebx_binary_input<R>,
    warnings: &mut rebx_input_binary_messages,
) {
    let mut readbuf = [0u8; REBX_BINARY_HEADER_SIZE];
    if inf.inf.read_exact(&mut readbuf).is_err() {
        *warnings |= REBX_INPUT_BINARY_ERROR_CORRUPT;
        return;
    }
    inf.rebx_detect_field_size();
    let curvbuf = format!("{}{}", rebx_binary_header_str, rebx_version_str);
    let end = match readbuf.iter().position(|b| *b == 0) {
        Some(i) => i,
        None => readbuf.len(),
    };
    if readbuf[..end] != *curvbuf.as_bytes() {
        *warnings |= REBX_INPUT_BINARY_WARNING_VERSION;
    }
}

/*****************************
 Public interface
 ****************************/

/// reboundx.h `rebx_input_skip_binary_field`.
pub fn rebx_input_skip_binary_field<R: Read + Seek>(
    inf: &mut rebx_binary_input<R>,
    field_size: i64,
) {
    inf.rebx_skip(field_size);
}

/// reboundx.h `rebx_input_read_binary_field`. Returns a zeroed field
/// (C: `struct rebx_binary_field empty = {0}`) if the read fails.
pub fn rebx_input_read_binary_field<R: Read + Seek>(
    inf: &mut rebx_binary_input<R>,
) -> rebx_binary_field {
    inf.rebx_fread_binary_field().unwrap_or_default()
}

/// reboundx.h `rebx_input_inspect_binary`. Opens the file, reads the
/// header, and hands back the stream positioned just after it.
///
/// The C returns a `FILE*` (NULL on failure); here it is an
/// `Option<rebx_binary_input<std::fs::File>>`, which also carries the
/// field layout detected from the header.
pub fn rebx_input_inspect_binary(
    filename: &str,
    warnings: &mut rebx_input_binary_messages,
) -> Option<rebx_binary_input<std::fs::File>> {
    let file = match std::fs::File::open(filename) {
        Ok(file) => file,
        Err(_) => {
            *warnings |= REBX_INPUT_BINARY_ERROR_NOFILE;
            return None;
        }
    };
    let mut inf = rebx_binary_input::new(file);
    rebx_input_read_header(&mut inf, warnings);
    Some(inf)
}

/// reboundx.h `rebx_init_extras_from_binary`. Loads a binary into the
/// REBOUNDx state already attached to `sim`.
///
/// The C takes `struct rebx_extras*` and reaches the simulation through
/// `rebx->sim`; here the extras live in `sim.extras`.
pub fn rebx_init_extras_from_binary(
    sim: &mut reb_simulation,
    filename: &str,
    warnings: &mut rebx_input_binary_messages,
) {
    // C: if (rebx->sim == NULL){ rebx_error(rebx, ""); return; }
    // rebx_with reports exactly that condition (no REBOUNDx state on
    // this simulation) and returns None without calling the closure.
    if rebx_with(sim, |_sim, _rebx| ()).is_none() {
        return;
    }

    let file = match std::fs::File::open(filename) {
        Ok(file) => file,
        Err(_) => {
            *warnings |= REBX_INPUT_BINARY_ERROR_NOFILE;
            return;
        }
    };
    let mut inf = rebx_binary_input::new(file);

    rebx_input_read_header(&mut inf, warnings);
    rebx_load_snapshot(sim, &mut inf, warnings);
}

/// reboundx.h `rebx_create_extras_from_binary`. Attaches a fresh
/// REBOUNDx state to `sim` and fills it from the binary file.
///
/// The C returns the new `struct rebx_extras*`; here the state lives in
/// `sim.extras` and is reached with `rebx_extras_mut` / `rebx_extras_ref`,
/// so there is nothing to hand back. As in the C, the state is created
/// with `rebx_initialize` rather than `rebx_attach`, so that the default
/// registered parameters are *not* registered — the file carries the
/// registered-parameter list that was in effect when it was written.
pub fn rebx_create_extras_from_binary(sim: &mut reb_simulation, filename: &str) {
    let mut warnings: rebx_input_binary_messages = REBX_INPUT_BINARY_WARNING_NONE;
    // create manually so that default registered parameters not loaded
    rebx_initialize(sim);
    rebx_init_extras_from_binary(sim, filename, &mut warnings);

    if warnings & REBX_INPUT_BINARY_ERROR_NOFILE != 0 {
        reb_simulation_error(sim, "REBOUNDx: Cannot open binary file. Check filename.");
    }
    if warnings & REBX_INPUT_BINARY_ERROR_CORRUPT != 0 {
        reb_simulation_error(sim, "REBOUNDx: Binary file is unreadable. Please open an issue on Github mentioning the version of REBOUND and REBOUNDx you are using and include the binary file.");
    }
    if warnings & REBX_INPUT_BINARY_ERROR_NO_MEMORY != 0 {
        reb_simulation_error(sim, "REBOUNDx: Ran out of system memory.");
    }
    if warnings & REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED != 0 {
        reb_simulation_error(sim, "REBOUNDx: REBOUNDx structure couldn't be loaded.");
    }
    if warnings & REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED != 0 {
        reb_simulation_error(sim, "REBOUNDx: At least one registered parameter was not loaded. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED != 0 {
        reb_simulation_warning(sim, "REBOUNDx: At least one force or operator parameter was not loaded from the binary file. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED != 0 {
        reb_simulation_warning(
            sim,
            "REBOUNDx: At least one particle's parameters were not loaded from the binary file.",
        );
    }
    if warnings & REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED != 0 {
        reb_simulation_warning(sim, "REBOUNDx: At least one force was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED != 0 {
        reb_simulation_warning(sim, "REBOUNDx: At least one operator was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED != 0 {
        reb_simulation_warning(
            sim,
            "REBOUNDx: At least one operator step was not loaded from the binary file.",
        );
    }
    if warnings & REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED != 0 {
        reb_simulation_warning(sim, "REBOUNDx: At least one force was not added to the simulation. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Unknown list in the REBOUNDx structure wasn't loaded. This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL != 0 {
        reb_simulation_warning(sim, "REBOUNDx: The value of at least one parameter was not loaded. This can happen if a custom structure was added by the user as a parameter. See Parameters.ipynb jupyter notebook example.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_VERSION != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed. Check that effects and parameters are loaded as expected.");
    }
    if warnings & REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED != 0 {
        reb_simulation_warning(sim, "REBOUNDx: A force parameter failed to load from the list of REBOUNDx implemented forces. Custom forces can't be saved to a REBOUNDx binary, and function points must be reset when a simulation is reloaded.");
    }
}
