//! output.rs — translation of REBOUNDx output.c
//! Writing the `rebx_extras` state to a REBOUNDx binary file.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx
//! 5.1.0 (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # The format (unchanged from the C)
//!
//! A 64-byte ASCII header, then a nested series of `rebx_binary_field`
//! records. Each record is `{type, size}`; abstract objects (whose
//! length depends on how much the user added) are terminated by an
//! `END` record, and their `size` counts every byte after their own
//! header up to *and including* that `END`, so a reader that does not
//! recognise the type can skip the whole object. Basic data fields
//! carry `size` = the number of payload bytes that follow immediately.
//!
//! ```text
//! SNAPSHOT
//!    REBX_STRUCTURE
//!        REGISTERED_PARAMETERS { REGISTERED_PARAM { PARAM_TYPE NAME END } ... END }
//!        ALLOCATED_FORCES      { FORCE { NAME PARAM_LIST { PARAM ... END } END } ... END }
//!        ALLOCATED_OPERATORS   { OPERATOR { NAME PARAM_LIST { ... } END } ... END }
//!        ADDITIONAL_FORCES     { ADDITIONAL_FORCE { NAME END } ... END }
//!        PRE_TIMESTEP_MODIFICATIONS  { STEP { NAME STEP_DT_FRACTION END } ... END }
//!        POST_TIMESTEP_MODIFICATIONS { STEP { NAME STEP_DT_FRACTION END } ... END }
//!    END
//!    PARTICLES { PARTICLE { PARTICLE_INDEX PARAM_LIST { PARAM ... END } END } ... END }
//! END
//! ```
//!
//! # Byte-format decisions
//!
//! * `struct rebx_binary_field { enum rebx_binary_field_type type;
//!   long size; }` is written with the layout a C compiler for **this
//!   same target** gives it, because that is the byte format of the
//!   `libreboundx` this crate has to interoperate with:
//!   - where C's `long` is 4 bytes (Windows, LLP64 — including the
//!     MSVC reference build used to validate this port): a 4-byte
//!     little-endian `int` immediately followed by a 4-byte
//!     little-endian `long`, **8 bytes** in total;
//!   - where C's `long` is 8 bytes (Linux, macOS, LP64): the same
//!     4-byte `int`, 4 bytes of padding, then an 8-byte little-endian
//!     `long`, **16 bytes** in total.
//!
//!   [`REBX_BINARY_FIELD_SIZE`] is that size, selected by `cfg`. The C
//!   has exactly the same platform dependence — a Windows-built
//!   `libreboundx` cannot read a Linux-written REBOUNDx binary either.
//!   The *reader* in the `input` module is more forgiving than the C
//!   and accepts both dialects (it detects which one a file uses from
//!   its opening `SNAPSHOT` record), so nothing that the C could read
//!   is lost.
//! * Names are written as C strings: `strlen(name) + 1` bytes,
//!   including the terminating NUL.
//! * Payload sizes come from core.c `rebx_sizeof`, exactly as the C:
//!   double = 8, int = 4, vec3d = 24, and **0** for every other type.
//!   REBOUNDx 5.1.0's `rebx_sizeof` has no case for `REBX_TYPE_STRING`,
//!   `REBX_TYPE_ORBIT` or `REBX_TYPE_UINT32`, so the C emits a
//!   zero-length `PARAM_VALUE` for those (and, for the latter two, the
//!   "Need to add new param type to switch statement in rebx_sizeof"
//!   error). That is a limitation of the C release, not of this
//!   translation; reproducing it is what keeps the files byte-identical.
//! * The C keeps the file open and `fseek`s backwards to patch each
//!   object's size once it is known. Here the same bytes are built in
//!   a `Vec<u8>` and patched in place, then written once — the file
//!   contents are identical.
//! * `rebx_githash_str` is deliberately not carried by this crate (see
//!   lib.rs), so the 27 githash bytes of the header are written as
//!   zeros. The C reader compares the header with `strcmp`, which stops
//!   at the NUL that terminates the version string, so it never looks
//!   at them.
//!
//! # Traversal order (load-bearing)
//!
//! `rebx_write_list` walks the C linked list from its **tail to its
//! head**:
//!
//! ```c
//! int N = rebx_len(list);
//! while (N > 0){ current = list; for(i=0;i<N-1;i++) current=current->next; ... N--; }
//! ```
//!
//! Because `rebx_add_node` prepends, the tail is the oldest node, so
//! objects go onto disk in *insertion* order and the reader — which
//! prepends each one it reads — rebuilds the original list exactly.
//! Here a list is a `Vec` whose index 0 is the head, so that traversal
//! is `(0..len).rev()`.
//!
//! The two exceptions are `allocated_forces` and `allocated_operators`:
//! those `Vec`s are *appended* to, because a force's index is its
//! identity (see the `core` module docs), so index 0 is the oldest
//! entry and the C's tail-to-head order is plain forward iteration.

use crate::core::{rebx_error, rebx_extras_mut, rebx_sizeof, rebx_with};
use crate::rebx_version_str;
use crate::types::rebx_param_type::*;
use crate::types::*;
use rebound_rs::{reb_orbit, reb_simulation, reb_vec3d};
use std::io::Write;

/*****************************
 reboundx.h `enum rebx_binary_field_type`
 ****************************/

/// reboundx.h `enum rebx_binary_field_type`. Carried as the C's
/// underlying `int` rather than a Rust enum so that a field type
/// written by a newer REBOUNDx survives a read unchanged (the C stores
/// whatever integer it finds and falls through to `default:`).
pub type rebx_binary_field_type = i32;

pub const REBX_BINARY_FIELD_TYPE_NONE: rebx_binary_field_type = 0;
pub const REBX_BINARY_FIELD_TYPE_OPERATOR: rebx_binary_field_type = 1;
pub const REBX_BINARY_FIELD_TYPE_PARTICLE: rebx_binary_field_type = 2;
pub const REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE: rebx_binary_field_type = 3;
pub const REBX_BINARY_FIELD_TYPE_PARAM: rebx_binary_field_type = 4;
pub const REBX_BINARY_FIELD_TYPE_NAME: rebx_binary_field_type = 5;
pub const REBX_BINARY_FIELD_TYPE_PARAM_TYPE: rebx_binary_field_type = 6;
pub const REBX_BINARY_FIELD_TYPE_PARAM_VALUE: rebx_binary_field_type = 7;
pub const REBX_BINARY_FIELD_TYPE_END: rebx_binary_field_type = 8;
pub const REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX: rebx_binary_field_type = 9;
pub const REBX_BINARY_FIELD_TYPE_REBX_INTEGRATOR: rebx_binary_field_type = 10;
pub const REBX_BINARY_FIELD_TYPE_FORCE_TYPE: rebx_binary_field_type = 11;
pub const REBX_BINARY_FIELD_TYPE_OPERATOR_TYPE: rebx_binary_field_type = 12;
pub const REBX_BINARY_FIELD_TYPE_STEP: rebx_binary_field_type = 13;
pub const REBX_BINARY_FIELD_TYPE_STEP_DT_FRACTION: rebx_binary_field_type = 14;
pub const REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM: rebx_binary_field_type = 15;
pub const REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE: rebx_binary_field_type = 16;
pub const REBX_BINARY_FIELD_TYPE_PARAM_LIST: rebx_binary_field_type = 17;
pub const REBX_BINARY_FIELD_TYPE_REGISTERED_PARAMETERS: rebx_binary_field_type = 18;
pub const REBX_BINARY_FIELD_TYPE_ALLOCATED_FORCES: rebx_binary_field_type = 19;
pub const REBX_BINARY_FIELD_TYPE_ALLOCATED_OPERATORS: rebx_binary_field_type = 20;
pub const REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCES: rebx_binary_field_type = 21;
pub const REBX_BINARY_FIELD_TYPE_PRE_TIMESTEP_MODIFICATIONS: rebx_binary_field_type = 22;
pub const REBX_BINARY_FIELD_TYPE_POST_TIMESTEP_MODIFICATIONS: rebx_binary_field_type = 23;
pub const REBX_BINARY_FIELD_TYPE_PARTICLES: rebx_binary_field_type = 24;
pub const REBX_BINARY_FIELD_TYPE_FORCE: rebx_binary_field_type = 25;
pub const REBX_BINARY_FIELD_TYPE_SNAPSHOT: rebx_binary_field_type = 26;

/// reboundx.h `struct rebx_binary_field` — precedes every object and
/// every data field in a binary file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct rebx_binary_field {
    /// Type of object (C: `enum rebx_binary_field_type type`).
    pub type_: rebx_binary_field_type,
    /// Size in bytes of the object data, not including this structure,
    /// so that it can be skipped (C: `long size`).
    pub size: i64,
}

/// `sizeof(struct rebx_binary_field)` where C's `long` is 4 bytes
/// (Windows, LLP64): the `long` follows the 4-byte enum directly.
pub const REBX_BINARY_FIELD_SIZE_LLP64: usize = 8;

/// `sizeof(struct rebx_binary_field)` where C's `long` is 8 bytes
/// (Linux, macOS, LP64): 4 bytes of padding sit between the enum and
/// the `long`.
pub const REBX_BINARY_FIELD_SIZE_LP64: usize = 16;

/// `sizeof(struct rebx_binary_field)` for **this** target — what a C
/// compiled alongside this crate would use, and therefore what the
/// writer emits. `REBX_BINARY_FIELD_SIZE / 2` is `offsetof(.., size)`
/// in both layouts. See the module docs.
#[cfg(windows)]
pub const REBX_BINARY_FIELD_SIZE: usize = REBX_BINARY_FIELD_SIZE_LLP64;
/// `sizeof(struct rebx_binary_field)` for this target. See above.
#[cfg(not(windows))]
pub const REBX_BINARY_FIELD_SIZE: usize = REBX_BINARY_FIELD_SIZE_LP64;

/// The fixed part of the binary file header (output.c / input.c
/// `const char str[]`).
pub const rebx_binary_header_str: &str = "REBOUNDx Binary File. Version: ";

/// Total size of the binary file header in bytes (output.c writes
/// 62 bytes of text plus two NULs; input.c `fread(readbuf, 1, 64, inf)`).
pub const REBX_BINARY_HEADER_SIZE: usize = 64;

/*****************************
 Low-level byte helpers
 ****************************/

/// C: `strlen(s) + 1` bytes — the string plus its terminating NUL.
///
/// A Rust `String` may contain interior NULs where a C `char*` cannot;
/// such a name would be truncated by the C's `strlen` on the next
/// write. Nothing in REBOUNDx creates one.
pub(crate) fn rebx_cstr_bytes(s: &str) -> Vec<u8> {
    let mut v = Vec::with_capacity(s.len() + 1);
    v.extend_from_slice(s.as_bytes());
    v.push(0);
    v
}

/// The 24 bytes of a `struct reb_vec3d`.
pub(crate) fn rebx_vec3d_bytes(v: &reb_vec3d) -> Vec<u8> {
    let mut d = Vec::with_capacity(24);
    d.extend_from_slice(&v.x.to_le_bytes());
    d.extend_from_slice(&v.y.to_le_bytes());
    d.extend_from_slice(&v.z.to_le_bytes());
    d
}

/// The 216 bytes of a `struct reb_orbit` (21 doubles followed by the
/// `hvec` and `evec` vectors), in rebound.h declaration order.
pub(crate) fn rebx_orbit_bytes(o: &reb_orbit) -> Vec<u8> {
    let mut d = Vec::with_capacity(REBX_ORBIT_RAW_SIZE);
    for x in [
        o.d, o.v, o.h, o.P, o.n, o.a, o.e, o.inc, o.Omega, o.omega, o.pomega, o.f, o.M, o.l,
        o.theta, o.T, o.rhill, o.pal_h, o.pal_k, o.pal_ix, o.pal_iy,
    ] {
        d.extend_from_slice(&x.to_le_bytes());
    }
    d.extend_from_slice(&rebx_vec3d_bytes(&o.hvec));
    d.extend_from_slice(&rebx_vec3d_bytes(&o.evec));
    d
}

/// `sizeof(struct reb_orbit)` on x86-64: 27 doubles.
pub const REBX_ORBIT_RAW_SIZE: usize = 216;

/// Write one `struct rebx_binary_field` (see the module docs for the
/// layout).
fn rebx_push_binary_field(of: &mut Vec<u8>, type_: rebx_binary_field_type, size: i64) {
    of.extend_from_slice(&type_.to_le_bytes());
    let bytes = size.to_le_bytes();
    match REBX_BINARY_FIELD_SIZE {
        // 4-byte `long`: it follows the enum directly. The low four
        // little-endian bytes are the same value narrowed to 32 bits.
        REBX_BINARY_FIELD_SIZE_LLP64 => of.extend_from_slice(&bytes[..4]),
        // 8-byte `long`: 4 bytes of struct padding, then the value.
        _ => {
            of.extend_from_slice(&[0u8; 4]);
            of.extend_from_slice(&bytes);
        }
    }
}

/// output.c `REBX_WRITE_DATA_FIELD` — a field header of `typename`
/// followed by `typesize` bytes of payload.
fn rebx_write_data_field(of: &mut Vec<u8>, type_: rebx_binary_field_type, value: &[u8]) {
    rebx_push_binary_field(of, type_, value.len() as i64);
    of.extend_from_slice(value);
}

/// output.c `REBX_START_OBJECT_FIELD` — writes a header with a
/// placeholder size and returns its position, which
/// `rebx_end_object_field` needs to patch it (the C caches the same
/// position in a local `long pos_start_header_<name>`).
fn rebx_start_object_field(of: &mut Vec<u8>, type_: rebx_binary_field_type) -> usize {
    let pos_start_header = of.len();
    rebx_push_binary_field(of, type_, 0);
    pos_start_header
}

/// output.c `REBX_END_OBJECT_FIELD` — writes the object's `END` field
/// and patches the header with the object's length. As in the C, the
/// recorded size spans everything after the object's own header up to
/// and including that `END` field.
fn rebx_end_object_field(of: &mut Vec<u8>, pos_start_header: usize) {
    rebx_push_binary_field(of, REBX_BINARY_FIELD_TYPE_END, 0);
    let pos_end = of.len();
    let pos_start = pos_start_header + REBX_BINARY_FIELD_SIZE;
    let size = (pos_end - pos_start) as i64;
    // C: `header_<name>.size = ...; fseek(of, pos_start_header_<name>,
    // SEEK_SET); fwrite(&header_<name>, ..)`. offsetof(.., size) is
    // half the struct size in both layouts (4 of 8, or 8 of 16).
    let off = pos_start_header + REBX_BINARY_FIELD_SIZE / 2;
    let bytes = size.to_le_bytes();
    match REBX_BINARY_FIELD_SIZE {
        REBX_BINARY_FIELD_SIZE_LLP64 => of[off..off + 4].copy_from_slice(&bytes[..4]),
        _ => of[off..off + 8].copy_from_slice(&bytes),
    }
}

/*****************************
 Objects
 ****************************/

/// The `size` bytes the C writes from `param->value`.
///
/// The C does `fwrite(param->value, rebx_sizeof(rebx, param->type), 1,
/// of)`, i.e. it copies the first `size` bytes of the malloc'd payload,
/// which is why a type whose `rebx_sizeof` is 0 contributes no bytes at
/// all. `None` stands for the C's `param->value == NULL`, which the C
/// would hand to `fwrite` — see `rebx_write_param`.
fn rebx_param_value_bytes(value: &rebx_param_value, size: usize) -> Option<Vec<u8>> {
    let full: Vec<u8> = match value {
        rebx_param_value::none => return None,
        rebx_param_value::double(v) => v.to_le_bytes().to_vec(),
        rebx_param_value::int(v) => v.to_le_bytes().to_vec(),
        rebx_param_value::uint32(v) => v.to_le_bytes().to_vec(),
        rebx_param_value::vec3d(v) => rebx_vec3d_bytes(v),
        rebx_param_value::orbit(v) => rebx_orbit_bytes(v),
        rebx_param_value::string(s) => rebx_cstr_bytes(s),
        // REBX_TYPE_FORCE never reaches here (rebx_write_force_param
        // handles it); REBX_TYPE_POINTER and REBX_TYPE_ODE are refused
        // by rebx_write_param before this point.
        rebx_param_value::force(_)
        | rebx_param_value::ode(_)
        | rebx_param_value::particles(_)
        | rebx_param_value::particle_index(_) => return None,
    };
    if size <= full.len() {
        Some(full[..size].to_vec())
    } else {
        Some(full)
    }
}

/// output.c `rebx_write_force_param`.
///
/// For `REBX_TYPE_FORCE` parameters the C agrees to store the force's
/// *name* in `PARAM_VALUE`, so that the reader can link the parameter
/// back up with the force it re-created from `ALLOCATED_FORCES`.
fn rebx_write_force_param(rebx: &rebx_extras, param: &rebx_param, of: &mut Vec<u8>) {
    let force_name = match &param.value {
        rebx_param_value::force(idx) => match rebx.allocated_forces.get(*idx) {
            Some(force) => force.name.clone(),
            // C: `struct rebx_force* force = param->value;` on a stale
            // pointer, then `strlen(force->name)`. Nothing sane to
            // write, so the parameter is left out of the file.
            None => return,
        },
        // C: the same dereference on a NULL / non-force value.
        _ => return,
    };
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PARAM);
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_PARAM_TYPE,
        &(param.type_ as i32).to_le_bytes(),
    );
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_NAME,
        &rebx_cstr_bytes(&param.name),
    );
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_PARAM_VALUE,
        &rebx_cstr_bytes(&force_name),
    );
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_param`.
///
/// `param` is a clone rather than the C's borrow into the list, because
/// `rebx_sizeof` needs the extras mutably (it reports unregistered and
/// unhandled types through `rebx_error`) while the parameter still has
/// to be read out of one of its lists.
fn rebx_write_param(rebx: &mut rebx_extras, param: &rebx_param, of: &mut Vec<u8>) {
    if param.type_ == REBX_TYPE_POINTER || param.type_ == REBX_TYPE_ODE {
        // Don't write pointers because we won't know how to load them
        // when we read the binary.
        return;
    }

    if param.type_ == REBX_TYPE_FORCE {
        // Force already written to the allocated_force list.
        rebx_write_force_param(rebx, param, of);
        return;
    }

    let size = rebx_sizeof(rebx, param.type_);
    let value = match rebx_param_value_bytes(&param.value, size) {
        Some(v) => v,
        // C: `param->value == NULL`, which it hands straight to fwrite.
        // Omitting the PARAM_VALUE field instead is exactly what the
        // reader's `param->value == NULL` path expects: the parameter
        // is dropped with REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL.
        None => return,
    };

    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PARAM);
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_PARAM_TYPE,
        &(param.type_ as i32).to_le_bytes(),
    );
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_NAME,
        &rebx_cstr_bytes(&param.name),
    );
    rebx_write_data_field(of, REBX_BINARY_FIELD_TYPE_PARAM_VALUE, &value);
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_registered_param`. Type and name only: a
/// registered parameter is a `rebx_param` with no value.
fn rebx_write_registered_param(param: &rebx_param, of: &mut Vec<u8>) {
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM);
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_PARAM_TYPE,
        &(param.type_ as i32).to_le_bytes(),
    );
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_NAME,
        &rebx_cstr_bytes(&param.name),
    );
    rebx_end_object_field(of, pos);
}

/// output.c's `REBX_WRITE_LIST_FIELD(PARAM_LIST, PARAM, <ap>)`.
///
/// See the module docs for why the `Vec` is walked backwards: the C
/// writes its linked list from tail to head.
fn rebx_write_param_list(rebx: &mut rebx_extras, sel: rebx_ap, of: &mut Vec<u8>) {
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PARAM_LIST);
    let N = rebx.ap(sel).len();
    for i in (0..N).rev() {
        let param = rebx.ap(sel)[i].clone();
        rebx_write_param(rebx, &param, of);
    }
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_force`. Takes the force's index into
/// `allocated_forces` — that index is the force's identity here, where
/// the C passes the `struct rebx_force*`.
fn rebx_write_force(rebx: &mut rebx_extras, force_idx: usize, of: &mut Vec<u8>) {
    let name = match rebx.allocated_forces.get(force_idx) {
        Some(force) => force.name.clone(),
        None => return,
    };
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_FORCE);
    // must write name first so that force can be loaded on read
    rebx_write_data_field(of, REBX_BINARY_FIELD_TYPE_NAME, &rebx_cstr_bytes(&name));
    rebx_write_param_list(rebx, rebx_ap::force(force_idx), of);
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_additional_force`. Same as a force, but only
/// holds the name for later loading, rather than the whole parameter
/// list.
fn rebx_write_additional_force(rebx: &rebx_extras, force_idx: usize, of: &mut Vec<u8>) {
    let name = match rebx.allocated_forces.get(force_idx) {
        Some(force) => force.name.clone(),
        None => return,
    };
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE);
    rebx_write_data_field(of, REBX_BINARY_FIELD_TYPE_NAME, &rebx_cstr_bytes(&name));
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_operator`.
fn rebx_write_operator(rebx: &mut rebx_extras, operator_idx: usize, of: &mut Vec<u8>) {
    let name = match rebx.allocated_operators.get(operator_idx) {
        Some(operator_) => operator_.name.clone(),
        None => return,
    };
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_OPERATOR);
    rebx_write_data_field(of, REBX_BINARY_FIELD_TYPE_NAME, &rebx_cstr_bytes(&name));
    rebx_write_param_list(rebx, rebx_ap::operator_(operator_idx), of);
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_step`.
fn rebx_write_step(rebx: &rebx_extras, step: &rebx_step, of: &mut Vec<u8>) {
    // Need operator name to load it from source when reading it back in
    let name = match rebx.allocated_operators.get(step.operator_) {
        Some(operator_) => operator_.name.clone(),
        None => return,
    };
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_STEP);
    rebx_write_data_field(of, REBX_BINARY_FIELD_TYPE_NAME, &rebx_cstr_bytes(&name));
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_STEP_DT_FRACTION,
        &step.dt_fraction.to_le_bytes(),
    );
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_particle`. The C passes the
/// `struct reb_particle*` and its index; here the parameter list is
/// reached by index alone (`rebx_extras::particle_params`).
fn rebx_write_particle(rebx: &mut rebx_extras, index: i32, of: &mut Vec<u8>) {
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PARTICLE);
    rebx_write_data_field(
        of,
        REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX,
        &index.to_le_bytes(),
    );
    rebx_write_param_list(rebx, rebx_ap::particle(index as usize), of);
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_rebx`.
///
/// The C's generic `rebx_write_list` switch is unrolled here, one loop
/// per list: the Rust lists are typed `Vec`s rather than a single
/// `struct rebx_node*` chain, so there is nothing left to dispatch on.
/// The emitted bytes and their order are unchanged.
fn rebx_write_rebx(rebx: &mut rebx_extras, of: &mut Vec<u8>) {
    let pos_rebx = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE);

    // REBX_WRITE_LIST_FIELD(REGISTERED_PARAMETERS, REGISTERED_PARAM, ..)
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_REGISTERED_PARAMETERS);
        for i in (0..rebx.registered_params.len()).rev() {
            let param = rebx.registered_params[i].clone();
            rebx_write_registered_param(&param, of);
        }
        rebx_end_object_field(of, pos);
    }

    // REBX_WRITE_LIST_FIELD(ALLOCATED_FORCES, FORCE, ..)
    // Forward: this Vec is appended to, so index 0 is the C list's tail.
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_ALLOCATED_FORCES);
        for i in 0..rebx.allocated_forces.len() {
            rebx_write_force(rebx, i, of);
        }
        rebx_end_object_field(of, pos);
    }

    // REBX_WRITE_LIST_FIELD(ALLOCATED_OPERATORS, OPERATOR, ..)
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_ALLOCATED_OPERATORS);
        for i in 0..rebx.allocated_operators.len() {
            rebx_write_operator(rebx, i, of);
        }
        rebx_end_object_field(of, pos);
    }

    // REBX_WRITE_LIST_FIELD(ADDITIONAL_FORCES, ADDITIONAL_FORCE, ..)
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCES);
        for i in (0..rebx.additional_forces.len()).rev() {
            let force_idx = rebx.additional_forces[i];
            rebx_write_additional_force(rebx, force_idx, of);
        }
        rebx_end_object_field(of, pos);
    }

    // REBX_WRITE_LIST_FIELD(PRE_TIMESTEP_MODIFICATIONS, STEP, ..)
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PRE_TIMESTEP_MODIFICATIONS);
        for i in (0..rebx.pre_timestep_modifications.len()).rev() {
            let step = rebx.pre_timestep_modifications[i];
            rebx_write_step(rebx, &step, of);
        }
        rebx_end_object_field(of, pos);
    }

    // REBX_WRITE_LIST_FIELD(POST_TIMESTEP_MODIFICATIONS, STEP, ..)
    {
        let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_POST_TIMESTEP_MODIFICATIONS);
        for i in (0..rebx.post_timestep_modifications.len()).rev() {
            let step = rebx.post_timestep_modifications[i];
            rebx_write_step(rebx, &step, of);
        }
        rebx_end_object_field(of, pos);
    }

    rebx_end_object_field(of, pos_rebx);
}

/// output.c `rebx_write_particles`. Writes a particle field for each
/// particle with a list of its parameters.
fn rebx_write_particles(sim: &reb_simulation, rebx: &mut rebx_extras, of: &mut Vec<u8>) {
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_PARTICLES);
    for i in 0..sim.N {
        rebx_write_particle(rebx, i as i32, of);
    }
    rebx_end_object_field(of, pos);
}

/// output.c `rebx_write_snapshot`. Could be extended to include time or
/// steps_done to make an archive.
fn rebx_write_snapshot(sim: &reb_simulation, rebx: &mut rebx_extras, of: &mut Vec<u8>) {
    let pos = rebx_start_object_field(of, REBX_BINARY_FIELD_TYPE_SNAPSHOT);
    rebx_write_rebx(rebx, of);
    rebx_write_particles(sim, rebx, of);
    rebx_end_object_field(of, pos);
}

/// The 64-byte file header (output.c `rebx_output_binary`).
///
/// `"REBOUNDx Binary File. Version: " + rebx_version_str + '\0'`, then
/// `62 - strlen(str) - strlen(version)` githash bytes, then a final NUL.
fn rebx_write_header(of: &mut Vec<u8>) {
    let str_ = rebx_binary_header_str;
    let lenheader = str_.len() + rebx_version_str.len();
    of.extend_from_slice(str_.as_bytes());
    of.extend_from_slice(rebx_version_str.as_bytes());
    of.push(0);
    // C: fwrite(rebx_githash_str, sizeof(char), 62-lenheader, of).
    // The githash is not carried by this crate and the C reader's
    // strcmp never reaches these bytes; they are written as zeros.
    let n_githash = 62usize.saturating_sub(lenheader);
    of.resize(of.len() + n_githash, 0);
    of.push(0);
}

/// reboundx.h `rebx_output_binary`.
///
/// The C takes `struct rebx_extras*` and reaches the simulation through
/// `rebx->sim`; here the extras live in `sim.extras`, so the simulation
/// is the argument and the extras are taken out of it for the duration
/// (`rebx_with`, see the `core` module docs).
pub fn rebx_output_binary(sim: &mut reb_simulation, filename: &str) {
    // C: FILE* of = fopen(filename,"wb");
    let of = std::fs::File::create(filename);
    if of.is_err() {
        if let Some(rebx) = rebx_extras_mut(sim) {
            rebx_error(
                rebx,
                "REBOUNDx error: Can not open file passed to rebx_output_binary.",
            );
        }
    }

    // C: if (rebx->sim == NULL){ rebx_error(rebx, ""); return; }
    // rebx_with reports the same condition (no REBOUNDx state attached
    // to this simulation) and returns None.
    let buf = match rebx_with(sim, |sim, rebx| {
        let mut buf: Vec<u8> = Vec::new();
        rebx_write_header(&mut buf);
        rebx_write_snapshot(sim, rebx, &mut buf);
        buf
    }) {
        Some(buf) => buf,
        None => return,
    };

    // The C carries on writing to a NULL FILE* when the fopen above
    // failed (undefined behaviour); there is nothing to write to here.
    let mut of = match of {
        Ok(of) => of,
        Err(_) => return,
    };
    let _ = of.write_all(&buf); // C: fwrite(..); fclose(of);
}
