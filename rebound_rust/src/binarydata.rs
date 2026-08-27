//! binarydata.rs — output, input and comparison of simulations in
//! binary format (from binarydata.c/h; (c) 2026 Hanno Rein).
//!
//! The C implementation is driven by `offsetof`-based field descriptor
//! tables and reads/writes raw struct memory. Rust cannot access
//! struct memory through offsets without `unsafe`, so this module
//! reproduces the SAME byte format with explicit per-field
//! serializers, in the same field order as the C descriptor tables.
//!
//! File format (identical to the C):
//! - 64-byte header ("REBOUND Binary File. Version: <v>\0<githash>...").
//! - Sequence of fields: 16-byte `reb_binarydata_field` (two u64 LE:
//!   size of name incl. NUL, size of data), the name bytes, the data.
//! - `struct reb_particle` payloads use the 112-byte x86-64 memory
//!   layout (11 doubles + 3 pointers). The Rust writer stores 0 for
//!   the `ap`/`sim` pointers, and a synthetic non-zero id in place of
//!   the `name` pointer (the C stores the actual heap pointer, which
//!   is only ever compared for equality against the pointers stored in
//!   `name_list` — the synthetic ids reproduce that protocol exactly,
//!   in both directions, including for files written by the C build).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use crate::simulation::reb_simulation_set_integrator;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;
use crate::{reb_githash_str, reb_version_str};
use std::io::{Read, Seek, SeekFrom};

/// binarydata.c `reb_binarydata_header` — corresponds to the first few
/// ASCII characters ("REBOUND ") in a binary file.
pub const reb_binarydata_header: u64 = 0x20444E554F424552;

// (rebound.h `REB_STRING_SIZE_MAX` is defined in types.rs.)

// binarydata.h `enum REB_BINARYDATA_ERROR_CODE` (bit flags).
pub type REB_BINARYDATA_ERROR_CODE = u32;
pub const REB_BINARYDATA_WARNING_NONE: u32 = 0;
pub const REB_BINARYDATA_ERROR_NOFILE: u32 = 1;
pub const REB_BINARYDATA_WARNING_VERSION: u32 = 2;
pub const REB_BINARYDATA_WARNING_POINTERS: u32 = 4;
pub const REB_BINARYDATA_WARNING_PARTICLES: u32 = 8;
pub const REB_BINARYDATA_ERROR_FILENOTOPEN: u32 = 16;
pub const REB_BINARYDATA_ERROR_OUTOFRANGE: u32 = 32;
pub const REB_BINARYDATA_ERROR_SEEK: u32 = 64;
pub const REB_BINARYDATA_WARNING_FIELD_UNKNOWN: u32 = 128;
pub const REB_BINARYDATA_ERROR_INTEGRATOR: u32 = 256;
pub const REB_BINARYDATA_WARNING_CORRUPTFILE: u32 = 512;
pub const REB_BINARYDATA_ERROR_OLD: u32 = 1024;
pub const REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR: u32 = 2048;

// binarydata.h `enum REB_BINARYDATA_OUTPUT`.
pub type REB_BINARYDATA_OUTPUT = i32;
pub const REB_BINARYDATA_OUTPUT_NONE: i32 = 0;
pub const REB_BINARYDATA_OUTPUT_PRINT: i32 = 1;
pub const REB_BINARYDATA_OUTPUT_STREAM: i32 = 2;
pub const REB_BINARYDATA_OUTPUT_BUFFER: i32 = 3;

/// binarydata.h `struct reb_binarydata_field` — precedes the actual
/// data of every field in a binary file.
#[derive(Clone, Copy, Debug, Default)]
pub struct reb_binarydata_field {
    /// Size of the name of the field including the \0 character.
    pub size_name: u64,
    /// Size of the data in the field.
    pub size_data: u64,
}

/// Size of the on-disk C structs (x86-64).
pub const REB_BINARYDATA_FIELD_SIZE: usize = 16;
/// sizeof(struct reb_particle) on x86-64 (11 doubles + name/ap/sim).
pub const REB_PARTICLE_RAW_SIZE: usize = 112;
/// sizeof(struct reb_variational_configuration) on x86-64.
pub const REB_VAR_CONFIG_RAW_SIZE: usize = 40;
/// sizeof(struct reb_simulationarchive_blob) (3 x int32).
pub const REB_SA_BLOB_SIZE: usize = 12;
/// sizeof(struct reb_particle_int) (6 x int64), integrator_janus.h.
pub const REB_PARTICLE_INT_RAW_SIZE: usize = 48;

// ---------------------------------------------------------------------
// Low-level byte helpers
// ---------------------------------------------------------------------

fn push_bytes(buf: &mut Vec<u8>, data: &[u8]) {
    buf.extend_from_slice(data);
}

/// Emit one field: header struct, name (with NUL), data.
fn push_field(buf: &mut Vec<u8>, name: &str, data: &[u8]) {
    let size_name = (name.len() + 1) as u64;
    let size_data = data.len() as u64;
    push_bytes(buf, &size_name.to_le_bytes());
    push_bytes(buf, &size_data.to_le_bytes());
    push_bytes(buf, name.as_bytes());
    buf.push(0);
    push_bytes(buf, data);
}

fn push_field_f64(buf: &mut Vec<u8>, name: &str, v: f64) {
    push_field(buf, name, &v.to_le_bytes());
}
fn push_field_i32(buf: &mut Vec<u8>, name: &str, v: i32) {
    push_field(buf, name, &v.to_le_bytes());
}
fn push_field_u32(buf: &mut Vec<u8>, name: &str, v: u32) {
    push_field(buf, name, &v.to_le_bytes());
}
fn push_field_i64(buf: &mut Vec<u8>, name: &str, v: i64) {
    push_field(buf, name, &v.to_le_bytes());
}
fn push_field_u64(buf: &mut Vec<u8>, name: &str, v: u64) {
    push_field(buf, name, &v.to_le_bytes());
}
fn push_field_usize(buf: &mut Vec<u8>, name: &str, v: usize) {
    push_field(buf, name, &(v as u64).to_le_bytes());
}
fn push_field_vec3d(buf: &mut Vec<u8>, name: &str, v: reb_vec3d) {
    let mut d = Vec::with_capacity(24);
    d.extend_from_slice(&v.x.to_le_bytes());
    d.extend_from_slice(&v.y.to_le_bytes());
    d.extend_from_slice(&v.z.to_le_bytes());
    push_field(buf, name, &d);
}
fn push_field_f64_slice(buf: &mut Vec<u8>, name: &str, v: &[f64]) {
    if v.is_empty() {
        return; // C: size_data == 0 -> field not written (REB_POINTER)
    }
    let mut d = Vec::with_capacity(8 * v.len());
    for x in v {
        d.extend_from_slice(&x.to_le_bytes());
    }
    push_field(buf, name, &d);
}

/// The synthetic "pointer" stored for a particle name (index into
/// `name_list`). Any non-zero value works; the C protocol only compares
/// these for equality with the values stored in the name_list blob.
fn name_fake_pointer(name: Option<usize>) -> u64 {
    match name {
        Some(i) => (i as u64) + 1,
        None => 0,
    }
}

/// Serialize one particle in the C x86-64 memory layout.
fn particle_to_raw(p: &reb_particle, out: &mut Vec<u8>) {
    out.extend_from_slice(&p.x.to_le_bytes());
    out.extend_from_slice(&p.y.to_le_bytes());
    out.extend_from_slice(&p.z.to_le_bytes());
    out.extend_from_slice(&p.vx.to_le_bytes());
    out.extend_from_slice(&p.vy.to_le_bytes());
    out.extend_from_slice(&p.vz.to_le_bytes());
    out.extend_from_slice(&p.ax.to_le_bytes());
    out.extend_from_slice(&p.ay.to_le_bytes());
    out.extend_from_slice(&p.az.to_le_bytes());
    out.extend_from_slice(&p.m.to_le_bytes());
    out.extend_from_slice(&p.r.to_le_bytes());
    out.extend_from_slice(&name_fake_pointer(p.name).to_le_bytes()); // name (char*)
    out.extend_from_slice(&0u64.to_le_bytes()); // ap (void*)
    out.extend_from_slice(&0u64.to_le_bytes()); // sim (struct reb_simulation*)
}

fn read_f64_le(d: &[u8], off: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[off..off + 8]);
    f64::from_le_bytes(b)
}
fn read_u64_le(d: &[u8], off: usize) -> u64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[off..off + 8]);
    u64::from_le_bytes(b)
}
fn read_i64_le(d: &[u8], off: usize) -> i64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[off..off + 8]);
    i64::from_le_bytes(b)
}
fn read_i32_le(d: &[u8], off: usize) -> i32 {
    let mut b = [0u8; 4];
    b.copy_from_slice(&d[off..off + 4]);
    i32::from_le_bytes(b)
}
fn read_u32_le(d: &[u8], off: usize) -> u32 {
    let mut b = [0u8; 4];
    b.copy_from_slice(&d[off..off + 4]);
    u32::from_le_bytes(b)
}

/// Deserialize one particle from the C memory layout; returns the
/// particle and the stored `name` pointer value (matched against the
/// name_list blob during the finish step).
fn particle_from_raw(d: &[u8]) -> (reb_particle, u64) {
    let mut p = reb_particle::default();
    p.x = read_f64_le(d, 0);
    p.y = read_f64_le(d, 8);
    p.z = read_f64_le(d, 16);
    p.vx = read_f64_le(d, 24);
    p.vy = read_f64_le(d, 32);
    p.vz = read_f64_le(d, 40);
    p.ax = read_f64_le(d, 48);
    p.ay = read_f64_le(d, 56);
    p.az = read_f64_le(d, 64);
    p.m = read_f64_le(d, 72);
    p.r = read_f64_le(d, 80);
    let name_ptr = read_u64_le(d, 88);
    (p, name_ptr)
}

fn f64_slice_from_raw(d: &[u8]) -> Vec<f64> {
    let n = d.len() / 8;
    let mut v = Vec::with_capacity(n);
    for i in 0..n {
        v.push(read_f64_le(d, 8 * i));
    }
    v
}

fn particles_from_raw(d: &[u8]) -> (Vec<reb_particle>, Vec<u64>) {
    let n = d.len() / REB_PARTICLE_RAW_SIZE;
    let mut ps = Vec::with_capacity(n);
    let mut ptrs = Vec::with_capacity(n);
    for i in 0..n {
        let (p, ptr) = particle_from_raw(&d[i * REB_PARTICLE_RAW_SIZE..]);
        ps.push(p);
        ptrs.push(ptr);
    }
    (ps, ptrs)
}

fn particles_to_raw(ps: &[reb_particle]) -> Vec<u8> {
    let mut out = Vec::with_capacity(ps.len() * REB_PARTICLE_RAW_SIZE);
    for p in ps {
        particle_to_raw(p, &mut out);
    }
    out
}

// ---------------------------------------------------------------------
// Serializer (binarydata.c `output_fields_from_list` +
// `reb_binarydata_simulation_to_stream`, with the field descriptor
// tables unrolled)
// ---------------------------------------------------------------------

/// Emit the integrator state fields (C: the per-integrator
/// `field_descriptor_list`s), with names prefixed
/// "integrator.<name>.".
fn output_integrator_fields(buf: &mut Vec<u8>, r: &reb_simulation) {
    let p = format!("integrator.{}", r.integrator.name());
    let n = |field: &str| format!("{}.{}", p, field);
    match &r.integrator {
        reb_integrator_state::none => {}
        reb_integrator_state::sei(s) => {
            push_field_f64(buf, &n("lastdt"), s.lastdt);
            push_field_f64(buf, &n("sindt"), s.sindt);
            push_field_f64(buf, &n("tandt"), s.tandt);
            push_field_f64(buf, &n("sindtz"), s.sindtz);
            push_field_f64(buf, &n("tandtz"), s.tandtz);
        }
        reb_integrator_state::leapfrog(s) => {
            push_field_u32(buf, &n("order"), s.order);
        }
        reb_integrator_state::ias15(s) => {
            push_field_f64(buf, &n("epsilon"), s.epsilon);
            push_field_f64(buf, &n("min_dt"), s.min_dt);
            push_field_u32(buf, &n("adaptive_mode"), s.adaptive_mode);
            push_field_u64(buf, &n("iterations_max_exceeded"), s.iterations_max_exceeded);
            push_field_f64_slice(buf, &n("at"), &s.at);
            push_field_f64_slice(buf, &n("x0"), &s.x0);
            push_field_f64_slice(buf, &n("v0"), &s.v0);
            push_field_f64_slice(buf, &n("a0"), &s.a0);
            push_field_f64_slice(buf, &n("csx"), &s.csx);
            push_field_f64_slice(buf, &n("csv"), &s.csv);
            push_field_f64_slice(buf, &n("csa0"), &s.csa0);
            push_field_f64_slice(buf, &n("g"), &s.g);
            push_field_f64_slice(buf, &n("b"), &s.b);
            push_field_f64_slice(buf, &n("csb"), &s.csb);
            push_field_f64_slice(buf, &n("e"), &s.e);
            push_field_f64_slice(buf, &n("br"), &s.br);
            push_field_f64_slice(buf, &n("er"), &s.er);
        }
        reb_integrator_state::whfast(s) => {
            push_field_u32(buf, &n("corrector"), s.corrector);
            push_field_u32(buf, &n("safe_mode"), s.safe_mode);
            push_field_i32(buf, &n("coordinates"), s.coordinates as i32);
            push_field_u32(buf, &n("corrector2"), s.corrector2);
            push_field_i32(buf, &n("kernel"), s.kernel as i32);
            push_field_u32(buf, &n("keep_unsynchronized"), s.keep_unsynchronized);
            if !s.p_jh.is_empty() {
                push_field(buf, &n("p_jh"), &particles_to_raw(&s.p_jh));
            }
            if !s.p_jh_var.is_empty() {
                push_field(buf, &n("p_jh_var"), &particles_to_raw(&s.p_jh_var));
            }
        }
        reb_integrator_state::saba(s) => {
            push_field_u32(buf, &n("safe_mode"), s.safe_mode);
            push_field_i32(buf, &n("type"), s.type_ as i32);
            push_field_u32(buf, &n("keep_unsynchronized"), s.keep_unsynchronized);
            if !s.p_jh.is_empty() {
                push_field(buf, &n("p_jh"), &particles_to_raw(&s.p_jh));
            }
        }
        reb_integrator_state::janus(s) => {
            push_field_f64(buf, &n("scale_pos"), s.scale_pos);
            push_field_f64(buf, &n("scale_vel"), s.scale_vel);
            push_field_u32(buf, &n("order"), s.order);
            push_field_u32(
                buf,
                &n("recalculate_integer_coordinates_this_timestep"),
                s.recalculate_integer_coordinates_this_timestep,
            );
            if !s.p_int.is_empty() {
                let mut d = Vec::with_capacity(s.p_int.len() * REB_PARTICLE_INT_RAW_SIZE);
                for pi in &s.p_int {
                    d.extend_from_slice(&pi.x.to_le_bytes());
                    d.extend_from_slice(&pi.y.to_le_bytes());
                    d.extend_from_slice(&pi.z.to_le_bytes());
                    d.extend_from_slice(&pi.vx.to_le_bytes());
                    d.extend_from_slice(&pi.vy.to_le_bytes());
                    d.extend_from_slice(&pi.vz.to_le_bytes());
                }
                push_field(buf, &n("p_int"), &d);
            }
        }
        reb_integrator_state::eos(s) => {
            push_field_i32(buf, &n("phi0"), s.phi0);
            push_field_i32(buf, &n("phi1"), s.phi1);
            push_field_u32(buf, &n("n"), s.n);
            push_field_u32(buf, &n("safe_mode"), s.safe_mode);
        }
        reb_integrator_state::mercurius(s) => {
            // (The `L` function pointer is REB_FUNCTIONPOINTER: not written.)
            push_field_f64(buf, &n("r_crit_hill"), s.r_crit_hill);
            push_field_u32(buf, &n("safe_mode"), s.safe_mode);
            push_field_f64_slice(buf, &n("dcrit"), &s.dcrit);
            push_field_vec3d(buf, &n("com_pos"), s.com_pos);
            push_field_vec3d(buf, &n("com_vel"), s.com_vel);
        }
        reb_integrator_state::bs(s) => {
            push_field_f64(buf, &n("eps_abs"), s.eps_abs);
            push_field_f64(buf, &n("eps_rel"), s.eps_rel);
            push_field_f64(buf, &n("min_dt"), s.min_dt);
            push_field_f64(buf, &n("max_dt"), s.max_dt);
            push_field_i32(buf, &n("first_or_last_step"), s.first_or_last_step);
            push_field_i32(buf, &n("previous_rejected"), s.previous_rejected);
            push_field_i32(buf, &n("target_iter"), s.target_iter);
        }
        reb_integrator_state::trace(s) => {
            // (S and S_peri are REB_FUNCTIONPOINTER: not written.)
            push_field_f64(buf, &n("r_crit_hill"), s.r_crit_hill);
            push_field_f64(buf, &n("peri_crit_eta"), s.peri_crit_eta);
            push_field_i32(buf, &n("peri_mode"), s.peri_mode);
        }
        reb_integrator_state::whfast512(s) => {
            push_field_u32(buf, &n("gr_potential"), s.gr_potential);
            push_field_u32(buf, &n("corrector"), s.corrector);
            push_field_u64(buf, &n("concatenate_steps"), s.concatenate_steps);
            push_field_u32(buf, &n("N_systems"), s.N_systems);
            // ("data" is REB_POINTER_ALIGNED with a NULL pointer on the
            // Windows reference build: size 0, not written.)
            push_field_f64(buf, &n("last_synchronization"), s.last_synchronization);
        }
    }
}

/// binarydata.c `reb_binarydata_simulation_to_stream` — serializes a
/// simulation to a buffer (the C uses out-parameters `bufp`/`sizep`;
/// Rust returns the buffer).
pub fn reb_binarydata_simulation_to_stream(r: &mut reb_simulation) -> Vec<u8> {
    if r.simulationarchive_version < 5 {
        reb_simulation_error(
            r,
            "Simulationarchives with version < 5 are no longer supported.\n",
        );
    }
    let mut buf: Vec<u8> = Vec::new();

    // Output header.
    let mut header = [0u8; 64];
    let s = format!("REBOUND Binary File. Version: {}", reb_version_str);
    header[..s.len()].copy_from_slice(s.as_bytes());
    // C: snprintf(header+cwritten+1, ...) — the githash follows the NUL.
    let g = reb_githash_str.as_bytes();
    let gstart = s.len() + 1;
    let gcap = 64usize.saturating_sub(gstart + 1); // snprintf keeps a trailing NUL
    let glen = std::cmp::min(g.len(), gcap);
    header[gstart..gstart + glen].copy_from_slice(&g[..glen]);
    push_bytes(&mut buf, &header);

    // Output all fields — same order as reb_binarydata_field_descriptor_list.
    push_field_f64(&mut buf, "t", r.t);
    push_field_f64(&mut buf, "G", r.G);
    push_field_f64(&mut buf, "softening", r.softening);
    push_field_f64(&mut buf, "dt", r.dt);
    push_field_usize(&mut buf, "N", r.N);
    push_field_usize(&mut buf, "N_var", r.N_var);
    push_field_usize(&mut buf, "N_active", r.N_active);
    push_field_i32(&mut buf, "testparticle_type", r.testparticle_type);
    push_field_f64(&mut buf, "opening_angle2", r.opening_angle2);
    push_field_i32(&mut buf, "status", r.status);
    push_field_i32(&mut buf, "exact_finish_time", r.exact_finish_time);
    push_field_u32(
        &mut buf,
        "force_is_velocity_dependent",
        r.force_is_velocity_dependent as u32,
    );
    push_field_u32(&mut buf, "gravity_ignore_terms", r.gravity_ignore_terms);
    push_field_f64(&mut buf, "output_timing_last", r.output_timing_last);
    push_field_i32(&mut buf, "save_messages", r.save_messages);
    push_field_f64(&mut buf, "exit_max_distance", r.exit_max_distance);
    push_field_f64(&mut buf, "exit_min_distance", r.exit_min_distance);
    push_field_f64(&mut buf, "usleep", r.usleep);
    push_field_i32(&mut buf, "track_energy_offset", r.track_energy_offset);
    push_field_f64(&mut buf, "energy_offset", r.energy_offset);
    push_field_f64(&mut buf, "root_size", r.root_size);
    push_field_usize(&mut buf, "N_root_x", r.N_root_x);
    push_field_usize(&mut buf, "N_root_y", r.N_root_y);
    push_field_usize(&mut buf, "N_root_z", r.N_root_z);
    push_field_i32(&mut buf, "N_ghost_x", r.N_ghost_x);
    push_field_i32(&mut buf, "N_ghost_y", r.N_ghost_y);
    push_field_i32(&mut buf, "N_ghost_z", r.N_ghost_z);
    push_field_f64(
        &mut buf,
        "minimum_collision_velocity",
        r.minimum_collision_velocity,
    );
    push_field_f64(&mut buf, "collisions_plog", r.collisions_plog);
    push_field_i64(&mut buf, "collisions_log_n", r.collisions_log_n);
    push_field_i32(&mut buf, "calculate_megno", r.calculate_megno);
    push_field_f64(&mut buf, "megno_Ys", r.megno_Ys);
    push_field_f64(&mut buf, "megno_Yss", r.megno_Yss);
    push_field_f64(&mut buf, "megno_cov_Yt", r.megno_cov_Yt);
    push_field_f64(&mut buf, "megno_var_t", r.megno_var_t);
    push_field_f64(&mut buf, "megno_mean_t", r.megno_mean_t);
    push_field_f64(&mut buf, "megno_mean_Y", r.megno_mean_Y);
    push_field_f64(&mut buf, "megno_initial_t", r.megno_initial_t);
    push_field_i64(&mut buf, "megno_n", r.megno_n);
    push_field_f64(
        &mut buf,
        "simulationarchive_auto_interval",
        r.simulationarchive_auto_interval,
    );
    push_field_f64(
        &mut buf,
        "simulationarchive_auto_walltime",
        r.simulationarchive_auto_walltime,
    );
    push_field_f64(&mut buf, "simulationarchive_next", r.simulationarchive_next);
    push_field_i32(&mut buf, "collision", r.collision as i32);
    {
        // "integrator.name" (REB_STRING): name bytes + NUL.
        let name = r.integrator.name();
        let mut d = Vec::with_capacity(name.len() + 1);
        d.extend_from_slice(name.as_bytes());
        d.push(0);
        push_field(&mut buf, "integrator.name", &d);
    }
    push_field_i32(&mut buf, "boundary", r.boundary as i32);
    push_field_i32(&mut buf, "gravity", r.gravity as i32);
    push_field_f64(&mut buf, "OMEGA", r.OMEGA);
    push_field_f64(&mut buf, "OMEGAZ", r.OMEGAZ);
    push_field_u32(&mut buf, "is_synchronized", r.is_synchronized);
    push_field_u32(&mut buf, "did_modify_particles", r.did_modify_particles);
    if r.N != 0 {
        push_field(&mut buf, "particles", &particles_to_raw(&r.particles[..r.N]));
    }
    if r.N_var != 0 {
        push_field(
            &mut buf,
            "particles_var",
            &particles_to_raw(&r.particles_var[..r.N_var]),
        );
    }
    if !r.var_config.is_empty() {
        // REB_POINTER with element size sizeof(struct reb_variational_configuration).
        let mut d = Vec::with_capacity(r.var_config.len() * REB_VAR_CONFIG_RAW_SIZE);
        for vc in &r.var_config {
            d.extend_from_slice(&0u64.to_le_bytes()); // sim (pointer)
            d.extend_from_slice(&vc.order.to_le_bytes());
            d.extend_from_slice(&(vc.index as i32).to_le_bytes());
            d.extend_from_slice(&vc.testparticle.to_le_bytes());
            d.extend_from_slice(&(vc.index_1st_order_a as i32).to_le_bytes());
            d.extend_from_slice(&(vc.index_1st_order_b as i32).to_le_bytes());
            d.extend_from_slice(&0u32.to_le_bytes()); // struct padding
            d.extend_from_slice(&vc.lrescale.to_le_bytes());
        }
        push_field(&mut buf, "var_config", &d);
    }
    push_field_i32(&mut buf, "simulationarchive_version", r.simulationarchive_version);
    push_field_f64(&mut buf, "walltime", r.walltime);
    push_field_f64(&mut buf, "walltime_last_steps", r.walltime_last_steps);
    push_field_u32(&mut buf, "python_unit_l", r.python_unit_l);
    push_field_u32(&mut buf, "python_unit_m", r.python_unit_m);
    push_field_u32(&mut buf, "python_unit_t", r.python_unit_t);
    push_field_u64(&mut buf, "simulationarchive_auto_step", r.simulationarchive_auto_step);
    push_field_u64(&mut buf, "simulationarchive_next_step", r.simulationarchive_next_step);
    push_field_u64(&mut buf, "steps_done", r.steps_done);
    push_field_f64(&mut buf, "dt_last_done", r.dt_last_done);
    push_field_u32(&mut buf, "rand_seed", r.rand_seed);
    push_field_i32(&mut buf, "testparticle_hidewarnings", r.testparticle_hidewarnings);
    // "display_settings" (REB_POINTER, SIZE_MAX): the display subsystem
    // is excluded from this port; the pointer is always NULL -> skipped.
    if !r.name_list.is_empty() {
        // REB_CHARP_LIST: each entry is the string incl. NUL, followed
        // by the (synthetic) pointer stored in the particles.
        let mut d = Vec::new();
        for (i, s) in r.name_list.iter().enumerate() {
            d.extend_from_slice(s.as_bytes());
            d.push(0);
            d.extend_from_slice(&name_fake_pointer(Some(i)).to_le_bytes());
        }
        push_field(&mut buf, "name_list", &d);
    }

    // Integrator state fields.
    output_integrator_fields(&mut buf, r);

    // Write function pointer warning flag.
    let functionpointersused: i32 = if r.coefficient_of_restitution.is_some()
        || r.collision_resolve.is_some()
        || r.additional_forces.is_some()
        || r.heartbeat.is_some()
        || r.post_timestep_modifications.is_some()
    {
        1
    } else {
        0
    };
    push_field_i32(&mut buf, "functionpointers", functionpointersused);

    // Write last field.
    push_field(&mut buf, "end", &[]);

    // Trailing zeroed reb_simulationarchive_blob.
    push_bytes(&mut buf, &[0u8; REB_SA_BLOB_SIZE]);
    buf
}

// ---------------------------------------------------------------------
// Reader (binarydata.c `reb_binarydata_input_fields`)
// ---------------------------------------------------------------------

fn read_exact_or_eof<R: Read>(inf: &mut R, out: &mut [u8]) -> bool {
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

fn read_field_header<R: Read>(inf: &mut R) -> Option<reb_binarydata_field> {
    let mut b = [0u8; REB_BINARYDATA_FIELD_SIZE];
    if !read_exact_or_eof(inf, &mut b) {
        return None;
    }
    Some(reb_binarydata_field {
        size_name: u64::from_le_bytes(b[0..8].try_into().unwrap_or([0; 8])),
        size_data: u64::from_le_bytes(b[8..16].try_into().unwrap_or([0; 8])),
    })
}

/// Apply one "integrator.<int>.<field>" data blob to the current
/// integrator state. Returns false if the field is unknown.
fn input_integrator_field(r: &mut reb_simulation, field: &str, d: &[u8]) -> bool {
    let mut state = std::mem::replace(&mut r.integrator, reb_integrator_state::none);
    let known = match &mut state {
        reb_integrator_state::none => false,
        reb_integrator_state::sei(s) => match field {
            "lastdt" => {
                s.lastdt = read_f64_le(d, 0);
                true
            }
            "sindt" => {
                s.sindt = read_f64_le(d, 0);
                true
            }
            "tandt" => {
                s.tandt = read_f64_le(d, 0);
                true
            }
            "sindtz" => {
                s.sindtz = read_f64_le(d, 0);
                true
            }
            "tandtz" => {
                s.tandtz = read_f64_le(d, 0);
                true
            }
            _ => false,
        },
        reb_integrator_state::leapfrog(s) => match field {
            "order" => {
                s.order = read_u32_le(d, 0);
                true
            }
            _ => false,
        },
        reb_integrator_state::ias15(s) => match field {
            "epsilon" => {
                s.epsilon = read_f64_le(d, 0);
                true
            }
            "min_dt" => {
                s.min_dt = read_f64_le(d, 0);
                true
            }
            "adaptive_mode" => {
                s.adaptive_mode = read_u32_le(d, 0);
                true
            }
            "iterations_max_exceeded" => {
                s.iterations_max_exceeded = read_u64_le(d, 0);
                true
            }
            "at" => {
                s.at = f64_slice_from_raw(d);
                s.N_allocated = s.at.len();
                true
            }
            "x0" => {
                s.x0 = f64_slice_from_raw(d);
                s.N_allocated = s.x0.len();
                true
            }
            "v0" => {
                s.v0 = f64_slice_from_raw(d);
                s.N_allocated = s.v0.len();
                true
            }
            "a0" => {
                s.a0 = f64_slice_from_raw(d);
                s.N_allocated = s.a0.len();
                true
            }
            "csx" => {
                s.csx = f64_slice_from_raw(d);
                s.N_allocated = s.csx.len();
                true
            }
            "csv" => {
                s.csv = f64_slice_from_raw(d);
                s.N_allocated = s.csv.len();
                true
            }
            "csa0" => {
                s.csa0 = f64_slice_from_raw(d);
                s.N_allocated = s.csa0.len();
                true
            }
            "g" => {
                s.g = f64_slice_from_raw(d);
                s.N_allocated = s.g.len() / 7;
                true
            }
            "b" => {
                s.b = f64_slice_from_raw(d);
                s.N_allocated = s.b.len() / 7;
                true
            }
            "csb" => {
                s.csb = f64_slice_from_raw(d);
                s.N_allocated = s.csb.len() / 7;
                true
            }
            "e" => {
                s.e = f64_slice_from_raw(d);
                s.N_allocated = s.e.len() / 7;
                true
            }
            "br" => {
                s.br = f64_slice_from_raw(d);
                s.N_allocated = s.br.len() / 7;
                true
            }
            "er" => {
                s.er = f64_slice_from_raw(d);
                s.N_allocated = s.er.len() / 7;
                true
            }
            _ => false,
        },
        reb_integrator_state::whfast(s) => match field {
            "corrector" => {
                s.corrector = read_u32_le(d, 0);
                true
            }
            "safe_mode" => {
                s.safe_mode = read_u32_le(d, 0);
                true
            }
            "coordinates" => {
                s.coordinates = read_i32_le(d, 0) as u32;
                true
            }
            "corrector2" => {
                s.corrector2 = read_u32_le(d, 0);
                true
            }
            "kernel" => {
                s.kernel = read_i32_le(d, 0) as u32;
                true
            }
            "keep_unsynchronized" => {
                s.keep_unsynchronized = read_u32_le(d, 0);
                true
            }
            "p_jh" => {
                s.p_jh = particles_from_raw(d).0;
                true
            }
            "p_jh_var" => {
                s.p_jh_var = particles_from_raw(d).0;
                true
            }
            _ => false,
        },
        reb_integrator_state::saba(s) => match field {
            "safe_mode" => {
                s.safe_mode = read_u32_le(d, 0);
                true
            }
            "type" => {
                s.type_ = read_i32_le(d, 0);
                true
            }
            "keep_unsynchronized" => {
                s.keep_unsynchronized = read_u32_le(d, 0);
                true
            }
            "p_jh" => {
                s.p_jh = particles_from_raw(d).0;
                true
            }
            _ => false,
        },
        reb_integrator_state::janus(s) => match field {
            "scale_pos" => {
                s.scale_pos = read_f64_le(d, 0);
                true
            }
            "scale_vel" => {
                s.scale_vel = read_f64_le(d, 0);
                true
            }
            "order" => {
                s.order = read_u32_le(d, 0);
                true
            }
            "recalculate_integer_coordinates_this_timestep" => {
                s.recalculate_integer_coordinates_this_timestep = read_u32_le(d, 0);
                true
            }
            "p_int" => {
                let n = d.len() / REB_PARTICLE_INT_RAW_SIZE;
                let mut v = Vec::with_capacity(n);
                for i in 0..n {
                    let o = i * REB_PARTICLE_INT_RAW_SIZE;
                    v.push(crate::integrator_janus::reb_particle_int {
                        x: read_i64_le(d, o),
                        y: read_i64_le(d, o + 8),
                        z: read_i64_le(d, o + 16),
                        vx: read_i64_le(d, o + 24),
                        vy: read_i64_le(d, o + 32),
                        vz: read_i64_le(d, o + 40),
                    });
                }
                s.p_int = v;
                true
            }
            _ => false,
        },
        reb_integrator_state::eos(s) => match field {
            "phi0" => {
                s.phi0 = read_i32_le(d, 0);
                true
            }
            "phi1" => {
                s.phi1 = read_i32_le(d, 0);
                true
            }
            "n" => {
                s.n = read_u32_le(d, 0);
                true
            }
            "safe_mode" => {
                s.safe_mode = read_u32_le(d, 0);
                true
            }
            _ => false,
        },
        reb_integrator_state::mercurius(s) => match field {
            "r_crit_hill" => {
                s.r_crit_hill = read_f64_le(d, 0);
                true
            }
            "safe_mode" => {
                s.safe_mode = read_u32_le(d, 0);
                true
            }
            "dcrit" => {
                s.dcrit = f64_slice_from_raw(d);
                true
            }
            "com_pos" => {
                s.com_pos = reb_vec3d {
                    x: read_f64_le(d, 0),
                    y: read_f64_le(d, 8),
                    z: read_f64_le(d, 16),
                };
                true
            }
            "com_vel" => {
                s.com_vel = reb_vec3d {
                    x: read_f64_le(d, 0),
                    y: read_f64_le(d, 8),
                    z: read_f64_le(d, 16),
                };
                true
            }
            _ => false,
        },
        reb_integrator_state::bs(s) => match field {
            "eps_abs" => {
                s.eps_abs = read_f64_le(d, 0);
                true
            }
            "eps_rel" => {
                s.eps_rel = read_f64_le(d, 0);
                true
            }
            "min_dt" => {
                s.min_dt = read_f64_le(d, 0);
                true
            }
            "max_dt" => {
                s.max_dt = read_f64_le(d, 0);
                true
            }
            "first_or_last_step" => {
                s.first_or_last_step = read_i32_le(d, 0);
                true
            }
            "previous_rejected" => {
                s.previous_rejected = read_i32_le(d, 0);
                true
            }
            "target_iter" => {
                s.target_iter = read_i32_le(d, 0);
                true
            }
            _ => false,
        },
        reb_integrator_state::trace(s) => match field {
            "r_crit_hill" => {
                s.r_crit_hill = read_f64_le(d, 0);
                true
            }
            "peri_crit_eta" => {
                s.peri_crit_eta = read_f64_le(d, 0);
                true
            }
            "peri_mode" => {
                s.peri_mode = read_i32_le(d, 0);
                true
            }
            _ => false,
        },
        reb_integrator_state::whfast512(s) => match field {
            "gr_potential" => {
                s.gr_potential = read_u32_le(d, 0);
                true
            }
            "corrector" => {
                s.corrector = read_u32_le(d, 0);
                true
            }
            "concatenate_steps" => {
                s.concatenate_steps = read_u64_le(d, 0);
                true
            }
            "N_systems" => {
                s.N_systems = read_u32_le(d, 0);
                true
            }
            "data" => true, // SIMD block: not carried on the Windows reference build
            "last_synchronization" => {
                s.last_synchronization = read_f64_le(d, 0);
                true
            }
            _ => false,
        },
    };
    r.integrator = state;
    known
}

fn collision_from_i32(v: i32) -> REB_COLLISION {
    match v {
        1 => REB_COLLISION::DIRECT,
        2 => REB_COLLISION::TREE,
        4 => REB_COLLISION::LINE,
        5 => REB_COLLISION::LINETREE,
        _ => REB_COLLISION::NONE,
    }
}
fn boundary_from_i32(v: i32) -> REB_BOUNDARY {
    match v {
        1 => REB_BOUNDARY::OPEN,
        2 => REB_BOUNDARY::PERIODIC,
        3 => REB_BOUNDARY::SHEAR,
        _ => REB_BOUNDARY::NONE,
    }
}
fn gravity_from_i32(v: i32) -> REB_GRAVITY {
    match v {
        1 => REB_GRAVITY::BASIC,
        2 => REB_GRAVITY::COMPENSATED,
        3 => REB_GRAVITY::TREE,
        5 => REB_GRAVITY::JACOBI,
        7 => REB_GRAVITY::CUSTOM,
        _ => REB_GRAVITY::NONE,
    }
}

/// binarydata.c `reb_binarydata_input_fields` — reads field data into
/// the simulation from a stream (File or `std::io::Cursor` over a
/// buffer; the latter replaces the C's non-portable fmemopen).
pub fn reb_binarydata_input_fields<R: Read + Seek>(
    r: &mut reb_simulation,
    inf: &mut R,
    warnings: &mut REB_BINARYDATA_ERROR_CODE,
) {
    // The stored particle name "pointers" and the name_list blob's
    // stored pointers, matched in the finish step (like the C).
    let mut particle_name_ptrs: Vec<u64> = Vec::new();
    let mut particle_var_name_ptrs: Vec<u64> = Vec::new();
    let mut name_list_ptrs: Vec<u64> = Vec::new();

    'fields: loop {
        let field = match read_field_header(inf) {
            Some(f) => f,
            None => break 'fields, // End of file
        };
        // Is this a real field or the header?
        if field.size_name == reb_binarydata_header {
            let bufsize = 64 - REB_BINARYDATA_FIELD_SIZE;
            let mut readbuf = vec![0u8; bufsize];
            let header = "REBOUND Binary File. Version: ";
            let curv = format!("{}{}", &header[REB_BINARYDATA_FIELD_SIZE..], reb_version_str);
            if !read_exact_or_eof(inf, &mut readbuf) {
                *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            } else {
                // Note: compares version, ignores githash.
                let mut matches = true;
                for (i, b) in curv.as_bytes().iter().enumerate() {
                    if i >= bufsize || readbuf[i] != *b {
                        matches = false;
                        break;
                    }
                }
                // strncmp also requires the NUL right after the version
                if matches && curv.len() < bufsize && readbuf[curv.len()] != 0 {
                    matches = false;
                }
                if !matches {
                    *warnings |= REB_BINARYDATA_WARNING_VERSION;
                }
            }
            continue 'fields;
        }
        // Try to get name of field
        if field.size_name as usize > REB_STRING_SIZE_MAX {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break 'fields;
        }
        let mut name_bytes = vec![0u8; field.size_name as usize];
        if !read_exact_or_eof(inf, &mut name_bytes) {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break 'fields;
        }
        let nul = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
        let name = String::from_utf8_lossy(&name_bytes[..nul]).to_string();

        // Fields that require special handling
        if name == "end" {
            break 'fields; // End of snapshot
        }
        let mut data = vec![0u8; field.size_data as usize];
        if name == "functionpointers" {
            // Warning for when function pointers were used.
            // No effect on simulation.
            if read_exact_or_eof(inf, &mut data) && !data.is_empty() && read_i32_le(&data, 0) != 0 {
                *warnings |= REB_BINARYDATA_WARNING_POINTERS;
            }
            continue 'fields;
        }
        if name == "integrator.name" {
            if !read_exact_or_eof(inf, &mut data) {
                *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
                break 'fields;
            }
            let nul = data.iter().position(|&b| b == 0).unwrap_or(data.len());
            let iname = String::from_utf8_lossy(&data[..nul]).to_string();
            reb_simulation_set_integrator(r, &iname);
            continue 'fields;
        }
        // Check that we only read integrator values that match the
        // currently set integrator.
        if let Some(rest) = name.strip_prefix("integrator.") {
            if let Some(dot) = rest.find('.') {
                let (iname, ifield) = rest.split_at(dot);
                if r.integrator.name() != iname {
                    *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
                    break 'fields;
                }
                if !read_exact_or_eof(inf, &mut data) {
                    *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
                    break 'fields;
                }
                if !input_integrator_field(r, &ifield[1..], &data) {
                    *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
                    break 'fields;
                }
                continue 'fields;
            }
        }

        if !read_exact_or_eof(inf, &mut data) {
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            break 'fields;
        }
        let d = &data[..];
        match name.as_str() {
            "t" => r.t = read_f64_le(d, 0),
            "G" => r.G = read_f64_le(d, 0),
            "softening" => r.softening = read_f64_le(d, 0),
            "dt" => r.dt = read_f64_le(d, 0),
            "N" => r.N = read_u64_le(d, 0) as usize,
            "N_var" => r.N_var = read_u64_le(d, 0) as usize,
            "N_active" => r.N_active = read_u64_le(d, 0) as usize,
            "testparticle_type" => r.testparticle_type = read_i32_le(d, 0),
            "opening_angle2" => r.opening_angle2 = read_f64_le(d, 0),
            "status" => r.status = read_i32_le(d, 0),
            "exact_finish_time" => r.exact_finish_time = read_i32_le(d, 0),
            "force_is_velocity_dependent" => {
                r.force_is_velocity_dependent = read_u32_le(d, 0) as i32
            }
            "gravity_ignore_terms" => r.gravity_ignore_terms = read_u32_le(d, 0),
            "output_timing_last" => r.output_timing_last = read_f64_le(d, 0),
            "save_messages" => r.save_messages = read_i32_le(d, 0),
            "exit_max_distance" => r.exit_max_distance = read_f64_le(d, 0),
            "exit_min_distance" => r.exit_min_distance = read_f64_le(d, 0),
            "usleep" => r.usleep = read_f64_le(d, 0),
            "track_energy_offset" => r.track_energy_offset = read_i32_le(d, 0),
            "energy_offset" => r.energy_offset = read_f64_le(d, 0),
            "root_size" => r.root_size = read_f64_le(d, 0),
            "N_root_x" => r.N_root_x = read_u64_le(d, 0) as usize,
            "N_root_y" => r.N_root_y = read_u64_le(d, 0) as usize,
            "N_root_z" => r.N_root_z = read_u64_le(d, 0) as usize,
            "N_ghost_x" => r.N_ghost_x = read_i32_le(d, 0),
            "N_ghost_y" => r.N_ghost_y = read_i32_le(d, 0),
            "N_ghost_z" => r.N_ghost_z = read_i32_le(d, 0),
            "minimum_collision_velocity" => r.minimum_collision_velocity = read_f64_le(d, 0),
            "collisions_plog" => r.collisions_plog = read_f64_le(d, 0),
            "collisions_log_n" => r.collisions_log_n = read_i64_le(d, 0),
            "calculate_megno" => r.calculate_megno = read_i32_le(d, 0),
            "megno_Ys" => r.megno_Ys = read_f64_le(d, 0),
            "megno_Yss" => r.megno_Yss = read_f64_le(d, 0),
            "megno_cov_Yt" => r.megno_cov_Yt = read_f64_le(d, 0),
            "megno_var_t" => r.megno_var_t = read_f64_le(d, 0),
            "megno_mean_t" => r.megno_mean_t = read_f64_le(d, 0),
            "megno_mean_Y" => r.megno_mean_Y = read_f64_le(d, 0),
            "megno_initial_t" => r.megno_initial_t = read_f64_le(d, 0),
            "megno_n" => r.megno_n = read_i64_le(d, 0),
            "simulationarchive_auto_interval" => {
                r.simulationarchive_auto_interval = read_f64_le(d, 0)
            }
            "simulationarchive_auto_walltime" => {
                r.simulationarchive_auto_walltime = read_f64_le(d, 0)
            }
            "simulationarchive_next" => r.simulationarchive_next = read_f64_le(d, 0),
            "collision" => r.collision = collision_from_i32(read_i32_le(d, 0)),
            "boundary" => r.boundary = boundary_from_i32(read_i32_le(d, 0)),
            "gravity" => r.gravity = gravity_from_i32(read_i32_le(d, 0)),
            "OMEGA" => r.OMEGA = read_f64_le(d, 0),
            "OMEGAZ" => r.OMEGAZ = read_f64_le(d, 0),
            "is_synchronized" => r.is_synchronized = read_u32_le(d, 0),
            "did_modify_particles" => r.did_modify_particles = read_u32_le(d, 0),
            "particles" => {
                let (ps, ptrs) = particles_from_raw(d);
                r.particles = ps;
                particle_name_ptrs = ptrs;
            }
            "particles_var" => {
                let (ps, ptrs) = particles_from_raw(d);
                r.particles_var = ps;
                particle_var_name_ptrs = ptrs;
            }
            "var_config" => {
                let n = d.len() / REB_VAR_CONFIG_RAW_SIZE;
                let mut v = Vec::with_capacity(n);
                for i in 0..n {
                    let o = i * REB_VAR_CONFIG_RAW_SIZE;
                    v.push(reb_variational_configuration {
                        order: read_i32_le(d, o + 8),
                        index: read_i32_le(d, o + 12) as usize,
                        testparticle: read_i32_le(d, o + 16),
                        index_1st_order_a: read_i32_le(d, o + 20) as usize,
                        index_1st_order_b: read_i32_le(d, o + 24) as usize,
                        lrescale: read_f64_le(d, o + 32),
                    });
                }
                r.var_config = v;
            }
            "simulationarchive_version" => r.simulationarchive_version = read_i32_le(d, 0),
            "walltime" => r.walltime = read_f64_le(d, 0),
            "walltime_last_steps" => r.walltime_last_steps = read_f64_le(d, 0),
            "python_unit_l" => r.python_unit_l = read_u32_le(d, 0),
            "python_unit_m" => r.python_unit_m = read_u32_le(d, 0),
            "python_unit_t" => r.python_unit_t = read_u32_le(d, 0),
            "simulationarchive_auto_step" => r.simulationarchive_auto_step = read_u64_le(d, 0),
            "simulationarchive_next_step" => r.simulationarchive_next_step = read_u64_le(d, 0),
            "steps_done" => r.steps_done = read_u64_le(d, 0),
            "dt_last_done" => r.dt_last_done = read_f64_le(d, 0),
            "rand_seed" => r.rand_seed = read_u32_le(d, 0),
            "testparticle_hidewarnings" => r.testparticle_hidewarnings = read_i32_le(d, 0),
            "display_settings" => {
                // Display subsystem excluded from this port; ignore payload.
            }
            "name_list" => {
                // Serialized as: string bytes + NUL + original pointer.
                r.name_list.clear();
                name_list_ptrs.clear();
                let mut pos = 0usize;
                while pos < d.len() {
                    let nul = match d[pos..].iter().position(|&b| b == 0) {
                        Some(p) => p,
                        None => break,
                    };
                    let s = String::from_utf8_lossy(&d[pos..pos + nul]).to_string();
                    let ptr_off = pos + nul + 1;
                    if ptr_off + 8 > d.len() {
                        break;
                    }
                    let ptr = read_u64_le(d, ptr_off);
                    r.name_list.push(s);
                    name_list_ptrs.push(ptr);
                    pos = ptr_off + 8;
                }
            }
            _ => {
                // Unknown field id (C: REB_FIELD_NOT_FOUND).
                *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
                // (Data already consumed; the C seeks past it and stops.)
                let _ = inf.seek(SeekFrom::Current(0));
                break 'fields;
            }
        }
    }

    // Some final initializations (the C's pointer restorations).
    // Restore particle names from the stored pointer values.
    for l in 0..r.particles.len() {
        let ptr = particle_name_ptrs.get(l).copied().unwrap_or(0);
        if ptr != 0 {
            let mut name_found = None;
            for (n, &stored) in name_list_ptrs.iter().enumerate() {
                if ptr == stored {
                    name_found = Some(n);
                }
            }
            if name_found.is_none() {
                reb_simulation_warning(
                    r,
                    "A name for a particle was not stored in the Simulationarchive.",
                );
            }
            r.particles[l].name = name_found;
        } else {
            r.particles[l].name = None;
        }
    }
    for l in 0..r.particles_var.len() {
        let ptr = particle_var_name_ptrs.get(l).copied().unwrap_or(0);
        let _ = ptr; // variational particles never carry names
        r.particles_var[l].name = None;
    }
}

/// binarydata.c `reb_binarydata_process_warnings`. The C frees the
/// simulation and returns NULL on fatal errors; here the return value
/// is 0 on success and -1 when the simulation should be discarded.
pub fn reb_binarydata_process_warnings(
    r: &mut reb_simulation,
    warnings: REB_BINARYDATA_ERROR_CODE,
) -> i32 {
    if warnings & REB_BINARYDATA_ERROR_NOFILE != 0 {
        reb_simulation_error(r, "Cannot read binary file. Check filename and file contents.");
        return -1;
    }
    if warnings & REB_BINARYDATA_WARNING_VERSION != 0 {
        reb_simulation_warning(r, "Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if warnings & REB_BINARYDATA_WARNING_POINTERS != 0 {
        reb_simulation_warning(r, "You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if warnings & REB_BINARYDATA_WARNING_PARTICLES != 0 {
        reb_simulation_warning(r, "Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if warnings & REB_BINARYDATA_ERROR_FILENOTOPEN != 0 {
        reb_simulation_error(r, "Error while reading binary file (file was not open).");
        return -1;
    }
    if warnings & REB_BINARYDATA_ERROR_OUTOFRANGE != 0 {
        reb_simulation_error(r, "Index out of range.");
        return -1;
    }
    if warnings & REB_BINARYDATA_ERROR_SEEK != 0 {
        reb_simulation_error(r, "Error while trying to seek file.");
        return -1;
    }
    if warnings & REB_BINARYDATA_WARNING_FIELD_UNKNOWN != 0 {
        reb_simulation_warning(r, "Unknown field found in binary file.");
    }
    if warnings & REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR != 0 {
        reb_simulation_warning(r, "Custom integrator encountered in Simulationarchive. Call reb_simulation_set_integrator after the simulation is loaded to reset function pointers and initialize data.");
    }
    if warnings & REB_BINARYDATA_ERROR_OLD != 0 {
        reb_simulation_error(r, "Reading old Simulationarchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3");
        return -1;
    }
    if warnings & REB_BINARYDATA_WARNING_CORRUPTFILE != 0 {
        reb_simulation_warning(r, "The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it.");
    }
    0
}

// ---------------------------------------------------------------------
// Diff (binarydata.c `reb_binarydata_diff`)
// ---------------------------------------------------------------------

struct FieldIter<'a> {
    buf: &'a [u8],
    pos: usize,
}

struct FieldRef<'a> {
    name: &'a str,
    data_pos: usize,
    size_data: usize,
}

impl<'a> FieldIter<'a> {
    fn new(buf: &'a [u8]) -> Self {
        FieldIter { buf, pos: 64 }
    }
    /// Read one field header+name; leaves pos at the data.
    fn next_field(&mut self) -> Option<FieldRef<'a>> {
        if self.pos + REB_BINARYDATA_FIELD_SIZE > self.buf.len() {
            return None;
        }
        let size_name = read_u64_le(self.buf, self.pos) as usize;
        let size_data = read_u64_le(self.buf, self.pos + 8) as usize;
        self.pos += REB_BINARYDATA_FIELD_SIZE;
        if self.pos + size_name > self.buf.len() {
            return None;
        }
        let name_bytes = &self.buf[self.pos..self.pos + size_name];
        self.pos += size_name;
        let nul = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
        let name = std::str::from_utf8(&name_bytes[..nul]).unwrap_or("");
        Some(FieldRef {
            name,
            data_pos: self.pos,
            size_data,
        })
    }
}

/// Compare the "particles" payloads like the C (`reb_particle_cmp`,
/// which compares the 11 double members and the name pointer, but not
/// `ap`/`sim`).
fn particles_payload_differ(d1: &[u8], d2: &[u8]) -> bool {
    let n = d1.len() / REB_PARTICLE_RAW_SIZE;
    let mut differ = false;
    for i in 0..n {
        let (p1, n1) = particle_from_raw(&d1[i * REB_PARTICLE_RAW_SIZE..]);
        let (p2, n2) = particle_from_raw(&d2[i * REB_PARTICLE_RAW_SIZE..]);
        differ = differ || (p1.x != p2.x);
        differ = differ || (p1.y != p2.y);
        differ = differ || (p1.z != p2.z);
        differ = differ || (p1.vx != p2.vx);
        differ = differ || (p1.vy != p2.vy);
        differ = differ || (p1.vz != p2.vz);
        differ = differ || (p1.ax != p2.ax);
        differ = differ || (p1.ay != p2.ay);
        differ = differ || (p1.az != p2.az);
        differ = differ || (p1.m != p2.m);
        differ = differ || (p1.r != p2.r);
        differ = differ || (n1 != n2); // name pointer comparison
    }
    differ
}

fn emit_diff_field(out: &mut Vec<u8>, buf: &[u8], f: &FieldRef) {
    let size_name = (f.name.len() + 1) as u64;
    out.extend_from_slice(&size_name.to_le_bytes());
    out.extend_from_slice(&(f.size_data as u64).to_le_bytes());
    out.extend_from_slice(f.name.as_bytes());
    out.push(0);
    out.extend_from_slice(&buf[f.data_pos..f.data_pos + f.size_data]);
}

/// binarydata.c `reb_binarydata_diff` — compares two serialized
/// simulations. Returns (are_different, diff_stream). Only the
/// REB_BINARYDATA_OUTPUT_STREAM and _NONE modes are implemented (the
/// PRINT/BUFFER human-readable modes are console conveniences); STREAM
/// output is what the Simulationarchive append path uses.
pub fn reb_binarydata_diff(
    buf1: &[u8],
    buf2: &[u8],
    output_option: REB_BINARYDATA_OUTPUT,
) -> (i32, Vec<u8>) {
    let mut out: Vec<u8> = Vec::new();
    if buf1.len() < 64 || buf2.len() < 64 {
        println!("Cannot read input buffers.");
        return (0, out);
    }
    let mut are_different = 0;

    // Header comparison (PRINT mode only in the C).
    if buf1[..64] != buf2[..64] && output_option == REB_BINARYDATA_OUTPUT_PRINT {
        println!("Header in binary files are different.");
    }

    // Pass 1: fields of buf1 vs buf2.
    let mut it1 = FieldIter::new(buf1);
    let mut it2 = FieldIter::new(buf2);
    loop {
        let f1 = match it1.next_field() {
            Some(f) => f,
            None => break,
        };
        if f1.name == "end" {
            break;
        }
        let mut f2 = match it2.next_field() {
            Some(f) => f,
            None => {
                it2 = FieldIter::new(buf2);
                match it2.next_field() {
                    Some(f) => f,
                    None => break,
                }
            }
        };
        // Fields might not be in the same order.
        if f1.name != f2.name {
            // Search for the element in buf2, from just past the header.
            it2 = FieldIter::new(buf2);
            let mut notfound = false;
            loop {
                match it2.next_field() {
                    None => {
                        notfound = true;
                        break;
                    }
                    Some(c) => {
                        if c.name == "end" {
                            notfound = true;
                            break;
                        }
                        if c.name == f1.name {
                            f2 = c;
                            break;
                        } else {
                            it2.pos += c.size_data; // skip
                        }
                    }
                }
            }
            if notfound {
                are_different = 1;
                if output_option == REB_BINARYDATA_OUTPUT_STREAM {
                    // The C writes only the header+name here (data
                    // pointer no longer valid on the other side).
                    let size_name = (f1.name.len() + 1) as u64;
                    out.extend_from_slice(&size_name.to_le_bytes());
                    out.extend_from_slice(&(f1.size_data as u64).to_le_bytes());
                    out.extend_from_slice(f1.name.as_bytes());
                    out.push(0);
                }
                // Set offsets for next search
                it2 = FieldIter::new(buf2);
                it1.pos += f1.size_data;
                continue;
            }
        }
        // Same names from here on.
        if f1.data_pos + f1.size_data > buf1.len() {
            println!("Corrupt binary file buf1.");
        }
        if f2.data_pos + f2.size_data > buf2.len() {
            println!("Corrupt binary file buf2.");
        }
        let mut fields_differ = false;
        if f1.size_data == f2.size_data {
            let d1 = &buf1[f1.data_pos..f1.data_pos + f1.size_data];
            let d2 = &buf2[f2.data_pos..f2.data_pos + f2.size_data];
            if f1.name == "particles" {
                fields_differ = particles_payload_differ(d1, d2);
            } else if d1 != d2 {
                fields_differ = true;
            }
        } else {
            fields_differ = true;
        }
        if fields_differ {
            if !f1.name.starts_with("walltime") {
                // Ignore all fields that start with walltime for the
                // return value.
                are_different = 1;
            }
            if output_option == REB_BINARYDATA_OUTPUT_STREAM {
                emit_diff_field(&mut out, buf2, &f2);
            }
        }
        it1.pos = f1.data_pos + f1.size_data;
        it2.pos = f2.data_pos + f2.size_data;
    }

    // Pass 2: fields present in buf2 but not in buf1.
    let mut it1 = FieldIter::new(buf1);
    let mut it2 = FieldIter::new(buf2);
    loop {
        let f2 = match it2.next_field() {
            Some(f) => f,
            None => break,
        };
        if f2.name == "end" {
            break;
        }
        let f1 = match it1.next_field() {
            Some(f) => f,
            None => {
                it1 = FieldIter::new(buf1);
                match it1.next_field() {
                    Some(f) => f,
                    None => break,
                }
            }
        };
        if f1.name == f2.name {
            // Not a new field. Skip.
            it1.pos = f1.data_pos + f1.size_data;
            it2.pos = f2.data_pos + f2.size_data;
            continue;
        }
        // Search for the element in buf1.
        it1 = FieldIter::new(buf1);
        let mut notfound = false;
        loop {
            match it1.next_field() {
                None => {
                    notfound = true;
                    break;
                }
                Some(c) => {
                    if c.name == "end" {
                        notfound = true;
                        break;
                    }
                    if c.name == f2.name {
                        break;
                    } else {
                        it1.pos += c.size_data;
                    }
                }
            }
        }
        if !notfound {
            // Not a new field. Skip.
            it1 = FieldIter::new(buf1);
            it2.pos = f2.data_pos + f2.size_data;
            continue;
        }
        are_different = 1;
        if output_option == REB_BINARYDATA_OUTPUT_STREAM {
            emit_diff_field(&mut out, buf2, &f2);
        }
        it1 = FieldIter::new(buf1);
        it2.pos = f2.data_pos + f2.size_data;
    }

    (are_different, out)
}
