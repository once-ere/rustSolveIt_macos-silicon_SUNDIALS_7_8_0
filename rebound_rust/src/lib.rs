//! rebound_rs — a pure-Rust translation of REBOUND 5.1.1
//! (github.com/hannorein/rebound @ dad5f978, "Patch (#931)").
//!
//! REBOUND is an open-source multi-purpose N-body code by Hanno Rein and
//! collaborators, licensed under the GNU General Public License v3 (or
//! later). This translation is a derivative work under the same license;
//! see the LICENSE file carried in this crate. Original authors and
//! copyright holders: Hanno Rein, Shangfei Liu, and the REBOUND
//! contributors.
//!
//! Translation rules (mirroring the sundials_rs porting discipline):
//! - zero `unsafe`, zero external dependencies, zero warnings;
//! - C function and struct names are preserved (`reb_simulation_create`,
//!   `reb_particle`, ...) — the crate root allows the C spellings;
//! - control flow, constants and arithmetic ORDER match the C source
//!   line for line (floating point is not associative);
//! - the glibc `rand_r` generator is reproduced exactly, so random
//!   initial conditions are bit-identical to the C build's;
//! - C's malloc'd pointer graphs become owned Rust containers: the
//!   particle array is a `Vec<reb_particle>`, the octree an index
//!   arena rebuilt each use exactly as the C rebuilds its cells.
//!
//! Deviations from C, all mechanical, are documented in
//! `rebound_rust.md` §"Deviation classes".

#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
//
// ---------------------------------------------------------------------
// Clippy waivers. Every lint below fires on a pattern that deliberately
// mirrors the C source. Applying clippy's suggestion would either change
// floating-point evaluation order (which is not associative, so it would
// change results) or destroy the line-for-line correspondence to the C
// that makes this port reviewable. Each is justified in rebound_rust.md
// section 17.
//
// READ THIS BEFORE "FIXING" ANYTHING BELOW: `neg_cmp_op_on_partial_ord`
// is load-bearing. The port previously hung forever because a C loop
// condition `while (a > b)` had been negated to `if a <= b { break }`.
// Those are NOT equivalent when either side is NaN — which really
// happens for near-rectilinear hyperbolic orbits in the WHFast Kepler
// solver. The negation must stay written as `!(a > b)`. See
// rebound_rust.md section 15.10, "Defect 1".
// ---------------------------------------------------------------------
#![allow(clippy::neg_cmp_op_on_partial_ord)] // MUST STAY: NaN semantics, see above
#![allow(clippy::excessive_precision)] // constants carry the C's exact digits
#![allow(clippy::identity_op)] // `m[0 + 4*0]` keeps C's stride arithmetic visible
#![allow(clippy::erasing_op)] // same: `0 * n` inside C-mirroring index expressions
#![allow(clippy::needless_range_loop)] // `for i in 0..N` mirrors `for(i=0;i<N;i++)`
#![allow(clippy::assign_op_pattern)] // `a = a + b` mirrors the C statement
#![allow(clippy::field_reassign_with_default)] // mirrors C's `struct X x = {0}; x.a = ..`
#![allow(clippy::too_many_arguments)] // C signatures preserved verbatim (HARD RULE 3)
#![allow(clippy::manual_range_contains)] // `x >= a || x < b` mirrors the C test
#![allow(clippy::manual_memcpy)] // explicit copy loops mirror the C's
#![allow(clippy::manual_swap)] // explicit 3-line swaps mirror the C's
#![allow(clippy::manual_div_ceil)] // mirrors the C's `(a + b - 1) / b`
#![allow(clippy::manual_is_multiple_of)] // mirrors the C's `x % n == 0`
#![allow(clippy::misrefactored_assign_op)] // `a = a * b + c` shape is the C's
#![allow(clippy::neg_multiply)] // `-1. * x` appears verbatim in the C
#![allow(clippy::collapsible_if)] // nested `if`s mirror the C's nesting
#![allow(clippy::collapsible_else_if)] // same
#![allow(clippy::needless_late_init)] // C declares at top of block, assigns later
#![allow(clippy::while_let_loop)] // mirrors the C's `for(;;)` with an inner break
#![allow(clippy::unnecessary_cast)] // explicit casts mirror the C's conversions
#![allow(clippy::ptr_arg)] // signature shape mirrors the C's array parameter
#![allow(clippy::seek_from_current)] // mirrors the C's `fseek(f, n, SEEK_CUR)`
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::drop_non_drop)] // explicit scope end mirrors the C's lifetime

pub mod types;
pub mod tools;
pub mod boundary;
pub mod tree;
pub mod gravity;
pub mod collision;
pub mod particle;
pub mod simulation;
pub mod output;
pub mod transformations;
pub mod rotations;
pub mod integrator_none;
pub mod integrator_sei;
pub mod integrator_leapfrog;
pub mod integrator_ias15;
pub mod integrator_whfast;
pub mod integrator_saba;
pub mod integrator_janus;
pub mod integrator_eos;
pub mod integrator_mercurius;
pub mod integrator_bs;
pub mod integrator_trace;
pub mod integrator_whfast512;
pub mod derivatives;
pub mod frequency_analysis;
pub mod binarydata;
pub mod simulationarchive;
pub mod server;

pub use types::*;
pub use tools::*;
pub use boundary::*;
pub use gravity::*;
pub use collision::*;
pub use particle::*;
pub use simulation::*;
pub use output::*;
pub use transformations::*;
pub use rotations::*;
pub use derivatives::*;
pub use frequency_analysis::*;
pub use binarydata::*;
pub use simulationarchive::*;
pub use server::*;

/// Version of the C release this crate translates (rebound.c `reb_version_str`).
pub const reb_version_str: &str = "5.1.1";
/// Git hash of the C source tree the translation was made from.
pub const reb_githash_str: &str = "dad5f97806ecbb408dcaff728851c64e67f9f6eb";
