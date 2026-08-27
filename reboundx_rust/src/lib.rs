//! reboundx_rs — a pure-Rust translation of REBOUNDx 5.1.0
//! (github.com/dtamayo/reboundx), the library for adding additional
//! forces to the REBOUND N-body integration package.
//!
//! REBOUNDx is an open-source library by Dan Tamayo, Hanno Rein and
//! collaborators, licensed under the GNU General Public License v3 (or
//! later). This translation is a derivative work under the same
//! license; see the LICENSE file carried in this crate. If you use it
//! for published science, cite **Tamayo, Rein, Shi & Hernandez 2019**
//! (MNRAS 491, 2885; arXiv:1908.05634) plus the papers for the
//! individual effects you enable. README.md in this crate lists the
//! paper to cite for every effect.
//!
//! Translation rules (identical to the sibling `rebound_rs` crate):
//! - zero `unsafe`, zero external dependencies, zero warnings;
//! - C function and struct names are preserved (`rebx_attach`,
//!   `rebx_add_force`, `rebx_set_param_double`, ...);
//! - control flow, constants and arithmetic ORDER match the C source
//!   expression for expression (floating point is not associative);
//! - the C pointer graph becomes owned Rust containers — see the
//!   `types` module docs for the three mechanical substitutions used.
//!
//! Deviations from C are all mechanical - parameters carry their
//! type instead of being a `void*` plus a tag, linked lists become
//! vectors whose index 0 is the head (preserving the C's prepend
//! order, which decides the order accelerations are summed),
//! forces and operators are referred to by index instead of by
//! pointer, and the simulation is passed explicitly rather than
//! reached through a `rebx->sim` back-pointer. None of them
//! changes a computed number. README.md lists them with the
//! reasoning under "How this differs from the C".

#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
// The lint set below is waived deliberately and repo-wide: every one of
// these fires on a pattern that mirrors the C source exactly, and
// "fixing" it would either change floating-point evaluation order or
// obscure the correspondence to the C that makes this port reviewable.
// The reason for each one is written beside it below, so a reviewer
// reads it here rather than having to find a configuration file.
#![allow(clippy::too_many_arguments)] // C signatures are preserved verbatim
#![allow(clippy::excessive_precision)] // physical constants carry the C's digits
#![allow(clippy::needless_range_loop)] // index loops mirror the C's for(i=0;...)
#![allow(clippy::identity_op)] // i*6+0 keeps the C's index arithmetic visible
#![allow(clippy::erasing_op)] // 0*n likewise
#![allow(clippy::assign_op_pattern)] // a = a + b mirrors the C statement
#![allow(clippy::collapsible_if)] // nested ifs mirror the C's nesting
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::manual_range_contains)] // (x < lo || x > hi) is the C's test
#![allow(clippy::field_reassign_with_default)] // mirrors `struct X x = {0}; x.a = ..`
#![allow(clippy::needless_late_init)] // mirrors `double t_wave; ... t_wave = ..`
#![allow(clippy::unnecessary_unwrap)] // mirrors `if (i>0 && p != NULL) { *p ... }`
#![allow(clippy::new_without_default)]
#![allow(clippy::type_complexity)]

pub mod types;

pub mod core;
pub mod rebxtools;

// Binary serialization of the REBOUNDx state.
pub mod input;
pub mod output;

pub mod interpolation;
pub mod steppers;

// REBOUNDx's own integrators (used by `integrate_force` and by effects
// that evolve their own ODEs).
pub mod integrator_euler;
pub mod integrator_implicit_midpoint;
pub mod integrator_rk2;
pub mod integrator_rk4;

// Forces.
pub mod central_force;
pub mod exponential_migration;
pub mod gas_damping_timescale;
pub mod gas_dynamical_friction;
pub mod gr;
pub mod gr_full;
pub mod gr_potential;
pub mod gravitational_harmonics;
pub mod inner_disk_edge;
pub mod lense_thirring;
pub mod modify_orbits_forces;
pub mod radiation_forces;
pub mod stochastic_forces;
pub mod tides_constant_time_lag;
pub mod tides_dynamical;
pub mod tides_spin;
pub mod type_I_migration;
pub mod yarkovsky_effect;

// Operators.
pub mod integrate_force;
pub mod modify_mass;
pub mod modify_orbits_direct;
pub mod track_min_distance;

pub use core::*;
pub use rebxtools::*;
pub use types::*;

pub use input::*;
pub use output::*;

pub use interpolation::*;
pub use steppers::*;

pub use integrator_euler::*;
pub use integrator_implicit_midpoint::*;
pub use integrator_rk2::*;
pub use integrator_rk4::*;

pub use central_force::*;
pub use exponential_migration::*;
pub use gas_damping_timescale::*;
pub use gas_dynamical_friction::*;
pub use gr::*;
pub use gr_full::*;
pub use gr_potential::*;
pub use gravitational_harmonics::*;
pub use inner_disk_edge::*;
pub use lense_thirring::*;
pub use modify_orbits_forces::*;
pub use radiation_forces::*;
pub use stochastic_forces::*;
pub use tides_constant_time_lag::*;
pub use tides_dynamical::*;
pub use tides_spin::*;
pub use type_I_migration::*;
pub use yarkovsky_effect::*;

pub use integrate_force::*;
pub use modify_mass::*;
pub use modify_orbits_direct::*;
pub use track_min_distance::*;

/// Version of the C release this crate translates
/// (core.c `rebx_version_str`).
pub const rebx_version_str: &str = "5.1.0";
