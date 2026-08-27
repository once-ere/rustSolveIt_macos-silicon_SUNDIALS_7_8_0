//! integrator_whfast512.rs — the AVX512-accelerated WHFast512
//! integrator (from integrator_whfast512.c/h; (c) 2026 Hanno Rein,
//! Rishit Dagli, Pejvak Javaheri).
//!
//! The computational core of WHFast512 is hand-written AVX-512 assembly
//! (integrator_whfast512.s) compiled only on 64-bit GCC/Clang targets.
//! The reference build for this port — REBOUND compiled with MSVC `cl`
//! on Windows — takes the `#else // Not 64 bit, Windows + cl` branch,
//! which contains only stubs that raise an error. This module mirrors
//! that reference build exactly: the state struct and its defaults are
//! carried (`create`), and `step`/`synchronize` emit the same error and
//! set the same status the Windows C build does.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use crate::tools::reb_simulation_error;
use crate::types::*;

/// integrator_whfast512.h `struct reb_integrator_whfast512_state`. The
/// internal `void* data` (the 64-byte-aligned AVX-512 SIMD block) is
/// never allocated on the Windows reference build and is not carried.
#[derive(Clone, Debug)]
pub struct reb_integrator_whfast512_state {
    /// 1: include general relativistic corrections (assumes AU, yr/2pi).
    pub gr_potential: u32,
    /// Number of independent planetary systems (1, 2, or 4).
    pub N_systems: u32,
    /// 17: use symplectic correctors (default 0).
    pub corrector: u32,
    /// Number of timesteps combined into one call (default 1e6).
    pub concatenate_steps: u64,
    // Internal use
    pub last_synchronization: f64,
}

impl Default for reb_integrator_whfast512_state {
    /// integrator_whfast512.c `reb_integrator_whfast512_create`.
    fn default() -> Self {
        reb_integrator_whfast512_state {
            gr_potential: 0,
            N_systems: 1,
            corrector: 0,
            concatenate_steps: 1e6 as u64,
            last_synchronization: 0.,
        }
    }
}

/// integrator_whfast512.c `reb_integrator_whfast512_step` — the
/// `#else // Not 64 bit, Windows + cl` branch of the reference build.
pub fn reb_integrator_whfast512_step(r: &mut reb_simulation) {
    reb_simulation_error(r, "AVX512 is not supported on your platform.");
    r.status = REB_STATUS_GENERIC_ERROR;
}

/// integrator_whfast512.c `reb_integrator_whfast512_synchronize` —
/// same stub branch.
pub fn reb_integrator_whfast512_synchronize(r: &mut reb_simulation) {
    reb_simulation_error(r, "AVX512 is not supported on your platform.");
    r.status = REB_STATUS_GENERIC_ERROR;
}
