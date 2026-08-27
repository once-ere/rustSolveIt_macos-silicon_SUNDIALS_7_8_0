//! integrator_none.rs — the dummy integrator (rebound.c
//! `reb_integrator_none_step`): advances time and nothing else.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::types::*;

pub fn reb_integrator_none_step(r: &mut reb_simulation) {
    r.t += r.dt;
    r.dt_last_done = r.dt;
}
