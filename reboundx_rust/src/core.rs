//! core.rs — translation of REBOUNDx core.c
//! Central internal functions for REBOUNDx (not called by user): the
//! lifecycle of the `rebx_extras` state, the force/operator library, the
//! parameter store, and the callbacks REBOUND invokes each timestep.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # What this module is
//!
//! `core.c` is the keystone of REBOUNDx. It holds
//!
//! * `rebx_attach` / `rebx_detach` / `rebx_free` — the lifecycle,
//! * `rebx_register_default_params` — the ~110 parameter names every
//!   effect in the library reads or writes, with their types,
//! * `rebx_load_force` / `rebx_load_operator` — the name → function
//!   tables that turn `"gr"` into `rebx_gr` with `REBX_FORCE_VEL`,
//! * `rebx_add_force` / `rebx_add_operator` / `rebx_add_operator_step`
//!   — installing effects into the simulation's callback slots,
//! * the `rebx_get_param_*` / `rebx_set_param_*` parameter interface,
//! * `rebx_additional_forces`, `rebx_pre_timestep_modifications` and
//!   `rebx_post_timestep_modifications` — the three functions REBOUND
//!   itself calls, which walk REBOUNDx's lists in order.
//!
//! # The three mechanical substitutions (see `types` module docs)
//!
//! 1. C linked lists become `Vec`s whose **index 0 is the head**. The C
//!    `rebx_add_node` prepends, so every list that the C prepends to is
//!    built here with `insert(0, ..)`: the traversal order — and hence
//!    the order in which accelerations are summed — is identical.
//!
//!    The one exception is `allocated_forces` / `allocated_operators`.
//!    There, the C's node address *is* the object's identity; here the
//!    **Vec index is the identity**, handed back by `rebx_create_force`
//!    and passed to every later call. Prepending would renumber every
//!    previously issued index, so those two vectors are appended to and
//!    the C's traversal order is restored where it is observable, by
//!    searching them in reverse (`rebx_get_force`, `rebx_get_operator`).
//!    Neither list's order affects any floating-point result: they exist
//!    in the C purely for `free()`ing and for name lookup.
//!
//! 2. The `void* value` + type tag pair becomes the single
//!    `rebx_param_value` enum, so a wrong-type read is impossible.
//!
//! 3. `struct rebx_node** apptr` (always the address of some object's
//!    `ap` member) becomes `rebx_ap`, which names the same list by index.
//!
//! # Taking the extras out of the simulation
//!
//! The C effect functions receive `sim` and reach REBOUNDx through
//! `sim->extras`; that is two live pointers into overlapping state,
//! which safe Rust forbids. Every entry point here that needs both
//! therefore *moves* the box out of `sim.extras`, works with the two
//! halves side by side, and puts it back. Every early return restores
//! it — a lost extras box would be a silent, hard-to-find bug.

use rebound_rs::{
    reb_integrator_state, reb_orbit, reb_particle, reb_simulation, reb_simulation_error,
    reb_simulation_warning, reb_vec3d,
};

use crate::types::rebx_force_type::*;
use crate::types::rebx_operator_type::*;
use crate::types::rebx_param_type::*;
use crate::types::*;

// core.c's `rebx_version_str` lives at the crate root (lib.rs), where
// `rebx_build_str` and `rebx_githash_str` — compile-time __DATE__ /
// __TIME__ / git-hash strings with no meaning after translation — are
// deliberately not carried.

/*****************************
 Attaching / detaching
 ****************************/

/// Borrow the REBOUNDx state held in `sim.extras` (C: `sim->extras`).
///
/// Returns `None` when `rebx_attach` was never called on this
/// simulation, which is the C's `sim->extras == NULL`.
pub fn rebx_extras_ref(sim: &reb_simulation) -> Option<&rebx_extras> {
    match sim.extras.as_ref() {
        None => None,
        Some(b) => b.downcast_ref::<rebx_extras>(),
    }
}

/// Mutably borrow the REBOUNDx state held in `sim.extras`.
pub fn rebx_extras_mut(sim: &mut reb_simulation) -> Option<&mut rebx_extras> {
    match sim.extras.as_mut() {
        None => None,
        Some(b) => b.downcast_mut::<rebx_extras>(),
    }
}

/// Move the extras box out of `sim` so that `sim` and the REBOUNDx
/// state can be borrowed mutably at the same time. Always paired with
/// `rebx_put`. If `sim.extras` holds something that is not a
/// `rebx_extras`, it is put straight back and `None` is returned.
fn rebx_take(sim: &mut reb_simulation) -> Option<Box<rebx_extras>> {
    match sim.extras.take() {
        None => None,
        Some(b) => match b.downcast::<rebx_extras>() {
            Ok(rebx) => Some(rebx),
            Err(other) => {
                sim.extras = Some(other);
                None
            }
        },
    }
}

/// Put back what `rebx_take` removed.
fn rebx_put(sim: &mut reb_simulation, rebx: Box<rebx_extras>) {
    sim.extras = Some(rebx);
}

/// Run `f` with the simulation and the REBOUNDx state both borrowed
/// mutably at once.
///
/// In C both are reachable at the same time because `sim->extras` is a
/// raw pointer. Safe Rust cannot hand out `&mut sim` and a `&mut` into
/// `sim.extras` simultaneously, so this helper moves the state out for
/// the duration of the call and puts it back afterwards — the same
/// take/use/put pattern the integrator states use.
///
/// Use it for the REBOUNDx functions that need both, such as
/// [`crate::rebxtools::rebx_tools_spin_angular_momentum`] and
/// [`crate::rebxtools::rebx_simulation_irotate`]:
///
/// ```ignore
/// let L = rebx_with(&mut sim, |sim, rebx| {
///     rebx_tools_spin_angular_momentum(sim, rebx)
/// }).unwrap();
/// ```
///
/// Returns `None` (without calling `f`) if no REBOUNDx state is
/// attached, which is the situation where the C would dereference a
/// NULL `sim->extras`.
pub fn rebx_with<R>(
    sim: &mut reb_simulation,
    f: impl FnOnce(&mut reb_simulation, &mut rebx_extras) -> R,
) -> Option<R> {
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => {
            rebx_error_detached();
            return None;
        }
    };
    let out = f(sim, &mut rebx);
    rebx_put(sim, rebx);
    Some(out)
}

/// C: `rebx_error` when `rebx->sim == NULL`. Reached here when the
/// simulation has no REBOUNDx state attached at all.
fn rebx_error_detached() {
    eprintln!(
        "\nError! REBOUNDx Error: A Simulation is no longer attached to this REBOUNDx extras \
         instance. Most likely the Simulation has been freed."
    );
}

/// reboundx.h `rebx_attach`.
///
/// The C mallocs a `rebx_extras`, calls `rebx_initialize` (which stores
/// it in `sim->extras`, NULLs every list, installs the
/// `free_particle_ap` / `extras_cleanup` callbacks and warns about
/// callback slots that are already in use), then registers the default
/// parameters. The order below is the C's.
///
/// The C returns the new `struct rebx_extras*`; here the state lives in
/// `sim.extras` and is reached with `rebx_extras_mut` / `rebx_extras_ref`,
/// so there is nothing to hand back.
pub fn rebx_attach(sim: &mut reb_simulation) {
    rebx_initialize(sim);

    if let Some(rebx) = rebx_extras_mut(sim) {
        rebx_register_default_params(rebx);
    }
}

/// core.c `rebx_initialize` — attaches an empty REBOUNDx state to the
/// simulation, *without* registering the default parameters.
///
/// The C takes the freshly malloc'd `struct rebx_extras*` and stores it
/// in `sim->extras`; here the state is created and stored in one step.
/// `rebx_attach` calls this and then registers the defaults;
/// `rebx_create_extras_from_binary` calls it alone, because a binary
/// file carries the registered-parameter list that was in effect when
/// it was written.
pub fn rebx_initialize(sim: &mut reb_simulation) {
    // C: sim->extras = rebx; all lists NULL.
    // `rebx_extras::default()` gives the same empty lists.
    //
    // sim->free_particle_ap and sim->extras_cleanup have no counterpart:
    // particle parameter lists live in `rebx_extras::particle_params`
    // (they are not owned by `reb_particle`), and Rust drops the extras
    // box with the simulation, so there is nothing for REBOUND to call
    // back into.
    sim.extras = Some(Box::new(rebx_extras::default()));

    if sim.additional_forces.is_some()
        || sim.pre_timestep_modifications.is_some()
        || sim.post_timestep_modifications.is_some()
    {
        reb_simulation_warning(sim, "REBOUNDx overwrites sim->additional_forces, sim->pre_timestep_modifications and sim->post_timestep_modifications whenever forces or operators that use them get added.  If you want to use REBOUNDx together with your own custom functions that use these callbacks, you should add them through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.");
    }
}

/// core.c `rebx_detach`. Clears the three callback slots if — and only
/// if — they still point at REBOUNDx's own functions, then drops the
/// extras.
///
/// The C additionally NULLs `sim->extras_cleanup` and
/// `sim->free_particle_ap`; neither exists here (see `rebx_attach`).
pub fn rebx_detach(sim: &mut reb_simulation) {
    if sim.extras.is_none() {
        return;
    }
    let rebx_af: fn(&mut reb_simulation) = rebx_additional_forces;
    let rebx_pre: fn(&mut reb_simulation) = rebx_pre_timestep_modifications;
    let rebx_post: fn(&mut reb_simulation) = rebx_post_timestep_modifications;

    if let Some(f) = sim.additional_forces {
        if std::ptr::fn_addr_eq(f, rebx_af) {
            sim.additional_forces = None;
        }
    }
    if let Some(f) = sim.pre_timestep_modifications {
        if std::ptr::fn_addr_eq(f, rebx_pre) {
            sim.pre_timestep_modifications = None;
        }
    }
    if let Some(f) = sim.post_timestep_modifications {
        if std::ptr::fn_addr_eq(f, rebx_post) {
            sim.post_timestep_modifications = None;
        }
    }
    sim.extras = None;
}

/// core.c `rebx_free`.
///
/// The C calls `rebx_free_pointers` (walking every list and `free()`ing
/// each node, parameter, force, operator and step) and then
/// `rebx_detach`. Here every one of those allocations is owned by the
/// `rebx_extras` box, so dropping the box frees all of it; only the
/// detaching — unhooking the callbacks — carries real work. Nothing is
/// faked: there simply is no manual free to perform.
pub fn rebx_free(sim: &mut reb_simulation) {
    rebx_detach(sim);
}

/// core.c `rebx_free_ap`. Empties one parameter list.
///
/// In the C this frees the nodes and their parameters; here `clear()`
/// drops them. As in the C, it is safe to call twice.
pub fn rebx_free_ap(rebx: &mut rebx_extras, sel: rebx_ap) {
    if let Some(ap) = rebx.ap_mut(sel) {
        ap.clear();
    }
}

/// core.c `rebx_free_particle_ap`. Takes the particle's index, since
/// particle parameter lists live in `rebx_extras::particle_params`.
pub fn rebx_free_particle_ap(rebx: &mut rebx_extras, particle_index: usize) {
    rebx_free_ap(rebx, rebx_ap::particle(particle_index));
}

/*****************************
 Default parameter registration
 ****************************/

/// core.c `rebx_register_default_params`. Every parameter name that any
/// effect in the library reads or writes, with its type, in the C's
/// order.
pub fn rebx_register_default_params(rebx: &mut rebx_extras) {
    rebx_register_param(rebx, "c", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gr_source", REBX_TYPE_INT);
    rebx_register_param(rebx, "tau_mass", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "force", REBX_TYPE_FORCE);
    rebx_register_param(rebx, "particle", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "Acentral", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gammacentral", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "max_iterations", REBX_TYPE_INT);
    rebx_register_param(rebx, "J2", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "J4", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "R_eq", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "coordinates", REBX_TYPE_INT);
    rebx_register_param(rebx, "p", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "d_factor", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "cs_coeff", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_coeff", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_a", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_e", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_inc", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_omega", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_Omega", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "em_tau_a", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "em_aini", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "em_afin", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "primary", REBX_TYPE_INT);
    rebx_register_param(rebx, "radiation_source", REBX_TYPE_INT);
    rebx_register_param(rebx, "kappa", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "kappa_x", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "kappa_y", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "kappa_z", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_kappa", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_kappa_x", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_kappa_y", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau_kappa_z", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "stochastic_force_r", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "stochastic_force_phi", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "stochastic_force_x", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "stochastic_force_y", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "stochastic_force_z", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "beta", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tides_primary", REBX_TYPE_INT);
    rebx_register_param(rebx, "R_tides", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tctl_k2", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tctl_tau", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "integrator", REBX_TYPE_INT);
    rebx_register_param(rebx, "im_ps_final", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "im_ps_prev", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "im_ps_avg", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "rk2_k2", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "rk4_k2", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "rk4_k3", REBX_TYPE_POINTER);
    rebx_register_param(rebx, "min_distance", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "min_distance_from", REBX_TYPE_STRING);
    rebx_register_param(rebx, "min_distance_orbit", REBX_TYPE_ORBIT);
    rebx_register_param(rebx, "luminosity", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ide_position", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ide_width", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tIm_flaring_index", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tIm_scale_height_1", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tIm_surface_density_1", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tIm_surface_density_exponent", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_c", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_body_density", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_lstar", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_flag", REBX_TYPE_INT);
    rebx_register_param(rebx, "ye_rotation_period", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_thermal_inertia", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_albedo", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_emissivity", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_k", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_stef_boltz", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_spin_axis_x", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_spin_axis_y", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ye_spin_axis_z", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "OmegaMag", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "Omega", REBX_TYPE_VEC3D);
    rebx_register_param(rebx, "k2", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "I", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "tau", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "ode", REBX_TYPE_ODE);
    rebx_register_param(rebx, "gas_df_rhog", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_alpha_rhog", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_cs", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_alpha_cs", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_xmin", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_hr", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gas_df_Qd", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_R_eq", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_Mom_I_fac", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_rot_rate", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_p_hatx", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_p_haty", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_p_hatz", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "lt_c", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_M_last", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_num_apoapsis", REBX_TYPE_INT);
    rebx_register_param(rebx, "td_c_imag", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_c_real", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_dP_hat", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_dP_crit", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_EB0", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_E_max", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_E_resid", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_dE_last", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_last_apoapsis", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_drag_coef", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "td_disruption_flag", REBX_TYPE_INT);
}

/// core.c `rebx_register_param`.
pub fn rebx_register_param(rebx: &mut rebx_extras, name: &str, type_: rebx_param_type) {
    // check registered_params for entry
    let reg_type = rebx_get_type(rebx, name);

    if reg_type != REBX_TYPE_NONE {
        let str = format!(
            "REBOUNDx Error: Parameter name '{}' already in registered list. Cannot add duplicates.\n",
            name
        );
        rebx_error(rebx, &str);
        return;
    }

    // Create new entry. These are just rebx_param structs without value populated
    let param = rebx_create_param(name, type_);
    rebx_add_param(&mut rebx.registered_params, param);
}

/**********************************************
 User Interface for adding forces and operators
 *********************************************/

/// reboundx.h `rebx_create_force`.
///
/// The C returns the new `struct rebx_force*`; here the identity of a
/// force is its index into `rebx_extras::allocated_forces`, which is
/// what is returned. Allocation cannot fail, so unlike the C there is
/// no NULL path.
pub fn rebx_create_force(rebx: &mut rebx_extras, name: &str) -> usize {
    let force = rebx_force {
        name: name.to_string(),
        ap: Vec::new(),
        force_type: REBX_FORCE_NONE,
        update_accelerations: None,
    };

    // Add force to allocated_forces list for later freeing. Appended,
    // not prepended: the index is the force's identity and must not
    // shift. See the module docs; rebx_get_force restores the C's
    // traversal order by searching in reverse.
    rebx.allocated_forces.push(force);

    rebx.allocated_forces.len() - 1
}

/// reboundx.h `rebx_load_force` — the name → (function, force_type)
/// table of every force in the REBOUNDx library.
pub fn rebx_load_force(sim: &mut reb_simulation, name: &str) -> Option<usize> {
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => {
            rebx_error_detached();
            return None;
        }
    };
    let force = rebx_create_force(&mut rebx, name);

    let (update_accelerations, force_type): (rebx_force_fn, rebx_force_type) = if name == "gr" {
        (crate::gr::rebx_gr, REBX_FORCE_VEL)
    } else if name == "central_force" {
        (crate::central_force::rebx_central_force, REBX_FORCE_POS)
    } else if name == "modify_orbits_forces" {
        (
            crate::modify_orbits_forces::rebx_modify_orbits_forces,
            REBX_FORCE_VEL,
        )
    } else if name == "gas_damping_timescale" {
        (
            crate::gas_damping_timescale::rebx_gas_damping_timescale,
            REBX_FORCE_VEL,
        )
    } else if name == "exponential_migration" {
        (
            crate::exponential_migration::rebx_exponential_migration,
            REBX_FORCE_VEL,
        )
    } else if name == "gr_full" {
        (crate::gr_full::rebx_gr_full, REBX_FORCE_VEL)
    } else if name == "gravitational_harmonics" {
        (
            crate::gravitational_harmonics::rebx_gravitational_harmonics,
            REBX_FORCE_POS,
        )
    } else if name == "gr_potential" {
        (crate::gr_potential::rebx_gr_potential, REBX_FORCE_POS)
    } else if name == "radiation_forces" {
        (
            crate::radiation_forces::rebx_radiation_forces,
            REBX_FORCE_VEL,
        )
    } else if name == "stochastic_forces" {
        (
            crate::stochastic_forces::rebx_stochastic_forces,
            REBX_FORCE_VEL,
        )
    } else if name == "tides_constant_time_lag" {
        (
            crate::tides_constant_time_lag::rebx_tides_constant_time_lag,
            REBX_FORCE_VEL,
        )
    } else if name == "type_I_migration" {
        (
            crate::type_I_migration::rebx_modify_orbits_with_type_I_migration,
            REBX_FORCE_VEL,
        )
    } else if name == "tides_spin" {
        reb_simulation_warning(sim, "tides_spin was updated in version 4.5.0 to halve the acceleration from the conservative piece of the tidal potential, reflecting a typo discovered in Eggleton et. al (1998). This warning will be removed in a future version.\n");
        (crate::tides_spin::rebx_tides_spin, REBX_FORCE_VEL)
    } else if name == "yarkovsky_effect" {
        (
            crate::yarkovsky_effect::rebx_yarkovsky_effect,
            REBX_FORCE_VEL,
        )
    } else if name == "gas_dynamical_friction" {
        (
            crate::gas_dynamical_friction::rebx_gas_dynamical_friction,
            REBX_FORCE_VEL,
        )
    } else if name == "lense_thirring" {
        (crate::lense_thirring::rebx_lense_thirring, REBX_FORCE_VEL)
    } else if name == "tides_dynamical" {
        (crate::tides_dynamical::rebx_tides_dynamical, REBX_FORCE_VEL)
    } else {
        let str = format!(
            "REBOUNDx error: Force '{}' not found in REBOUNDx library.\n",
            name
        );
        rebx_error(&mut rebx, &str);
        rebx_remove_force(&mut rebx, force); // Not free_force. Must remove from allocated_forces
        rebx_put(sim, rebx);
        return None;
    };

    rebx.allocated_forces[force].update_accelerations = Some(update_accelerations);
    rebx.allocated_forces[force].force_type = force_type;

    rebx_put(sim, rebx);
    Some(force)
}

/// reboundx.h `rebx_create_operator`. Returns the operator's index into
/// `rebx_extras::allocated_operators` (the C returns a pointer).
pub fn rebx_create_operator(rebx: &mut rebx_extras, name: &str) -> usize {
    let operator_ = rebx_operator {
        name: name.to_string(),
        ap: Vec::new(),
        operator_type: REBX_OPERATOR_NONE,
        step_function: None,
    };

    // Add operator to allocated_operators list for later freeing.
    // Appended for the same reason as allocated_forces above.
    rebx.allocated_operators.push(operator_);

    rebx.allocated_operators.len() - 1
}

/// reboundx.h `rebx_load_operator` — the name → (function,
/// operator_type) table of every operator in the REBOUNDx library.
pub fn rebx_load_operator(sim: &mut reb_simulation, name: &str) -> Option<usize> {
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => {
            rebx_error_detached();
            return None;
        }
    };
    let operator_ = rebx_create_operator(&mut rebx, name);

    let (step_function, operator_type): (rebx_operator_fn, rebx_operator_type) =
        if name == "modify_mass" {
            (crate::modify_mass::rebx_modify_mass, REBX_OPERATOR_UPDATER)
        } else if name == "integrate_force" {
            (
                crate::integrate_force::rebx_integrate_force,
                REBX_OPERATOR_UPDATER,
            )
        } else if name == "drift" {
            (crate::steppers::rebx_drift_step, REBX_OPERATOR_UPDATER)
        } else if name == "kick" {
            (crate::steppers::rebx_kick_step, REBX_OPERATOR_UPDATER)
        } else if name == "kepler" {
            (crate::steppers::rebx_kepler_step, REBX_OPERATOR_UPDATER)
        } else if name == "jump" {
            (crate::steppers::rebx_jump_step, REBX_OPERATOR_UPDATER)
        } else if name == "interaction" {
            (crate::steppers::rebx_interaction_step, REBX_OPERATOR_UPDATER)
        } else if name == "ias15" {
            (crate::steppers::rebx_ias15_step, REBX_OPERATOR_UPDATER)
        } else if name == "modify_orbits_direct" {
            (
                crate::modify_orbits_direct::rebx_modify_orbits_direct,
                REBX_OPERATOR_UPDATER,
            )
        } else if name == "track_min_distance" {
            (
                crate::track_min_distance::rebx_track_min_distance,
                REBX_OPERATOR_RECORDER,
            )
        } else {
            let str = format!(
                "REBOUNDx error: Operator '{}' not found in REBOUNDx library.\n",
                name
            );
            rebx_error(&mut rebx, &str);
            rebx_remove_operator(&mut rebx, operator_); // Not free_op. Must rm from allocated_forces
            rebx_put(sim, rebx);
            return None;
        };

    rebx.allocated_operators[operator_].step_function = Some(step_function);
    rebx.allocated_operators[operator_].operator_type = operator_type;

    rebx_put(sim, rebx);
    Some(operator_)
}

/// reboundx.h `rebx_add_force`. Returns 1 on success, 0 on failure,
/// exactly as the C.
///
/// Note the C does *not* refuse a force that has already been added:
/// adding the same force twice makes it contribute twice. That is
/// preserved here.
pub fn rebx_add_force(sim: &mut reb_simulation, force_idx: usize) -> i32 {
    if sim.integrator.name() == "whfast512" {
        reb_simulation_error(sim, "REBOUNDx Error: WHFast512 has been optimized for speed with options stripped out. This integrator will never be compatible with REBOUNDx.\n");
        return 0;
    }
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => {
            rebx_error_detached();
            return 0;
        }
    };

    // C: force == NULL. Here: an index that names no allocated force.
    if force_idx >= rebx.allocated_forces.len() {
        rebx_error(
            &mut rebx,
            "REBOUNDx error: Passed NULL pointer to rebx_add_force.\n",
        );
        rebx_put(sim, rebx);
        return 0;
    }

    if rebx.allocated_forces[force_idx].update_accelerations.is_none() {
        rebx_error(&mut rebx, "REBOUNDx error: Need to set update_accelerations function pointer on force before calling rebx_add_force. See custom effects example.\n");
        rebx_put(sim, rebx);
        return 0;
    }

    if rebx.allocated_forces[force_idx].force_type == REBX_FORCE_NONE {
        rebx_error(&mut rebx, "REBOUNDx error: Need to set force_type field on force before calling rebx_add_force. See custom effects example.\n");
        rebx_put(sim, rebx);
        return 0;
    }

    if rebx.allocated_forces[force_idx].force_type == REBX_FORCE_VEL {
        sim.force_is_velocity_dependent = 1;
    }

    // Could add logic based on different integrators
    // The C prepends, so the list traverses in reverse order of addition
    // and the accelerations sum in that order. insert(0, ..) keeps it.
    rebx.additional_forces.insert(0, force_idx);
    rebx_put(sim, rebx);

    let rebx_af: fn(&mut reb_simulation) = rebx_additional_forces;
    if let Some(f) = sim.additional_forces {
        if !std::ptr::fn_addr_eq(f, rebx_af) {
            reb_simulation_warning(sim, "REBOUNDx Warning: additional_forces was set and is being overwritten by REBOUNDx. To incorporate both, you can add your own custom effects through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.\n");
        }
    }
    sim.additional_forces = Some(rebx_af);

    1
}

/// reboundx.h `rebx_add_operator_step`.
pub fn rebx_add_operator_step(
    sim: &mut reb_simulation,
    operator_idx: usize,
    dt_fraction: f64,
    timing: rebx_timing,
) -> i32 {
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => {
            rebx_error_detached();
            return 0;
        }
    };
    // C: operator == NULL.
    if operator_idx >= rebx.allocated_operators.len() {
        rebx_error(
            &mut rebx,
            "REBOUNDx error: Passed NULL pointer to rebx_add_operator_step.\n",
        );
        rebx_put(sim, rebx);
        return 0;
    }
    if rebx.allocated_operators[operator_idx].step_function.is_none() {
        rebx_error(&mut rebx, "REBOUNDx error: Need to set step_function pointer on operator before adding to simulation. See custom effects example.\n");
        rebx_put(sim, rebx);
        return 0;
    }

    if rebx.allocated_operators[operator_idx].operator_type == REBX_OPERATOR_NONE {
        rebx_error(&mut rebx, "REBOUNDx error: Need to set operator_type field on operator before adding to simulation. See custom effects example.\n");
        rebx_put(sim, rebx);
        return 0;
    }

    let step = rebx_step {
        operator_: operator_idx,
        dt_fraction,
    };

    match timing {
        rebx_timing::REBX_TIMING_PRE => {
            rebx.pre_timestep_modifications.insert(0, step);
            rebx_put(sim, rebx);
            let rebx_pre: fn(&mut reb_simulation) = rebx_pre_timestep_modifications;
            if let Some(f) = sim.pre_timestep_modifications {
                if !std::ptr::fn_addr_eq(f, rebx_pre) {
                    reb_simulation_warning(sim, "REBOUNDx Warning: pre_timestep_modifications was set in the simulation and is being overwritten by REBOUNDx. To incorporate both, you can add your own custom effects through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.\n");
                }
            }
            sim.pre_timestep_modifications = Some(rebx_pre);
            1
        }
        rebx_timing::REBX_TIMING_POST => {
            rebx.post_timestep_modifications.insert(0, step);
            rebx_put(sim, rebx);
            let rebx_post: fn(&mut reb_simulation) = rebx_post_timestep_modifications;
            if let Some(f) = sim.post_timestep_modifications {
                if !std::ptr::fn_addr_eq(f, rebx_post) {
                    reb_simulation_warning(sim, "REBOUNDx Warning: post_timestep_modifications was set in the simulation and is being overwritten by REBOUNDx. To incorporate both, you can add your own custom effects through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.\n");
                }
            }
            sim.post_timestep_modifications = Some(rebx_post);
            1
        }
    }
    // The C has a trailing `return 0` for a timing value that is neither
    // PRE nor POST; the enum has exactly those two variants, so that
    // path is unreachable here.
}

/// reboundx.h `rebx_add_operator`. Picks the timing and dt_fraction that
/// suit the simulation's current integrator, then delegates to
/// `rebx_add_operator_step`.
pub fn rebx_add_operator(sim: &mut reb_simulation, operator_idx: usize) -> i32 {
    let operator_type = {
        let mut rebx = match rebx_take(sim) {
            Some(r) => r,
            None => {
                rebx_error_detached();
                return 0;
            }
        };
        // C: operator == NULL.
        if operator_idx >= rebx.allocated_operators.len() {
            rebx_error(
                &mut rebx,
                "REBOUNDx error: Passed NULL pointer to rebx_add_operator.\n",
            );
            rebx_put(sim, rebx);
            return 0;
        }
        let operator_type = rebx.allocated_operators[operator_idx].operator_type;
        rebx_put(sim, rebx);
        operator_type
    };

    if operator_type == REBX_OPERATOR_RECORDER {
        // Doesn't alter state. Add once after timestep.
        let dt_fraction = 1.;
        let success = rebx_add_operator_step(
            sim,
            operator_idx,
            dt_fraction,
            rebx_timing::REBX_TIMING_POST,
        );
        return success;
    }

    reb_simulation_warning(sim, "REBOUNDx Warning: Do not change the integrator after adding an operator (function assumed the current integrator when adding)");

    let integrator_name = sim.integrator.name();

    if integrator_name == "ias15" || integrator_name == "bs" {
        // don't add pre-timestep b/c don't know what IAS/BS will choose as dt
        let dt_fraction = 1.;
        let success = rebx_add_operator_step(
            sim,
            operator_idx,
            dt_fraction,
            rebx_timing::REBX_TIMING_POST,
        );
        return success;
    }
    if integrator_name == "whfast"
        || integrator_name == "saba"
        || integrator_name == "leapfrog"
        || integrator_name == "eos"
    {
        // half step pre and post
        let dt_fraction = 1. / 2.;
        let success1 = rebx_add_operator_step(
            sim,
            operator_idx,
            dt_fraction,
            rebx_timing::REBX_TIMING_PRE,
        );
        let success2 = rebx_add_operator_step(
            sim,
            operator_idx,
            dt_fraction,
            rebx_timing::REBX_TIMING_POST,
        );
        return if success1 != 0 && success2 != 0 { 1 } else { 0 };
    }
    if integrator_name == "mercurius"
        || integrator_name == "trace"
        || integrator_name == "janus"
        || integrator_name == "sei"
    {
        // TODO: Not yet implemented.
        if operator_type == REBX_OPERATOR_UPDATER {
            reb_simulation_error(sim, "REBOUNDx Error: Operators that affect particle trajectories are not supported with MERCURIUS, TRACE, SEI or Janus. Can only add forces.\n");
            return 0;
        }
    }
    if integrator_name == "whfast512" {
        reb_simulation_error(sim, "REBOUNDx Error: WHFast512 has been optimized for speed with options stripped out. This integrator will never be compatible with REBOUNDx.\n");
        return 0;
    }
    0 // didn't reach a successful outcome
}

/*****************************************************************
 User interface for setting parameter values
 *****************************************************************/

/// core.c `rebx_get_or_add_param`. Gets the parameter's position in the
/// list named by `sel` if it already exists, otherwise creates it (with
/// the type it was registered with) and prepends it, as the C does.
///
/// Returns the index of the parameter inside that list, or `None` where
/// the C returns NULL.
pub fn rebx_get_or_add_param(
    rebx: &mut rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<usize> {
    // C: `if (apptr == NULL)`. A rebx_ap always names a list, so the
    // only way to name nothing is an index past the end of
    // allocated_forces / allocated_operators.
    let bad_selector = match sel {
        rebx_ap::particle(_) => false,
        rebx_ap::force(i) => i >= rebx.allocated_forces.len(),
        rebx_ap::operator_(i) => i >= rebx.allocated_operators.len(),
    };
    if bad_selector {
        rebx_error(
            rebx,
            "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n",
        );
        return None;
    }

    // Check whether it already exists in linked list
    if let Some(i) = rebx.ap(sel).iter().position(|p| p.name == param_name) {
        return Some(i);
    }

    let type_ = rebx_get_type(rebx, param_name);
    if type_ == REBX_TYPE_NONE {
        let str = format!(
            "REBOUNDx Error: Need to register parameter name '{}' before using it. See examples.\n",
            param_name
        );
        rebx_error(rebx, &str);
        return None;
    }
    let param = rebx_create_param(param_name, type_);
    match rebx.ap_mut(sel) {
        Some(ap) => {
            rebx_add_param(ap, param);
            Some(0) // rebx_add_param prepends, so the new param is at the head
        }
        None => None,
    }
}

/// Shared tail of every `rebx_set_param_*`: store `value` on the
/// parameter named by `param_name` in the list named by `sel`,
/// creating it first if needed.
fn rebx_set_param_value(
    rebx: &mut rebx_extras,
    sel: rebx_ap,
    param_name: &str,
    value: rebx_param_value,
) {
    let idx = match rebx_get_or_add_param(rebx, sel, param_name) {
        Some(i) => i,
        None => return,
    };
    if let Some(ap) = rebx.ap_mut(sel) {
        ap[idx].value = value;
    }
}

/// reboundx.h `rebx_set_param_double`.
pub fn rebx_set_param_double(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: f64) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::double(val));
}

/// reboundx.h `rebx_set_param_int`.
pub fn rebx_set_param_int(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: i32) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::int(val));
}

/// reboundx.h `rebx_set_param_uint32`.
pub fn rebx_set_param_uint32(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: u32) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::uint32(val));
}

/// reboundx.h `rebx_set_param_vec3d`.
///
/// The C copies x, y and z field by field into the malloc'd vector; the
/// value is the same either way.
pub fn rebx_set_param_vec3d(
    rebx: &mut rebx_extras,
    sel: rebx_ap,
    param_name: &str,
    val: reb_vec3d,
) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::vec3d(val));
}

/// reboundx.h `rebx_set_param_orbit` (`REBX_TYPE_ORBIT`, set through
/// the generic pointer setter in the C).
pub fn rebx_set_param_orbit(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: reb_orbit) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::orbit(val));
}

/// reboundx.h `rebx_set_param_string`.
///
/// The C hands the string to REBOUND (`reb_simulation_register_name`)
/// so that REBOUND owns the memory; here the `String` is owned by the
/// parameter itself, so no registration is needed. The stored text is
/// the same.
pub fn rebx_set_param_string(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: &str) {
    rebx_set_param_value(
        rebx,
        sel,
        param_name,
        rebx_param_value::string(val.to_string()),
    );
}

/// reboundx.h `rebx_set_param_pointer` for `REBX_TYPE_FORCE`. `val` is
/// an index into `rebx_extras::allocated_forces`.
pub fn rebx_set_param_force(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: usize) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::force(val));
}

/// reboundx.h `rebx_set_param_pointer` for `REBX_TYPE_ODE`. `val` is the
/// id of an ODE in `reb_simulation::odes`.
pub fn rebx_set_param_ode(rebx: &mut rebx_extras, sel: rebx_ap, param_name: &str, val: usize) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::ode(val));
}

/// reboundx.h `rebx_set_param_pointer` for the internally-allocated
/// `struct reb_particle` buffers (`im_ps_final`, `rk4_k2`, ...).
pub fn rebx_set_param_particles(
    rebx: &mut rebx_extras,
    sel: rebx_ap,
    param_name: &str,
    val: Vec<reb_particle>,
) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::particles(val));
}

/// reboundx.h `rebx_set_param_pointer` for the `"particle"`
/// back-reference (`track_min_distance`): the C stores a
/// `struct reb_particle*` into the simulation's own array, which here
/// is that particle's index.
pub fn rebx_set_param_particle_index(
    rebx: &mut rebx_extras,
    sel: rebx_ap,
    param_name: &str,
    val: usize,
) {
    rebx_set_param_value(rebx, sel, param_name, rebx_param_value::particle_index(val));
}

/*******************************************************************
 User interface for getting REBOUNDx objects and parameters
 *******************************************************************/

/// Name of the payload a parameter is actually holding, for the
/// wrong-type diagnostic below.
fn rebx_param_value_type_name(value: &rebx_param_value) -> &'static str {
    match value {
        rebx_param_value::none => "unset",
        rebx_param_value::double(_) => "double",
        rebx_param_value::int(_) => "int",
        rebx_param_value::uint32(_) => "uint32",
        rebx_param_value::vec3d(_) => "vec3d",
        rebx_param_value::orbit(_) => "orbit",
        rebx_param_value::string(_) => "string",
        rebx_param_value::force(_) => "force",
        rebx_param_value::ode(_) => "ode",
        rebx_param_value::particles(_) => "particle buffer",
        rebx_param_value::particle_index(_) => "particle index",
    }
}

/// Reported when a parameter exists but holds a different type than the
/// getter asked for. The C would blindly cast the `void*` and read
/// whatever bytes were there; this is a safety improvement, and for
/// correct code (a parameter always read with the type it was
/// registered and set with) it never fires.
///
/// Takes no `&mut rebx_extras` because the getters borrow the extras
/// immutably, so it writes to stderr rather than to `rebx.messages`.
fn rebx_param_type_mismatch(param_name: &str, requested: &str, value: &rebx_param_value) {
    eprintln!(
        "\nError! REBOUNDx Error: Parameter '{}' was requested as {} but holds a {}. \
         Returning no value.",
        param_name,
        requested,
        rebx_param_value_type_name(value)
    );
}

/// core.c `rebx_get_param_struct`. Walks the list named by `sel` in
/// order — index 0 is the C list head — and returns the first parameter
/// whose name matches.
///
/// Returns `None` where the C returns NULL. As in the C, no error is
/// raised: optional parameters are looked up this way all the time.
pub fn rebx_get_param_struct<'a>(
    rebx: &'a rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<&'a rebx_param> {
    rebx.ap(sel).iter().find(|param| param.name == param_name)
    // name not found. Don't want warnings for optional parameters so don't reb_simulation_error
}

/// reboundx.h `rebx_get_param_double`. `None` exactly where the C
/// returns NULL: the name is absent, or the parameter exists but was
/// never given a value (`param->value == NULL`).
pub fn rebx_get_param_double(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<f64> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::double(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "double", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param_int`.
pub fn rebx_get_param_int(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<i32> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::int(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "int", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param_uint32`.
pub fn rebx_get_param_uint32(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<u32> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::uint32(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "uint32", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param_vec3d`.
pub fn rebx_get_param_vec3d(
    rebx: &rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<reb_vec3d> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::vec3d(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "vec3d", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param_orbit`.
pub fn rebx_get_param_orbit(
    rebx: &rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<reb_orbit> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::orbit(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "orbit", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param_string`. Returns an owned copy; the C
/// returns the pointer REBOUND owns.
pub fn rebx_get_param_string(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<String> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::string(v) => Some(v.clone()),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "string", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param` for `REBX_TYPE_FORCE`. The returned
/// `usize` indexes `rebx_extras::allocated_forces`.
pub fn rebx_get_param_force(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<usize> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::force(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "force", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param` for `REBX_TYPE_ODE`.
pub fn rebx_get_param_ode(rebx: &rebx_extras, sel: rebx_ap, param_name: &str) -> Option<usize> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::ode(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "ode", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param` for the internally-allocated
/// `struct reb_particle` buffers.
pub fn rebx_get_param_particles<'a>(
    rebx: &'a rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<&'a Vec<reb_particle>> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::particles(v) => Some(v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "particle buffer", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_param` for the `"particle"` back-reference:
/// returns the index of the particle the C stored a pointer to.
pub fn rebx_get_param_particle_index(
    rebx: &rebx_extras,
    sel: rebx_ap,
    param_name: &str,
) -> Option<usize> {
    let param = rebx_get_param_struct(rebx, sel, param_name)?;
    match &param.value {
        rebx_param_value::particle_index(v) => Some(*v),
        rebx_param_value::none => None,
        other => {
            rebx_param_type_mismatch(param_name, "particle index", other);
            None
        }
    }
}

/// reboundx.h `rebx_get_force`. Returns the index into
/// `allocated_forces` of the most recently created force with this
/// name, or `None` where the C returns NULL.
///
/// The C walks `allocated_forces` from its head, which — because
/// `rebx_add_node` prepends — is the newest force first. This vector is
/// appended to instead (the index is the force's identity), so the same
/// traversal order is obtained by iterating in reverse.
pub fn rebx_get_force(rebx: &rebx_extras, name: &str) -> Option<usize> {
    // rposition: the first match walking from the newest force backwards,
    // i.e. from the head of the C's list. None where the C returns NULL.
    rebx.allocated_forces
        .iter()
        .rposition(|force| force.name == name)
}

/// reboundx.h `rebx_get_operator`. Reverse traversal for the same
/// reason as `rebx_get_force`.
pub fn rebx_get_operator(rebx: &rebx_extras, name: &str) -> Option<usize> {
    rebx.allocated_operators
        .iter()
        .rposition(|operator_| operator_.name == name)
}

/*******************************************************************
 User interface for removing REBOUNDx objects
 *******************************************************************/

/// reboundx.h `rebx_remove_force`. Returns 1 if the force was removed
/// from the list that actually affects the simulation
/// (`additional_forces`), as the C does.
///
/// The C also frees the force and unlinks it from `allocated_forces`.
/// Here a force's index *is* its identity, so unlinking it from the
/// middle of `allocated_forces` would renumber every force created
/// after it. The entry is therefore dropped only when it is the last
/// one — the case `rebx_load_force` relies on to undo a failed load —
/// and otherwise left in place. It holds no resource that needs
/// releasing before the extras are dropped.
pub fn rebx_remove_force(rebx: &mut rebx_extras, force_idx: usize) -> i32 {
    if force_idx + 1 == rebx.allocated_forces.len() {
        rebx.allocated_forces.pop();
    }

    // success only cares about removal from add_forces that affects sim
    let mut success = 0;
    if let Some(i) = rebx.additional_forces.iter().position(|f| *f == force_idx) {
        rebx.additional_forces.remove(i);
        success = 1;
    }
    success
}

/// core.c `rebx_remove_step_node`. Removes the first step in `steps`
/// that uses `operator_idx`. Returns 1 if one was removed.
fn rebx_remove_step_node(steps: &mut Vec<rebx_step>, operator_idx: usize) -> i32 {
    if let Some(i) = steps.iter().position(|s| s.operator_ == operator_idx) {
        steps.remove(i);
        return 1;
    }
    0
}

/// reboundx.h `rebx_remove_operator`. Removes every pre- and
/// post-timestep step that uses this operator. Returns 1 if at least
/// one was removed.
///
/// `allocated_operators` is treated exactly as `allocated_forces` is in
/// `rebx_remove_force` — see the note there.
pub fn rebx_remove_operator(rebx: &mut rebx_extras, operator_idx: usize) -> i32 {
    if operator_idx + 1 == rebx.allocated_operators.len() {
        rebx.allocated_operators.pop();
    }

    // success only cares about removal from lists that actually do
    // something to sim below. Success if EITHER one successful.
    let mut success = 0;
    let mut keep_searching = 1;
    while keep_searching == 1 {
        // keep searching while steps are found
        keep_searching = rebx_remove_step_node(&mut rebx.pre_timestep_modifications, operator_idx);
        if keep_searching == 1 {
            // success if at least one step found
            success = 1;
        }
    }

    keep_searching = 1;
    while keep_searching == 1 {
        // keep searching while steps are found
        keep_searching = rebx_remove_step_node(&mut rebx.post_timestep_modifications, operator_idx);
        if keep_searching == 1 {
            // success if at least one step found
            success = 1;
        }
    }

    success
}

/**********************************************
 Internal Functions executing forces & pre/post timestep modifications each timestep
 *********************************************/

/// core.c `rebx_reset_accelerations`. Used by the REBOUNDx integrators
/// (`integrate_force` and friends) before re-evaluating a force on a
/// scratch particle buffer.
pub fn rebx_reset_accelerations(ps: &mut [reb_particle], N: usize) {
    for i in 0..N {
        ps[i].ax = 0.;
        ps[i].ay = 0.;
        ps[i].az = 0.;
    }
}

/// core.c `rebx_additional_forces` — installed into
/// `sim.additional_forces` and called by REBOUND every time it
/// evaluates the gravity.
///
/// Walks `rebx.additional_forces` from index 0 (the C list head) and
/// calls each force's `update_accelerations`. That order is the reverse
/// of the order the forces were added, and it decides the order in
/// which the accelerations are summed — so it is load-bearing for
/// bit-for-bit agreement.
///
/// Note that REBOUNDx 5.1.0 does *not* zero the accelerations here:
/// REBOUND has already filled them with the gravitational term and each
/// effect adds to them. `rebx_reset_accelerations` exists for the
/// REBOUNDx integrators, which evaluate a force in isolation.
pub fn rebx_additional_forces(sim: &mut reb_simulation) {
    if sim.N_var != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Variational particles have been added to the simulation but are not implemented in REBOUNDx and will not be evolved self-consistently.");
    }
    // Take the extras out so that both `sim` and the REBOUNDx state can
    // be handed to the force functions. Every exit below restores it.
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => return,
    };

    let mut current = 0;
    while current < rebx.additional_forces.len() {
        if sim.force_is_velocity_dependent != 0 && sim.integrator.name() == "whfast" {
            reb_simulation_warning(sim, "REBOUNDx: Passing a velocity-dependent force to WHFAST, will accumulate errors proportional to the force. If forces get big, consider using IAS15 or applying force as an operator. See REBOUNDx paper sec 5.1 and ipython_examples/IntegrateForce.ipynb.");
        }
        let force = rebx.additional_forces[current];
        let N = sim.N;
        // Read the function pointer out in its own statement: it must
        // not still be borrowed from `rebx` when `rebx` is handed to it.
        let update_accelerations = rebx.allocated_forces[force].update_accelerations;
        if let Some(update_accelerations) = update_accelerations {
            update_accelerations(sim, &mut rebx, force, N);
        }
        current += 1;
    }

    rebx_put(sim, rebx);
}

/// core.c `rebx_pre_timestep_modifications` — installed into
/// `sim.pre_timestep_modifications`. Each step is applied with
/// `dt = sim.dt * step.dt_fraction`.
pub fn rebx_pre_timestep_modifications(sim: &mut reb_simulation) {
    if sim.N_var != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Variational particles have been added to the simulation but are not implemented in REBOUNDx and will not be evolved self-consistently.");
    }
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => return,
    };
    let dt = sim.dt;

    let mut current = 0;
    while current < rebx.pre_timestep_modifications.len() {
        let step = rebx.pre_timestep_modifications[current];
        let operator_ = step.operator_;
        let operator_type = rebx.allocated_operators[operator_].operator_type;
        // C: `struct reb_integrator_ias15_state* ias15 = sim->integrator.state;`
        // Copied out in its own statement so that `sim` is free to be
        // borrowed mutably by reb_simulation_warning below.
        let ias15_epsilon = match &sim.integrator {
            reb_integrator_state::ias15(ias15) => Some(ias15.epsilon),
            _ => None,
        };
        if sim.integrator.name() == "ias15" && operator_type == REBX_OPERATOR_UPDATER {
            if let Some(epsilon) = ias15_epsilon {
                if epsilon != 0. {
                    reb_simulation_warning(sim, "REBOUNDx: Operators that affect particle trajectories with adaptive timesteps can give spurious results. Use sim.ri_ias15.epsilon=0 for fixed timestep with IAS, or use a different integrator.");
                }
            }
        }
        let step_function = rebx.allocated_operators[operator_].step_function;
        if let Some(step_function) = step_function {
            step_function(sim, &mut rebx, operator_, dt * step.dt_fraction);
        }
        current += 1;
    }

    rebx_put(sim, rebx);
}

/// core.c `rebx_post_timestep_modifications` — installed into
/// `sim.post_timestep_modifications`. Note the C takes the timestep
/// from `sim->dt_last_done` here, not `sim->dt`: the step that just
/// finished may not have been the one that was requested.
pub fn rebx_post_timestep_modifications(sim: &mut reb_simulation) {
    if sim.N_var != 0 {
        reb_simulation_warning(sim, "REBOUNDx: Variational particles have been added to the simulation but are not implemented in REBOUNDx and will not be evolved self-consistently.");
    }
    let mut rebx = match rebx_take(sim) {
        Some(r) => r,
        None => return,
    };
    let dt = sim.dt_last_done;

    let mut current = 0;
    while current < rebx.post_timestep_modifications.len() {
        let step = rebx.post_timestep_modifications[current];
        let operator_ = step.operator_;
        let operator_type = rebx.allocated_operators[operator_].operator_type;
        let ias15_epsilon = match &sim.integrator {
            reb_integrator_state::ias15(ias15) => Some(ias15.epsilon),
            _ => None,
        };
        if sim.integrator.name() == "ias15" && operator_type == REBX_OPERATOR_UPDATER {
            if let Some(epsilon) = ias15_epsilon {
                if epsilon != 0. {
                    reb_simulation_warning(sim, "REBOUNDx: Operators that affect particle trajectories with adaptive timesteps can give spurious results. Use sim.ri_ias15.epsilon=0 for fixed timestep with IAS, or use a different integrator.");
                }
            }
        }
        let step_function = rebx.allocated_operators[operator_].step_function;
        if let Some(step_function) = step_function {
            step_function(sim, &mut rebx, operator_, dt * step.dt_fraction);
        }
        current += 1;
    }

    rebx_put(sim, rebx);
}

/****************************************************************
 Internal functions for dealing with parameters
 ****************************************************************/

/// core.c `rebx_create_param`. The C mallocs a node and copies the
/// name; here the struct owns its `String`. `value` starts unset,
/// which is the C's `param->value = NULL`.
///
/// core.c `rebx_create_node` has no counterpart: a `Vec` element is its
/// own node.
pub fn rebx_create_param(name: &str, type_: rebx_param_type) -> rebx_param {
    rebx_param {
        name: name.to_string(),
        type_,
        value: rebx_param_value::none,
    }
}

/// core.c `rebx_add_param`. Prepends, exactly as `rebx_add_node` does,
/// so the list traverses in reverse order of addition like the C's.
///
/// The C returns 0 when the node allocation fails; allocation cannot
/// fail here, so this always succeeds and returns nothing.
pub fn rebx_add_param(ap: &mut Vec<rebx_param>, param: rebx_param) {
    ap.insert(0, param);
}

/// core.c `rebx_get_type`. The type a parameter name was registered
/// with, or `REBX_TYPE_NONE` if it was never registered.
pub fn rebx_get_type(rebx: &rebx_extras, name: &str) -> rebx_param_type {
    for param in rebx.registered_params.iter() {
        if param.name == name {
            return param.type_;
        }
    }

    REBX_TYPE_NONE // param not found
}

/// core.c `rebx_sizeof`. Size in bytes of the payload for a parameter
/// type, used by the (unported) binary I/O.
pub fn rebx_sizeof(rebx: &mut rebx_extras, type_: rebx_param_type) -> usize {
    match type_ {
        REBX_TYPE_DOUBLE => std::mem::size_of::<f64>(),
        REBX_TYPE_INT => std::mem::size_of::<i32>(),
        REBX_TYPE_FORCE => std::mem::size_of::<rebx_force>(),
        REBX_TYPE_VEC3D => std::mem::size_of::<reb_vec3d>(),
        REBX_TYPE_POINTER => 0,
        REBX_TYPE_STRING => 0,
        REBX_TYPE_NONE => {
            rebx_error(rebx, "REBOUNDx Error: Parameter name passed to rebx_sizeof was not registered. This should not happen. Please open issue on github.com/dtamayo/reboundx.\n");
            0
        }
        _ => {
            rebx_error(rebx, "REBOUNDx Error: Need to add new param type to switch statement in rebx_sizeof. Please open issue on github.com/dtamayo/reboundx.\n");
            0
        }
    }
}

/// core.c `rebx_error`.
///
/// The C routes the message to `reb_simulation_error(rebx->sim, msg)`,
/// which prints it to stderr (or stores it when `sim->save_messages` is
/// set). `rebx_extras` holds no back-pointer to the simulation here, so
/// the message goes to stderr in REBOUND's format and is also kept in
/// `rebx.messages` for callers that want to inspect it.
pub fn rebx_error(rebx: &mut rebx_extras, msg: &str) {
    eprintln!("\nError! {}", msg);
    rebx.messages.push(msg.to_string());
}
