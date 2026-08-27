//! Integration tests for the core_params group of reboundx_rs.
//! Part of reboundx_rs, GPL-3.0-or-later.
#![allow(non_snake_case)]
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::redundant_closure)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::too_many_arguments)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. Each waiver below carries its own reason; the same
// list and the rationale are in README.md under "Building and testing".
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::assign_op_pattern)]
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
#![allow(clippy::ptr_arg)]
#![allow(clippy::seek_from_current)]
#![allow(clippy::drop_non_drop)]
#![allow(clippy::approx_constant)]
#![allow(clippy::useless_vec)]
#![allow(clippy::type_complexity)]
use rebound_rs::*;
use reboundx_rs::*;

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Every state variable that a timestep can touch, as raw IEEE-754 bits.
/// Used for the bit-for-bit comparisons the port's prime directive asks
/// for (f64 equality would silently accept +0.0 == -0.0).
fn state_bits(sim: &reb_simulation) -> Vec<u64> {
    let mut v = Vec::with_capacity(1 + 7 * sim.N);
    v.push(sim.t.to_bits());
    for i in 0..sim.N {
        let p = sim.particles[i];
        v.push(p.x.to_bits());
        v.push(p.y.to_bits());
        v.push(p.z.to_bits());
        v.push(p.vx.to_bits());
        v.push(p.vy.to_bits());
        v.push(p.vz.to_bits());
        v.push(p.m.to_bits());
    }
    v
}

/// A star at the origin plus one bound, eccentric massless companion.
///
/// G = 1 and the companion is exactly massless, so the star feels no
/// acceleration at all and stays pinned at the origin; that removes the
/// centre-of-mass drift from every comparison below. With `star_m = 1`
/// the companion has v^2 = 1.21 at r = 1, i.e. vis-viva 1/a = 2 - 1.21 =
/// 0.79, a = 1.2658..., P = 2*pi*a^(3/2) = 8.95 — comfortably bound.
fn two_body(star_m: f64) -> reb_simulation {
    let mut sim = reb_simulation_create();
    sim.G = 1.0;
    sim.testparticle_hidewarnings = 1;
    reb_simulation_add(
        &mut sim,
        reb_particle {
            m: star_m,
            ..Default::default()
        },
    );
    reb_simulation_add(
        &mut sim,
        reb_particle {
            m: 0.0,
            x: 1.0,
            vy: 1.1,
            ..Default::default()
        },
    );
    sim
}

/// A stand-in for a user's own `additional_forces` callback.
fn foreign_additional_forces(_sim: &mut reb_simulation) {}

/// The 107 default parameter names, in the order `rebx_register_default_params`
/// registers them. Transcribed from REBOUNDx 5.1.0 core.c; used to check
/// both the contents and the (prepending) list order.
const DEFAULT_PARAM_NAMES: [&str; 107] = [
    "c",
    "gr_source",
    "tau_mass",
    "force",
    "particle",
    "Acentral",
    "gammacentral",
    "max_iterations",
    "J2",
    "J4",
    "R_eq",
    "coordinates",
    "p",
    "d_factor",
    "cs_coeff",
    "tau_coeff",
    "tau_a",
    "tau_e",
    "tau_inc",
    "tau_omega",
    "tau_Omega",
    "em_tau_a",
    "em_aini",
    "em_afin",
    "primary",
    "radiation_source",
    "kappa",
    "kappa_x",
    "kappa_y",
    "kappa_z",
    "tau_kappa",
    "tau_kappa_x",
    "tau_kappa_y",
    "tau_kappa_z",
    "stochastic_force_r",
    "stochastic_force_phi",
    "stochastic_force_x",
    "stochastic_force_y",
    "stochastic_force_z",
    "beta",
    "tides_primary",
    "R_tides",
    "tctl_k2",
    "tctl_tau",
    "integrator",
    "im_ps_final",
    "im_ps_prev",
    "im_ps_avg",
    "rk2_k2",
    "rk4_k2",
    "rk4_k3",
    "min_distance",
    "min_distance_from",
    "min_distance_orbit",
    "luminosity",
    "ide_position",
    "ide_width",
    "tIm_flaring_index",
    "tIm_scale_height_1",
    "tIm_surface_density_1",
    "tIm_surface_density_exponent",
    "ye_c",
    "ye_body_density",
    "ye_lstar",
    "ye_flag",
    "ye_rotation_period",
    "ye_thermal_inertia",
    "ye_albedo",
    "ye_emissivity",
    "ye_k",
    "ye_stef_boltz",
    "ye_spin_axis_x",
    "ye_spin_axis_y",
    "ye_spin_axis_z",
    "OmegaMag",
    "Omega",
    "k2",
    "I",
    "tau",
    "ode",
    "gas_df_rhog",
    "gas_df_alpha_rhog",
    "gas_df_cs",
    "gas_df_alpha_cs",
    "gas_df_xmin",
    "gas_df_hr",
    "gas_df_Qd",
    "lt_R_eq",
    "lt_Mom_I_fac",
    "lt_rot_rate",
    "lt_p_hatx",
    "lt_p_haty",
    "lt_p_hatz",
    "lt_c",
    "td_M_last",
    "td_num_apoapsis",
    "td_c_imag",
    "td_c_real",
    "td_dP_hat",
    "td_dP_crit",
    "td_EB0",
    "td_E_max",
    "td_E_resid",
    "td_dE_last",
    "td_last_apoapsis",
    "td_drag_coef",
    "td_disruption_flag",
];

// ---------------------------------------------------------------------
// rebx_attach: default parameter registration
// ---------------------------------------------------------------------

#[test]
fn attach_registers_the_107_default_params_in_C_order() {
    let mut sim = reb_simulation_create();
    assert!(
        rebx_extras_ref(&sim).is_none(),
        "extras must be absent before rebx_attach, got Some"
    );

    rebx_attach(&mut sim);
    let rebx = rebx_extras_ref(&sim).expect("rebx_attach must install the extras box");

    assert_eq!(
        rebx.registered_params.len(),
        DEFAULT_PARAM_NAMES.len(),
        "registered_params.len() = {}, expected {} (core.c registers 107 names)",
        rebx.registered_params.len(),
        DEFAULT_PARAM_NAMES.len()
    );

    // rebx_add_param PREPENDS (mirroring the C's rebx_add_node), so the
    // list reads back in reverse registration order: index 0 is the last
    // name registered.
    for (k, expected) in DEFAULT_PARAM_NAMES.iter().rev().enumerate() {
        assert_eq!(
            &rebx.registered_params[k].name, expected,
            "registered_params[{}].name = '{}', expected '{}' (list is prepended, \
             so it must read back in reverse registration order)",
            k, rebx.registered_params[k].name, expected
        );
    }

    // Every registered parameter starts with no value: the C's
    // `param->value = NULL`.
    for p in rebx.registered_params.iter() {
        assert_eq!(
            p.value,
            rebx_param_value::none,
            "registered param '{}' should carry no value, got {:?}",
            p.name,
            p.value
        );
    }
}

#[test]
fn registered_param_types_match_the_C_table() {
    let mut sim = reb_simulation_create();
    rebx_attach(&mut sim);
    let rebx = rebx_extras_ref(&sim).expect("extras attached");

    // One representative of each rebx_param_type used by core.c.
    let expected: [(&str, rebx_param_type); 8] = [
        ("c", rebx_param_type::REBX_TYPE_DOUBLE),
        ("gr_source", rebx_param_type::REBX_TYPE_INT),
        ("force", rebx_param_type::REBX_TYPE_FORCE),
        ("particle", rebx_param_type::REBX_TYPE_POINTER),
        ("min_distance_from", rebx_param_type::REBX_TYPE_STRING),
        ("min_distance_orbit", rebx_param_type::REBX_TYPE_ORBIT),
        ("Omega", rebx_param_type::REBX_TYPE_VEC3D),
        ("ode", rebx_param_type::REBX_TYPE_ODE),
    ];
    for (name, want) in expected.iter() {
        let got = rebx_get_type(rebx, name);
        assert_eq!(
            got, *want,
            "rebx_get_type('{}') = {:?}, expected {:?}",
            name, got, want
        );
    }

    let got = rebx_get_type(rebx, "definitely_not_a_registered_parameter");
    assert_eq!(
        got,
        rebx_param_type::REBX_TYPE_NONE,
        "rebx_get_type on an unregistered name must be REBX_TYPE_NONE, got {:?}",
        got
    );
}

#[test]
fn register_param_rejects_a_duplicate_name() {
    let mut sim = reb_simulation_create();
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    let before = rebx.registered_params.len();
    rebx_register_param(rebx, "c", rebx_param_type::REBX_TYPE_DOUBLE);
    assert_eq!(
        rebx.registered_params.len(),
        before,
        "registering the duplicate name 'c' must not grow registered_params: \
         len {} -> {}",
        before,
        rebx.registered_params.len()
    );
    assert!(
        rebx.messages.iter().any(|m| m.contains("already in registered list")),
        "expected an 'already in registered list' error for the duplicate, messages = {:?}",
        rebx.messages
    );

    // A genuinely new name is accepted and gets the type it was given.
    rebx_register_param(rebx, "core_params_test_u32", rebx_param_type::REBX_TYPE_UINT32);
    assert_eq!(
        rebx.registered_params.len(),
        before + 1,
        "registering a fresh name must grow registered_params: len {} -> {}",
        before,
        rebx.registered_params.len()
    );
    let t = rebx_get_type(rebx, "core_params_test_u32");
    assert_eq!(
        t,
        rebx_param_type::REBX_TYPE_UINT32,
        "freshly registered 'core_params_test_u32' has type {:?}, expected REBX_TYPE_UINT32",
        t
    );
}

// ---------------------------------------------------------------------
// rebx_attach / rebx_detach: the REBOUND callback hooks
// ---------------------------------------------------------------------

#[test]
fn attach_installs_no_hooks_until_effects_are_added() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);

    assert!(
        sim.additional_forces.is_none(),
        "rebx_attach alone must not set sim.additional_forces"
    );
    assert!(
        sim.pre_timestep_modifications.is_none(),
        "rebx_attach alone must not set sim.pre_timestep_modifications"
    );
    assert!(
        sim.post_timestep_modifications.is_none(),
        "rebx_attach alone must not set sim.post_timestep_modifications"
    );

    let f = rebx_load_force(&mut sim, "central_force").expect("central_force is in the library");
    let ok = rebx_add_force(&mut sim, f);
    assert_eq!(ok, 1, "rebx_add_force returned {}, expected 1 (success)", ok);
    assert!(
        sim.additional_forces.is_some(),
        "rebx_add_force must install sim.additional_forces"
    );

    reb_simulation_set_integrator(&mut sim, "whfast");
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("modify_mass is in the library");
    let ok = rebx_add_operator(&mut sim, op);
    assert_eq!(ok, 1, "rebx_add_operator returned {}, expected 1", ok);
    assert!(
        sim.pre_timestep_modifications.is_some() && sim.post_timestep_modifications.is_some(),
        "an updater operator under WHFast must install BOTH the pre- and \
         post-timestep hooks (pre set: {}, post set: {})",
        sim.pre_timestep_modifications.is_some(),
        sim.post_timestep_modifications.is_some()
    );

    rebx_detach(&mut sim);
    assert!(
        sim.extras.is_none(),
        "rebx_detach must drop the extras box"
    );
    assert!(
        sim.additional_forces.is_none()
            && sim.pre_timestep_modifications.is_none()
            && sim.post_timestep_modifications.is_none(),
        "rebx_detach must clear all three REBOUNDx hooks (af {}, pre {}, post {})",
        sim.additional_forces.is_some(),
        sim.pre_timestep_modifications.is_some(),
        sim.post_timestep_modifications.is_some()
    );
}

#[test]
fn detach_leaves_a_foreign_callback_alone() {
    // core.c's rebx_detach only clears a hook that still points at
    // REBOUNDx's own function.
    let mut sim = two_body(1.0);
    let foreign: fn(&mut reb_simulation) = foreign_additional_forces;
    sim.additional_forces = Some(foreign);

    rebx_attach(&mut sim);
    rebx_detach(&mut sim);

    match sim.additional_forces {
        None => panic!("rebx_detach wrongly cleared a foreign additional_forces callback"),
        Some(f) => assert!(
            std::ptr::fn_addr_eq(f, foreign),
            "sim.additional_forces no longer points at the caller's own function"
        ),
    }
}

#[test]
fn rebx_with_fails_cleanly_when_nothing_is_attached_and_restores_otherwise() {
    let mut sim = two_body(1.0);
    let r = rebx_with(&mut sim, |_sim, _rebx| 42i32);
    assert!(
        r.is_none(),
        "rebx_with on a simulation with no extras must return None, got {:?}",
        r
    );

    rebx_attach(&mut sim);
    let r = rebx_with(&mut sim, |sim, rebx| {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "c", 1234.5);
        sim.N
    });
    assert_eq!(r, Some(2), "rebx_with must return the closure's value, got {:?}", r);
    // The take/put pair must have put the box back.
    let v = rebx_extras_ref(&sim)
        .and_then(|rebx| rebx_get_param_double(rebx, rebx_ap::particle(0), "c"));
    assert_eq!(
        v,
        Some(1234.5),
        "after rebx_with the extras box must still hold the parameter written inside it, got {:?}",
        v
    );
}

// ---------------------------------------------------------------------
// set / get round-trips, one per payload type
// ---------------------------------------------------------------------

#[test]
fn double_param_round_trips_bit_for_bit() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    // Values chosen to catch anything that goes through a lossy path:
    // a negative, a subnormal, the largest finite double, and a value
    // whose decimal form is not exactly representable.
    let vals = [
        -3.25e-7_f64,
        f64::MIN_POSITIVE / 4.0,
        f64::MAX,
        0.1 + 0.2,
        -0.0,
    ];
    for (k, v) in vals.iter().enumerate() {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "tau_a", *v);
        let got = rebx_get_param_double(rebx, rebx_ap::particle(0), "tau_a")
            .unwrap_or_else(|| panic!("tau_a missing after set (case {})", k));
        assert_eq!(
            got.to_bits(),
            v.to_bits(),
            "double round-trip case {}: stored {:e} (bits {:016x}), read back {:e} (bits {:016x})",
            k,
            v,
            v.to_bits(),
            got,
            got.to_bits()
        );
    }

    // Re-setting overwrites in place; it must not create a second entry.
    let n = rebx.ap(rebx_ap::particle(0)).len();
    assert_eq!(
        n, 1,
        "after {} sets of the same name the list holds {} params, expected 1",
        vals.len(),
        n
    );
}

#[test]
fn int_and_uint32_params_round_trip() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    for v in [i32::MIN, -1, 0, 1, i32::MAX] {
        rebx_set_param_int(rebx, rebx_ap::particle(1), "gr_source", v);
        let got = rebx_get_param_int(rebx, rebx_ap::particle(1), "gr_source");
        assert_eq!(
            got,
            Some(v),
            "int round-trip: stored {}, read back {:?}",
            v,
            got
        );
    }

    // No default parameter is registered as UINT32, so register one —
    // this is exactly the custom-effect path in the C.
    rebx_register_param(rebx, "cp_u32", rebx_param_type::REBX_TYPE_UINT32);
    for v in [0u32, 1, 0x8000_0000, u32::MAX] {
        rebx_set_param_uint32(rebx, rebx_ap::particle(1), "cp_u32", v);
        let got = rebx_get_param_uint32(rebx, rebx_ap::particle(1), "cp_u32");
        assert_eq!(
            got,
            Some(v),
            "uint32 round-trip: stored {}, read back {:?}",
            v,
            got
        );
    }
}

#[test]
fn vec3d_param_round_trips_component_wise() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    let v = reb_vec3d {
        x: -1.0 / 3.0,
        y: 7.5e12,
        z: f64::MIN_POSITIVE,
    };
    rebx_set_param_vec3d(rebx, rebx_ap::particle(0), "Omega", v);
    let got = rebx_get_param_vec3d(rebx, rebx_ap::particle(0), "Omega").expect("Omega set above");
    assert_eq!(
        (got.x.to_bits(), got.y.to_bits(), got.z.to_bits()),
        (v.x.to_bits(), v.y.to_bits(), v.z.to_bits()),
        "vec3d round-trip: stored ({:e},{:e},{:e}), read back ({:e},{:e},{:e})",
        v.x,
        v.y,
        v.z,
        got.x,
        got.y,
        got.z
    );
}

#[test]
fn orbit_param_round_trips_field_for_field() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);

    // A real orbit rather than a made-up struct, so every field is a
    // non-trivial double.
    let orb = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    assert!(
        orb.a > 0.0 && orb.e < 1.0,
        "the fixture orbit must be bound: a = {:e}, e = {:e}",
        orb.a,
        orb.e
    );

    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");
    rebx_set_param_orbit(rebx, rebx_ap::particle(1), "min_distance_orbit", orb);
    let got = rebx_get_param_orbit(rebx, rebx_ap::particle(1), "min_distance_orbit")
        .expect("min_distance_orbit set above");

    let fields: [(&str, f64, f64); 12] = [
        ("d", orb.d, got.d),
        ("v", orb.v, got.v),
        ("h", orb.h, got.h),
        ("P", orb.P, got.P),
        ("n", orb.n, got.n),
        ("a", orb.a, got.a),
        ("e", orb.e, got.e),
        ("inc", orb.inc, got.inc),
        ("omega", orb.omega, got.omega),
        ("f", orb.f, got.f),
        ("M", orb.M, got.M),
        ("T", orb.T, got.T),
    ];
    for (name, want, have) in fields.iter() {
        assert_eq!(
            have.to_bits(),
            want.to_bits(),
            "orbit round-trip field '{}': stored {:e}, read back {:e}",
            name,
            want,
            have
        );
    }
    assert_eq!(
        got, orb,
        "the whole reb_orbit must round-trip unchanged: stored {:?}, read back {:?}",
        orb, got
    );
}

#[test]
fn string_force_ode_and_pointer_params_round_trip() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let gr = rebx_load_force(&mut sim, "gr").expect("gr is in the library");
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    rebx_set_param_string(rebx, rebx_ap::particle(1), "min_distance_from", "primary");
    let s = rebx_get_param_string(rebx, rebx_ap::particle(1), "min_distance_from");
    assert_eq!(
        s.as_deref(),
        Some("primary"),
        "string round-trip: stored \"primary\", read back {:?}",
        s
    );
    // Overwriting replaces the text.
    rebx_set_param_string(rebx, rebx_ap::particle(1), "min_distance_from", "com");
    let s = rebx_get_param_string(rebx, rebx_ap::particle(1), "min_distance_from");
    assert_eq!(
        s.as_deref(),
        Some("com"),
        "string overwrite: expected \"com\", read back {:?}",
        s
    );

    rebx_set_param_force(rebx, rebx_ap::particle(0), "force", gr);
    let f = rebx_get_param_force(rebx, rebx_ap::particle(0), "force");
    assert_eq!(
        f,
        Some(gr),
        "force round-trip: stored force index {}, read back {:?}",
        gr,
        f
    );

    rebx_set_param_ode(rebx, rebx_ap::particle(0), "ode", 3);
    let o = rebx_get_param_ode(rebx, rebx_ap::particle(0), "ode");
    assert_eq!(o, Some(3), "ode round-trip: stored 3, read back {:?}", o);

    // The "particle" back-reference (track_min_distance) is a particle
    // index here, where the C stores a reb_particle*.
    rebx_set_param_particle_index(rebx, rebx_ap::particle(0), "particle", 1);
    let pi = rebx_get_param_particle_index(rebx, rebx_ap::particle(0), "particle");
    assert_eq!(
        pi,
        Some(1),
        "particle-index round-trip: stored 1, read back {:?}",
        pi
    );

    // An internally-allocated reb_particle buffer (the REBOUNDx
    // integrators' scratch space).
    let buf = vec![
        reb_particle { m: 2.0, x: -4.0, ..Default::default() },
        reb_particle { m: 5.0, vz: 0.25, ..Default::default() },
    ];
    rebx_set_param_particles(rebx, rebx_ap::particle(0), "im_ps_final", buf.clone());
    let got = rebx_get_param_particles(rebx, rebx_ap::particle(0), "im_ps_final")
        .expect("im_ps_final set above");
    assert_eq!(
        got.len(),
        buf.len(),
        "particle-buffer round-trip: stored {} particles, read back {}",
        buf.len(),
        got.len()
    );
    for (k, (a, b)) in buf.iter().zip(got.iter()).enumerate() {
        assert_eq!(
            (a.m.to_bits(), a.x.to_bits(), a.vz.to_bits()),
            (b.m.to_bits(), b.x.to_bits(), b.vz.to_bits()),
            "particle-buffer element {}: stored (m {:e}, x {:e}, vz {:e}), \
             read back (m {:e}, x {:e}, vz {:e})",
            k,
            a.m,
            a.x,
            a.vz,
            b.m,
            b.x,
            b.vz
        );
    }
}

// ---------------------------------------------------------------------
// Where the C returns NULL, the getters return None
// ---------------------------------------------------------------------

#[test]
fn getters_return_none_for_a_name_that_was_never_set() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    // Registered, but never set on this particle: the C's
    // rebx_get_param_struct returns NULL.
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "tau_a"),
        None,
        "tau_a was never set on particle 0, so the getter must return None"
    );
    assert_eq!(
        rebx_get_param_int(rebx, rebx_ap::particle(0), "gr_source"),
        None,
        "gr_source was never set on particle 0, so the getter must return None"
    );
    assert_eq!(
        rebx_get_param_vec3d(rebx, rebx_ap::particle(0), "Omega"),
        None,
        "Omega was never set on particle 0, so the getter must return None"
    );
    assert!(
        rebx_get_param_struct(rebx, rebx_ap::particle(0), "tau_a").is_none(),
        "rebx_get_param_struct must also report the parameter as absent"
    );

    // Setting it on particle 0 must not make it appear on particle 1, on
    // a force, or on an operator.
    rebx_set_param_double(rebx, rebx_ap::particle(0), "tau_a", 9.0);
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(1), "tau_a"),
        None,
        "tau_a set on particle 0 must not be visible on particle 1"
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "tau_a"),
        Some(9.0),
        "tau_a on particle 0 should read back 9.0"
    );

    // A particle index far past any particle the setters ever touched has
    // an empty list, not a stale one.
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(99), "tau_a"),
        None,
        "an untouched particle index must have no parameters"
    );
}

#[test]
fn getter_returns_none_for_a_registered_but_unvalued_param() {
    // rebx_get_or_add_param creates the node with value == NULL in the C;
    // rebx_get_param then returns that NULL.
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    let idx = rebx_get_or_add_param(rebx, rebx_ap::particle(0), "min_distance")
        .expect("min_distance is a registered name, so the node must be created");
    assert_eq!(idx, 0, "the newly prepended param must sit at the list head, got {}", idx);

    // The node exists ...
    let p = rebx_get_param_struct(rebx, rebx_ap::particle(0), "min_distance")
        .expect("the node was just created");
    assert_eq!(
        p.value,
        rebx_param_value::none,
        "a freshly created param must hold no value, got {:?}",
        p.value
    );
    assert_eq!(
        p.type_,
        rebx_param_type::REBX_TYPE_DOUBLE,
        "the node must inherit the registered type, got {:?}",
        p.type_
    );

    // ... but the typed getter still reports NULL.
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "min_distance"),
        None,
        "a registered-but-unvalued param must read back as None (C: value == NULL)"
    );

    // Calling get_or_add again finds the existing node instead of adding
    // a second one.
    let again = rebx_get_or_add_param(rebx, rebx_ap::particle(0), "min_distance");
    assert_eq!(
        again,
        Some(0),
        "the second get_or_add must find the existing node at index 0, got {:?}",
        again
    );
    assert_eq!(
        rebx.ap(rebx_ap::particle(0)).len(),
        1,
        "two get_or_add calls must leave exactly one node, got {}",
        rebx.ap(rebx_ap::particle(0)).len()
    );
}

#[test]
fn getter_returns_none_when_the_stored_type_differs() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    rebx_set_param_double(rebx, rebx_ap::particle(0), "c", 1.0e4);
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "c"),
        Some(1.0e4),
        "the matching getter must return the stored value"
    );
    assert_eq!(
        rebx_get_param_int(rebx, rebx_ap::particle(0), "c"),
        None,
        "reading a double-valued param as int must yield None, not reinterpreted bytes"
    );
    assert_eq!(
        rebx_get_param_vec3d(rebx, rebx_ap::particle(0), "c"),
        None,
        "reading a double-valued param as vec3d must yield None"
    );
    assert_eq!(
        rebx_get_param_string(rebx, rebx_ap::particle(0), "c"),
        None,
        "reading a double-valued param as string must yield None"
    );
    // The value itself is untouched by the failed reads.
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "c"),
        Some(1.0e4),
        "a wrong-type read must not disturb the stored value"
    );
}

#[test]
fn setting_an_unregistered_name_is_refused() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    rebx_set_param_double(rebx, rebx_ap::particle(0), "not_a_registered_name", 1.0);
    assert_eq!(
        rebx.ap(rebx_ap::particle(0)).len(),
        0,
        "an unregistered name must not create a parameter, list len = {}",
        rebx.ap(rebx_ap::particle(0)).len()
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "not_a_registered_name"),
        None,
        "an unregistered name must read back as None"
    );
    assert!(
        rebx.messages
            .iter()
            .any(|m| m.contains("Need to register parameter name")),
        "expected a 'Need to register parameter name' error, messages = {:?}",
        rebx.messages
    );
}

#[test]
fn a_selector_naming_no_object_is_refused() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    // No forces or operators have been created, so index 0 names nothing;
    // this is the C's `apptr == NULL`.
    let r = rebx_get_or_add_param(rebx, rebx_ap::force(0), "tau_a");
    assert!(
        r.is_none(),
        "get_or_add on a force index that names nothing must return None, got {:?}",
        r
    );
    rebx_set_param_double(rebx, rebx_ap::operator_(0), "tau_a", 1.0);
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::operator_(0), "tau_a"),
        None,
        "setting on an operator index that names nothing must store nothing"
    );
    assert!(
        rebx.messages.iter().any(|m| m.contains("NULL apptr")),
        "expected a 'NULL apptr' error, messages = {:?}",
        rebx.messages
    );
}

// ---------------------------------------------------------------------
// Parameters attach independently to particles, forces and operators
// ---------------------------------------------------------------------

#[test]
fn params_on_particles_forces_and_operators_are_independent() {
    let mut sim = two_body(1.0);
    reb_simulation_set_integrator(&mut sim, "whfast");
    rebx_attach(&mut sim);
    let f = rebx_load_force(&mut sim, "modify_orbits_forces").expect("force in library");
    let op = rebx_load_operator(&mut sim, "modify_orbits_direct").expect("operator in library");
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    // Four different owners, same parameter name, four different values.
    let owners = [
        (rebx_ap::particle(0), 10.0_f64),
        (rebx_ap::particle(1), 20.0),
        (rebx_ap::force(f), 30.0),
        (rebx_ap::operator_(op), 40.0),
    ];
    for (sel, v) in owners.iter() {
        rebx_set_param_double(rebx, *sel, "tau_a", *v);
    }
    for (sel, v) in owners.iter() {
        let got = rebx_get_param_double(rebx, *sel, "tau_a");
        assert_eq!(
            got,
            Some(*v),
            "tau_a on {:?}: expected {}, got {:?}",
            sel,
            v,
            got
        );
    }

    // Overwriting one must leave the other three alone.
    rebx_set_param_double(rebx, rebx_ap::force(f), "tau_a", -1.0);
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::force(f), "tau_a"),
        Some(-1.0),
        "the force's tau_a should now be -1.0"
    );
    for (sel, v) in [
        (rebx_ap::particle(0), 10.0_f64),
        (rebx_ap::particle(1), 20.0),
        (rebx_ap::operator_(op), 40.0),
    ] {
        let got = rebx_get_param_double(rebx, sel, "tau_a");
        assert_eq!(
            got,
            Some(v),
            "writing the force's tau_a disturbed {:?}: expected {}, got {:?}",
            sel,
            v,
            got
        );
    }

    // rebx_free_ap empties exactly one list.
    rebx_free_ap(rebx, rebx_ap::particle(0));
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "tau_a"),
        None,
        "rebx_free_ap must empty particle 0's list"
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(1), "tau_a"),
        Some(20.0),
        "rebx_free_ap on particle 0 must not touch particle 1"
    );
    rebx_free_particle_ap(rebx, 1);
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(1), "tau_a"),
        None,
        "rebx_free_particle_ap must empty particle 1's list"
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::operator_(op), "tau_a"),
        Some(40.0),
        "freeing particle lists must not touch the operator's list"
    );
}

#[test]
fn a_param_list_reads_back_head_first() {
    // core.c's rebx_add_param prepends, so the most recently added
    // parameter is at the head of the list, and re-setting an existing
    // one does not move it.
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    for (name, v) in [("tau_a", 1.0_f64), ("tau_e", 2.0), ("tau_inc", 3.0)] {
        rebx_set_param_double(rebx, rebx_ap::particle(0), name, v);
    }
    let names: Vec<String> = rebx
        .ap(rebx_ap::particle(0))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    assert_eq!(
        names,
        vec!["tau_inc", "tau_e", "tau_a"],
        "the list must read back in reverse insertion order, got {:?}",
        names
    );

    rebx_set_param_double(rebx, rebx_ap::particle(0), "tau_a", 99.0);
    let names2: Vec<String> = rebx
        .ap(rebx_ap::particle(0))
        .iter()
        .map(|p| p.name.clone())
        .collect();
    assert_eq!(
        names2, names,
        "re-setting an existing parameter must not reorder the list: {:?} -> {:?}",
        names, names2
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "tau_a"),
        Some(99.0),
        "the re-set value should be 99.0"
    );
}

// ---------------------------------------------------------------------
// Loading effects by name
// ---------------------------------------------------------------------

#[test]
fn loading_an_unknown_force_fails_without_leaking_an_entry() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);

    let bad = rebx_load_force(&mut sim, "no_such_force_exists");
    assert!(
        bad.is_none(),
        "rebx_load_force on an unknown name must return None, got {:?}",
        bad
    );
    {
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        assert_eq!(
            rebx.allocated_forces.len(),
            0,
            "the failed load must undo its own rebx_create_force; \
             allocated_forces.len() = {}",
            rebx.allocated_forces.len()
        );
        assert!(
            rebx.messages
                .iter()
                .any(|m| m.contains("not found in REBOUNDx library")),
            "expected a 'not found in REBOUNDx library' error, messages = {:?}",
            rebx.messages
        );
    }

    // The next good load must therefore get index 0.
    let good = rebx_load_force(&mut sim, "gr_potential").expect("gr_potential is in the library");
    assert_eq!(
        good, 0,
        "after a failed load the next force must take index 0, got {}",
        good
    );
    let rebx = rebx_extras_ref(&sim).expect("extras attached");
    assert_eq!(
        rebx.allocated_forces[good].force_type,
        rebx_force_type::REBX_FORCE_POS,
        "gr_potential is derivable from a potential, so force_type must be REBX_FORCE_POS, got {:?}",
        rebx.allocated_forces[good].force_type
    );
    assert!(
        rebx.allocated_forces[good].update_accelerations.is_some(),
        "a loaded force must carry an update_accelerations function"
    );
}

#[test]
fn loading_an_unknown_operator_fails_without_leaking_an_entry() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);

    let bad = rebx_load_operator(&mut sim, "no_such_operator_exists");
    assert!(
        bad.is_none(),
        "rebx_load_operator on an unknown name must return None, got {:?}",
        bad
    );
    {
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        assert_eq!(
            rebx.allocated_operators.len(),
            0,
            "the failed load must undo its own rebx_create_operator; \
             allocated_operators.len() = {}",
            rebx.allocated_operators.len()
        );
    }

    let good = rebx_load_operator(&mut sim, "track_min_distance").expect("operator in library");
    let rebx = rebx_extras_ref(&sim).expect("extras attached");
    assert_eq!(
        rebx.allocated_operators[good].operator_type,
        rebx_operator_type::REBX_OPERATOR_RECORDER,
        "track_min_distance only records, so operator_type must be RECORDER, got {:?}",
        rebx.allocated_operators[good].operator_type
    );
    assert!(
        rebx.allocated_operators[good].step_function.is_some(),
        "a loaded operator must carry a step_function"
    );
}

#[test]
fn get_force_and_get_operator_find_things_by_name() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let gr = rebx_load_force(&mut sim, "gr").expect("gr in library");
    let cf = rebx_load_force(&mut sim, "central_force").expect("central_force in library");
    let mm = rebx_load_operator(&mut sim, "modify_mass").expect("modify_mass in library");

    let rebx = rebx_extras_ref(&sim).expect("extras attached");
    assert_eq!(
        rebx_get_force(rebx, "gr"),
        Some(gr),
        "rebx_get_force('gr') should return the index rebx_load_force handed out ({})",
        gr
    );
    assert_eq!(
        rebx_get_force(rebx, "central_force"),
        Some(cf),
        "rebx_get_force('central_force') should return {}",
        cf
    );
    assert_eq!(
        rebx_get_force(rebx, "modify_mass"),
        None,
        "modify_mass is an operator, so rebx_get_force must not find it"
    );
    assert_eq!(
        rebx_get_operator(rebx, "modify_mass"),
        Some(mm),
        "rebx_get_operator('modify_mass') should return {}",
        mm
    );
    assert_eq!(
        rebx_get_operator(rebx, "gr"),
        None,
        "gr is a force, so rebx_get_operator must not find it"
    );
    assert_eq!(
        rebx_get_force(rebx, "nothing_by_this_name"),
        None,
        "rebx_get_force on an absent name must be None"
    );
}

#[test]
fn get_force_returns_the_most_recently_created_namesake() {
    // The C's allocated_forces list is prepended to, so rebx_get_force
    // walks newest-first. Two forces with the same name must resolve to
    // the newer one.
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");

    let first = rebx_create_force(rebx, "twin");
    let second = rebx_create_force(rebx, "twin");
    assert!(
        second > first,
        "the second creation must get a later index: first {}, second {}",
        first,
        second
    );
    assert_eq!(
        rebx_get_force(rebx, "twin"),
        Some(second),
        "rebx_get_force must resolve to the newest namesake ({}), not the oldest ({})",
        second,
        first
    );

    let a = rebx_create_operator(rebx, "twin_op");
    let b = rebx_create_operator(rebx, "twin_op");
    assert_eq!(
        rebx_get_operator(rebx, "twin_op"),
        Some(b),
        "rebx_get_operator must resolve to the newest namesake ({}), not the oldest ({})",
        b,
        a
    );
}

// ---------------------------------------------------------------------
// Adding and removing effects
// ---------------------------------------------------------------------

#[test]
fn add_force_refuses_an_incomplete_or_nonexistent_force() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);

    // C: force == NULL.
    let ok = rebx_add_force(&mut sim, 7);
    assert_eq!(
        ok, 0,
        "rebx_add_force with an index naming no force must return 0, got {}",
        ok
    );

    // Created but with no update_accelerations and no force_type.
    let bare = rebx_extras_mut(&mut sim)
        .map(|rebx| rebx_create_force(rebx, "custom"))
        .expect("extras attached");
    let ok = rebx_add_force(&mut sim, bare);
    assert_eq!(
        ok, 0,
        "rebx_add_force on a force with no update_accelerations must return 0, got {}",
        ok
    );
    assert!(
        sim.additional_forces.is_none(),
        "a refused rebx_add_force must not install the REBOUND hook"
    );

    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");
    assert_eq!(
        rebx.additional_forces.len(),
        0,
        "a refused rebx_add_force must not grow additional_forces, len = {}",
        rebx.additional_forces.len()
    );
    assert!(
        rebx.messages
            .iter()
            .any(|m| m.contains("update_accelerations")),
        "expected an 'update_accelerations' error, messages = {:?}",
        rebx.messages
    );

    // Give it a function but still no force_type: still refused.
    rebx.allocated_forces[bare].update_accelerations = Some(crate_noop_force);
    let ok = rebx_add_force(&mut sim, bare);
    assert_eq!(
        ok, 0,
        "rebx_add_force with force_type == REBX_FORCE_NONE must return 0, got {}",
        ok
    );

    // With both fields set it is accepted.
    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");
    rebx.allocated_forces[bare].force_type = rebx_force_type::REBX_FORCE_POS;
    let ok = rebx_add_force(&mut sim, bare);
    assert_eq!(
        ok, 1,
        "a fully configured custom force must be accepted, got {}",
        ok
    );
}

fn crate_noop_force(
    _sim: &mut reb_simulation,
    _rebx: &mut rebx_extras,
    _force_idx: usize,
    _N: usize,
) {
}

#[test]
fn add_force_sets_the_velocity_dependent_flag_only_for_vel_forces() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    assert_eq!(
        sim.force_is_velocity_dependent, 0,
        "a fresh simulation must start with force_is_velocity_dependent = 0, got {}",
        sim.force_is_velocity_dependent
    );

    let cf = rebx_load_force(&mut sim, "central_force").expect("central_force in library");
    rebx_add_force(&mut sim, cf);
    assert_eq!(
        sim.force_is_velocity_dependent, 0,
        "central_force is REBX_FORCE_POS, so force_is_velocity_dependent must stay 0, got {}",
        sim.force_is_velocity_dependent
    );

    let gr = rebx_load_force(&mut sim, "gr").expect("gr in library");
    rebx_add_force(&mut sim, gr);
    assert_eq!(
        sim.force_is_velocity_dependent, 1,
        "gr is REBX_FORCE_VEL, so force_is_velocity_dependent must become 1, got {}",
        sim.force_is_velocity_dependent
    );
}

#[test]
fn additional_forces_list_is_traversed_newest_first() {
    // rebx_add_force prepends, so the acceleration summation order is the
    // reverse of the order the forces were added. That order is
    // load-bearing for bit-for-bit agreement with the C.
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let a = rebx_load_force(&mut sim, "central_force").expect("in library");
    let b = rebx_load_force(&mut sim, "gr_potential").expect("in library");
    let c = rebx_load_force(&mut sim, "gravitational_harmonics").expect("in library");
    for f in [a, b, c] {
        assert_eq!(rebx_add_force(&mut sim, f), 1, "rebx_add_force({}) failed", f);
    }
    let rebx = rebx_extras_ref(&sim).expect("extras attached");
    assert_eq!(
        rebx.additional_forces,
        vec![c, b, a],
        "additional_forces must read back newest-first, got {:?}, expected {:?}",
        rebx.additional_forces,
        vec![c, b, a]
    );
}

#[test]
fn remove_force_unhooks_it_from_the_simulation() {
    let mut sim = two_body(1.0);
    rebx_attach(&mut sim);
    let a = rebx_load_force(&mut sim, "central_force").expect("in library");
    let b = rebx_load_force(&mut sim, "gr_potential").expect("in library");
    rebx_add_force(&mut sim, a);
    rebx_add_force(&mut sim, b);

    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");
    let removed = rebx_remove_force(rebx, a);
    assert_eq!(
        removed, 1,
        "removing a force that is in additional_forces must return 1, got {}",
        removed
    );
    assert_eq!(
        rebx.additional_forces,
        vec![b],
        "after removing force {} only {} should remain, got {:?}",
        a,
        b,
        rebx.additional_forces
    );
    let removed_again = rebx_remove_force(rebx, a);
    assert_eq!(
        removed_again, 0,
        "removing the same force twice must return 0 the second time, got {}",
        removed_again
    );
}

#[test]
fn add_operator_picks_the_timing_that_suits_the_integrator() {
    // WHFast: an updater is split into a half step before and after.
    let mut sim = two_body(1.0);
    reb_simulation_set_integrator(&mut sim, "whfast");
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    assert_eq!(rebx_add_operator(&mut sim, op), 1, "rebx_add_operator failed under whfast");
    {
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        assert_eq!(
            (rebx.pre_timestep_modifications.len(), rebx.post_timestep_modifications.len()),
            (1, 1),
            "under WHFast an updater must add one pre and one post step, got ({}, {})",
            rebx.pre_timestep_modifications.len(),
            rebx.post_timestep_modifications.len()
        );
        for (label, s) in [
            ("pre", rebx.pre_timestep_modifications[0]),
            ("post", rebx.post_timestep_modifications[0]),
        ] {
            assert_eq!(
                s.dt_fraction, 0.5,
                "the {} step's dt_fraction must be 1/2 under WHFast, got {}",
                label, s.dt_fraction
            );
            assert_eq!(
                s.operator_, op,
                "the {} step must reference operator {}, got {}",
                label, op, s.operator_
            );
        }
    }

    // IAS15: only a full step afterwards, because the step size is not
    // known in advance.
    let mut sim = two_body(1.0);
    reb_simulation_set_integrator(&mut sim, "ias15");
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    assert_eq!(rebx_add_operator(&mut sim, op), 1, "rebx_add_operator failed under ias15");
    {
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        assert_eq!(
            rebx.pre_timestep_modifications.len(),
            0,
            "under IAS15 no pre-timestep step may be added, got {}",
            rebx.pre_timestep_modifications.len()
        );
        assert_eq!(
            rebx.post_timestep_modifications.len(),
            1,
            "under IAS15 exactly one post-timestep step is expected, got {}",
            rebx.post_timestep_modifications.len()
        );
        assert_eq!(
            rebx.post_timestep_modifications[0].dt_fraction, 1.0,
            "the IAS15 post step must use the full timestep, got {}",
            rebx.post_timestep_modifications[0].dt_fraction
        );
    }

    // A recorder is added once after the step whatever the integrator.
    let mut sim = two_body(1.0);
    reb_simulation_set_integrator(&mut sim, "whfast");
    rebx_attach(&mut sim);
    let rec = rebx_load_operator(&mut sim, "track_min_distance").expect("in library");
    assert_eq!(rebx_add_operator(&mut sim, rec), 1, "adding a recorder failed");
    {
        let rebx = rebx_extras_ref(&sim).expect("extras attached");
        assert_eq!(
            (rebx.pre_timestep_modifications.len(), rebx.post_timestep_modifications.len()),
            (0, 1),
            "a recorder must add exactly one post step and no pre step, got ({}, {})",
            rebx.pre_timestep_modifications.len(),
            rebx.post_timestep_modifications.len()
        );
        assert_eq!(
            rebx.post_timestep_modifications[0].dt_fraction, 1.0,
            "a recorder's step must use the full timestep, got {}",
            rebx.post_timestep_modifications[0].dt_fraction
        );
    }
}

#[test]
fn remove_operator_unhooks_every_step_that_uses_it() {
    let mut sim = two_body(1.0);
    reb_simulation_set_integrator(&mut sim, "whfast");
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    rebx_add_operator(&mut sim, op);

    let rebx = rebx_extras_mut(&mut sim).expect("extras attached");
    let removed = rebx_remove_operator(rebx, op);
    assert_eq!(
        removed, 1,
        "removing an installed operator must return 1, got {}",
        removed
    );
    assert_eq!(
        (rebx.pre_timestep_modifications.len(), rebx.post_timestep_modifications.len()),
        (0, 0),
        "rebx_remove_operator must drop BOTH the pre and post steps, got ({}, {})",
        rebx.pre_timestep_modifications.len(),
        rebx.post_timestep_modifications.len()
    );
    let removed_again = rebx_remove_operator(rebx, op);
    assert_eq!(
        removed_again, 0,
        "removing the same operator twice must return 0 the second time, got {}",
        removed_again
    );
}

// ---------------------------------------------------------------------
// The timestep hooks
// ---------------------------------------------------------------------

#[test]
fn pre_hook_uses_dt_and_post_hook_uses_dt_last_done() {
    // core.c takes `sim->dt` in rebx_pre_timestep_modifications but
    // `sim->dt_last_done` in rebx_post_timestep_modifications, and scales
    // each by the step's dt_fraction. modify_mass makes that visible:
    // it does exactly `p.m += p.m*dt/tau_mass`.
    let mut sim = two_body(1.0);
    sim.particles[1].m = 1.0e-3;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = 0.1;
    sim.dt_last_done = 0.0;

    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    assert_eq!(
        rebx_add_operator_step(&mut sim, op, 0.5, rebx_timing::REBX_TIMING_PRE),
        1,
        "adding the PRE step failed"
    );
    assert_eq!(
        rebx_add_operator_step(&mut sim, op, 0.5, rebx_timing::REBX_TIMING_POST),
        1,
        "adding the POST step failed"
    );

    let tau = -500.0_f64;
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_mass", tau);
    }

    let m0 = sim.particles[1].m;
    let star0 = sim.particles[0].m;

    // dt_last_done is still zero, so the post hook must be a no-op.
    rebx_post_timestep_modifications(&mut sim);
    assert_eq!(
        sim.particles[1].m.to_bits(),
        m0.to_bits(),
        "with dt_last_done = 0 the post hook must leave the mass untouched: \
         {:e} -> {:e}",
        m0,
        sim.particles[1].m
    );

    // The pre hook uses sim.dt * dt_fraction = 0.1 * 0.5.
    rebx_pre_timestep_modifications(&mut sim);
    let expect1 = m0 + m0 * (0.1_f64 * 0.5) / tau;
    assert_eq!(
        sim.particles[1].m.to_bits(),
        expect1.to_bits(),
        "pre hook: expected m = {:e} (bits {:016x}), got {:e} (bits {:016x})",
        expect1,
        expect1.to_bits(),
        sim.particles[1].m,
        sim.particles[1].m.to_bits()
    );

    // Now give dt_last_done a different value from dt; the post hook must
    // use dt_last_done.
    sim.dt_last_done = 0.04;
    let m1 = sim.particles[1].m;
    rebx_post_timestep_modifications(&mut sim);
    let expect2 = m1 + m1 * (0.04_f64 * 0.5) / tau;
    assert_eq!(
        sim.particles[1].m.to_bits(),
        expect2.to_bits(),
        "post hook must scale dt_last_done (0.04), not dt (0.1): expected {:e}, got {:e}",
        expect2,
        sim.particles[1].m
    );

    // tau_mass was never set on the star, so its mass is untouched.
    assert_eq!(
        sim.particles[0].m.to_bits(),
        star0.to_bits(),
        "the star has no tau_mass, so its mass must not change: {:e} -> {:e}",
        star0,
        sim.particles[0].m
    );
    // Mass loss (tau < 0) must be monotonic.
    assert!(
        sim.particles[1].m < m1 && m1 < m0,
        "with tau_mass < 0 the mass must decrease monotonically: {:e} -> {:e} -> {:e}",
        m0,
        m1,
        sim.particles[1].m
    );
}

// ---------------------------------------------------------------------
// The extras box survives a full integration
// ---------------------------------------------------------------------

#[test]
fn extras_survives_integration_with_a_force() {
    let mut sim = two_body(1.0);
    sim.particles[1].m = 1.0e-4;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = 0.01;
    rebx_attach(&mut sim);
    let cf = rebx_load_force(&mut sim, "central_force").expect("in library");
    assert_eq!(rebx_add_force(&mut sim, cf), 1, "rebx_add_force failed");
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", 1.0e-3);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "gammacentral", -2.0);
        rebx_set_param_double(rebx, rebx_ap::force(cf), "tau_a", 7.0);
    }

    reb_simulation_integrate(&mut sim, 5.0);

    let rebx = rebx_extras_ref(&sim).expect("the extras box must survive reb_simulation_integrate");
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(0), "Acentral"),
        Some(1.0e-3),
        "the particle parameter must be intact after integration"
    );
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::force(cf), "tau_a"),
        Some(7.0),
        "the force parameter must be intact after integration"
    );
    assert_eq!(
        rebx.additional_forces,
        vec![cf],
        "additional_forces must still hold the force after integration, got {:?}",
        rebx.additional_forces
    );
    assert!(
        sim.t >= 5.0,
        "the integration should have reached t >= 5.0, got t = {:e}",
        sim.t
    );
    assert!(
        sim.particles[1].x.is_finite() && sim.particles[1].y.is_finite(),
        "the companion's position must stay finite: ({:e}, {:e})",
        sim.particles[1].x,
        sim.particles[1].y
    );
}

#[test]
fn extras_survives_integration_with_pre_and_post_operators() {
    let mut sim = two_body(1.0);
    sim.particles[1].m = 1.0e-3;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = 0.01;
    rebx_attach(&mut sim);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    assert_eq!(rebx_add_operator(&mut sim, op), 1, "rebx_add_operator failed");

    let tau = -50.0_f64;
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_mass", tau);
    }
    let m0 = sim.particles[1].m;
    let star0 = sim.particles[0].m;

    reb_simulation_integrate(&mut sim, 2.0);

    let rebx = rebx_extras_ref(&sim).expect("the extras box must survive the operator hooks");
    assert_eq!(
        rebx_get_param_double(rebx, rebx_ap::particle(1), "tau_mass"),
        Some(tau),
        "tau_mass must be intact after integration"
    );

    // Exponential loss with e-folding time |tau| = 50 over t = 2 gives a
    // factor of exp(-2/50). The operator integrates that with the
    // first-order explicit form m *= (1 + h/tau), applied twice per
    // timestep (a half step before and after), so
    //     h = dt/2 = 0.005, x = |h/tau| = 1e-4, n = 2*(2/dt) = 400.
    // Since ln(1-x) < -x, the discrete product (1-x)^n falls BELOW
    // exp(-nx), by a relative deficit of n*x^2/2 = 2.0e-6 to leading
    // order. Both the sign and the size are checked.
    let m1 = sim.particles[1].m;
    let exact = m0 * (-2.0_f64 / 50.0).exp();
    assert!(
        m1 < m0,
        "with tau_mass < 0 the mass must decrease: {:e} -> {:e}",
        m0,
        m1
    );
    assert!(
        m1 < exact,
        "a first-order explicit decay must undershoot the exact exponential \
         (ln(1-x) < -x): m = {:e}, exp = {:e}",
        m1,
        exact
    );
    let rel = (exact - m1) / exact;
    let x = 0.005_f64 / 50.0;
    let n = 2.0 * (2.0 / 0.01);
    let predicted = n * x * x / 2.0;
    assert!(
        (rel / predicted - 1.0).abs() < 0.25,
        "the discretisation deficit should be n*x^2/2 = {:e} (n = {}, x = {:e}); \
         measured {:e} (m = {:e}, exp = {:e})",
        predicted,
        n,
        x,
        rel,
        m1,
        exact
    );
    assert_eq!(
        sim.particles[0].m.to_bits(),
        star0.to_bits(),
        "the star has no tau_mass and must keep its mass exactly: {:e} -> {:e}",
        star0,
        sim.particles[0].m
    );
}

#[test]
fn an_effectless_force_reproduces_plain_gravity_bit_for_bit() {
    // central_force is a no-op unless a particle carries BOTH Acentral
    // and gammacentral, and it is REBX_FORCE_POS so it does not flip
    // force_is_velocity_dependent. Installing it without parameters must
    // therefore leave every bit of the trajectory alone.
    let mut plain = two_body(1.0);
    plain.particles[1].m = 1.0e-4;
    reb_simulation_set_integrator(&mut plain, "whfast");
    plain.dt = 0.012_345;
    reb_simulation_integrate(&mut plain, 5.0);

    let mut hooked = two_body(1.0);
    hooked.particles[1].m = 1.0e-4;
    reb_simulation_set_integrator(&mut hooked, "whfast");
    hooked.dt = 0.012_345;
    rebx_attach(&mut hooked);
    let cf = rebx_load_force(&mut hooked, "central_force").expect("in library");
    assert_eq!(rebx_add_force(&mut hooked, cf), 1, "rebx_add_force failed");
    // Only Acentral, deliberately without gammacentral: the C requires
    // both before it does any work.
    if let Some(rebx) = rebx_extras_mut(&mut hooked) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", 1.0);
    }
    reb_simulation_integrate(&mut hooked, 5.0);

    assert_eq!(
        hooked.force_is_velocity_dependent, 0,
        "central_force is REBX_FORCE_POS, so the flag must stay 0, got {}",
        hooked.force_is_velocity_dependent
    );
    let a = state_bits(&plain);
    let b = state_bits(&hooked);
    assert_eq!(
        a, b,
        "a parameterless central_force must reproduce plain gravity bit-for-bit;\n\
         plain  = {:?}\nhooked = {:?}",
        a, b
    );
}

#[test]
fn central_force_can_stand_in_for_a_missing_central_mass() {
    // central_force with gamma = -2 gives prefac = A*r^-3 and adds
    // prefac*dx, i.e. A*dx/r^3. Newtonian gravity from a body of mass M
    // at the same place gives -G*M*dx/r^3. So A = -G*M reproduces the
    // gravity of an extra central mass M exactly. Splitting the star's
    // mass between real gravity and the effect must therefore give the
    // same orbit as the undivided star.
    let mut whole = two_body(1.0);
    reb_simulation_set_integrator(&mut whole, "ias15");
    reb_simulation_integrate(&mut whole, 6.0);

    let mut split = two_body(0.5);
    reb_simulation_set_integrator(&mut split, "ias15");
    rebx_attach(&mut split);
    let cf = rebx_load_force(&mut split, "central_force").expect("in library");
    assert_eq!(rebx_add_force(&mut split, cf), 1, "rebx_add_force failed");
    if let Some(rebx) = rebx_extras_mut(&mut split) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", -0.5);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "gammacentral", -2.0);
    }
    reb_simulation_integrate(&mut split, 6.0);

    // The companion is massless, so neither star ever moves.
    for (label, s) in [("whole", &whole), ("split", &split)] {
        assert!(
            s.particles[0].x.abs() < 1.0e-15 && s.particles[0].y.abs() < 1.0e-15,
            "the {} run's star should stay pinned at the origin (massless companion), \
             got ({:e}, {:e})",
            label,
            s.particles[0].x,
            s.particles[0].y
        );
    }

    let dx = whole.particles[1].x - split.particles[1].x;
    let dy = whole.particles[1].y - split.particles[1].y;
    let dz = whole.particles[1].z - split.particles[1].z;
    let sep = (dx * dx + dy * dy + dz * dz).sqrt();
    let r = (whole.particles[1].x.powi(2) + whole.particles[1].y.powi(2)).sqrt();
    assert!(
        sep / r < 1.0e-8,
        "half the star's mass moved into central_force must reproduce the same orbit \
         after t = 6: separation {:e} at r = {:e} (relative {:e}), tolerance 1e-8",
        sep,
        r,
        sep / r
    );

    // Control: without the effect, half a star is not enough to bind the
    // companion at all (v^2 = 1.21 > 2GM/r = 1.0), so it escapes.
    let mut halved = two_body(0.5);
    reb_simulation_set_integrator(&mut halved, "ias15");
    reb_simulation_integrate(&mut halved, 6.0);
    let r_esc = (halved.particles[1].x.powi(2) + halved.particles[1].y.powi(2)).sqrt();
    assert!(
        r_esc > 2.0 * r,
        "control: with no central_force the companion is unbound and should have run \
         far away — r = {:e}, bound run r = {:e}",
        r_esc,
        r
    );
}

// ---------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------

fn determinism_run() -> reb_simulation {
    let mut sim = two_body(1.0);
    sim.particles[1].m = 1.0e-3;
    reb_simulation_set_integrator(&mut sim, "whfast");
    sim.dt = 0.007;
    rebx_attach(&mut sim);

    let cf = rebx_load_force(&mut sim, "central_force").expect("in library");
    rebx_add_force(&mut sim, cf);
    let op = rebx_load_operator(&mut sim, "modify_mass").expect("in library");
    rebx_add_operator(&mut sim, op);
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::particle(0), "Acentral", 2.5e-3);
        rebx_set_param_double(rebx, rebx_ap::particle(0), "gammacentral", -2.0);
        rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_mass", -80.0);
    }
    reb_simulation_integrate(&mut sim, 3.0);
    sim
}

#[test]
fn the_same_setup_integrates_bit_identically_twice() {
    let a = determinism_run();
    let b = determinism_run();

    let ba = state_bits(&a);
    let bb = state_bits(&b);
    assert_eq!(
        ba.len(),
        bb.len(),
        "the two runs produced different particle counts: {} vs {}",
        a.N,
        b.N
    );
    for (k, (x, y)) in ba.iter().zip(bb.iter()).enumerate() {
        assert_eq!(
            x, y,
            "state word {} differs between two identical runs: {:016x} vs {:016x}",
            k, x, y
        );
    }

    // The parameters the effects read must also come back identical.
    let pa = rebx_extras_ref(&a)
        .and_then(|r| rebx_get_param_double(r, rebx_ap::particle(1), "tau_mass"))
        .expect("tau_mass survives run A");
    let pb = rebx_extras_ref(&b)
        .and_then(|r| rebx_get_param_double(r, rebx_ap::particle(1), "tau_mass"))
        .expect("tau_mass survives run B");
    assert_eq!(
        pa.to_bits(),
        pb.to_bits(),
        "tau_mass differs between runs: {:e} vs {:e}",
        pa,
        pb
    );

    // And the effect really did something: the run is not the trivial
    // one where nothing changed.
    assert!(
        a.particles[1].m < 1.0e-3,
        "modify_mass should have reduced the companion's mass below 1e-3, got {:e}",
        a.particles[1].m
    );
}
