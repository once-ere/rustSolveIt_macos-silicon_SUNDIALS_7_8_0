//! simulation.rs — simulation life cycle, the main timestep and the
//! integrate loop (from simulation.c; OpenGL/SERVER/MPI/EMSCRIPTEN
//! branches excluded, exactly the subsystems the C build can disable).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein and contributors. See crate root.

use crate::boundary::reb_boundary_check;
use crate::collision::reb_collision_search;
use crate::gravity::*;
use crate::particle::reb_particle_check_testparticles;
use crate::tools::{
    reb_simulation_error, reb_simulation_rescale_var, reb_simulation_warning,
    reb_tools_get_rand_seed,
};
use crate::types::*;

/// simulation.c `reb_simulation_create` — allocates and initializes a
/// simulation with the C default values.
pub fn reb_simulation_create() -> reb_simulation {
    let mut r = reb_simulation {
        t: 0.,
        G: 1.,
        softening: 0.,
        OMEGA: 0.,
        OMEGAZ: -1.0,
        dt: 0.001,
        dt_last_done: 0.,
        steps_done: 0,
        is_synchronized: 1,
        did_modify_particles: 0,
        N: 0,
        particles: Vec::new(),
        N_map: 0,
        map: None,
        N_var: 0,
        particles_var: Vec::new(),
        var_config: Vec::new(),
        odes: Vec::new(),
        ode_id_next: 1,
        N_active: usize::MAX,
        testparticle_type: 0,
        testparticle_hidewarnings: 0,
        name_list: Vec::new(),
        gravity_cs: Vec::new(),
        tree_root: Vec::new(),
        tree_cells: Vec::new(),
        opening_angle2: 0.25,
        status: 0,
        exact_finish_time: 1,
        force_is_velocity_dependent: 0,
        gravity_ignore_terms: REB_GRAVITY_IGNORE_TERMS_NONE,
        output_timing_last: -1.,
        save_messages: 0,
        messages: Vec::new(),
        messages_var_rescale_warning: 0,
        messages_timestep_warning: 0,
        exit_max_distance: 0.,
        exit_min_distance: 0.,
        usleep: 0.,
        track_energy_offset: 0,
        energy_offset: 0.,
        walltime: 0.,
        walltime_last_step: 0.,
        walltime_last_steps: 0.,
        walltime_last_steps_sum: 0.,
        walltime_last_steps_N: 0,
        root_size: -1.,
        N_root_x: 1,
        N_root_y: 1,
        N_root_z: 1,
        N_ghost_x: 0,
        N_ghost_y: 0,
        N_ghost_z: 0,
        collisions: Vec::new(),
        N_collisions: 0,
        N_targets: usize::MAX,
        minimum_collision_velocity: 0.,
        collisions_plog: 0.,
        collisions_log_n: 0,
        calculate_megno: 0,
        megno_Ys: 0.,
        megno_Yss: 0.,
        megno_cov_Yt: 0.,
        megno_var_t: 0.,
        megno_mean_t: 0.,
        megno_mean_Y: 0.,
        megno_initial_t: 0.,
        megno_n: 0,
        simulationarchive_version: 5,
        simulationarchive_auto_interval: 0.,
        simulationarchive_auto_walltime: 0.,
        simulationarchive_auto_step: 0,
        simulationarchive_next: 0.,
        simulationarchive_next_step: 0,
        simulationarchive_filename: None,
        python_unit_l: 0,
        python_unit_m: 0,
        python_unit_t: 0,
        rand_seed:reb_tools_get_rand_seed(),
        collision: REB_COLLISION::NONE,
        boundary: REB_BOUNDARY::NONE,
        gravity: REB_GRAVITY::BASIC,
        integrator: reb_integrator_state::none,
        gravity_custom: None,
        additional_forces: None,
        pre_timestep_modifications: None,
        post_timestep_modifications: None,
        heartbeat: None,
        coefficient_of_restitution: None,
        collision_resolve: None,
        server_data: None,
        extras: None,
    };
    reb_simulation_set_integrator(&mut r, "ias15");
    r
}

/// simulation.c `reb_simulation_stop`.
pub fn reb_simulation_stop(r: &mut reb_simulation) {
    r.status = REB_STATUS_USER;
}

/// simulation.c `reb_simulation_set_integrator`.
pub fn reb_simulation_set_integrator(r: &mut reb_simulation, name: &str) {
    if r.is_synchronized == 0 {
        reb_simulation_warning(
            r,
            "Changing integrators while simulation is not synchronized results in undefined behaviour.",
        );
    }
    match name {
        "none" => r.integrator = reb_integrator_state::none,
        "sei" => {
            r.integrator = reb_integrator_state::sei(
                crate::integrator_sei::reb_integrator_sei_state::default(),
            )
        }
        "leapfrog" => {
            r.integrator = reb_integrator_state::leapfrog(
                crate::integrator_leapfrog::reb_integrator_leapfrog_state::default(),
            )
        }
        "ias15" => {
            r.integrator = reb_integrator_state::ias15(
                crate::integrator_ias15::reb_integrator_ias15_state::default(),
            )
        }
        "whfast" => {
            r.integrator = reb_integrator_state::whfast(
                crate::integrator_whfast::reb_integrator_whfast_state::default(),
            )
        }
        "saba" => {
            r.integrator = reb_integrator_state::saba(
                crate::integrator_saba::reb_integrator_saba_state::default(),
            )
        }
        "janus" => {
            r.integrator = reb_integrator_state::janus(
                crate::integrator_janus::reb_integrator_janus_state::default(),
            )
        }
        "eos" => {
            r.integrator = reb_integrator_state::eos(
                crate::integrator_eos::reb_integrator_eos_state::default(),
            )
        }
        "mercurius" => {
            r.integrator = reb_integrator_state::mercurius(
                crate::integrator_mercurius::reb_integrator_mercurius_state::default(),
            )
        }
        "bs" => {
            r.integrator = reb_integrator_state::bs(
                crate::integrator_bs::reb_integrator_bs_state::default(),
            )
        }
        "trace" => {
            r.integrator = reb_integrator_state::trace(
                crate::integrator_trace::reb_integrator_trace_state::default(),
            )
        }
        "whfast512" => {
            r.integrator = reb_integrator_state::whfast512(
                crate::integrator_whfast512::reb_integrator_whfast512_state::default(),
            )
        }
        _ => reb_simulation_error(r, "Integrator not found."),
    }
}

/// Integrator `did_add_particle` hook dispatch (C: the
/// `did_add_particle` member of the integrator vtable; of the ported
/// integrators only MERCURIUS registers one).
pub fn reb_integrator_did_add_particle(r: &mut reb_simulation) {
    if matches!(r.integrator, reb_integrator_state::mercurius(_)) {
        crate::integrator_mercurius::reb_integrator_mercurius_did_add_particle(r);
    }
    if matches!(r.integrator, reb_integrator_state::trace(_)) {
        crate::integrator_trace::reb_integrator_trace_did_add_particle(r);
    }
}

/// Integrator `will_remove_particle` hook dispatch (same situation).
pub fn reb_integrator_will_remove_particle(r: &mut reb_simulation, index: usize) {
    if matches!(r.integrator, reb_integrator_state::mercurius(_)) {
        crate::integrator_mercurius::reb_integrator_mercurius_will_remove_particle(r, index);
    }
    if matches!(r.integrator, reb_integrator_state::trace(_)) {
        crate::integrator_trace::reb_integrator_trace_will_remove_particle(r, index);
    }
}

/// simulation.c `run_heartbeat` — heartbeat wrapper with exit checks.
fn run_heartbeat(r: &mut reb_simulation) {
    if let Some(hb) = r.heartbeat {
        hb(r);
    }
    if r.exit_max_distance != 0. {
        // Check for escaping particles
        let max2 = r.exit_max_distance * r.exit_max_distance;
        let N = r.N;
        for i in 0..N {
            let p = r.particles[i];
            let r2 = p.x * p.x + p.y * p.y + p.z * p.z;
            if r2 > max2 {
                r.status = REB_STATUS_ESCAPE;
            }
        }
    }
    if r.exit_min_distance != 0. {
        // Check for close encounters
        let min2 = r.exit_min_distance * r.exit_min_distance;
        let N = r.N;
        for i in 0..N {
            let pi = r.particles[i];
            for j in 0..i {
                let pj = r.particles[j];
                let x = pi.x - pj.x;
                let y = pi.y - pj.y;
                let z = pi.z - pj.z;
                let r2 = x * x + y * y + z * z;
                if r2 < min2 {
                    r.status = REB_STATUS_ENCOUNTER;
                }
            }
        }
    }
}

fn error_message_waiting(r: &reb_simulation) -> bool {
    for (t, _) in &r.messages {
        if *t == REB_MESSAGE_TYPE::ERROR {
            return true;
        }
    }
    false
}

/// simulation.c `reb_check_exit`.
fn reb_check_exit(r: &mut reb_simulation, tmax: f64, last_full_dt: &mut f64) -> REB_STATUS {
    if r.status <= REB_STATUS_SINGLE_STEP {
        if r.status == REB_STATUS_SINGLE_STEP {
            r.status = REB_STATUS_PAUSED;
        } else {
            // This allows an arbitrary number of steps before pausing
            r.status += 1;
        }
    }
    while r.status == REB_STATUS_PAUSED || r.status == REB_STATUS_SCREENSHOT {
        // Wait for user to disable paused simulation (the C server
        // thread flips the status directly; here the queued keyboard
        // commands are applied while waiting).
        crate::server::reb_server_update(r);
        std::thread::sleep(std::time::Duration::from_millis(1));
    }
    let dtsign = 1.0_f64.copysign(r.dt); // Used to determine integration direction
    if error_message_waiting(r) {
        r.status = REB_STATUS_GENERIC_ERROR;
    }
    if r.status >= 0 {
        // Exit now.
    } else if tmax != f64::INFINITY {
        if r.exact_finish_time == 1 {
            if (r.t + r.dt) * dtsign >= tmax * dtsign {
                // Next step would overshoot
                if r.t == tmax {
                    r.status = REB_STATUS_SUCCESS;
                } else if r.status == REB_STATUS_LAST_STEP {
                    let mut tscale = 1e-12 * tmax.abs(); // Find order of magnitude for time
                    if tscale < 1e-200 {
                        // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (r.t - tmax).abs() < tscale {
                        r.status = REB_STATUS_SUCCESS;
                    } else {
                        // not there yet, do another step.
                        reb_simulation_synchronize(r);
                        r.dt = tmax - r.t;
                    }
                } else {
                    r.status = REB_STATUS_LAST_STEP; // Do one small step, then exit.
                    reb_simulation_synchronize(r);
                    if r.dt_last_done != 0. {
                        // If first timestep is also last, do not use dt_last_done
                        *last_full_dt = r.dt_last_done;
                    }
                    r.dt = tmax - r.t;
                }
            } else if r.status == REB_STATUS_LAST_STEP {
                // An adaptive integrator reduced the timestep in what was
                // supposed to be the last timestep.
                r.status = REB_STATUS_RUNNING;
            }
        } else if r.t * dtsign >= tmax * dtsign {
            // Past the integration time
            r.status = REB_STATUS_SUCCESS; // Exit now.
        }
    }
    if r.N == 0 && r.odes.is_empty() {
        reb_simulation_warning(r, "No particles in simulation. Will exit.");
        r.status = REB_STATUS_NO_PARTICLES; // Exit now.
    }
    r.status
}

/// simulation.c `reb_simulation_step` — one bare timestep, without the
/// heartbeat.
pub fn reb_simulation_step(r: &mut reb_simulation) {
    let time_beginning = std::time::Instant::now();

    if r.pre_timestep_modifications.is_some() {
        reb_simulation_synchronize(r);
        (r.pre_timestep_modifications.unwrap())(r);
    }

    match r.integrator {
        reb_integrator_state::none => crate::integrator_none::reb_integrator_none_step(r),
        reb_integrator_state::sei(_) => crate::integrator_sei::reb_integrator_sei_step(r),
        reb_integrator_state::leapfrog(_) => {
            crate::integrator_leapfrog::reb_integrator_leapfrog_step(r)
        }
        reb_integrator_state::ias15(_) => crate::integrator_ias15::reb_integrator_ias15_step(r),
        reb_integrator_state::whfast(_) => crate::integrator_whfast::reb_integrator_whfast_step(r),
        reb_integrator_state::saba(_) => crate::integrator_saba::reb_integrator_saba_step(r),
        reb_integrator_state::janus(_) => crate::integrator_janus::reb_integrator_janus_step(r),
        reb_integrator_state::eos(_) => crate::integrator_eos::reb_integrator_eos_step(r),
        reb_integrator_state::mercurius(_) => {
            crate::integrator_mercurius::reb_integrator_mercurius_step(r)
        }
        reb_integrator_state::bs(_) => crate::integrator_bs::reb_integrator_bs_step(r),
        reb_integrator_state::trace(_) => crate::integrator_trace::reb_integrator_trace_step(r),
        reb_integrator_state::whfast512(_) => {
            crate::integrator_whfast512::reb_integrator_whfast512_step(r)
        }
    }

    // Integrate other ODEs (simulation.c: user ODEs are advanced with a
    // temporary BS state when the main integrator is not BS).
    if !r.odes.is_empty() && !matches!(r.integrator, reb_integrator_state::bs(_)) {
        let mut dt = r.dt_last_done;
        let mut t = r.t - r.dt_last_done; // Note: floating point inaccuracy
        let forward = if dt > 0. { 1. } else { -1. };
        let mut bs = crate::integrator_bs::reb_integrator_bs_state::default();
        while t * forward < r.t * forward && ((r.t - t) / (r.t.abs() + 1e-16)).abs() > 1e-15 {
            if bs.dt_proposed != 0. {
                let max_dt = (r.t - t).abs();
                dt = bs.dt_proposed.abs();
                if dt > max_dt {
                    // Don't overshoot N-body timestep
                    dt = max_dt;
                    bs.first_or_last_step = 1;
                }
                dt *= forward;
            }
            let success = crate::integrator_bs::reb_integrator_bs_step_odes(r, &mut bs, dt);
            if success != 0 {
                t += dt;
            }
        }
    }

    if r.post_timestep_modifications.is_some() {
        reb_simulation_synchronize(r);
        (r.post_timestep_modifications.unwrap())(r);
    }

    // Reset tainted particle flag
    r.did_modify_particles = 0;

    if r.N_var != 0 {
        reb_simulation_rescale_var(r);
    }

    reb_boundary_check(r);

    if r.collision != REB_COLLISION::NONE {
        reb_collision_search(r);
    }

    // Update walltime
    let elapsed = time_beginning.elapsed();
    r.walltime_last_step = elapsed.as_secs() as f64 + (elapsed.subsec_micros() as f64) / 1e6;
    r.walltime_last_steps_sum += r.walltime_last_step;
    r.walltime_last_steps_N += 1;
    if r.walltime_last_steps_sum > 0.1 {
        r.walltime_last_steps = r.walltime_last_steps_sum / (r.walltime_last_steps_N as f64);
        r.walltime_last_steps_sum = 0.;
        r.walltime_last_steps_N = 0;
    }
    r.walltime += r.walltime_last_step;

    // Update step counter
    r.steps_done += 1; // This also counts failed IAS15 steps
}

/// simulation.c `reb_simulation_integrate_raw` + `reb_simulation_integrate`
/// (the compute runs on the calling thread; the OpenGL display thread
/// split does not exist here).
pub fn reb_simulation_integrate(r: &mut reb_simulation, tmax: f64) -> REB_STATUS {
    if tmax != r.t {
        let dt_sign = if tmax > r.t { 1.0 } else { -1.0 };
        r.dt = r.dt.copysign(dt_sign);
    }

    let mut last_full_dt = r.dt; // stored in case dt is shrunk for exact_finish_time
    r.dt_last_done = 0.; // Reset in case first timestep attempt will fail

    if r.testparticle_hidewarnings == 0 && reb_particle_check_testparticles(r) {
        reb_simulation_warning(r, "At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }
    if r.status != REB_STATUS_PAUSED && r.status != REB_STATUS_SCREENSHOT {
        r.status = REB_STATUS_RUNNING;
    }
    run_heartbeat(r);
    while reb_check_exit(r, tmax, &mut last_full_dt) < 0 {
        if r.server_data.is_some() {
            // C: the integrate loop holds the server mutex around the
            // step; here the snapshot/key handshake runs between steps.
            crate::server::reb_server_update(r);
        }
        if r.simulationarchive_filename.is_some() {
            crate::simulationarchive::reb_simulationarchive_heartbeat(r);
        }
        reb_simulation_step(r);
        run_heartbeat(r);
        if r.usleep > 0. {
            std::thread::sleep(std::time::Duration::from_micros(r.usleep as u64));
        }
    }
    reb_simulation_synchronize(r);
    if r.exact_finish_time == 1 {
        // dt could have been shrunk; restore the last full timestep
        r.dt = last_full_dt;
    }
    r.status
}

/// simulation.c `reb_simulation_steps`.
pub fn reb_simulation_steps(r: &mut reb_simulation, N_steps: usize) -> REB_STATUS {
    r.status = REB_STATUS_RUNNING;
    run_heartbeat(r);
    let mut i = 0;
    while i < N_steps && r.status < 0 {
        reb_simulation_step(r);
        run_heartbeat(r);
        i += 1;
    }
    reb_simulation_synchronize(r);
    if r.status <= 0 {
        r.status = REB_STATUS_SUCCESS; // No error occurred. Success.
    }
    r.status
}

/// simulation.c `reb_simulation_synchronize`. Of the built-in
/// integrators only WHFast defines a synchronize callback in C.
pub fn reb_simulation_synchronize(r: &mut reb_simulation) {
    match r.integrator {
        reb_integrator_state::whfast(_) => {
            crate::integrator_whfast::reb_integrator_whfast_synchronize(r)
        }
        reb_integrator_state::saba(_) => {
            crate::integrator_saba::reb_integrator_saba_synchronize(r)
        }
        reb_integrator_state::janus(_) => {
            crate::integrator_janus::reb_integrator_janus_synchronize(r)
        }
        reb_integrator_state::eos(_) => {
            crate::integrator_eos::reb_integrator_eos_synchronize(r)
        }
        reb_integrator_state::mercurius(_) => {
            crate::integrator_mercurius::reb_integrator_mercurius_synchronize(r)
        }
        reb_integrator_state::whfast512(_) => {
            crate::integrator_whfast512::reb_integrator_whfast512_synchronize(r)
        }
        _ => {}
    }
}

/// simulation.c `reb_simulation_update_acceleration`.
pub fn reb_simulation_update_acceleration(r: &mut reb_simulation) {
    if r.gravity == REB_GRAVITY::CUSTOM {
        match r.gravity_custom {
            None => {
                reb_simulation_error(
                    r,
                    "REB_GRAVITY_CUSTOM selected, but r->gravity_custom function pointer not provided.",
                );
            }
            Some(f) => f(r),
        }
        return;
    }
    match r.gravity {
        REB_GRAVITY::NONE => {
            for j in 0..r.N {
                r.particles[j].ax = 0.;
                r.particles[j].ay = 0.;
                r.particles[j].az = 0.;
            }
        }
        REB_GRAVITY::TREE => reb_gravity_tree_calculate_acceleration(r),
        REB_GRAVITY::JACOBI => reb_gravity_jacobi_calculate_acceleration(r),
        REB_GRAVITY::BASIC => reb_gravity_basic_calculate_acceleration(r),
        REB_GRAVITY::COMPENSATED => reb_gravity_compensated_calculate_acceleration(r),
        REB_GRAVITY::CUSTOM => unreachable!(),
    }
    if r.N_var != 0 {
        match r.gravity {
            REB_GRAVITY::BASIC => reb_gravity_basic_calculate_acceleration_var(r),
            REB_GRAVITY::NONE => {}
            _ => {
                reb_simulation_error(r, "Variational gravity calculation not implemented in selected gravity module. Please use REB_GRAVITY_BASIC.");
                return;
            }
        }
    }

    if let Some(af) = r.additional_forces {
        af(r);
    }
}
