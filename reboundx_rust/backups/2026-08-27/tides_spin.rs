//! tides_spin.rs — translation of REBOUNDx tides_spin.c
//! Add self-consistent spin, tidal and dynamical equations of motion for bodies with structure.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Authors: Tiger Lu <tiger.lu@yale.edu>, Dan Tamayo <tamayo.daniel@gmail.com>,
//! Hanno Rein <hanno.rein@utoronto.ca>.
//!
//! # $Tides$
//!
//! ======================= ===============================================
//! Authors                 Tiger Lu, Hanno Rein, D. Tamayo, Sam Hadden, Rosemary Mardling, Sarah Millholland, Gregory Laughlin
//! Implementation Paper    `Lu et al., 2023 <https://arxiv.org/abs/2303.00006>`_.
//! Based on                `Eggleton et al. 1998 <https://ui.adsabs.harvard.edu/abs/1998ApJ...499..853E/abstract>`_.
//! C Example               :ref:`c_example_tides_spin_pseudo_synchronization`, :ref:`c_example_tides_spin_migration_driven_obliquity_tides`, :ref:`c_example_tides_spin_kozai`.
//! Python Example          `SpinsIntro.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/SpinsIntro.ipynb>`_, `TidesSpinPseudoSynchronization.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesSpinPseudoSynchronization.ipynb>`_, `TidesSpinEarthMoon.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesSpinEarthMoon.ipynb>`_.
//! ======================= ===============================================
//!
//! This effect consistently tracks both the spin and orbital evolution of bodies
//! under constant-time lag tides raised on both the primary and on the orbiting
//! bodies. In all cases, we need to set masses for all the particles that will
//! feel these tidal forces. Particles with only mass are point particles.
//!
//! Particles are assumed to have structure (i.e - physical extent & distortion
//! from spin) if the following parameters are set: physical radius
//! `particles[i].r`, potential Love number of degree 2 `k2` (Q/(1-Q) in Eggleton
//! 1998), and the spin angular rotation frequency vector `Omega`.
//! If we wish to evolve a body's spin components, the fully dimensional moment
//! of inertia `I` must be set as well. If this parameter is not set, the spin
//! components will be stationary. Note that if the body is a test particle, this
//! is assumed to be the specific moment of inertia.
//! Finally, if we wish to consider the effects of tides raised on a specific
//! body, we must set the constant time lag `tau` as well.
//!
//! For spins that are synchronized with a circular orbit, the constant time lag
//! can be related to the tidal quality factor Q as tau = 1/(2*n*tau), with n the
//! orbital mean motion. See Lu et. al (in review) and Eggleton et. al (1998)
//! above for discussion.
//!
//! **Effect Parameters**
//!
//! None
//!
//! **Particle Parameters**
//!
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
//! k2 (float)                   Yes         Potential Love number of degree 2.
//! Omega (reb_vec3d)            Yes         Angular rotation frequency (Omega_x, Omega_y, Omega_z)
//! I (float)                    No          Moment of inertia (for test particles, assumed to be the specific MoI I/m)
//! tau (float)                  No          Constant time lag. If not set, defaults to 0
//! ============================ =========== ==================================================================
//!
//! # How the C's spin ODE maps onto `rebound_rs`
//!
//! The C registers a `struct reb_ode` with REBOUND's Gragg-Bulirsch-Stoer
//! framework and reaches REBOUNDx from the callbacks through
//! `ode->ref == sim` and then `sim->extras`. Here the callbacks already
//! receive `&mut reb_simulation`, so `ode->ref` has no counterpart and the
//! REBOUNDx state is taken out of `sim.extras` with
//! [`crate::core::rebx_with`], which always puts it back.
//!
//! A `struct reb_ode*` is identified by its `id` in `sim.odes` rather than
//! by its address, so `rebx_spin_initialize_ode` stores that id under the
//! force's `"ode"` parameter (C: `rebx_set_param_pointer(.., "ode", spin_ode)`).

use rebound_rs::integrator_bs::{reb_ode, reb_ode_create, reb_ode_derivatives_fn, reb_ode_free};
use rebound_rs::{
    reb_particle, reb_simulation, reb_simulation_error, reb_simulation_warning, reb_vec3d,
};

use crate::core::{
    rebx_get_param_double, rebx_get_param_vec3d, rebx_set_param_ode, rebx_set_param_vec3d,
    rebx_with,
};
use crate::types::{rebx_ap, rebx_extras};

/// C: `struct reb_vec3d rebx_calculate_spin_orbit_accelerations(
/// struct reb_particle* source, struct reb_particle* target, const double G,
/// const double k2, const double sigma, const struct reb_vec3d Omega)`.
///
/// Neither particle is modified, so the two C pointers become shared
/// references.
pub fn rebx_calculate_spin_orbit_accelerations(
    source: &reb_particle,
    target: &reb_particle,
    G: f64,
    k2: f64,
    sigma: f64,
    Omega: reb_vec3d,
) -> reb_vec3d {
    // All quantities associated with SOURCE
    // This is the quadrupole potential/tides raised on the SOURCE
    let ms = source.m;
    let Rs = source.r;
    let mt = target.m;
    let mtot = ms + mt;
    let mu_ij = ms * mt / mtot; // have already checked for 0 and inf
    let big_a = k2 * (Rs * Rs * Rs * Rs * Rs);

    // distance vector FROM j TO i
    let dx = source.x - target.x;
    let dy = source.y - target.y;
    let dz = source.z - target.z;
    let d2 = dx * dx + dy * dy + dz * dz;
    let dr = d2.sqrt();

    // Velocity vector: i to j
    let dvx = source.vx - target.vx;
    let dvy = source.vy - target.vy;
    let dvz = source.vz - target.vz;
    //const double vel2 = dvx * dvx + dvy * dvy + dvz * dvz;
    //const double vr = sqrt(vel2);

    // C: `struct reb_vec3d tot_force = {0};`
    let mut tot_force = reb_vec3d {
        x: 0.,
        y: 0.,
        z: 0.,
    };

    if k2 != 0.0 {
        // Eggleton et. al 1998 quadrupole (equation 33)
        let quad_prefactor = mt * big_a / mu_ij;
        let omega_dot_d = Omega.x * dx + Omega.y * dy + Omega.z * dz;
        let omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;

        let t1 = 5. * omega_dot_d * omega_dot_d / (2. * (dr * dr * dr * dr * dr * dr * dr));
        let t2 = omega_squared / (2. * (dr * dr * dr * dr * dr));
        let t3 = omega_dot_d / (dr * dr * dr * dr * dr);
        let t4 = 3. * G * mt / (dr * dr * dr * dr * dr * dr * dr * dr);

        tot_force.x = quad_prefactor * ((t1 - t2 - t4) * dx - (t3 * Omega.x));
        tot_force.y = quad_prefactor * ((t1 - t2 - t4) * dy - (t3 * Omega.y));
        tot_force.z = quad_prefactor * ((t1 - t2 - t4) * dz - (t3 * Omega.z));

        if sigma != 0.0 {
            // Eggleton et. al 1998 tidal (equation 45)
            let d_dot_vel = dx * dvx + dy * dvy + dz * dvz;

            // first vector
            let vec1_x = 3. * d_dot_vel * dx;
            let vec1_y = 3. * d_dot_vel * dy;
            let vec1_z = 3. * d_dot_vel * dz;

            // h vector - EKH
            let hx = dy * dvz - dz * dvy;
            let hy = dz * dvx - dx * dvz;
            let hz = dx * dvy - dy * dvx;

            // h - r^2 Omega
            let comp_2_x = hx - d2 * Omega.x;
            let comp_2_y = hy - d2 * Omega.y;
            let comp_2_z = hz - d2 * Omega.z;

            // second vector
            let vec2_x = comp_2_y * dz - comp_2_z * dy;
            let vec2_y = comp_2_z * dx - comp_2_x * dz;
            let vec2_z = comp_2_x * dy - comp_2_y * dx;

            let prefactor =
                (-9. * sigma * mt * mt * big_a * big_a) / (2. * mu_ij * (d2 * d2 * d2 * d2 * d2));

            tot_force.x += prefactor * (vec1_x + vec2_x);
            tot_force.y += prefactor * (vec1_y + vec2_y);
            tot_force.z += prefactor * (vec1_z + vec2_z);
        }
    }

    tot_force
}

/// C: `static void rebx_spin_orbit_accelerations(struct reb_particle* source,
/// struct reb_particle* target, const double G, const double k2,
/// const double sigma, const struct reb_vec3d Omega)`.
///
/// The two C pointers become indices into the one particle array, so both
/// bodies can be written through safe Rust. The accelerations are written
/// in the C's order: target first (subtracted), then source (added).
fn rebx_spin_orbit_accelerations(
    particles: &mut [reb_particle],
    source: usize,
    target: usize,
    G: f64,
    k2: f64,
    sigma: f64,
    Omega: reb_vec3d,
) {
    // Input params all associated with source
    let ms = particles[source].m;
    let mt = particles[target].m;
    let mtot = ms + mt;

    // check if ODE is set here
    // (the helper reads positions/velocities/masses only, none of which is
    // written below, so the two copies are the same bytes the C dereferences)
    let source_p = particles[source];
    let target_p = particles[target];
    let tot_force =
        rebx_calculate_spin_orbit_accelerations(&source_p, &target_p, G, k2, sigma, Omega);

    particles[target].ax -= (ms / mtot) * tot_force.x;
    particles[target].ay -= (ms / mtot) * tot_force.y;
    particles[target].az -= (ms / mtot) * tot_force.z;

    particles[source].ax += (mt / mtot) * tot_force.x;
    particles[source].ay += (mt / mtot) * tot_force.y;
    particles[source].az += (mt / mtot) * tot_force.z;
}

/// Collects the spin vectors of every particle that has both a moment of
/// inertia `I` and a spin axis `Omega` set, flattened as
/// `[Omega_x, Omega_y, Omega_z, ...]` in particle order — i.e. exactly the
/// layout the C's `rebx_spin_sync_pre` writes into `ode->y`.
///
/// Factored out because `rebx_spin_initialize_ode` needs the very same
/// traversal (see the note on `rebx_spin_sync_pre`).
fn rebx_spin_collect(sim: &reb_simulation, rebx: &rebx_extras) -> Vec<f64> {
    let mut spins: Vec<f64> = Vec::new();
    let N = sim.N;
    for i in 0..N {
        let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        if let (Some(_I), Some(Omega)) = (I, Omega) {
            spins.push(Omega.x);
            spins.push(Omega.y);
            spins.push(Omega.z);
        }
    }
    spins
}

/// C: `static void rebx_spin_derivatives(struct reb_ode* const ode,
/// double* const yDot, const double* const y, const double t)`.
///
/// The C reaches REBOUNDx through `ode->ref == sim` and then `sim->extras`;
/// here the simulation arrives as an argument and the extras are taken out
/// of it (and always put back) by [`crate::core::rebx_with`].
///
/// `t` is unused in the C body as well.
fn rebx_spin_derivatives(
    sim: &mut reb_simulation,
    ode: &mut reb_ode,
    yDot: &mut [f64],
    y: &[f64],
    _t: f64,
) {
    // Read out before the extras are taken: the C reads `ode->length` at
    // the end of the function, but nothing in between can change it.
    let ode_length = ode.length;

    let _ = rebx_with(sim, |sim, rebx| {
        let mut Nspins: usize = 0; // C: unsigned int
        let N = sim.N;
        for i in 0..N {
            let pi = sim.particles[i]; // target particle
            let k2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "k2"); // This is slow
            let tau = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau");
            let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");

            // Particle MUST have k2 and moment of inertia to feel effects
            if let (Some(k2), Some(I)) = (k2, I) {
                // Tidal dissipation off by default. Check for non-zero tau here.
                let mut sigma_in = 0.0;
                if let Some(tau) = tau {
                    sigma_in = 4. * tau * sim.G / (3. * pi.r * pi.r * pi.r * pi.r * pi.r * k2);
                }
                // Set initial spin accelerations to 0
                yDot[3 * Nspins] = 0.;
                yDot[3 * Nspins + 1] = 0.;
                yDot[3 * Nspins + 2] = 0.;

                let Omega = reb_vec3d {
                    x: y[3 * Nspins],
                    y: y[3 * Nspins + 1],
                    z: y[3 * Nspins + 2],
                };
                for j in 0..N {
                    if i != j {
                        let pj = sim.particles[j];

                        let mi = pi.m;
                        let mj = pj.m;
                        // C declares `double mu_ij;` here; it is assigned and
                        // read only in the non-test-particle branch below, so
                        // it lives in that branch here.

                        if mj == 0. {
                            continue;
                        }

                        let I_specific = if mi == 0. {
                            // If test particle, assume I = specific moment of inertia
                            I
                        } else {
                            let mu_ij = (mi * mj) / (mi + mj);
                            I / mu_ij
                        };

                        // di - dj
                        let dx = pi.x - pj.x;
                        let dy = pi.y - pj.y;
                        let dz = pi.z - pj.z;

                        let tf = rebx_calculate_spin_orbit_accelerations(
                            &pi, &pj, sim.G, k2, sigma_in, Omega,
                        );
                        // Eggleton et. al 1998 spin EoM (equation 36)
                        yDot[3 * Nspins] += (dy * tf.z - dz * tf.y) / (-I_specific);
                        yDot[3 * Nspins + 1] += (dz * tf.x - dx * tf.z) / (-I_specific);
                        yDot[3 * Nspins + 2] += (dx * tf.y - dy * tf.x) / (-I_specific);
                    }
                }
                Nspins += 1;
            }
        }
        if ode_length != Nspins * 3 {
            reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
            std::process::exit(1); // C: exit(1)
        }
    });
}

/// C: `static void rebx_spin_sync_pre(struct reb_ode* const ode,
/// const double* const y0)`.
///
/// Copies each spin-tracked particle's `Omega` parameter into the ODE state
/// at the start of every BS step. The C writes straight into `ode->y` (its
/// `y0` argument is unused).
///
/// **Deviation, reported to the caller.** `rebound_rs`'s BS driver hands the
/// pre-timestep callback the state as an *immutable* slice and moves
/// `ode.y` out of the ODE for the duration of the call, restoring it
/// afterwards — so a store into `ode.y` made here is discarded. The
/// assignment below is the faithful translation and becomes effective the
/// moment `reb_ode_prepost_fn` takes `&mut [f64]` (the shape the C has).
/// Until then the initial load is performed once by
/// [`rebx_spin_initialize_ode`], which seeds `ode.y` with exactly these
/// values at creation time. That reproduces the C bit-for-bit for the
/// normal usage pattern (set `Omega`, call `rebx_spin_initialize_ode`,
/// integrate), because after the first successful step `rebx_spin_sync_post`
/// has already written `ode->y` back into the `Omega` parameters and this
/// function is a no-op. It diverges only if `Omega` is changed from outside
/// mid-integration without re-initializing the ODE.
///
/// The C indexes `ode->y` as it goes and only afterwards checks that the
/// length matches, which is a buffer overrun when it does not; the Vec built
/// here cannot overrun, and the same check follows.
fn rebx_spin_sync_pre(sim: &mut reb_simulation, ode: &mut reb_ode, _y0: &[f64]) {
    let ode_length = ode.length;

    let spins = rebx_with(sim, |sim, rebx| {
        let mut Nspins: usize = 0; // C: unsigned int
        let N = sim.N;
        let mut spins: Vec<f64> = Vec::new();
        for i in 0..N {
            let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");
            let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
            if let (Some(_I), Some(Omega)) = (I, Omega) {
                // C re-reads the same "Omega" parameter here into a shadowing
                // local; it is the identical value.
                spins.push(Omega.x);
                spins.push(Omega.y);
                spins.push(Omega.z);
                Nspins += 1;
            }
        }

        if ode_length != Nspins * 3 {
            reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
            std::process::exit(1); // C: exit(1)
        }
        spins
    });

    if let Some(spins) = spins {
        ode.y = spins;
    }
}

/// C: `static void rebx_spin_sync_post(struct reb_ode* const ode,
/// const double* const y0)`.
///
/// Writes the integrated spin state back into each tracked particle's
/// `Omega` parameter after a successful BS step.
fn rebx_spin_sync_post(sim: &mut reb_simulation, ode: &mut reb_ode, y0: &[f64]) {
    let ode_length = ode.length;

    let _ = rebx_with(sim, |sim, rebx| {
        let mut Nspins: usize = 0; // C: unsigned int
        let N = sim.N;
        for i in 0..N {
            let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");
            let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
            if let (Some(_I), Some(_Omega)) = (I, Omega) {
                rebx_set_param_vec3d(
                    rebx,
                    rebx_ap::particle(i),
                    "Omega",
                    reb_vec3d {
                        x: y0[3 * Nspins],
                        y: y0[3 * Nspins + 1],
                        z: y0[3 * Nspins + 2],
                    },
                );
                Nspins += 1;
            }
        }
        if ode_length != Nspins * 3 {
            reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
            std::process::exit(0); // C: exit(0) — note the C really does use 0 here
        }
    });
}

/// C: `void rebx_spin_initialize_ode(struct rebx_extras* const rebx,
/// struct rebx_force* const effect)`.
///
/// Must be called before integrating if the spin axes are to be evolved.
///
/// The C reaches the simulation through `rebx->sim`, which has no
/// counterpart here, so the simulation is passed explicitly (**sim first,
/// then rebx**, like its siblings in `rebxtools`) and the effect is named by
/// its index into `rebx_extras::allocated_forces`. The C's
/// `spin_ode->ref = sim` likewise has no counterpart: the callbacks receive
/// the simulation directly.
///
/// Call it as the C examples do, from inside [`crate::core::rebx_with`]:
///
/// ```ignore
/// rebx_with(&mut sim, |sim, rebx| {
///     rebx_spin_initialize_ode(sim, rebx, effect);
/// });
/// ```
///
/// The C identifies a previously-registered spin ODE by comparing
/// `ode->derivatives` against `rebx_spin_derivatives`; the same comparison
/// is made here with `std::ptr::fn_addr_eq`.
pub fn rebx_spin_initialize_ode(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
) {
    let mut Nspins: usize = 0; // C: unsigned int
    let N = sim.N;
    for i in 0..N {
        // Only track spin if particle has moment of inertia and valid spin axis set
        let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        if let (Some(_I), Some(_Omega)) = (I, Omega) {
            Nspins += 1;
        }
    }

    // Search for previous spin ode.
    // The C is `for (i=0; i<N_odes; i++){ if (..){ reb_ode_free(..); i--; } }`
    // — freeing shifts the array down and the `i--`/`i++` pair leaves the
    // cursor on the element that moved into the freed slot. The loop below
    // does exactly that by not advancing after a removal.
    let spin_derivatives: reb_ode_derivatives_fn = rebx_spin_derivatives;
    let mut i = 0;
    while i < sim.odes.len() {
        let is_spin_ode = match sim.odes[i].derivatives {
            Some(f) => std::ptr::fn_addr_eq(f, spin_derivatives),
            None => false,
        };
        if is_spin_ode {
            let id = sim.odes[i].id;
            reb_ode_free(sim, id);
        } else {
            i += 1;
        }
    }

    if Nspins > 0 {
        let spin_ode = reb_ode_create(sim, Nspins * 3);
        // C: spin_ode->ref = sim;  (no counterpart — see the module docs)
        if let Some(ode) = sim.odes.iter_mut().find(|o| o.id == spin_ode) {
            ode.derivatives = Some(rebx_spin_derivatives);
            ode.pre_timestep = Some(rebx_spin_sync_pre);
            ode.post_timestep = Some(rebx_spin_sync_post);
        }

        // Deviation from the C, forced by `rebound_rs`'s pre-timestep
        // callback signature: seed the freshly created state with the
        // spin vectors that `rebx_spin_sync_pre` would have loaded on the
        // first step. See the note on `rebx_spin_sync_pre`. The C leaves
        // `ode->y` zeroed here and lets `rebx_spin_sync_pre` fill it, so
        // the values integrated from are identical.
        let spins = rebx_spin_collect(sim, rebx);
        if let Some(ode) = sim.odes.iter_mut().find(|o| o.id == spin_ode) {
            ode.y = spins;
        }

        rebx_set_param_ode(rebx, rebx_ap::force(force_idx), "ode", spin_ode);
    }
}

/// C: `void rebx_tides_spin(struct reb_simulation* const sim,
/// struct rebx_force* const effect, struct reb_particle* const particles,
/// const int N)`.
///
/// The C `effect` pointer is unused in the body (this effect has no effect
/// parameters), so `_force_idx` is likewise unused here.
pub fn rebx_tides_spin(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    _force_idx: usize,
    N: usize,
) {
    let G = sim.G;

    // check if ODE is initialized
    // C: `struct reb_ode** ode = sim->odes; if (ode == NULL) ...` — the C
    // array pointer is NULL until the first `reb_ode_create`, which here is
    // an empty `sim.odes`.
    if sim.odes.is_empty() {
        reb_simulation_warning(
            sim,
            "Spin axes are not being evolved. Call rebx_spin_initialize_ode to evolve\n",
        );
    }

    for i in 0..N {
        // Particle must have a k2 set, otherwise we treat this body as a point particle
        let k2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "k2");
        let tau = rebx_get_param_double(rebx, rebx_ap::particle(i), "tau");
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");

        // Particle needs all three spin components and k2 to feel additional forces
        if let (Some(Omega), Some(k2)) = (Omega, k2) {
            // Tidal dissipation off by default. Check for non-zero tau here.
            let mut sigma_in = 0.0;
            if let Some(tau) = tau {
                let source_r = sim.particles[i].r;
                sigma_in = 4. * tau * sim.G
                    / (3. * source_r * source_r * source_r * source_r * source_r * k2);
            }

            for j in 0..N {
                if i == j {
                    continue;
                }
                // j raises tides on i
                if sim.particles[i].m == 0. || sim.particles[j].m == 0. {
                    continue;
                }

                rebx_spin_orbit_accelerations(&mut sim.particles, i, j, G, k2, sigma_in, Omega);
            }
        }
    }
}

/// Calculate potential of conservative piece of interaction between a point
/// mass target and a source with a tidally and rotationally induced
/// quadrupole. Equation 31 in Eggleton et. al (1998).
///
/// C: `static double rebx_calculate_spin_potential(struct reb_particle* source,
/// struct reb_particle* target, const double G, const double k2,
/// const struct reb_vec3d Omega)`.
fn rebx_calculate_spin_potential(
    source: &reb_particle,
    target: &reb_particle,
    G: f64,
    k2: f64,
    Omega: reb_vec3d,
) -> f64 {
    let Rs = source.r;
    let mt = target.m;

    let big_a = k2 * (Rs * Rs * Rs * Rs * Rs);

    // distance vector FROM j TO i
    let dx = source.x - target.x;
    let dy = source.y - target.y;
    let dz = source.z - target.z;
    let d2 = dx * dx + dy * dy + dz * dz;
    let dr = d2.sqrt();

    let omega_dot_d = Omega.x * dx + Omega.y * dy + Omega.z * dz;
    let omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;

    let t1 = -mt * big_a * omega_dot_d * omega_dot_d / (2. * d2 * d2 * dr);
    let t2 = mt * big_a * omega_squared / (6. * d2 * dr);
    let t3 = G * mt * mt * big_a / (2. * d2 * d2 * d2);

    -(t1 + t2 + t3)
}

/// C: `double rebx_tides_spin_energy(struct rebx_extras* const rebx)`.
///
/// The C reaches the simulation through the `rebx->sim` back-pointer and
/// bails out with `rebx_error` when it is NULL. `rebx_extras` carries no
/// back-pointer here, so the simulation is passed explicitly (**sim first,
/// then rebx**, like its siblings) and that NULL branch has no counterpart.
pub fn rebx_tides_spin_energy(sim: &reb_simulation, rebx: &rebx_extras) -> f64 {
    let N = sim.N;
    let particles = &sim.particles;
    let G = sim.G;
    let mut E = 0.;

    for i in 0..N {
        // Particle must have a k2, radius and mass set, otherwise we treat this body as a point particle
        let k2 = rebx_get_param_double(rebx, rebx_ap::particle(i), "k2");
        let Omegaptr = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        let k2 = match k2 {
            None => continue,
            Some(k2) => k2,
        };
        if particles[i].m == 0. || particles[i].r == 0. {
            continue;
        }
        // C: `struct reb_vec3d Omega = {0};`
        let mut Omega = reb_vec3d {
            x: 0.,
            y: 0.,
            z: 0.,
        };
        if let Some(Omegaptr) = Omegaptr {
            Omega = Omegaptr;
        }
        let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");
        if let Some(I) = I {
            let omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;
            E += 0.5 * I * omega_squared;
        }
        for j in 0..N {
            if i == j {
                continue;
            }
            // planet raising the tides on the star
            if particles[j].m > 0. {
                E += rebx_calculate_spin_potential(&particles[i], &particles[j], G, k2, Omega);
            }
        }
    }

    E
}
