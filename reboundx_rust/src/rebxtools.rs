//! rebxtools.rs — translation of REBOUNDx rebxtools.c
//! Helper functions shared by the REBOUNDx effects: centre-of-mass
//! bookkeeping for the three supported coordinate systems, Jacobi
//! masses, energy/angular-momentum diagnostics and whole-simulation
//! rotations that also carry the spin vectors along.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # Header (rebxtools.h)
//!
//! ```text
//! @file    rebxtools.h
//! @brief   Helper functions for reboundx
//! @author  Dan Tamayo <tamayo.daniel@gmail.com>
//! ```
//!
//! These helpers take no parameters of their own. The two parameters
//! they read off other objects are:
//!
//! ```text
//! ============================ =========== ==================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ==================================================================
//! primary (int)                No          Flag on the particle used as the reference body when an effect's
//!                                          `coordinates` is REBX_COORDINATES_PARTICLE. Only its presence is
//!                                          tested; the value is ignored.
//! Omega (reb_vec3d)            No          Spin angular rotation frequency vector of a particle. Read by
//!                                          rebx_tools_spin_angular_momentum, rebx_tools_spin_energy and
//!                                          rebx_simulation_irotate.
//! I (double)                   No          Moment of inertia of a particle. Read by
//!                                          rebx_tools_spin_angular_momentum and rebx_tools_spin_energy.
//! ============================ =========== ==================================================================
//! ```
//!
//! # Deviations from the C, all mechanical
//!
//! * The C passes `struct rebx_force*` / `struct rebx_operator*` and a
//!   `struct reb_particle*` array. Here a force/operator is an index
//!   into `rebx_extras::allocated_forces` / `allocated_operators`, the
//!   REBOUNDx state is passed explicitly as `&mut rebx_extras`, and the
//!   particles are `sim.particles` (the same array the C force functions
//!   are handed).
//! * The C callbacks reach a particle's parameter list through
//!   `p->ap`. `reb_particle` owns no list in this translation, so the
//!   callbacks additionally receive the particle's index `p_index`,
//!   which is exactly what `rebx_ap::particle(p_index)` names. The
//!   particle and the reference body are passed by value because
//!   `reb_particle` is `Copy` and the callbacks only read them.
//! * `enum REBX_COORDINATES` (declared in reboundx.h) is an `i32` alias
//!   with the C's three constants, mirroring how `rebx_integrator` is
//!   handled in `types`. It is stored as a `REBX_TYPE_INT` parameter
//!   named `"coordinates"`, so keeping it an integer both makes the
//!   assignment `coordinates = *ptr` direct and keeps the C's
//!   `default:` "Coordinates not supported" branch reachable — a real
//!   Rust enum would make that error path dead code.

use crate::core::{rebx_get_param_double, rebx_get_param_int, rebx_get_param_vec3d,
                  rebx_set_param_vec3d};
use crate::types::{rebx_ap, rebx_extras};
use rebound_rs::{
    reb_particle, reb_rotation, reb_simulation, reb_simulation_com, reb_simulation_error,
    reb_simulation_irotate, reb_vec3d, reb_vec3d_irotate,
};

/// reboundx.h `enum REBX_COORDINATES`. Carried as an `i32` because the
/// value travels through a `REBX_TYPE_INT` parameter named
/// `"coordinates"` (see the module docs).
pub type REBX_COORDINATES = i32;
/// Jacobi coordinates (default).
pub const REBX_COORDINATES_JACOBI: REBX_COORDINATES = 0;
/// Coordinates referenced to pos/vel of system's center of mass.
pub const REBX_COORDINATES_BARYCENTRIC: REBX_COORDINATES = 1;
/// Coordinates relative to pos/vel of a particular particle.
pub const REBX_COORDINATES_PARTICLE: REBX_COORDINATES = 2;

/// C: `struct reb_vec3d (*calculate_force)(struct reb_simulation* const
/// sim, struct rebx_force* const force, struct reb_particle* p,
/// struct reb_particle* source)`, the per-particle acceleration
/// callback handed to [`rebx_com_force`].
///
/// `force_idx` indexes `rebx_extras::allocated_forces` (C: the `force`
/// pointer) and `p_index` names the particle's parameter list
/// (C: `p->ap`).
pub type rebx_calculate_force_fn = fn(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    p_index: usize,
    p: reb_particle,
    source: reb_particle,
) -> reb_vec3d;

/// C: `struct reb_particle (*calculate_step)(struct reb_simulation*
/// const sim, struct rebx_operator* const operator, struct reb_particle*
/// p, struct reb_particle* source, const double dt)`, the per-particle
/// update callback handed to [`rebx_tools_com_ptm`].
pub type rebx_calculate_step_fn = fn(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    operator_idx: usize,
    p_index: usize,
    p: reb_particle,
    source: reb_particle,
    dt: f64,
) -> reb_particle;

/// rebxtools.c `rebx_get_com_without_particle`. Returns the
/// centre-of-mass `com` with particle `p` removed from it.
pub fn rebx_get_com_without_particle(com: reb_particle, p: reb_particle) -> reb_particle {
    let mut com = com;
    com.x = com.x * com.m - p.x * p.m;
    com.y = com.y * com.m - p.y * p.m;
    com.z = com.z * com.m - p.z * p.m;
    com.vx = com.vx * com.m - p.vx * p.m;
    com.vy = com.vy * com.m - p.vy * p.m;
    com.vz = com.vz * com.m - p.vz * p.m;
    com.ax = com.ax * com.m - p.ax * p.m;
    com.ay = com.ay * com.m - p.ay * p.m;
    com.az = com.az * com.m - p.az * p.m;
    com.m -= p.m;

    if com.m > 0. {
        com.x /= com.m;
        com.y /= com.m;
        com.z /= com.m;
        com.vx /= com.m;
        com.vy /= com.m;
        com.vz /= com.m;
        com.ax /= com.m;
        com.ay /= com.m;
        com.az /= com.m;
    }
    com
}

/// rebxtools.c `rebx_particle_minus` (C: `static inline`).
fn rebx_particle_minus(p1: reb_particle, p2: reb_particle) -> reb_particle {
    let mut p = reb_particle::default(); // C: struct reb_particle p = {0};
    p.m = p1.m - p2.m;
    p.x = p1.x - p2.x;
    p.y = p1.y - p2.y;
    p.z = p1.z - p2.z;
    p.vx = p1.vx - p2.vx;
    p.vy = p1.vy - p2.vy;
    p.vz = p1.vz - p2.vz;
    p.ax = p1.ax - p2.ax;
    p.ay = p1.ay - p2.ay;
    p.az = p1.az - p2.az;
    p
}

/// rebxtools.c `rebx_calculate_jacobi_masses`.
///
/// Jacobi masses are the reduced mass of a particle with the interior
/// masses; `m_j[0]` ends up holding the total mass.
pub fn rebx_calculate_jacobi_masses(ps: &[reb_particle], m_j: &mut [f64], N: usize) {
    let mut eta = ps[0].m;
    for i in 1..N {
        // jacobi masses are reduced mass of particle with interior masses
        m_j[i] = ps[i].m * eta;
        eta += ps[i].m;
        m_j[i] /= eta;
    }
    m_j[0] = eta;
}

/// rebxtools.c `rebx_Edot`. Rate of change of the kinetic energy given
/// the accelerations currently stored on the particles.
pub fn rebx_Edot(ps: &[reb_particle], N: usize) -> f64 {
    let mut Edot = 0.;
    for i in 0..N {
        Edot += ps[i].m * (ps[i].ax * ps[i].vx + ps[i].ay * ps[i].vy + ps[i].az * ps[i].vz);
    }
    Edot
}

/// rebxtools.c `rebx_com_force`.
///
/// Applies `calculate_force` to every particle (excluding the reference
/// body) in the requested coordinate system and adds the matching back
/// reaction, so that momentum is conserved.
///
/// Argument order follows the C — `(sim, force, coordinates,
/// back_reactions_inclusive, reference_name, calculate_force,
/// particles, N)` — with the C's `force` pointer replaced by
/// `(rebx, force_idx)` and `particles` dropped, since it is
/// `sim.particles`.
pub fn rebx_com_force(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    force_idx: usize,
    coordinates: REBX_COORDINATES,
    back_reactions_inclusive: i32,
    reference_name: &str,
    calculate_force: rebx_calculate_force_fn,
    N: usize,
) {
    // Start with full com for jacobi and barycentric coordinates.
    let mut com = reb_simulation_com(sim);

    let mut refindex: i64 = -1;
    if coordinates == REBX_COORDINATES_JACOBI {
        // There is no jacobi coordinate for the 0th particle, so set
        // refindex to skip it in loop below.
        refindex = 0;
    } else if coordinates == REBX_COORDINATES_PARTICLE {
        for i in 0..N {
            let reference = rebx_get_param_int(rebx, rebx_ap::particle(i), reference_name);
            if reference.is_some() {
                com = sim.particles[i];
                refindex = i as i64;
                break;
            }
            if i == N - 1 {
                let str = format!(
                    "Coordinates set to REBX_COORDINATES_PARTICLE, but {} param was not found in any particle.  Need to set parameter.\n",
                    reference_name
                );
                reb_simulation_error(sim, &str);
            }
        }
    }

    // Run through backwards so each iteration does not depend on
    // previous ones in Jacobi coordinates.
    for i in (0..N).rev() {
        if i as i64 == refindex {
            continue;
        }
        if coordinates == REBX_COORDINATES_JACOBI {
            com = rebx_get_com_without_particle(com, sim.particles[i]);
        }

        let p = sim.particles[i];
        let a = calculate_force(sim, rebx, force_idx, i, p, com);
        sim.particles[i].ax += a.x;
        sim.particles[i].ay += a.y;
        sim.particles[i].az += a.z;

        let massratio;
        if coordinates == REBX_COORDINATES_BARYCENTRIC {
            massratio = sim.particles[i].m / com.m;
            for j in 0..N {
                sim.particles[j].ax -= massratio * a.x;
                sim.particles[j].ay -= massratio * a.y;
                sim.particles[j].az -= massratio * a.z;
            }
        } else if coordinates == REBX_COORDINATES_JACOBI {
            if back_reactions_inclusive != 0 {
                massratio = sim.particles[i].m / (com.m + sim.particles[i].m);
            } else {
                massratio = sim.particles[i].m / com.m;
            }
            // stop at j=i if inclusive, at i-1 if not
            let jmax: i64 = i as i64 + back_reactions_inclusive as i64;
            let mut j: i64 = 0;
            while j < jmax {
                let ju = j as usize;
                sim.particles[ju].ax -= massratio * a.x;
                sim.particles[ju].ay -= massratio * a.y;
                sim.particles[ju].az -= massratio * a.z;
                j += 1;
            }
        } else if coordinates == REBX_COORDINATES_PARTICLE {
            if back_reactions_inclusive != 0 {
                massratio = sim.particles[i].m / (com.m + sim.particles[i].m);
                sim.particles[i].ax -= massratio * a.x;
                sim.particles[i].ay -= massratio * a.y;
                sim.particles[i].az -= massratio * a.z;
            } else {
                massratio = sim.particles[i].m / com.m;
            }
            // If no particle carried `reference_name`, refindex is still
            // -1 here. The C indexes particles[-1] (undefined behaviour);
            // safe Rust panics on the out-of-range index instead.
            let refi = refindex as usize;
            sim.particles[refi].ax -= massratio * a.x;
            sim.particles[refi].ay -= massratio * a.y;
            sim.particles[refi].az -= massratio * a.z;
        } else {
            reb_simulation_error(sim, "Coordinates not supported in REBOUNDx.\n");
        }
    }
}

/// rebxtools.c `rebx_subtract_posvel` (C: `static inline`).
fn rebx_subtract_posvel(p: &mut reb_particle, diff: &reb_particle, massratio: f64) {
    p.x -= massratio * diff.x;
    p.y -= massratio * diff.y;
    p.z -= massratio * diff.z;
    p.vx -= massratio * diff.vx;
    p.vy -= massratio * diff.vy;
    p.vz -= massratio * diff.vz;
}

/// rebxtools.c `rebx_tools_com_ptm` — the operator counterpart of
/// [`rebx_com_force`].
///
/// Only accepts one reference particle if
/// `coordinates == REBX_COORDINATES_PARTICLE`. The `calculate_step`
/// function should check for the edge case where the particle and the
/// reference are the same (could happen e.g. with barycentric
/// coordinates with test particles and a single massive body).
///
/// Argument order follows the C, with the `operator` pointer replaced by
/// `(rebx, operator_idx)`.
pub fn rebx_tools_com_ptm(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
    operator_idx: usize,
    coordinates: REBX_COORDINATES,
    back_reactions_inclusive: i32,
    reference_name: &str,
    calculate_step: rebx_calculate_step_fn,
    dt: f64,
) {
    let N = sim.N;
    // Start with full com for jacobi and barycentric coordinates.
    let mut com = reb_simulation_com(sim);

    let mut refindex: i64 = -1;
    if coordinates == REBX_COORDINATES_JACOBI {
        // There is no jacobi coordinate for the 0th particle, so should
        // skip index 0
        refindex = 0;
    } else if coordinates == REBX_COORDINATES_PARTICLE {
        for i in 0..N {
            let reference = rebx_get_param_int(rebx, rebx_ap::particle(i), reference_name);
            if reference.is_some() {
                com = sim.particles[i];
                refindex = i as i64;
                break;
            }
            if i == N - 1 {
                let str = format!(
                    "Coordinates set to REBX_COORDINATES_PARTICLE, but {} param was not found in any particle.  Need to set parameter.\n",
                    reference_name
                );
                reb_simulation_error(sim, &str);
            }
        }
    }

    // Run through backwards so each iteration does not depend on
    // previous ones in Jacobi coordinates.
    for i in (0..N).rev() {
        if i as i64 == refindex {
            continue;
        }
        if coordinates == REBX_COORDINATES_JACOBI {
            com = rebx_get_com_without_particle(com, sim.particles[i]);
        }

        let p = sim.particles[i];
        let modified_particle = calculate_step(sim, rebx, operator_idx, i, p, com, dt);
        let diff = rebx_particle_minus(modified_particle, p);
        sim.particles[i].x = modified_particle.x;
        sim.particles[i].y = modified_particle.y;
        sim.particles[i].z = modified_particle.z;
        sim.particles[i].vx = modified_particle.vx;
        sim.particles[i].vy = modified_particle.vy;
        sim.particles[i].vz = modified_particle.vz;

        let massratio;
        if coordinates == REBX_COORDINATES_BARYCENTRIC {
            massratio = sim.particles[i].m / com.m;
            for j in 0..N {
                rebx_subtract_posvel(&mut sim.particles[j], &diff, massratio);
            }
        } else if coordinates == REBX_COORDINATES_JACOBI {
            if back_reactions_inclusive != 0 {
                massratio = sim.particles[i].m / (com.m + sim.particles[i].m);
            } else {
                massratio = sim.particles[i].m / com.m;
            }
            // stop at j=i if inclusive, at i-1 if not
            let jmax: i64 = i as i64 + back_reactions_inclusive as i64;
            let mut j: i64 = 0;
            while j < jmax {
                rebx_subtract_posvel(&mut sim.particles[j as usize], &diff, massratio);
                j += 1;
            }
        } else if coordinates == REBX_COORDINATES_PARTICLE {
            if back_reactions_inclusive != 0 {
                massratio = sim.particles[i].m / (com.m + sim.particles[i].m);
                rebx_subtract_posvel(&mut sim.particles[i], &diff, massratio);
            } else {
                massratio = sim.particles[i].m / com.m;
            }
            // See the note in rebx_com_force: refindex is -1 here if no
            // particle carried `reference_name`.
            let refi = refindex as usize;
            rebx_subtract_posvel(&mut sim.particles[refi], &diff, massratio);
        } else {
            reb_simulation_error(sim, "Coordinates not supported in REBOUNDx.\n");
        }
    }
}

/// rebxtools.c `rebx_tools_spin_angular_momentum`.
///
/// Calculates the spin angular momentum in the simulation of any bodies
/// with spin parameters set (moment of inertia `I` and angular rotation
/// frequency vector `Omega`).
///
/// The C reads the simulation through `rebx->sim`; there is no such
/// back-pointer here, so the simulation is passed explicitly. Both
/// arguments are taken by `&mut` because the per-particle parameters
/// ("Omega", "I") live in `rebx_extras`, in the same order the sibling
/// [`rebx_simulation_irotate`] needs them: **sim first, then rebx**.
pub fn rebx_tools_spin_angular_momentum(
    sim: &mut reb_simulation,
    rebx: &mut rebx_extras,
) -> reb_vec3d {
    // Add spin angular momentum of any particles with spin parameters set
    let N = sim.N;
    let mut L = reb_vec3d { x: 0., y: 0., z: 0. }; // C: struct reb_vec3d L = {0.};
    for i in 0..N {
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");

        if let (Some(Omega), Some(I)) = (Omega, I) {
            L.x += I * (Omega.x);
            L.y += I * (Omega.y);
            L.z += I * (Omega.z);
        }
    }
    L
}

/// rebxtools.c `rebx_tools_spin_energy`.
///
/// Calculates the spin energy in the simulation of any bodies with spin
/// parameters set (moment of inertia `I` and angular rotation frequency
/// vector `Omega`). Takes **sim first, then rebx**, like its siblings.
pub fn rebx_tools_spin_energy(sim: &mut reb_simulation, rebx: &mut rebx_extras) -> f64 {
    // Add spin energy of any particles with spin parameters set
    let N = sim.N;
    let mut E = 0.;
    for i in 0..N {
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        let I = rebx_get_param_double(rebx, rebx_ap::particle(i), "I");

        if let (Some(Omega), Some(I)) = (Omega, I) {
            E += 0.5 * I * ((Omega.x) * (Omega.x) + (Omega.y) * (Omega.y) + (Omega.z) * (Omega.z));
        }
    }
    E
}

/// rebxtools.c `rebx_simulation_irotate`.
///
/// Rotates every particle's position and velocity by `q` (through
/// REBOUND's `reb_simulation_irotate`) and then rotates the spin vector
/// `Omega` of every particle that has one, so that spin angular momentum
/// is carried along with the orbits.
///
/// Modified from celmech `nbody_simulation_utilities.py` to include spin
/// angular momentum.
///
/// Takes **sim first, then rebx**: the C reaches the simulation through
/// `rebx->sim`, which does not exist here, and the "Omega" parameters it
/// rewrites live in `rebx_extras`.
pub fn rebx_simulation_irotate(sim: &mut reb_simulation, rebx: &mut rebx_extras, q: reb_rotation) {
    reb_simulation_irotate(sim, q); // rotate all the orbits first
    for i in 0..sim.N {
        // Rotate spins
        let Omega = rebx_get_param_vec3d(rebx, rebx_ap::particle(i), "Omega");
        if let Some(mut Omega) = Omega {
            reb_vec3d_irotate(&mut Omega, q);
            rebx_set_param_vec3d(rebx, rebx_ap::particle(i), "Omega", Omega);
        }
    }
}
