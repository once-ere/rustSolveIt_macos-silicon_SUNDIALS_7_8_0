//! boundary.rs — boundary conditions and ghost boxes (from boundary.c).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::particle::reb_simulation_remove_particle;
use crate::tools::reb_simulation_energy;
use crate::types::*;

/// boundary.c `reb_boundary_check`.
pub fn reb_boundary_check(r: &mut reb_simulation) {
    let mut N = r.N;
    let boxsize = reb_vec3d {
        x: r.root_size * (r.N_root_x as f64),
        y: r.root_size * (r.N_root_y as f64),
        z: r.root_size * (r.N_root_z as f64),
    };
    match r.boundary {
        REB_BOUNDARY::OPEN => {
            let mut i = 0usize;
            while i < N {
                let p = r.particles[i];
                let mut removep = false;
                if p.x > boxsize.x / 2. {
                    removep = true;
                }
                if p.x < -boxsize.x / 2. {
                    removep = true;
                }
                if p.y > boxsize.y / 2. {
                    removep = true;
                }
                if p.y < -boxsize.y / 2. {
                    removep = true;
                }
                if p.z > boxsize.z / 2. {
                    removep = true;
                }
                if p.z < -boxsize.z / 2. {
                    removep = true;
                }
                if removep {
                    if r.track_energy_offset != 0 {
                        let Ei = reb_simulation_energy(r);
                        reb_simulation_remove_particle(r, i);
                        r.energy_offset += Ei - reb_simulation_energy(r);
                    } else {
                        reb_simulation_remove_particle(r, i);
                    }
                    // i stays: need to recheck the particle that replaced
                    // the removed one (C decrements then re-increments).
                    N -= 1;
                } else {
                    i += 1;
                }
            }
        }
        REB_BOUNDARY::SHEAR => {
            // The offset of ghostcell is time dependent.
            let OMEGA = r.OMEGA;
            let offsetp1 =
                -(-1.5 * OMEGA * boxsize.x * r.t + boxsize.y / 2.) % boxsize.y - boxsize.y / 2.;
            let offsetm1 =
                -(1.5 * OMEGA * boxsize.x * r.t - boxsize.y / 2.) % boxsize.y + boxsize.y / 2.;
            for i in 0..N {
                // Radial
                while r.particles[i].x > boxsize.x / 2. {
                    r.particles[i].x -= boxsize.x;
                    r.particles[i].y += offsetp1;
                    r.particles[i].vy += 3. / 2. * OMEGA * boxsize.x;
                }
                while r.particles[i].x < -boxsize.x / 2. {
                    r.particles[i].x += boxsize.x;
                    r.particles[i].y += offsetm1;
                    r.particles[i].vy -= 3. / 2. * OMEGA * boxsize.x;
                }
                // Azimuthal
                while r.particles[i].y > boxsize.y / 2. {
                    r.particles[i].y -= boxsize.y;
                }
                while r.particles[i].y < -boxsize.y / 2. {
                    r.particles[i].y += boxsize.y;
                }
                // Vertical (there should be no boundary, but periodic
                // makes life easier)
                while r.particles[i].z > boxsize.z / 2. {
                    r.particles[i].z -= boxsize.z;
                }
                while r.particles[i].z < -boxsize.z / 2. {
                    r.particles[i].z += boxsize.z;
                }
            }
        }
        REB_BOUNDARY::PERIODIC => {
            for i in 0..N {
                while r.particles[i].x > boxsize.x / 2. {
                    r.particles[i].x -= boxsize.x;
                }
                while r.particles[i].x < -boxsize.x / 2. {
                    r.particles[i].x += boxsize.x;
                }
                while r.particles[i].y > boxsize.y / 2. {
                    r.particles[i].y -= boxsize.y;
                }
                while r.particles[i].y < -boxsize.y / 2. {
                    r.particles[i].y += boxsize.y;
                }
                while r.particles[i].z > boxsize.z / 2. {
                    r.particles[i].z -= boxsize.z;
                }
                while r.particles[i].z < -boxsize.z / 2. {
                    r.particles[i].z += boxsize.z;
                }
            }
        }
        REB_BOUNDARY::NONE => {}
    }
}

/// boundary.c `reb_boundary_get_ghostbox`.
pub fn reb_boundary_get_ghostbox(r: &reb_simulation, i: i32, j: i32, k: i32) -> reb_vec6d {
    let boxsize = reb_vec3d {
        x: r.root_size * (r.N_root_x as f64),
        y: r.root_size * (r.N_root_y as f64),
        z: r.root_size * (r.N_root_z as f64),
    };
    match r.boundary {
        REB_BOUNDARY::OPEN | REB_BOUNDARY::PERIODIC => reb_vec6d {
            x: boxsize.x * (i as f64),
            y: boxsize.y * (j as f64),
            z: boxsize.z * (k as f64),
            vx: 0.,
            vy: 0.,
            vz: 0.,
        },
        REB_BOUNDARY::SHEAR => {
            let OMEGA = r.OMEGA;
            let mut gb = reb_vec6d::default();
            // Ghostboxes have a finite velocity.
            gb.vx = 0.;
            gb.vy = -1.5 * (i as f64) * OMEGA * boxsize.x;
            gb.vz = 0.;
            // The shift in the y direction is time dependent.
            let shift;
            if i == 0 {
                shift = -(gb.vy * r.t) % boxsize.y;
            } else if i > 0 {
                shift = -(gb.vy * r.t - boxsize.y / 2.) % boxsize.y - boxsize.y / 2.;
            } else {
                shift = -(gb.vy * r.t + boxsize.y / 2.) % boxsize.y + boxsize.y / 2.;
            }
            gb.x = boxsize.x * (i as f64);
            gb.y = boxsize.y * (j as f64) - shift;
            gb.z = boxsize.z * (k as f64);
            gb
        }
        REB_BOUNDARY::NONE => reb_vec6d::default(),
    }
}

/// boundary.c `reb_boundary_particle_is_in_box`.
pub fn reb_boundary_particle_is_in_box(r: &reb_simulation, p: reb_particle) -> bool {
    match r.boundary {
        REB_BOUNDARY::OPEN | REB_BOUNDARY::SHEAR | REB_BOUNDARY::PERIODIC => {
            let boxsize = reb_vec3d {
                x: r.root_size * (r.N_root_x as f64),
                y: r.root_size * (r.N_root_y as f64),
                z: r.root_size * (r.N_root_z as f64),
            };
            if p.x > boxsize.x / 2. {
                return false;
            }
            if p.x < -boxsize.x / 2. {
                return false;
            }
            if p.y > boxsize.y / 2. {
                return false;
            }
            if p.y < -boxsize.y / 2. {
                return false;
            }
            if p.z > boxsize.z / 2. {
                return false;
            }
            if p.z < -boxsize.z / 2. {
                return false;
            }
            true
        }
        REB_BOUNDARY::NONE => true,
    }
}
