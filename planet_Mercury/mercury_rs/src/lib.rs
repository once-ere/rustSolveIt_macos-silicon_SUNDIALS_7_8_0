#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

//! mercury_rs — the historical tidal despinning and 3:2 spin-orbit resonant
//! capture of planet Mercury, as a two-body Sun–Mercury system.
//!
//! Model: the Sun is a point mass; Mercury is an extended, deformable,
//! slightly triaxial body. Five state variables y = [a, e, M, theta, Omega]
//! evolve under Kepler's law, the Hut (1981) constant-time-lag tidal torque
//! (with its back-reaction on a and e), and the permanent-triaxiality
//! ("handle") torque of Goldreich & Peale (1966).
//!
//! The ONLY integrator is the vendored pure-Rust SUNDIALS 7.8.0 CVODE
//! (BDF + Newton + dense linear solver) from the rustSolveIt engine.

pub mod params;
pub mod kepler;
pub mod hut;
pub mod rhs;
pub mod driver;
pub mod output;
pub mod test2;
