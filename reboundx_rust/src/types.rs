//! types.rs — the REBOUNDx data structures (translated from
//! reboundx.h, linkedlist.h and core.h).
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx
//! 5.1.0 (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! # How the C pointer graph becomes safe Rust
//!
//! REBOUNDx is built on three C constructs that safe Rust cannot copy
//! literally. Each has one mechanical, behavior-preserving replacement,
//! used consistently throughout this crate:
//!
//! 1. **`struct rebx_node` linked lists** (`ap`, `additional_forces`,
//!    `allocated_forces`, ...). The C `rebx_add_node` *prepends*:
//!    `node->next = *head; *head = node;`  So a C list iterates in
//!    reverse insertion order, and that order decides the order in which
//!    accelerations are summed — which changes floating-point results.
//!    Here a list is a `Vec` whose **index 0 is the head**, and the
//!    add helper does `insert(0, ..)`. Iteration order is therefore
//!    identical to the C, element for element.
//!
//! 2. **`void* value` parameters with a type tag.** C stores a type
//!    enum plus an untyped pointer and casts on read. Here the tag and
//!    the payload are one `enum rebx_param_value`, so a wrong-type read
//!    is impossible rather than undefined.
//!
//! 3. **Pointers naming *which* list to act on** (`&p->ap`,
//!    `&force->ap`). Here `rebx_ap` names the same thing by index:
//!    `rebx_ap::particle(0)` is exactly C's `&sim->particles[0].ap`.
//!    Particle parameter lists live in `rebx_extras::particle_params`
//!    rather than inside `reb_particle`, because `reb_particle` is
//!    `Copy` in the REBOUND translation and cannot own a heap list.

use rebound_rs::{reb_orbit, reb_particle, reb_simulation, reb_vec3d};

/// reboundx.h `enum rebx_param_type`. Kept for the binary format and
/// for the registered-parameter type checks; the payload itself is
/// carried by `rebx_param_value`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_param_type {
    REBX_TYPE_NONE = 0,
    REBX_TYPE_DOUBLE = 1,
    REBX_TYPE_INT = 2,
    REBX_TYPE_POINTER = 3,
    REBX_TYPE_FORCE = 4,
    REBX_TYPE_UINT32 = 5,
    REBX_TYPE_ORBIT = 6,
    REBX_TYPE_ODE = 7,
    REBX_TYPE_VEC3D = 8,
    REBX_TYPE_STRING = 9,
}

/// The value of one parameter. The C keeps `void* value` plus the type
/// tag above; this enum fuses the two so the pair can never disagree.
#[derive(Clone, Debug, PartialEq)]
pub enum rebx_param_value {
    /// C: `value == NULL` (registered but never set).
    none,
    double(f64),
    int(i32),
    uint32(u32),
    vec3d(reb_vec3d),
    orbit(reb_orbit),
    /// C: `char*` owned by REBOUND's name list.
    string(String),
    /// C: `struct rebx_force*` — an index into `allocated_forces`.
    force(usize),
    /// C: `struct reb_ode*` — an id of an ODE in `reb_simulation::odes`.
    ode(usize),
    /// C: `REBX_TYPE_POINTER` payloads. Every use of this type in
    /// REBOUNDx 5.1.0 is an internally-allocated `struct reb_particle`
    /// buffer belonging to one of the REBOUNDx integrators
    /// (`im_ps_final`, `im_ps_prev`, `im_ps_avg`, `rk2_k2`, `rk4_k2`,
    /// `rk4_k3`) or the `particle` back-reference used by
    /// `track_min_distance`.
    particles(Vec<reb_particle>),
    /// C: `REBX_TYPE_POINTER` used as a particle index back-reference.
    particle_index(usize),
}

/// reboundx.h `struct rebx_param`.
#[derive(Clone, Debug, PartialEq)]
pub struct rebx_param {
    /// For searching lists and informative errors.
    pub name: String,
    /// Needed to cast value in the C; carried for type checks and I/O.
    pub type_: rebx_param_type,
    /// Parameter value (C: `void* value`).
    pub value: rebx_param_value,
}

/// Names *which* parameter list (C: the `struct rebx_node** apptr`
/// argument, always the address of some object's `ap` member).
///
/// `rebx_ap::particle(i)` is C's `&sim->particles[i].ap`,
/// `rebx_ap::force(i)` is C's `&force->ap`, and
/// `rebx_ap::operator_(i)` is C's `&operator->ap`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_ap {
    particle(usize),
    force(usize),
    operator_(usize),
}

/// reboundx.h `enum rebx_timing`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_timing {
    REBX_TIMING_PRE = -1,
    REBX_TIMING_POST = 1,
}

/// reboundx.h `enum rebx_force_type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_force_type {
    REBX_FORCE_NONE = 0,
    /// Force derivable from a position-dependent potential.
    REBX_FORCE_POS = 1,
    /// Velocity (or pos and vel) dependent force.
    REBX_FORCE_VEL = 2,
}

/// reboundx.h `enum rebx_operator_type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_operator_type {
    REBX_OPERATOR_NONE = 0,
    /// Operator that modifies x, v or m.
    REBX_OPERATOR_UPDATER = 1,
    /// Operator that leaves state unchanged; just records.
    REBX_OPERATOR_RECORDER = 2,
}

/// reboundx.h `enum rebx_integrator`.
pub type rebx_integrator = i32;
pub const REBX_INTEGRATOR_NONE: i32 = -1;
pub const REBX_INTEGRATOR_IMPLICIT_MIDPOINT: i32 = 0;
pub const REBX_INTEGRATOR_RK4: i32 = 1;
pub const REBX_INTEGRATOR_EULER: i32 = 2;
pub const REBX_INTEGRATOR_RK2: i32 = 3;

/// reboundx.h `enum rebx_interpolation_type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum rebx_interpolation_type {
    REBX_INTERPOLATION_NONE = 0,
    REBX_INTERPOLATION_SPLINE = 1,
}

/// C: `void (*update_accelerations)(struct reb_simulation* const sim,
/// struct rebx_force* const force, struct reb_particle* const
/// particles, const int N)`.
///
/// The C reaches the REBOUNDx state through `sim->extras` and its own
/// parameters through the `force` pointer; both are passed explicitly
/// here. `force_idx` indexes `rebx_extras::allocated_forces`, and the
/// particles are `sim.particles` (the C passes the same array, except
/// inside a REBOUNDx integrator where it passes a scratch copy — those
/// call sites pass their buffer through `rebx_extras` instead).
pub type rebx_force_fn =
    fn(sim: &mut reb_simulation, rebx: &mut rebx_extras, force_idx: usize, N: usize);

/// C: `void (*step_function)(struct reb_simulation* sim,
/// struct rebx_operator* operator, const double dt)`.
pub type rebx_operator_fn =
    fn(sim: &mut reb_simulation, rebx: &mut rebx_extras, operator_idx: usize, dt: f64);

/// reboundx.h `struct rebx_force`.
#[derive(Clone, Debug)]
pub struct rebx_force {
    /// For searching lists and informative errors.
    pub name: String,
    /// Additional parameters list (C: `struct rebx_node* ap`).
    pub ap: Vec<rebx_param>,
    /// Force type for internal logic.
    pub force_type: rebx_force_type,
    /// Adds this effect's accelerations (C: function pointer).
    pub update_accelerations: Option<rebx_force_fn>,
}

/// reboundx.h `struct rebx_operator`.
#[derive(Clone, Debug)]
pub struct rebx_operator {
    pub name: String,
    pub ap: Vec<rebx_param>,
    pub operator_type: rebx_operator_type,
    pub step_function: Option<rebx_operator_fn>,
}

/// reboundx.h `struct rebx_step`.
#[derive(Clone, Copy, Debug)]
pub struct rebx_step {
    /// Index into `rebx_extras::allocated_operators`
    /// (C: `struct rebx_operator* operator`).
    pub operator_: usize,
    /// Fraction of sim.dt to use each time it is called.
    pub dt_fraction: f64,
}

/// reboundx.h `struct rebx_interpolator`.
#[derive(Clone, Debug)]
pub struct rebx_interpolator {
    pub interpolation: rebx_interpolation_type,
    pub times: Vec<f64>,
    pub values: Vec<f64>,
    pub Nvalues: i32,
    pub y2: Vec<f64>,
    pub klo: i32,
}

/// reboundx.h `struct rebx_extras` — the REBOUNDx state, stored in
/// `reb_simulation::extras` (C: `sim->extras`).
///
/// Every C member that was a `struct rebx_node*` list is a `Vec` whose
/// index 0 is the list head, so iteration order matches the C exactly
/// (see the module docs: the C prepends).
#[derive(Clone, Debug, Default)]
pub struct rebx_extras {
    /// Forces evaluated each timestep: indices into `allocated_forces`.
    pub additional_forces: Vec<usize>,
    /// Steps applied before each timestep.
    pub pre_timestep_modifications: Vec<rebx_step>,
    /// Steps applied after each timestep.
    pub post_timestep_modifications: Vec<rebx_step>,
    /// All parameter names registered with their type (for type safety).
    pub registered_params: Vec<rebx_param>,
    /// Owns every force ever created (C: for memory management).
    pub allocated_forces: Vec<rebx_force>,
    /// Owns every operator ever created.
    pub allocated_operators: Vec<rebx_operator>,
    /// Per-particle parameter lists; `particle_params[i]` is C's
    /// `sim->particles[i].ap`. Grown on demand by the setters.
    pub particle_params: Vec<Vec<rebx_param>>,
    /// Diagnostics emitted by `rebx_error` (the C prints to stderr and,
    /// where a simulation is attached, calls `reb_simulation_error`).
    pub messages: Vec<String>,
}

impl rebx_extras {
    /// Borrow the parameter list named by `sel` (C: dereferencing the
    /// `struct rebx_node** apptr` argument).
    pub fn ap(&self, sel: rebx_ap) -> &[rebx_param] {
        match sel {
            rebx_ap::particle(i) => match self.particle_params.get(i) {
                Some(v) => v,
                None => &[],
            },
            rebx_ap::force(i) => match self.allocated_forces.get(i) {
                Some(f) => &f.ap,
                None => &[],
            },
            rebx_ap::operator_(i) => match self.allocated_operators.get(i) {
                Some(o) => &o.ap,
                None => &[],
            },
        }
    }

    /// Mutably borrow the parameter list named by `sel`, growing the
    /// per-particle table if needed (the C allocates the node lazily).
    pub fn ap_mut(&mut self, sel: rebx_ap) -> Option<&mut Vec<rebx_param>> {
        match sel {
            rebx_ap::particle(i) => {
                if self.particle_params.len() <= i {
                    self.particle_params.resize(i + 1, Vec::new());
                }
                self.particle_params.get_mut(i)
            }
            rebx_ap::force(i) => self.allocated_forces.get_mut(i).map(|f| &mut f.ap),
            rebx_ap::operator_(i) => self.allocated_operators.get_mut(i).map(|o| &mut o.ap),
        }
    }
}
