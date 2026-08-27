//! integrator_saba.rs — the SABA integrator family (from
//! integrator_saba.c/h; Laskar & Robutel 2001, Blanes et al. 2013,
//! Rein, Tamayo & Brown 2019). SABA1-4, corrected SABACM/SABACL 1-4,
//! and the generalized-order methods SABA(10,4), (8,6,4), (10,6,4),
//! SABAH(8,4,4), (8,6,4), (10,6,4). Built on the WHFast machinery.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein and contributors. See crate root.

use crate::integrator_whfast::*;
use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::reb_simulation_error;
use crate::transformations::*;
use crate::types::*;

pub const REB_INTEGRATOR_SABA_TYPE_1: i32 = 0x0;
pub const REB_INTEGRATOR_SABA_TYPE_2: i32 = 0x1;
pub const REB_INTEGRATOR_SABA_TYPE_3: i32 = 0x2;
pub const REB_INTEGRATOR_SABA_TYPE_4: i32 = 0x3;
pub const REB_INTEGRATOR_SABA_TYPE_CM_1: i32 = 0x100;
pub const REB_INTEGRATOR_SABA_TYPE_CM_2: i32 = 0x101;
pub const REB_INTEGRATOR_SABA_TYPE_CM_3: i32 = 0x102;
pub const REB_INTEGRATOR_SABA_TYPE_CM_4: i32 = 0x103;
pub const REB_INTEGRATOR_SABA_TYPE_CL_1: i32 = 0x200;
pub const REB_INTEGRATOR_SABA_TYPE_CL_2: i32 = 0x201;
pub const REB_INTEGRATOR_SABA_TYPE_CL_3: i32 = 0x202;
pub const REB_INTEGRATOR_SABA_TYPE_CL_4: i32 = 0x203;
pub const REB_INTEGRATOR_SABA_TYPE_10_4: i32 = 0x4;
pub const REB_INTEGRATOR_SABA_TYPE_8_6_4: i32 = 0x5;
pub const REB_INTEGRATOR_SABA_TYPE_10_6_4: i32 = 0x6;
pub const REB_INTEGRATOR_SABA_TYPE_H_8_4_4: i32 = 0x7;
pub const REB_INTEGRATOR_SABA_TYPE_H_8_6_4: i32 = 0x8;
pub const REB_INTEGRATOR_SABA_TYPE_H_10_6_4: i32 = 0x9;

/// integrator_saba.h `struct reb_integrator_saba_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_saba_state {
    pub type_: i32,
    pub safe_mode: u32,
    pub keep_unsynchronized: u32,
    // Internal use
    pub p_jh: Vec<reb_particle>,
    pub p_temp: Vec<reb_particle>,
}

impl Default for reb_integrator_saba_state {
    /// integrator_saba.c `reb_integrator_saba_create`.
    fn default() -> Self {
        reb_integrator_saba_state {
            type_: REB_INTEGRATOR_SABA_TYPE_10_6_4,
            safe_mode: 1,
            keep_unsynchronized: 0,
            p_jh: Vec::new(),
            p_temp: Vec::new(),
        }
    }
}

/// integrator_saba.c `reb_saba_stages`.
fn reb_saba_stages(type_: i32) -> i32 {
    match type_ {
        REB_INTEGRATOR_SABA_TYPE_1 | REB_INTEGRATOR_SABA_TYPE_CM_1 | REB_INTEGRATOR_SABA_TYPE_CL_1 => 1,
        REB_INTEGRATOR_SABA_TYPE_2 | REB_INTEGRATOR_SABA_TYPE_CM_2 | REB_INTEGRATOR_SABA_TYPE_CL_2 => 2,
        REB_INTEGRATOR_SABA_TYPE_3 | REB_INTEGRATOR_SABA_TYPE_CM_3 | REB_INTEGRATOR_SABA_TYPE_CL_3 => 3,
        REB_INTEGRATOR_SABA_TYPE_4 | REB_INTEGRATOR_SABA_TYPE_CM_4 | REB_INTEGRATOR_SABA_TYPE_CL_4 => 4,
        REB_INTEGRATOR_SABA_TYPE_H_8_4_4 => 6,
        REB_INTEGRATOR_SABA_TYPE_10_4 | REB_INTEGRATOR_SABA_TYPE_8_6_4 => 7,
        REB_INTEGRATOR_SABA_TYPE_10_6_4 | REB_INTEGRATOR_SABA_TYPE_H_8_6_4 => 8,
        REB_INTEGRATOR_SABA_TYPE_H_10_6_4 => 9,
        _ => 0,
    }
}

// Some coefficients appear multiple times to simplify the loop structures.
const reb_saba_c: [[f64; 5]; 10] = [
    [0.5, 0., 0., 0., 0.], // SABA1
    [0.2113248654051871177454256097490212721762, 0.5773502691896257645091487805019574556476, 0., 0., 0.], // SABA2
    [0.1127016653792583114820734600217600389167, 0.3872983346207416885179265399782399610833, 0., 0., 0.], // SABA3
    [0.06943184420297371238802675555359524745214, 0.2605776340045981552106403648947824089476,
        0.3399810435848562648026657591032446872006, 0., 0.], // SABA4
    [0.04706710064597250612947887637243678556564, 0.1847569354170881069247376193702560968574,
        0.2827060056798362053243616565541452479160, -0.01453004174289681837857815229683813033908, 0.], // ABA(10,4)
    [0.0711334264982231177779387300061549964174, 0.241153427956640098736487795326289649618,
        0.521411761772814789212136078067994229991, -0.333698616227678005726562603400438876027, 0.], // ABA(8,6,4)
    [0.03809449742241219545697532230863756534060, 0.1452987161169137492940200726606637497442,
        0.2076276957255412507162056113249882065158, 0.4359097036515261592231548624010651844006,
        -0.6538612258327867093807117373907094120024], // ABA(10,6,4)
    [0.2741402689434018761640565440378637101205, -0.1075684384401642306251105297063236526845,
        -0.0480185025906016926911954171508475065370, 0.7628933441747280943044988056386148982021, 0.], // ABAH(8,4,4)
    [0.06810235651658372084723976682061164571212, 0.2511360387221033233072829580455350680082,
        -0.07507264957216562516006821767601620052338, -0.009544719701745007811488218957217113269121,
        0.5307579480704471776340674235341732001443], // ABAH(8,6,4)
    [0.04731908697653382270404371796320813250988, 0.2651105235748785159539480036185693201078,
        -0.009976522883811240843267468164812380613143, -0.05992919973494155126395247987729676004016,
        0.2574761120673404534492282264603316880356], // ABAH(10,6,4)
];
const reb_saba_d: [[f64; 5]; 10] = [
    [1., 0., 0., 0., 0.],
    [0.5, 0., 0., 0., 0.],
    [0.2777777777777777777777777777777777777778, 0.4444444444444444444444444444444444444444, 0., 0., 0.],
    [0.1739274225687269286865319746109997036177, 0.3260725774312730713134680253890002963823, 0., 0., 0.],
    [0.1188819173681970199453503950853885936957, 0.2410504605515015657441667865901651105675,
        -0.2732866667053238060543113981664559460630, 0.8267085775712504407295884329818044835997, 0.], // ABA(10,4)
    [0.183083687472197221961703757166430291072, 0.310782859898574869507522291054262796375,
        -0.0265646185119588006972121379164987592663, 0.0653961422823734184559721793911134363710, 0.], // ABA(8,6,4)
    [0.09585888083707521061077150377145884776921, 0.2044461531429987806805077839164344779763,
        0.2170703479789911017143385924306336714532, -0.01737538195906509300561788011852699719871, 0.], // ABA(10,6,4)
    [0.6408857951625127177322491164716010349386, -0.8585754489567828565881283246356000103664,
        0.7176896537942701388558792081639989754277, 0., 0.], // ABAH(8,4,4)
    [0.1684432593618954534310382697756917558148, 0.4243177173742677224300351657407231801453,
        -0.5858109694681756812309015355404036521923, 0.4930499927320125053698281000239887162321, 0.], // ABAH(8,6,4)
    [0.1196884624585322035312864297489892143852, 0.3752955855379374250420128537687503199451,
        -0.4684593418325993783650820409805381740605, 0.3351397342755897010393098942949569049275,
        0.2766711191210800975049457263356834696055], // ABAH(10,6,4)
];
const reb_saba_cc: [f64; 4] = [
    0.08333333333333333333333333333333333333333,  // SABAC1
    0.01116454968463011276968973577058865137738,  // SABAC2
    0.005634593363122809402267823769797538671562, // SABAC3
    0.003396775048208601331532157783492144,       // SABAC4
];

/// integrator_saba.c `reb_saba_corrector_step`.
fn reb_saba_corrector_step(
    r: &mut reb_simulation,
    saba: &mut reb_integrator_saba_state,
    cc: f64,
) {
    let N = r.N;
    let mut empty_var: Vec<reb_particle> = Vec::new();
    match saba.type_ / 0x100 {
        1 => {
            // modified kick: calculate normal kick
            {
                let masses = r.particles.clone();
                reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &saba.p_jh, &masses, N, N);
            }
            reb_simulation_update_acceleration(r);
            // Calculate jerk. p_jh used as temporary buffer
            reb_integrator_whfast_calculate_jerk(r, &mut saba.p_jh);

            for i in 0..N {
                let prefact = r.dt * r.dt;
                r.particles[i].ax = prefact * saba.p_jh[i].ax;
                r.particles[i].ay = prefact * saba.p_jh[i].ay;
                r.particles[i].az = prefact * saba.p_jh[i].az;
            }
            let ccdt = cc * r.dt;
            reb_integrator_whfast_interaction_step(
                r,
                &mut saba.p_jh,
                &mut empty_var,
                REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI,
                ccdt,
            );
        }
        2 => {
            // lazy corrector
            if saba.p_temp.len() != N {
                saba.p_temp.resize(N, reb_particle::default());
            }

            // Calculate normal kick
            {
                let masses = r.particles.clone();
                reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &saba.p_jh, &masses, N, N);
            }
            reb_simulation_update_acceleration(r);
            reb_transformations_inertial_to_jacobi_acc(&r.particles.clone(), &mut saba.p_jh, &r.particles, N, N);

            // make copy of original positions and accelerations
            saba.p_temp.copy_from_slice(&saba.p_jh);

            // WHT96 Eq 10.6
            let prefac1 = r.dt * r.dt / 12.;
            for i in 1..N {
                saba.p_jh[i].x += prefac1 * saba.p_temp[i].ax;
                saba.p_jh[i].y += prefac1 * saba.p_temp[i].ay;
                saba.p_jh[i].z += prefac1 * saba.p_temp[i].az;
            }

            // recalculate kick
            {
                let masses = r.particles.clone();
                reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &saba.p_jh, &masses, N, N);
            }
            reb_simulation_update_acceleration(r);
            reb_transformations_inertial_to_jacobi_acc(&r.particles.clone(), &mut saba.p_jh, &r.particles, N, N);

            let prefact = cc * r.dt * 12.;
            for i in 1..N {
                // Lazy implementer's commutator
                saba.p_jh[i].vx += prefact * (saba.p_jh[i].ax - saba.p_temp[i].ax);
                saba.p_jh[i].vy += prefact * (saba.p_jh[i].ay - saba.p_temp[i].ay);
                saba.p_jh[i].vz += prefact * (saba.p_jh[i].az - saba.p_temp[i].az);
                // reset positions
                saba.p_jh[i].x = saba.p_temp[i].x;
                saba.p_jh[i].y = saba.p_temp[i].y;
                saba.p_jh[i].z = saba.p_temp[i].z;
            }
        }
        _ => {}
    }
}

/// integrator_saba.c `reb_integrator_saba_step` (state-explicit).
pub fn reb_integrator_saba_step_state(
    r: &mut reb_simulation,
    saba: &mut reb_integrator_saba_state,
) {
    let type_ = saba.type_;
    let stages = reb_saba_stages(type_);
    let N = r.N;
    let mut empty_var: Vec<reb_particle> = Vec::new();
    if !r.var_config.is_empty() {
        reb_simulation_error(r, "Variational particles are not supported in the SABA integrator.");
        return;
    }
    if saba.keep_unsynchronized == 1 && saba.safe_mode == 1 {
        reb_simulation_error(r, "saba->keep_unsynchronized == 1 is not compatible with safe_mode. Must set saba->safe_mode = 0.");
    }
    let valid = matches!(
        type_,
        0x0 | 0x1 | 0x2 | 0x3 | 0x100 | 0x101 | 0x102 | 0x103 | 0x200 | 0x201 | 0x202 | 0x203
            | 0x4 | 0x5 | 0x6 | 0x7 | 0x8 | 0x9
    );
    if !valid {
        reb_simulation_error(r, "Invalid SABA integrator type used.");
        return;
    }
    if type_ >= 0x100 {
        // Force Jacobi terms in update_acceleration when corrector is used
        r.gravity = REB_GRAVITY::JACOBI;
    } else {
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1;
    }
    if saba.p_jh.len() != N {
        saba.p_jh.resize(N, reb_particle::default());
        r.did_modify_particles = 1;
    }

    // Only recalculate Jacobi coordinates if needed
    if saba.safe_mode != 0 || r.did_modify_particles != 0 {
        reb_integrator_whfast_from_inertial(
            r,
            &mut saba.p_jh,
            &mut empty_var,
            REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI,
        );
    }
    let row = (type_ % 0x100) as usize;
    if type_ >= 0x100 {
        // Correctors on
        if r.is_synchronized != 0 {
            reb_saba_corrector_step(r, saba, reb_saba_cc[row]);
        } else {
            reb_saba_corrector_step(r, saba, 2. * reb_saba_cc[row]);
        }
        // First half DRIFT step
        let d0 = reb_saba_c[row][0] * r.dt;
        reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);
        reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, d0);
    } else {
        // Correctors off
        if r.is_synchronized != 0 {
            // First half DRIFT step
            let d0 = reb_saba_c[row][0] * r.dt;
            reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);
            reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, d0);
        } else {
            // Combined DRIFT step
            let d0 = 2. * reb_saba_c[row][0] * r.dt;
            reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);
            reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, d0);
        }
    }

    reb_integrator_whfast_to_inertial(r, &saba.p_jh, &empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI);

    reb_simulation_update_acceleration(r);

    let d0 = reb_saba_d[row][0] * r.dt;
    reb_integrator_whfast_interaction_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);

    for j in 1..stages {
        {
            let mut i = j;
            if j > stages / 2 {
                i = stages - j;
            }
            let dc = reb_saba_c[row][i as usize] * r.dt;
            reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, dc);
            reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, dc);
        }
        {
            let mut i = j;
            if j > (stages - 1) / 2 {
                i = stages - j - 1;
            }
            {
                let masses = r.particles.clone();
                reb_transformations_jacobi_to_inertial_pos(&mut r.particles, &saba.p_jh, &masses, N, N);
            }
            reb_simulation_update_acceleration(r);
            let dd = reb_saba_d[row][i as usize] * r.dt;
            reb_integrator_whfast_interaction_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, dd);
        }
    }

    if saba.type_ >= 0x100 {
        // correctors on: always need to do drift step
        let d0 = reb_saba_c[row][0] * r.dt;
        reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);
        reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, d0);
    }

    r.is_synchronized = 0;
    if saba.safe_mode != 0 {
        reb_integrator_saba_synchronize_state(r, saba);
    }

    r.t += r.dt;
    r.dt_last_done = r.dt;
}

/// integrator_saba.c `reb_integrator_saba_synchronize` (state-explicit).
pub fn reb_integrator_saba_synchronize_state(
    r: &mut reb_simulation,
    saba: &mut reb_integrator_saba_state,
) {
    let type_ = saba.type_;
    let mut sync_pj: Option<Vec<reb_particle>> = None;
    if saba.keep_unsynchronized != 0 {
        sync_pj = Some(saba.p_jh.clone());
    }
    if r.is_synchronized == 0 {
        let N = r.N;
        let mut empty_var: Vec<reb_particle> = Vec::new();
        let row = (type_ % 0x100) as usize;
        if type_ >= 0x100 {
            // correctors on: drift already done, just need corrector
            reb_saba_corrector_step(r, saba, reb_saba_cc[row]);
        } else {
            let d0 = reb_saba_c[row][0] * r.dt;
            reb_integrator_whfast_kepler_step(r, &mut saba.p_jh, &mut empty_var, REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI, d0);
            reb_integrator_whfast_com_step(r, &mut saba.p_jh, &mut empty_var, d0);
        }
        {
            let masses = r.particles.clone();
            reb_transformations_jacobi_to_inertial_posvel(&mut r.particles, &saba.p_jh, &masses, N, N);
        }
        if let Some(saved) = sync_pj {
            saba.p_jh = saved;
        } else {
            r.is_synchronized = 1;
        }
    }
}

/// Step entry point for the dispatcher.
pub fn reb_integrator_saba_step(r: &mut reb_simulation) {
    let mut saba = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::saba(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_saba_step_state(r, &mut saba);
    r.integrator = reb_integrator_state::saba(saba);
}

/// Synchronize entry point for the dispatcher.
pub fn reb_integrator_saba_synchronize(r: &mut reb_simulation) {
    let mut saba = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
        reb_integrator_state::saba(s) => s,
        other => {
            r.integrator = other;
            return;
        }
    };
    reb_integrator_saba_synchronize_state(r, &mut saba);
    r.integrator = reb_integrator_state::saba(saba);
}
