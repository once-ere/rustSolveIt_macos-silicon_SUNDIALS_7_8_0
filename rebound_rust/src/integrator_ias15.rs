//! integrator_ias15.rs — the IAS15 integrator (from integrator_ias15.c
//! and integrator_ias15.h; Rein & Spiegel 2015, Everhart 1985, timestep
//! criterion Pham, Rein & Spiegel 2024).
//!
//! The C stores seven-segment coefficient arrays (`reb_dp7`) as one
//! malloc'd block of 7*N3 doubles addressed through seven pointers;
//! here each is a `Vec<f64>` of length 7*N3 with segment s at
//! `[s*N3 .. (s+1)*N3]` — the identical memory layout.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Dave Spiegel, Pasquale Tricarico and contributors.
//! See crate root.

use crate::simulation::reb_simulation_update_acceleration;
use crate::tools::{reb_simulation_warning, reb_tools_megno_deltad_delta, reb_tools_megno_update};
use crate::types::*;

pub const REB_IAS15_ADAPTIVEMODE_INDIVIDUAL: u32 = 0;
pub const REB_IAS15_ADAPTIVEMODE_GLOBAL: u32 = 1;
pub const REB_IAS15_ADAPTIVEMODE_PRS23: u32 = 2;
pub const REB_IAS15_ADAPTIVEMODE_AARSETH85: u32 = 3;

/// integrator_ias15.h `struct reb_integrator_ias15_state`.
#[derive(Clone, Debug)]
pub struct reb_integrator_ias15_state {
    pub epsilon: f64,
    pub min_dt: f64,
    pub adaptive_mode: u32,
    pub iterations_max_exceeded: u64,
    pub N_allocated: usize,
    pub at: Vec<f64>,
    pub x0: Vec<f64>,
    pub v0: Vec<f64>,
    pub a0: Vec<f64>,
    pub csx: Vec<f64>,
    pub csv: Vec<f64>,
    pub csa0: Vec<f64>,
    pub g: Vec<f64>,
    pub b: Vec<f64>,
    pub csb: Vec<f64>,
    pub e: Vec<f64>,
    pub br: Vec<f64>,
    pub er: Vec<f64>,
}

impl Default for reb_integrator_ias15_state {
    /// integrator_ias15.c `reb_integrator_ias15_create`.
    fn default() -> Self {
        reb_integrator_ias15_state {
            epsilon: 1e-9,
            min_dt: 0.0,
            adaptive_mode: REB_IAS15_ADAPTIVEMODE_PRS23,
            iterations_max_exceeded: 0,
            N_allocated: 0,
            at: Vec::new(),
            x0: Vec::new(),
            v0: Vec::new(),
            a0: Vec::new(),
            csx: Vec::new(),
            csv: Vec::new(),
            csa0: Vec::new(),
            g: Vec::new(),
            b: Vec::new(),
            csb: Vec::new(),
            e: Vec::new(),
            br: Vec::new(),
            er: Vec::new(),
        }
    }
}

/// Maximum increase/decrease of consecutive timesteps.
const safety_factor: f64 = 0.25;

/// Gauss-Radau spacings.
const h: [f64; 8] = [
    0.0,
    0.0562625605369221464656521910318,
    0.180240691736892364987579942780,
    0.352624717113169637373907769648,
    0.547153626330555383001448554766,
    0.734210177215410531523210605558,
    0.885320946839095768090359771030,
    0.977520613561287501891174488626,
];
const rr: [f64; 28] = [
    0.0562625605369221464656522,
    0.1802406917368923649875799,
    0.1239781311999702185219278,
    0.3526247171131696373739078,
    0.2963621565762474909082556,
    0.1723840253762772723863278,
    0.5471536263305553830014486,
    0.4908910657936332365357964,
    0.3669129345936630180138686,
    0.1945289092173857456275408,
    0.7342101772154105315232106,
    0.6779476166784883850575584,
    0.5539694854785181665356307,
    0.3815854601022408941493028,
    0.1870565508848551485217621,
    0.8853209468390957680903598,
    0.8290583863021736216247076,
    0.7050802551022034031027798,
    0.5326962297259261307164520,
    0.3381673205085403850889112,
    0.1511107696236852365671492,
    0.9775206135612875018911745,
    0.9212580530243653554255223,
    0.7972799218243951369035945,
    0.6248958964481178645172667,
    0.4303669872307321188897259,
    0.2433104363458769703679639,
    0.0921996667221917338008147,
];
const c: [f64; 21] = [
    -0.0562625605369221464656522,
    0.0101408028300636299864818,
    -0.2365032522738145114532321,
    -0.0035758977292516175949345,
    0.0935376952594620658957485,
    -0.5891279693869841488271399,
    0.0019565654099472210769006,
    -0.0547553868890686864408084,
    0.4158812000823068616886219,
    -1.1362815957175395318285885,
    -0.0014365302363708915424460,
    0.0421585277212687077072973,
    -0.3600995965020568122897665,
    1.2501507118406910258505441,
    -1.8704917729329500633517991,
    0.0012717903090268677492943,
    -0.0387603579159067703699046,
    0.3609622434528459832253398,
    -1.4668842084004269643701553,
    2.9061362593084293014237913,
    -2.7558127197720458314421588,
];
const d: [f64; 21] = [
    0.0562625605369221464656522,
    0.0031654757181708292499905,
    0.2365032522738145114532321,
    0.0001780977692217433881125,
    0.0457929855060279188954539,
    0.5891279693869841488271399,
    0.0000100202365223291272096,
    0.0084318571535257015445000,
    0.2535340690545692665214616,
    1.1362815957175395318285885,
    0.0000005637641639318207610,
    0.0015297840025004658189490,
    0.0978342365324440053653648,
    0.8752546646840910912297246,
    1.8704917729329500633517991,
    0.0000000317188154017613665,
    0.0002762930909826476593130,
    0.0360285539837364596003871,
    0.5767330002770787313544596,
    2.2485887607691597933926895,
    2.7558127197720458314421588,
];

/// Weights for integration of a first order differential equation.
const w: [f64; 8] = [
    0.03125,
    0.185358154802979278540728972807180754479812609,
    0.304130620646785128975743291458180383736715043,
    0.376517545389118556572129261157225608762708603,
    0.391572167452493593082499533303669362149363727,
    0.347014795634501068709955597003528601733139176,
    0.249647901329864963257869294715235590174262844,
    0.114508814744257199342353731044292225247093225,
];

/// integrator_ias15.c `sqrt7` — machine independent pow(a, 1/7).
fn sqrt7(mut a: f64) -> f64 {
    // Without scaling accurate for [1e-7, 1e2]; with scaling [1e-14, 1e8].
    let mut scale = 1.0;
    while a < 1e-7 && a.is_normal() {
        scale *= 0.1;
        a *= 1e7;
    }
    while a > 1e2 && a.is_normal() {
        scale *= 10.;
        a *= 1e-7;
    }
    let mut x = 1.;
    for _k in 0..20 {
        let x6 = x * x * x * x * x * x;
        x += (a / x6 - x) / 7.;
    }
    x * scale
}

/// integrator_ias15.c `add_cs` — one compensated-summation add.
#[inline]
fn add_cs(p: &mut f64, csp: &mut f64, inp: f64) {
    let y = inp - *csp;
    let t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

fn realloc_dp7(p: &mut Vec<f64>, N3: usize) {
    p.clear();
    p.resize(N3 * 7, 0.0);
}

/// integrator_ias15.c `reb_integrator_ias15_alloc`.
fn reb_integrator_ias15_alloc(r: &reb_simulation, ias15: &mut reb_integrator_ias15_state) {
    let N = if r.map.is_some() { r.N_map } else { r.N + r.N_var };
    let N3 = 3 * N;
    if N3 > ias15.N_allocated {
        realloc_dp7(&mut ias15.g, N3);
        realloc_dp7(&mut ias15.e, N3);
        realloc_dp7(&mut ias15.b, N3);
        realloc_dp7(&mut ias15.csb, N3);
        realloc_dp7(&mut ias15.er, N3);
        realloc_dp7(&mut ias15.br, N3);
        ias15.at.resize(N3, 0.0);
        ias15.x0.resize(N3, 0.0);
        ias15.v0.resize(N3, 0.0);
        ias15.a0.resize(N3, 0.0);
        ias15.csx.resize(N3, 0.0);
        ias15.csv.resize(N3, 0.0);
        ias15.csa0.resize(N3, 0.0);
        for i in 0..N3 {
            // Kill compensated summation coefficients
            ias15.csx[i] = 0.;
            ias15.csv[i] = 0.;
        }
        ias15.N_allocated = N3;
    }
}

/// The flattened `((double*)gravity_cs)[k]` read: element k of the
/// vec3d array viewed as doubles. `use_sim_cs` is the
/// REB_GRAVITY_COMPENSATED case; otherwise the C aliases `csa0`
/// (all zeros at that point).
#[inline]
fn gravity_cs_flat(r: &reb_simulation, csa0: &[f64], use_sim_cs: bool, k: usize) -> f64 {
    if use_sim_cs {
        let v = r.gravity_cs[k / 3];
        match k % 3 {
            0 => v.x,
            1 => v.y,
            _ => v.z,
        }
    } else {
        csa0[k]
    }
}

/// integrator_ias15.c `predict_next_step`. `src=true` predicts from the
/// retained er/br buffers, `src=false` from e/b in place (the C passes
/// the source and destination dp7s; in Rust the two cases are made
/// explicit to satisfy aliasing).
fn predict_next_step(
    ratio: f64,
    N3: usize,
    _e: &[f64],
    _b: &[f64],
    e: &mut [f64],
    b: &mut [f64],
) {
    if ratio > 20. {
        // Do not predict if stepsize increase is very large.
        for k in 0..N3 {
            for s in 0..7 {
                e[s * N3 + k] = 0.;
                b[s * N3 + k] = 0.;
            }
        }
    } else {
        let q1 = ratio;
        let q2 = q1 * q1;
        let q3 = q1 * q2;
        let q4 = q2 * q2;
        let q5 = q2 * q3;
        let q6 = q3 * q3;
        let q7 = q3 * q4;

        for k in 0..N3 {
            let _b0 = _b[0 * N3 + k];
            let _b1 = _b[1 * N3 + k];
            let _b2 = _b[2 * N3 + k];
            let _b3 = _b[3 * N3 + k];
            let _b4 = _b[4 * N3 + k];
            let _b5 = _b[5 * N3 + k];
            let _b6 = _b[6 * N3 + k];
            let be0 = _b0 - _e[0 * N3 + k];
            let be1 = _b1 - _e[1 * N3 + k];
            let be2 = _b2 - _e[2 * N3 + k];
            let be3 = _b3 - _e[3 * N3 + k];
            let be4 = _b4 - _e[4 * N3 + k];
            let be5 = _b5 - _e[5 * N3 + k];
            let be6 = _b6 - _e[6 * N3 + k];

            e[0 * N3 + k] = q1 * (_b6 * 7.0 + _b5 * 6.0 + _b4 * 5.0 + _b3 * 4.0 + _b2 * 3.0 + _b1 * 2.0 + _b0);
            e[1 * N3 + k] = q2 * (_b6 * 21.0 + _b5 * 15.0 + _b4 * 10.0 + _b3 * 6.0 + _b2 * 3.0 + _b1);
            e[2 * N3 + k] = q3 * (_b6 * 35.0 + _b5 * 20.0 + _b4 * 10.0 + _b3 * 4.0 + _b2);
            e[3 * N3 + k] = q4 * (_b6 * 35.0 + _b5 * 15.0 + _b4 * 5.0 + _b3);
            e[4 * N3 + k] = q5 * (_b6 * 21.0 + _b5 * 6.0 + _b4);
            e[5 * N3 + k] = q6 * (_b6 * 7.0 + _b5);
            e[6 * N3 + k] = q7 * _b6;

            b[0 * N3 + k] = e[0 * N3 + k] + be0;
            b[1 * N3 + k] = e[1 * N3 + k] + be1;
            b[2 * N3 + k] = e[2 * N3 + k] + be2;
            b[3 * N3 + k] = e[3 * N3 + k] + be3;
            b[4 * N3 + k] = e[4 * N3 + k] + be4;
            b[5 * N3 + k] = e[5 * N3 + k] + be5;
            b[6 * N3 + k] = e[6 * N3 + k] + be6;
        }
    }
}

/// integrator_ias15.c `predict_next_step` for the in-place case
/// (source == destination, C's final `predict_next_step(ratio, N3, e,
/// b, e, b)`): the C reads `_b`/`_e` element k before overwriting it,
/// so a per-k snapshot reproduces it exactly.
fn predict_next_step_inplace(ratio: f64, N3: usize, e: &mut [f64], b: &mut [f64]) {
    if ratio > 20. {
        for k in 0..N3 {
            for s in 0..7 {
                e[s * N3 + k] = 0.;
                b[s * N3 + k] = 0.;
            }
        }
    } else {
        let q1 = ratio;
        let q2 = q1 * q1;
        let q3 = q1 * q2;
        let q4 = q2 * q2;
        let q5 = q2 * q3;
        let q6 = q3 * q3;
        let q7 = q3 * q4;

        for k in 0..N3 {
            let _b0 = b[0 * N3 + k];
            let _b1 = b[1 * N3 + k];
            let _b2 = b[2 * N3 + k];
            let _b3 = b[3 * N3 + k];
            let _b4 = b[4 * N3 + k];
            let _b5 = b[5 * N3 + k];
            let _b6 = b[6 * N3 + k];
            let be0 = _b0 - e[0 * N3 + k];
            let be1 = _b1 - e[1 * N3 + k];
            let be2 = _b2 - e[2 * N3 + k];
            let be3 = _b3 - e[3 * N3 + k];
            let be4 = _b4 - e[4 * N3 + k];
            let be5 = _b5 - e[5 * N3 + k];
            let be6 = _b6 - e[6 * N3 + k];

            e[0 * N3 + k] = q1 * (_b6 * 7.0 + _b5 * 6.0 + _b4 * 5.0 + _b3 * 4.0 + _b2 * 3.0 + _b1 * 2.0 + _b0);
            e[1 * N3 + k] = q2 * (_b6 * 21.0 + _b5 * 15.0 + _b4 * 10.0 + _b3 * 6.0 + _b2 * 3.0 + _b1);
            e[2 * N3 + k] = q3 * (_b6 * 35.0 + _b5 * 20.0 + _b4 * 10.0 + _b3 * 4.0 + _b2);
            e[3 * N3 + k] = q4 * (_b6 * 35.0 + _b5 * 15.0 + _b4 * 5.0 + _b3);
            e[4 * N3 + k] = q5 * (_b6 * 21.0 + _b5 * 6.0 + _b4);
            e[5 * N3 + k] = q6 * (_b6 * 7.0 + _b5);
            e[6 * N3 + k] = q7 * _b6;

            b[0 * N3 + k] = e[0 * N3 + k] + be0;
            b[1 * N3 + k] = e[1 * N3 + k] + be1;
            b[2 * N3 + k] = e[2 * N3 + k] + be2;
            b[3 * N3 + k] = e[3 * N3 + k] + be3;
            b[4 * N3 + k] = e[4 * N3 + k] + be4;
            b[5 * N3 + k] = e[5 * N3 + k] + be5;
            b[6 * N3 + k] = e[6 * N3 + k] + be6;
        }
    }
}

/// integrator_ias15.c `copybuffers`.
fn copybuffers(a: &[f64], b: &mut [f64], N3: usize) {
    b[..N3 * 7].copy_from_slice(&a[..N3 * 7]);
}

#[inline]
fn map_index(r: &reb_simulation, k: usize) -> usize {
    match &r.map {
        Some(m) => m[k],
        None => k,
    }
}

/// integrator_ias15.c `reb_integrator_ias15_step_try`.
/// Returns true if the step was accepted.
fn reb_integrator_ias15_step_try(
    r: &mut reb_simulation,
    ias15: &mut reb_integrator_ias15_state,
) -> bool {
    reb_integrator_ias15_alloc(r, ias15);

    let (N, N_var) = if r.map.is_some() {
        (r.N_map, 0usize)
    } else {
        (r.N, r.N_var)
    };
    let N3 = 3 * (N + N_var);

    reb_simulation_update_acceleration(r);

    for k in 0..N {
        let mk = map_index(r, k);
        ias15.x0[3 * k] = r.particles[mk].x;
        ias15.x0[3 * k + 1] = r.particles[mk].y;
        ias15.x0[3 * k + 2] = r.particles[mk].z;
        ias15.v0[3 * k] = r.particles[mk].vx;
        ias15.v0[3 * k + 1] = r.particles[mk].vy;
        ias15.v0[3 * k + 2] = r.particles[mk].vz;
        ias15.a0[3 * k] = r.particles[mk].ax;
        ias15.a0[3 * k + 1] = r.particles[mk].ay;
        ias15.a0[3 * k + 2] = r.particles[mk].az;
    }
    for mk in 0..N_var {
        let k = mk + N;
        ias15.x0[3 * k] = r.particles_var[mk].x;
        ias15.x0[3 * k + 1] = r.particles_var[mk].y;
        ias15.x0[3 * k + 2] = r.particles_var[mk].z;
        ias15.v0[3 * k] = r.particles_var[mk].vx;
        ias15.v0[3 * k + 1] = r.particles_var[mk].vy;
        ias15.v0[3 * k + 2] = r.particles_var[mk].vz;
        ias15.a0[3 * k] = r.particles_var[mk].ax;
        ias15.a0[3 * k + 1] = r.particles_var[mk].ay;
        ias15.a0[3 * k + 2] = r.particles_var[mk].az;
    }
    let use_sim_cs = r.gravity == REB_GRAVITY::COMPENSATED;
    if use_sim_cs {
        for k in 0..(N + N_var) {
            let mk = map_index(r, k);
            ias15.csa0[3 * k] = r.gravity_cs[mk].x;
            ias15.csa0[3 * k + 1] = r.gravity_cs[mk].y;
            ias15.csa0[3 * k + 2] = r.gravity_cs[mk].z;
        }
    } else {
        for k in 0..N3 {
            ias15.csa0[k] = 0.;
        }
    }
    for k in 0..N3 {
        for s in 0..7 {
            ias15.csb[s * N3 + k] = 0.;
        }
    }

    {
        let g = &mut ias15.g;
        let b = &ias15.b;
        for k in 0..N3 {
            g[0 * N3 + k] = b[6 * N3 + k] * d[15] + b[5 * N3 + k] * d[10] + b[4 * N3 + k] * d[6] + b[3 * N3 + k] * d[3] + b[2 * N3 + k] * d[1] + b[1 * N3 + k] * d[0] + b[0 * N3 + k];
            g[1 * N3 + k] = b[6 * N3 + k] * d[16] + b[5 * N3 + k] * d[11] + b[4 * N3 + k] * d[7] + b[3 * N3 + k] * d[4] + b[2 * N3 + k] * d[2] + b[1 * N3 + k];
            g[2 * N3 + k] = b[6 * N3 + k] * d[17] + b[5 * N3 + k] * d[12] + b[4 * N3 + k] * d[8] + b[3 * N3 + k] * d[5] + b[2 * N3 + k];
            g[3 * N3 + k] = b[6 * N3 + k] * d[18] + b[5 * N3 + k] * d[13] + b[4 * N3 + k] * d[9] + b[3 * N3 + k];
            g[4 * N3 + k] = b[6 * N3 + k] * d[19] + b[5 * N3 + k] * d[14] + b[4 * N3 + k];
            g[5 * N3 + k] = b[6 * N3 + k] * d[20] + b[5 * N3 + k];
            g[6 * N3 + k] = b[6 * N3 + k];
        }
    }

    let mut integrator_megno_thisdt = 0.;
    let mut integrator_megno_thisdt_init = 0.;
    if r.calculate_megno != 0 {
        integrator_megno_thisdt_init =
            w[0] * (r.t - r.megno_initial_t) * reb_tools_megno_deltad_delta(r);
    }

    let t_beginning = r.t;
    let mut predictor_corrector_error: f64 = 1e300;
    let mut predictor_corrector_error_last: f64 = 2.;
    let mut iterations: usize = 0;
    // Predictor corrector loop
    loop {
        if predictor_corrector_error < 1e-16 {
            break;
        }
        if iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error {
            break;
        }
        if iterations >= 12 {
            ias15.iterations_max_exceeded += 1;
            let integrator_iterations_warning = 10;
            if ias15.iterations_max_exceeded == integrator_iterations_warning {
                reb_simulation_warning(r, "At least 10 predictor corrector loops in IAS15 did not converge. This is typically an indication of the timestep being too large.");
            }
            break; // Quit predictor corrector loop
        }
        predictor_corrector_error_last = predictor_corrector_error;
        predictor_corrector_error = 0.;
        iterations += 1;

        integrator_megno_thisdt = integrator_megno_thisdt_init;

        for n in 1..8usize {
            // Loop over interval using Gauss-Radau spacings
            r.t = t_beginning + r.dt * h[n];

            // Predict positions at interval n using b values
            for i in 0..(N + N_var) {
                let k0 = 3 * i;
                let k1 = 3 * i + 1;
                let k2 = 3 * i + 2;
                let b = &ias15.b;
                let xk0 = -ias15.csx[k0]
                    + ((((((((b[6 * N3 + k0] * 7. * h[n] / 9. + b[5 * N3 + k0]) * 3. * h[n] / 4.
                        + b[4 * N3 + k0]) * 5. * h[n] / 7.
                        + b[3 * N3 + k0]) * 2. * h[n] / 3.
                        + b[2 * N3 + k0]) * 3. * h[n] / 5.
                        + b[1 * N3 + k0]) * h[n] / 2.
                        + b[0 * N3 + k0]) * h[n] / 3.
                        + ias15.a0[k0]) * r.dt * h[n] / 2.
                        + ias15.v0[k0]) * r.dt * h[n];
                let xk1 = -ias15.csx[k1]
                    + ((((((((b[6 * N3 + k1] * 7. * h[n] / 9. + b[5 * N3 + k1]) * 3. * h[n] / 4.
                        + b[4 * N3 + k1]) * 5. * h[n] / 7.
                        + b[3 * N3 + k1]) * 2. * h[n] / 3.
                        + b[2 * N3 + k1]) * 3. * h[n] / 5.
                        + b[1 * N3 + k1]) * h[n] / 2.
                        + b[0 * N3 + k1]) * h[n] / 3.
                        + ias15.a0[k1]) * r.dt * h[n] / 2.
                        + ias15.v0[k1]) * r.dt * h[n];
                let xk2 = -ias15.csx[k2]
                    + ((((((((b[6 * N3 + k2] * 7. * h[n] / 9. + b[5 * N3 + k2]) * 3. * h[n] / 4.
                        + b[4 * N3 + k2]) * 5. * h[n] / 7.
                        + b[3 * N3 + k2]) * 2. * h[n] / 3.
                        + b[2 * N3 + k2]) * 3. * h[n] / 5.
                        + b[1 * N3 + k2]) * h[n] / 2.
                        + b[0 * N3 + k2]) * h[n] / 3.
                        + ias15.a0[k2]) * r.dt * h[n] / 2.
                        + ias15.v0[k2]) * r.dt * h[n];
                if i < N {
                    let mi = map_index(r, i);
                    r.particles[mi].x = xk0 + ias15.x0[k0];
                    r.particles[mi].y = xk1 + ias15.x0[k1];
                    r.particles[mi].z = xk2 + ias15.x0[k2];
                } else {
                    let mi = i - N;
                    r.particles_var[mi].x = xk0 + ias15.x0[k0];
                    r.particles_var[mi].y = xk1 + ias15.x0[k1];
                    r.particles_var[mi].z = xk2 + ias15.x0[k2];
                }
            }
            if r.calculate_megno != 0
                || (r.additional_forces.is_some() && r.force_is_velocity_dependent != 0)
            {
                // Predict velocities at interval n using b values
                for i in 0..(N + N_var) {
                    let k0 = 3 * i;
                    let k1 = 3 * i + 1;
                    let k2 = 3 * i + 2;
                    let b = &ias15.b;
                    let vk0 = -ias15.csv[k0]
                        + (((((((b[6 * N3 + k0] * 7. * h[n] / 8. + b[5 * N3 + k0]) * 6. * h[n] / 7.
                            + b[4 * N3 + k0]) * 5. * h[n] / 6.
                            + b[3 * N3 + k0]) * 4. * h[n] / 5.
                            + b[2 * N3 + k0]) * 3. * h[n] / 4.
                            + b[1 * N3 + k0]) * 2. * h[n] / 3.
                            + b[0 * N3 + k0]) * h[n] / 2.
                            + ias15.a0[k0]) * r.dt * h[n];
                    let vk1 = -ias15.csv[k1]
                        + (((((((b[6 * N3 + k1] * 7. * h[n] / 8. + b[5 * N3 + k1]) * 6. * h[n] / 7.
                            + b[4 * N3 + k1]) * 5. * h[n] / 6.
                            + b[3 * N3 + k1]) * 4. * h[n] / 5.
                            + b[2 * N3 + k1]) * 3. * h[n] / 4.
                            + b[1 * N3 + k1]) * 2. * h[n] / 3.
                            + b[0 * N3 + k1]) * h[n] / 2.
                            + ias15.a0[k1]) * r.dt * h[n];
                    let vk2 = -ias15.csv[k2]
                        + (((((((b[6 * N3 + k2] * 7. * h[n] / 8. + b[5 * N3 + k2]) * 6. * h[n] / 7.
                            + b[4 * N3 + k2]) * 5. * h[n] / 6.
                            + b[3 * N3 + k2]) * 4. * h[n] / 5.
                            + b[2 * N3 + k2]) * 3. * h[n] / 4.
                            + b[1 * N3 + k2]) * 2. * h[n] / 3.
                            + b[0 * N3 + k2]) * h[n] / 2.
                            + ias15.a0[k2]) * r.dt * h[n];
                    if i < N {
                        let mi = map_index(r, i);
                        r.particles[mi].vx = vk0 + ias15.v0[k0];
                        r.particles[mi].vy = vk1 + ias15.v0[k1];
                        r.particles[mi].vz = vk2 + ias15.v0[k2];
                    } else {
                        let mi = i - N;
                        r.particles_var[mi].vx = vk0 + ias15.v0[k0];
                        r.particles_var[mi].vy = vk1 + ias15.v0[k1];
                        r.particles_var[mi].vz = vk2 + ias15.v0[k2];
                    }
                }
            }

            reb_simulation_update_acceleration(r); // Calculate forces at interval n
            if r.calculate_megno != 0 {
                integrator_megno_thisdt +=
                    w[n] * (r.t - r.megno_initial_t) * reb_tools_megno_deltad_delta(r);
            }

            for k in 0..N {
                let mk = map_index(r, k);
                ias15.at[3 * k] = r.particles[mk].ax;
                ias15.at[3 * k + 1] = r.particles[mk].ay;
                ias15.at[3 * k + 2] = r.particles[mk].az;
            }
            for mk in 0..N_var {
                let k = mk + N;
                ias15.at[3 * k] = r.particles_var[mk].ax;
                ias15.at[3 * k + 1] = r.particles_var[mk].ay;
                ias15.at[3 * k + 2] = r.particles_var[mk].az;
            }
            match n {
                1 => {
                    for k in 0..N3 {
                        let tmp = ias15.g[0 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[0 * N3 + k] = gk / rr[0];
                        let delta = ias15.g[0 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], delta);
                    }
                }
                2 => {
                    for k in 0..N3 {
                        let mut tmp = ias15.g[1 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[1 * N3 + k] = (gk / rr[1] - ias15.g[0 * N3 + k]) / rr[2];
                        tmp = ias15.g[1 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[0]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp);
                    }
                }
                3 => {
                    for k in 0..N3 {
                        let mut tmp = ias15.g[2 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[2 * N3 + k] =
                            ((gk / rr[3] - ias15.g[0 * N3 + k]) / rr[4] - ias15.g[1 * N3 + k]) / rr[5];
                        tmp = ias15.g[2 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[1]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp * c[2]);
                        add_cs(&mut b[2 * N3 + k], &mut csb[2 * N3 + k], tmp);
                    }
                }
                4 => {
                    for k in 0..N3 {
                        let mut tmp = ias15.g[3 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[3 * N3 + k] = (((gk / rr[6] - ias15.g[0 * N3 + k]) / rr[7]
                            - ias15.g[1 * N3 + k]) / rr[8]
                            - ias15.g[2 * N3 + k]) / rr[9];
                        tmp = ias15.g[3 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[3]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp * c[4]);
                        add_cs(&mut b[2 * N3 + k], &mut csb[2 * N3 + k], tmp * c[5]);
                        add_cs(&mut b[3 * N3 + k], &mut csb[3 * N3 + k], tmp);
                    }
                }
                5 => {
                    for k in 0..N3 {
                        let mut tmp = ias15.g[4 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[4 * N3 + k] = ((((gk / rr[10] - ias15.g[0 * N3 + k]) / rr[11]
                            - ias15.g[1 * N3 + k]) / rr[12]
                            - ias15.g[2 * N3 + k]) / rr[13]
                            - ias15.g[3 * N3 + k]) / rr[14];
                        tmp = ias15.g[4 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[6]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp * c[7]);
                        add_cs(&mut b[2 * N3 + k], &mut csb[2 * N3 + k], tmp * c[8]);
                        add_cs(&mut b[3 * N3 + k], &mut csb[3 * N3 + k], tmp * c[9]);
                        add_cs(&mut b[4 * N3 + k], &mut csb[4 * N3 + k], tmp);
                    }
                }
                6 => {
                    for k in 0..N3 {
                        let mut tmp = ias15.g[5 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[5 * N3 + k] = (((((gk / rr[15] - ias15.g[0 * N3 + k]) / rr[16]
                            - ias15.g[1 * N3 + k]) / rr[17]
                            - ias15.g[2 * N3 + k]) / rr[18]
                            - ias15.g[3 * N3 + k]) / rr[19]
                            - ias15.g[4 * N3 + k]) / rr[20];
                        tmp = ias15.g[5 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[10]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp * c[11]);
                        add_cs(&mut b[2 * N3 + k], &mut csb[2 * N3 + k], tmp * c[12]);
                        add_cs(&mut b[3 * N3 + k], &mut csb[3 * N3 + k], tmp * c[13]);
                        add_cs(&mut b[4 * N3 + k], &mut csb[4 * N3 + k], tmp * c[14]);
                        add_cs(&mut b[5 * N3 + k], &mut csb[5 * N3 + k], tmp);
                    }
                }
                7 => {
                    let mut maxak: f64 = 0.0;
                    let mut maxb6ktmp: f64 = 0.0;
                    for k in 0..N3 {
                        let mut tmp = ias15.g[6 * N3 + k];
                        let mut gk = ias15.at[k];
                        let mut gk_cs = gravity_cs_flat(r, &ias15.csa0, use_sim_cs, k);
                        add_cs(&mut gk, &mut gk_cs, -ias15.a0[k]);
                        add_cs(&mut gk, &mut gk_cs, ias15.csa0[k]);
                        ias15.g[6 * N3 + k] = ((((((gk / rr[21] - ias15.g[0 * N3 + k]) / rr[22]
                            - ias15.g[1 * N3 + k]) / rr[23]
                            - ias15.g[2 * N3 + k]) / rr[24]
                            - ias15.g[3 * N3 + k]) / rr[25]
                            - ias15.g[4 * N3 + k]) / rr[26]
                            - ias15.g[5 * N3 + k]) / rr[27];
                        tmp = ias15.g[6 * N3 + k] - tmp;
                        let (b, csb) = (&mut ias15.b, &mut ias15.csb);
                        add_cs(&mut b[0 * N3 + k], &mut csb[0 * N3 + k], tmp * c[15]);
                        add_cs(&mut b[1 * N3 + k], &mut csb[1 * N3 + k], tmp * c[16]);
                        add_cs(&mut b[2 * N3 + k], &mut csb[2 * N3 + k], tmp * c[17]);
                        add_cs(&mut b[3 * N3 + k], &mut csb[3 * N3 + k], tmp * c[18]);
                        add_cs(&mut b[4 * N3 + k], &mut csb[4 * N3 + k], tmp * c[19]);
                        add_cs(&mut b[5 * N3 + k], &mut csb[5 * N3 + k], tmp * c[20]);
                        add_cs(&mut b[6 * N3 + k], &mut csb[6 * N3 + k], tmp);

                        // Monitor change in b6 relative to at. The
                        // predictor corrector scheme is converged if
                        // it is close to 0.
                        if ias15.adaptive_mode != REB_IAS15_ADAPTIVEMODE_INDIVIDUAL {
                            let ak = ias15.at[k].abs();
                            if ak.is_normal() && ak > maxak {
                                maxak = ak;
                            }
                            let b6ktmp = tmp.abs();
                            if b6ktmp.is_normal() && b6ktmp > maxb6ktmp {
                                maxb6ktmp = b6ktmp;
                            }
                        } else {
                            let ak = ias15.at[k];
                            let b6ktmp = tmp;
                            let errork = (b6ktmp / ak).abs();
                            if errork.is_normal() && errork > predictor_corrector_error {
                                predictor_corrector_error = errork;
                            }
                        }
                    }
                    if ias15.adaptive_mode != REB_IAS15_ADAPTIVEMODE_INDIVIDUAL {
                        predictor_corrector_error = maxb6ktmp / maxak;
                    }
                }
                _ => unreachable!(),
            }
        }
    }
    // Set time back to initial value (will be updated below)
    r.t = t_beginning;
    // Find new timestep
    let dt_done = r.dt;

    if ias15.epsilon > 0. {
        let mut dt_new;
        if ias15.adaptive_mode == REB_IAS15_ADAPTIVEMODE_INDIVIDUAL
            || ias15.adaptive_mode == REB_IAS15_ADAPTIVEMODE_GLOBAL
        {
            // Old adaptive timestepping methods
            let mut integrator_error: f64 = 0.0;
            if ias15.adaptive_mode == REB_IAS15_ADAPTIVEMODE_GLOBAL {
                let mut maxa: f64 = 0.0;
                let mut maxj: f64 = 0.0;
                for i in 0..N {
                    let mi = map_index(r, i);
                    let p = r.particles[mi];
                    let v2 = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
                    let x2 = p.x * p.x + p.y * p.y + p.z * p.z;
                    // Skip slowly varying accelerations
                    if (v2 * r.dt * r.dt / x2).abs() < 1e-16 {
                        continue;
                    }
                    for k in (3 * i)..(3 * (i + 1)) {
                        let ak = ias15.at[k].abs();
                        if ak.is_normal() && ak > maxa {
                            maxa = ak;
                        }
                        let b6k = ias15.b[6 * N3 + k].abs();
                        if b6k.is_normal() && b6k > maxj {
                            maxj = b6k;
                        }
                    }
                    integrator_error = maxj / maxa;
                }
            } else {
                for k in 0..N3 {
                    let ak = ias15.at[k];
                    let bk = ias15.b[6 * N3 + k];
                    let errork = (bk / ak).abs();
                    if errork.is_normal() && errork > integrator_error {
                        integrator_error = errork;
                    }
                }
            }
            // Use error estimate to predict new timestep
            if integrator_error.is_normal() {
                dt_new = sqrt7(ias15.epsilon / integrator_error) * dt_done;
            } else {
                // Error estimate is not finite: increase timestep a little
                dt_new = dt_done / safety_factor;
            }
        } else {
            // New adaptive timestepping method (default since Jan 2024)
            let mut min_timescale2 = f64::INFINITY; // factor dt_done^2 not included
            for i in 0..N {
                let mut a0i = 0.; // (acceleration at beginning)^2
                let mut y2 = 0.; // (acceleration at end)^2
                let mut y3 = 0.; // (jerk * dt_done)^2
                let mut y4 = 0.; // (snap * dt_done^2)^2
                let mut y5 = 0.; // (crackle * dt_done^3)^2
                for k in (3 * i)..(3 * (i + 1)) {
                    let b = &ias15.b;
                    a0i += ias15.a0[k] * ias15.a0[k];
                    let mut tmp = ias15.a0[k]
                        + b[0 * N3 + k] + b[1 * N3 + k] + b[2 * N3 + k] + b[3 * N3 + k]
                        + b[4 * N3 + k] + b[5 * N3 + k] + b[6 * N3 + k];
                    y2 += tmp * tmp;
                    tmp = b[0 * N3 + k] + 2. * b[1 * N3 + k] + 3. * b[2 * N3 + k]
                        + 4. * b[3 * N3 + k] + 5. * b[4 * N3 + k] + 6. * b[5 * N3 + k]
                        + 7. * b[6 * N3 + k];
                    y3 += tmp * tmp;
                    tmp = 2. * b[1 * N3 + k] + 6. * b[2 * N3 + k] + 12. * b[3 * N3 + k]
                        + 20. * b[4 * N3 + k] + 30. * b[5 * N3 + k] + 42. * b[6 * N3 + k];
                    y4 += tmp * tmp;
                    tmp = 6. * b[2 * N3 + k] + 24. * b[3 * N3 + k] + 60. * b[4 * N3 + k]
                        + 120. * b[5 * N3 + k] + 210. * b[6 * N3 + k];
                    y5 += tmp * tmp;
                }
                if !a0i.is_normal() {
                    // Skip particles with no or non-finite acceleration
                    continue;
                }
                let mut timescale2 = 0.;
                if ias15.adaptive_mode == REB_IAS15_ADAPTIVEMODE_PRS23 {
                    timescale2 = 2. * y2 / (y3 + (y4 * y2).sqrt()); // PRS23
                } else if ias15.adaptive_mode == REB_IAS15_ADAPTIVEMODE_AARSETH85 {
                    timescale2 = ((y2 * y4).sqrt() + y3) / ((y3 * y5).sqrt() + y4); // A85
                }

                if timescale2.is_normal() && timescale2 < min_timescale2 {
                    min_timescale2 = timescale2;
                }
            }
            if min_timescale2.is_normal() {
                // Numerical factor matches adaptive_mode GLOBAL with default epsilon
                dt_new = min_timescale2.sqrt() * dt_done * sqrt7(ias15.epsilon * 5040.0);
            } else {
                dt_new = dt_done / safety_factor; // increase timestep a little
            }
        }

        if dt_new.abs() < ias15.min_dt {
            dt_new = ias15.min_dt.copysign(dt_new);
        }

        if (dt_new / dt_done).abs() < safety_factor {
            // New timestep is significantly smaller. Reset particles.
            for k in 0..N {
                let mk = map_index(r, k);
                r.particles[mk].x = ias15.x0[3 * k];
                r.particles[mk].y = ias15.x0[3 * k + 1];
                r.particles[mk].z = ias15.x0[3 * k + 2];
                r.particles[mk].vx = ias15.v0[3 * k];
                r.particles[mk].vy = ias15.v0[3 * k + 1];
                r.particles[mk].vz = ias15.v0[3 * k + 2];
                r.particles[mk].ax = ias15.a0[3 * k];
                r.particles[mk].ay = ias15.a0[3 * k + 1];
                r.particles[mk].az = ias15.a0[3 * k + 2];
            }
            for mk in 0..N_var {
                let k = mk + N;
                r.particles_var[mk].x = ias15.x0[3 * k];
                r.particles_var[mk].y = ias15.x0[3 * k + 1];
                r.particles_var[mk].z = ias15.x0[3 * k + 2];
                r.particles_var[mk].vx = ias15.v0[3 * k];
                r.particles_var[mk].vy = ias15.v0[3 * k + 1];
                r.particles_var[mk].vz = ias15.v0[3 * k + 2];
                r.particles_var[mk].ax = ias15.a0[3 * k];
                r.particles_var[mk].ay = ias15.a0[3 * k + 1];
                r.particles_var[mk].az = ias15.a0[3 * k + 2];
            }
            r.dt = dt_new;
            if r.dt_last_done != 0. {
                // Do not predict next e/b values on the very first step
                let ratio = r.dt / r.dt_last_done;
                let (er, br) = (ias15.er.clone(), ias15.br.clone());
                predict_next_step(ratio, N3, &er, &br, &mut ias15.e, &mut ias15.b);
            }
            return false; // Step rejected. Do again.
        }
        if (dt_new / dt_done).abs() > 1.0 {
            // New timestep is larger.
            if dt_new / dt_done > 1. / safety_factor {
                // Don't increase the timestep by too much
                dt_new = dt_done / safety_factor;
            }
        }
        r.dt = dt_new;
    }

    // Find new position and velocity values at end of the sequence
    for k in 0..N3 {
        // Note: dt_done*dt_done is not precalculated to avoid biased
        // round-off errors when a fixed timestep is used.
        let b = &ias15.b;
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[6 * N3 + k] / 72. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[5 * N3 + k] / 56. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[4 * N3 + k] / 42. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[3 * N3 + k] / 30. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[2 * N3 + k] / 20. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[1 * N3 + k] / 12. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], b[0 * N3 + k] / 6. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], ias15.a0[k] / 2. * dt_done * dt_done);
        add_cs(&mut ias15.x0[k], &mut ias15.csx[k], ias15.v0[k] * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[6 * N3 + k] / 8. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[5 * N3 + k] / 7. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[4 * N3 + k] / 6. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[3 * N3 + k] / 5. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[2 * N3 + k] / 4. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[1 * N3 + k] / 3. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], b[0 * N3 + k] / 2. * dt_done);
        add_cs(&mut ias15.v0[k], &mut ias15.csv[k], ias15.a0[k] * dt_done);
    }

    r.t += dt_done;
    r.dt_last_done = dt_done;

    if r.calculate_megno != 0 {
        let dY = dt_done * integrator_megno_thisdt;
        reb_tools_megno_update(r, dY, dt_done);
    }

    // Swap particle buffers
    for k in 0..N {
        let mk = map_index(r, k);
        r.particles[mk].x = ias15.x0[3 * k];
        r.particles[mk].y = ias15.x0[3 * k + 1];
        r.particles[mk].z = ias15.x0[3 * k + 2];
        r.particles[mk].vx = ias15.v0[3 * k];
        r.particles[mk].vy = ias15.v0[3 * k + 1];
        r.particles[mk].vz = ias15.v0[3 * k + 2];
    }
    for mk in 0..N_var {
        let k = mk + N;
        r.particles_var[mk].x = ias15.x0[3 * k];
        r.particles_var[mk].y = ias15.x0[3 * k + 1];
        r.particles_var[mk].z = ias15.x0[3 * k + 2];
        r.particles_var[mk].vx = ias15.v0[3 * k];
        r.particles_var[mk].vy = ias15.v0[3 * k + 1];
        r.particles_var[mk].vz = ias15.v0[3 * k + 2];
    }
    copybuffers(&ias15.e, &mut ias15.er, N3);
    copybuffers(&ias15.b, &mut ias15.br, N3);
    let ratio = r.dt / dt_done;
    predict_next_step_inplace(ratio, N3, &mut ias15.e, &mut ias15.b);
    true // Success.
}

/// integrator_ias15.c `reb_integrator_ias15_step` (state-explicit; the
/// C `.step` callback takes `void* state`). MERCURIUS drives a private
/// IAS15 instance through this during close encounters.
pub fn reb_integrator_ias15_step_state(
    r: &mut reb_simulation,
    ias15: &mut reb_integrator_ias15_state,
) {
    r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    if r.N != 0 {
        // Try until a step was successful.
        while !reb_integrator_ias15_step_try(r, ias15) {}
    } else {
        r.t += r.dt;
        r.dt_last_done = r.dt;
    }
}

/// Step entry point for the dispatcher: takes the state out of the enum
/// for the duration of the step (the C passes `void* state` alongside
/// `r`; Rust needs the aliasing made explicit).
pub fn reb_integrator_ias15_step(r: &mut reb_simulation) {
    if r.N != 0 {
        let mut ias15 = match std::mem::replace(&mut r.integrator, reb_integrator_state::none) {
            reb_integrator_state::ias15(s) => s,
            other => {
                r.integrator = other;
                return;
            }
        };
        reb_integrator_ias15_step_state(r, &mut ias15);
        r.integrator = reb_integrator_state::ias15(ias15);
    } else {
        r.gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
        r.t += r.dt;
        r.dt_last_done = r.dt;
    }
}

/// integrator_ias15.c `reb_integrator_ias15_timescale` (PRS23).
pub fn reb_integrator_ias15_timescale(r: &mut reb_simulation) -> f64 {
    reb_simulation_update_acceleration(r);
    let N = if r.map.is_some() { r.N_map } else { r.N };

    let mut min_timescale2 = f64::INFINITY;

    for i in 0..N {
        let mi = map_index(r, i);
        let p_i = r.particles[mi];
        let y2 = p_i.ax * p_i.ax + p_i.ay * p_i.ay + p_i.az * p_i.az;
        let mut vec_y3 = reb_vec3d::default();
        let mut vec_y4 = reb_vec3d::default();

        if !y2.is_normal() {
            continue;
        }
        for j in 0..N {
            let mj = map_index(r, j);
            if mi == mj {
                continue;
            }
            let p_j = r.particles[mj];

            let rij_x = p_j.x - p_i.x;
            let rij_y = p_j.y - p_i.y;
            let rij_z = p_j.z - p_i.z;
            let vij_x = p_j.vx - p_i.vx;
            let vij_y = p_j.vy - p_i.vy;
            let vij_z = p_j.vz - p_i.vz;
            let aij_x = p_j.ax - p_i.ax;
            let aij_y = p_j.ay - p_i.ay;
            let aij_z = p_j.az - p_i.az;

            let r_sq = rij_x * rij_x + rij_y * rij_y + rij_z * rij_z;

            let r_mag = r_sq.sqrt();
            let r_cubed = r_sq * r_mag;
            let r_fifth = r_cubed * r_sq;
            let r_seventh = r_fifth * r_sq;

            let r_dot_v = rij_x * vij_x + rij_y * vij_y + rij_z * vij_z;
            let r_dot_a = rij_x * aij_x + rij_y * aij_y + rij_z * aij_z;
            let v_sq = vij_x * vij_x + vij_y * vij_y + vij_z * vij_z;

            // Jerk
            let jerk_factor1 = p_j.m / r_cubed;
            let jerk_factor2 = -3.0 * p_j.m * r_dot_v / r_fifth;

            vec_y3.x += jerk_factor1 * vij_x + jerk_factor2 * rij_x;
            vec_y3.y += jerk_factor1 * vij_y + jerk_factor2 * rij_y;
            vec_y3.z += jerk_factor1 * vij_z + jerk_factor2 * rij_z;

            // Snap
            let snap_c1 = p_j.m / r_cubed;
            let snap_c2 = -6.0 * p_j.m * r_dot_v / r_fifth;
            let snap_c3_rij = -3.0 * p_j.m * v_sq / r_fifth;
            let snap_c4_rij = -3.0 * p_j.m * r_dot_a / r_fifth;
            let snap_c5_rij = 15.0 * p_j.m * r_dot_v * r_dot_v / r_seventh;

            vec_y4.x += snap_c1 * aij_x + snap_c2 * vij_x + (snap_c3_rij + snap_c4_rij + snap_c5_rij) * rij_x;
            vec_y4.y += snap_c1 * aij_y + snap_c2 * vij_y + (snap_c3_rij + snap_c4_rij + snap_c5_rij) * rij_y;
            vec_y4.z += snap_c1 * aij_z + snap_c2 * vij_z + (snap_c3_rij + snap_c4_rij + snap_c5_rij) * rij_z;
        }
        vec_y3.x *= r.G;
        vec_y3.y *= r.G;
        vec_y3.z *= r.G;
        vec_y4.x *= r.G;
        vec_y4.y *= r.G;
        vec_y4.z *= r.G;
        let y3 = vec_y3.x * vec_y3.x + vec_y3.y * vec_y3.y + vec_y3.z * vec_y3.z;
        let y4 = vec_y4.x * vec_y4.x + vec_y4.y * vec_y4.y + vec_y4.z * vec_y4.z;
        let timescale2 = 2. * y2 / (y3 + (y4 * y2).sqrt()); // PRS23
        if timescale2.is_normal() && timescale2 < min_timescale2 {
            min_timescale2 = timescale2;
        }
    }
    min_timescale2.sqrt()
}
