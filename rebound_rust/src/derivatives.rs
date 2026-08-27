//! derivatives.rs — functions to calculate derivatives of Keplerian orbits
//! (from derivatives.c; (c) 2016 Hanno Rein, Dan Tamayp, Rejean Leblanc).
//! Arithmetic order matches the C source expression by expression.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Daniel Tamayo and contributors. See crate root.

use crate::tools::{reb_orbit_from_particle, reb_tools_particle_to_pal, reb_tools_solve_kepler_pal};
use crate::types::*;

pub fn reb_particle_derivative_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_dlambda = a * (dclp_dlambda + dp_dlambda / (2. - l) * h);
    let deta_dlambda = a * (dslp_dlambda - dp_dlambda / (2. - l) * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dlambda = deta_dlambda * ix - dxi_dlambda * iy;

    np.x = dxi_dlambda + 0.5 * iy * dW_dlambda;
    np.y = deta_dlambda - 0.5 * ix * dW_dlambda;
    np.z = 0.5 * iz * dW_dlambda;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dlambda + dq_dlambda / (2. - l) * h);
    let ddeta_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dlambda - dq_dlambda / (2. - l) * k);
    let ddW_dlambda = ddeta_dlambda * ix - ddxi_dlambda * iy;
    np.vx = ddxi_dlambda + 0.5 * iy * ddW_dlambda;
    np.vy = ddeta_dlambda - 0.5 * ix * ddW_dlambda;
    np.vz = 0.5 * iz * ddW_dlambda;

    np
}

pub fn reb_particle_derivative_h(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dxi_dh = a * (dclp_dh + dp_dh / (2. - l) * h + p / (2. - l) + p / ((2. - l) * (2. - l)) * dl_dh * h);
    let deta_dh = a * (dslp_dh - dp_dh / (2. - l) * k - p / ((2. - l) * (2. - l)) * k * dl_dh - 1.);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dh = deta_dh * ix - dxi_dh * iy;

    np.x = dxi_dh + 0.5 * iy * dW_dh;
    np.y = deta_dh - 0.5 * ix * dW_dh;
    np.z = 0.5 * iz * dW_dh;

    let dq_dh = 1. / (1. - q) * (slp - h);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dh = dq_dh * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l));
    let ddeta_dh = dq_dh * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k);
    let ddW_dh = ddeta_dh * ix - ddxi_dh * iy;

    np.vx = ddxi_dh + 0.5 * iy * ddW_dh;
    np.vy = ddeta_dh - 0.5 * ix * ddW_dh;
    np.vz = 0.5 * iz * ddW_dh;

    np
}

pub fn reb_particle_derivative_k(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dxi_dk = a * (dclp_dk + dp_dk / (2. - l) * h + p / ((2. - l) * (2. - l)) * dl_dk * h - 1.);
    let deta_dk = a * (dslp_dk - dp_dk / (2. - l) * k - p / (2. - l) - p / ((2. - l) * (2. - l)) * dl_dk * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dk = deta_dk * ix - dxi_dk * iy;

    np.x = dxi_dk + 0.5 * iy * dW_dk;
    np.y = deta_dk - 0.5 * ix * dW_dk;
    np.z = 0.5 * iz * dW_dk;

    let dq_dk = 1. / (1. - q) * (clp - k);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dk = dq_dk * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h);
    let ddeta_dk = dq_dk * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l));
    let ddW_dk = ddeta_dk * ix - ddxi_dk * iy;

    np.vx = ddxi_dk + 0.5 * iy * ddW_dk;
    np.vy = ddeta_dk - 0.5 * ix * ddW_dk;
    np.vz = 0.5 * iz * ddW_dk;

    np
}

pub fn reb_particle_derivative_k_k(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dl_dkk = 1. / (1. - h * h - k * k).sqrt()
        + (k * k) / ((1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt());
    let dp_dk = 1. / (1. - q) * (slp);
    let dq_dk = 1. / (1. - q) * (clp - k);
    let dp_dkk = dq_dk / ((1. - q) * (1. - q)) * (slp) + 1. / (1. - q) * (dslp_dk);
    let dq_dkk = dq_dk / ((1. - q) * (1. - q)) * (clp - k) + 1. / (1. - q) * (dclp_dk - 1.);
    let dclp_dkk = -dq_dk / ((1. - q) * (1. - q)) * (slp * slp) - 2. / (1. - q) * slp * dslp_dk;
    let dslp_dkk = -dq_dk / ((1. - q) * (1. - q)) * (-slp * clp) - 1. / (1. - q) * -slp * dclp_dk - 1. / (1. - q) * -dslp_dk * clp;

    let dxi_dkk = a * (dclp_dkk + dp_dkk / (2. - l) * h + dl_dk * dp_dk / ((2. - l) * (2. - l)) * h + dp_dk / ((2. - l) * (2. - l)) * dl_dk * h + 2. * dl_dk * p / ((2. - l) * (2. - l) * (2. - l)) * dl_dk * h + p / ((2. - l) * (2. - l)) * dl_dkk * h);
    let deta_dkk = a * (dslp_dkk - dp_dkk / (2. - l) * k - dl_dk * dp_dk / ((2. - l) * (2. - l)) * k - dp_dk / (2. - l) - dp_dk / (2. - l) - dl_dk * p / ((2. - l) * (2. - l))
        - dp_dk / ((2. - l) * (2. - l)) * dl_dk * k - 2. * dl_dk * p / ((2. - l) * (2. - l) * (2. - l)) * dl_dk * k - p / ((2. - l) * (2. - l)) * dl_dkk * k - p / ((2. - l) * (2. - l)) * dl_dk);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dkk = deta_dkk * ix - dxi_dkk * iy;

    np.x = dxi_dkk + 0.5 * iy * dW_dkk;
    np.y = deta_dkk - 0.5 * ix * dW_dkk;
    np.z = 0.5 * iz * dW_dkk;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dkk = dq_dkk * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h) + 2. * dq_dk * dq_dk * an / ((1. - q) * (1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dq_dk * an / ((1. - q) * (1. - q)) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h)
        + dq_dk * an / ((1. - q) * (1. - q)) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h)
        + an / (1. - q) * (-dslp_dkk + dq_dkk / (2. - l) * h + dl_dk * dq_dk / ((2. - l) * (2. - l)) * h
            + dl_dkk * q / ((2. - l) * (2. - l)) * h + dl_dk * dq_dk / ((2. - l) * (2. - l)) * h + 2. * dl_dk * dl_dk * q / ((2. - l) * (2. - l) * (2. - l)) * h);
    let ddeta_dkk = dq_dkk * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k) + 2. * dq_dk * dq_dk * an / ((1. - q) * (1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dq_dk * an / ((1. - q) * (1. - q)) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l))
        + dq_dk * an / ((1. - q) * (1. - q)) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l))
        + an / (1. - q) * (dclp_dkk - dq_dkk / (2. - l) * k - dq_dk * dl_dk / ((2. - l) * (2. - l)) * k - dq_dk / (2. - l)
            - dl_dkk * q / ((2. - l) * (2. - l)) * k - dl_dk * dq_dk / ((2. - l) * (2. - l)) * k - 2. * dl_dk * dl_dk * q / ((2. - l) * (2. - l) * (2. - l)) * k - dl_dk * q / ((2. - l) * (2. - l)) - dq_dk / (2. - l) - dl_dk * q / ((2. - l) * (2. - l)));
    let ddW_dkk = ddeta_dkk * ix - ddxi_dkk * iy;

    np.vx = ddxi_dkk + 0.5 * iy * ddW_dkk;
    np.vy = ddeta_dkk - 0.5 * ix * ddW_dkk;
    np.vz = 0.5 * iz * ddW_dkk;

    np
}

pub fn reb_particle_derivative_h_h(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dl_dhh = 1. / (1. - h * h - k * k).sqrt()
        + (h * h) / ((1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt());
    let dp_dh = 1. / (1. - q) * (-clp);
    let dq_dh = 1. / (1. - q) * (slp - h);
    let dq_dhh = 1. / ((1. - q) * (1. - q)) * dq_dh * (slp - h) + 1. / (1. - q) * (dslp_dh - 1.);
    let dp_dhh = 1. / ((1. - q) * (1. - q)) * dq_dh * (-clp) + 1. / (1. - q) * (-dclp_dh);
    let dclp_dhh = -1. / ((1. - q) * (1. - q)) * dq_dh * (-slp * clp) - 1. / (1. - q) * (-dslp_dh * clp) - 1. / (1. - q) * (-slp * dclp_dh);
    let dslp_dhh = -1. / ((1. - q) * (1. - q)) * dq_dh * (clp * clp) - 2. / (1. - q) * (clp * dclp_dh);

    let dxi_dhh = a * (dclp_dhh + (dp_dhh / (2. - l) * h + dl_dh * dp_dh / ((2. - l) * (2. - l)) * h + dp_dh / (2. - l)) + (dp_dh / (2. - l) + dl_dh * p / ((2. - l) * (2. - l)))
        + (dp_dh / ((2. - l) * (2. - l)) * dl_dh * h + 2. * p / ((2. - l) * (2. - l) * (2. - l)) * dl_dh * dl_dh * h + p / ((2. - l) * (2. - l)) * dl_dhh * h + p / ((2. - l) * (2. - l)) * dl_dh));
    let deta_dhh = a * (dslp_dhh + (-dp_dhh / (2. - l) * k - dl_dh * dp_dh / ((2. - l) * (2. - l)) * k) + (-dp_dh / ((2. - l) * (2. - l)) * k * dl_dh - 2. * p / ((2. - l) * (2. - l) * (2. - l)) * k * dl_dh * dl_dh - p / ((2. - l) * (2. - l)) * k * dl_dhh));

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dhh = deta_dhh * ix - dxi_dhh * iy;

    np.x = dxi_dhh + 0.5 * iy * dW_dhh;
    np.y = deta_dhh - 0.5 * ix * dW_dhh;
    np.z = 0.5 * iz * dW_dhh;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dhh = dq_dhh * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h) + 2. * dq_dh * dq_dh * an / ((1. - q) * (1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dq_dh * an / ((1. - q) * (1. - q)) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l))
        + dq_dh * an / ((1. - q) * (1. - q)) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l))
        + an / (1. - q) * (-dslp_dhh + (dq_dhh / (2. - l) * h + dl_dh * dq_dh / ((2. - l) * (2. - l)) * h + dq_dh / (2. - l))
            + (dl_dhh * q / ((2. - l) * (2. - l)) * h + dl_dh * dq_dh / ((2. - l) * (2. - l)) * h + 2. * dl_dh * dl_dh * q / ((2. - l) * (2. - l) * (2. - l)) * h + dl_dh * q / ((2. - l) * (2. - l))) + (dq_dh / (2. - l) + dl_dh * q / ((2. - l) * (2. - l))));
    let ddeta_dhh = dq_dhh * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k) + 2. * dq_dh * dq_dh * an / ((1. - q) * (1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dq_dh * an / ((1. - q) * (1. - q)) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k)
        + dq_dh * an / ((1. - q) * (1. - q)) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k)
        + an / (1. - q) * (dclp_dhh - dq_dhh / (2. - l) * k - dl_dh * dq_dh / ((2. - l) * (2. - l)) * k
            - dl_dhh * q / ((2. - l) * (2. - l)) * k - dl_dh * dq_dh / ((2. - l) * (2. - l)) * k - 2. * dl_dh * dl_dh * q / ((2. - l) * (2. - l) * (2. - l)) * k);

    let ddW_dhh = ddeta_dhh * ix - ddxi_dhh * iy;

    np.vx = ddxi_dhh + 0.5 * iy * ddW_dhh;
    np.vy = ddeta_dhh - 0.5 * ix * ddW_dhh;
    np.vz = 0.5 * iz * ddW_dhh;

    np
}

pub fn reb_particle_derivative_lambda_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);
    let dq_dlambdalambda = -dp_dlambda / (1. - q) - p / ((1. - q) * (1. - q)) * dq_dlambda;
    let dp_dlambdalambda = dq_dlambda / (1. - q) + q / ((1. - q) * (1. - q)) * dq_dlambda;

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;
    let dclp_dlambdalambda = -1. / ((1. - q) * (1. - q)) * dq_dlambda * slp - 1. / (1. - q) * dslp_dlambda;
    let dslp_dlambdalambda = 1. / ((1. - q) * (1. - q)) * dq_dlambda * clp + 1. / (1. - q) * dclp_dlambda;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_dlambdalambda = a * (dclp_dlambdalambda + dp_dlambdalambda / (2. - l) * h);
    let deta_dlambdalambda = a * (dslp_dlambdalambda - dp_dlambdalambda / (2. - l) * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dlambdalambda = deta_dlambdalambda * ix - dxi_dlambdalambda * iy;

    np.x = dxi_dlambdalambda + 0.5 * iy * dW_dlambdalambda;
    np.y = deta_dlambdalambda - 0.5 * ix * dW_dlambdalambda;
    np.z = 0.5 * iz * dW_dlambdalambda;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dlambdalambda = 2. * an / ((1. - q) * (1. - q) * (1. - q)) * dq_dlambda * dq_dlambda * (-slp + q / (2. - l) * h)
        + an / ((1. - q) * (1. - q)) * dq_dlambdalambda * (-slp + q / (2. - l) * h) + an / ((1. - q) * (1. - q)) * dq_dlambda * (-dslp_dlambda + dq_dlambda / (2. - l) * h)
        + an / ((1. - q) * (1. - q)) * dq_dlambda * (-dslp_dlambda + dq_dlambda / (2. - l) * h) + an / (1. - q) * (-dslp_dlambdalambda + dq_dlambdalambda / (2. - l) * h);
    let ddeta_dlambdalambda = 2. * an / ((1. - q) * (1. - q) * (1. - q)) * dq_dlambda * dq_dlambda * (clp - q / (2. - l) * k)
        + an / ((1. - q) * (1. - q)) * dq_dlambdalambda * (clp - q / (2. - l) * k) + an / ((1. - q) * (1. - q)) * dq_dlambda * (dclp_dlambda - dq_dlambda / (2. - l) * k)
        + an / ((1. - q) * (1. - q)) * dq_dlambda * (dclp_dlambda - dq_dlambda / (2. - l) * k) + an / (1. - q) * (dclp_dlambdalambda - dq_dlambdalambda / (2. - l) * k);

    let ddW_dlambdalambda = ddeta_dlambdalambda * ix - ddxi_dlambdalambda * iy;
    np.vx = ddxi_dlambdalambda + 0.5 * iy * ddW_dlambdalambda;
    np.vy = ddeta_dlambdalambda - 0.5 * ix * ddW_dlambdalambda;
    np.vz = 0.5 * iz * ddW_dlambdalambda;

    np
}

pub fn reb_particle_derivative_k_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dq_dk = 1. / (1. - q) * (clp - k);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);
    let dq_dklambda = -dp_dk / (1. - q) - p / ((1. - q) * (1. - q)) * dq_dk;
    let dp_dklambda = dq_dk / (1. - q) + q / ((1. - q) * (1. - q)) * dq_dk;
    let dclp_dklambda = -1. / (1. - q) * dslp_dk - 1. / ((1. - q) * (1. - q)) * dq_dk * slp;
    let dslp_dklambda = 1. / (1. - q) * dclp_dk + 1. / ((1. - q) * (1. - q)) * dq_dk * clp;

    let dxi_dklambda = a * (dclp_dklambda + dp_dklambda / (2. - l) * h + dp_dlambda / ((2. - l) * (2. - l)) * dl_dk * h);
    let deta_dklambda = a * (dslp_dklambda - dp_dklambda / (2. - l) * k - dp_dlambda / (2. - l) - dp_dlambda / ((2. - l) * (2. - l)) * dl_dk * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dklambda = deta_dklambda * ix - dxi_dklambda * iy;

    np.x = dxi_dklambda + 0.5 * iy * dW_dklambda;
    np.y = deta_dklambda - 0.5 * ix * dW_dklambda;
    np.z = 0.5 * iz * dW_dklambda;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dklambda = dq_dklambda * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h) + 2. * dq_dk * dq_dlambda * an / ((1. - q) * (1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dq_dk * an / ((1. - q) * (1. - q)) * (-dslp_dlambda + dq_dlambda / (2. - l) * h)
        + dq_dlambda * an / ((1. - q) * (1. - q)) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h)
        + an / (1. - q) * (-dslp_dklambda + dq_dklambda / (2. - l) * h + dl_dk * dq_dlambda / ((2. - l) * (2. - l)) * h);
    let ddeta_dklambda = dq_dklambda * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k) + 2. * dq_dk * dq_dlambda * an / ((1. - q) * (1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dq_dk * an / ((1. - q) * (1. - q)) * (dclp_dlambda - dq_dlambda / (2. - l) * k)
        + dq_dlambda * an / ((1. - q) * (1. - q)) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l))
        + an / (1. - q) * (dclp_dklambda - dq_dklambda / (2. - l) * k - dl_dk * dq_dlambda / ((2. - l) * (2. - l)) * k - dq_dlambda / (2. - l));
    let ddW_dklambda = ddeta_dklambda * ix - ddxi_dklambda * iy;

    np.vx = ddxi_dklambda + 0.5 * iy * ddW_dklambda;
    np.vy = ddeta_dklambda - 0.5 * ix * ddW_dklambda;
    np.vz = 0.5 * iz * ddW_dklambda;

    np
}

pub fn reb_particle_derivative_h_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dq_dh = 1. / (1. - q) * (slp - h);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);
    let dq_dhlambda = -dp_dh / (1. - q) - p / ((1. - q) * (1. - q)) * dq_dh;
    let dp_dhlambda = dq_dh / (1. - q) + q / ((1. - q) * (1. - q)) * dq_dh;
    let dclp_dhlambda = -1. / ((1. - q) * (1. - q)) * (-slp * clp) * dq_dlambda - 1. / (1. - q) * (-dslp_dlambda * clp) - 1. / (1. - q) * (-slp * dclp_dlambda);
    let dslp_dhlambda = -1. / ((1. - q) * (1. - q)) * (clp * clp) * dq_dlambda - 2. / (1. - q) * (clp * dclp_dlambda);

    let dxi_dhlambda = a * (dclp_dhlambda + dp_dhlambda / (2. - l) * h + dp_dlambda / (2. - l) + dp_dlambda / ((2. - l) * (2. - l)) * dl_dh * h);
    let deta_dhlambda = a * (dslp_dhlambda - dp_dhlambda / (2. - l) * k - dp_dlambda / ((2. - l) * (2. - l)) * k * dl_dh);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dhlambda = deta_dhlambda * ix - dxi_dhlambda * iy;

    np.x = dxi_dhlambda + 0.5 * iy * dW_dhlambda;
    np.y = deta_dhlambda - 0.5 * ix * dW_dhlambda;
    np.z = 0.5 * iz * dW_dhlambda;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dhlambda = dq_dhlambda * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h) + 2. * dq_dlambda * dq_dh * an / ((1. - q) * (1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dq_dh * an / ((1. - q) * (1. - q)) * (-dslp_dlambda + dq_dlambda / (2. - l) * h)
        + dq_dlambda * an / ((1. - q) * (1. - q)) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l))
        + an / (1. - q) * (-dslp_dhlambda + dq_dhlambda / (2. - l) * h + dl_dh * dq_dlambda / ((2. - l) * (2. - l)) * h + dq_dlambda / (2. - l));
    let ddeta_dhlambda = dq_dhlambda * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k) + 2. * dq_dh * dq_dlambda * an / ((1. - q) * (1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dq_dh * an / ((1. - q) * (1. - q)) * (dclp_dlambda - dq_dlambda / (2. - l) * k)
        + dq_dlambda * an / ((1. - q) * (1. - q)) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k)
        + an / (1. - q) * (dclp_dhlambda - dq_dhlambda / (2. - l) * k - dl_dh * dq_dlambda / ((2. - l) * (2. - l)) * k);
    let ddW_dhlambda = ddeta_dhlambda * ix - ddxi_dhlambda * iy;

    np.vx = ddxi_dhlambda + 0.5 * iy * ddW_dhlambda;
    np.vy = ddeta_dhlambda - 0.5 * ix * ddW_dhlambda;
    np.vz = 0.5 * iz * ddW_dhlambda;

    np
}

pub fn reb_particle_derivative_k_h(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dq_dh = 1. / (1. - q) * (slp - h);
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dq_dk = 1. / (1. - q) * (clp - k);
    let dl_dkh = k * h / ((1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt() * (1. - h * h - k * k).sqrt());
    let dp_dkh = 1. / ((1. - q) * (1. - q)) * dq_dh * (slp) + 1. / (1. - q) * (dslp_dh);
    let dq_dkh = 1. / ((1. - q) * (1. - q)) * dq_dh * (clp - k) + 1. / (1. - q) * (dclp_dh);
    let dclp_dkh = -1. / ((1. - q) * (1. - q)) * dq_dh * (slp * slp) - 2. / (1. - q) * (slp * dslp_dh);
    let dslp_dkh = -1. / ((1. - q) * (1. - q)) * dq_dh * (-slp * clp) - 1. / (1. - q) * (-dslp_dh * clp) - 1. / (1. - q) * (-slp * dclp_dh);

    let dxi_dkh = a * (dclp_dkh + dp_dkh / (2. - l) * h + dl_dh * dp_dk / ((2. - l) * (2. - l)) * h + dp_dk / (2. - l)
        + dp_dh / ((2. - l) * (2. - l)) * dl_dk * h + 2. * p / ((2. - l) * (2. - l) * (2. - l)) * dl_dk * dl_dh * h + p / ((2. - l) * (2. - l)) * dl_dkh * h + p / ((2. - l) * (2. - l)) * dl_dk);
    let deta_dkh = a * (dslp_dkh - dp_dkh / (2. - l) * k - dl_dh * dp_dk / ((2. - l) * (2. - l)) * k - dp_dh / (2. - l) - dl_dh * p / ((2. - l) * (2. - l))
        - dp_dh / ((2. - l) * (2. - l)) * dl_dk * k - p / ((2. - l) * (2. - l)) * dl_dkh * k - 2. * p / ((2. - l) * (2. - l) * (2. - l)) * dl_dk * dl_dh * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dkh = deta_dkh * ix - dxi_dkh * iy;

    np.x = dxi_dkh + 0.5 * iy * dW_dkh;
    np.y = deta_dkh - 0.5 * ix * dW_dkh;
    np.z = 0.5 * iz * dW_dkh;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dkh = dq_dkh * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h) + 2. * dq_dh * dq_dk * an / ((1. - q) * (1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dq_dk * an / ((1. - q) * (1. - q)) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l))
        + dq_dh * an / ((1. - q) * (1. - q)) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h)
        + an / (1. - q) * (-dslp_dkh + (dq_dkh / (2. - l) * h + dl_dh * dq_dk / ((2. - l) * (2. - l)) * h + dq_dk / (2. - l))
            + dl_dkh * q / ((2. - l) * (2. - l)) * h + dl_dk * dq_dh / ((2. - l) * (2. - l)) * h + 2. * dl_dh * dl_dk * q / ((2. - l) * (2. - l) * (2. - l)) * h + dl_dk * q / ((2. - l) * (2. - l)));
    let ddeta_dkh = dq_dkh * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k) + 2. * dq_dh * dq_dk * an / ((1. - q) * (1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dq_dk * an / ((1. - q) * (1. - q)) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k)
        + dq_dh * an / ((1. - q) * (1. - q)) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l))
        + an / (1. - q) * (dclp_dkh - dq_dkh / (2. - l) * k - dl_dh * dq_dk / ((2. - l) * (2. - l)) * k
            - dl_dkh * q / ((2. - l) * (2. - l)) * k - dl_dk * dq_dh / ((2. - l) * (2. - l)) * k - 2. * dl_dk * dl_dh * q / ((2. - l) * (2. - l) * (2. - l)) * k - dq_dh / (2. - l) - dl_dh * q / ((2. - l) * (2. - l)));
    let ddW_dkh = ddeta_dkh * ix - ddxi_dkh * iy;

    np.vx = ddxi_dkh + 0.5 * iy * ddW_dkh;
    np.vy = ddeta_dkh - 0.5 * ix * ddW_dkh;
    np.vz = 0.5 * iz * ddW_dkh;

    np
}

pub fn reb_particle_derivative_a(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_da = clp + p / (2. - l) * h - k;
    let deta_da = slp - p / (2. - l) * k - h;

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_da = deta_da * ix - dxi_da * iy;

    np.x = dxi_da + 0.5 * iy * dW_da;
    np.y = deta_da - 0.5 * ix * dW_da;
    np.z = 0.5 * iz * dW_da;

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddxi_da = dan_da / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_da = dan_da / (1. - q) * (clp - q / (2. - l) * k);

    let ddW_da = ddeta_da * ix - ddxi_da * iy;
    np.vx = ddxi_da + 0.5 * iy * ddW_da;
    np.vy = ddeta_da - 0.5 * ix * ddW_da;
    np.vz = 0.5 * iz * ddW_da;

    np
}

pub fn reb_particle_derivative_a_a(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_daa = 0.0; //clp + p/(2.-l)*h -k;
    let deta_daa = 0.0; //slp - p/(2.-l)*k -h;

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_daa = deta_daa * ix - dxi_daa * iy;

    np.x = dxi_daa + 0.5 * iy * dW_daa;
    np.y = deta_daa - 0.5 * ix * dW_daa;
    np.z = 0.5 * iz * dW_daa;

    let dan_daa = 0.75 * (G * (po.m + primary.m) / (a * a * a * a * a)).sqrt();
    let ddxi_daa = dan_daa / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_daa = dan_daa / (1. - q) * (clp - q / (2. - l) * k);

    let ddW_daa = ddeta_daa * ix - ddxi_daa * iy;
    np.vx = ddxi_daa + 0.5 * iy * ddW_daa;
    np.vy = ddeta_daa - 0.5 * ix * ddW_daa;
    np.vz = 0.5 * iz * ddW_daa;

    np
}

pub fn reb_particle_derivative_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let xi = a * (clp + p / (2. - l) * h - k);
    let eta = a * (slp - p / (2. - l) * k - h);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let W = eta * ix - xi * iy;
    let dW_dix = eta;

    np.x = 0.5 * iy * dW_dix;
    np.y = -0.5 * W - 0.5 * ix * dW_dix;
    np.z = 0.5 * diz_dix * W + 0.5 * iz * dW_dix;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let dxi = an / (1. - q) * (-slp + q / (2. - l) * h);
    let deta = an / (1. - q) * (clp - q / (2. - l) * k);
    let dW = deta * ix - dxi * iy;
    let ddW_dix = deta;

    np.vx = 0.5 * iy * ddW_dix;
    np.vy = -0.5 * dW - 0.5 * ix * ddW_dix;
    np.vz = 0.5 * diz_dix * dW + 0.5 * iz * ddW_dix;

    np
}

pub fn reb_particle_derivative_ix_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let xi = a * (clp + p / (2. - l) * h - k);
    let eta = a * (slp - p / (2. - l) * k - h);

    let iz = (4. - ix * ix - iy * iy).sqrt();
    let diz_dix = -ix / iz;
    let diz_dixix = -1. / iz - ix * ix / (iz * iz * iz);
    let W = eta * ix - xi * iy;
    let dW_dix = eta;
    let dW_dixix = 0.0;

    np.x = 0.5 * iy * dW_dixix;
    np.y = -dW_dix - 0.5 * ix * dW_dixix;
    np.z = 0.5 * diz_dixix * W + diz_dix * dW_dix;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let dxi = an / (1. - q) * (-slp + q / (2. - l) * h);
    let deta = an / (1. - q) * (clp - q / (2. - l) * k);
    let dW = deta * ix - dxi * iy;
    let ddW_dix = deta;
    let ddW_dixix = 0.0;

    np.vx = 0.5 * iy * ddW_dixix;
    np.vy = -ddW_dix - 0.5 * ix * ddW_dixix;
    np.vz = 0.5 * diz_dixix * dW + diz_dix * ddW_dix;

    np
}

pub fn reb_particle_derivative_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let xi = a * (clp + p / (2. - l) * h - k);
    let eta = a * (slp - p / (2. - l) * k - h);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let W = eta * ix - xi * iy;
    let dW_diy = -xi;

    np.x = 0.5 * W + 0.5 * iy * dW_diy;
    np.y = -0.5 * ix * dW_diy;
    np.z = 0.5 * diz_diy * W + 0.5 * iz * dW_diy;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let dxi = an / (1. - q) * (-slp + q / (2. - l) * h);
    let deta = an / (1. - q) * (clp - q / (2. - l) * k);
    let dW = deta * ix - dxi * iy;
    let ddW_diy = -dxi;

    np.vx = 0.5 * dW + 0.5 * iy * ddW_diy;
    np.vy = -0.5 * ix * ddW_diy;
    np.vz = 0.5 * diz_diy * dW + 0.5 * iz * ddW_diy;

    np
}

pub fn reb_particle_derivative_iy_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let xi = a * (clp + p / (2. - l) * h - k);
    let eta = a * (slp - p / (2. - l) * k - h);

    //double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diyiy = -1. / (4. - ix * ix - iy * iy).abs().sqrt()
        - iy * iy / ((4. - ix * ix - iy * iy).abs().sqrt() * (4. - ix * ix - iy * iy).abs().sqrt() * (4. - ix * ix - iy * iy).abs().sqrt());
    let W = eta * ix - xi * iy;
    let dW_diy = -xi;
    let dW_diyiy = 0.0;

    np.x = dW_diy + 0.5 * iy * dW_diyiy;
    np.y = -0.5 * ix * dW_diyiy;
    np.z = 0.5 * diz_diyiy * W + diz_diy * dW_diy;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let dxi = an / (1. - q) * (-slp + q / (2. - l) * h);
    let deta = an / (1. - q) * (clp - q / (2. - l) * k);
    let dW = deta * ix - dxi * iy;
    let ddW_diy = -dxi;
    let ddW_diyiy = 0.0;

    np.vx = ddW_diy - 0.5 * iy * ddW_diyiy;
    np.vy = -0.5 * ix * ddW_diyiy;
    np.vz = 0.5 * diz_diyiy * dW + diz_diy * ddW_diy;

    np
}

pub fn reb_particle_derivative_k_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dxi_dk = a * (dclp_dk + dp_dk / (2. - l) * h + p / ((2. - l) * (2. - l)) * dl_dk * h - 1.);
    let deta_dk = a * (dslp_dk - dp_dk / (2. - l) * k - p / (2. - l) - p / ((2. - l) * (2. - l)) * dl_dk * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dk = deta_dk * ix - dxi_dk * iy;
    let dW_dkix = deta_dk;

    np.x = 0.5 * iy * dW_dkix;
    np.y = -0.5 * dW_dk - 0.5 * ix * dW_dkix;
    np.z = 0.5 * diz_dix * dW_dk + 0.5 * iz * dW_dkix;

    let dq_dk = 1. / (1. - q) * (clp - k);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dk = dq_dk * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h);
    let ddeta_dk = dq_dk * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l));
    let ddW_dk = ddeta_dk * ix - ddxi_dk * iy;
    let ddW_dkix = ddeta_dk;

    np.vx = 0.5 * iy * ddW_dkix;
    np.vy = -0.5 * ddW_dk - 0.5 * ix * ddW_dkix;
    np.vz = 0.5 * diz_dix * ddW_dk + 0.5 * iz * ddW_dkix;

    np
}

pub fn reb_particle_derivative_h_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dxi_dh = a * (dclp_dh + dp_dh / (2. - l) * h + p / (2. - l) + p / ((2. - l) * (2. - l)) * dl_dh * h);
    let deta_dh = a * (dslp_dh - dp_dh / (2. - l) * k - p / ((2. - l) * (2. - l)) * k * dl_dh - 1.);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dh = deta_dh * ix - dxi_dh * iy;
    let dW_dhix = deta_dh;

    np.x = 0.5 * iy * dW_dhix;
    np.y = -0.5 * dW_dh - 0.5 * ix * dW_dhix;
    np.z = 0.5 * diz_dix * dW_dh + 0.5 * iz * dW_dhix;

    let dq_dh = 1. / (1. - q) * (slp - h);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dh = dq_dh * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l));
    let ddeta_dh = dq_dh * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k);
    let ddW_dh = ddeta_dh * ix - ddxi_dh * iy;
    let ddW_dhix = ddeta_dh;

    np.vx = 0.5 * iy * ddW_dhix;
    np.vy = -0.5 * ddW_dh - 0.5 * ix * ddW_dhix;
    np.vz = 0.5 * diz_dix * ddW_dh + 0.5 * iz * ddW_dhix;

    np
}

pub fn reb_particle_derivative_lambda_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_dlambda = a * (dclp_dlambda + dp_dlambda / (2. - l) * h);
    let deta_dlambda = a * (dslp_dlambda - dp_dlambda / (2. - l) * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dlambda = deta_dlambda * ix - dxi_dlambda * iy;
    let dW_dlambdaix = deta_dlambda;

    np.x = 0.5 * iy * dW_dlambdaix;
    np.y = -0.5 * dW_dlambda - 0.5 * ix * dW_dlambdaix;
    np.z = 0.5 * diz_dix * dW_dlambda + 0.5 * iz * dW_dlambdaix;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dlambda + dq_dlambda / (2. - l) * h);
    let ddeta_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dlambda - dq_dlambda / (2. - l) * k);
    let ddW_dlambda = ddeta_dlambda * ix - ddxi_dlambda * iy;
    let ddW_dlambdaix = ddeta_dlambda;
    np.vx = 0.5 * iy * ddW_dlambdaix;
    np.vy = -0.5 * ddW_dlambda - 0.5 * ix * ddW_dlambdaix;
    np.vz = 0.5 * diz_dix * ddW_dlambda + 0.5 * iz * ddW_dlambdaix;

    np
}

pub fn reb_particle_derivative_lambda_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_dlambda = a * (dclp_dlambda + dp_dlambda / (2. - l) * h);
    let deta_dlambda = a * (dslp_dlambda - dp_dlambda / (2. - l) * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dlambda = deta_dlambda * ix - dxi_dlambda * iy;
    let dW_dlambdaiy = -dxi_dlambda;
    np.x = 0.5 * dW_dlambda + 0.5 * iy * dW_dlambdaiy;
    np.y = -0.5 * ix * dW_dlambdaiy;
    np.z = 0.5 * diz_diy * dW_dlambda + 0.5 * iz * dW_dlambdaiy;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dlambda + dq_dlambda / (2. - l) * h);
    let ddeta_dlambda = an / ((1. - q) * (1. - q)) * dq_dlambda * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dlambda - dq_dlambda / (2. - l) * k);
    let ddW_dlambda = ddeta_dlambda * ix - ddxi_dlambda * iy;
    let ddW_dlambdaiy = -ddxi_dlambda;
    np.vx = 0.5 * ddW_dlambda + 0.5 * iy * ddW_dlambdaiy;
    np.vy = -0.5 * ix * ddW_dlambdaiy;
    np.vz = 0.5 * diz_diy * ddW_dlambda + 0.5 * iz * ddW_dlambdaiy;

    np
}

pub fn reb_particle_derivative_h_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dxi_dh = a * (dclp_dh + dp_dh / (2. - l) * h + p / (2. - l) + p / ((2. - l) * (2. - l)) * dl_dh * h);
    let deta_dh = a * (dslp_dh - dp_dh / (2. - l) * k - p / ((2. - l) * (2. - l)) * k * dl_dh - 1.);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dh = deta_dh * ix - dxi_dh * iy;
    let dW_dhiy = -dxi_dh;
    np.x = 0.5 * dW_dh + 0.5 * iy * dW_dhiy;
    np.y = -0.5 * ix * dW_dhiy;
    np.z = 0.5 * diz_diy * dW_dh + 0.5 * iz * dW_dhiy;

    let dq_dh = 1. / (1. - q) * (slp - h);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dh = dq_dh * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l));
    let ddeta_dh = dq_dh * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k);
    let ddW_dh = ddeta_dh * ix - ddxi_dh * iy;
    let ddW_dhiy = -ddxi_dh;
    np.vx = 0.5 * ddW_dh + 0.5 * iy * ddW_dhiy;
    np.vy = -0.5 * ix * ddW_dhiy;
    np.vz = 0.5 * diz_diy * ddW_dh + 0.5 * iz * ddW_dhiy;

    np
}

pub fn reb_particle_derivative_k_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dxi_dk = a * (dclp_dk + dp_dk / (2. - l) * h + p / ((2. - l) * (2. - l)) * dl_dk * h - 1.);
    let deta_dk = a * (dslp_dk - dp_dk / (2. - l) * k - p / (2. - l) - p / ((2. - l) * (2. - l)) * dl_dk * k);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dk = deta_dk * ix - dxi_dk * iy;
    let dW_dkiy = -dxi_dk;
    np.x = 0.5 * dW_dk + 0.5 * iy * dW_dkiy;
    np.y = -0.5 * ix * dW_dkiy;
    np.z = 0.5 * diz_diy * dW_dk + 0.5 * iz * dW_dkiy;

    let dq_dk = 1. / (1. - q) * (clp - k);

    let an = (G * (po.m + primary.m) / a).sqrt();
    let ddxi_dk = dq_dk * an / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + an / (1. - q) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h);
    let ddeta_dk = dq_dk * an / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + an / (1. - q) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l));
    let ddW_dk = ddeta_dk * ix - ddxi_dk * iy;
    let ddW_dkiy = -ddxi_dk;
    np.vx = 0.5 * ddW_dk + 0.5 * iy * ddW_dkiy;
    np.vy = -0.5 * ix * ddW_dkiy;
    np.vz = 0.5 * diz_diy * ddW_dk + 0.5 * iz * ddW_dkiy;

    np
}

pub fn reb_particle_derivative_ix_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let xi = a * (clp + p / (2. - l) * h - k);
    let eta = a * (slp - p / (2. - l) * k - h);

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dixiy = -ix * iy / ((4. - ix * ix - iy * iy).abs().sqrt() * (4. - ix * ix - iy * iy).abs().sqrt() * (4. - ix * ix - iy * iy).abs().sqrt());
    let W = eta * ix - xi * iy;
    let dW_dix = eta;
    let dW_diy = -xi;
    let dW_dixiy = 0.0;
    np.x = 0.5 * dW_dix + 0.5 * iy * dW_dixiy;
    np.y = -0.5 * dW_diy - 0.5 * ix * dW_dixiy;
    np.z = 0.5 * diz_dixiy * W + 0.5 * diz_dix * dW_diy + 0.5 * diz_diy * dW_dix + 0.5 * iz * dW_dixiy;

    let an = (G * (po.m + primary.m) / a).sqrt();
    let dxi = an / (1. - q) * (-slp + q / (2. - l) * h);
    let deta = an / (1. - q) * (clp - q / (2. - l) * k);
    let dW = deta * ix - dxi * iy;
    let ddW_dix = deta;
    let ddW_diy = -dxi;
    let ddW_dixiy = 0.0;

    np.vx = 0.5 * ddW_dix + 0.5 * iy * ddW_dixiy;
    np.vy = -0.5 * ddW_diy - 0.5 * ix * ddW_dixiy;
    np.vz = 0.5 * diz_dixiy * dW + 0.5 * diz_dix * ddW_diy + 0.5 * diz_diy * ddW_dix + 0.5 * iz * ddW_dixiy;

    np
}

pub fn reb_particle_derivative_a_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let deta_da = slp - p / (2. - l) * k - h;
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_daix = deta_da;
    let dxi_da = clp + p / (2. - l) * h - k;
    let dW_da = deta_da * ix - dxi_da * iy;

    np.x = 0.5 * iy * dW_daix;
    np.y = -0.5 * dW_da - 0.5 * ix * dW_daix;
    np.z = 0.5 * diz_dix * dW_da + 0.5 * iz * dW_daix;

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddeta_da = dan_da / (1. - q) * (clp - q / (2. - l) * k);
    let ddW_daix = ddeta_da;
    let ddxi_da = dan_da / (1. - q) * (-slp + q / (2. - l) * h);
    let ddW_da = ddeta_da * ix - ddxi_da * iy;

    np.vx = 0.5 * iy * ddW_daix;
    np.vy = -0.5 * ddW_da - 0.5 * ix * ddW_daix;
    np.vz = 0.5 * diz_dix * ddW_da + 0.5 * iz * ddW_daix;

    np
}

pub fn reb_particle_derivative_a_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();
    let dxi_da = clp + p / (2. - l) * h - k;
    let deta_da = slp - p / (2. - l) * k - h;
    let dW_da = deta_da * ix - dxi_da * iy;
    let dW_daiy = -dxi_da;
    np.x = 0.5 * dW_da + 0.5 * iy * dW_daiy;
    np.y = -0.5 * ix * dW_daiy;
    np.z = 0.5 * diz_diy * dW_da + 0.5 * iz * dW_daiy;

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddxi_da = dan_da / (1. - q) * (-slp + q / (2. - l) * h);
    let ddW_daiy = -ddxi_da;
    let ddeta_da = dan_da / (1. - q) * (clp - q / (2. - l) * k);
    let ddW_da = ddeta_da * ix - ddxi_da * iy;

    np.vx = 0.5 * ddW_da + 0.5 * iy * ddW_daiy;
    np.vy = -0.5 * ix * ddW_daiy;
    np.vz = 0.5 * diz_diy * ddW_da + 0.5 * iz * ddW_daiy;

    np
}

pub fn reb_particle_derivative_a_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();

    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);
    let dp_dlambda = q / (1. - q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dxi_dalambda = dclp_dlambda + dp_dlambda / (2. - l) * h;
    let deta_dalambda = dslp_dlambda - dp_dlambda / (2. - l) * k;

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dalambda = deta_dalambda * ix - dxi_dalambda * iy;

    np.x = dxi_dalambda + 0.5 * iy * dW_dalambda;
    np.y = deta_dalambda - 0.5 * ix * dW_dalambda;
    np.z = 0.5 * iz * dW_dalambda;

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddxi_dalambda = dan_da / ((1. - q) * (1. - q)) * dq_dlambda * (-slp + q / (2. - l) * h)
        + dan_da / (1. - q) * (-dslp_dlambda + dq_dlambda / (2. - l) * h);
    let ddeta_dalambda = dan_da / ((1. - q) * (1. - q)) * dq_dlambda * (clp - q / (2. - l) * k)
        + dan_da / (1. - q) * (dclp_dlambda - dq_dlambda / (2. - l) * k);
    let ddW_dalambda = ddeta_dalambda * ix - ddxi_dalambda * iy;
    np.vx = ddxi_dalambda + 0.5 * iy * ddW_dalambda;
    np.vy = ddeta_dalambda - 0.5 * ix * ddW_dalambda;
    np.vz = 0.5 * iz * ddW_dalambda;

    np
}

pub fn reb_particle_derivative_a_h(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let dp_dh = 1. / (1. - q) * (-clp);
    let dxi_dah = dclp_dh + dp_dh / (2. - l) * h + p / (2. - l) + p / ((2. - l) * (2. - l)) * dl_dh * h;
    let deta_dah = dslp_dh - dp_dh / (2. - l) * k - p / ((2. - l) * (2. - l)) * k * dl_dh - 1.;

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dah = deta_dah * ix - dxi_dah * iy;

    np.x = dxi_dah + 0.5 * iy * dW_dah;
    np.y = deta_dah - 0.5 * ix * dW_dah;
    np.z = 0.5 * iz * dW_dah;

    let dq_dh = 1. / (1. - q) * (slp - h);

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddxi_dah = dq_dh * dan_da / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dan_da / (1. - q) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l));
    let ddeta_dah = dq_dh * dan_da / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dan_da / (1. - q) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k);
    let ddW_dah = ddeta_dah * ix - ddxi_dah * iy;

    np.vx = ddxi_dah + 0.5 * iy * ddW_dah;
    np.vy = ddeta_dah - 0.5 * ix * ddW_dah;
    np.vz = 0.5 * iz * ddW_dah;

    np
}

pub fn reb_particle_derivative_a_k(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let dp_dk = 1. / (1. - q) * (slp);
    let dxi_dak = dclp_dk + dp_dk / (2. - l) * h + p / ((2. - l) * (2. - l)) * dl_dk * h - 1.;
    let deta_dak = dslp_dk - dp_dk / (2. - l) * k - p / (2. - l) - p / ((2. - l) * (2. - l)) * dl_dk * k;

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let dW_dak = deta_dak * ix - dxi_dak * iy;

    np.x = dxi_dak + 0.5 * iy * dW_dak;
    np.y = deta_dak - 0.5 * ix * dW_dak;
    np.z = 0.5 * iz * dW_dak;

    let dq_dk = 1. / (1. - q) * (clp - k);

    let dan_da = -0.5 * (G * (po.m + primary.m) / (a * a * a)).sqrt();
    let ddxi_dak = dq_dk * dan_da / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dan_da / (1. - q) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h);
    let ddeta_dak = dq_dk * dan_da / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dan_da / (1. - q) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l));
    let ddW_dak = ddeta_dak * ix - ddxi_dak * iy;

    np.vx = ddxi_dak + 0.5 * iy * ddW_dak;
    np.vy = ddeta_dak - 0.5 * ix * ddW_dak;
    np.vz = 0.5 * iz * ddW_dak;

    np
}

pub fn reb_particle_derivative_m(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    np.m = 1.;
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dm = dan_dm / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_dm = dan_dm / (1. - q) * (clp - q / (2. - l) * k);

    let ddW_dm = ddeta_dm * ix - ddxi_dm * iy;
    np.vx = ddxi_dm + 0.5 * iy * ddW_dm;
    np.vy = ddeta_dm - 0.5 * ix * ddW_dm;
    np.vz = 0.5 * iz * ddW_dm;

    np
}

pub fn reb_particle_derivative_m_a(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let l = 1. - (1. - h * h - k * k).sqrt();

    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dma = -0.5 * 0.5 * (G / (a * a * a * (po.m + primary.m))).sqrt();
    let ddxi_dma = dan_dma / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_dma = dan_dma / (1. - q) * (clp - q / (2. - l) * k);

    let ddW_dma = ddeta_dma * ix - ddxi_dma * iy;
    np.vx = ddxi_dma + 0.5 * iy * ddW_dma;
    np.vy = ddeta_dma - 0.5 * ix * ddW_dma;
    np.vz = 0.5 * iz * ddW_dma;

    np
}

pub fn reb_particle_derivative_m_lambda(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let dq_dlambda = -p / (1. - q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dlambda = -1. / (1. - q) * slp;
    let dslp_dlambda = 1. / (1. - q) * clp;

    let l = 1. - (1. - h * h - k * k).sqrt();
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dmlambda = dan_dm / ((1. - q) * (1. - q)) * dq_dlambda * (-slp + q / (2. - l) * h)
        + dan_dm / (1. - q) * (-dslp_dlambda + dq_dlambda / (2. - l) * h);
    let ddeta_dmlambda = dan_dm / ((1. - q) * (1. - q)) * dq_dlambda * (clp - q / (2. - l) * k)
        + dan_dm / (1. - q) * (dclp_dlambda - dq_dlambda / (2. - l) * k);
    let ddW_dmlambda = ddeta_dmlambda * ix - ddxi_dmlambda * iy;
    np.vx = ddxi_dmlambda + 0.5 * iy * ddW_dmlambda;
    np.vy = ddeta_dmlambda - 0.5 * ix * ddW_dmlambda;
    np.vz = 0.5 * iz * ddW_dmlambda;

    np
}

pub fn reb_particle_derivative_m_h(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dh = -1. / (1. - q) * (-slp * clp);
    let dslp_dh = -1. / (1. - q) * (clp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dh = 1. / (1. - h * h - k * k).sqrt() * h;
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dq_dh = 1. / (1. - q) * (slp - h);

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dmh = dq_dh * dan_dm / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dan_dm / (1. - q) * (-dslp_dh + dq_dh / (2. - l) * h + dl_dh * q / ((2. - l) * (2. - l)) * h + q / (2. - l));
    let ddeta_dmh = dq_dh * dan_dm / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dan_dm / (1. - q) * (dclp_dh - dq_dh / (2. - l) * k - dl_dh * q / ((2. - l) * (2. - l)) * k);
    let ddW_dmh = ddeta_dmh * ix - ddxi_dmh * iy;

    np.vx = ddxi_dmh + 0.5 * iy * ddW_dmh;
    np.vy = ddeta_dmh - 0.5 * ix * ddW_dmh;
    np.vz = 0.5 * iz * ddW_dmh;

    np
}

pub fn reb_particle_derivative_m_k(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();
    let dclp_dk = -1. / (1. - q) * (slp * slp);
    let dslp_dk = -1. / (1. - q) * (-slp * clp);

    let l = 1. - (1. - h * h - k * k).sqrt();
    let dl_dk = 1. / (1. - h * h - k * k).sqrt() * k;
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dq_dk = 1. / (1. - q) * (clp - k);

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dmk = dq_dk * dan_dm / ((1. - q) * (1. - q)) * (-slp + q / (2. - l) * h)
        + dan_dm / (1. - q) * (-dslp_dk + dq_dk / (2. - l) * h + dl_dk * q / ((2. - l) * (2. - l)) * h);
    let ddeta_dmk = dq_dk * dan_dm / ((1. - q) * (1. - q)) * (clp - q / (2. - l) * k)
        + dan_dm / (1. - q) * (dclp_dk - dq_dk / (2. - l) * k - dl_dk * q / ((2. - l) * (2. - l)) * k - q / (2. - l));
    let ddW_dmk = ddeta_dmk * ix - ddxi_dmk * iy;

    np.vx = ddxi_dmk + 0.5 * iy * ddW_dmk;
    np.vy = ddeta_dmk - 0.5 * ix * ddW_dmk;
    np.vz = 0.5 * iz * ddW_dmk;

    np
}

pub fn reb_particle_derivative_m_ix(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_dix = -ix / (4. - ix * ix - iy * iy).abs().sqrt();

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dm = dan_dm / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_dm = dan_dm / (1. - q) * (clp - q / (2. - l) * k);
    let ddW_dm = ddeta_dm * ix - ddxi_dm * iy;
    let ddW_dmix = ddeta_dm;

    np.vx = 0.5 * iy * ddW_dmix;
    np.vy = -0.5 * ddW_dm - 0.5 * ix * ddW_dmix;
    np.vz = 0.5 * diz_dix * ddW_dm + 0.5 * iz * ddW_dmix;

    np
}

pub fn reb_particle_derivative_m_iy(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);
    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    let diz_diy = -iy / (4. - ix * ix - iy * iy).abs().sqrt();

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dm = 0.5 * (G / (a * (po.m + primary.m))).sqrt();
    let ddxi_dm = dan_dm / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_dm = dan_dm / (1. - q) * (clp - q / (2. - l) * k);
    let ddW_dm = ddeta_dm * ix - ddxi_dm * iy;
    let ddW_dmiy = -ddxi_dm;

    np.vx = 0.5 * ddW_dm + 0.5 * iy * ddW_dmiy;
    np.vy = -0.5 * ix * ddW_dmiy;
    np.vz = 0.5 * diz_diy * ddW_dm + 0.5 * iz * ddW_dmiy;

    np
}

pub fn reb_particle_derivative_m_m(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let (mut a, mut lambda, mut k, mut h, mut ix, mut iy) = (0., 0., 0., 0., 0., 0.);
    reb_tools_particle_to_pal(G, po, primary, &mut a, &mut lambda, &mut k, &mut h, &mut ix, &mut iy);

    let mut np = reb_particle::default();
    let mut p = 0.;
    let mut q = 0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &mut p, &mut q);

    let slp = (lambda + p).sin();
    let clp = (lambda + p).cos();

    let l = 1. - (1. - h * h - k * k).sqrt();
    let iz = (4. - ix * ix - iy * iy).abs().sqrt();
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    let dan_dmm = -0.25 * (G / (a * (po.m + primary.m) * (po.m + primary.m) * (po.m + primary.m))).sqrt();
    let ddxi_dmm = dan_dmm / (1. - q) * (-slp + q / (2. - l) * h);
    let ddeta_dmm = dan_dmm / (1. - q) * (clp - q / (2. - l) * k);

    let ddW_dmm = ddeta_dmm * ix - ddxi_dmm * iy;
    np.vx = ddxi_dmm + 0.5 * iy * ddW_dmm;
    np.vy = ddeta_dmm - 0.5 * ix * ddW_dmm;
    np.vz = 0.5 * iz * ddW_dmm;

    np
}

pub fn reb_particle_derivative_e(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dr = -o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0 = (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y = dr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z = dr * (so * cf + co * sf) * si;

    p.vx = dv0 * ((o.e + cf) * (-ci * co * sO - cO * so) - sf * (co * cO - ci * so * sO));
    p.vy = dv0 * ((o.e + cf) * (ci * co * cO - sO * so) - sf * (co * sO + ci * so * cO));
    p.vz = dv0 * ((o.e + cf) * co * si - sf * si * so);

    p.vx += v0 * (-ci * co * sO - cO * so);
    p.vy += v0 * (ci * co * cO - sO * so);
    p.vz += v0 * (co * si);

    p
}

pub fn reb_particle_derivative_e_e(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let ddr = o.a * 2. * (cosf * cosf - 1.) / ((cosf * o.e + 1.) * (cosf * o.e + 1.) * (cosf * o.e + 1.));
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dv0 = o.e * v0 / (1. - o.e * o.e);
    let ddv0 = v0 / ((o.e * o.e - 1.) * (o.e * o.e - 1.)) * (2. * o.e * o.e + 1.);

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = ddr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y = ddr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z = ddr * (so * cf + co * sf) * si;

    p.vx = ddv0 * ((o.e + cf) * (-ci * co * sO - cO * so) - sf * (co * cO - ci * so * sO));
    p.vy = ddv0 * ((o.e + cf) * (ci * co * cO - sO * so) - sf * (co * sO + ci * so * cO));
    p.vz = ddv0 * ((o.e + cf) * co * si - sf * si * so);

    p.vx += 2. * dv0 * (-ci * co * sO - cO * so);
    p.vy += 2. * dv0 * (ci * co * cO - sO * so);
    p.vz += 2. * dv0 * (co * si);

    p
}

pub fn reb_particle_derivative_inc(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.x = r * (-sO * (so * cf + co * sf) * dci);
    p.y = r * (cO * (so * cf + co * sf) * dci);
    p.z = r * (so * cf + co * sf) * dsi;

    p.vx = v0 * ((o.e + cf) * (-dci * co * sO) - sf * (-dci * so * sO));
    p.vy = v0 * ((o.e + cf) * (dci * co * cO) - sf * (dci * so * cO));
    p.vz = v0 * ((o.e + cf) * co * dsi - sf * dsi * so);

    p
}

pub fn reb_particle_derivative_inc_inc(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ddci = -o.inc.cos();
    let ddsi = -o.inc.sin();

    p.x = r * (-sO * (so * cf + co * sf) * ddci);
    p.y = r * (cO * (so * cf + co * sf) * ddci);
    p.z = r * (so * cf + co * sf) * ddsi;

    p.vx = v0 * ((o.e + cf) * (-ddci * co * sO) - sf * (-ddci * so * sO));
    p.vy = v0 * ((o.e + cf) * (ddci * co * cO) - sf * (ddci * so * cO));
    p.vz = v0 * ((o.e + cf) * co * ddsi - sf * ddsi * so);

    p
}

pub fn reb_particle_derivative_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = r * (dcO * (co * cf - so * sf) - dsO * (so * cf + co * sf) * ci);
    p.y = r * (dsO * (co * cf - so * sf) + dcO * (so * cf + co * sf) * ci);
    p.z = 0.;

    p.vx = v0 * ((o.e + cf) * (-ci * co * dsO - dcO * so) - sf * (co * dcO - ci * so * dsO));
    p.vy = v0 * ((o.e + cf) * (ci * co * dcO - dsO * so) - sf * (co * dsO + ci * so * dcO));
    p.vz = 0.;

    p
}

pub fn reb_particle_derivative_Omega_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let ddcO = -o.Omega.cos();
    let ddsO = -o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = r * (ddcO * (co * cf - so * sf) - ddsO * (so * cf + co * sf) * ci);
    p.y = r * (ddsO * (co * cf - so * sf) + ddcO * (so * cf + co * sf) * ci);
    p.z = 0.;

    p.vx = v0 * ((o.e + cf) * (-ci * co * ddsO - ddcO * so) - sf * (co * ddcO - ci * so * ddsO));
    p.vy = v0 * ((o.e + cf) * (ci * co * ddcO - ddsO * so) - sf * (co * ddsO + ci * so * ddcO));
    p.vz = 0.;

    p
}

pub fn reb_particle_derivative_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = r * (cO * (dco * cf - dso * sf) - sO * (dso * cf + dco * sf) * ci);
    p.y = r * (sO * (dco * cf - dso * sf) + cO * (dso * cf + dco * sf) * ci);
    p.z = r * (dso * cf + dco * sf) * si;

    p.vx = v0 * ((o.e + cf) * (-ci * dco * sO - cO * dso) - sf * (dco * cO - ci * dso * sO));
    p.vy = v0 * ((o.e + cf) * (ci * dco * cO - sO * dso) - sf * (dco * sO + ci * dso * cO));
    p.vz = v0 * ((o.e + cf) * dco * si - sf * si * dso);

    p
}

pub fn reb_particle_derivative_omega_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let ddco = -o.omega.cos();
    let ddso = -o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = r * (cO * (ddco * cf - ddso * sf) - sO * (ddso * cf + ddco * sf) * ci);
    p.y = r * (sO * (ddco * cf - ddso * sf) + cO * (ddso * cf + ddco * sf) * ci);
    p.z = r * (ddso * cf + ddco * sf) * si;

    p.vx = v0 * ((o.e + cf) * (-ci * ddco * sO - cO * ddso) - sf * (ddco * cO - ci * ddso * sO));
    p.vy = v0 * ((o.e + cf) * (ci * ddco * cO - sO * ddso) - sf * (ddco * sO + ci * ddso * cO));
    p.vz = v0 * ((o.e + cf) * ddco * si - sf * si * ddso);

    p
}

pub fn reb_particle_derivative_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dr = o.a * (1. - o.e * o.e) / ((1. + o.e * o.f.cos()) * (1. + o.e * o.f.cos())) * o.e * o.f.sin();
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y = dr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z = dr * (so * cf + co * sf) * si;

    p.x += r * (cO * (co * dcf - so * dsf) - sO * (so * dcf + co * dsf) * ci);
    p.y += r * (sO * (co * dcf - so * dsf) + cO * (so * dcf + co * dsf) * ci);
    p.z += r * (so * dcf + co * dsf) * si;

    p.vx = v0 * (dcf * (-ci * co * sO - cO * so) - dsf * (co * cO - ci * so * sO));
    p.vy = v0 * (dcf * (ci * co * cO - sO * so) - dsf * (co * sO + ci * so * cO));
    p.vz = v0 * (dcf * co * si - dsf * si * so);

    p
}

pub fn reb_particle_derivative_f_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dr = o.a * (1. - o.e * o.e) / ((1. + o.e * o.f.cos()) * (1. + o.e * o.f.cos())) * o.e * o.f.sin();
    let ddr = 2. * o.a * (1. - o.e * o.e) / ((1. + o.e * o.f.cos()) * (1. + o.e * o.f.cos()) * (1. + o.e * o.f.cos())) * o.e * o.e * o.f.sin() * o.f.sin() + o.a * (1. - o.e * o.e) * o.e * o.f.cos() / ((1. + o.e * o.f.cos()) * (1. + o.e * o.f.cos()));
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let ddcf = -o.f.cos();
    let ddsf = -o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = ddr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y = ddr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z = ddr * (so * cf + co * sf) * si;

    p.x += 2. * dr * (cO * (co * dcf - so * dsf) - sO * (so * dcf + co * dsf) * ci);
    p.y += 2. * dr * (sO * (co * dcf - so * dsf) + cO * (so * dcf + co * dsf) * ci);
    p.z += 2. * dr * (so * dcf + co * dsf) * si;

    p.x += r * (cO * (co * ddcf - so * ddsf) - sO * (so * ddcf + co * ddsf) * ci);
    p.y += r * (sO * (co * ddcf - so * ddsf) + cO * (so * ddcf + co * ddsf) * ci);
    p.z += r * (so * ddcf + co * ddsf) * si;

    p.vx = v0 * (ddcf * (-ci * co * sO - cO * so) - ddsf * (co * cO - ci * so * sO));
    p.vy = v0 * (ddcf * (ci * co * cO - sO * so) - ddsf * (co * sO + ci * so * cO));
    p.vz = v0 * (ddcf * co * si - ddsf * si * so);

    p
}

pub fn reb_particle_derivative_a_e(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let ddr = -(cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0_da = -0.5 / (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt() * G * (po.m + primary.m) / (o.a * o.a) / (1. - o.e * o.e);

    let dv0_da_de = o.e * dv0_da / (1. - o.e * o.e);

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = ddr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y = ddr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z = ddr * (so * cf + co * sf) * si;

    p.vx = dv0_da_de * ((o.e + cf) * (-ci * co * sO - cO * so) - sf * (co * cO - ci * so * sO));
    p.vy = dv0_da_de * ((o.e + cf) * (ci * co * cO - sO * so) - sf * (co * sO + ci * so * cO));
    p.vz = dv0_da_de * ((o.e + cf) * co * si - sf * si * so);

    p.vx += dv0_da * (-ci * co * sO - cO * so);
    p.vy += dv0_da * (ci * co * cO - sO * so);
    p.vz += dv0_da * (co * si);

    p
}

pub fn reb_particle_derivative_a_inc(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dr = (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dv0 = -0.5 / (o.a * o.a * o.a).sqrt() * (G * (po.m + primary.m) / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.x = dr * (-sO * (so * cf + co * sf) * dci);
    p.y = dr * (cO * (so * cf + co * sf) * dci);
    p.z = dr * (so * cf + co * sf) * dsi;

    p.vx = dv0 * ((o.e + cf) * (-dci * co * sO) - sf * (-dci * so * sO));
    p.vy = dv0 * ((o.e + cf) * (dci * co * cO) - sf * (dci * so * cO));
    p.vz = dv0 * ((o.e + cf) * co * dsi - sf * dsi * so);

    p
}

pub fn reb_particle_derivative_a_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dr = (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dv0 = -0.5 / (o.a * o.a * o.a).sqrt() * (G * (po.m + primary.m) / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = dr * (dcO * (co * cf - so * sf) - dsO * (so * cf + co * sf) * ci);
    p.y = dr * (dsO * (co * cf - so * sf) + dcO * (so * cf + co * sf) * ci);

    p.vx = dv0 * ((o.e + cf) * (-ci * co * dsO - dcO * so) - sf * (co * dcO - ci * so * dsO));
    p.vy = dv0 * ((o.e + cf) * (ci * co * dcO - dsO * so) - sf * (co * dsO + ci * so * dcO));

    p
}

pub fn reb_particle_derivative_a_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dr = (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dv0 = -0.5 / (o.a * o.a * o.a).sqrt() * (G * (po.m + primary.m) / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (dco * cf - dso * sf) - sO * (dso * cf + dco * sf) * ci);
    p.y = dr * (sO * (dco * cf - dso * sf) + cO * (dso * cf + dco * sf) * ci);
    p.z = dr * (dso * cf + dco * sf) * si;

    p.vx = dv0 * ((o.e + cf) * (-ci * dco * sO - cO * dso) - sf * (dco * cO - ci * dso * sO));
    p.vy = dv0 * ((o.e + cf) * (ci * dco * cO - sO * dso) - sf * (dco * sO + ci * dso * cO));
    p.vz = dv0 * ((o.e + cf) * dco * si - sf * si * dso);

    p
}

pub fn reb_particle_derivative_a_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dr = (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let ddr = o.e * o.f.sin() * (1. - o.e * o.e) / (1. + o.e * o.f.cos()) / (1. + o.e * o.f.cos());
    let dv0 = -0.5 / (o.a * o.a * o.a).sqrt() * (G * (po.m + primary.m) / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (co * dcf - so * dsf) - sO * (so * dcf + co * dsf) * ci);
    p.y = dr * (sO * (co * dcf - so * dsf) + cO * (so * dcf + co * dsf) * ci);
    p.z = dr * (so * dcf + co * dsf) * si;

    p.x += ddr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y += ddr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z += ddr * (so * cf + co * sf) * si;

    p.vx = dv0 * (dcf * (-ci * co * sO - cO * so) - dsf * (co * cO - ci * so * sO));
    p.vy = dv0 * (dcf * (ci * co * cO - sO * so) - dsf * (co * sO + ci * so * cO));
    p.vz = dv0 * (dcf * co * si - dsf * si * so);

    p
}

pub fn reb_particle_derivative_e_inc(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dr = -o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0 = (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.x = dr * (-sO * (so * cf + co * sf) * dci);
    p.y = dr * (cO * (so * cf + co * sf) * dci);
    p.z = dr * (so * cf + co * sf) * dsi;

    p.vx = dv0 * ((o.e + cf) * (-dci * co * sO) - sf * (-dci * so * sO));
    p.vy = dv0 * ((o.e + cf) * (dci * co * cO) - sf * (dci * so * cO));
    p.vz = dv0 * ((o.e + cf) * co * dsi - sf * dsi * so);

    p.vx += v0 * (-dci * co * sO);
    p.vy += v0 * (dci * co * cO);
    p.vz += v0 * (co * dsi);

    p
}

pub fn reb_particle_derivative_e_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dr = -o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0 = (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = dr * (dcO * (co * cf - so * sf) - dsO * (so * cf + co * sf) * ci);
    p.y = dr * (dsO * (co * cf - so * sf) + dcO * (so * cf + co * sf) * ci);

    p.vx = dv0 * ((o.e + cf) * (-ci * co * dsO - dcO * so) - sf * (co * dcO - ci * so * dsO));
    p.vy = dv0 * ((o.e + cf) * (ci * co * dcO - dsO * so) - sf * (co * dsO + ci * so * dcO));

    p.vx += v0 * (-ci * co * dsO - dcO * so);
    p.vy += v0 * (ci * co * dcO - dsO * so);

    p
}

pub fn reb_particle_derivative_e_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dr = -o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0 = (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (dco * cf - dso * sf) - sO * (dso * cf + dco * sf) * ci);
    p.y = dr * (sO * (dco * cf - dso * sf) + cO * (dso * cf + dco * sf) * ci);
    p.z = dr * (dso * cf + dco * sf) * si;

    p.vx = dv0 * ((o.e + cf) * (-ci * dco * sO - cO * dso) - sf * (dco * cO - ci * dso * sO));
    p.vy = dv0 * ((o.e + cf) * (ci * dco * cO - sO * dso) - sf * (dco * sO + ci * dso * cO));
    p.vz = dv0 * ((o.e + cf) * dco * si - sf * si * dso);

    p.vx += v0 * (-ci * dco * sO - cO * dso);
    p.vy += v0 * (ci * dco * cO - sO * dso);
    p.vz += v0 * (dco * si);

    p
}

pub fn reb_particle_derivative_e_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let cosf = o.f.cos();
    let dr = -o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.));
    let ddr = -o.a * (-o.f.sin() * o.e * o.e - o.f.sin()) / ((cosf * o.e + 1.) * (cosf * o.e + 1.))
        - 2. * o.e * o.f.sin() * o.a * (cosf * o.e * o.e + cosf + 2. * o.e) / ((cosf * o.e + 1.) * (cosf * o.e + 1.) * (cosf * o.e + 1.));
    let dv0 = (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = dr * (cO * (co * dcf - so * dsf) - sO * (so * dcf + co * dsf) * ci);
    p.y = dr * (sO * (co * dcf - so * dsf) + cO * (so * dcf + co * dsf) * ci);
    p.z = dr * (so * dcf + co * dsf) * si;

    p.x += ddr * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci);
    p.y += ddr * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci);
    p.z += ddr * (so * cf + co * sf) * si;

    p.vx = dv0 * (dcf * (-ci * co * sO - cO * so) - dsf * (co * cO - ci * so * sO));
    p.vy = dv0 * (dcf * (ci * co * cO - sO * so) - dsf * (co * sO + ci * so * cO));
    p.vz = dv0 * (dcf * co * si - dsf * si * so);

    p
}

pub fn reb_particle_derivative_m_e(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dv0m = 0.5 * G / o.a / (1. - o.e * o.e) / (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();
    let dv0ea = 0.5 * G / o.a / (G * (po.m + primary.m) / o.a).sqrt() * o.e / ((1. - o.e * o.e) * (1. - o.e * o.e).sqrt());

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.vx = dv0ea * ((o.e + cf) * (-ci * co * sO - cO * so) - sf * (co * cO - ci * so * sO));
    p.vy = dv0ea * ((o.e + cf) * (ci * co * cO - sO * so) - sf * (co * sO + ci * so * cO));
    p.vz = dv0ea * ((o.e + cf) * co * si - sf * si * so);

    p.vx += dv0m * (-ci * co * sO - cO * so);
    p.vy += dv0m * (ci * co * cO - sO * so);
    p.vz += dv0m * (co * si);

    p
}

pub fn reb_particle_derivative_inc_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();

    p.x = r * (-dsO * (so * cf + co * sf) * dci);
    p.y = r * (dcO * (so * cf + co * sf) * dci);

    p.vx = v0 * ((o.e + cf) * (-dci * co * dsO) - sf * (-dci * so * dsO));
    p.vy = v0 * ((o.e + cf) * (dci * co * dcO) - sf * (dci * so * dcO));

    p
}

pub fn reb_particle_derivative_inc_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.x = r * (-sO * (dso * cf + dco * sf) * dci);
    p.y = r * (cO * (dso * cf + dco * sf) * dci);
    p.z = r * (dso * cf + dco * sf) * dsi;

    p.vx = v0 * ((o.e + cf) * (-dci * dco * sO) - sf * (-dci * dso * sO));
    p.vy = v0 * ((o.e + cf) * (dci * dco * cO) - sf * (dci * dso * cO));
    p.vz = v0 * ((o.e + cf) * dco * dsi - sf * dsi * dso);

    p
}

pub fn reb_particle_derivative_inc_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dr = o.e * o.f.sin() * o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos()) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.x = r * (-sO * (so * dcf + co * dsf) * dci);
    p.y = r * (cO * (so * dcf + co * dsf) * dci);
    p.z = r * (so * dcf + co * dsf) * dsi;

    p.x += dr * (-sO * (so * cf + co * sf) * dci);
    p.y += dr * (cO * (so * cf + co * sf) * dci);
    p.z += dr * (so * cf + co * sf) * dsi;

    p.vx = v0 * (dcf * (-dci * co * sO) - dsf * (-dci * so * sO));
    p.vy = v0 * (dcf * (dci * co * cO) - dsf * (dci * so * cO));
    p.vz = v0 * (dcf * co * dsi - dsf * dsi * so);

    p
}

pub fn reb_particle_derivative_m_inc(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dv0 = 0.5 / (po.m + primary.m).sqrt() * (G / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let dci = -o.inc.sin();
    let dsi = o.inc.cos();

    p.vx = dv0 * ((o.e + cf) * (-dci * co * sO) - sf * (-dci * so * sO));
    p.vy = dv0 * ((o.e + cf) * (dci * co * cO) - sf * (dci * so * cO));
    p.vz = dv0 * ((o.e + cf) * co * dsi - sf * dsi * so);

    p
}

pub fn reb_particle_derivative_omega_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = r * (dcO * (dco * cf - dso * sf) - dsO * (dso * cf + dco * sf) * ci);
    p.y = r * (dsO * (dco * cf - dso * sf) + dcO * (dso * cf + dco * sf) * ci);

    p.vx = v0 * ((o.e + cf) * (-ci * dco * dsO - dcO * dso) - sf * (dco * dcO - ci * dso * dsO));
    p.vy = v0 * ((o.e + cf) * (ci * dco * dcO - dsO * dso) - sf * (dco * dsO + ci * dso * dcO));

    p
}

pub fn reb_particle_derivative_Omega_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dr = o.e * o.f.sin() * o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos()) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.x = r * (dcO * (co * dcf - so * dsf) - dsO * (so * dcf + co * dsf) * ci);
    p.y = r * (dsO * (co * dcf - so * dsf) + dcO * (so * dcf + co * dsf) * ci);

    p.x += dr * (dcO * (co * cf - so * sf) - dsO * (so * cf + co * sf) * ci);
    p.y += dr * (dsO * (co * cf - so * sf) + dcO * (so * cf + co * sf) * ci);

    p.vx = v0 * ((dcf) * (-ci * co * dsO - dcO * so) - dsf * (co * dcO - ci * so * dsO));
    p.vy = v0 * ((dcf) * (ci * co * dcO - dsO * so) - dsf * (co * dsO + ci * so * dcO));

    p
}

pub fn reb_particle_derivative_m_Omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dv0 = 0.5 / (po.m + primary.m).sqrt() * (G / o.a / (1. - o.e * o.e)).sqrt();

    let dcO = -o.Omega.sin();
    let dsO = o.Omega.cos();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();

    p.vx = dv0 * ((o.e + cf) * (-ci * co * dsO - dcO * so) - sf * (co * dcO - ci * so * dsO));
    p.vy = dv0 * ((o.e + cf) * (ci * co * dcO - dsO * so) - sf * (co * dsO + ci * so * dcO));

    p
}

pub fn reb_particle_derivative_omega_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let r = o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos());
    let dr = o.e * o.f.sin() * o.a * (1. - o.e * o.e) / (1. + o.e * o.f.cos()) / (1. + o.e * o.f.cos());
    let v0 = (G * (po.m + primary.m) / o.a / (1. - o.e * o.e)).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.x = r * (cO * (dco * dcf - dso * dsf) - sO * (dso * dcf + dco * dsf) * ci);
    p.y = r * (sO * (dco * dcf - dso * dsf) + cO * (dso * dcf + dco * dsf) * ci);
    p.z = r * (dso * dcf + dco * dsf) * si;

    p.x += dr * (cO * (dco * cf - dso * sf) - sO * (dso * cf + dco * sf) * ci);
    p.y += dr * (sO * (dco * cf - dso * sf) + cO * (dso * cf + dco * sf) * ci);
    p.z += dr * (dso * cf + dco * sf) * si;

    p.vx = v0 * ((dcf) * (-ci * dco * sO - cO * dso) - dsf * (dco * cO - ci * dso * sO));
    p.vy = v0 * ((dcf) * (ci * dco * cO - sO * dso) - dsf * (dco * sO + ci * dso * cO));
    p.vz = v0 * ((dcf) * dco * si - dsf * si * dso);

    p
}

pub fn reb_particle_derivative_m_omega(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dv0 = 0.5 * (G / o.a / (1. - o.e * o.e)).sqrt() / (po.m + primary.m).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let dco = -o.omega.sin();
    let dso = o.omega.cos();
    let cf = o.f.cos();
    let sf = o.f.sin();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.vx = dv0 * ((o.e + cf) * (-ci * dco * sO - cO * dso) - sf * (dco * cO - ci * dso * sO));
    p.vy = dv0 * ((o.e + cf) * (ci * dco * cO - sO * dso) - sf * (dco * sO + ci * dso * cO));
    p.vz = dv0 * ((o.e + cf) * dco * si - sf * si * dso);

    p
}

pub fn reb_particle_derivative_m_f(G: f64, primary: reb_particle, po: reb_particle) -> reb_particle {
    let o = reb_orbit_from_particle(G, po, primary);
    let mut p = reb_particle::default();
    let dv0 = 0.5 * (G / o.a / (1. - o.e * o.e)).sqrt() / (po.m + primary.m).sqrt();

    let cO = o.Omega.cos();
    let sO = o.Omega.sin();
    let co = o.omega.cos();
    let so = o.omega.sin();
    let dcf = -o.f.sin();
    let dsf = o.f.cos();
    let ci = o.inc.cos();
    let si = o.inc.sin();

    p.vx = dv0 * (dcf * (-ci * co * sO - cO * so) - dsf * (co * cO - ci * so * sO));
    p.vy = dv0 * (dcf * (ci * co * cO - sO * so) - dsf * (co * sO + ci * so * cO));
    p.vz = dv0 * (dcf * co * si - dsf * si * so);

    p
}
