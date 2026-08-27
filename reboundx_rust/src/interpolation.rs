//! interpolation.rs — translation of REBOUNDx interpolation.c
//! Interpolate particle parameters from a passed dataset between timesteps in the simulation.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! Author: Stanley A. Baronett, Dan Tamayo.
//!
//! # Misc
//!
//! ======================= ===============================================
//! Authors                 S.A. Baronett, D. Tamayo, N. Ferich
//! Implementation Paper    `Baronett et al., 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.6001B/abstract>`_.
//! Based on                `Press et al., 1992 <https://ui.adsabs.harvard.edu/abs/1992nrca.book.....P/abstract>`_.
//! C Example               :ref:`c_example_parameter_interpolation`
//! Python Example          `ParameterInterpolation.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ParameterInterpolation.ipynb>`_.
//! ======================= ===============================================
//!
//! This isn't an effect that's loaded like the others, but an object that
//! facilitates machine-independent interpolation of parameters that can be
//! shared by both the C and Python versions. See the examples for how to use
//! them.
//!
//! **Effect Parameters**
//!
//! Not applicable. See examples.
//!
//! **Particle Parameters**
//!
//! Not applicable. See examples.
//!
//! # Translation notes
//!
//! * The C's `double* times`, `double* values` and `double* y2` become
//!   `Vec<f64>` owned by `rebx_interpolator`. The C's `y2 == NULL` (set
//!   whenever the interpolation type is not `REBX_INTERPOLATION_SPLINE`) is
//!   represented by an **empty** `y2` vector.
//! * `rebx_create_interpolator` returns the structure by value instead of a
//!   heap pointer; `rebx_free_interpolator` therefore consumes it.
//! * `rebx_malloc`'s only observable behaviour in the C is the error report
//!   on allocation failure, which cannot happen here, so the `rebx` argument
//!   of the constructors is unused. It is retained so the signatures match
//!   the C.
//! * The `klo` bisection cache is kept exactly as in the C: it is stored in
//!   the interpolator and updated in place on every call, so a sequence of
//!   calls with increasing `time` walks the table forward one interval at a
//!   time (and the same holds backwards).

use crate::core::rebx_error;
use crate::types::{rebx_extras, rebx_interpolation_type, rebx_interpolator};

/// Given a monotonic array `x[0..(n-1)]` and any array `y[0..(n-1)]`,
/// sets `y2[0..(n-1)]` with second-order derivatives of the
/// interpolating function at the tabulated points `x[i]`.
/// This routine assumes a "natural" spline, i.e. boundary
/// conditions with zero second derivatives at `y2[0]` and `y2[(n-1)]`.
///
/// Adapted from "Numerical Recipes for C," 2nd Ed., §3.3, p. 115.
fn rebx_spline(x: &[f64], y: &[f64], n: i32, y2: &mut [f64]) {
    // C: double p, qn, sig, un, u[n];  (u is a VLA)
    let mut u: Vec<f64> = vec![0.0; n as usize];

    y2[0] = 0.;
    u[0] = 0.0; // lower boundary condition is set to "natural"
    for i in 1..(n - 1) {
        let i = i as usize;
        // the decomposition loop of the tridiagonal algorithm.
        // y2 and u are used for temporary storage of the decompsed factors.
        let sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        let p = sig * y2[i - 1] + 2.;
        y2[i] = (sig - 1.) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6. * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    let qn = 0.;
    let un = 0.; // upper boundary condition is set to "natural"
    let nm1 = (n - 1) as usize;
    let nm2 = (n - 2) as usize;
    y2[nm1] = (un - qn * u[nm2]) / (qn * y2[nm2] + 1.);
    for k in (0..=(n - 2)).rev() {
        // backsubstitution loop of tridiagonal alg.
        let k = k as usize;
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
}

/// Given a monotonic array `xa[0..(n-1)]`, any array `ya[0..(n-1)]`, an array
/// of second derivatives `y2a[0..(n-1)]` outputted from `rebx_spline` above,
/// and a value of `x`, this returns a cubic-spline interpolated value `y`.
/// "Splint" comes from spl(ine)-int(erpolation).
///
/// Adapted from "Numerical Recipes for C," 2nd Ed., §3.3, p. 116.
fn rebx_splint(
    rebx: &mut rebx_extras,
    xa: &[f64],
    ya: &[f64],
    y2a: &[f64],
    x: f64,
    klo: &mut i32,
    n: i32,
) -> f64 {
    // C: double h, b, a;

    // since calls are generally sequential, find and update place for current
    // and future calls
    if xa[*klo as usize] > x {
        // backward case
        while xa[(*klo - 1) as usize] > x {
            *klo = *klo - 1;
        }
        if xa[(*klo - 1) as usize] <= x {
            *klo = *klo - 1; // back one more
        }
    } else {
        // forward case
        while xa[(*klo + 1) as usize] <= x && *klo + 1 != n - 1 {
            *klo = *klo + 1;
        }
    }
    let h = xa[(*klo + 1) as usize] - xa[*klo as usize];
    if h == 0.0 {
        // xa's must be distinct
        rebx_error(rebx, "Cubic spline run-time error...\n");
        rebx_error(rebx, "Bad xa input to routine splint\n");
        rebx_error(rebx, "...now exiting to system...\n");
        return 0.;
    }
    let a = (xa[(*klo + 1) as usize] - x) / h;
    let b = (x - xa[*klo as usize]) / h;
    // evaluate cubic spline
    a * ya[*klo as usize]
        + b * ya[(*klo + 1) as usize]
        + ((a * a * a - a) * y2a[*klo as usize] + (b * b * b - b) * y2a[(*klo + 1) as usize])
            * (h * h)
            / 6.
}

/// Takes an array of times and corresponding array of values and returns a
/// structure that allows interpolation of values at arbitrary times.
///
/// * `rebx` — the REBOUNDx extras instance.
/// * `Nvalues` — length of times and values arrays (must be equal for both).
/// * `times` — array of times at which the corresponding values are supplied.
/// * `values` — array of values at each corresponding time.
/// * `interpolation` — enum specifying the interpolation method.
///
/// Returns a `rebx_interpolator`. Call `rebx_interpolate` to get values.
///
/// C: returns a `struct rebx_interpolator*` from `rebx_malloc`; here the
/// structure is returned by value.
pub fn rebx_create_interpolator(
    rebx: &mut rebx_extras,
    Nvalues: i32,
    times: &[f64],
    values: &[f64],
    interpolation: rebx_interpolation_type,
) -> rebx_interpolator {
    // C: struct rebx_interpolator* interp = rebx_malloc(rebx, sizeof(*interp));
    // The malloc'd block is uninitialized; rebx_init_interpolator sets every
    // field, so a zeroed structure is an exact stand-in here.
    let mut interp = rebx_interpolator {
        interpolation: rebx_interpolation_type::REBX_INTERPOLATION_NONE,
        times: Vec::new(),
        values: Vec::new(),
        Nvalues: 0,
        y2: Vec::new(),
        klo: 0,
    };
    rebx_init_interpolator(rebx, &mut interp, Nvalues, times, values, interpolation);
    interp
}

/// Initializes an already-allocated `rebx_interpolator` (C:
/// `rebx_init_interpolator`).
pub fn rebx_init_interpolator(
    rebx: &mut rebx_extras,
    interp: &mut rebx_interpolator,
    Nvalues: i32,
    times: &[f64],
    values: &[f64],
    interpolation: rebx_interpolation_type,
) {
    // `rebx` is only consulted by rebx_malloc in the C, on the allocation
    // failure path that cannot occur here.
    let _ = &rebx;
    interp.Nvalues = Nvalues;
    interp.interpolation = interpolation;
    // C: calloc + memcpy of Nvalues doubles each.
    interp.times = times[..Nvalues as usize].to_vec();
    interp.values = values[..Nvalues as usize].to_vec();
    interp.y2 = Vec::new(); // C: interp->y2 = NULL;
    interp.klo = 0;
    if interpolation == rebx_interpolation_type::REBX_INTERPOLATION_SPLINE {
        // C: rebx_malloc(rebx, Nvalues*sizeof(*interp->y2)); every element is
        // written by rebx_spline below.
        let mut y2 = vec![0.0; Nvalues as usize];
        rebx_spline(&interp.times, &interp.values, interp.Nvalues, &mut y2);
        interp.y2 = y2;
    }
}

/// Frees the arrays owned by a `rebx_interpolator` (C:
/// `rebx_free_interpolator_pointers`). Here the owned `Vec`s are simply
/// emptied; the memory is released by `Vec`'s own `Drop`.
pub fn rebx_free_interpolator_pointers(interpolator: &mut rebx_interpolator) {
    interpolator.times = Vec::new();
    interpolator.values = Vec::new();
    if !interpolator.y2.is_empty() {
        // C: if (interpolator->y2 != NULL) free(interpolator->y2);
        interpolator.y2 = Vec::new();
    }
}

/// Frees the memory for a `rebx_interpolator` structure (C:
/// `rebx_free_interpolator`). Consumes the structure, since here it is owned
/// by value rather than by pointer.
pub fn rebx_free_interpolator(mut interpolator: rebx_interpolator) {
    rebx_free_interpolator_pointers(&mut interpolator);
    // C: free(interpolator); — the structure is dropped on return.
}

/// Interpolate value at arbitrary times.
///
/// Need to first call `rebx_create_interpolator` with an array of times and
/// corresponding values to interpolate between. See the parameter
/// interpolation examples.
///
/// * `rebx` — the REBOUNDx extras instance.
/// * `interpolator` — the `rebx_interpolator` structure to interpolate from
///   (taken by `&mut` because the C updates the `klo` cache in place).
/// * `time` — time at which to interpolate value.
///
/// Returns the interpolated value at the passed time.
///
/// C's `default:` branch of the switch ("Interpolation option not supported")
/// is unreachable here: `rebx_interpolation_type` is a Rust enum, so the two
/// cases below are exhaustive.
pub fn rebx_interpolate(
    rebx: &mut rebx_extras,
    interpolator: &mut rebx_interpolator,
    time: f64,
) -> f64 {
    match interpolator.interpolation {
        rebx_interpolation_type::REBX_INTERPOLATION_NONE => {
            0. // UPDATE
        }
        rebx_interpolation_type::REBX_INTERPOLATION_SPLINE => {
            // interpolate at passed time
            let mut klo = interpolator.klo;
            let result = rebx_splint(
                rebx,
                &interpolator.times,
                &interpolator.values,
                &interpolator.y2,
                time,
                &mut klo,
                interpolator.Nvalues,
            );
            interpolator.klo = klo;
            result
        }
    }
}
