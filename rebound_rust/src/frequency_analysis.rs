//! frequency_analysis.rs — Modified Fourier Transform of Laskar (1988)
//! and the Frequency Modified Fourier Transform of Sidlichovsky &
//! Nesvorny (1996), from frequency_analysis.c ((c) 2025 Hanno Rein,
//! based on David Nesvorny's FMFT code). Given a quasi-periodic complex
//! signal X + iY, estimates the frequencies, amplitudes and phases of
//! its decomposition.
//!
//! Deviation note: the final amplitude sort uses Rust's stable
//! `sort_by` where the C uses `qsort`; for distinct amplitudes (the
//! generic case) the resulting order is identical, only exactly-equal
//! amplitudes could tie-break differently.
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use std::f64::consts::PI as M_PI;

/// FMFT nominal precision.
const FMFT_TOL: f64 = 1.0e-10;
/// FMFT overlap exclusion parameter.
const FMFT_NEAR: f64 = 0.;
/// Maximum number of frequencies to remove from signal before giving up.
const FMFT_MAX_REMOVE: i32 = 64;

const TWOPI: f64 = 2. * M_PI;

/// rebound.h `enum REB_FREQUENCY_ANALYSIS_TYPE`.
pub type REB_FREQUENCY_ANALYSIS_TYPE = i32;
pub const REB_FREQUENCY_ANALYSIS_MFT: i32 = 0;
pub const REB_FREQUENCY_ANALYSIS_FMFT: i32 = 1;
pub const REB_FREQUENCY_ANALYSIS_FMFT2: i32 = 2;

/// frequency_analysis.c `reb_frequency_analysis`.
///
/// The output array needs to have space for 3*nfreq values. They are
/// stored as three blocks: frequencies, amplitudes, phases (each nfreq
/// long). The input array needs to be ndata*2 long with values
/// x(0), y(0), x(1), y(1), ... where x(t) and y(t) are the complex time
/// series to be analyzed. Returns 0 on success and a negative value for
/// various errors.
pub fn reb_frequency_analysis(
    output: &mut [f64],
    nfreq: usize,
    minfreq: f64,
    maxfreq: f64,
    type_: REB_FREQUENCY_ANALYSIS_TYPE,
    input: &[f64],
    ndata: usize,
) -> i32 {
    if minfreq >= maxfreq {
        println!("Frequency analysis error: minfreq must be smaller than maxfreq.");
        return -1;
    }
    if nfreq == 0 {
        println!("Frequency analysis error: nfreq must be larger than 0.");
        return -2;
    }
    if (ndata & (ndata.wrapping_sub(1))) != 0 {
        println!("Frequency analysis error: ndata must be power of 2.");
        return -3;
    }
    // (The C NULL-pointer checks for input/output do not apply: slices
    // are always valid.)

    /* ALLOCATION OF VARIABLES */

    let mut xdata = vec![0.0_f64; ndata];
    let mut ydata = vec![0.0_f64; ndata];
    let mut x = vec![0.0_f64; ndata];
    let mut y = vec![0.0_f64; ndata];
    let mut powsd = vec![0.0_f64; ndata];

    let mut freq = vec![0.0_f64; 3 * ((type_ as usize) + 1) * nfreq];
    let mut amp = vec![0.0_f64; 3 * ((type_ as usize) + 1) * nfreq];
    let mut phase = vec![0.0_f64; 3 * ((type_ as usize) + 1) * nfreq];

    let mut f = vec![0.0_f64; nfreq];
    let mut A = vec![0.0_f64; nfreq];
    let mut psi = vec![0.0_f64; nfreq];

    let mut Q = vec![0.0_f64; nfreq * nfreq];
    let mut alpha = vec![0.0_f64; nfreq * nfreq];
    let mut B = vec![0.0_f64; nfreq];

    /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

    for l in 0..=(type_ as usize) {
        if l == 0 {
            /* SEPARATE REAL AND IMAGINARY PARTS */
            for j in 0..ndata {
                xdata[j] = input[j * 2];
                ydata[j] = input[j * 2 + 1];
            }
        } else {
            /* GENERATE THE QUASIPERIODIC FUNCTION COMPUTED BY MFT */
            for i in 0..ndata {
                xdata[i] = 0.;
                ydata[i] = 0.;
                for k in 0..nfreq {
                    xdata[i] += amp[(l - 1) * nfreq + k]
                        * (freq[(l - 1) * nfreq + k] * (i as f64) + phase[(l - 1) * nfreq + k])
                            .cos();
                    ydata[i] += amp[(l - 1) * nfreq + k]
                        * (freq[(l - 1) * nfreq + k] * (i as f64) + phase[(l - 1) * nfreq + k])
                            .sin();
                }
            }
        }

        /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
        window(&mut x, &mut y, &xdata, &ydata, ndata);

        /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
        power(&mut powsd, &x, &y, ndata);

        let mut centerf: f64;

        if l == 0 {
            /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
            let mut frequencies_removed = 0;
            loop {
                centerf = bracket(&powsd, ndata);
                if !(centerf < minfreq || centerf > maxfreq) {
                    break;
                }
                /* IF NO, SUBTRACT IT FROM THE SIGNAL */
                f[0] = golden(centerf, TWOPI / (ndata as f64), &x, &y, ndata);

                {
                    let (a0, p0) = (&mut A[0], &mut psi[0]);
                    amph(a0, p0, f[0], &x, &y, ndata);
                }

                for j in 0..ndata {
                    xdata[j] -= A[0] * (f[0] * (j as f64) + psi[0]).cos();
                    ydata[j] -= A[0] * (f[0] * (j as f64) + psi[0]).sin();
                }

                window(&mut x, &mut y, &xdata, &ydata, ndata);

                power(&mut powsd, &x, &y, ndata);

                frequencies_removed += 1;
                if frequencies_removed > FMFT_MAX_REMOVE {
                    println!(
                        "Frequency analysis error: cannot find frequencies in range [minfreq, maxfreq]."
                    );
                    return -6;
                }
            }
        } else {
            centerf = freq[0];
        }

        /* DETERMINE THE FIRST FREQUENCY */
        f[0] = golden(centerf, TWOPI / (ndata as f64), &x, &y, ndata);

        /* COMPUTE AMPLITUDE AND PHASE */
        {
            let (mut a0, mut p0) = (0.0, 0.0);
            amph(&mut a0, &mut p0, f[0], &x, &y, ndata);
            A[0] = a0;
            psi[0] = p0;
        }

        /* SUBTRACT THE FIRST HARMONIC FROM THE SIGNAL */
        for j in 0..ndata {
            xdata[j] -= A[0] * (f[0] * (j as f64) + psi[0]).cos();
            ydata[j] -= A[0] * (f[0] * (j as f64) + psi[0]).sin();
        }

        /* HERE STARTS THE MAIN LOOP  *************************************/
        Q[0] = 1.;
        alpha[0] = 1.;

        for m in 1..nfreq {
            /* MULTIPLY SIGNAL BY WINDOW FUNCTION */
            window(&mut x, &mut y, &xdata, &ydata, ndata);

            /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
            power(&mut powsd, &x, &y, ndata);

            if l == 0 {
                let mut centerf = bracket(&powsd, ndata);
                f[m] = golden(centerf, TWOPI / (ndata as f64), &x, &y, ndata);

                /* CHECK WHETHER THE NEW FREQUENCY IS NOT TOO CLOSE TO ANY PREVIOUSLY
                DETERMINED ONE */
                let mut nearfreqflag = 0;
                for k in 0..m.saturating_sub(1) {
                    if (f[m] - f[k]).abs() < FMFT_NEAR * TWOPI / (ndata as f64) {
                        nearfreqflag = 1;
                    }
                }

                let mut frequencies_removed = 0;
                /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
                while f[m] < minfreq || f[m] > maxfreq || nearfreqflag == 1 {
                    /* IF NO, SUBTRACT IT FROM THE SIGNAL */
                    f[m] = golden(centerf, TWOPI / (ndata as f64), &x, &y, ndata);

                    {
                        let (mut am, mut pm) = (0.0, 0.0);
                        amph(&mut am, &mut pm, f[m], &x, &y, ndata);
                        A[m] = am;
                        psi[m] = pm;
                    }

                    for j in 0..ndata {
                        xdata[j] -= A[m] * (f[m] * (j as f64) + psi[m]).cos();
                        ydata[j] -= A[m] * (f[m] * (j as f64) + psi[m]).sin();
                    }

                    /* AND RECOMPUTE THE NEW ONE */
                    window(&mut x, &mut y, &xdata, &ydata, ndata);

                    power(&mut powsd, &x, &y, ndata);

                    centerf = bracket(&powsd, ndata);
                    f[m] = golden(centerf, TWOPI / (ndata as f64), &x, &y, ndata);

                    nearfreqflag = 0;
                    for k in 0..m.saturating_sub(1) {
                        if (f[m] - f[k]).abs() < FMFT_NEAR * TWOPI / (ndata as f64) {
                            nearfreqflag = 1;
                        }
                    }
                    frequencies_removed += 1;
                    if frequencies_removed > FMFT_MAX_REMOVE {
                        println!(
                            "Frequency analysis error: cannot find frequencies in range [minfreq, maxfreq]."
                        );
                        return -6;
                    }
                }
            } else {
                /* DETERMINE THE NEXT FREQUENCY */
                f[m] = golden(freq[m], TWOPI / (ndata as f64), &x, &y, ndata);
            }

            /* COMPUTE ITS AMPLITUDE AND PHASE */
            {
                let (mut am, mut pm) = (0.0, 0.0);
                amph(&mut am, &mut pm, f[m], &x, &y, ndata);
                A[m] = am;
                psi[m] = pm;
            }

            /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
            Q[m * nfreq + m] = 1.;
            for j in 0..m {
                let fac = (f[m] - f[j]) * ((ndata as f64) - 1.) / 2.;
                Q[m * nfreq + j] = fac.sin() / fac * M_PI * M_PI / (M_PI * M_PI - fac * fac);
                Q[j * nfreq + m] = Q[m * nfreq + j];
            }

            /* EQUATION (17) */
            for k in 0..m {
                B[k] = 0.;
                for j in 0..k {
                    B[k] += -alpha[k * nfreq + j] * Q[m * nfreq + j];
                }
            }

            /* EQUATION (18) */
            alpha[m * nfreq + m] = 1.;
            for j in 0..m {
                alpha[m * nfreq + m] -= B[j] * B[j];
            }
            alpha[m * nfreq + m] = 1. / alpha[m * nfreq + m].sqrt();

            /* EQUATION (19) */
            for k in 0..m {
                alpha[m * nfreq + k] = 0.;
                for j in k..m {
                    alpha[m * nfreq + k] += B[j] * alpha[j * nfreq + k];
                }
                alpha[m * nfreq + k] = alpha[m * nfreq + m] * alpha[m * nfreq + k];
            }

            /* EQUATION (22) */
            for i in 0..ndata {
                let mut xsum = 0.;
                let mut ysum = 0.;
                for j in 0..=m {
                    let fac = f[j] * (i as f64) + (f[m] - f[j]) * ((ndata as f64) - 1.) / 2.
                        + psi[m];
                    xsum += alpha[m * nfreq + j] * fac.cos();
                    ysum += alpha[m * nfreq + j] * fac.sin();
                }
                xdata[i] -= alpha[m * nfreq + m] * A[m] * xsum;
                ydata[i] -= alpha[m * nfreq + m] * A[m] * ysum;
            }
        }

        /* EQUATION (26) */
        for k in 0..nfreq {
            let mut xsum = 0.;
            let mut ysum = 0.;
            for j in k..nfreq {
                let fac = (f[j] - f[k]) * ((ndata as f64) - 1.) / 2. + psi[j];
                xsum += alpha[j * nfreq + j] * alpha[j * nfreq + k] * A[j] * fac.cos();
                ysum += alpha[j * nfreq + j] * alpha[j * nfreq + k] * A[j] * fac.sin();
            }
            A[k] = (xsum * xsum + ysum * ysum).sqrt();
            psi[k] = ysum.atan2(xsum);
        }

        /* REMEMBER THE COMPUTED VALUES FOR THE FMFT */
        for k in 0..nfreq {
            freq[l * nfreq + k] = f[k];
            amp[l * nfreq + k] = A[k];
            phase[l * nfreq + k] = psi[k];
        }
    }

    /* RETURN THE FINAL FREQUENCIES, AMPLITUDES AND PHASES */
    match type_ {
        REB_FREQUENCY_ANALYSIS_MFT => {
            for k in 0..nfreq {
                output[k] = freq[k];
                output[nfreq + k] = amp[k];
                output[2 * nfreq + k] = phase[k];
            }
        }
        REB_FREQUENCY_ANALYSIS_FMFT => {
            for k in 0..nfreq {
                output[k] = freq[k] + (freq[k] - freq[nfreq + k]);
                output[nfreq + k] = amp[k] + (amp[k] - amp[nfreq + k]);
                output[2 * nfreq + k] = phase[k] + (phase[k] - phase[nfreq + k]);
            }
        }
        REB_FREQUENCY_ANALYSIS_FMFT2 => {
            for k in 0..nfreq {
                output[k] = freq[k];
                let mut fac = freq[nfreq + k] - freq[2 * nfreq + k];
                if (fac / freq[nfreq + k]).abs() > FMFT_TOL {
                    let tmp = freq[k] - freq[nfreq + k];
                    output[k] += tmp * tmp / fac;
                } else {
                    output[k] += freq[k] - freq[nfreq + k];
                }
                output[nfreq + k] = amp[k];
                fac = amp[nfreq + k] - amp[2 * nfreq + k];
                if (fac / amp[nfreq + k]).abs() > FMFT_TOL {
                    let tmp = amp[k] - amp[nfreq + k];
                    output[nfreq + k] += tmp * tmp / fac;
                } else {
                    output[nfreq + k] += amp[k] - amp[nfreq + k];
                }
                output[2 * nfreq + k] = phase[k];
                fac = phase[nfreq + k] - phase[2 * nfreq + k];
                if (fac / phase[nfreq + k]).abs() > FMFT_TOL {
                    let tmp = phase[k] - phase[nfreq + k];
                    output[2 * nfreq + k] += tmp * tmp / fac;
                } else {
                    output[2 * nfreq + k] += phase[k] - phase[nfreq + k];
                }
            }
        }
        _ => {
            println!("REB_FREQUENCY_ANALYSIS_TYPE not implemented.");
        }
    }
    for k in 0..nfreq {
        if output[2 * nfreq + k] < 0.0 {
            output[2 * nfreq + k] += TWOPI;
        }
        if output[2 * nfreq + k] >= 2.0 * M_PI {
            output[2 * nfreq + k] -= TWOPI;
        }
    }

    // SORT THE FREQUENCIES IN DECREASING ORDER OF AMPLITUDE
    sort3(nfreq, output);

    0
}

/// frequency_analysis.c static `window` — Hanning window.
fn window(x: &mut [f64], y: &mut [f64], xdata: &[f64], ydata: &[f64], ndata: usize) {
    for j in 0..ndata {
        let window = (1. - (TWOPI * (j as f64) / ((ndata - 1) as f64)).cos()) * 0.5;
        x[j] = xdata[j] * window;
        y[j] = ydata[j] * window;
    }
}

/// frequency_analysis.c static `power` — rearranges data for the FFT,
/// calls the FFT and returns the power spectral density.
fn power(powsd: &mut [f64], x: &[f64], y: &[f64], ndata: usize) {
    let mut z = vec![0.0_f64; 2 * ndata];
    for j in 0..ndata {
        z[2 * j] = x[j];
        z[2 * j + 1] = y[j];
    }
    four1(&mut z, ndata);
    for j in 0..ndata {
        powsd[j] = z[2 * j] * z[2 * j] + z[2 * j + 1] * z[2 * j + 1];
    }
}

/// frequency_analysis.c static `four1` — in-place FFT (Numerical
/// Recipes style, 0-based); nn must be a power of 2.
fn four1(data: &mut [f64], nn: usize) {
    let n = nn << 1;
    let mut j: usize = 0;
    let mut i: usize = 0;
    while i < n - 1 {
        /* bit-reversal section */
        if j > i {
            data.swap(j, i);
            data.swap(j + 1, i + 1);
        }
        let mut m = n >> 1;
        while m >= 2 && j + 1 > m {
            j -= m;
            m >>= 1;
        }
        j += m;
        i += 2;
    }
    /* Danielson-Lanczos section */
    let mut mmax: usize = 2;
    while n > mmax {
        /* outer ln nn loop */
        let istep = mmax << 1;
        let theta = TWOPI / (mmax as f64); /* initialize */
        let mut wtemp = (0.5 * theta).sin();
        let wpr = -2.0 * wtemp * wtemp;
        let wpi = theta.sin();
        let mut wr = 1.0;
        let mut wi = 0.0;
        let mut m: usize = 0;
        while m < mmax {
            /* two inner loops */
            let mut i = m;
            while i < n {
                j = i + mmax; /* D-L formula */
                let tempr = wr * data[j] - wi * data[j + 1];
                let tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
                i += istep;
            }
            wtemp = wr;
            wr = wtemp * wpr - wi * wpi + wr; /* trig. recurrence */
            wi = wi * wpr + wtemp * wpi + wi;
            m += 2;
        }
        mmax = istep;
    }
}

/// frequency_analysis.c static `bracket` — finds the maximum of the
/// power spectral density.
fn bracket(powsd: &[f64], ndata: usize) -> f64 {
    let mut maxj: usize = 0;
    let mut maxpow = 0.;

    for j in 1..(ndata / 2 - 1) {
        // Changed end from -2 to -1 (comment carried from the C).
        if powsd[j] > powsd[j - 1] && powsd[j] > powsd[j + 1] && powsd[j] > maxpow {
            maxj = j;
            maxpow = powsd[j];
        }
    }

    for j in (ndata / 2 + 1)..(ndata - 1) {
        if powsd[j] > powsd[j - 1] && powsd[j] > powsd[j + 1] && powsd[j] > maxpow {
            maxj = j;
            maxpow = powsd[j];
        }
    }

    if powsd[0] > powsd[1] && powsd[0] > powsd[ndata - 1] && powsd[0] > maxpow {
        maxj = 0;
        maxpow = powsd[0];
    }

    if maxpow == 0. {
        println!("DFT has no maximum ...");
    }

    if maxj < ndata / 2 - 1 {
        -TWOPI * (maxj as f64) / (ndata as f64)
    } else {
        //  maxj > ndata/2-1
        -TWOPI * (((maxj as i64) - (ndata as i64)) as f64) / (ndata as f64)
    }
    /* negative signs and TWOPI compensate for the Numerical Recipes
    definition of the DFT */
}

/// frequency_analysis.c static `golden` — calculates the maximum of a
/// function bracketed by ax, bx and cx.
fn golden(bx: f64, width: f64, xdata: &[f64], ydata: &[f64], n: usize) -> f64 {
    let gold_r = 0.6180339887498948482_f64;
    let gold_c = 1.0 - gold_r;

    let ax = bx - width;
    let cx = bx + width;
    let mut x0 = ax;
    let mut x3 = cx;

    let mut x1: f64;
    let mut x2: f64;
    if (cx - bx).abs() > (bx - ax).abs() {
        x1 = bx;
        x2 = bx + gold_c * (cx - bx);
    } else {
        x2 = bx;
        x1 = bx - gold_c * (bx - ax);
    }

    let mut f1 = phisqr(x1, xdata, ydata, n);
    let mut f2 = phisqr(x2, xdata, ydata, n);

    while (x3 - x0).abs() > FMFT_TOL * (x1.abs() + x2.abs()) {
        if f2 > f1 {
            x0 = x1;
            x1 = x2;
            x2 = gold_r * x1 + gold_c * x3;
            f1 = f2;
            f2 = phisqr(x2, xdata, ydata, n);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = gold_r * x2 + gold_c * x0;
            f2 = f1;
            f1 = phisqr(x1, xdata, ydata, n);
        }
    }

    if f1 > f2 {
        x1
    } else {
        x2
    }
}

/// frequency_analysis.c static `amph` — calculates amplitude and phase.
fn amph(amp: &mut f64, phase: &mut f64, freq: f64, xdata: &[f64], ydata: &[f64], ndata: usize) {
    let mut xphi = 0.;
    let mut yphi = 0.;

    phifun(&mut xphi, &mut yphi, freq, xdata, ydata, ndata);

    *amp = (xphi * xphi + yphi * yphi).sqrt();
    *phase = yphi.atan2(xphi);
}

/// frequency_analysis.c static `phisqr` — square power of phi.
fn phisqr(freq: f64, xdata: &[f64], ydata: &[f64], ndata: usize) -> f64 {
    let mut xphi = 0.;
    let mut yphi = 0.;

    phifun(&mut xphi, &mut yphi, freq, xdata, ydata, ndata);

    xphi * xphi + yphi * yphi
}

/// frequency_analysis.c static `phifun` — computes the function phi.
fn phifun(xphi: &mut f64, yphi: &mut f64, freq: f64, xdata: &[f64], ydata: &[f64], n: usize) {
    let mut xdata2 = vec![0.0_f64; n];
    let mut ydata2 = vec![0.0_f64; n];

    xdata2[0] = xdata[0] / 2.;
    ydata2[0] = ydata[0] / 2.;
    xdata2[n - 1] = xdata[n - 1] / 2.;
    ydata2[n - 1] = ydata[n - 1] / 2.;

    for i in 1..(n - 1) {
        xdata2[i] = xdata[i];
        ydata2[i] = ydata[i];
    }

    let mut nn = n;
    while nn != 1 {
        nn /= 2;
        let c = (-freq * (nn as f64)).cos();
        let s = (-freq * (nn as f64)).sin();

        for i in 0..nn {
            let j = i + nn;
            xdata2[i] += c * xdata2[j] - s * ydata2[j];
            ydata2[i] += c * ydata2[j] + s * xdata2[j];
        }
    }

    *xphi = 2. * xdata2[0] / ((n - 1) as f64);
    *yphi = 2. * ydata2[0] / ((n - 1) as f64);
}

/// frequency_analysis.c static `sort3` — sorts the three output blocks
/// (freq, amp, phase) in decreasing order of amplitude. The C sorts an
/// index array with `qsort` on the amplitude block and permutes; Rust's
/// stable sort produces the identical permutation for distinct keys.
fn sort3(n: usize, output: &mut [f64]) {
    let mut idx: Vec<usize> = (0..n).collect();
    let amps: Vec<f64> = output[n..2 * n].to_vec();
    idx.sort_by(|&a, &b| {
        let (aa, ba) = (amps[a], amps[b]);
        if aa > ba {
            std::cmp::Ordering::Greater
        } else if aa < ba {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    });

    for block in 0..3 {
        let wksp: Vec<f64> = output[block * n..(block + 1) * n].to_vec();
        for j in 0..n {
            output[block * n + j] = wksp[idx[n - j - 1]];
        }
    }
}
