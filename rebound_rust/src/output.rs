//! output.rs — ASCII output routines (from output.c; screenshot and
//! MPI paths excluded, PROFILING excluded like a default C build).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1
//! (c) Hanno Rein, Shangfei Liu and contributors. See crate root.

use crate::tools::{reb_orbit_from_particle, reb_particle_com_of_pair, reb_simulation_error};
use crate::types::*;
use std::io::Write;

use std::f64::consts::PI as M_PI;

/// C's `printf("%- 9f", x)`: `%f` with precision 6, minimum width 9,
/// left-justified, space for the sign of positive numbers.
pub fn fmt_f_space_left9(x: f64) -> String {
    let body = if x.is_sign_negative() {
        format!("{:.6}", x)
    } else {
        format!(" {:.6}", x)
    };
    format!("{:<9}", body)
}

/// C's `printf("%- 9d", n)`.
pub fn fmt_d_space_left9(n: i64) -> String {
    let body = if n < 0 {
        format!("{}", n)
    } else {
        format!(" {}", n)
    };
    format!("{:<9}", body)
}

/// C's `printf("%e", x)`: precision 6, two-or-more exponent digits
/// (UCRT and glibc both print at least two).
pub fn fmt_e(x: f64) -> String {
    let s = format!("{:.6e}", x);
    // Rust prints e.g. "1.234000e5" / "1.234000e-5"; C prints
    // "1.234000e+05" / "1.234000e-05".
    if let Some(pos) = s.rfind('e') {
        let (mantissa, exp) = s.split_at(pos);
        let exp = &exp[1..];
        let (sign, digits) = if let Some(stripped) = exp.strip_prefix('-') {
            ("-", stripped)
        } else {
            ("+", exp)
        };
        if digits.len() < 2 {
            format!("{}e{}0{}", mantissa, sign, digits)
        } else {
            format!("{}e{}{}", mantissa, sign, digits)
        }
    } else {
        s
    }
}

/// output.c `reb_simulation_output_check_phase`.
pub fn reb_simulation_output_check_phase(r: &reb_simulation, interval: f64, phase: f64) -> bool {
    let shift = r.t + interval * phase;
    if (shift / interval).floor() != ((shift - r.dt) / interval).floor() {
        return true;
    }
    // Output at beginning
    if r.t == 0. {
        return true;
    }
    false
}

/// output.c `reb_simulation_output_check`.
pub fn reb_simulation_output_check(r: &reb_simulation, interval: f64) -> bool {
    reb_simulation_output_check_phase(r, interval, 0.)
}

fn now_seconds() -> f64 {
    match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(d) => d.as_secs() as f64 + (d.subsec_micros() as f64) / 1000000.0,
        Err(_) => 0.0,
    }
}

/// output.c `reb_simulation_output_timing` — the heartbeat status line.
pub fn reb_simulation_output_timing(r: &mut reb_simulation, tmax: f64) {
    let N_tot = r.N;
    let temp = now_seconds();
    if r.output_timing_last == -1. {
        r.output_timing_last = temp;
    } else {
        print!("\r");
    }
    print!("N_tot={}  ", fmt_d_space_left9(N_tot as i64));
    if r.OMEGA != 0. {
        print!("t={} [orb]  ", fmt_f_space_left9(r.t * r.OMEGA / 2. / M_PI));
    } else {
        print!("t={}  ", fmt_f_space_left9(r.t));
    }
    print!("dt={}  ", fmt_f_space_left9(r.dt));
    print!("cpu={} [s]  ", fmt_f_space_left9(temp - r.output_timing_last));
    if tmax > 0. {
        print!("t/tmax= {:5.2}%", r.t / tmax * 100.0);
    }
    let _ = std::io::stdout().flush();
    r.output_timing_last = temp;
}

/// output.c `reb_simulation_output_ascii` (appends to file).
pub fn reb_simulation_output_ascii(r: &mut reb_simulation, filename: &str) {
    let of = std::fs::OpenOptions::new().create(true).append(true).open(filename);
    let mut of = match of {
        Ok(f) => f,
        Err(_) => {
            reb_simulation_error(r, "Can not open file.");
            return;
        }
    };
    for i in 0..r.N {
        let p = r.particles[i];
        let _ = writeln!(
            of,
            "{}\t{}\t{}\t{}\t{}\t{}",
            fmt_e(p.x),
            fmt_e(p.y),
            fmt_e(p.z),
            fmt_e(p.vx),
            fmt_e(p.vy),
            fmt_e(p.vz)
        );
    }
}

/// output.c `reb_simulation_output_orbits` (appends to file).
pub fn reb_simulation_output_orbits(r: &mut reb_simulation, filename: &str) {
    let of = std::fs::OpenOptions::new().create(true).append(true).open(filename);
    let mut of = match of {
        Ok(f) => f,
        Err(_) => {
            reb_simulation_error(r, "Can not open file.");
            return;
        }
    };
    let mut com = r.particles[0];
    for i in 1..r.N {
        let o = reb_orbit_from_particle(r.G, r.particles[i], com);
        let _ = writeln!(
            of,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fmt_e(r.t),
            fmt_e(o.a),
            fmt_e(o.e),
            fmt_e(o.inc),
            fmt_e(o.Omega),
            fmt_e(o.omega),
            fmt_e(o.l),
            fmt_e(o.P),
            fmt_e(o.f)
        );
        com = reb_particle_com_of_pair(com, r.particles[i]);
    }
}

/// output.c `reb_simulation_output_velocity_dispersion`.
pub fn reb_simulation_output_velocity_dispersion(r: &mut reb_simulation, filename: &str) {
    let N = r.N;
    // Algorithm with reduced roundoff errors (see wikipedia)
    let mut A = reb_vec3d::default();
    let mut Q = reb_vec3d::default();
    for i in 0..N {
        let Aim1 = A;
        let p = r.particles[i];
        A.x = A.x + (p.vx - A.x) / ((i + 1) as f64);
        if r.OMEGA != 0. {
            A.y = A.y + (p.vy + 1.5 * r.OMEGA * p.x - A.y) / ((i + 1) as f64);
        } else {
            A.y = A.y + (p.vy - A.y) / ((i + 1) as f64);
        }
        A.z = A.z + (p.vz - A.z) / ((i + 1) as f64);
        Q.x = Q.x + (p.vx - Aim1.x) * (p.vx - A.x);
        if r.OMEGA != 0. {
            Q.y = Q.y
                + (p.vy + 1.5 * r.OMEGA * p.x - Aim1.y) * (p.vy + 1.5 * r.OMEGA * p.x - A.y);
        } else {
            Q.y = Q.y + (p.vy - Aim1.y) * (p.vy - A.y);
        }
        Q.z = Q.z + (p.vz - Aim1.z) * (p.vz - A.z);
    }
    let N_tot = N;
    let A_tot = A;
    let mut Q_tot = Q;
    Q_tot.x = (Q_tot.x / (N_tot as f64)).sqrt();
    Q_tot.y = (Q_tot.y / (N_tot as f64)).sqrt();
    Q_tot.z = (Q_tot.z / (N_tot as f64)).sqrt();
    let of = std::fs::OpenOptions::new().create(true).append(true).open(filename);
    let mut of = match of {
        Ok(f) => f,
        Err(_) => {
            reb_simulation_error(r, "Can not open file.");
            return;
        }
    };
    let _ = writeln!(
        of,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        fmt_e(r.t),
        fmt_e(A_tot.x),
        fmt_e(A_tot.y),
        fmt_e(A_tot.z),
        fmt_e(Q_tot.x),
        fmt_e(Q_tot.y),
        fmt_e(Q_tot.z)
    );
}
