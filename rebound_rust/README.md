# rebound_rust (`rebound_rs`) + reboundx_rust (`reboundx_rs`)

Pure-Rust translations of two astronomy libraries:

- **[REBOUND](https://github.com/hannorein/rebound) 5.1.1** — an open-source
  multi-purpose N-body code by **Hanno Rein** and collaborators. It calculates
  how objects move under gravity: planets, moons, ring particles, star clusters.
- **[REBOUNDx](https://github.com/dtamayo/reboundx) 5.1.0** — "REBOUND
  eXtras", a library by **Dan Tamayo**, Hanno Rein and collaborators that adds
  extra physics to a REBOUND simulation: general-relativistic precession,
  tides, spin evolution, radiation forces, migration, and more.

Both are translated line-for-line into safe Rust with **zero `unsafe`, zero
external dependencies (std only), and zero build warnings**, keeping the exact
C function names so existing C code reads across almost unchanged.

This is the **macOS / Apple Silicon edition** of the port. A sister edition
exists for Windows 11; where the two platforms measured different things, this
page reports the macOS numbers and labels the Windows ones as the sister
port's record.

**All of the science, and all of the credit, belongs to the original authors.**
This is a translation of their work. See [Attribution and how to cite](#attribution-and-how-to-cite).

---

## Is it actually the same?

Yes — verified bit-for-bit against the original C, compiled with Apple's
`clang` compiler (the one that comes with the free Xcode Command Line Tools),
with both builds run side by side on the same Apple Silicon Mac (an M5 Max
running macOS 26, Apple clang 21, Rust 1.94):

| Test | Result |
|---|---|
| All 21 maths functions (`sin`, `cos`, `tan`, `atan2`, `sqrt`, `fmod`, `exp`, `log`, `cbrt`, `pow`, …), 200,000 samples each | **all bit-identical — including `pow`** |
| 63 REBOUND integrator configurations, 500 steps each | **all bit-identical** |
| Shearing sheet: seed 42, 1,482 particles, 400 steps | **byte-identical dumps, equal SHA-256** (`418c864d…`) |
| All 65 orbital-derivative functions | **130/130 outputs bit-identical** |
| Frequency analysis (MFT, FMFT, FMFT2) | **bit-identical** |
| MEGNO / variational equations (`movetocom_var`) | **bit-identical** |
| Simulationarchive files written by C, continued in Rust — and vice versa | **6/6 directions bit-identical** |
| Web server: blob served by Rust, loaded by the C build | **bit-identical state** |
| REBOUNDx `tides_spin` examples (3 of them), short and long runs | **all 6 runs bit-identical, spins included** |
| REBOUNDx binary files written by C, read by Rust — and vice versa | **round-trips identical, 25/25 checks passed** |
| Automated test suites: 394 (REBOUND) + 137 (REBOUNDx) | **531 pass, 0 fail** |

"Bit-identical" means all 64 bits of every number match — not "agrees to 10
decimal places". A SHA-256 is a 64-character digital fingerprint of a file:
two files with the same fingerprint have exactly the same bytes. The REBOUNDx
binary files are both 10,784 bytes, and the only bytes that differ between the
C file and the Rust file are the 26 bytes of the version-stamp ("githash")
header field — a documented, deliberate deviation; every particle, parameter
and effect byte matches, and each side reads the other's file perfectly.

Full details, and every command needed to reproduce all of this, are in
**`rebound_rust.md`** (also typeset as `rebound_rust.pdf`).

### One platform quirk: `rand_r` on macOS

The shearing-sheet and MEGNO tests start from *random* initial conditions
made reproducible by a **seed** — a starting number that makes the "random"
sequence come out the same every run. REBOUND's C source carries its own copy
of the random-number routine `rand_r` (borrowed from glibc, the GNU C library
used on Linux) — but only switches it on when building for Windows. On macOS
that switch is off, so a plain C build silently uses Apple's own `rand_r`,
which is a *different* generator: same seed 42, different random stream, and
the C reference then builds a different shearing sheet (1,441 particles
instead of 1,482) — the comparison dies at particle 1 through no fault of
either port.

The fix is a **shim** — a small piece of code slotted in to make two things
fit: `porttest/macos_shim/rand_r_glibc.c` contains the same glibc algorithm
compiled as its own object file and linked into every C test harness, so the
C reference uses the glibc generator on macOS too. The upstream C source is
untouched. The Rust port implements the glibc algorithm on every platform, so
it needs no shim. (This is the macOS twin of the Windows edition's "VLA shim",
which Microsoft's compiler needed for the opposite kind of reason; Apple's
clang compiles those files unmodified.)

### The one known difference — there isn't one on macOS

On this platform, **no tested maths function differs**. All 21 functions —
`pow` (raise-to-a-power) included — are bit-identical between the C build and
the Rust build across 200,000 sample inputs each, because both clang-compiled
C and Rust resolve those calls to Apple's maths library.

The Windows sister port's record is different, and worth keeping for history:
there, `pow` was the single function where Rust's implementation and
Microsoft's C library disagreed — about **0.03% of inputs, by at most 2 ULP**
(a ULP, "unit in the last place", is the smallest possible difference between
two adjacent floating-point numbers). That was an artifact of comparing
against Microsoft's compiler, and it does not exist here. Relatedly, the
Windows shearing-sheet dump carried a different SHA-256 (`75bdaab7…`) than the
macOS one: each platform's C-and-Rust pair share *that platform's* maths
library for `exp` and `log`, so the dumps differ *between* platforms while
C-versus-Rust agreement *within* each platform — the actual acceptance
criterion — is exact on both.

---

## Quick start

You need two free tools, installed from the Terminal (macOS's command-line
app, in Applications → Utilities):

```bash
# Apple's developer tools: the clang C compiler and the linker Rust uses
xcode-select --install

# Rust itself, via rustup (the official installer from https://rustup.rs)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Suppose the repository is cloned at
`~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0` (`~` is your home folder)
and your own project sits beside it at `~/work/myproject` (create one with
`cargo new myproject`). Then in your project's `Cargo.toml` — the file that
lists what your project depends on:

```toml
[dependencies]
rebound_rs  = { path = "../rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rebound_rust" }
reboundx_rs = { path = "../rustSolveIt_macos-silicon_SUNDIALS_7_8_0/reboundx_rust" }   # only if you want the extra physics
```

A two-body simulation:

```rust
use rebound_rs::*;

fn main() {
    let mut sim = reb_simulation_create();
    let r = &mut sim;
    reb_simulation_set_integrator(r, "whfast");
    r.G = 1.0;
    r.dt = 0.01;

    let mut star = reb_particle::default();
    star.m = 1.0;
    reb_simulation_add(r, star);

    reb_simulation_add_fmt(r, "m a e", &[
        reb_fmt_arg::d(1e-3), reb_fmt_arg::d(1.0), reb_fmt_arg::d(0.1),
    ]);

    reb_simulation_move_to_com(r);
    reb_simulation_integrate(r, 100.0);
    println!("t = {}  E = {}", r.t, reb_simulation_energy(r));
}
```

Run it with `cargo run --release` from inside `~/work/myproject`.

Adding extra physics with REBOUNDx (general-relativistic precession):

```rust
use rebound_rs::*;
use reboundx_rs::*;

rebx_attach(&mut sim);
let gr = rebx_load_force(&mut sim, "gr_potential").unwrap();
rebx_add_force(&mut sim, gr);
if let Some(rebx) = rebx_extras_mut(&mut sim) {
    rebx_set_param_double(rebx, rebx_ap::force(gr), "c", 10065.32);
}
reb_simulation_integrate(&mut sim, 1000.0);
```

There are **17 runnable examples** across the two crates (13 REBOUND +
4 REBOUNDx), each run with

```bash
cargo run --release --example <name>
```

from the crate's folder, and **every one of the 17 has a companion Jupyter
notebook** (an interactive document mixing runnable code with explanation) in
`notebooks/` — including `tides_spin_pseudo`.

---

## Attribution and how to cite

### REBOUND

REBOUND is © **Hanno Rein**, Shangfei Liu, Dan Tamayo, David S. Spiegel,
Daniel Tamayo, Tiger Lu, Pejvak Javaheri, Rishit Dagli, Dave O'Hallaron,
Ernst Hairer and the REBOUND contributors. GPL-3.0-or-later.

From the REBOUND project:

> If you use this code or parts of this code for results presented in a
> scientific publication, we would greatly appreciate a citation. The simplest
> way to find the citations relevant to the specific setup of your REBOUND
> simulation is:
>
> ```python
> sim = rebound.Simulation()
> # -your setup-
> sim.cite()
> ```

### REBOUNDx

REBOUNDx is © **Dan Tamayo**, Hanno Rein and the REBOUNDx contributors.
GPL-3.0-or-later.

> REBOUNDx allows you to easily incorporate additional physics into your
> REBOUND simulations. For an overview of the technical details and some
> practical recommendations, see **Tamayo, Rein, Shi and Hernandez 2019**. The
> paper publication lines up with REBOUNDx version 3.0.0.

**Cite for any use of REBOUNDx:** Tamayo, Rein, Shi & Hernandez 2019, MNRAS
491, 2885 — [ADS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T) ·
[arXiv:1908.05634](https://arxiv.org/abs/1908.05634).
Per-effect attribution is listed at <https://reboundx.readthedocs.io>.

### `sim.cite()` in this port

`sim.cite()` lives in the upstream **Python** package
(`pip install rebound`), which this Rust port does not reimplement. The tables
below stand in for it: find the rows matching what your simulation uses, and
cite those papers. **At minimum, cite Rein & Liu 2012 for REBOUND, and
additionally Tamayo et al. 2019 if you use REBOUNDx.**

#### REBOUND

| If you use... | Cite |
|---|---|
| REBOUND at all | Rein & Liu 2012, A&A 537, A128 — [ADS](https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.128R) |
| SEI / shearing sheet | Rein & Tremaine 2011, MNRAS 415, 3168 — [ADS](https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.3168R) |
| IAS15 | Rein & Spiegel 2015, MNRAS 446, 1424 — [ADS](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.1424R) |
| IAS15 adaptive timestep (default) | Pham, Rein & Spiegel 2024, OJAp 7, 1 — [ADS](https://ui.adsabs.harvard.edu/abs/2024OJAp....7E...1P) |
| WHFast | Rein & Tamayo 2015, MNRAS 452, 376 — [ADS](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..376R); Wisdom & Holman 1991, AJ 102, 1528 |
| SABA, WH kernels | Rein, Tamayo & Brown 2019, MNRAS 489, 4632 — [ADS](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R) |
| WHFast512 | Javaheri, Rein & Tamayo 2023, OJAp 6, 29 — [ADS](https://ui.adsabs.harvard.edu/abs/2023OJAp....6E..29J) |
| JANUS | Rein & Tamayo 2018, MNRAS 473, 3351 — [ADS](https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3351R) |
| MERCURIUS | Rein et al. 2019, MNRAS 485, 5490 — [ADS](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R); Chambers 1999, MNRAS 304, 793 |
| EOS | Rein 2020, MNRAS 492, 5413 — [ADS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5413R) |
| TRACE | Lu, Hernandez & Rein 2024, MNRAS 533, 3708 — [ADS](https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.3708L); Hernandez & Dehnen 2023, MNRAS 522, 4639 |
| BS / ODE framework | Rein & Liu 2012; implementation follows Hairer, Nørsett & Wanner 1993 §II.9 via [Hipparchus](https://hipparchus.org) (© 2004 Ernst Hairer) |
| Simulationarchive | Rein & Tamayo 2017, MNRAS 467, 2377 — [ADS](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R) |
| Variational equations / MEGNO | Rein & Tamayo 2016, MNRAS 459, 2275 — [ADS](https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2275R) |
| Frequency analysis | Šidlichovský & Nesvorný 1996, CeMDA 65, 137; Laskar 1988, A&A 198, 341; after David Nesvorný's [FMFT](https://www2.boulder.swri.edu/~davidn/fmft/fmft.html) |
| Tree gravity | Rein & Liu 2012; Barnes & Hut 1986, Nature 324, 446 |

#### REBOUNDx

| If you use the effect... | Cite (in addition to Tamayo et al. 2019) |
|---|---|
| `tides_spin` | Lu, Hernandez & Rein 2023, MNRAS 526, 66; Eggleton, Kiseleva & Hut 1998, ApJ 499, 853; Hut 1981, A&A 99, 126 |
| `gr`, `gr_full`, `gr_potential` | Tamayo et al. 2019 §2; Anderson et al. 1975; Nobili & Roxburgh 1986; Newhall, Standish & Williams 1983 |
| `tides_constant_time_lag` | Hut 1981, A&A 99, 126; Bolmont et al. 2015, A&A 583, A116 |
| `radiation_forces` | Burns, Lamy & Soter 1979, Icarus 40, 1 |
| `modify_orbits_forces`, `modify_orbits_direct` | Papaloizou & Larwood 2000, MNRAS 315, 823; Kominami & Ida 2002, Icarus 157, 43 |
| `type_I_migration` | Cresswell & Nelson 2008, A&A 482, 677; Pichierri et al. 2018, CeMDA 130, 54 |
| `exponential_migration` | Hermosillo Ruiz, Lau & Malhotra 2023, MNRAS |
| `gravitational_harmonics` (J2, J4) | standard; see the REBOUNDx documentation |
| `stochastic_forces` | Rein & Papaloizou 2009, A&A 497, 595 |
| `yarkovsky_effect` | Veras, Higuchi & Ida 2019, MNRAS 485, 708 |
| `lense_thirring` | Park et al. 2017, AJ 153, 121 |
| `gas_dynamical_friction` | Ostriker 1999, ApJ 513, 252; Kim & Kim 2007, ApJ 665, 432 |
| `gas_damping_timescale` | Kominami & Ida 2002, Icarus 157, 43 |
| `tides_dynamical` | Press & Teukolsky 1977, ApJ 213, 183 |
| `central_force`, `track_min_distance`, `modify_mass`, `integrate_force`, `interpolation` | Tamayo et al. 2019 |

BibTeX for all of these is collected upstream at
<https://github.com/hannorein/rebound#papers>,
<https://rebound.readthedocs.io/en/latest/citations/> and
<https://reboundx.readthedocs.io>.

### A note on the upstream projects' AI policy

REBOUND's README states:

> REBOUND is a labour of love, created by people. Please refrain from
> submitting issues or pull requests that have been generated by an LLM or
> other fully-automated tools. […] You may of course use AI assistants for your
> own work with REBOUND. Just don't submit any AI generated code.

This translation was produced with AI assistance, and it is **exactly the
"your own work" case that policy permits**: it is a separate derivative work in
its own repository. Accordingly, **nothing from this port has been or will be
submitted to the upstream REBOUND or REBOUNDx repositories** — no issues, no
pull requests. Please do not file bug reports upstream about this port; the
upstream maintainers did not write it and should not be asked to support it.
Report problems here instead, and verify against the original C before
reporting anything upstream.

---

## Repository layout

Both crates live as **top-level folders inside the
`rustSolveIt_macos-silicon_SUNDIALS_7_8_0` repository**: `rebound_rust/`
(this folder) and `reboundx_rust/` beside it.

| Path | What it is |
|---|---|
| `src/` | the REBOUND translation (29 modules) |
| `examples/` | runnable examples — 13 here + 4 in `reboundx_rust/`, 17 in all, each with a companion Jupyter notebook |
| `notebooks/` | the executed notebooks, one per example (including `tides_spin_pseudo`) |
| `porttest/` | C reference harnesses, their raw-bit output, and the macOS `rand_r` shim in `macos_shim/` |
| `rebound_rust.md` / `.pdf` | **the master document**: complete instructions, provenance and verification for both ports |
| `shearing_sheet_port_test.md` | the shearing-sheet acceptance test, in full |
| `reboundx_port_test.md` | the REBOUNDx `tides_spin` acceptance tests, in full |
| `../reboundx_rust/` | the REBOUNDx translation — a sibling top-level folder in this same repository |
| `../rebound/rebound/`, `../reboundx/` | the original C source trees (REBOUND 5.1.1, REBOUNDx 5.1.0), cloned beside the ports for verification — **not committed** to this repository |

---

## Related projects

The companion crate is `reboundx_rust/` (`reboundx_rs` 5.1.0), the pure-Rust
translation of **REBOUNDx**, which depends on this crate and sits beside it in
the same repository — the same arrangement the C REBOUNDx requires of the C
REBOUND.

This is the macOS / Apple Silicon edition of the port; its Windows 11 sister
edition was produced alongside
[`rustSolveIt_Win11_SUNDIALS_7_8_0`](https://github.com/once-ere/rustSolveIt_Win11_SUNDIALS_7_8_0).
Both editions serve as the physics acceptance test of **rustSolveIt**, a
pure-Rust physics simulator built on a pure-Rust translation of SUNDIALS
7.8.0 — here, the surrounding `rustSolveIt_macos-silicon_SUNDIALS_7_8_0`
repository — and follow the same porting discipline: zero `unsafe`, zero
dependencies, zero warnings, C names preserved, and bit-for-bit verification
against the platform's C reference build.

## License

GPL-3.0-or-later — the same license as REBOUND and REBOUNDx. See `LICENSE`.
Both originals are free software; these translations are and remain free
software.
