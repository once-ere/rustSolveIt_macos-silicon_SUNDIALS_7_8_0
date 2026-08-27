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

**All of the science, and all of the credit, belongs to the original authors.**
This is a translation of their work. See [Attribution and how to cite](#attribution-and-how-to-cite).

---

## Is it actually the same?

Yes — verified bit-for-bit against the original C compiled with Microsoft's
C compiler on the same machine:

| Test | Result |
|---|---|
| 63 REBOUND integrator configurations, 500 steps each | **all bit-identical** |
| Shearing sheet: 1,482 particles, 400 steps, 102,533 collisions | **byte-identical dumps, equal SHA-256** |
| All 65 orbital-derivative functions | **130/130 outputs bit-identical** |
| Frequency analysis (MFT, FMFT, FMFT2) | **bit-identical** |
| Simulationarchive files written by C, continued in Rust — and vice versa | **bit-identical continuations** |
| Web server: blob served by Rust, loaded by the C build | **bit-identical state** |
| REBOUNDx `tides_spin` examples (3 of them), short and long runs | **all 6 runs bit-identical, spins included** |
| REBOUNDx binary files written by C, read by Rust — and vice versa | **round-trips identical** |
| Automated test suites: 394 (REBOUND) + 137 (REBOUNDx) | **531 pass, 0 fail** |

"Bit-identical" means all 64 bits of every number match — not "agrees to 10
decimal places". Full details, and every command needed to reproduce it, are in
**`rebound_rust.md`** (also typeset as `rebound_rust.pdf`).

### The one known difference

`pow` (raise-to-a-power) is the single maths function where Rust's own
implementation and Microsoft's C library disagree: about **0.03% of inputs, by
at most 2 ULP** (the smallest representable difference). Every other function
tested — `sin`, `cos`, `tan`, `atan2`, `sqrt`, `fmod`, `exp`, `log`, `cbrt` —
is bit-identical. This affects only code paths that call `pow` at run time.
It is characterised precisely in `rebound_rust.md`.

---

## Quick start

You need Rust with the MSVC toolchain (get it from <https://rustup.rs>).

```toml
[dependencies]
rebound_rs  = { path = "path/to/rebound_rust" }
reboundx_rs = { path = "path/to/reboundx_rust" }   # only if you want the extra physics
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

Every example ships with a companion Jupyter notebook in `notebooks/`.

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

| Path | What it is |
|---|---|
| `src/` | the REBOUND translation (29 modules) |
| `examples/` | runnable examples, one Jupyter notebook each in `notebooks/` |
| `porttest/` | C reference harnesses and their raw-bit output |
| `rebound_rust.md` / `.pdf` | **the master document**: complete instructions, provenance and verification for both ports |
| `shearing_sheet_port_test.md` | the shearing-sheet acceptance test, in full |
| `reboundx_port_test.md` | the REBOUNDx `tides_spin` acceptance tests, in full |
| [`once-ere/reboundx_rust`](https://github.com/once-ere/reboundx_rust) | the REBOUNDx translation — its own repository, cloned as a sibling of this one |

---

## Related projects

Its companion is [`once-ere/reboundx_rust`](https://github.com/once-ere/reboundx_rust) (`reboundx_rs` 5.1.0), the pure-Rust translation of **REBOUNDx**, which depends on this crate and expects to be cloned beside it — the same arrangement the C REBOUNDx requires of the C REBOUND.

This port was produced as the physics acceptance test of
[`rustSolveIt_Win11_SUNDIALS_7_8_0`](https://github.com/once-ere/rustSolveIt_Win11_SUNDIALS_7_8_0),
a pure-Rust physics simulator built on a pure-Rust translation of SUNDIALS
7.8.0, and follows the same porting discipline: zero `unsafe`, zero
dependencies, zero warnings, C names preserved, and bit-for-bit verification
against the platform's C reference build.

## License

GPL-3.0-or-later — the same license as REBOUND and REBOUNDx. See `LICENSE`.
Both originals are free software; these translations are and remain free
software.
