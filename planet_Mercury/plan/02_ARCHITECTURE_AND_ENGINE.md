# Planet Mercury Tidal-Locking Plan — Document 2 of 8: Architecture and Engine Use

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science, with every
technical term explained where it first appears. Everything needed to understand this
document is inside this document.
**Status:** PLAN (awaiting approval). Nothing has been built yet. Every claim below
about the engine's real code was verified against the engine's actual source files by
an independent audit agent.

This document specifies **how the software will be built**: the Rust compute crate,
the exact solver calls, the staged-integration design, the output files, the browser
display page, the Jupyter notebook, and how the finished work is mirrored into the
engine's git repository. (The physics equations themselves are restated here wherever
needed, so this document stands alone.)

Words expanded once for this document: **ODE** = ordinary differential equation (a
rule for how fast something changes); **CSV** = comma-separated values (a data table
saved as plain text); **SQL** = Structured Query Language (the standard mini-language
for asking a database questions); **HUD** = heads-up display (on-screen number
readouts); **API** = application programming interface (the set of functions a library
offers); **CLI** = command-line interface; **vendored** = a copy of a library stored
inside the project itself (nothing downloaded); **secular** = astronomers' word for
the slow, orbit-averaged part of a motion, as opposed to its fast wiggles.

---

## 1. The compute engine this project sits on

The engine is the repository
`/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/`
— a pure-Rust physics-simulation project whose **only** numerical-integration backend
is a vendored, read-only, pure-Rust port of SUNDIALS 7.8.0 (the Lawrence Livermore
National Laboratory ODE-solver suite), living in that repository's `sundials_rs/`
folder as its own cargo workspace of seven crates (all version 7.8.0):
`sundials_core`, `cvode_rs`, `cvodes_rs`, `kinsol_rs`, `ida_rs`, `idas_rs`,
`arkode_rs`.

Rules of that engine that this project inherits and will obey:

1. **All time integration goes through SUNDIALS.** No hand-written Euler/Runge-Kutta
   steppers anywhere, ever — not in code, not in tooling, not in display pages
   (a display page only *replays* recorded samples).
2. **Zero `unsafe` Rust, zero external dependencies, zero compiler warnings.** Every
   crate root carries `#![forbid(unsafe_code)]`, `#![deny(warnings)]`, and allows
   `non_snake_case`, `non_camel_case_types`, `non_upper_case_globals`. The lock file
   may list only local crates — nothing downloaded from the internet.
3. **`sundials_rs/` is read-only.** If a needed solver symbol is missing there, the
   build STOPS and reports exactly which symbol is missing — it is never reimplemented
   locally.
4. **Numbers printed for verification use the engine's C-style formatters** —
   `fmt_e(x, prec)` (like C's `%.*e`), `fmt_f`, `fmt_g`, and width variants `fmt_ew`,
   `fmt_fw`, `fmt_gw` — from the module `sundials_core::sundials_utils`. Rust's
   built-in `{:e}` formatting is forbidden (it prints exponents differently from C and
   breaks byte-for-byte reproducibility).
5. **Python tooling (notebooks, display-page bakers, test drivers) is standard-library
   Python only** — that is the engine's own convention for everything outside Rust.

---

## 2. The new Rust crate: `mercury_rs`

### 2.1 Location and independence

`planet_Mercury/mercury_rs/` is a standalone cargo package (its `Cargo.toml` carries a
bare `[workspace]` table so it does not attach itself to any parent workspace). Its
only dependencies are two of the engine's vendored SUNDIALS crates, by relative path:

```toml
[workspace]

[package]
name = "mercury_rs"
version = "0.1.0"
edition = "2021"
license = "BSD-3-Clause"

[dependencies]
sundials_core = { path = "../../rustSolveIt_macos-silicon_SUNDIALS_7_8_0/sundials_rs/crates/sundials_core" }
cvode_rs      = { path = "../../rustSolveIt_macos-silicon_SUNDIALS_7_8_0/sundials_rs/crates/cvode_rs" }
```

(When the finished project is mirrored into the engine repository — Section 8 — those
two paths become `../../sundials_rs/crates/sundials_core` and
`../../sundials_rs/crates/cvode_rs`; that one-line difference is the only difference
between the two copies, and the mirror is re-built and re-tested to prove it.)

### 2.2 Source layout

```
mercury_rs/
├── Cargo.toml
├── src/
│   ├── main.rs        command-line entry: pick a run, run it, self-check,
│   │                  print SUCCESS or FAILURE, exit nonzero on failure
│   ├── lib.rs         crate root (lint headers), module declarations
│   ├── params.rs      every physical constant and per-run configuration, in code
│   ├── kepler.rs      Kepler-equation solver (Newton's method) + true anomaly + radius
│   ├── hut.rs         the five Hut eccentricity polynomials f1..f5 + torque helpers
│   ├── rhs.rs         the 5-ODE right-hand side handed to CVODE
│   ├── driver.rs      CVODE setup, the output loop, stage switching, root-finding on
│   │                  the resonance crossing, restart save/load
│   └── output.rs      CSV writers (using fmt_e / fmt_ew) + the run-manifest writer
└── tests/
    └── physics.rs     analytic-expectation unit tests (listed in plan document 5)
```

Crate root lint headers (mandatory, engine convention):

```rust
#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
```

### 2.3 The five equations the RHS computes

State vector y = [a, e, M, θ, Ω] (orbit size in meters, ovalness, orbit clock angle in
radians, spin angle in radians, spin rate in radians/second). With
n = √(G(M☉+m)/a³) (the mean motion — the orbit's average angular speed),
K = 3·G·M☉²·R⁵·k₂·τ / a⁶ (the tidal-brake strength), C = 0.34·m·R² (Mercury's spin
moment of inertia — its resistance to changes in spin; A < B < C are Mercury's three
principal moments of inertia, and B−A measures how lopsided its equator is), and the
Hut (1981) polynomials

- f₁(e) = (1 + 3e² + (3/8)e⁴) / (1−e²)^(9/2)
- f₂(e) = (1 + (15/2)e² + (45/8)e⁴ + (5/16)e⁶) / (1−e²)⁶
- f₃(e) = (1 + (31/2)e² + (255/8)e⁴ + (185/16)e⁶ + (25/64)e⁸) / (1−e²)^(15/2)
- f₄(e) = (1 + (3/2)e² + (1/8)e⁴) / (1−e²)⁵
- f₅(e) = (1 + (15/4)e² + (15/8)e⁴ + (5/64)e⁶) / (1−e²)^(13/2)

the right-hand side is:

- da/dt = (2K/(m·n·a)) · [Ω·f₂ − n·f₃]
- de/dt = (9K·e/(m·n·a²)) · [(11/18)·Ω·f₄ − n·f₅]
- dM/dt = n
- dθ/dt = Ω
- dΩ/dt = ( T_tri + ⟨T_tidal⟩ ) / C, where
  - ⟨T_tidal⟩ = −K·[Ω·f₁ − n·f₂]  (the tidal brake, always on), and
  - T_tri = −(3/2)·G·M☉·(B−A)/r³ · sin(2(θ − f))  (the "handle" torque on Mercury's
    lopsided figure), **on only in the resonant stage** (Section 4). Here f is the
    true anomaly and r the Sun–Mercury distance, both recovered from (a, e, M) by
    solving Kepler's equation M = E − e·sin E with Newton's method (tolerance
    1×10⁻¹⁴ on E, at most 50 iterations, guarded against non-convergence by returning
    a named error), then
    tan(f/2) = √((1+e)/(1−e))·tan(E/2) and r = a(1−e²)/(1+e·cos f).

Constants (identical to the source specification):
G = 6.67430×10⁻¹¹ m³·kg⁻¹·s⁻², M☉ = 1.98847×10³⁰ kg, m = 3.3011×10²³ kg,
R = 2.4397×10⁶ m, C/(mR²) = 0.34, (B−A)/C = 1.0×10⁻⁴, k₂ = 0.12 (the Love number —
how squishy Mercury is), τ = 100 s (spec-literal tidal time lag) or
τ_eff = 1.0×10⁵ s (the documented 1000× time-compressed value used by the main runs).
Initial conditions: a₀ = 5.790905×10¹⁰ m, e₀ = 0.20563, M₀ = 0, θ₀ = 0 (plus the
per-branch sweep offset), Ω₀ = 1.5×10⁻⁴ rad/s.

Floating-point discipline: each formula is written once, in the order shown above, and
never algebraically "simplified" afterward — floating-point arithmetic is not
associative, and reproducibility depends on a fixed evaluation order.

### 2.4 Exact SUNDIALS calls (verified against the vendored engine's real source)

Plain-words preamble for this section: the vendored solver keeps C's calling style
inside safe Rust. "Handles" (solver objects) are shared reference-counted cells that
every call borrows immutably; constructors return `Option` (Rust's "maybe nothing"
type — a `None` becomes a named error, never a crash); vector data is reached through
a temporary borrow guard that must be released before the next solver call. The
"dense Jacobian" below is a 5×5 table of how sensitively each equation's rate responds
to each of the five tracked numbers — the solver estimates it automatically. The
driver uses exactly these symbols (each verified present in the engine):

```text
sundials_core::sundials_context::SUNContext_Create(SUN_COMM_NULL, &mut opt) -> i32 (0 = ok)
sundials_core::nvector_serial::N_VNew_Serial(5, &ctx) -> Option<N_Vector>
sundials_core::sundials_nvector::N_VGetArrayPointer(&v) -> Option<RefMut<Vec<f64>>>
        (wrapped in local with_data / with_data_mut helpers so the guard is scoped)
cvode_rs::cvode::CVodeCreate(CV_BDF, &ctx) -> Option<CVodeMem>
cvode_rs::cvode::CVodeInit(&cv, rhs, t0, &y) -> i32
cvode_rs::cvode::CVodeSVtolerances(&cv, 1.0e-12, &abstol_vector) -> i32
        abstol_vector = [1.0e-3, 1.0e-6, 1.0e-10, 1.0e-10, 1.0e-14]  (a, e, M, θ, Ω)
sundials_core::sunmatrix_dense::SUNDenseMatrix(5, 5, &ctx) -> Option<SUNMatrix>
sundials_core::sunlinsol_dense::SUNLinSol_Dense(&y, &A, &ctx) -> Option<SUNLinearSolver>
cvode_rs::cvode_ls::CVodeSetLinearSolver(&cv, &ls, Some(&A)) -> i32
        (attaching the dense linear solver gives Newton iteration with the internally
         estimated dense Jacobian — exactly the spec's "NEWTON + direct dense")
cvode_rs::cvode_io::CVodeSetUserData(&cv, Some(Box::new(params))) -> i32
        (user data is the engine's standard Option<Box<dyn Any>>; the RHS is a plain
         fn pointer that downcasts it back — i.e., recovers the typed parameters —
         returning -1 on failure)
cvode_rs::cvode_io::CVodeSetMaxNumSteps(&cv, 500_000_000) -> i32
cvode_rs::cvode_io::CVodeSetMaxStep(&cv, 864000.0) -> i32
cvode_rs::cvode_io::CVodeSetStopTime(&cv, t_end) -> i32
cvode_rs::cvode::CVodeRootInit(&cv, 1, Some(root_fn)) -> i32
        (VERIFIED PRESENT in the vendored engine, with the callback type CVRootFn and
         the result reader CVodeGetRootInfo(&cv, &mut rootsfound). The driver uses one
         root function g(t, y) = Ω − 1.5·n(a) so CVode stops exactly at the resonance
         crossing, returning CV_ROOT_RETURN — used for the restart save and the
         capture timestamp. No sample-scanning fallback is needed.)
cvode_rs::cvode::CVode(&cv, t_out, &y, &mut t_ret, CV_NORMAL) -> i32
        (CV_NORMAL integrates to each requested output time; return 0 = CV_SUCCESS,
         1 = CV_TSTOP_RETURN at the stop time, 2 = CV_ROOT_RETURN at a root — all
         fine; negative = error)
cvode_rs::cvode::CVodeReInit(&cv, t, &y) -> i32   (used at stage handovers and angle
        re-anchoring; solver statistics are harvested BEFORE each ReInit because
        ReInit zeroes the counters)
stats:  CVodeGetNumSteps, CVodeGetNumRhsEvals, CVodeGetNumErrTestFails (out-params)
teardown, in order: CVodeFree(&mut Some(cv)); SUNLinSolFree(Some(ls));
        SUNMatDestroy(A); N_VDestroy(y); N_VDestroy(abstol); SUNContext_Free(&mut Some(ctx))
```

The RHS callback signature (fixed by the engine):

```rust
fn rhs(t: f64, y: &N_Vector, ydot: &N_Vector, user_data: &mut Option<Box<dyn Any>>) -> i32
```

### 2.5 The command-line interface

`cargo run --release -- <subcommand>` (run from inside `mercury_rs/`):

| Subcommand | Run id / directory | What it does |
|---|---|---|
| `print-config` | — | prints every constant, tolerance, and run setting (the notebook shows this) |
| `run-a` | `A_spec_literal` | spec-literal run: τ = 100 s, 10 Myr, secular stage + a short full-model proof segment |
| `run-b` | `B_movie` | movie run: τ_eff = 1.0×10⁵ s, staged; integrates until Ω/n first falls to **1.6** (just above the 3:2 crossing), where it writes the restart state |
| `sweep --branches 64` | `C_sweep` | the phase sweep from the saved restart state; one branches.csv row per branch outcome |
| `run-b-final --branch <k>` | `B_final` | continues the chosen captured branch from the restart to the end of the 10-Myr window (the canonical locked history) |
| `run-d` | `D_high_e` | the e = 0.285 guaranteed-capture encore (crossing window only) |
| `run-e` | `E_seam` | the seam validation: secular vs. full model over the same window at Ω/n = 3 |

Every subcommand ends by printing explicit self-check lines and `SUCCESS` or
`FAILURE` (nonzero exit on failure) — the same pattern as the engine's own
self-checking physics examples. All numeric output goes through `fmt_e`/`fmt_ew`.

---

## 3. Output files: what the crate writes (the notebook's raw material)

Each run writes a directory `planet_Mercury/data/runs/<run_id>/`, where `<run_id>` is
one of the six strings in the table above (`A_spec_literal`, `B_movie`, `C_sweep`,
`B_final`, `D_high_e`, `E_seam`), containing:

**`samples.csv`** — one row per output sample. The column names match the database's
`sample` table exactly, one for one (units baked into the names: `_s` seconds,
`_m` meters, `_rad` radians, `_rad_s` radians/second, `_kgm2s` kg·m²/s of angular
momentum, `_j` joules):

```
t_s,a_m,e,M_rad,theta_rad,Omega_rad_s,n_rad_s,ratio,gamma_rad,P_orb_s,P_rot_s,L_spin_kgm2s,L_orb_kgm2s,L_tot_kgm2s,E_spin_j,E_orb_j,stage
0.000000000000e+00,5.790905000000e+10,2.056300000000e-01,...,S
```

(`ratio` = Ω/n; `gamma_rad` = the resonance angle γ = 2θ − 3M, reduced to (−π, π];
`L_spin_kgm2s` = C·Ω; `L_orb_kgm2s` = m·n·a²·√(1−e²); `L_tot_kgm2s` = their sum;
`E_spin_j` = ½CΩ²; `E_orb_j` = −G·M☉·m/(2a); `stage` = S for secular, R for
resonant-full. All floats are written with `fmt_e(x, 12)` so re-runs are
byte-identical.)

**`events.csv`** — one row per detected event, columns `t_s,event,value`, with the
event vocabulary exactly as the database expects it: `stage_handover`, `cross_5:2`,
`cross_2:1`, `cross_3:2`, `capture_detected`, `reanchor`, `restart_saved` (`value` is
the ratio Ω/n at that moment; for `reanchor` it is the running re-anchor count).

**`branches.csv`** (run C only) — one row per sweep branch, with the exact header
`branch_id,theta_offset_rad,captured,t_outcome_s,final_ratio,canonical` matching the
database's `branch` table (Section 5). Sweep branches do NOT write samples.csv rows —
their trajectories are summarized by this one file, and only the chosen canonical
branch is re-run in full as `B_final` (which does write samples).

**`manifest.json`** — the run's complete configuration echo (every constant, τ used,
compression factor, tolerances, branch offset if any, CVODE step statistics — summed
across the 64 branch integrations for `C_sweep` — the re-anchoring count, and the
crate's SUCCESS/FAILURE verdict) — hand-written JSON via the crate's own tiny
serializer (string concatenation with fixed field order; no external JSON crate).

**`restart.csv`** (run B only) — the exact 5-number state + time at the save point
(when Ω/n first falls to 1.6), written with full `fmt_e(x, 17)` precision, so sweep
branches restart bit-reproducibly.

Output cadence: the secular stage samples every 1,000 simulated years; the resonant
stage every 100 years; plus a dense window (every 2 years) for ±50,000 years around
the detected capture so the libration is beautifully resolved. Every cadence value is
recorded in the manifest.

---

## 4. Staged integration and angle re-anchoring (the two numerical designs)

### 4.1 Stages

- **Stage S (secular):** the handle torque T_tri is omitted from dΩ/dt (physically:
  far from resonance it averages to zero; plan document 1 works out why simulating its
  ~6-hour wiggle for 10 Myr is impossible and pointless). The remaining system is
  smooth, so CVODE takes enormous steps. Used while Ω/n > 2.2.
- **Stage R (resonant-full):** the complete five-equation system, handle torque on.
  Entered when Ω/n first reaches 2.2 (safely above the 2:1 crossing), used through
  capture and to the end. The handover is one `CVodeReInit` at the same (t, y).
- **Seam validation (run E):** starting from the same state at Ω/n = 3, integrate one
  window BOTH ways; the secular drift rates of a, e, Ω must agree to within the
  tolerance budget (the full model shows tiny wiggles around the secular trend; their
  window-average is compared). This proves the stage switch changes nothing physical.

### 4.2 Angle re-anchoring (keeping precision over 10 million years)

θ and M are angles that grow forever (θ reaches ~10¹⁰ radians over the run). A
relative-tolerance solver slowly loses absolute angle precision on such huge numbers.
The fix, applied at output boundaries in stage R: subtract 3π·j from θ and 2π·j from
M, where j is the whole number of completed 2π turns in M (that is, M/(2π) rounded
down to an integer), then `CVodeReInit`. This exact pair of subtractions leaves the
resonance angle γ = 2θ − 3M **mathematically unchanged**
(2·(−3πj) − 3·(−2πj) = 0), so the physics is untouched while the numbers stay small.
The manifest counts how many re-anchorings occurred; a unit test proves γ is preserved
across one re-anchoring to the last bit.

---

## 5. The database (full schema in plan document 4)

The notebook ingests every CSV into one SQLite file,
`planet_Mercury/data/mercury_orbit.sqlite3`, using Python's built-in `sqlite3` module
(standard library — no installation). Five tables: `run` (one row per run:
configuration echo + verdict), `sample` (every samples.csv row, keyed by run),
`event` (every events.csv row), `branch` (one row per sweep branch, keyed by run:
offset + outcome + whether it became the canonical branch), and `target` (the
observed-Mercury verification targets the notebook checks against). The notebook then
demonstrates **retrieval** with seven documented SQL queries — the final state of the
canonical run, the capture time, the libration amplitude decay, the orbital-period
history, the angular-momentum ledger, the sweep's capture fraction, and a
traceability query returning the solver settings that produced the canonical run —
whose printed results feed the acceptance checks.

---

## 6. The browser display page (`gui/mercury_orbit.html`)

One completely self-contained HTML file, following the engine's recorded-video-player
pattern (the engine ships thirteen such pages; the pattern is: all data embedded as
JavaScript constants, plain "vanilla" JavaScript — no frameworks — plus a 2-D canvas,
zero network fetches, dark fixed palette, keyboard + mouse controls). Specifics:

- **Data:** flat columnar arrays (`const T=[...], A=[...], E=[...], TH=[...],
  OM=[...], RATIO=[...], GAMMA=[...], PORB=[...], LSPIN=[...], LORB=[...]`),
  decimated — thinned to fewer samples — by the baker to a few thousand points plus
  the dense capture window, values rounded to what plots can resolve — total file
  size a few megabytes.
- **Left panel — the orbit view:** the Sun at a focus, the orbit ellipse drawn from
  the current (a, e), Mercury positioned from (a, e, M) by the same Kepler relations
  (display math on recorded state — no integration in the page, ever), and a **long-axis
  arrow through Mercury** rotated by θ, so the eye can see the spin ticking three
  times per two orbits after capture. A trail shows recent motion.
- **Right panel — the story plots:** stacked time-series canvases with a shared time
  cursor: Ω/n (log-time; horizontal rules at 1.5 and 1.256), γ libration, P_orb and
  P_rot histories, and the angular-momentum ledger (L_spin falling, L_orb rising,
  L_tot flat). Angles are unwrapped for plotting: when the stored angle jumps from
  just below +π to just above −π, the plot adds 2π so the curve stays continuous.
- **HUD readouts:** simulated year, Ω/n to 8 digits beside "target 1.5", current
  P_orb / P_rot / P_solar in Earth days beside the observed 87.969 / 58.646 / 175.938,
  and the locked/unlocked status.
- **Controls:** Play/Pause (Space), scrub slider, speed select, single-step arrows,
  and a "jump to capture" button.
- **Baker:** `gui/bake_page.py`, standard-library Python, deterministic (no
  timestamps, fixed float formatting), reads the SQLite database, writes the page by
  plain string-template substitution. Re-running the baker on the same database must
  reproduce the page byte-for-byte (a verification gate).

---

## 7. The Jupyter notebook (authoring spec in plan document 3)

`notebook/mercury_tidal_locking.ipynb` — a **Python 3 (ipykernel)** notebook, the
engine project's exact notebook convention (its 109 existing notebooks all work this
way: the notebook's Python cells drive the Rust engine as a subprocess; every code
cell is preceded by an explanatory markdown cell; all prose is self-contained with no
references to other files; the last cell is an interactive save-dialog cell that batch
runners skip). Cell-by-cell content, execution procedure, and the teacher/student
walk-through are fully specified in plan document 3. Companion tools
`notebook/run_notebook.py` (executes every code cell in order, embeds the real
outputs, fails loudly on any error) and `notebook/check_notebook.py` (audits the
structure rules) are standard-library-only Python, modeled on the engine's own
notebook machinery.

Headless rule: the notebook honors the environment variable `MERCURY_NO_BROWSER=1`
(skip opening the browser page — used by batch verification), mirroring the engine's
`POSIM_NO_BROWSER` convention.

---

## 8. Mirroring into the engine repository and pushing (the endgame)

Verified facts this procedure rests on: the engine repository's `origin` remote is
exactly `https://github.com/once-ere/rustSolveIt_macos-silicon_SUNDIALS_7_8_0.git`;
its branch `main` is clean and in sync; its root `Cargo.toml` keeps sub-projects out
of the cargo workspace via one line —
`exclude = ["sundials_rs", "vendor/spec_math", "rebound_rust", "reboundx_rust"]`;
its ignore rules exclude any `target/` build directory at any depth; and its
convention for added work is a root-level `*_PROVENANCE*.md` file.

Procedure (each step verified before the next):

1. Back up the engine root `Cargo.toml` into the engine's `.backups/<date>/` folder
   (the engine's own backup convention; that folder is never committed).
2. Copy the verified `planet_Mercury/` tree (minus any `target/` build directories and
   minus nothing else) to `<engine-repo>/planet_Mercury/`. The tree includes
   `SOURCE_SPECIFICATION.md` — a verbatim copy of the consolidated physics
   specification — so the pushed repository carries its own source document.
3. Apply the single documented difference: in the mirror's `mercury_rs/Cargo.toml`
   (and `mercury_crosscheck/Cargo.toml` if built), the SUNDIALS path-dependency lines
   change to `../../sundials_rs/crates/...` (and the crosscheck's rebound paths to
   `../../rebound_rust` / `../../reboundx_rust`).
4. Edit the engine root `Cargo.toml`: append `"planet_Mercury"` to the `exclude`
   array. Nothing else changes.
5. Re-build and re-test the mirror in place (`cargo build --release`, `cargo test`,
   plus one full run-E re-execution) — proving the mirror is functional, not just
   copied.
6. Re-run the engine's own gates to prove the engine is untouched:
   `cargo build --workspace --all-targets` (warning-free),
   `cargo test --workspace` (622 tests), `bash tools/macos_verify_physics.sh`.
7. Write `<engine-repo>/PLANET_MERCURY_PROVENANCE.md` (self-contained summary of what
   was added, from where, with verification results and the findings register).
8. Stage by explicit path only: `git add planet_Mercury Cargo.toml
   PLANET_MERCURY_PROVENANCE.md`; inspect `git status --porcelain` to confirm only
   intended paths are staged (no `target/`, no strays); commit with an import-style
   message; `git push` (upstream already configured).
9. Post-push verification: `git status -sb` shows `main...origin/main` in sync;
   `git ls-remote origin main` returns the new commit hash.

Safety rails: no `.gitignore` files are created or committed in `planet_Mercury/`
(assignment rule); the outer folder's `p*.txt` files are physically unstageable (they
live outside the repository work tree); `git add -A` is never used (the repository's
own ignore-file comments document a past incident); nothing under the read-only
reference trees is touched.
