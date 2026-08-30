# Mercury Tidal Locking — Complete Instructions for Teachers and Students

**What this is:** the full, self-contained guide to running (and re-creating) the
`planet_Mercury` project: a Jupyter notebook that simulates how the planet Mercury
became tidally locked to the Sun in its strange 3:2 spin-orbit resonance — spinning
exactly three times for every two trips around the Sun. Everything you need is in
this one document.

Words used below, explained once: **ODE** = ordinary differential equation (a rule
for how fast something changes); **CSV** = comma-separated values (a data table saved
as plain text); **SQL** = Structured Query Language (the standard mini-language for
asking a database questions); **BDF** = backward differentiation formulas (the solver
method used, built for stiff problems); **libration** = the gentle rocking of a locked
angle around its resting point — a librating resonance angle IS the lock.

---

## Part A — What you will have when you follow these instructions

A Jupyter notebook, `notebook/mercury_tidal_locking.ipynb`, that on one top-to-bottom
run:

1. **Explains** (inside the notebook itself, in plain English) how the Sun's tides
   slowly braked Mercury's spin, and why Mercury got caught at exactly 3 spins per
   2 orbits.
2. **Computes** that entire history by running a pure-Rust program (`mercury_rs`)
   that integrates five ODEs with the SUNDIALS 7.8.0 CVODE solver (BDF method, Newton
   iteration, dense linear solver — the rustSolveIt engine's only integration
   backend). Python never integrates anything.
3. **Stores** every result in a real database file (`data/mercury_orbit.sqlite3`) and
   **retrieves** results from it with seven documented SQL queries.
4. **Checks itself** against the observed Mercury: a year of 87.969 Earth days, a
   rotation of 58.646 Earth days (exactly 2/3 of the year), and a 175.938-day solar
   day (exactly two Mercury years, sunrise to sunrise). The notebook stops with an
   error if any check fails.
5. **Displays** the whole story as one interactive page in your web browser
   (`gui/mercury_orbit.html`): the orbit with a spinning Mercury and its long-axis
   "handle," a spin/orbit-ratio dial that settles on exactly 1.5, the libration plot,
   and the period and angular-momentum histories, with play/scrub controls.

**What the shipped, verified pipeline actually measured** (your re-run reproduces
these):

| Result | Value |
|---|---|
| Spec-literal braking after 10 Myr (run A) | spin barely moves: 181.4x → 178.9x the orbital rate — the real braking takes ~4.7 billion years, so the main runs use a documented 1000x time compression |
| Capture into 3:2 (canonical run) | at 4.6499 million movie-years (≈ 4.65 billion real-equivalent years) |
| Measured capture odds (64-branch phase sweep) | 10/64 = 15.6% ± 4.5% per crossing (Goldreich–Peale 1966 theory estimates ~7% — same order; the classic formula is approximate) |
| Settled spin/orbit ratio (time-averaged) | 1.5000000498 |
| Final year / rotation / solar day | 87.968 d / 58.642 d / 175.91 d (observed: 87.969 / 58.646 / 175.938; instantaneous samples carry a ~2×10⁻⁴ residual-libration wiggle) |
| Angular-momentum ledger | what the spin lost, the orbit gained: pre-capture total drift 3.9×10⁻¹⁰ (budget 10⁻⁹, where the model conserves exactly); locked-era drift 9.2×10⁻¹² vs the model's own predicted 1.1×10⁻⁹ secular leak (the spec's handle torque has no orbital back-reaction; the recorded value is further suppressed by 64-bit orbit quantization) |
| Libration decay after capture | the notebook's 12-bin table shows the swing falling 6.22 rad → 0.15 rad over the first 1.2 million years after capture (the crate's own check, spanning the whole locked era, sees 4.58 → 0.06 rad) — the lock settling in |
| Independent cross-check | REBOUNDx-port `tides_spin` despin rate agrees with the Hut formula to 0.25% — a separate one-command crate outside the notebook: `cd planet_Mercury/mercury_crosscheck && cargo run --release` prints the agreement and SUCCESS (an upstream REBOUNDx notice about a 4.5.0 tidal-potential fix also prints — expected, not an error) |

## Part B — One-time setup (teachers: once per classroom machine)

Requirements: a Mac with Apple Silicon (this project's engine is the macOS/Apple
Silicon port), the Rust toolchain (`cargo`), `git`, and Python 3. Check all three in
a terminal:

```
cargo --version
git --version
python3 --version
```

If `cargo` is missing, install Rust from https://rustup.rs (one command, no
administrator rights). macOS provides `git` and `python3` after a one-time "install
command line developer tools" prompt — click Install if it appears, wait, re-run.

**Step B1 — get the project.** Either you already have the workspace folder
`/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/` (which
contains `planet_Mercury/` beside the engine folder), or clone the published
repository — from whatever folder you keep projects in:

```
git clone https://github.com/once-ere/rustSolveIt_macos-silicon_SUNDIALS_7_8_0.git
cd rustSolveIt_macos-silicon_SUNDIALS_7_8_0
```

(in the clone, `planet_Mercury/` is a top-level folder inside the repository). The
compute crate finds the SUNDIALS solver by a relative path into the engine that
ships with it — specifically, `planet_Mercury/mercury_rs/Cargo.toml` points two
folders up and into `sundials_rs/crates/` (in the workspace layout that engine is
the nested `rustSolveIt_macos-silicon_SUNDIALS_7_8_0/` folder; in a clone it is
the repository root itself) — keep the folder layout exactly as delivered.

**Step B2 — build the compute program** (5–10 minutes the first time). From the
folder that contains `planet_Mercury/`:

```
cd planet_Mercury/mercury_rs
cargo build --release 2>&1 | tee /tmp/mercury_build.log
```

Zero errors and zero warnings expected. (Every command here uses
`2>&1 | tee <logfile>` so the full output is both shown and saved — a habit worth
teaching.)

**Step B3 — prove the program works before touching the notebook** (still inside
`planet_Mercury/mercury_rs/`):

```
cargo test --release 2>&1 | tee /tmp/mercury_test.log
cargo run --release -- print-config
```

`cargo test` must show a line reporting **15 passed; 0 failed** (the physics
suite; the other harness lines saying "0 passed; 0 failed" are normal empty
harnesses) — every test checks a pencil-and-paper
expectation (Kepler's law, the Hut polynomials, torque signs, exact
angular-momentum closure, the libration frequency…). `print-config` prints every
constant and solver setting; have students read it aloud.

**Step B4 — install JupyterLab** (the viewer only; the notebook's own code cells use
nothing outside Python's standard library). Modern macOS Python refuses system-wide
installs, so use a small private environment — from the folder that contains
`planet_Mercury/`:

```
python3 -m venv planet_Mercury/notebook/.venv
planet_Mercury/notebook/.venv/bin/pip install jupyterlab
```

(One-time, about a minute. Deleting the `.venv` folder uninstalls everything.)

## Part C — Running the notebook (students start here)

**Step C1 — open it.** From the folder that contains `planet_Mercury/`:

```
cd planet_Mercury/notebook
.venv/bin/jupyter lab mercury_tidal_locking.ipynb
```

A browser tab opens with the notebook. If asked to pick a kernel, choose
**Python 3 (ipykernel)**.

**Step C2 — run it, top to bottom.** Click the first cell, then press
**Shift+Enter** repeatedly (each press runs one cell and moves on), or use
Run → Run All Cells. The six simulation cells print progress lines and take minutes
each (the whole notebook is roughly an hour of computing — it integrates millions of
years of physics at 12-digit accuracy, including a 64-branch parallel sweep).

**Step C3 — the finale.** The last computational cell bakes the display page twice,
proves the two bakes byte-identical, and opens `gui/mercury_orbit.html` in your
browser. Press **Space** to play. Watch the spin/orbit dial: it starts near 181 (a
young Mercury spinning every 11.6 hours), falls for millions of years, hesitates at
the 3:2 line — and locks. Use "Jump to capture" and the scrub bar; the libration
plot below shows the resonance angle's circulation turning into a dying rocking —
tidal locking made visible.

Everything is also saved on disk: raw simulation output in
`planet_Mercury/data/runs/`, the database at
`planet_Mercury/data/mercury_orbit.sqlite3`, the page at
`planet_Mercury/gui/mercury_orbit.html`.

**Headless / grading re-run** (no browser, fails loudly on any error) — from
`planet_Mercury/notebook/`:

```
MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_tidal_locking.ipynb
python3 check_notebook.py mercury_tidal_locking.ipynb
```

Both must end happily (`ok … cells` and `all structure rules pass`).

## Part D — What each section of the notebook contains

The cells alternate: a markdown (text) cell explains, then a code cell does exactly
what was explained.

- **§1 What this notebook computes** — the story and the three observed numbers it
  must reproduce.
- **§2 How to run this notebook** — the Shift+Enter instructions, repeated inside
  the notebook so it stands alone.
- **§3 The words used in this notebook** — a one-line glossary (tide, torque,
  resonance, libration, ODE, tolerance, secular…).
- **§4 The physical situation** — the complete model in the notebook's own words:
  point-mass Sun; extended, deformable, lopsided Mercury; the five tracked numbers
  [a, e, M, θ, Ω]; every governing equation including all five Hut eccentricity
  factors; every constant with its meaning; and the honest deviations (the 1000x
  time compression, the two-stage integration, the capture coin-flip) with reasons.
  **§4.3** states the exact solver contract handed to SUNDIALS: CVODE, BDF, Newton,
  dense; relative tolerance 10⁻¹²; absolute tolerances [10⁻³, 10⁻⁶, 10⁻¹⁰, 10⁻¹⁰,
  10⁻¹⁴]; max step 864,000 s.
- **§5 How this notebook talks to the simulator** — the subprocess driver (finds the
  binary, streams progress, insists on the program's own SUCCESS verdict).
- **§6 Building and running the simulation** — one explained cell per run:
  6.1 config echo; 6.2 run A (spec-literal honesty run); 6.3 run B (compressed
  braking to the saved restart at ratio 1.6); 6.4 run C (the 64-branch sweep and the
  measured odds, with error bars); 6.5 run B-final (the first captured branch,
  continued to the locked present); 6.6 run D (guaranteed capture at eccentricity
  0.285 — why ovalness is the secret); 6.7 run E (proof the two integration stages
  agree, at the deliberately non-resonant checkpoint 2.7).
- **§7 Storing the results in a database** — creates the five tables (`run`,
  `sample`, `event`, `branch`, `target`) and bulk-loads every CSV; prints row counts
  and a referential-integrity check.
- **§8 Retrieving results from the database** — seven explained SQL queries: the
  present day; the capture time; the libration decay; the orbital-period history;
  the angular-momentum ledger; the measured capture odds; and the traceability query
  (which solver settings produced the canonical history).
- **§9 The verification gauntlet** — every acceptance target asserted in one place;
  one PASS line each, or the notebook stops.
- **§10 Baking and opening the display page** — deterministic double-bake plus the
  browser launch (skipped under MERCURY_NO_BROWSER=1).
- **§11 What we learned** — the story retold with the measured numbers, plus
  exercises (bigger sweeps; where does guaranteed capture stop?; the excluded
  physics — Jupiter's tugs and Einstein's perihelion correction — is the planned
  second test).
- **§12 Name and save this notebook** — an interactive save-dialog cell (batch
  runners skip it).

## Part E — For teachers: how the notebook was authored (make your own!)

1. **One notebook = one complete lesson.** Every explanation lives inside the
   notebook; never "see another file" — restate instead.
2. **Markdown before code, always** — every code cell has a ≥80-character
   explanation cell above it.
3. **Plain Python 3, standard library only** in code cells (`subprocess`, `sqlite3`,
   `csv`, `json`, `pathlib`, `webbrowser`). JupyterLab is only the viewer.
4. **The Rust program does all the physics.** Python runs it, loads, queries,
   checks, displays. New physics goes into `mercury_rs` (Rust), never the notebook.
5. **Real outputs are part of the notebook.** The shipped file was executed
   top-to-bottom by `run_notebook.py`, so every output you see is a genuine capture;
   `check_notebook.py` audits the structure rules mechanically.
6. **Determinism is sacred.** Same inputs → byte-identical outputs: no timestamps,
   fixed number formatting, a recorded phase for the capture coin-flip. The notebook
   is authored by `build_notebook.py` (also deterministic) — editing that builder
   and re-running it, then `run_notebook.py`, is the full authoring loop for a
   sibling notebook (a Venus? the Moon?).

## Part F — Troubleshooting (verified answers)

| Symptom | Cause and fix |
|---|---|
| `cargo: command not found` | Install Rust from https://rustup.rs, reopen the terminal |
| `python3`/`git` opens an "install command line developer tools" dialog | Normal on a fresh Mac — Install, wait, re-run |
| `pip install` refuses ("externally managed environment") | Use the Part B4 venv commands exactly |
| The driver cell cannot find the binary | Run the Part B2 build, or set `MERCURY_BIN` to the binary's full path |
| Build error about a missing `sundials_core`/`cvode_rs` path | The relative path from `planet_Mercury/mercury_rs/` to the engine's `sundials_rs/crates/` folder was broken by moving folders — restore the shipped layout |
| A cell raises `RuntimeError: … FAILURE` | The Rust program's self-check failed and its output names the exact check; re-run once, then `cargo test --release` and report |
| The browser page does not open | Open `planet_Mercury/gui/mercury_orbit.html` by double-clicking it; under MERCURY_NO_BROWSER=1 no page is expected |
| The notebook takes very long | Expected: ~an hour in release mode; if you built without `--release` it can take many hours — rebuild |
| Sweep reports zero captures | Impossible for the shipped, deterministic pipeline on this platform (it always finds 10/64); when re-creating the pipeline or running on another platform it is a real ~1% possibility — re-run section 6.4 once with `run("sweep", "--branches", "128")` (the documented contingency) |
| Numbers differ from this document in the last digit | They should not (the pipeline is deterministic on this platform); if it persists after a clean re-run, that IS a reportable finding |

## Part G — Test 2: Jupiter and Einstein (the second notebook)

The folder ships a **second, independent notebook**:
`notebook/mercury_test2_jupiter_gr.ipynb`. It adds the two pieces of physics the
first lesson deliberately left out, and it is just as self-contained — every
explanation is restated inside it.

**What's new in the physics.**

1. **Einstein's correction.** Very close to the Sun, gravity is a whisker
   stronger than Newton's law predicts. The orbit stays an ellipse, but the
   ellipse's long axis slowly swings around — 42.98 arcseconds per century, the
   number that made general relativity famous. The notebook makes the simulator
   reproduce it to better than one part in a thousand with everything else
   switched off.
2. **Jupiter's tugs** (in the classic averaged, "secular" form): they turn
   Mercury's ellipse another ~160 arcseconds per century AND make the orbit's
   ovalness breathe up and down by about 0.0045 on a ~809,000-year cycle. Both
   are checked against the pencil-and-paper Laplace–Lagrange values.

**The punchline your students should look for.** With the ellipse itself
turning, a locked Mercury cannot average EXACTLY 1.5 spins per orbit any more —
the lock tracks the moving perihelion, so the settled average is
1.5 + (perihelion turn rate)/(orbital rate) ≈ 1.5000004. The notebook measures
that four-parts-in-ten-million offset and checks it against the prediction:
**the lock follows the precessing ellipse, not the stars.** A genuinely
relativistic fingerprint, resolved by a 12-digit integration your class can run.

**How to run it.** Identical to Parts A–B: same build, same JupyterLab
environment, then open `mercury_test2_jupiter_gr.ipynb` and Shift+Enter top to
bottom (or batch:
`MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_test2_jupiter_gr.ipynb`).
It writes its own database (`data/mercury_test2.sqlite3`) and bakes its own
display page (`gui/mercury_test2.html`) — nothing it does touches the first
lesson's files, so the two notebooks can be re-run in any order.

**One honest difference to point out.** The first lesson's angular-momentum
ledger is deliberately ABSENT from test 2: Jupiter's averaged tugs exchange
angular momentum with Jupiter itself, which the model does not track, so a
two-body "books must balance" check would be checking the wrong law. The
notebook (and the simulator's own output) say so — auditing what a model can
and cannot promise is itself the lesson there.
