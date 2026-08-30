# Planet Mercury Tidal-Locking Plan — Document 3 of 8: The Notebook, and the Complete Instructions for Teachers and Students

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science. Everything
needed is inside this document — no other file must be consulted to follow it.
**Status:** PLAN (awaiting approval). This document is the full draft of the
instructions file that will ship with the finished project (as
`INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md`); during the build, every step below is
verified by actually performing it, and the shipped copy is updated with the real
captured outputs.

Words expanded once for this document: **ODE** = ordinary differential equation (a
rule for how fast something changes); **CSV** = comma-separated values (a data table
saved as plain text); **SQL** = Structured Query Language (the standard mini-language
for asking a database questions); **BDF** = backward differentiation formulas (the
solver method used — built for problems mixing fast wiggles with very slow drifts).

---

## Part A — What you will have when you follow these instructions

A Jupyter notebook, `mercury_tidal_locking.ipynb`, that on one top-to-bottom run:

1. **Explains** (in plain English, inside the notebook itself) how the Sun's tides
   slowly braked Mercury's spin, and why Mercury got caught spinning exactly three
   times per two orbits (the "3:2 spin-orbit resonance" — Mercury's form of tidal
   locking).
2. **Computes** that entire history by running a pure-Rust program (`mercury_rs`)
   that integrates five ODEs with the SUNDIALS CVODE solver — the same solver engine
   used by the whole rustSolveIt project. Nothing else does any integration.
3. **Stores** every result in a real database file (`mercury_orbit.sqlite3`) and then
   **retrieves** results from it with documented SQL queries.
4. **Checks itself**: the final simulated Mercury must show a year of 87.969 Earth
   days, a rotation period of 58.646 Earth days (exactly 2/3 of the year), and a
   sunrise-to-sunrise "solar day" of 175.938 Earth days (exactly two Mercury years).
   The notebook stops with an error if any check fails.
5. **Displays** the whole story as an interactive page in your web browser: the orbit,
   a spinning Mercury with a visible long axis, a spin/orbit-ratio dial that settles
   on exactly 1.5, and history plots of the orbital period and the angular momenta.

---

## Part B — One-time setup (teachers: do this once per classroom machine)

Requirements: a Mac with Apple Silicon (this project's engine is the macOS/Apple
Silicon port), the Rust toolchain (`cargo`), git, and Python 3. Check all three in a
terminal:

```
cargo --version
git --version
python3 --version
```

If `cargo` is missing, install Rust from https://rustup.rs (one command, no
administrator rights needed). macOS provides `git` and `python3` after a one-time
prompt to install the "command line developer tools" — if a dialog pops up, click
Install, wait for it to finish, then re-run the command.

Step B1 — get the project. Either you already have the workspace folder

```
/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/
```

(which contains `planet_Mercury/` beside the engine folder), or you clone the
published repository — from whatever folder you keep projects in:

```
git clone https://github.com/once-ere/rustSolveIt_macos-silicon_SUNDIALS_7_8_0.git
cd rustSolveIt_macos-silicon_SUNDIALS_7_8_0
```

(in the clone, `planet_Mercury/` is a top-level folder inside the repository). The
compute crate finds the SUNDIALS solver by a relative path into the engine folder
that ships with it — so keep the folder layout exactly as delivered; do not move
`planet_Mercury/` somewhere else on its own.

Step B2 — build the compute program once (5–10 minutes the first time). From the
folder that contains `planet_Mercury/` (the workspace folder, or the cloned
repository folder):

```
cd planet_Mercury/mercury_rs
cargo build --release 2>&1 | tee /tmp/mercury_build.log
```

The build must end with zero errors and zero warnings. (Every terminal command in
this project is written with `2>&1 | tee <logfile>` so the full output is both shown
and saved — a project habit worth teaching.)

Step B3 — prove the program works before touching the notebook (still inside
`planet_Mercury/mercury_rs/`):

```
cargo test 2>&1 | tee /tmp/mercury_test.log
cargo run --release -- print-config
```

`cargo test` must report every test passing. `print-config` prints every physical
constant and solver setting the simulation will use — have students read it aloud:
the Sun's mass, Mercury's mass and radius, the tidal "squishiness" k₂ = 0.12, the
solver's relative tolerance 1.0×10⁻¹², and so on.

Step B4 — install JupyterLab (the notebook viewer; the notebook's code cells
themselves use nothing outside Python's standard library). Modern macOS Python
refuses system-wide pip installs (it is "externally managed"), so use a small
private environment — from the folder that contains `planet_Mercury/`:

```
python3 -m venv planet_Mercury/notebook/.venv
planet_Mercury/notebook/.venv/bin/pip install jupyterlab
```

(One-time, ~a minute. Nothing outside this folder is touched; deleting the `.venv`
folder uninstalls everything.)

---

## Part C — Running the notebook (students start here)

Step C1 — open the notebook. From the folder that contains `planet_Mercury/`:

```
cd planet_Mercury/notebook
.venv/bin/jupyter lab mercury_tidal_locking.ipynb
```

A browser tab opens showing the notebook. If asked to pick a kernel, choose
**Python 3 (ipykernel)**.

Step C2 — run it, top to bottom. Click the first cell, then press **Shift+Enter**
repeatedly (each press runs one cell and moves to the next), or use the menu
Run → Run All Cells. Expect the big simulation cell to run for a while (minutes, not
seconds — it is integrating millions of years of physics at 12-digit accuracy; the
cell prints progress lines as it goes).

Step C3 — at the end, the notebook bakes and opens the interactive display page
(`mercury_orbit.html`) in your browser. Press **Space** to play. Watch the ratio dial:
it starts near 181 (a young Mercury spinning once every 11.6 hours), falls for
millions of years, hesitates at the 3:2 line — and locks. Use the scrub bar to jump
to the capture moment; the "libration" plot below shows the resonance angle rocking
like a settling pendulum, the fingerprint of tidal locking.

Everything the notebook prints is also saved: the raw simulation output in
`planet_Mercury/data/runs/`, the database in
`planet_Mercury/data/mercury_orbit.sqlite3`, the display page in
`planet_Mercury/gui/mercury_orbit.html`.

Headless / classroom-batch note: to re-run the whole notebook non-interactively (for
grading, or to re-verify a machine), use the standard-library runner — from
`planet_Mercury/notebook/`:

```
MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_tidal_locking.ipynb
python3 check_notebook.py mercury_tidal_locking.ipynb
```

The first command executes every code cell in order and fails loudly if any cell
errors; the second audits the notebook's structure rules (Part E). Both must pass.
(These two commands use plain `python3` — no venv needed, since they use only the
standard library.)

---

## Part D — What each section of the notebook contains (the walk-through)

The notebook's cells alternate: a markdown (text) cell that explains, then a code
cell that does exactly what was just explained. The sections:

**§1 Title and abstract.** What the notebook computes and the three numbers it must
reproduce (87.969 d / 58.646 d / 175.938 d).

**§2 How to run this notebook.** The Shift+Enter instructions (repeated inside the
notebook so it is self-contained).

**§3 The words used in this notebook.** A glossary: tide, torque, resonance,
eccentricity, mean anomaly, libration, ODE, tolerance — each in one plain sentence.

**§4 The physical situation.** The full model, in the notebook's own words: the Sun
as a point mass; Mercury as an extended, deformable, slightly lopsided body; the five
tracked quantities [a, e, M, θ, Ω]; every governing equation written out (Kepler's
mean motion; the Hut tidal brake with its five eccentricity polynomials; the triaxial
"handle" torque; the two slow orbit-reshaping equations); every constant with its
value and meaning; and the two documented deviations from the source specification
(the 1000× time-compression of the tidal strength, and the two-stage integration far
from vs. near resonance) with the reasons why, spelled out honestly.

**§4.3 The first-order system actually handed to SUNDIALS.** The five equations
gathered as dy/dt = F(t, y), plus the solver contract: CVODE, BDF method, Newton
iteration, dense linear solver, relative tolerance 1.0×10⁻¹², absolute tolerances
[1.0×10⁻³, 1.0×10⁻⁶, 1.0×10⁻¹⁰, 1.0×10⁻¹⁰, 1.0×10⁻¹⁴], max step 864,000 s.

**§5 How this notebook talks to the simulator.** The driver code: a small Python
helper that locates the built `mercury_rs` binary (environment variable
`MERCURY_BIN`, else `../mercury_rs/target/release/mercury_rs`, else it tells you to
run the Part B build command), runs a subcommand with `subprocess.run`, streams the
program's printed progress into the cell output, and raises an error if the program
does not end with `SUCCESS`.

**§6 Building and running the simulation.** One explained code cell per run:
- 6.1 `print-config` — read the constants back from the program itself.
- 6.2 Run A, spec-literal — shows honestly that with the specification's own tidal
  strength the spin barely budges in 10 million years (the braking really takes
  billions of years — that IS the finding), motivating the compression.
- 6.3 Run B — the time-compressed braking movie: fast spin, the 2:1 pass, and the
  restart state saved when the spin ratio first reaches 1.6.
- 6.4 Run C, the 64-branch phase sweep — capture is a game of chance; the cell prints
  the measured capture fraction beside the ≈ 7% theoretical estimate (with error
  bars), and identifies the first captured branch.
- 6.5 Run B-final — continue that captured branch to the end of the window: Mercury
  locked, forever.
- 6.6 Run D, the encore — at eccentricity 0.285 capture is guaranteed; why ovalness
  is the secret ingredient.
- 6.7 Run E, the seam validation — the two integration stages agree where they meet.

**§7 Storing the results in a database.** Creates `mercury_orbit.sqlite3` with
Python's built-in `sqlite3` module, creates the five tables (`run`, `sample`,
`event`, `branch`, `target` — the full schema with every column explained appears in
the notebook text itself), and bulk-loads every CSV the runs produced. Prints row
counts.

**§8 Retrieving results from the database.** Seven documented SQL queries, each
explained then executed: (1) the final state of the canonical run; (2) the capture
time; (3) the libration amplitude over time (proving it decays); (4) the
orbital-period history; (5) the angular-momentum ledger (spin momentum lost ≈ orbit
momentum gained); (6) the sweep's capture fraction; (7) the traceability query — the
exact solver settings and step counts that produced the canonical history.

**§9 The verification gauntlet.** Assertions, all must pass: final Ω/n within 10⁻⁴
of 1.5; final P_orb / P_rot / P_solar within 0.1% of 87.969 / 58.646 / 175.938 Earth
days AND P_rot/P_orb = 2/3 to 5 significant figures; γ librating with decaying
amplitude; total angular momentum |L_tot(t) − L_tot(0)|/L_tot(0) below the documented
budget of 1×10⁻⁹; run A's verdict matches Finding F1 (the spec-literal tidal strength
cannot capture in 10 Myr). The cell prints one PASS line per check.

**§10 Baking and opening the display page.** Runs `gui/bake_page.py` on the
database, verifies the output byte-for-byte against a second bake, and opens it in
the browser (skipped when `MERCURY_NO_BROWSER=1`).

**§11 What we learned.** The story retold with the measured numbers filled in, plus
two "explore further" exercises (change the phase-sweep size; try other
eccentricities) and the explicit list of what was excluded (Jupiter, general
relativity) and why those belong to the second test.

**§12 Name and save this notebook.** An interactive final cell (standard project
convention) that asks for a name and opens a save-file dialog; batch runners skip
this cell automatically.

---

## Part E — For teachers: how the notebook is authored (and how to make your own)

These are the authoring rules the build follows, and that a class can follow to
create a sibling notebook (for example, for Venus or for the Moon):

1. **One notebook = one complete lesson.** Every explanation lives inside the
   notebook. Never write "see another file" — restate instead. (The engine project
   enforces the same rule on all 109 of its own notebooks.)
2. **Markdown before code, always.** Every code cell is preceded by a markdown cell
   of at least 80 characters that says what is about to happen and why.
3. **Plain Python 3, standard library only** in code cells: `subprocess` to run the
   Rust program, `sqlite3` for the database, `json`, `pathlib`, `csv`, `webbrowser`.
   No pip packages in cells (JupyterLab itself is only the viewer).
4. **The Rust program does all the physics.** Python never integrates, never
   "fixes up" numbers; it runs the program, loads its output, queries it, checks it,
   displays it. If you need new physics, extend `mercury_rs` (in Rust) — not the
   notebook.
5. **Real outputs are part of the notebook.** After authoring, execute top-to-bottom
   with `run_notebook.py` so the committed file contains genuine captured outputs;
   `check_notebook.py` then audits rules 1–3 mechanically (it checks: the how-to-run
   text exists; every code cell has its ≥80-character markdown lead-in; no
   cross-references to other notebooks; the required section headings exist; every
   non-interactive code cell was actually executed).
6. **Determinism is sacred.** Same inputs → byte-identical outputs. No timestamps, no
   random seeds without recording them, all floating-point printing through the fixed
   formatters. This is what makes "check and verify" meaningful.

---

## Part F — Troubleshooting (verified answers to the likely failures)

| Symptom | Cause and fix |
|---|---|
| `cargo: command not found` | Rust is not installed — install from https://rustup.rs, reopen the terminal |
| Running `python3` or `git` opens an "install command line developer tools" dialog | Normal on a fresh Mac — click Install, wait for it to finish, then re-run the command |
| `pip install` refuses with an "externally managed environment" message | You skipped the venv — use the two Part B4 commands exactly (they create a private environment and install into it) |
| The notebook's driver cell says it cannot find the binary | Run the Part B build (`cargo build --release` inside `planet_Mercury/mercury_rs/`), or set `MERCURY_BIN` to the binary's full path |
| Build error mentioning a missing `sundials_core`/`cvode_rs` path | The folder layout was changed — `planet_Mercury/` must sit beside (or inside) the engine repository exactly as shipped |
| A cell raises `RuntimeError: ... FAILURE` | The Rust program's self-check failed; read the cell's captured program output — it names the exact check. Re-run the cell once; if it persists, re-run `cargo test` and report what fails |
| The browser page does not open | Open `planet_Mercury/gui/mercury_orbit.html` yourself by double-clicking it; in batch runs set `MERCURY_NO_BROWSER=1` and no page is expected |
| The big run takes very long | Expected: minutes in release mode. If you accidentally built without `--release`, it can take hours — rebuild with `--release` |
| Numbers differ in the last digit from a friend's machine | Should not happen on the same platform (the engine's math is host-independent by design); re-run `run_notebook.py` and compare again — and if it truly persists, that IS a reportable finding |
