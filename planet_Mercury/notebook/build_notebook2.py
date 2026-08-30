#!/usr/bin/env python3
"""Author mercury_test2_jupiter_gr.ipynb deterministically.

Test 2 of the planet_Mercury project: the 3:2 spin-orbit capture WITH
Jupiter's secular forcing and Einstein's general-relativistic perihelion
correction. Builds the un-executed notebook JSON (nbformat 4, Python 3
ipykernel); run_notebook.py then executes it to embed the real captured
outputs. Standard library only. Re-running this builder always produces the
identical file.
"""

import json
from pathlib import Path


def md(text):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.strip("\n").splitlines(keepends=True),
    }


def code(text, interactive=False):
    cell = {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": text.strip("\n").splitlines(keepends=True),
    }
    if interactive:
        cell["metadata"]["tags"] = ["interactive"]
    return cell


CELLS = []

CELLS.append(md("""
# Mercury with Jupiter and Einstein — the 3:2 lock follows the precessing ellipse

## 1. What this notebook computes

This is **test 2** of the planet_Mercury project. Test 1 established, with a
pure two-body model, that the Sun's tides braked Mercury's spin until the
Sun's grip on Mercury's slightly lopsided shape snapped it into the strange
**3:2 spin-orbit resonance** — three spins for every two trips around the Sun.
Test 2 now adds the two pieces of physics test 1 deliberately excluded:

1. **Einstein's general-relativistic correction.** Near the Sun, gravity is
   very slightly stronger than Newton's law says. The orbit stays an ellipse,
   but the ellipse's long axis slowly swings around — the famous
   **43 arcseconds per century** of extra perihelion advance that Newtonian
   physics could not explain and Einstein's 1915 theory nailed exactly.
2. **Jupiter's perturbations.** Jupiter's steady pull, averaged over many
   orbits ("secular" forcing, in the classic Laplace-Lagrange form), does two
   things: it turns Mercury's ellipse an additional ~160 arcseconds per
   century, and it makes Mercury's eccentricity **breathe** — rising and
   falling by about 0.0045 in a roughly 800,000-year cycle.

The headline result this notebook must demonstrate: with the ellipse itself
turning, the locked spin no longer averages EXACTLY 1.5 orbits-per-orbit.
The lock tracks the **moving perihelion**, so the settled mean spin ratio is

    ratio = 3/2 + (perihelion turn rate) / (orbital rate)  =  1.5000004...

— a shift of about four parts in ten million, resolved cleanly by the solver
at its 12-digit accuracy. **The lock follows the precessing ellipse, not the
stars.** All computation is done by the same pure-Rust program (`mercury_rs`)
using the SUNDIALS 7.8.0 CVODE solver, the rustSolveIt engine's only
integration backend; this notebook's Python cells (standard library only) run
that program, store the results in a real SQLite database, query them back,
check every acceptance target, and bake a browser display page.

| Quantity | Value this notebook must reproduce |
|---|---|
| Einstein's extra perihelion advance | 42.98 arcsec/century |
| Jupiter's secular perihelion advance | ~160 arcsec/century (computed, then verified) |
| Jupiter's eccentricity cycle | amplitude ~0.0045, period ~809,000 years |
| Locked mean spin ratio | 1.5 + (perihelion rate)/(orbital rate), NOT 1.5 |
| Year / rotation / solar day today | 87.969 / 58.646 / 175.938 Earth days |

The notebook stops with an error if any check fails.
"""))

CELLS.append(md("""
## 2. How to run this notebook

One-time setup, in a terminal (5-10 minutes the first time):

1. Build the compute program — from the folder that contains `planet_Mercury/`:
   `cd planet_Mercury/mercury_rs` then `cargo build --release 2>&1 | tee /tmp/mercury_build.log`
   (zero errors, zero warnings expected; install Rust from https://rustup.rs if
   `cargo` is missing).
2. Install JupyterLab into a small private environment (modern macOS Python
   refuses system-wide installs): from the folder that contains `planet_Mercury/`:
   `python3 -m venv planet_Mercury/notebook/.venv` then
   `planet_Mercury/notebook/.venv/bin/pip install jupyterlab`.
3. Open this notebook: `cd planet_Mercury/notebook` then
   `.venv/bin/jupyter lab mercury_test2_jupiter_gr.ipynb`.
   If asked to pick a kernel, choose **Python 3 (ipykernel)**.

Then run the cells top to bottom: click the first cell and press
**Shift+Enter** repeatedly (each press runs one cell and moves on), or use
Run -> Run All Cells. The big simulation cells print progress as they go and
take minutes each — they are integrating millions of years of physics at
12-digit accuracy. Batch re-verification without a browser:
`MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_test2_jupiter_gr.ipynb`.
"""))

CELLS.append(md("""
## 3. The words used in this notebook

- **Tide**: the stretching a body feels because gravity pulls its near side
  harder than its far side. Flexing under that stretch wastes energy as heat.
- **Torque**: a twisting push that changes how fast something spins.
- **Eccentricity (e)**: how oval an orbit is; 0 is a perfect circle. Mercury's
  is 0.20563 — noticeably stretched.
- **Semi-major axis (a)**: the orbit's size — half its long diameter.
- **Perihelion**: the point of the orbit closest to the Sun. The angle of that
  point, measured around the Sun, is called **pomega** (written as the symbol
  curly-pi in textbooks) and is the sixth variable of this test.
- **Precession**: the slow turning of the ellipse's long axis — the perihelion
  swinging around the Sun like the slow hand of a clock.
- **Arcsecond**: 1/3600 of a degree. 43 arcseconds is roughly the width of a
  human hair seen from 30 meters away — per CENTURY. Tiny, and measurable.
- **Mean anomaly (M)**: a "clock angle" marking where a planet is along its
  orbit; it advances at a perfectly steady rate.
- **Mean motion (n)**: the orbit's average angular speed; the year is 2*pi/n.
- **Spin rate (Omega)**: how fast the planet rotates, in radians per second.
- **Resonance**: a lock between two rhythms — here, 3 spins per 2 orbits.
- **Libration**: the gentle rocking of a locked angle around its resting
  point, like a settling pendulum. A librating resonance angle IS the lock.
- **Secular**: astronomers' word for the slow, averaged-over-many-orbits part
  of a motion, as opposed to its fast wiggles. Jupiter's influence here is
  purely secular — the classic Laplace-Lagrange averaged form.
- **Laplace coefficient**: a standard number b(alpha) from celestial
  mechanics that measures how strongly two orbits of size ratio alpha talk to
  each other; the program computes it by numerical integration and this
  notebook's test suite verifies it against the textbook power series.
- **General relativity (GR)**: Einstein's 1915 theory of gravity. Its only
  effect on this problem is the extra 42.98 arcsec/century perihelion drift.
- **ODE / tolerance**: an ordinary differential equation is a rule for how
  fast something changes; the tolerance is the accuracy the solver must
  maintain at every single step (relative 1e-12 here).
"""))

CELLS.append(md("""
## 4. The physical situation

### 4.1 The stage, restated in full

Exactly two bodies are integrated. The **Sun** is a point mass. **Mercury**
is an extended, deformable, nearly spherical body: it flexes under the Sun's
tide (squishiness given by its Love number k2 = 0.12 and a constant time lag
tau), and it carries a tiny permanent "handle" — an equator longer in one
direction than the other by one part in ten thousand, (B-A)/C = 1e-4.

Two torques act on the spin. The **tidal torque** (Hut 1981,
constant-time-lag form) always drags the spin toward 1.256x the orbital rate
(the "pseudo-synchronous" rate for Mercury's oval orbit) and slowly feeds the
lost spin momentum into the orbit, changing a and e as it does. The **handle
torque** (Goldreich & Peale 1966) is the Sun's grip on the lopsided equator;
it averages to nearly zero while the spin is fast, but at a resonance
crossing it can snap the spin into step. Because the spin must fall THROUGH
3:2 on its way to 1.256, and because each crossing is a phase-dependent coin
flip, Mercury had its chance to be caught — and was.

As in test 1, the honest finding stands: with the specification's real tidal
numbers the braking takes billions of years, so the "movie" runs compress the
TIDAL strength by a documented factor of 1000 to fit the 10-million-year
integration window. Every compressed result says so.

### 4.2 The two new pushes, and what they change

**Einstein's push.** GR adds a steady perihelion drift
d(pomega)/dt = 3 n G M_sun / (c^2 a (1 - e^2)) — at Mercury's orbit,
42.98 arcsec/century. It changes no other variable.

**Jupiter's push (Laplace-Lagrange secular form).** Averaged over both
orbits, Jupiter contributes two coupled effects governed by two computed
rates A11 and A12 (built from Jupiter's mass, the orbit size ratio
alpha = a_Mercury / a_Jupiter, and Laplace coefficients):

- d(pomega)/dt gains A11 + A12 (e_J / e) cos(pomega - pomega_J): a further
  ~160 arcsec/century of turning, slightly modulated;
- de/dt gains A12 e_J sin(pomega - pomega_J): the eccentricity breathes with
  amplitude |A12/A11| * e_J (about 0.0045) as pomega circulates past
  Jupiter's own (held fixed at 0), with period 2*pi/A11 (about 809,000 yr).

**The crucial geometric consequence.** The resonance is a lock between the
spin and the ORBIT'S SHAPE — the resonance angle is now measured from the
moving perihelion: gamma2 = 2*theta - 3*M - 2*pomega. When gamma2 librates
(the lock), the spin's long-term average must satisfy
Omega = 1.5 n + d(pomega)/dt: the spin runs very slightly FASTER than
exactly 3:2, by precisely the perihelion turn rate. That offset — about
+3.8e-7 in the ratio — is the headline check of this whole test.

**What the movie compression does and does not touch.** Only the tidal
strength is compressed 1000x. Einstein's and Jupiter's rates are REAL,
uncompressed rates — so the braking movie shows about 12 Jupiter
eccentricity cycles in its 10-million-year window (the real, uncompressed
history would have run through thousands). The headline lock-offset check is
unaffected: it compares the measured mean spin ratio against the same real
perihelion rate that acted in the integration.

### 4.3 The first-order system actually handed to SUNDIALS

Six state variables y = [a, e, M, pomega, theta, Omega]:

1. da/dt — tidal back-reaction on the orbit size (Hut), 1000x in movie runs.
2. de/dt — tidal circularization (Hut, 1000x) **plus** Jupiter's
   A12 e_J sin(pomega - pomega_J).
3. dM/dt = n(a) — Kepler's clock ticks at the mean motion.
4. d(pomega)/dt = GR term + Jupiter's A11 + A12 (e_J/e) cos(pomega - pomega_J).
5. d(theta)/dt = Omega — the rotation angle integrates the spin rate.
6. d(Omega)/dt — tidal braking torque (1000x) plus, when the handle stage is
   active, the handle torque -(3/2) (G M_sun (B-A)/C / r^3) sin(2(theta - f - pomega)),
   where f is the true anomaly solved from Kepler's equation at every
   evaluation and the perihelion angle now appears because the handle's
   reference direction rides on the ellipse.

Solver: CVODE (BDF multistep + Newton iteration + dense linear solver) from
the pure-Rust SUNDIALS 7.8.0 port, relative tolerance 1e-12, per-variable
absolute tolerances [1e-3, 1e-6, 1e-10, 1e-10, 1e-10, 1e-14], maximum step
10 days. CVODE's exact root-finder stops runs precisely at chosen spin/orbit
ratios (the 2.2 stage handover, the 1.6 restart save, the 1.5 crossing). The
same two-stage strategy as test 1 is used (handle torque averaged out while
the spin is far above resonance, full model below 2.2), and the same angle
re-anchoring keeps M and theta small — a step that leaves gamma2 EXACTLY
unchanged, which the program's unit tests prove.
"""))

CELLS.append(md("""
## 5. How this notebook talks to the simulator

The next cell defines the driver: it locates the built `mercury_rs` binary
(environment variable `MERCURY_BIN` first, then
`../mercury_rs/target/release/mercury_rs`, else it tells you to build), runs
one subcommand at a time with `subprocess`, streams the program's printed
progress into the cell output, and raises an error unless the program ends
with its own `SUCCESS` verdict. All simulation output lands in
`planet_Mercury/data/runs/<run_id>/` as plain CSV (comma-separated values)
files plus a `manifest.json` configuration echo; test 2's sample files carry
one extra column, `pomega_rad`, the perihelion angle.
"""))

CELLS.append(code("""
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

NB_DIR = Path.cwd()
BASE = NB_DIR.parent if NB_DIR.name == "notebook" else NB_DIR
MR = BASE / "mercury_rs"
DATA = BASE / "data" / "runs"
DB_PATH = BASE / "data" / "mercury_test2.sqlite3"
YEAR = 3.15576e7

def find_binary():
    env = os.environ.get("MERCURY_BIN", "")
    if env and Path(env).exists():
        return Path(env)
    p = MR / "target" / "release" / "mercury_rs"
    if p.exists():
        return p
    raise RuntimeError(
        "mercury_rs binary not found - build it first: "
        "cd planet_Mercury/mercury_rs && cargo build --release")

def run(*args):
    binary = find_binary()
    proc = subprocess.Popen([str(binary), *args], cwd=str(MR),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    lines = []
    for line in proc.stdout:
        print(line, end="")
        lines.append(line.rstrip("\\n"))
    proc.wait()
    if proc.returncode != 0 or (lines and lines[-1] != "SUCCESS"):
        raise RuntimeError(f"mercury_rs {' '.join(args)} ended in FAILURE "
                           f"(exit {proc.returncode}) - read the output above")
    return lines

print("driver ready; simulator:", find_binary().name, "| data ->", DATA)
"""))

CELLS.append(md("""
## 6. Building and running the simulation

### 6.1 Read the configuration back from the program itself

Before trusting any result, we make the simulator recite its own shared
constants — masses, radius, squishiness, tolerances, the compression factor —
so what you verify below is what actually ran (never a number remembered by
hand). Test 2's OWN constants — Jupiter's mass, orbit and eccentricity, the
computed A11/A12 rates, and the GR rate — are echoed by each test-2 run into
its `manifest.json` and land in this notebook's database, where section 8.7
queries them back as provenance.
"""))

CELLS.append(code('run("print-config")\n'))

CELLS.append(md("""
### 6.2 Gate A — Einstein alone must give the famous 43

First, physics in isolation. With tides, handle, and Jupiter all switched
off, ONLY the GR term drives pomega. A thousand years of integration must
reproduce the analytic rate — 42.98 arcseconds per century, the number that
made Einstein's reputation — to better than one part in a thousand (the
program checks itself and prints its verdict).
"""))

CELLS.append(code('run("t2-gr-check")\n'))

CELLS.append(md("""
### 6.3 Gate B — Jupiter alone must breathe the eccentricity

Second isolation test: ONLY Jupiter's secular terms active. Two million
years of integration must show the eccentricity oscillating with the
Laplace-Lagrange amplitude |A12/A11|*e_J and period 2*pi/A11 — both to 5%.
This proves the coupled e-pomega secular system is integrated correctly
before it is mixed with tides and capture physics.
"""))

CELLS.append(code('run("t2-jupiter-check")\n'))

CELLS.append(md("""
### 6.4 The braking movie with all the physics on

Now everything at once: tides 1000x, GR, Jupiter, and (below spin ratio 2.2)
the handle torque. The spin falls from 181x the orbital rate through the 5:2
and 2:1 crossings while the perihelion steadily turns and the eccentricity
breathes; CVODE's root-finder stops exactly at ratio 1.6 and saves a
full-precision restart state — the common starting line for every sweep
branch below.
"""))

CELLS.append(code('run("t2-movie")\n'))

CELLS.append(md("""
### 6.5 The 16-branch phase sweep — is capture still a coin flip?

Whether a 3:2 crossing captures depends on the spin's phase angle on
arrival. Sixteen branches restart from the saved state, each with its spin
angle nudged by k*pi/16, and each runs to a definite outcome: CAPTURED into
the (now precessing) 3:2 lock, or passed through. Jupiter's breathing
eccentricity is in play the whole time, wiggling each branch's capture odds —
the point is that capture still happens robustly with the full physics on.
"""))

CELLS.append(code("""
run("t2-sweep", "--branches", "16")
rows = list(csv.DictReader(open(DATA / "T2_sweep" / "branches.csv", encoding="utf-8")))
captured = [int(r["branch_id"]) for r in rows if r["captured"] == "1"]
frac = len(captured) / len(rows)
se = math.sqrt(frac * (1.0 - frac) / len(rows))
print(f"measured capture fraction: {len(captured)}/{len(rows)} = "
      f"{frac:.3f} +/- {se:.3f} (68% confidence)")
assert captured, "zero captures: re-sweep with a finer phase grid"
CANONICAL = captured[0]
print("canonical branch (first captured):", CANONICAL)
"""))

CELLS.append(md("""
### 6.6 The canonical history — capture and lock with Jupiter and Einstein

The first captured branch is continued from the restart to the end of the
10-million-year window. Watch for the capture event, then the program's own
checks: the settled mean spin ratio equal to 1.5 PLUS the perihelion rate
over the orbital rate (the headline), the eccentricity still breathing on
Jupiter's cycle right through capture and lock, the total perihelion advance
matching prediction, and today's 87.969-day year. Note what is deliberately
absent: test 1's angular-momentum ledger. Jupiter's secular terms EXCHANGE
angular momentum with Jupiter itself, which this model does not track, so a
two-body ledger would be physically wrong here — the program prints the same
explanation.
"""))

CELLS.append(code('run("t2-final", "--branch", str(CANONICAL))\n'))

CELLS.append(md("""
## 7. Storing the results in a database

Everything the five runs produced goes into one SQLite database file —
`data/mercury_test2.sqlite3` — via Python's built-in `sqlite3` module.
Six tables (the documented schema):

- `run` — one row per run: the full configuration echo (masses, tolerances,
  initial conditions, solver identity, step counts) and the verdict.
- `run_extra` — key/value numbers specific to test 2, straight from each
  run's manifest: Jupiter's mass/orbit/eccentricity, the computed A11 and
  A12 secular rates, the GR rate, the starting perihelion angle, sweep
  branch bookkeeping. This is the notebook's provenance for every new
  constant.
- `sample` — every time-history row, keyed by (run_id, idx): time, orbit
  (a, e, M, pomega), spin (theta, Omega), derived columns (mean motion,
  spin/orbit ratio, the resonance angle gamma2 = 2*theta - 3*M - 2*pomega,
  periods, angular momenta, energies) and the integration stage tag.
- `event` — stage handovers, resonance crossings, the capture, restart saves.
- `branch` — each sweep branch's outcome (captured or passed, when, final
  ratio, and which one is canonical).
- `target` — the observed-Mercury numbers and Einstein's 42.98, with
  sources, so every later check compares against a stored, attributed value.

The cell prints row counts and proves referential integrity (no orphan
rows) before anything else is trusted.
"""))

CELLS.append(code("""
import sqlite3

SCHEMA = '''
CREATE TABLE run (
  run_id TEXT PRIMARY KEY, description TEXT NOT NULL, k2 REAL NOT NULL,
  tau_lag_s REAL NOT NULL, compression REAL NOT NULL, a0_m REAL NOT NULL,
  e0 REAL NOT NULL, M0_rad REAL NOT NULL, theta0_rad REAL NOT NULL,
  Omega0_rad_s REAL NOT NULL, t_final_s REAL NOT NULL, rel_tol REAL NOT NULL,
  abs_tol_a REAL NOT NULL, abs_tol_e REAL NOT NULL, abs_tol_M REAL NOT NULL,
  abs_tol_theta REAL NOT NULL, abs_tol_Omega REAL NOT NULL,
  max_step_s REAL NOT NULL, solver TEXT NOT NULL, n_steps INTEGER NOT NULL,
  n_rhs_evals INTEGER NOT NULL, n_reanchor INTEGER NOT NULL,
  verdict TEXT NOT NULL, engine TEXT NOT NULL);
CREATE TABLE run_extra (
  run_id TEXT NOT NULL REFERENCES run(run_id), key TEXT NOT NULL,
  value REAL NOT NULL, PRIMARY KEY (run_id, key));
CREATE TABLE sample (
  run_id TEXT NOT NULL REFERENCES run(run_id), idx INTEGER NOT NULL,
  t_s REAL NOT NULL, a_m REAL NOT NULL, e REAL NOT NULL, M_rad REAL NOT NULL,
  theta_rad REAL NOT NULL, Omega_rad_s REAL NOT NULL, n_rad_s REAL NOT NULL,
  ratio REAL NOT NULL, gamma2_rad REAL NOT NULL, P_orb_s REAL NOT NULL,
  P_rot_s REAL NOT NULL, L_spin_kgm2s REAL NOT NULL, L_orb_kgm2s REAL NOT NULL,
  L_tot_kgm2s REAL NOT NULL, E_spin_j REAL NOT NULL, E_orb_j REAL NOT NULL,
  stage TEXT NOT NULL, pomega_rad REAL NOT NULL, PRIMARY KEY (run_id, idx));
CREATE TABLE event (
  run_id TEXT NOT NULL REFERENCES run(run_id), t_s REAL NOT NULL,
  event TEXT NOT NULL, value REAL NOT NULL);
CREATE TABLE branch (
  run_id TEXT NOT NULL REFERENCES run(run_id), branch_id INTEGER NOT NULL,
  theta_offset_rad REAL NOT NULL, captured INTEGER NOT NULL,
  t_outcome_s REAL NOT NULL, final_ratio REAL NOT NULL,
  canonical INTEGER NOT NULL, PRIMARY KEY (run_id, branch_id));
CREATE TABLE target (
  name TEXT PRIMARY KEY, value REAL NOT NULL, source TEXT NOT NULL);
CREATE INDEX idx_sample_time ON sample(run_id, t_s);
'''

RUN_IDS = ["T2_gr_check", "T2_jupiter", "T2_movie", "T2_sweep", "T2_final"]
if DB_PATH.exists():
    DB_PATH.unlink()
con = sqlite3.connect(DB_PATH)
cur = con.cursor()
cur.executescript(SCHEMA)
for rid in RUN_IDS:
    d = DATA / rid
    man = json.loads((d / "manifest.json").read_text(encoding="utf-8"))
    cur.execute("INSERT INTO run VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
        man["run_id"], man["description"], man["k2"], man["tau_lag_s"],
        man["compression"], man["a0_m"], man["e0"], man["M0_rad"],
        man["theta0_rad"], man["Omega0_rad_s"], man["t_final_s"],
        man["rel_tol"], *man["abs_tol"], man["max_step_s"], man["solver"],
        man["n_steps"], man["n_rhs_evals"], man["n_reanchor"],
        man["verdict"], man["engine"]))
    cur.executemany("INSERT INTO run_extra VALUES (?,?,?)",
                    ((rid, k, float(v)) for k, v in man["extras"].items()))
    spath = d / "samples.csv"
    if spath.exists():
        with open(spath, encoding="utf-8") as f:
            r = csv.reader(f)
            next(r)
            cur.executemany(
                "INSERT INTO sample VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                ((rid, i, *[float(x) for x in row[:16]], row[16], float(row[17]))
                 for i, row in enumerate(r)))
    with open(d / "events.csv", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r)
        cur.executemany("INSERT INTO event VALUES (?,?,?,?)",
                        ((rid, float(row[0]), row[1], float(row[2])) for row in r))
    bpath = d / "branches.csv"
    if bpath.exists():
        with open(bpath, encoding="utf-8") as f:
            r = csv.reader(f)
            next(r)
            cur.executemany("INSERT INTO branch VALUES (?,?,?,?,?,?,?)",
                            ((rid, int(row[0]), float(row[1]), int(row[2]),
                              float(row[3]), float(row[4]), int(row[5])) for row in r))
cur.executemany("INSERT INTO target VALUES (?,?,?)", [
    ("P_orb_days", 87.969, "radar/ephemerides; Pettengill & Dyce 1965 lineage"),
    ("P_rot_days", 58.646, "radar/ephemerides; exactly 2/3 of the year"),
    ("P_solar_days", 175.938, "derived: exactly two Mercury years in a 3:2 lock"),
    ("ratio", 1.5, "the 3:2 spin-orbit resonance (test 2 measures the offset ABOVE this)"),
    ("gr_arcsec_cy", 42.98, "Einstein 1915; Le Verrier's unexplained remainder"),
])
con.commit()
for table in ["run", "run_extra", "sample", "event", "branch", "target"]:
    n = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    print(f"table {table:9s}: {n} rows")
orphans = cur.execute(
    "SELECT (SELECT COUNT(*) FROM sample WHERE run_id NOT IN (SELECT run_id FROM run))"
    " + (SELECT COUNT(*) FROM event  WHERE run_id NOT IN (SELECT run_id FROM run))"
    " + (SELECT COUNT(*) FROM branch WHERE run_id NOT IN (SELECT run_id FROM run))"
    " + (SELECT COUNT(*) FROM run_extra WHERE run_id NOT IN (SELECT run_id FROM run))"
).fetchone()[0]
print("referential integrity: every stored row points at a real run ->",
      "OK" if orphans == 0 else f"{orphans} ORPHANS")
assert orphans == 0
"""))

CELLS.append(md("""
## 8. Retrieving results from the database

Storage is only half the assignment — now we ask the database questions.
Seven documented SQL queries follow (unit conversions used throughout:
3.15576e13 seconds = 1 million years; 86400 seconds = 1 day; multiplying a
rate in rad/s by 6.50925...e14 gives arcseconds per century).

### 8.1 Query 1 — the present day (the final state of the canonical run)

The last sample of run T2_final is "today": the spin/orbit ratio near 1.5,
the year and rotation matching the observed Mercury, and a perihelion angle
that has turned through many full circles since the run began.
"""))

CELLS.append(code("""
ARCSEC_CY = 206264.80624709636 * 3.15576e9
q1 = cur.execute('''
SELECT t_s/3.15576e13, ratio, P_orb_s/86400.0, P_rot_s/86400.0, pomega_rad
FROM sample WHERE run_id='T2_final' ORDER BY idx DESC LIMIT 1''').fetchone()
print(f"t = {q1[0]:.3f} Myr | ratio = {q1[1]:.8f} | "
      f"P_orb = {q1[2]:.4f} d | P_rot = {q1[3]:.4f} d")
print(f"perihelion angle = {q1[4]:.2f} rad = "
      f"{q1[4]/(2*math.pi):.2f} full turns of the ellipse since t = 0")
FINAL = {"ratio": q1[1], "p_orb_d": q1[2], "p_rot_d": q1[3]}
"""))

CELLS.append(md("""
### 8.2 Query 2 — when was Mercury captured?

The capture event pins the moment the resonance angle gamma2 (measured from
the MOVING perihelion) stopped circulating and began librating. (Remember
the movie clock: the tides ran 1000x strong, so multiply the braking
timescale by ~1000 for the real history.)
"""))

CELLS.append(code("""
q2 = cur.execute('''
SELECT t_s/3.15576e13, value FROM event
WHERE run_id='T2_final' AND event='capture_detected'
ORDER BY t_s LIMIT 1''').fetchone()
print(f"capture at t = {q2[0]:.4f} Myr of movie time (ratio locked near {q2[1]})")
T_CAPTURE_S = q2[0] * 3.15576e13
"""))

CELLS.append(md("""
### 8.3 Query 3 — the headline: the lock follows the precessing ellipse

The settled-era mean spin ratio (every sample from one million years after
capture onward, averaged in Python over a FIXED iteration order so the
printed digits cannot depend on the database engine's summation plan) is
compared against 1.5 + (GR rate + Jupiter's A11)/n. The two must agree to
about a part in ten million — and BOTH must sit measurably ABOVE exactly
1.5. The perihelion rates come from the run's own manifest echo
(`run_extra`), not from numbers typed here.
"""))

CELLS.append(code("""
xt = dict(cur.execute(
    "SELECT key, value FROM run_extra WHERE run_id='T2_final'").fetchall())
G = 6.67430e-11; M_SUN = 1.98847e30; M_MERC = 3.3011e23
a0 = cur.execute("SELECT a0_m FROM run WHERE run_id='T2_final'").fetchone()[0]
n0 = math.sqrt(G * (M_SUN + M_MERC) / a0**3)
PW_DOT_PRED = xt["gr_rate_rad_s"] + xt["ll_A11_rad_s"]
EXPECTED_RATIO = 1.5 + PW_DOT_PRED / n0
rows_r = cur.execute(
    "SELECT ratio FROM sample WHERE run_id='T2_final' AND t_s > ? ORDER BY idx",
    (T_CAPTURE_S + 3.15576e13,)).fetchall()
MEAN_RATIO = sum(r[0] for r in rows_r) / len(rows_r)
print(f"settled-era mean spin ratio: {MEAN_RATIO:.10f}  ({len(rows_r)} samples)")
print(f"predicted 1.5 + pomega_dot/n: {EXPECTED_RATIO:.10f}")
print(f"offset above exactly 1.5:    {MEAN_RATIO - 1.5:.3e} "
      f"(prediction {EXPECTED_RATIO - 1.5:.3e})")
print("the lock follows the precessing ellipse, not the stars")
"""))

CELLS.append(md("""
### 8.4 Query 4 — Jupiter's breathing eccentricity, through thick and thin

Mercury's eccentricity must keep oscillating on Jupiter's cycle right
through braking, capture, and lock. The swing observed in the canonical run
is compared with the Laplace-Lagrange amplitude |A12/A11|*e_J; the cycle
period is measured from the clean Jupiter-only gate run by locating
successive eccentricity maxima.
"""))

CELLS.append(code("""
q4 = cur.execute('''SELECT MIN(e), MAX(e) FROM sample
WHERE run_id='T2_final' ''').fetchone()
AMP_MEASURED = 0.5 * (q4[1] - q4[0])
AMP_PRED = abs(xt["ll_A12_rad_s"] / xt["ll_A11_rad_s"]) * xt["e_jupiter"]
print(f"e range in the canonical run: {q4[0]:.5f} .. {q4[1]:.5f} "
      f"-> swing amplitude {AMP_MEASURED:.3e} (LL prediction {AMP_PRED:.3e})")
ej = cur.execute('''SELECT t_s, e FROM sample WHERE run_id='T2_jupiter'
ORDER BY idx''').fetchall()
maxima = [ej[i][0] for i in range(1, len(ej) - 1)
          if ej[i][1] > ej[i-1][1] and ej[i][1] > ej[i+1][1]]
PERIOD_MEASURED = (maxima[-1] - maxima[0]) / (len(maxima) - 1)
PERIOD_PRED = 2 * math.pi / xt["ll_A11_rad_s"]
print(f"cycle period (Jupiter-only gate): {PERIOD_MEASURED/3.15576e10:.1f} kyr "
      f"(LL prediction {PERIOD_PRED/3.15576e10:.1f} kyr)")
"""))

CELLS.append(md("""
### 8.5 Query 5 — the perihelion odometer

The total distance the perihelion traveled in the canonical run, divided by
the elapsed time, must equal the predicted GR + Jupiter rate to 5% — with
Einstein's 42.98 arcseconds per century of it separately attributed. (Test
1's angular-momentum ledger is deliberately absent here: Jupiter's secular
terms exchange angular momentum with Jupiter itself, which this model does
not track, so a two-body ledger would be checking the wrong conservation
law. The program's own output says the same.)
"""))

CELLS.append(code("""
q5a = cur.execute('''SELECT t_s, pomega_rad FROM sample
WHERE run_id='T2_final' ORDER BY idx LIMIT 1''').fetchone()
q5b = cur.execute('''SELECT t_s, pomega_rad FROM sample
WHERE run_id='T2_final' ORDER BY idx DESC LIMIT 1''').fetchone()
PW_RATE = (q5b[1] - q5a[1]) / (q5b[0] - q5a[0])
print(f"perihelion advanced {q5b[1] - q5a[1]:.2f} rad over "
      f"{(q5b[0]-q5a[0])/3.15576e13:.3f} Myr")
print(f"measured rate  : {PW_RATE*ARCSEC_CY:8.2f} arcsec/century")
print(f"predicted rate : {PW_DOT_PRED*ARCSEC_CY:8.2f} arcsec/century "
      f"(GR {xt['gr_rate_rad_s']*ARCSEC_CY:.2f} + Jupiter {xt['ll_A11_rad_s']*ARCSEC_CY:.2f})")
"""))

CELLS.append(md("""
### 8.6 Query 6 — the measured capture odds with the full physics on

The sweep's record, straight from the branch table: how many of the 16
phase branches locked while Jupiter breathed the eccentricity and Einstein
turned the ellipse. Sixteen coin flips give rough odds — the robust point
is that capture happens at all under the full physics.
"""))

CELLS.append(code("""
q6 = cur.execute('''
SELECT AVG(captured), COUNT(*), SUM(captured) FROM branch
WHERE run_id='T2_sweep' ''').fetchone()
print(f"capture fraction = {q6[2]}/{q6[1]} = {q6[0]:.4f} "
      f"(test 1's two-body sweep measured ~0.16; Goldreich-Peale theory ~0.07)")
CAPTURE_FRACTION = q6[0]
"""))

CELLS.append(md("""
### 8.7 Query 7 — traceability: what produced the canonical history?

Every number above traces back to a run row plus its `run_extra` block
recording the exact solver settings, step counts, Jupiter/GR constants, and
verdict that produced it — the database is its own provenance record.
"""))

CELLS.append(code("""
q7 = cur.execute('''
SELECT solver, rel_tol, abs_tol_Omega, max_step_s, n_steps, n_reanchor, verdict
FROM run WHERE run_id='T2_final' ''').fetchone()
print(f"solver {q7[0]} | rel_tol {q7[1]:.0e} | abs_tol_Omega {q7[2]:.0e}")
print(f"max_step {q7[3]:.0f} s | {q7[4]} CVODE steps | "
      f"{q7[5]} re-anchorings | verdict {q7[6]}")
print("test-2 constants echoed by the run itself:")
for k in ["m_jupiter_kg", "a_jupiter_m", "e_jupiter",
          "ll_A11_rad_s", "ll_A12_rad_s", "gr_rate_rad_s"]:
    print(f"  {k:14s} = {xt[k]:.6e}")
"""))

CELLS.append(md("""
## 9. The verification gauntlet

Every acceptance target, asserted in one place. If any line fails, the
notebook stops here with an error — there is no "mostly passed."
"""))

CELLS.append(code("""
def gauntlet(name, ok, detail):
    print(("PASS - " if ok else "FAIL - ") + name + ": " + detail)
    assert ok, name

verdicts = cur.execute("SELECT COUNT(*) FROM run WHERE verdict='SUCCESS'").fetchone()[0]
gauntlet("all_run_verdicts", verdicts == 5, f"{verdicts}/5 runs ended SUCCESS")
gauntlet("capture_detected", T_CAPTURE_S > 0,
         f"capture at {T_CAPTURE_S/3.15576e13:.4f} Myr with Jupiter + GR active")
gr_as = xt["gr_rate_rad_s"] * ARCSEC_CY
tgt_gr = cur.execute("SELECT value FROM target WHERE name='gr_arcsec_cy'").fetchone()[0]
gauntlet("gr_is_einsteins_number", abs(gr_as - tgt_gr) < 0.05,
         f"the run's own GR rate {gr_as:.3f} arcsec/cy vs stored target {tgt_gr}")
gauntlet("lock_follows_the_precessing_ellipse",
         abs(MEAN_RATIO - EXPECTED_RATIO) <= 1.5e-7,
         f"settled mean ratio {MEAN_RATIO:.10f} vs 1.5 + pomega_dot/n = "
         f"{EXPECTED_RATIO:.10f} (|diff| = {abs(MEAN_RATIO-EXPECTED_RATIO):.2e})")
gauntlet("shift_off_exact_three_halves", (MEAN_RATIO - 1.5) >= 2e-7,
         f"mean ratio sits {MEAN_RATIO-1.5:.3e} ABOVE exactly 1.5 - the lock "
         f"tracks the moving perihelion, not the fixed stars")
gauntlet("P_orb", abs(FINAL["p_orb_d"] - 87.969) / 87.969 < 1e-3,
         f"{FINAL['p_orb_d']:.4f} d vs observed 87.969 d")
gauntlet("P_rot", abs(FINAL["p_rot_d"] - 58.646) / 58.646 < 1e-3,
         f"{FINAL['p_rot_d']:.4f} d vs observed 58.646 d")
p_solar = 1.0 / abs(1.0 / FINAL["p_rot_d"] - 1.0 / FINAL["p_orb_d"])
gauntlet("P_solar", abs(p_solar - 175.938) / 175.938 < 1e-3,
         f"{p_solar:.4f} d vs observed 175.938 d (two Mercury years)")
gauntlet("jupiter_cycle_amplitude", abs(AMP_MEASURED - AMP_PRED) / AMP_PRED < 0.3,
         f"e swing {AMP_MEASURED:.3e} vs Laplace-Lagrange {AMP_PRED:.3e} "
         f"(persists through braking, capture, and lock)")
gauntlet("jupiter_cycle_period",
         abs(PERIOD_MEASURED - PERIOD_PRED) / PERIOD_PRED < 0.05,
         f"{PERIOD_MEASURED/3.15576e10:.1f} kyr vs 2*pi/A11 = "
         f"{PERIOD_PRED/3.15576e10:.1f} kyr")
gauntlet("perihelion_advance", abs(PW_RATE - PW_DOT_PRED) / PW_DOT_PRED < 0.05,
         f"measured {PW_RATE*ARCSEC_CY:.2f} vs predicted "
         f"{PW_DOT_PRED*ARCSEC_CY:.2f} arcsec/century")
gauntlet("sweep_measured_odds", 0.0 < CAPTURE_FRACTION < 1.0,
         f"fraction {CAPTURE_FRACTION:.4f} (16-flip error bars apply)")
print("ALL CHECKS PASSED")
"""))

CELLS.append(md("""
## 10. Baking and opening the display page

The last computational step turns the database into one self-contained
browser page (`gui/mercury_test2.html`): the spin-ratio descent into the
lock, the breathing eccentricity, and the perihelion odometer climbing
through sixteen-odd full turns, with the headline numbers on top. The cell
bakes the page TWICE and proves the two bakes are byte-identical
(determinism), then opens it in your browser — set MERCURY_NO_BROWSER=1 to
skip opening (batch runs do).
"""))

CELLS.append(code("""
import hashlib
import webbrowser

page = BASE / "gui" / "mercury_test2.html"
for i in (1, 2):
    r = subprocess.run([sys.executable, str(BASE / "gui" / "bake_page2.py")],
                       capture_output=True, text=True, cwd=str(BASE))
    print(r.stdout, end="")
    if r.returncode != 0:
        print(r.stderr, end="")
        raise RuntimeError("bake_page2.py failed")
    if i == 1:
        h1 = hashlib.sha256(page.read_bytes()).hexdigest()
h2 = hashlib.sha256(page.read_bytes()).hexdigest()
assert h1 == h2, "two bakes differed - determinism broken"
print("bake determinism: byte-identical (sha256 " + h1[:16] + "...)")
if os.environ.get("MERCURY_NO_BROWSER"):
    print("MERCURY_NO_BROWSER set - not opening a browser; the page is at:", page)
else:
    webbrowser.open(page.as_uri())
    print("opened in your browser:", page)
"""))

CELLS.append(md("""
## 11. What we learned

- **Einstein's 43 arcseconds are in the machine.** With everything else
  switched off, the integrated perihelion advanced at 42.98 arcsec/century —
  matching the analytic general-relativistic rate to a part in a thousand,
  the same number that historically separated Newton from Einstein.
- **Jupiter breathes Mercury's orbit.** The secular (averaged) pull makes
  the eccentricity oscillate by ~0.0045 on a ~809,000-year cycle and turns
  the ellipse ~160 arcsec/century — both reproduced against the
  Laplace-Lagrange pencil-and-paper values, and both persisting right
  through braking, capture, and lock.
- **The lock follows the precessing ellipse, not the stars.** The resonance
  angle is measured from the MOVING perihelion, so the settled spin averages
  1.5 + pomega_dot/n — about four parts in ten million ABOVE exactly 3:2.
  The measured offset matched the prediction at the part-in-ten-million
  level: a genuinely relativistic-plus-planetary fingerprint on Mercury's
  clock, resolved by a 12-digit integration.
- **Capture survives the full physics.** With the eccentricity breathing and
  the ellipse turning, the 16-branch sweep still locked a healthy fraction
  of its phase branches — Mercury's trap works in the real, messy solar
  system, not just in a clean two-body toy.
- **Honest bookkeeping.** Test 1's angular-momentum ledger is deliberately
  absent: Jupiter's secular terms exchange angular momentum with Jupiter,
  which this model does not track. Checking a two-body conservation law here
  would be checking the wrong thing — so the notebook says so instead.

**Explore further (exercises):** (1) rerun section 6.5 with more branches
for tighter odds; (2) set Jupiter's eccentricity to zero in the simulator
source and watch the breathing stop while the 160 arcsec/century advance
survives; (3) double Jupiter's mass and predict — before running — how the
lock offset and the cycle period must change.
"""))

CELLS.append(md("""
## 12. Name and save this notebook

Run the final cell to save a named copy of this executed notebook wherever
you like (a save-file dialog opens; if the dialog cannot open, you will be
asked to type a folder path instead). Batch runners skip this cell.
"""))

CELLS.append(code("""
name = input("Name for this notebook copy [mercury_test2_jupiter_gr]: ").strip() \\
       or "mercury_test2_jupiter_gr"
dest = ""
try:
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk()
    root.withdraw()
    dest = filedialog.asksaveasfilename(defaultextension=".ipynb",
                                        initialfile=name + ".ipynb")
    root.destroy()
except Exception:
    print("Falling back to a typed folder path")
    folder = input("Folder to save into: ").strip()
    if folder:
        dest = str(Path(folder) / (name + ".ipynb"))
if dest:
    import shutil
    shutil.copyfile("mercury_test2_jupiter_gr.ipynb", dest)
    print("saved:", dest)
else:
    print("not saved")
""", interactive=True))


def main():
    # Deterministic cell ids (nbformat 4.5 requires ids; random ones would
    # make Jupyter-saved and builder-built files diverge for no reason).
    for i, c in enumerate(CELLS):
        c["id"] = f"cell-{i:03d}"
    nb = {
        "cells": CELLS,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3 (ipykernel)",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
            "mercury": {"pairs_with": "mercury_rs/src/test2.rs"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out = Path(__file__).resolve().parent / "mercury_test2_jupiter_gr.ipynb"
    out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
    print(f"built {out} ({len(CELLS)} cells)")


if __name__ == "__main__":
    main()
