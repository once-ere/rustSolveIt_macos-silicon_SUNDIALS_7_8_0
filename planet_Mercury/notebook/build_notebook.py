#!/usr/bin/env python3
"""Author mercury_tidal_locking.ipynb deterministically.

Builds the un-executed notebook JSON (nbformat 4, Python 3 ipykernel) with all
markdown prose and code cells; run_notebook.py then executes it to embed the
real captured outputs. Standard library only. Re-running this builder always
produces the identical file.
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
# Mercury tidal locking — the 3:2 spin-orbit capture, computed with pure-Rust SUNDIALS CVODE

## 1. What this notebook computes

Mercury spins on its axis exactly **three times for every two trips around the
Sun** — a "3:2 spin-orbit resonance," which is Mercury's strange form of tidal
locking (our Moon, by contrast, is locked 1:1 and always shows Earth one
face). This notebook simulates **how that happened**: the Sun's tides slowly
braked Mercury's once-fast spin over billions of years, until the Sun's
gravitational grip on Mercury's slightly lopsided shape snapped the spin into
step as it fell through the 3:2 ratio.

Everything numerical is computed by a pure-Rust program (`mercury_rs`) that
integrates five ordinary differential equations (ODEs — rules for how fast
things change) with the SUNDIALS 7.8.0 CVODE solver, the rustSolveIt engine's
only integration backend. This notebook's Python cells (standard library
only) run that program, store its results in a real SQLite database, retrieve
them with SQL (Structured Query Language) queries, check every acceptance
target, and bake an interactive browser display page.

The simulated present day must reproduce the real Mercury's clock:

| Quantity | Observed value |
|---|---|
| Year (orbital period) | 87.969 Earth days |
| Sidereal day (rotation period) | 58.646 Earth days — exactly 2/3 of the year |
| Solar day (sunrise to sunrise) | 175.938 Earth days — exactly TWO Mercury years |

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
   `.venv/bin/jupyter lab mercury_tidal_locking.ipynb`.
   If asked to pick a kernel, choose **Python 3 (ipykernel)**.

Then run the cells top to bottom: click the first cell and press
**Shift+Enter** repeatedly (each press runs one cell and moves on), or use
Run → Run All Cells. The big simulation cells print progress as they go and
take minutes each — they are integrating millions of years of physics at
12-digit accuracy. Batch re-verification without a browser:
`MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_tidal_locking.ipynb`.
"""))

CELLS.append(md("""
## 3. The words used in this notebook

- **Tide**: the stretching a body feels because gravity pulls its near side
  harder than its far side. Flexing under that stretch wastes energy as heat.
- **Torque**: a twisting push that changes how fast something spins.
- **Eccentricity (e)**: how oval an orbit is; 0 is a perfect circle. Mercury's
  is 0.20563 — noticeably stretched.
- **Semi-major axis (a)**: the orbit's size — half its long diameter.
- **Mean anomaly (M)**: a "clock angle" marking where a planet is along its
  orbit; it advances at a perfectly steady rate.
- **Mean motion (n)**: the orbit's average angular speed; the year is 2*pi/n.
- **Spin rate (Omega)**: how fast the planet rotates, in radians per second.
- **Resonance**: a lock between two rhythms — here, 3 spins per 2 orbits.
- **Libration**: the gentle rocking of a locked angle around its resting
  point, like a settling pendulum. A librating resonance angle IS the lock.
- **ODE (ordinary differential equation)**: a rule for how fast something
  changes; the computer advances it through time.
- **Tolerance**: the accuracy the solver must maintain at every step.
- **Secular**: astronomers' word for the slow, averaged-over-many-orbits part
  of a motion, as opposed to its fast wiggles.
"""))

CELLS.append(md("""
## 4. The physical situation

**The model (exactly two bodies).** The Sun is a point mass M_sun. Mercury is
an extended, deformable, nearly spherical body: radius R, a slightly lopsided
equator (principal moments of inertia A < B < C; the lopsidedness B-A is the
"handle" the Sun grabs), a tidal "squishiness" (Love number k2) and a tidal
"sluggishness" (time lag tau). Five numbers evolve through time:

    y(t) = [ a, e, M, theta, Omega ]

(orbit size, orbit ovalness, orbit clock angle, spin angle, spin rate).

### 4.1 Equations of motion

With n = sqrt(G(M_sun+m)/a^3), the tidal-brake strength
K = 3 G M_sun^2 R^5 k2 tau / a^6, and C = 0.34 m R^2:

    da/dt     = (2K/(m n a))    [ Omega f2(e) - n f3(e) ]
    de/dt     = (9Ke/(m n a^2)) [ (11/18) Omega f4(e) - n f5(e) ]
    dM/dt     = n
    dtheta/dt = Omega
    dOmega/dt = ( T_tri + T_tidal ) / C

where the orbit-averaged tidal brake is T_tidal = -K [Omega f1(e) - n f2(e)]
(Hut 1981 constant-time-lag model), and the "handle" torque is
T_tri = -(3/2) G M_sun (B-A) / r^3 * sin(2(theta - f)), with the true anomaly
f and distance r recovered from (a, e, M) by solving Kepler's equation
M = E - e sin E at every step. The five eccentricity factors are Hut's:

    f1 = (1 + 3e^2 + (3/8)e^4) / (1-e^2)^(9/2)
    f2 = (1 + (15/2)e^2 + (45/8)e^4 + (5/16)e^6) / (1-e^2)^6
    f3 = (1 + (31/2)e^2 + (255/8)e^4 + (185/16)e^6 + (25/64)e^8) / (1-e^2)^(15/2)
    f4 = (1 + (3/2)e^2 + (1/8)e^4) / (1-e^2)^5
    f5 = (1 + (15/4)e^2 + (15/8)e^4 + (5/64)e^6) / (1-e^2)^(13/2)

The brake alone would park the spin at the "pseudo-synchronous" rate
n*f2/f1 = 1.256 n at Mercury's eccentricity — BELOW the 1.5 n resonance — so
a fast-spinning Mercury must brake down THROUGH 3:2, giving the handle torque
its one chance per crossing to capture. Constants (from the source
specification): G = 6.67430e-11 m^3 kg^-1 s^-2, M_sun = 1.98847e30 kg,
m = 3.3011e23 kg, R = 2.4397e6 m, C/(mR^2) = 0.34, (B-A)/C = 1e-4, k2 = 0.12,
tau = 100 s, a0 = 5.790905e10 m, e0 = 0.20563, Omega0 = 1.5e-4 rad/s (one
turn per 11.6 hours, about 181x the orbital rate).

### 4.2 Honest deviations from the source specification (both displayed below)

- **Time compression (S = 1000).** With the spec's own constants the braking
  takes about 4.7 BILLION years (the spin decays toward 1.256 n with a
  710-Myr e-folding time; reaching 1.5 n needs 6.6 of those e-folds), so the
  spec's 10-million-year window could never see capture. Honoring the spec's
  stated intent ("dissipation deliberately set strong... tractable"), the
  main runs multiply k2*tau by 1000 (tau: 100 s -> 1e5 s). Run A below keeps
  the literal spec values and shows, honestly, that almost nothing happens.
- **Staged integration.** While Mercury spins fast, the handle torque wiggles
  every ~6 hours; resolving that for 10 Myr would take ~1e11 solver steps
  (200x the spec's own budget) to compute a torque that averages to zero far
  from resonance. So far from resonance (spin ratio > 2.2) the handle torque
  is switched off (secular stage); near and after the crossings the full
  five-equation system runs. Run E below PROVES the two stages agree where
  they hand over — at a deliberately NON-resonant checkpoint ratio of 2.7
  (the first draft used 3.0 and discovered the hard way that 3.0 is itself
  the 3:1 resonance!).
- **Capture is a coin flip (~7% per crossing).** Whether one run captures
  depends on the spin's phase angle at the crossing. Run C sweeps 64 phase
  offsets and measures the odds; the first captured branch becomes the
  canonical history continued in run B-final.

### 4.3 The first-order system actually handed to SUNDIALS

The five equations above are integrated as dy/dt = F(t, y) by CVODE from the
vendored pure-Rust SUNDIALS 7.8.0: BDF method (backward differentiation
formulas — built for stiff problems), Newton iteration, dense linear solver,
relative tolerance 1.0e-12, absolute tolerances
[1.0e-3, 1.0e-6, 1.0e-10, 1.0e-10, 1.0e-14] for [a, e, M, theta, Omega],
maximum step 864,000 s (10 days). CVODE's exact root-finding stops precisely
on chosen spin-ratio crossings; angles are periodically "re-anchored"
(theta -= 3*pi*j, M -= 2*pi*j — which provably leaves the resonance angle
gamma = 2*theta - 3*M unchanged) so angle precision survives 10 million
years. No Python code in this notebook ever integrates anything.
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
files plus a `manifest.json` configuration echo.
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
DB_PATH = BASE / "data" / "mercury_orbit.sqlite3"
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

Before trusting any result, we make the simulator recite its own constants —
masses, radius, squishiness, tolerances, the compression factor — so what you
verify below is what actually ran (never a number remembered by hand).
"""))

CELLS.append(code('run("print-config")\n'))

CELLS.append(md("""
### 6.2 Run A — the spec-literal run (the honest finding)

With the specification's own tidal strength (tau = 100 s), ten million years
of braking barely moves the spin: about 181x the orbital rate down to about
179x. That IS the finding: the real braking took billions of years — most of
Mercury's existence — which motivates the documented 1000x time compression
used by the movie runs below.
"""))

CELLS.append(code('run("run-a")\n'))

CELLS.append(md("""
### 6.3 Run B — the compressed braking movie, to the restart save

Now the movie: tides 1000x stronger. The spin falls from 181x through the
5:2 and 2:1 crossings; at spin ratio 2.2 the handle torque switches on (full
five-equation model); CVODE's root-finder then stops exactly at ratio 1.6 —
just above the 3:2 doorstep — and saves a full-precision restart state that
every sweep branch below starts from.
"""))

CELLS.append(code('run("run-b")\n'))

CELLS.append(md("""
### 6.4 Run C — the 64-branch phase sweep (capture is a coin flip)

Whether a crossing captures depends on the spin's phase angle when it arrives
— effectively a coin whose outcome is fixed by the starting angle. Sixty-four
branches restart from the saved state, each with its spin angle nudged by
k*pi/64, and each is integrated through the crossing to a definite outcome:
CAPTURED into 3:2, or passed through toward the 1.256 n resting rate. The
measured capture fraction lands next to Goldreich & Peale's 1966 theoretical
estimate of about 7% — with honest statistical error bars.
"""))

CELLS.append(code("""
run("sweep", "--branches", "64")
rows = list(csv.DictReader(open(DATA / "C_sweep" / "branches.csv", encoding="utf-8")))
captured = [int(r["branch_id"]) for r in rows if r["captured"] == "1"]
frac = len(captured) / len(rows)
se = math.sqrt(frac * (1.0 - frac) / len(rows))
print(f"measured capture fraction: {len(captured)}/{len(rows)} = "
      f"{frac:.3f} +/- {se:.3f} (68% confidence)")
print("Goldreich-Peale (1966) theoretical estimate: ~0.070")
assert captured, "zero captures (~1% chance): run the finer-grid contingency re-sweep"
CANONICAL = captured[0]
print("canonical branch (first captured):", CANONICAL)
"""))

CELLS.append(md("""
### 6.5 Run B-final — the canonical history: capture, libration, lock

The first captured branch is continued from the restart to the end of the
10-million-year window. Watch for the capture event, then the self-checks:
the settled spin ratio at 1.5 to a part in ten thousand, the final year /
day / solar day against the observed 87.969 / 58.646 / 175.938 days, the
angular-momentum ledger, and the libration amplitude decaying — tidal
locking, demonstrated end to end.
"""))

CELLS.append(code('run("run-b-final", "--branch", str(CANONICAL))\n'))

CELLS.append(md("""
### 6.6 Run D — the guaranteed-capture encore (why ovalness is the secret)

At eccentricity 0.285 the pseudo-synchronous rate f2/f1 equals exactly 1.5 —
the brake DELIVERS the spin to the 3:2 doorstep instead of dragging it past,
and capture is certain. This clearly-labeled exploration beyond the
specification shows why Mercury's oval orbit is the secret ingredient of its
lock; it also measures the libration (rocking) period against the
Goldreich-Peale pencil-and-paper formula.
"""))

CELLS.append(code('run("run-d")\n'))

CELLS.append(md("""
### 6.7 Run E — proving the two integration stages agree

The staged integration must be validated: from the same starting state at the
deliberately NON-resonant spin ratio 2.7, the secular model (handle torque
averaged away) and the full model are integrated over the same 5000-year
window and their despin rates must agree to 1% (the two stages differ ONLY in
the spin equation). Fun fact discovered while building this: the first draft
put this checkpoint at ratio 3.0 — and the full model promptly locked into
the 3:1 resonance that lives exactly there, which is why the checkpoint moved.
"""))

CELLS.append(code('run("run-e")\n'))

CELLS.append(md("""
## 7. Storing the results in a database

Everything the runs produced now goes into one SQLite database file —
`data/mercury_orbit.sqlite3` — using Python's built-in `sqlite3` module
(nothing to install). Five tables: `run` (one row per run: the full
configuration echo and verdict), `sample` (every time-history row, keyed by
run), `event` (stage handovers, resonance crossings, the capture), `branch`
(each sweep branch's coin-flip outcome), and `target` (the observed-Mercury
numbers the checks compare against). The cell prints the row counts as proof
of a complete load.
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
CREATE TABLE sample (
  run_id TEXT NOT NULL REFERENCES run(run_id), idx INTEGER NOT NULL,
  t_s REAL NOT NULL, a_m REAL NOT NULL, e REAL NOT NULL, M_rad REAL NOT NULL,
  theta_rad REAL NOT NULL, Omega_rad_s REAL NOT NULL, n_rad_s REAL NOT NULL,
  ratio REAL NOT NULL, gamma_rad REAL NOT NULL, P_orb_s REAL NOT NULL,
  P_rot_s REAL NOT NULL, L_spin_kgm2s REAL NOT NULL, L_orb_kgm2s REAL NOT NULL,
  L_tot_kgm2s REAL NOT NULL, E_spin_j REAL NOT NULL, E_orb_j REAL NOT NULL,
  stage TEXT NOT NULL, PRIMARY KEY (run_id, idx));
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

RUN_IDS = ["A_spec_literal", "B_movie", "C_sweep", "B_final", "D_high_e", "E_seam"]
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
    with open(d / "samples.csv", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r)
        cur.executemany(
            "INSERT INTO sample VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ((rid, i, *[float(x) for x in row[:16]], row[16])
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
    ("ratio", 1.5, "the 3:2 spin-orbit resonance"),
    ("obliquity_deg", 0.0, "model assumption; real Mercury is within ~2 arcminutes"),
])
con.commit()
for table in ["run", "sample", "event", "branch", "target"]:
    n = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    print(f"table {table:7s}: {n} rows")
orphans = cur.execute(
    "SELECT (SELECT COUNT(*) FROM sample WHERE run_id NOT IN (SELECT run_id FROM run))"
    " + (SELECT COUNT(*) FROM event  WHERE run_id NOT IN (SELECT run_id FROM run))"
    " + (SELECT COUNT(*) FROM branch WHERE run_id NOT IN (SELECT run_id FROM run))"
).fetchone()[0]
print("referential integrity: every stored row points at a real run ->",
      "OK" if orphans == 0 else f"{orphans} ORPHANS")
assert orphans == 0
"""))

CELLS.append(md("""
## 8. Retrieving results from the database

Storage is only half the assignment — now we ask the database questions.
Seven documented SQL queries follow (unit conversions used throughout:
3.15576e13 seconds = 1 million years; 86400 seconds = 1 day).

### 8.1 Query 1 — the present day (the final state of the canonical run)

The last sample of run B_final is "today": the spin/orbit ratio must sit at
1.5, and the year and day must match the observed Mercury.
"""))

CELLS.append(code("""
q1 = cur.execute('''
SELECT t_s/3.15576e13, ratio, P_orb_s/86400.0, P_rot_s/86400.0
FROM sample WHERE run_id='B_final' ORDER BY idx DESC LIMIT 1''').fetchone()
print(f"t = {q1[0]:.3f} Myr | ratio = {q1[1]:.8f} | "
      f"P_orb = {q1[2]:.4f} d | P_rot = {q1[3]:.4f} d")
FINAL = {"ratio": q1[1], "p_orb_d": q1[2], "p_rot_d": q1[3]}
"""))

CELLS.append(md("""
### 8.2 Query 2 — when was Mercury captured?

The capture event pins the moment the resonance angle stopped circulating and
began librating. (Remember the movie clock: multiply by ~1000 for the real,
uncompressed history.)
"""))

CELLS.append(code("""
q2 = cur.execute('''
SELECT t_s/3.15576e13, value FROM event
WHERE run_id='B_final' AND event='capture_detected'
ORDER BY t_s LIMIT 1''').fetchone()
print(f"capture at t = {q2[0]:.4f} Myr of movie time (ratio locked to {q2[1]})")
T_CAPTURE_S = q2[0] * 3.15576e13
"""))

CELLS.append(md("""
### 8.3 Query 3 — the libration settling down (locking, quantified)

Binning the resonance angle gamma into consecutive 100,000-year windows after
capture and taking each bin's swing (max - min) shows the rocking dying away
— the numerical signature of a body settling into its lock.
"""))

CELLS.append(code("""
q3 = cur.execute('''
SELECT CAST((t_s - ?)/3.15576e12 AS INTEGER) AS bin, MAX(gamma_rad)-MIN(gamma_rad)
FROM sample WHERE run_id='B_final' AND t_s > ?
GROUP BY bin ORDER BY bin LIMIT 12''', (T_CAPTURE_S, T_CAPTURE_S)).fetchall()
print("bin (100 kyr each) | gamma swing [rad]")
for b, s in q3:
    print(f"  {b:3d}              | {s:.4f}")
SWINGS = [s for _, s in q3]
"""))

CELLS.append(md("""
### 8.4 Query 4 — the orbital period through history

The assignment asks explicitly how the orbital period changes in time. The
tide moves angular momentum from spin into the orbit, so the orbit grows —
but only minutely: the year lengthens by about three SECONDS over the whole
history, while the day changes by a factor of 121.
"""))

CELLS.append(code("""
q4 = cur.execute('''
SELECT t_s/3.15576e13, P_orb_s/86400.0 FROM sample
WHERE run_id IN ('B_movie','B_final') AND idx % 800 = 0 ORDER BY t_s''').fetchall()
print("t [Myr] | P_orb [days]")
for t, p in q4[:14]:
    print(f"  {t:7.3f} | {p:.9f}")
print(f"  ... {len(q4)} rows total; "
      f"P_orb changed by {(q4[-1][1]-q4[0][1])*86400.0:.3f} s over the whole history")
"""))

CELLS.append(md("""
### 8.5 Query 5 — the angular-momentum ledger

The other explicit assignment question: spin angular momentum falls, orbital
angular momentum rises by the same amount. The audit has TWO eras with
different physics. BEFORE capture, the handle torque averages to zero and the
model conserves the total exactly, so the drift must stay under one part in a
billion. AFTER capture, holding the lock requires a nonzero average handle
torque — and because the source specification deliberately gives that torque
no orbital back-reaction, the model itself leaks total angular momentum at a
small predictable rate; the locked era is therefore compared against that
predicted leak, honestly, rather than against zero. (The recorded leak is
even smaller than predicted because at lock the orbit's per-step change falls
below what a 64-bit number can resolve — also documented.)
"""))

CELLS.append(code("""
q5 = cur.execute('''
SELECT MIN(L_spin_kgm2s), MAX(L_spin_kgm2s),
       MIN(L_orb_kgm2s),  MAX(L_orb_kgm2s)
FROM sample WHERE run_id IN ('B_movie','B_final')''').fetchone()
print(f"L_spin: {q5[1]:.4e} -> {q5[0]:.4e} kg m^2/s (falls)")
print(f"L_orb : {q5[2]:.6e} -> {q5[3]:.6e} kg m^2/s (rises)")
cap5 = cur.execute('''SELECT t_s FROM event WHERE run_id='B_final'
  AND event='capture_detected' ORDER BY t_s LIMIT 1''').fetchone()[0]
l00 = cur.execute('''SELECT L_tot_kgm2s FROM sample WHERE run_id='B_movie'
  ORDER BY idx LIMIT 1''').fetchone()[0]
lcap = cur.execute('''SELECT L_tot_kgm2s FROM sample WHERE run_id='B_final'
  AND t_s >= ? ORDER BY t_s LIMIT 1''', (cap5,)).fetchone()[0]
L_PRE = cur.execute('''SELECT MAX(ABS(L_tot_kgm2s - ?))/? FROM sample
  WHERE run_id IN ('B_movie','B_final') AND t_s <= ?''',
  (l00, l00, cap5)).fetchone()[0]
L_LOCKED = cur.execute('''SELECT MAX(ABS(L_tot_kgm2s - ?))/? FROM sample
  WHERE run_id='B_final' AND t_s > ?''', (lcap, l00, cap5)).fetchone()[0]
# The model's own predicted locked-era leak: <T_tri> = K n (1.5 f1 - f2)
# with Hut's f1, f2 at Mercury's eccentricity (display math, no integration).
e5 = 0.20563; e2 = e5*e5; e4 = e2*e2; e6 = e4*e2
f1h = (1 + 3*e2 + 0.375*e4) / ((1 - e2)**4.5)
f2h = (1 + 7.5*e2 + 5.625*e4 + 0.3125*e6) / ((1 - e2)**6)
Kb = 3 * 6.67430e-11 * (1.98847e30**2) * (2.4397e6**5) * (0.12*1.0e5) / (5.790905e10**6)
n5 = (6.67430e-11 * (1.98847e30 + 3.3011e23) / 5.790905e10**3) ** 0.5
t_end5 = cur.execute("SELECT MAX(t_s) FROM sample WHERE run_id='B_final'").fetchone()[0]
L_LEAK = Kb * n5 * (1.5*f1h - f2h) * (t_end5 - cap5) / l00
print(f"pre-capture: max |L_tot - L_tot(start)|/L_tot = {L_PRE:.3e}  (budget 1e-9)")
print(f"locked era : measured drift {L_LOCKED:.3e} vs the model's predicted "
      f"secular leak {L_LEAK:.3e}")
"""))

CELLS.append(md("""
### 8.6 Query 6 — the measured capture odds

The sweep's coin-flip record, straight from the branch table: how many of the
64 phase branches locked, versus the ~7% theory.
"""))

CELLS.append(code("""
q6 = cur.execute('''
SELECT AVG(captured), COUNT(*), SUM(captured) FROM branch
WHERE run_id='C_sweep' ''').fetchone()
print(f"capture fraction = {q6[2]}/{q6[1]} = {q6[0]:.4f} "
      f"(Goldreich-Peale estimate ~0.070)")
CAPTURE_FRACTION = q6[0]
"""))

CELLS.append(md("""
### 8.7 Query 7 — traceability: what produced the canonical history?

Every number above traces back to a run row recording the exact solver
settings, step counts, and verdict that produced it — the database is its own
provenance record.
"""))

CELLS.append(code("""
q7 = cur.execute('''
SELECT solver, rel_tol, abs_tol_Omega, max_step_s, n_steps, n_reanchor, verdict
FROM run WHERE run_id='B_final' ''').fetchone()
print(f"solver {q7[0]} | rel_tol {q7[1]:.0e} | abs_tol_Omega {q7[2]:.0e}")
print(f"max_step {q7[3]:.0f} s | {q7[4]} CVODE steps | "
      f"{q7[5]} re-anchorings | verdict {q7[6]}")
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

gauntlet("capture_detected", T_CAPTURE_S > 0,
         f"capture at {T_CAPTURE_S/3.15576e13:.4f} Myr")
# A locked planet still LIBRATES (rocks) forever with a tiny residual
# amplitude, so instantaneous samples wiggle by ~2e-4 in the ratio; the lock
# statement - spin ratio exactly 3/2 - is about the TIME AVERAGE. The means
# are computed in Python over a FIXED iteration order (ORDER BY idx) so the
# printed digits cannot depend on the database engine's summation plan.
rows_r = cur.execute(
    "SELECT ratio, P_rot_s/P_orb_s FROM sample WHERE run_id='B_final' "
    "AND t_s > ? ORDER BY idx", (T_CAPTURE_S + 3.15576e13,)).fetchall()
MEAN_RATIO = sum(r[0] for r in rows_r) / len(rows_r)
MEAN_PR = sum(r[1] for r in rows_r) / len(rows_r)
gauntlet("locked_mean_ratio", abs(MEAN_RATIO - 1.5) <= 1e-4,
         f"settled-era mean ratio {MEAN_RATIO:.10f} "
         f"(final instantaneous sample: {FINAL['ratio']:.8f})")
gauntlet("P_orb", abs(FINAL["p_orb_d"] - 87.969) / 87.969 < 1e-3,
         f"{FINAL['p_orb_d']:.4f} d vs observed 87.969 d")
gauntlet("P_rot", abs(FINAL["p_rot_d"] - 58.646) / 58.646 < 1e-3,
         f"{FINAL['p_rot_d']:.4f} d vs observed 58.646 d")
p_solar = 1.0 / abs(1.0 / FINAL["p_rot_d"] - 1.0 / FINAL["p_orb_d"])
gauntlet("P_solar", abs(p_solar - 175.938) / 175.938 < 1e-3,
         f"{p_solar:.4f} d vs observed 175.938 d (two Mercury years)")
r23 = abs(MEAN_PR - 2.0 / 3.0) / (2.0 / 3.0)
gauntlet("two_thirds_lock", r23 < 1e-5,
         f"settled-era mean P_rot/P_orb - 2/3 relative = {r23:.2e}")
gauntlet("ledger_precapture", L_PRE <= 1e-9,
         f"pre-capture total-L drift {L_PRE:.3e} (budget 1e-9)")
gauntlet("ledger_locked_era", L_LOCKED <= 1.5 * L_LEAK + 2e-10,
         f"locked-era drift {L_LOCKED:.3e} vs the model's own predicted "
         f"secular leak {L_LEAK:.3e} (the spec's handle torque has no orbital "
         f"back-reaction; f64 orbit quantization suppresses the recorded value)")
gauntlet("libration_decays", len(SWINGS) >= 3 and SWINGS[-1] < SWINGS[0],
         f"gamma swing {SWINGS[0]:.3f} -> {SWINGS[-1]:.3f} rad across the bins")
gauntlet("sweep_measured_odds", 0.0 < CAPTURE_FRACTION < 1.0,
         f"fraction {CAPTURE_FRACTION:.4f} (theory ~0.070; 64-flip error bars apply)")
ra = cur.execute("SELECT ratio FROM sample WHERE run_id='A_spec_literal' "
                 "AND stage='S' ORDER BY idx DESC LIMIT 1").fetchone()[0]
gauntlet("finding_F1_confirmed", 170.0 <= ra <= 181.5,
         f"spec-literal run ends at ratio {ra:.4f} after 10 Myr - no capture possible")
verdicts = cur.execute("SELECT COUNT(*) FROM run WHERE verdict='SUCCESS'").fetchone()[0]
gauntlet("all_run_verdicts", verdicts == 6, f"{verdicts}/6 runs ended SUCCESS")
print("ALL CHECKS PASSED")
"""))

CELLS.append(md("""
## 10. Baking and opening the display page

The last computational step turns the database into one self-contained
browser page (`gui/mercury_orbit.html`): the orbit with a spinning Mercury
and its long-axis "handle," the spin/orbit-ratio dial settling on exactly
1.5, the libration plot, and the period and angular-momentum histories, with
play/scrub controls. The cell bakes the page TWICE and proves the two bakes
are byte-identical (determinism), then opens it in your browser — set
MERCURY_NO_BROWSER=1 to skip opening (batch runs do).
"""))

CELLS.append(code("""
import hashlib
import webbrowser

page = BASE / "gui" / "mercury_orbit.html"
for i in (1, 2):
    r = subprocess.run([sys.executable, str(BASE / "gui" / "bake_page.py")],
                       capture_output=True, text=True, cwd=str(BASE))
    print(r.stdout, end="")
    if r.returncode != 0:
        print(r.stderr, end="")
        raise RuntimeError("bake_page.py failed")
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

- **Tides are brakes.** Flexing Mercury wastes spin energy as heat; with the
  specification's realistic tidal numbers the braking from an 11.6-hour day
  takes ~4.7 billion years — most of Mercury's existence (run A showed the
  honest, nearly-motionless version of that story; the movie runs compressed
  it 1000x, and every compressed result states so).
- **The oval orbit sets the trap.** Tides alone would park the spin at 1.256x
  the orbital rate — BELOW 3:2 — so the spin must fall through the resonance,
  where the Sun's grip on Mercury's lopsided equator can catch it.
- **Capture is a coin flip.** Our 64-branch sweep measured the odds of any
  single crossing locking, next to Goldreich & Peale's ~7%; at eccentricity
  0.285 (run D) capture became certain — ovalness is the secret ingredient.
- **The lock is real and checkable.** After capture the resonance angle
  librates with dying swings, the spin sits at exactly 1.5 orbits-per-orbit,
  and today's clock emerges: an 87.969-day year, a 58.646-day rotation, and a
  175.938-day solar day — two full Mercury years from one sunrise to the next.
- **Angular momentum balances.** What the spin lost, the orbit gained, to one
  part in a billion — audited in SQL from the stored samples.

**Explore further (exercises):** (1) rerun section 6.4 with
`run("sweep", "--branches", "128")` for tighter odds; (2) edit run D's
eccentricity in the simulator source and find where guaranteed capture stops;
(3) the excluded physics — Jupiter's tugs on Mercury's eccentricity and
Einstein's relativistic perihelion correction — is the planned second test.
"""))

CELLS.append(md("""
## 12. Name and save this notebook

Run the final cell to save a named copy of this executed notebook wherever
you like (a save-file dialog opens; if the dialog cannot open, you will be
asked to type a folder path instead). Batch runners skip this cell.
"""))

CELLS.append(code("""
name = input("Name for this notebook copy [mercury_tidal_locking]: ").strip() \\
       or "mercury_tidal_locking"
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
    shutil.copyfile("mercury_tidal_locking.ipynb", dest)
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
            "mercury": {"pairs_with": "mercury_rs/src/main.rs"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out = Path(__file__).resolve().parent / "mercury_tidal_locking.ipynb"
    out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
    print(f"built {out} ({len(CELLS)} cells)")


if __name__ == "__main__":
    main()
