# Planet Mercury Tidal-Locking Plan — Document 4 of 8: The Database Schema, Completely Documented

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science, including a
reader who has never used a database. Everything needed is inside this document.
**Status:** PLAN (awaiting approval). The schema below is what the build will create.

Words expanded once for this document: **CSV** = comma-separated values (a data table
saved as plain text); **SQL** = Structured Query Language — the standard mini-language
for asking a database questions; every query below can be pasted and run as-is;
**SI** = the international metric system (meters, kilograms, seconds). Symbols used in
column formulas: G is Newton's gravitational constant, Msun the Sun's mass, m
Mercury's mass, a the orbit size, e the orbit ovalness, n the orbit's average angular
speed, Omega the spin rate, theta the spin angle, M the orbit clock angle, and C
Mercury's spin moment of inertia (its resistance to changes in spin). The "Love
number" k2 measures how squishy Mercury is; "obliquity" is the tilt of the spin axis
(held at 0 in this model).

---

## 1. What the database is, and why we use one

The simulation program (`mercury_rs`, pure Rust) writes its results as plain CSV text
files — simple, portable, and byte-reproducible. The Jupyter notebook then loads all
of those files into **one SQLite database file**:

```
planet_Mercury/data/mercury_orbit.sqlite3
```

SQLite is a complete database engine that lives in a single ordinary file and is built
into Python (`import sqlite3` — nothing to install). A database gives us the two
things CSV files alone cannot:

- **Retrieval by question, not by file:** "What was the orbital period 2 million years
  in?" or "Which sweep branches got captured?" become one-line SQL queries.
- **Cross-run bookkeeping:** every stored row is traceable to the run that produced
  it, and every run row records the exact configuration that produced it — so a
  number can always be traced back to the settings that made it.

The notebook demonstrates both directions, storing AND retrieving, as the assignment
requires.

## 2. Conventions used in every table

- **Units are SI and are baked into column names:** `_s` seconds, `_m` meters,
  `_rad` radians, `_rad_s` radians per second, `_kgm2s` kilogram-meter²-per-second
  (angular momentum), `_j` joules. Time `t_s` is **simulated** seconds since the run's
  start (t = 0 is the start of the simulated history; there are no calendar dates
  anywhere in the data).
- **REAL** columns are 64-bit floating-point numbers; **TEXT** are strings;
  **INTEGER** are whole numbers. SQLite stores them exactly as given.
- Simulation-output tables (`sample`, `event`, `branch`) carry a `run_id` column
  tracing every row to its run; `run_id` matches the folder name of the raw CSVs
  (`planet_Mercury/data/runs/<run_id>/`). The one deliberate exception is `target`,
  which holds *observed real-Mercury* values, not simulation output, so it carries no
  `run_id`.
- The database is rebuilt from the CSVs by the notebook at any time; the CSVs remain
  the raw source of truth. Rebuilding twice from the same CSVs yields identical
  content (a verification gate). The CSV column headers match the table columns
  below one-for-one, by name.

## 3. The complete schema (the exact SQL the notebook executes)

```sql
-- Table 1: run — one row per simulation run (the configuration echo + verdict)
CREATE TABLE run (
    run_id            TEXT PRIMARY KEY,  -- 'A_spec_literal' | 'B_movie' | 'C_sweep'
                                         -- | 'B_final' | 'D_high_e' | 'E_seam'
    description       TEXT NOT NULL,     -- one plain-English sentence
    k2                REAL NOT NULL,     -- Love number used (0.12)
    tau_lag_s         REAL NOT NULL,     -- tidal time lag used [s] (100.0 or 1.0e5)
    compression       REAL NOT NULL,     -- time-compression factor S (1 or 1000)
    a0_m              REAL NOT NULL,     -- starting semi-major axis [m]
    e0                REAL NOT NULL,     -- starting eccentricity
    M0_rad            REAL NOT NULL,     -- starting mean anomaly [rad]
    theta0_rad        REAL NOT NULL,     -- starting spin angle [rad]
    Omega0_rad_s      REAL NOT NULL,     -- starting spin rate [rad/s]
    t_final_s         REAL NOT NULL,     -- requested end time [s]
    rel_tol           REAL NOT NULL,     -- CVODE relative tolerance (1.0e-12)
    abs_tol_a         REAL NOT NULL,     -- CVODE absolute tolerances, one per state:
    abs_tol_e         REAL NOT NULL,     --   [a, e, M, theta, Omega] =
    abs_tol_M         REAL NOT NULL,     --   [1e-3, 1e-6, 1e-10, 1e-10, 1e-14]
    abs_tol_theta     REAL NOT NULL,
    abs_tol_Omega     REAL NOT NULL,
    max_step_s        REAL NOT NULL,     -- CVODE maximum step (864000 s = 10 days)
    solver            TEXT NOT NULL,     -- 'CVODE_BDF_NEWTON_DENSE' (always)
    n_steps           INTEGER NOT NULL,  -- CVODE internal steps taken (for C_sweep:
                                         -- summed over all 64 branch integrations)
    n_rhs_evals       INTEGER NOT NULL,  -- right-hand-side evaluations (same rule)
    n_reanchor        INTEGER NOT NULL,  -- angle re-anchoring events (see Table 2 note)
    verdict           TEXT NOT NULL,     -- 'SUCCESS' or 'FAILURE' from the program
    engine            TEXT NOT NULL      -- 'sundials_rs 7.8.0 (pure Rust, macOS arm64)'
);

-- Table 2: sample — the time history (one row per output sample of a run).
-- Sweep branches do NOT write samples; the sweep's outcome lives in Table 4, and the
-- one branch that becomes the canonical history is re-run in full as run 'B_final'.
CREATE TABLE sample (
    run_id        TEXT    NOT NULL REFERENCES run(run_id),
    idx           INTEGER NOT NULL,      -- 0, 1, 2, ... within the run
    t_s           REAL NOT NULL,         -- simulated time [s]
    a_m           REAL NOT NULL,         -- semi-major axis (orbit size) [m]
    e             REAL NOT NULL,         -- eccentricity (orbit ovalness)
    M_rad         REAL NOT NULL,         -- mean anomaly, re-anchored [rad]
    theta_rad     REAL NOT NULL,         -- spin angle, re-anchored [rad]
    Omega_rad_s   REAL NOT NULL,         -- spin rate [rad/s]
    n_rad_s       REAL NOT NULL,         -- mean motion n = sqrt(G(Msun+m)/a^3) [rad/s]
    ratio         REAL NOT NULL,         -- Omega/n  (the star of the show; 1.5 = locked)
    gamma_rad     REAL NOT NULL,         -- resonance angle gamma = 2*theta - 3*M,
                                         -- reduced to (-pi, pi]; librates after capture
    P_orb_s       REAL NOT NULL,         -- orbital period 2*pi/n [s]
    P_rot_s       REAL NOT NULL,         -- rotation period 2*pi/Omega [s]
    L_spin_kgm2s  REAL NOT NULL,         -- spin angular momentum C*Omega
    L_orb_kgm2s   REAL NOT NULL,         -- orbital angular momentum m*n*a^2*sqrt(1-e^2)
    L_tot_kgm2s   REAL NOT NULL,         -- L_spin + L_orb (the conservation audit)
    E_spin_j      REAL NOT NULL,         -- spin kinetic energy (1/2)*C*Omega^2
    E_orb_j       REAL NOT NULL,         -- orbital energy -G*Msun*m/(2a)
    stage         TEXT NOT NULL,         -- 'S' secular stage | 'R' resonant-full stage
    PRIMARY KEY (run_id, idx)
);
-- Note on "re-anchored": theta and M are angles that would grow to ~10^10 radians
-- over a run; to preserve numerical precision the program periodically subtracts
-- 3*pi*j from theta and 2*pi*j from M (j = the whole number of completed 2*pi turns
-- in M), which leaves gamma = 2*theta - 3*M mathematically unchanged. Stored values
-- are the re-anchored (small) angles.

-- Table 3: event — notable moments detected during a run
CREATE TABLE event (
    run_id   TEXT NOT NULL REFERENCES run(run_id),
    t_s      REAL NOT NULL,              -- when it happened [simulated s]
    event    TEXT NOT NULL,              -- 'stage_handover' | 'cross_5:2' | 'cross_2:1'
                                         -- | 'cross_3:2' | 'capture_detected'
                                         -- | 'reanchor' | 'restart_saved'
    value    REAL NOT NULL               -- the ratio Omega/n at that moment (for
                                         -- 'reanchor': the running re-anchor count)
);

-- Table 4: branch — the phase-sweep bookkeeping (one row per branch per sweep run)
CREATE TABLE branch (
    run_id           TEXT NOT NULL REFERENCES run(run_id),  -- 'C_sweep' (a finer-grid
                                          -- contingency re-sweep, if one is ever
                                          -- needed, loads as 'C_sweep2')
    branch_id        INTEGER NOT NULL,    -- 0..63 within its sweep
    theta_offset_rad REAL NOT NULL,       -- the spin-phase offset added at the restart
                                          -- point (first sweep: branch_id * pi/64 —
                                          -- the handle torque repeats every pi in
                                          -- theta, so 64 such offsets tile one full
                                          -- period; a re-sweep records its own grid)
    captured         INTEGER NOT NULL,    -- 1 = locked into 3:2, 0 = sailed through
    t_outcome_s      REAL NOT NULL,       -- when the outcome was decided [s]
    final_ratio      REAL NOT NULL,       -- Omega/n at the branch's end
    canonical        INTEGER NOT NULL,    -- 1 = this branch was continued as 'B_final'
    PRIMARY KEY (run_id, branch_id)
);

-- Table 5: target — the observed-Mercury numbers the simulation must reproduce
-- (observed values, not simulation output — deliberately no run_id; see Section 2)
CREATE TABLE target (
    name     TEXT PRIMARY KEY,   -- 'P_orb_days' | 'P_rot_days' | 'P_solar_days'
                                 -- | 'ratio' | 'obliquity_deg'
    value    REAL NOT NULL,      -- 87.969 | 58.646 | 175.938 | 1.5 | 0.0
    source   TEXT NOT NULL       -- one-sentence citation, e.g. 'radar, Pettengill &
                                 -- Dyce 1965; modern ephemerides'
);

CREATE INDEX idx_sample_time ON sample(run_id, t_s);
```

## 4. Key columns explained once more, in plain words

- **run.compression** — how much faster than reality the tidal friction was set, so a
  multi-billion-year braking story fits in a 10-million-year simulation. `1` means
  "exactly the source specification's strength"; `1000` means "1000× stronger, the
  documented movie speed-up."
- **sample.ratio** — the single most important number in the project: how fast Mercury
  spins compared with its orbital rate. It starts near 181 (a spin every 11.6 hours),
  and Mercury is tidally locked when it stays at exactly 1.5 forever.
- **sample.gamma_rad** — the "lock dial." Before capture, γ spins round and round
  (circulates). After capture, it rocks gently back and forth around a fixed value
  (librates) with shrinking swings. A librating γ IS what "locked in 3:2 resonance"
  means, made visible.
- **sample.L_spin_kgm2s / L_orb_kgm2s / L_tot_kgm2s** — the angular-momentum ledger.
  The tide moves angular momentum from Mercury's spin into the orbit; L_tot must stay
  flat to within the documented budget of 1×10⁻⁹ (relative to its starting value) —
  the leftover covers solver tolerance plus one deliberately unmodeled, tiny,
  oscillatory term (the handle torque's back-reaction on the orbit), stated honestly.
- **sample.stage** — which model produced the row: `S` = secular (handle torque
  averaged away, far from resonance), `R` = full five-equation model (near and after
  the resonance crossings).
- **branch.theta_offset_rad** — capture is a game of chance decided by the spin's
  phase at the crossing; each branch nudges that phase by one tick and records the
  coin-flip outcome. **branch.canonical** marks the one branch whose continuation
  became the locked history.

## 5. Worked example queries (each appears, explained, in the notebook)

```sql
-- Unit conversions used below: 3.15576e13 s = 1 million years (Myr);
-- 3.15576e12 s = 100,000 years; 86400 s = 1 day; 'idx % 100 = 0' keeps every
-- 100th sample (% is the remainder operator).

-- Q1: The present day — the final state of the canonical locked run.
SELECT t_s/3.15576e13 AS t_Myr, ratio, P_orb_s/86400.0 AS P_orb_days,
       P_rot_s/86400.0 AS P_rot_days
FROM sample WHERE run_id = 'B_final' ORDER BY idx DESC LIMIT 1;
-- expect: ratio within 1e-4 of 1.5; periods within 0.1% of 87.969 / 58.646 days

-- Q2: When was Mercury captured?
SELECT t_s/3.15576e13 AS t_Myr, value AS ratio_at_capture
FROM event WHERE run_id = 'B_final' AND event = 'capture_detected';

-- Q3: The libration settling down (tidal locking, quantified):
--     the swing of gamma in consecutive 100,000-year bins after capture.
SELECT CAST(t_s/3.15576e12 AS INTEGER) AS bin_100kyr,
       MAX(gamma_rad) - MIN(gamma_rad) AS swing_rad
FROM sample WHERE run_id = 'B_final'
  AND t_s > (SELECT t_s FROM event WHERE run_id='B_final' AND event='capture_detected')
GROUP BY bin_100kyr ORDER BY bin_100kyr;
-- expect: swing_rad shrinking bin over bin

-- Q4: The orbital period through history (the assignment's explicit question).
SELECT t_s/3.15576e13 AS t_Myr, P_orb_s/86400.0 AS P_orb_days
FROM sample WHERE run_id = 'B_final' AND idx % 100 = 0 ORDER BY idx;

-- Q5: The angular-momentum ledger (the assignment's other explicit question).
SELECT t_s/3.15576e13 AS t_Myr, L_spin_kgm2s, L_orb_kgm2s, L_tot_kgm2s
FROM sample WHERE run_id = 'B_final' AND idx % 100 = 0 ORDER BY idx;
-- expect: L_spin falls, L_orb rises, L_tot flat to within the 1e-9 budget

-- Q6: The measured capture odds from the 64-branch sweep.
SELECT AVG(captured) AS capture_fraction, COUNT(*) AS branches
FROM branch WHERE run_id = 'C_sweep';
-- expect: a fraction in the neighborhood of the ~7% theoretical estimate
-- (binomial statistics: 64 flips at 7% typically give 2-8 captures)

-- Q7: Traceability — which solver settings produced the canonical history?
SELECT solver, rel_tol, abs_tol_Omega, max_step_s, n_steps, verdict
FROM run WHERE run_id = 'B_final';
```

## 6. How the data flows in and out (the lifecycle)

1. `mercury_rs` (Rust) runs a simulation and writes
   `data/runs/<run_id>/samples.csv`, `events.csv`, `manifest.json` — plus
   `restart.csv` for run B and `branches.csv` for run C. All numbers printed via the
   engine's C-style formatters, so re-runs are byte-identical.
2. The notebook (Python stdlib `sqlite3` + `csv` + `json`) creates the tables of
   Section 3, loads every CSV/manifest, and prints the row counts as proof
   (expected magnitudes: `run` 6 rows; `sample` ≈ 60,000–80,000 rows total —
   ~10,000 secular + ~25,000 resonant + ~25,000 dense capture-window + the short
   runs; `event` ≈ 20 rows; `branch` 64; `target` 5).
3. The notebook retrieves (Section 5's queries), checks the acceptance targets, and
   hands the database to `gui/bake_page.py`, which reads it (Python stdlib
   `sqlite3`) and bakes the self-contained browser page.
4. Anyone can later open the database directly with any SQLite viewer, or with three
   lines of Python — the file is the permanent, queryable record of the simulated
   history of Mercury.
