# PLANET_MERCURY_PROVENANCE.md — the planet_Mercury sub-project import

This document records, in the repository's provenance convention, what the
top-level `planet_Mercury/` directory is, where it came from, how it was
verified, and every deviation made along the way. It is self-contained.

## 1. What was added

`planet_Mercury/` — a complete, verified teaching-and-research sub-project that
simulates the historical tidal despinning and 3:2 spin-orbit resonant capture of
planet Mercury as a two-body Sun–Mercury system (point-mass Sun; Mercury as an
extended, deformable, triaxial body), using THIS repository's vendored pure-Rust
SUNDIALS 7.8.0 CVODE as the only integrator. Contents:

- `mercury_rs/` — pure-Rust compute crate (BSD-3-Clause; `#![forbid(unsafe_code)]`,
  `#![deny(warnings)]`, zero external dependencies; path-depends on
  `../../sundials_rs/crates/{sundials_core,cvode_rs}`). Five-state ODE
  y = [a, e, M, θ, Ω]: Kepler's law, the Hut (1981) constant-time-lag tidal
  torque with its orbital back-reaction, and the Goldreich–Peale triaxial
  "handle" torque. CVODE BDF + Newton + dense at rtol 1e-12,
  atol [1e-3, 1e-6, 1e-10, 1e-10, 1e-14], max step 864000 s, with exact
  root-finding on spin-ratio crossings. 15 analytic-expectation unit tests.
- `mercury_crosscheck/` — small GPL-3.0-or-later crate (links this repository's
  GPL `rebound_rust`/`reboundx_rust` ports, exactly like those imports
  themselves): independently reproduces the tidal despin rate with the
  REBOUNDx `tides_spin` (Eggleton/Lu) formulation — agreement 0.25%, after
  documenting that REBOUNDx's "constant time lag" is exactly HALF of Hut's τ
  (σ = 4τG/(3R⁵k₂); derived by hand in the source).
- `notebook/` — the executed Jupyter deliverable `mercury_tidal_locking.ipynb`
  (Python 3 ipykernel, standard-library-only cells driving the Rust binary by
  subprocess, per this repository's own notebook convention) plus its
  deterministic builder, batch runner, and structure auditor.
- `gui/mercury_orbit.html` + `gui/bake_page.py` — a self-contained browser
  display page in this repository's recorded-player style (all data embedded,
  vanilla JS + canvas, zero network requests): orbit animation with the spin
  "handle," a spin/orbit-ratio dial settling on exactly 1.5, libration,
  period, and angular-momentum-ledger plots, play/scrub controls.
- `plan/` (eight planning documents), `SOURCE_SPECIFICATION.md` (verbatim copy
  of the consolidated physics specification), and
  `INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md` (complete self-contained guide).
- Deliberately NOT committed: `planet_Mercury/data/` (~267 MB of CSV/SQLite the
  notebook regenerates deterministically, bit-for-bit, on one run) and build
  `target/` dirs. No `.gitignore` files were added anywhere in the sub-project.

## 2. Where the physics came from

The model implements *"Mercury 3:2 Spin-Orbit Resonant Capture Provenance
Specification"* (Patrick Nash, 2026-08-24), as consolidated and corrected in
`planet_Mercury/SOURCE_SPECIFICATION.md` (2026-08-25): Hut (1981) CTL tides,
Goldreich & Peale (1966) capture mechanics, spec constants verbatim
(k₂ = 0.12, τ = 100 s, (B−A)/C = 1e-4, C = 0.34mR², a₀ = 5.790905e10 m,
e₀ = 0.20563, Ω₀ = 1.5e-4 rad/s). The specification's four inherited equation
errata (E1–E4, documented in that file) are corrected; the corrected trio
provably conserves C·Ω + m·n·a²·√(1−e²) exactly in the secular part (verified
symbolically during planning and enforced by unit test).

## 3. Verification (all measured on this machine, Apple Silicon, 2026-08-30)

- Engine untouched: after the one-line root `Cargo.toml` change (appending
  `"planet_Mercury"` to the workspace `exclude` list; prior copy in
  `.backups/2026-08-30/`), `cargo metadata` shows the original four workspace
  packages only; `cargo build --workspace --all-targets` warning-free; all 622
  workspace tests green; `bash tools/macos_verify_physics.sh` exit 0.
- Crate gates: `cargo build --release` warning-free and `cargo test --release`
  15/15 in BOTH the workspace original and this in-repo mirror; the mirror's
  run-E seam validation reproduces the original bit-for-bit (despin-rate
  agreement between the secular and full models: 0.12%).
- Science results (the executed notebook's own gauntlet asserts all of these):
  spec-literal run A: Ω/n falls only 181.4 → 178.9 in 10 Myr — the spec's own
  constants imply ~4.7 Gyr of braking (Finding F1), motivating the documented
  1000× time compression of the movie runs; canonical capture into 3:2 at
  4.6499 Myr movie-time; 64-branch phase sweep: 10/64 = 15.6% ± 4.5% captured
  (Goldreich–Peale's simple estimate ~7%); settled-era mean spin ratio
  1.5000000498; final periods 87.968 d / 58.642 d / 175.91 d vs observed
  87.969 / 58.646 / 175.938 (instantaneous samples carry the residual
  ~2e-4 libration wiggle; the time-averaged 2/3 lock holds to 3.3e-8);
  libration swing decays 6.22 → 0.15 rad across the notebook's 12 bins;
  angular-momentum ledger: pre-capture drift ≤ 1e-9 budget, locked era checked
  against the model's own predicted secular leak (the spec's handle torque has
  no orbital back-reaction — see deviation DEV-7 below).
- Notebook: executed end-to-end headless twice; the two executed files are
  byte-identical (determinism gate); the structure auditor passes all rules.
- Display page: baked twice byte-identically; zero network references; loads
  with zero console errors; capture-jump, scrub, HUD readouts and all four
  plots verified in a real browser.

## 4. Deviations register (every one documented and justified)

Approved-plan decisions: D1 build-then-mirror layout (this import); D2 time
compression S = 1000 (spec's 10-Myr window vs its own ~4.7-Gyr braking
timescale); D3 staged integration (the handle torque's ~6-hour wiggle would
need ~1e11 steps; it averages to zero far from resonance, switched on at
Ω/n = 2.2, seam-validated); D4 64-branch phase sweep (capture is genuinely
probabilistic); D5 notebook/display technology per this repository's own
conventions; D6 no `.gitignore` files added.

Build-time findings (full detail in `planet_Mercury/plan/07_PROVENANCE_AND_DEVIATIONS.md`):
DEV-1 the seam checkpoint was moved from Ω/n = 3.0 — itself the 3:1 resonance,
which the validation run caught — to the non-resonant 2.7; DEV-2 run D's
libration period is checked with the exact large-amplitude pendulum formula
(fresh capture librates near the separatrix); DEV-3 sweep branches and the
canonical continuation use numerically identical solver paths (capture is
bit-level path-sensitive); DEV-4 the 2/3-lock statement is checked on the
time-averaged spin (a locked planet librates forever); DEV-5 the REBOUNDx τ
convention is half of Hut's (mapping documented in mercury_crosscheck);
DEV-6 regenerable data not committed; DEV-7 in lock the handle torque's mean
is nonzero by necessity, so the model leaks total angular momentum secularly
(~2e-10 per compressed Myr) — the ledger check is split into the two eras
honestly, and the f64 quantization that freezes the locked orbit's recorded
drift is documented; DEV-8 two reviewer notes documented-not-changed (θ's
magnitude-degraded absolute tolerance; the γ-unwrap margin at the dense
trigger), each harmless as configured and pinned with comments; DEV-9 a
review-driven manifest-clearing hardening briefly deleted a sibling run's
manifest on read access — caught by the end-to-end pipeline and fixed by
splitting read-safe `run_dir` from owning-run `fresh_run_dir`.

## 5. Reproduction

From a clone of this repository:
`cd planet_Mercury/mercury_rs && cargo build --release && cargo test --release`
(15 tests), then follow `planet_Mercury/INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md`
— one top-to-bottom notebook run (~1 hour) regenerates every dataset, database,
check, and the display page, deterministically.
