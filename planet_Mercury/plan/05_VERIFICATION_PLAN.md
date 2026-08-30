# Planet Mercury Tidal-Locking Plan — Document 5 of 8: The Verification Plan (Every Check That Must Pass)

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science. Everything
needed is inside this document.
**Status:** PLAN (awaiting approval). These are the gates the build must pass; none
have been run yet.

Words expanded once for this document: **CSV** = comma-separated values (a data table
saved as plain text); **SQL** = Structured Query Language (database queries);
**HUD** = heads-up display (the on-screen number readouts of the display page);
**BDF** = backward differentiation formulas (the solver method); an **e-fold** is one
shrink of an exponentially decaying quantity by the factor e ≈ 2.718.

Philosophy (inherited from the engine project): every check compares against an
**analytic expectation** (a number derivable with pencil and paper) or a **pinned
observed value** — never against "whatever the program printed last time." Every
command is run with `2>&1 | tee <log>` and the log is read before proceeding. A gate
either passes completely or the build stops and the failure is reported verbatim.
The one deliberately statistical check (the capture-odds sweep) states its chance
element and its contingency explicitly (gate 2.3).

---

## Gate 0 — Baseline: the engine is green BEFORE we start

Run from the engine repository
(`/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/`):

| # | Check | Command | Must show |
|---|---|---|---|
| 0.1 | Workspace builds warning-free | `cargo build --workspace --all-targets 2>&1 \| tee /tmp/g0_build.log` | zero warnings, zero errors |
| 0.2 | All engine tests pass | `cargo test --workspace 2>&1 \| tee /tmp/g0_test.log` | 622 tests, all green |
| 0.3 | Physics byte-identity gate | `bash tools/macos_verify_physics.sh 2>&1 \| tee /tmp/g0_phys.log` | exit 0 |
| 0.4 | Git tree clean & synced | `git status -sb` | `## main...origin/main`, nothing else |

## Gate 1 — The `mercury_rs` crate (unit correctness)

Run from `planet_Mercury/mercury_rs/`. Build: `cargo build --release --all-targets`
(zero warnings — the crate root carries `#![deny(warnings)]`, so a warning IS a build
failure). Then `cargo test`, which must include ALL of these named tests, each with
its analytic expectation written into the test (every constant needed to reproduce
them by hand is in the numbers wall at the end of this document):

| # | Test | Analytic expectation |
|---|---|---|
| 1.1 | `mean_motion_matches_keplers_third_law` | n(a₀) = 8.2669×10⁻⁷ rad/s from n = √(G(M☉+m)/a₀³); P_orb = 2π/n = 87.97 days |
| 1.2 | `kepler_solver_inverts_keplers_equation` | for a grid of (e, E) pairs, solving M = E − e·sinE back for E recovers E to 1×10⁻¹³; f and r match the closed forms |
| 1.3 | `kepler_solver_circular_orbit_is_identity` | at e = 0: E = M, f = M, r = a exactly |
| 1.4 | `hut_polynomials_at_zero_eccentricity_are_one` | f₁(0)=f₂(0)=f₃(0)=f₄(0)=f₅(0)=1 exactly |
| 1.5 | `pseudo_synchronous_ratio_at_mercury_e` | f₂/f₁ at e = 0.20563 = 1.2560 ± 0.0005 (the "tides alone park at 1.256 n" number) |
| 1.6 | `pseudo_synchronous_equals_three_halves_at_e285` | f₂/f₁ = 1.5000 ± 0.0005 at e = 0.285 (run D's foundation) |
| 1.7 | `tidal_torque_sign_brakes_fast_spin_and_spins_up_slow` | ⟨T⟩ < 0 for Ω = 2n; ⟨T⟩ > 0 for Ω = n (since Ω_eq = 1.256 n > n) |
| 1.8 | `despin_rate_matches_linearized_prediction` | a short CVODE run at Ω = 100n, secular stage: measured dΩ/dt within 0.1% of −K·(Ω·f₁ − n·f₂)/C computed by hand |
| 1.9 | `angular_momentum_ledger_closes_in_secular_stage` | over a secular run long enough that the transferred spin momentum Δ(C·Ω) is well resolved, \|Δ(C·Ω) + ΔL_orb\| ≤ 1×10⁻⁸ · \|Δ(C·Ω)\| — i.e., the ledger closes to one part in 10⁸ **of the momentum actually moved** (the sharp baseline; the corrected Hut equations conserve the total exactly, proven symbolically during the plan audit, so this measures only solver error — and the uncorrected source equations would fail it catastrophically, proving errata E2/E3 matter) |
| 1.10 | `reanchoring_preserves_gamma_exactly` | γ = 2θ − 3M identical to the last bit before/after the (−3πj, −2πj) re-anchoring |
| 1.11 | `triaxial_torque_zero_when_aligned` | T_tri = 0 at θ = f and maximal at θ − f = π/4 |
| 1.12 | `csv_row_roundtrips_through_fmt_e` | a written sample row re-parses to the same f64s |
| 1.13 | `restart_reproduces_bit_identical_state` | run to t₁, save restart, reload, run to t₂ — matches running straight to t₂ at the sample points within tolerance; the reload itself is bit-exact |
| 1.14 | `libration_frequency_matches_goldreich_peale` | in run-D conditions, measured libration period ≈ 2π/(n·√(3(B−A)/C·H(e))) within 5% (H(e) = (7/2)e − (123/16)e³) |
| 1.15 | `root_function_stops_at_the_crossing` | the CVodeRootInit root on (Ω − 1.5·n) returns CV_ROOT_RETURN with \|Ω/n − 1.5\| < 10⁻⁹ at the reported time |

## Gate 2 — The science runs (physics acceptance)

Each run's binary self-checks print PASS lines and end `SUCCESS` (nonzero exit
otherwise). The specific acceptance criteria:

| # | Run | Acceptance |
|---|---|---|
| 2.1 | A (spec-literal, τ = 100 s) | final Ω/n between 170 and 181 (barely moved — Finding F1 confirmed numerically); NO capture event; verdict SUCCESS (success = honest confirmation, not capture) |
| 2.2 | B (movie, τ_eff = 1.0×10⁵ s) | despin passes 2:1 without sticking (or the event log documents a 2:1 capture and the sweep restarts above 2:1 — decided by the data, recorded either way); restart state saved by the root-finder when Ω/n first reaches 1.6, at ≈ 4.6 Myr (±20%) |
| 2.3 | C (64-branch sweep) | every branch terminates with a definite outcome; **pass condition: at least one branch captures, allowing at most one documented finer-grid contingency re-sweep** (at ≈ 7% odds per branch there is ≈ 1% chance a 64-branch sweep gets zero captures). The measured fraction is then **reported with its 95% binomial uncertainty band** next to the 7.0% theoretical estimate (64 flips at 7% typically give 2–8 captures) — the comparison to theory is informational, not pass/fail, and the notebook says so honestly |
| 2.4 | B-final (canonical) | capture event present at ≈ 4.7 Myr (±20%); thereafter mean Ω/n − 1.5 within 1×10⁻⁴; γ librates with decaying swing; final P_orb, P_rot, P_solar within 0.1% of 87.969 / 58.646 / 175.938 Earth days; P_rot/P_orb = 2/3 to ≥5 significant figures |
| 2.5 | D (e = 0.285 encore) | capture occurs (the guaranteed-arrival mechanism); libration confirmed |
| 2.6 | E (seam validation) | window-averaged secular drift rates of a, e, Ω from the full model match the secular model within 1% each |
| 2.7 | Cross-run sanity | run B's secular-stage despin e-folding time ≈ 710 Myr / 1000 = 0.71 Myr within 10% (matches the pencil-and-paper C/(K·f₁) with the movie's k₂τ) |
| 2.8 | Optional cross-check (`mercury_crosscheck`, if built) | the REBOUNDx-port `tides_spin` despin rate for the same Sun–Mercury parameters agrees with `mercury_rs`'s secular dΩ/dt within a documented tolerance (an independent implementation of the same Hut physics) |

## Gate 3 — The database

| # | Check |
|---|---|
| 3.1 | Ingest completes; row counts printed match the CSV line counts exactly (sweep branches contribute rows only to `branch`, by design) |
| 3.2 | `run` has one row per executed run, each `verdict = 'SUCCESS'` |
| 3.3 | Rebuild determinism: deleting the .sqlite3 and re-ingesting yields identical query results (spot-checked via the notebook's seven queries Q1–Q7) |
| 3.4 | Referential integrity (every stored row points at a run that really exists): every `sample.run_id`, `event.run_id`, and `branch.run_id` exists in `run` (SQL join check, zero orphans) |
| 3.5 | Ledger invariant in SQL, two eras (refined by the build's adversarial review): PRE-CAPTURE drift \|L_tot − L_tot(0)\|/L_tot(0) ≤ **1×10⁻⁹** (the model conserves exactly there); LOCKED-ERA drift compared against the model's own predicted secular leak K·n·(1.5f₁−f₂)·Δt/L_tot (the spec's handle torque has no orbital back-reaction, so in lock L_tot leaks at ~2×10⁻¹⁰ per compressed Myr by design; recorded values are further suppressed by f64 orbit quantization — both stated honestly) |

## Gate 4 — The display page

| # | Check |
|---|---|
| 4.1 | Bake determinism, spelled out (run from `planet_Mercury/`): `python3 gui/bake_page.py` then `cp gui/mercury_orbit.html /tmp/bake1.html`, then `python3 gui/bake_page.py` again, then `cmp /tmp/bake1.html gui/mercury_orbit.html` — `cmp` must print nothing (byte-identical) |
| 4.2 | Self-containedness: a scan of `mercury_orbit.html` finds no `fetch(`, no `XMLHttpRequest`, no `WebSocket`, no external `src=`/`href=` resources (plain-text mentions of citations excepted, verified by eye) |
| 4.3 | No integration in the page: code review confirms playback only replays embedded samples |
| 4.4 | Browser smoke test (real browser, automated): page loads with zero console errors; Play advances the year readout; scrub to the capture time shows the ratio dial at 1.500…; the libration plot shows a decaying rocking curve; the HUD's three periods match the observed 87.969 / 58.646 / 175.938 within 0.1% |
| 4.5 | File size sane (< 10 MB) |

## Gate 5 — The notebook

| # | Check |
|---|---|
| 5.1 | `MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_tidal_locking.ipynb` → every non-interactive cell executes without error; real outputs embedded |
| 5.2 | `python3 check_notebook.py mercury_tidal_locking.ipynb` → all structure rules pass (how-to-run text present; every code cell has a ≥80-character markdown lead-in; no references to any other file for required information; required section headings present; every non-interactive code cell carries a real execution count) |
| 5.3 | The §9 verification-gauntlet cell prints one PASS per acceptance target and no FAILs |
| 5.4 | Re-execution determinism: running 5.1 twice produces byte-identical .ipynb files |
| 5.5 | Fresh-eyes walk-through: an agent with no build context follows INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md literally in a clean shell, from Part B step 1, and reaches the end with zero improvisation (any needed improvisation = a doc bug to fix) |

## Gate 6 — Adversarial review (multi-agent, before any commit)

Independent reviewer agents, each tasked to REFUTE, not confirm:

| # | Reviewer | Attack surface |
|---|---|---|
| 6.1 | Physics reviewer | every formula in `hut.rs`/`rhs.rs` vs. the corrected Hut/Goldreich–Peale forms; arithmetic order preserved; unit consistency |
| 6.2 | Rust reviewer | borrow-guard discipline around N_VGetArrayPointer; no `unwrap` on solver `Option`s; error style `Result<_, String>` with actionable messages; zero `unsafe`; no external deps in Cargo.lock |
| 6.3 | Determinism reviewer | hunts timestamps, HashMap iteration order, platform-dependent formatting, uncontrolled randomness |
| 6.4 | Docs reviewer | every number quoted in the markdown files and notebook prose is measured from the actual runs (never remembered); the self-containedness rule (no "see other file" for required information) holds in every shipped document |

## Gate 7 — Mirror, commit, push

| # | Check |
|---|---|
| 7.1 | Engine root `Cargo.toml` backed up to `.backups/<date>/` before editing |
| 7.2 | Mirror copy re-builds and re-tests green in place inside the engine repo (Gate 1 repeated there); run E re-executed as a smoke test |
| 7.3 | Engine gates re-run AFTER the root Cargo.toml edit: workspace build warning-free; 622 tests green; `tools/macos_verify_physics.sh` exit 0; `cargo metadata` shows `planet_Mercury` absent from the workspace package list (the `exclude` works) |
| 7.4 | `git status --porcelain` shows ONLY: `planet_Mercury/` files, the one-line `Cargo.toml` change, and `PLANET_MERCURY_PROVENANCE.md` — nothing else; zero `.gitignore` files anywhere under `planet_Mercury/` (assignment rule); no `target/` content staged |
| 7.5 | Secret scan over every staged file (pattern scan for keys/tokens/credentials) comes back clean |
| 7.6 | Commit message written; `git push` succeeds; `git ls-remote origin main` equals the local HEAD hash |
| 7.7 | Post-push clone test: `git clone` into a scratch directory; inside the clone, `planet_Mercury/mercury_rs` builds and its tests pass (proving the pushed artifact is complete and self-sufficient) |

## The numbers wall (single, self-sufficient source for every target used above)

Constants: G = 6.67430×10⁻¹¹ m³·kg⁻¹·s⁻² (Newton's constant), M☉ = 1.98847×10³⁰ kg
(Sun), m = 3.3011×10²³ kg (Mercury), R = 2.4397×10⁶ m (Mercury's radius),
a₀ = 5.790905×10¹⁰ m, e₀ = 0.20563, Ω₀ = 1.5×10⁻⁴ rad/s, C = 0.34·m·R²
(moment of inertia), (B−A)/C = 1.0×10⁻⁴ (equator lopsidedness), k₂ = 0.12,
τ = 100 s (spec) / 1.0×10⁵ s (movie).

The five Hut (1981) eccentricity polynomials used everywhere above:
f₁(e) = (1 + 3e² + (3/8)e⁴)/(1−e²)^(9/2); f₂(e) = (1 + (15/2)e² + (45/8)e⁴ +
(5/16)e⁶)/(1−e²)⁶; f₃(e) = (1 + (31/2)e² + (255/8)e⁴ + (185/16)e⁶ +
(25/64)e⁸)/(1−e²)^(15/2); f₄(e) = (1 + (3/2)e² + (1/8)e⁴)/(1−e²)⁵;
f₅(e) = (1 + (15/4)e² + (15/8)e⁴ + (5/64)e⁶)/(1−e²)^(13/2).
At e₀ = 0.20563: f₁ = 1.3695, f₂ = 1.7200.

| Quantity | Value | Where it comes from |
|---|---|---|
| n at a₀ | 8.2669×10⁻⁷ rad/s | √(G(M☉+m)/a₀³) |
| P_orb target | 87.969 Earth days | 2π/n; matches radar/ephemeris Mercury |
| P_rot target | 58.646 Earth days | (2/3)·P_orb — the 3:2 lock |
| P_solar target | 175.938 Earth days | 1/\|1/P_rot − 1/P_orb\| = 2·P_orb exactly in a 3:2 lock |
| Ω_eq/n at e₀ | 1.256 | f₂(e₀)/f₁(e₀) |
| C | 6.68×10³⁵ kg·m² | 0.34·m·R² |
| K (spec-literal) | 2.18×10¹⁹ kg·m²/s | 3GM☉²R⁵k₂τ/a₀⁶ at k₂τ = 12 s |
| Despin e-fold, spec-literal | ≈ 710 Myr | C/(K·f₁) |
| E-folds from 181.4 n to 1.5 n | 6.6 | ln((181.4 − 1.256)/(1.5 − 1.256)) — decay is TOWARD Ω_eq = 1.256 n, not toward zero |
| Spec-literal braking total | ≈ 4.7 Gyr | 6.6 × 710 Myr (hence Finding F1: the 10-Myr window is ~470× too short) |
| Despin e-fold, movie (S = 1000) | ≈ 0.71 Myr | 710 Myr / 1000 |
| Movie 3:2 crossing time | ≈ 4.7 Myr | 6.6 × 0.71 Myr (stage handover ≈ 3.7 Myr, 2:1 crossing ≈ 3.9 Myr) |
| Libration frequency at 3:2 | ω_lib ≈ 0.014·n ≈ 1.16×10⁻⁸ rad/s | n·√(3·(B−A)/C·H(e)), H = (7/2)e − (123/16)e³ ≈ 0.653 |
| Capture probability estimate | ≈ 7.0% | Goldreich & Peale: P = 2/(1 + πV/(2W)), V = K·n·(1.5f₁−f₂) ≈ 6.0×10¹² N·m, W = K·f₁·ω_lib ≈ 3.5×10¹¹ N·m (ratio independent of the compression S) |
| Angular-momentum budget | 1×10⁻⁹ of L_tot(0) | gate 3.5 / notebook §9; predicted physical wobble ~10⁻¹¹, rest is solver headroom |
