# Planet Mercury Tidal-Locking Plan — Document 7 of 8: Provenance and Deviations Register

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science. Everything
needed is inside this document.
**Status:** PLAN (awaiting approval).

**What this document is:** the assignment requires that any issue with following the
instructions exactly be explained in a provenance markdown file. This is that file.
**Bottom line: the work CAN proceed — nothing here forces a stop.** Every issue found
during planning is listed below with its evidence and its resolution.

Symbol key for the formulas quoted below: a = orbit size (semi-major axis), e = orbit
ovalness (eccentricity), M = orbit clock angle (mean anomaly), θ = spin angle,
Ω = spin rate, n = the orbit's average angular speed (mean motion), m = Mercury's
mass, M☉ = the Sun's mass, R = Mercury's radius, C = Mercury's spin moment of inertia
(resistance to spin changes), B−A = how lopsided Mercury's equator is, k₂ = the Love
number (squishiness), τ = the tidal time lag (sluggishness), K = the tidal-brake
strength 3·G·M☉²·R⁵·k₂·τ/a⁶, f₁…f₅ = standard eccentricity factors from Hut (1981).
An **e-fold** is one shrink of an exponentially decaying quantity by the factor
e ≈ 2.718; **BDF** = backward differentiation formulas, the solver method used.

## Part 0 — The decisions submitted for approval, restated in full

| # | Title | Substance |
|---|---|---|
| D1 | Build location vs. push target | Build and verify at the required outer path `…/planet_Mercury/`; mirror the verified tree into the engine git repository as a top-level `planet_Mercury/` directory (that repo's exact precedent for sub-projects), adjust one dependency-path line, re-verify in place, commit by explicit path, push to `origin` — which is exactly the required URL |
| D2 | Time compression | Multiply the tidal-brake strength k₂τ by S = 1000 (τ: 100 s → 1.0×10⁵ s) so the ~4.7-billion-year braking story fits the 10-million-year window, with capture landing at ≈ 4.7 Myr; the spec-literal run is kept alongside for honesty |
| D3 | Staged integration | Switch the oscillating "handle" torque off far from resonance (where it averages to zero) and on for the resonant era; both stages are CVODE BDF with Newton iteration and a dense linear solver at the spec tolerances (relative 1.0×10⁻¹²; absolute [1.0×10⁻³, 1.0×10⁻⁶, 1.0×10⁻¹⁰, 1.0×10⁻¹⁰, 1.0×10⁻¹⁴] for [a, e, M, θ, Ω]; max step 864,000 s); a seam-validation run proves the two stages agree where they hand over |
| D4 | The capture coin-flip | Run a 64-branch spin-phase sweep at the 3:2 crossing (one finer-grid re-sweep if no branch captures — about a 1% possibility); the first captured branch, its phase offset recorded, becomes the canonical locked history; the measured capture fraction is displayed beside the ≈ 7% theory with honest error bars |
| D5 | Notebook & display technology | The deliverable notebook is a Python 3 (ipykernel) Jupyter notebook whose standard-library-only code cells drive the pure-Rust `mercury_rs` program as a subprocess (the engine project's exact convention for all 109 of its own notebooks); the display is one self-contained browser HTML page in the engine's recorded-player style (all data embedded, no frameworks, no network) |
| D6 | .gitignore rule | No `.gitignore` files are created or committed anywhere under `planet_Mercury/` (the assignment forbids pushing them; the repository's existing ignore rule for `target/` build folders already keeps build output out of git) |

## Part 1 — Where the problem specification came from (chain of custody)

1. The physics problem originates in *"Mercury 3:2 Spin-Orbit Resonant Capture
   Provenance Specification"* (Patrick Nash, self-addressed email, 2026-08-24),
   fragmented by a Gemini agent into 15 files, then consolidated, deduplicated,
   verified, and corrected into the workspace file `Planet Mercury by Fable.md`
   (2026-08-25). That consolidated document was read in full, is the authoritative
   physics source for this plan, and will be copied verbatim into the deliverable as
   `planet_Mercury/SOURCE_SPECIFICATION.md` so the pushed repository carries its own
   source.
2. The compute/display engine is the repository
   `rustSolveIt_macos-silicon_SUNDIALS_7_8_0` (macOS / Apple Silicon port), whose only
   integration backend is a vendored (stored-inside-the-project), read-only, pure-Rust
   port of SUNDIALS 7.8.0. Its conventions (zero unsafe / zero external dependencies /
   zero warnings; notebooks as Python-3-driving-Rust; self-contained browser display
   pages; C-style number formatting; provenance files for imported work) were
   established by direct inspection during planning and are adopted wholesale.
3. The tidal model is Hut (1981) constant-time-lag theory; the capture theory is
   Goldreich & Peale (1966); the specification itself cross-references the REBOUND /
   REBOUNDx frameworks, whose pure-Rust ports also live in the engine repository and
   are used only as an optional independent cross-check (they contain the tidal
   physics but no permanent-triaxiality torque, so they cannot reproduce the 3:2
   capture itself — verified by direct source inspection).
4. This plan itself was adversarially audited before presentation: four independent
   reviewer agents recomputed every physics number, verified every engine claim
   against real source files, and hunted contradictions; their confirmed findings
   (including a wrong e-fold count in the first draft of the timescale audit, fixed
   to the numbers quoted below) are incorporated.

## Part 2 — Errors inherited from the source specification (all corrected)

These four errors were found in the source document by checking against Hut (1981) —
the paper it cites — and by dimensional analysis; an independent audit re-verified all
four corrections, including a symbolic proof that the corrected equations conserve
total angular momentum exactly.

| # | What was wrong | The correction |
|---|---|---|
| E1 | The tidal-torque factor f₁(e) was written with denominator (1−e²)^(15/2) | (1−e²)^(9/2). With the wrong exponent the "tides-alone" equilibrium spin would sit below the orbital rate at every eccentricity, making 3:2 capture impossible — contradicting the document's own goal |
| E2 | da/dt had a dimensionally wrong prefactor (−2aK/(G·M☉·m)) and the wrong sign | da/dt = +(2K/(m·n·a))·[Ω·f₂ − n·f₃] — a fast-spinning planet pushes its orbit outward (same physics as the Moon receding from Earth); reduces to the classical circular-orbit formula at e = 0 |
| E3 | de/dt scrambled Hut's equation: the 11/18 factor on the wrong term, factors swapped, one factor matching nothing in Hut, and E2's dimensional error repeated | de/dt = +(9Ke/(m·n·a²))·[(11/18)·Ω·f₄ − n·f₅] with f₄ = (1+(3/2)e²+(1/8)e⁴)/(1−e²)⁵ and f₅ = (1+(15/4)e²+(15/8)e⁴+(5/64)e⁶)/(1−e²)^(13/2) — gives eccentricity damping for slow rotators, as physics requires |
| E4 | The initial spin Ω₀ = 1.5×10⁻⁴ rad/s was described as a "~1.16 day" period | 2π/Ω₀ = 41,888 s ≈ 11.6 hours; the numeric value is kept authoritative (the outcome is insensitive to this choice) |

## Part 3 — New findings made during this planning (the deviations, each justified)

| # | Finding (with evidence) | Resolution (decision #) |
|---|---|---|
| F1 | **The spec cannot capture in its own window.** Its constants give a tidal-braking e-folding time of C/(K·f₁) ≈ 710 million years (C = 6.68×10³⁵ kg·m², K = 2.18×10¹⁹ kg·m²/s at k₂τ = 12 s, f₁ = 1.3695). The spin decays toward the tides-only resting rate 1.256·n — not toward zero — so reaching the 1.5·n mark from 181·n takes ln((181.4−1.256)/(1.5−1.256)) ≈ 6.6 e-folds ≈ **4.7 billion years**, while final_time is 10 million years (~470× too short). The spec's own checks verified periods and dimensions, never this timescale | Honor the spec's stated intent ("dissipation deliberately set strong… tractable") with a documented **1000× compression** of k₂τ, putting the 3:2 crossing at ≈ 4.7 Myr with >5 Myr of locked-state display left in the window; keep a spec-literal run that shows the honest slow result (**D2**) |
| F2 | **The spec run is computationally impossible as written.** The triaxial torque oscillates with a ~6-hour period at the initial spin; resolving it at the spec tolerances for 10 Myr implies ~10¹¹ solver steps, versus the spec's own max_steps = 5×10⁸ | **Staged integration**: the oscillating torque is switched off far from resonance, where its orbit-average is essentially zero (standard practice in the research literature), and the full system runs near/after the crossings; a seam-validation run proves the two agree where they hand over; every stage is still CVODE BDF at the spec tolerances listed under D3 above (**D3**) |
| F3 | **Capture is probabilistic, not guaranteed.** Goldreich–Peale theory gives ≈ 7.0% capture odds per 3:2 crossing at e = 0.20563 for this tidal model (independently re-derived during the audit; the odds do not depend on the compression factor, which cancels out) — a single run with an arbitrary phase usually sails through; the spec implicitly assumes a capturing phase | A documented **64-branch spin-phase sweep** (~1% chance of zero captures → one documented finer-grid re-sweep); the measured capture fraction is displayed beside the theory with binomial error bars; the first captured branch becomes the canonical history, its phase offset recorded for bit-reproducibility (**D4**) |
| F4 | **The build location and the push target belong to different trees.** The required build folder (workspace level) is not a git repository (verified: a standard git query fails there); the required push URL is already the `origin` of the engine repository nested beside it (verified) | Build at the required outer location; **mirror** the verified result into the engine repository as top-level `planet_Mercury/` (that repo's exact precedent for hosting sub-projects), adjust one dependency path line, re-verify in place, commit by explicit path, push (**D1**) |
| F5 | **Two instruction sets disagree about `.gitignore`.** The engine repository's precedent tracks per-crate `.gitignore` files, but the assignment says "do not push .gitignore files" | The assignment wins: **no `.gitignore` files are created or committed** anywhere under `planet_Mercury/`; the repository's existing ignore rule for any `target/` build folder already keeps build output out of git (**D6**) |

## Part 4 — Smaller notes for the record (no action needed)

| # | Note |
|---|---|
| N1 | The assignment's HARD RULES name the C reference as `/home/nsh/Developer/sundials-7.8.0/` — a Linux path that does not exist on this macOS machine (verified). The actual reference tree exists at `/Users/nsh/Developer/sundials-7.8.0/` (verified), which is also where the engine's own documentation points. The macOS path is used. |
| N2 | "Tidal locking" for Mercury means the **3:2 spin-orbit resonance** (three spins per two orbits), not the Moon-style 1:1 lock; the source specification, the physics literature, and this plan all use it that way. The plan's displays make the distinction explicit for students. |
| N3 | "rustSolveIt Jupyter notebook" is implemented the way every one of the engine project's own 109 notebooks is implemented — see decision D5 in Part 0 for the full statement. All physics computation is Rust; the notebook's standard-library Python is the conductor, database loader, and verifier. This is the project's native, verified notebook form. |
| N4 | The workspace's global coding-rules file for TypeScript does not apply: this project contains no TypeScript. The display page's small amount of JavaScript follows the engine's own established page conventions (plain "vanilla" JavaScript, no frameworks, no network). |
| N5 | The optional cross-check crate links the workspace's REBOUND/REBOUNDx pure-Rust ports, which are GPL-3.0-or-later licensed (the engine workspace itself is BSD-3-Clause and already hosts those GPL crates as excluded sub-projects). The cross-check crate is therefore GPL-licensed and stands alone, exactly like the precedent. It validates only the tidal-braking rate (the ports contain no triaxial torque — verified by source inspection — so they cannot reproduce capture, and are not used for it). If it is dropped for time, the drop is recorded here. |
| N6 | The outer workspace's `p*.txt` files cannot be pushed even by accident: they live outside the engine repository's work tree, and git refuses to stage paths outside the work tree. No secrets exist in any planned artifact; a secret scan is still run before the commit as a gate. |
| N7 | Jupiter's perturbations and Einstein's general-relativistic perihelion correction are excluded by explicit assignment instruction; both are named in the plan as the content of the second test. (For the record: the engine's REBOUNDx port already contains `gr_potential` — the standard tool for the GR correction — so the second test has a ready ingredient.) |
| N8 | Mercury's obliquity (spin-axis tilt) is fixed at 0° per the specification; the real value is within ~2 arcminutes of that. |
| N9 | The consolidated specification file lives at the workspace root, outside the deliverable tree; to keep the published repository self-explanatory, a verbatim copy ships inside the project as `SOURCE_SPECIFICATION.md` (see Part 1, item 1). |

## Part 4b — As-built addendum: deviations discovered DURING the build (all verified)

The build itself surfaced five more findings; each was fixed, verified, and is
recorded here. None changes the physics model or the solver contract.

| # | Finding during the build | Resolution |
|---|---|---|
| DEV-1 | The plan put the staging-seam validation (run E) at spin ratio Ω/n = 3.0 — which is itself the **3:1 spin-orbit resonance** (resonances sit at every half-integer ratio). The full model locked there and the drift comparison was meaningless — the seam test caught a real plan bug | Run E moved to the non-resonant ratio **2.7**; and since the two stages differ ONLY in the spin equation, the seam gate compares the despin rate (measured agreement: 0.12%, gate 1%) plus end-state sanity on a and e |
| DEV-2 | The plan checked run D's libration period against the small-amplitude Goldreich–Peale formula — but a **fresh capture librates at large amplitude** near the separatrix (measured 2.73 rad), where the pendulum period is much longer | The check uses the exact large-amplitude pendulum period T = (4/ω_lib)·K(sin(γmax/2)) with K the complete elliptic integral (arithmetic–geometric mean), amplitude measured from the data. Result: measured 27.52 yr vs predicted 29.25 yr — 5.9% agreement (gate 15%); the small-amplitude formula itself is separately unit-tested at small amplitude to 5% |
| DEV-3 | A sweep branch's capture verdict did NOT initially reproduce in the canonical continuation: capture is **bit-level path-sensitive**, and the sweep and the continuation used slightly different solver paths (different root arming and cadence triggers change CVODE's step sequence — a different coin flip) | The sweep branches and run B-final now use the numerically **identical configuration** (same cadence schedule, same armed 1.5-root, same detector windows); the branch's decision then reproduces bit-for-bit — verified: branch 6 captured at t = 4.6499 Myr in both |
| DEV-4 | The plan's "P_rot/P_orb = 2/3 to 5 significant figures" was first checked on the final **instantaneous** sample — but a locked planet librates forever with tiny residual amplitude (~2×10⁻⁴ wiggle in the instantaneous ratio) | The lock statement is about the **time average**: the settled-era mean spin ratio measured 1.5000000498 (3.3×10⁻⁸ from exactly 3/2); the check now uses the mean, and reports the instantaneous wiggle alongside honestly |
| DEV-5 | The REBOUNDx cross-check first disagreed by exactly **2.005×** — the Eggleton (1998) σ-convention used by REBOUNDx defines its "constant time lag" as HALF of Hut's τ (derived by hand from σ = 4τG/(3R⁵k₂) in the circular-orbit limit) | With the mapping τ_reboundx = τ_Hut/2, the two independent formulations agree to **0.25%** — exactly the convention trap an independent cross-check exists to catch, now documented in the cross-check source |
| DEV-6 | Committing ~267 MB of regenerable CSV samples and the SQLite file would multiply the repository's tracked size for data any clone reproduces deterministically by running the notebook | The mirror pushed to the repository ships all code, documents, the executed notebook (real captured outputs), and the display page (with its embedded decimated dataset), but **not** `data/runs/` or the .sqlite3 — the notebook regenerates them bit-for-bit; recorded here and in the repository's provenance file |
| DEV-7 | The adversarial review of the BUILT code found the planned angular-momentum budget physically wrong for the locked era: holding the lock requires a **nonzero mean handle torque** ⟨T_tri⟩ = −⟨T_tidal⟩, and since the spec gives that torque no orbital back-reaction, the model itself leaks L_tot secularly in lock (~2×10⁻¹⁰ per compressed Myr — about 1.1×10⁻⁹ over the locked era, exceeding the flat 10⁻⁹ budget in exact arithmetic); the earlier PASS rested on the recorded orbit being frozen below f64 resolution at lock | The ledger check was split honestly: pre-capture era against the strict 10⁻⁹ budget (the model conserves exactly there); locked era against the model's own predicted secular leak, with the f64 orbit-quantization suppression stated in the check's own output, the notebook, and the plan documents |
| DEV-9 | The review-driven "clear stale manifests" hardening itself introduced a bug: the clearing lived in the shared directory accessor, so the sweep and the canonical continuation — which merely READ B_movie's restart file — deleted B_movie's manifest as a side effect, and the notebook's database cell then failed. The end-to-end pipeline caught it (and a shell pipe in the runner chain briefly masked the failure's exit code, costing one verification cycle) | Directory access was split into read-safe `run_dir` and owning-run `fresh_run_dir` (only the run that writes a directory clears its stale manifest); a targeted smoke test proved a sweep leaves B_movie's manifest untouched; the runner chain was rebuilt without the exit-masking pipe |
| DEV-8 | Review notes documented rather than changed (each would alter the bit-level trajectory and therefore every measured number, for no physical gain): (a) θ is re-anchored only jointly with M (3πj/2πj), so it sits at ~3.3×10⁹ rad through the resonance era — its effective absolute tolerance is rtol·θ ≈ 3×10⁻³ rad rather than the spec's 10⁻¹⁰, harmless because θ enters only 2π-periodic expressions, step size is governed by Ω's 10⁻¹⁴ tolerance, and capture phases are swept; (b) the dense-cadence γ-unwrap margin at the 1.53 trigger is 3.13 rad vs the π aliasing limit (and briefly beyond under run D's 1.535 trigger) — provably harmless as configured since aliased samples age out of the decision window and the "passed" decision never uses γ | Both documented in the code comments and here; any future change of trigger, cadence, or eccentricity must re-derive the margin and re-run the sweep |

Key measured results of the verified build: spec-literal run A after 10 Myr: spin
181.4→178.9× the orbital rate (Finding F1 confirmed); canonical capture at
4.6499 Myr movie time; sweep odds 10/64 = 15.6% ± 4.5% vs ≈7% theory; settled mean
spin ratio 1.5000000498; final periods 87.968 / 58.642 / 175.91 days (observed
87.969 / 58.646 / 175.938; instantaneous-libration wiggle ~2×10⁻⁴); angular-momentum ledger:
pre-capture drift 3.9×10⁻¹⁰ (budget 10⁻⁹), locked-era drift 9.2×10⁻¹² vs the
model's predicted 1.1×10⁻⁹ secular leak (DEV-7); libration swing decayed
6.22 → 0.15 rad across the notebook's bins (4.58 → 0.06 in the crate's
whole-locked-era check).

## Part 5 — Why no QUIT was necessary

The assignment requires quitting only if an instruction cannot be followed and the
issue cannot be resolved. Every issue above has a concrete, documented resolution that
preserves the assignment's intent: the physics model is exactly the specified one
(with the specification's own arithmetic errors corrected against its cited source);
the solver is exactly the specified one (CVODE BDF, Newton iteration, dense linear
solver, at the specified tolerances); the deviations (compression, staging, phase
sweep, mirror-then-push) are each forced by arithmetic or by the machine's real
directory/git layout, are documented here, and are submitted for approval as
decisions D1–D6 (restated in full in Part 0) before anything is built. Approval of
the master plan is approval of these resolutions.

## Part 6 — TEST 2 addendum (Jupiter + Einstein), as built

Test 2 was commissioned after test 1 shipped, with the instruction: "now run your
second test: add Jupiter and Einstein's GR correction". It adds exactly the two
pieces of physics the assignment reserved for the second test, on top of the
approved test-1 model, in a separate six-variable module (`mercury_rs/src/test2.rs`)
so that every test-1 code path, number, and verdict is untouched.

**The model additions (restated in full).**

- **Einstein's correction**: a constant-in-form apsidal drift
  dϖ/dt = 3 n G M_sun / (c² a (1 − e²)) added to a NEW sixth state variable ϖ
  (the perihelion longitude). At the spec's a₀, e₀ this is 42.98 arcseconds per
  century — the historical GR test value.
- **Jupiter's perturbations**, in the classic **Laplace–Lagrange secular** form
  (orbit-averaged; Jupiter itself on a fixed orbit, ϖ_J = 0, e_J = 0.0489,
  m_J = 1.89813×10²⁷ kg, a_J = 7.78479×10¹¹ m): de/dt += A12 e_J sin(ϖ − ϖ_J)
  and dϖ/dt += A11 + A12 (e_J/e) cos(ϖ − ϖ_J), with A11 = +n/4 (m_J/M_sun) α²
  b⁽¹⁾₃⁄₂(α) ≈ 160 ″/cy and A12 = −n/4 (m_J/M_sun) α² b⁽²⁾₃⁄₂(α); the Laplace
  coefficients are computed by a deterministic 4096-point trapezoid rule and
  unit-tested against the textbook power series to 10⁻¹⁰.
- **The resonance angle moves house**: the handle torque's argument becomes
  2(θ − f − ϖ) and the libration angle γ₂ = 2θ − 3M − 2ϖ, because the lock is to
  the ORBIT'S SHAPE and the shape now turns. Test 1's re-anchoring
  (M −= 2πj, θ −= 3πj) leaves γ₂ exactly invariant (unit-tested).

**What is deliberately NOT checked in test 2, and why.** Test 1's
angular-momentum ledger is absent: the Laplace–Lagrange terms exchange angular
momentum with Jupiter, which the model does not track, so a two-body conservation
check would be checking the wrong law. The program's output, the notebook, and
this document all state this.

**Movie-compression scope (documented deviation, same spirit as D2).** Only the
TIDAL strength is compressed 1000×; Einstein's and Jupiter's rates run at their
real values. Consequence: the braking movie contains ~12 Jupiter eccentricity
cycles instead of the thousands of the real history. The headline lock-offset
check is unaffected because it compares the measured mean ratio against the same
real ϖ̇ that acted in the integration.

**The headline acceptance target.** With the ellipse precessing at
ϖ̇ = GR + A11, a locked spin must average Ω = 1.5 n + ϖ̇, i.e. a mean spin ratio
of 1.5 + ϖ̇/n ≈ 1.5000003772 — about 3.8×10⁻⁷ ABOVE exactly 3/2. The build must
measure that offset to |Δ| ≤ 1.5×10⁻⁷ and show it exceeds 2×10⁻⁷: **the lock
follows the precessing ellipse, not the stars.**

**Verification gates (all must print SUCCESS).** Gate T2-A: GR alone reproduces
42.98 ″/cy to 10⁻³ (measured: rel. diff 7×10⁻¹³). Gate T2-B: Jupiter alone
reproduces the LL forced-eccentricity amplitude |A12/A11| e_J and period 2π/A11
to 5% (measured: 4.5438×10⁻³ vs 4.5438×10⁻³; 808.8 kyr vs 808.8 kyr). Then the
full chain: T2_movie (braking to the 1.6 restart), T2_sweep (16 phase branches),
T2_final (canonical capture + lock with all checks). Four new analytic unit tests
join the crate's suite (19 total). A second notebook
(`notebook/mercury_test2_jupiter_gr.ipynb`, authored by `build_notebook2.py`,
audited by the same structure rules, executed for real, byte-determinism proven
by double execution) stores everything in its own documented database
(`data/mercury_test2.sqlite3`, with a `run_extra` provenance table and a
`pomega_rad` sample column) and bakes its own display page
(`gui/mercury_test2.html` via `gui/bake_page2.py`).
