# Planet Mercury Tidal-Locking Plan — Document 0 of 8: The Master Plan

**Project:** `planet_Mercury` — simulate, store, display, and explore how the planet
Mercury became tidally locked to the Sun in a 3:2 spin-orbit resonance, using the
pure-Rust rustSolveIt / SUNDIALS 7.8.0 engine as the only compute engine, with the
result delivered as a Jupyter notebook plus a browser display page.
**Audience:** written for a reader with U.S. high-school math and science. Everything
needed to understand this document is inside this document.
**Status:** PLAN (awaiting approval). Nothing has been built yet. After approval, the
build is carried out by an orchestrated multi-agent ("ultracode") workflow.
**Quality note:** this plan was adversarially audited before being presented — four
independent reviewer agents recomputed every physics number, checked every claim about
the engine against its actual source code, and hunted contradictions; all confirmed
findings are already fixed in what you are reading.

Words used throughout this plan, expanded once here and re-explained in each document
where they matter: **ODE** = ordinary differential equation (a rule for how fast
something changes); **BDF** = backward differentiation formulas (the solver's method
for stiff problems — problems mixing fast wiggles with very slow drifts); **CSV** =
comma-separated values (a data table saved as plain text); **CLI** = command-line
interface (running a program by typing its name and options in a terminal); **GUI** =
graphical user interface; **HUD** = heads-up display (the on-screen number readouts on
the display page); **SQL** = Structured Query Language (the standard mini-language for
asking a database questions).

---

## 1. The problem, in three sentences

Mercury spins on its axis exactly three times for every two trips around the Sun — a
"3:2 spin-orbit resonance," which is Mercury's strange form of tidal locking (our Moon,
by contrast, is locked 1:1 and always shows Earth one face). This happened because the
Sun's tides slowly braked Mercury's once-fast spin over billions of years of its
4.5-billion-year existence, until the Sun's gravitational grip on Mercury's slightly
lopsided shape snapped the spin into step as it fell through the 3:2 ratio. This
project simulates that entire history as a two-body problem — the Sun as a point mass,
Mercury as an extended, deformable, slightly lopsided body — by integrating five
coupled ordinary differential equations (ODEs) with the SUNDIALS CVODE solver (BDF
method, Newton iteration, dense linear solver), then stores every result in a
documented database and displays the whole story interactively in a web browser.

Two things are deliberately **excluded** in this first test, per the assignment: the
gravitational tugs of Jupiter and the other planets, and Einstein's general-relativity
perihelion correction. Both are reserved for a planned second test.

---

## 2. What will exist when the build is done (the deliverable inventory)

Everything below is created inside the working folder
`/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/planet_Mercury/`:

```
planet_Mercury/
├── plan/                            ← these eight planning documents (already written)
│   ├── 00_PLAN_OVERVIEW.md          (this file — the master plan)
│   ├── 01_PHYSICS_AND_MATH.md       (every equation, constant, and timescale)
│   ├── 02_ARCHITECTURE_AND_ENGINE.md(the Rust crate, solver calls, staging, GUI page)
│   ├── 03_NOTEBOOK_INSTRUCTIONS.md  (teacher/student instructions, complete)
│   ├── 04_DATABASE_SCHEMA.md        (every table, column, and query, documented)
│   ├── 05_VERIFICATION_PLAN.md      (every check that must pass, with commands)
│   ├── 06_BUILD_ORCHESTRATION.md    (how the multi-agent build will run)
│   └── 07_PROVENANCE_AND_DEVIATIONS.md (every issue found + how each is resolved)
├── SOURCE_SPECIFICATION.md          ← a verbatim copy of the consolidated physics
│                                      specification this project implements, so the
│                                      published repository carries its own source
├── mercury_rs/                      ← the pure-Rust compute crate (zero unsafe,
│   │                                  zero external dependencies, zero warnings)
│   ├── Cargo.toml                   (standalone; path-depends only on the engine's
│   │                                  vendored pure-Rust SUNDIALS crates)
│   ├── src/                         (params, Kepler solver, the 5-ODE right-hand
│   │                                  side, the CVODE driver, CSV output)
│   └── tests/                       (analytic-expectation unit tests)
├── mercury_crosscheck/              ← OPTIONAL small GPL-licensed crate that uses the
│                                      workspace's pure-Rust REBOUNDx port (tides_spin)
│                                      to independently check the tidal braking rate
├── notebook/
│   ├── mercury_tidal_locking.ipynb  ← THE deliverable Jupyter notebook (executed,
│   │                                  with real captured outputs)
│   ├── run_notebook.py              (stdlib-only batch executor / re-verifier)
│   └── check_notebook.py            (stdlib-only structure auditor)
├── data/
│   ├── runs/                        (CSV files written by mercury_rs, one dir per run)
│   └── mercury_orbit.sqlite3        (the database, built by the notebook from the CSVs)
├── gui/
│   ├── bake_page.py                 (stdlib-only Python that embeds the decimated —
│   │                                  i.e. thinned-down — dataset into a
│   │                                  self-contained HTML player)
│   └── mercury_orbit.html           (the browser display page — orbit animation, spin
│                                      arrow, spin/orbit-ratio dial, libration plot,
│                                      period & angular-momentum history plots)
└── INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md
                                     ← the complete authoring + running instructions
                                       (a full standalone copy, as required)
```

In addition, after everything above is verified, a **mirror copy** of the finished
sub-project is committed into the engine's git repository and pushed to GitHub
(Section 5 explains exactly how and why).

---

## 3. The scientific plan in brief (fully worked in plan document 1)

The state of the system is five numbers: orbit size a, orbit ovalness e, orbit clock
angle M, spin angle θ, and spin rate Ω. Five coupled ODEs evolve them: Kepler's law
drives M; the Hut (1981) "constant time lag" tidal torque brakes Ω and gently reshapes
a and e; the "handle" torque on Mercury's lopsided figure (the triaxial torque) is what
captures the spin at resonance. All integration is SUNDIALS CVODE: BDF method, Newton
iteration, dense linear solver, relative tolerance 1.0×10⁻¹², absolute tolerances
[1.0×10⁻³, 1.0×10⁻⁶, 1.0×10⁻¹⁰, 1.0×10⁻¹⁰, 1.0×10⁻¹⁴] for [a, e, M, θ, Ω], maximum
step 864,000 s.

Planning-stage arithmetic (worked in full in plan document 1, and re-verified by an
independent audit agent) produced two findings that shape the run plan:

- **Finding F1:** with the specification's own tidal-strength constants, the braking
  from the fast newborn spin down to the 3:2 mark takes about **4.7 billion years**
  (the spin decays exponentially toward its tides-only resting rate of 1.256× the
  orbital rate, which takes 6.6 shrink-steps of factor e ≈ 2.718 each ≈ 710 million
  years long). The spec's 10-million-year window is therefore ~470× too short to reach
  capture. Fix: a documented **time-compression factor S = 1000** on the tidal
  strength (the spec itself says dissipation is "deliberately set strong" to make the
  story tractable; it just was not set strong enough). At S = 1000 the capture lands
  at ≈ 4.7 million years — comfortably inside the window, with >5 million years left
  to display the locked state.
- **Finding F2:** the handle torque wiggles every ~6 hours while Mercury spins fast,
  which would force ~10¹¹ solver steps (hundreds of times the spec's own step budget).
  Fix: **staged integration** — far from resonance the wiggle averages to zero and is
  switched off (the solver then takes giant steps); near the resonances the full
  five-equation system runs. Both stages use CVODE; a validation run proves the two
  models agree where they hand over.

And one piece of real science the plan turns into a feature: capture into 3:2 at
Mercury's eccentricity is **probabilistic** (≈ 7.0% per crossing — Goldreich & Peale
1966; re-derived independently during the audit). A 64-branch sweep of the spin phase
at the crossing measures the capture odds, and the first captured branch becomes the
canonical "history of Mercury" trajectory.

The six runs: **A** spec-literal (documents F1 honestly), **B** the time-compressed
braking movie up to a saved restart point just above the 3:2 crossing, **C** the
64-branch phase sweep from that restart, **B-final** the first captured branch
continued to the end (the canonical locked history), **D** a guaranteed-capture
high-eccentricity encore, **E** the staging-seam validation.

---

## 4. How each piece of the assignment is satisfied

Symbols used in this table: G (Newton's gravitational constant), M☉ (the Sun's mass),
m (Mercury's mass), R (Mercury's radius), a (orbit size), e (orbit ovalness),
n (the orbit's average angular speed, n = √(G·(M☉+m)/a³)), Ω (spin rate),
C (Mercury's spin moment of inertia — its resistance to changes in spin,
C = 0.34·m·R²), θ (spin angle), M (orbit clock angle).

| Assignment requirement | How the plan satisfies it |
|---|---|
| "history of, and present day, orbit of planet Mercury" | Runs B + B-final evolve a, e, M, θ, Ω from Mercury's fast-spinning youth to the locked present; the final state reproduces today's measured periods (year 87.969 d, day 58.646 d, solar day 175.938 d) |
| "clearly displays and explores tidal locking" | The browser page animates the orbit + a spin arrow with a spin/orbit-ratio dial that settles on exactly 1.5; a libration plot shows the resonance angle γ = 2θ − 3M rocking after capture; the notebook explores the capture odds and the eccentricity effect |
| "changes in time of the orbital period" | P_orb(t) = 2π/n(t) is stored for every sample and plotted |
| "changes in time of the angular momenta" | Spin angular momentum C·Ω, orbital angular momentum ≈ m·n·a²·√(1−e²), and their sum (the conservation audit) are stored for every sample and plotted |
| "over the entire (time) existence of planet Mercury" | The time-compressed run covers the full braking story (a faithful, documented, 1000×-sped-up movie of the multi-billion-year history) |
| "two-body Sun-Mercury; Sun = point mass; Mercury = extended, deformable, nearly spherical" | Exactly the model: point-mass Sun; Mercury with radius, moment of inertia C = 0.34·m·R², triaxial asymmetry (B−A)/C = 10⁻⁴ (B−A measures how lopsided Mercury's equator is), Love number k₂ (squishiness), tidal lag τ (sluggishness) |
| "do NOT include Jupiter … or GR" | Excluded; explicitly reserved for the second test |
| "rustSolveIt … as the sole compute engine and display engine" | All integration is the engine's vendored pure-Rust SUNDIALS CVODE; the display page follows the engine's established self-contained browser-player pattern exactly |
| "pure rust" code | The compute crate is pure Rust: `#![forbid(unsafe_code)]`, zero external dependencies, zero warnings (the notebook/GUI tooling is Python-stdlib by the engine's own convention — Python tooling is explicitly outside the Rust constraints in that project, and the Jupyter notebook format itself requires it) |
| "correct, legal, rust SolveIt Jupyter Notebook" | The notebook follows the engine project's exact notebook conventions: Python 3 (ipykernel), standard-library-only code cells, every code cell preceded by an explaining markdown cell, no cross-references to other documents, real captured outputs, a final interactive save cell |
| "complete, correct, precise and verified set of clear instructions for teachers and students" | `INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md` (plan document 3 is its full draft) — self-contained, step-by-step, verified by actually following it during the build |
| "completely and clearly documented database schema" | An SQLite database with five documented tables; plan document 4 contains the complete schema with every column explained and worked example queries |
| "check and verify" everything | Plan document 5 defines 40+ concrete checks: unit tests with analytic expectations, notebook re-execution, byte-determinism of the baked page, physics acceptance targets, and the engine's own untouched-ness gates |
| "push the repo" to the GitHub URL | Section 5 below — mirror into the engine repo, explicit-path staging, push to `origin` (which is exactly the required URL); the pushed tree includes a verbatim copy of the source physics specification |

---

## 5. Where the work lives, and how the push works (Decision D1 — please approve)

Plain-words git glossary for this section: a **repository** ("repo") is a folder whose
history git tracks; the **work tree** is that folder's files; **staging** means telling
git exactly which files go into the next saved snapshot (commit); **origin** is the
nickname git gives the online copy it pushes to.

Two facts discovered during planning make this the one place where the assignment's
letter needs a small, explicit resolution:

1. The required build location
   `/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/planet_Mercury/`
   sits in the **outer workspace folder**, which is **not a git repository** (a
   standard git command confirms that folder is not itself a git project — verified).
2. The required push target
   `https://github.com/once-ere/rustSolveIt_macos-silicon_SUNDIALS_7_8_0.git` is
   already the `origin` remote of the **inner engine repository**
   `/Users/nsh/Developer/github/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/`
   (branch `main`, clean, in sync with origin — verified during planning).

**Default resolution (D1):** build and verify everything at the required outer
location, then **mirror** the finished sub-project into the inner repository as a
top-level `planet_Mercury/` directory — exactly the way that repository already hosts
its `rebound_rust/` and `reboundx_rust/` sub-projects: committed in full, but kept out
of the root cargo workspace by adding `"planet_Mercury"` to the one `exclude = [...]`
line in the root `Cargo.toml` (a backup of that file is made first). The only file that
differs between the outer original and the inner mirror is one relative path in
`mercury_rs/Cargo.toml` (the vendored SUNDIALS crates sit at a different relative depth
from each location); the mirror is re-built and re-tested in place to prove it. Then a
root-level `PLANET_MERCURY_PROVENANCE.md` is added (that repository's naming convention
for documenting added work), changes are staged **by explicit path only** (never
`git add -A` — that repository's own history documents why), committed, and pushed.

Safety facts, verified during planning: the outer `p*.txt` files physically cannot be
pushed (they live outside the repository work tree; git refuses to stage anything
outside it); no secrets are involved anywhere; per the assignment's rule, **no
`.gitignore` files are created or committed** in `planet_Mercury/` (the repository's
existing ignore rule for any `target/` build folder already keeps build output out of
git).

---

## 6. Architecture in one picture

```
                     ┌───────────────────────────────────────────────┐
                     │  mercury_rs  (pure Rust, zero deps, no unsafe)│
 constants, run      │  ┌─────────┐   ┌──────────────────────────┐  │
 selection (CLI) ───▶│  │ params  │──▶│ 5-ODE right-hand side    │  │
                     │  └─────────┘   │ (Kepler solve, tidal     │  │
                     │                │  brake, handle torque)   │  │
                     │                └───────────┬──────────────┘  │
                     │                            ▼                 │
                     │        vendored pure-Rust SUNDIALS 7.8.0     │
                     │        CVODE: BDF + Newton + dense solver    │
                     │                            │                 │
                     │                            ▼                 │
                     │        CSV sample files + run manifest       │
                     └────────────────────────────┬──────────────────┘
                                                  ▼
      ┌──────────────────────────────────────────────────────────────┐
      │  mercury_tidal_locking.ipynb  (Python 3, stdlib only)        │
      │  • explains the physics (self-contained prose)               │
      │  • builds & runs mercury_rs via subprocess (cargo run)       │
      │  • ingests CSVs into data/mercury_orbit.sqlite3 (sqlite3)    │
      │  • queries the database back; prints numeric tables         │
      │  • asserts every acceptance target (capture, periods, L)    │
      │  • bakes gui/mercury_orbit.html and opens it in the browser │
      └──────────────────────────────┬───────────────────────────────┘
                                     ▼
      ┌──────────────────────────────────────────────────────────────┐
      │  mercury_orbit.html  (one self-contained file, vanilla JS)   │
      │  orbit animation • spin arrow • ratio dial → 1.5             │
      │  libration plot • P_orb & angular-momentum history plots     │
      │  play/pause • scrub • speed • HUD readouts                   │
      └──────────────────────────────────────────────────────────────┘
```

---

## 7. The build, in phases (fully detailed in plan document 6)

1. **Phase 0 — baseline:** prove the engine's own gates are green before touching
   anything (workspace build + 622 tests + physics byte-identity script).
2. **Phase 1 — mercury_rs crate:** implement params → Kepler → RHS → CVODE driver →
   CSV output; unit tests with analytic expectations at every step; adversarial code
   review before the phase closes.
3. **Phase 2 — science runs:** execute runs A, B, C, B-final, D, E in that order;
   select the canonical captured branch; write all CSVs and the run manifests.
   (Optional Phase 2b: the REBOUNDx-port cross-check of the braking rate.)
4. **Phase 3 — database:** build the SQLite database from the CSVs; run the
   documented example queries; verify row counts and invariants.
5. **Phase 4 — display page:** bake the self-contained HTML player; verify
   byte-determinism and zero external fetches; visual checks in a real browser.
6. **Phase 5 — notebook:** author + execute the Jupyter notebook end-to-end
   (headless); audit its structure; capture real outputs.
7. **Phase 6 — instructions & docs:** finalize the teacher/student instructions
   (verified by following them literally in a clean shell), copy the source
   specification into the project as `SOURCE_SPECIFICATION.md`, and complete the
   provenance file.
8. **Phase 7 — mirror, commit, push:** Section 5's procedure, with the engine's gates
   re-run after the root `Cargo.toml` edit, then push and post-push verification.

Each phase ends with its own verification gate (plan document 5); the multi-agent
workflow runs adversarial reviewers over the physics code and the notebook before
anything is committed.

---

## 8. Decisions submitted for approval

| # | Decision | Default in this plan |
|---|---|---|
| D1 | Where the pushed copy lives | Build at the required outer path; mirror into the engine repo top level (rebound_rust precedent); push to origin = the required URL |
| D2 | Spec's 10-Myr window vs. its ~4.7-Gyr braking timescale (Finding F1) | Time-compression factor **S = 1000** on tidal strength k₂τ (capture at ≈ 4.7 Myr, >5 Myr of locked display), fully documented; spec-literal run A kept alongside for honesty |
| D3 | Spec's step budget vs. the oscillating torque (Finding F2) | Staged integration (handle torque off far from resonance), seam-validated; all stages CVODE BDF at spec tolerances |
| D4 | Capture is probabilistic (≈ 7%) | 64-branch spin-phase sweep (with one finer-grid contingency re-sweep if no branch captures — about a 1% possibility); first captured branch = canonical history; measured odds displayed beside theory with honest error bars |
| D5 | Notebook & display technology | Python 3 (ipykernel) stdlib-only notebook driving the Rust crate by subprocess (the engine project's exact convention); display = one self-contained HTML player page (the engine's videos-page pattern) |
| D6 | ".gitignore files must not be pushed" | No `.gitignore` created or committed anywhere in planet_Mercury/ |

**To approve the plan as-is, reply "approved."** To change any decision, name the
number (for example: "approved, but D4: use 32 branches"). After approval, the
orchestrated build begins immediately and runs through Phase 7 (including the final
push) without further pauses, reporting progress as it goes.
