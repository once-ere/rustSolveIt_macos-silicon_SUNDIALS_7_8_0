# Planet Mercury Tidal-Locking Plan — Document 6 of 8: How the Build Will Be Orchestrated

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how
Mercury became locked to the Sun in a 3:2 spin-orbit resonance.
**Audience:** written for a reader with U.S. high-school math and science. Everything
needed is inside this document.
**Status:** PLAN (awaiting approval). This document describes how the approved plan
will be EXECUTED by an orchestrated team of AI build agents ("ultracode" mode), with
the human-readable checkpoints along the way.

Words expanded once for this document: **CSV** = comma-separated values (a data table
saved as plain text); **SQL** = Structured Query Language (database queries);
**referential integrity** = every stored row points at a run that really exists.

---

## 1. Ground rules every build agent obeys

1. **Backups before modification:** before any pre-existing file is edited, a copy
   goes into the engine's `.backups/<date>/` folder (which is never committed). The
   only pre-existing file this project modifies at all is the engine repository's
   root `Cargo.toml` (one line appended to its `exclude` list at mirror time).
2. **Logs before edits:** every build/test/run command is executed as
   `<command> 2>&1 | tee <logfile>`, and the log is read before the next action. No
   command that produced no visible output is ever blindly re-run.
3. **Two strikes, then rethink:** at most 2 retries of any failing command before the
   strategy is changed (never grind the same failure).
4. **Read-only zones:** the engine's `sundials_rs/` (vendored solver), the reference
   trees, and everything not listed as a deliverable are never modified. A missing
   solver symbol stops the build with a report naming the exact symbol — it is never
   stubbed or reimplemented.
5. **The run outputs, the database, the page, and the notebook are produced by the
   code being shipped** — no hand-edited numbers anywhere. Every number quoted in
   documentation is measured from an actual logged run.
6. **Verification gates between phases.** The full gate definitions (all commands and
   expected values) are plan document 5; for self-containedness each phase below also
   states its closing check inline, in one line.

## 2. The phases, their agents, and their checkpoints

**Phase 0 — Baseline (1 agent).** Re-run the engine's own health gates. Nothing of
ours is built until the foundation is proven green.
*Closes when:* `cargo build --workspace --all-targets` is warning-free,
`cargo test --workspace` shows 622 passing tests, `bash tools/macos_verify_physics.sh`
exits 0, and `git status -sb` shows a clean tree in sync with origin (Gate 0).

**Phase 1 — The compute crate (parallel build + adversarial verify).**
- Builder agents implement `mercury_rs` module by module in this order:
  `params` → `kepler` → `hut` → `rhs` → `output` → `driver` → `main`, each module
  landing together with its unit tests (the fifteen named analytic tests of Gate 1 —
  Kepler inversion, Hut-polynomial values, torque signs, the linearized despin rate,
  exact angular-momentum closure, re-anchoring exactness, bit-identical restarts,
  libration frequency, root-finding precision, and format round-trips).
- After the crate is green, independent reviewer agents attack it in parallel:
  a physics reviewer (formulas vs. the corrected Hut/Goldreich–Peale forms), a Rust
  reviewer (borrow-guard discipline, no unwraps on solver options, zero unsafe, no
  external dependencies), and a determinism reviewer (no timestamps, no
  platform-dependent output). Every confirmed finding is fixed and re-tested.
*Closes when:* `cargo build --release --all-targets` is warning-free, `cargo test`
is fully green including all fifteen named tests, and every confirmed review finding
is fixed (Gate 1 + Gate 6 crate portion).

**Phase 2 — The science runs (sequential, data-dependent).** Runs execute in
dependency order: **A** (spec-literal honesty run), **B** (movie to the restart point
saved at spin ratio 1.6), **C** (the 64-branch sweep — branches run in parallel),
**B-final** (first captured branch continued to the end), **D** (high-eccentricity
encore), **E** (seam validation). After each run, its self-check verdict and its log
are read before the next run starts. Contingency, aligned with the verification
plan's gate 2.3: if the sweep produces zero captures (about a 1% possibility at ≈ 7%
odds per branch), exactly one finer-grid re-sweep is run (stored separately as
`C_sweep2`), and the outcome either way is recorded in the provenance file.
*Closes when:* all six runs print SUCCESS; run B-final shows a capture event at
≈ 4.7 Myr (±20%) with final periods within 0.1% of 87.969 / 58.646 / 175.938 Earth
days; the sweep has at least one capture after at most one re-sweep (Gate 2).

**Phase 2b — Optional cross-check (1 agent, may run in parallel with Phase 3).**
Build the small GPL-licensed `mercury_crosscheck` crate against the workspace's
pure-Rust REBOUNDx port, configure its `tides_spin` effect with the same Sun–Mercury
parameters, and compare its tidal despin rate over a short window against
`mercury_rs`'s secular rate. If time or an unexpected obstacle prevents this, the
crate is dropped and the drop is recorded in the provenance file — it is a
nice-to-have independent check, not a required deliverable.
*Closes when:* the two despin rates agree within the documented tolerance, or the
drop is recorded (Gate 2.8).

**Phase 3 — The database (1 agent).** Write the ingest logic, build
`mercury_orbit.sqlite3`, run the seven documented queries, check row counts,
referential integrity, and the angular-momentum ledger in SQL.
*Closes when:* row counts match the CSV line counts, zero orphan rows, every run row
says SUCCESS, and the ledger stays within its 1×10⁻⁹ budget (Gate 3).

**Phase 4 — The display page (build + browser verify).** Implement `bake_page.py`,
bake `mercury_orbit.html`, prove byte-determinism (bake, copy aside, bake again,
compare with `cmp`), scan for self-containedness (no network use of any kind), then
drive a real browser over the page.
*Closes when:* the double-bake compares byte-identical, the page makes zero network
requests and zero console errors, playback advances, the capture moment shows the
ratio dial at 1.5, and the libration plot rocks and decays (Gate 4).

**Phase 5 — The notebook (author + execute + audit).** Author
`mercury_tidal_locking.ipynb` to the structure fixed in plan document 3 — in brief:
title/abstract; how-to-run; glossary; the physical situation with every equation and
both documented deviations; the solver contract; the subprocess driver; one explained
cell per run (config echo, A, B, C, B-final, D, E); database build; the seven
retrieval queries; the verification gauntlet; page bake + browser open; conclusions;
interactive save cell. Implement `run_notebook.py` and `check_notebook.py`; execute
the notebook headless end-to-end so every committed output is real; audit the
structure; re-execute to prove byte-determinism.
*Closes when:* the headless run executes every cell clean, the structure audit
passes, the gauntlet cell prints all PASS, and a second run is byte-identical
(Gate 5).

**Phase 6 — Documentation freeze (writer + fresh-eyes agent).** Finalize
`INSTRUCTIONS_FOR_TEACHERS_AND_STUDENTS.md` with the measured numbers filled in; copy
the consolidated physics specification into the project as
`SOURCE_SPECIFICATION.md` (so the pushed repository carries its own source); complete
the provenance file with every finding and its resolution. A fresh-eyes agent (no
build context) then follows the instructions literally, from a clean shell, and must
reach the end without improvising — any gap found is a documentation bug, fixed and
re-walked.
*Closes when:* the fresh-eyes walk-through (check 5.5 of the verification plan)
succeeds with zero improvisation and the docs reviewer confirms every quoted number
is measured (Gate 5.5 + Gate 6 docs portion).

**Phase 7 — Mirror, commit, push (1 careful agent, serialized).** Exactly the
procedure fixed in plan document 2: back up the engine root `Cargo.toml` into
`.backups/<date>/`; copy the verified tree into the engine repository as top-level
`planet_Mercury/`; apply the one-line dependency-path change in the mirror's
`Cargo.toml`(s); append `"planet_Mercury"` to the engine root `Cargo.toml` workspace
`exclude` list; re-build/re-test the mirror in place; re-run the engine's own gates;
write `PLANET_MERCURY_PROVENANCE.md` at the engine-repo root; stage by explicit path
only; inspect the staged list; commit; push to `origin` (which is the required URL
`https://github.com/once-ere/rustSolveIt_macos-silicon_SUNDIALS_7_8_0.git`); verify
the remote hash; clone-test the pushed repository in a scratch directory.
*Closes when:* the staged list contains only the intended paths (and zero
`.gitignore` files under `planet_Mercury/`), the engine's 622 tests and physics gate
are still green, the push lands, the remote hash matches local HEAD, and the fresh
clone builds and tests green (Gate 7).

## 3. Commit discipline

Commits happen only in Phase 7 (the outer build folder is not itself a git
repository; the engine repository receives the work as one reviewed, verified
import — the same way it received its rebound_rust and reboundx_rust sub-projects,
each as a single clean import commit with a root-level provenance document). If
review findings after the first push ever require fixes, each fix lands as its own
small, described commit.

## 4. Estimated effort and runtime

- Compute: run B-final's resonant stage is the big one — order 10⁸ CVODE steps,
  minutes to tens of minutes in release mode on this machine; the 64 sweep branches
  are short and run in parallel. Total simulation wall-time budget: under two hours.
- The full orchestrated build (all phases, all reviews, all gates): a working
  session on the order of hours, reported live as it proceeds.

## 5. What the human sees along the way

At each phase boundary: a short progress report naming the gate that closed and any
finding that changed anything. At the end: a final report with (1) every gate's
result, (2) the measured science numbers (capture time, capture fraction with error
bars, final periods, ledger drift), (3) the commit hash pushed, and (4) the complete
file inventory of what was delivered.
