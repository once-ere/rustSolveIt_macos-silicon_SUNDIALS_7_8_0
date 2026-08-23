# MACOS_PORT_COMPLETE — summary of the finished port

**rustSolveIt_macos-silicon_SUNDIALS_7_8_0 is complete and verified on
macOS running on Apple Silicon.** This page summarizes what was
delivered and carries the complete set of commands to reproduce every
result — you do not need any other document to run anything named
here.

## What this project is

A pure-Rust physics simulator for macOS on Apple Silicon (arm64): a
notebook command language (`posim`) over a rigid-body/particle physics
library (`physical_object`), with quantum-mechanics and
special-function command families, a live browser scene window, 13
live GUI web pages, 13 recorded browser videos, and 109 executed
Jupyter notebooks — every example in the repository has exactly one
concomitant Jupyter notebook that defines, implements, executes and
displays it. All numerical integration is done by a **vendored
pure-Rust translation of SUNDIALS 7.8.0** (CVODE, CVODES, IDA, IDAS,
KINSOL, ARKODE) whose pure-Rust glibc-translated math library makes
the physics **byte-identical** to the Linux release's recorded
evidence. Zero `unsafe`, zero crates.io dependencies, zero warnings,
zero network access at build time.

It was produced by porting
`rustSolveIt_Using_SUNDIALS_7_8_0/version-7.8.0` (Linux) to macOS on
Apple Silicon, keeping the vendored engine byte-identical. Machine:
Apple **M5 Max**, 128 GB RAM, macOS Tahoe 26 (26.6.1). Toolchain:
rustc/cargo 1.94.0, `aarch64-apple-darwin`, linking with Apple ld from
the Xcode Command Line Tools (the project is pure Rust — no
C/C++/Fortran compiler compiled anything). Apple Silicon executes FMA
natively, which the engine's deterministic math library uses.

## The scoreboard (every number measured on this machine, 2026-08-23)

| gate | result |
|---|---|
| `cargo build --workspace --all-targets` | 0 errors, **0 warnings** |
| `cargo test --workspace` | **622 passed, 0 failed** (49 lib + 19 collision + 9 conservation + 42 DAE/equilibrium/sensitivity + 112 posim + 92 quantum + 233 special_functions + 11 vendored identities + 55 doctests) |
| 6 self-checking physics examples | 6/6 **SUCCESS** (Pluto anchor, step counts and drifts byte-identical to Linux) |
| 12 collision scripts | 12/12 exit 0; combined output **byte-identical** to the Linux evidence |
| 19 SolveIt worked-example scripts | 19/19 exit 0; 01–16 **byte-identical** to the Linux evidence |
| 59 dynamic notebooks | 59/59 exit 0; 57 **byte-identical** to Linux, 2 (the quantum pair) differ only in the last printed digit — pinned in `evidence/macos/accepted-divergences-dynamic.diff` |
| 13 recorded browser videos | 12/13 **byte-identical** to the Linux recordings; `rack_and_pinion` re-recorded here (host-libm `sin` in the GEAR/RACK residual); `record_all.py --check`: **all 13 reproduce byte for byte** |
| 109 Jupyter notebooks | executed: **109 ok, 0 failed**; checker: **109/109 pass all seven requirements** |
| 13 live GUI servers | **13/13 pass** the automated HTTP smoke test (page, live state, Start advances, Reset returns bit-exactly to t = 0) |
| index of 6,309 entities | **1177/1177** runnable posim/machine examples pass; **258/258** Rust snippets compile |
| Jupyter wire protocol + kernel | all protocol checks passed; all 7 kernel cells ok |
| engine sanity (sibling Apple-Silicon SUNDIALS port checkout) | build 0 warnings; 25/25 lib tests |
| documentation PDFs | `grammar.pdf` (66 pages, 18 worked examples) and `SolveIt.pdf` (24 pages, 16 worked examples) rebuilt from their `.tex` |

The port changed **no Rust source file**. It fixed the notebook
tooling's repository-root paths for the flattened layout (3 Python
files + 8 spec files), added two tools (`tools/macos_verify_physics.sh`,
`tools/gui_smoke.py`), re-recorded one host-libm-sensitive video, and
re-executed all 109 notebooks on this machine. Not one physics
formula, solver call, tolerance or heuristic moved — and the
byte-identity gate is what proves it.

## Complete commands

All from the repository root, in Terminal (zsh or bash).

Build and test:

```bash
cargo build --workspace --all-targets
cargo build --release --workspace --all-targets
cargo test --workspace
```

The interactive notebook (type `HELP`; `scene create` opens the live
browser window with Start/Pause/Reverse/Reset and live E, P, L
readouts):

```bash
cargo run
```

The six self-checking physics examples:

```bash
cargo run -p physical_object --release --example kepler_orbit
cargo run -p physical_object --release --example outer_solar_system
cargo run -p physical_object --release --example tumbling_body
cargo run -p physical_object --release --example charged_in_b_field
cargo run -p physical_object --release --example newtons_cradle
cargo run -p physical_object --release --example bouncing_ball_restitution
```

Any example script — the 12 collision scripts
`scripts/collisions/01_head_on_exchange.posim` …
`12_two_dumbbells.posim`, the 19 SolveIt scripts
`scripts/solveit/01_elastic_head_on.posim` … `19_hinged_door.posim`,
and the 59 dynamic notebooks `dynamic_notebooks/*.posim` — runs as
(set `POSIM_NO_BROWSER=1` first for headless):

```bash
cargo run -p posim --release -- --script scripts/collisions/01_head_on_exchange.posim
```

A dynamic notebook, interactively with its live scene window:

```bash
tools/posim_notebook kepler_orbit
tools/posim_notebook --list
```

The physics byte-identity gate against the Linux evidence:

```bash
bash tools/macos_verify_physics.sh
```

The Jupyter notebooks — regenerate, execute all 109, check all 109:

```bash
python3 notebooks/_build/regen.py
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/*.ipynb
python3 notebooks/_build/nbcheck.py notebooks/*.ipynb
```

Open the notebooks in JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/
```

The wire protocol and the JupyterLab wrapper kernel:

```bash
cargo build --release
python3 jupyter/test_protocol.py
uv venv jupyter/.venv
uv pip install -p jupyter/.venv/bin/python ipykernel jupyter_client
POSIM_NO_BROWSER=1 jupyter/.venv/bin/python jupyter/test_kernel.py
```

The live GUI pages (one server per scene, fixed ports 8895–8907; open
the printed URL in your browser), and the automated pass over all 13:

```bash
python3 gui/kepler_ellipse/server.py
python3 tools/gui_smoke.py
```

The recorded browser videos — watch, re-record, verify byte-identity:

```bash
open videos/kepler_ellipse.html
python3 recorder/src/record_all.py
python3 recorder/src/record_all.py --check
python3 recorder/tests/test_units.py
python3 recorder/tests/test_end_to_end.py
```

Record a new video from any scene script:

```bash
cargo build --release -p posim
python3 recorder/src/record_video.py videos/scenes/kepler_ellipse.posim \
     -o mine.html --frames 360 --dt 0.02 --title "Kepler orbit, e = 0.6"
```

The entity index and its verifiers:

```bash
open index_of_entities.html
python3 tools/verify_index_examples.py
python3 tools/verify_tierb_examples.py
```

The documentation PDFs (each twice — the table of contents needs the
second pass):

```bash
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode SolveIt.tex
pdflatex -interaction=nonstopmode SolveIt.tex
```

## Where the deliverables live

- `grammar.md` / `grammar.tex` / `grammar.pdf` — the complete command
  language: lexer, full EBNF, type system, every command, the stack
  machine, the SUNDIALS 7.8.0 engine, the browser videos, and 18 fully
  worked examples.
- `SolveIt.md` / `SolveIt.tex` / `SolveIt.pdf` — the complete solution
  guide for a first-time reader, with 16 fully documented worked
  examples (scripts in `scripts/solveit/`).
- `ARCHITECTURE.md` — module responsibilities and pinned cross-module
  contracts. `CLAUDE.md` — working rules for contributors and agents,
  including the macOS platform notes.
- `PORT_MACOS_PROVENANCE.md` — this port's full provenance: machine,
  toolchain, the complete change list, and every verification command
  with its measured result. `VERIFICATION_MACOS.md` — the per-suite
  verification matrix and every documented divergence.
- `REBOUND_PROVENANCE_*.md` — one self-contained provenance page per
  rebound-physics (restitution) notebook, each carrying the complete
  commands to test, verify, execute and display it.
- `notebooks/` — the 109 executed Jupyter notebooks (one per example);
  `gui/` — the 13 live GUI pages; `videos/` — the 13 recorded browser
  videos; `evidence/macos/` — the macOS gate logs and the pinned
  divergence diffs.
