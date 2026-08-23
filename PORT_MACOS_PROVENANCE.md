# PORT_MACOS_PROVENANCE — the macOS / Apple Silicon port, in full

This document records how **rustSolveIt_macos-silicon_SUNDIALS_7_8_0**
was produced, on what machine, what changed relative to the donor, and
every verification command with its measured result. It is
self-contained: every command needed to reproduce every claim on this
page appears **on this page**, run from the repository root.

## 1. What this project is

A pure-Rust physics simulator for macOS on Apple Silicon (arm64): a
notebook command language (`posim`) over a rigid-body/particle physics
library (`physical_object`), with quantum-mechanics and
special-function command families, a live browser scene window, 13
live GUI web pages, 13 recorded browser videos, and 109 executed
Jupyter notebooks — every example in the repository has exactly one
concomitant Jupyter notebook that defines, implements, executes and
displays it. All numerical integration is done by a **vendored
pure-Rust translation of SUNDIALS 7.8.0** (CVODE, CVODES, IDA, IDAS,
KINSOL, ARKODE). Zero `unsafe`, zero crates.io dependencies, zero
warnings, zero network access at build time.

## 2. Lineage

It was produced by porting
`once-ere/rustSolveIt_Using_SUNDIALS_7_8_0/version-7.8.0` (the Linux
release) to macOS on Apple Silicon. The layout is flattened: this
repository's root **is** the workspace (the donor nested it under
`version-7.8.0/`).

The vendored engine `sundials_rs/` was carried over **byte-identical**
to the donor's, which is itself byte-identical to
`once-ere/SUNDIALS_7_8_Rust_port_for_Linux@780b916` (2,929 files).
That engine's own pure-Rust glibc-translated math library
(`sundials_libm`, used by default; the `host-libm` cargo feature that
bypasses it is **off**) makes the solver numerics host-independent —
which is why this port needed **no engine changes** and reproduces the
Linux physics byte for byte (section 6).

A second engine lineage exists on this machine:
`SUNDIALS_7_8_Rust_port_for_AppleSilicon_macos`, an independent
Apple-Silicon-verified port of the same SUNDIALS 7.8.0 release (199
example variants verified against upstream C on this platform). It was
used as platform evidence that pure-Rust SUNDIALS runs correctly on
Apple Silicon, and its build + unit-test gate was re-run first-hand
during this port (7 crates, zero warnings, 25/25 lib tests). It was
**not** vendored here, because it predates the `sundials_libm`
host-independent math library that the byte-identity results below
depend on, and the simulator was written against the later engine's
API.

## 3. Machine and toolchain (as measured, 2026-08-23)

| item | value |
|---|---|
| machine | Apple **M5 Max**, 128 GB RAM |
| OS | macOS Tahoe 26 (`sw_vers`: macOS 26.6.1, build 25G76) |
| architecture | arm64 (`uname -m`) |
| Rust | rustc 1.94.0 / cargo 1.94.0, `aarch64-apple-darwin` |
| Python | 3.14.7 (Homebrew), used for tooling only — no Rust dependency |
| LaTeX | TeX Live 2026 `pdflatex` |
| linker | Apple ld from the Xcode Command Line Tools |

No C, C++ or Fortran compiler compiled anything: the project is pure
Rust. Apple Silicon executes FMA natively, which the vendored engine's
deterministic `pow` requires.

## 4. The complete change list

The port changed **no Rust source file**. Every change is tooling,
documentation, or platform evidence:

1. `notebooks/_build/nbcheck.py` — the repository root is
   `parents[2]` of the checker, not `parents[3]` (the donor had one
   extra `version-7.8.0/` nesting level).
2. `notebooks/_build/regen.py` — the `pairs_with` metadata written
   into each notebook is now repo-root-relative (the donor hardcoded a
   `version-7.8.0/` prefix).
3. `notebooks/_build/nbtext.py` — the launch instructions embedded in
   every notebook clone **this** repository and build from the
   repository root.
4. `notebooks/_build/specs/rust_*.json` (6 files),
   `specs/video_piston_crankshaft.json`,
   `specs/video_piston_crankshaft.handwritten.py` — the hand-maintained
   spec paths lose their `version-7.8.0/` prefix.
5. `tools/macos_verify_physics.sh` — **new**: the macOS byte-identity
   gate (section 6).
6. `tools/gui_smoke.py` — **new**: the automated HTTP smoke test over
   all 13 live GUI servers (portable stdlib Python).
7. `videos/rack_and_pinion.html` — re-recorded on this machine
   (section 6.5 explains why this one video is host-sensitive).
8. All 109 notebooks under `notebooks/` — regenerated and re-executed
   on this machine, so their recorded outputs are this machine's.
9. Documentation — `README.md`, `CLAUDE.md`, `ARCHITECTURE.md`,
   `SolveIt.md`/`.tex`, `grammar.md`/`.tex` re-scoped for this
   platform; PDFs rebuilt; this file, `VERIFICATION_MACOS.md`,
   `MACOS_PORT_COMPLETE.md` and the `REBOUND_PROVENANCE_*.md` files
   added; macOS gate evidence added under `evidence/macos/`.

Not one physics formula, solver call, tolerance or heuristic moved —
and the byte-identity gate is what proves it.

## 5. Build and test (commands + measured results)

```bash
cargo build --workspace --all-targets
```
Result: all 12 crates (4 workspace + 7 engine + vendored `spec_math`)
compile with **zero warnings**.

```bash
cargo build --release --workspace --all-targets
```
Result: zero warnings.

```bash
cargo test --workspace
```
Result: **622 passed, 0 failed** — 49 physical_object lib +
19 collision + 9 conservation + 42 constrained/DAE/equilibrium/
sensitivity + 112 posim + 92 quantum + 233 special_functions +
11 vendored identities + 55 doctests.

The engine's own gate, run in a sibling checkout of
`SUNDIALS_7_8_Rust_port_for_AppleSilicon_macos`:

```bash
cargo build --workspace          # zero warnings
cargo test --workspace --lib     # 25 passed, 0 failed
```

## 6. Physics byte-identity against the Linux evidence

The gate re-runs, on this machine, exactly the three suites whose
Linux outputs are recorded in `evidence/port-7.8.0/`, and diffs the
results byte for byte (scene ports and the current directory are
normalised, as in the originals):

```bash
cargo build --release -p posim
bash tools/macos_verify_physics.sh
```

Measured result (exit 0):

```
IDENTICAL-MODULO-DOCUMENTED  examples-macos.log vs examples-7.8.0.log (pinned: accepted-divergences-examples.diff)
IDENTICAL  collision-scripts-macos.log == collision-scripts-7.8.0.log
IDENTICAL-MODULO-DOCUMENTED  dynamic-notebooks-macos.log vs dynamic-notebooks-7.8.0.log (pinned: accepted-divergences-dynamic.diff)
```

### 6.1 The six self-checking examples

All six print SUCCESS and exit 0; output is **byte-identical** to the
Linux evidence except one pinned blank line:

```bash
cargo run -p physical_object --release --example kepler_orbit
cargo run -p physical_object --release --example outer_solar_system
cargo run -p physical_object --release --example tumbling_body
cargo run -p physical_object --release --example charged_in_b_field
cargo run -p physical_object --release --example newtons_cradle
cargo run -p physical_object --release --example bouncing_ball_restitution
```

Anchors, identical to Linux to the last digit: Pluto
`x = 31.78592516  y = 38.63618957  z = 3.19279415` at t = 500,000
days, energy drift `7.835809e-07`, `12581` internal steps; Kepler
`|dA|/|A| = 1.131858e-07`; tumbling body `|dL|/|L| = 0.000000e+00`
exactly.

**The pinned examples divergence** (2 diff lines): one inserted blank
line after the `outer_solar_system` header. The donor added a leading
`"\n"` to that example's first `println!` *after* its evidence log was
recorded — the donor's own current source cannot reproduce its
evidence on Linux either. Every number is identical. (The Windows port
pinned the same divergence.)

### 6.2 The twelve collision scripts

```bash
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/collisions/01_head_on_exchange.posim
# … through …
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/collisions/12_two_dumbbells.posim
```

Result: 12/12 exit 0; concatenated output **byte-identical** to the
Linux evidence — zero divergent lines.

### 6.3 The fifty-nine dynamic notebooks

```bash
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script dynamic_notebooks/<name>.posim
```

Result: 59/59 exit 0; **57 byte-identical** to the Linux evidence. The
two divergent notebooks are `double_slit` (1 line) and `tunneling`
(10 lines), every difference last-digit or 1e-13-scale (e.g. norm
drift `6.393e-13` vs `6.399e-13`). Cause: the **quantum crate** calls
the host libm (Apple libm here, glibc there); the physics engine
routes through the vendored deterministic libm and is byte-identical,
as the other 57 notebooks and all 12 collision scripts show. Pinned
byte-for-byte in `evidence/macos/accepted-divergences-dynamic.diff`;
the gate fails on anything beyond that pinned diff.

### 6.4 The nineteen SolveIt scripts

```bash
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/solveit/01_elastic_head_on.posim
# … through …
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/solveit/19_hinged_door.posim
```

Result: 19/19 exit 0. Scripts 01–16 are **byte-identical** to the
Linux evidence log (zero divergent lines). Scripts 17–18 differ from
their evidence log only because that log is **stale upstream
evidence**: it was recorded before the donor's joint-family upgrade,
and its `CONSTRAIN` status format (`held at 1`) predates the donor's
own current source (`rod … 1 row(s)`) — the donor itself cannot
reproduce it on Linux. Script 19 postdates all evidence logs. The DAE
path is verified byte-identically through the dynamic-notebook suite
(hinges, ball chains, rods, gears) and by the 42 DAE unit tests.

### 6.5 The thirteen recorded browser videos

```bash
cargo build --release -p posim
python3 recorder/src/record_all.py --check
```

Measured result: **all 13 recordings reproduce byte for byte** on this
machine. Against the donor's Linux recordings, **12 of 13 were already
byte-identical**; `rack_and_pinion` was re-recorded here. Cause,
measured: the GEAR/RACK joint residual is `g = sin(qθᵢ + pθⱼ)`
(`physical_object/src/constrain.rs`), evaluated with **host-libm**
`sin`/`cos` on every IDA residual call; Apple libm's one-ulp
differences fork the rounding noise (trajectory agreement ~1e-12; sign
flips only on 1e-26-scale values around zero). The Windows port could
not reproduce the Linux original either. SHA-256 of both recordings:
`evidence/macos/rack_and_pinion-recording-shas.txt`.

Recorder's own tests:

```bash
python3 recorder/tests/test_units.py        # Ran 17 tests … OK
python3 recorder/tests/test_end_to_end.py   # Ran 9 tests … OK
```

## 7. The 109 Jupyter notebooks

```bash
python3 notebooks/_build/regen.py
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/*.ipynb
python3 notebooks/_build/nbcheck.py notebooks/*.ipynb
```

Measured results: `regenerated`; **109 ok, 0 failed**;
**109/109 notebooks pass all seven requirements**.

Open them in JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/
```

## 8. The Jupyter wire protocol and wrapper kernel

```bash
cargo build --release
python3 jupyter/test_protocol.py
```
Measured: `all protocol checks passed`.

```bash
uv venv jupyter/.venv
uv pip install -p jupyter/.venv/bin/python ipykernel jupyter_client
POSIM_NO_BROWSER=1 jupyter/.venv/bin/python jupyter/test_kernel.py
```
Measured: `all kernel checks passed — JupyterLab can drive this kernel`
(7/7 cells ok).

## 9. The thirteen live GUI pages

One server per recorded scene, fixed ports 8895–8907. Run one and open
the printed URL in your browser:

```bash
python3 gui/kepler_ellipse/server.py
```

The automated pass over all 13 (page served with canvas, live state,
Start advances t, Reset returns bit-exactly to t = 0):

```bash
python3 tools/gui_smoke.py
```
Measured: **13 of 13 GUIs pass**.

## 10. The entity index

```bash
open index_of_entities.html
python3 tools/verify_index_examples.py
python3 tools/verify_tierb_examples.py
```
Measured: **PASS 1177/1177 (100.0%)** runnable posim/machine
fragments; **COMPILED 258/258 (100.0%)** Rust snippets.

## 11. The documentation PDFs

```bash
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode SolveIt.tex
pdflatex -interaction=nonstopmode SolveIt.tex
```
Measured: `grammar.pdf` 66 pages, `SolveIt.pdf` 24 pages, zero LaTeX
errors.

## 12. Dependency posture, re-checked here

```bash
grep -c "registry" Cargo.lock sundials_rs/Cargo.lock
# measured: 0 and 0 — no crates.io entry anywhere

grep -rn "unsafe {" --include="*.rs" physical_object/src posim/src quantum/src special_functions/src | wc -l
# measured: 0 — and every crate root carries #![forbid(unsafe_code)],
# so the compiler enforces it:
grep -rn "forbid(unsafe_code)" physical_object/src/lib.rs posim/src/main.rs quantum/src/lib.rs special_functions/src/lib.rs
# measured: all four crate roots
```

`Cargo.lock` lists exactly the local crates: `arkode_rs 7.8.0`,
`cvode_rs 7.8.0`, `sundials_core 7.8.0`, `physical_object 0.1.0`,
`posim 0.1.0`, `quantum 0.1.0`, `spec_math 0.1.6`,
`special_functions 0.1.0`. Nothing is downloaded at build time.
