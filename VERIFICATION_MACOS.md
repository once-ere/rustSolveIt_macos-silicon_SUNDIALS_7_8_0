# VERIFICATION_MACOS — every suite, every verdict, on this machine

Platform: **macOS Tahoe 26 (26.6.1, build 25G76) on Apple Silicon
(Apple M5 Max, arm64, 128 GB RAM)**; rustc/cargo 1.94.0
(`aarch64-apple-darwin`); Python 3.14.7; TeX Live 2026. Date:
2026-08-23. Every command below runs from the repository root, and
every number on this page came out of a run on this machine.

> **Platform scope.** "Byte-identical" on this page means identical to
> the Linux release's recorded evidence (`evidence/port-7.8.0/`),
> reproduced on macOS/arm64. The engine's vendored deterministic math
> library makes the solver numerics host-independent; the only
> divergences are the three documented host-libm cases in §4, each
> pinned byte-for-byte so the gate fails on anything new.

## 1. The scoreboard

| gate | command | result |
|---|---|---|
| build, debug | `cargo build --workspace --all-targets` | 0 errors, **0 warnings** |
| build, release | `cargo build --release --workspace --all-targets` | 0 errors, **0 warnings** |
| unit + integration tests | `cargo test --workspace` | **622 passed, 0 failed** |
| engine sanity (sibling checkout of the Apple-Silicon SUNDIALS port) | `cargo build --workspace && cargo test --workspace --lib` | 0 warnings; **25 passed, 0 failed** |
| 6 self-checking physics examples | `bash tools/macos_verify_physics.sh` | 6/6 SUCCESS; byte-identical modulo one pinned blank line (§4.1) |
| 12 collision scripts | `bash tools/macos_verify_physics.sh` | 12/12 exit 0; **byte-identical, zero divergent lines** |
| 59 dynamic notebooks | `bash tools/macos_verify_physics.sh` | 59/59 exit 0; **57 byte-identical**, quantum pair pinned (§4.2) |
| 19 SolveIt scripts | `POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/solveit/NN_name.posim` | 19/19 exit 0; **01–16 byte-identical**; 17–19 see §4.4 |
| 13 recorded browser videos | `python3 recorder/src/record_all.py --check` | **all 13 reproduce byte for byte**; 12/13 byte-identical to the Linux recordings (§4.3) |
| recorder tests | `python3 recorder/tests/test_units.py` and `python3 recorder/tests/test_end_to_end.py` | 17 OK; 9 OK |
| 109 Jupyter notebooks, executed | `POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/*.ipynb` | **109 ok, 0 failed** |
| 109 Jupyter notebooks, checked | `python3 notebooks/_build/nbcheck.py notebooks/*.ipynb` | **109/109 pass all seven requirements** |
| 13 live GUI servers | `python3 tools/gui_smoke.py` | **13/13 pass** (page + canvas, state advances, Reset returns bit-exactly to t = 0) |
| entity index, posim fragments | `python3 tools/verify_index_examples.py` | **PASS 1177/1177 (100.0%)** |
| entity index, Rust snippets | `python3 tools/verify_tierb_examples.py` | **COMPILED 258/258 (100.0%)** |
| Jupyter wire protocol | `python3 jupyter/test_protocol.py` | all protocol checks passed |
| JupyterLab wrapper kernel | `POSIM_NO_BROWSER=1 jupyter/.venv/bin/python jupyter/test_kernel.py` | all kernel checks passed (7/7 cells) |
| documentation PDFs | `pdflatex -interaction=nonstopmode grammar.tex` (×2), same for `SolveIt.tex` | grammar.pdf 66 pp., SolveIt.pdf 24 pp., zero errors |

## 2. The byte-identity gate

```bash
cargo build --release -p posim
bash tools/macos_verify_physics.sh
```

re-runs the three suites whose Linux outputs are recorded in
`evidence/port-7.8.0/`, writes this machine's outputs to
`evidence/macos/`, and byte-diffs them (the OS-assigned scene port and
the current directory are normalised, exactly as the originals were).
Measured verdict, exit 0:

```
IDENTICAL-MODULO-DOCUMENTED  examples-macos.log vs examples-7.8.0.log (pinned: accepted-divergences-examples.diff)
IDENTICAL  collision-scripts-macos.log == collision-scripts-7.8.0.log
IDENTICAL-MODULO-DOCUMENTED  dynamic-notebooks-macos.log vs dynamic-notebooks-7.8.0.log (pinned: accepted-divergences-dynamic.diff)
```

A divergence is ACCEPTED only if it matches its pinned diff file in
`evidence/macos/` **byte for byte**; anything else fails the gate.

## 3. Regression anchors, re-measured here (identical to Linux)

- Outer solar system: Pluto `x = 31.78592516  y = 38.63618957
  z = 3.19279415` at t = 500,000 days; energy drift `7.835809e-07`;
  `12581` internal steps.
- Kepler e = 0.6: `|dE/E| = 2.430222e-09`, `|dL|/|L| = 5.144452e-08`,
  `|dA|/|A| = 1.131858e-07` (Runge–Lenz).
- Tumbling body: `|dL|/|L| = 0.000000e+00` — exact.
- Ball on plate TOI `0.8944271909999157` vs analytic `sqrt(0.8)` —
  1.24e-16 relative.
- Sixteen SolveIt worked-example scripts: byte-identical transcripts.

## 4. The documented divergences — all of them

### 4.1 examples: one blank line (pinned: `accepted-divergences-examples.diff`, 2 lines)

The Linux evidence for `outer_solar_system` lacks a blank line that
the donor later added (a leading `"\n"` in the example's first
`println!`). The donor's own current source cannot reproduce its
evidence on Linux either; every number is identical. The Windows port
pinned the same divergence.

### 4.2 dynamic notebooks: the quantum pair (pinned: `accepted-divergences-dynamic.diff`, 44 lines)

11 output lines across `double_slit` (1) and `tunneling` (10) differ
in the last printed digit or at 1e-13 scale (e.g. norm drift
`6.393e-13` → `6.399e-13`). The **quantum crate** calls the host libm
(Apple libm here, glibc on the Linux host); the physics engine routes
through the vendored deterministic libm and is byte-identical, as the
other 57 notebooks and all 12 collision scripts show. The Windows port
recorded the same two notebooks as its only dynamic divergence.

### 4.3 videos: `rack_and_pinion` re-recorded on this machine

12 of 13 Linux-recorded videos reproduce **byte-for-byte** on Apple
Silicon. The rack-and-pinion scene alone is host-sensitive: the
GEAR/RACK joint residual `g = sin(qθᵢ + pθⱼ)`
(`physical_object/src/constrain.rs`) evaluates host-libm `sin`/`cos`
inside every IDA residual call, and a one-ulp libm difference forks
the step-size trajectory's rounding noise. Measured: trajectory
agreement ~1e-12 relative; the only sign changes are on 1e-26-scale
values around zero; final energy differs by ~1e-15 relative. The
Windows port could not reproduce the Linux original either. The video
was re-recorded here; SHA-256 of both recordings is pinned in
`evidence/macos/rack_and_pinion-recording-shas.txt`, and
`python3 recorder/src/record_all.py --check` now reports **all 13
reproduce byte for byte** on this machine.

### 4.4 SolveIt scripts 17–19: stale upstream evidence, not drift

Scripts 17–18 differ from `evidence/port-7.8.0/dae-examples.log` only
because that log predates the donor's joint-family upgrade: its
`CONSTRAIN` status format (`constraint0: obj0 <-> obj1 held at 1`)
does not match the donor's own **current** source (`constraint0: rod
obj0 <-> obj1, 1 row(s)`), so the log is unreproducible on Linux too.
Script 19 postdates every evidence log. All three exit 0, and the DAE
path they exercise is verified byte-identically by the dynamic
notebook suite (rod chains, hinges, ball joints, gears) and by the 42
constrained/DAE unit tests.

## 5. Reproduce everything, in order

```bash
cargo build --workspace --all-targets
cargo build --release --workspace --all-targets
cargo test --workspace
bash tools/macos_verify_physics.sh
for f in scripts/solveit/*.posim; do POSIM_NO_BROWSER=1 target/release/posim --script "$f" || echo "FAILED $f"; done
python3 recorder/src/record_all.py --check
python3 recorder/tests/test_units.py
python3 recorder/tests/test_end_to_end.py
python3 notebooks/_build/regen.py
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/*.ipynb
python3 notebooks/_build/nbcheck.py notebooks/*.ipynb
python3 tools/gui_smoke.py
python3 tools/verify_index_examples.py
python3 tools/verify_tierb_examples.py
python3 jupyter/test_protocol.py
uv venv jupyter/.venv
uv pip install -p jupyter/.venv/bin/python ipykernel jupyter_client
POSIM_NO_BROWSER=1 jupyter/.venv/bin/python jupyter/test_kernel.py
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode grammar.tex
pdflatex -interaction=nonstopmode SolveIt.tex
pdflatex -interaction=nonstopmode SolveIt.tex
```
