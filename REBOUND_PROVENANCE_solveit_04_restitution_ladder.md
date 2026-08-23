# REBOUND_PROVENANCE — notebooks/solveit_04_restitution_ladder.ipynb

Rebound physics, verified on **macOS / Apple Silicon (Apple M5 Max,
arm64, macOS Tahoe 26)**, 2026-08-23. This page is self-contained:
every command needed to test, verify, execute and display this
notebook appears on this page, run from the repository root.

## What this notebook is

The Jupyter notebook paired with the SolveIt worked example
`scripts/solveit/04_restitution_ladder.posim` (worked example 4 of the
SolveIt solution guide). A ball with restitution e = 0.7 is dropped so
its centre starts at 5.5 (radius 0.5) and falls 5.0 to the floor's top
face under g = 10. Each **rebound** reaches `e^(2n)` of the drop
height on bounce n — the restitution ladder:

- bounce 1 at t = 1.0, rebound speed 7.00, apex at t = 1.7 —
  height above floor `5·e² = 2.45`, centre at 2.95
- bounce 2 at t = 2.4, rebound speed 4.90, apex at t = 2.89 —
  height `5·e⁴ = 1.2005`, centre at 1.7005
- bounce 3 at t = 3.38, rebound speed 3.43, apex at t = 3.723 —
  height `5·e⁶`, centre at 1.088245

The notebook derives the ladder, runs the three bounces, and compares
each measured apex against its closed form.

## Measured on this machine

| bounce | measured centre apex | analytic centre apex |
|---|---|---|
| 1 | `2.949999999869974` | `2.9499999999999997` |
| 2 | `1.7004999998062575` | `1.7004999999999997` |
| 3 | `1.0882449997748824` | `1.0882449999999997` |

Three collisions counted; every apex agrees with its closed form to
~1e-10 (the apex is sampled on the output grid, which is what limits
it — the impact times themselves are solver rootfinding events). The
script's transcript is **byte-identical** to the Linux release's
recorded evidence (SolveIt scripts 01–16 diff clean). Notebook
execution: ok, all cells; checker: passes all seven requirements.

## Complete commands

Build the simulator, once:

```bash
cargo build --release -p posim
```

Run the paired worked-example script headlessly and read its
transcript (also documented, with captured output, as example 4 in the
SolveIt guide):

```bash
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script scripts/solveit/04_restitution_ladder.posim
```

Execute the Jupyter notebook headlessly and write its outputs back
into it, then check it (stdlib-Python tools vendored here — nothing to
install):

```bash
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/solveit_04_restitution_ladder.ipynb
python3 notebooks/_build/nbcheck.py notebooks/solveit_04_restitution_ladder.ipynb
```

Display the notebook in a web browser with JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/solveit_04_restitution_ladder.ipynb
```

Display a bouncing-ball rebound live in this project's browser GUI —
the scene window opens in your browser with Start/Pause/Reverse/Reset
and live E, P, L readouts (the dynamic-notebook variant of the same
physics; press Start in the window):

```bash
tools/posim_notebook bouncing_ball_restitution
```

Verify byte-identity of the SolveIt scripts against the Linux
evidence — run all nineteen and confirm each exits 0 (01–16 are
byte-compared by the record in `evidence/`):

```bash
for f in scripts/solveit/*.posim; do POSIM_NO_BROWSER=1 target/release/posim --script "$f" || echo "FAILED $f"; done
```
