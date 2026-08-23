# REBOUND_PROVENANCE — notebooks/dynamic_routh_p1_geometric_progression.ipynb

Rebound physics, verified on **macOS / Apple Silicon (Apple M5 Max,
arm64, macOS Tahoe 26)**, 2026-08-23. This page is self-contained:
every command needed to test, verify, execute and display this
notebook appears on this page, run from the repository root.

## What this notebook is

The Jupyter notebook paired with
`dynamic_notebooks/routh_p1_geometric_progression.posim` — E. J.
Routh, *A Treatise on Dynamics of a Particle* (Cambridge, 1898),
Part I, Art. 88: perfectly elastic balls in a line, each impinging
directly on the next; if their masses form a geometric progression of
common ratio 2, the velocities after impact form a geometric
progression. Four spheres of masses 1, 2, 4, 8 (restitution 1, no
gravity) realize the theorem: each impact passes 2/3 of the incoming
speed to the struck ball and **rebounds** the striker backwards, so
the struck-ball speeds run 2/3, (2/3)², (2/3)³ of the original — the
geometric progression Routh asks for. The notebook derives the impact
equations of Arts. 83–84 in full and checks the measured ratios.

## Measured on this machine

- Final velocity of the last (mass-8) ball:
  `0.2962962962962963 = 8/27` — the measured ratio against the exact
  value prints `1`.
- Rebound of the third ball after it strikes the fourth:
  `-0.14814814814814814` (exactly −4/27 — it bounces backwards).
- `3 collision(s)` through the chain; total momentum `[1, 0, 0]`
  conserved exactly; total energy `0.49999999999999994` against the
  initial 0.5; at tight tolerance the run reports
  `|dE/E| = 1.110e-16` — one machine epsilon.
- The script's transcript is **byte-identical** to the Linux release's
  recorded evidence. Notebook execution: ok, all cells; checker:
  passes all seven requirements.

## Complete commands

Build the simulator, once:

```bash
cargo build --release -p posim
```

Run the paired script headlessly and read its transcript:

```bash
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script dynamic_notebooks/routh_p1_geometric_progression.posim
```

Execute the Jupyter notebook headlessly and write its outputs back
into it, then check it (stdlib-Python tools vendored here — nothing to
install):

```bash
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/dynamic_routh_p1_geometric_progression.ipynb
python3 notebooks/_build/nbcheck.py notebooks/dynamic_routh_p1_geometric_progression.ipynb
```

Display the notebook in a web browser with JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/dynamic_routh_p1_geometric_progression.ipynb
```

Display the rebound chain live in this project's browser GUI — the
scene window opens in your browser; press Start (or type
`SCENE START`) to watch the impulse run down the chain of doubling
masses; press C for the golden contact arrows at each impact:

```bash
tools/posim_notebook routh_p1_geometric_progression
```

Verify byte-identity of the whole 59-script dynamic-notebook suite
against the Linux evidence (this script is one of them; exit 0 =
pass):

```bash
bash tools/macos_verify_physics.sh
```
