# REBOUND_PROVENANCE — notebooks/dynamic_bouncing_ball_restitution.ipynb

Rebound physics, verified on **macOS / Apple Silicon (Apple M5 Max,
arm64, macOS Tahoe 26)**, 2026-08-23. This page is self-contained:
every command needed to test, verify, execute and display this
notebook appears on this page, run from the repository root.

## What this notebook is

The Jupyter notebook paired with the dynamic notebook script
`dynamic_notebooks/bouncing_ball_restitution.posim`. A ball (radius
0.5, mass 1) falls from rest at y = 5 onto a static floor slab under
uniform gravity g = 10, with restitution e = 0.8. The centre drops
h = 4.5 to its rest height 0.5, so the first impact lands at
`t = sqrt(2h/g) = sqrt(0.9) ≈ 0.9487`; the **rebound** speed is 0.8 of
the impact speed, so the next apex is `0.5 + e²·4.5 = 3.38`, and each
later apex is lower by the same geometric factor. The Jupyter notebook
drives the simulator over its JSON machine protocol, derives the
closed forms in its own text, and records the real outputs.

## Measured on this machine

- Time of first impact: `0.9486832980505138` — the analytic
  `sqrt(0.9)` to the last digit.
- Rebound apex of the centre: `3.38` — exactly `0.5 + e²·4.5`.
- Initial energy 50 (pure potential, m·g·y), dropping at every impact.
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
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script dynamic_notebooks/bouncing_ball_restitution.posim
```

Execute the Jupyter notebook headlessly and write its outputs back
into it, then check it (stdlib-Python tools vendored here — nothing to
install):

```bash
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/dynamic_bouncing_ball_restitution.ipynb
python3 notebooks/_build/nbcheck.py notebooks/dynamic_bouncing_ball_restitution.ipynb
```

Display the notebook in a web browser with JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/dynamic_bouncing_ball_restitution.ipynb
```

Display the rebound live in this project's browser GUI — the scene
window opens in your browser with Start/Pause/Reverse/Reset and live
E, P, L readouts; press Start (or type `SCENE START`) and watch the
ball climb to just below 3.4, then geometrically lower each bounce;
press C for the golden action–reaction contact arrows:

```bash
tools/posim_notebook bouncing_ball_restitution
```

Verify byte-identity of the whole 59-script dynamic-notebook suite
against the Linux evidence (this script is one of them; exit 0 =
pass):

```bash
bash tools/macos_verify_physics.sh
```
