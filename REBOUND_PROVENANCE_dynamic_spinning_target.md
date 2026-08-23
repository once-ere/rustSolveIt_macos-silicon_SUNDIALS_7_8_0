# REBOUND_PROVENANCE — notebooks/dynamic_spinning_target.ipynb

Rebound physics, verified on **macOS / Apple Silicon (Apple M5 Max,
arm64, macOS Tahoe 26)**, 2026-08-23. This page is self-contained:
every command needed to test, verify, execute and display this
notebook appears on this page, run from the repository root.

## What this notebook is

The Jupyter notebook paired with the dynamic notebook script
`dynamic_notebooks/spinning_target.posim` — a **rebound off a moving
surface**. A cuboid (m = 2, half-extents 0.6) spins at ω_z = 2; a
sphere (m = 1, r = 0.3) flies in at [2, 0, 0], aimed above centre
(y = 0.4). The surface point the sphere meets is itself moving
(ω × r), so the approach speed at contact is 2.40, not 2 — the
K-matrix's angular terms handle exactly this. The t ≈ 0.97 impact
trades the target's spin for the sphere's rebound, yet the conserved
totals stay pinned through the collision.

## Measured on this machine

- Cuboid's own spin before impact: `I_z = (2/3)(0.36+0.36) = 0.48`,
  so `L_z = 0.96` — printed as `[0, 0, 0.96]`.
- Conserved totals, held through the impact: `L = [0, 0, 0.16]`,
  `P = [2, 0, 0]`, `E = 2.96` (0.96 spin + 2 translational).
- The impact trades spin for rebound: cuboid `L_z` 0.96 → 0.31.
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
POSIM_NO_BROWSER=1 cargo run -p posim --release -- --script dynamic_notebooks/spinning_target.posim
```

Execute the Jupyter notebook headlessly and write its outputs back
into it, then check it (stdlib-Python tools vendored here — nothing to
install):

```bash
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/dynamic_spinning_target.ipynb
python3 notebooks/_build/nbcheck.py notebooks/dynamic_spinning_target.ipynb
```

Display the notebook in a web browser with JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/dynamic_spinning_target.ipynb
```

Display the rebound live in this project's browser GUI — the scene
window opens in your browser; press Start (or type `SCENE START`) to
watch the sphere strike the turning face and rebound while the hud's
E, P, L readouts hold their values through the impact; press C for the
golden contact arrows:

```bash
tools/posim_notebook spinning_target
```

Verify byte-identity of the whole 59-script dynamic-notebook suite
against the Linux evidence (this script is one of them; exit 0 =
pass):

```bash
bash tools/macos_verify_physics.sh
```
