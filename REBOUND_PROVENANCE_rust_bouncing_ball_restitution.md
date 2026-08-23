# REBOUND_PROVENANCE — notebooks/rust_bouncing_ball_restitution.ipynb

Rebound physics, verified on **macOS / Apple Silicon (Apple M5 Max,
arm64, macOS Tahoe 26)**, 2026-08-23. This page is self-contained:
every command needed to test, verify, execute and display this
notebook appears on this page, run from the repository root.

## What this notebook is

The Jupyter notebook paired with the compiled, self-checking Rust
example `physical_object/examples/bouncing_ball_restitution.rs`. A
ball (radius 0.5, mass 1) is dropped so its centre falls h = 4.5 onto
a static floor under gravity g = 10, with coefficient of restitution
e = 0.8. The physics of a **rebound**: the impact must land at the
analytic time of impact `t = sqrt(2h/g) = sqrt(0.9)`, and the rebound
apex must be `e² × h` above the rest height, because the impulse
response scales the separation speed by exactly e. The example checks
both against closed forms and prints SUCCESS or FAILURE; the notebook
compiles and runs that example, asserts SUCCESS, and reproduces the
same physics interactively over the simulator's JSON machine protocol.

## Measured on this machine

```
Bouncing ball: h = 4.5, g = 10, e = 0.8
  impact time  : 0.9486832980505138 (analytic 0.9486832980505138)
  impact normal: [-0, 1, -0] (floor → ball)
  rebound apex : 3.37997191386479 (analytic 3.380000000000001)
SUCCESS: rebound apex = e^2 x drop height, impact at the exact TOI
```

The impact time matches the analytic value **to the last digit**, and
this output is **byte-identical** to the Linux release's recorded
evidence. Notebook execution: ok, all cells; checker: passes all seven
requirements.

## Complete commands

Build the simulator and the examples, once:

```bash
cargo build --release --workspace --all-targets
```

Run the self-checking example directly (prints the block quoted above,
exits 0 on SUCCESS and nonzero on FAILURE):

```bash
cargo run -p physical_object --release --example bouncing_ball_restitution
```

Execute the notebook headlessly and write its outputs back into it,
then check it (both are stdlib-Python tools vendored in this
repository — nothing to install):

```bash
POSIM_NO_BROWSER=1 python3 notebooks/_build/nbrun.py notebooks/rust_bouncing_ball_restitution.ipynb
python3 notebooks/_build/nbcheck.py notebooks/rust_bouncing_ball_restitution.ipynb
```

Display the notebook in a web browser with JupyterLab:

```bash
python3 -m pip install --user jupyterlab
jupyter lab notebooks/rust_bouncing_ball_restitution.ipynb
```

Display the same rebound live in this project's browser GUI — the
scene window with Start/Pause/Reverse/Reset and live E, P, L readouts
(the dynamic-notebook twin of this example; press Start in the window):

```bash
tools/posim_notebook bouncing_ball_restitution
```

Verify byte-identity of the whole examples suite against the Linux
evidence (this example is one of the six it covers; exit 0 = pass):

```bash
bash tools/macos_verify_physics.sh
```
