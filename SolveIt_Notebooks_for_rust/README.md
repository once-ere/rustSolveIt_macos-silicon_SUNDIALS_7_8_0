# SolveIt_Notebooks_for_rust — the notebook encyclopedia

**128 Jupyter notebooks** (128/128 executed green on macOS/Apple Silicon),
organized into seven major topics, each in its own tagged subfolder with a
companion `player_<tag>.py`, plus the encyclopedia itself:

- `ENCYCLOPEDIA.pdf` / `ENCYCLOPEDIA.md` / `ENCYCLOPEDIA.html` — one entry
  per notebook: description, execution verdict, GUI image, and every
  pointer (source, paired scene, movie, database, player).
- `encyclopedia/images/` — 105 captured GUI screenshots (the 6 `rust_*`
  and 17 rebound notebooks have textual verdicts, no browser GUI).
- `_tools/` — the deterministic pipeline: `gen.py` (copies + header cells
  + players), `run_copy.py` (family-aware executor), `shoot.py` (GUI
  capture), `encyclo.py` (encyclopedia builder), `verify.py` (six
  mechanical gates), `manifest.json`, and the execution logs the
  encyclopedia's verdicts come from.

Every notebook's **first cell** is its encyclopedia header: the exact
commands to execute it, display its browser GUI, create/access its movies,
access its database/data — and the full player script inline, byte-identical
to the shipped `player_<tag>.py`.

## Two layouts, one rule

The paths in the header cells and players resolve automatically in both
places this tree lives:

1. **The original workspace** — a folder containing both
   `SolveIt_Notebooks_for_rust/` and the engine repository directory
   `rustSolveIt_macos-silicon_SUNDIALS_7_8_0/`. Every header command runs
   verbatim from that workspace root.
2. **A fresh clone of the repository** — here the clone root IS the engine
   (it contains `posim/`, `notebooks/`, `planet_Mercury/`, and this
   folder). The tools detect this and resolve engine paths to the clone
   root. Run commands **from the clone root**, with two substitutions:
   the one-time engine build is simply
   `cargo build --release -p posim`
   (and `cd planet_Mercury/mercury_rs && cargo build --release` for the
   Mercury pair); every command that starts with
   `SolveIt_Notebooks_for_rust/...` runs verbatim.

## Reproduce everything

```bash
# one-time builds (from the engine root):
cargo build --release -p posim
cd planet_Mercury/mercury_rs && cargo build --release && cd ../..

# execute any notebook copy (writes its real outputs back):
python3 SolveIt_Notebooks_for_rust/_tools/run_copy.py \
        SolveIt_Notebooks_for_rust/<topic>/<tag>/<tag>.ipynb

# per-notebook player: open GUI, jump to the key event, capture, export
python3 SolveIt_Notebooks_for_rust/<topic>/<tag>/player_<tag>.py gui
python3 SolveIt_Notebooks_for_rust/<topic>/<tag>/player_<tag>.py jump
python3 SolveIt_Notebooks_for_rust/<topic>/<tag>/player_<tag>.py capture out.png
python3 SolveIt_Notebooks_for_rust/<topic>/<tag>/player_<tag>.py data exported

# rebuild the whole encyclopedia:
python3 SolveIt_Notebooks_for_rust/_tools/shoot.py
python3 SolveIt_Notebooks_for_rust/_tools/encyclo.py
python3 SolveIt_Notebooks_for_rust/_tools/verify.py
```

Headless GUI capture and the PDF build use Google Chrome; everything else
is standard-library Python 3 plus the pure-Rust engine. The rebound family
additionally needs `jupyterlab` (and `matplotlib` for its two plotting
notebooks) in `planet_Mercury/notebook/.venv`.

## Deliberately not committed

Two notebooks bake large regenerable artifacts beside themselves when
executed — `dynamic_double_slit/double_slit.html` (~49 MB) and
`dynamic_tunneling/scatter.html` — one run of each notebook recreates
them. The Mercury notebooks likewise regenerate their SQLite databases,
CSV runs, and player pages under `planet_Mercury/` (documented there).
