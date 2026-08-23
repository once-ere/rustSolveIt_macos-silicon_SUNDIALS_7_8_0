#!/usr/bin/env python3
"""Regenerate all 109 notebook specs and render them.

The piston-crankshaft spec is hand-written (specs/*.handwritten.py) and
the six rust_* specs carry hand-written pairing text; both are refreshed
by their own generators, so this script rebuilds the other 102 from
their .posim sources and then re-renders every spec present.
"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.parents[1]          # the repository root

def sh(*args):
    r = subprocess.run([sys.executable, *map(str, args)], cwd=ROOT)
    if r.returncode: sys.exit(r.returncode)

SETS = [("videos/scenes", "video"), ("scripts/collisions", "collision"),
        ("scripts/solveit", "solveit"), ("dynamic_notebooks", "dynamic")]
HAND = {"video_piston_crankshaft"}    # hand-written spec: do not overwrite

for d, cat in SETS:
    for f in sorted((ROOT / d).glob("*.posim")):
        key = f"{cat}_{f.stem}"
        if key in HAND: continue
        # Repo-root-relative POSIX path: this repository IS the workspace
        # (no version-7.8.0/ nesting), so the pairs_with metadata carries
        # the path exactly as a reader sees it after cloning.
        rel = f.relative_to(ROOT).as_posix()
        sh(HERE / "nbgen.py", rel, key, cat, rel,
           HERE / "specs" / f"{key}.json")

for spec in sorted((HERE / "specs").glob("*.json")):
    sh(HERE / "nbbuild.py", spec, ROOT / "notebooks")
print("regenerated")
