#!/usr/bin/env python3
"""verify.py — mechanical gate over the whole encyclopedia tree.

Checks (read-only):
  V1  manifest has 128 entries; topic counts are 2/19/59/12/13/6/17
  V2  every copy exists, parses, and its FIRST cell is the tagged header
      with all five required sections
  V3  every player_<tag>.py exists and compiles; the header cell embeds
      the identical script text
  V4  every code cell of every copy carries an execution_count
      (interactive-tagged cells exempt)
  V5  every GUI-family entry has a non-trivial PNG; rust/rebound have none
  V6  ENCYCLOPEDIA.md exists and contains every tag exactly once
Exit 0 only if all pass.
"""

import json
import re
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[2]
ROOT = WS / "SolveIt_Notebooks_for_rust"

REQUIRED_SECTIONS = [
    "## 1. Execute this notebook",
    "## 2. Display the browser GUI",
    "## 3. Movies",
    "## 4. Database / data access",
    "## 5. The full player",
]
EXPECT = {"01_planet_mercury_tidal_locking": 2, "02_solveit_worked_examples": 19,
          "03_dynamics_and_routh": 59, "04_collisions": 12,
          "05_mechanism_videos": 13, "06_rust_compiled_examples": 6,
          "07_nbody_rebound_rust": 17}

problems = []


def check(ok, label):
    if not ok:
        problems.append(label)
    return ok


def main():
    entries = json.loads((ROOT / "_tools" / "manifest.json").read_text())
    check(len(entries) == 128, f"V1: {len(entries)} entries, expected 128")
    from collections import Counter
    counts = Counter(e["topic"] for e in entries)
    for t, n in EXPECT.items():
        check(counts[t] == n, f"V1: {t} has {counts[t]}, expected {n}")

    md = (ROOT / "ENCYCLOPEDIA.md").read_text(encoding="utf-8") \
        if (ROOT / "ENCYCLOPEDIA.md").exists() else ""
    check(bool(md), "V6: ENCYCLOPEDIA.md missing")

    for e in entries:
        tag = e["tag"]
        copy = WS / e["copy"]
        if not check(copy.exists(), f"V2 {tag}: copy missing"):
            continue
        nb = json.loads(copy.read_text(encoding="utf-8"))
        first = nb["cells"][0]
        htext = "".join(first["source"])
        check(first["cell_type"] == "markdown" and
              "encyclopedia-header" in first["metadata"].get("tags", []),
              f"V2 {tag}: first cell is not the tagged header")
        for s in REQUIRED_SECTIONS:
            check(s in htext, f"V2 {tag}: header lacks '{s}'")
        player = WS / e["player"]
        if check(player.exists(), f"V3 {tag}: player missing"):
            ptext = player.read_text(encoding="utf-8")
            try:
                compile(ptext, str(player), "exec")
            except SyntaxError as exc:
                check(False, f"V3 {tag}: player does not compile: {exc}")
            m = re.search(r"```python\n(.*?)\n```", htext, re.S)
            check(m is not None and m.group(1).strip() == ptext.strip(),
                  f"V3 {tag}: header-embedded script != shipped player")
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] != "code":
                continue
            if "interactive" in c.get("metadata", {}).get("tags", []):
                continue
            check(c.get("execution_count") is not None,
                  f"V4 {tag}: code cell {i} never executed")
        if e["family"] in ("solveit", "dynamic", "collision", "video", "mercury"):
            img = ROOT / "encyclopedia" / "images" / f"{tag}.png"
            check(img.exists() and img.stat().st_size > 10000,
                  f"V5 {tag}: GUI image missing/trivial")
        else:
            check(e.get("image") is None, f"V5 {tag}: unexpected image record")
        check(md.count(f"### `{tag}`") == 1,
              f"V6 {tag}: appears {md.count(f'### `{tag}`')}x in ENCYCLOPEDIA.md")

    if problems:
        print(f"FAIL — {len(problems)} problems:")
        for p in problems[:40]:
            print("  -", p)
        return 1
    print(f"ALL VERIFY GATES PASS ({len(entries)} notebooks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
