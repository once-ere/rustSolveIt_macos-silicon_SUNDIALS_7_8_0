#!/usr/bin/env python3
"""Audit the structure rules of the planet_Mercury notebook (Gate 5.2),
mirroring the rustSolveIt engine's own notebook auditor. Standard library
only. Exit 0 only if every rule passes.

Rules:
  R1 the how-to-run instructions are present (needles below)
  R2 no cross-references ("see ... notebook", any *.ipynb name in markdown)
  R3 every code cell is preceded by a markdown cell of >= 80 characters
  R4 the required section headings exist
  R5 the interactive save cell exists (tkinter + asksaveasfilename + fallback)
  R6 every non-interactive code cell carries a non-null execution_count
"""

import json
import re
import sys
from pathlib import Path

NEEDLES_R1 = [
    "jupyter lab",
    "Shift+Enter",
    "Python 3 (ipykernel)",
    "cargo build --release",
]
HEADINGS_R4 = [
    "## 1. What this notebook computes",
    "## 2. How to run this notebook",
    "## 3. The words used in this notebook",
    "## 4. The physical situation",
    "### 4.3 The first-order system actually handed to SUNDIALS",
    "## 5. How this notebook talks to the simulator",
    "## 6. Building and running the simulation",
    "## 7. Storing the results in a database",
    "## 8. Retrieving results from the database",
    "## 9. The verification gauntlet",
    "## 10. Baking and opening the display page",
    "## 11. What we learned",
    "## 12. Name and save this notebook",
]
NEEDLES_R5 = ["tkinter", "asksaveasfilename", "Falling back to a typed folder path"]


def audit(path: Path) -> list:
    problems = []
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return [f"unreadable notebook JSON: {exc!r}"]
    cells = nb.get("cells", [])
    all_md = "\n".join(
        "".join(c.get("source", [])) for c in cells if c.get("cell_type") == "markdown"
    )
    all_text = "\n".join("".join(c.get("source", [])) for c in cells)

    for needle in NEEDLES_R1:
        if needle not in all_text:
            problems.append(f"R1: missing how-to-run needle {needle!r}")
    if re.search(r"see\s+(?:the\s+)?[\w\- ]*notebook\b", all_md, re.IGNORECASE):
        problems.append("R2: markdown refers the reader to another notebook")
    for m in re.findall(r"\b[\w\-]+\.ipynb\b", all_md):
        # Naming ITSELF (launch/batch instructions) is self-reference, fine;
        # naming any other notebook is a forbidden cross-reference.
        if m != path.name:
            problems.append(f"R2: markdown names another notebook file: {m}")
    prev = None
    for i, c in enumerate(cells):
        if c.get("cell_type") == "code":
            tags = c.get("metadata", {}).get("tags", [])
            if prev is None or prev.get("cell_type") != "markdown" or len(
                "".join(prev.get("source", [])).strip()
            ) < 80:
                problems.append(f"R3: code cell {i} lacks a >=80-char markdown lead-in")
            if "interactive" not in tags and c.get("execution_count") is None:
                problems.append(f"R6: code cell {i} was never executed")
        prev = c
    for h in HEADINGS_R4:
        if h not in all_md:
            problems.append(f"R4: missing heading {h!r}")
    for needle in NEEDLES_R5:
        if needle not in all_text:
            problems.append(f"R5: missing save-cell needle {needle!r}")
    return problems


def main() -> int:
    paths = [Path(p) for p in sys.argv[1:]]
    if not paths:
        print("usage: check_notebook.py <notebook.ipynb> [...]")
        return 2
    bad = 0
    for p in paths:
        problems = audit(p)
        if problems:
            bad += 1
            print(f"FAIL {p}:")
            for x in problems:
                print(f"  - {x}")
        else:
            print(f"ok {p}: all structure rules pass")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
