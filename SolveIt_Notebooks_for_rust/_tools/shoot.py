#!/usr/bin/env python3
"""shoot.py — capture one GUI PNG per encyclopedia notebook.

Calls each notebook's own companion player in `capture` mode (scene replay
for the posim families, served page for mercury/video — mercury pages are
captured jumped to the capture moment). Families with no browser GUI
(rust_*, rebound) are recorded honestly as image=None.

Usage: python3 shoot.py [--workers N] [--only PAT]
Updates manifest.json's "image" fields. Exit 0 if every capturable
notebook produced a non-empty PNG.
"""

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

WS = Path(__file__).resolve().parents[2]
ROOT = WS / "SolveIt_Notebooks_for_rust"
MANIFEST = ROOT / "_tools" / "manifest.json"
IMAGES = ROOT / "encyclopedia" / "images"

GUI_FAMILIES = {"solveit", "dynamic", "collision", "video", "mercury"}


def capture(entry):
    tag = entry["tag"]
    out = IMAGES / f"{tag}.png"
    player = WS / entry["player"]
    try:
        r = subprocess.run(
            [sys.executable, str(player), "capture", str(out)],
            capture_output=True, text=True, timeout=300, cwd=str(player.parent))
    except subprocess.TimeoutExpired:
        return tag, None, "timeout"
    if r.returncode != 0 or not out.exists() or out.stat().st_size < 1000:
        return tag, None, (r.stderr or r.stdout)[-300:]
    return tag, str(out.relative_to(WS)), None


def main():
    workers = 4
    only = None
    args = sys.argv[1:]
    if "--workers" in args:
        workers = int(args[args.index("--workers") + 1])
    if "--only" in args:
        only = args[args.index("--only") + 1]
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    todo = [e for e in entries if e["family"] in GUI_FAMILIES
            and (only is None or only in e["tag"])]
    IMAGES.mkdir(parents=True, exist_ok=True)
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for tag, image, err in ex.map(capture, todo):
            for e in entries:
                if e["tag"] == tag:
                    e["image"] = image
            if image:
                print(f"ok   {tag}")
            else:
                failures.append(tag)
                print(f"FAIL {tag}: {err}")
    MANIFEST.write_text(json.dumps(entries, indent=1) + "\n", encoding="utf-8")
    print(f"{len(todo) - len(failures)} captured, {len(failures)} failed"
          + (f": {failures}" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
