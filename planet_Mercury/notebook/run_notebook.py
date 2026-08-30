#!/usr/bin/env python3
"""Execute every code cell of a Jupyter notebook in order, embedding the real
captured outputs, exactly like the rustSolveIt engine's own notebook runner:
cells are exec()'d in one shared namespace with stdout/stderr captured; cells
tagged "interactive" are skipped (their outputs blanked); the file is written
back ONLY if every executed cell succeeded. Standard library only.

Usage:  MERCURY_NO_BROWSER=1 python3 run_notebook.py mercury_tidal_locking.ipynb
"""

import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


def run(path: Path) -> bool:
    nb = json.loads(path.read_text(encoding="utf-8"))
    ns: dict = {"__name__": "__main__"}
    count = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        tags = cell.get("metadata", {}).get("tags", [])
        src = "".join(cell.get("source", []))
        if "interactive" in tags:
            cell["outputs"] = []
            cell["execution_count"] = None
            continue
        count += 1
        buf = io.StringIO()
        try:
            with redirect_stdout(buf), redirect_stderr(buf):
                exec(compile(src, f"<cell {count}>", "exec"), ns)
        except BaseException as exc:  # noqa: BLE001 - report and fail loudly
            text = buf.getvalue()
            print(f"FAIL {path} (cell {count}): {exc!r}")
            if text:
                print("--- captured output of the failing cell ---")
                print(text)
            return False
        cell["execution_count"] = count
        text = buf.getvalue()
        cell["outputs"] = (
            [{"output_type": "stream", "name": "stdout", "text": text.splitlines(keepends=True)}]
            if text
            else []
        )
    sim = ns.get("sim")
    if sim is not None and hasattr(sim, "close"):
        try:
            sim.close()
        except Exception:
            pass
    path.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
    print(f"ok {path} ({count} cells)")
    return True


def main() -> int:
    paths = [Path(p) for p in sys.argv[1:]]
    if not paths:
        print("usage: run_notebook.py <notebook.ipynb> [...]")
        return 2
    ok = 0
    for p in paths:
        if run(p):
            ok += 1
    print(f"{ok} ok, {len(paths) - ok} failed")
    return 0 if ok == len(paths) else 1


if __name__ == "__main__":
    sys.exit(main())
