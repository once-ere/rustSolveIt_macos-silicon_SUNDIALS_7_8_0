#!/usr/bin/env python3
"""run_copy.py — execute one encyclopedia notebook copy, family-aware.

Runs every non-interactive code cell in order in one shared namespace
(exactly the engine's own miniature-kernel discipline), captures the real
outputs, and writes them back into the copy ONLY if every cell succeeded.
The rebound family is executed through nbclient (its notebooks were written
for a real Jupyter kernel), with the kernel's working directory set to the
canonical rebound_rust/notebooks folder so their relative paths resolve.

Usage:  python3 run_copy.py <copy.ipynb> [...]
Exit 0 only if every notebook passed.
"""

import contextlib
import io
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path

WS = Path(__file__).resolve().parents[2]
ENGINE = WS / "rustSolveIt_macos-silicon_SUNDIALS_7_8_0"
if not ENGINE.exists():
    ENGINE = WS                # fresh clone: the repository root IS the engine


def family_of(path: Path) -> str:
    topic = path.resolve().parents[1].name
    return {
        "01_planet_mercury_tidal_locking": "mercury",
        "02_solveit_worked_examples": "posim",
        "03_dynamics_and_routh": "posim",
        "04_collisions": "posim",
        "05_mechanism_videos": "posim",
        "06_rust_compiled_examples": "posim",
        "07_nbody_rebound_rust": "rebound",
    }[topic]


def stream(text):
    return {"output_type": "stream", "name": "stdout",
            "text": text.splitlines(keepends=True)}


def run_exec(path: Path, cwd: Path, env: dict):
    """The engine's nbrun discipline: exec cells in one namespace."""
    nb = json.loads(path.read_text(encoding="utf-8"))
    ns = {"__name__": "__main__"}
    count = 0
    ok = True
    old_cwd = os.getcwd()
    old_env = {k: os.environ.get(k) for k in env}
    os.chdir(cwd)
    os.environ.update(env)
    try:
        for cell in nb["cells"]:
            if cell["cell_type"] != "code":
                continue
            if "interactive" in cell.get("metadata", {}).get("tags", []):
                cell["outputs"] = []
                cell["execution_count"] = None
                continue
            src = "".join(cell["source"])
            count += 1
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                    exec(compile(src, f"{path.name}:cell{count}", "exec"), ns)
            except BaseException:
                ok = False
                print(f"FAIL {path} cell {count}", file=sys.stderr)
                print(buf.getvalue(), file=sys.stderr)
                traceback.print_exc()
                break
            cell["execution_count"] = count
            text = buf.getvalue()
            cell["outputs"] = [stream(text)] if text else []
    finally:
        sim = ns.get("sim")
        if sim is not None:
            try:
                sim.close()
            except Exception:
                pass
        os.chdir(old_cwd)
        for k, v in old_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    if ok:
        path.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
    return ok, count


def run_rebound(path: Path):
    """nbclient execution with the canonical folder as kernel cwd."""
    venv_py = WS / "planet_Mercury" / "notebook" / ".venv" / "bin" / "python"
    py = str(venv_py) if venv_py.exists() else sys.executable
    canon = ENGINE / "rebound_rust" / "notebooks"
    code = (
        "import sys, nbformat\n"
        "from nbclient import NotebookClient\n"
        "p = sys.argv[1]\n"
        "nb = nbformat.read(p, as_version=4)\n"
        "NotebookClient(nb, timeout=1800, kernel_name='python3',\n"
        "               resources={'metadata': {'path': sys.argv[2]}}).execute()\n"
        "nbformat.write(nb, p)\n"
        "print('ok', p)\n"
    )
    r = subprocess.run([py, "-c", code, str(path), str(canon)],
                       capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-4000:])
    nb = json.loads(path.read_text(encoding="utf-8"))
    n = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
    return r.returncode == 0, n


def main() -> int:
    good, bad = [], []
    for arg in sys.argv[1:]:
        path = Path(arg).resolve()
        fam = family_of(path)
        if fam == "posim":
            # rust_* cells invoke cargo, which walks parent dirs for the
            # workspace manifest — run them from the engine root so it finds
            # the real one (a stray file above the copy tree breaks it).
            cwd = (ENGINE if path.parents[1].name == "06_rust_compiled_examples"
                   else path.parent)
            ok, n = run_exec(path, cwd=cwd,
                             env={"POSIM_BIN": str(ENGINE / "target" / "release" / "posim"),
                                  "POSIM_NO_BROWSER": "1"})
        elif fam == "mercury":
            ok, n = run_exec(path, cwd=WS / "planet_Mercury" / "notebook",
                             env={"MERCURY_NO_BROWSER": "1"})
        else:
            ok, n = run_rebound(path)
        (good if ok else bad).append(str(path))
        print(f"{'ok  ' if ok else 'FAIL'} {path.relative_to(WS)}  ({n} cells)")
    print(f"{len(good)} ok, {len(bad)} failed")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
