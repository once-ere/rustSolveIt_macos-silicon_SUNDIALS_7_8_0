#!/usr/bin/env python3
"""run_notebooks.py - execute every notebook in this folder and report.

Each notebook is run top to bottom in a fresh kernel, exactly as a reader
would experience it by choosing "Run All". The executed copy (with its
outputs) is written back in place, so the notebooks in the repository are
the ones that actually ran.

Usage:
    python run_notebooks.py                # run all of them
    python run_notebooks.py shearing_sheet # run only matching ones
    python run_notebooks.py --timeout 3600 # per-cell timeout, seconds

Part of the rebound_rs documentation toolchain, GPL-3.0-or-later.
"""

import io
import os
import sys
import time

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    args = [a for a in sys.argv[1:]]

    timeout = 1800
    if "--timeout" in args:
        k = args.index("--timeout")
        timeout = int(args[k + 1])
        del args[k:k + 2]

    names = sorted(f for f in os.listdir(here) if f.endswith(".ipynb"))
    if args:
        names = [f for f in names if any(a in f for a in args)]

    if not names:
        print("no notebooks matched")
        return 1

    print("Executing %d notebook(s), per-cell timeout %ds\n" % (len(names), timeout))

    failures = []
    for name in names:
        path = os.path.join(here, name)
        nb = nbformat.read(path, as_version=4)
        client = NotebookClient(
            nb,
            timeout=timeout,
            kernel_name="python3",
            resources={"metadata": {"path": here}},
            allow_errors=False,
        )
        t0 = time.time()
        sys.stdout.write("  %-34s " % name)
        sys.stdout.flush()
        try:
            client.execute()
            nbformat.write(nb, path)
            print("ok      (%5.1fs)" % (time.time() - t0))
        except CellExecutionError as e:
            print("FAILED  (%5.1fs)" % (time.time() - t0))
            failures.append((name, str(e).strip().split("\n")[-1]))
            nbformat.write(nb, path)
        except Exception as e:  # timeout, dead kernel, ...
            print("ERROR   (%5.1fs)  %s" % (time.time() - t0, type(e).__name__))
            failures.append((name, "%s: %s" % (type(e).__name__, e)))

    print("")
    if failures:
        print("%d of %d notebooks failed:" % (len(failures), len(names)))
        for name, msg in failures:
            print("  %-34s %s" % (name, msg[:110]))
        return 1

    print("All %d notebooks executed with no errors." % len(names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
