#!/usr/bin/env python3
"""player_video_double_pendulum_hinges.py — full player for video_double_pendulum_hinges.ipynb (standard library only).

Modes
  python3 player_video_double_pendulum_hinges.py run              re-execute this notebook copy
  python3 player_video_double_pendulum_hinges.py gui              open the browser GUI
  python3 player_video_double_pendulum_hinges.py jump             open the GUI jumped to the key
                                           event (the recorded scene's final state)
  python3 player_video_double_pendulum_hinges.py capture out.png  save a PNG of the GUI
                                           (headless Chrome)
  python3 player_video_double_pendulum_hinges.py data outdir      export this notebook's data

Family: video.  Everything runs against recorded, verified artifacts of
the pure-Rust SUNDIALS engine; the player never computes physics itself.
"""

import http.server
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent
WS = HERE.parents[2]                       # workspace root
ENGINE = WS / "rustSolveIt_macos-silicon_SUNDIALS_7_8_0"
if not ENGINE.exists():
    ENGINE = WS                # fresh clone: the repository root IS the engine
TAG = "video_double_pendulum_hinges"
FAMILY = "video"
GUI_PAGE = None                      # main browser-GUI artifact (or None)
MOVIE = ENGINE / "videos/double_pendulum_hinges.html"                            # recorded movie page (or None)
GUI_SERVER = ENGINE / "gui/double_pendulum_hinges/server.py"                  # live-GUI server.py (or None)
POSIM_SCRIPT = ENGINE / "videos/scenes/double_pendulum_hinges.posim"              # paired .posim scene (or None)
DBS = []                                # SQLite database files
CAPTURE_BUTTON = None          # id of the page's jump button


def chrome():
    for c in ("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
              shutil.which("google-chrome") or "",
              shutil.which("chromium") or ""):
        if c and Path(c).exists():
            return c
    raise SystemExit("headless capture needs Google Chrome installed")


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def serve(directory, page_name, jump=False):
    """Serve `directory` on an ephemeral port; return the page URL.
    With jump=True, serve an augmented copy that clicks the page's own
    jump-to-capture button after load (the original file is untouched)."""
    import functools
    port = free_port()
    workdir = directory
    name = page_name
    if jump and CAPTURE_BUTTON:
        tmp = Path(tempfile.mkdtemp(prefix="player_" + TAG + "_"))
        body = (Path(directory) / page_name).read_text(encoding="utf-8")
        body += ('<script>window.addEventListener("load",function(){'
                 'var b=document.getElementById("' + CAPTURE_BUTTON + '");'
                 'if(b){b.click();}});</script>')
        (tmp / name).write_text(body, encoding="utf-8")
        workdir = tmp
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=str(workdir))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", port), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return f"http://127.0.0.1:{port}/{name}", httpd


def scene_replay_url():
    """Replay the paired .posim scene into a fresh posim child and open its
    live scene window server; returns (url, process)."""
    binary = os.environ.get("POSIM_BIN", str(ENGINE / "target" / "release" / "posim"))
    if not Path(binary).is_file():
        raise SystemExit("build the engine first: cd rustSolveIt_macos-silicon_"
                         "SUNDIALS_7_8_0 && cargo build --release -p posim")
    port = free_port()
    proc = subprocess.Popen([binary, "--machine"],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            text=True, bufsize=1,
                            env=dict(os.environ, POSIM_NO_BROWSER="1"))
    def rpc(cmd):
        proc.stdin.write(json.dumps({"op": "exec", "code": cmd}) + "\n")
        proc.stdin.flush()
        while True:
            line = proc.stdout.readline()
            if not line:
                raise RuntimeError("posim closed the connection")
            r = json.loads(line)
            if "event" in r:          # asynchronous scene notice, not a reply
                continue
            if not r.get("ok"):
                raise RuntimeError("posim refused %r: %s" % (cmd, r.get("error")))
            return r
    # Commands are line-based, but DEF blocks span lines: accumulate until
    # braces balance and send each complete statement as one code string.
    pending = []
    depth = 0
    for raw in Path(POSIM_SCRIPT).read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not pending and (not stripped or stripped.startswith("#")):
            continue
        pending.append(raw)
        depth += raw.count("{") - raw.count("}")
        if depth <= 0:
            rpc("\n".join(pending))
            pending = []
            depth = 0
    if pending:
        rpc("\n".join(pending))
    rpc(f"scene create {port}")
    time.sleep(0.5)
    return f"http://127.0.0.1:{port}/", proc


def open_gui(jump=False):
    if FAMILY in ("solveit", "dynamic", "collision") or (
            FAMILY == "video" and GUI_PAGE is None and MOVIE is None):
        if not POSIM_SCRIPT:
            print("this notebook has no paired .posim scene; its verified")
            print("artifact is the executed notebook itself — open it with:")
            print("  jupyter lab", TAG + ".ipynb")
            return
        url, proc = scene_replay_url()
        print("live posim scene window:", url)
        print("(the paired scene", Path(POSIM_SCRIPT).name,
              "was replayed; Ctrl-C to quit)")
        webbrowser.open(url)
        try:
            proc.wait()
        except KeyboardInterrupt:
            proc.terminate()
        return
    if GUI_SERVER and not jump:
        print("starting the live GUI server (Ctrl-C to quit):", GUI_SERVER)
        subprocess.run([sys.executable, str(GUI_SERVER)], cwd=str(ENGINE))
        return
    page = GUI_PAGE or MOVIE
    if page is None:
        print("this notebook family has no browser GUI; its verified artifact")
        print("is the executed notebook itself — open it with: jupyter lab",
              TAG + ".ipynb")
        return
    url, httpd = serve(str(Path(page).parent), Path(page).name, jump=jump)
    print(("jumped to " + "the recorded scene's final state" + ": ") if jump else "GUI: ", url)
    webbrowser.open(url)
    print("serving; Ctrl-C to quit")
    try:
        while True:
            time.sleep(3600)
    except KeyboardInterrupt:
        httpd.shutdown()


def capture(out_png, jump=True):
    c = chrome()
    if FAMILY in ("solveit", "dynamic", "collision") or (
            FAMILY == "video" and GUI_PAGE is None and MOVIE is None):
        if not POSIM_SCRIPT:
            raise SystemExit("no paired .posim scene recorded for this notebook")
        url, proc = scene_replay_url()
        try:
            subprocess.run([c, "--headless=new", "--disable-gpu",
                            "--window-size=1500,950",
                            "--screenshot=" + str(out_png), url],
                           check=True, capture_output=True, timeout=120)
        finally:
            proc.terminate()
    else:
        page = GUI_PAGE or MOVIE
        if page is None:
            raise SystemExit("no browser GUI to capture for this notebook")
        url, httpd = serve(str(Path(page).parent), Path(page).name, jump=jump)
        time.sleep(0.3)
        subprocess.run([c, "--headless=new", "--disable-gpu",
                        "--window-size=1500,950", "--virtual-time-budget=8000",
                        "--screenshot=" + str(out_png), url],
                       check=True, capture_output=True, timeout=120)
        httpd.shutdown()
    print("captured:", out_png)


def export_data(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    for db in DBS:
        import csv
        import sqlite3
        con = sqlite3.connect(db)
        cur = con.cursor()
        tables = [r[0] for r in cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
        for t in tables:
            cur.execute(f"SELECT * FROM {t}")
            cols = [d[0] for d in cur.description]
            rows = cur.fetchall()
            with open(out / f"{Path(db).stem}_{t}.csv", "w", newline="",
                      encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(cols)
                w.writerows(rows)
        con.close()
        print("exported", len(tables), "tables from", Path(db).name)
    nb = json.loads((HERE / (TAG + ".ipynb")).read_text(encoding="utf-8"))
    with open(out / "transcript.txt", "w", encoding="utf-8") as f:
        for cell in nb["cells"]:
            if cell["cell_type"] != "code":
                continue
            f.write("### In:\n" + "".join(cell["source"]) + "\n### Out:\n")
            for o in cell.get("outputs", []):
                f.write("".join(o.get("text", [])))
            f.write("\n" + "=" * 70 + "\n")
    print("exported executed transcript ->", out / "transcript.txt")
    if MOVIE:
        shutil.copyfile(MOVIE, out / Path(MOVIE).name)
        print("copied movie page ->", out / Path(MOVIE).name)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "gui"
    if mode == "run":
        sys.exit(subprocess.run([sys.executable,
                                 str(WS / "SolveIt_Notebooks_for_rust" / "_tools" / "run_copy.py"),
                                 str(HERE / (TAG + ".ipynb"))]).returncode)
    elif mode == "gui":
        open_gui(jump=False)
    elif mode == "jump":
        open_gui(jump=True)
    elif mode == "capture":
        capture(sys.argv[2] if len(sys.argv) > 2 else TAG + ".png")
    elif mode == "data":
        export_data(sys.argv[2] if len(sys.argv) > 2 else "exported_data")
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
