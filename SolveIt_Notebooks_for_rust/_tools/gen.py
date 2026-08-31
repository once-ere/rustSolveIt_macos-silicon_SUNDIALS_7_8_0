#!/usr/bin/env python3
"""gen.py — build the SolveIt_Notebooks_for_rust encyclopedia tree.

For every canonical project-authored Jupyter notebook (128 of them):
  1. copy it into  SolveIt_Notebooks_for_rust/<topic>/<tag>/<tag>.ipynb
  2. prepend ONE header markdown cell: exact execute commands, exact
     browser-GUI display commands, movie access, database access, and the
     complete companion player script (inline, verbatim);
  3. write the companion  player_<tag>.py  next to it.

Canonical sources (duplicates elsewhere are deliberately excluded):
  - rustSolveIt_macos-silicon_SUNDIALS_7_8_0/notebooks/*.ipynb        (109)
  - rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rebound_rust/notebooks/  (17)
  - planet_Mercury/notebook/*.ipynb                                   (2)

Standard library only. Deterministic. Usage:
    python3 gen.py            # generate everything
    python3 gen.py --only PAT # only tags containing PAT (pilot runs)
"""

import json
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[2]          # workspace root
ENGINE = WS / "rustSolveIt_macos-silicon_SUNDIALS_7_8_0"
ROOT = WS / "SolveIt_Notebooks_for_rust"
MANIFEST = ROOT / "_tools" / "manifest.json"

TOPICS = {
    "mercury": "01_planet_mercury_tidal_locking",
    "solveit": "02_solveit_worked_examples",
    "dynamic": "03_dynamics_and_routh",
    "collision": "04_collisions",
    "video": "05_mechanism_videos",
    "rust": "06_rust_compiled_examples",
    "rebound": "07_nbody_rebound_rust",
}


def family_of(src: Path) -> str:
    if "planet_Mercury" in str(src):
        return "mercury"
    if "rebound_rust" in str(src):
        return "rebound"
    stem = src.stem
    for fam in ("solveit", "dynamic", "collision", "video", "rust"):
        if stem.startswith(fam + "_"):
            return fam
    raise SystemExit(f"unclassifiable notebook: {src}")


def collect():
    out = []
    for p in sorted((ENGINE / "notebooks").glob("*.ipynb")):
        out.append(p)
    for p in sorted((ENGINE / "rebound_rust" / "notebooks").glob("*.ipynb")):
        out.append(p)
    for p in sorted((WS / "planet_Mercury" / "notebook").glob("*.ipynb")):
        out.append(p)
    return out


def paired_assets(fam: str, tag: str):
    """Return dict of the notebook's paired scripts / GUIs / movies / data."""
    a = {"script": None, "movie": None, "gui_server": None, "gui_pages": [],
         "db": [], "scene_replay": False}
    if fam == "solveit":
        s = ENGINE / "scripts" / "solveit" / (tag[len("solveit_"):] + ".posim")
        a["script"] = s if s.exists() else None
        a["scene_replay"] = a["script"] is not None
    elif fam == "collision":
        s = ENGINE / "scripts" / "collisions" / (tag[len("collision_"):] + ".posim")
        a["script"] = s if s.exists() else None
        a["scene_replay"] = a["script"] is not None
    elif fam == "dynamic":
        s = ENGINE / "dynamic_notebooks" / (tag[len("dynamic_"):] + ".posim")
        a["script"] = s if s.exists() else None
        a["scene_replay"] = a["script"] is not None
    elif fam == "video":
        name = tag[len("video_"):]
        s = ENGINE / "videos" / "scenes" / (name + ".posim")
        a["script"] = s if s.exists() else None
        m = ENGINE / "videos" / (name + ".html")
        a["movie"] = m if m.exists() else None
        g = ENGINE / "gui" / name / "server.py"
        a["gui_server"] = g if g.exists() else None
        a["scene_replay"] = a["script"] is not None
    elif fam == "rust":
        s = ENGINE / "physical_object" / "examples" / (tag[len("rust_"):] + ".rs")
        a["script"] = s if s.exists() else None
    elif fam == "mercury":
        if tag == "mercury_tidal_locking":
            a["gui_pages"] = [WS / "planet_Mercury" / "gui" / "mercury_orbit.html"]
            a["db"] = [WS / "planet_Mercury" / "data" / "mercury_orbit.sqlite3"]
        else:
            a["gui_pages"] = [WS / "planet_Mercury" / "gui" / "mercury_test2.html"]
            a["db"] = [WS / "planet_Mercury" / "data" / "mercury_test2.sqlite3"]
    return a


# --------------------------------------------------------------------------
# The companion player. @TOKENS@ are substituted per notebook; the SAME text
# (after substitution) is embedded in the header cell, verbatim.
# --------------------------------------------------------------------------
PLAYER_TEMPLATE = r'''#!/usr/bin/env python3
"""player_@TAG@.py — full player for @TAG@.ipynb (standard library only).

Modes
  python3 player_@TAG@.py run              re-execute this notebook copy
  python3 player_@TAG@.py gui              open the browser GUI
  python3 player_@TAG@.py jump             open the GUI jumped to the key
                                           event (@KEY_EVENT@)
  python3 player_@TAG@.py capture out.png  save a PNG of the GUI
                                           (headless Chrome)
  python3 player_@TAG@.py data outdir      export this notebook's data

Family: @FAMILY@.  Everything runs against recorded, verified artifacts of
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
TAG = "@TAG@"
FAMILY = "@FAMILY@"
GUI_PAGE = @GUI_PAGE@                      # main browser-GUI artifact (or None)
MOVIE = @MOVIE@                            # recorded movie page (or None)
GUI_SERVER = @GUI_SERVER@                  # live-GUI server.py (or None)
POSIM_SCRIPT = @POSIM_SCRIPT@              # paired .posim scene (or None)
DBS = @DBS@                                # SQLite database files
CAPTURE_BUTTON = @CAPTURE_BUTTON@          # id of the page's jump button


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
    print(("jumped to " + "@KEY_EVENT@" + ": ") if jump else "GUI: ", url)
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
'''


def build_player(fam, tag, assets):
    def lit(p):
        if p is None:
            return "None"
        r = str(Path(p).relative_to(WS)).replace("\\", "/")
        # Engine-relative paths go through ENGINE so a fresh clone (where
        # the repo root IS the engine) resolves them too.
        prefix = "rustSolveIt_macos-silicon_SUNDIALS_7_8_0/"
        if r.startswith(prefix):
            return 'ENGINE / "' + r[len(prefix):] + '"'
        return 'WS / "' + r + '"'
    key_event = {
        "mercury": "the 3:2 resonance capture moment",
        "video": "the recorded scene's final state",
        "solveit": "the worked example's final state",
        "dynamic": "the dynamic scene's final state",
        "collision": "the collision outcome",
        "rust": "the compiled example's SUCCESS verdict",
        "rebound": "the N-body run's final state",
    }[fam]
    txt = PLAYER_TEMPLATE
    txt = txt.replace("@TAG@", tag).replace("@FAMILY@", fam)
    txt = txt.replace("@KEY_EVENT@", key_event)
    txt = txt.replace("@GUI_PAGE@", lit(assets["gui_pages"][0]) if assets["gui_pages"] else "None")
    txt = txt.replace("@MOVIE@", lit(assets["movie"]))
    txt = txt.replace("@GUI_SERVER@", lit(assets["gui_server"]))
    txt = txt.replace("@POSIM_SCRIPT@", lit(assets["script"] if assets["scene_replay"] else None))
    txt = txt.replace("@DBS@", "[" + ", ".join(lit(d) for d in assets["db"]) + "]")
    txt = txt.replace("@CAPTURE_BUTTON@", '"bt-cap"' if fam == "mercury" else "None")
    return txt


def rel(p):
    return str(Path(p).relative_to(WS)) if p else None


def header_cell(fam, tag, topic, assets, player_text):
    e = "rustSolveIt_macos-silicon_SUNDIALS_7_8_0"
    L = []
    L.append(f"# Encyclopedia header — `{tag}` ({topic})")
    L.append("")
    L.append("This is the **SolveIt_Notebooks_for_rust encyclopedia copy** of this "
             "notebook. Everything below is exact and copy-paste runnable from the "
             "**workspace root** (the folder that contains both "
             f"`{e}/` and `SolveIt_Notebooks_for_rust/`).")
    L.append("")
    L.append("## 1. Execute this notebook")
    L.append("")
    L.append("One-time build of the pure-Rust engine behind it:")
    L.append("```bash")
    if fam in ("solveit", "dynamic", "collision", "video", "rust"):
        L.append(f"cd {e} && cargo build --release -p posim")
    elif fam == "rebound":
        L.append(f"cd {e}/reboundx_rust && cargo build --release")
    else:
        L.append("cd planet_Mercury/mercury_rs && cargo build --release 2>&1 | tee /tmp/mercury_build.log")
    L.append("```")
    L.append("Batch execution of this copy (writes real outputs back into it):")
    L.append("```bash")
    L.append(f"python3 SolveIt_Notebooks_for_rust/_tools/run_copy.py "
             f"SolveIt_Notebooks_for_rust/{topic}/{tag}/{tag}.ipynb")
    L.append("```")
    L.append("Interactive execution (Shift+Enter through every cell, "
             "Python 3 (ipykernel) kernel):")
    L.append("```bash")
    if fam == "mercury":
        L.append("python3 -m venv planet_Mercury/notebook/.venv  "
                 "# once; then: planet_Mercury/notebook/.venv/bin/pip install jupyterlab")
        L.append(f"cd SolveIt_Notebooks_for_rust/{topic}/{tag} && "
                 f"../../../planet_Mercury/notebook/.venv/bin/jupyter lab {tag}.ipynb")
    else:
        L.append("python3 -m pip install --user jupyterlab   # once")
        L.append(f"cd SolveIt_Notebooks_for_rust/{topic}/{tag} && jupyter lab {tag}.ipynb")
    L.append("```")
    if fam in ("solveit", "dynamic", "collision", "video"):
        L.append(f"(The notebook finds the simulator itself; if needed set "
                 f"`POSIM_BIN=$PWD/{e}/target/release/posim`.)")
    L.append("")
    L.append("## 2. Display the browser GUI")
    L.append("")
    if fam == "mercury":
        page = rel(assets["gui_pages"][0])
        L.append("This notebook bakes a self-contained animated player page and its "
                 "own database; after one execution:")
        L.append("```bash")
        L.append(f"open {page}")
        L.append("```")
        L.append("or, with the companion player (serves it and can auto-jump):")
        L.append("```bash")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py gui")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py jump   # straight to the capture moment")
        L.append("```")
    elif fam == "video":
        if assets["movie"]:
            L.append("The recorded **movie page** (self-contained; double-click also works):")
            L.append("```bash")
            L.append(f"open {rel(assets['movie'])}")
            L.append("```")
        if assets["gui_server"]:
            L.append("The **live GUI** (a posim child integrating in real time behind a "
                     "canvas page with Start/Pause/Reset and live physics readouts):")
            L.append("```bash")
            L.append(f"python3 {rel(assets['gui_server'])}   # prints its URL, Ctrl-C to quit")
            L.append("```")
    elif fam in ("solveit", "dynamic", "collision"):
        L.append("This notebook's GUI is posim's own live **scene window**. Replay the "
                 "paired scene and open it:")
        L.append("```bash")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py gui")
        L.append("```")
        if assets["script"]:
            L.append("Manually, in the interactive REPL "
                     f"(`{e}/target/release/posim`): type the commands of "
                     f"`{rel(assets['script'])}`, then `SCENE CREATE 8900` and open "
                     "http://127.0.0.1:8900/ in a browser.")
    elif fam == "rust":
        L.append("This family pairs with a **compiled, self-checking Rust example** — "
                 "its verdict is textual (SUCCESS), and the notebook itself is the "
                 "display artifact. Run the example directly:")
        L.append("```bash")
        L.append(f"cd {e} && cargo run --release -p physical_object --example {tag[len('rust_'):]}")
        L.append("```")
    else:
        L.append("This family's artifacts are the executed notebook and the files it "
                 "writes (see section 4); where a browser GUI exists the companion "
                 "player opens it:")
        L.append("```bash")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py gui")
        L.append("```")
    L.append("")
    L.append("## 3. Movies")
    L.append("")
    if fam == "video":
        L.append(f"The movie already exists: `{rel(assets['movie'])}` is the recorded "
                 "browser video page for this exact scene (open it with the command "
                 "above). To re-record it: run the live GUI server (section 2) and "
                 "use any screen recorder (macOS: Cmd+Shift+5), or re-bake the scene "
                 "page with the engine's own tooling in `videos/`.")
    elif fam == "mercury":
        L.append("The player page IS the movie: press Play (Space) and it replays "
                 "the recorded CVODE samples — braking, capture, and lock — with a "
                 "scrub bar and speed control. To export a conventional video file, "
                 "screen-record the playing page (macOS: Cmd+Shift+5). The "
                 "`Jump to capture` button (or `player jump`) goes straight to the "
                 "capture moment.")
    else:
        L.append("No pre-recorded movie pairs with this notebook. To watch it as a "
                 "movie: open the live scene window (section 2) — the scene plays in "
                 "real time — and screen-record it (macOS: Cmd+Shift+5). The 13 "
                 "`video_*` notebooks in `05_mechanism_videos/` are the family with "
                 "recorded movie pages.")
    L.append("")
    L.append("## 4. Database / data access")
    L.append("")
    if fam == "mercury":
        db = rel(assets["db"][0])
        L.append(f"Executing the notebook builds `{db}` (SQLite, schema documented "
                 "inside the notebook's own section 7: run, sample, event, branch, "
                 "target" + (", run_extra" if "test2" in tag else "") +
                 " tables). Open it:")
        L.append("```bash")
        L.append(f"sqlite3 {db} \".tables\"")
        L.append(f"sqlite3 {db} \"SELECT * FROM run;\"")
        L.append("```")
        L.append("Or export every table to CSV plus the executed transcript:")
        L.append("```bash")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py data exported_data")
        L.append("```")
    else:
        L.append("This notebook's data is its **executed transcript** (every command "
                 "sent to the pure-Rust engine and every reply, embedded as real "
                 "cell outputs) — there is no SQLite database in this family. "
                 "Export the transcript (and any paired artifacts) with:")
        L.append("```bash")
        L.append(f"python3 SolveIt_Notebooks_for_rust/{topic}/{tag}/player_{tag}.py data exported_data")
        L.append("```")
        if assets["script"]:
            L.append(f"The paired source scene/script is `{rel(assets['script'])}`.")
    L.append("")
    L.append("## 5. The full player (companion script, shipped in this folder)")
    L.append("")
    L.append(f"`player_{tag}.py` sits next to this notebook. It re-executes the "
             "notebook (`run`), opens the browser GUI (`gui`), **jumps to the key "
             "event** (`jump`), captures a PNG of the GUI headlessly (`capture "
             "out.png`, needs Google Chrome), and exports the data (`data outdir`). "
             "It is standard-library Python and never computes physics — it drives "
             "and displays the recorded, verified artifacts. The complete script, "
             "verbatim:")
    L.append("")
    L.append("```python")
    L.append(player_text.rstrip("\n"))
    L.append("```")
    src = "\n".join(L)
    return {"cell_type": "markdown",
            "metadata": {"tags": ["encyclopedia-header"]},
            "source": src.splitlines(keepends=True)}


def main():
    only = None
    refresh = "--refresh" in sys.argv
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    entries = []
    made = 0
    for src in collect():
        fam = family_of(src)
        tag = src.stem
        if only and only not in tag:
            continue
        topic = TOPICS[fam]
        dest_dir = ROOT / topic / tag
        dest_dir.mkdir(parents=True, exist_ok=True)
        assets = paired_assets(fam, tag)
        player_text = build_player(fam, tag, assets)
        (dest_dir / f"player_{tag}.py").write_text(player_text, encoding="utf-8")
        # --refresh: rebuild player + header IN the existing executed copy,
        # preserving every executed output cell.
        read_from = dest_dir / f"{tag}.ipynb" if (
            refresh and (dest_dir / f"{tag}.ipynb").exists()) else src
        nb = json.loads(read_from.read_text(encoding="utf-8"))
        nb["cells"] = ([header_cell(fam, tag, topic, assets, player_text)]
                       + [c for c in nb["cells"]
                          if "encyclopedia-header" not in
                          c.get("metadata", {}).get("tags", [])])
        ids = {c.get("id") for c in nb["cells"] if c.get("id")}
        if ids:
            nb["cells"][0]["id"] = "cell-enc-header"
        (dest_dir / f"{tag}.ipynb").write_text(
            json.dumps(nb, indent=1) + "\n", encoding="utf-8")
        entries.append({
            "tag": tag, "family": fam, "topic": topic,
            "source": rel(src),
            "copy": str((dest_dir / (tag + ".ipynb")).relative_to(WS)),
            "player": str((dest_dir / f"player_{tag}.py").relative_to(WS)),
            "script": rel(assets["script"]), "movie": rel(assets["movie"]),
            "gui_server": rel(assets["gui_server"]),
            "gui_pages": [rel(p) for p in assets["gui_pages"]],
            "db": [rel(d) for d in assets["db"]],
            "scene_replay": assets["scene_replay"],
            "executed": False, "image": None,
        })
        made += 1
    if only is None:
        MANIFEST.write_text(json.dumps(entries, indent=1) + "\n", encoding="utf-8")
        print(f"manifest: {MANIFEST} ({len(entries)} entries)")
    print(f"generated {made} notebook folders under {ROOT}")


if __name__ == "__main__":
    main()
