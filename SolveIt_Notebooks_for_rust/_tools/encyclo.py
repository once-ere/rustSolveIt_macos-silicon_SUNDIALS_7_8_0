#!/usr/bin/env python3
"""encyclo.py — build the encyclopedia (Markdown + HTML + PDF).

Reads manifest.json, the executed notebook copies (for each notebook's own
opening description), the execution logs (for verdicts), and the captured
GUI images; emits:

    SolveIt_Notebooks_for_rust/ENCYCLOPEDIA.md
    SolveIt_Notebooks_for_rust/ENCYCLOPEDIA.html
    SolveIt_Notebooks_for_rust/ENCYCLOPEDIA.pdf   (headless Chrome print)

Standard library only (Chrome only for the PDF). Deterministic ordering.
"""

import html
import json
import re
import subprocess
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[2]
ENGINE = WS / "rustSolveIt_macos-silicon_SUNDIALS_7_8_0"
ROOT = WS / "SolveIt_Notebooks_for_rust"
MANIFEST = ROOT / "_tools" / "manifest.json"

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

TOPIC_INTROS = {
    "01_planet_mercury_tidal_locking": (
        "Planet Mercury's 3:2 spin-orbit capture, computed end to end by the "
        "pure-Rust SUNDIALS 7.8.0 CVODE solver: test 1 is the pure two-body "
        "story (tidal braking, resonance capture, lock), test 2 adds "
        "Einstein's general-relativistic perihelion precession and Jupiter's "
        "Laplace-Lagrange secular forcing — the lock follows the precessing "
        "ellipse. Each notebook builds its own SQLite database and bakes an "
        "animated browser player with a Jump-to-capture button."),
    "02_solveit_worked_examples": (
        "The SolveIt worked examples: each notebook drives the posim "
        "simulator (pure-Rust SUNDIALS backend) through one classic physics "
        "problem with closed-form answers checked in the transcript."),
    "03_dynamics_and_routh": (
        "The dynamic notebooks — free rigid-body motion, charged particles "
        "in fields, orbital mechanics, and the 34 Routh rigid-body problems "
        "— each a machine-mode posim session with its analytic checks."),
    "04_collisions": (
        "The collision-detection walkthroughs: elastic exchanges, "
        "restitution ladders, billiards, dumbbells, and the box-of-shapes "
        "stress test, all integrated by the pure-Rust engine."),
    "05_mechanism_videos": (
        "The thirteen recorded mechanism scenes: each notebook pairs with a "
        "recorded browser movie page in videos/ AND a live GUI server in "
        "gui/ that integrates the same scene in real time behind canvas "
        "controls with live physics readouts."),
    "06_rust_compiled_examples": (
        "The compiled, self-checking Rust examples: each notebook builds and "
        "runs one physical_object example and asserts its SUCCESS verdict — "
        "the verdict is textual, so these entries carry no GUI image, "
        "honestly."),
    "07_nbody_rebound_rust": (
        "The N-body verification notebooks of the rebound/reboundx pure-Rust "
        "port: integrator equivalence, simulation archives, tides and spin, "
        "shearing sheets. Verdicts are textual (no browser GUI in this "
        "family)."),
}


def first_summary(copy_path: Path) -> str:
    """First original markdown paragraph of the notebook (past our header)."""
    nb = json.loads((WS / copy_path).read_text(encoding="utf-8"))
    for cell in nb["cells"]:
        if cell["cell_type"] != "markdown":
            continue
        if "encyclopedia-header" in cell.get("metadata", {}).get("tags", []):
            continue
        text = "".join(cell["source"])
        text = re.sub(r"^#.*$", "", text, flags=re.M)          # strip headings
        text = re.sub(r"`([^`]*)`", r"\1", text)
        text = re.sub(r"\*\*([^*]*)\*\*", r"\1", text)
        text = re.sub(r"\|.*\|", "", text)                      # strip tables
        paras = [p.strip().replace("\n", " ")
                 for p in text.split("\n\n") if len(p.strip()) > 60]
        if paras:
            p = paras[0]
            # A paragraph ending in ':' introduces an enumerated list —
            # pull the list in so the summary does not dangle mid-thought.
            if p.endswith(":") and len(paras) > 1:
                p = p + " " + paras[1]
            return (p[:700] + "…") if len(p) > 700 else p
    return "(no descriptive markdown found)"


# Editorial corrections where a notebook's own intro prose carries claims
# from an earlier (Windows/MSVC) lineage that THIS machine's embedded run
# contradicts — the quoted summary stays faithful, the note sets the record.
EDITORIAL_NOTES = {
    "libm_diff": (
        "Encyclopedia note: the intro above is the notebook's own "
        "Windows-lineage description; the executed run embedded in THIS "
        "copy reports 'functions that differ: NONE - all bit-identical' "
        "on macOS/Apple Silicon."),
    "bs_pow_diff": (
        "Encyclopedia note: the intro above describes the Windows/MSVC "
        "measurement; the executed run embedded in THIS copy reports zero "
        "pow divergences ('ULP distribution: none') on macOS/Apple "
        "Silicon."),
    "integrators_test": (
        "Encyclopedia note: the bit-identical-vs-MSVC claim is history "
        "from the port's Windows lineage; the verdict embedded in THIS "
        "copy is the macOS/Apple Silicon run's own."),
}


def verdicts():
    """tag -> ('ok'|'FAIL', cells) from the execution logs."""
    v = {}
    for log in sorted((ROOT / "_tools").glob("exec_*.log")):
        for line in log.read_text(encoding="utf-8").splitlines():
            m = re.match(r"(ok|FAIL)\s+.*/([A-Za-z0-9_.]+)\.ipynb\s+\((\d+) cells\)",
                         line.strip())
            if m:
                v[m.group(2)] = (m.group(1), int(m.group(3)))
    return v


def entry_md(e, v):
    tag = e["tag"]
    verdict, cells = v.get(tag, ("NOT RUN", 0))
    L = [f"### `{tag}`", ""]
    L.append(first_summary(e["copy"]))
    if tag in EDITORIAL_NOTES:
        L.append("")
        L.append("*[" + EDITORIAL_NOTES[tag] + "]*")
    L.append("")
    L.append(f"- **Executed**: {'PASS — every cell green' if verdict == 'ok' else verdict}"
             f" ({cells} code cells)")
    L.append(f"- **Copy**: `{e['copy']}` (header cell prepended: run commands, "
             f"GUI commands, movies, data access, full player script)")
    L.append(f"- **Canonical source**: `{e['source']}`")
    if e.get("script"):
        L.append(f"- **Paired scene/source**: `{e['script']}`")
    if e.get("movie"):
        L.append(f"- **Movie page**: `{e['movie']}`")
    if e.get("gui_server"):
        L.append(f"- **Live GUI**: `python3 {e['gui_server']}`")
    for db in e.get("db", []):
        L.append(f"- **Database**: `{db}` (SQLite)")
    L.append(f"- **Player**: `python3 {e['player']} "
             f"run|gui|jump|capture out.png|data outdir`")
    L.append("")
    if e.get("image"):
        img = e["image"].split("SolveIt_Notebooks_for_rust/")[-1]
        L.append(f"![{tag} GUI]({img})")
    else:
        L.append("*(no browser GUI in this family — the verdict is textual, "
                 "embedded in the executed notebook)*")
    L.append("")
    return "\n".join(L)


def upstream_appendix():
    names = set()
    for d in [ENGINE / "rebound" / "rebound" / "ipython_examples",
              ENGINE / "reboundx" / "ipython_examples"]:
        if d.exists():
            names.update(p.name for p in d.glob("*.ipynb"))
    return sorted(names)


def build_md(entries, v):
    topics = {}
    for e in entries:
        topics.setdefault(e["topic"], []).append(e)
    n_ok = sum(1 for e in entries if v.get(e["tag"], ("", 0))[0] == "ok")
    L = []
    L.append("# The SolveIt Notebooks for Rust — an encyclopedia")
    L.append("")
    L.append(f"**{len(entries)} Jupyter notebooks**, every one project-authored, "
             f"every one executed for real against the pure-Rust SUNDIALS 7.8.0 "
             f"engine on this machine ({n_ok}/{len(entries)} green), organized "
             "into seven major topics. Each notebook lives in its own tagged "
             "subfolder together with its companion player script "
             "(`player_<tag>.py`), and each carries a prepended header cell "
             "with the exact commands to execute it, display its browser GUI, "
             "create/access its movies, and access its data — plus the full "
             "player script inline.")
    L.append("")
    L.append("Every simulation below was integrated by the vendored pure-Rust "
             "port of SUNDIALS 7.8.0 CVODE (BDF + Newton + dense solver) — no "
             "C, no FFI, no unsafe code — and every notebook checks its own "
             "physics against closed forms before it is allowed to say ok.")
    L.append("")
    L.append("## Contents")
    L.append("")
    for t in sorted(topics):
        L.append(f"- **{t}** — {len(topics[t])} notebooks")
    L.append("- Appendix A — excluded duplicates")
    L.append("- Appendix B — vendored upstream reference notebooks")
    L.append("- Appendix C — reproduction commands")
    L.append("")
    for t in sorted(topics):
        L.append(f"## {t}")
        L.append("")
        L.append(TOPIC_INTROS.get(t, ""))
        L.append("")
        for e in sorted(topics[t], key=lambda x: x["tag"]):
            L.append(entry_md(e, v))
    L.append("## Appendix A — excluded duplicates")
    L.append("")
    L.append("The workspace holds several parallel copies of the same notebooks; "
             "the encyclopedia uses one canonical edition of each and excludes, "
             "deliberately: the 109 earlier editions under "
             "`rustSolveIt_Using_SUNDIALS_7_8_0/version-7.8.0/notebooks/` (the "
             "pre-macOS lineage of the same notebooks), the two engine-mirror "
             "copies of the planet-Mercury notebooks under "
             "`rustSolveIt_macos-silicon_SUNDIALS_7_8_0/planet_Mercury/notebook/`, "
             "two strays at the workspace root (one byte-identical to its "
             "canonical copy; the other an earlier edition of the test-2 "
             "notebook superseded by the canonical one), and the "
             "duplicated `docs/ipython_examples` copy inside the vendored "
             "rebound sources. `StageNbks/*.posim` are posim REPL sessions, "
             "not Jupyter notebooks, and are documented in the engine's own "
             "STAGENBKS_PROVENANCE.md.")
    L.append("")
    L.append("## Appendix B — vendored upstream reference notebooks (not executed)")
    L.append("")
    L.append("The vendored third-party rebound / REBOUNDx sources ship their own "
             "example notebooks (they exercise the upstream C libraries, not "
             "this project's Rust engine, and so are catalogued here as "
             "reference only):")
    L.append("")
    for n in upstream_appendix():
        L.append(f"- `{n}`")
    L.append("")
    L.append("## Appendix C — reproduction commands")
    L.append("")
    L.append("```bash")
    L.append("# from the workspace root:")
    L.append("cd rustSolveIt_macos-silicon_SUNDIALS_7_8_0 && cargo build --release -p posim && cd ..")
    L.append("cd planet_Mercury/mercury_rs && cargo build --release && cd ../..")
    L.append("python3 SolveIt_Notebooks_for_rust/_tools/gen.py       # regenerate copies+players")
    L.append("python3 SolveIt_Notebooks_for_rust/_tools/run_copy.py SolveIt_Notebooks_for_rust/<topic>/<tag>/<tag>.ipynb")
    L.append("python3 SolveIt_Notebooks_for_rust/_tools/shoot.py     # re-capture every GUI image")
    L.append("python3 SolveIt_Notebooks_for_rust/_tools/encyclo.py   # rebuild this encyclopedia")
    L.append("```")
    L.append("")
    return "\n".join(L)


def md_to_html_body(md_text):
    """A small, honest renderer for exactly the Markdown this file emits."""
    out = []
    in_code = False
    in_list = False
    for line in md_text.splitlines():
        if line.startswith("```"):
            if in_list:
                out.append("</ul>")
                in_list = False
            out.append("<pre>" if not in_code else "</pre>")
            in_code = not in_code
            continue
        if in_code:
            out.append(html.escape(line))
            continue
        img = re.match(r"!\[(.*)\]\((.*)\)", line)
        if img:
            out.append(f'<img src="{img.group(2)}" alt="{html.escape(img.group(1))}">')
            continue
        esc = html.escape(line)
        esc = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", esc)
        esc = re.sub(r"`([^`]+)`", r"<code>\1</code>", esc)
        esc = re.sub(r"\*\(([^)]*)\)\*", r"<i>(\1)</i>", esc)
        if line.startswith("- "):
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append("<li>" + esc[2:] + "</li>")
            continue
        if in_list:
            out.append("</ul>")
            in_list = False
        if line.startswith("### "):
            out.append("<h3>" + esc[4:] + "</h3>")
        elif line.startswith("## "):
            out.append("<h2>" + esc[3:] + "</h2>")
        elif line.startswith("# "):
            out.append("<h1>" + esc[2:] + "</h1>")
        elif line.strip() == "":
            out.append("")
        else:
            out.append("<p>" + esc + "</p>")
    if in_list:
        out.append("</ul>")
    return "\n".join(out)


def main():
    entries = json.loads(MANIFEST.read_text(encoding="utf-8"))
    v = verdicts()
    md = build_md(entries, v)
    (ROOT / "ENCYCLOPEDIA.md").write_text(md, encoding="utf-8")
    body = md_to_html_body(md)
    html_doc = (
        "<!doctype html><meta charset='utf-8'>"
        "<title>The SolveIt Notebooks for Rust</title><style>"
        "body{font:13px/1.55 -apple-system,'Helvetica Neue',sans-serif;"
        "max-width:900px;margin:24px auto;padding:0 16px;color:#1a1e26}"
        "h1{font-size:26px;border-bottom:3px solid #2b5aa0;padding-bottom:8px}"
        "h2{font-size:20px;color:#2b5aa0;border-bottom:1px solid #ccd;"
        "padding-bottom:4px;margin-top:34px;page-break-before:always}"
        "h3{font-size:15px;margin:22px 0 6px;color:#333}"
        "code{background:#eef1f6;padding:1px 4px;border-radius:3px;"
        "font:11.5px ui-monospace,Menlo,monospace}"
        "pre{background:#f4f6fa;border:1px solid #dde;padding:10px;"
        "border-radius:6px;font:11.5px ui-monospace,Menlo,monospace;"
        "overflow-x:auto}"
        "img{max-width:100%;border:1px solid #ccd;border-radius:6px;"
        "margin:6px 0;page-break-inside:avoid}"
        "ul{margin:6px 0}li{margin:2px 0}"
        "</style>" + body)
    (ROOT / "ENCYCLOPEDIA.html").write_text(html_doc, encoding="utf-8")
    r = subprocess.run(
        [CHROME, "--headless=new", "--disable-gpu",
         "--print-to-pdf=" + str(ROOT / "ENCYCLOPEDIA.pdf"),
         "--no-pdf-header-footer",
         (ROOT / "ENCYCLOPEDIA.html").as_uri()],
        capture_output=True, text=True, timeout=600)
    ok = (ROOT / "ENCYCLOPEDIA.pdf").exists() and \
        (ROOT / "ENCYCLOPEDIA.pdf").stat().st_size > 10000
    print("ENCYCLOPEDIA.md :", len(md), "chars")
    print("ENCYCLOPEDIA.html:", len(html_doc), "chars")
    print("ENCYCLOPEDIA.pdf:",
          (ROOT / "ENCYCLOPEDIA.pdf").stat().st_size if ok else "FAILED",
          "" if ok else r.stderr[-400:])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
