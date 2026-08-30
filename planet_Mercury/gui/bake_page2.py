#!/usr/bin/env python3
"""Bake gui/mercury_test2.html from data/mercury_test2.sqlite3.

Test 2's display page: the spin-ratio descent into the 3:2 lock, Jupiter's
breathing eccentricity, and the perihelion odometer, with the headline
numbers on top. Standard library only, fully deterministic (no timestamps,
fixed float formatting), self-contained output (inline CSS + SVG, no
external resources). Run from the planet_Mercury folder or the repo root:

    python3 gui/bake_page2.py
"""

import math
import sqlite3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
DB = BASE / "data" / "mercury_test2.sqlite3"
OUT = HERE / "mercury_test2.html"

MYR = 3.15576e13
ARCSEC_CY = 206264.80624709636 * 3.15576e9


def fetch_series(cur, expr):
    n = cur.execute(
        "SELECT COUNT(*) FROM sample WHERE run_id IN ('T2_movie','T2_final')"
    ).fetchone()[0]
    stride = max(1, (n + 3999) // 4000)
    return cur.execute(
        f"SELECT t_s/{MYR}, {expr} FROM sample "
        "WHERE run_id IN ('T2_movie','T2_final') AND idx % ? = 0 "
        "ORDER BY t_s", (stride,)).fetchall()


def svg_chart(points, title, y_label, w=920, h=240):
    pad_l, pad_r, pad_t, pad_b = 64, 16, 28, 34
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if y1 - y0 < 1e-12:
        y0, y1 = y0 - 1.0, y1 + 1.0
    def px(x):
        return pad_l + (x - x0) / (x1 - x0) * (w - pad_l - pad_r)
    def py(y):
        return h - pad_b - (y - y0) / (y1 - y0) * (h - pad_t - pad_b)
    d = " ".join(
        ("M" if i == 0 else "L") + f"{px(x):.2f},{py(y):.2f}"
        for i, (x, y) in enumerate(points))
    return f"""<figure>
<figcaption>{title}</figcaption>
<svg viewBox="0 0 {w} {h}" role="img" aria-label="{title}">
<rect x="{pad_l}" y="{pad_t}" width="{w-pad_l-pad_r}" height="{h-pad_t-pad_b}"
 class="plot-bg"/>
<path d="{d}" class="line"/>
<text x="{pad_l-6}" y="{py(y1)+4:.2f}" class="ax" text-anchor="end">{y1:.6g}</text>
<text x="{pad_l-6}" y="{py(y0)+4:.2f}" class="ax" text-anchor="end">{y0:.6g}</text>
<text x="{pad_l}" y="{h-10}" class="ax">{x0:.3g} Myr</text>
<text x="{w-pad_r}" y="{h-10}" class="ax" text-anchor="end">{x1:.3g} Myr</text>
<text x="14" y="{(pad_t+h-pad_b)/2:.0f}" class="ax"
 transform="rotate(-90 14 {(pad_t+h-pad_b)/2:.0f})" text-anchor="middle">{y_label}</text>
</svg>
</figure>"""


def main():
    if not DB.exists():
        print(f"database not found: {DB} - run the test-2 notebook first",
              file=sys.stderr)
        return 1
    con = sqlite3.connect(DB)
    cur = con.cursor()

    xt = dict(cur.execute(
        "SELECT key, value FROM run_extra WHERE run_id='T2_final'").fetchall())
    a0 = cur.execute("SELECT a0_m FROM run WHERE run_id='T2_final'").fetchone()[0]
    n0 = math.sqrt(6.67430e-11 * (1.98847e30 + 3.3011e23) / a0**3)
    pw_pred = xt["gr_rate_rad_s"] + xt["ll_A11_rad_s"]
    expected = 1.5 + pw_pred / n0
    t_cap = cur.execute(
        "SELECT t_s FROM event WHERE run_id='T2_final' "
        "AND event='capture_detected' ORDER BY t_s LIMIT 1").fetchone()[0]
    rows = cur.execute(
        "SELECT ratio FROM sample WHERE run_id='T2_final' AND t_s > ? "
        "ORDER BY idx", (t_cap + MYR,)).fetchall()
    mean_ratio = sum(r[0] for r in rows) / len(rows)
    fin = cur.execute(
        "SELECT P_orb_s/86400.0, P_rot_s/86400.0, pomega_rad FROM sample "
        "WHERE run_id='T2_final' ORDER BY idx DESC LIMIT 1").fetchone()
    nb = cur.execute(
        "SELECT SUM(captured), COUNT(*) FROM branch WHERE run_id='T2_sweep'"
    ).fetchone()

    ratio_pts = [(t, math.log10(v)) for t, v in fetch_series(cur, "ratio")]
    e_pts = fetch_series(cur, "e")
    pw_pts = fetch_series(cur, "pomega_rad")

    cards = [
        ("Einstein's rate", f"{xt['gr_rate_rad_s']*ARCSEC_CY:.2f}&Prime;/cy",
         "GR perihelion advance"),
        ("Jupiter's rate", f"{xt['ll_A11_rad_s']*ARCSEC_CY:.2f}&Prime;/cy",
         "Laplace-Lagrange A11"),
        ("Settled mean spin ratio", f"{mean_ratio:.9f}",
         f"predicted 1.5 + &#x03D6;&#x0307;/n = {expected:.9f}"),
        ("Capture time", f"{t_cap/MYR:.3f} Myr",
         "movie clock (tides 1000&times;)"),
        ("Today's year / rotation", f"{fin[0]:.4f} / {fin[1]:.4f} d",
         "observed: 87.969 / 58.646 d"),
        ("Sweep captures", f"{nb[0]}/{nb[1]}",
         "phase branches locked with full physics"),
        ("Ellipse turns in the run", f"{fin[2]/(2*math.pi):.2f}",
         "the perihelion odometer's total"),
    ]
    cards_html = "\n".join(
        f'<div class="card"><div class="v">{v}</div>'
        f'<div class="k">{k}</div><div class="s">{s}</div></div>'
        for k, v, s in cards)

    charts = "\n".join([
        svg_chart(ratio_pts,
                  "The braking descent into the 3:2 lock (log scale: "
                  "181&times; the orbital rate down to 1.5)",
                  "log10(&Omega;/n)"),
        svg_chart(e_pts,
                  "Jupiter breathes Mercury's eccentricity — through braking, "
                  "capture, and lock",
                  "eccentricity e"),
        svg_chart(pw_pts,
                  "The perihelion odometer: Einstein + Jupiter turn the "
                  "ellipse (radians)",
                  "&#x03D6; [rad]"),
    ])

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mercury test 2 — Jupiter and Einstein</title>
<style>
body {{ font: 15px/1.5 -apple-system, "Segoe UI", Roboto, sans-serif;
  margin: 0; background: #0b0f1a; color: #dce3f0; }}
main {{ max-width: 980px; margin: 0 auto; padding: 24px 16px 48px; }}
h1 {{ font-size: 1.5rem; margin: 8px 0 2px; }}
.sub {{ color: #8fa0bd; margin: 0 0 18px; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 10px; margin: 0 0 22px; }}
.card {{ background: #141b2d; border: 1px solid #24304a; border-radius: 10px;
  padding: 12px 14px; }}
.card .v {{ font-size: 1.15rem; font-weight: 600; color: #ffd479; }}
.card .k {{ font-size: .85rem; color: #aebadd; margin-top: 2px; }}
.card .s {{ font-size: .75rem; color: #7684a3; }}
figure {{ margin: 0 0 22px; background: #141b2d; border: 1px solid #24304a;
  border-radius: 10px; padding: 12px; }}
figcaption {{ font-size: .9rem; color: #aebadd; margin-bottom: 6px; }}
svg {{ width: 100%; height: auto; display: block; }}
.plot-bg {{ fill: #0e1424; stroke: #24304a; }}
.line {{ fill: none; stroke: #ffd479; stroke-width: 1.4; }}
.ax {{ fill: #7684a3; font-size: 11px; }}
footer {{ color: #7684a3; font-size: .8rem; border-top: 1px solid #24304a;
  padding-top: 12px; }}
</style>
</head>
<body>
<main>
<h1>Mercury, Jupiter, and Einstein — the lock follows the precessing ellipse</h1>
<p class="sub">Test 2 of the planet_Mercury project: the 3:2 spin-orbit capture
with general relativity and Jupiter's secular forcing, integrated by the
pure-Rust SUNDIALS 7.8.0 CVODE solver at relative tolerance 1e-12. Every
number below was queried from <code>data/mercury_test2.sqlite3</code>.</p>
<div class="cards">
{cards_html}
</div>
{charts}
<footer>Baked deterministically by gui/bake_page2.py from the test-2 database.
The settled spin ratio sits above exactly 3/2 by the perihelion turn rate over
the mean motion — the relativistic-plus-planetary fingerprint on Mercury's
clock. Tidal strength in the movie runs is compressed 1000&times; (documented);
Einstein's and Jupiter's rates are real, uncompressed rates.</footer>
</main>
</body>
</html>
"""
    OUT.write_text(html, encoding="utf-8")
    print(f"baked {OUT} ({OUT.stat().st_size} bytes, "
          f"{len(ratio_pts)}+{len(e_pts)}+{len(pw_pts)} chart points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
