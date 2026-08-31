#!/usr/bin/env python3
"""Bake the self-contained test-2 display page (Jupiter + Einstein).

Reads the SQLite database the test-2 notebook builds
(data/mercury_test2.sqlite3), decimates the canonical history to
plot-resolution arrays, and writes ONE self-contained HTML file
(gui/mercury_test2.html) in the rustSolveIt engine's recorded-player style:
all data embedded as JavaScript constants, vanilla JavaScript + 2-D canvas,
zero network requests, dark fixed palette.

The page only REPLAYS recorded samples (no integration in the page — the
orbit panel's moving dot is display math on the recorded orbit shape, the
recorded perihelion angle, and the recorded spin/orbit ratio). What is new
versus test 1's player: the ellipse itself visibly PRECESSES (the recorded
pomega turns it), the eccentricity breathes on Jupiter's cycle, and the HUD
tracks the perihelion odometer. Deterministic: same database ->
byte-identical page. Standard library only.

Usage: python3 gui/bake_page2.py   (from planet_Mercury/, or pass paths)
"""

import math
import sqlite3
import sys
from pathlib import Path

YEAR = 3.15576e7
ARCSEC_CY = 206264.80624709636 * 3.15576e9


def fnum(x, fmt):
    s = fmt % x
    # Normalize negative zero for determinism across query paths.
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    return s


def js_array(values, fmt):
    return "[" + ",".join(fnum(v, fmt) for v in values) + "]"


def decimate(rows, cap):
    if len(rows) <= cap:
        return rows
    step = (len(rows) + cap - 1) // cap
    out = rows[::step]
    if out[-1] != rows[-1]:
        out.append(rows[-1])
    return out


def bake(db_path: Path, out_path: Path) -> None:
    con = sqlite3.connect(str(db_path))
    cur = con.cursor()

    hist = cur.execute(
        "SELECT t_s, a_m, e, ratio, gamma2_rad, P_orb_s, P_rot_s, pomega_rad"
        " FROM sample WHERE run_id IN ('T2_movie','T2_final') ORDER BY t_s"
    ).fetchall()
    if not hist:
        raise SystemExit(
            "no T2_movie/T2_final samples in the database — run the test-2 notebook first")

    cap_row = cur.execute(
        "SELECT t_s FROM event WHERE run_id='T2_final'"
        " AND event='capture_detected' ORDER BY t_s LIMIT 1"
    ).fetchone()
    t_capture = cap_row[0] if cap_row else None

    dense = []
    if t_capture is not None:
        dense = cur.execute(
            "SELECT t_s, gamma2_rad, ratio FROM sample WHERE run_id='T2_final'"
            " AND t_s >= ? AND t_s <= ? ORDER BY t_s",
            (t_capture - 120_000 * YEAR, t_capture + 50_000 * YEAR),
        ).fetchall()

    sweep = cur.execute(
        "SELECT COUNT(*), SUM(captured) FROM branch WHERE run_id='T2_sweep'"
    ).fetchone()
    n_branches = sweep[0] or 0
    n_captured = sweep[1] or 0

    events = cur.execute(
        "SELECT run_id, t_s, event, value FROM event"
        " WHERE run_id IN ('T2_movie','T2_final') ORDER BY t_s"
    ).fetchall()

    run_bf = cur.execute(
        "SELECT tau_lag_s, compression, a0_m FROM run WHERE run_id='T2_final'"
    ).fetchone()
    xt = dict(cur.execute(
        "SELECT key, value FROM run_extra WHERE run_id='T2_final'").fetchall())

    settled = []
    if t_capture is not None:
        settled = cur.execute(
            "SELECT ratio FROM sample WHERE run_id='T2_final' AND t_s > ?"
            " ORDER BY idx", (t_capture + 1.0e6 * YEAR,)).fetchall()
    con.close()

    n0 = math.sqrt(6.67430e-11 * (1.98847e30 + 3.3011e23) / run_bf[2] ** 3)
    pw_dot = xt["gr_rate_rad_s"] + xt["ll_A11_rad_s"]
    expected_ratio = 1.5 + pw_dot / n0
    mean_ratio = (sum(r[0] for r in settled) / len(settled)) if settled else 0.0
    e_min = min(r[2] for r in hist)
    e_max = max(r[2] for r in hist)
    turns = hist[-1][7] / (2.0 * math.pi)

    hist = decimate(hist, 6000)
    dense = decimate(dense, 20000)

    last = hist[-1]
    p_orb_d = last[5] / 86400.0
    p_rot_d = last[6] / 86400.0
    p_solar_d = 1.0 / abs(1.0 / p_rot_d - 1.0 / p_orb_d)

    data_js = []
    data_js.append(
        "const H={"
        + "t:" + js_array([r[0] / YEAR for r in hist], "%.1f")
        + ",a:" + js_array([r[1] for r in hist], "%.8e")
        + ",e:" + js_array([r[2] for r in hist], "%.8f")
        + ",ratio:" + js_array([r[3] for r in hist], "%.8f")
        + ",gamma2:" + js_array([r[4] for r in hist], "%.6f")
        + ",porb:" + js_array([r[5] / 86400.0 for r in hist], "%.6f")
        + ",prot:" + js_array([r[6] / 86400.0 for r in hist], "%.6f")
        + ",pomega:" + js_array([r[7] for r in hist], "%.7f")
        + "};"
    )
    data_js.append(
        "const DENSE={"
        + "t:" + js_array([r[0] / YEAR for r in dense], "%.2f")
        + ",gamma2:" + js_array([r[1] for r in dense], "%.6f")
        + ",ratio:" + js_array([r[2] for r in dense], "%.8f")
        + "};"
    )
    ev_js = ",".join(
        '{run:"%s",t:%s,name:"%s",v:%s}'
        % (r[0], fnum(r[1] / YEAR, "%.1f"), r[2], fnum(r[3], "%.6f"))
        for r in events
    )
    data_js.append("const EVENTS=[" + ev_js + "];")
    meta = {
        "t_capture_yr": fnum(t_capture / YEAR, "%.1f") if t_capture else "null",
        "branches": str(n_branches),
        "captured": str(n_captured),
        "p_orb_d": fnum(p_orb_d, "%.5f"),
        "p_rot_d": fnum(p_rot_d, "%.5f"),
        "p_solar_d": fnum(p_solar_d, "%.5f"),
        "final_ratio": fnum(last[3], "%.8f"),
        "compression": fnum(run_bf[1] if run_bf else 1000.0, "%.0f"),
        "gr_as": fnum(xt["gr_rate_rad_s"] * ARCSEC_CY, "%.2f"),
        "jup_as": fnum(xt["ll_A11_rad_s"] * ARCSEC_CY, "%.2f"),
        "mean_ratio": fnum(mean_ratio, "%.9f"),
        "expected_ratio": fnum(expected_ratio, "%.9f"),
        "e_min": fnum(e_min, "%.5f"),
        "e_max": fnum(e_max, "%.5f"),
        "turns": fnum(turns, "%.2f"),
    }
    data_js.append(
        "const META={t_capture_yr:%(t_capture_yr)s,branches:%(branches)s,"
        "captured:%(captured)s,p_orb_d:%(p_orb_d)s,p_rot_d:%(p_rot_d)s,"
        "p_solar_d:%(p_solar_d)s,final_ratio:%(final_ratio)s,"
        "compression:%(compression)s,gr_as:%(gr_as)s,jup_as:%(jup_as)s,"
        "mean_ratio:%(mean_ratio)s,expected_ratio:%(expected_ratio)s,"
        "e_min:%(e_min)s,e_max:%(e_max)s,turns:%(turns)s};" % meta
    )

    caption = (
        "Test 2: Mercury's 3:2 capture with Einstein's relativity (%s arcsec/century) and "
        "Jupiter's secular forcing (%s arcsec/century; the eccentricity breathes %s-%s on a "
        "~809,000-year cycle). Capture at %s million simulated years (tides compressed x%s; "
        "the GR and Jupiter rates are real). The ellipse turns %s times in the window, and the "
        "settled mean spin ratio is %s vs predicted 1.5 + pomega_dot/n = %s: the lock follows "
        "the PRECESSING ellipse, not the stars. Sweep: %s of %s phase branches captured. "
        "Playback replays recorded SUNDIALS CVODE samples at display speed (not real time)."
    ) % (
        meta["gr_as"], meta["jup_as"], meta["e_min"], meta["e_max"],
        fnum(t_capture / YEAR / 1e6, "%.2f") if t_capture else "?",
        meta["compression"], meta["turns"], meta["mean_ratio"],
        meta["expected_ratio"], meta["captured"], meta["branches"],
    )

    html = PAGE.replace("__DATA__", "\n".join(data_js)).replace("__CAPTION__", caption)
    out_path.write_text(html, encoding="utf-8")
    print(
        "baked %s: %d history + %d dense points, %.2f MB"
        % (out_path, len(hist), len(dense), len(html) / 1e6)
    )


PAGE = r"""<!doctype html>
<meta charset="utf-8">
<title>Mercury test 2 — Jupiter and Einstein</title>
<style>
:root{--bg:#0d1117;--panel:#161b22;--line:#30363d;--fg:#e6edf3;--dim:#8b949e;
--accent:#5d84a8;--gold:#d9a441;--green:#7fd18b;--red:#e06c75;
font:14px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--bg);color:var(--fg);display:flex;flex-direction:column;height:100vh}
header{padding:8px 14px;border-bottom:1px solid var(--line)}
header h1{font-size:16px}
header p{color:var(--dim);font-size:11.5px;max-width:1200px}
#main{flex:1;display:flex;min-height:0}
#left{flex:1.1;position:relative;border-right:1px solid var(--line)}
#orbit{width:100%;height:100%;display:block}
#hud{position:absolute;top:8px;left:8px;background:rgba(13,17,23,.85);
border:1px solid var(--line);border-radius:6px;padding:8px 10px;font-size:.76rem;
pointer-events:none;line-height:1.65}
#hud b{color:var(--fg)} #hud .k{color:var(--dim)} #hud .good{color:var(--green)}
#right{flex:1;display:flex;flex-direction:column;min-width:0}
#right canvas{flex:1;width:100%;min-height:0;border-bottom:1px solid var(--line)}
footer{display:flex;gap:8px;align-items:center;padding:8px 12px;border-top:1px solid var(--line);background:var(--panel)}
button,select{background:#23323f;color:var(--fg);border:1px solid #3c4f61;border-radius:4px;
padding:4px 10px;font:inherit;cursor:pointer}
button.primary{background:#1f5131}
#scrub{flex:1}
</style>
<header>
 <h1>Mercury test 2 &mdash; Jupiter and Einstein: the lock follows the precessing ellipse</h1>
 <p>__CAPTION__</p>
</header>
<div id="main">
 <div id="left"><canvas id="orbit"></canvas>
  <div id="hud"></div>
 </div>
 <div id="right">
  <canvas id="p1"></canvas><canvas id="p2"></canvas><canvas id="p3"></canvas><canvas id="p4"></canvas>
 </div>
</div>
<footer>
 <button id="bt-play" class="primary">&#9654; Play</button>
 <button id="bt-back">&#9664;&#9664;</button>
 <button id="bt-fwd">&#9654;&#9654;</button>
 <button id="bt-cap">Jump to capture</button>
 <input id="scrub" type="range" min="0" max="0" value="0">
 <select id="speed"><option>1</option><option>2</option><option selected>5</option><option>15</option><option>40</option></select>
 <span style="color:var(--dim);font-size:11px">Space = play/pause, arrows = step</span>
</footer>
<script>
"use strict";
__DATA__
const N=H.t.length;
let idx=0,playing=false,speed=5,animM=0,animTheta=0,trail=[];
const $=id=>document.getElementById(id);
const cO=$("orbit"),cP=[$("p1"),$("p2"),$("p3"),$("p4")];
function fit(c){const d=window.devicePixelRatio||1,r=c.getBoundingClientRect();
 c.width=Math.max(50,r.width*d);c.height=Math.max(50,r.height*d);
 const x=c.getContext("2d");x.setTransform(d,0,0,d,0,0);return x;}
// Display-only Kepler solve (recorded state -> screen position; no physics
// integration happens in this page).
function keplerPos(a,e,M){let E=M;for(let i=0;i<12;i++){E-=(E-e*Math.sin(E)-M)/(1-e*Math.cos(E));}
 const f=2*Math.atan2(Math.sqrt(1+e)*Math.sin(E/2),Math.sqrt(1-e)*Math.cos(E/2));
 const r=a*(1-e*e)/(1+e*Math.cos(f));return {x:r*Math.cos(f),y:r*Math.sin(f),f:f};}
function logt(t){return Math.log10(Math.max(t,1000));}
const TX0=logt(H.t[0]),TX1=logt(H.t[N-1]);
const TMYR1=H.t[N-1]/1e6;
function drawPlot(ctx,c,title,opts){
 const W=c.getBoundingClientRect().width,Hh=c.getBoundingClientRect().height;
 ctx.clearRect(0,0,W,Hh);
 const mL=56,mR=10,mT=20,mB=22,pw=W-mL-mR,ph=Hh-mT-mB;
 ctx.strokeStyle="#30363d";ctx.strokeRect(mL,mT,pw,ph);
 ctx.fillStyle="#8b949e";ctx.font="11px ui-monospace";ctx.fillText(title,mL,13);
 const x0=opts.x0,x1=opts.x1,y0=opts.y0,y1=opts.y1;
 const X=v=>mL+pw*(v-x0)/(x1-x0), Y=v=>mT+ph*(1-(v-y0)/(y1-y0));
 (opts.rules||[]).forEach(rv=>{const yy=Y(opts.ymap(rv));if(yy>mT&&yy<mT+ph){
   ctx.strokeStyle="#d9a44166";ctx.beginPath();ctx.moveTo(mL,yy);ctx.lineTo(mL+pw,yy);ctx.stroke();
   ctx.fillStyle="#d9a441";ctx.fillText(String(rv),mL+pw-46,yy-3);}});
 for(const s of opts.series){ctx.strokeStyle=s.color;ctx.beginPath();let started=false;
  for(let i=0;i<s.xs.length;i++){const xx=X(s.xmap(s.xs[i])),yy=Y(s.ymap(s.ys[i]));
   if(!isFinite(xx)||!isFinite(yy))continue;
   if(!started){ctx.moveTo(xx,yy);started=true;}else ctx.lineTo(xx,yy);}
  ctx.stroke();}
 ctx.fillStyle="#8b949e";
 ctx.fillText(opts.xlabel,mL+pw/2-30,Hh-6);
 (opts.yticks||[]).forEach(tv=>{const yy=Y(opts.ymap(tv));if(yy>mT-1&&yy<mT+ph+1){
   ctx.fillText(String(tv),4,yy+4);}});
 if(opts.cursor!==undefined){const xx=X(opts.cursor);if(xx>=mL&&xx<=mL+pw){
   ctx.strokeStyle="#e06c75";ctx.beginPath();ctx.moveTo(xx,mT);ctx.lineTo(xx,mT+ph);ctx.stroke();}}
 opts.series.forEach((s,i)=>{ctx.fillStyle=s.color;ctx.fillText(s.name,mL+8+i*150,mT+12);});
}
function draw(){
 const t=H.t[idx],a=H.a[idx],e=H.e[idx],ratio=H.ratio[idx],pw=H.pomega[idx];
 // ---- orbit panel: the Sun pinned at the focus, the ellipse turned by the
 // recorded perihelion angle pomega (Einstein + Jupiter turn it for real) ----
 const ctx=fit(cO),W=cO.getBoundingClientRect().width,Hh=cO.getBoundingClientRect().height;
 ctx.clearRect(0,0,W,Hh);
 const SC=Math.min(W,Hh)/(2.6*a);
 const sx=W*0.5,sy=Hh*0.52; // the Sun (the focus) pinned at screen center
 const b=a*Math.sqrt(1-e*e);
 const cx=sx-a*e*Math.cos(pw)*SC, cy=sy+a*e*Math.sin(pw)*SC; // ellipse center
 ctx.strokeStyle="#3c4f61";ctx.beginPath();
 ctx.ellipse(cx,cy,a*SC,b*SC,-pw,0,2*Math.PI);ctx.stroke();
 // apsidal line: aphelion through the Sun to perihelion (it precesses!)
 const px=sx+a*(1-e)*Math.cos(pw)*SC, py=sy-a*(1-e)*Math.sin(pw)*SC;
 const qx=sx-a*(1+e)*Math.cos(pw)*SC, qy=sy+a*(1+e)*Math.sin(pw)*SC;
 ctx.strokeStyle="#d9a44155";ctx.beginPath();ctx.moveTo(qx,qy);ctx.lineTo(px,py);ctx.stroke();
 ctx.fillStyle="#d9a441";ctx.beginPath();ctx.arc(px,py,3,0,2*Math.PI);ctx.fill();
 ctx.font="10px ui-monospace";ctx.fillText("perihelion",px+6,py-6);
 ctx.fillStyle="#d9a441";ctx.beginPath();ctx.arc(sx,sy,9,0,2*Math.PI);ctx.fill();
 const p=keplerPos(a,e,animM);
 const ang=p.f+pw, r=Math.hypot(p.x,p.y);
 const mx=sx+r*Math.cos(ang)*SC, my=sy-r*Math.sin(ang)*SC;
 trail.push([mx,my]);if(trail.length>240)trail.shift();
 ctx.strokeStyle="#5d84a866";ctx.beginPath();
 trail.forEach((q,i)=>{i?ctx.lineTo(q[0],q[1]):ctx.moveTo(q[0],q[1]);});ctx.stroke();
 ctx.fillStyle="#5d84a8";ctx.beginPath();ctx.arc(mx,my,7,0,2*Math.PI);ctx.fill();
 // long-axis arrow (the "handle") at the display spin angle
 ctx.strokeStyle="#e6edf3";ctx.lineWidth=2;ctx.beginPath();
 ctx.moveTo(mx-14*Math.cos(animTheta),my+14*Math.sin(animTheta));
 ctx.lineTo(mx+14*Math.cos(animTheta),my-14*Math.sin(animTheta));ctx.stroke();ctx.lineWidth=1;
 ctx.fillStyle="#e06c75";ctx.beginPath();
 ctx.arc(mx+14*Math.cos(animTheta),my-14*Math.sin(animTheta),3,0,2*Math.PI);ctx.fill();
 ctx.fillStyle="#8b949e";ctx.font="11px ui-monospace";
 ctx.fillText("spin/orbit ratio dial",12,Hh-34);
 ctx.font="26px ui-monospace";
 ctx.fillStyle=Math.abs(ratio-1.5)<1e-3?"#7fd18b":"#e6edf3";
 ctx.fillText(ratio.toFixed(7),12,Hh-10);
 // ---- HUD ----
 const locked=Math.abs(ratio-1.5)<5e-3&&META.t_capture_yr&&t>=META.t_capture_yr;
 $("hud").innerHTML=
  '<span class="k">simulated year</span> <b>'+t.toLocaleString("en-US",{maximumFractionDigits:0})+'</b> (tides x'+META.compression+'; GR + Jupiter real)<br>'+
  '<span class="k">spin/orbit ratio</span> <b>'+ratio.toFixed(8)+'</b> <span class="k">locked mean '+META.mean_ratio+'</span><br>'+
  '<span class="k">eccentricity e</span> <b>'+e.toFixed(5)+'</b> <span class="k">breathes '+META.e_min+'-'+META.e_max+'</span><br>'+
  '<span class="k">ellipse turns</span> <b>'+(pw/(2*Math.PI)).toFixed(3)+'</b> <span class="k">of '+META.turns+' total ('+META.gr_as+'&Prime; GR + '+META.jup_as+'&Prime; Jupiter /cy)</span><br>'+
  '<span class="k">year P_orb</span> <b>'+H.porb[idx].toFixed(3)+' d</b> <span class="k">obs 87.969</span><br>'+
  '<span class="k">status</span> '+(locked?'<b class="good">LOCKED 3:2 &mdash; riding the precessing ellipse</b>':'<b>braking...</b>');
 // ---- plots ----
 const tmap=logt, tmyr=v=>v/1e6;
 drawPlot(fit(cP[0]),cP[0],"spin/orbit ratio Omega/n vs time (log-log)",{
  x0:TX0,x1:TX1,y0:Math.log10(1.2),y1:Math.log10(200),
  xlabel:"log10 t [yr]",cursor:tmap(t),
  ymap:v=>Math.log10(v),yticks:[1.5,3,10,30,100],rules:[1.5,1.256],
  series:[{name:"Omega/n",color:"#5d84a8",xs:H.t,ys:H.ratio,xmap:tmap,ymap:v=>Math.log10(v)}]});
 const dcap=META.t_capture_yr||0;
 drawPlot(fit(cP[1]),cP[1],"resonance angle gamma2 = 2 theta - 3 M - 2 pomega: circulation -> libration",{
  x0:(DENSE.t[0]||0)-dcap,x1:(DENSE.t[DENSE.t.length-1]||1)-dcap,y0:-3.4,y1:3.4,
  xlabel:"t - t_capture [yr]",cursor:t-dcap,
  ymap:v=>v,yticks:[-3,0,3],rules:[],
  series:[{name:"gamma2 [rad]",color:"#7fd18b",xs:DENSE.t.map(v=>v-dcap),ys:DENSE.gamma2,xmap:v=>v,ymap:v=>v}]});
 drawPlot(fit(cP[2]),cP[2],"Jupiter breathes the eccentricity (through braking, capture, and lock)",{
  x0:0,x1:TMYR1,y0:0.195,y1:0.208,
  xlabel:"t [Myr]",cursor:tmyr(t),
  ymap:v=>v,yticks:[0.196,0.2,0.204],rules:[0.20563],
  series:[{name:"e (LL cycle ~809 kyr)",color:"#d9a441",xs:H.t,ys:H.e,xmap:tmyr,ymap:v=>v}]});
 drawPlot(fit(cP[3]),cP[3],"the perihelion odometer: Einstein + Jupiter turn the ellipse",{
  x0:0,x1:TMYR1,y0:0,y1:105,
  xlabel:"t [Myr]",cursor:tmyr(t),
  ymap:v=>v,yticks:[0,31.4,62.8,94.2],rules:[],
  series:[{name:"pomega [rad] (2 pi per turn)",color:"#7fd18b",xs:H.t,ys:H.pomega,xmap:tmyr,ymap:v=>v}]});
 $("scrub").value=idx;
}
function tick(){
 if(playing){idx=Math.min(N-1,idx+speed);if(idx>=N-1){playing=false;$("bt-play").innerHTML="&#9654; Play";}}
 const dM=0.11;animM+=dM;animTheta+=H.ratio[idx]*dM;
 draw();requestAnimationFrame(tick);}
$("bt-play").onclick=()=>{if(idx>=N-1)idx=0;playing=!playing;
 $("bt-play").innerHTML=playing?"&#10074;&#10074; Pause":"&#9654; Play";};
$("bt-back").onclick=()=>{idx=Math.max(0,idx-1);};
$("bt-fwd").onclick=()=>{idx=Math.min(N-1,idx+1);};
$("bt-cap").onclick=()=>{if(META.t_capture_yr){for(let i=0;i<N;i++)if(H.t[i]>=META.t_capture_yr){idx=i;break;}}};
$("scrub").max=N-1;
$("scrub").oninput=e=>{idx=+e.target.value;};
$("speed").onchange=e=>{speed=+e.target.value;};
window.addEventListener("keydown",e=>{
 if(e.code==="Space"){e.preventDefault();$("bt-play").onclick();}
 else if(e.key==="ArrowLeft")$("bt-back").onclick();
 else if(e.key==="ArrowRight")$("bt-fwd").onclick();});
requestAnimationFrame(tick);
</script>
"""


def main() -> int:
    base = Path(__file__).resolve().parent.parent
    db = Path(sys.argv[1]) if len(sys.argv) > 1 else base / "data" / "mercury_test2.sqlite3"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else base / "gui" / "mercury_test2.html"
    bake(db, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
