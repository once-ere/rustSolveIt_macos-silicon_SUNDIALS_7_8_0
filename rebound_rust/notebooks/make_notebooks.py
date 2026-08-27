"""Generates one Jupyter notebook per example program, for BOTH crates.

Golden rule of this project: every example has a concomitant notebook.
Each notebook builds its example with cargo, runs it, and shows the output
(plotting or decoding the raw-bit dumps where that is useful).

Run:  python make_notebooks.py
Part of the rebound_rs / reboundx_rs port. GPL-3.0-or-later.
"""
import json
import os

# Derived, never hard-coded: this file lives at
# <workspace>/rebound_rust/notebooks/make_notebooks.py, so two levels
# up is the workspace holding both crates. Hard-coding an absolute
# path here would bake one machine's layout into every notebook.
_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(_HERE))
REBOUND = os.path.join(ROOT, "rebound_rust")
REBOUNDX = os.path.join(ROOT, "reboundx_rust")
OUT = os.path.join(REBOUND, "notebooks")

# Helper source injected into notebooks that decode raw-bit dumps.
DECODE = '''import struct

def unbits(h):
    """Turn a 16-hex-digit IEEE-754 bit pattern back into a float."""
    return struct.unpack("<d", int(h, 16).to_bytes(8, "little"))[0]

def read_state(path):
    """Read one of the raw-bit state dumps into {label: [floats]}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            key, rest = parts[0], parts[1:]
            vals = []
            for tok in rest:
                if len(tok) == 16:
                    try:
                        vals.append(unbits(tok))
                        continue
                    except ValueError:
                        pass
                vals.append(tok)
            out.setdefault(key, []).append(vals)
    return out
'''

COMPARE = '''def compare(a, b, label_a="C", label_b="Rust"):
    """Byte-compare two dump files and report."""
    ta = open(a, "rb").read().replace(b"\\r\\n", b"\\n")
    tb = open(b, "rb").read().replace(b"\\r\\n", b"\\n")
    if ta == tb:
        print(f"BIT-IDENTICAL: {label_a} and {label_b} agree on every bit")
        return True
    print(f"MISMATCH between {label_a} and {label_b}")
    la, lb = ta.decode().splitlines(), tb.decode().splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            print(f"  line {i}:\\n    {label_a}: {x}\\n    {label_b}: {y}")
    return False
'''

# name -> (crate, title, description, argv, extra cells)
EX = {}


def add(crate, name, title, desc, args=(), post=None, plot=None, runs=None):
    EX[name] = dict(crate=crate, title=title, desc=desc, args=list(args),
                    post=post, plot=plot, runs=runs)


# ----------------------------- REBOUND --------------------------------
add("rebound", "shearing_sheet",
    "shearing_sheet — a patch of Saturn's rings",
    "A straight port of REBOUND's flagship example: a small sheared box of "
    "colliding ice particles, using the SEI integrator, octree self-gravity, "
    "tree collision search and shear-periodic boundaries. "
    "**Note:** like the stock C example, `shearing_sheet` integrates to "
    "infinity — you stop it with Ctrl+C. A notebook cannot run that, so this "
    "notebook builds it (to prove it compiles) and then runs the terminating "
    "seeded harness `shearing_sheet_test` for 400 steps to produce the plot. "
    "The two set up identical physics; the harness just fixes the random seed "
    "and stops.",
    args=["400"], plot="shearing", runs="shearing_sheet_test")

add("rebound", "shearing_sheet_test",
    "shearing_sheet_test — the byte-identity acceptance test",
    "Runs the seeded 400-step shearing sheet and, if the C reference dump is "
    "present, verifies the two agree on every bit and prints both SHA-256 "
    "fingerprints. This is the project's headline acceptance test: 1,482 "
    "particles and ~102,500 collisions.",
    args=["400"], post="sha")

add("rebound", "integrators_test",
    "integrators_test — the 63-configuration matrix",
    "Runs the fixed three-body problem under a chosen integrator "
    "configuration and dumps the final state as raw IEEE-754 bits. "
    "63 configurations of this harness were verified bit-identical against "
    "the MSVC C build. Change the arguments below to try others "
    "(e.g. 'ias15', 'saba-h1064', 'mercurius', 'trace').",
    args=["whfast", "2", "500"], post="state")

add("rebound", "libm_diff",
    "libm_diff — which maths functions agree with C?",
    "Samples 200,000 inputs per maths function and dumps raw bits, so the C "
    "and Rust results can be compared exactly. This is the measurement that "
    "underpins the whole project: sin, cos, tan, atan2, sqrt, fmod, exp, log "
    "and cbrt are bit-identical to Microsoft's C library, and pow is the one "
    "exception.",
    post="libm")

add("rebound", "bs_pow_diff",
    "bs_pow_diff — isolating the one library difference",
    "The BS integrator picks its own step size using pow(). This evaluates "
    "exactly those pow() calls (200,000 samples) so the disagreement can be "
    "measured precisely: 56 differ, every one by exactly 1 ULP — the smallest "
    "difference two numbers can have.",
    post="powdiff")

add("rebound", "kepler_rectilinear",
    "kepler_rectilinear — a real bug this test suite caught",
    "Calls the Kepler solver with (near-)rectilinear hyperbolic motion — the "
    "regime that strains its iteration hardest. This probe found a genuine "
    "port defect: with a large timestep the Rust looped forever where the C "
    "returned, because `while (a > b)` and `if (a <= b) break` differ when a "
    "value is NaN. Now fixed; the two agree bit-for-bit.",
    args=["1e-12", "1", "2.0"], post="plain")

add("rebound", "movetocom_var",
    "movetocom_var — the variational centre-of-mass fix",
    "Probes reb_simulation_move_to_com with MEGNO variational particles. An "
    "audit found the first-order shift summed the wrong particle array, which "
    "silently changed every MEGNO/Lyapunov result. Now fixed and verified "
    "bit-identical against C.",
    post="movetocom")

add("rebound", "movetocom_var_test",
    "movetocom_var_test — the audit's own probe",
    "The probe written during the audit to demonstrate the variational "
    "centre-of-mass defect, kept as a regression check. It prints both dm "
    "accumulators so you can see why the wrong array mattered.",
    post="plain")

add("rebound", "derivatives_test",
    "derivatives_test — 65 orbital derivative functions",
    "Evaluates all 65 reb_particle_derivative_* functions on two "
    "configurations and dumps raw bits; verified 130/130 lines bit-identical "
    "against the C build.",
    post="head")

add("rebound", "frequency_test",
    "frequency_test — MFT / FMFT / FMFT2 frequency analysis",
    "Runs REBOUND's frequency analysis in all three modes on a synthetic "
    "three-frequency signal (true frequencies 0.30, 0.55, 0.11 radians per "
    "sample) and prints the recovered frequencies, amplitudes and phases.",
    post="freq")

add("rebound", "archive_test",
    "archive_test — Simulationarchive round trip",
    "Writes a 3-snapshot Simulationarchive. The same file loads in the C "
    "build and continues bit-identically, and vice versa — the formats are "
    "fully interchangeable.",
    args=["whfast-usafe", "write"], post="archive")

add("rebound", "server_test",
    "server_test — the REBOUND web server in pure Rust",
    "Starts the ported web server, fetches the /simulation binary over HTTP, "
    "then shuts it down with the 'Q' key endpoint. The served blob is a valid "
    "REBOUND binary that the C build loads to the identical state.",
    post="server")

add("rebound", "addfmt_test",
    "addfmt_test — creating particles from orbital elements",
    "Adds the built-in solar system plus particles specified by orbital "
    "elements, by Pal coordinates, and by orbital period. Verified "
    "bit-identical against C.",
    post="head")

# ----------------------------- REBOUNDx -------------------------------
add("reboundx", "tides_spin_pseudo",
    "tides_spin_pseudo — pseudo-synchronisation of a hot Jupiter",
    "A giant planet orbiting very close to its star, starting slightly "
    "elliptical, tilted 30 degrees and spinning fast. Tides should circularise "
    "the orbit, damp the tilt, and settle the spin. Exercises WHFast plus "
    "REBOUNDx's tides_spin force and its spin differential equations. "
    "Verified bit-identical to the C REBOUNDx.",
    args=["62.83185307179586"], post="rxstate")

add("reboundx", "tides_spin_kozai",
    "tides_spin_kozai — a Kozai cycle with tides and relativity",
    "A planet with a distant stellar companion that periodically drives its "
    "orbit to high eccentricity. Uses the ADAPTIVE IAS15 integrator plus two "
    "REBOUNDx forces at once. Matching the C bit-for-bit here means both "
    "programs chose the identical sequence of thousands of adaptive steps.",
    args=["1000.0"], post="rxstate")

add("reboundx", "tides_spin_migration",
    "tides_spin_migration — migration-driven obliquity tides",
    "Two Earth-sized planets migrating through a gas disk while tides evolve "
    "their spins. Exercises two REBOUNDx forces simultaneously and a "
    "parameter changed in the middle of the run (migration is switched off "
    "half-way). Verified bit-identical to the C.",
    args=["62.83185307179586"], post="rxstate")

add("reboundx", "rebx_binary_roundtrip",
    "rebx_binary_roundtrip — REBOUNDx binary save/load",
    "Serialises a REBOUNDx state (forces and parameters of several types), "
    "reads it back into a fresh simulation, and checks every value returns "
    "with identical bits.",
    post="plain")


# Every cell needs a unique "id" (nbformat 4.5 and later). Without one,
# nbformat warns on every read and says it will become a hard error. The
# ids must be stable across regeneration so that re-running the generator
# does not produce a spurious diff in git, so they are simple counters
# rather than random strings.
_CELL_N = [0]


def _next_id():
    _CELL_N[0] += 1
    return "cell-%03d" % _CELL_N[0]


def code(src):
    return {"cell_type": "code", "id": _next_id(), "execution_count": None,
            "metadata": {}, "outputs": [],
            "source": src.splitlines(keepends=True)}


def md(src):
    return {"cell_type": "markdown", "id": _next_id(), "metadata": {},
            "source": src.splitlines(keepends=True)}


POSTS = {
    "plain": 'print("(see the program output above)")',
    "head": '''p = os.path.join(WORK, OUTFILE) if OUTFILE else None
if p and os.path.exists(p):
    lines = open(p).read().splitlines()
    print(f"{len(lines)} result lines; first 6:")
    print("\\n".join(lines[:6]))
''',
    "state": '''p = os.path.join(WORK, "state_rust_final.txt")
if os.path.exists(p):
    print(open(p).read())
''',
    "rxstate": '''name = EXAMPLE.replace("tides_spin_", "")
rust = os.path.join(WORK, f"state_{name}_rust.txt")
cref = os.path.join(WORK, f"state_{name}_c.txt")
if os.path.exists(rust):
    st = read_state(rust)
    print("--- final state (decoded from raw bits) ---")
    for p in st.get("p", []):
        print(f"  particle {p[0]}: x={p[1]:+.9e} y={p[2]:+.9e} z={p[3]:+.9e}")
    for o in st.get("Omega", []):
        if len(o) > 1 and not isinstance(o[1], str):
            mag = (o[1]**2 + o[2]**2 + o[3]**2) ** 0.5
            print(f"  spin {o[0]}: |Omega| = {mag:.9e}")
if os.path.exists(cref):
    print()
    compare(cref, rust)
else:
    print("\\n(no C reference dump present - build and run the C harness to compare)")
''',
    "sha": '''import hashlib
def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest().upper()
c  = os.path.join(WORK, "state_c_final.txt")
rs = os.path.join(WORK, "state_rust_final.txt")
print("Rust SHA-256:", sha(rs))
if os.path.exists(c):
    print("C    SHA-256:", sha(c))
    compare(c, rs)
else:
    print("(C reference dump not present; run porttest\\\\rebound_test.exe 400 to compare)")
''',
    "libm": '''c  = os.path.join(WORK, "libm_c.txt")
rs = os.path.join(WORK, "libm_rust.txt")
if os.path.exists(c) and os.path.exists(rs):
    lc, lr = open(c).read().splitlines(), open(rs).read().splitlines()
    from collections import Counter
    bad = Counter()
    for a, b in zip(lc, lr):
        if a != b:
            bad[a.split()[0]] += 1
    print(f"compared {len(lc)} samples")
    print("functions that differ:", dict(bad) if bad else "NONE - all bit-identical")
else:
    print("(run both libm_diff programs in porttest/ to compare)")
''',
    "powdiff": '''c  = os.path.join(WORK, "bs_pow_c.txt")
rs = os.path.join(WORK, "bs_pow_rust.txt")
if os.path.exists(c) and os.path.exists(rs):
    lc, lr = open(c).read().splitlines(), open(rs).read().splitlines()
    ulps = {}
    n = 0
    for a, b in zip(lc, lr):
        if a != b:
            n += 1
            d = abs(int(a.split()[-1], 16) - int(b.split()[-1], 16))
            ulps[d] = ulps.get(d, 0) + 1
    print(f"samples: {len(lc)}   mismatches: {n}  ({100*n/len(lc):.4f}%)")
    print("ULP distribution:", dict(sorted(ulps.items())) or "none")
else:
    print("(run both bs_pow_diff programs in porttest/ to compare)")
''',
    "freq": '''p = os.path.join(WORK, "frequency_rust.txt")
vals, mode = [], None
for line in open(p).read().splitlines():
    parts = line.split()
    if len(parts) == 3 and parts[1] == "ret":
        mode, vals = parts[0], []
        print(f"--- {mode} (ret {parts[2]}) ---")
    elif len(parts) == 2:
        vals.append(unbits(parts[1]))
        if len(vals) == 9:
            print("  frequencies:", [round(v, 6) for v in vals[0:3]])
            print("  amplitudes :", [round(v, 6) for v in vals[3:6]])
            print("  phases     :", [round(v, 6) for v in vals[6:9]])
''',
    "archive": '''import glob
for f in sorted(glob.glob(os.path.join(WORK, "archive_rust_*.bin"))):
    print(f"{os.path.basename(f)}: {os.path.getsize(f)} bytes")
p = os.path.join(WORK, "archive_state_rust.txt")
if os.path.exists(p):
    print(open(p).read())
''',
    "server": None,  # custom run cell
    "movetocom": '''c  = os.path.join(WORK, "movetocom_var_c.txt")
rs = os.path.join(WORK, "movetocom_var_rust.txt")
if os.path.exists(rs):
    print(open(rs).read())
if os.path.exists(c):
    compare(c, rs)
''',
}

SERVER_RUN = '''import subprocess, time, urllib.request
exe = os.path.join(CRATE, "target", "release", "examples", "server_test.exe")
proc = subprocess.Popen([exe], cwd=WORK, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, text=True)
time.sleep(2.0)
blob = urllib.request.urlopen("http://localhost:12873/simulation", timeout=10).read()
print(f"/simulation returned {len(blob)} bytes")
print(f"header: {blob[:32]!r}")
urllib.request.urlopen("http://localhost:12873/keyboard/81", timeout=10).read()
proc.wait(timeout=15)
print("server exited cleanly")
'''

SHEARING_PLOT = '''import matplotlib.pyplot as plt
p = os.path.join(WORK, "state_rust_final.txt")
xs, ys, rs = [], [], []
for line in open(p):
    parts = line.split()
    if len(parts) == 7 and parts[0].isdigit():
        xs.append(unbits(parts[1]))
        ys.append(unbits(parts[2]))
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(xs, ys, s=3, alpha=0.7)
ax.set_xlabel("x  [m]")
ax.set_ylabel("y  [m]")
ax.set_title(f"Shearing sheet: {len(xs)} ring particles after 400 steps")
ax.set_aspect("equal")
plt.show()
'''


def build(name, spec):
    crate_dir = REBOUND if spec["crate"] == "rebound" else REBOUNDX
    crate_name = "rebound_rust" if spec["crate"] == "rebound" else "reboundx_rust"
    work = os.path.join(crate_dir, "porttest")
    # Which example actually gets built and run. `runs` overrides the
    # notebook's own name for the one case where the example named in the
    # title never terminates on its own (shearing_sheet).
    runs = spec.get("runs") or name
    exe_rel = os.path.join("..", "target", "release", "examples", runs + ".exe")
    cmdline = exe_rel + "".join(" " + a for a in spec["args"])
    cells = [
        md(f"# {spec['title']}\n\n{spec['desc']}\n\n"
           f"This notebook is self-contained: it builds the example with cargo, "
           f"runs it, and shows the result. Everything it does can also be done "
           f"by hand in a terminal:\n\n"
           f"```\n"
           f"cd {crate_name}\n"
           f"cargo build --release --example {runs}\n"
           # Relative, because the block has already cd'd into the crate.
           # An absolute path here would bake this machine's directory
           # layout, and its username, into every generated notebook -
           # which is exactly what happened before, and was only caught
           # because regenerating reintroduced a path a text sweep had
           # already removed from the committed files.
           f"cd porttest\n"
           f"{cmdline}\n"
           f"```"),
        code(
            f'import os, subprocess\n'
            # Derived at run time from the notebook's own location so the
            # notebook works in any checkout, and carries no absolute path.
            f'NB_DIR  = os.getcwd()                       # <crate>/notebooks\n'
            f'ROOT    = os.path.dirname(os.path.dirname(NB_DIR))\n'
            f'CRATE   = os.path.join(ROOT, "{crate_name}")\n'
            f'WORK    = os.path.join(CRATE, "porttest")\n'
            f'if not os.path.exists(os.path.join(CRATE, "Cargo.toml")):\n'
            f'    raise SystemExit(\n'
            f'        "Could not find the crate. Run this notebook from "\n'
            f'        "the notebooks folder of a full checkout: " + CRATE)\n'
            f'# The example that is BUILT and RUN. It is usually the one the\n'
            f'# notebook is named after; where it differs (the stock\n'
            f'# shearing_sheet integrates forever by design) the terminating\n'
            f'# variant is used instead, and the note above says so.\n'
            f'EXAMPLE = "{runs}"\n'
            f'OUTFILE = None\n'
            f'os.makedirs(WORK, exist_ok=True)\n'
            f'res = subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE],\n'
            f'                     cwd=CRATE, capture_output=True, text=True)\n'
            f'print(res.stderr.strip()[-400:] or "build ok")'),
    ]
    post = spec["post"]
    if post in ("rxstate", "sha", "libm", "powdiff", "movetocom"):
        cells.append(code(DECODE + "\n" + COMPARE))
    elif post in ("freq",):
        cells.append(code(DECODE))

    if post == "server":
        cells.append(code(SERVER_RUN))
    else:
        arglist = "".join(f', "{a}"' for a in spec["args"])
        cells.append(code(
            f'exe = os.path.join(CRATE, "target", "release", "examples", EXAMPLE + ".exe")\n'
            f'res = subprocess.run([exe{arglist}], cwd=WORK, capture_output=True, text=True)\n'
            f'print(res.stdout)\n'
            f'if res.returncode != 0:\n'
            f'    print("STDERR:", res.stderr[-2000:])'))

    if post and post != "server" and POSTS.get(post):
        cells.append(code(POSTS[post]))
    if spec["plot"] == "shearing":
        cells.append(code(DECODE))
        cells.append(code(SHEARING_PLOT))

    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                           "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = os.path.join(OUT, f"{name}.ipynb")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(nb, fh, indent=1)
    return path


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for n, s in EX.items():
        print("wrote", build(n, s))
    print(f"\n{len(EX)} notebooks written to {OUT}")
