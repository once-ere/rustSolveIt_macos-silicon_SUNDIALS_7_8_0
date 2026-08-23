#!/usr/bin/env python3
"""Smoke-test all thirteen live GUI servers (macOS/POSIX).

For each gui/<name>/server.py, in sequence:

  1. launch it (it owns a `posim --machine` child; POSIM_NO_BROWSER has
     no role here — the servers never open a browser themselves),
  2. GET  /            and check the page arrives with its canvas,
  3. GET  /api/state   and record t,
  4. POST /api/start, wait ~1.5 s of wall clock,
  5. GET  /api/state   and check t advanced (the engine really stepped),
  6. POST /api/stop, POST /api/reset, GET /api/state and check t == t0,
  7. terminate the server (its posim child dies with it).

Every check talks HTTP to the same fixed port the README documents, so a
pass here is exactly what a reader gets when they run
`python gui/<name>/server.py` and open the printed URL.

Stdlib only. Exit 0 only if every GUI passes every check.
"""
import json
import pathlib
import subprocess
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent

GUIS = {
    "piston_crankshaft": 8895,
    "rack_and_pinion": 8896,
    "gyroscope_gimbal": 8897,
    "cardan_compass": 8898,
    "universal_joint": 8899,
    "spinning_top": 8900,
    "ball_joint_chain": 8901,
    "cardan_gear": 8902,
    "rod_pendulum_chain": 8903,
    "double_pendulum_hinges": 8904,
    "tumbling_racket": 8905,
    "kepler_ellipse": 8906,
    "box_of_shapes": 8907,
}


def http(method, url, timeout=10):
    req = urllib.request.Request(url, method=method)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def wait_for(url, tries=60):
    for _ in range(tries):
        try:
            return http("GET", url)
        except OSError:
            time.sleep(0.5)
    raise RuntimeError(f"server never answered at {url}")


def state(port):
    return json.loads(http("GET", f"http://127.0.0.1:{port}/api/state"))


def run_one(name, port):
    proc = subprocess.Popen(
        [sys.executable, str(ROOT / "gui" / name / "server.py")],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=ROOT,
    )
    try:
        page = wait_for(f"http://127.0.0.1:{port}/").decode("utf-8")
        assert "<canvas" in page, f"{name}: page has no canvas"
        s0 = state(port)
        t0 = s0.get("t", 0.0)
        http("POST", f"http://127.0.0.1:{port}/api/start")
        time.sleep(1.5)
        s1 = state(port)
        assert s1.get("t", 0.0) > t0, f"{name}: t did not advance ({t0} -> {s1.get('t')})"
        http("POST", f"http://127.0.0.1:{port}/api/stop")
        http("POST", f"http://127.0.0.1:{port}/api/reset")
        time.sleep(0.5)
        s2 = state(port)
        assert abs(s2.get("t", -1.0) - t0) < 1e-12, \
            f"{name}: reset did not return to t = {t0} (got {s2.get('t')})"
        return s1.get("t", 0.0)
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=10)


def main():
    failed = []
    for name, port in GUIS.items():
        try:
            t = run_one(name, port)
            print(f"ok    {name} (port {port}): page served, ran to t = {t:.3f}, reset to 0")
        except Exception as e:
            print(f"FAIL  {name} (port {port}): {e}")
            failed.append(name)
    print(f"\n{len(GUIS) - len(failed)} of {len(GUIS)} GUIs pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
