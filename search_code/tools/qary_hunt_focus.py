"""Focused hunt: open truly-mixed cells of data/brouwer/mixed23.json whose live lower bound is a product-type number
2^a * 3^b >= 24 (likely never searched directly), small-group monomial families only (mono2, mono), larger orbit target so that
small groups are allowed, larger clique graph cap. Same progress-file format as qary_hunt.py so the watcher works.
    python tools/qary_hunt_focus.py [--procs 8] [--secs 600] [--seed0 300]
Output: runs/qary_focus/ (HIT_*.json, hunt_progress.txt, one log per job)."""
import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs", "qary_focus")


def smooth(x):
    for p in (2, 3):
        while x % p == 0 and x > 1:
            x //= p
    return x == 1


def cells():
    out = []
    for r in json.load(open(os.path.join(ROOT, "data", "brouwer", "mixed23.json"), encoding="utf-8")):
        if r["lower"] and not r["exact"] and r["m_bin"] > 0 and r["n_ter"] > 0 and r["lower"] >= 24 and smooth(r["lower"]):
            out.append((r["lower"] / r["upper"], f"M:{r['m_bin']},{r['n_ter']},{r['d']}", r["lower"]))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--secs", type=float, default=600)
    ap.add_argument("--seed0", type=int, default=300)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    jobs = [(fam, c) for c in cells() for fam in ("mono2", "mono")]
    env = dict(os.environ, NUMBA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    running = []
    prog = open(os.path.join(OUT, "hunt_progress.txt"), "a", buffering=1)
    prog.write(f"{time.ctime()} start: {len(jobs)} jobs, procs {a.procs}, secs {a.secs}\n")
    while jobs or running:
        for p, tag in list(running):
            if p.poll() is not None:
                running.remove((p, tag))
                prog.write(f"{time.ctime()} done {tag} rc={p.returncode}\n")
        while jobs and len(running) < a.procs:
            fam, (soft, spec, lower) = jobs.pop(0)
            tag = spec.replace(":", "_").replace(",", "_") + f"_{fam}"
            cmd = [sys.executable, os.path.join(ROOT, "tools", "qary_search.py"), "--cell", spec,
                   "--secs", str(a.secs - 30), "--refill", "30", "--seed", str(a.seed0 + len(running)),
                   "--fam", fam, "--orb", "250000", "--vmax", "14000", "--out", OUT]
            log = open(os.path.join(OUT, f"hunt_{tag}.txt"), "w")
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=env), tag))
            prog.write(f"{time.ctime()} start {tag} (live lower {lower}, softness {soft:.2f})\n")
        time.sleep(5)
    prog.write(f"{time.ctime()} all done\n")


if __name__ == "__main__":
    main()
