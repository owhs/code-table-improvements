"""Detached hunt driver: run tools/qary_search.py over all open cells of data/brouwer/{ternary,mixed23}.json,
softest first (lower/upper ascending; truly mixed + ternary), at most --procs concurrent single-threaded workers.
Each worker compares with the parsed live lower bound and writes runs/qary/HIT_<table>_<params>.json on a beat.
    python tools/qary_hunt.py [--procs 6] [--secs 900] [--rounds 3]
Progress: runs/qary/hunt_progress.txt; per-cell logs runs/qary/hunt_<cell>_r<round>.txt"""
import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs", "qary")


def cells():
    out = []
    for r in json.load(open(os.path.join(ROOT, "data", "brouwer", "ternary.json"), encoding="utf-8")):
        if not r["exact"] and r["n"] >= 7:
            out.append((r["lower"] / r["upper"], f"T:{r['n']},{r['d']}", r["lower"]))
    for r in json.load(open(os.path.join(ROOT, "data", "brouwer", "mixed23.json"), encoding="utf-8")):
        if r["lower"] and not r["exact"] and r["m_bin"] > 0 and r["n_ter"] > 0:
            out.append((r["lower"] / r["upper"], f"M:{r['m_bin']},{r['n_ter']},{r['d']}", r["lower"]))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--secs", type=float, default=900)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--seed0", type=int, default=100)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    jobs = [(rd, c) for rd in range(a.rounds) for c in cells()]
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
            rd, (soft, spec, lower) = jobs.pop(0)
            tag = spec.replace(":", "_").replace(",", "_") + f"_r{rd}"
            refill = min(300, a.secs / 4)
            cmd = [sys.executable, os.path.join(ROOT, "tools", "qary_search.py"), "--cell", spec,
                   "--secs", str(a.secs - refill), "--refill", str(refill), "--seed", str(a.seed0 + rd),
                   "--orb", "30000", "--vmax", "8000", "--out", OUT]
            log = open(os.path.join(OUT, f"hunt_{tag}.txt"), "w")
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=env), tag))
            prog.write(f"{time.ctime()} start {tag} (live lower {lower}, softness {soft:.2f})\n")
        time.sleep(5)
    prog.write(f"{time.ctime()} all done\n")


if __name__ == "__main__":
    main()
