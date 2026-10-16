"""Run tools/cw_improve.py over Brouwer's downloaded record code files (runs/brouwer_codes/a<n>.<d>.<w>.<size>[H]), a few at a time.
Progress file runs/cwimprove/progress.txt uses the same 'start:' / 'done' format as the other hunt drivers (for the watcher).
    python tools/cw_improve_all.py [--procs 3]"""
import argparse, glob, os, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs", "cwimprove")
ap = argparse.ArgumentParser()
ap.add_argument("--procs", type=int, default=3)
a = ap.parse_args()
os.makedirs(OUT, exist_ok=True)
jobs = []
for f in glob.glob(os.path.join(ROOT, "runs", "brouwer_codes", "a*")):
    m = re.match(r"a(\d+)\.(\d+)\.(\d+)\.(\d+)", os.path.basename(f))
    if m:
        n, d, w, size = (int(x) for x in m.groups())
        jobs.append((n, d, w, f))
jobs.sort(key=lambda j: (j[0], j[2]))  # smaller scans first
prog = open(os.path.join(OUT, "progress.txt"), "a", buffering=1)
prog.write(f"{time.ctime()} start: {len(jobs)} jobs, procs {a.procs}\n")
running = []
env = dict(os.environ, NUMBA_NUM_THREADS="1")
while jobs or running:
    for p, tag in list(running):
        if p.poll() is not None:
            running.remove((p, tag))
            prog.write(f"{time.ctime()} done {tag} rc={p.returncode}\n")
    while jobs and len(running) < a.procs:
        n, d, w, f = jobs.pop(0)
        tag = f"A_{n}_{d}_{w}"
        log = open(os.path.join(OUT, f"log_{tag}.txt"), "w")
        running.append((subprocess.Popen([sys.executable, "-I", os.path.join(ROOT, "tools", "cw_improve.py"), f, str(n), str(d), str(w)],
                                         stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=env), tag))
        prog.write(f"{time.ctime()} start {tag}\n")
    time.sleep(5)
prog.write(f"{time.ctime()} all done\n")
