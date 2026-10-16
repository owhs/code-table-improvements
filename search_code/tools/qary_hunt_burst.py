"""Burst hunt: many SHORT randomised runs per cell instead of a few long ones.
Reason: in the long hunts most winning codes appeared in the first minutes of a run (7 of 12 within ~2 min, hits then go quiet),
so compute is better spent on many fresh starts than on polishing one. Every job draws a new seed, a random group family
(biased to mono/mono2, which produced 11 of the first 17 hits) and a random orbit target (30k / 100k / 250k, which sets the
group order range), and stops as soon as it beats the best size known so far for the cell (--stop_at).
Open truly-mixed cells of data/brouwer/mixed23.json, shuffled each round.
    python tools/qary_hunt_burst.py [--procs 8] [--secs 120] [--rounds 2] [--seed0 1000]
Output: runs/qary_burst (HIT_*.json when the live table bound is beaten, hunt_progress.txt in the usual format)."""
import argparse, glob, json, os, random, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs", "qary_burst")
FAMS = ["mono"] * 7 + ["mono2"] * 7 + ["mono+K"] * 2 + ["struct"] * 2 + ["any"] * 2
ORBS = [30000, 100000, 250000]


def ours():
    best = {}
    for f in glob.glob(os.path.join(ROOT, "runs", "qary_hits", "HIT_mixed23_*.json")):
        try:
            j = json.load(open(f))
            best[(j["m_bin"], j["n_ter"], j["d"])] = len(j["words"])
        except Exception:
            pass
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--secs", type=float, default=120)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--seed0", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    known = ours()
    cells = [(r["m_bin"], r["n_ter"], r["d"], r["lower"]) for r in
             json.load(open(os.path.join(ROOT, "data", "brouwer", "mixed23.json"), encoding="utf-8"))
             if r["lower"] and not r["exact"] and r["m_bin"] > 0 and r["n_ter"] > 0]
    jobs = []
    for rd in range(a.rounds):
        rng = random.Random(a.seed0 + rd)
        order = cells[:]
        rng.shuffle(order)
        for i, c in enumerate(order):
            jobs.append((rd, i, c, rng.choice(FAMS), rng.choice(ORBS)))
    prog = open(os.path.join(OUT, "hunt_progress.txt"), "a", buffering=1)
    prog.write(f"{time.ctime()} start: {len(jobs)} jobs, procs {a.procs}, secs {a.secs}\n")
    env = dict(os.environ, NUMBA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    running = []
    while jobs or running:
        for p, tag in list(running):
            if p.poll() is not None:
                running.remove((p, tag))
                prog.write(f"{time.ctime()} done {tag} rc={p.returncode}\n")
        while jobs and len(running) < a.procs:
            rd, i, (m, n, d, lower), fam, orb = jobs.pop(0)
            tag = f"M_{m}_{n}_{d}_r{rd}_{fam.replace('+', 'K')}_{orb // 1000}k"
            stop = max(lower, known.get((m, n, d), 0)) + 1
            cmd = [sys.executable, os.path.join(ROOT, "tools", "qary_search.py"), "--cell", f"M:{m},{n},{d}",
                   "--secs", str(a.secs - 15), "--refill", "15", "--seed", str(a.seed0 + rd * 100000 + i),
                   "--fam", fam, "--orb", str(orb), "--vmax", "10000", "--stop_at", str(stop), "--out", OUT]
            log = open(os.path.join(OUT, f"hunt_{tag}.txt"), "w")
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=env), tag))
            prog.write(f"{time.ctime()} start {tag}\n")
        time.sleep(3)
    prog.write(f"{time.ctime()} all done\n")


if __name__ == "__main__":
    main()
