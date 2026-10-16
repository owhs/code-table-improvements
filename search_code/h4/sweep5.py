"""Sweep 3.0 for constant-weight codes A(n,d,w): the sweep4 pipeline (exact orbit pool, compiled tabu clique, exact clique on small
pools, certify) with a much wider group library, because sweep4 (affine groups of a prime p, n = p..p+2) could not reach the
records for w >= 13 (e.g. A(32,10,13): 1302 vs 1984).
New groups: block cyclic groups. n = c*b + k (k fixed points, c blocks of length b); generators: shift every block together,
multiply positions by a unit m of Z_b (optionally also rotate the blocks, giving Z_b x Z_c); b composite is allowed.
Groups are ranked by closeness of |G| to the order that makes the orbit pool about POOL_TARGET orbits (so pools are complete when
possible), and the best few are run per cell.
    python -m h4.sweep5 <table json> <workers> <minutes> <nmin> <nmax> <dmin> <dmax> [wmin wmax]
Output: runs/sweep5/results.jsonl, runs/sweep5/hit_<n>_<d>_<w>.json (certified improvements)."""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
POOL_TARGET = int(os.environ.get("POOL_TARGET", 20000))
MAX_GROUPS = int(os.environ.get("MAX_GROUPS", 10))
MAXCOMB = float(os.environ.get("MAXCOMB", 1e18))   # skip cells with C(n,w) above this (orbit enumeration cost)
MAXJOBS = int(os.environ.get("MAXJOBS", 10 ** 6))  # keep only the softest MAXJOBS cells (lower/upper ascending)


def units_cyclic_subgroups(b):
    """one generator per distinct cyclic subgroup of the unit group of Z_b (with its order)."""
    seen, out = set(), []
    for m in range(1, b):
        if math.gcd(m, b) != 1:
            continue
        orb, x = [], 1
        while True:
            x = (x * m) % b
            orb.append(x)
            if x == 1:
                break
        key = frozenset(orb)
        if key not in seen:
            seen.add(key)
            out.append((m, len(orb)))
    return out


def block_group(n, c, b, k, m, rotate):
    """generators (as permutations of range(n)) of the block group; order is computed by the closure later."""
    shift, mult = list(range(n)), list(range(n))
    for blk in range(c):
        for i in range(b):
            shift[blk * b + i] = blk * b + (i + 1) % b
            mult[blk * b + i] = blk * b + (m * i) % b
    gens = [shift, mult]
    if rotate and c > 1:
        rot = list(range(n))
        for blk in range(c):
            for i in range(b):
                rot[blk * b + i] = ((blk + 1) % c) * b + i
        gens.append(rot)
    return gens


def groups_for(n, w):
    target = max(8, math.comb(n, w) / POOL_TARGET)
    cand = []
    for k in range(0, 4):
        m_pts = n - k
        for c in range(1, 5):
            if m_pts % c:
                continue
            b = m_pts // c
            if b < 7:
                continue
            for mm, h in units_cyclic_subgroups(b):
                for rotate in ((False, True) if c > 1 else (False,)):
                    order = b * h * (c if rotate else 1)
                    if order < 24 or order > 20000:
                        continue
                    cand.append((abs(math.log(order / target)), order, {
                        "type": "generators", "n": n, "gens": block_group(n, c, b, k, mm, rotate),
                        "tag": f"c{c}b{b}k{k}m{mm}{'R' if rotate else ''} order {order}"}))
    cand.sort(key=lambda x: x[0])
    seen, out = set(), []
    for _, order, g in cand:
        key = (order, g["tag"].split(" order")[0])
        if key in seen:
            continue
        seen.add(key)
        out.append(g)
        if len(out) >= MAX_GROUPS:
            break
    return out


def run_entry(job):
    from engine.problems.cwbase import make
    from h4 import clique as CL
    from h4 import cw_ops as O
    n, d, w, rec = job["n"], job["d"], job["w"], job["record"]
    lam = w - d // 2
    prob = make(n, d, w, rec)
    t0 = time.time()
    best, best_words, log = 0, None, []
    for g in groups_for(n, w):
        tag = g.pop("tag")
        if time.time() - t0 > job["minutes"] * 60:
            log.append(f"{tag}: skipped (time)")
            continue
        try:
            orbs, complete = O.orbit_pool(n, w, lam, {"group": g, "method": "exact", "max_iter": 2 * 10 ** 9})
        except Exception as e:
            log.append(f"{tag}: error {e!r}"[:160])
            continue
        if not orbs:
            log.append(f"{tag}: 0 orbits")
            continue
        if len(orbs) > POOL_TARGET + 5000:  # random sample, weighted toward long orbits, instead of keeping only the longest
            rng = np.random.default_rng(1)
            lens = np.array([len(o) for o in orbs], dtype=float)
            pick = rng.choice(len(orbs), POOL_TARGET + 5000, replace=False, p=lens / lens.sum())
            orbs = [orbs[i] for i in pick]
            complete = False
        orbs_list, C = O.compat([orbs], lam, same_group=True)
        wts = np.array([len(o) for o in orbs_list], dtype=np.int64)
        idx, size = CL.tabu(C, wts, 100000, 0)
        how = "tabu"
        if len(orbs_list) <= 8000 and C.sum() // 2 <= 2_000_000:
            idx2, size2, how2 = O.select(orbs_list, C, {"solver": "exact"}, 1.0)
            if size2 >= size:
                idx, size, how = idx2, size2, ("exact" if complete else "exact-on-partial-pool")
        log.append(f"{tag}: {len(orbs_list)} orbits{'' if complete else ' (partial)'}, {how} {size}")
        if size > best:
            best = size
            best_words = [[b for b in range(n) if (int(m) >> b) & 1] for i in idx for m in orbs_list[i]]
    res = {"n": n, "d": d, "w": w, "record": rec, "best": best, "seconds": round(time.time() - t0), "log": log}
    if best_words and best > rec:
        cert = prob["certify"]({"words": best_words})
        res["certified"] = bool(cert["ok"])
        os.makedirs(os.path.join(ROOT, "runs", "sweep5"), exist_ok=True)
        json.dump({"n": n, "d": d, "w": w, "record": rec, "size": best, "words": best_words, "log": log},
                  open(os.path.join(ROOT, "runs", "sweep5", f"hit_{n}_{d}_{w}.json"), "w"))
    return res


from functools import lru_cache


@lru_cache(maxsize=None)
def johnson(n, d, w):
    """Johnson bound A(n,d,w) <= floor(n/w * A(n-1,d,w-1)); the tables use it where no better upper bound is given."""
    if 2 * w < d:
        return 1
    if 2 * w == d:
        return n // w
    return (n * johnson(n - 1, d, w - 1)) // w


def main():
    table = json.load(open(sys.argv[1], encoding="utf-8-sig"))
    for r in table:  # cells without a tabulated upper bound: use the Johnson bound
        if r["upper"] is None:
            r["upper"] = johnson(r["n"], r["d"], min(r["w"], r["n"] - r["w"]))
    workers, minutes = int(sys.argv[2]), float(sys.argv[3])
    nmin, nmax, dmin, dmax = (int(x) for x in sys.argv[4:8])
    wmin, wmax = (int(sys.argv[8]), int(sys.argv[9])) if len(sys.argv) > 9 else (3, 99)
    jobs = [{"n": r["n"], "d": r["d"], "w": r["w"], "record": r["lower"], "minutes": minutes}
            for r in table if nmin <= r["n"] <= nmax and dmin <= r["d"] <= dmax and wmin <= r["w"] <= wmax
            and 2 * r["w"] <= r["n"] and r["upper"] is not None and r["lower"] < r["upper"]]
    up = {(r["n"], r["d"], r["w"]): r["upper"] for r in table}
    jobs = [j for j in jobs if math.comb(j["n"], j["w"]) <= MAXCOMB]
    jobs.sort(key=lambda j: j["record"] / up[(j["n"], j["d"], j["w"])])
    jobs = jobs[:MAXJOBS]
    os.makedirs(os.path.join(ROOT, "runs", "sweep5"), exist_ok=True)
    out = os.path.join(ROOT, "runs", "sweep5", "results.jsonl")
    print(f"{len(jobs)} entries", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        futs = {ex.submit(run_entry, j): j for j in jobs}
        for f in as_completed(futs):
            j = futs[f]
            try:
                r = f.result()
            except Exception as e:
                r = {**j, "error": repr(e)[:200]}
            open(out, "a").write(json.dumps(r) + "\n")
            flag = "  <<< IMPROVEMENT" + (" CERTIFIED" if r.get("certified") else "") if r.get("best", 0) > j["record"] else ""
            print(f"A({j['n']},{j['d']},{j['w']}): best {r.get('best')} vs record {j['record']} [{r.get('seconds')}s]{flag}", flush=True)


if __name__ == "__main__":
    main()
