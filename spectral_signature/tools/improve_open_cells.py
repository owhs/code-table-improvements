"""Strengthen lower bounds for cells where the integer program hit its time limit: randomised min-degree greedy with restarts + (1,2)-swap local search on the conflict graph,
then (if the heuristic beats the incumbent) replace the stored code. Upper bounds are left as proven by the solver. Updates _build_data.json in place."""
import itertools, json, os, sys, time
import numpy as np

OUT = r"C:\projects\2026\reframe\publish\Spectral_codes_2026"
P_ = os.path.join(OUT, "_build_data.json")
D = json.load(open(P_))
cells = D["cells"]


def classes(n, w):
    def canon(c):
        best = None
        for r in range(n):
            for s in (1, -1):
                t = tuple(sorted((s * i + r) % n for i in c))
                if best is None or t < best:
                    best = t
        return best
    reps = sorted({canon(c) for c in itertools.combinations(range(n), w)})
    Wb = np.zeros((len(reps), n))
    for i, c in enumerate(reps):
        Wb[i, list(c)] = 1
    keys = {}
    for i, c in enumerate(reps):
        b = Wb[i].astype(int)
        keys.setdefault(tuple(int(sum(b[j] * b[(j + k) % n] for j in range(n))) for k in range(1, n // 2 + 1)), i)
    first = sorted(keys.values())
    P = np.abs(np.fft.rfft(Wb[first], axis=1)[:, 1:]) ** 2
    return [reps[i] for i in first], P


def heuristic(conf, start, seconds, rng):
    """maximum independent set on the conflict graph `conf` (bool NxN): greedy min-degree restarts + plateau (1,2)-swaps"""
    N = len(conf)
    best = list(start)
    t0 = time.time()
    deg0 = conf.sum(1)
    while time.time() - t0 < seconds:
        alive = np.ones(N, bool)
        sel = []
        while alive.any():
            idx = np.where(alive)[0]
            d = conf[np.ix_(idx, idx)].sum(1) + rng.random(len(idx)) * 1.5
            v = idx[int(np.argmin(d))]
            sel.append(v)
            alive &= ~conf[v]
            alive[v] = False
        # swap improvement: remove one selected vertex, try to add two
        improved = True
        while improved and time.time() - t0 < seconds:
            improved = False
            selset = set(sel)
            cnt = conf[:, sel].sum(1)                         # number of selected neighbours per vertex
            for v in rng.permutation(sel):
                cand = [u for u in np.where((cnt - conf[:, v]) == 0)[0] if u not in selset or u == v]
                cand = [u for u in cand if u != v]
                if len(cand) >= 2:
                    sub = conf[np.ix_(cand, cand)]
                    for a, b in itertools.combinations(range(len(cand)), 2):
                        if not sub[a, b]:
                            sel = [x for x in sel if x != v] + [cand[a], cand[b]]
                            improved = True
                            break
                if improved:
                    break
        if len(sel) > len(best):
            best = list(sel)
    return best


rng = np.random.default_rng(11)
changed = 0
for key, c in cells.items():
    if c.get("skipped"):
        continue
    open_th = [th for th, cd in c["codes"].items() if not cd["optimal"]]
    if not open_th:
        continue
    n, w = c["n"], c["w"]
    reps, P = classes(n, w)
    N = len(reps)
    Dm = np.abs(P[:, None, :] - P[None, :, :]).sum(2) / (w * w)
    for th in open_th:
        cd = c["codes"][th]
        if cd["upper"] is None:
            cd["upper"] = N                                   # trivial bound: at most one word per class
        conf = (Dm < float(th)) & ~np.eye(N, dtype=bool)
        cur = [reps.index(tuple(x)) for x in cd["words"]]
        new = heuristic(conf, cur, 40, rng)
        # validate before storing
        ok = all(Dm[i, j] >= float(th) - 1e-12 for i, j in itertools.combinations(new, 2))
        if ok and len(new) > cd["size"]:
            print(f"n={n} w={w} theta={th}: {cd['size']} -> {len(new)} (upper {cd['upper']})", flush=True)
            cd["size"] = len(new)
            cd["words"] = [list(reps[i]) for i in sorted(new)]
            if cd["size"] >= cd["upper"]:
                cd["optimal"] = True
            changed += 1
        else:
            print(f"n={n} w={w} theta={th}: kept {cd['size']} (upper {cd['upper']})", flush=True)
json.dump(D, open(P_, "w"))
print("updated", changed, "cells")
