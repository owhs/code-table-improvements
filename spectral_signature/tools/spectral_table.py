"""Table for the new category: spectral-signature codes. Normalised spectral distance  d(x,y) = sum_{k=1}^{floor(n/2)} | |X_x(k)|^2 - |X_y(k)|^2 | / w^2.
For each (n,w): number of distinct power-spectrum classes (the ceiling), and the maximum code size at thresholds d >= 0.05, 0.10, 0.20, 0.40
(proven optimal by integer programming where the status says so, else best found with an upper bound)."""
import itertools, json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix


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
    P = np.abs(np.fft.rfft(Wb, axis=1)[:, 1:]) ** 2
    _, first = np.unique(np.round(P, 6), axis=0, return_index=True)
    return P[np.sort(first)]


def solve(D, theta, limit):
    N = len(D)
    bad = np.argwhere(np.triu(D < theta, 1))
    cons = []
    if len(bad):
        A = coo_matrix((np.ones(2 * len(bad)), (np.repeat(np.arange(len(bad)), 2), bad.ravel())), shape=(len(bad), N)).tocsr()
        cons = [LinearConstraint(A, -np.inf, 1)]
    r = milp(-np.ones(N), constraints=cons, integrality=np.ones(N), bounds=Bounds(0, 1), options=dict(time_limit=limit, disp=False))
    size = int(round(-r.fun)) if r.x is not None else 0
    ub = int(np.floor(-r.mip_dual_bound + 1e-6)) if getattr(r, "mip_dual_bound", None) is not None else None
    return size, (r.status == 0), ub


res = []
print(f"{'n':>2} {'w':>2} {'classes':>8} | " + " ".join(f"d>={t:<5}" for t in (0.05, 0.10, 0.20, 0.40)), flush=True)
for n, ws in ((10, (3, 5)), (12, (4, 6)), (13, (4, 6)), (14, (5, 7)), (15, (5, 7)), (16, (5, 8)), (18, (6, 9))):
    for w in ws:
        P = classes(n, w)
        N = len(P)
        if N > 3500:
            continue
        D = np.abs(P[:, None, :] - P[None, :, :]).sum(2) / (w * w)
        row = []
        for th in (0.05, 0.10, 0.20, 0.40):
            s, opt, ub = solve(D, th, 40)
            row.append((s, opt, ub))
        res.append(dict(n=n, w=w, classes=N, cells=row))
        print(f"{n:2d} {w:2d} {N:8d} | " + " ".join((f"{s}" if opt else f"{s}(<= {ub})").ljust(8) for s, opt, ub in row), flush=True)
json.dump(res, open("research/spectral_table.json", "w"), indent=1)
