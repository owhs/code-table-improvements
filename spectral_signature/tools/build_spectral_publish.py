"""Build the data for publish/Spectral_codes_2026 (run: python -u tools/build_spectral_publish.py).  Writes codes.json, CSV tables, code and matrix files.
Everything is recomputed from scratch here (class enumeration, spectra, integer programs); the weight-resolved class COUNTS for n up to 30 come from
theory/homo-seq_fast_out_n*.json (exact integer autocorrelation keys) and are cross-checked against this script's own enumeration for every cell it enumerates."""
import itertools, json, os, sys, time, glob, math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

OUT = r"C:\projects\2026\reframe\publish\Spectral_codes_2026"
THETAS = (0.05, 0.10, 0.20, 0.40)
MAXCLASSES = 3500
LIMIT = 30
rng = np.random.default_rng(7)

counts = {}
for f in glob.glob("theory/homo-seq_fast_out_n*.json"):
    j = json.load(open(f))
    counts[j["n"]] = {r["w"]: (r["bracelets"], r["classes"]) for r in j["per_weight"]}


def enum_classes(n, w):
    def canon(c):
        best = None
        for r in range(n):
            for s in (1, -1):
                t = tuple(sorted((s * i + r) % n for i in c))
                if best is None or t < best:
                    best = t
        return best
    reps = sorted({canon(c) for c in itertools.combinations(range(n), w)})
    br = len(reps)
    Wb = np.zeros((br, n))
    for i, c in enumerate(reps):
        Wb[i, list(c)] = 1
    P = np.abs(np.fft.rfft(Wb, axis=1)[:, 1:]) ** 2
    # exact class keys via integer autocorrelation
    keys = {}
    for i, c in enumerate(reps):
        b = Wb[i].astype(int)
        key = tuple(int(sum(b[j] * b[(j + k) % n] for j in range(n))) for k in range(1, n // 2 + 1))
        keys.setdefault(key, i)
    first = sorted(keys.values())
    return [reps[i] for i in first], P[first], br


def solve(D, theta):
    N = len(D)
    bad = np.argwhere(np.triu(D < theta, 1))
    cons = []
    if len(bad):
        A = coo_matrix((np.ones(2 * len(bad)), (np.repeat(np.arange(len(bad)), 2), bad.ravel())), shape=(len(bad), N)).tocsr()
        cons = [LinearConstraint(A, -np.inf, 1)]
    r = milp(-np.ones(N), constraints=cons, integrality=np.ones(N), bounds=Bounds(0, 1), options=dict(time_limit=LIMIT, disp=False))
    chosen = sorted(np.where(r.x > 0.5)[0].tolist()) if r.x is not None else []
    ub = int(np.floor(-r.mip_dual_bound + 1e-6)) if getattr(r, "mip_dual_bound", None) is not None else None
    return chosen, (r.status == 0), ub


cells = {}
t0 = time.time()
for n in range(6, 21):
    for w in range(3, n // 2 + 1):
        reps, P, br = enum_classes(n, w)
        N = len(reps)
        ref = counts.get(n, {}).get(w)
        if ref is not None and (ref[0] != br or ref[1] != N):
            print(f"!! MISMATCH n={n} w={w}: own (bracelets {br}, classes {N}) vs agent {ref}", flush=True)
        if N > MAXCLASSES:
            cells[f"{n},{w}"] = dict(n=n, w=w, bracelets=br, classes=N, skipped=True)
            print(f"n={n:2d} w={w:2d} classes {N:6d} skipped (> {MAXCLASSES})", flush=True)
            continue
        D = np.abs(P[:, None, :] - P[None, :, :]).sum(2) / (w * w)
        iu = np.triu_indices(N, 1)
        dmin = float(D[iu].min()) if N > 1 else None
        rec = dict(n=n, w=w, bracelets=br, classes=N, delta_min=dmin, codes={})
        for th in THETAS:
            if dmin is not None and th <= dmin * (1 - 1e-9):
                chosen, opt, ub = list(range(N)), True, N
            else:
                chosen, opt, ub = solve(D, th)
            rec["codes"][str(th)] = dict(size=len(chosen), optimal=bool(opt), upper=ub if not opt else len(chosen), words=[list(reps[i]) for i in chosen])
        cells[f"{n},{w}"] = rec
        print(f"n={n:2d} w={w:2d} bracelets {br:5d} classes {N:5d} delta_min {dmin if dmin is None else round(dmin,3)} | " +
              " ".join(f"{rec['codes'][str(t)]['size']}{'' if rec['codes'][str(t)]['optimal'] else '(<=%s)' % rec['codes'][str(t)]['upper']}" for t in THETAS) +
              f"   [{time.time()-t0:.0f}s]", flush=True)

# ------------- noise demos
def spec_of(code, n):
    M = np.zeros((len(code), n))
    for i, c in enumerate(code):
        M[i, list(c)] = 1
    return np.abs(np.fft.rfft(M, axis=1)[:, 1:]) ** 2, M


def ser(P, M, n, sigma, T=3000):
    err = 0
    for _ in range(T):
        k = int(rng.integers(len(P)))
        y = M[k] + sigma * rng.standard_normal(n)
        Py = np.abs(np.fft.rfft(y)[1:]) ** 2
        err += int(np.argmin(np.abs(P - Py).sum(1)) != k)
    return err / T


demos = []
for (n, w, K) in ((16, 6, 12), (13, 4, 13), (18, 6, 16)):
    cell = cells[f"{n},{w}"]
    best_th = None
    for th in THETAS[::-1]:
        if cell["codes"][str(th)]["size"] >= K:
            best_th = th
            break
    spec = [tuple(x) for x in cell["codes"][str(best_th)]["words"][:K]]
    allw = list(itertools.combinations(range(n), w))
    ham = []
    for dmin in range(2 * w, 1, -2):
        ham = []
        for i in rng.permutation(len(allw)):
            c = allw[i]
            if all(len(set(c) ^ set(h)) >= dmin for h in ham):
                ham.append(c)
            if len(ham) >= K:
                break
        if len(ham) >= K:
            break
    base = spec[0]
    cyc = [tuple(sorted((i + r) % n for i in base)) for r in range(min(K, n))]
    dh = min(len(set(a) ^ set(b)) for a, b in itertools.combinations(ham, 2))
    ps, ms = spec_of(spec, n); ph, mh = spec_of(ham, n); pc, mc = spec_of(cyc, n)
    dsp = min(np.abs(a - b).sum() for a, b in itertools.combinations(ps, 2)) / w ** 2
    dhs = min(np.abs(a - b).sum() for a, b in itertools.combinations(ph, 2)) / w ** 2
    rows = []
    for sg in (0.1, 0.2, 0.3, 0.5):
        rows.append(dict(sigma=sg, spectral=ser(ps, ms, n, sg), hamming=ser(ph, mh, n, sg), cyclic=ser(pc, mc, n, sg)))
    demos.append(dict(n=n, w=w, K=K, threshold_used=best_th, spectral_min_dist=float(dsp), hamming_min_hamming=int(dh), hamming_min_spectral=float(dhs),
                      chance=1 - 1 / K, rows=rows, cyclic_size=len(cyc)))
    print(f"demo n={n} w={w} K={K}: spectral min d {dsp:.3f} | Hamming code d_H={dh}, its min spectral d {dhs:.3f}", flush=True)
    for r in rows:
        print("   ", r, flush=True)

os.makedirs(OUT, exist_ok=True)
json.dump(dict(counts={str(n): {str(w): v for w, v in c.items()} for n, c in counts.items()}, cells=cells, demos=demos, thetas=list(THETAS)),
          open(os.path.join(OUT, "_build_data.json"), "w"))
print("build data written", flush=True)
