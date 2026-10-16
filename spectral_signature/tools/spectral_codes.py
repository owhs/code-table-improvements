"""NEW CATEGORY: spectral-signature codes (codes for intensity-only / phase-blind receivers).
Words: binary length n, weight w.  Receiver sees only P(x)_k = |sum_j x_j e^{-2 pi i jk/n}|^2, k = 1..floor(n/2) (+ noise).
Words in the same dihedral class (shift / reflection) have IDENTICAL P, so a code can hold at most one word per class, and homometric classes also collide.
Spectral distance  D(x,y) = sum_k |P(x)_k - P(y)_k|.   A (n,w,theta) spectral code = set of classes with pairwise D >= theta.
Exact maximum size = maximum clique of the graph on classes (edge iff D >= theta).   Then a USE demo: symbol error rate of a phase-blind receiver
for (a) a spectral code, (b) a Hamming-distance code of the same size, (c) the cyclic-shift family of one word."""
import sys, itertools, json, time
import numpy as np
import networkx as nx

n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
w = int(sys.argv[2]) if len(sys.argv) > 2 else 6


def canon(c):
    best = None
    for r in range(n):
        for s in (1, -1):
            t = tuple(sorted((s * i + r) % n for i in c))
            if best is None or t < best:
                best = t
    return best


t0 = time.time()
reps = sorted({canon(c) for c in itertools.combinations(range(n), w)})
Wb = np.zeros((len(reps), n))
for i, c in enumerate(reps):
    Wb[i, list(c)] = 1
P = np.abs(np.fft.rfft(Wb, axis=1)[:, 1:]) ** 2          # (classes, n//2)
# merge homometric classes (identical P)
key = np.round(P, 6)
_, first = np.unique(key, axis=0, return_index=True)
first = np.sort(first)
reps = [reps[i] for i in first]; P = P[first]; Wb = Wb[first]
N = len(reps)
print(f"n={n} w={w}: {comb_total if (comb_total:=0) else ''}", end="")
print(f"{N} power-spectrum classes (the maximum alphabet size for ANY phase-blind separation > 0)  [{time.time()-t0:.1f}s]", flush=True)
D = np.abs(P[:, None, :] - P[None, :, :]).sum(2)
iu = np.triu_indices(N, 1)
dvals = D[iu]
print(f"spectral distance: min {dvals.min():.3f}, median {np.median(dvals):.2f}, max {dvals.max():.2f}")


LAST = {}


def max_code(theta, limit=60):
    """maximum set of classes with pairwise D >= theta: integer program over conflict pairs (HiGHS); LAST records proof status and bound"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    bad = np.argwhere(np.triu(D < theta, 1))
    if len(bad):
        rows = np.repeat(np.arange(len(bad)), 2)
        A = coo_matrix((np.ones(2 * len(bad)), (rows, bad.ravel())), shape=(len(bad), N)).tocsr()
        cons = [LinearConstraint(A, -np.inf, 1)]
    else:
        cons = []
    r = milp(-np.ones(N), constraints=cons, integrality=np.ones(N), bounds=Bounds(0, 1), options=dict(time_limit=limit, disp=False))
    LAST["optimal"] = (r.status == 0)
    LAST["bound"] = int(np.floor(-r.mip_dual_bound + 1e-6)) if getattr(r, "mip_dual_bound", None) is not None else None
    return sorted(np.where(r.x > 0.5)[0].tolist()) if r.x is not None else []


print("\nexact maximum spectral code size by threshold (fraction of the median distance):")
table = {}
for frac in (0.25, 0.5, 0.75, 1.0, 1.25):
    th = frac * np.median(dvals)
    t1 = time.time()
    cl = max_code(th)
    table[frac] = (th, len(cl))
    tag = "proven optimal" if LAST["optimal"] else f"best found, upper bound {LAST['bound']}"
    print(f"  theta = {th:7.2f} ({frac:.2f} x median): max size {len(cl):4d} ({tag})   [{time.time()-t1:.1f}s]", flush=True)

# ---------------- use demo: phase-blind receiver SER
K = 12
th = 0.75 * np.median(dvals)
cl = max_code(th)
spec_code = cl[:K] if len(cl) >= K else cl
K = len(spec_code)
# Hamming code of size K: greedy min-distance words (all words, not classes), pick a code with large Hamming distance
all_words = [c for c in itertools.combinations(range(n), w)]
rng = np.random.default_rng(1)
perm = rng.permutation(len(all_words))
ham = []
for dmin in range(2 * w, 1, -2):
    ham = []
    for i in perm:
        c = all_words[i]
        if all(len(set(c) ^ set(h)) >= dmin for h in ham):
            ham.append(c)
        if len(ham) >= K:
            break
    if len(ham) >= K:
        break
ham = ham[:K]
dmin_h = min(len(set(a) ^ set(b)) for a, b in itertools.combinations(ham, 2))
# cyclic family of one word
base = reps[0]
cyc = [tuple(sorted((i + r) % n for i in base)) for r in range(min(K, n))]


def spec_of(code):
    M = np.zeros((len(code), n))
    for i, c in enumerate(code):
        M[i, list(c)] = 1
    return np.abs(np.fft.rfft(M, axis=1)[:, 1:]) ** 2, M


Psc, Msc = spec_of([reps[i] for i in spec_code])
Psh, Msh = spec_of(ham)
Psy, Msy = spec_of(cyc)
print(f"\nUSE DEMO: K={K} words; Hamming code min Hamming distance {dmin_h}; spectral code min spectral distance {min(np.abs(a-b).sum() for a,b in itertools.combinations(Psc,2)):.2f}; "
      f"Hamming code min spectral distance {min(np.abs(a-b).sum() for a,b in itertools.combinations(Psh,2)):.3f}")
print("phase-blind receiver: observes |FFT(x + noise)|^2, decodes by nearest signature; symbol error rate (chance = %.2f)" % (1 - 1 / K))
print(f"{'sigma':>6} | {'spectral code':>13} {'Hamming code':>13} {'cyclic family':>14}")
res = {}
for sigma in (0.05, 0.1, 0.2, 0.3, 0.5):
    out = []
    for Pc, Mc in ((Psc, Msc), (Psh, Msh), (Psy, Msy)):
        err = 0; T = 2000
        for _ in range(T):
            k = int(rng.integers(len(Pc)))
            y = Mc[k] + sigma * rng.standard_normal(n)
            Py = np.abs(np.fft.rfft(y)[1:]) ** 2
            err += int(np.argmin(np.abs(Pc - Py).sum(1)) != k)
        out.append(err / T)
    res[sigma] = out
    print(f"{sigma:6.2f} | {out[0]:13.3f} {out[1]:13.3f} {out[2]:14.3f}", flush=True)
json.dump(dict(n=n, w=w, classes=N, table={str(k): v for k, v in table.items()}, ser={str(k): v for k, v in res.items()}), open(f"research/spectral_codes_{n}_{w}.json", "w"), indent=1, default=float)
