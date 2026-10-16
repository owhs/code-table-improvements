"""Bispectral-signature codes vs power-spectrum codes (constant-weight binary words of length n, weight w).
Three receivers, each with its own invariances:
  P   : sees |DFT|^2 (k=1..n/2).  Invariant to shift AND reflection AND homometry.  d_P(x,y) = sum_k |P_x(k)-P_y(k)| / w^2     (Spectral_codes_2026 convention)
  T   : sees the triple correlation T(a,b)=sum_t x_t x_{t+a} x_{t+b} (= bispectrum info).  Invariant to shift only (splits mirror pairs).  d_T = sum_{a,b}|T_x-T_y| / w^3
  Tm  : T receiver that is also blind to reflection (signature of x merged with its mirror): d_Tm = min(d_T(x,y), d_T(x,mirror y)).
Capacity S(theta) = max set of classes with pairwise d >= theta (exact MILP, HiGHS, time-limited, bound reported).
Then a noise demo: K-subset by randomised farthest-point search, nearest-signature decoder, empirical signature of y = x + sigma*noise.
    python tools/bispectral_codes.py n w"""
import sys, itertools, json, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

n = int(sys.argv[1]); w = int(sys.argv[2])
THETAS = (0.05, 0.10, 0.20, 0.40)
t0 = time.time()

# translation classes
S = np.array(list(itertools.combinations(range(n), w)), dtype=np.int16)
X = np.zeros((len(S), n), dtype=np.int8); np.put_along_axis(X, S, 1, axis=1)
pw = (1 << np.arange(n, dtype=np.int64))
codes = np.stack([(np.roll(X, k, axis=1).astype(np.int64) * pw).sum(1) for k in range(n)]).min(0)
_, first = np.unique(codes, return_index=True)
X = X[np.sort(first)].astype(np.int64)
C = len(X)


def triple(Z):
    out = np.empty((len(Z), n, n), dtype=np.int16)
    for a in range(n):
        xa = np.roll(Z, -a, axis=1)
        for b in range(n):
            out[:, a, b] = (Z * xa * np.roll(Z, -b, axis=1)).sum(1)
    return out


T = triple(X)                                              # (C,n,n)
Tref = T[:, (-np.arange(n)) % n][:, :, (-np.arange(n)) % n]   # T of the mirror image: T(-a,-b)
P = np.abs(np.fft.rfft(X, axis=1)[:, 1:n // 2 + 1]) ** 2
auto = np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(1, n // 2 + 1)], axis=1)
print(f"n={n} w={w}: {C} translation classes  [{time.time()-t0:.1f}s]", flush=True)


def reps_of(keys):
    _, idx = np.unique(keys, axis=0, return_index=True)
    return np.sort(idx)


iP = reps_of(auto)                                          # P classes (shift+mirror+homometric merged)
Tf = T.reshape(C, -1); Tr = Tref.reshape(C, -1)
# mirror-merged class key: lexicographically smaller of (T, Tref)
less = np.zeros(C, dtype=bool)
for i in range(C):
    a, b = Tf[i].tolist(), Tr[i].tolist()
    less[i] = a <= b
Km = np.where(less[:, None], Tf, Tr)
iTm = reps_of(Km)
iT = np.arange(C)
print(f"classes: P (power spectrum) {len(iP)} | Tm (triple, mirror-blind) {len(iTm)} | T (triple, shift-only) {len(iT)}  [{time.time()-t0:.1f}s]", flush=True)


def dist_P(idx):
    p = P[idx]
    return np.abs(p[:, None, :] - p[None, :, :]).sum(2) / w ** 2


def dist_T(idx, mirror=False):
    t = Tf[idx].astype(np.int32); r = Tr[idx].astype(np.int32)
    D = np.empty((len(idx), len(idx)), dtype=np.float64)
    for i in range(len(idx)):
        d1 = np.abs(t - t[i]).sum(1)
        if mirror:
            d1 = np.minimum(d1, np.abs(r - t[i]).sum(1))
        D[i] = d1
    return D / w ** 3


def max_code(D, theta, limit=60):
    N = len(D)
    bad = np.argwhere(np.triu(D < theta - 1e-12, 1))
    if len(bad) == 0:
        return list(range(N)), True, N
    rows = np.repeat(np.arange(len(bad)), 2)
    A = coo_matrix((np.ones(2 * len(bad)), (rows, bad.ravel())), shape=(len(bad), N)).tocsr()
    r = milp(-np.ones(N), constraints=[LinearConstraint(A, -np.inf, 1)], integrality=np.ones(N), bounds=Bounds(0, 1), options=dict(time_limit=limit, disp=False))
    sel = sorted(np.where(r.x > 0.5)[0].tolist()) if r.x is not None else []
    ub = int(np.floor(-r.mip_dual_bound + 1e-6)) if getattr(r, "mip_dual_bound", None) is not None else N
    return sel, r.status == 0, ub


rec = {}
Ds = {"P": (dist_P(iP), iP), "Tm": (dist_T(iTm, True), iTm), "T": (dist_T(iT, False), iT)}
print(f"\n{'theta':>6} | " + " | ".join(f"{k:>22}" for k in Ds), flush=True)
for th in THETAS:
    cells = []
    for k, (D, idx) in Ds.items():
        sel, opt, ub = max_code(D, th)
        cells.append(f"{len(sel):5d} {'opt' if opt else f'<= {ub}':>6} (of {len(idx)})")
        rec.setdefault(k, {})[th] = dict(size=len(sel), optimal=bool(opt), upper=ub, classes=len(idx))
    print(f"{th:6.2f} | " + " | ".join(f"{c:>22}" for c in cells), flush=True)
for k, (D, idx) in Ds.items():
    iu = np.triu_indices(len(D), 1)
    print(f"  {k}: distance median {np.median(D[iu]):.3f}, max {D[iu].max():.3f}")

# ---- noise demo: K words by farthest-point search, nearest-signature decoder on empirical signatures
rng = np.random.default_rng(7)
K = 12


def farthest(D, K, tries=40):
    best = None; bm = -1
    for _ in range(tries):
        s = [int(rng.integers(len(D)))]
        while len(s) < K:
            s.append(int(np.argmax(D[:, s].min(1))))
        m = D[np.ix_(s, s)][np.triu_indices(K, 1)].min()
        if m > bm:
            bm, best = m, s
    return best, bm


def sig(kind, y):
    if kind == "P":
        return np.abs(np.fft.rfft(y)[1:n // 2 + 1]) ** 2 / w ** 2
    t = np.array([[(y * np.roll(y, -a) * np.roll(y, -b)).sum() for b in range(n)] for a in range(n)]).ravel()
    return t / w ** 3


def sigclean(kind, i):
    return (P[i] / w ** 2) if kind == "P" else (Tf[i] / w ** 3)


print(f"\nNOISE DEMO (K={K}, farthest-point codes, nearest-signature decoding on the empirical signature of x + sigma*noise; chance {1-1/K:.2f}):")
print(f"{'sigma':>6} | {'P':>8} {'Tm':>8} {'T':>8}   (code min distance: " + ", ".join(f"{k} {farthest(D, K)[1]:.3f}" for k, (D, _) in Ds.items()) + ")", flush=True)
codes_ = {}
for k, (D, idx) in Ds.items():
    s, m = farthest(D, K); codes_[k] = idx[s]
ser = {}
for sigma in (0.05, 0.1, 0.2, 0.3):
    row = []
    for k in ("P", "Tm", "T"):
        ids = codes_[k]
        if k == "Tm":
            # mirror-blind decoding: compare against signature and mirrored signature
            refs = np.stack([Tf[i] / w ** 3 for i in ids]); refm = np.stack([Tr[i] / w ** 3 for i in ids])
        elif k == "T":
            refs = np.stack([Tf[i] / w ** 3 for i in ids]); refm = None
        else:
            refs = np.stack([P[i] / w ** 2 for i in ids]); refm = None
        err = 0; trials = 400
        for _ in range(trials):
            j = int(rng.integers(K))
            y = X[ids[j]] + sigma * rng.standard_normal(n)
            s_ = sig("P" if k == "P" else "T", y)
            d = np.abs(refs - s_).sum(1)
            if refm is not None:
                d = np.minimum(d, np.abs(refm - s_).sum(1))
            err += int(np.argmin(d) != j)
        row.append(err / trials)
    ser[sigma] = row
    print(f"{sigma:6.2f} | {row[0]:8.3f} {row[1]:8.3f} {row[2]:8.3f}", flush=True)
json.dump(dict(n=n, w=w, translation_classes=C, classes={k: len(v[1]) for k, v in Ds.items()}, capacity=rec, ser={str(k): v for k, v in ser.items()}),
          open(f"research/bispectral_{n}_{w}.json", "w"), indent=1, default=float)
print(f"done [{time.time()-t0:.1f}s]")
