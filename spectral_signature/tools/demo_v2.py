"""Fair demo: phase-blind receiver (|FFT|^2 nearest-signature decoding of word + Gaussian noise).  Codebooks of K words:
  spectral-opt : K classes chosen to maximise the minimum spectral distance (best of 300 randomised farthest-point selections)
  distinct     : K random classes (distinct signatures, no distance optimisation)           <- separates 'distinct' from 'well separated'
  Hamming x5   : greedy codes of maximal Hamming distance from 5 different random draws (mean, plus min/max over draws)
  cyclic       : K cyclic shifts of one word
3000 trials per cell. Writes _demo_v2.json in the publish folder."""
import itertools, json, os
import numpy as np

OUT = r"C:\projects\2026\reframe\publish\Spectral_codes_2026"
rng = np.random.default_rng(21)


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
    for i in range(len(reps)):
        b = Wb[i].astype(int)
        keys.setdefault(tuple(int(sum(b[j] * b[(j + k) % n] for j in range(n))) for k in range(1, n // 2 + 1)), i)
    first = sorted(keys.values())
    return [reps[i] for i in first]


def sigs(code, n):
    M = np.zeros((len(code), n))
    for i, c in enumerate(code):
        M[i, list(c)] = 1
    return np.abs(np.fft.rfft(M, axis=1)[:, 1:]) ** 2, M


def mind(P, w):
    return min(np.abs(a - b).sum() for a, b in itertools.combinations(P, 2)) / w ** 2


def ser(P, M, n, sigma, T=3000):
    err = 0
    for _ in range(T):
        k = int(rng.integers(len(P)))
        y = M[k] + sigma * rng.standard_normal(n)
        err += int(np.argmin(np.abs(P - np.abs(np.fft.rfft(y)[1:]) ** 2).sum(1)) != k)
    return err / T


def farthest(Dm, K, tries=300):
    N = len(Dm)
    best, bd = None, -1
    for _ in range(tries):
        sel = [int(rng.integers(N))]
        md = Dm[sel[0]].copy()
        while len(sel) < K:
            v = int(np.argmax(md + 1e-9 * rng.random(N)))
            sel.append(v); md = np.minimum(md, Dm[v])
        d = min(Dm[i, j] for i, j in itertools.combinations(sel, 2))
        if d > bd:
            best, bd = sel, d
    return best, bd


SIG = (0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
res = []
for (n, w, K) in ((16, 6, 12), (13, 4, 13), (18, 6, 16), (14, 5, 20)):
    reps = classes(n, w)
    Pc, _ = sigs(reps, n)
    Dm = np.abs(Pc[:, None, :] - Pc[None, :, :]).sum(2) / w ** 2
    sel, bd = farthest(Dm, K)
    spec = [reps[i] for i in sel]
    dist = [reps[i] for i in rng.choice(len(reps), K, replace=False)]
    allw = list(itertools.combinations(range(n), w))
    hams = []
    for seed in range(5):
        for dmin in range(2 * w, 1, -2):
            h = []
            for i in rng.permutation(len(allw)):
                c = allw[i]
                if all(len(set(c) ^ set(x)) >= dmin for x in h):
                    h.append(c)
                if len(h) >= K:
                    break
            if len(h) >= K:
                break
        hams.append(h[:K])
    base = spec[0]
    cyc = [tuple(sorted((i + r) % n for i in base)) for r in range(min(K, n))]
    sp, ms = sigs(spec, n); dp, md_ = sigs(dist, n); cp, mc = sigs(cyc, n)
    hs = [sigs(h, n) for h in hams]
    info = dict(n=n, w=w, K=K, classes=len(reps), chance=1 - 1 / K,
                spectral_min_d=float(bd), distinct_min_d=float(mind(dp, w)),
                hamming_min_d_h=[int(min(len(set(a) ^ set(b)) for a, b in itertools.combinations(h, 2))) for h in hams],
                hamming_min_spectral=[float(mind(p, w)) for p, _ in hs], rows=[])
    for sg in SIG:
        r = dict(sigma=sg, spectral_opt=ser(sp, ms, n, sg), distinct=ser(dp, md_, n, sg), cyclic=ser(cp, mc, n, sg))
        hv = [ser(p, m, n, sg, 1500) for p, m in hs]
        r.update(hamming_mean=float(np.mean(hv)), hamming_min=float(min(hv)), hamming_max=float(max(hv)))
        info["rows"].append(r)
    res.append(info)
    print(f"n={n} w={w} K={K}: spectral-opt min d {bd:.3f}; random-distinct min d {info['distinct_min_d']:.3f}; Hamming draws min spectral d {np.round(info['hamming_min_spectral'],3)}", flush=True)
    for r in info["rows"]:
        print(f"   sigma {r['sigma']:.2f}: spectral-opt {r['spectral_opt']:.3f} | distinct {r['distinct']:.3f} | Hamming mean {r['hamming_mean']:.3f} [{r['hamming_min']:.3f}-{r['hamming_max']:.3f}] | cyclic {r['cyclic']:.3f}", flush=True)
json.dump(res, open(os.path.join(OUT, "_demo_v2.json"), "w"), indent=1)
print("demo v2 written")
