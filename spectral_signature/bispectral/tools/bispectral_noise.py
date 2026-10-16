"""Does a third-order (triple-correlation) receiver beat a power-spectrum receiver in UNKNOWN COLOURED Gaussian noise?
Theory: for zero-mean Gaussian noise with autocovariance R (circulant, PSD S=DFT(R)):
    E|Y(k)|^2        = P_x(k) + S(k)                          (noise floor; needs S to correct)
    E[T_y(a,b)]      = T_x(a,b) + w*(R(a)+R(b)+R(b-a))        (noise bias, same for every codeword of weight w)
so both signatures have a codeword-independent additive bias that depends on the unknown noise; the fluctuations differ.
Modes: 'blind' (no noise knowledge, uses raw signature), 'oracle' (true bias subtracted).  M snapshots averaged (signatures are shift invariant, so each snapshot may have a different shift).
Codes: K=12 by randomised farthest-point search in each receiver's own distance (P: dihedral/homometric classes; Tm: triple correlation merged with mirror).
    python tools/bispectral_noise.py n w [trials]"""
import sys, itertools, json, time
import numpy as np

n = int(sys.argv[1]); w = int(sys.argv[2]); TR = int(sys.argv[3]) if len(sys.argv) > 3 else 400
rng = np.random.default_rng(11)
K = 12
S_ = np.array(list(itertools.combinations(range(n), w)), dtype=np.int16)
X = np.zeros((len(S_), n), dtype=np.int8); np.put_along_axis(X, S_, 1, axis=1)
pw = (1 << np.arange(n, dtype=np.int64))
cod = np.stack([(np.roll(X, k, axis=1).astype(np.int64) * pw).sum(1) for k in range(n)]).min(0)
_, first = np.unique(cod, return_index=True)
X = X[np.sort(first)].astype(np.float64)
C = len(X)
ar = np.arange(n)


def triple_batch(Z):
    out = np.empty((len(Z), n, n))
    for a in range(n):
        xa = np.roll(Z, -a, axis=1)
        for b in range(n):
            out[:, a, b] = (Z * xa * np.roll(Z, -b, axis=1)).sum(1)
    return out


def triple1(y):
    R = np.stack([np.roll(y, -a) for a in range(n)])        # R[a,t]=y[t+a]
    return np.einsum('t,at,bt->ab', y, R, R)


T = triple_batch(X); Tref = T[:, (-ar) % n][:, :, (-ar) % n]
Tf = T.reshape(C, -1) / w ** 3; Tr = Tref.reshape(C, -1) / w ** 3
P = (np.abs(np.fft.rfft(X, axis=1)[:, 1:n // 2 + 1]) ** 2) / w ** 2
auto = np.round(np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(1, n // 2 + 1)], axis=1)).astype(int)
_, iP = np.unique(auto, axis=0, return_index=True); iP = np.sort(iP)
less = np.array([Tf[i].tolist() <= Tr[i].tolist() for i in range(C)])
Km = np.where(less[:, None], Tf, Tr)
_, iTm = np.unique(np.round(Km, 9), axis=0, return_index=True); iTm = np.sort(iTm)
print(f"n={n} w={w}: classes P {len(iP)}, Tm {len(iTm)}", flush=True)


def dP(idx):
    p = P[idx]; return np.abs(p[:, None, :] - p[None, :, :]).sum(2)


def dTm(idx):
    t = Tf[idx]; r = Tr[idx]
    D = np.empty((len(idx), len(idx)))
    for i in range(len(idx)):
        D[i] = np.minimum(np.abs(t - t[i]).sum(1), np.abs(r - t[i]).sum(1))
    return D


def farthest(D, K, tries=60):
    best = None; bm = -1
    for _ in range(tries):
        s = [int(rng.integers(len(D)))]
        while len(s) < K:
            s.append(int(np.argmax(D[:, s].min(1))))
        m = D[np.ix_(s, s)][np.triu_indices(K, 1)].min()
        if m > bm:
            bm, best = m, s
    return best, bm


sP, mP = farthest(dP(iP), K); idP = iP[sP]
sT, mT = farthest(dTm(iTm), K); idT = iTm[sT]
print(f"code min distances: P {mP:.3f}  Tm {mT:.3f}", flush=True)


def coloured_psd():
    k0 = rng.uniform(1.0, n / 4)
    al = rng.uniform(0.5, 2.5)
    k = np.minimum(ar, n - ar).astype(float)
    S = 1.0 / (1.0 + (k / k0) ** 2) ** al
    return S / S.mean()                     # mean PSD = 1 -> per-sample variance sigma^2 after scaling


def noise(S, sigma):
    z = np.fft.fft(rng.standard_normal(n))
    return np.real(np.fft.ifft(z * np.sqrt(S))) * sigma


def run(sigma, M, mode, coloured=True):
    errs = {"P": 0, "Tm": 0}
    for _ in range(TR):
        j = int(rng.integers(K))
        S = coloured_psd() if coloured else np.ones(n)
        Rn = np.real(np.fft.ifft(S)) * sigma ** 2             # autocovariance R(tau)
        Pm = np.zeros(n // 2); Tm_ = np.zeros((n, n))
        for _m in range(M):
            xs = np.roll(X[idP[j]], int(rng.integers(n)))
            y = xs + noise(S, sigma)
            Pm += (np.abs(np.fft.rfft(y)[1:n // 2 + 1]) ** 2) / w ** 2
        Pm /= M
        for _m in range(M):
            xs = np.roll(X[idT[j]], int(rng.integers(n)))
            y = xs + noise(S, sigma)
            Tm_ += triple1(y) / w ** 3
        Tm_ /= M
        if mode == "oracle":
            Sk = np.real(np.fft.fft(Rn))                      # PSD of the noise = DFT of R
            Pm = Pm - Sk[1:n // 2 + 1] / w ** 2
            bias = np.add.outer(Rn, Rn) + Rn[(ar[None, :] - ar[:, None]) % n]
            Tm_ = Tm_ - w * bias / w ** 3
        dpp = np.abs(P[idP] - Pm).sum(1)
        errs["P"] += int(np.argmin(dpp) != j)
        tt = Tm_.ravel()
        dtt = np.minimum(np.abs(Tf[idT] - tt).sum(1), np.abs(Tr[idT] - tt).sum(1))
        errs["Tm"] += int(np.argmin(dtt) != j)
    return errs["P"] / TR, errs["Tm"] / TR


res = {}
print(f"\nsymbol error rate, K={K}, {TR} trials per cell (chance {1-1/K:.2f}); noise = random coloured PSD (unknown) of per-sample std sigma\n")
print(f"{'noise':>9} {'mode':>7} {'M':>3} {'sigma':>6} | {'P':>7} {'Tm':>7}", flush=True)
for col in (False, True):
    for mode in ("blind", "oracle"):
        for M in (1, 8):
            for sigma in (0.3, 0.6, 1.0):
                a, b = run(sigma, M, mode, col)
                res[f"{'col' if col else 'white'}|{mode}|{M}|{sigma}"] = (a, b)
                print(f"{'coloured' if col else 'white':>9} {mode:>7} {M:3d} {sigma:6.2f} | {a:7.3f} {b:7.3f}", flush=True)
json.dump(dict(n=n, w=w, K=K, trials=TR, classes=[len(iP), len(iTm)], codemin=[mP, mT], res=res), open(f"research/bispnoise_{n}_{w}.json", "w"), indent=1)
