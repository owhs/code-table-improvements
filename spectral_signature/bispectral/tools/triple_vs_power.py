"""Phase-blind sensor (power spectrum = autocorrelation) vs one that sees third-order structure (triple correlation = bispectrum), on w-subsets of Z_n.
Classes: translation orbits; translation+reflection orbits; autocorrelation (|DFT|^2) classes; triple-correlation classes.
Homometric (spectrum-blind) pairs that are NOT mirror images are the genuinely hidden ones; the question is whether third-order structure resolves them."""
import itertools, sys
import numpy as np

def run(n):
    out = []
    for w in range(2, n // 2 + 1):
        S = np.array(list(itertools.combinations(range(n), w)), dtype=np.int16)
        M = len(S)
        X = np.zeros((M, n), dtype=np.int8)
        np.put_along_axis(X, S, 1, axis=1)
        # translation canonical form: lexicographically min rotation, as bytes
        def canon(Y):
            best = None
            keys = []
            R = [np.roll(Y, k, axis=1) for k in range(n)]
            R = np.stack(R)  # n,M,n
            # encode as int via dot with powers (n<=22 fits int64)
            pw = (1 << np.arange(n, dtype=np.int64))
            codes = (R.astype(np.int64) * pw).sum(2)  # n,M
            return codes.min(0)
        ct = canon(X)
        cm = np.minimum(ct, canon(X[:, ::-1].copy()))
        # autocorrelation
        A = np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(n)], axis=1).astype(np.int16)
        # triple correlation T[a,b]=sum_t x(t)x(t+a)x(t+b), a<=b suffices up to symmetry; use all a,b
        T = np.empty((M, n * n), dtype=np.int8)
        for a in range(n):
            xa = np.roll(X, -a, axis=1)
            for b in range(n):
                T[:, a * n + b] = (X * xa * np.roll(X, -b, axis=1)).sum(1)
        def ncls(Z):
            Zc = np.ascontiguousarray(Z)
            v = Zc.view(np.dtype((np.void, Zc.dtype.itemsize * Zc.shape[1])))
            return len(np.unique(v))
        nt, nm = len(np.unique(ct)), len(np.unique(cm))
        na, n3 = ncls(A), ncls(T)
        out.append((n, w, M, nt, nm, na, n3))
        print(f"n={n:2d} w={w:2d}  subsets {M:7d} | translation classes {nt:6d} | +mirror {nm:6d} | autocorr (power-spectrum) classes {na:6d} | triple-corr classes {n3:6d}"
              f" | hidden from spectrum (beyond mirror): {nm-na:5d} | resolved by triple: {'ALL' if n3==nt else f'{n3} of {nt}'}", flush=True)
    return out

if __name__ == "__main__":
    for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        run(n)
