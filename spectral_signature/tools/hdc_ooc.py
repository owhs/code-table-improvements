"""HDC-style bundling memory with sparse binary (constant-weight) codewords:  an item = (codeword m, role = cyclic shift s).  Memory = SUM of T shifted codewords
(integer vector of length n).  Readout (resonance): score(m,s) = <memory, shift_s(word_m)>; pick the T best pairs.  Compares codebooks:
  OOC      greedy set with all auto/cross-correlation <= 1 (an optical orthogonal code: every pair of shifted words overlaps in at most 1 position)
  random   random weight-w words
  Hamming  greedy set with large Hamming distance (the usual 'good constant-weight code', no correlation constraint)
Metric: fraction of the T stored items recovered exactly.  Also report how many (m,s) pairs there are (capacity of the address space)."""
import itertools, sys, json
import numpy as np

rng = np.random.default_rng(5)
n, w = int(sys.argv[1]) if len(sys.argv) > 1 else 127, int(sys.argv[2]) if len(sys.argv) > 2 else 4


def diffs(S):
    return {(a - b) % n for a in S for b in S if a != b}


def ooc_book(M, tries=200000):
    book, used = [], set()
    for _ in range(tries):
        S = tuple(sorted(rng.choice(n, w, replace=False).tolist()))
        d = [(a - b) % n for a in S for b in S if a != b]
        if len(set(d)) == len(d) and not (set(d) & used):
            book.append(S); used |= set(d)
            if len(book) == M:
                break
    return book


def rand_book(M):
    return [tuple(sorted(rng.choice(n, w, replace=False).tolist())) for _ in range(M)]


def ham_book(M, tries=20000):
    book = []
    best_d = 2 * w
    while len(book) < M and best_d >= 2:
        book = []
        for _ in range(tries):
            S = tuple(sorted(rng.choice(n, w, replace=False).tolist()))
            if all(len(set(S) ^ set(B)) >= best_d for B in book):
                book.append(S)
                if len(book) == M:
                    break
        if len(book) < M:
            best_d -= 2
    return book


def mat(book):
    W = np.zeros((len(book), n))
    for i, S in enumerate(book):
        W[i, list(S)] = 1
    return W


def trial(W, T):
    M = len(W)
    pairs = set()
    while len(pairs) < T:
        pairs.add((int(rng.integers(M)), int(rng.integers(n))))
    mem = np.zeros(n)
    for m, s in pairs:
        mem += np.roll(W[m], s)
    # resonance scores for all (m, s): circular cross-correlation of memory with each word
    F = np.fft.rfft(mem)
    sc = np.fft.irfft(F[None, :] * np.conj(np.fft.rfft(W, axis=1)), n=n, axis=1)   # sc[m, s] = sum_j mem[j + s] W[m, j]  (shift convention: roll by -s)
    # score of hypothesis shift s means word rolled by s: <mem, roll(W,s)> = sum_j mem[j] W[(j - s) % n]  -> corr at lag s
    sc = np.fft.irfft(F[None, :] * np.conj(np.fft.rfft(W, axis=1)), n=n, axis=1)
    top = np.argsort(-sc.ravel())[:T]
    got = {(int(i // n), int(i % n)) for i in top}
    return len(got & pairs) / T


M = 10
books = {"OOC (corr<=1)": ooc_book(M), "random": rand_book(M), "Hamming code": ham_book(M)}
print(f"n={n} w={w}; codebook sizes {[len(b) for b in books.values()]}; address space {M*n} (word, shift) pairs")
for name, b in books.items():
    ac = max(len(set(S) & {(s + 1) % n for s in S}) for S in b)  # placeholder crude value not used
    W = mat(b)
    maxcross = 0
    for i in range(len(b)):
        for j in range(len(b)):
            for s in range(n):
                if i == j and s == 0:
                    continue
                maxcross = max(maxcross, int((W[i] * np.roll(W[j], s)).sum()))
    print(f"  {name:14s} max auto/cross-correlation over all shifts: {maxcross}")
print(f"{'T items':>8} | " + " ".join(f"{k:>14s}" for k in books))
res = {}
for T in (2, 3, 4, 6, 8, 12, 16, 24):
    row = {}
    for k, b in books.items():
        W = mat(b)
        row[k] = float(np.mean([trial(W, T) for _ in range(400)]))
    res[T] = row
    print(f"{T:8d} | " + " ".join(f"{row[k]:14.3f}" for k in books), flush=True)
json.dump(res, open(f"research/hdc_ooc_{n}_{w}.json", "w"), indent=1)
