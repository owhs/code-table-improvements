"""Dependency-free verifier (standard library only) for the A_5(10,4) >= 18750 code in this folder.
    python check_code.py            # structural proof + the file equals the structure   (about 10 s)
    python check_code.py --brute    # additionally every pair of code_18750.txt (175 million pairs; slow in pure Python, uses numpy if available)

Structure. Let D = the 30 vectors of D_30.txt (elements of Z_5^5). The code is
    C = { (u | v) : u, v in Z_5^5,  v - u in D,  sum(u) + sum(v) = 1 (mod 5) }        (30 * 625 = 18750 words, coordinates 0-4 = u, 5-9 = v).
Equivalently C is a union of 30 cosets of the linear [10,4] code K' = { (a | a) : a in Z_5^5, sum(a) = 0 }.
Distance proof: take (u|u+d1), (u'|u'+d2) in C, w = u - u', e = d1 - d2. Their distance is wt(w) + wt(w + e) and the layer rule gives
    2 sum(w) + sum(e) = 0, i.e. sum(w) = -sum(e)/2 = 2 sum(e) (mod 5)    [w = every vector of that sum class occurs].
For d1 = d2 the same formula holds with w != 0 (and sum(w) = 0). So the minimum distance is the minimum over d1, d2 in D of
    min { wt(w) + wt(w + e) : sum(w) = 2 sum(e), and w != 0 if e = 0 }.
That is computed below for all 30 * 31 / 2 pairs (about 1.45 million evaluations)."""
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
q, m, d = 5, 5, 4
D = [tuple(int(c) for c in ln.strip()) for ln in open(os.path.join(HERE, "D_30.txt")) if ln.strip()]
assert len(D) == 30 and len(set(D)) == 30
wt = lambda v: sum(1 for x in v if x)
vecs = list(itertools.product(range(q), repeat=m))
by_sum = {s: [w for w in vecs if sum(w) % q == s] for s in range(q)}
md = 99
for i in range(len(D)):
    for j in range(i, len(D)):
        e = tuple((a - b) % q for a, b in zip(D[i], D[j]))
        need = (2 * sum(e)) % q
        for w in by_sum[need]:
            if i == j and not any(w):
                continue
            val = wt(w) + wt(tuple((x + y) % q for x, y in zip(w, e)))
            if val < md:
                md = val
print("structural minimum distance:", md, "(need >=", d, ")")
words = set()
for dl in D:
    for u in vecs:
        if (2 * sum(u) + sum(dl)) % q == 1:
            words.add(u + tuple((a + b) % q for a, b in zip(u, dl)))
print("structure words:", len(words))
fileset = {tuple(int(c) for c in ln.strip()) for ln in open(os.path.join(HERE, "code_18750.txt")) if ln.strip()}
print(f"file words: {len(fileset)}; file == structure: {fileset == words}")
ok = md >= d and fileset == words and len(words) == 18750
if "--brute" in sys.argv:
    try:
        import numpy as np
        A = np.array(sorted(fileset), dtype=np.int8)
        mdb = A.shape[1] + 1
        for r0 in range(0, len(A), 150):
            blk = A[r0:r0 + 150]
            dist = np.zeros((len(blk), len(A)), dtype=np.int16)
            for c in range(A.shape[1]):
                dist += blk[:, c][:, None] != A[:, c][None, :]
            for i in range(len(blk)):
                dist[i, r0 + i] = 99
            mdb = min(mdb, int(dist.min()))
        print("brute-force minimum distance (numpy):", mdb)
        ok = ok and mdb >= d
    except ImportError:
        print("numpy not available: brute force skipped")
print("RESULT:", "VALID A_5(10,4) >= 18750" if ok else "FAILED")
sys.exit(0 if ok else 1)
