"""Dependency-free verifier (standard library only) for the A_5(9,4) >= 4500 code in this folder.
    python check_code.py          # rebuilds the code from D_36.txt, compares it with code_4500.txt, then checks every pair (10 million pairs, about 1-2 minutes)
Structure. Write a word as (u | v | z) with u, v in Z_5^4 and z in Z_5 (coordinates 0-3 = u, 4-7 = v, 8 = z). D is the set of 36 pairs (delta, z) in D_36.txt (delta in Z_5^4, z in Z_5). The code is
    C = { (u | u + delta | z) : (delta, z) in D,  sum(u) + sum(u + delta) = 1 (mod 5) }        36 * 125 = 4500 words
(each (delta, z) contributes the 125 vectors u with sum(2u + delta) = 1)."""
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
q = 5
D = []
for ln in open(os.path.join(HERE, "D_36.txt")):
    if ln.strip():
        d, z = ln.split()
        D.append((tuple(int(c) for c in d), tuple(int(c) for c in z)))
assert len(D) == 36
vecs = list(itertools.product(range(q), repeat=4))
words = set()
for d, z in D:
    for u in vecs:
        v = tuple((a + b) % q for a, b in zip(u, d))
        if (sum(u) + sum(v)) % q == 1:
            words.add(u + v + z)
print("words from the structure:", len(words))
fileset = {tuple(int(c) for c in ln.strip()) for ln in open(os.path.join(HERE, "code_4500.txt")) if ln.strip()}
print("file words:", len(fileset), "; file == structure:", fileset == words)
ok = fileset == words and len(words) == 4500
W = sorted(fileset)
md = 99
for i in range(len(W)):
    a = W[i]
    for j in range(i + 1, len(W)):
        b = W[j]
        dist = (a[0] != b[0]) + (a[1] != b[1]) + (a[2] != b[2]) + (a[3] != b[3]) + (a[4] != b[4]) + (a[5] != b[5]) + (a[6] != b[6]) + (a[7] != b[7]) + (a[8] != b[8])
        if dist < md:
            md = dist
            if md < 4:
                break
    if md < 4:
        break
print("minimum distance over all pairs:", md, "(need >= 4)")
ok = ok and md >= 4
print("RESULT:", "VALID A_5(9,4) >= 4500" if ok else "FAILED")
sys.exit(0 if ok else 1)
