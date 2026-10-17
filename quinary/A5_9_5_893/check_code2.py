"""Check 2 of code_893.txt, a different method: agreements between all pairs as a one-hot matrix product (distance = length - agreements).
Needs numpy.  python check_code2.py"""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
rows = [ln.strip() for ln in open(os.path.join(here, "code_893.txt")) if ln.strip()]
n = len(rows[0])
assert all(len(r) == n and set(r) <= set("01234") for r in rows)
H = np.zeros((len(rows), 5 * n), dtype=np.float32)
for i, r in enumerate(rows):
    for k, ch in enumerate(r):
        H[i, 5 * k + int(ch)] = 1.0
A = H @ H.T  # A[i, j] = number of positions where words i and j agree
np.fill_diagonal(A, -1)
agree_max = int(A.max())
print("words", len(rows), "length", n, "largest agreement between two different words", agree_max, "-> minimum distance", n - agree_max)
print("valid A_5(9,5) >= 893" if n - agree_max >= 5 and len(set(rows)) == len(rows) and len(rows) == 893 else "NOT valid")
