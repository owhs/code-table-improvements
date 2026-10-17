"""Check 1 of code_893.txt: a 5-ary code of length 9 with 893 words and minimum distance >= 5, by direct comparison of all pairs.
Needs numpy.  python check_code.py"""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
W = np.array([[int(c) for c in ln.strip()] for ln in open(os.path.join(here, "code_893.txt")) if ln.strip()], dtype=np.int16)
print("words", len(W), "length", W.shape[1], "symbols", sorted(set(W.ravel().tolist())))
print("distinct", len({tuple(w) for w in W.tolist()}) == len(W))
D = (W[:, None, :] != W[None, :, :]).sum(2)
np.fill_diagonal(D, 99)
md = int(D.min())
print("minimum distance", md, "-> valid A_5(9,5) >= 893" if md >= 5 and len(W) == 893 else "-> NOT valid")
