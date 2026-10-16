"""Independent plain-Python verifier for q-ary / mixed binary-ternary codes (no numpy, no numba, no project imports).
    python -I tools/qary_verify.py code.json [more.json ...]
JSON keys: "qs" (alphabet size per coordinate) or ("m_bin", "n_ter") meaning qs = [2]*m_bin + [3]*n_ter;
"d" (required minimum distance); "words" (list of digit lists). Checks lengths, symbol ranges, distinctness and every
pairwise Hamming distance by direct symbol comparison. Prints size and exact minimum distance; exit code 1 on failure."""
import json
import sys
from itertools import combinations
from operator import ne


def verify(path):
    c = json.load(open(path))
    if "qs" in c:
        qs = [int(q) for q in c["qs"]]
    else:
        qs = [2] * int(c["m_bin"]) + [3] * int(c["n_ter"])
    d = int(c["d"])
    words = [tuple(int(v) for v in w) for w in c["words"]]
    n = len(qs)
    for w in words:
        if len(w) != n:
            return False, f"word of length {len(w)} != {n}"
        for v, q in zip(w, qs):
            if not 0 <= v < q:
                return False, f"symbol {v} out of range for alphabet {q}"
    if len(set(words)) != len(words):
        return False, "duplicate words"
    mind = n + 1
    for a, b in combinations(words, 2):
        dd = sum(map(ne, a, b))
        if dd < mind:
            mind = dd
            if mind < d:
                return False, f"distance {dd} < {d} between {a} and {b}"
    claimed = c.get("size")
    if claimed is not None and int(claimed) != len(words):
        return False, f"claimed size {claimed} != {len(words)} words"
    return True, f"OK: {len(words)} words, length {n} (qs: {qs.count(2)} binary + {qs.count(3)} ternary), " \
                 f"min distance {mind if len(words) > 1 else 'n/a'} >= {d}"


def main():
    bad = 0
    for p in sys.argv[1:]:
        ok, msg = verify(p)
        print(f"{p}: {msg}")
        bad += not ok
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
