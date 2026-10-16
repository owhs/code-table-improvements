"""Dependency-free checker for the spectral-signature codes in this folder.
    python check_codes.py            # verify every codes/*.txt (weight, length, distinctness, signature separation) and compare with codes.json
    python check_codes.py --recount  # additionally recount c(n,w) for n <= 14 by brute force (exact integer autocorrelation keys) against table_classes.csv
Only the standard library is used. The signature of a word with support S is P(k) = sum_{a,b in S} cos(2 pi k (a-b)/n), k = 1..floor(n/2) (= |DFT|^2);
the distance is sum_k |P_x(k) - P_y(k)| / w^2. A code passes if all pairwise distances are >= theta - 1e-9 and all integer autocorrelations are pairwise distinct."""
import csv
import itertools
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def signature(S, n):
    return [math.fsum(math.cos(2 * math.pi * k * (a - b) / n) for a in S for b in S) for k in range(1, n // 2 + 1)]


def autocorr(S, n):
    s = set(S)
    return tuple(sum(1 for a in S if (a + k) % n in s) for k in range(1, n // 2 + 1))


def check_file(path, claimed=None):
    m = re.fullmatch(r"n(\d+)_w(\d+)_d([\d.]+)\.txt", os.path.basename(path))
    n, w, theta = int(m.group(1)), int(m.group(2)), float(m.group(3))
    rows = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    supports = []
    for r in rows:
        if len(r) != n or set(r) - {"0", "1"}:
            return False, f"bad row {r!r}"
        S = [i for i, c in enumerate(r) if c == "1"]
        if len(S) != w:
            return False, f"weight {len(S)} != {w}"
        supports.append(S)
    if len(set(map(tuple, supports))) != len(supports):
        return False, "duplicate words"
    sigs = [signature(S, n) for S in supports]
    acs = [autocorr(S, n) for S in supports]
    if len(set(acs)) != len(acs):
        return False, "two words share the same autocorrelation"
    md = min((sum(abs(a - b) for a, b in zip(sigs[i], sigs[j])) / w ** 2 for i, j in itertools.combinations(range(len(rows)), 2)), default=float("inf"))
    ok = md >= theta - 1e-9
    if claimed is not None and claimed != len(rows):
        return False, f"size {len(rows)} != claimed {claimed}"
    return ok, f"size {len(rows)}, min spectral distance {md:.4f} (theta {theta}, margin {md - theta:+.4f})"


def main():
    cj = json.load(open(os.path.join(HERE, "codes.json")))
    bad = 0
    files = sorted(os.listdir(os.path.join(HERE, "codes")))
    for fn in files:
        key = fn[:-4]
        ok, msg = check_file(os.path.join(HERE, "codes", fn), cj[key]["size"] if key in cj else None)
        if not ok:
            bad += 1
        print(("OK   " if ok else "FAIL ") + fn, msg)
    print(f"\n{len(files)} files, {bad} failures")
    if "--recount" in sys.argv:
        tab = {(int(r["n"]), int(r["w"])): int(r["classes"]) for r in csv.DictReader(open(os.path.join(HERE, "table_classes.csv")))}
        mism = 0
        for n in range(4, 15):
            for w in range(2, n // 2 + 1):
                keys = {autocorr([i for i in range(n) if c >> i & 1], n) for c in (sum(1 << i for i in comb) for comb in itertools.combinations(range(n), w))}
                if len(keys) != tab[(n, w)]:
                    mism += 1
                    print(f"RECOUNT MISMATCH n={n} w={w}: {len(keys)} vs table {tab[(n, w)]}")
        print(f"recount n<=14: {mism} mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
