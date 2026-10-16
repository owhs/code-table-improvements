"""Word-level local search starting from a known record constant-weight code (Brouwer's code files, hex words after '$BASE=16').
For every weight-w word outside the code we count the code words it conflicts with (overlap > lam = w - d/2), stopping at 3.
  0 conflicts: the code is not maximal, add the word.
  1 conflict : candidate for a 1-swap; two mutually compatible candidates sharing the same single conflict give +1.
  2 conflicts: candidates for 2-swaps; three mutually compatible candidates whose conflicts lie in {c1,c2} give +1.
Repeats until no gain (max --rounds). Certified gains are written to runs/cwimprove/hit_<n>_<d>_<w>.json.
usage: python -I tools/cw_improve.py <code file> <n> <d> <w> [--rounds 5]
Needs n <= 32 (one scan costs C(n,w) * |code| popcounts in the worst case)."""
import argparse, json, os, sys, time
from itertools import combinations

import numpy as np
from numba import njit

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


@njit(cache=True)
def pc(x):
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    return (x * 0x0101010101010101) >> 56


@njit(cache=True)
def scan(n, w, lam, code, cap):
    ncode = code.shape[0]
    words = np.zeros(cap, dtype=np.int64)
    c1 = np.full(cap, -1, dtype=np.int64)
    c2 = np.full(cap, -1, dtype=np.int64)
    k = 0
    m = (np.int64(1) << w) - 1
    lim = np.int64(1) << n
    while m < lim:
        conf = 0
        a = -1
        b = -1
        for ci in range(ncode):
            if pc(m & code[ci]) > lam:
                conf += 1
                if conf == 1:
                    a = ci
                elif conf == 2:
                    b = ci
                else:
                    break
        if conf <= 2 and k < cap:
            if not (conf == 1 and code[a] == m):
                words[k] = m
                c1[k] = a
                c2[k] = b
                k += 1
        c = m & -m
        r = m + c
        m = (((r ^ m) >> 2) // c) | r
    return words[:k], c1[:k], c2[:k]


def load(path, w):
    ws = []
    for line in open(path):
        line = line.strip()
        if line and not line.startswith("$"):
            ws.append(int(line, 16))
    assert all(bin(x).count("1") == w for x in ws), "weights differ from w"
    return ws


def valid(code, lam):
    a = np.array(code, dtype=np.int64)
    for i in range(len(a) - 1):
        ov = np.array([bin(int(x)).count("1") for x in (a[i + 1:] & a[i])])
        if len(ov) and ov.max() > lam:
            return False
    return True


def comp(x, y, lam):
    return bin(x & y).count("1") <= lam


def independent(cands, need, lam):
    """is there a set of `need` mutually compatible words (small backtracking)?"""
    def rec(chosen, rest):
        if len(chosen) == need:
            return chosen
        for i, x in enumerate(rest):
            if all(comp(x, y, lam) for y in chosen):
                r = rec(chosen + [x], rest[i + 1:])
                if r:
                    return r
        return None
    return rec([], cands)


def improve_once(code, n, w, lam, cap=6_000_000):
    arr = np.array(code, dtype=np.int64)
    words, c1, c2 = scan(n, w, lam, arr, cap)
    z = [int(x) for x, a in zip(words, c1) if a < 0]
    stats = {"free": len(z), "one": int((c1 >= 0).sum() - (c2 >= 0).sum()), "two": int((c2 >= 0).sum()), "capped": len(words) >= cap}
    if z:  # not maximal: greedy add
        added = []
        for x in z:
            if all(comp(x, y, lam) for y in added):
                added.append(x)
        return code + added, stats
    by1, by2 = {}, {}
    for x, a, b in zip(words, c1, c2):
        if b < 0:
            by1.setdefault(int(a), []).append(int(x))
        else:
            by2.setdefault((int(a), int(b)), []).append(int(x))
    for a, L in by1.items():
        r = independent(L[:400], 2, lam)
        if r:
            return [x for i, x in enumerate(code) if i != a] + r, stats
    for (a, b), L in by2.items():
        L = L + by1.get(a, []) + by1.get(b, [])
        if len(L) >= 3:
            r = independent(L[:300], 3, lam)
            if r:
                return [x for i, x in enumerate(code) if i not in (a, b)] + r, stats
    return None, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file"); ap.add_argument("n", type=int); ap.add_argument("d", type=int); ap.add_argument("w", type=int)
    ap.add_argument("--rounds", type=int, default=5)
    a = ap.parse_args()
    lam = a.w - a.d // 2
    code = load(a.file, a.w)
    rec = len(code)
    print(f"A({a.n},{a.d},{a.w}): loaded {rec} words, valid={valid(code, lam)}", flush=True)
    for rd in range(a.rounds):
        t = time.time()
        new, st = improve_once(code, a.n, a.w, lam)
        print(f"  round {rd}: {st} {time.time() - t:.0f}s -> {'IMPROVED to ' + str(len(new)) if new else 'no gain'}", flush=True)
        if not new:
            break
        code = new
    if len(code) > rec:
        from engine.problems.cwbase import make
        words = [[b for b in range(a.n) if (int(m) >> b) & 1] for m in code]
        cert = make(a.n, a.d, a.w, rec)["certify"]({"words": words})
        os.makedirs(os.path.join(ROOT, "runs", "cwimprove"), exist_ok=True)
        json.dump({"n": a.n, "d": a.d, "w": a.w, "record": rec, "size": len(code), "certified": bool(cert["ok"]), "words": words},
                  open(os.path.join(ROOT, "runs", "cwimprove", f"hit_{a.n}_{a.d}_{a.w}.json"), "w"))
        print(f"  >>> A({a.n},{a.d},{a.w}) >= {len(code)} (record {rec}) certified={cert['ok']}", flush=True)


if __name__ == "__main__":
    main()
