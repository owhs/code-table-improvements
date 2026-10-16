"""Root map of a table of lower bounds (new file; plain python). Cells with equal values are often one code seen through puncturing/shortening.
For a cell c with lower bound L(c): c is INHERITED if some neighbour gives it exactly L(c) by a standard rule; otherwise it is a ROOT (the record really lives there).
Rules (q-ary A(n,d)):  A(n,d) >= A(n+1,d+1) (puncture);  A(n,d) >= ceil(A(n+1,d)/q) (shorten);  A(n,d) >= A(n-1,d) (lengthen);  A(n,d) >= A(n,d+1).
Mixed N(m,n,d) (m binary, n ternary coordinates): puncture a binary or a ternary coordinate: N(m,n,d) >= N(m+1,n,d+1), N(m,n+1,d+1); shorten: >= ceil(N(m+1,n,d)/2), ceil(N(m,n+1,d)/3);
lengthen: >= N(m-1,n,d), N(m,n-1,d); relax distance: >= N(m,n,d+1); alphabet: N(m,n,d) >= N(m+1,n-1,d) (a binary coordinate can be read as ternary).
Descendants of a root = cells whose value is reached from it through such steps. High leverage open roots (not proven optimal) with many descendants are the targets.
    python tools/root_map.py [--top 15]"""
import argparse
import json
import math
import os

ap = argparse.ArgumentParser()
ap.add_argument("--top", type=int, default=15)
a = ap.parse_args()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(fn):
    return json.load(open(os.path.join(ROOT, "data", "brouwer", fn), encoding="utf-8-sig"))


def qary_rules(q):
    def f(c, L):
        n, d = c
        out = []
        for key, val in ((  (n + 1, d + 1), lambda v: v), ((n + 1, d), lambda v: math.ceil(v / q)), ((n - 1, d), lambda v: v), ((n, d + 1), lambda v: v)):
            if key in L:
                out.append((key, val(L[key])))
        return out
    return f


def mixed_rules(c, L):
    m, n, d = c
    out = []
    for key, fn in (((m + 1, n, d + 1), lambda v: v), ((m, n + 1, d + 1), lambda v: v), ((m + 1, n, d), lambda v: math.ceil(v / 2)), ((m, n + 1, d), lambda v: math.ceil(v / 3)),
                    ((m - 1, n, d), lambda v: v), ((m, n - 1, d), lambda v: v), ((m, n, d + 1), lambda v: v), ((m + 1, n - 1, d), lambda v: v)):
        if key in L:
            out.append((key, fn(L[key])))
    return out


def analyse(name, cells, rules):
    L = {k: v["lower"] for k, v in cells.items()}
    parent = {}
    for c, lo in L.items():
        srcs = [k for k, val in rules(c, L) if val == lo and L[k] > 0]
        if srcs:
            parent[c] = srcs
    roots = [c for c in L if c not in parent]
    desc = {r: set() for r in roots}
    # propagate: a cell is a descendant of every root reachable by following parents
    def reach(c, seen):
        for p in parent.get(c, []):
            if p in seen:
                continue
            seen.add(p)
            reach(p, seen)
        return seen
    for c in L:
        if c in parent:
            for r in reach(c, set()):
                if r in desc:
                    desc[r].add(c)
    open_roots = [r for r in roots if not cells[r]["exact"]]
    ranked = sorted(open_roots, key=lambda r: (-len(desc[r]), cells[r]["lower"] / max(cells[r]["upper"] or 1, 1)))
    print(f"\n== {name}: {len(L)} cells, {len(parent)} inherited, {len(roots)} roots ({len(open_roots)} open)")
    print("   open roots by number of descendant cells (leverage), then softness lower/upper:")
    for r in ranked[:a.top]:
        up = cells[r]["upper"]
        print(f"   {str(r):14} lower {cells[r]['lower']:>9} upper {str(up):>9}  soft {cells[r]['lower'] / up if up else 0:4.2f}  descendants {len(desc[r]):3d}  e.g. {sorted(desc[r])[:4]}")


T = {}
for q, fn, name in ((3, "ternary.json", "ternary A_3(n,d)"), (4, "quaternary.json", "quaternary A_4(n,d)"), (5, "5ary.json", "5-ary A_5(n,d)")):
    cells = {(r["n"], r["d"]): r for r in load(fn) if r.get("lower")}
    analyse(name, cells, qary_rules(q))
M = load("mixed23.json")
cells = {}
for r in M:
    m = r.get("m", r.get("m_bin"))
    n = r.get("n", r.get("n_ter"))
    if m is None or n is None or not r.get("lower"):
        continue
    cells[(m, n, r["d"])] = r
analyse("mixed binary/ternary N(m,n,d)", cells, mixed_rules)
