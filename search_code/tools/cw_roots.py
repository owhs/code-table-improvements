"""Root map for Brouwer's constant-weight table A(n,d,w) (data/brouwer/table.json, parsed from Andw.html).
A cell is INHERITED when its lower bound equals what a neighbour gives for free:
  A(n,d,w) >= A(n-1,d,w)      (append a 0)       A(n,d,w) >= A(n-1,d,w-1)   (append a 1)
  A(n,d,w) >= A(n,d+2,w)      (weaker distance)  A(n,d,w) = A(n,d,n-w)      (complement)
otherwise it is a ROOT (a record really constructed there). Open roots (upper bound known and > lower) are ranked by the number
of descendant cells, then softness lower/upper. Only n >= --nmin is listed.
usage: python -I tools/cw_roots.py [--nmin 29] [--top 30]"""
import argparse, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--nmin", type=int, default=29)
ap.add_argument("--top", type=int, default=30)
a = ap.parse_args()
T = json.load(open(os.path.join(ROOT, "data", "brouwer", "table.json"), encoding="utf-8-sig"))
cell = {}
for r in T:
    cell[(r["n"], r["d"], min(r["w"], r["n"] - r["w"]))] = r


def key(n, d, w):
    return (n, d, min(w, n - w))


def parents(k):
    n, d, w = k
    out = [key(n - 1, d, w), key(n - 1, d, w - 1), key(n, d + 2, w)]
    return [p for p in out if p in cell and p != k]


inh, par = {}, {}
for k, r in cell.items():
    best = [(cell[p]["lower"], p) for p in parents(k)]
    ok = [p for v, p in best if v == r["lower"]]
    inh[k] = bool(ok)
    par[k] = ok
desc = {k: set() for k in cell}
for k, ps in par.items():
    for p in ps:
        desc[p].add(k)


def all_desc(k, seen=None):
    seen = seen if seen is not None else set()
    for c in desc[k]:
        if c not in seen:
            seen.add(c)
            all_desc(c, seen)
    return seen


rows = []
for k, r in cell.items():
    if inh[k] or k[0] < a.nmin or r["upper"] is None or r["lower"] >= r["upper"]:
        continue
    rows.append((len(all_desc(k)), r["lower"] / r["upper"], k, r))
rows.sort(key=lambda x: (-x[0], x[1]))
print(f"{len(cell)} cells, {sum(inh.values())} inherited; open roots with n>={a.nmin}: {len(rows)}")
for nd, soft, k, r in rows[:a.top]:
    print(f"A{k}: lower {r['lower']} upper {r['upper']} soft {soft:.2f} desc {nd} source '{r['label']}'")
