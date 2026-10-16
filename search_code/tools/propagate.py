"""Which cells of Brouwer's mixed binary/ternary table (data/brouwer/mixed23.json, N(m,n,d), m binary + n ternary
coordinates, m=0 is the ternary table) are improved for free by standard operations on OUR new codes?
Rules (each is an explicit operation on a code, so every derived value is a real code):
  keep size, d -> d-1 (weaker distance)        N(m,n,d-1)   >= N(m,n,d)
  binary coordinate used as ternary            N(m-1,n+1,d) >= N(m,n,d)
  puncture a coordinate (distance drops by 1)  N(m-1,n,d-1) >= N(m,n,d),  N(m,n-1,d-1) >= N(m,n,d)
  shorten (keep the largest fibre)             N(m-1,n,d)   >= ceil(N(m,n,d)/2),  N(m,n-1,d) >= ceil(N(m,n,d)/3)
usage: python -I tools/propagate.py runs/qary_hits"""
import glob, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
table = {}
for r in json.load(open(os.path.join(ROOT, "data", "brouwer", "mixed23.json"), encoding="utf-8")):
    table[(r["m_bin"], r["n_ter"], r["d"])] = (r["lower"], r["upper"])
ours = {}
for f in glob.glob(os.path.join(sys.argv[1], "HIT_*.json")):
    j = json.load(open(f))
    key = (j["m_bin"], j["n_ter"], j["d"]) if j["table"] == "mixed23" else (0, len(j["qs"]), j["d"])
    ours[key] = max(ours.get(key, 0), len(j["words"]))
val = dict(ours)
src = {k: "found" for k in ours}
changed = True
while changed:
    changed = False
    for (m, n, d), s in list(val.items()):
        cand = []
        if d > 1: cand += [((m, n, d - 1), s, "d-1")]
        if m > 0: cand += [((m - 1, n + 1, d), s, "bin->ter")]
        if d > 1 and m > 0: cand += [((m - 1, n, d - 1), s, "puncture bin")]
        if d > 1 and n > 0: cand += [((m, n - 1, d - 1), s, "puncture ter")]
        if m > 0: cand += [((m - 1, n, d), -(-s // 2), "shorten bin")]
        if n > 0: cand += [((m, n - 1, d), -(-s // 3), "shorten ter")]
        for k, v, how in cand:
            if v > val.get(k, 0):
                val[k], src[k] = v, f"{how} of N{(m, n, d)}={s}"
                changed = True
print("cells where a derived value beats the live table lower bound (and is not one of our found cells):")
rows = []
for k, v in sorted(val.items()):
    if k in ours or k not in table:
        continue
    lo, up = table[k]
    if lo and v > lo:
        rows.append((k, lo, v, up, src[k]))
for k, lo, v, up, how in rows:
    print(f"N{k}: table {lo} (ub {up}) -> {v}   [{how}]")
print(len(rows), "derived improvements")
