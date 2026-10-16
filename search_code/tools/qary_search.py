"""Prescribed-automorphism search for unrestricted ternary / mixed binary-ternary codes (Brouwer tables).
    python tools/qary_search.py --cell T:16,7 [--secs 600] [--fam any] [--seed 1] [--orb 30000] [--vmax 12000]
                                [--refill 120] [--out runs/qary] [--threads 4]
    --cell T:n,d   A_3(n,d)            (data/brouwer/ternary.json)
    --cell M:m,n,d N(m,n,d), m binary + n ternary coordinates (data/brouwer/mixed23.json); binary coordinates first
Loop: sample a group G of isometries (coordinate permutations within each alphabet block, affine symbol maps
x -> a x + b per coordinate, translations by a random additive code K = union-of-cosets reduction), enumerate all
orbits exactly, keep the internally valid ones, build the exact orbit compatibility graph, solve max-weight clique
(tabu), augment with every compatible word. The best code is polished by remove-and-refill. Results:
    <out>/best_<T|M>_<params>.json  (always)   <out>/HIT_<table>_<params>.json  (if it beats the live lower bound)
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from engine.problems import qary as QA  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_cell(spec):
    kind, rest = spec.split(":")
    nums = [int(x) for x in rest.split(",")]
    if kind == "T":
        n, d = nums
        rows = json.load(open(os.path.join(ROOT, "data", "brouwer", "ternary.json"), encoding="utf-8"))
        r = next((r for r in rows if r["n"] == n and r["d"] == d), None)
        qs = [3] * n
        name = f"ternary_{n}_{d}"
    else:
        m, n, d = nums
        rows = json.load(open(os.path.join(ROOT, "data", "brouwer", "mixed23.json"), encoding="utf-8"))
        r = next((r for r in rows if r["m_bin"] == m and r["n_ter"] == n and r["d"] == d), None)
        qs = [2] * m + [3] * n
        name = f"mixed23_{m}_{n}_{d}"
    lower = r["lower"] if r else None
    upper = r["upper"] if r else None
    return qs, d, lower, upper, name


# ----------------------------------------------------------------------------------------------- group sampling
def _blocks(qs):
    b2 = [i for i, q in enumerate(qs) if q == 2]
    b3 = [i for i, q in enumerate(qs) if q == 3]
    return [b for b in (b2, b3) if b]


def _primes_upto(m):
    return [p for p in range(2, m + 1) if all(p % k for k in range(2, int(p ** 0.5) + 1))]


def _prim_root(p):
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in _primes_upto(p - 1) if (p - 1) % q == 0):
            return g
    return 1


def _struct_gens(rng, L):
    """coordinate permutation generators on points 0..L-1 (structured groups); returns list of perms"""
    opts = []
    opts.append([[(x + 1) % L for x in range(L)]])                       # cyclic
    for p in _primes_upto(L):
        if p >= 5 and L - p <= 2:
            g = _prim_root(p)
            sub = [h for h in range(1, p) if (p - 1) % h == 0]
            h = int(rng.choice(sub))
            a = pow(g, (p - 1) // h, p)
            t = [(x + 1) % p if x < p else x for x in range(L)]
            mlt = [(a * x) % p if x < p else x for x in range(L)]
            opts.append([t, mlt])                                         # AGL-type on Z_p (+ fixed points)
            opts.append([mlt])                                            # multiplier only
        if p >= 5 and L - p - 1 in (0, 1):                                # PSL(2,p) on p+1 points (+ fixed)
            INF = p
            t = [(x + 1) % p if x < p else x for x in range(L)]
            s = []
            for x in range(L):
                if x == INF:
                    s.append(0)
                elif x == 0:
                    s.append(INF)
                elif x < p:
                    s.append((-pow(x, p - 2, p)) % p)
                else:
                    s.append(x)
            opts.append([t, s])
    for k in (3, 4):
        if L in (2 ** k, 2 ** k + 1):                                     # AGL(1, 2^k) on GF(2^k)
            poly = {3: 0b1011, 4: 0b10011}[k]
            def mul2(x):
                x <<= 1
                return x ^ poly if x >> k else x
            t = [x ^ 1 if x < 2 ** k else x for x in range(L)]
            mlt = [mul2(x) if x < 2 ** k else x for x in range(L)]
            opts.append([t, mlt])
            opts.append([mlt])
    return opts[int(rng.integers(len(opts)))]


def _rand_cycle_perm(rng, L):
    o = int(rng.integers(2, L + 1))
    pts = list(rng.permutation(L))
    p = list(range(L))
    pos = 0
    while pos < L:
        divs = [k for k in range(1, o + 1) if o % k == 0 and k <= L - pos]
        c = o if (o <= L - pos and rng.random() < 0.7) else int(rng.choice(divs))
        cyc = pts[pos:pos + c]
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            p[int(a)] = int(b)
        pos += c
    return p


def _decorate(rng, qs, perm, p_scale, p_shift):
    n = len(qs)
    scale = [int(rng.choice([1, 2])) if (qs[i] == 3 and rng.random() < p_scale) else 1 for i in range(n)]
    shift = [int(rng.integers(qs[i])) if rng.random() < p_shift else 0 for i in range(n)]
    return [list(map(int, perm)), scale, shift]


def _block_perm(rng, qs, mode):
    n = len(qs)
    perm = list(range(n))
    for b in _blocks(qs):
        L = len(b)
        if L < 2 or rng.random() < 0.15:
            continue
        if mode == "struct":
            gl = _struct_gens(rng, L)
            lp = gl[int(rng.integers(len(gl)))]
        else:
            lp = _rand_cycle_perm(rng, L)
        for a in range(L):
            perm[b[a]] = b[lp[a]]
    return perm


def _struct_perms(rng, qs):
    """a list of whole-space permutations (one structured group per block, generators aligned)"""
    n = len(qs)
    gens = []
    for b in _blocks(qs):
        L = len(b)
        if L < 3 or rng.random() < 0.15:
            continue
        for lp in _struct_gens(rng, L):
            perm = list(range(n))
            for a in range(L):
                perm[b[a]] = b[lp[a]]
            gens.append(perm)
    return gens


def _rand_translation(rng, qs, d):
    n = len(qs)
    while True:
        w = int(rng.integers(min(d, n), n + 1))
        sup = rng.choice(n, size=w, replace=False)
        c = [0] * n
        for i in sup:
            c[int(i)] = int(rng.integers(1, qs[int(i)]))
        return c


def _coord_orbits(n, perms):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for p in perms:
        for i in range(n):
            ra, rb = find(i), find(p[i])
            if ra != rb:
                parent[ra] = rb
    orbs = {}
    for i in range(n):
        orbs.setdefault(find(i), []).append(i)
    return list(orbs.values())


def _invariant_translation(rng, qs, d, perms):
    """vector constant on a random set of coordinate orbits of the permutation parts (weight >= d if possible)"""
    n = len(qs)
    orbs = _coord_orbits(n, perms)
    for _ in range(50):
        c = [0] * n
        for o in orbs:
            if rng.random() < 0.6:
                val = int(rng.integers(1, qs[o[0]]))
                for i in o:
                    c[i] = val
        if sum(1 for x in c if x) >= d:
            return c
    return None


def _deco(rng, qs, perm, ps, pt):
    mode = rng.choice(["none", "neg", "one", "negone", "rand"], p=[0.25, 0.2, 0.2, 0.15, 0.2])
    n = len(qs)
    if mode == "rand":
        return _decorate(rng, qs, perm, ps, pt)
    scale = [2 if (qs[i] == 3 and mode in ("neg", "negone")) else 1 for i in range(n)]
    shift = [1 if mode in ("one", "negone") else 0 for i in range(n)]
    return [list(map(int, perm)), scale, shift]


def _commuting_perm(rng, n, perm, qs=None):
    """a permutation commuting with perm: rotate each cycle independently, or swap two equal-length cycles"""
    seen, cycles = set(), []
    for i in range(n):
        if i not in seen:
            c = [i]; seen.add(i)
            j = perm[i]
            while j != i:
                c.append(j); seen.add(j); j = perm[j]
            cycles.append(c)
    q = list(range(n))
    if rng.random() < 0.5:
        for c in cycles:
            r = int(rng.integers(len(c)))
            for k, x in enumerate(c):
                q[x] = c[(k + r) % len(c)]
    else:
        bylen = {}
        for c in cycles:
            bylen.setdefault((len(c), None if qs is None else qs[c[0]]), []).append(c)
        pairs = [v for v in bylen.values() if len(v) >= 2]
        if not pairs:
            return None
        grp = pairs[int(rng.integers(len(pairs)))]
        a, b = rng.choice(len(grp), size=2, replace=False)
        ca, cb = grp[int(a)], grp[int(b)]
        r = int(rng.integers(len(ca)))
        for k in range(len(ca)):
            q[ca[k]] = cb[(k + r) % len(cb)]
            q[cb[k]] = ca[(k - r) % len(ca)]
    return q


def sample_group(rng, qs, d, min_order, cap, fam="any"):
    fams = ["mono", "mono+K", "struct", "struct+K", "K", "mono2"]
    f = fams[int(rng.integers(len(fams)))] if fam == "any" else fam
    gens = []
    ps, pt = float(rng.random()), float(rng.random()) * 0.5
    n = len(qs)
    if f.startswith("mono"):
        p1 = _block_perm(rng, qs, "rand")
        gens.append(_deco(rng, qs, p1, ps, pt))
        if f == "mono2":
            p2 = _commuting_perm(rng, n, p1, qs)
            if p2 is not None:
                gens.append(_deco(rng, qs, p2, ps, pt))
    if f.startswith("struct"):
        for perm in _struct_perms(rng, qs):
            gens.append(_deco(rng, qs, perm, ps, pt))
        if not gens:
            gens.append(_deco(rng, qs, _block_perm(rng, qs, "rand"), ps, pt))
    if f.endswith("K"):
        k = int(rng.integers(1, 4))
        perms = [g[0] for g in gens]
        for _ in range(k):
            c = _invariant_translation(rng, qs, d, perms)
            if c is not None:
                gens.append([list(range(n)), [1] * n, c])
    if f == "K":
        k = int(rng.integers(1, 6))
        for _ in range(k):
            gens.append([list(range(n)), [1] * n, _rand_translation(rng, qs, d)])
    els = [QA.elem(qs, *g) for g in gens]
    G = QA.closure(qs, els, cap=cap)
    if G is None:
        return None, gens, f, None
    order = G[0].shape[0]
    if order < min_order:
        return None, gens, f, order
    return G, gens, f, order


# ----------------------------------------------------------------------------------------------- one group
def run_group(sp, d, G, vmax, iters, rng):
    t0 = time.time()
    oid, reps, sizes, info = QA.enum_orbits(sp, G[0], G[1], d, vcap=4 * vmax)
    V = len(reps)
    if info["overflow"]:
        info["note"] = "valid-orbit overflow"
    if V == 0:
        return None, info
    if V > vmax:
        keep = np.sort(rng.choice(V, size=vmax, replace=False))
        remap = np.full(V, -2, dtype=np.int32)
        remap[keep] = np.arange(vmax, dtype=np.int32)
        m = oid >= 0
        oid[m] = remap[oid[m]]
        reps, sizes = reps[keep], sizes[keep]
        V = vmax
        info["subsampled"] = True
    t1 = time.time()
    Cf = QA.orbit_conflicts(sp, oid, reps, sizes, d)
    C = ~Cf
    np.fill_diagonal(C, False)
    t2 = time.time()
    idx, w = QA.max_clique(C, sizes, iters=iters, seeds=2, exact_limit=60)
    chosen = np.zeros(V, dtype=bool)
    chosen[idx] = True
    m = oid >= 0
    words_idx = np.nonzero(m)[0][chosen[oid[m]]]
    info.update(V=int(V), orbit_size_max=int(sizes.max()), enum_s=round(t1 - t0, 2), graph_s=round(t2 - t1, 2),
                clique_s=round(time.time() - t2, 2), clique=int(w), reps=[sp.digits(int(reps[i])) for i in idx])
    return words_idx, info


def save(path, sp, d, cell, words, construction):
    ok, md, msg = QA.check(sp, d, words)
    out = {"table": cell["table"], "cell": cell["spec"], "qs": [int(q) for q in sp.qs], "d": d,
           "m_bin": int((sp.qs == 2).sum()), "n_ter": int((sp.qs == 3).sum()), "size": len(words),
           "live_lower": cell["lower"], "live_upper": cell["upper"], "checker1": msg, "checker1_ok": bool(ok),
           "construction": construction, "words": words}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--secs", type=float, default=600)
    ap.add_argument("--fam", default="any")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--orb", type=float, default=30000, help="target max number of orbits (sets min |G|)")
    ap.add_argument("--vmax", type=int, default=12000)
    ap.add_argument("--cap", type=int, default=20000)
    ap.add_argument("--iters", type=int, default=60000)
    ap.add_argument("--refill", type=float, default=120)
    ap.add_argument("--rho", type=int, default=None)
    ap.add_argument("--out", default=os.path.join(ROOT, "runs", "qary"))
    ap.add_argument("--stop_at", type=int, default=None, help="stop the group loop once this size is reached")
    a = ap.parse_args()
    qs, d, lower, upper, cname = load_cell(a.cell)
    name = f"{cname}_s{a.seed}"
    sp = QA.Space(qs)
    cell = {"spec": a.cell, "table": "ternary" if a.cell.startswith("T") else "mixed23", "lower": lower, "upper": upper}
    rng = np.random.default_rng(a.seed)
    log = lambda s: print(s, flush=True)  # noqa: E731
    log(f"{name}: Q={sp.Q} d={d} live [{lower},{upper}] ball(d-1)={sp.ball_size(d - 1)}")
    min_order = max(1, int(math.ceil(sp.Q / a.orb)))
    best = (0, None, None)
    t0 = time.time(); tried = 0
    stop_at = a.stop_at
    while time.time() - t0 < a.secs:
        try:
            G, gens, fam, order = sample_group(rng, qs, d, min_order, a.cap, a.fam)
        except AssertionError:
            continue
        if G is None:
            continue
        tried += 1
        words_idx, info = run_group(sp, d, G, a.vmax, a.iters, rng)
        if words_idx is None:
            continue
        size0 = len(words_idx)
        added = 0
        if size0 >= 0.9 * best[0]:
            cov = QA.Cover(sp, d, words_idx)
            added = QA.augment(cov)
            words_idx = np.array(sorted(cov.code), dtype=np.int64)
            del cov
        size = len(words_idx)
        log(f"[{time.time() - t0:6.0f}s] fam={fam:8s} |G|={order:6d} orbits={info['orbits']} V={info['V']} "
            f"clique={size0} +aug {added} = {size}  (enum {info['enum_s']}s graph {info['graph_s']}s "
            f"clique {info['clique_s']}s){' SUB' if info.get('subsampled') else ''}")
        if size > best[0]:
            best = (size, words_idx, {"family": fam, "group_order": order, "generators": gens,
                                      "generator_semantics": "y[perm[i]] = (scale[i]*x[i] + shift[i]) mod q_i",
                                      "orbit_representatives": info["reps"], "clique_size": size0,
                                      "augmented": added, "seed": a.seed})
            words = [sp.digits(int(i)) for i in words_idx]
            save(os.path.join(a.out, f"best_{name}.json"), sp, d, cell, words, best[2])
            log(f"   NEW BEST {size}")
            if stop_at and size >= stop_at:
                break
    log(f"groups tried: {tried}; best {best[0]}")
    if best[1] is None:
        return
    if a.refill > 0:
        cov = QA.Cover(sp, d, best[1])
        before = len(cov.code)
        gains, moves = QA.refill(cov, a.refill, rho=a.rho, rng=rng, log=log)
        log(f"refill: {moves} moves, {before} -> {len(cov.code)}")
        if len(cov.code) > before:
            words = [sp.digits(int(i)) for i in sorted(cov.code)]
            c = dict(best[2]); c["refill"] = f"remove-and-refill {a.refill}s: {before} -> {len(cov.code)}"
            best = (len(cov.code), None, c)
            save(os.path.join(a.out, f"best_{name}.json"), sp, d, cell, words, c)
    final = json.load(open(os.path.join(a.out, f"best_{name}.json")))
    log(f"FINAL {name}: {final['size']} ({final['checker1']}) live lower {lower}")
    if lower is not None and final["size"] > lower and final["checker1_ok"]:
        hp = os.path.join(a.out, f"HIT_{cname}.json")
        prev = json.load(open(hp))["size"] if os.path.exists(hp) else 0
        if final["size"] > prev:
            with open(hp, "w") as f:
                json.dump(final, f)
        log(f"*** HIT {name}: {final['size']} > {lower} -> {hp}")


if __name__ == "__main__":
    main()
