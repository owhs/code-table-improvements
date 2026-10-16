"""Harness 4.0 operator library for set-system codes (constant-weight codes A(n,d,w)). Every heavy result is cached.

GROUPS (plan key "group"):
  {"type": "affine", "p": P, "h": H}            x -> a*x + b mod P, a in the order-H subgroup of GF(P)^* (H | P-1)
  {"type": "multiplier", "p": P, "h": H}        x -> a*x mod P (fixes 0)
  {"type": "psl2", "q": Q}                      PSL(2,Q) on Q+1 points (Q prime)
  {"type": "generators", "n": N, "gens": [[perm], ...]}   closure of explicit permutations (capped at 20000)
  {"type": "pg25", "name": NAME}                31 points of the projective plane PG(2,5) (NOT affine maps):
        NAME in singer (31), singer_frob (93), sl_diag_mon (96), flag_stab (2000), point_stab (12000), line_stab (12000)
  optional "fixed": K                           add K extra points fixed by every element (n = base + K)
ORBIT POOLS (plan key "pools": list of {"group": ..., "method": "exact"|"sample", "min_size": m, "max_orbits": k,
  "samples": s}): all internally valid orbits of weight-w words (exact pruned enumeration, or random sampling).
SELECT ({"solver": "exact"|"tabu"|"greedy", "iters": I}): best union of mutually compatible orbits from the pools.
  exact = networkx max-weight clique (proven optimum for the pool; only when the pool graph is small);
  tabu  = compiled add/swap/drop tabu search (any size).
POLISH ({"augment": true, "refill": K}): add every compatible word (exact over all C(n,w) words), then remove K
  orbits at a time and re-solve the hole exactly.
"""
import itertools
import time

import numpy as np

from engine.problems import cwexact as E
from engine.problems import cwgeneric as X
from engine.problems import cwgroups as G
from engine.problems.cwbase import popcount64, to_masks
from h4 import cache
from h4 import clique as CL


def _close(gens, n, cap=20000):
    gens = [tuple(g) for g in gens]
    ident = tuple(range(n))
    seen = {ident}; frontier = [ident]
    while frontier and len(seen) < cap:
        g = frontier.pop()
        for s in gens:
            h = tuple(s[g[i]] for i in range(n))
            if h not in seen:
                seen.add(h); frontier.append(h)
    return np.array(sorted(seen), dtype=np.int64)


def group(spec):
    t = spec["type"]
    if t == "affine":
        P = G.affine_group(spec["p"], spec["h"])
    elif t == "multiplier":
        P = G.multiplier_group(spec["p"], spec["h"])
    elif t == "psl2":
        P = G.psl2_group(spec["q"])
    elif t == "generators":
        P = _close(spec["gens"], spec["n"])
    elif t == "pg25":                                   # collineation subgroups of the projective plane PG(2,5)
        from h4 import pg25
        P = pg25.group(spec["name"])
    else:
        raise ValueError(f"unknown group type {t}")
    if spec.get("fixed"):
        P = G.add_fixed_points(P, int(spec["fixed"]))
    return P


def orbit_pool(n, w, lam, pspec):
    gspec = pspec["group"]
    params = {"n": n, "w": w, "lam": lam, "group": gspec, "method": pspec.get("method", "exact"),
              "samples": pspec.get("samples", 200000), "max_iter": pspec.get("max_iter", 10 ** 11)}

    def compute():
        perms = group(gspec)
        if perms.shape[1] != n:
            raise ValueError(f"group acts on {perms.shape[1]} points, problem has n={n}")
        info = {}
        if params["method"] == "sample":
            orbs = G.random_valid_orbits(n, w, lam, perms, tries=params["samples"], seed=1)
            info["complete"] = False
        elif gspec["type"] == "affine" and not gspec.get("fixed"):
            orbs = E.valid_orbit_reps(gspec["p"], w, lam, gspec["h"], cap=400000, max_iter=params["max_iter"], info=info)
        else:
            orbs = X.valid_orbits(n, w, lam, perms, cap=400000, max_iter=params["max_iter"], info=info)
        return {"orbits": orbs, "complete": info.get("complete", False)}
    res = cache.cached("orbit_pool", params, compute, version=2)   # v2: enumerators store canonical reps only
    orbs = [o for o in res["orbits"] if len(o) >= pspec.get("min_size", 0)]
    if pspec.get("max_orbits") and len(orbs) > pspec["max_orbits"]:
        orbs = sorted(orbs, key=len, reverse=True)[:pspec["max_orbits"]]
    return orbs, res["complete"]


def compat(orbit_lists, lam, same_group, group_specs=None):
    """compatibility matrix over the concatenated pools: fast rep-vs-all inside one group; across affine groups on
    the same Z_p, representatives of the orbits of the COMMON subgroup suffice (exact); otherwise a full check."""
    orbs = [o for lst in orbit_lists for o in lst]
    if same_group:
        return orbs, G.compat_matrix(orbs, lam)
    if group_specs and all(g["type"] == "affine" and not g.get("fixed") for g in group_specs) \
            and len({g["p"] for g in group_specs}) == 1:
        import math
        p = group_specs[0]["p"]
        hk = 0
        for g in group_specs:
            hk = math.gcd(hk, g["h"])
        K = G.affine_group(p, hk)
        reps_flat, rstarts, rlens = [], [], []
        for o in orbs:                                   # split each orbit into K-orbits, keep one rep of each
            seen, reps = set(), []
            for m in o:
                m = int(m)
                if m in seen:
                    continue
                ko = G.orbit(m, K)
                seen.update(int(x) for x in ko)
                reps.append(m)
            rstarts.append(len(reps_flat)); rlens.append(len(reps)); reps_flat += reps
        lens = np.array([len(o) for o in orbs], dtype=np.int64)
        starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
        flat = np.concatenate(orbs).astype(np.uint64)
        return orbs, _krep_compat(flat, starts, lens, np.array(reps_flat, dtype=np.uint64),
                                  np.array(rstarts, dtype=np.int64), np.array(rlens, dtype=np.int64), lam)
    lens = np.array([len(o) for o in orbs], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
    flat = np.concatenate(orbs).astype(np.uint64)
    return orbs, _full_compat(flat, starts, lens, lam)


from numba import njit  # noqa: E402


@njit(cache=True)
def _krep_compat(flat, starts, lens, reps, rstarts, rlens, lam):
    k = lens.shape[0]
    C = np.zeros((k, k), dtype=np.bool_)
    for i in range(k):
        for j in range(i + 1, k):
            ok = True
            for a in range(rstarts[i], rstarts[i] + rlens[i]):
                x = reps[a]
                for b in range(starts[j], starts[j] + lens[j]):
                    if popcount64(x & flat[b]) > lam:
                        ok = False
                        break
                if not ok:
                    break
            C[i, j] = ok; C[j, i] = ok
    return C


@njit(cache=True)
def _full_compat(flat, starts, lens, lam):
    k = lens.shape[0]
    C = np.zeros((k, k), dtype=np.bool_)
    for i in range(k):
        for j in range(i + 1, k):
            ok = True
            for a in range(starts[i], starts[i] + lens[i]):
                x = flat[a]
                for b in range(starts[j], starts[j] + lens[j]):
                    if popcount64(x & flat[b]) > lam:
                        ok = False
                        break
                if not ok:
                    break
            C[i, j] = ok; C[j, i] = ok
    return C


def select(orbs, C, sspec, effort):
    wts = np.array([len(o) for o in orbs], dtype=np.int64)
    solver = sspec.get("solver", "tabu")
    edges = int(C.sum() // 2)
    if solver == "exact" and len(orbs) <= 12000 * effort and edges <= 4_000_000 * effort:
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from((i, {"weight": int(wts[i])}) for i in range(len(orbs)))
        ii, jj = np.nonzero(np.triu(C, 1))
        g.add_edges_from(zip(ii.tolist(), jj.tolist()))
        idx, wsum = nx.max_weight_clique(g, weight="weight")
        return list(idx), int(wsum), "exact"
    if solver == "greedy":
        order = np.argsort(-wts)
        allowed = np.ones(len(orbs), dtype=bool); idx = []
        for i in order:
            if allowed[i]:
                idx.append(int(i)); allowed &= C[i]
        return idx, int(wts[idx].sum()), "greedy"
    iters = int(sspec.get("iters", 50000) * effort)
    best = ([], 0)
    for s in range(max(1, int(effort))):
        idx, wsum = CL.tabu(C, wts, iters, seed=s)
        if wsum > best[1]:
            best = (idx, wsum)
    return best[0], best[1], "tabu"


def augment(n, w, lam, masks):
    from tools.augment import _compatible_words
    extra = _compatible_words(n, w, lam, np.asarray(masks, dtype=np.uint64), 200000)
    if len(extra) == 0:
        return list(masks)
    k = len(extra)
    C = np.zeros((k, k), dtype=bool)
    for i in range(k):
        for j in range(i + 1, k):
            C[i, j] = C[j, i] = popcount64(extra[i] & extra[j]) <= lam
    idx, _ = CL.tabu(C, np.ones(k, dtype=np.int64), 20000, 0)
    return list(masks) + [int(extra[i]) for i in idx]


def refill(n, w, lam, blocks, k, budget_s):
    """remove k blocks at a time, refill the hole exactly over all compatible words (tabu on the hole graph)"""
    from tools.augment import _compatible_words
    t0 = time.time()
    for drop in itertools.combinations(range(len(blocks)), k):
        if time.time() - t0 > budget_s:
            break
        keep = [m for i, b in enumerate(blocks) if i not in drop for m in b]
        removed = sum(len(blocks[i]) for i in drop)
        cands = _compatible_words(n, w, lam, np.array(keep, dtype=np.uint64), 400000)
        if len(cands) <= removed:
            continue
        kk = len(cands)
        if kk > 20000:
            continue
        C = np.zeros((kk, kk), dtype=bool)
        for i in range(kk):
            x = cands[i]
            for j in range(i + 1, kk):
                C[i, j] = C[j, i] = popcount64(x & cands[j]) <= lam
        idx, got = CL.tabu(C, np.ones(kk, dtype=np.int64), 100000, 0)
        if got > removed:
            new_blocks = [b for i, b in enumerate(blocks) if i not in drop] + [[int(cands[i]) for i in idx]]
            return new_blocks, True
    return blocks, False
