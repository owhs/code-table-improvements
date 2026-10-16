"""Mixed-alphabet (q-ary / binary-ternary) unrestricted codes: domain adapter for the orbit + clique machinery.

SPACE  qs = alphabet size per coordinate (2 or 3), e.g. [2]*m + [3]*n for Brouwer's N(m,n,d) (binary coordinates first)
       or [3]*n for A_3(n,d). A word is a digit vector x (x_i in 0..q_i-1); its INDEX is sum x_i * radix_i
       (radix_i = prod_{j<i} q_j), its PACKED form is a uint64 with bit i set if x_i == 1 and bit 32+i if x_i == 2, so
       Hamming distance = popcount((a ^ b | (a ^ b) >> 32) & 0xffffffff).  n <= 32, prod(q) <= ~6e7 (index arrays).
GROUP  isometries g = (perm, sym): y_{perm[i]} = sym[i, x_i]; perm maps coordinates to coordinates of the same
       alphabet size, sym[i] is a permutation of the alphabet. Arrays perms (|G|, n) and syms (|G|, n, 3).
       elem(qs, perm, scale, shift) builds the affine element y_{perm[i]} = scale_i x_i + shift_i (mod q_i);
       closure(qs, gens, cap) -> (perms, syms) or None if the group is larger than cap.
       Translations by a code K (x -> x + c) give orbits = cosets of K (union-of-cosets constructions);
       affine symbol maps + coordinate permutations normalising K give K semidirect P.
ORBITS enum_orbits(sp, perms, syms, d, vcap) -> oid (int32 per word: >= 0 valid orbit id, -2 invalid orbit),
       reps (indices), sizes. Valid = every two distinct words of the orbit at distance >= d (rep vs its images).
GRAPH  orbit_conflicts(sp, oid, reps, d) -> bool matrix of incompatible valid-orbit pairs (exact: rep_i against all words
       of orbit j, by Hamming-ball lookup or by a flat pairwise scan, whichever is cheaper).
POLISH Cover (word-level): cnt[idx] = number of codewords within distance d-1; free words have cnt == 0.
       augment(sp, d, code) adds a max clique of mutually compatible free words; refill(...) removes the codewords
       around a random centre and re-solves the hole (accepting gains, and equal-size changes as plateau moves).
CHECK  check(sp, d, words) -> (ok, min distance, message): compiled pairwise check (checker 1; tools/qary_verify.py is
       the independent plain-Python checker 2).
"""
import math
import time

import numpy as np
from numba import njit, prange

from engine.problems.cwbase import popcount64

LOW32 = np.uint64(0xFFFFFFFF)


# ----------------------------------------------------------------------------------------------- space
class Space:
    def __init__(self, qs):
        self.qs = np.array(qs, dtype=np.int64)
        self.n = len(qs)
        assert self.n <= 32 and all(q in (2, 3) for q in qs)
        self.radix = np.ones(self.n, dtype=np.int64)
        for i in range(1, self.n):
            self.radix[i] = self.radix[i - 1] * self.qs[i - 1]
        self.Q = int(np.prod(self.qs))

    def index(self, word):
        return int(sum(int(x) * int(r) for x, r in zip(word, self.radix)))

    def digits(self, idx):
        out = []
        for q in self.qs:
            out.append(int(idx % q)); idx //= int(q)
        return out

    def ball_size(self, r):
        # number of words within distance r (coordinates of mixed alphabets: polynomial product)
        poly = [1]
        for q in self.qs:
            new = [0] * (len(poly) + 1)
            for k, c in enumerate(poly):
                new[k] += c; new[k + 1] += c * (int(q) - 1)
            poly = new
        return sum(poly[:r + 1])


@njit(cache=True)
def qdist(a, b):
    z = a ^ b
    return popcount64((z | (z >> np.uint64(32))) & LOW32)


@njit(cache=True)
def _decode(idx, qs, x):
    t = idx
    for i in range(qs.shape[0]):
        x[i] = t % qs[i]
        t //= qs[i]


@njit(cache=True)
def _pack_digits(x):
    pk = np.uint64(0)
    for i in range(x.shape[0]):
        if x[i] == 1:
            pk |= np.uint64(1) << np.uint64(i)
        elif x[i] == 2:
            pk |= np.uint64(1) << np.uint64(32 + i)
    return pk


@njit(cache=True)
def pack_indices(idxs, qs):
    out = np.empty(idxs.shape[0], dtype=np.uint64)
    x = np.empty(qs.shape[0], dtype=np.int64)
    for k in range(idxs.shape[0]):
        _decode(idxs[k], qs, x)
        out[k] = _pack_digits(x)
    return out


# ----------------------------------------------------------------------------------------------- groups
def elem(qs, perm, scale=None, shift=None):
    n = len(qs)
    perm = list(perm)
    scale = [1] * n if scale is None else list(scale)
    shift = [0] * n if shift is None else list(shift)
    sym = np.zeros((n, 3), dtype=np.int64)
    for i in range(n):
        q = qs[i]
        assert qs[perm[i]] == q, "permutation must preserve alphabet sizes"
        for a in range(q):
            sym[i, a] = (scale[i] * a + shift[i]) % q
        assert sorted(sym[i, :q]) == list(range(q))
    return np.array(perm, dtype=np.int64), sym


def translation(qs, c):
    return elem(qs, range(len(qs)), None, c)


def _key(p, s):
    return p.tobytes() + s.tobytes()


def _compose(h, g, binmask=None):
    """h after g"""
    ph, sh = h; pg, sg = g
    p = ph[pg]
    s = np.take_along_axis(sh[pg], sg, axis=1)
    if binmask is not None:
        s[binmask, 2] = 0                 # unused column of binary coordinates: keep canonical
    return p, s


def closure(qs, gens, cap=50000):
    n = len(qs)
    binmask = np.array([q == 2 for q in qs])
    ident = elem(qs, range(n))
    seen = {_key(*ident): ident}
    frontier = [ident]
    while frontier:
        g = frontier.pop()
        for s in gens:
            h = _compose(s, g, binmask)
            k = _key(*h)
            if k not in seen:
                seen[k] = h
                if len(seen) > cap:
                    return None
                frontier.append(h)
    els = list(seen.values())
    return np.array([e[0] for e in els], dtype=np.int64), np.array([e[1] for e in els], dtype=np.int64)


def cyclic_perm(n, cycles):
    """permutation of range(n) with the given cycles (lists of points); other points fixed"""
    p = list(range(n))
    for c in cycles:
        for a, b in zip(c, c[1:] + c[:1]):
            p[a] = b
    return p


# ----------------------------------------------------------------------------------------------- orbits
@njit(cache=True)
def _enum_orbits(Q, qs, radix, perms, syms, d, oid, vcap, reps, sizes):
    n = qs.shape[0]
    G = perms.shape[0]
    tmp = np.empty(G, dtype=np.int64)
    x = np.empty(n, dtype=np.int64)
    nv = 0; norb = 0; overflow = 0
    for w in range(Q):
        if oid[w] != -1:
            continue
        _decode(w, qs, x)
        pw = _pack_digits(x)
        k = 0; mind = n + 1
        for g in range(G):
            idx = 0
            pk = np.uint64(0)
            for i in range(n):
                v = syms[g, i, x[i]]
                j = perms[g, i]
                idx += v * radix[j]
                if v == 1:
                    pk |= np.uint64(1) << np.uint64(j)
                elif v == 2:
                    pk |= np.uint64(1) << np.uint64(32 + j)
            if oid[idx] == -1:
                oid[idx] = -3
                tmp[k] = idx; k += 1
                if idx != w:
                    dd = qdist(pw, pk)
                    if dd < mind:
                        mind = dd
        norb += 1
        lab = -2
        if mind >= d:
            if nv < vcap:
                lab = nv
                reps[nv] = w; sizes[nv] = k
                nv += 1
            else:
                overflow += 1
        for t in range(k):
            oid[tmp[t]] = lab
    return nv, norb, overflow


def enum_orbits(sp, perms, syms, d, vcap=200000):
    oid = np.full(sp.Q, -1, dtype=np.int32)
    reps = np.zeros(vcap, dtype=np.int64)
    sizes = np.zeros(vcap, dtype=np.int64)
    nv, norb, overflow = _enum_orbits(sp.Q, sp.qs, sp.radix, perms, syms, d, oid, vcap, reps, sizes)
    return oid, reps[:nv].copy(), sizes[:nv].copy(), {"orbits": int(norb), "valid": int(nv), "overflow": int(overflow)}


# ----------------------------------------------------------------------------------------------- Hamming-ball walks
@njit(cache=True)
def _ball_conflicts_row(i, idx0, qs, radix, r, oid, row, x, cp, ca, ix):
    n = qs.shape[0]
    _decode(idx0, qs, x)
    k = 0; cp[0] = 0; ca[0] = 0; ix[0] = idx0
    while k >= 0:
        p = cp[k]; a = ca[k] + 1
        if a >= qs[p]:
            p += 1; a = 1
        if p >= n:
            k -= 1
            continue
        cp[k] = p; ca[k] = a
        nd = (x[p] + a) % qs[p]
        idx = ix[k] + (nd - x[p]) * radix[p]
        j = oid[idx]
        if j >= 0 and j != i:
            row[j] = True
        if k + 1 < r and p + 1 < n:
            ix[k + 1] = idx; cp[k + 1] = p + 1; ca[k + 1] = 0; k += 1


@njit(cache=True, parallel=True)
def _conflicts_ball(reps, qs, radix, r, oid, Cf):
    V = reps.shape[0]
    n = qs.shape[0]
    for i in prange(V):
        x = np.empty(n, dtype=np.int64)
        cp = np.zeros(r + 1, dtype=np.int64); ca = np.zeros(r + 1, dtype=np.int64); ix = np.zeros(r + 1, dtype=np.int64)
        _ball_conflicts_row(i, reps[i], qs, radix, r, oid, Cf[i], x, cp, ca, ix)


@njit(cache=True, parallel=True)
def _conflicts_pair(rep_pk, flat_pk, owner, d, Cf):
    V = rep_pk.shape[0]
    for i in prange(V):
        a = rep_pk[i]
        for t in range(flat_pk.shape[0]):
            j = owner[t]
            if j != i and qdist(a, flat_pk[t]) < d:
                Cf[i, j] = True


def orbit_conflicts(sp, oid, reps, sizes, d):
    V = len(reps)
    Cf = np.zeros((V, V), dtype=np.bool_)
    W = int(sizes.sum())
    if V == 0:
        return Cf
    ball = sp.ball_size(d - 1)
    if W <= ball:
        flat = np.nonzero(oid >= 0)[0].astype(np.int64)
        owner = oid[flat].astype(np.int64)
        _conflicts_pair(pack_indices(reps, sp.qs), pack_indices(flat, sp.qs), owner, d, Cf)
    else:
        _conflicts_ball(reps, sp.qs, sp.radix, d - 1, oid, Cf)
    Cf |= Cf.T
    return Cf


# ----------------------------------------------------------------------------------------------- word-level cover
@njit(cache=True)
def _ball_add(idx0, qs, radix, r, cnt, delta):
    n = qs.shape[0]
    x = np.empty(n, dtype=np.int64)
    cp = np.zeros(r + 1, dtype=np.int64); ca = np.zeros(r + 1, dtype=np.int64); ix = np.zeros(r + 1, dtype=np.int64)
    _decode(idx0, qs, x)
    cnt[idx0] += delta
    if r == 0:
        return
    k = 0; cp[0] = 0; ca[0] = 0; ix[0] = idx0
    while k >= 0:
        p = cp[k]; a = ca[k] + 1
        if a >= qs[p]:
            p += 1; a = 1
        if p >= n:
            k -= 1
            continue
        cp[k] = p; ca[k] = a
        nd = (x[p] + a) % qs[p]
        idx = ix[k] + (nd - x[p]) * radix[p]
        cnt[idx] += delta
        if k + 1 < r and p + 1 < n:
            ix[k + 1] = idx; cp[k + 1] = p + 1; ca[k + 1] = 0; k += 1


@njit(cache=True)
def _ball_free(idx0, qs, radix, r, cnt, out, m):
    """append to out[m:] the words of the ball around idx0 with cnt == 0; returns new m"""
    n = qs.shape[0]
    x = np.empty(n, dtype=np.int64)
    cp = np.zeros(r + 1, dtype=np.int64); ca = np.zeros(r + 1, dtype=np.int64); ix = np.zeros(r + 1, dtype=np.int64)
    _decode(idx0, qs, x)
    if cnt[idx0] == 0 and m < out.shape[0]:
        out[m] = idx0; m += 1
    if r == 0:
        return m
    k = 0; cp[0] = 0; ca[0] = 0; ix[0] = idx0
    while k >= 0:
        p = cp[k]; a = ca[k] + 1
        if a >= qs[p]:
            p += 1; a = 1
        if p >= n:
            k -= 1
            continue
        cp[k] = p; ca[k] = a
        nd = (x[p] + a) % qs[p]
        idx = ix[k] + (nd - x[p]) * radix[p]
        if cnt[idx] == 0 and m < out.shape[0]:
            out[m] = idx; m += 1
        if k + 1 < r and p + 1 < n:
            ix[k + 1] = idx; cp[k + 1] = p + 1; ca[k + 1] = 0; k += 1
    return m


@njit(cache=True)
def _build_cover(code, qs, radix, r, cnt):
    for c in range(code.shape[0]):
        _ball_add(code[c], qs, radix, r, cnt, 1)


@njit(cache=True)
def _pair_compat(pk, d):
    k = pk.shape[0]
    C = np.zeros((k, k), dtype=np.bool_)
    for i in range(k):
        for j in range(i + 1, k):
            if qdist(pk[i], pk[j]) >= d:
                C[i, j] = True; C[j, i] = True
    return C


@njit(cache=True)
def _min_dist(pk):
    m = 64
    for i in range(pk.shape[0]):
        for j in range(i + 1, pk.shape[0]):
            dd = qdist(pk[i], pk[j])
            if dd < m:
                m = dd
    return m


def check(sp, d, words):
    idx = np.array([sp.index(w) for w in words], dtype=np.int64)
    for w in words:
        if len(w) != sp.n or any(not (0 <= int(v) < int(q)) for v, q in zip(w, sp.qs)):
            return False, None, "bad word"
    if len(set(idx.tolist())) != len(idx):
        return False, 0, "duplicate words"
    md = int(_min_dist(pack_indices(idx, sp.qs))) if len(idx) > 1 else 64
    return md >= d, md, f"{len(idx)} words, min distance {md}"


# ----------------------------------------------------------------------------------------------- clique helper
def max_clique(C, wts, iters=200000, seeds=3, exact_limit=80, time_limit=None):
    from h4 import clique as CL
    k = len(wts)
    if k == 0:
        return [], 0
    if k <= exact_limit:
        idx, w, proven = CL.exact(C, wts, max_nodes=5_000_000)
        if proven:
            return idx, w
    best = ([], 0)
    t0 = time.time()
    for s in range(seeds):
        idx, w = CL.tabu_fast(C, wts, iters, s)
        if w > best[1]:
            best = (idx, w)
        if time_limit and time.time() - t0 > time_limit:
            break
    return best


# ----------------------------------------------------------------------------------------------- polish
class Cover:
    def __init__(self, sp, d, code_idx):
        self.sp, self.d = sp, d
        self.r = d - 1
        self.cnt = np.zeros(sp.Q, dtype=np.int32)
        self.code = set(int(c) for c in code_idx)
        _build_cover(np.array(sorted(self.code), dtype=np.int64), sp.qs, sp.radix, self.r, self.cnt)

    def add(self, idx):
        self.code.add(int(idx)); _ball_add(int(idx), self.sp.qs, self.sp.radix, self.r, self.cnt, 1)

    def remove(self, idx):
        self.code.discard(int(idx)); _ball_add(int(idx), self.sp.qs, self.sp.radix, self.r, self.cnt, -1)

    def free_all(self):
        return np.nonzero(self.cnt == 0)[0].astype(np.int64)

    def free_near(self, centres):
        cap = max(1000, len(centres) * self.sp.ball_size(self.r))
        cap = min(cap, 20_000_000)
        out = np.empty(cap, dtype=np.int64); m = 0
        for c in centres:
            m = _ball_free(int(c), self.sp.qs, self.sp.radix, self.r, self.cnt, out, m)
        return np.unique(out[:m])

    def words(self):
        return [self.sp.digits(i) for i in sorted(self.code)]


def augment(cov, max_cands=6000, iters=100000):
    F = cov.free_all()
    if len(F) == 0:
        return 0
    if len(F) > max_cands:
        # greedy in random order
        rng = np.random.default_rng(len(F))
        added = 0
        for i in rng.permutation(F):
            if cov.cnt[i] == 0:
                cov.add(i); added += 1
        return added
    C = _pair_compat(pack_indices(F, cov.sp.qs), cov.d)
    idx, w = max_clique(C, np.ones(len(F), dtype=np.int64), iters=iters, seeds=2)
    for i in idx:
        cov.add(F[i])
    return len(idx)


def refill(cov, seconds, rho=None, rng=None, plateau=0.5, max_cands=3000, log=None):
    """remove the codewords within distance rho of a random codeword, solve the hole exactly-ish, keep gains"""
    rng = rng or np.random.default_rng(0)
    sp, d = cov.sp, cov.d
    t0 = time.time(); moves = 0; gains = 0
    rho = rho if rho is not None else d
    while time.time() - t0 < seconds:
        code = np.array(sorted(cov.code), dtype=np.int64)
        pk = pack_indices(code, sp.qs)
        c0 = pk[rng.integers(len(code))]
        rr = int(rng.integers(max(1, rho - 1), rho + 2))
        dist = np.array([qdist(c0, p) for p in pk])
        R = code[dist <= rr]
        if len(R) == 0:
            continue
        for c in R:
            cov.remove(c)
        F = cov.free_near(R)
        if len(F) > max_cands:
            for c in R:
                cov.add(c)
            rho = max(1, rr - 1)
            continue
        C = _pair_compat(pack_indices(F, sp.qs), d)
        idx, w = max_clique(C, np.ones(len(F), dtype=np.int64), iters=20000, seeds=1, exact_limit=60)
        moves += 1
        newset = set(int(F[i]) for i in idx)
        if w > len(R) or (w == len(R) and newset != set(int(c) for c in R) and rng.random() < plateau):
            for c in newset:
                cov.add(c)
            if w > len(R):
                gains += w - len(R)
                if log:
                    log(f"  refill gain +{w - len(R)} -> {len(cov.code)} (removed {len(R)}, cands {len(F)})")
        else:
            for c in R:
                cov.add(c)
        # extra free words may have appeared elsewhere (never, but cheap to absorb)
    return gains, moves
