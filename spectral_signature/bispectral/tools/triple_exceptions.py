"""Which translation classes share a triple correlation? Test whether each exception is a mirror pair or a set/complement pair (w = n/2).
Result for (12,6), (14,7), (16,8), (18,9): 1, 1, 4, 7 pairs, all mirror pairs, no complement pairs.
(An earlier version tested only whether the complement of one member lay in the group, which is also true when a set is a
translate of its own complement; that test wrongly reported the pairs as complements.)"""
import itertools, sys, collections
import numpy as np

def trip(x, n):
    return tuple(int((x * np.roll(x, -a) * np.roll(x, -b)).sum()) for a in range(n) for b in range(n))

def canon(mask, n):
    full = (1 << n) - 1
    return min(((mask << k) | (mask >> (n - k))) & full for k in range(n))

for n, w in [(12, 6), (14, 7), (16, 8), (18, 9)]:
    groups = collections.defaultdict(set)
    for S in itertools.combinations(range(n), w):
        m = sum(1 << s for s in S)
        c = canon(m, n)
        if c in {c2 for g in groups.values() for c2 in g} and False:
            continue
        x = np.zeros(n, dtype=np.int64); x[list(S)] = 1
        groups[trip(x, n)].add(c)
    bad = [g for g in groups.values() if len(g) > 1]
    full = (1 << n) - 1
    mirror = lambda m: sum(1 << ((-i) % n) for i in range(n) if m >> i & 1)  # noqa: E731
    comp = sum(1 for g in bad if len(g) == 2 and canon(full ^ min(g), n) == max(g))
    mir = sum(1 for g in bad if len(g) == 2 and canon(mirror(min(g)), n) == max(g))
    print(f"n={n} w={w}: triple-correlation collisions between distinct translation classes: {len(bad)} groups; "
          f"mirror pairs: {mir}; set/complement pairs: {comp}; other: {len(bad) - mir - comp}", flush=True)
