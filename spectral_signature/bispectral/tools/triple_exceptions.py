"""Which translation classes share a triple correlation? Test: are the exceptions exactly set/complement pairs (w = n/2)?"""
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
    comp = sum(1 for g in bad if len(g) == 2 and canon(full ^ next(iter(g)), n) in g)
    print(f"n={n} w={w}: triple-correlation collisions between distinct translation classes: {len(bad)} groups; of these set/complement pairs: {comp}; other: {len(bad)-comp}", flush=True)
