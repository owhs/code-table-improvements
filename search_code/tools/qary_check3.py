"""Third, dependency-free checker for q-ary / mixed HIT files.
usage: python -I tools/qary_check3.py runs/qary_hits/HIT_*.json"""
import json, sys


def check(path, j=None):
    j = j or json.load(open(path))
    qs, d, words = j["qs"], j["d"], j["words"]
    n = len(qs)
    if any(len(w) != n for w in words):
        return "BAD length"
    for w in words:
        for x, q in zip(w, qs):
            if not isinstance(x, int) or not 0 <= x < q:
                return "BAD symbol range"
    if len({tuple(w) for w in words}) != len(words):
        return "BAD duplicates"
    if len(words) != j["size"]:
        return "BAD size field"
    worst = n
    for i in range(len(words)):
        a = words[i]
        for k in range(i + 1, len(words)):
            b = words[k]
            dist = 0
            for x, y in zip(a, b):
                if x != y:
                    dist += 1
            if dist < worst:
                worst = dist
                if worst < d:
                    return f"BAD distance {worst} < {d} between words {i},{k}"
    nb = sum(1 for q in qs if q == 2)
    return f"OK size={len(words)} d>={d} (min {worst}) binary={nb} ternary={sum(1 for q in qs if q == 3)} live_lower={j.get('live_lower')}"


if len(sys.argv) > 1:
    for p in sys.argv[1:]:
        print(p.split("/")[-1].split("\\")[-1], check(p))
else:  # published folder: check every code in codes.json next to this file
    import os
    allc = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "codes.json")))
    for name, j in allc.items():
        print(name, check(name, j))
