"""Independent verifier for a constant-weight code, written separately from the engine (plain Python, no numpy, no
engine imports). Usage: python -I tools/verify_cw.py <code.json>
Checks: every word has exactly w distinct positions in 0..n-1, no duplicates, and every pair of words is at Hamming
distance >= d (computed as the size of the symmetric difference)."""
import itertools
import json
import sys


def main():
    c = json.load(open(sys.argv[1]))
    n, d, w = c["n"], c["d"], c["w"]
    words = [frozenset(x) for x in c["words"]]
    assert all(len(x) == w for x in words), "a word has the wrong weight"
    assert all(0 <= p < n for x in words for p in x), "a position is out of range"
    assert len(set(words)) == len(words), "duplicate words"
    worst = min(len(a ^ b) for a, b in itertools.combinations(words, 2))
    print(f"n={n} d={d} w={w}: {len(words)} words, minimum pairwise distance {worst} -> "
          f"{'VALID' if worst >= d else 'INVALID'}")


if __name__ == "__main__":
    main()
