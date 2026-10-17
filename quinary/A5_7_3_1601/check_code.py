"""Check the 5-ary code in this folder (standard library only): distinct words, alphabet 0-4, length, minimum distance.
    python check_codes.py
"""
import itertools
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CODES = (("code_1601.txt", 7, 3, 1601),)

ok_all = True
for name, n, d, size in CODES:
    words = [line.strip() for line in open(os.path.join(HERE, name)) if line.strip()]
    ok = (len(words) == size and len(set(words)) == size
          and all(len(w) == n and set(w) <= set("01234") for w in words))
    md = min(sum(a != b for a, b in zip(x, y)) for x, y in itertools.combinations(words, 2))
    ok = ok and md >= d
    ok_all &= ok
    print(f"{name}: {len(words)} words, length {n}, minimum distance {md} (need {d}) -> {'VALID' if ok else 'INVALID'}")
print("RESULT:", "all valid" if ok_all else "FAILED")
