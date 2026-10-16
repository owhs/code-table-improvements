"""Write the 961 codewords of code_961.json as a plain 0-1 matrix (one 31-character row per codeword)
and re-check it: 961 distinct rows, each of weight 12, pairwise Hamming distance >= 10."""
import json, os

here = os.path.dirname(os.path.abspath(__file__))
code = json.load(open(os.path.join(here, "code_961.json")))
n = code["n"]
words = code.get("codewords") or code.get("words")
rows = sorted("".join("1" if i in set(w) else "0" for i in range(n)) for w in words)

assert len(rows) == 961 and len(set(rows)) == 961
assert all(r.count("1") == 12 for r in rows)
ints = [int(r, 2) for r in rows]
dmin = min(bin(a ^ b).count("1") for i, a in enumerate(ints) for b in ints[i + 1:])
assert dmin >= 10, dmin

with open(os.path.join(here, "code_961_matrix.txt"), "w", newline="\n") as f:
    f.write("\n".join(rows) + "\n")
print("rows", len(rows), "weight 12, min distance", dmin)
