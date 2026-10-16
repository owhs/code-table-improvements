"""Check a binary code file (one word per line, 0/1): equal lengths, no duplicates, minimum Hamming distance. Standard library only.
    python check_binary.py A18_4_6144.txt 4
    python check_binary.py A17_4_3072.txt 4
    python check_binary.py A16_3_3072.txt 3"""
import sys

path, d = sys.argv[1], int(sys.argv[2])
words = [l.strip() for l in open(path) if l.strip()]
n = len(words[0])
assert all(len(w) == n and set(w) <= {"0", "1"} for w in words), "bad word"
x = [int(w, 2) for w in words]
assert len(set(x)) == len(x), "duplicate words"
md = min(bin(x[i] ^ x[j]).count("1") for i in range(len(x)) for j in range(i + 1, len(x)))
print(f"{path}: {len(x)} words of length {n}, minimum distance {md}: {'OK' if md >= d else 'FAIL'} (required {d})")
