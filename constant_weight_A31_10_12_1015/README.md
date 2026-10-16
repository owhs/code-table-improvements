# A(31,10,12) >= 1015

A binary constant-weight code of length 31, weight 12, minimum distance 10, with 1015 words. Brouwer's page (read 2026-10-16) has 961, which
is the earlier code in `../constant_weight_A31_10_12/`.

- `code_1015.json`: the words as lists of the 12 positions holding a 1.
- `code_1015_matrix.txt`: the same code as a 0/1 matrix, one word per line.
- `verify_cw.py`: plain-Python check of weights, duplicates and all pairwise distances (`python verify_cw.py code_1015.json`).

## Structure
The code is a union of 5 orbits of 203 words under the group Z_29 x| Z_7 (x -> 7^k x + b on coordinates 0-28), with coordinates 29 and 30
fixed. By the fixed coordinates: 406 words contain neither, 203 contain only 29, 203 contain only 30, and 203 contain both.
The code is maximal: an exhaustive check of all C(31,12) words found none that can be added.

## How it was found and checked
Prescribed-automorphism search with this group, orbits chosen by a tabu clique search on a randomly sampled orbit pool. Two runs with different
pool seeds on two machines gave the same code. Also checked by an engine certifier and an all-pairs numpy check (maximum overlap 7, as required
for d = 10).
