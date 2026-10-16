# A(32,10,10) >= 558

A binary constant-weight code of length 32, weight 10, minimum distance 10, with 558 words. The table value on Brouwer's page (read 2026-10-15 and 2026-10-16) is 500.
The code is maximal: an exhaustive check of all C(32,10) words found none that can be added.

- `code_558.json`: the words as lists of the 10 positions holding a 1.
- `code_558_matrix.txt`: the same code as a 0/1 matrix, one word per line.
- `verify_cw.py`: plain-Python check of weights, duplicates and all pairwise distances (`python verify_cw.py code_558.json`).

Found by a prescribed-automorphism search: the group Z_31 x| Z_5 acting on 31 coordinates with one fixed point (order 155), orbits chosen by a tabu clique
search on a sampled orbit pool. Also checked by an engine certifier and an all-pairs numpy check (maximum overlap 5, as required for d = 10).
