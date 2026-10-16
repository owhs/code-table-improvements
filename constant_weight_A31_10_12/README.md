# A(31,10,12) ≥ 961

Superseded by A(31,10,12) ≥ 1015 (`../constant_weight_A31_10_12_1015/`). Brouwer's page has listed 961 (label "San") since 2026-10-14.

A binary constant-weight code with length 31, minimum distance 10, weight 12 and **961 codewords**. The previous best
known lower bound was 930 (A. E. Brouwer, *Bounds for binary constant weight codes*,
https://aeb.win.tue.nl/codes/Andw.html, label "g", checked 2026-10-13 and again 2026-10-14; also 930 in
S. Polak, arXiv:1703.05171, 2017). By complementation, A(31,10,19) ≥ 961 too.

## Construction
Positions are Z_31. G is the affine group {x ↦ a·x + b mod 31 : a ∈ H, b ∈ Z_31}, where H is the subgroup of order 6
of GF(31)^* (so |G| = 186). The code is the union of the G-orbits of these 11 base words:

| orbit size | base word (positions) |
|---|---|
| 186 | 0 2 3 4 6 8 9 12 15 17 20 24 |
| 186 | 0 1 3 4 5 6 11 13 15 17 19 22 |
| 93 | 0 1 2 5 7 9 13 15 17 20 21 22 |
| 93 | 0 2 3 6 7 8 10 11 12 15 16 18 |
| 93 | 0 2 3 6 9 10 13 14 17 20 21 23 |
| 93 | 0 4 5 6 9 11 12 14 17 18 19 23 |
| 62 | 0 1 2 3 4 6 8 12 13 18 23 25 |
| 62 | 0 1 2 6 10 13 14 16 17 19 22 24 |
| 31 | 0 1 4 8 10 12 13 15 17 21 24 25 |
| 31 | 0 1 3 4 8 10 14 15 17 18 23 26 |
| 31 | 0 1 2 4 6 9 11 13 14 15 21 25 |

Total 2·186 + 4·93 + 2·62 + 3·31 = 961. Two weight-12 words are at distance ≥ 10 iff they share ≤ 7 positions.

## How it was found
All internally valid G-orbits (9,540 of them) were enumerated exactly by a pruned search, a compatibility graph
between orbits was built (an orbit pair is compatible when no two words share more than 7 positions), and an
**exact maximum-weight clique** was computed (networkx). The enumeration is complete (each orbit stored once, no
cap reached), so 961 is the largest code invariant under G. The same method proves that the larger affine subgroups
with |H| = 10 and |H| = 15 cannot exceed 930, and that unions of G-orbits and orbits of the |H| = 10 group together
(12,194 orbits) also cannot exceed 961.

## Verify
`code_961.json` lists all 961 words. `python -I verify_cw.py code_961.json` (plain Python, no dependencies) checks
weights, range, duplicates and every one of the 461,280 pairs, and prints the minimum distance (10).
