# Improved lower bounds for ternary and mixed binary/ternary codes (2026-10-14)

Codes beating the lower bounds on A. E. Brouwer's tables (https://aeb.win.tue.nl/codes/ternary-1.html and
https://aeb.win.tue.nl/codes/23codes.html, fetched and compared 2026-10-14, and again on 2026-10-16 with no change in
these cells). N(m,n,d) is the maximum size of a code with
m binary and n ternary coordinates and minimum Hamming distance d. In the files the binary coordinates come first.

| cell | table lower bound | table upper bound | new size | gain | table page |
|---|---:|---:|---:|---:|---|
| A_3(14,8) | 108 | 237 | **123** | +15 | ternary-1.html |
| N(1,13,8) | 81 | 179 | **86** | +5 | 23codes.html |
| N(2,9,7) | 18 | 24 | **19** | +1 | 23codes.html |
| N(2,11,7) | 75 | 144 | **77** | +2 | 23codes.html |
| N(3,9,7) | 28 | 48 | **30** | +2 | 23codes.html |
| N(3,10,5) | 729 | 1519 | **736** | +7 | 23codes.html |
| N(3,10,8) | 24 | 36 | **26** | +2 | 23codes.html |
| N(5,8,5) | 486 | 791 | **513** | +27 | 23codes.html |
| N(5,8,7) | 36 | 72 | **39** | +3 | 23codes.html |
| N(5,9,6) | 486 | 744 | **513** | +27 | 23codes.html |
| N(5,9,7) | 69 | 174 | **73** | +4 | 23codes.html |
| N(5,9,8) | 30 | 52 | **32** | +2 | 23codes.html |
| N(6,8,8) | 28 | 44 | **32** | +4 | 23codes.html |
| N(7,4,5) | 44 | 56 | **45** | +1 | 23codes.html |
| N(7,5,5) | 99 | 144 | **104** | +5 | 23codes.html |
| N(7,6,7) | 24 | 41 | **27** | +3 | 23codes.html |
| N(8,4,5) | 84 | 112 | **85** | +1 | 23codes.html |
| N(9,4,4) | 548 | 768 | **573** | +25 | 23codes.html |
| N(9,4,5) | 136 | 216 | **144** | +8 | 23codes.html |
| N(10,3,3) | 1152 | 1492 | **1176** | +24 | 23codes.html |
| N(10,4,4) | 1050 | 1480 | **1104** | +54 | 23codes.html |
| N(10,4,5) | 234 | 390 | **235** | +1 | 23codes.html |
| N(10,4,6) | 88 | 144 | **96** | +8 | 23codes.html |
| N(12,2,3) | 1536 | 2011 | **1552** | +16 | 23codes.html |
| N(6,8,6) | 324 | 527 | **342** | +18 | 23codes.html (derived) |
| N(7,7,6) | 216 | 375 | **228** | +12 | 23codes.html (derived) |

The two derived codes come from N(5,9,6) >= 513: keep the words whose ternary coordinate 5 takes its two most frequent
symbols and relabel that coordinate as binary (342 words); doing the same on coordinate 6 as well gives 228 words. Hamming
distances are unchanged by this step. Both were found on 2026-10-16 by applying shortening, puncturing and relabelling
to all codes here and comparing with the live table.

## Method
Each code is a union of orbits of a group of coordinate isometries (permutations, symbol scalings and shifts), found by
orbit enumeration and a weighted clique search on the orbit conflict graph, then extended by word-level augment and refill.

## Files and verification
- codes/<cell>.txt: one codeword per line, one digit per coordinate; codes.json: all codes with alphabets and constructions.
- check_codes.py: dependency-free checker (symbol ranges, distinctness, pairwise Hamming distance). `python check_codes.py`.
- Every code was also checked by `search_code/tools/qary_verify.py` and `search_code/tools/qary_check3.py`, which do not share code with the search.

Status: not yet submitted to the maintainer. Literature beyond the tables has not been checked.
