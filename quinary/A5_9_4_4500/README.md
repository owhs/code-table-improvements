# A_5(9,4) >= 4500 (2026-10-15)

A 5-ary code of length 9, minimum Hamming distance 4, with **4500 words**.

| | value | source |
|---|---:|---|
| previous lower bound | 4375 | Laaksonen & Ostergard, "New lower bounds on q-ary error-correcting codes", Cryptogr. Commun. 11(5) (2019) 881-889, Table 3 (row n = 9, d = 4: upper 9672, lower 4375); Brouwer's 5-ary page lists 4375 too |
| this code | **4500** | found 2026-10-15 |
| upper bound | 9672 | same table |

## Structure
Write a word as (u | v | z) with u, v in Z_5^4 and z in Z_5 (coordinates 0-3, 4-7, 8). `D_36.txt` lists 36 pairs (delta, z), delta in Z_5^4, z in Z_5. The code is

    C = { (u | u + delta | z) : (delta, z) in D,  sum(u) + sum(u + delta) = 1 (mod 5) }        36 * 125 = 4500 words,

i.e. 36 cosets of the linear [9,3] code { (a | a | 0) : a in Z_5^4, sum(a) = 0 } (125 words, minimum weight 4). This is the same construction as in `../A5_10_4_18750/` (generalised (u | u+v), one layer) with one extra free coordinate z so that the length is odd.

## Files and checks
- `code_4500.txt`: the 4500 words (digits 0-4). `D_36.txt`: the 36 pairs. `check_code.py`: standard library only; rebuilds the code from D_36.txt, checks it equals the file and checks all 10,124,250 pairs (about 10 s): minimum distance 4.
- Also checked with a numpy all-pairs check inside the search script and with the plain-Python verifier `search_code/tools/qary_verify.py`: OK.
- Live tables read 2026-10-14, 15 and 16: Brouwer's 5-ary page lists 4375; the 2019 Table 3 (PDF) has 4375.

## How it was found
The standard prescribed-automorphism approach (Laaksonen and Ostergard 2017, 2019; Braun, Kohnert and Wassermann; Kramer and Mesner) is the
background. The layered search with one free coordinate reproduces the published 4375 (35 layers); the same search with many random
relabellings of the vertices then found 36 layers (6 of 10 runs of 20 minutes). The layered family is a generalised (u | u+v) construction
(Plotkin; generalised concatenation); no new principle is claimed. The same family gives at most 225 for A_5(7,4) (proven by an exact
solve), below the published 257, and stays below the tables at ternary odd lengths.

## Status
- Not submitted anywhere. Literature after 2019 searched by web search only (no improvement of A_5(9,4) found); not exhaustive.
- 36 layers is the best found; 37 was not seen in 10 runs (not a proof of a maximum).
