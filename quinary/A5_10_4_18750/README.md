# A_5(10,4) >= 18750 (2026-10-15)

A 5-ary code of length 10, minimum Hamming distance 4, with **18750 words**.

| | value | source |
|---|---:|---|
| previous lower bound | 17500 | Laaksonen and Ostergard, "New lower bounds on q-ary error-correcting codes", Cryptogr. Commun. 11(5) (2019) 881-889, Table 3; also Brouwer's 5-ary page (read 2026-10-16) and Laaksonen's results page |
| this code | **18750** | found 2026-10-15 |
| upper bound | 44642 | same table |

## Structure (a short proof of validity)
Write a word as (u | v) with u, v in Z_5^5 (coordinates 0-4 are u, 5-9 are v). Let D be the 30 vectors of `D_30.txt`. The code is

    C = { (u | v) : v - u in D,  sum(u) + sum(v) = 1 (mod 5) }        30 * 625 = 18750 words.

Equivalently, C is a union of 30 cosets of the linear [10,4] code K' = { (a | a) : a in Z_5^5, sum(a) = 0 }, which has minimum weight 4.
For two words with differences d1, d2 in D, put w = u - u', e = d1 - d2; the distance is wt(w) + wt(w + e), and the layer rule forces
sum(w) = 2 sum(e) (mod 5); every w of that sum class occurs. So the minimum distance is a minimum over the 30*31/2 pairs of D and the 625
vectors of one sum class: about 1.45 million evaluations. `check_code.py` does exactly that, and checks that `code_18750.txt` equals the
structure.

## Files
- `code_18750.txt`: the 18750 codewords, one per line, digits 0-4.
- `D_30.txt`: the 30 vectors of D.
- `check_code.py`: verifier, standard library only (about 10 s for the structural proof; `--brute` adds an all-pairs numpy check).

## Verification
- `check_code.py` (structural proof): minimum distance 4, file equals structure.
- `check_code.py --brute`: all 175,771,875 pairs, minimum distance 4 (numpy, per-coordinate comparison).
- A separate rebuild of the code from D (no code shared with the search): all pairs, minimum distance 4, 18750 distinct words.
- `search_code/tools/qary_verify.py` (plain Python, no numpy): 18750 words, minimum distance 4.

## How it was found
The method is the standard prescribed-automorphism search (Laaksonen and Ostergard 2017, 2019; Braun, Kohnert and Wassermann 2005;
Kramer and Mesner 1976). A first code with 18125 words came out of such a search. It turned out to be a union of cosets of a doubled code
(a | a), a generalised (u | u + v) construction (Plotkin; generalised concatenation). Restricting to one layer of that structure turns the
problem into an independent-set problem on the 3125 vectors v - u in Z_5^5, which a local search solves in minutes; one solution has 30
vectors. No new principle is claimed.

The published 17500 code (file from Laaksonen's results page) has the same form: it is a union of 28 cosets of a 625-word linear code of
the form (a | a), with coordinates paired as (i, i+5). Section 2 of the 2019 paper limits Cliquer to graphs of at most 5000 vertices; the
coset graph for this family has 15625 vertices. Whether that is why 30 cosets were not reached there is not known.

## Status
- Not submitted anywhere. Literature after 2019 searched by web search only (no improvement of A_5(10,4) found); not exhaustive.
- 30 vectors is the best found, not a proven maximum for this family.
