# Independent reproductions of two 5-ary values of Laaksonen and Ostergard (2019)

The values are **not** ours. A. Laaksonen and P. R. J. Ostergard, *New lower bounds on q-ary error-correcting codes*,
Cryptogr. Commun. 11(5) (2019) 881-889, Table 3, gives A_5(7,4) >= 257, A_5(7,5) >= 57, A_5(8,5) >= 257,
A_5(9,5) >= 857 and A_5(9,6) >= 157. Laaksonen's results page (https://www.cs.helsinki.fi/u/ahslaaks/codes/) lists the
same five values as additional unpublished results and has the codes for the last four. Brouwer's 5-ary page, read
2026-10-16, still shows 250, 53, 165, 725 and 135 for these cells.

Our search found codes of the same sizes for two of the cells, independently of their files:

| file | cell | words | minimum distance |
|---|---|---:|---:|
| `code_A5_8_5_257.txt` | A_5(8,5) | 257 | 5 |
| `code_A5_7_5_57.txt` | A_5(7,5) | 57 | 5 |

- The 257 code is the zero word plus 8 sub-codes, one for each coordinate that is zero, each an optimal A_4(7,5) = 32
  code. Any code built this way has at most 1 + 8 * 32 = 257 words.
- The 57 code was found by the same prescribed-automorphism search as our other codes.
- For A_5(7,4) our best is 254, below 257. We have not seen a code for A_5(7,4) >= 257.

`check_codes.py` checks both files with the Python standard library only. Laaksonen's four published codes were also
checked here with an all-pairs distance check (57, 257, 857 and 157 words; minimum distances 5, 5, 5 and 6).
