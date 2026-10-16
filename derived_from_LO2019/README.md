# Mixed binary/ternary bounds that follow from Laaksonen and Ostergard's ternary codes

These five values are **consequences of published codes**, not new searches. The ternary codes A_3(13,4) >= 13122 and A_3(15,4) >= 83106
are due to A. Laaksonen and P. R. J. Ostergard (Castle Meeting 2017; Cryptogr. Commun. 11(5) (2019)); the code files are from Laaksonen's
results page, https://www.cs.helsinki.fi/u/ahslaaks/codes/. Brouwer's ternary page has these values, but his mixed binary/ternary page
(https://aeb.win.tue.nl/codes/23codes.html, read 2026-10-16) does not yet have what follows from them.

The step used: restrict one ternary coordinate to its two most frequent symbols and relabel it as binary. Hamming distances are unchanged and
at least 2/3 of the words remain.

| cell | Brouwer's page (lower-upper) | this code | construction |
|---|---|---:|---|
| N(1,12,4) | 8019-13531 | **8748** | A_3(13,4) = 13122, one ternary coordinate made binary |
| N(1,13,4) | 16767-37714 | **18711** | A_3(15,4) = 83106 shortened at one coordinate (largest fibre), then one ternary coordinate made binary |
| N(2,11,4) | 5589-9450 | **5832** | A_3(13,4) = 13122, two coordinates made binary |
| N(3,10,4) | 3726-6581 | **3888** | A_3(13,4) = 13122, three coordinates made binary |
| N(4,9,4) | 2484-4590 | **2592** | A_3(13,4) = 13122, four coordinates made binary |

The check that found them: the live mixed table was tested for closure under weaker distance, puncturing, shortening, binary-as-ternary and
this ternary-to-binary step (lower bound ceil(2N/3)). These five cells are the only ones that move; the quaternary and 5-ary tables are
closed under shortening and puncturing. The codes here were then built explicitly from the code files (N(1,13,4) gives 18711, more than
the guaranteed 18468).

## Files and checks
- `codes/<cell>.txt`: one word per line, binary coordinates first. `codes.json`: the same codes with the construction.
- `check_codes.py`: standard-library checker (`python check_codes.py`). All five also pass `search_code/tools/qary_verify.py` and
  `search_code/tools/qary_check3.py`.
