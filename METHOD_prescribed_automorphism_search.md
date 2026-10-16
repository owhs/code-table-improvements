# Method: prescribed-automorphism ("orbit + clique") search

## The search
A symmetry group G acting on the words is fixed (coordinate permutations, symbol scalings and shifts, translations by a linear code).
Whole orbits of G are then chosen instead of single words:
1. **Enumerate the orbits** of G exactly and discard orbits that are not valid on their own (two words of the orbit closer than d).
2. **Build the conflict graph** on the remaining orbits: two orbits conflict if some pair of words, one from each, is closer than d.
   Checking one representative against the other orbit is enough, because G acts by isometries.
3. **Choose a maximum-weight clique** of compatible orbits (weight = orbit size): exact when the graph is small, tabu search otherwise.
   For example, a 5-ary code of length 10 has 9.8 million words, but under a group of order 1250 there are about 8 thousand orbits.
4. **Polish at the word level**: add every word still compatible with the code, then remove and refill the words around a random codeword.
5. **Repeat with many groups.** Groups are sampled from cyclic, dihedral and affine maps x -> ax+b on Z_p, multiplier subgroups,
   PSL(2,p), AGL(1,2^k), block-cyclic groups with fixed points (constant weight), per-coordinate symbol maps, and translations by a code K
   (so that the code is a union of cosets of K).

Every code was checked by programs that do not share code with the search, and compared with the live table before being listed.

## Prior work
The method is standard; it is not claimed as new.
- E. S. Kramer and D. M. Mesner, *t-designs on hypergraphs*, Discrete Math. 15(3) (1976) 263-296: designs as unions of orbits of a
  prescribed group, chosen by solving a 0-1 system.
- M. Braun, A. Kohnert, A. Wassermann, *Optimal linear codes from matrix groups*, IEEE Trans. Inform. Theory 51(12) (2005) 4247-4251;
  A. Kohnert, *Constructing two-weight codes with prescribed groups of automorphisms*, Discrete Appl. Math. 155 (2007) 1451-1457.
- P. Kaski and P. R. J. Ostergard, *Classification Algorithms for Codes and Designs*, Springer 2006 (prescribed automorphism groups).
- Exact clique search: Cliquer (S. Niskanen and P. R. J. Ostergard); P. R. J. Ostergard, *A fast algorithm for the maximum clique problem*,
  Discrete Appl. Math. 120 (2002), doi:10.1016/S0166-218X(01)00290-6.
- A. Laaksonen and P. R. J. Ostergard, *Constructing error-correcting binary codes using transitive permutation groups*,
  Discrete Appl. Math. 233 (2017) 65-70 (arXiv:1604.06022): orbits chosen as a maximum-weight clique, for binary codes.
- A. Laaksonen and P. R. J. Ostergard, *New lower bounds on error-correcting ternary, quaternary and quinary codes*, Proc. 5th International
  Castle Meeting on Coding Theory and Applications (2017) 228-237, and *New lower bounds on q-ary error-correcting codes*,
  Cryptogr. Commun. 11(5) (2019) 881-889: the same pipeline for q-ary codes. Our 5-ary results are directly comparable with their Table 3.
- The tables: A. E. Brouwer, J. B. Shearer, N. J. A. Sloane, W. D. Smith, *A new table of constant weight codes*, IEEE Trans. Inform.
  Theory 36 (1990) 1334-1380; W. D. Smith, S. Hughes, S. Perkins, *A new table of constant weight codes of length greater than 28*,
  Electron. J. Combin. 13 (2006) #A2; A. E. Brouwer, H. O. Hamalainen, P. R. J. Ostergard, N. J. A. Sloane, *Bounds on mixed binary/ternary
  codes*, IEEE Trans. Inform. Theory 44 (1998) 140-161; G. T. Bogdanova and P. R. J. Ostergard, *Bounds on codes over an alphabet of five
  elements*, Discrete Math. 240 (2001) 13-19; and the live pages at https://aeb.win.tue.nl/codes/.

## What differs here
These are practical differences, not theoretical advances.
1. **Group sampling.** Groups are drawn automatically from a wide library instead of being chosen by hand, with many short runs instead of a
   few long ones.
2. **One program for several alphabets**: binary, ternary, quaternary, 5-ary and mixed binary/ternary words.
3. **Structure read off the results.** The 5-ary records are unions of cosets of a linear code of the doubled form (a | a):
   18750 = 30 cosets of a [10,4] code and 4500 = 36 cosets of a [9,3] code. The published codes for these cells have the same structure:
   Laaksonen's 17500 is 28 cosets of a 625-word code of the form (a | a) (checked on the code file from his page), and 4375 is 35 cosets.

## Limits
- The method is only as good as its groups. In 15-minute single-core runs it did not reach several published values,
  for example A_5(9,5) (625 against 857) and A_5(9,6) (125 against 157).
- For A_5(7,4) our best is 254, below the 257 of Laaksonen and Ostergard (2019).
- Literature after the table sources was searched by web search only.
