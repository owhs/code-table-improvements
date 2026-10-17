# A_5(9,5) >= 893

5-ary code, length 9, minimum distance 5, 893 words (`code_893.txt`, one word per line, symbols 0..4).
Brouwer's table lists 725 (Laaksonen and Ostergard, ICMCTA 2017). Their 2019 paper reports 857, and this code has 36 words more than that.

Checks: `python check_code.py` (direct all-pairs comparison) and `python check_code2.py` (one-hot matrix product). Both need numpy and print the word count and the minimum distance.

How it was found: start from the 857-word code of Laaksonen and Ostergard (2019), and run a large-neighbourhood search over single words (the trivial symmetry group): remove a random set of words near a
random codeword, list the words that become free, and solve the independent-set problem on them with penalty local search, keeping the result if it is no smaller. This is a standard method, not a new one.
It is the same procedure that raised A_3(15,5) from 7812 to 7929. The 857 code has an automorphism group of order 216 and only 13 orbits; at that granularity the search found nothing, and at word level
(857 pieces) removing about 171 words per move gained 36. Two continuations from 893 (removing 133 to 532 words per move, 60 moves per size, none capped) found no further gain, so 893 is a plateau of this search, not a proven optimum.

The 857 code is Laaksonen and Ostergard's published code, used as the starting point. Tool: `tools/orbit_lns.py` in the reframe harness.
