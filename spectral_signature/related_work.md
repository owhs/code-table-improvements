# Related work and what is (not) new

Sources below were found by web search on 2026-10-14. Where only an abstract or a search summary was read, it says so. Nothing here was read in full text.

## Closest prior art (the channel, not the code design)
- Wing Hong Wong, Yifei Lou, Stefano Marchesini, Tieyong Zeng, "One-Dimensional Phase Retrieval: Regularization, Box Relaxation and Uniqueness", arXiv:1904.10157 (2019, revised 2020).
  Abstract read: binary signals can be recovered from their Fourier magnitude uniquely up to trivial ambiguities under certain conditions. The conditions and the list of trivial
  ambiguities were NOT read. (A search summary mentioned shift, reflection and complement; for fixed weight the complement is a different word unless w = n/2.)
- Veit Elser and Stefan Eisebitt, "Uniqueness transition in noisy phase retrieval", arXiv:1009.4693 (2010). Abstract read: a mutual-information criterion for noisy phase
  retrieval from intensities. Not a code construction.
These treat recovery of ONE signal from intensities. The question here is different: which ALPHABET of weight-w words can a phase-blind receiver tell apart, how large can it be
at a given separation, and how much better is it than a Hamming-distance code.

## Non-coherent coding (a different invariance and metric)
- Grassmannian packings: codewords are subspaces separated in chordal distance, invariant to an unknown channel matrix; e.g. A. Ashikhmin and R. Calderbank,
  "Grassmannian packings from operator Reed-Muller codes", IEEE Trans. Inform. Theory 56(11), 2010, pp. 5689-5714 (listing seen in search results, not read).
  Our invariance is to cyclic shifts, reflections and Fourier phase of a discrete 0/1 word, with an L1 distance between power spectra.
- Constant-weight codes closed under cyclic shift for pulse position modulation, arXiv:1301.3369 ("Self-synchronizing pulse position modulation with error tolerance").
  These are exactly the shift-closed families that are BLIND to intensity-only detection (all shifts share one signature); see section 6 of the page.
- Optical orthogonal codes (OOCs): constant-weight binary sequences with bounded autocorrelation and cross-correlation (Chung-Salehi-Wei; Chung-Kumar 1990 "Optical orthogonal
  codes - new bounds and an optimal construction"; many later papers, e.g. arXiv:1804.07638, arXiv:1511.09289, arXiv:1702.06455 found in search, abstracts only). OOCs bound
  correlation values; spectral codes require DISTINCT signatures. Complementary, not the same constraint.

## Homometry, difference multisets, counting
- A. Rosenblatt and P. D. Seymour, "The structure of homometric sets", SIAM J. Algebraic Discrete Methods 3 (1982) 343-350 (via search snippets): A and B homometric iff
  A = U + V and B = U - V for suitable multisets; the structure theorem behind the merging of signatures.
- Erickson and Jones, "Homometric subsets of Z_n with cardinality 5: classification and enumeration", arXiv:2412.08997. Enumerates homometric 5-subsets; our c(n,5) can be compared.
- Lemke, Skiena, Smith ("Reconstructing sets from interpoint distances"), hexachordal theorem: on a cycle of length 2n every n-subset is homometric to its complement
  (via a search summary). This is why classes merge at w = n/2 in the tables.
- Music theory: Z-relation = homometry (Lewin; Mandereau, Ghisi, Amiot, "Z-relation and homometry in musical distributions" and Amiot's "The structure of Z-related sets",
  arXiv:1304.6608; both seen in search results).
- OEIS A005434 relates to the number of distinct autocorrelations of length-n words (aperiodic, as described by a LIRMM page); the cyclic, weight-resolved counts here are a different
  sequence. A search for the total counts of the class numbers returned no match on 2026-10-14.

## Literature check of 2026-10-16
Five searches (web search, about 50 queries in all), followed by reading the key sources directly. Findings:
- **Counts c(n,w).** No OEIS sequence contains our values. F. Jedrzejewski and T. Johnson, *The structure of Z-related sets*
  (arXiv:1304.6608; Mathematics and Computation in Music 2013) enumerate the Z-related (homometric) sets in Z_N for N up to 24, counted up to
  translation and inversion. Together with the number of bracelets this determines c(n,w) for n <= 24. Checked against our table:
  (12,6) 50 - 15 = 35, (16,6) 280 - 34 = 246, and the whole w = 4 column for n = 8 to 19. All agree. So c(n,w) is not new for n <= 24;
  the values for 25 <= n <= 30 were not found in print (the paper credits J. Wild with tables for larger N, not seen).
- **First merge.** The first homometric pair is at n = 8, w = 4: {0,1,2,5} and {0,1,3,4}.
- **Code sizes S(n,w,theta).** Not found. Searched: coding for the turnpike and beltway problems and polymer-based storage (Pattabiraman,
  Gabrys, Milenkovic, arXiv:1904.09280 and arXiv:2110.02352: linear prefix/suffix compositions, not cyclic autocorrelation);
  phaseless and intensity-only codebooks (Tasbihi and Kschischang, arXiv:2204.13806: direct detection, not binary constant-weight);
  spectral-amplitude-coding OCDMA (cross-correlation constraints, not equal spectra); cyclically permutable codes (Bitan and Etzion 1995),
  optical orthogonal codes (Chung, Salehi, Wei 1989) and comma-free codes. Those use Hamming or correlation constraints over shifts and
  keep homometric words apart; none tabulates codes under a spectral distance.
- **Name.** "Spectral-signature code" is not used in coding theory as far as found; "signature codes" is an existing term for another
  object (fingerprinting), so the full name is kept.
- **Bispectrum.** See `bispectral/README.md`: completeness of the 3-deck in Z_n is due to Grunbaum-Moore, Pebody and Keleti-Kolountzakis;
  bispectral codeword encoding for radio exists (Molchanov and Totsky 2008, search summaries only); separation-optimised codes and class
  counts were not found.

## What is new
- The packaging as a coding problem: weight-w words, spectral (|DFT|^2, equivalently autocorrelation) distance, exact optimal sizes S(n,w,theta) by integer programming, and files of optimal codes.
- The tables of S(n,w,theta) and the full-alphabet separation delta(n,w) (computed, checker in this folder). The table of c(n,w) is
  a convenience: for n <= 24 it follows from published enumerations (see above).
- The measured design rule: shift-closed (cyclic) constant-weight codes, including those from difference sets, are useless for intensity-only detection, while a spectral code of the
  same size has a lower error rate than a Hamming-distance code at moderate noise.
Not claimed: any theorem, any asymptotics, any new bound beyond the computed cells, or novelty relative to literature that was not found by the searches above.
