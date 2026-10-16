# Spectral-signature codes (2026-10-14)

Tables, optimal codes and a checker for **codes for intensity-only (phase-blind) receivers**: sets of constant-weight binary words that a receiver which only measures
|Fourier|^2 (equivalently the cyclic autocorrelation, i.e. the multiset of cyclic differences of the support) can tell apart.

`brouwer.html` is the table page laid out like Brouwer's code tables (n down, w across, each entry linking to a file of codewords);
`index.html` is the longer page with the decoding experiments. Both are served at https://owhs.github.io/code-table-improvements/spectral_signature/.

## The problem in four lines
- Word x: binary, length n, weight w. Signature P_x(k) = |sum_j x_j e^{-2 pi i jk/n}|^2, k = 1..floor(n/2).
- Cyclic shifts, reflections and homometric (Z-related) partners have identical signatures; for w = n/2 so does the complement.
- Spectral distance d(x,y) = sum_k |P_x(k) - P_y(k)| / w^2. A spectral code C(n,w,theta) is a set of weight-w words with pairwise d >= theta.
- c(n,w) = number of distinct signatures (ceiling on any alphabet); S(n,w,theta) = largest spectral code (exact by integer programming where the table says so).

## What is in this folder
| file | content |
|---|---|
| index.html | the table page (c(n,w) for n <= 30; S(n,w,theta) for n <= 20 at theta = 0.05, 0.10, 0.20, 0.40; separation delta(n,w); matrices; demos) |
| codes/n<N>_w<W>_d<theta>.txt | one codeword per line as a 0/1 matrix row, for every computed cell and threshold |
| codes.json | the same codes as position lists, with size and proof status |
| table_classes.csv, table_maxcode.csv | the tables in machine-readable form |
| check_codes.py | dependency-free checker: weight, length, distinctness, signature separation (float and exact integer autocorrelation); `--recount` recounts c(n,w) for n <= 14 |
| hdc_demo.json | the bundling / readout experiment of section 7 of the page |
| related_work.md | prior art and what is and is not new |
| tools/ | the scripts that produced the data (if present) |

## Headline results (all produced by the scripts here)
1. **Class counts.** c(n,w) for all n <= 30, w <= n/2, from exact integer autocorrelation keys; cross-checked against an independent enumeration for every cell computed for the code tables (n <= 20), no mismatches. For n <= 24 these counts also follow from the enumeration of Z-related sets by Jedrzejewski and Johnson (2013); see `related_work.md`.
2. **Optimal codes.** S(n,w,theta) for n <= 20 at four thresholds: 239 of the 252 entries are proven optimal by integer programming; the other 13 are shown as lb-ub with a validated code for lb. (n,w) = (20,9) was not computed (3651 classes, above the 3500 limit used).
3. **Spectral codes against Hamming-distance codes under intensity-only detection.** Shift-closed families (including difference-set / projective-plane families that are optimal for ordinary Hamming codes) have a single signature: their error under
   intensity-only decoding is chance (0.91-0.94 at every noise level tested). Against Hamming-distance codes of the same size, a K-subset chosen to maximise the minimum spectral distance cut the symbol error rate by
   13x at sigma = 0.2 (n = 16, w = 6, K = 12: 0.009 vs 0.117) and 19x (n = 18, w = 6, K = 16: 0.006 vs 0.113), and still 3.3x at sigma = 0.3 (n = 16: 0.111 vs 0.361). Greedy Hamming codes of the same size contained two
   words with identical signatures in 5 of 5 draws at n = 13 (w = 4, K = 13) and 3 of 5 at n = 14 (w = 5, K = 20). Choosing K random distinct signatures (no distance optimisation) is clearly worse than optimising the
   minimum distance (e.g. 0.080 vs 0.009 at n = 16, sigma = 0.2). At sigma = 0.5 the advantage shrinks (0.474 vs 0.666 at n = 16).
4. **A link to hyperdimensional computing.** Bundles of shifted constant-weight codewords read out by correlation: an optical orthogonal code (all shifted-word overlaps <= 1) recovers every item exactly up to w - 1 = 3 bundled items
   (n = 127, w = 4) and 98% at 4 items, against 96% and 89% for random words (more at larger bundles); see the page.

## Verify
```
python check_codes.py            # every code file
python check_codes.py --recount  # plus brute-force recount of c(n,w) for n <= 14 against table_classes.csv
```

## Status
- Not submitted anywhere. Literature: see `related_work.md` (mostly abstracts and search summaries; the Keleti-Kolountzakis and Jedrzejewski-Johnson texts were read).
- The distance is a convenient L1 distance between noise-free signatures; it is not the optimal criterion for every noise model (the demos use additive Gaussian noise on the 0/1 positions before |FFT|^2).
- An earlier version of the demo took the first K words of a threshold code and showed a mixed result; the demos here use a proper maximum-minimum-distance K-subset, and the page says how the comparison was built.
- The mixed or worse behaviour at high noise (sigma >= 0.5) is real; nothing here claims an optimal code for a specific channel.
