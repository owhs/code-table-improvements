# Bispectral-signature codes vs power-spectrum codes (extension of the spectral-signature table in the parent folder)

Constant-weight binary words, length n, weight w. Three receivers, each with its own invariances. Nothing in the parent folder was edited.

| receiver | sees | invariant to shift | invariant to mirror | homometric partners | complement at w=n/2 |
|---|---|---|---|---|---|
| **P** (baseline) | power spectrum, `d_P = sum_{k=1..n/2} abs(P_x(k)-P_y(k)) / w^2` | yes | yes | merged (same signature) | merged (hexachordal theorem) |
| **Tm** | triple correlation `T(a,b)=sum_t x_t x_{t+a} x_{t+b}`, signature merged with its mirror image, `d = min(d_T(x,y), d_T(x,mirror y))`, `d_T = sum_{a,b} abs(dT) / w^3` | yes | yes (by construction) | separated | separated |
| **T** | the same triple correlation, no mirror merging | yes | **no** (splits most mirror pairs; see below) | separated | separated |

Tm is the like-for-like comparison with P (both are blind to shift and reflection). T is shown to make the size of the mirror-splitting effect visible.

## Number of classes (normalisation-free)

| (n,w) | P | Tm | T | Tm/P | T/Tm (mirror splitting) |
|---|---|---|---|---|---|
| (12,6) | 35 | 50 | 80 | 1.43 | 1.60 |
| (14,7) | 85 | 133 | 246 | 1.56 | 1.85 |
| (16,6) | 246 | 280 | 504 | 1.14 | 1.80 |
| (16,8) | 254 | 440 | 810 | 1.73 | 1.84 |
| (17,8) | 698 | 750 | 1430 | 1.07 | 1.91 |
| (18,9) | 701 | 1387 | 2704 | 1.98 | 1.95 |

**Where the extra classes come from:** roughly a factor 1.6 to 1.95 of T's gain over P is mirror splitting (T/Tm), and the remaining factor 1.07 to 1.98 (Tm/P) is homometric partners becoming distinguishable. For n=16, w=8, T has 3.19 times as many classes as P; 1.84 of that is mirror splitting and 1.73 is resolving homometric partners. The T column alone overstates the benefit and should not be quoted as a like-for-like gain.

## Maximum code size S(theta) (exact MILP, HiGHS, 60 s limit; "opt" = proven optimal, `<= u` = best found with upper bound u)

| (n,w) | theta | P | Tm | T |
|---|---|---|---|---|
| (12,6) | 0.20 | 29 opt | 46 opt | 65 opt |
| (12,6) | 0.40 | 11 opt | 15 opt | 16 opt |
| (14,7) | 0.20 | 85 opt | 121 opt | 152 opt |
| (14,7) | 0.40 | 20 opt | 20 opt | 25 opt |
| (16,6) | 0.20 | 246 opt | 277 opt | 497 opt |
| (16,6) | 0.40 | 108 opt | 116 opt | 177 opt |
| (16,8) | 0.20 | 172 opt | 267 opt | 427 opt |
| (16,8) | 0.40 | 36 opt | 35 opt | 42 (<= 46) |
| (17,8) | 0.20 | 698 opt | **598 opt** | 1006 opt |
| (17,8) | 0.40 | 46 (<= 61) | 51 (<= 67) | unresolved |
| (18,9) | 0.20 | 300 opt | 592 opt | 960 opt |
| (18,9) | 0.40 | 45 (<= 57) | unresolved | unresolved |

Thresholds 0.05 and 0.10 are in `results/*.txt` (there almost every class is separable, so the table is about class counts).

**Cautions that matter:**
- At theta 0.40 for (17,8) and (18,9) the solver returned no incumbent within 60 s for T and Tm (printed as "0 <= N" in the raw files). That means unresolved, not zero.
- The two distances use different normalisations (`/w^2` vs `/w^3`), so the same theta does not mean the same noise margin. The class counts above are the normalisation-free comparison. At (17,8), Tm at theta 0.20 (598) is below P (698), because the triple-correlation distances are smaller. A bigger signature space does not automatically mean a bigger code at a fixed theta.
- At large theta the gain over P largely vanishes (for example 36 / 35 / 42 at (16,8), theta 0.40).

## Noise demo (K=12 words from a randomised farthest-point search, nearest-signature decoding on the empirical signature of x + sigma*noise)

Error rates are about the same for all three receivers at every size tried. At sigma 0.3: (12,6) 0.225 / 0.220 / 0.185 (P / Tm / T), (14,7) 0.117 / 0.152 / 0.102, (16,8) 0.120 / 0.090 / 0.105, (16,6) 0.110 / 0.052 / 0.048, (17,8) 0.045 / 0.060 / 0.043, (18,9) 0.025 / 0.050 / 0.040. 400 trials per cell, so each entry is good to about +-0.03 to +-0.05; the only difference clearly outside that is (16,6) at sigma 0.3. The extra classes do not buy noise robustness at these sizes.

## Noise test with averaging and coloured noise (`tools/bispectral_noise.py`, `results/bispnoise_*`)

Setup: K=12 words per receiver (farthest-point in that receiver's own distance); each observation is the word at a random shift plus zero-mean Gaussian noise; M=1 or 8 snapshots averaged; noise either white or with a random unknown smooth PSD (mean PSD 1); 'blind' = raw signature, 'oracle' = true noise floor / bias subtracted (E|Y(k)|^2 = P_x(k) + S(k); E T_y(a,b) = T_x(a,b) + w(R(a)+R(b)+R(b-a)), the same for every codeword). 500 trials per cell, so each error rate is good to about +-0.02. Error rate (P vs Tm), M=8, sigma = 0.6 and 1.0, blind mode:

| (n,w) | white 0.6 | white 1.0 | coloured 0.6 | coloured 1.0 |
|---|---|---|---|---|
| (12,6) | 0.414 vs 0.126 | 0.842 vs 0.556 | 0.096 vs 0.380 | 0.440 vs 0.786 |
| (16,6) | 0.174 vs 0.042 | 0.824 vs 0.428 | 0.132 vs 0.310 | 0.430 vs 0.720 |
| (16,8) | 0.178 vs 0.056 | 0.804 vs 0.392 | 0.046 vs 0.266 | 0.456 vs 0.666 |
| (17,8) | 0.082 vs 0.022 | 0.832 vs 0.318 | 0.054 vs 0.240 | 0.642 vs 0.702 |
| (18,9) | 0.092 vs 0.006 | 0.716 vs 0.276 | 0.012 vs 0.198 | 0.206 vs 0.674 |

- **White noise: the triple-correlation receiver (Tm) wins clearly** in all five cells once snapshots are averaged. Oracle bias removal changes little for either receiver (white P blind vs oracle: for example (16,8), sigma 0.6: 0.178 vs 0.144), so the gain is lower estimator variance, not bias.
- **Coloured noise: the power-spectrum receiver (P) wins** in blind mode in every cell, and in oracle mode too, with a smaller margin (oracle, M=8: P wins every cell at sigma 0.6, for example (16,8) 0.048 vs 0.098; at sigma 1.0 P wins four cells, nearly tied at (16,8) 0.446 vs 0.472, and Tm wins only (17,8), 0.624 vs 0.512). The hypothesis that a bispectral receiver is the better choice when the noise colour is unknown was **not** supported.
- **Likely reason, not tested:** my random PSD puts most noise power at low frequencies, leaving high-k power-spectrum bins nearly clean, while the triple correlation mixes all bins. Other noise colours or non-Gaussian noise could change this.
- **Single snapshot (M=1):** no consistent winner in the pilot.
- **Caveats:** the two receivers use different code sets and distance scales; the gain in white noise could partly reflect Tm having more classes and more signature entries, not a property of third-order statistics as such. One noise family, one K, no theory behind the sign of the effect.

## Prior art (so nothing here is claimed as new beyond the capacity view)

- Completeness of the triple correlation for subsets of finite cyclic groups: Grunbaum and Moore (3-deck), and Keleti and Kolountzakis, arXiv:math/0603415 (which n make the 3-deck determine every set, plus the random-set result). I read the arXiv text (introduction and results summary): the 3-deck (triple correlation) determines every subset of Z_n up to translation for n a power of an odd prime, a product of two odd primes (Pebody, Gruenbaum-Moore), p^2 q or pqr with distinct odd primes (their Theorem 2.1); it fails for n = 2k with k >= 6 (counterexamples given) and for n = pqrd with p, q distinct primes and r, d > 1. That matches our enumeration: n = 15 and 17 fully resolved, even n = 12, 14, 16, 18 show collisions. Their Theorem 2.22 (read in the full text) gives, for n = 2k with k >= 6, the pair E = {0, 3, 4, ..., k-1, k+1, k+2} and F = {0, 1, 3, ..., k-1, k+2}: both of size n/2, same 3-deck, and F a translate of -E (a mirror pair).
- Bispectrum as a complete invariant for generic signals: Kakarala (thesis 1992; arXiv:0902.0196); the skew spectrum (arXiv:0712.4259). Kondor (arXiv:cs/0701127) proposed related invariants for images; the arXiv comments say its uniqueness claim was dropped, so it is not cited for completeness.
- Bispectrum samples have been used to encode codewords for radio links (Molchanov and Totsky 2008; Naumenko et al. 2015; seen in search summaries only, not read). Those works do not optimise a minimum bispectral distance or count classes.
- Power spectrum at w = n/2: a set and its complement have the same power spectrum (hexachordal theorem; Babbitt; Kassler 1961, Regener 1974, Wilcox 1983; Lemke-Skiena-Smith). This is a property of P, not of T.
- Our own checks, exhaustive for (12,6), (14,7), (16,8), (18,9): the triple correlation separates every translation class except 1, 1, 4 and 7 pairs, and every one of these pairs is a set and its mirror image (no set/complement pair). At (12,6) and (14,7) the pair is exactly the Keleti-Kolountzakis pair above. An earlier version of this README called these pairs set/complement pairs; that came from a faulty test in `tools/triple_exceptions.py`, now fixed.

The only part not found in the literature search is the code-capacity question (maximum signature-separated code size at a distance threshold), and a thin search can miss things.

## Files
- `tools/bispectral_codes.py n w` - class counts, S(theta) table, noise demo. `tools/triple_vs_power.py a b` and `tools/triple_exceptions.py` - the class-count and complement-collision checks.
- `results/bispectral_n_w.txt|json` - raw outputs. (12,6) has only the json.
- Scripts expect to write to a `research/` folder in the current directory (create it first).
