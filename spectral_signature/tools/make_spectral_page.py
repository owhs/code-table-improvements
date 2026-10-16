"""Generate publish/Spectral_codes_2026: code files, matrices, CSV tables and index.html from _build_data.json."""
import json, os, html

OUT = r"C:\projects\2026\reframe\publish\Spectral_codes_2026"
D = json.load(open(os.path.join(OUT, "_build_data.json")))
cells, counts, demos, THETAS = D["cells"], D["counts"], D["demos"], D["thetas"]
os.makedirs(os.path.join(OUT, "codes"), exist_ok=True)

# ------------------------------------------------------------------ code files + codes.json
allcodes = {}
for key, c in cells.items():
    if c.get("skipped"):
        continue
    n, w = c["n"], c["w"]
    for th, cd in c["codes"].items():
        name = f"n{n}_w{w}_d{float(th):.2f}"
        rows = []
        for pos in cd["words"]:
            v = ["0"] * n
            for p in pos:
                v[p] = "1"
            rows.append("".join(v))
        status = "optimal" if cd["optimal"] else f"best found; upper bound {cd['upper']}"
        header = [f"# spectral-signature code, n={n} w={w}, normalised spectral distance >= {th}: {len(rows)} words ({status})",
                  "# one codeword per line, bit i = position i; words are class representatives",
                  "# distance d(x,y) = sum_{k=1}^{floor(n/2)} | |X_x(k)|^2 - |X_y(k)|^2 | / w^2"]
        with open(os.path.join(OUT, "codes", name + ".txt"), "w") as f:
            f.write("\n".join(header + rows) + "\n")
        allcodes[name] = dict(n=n, w=w, theta=float(th), size=len(rows), optimal=cd["optimal"], upper=cd["upper"], words=cd["words"])
json.dump(allcodes, open(os.path.join(OUT, "codes.json"), "w"))

# ------------------------------------------------------------------ CSV
with open(os.path.join(OUT, "table_classes.csv"), "w") as f:
    f.write("n,w,bracelets,classes\n")
    for n in sorted(counts, key=int):
        for w in sorted(counts[n], key=int):
            b, c = counts[n][w]
            if 2 <= int(w) <= int(n) // 2:
                f.write(f"{n},{w},{b},{c}\n")
with open(os.path.join(OUT, "table_maxcode.csv"), "w") as f:
    f.write("n,w,classes,delta_min," + ",".join(f"size_d{t:.2f},optimal_d{t:.2f},upper_d{t:.2f}" for t in THETAS) + "\n")
    for key, c in cells.items():
        if c.get("skipped"):
            continue
        row = [c["n"], c["w"], c["classes"], round(c["delta_min"], 4) if c["delta_min"] is not None else ""]
        for t in THETAS:
            cd = c["codes"][str(t)]
            row += [cd["size"], int(cd["optimal"]), cd["upper"]]
        f.write(",".join(map(str, row)) + "\n")

# ------------------------------------------------------------------ HTML
def cell_class(n, w):
    r = counts.get(str(n), {}).get(str(w))
    if r is None:
        return ""
    b, c = r
    return f"{c}" if b == c else f'<span title="{b} bracelets, {b - c} merged by homometry">{c}*</span>'

ns = sorted(int(n) for n in counts if 4 <= int(n) <= 30)
t1 = ['<TABLE BORDER="1"><TR><TD>n \\ w</TD>' + "".join(f"<TD> {w} </TD>" for w in range(2, 16)) + "</TR>"]
for n in ns:
    t1.append(f"<TR><TD>{n}</TD>" + "".join(f"<TD>{cell_class(n, w) if w <= n // 2 else ''}</TD>" for w in range(2, 16)) + "</TR>")
t1.append("</TABLE>")

def bound_cell(cd):
    return f"{cd['size']}" if cd["optimal"] else f"{cd['size']}-{cd['upper']}"

maxtabs = []
for t in THETAS:
    rows = ['<TABLE BORDER="1"><TR><TD>n \\ w</TD>' + "".join(f"<TD> {w} </TD>" for w in range(3, 11)) + "</TR>"]
    for n in range(6, 21):
        tds = []
        for w in range(3, 11):
            c = cells.get(f"{n},{w}")
            if c is None or w > n // 2:
                tds.append("<TD></TD>")
            elif c.get("skipped"):
                tds.append("<TD>&mdash;</TD>")
            else:
                cd = c["codes"][str(t)]
                full = cd["size"] == c["classes"]
                label = bound_cell(cd)
                link = f'<A HREF="codes/n{n}_w{w}_d{t:.2f}.txt">{label}</A>'
                tds.append(f"<TD>{'<B>' + link + '</B>' if not full else link}</TD>")
        rows.append(f"<TR><TD>{n}</TD>" + "".join(tds) + "</TR>")
    rows.append("</TABLE>")
    maxtabs.append((t, "\n".join(rows)))

dm = ['<TABLE BORDER="1"><TR><TD>n \\ w</TD>' + "".join(f"<TD> {w} </TD>" for w in range(3, 11)) + "</TR>"]
for n in range(6, 21):
    tds = []
    for w in range(3, 11):
        c = cells.get(f"{n},{w}")
        tds.append("<TD></TD>" if c is None or w > n // 2 or c.get("skipped") or c["delta_min"] is None else f"<TD>{c['delta_min']:.3f}</TD>")
    dm.append(f"<TR><TD>{n}</TD>" + "".join(tds) + "</TR>")
dm.append("</TABLE>")

demo_html = []
DV2 = json.load(open(os.path.join(OUT, "_demo_v2.json")))
for d in DV2:
    coll = sum(1 for x in d["hamming_min_spectral"] if x < 1e-9)
    demo_html.append(f"<P>n={d['n']}, w={d['w']}, K={d['K']} words (chance error {d['chance']:.3f}; {d['classes']} signatures exist). "
                     f"<B>spectral-opt</B>: the K classes found by randomised farthest-point search with the largest minimum spectral distance ({d['spectral_min_d']:.3f}). "
                     f"<B>distinct</B>: K random classes (min distance {d['distinct_min_d']:.3f}). "
                     f"<B>Hamming</B>: greedy codes of maximal Hamming distance, 5 random draws (mean, [min-max] of error); "
                     f"{coll} of 5 draws contain two words with identical signatures. <B>cyclic</B>: K shifts of one word.</P>")
    demo_html.append('<TABLE BORDER="1"><TR><TD>noise sigma</TD><TD>spectral-opt</TD><TD>distinct</TD><TD>Hamming mean [min-max]</TD><TD>cyclic family</TD></TR>')
    for r in d["rows"]:
        demo_html.append(f"<TR><TD>{r['sigma']}</TD><TD><B>{r['spectral_opt']:.3f}</B></TD><TD>{r['distinct']:.3f}</TD>"
                         f"<TD>{r['hamming_mean']:.3f} [{r['hamming_min']:.3f}-{r['hamming_max']:.3f}]</TD><TD>{r['cyclic']:.3f}</TD></TR>")
    demo_html.append("</TABLE>")

HDC = json.load(open(os.path.join(OUT, "hdc_demo.json")))
hdc_rows = ['<TABLE BORDER="1"><TR><TD>items bundled T</TD><TD>OOC codebook (corr &le; 1)</TD><TD>random weight-4 words</TD><TD>Hamming code</TD></TR>']
for T, r in HDC.items():
    hdc_rows.append(f"<TR><TD>{T}</TD><TD><B>{r['OOC (corr<=1)']:.3f}</B></TD><TD>{r['random']:.3f}</TD><TD>{r['Hamming code']:.3f}</TD></TR>")
hdc_rows.append("</TABLE>")

# matrices
singer = [[1 if ((j - i) % 13) in (0, 1, 3, 9) else 0 for j in range(13)] for i in range(13)]
singer_txt = "\n".join("".join(map(str, r)) for r in singer)
ex = allcodes.get("n12_w6_d0.20")
ex_txt = "\n".join(l for l in open(os.path.join(OUT, "codes", "n12_w6_d0.20.txt")).read().splitlines() if not l.startswith("#")) if ex else ""

page = f"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<HTML>
<HEAD>
 <META http-equiv=Content-Type content="text/html; charset=utf-8">
 <TITLE>Spectral-signature codes: codes for intensity-only receivers</TITLE>
</HEAD>
<BODY>
<H2>Spectral-signature codes: tables of bounds</H2>
<P>Computed 2026-10-14 (exact where marked). Files: <A HREF="README.md">README</A>, <A HREF="codes.json">codes.json</A>, <A HREF="table_classes.csv">table_classes.csv</A>,
<A HREF="table_maxcode.csv">table_maxcode.csv</A>, <A HREF="check_codes.py">check_codes.py</A>, <A HREF="related_work.md">related work</A>.</P>

<H3>1. The problem</H3>
<P>A word x is a binary vector of length n and weight w. A receiver that can only measure intensities sees the power signature
P<sub>x</sub>(k) = |&sum;<sub>j</sub> x<sub>j</sub> e<sup>-2&pi;ijk/n</sup>|<sup>2</sup>, k = 1..&lfloor;n/2&rfloor;
(equivalently the cyclic autocorrelation, i.e. the multiset of cyclic differences, of the support). Cyclic shifts, reflections and homometric (Z-related)
partners have identical signatures; for w = n/2 a word and its complement do too. Define the normalised spectral distance
d(x,y) = &sum;<sub>k</sub> |P<sub>x</sub>(k) - P<sub>y</sub>(k)| / w<sup>2</sup>.
A <I>spectral code</I> C(n,w,&theta;) is a set of weight-w words with pairwise d &ge; &theta;. Let S(n,w,&theta;) be the largest size.
S(n,w,&theta;) &le; c(n,w), the number of distinct signatures; equality holds exactly for &theta; &le; &delta;(n,w), the smallest distance between two classes.
The ordinary constant-weight code problem A(n,d,w) uses Hamming distance and counts shifts as different words; here a shift-closed code of any size is worthless.</P>

<H3>2. Table of c(n,w): the number of distinct power signatures (weight w, length n)</H3>
<P>Values with * are smaller than the number of bracelets (binary necklaces up to rotation and reflection) because of homometric merging; hover for the count.
c(n,w) = c(n,n-w); only w &le; n/2 is shown. Exact integer autocorrelation keys. w = 2: n/2 classes (floor).</P>
{chr(10).join(t1)}

<H3>3. Tables of S(n,w,&theta;): the largest spectral code</H3>
<P>Entries are exact (proven by integer programming) unless written lb-ub. Bold = the code is smaller than the full alphabet c(n,w) (the threshold bites).
Each entry links to a file with one codeword per line (a 0/1 matrix). &mdash; = more than {3500} classes, not computed.</P>
{"".join(f"<H4>&theta; = {t:.2f}</H4>{tb}" for t, tb in maxtabs)}

<H3>4. Table of &delta;(n,w): separation of the full alphabet</H3>
<P>The largest &theta; for which all c(n,w) classes form a code. Large &delta; = every phase-blind-distinguishable word is also well separated.</P>
{chr(10).join(dm)}

<H3>5. Matrices</H3>
<P>The 13 shifts of the (13,4,1) difference set {{0,1,3,9}} (a circulant; optimal for the Hamming problem A(13,6,4) = 13) have ONE signature: S = 1 for these 13 words.</P>
<PRE>{singer_txt}</PRE>
<P>A spectral code for n=12, w=6, &theta; = 0.20 (the proven maximum, {ex['size'] if ex else '?'} words):</P>
<PRE>{ex_txt}</PRE>

<H3>6. Use: an intensity-only receiver in additive noise</H3>
<P>Each position of the 0/1 word is observed with Gaussian noise sigma; the receiver computes |FFT|<sup>2</sup> and decodes by the nearest stored signature (3000 trials per cell).</P>
{chr(10).join(demo_html)}

<H3>7. Link to hyperdimensional computing (HDC) and optical orthogonal codes</H3>
<P>In HDC an item is bound to a role by a cyclic shift (permutation) or circular convolution, items are bundled by superposition, and the readout correlates the bundle against
a codebook (a resonator network iterates this). With sparse binary constant-weight vectors the bundle is the integer sum of shifted codewords and the readout scores are
the circular cross-correlations with the codebook: the same interference as optical CDMA. Three consequences, one of them measured below.</P>
<UL>
<LI>Codebooks whose shifted words pairwise overlap in at most one position (optical orthogonal codes, correlation &le; 1) have a counting guarantee: a true item scores at least w, a false
 hypothesis at most T, so T &le; w-1 bundled items are read out exactly. The size of such a codebook is limited by the Johnson-type bound floor((n-1)/(w(w-1))).</LI>
<LI>Rosenblatt-Seymour (A = U+V against B = U-V) says binding with V and with the reversal of V give words with the same power spectrum: a magnitude-only readout cannot
 distinguish them. This is the homometric ambiguity of section 1 seen as an unbinding ambiguity.</LI>
<LI>Difference sets are the binary analogue of the unit-magnitude (FHRR) role vectors: perfect autocorrelation means the shifts of one word are quasi-orthogonal, which is ideal for roles
 but, by section 6, worthless as an alphabet for an intensity-only receiver (all shifts share one signature). Roles want flat spectra; alphabets for phase-blind readout want spread-out spectra.</LI>
</UL>
<P>Measured (n = 127, w = 4, 10 words, bundle of T random (word, shift) pairs, readout = top-T correlation peaks, 400 trials per cell; fraction of items recovered):</P>
{chr(10).join(hdc_rows)}
<P>Counted: the OOC codebook has maximum auto/cross-correlation 1, the random and Hamming codebooks 3 and 2. Guarantee holds for T &le; 3 (exact, 1.000).</P>

<H3>8. Notes</H3>
<UL>
<LI>All sizes in sections 3 are lower bounds by explicit codes (files) and, where marked exact, upper bounds by the integer program. One cell may be open (lb-ub).</LI>
<LI>c(n,w) agrees with the number of bracelets for w &le; 3 in every cell computed (n &le; 30) and falls below it first for w = 4 (n = 12: the pair {{0,1,4,6}}, {{0,1,3,7}}).</LI>
<LI>Not literature-checked beyond the items in <A HREF="related_work.md">related_work.md</A>. Nothing here is submitted anywhere.</LI>
</UL>
</BODY>
</HTML>
"""
open(os.path.join(OUT, "index.html"), "w", encoding="utf-8").write(page)
print("page written;", len(allcodes), "code files")
