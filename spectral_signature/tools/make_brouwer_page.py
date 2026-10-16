"""Write brouwer.html in the layout of A. E. Brouwer's page https://aeb.win.tue.nl/codes/Andw.html
(HTML 4.01 Transitional; h1 title, definition paragraph, 'Jump to', 'Sources' with blockquoted references,
one 'Bounds on ...' h1 + TABLE BORDER=1 per parameter, <th>n</th> row headers, cells 'lb-ub' with the lower bound linking to
a code file, exact values as a single number). Reads table_maxcode.csv and table_classes.csv. Does not touch index.html.
usage: python -I tools/make_brouwer_page.py"""
import csv, os

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
THETAS = ["0.05", "0.10", "0.20", "0.40"]
WS = list(range(3, 11))
NS = list(range(6, 21))
mx, cls = {}, {}
for r in csv.DictReader(open(os.path.join(OUT, "table_maxcode.csv"))):
    mx[(int(r["n"]), int(r["w"]))] = r
for r in csv.DictReader(open(os.path.join(OUT, "table_classes.csv"))):
    cls[(int(r["n"]), int(r["w"]))] = int(r["classes"])


def th_row(ws):
    return "<tr>\n<th>n\\w</th>\n" + "\n".join(f"<th>{w}</th>" for w in ws) + "\n</tr>\n"


def table_classes():
    ws = list(range(3, 16))
    out = ['<table border="1">', th_row(ws)]
    for n in range(6, 31):
        cells = [f"<td>{cls[(n, w)]}</td>" if (n, w) in cls and w <= n // 2 else None for w in ws]
        if any(cells):
            out.append(f"<tr>\n<th>{n}</th>" + "".join(c for c in cells if c) + "\n</tr>\n")
    out.append("</table>")
    return "\n".join(out)


def table_theta(th):
    out = ['<table border="1">', th_row(WS)]
    for n in NS:
        cells, any_ = [], False
        for w in WS:
            if w > n // 2:
                continue
            r = mx.get((n, w))
            if r is None or not r.get(f"size_d{th}"):
                cells.append("<td></td>")
                continue
            any_ = True
            s, u, opt = int(r[f"size_d{th}"]), int(r[f"upper_d{th}"]), r[f"optimal_d{th}"] == "1"
            link = f'<a href="codes/n{n}_w{w}_d{th}.txt">{s}</a>'
            cells.append(f"<td>{link}</td>" if opt else f"<td>{link}-{u}</td>")
        if any_:
            out.append(f"<tr>\n<th>{n}</th>" + "".join(cells) + "\n</tr>\n")
    out.append("</table>")
    return "\n".join(out)


jump = ", ".join(f'<a href="#t{th}">&theta;={th}</a>' for th in THETAS)
page = f"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>Bounds for spectral-signature codes</title>
<meta name="description" content="Tables of the sizes of binary constant weight codes whose words have pairwise different power spectra. Also actual codes are given.">
</head>
<body>
<h1>Bounds for spectral-signature codes</h1>
A word <i>x</i> of length <i>n</i> and weight <i>w</i> has the power signature
P<sub><i>x</i></sub>(<i>k</i>) = |&sum;<sub><i>j</i></sub> <i>x</i><sub><i>j</i></sub> exp(-2&pi;<i>ijk</i>/<i>n</i>)|<sup>2</sup>,
<i>k</i> = 1,...,&lfloor;<i>n</i>/2&rfloor;
(the same information as the cyclic autocorrelation of the support, that is, the multiset of cyclic differences).
Cyclic shifts, reflections and homometric words have the same signature. For <i>w</i> = <i>n</i>/2 a word and its complement have the same signature.
Put <i>d</i>(<i>x</i>,<i>y</i>) = &sum;<sub><i>k</i></sub> |P<sub><i>x</i></sub>(<i>k</i>) - P<sub><i>y</i></sub>(<i>k</i>)| / <i>w</i><sup>2</sup>.
A spectral code <i>C</i>(<i>n</i>,<i>w</i>,&theta;) is a set of binary words of length <i>n</i> and weight <i>w</i>
with <i>d</i>(<i>x</i>,<i>y</i>) &ge; &theta; for all <i>x</i> &ne; <i>y</i>.
<i>S</i>(<i>n</i>,<i>w</i>,&theta;) is the maximum size of such a code.
The tables below give lower (and upper) bounds for this function for small parameters.
A lower bound is demonstrated by giving an explicit code.
Contact: GitHub user <a href="https://github.com/owhs">owhs</a>.
<p>
Jump to <a href="#sig">the number of signatures</a>, {jump}.
<p>

<h1>Sources</h1>
All bounds in these tables were first computed here. Lower bounds are explicit codes (the file is linked from the table, one word per line
as a 0/1 matrix); an entry <i>a</i>-<i>b</i> means <i>a</i> &le; <i>S</i> &le; <i>b</i>, and a single number is exact (proved by integer programming).
<blockquote>
San: Oscar Sanchez,
<br>search harness built with Claude (Anthropic), 2026-10-14.
<a href="README.md">(README)</a> <a href="check_codes.py">(checker)</a> <a href="codes.json">(codes)</a>
<a href="index.html">(longer page with decoding experiments)</a>
</blockquote>
<p>
The ordinary constant weight codes <i>A</i>(<i>n</i>,<i>d</i>,<i>w</i>) are tabulated in
<blockquote>
BSSS:
A. E. Brouwer, J. B. Shearer, N. J. A. Sloane &amp; W. D. Smith,
<br><i>A new table of constant weight codes</i>,
<br>IEEE Trans. Inform. Theory <b>36</b> (1990) 1334-1380.
</blockquote>
There shifts of a word are different words and the distance is the Hamming distance; here two shifts of the same word cannot both be in a code.
The 13 shifts of the (13,4,1) difference set {{0,1,3,9}} give <i>A</i>(13,6,4) = 13 but a single signature, so <i>S</i> = 1 for that family.
Recovery of one binary signal from its Fourier magnitude is treated in
<blockquote>
W. H. Wong, Y. Lou, S. Marchesini &amp; T. Zeng,
<br><i>One-dimensional phase retrieval: regularization, box relaxation and uniqueness</i>,
<br>arXiv:1904.10157.
</blockquote>
and homometric sets in
<blockquote>
A. Rosenblatt &amp; P. D. Seymour,
<br><i>The structure of homometric sets</i>,
<br>SIAM J. Algebraic Discrete Methods <b>3</b> (1982) 343-350.
</blockquote>
<p>
For <i>n</i> = 16, <i>w</i> = 6, 12 words and Gaussian noise 0.2 on the 0/1 positions, a code of maximal minimum spectral distance has symbol
error rate 0.009, against 0.117 for greedy Hamming codes of the same size (for <i>n</i> = 18, <i>w</i> = 6, 16 words: 0.006 against 0.113).
At noise 0.5 the gain is small (0.474 against 0.666 for <i>n</i> = 16). The distance used is the L1 distance between noise-free signatures.
No theorem is claimed; the literature was searched but not read in full (<a href="related_work.md">related_work.md</a>).
<p>

<hr>
<h1><a name="sig">The number of signatures <i>c</i>(<i>n</i>,<i>w</i>)</a></h1>
<i>S</i>(<i>n</i>,<i>w</i>,&theta;) &le; <i>c</i>(<i>n</i>,<i>w</i>), with equality if and only if &theta; is at most the smallest distance between two classes.
Since <i>c</i>(<i>n</i>,<i>w</i>) = <i>c</i>(<i>n</i>,<i>n-w</i>) only <i>w</i> &le; <i>n</i>/2 is shown. For <i>w</i> &le; 3 the number equals the number of bracelets;
the first difference is <i>w</i> = 4, <i>n</i> = 12, where {{0,1,4,6}} and {{0,1,3,7}} have the same signature.
{table_classes()}
""" + "\n".join(f'<p>\n<hr>\n<h1><a name="t{th}">Bounds on <i>S</i>(<i>n</i>,<i>w</i>,{th})</a></h1>\n{table_theta(th)}' for th in THETAS) + """
<p>
<hr>
</body>
</html>
"""
open(os.path.join(OUT, "brouwer.html"), "w", encoding="utf-8").write(page)
print("wrote", os.path.join(OUT, "brouwer.html"), len(page), "bytes")
