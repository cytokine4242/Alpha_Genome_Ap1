#!/usr/bin/env python3
"""Method diagrams for the fibroblast lab-meeting deck, drawn in the same idiom as the HTML report's diagrams.

Same primitives as Scripts/report_diagrams.py (rounded rects, hairline strokes, arrow markers, orange = edited), but
emitted as standalone SVG and rasterised for slides rather than inlined in a page. The report's blue is deliberately
not used: on these slides blue already means Stable and crimson means Dynamic, and a method diagram has no classes in
it, so native/unchanged is drawn in ink and the only strong colours are the edit (orange), the matched control
(green) and the AP-1 motif itself (gold, as in the data figures).

One diagram per design, so each can sit with the result it produced rather than arriving all at once:

  m1_pipeline     how any effect on these slides is measured at all
  m2_transplant   an element dropped into chromatin that has none          -> the transplant result
  m3_swap         every element written into every other location          -> the element-or-place result
  m4_shuffle      the site's own AP-1 motif, and one partner at a time     -> the motif-specificity results
  m5_scan         the tile scan and single-base saturation                 -> the importance-density result
  m6_titration    the AP-1 ladder and the spacing designs, base by base    -> the titration result

Writes Analysis/slide_figures/m_*.svg and .png (PNG via headless Chrome at 3x, which is the only rasteriser on this
machine - no rsvg-convert, inkscape or cairosvg).

  python3 Scripts/slides/perturbation_diagrams.py [--no-png]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Scripts"))
import report_diagrams as _rd  # noqa: E402
from report_diagrams import L, P, R  # noqa: E402,F401

FS = 1.18   # these diagrams are read on a projector, not on a page


def T(x, y, s, size=12, *a, **k):
    return _rd.T(x, y, s, size * FS, *a, **k)


def TT(x, y, lines, size=12, anchor="middle", lh=15, *a, **k):
    return _rd.TT(x, y, lines, size * FS, anchor, lh * FS, *a, **k)

OUT = ROOT / "Analysis" / "slide_figures"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

INK = "#111111"
EDIT = "#EB6834"        # the report's "edited" orange
CTRL = "#1BAF7A"        # the report's "measured / control" green
GOLD = "#D9A32B"        # the AP-1 motif, as in the data figures
HAIR = "#9A9993"


def defs(uid):
    """Arrowheads in each colour the diagrams use."""
    out = ["<defs>"]
    for suf, col in (("", "currentColor"), ("e", EDIT), ("c", CTRL), ("g", GOLD)):
        out.append(f'<marker id="{uid}{suf}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7"'
                   f' orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
    return "".join(out) + "</defs>"


def svg(inner, w, h, uid, label):
    return (f'<svg viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" aria-label="{label}" '
            f'xmlns="http://www.w3.org/2000/svg">{defs(uid)}{inner}</svg>')


def bar(x, y, w, h=20, col="currentColor", fo=0.055, sw=1.3, rx=4, dash=None):
    return R(x, y, w, h, col, col, fo, sw, rx=rx, dash=dash)


def block(x, y, w, h=20, col=EDIT, fo=0.55, sw=1.2, rx=3, dash=None):
    return R(x, y, w, h, col, col, fo, sw, rx=rx, dash=dash)


def lane_label(x, y, name, sub=None):
    s = T(x, y, name, 13.5, weight="bold")
    if sub:
        s += T(x, y + 17, sub, 11.5, op=0.72)
    return s


# ==========================================================================================================
def m1_pipeline():
    """One site, predicted twice, differenced. Everything else on the slides is a variation on this."""
    s, uid = [], "p1"
    s.append(R(20, 56, 212, 74, "currentColor", "currentColor", 0.03, 1.5, rx=9))
    s.append(TT(126, 82, ["One fibroblast AP-1 site", "TOBIAS footprint, 11 bp", "5,018 sites"], 12, lh=15.5))
    s.append(R(20, 198, 212, 74, "currentColor", "currentColor", 0.03, 1.5, rx=9))
    s.append(TT(126, 224, ["131,072 bp of mm10", "centred on the motif", "fetched once"], 12, lh=15.5))
    s.append(L(126, 130, 126, 196, arrow=uid))
    s.append(T(136, 168, "fetch once", 11, op=0.7))

    s.append(R(262, 48, 212, 74, INK, INK, 0.05, 1.7, rx=9))
    s.append(TT(368, 74, ["NATIVE sequence", "unchanged"], 12.5, weight="bold", lh=17))
    s.append(R(262, 190, 212, 92, EDIT, EDIT, 0.08, 1.7, rx=9))
    s.append(TT(368, 214, ["EDITED sequence"], 12.5, weight="bold"))
    s.append(TT(368, 234, ["shuffle · mutate · insert · swap", "10 replicate edits per site"], 11.5, lh=15))
    s.append(P([(232, 235), (248, 235), (248, 85), (260, 85)], arrow=uid))
    s.append(P([(232, 235), (260, 235)], arrow=uid))

    s.append(R(518, 48, 196, 234, "currentColor", "currentColor", 0.045, 1.7, rx=9))
    s.append(TT(616, 76, ["AlphaGenome"], 13, weight="bold"))
    s.append(TT(616, 100, ["fine-tuned 1 bp ATAC head", "on our own fibroblast Tn5 data"], 11.5, lh=15))
    s.append(L(540, 128, 692, 128, "currentColor", 1, op=0.25))
    s.append(TT(616, 150, ["held-out Spearman 0.68", "(0.27 for the built-in head)"], 11.5, lh=15))
    s.append(L(540, 174, 692, 174, "currentColor", 1, op=0.25))
    s.append(TT(616, 196, ["trackset FB = the mean of", "the four conditions", "(young/old × sed/exercise)"], 11.5, lh=15))
    s.append(TT(616, 258, ["same weights for both lanes"], 11, op=0.7))
    s.append(L(474, 85, 516, 85, INK, 1.7, arrow=uid))
    s.append(L(474, 235, 516, 235, EDIT, 1.7, arrow=uid + "e"))

    s.append(R(758, 48, 206, 74, INK, "none", 0, 1.7, rx=9))
    s.append(TT(861, 72, ["native ATAC"], 12.5, weight="bold"))
    s.append(TT(861, 92, ["mean over the motif,", "the 170 bp and the peak"], 11.5, lh=15))
    s.append(R(758, 190, 206, 92, EDIT, "none", 0, 1.7, rx=9))
    s.append(TT(861, 222, ["edited ATAC"], 12.5, weight="bold", fill=EDIT))
    s.append(TT(861, 244, ["the same regions,", "the same tracks"], 11.5, lh=15))
    s.append(L(714, 85, 756, 85, INK, 1.7, arrow=uid))
    s.append(L(714, 235, 756, 235, EDIT, 1.7, arrow=uid + "e"))

    s.append(R(518, 320, 422, 58, "currentColor", CTRL, 0.10, 1.7, rx=9))
    s.append(TT(729, 344, ["effect = log₂( (edited + c) ⁄ (native + c) )"], 13.5, weight="bold"))
    s.append(TT(729, 364, ["averaged over the replicate edits — one number per site"], 11.5))
    s.append(P([(964, 85), (980, 85), (980, 349), (942, 349)], INK, arrow=uid))
    s.append(P([(861, 282), (861, 318)], EDIT, arrow=uid + "e"))

    s.append(TT(20, 318, ["The two predictions come from the same 131 kb", "window, so a change can only come from the edit.",
                          "c is 5% of the median native value, which keeps", "the ratio stable near zero."],
                11.5, anchor="start", lh=16, op=0.8))
    return svg("".join(s), 1000, 398, uid, "How an effect is measured"), 1000, 398


# ==========================================================================================================
def m2_transplant():
    """An element lifted out and dropped into chromatin that has none, against the same letters shuffled."""
    s, uid = [], "p2"
    s.append(lane_label(20, 62, "The element", "its whole peak, lifted out"))
    s.append(bar(250, 48, 320, 26))
    s.append(block(400, 44, 18, 34, GOLD, 0.75, 1.3))
    s.append(T(409, 36, "AP-1", 11, "middle", weight="bold", fill="#8A6400"))
    s.append(L(409, 82, 409, 118, arrow=uid))

    for x, lab, sub, shuf in ((250, "the element", "as it is", False),
                              (490, "the same letters", "shuffled", True)):
        s.append(R(x, 126, 200, 108, EDIT if shuf else INK, EDIT if shuf else INK, 0.045, 1.5, rx=9))
        s.append(T(x + 100, 150, lab, 12.5, "middle", weight="bold", fill=EDIT if shuf else INK))
        s.append(T(x + 100, 168, sub, 11.5, "middle", op=0.75))
        for j, op_ in enumerate((0.22, 0.45, 1.0)):      # the same pair goes into 20 hosts, drawn as a short stack
            s.append(L(x + 16, 194 + j * 10, x + 184, 194 + j * 10, HAIR, 2.4, op=op_))
        if shuf:
            s.append(block(x + 82, 202, 36, 24, EDIT, 0.45, 1.1))
        else:
            s.append(bar(x + 70, 201, 60, 26, INK, 0.10, 1.3))
            s.append(block(x + 93, 197, 14, 34, GOLD, 0.75, 1.2))
    s.append(R(458, 166, 28, 28, "currentColor", "#ffffff", 1.0, 1.3, rx=14))
    s.append(T(472, 185, "vs", 12, "middle", weight="bold"))
    s.append(T(470, 254, "each into the same 20 closed backgrounds", 11.5, "middle", op=0.82))

    s.append(R(250, 278, 440, 62, "currentColor", CTRL, 0.10, 1.6, rx=9))
    s.append(TT(470, 302, ["effect = log₂( element ⁄ same sequence shuffled )"], 13, weight="bold"))
    s.append(TT(470, 323, ["the host contributes to both, so it cancels"], 11.5, op=0.85))

    s.append(TT(752, 150, ["50 donor elements per class", "20 backgrounds each",
                           "scored at the whole-peak width"], 12.5, anchor="start", lh=23))
    s.append(TT(752, 240, ["A background is a genomic window with no",
                           "called ATAC peak, so there is nothing there",
                           "for the element to add to."], 11.5, anchor="start", lh=17, op=0.8))
    return svg("".join(s), 1100, 366, uid, "Transplant into closed chromatin"), 1100, 366


# ==========================================================================================================
def m3_swap():
    """Every element written into every other element's location, at two widths."""
    s, uid = [], "p3"
    BX0, BW = 250, 420
    rows = (("Donor", "where the sequence comes from", "donor"),
            ("Recipient", "the location it is written into", "recip"),
            ("Result", "donor sequence, recipient location", "res"))
    for i, (name, sub, kind) in enumerate(rows):
        y = 50 + i * 82
        s.append(lane_label(20, y + 16, name, sub))
        if kind == "donor":
            s.append(bar(BX0, y, BW, 26, EDIT, 0.10))
            s.append(block(BX0 + 200, y - 4, 18, 34, GOLD, 0.75, 1.3))
        elif kind == "recip":
            s.append(bar(BX0, y, BW, 26))
            s.append(block(BX0 + 200, y - 4, 18, 34, GOLD, 0.30, 1.2))
        else:
            s.append(bar(BX0, y, BW, 26, "currentColor", 0.035))
            s.append(block(BX0, y - 4, BW, 34, EDIT, 0.26, 1.3))
            s.append(block(BX0 + 200, y - 4, 18, 34, GOLD, 0.75, 1.3))
    s.append(P([(BX0 + BW + 12, 63), (BX0 + BW + 36, 63), (BX0 + BW + 36, 145), (BX0 + BW + 14, 145)],
               EDIT, 1.7, arrow=uid + "e"))

    s.append(R(752, 48, 328, 72, "currentColor", "currentColor", 0.04, 1.5, rx=9))
    s.append(TT(916, 72, ["Two widths"], 12.5, weight="bold"))
    s.append(TT(916, 92, ["the 170 bp window,", "and the whole peak"], 11.5, lh=16))
    s.append(R(752, 136, 328, 72, "currentColor", "currentColor", 0.04, 1.5, rx=9))
    s.append(TT(916, 160, ["62,250 swaps"], 12.5, weight="bold"))
    s.append(TT(916, 180, ["every element into every other", "location; self-swaps excluded"], 11.5, lh=16))

    s.append(R(250, 258, 830, 60, "currentColor", CTRL, 0.10, 1.6, rx=9))
    s.append(TT(665, 281, ["regress the swapped site's accessibility on both:"], 13, weight="bold"))
    s.append(TT(665, 301, ["what the donor element was, and where the recipient location was"], 11.5, op=0.85))
    s.append(TT(250, 340, ["Each donor appears in about 250 of those rows, so the confidence intervals come from",
                           "a cluster bootstrap over donors, not from the 62,250 swaps."],
                11.5, anchor="start", lh=17, op=0.8))
    return svg("".join(s), 1100, 380, uid, "The element swap"), 1100, 380


# ==========================================================================================================
X0, X1 = 250, 1080                      # the drawn peak for the two in-place diagrams, ~1 kb across
CX = (X0 + X1) / 2
PXBP = (X1 - X0) / 1000.0
FOOTPRINTS = (-260, -150, 95, 210)      # illustrative partner footprints


def site_row(s, y, uid, show_footprints=True, scale=True):
    """The reference bar both in-place diagrams are read against, drawn to scale."""
    s.append(lane_label(20, y + 16, "The site", "ATAC peak, median 0.8 kb"))
    s.append(bar(X0, y, X1 - X0, 26))
    if show_footprints:
        for d in FOOTPRINTS:
            s.append(block(CX + d * PXBP - 10, y + 4, 20, 18, CTRL, 0.30, 1.1))
    s.append(block(CX - 8, y - 4, 16, 34, GOLD, 0.75, 1.3))
    s.append(T(CX, y - 12, "AP-1, 11 bp", 12, "middle", weight="bold", fill="#8A6400"))
    if scale:
        for bp in (-400, -200, 0, 200, 400):
            x = CX + bp * PXBP
            s.append(L(x, y + 26, x, y + 33, HAIR, 1))
            s.append(T(x, y + 48, "0" if bp == 0 else f"{bp:+d}", 11, "middle", op=0.7))
        s.append(T(X0 - 12, y + 48, "bp", 11, "end", op=0.7))


def m4_shuffle():
    """The site's own AP-1 motif, and one partner motif at a time against a matched control."""
    s, uid = [], "p4"
    site_row(s, 50, uid)
    s.append(T(X1, 14, "green = a partner motif footprint", 11.5, "end", fill=CTRL, op=0.95))
    s.append(L(20, 132, X1, 132, HAIR, 1, dash="4 4"))

    y = 172
    s.append(lane_label(20, y + 14, "Its own AP-1 motif", "the 11 bp footprint, letters permuted"))
    s.append(bar(X0, y, X1 - X0, 24, "currentColor", 0.035))
    s.append(block(CX - 8, y - 4, 16, 32))
    s.append(T(CX + 18, y + 16, "11 bp", 11.5, fill=EDIT, weight="bold"))

    y = 268
    s.append(lane_label(20, y + 14, "One partner motif", "scored against a matched control"))
    s.append(bar(X0, y, X1 - X0, 24, "currentColor", 0.035))
    for d in (-260, -150, 210):
        s.append(block(CX + d * PXBP - 10, y + 3, 20, 18, CTRL, 0.10, 1.0))
    px, cx_ = CX + 95 * PXBP, CX - 95 * PXBP
    s.append(block(px - 10, y - 4, 20, 32))
    s.append(T(px, y - 12, "shuffled", 11.5, "middle", fill=EDIT, weight="bold"))
    s.append(block(cx_ - 10, y - 4, 20, 32, CTRL, 0.12, 1.3, dash="3 3"))
    s.append(T(cx_, y - 12, "matched control", 11.5, "middle", weight="bold", fill=CTRL))
    s.append(L(cx_ + 12, y + 40, px - 12, y + 40, "currentColor", 1, dash="4 4", op=0.45))
    s.append(T(CX, y + 56, "same width, mirrored distance, no footprint there", 11.5, "middle", op=0.82))

    s.append(TT(20, 356, ["Both effects are measured inside one site, so the accessibility gap between the classes "
                          "cannot reach the comparison. Repeated for every partner family the site carries."],
                11.5, anchor="start", op=0.8))
    return svg("".join(s), 1100, 376, uid, "Shuffling one motif at a time"), 1100, 376


# ==========================================================================================================
def m5_scan():
    """Where along the sequence the dependence sits: tiles, then single bases."""
    s, uid = [], "p5"
    site_row(s, 50, uid, show_footprints=False)
    s.append(L(20, 132, X1, 132, HAIR, 1, dash="4 4"))

    y = 176
    s.append(lane_label(20, y + 14, "Tile scan", "each 10 bp tile replaced in turn"))
    s.append(bar(X0, y, X1 - X0, 24, "currentColor", 0.035))
    tw = 10 * PXBP
    for k in range(-30, 30):
        x = CX + k * tw
        on = k == -8
        s.append(R(x + 0.6, y - 3, tw - 1.2, 30, EDIT, EDIT, 0.55 if on else 0.0,
                   1.0 if on else 0.7, rx=1.5, so=1.0 if on else 0.5))
    s.append(L(CX - 30 * tw, y + 40, CX + 30 * tw, y + 40, "currentColor", 1, op=0.4))
    s.append(T(CX, y + 56, "60 tiles from −300 to +290 bp, one prediction each", 11.5, "middle", op=0.82))

    y = 282
    s.append(lane_label(20, y + 14, "Saturation", "every base, one at a time"))
    s.append(bar(X0, y, X1 - X0, 24, "currentColor", 0.035))
    sw_ = 30 * PXBP
    s.append(block(CX - sw_, y - 3, 2 * sw_, 30))
    for k in range(-5, 6):
        s.append(L(CX + k * sw_ / 5.5, y - 3, CX + k * sw_ / 5.5, y + 27, "#ffffff", 0.9, op=0.85))
    s.append(L(CX - sw_, y + 30, CX - 150, y + 48, HAIR, 0.9, dash="3 3"))
    s.append(L(CX + sw_, y + 30, CX + 150, y + 48, HAIR, 0.9, dash="3 3"))
    s.append(T(CX, y + 62, "61 positions × 3 alternative bases — the only row not drawn to scale",
               11.5, "middle", op=0.82))
    return svg("".join(s), 1100, 376, uid, "Tile scan and saturation"), 1100, 376


# ==========================================================================================================
def m6_titration():
    """The AP-1 affinity ladder and the spacing experiments."""
    s, uid = [], "p6"
    native = "ACGTTCCTGAGTCACCGTAGA"        # motif at 5..15, core TGAGTCA at 7..13
    rows = [
        ("native", native, None, "the site as it is"),
        ("1 base of the core changed", "ACGTTCCTCAGTCACCGTAGA", [8], "the gentlest rung"),
        ("3 bases changed", "ACGTTCCTCATTAACCGTAGA", [8, 10, 12], "the ladder runs from 1 to 7 bases"),
        ("whole core changed", "ACGTTCCCATCAGGCCGTAGA", range(7, 14), "no AP-1 half-site left"),
        ("half site", "ACGTTCCTGACCTACCGTAGA", [10, 11, 12], "one half kept, the other destroyed"),
        ("reverse complement", "ACGTTGGTGACTCAGGGTAGA", [5, 6, 10, 14, 15],
         "same motif, other strand — only 5 of the 11 bases move"),
    ]
    x0, lw = 268, 25
    s.append(T(20, 22, "the same 21 bp excerpt under each edit;  orange = changed bases", 12, op=0.8))
    s.append(L(x0 + 5 * lw, 50, x0 + 16 * lw - 3, 50, "currentColor", 1.4))
    s.append(T(x0 + 10.5 * lw, 44, "11 bp motif", 11.5, "middle"))
    s.append(L(x0 + 7 * lw, 60, x0 + 14 * lw - 3, 60, GOLD, 2.2))
    s.append(T(x0 + 10.5 * lw, 74, "7 bp core  TGA[C/G]TCA", 11, "middle", fill="#8A6400", weight="bold"))
    for r, (name, seq, chg, desc) in enumerate(rows):
        y = 86 + r * 32
        s.append(T(20, y + 17, name, 12, weight="bold" if r == 0 else "normal"))
        for i, ch in enumerate(seq):
            hot = chg is not None and i in chg
            col = EDIT if hot else "currentColor"
            s.append(R(x0 + i * lw, y, lw - 2, 25, col, col, 0.42 if hot else 0.055, 1.1, rx=3))
            s.append(T(x0 + i * lw + (lw - 2) / 2, y + 17, ch, 12.5, "middle", weight="bold"))
        s.append(T(x0 + 21 * lw + 12, y + 17, desc, 11.5, op=0.88))

    # the core upgrade only makes sense as a before and after, since this excerpt's core is already perfect
    uy = 86 + len(rows) * 32 + 14
    s.append(L(20, uy - 14, 1180, uy - 14, HAIR, 1, dash="4 4"))
    s.append(T(20, uy + 17, "core upgraded", 12))
    for j, (bx, seq7, hot, cap) in enumerate(((x0, "TGAGTAA", None, "an imperfect core"),
                                              (x0 + 250, "TGAGTCA", [5], "made TGAGTCA"))):
        for i, ch in enumerate(seq7):
            on = hot is not None and i in hot
            col = EDIT if on else "currentColor"
            s.append(R(bx + i * lw, uy, lw - 2, 25, col, col, 0.42 if on else 0.055, 1.1, rx=3))
            s.append(T(bx + i * lw + (lw - 2) / 2, uy + 17, ch, 12.5, "middle", weight="bold"))
        s.append(T(bx + 3.5 * lw, uy + 40, cap, 11, "middle", op=0.8))
    s.append(L(x0 + 7 * lw + 12, uy + 12, x0 + 250 - 12, uy + 12, "currentColor", 1.3, arrow=uid))
    s.append(T(x0 + 21 * lw + 12, uy + 17, "run only where the core is not already perfect — 1 or 2 bases",
               11.5, op=0.88))

    # spacing block
    ty = uy + 78
    s.append(T(20, ty, "Spacing  —  a second copy, or the same copy moved, at a set distance", 13.5, weight="bold"))
    BX0, BX1 = 404, 1044
    CX = (BX0 + BX1) / 2
    for i, (name, sub, kind) in enumerate((("second copy at ±d", "and a width-matched insert as the control", "dup"),
                                           ("motif moved to ±d", "and the same move to a footprint-free spot", "move"))):
        y = ty + 28 + i * 86
        s.append(lane_label(20, y + 14, name, sub))
        s.append(bar(BX0, y, BX1 - BX0, 22))
        if kind == "dup":
            s.append(block(CX - 8, y - 3, 16, 28, GOLD, 0.75, 1.3))
            s.append(block(CX + 150, y - 3, 16, 28, GOLD, 0.75, 1.3))
            s.append(block(CX - 158, y - 3, 16, 28, CTRL, 0.12, 1.3, dash="3 3"))
            s.append(T(CX + 158, y - 10, "second AP-1", 10.5, "middle", weight="bold", fill="#8A6400"))
            s.append(T(CX - 150, y - 10, "scrambled control, same width", 10.5, "middle", weight="bold", fill=CTRL))
            s.append(L(CX, y + 36, CX - 150, y + 36, "currentColor", 1.2, op=0.55, arrow=uid))
        else:
            s.append(block(CX - 8, y - 3, 16, 28, EDIT, 0.30, 1.3, dash="3 3"))
            s.append(block(CX + 150, y - 3, 16, 28, GOLD, 0.75, 1.3))
            s.append(T(CX, y - 10, "removed", 10.5, "middle", weight="bold", fill=EDIT))
            s.append(T(CX + 158, y - 10, "placed here", 10.5, "middle", weight="bold", fill="#8A6400"))
        s.append(L(CX, y + 36, CX + 158, y + 36, "currentColor", 1.2, op=0.55, arrow=uid))
        s.append(T(CX + 79, y + 50, "d", 11, "middle", op=0.8))
    s.append(T(CX, ty + 182, "d = 20, 40, 60, 100 or 150 bp, either side — the control version is run separately",
               11, "middle", op=0.8))
    return svg("".join(s), 1220, ty + 206, uid, "The AP-1 titration ladder"), 1220, ty + 206


# ==========================================================================================================
DIAGRAMS = {"m1_pipeline": m1_pipeline, "m2_transplant": m2_transplant, "m3_swap": m3_swap,
            "m4_shuffle": m4_shuffle, "m5_scan": m5_scan, "m6_titration": m6_titration}


def rasterise(name, w, h, scale=3):
    """Headless Chrome is the only rasteriser here; it honours the SVG's own pixel size at the given scale."""
    html = OUT / f"{name}.html"
    html.write_text(
        "<!doctype html><meta charset=utf-8><style>html,body{margin:0;padding:0;background:#fff}"
        "svg{display:block;color:#111111;"
        "font-family:Arial,Helvetica,sans-serif}</style><body>"
        + (OUT / f"{name}.svg").read_text() + "</body>")
    r = subprocess.run([CHROME, "--headless", "--disable-gpu", "--hide-scrollbars",
                        f"--force-device-scale-factor={scale}", f"--window-size={w},{h}",
                        f"--screenshot={OUT / (name + '.png')}", str(html)],
                       capture_output=True, text=True)
    html.unlink(missing_ok=True)
    if r.returncode != 0:
        print(f"  ! chrome failed for {name}: {r.stderr.strip()[:200]}")
    return r.returncode == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-png", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in DIAGRAMS.items():
        body, w, h = fn()
        (OUT / f"{name}.svg").write_text(body)
        ok = True if a.no_png else rasterise(name, w, h)
        print(f"  wrote {OUT / (name + '.svg')}  ({w}x{h}){'  + png' if ok and not a.no_png else ''}")
    print("done")
