#!/usr/bin/env python3
"""Method diagrams for the fibroblast deck, drawn to graphical-abstract standard.

An earlier version of this file drew boxes-and-arrows flowcharts in the HTML report's idiom. These are built
instead from a small set of drawn primitives - a locus, a predicted ATAC track, a bZIP dimer on its motif, a
scrambled stretch of sequence - so the mechanism is visible rather than described. One diagram per design, each
sitting with the result it produced.

  d1_pipeline     one window predicted twice and differenced
  d2_transplant   an element dropped into chromatin that has none
  d3_swap         the donor's sequence written into the recipient's place
  d4_shuffle      one motif scrambled at a time, against an empty window of the same width
  d5_scan         a 10 bp window walked across the locus
  d6_titration    the motif weakened one base at a time, and a second copy at a set distance

The grammar is the same in all six, so the audience learns it once:
  left to right   what we do, then what the model predicts
  gold            the AP-1 motif and the protein on it
  orange          sequence we changed
  slate fill      predicted ATAC signal
  dashed grey     a control: the same thing done where there is no motif

Writes Analysis/slide_figures/d*.svg and .png (PNG via headless Chrome at 3x - no rsvg/inkscape/cairosvg here).
  python3 Scripts/slides/perturbation_diagrams.py [--no-png] [--only NAME]
"""
from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Analysis" / "slide_figures"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

INK = "#15181C"
MUTED = "#5A6471"
HAIR = "#AFB5BD"
SIG = "#55637A"          # predicted ATAC signal
GOLD = "#D9A32B"         # the AP-1 motif and its protein
GOLD_D = "#8A6400"
EDIT = "#E2622F"         # sequence we changed
EDIT_D = "#A2380F"
PART = "#4E8E8A"         # a partner transcription factor


# ---------------------------------------------------------------------------------------------- primitives
def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


FS = 1.22   # these are read across a room: a 12.5 unit label lands near 12 pt at slide width


def T(x, y, s, size=13, anchor="start", weight="400", fill=INK, op=1.0, ls=0):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size * FS:.2f}" text-anchor="{anchor}" '
            f'font-weight="{weight}" fill="{fill}" fill-opacity="{op}" letter-spacing="{ls}">{esc(s)}</text>')


def tag(x, y, s, col=MUTED, size=11.5):
    """A count or a unit, set apart from prose."""
    return T(x, y, s, size, fill=col, op=0.95, ls=0.2)


def rule(x0, x1, y, col=HAIR, w=1, dash=None, op=1.0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="{w}"'
            f' stroke-opacity="{op}"{d}/>')


def vrule(x, y0, y1, col=HAIR, w=1, dash=None, op=1.0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x:.1f}" y1="{y0:.1f}" x2="{x:.1f}" y2="{y1:.1f}" stroke="{col}" stroke-width="{w}"'
            f' stroke-opacity="{op}"{d}/>')


def arrow(x0, y0, x1, y1, uid, col=INK, w=1.6, op=0.85):
    return (f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" stroke="{col}" stroke-width="{w}"'
            f' stroke-opacity="{op}" marker-end="url(#{uid}-{col.lstrip("#")})"/>')


def dna(x0, x1, y, col=INK, op=0.55):
    """The locus itself: one quiet line. Detail is saved for the motif."""
    return rule(x0, x1, y, col, 2.0, op=op)


def _noise(x, seed):
    return (math.sin(x * 0.11 + seed) * 0.5 + math.sin(x * 0.29 + seed * 2.1) * 0.3
            + math.sin(x * 0.53 + seed * 3.7) * 0.2)


def track(x0, x1, ybase, peaks, floor=2.0, seed=1.0, col=SIG, fo=0.16, n=200):
    """A predicted ATAC track: gaussians plus a little texture, so it reads as signal rather than a bell curve."""
    pts = []
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        h = floor + 1.6 * abs(_noise(x, seed))
        for cx, w, a in peaks:
            h += a * math.exp(-((x - cx) ** 2) / (2.0 * w * w)) * (1 + 0.10 * _noise(x, seed + 5))
        pts.append((x, ybase - max(h, 0.0)))
    area = (f'<path d="M {x0:.1f} {ybase:.1f} ' + " ".join(f"L {a:.1f} {b:.1f}" for a, b in pts)
            + f' L {x1:.1f} {ybase:.1f} Z" fill="{col}" fill-opacity="{fo}" stroke="none"/>')
    line = ('<polyline points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)
            + f'" fill="none" stroke="{col}" stroke-width="1.5" stroke-opacity="0.9" stroke-linejoin="round"/>')
    return area + line


def dimer(x, y, col=GOLD, s=1.0, op=1.0, tilt=0.0):
    """A bZIP dimer clamping the DNA where its motif is. The deck's one piece of iconography."""
    w, h, out = 8.0 * s, 32 * s, []
    for a in (-17 + tilt, 17 + tilt):
        out.append(f'<g transform="translate({x:.1f},{y:.1f}) rotate({a})" opacity="{op}">'
                   f'<rect x="{-w / 2:.1f}" y="{-h:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{w / 2:.1f}" '
                   f'fill="{col}" fill-opacity="0.80" stroke="{col}" stroke-width="1.3"/></g>')
    return "".join(out)


def blob(x, y, col=PART, s=1.0, op=1.0):
    """A partner factor: deliberately not a dimer, so it never reads as AP-1."""
    return (f'<ellipse cx="{x:.1f}" cy="{y - 9 * s:.1f}" rx="{9 * s:.1f}" ry="{8 * s:.1f}" fill="{col}" '
            f'fill-opacity="{0.70 * op}" stroke="{col}" stroke-width="1.2" stroke-opacity="{op}"/>')


def motif(x, y, w=20, col=GOLD, op=1.0):
    """The motif's footprint on the DNA."""
    return (f'<rect x="{x - w / 2:.1f}" y="{y - 4:.1f}" width="{w:.1f}" height="8" rx="2" fill="{col}" '
            f'fill-opacity="{0.85 * op}" stroke="none"/>')


def scramble(x, y, w, col=EDIT, h=9.0):
    """Changed sequence: a band with a scribble through it, which needs no label."""
    n = max(4, int(w / 5))
    pts = " ".join(f"{x - w / 2 + w * i / n:.1f},{y + (h / 2 if i % 2 else -h / 2):.1f}" for i in range(n + 1))
    return (f'<rect x="{x - w / 2:.1f}" y="{y - h / 2 - 1:.1f}" width="{w:.1f}" height="{h + 2:.1f}" rx="2" '
            f'fill="{col}" fill-opacity="0.18" stroke="none"/>'
            f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.6" stroke-linejoin="round"/>')


def ghost_window(x, y, w, col=MUTED, h=13):
    """A control: the same stretch of sequence, where no motif sits. Drawn to be found, not to hide."""
    n = max(4, int(w / 5))
    pts = " ".join(f"{x - w / 2 + w * i / n:.1f},{y + (4 if i % 2 else -4):.1f}" for i in range(n + 1))
    return (f'<rect x="{x - w / 2 - 3:.1f}" y="{y - h / 2:.1f}" width="{w + 6:.1f}" height="{h:.1f}" rx="2.5" '
            f'fill="{col}" fill-opacity="0.07" stroke="{col}" stroke-width="1.3" stroke-dasharray="3.5 2.5" '
            f'stroke-opacity="0.9"/>'
            f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.5" stroke-opacity="0.75" '
            f'stroke-linejoin="round"/>')


def delta(x, y0, y1, label, col=INK, side="right", size=12.5):
    """The quantity a panel actually produces, bracketed between two heights."""
    k = 7 if side == "right" else -7
    s = [vrule(x, y0, y1, col, 1.4, op=0.9),
         rule(x, x + k, y0, col, 1.4, op=0.9), rule(x, x + k, y1, col, 1.4, op=0.9),
         T(x + (12 if side == "right" else -12), (y0 + y1) / 2 + 4, label, size,
           "start" if side == "right" else "end", "600", col)]
    return "".join(s)


def defs(uid):
    out = ["<defs>"]
    for col in (INK, EDIT, GOLD_D, MUTED):
        out.append(f'<marker id="{uid}-{col.lstrip("#")}" viewBox="0 0 10 10" refX="8.5" refY="5" '
                   f'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
                   f'<path d="M0,0.8 L9,5 L0,9.2 z" fill="{col}"/></marker>')
    return "".join(out) + "</defs>"


def svg(inner, w, h, uid, label):
    return (f'<svg viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" aria-label="{esc(label)}" '
            f'xmlns="http://www.w3.org/2000/svg" font-family="Arial, Helvetica, sans-serif">'
            f'{defs(uid)}{inner}</svg>')


def heading(x, y, s):
    return T(x, y, s, 14, weight="700", fill=INK, ls=0.3)


def caption(x, y, s, anchor="start"):
    return T(x, y, s, 12.5, anchor, fill=MUTED)


# ================================================================================================ d1
def d1_pipeline():
    """The same window, twice, differenced."""
    s, uid = [], "d1"
    LX, LW = 124, 236
    PX, PW = 614, 286
    yA, yB = 108, 226

    s.append(heading(20, yA - 38, "NATIVE"))
    s.append(caption(20, yA - 20, "as it is"))
    s.append(dna(LX, LX + LW, yA))
    s.append(motif(LX + LW / 2, yA))
    s.append(dimer(LX + LW / 2, yA - 4))

    s.append(heading(20, yB - 38, "EDITED"))
    s.append(caption(20, yB - 20, "one thing changed"))
    s.append(dna(LX, LX + LW, yB))
    s.append(scramble(LX + LW / 2, yB, 26))

    s.append(f'<rect x="436" y="{yA - 58}" width="126" height="{yB - yA + 96}" rx="14" fill="{INK}" '
             f'fill-opacity="0.045" stroke="{INK}" stroke-opacity="0.32" stroke-width="1.2"/>')
    s.append(T(499, yA + 46, "Alpha", 14, "middle", "700"))
    s.append(T(499, yA + 63, "Genome", 14, "middle", "700"))
    s.append(T(499, yA + 86, "1 bp ATAC head,", 11.5, "middle", fill=MUTED))
    s.append(T(499, yA + 101, "fine-tuned on", 11.5, "middle", fill=MUTED))
    s.append(T(499, yA + 116, "our own Tn5", 11.5, "middle", fill=MUTED))
    s.append(arrow(LX + LW + 16, yA, 430, yA, uid))
    s.append(arrow(LX + LW + 16, yB, 430, yB, uid, col=EDIT))
    s.append(arrow(568, yA, PX - 14, yA, uid))
    s.append(arrow(568, yB, PX - 14, yB, uid, col=EDIT))

    s.append(track(PX, PX + PW, yA, [(PX + PW / 2, 26, 48)], seed=1.3))
    s.append(dna(PX, PX + PW, yA))
    s.append(track(PX, PX + PW, yB, [(PX + PW / 2, 26, 17)], seed=1.3, col=EDIT, fo=0.14))
    s.append(dna(PX, PX + PW, yB))
    s.append(caption(PX, yA - 62, "predicted accessibility"))

    s.append(delta(PX + PW + 36, yA - 48, yB - 17, "log₂ ratio"))
    s.append(T(PX + PW + 48, yB + 6, "one number", 12, fill=MUTED))
    s.append(T(PX + PW + 48, yB + 22, "per site", 12, fill=MUTED))
    s.append(rule(20, 1060, 292, HAIR, 1, dash="4 4"))
    s.append(tag(20, 318, "both predictions come from the same 131 kb window, so only the edit can move them"))
    return svg("".join(s), 1100, 338, uid, "One window predicted twice and differenced"), 1100, 338


# ================================================================================================ d2
def d2_transplant():
    """Three states, left to right: a real site, an empty place, and what the element does to it."""
    s, uid = [], "d2"
    W = 258
    C = (56, 404, 752)                       # three columns, one story
    yT, yB = 134, 268

    s.append(heading(C[0], 46, "1   A REAL SITE"))
    s.append(track(C[0], C[0] + W, yT, [(C[0] + W / 2, 32, 52)], seed=2.2))
    s.append(dna(C[0], C[0] + W, yT))
    s.append(motif(C[0] + W / 2, yT)); s.append(dimer(C[0] + W / 2, yT - 4))
    s.append(f'<rect x="{C[0] + W / 2 - 86}" y="{yT + 8}" width="172" height="18" rx="4" fill="none" '
             f'stroke="{INK}" stroke-width="1.3" stroke-dasharray="4 3" stroke-opacity="0.6"/>')
    s.append(caption(C[0] + W / 2, yT + 44, "lift out its whole peak", "middle"))

    s.append(heading(C[1], 46, "2   AN EMPTY PLACE"))
    s.append(track(C[1], C[1] + W, yT, [], seed=4.1))
    s.append(dna(C[1], C[1] + W, yT))
    s.append(caption(C[1] + W / 2, yT + 26, "a closed background", "middle"))
    s.append(caption(C[1] + W / 2, yT + 44, "with no called peak", "middle"))

    s.append(heading(C[2], 46, "3   WHAT IT BUILDS THERE"))
    for yy, lab, amp, col, shuf in ((yT, "the element", 46, SIG, False),
                                    (yB, "the same letters, shuffled", 6, EDIT, True)):
        s.append(track(C[2], C[2] + W, yy, [(C[2] + W / 2, 28, amp)], seed=2.6, col=col,
                       fo=0.13 if shuf else 0.16))
        s.append(dna(C[2], C[2] + W, yy))
        if shuf:
            s.append(scramble(C[2] + W / 2, yy, 28))
        else:
            s.append(motif(C[2] + W / 2, yy)); s.append(dimer(C[2] + W / 2, yy - 4))
        s.append(T(C[2], yy + 22, lab, 12.5, fill=EDIT_D if shuf else INK, weight="400" if shuf else "600"))
    s.append(delta(C[2] + W + 26, yT - 46, yB - 8, "effect"))
    s.append(caption(C[2], yB + 46, "the host is in both, so it cancels"))

    s.append(arrow(C[0] + W + 20, yT - 16, C[1] - 20, yT - 16, uid, col=MUTED, w=1.4, op=0.6))
    s.append(arrow(C[1] + W + 20, yT - 16, C[2] - 20, yT - 16, uid, col=MUTED, w=1.4, op=0.6))

    s.append(rule(20, 1080, 320, HAIR, 1, dash="4 4"))
    s.append(tag(20, 346, "20 backgrounds per element   ·   50 elements per class   ·   whole-peak width"))
    return svg("".join(s), 1100, 366, uid, "Transplant into closed chromatin"), 1100, 366


# ================================================================================================ d3
def d3_swap():
    """Is a site's level set by its sequence, or by where it sits?"""
    s, uid = [], "d3"
    AX, AW = 50, 238
    BX, BW = 660, 270
    yD, yR = 112, 246

    s.append(heading(AX, 44, "DONOR"))
    s.append(caption(AX, 62, "a strong element"))
    s.append(track(AX, AX + AW, yD, [(AX + AW / 2, 28, 50)], seed=3.1))
    s.append(dna(AX, AX + AW, yD))
    s.append(motif(AX + AW / 2, yD)); s.append(dimer(AX + AW / 2, yD - 4))

    s.append(heading(AX, yR - 44, "RECIPIENT"))
    s.append(caption(AX, yR - 26, "a weak location"))
    s.append(track(AX, AX + AW, yR, [(AX + AW / 2, 28, 15)], seed=5.4))
    s.append(dna(AX, AX + AW, yR))
    s.append(motif(AX + AW / 2, yR, op=0.45))

    s.append(f'<path d="M {AX + AW + 14} {yD - 12} C {AX + AW + 86} {yD - 12}, {AX + AW + 86} {yR - 12}, '
             f'{AX + AW + 16} {yR - 12}" fill="none" stroke="{EDIT}" stroke-width="1.8" '
             f'marker-end="url(#{uid}-{EDIT.lstrip("#")})"/>')
    s.append(T(AX + AW + 98, (yD + yR) / 2 - 9, "write the donor's sequence", 12.5, fill=EDIT_D, weight="600"))
    s.append(T(AX + AW + 98, (yD + yR) / 2 + 9, "into the recipient's place", 12.5, fill=EDIT_D, weight="600"))

    s.append(vrule(612, 34, 300, HAIR, 1, dash="5 5"))
    s.append(heading(BX, 44, "WHICH HEIGHT DOES IT TAKE?"))
    # the two possibilities, as outlines: neither has happened yet
    for yy, amp, col, lab, sub in ((yD + 24, 50, SIG, "the donor's", "the sequence wins"),
                                   (yR + 18, 15, EDIT, "the recipient's", "the place wins")):
        s.append(track(BX, BX + BW, yy, [(BX + BW / 2, 26, amp)], seed=3.1 if amp > 30 else 5.4,
                       col=col, fo=0.07))
        s.append(dna(BX, BX + BW, yy))
        s.append(motif(BX + BW / 2, yy)); s.append(dimer(BX + BW / 2, yy - 4, op=0.85))
        s.append(T(BX + BW + 22, yy - 6, lab, 12.5, fill=col if col != SIG else SIG, weight="600"))
        s.append(T(BX + BW + 22, yy + 10, sub, 12, fill=MUTED))
    s.append(T(BX + BW / 2, (yD + yR) / 2 + 34, "?", 42, "middle", "700", MUTED, op=0.35))

    s.append(rule(20, 1080, 312, HAIR, 1, dash="4 4"))
    s.append(tag(20, 338, "62,250 swaps   ·   two widths: 170 bp and the whole peak   ·   "
                          "intervals from a cluster bootstrap over the 250 donors"))
    return svg("".join(s), 1100, 358, uid, "The element swap"), 1100, 358


# ================================================================================================ d4
def d4_shuffle():
    """One motif at a time, against the same stretch of sequence where no motif sits."""
    s, uid = [], "d4"
    X0, W = 232, 556
    CX = X0 + W / 2
    PA, PB = CX - 146, CX + 146
    CTRLX = CX - 60
    y0 = 118

    s.append(heading(20, 46, "ONE SITE"))
    s.append(caption(20, 64, "its own motif,"))
    s.append(caption(20, 80, "and its partners"))
    s.append(track(X0, X0 + W, y0, [(CX, 26, 46), (PA, 18, 12), (PB, 18, 10)], seed=6.3))
    s.append(dna(X0, X0 + W, y0))
    s.append(motif(CX, y0)); s.append(dimer(CX, y0 - 4))
    s.append(blob(PA, y0)); s.append(blob(PB, y0))
    s.append(T(CX, y0 + 24, "AP-1", 12, "middle", "700", GOLD_D))
    s.append(T(PA, y0 + 24, "partner", 12, "middle", fill=PART))
    s.append(T(PB, y0 + 24, "partner", 12, "middle", fill=PART))
    s.append(rule(20, 1080, 156, HAIR, 1, dash="4 4"))

    rows = (("scramble AP-1", "ap1", "big drop"),
            ("scramble a partner", "partner", "small drop"),
            ("scramble an empty window", "control", "no drop"))
    for i, (lab, kind, out) in enumerate(rows):
        y = 214 + i * 66
        ctrl = kind == "control"
        s.append(T(20, y + 4, lab, 12.5, weight="400" if ctrl else "600", fill=MUTED if ctrl else INK))
        peaks = {"ap1": [(CX, 26, 14), (PA, 18, 12), (PB, 18, 10)],
                 "partner": [(CX, 26, 46), (PA, 18, 2), (PB, 18, 10)],
                 "control": [(CX, 26, 46), (PA, 18, 12), (PB, 18, 10)]}[kind]
        s.append(track(X0, X0 + W, y, peaks, seed=6.3, floor=1.4, col=SIG if ctrl else EDIT, fo=0.13))
        s.append(dna(X0, X0 + W, y))
        if kind == "ap1":
            s.append(scramble(CX, y, 26)); s.append(blob(PA, y, op=0.9)); s.append(blob(PB, y, op=0.9))
        elif kind == "partner":
            s.append(scramble(PA, y, 24)); s.append(dimer(CX, y - 4)); s.append(blob(PB, y, op=0.9))
        else:
            s.append(scramble(CTRLX, y, 24, col=HAIR)); s.append(ghost_window(CTRLX, y, 30))
            s.append(dimer(CX, y - 4)); s.append(blob(PA, y, op=0.9)); s.append(blob(PB, y, op=0.9))
        s.append(T(X0 + W + 24, y + 4, out, 12.5, fill=MUTED if ctrl else EDIT_D,
                   weight="400" if ctrl else "600"))

    bx = X0 + W + 132
    s.append(f'<path d="M {bx} 274 L {bx + 11} 274 L {bx + 11} 350 L {bx} 350" fill="none" '
             f'stroke="{INK}" stroke-width="1.4" stroke-opacity="0.85"/>')
    s.append(T(bx + 20, 300, "a partner's effect", 12, fill=MUTED))
    s.append(T(bx + 20, 316, "is the difference", 12.5, fill=INK, weight="700"))

    s.append(rule(20, 1080, 376, HAIR, 1, dash="4 4"))
    s.append(tag(20, 402, "the empty window is the same width, at the mirrored distance, with no footprint in it "
                          "— so the comparison never leaves the site"))
    return svg("".join(s), 1100, 422, uid, "Shuffling one motif at a time"), 1100, 422


# ================================================================================================ d5
def d5_scan():
    """Walk a window across the locus; then zoom in and change every base."""
    s, uid = [], "d5"
    X0, W = 176, 740
    CX = X0 + W / 2
    y0 = 124

    s.append(heading(20, 46, "WALK A 10 bp WINDOW ACROSS THE LOCUS"))
    s.append(track(X0, X0 + W, y0, [(CX, 26, 48)], seed=7.7))
    s.append(dna(X0, X0 + W, y0))
    s.append(motif(CX, y0)); s.append(dimer(CX, y0 - 4))
    s.append(T(CX, y0 - 62, "AP-1", 12, "middle", "700", GOLD_D))

    tw = W / 60.0
    for k in range(60):
        x = X0 + (k + 0.5) * tw
        live = k == 21
        s.append(f'<rect x="{x - tw / 2 + 0.7:.1f}" y="{y0 + 12:.1f}" width="{tw - 1.4:.1f}" height="16" rx="1.6" '
                 f'fill="{EDIT}" fill-opacity="{0.78 if live else 0.0}" stroke="{EDIT}" '
                 f'stroke-width="{1.0 if live else 0.7}" stroke-opacity="{1.0 if live else 0.30}"/>')
    s.append(rule(X0, X0 + W, y0 + 42, HAIR, 1))
    s.append(caption(X0, y0 + 62, "−300 bp"))
    s.append(caption(X0 + W, y0 + 62, "+290 bp", "end"))
    s.append(caption(CX, y0 + 62, "60 tiles, one prediction each", "middle"))

    s.append(rule(20, 1080, 230, HAIR, 1, dash="4 4"))
    s.append(heading(20, 268, "THEN ZOOM IN AND CHANGE EVERY BASE"))

    seq, bw, yS, hi = "ACGTTCCTGAGTCACCGT", 24.0, 318, 9
    sx = CX - len(seq) * bw / 2
    for i, ch in enumerate(seq):
        x = sx + i * bw
        on = i == hi
        s.append(f'<rect x="{x + 1.5:.1f}" y="{yS - 11:.1f}" width="{bw - 3:.1f}" height="23" rx="3" '
                 f'fill="{EDIT if on else INK}" fill-opacity="{0.20 if on else 0.045}" '
                 f'stroke="{EDIT if on else HAIR}" stroke-width="{1.2 if on else 0.9}"/>')
        s.append(T(x + bw / 2, yS + 6, ch, 12.5, "middle", "700" if on else "400", EDIT_D if on else INK))
    s.append(f'<path d="M {sx:.1f} {yS - 26:.1f} L {X0 + 290:.1f} {y0 + 30:.1f}" fill="none" stroke="{HAIR}" '
             f'stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(f'<path d="M {sx + len(seq) * bw:.1f} {yS - 26:.1f} L {X0 + 300:.1f} {y0 + 30:.1f}" fill="none" '
             f'stroke="{HAIR}" stroke-width="1" stroke-dasharray="3 3"/>')

    hx = sx + hi * bw + bw / 2
    s.append(vrule(hx, yS + 14, yS + 28, EDIT, 1.2, op=0.8))
    for j, alt in enumerate("ACG"):
        x = hx + (j - 1) * 30
        s.append(f'<rect x="{x - 11:.1f}" y="{yS + 30:.1f}" width="22" height="22" rx="3" fill="{EDIT}" '
                 f'fill-opacity="0.12" stroke="{EDIT}" stroke-width="1.1"/>')
        s.append(T(x, yS + 46, alt, 12.5, "middle", "700", EDIT_D))
    s.append(rule(hx - 41, hx + 41, yS + 26, EDIT, 1.1, op=0.6))
    s.append(caption(hx + 72, yS + 46, "each of the other three bases, one prediction each"))
    s.append(caption(sx - 14, yS + 6, "61 positions", "end"))

    s.append(rule(20, 1080, 398, HAIR, 1, dash="4 4"))
    s.append(tag(20, 424, "each site's own noise floor is removed before its profile is normalised, because a less "
                          "accessible site is a noisier one"))
    return svg("".join(s), 1100, 444, uid, "Tile scan and saturation"), 1100, 444


# ================================================================================================ d6
def d6_titration():
    """Weaken the motif a base at a time; or add a second copy at a set distance."""
    s, uid = [], "d6"
    y = 142
    xs = (132, 352, 572, 792)
    core = "TGAGTCA"
    states = ((core, "perfect core", 1.00, 0.0),
              ("TGACTCA", "1 base changed", 0.78, 5.0),
              ("TCACTTA", "3 bases changed", 0.48, 12.0),
              ("CATCAGG", "core gone", 0.0, 0.0))

    s.append(heading(20, 46, "WEAKEN THE MOTIF, ONE BASE AT A TIME"))
    s.append(arrow(xs[0] + 96, y - 56, xs[-1] - 96, y - 56, uid, col=MUTED, w=1.2, op=0.3))
    s.append(T((xs[0] + xs[-1]) / 2, y - 64, "weaker grip", 12.5, "middle", fill=MUTED))
    for (alt, lab, grip, tilt), x in zip(states, xs):
        s.append(dna(x - 78, x + 78, y))
        if grip > 0:
            s.append(motif(x, y, 46, op=0.35 + 0.65 * grip))
            s.append(dimer(x, y - 4, op=0.32 + 0.68 * grip, tilt=tilt))
        else:
            s.append(scramble(x, y, 46))
        for i, ch in enumerate(alt):
            cx = x - 42 + i * 14
            hot = ch != core[i]
            s.append(f'<rect x="{cx - 6.2:.1f}" y="{y + 20:.1f}" width="12.6" height="17" rx="2.5" '
                     f'fill="{EDIT if hot else INK}" fill-opacity="{0.22 if hot else 0.05}" '
                     f'stroke="{EDIT if hot else HAIR}" stroke-width="{1.1 if hot else 0.9}"/>')
            s.append(T(cx, y + 33, ch, 11.5, "middle", "700" if hot else "400", EDIT_D if hot else INK))
        s.append(caption(x, y + 56, lab, "middle"))
    s.append(T(916, y - 16, "also run", 12.5, fill=MUTED, ls=0.2))
    for i, t in enumerate(("upgrade the core", "half site", "reverse complement")):
        s.append(T(916, y + 6 + i * 19, "·  " + t, 12.5, fill=INK))

    s.append(rule(20, 1080, 240, HAIR, 1, dash="4 4"))
    yy = 320
    s.append(heading(20, 278, "OR ADD A SECOND COPY AT A SET DISTANCE"))
    BX, BW = 250, 520
    s.append(dna(BX, BX + BW, yy))
    s.append(motif(BX + 160, yy)); s.append(dimer(BX + 160, yy - 4))
    s.append(motif(BX + 360, yy)); s.append(dimer(BX + 360, yy - 4))
    s.append(rule(BX + 160, BX + 360, yy + 26, MUTED, 1.2, op=0.7))
    s.append(vrule(BX + 160, yy + 21, yy + 31, MUTED, 1.2, op=0.7))
    s.append(vrule(BX + 360, yy + 21, yy + 31, MUTED, 1.2, op=0.7))
    s.append(T(BX + 260, yy + 20, "d", 12.5, "middle", "600", MUTED))
    s.append(T(BX + BW + 26, yy + 4, "a second AP-1", 12.5, fill=GOLD_D, weight="600"))

    s.append(dna(BX, BX + BW, yy + 76))
    s.append(motif(BX + 160, yy + 76)); s.append(dimer(BX + 160, yy + 72))
    s.append(scramble(BX + 360, yy + 76, 24, col=HAIR)); s.append(ghost_window(BX + 360, yy + 76, 30))
    s.append(T(20, yy + 80, "control", 12.5, fill=MUTED))
    s.append(T(BX + BW + 26, yy + 80, "or the same width of scramble", 12.5, fill=MUTED))
    s.append(caption(BX + 260, yy + 112, "d = 20 to 150 bp, either side", "middle"))
    return svg("".join(s), 1100, 452, uid, "The AP-1 ladder and the spacing designs"), 1100, 452


# ================================================================================================
DIAGRAMS = {"d1_pipeline": d1_pipeline, "d2_transplant": d2_transplant, "d3_swap": d3_swap,
            "d4_shuffle": d4_shuffle, "d5_scan": d5_scan, "d6_titration": d6_titration}


def rasterise(name, w, h, scale=3):
    html = OUT / f"{name}.html"
    html.write_text("<!doctype html><meta charset=utf-8><style>html,body{margin:0;padding:0;background:#fff}"
                    "svg{display:block}</style><body>" + (OUT / f"{name}.svg").read_text() + "</body>")
    r = subprocess.run([CHROME, "--headless", "--disable-gpu", "--hide-scrollbars",
                        f"--force-device-scale-factor={scale}", f"--window-size={w},{h}",
                        f"--screenshot={OUT / (name + '.png')}", str(html)], capture_output=True, text=True)
    html.unlink(missing_ok=True)
    if r.returncode != 0:
        print(f"  ! chrome failed for {name}: {r.stderr.strip()[:200]}")
    return r.returncode == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-png", action="store_true")
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in DIAGRAMS.items():
        if a.only and a.only not in name:
            continue
        body, w, h = fn()
        (OUT / f"{name}.svg").write_text(body)
        ok = True if a.no_png else rasterise(name, w, h)
        print(f"  {name}  {w}x{h}{'  + png' if ok and not a.no_png else ''}")
    print("done")
