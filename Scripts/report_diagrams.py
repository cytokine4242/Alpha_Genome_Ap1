"""Inline-SVG diagrams for the HTML report. Each function returns a <figure> with an <svg> (role=img, aria-label) and a caption.
Strokes and text use currentColor so they follow the page theme; blue = native / Stable, orange = edited / Dynamic, green = measured data."""
from __future__ import annotations

BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"


def T(x, y, s, size=12, anchor="start", weight="normal", fill="currentColor", op=1.0):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" fill-opacity="{op}" font-family="inherit">{s}</text>'


def TT(x, y, lines, size=12, anchor="middle", lh=15, weight="normal", fill="currentColor", op=1.0):
    return "".join(T(x, y + i * lh, ln, size, anchor, weight, fill, op) for i, ln in enumerate(lines))


def R(x, y, w, h, stroke="currentColor", fill="none", fo=0.0, sw=1.3, rx=7, dash=None, so=1.0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" fill-opacity="{fo}" stroke="{stroke}" stroke-opacity="{so}" stroke-width="{sw}"{d}/>'


def L(x1, y1, x2, y2, stroke="currentColor", sw=1.4, dash=None, arrow=None, op=1.0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{arrow})"' if arrow else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-opacity="{op}" stroke-width="{sw}"{d}{m}/>'


def P(points, stroke="currentColor", sw=1.4, arrow=None, dash=None, fill="none", fo=0.0):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{arrow})"' if arrow else ""
    pts = " ".join(f"{a},{b}" for a, b in points)
    return f'<polyline points="{pts}" fill="{fill}" fill-opacity="{fo}" stroke="{stroke}" stroke-width="{sw}"{d}{m}/>'


def defs(uid):
    return (f'<defs><marker id="{uid}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker>'
            f'<marker id="{uid}b" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{BLUE}"/></marker>'
            f'<marker id="{uid}o" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{ORANGE}"/></marker></defs>')


def fig(svg_inner, w, h, uid, label, caption):
    return (f'<figure class="dia"><svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">{defs(uid)}{svg_inner}</svg>'
            f'<figcaption>{caption}</figcaption></figure>')


# ---------------------------------------------------------------- D1: pipeline
def d1_pipeline():
    s = []
    # column 1: inputs
    s.append(R(20, 40, 190, 70)); s.append(TT(115, 66, ["Fos site (TOBIAS)", "bound in YS / MAS / MAE?"], 12.5))
    s.append(R(20, 170, 190, 70)); s.append(TT(115, 196, ["mm10 genome window", "131,072 bp around the motif"], 12.5))
    s.append(L(115, 110, 115, 168, arrow="a1")); s.append(T(123, 144, "fetch", 11, op=0.8))
    # lanes
    s.append(R(265, 40, 190, 70, BLUE, BLUE, 0.08, 1.8)); s.append(TT(360, 66, ["NATIVE sequence", "unchanged"], 12.5, weight="bold"))
    s.append(R(265, 170, 190, 70, ORANGE, ORANGE, 0.08, 1.8)); s.append(TT(360, 190, ["EDITED sequence", "shuffle x10 / mutate / insert", "(see 'the experiments')"], 12, lh=14.5))
    s.append(P([(210, 205), (238, 205), (238, 75), (263, 75)], arrow="a1"))
    s.append(P([(210, 205), (263, 205)], arrow="a1"))
    # model
    s.append(R(515, 40, 170, 200, "currentColor", "currentColor", 0.05, 1.8))
    s.append(TT(600, 70, ["AlphaGenome", "(mouse, all-folds", "weights)"], 12.5, weight="bold"))
    s.append(TT(600, 135, ["one-hot DNA in", "predicted tracks out:", "ATAC at 1 bp", "RNA at 128 bp"], 12, lh=16))
    s.append(L(455, 75, 513, 75, BLUE, 1.8, arrow="a1b")); s.append(L(455, 205, 513, 205, ORANGE, 1.8, arrow="a1o"))
    # readouts
    s.append(R(745, 40, 200, 70, BLUE, "none", 0, 1.8)); s.append(TT(845, 62, ["native readout", "mean ATAC over motif,", "170 bp, peak; RNA nearby"], 12, lh=14.5))
    s.append(R(745, 170, 200, 70, ORANGE, "none", 0, 1.8)); s.append(TT(845, 192, ["edited readout", "same regions,", "same tracks"], 12, lh=14.5))
    s.append(L(685, 75, 743, 75, BLUE, 1.8, arrow="a1b")); s.append(L(685, 205, 743, 205, ORANGE, 1.8, arrow="a1o"))
    # effect
    s.append(R(540, 300, 420, 64, "currentColor", GREEN, 0.10, 1.8))
    s.append(TT(750, 325, ["effect per site = log2( (edited + c) / (native + c) )", "mean over the 10 replicates, matched readout"], 13, weight="bold", lh=18))
    s.append(P([(945, 75), (972, 75), (972, 332), (962, 332)], BLUE, arrow="a1b")); s.append(P([(845, 240), (845, 298)], ORANGE, arrow="a1o"))
    s.append(T(855, 272, "compared", 11, op=0.8))
    s.append(R(540, 400, 420, 56, "currentColor", "none", 0, 1.3)); s.append(TT(750, 424, ["classes compared: median effect per class, Mann-Whitney", "tests, adjusted models, confidence tiers"], 12.5, lh=17))
    s.append(L(750, 364, 750, 398, arrow="a1")); s.append(T(760, 386, "one number per site", 11, op=0.8))
    s.append(TT(20, 330, ["The comparison that matters is NATIVE vs EDITED through the", "same model; both come from the same window, so a change in", "the prediction can only come from the edit."], 12.5, lh=18, anchor="start", op=0.9))
    s.append(TT(20, 410, ["c = a small constant (5% of the median native value)", "that keeps the ratio stable when a value is near zero."], 11.5, lh=16, anchor="start", op=0.75))
    return fig("".join(s), 980, 480, "a1", "Pipeline from a Fos site to an effect size",
               "Pipeline for every experiment: one 131 kb window per site is predicted twice, unedited and edited, by the same model. Both predictions are summarised over the same regions and compared as a log2 ratio, one number per site, which feeds the class statistics.")


# ---------------------------------------------------------------- D2: window anatomy
def d2_window():
    s = []
    x0, x1 = 40, 940  # 131 kb bar
    c = (x0 + x1) // 2
    kb = (x1 - x0) / 131.072
    s.append(T(x0, 16, "Whole window: 131,072 bp, predicted in one pass", 12.5, weight="bold"))
    # +-50 kb bracket above the bar
    b0, b1 = c - 50 * kb, c + 50 * kb
    s.append(L(b0, 34, b1, 34, GREEN, 2.2)); s.append(L(b0, 28, b0, 40, GREEN, 2.2)); s.append(L(b1, 28, b1, 40, GREEN, 2.2))
    s.append(T(b1 - 6, 26, "RNA: mean and max over +-50 kb around the site", 12, "end", fill=GREEN, weight="bold"))
    s.append(R(x0, 46, x1 - x0, 22, "currentColor", "currentColor", 0.06, 1.4, rx=4))
    tss = c + 205
    s.append(R(tss - 2 * kb, 46, 4 * kb, 22, GREEN, GREEN, 0.4, 1.2, rx=2))
    s.append(L(tss, 96, tss, 70, "currentColor", 1.6, arrow="a2"))
    s.append(T(tss + 10, 96, "nearest TSS; green block = RNA +-2 kb around it", 11.5, op=0.95))
    s.append(R(c - 3, 44, 6, 26, ORANGE, ORANGE, 0.8, 1, rx=1))
    # zoom 1: peak
    y1 = 160
    s.append(P([(c - 3, 72), (150, y1)], "currentColor", 0.8, dash="3 3")); s.append(P([(c + 3, 72), (830, y1)], "currentColor", 0.8, dash="3 3"))
    s.append(T(20, y1 + 11, "Zoom 1", 12, weight="bold")); s.append(T(20, y1 + 26, "about 1.2 kb", 11.5, op=0.85))
    s.append(R(150, y1, 680, 22, "currentColor", "currentColor", 0.04, 1.2, rx=3))
    s.append(R(235, y1, 510, 22, ORANGE, ORANGE, 0.22, 1.6, rx=3)); s.append(T(490, y1 + 16, "PEAK (median 0.85 kb): peak shuffle edits all of it; ATAC readout = mean over it", 11.5, "middle"))
    # zoom 2: 170 bp
    y2 = 262
    s.append(P([(441, y1 + 22), (150, y2)], "currentColor", 0.8, dash="3 3")); s.append(P([(539, y1 + 22), (830, y2)], "currentColor", 0.8, dash="3 3"))
    s.append(T(20, y2 + 11, "Zoom 2", 12, weight="bold")); s.append(T(20, y2 + 26, "170 bp window", 11.5, op=0.85))
    s.append(R(150, y2, 680, 22, ORANGE, ORANGE, 0.22, 1.6, rx=3)); s.append(T(490, y2 + 16, "170 bp = motif midpoint +-85 bp: 170 bp shuffle edits this; ATAC readout = mean over it", 11.5, "middle"))
    # zoom 3: motif letters
    y3 = 372
    s.append(P([(470, y2 + 22), (330, y3)], "currentColor", 0.8, dash="3 3")); s.append(P([(510, y2 + 22), (640, y3)], "currentColor", 0.8, dash="3 3"))
    s.append(T(20, y3 + 12, "Zoom 3", 12, weight="bold")); s.append(T(20, y3 + 27, "the 11 bp motif", 11.5, op=0.85)); s.append(T(20, y3 + 42, "(example sequence)", 11.5, op=0.85))
    seq = "CCTGAGTCACC"
    for i, ch in enumerate(seq):
        x = 330 + i * 28
        core = 2 <= i <= 8
        s.append(R(x, y3, 27, 30, ORANGE if core else "currentColor", ORANGE if core else "currentColor", 0.30 if core else 0.06, 1.3, rx=3))
        s.append(T(x + 13.5, y3 + 21, ch, 15, "middle", weight="bold"))
    s.append(L(330 + 2 * 28, y3 + 40, 330 + 9 * 28 - 1, y3 + 40, ORANGE, 2)); s.append(T(330 + 5.5 * 28, y3 + 57, "AP-1 core TGA[C/G]TCA (7 bp)", 12, "middle", fill=ORANGE, weight="bold"))
    s.append(T(330 + 5.5 * 28, y3 + 76, "motif shuffle edits all 11 bp; ATAC readout = mean over them", 11.5, "middle", op=0.9))
    s.append(T(330 + 5.5 * 28, y3 + 94, "mutation arm: core edit changes the 7 bp core; flank edit changes bases outside it", 11.5, "middle", op=0.9))
    return fig("".join(s), 980, 490, "a2", "Anatomy of one prediction window: peak, 170 bp window, motif, core, TSS and RNA windows",
               "Nested regions around one site (not to scale at the lower zooms). Each shuffle edits one region and the matching ATAC readout is the mean predicted signal over that same region. RNA readouts use windows around the site (+-50 kb) and around the nearest TSS (+-2 kb).")


# ---------------------------------------------------------------- D3: experiments
def seqbar(x, y, w, col, seg=None, h=14, fo=0.12):
    out = R(x, y, w, h, col, col, fo, 1.3, rx=3)
    if seg:
        sx, sw = seg
        out += R(x + sx, y, sw, h, ORANGE, ORANGE, 0.55, 1.2, rx=2)
    return out


def d3_experiments():
    s = []

    def card(x, y, w, h, title, status, okc):
        o = R(x, y, w, h, "currentColor", "currentColor", 0.03, 1.4, rx=10)
        o += T(x + 14, y + 24, title, 14, weight="bold")
        o += R(x + w - 128, y + 9, 116, 22, okc, okc, 0.15, 1.3, rx=11) + T(x + w - 70, y + 24, status, 11.5, "middle", fill=okc, weight="bold")
        return o

    # card 1: shuffle
    x, y = 20, 20
    s.append(card(x, y, 470, 250, "1  Shuffle arm", "done (CM, FB)", GREEN))
    s.append(T(x + 14, y + 52, "CHANGE", 11, weight="bold", op=0.7))
    s.append(T(x + 14, y + 74, "native", 11.5)); s.append(seqbar(x + 115, y + 63, 335, BLUE))
    for i, (nm, seg) in enumerate((("motif", (165, 4)), ("170 bp", (148, 38)), ("peak", (70, 190)))):
        yy = y + 96 + i * 24
        s.append(T(x + 14, yy + 11, nm + " shuffle", 11.5)); s.append(seqbar(x + 115, yy, 335, ORANGE, seg))
    s.append(T(x + 14, y + 176, "10 random shuffles of that region (letter composition kept)", 11.5, op=0.9))
    s.append(T(x + 14, y + 198, "MEASURED", 11, weight="bold", op=0.7))
    s.append(TT(x + 14, y + 216, ["predicted ATAC over the shuffled region (and RNA near the TSS);", "effect = log2(shuffled / native), averaged over the 10 shuffles"], 11.5, "start", lh=15))
    # card 2: mutation
    x, y = 510, 20
    s.append(card(x, y, 470, 250, "2  Mutation arm", "done (CM, FB)", GREEN))
    s.append(T(x + 14, y + 52, "CHANGE", 11, weight="bold", op=0.7))
    rows = (("random", "every base replaced by a different one", (0, 360)), ("dinucleotide", "order shuffled, dinucleotides kept", (160, 40)),
            ("core", "TGA[C/G]TCA 7 bp disrupted", (178, 4)), ("flanks", "motif kept, flanking bases mutated", (160, 40)))
    for i, (nm, desc, seg) in enumerate(rows):
        yy = y + 62 + i * 30
        s.append(T(x + 14, yy + 11, nm, 11.5)); s.append(seqbar(x + 100, yy, 150, ORANGE if nm != "flanks" else ORANGE, (60, 10) if nm in ("core", "random") else (40, 70) if nm == "dinucleotide" else (10, 40)))
        s.append(T(x + 262, yy + 11, desc, 11, op=0.9))
    s.append(T(x + 14, y + 184, "each edit is applied to the motif, the 170 bp or the peak region", 11, op=0.85))
    s.append(T(x + 14, y + 204, "MEASURED", 11, weight="bold", op=0.7))
    s.append(TT(x + 14, y + 222, ["same as the shuffle arm; reported as dependency = -log2FC", "(random 5x, dinucleotide 5x, core 3x, flank 5x per site)"], 11.5, "start", lh=15))
    # card 3: insertion
    x, y = 20, 290
    s.append(card(x, y, 470, 290, "3  Insertion experiments", "done (CM)", GREEN))
    s.append(T(x + 14, y + 50, "A  PEAK SWAP", 11, weight="bold", op=0.7))
    s.append(T(x + 14, y + 72, "recipient", 11.5)); s.append(seqbar(x + 90, y + 61, 330, BLUE, (130, 70)))
    s.append(T(x + 14, y + 100, "donor peak", 11.5)); s.append(seqbar(x + 90, y + 89, 70, ORANGE, None))
    s.append(P([(x + 125, y + 89), (x + 125, y + 80), (x + 255, y + 80), (x + 255, y + 77)], ORANGE, arrow="a3o"))
    s.append(TT(x + 14, y + 128, ["4 recipient windows x 73 donors: the recipient's peak is replaced by the donor's;", "ATAC over the inserted piece is compared with the donor's own native", "ATAC (sequence effect) and the recipient's (location effect)"], 11, "start", lh=14.5, op=0.95))
    s.append(T(x + 14, y + 178, "B  MOTIF INSERTION", 11, weight="bold", op=0.7))
    s.append(T(x + 14, y + 200, "random 131 kb", 11.5)); s.append(seqbar(x + 110, y + 189, 310, BLUE))
    s.append(T(x + 14, y + 224, "shuffled", 11.5)); s.append(seqbar(x + 110, y + 213, 310, ORANGE, None, fo=0.2))
    s.append(R(x + 255, y + 211, 18, 18, ORANGE, ORANGE, 0.7, 1.2, rx=2)); s.append(T(x + 282, y + 225, "+ TGAGTCA (or a shuffled copy)", 11.5))
    s.append(TT(x + 14, y + 254, ["50 backgrounds: ATAC at the centre with the motif vs with a shuffled-motif control", "tests whether the motif alone can create accessibility"], 11, "start", lh=14.5, op=0.95))
    # card 4: MPRA
    x, y = 510, 290
    s.append(card(x, y, 470, 290, "4  MPRA comparison", "done", GREEN))
    s.append(T(x + 14, y + 50, "MPRA (measured, fibroblast, episomal)", 11.5, weight="bold", fill=GREEN))
    for i, nm in enumerate(("native 170 bp", "core motif scrambled", "fully scrambled")):
        yy = y + 62 + i * 22
        s.append(T(x + 14, yy + 11, nm, 11.5)); s.append(seqbar(x + 150, yy, 120, GREEN if i == 0 else ORANGE, (50, 12) if i == 1 else (0, 120) if i == 2 else None))
    s.append(T(x + 285, y + 100, "activity = log2(RNA / DNA)", 11.5, fill=GREEN))
    s.append(T(x + 14, y + 140, "AlphaGenome (predicted, genomic context, same sequences)", 11.5, weight="bold", fill=BLUE))
    s.append(TT(x + 14, y + 160, ["predicted ATAC for native / motif-shuffled / fully shuffled", "(10 shuffles each) in the element's own 131 kb window"], 11.5, "start", lh=15))
    s.append(L(x + 235, y + 190, x + 235, y + 210, arrow="a3"))
    s.append(R(x + 14, y + 212, 442, 62, "currentColor", "none", 0, 1.3, rx=8))
    s.append(TT(x + 235, y + 232, ["compare per element: MPRA change (native - scrambled)", "vs predicted accessibility change; and native vs native", "Spearman rho: about 0 for the changes, 0.2 for native levels"], 11.5, lh=15))
    return fig("".join(s), 1000, 600, "a3", "The four experiments: what each one changes and what it measures",
               "What each experiment edits and what is read out. Orange marks the edited sequence, blue the unedited one, green the measured MPRA data. All four have been run; the shuffle and mutation arms cover cardiomyocyte and fibroblast sites.")


# ---------------------------------------------------------------- D4: classes
def d4_classes():
    s = []
    cols = (("YS", 260), ("MAS", 330), ("MAE", 400))
    s.append(T(20, 24, "Binding call at the site (filled = bound)", 12.5, weight="bold"))
    for nm, x in cols:
        s.append(T(x + 25, 52, nm, 12.5, "middle", weight="bold"))
    s.append(P([(285, 40), (285, 30), (355, 30), (355, 40)], "currentColor", 1.2)); s.append(T(320, 22, "AGE: YS vs MAS", 11.5, "middle", fill=BLUE, weight="bold"))
    rows = (("111", 1, 1, 1, "Stable", "Stable", 2907), ("110", 1, 1, 0, "Stable", "Exercise Down", 90), ("011", 0, 1, 1, "Age Up", "Stable", 210), ("010", 0, 1, 0, "Age Up", "Exercise Down", 186),
            ("101", 1, 0, 1, "Age Down", "Exercise Up", 79), ("100", 1, 0, 0, "Age Down", "YS only (neither)", 100), ("001", 0, 0, 1, "MAE only (neither)", "Exercise Up", 107))
    s.append(T(540, 52, "Age class", 12.5, weight="bold", fill=BLUE)); s.append(T(690, 52, "Exercise class", 12.5, weight="bold", fill=ORANGE)); s.append(T(870, 52, "Sites", 12.5, "middle", weight="bold"))
    s.append(T(20, 52, "Pattern", 12.5, weight="bold"))
    for i, (pat, a, b, c, age, ex, n) in enumerate(rows):
        y = 64 + i * 34
        s.append(T(20, y + 21, pat + ("   -> Stable site" if i == 0 else "   -> Dynamic site"), 12.5, weight="bold" if i == 0 else "normal"))
        for (nm, x), v in zip(cols, (a, b, c)):
            s.append(R(x, y, 50, 26, BLUE if v else "currentColor", BLUE if v else "none", 0.55 if v else 0, 1.2 if v else 0.8, rx=5, so=1 if v else 0.35))
        s.append(T(540, y + 18, age, 12.5, fill=BLUE if age != "Stable" else "currentColor")); s.append(T(690, y + 18, ex, 12.5, fill=ORANGE if ex != "Stable" else "currentColor"))
        s.append(T(870, y + 18, f"{n:,}", 12.5, "middle"))
    s.append(P([(355, 322), (355, 332), (425, 332), (425, 322)], "currentColor", 1.2)); s.append(T(390, 350, "EXERCISE: MAS vs MAE", 11.5, "middle", fill=ORANGE, weight="bold"))
    s.append(TT(540, 322, ["Stable = bound in all three.", "Age and exercise labels come from the", "same three calls, so they are not", "independent (only 7 patterns)."], 11.5, "start", lh=15, op=0.9))
    return fig("".join(s), 980, 372, "a4", "How binding patterns map to Stable, age and exercise classes",
               "Each Fos site has a bound / not-bound call in young sedentary (YS), aged sedentary (MAS) and aged exercised (MAE) hearts. Age class compares the first two calls, exercise class the last two; a site can therefore carry an age label and an exercise label at once.")


# ---------------------------------------------------------------- D5: worked example
def d5_example():
    s = []
    s.append(T(20, 26, "One site, one shuffle width (illustrative numbers)", 12.5, weight="bold"))
    # chart area
    s.append(L(70, 60, 70, 250, "currentColor", 1.2, op=0.7)); s.append(L(70, 250, 440, 250, "currentColor", 1.2, op=0.7))
    s.append(T(30, 150, "predicted ATAC", 11.5, "middle", op=0.8).replace('<text ', '<text transform="rotate(-90 30 150)" '))
    s.append(R(100, 250 - 160, 70, 160, BLUE, BLUE, 0.35, 1.6, rx=3)); s.append(T(135, 276, "native", 12, "middle")); s.append(T(135, 250 - 166, "0.80", 12, "middle", weight="bold"))
    vals = [0.14, 0.2, 0.22, 0.25, 0.27, 0.3, 0.19, 0.24, 0.33, 0.16]
    for i, v in enumerate(vals):
        cx = 232 + (i % 5) * 18 + (9 if i >= 5 else 0)
        s.append(f'<circle cx="{cx}" cy="{250 - v * 200}" r="5" fill="{ORANGE}" fill-opacity="0.85"/>')
    mean = sum(vals) / len(vals)
    s.append(L(215, 250 - mean * 200, 345, 250 - mean * 200, ORANGE, 2, dash="5 3")); s.append(T(350, 250 - mean * 200 + 4, f"mean {mean:.2f}", 11.5, fill=ORANGE))
    s.append(T(290, 276, "10 shuffles", 12, "middle")); s.append(L(100, 250 - 160, 440, 250 - 160, BLUE, 1, dash="2 4", op=0.6))
    # effect
    s.append(R(500, 70, 460, 82, "currentColor", GREEN, 0.10, 1.6, rx=9))
    s.append(TT(730, 96, ["effect = log2( (0.23 + c) / (0.80 + c) )  =  -1.8", "negative = this site lost accessibility when shuffled", "(c is 5% of the median native value, here about 0.01)"], 12.5, lh=19, weight="bold"))
    s.append(L(730, 152, 730, 184, arrow="a5")); s.append(T(742, 172, "repeat for every site in a class", 11.5, op=0.8))
    s.append(R(500, 186, 460, 56, "currentColor", "none", 0, 1.4, rx=9)); s.append(TT(730, 210, ["class result = median effect over its sites, with a 95% CI", "(e.g. Stable -1.75, Dynamic -1.23 for the 170 bp shuffle)"], 12.5, lh=18))
    s.append(L(730, 242, 730, 274, arrow="a5")); s.append(T(742, 262, "compare classes", 11.5, op=0.8))
    s.append(R(500, 276, 460, 70, "currentColor", "none", 0, 1.4, rx=9))
    s.append(TT(730, 298, ["Mann-Whitney test vs the reference class (Stable), corrected for many", "tests (FDR); effect size reported next to the stars; adjusted model asks", "if the gap survives after accounting for accessibility and peak width"], 12, lh=17))
    return fig("".join(s), 980, 366, "a5", "A worked example of how one site's effect is calculated and compared between classes",
               "How one number per site is built and used. Ten shuffled predictions are averaged against the native prediction on a log2 scale; sites are then summarised per class and classes are compared with significance tests and adjusted models.")


# ---------------------------------------------------------------- D6: evidence map
def d6_evidence():
    s = []

    def node(x, y, w, h, title, sub, col):
        return R(x, y, w, h, col, col, 0.10, 1.8, rx=10) + T(x + w / 2, y + 26, title, 13.5, "middle", weight="bold") + TT(x + w / 2, y + 46, sub, 11.5, "middle", lh=14.5, op=0.9)

    s.append(node(40, 20, 250, 90, "snATAC (measured)", ["pseudobulk by condition,", "CPM per kb at each peak"], GREEN))
    s.append(node(690, 20, 250, 90, "TOBIAS binding calls", ["bound in YS / MAS / MAE,", "footprint score; site classes"], ORANGE))
    s.append(node(360, 170, 260, 100, "AlphaGenome (predicted)", ["native and edited sequences:", "accessibility and RNA tracks"], BLUE))
    s.append(node(40, 340, 250, 90, "MPRA (measured)", ["fibroblast, episomal activity of", "native / scrambled 170 bp"], GREEN))
    # comparison between measured and classes (line, no arrowheads) and derivation (arrow)
    s.append(L(292, 46, 688, 46, GREEN, 1.8)); s.append(T(490, 38, "A4: do the class labels match the measured accessibility change?", 11.5, "middle", op=0.95))
    s.append(L(292, 92, 686, 92, "currentColor", 1.4, dash="5 4", arrow="a6")); s.append(T(490, 84, "calls are computed FROM this ATAC data", 11.5, "middle", op=0.95))
    # classes define the edited sites
    s.append(P([(815, 112), (815, 215), (622, 215)], ORANGE, 1.8, arrow="a6o")); s.append(T(823, 150, "defines which", 11.5, op=0.95)); s.append(T(823, 165, "sites are edited", 11.5, op=0.95)); s.append(T(823, 180, "and how they", 11.5, op=0.95)); s.append(T(823, 195, "are grouped", 11.5, op=0.95))
    s.append(P([(905, 112), (905, 252), (622, 252)], "currentColor", 1.4, dash="3 3")); s.append(T(898, 272, "C4: footprint score vs predicted", 11.5, "end", op=0.95)); s.append(T(898, 287, "dependence (partly circular)", 11.5, "end", op=0.95))
    # prediction vs measurement
    s.append(P([(165, 112), (165, 215), (358, 215)], BLUE, 1.8)); s.append(T(173, 168, "B1 / A5: predicted native", 11.5, op=0.95)); s.append(T(173, 183, "vs measured accessibility", 11.5, op=0.95)); s.append(T(173, 198, "(rho 0.40)", 11.5, op=0.95))
    s.append(P([(110, 338), (110, 255), (358, 255)], BLUE, 1.8)); s.append(T(118, 232, "C7: predicted native level", 11.5, op=0.95)); s.append(T(118, 247, "vs MPRA activity (rho 0.2)", 11.5, op=0.95))
    s.append(P([(250, 338), (250, 282), (358, 282)], ORANGE, 1.8)); s.append(T(262, 302, "C7: predicted shuffle effect", 11.5, op=0.95)); s.append(T(262, 317, "vs MPRA change (rho about 0)", 11.5, op=0.95))
    s.append(TT(500, 340, ["Lines = comparisons between two things.", "Arrows = one thing computed from / defining another.", "Dashed = the two ends are not independent."], 11.5, "start", lh=16, op=0.85))
    return fig("".join(s), 980, 450, "a6", "Which data sources support which comparisons",
               "Evidence map. Lines are comparisons between a prediction and a measurement (or between two measurements); arrows show what is computed from what. Footprint calls come from the same ATAC data used as the accessibility check, so footprint-score analyses (C4) are partly circular.")


# ---------------------------------------------------------------- D7: controls + titration
def d7_controls():
    s = []
    cx, sc = 490, 0.5  # 0.5 px per bp
    s.append(T(20, 22, "Where the shuffled region sits, relative to the motif (to scale, +-800 bp shown)", 12.5, weight="bold"))
    s.append(R(cx - 422 / 1.0, 40, 844, 22, "currentColor", "currentColor", 0.04, 1.2, rx=3))
    s.append(T(cx - 200, 56, "snATAC peak (median 0.85 kb)", 11.5, "middle", op=0.8))
    s.append(R(cx - 3, 36, 6, 30, ORANGE, ORANGE, 0.9, 1, rx=1)); s.append(T(cx, 80, "motif", 11.5, "middle", fill=ORANGE, weight="bold"))
    for off, lab in ((-600, "-600"), (-300, "-300"), (300, "+300"), (600, "+600")):
        x = cx + off * sc - 42.5
        s.append(R(x, 96, 85, 22, "currentColor", "currentColor", 0.12, 1.4, rx=3)); s.append(T(x + 42.5, 111, lab, 11.5, "middle"))
    s.append(R(cx - 42.5, 96, 85, 22, ORANGE, ORANGE, 0.35, 1.8, rx=3)); s.append(T(cx, 111, "centred", 11.5, "middle", weight="bold"))
    s.append(T(20, 111, "170 bp shuffles", 12, weight="bold"))
    s.append(TT(20, 136, ["Control question: is the loss specific to the AP-1 site, or does shuffling ANY 170 bp of the peak do this?", "Each is read out as ATAC at the motif; centred = the standard 170 bp shuffle, the others are matched-size controls (5 shuffles each)."], 11.5, "start", lh=15, op=0.9))
    s.append(T(20, 190, "Width titration: centred shuffles of increasing width", 12.5, weight="bold"))
    widths = (7, 11, 21, 41, 85, 170, 340, 680)
    for i, w in enumerate(widths):
        y = 204 + i * 20
        bw = max(w * sc, 3)
        s.append(R(cx - bw / 2, y, bw, 14, ORANGE, ORANGE, 0.15 + 0.07 * i, 1.4, rx=2))
        s.append(T(cx + bw / 2 + 8, y + 12, f"{w} bp" + ("  (core only)" if w == 7 else "  (whole motif)" if w == 11 else "  (the standard shuffle)" if w == 170 else ""), 11.5))
    s.append(R(cx - 3, 200, 6, 168, ORANGE, "none", 0, 0, rx=0)); s.append(L(cx, 196, cx, 366, ORANGE, 1, dash="3 3"))
    s.append(TT(20, 236, ["Question: how does the loss grow", "with the shuffled width, and", "how much does the core alone", "account for?", "", "Each width is read out as ATAC", "at the motif, the 170 bp and", "the peak; 5 shuffles each."], 11.5, "start", lh=15, op=0.9))
    return fig("".join(s), 980, 384, "a7", "Design of the matched-control shuffles and the width titration",
               "Controls and titration. Top: same-size (170 bp) shuffles centred on the motif or shifted by -600, -300, +300 or +600 bp, which tests whether the effect is specific to the site. Bottom: centred shuffles of eight widths from the 7 bp core to 680 bp, which measures how the effect scales.")


# ---------------------------------------------------------------- D8: edit types on sequence
def d8_edits():
    s = []
    native = "ACGTTCCTGAGTCACCGTAGA"  # motif = index 5..15 (CCTGAGTCACC), core = 7..13 (TGAGTCA)
    rows = [
        ("native", native, None, "unedited reference"),
        ("shuffle", "ACGTTCGCACTACTGCGTAGA", range(5, 16), "shuffle arm: letters permuted, composition kept"),
        ("random", "ACGTTTAGCGTAACTAGTAGA", range(5, 16), "mutation arm: every base replaced by a different one"),
        ("dinucleotide", "ACGTTCCACTGAGTCCGTAGA", range(5, 16), "mutation arm: shuffled, dinucleotide counts kept"),
        ("core", "ACGTTCCCATCAGGCCGTAGA", range(7, 14), "mutation arm: 7 bp core replaced (no AP-1 half-site left)"),
        ("flank (motif)", "ACGTTTATGAGTCAGTGTAGA", [5, 6, 14, 15], "mutation arm: 4 motif flank bases mutated, core kept"),
        ("flank (170 bp)", "GACACCCTGAGTCACCACTAC", list(range(0, 5)) + list(range(16, 21)), "mutation arm: surroundings mutated, whole motif kept"),
    ]
    x0, lw = 190, 26
    s.append(T(20, 24, "The same 21 bp excerpt (illustrative letters) under each edit; orange = changed bases", 12.5, weight="bold"))
    # motif and core brackets over the native row
    s.append(L(x0 + 5 * lw, 52, x0 + 16 * lw - 2, 52, "currentColor", 1.4)); s.append(T(x0 + 10.5 * lw, 46, "11 bp motif", 11.5, "middle"))
    s.append(L(x0 + 7 * lw, 62, x0 + 14 * lw - 2, 62, ORANGE, 2)); s.append(T(x0 + 10.5 * lw, 76, "7 bp core TGA[C/G]TCA", 11, "middle", fill=ORANGE, weight="bold"))
    for r, (name, seq, chg, desc) in enumerate(rows):
        y = 88 + r * 34
        s.append(T(20, y + 18, name, 12.5, weight="bold" if r == 0 else "normal"))
        for i, ch in enumerate(seq):
            ch_ = chg is not None and i in chg
            s.append(R(x0 + i * lw, y, lw - 2, 26, ORANGE if ch_ else "currentColor", ORANGE if ch_ else "currentColor", 0.40 if ch_ else 0.06, 1.2, rx=3))
            s.append(T(x0 + i * lw + (lw - 2) / 2, y + 18, ch, 13, "middle", weight="bold"))
        s.append(T(x0 + 21 * lw + 14, y + 18, desc, 11.5, op=0.9))
    return fig("".join(s), 1160, 340, "a8", "What each edit does to the sequence",
               "Edit types on one illustrative excerpt: the motif (11 bp, with its 7 bp core) and 5 bp of surrounding sequence either side. Shuffles keep the letters, random edits replace every base, dinucleotide shuffles keep adjacent pairs, the core and flank edits change only part of the motif, and the 170 bp flank edit changes everything but the motif. At the 170 bp and peak widths the same edits apply to those regions.")


# ---------------------------------------------------------------- D9: neutral-context MPRA
def d9_neutral():
    s = []

    def lane(y, col, title, box1, box2, box3):
        o = T(20, y + 2, title, 13, weight="bold", fill=col)
        o += R(20, y + 14, 250, 70, col, col, 0.08, 1.8) + TT(145, y + 40, box1, 12, lh=15)
        o += L(270, y + 49, 330, y + 49, col, 1.8, arrow="a9" + ("b" if col == BLUE else "o"))
        o += R(330, y + 14, 150, 70, "currentColor", "currentColor", 0.05, 1.8) + TT(405, y + 42, box2, 12, lh=16, weight="bold")
        o += L(480, y + 49, 540, y + 49, col, 1.8, arrow="a9" + ("b" if col == BLUE else "o"))
        o += R(540, y + 14, 210, 70, col, "none", 0, 1.8) + TT(645, y + 40, box3, 12, lh=15)
        return o

    s.append(lane(24, BLUE, "Genomic context (done earlier)", ["the 170 bp element sits in its own", "131 kb mouse window"], ["AlphaGenome", "predicted ATAC"], ["accessibility over the", "element, native and", "scrambled copies"]))
    s.append(lane(150, ORANGE, "Neutral context (overnight)", ["the same 170 bp element is placed at", "the centre of a dinucleotide-shuffled", "random 131 kb window (2 backgrounds)"], ["AlphaGenome", "predicted ATAC"], ["accessibility over the", "insert, native and", "scrambled copies"]))
    s.append(R(790, 40, 190, 160, GREEN, GREEN, 0.10, 1.8)); s.append(TT(885, 66, ["MPRA (measured)", "fibroblast, episomal"], 12.5, weight="bold", lh=16))
    s.append(TT(885, 108, ["native activity", "log2 alpha", "native - scrambled", "(motif / whole)"], 12, lh=16))
    s.append(P([(750, 73), (788, 73)], BLUE, 1.6, arrow="a9b")); s.append(P([(750, 199), (770, 199), (770, 150), (788, 150)], ORANGE, 1.6, arrow="a9o"))
    s.append(TT(20, 270, ["Question: the native ATAC head did not track the MPRA changes. Was that because the model also sees the genomic", "neighbourhood, which a plasmid reporter does not have? Placing the element in a neutral window removes the neighbourhood."], 12, "start", lh=17, op=0.9))
    s.append(TT(20, 320, ["Result: native levels still correlate weakly with the MPRA (rho 0.1 to 0.3) and the shuffle effects still do not (rho about 0)", "in either context, so context is not the explanation."], 12, "start", lh=17, weight="bold"))
    return fig("".join(s), 980, 345, "a9", "Design of the genomic-context versus neutral-context MPRA comparison",
               "Each MPRA element is predicted twice: in its own genomic window and inserted into a neutral background, then compared with the measured MPRA activity and with the MPRA change on scrambling. If genomic context were responsible for the earlier mismatch, the neutral-context prediction should agree better.")


def part_d_diagrams():
    return d7_controls(), d8_edits(), d9_neutral()


def design_section():
    return f"""
<h2 id="design">Design at a glance</h2>
<p>Six diagrams explain what was done and what each number means. Blue is the unedited (native) sequence or Stable sites, orange the edited sequence or Dynamic sites, green measured data.</p>
<h3>1. From a site to a number</h3>{d1_pipeline()}
<h3>2. What a window contains, and which edit and readout belong to which region</h3>{d2_window()}
<h3>3. The four experiments</h3>{d3_experiments()}
<h3>4. How the site labels are defined</h3>{d4_classes()}
<h3>5. How an effect is calculated and compared</h3>{d5_example()}
<h3>6. Where each comparison comes from, and which ones are independent</h3>{d6_evidence()}
"""
