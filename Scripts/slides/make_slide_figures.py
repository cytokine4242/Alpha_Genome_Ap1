#!/usr/bin/env python3
"""Presentation figures: simplified, purpose-built versions of the analysis plots, plus method schematics.
The dense analysis figures stay in the reports; these are built for a projector — few series, large type, one message each.
Writes into Analysis/slide_figures/.   python3 Scripts/slides/make_slide_figures.py"""
from __future__ import annotations

import glob
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.stats import binomtest, fisher_exact, mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Analysis" / "slide_figures"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(50)
# Same visual language as the analysis figures in the HTML reports: an off-white panel, a recessive grid, muted axis
# furniture, ink only where the data is. Type is larger here because these are read from a projector, not a page.
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df"
BLUE, ORANGE, GREY, GREEN = "#1F6FB2", "#C2551C", "#9a9993", "#1baf7a"
STABLE_GREY, GREEN_D, PURPLE_D = "#4F4F4F", "#2E7D32", "#7B3FA0"
CTRL_GREY = "#B3B2AE"          # the matched controls sit behind the classes, so they are the lighter mark
plt.rcParams.update({"font.size": 12.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7, "axes.axisbelow": True,
                     "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]


def boot(v, n=600):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5]) if len(v) > 5 else (np.nan, np.nan)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor=SURF)
    plt.close(fig)
    print("  wrote", name)


# ---------------------------------------------------------------- shared tables
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
SCF = pd.read_csv(ROOT / "Analysis/scan_features/tables/scan_features_per_site.csv").set_index("site_id")
S = sl.join(scl[["age_score_class"]]).join(meas.rename("measured"))
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))


def pairs(frame, col, a_, b_):
    m_ = frame[frame[col].isin([a_, b_])].dropna(subset=MCOLS).copy()
    Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
    avail, out = set(m_.index[m_[col] == b_]), []
    for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < 0.5:
            out.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(out, columns=["case", "control"])


# ================================================================ 1. method schematic
def eyebrow(ax, x, y, text):
    """The small letterspaced caps the report uses to label a panel without competing with the slide title."""
    ax.text(x, y, text.upper(), fontsize=11, color=MUTED, fontweight="bold", va="bottom", ha="left",
            fontstretch="expanded", transform=ax.transAxes)


SEQ, RAIL = "#d8d7d2", "#3f3e3b"

# ---------------------------------------------------------------- the five classes, and one matched control set each
GROUPS = [("Dynamic (all)", "stable_dynamic", "Dynamic", "#4F4F4F"), ("Age Up", "age_class", "Age Up", ORANGE),
          ("Age Down", "age_class", "Age Down", BLUE), ("Exercise Up", "ex_class", "Exercise Up", GREEN_D),
          ("Exercise Down", "ex_class", "Exercise Down", PURPLE_D)]
_STABLE = set(S.index[S.stable_dynamic == "Stable"])
_PAIRS = {}


def class_pairs(nm):
    """Each class matched to Stable on its own, so a gain and a loss can never cancel and every class gets
    controls it could actually be paired with. Cached: the matching is the slow part."""
    if nm not in _PAIRS:
        col, val = next((c, v) for n, c, v, _ in GROUPS if n == nm)
        f = S.copy()
        f["grp"] = np.where(f.index.isin(_STABLE), "Stable", np.where(f[col] == val, nm, ""))
        _PAIRS[nm] = pairs(f[f.grp != ""], "grp", nm, "Stable")
    return _PAIRS[nm]


KEY_DOT = ("grey line and band: all 2,907 Stable sites, median with a 95% bootstrap CI.   "
           "Each class is every site carrying that call, tested against that reference by Mann\u2013Whitney, "
           "BH across all {n} tests on this figure.   *** p<0.001   ** p<0.01   * p<0.05   ns \u2265 0.05")
KEY_PROP = ("grey line and band: all 2,907 Stable sites with a 95% bootstrap CI.   Fisher's exact test against that "
            "reference, BH across all {n} tests on this figure.   *** p<0.001   ** p<0.01   * p<0.05   ns \u2265 0.05")


PAN = OUT / "panels"
PAN.mkdir(parents=True, exist_ok=True)


def emit_class_figure(name, measures, panel_fn, test_note):
    """Emit one bare PNG per measure x mode plus a JSON manifest. Nothing but axes, ticks, axis labels and the
    significance marks goes in the image; every title, legend and caption is left for the slide to place as real text,
    so a wording change never means regenerating a figure."""
    MODES = ("all",)                      # accessibility is part of the biology, so the deck shows the uncontrolled view
    cols, panels, made, allps = [], [], [], []
    for ci, (v, title, sub, ylab) in enumerate(measures):
        cols.append({"title": title, "subtitle": sub})
        rowfiles = []
        for ri, mode in enumerate(MODES):
            fig, ax = plt.subplots(figsize=(3.3, 3.0))
            allps.append(panel_fn(ax, v, ylab, mode))
            f = f"{name}_c{ci}_r{ri}.png"
            made.append((fig, ax, f))
            rowfiles.append(f)
        panels.append(rowfiles)
    mark_fdr([a for _, a, _ in made], allps)
    for fig, _ax, f in made:
        fig.subplots_adjust(left=0.255, right=0.99, top=0.965, bottom=0.28)
        fig.savefig(PAN / f, facecolor=SURF)
        plt.close(fig)
    man = {"name": name, "cols": cols,
           "rows": ["every site of that class against all Stable sites"],
           "panels": panels,
           "legend": [[nm, col] for nm, _, _, col in GROUPS] + [["Stable (all 2,907 sites)", CTRL_GREY]],
           "key": test_note.format(n=sum(len(d) for d in allps))}
    (OUT / f"{name}.json").write_text(json.dumps(man, indent=1))
    print("  wrote", name, f"({sum(len(r) for r in panels)} panels)")


def class_dot_panel(ax, v, ylab, mode, minn=15, horizontal=False):
    """One measure, one mode. mode="all" compares every class against all Stable sites, nothing controlled.
    horizontal=True puts the classes on the y-axis, which stacks better in a narrow column."""
    ps = {}
    base_all = v.reindex(sorted(_STABLE)).dropna().values
    lo0, hi0 = boot(base_all, 400)
    span, line = (ax.axvspan, ax.axvline) if horizontal else (ax.axhspan, ax.axhline)
    if mode == "all":
        span(lo0, hi0, color=CTRL_GREY, alpha=0.30, lw=0, zorder=1)
        line(np.median(base_all), color=STABLE_GREY, lw=1.4, zorder=2)
    for k, (nm, col_, val_, col) in enumerate(GROUPS):
        if mode == "all":
            a = v.reindex(sorted(S.index[S[col_] == val_])).dropna().values
            if len(a) < minn:
                continue
            ps[nm] = mannwhitneyu(a, base_all).pvalue
            series = [(a, col, 0.0)]
        else:
            P = class_pairs(nm)
            x_, y_ = v.reindex(P.case).values, v.reindex(P.control).values
            ok = ~(np.isnan(x_) | np.isnan(y_))
            if ok.sum() < minn:
                continue
            a, b = x_[ok], y_[ok]
            ps[nm] = wilcoxon(a - b).pvalue
            series = [(b, CTRL_GREY, 0.19), (a, col, -0.19)]
        for vals, c_, off in series:
            lo, hi = boot(vals, 400)
            m_ = np.median(vals)
            if horizontal:
                ax.errorbar(m_, k + off, xerr=[[m_ - lo], [hi - m_]], fmt="o", ms=8, color=c_, lw=2.0, capsize=0, zorder=4)
            else:
                ax.errorbar(k + off, m_, yerr=[[m_ - lo], [hi - m_]], fmt="o", ms=9, color=c_, lw=2.2, capsize=0, zorder=4)
    labs = [n.replace("Exercise ", "Ex ").replace(" (all)", "") for n, _, _, _ in GROUPS]
    if horizontal:
        ax.set_yticks(range(len(GROUPS)), labs, fontsize=9.5)
        ax.set_ylim(len(GROUPS) - 0.45, -0.55)                   # Dynamic at the top, reading downward
        lo_, hi_ = ax.get_xlim()
        ax.set_xlim(lo_ - (hi_ - lo_) * 0.04, hi_ + (hi_ - lo_) * 0.12)
        if ylab:
            ax.set_xlabel(ylab, fontsize=9.5)
    else:
        lo_, hi_ = ax.get_ylim()
        ax.set_ylim(lo_, hi_ + (hi_ - lo_) * 0.16)
        ax.set_xticks(range(len(GROUPS)), labs, fontsize=10.5, rotation=20, ha="right")
        ax.set_xlim(-0.6, len(GROUPS) - 0.4)
        if ylab:
            ax.set_ylabel(ylab, fontsize=11)
    return ps


def class_prop_panel(ax, flag, ylab, mode, minn=15):
    """Same layout as class_dot_panel but for a yes/no site property: proportions with bootstrap CIs,
    Fisher on all sites and McNemar on the matched pairs, since there the two arms are paired."""
    ps = {}
    base_all = flag.reindex(sorted(_STABLE)).dropna().values
    if mode == "all":
        lo, hi = np.percentile([rng.choice(base_all, len(base_all)).mean() for _ in range(600)], [2.5, 97.5])
        ax.axhspan(lo * 100, hi * 100, color=CTRL_GREY, alpha=0.30, lw=0, zorder=1)
        ax.axhline(base_all.mean() * 100, color=STABLE_GREY, lw=1.4, zorder=2)
    for k, (nm, col_, val_, col) in enumerate(GROUPS):
        if mode == "all":
            a = flag.reindex(sorted(S.index[S[col_] == val_])).dropna().values
            b = base_all
            if len(a) < minn:
                continue
            ps[nm] = fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]])[1]
            series = [(a, col, 0.0)]
        else:
            P = class_pairs(nm)
            x_, y_ = flag.reindex(P.case).values, flag.reindex(P.control).values
            ok = ~(np.isnan(x_) | np.isnan(y_))
            if ok.sum() < minn:
                continue
            a, b = x_[ok], y_[ok]
            b01, b10 = int(((a == 1) & (b == 0)).sum()), int(((a == 0) & (b == 1)).sum())
            ps[nm] = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 else 1.0   # exact McNemar
            series = [(b, CTRL_GREY, 0.19), (a, col, -0.19)]
        for vals, c_, off in series:
            lo, hi = np.percentile([rng.choice(vals, len(vals)).mean() for _ in range(600)], [2.5, 97.5])
            ax.errorbar(k + off, vals.mean() * 100, yerr=[[(vals.mean() - lo) * 100], [(hi - vals.mean()) * 100]],
                        fmt="o", ms=9, color=c_, lw=2.2, capsize=0, zorder=4)
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(max(0, lo_), hi_ + (hi_ - lo_) * 0.16)
    ax.set_xticks(range(len(GROUPS)), [n.replace("Exercise ", "Ex ").replace(" (all)", "") for n, _, _, _ in GROUPS],
                  fontsize=10.5, rotation=20, ha="right")
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    if ylab:
        ax.set_ylabel(ylab, fontsize=11)
    return ps


TEST_KEY = ("grey line and band in the top row: the median of all Stable sites with its 95% bootstrap CI.   "
            "grey dots in the bottom row: that class's own matched Stable controls.\n"
            "stars test the class against that reference \u2014 top row Mann\u2013Whitney on all sites, bottom row paired "
            "Wilcoxon on matched pairs \u2014 BH-corrected across all {n} tests on this figure.   "
            "*** p < 0.001    ** p < 0.01    * p < 0.05    ns \u2265 0.05")


def mark_fdr(axes, allps):
    """BH across every class-by-measure test on the figure, drawn in each panel's top margin."""
    flat = [(i, nm, p) for i, d in enumerate(allps) for nm, p in d.items()]
    order = np.argsort([f[2] for f in flat])
    fdr = np.empty(len(flat))
    fdr[order] = np.minimum.accumulate(
        (np.array([flat[j][2] for j in order]) * len(flat) / np.arange(1, len(flat) + 1))[::-1])[::-1]
    for (i, nm, _p), q in zip(flat, fdr):
        k = [n for n, _, _, _ in GROUPS].index(nm)
        ax = axes[i]
        star = "***" if q < 1e-3 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"
        col_ = INK if q < 0.05 else MUTED
        bold = "bold" if q < 0.05 else "normal"
        if getattr(ax, "_horizontal", False):
            ax.text(0.985, k, star, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=10.5, color=col_, fontweight=bold)
        else:
            ax.text(k, 0.985, star, transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=12, color=col_, fontweight=bold)



def element_panel():
    """The setup slide: what one element looks like and what gets read off it. Nothing is edited here yet."""
    fig, ax = plt.subplots(figsize=(13.2, 4.3))
    ax.set_xlim(-345, 345)
    ax.set_ylim(-0.72, 1.12)
    ax.axis("off")
    x = np.linspace(-345, 345, 600)
    peak = 0.74 * np.exp(-(x ** 2) / (2 * 155 ** 2))
    ax.fill_between(x, 0.02, peak, color=BLUE, alpha=0.13, lw=0, zorder=1)
    ax.plot(x, peak, color=BLUE, lw=2.2, zorder=3, solid_capstyle="round")
    ax.text(-345, 0.90, "a cardiomyocyte ATAC peak, with an AP-1 footprint in it", fontsize=13.5, color=INK, va="bottom")
    ax.annotate("", xy=(-85, 0.22), xytext=(85, 0.22), arrowprops=dict(arrowstyle="|-|,widthA=0.4,widthB=0.4", color=MUTED, lw=1.3))
    ax.text(0, 0.265, "170 bp \u2014 the window every effect is read over", ha="center", fontsize=12, color=MUTED)
    ax.plot([-345, 345], [0.02, 0.02], color=RAIL, lw=1.5, zorder=4)
    for xx, lab, c, h in ((-215, "MEF2", GREY, 15), (0, "AP-1", ORANGE, 11), (120, "GATA", GREY, 14)):
        ax.add_patch(Rectangle((xx - h / 2 - 5, -0.085), h + 10, 0.21, color=c, ec="none", zorder=5))
        ax.text(xx, -0.15, lab, ha="center", va="top", fontsize=12.5, color=c if c == ORANGE else MUTED,
                fontweight="bold" if c == ORANGE else "normal")
    for xx in (-300, -150, 0, 150, 300):
        ax.plot([xx, xx], [0.02, -0.02], color=MUTED, lw=1, zorder=4)
        ax.text(xx, -0.33, f"{xx:+d} bp" if xx else "0", ha="center", va="top", fontsize=11, color=MUTED)
    ax.text(0, -0.56, "AlphaGenome reads 131 kb of sequence and predicts accessibility.\n"
                      "We change the middle of it and read what the prediction loses.",
            ha="center", va="top", fontsize=13, color=INK, linespacing=1.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    save(fig, "s0_element")


def edits_diagram():
    """The four edits of slide 4, drawn to scale on one rail: what is replaced, and how much of it."""
    fig, ax = plt.subplots(figsize=(7.0, 2.0))
    ax.set_xlim(-620, 620)
    ax.set_ylim(-0.42, 4.25)
    ax.axis("off")
    ROWS = [(3.55, 11, 0, "Replace the 11 bp motif", ORANGE),
            (2.55, 170, 0, "Replace the 170 bp element", ORANGE),
            (1.55, 900, 0, "Replace the whole peak", ORANGE),
            (0.55, 170, 300, "Replace 170 bp, 300 bp away \u2014 the motif untouched", STABLE_GREY)]
    for y, width, centre, lab, c in ROWS:
        ax.plot([-600, 600], [y, y], color=RAIL, lw=1.1, zorder=2)
        ax.add_patch(Rectangle((centre - width / 2, y - 0.13), max(width, 14), 0.26, color=c, ec="none", zorder=3))
        if centre:                                              # only the control spares the motif, so only it shows one
            ax.add_patch(Rectangle((-5.5, y - 0.11), 11, 0.22, color="#C9961F", ec=INK, lw=0.8, zorder=4))
        ax.text(-600, y + 0.26, lab, fontsize=10, color=INK, va="bottom",
                fontweight="bold" if c is ORANGE else "normal")
    for x_ in (-600, -300, 0, 300, 600):
        ax.plot([x_, x_], [0.18, 0.10], color=MUTED, lw=0.8)
        ax.text(x_, 0.02, f"{x_:+d}" if x_ else "0", ha="center", va="top", fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(PAN / "d_edits.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)
    print("  wrote d_edits")


def _band(kind):
    """A short, wide strip introducing one experiment, meant to sit under a slide title above that experiment's result."""
    fig, ax = plt.subplots(figsize=(10.4, 1.30))
    ax.set_xlim(0, 10.4)
    ax.set_ylim(0, 1.30)
    ax.axis("off")
    X0, W, Y = 3.35, 6.85, 0.62

    def head(title, sub):
        ax.text(0.0, 0.86, title, fontsize=13, fontweight="bold", va="center")
        ax.text(0.0, 0.42, sub, fontsize=11, color=MUTED, va="center")

    if kind == "block":
        head("Replace a whole region", "the motif, the 170 bp element or the peak,\neach swapped for random bases")
        w_ = 4.55                                              # three stacked rails, so the three widths can be compared
        for k, (frac, lab) in enumerate(((0.035, "the 11 bp motif"), (0.34, "the 170 bp element"), (1.0, "the whole peak"))):
            yy = 0.98 - k * 0.33
            ax.plot([X0, X0 + w_], [yy, yy], color=RAIL, lw=1.2, zorder=2)
            bw = max(w_ * frac, 0.10)
            ax.add_patch(Rectangle((X0 + w_ / 2 - bw / 2, yy - 0.085), bw, 0.17, color=ORANGE, ec="none", zorder=3))
            ax.text(X0 + w_ + 0.18, yy, lab, ha="left", va="center", fontsize=11, color=INK)
    elif kind == "tile":
        head("Sliding window", "every 10 bp from \u2212300 to +300 bp,\nreplaced by random bases in turn")
        ax.plot([X0, X0 + W], [Y, Y], color=RAIL, lw=1.3, zorder=2)
        for k in range(12):
            xx = X0 + k * W / 12
            ax.add_patch(Rectangle((xx + 0.015, Y - 0.12), W / 12 - 0.03, 0.24,
                                   color=ORANGE if k == 5 else SEQ, ec="none", zorder=3))
        ax.annotate("", xy=(X0 + W * 0.63, Y + 0.33), xytext=(X0 + W * 0.45, Y + 0.33),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.5, mutation_scale=13))
        ax.text(X0 + W * 0.43, Y + 0.33, "one window at a time  ", ha="right", va="center", fontsize=10.5, color=ORANGE)
    elif kind == "base":
        head("Single bases", "all three alternatives at every base,\n\u00b130 bp around the motif")
        letters = "C A G T G A C T C A G"
        for k, ch in enumerate(letters.split(" ")):
            xx = X0 + (k + 0.5) * W / 11
            ax.text(xx, Y + 0.10, ch, ha="center", va="center", fontsize=12.5, family="monospace",
                    color=ORANGE if k == 5 else MUTED, fontweight="bold" if k == 5 else "normal", zorder=4)
        xx = X0 + 5.5 * W / 11
        ax.add_patch(Rectangle((xx - W / 23, Y - 0.04), W / 11.5, 0.28, fc="none", ec=ORANGE, lw=1.4, zorder=3))
        for jj, alt in enumerate("ACT"):
            ax.text(xx + (jj - 1) * 0.26, Y - 0.32, alt, ha="center", va="center", fontsize=11.5,
                    family="monospace", color=ORANGE, zorder=4)
        ax.annotate("", xy=(xx, Y - 0.22), xytext=(xx, Y - 0.08), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.2, mutation_scale=11))
        ax.text(X0 + W + 0.06, Y + 0.10, "\u00d7 every position", ha="left", va="center", fontsize=10.5, color=MUTED)
    else:
        head("Motif titration", "1\u20137 core bases changed, and the motif\nmoved, duplicated or removed")
        ax.text(X0, Y + 0.30, "TGACTCA", fontsize=12, family="monospace", color=MUTED, va="center")
        ax.annotate("", xy=(X0 + 0.92, Y + 0.30), xytext=(X0 + 0.70, Y + 0.30),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.1, mutation_scale=11))
        for k, ch in enumerate("TGAgTCA"):
            ax.text(X0 + 1.02 + k * 0.125, Y + 0.30, ch, fontsize=12, family="monospace",
                    color=ORANGE if ch.islower() else MUTED, fontweight="bold" if ch.islower() else "normal", va="center")
        ax.text(X0 + 2.05, Y + 0.30, "   1 of 7 bases", fontsize=10.5, color=MUTED, va="center")
        y2 = Y - 0.26
        ax.plot([X0, X0 + W], [y2, y2], color=RAIL, lw=1.3, zorder=2)
        mw = W * 0.085
        ax.add_patch(Rectangle((X0 + W * 0.30, y2 - 0.10), mw, 0.20, color=ORANGE, ec="none", zorder=3))
        for frac in (0.06, 0.58, 0.80):
            ax.add_patch(Rectangle((X0 + W * frac, y2 - 0.10), mw, 0.20, fc="none", ec=ORANGE, lw=1.2, ls=(0, (2.5, 1.8)), zorder=3))
        ax.annotate("", xy=(X0 + W * 0.58, y2 + 0.13), xytext=(X0 + W * 0.355, y2 + 0.13),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.4, mutation_scale=12, connectionstyle="arc3,rad=-0.5"))
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
    save(fig, f"d_{kind}")


def bands():
    for k in ("block", "tile", "base", "titr"):
        _band(k)


# ================================================================ 2. the matched comparison, raw and matched together
def matched_null():
    KEY = [("replace motif, read over motif", "Replace the 11 bp motif", "the motif alone", "Effect (log2FC)"),
           ("replace 170 bp, read over 170 bp", "Replace the 170 bp element", "the element around it", "Effect (log2FC)"),
           ("replace peak, read over peak", "Replace the whole peak", "everything, read over the peak", "Effect (log2FC)"),
           ("170 bp shuffle \u00b1300 bp away, read over motif", "Edit 300 bp away", "the control: nothing at stake", "Effect (log2FC)")]
    d = S.join(FD[[k for k, _, _, _ in KEY]])
    emit_class_figure("s2_matched_null", [(d[c], t, sub, y) for c, t, sub, y in KEY], class_dot_panel, KEY_DOT)


# ================================================================ 3. dependence shape, raw and matched
def _profile_schematic(ax):
    """What the measures are read off. The motif is 11 bp -- a hairline at this scale -- so it gets its own narrow
    marker and a leader line; the shaded band is the \u00b185 bp readout window and is labelled as such."""
    ax.axis("off")
    ax.set_xlim(-310, 310)
    ax.set_ylim(-1.95, 1.05)
    rng2 = np.random.default_rng(7)
    c = np.arange(-295, 300, 10)
    prof = -(0.95 * np.exp(-(c ** 2) / (2 * 55 ** 2)) + 0.30 * np.exp(-((c - 150) ** 2) / (2 * 40 ** 2))
             + 0.22 * np.exp(-((c + 190) ** 2) / (2 * 45 ** 2))) - rng2.uniform(0, 0.05, len(c))
    ax.axvspan(-85, 85, color=ORANGE, alpha=0.09, lw=0, zorder=1)
    ax.bar(c[np.abs(c) <= 85], prof[np.abs(c) <= 85], width=9, color=ORANGE, zorder=3)
    ax.bar(c[np.abs(c) > 85], prof[np.abs(c) > 85], width=9, color=CTRL_GREY, zorder=3)
    ax.axhline(0, color=RAIL, lw=1.2, zorder=4)
    ax.add_patch(Rectangle((-6, -0.07), 12, 0.14, color="#C9961F", ec=INK, lw=0.9, zorder=6))
    ax.annotate("AP-1 motif, 11 bp", xy=(6, 0.02), xytext=(118, 0.52), fontsize=9.5, color="#9a7d0a",
                arrowprops=dict(arrowstyle="->", color="#9a7d0a", lw=1.0))
    ax.text(0, 0.90, "one site's sliding-window profile", ha="center", fontsize=11, color=INK, fontweight="bold")
    ax.text(0, 0.66, "each bar = one 10 bp window replaced", ha="center", fontsize=9.5, color=MUTED)
    ax.annotate("", xy=(-85, -1.34), xytext=(85, -1.34), arrowprops=dict(arrowstyle="|-|,widthA=0.3,widthB=0.3", color=ORANGE, lw=1.2))
    ax.text(0, -1.42, "near: within \u00b185 bp of the motif", ha="center", va="top", fontsize=9.5, color=ORANGE, fontweight="bold")
    ax.text(-215, 0.30, "far", ha="center", fontsize=10, color=MUTED)
    ax.text(215, 0.30, "far", ha="center", fontsize=10, color=MUTED)
    ax.annotate("", xy=(0, -1.78), xytext=(-185, -1.78), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.0))
    ax.text(-92, -1.74, "spread", ha="center", va="bottom", fontsize=9.5, color=MUTED)


def dependence():
    SCF2 = SCF.assign(near_abs_170=SCF.total_dependence_170 * SCF.share_within_85_170,
                      far_abs_170=SCF.total_dependence_170 * (1 - SCF.share_within_85_170))
    M = [("near_abs_170", "Dependence near the motif", "loss summed over windows within \u00b185 bp", "Summed loss within \u00b185 bp (log2)"),
         ("far_abs_170", "Dependence out in the context", "loss summed over windows beyond that", "Summed loss beyond \u00b185 bp (log2)"),
         ("spread_170", "Spread along the element", "SD of window position, weighted by loss", "Spread of the loss (bp)")]
    d = S.join(SCF2[[m[0] for m in M]], how="inner")
    cols, panels, made, allps = [], [], [], []
    for ci, (col, title, sub, xlab) in enumerate(M):
        cols.append({"title": title, "subtitle": sub})
        fig, ax = plt.subplots(figsize=(4.4, 1.40))
        ax._horizontal = True
        allps.append(class_dot_panel(ax, d[col], xlab, "all", horizontal=True))
        f = f"s3_dependence_c{ci}_r0.png"
        made.append((fig, ax, f))
        panels.append([f])
    mark_fdr([a_ for _, a_, _ in made], allps)
    for fig, _ax, f in made:
        fig.subplots_adjust(left=0.225, right=0.99, top=0.95, bottom=0.33)
        fig.savefig(PAN / f, facecolor=SURF)
        plt.close(fig)
    (OUT / "s3_dependence.json").write_text(json.dumps(
        {"name": "s3_dependence", "cols": cols, "rows": ["every site of that class against all Stable sites"],
         "panels": panels, "legend": [[nm, col] for nm, _, _, col in GROUPS] + [["Stable (all 2,907 sites)", CTRL_GREY]],
         "key": KEY_DOT.format(n=sum(len(x) for x in allps))}, indent=1))
    print("  wrote s3_dependence (3 stacked panels)")


# ================================================================ 4. titration, everything against native
def titration():
    tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/titrate*_cm_chunk*.csv")))], ignore_index=True)
    tt = tt[tt.trackset == "CM"]
    nat = tt[tt.kind == "native"].set_index("site_id").atac_170
    pc = 0.05 * nat.median()
    e = tt[tt.kind != "native"].copy()
    e["eff"] = np.log2((e.atac_170 + pc) / (e.site_id.map(nat) + pc))
    TF = pd.read_csv(ROOT / "Analysis/titration/tables/titration_features_per_site.csv").set_index("site_id")
    # the 1-base and 7-base results repeat the motif shuffle shown earlier, so only copy number and position are kept here
    M = [("add_copy_net_170", "A second copy added", "net of a shuffled-sequence control", "Effect (log2FC)"),
         ("relocate_net_170", "The motif moved elsewhere", "the motif destroyed and re-inserted, net of the same control", "Effect (log2FC)")]
    d = S.join(TF[[m[0] for m in M]], how="left")
    emit_class_figure("s4_titration", [(d[c], t, sub, y) for c, t, sub, y in M], class_dot_panel, KEY_DOT)
    CURVE = [("Stable", None, None, CTRL_GREY)] + [(n, c, v, col) for n, c, v, col in GROUPS if n != "Dynamic (all)"]
    members = lambda nm, c_, v_: set(S.index[S.stable_dynamic == "Stable"]) if nm == "Stable" else set(S.index[S[c_] == v_])
    for fname, kind, xlab, ylab, zero in (("d_addcopy", "dup", "Position of the added copy (bp)", "Change vs native (log2FC)", True),
                                          ("d_relocate", "move", "Position it was moved to (bp)", "Change vs native (log2FC)", True)):
        fig, ax = plt.subplots(figsize=(4.0, 2.8))
        for nm, c_, v_, col in CURVE:
            w = e[(e.kind == kind) & (e.site_id.isin(members(nm, c_, v_)))].groupby("pos").eff.median()
            ax.plot(w.index, w.values, color=col, lw=2.0, marker="o", ms=4)
        if zero:
            ax.axhline(0, color=INK, lw=1.0)
        ax.set_xlabel(xlab, fontsize=10.5)
        ax.set_ylabel(ylab, fontsize=10.5)
        fig.subplots_adjust(left=0.235, right=0.985, top=0.965, bottom=0.21)
        fig.savefig(PAN / f"{fname}.png", facecolor=SURF, bbox_inches="tight")
        plt.close(fig)


# ================================================================ 5. the age-gained result
def age_gained():
    """One bare panel: the three titration quantities on which age-gained (\u0394score) sites differ."""
    TI = pd.read_csv(ROOT / "Analysis/titration/tables/titration_by_class.csv")
    ti = lambda f: TI[(TI.contrast == "Age Up vs Unchanged (\u0394score)") & (TI.feature.str.startswith(f))].iloc[0]
    items = [("net value of a second", "Gains more from\nan added motif copy"),
             ("does a second copy help", "Benefits more when\nthat copy is near"),
             ("effect of changing 1 core", "Loses more from one\ndamaged core base")]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ys = np.arange(len(items))[::-1]
    for y, (key, lab) in zip(ys, items):
        r = ti(key)
        ax.scatter(r.median_reference, y, s=110, color=CTRL_GREY, zorder=3)
        ax.scatter(r.median_class, y, s=110, color=ORANGE, zorder=4)
        ax.plot([r.median_reference, r.median_class], [y, y], color=MUTED, lw=1.5, zorder=2)
        ax.text(max(r.median_class, r.median_reference) + 0.06, y,
                f"FDR {r.fdr:.3f}" if r.fdr >= 1e-3 else "FDR < 0.001", va="center", fontsize=10.5, color=MUTED)
    ax.set_yticks(ys, [l for _, l in items], fontsize=11)
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_xlim(-0.80, 0.80)
    ax.set_xlabel("Effect on predicted accessibility (log2FC)", fontsize=10.5)
    fig.subplots_adjust(left=0.265, right=0.99, top=0.97, bottom=0.165)
    fig.savefig(PAN / "s5_age_gained.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)
    print("  wrote s5 panel")


# ================================================================ 6. MPRA against the model, same axes
def mpra():
    """Two bare panels: the motif's share of local dependence in each assay, and the class direction, where they disagree."""
    FDl = FD["replace 170 bp, read over 170 bp"]
    agm = (FDl - FD["replace surroundings (motif kept), read over 170 bp"]) / FDl
    ag = pd.DataFrame({"share": agm.where(FDl < -0.5)}).join(sl[["stable_dynamic"]]).dropna()
    m = pd.read_csv(ROOT / "Data/MPRA_Results_with_TSS_distance.csv")
    m = m[m.replicated == True].copy()
    m["share"] = (m.delta_vs_MotifShuff / m.delta_vs_FullShuff).replace([np.inf, -np.inf], np.nan)
    mp = m[m.delta_vs_FullShuff > 0.1]

    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    for k, (v, c) in enumerate(((ag.share.values, BLUE), (mp.share.dropna().values, ORANGE))):
        v = v[(v > -0.5) & (v < 1.5)]
        lo, hi = boot(v)
        ax.bar(k, np.median(v), width=0.52, color=c, alpha=0.85)
        ax.errorbar(k, np.median(v), yerr=[[np.median(v) - lo], [hi - np.median(v)]], color=INK, lw=1.8, capsize=5)
        ax.text(k, np.median(v) + 0.028, f"{np.median(v):.0%}", ha="center", fontsize=15, fontweight="bold", color=c)
        ax.text(k, -0.026, f"n = {len(v):,}", ha="center", fontsize=10, color=MUTED)
    ax.set_xticks([0, 1], ["AlphaGenome\n(CM, in genome)", "MPRA\n(fibroblast)"], fontsize=10.5)
    ax.set_ylim(-0.05, 0.19)
    ax.set_ylabel("Motif's share of the element's\nsequence dependence", fontsize=10.5)
    fig.subplots_adjust(left=0.235, right=0.985, top=0.97, bottom=0.185)
    fig.savefig(PAN / "s6_share.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    agS = (-FDl.reindex(sl.index[sl.stable_dynamic == "Stable"]).dropna()).values
    agD = (-FDl.reindex(sl.index[sl.stable_dynamic == "Dynamic"]).dropna()).values
    mpS = m[m.Group == "Stable"].delta_vs_FullShuff.dropna().values
    mpD = m[m.Group != "Stable"].delta_vs_FullShuff.dropna().values
    for k, (vS, vD, c) in enumerate(((agS, agD, BLUE), (mpS, mpD, ORANGE))):
        r = np.median(vD) / np.median(vS)
        bs = [np.median(rng.choice(vD, len(vD))) / np.median(rng.choice(vS, len(vS))) for _ in range(800)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        ax.errorbar(k, r, yerr=[[max(r - max(lo, 0.4), 0)], [max(min(hi, 7.0) - r, 0)]], fmt="o", ms=13, color=c, lw=2.4, capsize=0, zorder=4)
        ax.text(k + 0.16, r, f"{r:.2f}\u00d7", ha="left", va="center", fontsize=14, fontweight="bold", color=c)
    ax.axhline(1, color=INK, lw=1.3, ls=":")
    ax.set_yscale("log")
    ax.set_yticks([0.5, 1, 2, 4], ["0.5\u00d7", "1\u00d7", "2\u00d7", "4\u00d7"])
    ax.set_ylim(0.38, 7.5)
    ax.set_xticks([0, 1], ["AlphaGenome\n(accessibility)", "MPRA\n(reporter activity)"], fontsize=10.5)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Responsive elements, relative to\nStable in the same assay", fontsize=10.5)
    fig.subplots_adjust(left=0.235, right=0.985, top=0.97, bottom=0.185)
    fig.savefig(PAN / "s6_direction.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)
    print("  wrote s6 panels")


# ================================================================ 7. spacing, with the test on it
def spacing():
    CL = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_peak.csv").set_index("family")
    FINE = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_profile.csv")
    IMP = pd.read_csv(ROOT / "Analysis/distance_rules/tables/g_importance_by_distance.csv")
    SER = (("MEF2", ORANGE, 2.6), ("KLF / SP", GREY, 2.2))

    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for fam, c, lw in SER:
        w = FINE[FINE.family == fam]
        ax.plot(w.distance, w.share.rolling(3, center=True, min_periods=1).mean() * 100, color=c, lw=lw)
    ax.axvspan(0, 5.5, color="#EDEDED", lw=0)
    ax.axvspan(11, 25, color=ORANGE, alpha=0.10, lw=0)
    ax.set_xlabel("Distance from the AP-1 motif (bp)", fontsize=10.5)
    ax.set_ylabel("Share of that family's footprints (%)", fontsize=10)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 2.1)
    fig.subplots_adjust(left=0.155, right=0.985, top=0.97, bottom=0.165)
    fig.savefig(PAN / "s7_occurrence.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for fam, c, lw in SER:
        w = IMP[IMP.family == fam]
        ax.plot(w.mid, w.importance, color=c, lw=lw, marker="o", ms=6)
    ax.axhline(0, color=MUTED, lw=1.0)
    ax.set_xlabel("Distance from the AP-1 motif (bp)", fontsize=10.5)
    ax.set_ylabel("Extra loss vs footprint-free tiles (log2FC)", fontsize=10)
    ax.set_xlim(-10, 270)
    fig.subplots_adjust(left=0.155, right=0.985, top=0.97, bottom=0.165)
    fig.savefig(PAN / "s7_importance.png", facecolor=SURF, bbox_inches="tight")
    plt.close(fig)
    print("  wrote s7 panels   MEF2 near/far ratio %.2f (FDR %.3f)" % (CL.loc["MEF2", "density_ratio_near_over_far"], CL.loc["MEF2", "fdr"]))


# ================================================================ 8. partners, ranked AND split by class
def partners():
    """Two bare panels: the ranking with every class on it, all sites and matched. Titles and key live on the slide."""
    FI = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    FI = FI[FI.family != "Other"]
    cl = S[["stable_dynamic", "age_class", "ex_class"]]
    stable_ids = set(cl.index[cl.stable_dynamic == "Stable"])
    FAM = ["MEF2", "NFI (Nfia/b/c/x)", "GATA", "TEAD", "AP-1 (Fos/Jun/Atf/Batf)", "KLF / SP"]
    matched = {nm: class_pairs(nm) for nm, _, _, _ in GROUPS if nm != "Dynamic (all)"}
    pooled_ctrl = set().union(*[set(P.control) for P in matched.values()])
    tests = []
    for mode, fname in (("all", "s8_partners_all.png"),):
        fig, ax = plt.subplots(figsize=(5.0, 3.6))
        ys = np.arange(len(FAM))[::-1]
        lim = 0.0
        for y, fam in zip(ys, FAM):
            g = FI[FI.family == fam].set_index("site_id").vs_free
            base = g.reindex(sorted(stable_ids if mode == "all" else pooled_ctrl)).dropna().values
            series = [("Stable", CTRL_GREY, base)]
            for nm, col_, val_, c in GROUPS:
                if nm == "Dynamic (all)":
                    continue
                ids = set(cl.index[cl[col_] == val_]) if mode == "all" else set(matched[nm].case)
                series.append((nm, c, g.reindex(sorted(ids)).dropna().values))
            for k, (nm, c, v) in enumerate(series):
                if len(v) < 18:
                    continue
                lo, hi = boot(v, 400)
                lim = min(lim, lo)
                ax.errorbar(np.median(v), y + (2 - k) * 0.155, xerr=[[np.median(v) - lo], [hi - np.median(v)]],
                            fmt="o", ms=6, color=c, lw=1.7, capsize=0, zorder=4)
                if nm != "Stable":
                    tests.append(mannwhitneyu(base, v).pvalue)
        ax.axvline(0, color=MUTED, lw=1.0, ls=":")
        ax.set_yticks(ys, [f.split(" (")[0].split(" /")[0] for f in FAM], fontsize=11)
        ax.set_xlim(lim * 1.05, 0.12)
        ax.set_ylim(-0.62, len(FAM) - 0.38)
        ax.set_xlabel("Loss when that tile is replaced (log2FC)", fontsize=10)
        fig.subplots_adjust(left=0.215, right=0.99, top=0.97, bottom=0.145)
        fig.savefig(PAN / fname, facecolor=SURF)
        plt.close(fig)
    q = np.sort(np.array(tests))
    q = float(np.min(q * len(q) / np.arange(1, len(q) + 1)))
    (OUT / "s8_partners.json").write_text(json.dumps({"n_tests": len(tests), "min_fdr": round(q, 3)}, indent=1))
    print(f"  wrote s8 panels   {len(tests)} tests, smallest FDR {q:.2f}")


# ================================================================ 9. the MEF2 deficit, by class, two ways
def mef2_by_class():
    FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
    FP = FP[~FP.own_ap1].copy()
    FP["mid"] = (FP.rel_start + FP.rel_end) / 2
    FI = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    imp = FI[FI.family == "MEF2"].set_index("site_id").vs_free
    near = set(FP[(FP.family == "MEF2") & (FP.mid.abs() <= 25)].site_id)
    leaned = {i for i in near if imp.get(i, 0) < -0.2}
    M = [(pd.Series({i: float(i in near) for i in S.index}), "A MEF2 footprint within 25 bp",
          "counting footprints, used by the model or not", "Sites with one (%)"),
         (pd.Series({i: float(i in leaned) for i in S.index}), "...one the model actually leans on",
          "restricted to tiles worth more than 0.2 log2", "Sites with one (%)")]
    emit_class_figure("s9_mef2_by_class", M, class_prop_panel, KEY_PROP)


# ================================================================ 10. presence vs rank, from the top-window odds
def top_tile_odds():
    """Bare forest panels. Is the MEF2 effect about having one, or about how it is ranked when there is one --
    asked against all Stable sites and against matched controls, since the two carry different information."""
    D = pd.read_csv(ROOT / "Analysis/top_tile_odds/tables/top3_tile_odds_by_class.csv")
    order = [nm for nm, _, _, _ in GROUPS]

    def forest(fname, rows, xlim):
        """rows: (filter dict, colour, offset) triples, each drawn as one series down the five classes."""
        fig, ax = plt.subplots(figsize=(4.3, 3.2))
        ys = np.arange(len(order))[::-1]
        for filt, col, off in rows:
            w = D[np.logical_and.reduce([D[k] == v for k, v in filt.items()])].set_index("class")
            for y, nm in zip(ys, order):
                if nm not in w.index:
                    continue
                r = w.loc[nm]
                lor = np.log2(r.odds_ratio)                      # log2, to match the fold changes used everywhere else
                ax.errorbar(lor, y + off, xerr=[[max(lor - np.log2(r.or_lo), 0)], [max(np.log2(r.or_hi) - lor, 0)]],
                            fmt="o", ms=7, color=col, lw=1.8, capsize=0, zorder=4)
        ax.axvline(0, color=INK, lw=1.2, ls=":")
        ax.set_xticks([t for t in (-2, -1, 0, 1, 2) if xlim[0] <= t <= xlim[1]])
        ax.set_xlim(*xlim)
        ax.set_yticks(ys, [n.replace("Exercise ", "Ex ").replace(" (all)", "") for n in order], fontsize=10.5)
        ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_xlabel("log2 odds ratio (95% CI)", fontsize=10.5)
        fig.subplots_adjust(left=0.265, right=0.985, top=0.97, bottom=0.165)
        fig.savefig(PAN / fname, facecolor=SURF, bbox_inches="tight")
        plt.close(fig)

    for mode, tag in (("all sites", "all"),):
        forest(f"s10_mef2_{tag}.png",
               [({"family": "MEF2", "mode": mode, "conditioning": "unconditional"}, ORANGE, -0.17),
                ({"family": "MEF2", "mode": mode, "conditioning": "given present"}, STABLE_GREY, 0.17)],
               (-2.6, 2.6))
    forest("s10_ap1.png",
           [({"family": "AP-1 (Fos/Jun/Atf/Batf)", "mode": "all sites", "conditioning": "unconditional"}, ORANGE, -0.17),
            ({"family": "AP-1 (Fos/Jun/Atf/Batf)", "mode": "all sites", "conditioning": "given present"}, STABLE_GREY, 0.17)],
           (-1.3, 1.8))
    print("  wrote s10 odds panels")


print("building slide figures")
element_panel(); bands(); edits_diagram(); matched_null(); dependence(); titration(); age_gained(); mpra(); spacing(); partners(); mef2_by_class(); top_tile_odds()
print("done ->", OUT)
