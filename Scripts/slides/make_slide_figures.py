#!/usr/bin/env python3
"""Presentation figures: simplified, purpose-built versions of the analysis plots, plus method schematics.
The dense analysis figures stay in the reports; these are built for a projector — few series, large type, one message each.
Writes into Analysis/slide_figures/.   python3 Scripts/slides/make_slide_figures.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon

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
    """Both classes on the same axis, all sites beside context-matched pairs \u2014 the idiom of the report's profile figure,
    cut down to the four edits that carry the argument."""
    KEY = [("replace motif, read over motif", "Replace the 11 bp motif"),
           ("replace 170 bp, read over 170 bp", "Replace the 170 bp element"),
           ("replace peak, read over peak", "Replace the whole peak"),
           ("170 bp shuffle \u00b1300 bp away, read over motif", "Edit 300 bp away  (control)")]
    d = S.join(FD[[k for k, _ in KEY]])
    P = pairs(d, "stable_dynamic", "Dynamic", "Stable")
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8), sharey=True)
    ys = np.arange(len(KEY))[::-1]
    for ax, mode in zip(axes, ("all", "matched")):
        for y, (col, _) in zip(ys, KEY):
            if mode == "all":
                vS, vD = d[d.stable_dynamic == "Stable"][col].dropna().values, d[d.stable_dynamic == "Dynamic"][col].dropna().values
                p = mannwhitneyu(vS, vD).pvalue
            else:
                a, b = d[col].reindex(P.case).values, d[col].reindex(P.control).values
                ok = ~(np.isnan(a) | np.isnan(b))
                vS, vD = b[ok], a[ok]
                p = wilcoxon(vD - vS).pvalue
            for v, c, off in ((vS, BLUE, 0.17), (vD, ORANGE, -0.17)):
                lo, hi = boot(v)
                ax.errorbar(np.median(v), y + off, xerr=[[np.median(v) - lo], [hi - np.median(v)]],
                            fmt="o", ms=10, color=c, lw=2.4, capsize=0, zorder=4)
            ax.text(0.30, y, stars(p), ha="left", va="center", fontsize=14,
                    color=INK if p < 0.05 else MUTED, fontweight="bold" if p < 0.05 else "normal")
        ax.axvline(0, color=MUTED, lw=1.0)
        ax.set_xlim(-2.25, 0.62)
        ax.set_xticks([-2.0, -1.5, -1.0, -0.5, 0.0])
        ax.set_ylim(-0.62, len(KEY) - 0.38)
    axes[0].set_yticks(ys, [l for _, l in KEY], fontsize=13.5)
    axes[0].set_title("All sites", fontsize=14.5, loc="left", pad=12)
    axes[1].set_title(f"Context-matched pairs  (n = {len(P)})", fontsize=14.5, loc="left", pad=12)
    h = [axes[0].plot([], [], "o", color=c, ms=10, label=l)[0] for c, l in ((BLUE, "Stable"), (ORANGE, "Dynamic"))]
    fig.subplots_adjust(left=0.205, right=0.985, top=0.90, bottom=0.255, wspace=0.09)
    fig.text(0.595, 0.155, "Median effect (log2FC);  more negative = more sensitive", ha="center", fontsize=13)
    fig.legend(handles=h, frameon=False, fontsize=13, loc="lower center", ncol=2, bbox_to_anchor=(0.595, 0.055))
    fig.text(0.595, 0.012, "each Dynamic site paired with a Stable site of the same measured and predicted accessibility, peak width and gene distance",
             ha="center", fontsize=11.5, color=MUTED)
    save(fig, "s2_matched_null")


# ================================================================ 3. dependence shape, raw and matched
def dependence():
    """Three measures of where a site's sequence dependence sits, each defined on the figure: a slide is no place for a
    name the audience has to take on trust."""
    M = [("share_within_85_170", "Concentrated near the motif",
          "of all the accessibility a site loses across the 60 sliding\nwindows, the fraction lost in windows within \u00b185 bp",
          "Share of the site's total loss\nfalling within \u00b185 bp", "{:.2f}"),
         ("spread_170", "Spread along the element",
          "the standard deviation of window position, weighting each\nwindow by how much of the loss it accounts for",
          "Spread of the loss (bp)", "{:.0f}"),
         ("ism_max", "Depth of the single best base",
          "the largest effect any single base change has, averaged\nover its three alternatives, anywhere in \u00b130 bp",
          "Largest single-base effect\n(|log2FC|)", "{:.2f}")]
    d = S.join(SCF[[m[0] for m in M]], how="inner")
    P = pairs(d, "stable_dynamic", "Dynamic", "Stable")
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 5.4))
    for ax, (col, lab, defn, ylab, _f) in zip(axes, M):
        for j_, mode in enumerate(("all sites", "matched")):
            if mode == "all sites":
                vs = [d[d.stable_dynamic == g][col].dropna().values for g in ("Stable", "Dynamic")]
            else:
                vs = [d[col].reindex(P.control).dropna().values, d[col].reindex(P.case).dropna().values]
            for i_, (v, c) in enumerate(zip(vs, (BLUE, ORANGE))):
                lo, hi = boot(v)
                ax.errorbar(j_ + (i_ - 0.5) * 0.26, np.median(v), yerr=[[np.median(v) - lo], [hi - np.median(v)]],
                            fmt="o", ms=11, color=c, lw=2.4, capsize=0)
        a, b = d[col].reindex(P.case).values, d[col].reindex(P.control).values
        ok = ~(np.isnan(a) | np.isnan(b))
        p_ = wilcoxon(a[ok] - b[ok]).pvalue
        ax.set_xticks([0, 1], ["all sites", "context-\nmatched"], fontsize=12)
        ax.set_xlim(-0.45, 1.45)
        ax.set_ylabel(ylab, fontsize=12.5)
        ax.set_title(lab, fontsize=13.5, loc="left", pad=30)
        ax.text(0.0, 1.015, defn, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5,
                color=MUTED, linespacing=1.45)
        ax.text(0.5, -0.195, f"matched p = {p_:.1g}" if p_ >= 1e-4 else "matched p < 0.0001",
                transform=ax.transAxes, ha="center", fontsize=12, color=MUTED)
    h = [axes[0].plot([], [], "o", color=c, ms=10, label=l)[0] for c, l in ((BLUE, "Stable"), (ORANGE, "Dynamic"))]
    fig.subplots_adjust(left=0.075, right=0.99, top=0.80, bottom=0.25, wspace=0.42)
    fig.legend(handles=h, frameon=False, fontsize=13, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.005))
    save(fig, "s3_dependence")


# ================================================================ 4. titration, everything against native
def titration():
    """Mechanism, but split by class throughout: the dose curve and the position curve are drawn for Stable and for
    Dynamic separately, and the third panel asks whether any titration feature separates the directional classes."""
    tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/titrate*_cm_chunk*.csv")))], ignore_index=True)
    tt = tt[tt.trackset == "CM"]
    nat = tt[tt.kind == "native"].set_index("site_id").atac_170
    pc = 0.05 * nat.median()
    e = tt[tt.kind != "native"].copy()
    e["eff"] = np.log2((e.atac_170 + pc) / (e.site_id.map(nat) + pc))
    e["cls"] = e.site_id.map(S.stable_dynamic)
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0), gridspec_kw={"width_ratios": [1, 1.15, 1.25]})

    ax = axes[0]
    for c_, col in (("Stable", STABLE_GREY), ("Dynamic", ORANGE)):
        mm = e[(e.kind == "core_mm") & (e.cls == c_)].groupby("pos").eff
        med = mm.median()
        lo = mm.apply(lambda v: boot(v.values, 300)[0])
        hi = mm.apply(lambda v: boot(v.values, 300)[1])
        ax.plot(med.index, med.values, color=col, lw=2.6, marker="o", ms=8, label=c_)
        ax.fill_between(med.index, lo.values, hi.values, color=col, alpha=0.16, lw=0)
    ax.set_xlabel("Core bases changed (of 7)", fontsize=12.5)
    ax.set_ylabel("Loss of predicted accessibility (log2FC)", fontsize=12)
    ax.set_title("One wrong base does almost all the damage", fontsize=13.5, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=11.5, loc="upper right")
    ax.text(0.5, -0.235, "the curve has the same shape in both classes", transform=ax.transAxes,
            ha="center", fontsize=11, color=MUTED)

    ax = axes[1]
    for kind, lab, ls in (("move", "motif destroyed, re-inserted elsewhere", "-"),
                          ("dup", "original kept, a second copy added", (0, (4, 2)))):
        for c_, col in (("Stable", STABLE_GREY), ("Dynamic", ORANGE)):
            w = e[(e.kind == kind) & (e.cls == c_)].groupby("pos").eff.median()
            ax.plot(w.index, w.values, color=col, lw=2.4, ls=ls, marker="o", ms=6)
    ax.axhline(0, color=INK, lw=1.2)
    ax.text(152, 0.03, "native", fontsize=11.5, color=INK, ha="right", va="bottom")
    ax.set_xlabel("Position of the inserted copy, relative to the motif (bp)", fontsize=12.5)
    ax.set_ylabel("Change vs native (log2FC)", fontsize=12)
    ax.set_title("Position is a dial \u2014 and the dial reads the same", fontsize=13.5, loc="left", pad=8)
    hh = [plt.Line2D([], [], color=STABLE_GREY, lw=2.4), plt.Line2D([], [], color=ORANGE, lw=2.4),
          plt.Line2D([], [], color=INK, lw=2.0, ls="-"), plt.Line2D([], [], color=INK, lw=2.0, ls=(0, (4, 2)))]
    ax.legend(hh, ["Stable", "Dynamic", "destroyed, re-inserted elsewhere", "original kept, a copy added"],
              frameon=False, fontsize=10, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=1, handlelength=2.4)
    ax.set_ylim(-0.92, 0.33)

    ax = axes[2]
    TI = pd.read_csv(ROOT / "Analysis/titration/tables/titration_by_class.csv")
    FEAT = [("effect of changing 1 core base", "one damaged\ncore base"),
            ("effect of changing all 7 core bases", "all seven\ncore bases"),
            ("effect of keeping only the first half-site", "half-site\nonly"),
            ("net value of a second motif copy", "a second\ncopy added"),
            ("net effect of relocating the motif", "motif moved\nelsewhere")]
    CON = [("Age Up vs Stable (threshold)", "Age Up", ORANGE), ("Age Down vs Stable (threshold)", "Age Down", BLUE),
           ("Ex Up vs Stable (threshold)", "Exercise Up", GREEN_D), ("Ex Down vs Stable (threshold)", "Exercise Down", PURPLE_D)]
    ys = np.arange(len(FEAT))[::-1]
    for y, (feat, lab) in zip(ys, FEAT):
        for k, (con, cl, col) in enumerate(CON):
            r = TI[(TI.contrast == con) & (TI.feature.str.startswith(feat))]
            if not len(r):
                continue
            r = r.iloc[0]
            off = (1.5 - k) * 0.19
            ax.scatter(r.mean_difference, y + off, s=58, color=col, zorder=4,
                       edgecolors=INK if r.fdr < 0.05 else "none", lw=1.4)
    ax.axvline(0, color=INK, lw=1.3)
    ax.set_yticks(ys, [l for _, l in FEAT], fontsize=11.5)
    ax.set_xlim(-0.30, 0.30)
    ax.set_ylim(-0.62, len(FEAT) - 0.38)
    ax.set_xlabel("Difference from that class's own\nmatched Stable controls (log2FC)", fontsize=12)
    ax.set_title("...and no directional class differs either", fontsize=13.5, loc="left", pad=8)
    ax.text(0.5, -0.235, "nothing survives correction (every FDR \u2265 0.08)", transform=ax.transAxes,
            ha="center", fontsize=11, color=MUTED)
    for _, cl, col in CON:
        ax.scatter([], [], s=58, color=col, label=cl)
    ax.legend(frameon=False, fontsize=10.5, loc="lower right", bbox_to_anchor=(1.02, -0.02))
    fig.subplots_adjust(left=0.065, right=0.99, top=0.90, bottom=0.20, wspace=0.37)
    save(fig, "s4_titration")


# ================================================================ 5. the age-gained result
def age_gained():
    TI = pd.read_csv(ROOT / "Analysis/titration/tables/titration_by_class.csv")
    ti = lambda f: TI[(TI.contrast == "Age Up vs Unchanged (Δscore)") & (TI.feature.str.startswith(f))].iloc[0]
    reg = pd.read_csv(ROOT / "Analysis/hypotheses/tables/hypothesis_register.csv").set_index("id")
    items = [("net value of a second", "Gains more from\nan added motif copy", "+0.38", "+0.32"),
             ("does a second copy help", "Benefits more when\nthat copy is near", "+0.21", "+0.14"),
             ("effect of changing 1 core", "Loses more from one\ndamaged core base", "−0.60", "−0.51")]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7), gridspec_kw={"width_ratios": [1.45, 1]})
    ax = axes[0]
    ys = np.arange(len(items))[::-1]
    for y, (key, lab, a_, b_) in zip(ys, items):
        r = ti(key)
        ax.scatter(r.median_reference, y, s=150, color=GREY, zorder=3, label="Unchanged" if y == ys[0] else None)
        ax.scatter(r.median_class, y, s=150, color=ORANGE, zorder=4, label="Age-gained (Δscore)" if y == ys[0] else None)
        ax.plot([r.median_reference, r.median_class], [y, y], color=MUTED, lw=1.6, zorder=2)
        ax.text(max(r.median_class, r.median_reference) + 0.075, y, f"FDR {r.fdr:.3f}" if r.fdr >= 1e-3 else "FDR < 0.001",
                va="center", fontsize=11.5, color=MUTED)
    ax.set_yticks(ys, [i[1] for i in items], fontsize=12.5)
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_xlim(-0.80, 0.78)
    ax.set_xlabel("Effect on predicted accessibility (log2FC)", fontsize=12.5)
    ax.set_title("Age-gained sites respond more to AP-1 dosage, both ways", fontsize=14, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=12, loc="upper left", bbox_to_anchor=(0.0, 0.42))
    ax = axes[1]
    vals = [84.5, 78.6]
    bars = ax.bar([0, 1], vals, width=0.58, color=[ORANGE, GREY])
    for x_, v in zip([0, 1], vals):
        ax.text(x_, v + 1.4, f"{v:.1f}%", ha="center", fontsize=16, fontweight="bold", color=ORANGE if x_ == 0 else MUTED)
    ax.set_xticks([0, 1], ["Age-gained\n(Δscore)", "Matched\nUnchanged"], fontsize=12.5)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Carry a perfect TGA(C/G)TCA core", fontsize=12.5)
    ax.set_title("...and they start from better motifs", fontsize=14, loc="left", pad=8)
    ax.text(0.5, -0.20, "p = 0.02   \u00b7   all 3,678 sites, read from the genome sequence", transform=ax.transAxes,
            ha="center", fontsize=11.5, color=MUTED)
    fig.tight_layout()
    save(fig, "s5_age_gained")


# ================================================================ 6. MPRA against the model, same axes
def mpra():
    MSA = pd.read_csv(ROOT / "Analysis/mpra_scan/tables/motif_share_both_assays.csv").set_index("dataset")
    AGD, MPD = "AlphaGenome (CM, genomic context)", "MPRA (fibroblast, episomal)"
    FDl = FD["replace 170 bp, read over 170 bp"]
    agm = (FDl - FD["replace surroundings (motif kept), read over 170 bp"]) / FDl
    ag = pd.DataFrame({"share": agm.where(FDl < -0.5)}).join(sl[["stable_dynamic"]]).dropna()
    m = pd.read_csv(ROOT / "Data/MPRA_Results_with_TSS_distance.csv")
    m = m[m.replicated == True].copy()
    m["share"] = (m.delta_vs_MotifShuff / m.delta_vs_FullShuff).replace([np.inf, -np.inf], np.nan)
    mp = m[m.delta_vs_FullShuff > 0.1]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.8), gridspec_kw={"width_ratios": [1, 1.2]})
    ax = axes[0]
    for k, (v, lab, c) in enumerate(((ag.share.values, "AlphaGenome\npredicted accessibility\n(cardiomyocyte, in genome)", BLUE),
                                     (mp.share.dropna().values, "MPRA\nreporter activity\n(fibroblast, episomal)", ORANGE))):
        v = v[(v > -0.5) & (v < 1.5)]
        lo, hi = boot(v)
        ax.bar(k, np.median(v), width=0.52, color=c, alpha=0.85)
        ax.errorbar(k, np.median(v), yerr=[[np.median(v) - lo], [hi - np.median(v)]], color=INK, lw=2, capsize=6)
        ax.text(k, np.median(v) + 0.035, f"{np.median(v):.0%}", ha="center", fontsize=19, fontweight="bold", color=c)
        ax.text(k, -0.028, f"n = {len(v):,}", ha="center", fontsize=11.5, color=MUTED)
    ax.set_xticks([0, 1], ["AlphaGenome\n(CM, in genome)", "MPRA\n(fibroblast, episomal)"], fontsize=12)
    ax.set_ylim(-0.05, 0.19)
    ax.set_ylabel("Share of the element's sequence dependence\nattributable to the AP-1 motif", fontsize=12)
    ax.set_title("Same answer: the motif is a small part", fontsize=14, loc="left", pad=8)
    ax = axes[1]
    agS = (-FDl.reindex(sl.index[sl.stable_dynamic == "Stable"]).dropna()).values
    agD = (-FDl.reindex(sl.index[sl.stable_dynamic == "Dynamic"]).dropna()).values
    mpS = m[m.Group == "Stable"].delta_vs_FullShuff.dropna().values
    mpD = m[m.Group != "Stable"].delta_vs_FullShuff.dropna().values
    panels = [("AlphaGenome\nloss of predicted\naccessibility", agS, agD, BLUE),
              ("MPRA\nloss of reporter\nactivity", mpS, mpD, ORANGE)]
    for k, (lab, vS, vD, c) in enumerate(panels):
        r = np.median(vD) / np.median(vS)
        bs = [np.median(rng.choice(vD, len(vD))) / np.median(rng.choice(vS, len(vS))) for _ in range(800)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        p_ = mannwhitneyu(vS, vD).pvalue
        lo, hi = max(lo, 0.4), min(hi, 7.0)
        ax.errorbar(k, r, yerr=[[max(r - lo, 0)], [max(hi - r, 0)]], fmt="o", ms=15, color=c, lw=2.8, capsize=0, zorder=4)
        ax.text(k + 0.17, r, f"{r:.2f}\u00d7", ha="left", va="center", fontsize=17, fontweight="bold", color=c)
        ax.text(k, 0.43, "p < 0.001" if p_ < 1e-3 else f"p = {p_:.3f}", ha="center", fontsize=11.5, color=MUTED)
    ax.axhline(1, color=INK, lw=1.4, ls=":")
    ax.text(1.46, 1.06, "no difference", ha="right", fontsize=11.5, color=MUTED)
    ax.set_yscale("log")
    ax.set_yticks([0.5, 1, 2, 4], ["0.5\u00d7", "1\u00d7", "2\u00d7", "4\u00d7"])
    ax.set_ylim(0.38, 7.5)
    ax.set_xticks([0, 1], [p[0] for p in panels], fontsize=12)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Responsive elements, relative to Stable\nin the same assay", fontsize=12)
    ax.set_title("But they disagree on which elements depend more", fontsize=14, loc="left", pad=8)
    ax.text(0.5, -0.30, "below 1\u00d7: responsive lose less   \u00b7   above 1\u00d7: responsive lose more",
            transform=ax.transAxes, ha="center", fontsize=11.5, color=MUTED)
    fig.tight_layout(w_pad=4.0)
    save(fig, "s6_mpra")


# ================================================================ 7. spacing, with the test on it
def spacing():
    """Three questions kept apart, because the old single panel ran them together.
    Left: how OFTEN a partner sits at each distance (occurrence). Middle: how much it MATTERS at that distance
    (importance, against footprint-free tiles in the same distance bin). Right: whether the classes differ in either."""
    CL = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_peak.csv").set_index("family")
    FINE = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_profile.csv")
    IMP = pd.read_csv(ROOT / "Analysis/distance_rules/tables/g_importance_by_distance.csv")
    FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
    FP = FP[~FP.own_ap1].copy()
    FP["mid"] = (FP.rel_start + FP.rel_end) / 2

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 5.9), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    SER = (("MEF2", "MEF2", ORANGE, 2.8), ("KLF / SP", "KLF/SP", GREY, 2.4))

    ax = axes[0]
    for fam, lab, c, lw in SER:
        w = FINE[FINE.family == fam]
        ax.plot(w.distance, w.share.rolling(3, center=True, min_periods=1).mean() * 100, color=c, lw=lw, label=lab)
    ax.axvspan(0, 5.5, color="#EDEDED", lw=0)
    ax.axvspan(11, 25, color=ORANGE, alpha=0.10, lw=0)
    ax.text(18, 0.12, "abutting", ha="center", va="bottom", fontsize=11.5, color=ORANGE, fontweight="bold")
    ax.set_xlabel("Distance from the AP-1 motif (bp, centre to centre)", fontsize=12)
    ax.set_ylabel("Share of that family's footprints (%)", fontsize=12)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 2.1)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    ax.set_title("How often a partner sits there", fontsize=13.5, loc="left", pad=30)
    ax.text(0.0, 1.015, f"every footprint of that family, important or not.\nMEF2 {CL.loc['MEF2', 'density_ratio_near_over_far']:.2f}\u00d7 enriched at 11\u201325 bp vs {CL.loc['KLF / SP', 'density_ratio_near_over_far']:.2f}\u00d7 for KLF/SP (FDR {CL.loc['MEF2', 'fdr']:.3f})",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5, color=MUTED, linespacing=1.45)

    ax = axes[1]
    for fam, lab, c, lw in SER:
        w = IMP[IMP.family == fam]
        ax.plot(w.mid, w.importance, color=c, lw=lw, marker="o", ms=7, label=lab)
    ax.axhline(0, color=MUTED, lw=1.0)
    ax.set_xlabel("Distance from the AP-1 motif (bp)", fontsize=12)
    ax.set_ylabel("Extra loss when that tile is replaced (log2FC)", fontsize=11.5)
    ax.set_title("How much it matters there", fontsize=13.5, loc="left", pad=30)
    ax.text(0.0, 1.015, "a MEF2 neighbour is worth most within ~50 bp and still\nworth \u22120.30 at 250 bp; KLF/SP is worth nothing anywhere",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5, color=MUTED, linespacing=1.45)
    ax.set_xlim(-10, 270)

    ax = axes[2]
    CLS = [("Age Up", "age_class", "Age Up", ORANGE), ("Age Down", "age_class", "Age Down", BLUE),
           ("Exercise Up", "ex_class", "Exercise Up", GREEN_D), ("Exercise Down", "ex_class", "Exercise Down", PURPLE_D)]
    near = set(FP[(FP.family == "MEF2") & (FP.mid.abs() <= 25)].site_id)
    stable_ids = S.index[S.stable_dynamic == "Stable"]
    ps, labs = [], []
    for k, (nm, col, val, c) in enumerate(CLS):
        f = S.copy()
        f["grp"] = np.where(f.index.isin(stable_ids), "Stable", np.where(f[col] == val, nm, ""))
        P = pairs(f[f.grp != ""], "grp", nm, "Stable")
        for ids, cc, off in ((P.case, c, -0.16), (P.control, STABLE_GREY, 0.16)):
            v = np.array([1.0 if i_ in near else 0.0 for i_ in ids])
            lo, hi = np.percentile([rng.choice(v, len(v)).mean() for _ in range(700)], [2.5, 97.5])
            ax.errorbar(k + off, v.mean() * 100, yerr=[[(v.mean() - lo) * 100], [(hi - v.mean()) * 100]],
                        fmt="o", ms=10, color=cc, lw=2.2, capsize=0)
        a = sum(i_ in near for i_ in P.case)
        b = sum(i_ in near for i_ in P.control)
        ps.append(fisher_exact([[a, len(P) - a], [b, len(P) - b]])[1])
        labs.append(nm.replace("Exercise ", "Ex "))
    q = np.argsort(ps)
    fdr = np.empty(len(ps))
    fdr[q] = np.minimum.accumulate((np.array(ps)[q] * len(ps) / np.arange(1, len(ps) + 1))[::-1])[::-1]
    for k, f_ in enumerate(fdr):
        ax.text(k, ax.get_ylim()[1] * 0.985, f"FDR {f_:.2f}" if f_ < 0.1 else "ns",
                ha="center", va="top", fontsize=10.5, color=INK if f_ < 0.1 else MUTED,
                fontweight="bold" if f_ < 0.1 else "normal")
    ax.set_xticks(range(len(CLS)), labs, fontsize=11.5, rotation=18, ha="right")
    ax.set_ylabel("Sites with a MEF2 footprint\nwithin 25 bp of the motif (%)", fontsize=11.5)
    ax.set_title("And does it differ by class? Yes", fontsize=13.5, loc="left", pad=30)
    ax.text(0.0, 1.015, "coloured = the class, grey = its own context-matched\nStable controls. Age-gained sites are short of a close MEF2.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5, color=MUTED, linespacing=1.45)
    ax.set_xlim(-0.6, len(CLS) - 0.4)
    ax.set_ylim(0, 21)
    fig.subplots_adjust(left=0.062, right=0.99, top=0.775, bottom=0.20, wspace=0.42)
    save(fig, "s7_spacing")


# ================================================================ 8. partners, ranked AND split by class
def partners():
    """The ranking is validation; the class split is the question. Dynamic is broken into its four directional classes,
    because a gain and a loss could cancel in the pooled comparison, and shown both raw and context-matched."""
    FI = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    FI = FI[FI.family != "Other"]
    cl = S[["stable_dynamic", "age_class", "ex_class"]].copy()
    CLASSES = [("Age Up", "age_class", "Age Up", ORANGE), ("Age Down", "age_class", "Age Down", BLUE),
               ("Exercise Up", "ex_class", "Exercise Up", GREEN_D), ("Exercise Down", "ex_class", "Exercise Down", PURPLE_D)]
    STCOL = STABLE_GREY
    stable_ids = set(cl.index[cl.stable_dynamic == "Stable"])
    FAM = ["MEF2", "NFI (Nfia/b/c/x)", "GATA", "TEAD", "AP-1 (Fos/Jun/Atf/Batf)", "KLF / SP"]

    # one matched control set per class, so a class is never compared against sites it could not be matched to
    matched = {}
    for nm, col, val, _c in CLASSES:
        f = S.copy()
        f["grp"] = np.where(f.index.isin(stable_ids), "Stable", np.where(f[col] == val, nm, ""))
        P = pairs(f[f.grp != ""], "grp", nm, "Stable")
        matched[nm] = (set(P.case), set(P.control), len(P))
    pooled_ctrl = set().union(*[m[1] for m in matched.values()])

    fig, axes = plt.subplots(1, 2, figsize=(15.6, 5.7), sharey=True)
    ys = np.arange(len(FAM))[::-1]
    lim = 0.0
    tests, dropped = [], 0
    for ax, mode in zip(axes, ("all", "matched")):
        for y, fam in zip(ys, FAM):
            g = FI[FI.family == fam].set_index("site_id").vs_free
            series = [("Stable", STCOL, g.reindex(sorted(stable_ids if mode == "all" else pooled_ctrl)).dropna().values)]
            for nm, col, val, c in CLASSES:
                ids = set(cl.index[cl[col] == val]) if mode == "all" else matched[nm][0]
                series.append((nm, c, g.reindex(sorted(ids)).dropna().values))
            for k, (nm, c, v) in enumerate(series):
                if len(v) < 18:
                    dropped += 1
                    continue
                lo, hi = boot(v, 400)
                lim = min(lim, lo)
                ax.errorbar(np.median(v), y + (2 - k) * 0.155, xerr=[[np.median(v) - lo], [hi - np.median(v)]],
                            fmt="o", ms=7.5, color=c, lw=1.9, capsize=0, zorder=4)
        ax.axvline(0, color=MUTED, lw=1.0, ls=":")
        ax.set_ylim(-0.62, len(FAM) - 0.38)
    axes[0].set_yticks(ys, [f.split(" (")[0].split(" /")[0] for f in FAM], fontsize=13.5)
    axes[0].set_title("All sites", fontsize=14, loc="left", pad=10)
    axes[1].set_title("Context-matched to Stable, class by class", fontsize=14, loc="left", pad=10)
    for ax in axes:
        ax.set_xlim(lim * 1.05, 0.12)
    h = [axes[0].plot([], [], "o", color=c, ms=9, label=l)[0]
         for l, c in [("Stable", STCOL)] + [(nm, c) for nm, _, _, c in CLASSES]]
    fig.subplots_adjust(left=0.135, right=0.99, top=0.92, bottom=0.23, wspace=0.07)
    fig.text(0.56, 0.135, "How much accessibility is lost when a 10 bp tile over this family's footprint is replaced (log2FC)",
             ha="center", fontsize=12.5)
    fig.legend(handles=h, frameon=False, fontsize=12, loc="lower center", ncol=5, bbox_to_anchor=(0.56, 0.045))
    for mode in ("all", "matched"):
        for fam in sorted(FI.family.unique()):
            g = FI[FI.family == fam].set_index("site_id").vs_free
            base = g.reindex(sorted(stable_ids if mode == "all" else pooled_ctrl)).dropna().values
            if len(base) < 15:
                continue
            for nm, col, val, _c in CLASSES:
                ids = set(cl.index[cl[col] == val]) if mode == "all" else matched[nm][0]
                v = g.reindex(sorted(ids)).dropna().values
                if len(v) < 15:
                    continue
                tests.append({"mode": mode, "family": fam, "cls": nm, "p": mannwhitneyu(base, v).pvalue})
    TT = pd.DataFrame(tests).sort_values("p")
    TT["fdr"] = (TT.p * len(TT) / np.arange(1, len(TT) + 1)).clip(upper=1)
    best = {m: TT[TT["mode"] == m].iloc[0] for m in ("all", "matched")}
    bm = best["matched"]
    fig.text(0.56, 0.008, "six families shown; tests run over all "
                          + f"{int((TT['mode'] == 'matched').sum())} class-by-family comparisons and none survives correction "
                          + f"(smallest matched FDR {bm.fdr:.2f}, {bm.family.split(' (')[0].split(' /')[0]} in {bm.cls}, p = {bm.p:.3f})"
                          + "  \u00b7  matched pairs per class: " + ", ".join(f"{nm} {matched[nm][2]}" for nm, _, _, _ in CLASSES),
             ha="center", fontsize=10.5, color=MUTED)
    save(fig, "s8_partners")


print("building slide figures")
element_panel(); bands(); matched_null(); dependence(); titration(); age_gained(); mpra(); spacing(); partners()
print("done ->", OUT)
