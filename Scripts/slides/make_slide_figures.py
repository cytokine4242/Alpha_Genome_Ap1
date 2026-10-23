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


def class_dot_panel(ax, v, ylab, mode, minn=15):
    """One measure, one mode. mode="all": every class against all Stable sites, nothing controlled -- if the
    accessibility difference is real biology it lives here. mode="matched": every class against its own matched
    Stable controls, paired. Both are shown for every result; neither is the privileged one."""
    ps = {}
    base_all = v.reindex(sorted(_STABLE)).dropna().values
    if mode == "all":
        lo, hi = boot(base_all, 400)
        ax.axhspan(lo, hi, color=CTRL_GREY, alpha=0.30, lw=0, zorder=1)
        ax.axhline(np.median(base_all), color=STABLE_GREY, lw=1.4, zorder=2)
    for k, (nm, col_, val_, col) in enumerate(GROUPS):
        if mode == "all":
            a = v.reindex(sorted(S.index[S[col_] == val_])).dropna().values
            b = base_all
            if len(a) < minn:
                continue
            ps[nm] = mannwhitneyu(a, b).pvalue
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
            ax.errorbar(k + off, np.median(vals), yerr=[[np.median(vals) - lo], [hi - np.median(vals)]],
                        fmt="o", ms=9, color=c_, lw=2.2, capsize=0, zorder=4)
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_, hi_ + (hi_ - lo_) * 0.16)
    ax.set_xticks(range(len(GROUPS)), [n.replace("Exercise ", "Ex ").replace(" (all)", "") for n, _, _, _ in GROUPS],
                  fontsize=10.5, rotation=20, ha="right")
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
        ax.text(k, 0.985, ("***" if q < 1e-3 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"),
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=12,
                color=INK if q < 0.05 else MUTED, fontweight="bold" if q < 0.05 else "normal")



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
    """The four edits that carry the argument, every directional class, shown both ways: against all Stable sites
    and against each class's own matched controls. The scale-specific residual is the point -- the motif and the
    whole peak close up under matching, the 170 bp element does not."""
    KEY = [("replace motif, read over motif", "Replace the 11 bp motif", "the motif alone"),
           ("replace 170 bp, read over 170 bp", "Replace the 170 bp element", "the element around it"),
           ("replace peak, read over peak", "Replace the whole peak", "everything, read over the peak"),
           ("170 bp shuffle \u00b1300 bp away, read over motif", "Edit 300 bp away", "the control: nothing at stake")]
    d = S.join(FD[[k for k, _, _ in KEY]])
    fig, axes = plt.subplots(2, 4, figsize=(16.0, 7.4), sharex=True)
    allps, flat = [], []
    for c_, (col, lab, defn) in enumerate(KEY):
        for r_, mode in enumerate(("all", "matched")):
            ax = axes[r_, c_]
            allps.append(class_dot_panel(ax, d[col], "Effect on predicted\naccessibility (log2FC)" if c_ == 0 else "", mode))
            flat.append(ax)
            if r_ == 0:
                ax.set_title(lab, fontsize=13, loc="left", pad=26)
                ax.text(0.0, 1.015, defn, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5, color=MUTED)
    mark_fdr(flat, allps)
    for r_, txt in enumerate(("all sites, nothing controlled\n(grey band = all Stable sites)",
                              "each class against its own\ncontext-matched Stable controls")):
        axes[r_, 0].text(-0.42, 0.5, txt, transform=axes[r_, 0].transAxes, rotation=90, ha="center", va="center",
                         fontsize=11, color=INK, linespacing=1.5)
    h = [axes[0, 0].plot([], [], "o", color=col, ms=9, label=nm)[0] for nm, _, _, col in GROUPS]
    h.append(axes[0, 0].plot([], [], "o", color=CTRL_GREY, ms=9, label="Stable (all, or matched controls)")[0])
    fig.subplots_adjust(left=0.115, right=0.99, top=0.875, bottom=0.265, wspace=0.34, hspace=0.26)
    fig.legend(handles=h, frameon=False, fontsize=11.5, loc="lower center", ncol=6, bbox_to_anchor=(0.5, 0.008))
    fig.text(0.55, 0.145, TEST_KEY.format(n=sum(len(d) for d in allps)), ha="center", va="top",
             fontsize=10, color=MUTED, linespacing=1.6)
    save(fig, "s2_matched_null")


# ================================================================ 3. dependence shape, raw and matched
def _profile_schematic(ax):
    """What the four measures are read off: one site's sliding-window profile, with the near window marked."""
    ax.axis("off")
    ax.set_xlim(-310, 310)
    ax.set_ylim(-1.95, 0.80)
    rng2 = np.random.default_rng(7)
    c = np.arange(-295, 300, 10)
    prof = -(0.95 * np.exp(-(c ** 2) / (2 * 55 ** 2)) + 0.30 * np.exp(-((c - 150) ** 2) / (2 * 40 ** 2))
             + 0.22 * np.exp(-((c + 190) ** 2) / (2 * 45 ** 2))) - rng2.uniform(0, 0.05, len(c))
    ax.axvspan(-85, 85, color=ORANGE, alpha=0.10, lw=0, zorder=1)
    ax.bar(c[np.abs(c) <= 85], prof[np.abs(c) <= 85], width=9, color=ORANGE, zorder=3)
    ax.bar(c[np.abs(c) > 85], prof[np.abs(c) > 85], width=9, color=CTRL_GREY, zorder=3)
    ax.axhline(0, color=RAIL, lw=1.2, zorder=4)
    ax.text(0, 0.60, "one site's sliding-window profile", ha="center", fontsize=11.5, color=INK, fontweight="bold")
    ax.text(0, 0.30, "each bar = one 10 bp window replaced", ha="center", fontsize=10, color=MUTED)
    ax.annotate("", xy=(-85, -1.32), xytext=(85, -1.32), arrowprops=dict(arrowstyle="|-|,widthA=0.3,widthB=0.3", color=ORANGE, lw=1.2))
    ax.text(0, -1.40, "near: \u00b185 bp", ha="center", va="top", fontsize=10.5, color=ORANGE, fontweight="bold")
    ax.text(-200, 0.40, "far", ha="center", fontsize=10.5, color=MUTED)
    ax.text(200, 0.40, "far", ha="center", fontsize=10.5, color=MUTED)
    ax.annotate("", xy=(0, -1.72), xytext=(-185, -1.72), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.0))
    ax.text(-92, -1.68, "spread", ha="center", va="bottom", fontsize=10, color=MUTED)
    ax.text(0, 0.06, "AP-1 motif", ha="center", va="bottom", fontsize=9.5, color=ORANGE)


def dependence():
    """How much of a site's sequence dependence sits near its motif and how much sits out in the surrounding context,
    in absolute log2 units rather than as a share, so 'more distributed' can be read directly. Every class, both ways."""
    SCF2 = SCF.assign(near_abs_170=SCF.total_dependence_170 * SCF.share_within_85_170,
                      far_abs_170=SCF.total_dependence_170 * (1 - SCF.share_within_85_170))
    M = [("near_abs_170", "Dependence near the motif",
          "loss summed over windows within \u00b185 bp", "Summed loss within\n\u00b185 bp (log2)"),
         ("far_abs_170", "Dependence out in the context",
          "loss summed over windows beyond that", "Summed loss beyond\n\u00b185 bp (log2)"),
         ("spread_170", "Spread along the element",
          "SD of window position, weighted by loss", "Spread of the loss (bp)"),
         ("ism_max", "Depth of the single best base",
          "largest effect of any one base change, \u00b130 bp", "Largest single-base\neffect (|log2FC|)")]
    d = S.join(SCF2[[m[0] for m in M]], how="inner")
    fig = plt.figure(figsize=(16.6, 7.4))
    gs = fig.add_gridspec(2, 6, width_ratios=[1.2, 0.22, 1, 1, 1, 1], wspace=0.42, hspace=0.26,
                          left=0.012, right=0.99, top=0.855, bottom=0.265)
    _profile_schematic(fig.add_subplot(gs[:, 0]))
    allps, flat = [], []
    for c_, (col, lab, defn, ylab) in enumerate(M):
        for r_, mode in enumerate(("all", "matched")):
            ax = fig.add_subplot(gs[r_, c_ + 2])
            allps.append(class_dot_panel(ax, d[col], ylab, mode))
            flat.append(ax)
            if r_ == 0:
                ax.set_title(lab, fontsize=12.5, loc="left", pad=26)
                ax.text(0.0, 1.015, defn, transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color=MUTED)
                ax.tick_params(labelbottom=False)
    mark_fdr(flat, allps)
    for y_, txt in ((0.665, "all sites,\nnothing controlled"), (0.40, "each class vs its own\nmatched controls")):
        fig.text(0.205, y_, txt, rotation=90, ha="center", va="center", fontsize=10.5, color=INK, linespacing=1.5)
    h = [flat[0].plot([], [], "o", color=col, ms=9, label=nm)[0] for nm, _, _, col in GROUPS]
    h.append(flat[0].plot([], [], "o", color=CTRL_GREY, ms=9, label="Stable (all, or matched controls)")[0])
    fig.legend(handles=h, frameon=False, fontsize=11.5, loc="lower center", ncol=6, bbox_to_anchor=(0.55, 0.008))
    fig.text(0.55, 0.145, TEST_KEY.format(n=sum(len(d_) for d_ in allps)), ha="center", va="top",
             fontsize=10, color=MUTED, linespacing=1.6)
    save(fig, "s3_dependence")


# ================================================================ 4. titration, everything against native
def titration():
    """Mechanism split by class throughout. Left column: the dose curve and the position curve drawn once per class.
    Right: four titration quantities per site, every class, both against all Stable sites and against matched controls."""
    tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/titrate*_cm_chunk*.csv")))], ignore_index=True)
    tt = tt[tt.trackset == "CM"]
    nat = tt[tt.kind == "native"].set_index("site_id").atac_170
    pc = 0.05 * nat.median()
    e = tt[tt.kind != "native"].copy()
    e["eff"] = np.log2((e.atac_170 + pc) / (e.site_id.map(nat) + pc))
    TF = pd.read_csv(ROOT / "Analysis/titration/tables/titration_features_per_site.csv").set_index("site_id")
    M = [("mm1_170", "One damaged core base", "TGACTCA with a single base changed", "Effect (log2FC)"),
         ("mm7_170", "All seven core bases", "the motif destroyed outright", "Effect (log2FC)"),
         ("add_copy_net_170", "A second copy added", "net of a shuffled-sequence control", "Effect (log2FC)"),
         ("relocate_net_170", "The motif moved elsewhere", "net of the same control", "Effect (log2FC)")]
    d = S.join(TF[[m[0] for m in M]], how="left")
    fig = plt.figure(figsize=(16.6, 7.4))
    gs = fig.add_gridspec(2, 6, width_ratios=[1.35, 0.26, 1, 1, 1, 1], wspace=0.42, hspace=0.46,
                          left=0.045, right=0.99, top=0.855, bottom=0.265)
    CURVE = [("Stable", None, None, CTRL_GREY)] + [(n, c, v, col) for n, c, v, col in GROUPS if n != "Dynamic (all)"]

    def curve_members(nm, col_, val_):
        return set(S.index[S.stable_dynamic == "Stable"]) if nm == "Stable" else set(S.index[S[col_] == val_])

    ax = fig.add_subplot(gs[0, 0])
    for nm, col_, val_, c in CURVE:
        ids = curve_members(nm, col_, val_)
        w = e[(e.kind == "core_mm") & (e.site_id.isin(ids))].groupby("pos").eff.median()
        ax.plot(w.index, w.values, color=c, lw=2.2, marker="o", ms=5, label=nm.replace("Exercise ", "Ex "))
    ax.set_xlabel("Core bases changed (of 7)", fontsize=11)
    ax.set_ylabel("Loss of predicted\naccessibility (log2FC)", fontsize=11)
    ax.set_title("The dose curve, class by class", fontsize=12.5, loc="left", pad=6)
    ax.text(0.99, 0.04, "colours as in the legend below", transform=ax.transAxes, ha="right", fontsize=9.5, color=MUTED)

    ax = fig.add_subplot(gs[1, 0])
    for nm, col_, val_, c in CURVE:
        ids = curve_members(nm, col_, val_)
        w = e[(e.kind == "dup") & (e.site_id.isin(ids))].groupby("pos").eff.median()
        ax.plot(w.index, w.values, color=c, lw=2.2, marker="o", ms=4)
    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xlabel("Position of the added copy (bp)", fontsize=11)
    ax.set_ylabel("Change vs native (log2FC)", fontsize=11)
    ax.set_title("Where a second copy helps", fontsize=12.5, loc="left", pad=6)

    allps, flat = [], []
    for c_, (col, lab, defn, ylab) in enumerate(M):
        for r_, mode in enumerate(("all", "matched")):
            ax = fig.add_subplot(gs[r_, c_ + 2])
            allps.append(class_dot_panel(ax, d[col], ylab if c_ == 0 else "", mode))
            flat.append(ax)
            if r_ == 0:
                ax.set_title(lab, fontsize=12.5, loc="left", pad=26)
                ax.text(0.0, 1.015, defn, transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color=MUTED)
                ax.tick_params(labelbottom=False)
    mark_fdr(flat, allps)
    for y_, txt in ((0.665, "all sites,\nnothing controlled"), (0.40, "each class vs its own\nmatched controls")):
        fig.text(0.238, y_, txt, rotation=90, ha="center", va="center", fontsize=10.5, color=INK, linespacing=1.5)
    h = [flat[0].plot([], [], "o", color=col, ms=9, label=nm)[0] for nm, _, _, col in GROUPS]
    h.append(flat[0].plot([], [], "o", color=CTRL_GREY, ms=9, label="Stable (all, or matched controls)")[0])
    fig.legend(handles=h, frameon=False, fontsize=11.5, loc="lower center", ncol=6, bbox_to_anchor=(0.55, 0.008))
    fig.text(0.55, 0.145, TEST_KEY.format(n=sum(len(d_) for d_ in allps)), ha="center", va="top",
             fontsize=10, color=MUTED, linespacing=1.6)
    save(fig, "s4_titration")


# ================================================================ 5. the age-gained result
def age_gained():
    TI = pd.read_csv(ROOT / "Analysis/titration/tables/titration_by_class.csv")
    ti = lambda f: TI[(TI.contrast == "Age Up vs Unchanged (Δscore)") & (TI.feature.str.startswith(f))].iloc[0]
    reg = pd.read_csv(ROOT / "Analysis/hypotheses/tables/hypothesis_register.csv").set_index("id")
    items = [("net value of a second", "Gains more from\nan added motif copy", "+0.38", "+0.32"),
             ("does a second copy help", "Benefits more when\nthat copy is near", "+0.21", "+0.14"),
             ("effect of changing 1 core", "Loses more from one\ndamaged core base", "−0.60", "−0.51")]
    fig, ax = plt.subplots(figsize=(10.4, 4.4))
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

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), gridspec_kw={"width_ratios": [1.15, 1]})
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

    fig.subplots_adjust(left=0.085, right=0.985, top=0.80, bottom=0.165, wspace=0.34)
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


# ================================================================ 9. the MEF2 deficit, by class, two ways
def mef2_by_class():
    """The close-MEF2 result given its own figure. Two definitions of 'has one' -- a footprint at all, and a footprint
    the model actually leans on -- each shown against all Stable sites and against each class's own matched controls."""
    FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
    FP = FP[~FP.own_ap1].copy()
    FP["mid"] = (FP.rel_start + FP.rel_end) / 2
    FI = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    imp = FI[FI.family == "MEF2"].set_index("site_id").vs_free
    near = set(FP[(FP.family == "MEF2") & (FP.mid.abs() <= 25)].site_id)
    leaned = {i for i in near if imp.get(i, 0) < -0.2}
    M = [(pd.Series({i: float(i in near) for i in S.index}), "A MEF2 footprint within 25 bp",
          "counting footprints, whether or not the model uses them",
          "Sites with a MEF2 footprint\nwithin 25 bp (%)"),
         (pd.Series({i: float(i in leaned) for i in S.index}), "...one the model actually leans on",
          "the same, restricted to sites where replacing that tile\ncosts more than 0.2 log2 against footprint-free tiles",
          "Sites with a close MEF2\nthe model leans on (%)")]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6), sharex=True)
    allps, flat = [], []
    for c_, (flag, lab, defn, ylab) in enumerate(M):
        for r_, mode in enumerate(("all", "matched")):
            ax = axes[r_, c_]
            allps.append(class_prop_panel(ax, flag, ylab, mode))
            flat.append(ax)
            if r_ == 0:
                ax.set_title(lab, fontsize=13.5, loc="left", pad=40)
                ax.text(0.0, 1.025, defn, transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=10.5, color=MUTED, linespacing=1.45)
    mark_fdr(flat, allps)
    for r_, txt in enumerate(("all sites, nothing controlled\n(grey band = all Stable sites)",
                              "each class against its own\ncontext-matched Stable controls")):
        axes[r_, 0].text(-0.33, 0.5, txt, transform=axes[r_, 0].transAxes, rotation=90, ha="center", va="center",
                         fontsize=11, color=INK, linespacing=1.5)
    h = [axes[0, 0].plot([], [], "o", color=col, ms=9, label=nm)[0] for nm, _, _, col in GROUPS]
    h.append(axes[0, 0].plot([], [], "o", color=CTRL_GREY, ms=9, label="Stable (all, or matched controls)")[0])
    fig.subplots_adjust(left=0.155, right=0.985, top=0.835, bottom=0.235, wspace=0.36, hspace=0.26)
    fig.legend(handles=h, frameon=False, fontsize=11, loc="lower center", ncol=3, bbox_to_anchor=(0.56, 0.055))
    fig.text(0.56, 0.045, "grey line and band: all Stable sites with a 95% bootstrap CI.   grey dots: that class's own matched controls.\n"
                          "Fisher on all sites, exact McNemar on the matched pairs, BH across all "
                          f"{sum(len(d) for d in allps)} tests.    *** p < 0.001    ** p < 0.01    * p < 0.05    ns \u2265 0.05",
             ha="center", va="top", fontsize=10, color=MUTED, linespacing=1.6)
    save(fig, "s9_mef2_by_class")


print("building slide figures")
element_panel(); bands(); matched_null(); dependence(); titration(); age_gained(); mpra(); spacing(); partners(); mef2_by_class()
print("done ->", OUT)
