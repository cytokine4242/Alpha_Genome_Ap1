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
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Analysis" / "slide_figures"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(50)
# Same visual language as the analysis figures in the HTML reports: an off-white panel, a recessive grid, muted axis
# furniture, ink only where the data is. Type is larger here because these are read from a projector, not a page.
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, GREY, GREEN = "#1F6FB2", "#C2551C", "#9a9993", "#1baf7a"
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


def schematic():
    """One accent colour for 'this is the part we edited', grey for everything else, and real axis furniture.
    Every row below shares the same rail and the same chip geometry, so the three experiments read as variations of one idea."""
    fig = plt.figure(figsize=(13.6, 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.5], hspace=0.30)
    SEQ, RAIL = "#d8d7d2", "#3f3e3b"

    # ---- A. the element, drawn as a browser track
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(-345, 345)
    ax.set_ylim(-0.62, 1.22)
    ax.axis("off")
    x = np.linspace(-345, 345, 600)
    peak = 0.70 * np.exp(-(x ** 2) / (2 * 155 ** 2))
    ax.fill_between(x, 0.02, peak, color=BLUE, alpha=0.13, lw=0, zorder=1)
    ax.plot(x, peak, color=BLUE, lw=2.0, zorder=3, solid_capstyle="round")
    eyebrow(ax, 0.0, 0.90, "the element")
    ax.text(-345, 0.86, "a cardiomyocyte ATAC peak with an AP-1 footprint in it", fontsize=12.5, color=INK, va="bottom")
    ax.annotate("", xy=(-85, 0.20), xytext=(85, 0.20), arrowprops=dict(arrowstyle="|-|,widthA=0.4,widthB=0.4", color=MUTED, lw=1.2))
    ax.text(0, 0.235, "170 bp — the window every effect is read over", ha="center", fontsize=11, color=MUTED)
    ax.plot([-345, 345], [0.02, 0.02], color=RAIL, lw=1.4, zorder=4)
    for xx, lab, c, h in ((-215, "MEF2", GREY, 15), (0, "AP-1", ORANGE, 11), (120, "GATA", GREY, 14)):
        ax.add_patch(Rectangle((xx - h / 2 - 4, -0.075), h + 8, 0.19, color=c, ec="none", zorder=5))
        ax.text(xx, -0.135, lab, ha="center", va="top", fontsize=11.5, color=c if c == ORANGE else MUTED,
                fontweight="bold" if c == ORANGE else "normal")
    for xx in (-300, -150, 0, 150, 300):
        ax.plot([xx, xx], [0.02, -0.02], color=MUTED, lw=1, zorder=4)
        ax.text(xx, -0.30, f"{xx:+d} bp" if xx else "0", ha="center", va="top", fontsize=10.5, color=MUTED)
    ax.text(0, -0.50, "AlphaGenome reads 131 kb of sequence and predicts accessibility; we change the middle of it and read what the prediction loses",
            ha="center", va="top", fontsize=12, color=MUTED)

    # ---- B. the three edits, one shared rail
    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.42, 3.50)
    ax.axis("off")
    eyebrow(ax, 0.0, 0.955, "what we changed")
    X0, W = 4.25, 5.55
    rows = [(2.50, "Sliding window", "every 10 bp from −300 to +300 bp, replaced by random bases", "tile"),
            (1.35, "Single bases", "all three alternatives at every base, ±30 bp", "ism"),
            (0.25, "Motif titration", "1–7 core bases changed; the motif moved, duplicated or removed", "titr")]
    for i, (y, head, sub, kind) in enumerate(rows):
        ax.text(0.0, y + 0.50, f"{i + 1}", fontsize=13, color=SURF, ha="center", va="center", fontweight="bold",
                bbox=dict(boxstyle="circle,pad=0.34", fc=RAIL, ec="none"))
        ax.text(0.30, y + 0.63, head, fontsize=13.5, fontweight="bold", va="top")
        ax.text(0.30, y + 0.29, sub, fontsize=11.5, color=MUTED, va="top")
        if kind != "ism":                                     # the base row is drawn as sequence, so it needs no rail under it
            ax.plot([X0, X0 + W], [y + 0.40, y + 0.40], color=RAIL, lw=1.2, zorder=2)
        if kind == "tile":
            for k in range(12):
                xx = X0 + k * W / 12
                ax.add_patch(Rectangle((xx + 0.012, y + 0.28), W / 12 - 0.024, 0.24,
                                       color=ORANGE if k == 5 else SEQ, ec="none", zorder=3))
            ax.annotate("", xy=(X0 + W * 0.62, y + 0.74), xytext=(X0 + W * 0.44, y + 0.74),
                        arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.6, mutation_scale=14))
            ax.text(X0 + W * 0.42, y + 0.74, "one window at a time  ", ha="right", va="center", fontsize=11, color=ORANGE)
        elif kind == "ism":
            letters = "C A G T G A C T C A G"
            for k, ch in enumerate(letters.split(" ")):
                xx = X0 + (k + 0.5) * W / 11
                ax.text(xx, y + 0.40, ch, ha="center", va="center", fontsize=12, family="monospace",
                        color=ORANGE if k == 5 else MUTED, fontweight="bold" if k == 5 else "normal", zorder=4)
            xx = X0 + 5.5 * W / 11
            ax.add_patch(Rectangle((xx - W / 22, y + 0.26), W / 11, 0.28, fc="none", ec=ORANGE, lw=1.4, zorder=3))
            for j, alt in enumerate("ACT"):
                ax.text(xx + (j - 1) * 0.30, y + 0.03, alt, ha="center", va="center", fontsize=11,
                        family="monospace", color=ORANGE, zorder=4)
            ax.annotate("", xy=(xx, y + 0.14), xytext=(xx, y + 0.24), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.2, mutation_scale=11))
            ax.text(X0 + W + 0.06, y + 0.40, "× every position", ha="left", va="center", fontsize=11, color=MUTED)
        else:
            mw = W * 0.105
            ax.add_patch(Rectangle((X0 + W * 0.30, y + 0.28), mw, 0.24, color=ORANGE, ec="none", zorder=3))
            for frac in (0.06, 0.60, 0.82):
                ax.add_patch(Rectangle((X0 + W * frac, y + 0.28), mw, 0.24, fc="none", ec=ORANGE, lw=1.3, ls=(0, (2.5, 1.8)), zorder=3))
            ax.annotate("", xy=(X0 + W * 0.60, y + 0.24), xytext=(X0 + W * 0.365, y + 0.24),
                        arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.5, mutation_scale=13, connectionstyle="arc3,rad=0.45"))
            ax.text(X0 + W * 0.49, y + 0.62, "moved, or a second copy added", ha="center", va="bottom", fontsize=11, color=ORANGE)
    ax.text(0.30, -0.26, "3,678 sites · 364 edits each · 3 replicates · every effect is log2(edited / native)",
            fontsize=11.5, color=MUTED, va="top")
    fig.subplots_adjust(left=0.012, right=0.988, top=0.97, bottom=0.03)
    save(fig, "s1_method_schematic")


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
    M = [("share_within_85_170", "Share of the loss\nwithin ±85 bp"), ("spread_170", "Spread of the loss\n(bp)"),
         ("ism_max", "Largest single-base\neffect (log2FC)")]
    d = S.join(SCF[[m for m, _ in M]], how="inner")
    P = pairs(d, "stable_dynamic", "Dynamic", "Stable")
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.6))
    for ax, (col, lab) in zip(axes, M):
        xs, labels = [], []
        for mode in ("raw", "matched"):
            if mode == "raw":
                vs = [d[d.stable_dynamic == g][col].dropna().values for g in ("Stable", "Dynamic")]
            else:
                vs = [d[col].reindex(P.control).dropna().values, d[col].reindex(P.case).dropna().values]
            xs.append(vs)
        for j, (vs, mode) in enumerate(zip(xs, ("all sites", "matched"))):
            for i, (v, c) in enumerate(zip(vs, (BLUE, ORANGE))):
                lo, hi = boot(v)
                ax.errorbar(j + (i - 0.5) * 0.26, np.median(v), yerr=[[np.median(v) - lo], [hi - np.median(v)]],
                            fmt="o", ms=11, color=c, lw=2.4, capsize=0)
        a, b = d[col].reindex(P.case).values, d[col].reindex(P.control).values
        ok = ~(np.isnan(a) | np.isnan(b))
        p = wilcoxon(a[ok] - b[ok]).pvalue
        ax.set_xticks([0, 1], ["all sites", "matched"], fontsize=12)
        ax.set_xlim(-0.45, 1.45)
        ax.set_title(lab, fontsize=13.5, loc="left", pad=8)
        ax.text(0.5, -0.22, f"matched p = {p:.1g}" if p >= 1e-4 else "matched p < 0.0001",
                transform=ax.transAxes, ha="center", fontsize=12, color=MUTED)
    h = [axes[0].plot([], [], "o", color=c, ms=10, label=l)[0] for c, l in ((BLUE, "Stable"), (ORANGE, "Dynamic"))]
    fig.legend(handles=h, frameon=False, fontsize=13, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    save(fig, "s3_dependence")


# ================================================================ 4. titration, everything against native
def titration():
    tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/titrate*_cm_chunk*.csv")))], ignore_index=True)
    tt = tt[tt.trackset == "CM"]
    nat = tt[tt.kind == "native"].set_index("site_id").atac_170
    pc = 0.05 * nat.median()
    e = tt[tt.kind != "native"].copy()
    e["eff"] = np.log2((e.atac_170 + pc) / (e.site_id.map(nat) + pc))
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.9), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    mm = e[e.kind == "core_mm"].groupby("pos").eff
    med = mm.median()
    lo = mm.apply(lambda v: boot(v.values)[0])
    hi = mm.apply(lambda v: boot(v.values)[1])
    ax.plot(med.index, med.values, color=ORANGE, lw=2.6, marker="o", ms=9)
    ax.fill_between(med.index, lo.values, hi.values, color=ORANGE, alpha=0.18, lw=0)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlabel("Core bases changed (of 7)", fontsize=13)
    ax.set_ylabel("Loss of predicted accessibility (log2FC)", fontsize=12.5)
    ax.set_title("One wrong base does almost all the damage", fontsize=14, loc="left", pad=8)
    ax.annotate(f"{med.iloc[0]:.2f} after one base", xy=(1, med.iloc[0]), xytext=(2.4, med.iloc[0] + 0.07),
                fontsize=12, color=ORANGE, arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.5))
    ax = axes[1]
    for kind, lab, c, ls in (("ctrlmove", "motif destroyed, nothing put back", GREY, "--"),
                             ("move", "motif destroyed, re-inserted elsewhere", ORANGE, "-"),
                             ("dup", "original kept, a second copy added", GREEN, "-")):
        w = e[e.kind == kind].groupby("pos").eff.median()
        ax.plot(w.index, w.values, color=c, lw=2.6, ls=ls, marker="o", ms=8, label=lab)
    ax.axhline(0, color=INK, lw=1.3)
    ax.text(152, 0.045, "native", fontsize=12, color=INK, ha="right")
    ax.set_xlabel("Position of the inserted copy, relative to the motif (bp)", fontsize=13)
    ax.set_ylabel("Change vs native (log2FC)", fontsize=12.5)
    ax.set_title("A motif put back elsewhere never fully replaces the original", fontsize=14, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=11.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, handlelength=2.4)
    ax.set_ylim(-0.92, 0.33)
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
    """One panel: the important partner sits just outside the motif, the unimportant one does not.
    The helical-phasing test is a report detail (MEF2 raw p = 0.047, FDR 0.47 across families) and stays in the notes."""
    CL = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_peak.csv").set_index("family")
    MD = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_modal_distance.csv").set_index("family")
    NE = pd.read_csv(ROOT / "Analysis/distance_rules/tables/b_distance_summary.csv").set_index("family")
    FINE = pd.read_csv(ROOT / "Analysis/distance_rules/tables/f_close_range_profile.csv")
    fig, ax = plt.subplots(figsize=(11.8, 5.6))
    for fam, lab, c, lw in (("MEF2", "MEF2 \u2014 the partner that matters most", ORANGE, 3.0),
                            ("KLF / SP", "KLF/SP \u2014 the commonest, and matters least", GREY, 2.4)):
        w = FINE[FINE.family == fam]
        ax.plot(w.distance, w.share.rolling(3, center=True, min_periods=1).mean() * 100, color=c, lw=lw, label=lab)
    ax.axvspan(0, 5.5, color="#EDEDED", lw=0)
    ax.text(3.2, 1.88, "overlap\n(impossible)", ha="center", va="top", fontsize=11, color=MUTED)
    ax.axvspan(11, 25, color=ORANGE, alpha=0.10, lw=0)
    ax.text(18, 1.88, "abutting", ha="center", va="top", fontsize=12.5, color=ORANGE, fontweight="bold")
    ax.set_xlabel("Distance between the AP-1 motif and its neighbour (bp, centre to centre)", fontsize=13)
    ax.set_ylabel("Share of that family's sites (%)", fontsize=13)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 2.5)
    ax.legend(frameon=False, fontsize=12.5, loc="upper right", bbox_to_anchor=(1.0, 1.02))
    r = CL.loc["MEF2"]
    fig.subplots_adjust(left=0.085, right=0.985, top=0.97, bottom=0.235)
    fig.text(0.535, 0.045, f"MEF2 {r.density_ratio_near_over_far:.2f}\u00d7 enriched at 11\u201325 bp against {CL.loc['KLF / SP', 'density_ratio_near_over_far']:.2f}\u00d7 for KLF/SP (FDR {r.fdr:.3f})"
                           f"   \u00b7   most common spacing {MD.loc['MEF2', 'modal_distance_9_40bp']:.0f} bp"
                           f"   \u00b7   median distance MEF2 {NE.loc['MEF2', 'median_distance']:.0f} bp vs KLF/SP {NE.loc['KLF / SP', 'median_distance']:.0f} bp",
             ha="center", va="bottom", fontsize=12, color=INK)
    save(fig, "s7_spacing")


# ================================================================ 8. partners, ranked AND split by class
def partners():
    """The partner ranking is validation; the class split is the question. Both read off one panel: families ordered by how
    much the model leans on them, each shown for Dynamic and for its context-matched Stable partner."""
    FI = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    FI = FI[FI.family != "Other"]
    P = pairs(S, "stable_dynamic", "Dynamic", "Stable")
    cls = pd.Series("", index=S.index)
    cls.loc[P.case] = "Dynamic"
    cls.loc[P.control] = "Stable"
    FI["cls"] = FI.site_id.map(cls).fillna("")
    w = FI[FI.cls != ""]
    keep = [f for f, g in w.groupby("family") if (g.cls == "Dynamic").sum() >= 35 and (g.cls == "Stable").sum() >= 35]
    rank = FI[FI.family.isin(keep)].groupby("family").vs_free.median().sort_values()
    PIN = ["AP-1 (Fos/Jun/Atf/Batf)", "KLF / SP"]                 # the two that make the point by mattering least
    order = [f for f in rank.index if f not in PIN][:8] + [f for f in PIN if f in rank.index]
    fig, ax = plt.subplots(figsize=(12.4, 5.8))
    ys = np.arange(len(order))[::-1]
    lim = 0.0
    for y, fam in zip(ys, order):
        g = w[w.family == fam]
        ps = []
        for c_, col, off in (("Stable", BLUE, 0.19), ("Dynamic", ORANGE, -0.19)):
            v = g[g.cls == c_].vs_free.values
            lo, hi = boot(v)
            lim = min(lim, lo)
            ax.errorbar(np.median(v), y + off, xerr=[[np.median(v) - lo], [hi - np.median(v)]],
                        fmt="o", ms=10, color=col, lw=2.4, capsize=0, zorder=4)
            ps.append(v)
        p = mannwhitneyu(*ps).pvalue
        ax.text(0.12, y, stars(p), ha="left", va="center", fontsize=13.5, color=MUTED)
    ax.axvline(0, color=MUTED, lw=1.0, ls=":")
    ax.set_yticks(ys, [f.split(" (")[0].split(" /")[0] for f in order], fontsize=14)
    ax.set_xlim(lim * 1.06, 0.26)
    ax.set_ylim(-0.62, len(order) - 0.38)
    ax.set_xlabel("How much accessibility is lost when a 10 bp tile over this family's footprint is replaced (log2FC)", fontsize=12.5)
    h = [ax.plot([], [], "o", color=c, ms=10, label=l)[0] for c, l in ((BLUE, "Stable"), (ORANGE, "Dynamic (context-matched)"))]
    fig.subplots_adjust(left=0.165, right=0.985, top=0.965, bottom=0.225)
    fig.legend(handles=h, frameon=False, fontsize=13, loc="lower center", ncol=2, bbox_to_anchor=(0.575, 0.06))
    fig.text(0.575, 0.012, f"{len(P)} matched pairs · ns = the two classes do not differ for that family (Mann–Whitney, uncorrected; nothing survives correction)",
             ha="center", fontsize=11.5, color=MUTED)
    save(fig, "s8_partners")


print("building slide figures")
schematic(); matched_null(); dependence(); titration(); age_gained(); mpra(); spacing(); partners()
print("done ->", OUT)
