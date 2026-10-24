#!/usr/bin/env python3
"""Two lab-meeting figures from the fibroblast runs (fine-tuned 1 bp ATAC head), one message each.

  fb1_autonomy      Transplanting elements into closed chromatin: Stable AP-1 elements open it on their own,
                    Dynamic ones much less, and the information that does this lives at peak scale, not in 170 bp.
                    Inputs Analysis/fb_swap/tables/{host_insertion_lift,host_insertion_class_tests,sequence_vs_location}.csv
  fb2_ap1_specific  Shuffling the site's own AP-1 motif and each partner motif inside the same window: only the
                    AP-1-losing sites are disengaged, and only from AP-1.
                    Inputs Analysis/fb_partner_shuffle/tables/{own_vs_partners_per_site,net_effect_by_family}.csv

Writes into Analysis/slide_figures/.   python3 Scripts/slides/make_fb_slide_figures.py
"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Analysis" / "slide_figures"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(54)
# same visual language as the other slide figures: off-white panel, recessive grid, ink only where the data is
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df"
STABLE_GREY, ORANGE, BLUE, RED, CTRL_GREY = "#4F4F4F", "#C2551C", "#1F6FB2", "#b3261e", "#B3B2AE"
plt.rcParams.update({"font.size": 12.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.7, "axes.axisbelow": True, "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def boot(v, n=800):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5]) if len(v) > 5 else (np.nan, np.nan)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor=SURF)
    plt.close(fig)
    print("  wrote", OUT / f"{name}.png")


# ================================================================ 1. sequence autonomy
SW = ROOT / "Analysis" / "fb_swap" / "tables"
lift = pd.read_csv(SW / "host_insertion_lift.csv")
tests = pd.read_csv(SW / "host_insertion_class_tests.csv")
sl = pd.read_csv(SW / "sequence_vs_location.csv")

ORDER = ["Stable", "Age Down", "Age Up", "Exercise Down", "Exercise Up"]
LAB = {"Stable": "Stable", "Age Down": "AP-1 lost with age", "Age Up": "AP-1 gained with age",
       "Exercise Down": "AP-1 lost with exercise", "Exercise Up": "AP-1 gained with exercise"}
p_vs_stable = {}
for r in tests[tests.width == "peak"].itertuples():
    if r.class_a == "Stable":
        p_vs_stable[r.class_b] = r.p

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 4.9), gridspec_kw={"width_ratios": [1.45, 1]})
pk = lift[lift.width == "peak"].set_index("donor_class")
y = np.arange(len(ORDER))[::-1]
for yi, cls in zip(y, ORDER):
    r = pk.loc[cls]
    col = STABLE_GREY if cls == "Stable" else ORANGE
    ax.barh(yi, r["median"], height=0.62, color=col, alpha=0.9 if cls == "Stable" else 0.78, zorder=3)
    ax.plot([r.q25, r.q75], [yi, yi], color=INK, lw=1.6, zorder=4, solid_capstyle="butt")
    ax.text(r.q75 + 0.14, yi, f"{2 ** r['median']:.0f}$\\times$", va="center", ha="left", fontsize=12.5,
            color=INK, fontweight="bold", zorder=5)
    if cls != "Stable":
        ax.text(0.12, yi, stars(p_vs_stable[cls]), va="center", ha="left", fontsize=11, color="#ffffff", zorder=6)
ax.set_yticks(y, [LAB[c] for c in ORDER])
ax.set_xlabel("opening of a closed site after the element is moved in  (log$_2$)")
ax.set_xlim(0, 10.1)
ax.set_xticks(np.arange(0, 6.1, 1))
ax.set_ylim(-0.62, len(ORDER) - 0.42)
ax.set_title("A   Stable elements open closed chromatin on their own", loc="left", fontsize=13.5, pad=10)
# how much of the gap is the donors' own accessibility rather than their sequence (53_fb_host_insertion_matched.py)
adj = pd.read_csv(SW / "host_insertion_accessibility_matched_measured.csv")
adj = adj[adj.width == "peak"].iloc[0]
ax.text(9.95, 1.0, "Stable donors also sit in\nmore accessible loci. Hold\nthat constant and a Stable\n"
        f"element still gains {adj.gap_adjusted:+.1f} log$_2$\n({2 ** adj.gap_adjusted:.1f}$\\times$), "
        f"p = {adj.p_adjusted:.3f}",
        ha="right", va="center", fontsize=10.5, color=INK, linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f4f2ee", edgecolor=GRID))

w = 0.36
xs = np.arange(2)
seq = [sl[sl.width == "170"].beta_donor_sequence.iloc[0], sl[sl.width == "peak"].beta_donor_sequence.iloc[0]]
loc = [sl[sl.width == "170"].beta_recipient_location.iloc[0], sl[sl.width == "peak"].beta_recipient_location.iloc[0]]
ax2.bar(xs - w / 2, seq, w, color=BLUE, label="what the element is", zorder=3)
ax2.bar(xs + w / 2, loc, w, color=CTRL_GREY, label="where it sits", zorder=3)
for x, v in zip(xs - w / 2, seq):
    ax2.text(x, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=11.5, color=INK)
for x, v in zip(xs + w / 2, loc):
    ax2.text(x, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=11.5, color=MUTED)
ax2.set_xticks(xs, ["170 bp element\n(the MPRA window)", "whole peak\n(median 0.8 kb)"])
ax2.set_ylabel("weight on the swapped site's accessibility")
ax2.set_ylim(0, 0.95)
ax2.legend(frameon=False, loc="upper left", fontsize=11)
ax2.set_title("B   Identity lives at peak scale, not in 170 bp", loc="left", fontsize=13.5, pad=10)
ax2.set_yticks(np.arange(0, 0.81, 0.2))
fig.text(0.005, -0.045,
         "A: 1,000 transplants per class into 20 closed backgrounds; bar = median, line = IQR; $\\times$ = fold rise over the same sequence shuffled; "
         "*** p < 0.001 vs Stable.\nB: regression of the swapped site's accessibility on the donor element and the recipient location, over 62,250 swaps "
         "between fibroblast AP-1 sites.", ha="left", va="top", fontsize=10, color=MUTED)
save(fig, "fb1_autonomy")

# ================================================================ 2. is the AP-1 deficit specific to AP-1?
PS = ROOT / "Analysis" / "fb_partner_shuffle" / "tables"
J = pd.read_csv(PS / "own_vs_partners_per_site.csv")
F = pd.read_csv(PS / "net_effect_by_family.csv")
GR = ["Stable", "AP-1 gain", "other Dynamic", "AP-1 loss"]
GLAB = {"Stable": "Stable", "AP-1 gain": "AP-1 gained\nwith age", "other Dynamic": "other\nDynamic",
        "AP-1 loss": "AP-1 lost\nwith age"}

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14.2, 5.4), gridspec_kw={"width_ratios": [1, 1.1], "wspace": 0.42})
x = np.arange(len(GR))
w = 0.36
for i, (field, col, lab) in enumerate((("own_ap1", ORANGE, "the site's own AP-1 motif"),
                                       ("partners_mean_net", CTRL_GREY, "its partner motifs, net of a matched control"))):
    off = (-w / 2) if i == 0 else (w / 2)
    med = np.array([J[J.grp == g][field].median() for g in GR])
    lo = np.array([boot(J[J.grp == g][field])[0] for g in GR])
    hi = np.array([boot(J[J.grp == g][field])[1] for g in GR])
    ax.bar(x + off, med, w, color=col, zorder=3, label=lab)
    ax.errorbar(x + off, med, yerr=[med - lo, hi - med], fmt="none", ecolor=INK, lw=1.4, capsize=3, zorder=4)
    for xi, v in zip(x + off, med):
        inside = i == 0
        ax.text(xi, v + (0.06 if inside else -0.10), f"{v:.2f}", ha="center", va="bottom" if inside else "top",
                fontsize=11.5, color="#ffffff" if inside else MUTED, zorder=6)
ax.axhline(0, color=MUTED, lw=1)
ax.set_xticks(x, [GLAB[g] for g in GR])
ax.set_ylabel("effect of shuffling the motif  (log$_2$)")
ax.set_ylim(-1.45, 0.78)
ax.set_yticks(np.arange(-1.25, 0.01, 0.25))
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.03), fontsize=10.5, handlelength=1.3)
ax.set_title("A   Only AP-1-losing sites let go — and only of AP-1", loc="left", fontsize=13.5, pad=10)
# every Dynamic class against Stable, on each of the two motif sets
ST = pd.read_csv(PS / "specificity_tests.csv").set_index("comparison")
for xi, g in zip(x, GR):
    if g == "Stable":
        continue
    row = ST.loc[f"Stable vs {g}"] if f"Stable vs {g}" in ST.index else None
    for off, field, col in ((-w / 2, "p_own_ap1", ORANGE), (w / 2, "p_partners", MUTED)):
        if row is None:
            from scipy.stats import mannwhitneyu
            f2 = "own_ap1" if field == "p_own_ap1" else "partners_mean_net"
            p = mannwhitneyu(J[J.grp == "Stable"][f2].dropna(), J[J.grp == g][f2].dropna()).pvalue
        else:
            p = row[field]
        ax.text(xi + off, 0.06, stars(p), ha="center", va="bottom", fontsize=12, color=col, fontweight="bold")
ax.text(0.0, 0.33, "vs Stable", fontsize=10.5, color=MUTED, ha="center", va="bottom", style="italic")

SHORT = {"AP-1 (Fos/Jun/Atf/Batf)": "AP-1, 2nd copy", "Heat shock / stress": "heat shock",
         "CTCF / cohesin": "CTCF", "E-box bHLH (Myod/Tcf/Neuro)": "E-box bHLH"}
top = F.sort_values("net").head(12).iloc[::-1]
yy = np.arange(len(top))
ax2.barh(yy, top.net, height=0.66, color=CTRL_GREY, zorder=3)
ax2.set_yticks(yy, [SHORT.get(f, f.split(" (")[0]) for f in top.family])
own = J.own_ap1.median()
ax2.axvline(own, color=ORANGE, lw=2.4, zorder=5)
ax2.text(own + 0.03, -0.95, "the site's own AP-1 motif", color=ORANGE, fontsize=11.5, va="center", ha="left",
         fontweight="bold")
for yi, (v, n) in enumerate(zip(top.net, top.n_sites)):
    ax2.text(0.02, yi, f"{n}", ha="left", va="center", fontsize=10, color=MUTED)
ax2.text(0.02, len(top) - 0.1, "sites", ha="left", va="center", fontsize=10, color=MUTED, style="italic")
ax2.set_xlabel("net effect of shuffling the partner motif  (log$_2$)")
ax2.set_xlim(-1.30, 0.17)
ax2.set_ylim(-1.5, len(top) - 0.3)
ax2.set_title("B   C/EBP is the strongest partner — still half of AP-1", loc="left", fontsize=13.5, pad=10)
fig.text(0.005, -0.03,
         "2,352 fibroblast AP-1 sites. Within each window the site's own AP-1 motif and each footprinted partner are shuffled separately; every partner is "
         "scored against a footprint-free\nwindow of the same width at the mirrored distance, so the comparison is within one site and needs no accessibility "
         "matching. Bars = median, 95% CI; *** p < 0.001, ns = not\nsignificant (Mann-Whitney vs Stable). AP-1 lost with age vs Stable: own AP-1 p = 2.3×10$^{-5}$, partner motifs p = 0.66.",
         ha="left", va="top", fontsize=10, color=MUTED)
save(fig, "fb2_ap1_specific")
print("done")
