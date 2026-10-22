#!/usr/bin/env python3
"""Heatmaps of which TF families the shuffle screen marks as important, per class, and which ones separate the classes in each direction.
Importance = median effect of replacing a 10 bp tile overlapping that family's footprint, measured against the footprint-free tiles of the SAME site, one value per site per family.
Most per-class cells are small, so the difference panels are explicitly exploratory: the top 5 families in each direction are labelled whether or not they reach significance.
Outputs Analysis/motif_heatmap/{tables,figures}.   python3 45_motif_importance_heatmap.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "motif_heatmap"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
INK, MUTED, SURF, GRID, BLUE, ORANGE, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "savefig.dpi": 200})
MIN_SITES = 8                                        # a family needs this many sites in BOTH arms of a contrast to be shown

D = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
GRP = {"dynamic": "Dynamic (thr.)", "matched_stable": "Stable, matched (thr.)", "high_stable": "Stable, high acc."}
AGE = {"Age Up (score)": "Age Up (\u0394)", "Age Down (score)": "Age Down (\u0394)", "Unchanged": "Unchanged (\u0394)"}
D["cls_group"] = D.group.map(GRP)
D["cls_age"] = D.age_score_class.map(AGE)

# ---------------------------------------------------------------- absolute importance per class
panels = []
for col, order in (("cls_group", list(GRP.values())), ("cls_age", list(AGE.values()))):
    piv = D.pivot_table(index="family", columns=col, values="vs_free", aggfunc="median").reindex(columns=order)
    cnt = D.pivot_table(index="family", columns=col, values="vs_free", aggfunc="size").reindex(columns=order).fillna(0)
    piv = piv.where(cnt >= MIN_SITES)
    panels.append((col, piv, cnt))
ABS = panels[0][1].join(panels[1][1], how="outer")
ABS = ABS[ABS.notna().sum(1) >= 4]                   # drop families too sparse to compare across classes
ABS = ABS.loc[ABS.mean(1).sort_values().index]
ABS.to_csv(OUT / "tables" / "importance_by_class.csv")

# ---------------------------------------------------------------- differences, with the top 5 each way
CONTRASTS = [("Dynamic vs\nStable (thr.)", "group", "dynamic", "matched_stable"),
             ("Dynamic vs\nhigh-acc. Stable", "group", "dynamic", "high_stable"),
             ("Age Up vs\nUnchanged (Δ)", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down vs\nUnchanged (Δ)", "age_score_class", "Age Down (score)", "Unchanged")]
rows = []
for lab, col, a_, b_ in CONTRASTS:
    for fam, w in D.groupby("family"):
        x, y = w[w[col] == a_].vs_free.dropna(), w[w[col] == b_].vs_free.dropna()
        if min(len(x), len(y)) < MIN_SITES:
            continue
        rows.append({"contrast": lab, "family": fam, "n_a": len(x), "n_b": len(y),
                     "median_a": x.median(), "median_b": y.median(), "difference": x.median() - y.median(),
                     "p": mannwhitneyu(x, y).pvalue})
DIFF = pd.DataFrame(rows)
DIFF["rank_in_contrast"] = DIFF.groupby("contrast").difference.rank()
tops = []
for lab, w in DIFF.groupby("contrast"):
    w = w.sort_values("difference")
    for d_, lab2 in ((w.head(5), "more important in the first class"), (w.tail(5).iloc[::-1], "more important in the second class")):
        for r in d_.itertuples():
            tops.append({"contrast": lab, "direction": lab2, "family": r.family, "n_first": r.n_a, "n_second": r.n_b,
                         "median_first": r.median_a, "median_second": r.median_b, "difference": r.difference, "p": r.p})
TOP = pd.DataFrame(tops)
TOP.to_csv(OUT / "tables" / "top5_each_direction.csv", index=False)
DIFF.to_csv(OUT / "tables" / "importance_differences.csv", index=False)

# ---------------------------------------------------------------- figure
fams = ABS.index.tolist()
fig = plt.figure(figsize=(12.4, max(6.0, 0.31 * len(fams) + 3.0)))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.06)
ax = fig.add_subplot(gs[0])
v = ABS.values
norm = TwoSlopeNorm(vmin=np.nanmin(v), vcenter=0, vmax=max(np.nanmax(v), 0.05))
im = ax.imshow(v, cmap="RdBu_r", norm=norm, aspect="auto")
ax.set_xticks(range(ABS.shape[1]), ABS.columns, fontsize=7.5, rotation=40, ha="right")
ax.set_yticks(range(len(fams)), fams, fontsize=7.5)
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        if np.isfinite(v[i, j]):
            ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center", fontsize=6.2,
                    color="white" if abs(v[i, j]) > 0.45 * max(abs(np.nanmin(v)), abs(np.nanmax(v))) else INK)
ax.set_title("Importance of each factor, by class\nmedian log2FC of a tile over its footprint", fontsize=9, loc="left")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.16, location="bottom", shrink=0.9)
cb.set_label("more important  ←      →  less important", fontsize=7.5)
cb.ax.tick_params(labelsize=7)

ax = fig.add_subplot(gs[1])
cs = [c[0] for c in CONTRASTS]
MD = DIFF.pivot_table(index="family", columns="contrast", values="difference").reindex(index=fams, columns=cs)
m = MD.values
lim = np.nanmax(np.abs(m)) if np.isfinite(m).any() else 1
im2 = ax.imshow(m, cmap="PuOr", vmin=-lim, vmax=lim, aspect="auto")
ax.set_xticks(range(len(cs)), [c.replace(" vs ", "\nvs ") for c in cs], fontsize=7.5)
ax.set_yticks(range(len(fams)), ["" for _ in fams])
top_set = {(r.contrast, r.family) for r in TOP.itertuples()}
for i, fam in enumerate(fams):
    for j, c in enumerate(cs):
        if not np.isfinite(m[i, j]):
            continue
        txt = f"{m[i, j]:+.2f}"
        is_top = (c, fam) in top_set
        ax.text(j, i, txt, ha="center", va="center", fontsize=6.2, fontweight="bold" if is_top else "normal",
                color="white" if abs(m[i, j]) > 0.55 * lim else INK)
        if is_top:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.1))
ax.set_title("Difference between classes (first minus second)\nboxed = top 5 each way", fontsize=9, loc="left")
cb2 = fig.colorbar(im2, ax=ax, fraction=0.04, pad=0.16, location="bottom", shrink=0.9)
cb2.set_label("more important in the FIRST class  ←     →  in the SECOND", fontsize=7.5)
cb2.ax.tick_params(labelsize=7)
fig.suptitle("Which factors does the shuffle screen mark as important, and for which class?\nEach value is measured against that site's own tiles with no footprint. Negative = replacing the tile closes the site, so the factor matters.",
             fontsize=10.5, x=0.012, ha="left", y=0.998)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig(OUT / "figures" / "x1_motif_importance_heatmap.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
print("== importance by class (blank = fewer than %d sites)\n" % MIN_SITES, ABS.round(3).to_string())
print("\n== top 5 each direction\n", TOP.round(3).to_string(index=False))
