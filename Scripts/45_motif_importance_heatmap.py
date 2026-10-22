#!/usr/bin/env python3
"""Heatmaps of which TF families the shuffle screen marks as important, broken down by DIRECTION of change rather than by Stable/Dynamic.
Two parallel sets of classes, kept apart on purpose:
  threshold (binary call): Stable, Age Up, Age Down, Exercise Up, Exercise Down   -- from in_YS / in_MAS / in_MAE
  delta score (graded):    Unchanged, Age Up, Age Down                            -- top and bottom fifth of the score change
Importance = median effect of replacing a 10 bp tile overlapping that family's footprint, against the footprint-free tiles of the SAME site; one value per site per family.
The directional classes are small (Age Down 29 sites, Exercise Up 22), so every cell carries its n and the whole figure is exploratory.
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
MIN_SITES = 5                                        # a cell needs this many sites to be shown at all

D = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
D = D.join(sl[["age_class", "ex_class"]], on="site_id")
D["age_class"] = D.age_class.replace({"Neither YS nor MAS (= MAE only)": np.nan})
D["ex_class"] = D.ex_class.replace({"Neither MAS nor MAE (= YS only)": np.nan})

COLS = [("Stable\n(thr.)", "age_class", "Stable"),
        ("Age Up\n(thr.)", "age_class", "Age Up"),
        ("Age Down\n(thr.)", "age_class", "Age Down"),
        ("Ex Up\n(thr.)", "ex_class", "Exercise Up"),
        ("Ex Down\n(thr.)", "ex_class", "Exercise Down"),
        ("Unchanged\n(Δ)", "age_score_class", "Unchanged"),
        ("Age Up\n(Δ)", "age_score_class", "Age Up (score)"),
        ("Age Down\n(Δ)", "age_score_class", "Age Down (score)")]
med = pd.DataFrame({lab: D[D[col] == val].groupby("family").vs_free.median() for lab, col, val in COLS})
cnt = pd.DataFrame({lab: D[D[col] == val].groupby("family").vs_free.size() for lab, col, val in COLS}).fillna(0)
med = med.where(cnt >= MIN_SITES)
keep = med.notna().sum(1) >= 4
med, cnt = med[keep], cnt[keep]
med = med.loc[med.mean(1).sort_values().index]
cnt = cnt.loc[med.index]
med.to_csv(OUT / "tables" / "importance_by_class.csv")
cnt.to_csv(OUT / "tables" / "site_counts_by_class.csv")

CONTRASTS = [("Age Up vs Stable\n(threshold)", "age_class", "Age Up", "Stable"),
             ("Age Down vs Stable\n(threshold)", "age_class", "Age Down", "Stable"),
             ("Ex Up vs Stable\n(threshold)", "ex_class", "Exercise Up", "Stable"),
             ("Ex Down vs Stable\n(threshold)", "ex_class", "Exercise Down", "Stable"),
             ("Age Up vs Unchanged\n(Δscore)", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down vs Unchanged\n(Δscore)", "age_score_class", "Age Down (score)", "Unchanged")]
rows = []
for lab, col, a_, b_ in CONTRASTS:
    for fam, w in D.groupby("family"):
        x, y = w[w[col] == a_].vs_free.dropna(), w[w[col] == b_].vs_free.dropna()
        if len(x) < MIN_SITES or len(y) < MIN_SITES:
            continue
        rows.append({"contrast": lab.replace("\n", " "), "family": fam, "n_class": len(x), "n_reference": len(y),
                     "median_class": x.median(), "median_reference": y.median(), "difference": x.median() - y.median(),
                     "p": mannwhitneyu(x, y).pvalue})
DIFF = pd.DataFrame(rows)
DIFF.to_csv(OUT / "tables" / "importance_differences.csv", index=False)
tops = []
for lab, w in DIFF.groupby("contrast"):
    w = w.sort_values("difference")
    for d_, dlab in ((w.head(5), "more important in this class"), (w.tail(5).iloc[::-1], "more important in the reference")):
        for r in d_.itertuples():
            tops.append({"contrast": lab, "direction": dlab, "family": r.family, "n_class": r.n_class, "n_reference": r.n_reference,
                         "median_class": r.median_class, "median_reference": r.median_reference, "difference": r.difference, "p": r.p})
TOP = pd.DataFrame(tops)
TOP.to_csv(OUT / "tables" / "top5_each_direction.csv", index=False)

# ---------------------------------------------------------------- figure
fams = med.index.tolist()
fig = plt.figure(figsize=(14.6, max(6.0, 0.33 * len(fams) + 3.2)))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.06)

ax = fig.add_subplot(gs[0])
v = med.values
vlo = max(np.nanpercentile(v, 2), -0.8)                # clip the scale: one thin cell (n = 6) otherwise washes out the rest
norm = TwoSlopeNorm(vmin=vlo, vcenter=0, vmax=max(np.nanmax(v), 0.05))
im = ax.imshow(v, cmap="RdBu_r", norm=norm, aspect="auto")
ax.set_xticks(range(med.shape[1]), med.columns, fontsize=7.5, rotation=0)
ax.set_yticks(range(len(fams)), fams, fontsize=7.5)
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        if np.isfinite(v[i, j]):
            big = v[i, j] < 0.55 * vlo or v[i, j] > 0.55 * max(np.nanmax(v), 0.05)
            ax.text(j, i, f"{v[i, j]:.2f}\n(n={int(cnt.values[i, j])})", ha="center", va="center", fontsize=5.6,
                    color="white" if big else INK, linespacing=1.0)
ax.axvline(4.5, color=INK, lw=1.6)
ax.set_title("Importance of each factor, by class\nmedian log2FC of a tile over its footprint", fontsize=9, loc="left")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.14, location="bottom", shrink=0.9)
cb.set_label("more important  ←      →  less important   (colour clipped; printed values are exact)", fontsize=7.5)
cb.ax.tick_params(labelsize=7)

ax = fig.add_subplot(gs[1])
cs = [c[0] for c in CONTRASTS]
MD = DIFF.pivot_table(index="family", columns="contrast", values="difference")
MD = MD.reindex(index=fams, columns=[c.replace("\n", " ") for c in cs])
m = MD.values
lim = min(np.nanpercentile(np.abs(m), 96), 0.6) if np.isfinite(m).any() else 1      # same clipping on the difference panel
im2 = ax.imshow(m, cmap="PuOr", vmin=-lim, vmax=lim, aspect="auto")
ax.set_xticks(range(len(cs)), cs, fontsize=7, rotation=20, ha="right")
ax.set_yticks(range(len(fams)), ["" for _ in fams])
top_set = {(r.contrast, r.family) for r in TOP.itertuples()}
sig_set = {(r.contrast, r.family) for r in DIFF.itertuples() if r.p < 0.05}
for i, fam in enumerate(fams):
    for j, c in enumerate([c.replace("\n", " ") for c in cs]):
        if not np.isfinite(m[i, j]):
            continue
        is_top, is_sig = (c, fam) in top_set, (c, fam) in sig_set
        ax.text(j, i, f"{m[i, j]:+.2f}" + ("*" if is_sig else ""), ha="center", va="center", fontsize=6.2,
                fontweight="bold" if is_top else "normal", color="white" if abs(m[i, j]) > 0.55 * lim else INK)
        if is_top:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.1))
ax.axvline(3.5, color=INK, lw=1.6)
ax.set_title("Difference from the reference class\nboxed = top 5 each way, * = p < 0.05 uncorrected", fontsize=9, loc="left")
cb2 = fig.colorbar(im2, ax=ax, fraction=0.04, pad=0.14, location="bottom", shrink=0.9)
cb2.set_label("more important in this class  ←     →  in the reference   (colour clipped)", fontsize=7.5)
cb2.ax.tick_params(labelsize=7)

fig.suptitle("Which factors the shuffle screen marks as important, by DIRECTION of change\n"
             "Left of each black line: threshold (binary call) classes against Stable. Right: Δscore classes against Unchanged. "
             "Negative = replacing that tile closes the site, so the factor matters.", fontsize=10.5, x=0.012, ha="left", y=0.999)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(OUT / "figures" / "x1_motif_importance_heatmap.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
u = D.drop_duplicates("site_id")
print("== scanned sites per class\n", u.age_class.value_counts().to_string(), "\n", u.ex_class.value_counts().to_string(),
      "\n", u.age_score_class.value_counts().to_string())
print("\n== importance by class (blank = fewer than %d sites)\n" % MIN_SITES, med.round(3).to_string())
print("\n== top 5 each direction\n", TOP.round(3).to_string(index=False))
