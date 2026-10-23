#!/usr/bin/env python3
"""The density of importance across the 60 sliding windows, class by class.

The summary measures (loss within +-85 bp, loss beyond it, spread) are read off this curve; this is the curve.
Two versions, because they answer different questions:
  absolute   -- median loss at each window position, so depth and shape are mixed
  normalised -- each site's profile divided by its own total loss, so only the SHAPE is compared
Outputs Analysis/positional_profile/ and bare panels for the deck.
  python3 53_positional_profile.py
"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "positional_profile"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
PAN = ROOT / "Analysis" / "slide_figures" / "panels"
PAN.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(53)
INK, MUTED, SURF, GRID, CTRL = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df", "#B3B2AE"
ORANGE, BLUE, GREEN_D, PURPLE_D, DARK = "#C2551C", "#1F6FB2", "#2E7D32", "#7B3FA0", "#4F4F4F"
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 300})

sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
sc = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "pos", "atac_170", "trackset"])
                for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[sc.trackset == "CM"]
nat = sc[sc.kind == "native"].groupby("site_id").atac_170.median()
pc = 0.05 * nat.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
T = t.groupby(["site_id", "pos"], as_index=False).eff.mean()
T["centre"] = T.pos + 5
T["loss"] = np.clip(-T.eff, 0, None)
tot = T.groupby("site_id").loss.sum()
T["share"] = T.loss / T.site_id.map(tot).replace(0, np.nan)

GROUPS = [("Stable", "stable_dynamic", "Stable", CTRL), ("Dynamic (all)", "stable_dynamic", "Dynamic", DARK),
          ("Age Up", "age_class", "Age Up", ORANGE), ("Age Down", "age_class", "Age Down", BLUE),
          ("Exercise Up", "ex_class", "Exercise Up", GREEN_D), ("Exercise Down", "ex_class", "Exercise Down", PURPLE_D)]
rows = []
for nm, col, val, _c in GROUPS:
    ids = set(sl.index[sl[col] == val])
    w = T[T.site_id.isin(ids)]
    for metric in ("loss", "share"):
        # the MEAN, not the median: summed over the windows it equals the per-site totals the dot panels show,
        # so the area under the curve and the summary measures are the same quantity
        g = w.groupby("centre")[metric].mean()
        rows += [{"class": nm, "metric": metric, "centre": c, "mean": v, "n_sites": len(ids)} for c, v in g.items()]
P = pd.DataFrame(rows)
P.to_csv(OUT / "tables" / "positional_profile_by_class.csv", index=False)

# the core-vs-flank test the curve is a picture of
core, flank = T.centre.abs() <= 15, T.centre.abs() > 85
stat = []
for nm, col, val, _c in GROUPS[1:]:
    ids = set(sl.index[sl[col] == val])
    sids = set(sl.index[sl.stable_dynamic == "Stable"])
    for lab, m in (("core, |x| <= 15 bp", core), ("flank, |x| > 85 bp", flank)):
        for metric in ("loss", "share"):
            a_ = T[m & T.site_id.isin(ids)].groupby("site_id")[metric].mean()
            b_ = T[m & T.site_id.isin(sids)].groupby("site_id")[metric].mean()
            stat.append({"class": nm, "region": lab, "metric": metric, "class_median": a_.median(),
                         "stable_median": b_.median(), "ratio": a_.median() / b_.median(),
                         "p": mannwhitneyu(a_.dropna(), b_.dropna()).pvalue})
S_ = pd.DataFrame(stat)
o = np.argsort(S_.p.values); r = np.empty(len(S_)); r[o] = np.arange(1, len(S_) + 1)
S_["fdr"] = np.minimum(S_.p.values * len(S_) / r, 1)
S_.to_csv(OUT / "tables" / "core_vs_flank_by_class.csv", index=False)

for metric, ylab, fname in (("loss", "Mean loss per window (log2FC)", "s12_profile_absolute.png"),
                            ("share", "Mean share of the site's total", "s12_profile_shape.png")):
    fig, ax = plt.subplots(figsize=(7.8, 2.45))
    for nm, _col, _val, c in GROUPS:
        g = P[(P["class"] == nm) & (P.metric == metric)].sort_values("centre")
        y = g["mean"].rolling(3, center=True, min_periods=1).mean()     # 3-window smoothing; the tails are noisy
        ax.plot(g.centre, y, color=c, lw=2.4 if nm == "Stable" else 1.9, zorder=4 if nm == "Stable" else 3)
    ax.axvspan(-15, 15, color="#C9961F", alpha=0.13, lw=0, zorder=1)
    ax.set_xlim(-305, 305)
    ax.set_xlabel("Window position relative to the motif (bp)", fontsize=9.5)
    ax.set_ylabel(ylab, fontsize=9.5)
    fig.subplots_adjust(left=0.135, right=0.975, top=0.97, bottom=0.195)
    fig.savefig(PAN / fname, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)

pd.set_option("display.width", 210)
print(S_.round(4).to_string(index=False))
print("\nwrote", OUT / "tables", "and 2 panels")
