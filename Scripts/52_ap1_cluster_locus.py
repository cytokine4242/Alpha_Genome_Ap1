#!/usr/bin/env python3
"""A worked example of the top-window odds ratio: which 10 bp window does the model actually lean on?

The statistic being illustrated is "does one of this element's three most damaging windows cover a NEIGHBOURING
AP-1 footprint". So the figure has to make four things readable without a caption: where the called motif is,
where the other AP-1 footprints are, which windows rank 1-2-3, and whether those windows land on the called motif
or on a neighbour. Every footprint is labelled in the lane, each top window is badged with its rank, and a guide
line is drawn from each top window down to whatever sits beneath it, so the link is drawn rather than implied.
Both panels share a y-axis so the bar heights are comparable.
  python3 52_ap1_cluster_locus.py [--stable ID --dynamic ID]
"""
from __future__ import annotations

import argparse
import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
PAN = ROOT / "Analysis" / "slide_figures" / "panels"
PAN.mkdir(parents=True, exist_ok=True)
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df"
BLUE, ORANGE, GREY, GOLD = "#1F6FB2", "#C2551C", "#B3B2AE", "#C9961F"
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 300})
ap = argparse.ArgumentParser()
ap.add_argument("--stable", default="chr2-110186562-110186573")
ap.add_argument("--dynamic", default="chr13-35395938-35395949")
a = ap.parse_args()
SITES = [(a.stable, "stable", BLUE), (a.dynamic, "dynamic", ORANGE)]

sc = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "pos", "atac_170", "trackset"])
                for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[(sc.trackset == "CM") & sc.site_id.isin([s for s, _, _ in SITES])]
pc = 0.05 * pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "atac_170", "trackset"])
                       for f in sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))[:4]]
                      ).query("kind=='native' and trackset=='CM'").atac_170.median()
nat = sc[sc.kind == "native"].groupby("site_id").atac_170.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
T = t.groupby(["site_id", "pos"], as_index=False).eff.mean()
T["centre"] = T.pos + 5
FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
YLIM = (T.eff.min() * 1.12, 0.62)                               # shared, with headroom above zero for the rank badges

for sid, tag, col in SITES:
    w = T[T.site_id == sid].sort_values("centre")
    ranked = w.nsmallest(3, "eff").reset_index(drop=True)        # 1 = most damaging
    rank_of = {r.pos: i + 1 for i, r in ranked.iterrows()}
    ap1 = FP[(FP.site_id == sid) & (FP.family == "AP-1 (Fos/Jun/Atf/Batf)") & (~FP.own_ap1)]
    spans = []                                                  # overlapping calls merge, so a cluster reads as one
    for st, en in sorted(zip(ap1.rel_start, ap1.rel_end)):
        if spans and st <= spans[-1][1] + 6:
            spans[-1] = (spans[-1][0], max(spans[-1][1], en), spans[-1][2] + 1)
        else:
            spans.append((st, en, 1))

    fig, (ax, axf) = plt.subplots(2, 1, figsize=(12.6, 2.1), sharex=True,
                                  gridspec_kw={"height_ratios": [3.0, 0.9], "hspace": 0.0})
    is_top = w.pos.isin(rank_of).values
    ax.bar(w.centre[~is_top], w.eff[~is_top], width=9, color=GREY, zorder=3)
    ax.bar(w.centre[is_top], w.eff[is_top], width=9, color=col, zorder=4)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    for r in ranked.itertuples():                               # badges ride above zero, so they never collide with a bar
        ax.plot([r.centre, r.centre], [r.eff, YLIM[0]], color=col, lw=1.0, ls=(0, (2, 2)), alpha=0.75, zorder=1)
        ax.text(r.centre, 0.33, str(rank_of[r.pos]), ha="center", va="center", fontsize=9.5,
                color="white", fontweight="bold", zorder=6,
                bbox=dict(boxstyle="circle,pad=0.24", fc=col, ec="none"))
    ax.set_ylim(*YLIM)
    ax.set_xlim(-305, 305)
    ax.set_ylabel("Loss when replaced\n(log2FC)", fontsize=9.5)
    ax.tick_params(labelbottom=False)

    for st, en, n in spans:                                     # neighbours, labelled in place
        on = any(st - 10 < c < en + 10 for c in ranked.centre)
        x0, wd = st - 3, max(en - st, 12) + 6
        axf.add_patch(plt.Rectangle((x0, 0.18), wd, 0.46, color=col if on else GREY,
                                    ec=INK if on else "none", lw=1.1, zorder=3))
        axf.text(x0 + wd / 2, -0.04, "AP-1" + (f" ×{n}" if n > 1 else ""), ha="center", va="top",
                 fontsize=8.5, color=col if on else MUTED, fontweight="bold" if on else "normal")
    axf.add_patch(plt.Rectangle((-8, 0.18), 16, 0.46, color=GOLD, ec=INK, lw=1.2, zorder=4))
    axf.text(0, -0.04, "called AP-1 motif", ha="center", va="top", fontsize=8.5, color=GOLD, fontweight="bold")
    axf.set_ylim(-0.75, 0.80)
    axf.set_yticks([])
    axf.grid(False)
    axf.set_xlabel("Position relative to the called AP-1 motif (bp)", fontsize=9.5)
    for sp in ("left", "bottom"):
        axf.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.085, right=0.995, top=0.985, bottom=0.21)
    fig.savefig(PAN / f"s11_cluster_{tag}.png", facecolor=SURF)
    plt.close(fig)
    own = w.loc[w.centre.abs() <= 6, "eff"].min()
    best = w.loc[w.centre.abs() > 6, "eff"].min()
    where = ["the called motif" if abs(c) <= 10 else "a neighbouring AP-1 site" if any(st - 10 < c < en + 10 for st, en, _ in spans) else "neither"
             for c in ranked.centre]
    print(f"{tag:8s} {sid:26s} own {own:+.2f}  best neighbour {best:+.2f}  {len(ap1)} neighbouring AP-1 calls in {len(spans)} spans")
    print(f"         top 3 windows: " + ", ".join(f"#{i+1} at {int(c):+d} bp on {wh}" for i, (c, wh) in enumerate(zip(ranked.centre, where))))
