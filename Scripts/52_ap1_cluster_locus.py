#!/usr/bin/env python3
"""A worked example of the top-window result: which 10 bp window does the model actually lean on?

Two elements with a comparable number of neighbouring AP-1 footprints. At the Stable one the strongest window is
its own motif; at the age-gained one it is a neighbouring AP-1 site. Emits two bare track panels (no titles, no
legend) into Analysis/slide_figures/panels/ so the slide can carry the words.
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
BLUE, ORANGE, GREY, GOLD = "#1F6FB2", "#C2551C", "#B3B2AE", "#E8B33C"
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
allnat = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "atac_170", "trackset"])
                    for f in sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))[:4]])
pc = 0.05 * allnat.query("kind=='native' and trackset=='CM'").atac_170.median()
nat = sc[sc.kind == "native"].groupby("site_id").atac_170.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
T = t.groupby(["site_id", "pos"], as_index=False).eff.mean()
T["centre"] = T.pos + 5
FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
ymin = T.eff.min() * 1.12

for sid, tag, col in SITES:
    w = T[T.site_id == sid].sort_values("centre")
    top3 = set(w.nsmallest(3, "eff").pos)
    ap1 = FP[(FP.site_id == sid) & (FP.family == "AP-1 (Fos/Jun/Atf/Batf)") & (~FP.own_ap1)]
    fig, (ax, axf) = plt.subplots(2, 1, figsize=(12.2, 2.25), sharex=True,
                                  gridspec_kw={"height_ratios": [3.4, 1.0], "hspace": 0.10})
    is_top = w.pos.isin(top3).values
    ax.bar(w.centre[~is_top], w.eff[~is_top], width=9, color=GREY, zorder=3)
    ax.bar(w.centre[is_top], w.eff[is_top], width=9, color=col, zorder=4)
    ax.axvspan(-5.5, 5.5, color=GOLD, alpha=0.45, lw=0, zorder=1)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    ax.set_xlim(-305, 305)
    ax.set_ylim(ymin, max(0.1, w.eff.max() * 1.3))
    ax.set_ylabel("Loss (log2FC)", fontsize=10)
    ax.tick_params(labelbottom=False)
    spans = []                                              # overlapping calls merge, so a cluster reads as one span
    for st, en in sorted(zip(ap1.rel_start, ap1.rel_end)):
        if spans and st <= spans[-1][1] + 6:
            spans[-1] = (spans[-1][0], max(spans[-1][1], en), spans[-1][2] + 1)
        else:
            spans.append((st, en, 1))
    for st, en, n in spans:
        on = any(st - 10 < p + 5 < en + 10 for p in top3)
        axf.add_patch(plt.Rectangle((st - 3, -0.42), max(en - st, 12) + 6, 0.62,
                                    color=col if on else GREY, ec=INK if on else "none", lw=1.1, zorder=3))
        if n > 1:
            axf.text((st + en) / 2, -0.11, str(n), ha="center", va="center", fontsize=8,
                     color="white", fontweight="bold", zorder=5)
    axf.add_patch(plt.Rectangle((-7, -0.42), 14, 0.62, color=GOLD, ec=INK, lw=1.1, zorder=4))
    axf.set_ylim(-0.75, 0.35)
    axf.set_yticks([])
    axf.grid(False)
    axf.set_ylabel("AP-1\nfootprints", fontsize=8.5, color=MUTED)
    axf.set_xlabel("Position relative to the called AP-1 motif (bp)", fontsize=9.5)
    for sp in ("left", "bottom"):
        axf.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.072, right=0.995, top=0.985, bottom=0.185)
    fig.savefig(PAN / f"s11_cluster_{tag}.png", facecolor=SURF)
    plt.close(fig)
    own = w.loc[w.centre.abs() <= 6, "eff"].min()
    best = w.loc[w.centre.abs() > 6, "eff"].min()
    print(f"{tag:8s} {sid:26s} own motif {own:+.2f}   strongest neighbour {best:+.2f}   "
          f"neighbouring AP-1 footprints {len(ap1)}   top-3 windows at {sorted(p + 5 for p in top3)}")
