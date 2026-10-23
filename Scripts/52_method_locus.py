#!/usr/bin/env python3
"""A demonstration of the method at one real element: the sliding-window scan and the footprints it is read against.

Not an argument -- there is no class contrast here. One locus, the called AP-1 motif, what the prediction loses as
each 10 bp window is replaced, and every footprinted motif in the element labelled by family. All prose lives on the
slide, so the figure carries only axes, data and the footprint labels.
  python3 52_method_locus.py [--site ID]
Writes Analysis/slide_figures/panels/s11_method_locus.png
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
INK, MUTED, SURF, GRID, GREY, GOLD = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df", "#B3B2AE", "#C9961F"
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 300})
FCOL = {"MEF2": "#C0392B", "GATA": "#8E44AD", "NFI (Nfia/b/c/x)": "#16A085", "TEAD": "#D68910",
        "NKX / homeobox cardiac": "#2980B9", "TBX": "#7F8C8D", "CTCF / cohesin": "#2C3E50",
        "SRF / MADS": "#E67E22", "AP-1 (Fos/Jun/Atf/Batf)": "#555555", "RUNX": "#95A5A6",
        "ETS": "#95A5A6", "KLF / SP": "#95A5A6"}
ap = argparse.ArgumentParser()
ap.add_argument("--site", default="chr14-77485003-77485014")
a = ap.parse_args()

sc = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "pos", "atac_170", "trackset"])
                for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[(sc.trackset == "CM") & (sc.site_id == a.site)]
pc = 0.05 * pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "atac_170", "trackset"])
                       for f in sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))[:4]]
                      ).query("kind=='native' and trackset=='CM'").atac_170.median()
nat = sc[sc.kind == "native"].atac_170.median()
w = sc[sc.kind == "tile"].copy()
w["eff"] = np.log2((w.atac_170 + pc) / (nat + pc))
w = w.groupby("pos", as_index=False).eff.mean()
w["centre"] = w.pos + 5
FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
fp = FP[(FP.site_id == a.site) & (~FP.own_ap1) & (FP.family.isin(FCOL)) &
        (((FP.rel_start + FP.rel_end) / 2).abs() < 290)].copy()

fig, (ax, axf) = plt.subplots(2, 1, figsize=(12.0, 3.1), sharex=True,
                              gridspec_kw={"height_ratios": [3.0, 1.25], "hspace": 0.0})
on_motif = w.centre.abs() <= 6
ax.bar(w.centre[~on_motif], w.eff[~on_motif], width=9, color=GREY, zorder=3)
ax.bar(w.centre[on_motif], w.eff[on_motif], width=9, color=GOLD, ec=INK, lw=0.8, zorder=4)
ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
ax.set_xlim(-305, 305)
ax.set_ylim(w.eff.min() * 1.14, max(0.12, w.eff.max() * 1.5))
ax.set_ylabel("Loss of predicted\naccessibility (log2FC)", fontsize=10)
ax.tick_params(labelbottom=False)

lanes: list[float] = []                                        # stack footprints that would otherwise overlap
for r in fp.sort_values("rel_start").itertuples():
    k = next((i for i, end in enumerate(lanes) if r.rel_start > end + 42), len(lanes))
    nm = r.family.split(" (")[0].split(" /")[0]
    wid = max(r.rel_end - r.rel_start, 7.6 * len(nm) + 10)
    cx = (r.rel_start + r.rel_end) / 2
    if k == len(lanes):
        lanes.append(cx + wid / 2)
    else:
        lanes[k] = cx + wid / 2
    c = FCOL[r.family]
    axf.add_patch(plt.Rectangle((cx - wid / 2, -k - 0.40), wid, 0.62, color=c, ec="none", zorder=3))
    axf.text(cx, -k - 0.09, nm, ha="center", va="center", fontsize=8.2, color="white", fontweight="bold", zorder=4)
axf.add_patch(plt.Rectangle((-26, 0.35), 52, 0.56, color=GOLD, ec=INK, lw=1.0, zorder=4))
axf.text(0, 0.63, "AP-1", ha="center", va="center", fontsize=8.6, color=INK, fontweight="bold", zorder=5)
axf.set_ylim(-max(len(lanes), 1) - 0.1, 1.05)
axf.set_yticks([])
axf.grid(False)
axf.set_ylabel("footprinted\nmotifs", fontsize=9, color=MUTED)
axf.set_xlabel("Position relative to the AP-1 motif (bp)", fontsize=10.5)
for sp in ("left", "bottom"):
    axf.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.088, right=0.995, top=0.985, bottom=0.155)
fig.savefig(PAN / "s11_method_locus.png", facecolor=SURF)
plt.close(fig)
print(f"{a.site}: {len(w)} windows, motif window {w.eff[on_motif].min():+.2f}, "
      f"{len(fp)} footprints of {fp.family.nunique()} families in {len(lanes)} lanes")
