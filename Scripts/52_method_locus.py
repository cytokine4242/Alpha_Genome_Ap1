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
# every named family is shown; only the two catch-all bins are dropped, so the lane is the real footprint density
DROP = {"Other", "Zinc finger, other"}
FCOL = {"MEF2": "#C0392B", "GATA": "#8E44AD", "NFI (Nfia/b/c/x)": "#16A085", "TEAD": "#D68910",
        "NKX / homeobox cardiac": "#2980B9", "TBX": "#7F8C8D", "CTCF / cohesin": "#2C3E50",
        "SRF / MADS": "#E67E22", "AP-1 (Fos/Jun/Atf/Batf)": "#555555", "RUNX": "#8E6E53",
        "ETS": "#B9770E", "KLF / SP": "#7D6608", "MYC / MAX": "#117A65", "NF-kB / Rel": "#A04000",
        "E-box bHLH (Myod/Tcf/Neuro)": "#5B2C6F", "STAT / IRF": "#1A5276", "SOX": "#148F77",
        "FOX": "#7E5109", "RFX": "#4A235A", "NFAT": "#873600", "Smad / Tgfb": "#196F3D",
        "TALE homeodomain (Meis/Pknox/Tgif/Pbx)": "#633974", "Nuclear receptor (Esrr/Nr/Ppar/Thr)": "#145A32",
        "Heat shock / stress": "#6E2C00", "CEBP": "#512E5F", "YY / Gli / Pax": "#154360"}
ap = argparse.ArgumentParser()
ap.add_argument("--site", default="chr2-167698004-167698015")
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
fp = FP[(FP.site_id == a.site) & (~FP.own_ap1) & (~FP.family.isin(DROP)) &
        (((FP.rel_start + FP.rel_end) / 2).abs() < 290)].copy()
merged = []                                                     # overlapping calls of one family are one footprint
for fam, g in fp.groupby("family"):
    for st, en in sorted(zip(g.rel_start, g.rel_end)):
        if merged and merged[-1][0] == fam and st <= merged[-1][2] + 4:
            merged[-1] = (fam, merged[-1][1], max(merged[-1][2], en), merged[-1][3] + 1)
        else:
            merged.append((fam, st, en, 1))
fp = pd.DataFrame(merged, columns=["family", "rel_start", "rel_end", "n_calls"])

fig, (ax, axf) = plt.subplots(2, 1, figsize=(12.0, 3.5), sharex=True,
                              gridspec_kw={"height_ratios": [2.7, 1.6], "hspace": 0.0})
on_motif = w.centre.abs() <= 6
ax.bar(w.centre[~on_motif], w.eff[~on_motif], width=9, color=GREY, zorder=3)
ax.bar(w.centre[on_motif], w.eff[on_motif], width=9, color=GOLD, ec=INK, lw=0.8, zorder=4)
ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
ax.set_xlim(-305, 305)
ax.set_ylim(w.eff.min() * 1.14, max(0.12, w.eff.max() * 1.5))
ax.set_ylabel("Loss of predicted\naccessibility (log2FC)", fontsize=10)
ax.tick_params(labelbottom=False)

# Footprints are drawn at their true extent -- 6 to 39 bp, median 11 -- so nothing is inflated to fit its label.
# Labels sit in stacked rows beneath the rail with a leader line back to the chip they belong to.
RAILY, ROW = 0.0, 0.30
label_rows: list[float] = []
for r in fp.sort_values(["rel_start", "rel_end"]).itertuples():
    nm = r.family.split(" (")[0].split(" /")[0] + (f" \u00d7{r.n_calls}" if r.n_calls > 1 else "")
    cx = (r.rel_start + r.rel_end) / 2
    half = 3.9 * len(nm)                                        # half the label's width in bp at this scale
    k = next((i for i, end in enumerate(label_rows) if cx - half > end + 8), len(label_rows))
    if k == len(label_rows):
        label_rows.append(cx + half)
    else:
        label_rows[k] = cx + half
    c = FCOL[r.family]
    y = RAILY - 0.22 - k * ROW
    axf.add_patch(plt.Rectangle((r.rel_start, RAILY - 0.07), max(r.rel_end - r.rel_start, 3), 0.20,
                                color=c, ec="none", zorder=4))
    axf.plot([cx, cx], [RAILY - 0.08, y + 0.07], color=c, lw=0.7, alpha=0.8, zorder=3)
    axf.text(cx, y, nm, ha="center", va="top", fontsize=7.8, color=c, fontweight="bold", zorder=5)
axf.plot([-305, 305], [RAILY, RAILY], color=MUTED, lw=0.9, zorder=2)
mw = 11                                                         # the called motif is 11 bp, drawn to scale
axf.add_patch(plt.Rectangle((-mw / 2, RAILY + 0.17), mw, 0.22, color=GOLD, ec=INK, lw=0.9, zorder=5))
axf.annotate("AP-1 motif, 11 bp", xy=(mw / 2, RAILY + 0.28), xytext=(64, RAILY + 0.62), fontsize=8.2,
             color="#9a7d0a", fontweight="bold", arrowprops=dict(arrowstyle="->", color="#9a7d0a", lw=0.9))
axf.set_ylim(RAILY - 0.22 - max(len(label_rows), 1) * ROW - 0.08, RAILY + 0.86)
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
      f"{len(fp)} merged footprints of {fp.family.nunique()} families")
