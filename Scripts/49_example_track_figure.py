#!/usr/bin/env python3
"""A presentation figure: the sliding-window mutation scan across one Stable and one Dynamic AP-1 element, drawn like a genome browser track.
Each bar is one 10 bp window replaced by random bases (3 replicates, mean); the bar height is the predicted loss of accessibility over the
surrounding 170 bp. Beneath each track are the footprinted partner motifs. The pair is matched on predicted accessibility and peak width, and
was picked to show the trend the whole dataset carries: Stable elements concentrate their dependence on the AP-1 motif, Dynamic elements spread it out.
  python3 49_example_track_figure.py [--dynamic SITE_ID --stable SITE_ID]
Outputs Analysis/example_track/."""
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
OUT = ROOT / "Analysis" / "example_track"
OUT.mkdir(parents=True, exist_ok=True)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#ffffff", "#e4e3df", "#1F6FB2", "#C2551C", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 11, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 220})
ap = argparse.ArgumentParser()
ap.add_argument("--stable", default="chr11-50211071-50211082")       # Sqstm1 promoter, bound in all three conditions
ap.add_argument("--dynamic", default="chr11-80428854-80428865")      # Psmd11 promoter, AP-1 gained in both aged states
a = ap.parse_args()

# ---------------------------------------------------------------- data
sc = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[(sc.trackset == "CM") & sc.site_id.isin([a.stable, a.dynamic])]
nat = sc[sc.kind == "native"].set_index("site_id").atac_170
pc = 0.05 * pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "atac_170", "trackset"]) for f in sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))[:4]]).query("kind=='native' and trackset=='CM'").atac_170.median()
e = sc[sc.kind != "native"].copy()
e["eff"] = np.log2((e.atac_170 + pc) / (e.site_id.map(nat) + pc))
TILE = e[e.kind == "tile"].groupby(["site_id", "pos"], as_index=False).eff.mean()
TILE["centre"] = TILE.pos + 5
SF = pd.read_csv(ROOT / "Analysis/scan_features/tables/scan_features_per_site.csv").set_index("site_id")
M = pd.read_csv(ROOT / "Analysis/distance_rules/tables/merged_partner_instances.csv.gz")
SHOW = ["MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD", "NKX / homeobox cardiac", "TBX", "CTCF / cohesin", "SRF / MADS", "AP-1 (Fos/Jun/Atf/Batf)", "MYC / MAX", "RUNX", "ETS"]
FCOL = {"MEF2": "#c0392b", "GATA": "#8e44ad", "NFI (Nfia/b/c/x)": "#16a085", "TEAD": "#d68910", "NKX / homeobox cardiac": "#2980b9",
        "TBX": "#7f8c8d", "CTCF / cohesin": "#2c3e50", "SRF / MADS": "#e67e22", "AP-1 (Fos/Jun/Atf/Batf)": "#555555", "MYC / MAX": "#95a5a6", "RUNX": "#95a5a6", "ETS": "#95a5a6"}

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(16.4, 7.4))
gs = fig.add_gridspec(5, 1, height_ratios=[3.0, 0.9, 0.55, 3.0, 0.9], hspace=0.16)
GE = pd.read_csv(ROOT / "Data/derived/cm_sites_nearest_gene_mast.csv").set_index("site_id")
SL = pd.read_csv(ROOT / "Data/derived/cm_sites_nearest_gene_mast.csv").set_index("site_id")
CAT = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
PATTERN = {"YS+MAS+MAE": "AP-1 bound in all three conditions", "MAS+MAE": "no AP-1 in young, bound in both aged states",
           "MAS only": "no AP-1 in young, bound in aged sedentary only", "MAE only": "AP-1 only after exercise",
           "YS only": "AP-1 in young only, lost with age", "YS+MAS": "AP-1 lost with exercise", "YS+MAE": "AP-1 lost in aged sedentary"}
panels = [(a.stable, "Stable AP-1 site", BLUE, 0), (a.dynamic, "Dynamic AP-1 site \u2014 gains binding with age", ORANGE, 3)]
ymin = min(TILE.eff.min(), -0.2) * 1.15
for sid, lab, col, row in panels:
    ax = fig.add_subplot(gs[row])
    w = TILE[TILE.site_id == sid].sort_values("centre")
    f = SF.loc[sid]
    heavy = w.eff < 0.35 * w.eff.min()
    ax.bar(w.centre[~heavy], w.eff[~heavy], width=9, color=GREY, alpha=0.55, zorder=3)
    ax.bar(w.centre[heavy], w.eff[heavy], width=9, color=col, zorder=3)
    ax.axvspan(-5.5, 5.5, color="#f2c14e", alpha=0.45, lw=0, zorder=1)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    ax.set_xlim(-305, 305)
    ax.set_ylim(ymin, max(0.12, w.eff.max() * 1.55))
    ax.set_ylabel("Loss of predicted\naccessibility (log2FC)", fontsize=10.5)
    chrom, s_, e_ = sid.split("-")
    g = GE.loc[sid]
    gname, gdist, grna_a, grna_e = g.gene, g.dist_gene / 1000, g.rna_age, g.rna_ex
    ax.set_title(f"{lab}   ·   {chrom}:{int(s_):,}   ·   {gname} promoter, {gdist * 1000:.0f} bp from the TSS",
                 fontsize=12.5, loc="left", color=col, fontweight="bold", pad=6)
    fp = "AP-1 absent in young, gained in both aged states" if sid == a.dynamic else "AP-1 present in all three conditions"
    ax.text(0.004, 0.96, f"{fp}    \u00b7    {gname} RNA {grna_a:+.2f} with age, {grna_e:+.2f} with exercise",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=col, fontweight="bold")
    ax.text(0.996, 0.035, f"AP-1 motif window {f.own_motif_170:.2f}   \u00b7   strongest other {f.best_partner_170:.2f}   \u00b7   "
                          f"{f.share_within_85_170:.0%} of the loss within \u00b185 bp",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5, color=MUTED)
    ax.annotate("AP-1 motif", xy=(0, ymin * 0.14), xytext=(78, ymin * 0.42), fontsize=10, color="#9a7d0a",
                arrowprops=dict(arrowstyle="->", color="#9a7d0a", lw=1.3))
    ax.tick_params(labelbottom=False)

    axf = fig.add_subplot(gs[row + 1], sharex=ax)
    w2 = M[(M.site_id == sid) & (M.family.isin(SHOW)) & (M.rs.abs() < 320)]
    lanes = []
    for r in w2.sort_values("rs").itertuples():
        k = next((i for i, end in enumerate(lanes) if r.rs > end + 48), len(lanes))
        if k == len(lanes):
            lanes.append(r.re)
        else:
            lanes[k] = r.re
        c = FCOL.get(r.family, GREY)
        key = r.family in ("MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD", "SRF / MADS", "NKX / homeobox cardiac")
        nm0 = r.family.split(" (")[0].split(" /")[0]
        wmin = 7.5 * len(nm0) + 10                       # wide enough for the label, so chips never clip their own text
        cx = (r.rs + r.re) / 2
        wid = max(r.re - r.rs, wmin)
        axf.add_patch(plt.Rectangle((cx - wid / 2, -k - 0.42), wid, 0.62, color=c, alpha=0.95 if key else 0.45,
                                    ec=INK if key else "none", lw=1.1 if key else 0))
        axf.text(cx, -k - 0.11, nm0, ha="center", va="center", fontsize=8.2 if key else 7.2,
                 color="white" if key else INK, fontweight="bold" if key else "normal")
    axf.set_ylim(-max(len(lanes), 1) - 0.1, 0.35)
    axf.set_xlim(-305, 305)
    axf.set_yticks([])
    axf.grid(False)
    axf.set_ylabel("footprinted\npartners", fontsize=9, color=MUTED)
    for sp in ("left", "bottom"):
        axf.spines[sp].set_visible(False)
    if row == 3:
        axf.set_xlabel("Position of the replaced 10 bp window, relative to the AP-1 motif (bp)", fontsize=11.5)
        axf.tick_params(labelbottom=True, labelsize=10)
    else:
        axf.tick_params(labelbottom=False)

fig.suptitle("Sliding-window mutagenesis across two AP-1 elements matched for chromatin context\n"
             "each bar is one 10 bp window replaced by random bases; height is the predicted loss of accessibility over the surrounding 170 bp",
             fontsize=12.5, x=0.012, ha="left", y=0.998)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT / "example_scan_track.png", bbox_inches="tight")
fig.savefig(OUT / "example_scan_track.pdf", bbox_inches="tight")
plt.close(fig)

info = SF.loc[[a.stable, a.dynamic], ["own_motif_170", "best_partner_170", "share_within_85_170", "spread_170", "gini_170", "total_dependence_170", "n_tiles_matter_170"]]
info.to_csv(OUT / "example_sites.csv")
pd.set_option("display.width", 220)
print(info.round(3).to_string())
print(f"\nwrote {OUT / 'example_scan_track.png'} and .pdf")
