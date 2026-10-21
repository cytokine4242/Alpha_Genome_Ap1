#!/usr/bin/env python3
"""Insertion experiments (Katana/results/insertion_*.csv): peak swap (swapped ATAC vs donor and recipient native ATAC) and AP-1 motif
insertion into dinucleotide-shuffled backgrounds (log2FC vs shuffled background, vs shuffled-motif control, vs original background).
Tables and figures go to Analysis/insertion/.   python3 23_insertion_analysis.py"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "insertion"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
INK, MUTED, SURF, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
R = ROOT / "Katana" / "results"

ps = pd.read_csv(R / "insertion_peakswap.csv")
print(f"=== peak swap: {len(ps)} pairs, {ps.recipient_id.nunique()} recipients x {ps.donor_id.nunique()} donors")
rows = []
for name, x in (("donor native (sequence)", ps.donor_native_atac), ("recipient native (location)", ps.recipient_native_atac)):
    r, p = spearmanr(ps.swapped_atac, x)
    rows.append({"compared_with": name, "rho": r, "p": p, "n": len(ps)})
    for sd in ("Stable", "Dynamic"):
        g = ps[ps.donor_stable_dynamic == sd]
        r2, p2 = spearmanr(g.swapped_atac, g.donor_native_atac if "donor" in name else g.recipient_native_atac)
        rows.append({"compared_with": f"{name}, {sd} donors", "rho": r2, "p": p2, "n": len(g)})
t = pd.DataFrame(rows)
t.to_csv(OUT / "tables" / "peak_swap.csv", index=False)
print(t.round(3).to_string(index=False))
# within-recipient: does the donor's own accessibility still rank the swapped values?
w = ps.groupby("recipient_id").apply(lambda g: spearmanr(g.swapped_atac, g.donor_native_atac)[0])
print("within each recipient, swapped vs donor native rho:", w.round(2).to_dict())

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
for ax, col, xl in ((axes[0], "donor_native_atac", "Donor native ATAC (sequence)"), (axes[1], "recipient_native_atac", "Recipient native ATAC (location)")):
    for sd, c in (("Stable", BLUE), ("Dynamic", ORANGE)):
        g = ps[ps.donor_stable_dynamic == sd]
        ax.scatter(g[col], g.swapped_atac, s=14, alpha=0.6, color=c, linewidths=0, label=sd)
    r = spearmanr(ps.swapped_atac, ps[col])[0]
    ax.set_xlabel(xl)
    ax.set_title(f"rho = {r:.2f}", fontsize=10, loc="left")
axes[0].set_ylabel("Swapped ATAC (inserted peak)")
axes[0].legend(frameon=False, title="Donor")
fig.tight_layout()
fig.savefig(OUT / "figures" / "fig_peak_swap.png", bbox_inches="tight")
plt.close(fig)

mo = pd.read_csv(R / "insertion_motif.csv")
print("\n=== motif insertion: log2FC of ATAC at the insertion site (centre 7 bp / 128 bp / 170 bp), 50 backgrounds")
rows = []
fig, axes = plt.subplots(1, 3, figsize=(10, 3.4), sharey=True)
for ts in ("CM", "FB"):
    m = mo[mo.trackset == ts]
    for rd in ("atac_motif", "atac_peak", "atac_170"):
        p = m.pivot(index="background_id", columns="variant", values=rd)
        pc = 0.05 * p.median().median()
        lf = lambda a, b: np.log2((p[a] + pc) / (p[b] + pc))
        for label, a, b in (("+motif vs shuffled background", "shuffled_plus_motif", "background_shuffled"),
                            ("+motif vs +shuffled motif (control)", "shuffled_plus_motif", "shuffled_plus_shuffled_motif"),
                            ("+motif vs original background", "shuffled_plus_motif", "background_native")):
            x = lf(a, b)
            rows.append({"trackset": ts, "readout": rd, "comparison": label, "mean_log2fc": x.mean(), "median_log2fc": x.median(), "frac_positive": (x > 0).mean(),
                         "wilcoxon_p": wilcoxon(x).pvalue, "n": len(x)})
t = pd.DataFrame(rows)
t.to_csv(OUT / "tables" / "motif_insertion.csv", index=False)
print(t[t.trackset == "CM"].round(3).to_string(index=False))
m = mo[mo.trackset == "CM"]
for ax, (rd, ttl) in zip(axes, (("atac_motif", "Centre 7 bp"), ("atac_peak", "Centre 128 bp"), ("atac_170", "Centre 170 bp"))):
    p = m.pivot(index="background_id", columns="variant", values=rd)
    pc = 0.05 * p.median().median()
    for k, (a, b, c, lab) in enumerate((("shuffled_plus_motif", "background_shuffled", BLUE, "+AP-1 motif vs shuffled background"),
                                        ("shuffled_plus_motif", "shuffled_plus_shuffled_motif", ORANGE, "+AP-1 motif vs shuffled-motif control"))):
        x = np.log2((p[a] + pc) / (p[b] + pc))
        ax.scatter(np.full(len(x), k) + np.random.default_rng(1).uniform(-0.12, 0.12, len(x)), x, s=10, alpha=0.6, color=c, linewidths=0)
        ax.hlines(x.median(), k - 0.25, k + 0.25, color=INK, lw=1.5)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks([0, 1], ["vs shuffled\nbackground", "vs shuffled-\nmotif control"])
    ax.set_title(ttl, fontsize=10, loc="left")
axes[0].set_ylabel("log2FC ATAC (heart track), n = 50 backgrounds")
fig.tight_layout()
fig.savefig(OUT / "figures" / "fig_motif_insertion.png", bbox_inches="tight")
plt.close(fig)
print("\nbackground ATAC (CM heart track, mean, 50 backgrounds):", m.groupby("variant").atac_peak.mean().round(4).to_dict())
