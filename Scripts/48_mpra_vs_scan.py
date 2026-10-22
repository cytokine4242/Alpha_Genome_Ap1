#!/usr/bin/env python3
"""Does the fibroblast MPRA agree with what the mutation scan says about how much the AP-1 motif is worth?
The comparable quantity is the motif's share of the LOCAL sequence dependence, which both experiments measure in their own way:
  MPRA (wet lab, episomal 170 bp reporter, fibroblast):  (Endo - MotifShuff) / (Endo - FullShuff)
  AlphaGenome (in silico, genomic context, cardiomyocyte): (replace 170 bp - replace surroundings with the motif kept) / replace 170 bp
Both are ratios within the same element, so cell type, assay scale and dynamic range cancel; only the apportioning is compared.
Also tests the class direction, where the two datasets disagree.
Outputs Analysis/mpra_scan/{tables,figures}.   python3 48_mpra_vs_scan.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "mpra_scan"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(48)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})


def boot_med(v, n=1000):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5]) if len(v) > 5 else (np.nan, np.nan)


# ---------------------------------------------------------------- AlphaGenome side
FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
den_ag = FD["replace 170 bp, read over 170 bp"]
num_ag = den_ag - FD["replace surroundings (motif kept), read over 170 bp"]
AG = pd.DataFrame({"share": (num_ag / den_ag).where(den_ag < -0.5), "whole": den_ag, "motif_only": num_ag}).join(sl[["stable_dynamic"]])

# ---------------------------------------------------------------- MPRA side
m = pd.read_csv(ROOT / "Data/MPRA_Results_with_TSS_distance.csv")
m = m[m.replicated == True].copy()
m["share"] = (m.delta_vs_MotifShuff / m.delta_vs_FullShuff).replace([np.inf, -np.inf], np.nan)
MP = m[m.delta_vs_FullShuff > 0.1].copy()                     # elements whose activity actually depends on their sequence

rows = [{"dataset": "AlphaGenome (CM, genomic context)", "n": int(AG.share.notna().sum()), "median_share": AG.share.median(),
         "ci_lo": boot_med(AG.share.values)[0], "ci_hi": boot_med(AG.share.values)[1],
         "whole_element_effect": AG.whole.median(), "motif_only_effect": AG.motif_only.median()},
        {"dataset": "MPRA (fibroblast, episomal)", "n": int(MP.share.notna().sum()), "median_share": MP.share.median(),
         "ci_lo": boot_med(MP.share.values)[0], "ci_hi": boot_med(MP.share.values)[1],
         "whole_element_effect": MP.delta_vs_FullShuff.median(), "motif_only_effect": MP.delta_vs_MotifShuff.median()}]
A = pd.DataFrame(rows)
A.to_csv(OUT / "tables" / "motif_share_both_assays.csv", index=False)

# ---------------------------------------------------------------- class direction, where they disagree
rep = m
st, dy = rep[rep.Group == "Stable"], rep[rep.Group != "Stable"]
rows = []
for c, lab in (("delta_vs_FullShuff", "whole element shuffled"), ("delta_vs_MotifShuff", "AP-1 motif shuffled"), ("log2_alpha_Endo", "baseline activity")):
    u = mannwhitneyu(st[c].dropna(), dy[c].dropna())
    rows.append({"assay": "MPRA (fibroblast)", "measure": lab, "stable": st[c].median(), "dynamic": dy[c].median(),
                 "difference": dy[c].median() - st[c].median(), "p": u.pvalue, "n_stable": len(st), "n_dynamic": len(dy)})
agd = AG.dropna(subset=["stable_dynamic"])
for c, lab in (("whole", "whole element shuffled"), ("motif_only", "AP-1 motif only")):
    x, y = agd[agd.stable_dynamic == "Stable"][c].dropna(), agd[agd.stable_dynamic == "Dynamic"][c].dropna()
    rows.append({"assay": "AlphaGenome (CM)", "measure": lab, "stable": x.median(), "dynamic": y.median(),
                 "difference": y.median() - x.median(), "p": mannwhitneyu(x, y).pvalue, "n_stable": len(x), "n_dynamic": len(y)})
B = pd.DataFrame(rows)
B["kruskal_5_groups_mpra"] = [kruskal(*[g.delta_vs_FullShuff.dropna() for _, g in rep.groupby("Group")]).pvalue if r.assay.startswith("MPRA") and r.measure == "whole element shuffled" else np.nan for r in B.itertuples()]
B.to_csv(OUT / "tables" / "class_direction_both_assays.csv", index=False)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2), gridspec_kw={"width_ratios": [1, 1.15, 1.1]})
ax = axes[0]
for k, (v, lab, c) in enumerate(((AG.share.dropna().values, "AlphaGenome\n(cardiomyocyte,\ngenomic context)", BLUE),
                                 (MP.share.dropna().values, "MPRA\n(fibroblast,\nepisomal reporter)", ORANGE))):
    v = v[(v > -0.5) & (v < 1.5)]
    bp = ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, patch_artist=True,
                    boxprops={"facecolor": c, "alpha": 0.55}, medianprops={"color": INK, "lw": 1.8})
    ax.text(k, np.median(v) + 0.06, f"{np.median(v):.0%}\nn = {len(v):,}", ha="center", fontsize=8.5)
ax.set_xticks([0, 1], ["AlphaGenome\n(CM, in genome)", "MPRA\n(fibroblast, episomal)"], fontsize=8)
ax.axhline(0, color=MUTED, lw=0.9)
ax.set_ylabel("Share of the element's sequence dependence\nattributable to the AP-1 motif")
ax.set_title("Two independent assays, same answer", fontsize=10, loc="left")

ax = axes[1]
w = B[B.assay.str.startswith("MPRA")]
ys = np.arange(len(w))
for y_, r in zip(ys, w.itertuples()):
    ax.plot([r.stable, r.dynamic], [y_, y_], color=GREY, lw=1.2, zorder=2)
    ax.scatter(r.stable, y_, s=48, color=BLUE, zorder=3, label="Stable" if y_ == 0 else None)
    ax.scatter(r.dynamic, y_, s=48, color=ORANGE, zorder=3, label="dynamic groups" if y_ == 0 else None)
    if r.p < 0.05:
        ax.text(max(r.stable, r.dynamic) + 0.012, y_, f"p = {r.p:.3g}", va="center", fontsize=7.5)
ax.set_yticks(ys, w.measure, fontsize=8)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_xlabel("MPRA effect, log2 (Endo − shuffled)")
ax.set_title("MPRA: dynamic elements depend MORE on sequence", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.set_xlim(-0.02, 0.32)

ax = axes[2]
g = MP.groupby("Group").share.median().sort_values()
lo = [boot_med(MP[MP.Group == k].share.values)[0] for k in g.index]
hi = [boot_med(MP[MP.Group == k].share.values)[1] for k in g.index]
ax.errorbar(g.values, np.arange(len(g)), xerr=[g.values - np.array(lo), np.array(hi) - g.values], fmt="o", color=ORANGE, ms=6, lw=1.4, capsize=0)
ax.axvline(AG.share.median(), color=BLUE, lw=1.6, ls="--", label=f"AlphaGenome, all CM sites ({AG.share.median():.0%})")
ax.set_yticks(range(len(g)), [f"{k} (n = {int((MP.Group == k).sum())})" for k in g.index], fontsize=8)
ax.set_xlabel("MPRA: motif's share of the element's dependence")
ax.set_title("By MPRA group, against the model's estimate", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_xlim(0.0, 0.22)
fig.tight_layout()
fig.savefig(OUT / "figures" / "m1_mpra_vs_scan.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 240)
print("== how much of the local dependence is the AP-1 motif worth?\n", A.round(3).to_string(index=False))
print("\n== class direction in each assay\n", B.round(4).to_string(index=False))
print("\n== MPRA motif share by group\n", MP.groupby("Group").share.agg(["size", "median"]).round(3).to_string())
