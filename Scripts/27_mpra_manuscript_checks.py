#!/usr/bin/env python3
"""Re-derives the manuscript's MPRA statements from Table S27 (Katana/data/MPRA_Results_with_TSS_distance.csv) and compares them with AlphaGenome's
native ATAC-head predictions on the same elements (Katana/results/mpra_chunk*.csv).   python3 27_mpra_manuscript_checks.py > Analysis/mpra_manuscript_checks.txt
core loss = log2_alpha_Endo - log2_alpha_MotifShuff ; flank loss = log2_alpha_MotifShuff - log2_alpha_FullShuff (activity lost by shuffling the flanks)."""
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
mt = pd.read_csv(ROOT / "Katana" / "data" / "MPRA_Results_with_TSS_distance.csv")
mt["core_loss"] = mt.log2_alpha_Endo - mt.log2_alpha_MotifShuff
mt["flank_loss"] = mt.log2_alpha_MotifShuff - mt.log2_alpha_FullShuff
mt["dyn"] = np.where(mt.Group == "Stable", "Stable", "Dynamic")


def cohen(a, b):
    s = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / s


print(f"rows {len(mt)}, unique elements {mt.Coords.nunique()}, Stable rows {(mt.Group == 'Stable').sum()}")
print("\n(a) MPRA endogenous activity (log2 alpha) by class vs Stable")
st = mt[mt.Group == "Stable"].log2_alpha_Endo
for g in ["Age Up", "Age Down", "Ex Up", "Ex Down"]:
    x = mt[mt.Group == g].log2_alpha_Endo
    print(f"  {g:9s} n={len(x):4d} median {x.median():+.3f} vs Stable {st.median():+.3f} | MWU p={mannwhitneyu(x, st).pvalue:.3g}")
print("\n(b) core vs flank contribution (manuscript: flank loss dynamic 0.128 vs stable 0.016, d=0.144, FDR=0.034)")
for lab, sel in (("all dynamic rows", mt.dyn == "Dynamic"), ("age-gained (Age Up)", mt.Group == "Age Up")):
    for nm, col in (("flank", "flank_loss"), ("core", "core_loss")):
        a, b = mt[sel][col], mt[mt.Group == "Stable"][col]
        print(f"  {nm:5s} loss {lab:20s} mean {a.mean():+.3f} vs Stable {b.mean():+.3f} | d={cohen(a, b):.3f} | MWU p={mannwhitneyu(a, b).pvalue:.3g}")
print(f"  overall mean core loss {mt.core_loss.mean():.3f}, flank loss {mt.flank_loss.mean():.3f}")
mp = pd.concat([pd.read_csv(f) for f in glob.glob(str(ROOT / "Katana" / "results" / "mpra_chunk*.csv"))]).drop_duplicates(["mpra_id", "trackset", "condition", "rep"])
print("\n(c) AlphaGenome-predicted (native ATAC head, genomic context, 170 bp readout), same elements")
for ts, lab in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
    g = mp[mp.trackset == ts]
    nat = g[g.condition == "native"].groupby("mpra_id").atac_170.mean()
    pc = 0.05 * nat.median()
    res = {}
    for c in ("motif_shuffle", "full_shuffle"):
        s = g[g.condition == c].groupby("mpra_id").atac_170.mean()
        res[c] = np.log2((s + pc) / (nat.reindex(s.index) + pc))
    P = pd.DataFrame({"native": nat, "core_pred": -res["motif_shuffle"], "full_pred": -res["full_shuffle"]})
    P["flank_pred"] = P.full_pred - P.core_pred
    M = mt.merge(P, left_on="Coords", right_index=True)
    print(f"  [{lab}] n rows {len(M)}")
    for c_, nm in (("native", "predicted native ATAC"), ("core_pred", "predicted core loss"), ("flank_pred", "predicted flank loss")):
        a, b = M[M.dyn == "Dynamic"][c_].dropna(), M[M.dyn == "Stable"][c_].dropna()
        print(f"     {nm:24s} Dynamic mean {a.mean():+.3f} vs Stable {b.mean():+.3f} | d={cohen(a, b):+.3f} | MWU p={mannwhitneyu(a, b).pvalue:.3g}")
    print("     Spearman across elements: predicted vs MPRA flank loss rho=%.3f | core loss rho=%.3f" % (spearmanr(M.flank_pred, M.flank_loss, nan_policy="omit")[0], spearmanr(M.core_pred, M.core_loss, nan_policy="omit")[0]))
    print("     Spearman(predicted native ATAC, log2_alpha_Endo) = %.3f" % spearmanr(M.native, M.log2_alpha_Endo, nan_policy="omit")[0])
print("\n(d) GC / TSS distance in the MPRA table: Stable vs Dynamic")
for c in ("GC_content", "dist_to_TSS"):
    a, b = mt[mt.dyn == "Dynamic"][c], mt[mt.dyn == "Stable"][c]
    print(f"  {c:12s} Dynamic median {a.median():.3f} vs Stable {b.median():.3f} | MWU p={mannwhitneyu(a, b).pvalue:.3g}")
print("  Spearman(GC content, log2_alpha_Endo) =", round(spearmanr(mt.GC_content, mt.log2_alpha_Endo)[0], 3))
