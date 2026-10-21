#!/usr/bin/env python3
"""Summary of the pilot results (Katana/results/pilot_*): sanity checks, native vs measured ATAC, shuffle / mutation effects
against the slide numbers, ATAC-RNA-CAGE coupling, and the fibroblast MPRA comparison. Prints tables; no plots.

  python3 20_pilot_summary.py [results_dir] [prefix]      (default: Katana/results pilot_)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Katana" / "results"
PRE = sys.argv[2] if len(sys.argv) > 2 else "pilot_"
pd.set_option("display.width", 200, "display.max_columns", 30)


def log2fc(df, col, pc_frac=0.05):
    """Per site x condition mean over replicates of log2((edit + pc) / (native + pc)); pc = 5% of the median native value."""
    nat = df[df.condition == "native"].set_index(["site_id", "trackset"])[col]
    pc = pc_frac * nat.median()
    e = df[df.condition != "native"].copy()
    e["nat"] = e.set_index(["site_id", "trackset"]).index.map(nat)
    e["lfc"] = np.log2((e[col] + pc) / (e["nat"] + pc))
    return e.groupby(["site_id", "trackset", "condition", "stable_dynamic"], as_index=False)["lfc"].mean()


def rho(a, b):
    m = a.notna() & b.notna()
    if m.sum() < 5:
        return "n<5"
    r, p = spearmanr(a[m], b[m])
    return f"rho={r:+.2f} (p={p:.2g}, n={m.sum()})"


cm = pd.read_csv(RES / f"{PRE}cm_chunk0.csv")
fb = pd.read_csv(RES / f"{PRE}fb_chunk0.csv")
mp = pd.read_csv(RES / f"{PRE}mpra_chunk0.csv")
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
mpra = pd.read_csv(ROOT / "Katana" / "data" / "MPRA_Results_with_TSS_distance.csv")

print("=== 1. native ATAC prediction vs measured snATAC (CPM/kb of the peak containing the motif)")
for name, df, tab, ts, col in (("CM sites, heart ATAC", cm, "cm", "CM", "cm_cpm_per_kb"), ("FB sites, Patski ATAC", fb, "fb", "FB", "fb_cpm_per_kb"),
                               ("FB sites, heart ATAC", fb, "fb", "CM", "fb_cpm_per_kb"), ("CM sites, Patski ATAC", cm, "cm", "FB", "cm_cpm_per_kb")):
    n = df[(df.condition == "native") & (df.trackset == ts)].merge(meas[meas.site_table == tab][["site_id", col]], on="site_id")
    print(f"{name:24s} atac_peak {rho(n.atac_peak, n[col])} | atac_motif {rho(n.atac_motif, n[col])} | atac_170 {rho(n.atac_170, n[col])}")

READ = ["atac_motif", "atac_170", "atac_peak"]
for label, df, ts in (("CM sites, CM tracks", cm, "CM"), ("FB sites, FB tracks", fb, "FB")):
    print(f"\n=== 2. shuffle effects, {label}: median log2(shuffled/native) per site (rows = region shuffled, cols = ATAC readout)")
    parts = {c: log2fc(df[df.trackset == ts], c) for c in READ}
    for sd in ("Stable", "Dynamic"):
        t = pd.DataFrame({c: p[(p.stable_dynamic == sd) & p.condition.str.startswith("shuffle")].groupby("condition").lfc.median()
                          for c, p in parts.items()}).reindex(["shuffle@motif", "shuffle@170", "shuffle@peak"])
        print(f"{sd} (n sites = {df[df.stable_dynamic == sd].site_id.nunique()})\n{t.round(2).to_string()}")
print("\nslides (CM, Stable, fine-tuned model): motif shuffle -0.12/-0.26/-0.11, 170 bp shuffle -2.25/-1.84/-0.75, peak shuffle -2.76/-2.53/-1.53 (motif/170/peak readouts)")

print("\n=== 3. mutation dependency = -log2FC (mean over sites), CM sites, readouts motif / 170 / peak")
parts = {c: log2fc(cm[cm.trackset == "CM"], c) for c in READ}
t = pd.DataFrame({c: -p[~p.condition.str.startswith("shuffle")].groupby("condition").lfc.mean() for c, p in parts.items()})
print(t.round(2).to_string())
print("slides (CM, Stable, ATAC motif/170/peak): random@motif 0.47/0.38/0.19, dinuc@motif 0.29/0.31/0.15, targeted core 1.18/0.99/0.61")

print("\n=== 4. replicate noise: SD across the 10 shuffle replicates vs the effect (CM tracks, atac_peak, median over sites)")
s = cm[(cm.trackset == "CM") & cm.condition.str.startswith("shuffle")]
nat = cm[(cm.condition == "native") & (cm.trackset == "CM")].set_index("site_id").atac_peak
pc = 0.05 * nat.median()
s = s.assign(l=np.log2((s.atac_peak + pc) / (s.site_id.map(nat) + pc)))
g = s.groupby(["condition", "site_id"]).l.agg(["mean", "std"]).groupby("condition").median()
print(g.round(3).to_string())

print("\n=== 5. native coupling across sites (Spearman), CM tracks | FB tracks")
for df, ts in ((cm, "CM"), (fb, "FB")):
    n = df[(df.condition == "native") & (df.trackset == ts)]
    print(ts, "| ATAC peak vs RNA prox", rho(n.atac_peak, n.rna_prox2kb), "| ATAC peak vs CAGE prox", rho(n.atac_peak, n.cage_prox2kb),
          "| RNA prox vs CAGE prox", rho(n.rna_prox2kb, n.cage_prox2kb))
    print("   RNA mean50 vs CAGE mean50", rho(n.rna_mean50kb, n.cage_mean50kb), "| ATAC peak vs CAGE mean50", rho(n.atac_peak, n.cage_mean50kb))

print("\n=== 6. peak/170 bp shuffle effect on RNA / CAGE (median log2FC over sites, same-cell-type tracks)")
for label, df, ts in (("CM", cm, "CM"), ("FB", fb, "FB")):
    d = df[df.trackset == ts]
    out = {}
    for c in ("atac_peak", "rna_prox2kb", "rna_mean50kb", "cage_prox2kb", "cage_mean50kb"):
        p = log2fc(d, c)
        out[c] = p[p.condition.isin(["shuffle@motif", "shuffle@170", "shuffle@peak"])].groupby("condition").lfc.median()
    print(label); print(pd.DataFrame(out).round(2).to_string())
    a = log2fc(d, "atac_peak"); a = a[a.condition == "shuffle@peak"].set_index("site_id").lfc
    for c in ("rna_prox2kb", "cage_prox2kb"):
        p = log2fc(d, c); p = p[p.condition == "shuffle@peak"].set_index("site_id").lfc
        print(f"   peak-shuffle effect, ATAC vs {c}: {rho(a, p)}")

print("\n=== 7. fibroblast MPRA: predicted vs measured (per element)")
nat = mp[mp.condition == "native"].groupby(["mpra_id", "trackset"]).mean(numeric_only=True)
rep = mp[mp.condition != "native"].groupby(["mpra_id", "trackset", "condition"]).mean(numeric_only=True).reset_index()
# 359 elements are listed under two MPRA groups with identical values: collapse to unique elements
groups = mpra.groupby("Coords").Group.agg(lambda s: "+".join(sorted(set(s))))
m = mpra.drop_duplicates("Coords").set_index("Coords")
m["Group"] = groups
mp = mp.drop_duplicates(["mpra_id", "trackset", "condition", "rep"])
for ts in ("FB", "CM"):
    n = nat.xs(ts, level="trackset").join(m[["log2_alpha_Endo", "delta_vs_MotifShuff", "delta_vs_FullShuff", "dist_to_TSS", "Group"]])
    print(f"-- {ts} tracks, n elements {len(n)}")
    for c in ("atac_170", "cage_prox2kb", "rna_prox2kb"):
        print(f"   native {c:13s} vs log2_alpha_Endo: {rho(n[c], n.log2_alpha_Endo)}")
    for cond, mc in (("motif_shuffle", "delta_vs_MotifShuff"), ("full_shuffle", "delta_vs_FullShuff")):
        r = rep[(rep.trackset == ts) & (rep.condition == cond)].set_index("mpra_id")
        for c in ("atac_170", "cage_prox2kb"):
            pc = 0.05 * n[c].median()
            l = np.log2((r[c] + pc) / (n[c].reindex(r.index) + pc))
            print(f"   predicted log2FC {cond:13s} {c:13s} vs MPRA {mc}: {rho(l, m[mc].reindex(r.index))} | median pred {l.median():+.2f}")
print("\nMPRA motif source:", mp.drop_duplicates("mpra_id").motif_source.value_counts().to_dict())
