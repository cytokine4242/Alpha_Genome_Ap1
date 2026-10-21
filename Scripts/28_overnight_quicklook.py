#!/usr/bin/env python3
"""Quick look at the overnight results: controls + width titration (CM), mutation arm (CM, FB), FB shuffle arm, neutral-context MPRA.
All effects = per-site mean over replicates of log2((edited+pc)/(native+pc)), pc = 5% of the median native value of that readout and trackset.
  python3 28_overnight_quicklook.py"""
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

warnings.filterwarnings("ignore")
R = Path(__file__).resolve().parents[1] / "Katana" / "results"
pd.set_option("display.width", 220, "display.max_columns", 30)
load = lambda p: pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{p}_chunk*.csv")))], ignore_index=True)


def lfc(d, ts, col):
    nat = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[(d.condition != "native") & (d.trackset == ts)].copy()
    e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "condition"]).l.mean().unstack()


def by_class(t, d, cols):
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    t = t.join(sd)
    return pd.concat({k: t[t.stable_dynamic == k][cols].median() for k in ("Stable", "Dynamic")}, axis=1)


# ---------- A. controls and titration (CM sites, heart ATAC)
co = load("controls_cm")
print("=== A. controls + width titration, CM sites, heart ATAC: median log2FC at each readout")
for rd in ("atac_motif", "atac_170", "atac_peak"):
    t = lfc(co, "CM", rd)
    conds = ["titrate@7", "titrate@11", "titrate@21", "titrate@41", "titrate@85", "titrate@170", "titrate@340", "titrate@680", "offset@-600", "offset@-300", "offset@+300", "offset@+600"]
    print(f"\nreadout = {rd}"); print(by_class(t, co, conds).round(2).T.to_string())
# is the centred 170 bp effect larger than a same-size region 300 bp away (paired, per site)?
t = lfc(co, "CM", "atac_motif")
off = t[["offset@-300", "offset@+300"]].mean(1)
cen = t["titrate@170"]
sd = co.drop_duplicates("site_id").set_index("site_id").stable_dynamic
for k in ("Stable", "Dynamic"):
    m = sd[sd == k].index
    print(f"\n{k}: motif-readout effect, centred 170 bp shuffle {cen[m].median():+.2f} vs offset +-300 bp {off[m].median():+.2f}; share of sites where centred is more negative: {(cen[m] < off[m]).mean():.0%}")

# ---------- B. mutation arm
print("\n=== B. mutation arm: dependency = -log2FC (mean over sites); columns = class / readout region")
for lab, f, ts in (("CM sites, heart ATAC", "mut_cm", "CM"), ("FB sites, Patski ATAC", "mut_fb", "FB")):
    d = load(f)
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    cols = {}
    for rd in ("atac_motif", "atac_170", "atac_peak"):
        t_ = -lfc(d, ts, rd)
        for k in ("Stable", "Dynamic"):
            cols[(k, rd.replace("atac_", ""))] = t_.loc[sd[sd == k].index].mean()
    print(f"\n{lab}"); print(pd.DataFrame(cols).round(2).to_string())
print("slides (CM, Stable, motif/170/peak readouts): random@motif 0.47/0.38/0.19 | dinuc@motif 0.29/0.31/0.15 | targeted core 1.18/0.99/0.61")

# ---------- C. FB shuffle arm
print("\n=== C. FB sites, shuffle arm: median matched effect (shuffle width -> readout of the same region), Patski ATAC and heart ATAC")
fb = load("shuffle_fb")
sdf = fb.drop_duplicates("site_id").set_index("site_id").stable_dynamic
for ts, lab in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
    res = {}
    for w, rd in (("motif", "atac_motif"), ("170", "atac_170"), ("peak", "atac_peak")):
        e = lfc(fb, ts, rd)[f"shuffle@{w}"]
        s, dy = e[sdf[sdf == "Stable"].index].dropna(), e[sdf[sdf == "Dynamic"].index].dropna()
        res[w] = {"Stable": s.median(), "Dynamic": dy.median(), "diff": dy.median() - s.median(), "MWU p": mannwhitneyu(dy, s).pvalue, "n_s": len(s), "n_d": len(dy)}
    print(f"\n{lab}"); print(pd.DataFrame(res).T.round(3).to_string())

# ---------- D. neutral-context MPRA
print("\n=== D. MPRA elements in a neutral context vs MPRA (n elements = 1653)")
mn = load_m = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / "mpra_neutral_chunk*.csv")))], ignore_index=True)
mt = pd.read_csv(R.parent / "data" / "MPRA_Results_with_TSS_distance.csv").drop_duplicates("Coords").set_index("Coords")
for ts, lab in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
    g = mn[mn.trackset == ts]
    nat = g[g.condition == "native"].groupby("mpra_id").atac_170.mean()
    pc = 0.05 * nat.median()
    res = {"native": nat}
    for c in ("motif_shuffle", "full_shuffle"):
        s = g[g.condition == c].groupby(["mpra_id", "background"]).atac_170.mean().groupby("mpra_id").mean()
        res[c] = np.log2((s + pc) / (nat.reindex(s.index) + pc))
    for name, x, y in (("native predicted ATAC vs log2 alpha", res["native"], mt.log2_alpha_Endo), ("motif-shuffle effect vs MPRA native - motif-shuffled", res["motif_shuffle"], mt.delta_vs_MotifShuff), ("full-shuffle effect vs MPRA native - fully shuffled", res["full_shuffle"], mt.delta_vs_FullShuff)):
        ids = [i for i in x.index if i in y.index]
        xx, yy = x.reindex(ids), y.reindex(ids)
        m = xx.notna() & yy.notna()
        r, p = spearmanr(xx[m], yy[m])
        print(f"  [{lab}] {name:55s} rho = {r:+.3f} (p = {p:.2g}, n = {int(m.sum())})")
