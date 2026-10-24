#!/usr/bin/env python3
"""Do Stable and Dynamic fibroblast AP-1 sites USE their partner motifs differently?

The earlier partner analysis compared all sites, so it mixed two things: Dynamic sites have fewer detectable footprints
(they are less accessible, and footprint calling needs signal) and they are less sensitive overall. This asks the
conditional question instead - given that a site HAS a footprint of family F, how much does that footprint matter? -
and compares classes at matched accessibility.

Per site and family: importance = mean over the tiles overlapping any footprint of F of (|log2 effect| - the site's own
noise floor, from tiles at |pos| >= 200 bp). Two readouts:
  absolute    the importance itself
  share       that importance divided by the site's total tile signal, so a site that is weaker everywhere is not
              counted as using its partners less
Each Dynamic site is paired with an unused Stable site of similar native accessibility (caliper 0.1 log10) and the two
are compared with a paired test; the unmatched comparison is reported alongside so the effect of matching is visible.
For every family the smallest difference detectable at 80% power is reported, so a null can be read as a null.

Inputs: Katana/results_ft/scan_fb_chunk*.csv and Analysis/partners_fb/tables/footprints_near_sites.csv (from
38_tfbs_partner_mapping.py --cell fb).   Outputs Analysis/partners_fb/tables/utilisation_*.csv
  python3 51_fb_motif_utilisation.py
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--min-sites", type=int, default=25, help="minimum sites per class carrying the family")
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
PART = ROOT / "Analysis" / "partners_fb" / "tables"
TILE, HALF = 10, 300


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ---------------------------------------------------------------- tile effects, noise floor removed
cols = ["site_id", "group", "stable_dynamic", "kind", "pos", "rep", "trackset", a.readout]
sc = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in sorted(glob.glob(str(RES / "scan_fb_chunk*.csv")))],
               ignore_index=True)
sc = sc[sc.trackset == "FB"]
nat = sc[sc.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t[a.readout] + pc) / (t.site_id.map(nat) + pc))
t = t.groupby(["site_id", "group", "stable_dynamic", "pos"], as_index=False).eff.mean()
floor = t[t.pos.abs() >= 200].assign(x=lambda z: z.eff.abs()).groupby("site_id").x.median().rename("floor")
t = t.join(floor, on="site_id")
t["signal"] = (t.eff.abs() - t.floor).clip(lower=0)
total = t.groupby("site_id").signal.sum().rename("total")
meta = t.drop_duplicates("site_id").set_index("site_id")[["group", "stable_dynamic"]]
meta["native"] = nat.reindex(meta.index)
meta["total"] = total.reindex(meta.index)
print(f"{len(meta)} sites: " + ", ".join(f"{k} {v}" for k, v in meta.stable_dynamic.value_counts().items()))

# ---------------------------------------------------------------- which tiles overlap which family
H = pd.read_csv(PART / "footprints_near_sites.csv")
H = H[~H.own_ap1]                                        # the site's own AP-1 motif is not a partner
tiles = np.arange(-HALF, HALF, TILE)
rows = []
for (sid, fam), g in H.groupby(["site_id", "family"], observed=True):
    hit = np.zeros(len(tiles), bool)
    for s_, e_ in zip(g.rel_start.values, g.rel_end.values):
        hit |= (tiles + TILE > s_) & (tiles < e_)
    if hit.any():
        rows.append(pd.DataFrame({"site_id": sid, "family": fam, "pos": tiles[hit]}))
M = pd.concat(rows, ignore_index=True)
u = M.merge(t[["site_id", "pos", "signal"]], on=["site_id", "pos"], how="left").dropna(subset=["signal"])
per = u.groupby(["site_id", "family"], observed=True).signal.mean().rename("importance").reset_index()
per = per.join(meta, on="site_id")
per["share"] = per.importance / per.total.replace(0, np.nan)
per.to_csv(PART / "utilisation_per_site.csv", index=False)

# ---------------------------------------------------------------- accessibility-matched pairing (once, reused for every family)
dyn = meta[meta.stable_dynamic == "Dynamic"].copy()
sta = meta[meta.stable_dynamic == "Stable"].copy()
dyn["l"], sta["l"] = np.log10(dyn.native.clip(lower=1e-3)), np.log10(sta.native.clip(lower=1e-3))
sta = sta.sort_values("l")
used, pair_of = set(), {}
for sid, r in dyn.sort_values("l").iterrows():
    free = sta[~sta.index.isin(used)]
    if free.empty:
        break
    j = (free.l - r.l).abs().idxmin()
    if abs(sta.loc[j, "l"] - r.l) <= 0.1:
        used.add(j)
        pair_of[sid] = j
print(f"accessibility-matched pairs: {len(pair_of)}  "
      f"(Dynamic {dyn.loc[list(pair_of), 'native'].median():.2f} vs Stable {sta.loc[list(pair_of.values()), 'native'].median():.2f})")

imp = {(r.site_id, r.family): (r.importance, r.share) for r in per.itertuples()}
fams = sorted(per.family.unique())
rows = []
for fam in fams:
    d_un = per[(per.family == fam) & (per.stable_dynamic == "Dynamic")]
    s_un = per[(per.family == fam) & (per.stable_dynamic == "Stable")]
    if min(len(d_un), len(s_un)) < a.min_sites:
        continue
    # matched: keep pairs where BOTH partners carry this family, so presence is controlled as well as accessibility
    pa, pb = [], []
    for ds, ss in pair_of.items():
        if (ds, fam) in imp and (ss, fam) in imp:
            pa.append(imp[(ds, fam)])
            pb.append(imp[(ss, fam)])
    row = {"family": fam, "n_dynamic_unmatched": len(d_un), "n_stable_unmatched": len(s_un),
           "importance_dynamic": d_un.importance.median(), "importance_stable": s_un.importance.median(),
           "p_importance_unmatched": mannwhitneyu(d_un.importance, s_un.importance).pvalue,
           "share_dynamic": d_un.share.median(), "share_stable": s_un.share.median(),
           "p_share_unmatched": mannwhitneyu(d_un.share.dropna(), s_un.share.dropna()).pvalue, "n_pairs": len(pa)}
    if len(pa) >= a.min_sites:
        A, B = np.array(pa), np.array(pb)
        for k, name in ((0, "importance"), (1, "share")):
            diff = A[:, k] - B[:, k]
            ok = ~np.isnan(diff)
            row[f"{name}_paired_difference"] = float(np.median(diff[ok]))
            row[f"p_{name}_matched"] = float(wilcoxon(diff[ok]).pvalue) if ok.sum() > 10 and np.any(diff[ok] != 0) else np.nan
            # smallest difference this many pairs could detect at 80% power, 5% two-sided
            row[f"{name}_min_detectable_80pct"] = float(2.8 * np.nanstd(diff[ok]) / np.sqrt(ok.sum()))
    rows.append(row)
R = pd.DataFrame(rows)
for c in ("p_importance_matched", "p_share_matched", "p_importance_unmatched", "p_share_unmatched"):
    if c in R:
        R[c.replace("p_", "fdr_")] = bh(R[c].fillna(1).values)
R = R.sort_values("importance_dynamic")
R.to_csv(PART / "utilisation_by_family.csv", index=False)

pd.set_option("display.width", 250)
show = ["family", "n_pairs", "importance_dynamic", "importance_stable", "importance_paired_difference",
        "p_importance_matched", "fdr_importance_matched", "importance_min_detectable_80pct"]
print("\n=== ABSOLUTE importance of each partner family, Dynamic vs accessibility-matched Stable")
print(R[[c for c in show if c in R]].round(4).to_string(index=False))
show2 = ["family", "n_pairs", "share_dynamic", "share_stable", "share_paired_difference", "p_share_matched",
         "fdr_share_matched", "share_min_detectable_80pct"]
print("\n=== SHARE of the site's own total dependence (utilisation, controls for overall sensitivity)")
print(R[[c for c in show2 if c in R]].round(4).to_string(index=False))
if "fdr_share_matched" in R:
    k = R[R.fdr_share_matched < 0.05]
    print(f"\nfamilies whose SHARE differs between classes at FDR<0.05: {len(k)} of {len(R)}"
          + ("" if k.empty else ": " + ", ".join(k.family)))
    k2 = R[R.fdr_importance_matched < 0.05]
    print(f"families whose ABSOLUTE importance differs: {len(k2)} of {len(R)}"
          + ("" if k2.empty else ": " + ", ".join(k2.family)))
print("\noutputs in", PART)
