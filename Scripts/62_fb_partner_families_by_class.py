#!/usr/bin/env python3
"""Besides AP-1, which motif families matter most at each class of fibroblast AP-1 site?

A note on where this comes from, because it is easy to expect it from the wrong experiment. The tile and sliding
scans are POSITIONAL: they say which 10 bp of sequence the prediction depends on, not what is written there. They
cannot name a motif. The partner shuffle can, because it shuffles each TOBIAS-called footprint individually and
scores it against a width- and distance-matched footprint-free window in the same site, so every number is already
net of "disturbing that many bases at that distance".

Per class and per family:
  net      median of (partner shuffled - its matched control), log2. More negative = the site depends on it more.
  share    the fraction of that class's sites that carry a called footprint of the family at all
Families are ranked within each class, AP-1 itself excluded. Each Dynamic class is tested against Stable on the
same family, BH-adjusted across the families tested within that class.

Detection is not independent of accessibility - a weaker site yields fewer called footprints - so `share` is
reported alongside, and a family present in fewer than --min-sites sites of a class is not ranked for it.

Outputs Analysis/fb_partner_shuffle/tables/family_rank_by_class.csv.
  python3 62_fb_partner_families_by_class.py [--min-sites 25]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--min-sites", type=int, default=25)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=62)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
TAB = ROOT / "Analysis" / "fb_partner_shuffle" / "tables"
CLS5 = ["Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down"]

F = pd.read_csv(TAB / "net_effect_per_site_family.csv")
S = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                usecols=["site_id", "stable_dynamic", "age_class", "ex_class"]).drop_duplicates("site_id")
S["cls"] = np.where(S.stable_dynamic == "Stable", "Stable",
            np.where(S.age_class == "Age Up", "Age Up",
            np.where(S.age_class == "Age Down", "Age Down",
            np.where(S.ex_class == "Exercise Up", "Exercise Up",
            np.where(S.ex_class == "Exercise Down", "Exercise Down", None)))))
d = F.merge(S[["site_id", "cls"]], on="site_id", how="inner")
d = d[d.cls.notna()]
d = d[~d.family.str.startswith("AP-1")]            # the site's own family is the subject of its own slide
n_cls = S[S.cls.notna()].groupby("cls").size()
print(f"{d.site_id.nunique()} sites, {d.family.nunique()} partner families, "
      f"{len(d)} site x family observations")

rows = []
for cls in CLS5:
    sub = d[d.cls == cls]
    for fam, g in sub.groupby("family"):
        if len(g) < a.min_sites:
            continue
        v = g.net.to_numpy()
        b = np.array([np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(a.boot)])
        r = {"class": cls, "family": fam, "n_sites": len(g),
             "share_of_class": len(g) / n_cls[cls], "net": float(np.median(v)),
             "lo": float(np.percentile(b, 2.5)), "hi": float(np.percentile(b, 97.5))}
        if cls != "Stable":
            st = d[(d.cls == "Stable") & (d.family == fam)].net
            if len(st) >= a.min_sites:
                r["p_vs_stable"] = mannwhitneyu(st, v).pvalue
        rows.append(r)
T = pd.DataFrame(rows)
for cls, g in T.groupby("class"):
    m = g.p_vs_stable.notna() if "p_vs_stable" in g else pd.Series(False, index=g.index)
    if m.any():
        T.loc[g.index[m], "fdr_vs_stable"] = multipletests(g.loc[m, "p_vs_stable"], method="fdr_bh")[1]
    T.loc[g.index, "rank_in_class"] = g.net.rank().astype(int)      # most negative = rank 1
T = T.sort_values(["class", "rank_in_class"])
T.to_csv(TAB / "family_rank_by_class.csv", index=False)

print("\n=== top partner families per class (most negative net = most depended on), AP-1 excluded")
for cls in CLS5:
    g = T[T["class"] == cls].head(5)
    if g.empty:
        continue
    print(f"\n  {cls}  (n = {n_cls[cls]} sites)")
    for r in g.itertuples():
        extra = ""
        if not np.isnan(getattr(r, "fdr_vs_stable", np.nan)):
            extra = f", vs Stable FDR {r.fdr_vs_stable:.2g}"
        print(f"    {r.rank_in_class}. {r.family:<26} net {r.net:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]  "
              f"in {100 * r.share_of_class:4.1f}% of sites (n={r.n_sites}){extra}")

sig = T[T.get("fdr_vs_stable", pd.Series(dtype=float)) < 0.05] if "fdr_vs_stable" in T else T.iloc[:0]
print(f"\n{len(sig)} family x class comparisons differ from Stable at FDR < 0.05"
      + ("" if len(sig) else "  <- no partner family is used differently by any Dynamic class"))
if len(sig):
    print(sig[["class", "family", "net", "fdr_vs_stable"]].round(3).to_string(index=False))
print("\nwrote", TAB / "family_rank_by_class.csv")
