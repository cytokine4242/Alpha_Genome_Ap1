#!/usr/bin/env python3
"""Fixed cardiomyocyte site subset for the tile-scan / saturation-mutagenesis and AP-1 titration experiments. No TOBIAS footprint scores are used;
groups are defined from AlphaGenome's native predictions and the peak geometry only (seed 34):
  dynamic            500 of the 772 Dynamic sites
  matched_stable     for each of those, the Stable site closest in native predicted accessibility (170 bp and peak), log peak width and log TSS distance
                     (1:1 without replacement, caliper 0.5 SD): the context-matched comparison
  high_stable        500 Stable sites from the top quartile of native predicted 170 bp accessibility (typical strongly accessible Stable sites)
Rows are shuffled so every PBS chunk gets all three groups. Writes Katana/data/cm_subset_dynamics.csv.   python3 34_make_dynamics_subsets.py"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N_PER_GROUP, SEED = 500, 34
rng = np.random.default_rng(SEED)
s = pd.read_csv(ROOT / "Data" / "derived" / "cm_fos_sites.csv").set_index("site_id")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
s = s.drop(columns=["YS_score", "MAS_score", "MAE_score"]).join(sl[["native_170", "native_peak", "abs_tss"]]).dropna(subset=["native_170", "native_peak"])
s["l_170"], s["l_peak"] = np.log10(s.native_170.clip(lower=1e-3)), np.log10(s.native_peak.clip(lower=1e-3))
s["l_pw"], s["l_tss"] = np.log10(s.peak_width), np.log10(s.abs_tss + 1)
cols = ["l_170", "l_peak", "l_pw", "l_tss"]
Z = (s[cols] - s[cols].mean()) / s[cols].std()
dyn = s[s.stable_dynamic == "Dynamic"].sample(n=N_PER_GROUP, random_state=SEED)
avail = set(s.index[s.stable_dynamic == "Stable"])
matched = []
for sid in rng.permutation(dyn.index.values):
    cand = np.array(sorted(avail))
    d = np.sqrt(((Z.loc[cand].values - Z.loc[sid].values) ** 2).sum(1))
    j = int(np.argmin(d))
    if d[j] < 0.5:
        matched.append(cand[j]); avail.discard(cand[j])
ms = s.loc[matched]
rest = s.loc[sorted(avail)]
high = rest[rest.native_170 >= s[s.stable_dynamic == "Stable"].native_170.quantile(0.75)]
high = high.sample(n=min(N_PER_GROUP, len(high)), random_state=SEED)
parts = [dyn.assign(group="dynamic"), ms.assign(group="matched_stable"), high.assign(group="high_stable")]
sub = pd.concat(parts).drop(columns=cols).reset_index().rename(columns={"index": "site_id"})
sub = sub.sample(frac=1, random_state=SEED).reset_index(drop=True)
sub.to_csv(ROOT / "Katana" / "data" / "cm_subset_dynamics.csv", index=False)
print(sub.groupby("group").agg(n=("site_id", "size"), native_170=("native_170", "median"), native_peak=("native_peak", "median"), peak_width=("peak_width", "median"), abs_tss=("abs_tss", "median")).round(3))
print("matched", len(matched), "of", N_PER_GROUP, "| columns:", [c for c in sub.columns][:12], "...")
