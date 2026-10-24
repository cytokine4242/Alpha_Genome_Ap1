#!/usr/bin/env python3
"""Confidence intervals for the sequence-versus-location regression (47_fb_swap_analysis.py reports the point estimates only).

Each swap writes one donor element into one recipient location, so the 62,250 swaps are not independent: every donor
appears in ~250 of them. A plain OLS standard error would therefore be far too small. This resamples DONORS with
replacement (a cluster bootstrap), refits on each resample, and reports the percentile interval, which is what the
slide error bars need.

Outputs Analysis/fb_swap/tables/sequence_vs_location_ci.csv.
  python3 55_fb_swap_regression_ci.py [--res DIR] [--boot 2000]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=55)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_swap" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)


def lg(x):
    return np.log2(np.asarray(x, np.float64).clip(min=1e-3))


cols = ["donor_id", "recipient_id", "width", "trackset", "self_swap",
        "atac_insert", "donor_native_atac_insert", "recipient_native_atac_insert"]
files = sorted(glob.glob(str(RES / "swap_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no swap_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[(d.trackset == a.track) & (~d.self_swap.astype(bool))]
print(f"{len(files)} chunks, {len(d)} swaps, {d.donor_id.nunique()} donors, {d.recipient_id.nunique()} locations")

rows = []
for w, g in d.groupby("width"):
    y = lg(g.atac_insert)
    X = np.column_stack([np.ones(len(g)), lg(g.donor_native_atac_insert), lg(g.recipient_native_atac_insert)])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

    # cluster bootstrap over donors: redraw the set of donor elements, keep every swap each drawn donor made
    idx_of = {k: np.flatnonzero(g.donor_id.values == k) for k in g.donor_id.unique()}
    donors = np.array(list(idx_of))
    draws = np.empty((a.boot, 2))
    for b in range(a.boot):
        take = np.concatenate([idx_of[k] for k in rng.choice(donors, len(donors), replace=True)])
        draws[b] = np.linalg.lstsq(X[take], y[take], rcond=None)[0][1:]
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
    # how often does the bootstrap put the element above the location, and vice versa?
    p_seq_gt_loc = float((draws[:, 0] > draws[:, 1]).mean())
    rows.append({"width": w, "n_swaps": len(g), "n_donors": len(donors), "r2": r2,
                 "beta_donor_sequence": coef[1], "seq_lo": lo[0], "seq_hi": hi[0],
                 "beta_recipient_location": coef[2], "loc_lo": lo[1], "loc_hi": hi[1],
                 "boot_frac_sequence_gt_location": p_seq_gt_loc,
                 "p_sequence_vs_location": 2 * min(p_seq_gt_loc, 1 - p_seq_gt_loc)})
R = pd.DataFrame(rows)
R.to_csv(OUT / "sequence_vs_location_ci.csv", index=False)
pd.set_option("display.width", 220)
print(f"\ncluster bootstrap over donors, {a.boot} resamples:")
print(R.round(4).to_string(index=False))
print("\nwrote", OUT / "sequence_vs_location_ci.csv")
