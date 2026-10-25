#!/usr/bin/env python3
"""Does the element-versus-location balance differ BY CLASS, or is it the same for every site?

55_fb_swap_regression_ci.py fits one regression over all 62,250 swaps and reports two weights per width. That says
something about elements in general and nothing about Stable against Dynamic, which is the only contrast this deck
is actually about. If the balance is the same in every class, the swap result is a statement about assay width and
does not belong on a class slide.

Per donor class, and per width:
  beta_seq   how much of the swapped site's accessibility follows the DONOR element's own native level
  beta_loc   how much follows the RECIPIENT location's native level
  transfer   beta_seq / (beta_seq + beta_loc), the share of the level the element carries with it

Each donor appears in ~250 swaps, so every interval and every test here resamples DONORS with replacement
(a cluster bootstrap) rather than rows. The class tests are on the bootstrap difference against Stable, BH-adjusted
across the four Dynamic classes within a width.

Self-swaps are excluded. Readout atac_insert, trackset FB, fine-tuned head.

Outputs Analysis/fb_swap/tables/sequence_vs_location_by_class.csv.
  python3 60_fb_swap_by_class.py [--res DIR] [--boot 2000]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=60)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_swap" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)

CLS5 = ["Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down"]
cols = ["donor_id", "donor_class", "recipient_id", "recipient_class", "width", "self_swap", "trackset",
        "atac_insert", "donor_native_atac_insert", "recipient_native_atac_insert"]

files = sorted(glob.glob(str(RES / "swap_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no swap_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[(d.trackset == a.track) & (~d.self_swap.astype(bool))]
d = d.dropna(subset=["atac_insert", "donor_native_atac_insert", "recipient_native_atac_insert"])

pc = 0.05 * d.atac_insert.median()
lg = lambda v: np.log2(np.asarray(v, float) + pc)
d["y"] = lg(d.atac_insert)
d["dn"] = lg(d.donor_native_atac_insert)
d["rn"] = lg(d.recipient_native_atac_insert)
print(f"{len(files)} chunks, {len(d)} swaps, {d.donor_id.nunique()} donors, widths {sorted(d.width.unique())}")


def fit(sub):
    """OLS of y on the donor's own level and the recipient location's level. Returns (beta_seq, beta_loc, R2)."""
    X = np.column_stack([np.ones(len(sub)), sub.dn.to_numpy(), sub.rn.to_numpy()])
    y = sub.y.to_numpy()
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    ss = ((y - y.mean()) ** 2).sum()
    return b[1], b[2], 1 - (resid ** 2).sum() / ss if ss > 0 else np.nan


rows, draws = [], {}
for w, dw in d.groupby("width"):
    for cls in CLS5:
        sub = dw[dw.donor_class == cls]
        if len(sub) < 50:
            continue
        bs, bl, r2 = fit(sub)
        # resample DONORS, not rows: the ~250 swaps a donor appears in are not independent observations
        donors = sub.donor_id.unique()
        by = {k: g for k, g in sub.groupby("donor_id")}
        B = np.empty((a.boot, 3))
        for i in range(a.boot):
            pick = rng.choice(donors, len(donors), replace=True)
            rs = pd.concat([by[p] for p in pick], ignore_index=True)
            s_, l_, _ = fit(rs)
            B[i] = (s_, l_, s_ / (s_ + l_) if (s_ + l_) != 0 else np.nan)
        draws[(w, cls)] = B
        lo, hi = np.nanpercentile(B, [2.5, 97.5], axis=0)
        rows.append({"width": w, "class": cls, "n_swaps": len(sub), "n_donors": len(donors),
                     "beta_seq": bs, "beta_seq_lo": lo[0], "beta_seq_hi": hi[0],
                     "beta_loc": bl, "beta_loc_lo": lo[1], "beta_loc_hi": hi[1],
                     "transfer": bs / (bs + bl), "transfer_lo": lo[2], "transfer_hi": hi[2], "r2": r2})
T = pd.DataFrame(rows)

# ---- the question the slide has to answer: is any Dynamic class different from Stable?
tests = []
for w in T.width.unique():
    for cls in CLS5[1:]:
        if (w, cls) not in draws or (w, "Stable") not in draws:
            continue
        for j, what in ((0, "beta_seq"), (2, "transfer")):
            dif = draws[(w, "Stable")][:, j] - draws[(w, cls)][:, j]
            dif = dif[np.isfinite(dif)]
            # two-sided bootstrap p: how often the difference crosses zero
            p = 2 * min((dif <= 0).mean(), (dif >= 0).mean())
            tests.append({"width": w, "quantity": what, "class": cls,
                          "stable_minus_class": float(np.median(dif)),
                          "lo": float(np.percentile(dif, 2.5)), "hi": float(np.percentile(dif, 97.5)),
                          "p": max(p, 1.0 / len(dif))})
E = pd.DataFrame(tests)
for (w, q), g in E.groupby(["width", "quantity"]):
    E.loc[g.index, "fdr"] = multipletests(g.p, method="fdr_bh")[1]

T.to_csv(OUT / "sequence_vs_location_by_class.csv", index=False)
E.to_csv(OUT / "sequence_vs_location_class_tests.csv", index=False)
print("\n=== weights by donor class")
print(T.round(3).to_string(index=False))
print("\n=== each Dynamic class against Stable (cluster bootstrap, BH within width and quantity)")
print(E.round(4).to_string(index=False))
sig = E[E.fdr < 0.05]
print(f"\n{len(sig)} of {len(E)} class comparisons reach FDR < 0.05"
      + ("" if len(sig) else "  <- no class difference: the swap result is about assay width, not about class"))
print("\nwrote", OUT / "sequence_vs_location_by_class.csv", "| _class_tests.csv")
