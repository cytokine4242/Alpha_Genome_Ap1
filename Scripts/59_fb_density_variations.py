#!/usr/bin/env python3
"""Three other ways of showing the same importance-density result, because two near-identical curves hide it.

The profile in 57 is honest but hard to read: Stable and Dynamic differ by about a percentage point per position,
which is invisible when both are drawn at full height. Nothing here changes the measurement. What changes is what
is plotted:

  1  cumulative   the share of a site's importance within +-x bp, as a function of x. Per-position differences of
                  a point each accumulate, so a real gap becomes a visible one. The quantity at x = 25 is exactly
                  the number 57 already reports, so the two cannot disagree.
  2  difference   Stable minus Dynamic per position, on its own axis, with a band from resampling sites within
                  each class. One curve instead of two, and the band says which part of it is real.
  3  by class     the same two views with Dynamic split into the age and exercise classes as given, which is the
                  only one of the three that shows anything the pooled comparison cannot.

An axis that starts above zero or is normalised can mislead, so every output carries the raw per-site numbers
alongside, and the figures built from them print the effect size rather than only the shape.

Reads the tile scan by default, or the sliding scan with --kind slide once that exists.
Outputs Analysis/fb_scan/tables/density_{cumulative,difference,profile_5class,concentration_5class}.csv.
  python3 59_fb_density_variations.py [--kind scan|slide] [--boot 2000]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--kind", default="scan", choices=["scan", "slide"])
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=59)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_scan" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)

EDIT_KIND = "tile" if a.kind == "scan" else "slide"
BASE_FROM = 200 if a.kind == "scan" else 300
CLS5 = ["Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down"]

cols = ["site_id", "stable_dynamic", "kind", "pos", "rep", "trackset", a.readout]
files = sorted(glob.glob(str(RES / f"{a.kind}_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no {a.kind}_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
sites = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                    usecols=["site_id", "age_class", "ex_class"]).drop_duplicates("site_id")
d = d.merge(sites, on="site_id", how="left")

nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
t = d[d.kind == EDIT_KIND].copy()
t["eff"] = np.log2((t[a.readout] + pc) / (t.site_id.map(nat) + pc))
t = t.groupby(["site_id", "stable_dynamic", "age_class", "ex_class", "pos"], as_index=False).eff.mean()
t["imp"] = t.eff.abs()

floor = t[t.pos.abs() >= BASE_FROM].groupby("site_id").imp.median().rename("floor")
t = t.join(floor, on="site_id")
t["signal"] = (t.imp - t.floor).clip(lower=0)
tot = t.groupby("site_id").signal.sum().rename("total")
t = t.join(tot, on="site_id")
t = t[t.total > 0].copy()
t["density"] = t.signal / t.total

# the five classes as given; a site dynamic with both age and exercise is counted under its age class
t["cls"] = np.where(t.stable_dynamic == "Stable", "Stable",
            np.where(t.age_class == "Age Up", "Age Up",
            np.where(t.age_class == "Age Down", "Age Down",
            np.where(t.ex_class == "Exercise Up", "Exercise Up",
            np.where(t.ex_class == "Exercise Down", "Exercise Down", None)))))
t = t[t.cls.notna()]
POS = np.sort(t.pos.unique())
step = int(np.median(np.diff(POS))) if len(POS) > 1 else 10
print(f"{len(files)} chunks, {t.site_id.nunique()} sites, {len(POS)} windows "
      f"({POS.min()} to {POS.max()} bp, step {step}), edit kind '{EDIT_KIND}'")
print("  sites per class: " + ", ".join(f"{k} {v}" for k, v in
                                        t.drop_duplicates("site_id").cls.value_counts().items()))

MAT = {c: g.pivot(index="site_id", columns="pos", values="density").reindex(columns=POS).to_numpy(np.float64)
       for c, g in t.groupby("cls")}
SD = {"Stable": MAT["Stable"],
      "Dynamic": np.vstack([MAT[c] for c in CLS5[1:] if c in MAT])}


def boot_median(M, boot):
    med = np.nanmedian(M, axis=0)
    draws = np.empty((boot, M.shape[1]))
    for b in range(boot):
        draws[b] = np.nanmedian(M[rng.integers(0, M.shape[0], M.shape[0])], axis=0)
    return med, draws


# ---------------------------------------------------------------- 1. cumulative share within +-x
XS = np.array([x for x in POS if x >= 0])
rows = []
for name, M in {**SD, **MAT}.items():
    C = np.stack([np.nansum(M[:, np.abs(POS) <= x], axis=1) for x in XS], axis=1)   # sites x radii
    med, draws = boot_median(C, a.boot)
    lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    for j, x in enumerate(XS):
        rows.append({"grouping": "pooled" if name in SD else "class", "class": name, "radius_bp": int(x),
                     "n_sites": M.shape[0], "median": med[j], "lo": lo[j], "hi": hi[j]})
CU = pd.DataFrame(rows)
CU.to_csv(OUT / "density_cumulative.csv", index=False)
at25 = CU[(CU.radius_bp == 25)].set_index("class")["median"]
print("\n=== cumulative share of a site's importance within +-25 bp")
print(at25.round(4).to_string())

# ---------------------------------------------------------------- 2. the difference, with its own band
mS, dS = boot_median(SD["Stable"], a.boot)
mD, dD = boot_median(SD["Dynamic"], a.boot)
diff = 100 * (mS - mD)
dd = 100 * (dS - dD)
lo, hi = np.nanpercentile(dd, [2.5, 97.5], axis=0)
pv = [mannwhitneyu(SD["Stable"][:, j][~np.isnan(SD["Stable"][:, j])],
                   SD["Dynamic"][:, j][~np.isnan(SD["Dynamic"][:, j])]).pvalue for j in range(len(POS))]
fdr = multipletests(pv, method="fdr_bh")[1]
DF = pd.DataFrame({"pos": POS.astype(int), "diff_pp": diff, "lo": lo, "hi": hi,
                   "p_stable_vs_dynamic": pv, "fdr": fdr})
DF.to_csv(OUT / "density_difference.csv", index=False)
sig = POS[fdr < 0.05]
print(f"\n=== Stable − Dynamic per position: {len(sig)} of {len(POS)} at FDR < 0.05"
      + (f", spanning {sig.min():+d} to {sig.max():+d} bp" if len(sig) else ""))
print(f"    largest difference {diff.max():+.2f} points at {POS[int(np.argmax(diff))]:+d} bp")

# ---------------------------------------------------------------- 3. the five classes
rows = []
for c in CLS5:
    if c not in MAT:
        continue
    med, draws = boot_median(MAT[c], a.boot)
    lo_, hi_ = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    for j, p_ in enumerate(POS):
        rows.append({"class": c, "pos": int(p_), "n_sites": MAT[c].shape[0],
                     "median": med[j], "lo": lo_[j], "hi": hi_[j]})
pd.DataFrame(rows).to_csv(OUT / "density_profile_5class.csv", index=False)

near = {c: np.nansum(M[:, np.abs(POS) <= 25], axis=1) for c, M in MAT.items()}
rows = []
for c in CLS5:
    if c not in near:
        continue
    v = near[c]
    b = np.array([np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(a.boot)])
    rows.append({"class": c, "n_sites": len(v), "share_within_25bp": float(np.median(v)),
                 "lo": float(np.percentile(b, 2.5)), "hi": float(np.percentile(b, 97.5)),
                 "p_vs_stable": np.nan if c == "Stable" else
                 mannwhitneyu(near["Stable"], v).pvalue})
C5 = pd.DataFrame(rows)
m = C5.p_vs_stable.notna()
C5.loc[m, "fdr_vs_stable"] = multipletests(C5.loc[m, "p_vs_stable"], method="fdr_bh")[1]
C5.to_csv(OUT / "density_concentration_5class.csv", index=False)
print("\n=== concentration within +-25 bp, by class")
print(C5.round(4).to_string(index=False))
print("\nwrote", OUT / "density_cumulative.csv", "| _difference.csv | _profile_5class.csv | _concentration_5class.csv")
