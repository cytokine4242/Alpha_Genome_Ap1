#!/usr/bin/env python3
"""The sliding-window version of the importance density: a 10 bp edit stepped 2 bp from -400 to +390 bp.

57 did this on the 60-tile scan, which samples the locus once per 10 bp. At that resolution a difference confined
to the motif can only ever appear as one or two points, and the width of anything is unmeasurable. This reads the
sliding run instead (396 windows a site) and asks the two questions the tile version could not:

  how wide is the region the classes differ over, and where exactly does it start and stop?
  how sharp is the motif peak itself - does the dependence fall off over 10 bp or over 50?

Everything else matches 57 so the two are comparable: each site's own noise floor is subtracted before its profile
is normalised, because a less accessible site is a noisier one and noise is flat across position.

  1  per site: signal = max(|effect| - its own floor, 0), then density = signal / the site's total
  2  per position and class: median density, bootstrap band, Mann-Whitney and BH across positions
  3  the shape: full width at half maximum of the motif peak, per site and per class
  4  the difference: where Stable - Dynamic is significant, and how wide that run is
  5  concentration within +-25 bp, unmatched and accessibility-matched, as in 57

Outputs Analysis/fb_scan/tables/slide_density_*.csv.
  python3 58_fb_slide_density.py [--res DIR] [--boot 2000]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--caliper", type=float, default=0.1)
ap.add_argument("--seed", type=int, default=58)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_scan" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)

BASE_FROM = 300          # |pos| at or beyond which the profile is taken to be noise
NEAR_BP = 25             # the concentration window, as in 49/57

cols = ["site_id", "group", "stable_dynamic", "kind", "pos", "rep", "trackset", a.readout]
files = sorted(glob.glob(str(RES / "slide_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no slide_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
sites = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                    usecols=["site_id", "age_class", "ex_class"]).drop_duplicates("site_id")
d = d.merge(sites, on="site_id", how="left")

nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
t = d[d.kind == "slide"].copy()
t["eff"] = np.log2((t[a.readout] + pc) / (t.site_id.map(nat) + pc))
key = ["site_id", "stable_dynamic", "age_class"]
t = t.groupby(key + ["pos"], as_index=False).eff.mean()
t["imp"] = t.eff.abs()

floor = t[t.pos.abs() >= BASE_FROM].groupby("site_id").imp.median().rename("floor")
t = t.join(floor, on="site_id")
t["signal"] = (t.imp - t.floor).clip(lower=0)
tot = t.groupby("site_id").signal.sum().rename("total")
t = t.join(tot, on="site_id")
t = t[t.total > 0].copy()
t["density"] = t.signal / t.total

POS = np.sort(t.pos.unique())
step = int(np.median(np.diff(POS)))
print(f"{len(files)} chunks, {t.site_id.nunique()} sites, {len(POS)} windows "
      f"({POS.min()} to {POS.max()} bp, step {step}), readout {a.readout}")
fl = t.drop_duplicates("site_id").groupby("stable_dynamic").floor.median()
print("  per-site noise floor beyond +-%d bp: %s  <- subtracted before normalising"
      % (BASE_FROM, ", ".join(f"{k} {v:.4f}" for k, v in fl.items())))


def band(M, boot):
    med = np.nanmedian(M, axis=0)
    draws = np.empty((boot, M.shape[1]))
    for b in range(boot):
        draws[b] = np.nanmedian(M[rng.integers(0, M.shape[0], M.shape[0])], axis=0)
    lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    return med, lo, hi


# ---------------------------------------------------------------- 2. profile by class
rows = []
for cls, sub in t.groupby("stable_dynamic"):
    M = sub.pivot(index="site_id", columns="pos", values="density").reindex(columns=POS).to_numpy(np.float64)
    med, lo, hi = band(M, a.boot)
    for j, p_ in enumerate(POS):
        rows.append({"class": cls, "pos": int(p_), "n_sites": M.shape[0],
                     "median": med[j], "lo": lo[j], "hi": hi[j]})
P = pd.DataFrame(rows)

S, D = t[t.stable_dynamic == "Stable"], t[t.stable_dynamic == "Dynamic"]
Sp = S.pivot(index="site_id", columns="pos", values="density").reindex(columns=POS)
Dp = D.pivot(index="site_id", columns="pos", values="density").reindex(columns=POS)
pv = [mannwhitneyu(Sp[p_].dropna(), Dp[p_].dropna()).pvalue for p_ in POS]
fdr = multipletests(pv, method="fdr_bh")[1]
P = P.merge(pd.DataFrame({"pos": POS.astype(int), "p_stable_vs_dynamic": pv, "fdr": fdr}), on="pos", how="left")
P.to_csv(OUT / "slide_density_profile.csv", index=False)

sig = POS[fdr < 0.05]
print(f"\n=== per-position Stable vs Dynamic: {len(sig)} of {len(POS)} windows at FDR < 0.05")
if len(sig):
    # the longest unbroken run of significant windows is the region the classes actually differ over
    runs, cur = [], [sig[0]]
    for x in sig[1:]:
        if x - cur[-1] <= step:
            cur.append(x)
        else:
            runs.append(cur); cur = [x]
    runs.append(cur)
    runs.sort(key=len, reverse=True)
    print("    significant windows span %+d to %+d bp" % (sig.min(), sig.max()))
    for r in runs[:3]:
        print(f"    run: {r[0]:+d} to {r[-1]:+d} bp  ({len(r)} windows, {r[-1] - r[0] + 10} bp of sequence)")
    pd.DataFrame([{"run": i, "start": r[0], "end": r[-1], "n_windows": len(r),
                   "bp_covered": r[-1] - r[0] + 10} for i, r in enumerate(runs)]
                 ).to_csv(OUT / "slide_density_sig_runs.csv", index=False)


# ---------------------------------------------------------------- 3. how sharp is the motif peak
def fwhm(vals, pos):
    """Full width at half maximum of the central peak, read off the sampled profile."""
    v = np.asarray(vals, float)
    k = int(np.nanargmax(np.where(np.abs(pos) <= 60, v, -np.inf)))
    half = v[k] / 2.0
    lo = k
    while lo > 0 and v[lo] > half:
        lo -= 1
    hi = k
    while hi < len(v) - 1 and v[hi] > half:
        hi += 1
    return pos[hi] - pos[lo], pos[k]


wrows = []
for cls in ("Stable", "Dynamic"):
    M = (Sp if cls == "Stable" else Dp)
    w, apex = fwhm(np.nanmedian(M.to_numpy(np.float64), axis=0), POS)
    per = []
    for _, r in M.iterrows():
        vals = r.to_numpy(np.float64)
        if np.isfinite(vals).all() and vals.max() > 0:
            per.append(fwhm(vals, POS)[0])
    per = np.array(per, float)
    wrows.append({"class": cls, "fwhm_of_class_median_bp": w, "apex_bp": apex,
                  "fwhm_per_site_median_bp": np.median(per), "n_sites": len(per)})
W = pd.DataFrame(wrows)
pw = mannwhitneyu(
    [fwhm(r.to_numpy(np.float64), POS)[0] for _, r in Sp.iterrows() if np.isfinite(r).all()],
    [fwhm(r.to_numpy(np.float64), POS)[0] for _, r in Dp.iterrows() if np.isfinite(r).all()]).pvalue
W["p_stable_vs_dynamic"] = pw
W.to_csv(OUT / "slide_density_width.csv", index=False)
print("\n=== how sharp the motif peak is")
print(W.round(3).to_string(index=False))

# ---------------------------------------------------------------- 5. concentration, as in 57
near = t[t.pos.abs() <= NEAR_BP].groupby("site_id").density.sum().rename("near")
conc = t.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic"]].join(near)
conc["native"] = conc.index.map(nat)
p_near = mannwhitneyu(conc.near[conc.stable_dynamic == "Stable"].dropna(),
                      conc.near[conc.stable_dynamic == "Dynamic"].dropna()).pvalue
print(f"\n=== share of the density within +-{NEAR_BP} bp")
print(conc.groupby("stable_dynamic").near.median().round(4).to_string())
print(f"    Stable vs Dynamic p = {p_near:.3g}")

dm = conc[conc.stable_dynamic == "Dynamic"].dropna(subset=["native", "near"]).copy()
sm = conc[conc.stable_dynamic == "Stable"].dropna(subset=["native", "near"]).copy()
dm["l"], sm["l"] = np.log10(dm.native.clip(lower=1e-3)), np.log10(sm.native.clip(lower=1e-3))
avail = sm.sort_values("l").reset_index(drop=True)
used, pairs = set(), []
for r in dm.sort_values("l").itertuples():
    free = avail[~avail.index.isin(used)]
    if free.empty:
        break
    j = (free.l - r.l).abs().idxmin()
    if abs(avail.loc[j, "l"] - r.l) <= a.caliper:
        used.add(j)
        pairs.append((r.near, avail.loc[j, "near"]))
print(f"\n=== accessibility-matched (caliper {a.caliper} log10): {len(pairs)} pairs")
row = {"n_pairs": len(pairs), "p_unmatched": p_near}
if len(pairs) > 30:
    dn = np.array([s - dd for dd, s in pairs])
    p_m = wilcoxon(dn).pvalue
    print(f"    Dynamic {np.median([x[0] for x in pairs]):.4f}  Stable {np.median([x[1] for x in pairs]):.4f}  "
          f"paired {np.median(dn):+.4f}, Wilcoxon p = {p_m:.3g}")
    row.update({"near_dynamic": float(np.median([x[0] for x in pairs])),
                "near_stable": float(np.median([x[1] for x in pairs])),
                "paired_difference": float(np.median(dn)), "p_matched": p_m})
pd.DataFrame([row]).to_csv(OUT / "slide_density_concentration.csv", index=False)
print("\nwrote", OUT / "slide_density_profile.csv", "| _width.csv | _concentration.csv | _sig_runs.csv")
