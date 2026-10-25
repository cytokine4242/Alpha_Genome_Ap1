#!/usr/bin/env python3
"""Importance as a density along the sequence: where, across 600 bp, does a fibroblast AP-1 site's dependence sit?

49 writes the tile profile already, but as a bare per-class median with no interval, and its `shape` column divides
each site by its own total WITHOUT subtracting the noise floor first. That matters here more than anywhere else:
Dynamic sites are about half as accessible as Stable ones, so their log2 ratios are noisier, and noise is flat across
position. An uncorrected density therefore makes Dynamic profiles look flatter - more "distributed" - whether or not
they are. Everything below is computed on the floor-corrected signal, the uncorrected version is written alongside so
the two can be compared, and the headline concentration is also reported on accessibility-matched sites.

  1  per site: signal = max(|effect| - the site's own noise floor, 0), then density = signal / the site's total
  2  per position and class: median density, a bootstrap 95% band over sites, Mann-Whitney and BH-FDR across positions
  3  concentration: the share of a site's density within +-25 bp, and the half-width, unmatched and matched
Both groupings are written: Stable vs Dynamic, and Dynamic split into the age classes as given.

Outputs Analysis/fb_scan/tables/importance_density_*.csv.
  python3 57_fb_importance_density.py [--res DIR] [--boot 2000]
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
ap.add_argument("--seed", type=int, default=57)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_scan" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)

BASE_FROM = 200          # |pos| at or beyond which the profile is taken to be noise, as in 49
NEAR_BP = 25             # the concentration window, as in 49

# ---------------------------------------------------------------- load
cols = ["site_id", "group", "binding_category", "stable_dynamic", "kind", "pos", "rep", "trackset", a.readout]
files = sorted(glob.glob(str(RES / "scan_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no scan_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
sites = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                    usecols=["site_id", "age_class", "ex_class"]).drop_duplicates("site_id")
d = d.merge(sites, on="site_id", how="left")

nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
t = d[d.kind == "tile"].copy()
t["eff"] = np.log2((t[a.readout] + pc) / (t.site_id.map(nat) + pc))
key = ["site_id", "stable_dynamic", "age_class"]
t = t.groupby(key + ["pos"], as_index=False).eff.mean()
t["imp"] = t.eff.abs()

# ---------------------------------------------------------------- 1. per-site density
floor = t[t.pos.abs() >= BASE_FROM].groupby("site_id").imp.median().rename("floor")
t = t.join(floor, on="site_id")
t["signal"] = (t.imp - t.floor).clip(lower=0)
tot = t.groupby("site_id").signal.sum().rename("total")
t = t.join(tot, on="site_id")
t = t[t.total > 0].copy()
t["density"] = t.signal / t.total
raw_tot = t.groupby("site_id").imp.sum().rename("raw_total")
t = t.join(raw_tot, on="site_id")
t["density_raw"] = t.imp / t.raw_total

t["grp4"] = np.where(t.stable_dynamic == "Stable", "Stable",
                     np.where(t.age_class == "Age Down", "AP-1 loss",
                              np.where(t.age_class == "Age Up", "AP-1 gain", "other Dynamic")))
POS = np.sort(t.pos.unique())
print(f"{len(files)} chunks, {t.site_id.nunique()} sites kept, {len(POS)} positions "
      f"({POS.min()} to {POS.max()} bp), readout {a.readout}")
fl = t.drop_duplicates("site_id").groupby("stable_dynamic").floor.median()
print("  per-site noise floor (median |effect| beyond +-%d bp): %s  <- subtracted before normalising"
      % (BASE_FROM, ", ".join(f"{k} {v:.4f}" for k, v in fl.items())))


def matrix(sub, col):
    """sites x positions, in POS order - the shape every bootstrap below works on."""
    m = sub.pivot(index="site_id", columns="pos", values=col).reindex(columns=POS)
    return m.to_numpy(np.float64), m.index.to_numpy()


def band(M, boot):
    """median density per position, with a percentile bootstrap band over sites."""
    med = np.nanmedian(M, axis=0)
    draws = np.empty((boot, M.shape[1]))
    for b in range(boot):
        draws[b] = np.nanmedian(M[rng.integers(0, M.shape[0], M.shape[0])], axis=0)
    lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    return med, lo, hi


# ---------------------------------------------------------------- 2. the profile, by class
rows = []
for field in ("stable_dynamic", "grp4"):
    for cls, sub in t.groupby(field):
        if sub.site_id.nunique() < 30:
            continue
        for col, tag in (("density", "corrected"), ("density_raw", "raw")):
            M, _ = matrix(sub, col)
            med, lo, hi = band(M, a.boot if tag == "corrected" else 400)
            for j, p_ in enumerate(POS):
                rows.append({"grouping": field, "class": cls, "version": tag, "pos": int(p_),
                             "n_sites": M.shape[0], "median": med[j], "lo": lo[j], "hi": hi[j]})
P = pd.DataFrame(rows)

# per-position Stable vs Dynamic, on the corrected density
S = t[t.stable_dynamic == "Stable"]
D = t[t.stable_dynamic == "Dynamic"]
pv = []
for p_ in POS:
    x = S.density[S.pos == p_].to_numpy()
    y = D.density[D.pos == p_].to_numpy()
    pv.append(mannwhitneyu(x, y).pvalue)
fdr = multipletests(pv, method="fdr_bh")[1]
PP = pd.DataFrame({"pos": POS.astype(int), "p_stable_vs_dynamic": pv, "fdr": fdr})
P = P.merge(PP, on="pos", how="left")
P.to_csv(OUT / "importance_density_profile.csv", index=False)
print(f"\n=== per-position Stable vs Dynamic (corrected density): {int((fdr < 0.05).sum())} of {len(POS)} "
      f"positions at FDR < 0.05")
sig = PP[PP.fdr < 0.05]
if len(sig):
    print("    " + ", ".join(f"{int(r.pos)}" for r in sig.itertuples()))

# ---------------------------------------------------------------- 3. concentration
near = t[t.pos.abs() <= NEAR_BP].groupby("site_id").density.sum().rename("near")
conc = t.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic", "grp4"]].join(near)
conc["native"] = conc.index.map(nat)
t["absdist"] = t.pos.abs()


def half_width(g):
    s = g.sort_values("absdist")
    c = s.density.cumsum()
    hit = s.absdist[c >= 0.5]
    return hit.iloc[0] if len(hit) else np.nan


conc["half_width_bp"] = t.groupby("site_id", group_keys=False)[["absdist", "density"]].apply(half_width)
conc.to_csv(OUT / "importance_density_per_site.csv")

rows = []
for field in ("stable_dynamic", "grp4"):
    for cls, g in conc.groupby(field):
        if len(g) < 30:
            continue
        rows.append({"grouping": field, "class": cls, "n_sites": len(g),
                     "share_within_25bp": g.near.median(), "half_width_bp": g.half_width_bp.median(),
                     "native_median": g.native.median()})
C = pd.DataFrame(rows)
p_near = mannwhitneyu(conc.near[conc.stable_dynamic == "Stable"],
                      conc.near[conc.stable_dynamic == "Dynamic"]).pvalue
p_hw = mannwhitneyu(conc.half_width_bp[conc.stable_dynamic == "Stable"].dropna(),
                    conc.half_width_bp[conc.stable_dynamic == "Dynamic"].dropna()).pvalue
print("\n=== concentration of the corrected density")
print(C.round(4).to_string(index=False))
print(f"    Stable vs Dynamic: share within +-{NEAR_BP} bp p = {p_near:.3g}, half-width p = {p_hw:.3g}")

# accessibility-matched, the check the ratio measures in 50 failed
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
        pairs.append((r, avail.loc[j]))
print(f"\n=== accessibility-matched Stable (caliper {a.caliper} log10): {len(pairs)} pairs of {len(dm)} Dynamic sites")
mrow = {}
if len(pairs) > 30:
    PD = pd.DataFrame([{"near": x.near, "hw": x.half_width_bp, "native": x.native} for x, _ in pairs])
    PS = pd.DataFrame([{"near": y.near, "hw": y.half_width_bp, "native": y.native} for _, y in pairs])
    dn = (PS.near - PD.near).dropna()
    p_m = wilcoxon(dn).pvalue
    print(f"    native after matching: Dynamic {PD.native.median():.2f}  Stable {PS.native.median():.2f}")
    print(f"    share within +-{NEAR_BP} bp: Dynamic {PD.near.median():.4f}  Stable {PS.near.median():.4f}  "
          f"paired difference {dn.median():+.4f}, Wilcoxon p = {p_m:.3g}")
    mrow = {"n_pairs": len(pairs), "near_dynamic": PD.near.median(), "near_stable": PS.near.median(),
            "paired_difference": dn.median(), "p_wilcoxon": p_m,
            "native_median_dynamic": PD.native.median(), "native_median_stable": PS.native.median()}
    pd.DataFrame([mrow]).to_csv(OUT / "importance_density_matched.csv", index=False)

C["p_stable_vs_dynamic_share"] = p_near
C["p_stable_vs_dynamic_halfwidth"] = p_hw
C["p_share_matched"] = mrow.get("p_wilcoxon", np.nan)
C.to_csv(OUT / "importance_density_concentration.csv", index=False)
print("\nwrote", OUT / "importance_density_profile.csv", "| _concentration.csv | _per_site.csv | _matched.csv")
