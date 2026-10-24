#!/usr/bin/env python3
"""Does pooling every Dynamic site hide anything? The fibroblast tile scan split by the binding categories as given in the
site table, rather than Stable vs Dynamic (Katana/results_ft/scan_fb_chunk*.csv, fine-tuned 1 bp fibroblast head).

For every class it reports the two quantities the AP-1-centric / distributed-context question turns on:
  motif dependence   the tiles covering the 11 bp motif (positions -10 and 0)
  collar dependence  the tiles from 10 to 50 bp out, where the Stable/Dynamic difference was found to live
  motif share        motif / (motif + collar): high = the site's dependence is concentrated on AP-1,
                     low = it is spread into the surrounding sequence
Each site's own noise floor (median |effect| at |pos| >= 200 bp) is subtracted first, because Dynamic sites are less
accessible and therefore noisier, which would otherwise inflate their apparent context dependence.
Outputs Analysis/fb_scan/tables/per_class_*.csv.   python3 50_fb_scan_per_class.py [--res DIR]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_scan" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
MOTIF_TILES, COLLAR = (-10, 0), (10, 50)

files = sorted(glob.glob(str(RES / "scan_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no scan_fb_chunk*.csv in {RES}")
cols = ["site_id", "group", "binding_category", "stable_dynamic", "condition", "kind", "pos", "rep", "trackset", a.readout]
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
sites = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                    usecols=["site_id", "age_class", "ex_class"]).drop_duplicates("site_id")
d = d.merge(sites, on="site_id", how="left")

nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
t = d[d.kind == "tile"].copy()
t["eff"] = np.log2((t[a.readout] + pc) / (t.site_id.map(nat) + pc))
t = t.groupby(["site_id", "binding_category", "stable_dynamic", "age_class", "ex_class", "pos"], as_index=False).eff.mean()
floor = t[t.pos.abs() >= 200].assign(x=lambda z: z.eff.abs()).groupby("site_id").x.median().rename("floor")
t = t.join(floor, on="site_id")
t["signal"] = (t.eff.abs() - t.floor).clip(lower=0)

key = ["site_id", "binding_category", "stable_dynamic", "age_class", "ex_class"]
motif = t[t.pos.isin(MOTIF_TILES)].groupby(key).signal.sum().rename("motif")
collar = t[(t.pos.abs() >= COLLAR[0]) & (t.pos.abs() < COLLAR[1])].groupby(key).signal.sum().rename("collar")
per = pd.concat([motif, collar], axis=1).reset_index()
per["motif_share"] = per.motif / (per.motif + per.collar).replace(0, np.nan)
per["native"] = per.site_id.map(nat)
per.to_csv(OUT / "per_class_per_site.csv", index=False)


def table(field, order=None):
    rows = []
    for cls, g in per.groupby(field):
        if len(g) < 15:
            continue
        rows.append({"grouping": field, "class": cls, "n_sites": len(g), "native_median": g.native.median(),
                     "motif": g.motif.median(), "collar": g.collar.median(), "motif_share": g.motif_share.median()})
    r = pd.DataFrame(rows)
    if order:
        r["_o"] = r["class"].map({k: i for i, k in enumerate(order)}).fillna(99)
        r = r.sort_values("_o").drop(columns="_o")
    return r.sort_values("n_sites", ascending=False) if not order else r


print(f"{len(files)} chunks, {per.site_id.nunique()} sites, readout {a.readout}, trackset {a.track}")
print(f"motif tiles {MOTIF_TILES}, collar {COLLAR[0]}-{COLLAR[1]} bp, noise floor subtracted per site\n")
out = []
for field, order in (("stable_dynamic", ["Stable", "Dynamic"]), ("binding_category", None),
                     ("age_class", ["Stable", "Age Up", "Age Down", "Unbound"]),
                     ("ex_class", ["Stable", "Exercise Up", "Exercise Down", "Unbound"])):
    r = table(field, order)
    if r.empty:
        continue
    out.append(r)
    print(f"=== by {field}")
    print(r.round(4).to_string(index=False))
    # does this grouping separate the sites at all, beyond the pooled Stable/Dynamic split?
    groups = [g.motif_share.dropna().values for _, g in per.groupby(field) if len(g) >= 15]
    if len(groups) > 2:
        print(f"    Kruskal-Wallis across classes, motif share: p={kruskal(*groups).pvalue:.4f}")
    gc = [g.collar.dropna().values for _, g in per.groupby(field) if len(g) >= 15]
    if len(gc) > 2:
        print(f"    Kruskal-Wallis across classes, collar:      p={kruskal(*gc).pvalue:.4f}")
    print()
pd.concat(out).to_csv(OUT / "per_class_summary.csv", index=False)

print("=== do the Dynamic sub-classes differ from each other? (Dynamic sites only)")
dyn = per[per.stable_dynamic == "Dynamic"]
sub = [(c, g) for c, g in dyn.groupby("binding_category") if len(g) >= 15]
for name, metric in (("motif", "motif"), ("collar", "collar"), ("motif share", "motif_share")):
    vals = [g[metric].dropna().values for _, g in sub]
    if len(vals) > 2:
        print(f"  {name:12s} Kruskal-Wallis across {len(vals)} Dynamic categories: p={kruskal(*vals).pvalue:.4f}")
print("  medians per Dynamic category:")
print(pd.DataFrame([{"category": c, "n": len(g), "motif": g.motif.median(), "collar": g.collar.median(),
                     "motif_share": g.motif_share.median()} for c, g in sub]).round(4).to_string(index=False))

# ---------------------------------------------------------------- accessibility-matched comparison
# The Stable sites here were drawn at random, so they are about twice as accessible as the Dynamic ones. The earlier
# cardiomyocyte work found that the concentration differences largely disappeared once Stable sites were matched on
# accessibility, so the same check is run here: each Dynamic site is paired with the unused Stable site closest in
# log10 native accessibility (caliper 0.1 log10), and the pairs are compared with a paired test.
from scipy.stats import wilcoxon

dyn_m = per[per.stable_dynamic == "Dynamic"].dropna(subset=["native", "motif_share"]).copy()
sta_m = per[per.stable_dynamic == "Stable"].dropna(subset=["native", "motif_share"]).copy()
dyn_m["l"], sta_m["l"] = np.log10(dyn_m.native.clip(lower=1e-3)), np.log10(sta_m.native.clip(lower=1e-3))
avail = sta_m.sort_values("l").reset_index(drop=True)
used, pairs = set(), []
for r in dyn_m.sort_values("l").itertuples():
    free = avail[~avail.index.isin(used)]
    if free.empty:
        break
    j = (free.l - r.l).abs().idxmin()
    if abs(avail.loc[j, "l"] - r.l) <= 0.1:
        used.add(j)
        pairs.append((r, avail.loc[j]))
print(f"\n=== accessibility-matched Stable (caliper 0.1 log10): {len(pairs)} pairs of {len(dyn_m)} Dynamic sites")
if len(pairs) > 30:
    D = pd.DataFrame([{"motif": a_.motif, "collar": a_.collar, "motif_share": a_.motif_share, "native": a_.native} for a_, _ in pairs])
    S = pd.DataFrame([{"motif": b_.motif, "collar": b_.collar, "motif_share": b_.motif_share, "native": b_.native} for _, b_ in pairs])
    print(f"    native accessibility after matching: Dynamic {D.native.median():.2f}  Stable {S.native.median():.2f}")
    rows = []
    for metric in ("motif", "collar", "motif_share"):
        diff = (D[metric] - S[metric]).dropna()
        rows.append({"metric": metric, "n_pairs": len(diff), "median_dynamic": D[metric].median(),
                     "median_stable_matched": S[metric].median(), "paired_median_difference": diff.median(),
                     "p_wilcoxon": wilcoxon(diff).pvalue if (diff != 0).any() else np.nan})
    M = pd.DataFrame(rows)
    M.to_csv(OUT / "per_class_accessibility_matched.csv", index=False)
    print(M.round(4).to_string(index=False))

print("\n=== pooled Stable vs Dynamic, for reference")
s, dy = per[per.stable_dynamic == "Stable"], per[per.stable_dynamic == "Dynamic"]
for metric in ("motif", "collar", "motif_share"):
    print(f"  {metric:12s} Stable {s[metric].median():.4f}  Dynamic {dy[metric].median():.4f}  "
          f"p={mannwhitneyu(s[metric].dropna(), dy[metric].dropna()).pvalue:.2e}")
print("\noutputs in", OUT)
