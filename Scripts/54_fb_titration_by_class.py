#!/usr/bin/env python3
"""The fibroblast AP-1 affinity titration, by binding class (Katana/results_ft/titrate_fb_chunk*.csv, fine-tuned 1 bp head).

The scan and shuffle arms ask what a site loses when sequence is destroyed. The titration asks the complementary question:
what does a site do when the AP-1 motif is made weaker, stronger, reversed, duplicated or moved? Conditions come from
15_scan_titrate.py:
  core_mm1..7     the 7 bp core weakened by 1-7 bases      -> per-site dose-response slope and the plateau it saturates at
  core_upgrade    an imperfect core made perfect           -> only at sites that do not already have TGAGTCA
  core_halfsite   only the first half-site kept
  motif_rc        the whole 11 bp motif reverse-complemented
  dup@+-d         a second copy of the motif at +-d bp, against its own shuffled control ctrl@+-d
  move@+-d        the motif relocated to +-d bp,          against its own shuffled control ctrlmove@+-d
The duplication and move arms are the ones that bear on the AP-1-centric / distributed-context question: a site whose
output is already set by its own AP-1 motif has less to gain from a second copy than one running on distributed context.

Every class contrast is reported twice - all sites, and with each Dynamic site paired to an unused Stable site of similar
predicted native accessibility (caliper --caliper log10) - because Dynamic sites are about half as accessible and that
alone changes how much room a site has to move. BH-FDR within each contrast family.
Outputs Analysis/fb_titration/tables.   python3 54_fb_titration_by_class.py [--res DIR]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--caliper", type=float, default=0.1)
ap.add_argument("--min-sites", type=int, default=25)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_titration" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
OFFSETS = [-150, -80, -40, -20, 20, 40, 80, 150]
pd.set_option("display.width", 220)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ---------------------------------------------------------------- effects per site and condition
files = sorted(glob.glob(str(RES / f"titrate_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no titrate_fb_chunk*.csv in {RES}")
cols = ["site_id", "binding_category", "stable_dynamic", "core_found", "condition", "kind", "rep", "trackset", a.readout]
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
sites = pd.read_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv",
                    usecols=["site_id", "age_class", "ex_class"]).drop_duplicates("site_id")
d = d.merge(sites, on="site_id", how="left")

nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
e = d[d.kind != "native"].copy()
e["eff"] = np.log2((e[a.readout] + pc) / (e.site_id.map(nat) + pc))
KEY = ["site_id", "binding_category", "stable_dynamic", "age_class", "ex_class", "core_found"]
per = e.groupby(KEY + ["condition", "kind"], as_index=False).eff.mean()
meta = per.drop_duplicates("site_id").set_index("site_id")[KEY[1:]]
meta["native"] = nat.reindex(meta.index)
print(f"{len(files)} chunks, {per.site_id.nunique()} sites, readout {a.readout}, trackset {a.track}")
print("  " + ", ".join(f"{k} {v}" for k, v in meta.stable_dynamic.value_counts().items())
      + f";  {int(meta.core_found.sum())} sites already carry a perfect TGAGTCA core")

W = per.pivot_table(index="site_id", columns="condition", values="eff")


def net(prefix, ctrl):
    """dup/move minus its own shuffled control, per site per offset, and averaged over offsets."""
    out = {}
    for o in OFFSETS:
        c, k = f"{prefix}@{o:+d}", f"{ctrl}@{o:+d}"
        if c in W and k in W:
            out[o] = W[c] - W[k]
    N = pd.DataFrame(out)
    return N, N.mean(axis=1).rename(f"{prefix}_net")


DUP, dup_mean = net("dup", "ctrl")
MOVE, move_mean = net("move", "ctrlmove")

# per-site dose-response over the 1-7 base core weakening
mm = [f"core_mm{i}" for i in range(1, 8) if f"core_mm{i}" in W]
X = np.array([int(c[-1]) for c in mm], float)
Y = W[mm].values
ok = ~np.isnan(Y).any(axis=1)
slope = np.full(len(W), np.nan)
slope[ok] = np.polyfit(X, Y[ok].T, 1)[0]
S = pd.DataFrame({"mm_slope": slope, "mm_plateau": W[[c for c in mm if int(c[-1]) >= 5]].mean(axis=1),
                  "mm1": W.get("core_mm1"), "halfsite": W.get("core_halfsite"), "motif_rc": W.get("motif_rc"),
                  "core_upgrade": W.get("core_upgrade"), "dup_net": dup_mean, "move_net": move_mean}, index=W.index)
S = S.join(meta)
S.to_csv(OUT / "titration_per_site.csv")

METRICS = [("mm_slope", "dose-response slope, log2 per base of core weakened"),
           ("mm1", "one base of the core changed"),
           ("mm_plateau", "plateau: 5-7 bases changed (the motif is gone)"),
           ("halfsite", "only the first half-site kept"),
           ("motif_rc", "the 11 bp motif reverse-complemented"),
           ("core_upgrade", "imperfect core made perfect (sites without TGAGTCA)"),
           ("dup_net", "a second copy of the motif, net of a shuffled control"),
           ("move_net", "the motif relocated, net of a shuffled control")]

# ---------------------------------------------------------------- 1. medians by class
print(f"\n{'=' * 118}\n=== 1. every condition by class (median log2 effect; negative = the edit costs the site accessibility)")
rows = []
for field, lab in METRICS:
    if field not in S:
        continue
    sub = S.dropna(subset=[field])
    if field == "core_upgrade":
        sub = sub[~sub.core_found.astype(bool)]
    for grouping in ("stable_dynamic", "age_class", "ex_class"):
        for cls, g in sub.groupby(grouping):
            if len(g) >= a.min_sites:
                rows.append({"metric": field, "description": lab, "grouping": grouping, "class": cls,
                             "n": len(g), "median": g[field].median()})
M = pd.DataFrame(rows)
M.to_csv(OUT / "titration_by_class.csv", index=False)
piv = M[M.grouping == "stable_dynamic"].pivot_table(index=["metric", "description"], columns="class", values="median")
cnt = M[M.grouping == "stable_dynamic"].pivot_table(index="metric", columns="class", values="n")
piv = piv.reindex([(f, l) for f, l in METRICS if f in piv.index.get_level_values(0)])
print(piv.round(4).to_string())
print(f"  n sites: Stable {int(cnt.Stable.max())}, Dynamic {int(cnt.Dynamic.max())} "
      f"(fewer for the core upgrade, which only applies where there is no perfect core)")

# ---------------------------------------------------------------- 2. Stable vs Dynamic, all sites and matched
dyn = S[(S.stable_dynamic == "Dynamic")].dropna(subset=["native"]).copy()
sta = S[(S.stable_dynamic == "Stable")].dropna(subset=["native"]).copy()
dyn["l"], sta["l"] = np.log10(dyn.native.clip(lower=1e-3)), np.log10(sta.native.clip(lower=1e-3))
used, pair_of = set(), {}
for sid, r in dyn.sort_values("l").iterrows():
    free = sta[~sta.index.isin(used)]
    if free.empty:
        break
    j = (free.l - r.l).abs().idxmin()
    if abs(sta.loc[j, "l"] - r.l) <= a.caliper:
        used.add(j)
        pair_of[sid] = j
print(f"\n{'=' * 118}\n=== 2. Stable vs Dynamic, unmatched and matched on predicted native accessibility "
      f"(caliper {a.caliper} log10): {len(pair_of)} pairs")
print(f"    native after matching: Dynamic {dyn.loc[list(pair_of), 'native'].median():.2f} vs "
      f"Stable {sta.loc[list(pair_of.values()), 'native'].median():.2f} "
      f"(before: {dyn.native.median():.2f} vs {sta.native.median():.2f})")
rows = []
for field, lab in METRICS:
    if field not in S:
        continue
    keep = (~S.core_found.astype(bool)) if field == "core_upgrade" else pd.Series(True, index=S.index)
    du, su = dyn[keep.reindex(dyn.index, fill_value=False)][field].dropna(), sta[keep.reindex(sta.index, fill_value=False)][field].dropna()
    if min(len(du), len(su)) < a.min_sites:
        continue
    r = {"metric": field, "description": lab, "n_dynamic": len(du), "n_stable": len(su),
         "dynamic": du.median(), "stable": su.median(), "difference": du.median() - su.median(),
         "p_unmatched": mannwhitneyu(du, su).pvalue}
    pa = [(S.at[ds, field], S.at[ss, field]) for ds, ss in pair_of.items()
          if keep.get(ds, False) and keep.get(ss, False) and not (np.isnan(S.at[ds, field]) or np.isnan(S.at[ss, field]))]
    if len(pa) >= a.min_sites:
        A = np.array(pa)
        diff = A[:, 0] - A[:, 1]
        r.update({"n_pairs": len(diff), "paired_difference": float(np.median(diff)),
                  "p_matched": float(wilcoxon(diff).pvalue) if np.any(diff != 0) else np.nan,
                  "min_detectable_80pct": float(2.8 * diff.std() / np.sqrt(len(diff)))})
    rows.append(r)
T = pd.DataFrame(rows)
for c in ("p_unmatched", "p_matched"):
    if c in T:
        T[c.replace("p_", "fdr_")] = bh(T[c].fillna(1).values)
T.to_csv(OUT / "titration_stable_vs_dynamic.csv", index=False)
show = ["metric", "n_pairs", "stable", "dynamic", "difference", "p_unmatched", "paired_difference", "p_matched",
        "fdr_matched", "min_detectable_80pct"]
print(T[[c for c in show if c in T]].round(4).to_string(index=False))
print("  difference = Dynamic - Stable; a positive number means the edit costs the Dynamic site LESS.")

# ---------------------------------------------------------------- 3. per binding class
print(f"\n{'=' * 118}\n=== 3. by age and exercise class (all sites), and whether the classes separate at all")
for grouping, order in (("age_class", ["Stable", "Age Up", "Age Down", "Unbound"]),
                        ("ex_class", ["Stable", "Exercise Up", "Exercise Down", "Unbound"])):
    sub = M[M.grouping == grouping].pivot_table(index="metric", columns="class", values="median")
    sub = sub[[c for c in order if c in sub.columns]].reindex([f for f, _ in METRICS if f in sub.index])
    print(f"\n  by {grouping}")
    print(sub.round(4).to_string())
    for field, _ in METRICS:
        if field not in S:
            continue
        keep = S[~S.core_found.astype(bool)] if field == "core_upgrade" else S
        gs = [g[field].dropna().values for _, g in keep.groupby(grouping) if g[field].notna().sum() >= a.min_sites]
        if len(gs) > 2:
            p = kruskal(*gs).pvalue
            if p < 0.05:
                print(f"      {field:14s} Kruskal-Wallis across classes: p={p:.4g}")

# ---------------------------------------------------------------- 4. the duplication arm by distance
print(f"\n{'=' * 118}\n=== 4. the value of a second AP-1 copy, by distance (net of a shuffled control at the same place)")
rows = []
for o in OFFSETS:
    if o not in DUP:
        continue
    v = DUP[o].rename("net").to_frame().join(S[["stable_dynamic"]])
    sv, dv = v[v.stable_dynamic == "Stable"].net.dropna(), v[v.stable_dynamic == "Dynamic"].net.dropna()
    pr = [(DUP.at[ds, o], DUP.at[ss, o]) for ds, ss in pair_of.items()
          if ds in DUP.index and ss in DUP.index and not (np.isnan(DUP.at[ds, o]) or np.isnan(DUP.at[ss, o]))]
    row = {"offset_bp": o, "stable": sv.median(), "dynamic": dv.median(),
           "p_unmatched": mannwhitneyu(dv, sv).pvalue if min(len(sv), len(dv)) >= a.min_sites else np.nan}
    if len(pr) >= a.min_sites:
        diff = np.array(pr)[:, 0] - np.array(pr)[:, 1]
        row.update({"n_pairs": len(diff), "paired_difference": float(np.median(diff)),
                    "p_matched": float(wilcoxon(diff).pvalue)})
    rows.append(row)
Dd = pd.DataFrame(rows)
if "p_matched" in Dd:
    Dd["fdr_matched"] = bh(Dd.p_matched.fillna(1).values)
Dd.to_csv(OUT / "duplication_by_distance.csv", index=False)
print(Dd.round(4).to_string(index=False))
print("\noutputs in", OUT)
