#!/usr/bin/env python3
"""Is AP-1 special, and are AP-1-loss sites selectively or generally disengaged?
Analysis of Katana/results_ft/pshuf_fb_chunk*.csv (19_partner_shuffle.py, fine-tuned 1 bp fibroblast head).

Each site contributes its own AP-1 motif shuffled, each footprinted partner shuffled, and for every partner a
footprint-free window of the SAME width at the mirrored distance, also shuffled. The control is the point: the net
effect of a motif is (partner - control), which removes "shuffling that many bases at that distance would matter
anyway". Comparisons are within a site, so the accessibility difference between classes cannot drive them.

  1  net effect by family: which motifs matter beyond their own background, ranked
  2  the site's own AP-1 motif against its partners, by class
  3  the decisive test: at AP-1-loss sites, is the deficit specific to AP-1 or shared by every family?
     Per site: own AP-1 effect, and the mean net effect of its partners. If loss sites are weak only on AP-1 the first
     drops and the second does not; if they are weak elements generally, both drop together.
Outputs Analysis/fb_partner_shuffle/tables.   python3 52_fb_partner_shuffle_analysis.py [--res DIR]
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
ap.add_argument("--prefix", default="pshuf_fb")
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--min-sites", type=int, default=20)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_partner_shuffle" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


files = sorted(glob.glob(str(RES / f"{a.prefix}_chunk*.csv")))
if not files:
    raise SystemExit(f"no {a.prefix}_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[a.readout]
pc = 0.05 * nat.median()
e = d[d.kind != "native"].copy()
e["eff"] = np.log2((e[a.readout] + pc) / (e.site_id.map(nat) + pc))
per = e.groupby(["site_id", "stable_dynamic", "age_class", "kind", "family", "edit_width"], as_index=False).eff.mean()
per["grp"] = np.where(per.stable_dynamic == "Stable", "Stable",
                      np.where(per.age_class == "Age Down", "AP-1 loss",
                               np.where(per.age_class == "Age Up", "AP-1 gain", "other Dynamic")))
GR = ("Stable", "AP-1 loss", "AP-1 gain", "other Dynamic")
print(f"{len(files)} chunks, {per.site_id.nunique()} sites, readout {a.readout}")
print("  sites per class: " + ", ".join(f"{g} {per[per.grp == g].site_id.nunique()}" for g in GR))

# ---------------------------------------------------------------- net effect per site and family
pr = per[per.kind == "partner"][["site_id", "grp", "family", "edit_width", "eff"]].rename(columns={"eff": "partner"})
ct = per[per.kind == "ctrl"][["site_id", "family", "eff"]].rename(columns={"eff": "control"})
net = pr.merge(ct, on=["site_id", "family"], how="inner")
net["net"] = net.partner - net.control
own = per[per.kind == "own_ap1"][["site_id", "grp", "eff"]].rename(columns={"eff": "own_ap1"})
net.to_csv(OUT / "net_effect_per_site_family.csv", index=False)

rows = []
for fam, g in net.groupby("family"):
    if g.site_id.nunique() < a.min_sites:
        continue
    rows.append({"family": fam, "n_sites": g.site_id.nunique(), "median_width": g.edit_width.median(),
                 "partner": g.partner.median(), "control": g.control.median(), "net": g.net.median(),
                 "p_vs_control": wilcoxon(g.net).pvalue if (g.net != 0).any() else np.nan})
F = pd.DataFrame(rows)
F["fdr"] = bh(F.p_vs_control.fillna(1).values)
F = F.sort_values("net")
F.to_csv(OUT / "net_effect_by_family.csv", index=False)
print("\n=== 1. net effect of shuffling each motif, beyond a matched footprint-free window")
print(F.round(4).to_string(index=False))
print(f"   own AP-1 motif (no matched control is possible at distance 0): median {own.own_ap1.median():.4f}, n={len(own)}")

# ---------------------------------------------------------------- own AP-1 vs partners, by class
print("\n=== 2. the site's own AP-1 motif, by class")
for g in GR:
    w = own[own.grp == g]
    if len(w) >= a.min_sites:
        print(f"   {g:14s} n={len(w):4d}  own AP-1 {w.own_ap1.median():+.4f}")
for x, y in (("Stable", "AP-1 loss"), ("Stable", "AP-1 gain"), ("AP-1 loss", "AP-1 gain")):
    p1, p2 = own[own.grp == x].own_ap1, own[own.grp == y].own_ap1
    if min(len(p1), len(p2)) >= a.min_sites:
        print(f"     {x} vs {y}: p={mannwhitneyu(p1, p2).pvalue:.4f}")

# ---------------------------------------------------------------- the decisive test
pm = net.groupby(["site_id", "grp"]).net.mean().rename("partners_mean_net").reset_index()
J = own.merge(pm[["site_id", "partners_mean_net"]], on="site_id", how="inner")
J["ap1_minus_partners"] = J.own_ap1 - J.partners_mean_net
J.to_csv(OUT / "own_vs_partners_per_site.csv", index=False)
print("\n=== 3. is the AP-1-loss deficit specific to AP-1, or shared by every motif?")
print(f"   {'class':14s} {'n':>5s} {'own AP-1':>10s} {'partners (net)':>15s} {'difference':>11s}")
for g in GR:
    w = J[J.grp == g]
    if len(w) >= a.min_sites:
        print(f"   {g:14s} {len(w):5d} {w.own_ap1.median():10.4f} {w.partners_mean_net.median():15.4f} {w.ap1_minus_partners.median():11.4f}")
rows = []
for x, y in (("Stable", "AP-1 loss"), ("Stable", "AP-1 gain"), ("AP-1 loss", "AP-1 gain")):
    A, B = J[J.grp == x], J[J.grp == y]
    if min(len(A), len(B)) < a.min_sites:
        continue
    rows.append({"comparison": f"{x} vs {y}",
                 "p_own_ap1": mannwhitneyu(A.own_ap1, B.own_ap1).pvalue,
                 "p_partners": mannwhitneyu(A.partners_mean_net, B.partners_mean_net).pvalue,
                 "p_difference": mannwhitneyu(A.ap1_minus_partners, B.ap1_minus_partners).pvalue,
                 "own_ap1_gap": B.own_ap1.median() - A.own_ap1.median(),
                 "partners_gap": B.partners_mean_net.median() - A.partners_mean_net.median()})
T = pd.DataFrame(rows)
T.to_csv(OUT / "specificity_tests.csv", index=False)
print(T.round(4).to_string(index=False))
print("\n   Reading: a class with a large own_ap1_gap and a small partners_gap is selectively disengaged from AP-1;")
print("   if both gaps move together the sites are simply less sequence-dependent everywhere.")

# per-family version of the same question, so one family cannot drive the mean
rows = []
for fam, g in net.groupby("family"):
    A, B = g[g.grp == "Stable"], g[g.grp == "AP-1 loss"]
    if min(A.site_id.nunique(), B.site_id.nunique()) < a.min_sites:
        continue
    rows.append({"family": fam, "n_stable": A.site_id.nunique(), "n_loss": B.site_id.nunique(),
                 "net_stable": A.net.median(), "net_loss": B.net.median(),
                 "p": mannwhitneyu(A.net, B.net).pvalue})
P = pd.DataFrame(rows)
if not P.empty:
    P["fdr"] = bh(P.p.values)
    P.to_csv(OUT / "per_family_stable_vs_loss.csv", index=False)
    print("\n=== 4. per family: Stable vs AP-1-loss net effect (does any single family carry the difference?)")
    print(P.sort_values("p").round(4).to_string(index=False))
# ---------------------------------------------------------------- 5. how AP-1-centric is each class?
# The share of a site's measured motif dependence that sits on its own AP-1 motif rather than on its partners. Every
# term is an AlphaGenome effect; TOBIAS only says where the partner spans are. Because partner detection needs signal,
# Stable sites carry more detected partners, so the test is repeated inside each stratum of partner count.
J["ap1_share"] = J.own_ap1.abs() / (J.own_ap1.abs() + J.partners_mean_net.abs())
J["n_partners"] = J.site_id.map(net.groupby("site_id").family.nunique())
J["sd"] = np.where(J.grp == "Stable", "Stable", "Dynamic")
J.to_csv(OUT / "own_vs_partners_per_site.csv", index=False)
print("\n=== 5. AP-1 share = |own AP-1| / (|own AP-1| + |partners, net|): how AP-1-centric is the site?")
for g in GR:
    w = J[J.grp == g].ap1_share.dropna()
    if len(w) >= a.min_sites:
        print(f"   {g:14s} n={len(w):5d}  {w.median():.3f}  (IQR {w.quantile(.25):.3f}-{w.quantile(.75):.3f})")
S_, D_ = J[J.sd == "Stable"].ap1_share.dropna(), J[J.sd == "Dynamic"].ap1_share.dropna()
print(f"   pooled Stable {S_.median():.3f} vs Dynamic {D_.median():.3f}  p={mannwhitneyu(S_, D_).pvalue:.3g}")
print(f"   partners detected per site: Stable {J[J.sd == 'Stable'].n_partners.mean():.2f} vs "
      f"Dynamic {J[J.sd == 'Dynamic'].n_partners.mean():.2f} - so the same test within each stratum:")
rows = []
for n, g in J.groupby("n_partners"):
    x, y = g[g.sd == "Stable"].ap1_share.dropna(), g[g.sd == "Dynamic"].ap1_share.dropna()
    if min(len(x), len(y)) >= 20:
        rows.append({"n_partners": int(n), "n_stable": len(x), "n_dynamic": len(y), "stable": x.median(),
                     "dynamic": y.median(), "difference": y.median() - x.median(), "p": mannwhitneyu(x, y).pvalue})
A = pd.DataFrame(rows)
if not A.empty:
    A.to_csv(OUT / "ap1_share_by_partner_count.csv", index=False)
    print(A.round(4).to_string(index=False))

print("\noutputs in", OUT)
