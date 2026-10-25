#!/usr/bin/env python3
"""Where does a fibroblast AP-1 site's single most important 10 bp tile sit - on the AP-1 core, or somewhere else?

49/50 ask how much of a site's tile importance falls near the motif, which is a ratio of magnitudes and therefore
scales with how accessible the site is; 50's motif_share does not survive accessibility matching (p = 0.11). This asks
the rank question instead: of the 60 tiles scanned across the site, which one costs the most to replace? A rank is
invariant to any per-site multiplicative or additive change in signal, so the accessibility gap between the classes
cannot move it on its own - which is the point of running it.

Tiles are 10 bp, laid from -300 to +290 bp around the 11 bp AP-1 core, so each tile's midpoint is pos + 5. Two of the
60 tiles overlap the core (pos -10 and 0, the same MOTIF_TILES as 50), eight more lie inside +-50 bp, and the remaining
50 sit further out. Those counts give the chance levels the figure is read against: 3.3%, 13.3%, 83.3%.

  1  per site: the zone of the top tile, whether a core tile is in the top 3, and the rank of the better core tile
  2  by class - Stable vs Dynamic, and Dynamic split into the age classes as given (Age Up / Age Down / the rest)
  3  the same comparison on Stable sites matched to the Dynamic ones on native accessibility (caliper 0.1 log10)

Outputs Analysis/fb_scan/tables/top_tile_*.csv.
  python3 56_fb_top_tile_rank.py [--res DIR] [--readout atac_170]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2_contingency, fisher_exact, mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170")
ap.add_argument("--caliper", type=float, default=0.1)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_scan" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

MOTIF_TILES = (-10, 0)          # the two tiles that overlap the 11 bp core, as in 50_fb_scan_per_class.py
NEAR_BP = 50                    # the flank zone, by tile midpoint
ZONES = ("AP-1 core", "Flank, within 50 bp", "Beyond 50 bp")


def wilson(k, n, z=1.96):
    """Score interval for a proportion - behaves at the small counts the per-class split produces."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


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
key = ["site_id", "binding_category", "stable_dynamic", "age_class", "ex_class"]
t = t.groupby(key + ["pos"], as_index=False).eff.mean()
# importance is the size of the effect, as in 49/50. The per-site noise floor those scripts subtract is a constant
# within a site, so it cannot change the order of a site's own tiles; ranking is done on the raw magnitude.
t["imp"] = t.eff.abs()
t["mid"] = t.pos + 5
t["zone"] = np.where(t.pos.isin(MOTIF_TILES), ZONES[0],
                     np.where(t.mid.abs() <= NEAR_BP, ZONES[1], ZONES[2]))
n_tiles = int(t.groupby("site_id").pos.nunique().median())
per_zone = t.drop_duplicates("pos").zone.value_counts()
print(f"{len(files)} chunks, {t.site_id.nunique()} sites, {n_tiles} tiles per site, readout {a.readout}")
print("  tiles per zone, i.e. the chance level each percentage is read against: "
      + ", ".join(f"{z} {per_zone.get(z, 0)}/{n_tiles} = {100 * per_zone.get(z, 0) / n_tiles:.1f}%" for z in ZONES))
pd.DataFrame([{"zone": z, "n_tiles": int(per_zone.get(z, 0)), "chance_pct": 100 * per_zone.get(z, 0) / n_tiles}
              for z in ZONES]).to_csv(OUT / "top_tile_chance.csv", index=False)

# ---------------------------------------------------------------- 1. per site
t["rank"] = t.groupby("site_id").imp.rank(ascending=False, method="min")
top = t.sort_values("rank").groupby(key, as_index=False).first()[key + ["pos", "mid", "zone", "imp"]]
top = top.rename(columns={"pos": "top_pos", "mid": "top_mid", "zone": "top_zone", "imp": "top_imp"})
core_rank = t[t.pos.isin(MOTIF_TILES)].groupby("site_id")["rank"].min().rename("core_rank")
top = top.join(core_rank, on="site_id")
top["core_is_top"] = top.top_zone == ZONES[0]
top["core_in_top3"] = top.core_rank <= 3
top["native"] = top.site_id.map(nat)
top["grp4"] = np.where(top.stable_dynamic == "Stable", "Stable",
                       np.where(top.age_class == "Age Down", "AP-1 loss",
                                np.where(top.age_class == "Age Up", "AP-1 gain", "other Dynamic")))
top.to_csv(OUT / "top_tile_per_site.csv", index=False)

# ---------------------------------------------------------------- 2. by class
rows = []
for field, order in (("stable_dynamic", ("Stable", "Dynamic")),
                     ("grp4", ("Stable", "AP-1 gain", "AP-1 loss", "other Dynamic"))):
    for cls in order:
        g = top[top[field] == cls]
        if len(g) < 15:
            continue
        for z in ZONES:
            k = int((g.top_zone == z).sum())
            lo, hi = wilson(k, len(g))
            rows.append({"grouping": field, "class": cls, "n_sites": len(g), "zone": z,
                         "k": k, "pct": 100 * k / len(g), "lo": 100 * lo, "hi": 100 * hi,
                         "pct_core_in_top3": 100 * g.core_in_top3.mean(),
                         "median_core_rank": g.core_rank.median(),
                         "native_median": g.native.median()})
Z = pd.DataFrame(rows)
Z.to_csv(OUT / "top_tile_zone_by_class.csv", index=False)
pd.set_option("display.width", 200)
print("\n=== where the most important tile sits, by class")
print(Z.pivot_table(index=["grouping", "class", "n_sites"], columns="zone", values="pct").round(1).to_string())
print("\n=== an AP-1 core tile in the site's top 3, and the core tile's median rank of 60")
print(Z.drop_duplicates(["grouping", "class"])[["grouping", "class", "n_sites", "pct_core_in_top3",
                                                "median_core_rank", "native_median"]].round(2).to_string(index=False))

# ---------------------------------------------------------------- 3. tests
S = top[top.stable_dynamic == "Stable"]
D = top[top.stable_dynamic == "Dynamic"]
tab = [[int(S.core_is_top.sum()), int((~S.core_is_top).sum())],
       [int(D.core_is_top.sum()), int((~D.core_is_top).sum())]]
orv, p_fisher = fisher_exact(tab)
p_rank = mannwhitneyu(S.core_rank.dropna(), D.core_rank.dropna()).pvalue
p_top3 = fisher_exact([[int(S.core_in_top3.sum()), int((~S.core_in_top3).sum())],
                       [int(D.core_in_top3.sum()), int((~D.core_in_top3).sum())]])[1]
ct = pd.crosstab(top.grp4, top.top_zone)
p_chi4 = chi2_contingency(ct)[1]
tests = [{"comparison": "top tile is the AP-1 core, Stable vs Dynamic", "statistic": "Fisher odds ratio",
          "value": orv, "p": p_fisher, "n": len(top), "matched": False},
         {"comparison": "AP-1 core tile in the top 3, Stable vs Dynamic", "statistic": "Fisher",
          "value": np.nan, "p": p_top3, "n": len(top), "matched": False},
         {"comparison": "rank of the AP-1 core tile, Stable vs Dynamic", "statistic": "Mann-Whitney",
          "value": S.core_rank.median() - D.core_rank.median(), "p": p_rank, "n": len(top), "matched": False},
         {"comparison": "zone of the top tile across the four classes", "statistic": "chi-square",
          "value": np.nan, "p": p_chi4, "n": len(top), "matched": False}]

# ---------------------------------------------------------------- 4. the same thing, accessibility-matched
# A rank should already be immune to the accessibility gap; this is the check that it is, run exactly as 50 runs it.
dm = D.dropna(subset=["native"]).copy()
sm = S.dropna(subset=["native"]).copy()
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
if len(pairs) > 30:
    PD = pd.DataFrame([{"core_is_top": x.core_is_top, "core_rank": x.core_rank, "native": x.native} for x, _ in pairs])
    PS = pd.DataFrame([{"core_is_top": y.core_is_top, "core_rank": y.core_rank, "native": y.native} for _, y in pairs])
    print(f"    native accessibility after matching: Dynamic {PD.native.median():.2f}  Stable {PS.native.median():.2f}")
    print(f"    top tile is the AP-1 core: Dynamic {100 * PD.core_is_top.mean():.1f}%  "
          f"Stable {100 * PS.core_is_top.mean():.1f}%")
    # discordant pairs only, which is what a paired binary comparison actually rests on (McNemar, exact)
    b = int((PS.core_is_top & ~PD.core_is_top).sum())
    c = int((~PS.core_is_top & PD.core_is_top).sum())
    p_mc = binomtest(b, b + c, 0.5).pvalue if (b + c) > 0 else np.nan
    dr = (PS.core_rank - PD.core_rank).dropna()
    p_wr = wilcoxon(dr).pvalue if len(dr) > 10 and dr.abs().sum() > 0 else np.nan
    print(f"    McNemar on the {b + c} discordant pairs ({b} Stable-only, {c} Dynamic-only): p = {p_mc:.3g}")
    print(f"    paired difference in core-tile rank: median {dr.median():.1f}, Wilcoxon p = {p_wr:.3g}")
    tests += [{"comparison": "top tile is the AP-1 core, matched", "statistic": f"McNemar ({b + c} discordant)",
               "value": 100 * (PS.core_is_top.mean() - PD.core_is_top.mean()), "p": p_mc,
               "n": len(pairs), "matched": True},
              {"comparison": "rank of the AP-1 core tile, matched", "statistic": "Wilcoxon",
               "value": dr.median(), "p": p_wr, "n": len(pairs), "matched": True}]
    pd.DataFrame([{"n_pairs": len(pairs), "pct_core_top_dynamic": 100 * PD.core_is_top.mean(),
                   "pct_core_top_stable": 100 * PS.core_is_top.mean(),
                   "native_median_dynamic": PD.native.median(), "native_median_stable": PS.native.median(),
                   "p_mcnemar": p_mc, "median_rank_difference": dr.median(), "p_wilcoxon_rank": p_wr}]
                 ).to_csv(OUT / "top_tile_matched.csv", index=False)

T = pd.DataFrame(tests)
T.to_csv(OUT / "top_tile_tests.csv", index=False)
print("\n=== tests")
print(T.round(4).to_string(index=False))
print("\nwrote", OUT / "top_tile_zone_by_class.csv", "| top_tile_per_site.csv | top_tile_tests.csv")
