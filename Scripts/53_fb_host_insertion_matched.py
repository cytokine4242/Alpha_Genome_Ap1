#!/usr/bin/env python3
"""Does the Stable advantage in the host insertions survive matching on the donor's own accessibility?

The host-insertion result (47_fb_swap_analysis.py) compares classes with the host context fully controlled - every element
goes into the SAME 20 closed windows and is scored against its own dinucleotide-shuffled copy - but NOT with the donor
classes matched: Stable donors come from about twice as accessible loci as Dynamic ones. This asks whether a Stable element
still opens a closed site more than a Dynamic element taken from an equally accessible locus.

Three readings, because 1:1 matching can only use 50 of the 200 Dynamic donors:
  1  unmatched, as reported on the slide
  2  how the lift tracks the donor's own native accessibility at all (Spearman, within and across classes)
  3  matched: each Dynamic donor paired with the nearest unused Stable donor in log10 native accessibility (caliper
     --caliper), paired test; and, using every donor, an analysis of covariance - lift ~ class + log10 native - so the
     class term is the part of the gap that the accessibility difference does not explain
Lift is log2(element in the host) - log2(the same sequence shuffled in the same host), per donor the median over hosts.
Outputs Analysis/fb_swap/tables/host_insertion_accessibility_matched.csv.
  python3 53_fb_host_insertion_matched.py [--width peak] [--caliper 0.1]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--caliper", type=float, default=0.1, help="matching tolerance in log10 native accessibility")
ap.add_argument("--match-on", default="measured", choices=["measured", "predicted"],
                help="measured = the donor site's own fibroblast ATAC CPM/kb, which the model never saw; "
                     "predicted = the model's native call over the element, which partly IS the thing being measured")
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_swap" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def lg(x):
    v = np.log2(np.asarray(x, dtype=np.float64).clip(min=1e-3))
    return pd.Series(v, index=x.index) if isinstance(x, pd.Series) else v


files = sorted(glob.glob(str(RES / "bg_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no bg_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d = d[d.trackset == a.track]
piv = d.pivot_table(index=["donor_id", "donor_class", "host_id", "width"], columns="kind",
                    values=["atac_insert", "donor_native_atac_insert"]).reset_index()
piv.columns = [c[0] if not c[1] else f"{c[0]}__{c[1]}" for c in piv.columns]
piv["lift"] = lg(piv["atac_insert__element"]) - lg(piv["atac_insert__shuffled"])
piv["native"] = piv["donor_native_atac_insert__element"]
if a.match_on == "measured":
    meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
    meas = meas[meas.site_table == "fb"].drop_duplicates("site_id").set_index("site_id").fb_cpm_per_kb
    piv["native"] = piv.donor_id.map(meas)
    gone = piv[piv.native.isna()].drop_duplicates(["donor_id", "donor_class"]).donor_class.value_counts()
    if len(gone):
        print("donors with no measured fibroblast peak, dropped from the matched analysis: "
              + ", ".join(f"{k} {v}" for k, v in gone.items()))
    piv = piv[piv.native.notna()]
print(f"matching on the {a.match_on} accessibility of the donor site")

rows = []
for width, W in piv.groupby("width"):
    per = W.groupby(["donor_id", "donor_class"], as_index=False).agg(lift=("lift", "median"), native=("native", "median"),
                                                                     n_hosts=("lift", "size"))
    per["dyn"] = np.where(per.donor_class == "Stable", "Stable", "Dynamic")
    per["l"] = np.log10(per.native.clip(lower=1e-3))
    S, D = per[per.dyn == "Stable"], per[per.dyn == "Dynamic"]
    print(f"\n{'=' * 96}\n=== width {width}: {len(S)} Stable donors, {len(D)} Dynamic donors, {per.n_hosts.iloc[0]} hosts each")
    print(f"    donor native accessibility: Stable {S.native.median():.2f}  Dynamic {D.native.median():.2f}  "
          f"(ratio {S.native.median() / D.native.median():.2f}x)")

    # 1 ---------------------------------------------------------- unmatched, as on the slide
    p_un = mannwhitneyu(S.lift, D.lift).pvalue
    print(f"\n 1. unmatched   Stable {S.lift.median():+.3f}  Dynamic {D.lift.median():+.3f}  "
          f"gap {S.lift.median() - D.lift.median():+.3f} log2  p={p_un:.3g}")
    print("    by class: " + "  ".join(f"{c} {g.lift.median():+.2f}" for c, g in per.groupby("donor_class")))

    # 2 ---------------------------------------------------------- does the lift track accessibility at all?
    rho_all, p_all = spearmanr(per.l, per.lift)
    rho_s, _ = spearmanr(S.l, S.lift)
    rho_d, _ = spearmanr(D.l, D.lift)
    print(f"\n 2. lift vs the donor's own native accessibility: all donors rho={rho_all:+.3f} (p={p_all:.3g}), "
          f"within Stable {rho_s:+.3f}, within Dynamic {rho_d:+.3f}")

    # 3a --------------------------------------------------------- 1:1 matching on log10 native accessibility
    used, pairs = set(), []
    for r in D.sort_values("l").itertuples():
        free = S[~S.donor_id.isin(used)]
        if free.empty:
            break
        j = (free.l - r.l).abs().idxmin()
        if abs(S.loc[j, "l"] - r.l) <= a.caliper:
            used.add(S.loc[j, "donor_id"])
            pairs.append((r.lift, S.loc[j, "lift"], r.native, S.loc[j, "native"]))
    p_matched = np.nan
    if len(pairs) >= 10:
        P = pd.DataFrame(pairs, columns=["dyn_lift", "sta_lift", "dyn_native", "sta_native"])
        diff = P.dyn_lift - P.sta_lift
        p_matched = wilcoxon(diff).pvalue
        print(f"\n 3. matched (caliper {a.caliper} log10): {len(P)} pairs; native after matching "
              f"Dynamic {P.dyn_native.median():.2f} vs Stable {P.sta_native.median():.2f}")
        print(f"    Stable {P.sta_lift.median():+.3f}  Dynamic {P.dyn_lift.median():+.3f}  "
              f"(difference of medians {P.sta_lift.median() - P.dyn_lift.median():+.3f})")
        print(f"    paired gap: median of the per-pair differences {-diff.median():+.3f}, mean {-diff.mean():+.3f} log2  "
              f"p={p_matched:.3g}; Stable wins {(diff < 0).sum()}/{len(diff)} pairs")
        print(f"    smallest gap these pairs could detect at 80% power: {2.8 * diff.std() / np.sqrt(len(diff)):.3f} log2")
    else:
        print(f"\n 3. matched: only {len(pairs)} pairs within caliper {a.caliper} - too few to test")

    # 3b --------------------------------------------------------- covariance adjustment, every donor kept
    X = np.column_stack([np.ones(len(per)), (per.dyn == "Stable").astype(float), per.l.values])
    beta, *_ = np.linalg.lstsq(X, per.lift.values, rcond=None)
    resid = per.lift.values - X @ beta
    dof = len(per) - X.shape[1]
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * (resid @ resid) / dof)
    from scipy.stats import t as tdist
    p_adj = 2 * tdist.sf(abs(beta[1] / se[1]), dof)
    raw_mean = S.lift.mean() - D.lift.mean()
    print(f"\n 3b. lift ~ class + log10 native, all {len(per)} donors (means, to match the regression):")
    print(f"    raw mean gap {raw_mean:+.3f} log2  ->  Stable term after adjustment {beta[1]:+.3f} +- {se[1]:.3f} "
          f"(p={p_adj:.3g});  accessibility term {beta[2]:+.3f} per log10")
    print(f"    so {100 * (1 - beta[1] / raw_mean):.0f}% of the raw gap is explained by the donors' own accessibility")

    rows.append({"width": width, "n_stable": len(S), "n_dynamic": len(D),
                 "native_stable": S.native.median(), "native_dynamic": D.native.median(),
                 "lift_stable_unmatched": S.lift.median(), "lift_dynamic_unmatched": D.lift.median(),
                 "gap_unmatched": S.lift.median() - D.lift.median(), "p_unmatched": p_un,
                 "gap_unmatched_mean": S.lift.mean() - D.lift.mean(),
                 "rho_lift_vs_native": rho_all, "n_pairs": len(pairs),
                 "gap_matched": (P.sta_lift.median() - P.dyn_lift.median()) if len(pairs) >= 10 else np.nan,
                 "p_matched": p_matched, "gap_adjusted": beta[1], "se_adjusted": se[1], "p_adjusted": p_adj,
                 "match_on": a.match_on})

R = pd.DataFrame(rows)
R.to_csv(OUT / f"host_insertion_accessibility_matched_{a.match_on}.csv", index=False)
print(f"\n{'=' * 96}\nwrote {OUT / f'host_insertion_accessibility_matched_{a.match_on}.csv'}")
