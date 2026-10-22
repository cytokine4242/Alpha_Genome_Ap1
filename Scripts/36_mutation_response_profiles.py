#!/usr/bin/env python3
"""Do Dynamic and Stable AP-1 sites respond differently to different kinds of mutation? Cardiomyocyte sites, heart ATAC track. No TOBIAS footprint scores are used;
the classes are the TOBIAS binding labels, and Dynamic sites are compared with Stable sites matched on predicted accessibility, peak width and TSS distance
(pairs from 33_dynamic_vs_stable_deep_dive.py, tables/matched_pair_ids.csv).
 A  response SHAPE: each site's 21 mutation effects divided by its overall size (L2 norm), so a class difference in the mix of mutation types is not hidden by a
    difference in how strongly the site responds overall. Per mutation type (paired Wilcoxon) and all at once (cross-validated classifier, folds grouped by peak).
 B  ratios with a plain meaning, at sites that respond (denominator below -0.5): motif share of the local dependence, core share of the motif, composition vs order, reach
 C  non-responders: sites where destroying the motif changes accessibility by less than 10%
 D  tile scan (subset of 1,470 sites, 34_make_dynamics_subsets.py): does a site depend on its own AP-1 motif or on a neighbouring element ("partner") within ±300 bp?
Outputs Analysis/dynamic_profiles/{tables,figures}.   python3 36_mutation_response_profiles.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "dynamic_profiles"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(36)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_med(d, n=1000):
    d = np.asarray(d, float)
    return np.percentile([np.median(rng.choice(d, len(d))) for _ in range(n)], [2.5, 97.5])


FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
P = pd.read_csv(ROOT / "Analysis/dynamic/tables/matched_pair_ids.csv")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
prof = pd.read_csv(ROOT / "Analysis/dynamic/tables/mutation_type_profile_dynamic_vs_stable.csv")
TYPES = prof.feature.tolist()                       # the 21 mutation types of the earlier profile, same order and labels
SHORT = dict(zip(prof.feature, prof.mutation))
GROUP = dict(zip(prof.feature, prof.group))
peak_of = sl.peak_coord if "peak_coord" in sl else pd.Series(sl.index, index=sl.index)

# ============================================================ A. response shape
V = FD[TYPES].dropna()
S = V.div(np.sqrt((V ** 2).sum(1)), axis=0)         # unit-length profile: the mix of mutation types, not the overall size
pr = P[P.dyn.isin(S.index) & P.stable.isin(S.index)]
rows = []
for c in TYPES:
    d = S.loc[pr.dyn, c].values - S.loc[pr.stable, c].values
    lo, hi = boot_med(d, 500)
    rows.append({"group": GROUP[c], "mutation": SHORT[c], "feature": c, "n_pairs": len(d), "share_stable": S.loc[pr.stable, c].median(), "share_dynamic": S.loc[pr.dyn, c].median(),
                 "paired_median_diff": np.median(d), "ci_lo": lo, "ci_hi": hi, "p_wilcoxon": wilcoxon(d).pvalue})
A = pd.DataFrame(rows)
A["fdr"] = bh(A.p_wilcoxon)
A.to_csv(OUT / "tables" / "a_response_shape_by_mutation_type.csv", index=False)


def cv_auc(X, y, groups, n_rep=5):
    aucs = []
    ug = np.unique(groups)
    for r in range(n_rep):                          # each repeat assigns peaks to 5 folds at random (sites of one peak stay together)
        fold_of = dict(zip(ug, np.random.default_rng(r).integers(0, 5, len(ug))))
        f = np.array([fold_of[g] for g in groups])
        pred = np.zeros(len(y))
        for k in range(5):
            tr, te = f != k, f == k
            m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.5)).fit(X[tr], y[tr])
            pred[te] = m.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y, pred))
    return float(np.mean(aucs)), float(np.std(aucs))


ids = np.r_[pr.dyn.values, pr.stable.values]
y = np.r_[np.ones(len(pr)), np.zeros(len(pr))].astype(int)
grp = peak_of.reindex(ids).fillna(pd.Series(ids, index=ids)).values
crows = []
for lab, X in (("shape only (each site's profile scaled to unit length)", S.loc[ids].values), ("size and shape (raw effects)", V.loc[ids].values),
               ("overall size only (length of the profile)", np.sqrt((V.loc[ids] ** 2).sum(1)).values[:, None])):
    auc, sd = cv_auc(X, y, grp)
    null = []
    for k in range(100):
        null.append(cv_auc(X, np.random.default_rng(1000 + k).permutation(y), grp, n_rep=1)[0])
    crows.append({"features": lab, "n_sites": len(y), "auc": auc, "auc_sd": sd, "null_auc_95th": float(np.percentile(null, 95)), "p_perm": (np.sum(np.array(null) >= auc) + 1) / (len(null) + 1)})
C = pd.DataFrame(crows)
C.to_csv(OUT / "tables" / "a_classifier_matched_pairs.csv", index=False)

# ============================================================ B. ratios with a plain meaning (responding sites only)
RAT = {
    "motif share of the local dependence: replace motif / centred 680 bp shuffle (read over motif)": ("replace motif, read over motif", "centred shuffle of 680 bp, read over motif"),
    "core share of the motif: replace 7 bp core / replace 11 bp motif (read over motif)": ("replace core, read over motif", "replace motif, read over motif"),
    "composition kept: dinucleotide-shuffle / replace, 170 bp (read over 170 bp)": ("dinucleotide-shuffle 170 bp, read over 170 bp", "replace 170 bp, read over 170 bp"),
    "surroundings share: replace the 159 bp around the motif / replace 170 bp (read over 170 bp)": ("replace surroundings (motif kept), read over 170 bp", "replace 170 bp, read over 170 bp"),
    "reach: centred 41 bp shuffle / centred 680 bp shuffle (read over motif)": ("centred shuffle of 41 bp, read over motif", "centred shuffle of 680 bp, read over motif"),
}
rows = []
for lab, (num, den) in RAT.items():
    ok = (FD[den] < -0.5)
    r = (FD[num] / FD[den]).where(ok)
    q = pr[r.reindex(pr.dyn).notna().values & r.reindex(pr.stable).notna().values]
    d = r.loc[q.dyn].values - r.loc[q.stable].values
    lo, hi = boot_med(d, 500) if len(d) > 20 else (np.nan, np.nan)
    rows.append({"ratio": lab, "n_pairs_both_responding": len(d), "median_stable": np.median(r.loc[q.stable]), "median_dynamic": np.median(r.loc[q.dyn]),
                 "paired_median_diff": np.median(d), "ci_lo": lo, "ci_hi": hi, "p_wilcoxon": wilcoxon(d).pvalue if len(d) > 20 else np.nan})
B = pd.DataFrame(rows)
B["fdr"] = bh(B.p_wilcoxon.fillna(1))
B.to_csv(OUT / "tables" / "b_ratios_matched_pairs.csv", index=False)

# ============================================================ C. non-responders
NR = "replace motif, read over motif"
nr = FD[NR] > np.log2(0.9)                          # destroying the motif changes accessibility by less than 10%
both = pr[pr.dyn.isin(FD.index) & pr.stable.isin(FD.index)]
a, b = nr.loc[both.dyn].values, nr.loc[both.stable].values
n10, n01 = int((a & ~b).sum()), int((~a & b).sum())                       # discordant pairs (McNemar, exact)
from scipy.stats import binomtest
dyn_all = sl.index[sl.stable_dynamic == "Dynamic"].intersection(FD.index)
st_all = sl.index[sl.stable_dynamic == "Stable"].intersection(FD.index)
odds, p_raw = fisher_exact([[nr.loc[dyn_all].sum(), (~nr.loc[dyn_all]).sum()], [nr.loc[st_all].sum(), (~nr.loc[st_all]).sum()]])
Cn = pd.DataFrame([{"definition": "replace the 11 bp motif changes accessibility over the motif by less than 10%", "share_dynamic_all": nr.loc[dyn_all].mean(), "share_stable_all": nr.loc[st_all].mean(),
                    "p_fisher_all": p_raw, "n_pairs": len(both), "share_dynamic_matched": a.mean(), "share_stable_matched": b.mean(), "discordant_dyn_only": n10, "discordant_stable_only": n01,
                    "p_mcnemar_matched": binomtest(n10, n10 + n01).pvalue if n10 + n01 else np.nan}])
Cn.to_csv(OUT / "tables" / "c_non_responders.csv", index=False)

# ============================================================ D. tile scan: own motif vs neighbouring element
GROUPS = [("dynamic", "Dynamic", ORANGE), ("matched_stable", "Stable, context-matched", BLUE), ("high_stable", "Stable, high accessibility", AQUA)]
fs = sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))
D = None
if fs:
    sc = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    sc = sc[sc.trackset == "CM"]
    n_pass = sc.groupby("site_id").size()
    sc = sc[sc.site_id.isin(n_pass.index[n_pass == n_pass.max()])]              # complete sites only
    MOTIF_TILES = {-10, 0}                                                      # tiles -10..-1 and 0..9 cover the 11 bp motif (-5..+5)
    per, thr_rows = [], []
    for rd, col in (("170 bp", "atac_170"), ("motif", "atac_motif")):
        nat = sc[sc.kind == "native"].set_index("site_id")[col]
        pc = 0.05 * nat.median()
        t = sc[sc.kind == "tile"].copy()
        t["eff"] = np.log2((t[col] + pc) / (t.site_id.map(nat) + pc))
        tile = t.groupby(["site_id", "group", "pos"]).eff.mean().reset_index()   # mean of the 3 replicate tiles
        out_t = tile[(tile.pos <= -20) | (tile.pos >= 10)]
        thr = float(np.percentile(tile[tile.pos.abs() >= 250].eff, 1))          # 'strong' = beyond the 1st percentile of tiles 250-300 bp away
        thr_rows.append({"readout": rd, "strong_threshold_log2": thr, "replicate_sd_median": t.groupby(["site_id", "pos"]).eff.std().median()})
        for (sid, g), x in tile.groupby(["site_id", "group"]):
            own = x[x.pos.isin(MOTIF_TILES)].eff.min()
            o = out_t[out_t.site_id == sid]
            j = o.eff.idxmin()
            per.append({"readout": rd, "site_id": sid, "group": g, "own_motif_tile": own, "strongest_partner_tile": o.eff.min(), "partner_pos": int(o.loc[j, "pos"]) + 5,
                        "partner_beyond_85bp": abs(int(o.loc[j, "pos"]) + 5) > 85, "n_strong_partner_tiles": int((o.eff < thr).sum()), "partner_dominates": bool(o.eff.min() < own)})
    DA = pd.DataFrame(per)
    DA.to_csv(OUT / "tables" / "d_tile_scan_per_site.csv", index=False)
    D = DA[DA.readout == "170 bp"]
    rows = []
    for (rd, g), x in DA.groupby(["readout", "group"]):
        rows.append({"readout": rd, "group": dict((a, b) for a, b, _ in GROUPS)[g], "n_sites": len(x), "share_partner_beyond_85bp": x.partner_beyond_85bp.mean(), "own_motif_tile_median": x.own_motif_tile.median(), "strongest_partner_median": x.strongest_partner_tile.median(),
                     "share_partner_dominates": x.partner_dominates.mean(), "share_with_strong_partner": (x.n_strong_partner_tiles > 0).mean(),
                     "median_strong_partner_tiles": x.n_strong_partner_tiles.median()})
    Dsum = pd.DataFrame(rows)
    tests = []
    for (rd, a_, b_) in [(rd, a_, b_) for rd in ("170 bp", "motif") for a_, b_ in (("dynamic", "matched_stable"), ("dynamic", "high_stable"))]:
        xa, xb = DA[(DA.readout == rd) & (DA.group == a_)], DA[(DA.readout == rd) & (DA.group == b_)]
        for col in ("own_motif_tile", "strongest_partner_tile", "n_strong_partner_tiles"):
            tests.append({"readout": rd, "comparison": f"{a_} vs {b_}", "measure": col, "median_a": xa[col].median(), "median_b": xb[col].median(), "p_mannwhitney": mannwhitneyu(xa[col], xb[col]).pvalue})
        tab = [[xa.partner_dominates.sum(), (~xa.partner_dominates).sum()], [xb.partner_dominates.sum(), (~xb.partner_dominates).sum()]]
        tests.append({"readout": rd, "comparison": f"{a_} vs {b_}", "measure": "partner_dominates (share)", "median_a": xa.partner_dominates.mean(), "median_b": xb.partner_dominates.mean(), "p_mannwhitney": fisher_exact(tab)[1]})
    Dt = pd.DataFrame(tests)
    Dt["fdr"] = bh(Dt.p_mannwhitney)
    Dsum.to_csv(OUT / "tables" / "d_tile_scan_summary.csv", index=False)
    Dt.to_csv(OUT / "tables" / "d_tile_scan_tests.csv", index=False)
    pd.DataFrame(thr_rows).assign(n_complete_sites=D.site_id.nunique()).to_csv(OUT / "tables" / "d_tile_scan_thresholds.csv", index=False)

# ============================================================ figures
# z1: response shape per mutation type (paired difference in share, with CI)
fig, ax = plt.subplots(figsize=(8.2, 6.6))
A_ = A.iloc[::-1].reset_index(drop=True)
ys = np.arange(len(A_))
sig = A_.fdr < 0.05
ax.errorbar(A_.paired_median_diff, ys, xerr=[A_.paired_median_diff - A_.ci_lo, A_.ci_hi - A_.paired_median_diff], fmt="none", ecolor=GREY, elinewidth=1.4)
ax.scatter(A_.paired_median_diff[~sig], ys[~sig], s=36, color=GREY, zorder=3, label="not significant (FDR ≥ 0.05)")
ax.scatter(A_.paired_median_diff[sig], ys[sig], s=46, color=ORANGE, zorder=4, label="FDR < 0.05")
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(ys, [f"{m}  ·  {g.split(' (')[0].lower()}" for m, g in zip(A_.mutation, A_.group)], fontsize=8)
ax.set_xlabel("Dynamic minus matched Stable: share of the response\ntaken by this mutation type (> 0 = relatively less at Dynamic sites)")
ax.set_title("Response shape (each site's 21 effects scaled to unit length)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "z1_response_shape.png", bbox_inches="tight")
plt.close(fig)

# z2: classifier AUCs (shape, size, both) against the permutation null
fig, ax = plt.subplots(figsize=(7.2, 2.8))
for i, r in C.iterrows():
    ax.barh(i, r.auc - 0.5, left=0.5, height=0.55, color=[ORANGE, BLUE, GREY][i], alpha=0.85)
    ax.errorbar(r.auc, i, xerr=r.auc_sd, color=INK, capsize=3, lw=1)
    ax.plot([r.null_auc_95th] * 2, [i - 0.33, i + 0.33], color=INK, lw=1.4, ls=":")
    ax.text(max(r.auc + r.auc_sd, r.null_auc_95th) + 0.005, i, f"AUC {r.auc:.2f}, permutation p = {r.p_perm:.2g}", va="center", fontsize=8)
ax.set_yticks(range(len(C)), C.features, fontsize=8.5)
ax.set_xlim(0.45, 0.75)
ax.axvline(0.5, color=MUTED, lw=0.9)
ax.set_xlabel("Cross-validated AUC, Dynamic vs matched Stable (dotted line = 95th percentile of shuffled labels)")
ax.set_title("Can the mutation responses tell matched Dynamic and Stable sites apart?", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "z2_classifier.png", bbox_inches="tight")
plt.close(fig)

if D is not None:
    # z3: own motif vs strongest neighbouring element, per site, and share where the neighbour dominates
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    for g, lab, c in GROUPS[:2]:
        x = D[D.group == g]
        ax.scatter(x.own_motif_tile, x.strongest_partner_tile, s=16, color=c, alpha=0.7, label=f"{lab} (n = {len(x)})", edgecolors="none")
    lim = [min(D.own_motif_tile.min(), D.strongest_partner_tile.min()) - 0.1, 0.4]
    ax.plot(lim, lim, color=MUTED, lw=0.9, ls="--")
    ax.text(lim[0] + 0.1, lim[0] + 0.35, "neighbour matters more\n(below the line)", fontsize=8, color=MUTED)
    ax.set_xlabel("Effect of replacing the AP-1 motif tile (log2FC, read over 170 bp)")
    ax.set_ylabel("Effect of the strongest other 10 bp tile within ±300 bp")
    ax.set_title("Own motif vs strongest neighbouring element, per site", fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    for k, (g, lab, c) in enumerate(GROUPS):
        x = D[D.group == g].partner_dominates.values
        lo, hi = np.percentile([rng.choice(x, len(x)).mean() for _ in range(1000)], [2.5, 97.5])
        ax.bar(k, x.mean(), color=c, alpha=0.85, width=0.6)
        ax.errorbar(k, x.mean(), yerr=[[x.mean() - lo], [hi - x.mean()]], color=INK, capsize=3, lw=1)
        ax.text(k, hi + 0.02, f"{x.mean():.0%}\nn = {len(x)}", ha="center", fontsize=8)
    ax.set_xticks(range(3), [g[1] for g in GROUPS], fontsize=8, rotation=12)
    ax.set_ylabel("Share of sites where a neighbouring tile\nmatters more than the AP-1 motif")
    ax.set_ylim(0, 1)
    ax.set_title("How often does a neighbour dominate?", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "z3_own_motif_vs_partner.png", bbox_inches="tight")
    plt.close(fig)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 95)
print("== A shape\n", A[["group", "mutation", "share_stable", "share_dynamic", "paired_median_diff", "ci_lo", "ci_hi", "p_wilcoxon", "fdr"]].round(3).to_string())
print("== A classifier\n", C.round(3).to_string())
print("== B ratios\n", B.round(3).to_string())
print("== C non-responders\n", Cn.round(3).T.to_string())
if D is not None:
    print("== D summary\n", Dsum.round(3).to_string())
    print("== D tests\n", Dt.round(4).to_string())
    print(pd.read_csv(OUT / "tables" / "d_tile_scan_thresholds.csv").round(3).to_string())
