#!/usr/bin/env python3
"""What makes Dynamic sites different from Stable sites? Analyses of the existing results (cardiomyocyte sites, heart ATAC track unless stated):
 A  how well can model-derived and measured features tell Dynamic from Stable (cross-validated AUC, grouped by peak)?
 B  condition-specific accessibility: measured accessibility and the model's over-prediction in the conditions where a Dynamic site is bound vs unbound
 C  footprint strength in the bound conditions
 D  Stable and Dynamic sites in the SAME peak (paired)
 E  position of the site within its peak
 F  does a site's class carry over between cell types (same sequence, cardiomyocyte vs fibroblast)?
 G  matched pairs: each Dynamic site against the Stable site closest in accessibility, peak width and TSS distance
Outputs Analysis/dynamic/{tables,figures}.   python3 33_dynamic_vs_stable_deep_dive.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "dynamic"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(33)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id")
fb_sites = pd.read_csv(ROOT / "Data/derived/fb_fos_sites.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id")
bycond = pd.read_csv(ROOT / "Katana/data/measured_by_condition.csv")
bycond = bycond[bycond.site_table == "cm"].set_index("site_id")
df = sl.join(sites[["motif_mid", "peak_start", "peak_end", "in_YS", "in_MAS", "in_MAE"]]).join(meas[["cm_cpm_per_kb"]]).join(bycond[["cm_ys", "cm_mas", "cm_mae"]])
df["dyn"] = (df.stable_dynamic == "Dynamic").astype(int)
df["l_peak"], df["l_170"], df["l_motif"] = [np.log10(df[c].clip(lower=1e-3)) for c in ("native_peak", "native_170", "native_motif")]
df["l_pw"], df["l_tss"] = np.log10(df.peak_width), np.log10(df.abs_tss + 1)
df["l_meas"] = np.log10(df.cm_cpm_per_kb.clip(lower=1e-3))
df["pos_in_peak"] = (df.motif_mid - (df.peak_start + df.peak_end) / 2).abs() / (df.peak_width / 2)
df["resid"] = df.native_peak.rank(pct=True) - df.cm_cpm_per_kb.rank(pct=True)

# ---------------------------------------------------------------- A. classifier
FEATS = {
    "native predicted accessibility only (motif, 170 bp, peak)": ["l_motif", "l_170", "l_peak"],
    "peak geometry and position (width, TSS distance, motifs in peak, position in peak)": ["l_pw", "l_tss", "motifs_in_peak", "pos_in_peak"],
    "model effects only (motif, 170 bp, peak shuffle)": ["eff_motif", "eff_170", "eff_peak"],
    "context: native accessibility + geometry": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak"],
    "context + model effects": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak", "eff_motif", "eff_170", "eff_peak"],
    "measured accessibility only": ["l_meas"],
    "context + effects + measured accessibility + prediction error": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak", "eff_motif", "eff_170", "eff_peak", "l_meas", "resid"],
    "footprint score only (partly defines the classes)": ["fp_score"],
}
rows = []
d_ = df.dropna(subset=sum(FEATS.values(), [])).copy()
groups = pd.factorize(d_.peak_coord)[0]
for name, cols in FEATS.items():
    X, y = d_[cols].values, d_.dyn.values
    aucs = {"logistic": [], "boosted trees": []}
    for tr, te in GroupKFold(5).split(X, y, groups):
        lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X[tr], y[tr])
        gb = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.06, max_iter=200, random_state=0).fit(X[tr], y[tr])
        aucs["logistic"].append(roc_auc_score(y[te], lr.predict_proba(X[te])[:, 1]))
        aucs["boosted trees"].append(roc_auc_score(y[te], gb.predict_proba(X[te])[:, 1]))
    rows.append({"features": name, "n_features": len(cols), "auc_logistic": np.mean(aucs["logistic"]), "auc_logistic_sd": np.std(aucs["logistic"]),
                 "auc_trees": np.mean(aucs["boosted trees"]), "auc_trees_sd": np.std(aucs["boosted trees"]), "n_sites": len(d_)})
A = pd.DataFrame(rows)
A.to_csv(OUT / "tables" / "classifier_auc.csv", index=False)
cols = FEATS["context + effects + measured accessibility + prediction error"]
gb = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.06, max_iter=200, random_state=0)
pi = []
for tr, te in GroupKFold(5).split(d_[cols].values, d_.dyn.values, groups):
    gb.fit(d_[cols].values[tr], d_.dyn.values[tr])
    r = permutation_importance(gb, d_[cols].values[te], d_.dyn.values[te], scoring="roc_auc", n_repeats=5, random_state=0)
    pi.append(r.importances_mean)
PI = pd.DataFrame({"feature": cols, "auc_drop": np.mean(pi, 0), "sd": np.std(pi, 0)}).sort_values("auc_drop", ascending=False)
PI.to_csv(OUT / "tables" / "classifier_permutation_importance.csv", index=False)

# ---------------------------------------------------------------- B. condition-specific accessibility and prediction error
cond = {"YS": "cm_ys", "MAS": "cm_mas", "MAE": "cm_mae"}
inb = {"YS": "in_YS", "MAS": "in_MAS", "MAE": "in_MAE"}
mm = df.dropna(subset=list(cond.values()) + ["native_peak"]).copy()
for k, c in cond.items():
    mm[f"m_{k}"] = np.log2(mm[c].clip(lower=1e-3))
    mm[f"p_{k}"] = mm[c].rank(pct=True)
    mm[f"r_{k}"] = mm.native_peak.rank(pct=True) - mm[f"p_{k}"]
pat_rows = []
dyn = mm[mm.dyn == 1].copy()
bound_m, unb_m, bound_r, unb_r = [], [], [], []
for sid, r in dyn.iterrows():
    b = [k for k in cond if r[inb[k]]]
    u = [k for k in cond if not r[inb[k]]]
    bound_m.append(np.mean([r[f"m_{k}"] for k in b])); unb_m.append(np.mean([r[f"m_{k}"] for k in u]))
    bound_r.append(np.mean([r[f"r_{k}"] for k in b])); unb_r.append(np.mean([r[f"r_{k}"] for k in u]))
dyn["m_bound"], dyn["m_unbound"], dyn["r_bound"], dyn["r_unbound"] = bound_m, unb_m, bound_r, unb_r
for pat, g in list(dyn.groupby("binding_category")) + [("All Dynamic", dyn)]:
    pat_rows.append({"binding_pattern": pat, "n": len(g), "measured_log2_bound_minus_unbound": (g.m_bound - g.m_unbound).mean(), "p_measured": wilcoxon(g.m_bound - g.m_unbound).pvalue,
                     "overprediction_bound": g.r_bound.mean(), "overprediction_unbound": g.r_unbound.mean(), "overprediction_diff_unbound_minus_bound": (g.r_unbound - g.r_bound).mean(),
                     "p_overprediction": wilcoxon(g.r_unbound - g.r_bound).pvalue})
B = pd.DataFrame(pat_rows)
B["fdr_overprediction"] = np.nan
B.loc[B.binding_pattern != "All Dynamic", "fdr_overprediction"] = bh(B[B.binding_pattern != "All Dynamic"].p_overprediction.values)
B.to_csv(OUT / "tables" / "condition_specific_accessibility.csv", index=False)
mm["m_sd"] = mm[["m_YS", "m_MAS", "m_MAE"]].std(axis=1)
Bv = {"stable_sd": mm[mm.dyn == 0].m_sd.median(), "dynamic_sd": mm[mm.dyn == 1].m_sd.median(), "p": mannwhitneyu(mm[mm.dyn == 0].m_sd, mm[mm.dyn == 1].m_sd).pvalue}
pd.DataFrame([Bv]).to_csv(OUT / "tables" / "accessibility_variability_across_conditions.csv", index=False)
# does the prediction track the condition in which a site is bound?
cor_rows = []
for pat, g in list(mm[mm.dyn == 1].groupby("binding_category")) + [("Stable (YS+MAS+MAE)", mm[mm.dyn == 0])]:
    r = {"binding_pattern": pat, "n": len(g)}
    for k in cond:
        r[f"rho_pred_vs_{k}"] = spearmanr(g.native_peak, g[cond[k]])[0]
    cor_rows.append(r)
COR = pd.DataFrame(cor_rows)
COR.to_csv(OUT / "tables" / "prediction_vs_each_condition.csv", index=False)

# ---------------------------------------------------------------- C. footprint strength in the bound conditions
sc = {"YS": "YS_score", "MAS": "MAS_score", "MAE": "MAE_score"}
fp_bound = []
for sid, r in df.iterrows():
    b = [r[sc[k]] for k in sc if r[inb[k]] and pd.notna(r[sc[k]])]
    fp_bound.append(np.mean(b) if b else np.nan)
df["fp_bound"] = fp_bound
C = {"stable_median": df[df.dyn == 0].fp_bound.median(), "dynamic_median": df[df.dyn == 1].fp_bound.median(), "p": mannwhitneyu(df[df.dyn == 0].fp_bound.dropna(), df[df.dyn == 1].fp_bound.dropna()).pvalue}
Cc = pd.DataFrame([C])
Cc.to_csv(OUT / "tables" / "footprint_strength_when_bound.csv", index=False)

# ---------------------------------------------------------------- D. same-peak pairs
pk = df.groupby("peak_coord").dyn.agg(["min", "max", "size"])
mixed = pk[(pk["min"] == 0) & (pk["max"] == 1)].index
rows = []
for c, nm in (("eff_170", "170 bp shuffle effect"), ("eff_motif", "motif shuffle effect"), ("eff_peak", "peak shuffle effect"), ("l_170", "native accessibility, 170 bp (log10)"), ("l_motif", "native accessibility, motif (log10)")):
    g = df[df.peak_coord.isin(mixed)].groupby(["peak_coord", "dyn"])[c].mean().unstack().dropna()
    d = g[1] - g[0]
    rows.append({"quantity": nm, "n_peaks": len(g), "median_dynamic_minus_stable": d.median(), "mean": d.mean(), "ci_lo": np.percentile([np.median(rng.choice(d.values, len(d))) for _ in range(1000)], 2.5),
                 "ci_hi": np.percentile([np.median(rng.choice(d.values, len(d))) for _ in range(1000)], 97.5), "p_wilcoxon": wilcoxon(d).pvalue})
D = pd.DataFrame(rows)
D.to_csv(OUT / "tables" / "same_peak_pairs.csv", index=False)

# ---------------------------------------------------------------- E. position in peak
E = {"stable_median": df[df.dyn == 0].pos_in_peak.median(), "dynamic_median": df[df.dyn == 1].pos_in_peak.median(), "p": mannwhitneyu(df[df.dyn == 0].pos_in_peak.dropna(), df[df.dyn == 1].pos_in_peak.dropna()).pvalue,
     "rho_pos_vs_eff170": spearmanr(df.pos_in_peak, df.eff_170, nan_policy="omit")[0]}
pd.DataFrame([E]).to_csv(OUT / "tables" / "position_in_peak.csv", index=False)

# ---------------------------------------------------------------- F. cross-cell-type concordance
shared = sites.index.intersection(fb_sites.index)
x = pd.DataFrame({"cm": sites.loc[shared, "stable_dynamic"], "fb": fb_sites.loc[shared, "stable_dynamic"]})
ct = pd.crosstab(x.cm, x.fb)
po = (x.cm == x.fb).mean()
pe = sum((x.cm == k).mean() * (x.fb == k).mean() for k in ("Stable", "Dynamic"))
kappa = (po - pe) / (1 - pe)
F = {"shared_sites": len(shared), "stable_both": int(ct.loc["Stable", "Stable"]), "dynamic_both": int(ct.loc["Dynamic", "Dynamic"]), "stable_cm_dynamic_fb": int(ct.loc["Stable", "Dynamic"]), "dynamic_cm_stable_fb": int(ct.loc["Dynamic", "Stable"]),
     "agreement": po, "cohens_kappa": kappa, "p_fb_stable_given_cm_stable": (x[x.cm == "Stable"].fb == "Stable").mean(), "p_fb_stable_given_cm_dynamic": (x[x.cm == "Dynamic"].fb == "Stable").mean()}
pd.DataFrame([F]).to_csv(OUT / "tables" / "cross_cell_type_concordance.csv", index=False)

# ---------------------------------------------------------------- G. matched pairs
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss"]
m_ = df.dropna(subset=MCOLS + ["eff_170", "eff_motif", "eff_peak"]).copy()
Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
dyn_ids = m_.index[m_.dyn == 1].tolist()
st_ids = m_.index[m_.dyn == 0].tolist()
avail = set(st_ids)
pairs = []
order = rng.permutation(dyn_ids)
for d_id in order:
    cand = np.array(sorted(avail))
    dist = np.sqrt(((Z.loc[cand].values - Z.loc[d_id].values) ** 2).sum(1))
    j = int(np.argmin(dist))
    if dist[j] < 0.5:                       # caliper in standard deviations
        pairs.append((d_id, cand[j])); avail.discard(cand[j])
P = pd.DataFrame(pairs, columns=["dyn", "stable"])
bal = []
for c in MCOLS:
    pre = (m_[m_.dyn == 1][c].mean() - m_[m_.dyn == 0][c].mean()) / m_[c].std()
    post = (m_.loc[P.dyn, c].mean() - m_.loc[P.stable, c].mean()) / m_[c].std()
    bal.append({"covariate": c, "std_diff_before": pre, "std_diff_after": post})
pd.DataFrame(bal).to_csv(OUT / "tables" / "matching_balance.csv", index=False)
rows = []
for c, nm in (("eff_motif", "motif shuffle"), ("eff_170", "170 bp shuffle"), ("eff_peak", "peak shuffle")):
    raw = m_[m_.dyn == 1][c].median() - m_[m_.dyn == 0][c].median()
    d = m_.loc[P.dyn, c].values - m_.loc[P.stable, c].values
    bs = [np.median(rng.choice(d, len(d))) for _ in range(1000)]
    rows.append({"effect": nm, "n_pairs": len(P), "n_dynamic_total": len(dyn_ids), "raw_gap": raw, "matched_gap_median": np.median(d), "matched_ci_lo": np.percentile(bs, 2.5), "matched_ci_hi": np.percentile(bs, 97.5),
                 "p_wilcoxon": wilcoxon(d).pvalue})
G = pd.DataFrame(rows)
G.to_csv(OUT / "tables" / "matched_pairs.csv", index=False)

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(14, 4.4), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axes[0]
Ap = A.iloc[::-1].reset_index(drop=True)
ax.barh(range(len(Ap)), Ap.auc_trees, xerr=Ap.auc_trees_sd, color=[GREY if "footprint" in f else BLUE for f in Ap.features], alpha=0.85, ecolor=INK)
ax.axvline(0.5, color=MUTED, lw=0.8, ls=":")
ax.set_yticks(range(len(Ap)), [f.replace(" (", "\n(") for f in Ap.features], fontsize=8)
for i, v in enumerate(Ap.auc_trees):
    ax.text(v + 0.02, i, f"{v:.2f}", va="center", fontsize=8.5)
ax.set_xlim(0.4, 1.0)
ax.set_xlabel("Cross-validated AUC for telling Dynamic from Stable sites (grouped by peak)")
ax.set_title("How separable are the classes?", fontsize=10, loc="left")
ax = axes[1]
PIp = PI.iloc[::-1]
LAB = {"l_motif": "native accessibility, motif", "l_170": "native accessibility, 170 bp", "l_peak": "native accessibility, peak", "l_pw": "peak width", "l_tss": "TSS distance", "motifs_in_peak": "motifs in peak",
       "pos_in_peak": "position in peak", "eff_motif": "motif shuffle effect", "eff_170": "170 bp shuffle effect", "eff_peak": "peak shuffle effect", "l_meas": "measured accessibility", "resid": "prediction error (pred - meas rank)"}
ax.barh(range(len(PIp)), PIp.auc_drop, xerr=PIp["sd"], color=ORANGE, alpha=0.85, ecolor=INK)
ax.set_yticks(range(len(PIp)), [LAB[f] for f in PIp.feature], fontsize=8.5)
ax.set_xlabel("Drop in AUC when the feature is scrambled")
ax.set_title("Which features carry the separation (full model)", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y1_classifier.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
Bp = B[B.binding_pattern != "All Dynamic"].reset_index(drop=True)
for ax, col, ttl, xl in ((axes[0], "measured_log2_bound_minus_unbound", "Measured: more open when bound?", "Measured accessibility, bound minus unbound conditions (log2)"),
                         (axes[1], "overprediction_diff_unbound_minus_bound", "Model: over-predicts more when unbound?", "Over-prediction in unbound minus bound conditions (rank units)")):
    ax.barh(range(len(Bp)), Bp[col], color=ORANGE, alpha=0.85)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(range(len(Bp)), [f"{p} (n={n})" for p, n in zip(Bp.binding_pattern, Bp.n)], fontsize=9)
    ax.set_xlabel(xl, fontsize=8.5)
    ax.set_title(ttl, fontsize=10, loc="left")
for i, r in Bp.iterrows():
    axes[1].text(r.overprediction_diff_unbound_minus_bound + (0.004 if r.overprediction_diff_unbound_minus_bound >= 0 else -0.004), i, "***" if r.fdr_overprediction < 0.001 else "**" if r.fdr_overprediction < 0.01 else "*" if r.fdr_overprediction < 0.05 else "ns",
                 va="center", ha="left" if r.overprediction_diff_unbound_minus_bound >= 0 else "right", fontsize=8)
ax = axes[2]
for k, (nm, col) in enumerate((("Stable (all three conditions)", BLUE), ("Dynamic (bound in some)", ORANGE))):
    ax.hist(mm[mm.dyn == k].m_sd, bins=np.linspace(0, 1.5, 31), density=True, histtype="step", color=col, lw=1.8, label=nm)
ax.set_xlabel("SD of measured accessibility across YS, MAS, MAE (log2)")
ax.set_ylabel("Density")
ax.legend(frameon=False, fontsize=8.5)
ax.set_title("Does measured accessibility change more at Dynamic sites?", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y2_condition_specific.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.4, 3.6))
raw = G[G.effect == "170 bp shuffle"].iloc[0]
rows_ = [("Raw difference of medians (all sites)", raw.raw_gap, np.nan, np.nan),
         (f"Matched pairs ({int(raw.n_pairs)} of {int(raw.n_dynamic_total)} Dynamic sites)", raw.matched_gap_median, raw.matched_ci_lo, raw.matched_ci_hi),
         (f"Same-peak pairs ({int(D.iloc[0].n_peaks)} peaks with both classes)", D.iloc[0].median_dynamic_minus_stable, D.iloc[0].ci_lo, D.iloc[0].ci_hi)]
for i, (nm, v, lo, hi) in enumerate(rows_[::-1]):
    ax.errorbar(v, i, xerr=None if np.isnan(lo) else [[v - lo], [hi - v]], fmt="o", color=BLUE, ms=7, lw=1.5, capsize=0)
ax.set_yticks(range(3), [r[0] for r in rows_[::-1]])
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_xlabel("Dynamic minus Stable, 170 bp shuffle effect (log2FC); positive = Dynamic loses less")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y3_matched_and_same_peak.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- H. is Dynamic mostly a threshold effect on a continuous footprint score?
import statsmodels.formula.api as smf
THR = 0.542                                                  # smallest footprint score called "bound" in the TOBIAS table
df["fp_mean"] = df[["YS_score", "MAS_score", "MAE_score"]].mean(axis=1)
df["fp_sd"] = df[["YS_score", "MAS_score", "MAE_score"]].std(axis=1)
b_scores, u_scores = [], []
for sid, r in df[df.dyn == 1].iterrows():
    for k in sc:
        (b_scores if r[inb[k]] else u_scores).append(r[sc[k]])
H1 = {"dynamic_bound_median": np.median(b_scores), "dynamic_bound_p25": np.percentile(b_scores, 25), "dynamic_bound_p75": np.percentile(b_scores, 75),
      "dynamic_unbound_median": np.median(u_scores), "dynamic_unbound_p25": np.percentile(u_scores, 25), "dynamic_unbound_p75": np.percentile(u_scores, 75),
      "threshold": THR, "stable_median_score": df[df.dyn == 0][["YS_score", "MAS_score", "MAE_score"]].stack().median(),
      "share_dynamic_bound_calls_within_0.15_of_threshold": float(np.mean(np.abs(np.array(b_scores) - THR) < 0.15)),
      "share_dynamic_unbound_calls_within_0.15_of_threshold": float(np.mean(np.abs(np.array(u_scores) - THR) < 0.15)),
      "median_score_sd_across_conditions_stable": df[df.dyn == 0].fp_sd.median(), "median_score_sd_across_conditions_dynamic": df[df.dyn == 1].fp_sd.median()}
pd.DataFrame([H1]).to_csv(OUT / "tables" / "footprint_threshold.csv", index=False)
df["fp_bin"] = pd.qcut(df.fp_mean, 8, duplicates="drop")
HB = df.groupby("fp_bin").agg(n=("dyn", "size"), share_dynamic=("dyn", "mean"), fp_mean=("fp_mean", "median")).reset_index(drop=True)
HB.to_csv(OUT / "tables" / "dynamic_share_by_footprint_strength.csv", index=False)
rows = []
for q, g in df.dropna(subset=["eff_170"]).groupby("fp_bin"):
    a, b = g[g.dyn == 1].eff_170, g[g.dyn == 0].eff_170
    rows.append({"footprint_bin_median": g.fp_mean.median(), "n_stable": len(b), "n_dynamic": len(a), "median_stable": b.median() if len(b) else np.nan, "median_dynamic": a.median() if len(a) else np.nan,
                 "gap_dynamic_minus_stable": a.median() - b.median() if len(a) >= 10 and len(b) >= 10 else np.nan})
HG = pd.DataFrame(rows)
HG.to_csv(OUT / "tables" / "gap_within_footprint_bins.csv", index=False)
dm = df.dropna(subset=["eff_170", "fp_mean", "l_170", "l_peak", "l_pw"]).copy()
dm["grp"] = pd.factorize(dm.peak_coord)[0]
mods = [("raw", "eff_170 ~ dyn"), ("+ mean footprint score", "eff_170 ~ dyn + fp_mean"), ("+ native accessibility, peak width", "eff_170 ~ dyn + l_170 + l_peak + l_pw"),
        ("+ native accessibility, peak width, mean footprint score", "eff_170 ~ dyn + l_170 + l_peak + l_pw + fp_mean")]
rows = []
for nm, f in mods:
    r = smf.ols(f, dm).fit(cov_type="cluster", cov_kwds={"groups": dm.grp})
    ci = r.conf_int().loc["dyn"]
    rows.append({"model": nm, "gap_dynamic_minus_stable": r.params["dyn"], "ci_lo": ci[0], "ci_hi": ci[1], "p": r.pvalues["dyn"], "n": len(dm)})
HM = pd.DataFrame(rows)
HM.to_csv(OUT / "tables" / "gap_models_with_mean_footprint.csv", index=False)
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
ax = axes[0]
bins = np.linspace(0.2, 1.6, 36)
ax.hist(np.array(u_scores), bins=bins, color=GREY, alpha=0.75, label="Dynamic sites, conditions called unbound")
ax.hist(np.array(b_scores), bins=bins, color=ORANGE, alpha=0.75, label="Dynamic sites, conditions called bound")
ax.axvline(THR, color=INK, lw=1.2, ls="--")
ax.text(THR - 0.02, ax.get_ylim()[1] * 0.97, "binding threshold", fontsize=8.5, ha="right")
ax.set_xlabel("TOBIAS footprint score"); ax.set_ylabel("Site-conditions")
ax.legend(frameon=False, fontsize=8, loc="center right")
ax.set_title("Dynamic calls sit at the threshold", fontsize=10, loc="left")
ax = axes[1]
ax.plot(HB.fp_mean, HB.share_dynamic, "o-", color=BLUE, lw=1.8)
ax.set_xlabel("Mean footprint score over the three conditions (bin median)"); ax.set_ylabel("Share of sites that are Dynamic")
ax.set_title("Weak footprints are the Dynamic ones", fontsize=10, loc="left")
ax = axes[2]
ok = HG.dropna(subset=["gap_dynamic_minus_stable"])
ax.plot(HG.footprint_bin_median, HG.median_stable, "o-", color=BLUE, label="Stable")
ax.plot(HG.footprint_bin_median, HG.median_dynamic, "o-", color=ORANGE, label="Dynamic")
ax.set_xlabel("Mean footprint score (bin median)"); ax.set_ylabel("Median 170 bp shuffle effect (log2FC)")
ax.legend(frameon=False)
ax.set_title("At equal footprint strength, classes look alike", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y4_footprint_threshold.png", bbox_inches="tight")
plt.close(fig)
print("\n== H threshold\n", pd.DataFrame([H1]).round(3).T.to_string(), "\n", HB.round(3).to_string(index=False), "\n", HG.round(3).to_string(index=False), "\n", HM.round(3).to_string(index=False))

pd.set_option("display.width", 230, "display.max_columns", 30)
for nm, t in (("A classifier", A), ("A permutation importance", PI), ("B condition-specific", B), ("B variability", pd.DataFrame([Bv])), ("B prediction vs each condition", COR), ("C footprint when bound", Cc),
              ("D same peak", D), ("E position", pd.DataFrame([E])), ("F cross cell type", pd.DataFrame([F])), ("G balance", pd.DataFrame(bal)), ("G matched", G)):
    print(f"\n== {nm}\n", t.round(4).to_string(index=False))
print("\ncontingency CM x FB\n", ct)
