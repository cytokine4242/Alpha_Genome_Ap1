#!/usr/bin/env python3
"""What in the AlphaGenome results distinguishes Dynamic from Stable sites? Cardiomyocyte sites, heart ATAC track. No TOBIAS footprint scores are used.
 A  can the model's outputs tell the classes apart (cross-validated AUC, folds grouped by peak)?
 B  condition-specific accessibility: measured accessibility and the model's over-prediction where a Dynamic site is bound vs unbound
 D  Stable and Dynamic sites in the SAME peak (paired)
 E  position of the site within its peak
 F  does a site's class carry over between cell types (same sequence, cardiomyocyte vs fibroblast)?
 G  matched pairs: each Dynamic site against the Stable site closest in predicted accessibility, peak width and TSS distance
 I  scan of every AlphaGenome output (native predictions, shuffle/mutation/control effects, RNA) for differences between Dynamic sites and their matched Stable sites
Outputs Analysis/dynamic/{tables,figures}.   python3 33_dynamic_vs_stable_deep_dive.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from sklearn.ensemble import HistGradientBoostingClassifier
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
R = ROOT / "Katana" / "results"
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


def load(prefix):
    d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{prefix}_chunk*.csv")))], ignore_index=True)
    return d[d.trackset == "CM"]


def eff_table(d, col):
    """per-site mean log2((edit + pc) / (native + pc)) for every condition; pc = 5% of the median native value"""
    nat = d[d.condition == "native"].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[d.condition != "native"].copy()
    e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "condition"]).l.mean().unstack(), nat


sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id")
fb_sites = pd.read_csv(ROOT / "Data/derived/fb_fos_sites.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id")
bycond = pd.read_csv(ROOT / "Katana/data/measured_by_condition.csv")
bycond = bycond[bycond.site_table == "cm"].set_index("site_id")
drop_scores = [c for c in sl.columns if c in ("YS_score", "MAS_score", "MAE_score", "fp_score")]       # TOBIAS scores are not used anywhere below
df = sl.drop(columns=drop_scores).join(sites[["motif_mid", "peak_start", "peak_end", "in_YS", "in_MAS", "in_MAE"]]).join(meas[["cm_cpm_per_kb"]]).join(bycond[["cm_ys", "cm_mas", "cm_mae"]])
df["dyn"] = (df.stable_dynamic == "Dynamic").astype(int)
df["l_peak"], df["l_170"], df["l_motif"] = [np.log10(df[c].clip(lower=1e-3)) for c in ("native_peak", "native_170", "native_motif")]
df["l_pw"], df["l_tss"] = np.log10(df.peak_width), np.log10(df.abs_tss + 1)
df["l_meas"] = np.log10(df.cm_cpm_per_kb.clip(lower=1e-3))
df["pos_in_peak"] = (df.motif_mid - (df.peak_start + df.peak_end) / 2).abs() / (df.peak_width / 2)
df["resid"] = df.native_peak.rank(pct=True) - df.cm_cpm_per_kb.rank(pct=True)

# ---------------------------------------------------------------- per-site table of every AlphaGenome output
sh, mu, ct = load("shuffle_cm"), load("mut_cm"), load("controls_cm")
F, fam = {}, {}


def add(name, s, family):
    F[name] = s
    fam[name] = family


RDN = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
RDL = {"motif": "motif", "170": "170 bp", "peak": "peak"}
tabs = {rd: eff_table(sh, c) for rd, c in RDN.items()}
for w in RDN:
    for rd in RDN:
        add(f"shuffle {RDL[w]}, read over {RDL[rd]}", tabs[rd][0][f"shuffle@{w}"], "shuffle effects")
mtabs = {rd: eff_table(mu, c)[0] for rd, c in RDN.items()}
for cond, rd, lab in (("random@motif", "motif", "replace motif"), ("dinuc@motif", "motif", "dinucleotide-shuffle motif"), ("core@motif", "motif", "replace core"), ("flank@motif", "motif", "replace motif flanks"),
                      ("random@170", "170", "replace 170 bp"), ("dinuc@170", "170", "dinucleotide-shuffle 170 bp"), ("flank@170", "170", "replace surroundings (motif kept)"),
                      ("random@peak", "peak", "replace peak"), ("dinuc@peak", "peak", "dinucleotide-shuffle peak")):
    add(f"{lab}, read over {RDL[rd]}", mtabs[rd][cond], "mutation effects")
ctabs = {rd: eff_table(ct, c)[0] for rd, c in RDN.items()}
for rd in ("motif", "peak"):
    add(f"170 bp shuffle ±300 bp away, read over {RDL[rd]}", (ctabs[rd]["offset@-300"] + ctabs[rd]["offset@+300"]) / 2, "control and titration effects")
    add(f"170 bp shuffle ±600 bp away, read over {RDL[rd]}", (ctabs[rd]["offset@-600"] + ctabs[rd]["offset@+600"]) / 2, "control and titration effects")
add("170 bp shuffle −300 minus +300 bp (left-right asymmetry), read over motif", ctabs["motif"]["offset@-300"] - ctabs["motif"]["offset@+300"], "control and titration effects")
for w in (7, 11, 21, 41, 85, 340, 680):
    add(f"centred shuffle of {w} bp, read over motif", ctabs["motif"][f"titrate@{w}"], "control and titration effects")
for col, lab in (("rna_prox2kb", "RNA within 2 kb of nearest TSS"), ("rna_mean50kb", "RNA, mean ±50 kb")):
    nat = sh[sh.condition == "native"].set_index("site_id")[col]
    add(f"native predicted {lab} (log10)", np.log10(nat.clip(lower=1e-6)), "RNA predictions")
    t_, _ = eff_table(sh, col)
    for w in RDN:
        add(f"{RDL[w]} shuffle: effect on {lab}", t_[f"shuffle@{w}"], "RNA predictions")
add("native predicted accessibility, motif (log10)", df.l_motif, "native accessibility")
add("native predicted accessibility, 170 bp (log10)", df.l_170, "native accessibility")
add("native predicted accessibility, peak (log10)", df.l_peak, "native accessibility")
FD = pd.DataFrame(F)
FD["motif share: 170 bp shuffle minus motif shuffle, read over motif"] = FD["shuffle 170 bp, read over motif"] - FD["shuffle motif, read over motif"]
fam["motif share: 170 bp shuffle minus motif shuffle, read over motif"] = "derived contrasts"
FD["what the motif adds: replace 170 bp minus replace surroundings, read over 170 bp"] = FD["replace 170 bp, read over 170 bp"] - FD["replace surroundings (motif kept), read over 170 bp"]
fam["what the motif adds: replace 170 bp minus replace surroundings, read over 170 bp"] = "derived contrasts"
FD["what the core adds: replace motif minus replace core, read over motif"] = FD["replace motif, read over motif"] - FD["replace core, read over motif"]
fam["what the core adds: replace motif minus replace core, read over motif"] = "derived contrasts"
FD["locality: centred 170 bp shuffle minus ±300 bp shuffle, read over motif"] = FD["shuffle 170 bp, read over motif"] - FD["170 bp shuffle ±300 bp away, read over motif"]
fam["locality: centred 170 bp shuffle minus ±300 bp shuffle, read over motif"] = "derived contrasts"
FD["reach: 170 bp shuffle read over peak minus read over 170 bp"] = FD["shuffle 170 bp, read over peak"] - FD["shuffle 170 bp, read over 170 bp"]
fam["reach: 170 bp shuffle read over peak minus read over 170 bp"] = "derived contrasts"
FD.to_csv(OUT / "tables" / "per_site_alphagenome_features.csv")
df = df.join(FD, how="left", rsuffix="_f") if False else df

# ---------------------------------------------------------------- G. matched pairs (also the sample for the feature scan)
# measured accessibility belongs in the match: matching on the prediction alone leaves Dynamic sites 0.8 SD lower in the
# measured data, which is itself most of what makes a wide edit cost less. With it, all five covariates balance within 0.06 SD.
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]
m_ = df.dropna(subset=MCOLS + ["eff_170", "eff_motif", "eff_peak"]).copy()
Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
dyn_ids = m_.index[m_.dyn == 1].tolist()
avail = set(m_.index[m_.dyn == 0])
pairs = []
for d_id in rng.permutation(dyn_ids):
    cand = np.array(sorted(avail))
    dist = np.sqrt(((Z.loc[cand].values - Z.loc[d_id].values) ** 2).sum(1))
    j = int(np.argmin(dist))
    if dist[j] < 0.5:
        pairs.append((d_id, cand[j])); avail.discard(cand[j])
P = pd.DataFrame(pairs, columns=["dyn", "stable"])
P.to_csv(OUT / "tables" / "matched_pair_ids.csv", index=False)
bal = []
for c in MCOLS:
    bal.append({"covariate": c, "std_diff_before": (m_[m_.dyn == 1][c].mean() - m_[m_.dyn == 0][c].mean()) / m_[c].std(), "std_diff_after": (m_.loc[P.dyn, c].mean() - m_.loc[P.stable, c].mean()) / m_[c].std()})
pd.DataFrame(bal).to_csv(OUT / "tables" / "matching_balance.csv", index=False)
rows = []
for c, nm in (("eff_motif", "motif shuffle"), ("eff_170", "170 bp shuffle"), ("eff_peak", "peak shuffle")):
    raw = m_[m_.dyn == 1][c].median() - m_[m_.dyn == 0][c].median()
    d = m_.loc[P.dyn, c].values - m_.loc[P.stable, c].values
    bs = [np.median(rng.choice(d, len(d))) for _ in range(1000)]
    rows.append({"effect": nm, "n_pairs": len(P), "n_dynamic_total": len(dyn_ids), "raw_gap": raw, "matched_gap_median": np.median(d), "matched_ci_lo": np.percentile(bs, 2.5), "matched_ci_hi": np.percentile(bs, 97.5), "p_wilcoxon": wilcoxon(d).pvalue})
G = pd.DataFrame(rows)
G.to_csv(OUT / "tables" / "matched_pairs.csv", index=False)

# ---------------------------------------------------------------- I. scan of every AlphaGenome output in matched pairs
rows = []
dyn_all, st_all = FD.loc[FD.index.intersection(df.index[df.dyn == 1])], FD.loc[FD.index.intersection(df.index[df.dyn == 0])]
for c in FD.columns:
    a, b = FD.loc[P.dyn, c].values, FD.loc[P.stable, c].values
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 30:
        continue
    d = a[ok] - b[ok]
    sd = np.nanstd(FD.loc[df.index[df.dyn == 0], c])
    bs = [np.median(rng.choice(d, len(d))) for _ in range(500)]
    rows.append({"feature": c, "family": fam[c], "n_pairs": int(ok.sum()), "median_stable": np.median(b[ok]), "median_dynamic": np.median(a[ok]), "matched_median_diff": np.median(d), "ci_lo": np.percentile(bs, 2.5), "ci_hi": np.percentile(bs, 97.5),
                 "matched_std_diff": np.mean(d) / sd, "p_wilcoxon": wilcoxon(d).pvalue,
                 "raw_std_diff": (np.nanmean(dyn_all[c]) - np.nanmean(st_all[c])) / sd})
SC = pd.DataFrame(rows)
SC["fdr"] = bh(SC.p_wilcoxon.values)
SC = SC.sort_values("matched_std_diff", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
SC.to_csv(OUT / "tables" / "feature_scan_matched_pairs.csv", index=False)

# classifiers on the matched sample: which sets of model outputs can tell Dynamic from matched Stable?
ids = list(P.dyn) + list(P.stable)
y = np.array([1] * len(P) + [0] * len(P))
Xall = FD.loc[ids].copy()
Xall["l_pw"], Xall["l_tss"] = df.loc[ids, "l_pw"].values, df.loc[ids, "l_tss"].values
grp = pd.factorize(df.loc[ids, "peak_coord"])[0]
SETS = {
    "matching variables (native accessibility, peak width, TSS distance)": ["native predicted accessibility, 170 bp (log10)", "native predicted accessibility, peak (log10)", "l_pw", "l_tss"],
    "shuffle effects (9)": [c for c in FD.columns if fam[c] == "shuffle effects"],
    "mutation effects (9)": [c for c in FD.columns if fam[c] == "mutation effects"],
    "control and titration effects": [c for c in FD.columns if fam[c] == "control and titration effects"],
    "RNA predictions": [c for c in FD.columns if fam[c] == "RNA predictions"],
    "all effects and contrasts (no native values)": [c for c in FD.columns if fam[c] in ("shuffle effects", "mutation effects", "control and titration effects", "derived contrasts")],
    "everything AlphaGenome outputs": list(FD.columns) + ["l_pw", "l_tss"],
}
rows = []
for nm, cols in SETS.items():
    X = Xall[cols].copy()
    X = X.fillna(X.median())
    X = X.values
    aucs = {"logistic": [], "trees": []}
    for tr, te in GroupKFold(5).split(X, y, grp):
        lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, C=0.3)).fit(X[tr], y[tr])
        gb = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=150, random_state=0).fit(X[tr], y[tr])
        aucs["logistic"].append(roc_auc_score(y[te], lr.predict_proba(X[te])[:, 1]))
        aucs["trees"].append(roc_auc_score(y[te], gb.predict_proba(X[te])[:, 1]))
    rows.append({"feature_set": nm, "n_features": len(cols), "auc_logistic": np.mean(aucs["logistic"]), "auc_logistic_sd": np.std(aucs["logistic"]), "auc_trees": np.mean(aucs["trees"]), "auc_trees_sd": np.std(aucs["trees"]), "n_sites": len(y)})
CM_AUC = pd.DataFrame(rows)
CM_AUC.to_csv(OUT / "tables" / "classifier_auc_matched_sample.csv", index=False)

# ---------------------------------------------------------------- A. classifier on all sites (model outputs and prediction error only)
FEATS = {
    "native predicted accessibility only (motif, 170 bp, peak)": ["l_motif", "l_170", "l_peak"],
    "peak geometry and position (width, TSS distance, motifs in peak, position in peak)": ["l_pw", "l_tss", "motifs_in_peak", "pos_in_peak"],
    "model effects only (motif, 170 bp, peak shuffle)": ["eff_motif", "eff_170", "eff_peak"],
    "context: native accessibility + geometry": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak"],
    "context + model effects": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak", "eff_motif", "eff_170", "eff_peak"],
    "measured accessibility only": ["l_meas"],
    "context + effects + measured accessibility + prediction error": ["l_motif", "l_170", "l_peak", "l_pw", "l_tss", "motifs_in_peak", "pos_in_peak", "eff_motif", "eff_170", "eff_peak", "l_meas", "resid"],
}
rows = []
d_ = df.dropna(subset=sum(FEATS.values(), [])).copy()
groups = pd.factorize(d_.peak_coord)[0]
for name, cols in FEATS.items():
    X, yy = d_[cols].values, d_.dyn.values
    aucs = {"logistic": [], "boosted trees": []}
    for tr, te in GroupKFold(5).split(X, yy, groups):
        lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X[tr], yy[tr])
        gb = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.06, max_iter=200, random_state=0).fit(X[tr], yy[tr])
        aucs["logistic"].append(roc_auc_score(yy[te], lr.predict_proba(X[te])[:, 1]))
        aucs["boosted trees"].append(roc_auc_score(yy[te], gb.predict_proba(X[te])[:, 1]))
    rows.append({"features": name, "n_features": len(cols), "auc_logistic": np.mean(aucs["logistic"]), "auc_logistic_sd": np.std(aucs["logistic"]), "auc_trees": np.mean(aucs["boosted trees"]), "auc_trees_sd": np.std(aucs["boosted trees"]), "n_sites": len(d_)})
A = pd.DataFrame(rows)
A.to_csv(OUT / "tables" / "classifier_auc.csv", index=False)

# ---------------------------------------------------------------- B. condition-specific accessibility and prediction error
cond = {"YS": "cm_ys", "MAS": "cm_mas", "MAE": "cm_mae"}
inb = {"YS": "in_YS", "MAS": "in_MAS", "MAE": "in_MAE"}
mm = df.dropna(subset=list(cond.values()) + ["native_peak"]).copy()
for k, c in cond.items():
    mm[f"m_{k}"] = np.log2(mm[c].clip(lower=1e-3))
    mm[f"p_{k}"] = mm[c].rank(pct=True)
    mm[f"r_{k}"] = mm.native_peak.rank(pct=True) - mm[f"p_{k}"]
dyn = mm[mm.dyn == 1].copy()
bm, um, br, ur = [], [], [], []
for sid, r in dyn.iterrows():
    b = [k for k in cond if r[inb[k]]]
    u = [k for k in cond if not r[inb[k]]]
    bm.append(np.mean([r[f"m_{k}"] for k in b])); um.append(np.mean([r[f"m_{k}"] for k in u]))
    br.append(np.mean([r[f"r_{k}"] for k in b])); ur.append(np.mean([r[f"r_{k}"] for k in u]))
dyn["m_bound"], dyn["m_unbound"], dyn["r_bound"], dyn["r_unbound"] = bm, um, br, ur
pat_rows = []
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
cor_rows = []
for pat, g in list(mm[mm.dyn == 1].groupby("binding_category")) + [("Stable (YS+MAS+MAE)", mm[mm.dyn == 0])]:
    r = {"binding_pattern": pat, "n": len(g)}
    for k in cond:
        r[f"rho_pred_vs_{k}"] = spearmanr(g.native_peak, g[cond[k]])[0]
    cor_rows.append(r)
pd.DataFrame(cor_rows).to_csv(OUT / "tables" / "prediction_vs_each_condition.csv", index=False)

# ---------------------------------------------------------------- D. same-peak pairs, E. position, F. cross cell type
pk = df.groupby("peak_coord").dyn.agg(["min", "max", "size"])
mixed = pk[(pk["min"] == 0) & (pk["max"] == 1)].index
rows = []
for c, nm in (("eff_170", "170 bp shuffle effect"), ("eff_motif", "motif shuffle effect"), ("eff_peak", "peak shuffle effect"), ("l_170", "native accessibility, 170 bp (log10)"), ("l_motif", "native accessibility, motif (log10)")):
    g = df[df.peak_coord.isin(mixed)].groupby(["peak_coord", "dyn"])[c].mean().unstack().dropna()
    d = g[1] - g[0]
    bs = [np.median(rng.choice(d.values, len(d))) for _ in range(1000)]
    rows.append({"quantity": nm, "n_peaks": len(g), "median_dynamic_minus_stable": d.median(), "mean": d.mean(), "ci_lo": np.percentile(bs, 2.5), "ci_hi": np.percentile(bs, 97.5), "p_wilcoxon": wilcoxon(d).pvalue})
D = pd.DataFrame(rows)
D.to_csv(OUT / "tables" / "same_peak_pairs.csv", index=False)
E = {"stable_median": df[df.dyn == 0].pos_in_peak.median(), "dynamic_median": df[df.dyn == 1].pos_in_peak.median(), "p": mannwhitneyu(df[df.dyn == 0].pos_in_peak.dropna(), df[df.dyn == 1].pos_in_peak.dropna()).pvalue,
     "rho_pos_vs_eff170": spearmanr(df.pos_in_peak, df.eff_170, nan_policy="omit")[0]}
pd.DataFrame([E]).to_csv(OUT / "tables" / "position_in_peak.csv", index=False)
shared = sites.index.intersection(fb_sites.index)
x = pd.DataFrame({"cm": sites.loc[shared, "stable_dynamic"], "fb": fb_sites.loc[shared, "stable_dynamic"]})
ct_ = pd.crosstab(x.cm, x.fb)
po = (x.cm == x.fb).mean()
pe = sum((x.cm == k).mean() * (x.fb == k).mean() for k in ("Stable", "Dynamic"))
Fx = {"shared_sites": len(shared), "stable_both": int(ct_.loc["Stable", "Stable"]), "dynamic_both": int(ct_.loc["Dynamic", "Dynamic"]), "stable_cm_dynamic_fb": int(ct_.loc["Stable", "Dynamic"]), "dynamic_cm_stable_fb": int(ct_.loc["Dynamic", "Stable"]),
      "agreement": po, "cohens_kappa": (po - pe) / (1 - pe), "p_fb_stable_given_cm_stable": (x[x.cm == "Stable"].fb == "Stable").mean(), "p_fb_stable_given_cm_dynamic": (x[x.cm == "Dynamic"].fb == "Stable").mean()}
pd.DataFrame([Fx]).to_csv(OUT / "tables" / "cross_cell_type_concordance.csv", index=False)

# ---------------------------------------------------------------- figures
FAMCOL = {"shuffle effects": BLUE, "mutation effects": ORANGE, "control and titration effects": AQUA, "RNA predictions": "#8a5cc2", "derived contrasts": "#d1a000", "native accessibility": GREY}
top = SC.head(22).iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11.5, 7.2))
for i, r in top.iterrows():
    sig = r.fdr < 0.05
    ax.scatter(r.matched_std_diff, i, s=55, color=FAMCOL[r.family], zorder=3, facecolor=FAMCOL[r.family] if sig else "none", linewidths=1.6)
    ax.scatter(r.raw_std_diff, i, s=28, marker="D", color=GREY, zorder=2, facecolor="none", linewidths=1)
ax.set_yticks(range(len(top)), top.feature, fontsize=8)
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_xlabel("Dynamic minus Stable in standard deviations of the Stable sites (positive = larger/less negative in Dynamic)")
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker="o", ls="", color=c, label=f) for f, c in FAMCOL.items()] + [Line2D([], [], marker="o", ls="", color=INK, mfc="none", label="matched pairs, FDR ≥ 0.05"), Line2D([], [], marker="D", ls="", color=GREY, mfc="none", label="all sites (unmatched)")]
ax.legend(handles=leg, frameon=False, fontsize=8, loc="lower right")
ax.set_title("AlphaGenome outputs that differ most between Dynamic sites and context-matched Stable sites (22 largest of %d)" % len(SC), fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y4_feature_scan.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axes[0]
Ap = A.iloc[::-1].reset_index(drop=True)
ax.barh(range(len(Ap)), Ap.auc_trees, xerr=Ap.auc_trees_sd, color=BLUE, alpha=0.85, ecolor=INK)
ax.axvline(0.5, color=MUTED, lw=0.8, ls=":")
ax.set_yticks(range(len(Ap)), [f.replace(" (", "\n(") for f in Ap.features], fontsize=8)
for i, v in enumerate(Ap.auc_trees):
    ax.text(v + 0.02, i, f"{v:.2f}", va="center", fontsize=8.5)
ax.set_xlim(0.4, 1.0)
ax.set_xlabel("Cross-validated AUC, Dynamic vs Stable (all sites)")
ax.set_title("All sites", fontsize=10, loc="left")
ax = axes[1]
Cp = CM_AUC.iloc[::-1].reset_index(drop=True)
ax.barh(range(len(Cp)), Cp.auc_logistic, xerr=Cp.auc_logistic_sd, color=ORANGE, alpha=0.85, ecolor=INK)
ax.axvline(0.5, color=MUTED, lw=0.8, ls=":")
ax.set_yticks(range(len(Cp)), [f.replace(" (", "\n(") for f in Cp.feature_set], fontsize=8)
for i, v in enumerate(Cp.auc_logistic):
    ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=8.5)
ax.set_xlim(0.4, 0.8)
ax.set_xlabel("Cross-validated AUC, Dynamic vs context-matched Stable")
ax.set_title("Matched sample", fontsize=10, loc="left")
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
    axes[1].text(r.overprediction_diff_unbound_minus_bound + 0.004, i, "***" if r.fdr_overprediction < 0.001 else "**" if r.fdr_overprediction < 0.01 else "*" if r.fdr_overprediction < 0.05 else "ns", va="center", fontsize=8)
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
rows_ = [("Raw difference of medians (all sites)", raw.raw_gap, np.nan, np.nan), (f"Matched pairs ({int(raw.n_pairs)} of {int(raw.n_dynamic_total)} Dynamic sites)", raw.matched_gap_median, raw.matched_ci_lo, raw.matched_ci_hi),
         (f"Same-peak pairs ({int(D.iloc[0].n_peaks)} peaks with both classes)", D.iloc[0].median_dynamic_minus_stable, D.iloc[0].ci_lo, D.iloc[0].ci_hi)]
for i, (nm, v, lo, hi) in enumerate(rows_[::-1]):
    ax.errorbar(v, i, xerr=None if np.isnan(lo) else [[v - lo], [hi - v]], fmt="o", color=BLUE, ms=7, lw=1.5, capsize=0)
ax.set_yticks(range(3), [r[0] for r in rows_[::-1]])
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_xlabel("Dynamic minus Stable, 170 bp shuffle effect (log2FC); positive = Dynamic loses less")
fig.tight_layout()
fig.savefig(OUT / "figures" / "y3_matched_and_same_peak.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- J. sensitivity to each type of mutation: Dynamic vs Stable, and each binding pattern vs Stable
import statsmodels.formula.api as smf
import patsy
bs = patsy.bs        # (an earlier loop reused the name bs for a bootstrap list; the formulas below need the spline basis)
PROFILE = [
    ("Within the 11 bp motif (read over the motif)", [("shuffle motif, read over motif", "shuffle the motif"), ("replace motif, read over motif", "replace the motif"), ("dinucleotide-shuffle motif, read over motif", "dinucleotide-shuffle the motif"),
                                                       ("replace core, read over motif", "replace the 7 bp core"), ("replace motif flanks, read over motif", "replace the 4 flank bases")]),
    ("Within the 170 bp (read over the 170 bp)", [("shuffle 170 bp, read over 170 bp", "shuffle the 170 bp"), ("replace 170 bp, read over 170 bp", "replace the 170 bp"), ("dinucleotide-shuffle 170 bp, read over 170 bp", "dinucleotide-shuffle the 170 bp"),
                                                   ("replace surroundings (motif kept), read over 170 bp", "replace the 159 bp around the motif")]),
    ("Whole peak (read over the peak)", [("shuffle peak, read over peak", "shuffle the peak"), ("replace peak, read over peak", "replace the peak"), ("dinucleotide-shuffle peak, read over peak", "dinucleotide-shuffle the peak")]),
    ("Edits far from the motif (read over the motif)", [("170 bp shuffle ±300 bp away, read over motif", "170 bp shuffle ±300 bp away"), ("170 bp shuffle ±600 bp away, read over motif", "170 bp shuffle ±600 bp away")]),
    ("Centred shuffles of increasing width (read over the motif)", [(f"centred shuffle of {w} bp, read over motif", f"centred shuffle, {w} bp") for w in (7, 11, 21, 41, 85, 340, 680)]),
]
FEATS_J = [(g, f, lab) for g, items in PROFILE for f, lab in items]
NATCOV = lambda f: "l_motif" if f.endswith("read over motif") else "l_peak" if f.endswith("read over peak") else "l_170"
dd = df[["dyn", "binding_category", "l_motif", "l_170", "l_peak", "l_pw", "l_tss"]].join(FD[[f for _, f, _ in FEATS_J]])
rows = []
for g, f, lab in FEATS_J:
    a_all, b_all = dd[dd.dyn == 1][f].dropna(), dd[dd.dyn == 0][f].dropna()
    a_m, b_m = FD.loc[P.dyn, f].values, FD.loc[P.stable, f].values
    ok = ~(np.isnan(a_m) | np.isnan(b_m))
    ci = lambda v: tuple(np.percentile([np.median(rng.choice(v, len(v))) for _ in range(300)], [2.5, 97.5]))
    rows.append({"group": g, "mutation": lab, "feature": f, "stable_median": b_all.median(), "stable_lo": ci(b_all.values)[0], "stable_hi": ci(b_all.values)[1], "dynamic_median": a_all.median(), "dynamic_lo": ci(a_all.values)[0], "dynamic_hi": ci(a_all.values)[1],
                 "p_raw": mannwhitneyu(a_all, b_all).pvalue, "matched_stable_median": np.median(b_m[ok]), "matched_stable_lo": ci(b_m[ok])[0], "matched_stable_hi": ci(b_m[ok])[1], "matched_dynamic_median": np.median(a_m[ok]),
                 "matched_dynamic_lo": ci(a_m[ok])[0], "matched_dynamic_hi": ci(a_m[ok])[1], "p_matched": wilcoxon(a_m[ok] - b_m[ok]).pvalue})
J = pd.DataFrame(rows)
J["fdr_raw"], J["fdr_matched"] = bh(J.p_raw.values), bh(J.p_matched.values)
J.to_csv(OUT / "tables" / "mutation_type_profile_dynamic_vs_stable.csv", index=False)
stars_ = lambda q: "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"
fig, axes = plt.subplots(1, 2, figsize=(13.5, 7.4), sharey=True)
ypos, y, sep = [], 0, []
for g, items in PROFILE:
    for f, lab in items:
        ypos.append(y); y += 1
    sep.append(y - 0.5); y += 0.7
yy_ = {}
k = 0
for g, items in PROFILE:
    for f, lab in items:
        yy_[f] = ypos[k]; k += 1
for ax, pre, ttl, qcol in ((axes[0], "", "All sites", "fdr_raw"), (axes[1], "matched_", "Context-matched pairs", "fdr_matched")):
    for _, r in J.iterrows():
        y_ = yy_[r.feature]
        ax.errorbar(r[f"{pre}stable_median"], y_ - 0.12, xerr=[[r[f"{pre}stable_median"] - r[f"{pre}stable_lo"]], [r[f"{pre}stable_hi"] - r[f"{pre}stable_median"]]], fmt="o", color=BLUE, ms=5, lw=1.2, capsize=0, label="Stable" if r.name == 0 else None)
        ax.errorbar(r[f"{pre}dynamic_median"], y_ + 0.12, xerr=[[r[f"{pre}dynamic_median"] - r[f"{pre}dynamic_lo"]], [r[f"{pre}dynamic_hi"] - r[f"{pre}dynamic_median"]]], fmt="o", color=ORANGE, ms=5, lw=1.2, capsize=0, label="Dynamic" if r.name == 0 else None)
        ax.text(1.01, y_, stars_(r[qcol]), transform=ax.get_yaxis_transform(), va="center", fontsize=8.5)
    for s_ in sep[:-1]:
        ax.axhline(s_ + 0.35, color=GRID, lw=1)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_title(ttl, fontsize=10, loc="left")
    ax.set_xlabel("Median effect (log2FC); more negative = more sensitive")
axes[0].set_yticks([yy_[f] for _, f, _ in FEATS_J], [lab for _, _, lab in FEATS_J], fontsize=8.5)
axes[0].invert_yaxis()
axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.93))
ytop = 0
for g, items in PROFILE:
    axes[0].text(-0.02, yy_[items[0][0]] - 0.65, g, transform=axes[0].get_yaxis_transform(), ha="right", va="bottom", fontsize=8, fontweight="bold", color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "y5_mutation_type_profile.png", bbox_inches="tight")
plt.close(fig)

# each binding pattern against Stable, adjusted for native accessibility at the readout and peak width
pats = sorted(dd[dd.dyn == 1].binding_category.unique())
rows = []
for g, f, lab in FEATS_J:
    cov = NATCOV(f)
    base = dd[dd.dyn == 0]
    for pat in pats + ["All Dynamic"]:
        sub = dd[dd.dyn == 1] if pat == "All Dynamic" else dd[dd.binding_category == pat]
        x = pd.concat([base.assign(g=0), sub.assign(g=1)]).dropna(subset=[f, cov, "l_peak", "l_pw", "l_tss"])
        terms = [(cov, 5)] + [(c_, 4) for c_ in ("l_peak", "l_pw", "l_tss") if c_ != cov]
        r = smf.ols(f"Q('{f}') ~ g + " + " + ".join(f"bs({c_}, df={d_})" for c_, d_ in terms), x).fit(cov_type="HC3")
        rows.append({"mutation": lab, "feature": f, "pattern": pat, "n": int((x.g == 1).sum()), "raw_median_diff": x[x.g == 1][f].median() - x[x.g == 0][f].median(), "adj_diff": r.params["g"], "adj_lo": r.conf_int().loc["g", 0], "adj_hi": r.conf_int().loc["g", 1], "p": r.pvalues["g"]})
PJ = pd.DataFrame(rows)
sel = PJ.pattern != "All Dynamic"
PJ["fdr"] = np.nan
PJ.loc[sel, "fdr"] = bh(PJ[sel].p.values)
PJ.to_csv(OUT / "tables" / "mutation_type_by_binding_pattern.csv", index=False)
order_m = [lab for _, _, lab in FEATS_J]
fig, axes = plt.subplots(1, 2, figsize=(13.5, 7.2), sharey=True)
for ax, col, ttl in ((axes[0], "raw_median_diff", "Raw difference of medians from Stable"), (axes[1], "adj_diff", "Adjusted for native accessibility (readout and peak), peak width and TSS distance")):
    M = PJ.pivot(index="mutation", columns="pattern", values=col).reindex(order_m)[pats + ["All Dynamic"]]
    vmax = np.nanmax(np.abs(M.values))
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(M.shape[1]), M.columns, rotation=40, ha="right", fontsize=8.5)
    ax.set_title(ttl, fontsize=10, loc="left")
    if col == "adj_diff":
        Q = PJ.pivot(index="mutation", columns="pattern", values="fdr").reindex(order_m)[pats + ["All Dynamic"]]
        Pp = PJ.pivot(index="mutation", columns="pattern", values="p").reindex(order_m)[pats + ["All Dynamic"]]
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            txt = f"{M.values[i, j]:+.2f}"
            if col == "adj_diff":
                q = Q.values[i, j] if M.columns[j] != "All Dynamic" else Pp.values[i, j]
                txt += "\n" + ("***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.8, color="white" if abs(M.values[i, j]) > 0.6 * vmax else INK)
axes[0].set_yticks(range(len(order_m)), order_m, fontsize=8.5)
fig.colorbar(im, ax=axes, shrink=0.6, label="Dynamic pattern minus Stable (log2FC); positive = loses less")
fig.savefig(OUT / "figures" / "y6_mutation_type_by_pattern.png", bbox_inches="tight")
plt.close(fig)
print("\n== J profile\n", J[["mutation", "stable_median", "dynamic_median", "fdr_raw", "matched_stable_median", "matched_dynamic_median", "fdr_matched"]].round(3).to_string(index=False))
print("\n== J patterns significant after BH:", int((PJ.fdr < 0.05).sum()), "of", int(sel.sum()))
print(PJ[PJ.fdr < 0.05][["mutation", "pattern", "n", "raw_median_diff", "adj_diff", "adj_lo", "adj_hi", "fdr"]].round(3).to_string(index=False))

pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 80)
print("== matching balance\n", pd.DataFrame(bal).round(3).to_string(index=False), "\n== matched\n", G.round(3).to_string(index=False))
print("\n== feature scan (top 25)\n", SC.head(25)[["feature", "family", "n_pairs", "median_stable", "median_dynamic", "matched_median_diff", "matched_std_diff", "fdr", "raw_std_diff"]].round(3).to_string(index=False))
print("\nfeatures with FDR<0.05:", int((SC.fdr < 0.05).sum()), "of", len(SC))
print("\n== classifier, matched sample\n", CM_AUC.round(3).to_string(index=False), "\n== classifier all sites\n", A.round(3).to_string(index=False))
print("\n== B\n", B.round(4).to_string(index=False), "\n== same peak\n", D.round(3).to_string(index=False), "\n== cross\n", pd.DataFrame([Fx]).round(3).to_string(index=False))
