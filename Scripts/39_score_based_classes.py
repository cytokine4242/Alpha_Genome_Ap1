#!/usr/bin/env python3
"""A second, parallel definition of the age and exercise classes: the CHANGE IN TOBIAS BINDING SCORE between conditions, instead of crossing the bound/unbound threshold.
The threshold definitions in Data/derived/cm_fos_sites.csv are left untouched; these are added alongside for comparison.
 A  build the score-based classes (age: MAS - YS, exercise: MAE - MAS) and compare them with the threshold classes
 B  which definition better matches the MEASURED change in snATAC accessibility between conditions?
 C  do the score-based classes differ in their response to the 21 mutation types (raw, and matched on predicted context)?
 D  continuous version: partial correlation of the score change with each mutation effect, adjusted for predicted accessibility
Outputs Analysis/score_classes/{tables,figures}.   python3 39_score_based_classes.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "score_classes"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(39)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
Q = 0.20                                            # top and bottom fifth of the score change


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_med(d, n=500):
    d = np.asarray(d, float)
    return np.percentile([np.median(rng.choice(d, len(d))) for _ in range(n)], [2.5, 97.5])


def rank_resid(y, x):
    """y with the rank-linear part of x removed: residual of rank(y) on rank(x)."""
    ok = ~(y.isna() | x.isna())
    ry, rx = y[ok].rank(pct=True), x[ok].rank(pct=True)
    b = np.polyfit(rx, ry, 1)
    out = pd.Series(np.nan, index=y.index)
    out[ok] = ry - np.polyval(b, rx)
    return out


# ---------------------------------------------------------------- inputs
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
prof = pd.read_csv(ROOT / "Analysis/dynamic/tables/mutation_type_profile_dynamic_vs_stable.csv")
TYPES, SHORT, GRP = prof.feature.tolist(), dict(zip(prof.feature, prof.mutation)), dict(zip(prof.feature, prof.group))
bycond = pd.read_csv(ROOT / "Katana/data/measured_by_condition.csv")
bycond = bycond[bycond.site_table == "cm"].set_index("site_id")

d = sl.join(bycond[["cm_ys", "cm_mas", "cm_mae"]])
d["age_delta"] = d.MAS_score - d.YS_score                       # score change, young sedentary -> aged sedentary
d["ex_delta"] = d.MAE_score - d.MAS_score                       # aged sedentary -> aged + exercise
d["l_170"], d["l_peak"] = np.log10(d.native_170.clip(lower=1e-3)), np.log10(d.native_peak.clip(lower=1e-3))
d["l_pw"], d["l_tss"] = np.log10(d.peak_width), np.log10(d.abs_tss + 1)
pc = 0.05 * d[["cm_ys", "cm_mas", "cm_mae"]].stack().median()
d["meas_age"] = np.log2((d.cm_mas + pc) / (d.cm_ys + pc))       # measured accessibility change, same two contrasts
d["meas_ex"] = np.log2((d.cm_mae + pc) / (d.cm_mas + pc))


def by_quantile(v, up, down):
    lo, hi = v.quantile(Q), v.quantile(1 - Q)
    return pd.Series(np.where(v >= hi, up, np.where(v <= lo, down, "Unchanged")), index=v.index).where(v.notna())


d["age_score_class"] = by_quantile(d.age_delta, "Age Up (score)", "Age Down (score)")
d["ex_score_class"] = by_quantile(d.ex_delta, "Exercise Up (score)", "Exercise Down (score)")
d["score_dynamic"] = np.where((d.age_delta.abs().rank(pct=True) > 1 - Q) | (d.ex_delta.abs().rank(pct=True) > 1 - Q), "Changing (score)", "Steady (score)")
d[["stable_dynamic", "age_class", "ex_class", "age_delta", "ex_delta", "age_score_class", "ex_score_class", "score_dynamic", "meas_age", "meas_ex"]].to_csv(OUT / "tables" / "score_based_classes.csv")

# ---------------------------------------------------------------- A. agreement with the threshold definitions
ct_age = pd.crosstab(d.age_class, d.age_score_class)
ct_ex = pd.crosstab(d.ex_class, d.ex_score_class)
ct_sd = pd.crosstab(d.stable_dynamic, d.score_dynamic)
ct_age.to_csv(OUT / "tables" / "a_crosstab_age.csv")
ct_ex.to_csv(OUT / "tables" / "a_crosstab_exercise.csv")
ct_sd.to_csv(OUT / "tables" / "a_crosstab_stable_dynamic.csv")
rows = []
for lab, thr, sco, up, dn in (("age", "age_class", "age_score_class", "Age Up", "Age Down"), ("exercise", "ex_class", "ex_score_class", "Exercise Up", "Exercise Down")):
    w = d.dropna(subset=[sco])
    rows.append({"contrast": lab, "n": len(w), "threshold_Up_also_score_Up": (w[w[thr] == up][sco] == f"{up} (score)").mean(),
                 "threshold_Down_also_score_Down": (w[w[thr] == dn][sco] == f"{dn} (score)").mean(),
                 "score_Up_that_crossed_threshold": (w[w[sco] == f"{up} (score)"][thr] == up).mean(),
                 "threshold_Stable_called_changing": (w[w[thr] == "Stable"][sco] != "Unchanged").mean(),
                 "rho_delta_vs_threshold_direction": spearmanr(w["age_delta" if lab == "age" else "ex_delta"], w[thr].map({up: 1, dn: -1}).fillna(0)).statistic})
A = pd.DataFrame(rows)
A.to_csv(OUT / "tables" / "a_agreement.csv", index=False)

# ---------------------------------------------------------------- B. which definition tracks the measured change better?
rows = []
for lab, delta, meas, thr, up, dn in (("age (MAS vs YS)", "age_delta", "meas_age", "age_class", "Age Up", "Age Down"),
                                      ("exercise (MAE vs MAS)", "ex_delta", "meas_ex", "ex_class", "Exercise Up", "Exercise Down")):
    w = d.dropna(subset=[meas, delta])
    r_sc = spearmanr(w[delta], w[meas])
    dirn = w[thr].map({up: 1, dn: -1}).fillna(0)                 # Stable and the 'neither condition' class both count as no change
    r_th = spearmanr(dirn, w[meas])
    bs_sc = [spearmanr(*w.sample(len(w), replace=True, random_state=k)[[delta, meas]].values.T).statistic for k in range(300)]
    bs_th = []
    for k in range(300):
        s_ = w.sample(len(w), replace=True, random_state=1000 + k)
        bs_th.append(spearmanr(s_[thr].map({up: 1, dn: -1}).fillna(0), s_[meas]).statistic)
    rows.append({"contrast": lab, "n": len(w), "rho_score_change_vs_measured": r_sc.statistic, "ci_lo_score": np.percentile(bs_sc, 2.5), "ci_hi_score": np.percentile(bs_sc, 97.5), "p_score": r_sc.pvalue,
                 "rho_threshold_class_vs_measured": r_th.statistic, "ci_lo_threshold": np.percentile(bs_th, 2.5), "ci_hi_threshold": np.percentile(bs_th, 97.5), "p_threshold": r_th.pvalue})
B = pd.DataFrame(rows)
B.to_csv(OUT / "tables" / "b_vs_measured_change.csv", index=False)

# measured change by class, both definitions
rows = []
for defn, col, groups in (("threshold", "age_class", ["Age Up", "Age Down", "Stable"]), ("score", "age_score_class", ["Age Up (score)", "Age Down (score)", "Unchanged"]),
                          ("threshold", "ex_class", ["Exercise Up", "Exercise Down", "Stable"]), ("score", "ex_score_class", ["Exercise Up (score)", "Exercise Down (score)", "Unchanged"])):
    meas = "meas_age" if "age" in col else "meas_ex"
    for g in groups:
        w = d[d[col] == g][meas].dropna()
        if len(w) < 10:
            continue
        lo, hi = boot_med(w.values, 300)
        rows.append({"contrast": "age" if "age" in col else "exercise", "definition": defn, "class": g, "n": len(w), "median_measured_log2": w.median(), "ci_lo": lo, "ci_hi": hi})
BM = pd.DataFrame(rows)
BM.to_csv(OUT / "tables" / "b_measured_change_by_class.csv", index=False)

# ---------------------------------------------------------------- C. mutation responses under the score-based classes
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss"]


def matched_pairs2(frame, is_case, cols):
    """1:1 nearest-neighbour matching of cases to controls on `cols` (caliper 0.5 SD)."""
    m_ = frame.dropna(subset=cols).copy()
    Z = (m_[cols] - m_[cols].mean()) / m_[cols].std()
    case_ids = m_.index[is_case.reindex(m_.index).fillna(False)].tolist()
    avail = set(m_.index) - set(case_ids)
    pairs = []
    for cid in rng.permutation(case_ids):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        dist = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(dist))
        if dist[j] < 0.5:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


def matched_pairs(frame, is_case):
    return matched_pairs2(frame, is_case, MCOLS)


CONTRASTS = [("Age Up (score)", "age_score_class", "Unchanged"), ("Age Down (score)", "age_score_class", "Unchanged"),
             ("Exercise Up (score)", "ex_score_class", "Unchanged"), ("Exercise Down (score)", "ex_score_class", "Unchanged"),
             ("Changing (score)", "score_dynamic", "Steady (score)")]
rows = []
for case, col, ctrl in CONTRASTS:
    frame = d[d[col].isin([case, ctrl])]
    P = matched_pairs(frame, frame[col] == case)
    for c in TYPES:
        a, b_ = FD[c].reindex(P.case).values, FD[c].reindex(P.control).values
        ok = ~(np.isnan(a) | np.isnan(b_))
        if ok.sum() < 30:
            continue
        diff = a[ok] - b_[ok]
        raw_case, raw_ctrl = FD[c].reindex(frame.index[frame[col] == case]), FD[c].reindex(frame.index[frame[col] == ctrl])
        lo, hi = boot_med(diff, 300)
        rows.append({"contrast": case, "group": GRP[c], "mutation": SHORT[c], "feature": c, "n_pairs": int(ok.sum()),
                     "raw_case_median": raw_case.median(), "raw_control_median": raw_ctrl.median(), "raw_diff": raw_case.median() - raw_ctrl.median(),
                     "p_raw": mannwhitneyu(raw_case.dropna(), raw_ctrl.dropna()).pvalue,
                     "matched_diff": np.median(diff), "ci_lo": lo, "ci_hi": hi, "p_matched": wilcoxon(diff).pvalue})
C = pd.DataFrame(rows)
C["fdr_raw"] = C.groupby("contrast").p_raw.transform(lambda p: bh(p))
C["fdr_matched"] = C.groupby("contrast").p_matched.transform(lambda p: bh(p))
C.to_csv(OUT / "tables" / "c_mutation_response_by_score_class.csv", index=False)

# ---------------------------------------------------------------- D. continuous: score change vs each mutation effect
rows = []
for lab, delta in (("age", "age_delta"), ("exercise", "ex_delta")):
    adj_delta = rank_resid(d[delta], d.l_170)
    for c in TYPES:
        y = FD[c].reindex(d.index)
        r = spearmanr(d[delta], y, nan_policy="omit")
        ra = spearmanr(adj_delta, rank_resid(y, d.l_170), nan_policy="omit")
        rows.append({"contrast": lab, "group": GRP[c], "mutation": SHORT[c], "n": int((~(d[delta].isna() | y.isna())).sum()),
                     "rho_raw": r.statistic, "p_raw": r.pvalue, "rho_adjusted": ra.statistic, "p_adjusted": ra.pvalue})
D = pd.DataFrame(rows)
D["fdr_adjusted"] = D.groupby("contrast").p_adjusted.transform(lambda p: bh(p))
D.to_csv(OUT / "tables" / "d_continuous_correlations.csv", index=False)

# ---------------------------------------------------------------- E. control: add the absolute footprint score to the matching
d["abs_score"] = d[["YS_score", "MAS_score", "MAE_score"]].mean(1)
d["l_score"] = np.log10(d.abs_score)
rows = []
for case, col, ctrl in CONTRASTS:
    frame = d[d[col].isin([case, ctrl])]
    for cols, lab in ((MCOLS, "context only"), (MCOLS + ["l_score"], "context + absolute footprint score")):
        globals()["MCOLS_USE"] = cols
        P = matched_pairs2(frame, frame[col] == case, cols)
        for c in TYPES:
            a, b_ = FD[c].reindex(P.case).values, FD[c].reindex(P.control).values
            ok = ~(np.isnan(a) | np.isnan(b_))
            if ok.sum() < 30:
                continue
            diff = a[ok] - b_[ok]
            rows.append({"contrast": case, "matching": lab, "group": GRP[c], "mutation": SHORT[c], "n_pairs": int(ok.sum()),
                         "matched_diff": np.median(diff), "p_matched": wilcoxon(diff).pvalue})
E = pd.DataFrame(rows)
E["fdr"] = E.groupby(["contrast", "matching"]).p_matched.transform(lambda p: bh(p))
E.to_csv(OUT / "tables" / "e_matching_with_score.csv", index=False)
bal = []
for case, col, ctrl in CONTRASTS:
    frame = d[d[col].isin([case, ctrl])]
    for cols, lab in ((MCOLS, "context only"), (MCOLS + ["l_score"], "context + absolute footprint score")):
        P = matched_pairs2(frame, frame[col] == case, cols)
        for v in MCOLS + ["l_score"]:
            bal.append({"contrast": case, "matching": lab, "covariate": v,
                        "std_diff_after": (d.loc[P.case, v].mean() - d.loc[P.control, v].mean()) / d[v].std()})
pd.DataFrame(bal).to_csv(OUT / "tables" / "e_matching_balance.csv", index=False)

# ---------------------------------------------------------------- F. is the difference specific to the motif?
SPEC = {"motif share of the local dependence": ("replace motif, read over motif", "centred shuffle of 680 bp, read over motif"),
        "core share of the motif": ("replace core, read over motif", "replace motif, read over motif")}
rows = []
for case, col, ctrl in CONTRASTS:
    frame = d[d[col].isin([case, ctrl])]
    P = matched_pairs2(frame, frame[col] == case, MCOLS)
    for lab, (num, den) in SPEC.items():
        r = (FD[num] / FD[den]).where(FD[den] < -0.5)
        a, b_ = r.reindex(P.case).values, r.reindex(P.control).values
        ok = ~(np.isnan(a) | np.isnan(b_))
        if ok.sum() < 30:
            continue
        diff = a[ok] - b_[ok]
        rows.append({"contrast": case, "measure": lab, "n_pairs": int(ok.sum()), "median_case": np.nanmedian(a[ok]), "median_control": np.nanmedian(b_[ok]),
                     "paired_diff": np.median(diff), "p": wilcoxon(diff).pvalue})
    for c, lab in (("170 bp shuffle ±300 bp away, read over motif", "far control, ±300 bp"), ("170 bp shuffle ±600 bp away, read over motif", "far control, ±600 bp")):
        a, b_ = FD[c].reindex(P.case).values, FD[c].reindex(P.control).values
        ok = ~(np.isnan(a) | np.isnan(b_))
        diff = a[ok] - b_[ok]
        rows.append({"contrast": case, "measure": lab, "n_pairs": int(ok.sum()), "median_case": np.median(a[ok]), "median_control": np.median(b_[ok]),
                     "paired_diff": np.median(diff), "p": wilcoxon(diff).pvalue})
F = pd.DataFrame(rows)
F["fdr"] = bh(F.p)
F.to_csv(OUT / "tables" / "f_motif_specificity.csv", index=False)

# ---------------------------------------------------------------- G. do Age Up (score) sites carry a better AP-1 core? (scan subset, which records the core sequence)
import glob as _glob
from scipy.stats import fisher_exact
G = None
_f = sorted(_glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))
if _f:
    s_ = pd.concat([pd.read_csv(x, usecols=["site_id", "core_found", "trackset"]) for x in _f], ignore_index=True)
    s_ = s_[s_.trackset == "CM"].drop_duplicates("site_id").set_index("site_id").core_found
    w = d.join(s_.rename("core_found")).dropna(subset=["core_found"])
    rows = []
    for col, up, dn, ctrl in (("age_score_class", "Age Up (score)", "Age Down (score)", "Unchanged"), ("ex_score_class", "Exercise Up (score)", "Exercise Down (score)", "Unchanged")):
        for g in (up, dn):
            a, b_ = w[w[col] == g].core_found, w[w[col] == ctrl].core_found
            if min(len(a), len(b_)) < 10:
                continue
            rows.append({"class": g, "n": len(a), "share_with_exact_core": a.mean(), "n_control": len(b_), "share_control": b_.mean(),
                         "p_fisher": fisher_exact([[a.sum(), (~a.astype(bool)).sum()], [b_.sum(), (~b_.astype(bool)).sum()]])[1]})
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "tables" / "g_core_sequence_by_score_class.csv", index=False)

# ---------------------------------------------------------------- H. why the contrasts behave differently: a site can only lose score if it had score to lose
prof_rows = []
for col in ("age_score_class", "ex_score_class"):
    for g, w in d.groupby(col):
        prof_rows.append({"contrast": "age" if col.startswith("age") else "exercise", "class": g, "n": len(w),
                          "median_absolute_score": w.abs_score.median(), "median_start_score": w[("YS_score" if col.startswith("age") else "MAS_score")].median(),
                          "median_end_score": w[("MAS_score" if col.startswith("age") else "MAE_score")].median(),
                          "median_native_170": w.native_170.median(), "median_measured_change": w[("meas_age" if col.startswith("age") else "meas_ex")].median()})
HS = pd.DataFrame(prof_rows)
HS.to_csv(OUT / "tables" / "h_score_profile_by_class.csv", index=False)
print("== H score profile by class (regression to the mean)\n", HS.round(3).to_string())

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axes[0]
for lab, delta, c in (("Age (MAS − YS)", "age_delta", ORANGE), ("Exercise (MAE − MAS)", "ex_delta", BLUE)):
    v = d[delta].dropna()
    ax.hist(v, bins=80, histtype="step", color=c, lw=1.6, label=f"{lab}  (n = {len(v):,})")
for q_, s_ in ((Q, "-"), (1 - Q, "-")):
    ax.axvline(d.age_delta.quantile(q_), color=ORANGE, lw=0.9, ls=":")
    ax.axvline(d.ex_delta.quantile(q_), color=BLUE, lw=0.9, ls=":")
ax.set_xlabel("Change in TOBIAS binding score between conditions")
ax.set_ylabel("Sites")
ax.set_title(f"Score change, with the {Q:.0%} / {1 - Q:.0%} cut points (dotted)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
ax = axes[1]
w = BM[BM.contrast == "age"]
ys = np.arange(len(w))
cols = [ORANGE if "Up" in g else BLUE if "Down" in g else GREY for g in w["class"]]
ax.errorbar(w.median_measured_log2, ys, xerr=[w.median_measured_log2 - w.ci_lo, w.ci_hi - w.median_measured_log2], fmt="none", ecolor=GREY, elinewidth=1.3)
ax.scatter(w.median_measured_log2, ys, s=46, color=cols, zorder=3)
ax.set_yticks(ys, [f"{g}  ({dd}, n = {n})" for g, dd, n in zip(w["class"], w.definition, w.n)], fontsize=8)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_xlabel("Measured snATAC change, log2(aged sedentary / young sedentary)")
ax.set_title("Does each class show the measured change it is named for?", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "q1_score_classes.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.4, 6.4))
w = C[C.contrast == "Changing (score)"].reset_index(drop=True)[::-1]
ys = np.arange(len(w))
sig = (w.fdr_matched < 0.05).values
ax.errorbar(w.matched_diff, ys, xerr=[w.matched_diff - w.ci_lo, w.ci_hi - w.matched_diff], fmt="none", ecolor=GREY, elinewidth=1.3)
ax.scatter(w.raw_diff, ys, s=40, marker="D", facecolors="none", edgecolors=ORANGE, lw=1.2, label="raw difference")
ax.scatter(w.matched_diff[~sig], ys[~sig], s=38, color=GREY, zorder=3, label="matched, not significant")
ax.scatter(w.matched_diff[sig], ys[sig], s=46, color=BLUE, zorder=4, label="matched, FDR < 0.05")
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(ys, w.mutation, fontsize=8)
ax.set_xlabel("Changing (score) minus Steady (score), log2FC\nraw (open diamonds) and matched on predicted accessibility, peak width and TSS distance (filled)")
ax.set_title("Mutation responses under the score-based definition", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "q2_mutation_response_score_classes.png", bbox_inches="tight")
plt.close(fig)

# q3: what the score definition finds that the threshold definition does not
thr_prof = prof.set_index("feature")
fig, ax = plt.subplots(figsize=(8.6, 6.4))
w = C[C.contrast == "Age Up (score)"].set_index("feature").reindex(TYPES)
ys = np.arange(len(TYPES))[::-1]
ax.scatter(thr_prof.loc[TYPES, "matched_dynamic_median"] - thr_prof.loc[TYPES, "matched_stable_median"], ys, s=40, marker="s", color=GREY, label="threshold: Dynamic vs Stable (none significant)")
sig = (w.fdr_matched < 0.05).values
ax.errorbar(w.matched_diff, ys, xerr=[w.matched_diff - w.ci_lo, w.ci_hi - w.matched_diff], fmt="none", ecolor=AQUA, elinewidth=1.2, alpha=0.7)
ax.scatter(w.matched_diff[~sig], ys[~sig], s=40, color=AQUA, alpha=0.45, label="score: Age Up vs Unchanged, ns")
ax.scatter(w.matched_diff[sig], ys[sig], s=52, color=ORANGE, zorder=4, label="score: Age Up vs Unchanged, FDR < 0.05")
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(ys, [SHORT[c] for c in TYPES], fontsize=8)
ax.set_xlabel("Matched difference in effect, log2FC (negative = the responsive class loses more)")
ax.set_title("What the score-based definition finds that the threshold definition misses", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "q3_definitions_compared.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
print("== A agreement\n", A.round(3).to_string())
print("\nage: threshold (rows) vs score (columns)\n", ct_age.to_string(), "\n\nStable/Dynamic vs score\n", ct_sd.to_string())
print("== B vs measured change\n", B.round(3).to_string())
print("== B measured change by class\n", BM.round(3).to_string())
print("== C significant after matching\n", C[C.fdr_matched < 0.05].round(3).to_string())
print("== C summary per contrast\n", C.groupby("contrast").agg(n_types=("feature", "size"), n_sig_raw=("fdr_raw", lambda v: (v < 0.05).sum()), n_sig_matched=("fdr_matched", lambda v: (v < 0.05).sum()),
                                                               max_abs_matched=("matched_diff", lambda v: v.abs().max()), n_pairs=("n_pairs", "max")).round(3).to_string())
print("== E matching with the absolute score\n", E[(E.fdr < 0.05)].pivot_table(index=["contrast", "mutation"], columns="matching", values="matched_diff").round(3).to_string())
print("== E n significant per contrast and matching\n", E.groupby(["contrast", "matching"]).fdr.apply(lambda v: (v < 0.05).sum()).to_string())
print("== F motif specificity\n", F.round(3).to_string())
print("== G core sequence\n", "" if G is None else G.round(3).to_string())
print("== D strongest adjusted correlations\n", D.reindex(D.rho_adjusted.abs().sort_values(ascending=False).index).head(8).round(3).to_string())
