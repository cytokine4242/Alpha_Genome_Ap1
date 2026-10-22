#!/usr/bin/env python3
"""Four hypotheses about what else might differ between the AP-1 site classes, all tested on ALL 3,679 cardiomyocyte sites (not the scanned subset) so they are properly powered.
Footprint detection needs ATAC signal, so every comparison is also stratified by MEASURED accessibility; the unstratified version is reported alongside as the thing to distrust.
 A  repertoire: a cardiac-enhancer score (MEF2/GATA/NKX/TBX/TEAD) and a promoter-like score (KLF-SP/ETS/MYC/YY) per site, by class
 B  mass action: do the age-GAINED sites (score) carry fewer or more distant partners, as an affinity-driven site should?
 C  homotypic clustering: do Dynamic sites sit in peaks with more AP-1 motifs, making any single one dispensable?
 D  the model's error: is over-prediction (predicted high, measured low) explained by partner content?
Outputs Analysis/class_context/{tables,figures}.   python3 41_class_context_hypotheses.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, mannwhitneyu, spearmanr

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "class_context"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(41)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
CARDIAC = ["MEF2", "GATA", "NKX / homeobox cardiac", "TBX", "TEAD"]
PROMOTER = ["KLF / SP", "ETS", "MYC / MAX", "YY / Gli / Pax"]
NEAR = 100                                        # 'near' = within 100 bp, about where importance has halved


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def strat_test(frame, col, a_, b_, value, n_strata=4):
    """Mann-Whitney on the whole set, and a stratified version: the median difference within each quartile of measured accessibility,
    combined by van Elteren weighting, with a permutation p-value that shuffles the label inside each stratum."""
    w = frame[frame[col].isin([a_, b_])].dropna(subset=[value, "measured"])
    if min((w[col] == a_).sum(), (w[col] == b_).sum()) < 25:
        return None
    w = w.copy()
    w["q"] = pd.qcut(w.measured, n_strata, labels=False, duplicates="drop")
    obs = np.mean([w[(w.q == q) & (w[col] == a_)][value].median() - w[(w.q == q) & (w[col] == b_)][value].median() for q in sorted(w.q.unique())])
    null = []
    for k in range(2000):
        s = w.copy()
        s[col] = s.groupby("q")[col].transform(lambda v: rng.permutation(v.values))
        null.append(np.mean([s[(s.q == q) & (s[col] == a_)][value].median() - s[(s.q == q) & (s[col] == b_)][value].median() for q in sorted(s.q.unique())]))
    null = np.array(null)
    return {"n_a": int((w[col] == a_).sum()), "n_b": int((w[col] == b_).sum()),
            "median_a": w[w[col] == a_][value].median(), "median_b": w[w[col] == b_][value].median(),
            "p_unstratified": mannwhitneyu(w[w[col] == a_][value], w[w[col] == b_][value]).pvalue,
            "stratified_diff": obs, "p_stratified": (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)}


# ---------------------------------------------------------------- per-site context table
M = pd.read_csv(ROOT / "Analysis/distance_rules/tables/merged_partner_instances.csv.gz")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb

S = sl.join(scl[["age_score_class", "ex_score_class", "age_delta", "meas_age"]]).join(meas.rename("measured"))
near = M[M.dist <= NEAR]
for lab, fams in (("cardiac", CARDIAC), ("promoter", PROMOTER)):
    S[f"n_{lab}"] = near[near.family.isin(fams)].groupby("site_id").family.nunique().reindex(S.index).fillna(0)
S["n_families_near"] = near.groupby("site_id").family.nunique().reindex(S.index).fillna(0)
S["cardiac_share"] = np.where(S.n_families_near > 0, S.n_cardiac / S.n_families_near, np.nan)
S["has_mef2_near"] = S.index.isin(near[near.family == "MEF2"].site_id)
S["mef2_distance"] = M[M.family == "MEF2"].groupby("site_id").dist.min().reindex(S.index)
S["nkx_distance"] = M[M.family == "NKX / homeobox cardiac"].groupby("site_id").dist.min().reindex(S.index)
S["ctcf_distance"] = M[M.family == "CTCF / cohesin"].groupby("site_id").dist.min().reindex(S.index)
S["nearest_partner"] = M[M.family.isin(CARDIAC)].groupby("site_id").dist.min().reindex(S.index)
S["l_meas"] = np.log10(S.measured.clip(lower=1e-3))
S["overprediction"] = S.native_peak.rank(pct=True) - S.measured.rank(pct=True)
S.to_csv(OUT / "tables" / "per_site_context.csv")
print(f"{len(S):,} sites; median families within {NEAR} bp: {S.n_families_near.median():.0f} "
      f"(cardiac {S.n_cardiac.median():.0f}, promoter-like {S.n_promoter.median():.0f}); MEF2 within {NEAR} bp at {S.has_mef2_near.mean():.0%}")

# ---------------------------------------------------------------- matched pairs: predicted context AND measured accessibility, so footprint detection is matched by construction
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1)


def make_pairs(col, a_, b_):
    m_ = S[S[col].isin([a_, b_])].dropna(subset=MCOLS).copy()
    Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
    cases = m_.index[m_[col] == a_].tolist()
    avail = set(m_.index[m_[col] == b_])
    pairs = []
    for cid in rng.permutation(cases):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < 0.5:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


# ---------------------------------------------------------------- A/B. repertoire and spacing by class
CONTRASTS = [("Dynamic vs Stable", "stable_dynamic", "Dynamic", "Stable"),
             ("Age Up (score) vs Unchanged", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down (score) vs Unchanged", "age_score_class", "Age Down (score)", "Unchanged"),
             ("Exercise Up (score) vs Unchanged", "ex_score_class", "Exercise Up (score)", "Unchanged")]
MEASURES = [("n_families_near", "distinct families within 100 bp"), ("n_cardiac", "cardiac families within 100 bp"), ("n_promoter", "promoter-like families within 100 bp"),
            ("cardiac_share", "share of nearby families that are cardiac"), ("mef2_distance", "distance to the nearest MEF2 (bp)"),
            ("nkx_distance", "distance to the nearest NKX (bp)"), ("ctcf_distance", "distance to the nearest CTCF (bp)"),
            ("nearest_partner", "distance to the nearest cardiac partner (bp)"), ("motifs_in_peak", "AP-1 motifs in the same peak")]
from scipy.stats import wilcoxon
rows, PAIRS = [], {}
for lab, col, a_, b_ in CONTRASTS:
    P = make_pairs(col, a_, b_)
    PAIRS[lab] = P
    bal = {c: (S.loc[P.case, c].mean() - S.loc[P.control, c].mean()) / S[c].std() for c in MCOLS}
    for v, vlab in MEASURES:
        r = strat_test(S, col, a_, b_, v)
        x, y = S[v].reindex(P.case).values, S[v].reindex(P.control).values
        ok = ~(np.isnan(x) | np.isnan(y))
        pr = {"n_pairs": int(ok.sum()), "matched_median_diff": np.nan, "p_matched": np.nan}
        if ok.sum() >= 25 and np.any(x[ok] != y[ok]):
            pr = {"n_pairs": int(ok.sum()), "matched_median_diff": float(np.median(x[ok] - y[ok])), "p_matched": float(wilcoxon(x[ok] - y[ok]).pvalue)}
        if r:
            rows.append({"comparison": lab, "measure": vlab, **r, **pr, "worst_balance_sd": max(abs(v_) for v_ in bal.values())})
A = pd.DataFrame(rows)
A["fdr_matched"] = A.groupby("comparison").p_matched.transform(lambda p: bh(p.fillna(1)))
A["fdr_unstratified"] = A.groupby("comparison").p_unstratified.transform(lambda p: bh(p))
A["fdr_stratified"] = A.groupby("comparison").p_stratified.transform(lambda p: bh(p))
A.to_csv(OUT / "tables" / "a_context_by_class.csv", index=False)

# MEF2 presence specifically, with a Mantel-Haenszel test across accessibility quartiles
rows = []
for lab, col, a_, b_ in CONTRASTS:
    w = S[S[col].isin([a_, b_])].dropna(subset=["measured"]).copy()
    w["q"] = pd.qcut(w.measured, 4, labels=False, duplicates="drop")
    tabs = [np.array([[int((g[col] == a_).mul(g.has_mef2_near).sum()), int((g[col] == a_).mul(~g.has_mef2_near).sum())],
                      [int((g[col] == b_).mul(g.has_mef2_near).sum()), int((g[col] == b_).mul(~g.has_mef2_near).sum())]]) for _, g in w.groupby("q")]
    num = sum(tb[0, 0] * tb[1, 1] / tb.sum() for tb in tabs if tb.sum())
    den = sum(tb[0, 1] * tb[1, 0] / tb.sum() for tb in tabs if tb.sum())
    obs = sum(tb[0, 0] for tb in tabs)
    exp = sum(tb[0].sum() * tb[:, 0].sum() / tb.sum() for tb in tabs if tb.sum())
    var = sum(tb[0].sum() * tb[1].sum() * tb[:, 0].sum() * tb[:, 1].sum() / (tb.sum() ** 2 * (tb.sum() - 1)) for tb in tabs if tb.sum() > 1)
    P = PAIRS[lab]
    xa, xb = S.has_mef2_near.reindex(P.case).values, S.has_mef2_near.reindex(P.control).values
    n10, n01 = int((xa & ~xb).sum()), int((~xa & xb).sum())
    from scipy.stats import binomtest
    rows.append({"comparison": lab, "share_a": w[w[col] == a_].has_mef2_near.mean(), "share_b": w[w[col] == b_].has_mef2_near.mean(),
                 "n_pairs": len(P), "matched_share_a": xa.mean(), "matched_share_b": xb.mean(),
                 "p_matched_mcnemar": binomtest(n10, n10 + n01).pvalue if n10 + n01 else np.nan,
                 "p_unstratified": fisher_exact([[int(w[w[col] == a_].has_mef2_near.sum()), int((~w[w[col] == a_].has_mef2_near).sum())],
                                                 [int(w[w[col] == b_].has_mef2_near.sum()), int((~w[w[col] == b_].has_mef2_near).sum())]])[1],
                 "odds_ratio_stratified": num / den if den else np.nan,
                 "p_stratified": float(chi2.sf(max(abs(obs - exp) - 0.5, 0) ** 2 / var, 1)) if var > 0 else np.nan})
B = pd.DataFrame(rows)
B["fdr_stratified"] = bh(B.p_stratified)
B.to_csv(OUT / "tables" / "b_mef2_presence_by_class.csv", index=False)

# ---------------------------------------------------------------- C. is the partner deficit at Dynamic sites still detection?
# Decisive check: repeat inside the top quartile of measured accessibility, where footprint calling is least limited by signal.
from scipy.stats import binomtest as _bt
rows = []
for lab, col, a_, b_ in CONTRASTS[:1] + CONTRASTS[1:3]:
    for q_lab, sub in (("all sites", S), ("top quartile of measured accessibility", S[S.measured >= S.measured.quantile(0.75)]),
                       ("bottom half of measured accessibility", S[S.measured < S.measured.median()])):
        keep = sub.index
        m_ = S.loc[keep]
        m_ = m_[m_[col].isin([a_, b_])].dropna(subset=MCOLS)
        if min((m_[col] == a_).sum(), (m_[col] == b_).sum()) < 40:
            continue
        Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
        avail, pairs = set(m_.index[m_[col] == b_]), []
        for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
            cand = np.array(sorted(avail))
            if not len(cand):
                break
            d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
            j = int(np.argmin(d_))
            if d_[j] < 0.5:
                pairs.append((cid, cand[j]))
                avail.discard(cand[j])
        P = pd.DataFrame(pairs, columns=["case", "control"])
        if len(P) < 40:
            continue
        fa, fb = S.n_families_near.reindex(P.case).values, S.n_families_near.reindex(P.control).values
        ma, mb = S.has_mef2_near.reindex(P.case).values, S.has_mef2_near.reindex(P.control).values
        n10, n01 = int((ma & ~mb).sum()), int((~ma & mb).sum())
        rows.append({"comparison": lab, "subset": q_lab, "n_pairs": len(P), "median_families_a": np.median(fa), "median_families_b": np.median(fb),
                     "paired_diff_families": float(np.median(fa - fb)), "p_families": float(wilcoxon(fa - fb).pvalue) if np.any(fa != fb) else np.nan,
                     "mef2_share_a": ma.mean(), "mef2_share_b": mb.mean(), "p_mef2": _bt(n10, n10 + n01).pvalue if n10 + n01 else np.nan})
C = pd.DataFrame(rows)
C.to_csv(OUT / "tables" / "c_detection_sensitivity.csv", index=False)
print("== C is it detection? matched pairs inside accessibility strata\n", C.round(3).to_string())

# ---------------------------------------------------------------- D. does partner content explain the model's error?
rows = []
for v, vlab in (("n_families_near", "families within 100 bp"), ("n_cardiac", "cardiac families within 100 bp"), ("mef2_distance", "distance to nearest MEF2"), ("motifs_in_peak", "AP-1 motifs in the peak")):
    for grp, w in [("all sites", S)] + [(g, S[S.stable_dynamic == g]) for g in ("Stable", "Dynamic")]:
        ok = w.dropna(subset=[v, "overprediction"])
        r = spearmanr(ok[v], ok.overprediction)
        rp = spearmanr(ok[v], ok.native_peak.rank(pct=True))
        rm = spearmanr(ok[v], ok.measured.rank(pct=True))
        rows.append({"measure": vlab, "sites": grp, "n": len(ok), "rho_vs_overprediction": r.statistic, "p": r.pvalue,
                     "rho_vs_predicted": rp.statistic, "rho_vs_measured": rm.statistic})
D = pd.DataFrame(rows)
D["fdr"] = bh(D.p)
D.to_csv(OUT / "tables" / "d_partner_content_vs_error.csv", index=False)

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.2))
ax = axes[0]
sub = A[A.measure.isin([m[1] for m in MEASURES])].copy()
ordm = [m[1] for m in MEASURES]
for k, (lab, c) in enumerate((("Dynamic vs Stable", ORANGE), ("Age Up (score) vs Unchanged", BLUE), ("Age Down (score) vs Unchanged", AQUA))):
    w = sub[sub.comparison == lab].set_index("measure").reindex(ordm)
    z = (w.median_a - w.median_b) / np.array([S[m].std() for m, _ in MEASURES])
    ys = np.arange(len(ordm)) + (k - 1) * 0.24
    sig = (w.fdr_stratified < 0.05).values
    ax.scatter(z[~sig], ys[~sig], s=34, color=c, alpha=0.4)
    ax.scatter(z[sig], ys[sig], s=52, color=c, edgecolors=INK, lw=0.8, zorder=4)
    ax.scatter([], [], s=40, color=c, label=lab)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(np.arange(len(ordm)), ordm, fontsize=8)
ax.set_xlabel("Standardised difference (responsive class minus reference)\nring = survives stratification by measured accessibility, FDR < 0.05")
ax.set_title("Context of the AP-1 site, by class", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax = axes[1]
for g, c in (("Stable", BLUE), ("Dynamic", ORANGE)):
    w = S[S.stable_dynamic == g]
    h = w.n_families_near.value_counts(normalize=True).sort_index()
    ax.plot(h.index[:15], h.values[:15], color=c, lw=1.8, marker="o", ms=4, label=f"{g} (median {w.n_families_near.median():.0f})")
ax.set_xlabel(f"Distinct TF families footprinted within {NEAR} bp")
ax.set_ylabel("Share of sites")
ax.set_title("How crowded is the neighbourhood?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "figures" / "t1_context_by_class.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
print("== A context by class: matched pairs (predicted context + measured accessibility) vs stratified vs raw\n",
      A[["comparison", "measure", "n_a", "n_b", "median_a", "median_b", "fdr_unstratified", "fdr_stratified", "n_pairs", "matched_median_diff", "fdr_matched", "worst_balance_sd"]].round(3).to_string())
print("== B MEF2 presence within 100 bp\n", B.round(3).to_string())
print("== D partner content vs the model's error\n", D.round(3).to_string())
