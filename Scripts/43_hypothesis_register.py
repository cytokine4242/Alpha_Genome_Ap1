#!/usr/bin/env python3
"""A register of hypotheses about what distinguishes the AP-1 site classes: each one stated, tested, and recorded with its verdict, whether or not it held.
Every comparison uses 1:1 matching on predicted accessibility (motif and peak), peak width, TSS distance and MEASURED accessibility, unless stated otherwise.
New inputs: Data/derived/core_quality.csv (the AP-1 core of all 3,679 sites, from a scan of mm10) and Data/derived/motif_counts.csv (consensus matches nearby).
Outputs Analysis/hypotheses/{tables,figures}, including hypothesis_register.csv.   python3 43_hypothesis_register.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, mannwhitneyu, spearmanr, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "hypotheses"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(43)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY, RED = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993", "#b3261e"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]
REG = []


def record(hid, hypothesis, test, n, result, verdict, detail=""):
    REG.append({"id": hid, "hypothesis": hypothesis, "test": test, "n": n, "result": result, "verdict": verdict, "detail": detail})
    print(f"[{hid}] {verdict:<12} {hypothesis}\n     {result}")


# ---------------------------------------------------------------- data
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
core = pd.read_csv(ROOT / "Data/derived/core_quality.csv").set_index("site_id")
mc = pd.read_csv(ROOT / "Data/derived/motif_counts.csv").set_index("site_id")
PWM = ROOT / "Data/derived/pwm_counts.csv"
if PWM.exists():
    mc = mc.join(pd.read_csv(PWM).set_index("site_id"))
MEF2_MOTIF = "pwm_MEF2_100" if "pwm_MEF2_100" in mc.columns else "MEF2_100"
AP1_MOTIF = "pwm_AP1_100" if "pwm_AP1_100" in mc.columns else "AP1_100"
print(f"motif presence from: {MEF2_MOTIF} (PWM scan)" if MEF2_MOTIF.startswith("pwm") else f"motif presence from: {MEF2_MOTIF} (strict consensus; PWM scan not yet available)")
FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id")
S = sl.join(scl[["age_score_class", "ex_score_class"]]).join(meas.rename("measured")).join(core).join(mc).join(sites[["in_YS", "in_MAS", "in_MAE"]], rsuffix="_site")
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))
M = pd.read_csv(ROOT / "Analysis/distance_rules/tables/merged_partner_instances.csv.gz")


def make_pairs(frame, col, a_, b_, cols=MCOLS, caliper=0.5):
    m_ = frame[frame[col].isin([a_, b_])].dropna(subset=cols).copy()
    Z = (m_[cols] - m_[cols].mean()) / m_[cols].std()
    avail, pairs = set(m_.index[m_[col] == b_]), []
    for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < caliper:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


def paired_binary(frame, P, v):
    a, b = frame[v].reindex(P.case).values.astype(float), frame[v].reindex(P.control).values.astype(float)
    ok = ~(np.isnan(a) | np.isnan(b))
    a, b = a[ok].astype(bool), b[ok].astype(bool)
    n10, n01 = int((a & ~b).sum()), int((~a & b).sum())
    return a.mean(), b.mean(), (binomtest(n10, n10 + n01).pvalue if n10 + n01 else np.nan), int(ok.sum())


def paired_num(frame, P, v):
    a, b = frame[v].reindex(P.case).values.astype(float), frame[v].reindex(P.control).values.astype(float)
    ok = ~(np.isnan(a) | np.isnan(b))
    d = a[ok] - b[ok]
    return np.median(a[ok]), np.median(b[ok]), (wilcoxon(d).pvalue if np.any(d != 0) else np.nan), int(ok.sum()), float(np.mean(d))


# ================================================================ H1. do age-gained sites carry better AP-1 cores? (genome-wide, was 78 sites)
P = make_pairs(S, "age_score_class", "Age Up (score)", "Unchanged")
ca, cb, p1, n1 = paired_binary(S, P, "exact_core")
mm_a, mm_b, p1b, _, mm_mean = paired_num(S, P, "core_mismatches")
record("H1", "Age Up (delta score): sites in the top fifth of the CHANGE in AP-1 binding score with age carry better AP-1 cores.",
       "Exact TGA[C/G]TCA core, matched pairs, genome-wide scan of mm10", n1,
       f"{ca:.1%} of Age Up (delta score) sites vs {cb:.1%} of matched Unchanged have an exact core (McNemar p = {p1:.2g}); mean difference in mismatch count {mm_mean:+.3f} (p = {p1b:.2g}).",
       "SUPPORTED" if p1 < 0.05 and ca > cb else ("NOT SUPPORTED" if p1 >= 0.05 else "CONTRADICTED"),
       "Earlier evidence was 78 scanned sites (90% vs 77%); this is the whole set.")

Pd = make_pairs(S, "age_score_class", "Age Down (score)", "Unchanged")
da, db, p1c, n1c = paired_binary(S, Pd, "exact_core")
record("H1b", "Age Down (delta score): the same holds for sites LOSING binding score, i.e. it is about change, not direction.",
       "Exact core, matched pairs", n1c,
       f"{da:.1%} of Age Down (delta score) sites vs {db:.1%} of matched Unchanged (McNemar p = {p1c:.2g}).",
       "NOT SUPPORTED" if p1c >= 0.05 else "SUPPORTED", "If this were significant too, the H1 result would be about lability, not gain.")

# ================================================================ H2. do Dynamic sites carry worse cores?
P2 = make_pairs(S, "stable_dynamic", "Dynamic", "Stable")
ea, eb, p2, n2 = paired_binary(S, P2, "exact_core")
record("H2", "Dynamic (threshold): sites whose binding CALL changes in any condition carry worse AP-1 cores than context-matched Stable sites.",
       "Exact core, matched pairs", n2,
       f"{ea:.1%} vs {eb:.1%} have an exact core (McNemar p = {p2:.2g}).",
       "SUPPORTED" if p2 < 0.05 and ea < eb else "NOT SUPPORTED",
       "Tests whether the threshold classes have a sequence signature of their own.")

# ================================================================ H3. is the MEF2 footprint deficit at Dynamic sites detection, or occupancy?
near_mef2 = set(M[(M.family == "MEF2") & (M.dist <= 100)].site_id)
S["mef2_footprint"] = S.index.isin(near_mef2)
S["mef2_motif"] = S[MEF2_MOTIF] > 0
sub = S[S.mef2_motif & (S.measured >= S.measured.quantile(0.5))]
print(f"   sites with a MEF2 motif: {int(S.mef2_motif.sum()):,}; of those, above median accessibility: {len(sub):,}")            # motif present AND enough signal to call a footprint
P3 = make_pairs(sub, "stable_dynamic", "Dynamic", "Stable")
fa, fb, p3, n3 = paired_binary(sub, P3, "mef2_footprint")
P3all = make_pairs(S[S.mef2_motif], "stable_dynamic", "Dynamic", "Stable")
fa2, fb2, p3b, n3b = paired_binary(S[S.mef2_motif], P3all, "mef2_footprint")
record("H3", "Where a MEF2 motif exists in the DNA, it is less often OCCUPIED at Dynamic sites (not merely less often detected).",
       "MEF2 footprint called, among sites with a MEF2 consensus match within 100 bp and measured accessibility above the median; matched pairs", n3,
       f"{fa:.1%} of Dynamic vs {fb:.1%} of matched Stable carry a MEF2 footprint (McNemar p = {p3:.2g}); over all motif-bearing sites {fa2:.1%} vs {fb2:.1%} (p = {p3b:.2g}).",
       "SUPPORTED" if p3 < 0.05 and fa < fb else "NOT SUPPORTED",
       "Restricting to signal-rich sites removes most of the detection bias; a surviving deficit points to occupancy.")

# ================================================================ H4. does partner occupancy move with the AP-1 site across conditions?
tob = pd.read_csv(ROOT / "Data/derived/tobias_cm_all_motifs.csv.gz", usecols=["seqnames", "start", "end", "in_YS", "in_MAS", "in_MAE", "name"], dtype={"name": "category"})
import re
tob = tob[tob.name.astype(str).str.match(r"^(Mef2|Borcs8)", case=False, na=False)]
idx = {c: g.sort_values("start").reset_index(drop=True) for c, g in tob.groupby("seqnames", observed=True)}
rows = []
for sid, r in sites.drop(columns=[c for c in ("stable_dynamic",) if c in sites]).join(S[["stable_dynamic"]], how="inner").iterrows():
    g = idx.get(r.seqnames)
    if g is None:
        continue
    lo, hi = r.motif_mid - 100, r.motif_mid + 100
    w = g.iloc[max(np.searchsorted(g.start.values, lo - 60), 0):np.searchsorted(g.start.values, hi)]
    w = w[w.end.values > lo]
    if not len(w):
        continue
    rows.append({"site_id": sid, "stable_dynamic": r.stable_dynamic,
                 "ap1": (int(r.in_YS), int(r.in_MAS), int(r.in_MAE)),
                 "mef2": (int(w.in_YS.any()), int(w.in_MAS.any()), int(w.in_MAE.any()))})
CO = pd.DataFrame(rows)
CO["ap1_dynamic"] = CO.ap1.map(lambda t: len(set(t)) > 1)
CO["mef2_dynamic"] = CO.mef2.map(lambda t: len(set(t)) > 1)
CO["same_pattern"] = CO.ap1 == CO.mef2
obs = CO[CO.ap1_dynamic].mef2_dynamic.mean()
exp = CO[~CO.ap1_dynamic].mef2_dynamic.mean()
tab = [[int(CO[CO.ap1_dynamic].mef2_dynamic.sum()), int((~CO[CO.ap1_dynamic].mef2_dynamic).sum())],
       [int(CO[~CO.ap1_dynamic].mef2_dynamic.sum()), int((~CO[~CO.ap1_dynamic].mef2_dynamic).sum())]]
orv, p4 = fisher_exact(tab)
# do they change in the SAME direction, beyond chance? compare with MEF2 patterns shuffled across sites
null = [np.mean(CO.ap1.values == rng.permutation(CO.mef2.values)) for _ in range(500)]
p4b = (np.sum(np.array(null) >= CO.same_pattern.mean()) + 1) / 501
CO.to_csv(OUT / "tables" / "h4_ap1_mef2_cooccupancy.csv", index=False)
record("H4", "A condition-variable (threshold) AP-1 site sits beside a MEF2 footprint that is itself condition-variable, i.e. they are co-regulated.",
       "Among AP-1 sites with a MEF2 footprint within 100 bp: is MEF2 condition-variable more often when AP-1 is? (Fisher); and do the two share the exact YS/MAS/MAE pattern more than chance (permutation)?", len(CO),
       f"MEF2 is condition-variable at {obs:.1%} of condition-variable AP-1 sites vs {exp:.1%} of constitutive ones (OR {orv:.2f}, p = {p4:.2g}); identical three-condition pattern in {CO.same_pattern.mean():.1%} of loci vs {np.mean(null):.1%} expected (p = {p4b:.3f}).",
       "SUPPORTED" if (p4 < 0.05 and orv > 1) or p4b < 0.05 else "NOT SUPPORTED",
       "Co-variation across conditions cannot be produced by a static detection bias.")

# ================================================================ H5. homotypic AP-1 clustering
a5, b5, p5, n5, mean5 = paired_num(S, P2, AP1_MOTIF)
record("H5", "Dynamic (threshold) sites sit among more AP-1 motifs, so no single one is essential.",
       "Consensus TGA[C/G]TCA matches within 100 bp (sequence), matched pairs", n5,
       f"median {a5:.1f} vs {b5:.1f} matches, mean difference {mean5:+.2f} (paired Wilcoxon p = {p5:.2g}).",
       "SUPPORTED" if p5 < 0.05 and mean5 > 0 else ("CONTRADICTED" if p5 < 0.05 else "NOT SUPPORTED"),
       "A significant difference in the opposite direction counts as contradicted, not merely unsupported.")

# ================================================================ H6. positive control: does core quality predict the model's motif dependence?
d6 = S.join(FD[["replace motif, read over motif", "shuffle motif, read over motif"]])
r6 = mannwhitneyu(d6[d6.exact_core == 1]["replace motif, read over motif"].dropna(), d6[d6.exact_core == 0]["replace motif, read over motif"].dropna())
m1, m0 = d6[d6.exact_core == 1]["replace motif, read over motif"].median(), d6[d6.exact_core == 0]["replace motif, read over motif"].median()
rho6 = spearmanr(d6.core_mismatches, d6["replace motif, read over motif"], nan_policy="omit")
record("H6", "Control: AlphaGenome's motif dependence tracks the quality of the AP-1 core in the DNA.",
       "Effect of replacing the motif, exact core vs not; and correlation with mismatch count", int(d6["replace motif, read over motif"].notna().sum()),
       f"exact core {m1:.2f} vs imperfect {m0:.2f} log2 (p = {r6.pvalue:.2g}); rho with mismatches {rho6.statistic:+.3f} (p = {rho6.pvalue:.2g}).",
       "SUPPORTED" if r6.pvalue < 0.05 and m1 < m0 else "NOT SUPPORTED",
       "If this failed, the sequence scan and the model would be measuring different things.")

# ================================================================ H7. is core quality the MECHANISM behind the Age Up motif sensitivity?
TYPES = ["shuffle motif, read over motif", "replace motif, read over motif", "dinucleotide-shuffle motif, read over motif", "centred shuffle of 7 bp, read over motif"]
d7 = S.join(FD[TYPES])
rows = []
for lab, cols in (("context only", MCOLS), ("context + core quality", MCOLS + ["core_mismatches"])):
    P7 = make_pairs(d7, "age_score_class", "Age Up (score)", "Unchanged", cols)
    for v in TYPES:
        a, b, p_, n_, mean_d = paired_num(d7, P7, v)
        rows.append({"matching": lab, "mutation": v, "n_pairs": n_, "median_case": a, "median_control": b, "mean_diff": mean_d, "p": p_})
H7 = pd.DataFrame(rows)
H7.to_csv(OUT / "tables" / "h7_core_mediation.csv", index=False)
n_sig = {k: int((w.p < 0.05).sum()) for k, w in H7.groupby("matching")}
record("H7", "The extra motif sensitivity of Age Up (delta score) sites is explained by their better cores (core quality is the mechanism, not a correlate).",
       "Repeat the Age Up vs Unchanged matched comparison on the four motif-level edits, with and without core mismatch count added to the matching", int(H7.n_pairs.max()),
       f"significant motif edits: {n_sig.get('context only', 0)} of {len(TYPES)} matching on context only, {n_sig.get('context + core quality', 0)} of {len(TYPES)} after also matching on core quality.",
       "SUPPORTED" if n_sig.get("context + core quality", 0) < n_sig.get("context only", 0) else "NOT SUPPORTED",
       "A drop to zero means core quality carries the effect; no change means the two are separate.")

# ================================================================ H4b. is the AP-1/MEF2 co-variation specific, or just shared peak signal?
# Control: the same concordance between AP-1 and every OTHER factor footprinted in the same window. If a peak's signal drives everything,
# AP-1 should agree with all of its neighbours equally well, and MEF2 should not stand out.
tob2 = pd.read_csv(ROOT / "Data/derived/tobias_cm_all_motifs.csv.gz", usecols=["seqnames", "start", "end", "in_YS", "in_MAS", "in_MAE", "name"], dtype={"name": "category"})
fam_of = {}
import re as _re
FAMPAT = [("MEF2", r"^Mef2|^Borcs8"), ("GATA", r"^Gata"), ("NFI", r"^Nfi[abcx]"), ("TEAD", r"^Tead"), ("KLF / SP", r"^Klf|^Sp[0-9]"), ("CTCF", r"^Ctcf")]
def _fam(n):
    for lab, rx in FAMPAT:
        if _re.match(rx, str(n), _re.I):
            return lab
    return "other"
tob2["family"] = tob2.name.map(_fam)
idx2 = {c: g.sort_values("start").reset_index(drop=True) for c, g in tob2.groupby("seqnames", observed=True)}
rows = []
base = sites.drop(columns=[c for c in ("stable_dynamic",) if c in sites]).join(S[["stable_dynamic"]], how="inner")
for sid, r in base.iterrows():
    g = idx2.get(r.seqnames)
    if g is None:
        continue
    lo, hi = r.motif_mid - 100, r.motif_mid + 100
    w = g.iloc[max(np.searchsorted(g.start.values, lo - 60), 0):np.searchsorted(g.start.values, hi)]
    w = w[w.end.values > lo]
    if not len(w):
        continue
    ap1 = (int(r.in_YS), int(r.in_MAS), int(r.in_MAE))
    for fam_, ww in w.groupby("family", observed=True):
        pat = (int(ww.in_YS.any()), int(ww.in_MAS.any()), int(ww.in_MAE.any()))
        rows.append({"site_id": sid, "family": fam_, "same_pattern": pat == ap1, "partner_variable": len(set(pat)) > 1, "ap1_variable": len(set(ap1)) > 1})
CC = pd.DataFrame(rows)
CC.to_csv(OUT / "tables" / "h4b_concordance_by_family.csv", index=False)
cc = CC.groupby("family").agg(n=("same_pattern", "size"), concordance=("same_pattern", "mean"),
                              var_given_ap1_var=("partner_variable", lambda v: np.nan)).reset_index()
cc["concordance_when_ap1_variable"] = [CC[(CC.family == f) & CC.ap1_variable].same_pattern.mean() for f in cc.family]
cc["partner_variable_when_ap1_variable"] = [CC[(CC.family == f) & CC.ap1_variable].partner_variable.mean() for f in cc.family]
cc["partner_variable_when_ap1_constant"] = [CC[(CC.family == f) & ~CC.ap1_variable].partner_variable.mean() for f in cc.family]
cc = cc.drop(columns=["var_given_ap1_var"]).sort_values("concordance", ascending=False)
cc.to_csv(OUT / "tables" / "h4b_concordance_summary.csv", index=False)
mef2_c = cc[cc.family == "MEF2"].iloc[0]
others = cc[cc.family != "MEF2"]
record("H4b", "The AP-1/MEF2 co-variation is SPECIFIC to MEF2, rather than every factor in the peak moving together.",
       "Concordance of the three-condition pattern between AP-1 and each neighbouring family, MEF2 against the others", int(cc.n.sum()),
       f"MEF2 agrees with AP-1 at {mef2_c.concordance:.1%} of loci; other families {others.concordance.min():.1%} to {others.concordance.max():.1%} "
       f"(KLF/SP {cc[cc.family == 'KLF / SP'].concordance.iloc[0]:.1%}, other {cc[cc.family == 'other'].concordance.iloc[0]:.1%}). "
       f"When AP-1 is condition-variable, the neighbouring factor is too in {mef2_c.partner_variable_when_ap1_variable:.0%} of MEF2 cases and {others.partner_variable_when_ap1_variable.mean():.0%} on average for the others.",
       "NOT SUPPORTED" if mef2_c.concordance <= others.concordance.max() else "SUPPORTED",
       "If every family agrees with AP-1 equally, the co-variation is peak-level signal, not co-regulation.")

# ================================================================ H8. does the Age Up motif sensitivity survive matching on MEASURED accessibility?
rows = []
for lab, cols in (("predicted context only (4 variables)", ["l_170", "l_peak", "l_pw", "l_tss"]), ("+ measured accessibility (5 variables)", MCOLS)):
    P8 = make_pairs(d7, "age_score_class", "Age Up (score)", "Unchanged", cols)
    for v in TYPES:
        a, b, p_, n_, md = paired_num(d7, P8, v)
        rows.append({"matching": lab, "mutation": v.split(",")[0], "n_pairs": n_, "mean_diff": md, "p": p_})
H8 = pd.DataFrame(rows)
H8.to_csv(OUT / "tables" / "h8_measured_accessibility_matching.csv", index=False)
n8 = {k: int((w.p < 0.05).sum()) for k, w in H8.groupby("matching")}
record("H8", "The extra motif sensitivity of Age Up (delta score) sites is robust to matching on MEASURED accessibility, not only predicted.",
       "Age Up vs Unchanged on four motif-level edits, matched on 4 predicted-context variables vs the same plus measured accessibility", int(H8.n_pairs.max()),
       f"significant motif edits: {n8.get('predicted context only (4 variables)', 0)} of {len(TYPES)} with predicted context only, "
       f"{n8.get('+ measured accessibility (5 variables)', 0)} of {len(TYPES)} once measured accessibility is matched too.",
       "SUPPORTED" if n8.get("+ measured accessibility (5 variables)", 0) > 0 else "NOT SUPPORTED",
       "Measured accessibility is partly a consequence of the binding change being studied, so matching on it may be over-adjustment; this records the sensitivity either way.")

# ================================================================ H9. are the non-responders simply the sites without a perfect core?
d9 = S.join(FD[["replace motif, read over motif"]])
d9 = d9.dropna(subset=["replace motif, read over motif"])
d9["non_responder"] = d9["replace motif, read over motif"] > np.log2(0.9)
tab9 = pd.crosstab(d9.exact_core, d9.non_responder)
or9, p9 = fisher_exact(tab9.values)
share_imperfect = d9[d9.exact_core == 0].non_responder.mean()
share_exact = d9[d9.exact_core == 1].non_responder.mean()
record("H9", "The AP-1 sites that do not respond to destroying their motif are the ones without a perfect core.",
       "Non-responder (replacing the motif changes accessibility by less than 10%) vs exact core, all sites", len(d9),
       f"{share_imperfect:.0%} of sites with an imperfect core are non-responders vs {share_exact:.0%} with a perfect core (odds ratio {1 / or9:.1f}, p = {p9:.2g}); "
       f"imperfect cores are {d9.exact_core.eq(0).mean():.0%} of sites and {d9[d9.non_responder].exact_core.eq(0).mean():.0%} of all non-responders.",
       "SUPPORTED" if p9 < 0.05 and share_imperfect > share_exact else "NOT SUPPORTED")

for e in REG:                                      # H4b is the control for H4: if the co-variation is not specific, H4 does not stand
    if e["id"] == "H4" and any(x["id"] == "H4b" and x["verdict"] == "NOT SUPPORTED" for x in REG):
        e["verdict"] = "RETRACTED"
        e["detail"] += " RETRACTED by H4b: every family in the peak co-varies with AP-1 about equally, so this is peak-level signal, not co-regulation."
R = pd.DataFrame(REG)
R.to_csv(OUT / "tables" / "hypothesis_register.csv", index=False)

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
ax = axes[0]
groups = [("Age Up\n(\u0394score)", ca, cb, p1), ("Age Down\n(\u0394score)", da, db, p1c), ("Dynamic\n(threshold)", ea, eb, p2)]
for k, (lab, a_, b_, p_) in enumerate(groups):
    ax.bar(k - 0.19, a_, width=0.36, color=ORANGE, alpha=0.9)
    ax.bar(k + 0.19, b_, width=0.36, color=BLUE, alpha=0.9)
    ax.text(k, max(a_, b_) + 0.02, ("p = %.1g" % p_) if p_ >= 0.001 else "p < 0.001", ha="center", fontsize=8)
ax.set_xticks(range(3), [g[0] for g in groups], fontsize=8.5)
ax.set_ylim(0, 1)
ax.set_ylabel("Share with an exact TGA[C/G]TCA core")
ax.set_title("H1, H1b, H2: core quality by class definition", fontsize=10, loc="left")
ax.bar(0, 0, color=ORANGE, label="class")
ax.bar(0, 0, color=BLUE, label="matched control")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax = axes[1]
v = d6.dropna(subset=["replace motif, read over motif"])
for k, (mmv, c) in enumerate(((0, AQUA), (1, ORANGE), (2, RED))):
    w = v[v.core_mismatches == mmv]["replace motif, read over motif"]
    if len(w) < 5:
        continue
    bp = ax.boxplot(w, positions=[k], widths=0.55, showfliers=False, patch_artist=True,
                    boxprops={"facecolor": c, "alpha": 0.55}, medianprops={"color": INK})
    ax.text(k, w.median() + 0.08, f"n = {len(w):,}", ha="center", fontsize=7.5)
ax.set_xticks(range(3), ["perfect core", "1 mismatch", "2 mismatches"], fontsize=8.5)
ax.set_ylabel("Effect of replacing the motif (log2FC)")
ax.set_title("H6: does the model see core quality?", fontsize=10, loc="left")
ax = axes[2]
w = H7.pivot_table(index="mutation", columns="matching", values="mean_diff")
ys = np.arange(len(w))
ax.scatter(w["context only"], ys, s=46, color=ORANGE, label="matched on context", zorder=3)
ax.scatter(w["context + core quality"], ys, s=46, color=BLUE, label="+ core quality", zorder=3)
for y_, a_, b_ in zip(ys, w["context only"], w["context + core quality"]):
    ax.plot([a_, b_], [y_, y_], color=GREY, lw=1.1, zorder=2)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(ys, [i.split(",")[0] for i in w.index], fontsize=8)
ax.set_xlabel("Age Up (\u0394score) minus Unchanged (log2FC)")
ax.set_title("H7: core quality as the mechanism (\u0394score classes)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "v1_hypotheses.png", bbox_inches="tight")
plt.close(fig)

# second figure: the register itself, and the two controls that changed a verdict
fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2), gridspec_kw={"width_ratios": [1.5, 1, 1]})
ax = axes[0]
COL = {"SUPPORTED": AQUA, "NOT SUPPORTED": GREY, "CONTRADICTED": RED, "RETRACTED": RED}
ys = np.arange(len(R))[::-1]
for y_, r_ in zip(ys, R.itertuples()):
    ax.barh(y_, 1, color=COL.get(r_.verdict, GREY), alpha=0.85, height=0.7)
    ax.text(0.03, y_, f"{r_.id}  {r_.verdict}", va="center", fontsize=8.5, color="white" if r_.verdict != "NOT SUPPORTED" else INK, fontweight="bold")
ax.set_yticks(ys, [h[:62] + ("..." if len(h) > 62 else "") for h in R.hypothesis], fontsize=7.5)
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.grid(False)
ax.set_title("Hypothesis register", fontsize=10, loc="left")
ax = axes[1]
cc = pd.read_csv(OUT / "tables" / "h4b_concordance_summary.csv").sort_values("concordance")
ax.barh(np.arange(len(cc)), cc.concordance, color=[ORANGE if f == "MEF2" else GREY for f in cc.family], alpha=0.9, height=0.62)
ax.set_yticks(np.arange(len(cc)), cc.family, fontsize=8.5)
ax.set_xlim(0.5, 1)
ax.set_xlabel("Share of loci where the partner's\nYS/MAS/MAE pattern equals AP-1's")
ax.set_title("H4b: is the co-variation MEF2-specific?", fontsize=10, loc="left")
ax = axes[2]
h8 = pd.read_csv(OUT / "tables" / "h8_measured_accessibility_matching.csv")
w = h8.pivot_table(index="mutation", columns="matching", values="mean_diff")
ys = np.arange(len(w))
c4 = [c for c in w.columns if c.startswith("predicted")][0]
c5 = [c for c in w.columns if c.startswith("+")][0]
ax.scatter(w[c4], ys, s=48, color=ORANGE, label="predicted context only", zorder=3)
ax.scatter(w[c5], ys, s=48, color=BLUE, label="+ measured accessibility", zorder=3)
for y_, a_, b_ in zip(ys, w[c4], w[c5]):
    ax.plot([a_, b_], [y_, y_], color=GREY, lw=1.1, zorder=2)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_yticks(ys, [i for i in w.index], fontsize=8)
ax.set_xlabel("Age Up (\u0394score) minus Unchanged (log2FC)")
ax.set_title("H8: robustness of the Age Up (\u0394score) effect", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_ylim(-0.6, len(w) - 0.2)
fig.tight_layout()
fig.savefig(OUT / "figures" / "v2_register_and_controls.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 120)
print("\n================ REGISTER ================")
print(R[["id", "verdict", "hypothesis"]].to_string(index=False))
print("\n== H7 detail\n", H7.round(3).to_string())
