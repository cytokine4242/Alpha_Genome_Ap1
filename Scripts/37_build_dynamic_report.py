#!/usr/bin/env python3
"""Standalone HTML document for one question: in the AlphaGenome mutation results, do Dynamic AP-1 sites respond differently from Stable ones?
Reads only saved tables and figures (Analysis/dynamic, Analysis/dynamic_profiles, Analysis/scan_titrate); images are embedded, so the file is self-contained.
  python3 37_build_dynamic_report.py          ->  Analysis/report_dynamic_vs_stable.html"""
from __future__ import annotations

import base64
import html
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "Analysis"
DY, DPR, ST = AN / "dynamic", AN / "dynamic_profiles", AN / "scan_titrate"
SVC, MH, SF, TI, MS = AN / "stable_vs_changing", AN / "motif_heatmap", AN / "scan_features", AN / "titration", AN / "mpra_scan"
PT, SCL, DR, CX, DEP, HYP = AN / "partners", AN / "score_classes", AN / "distance_rules", AN / "class_context", AN / "dependence", AN / "hypotheses"

IMGS: dict[int, Path] = {}


def img(path, alt):
    k = len(IMGS)
    IMGS[k] = Path(path)
    return f'<img src="@@IMG{k}@@" alt="{html.escape(alt)}" style="max-width:100%">'


def fig(path, alt, caption):
    return f"<figure>{img(path, alt)}<figcaption>{caption}</figcaption></figure>"


def table(df, fmt=None):
    d = df.copy()
    for c in d.columns:
        if fmt and c in fmt:
            d[c] = d[c].map(lambda v, f=fmt[c]: "" if pd.isna(v) else f.format(v))
        elif d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.2f}")
    return d.to_html(index=False, escape=False, classes="t", border=0)


def U(*items):
    return "<ul class='pts'>" + "".join(f"<li>{x}</li>" for x in items) + "</ul>"


def badge(kind):
    lab = {"ok": "&#10003; High confidence", "mid": "~ Moderate confidence", "inc": "? Inconclusive",
           "yes": "&#10003; Supported", "no": "&ndash; Not supported", "bad": "&#10007; Contradicted", "ret": "&#10007; Retracted"}[kind]
    return f'<span class="badge {kind}">{lab}</span>'


def bio(text):
    return f'<div class="bio"><strong>What it means</strong> &middot; {text}</div>'


# ---------------------------------------------------------------- tables
dr = lambda n: pd.read_csv(DY / "tables" / n)
pr = lambda n: pd.read_csv(DPR / "tables" / n)
PROF = dr("mutation_type_profile_dynamic_vs_stable.csv").set_index("mutation")
MATCH, BAL = dr("matched_pairs.csv").set_index("effect"), dr("matching_balance.csv")
SCAN, MAUC, AUC = dr("feature_scan_matched_pairs.csv"), dr("classifier_auc_matched_sample.csv").set_index("feature_set"), dr("classifier_auc.csv").set_index("features")
PAT, FBR = dr("mutation_type_by_binding_pattern.csv"), dr("fibroblast_replication_by_pattern.csv")
COND, VAR, SAME, CROSS = dr("condition_specific_accessibility.csv").set_index("binding_pattern"), dr("accessibility_variability_across_conditions.csv").iloc[0], dr("same_peak_pairs.csv").set_index("quantity"), dr("cross_cell_type_concordance.csv").iloc[0]
A, CLF, B, N = pr("a_response_shape_by_mutation_type.csv"), pr("a_classifier_matched_pairs.csv").set_index("features"), pr("b_ratios_matched_pairs.csv"), pr("c_non_responders.csv").iloc[0]
DS_, DT_, DTH = pr("d_tile_scan_summary.csv"), pr("d_tile_scan_tests.csv"), pr("d_tile_scan_thresholds.csv")

NP = int(MATCH.loc["170 bp shuffle", "n_pairs"])
ND = int(MATCH.loc["170 bp shuffle", "n_dynamic_total"])
sig = A[A.fdr < 0.05]
clf = lambda k: CLF.loc[[i for i in CLF.index if i.startswith(k)][0]]
ds = lambda rd, g: DS_[(DS_.readout == rd) & (DS_.group == g)].iloc[0]
dt = lambda rd, m: DT_[(DT_.readout == rd) & (DT_.comparison == "dynamic vs matched_stable") & (DT_.measure == m)].iloc[0]
DYN, MST, HST = ds("170 bp", "Dynamic"), ds("170 bp", "Stable, context-matched"), ds("170 bp", "Stable, high accessibility")
mm = PAT[(PAT.pattern == "MAS+MAE") & (PAT.fdr < 0.05)]
fb = lambda mut: FBR[(FBR.track == "FB sites, Patski") & (FBR.mutation == mut) & (FBR.pattern == "MAS+MAE")].iloc[0]

ptr = lambda n: pd.read_csv(PT / "tables" / n)
scr = lambda n: pd.read_csv(SCL / "tables" / n)
FAMI = ptr("a_tile_importance_by_family.csv").set_index("family")
FAMB = ptr("b_base_importance_by_family.csv").set_index("family")
FAMC = ptr("c_dynamic_vs_stable_by_family.csv").set_index("family")
PI = ptr("i_class_differences_site_level.csv")
PCT, PPRV = ptr("i_ctcf_site_level.csv").set_index("group"), ptr("i_prevalence_stratified.csv").set_index("family")
pi = lambda cmp_, fam_: PI[(PI.comparison == cmp_) & (PI.family == fam_)].iloc[0]
PCNT, PDET, PPOS = ptr("c_partner_counts_by_group.csv").set_index("group"), ptr("e_partner_count_vs_detection.csv"), ptr("f_partner_positions.csv").set_index("family")
FREE = float((FAMI.median_adj_effect - FAMI.diff_vs_footprint_free).median())
SAG, SVM, SBM = scr("a_agreement.csv").set_index("contrast"), scr("b_vs_measured_change.csv").set_index("contrast"), scr("b_measured_change_by_class.csv")
SHS = scr("h_score_profile_by_class.csv")
hs = lambda cls, col: SHS[SHS["class"] == cls].iloc[0][col]
SMU, SEM, SFM, SCO = scr("c_mutation_response_by_score_class.csv"), scr("e_matching_with_score.csv"), scr("f_motif_specificity.csv"), scr("g_core_sequence_by_score_class.csv").set_index("class")
AGEUP = SMU[SMU.contrast == "Age Up (score)"].set_index("mutation")
AGEDN = SMU[SMU.contrast == "Age Down (score)"].set_index("mutation")
_ns = SEM.groupby(["contrast", "matching"]).fdr.apply(lambda v: int((v < 0.05).sum())).unstack()
bm = lambda defn, cls, col: SBM[(SBM.definition == defn) & (SBM["class"] == cls)].iloc[0][col]


drr = lambda n: pd.read_csv(DR / "tables" / n)
DNEAR, DPER, DCLOSE = drr("b_distance_summary.csv").set_index("family"), drr("c_periodicity.csv").set_index("family"), drr("f_close_range_peak.csv").set_index("family")
DMODE, DRHO, DIMP, DCLS = drr("f_modal_distance.csv").set_index("family"), drr("d_distance_importance_correlation.csv").set_index("family"), drr("d_importance_by_distance.csv"), drr("e_class_differences_in_spacing.csv")
di = lambda fam_, b_: DIMP[(DIMP.family == fam_) & (DIMP.distance_bin == b_)].median_importance_vs_free.iloc[0]
dc = lambda fam_, cmp_: DCLS[(DCLS.family == fam_) & (DCLS.comparison == cmp_)].iloc[0]


CXA = pd.read_csv(CX / "tables" / "a_context_by_class.csv")
CXB = pd.read_csv(CX / "tables" / "b_mef2_presence_by_class.csv").set_index("comparison")
SEQ = pd.read_csv(DEP / "tables" / "a_sequence_motif_counts_by_class.csv")
SHP = pd.read_csv(DEP / "tables" / "b_dependence_shape_by_class.csv")
cx = lambda c, m: CXA[(CXA.comparison == c) & (CXA.measure == m)].iloc[0]
sq = lambda c, m: SEQ[(SEQ.comparison == c) & (SEQ.measure == m)].iloc[0]
sh = lambda c, m: SHP[(SHP.comparison == c) & (SHP.measure == m)].iloc[0]


HREG = pd.read_csv(HYP / "tables" / "hypothesis_register.csv").set_index("id")
HCC = pd.read_csv(HYP / "tables" / "h4b_concordance_summary.csv").set_index("family")
H8T = pd.read_csv(HYP / "tables" / "h8_measured_accessibility_matching.csv")
VCOL = {"SUPPORTED": "yes", "NOT SUPPORTED": "no", "CONTRADICTED": "bad", "RETRACTED": "ret"}


def register_table():
    rows = []
    for i, r in HREG.iterrows():
        rows.append(f"<tr><td><b>{i}</b></td><td>{badge(VCOL.get(r.verdict, 'inc'))}</td><td>{html.escape(r.hypothesis)}</td>"
                    f"<td>n = {int(r.n):,}</td><td>{html.escape(str(r.result))}</td></tr>")
    return ("<table class='t reg'><tr><th>ID</th><th>Verdict</th><th>Hypothesis</th><th>n</th><th>Result</th></tr>" + "".join(rows) + "</table>")


SVCR = pd.read_csv(SVC / "tables" / "stable_vs_changing.csv")
MSA = pd.read_csv(MS / "tables" / "motif_share_both_assays.csv").set_index("dataset")
MSB = pd.read_csv(MS / "tables" / "class_direction_both_assays.csv")
msb = lambda a, m: MSB[(MSB.assay.str.startswith(a)) & (MSB.measure == m)].iloc[0]
AGD, MPD = "AlphaGenome (CM, genomic context)", "MPRA (fibroblast, episomal)"
SFR = pd.read_csv(SF / "tables" / "scan_features_by_class.csv")
TIR = pd.read_csv(TI / "tables" / "titration_by_class.csv")
TIO = pd.read_csv(TI / "tables" / "titration_overall.csv").set_index("feature")
TIV = pd.read_csv(TI / "tables" / "validity_checks.csv")
sf = lambda c, f: SFR[(SFR.contrast == c) & (SFR.feature.str.startswith(f))].iloc[0]
ti = lambda c, f: TIR[(TIR.contrast == c) & (TIR.feature.str.startswith(f))].iloc[0]
nsf = lambda c: int((SFR[SFR.contrast == c].detects == "YES").sum())
nti = lambda c: int((TIR[TIR.contrast == c].detects == "YES").sum())
GINI = "concentration of dependence (Gini)"
C_AU_D, C_AD_D = "Age Up vs Unchanged (\u0394score)", "Age Down vs Unchanged (\u0394score)"
MHT = pd.read_csv(MH / "tables" / "top5_each_direction.csv")
MHA = pd.read_csv(MH / "tables" / "importance_by_class.csv").set_index("family")


def top5(contrast, direction, n=5):
    w = MHT[(MHT.contrast == contrast) & (MHT.direction == direction)].head(n)
    return ", ".join(f"<b>{r.family.split(' (')[0]}</b> ({r.difference:+.2f}, n = {int(r.n_class)}{', p = %.3f' % r.p if r.p < 0.05 else ''})" for r in w.itertuples())


AGEUP_T, AGEDN_T = "Age Up vs Stable (threshold)", "Age Down vs Stable (threshold)"
EXUP_T, EXDN_T = "Ex Up vs Stable (threshold)", "Ex Down vs Stable (threshold)"
MORE, LESS = "more important in this class", "more important in the reference"
sv = lambda c, m: SVCR[(SVCR.contrast == c) & (SVCR.measure == m)].iloc[0]
TH, CH, AU, AD = "Threshold: Dynamic vs Stable", "\u0394score: Changing vs Steady", "\u0394score: Age Up vs Unchanged", "\u0394score: Age Down vs Unchanged"
nyes = lambda c: int((SVCR[SVCR.contrast == c].detects == "YES").sum())


# ---------------------------------------------------------------- sections
SEC = []


def sec(key, title, kind, why, points, extra=""):
    SEC.append((key, title, f'<h2 id="{key}"><span class="num">{len(SEC) + 1}</span>{title}</h2>'
                            f'<p class="conf">{badge(kind)} {why}</p>{points}{extra}'))


sec("design", "What was compared", "ok",
    "Definitions come from the TOBIAS binding calls; the matching uses AlphaGenome's own predictions and peak geometry only.",
    U("<b>Stable</b> = the AP-1 footprint is called in all three conditions (young sedentary, aged sedentary, aged + exercise). <b>Dynamic</b> = anything else. 3,679 cardiomyocyte Fos sites: 2,907 Stable, 772 Dynamic.",
      "<b>Effect</b> = log2(edited / native) predicted accessibility, averaged over replicate edits. Negative = the edit closes the site. &minus;1 means half the accessibility is lost.",
      "<b>Readout</b>: each edit is read over the region named with it (the 11 bp motif, the 170 bp around it, or the whole peak). Reading an edit over a region it did not touch measures spread, not dependence.",
      f"<b>Context matching</b>: each Dynamic site is paired with the Stable site closest in <i>predicted</i> accessibility (170&nbsp;bp and peak), peak width, distance to the nearest TSS and <i>measured</i> snATAC accessibility. {NP} of {ND} Dynamic sites find a partner within the caliper; after matching the five covariates differ by at most {BAL.std_diff_after.abs().max():.2f} SD.",
      "<b>Correction to an earlier version of this report.</b> The matching originally used the four predicted/geometric covariates only. That leaves measured accessibility unbalanced by 0.81&nbsp;SD &mdash; Dynamic sites really do sit in weaker chromatin even when the model predicts the same opening &mdash; and the unbalanced covariate was suppressing the very difference the matching was meant to isolate. With measured accessibility matched as well, the 170&nbsp;bp difference no longer vanishes. Everything below uses the five-covariate matching.",
      "<b>Two different definitions of 'responsive' are used below, and they pick out different sites. Keep them apart:</b>",
      "<b>Age Up (threshold)</b> &mdash; a <i>discrete</i> change in the binding call: no AP-1 footprint in young sedentary, a footprint in aged sedentary. 396 cardiomyocyte sites. This is the original TOBIAS class and the one the Stable/Dynamic split is built from.",
      "<b>Age Up (score)</b> &mdash; a <i>graded</i> change: the top fifth of the change in TOBIAS binding score (MAS &minus; YS), whether or not the call ever crossed the threshold. 736 sites. Only 22% of these ever cross the threshold, and 40% of threshold-Stable sites fall in a &Delta;score class.",
      "Where a result belongs to only one of them it is named in full. The core-quality and motif-sensitivity findings are <b>&Delta;score</b> results; the 2.2:1 excess of gains over losses with age and the Stable/Dynamic contrast are <b>threshold</b> results."),
    bio("Dynamic sites sit in weaker chromatin than Stable ones. Any comparison that does not hold that constant measures the chromatin, not the sequence. Every result below is reported raw and context-matched."))

sec("raw", "Most of the raw difference is chromatin &mdash; but not all of it", "ok",
    "Large raw differences, most of which matching removes; what survives is confined to one scale.",
    U(f"Raw, Dynamic sites lose less to every wide edit: replace the 170 bp {PROF.loc['replace the 170 bp', 'stable_median']:.2f} (Stable) vs {PROF.loc['replace the 170 bp', 'dynamic_median']:.2f} (Dynamic); shuffle the 170 bp {PROF.loc['shuffle the 170 bp', 'stable_median']:.2f} vs {PROF.loc['shuffle the 170 bp', 'dynamic_median']:.2f}; centred shuffle of 340 bp {PROF.loc['centred shuffle, 340 bp', 'stable_median']:.2f} vs {PROF.loc['centred shuffle, 340 bp', 'dynamic_median']:.2f} (all FDR &lt; 0.001).",
      f"Motif-level edits differ by only about 0.1: replace the motif {PROF.loc['replace the motif', 'stable_median']:.2f} vs {PROF.loc['replace the motif', 'dynamic_median']:.2f}; replace the 7 bp core {PROF.loc['replace the 7 bp core', 'stable_median']:.2f} vs {PROF.loc['replace the 7 bp core', 'dynamic_median']:.2f}.",
      "Edits far from the motif (&plusmn;300 and &plusmn;600 bp away) do not differ at all, which is the control: the two classes respond identically where nothing is at stake.",
      f"<b>Matching removes most of it.</b> The 170 bp gap falls from {MATCH.loc['170 bp shuffle', 'raw_gap']:+.2f} raw to {MATCH.loc['170 bp shuffle', 'matched_gap_median']:+.2f} (CI {MATCH.loc['170 bp shuffle', 'matched_ci_lo']:+.2f} to {MATCH.loc['170 bp shuffle', 'matched_ci_hi']:+.2f}) &mdash; about {1 - MATCH.loc['170 bp shuffle', 'matched_gap_median'] / MATCH.loc['170 bp shuffle', 'raw_gap']:.0%} of the raw difference was chromatin.",
      f"<b>What survives is confined to one scale.</b> Edits of the 11 bp motif no longer differ (replace the motif {PROF.loc['replace the motif', 'matched_stable_median']:.2f} vs {PROF.loc['replace the motif', 'matched_dynamic_median']:.2f}, FDR {PROF.loc['replace the motif', 'fdr_matched']:.2f}; the 7 bp core FDR {PROF.loc['replace the 7 bp core', 'fdr_matched']:.2f}), and neither does any edit of the whole peak (shuffle {PROF.loc['shuffle the peak', 'fdr_matched']:.2f}, replace {PROF.loc['replace the peak', 'fdr_matched']:.2f}, dinucleotide {PROF.loc['dinucleotide-shuffle the peak', 'fdr_matched']:.2f}). But every edit of the 170 bp does: replace {PROF.loc['replace the 170 bp', 'matched_stable_median']:.2f} vs {PROF.loc['replace the 170 bp', 'matched_dynamic_median']:.2f}, shuffle {PROF.loc['shuffle the 170 bp', 'matched_stable_median']:.2f} vs {PROF.loc['shuffle the 170 bp', 'matched_dynamic_median']:.2f}, dinucleotide-shuffle {PROF.loc['dinucleotide-shuffle the 170 bp', 'matched_stable_median']:.2f} vs {PROF.loc['dinucleotide-shuffle the 170 bp', 'matched_dynamic_median']:.2f} (all FDR &lt; 0.001).",
      f"<b>The width series says the same thing.</b> Walking the centred shuffle out from 7 bp to 680 bp, the matched gap is near zero at motif width ({PROF.loc['centred shuffle, 11 bp', 'matched_dynamic_median'] - PROF.loc['centred shuffle, 11 bp', 'matched_stable_median']:+.2f}), peaks around the element ({PROF.loc['centred shuffle, 85 bp', 'matched_dynamic_median'] - PROF.loc['centred shuffle, 85 bp', 'matched_stable_median']:+.2f} at 85 bp), and closes again once the edit is wide enough to cover the whole neighbourhood ({PROF.loc['centred shuffle, 680 bp', 'matched_dynamic_median'] - PROF.loc['centred shuffle, 680 bp', 'matched_stable_median']:+.2f} at 680 bp).",
      f"The far controls still hold: &plusmn;300 bp FDR {PROF.loc['170 bp shuffle ±300 bp away', 'fdr_matched']:.2f}, &plusmn;600 bp FDR {PROF.loc['170 bp shuffle ±600 bp away', 'fdr_matched']:.2f}."),
    fig(DY / "figures" / "y5_mutation_type_profile.png", "Sensitivity to each mutation type",
        "Median effect of each mutation type for Stable (blue) and Dynamic (orange) sites with 95% bootstrap CI. Left marks: all sites (Mann&ndash;Whitney, BH). Right marks: context-matched pairs (paired Wilcoxon, BH).")
    + bio("Most of the raw gap is a chromatin effect: a site in weaker chromatin has less accessibility to lose, so every broad edit takes less from it. What is left after matching is scale-specific, and it is informative. The motif is worth the same at both classes, and so is the whole peak; only the 170 bp element around the motif is worth less at a Dynamic site. "
          "That is what you would expect if a Dynamic site's dependence is spread more thinly across its peak &mdash; a 170 bp edit catches a smaller share of it, while a peak-wide edit catches all of it either way. The dependence-shape section below measures that spread directly, and it agrees."))

sec("shape", "Same shape: the mix of mutation types at equal sensitivity", "mid",
    "Matched pairs, paired tests and a permutation-checked classifier. A null result, so a small difference in shape could be missed.",
    U("The matched comparison above could still hide a difference in <i>proportions</i>: perhaps Dynamic sites lose less overall but depend relatively more on their motif. To test that, each site's 21 effects were divided by their overall size, so every site has the same total sensitivity and only the mix is compared.",
      f"The mix is the same. {len(sig)} of {len(A)} mutation types differ at FDR &lt; 0.05, all four of them edits of the 170 bp, and all by {sig.paired_median_diff.min():+.3f} to {sig.paired_median_diff.max():+.3f} on shares of about {sig.share_stable.abs().mean():.2f} &mdash; around 5% of each share.",
      f"A classifier given the whole profile cannot tell matched Dynamic from Stable sites apart: shape alone AUC {clf('shape only').auc:.2f}, raw effects {clf('size and shape').auc:.2f}, overall size {clf('overall size').auc:.2f}. Shuffling the labels gives about {CLF.null_auc_95th.mean():.2f}, so these are at chance (permutation p &ge; {CLF.p_perm.min():.2f}).",
      "Five ratios with a plain meaning were compared at sites that actually respond. None differs: " + "; ".join(f"<b>{r.ratio.split(':')[0]}</b> {r.median_stable:.2f} vs {r.median_dynamic:.2f}" for r in B.itertuples()) + ".",
      f"Nor are Dynamic sites a different mixture of AP-1-dependent and AP-1-independent sites: {N.share_dynamic_matched:.0%} of Dynamic and {N.share_stable_matched:.0%} of matched Stable sites are non-responders, where replacing the motif changes accessibility by less than 10% (McNemar p = {N.p_mcnemar_matched:.2f}; {N.share_dynamic_all:.0%} vs {N.share_stable_all:.0%} in all sites)."),
    fig(DPR / "figures" / "z1_response_shape.png", "Response shape by mutation type",
        "Paired difference (Dynamic minus matched Stable) in the share of each site's response taken by each mutation type, median with 95% bootstrap CI. Orange = FDR &lt; 0.05 (paired Wilcoxon, BH). Each site's 21 effects are scaled to unit length first.")
    + f"<details><summary>Classifier and the five ratios</summary>"
      + fig(DPR / "figures" / "z2_classifier.png", "Classifier on response profiles",
            "Cross-validated AUC in the matched sample, folds grouped by peak, averaged over 5 random fold assignments. Dotted line = 95th percentile of 100 label permutations.")
      + table(B[["ratio", "n_pairs_both_responding", "median_stable", "median_dynamic", "paired_median_diff", "p_wilcoxon", "fdr"]].rename(
          columns={"ratio": "Ratio", "n_pairs_both_responding": "n pairs", "median_stable": "Stable", "median_dynamic": "Dynamic", "paired_median_diff": "Paired diff", "p_wilcoxon": "p", "fdr": "FDR"}))
      + "</details>"
    + bio("Responsive and constitutive AP-1 sites are built the same way. They depend on the same parts of their sequence in the same proportions, and differ only in how much accessibility is at stake. This is the clearest answer the perturbation data give to your question."))

sec("patterns", "Do the individual binding patterns differ?", "inc",
    "126 comparisons; the one pattern that stood out does not replicate in fibroblast sites.",
    U("Stable/Dynamic is a coarse split, so each of the six Dynamic binding patterns was also compared with Stable separately, adjusted for accessibility at the readout, peak accessibility, peak width and TSS distance (splines).",
      f"{int((PAT.fdr < 0.05).sum())} of {int((PAT.pattern != 'All Dynamic').sum())} pattern-by-mutation comparisons are significant, {len(mm)} of them for a single pattern: MAS+MAE (bound in both aged states, not in young), which is less sensitive to the 170 bp and wider edits by {mm[mm.adj_diff > 0].adj_diff.min():+.2f} to {mm[mm.adj_diff > 0].adj_diff.max():+.2f} log2.",
      f"It does not replicate. In fibroblast sites with the Patski track, MAS+MAE shows no 170 bp difference ({fb('shuffle 170→170').adj_diff:+.2f}) and is slightly <i>more</i> sensitive to motif and core edits ({fb('shuffle motif→motif').adj_diff:+.2f}, {fb('replace core→motif').adj_diff:+.2f}).",
      "The remaining significant cells are isolated and inconsistent between tracks."),
    f"<details><summary>Full pattern-by-mutation figure</summary>"
    + fig(DY / "figures" / "y6_mutation_type_by_pattern.png", "Sensitivity by binding pattern and mutation type",
          "Difference from Stable for each binding pattern (columns) and mutation type (rows): raw (left) and adjusted (right). Marks are BH-corrected over all pattern-by-mutation cells. Positive = the pattern loses less.")
    + "</details>"
    + bio("No binding pattern carries a specific vulnerability or resistance. MAS+MAE is a hypothesis worth revisiting with more sites or a fine-tuned head, not a finding."))

sec("else", "Anything else in the model's outputs?", "mid",
    "Cross-validated; the matched sample has 1,440 sites, so small differences could be missed.",
    U(f"Every AlphaGenome output was scanned, not only the mutation effects: native predictions, the effect of every edit at every readout, RNA predictions and derived contrasts ({len(SCAN)} in all). {int((SCAN.fdr < 0.05).sum())} differ at FDR &lt; 0.05 between matched pairs, and the largest is {SCAN.matched_std_diff.abs().max():.2f} SD of the Stable sites.",
      f"A classifier given any subset is at chance in the matched sample: all outputs AUC {MAUC.loc['everything AlphaGenome outputs', 'auc_logistic']:.2f}, shuffle effects {MAUC.loc['shuffle effects (9)', 'auc_logistic']:.2f}, mutation effects {MAUC.loc['mutation effects (9)', 'auc_logistic']:.2f}, RNA predictions {MAUC.loc['RNA predictions', 'auc_logistic']:.2f}.",
      f"Across all sites, predicted accessibility alone separates the classes at AUC {AUC.loc['native predicted accessibility only (motif, 170 bp, peak)', 'auc_logistic']:.2f}, and adding every mutation effect changes nothing ({AUC.loc['context + model effects', 'auc_logistic']:.2f}). Only measured accessibility and the prediction error lift it, to {AUC.loc['context + effects + measured accessibility + prediction error', 'auc_logistic']:.2f}.",
      "Predicted RNA near the TSS is slightly lower at Dynamic sites (matched difference &minus;0.13 log10, FDR 0.053). Too weak to rely on."),
    f"<details><summary>Feature scan and classifiers</summary>"
    + fig(DY / "figures" / "y4_feature_scan.png", "Outputs that differ most between Dynamic and matched Stable sites",
          "The 22 largest standardised differences between Dynamic sites and their matched Stable partners (filled = FDR &lt; 0.05; grey diamonds = all sites, unmatched).")
    + fig(DY / "figures" / "y1_classifier.png", "Classifiers",
          "Cross-validated AUC for telling Dynamic from Stable with different feature sets, in all sites (left) and the matched sample (right). Folds grouped by peak.")
    + "</details>"
    + bio("The model's sequence-level outputs carry no signature of a Dynamic site. What separates the classes is the accessibility the model predicts, and the part of the measurement it fails to predict."))

sec("differs", "What does differ: context, and accessibility that changes with condition", "mid",
    "Small effects, but consistent across all six binding patterns; one pooled measurement per condition.",
    U(f"Measured accessibility at Dynamic sites is higher in the conditions where they are called bound: {COND.loc['All Dynamic', 'measured_log2_bound_minus_unbound']:+.2f} log2 on average, about 11%, in every pattern (p &lt; 0.001). It also varies more across conditions than at Stable sites (median SD {VAR.dynamic_sd:.2f} vs {VAR.stable_sd:.2f} log2).",
      f"The model predicts one baseline value per site and over-predicts Dynamic sites more in the conditions where they are unbound ({COND.loc['All Dynamic', 'overprediction_diff_unbound_minus_bound']:+.2f} rank units, every pattern FDR &lt; 0.001): the sequence looks permissive, and the site is open mostly when bound.",
      f"Within a single peak holding both classes ({int(SAME.loc['170 bp shuffle effect', 'n_peaks'])} peaks), the Dynamic site sits at lower accessibility ({SAME.loc['native accessibility, 170 bp (log10)', 'median_dynamic_minus_stable']:+.2f} log10) and loses less ({SAME.loc['170 bp shuffle effect', 'median_dynamic_minus_stable']:+.2f}): the effect tracks local accessibility even inside one peak.",
      f"Class is not a fixed property of the DNA. Of {int(CROSS.shared_sites)} sites present in both cell-type tables, {CROSS.agreement:.0%} keep their class, but Cohen's kappa is only {CROSS.cohens_kappa:.2f}, and only {int(CROSS.dynamic_both)} sites are Dynamic in both."),
    f"<details><summary>Condition-specific accessibility and paired comparisons</summary>"
    + fig(DY / "figures" / "y2_condition_specific.png", "Condition-specific accessibility at Dynamic sites",
          "Left: measured accessibility, bound minus unbound conditions, by binding pattern. Middle: the model's over-prediction in unbound minus bound conditions. Right: variability across the three conditions.")
    + fig(DY / "figures" / "y3_matched_and_same_peak.png", "Matched pairs and same-peak pairs",
          "The 170 bp gap in all sites, in context-matched pairs and in same-peak pairs (median, 95% bootstrap CI).")
    + "</details>"
    + bio("Dynamic sites are genuinely condition-specific in the measured data, and a baseline sequence model cannot see it: it predicts one number per site, with no notion of condition. That is the gap a fine-tuned, condition-aware head would be built to fill."))

sec("scoredef", "A second definition: change in binding score, not crossing the threshold", "mid",
    "A parallel definition, reported alongside the threshold classes rather than replacing them; one of the five contrasts survives every control.",
    U(f"The classes so far are threshold calls: a site is Age Up if its footprint is absent in young and called in aged. An alternative is the <b>size of the change in TOBIAS binding score</b>: age = MAS &minus; YS, exercise = MAE &minus; MAS, taking the top and bottom fifth of each.",
      f"<b>The two definitions disagree a great deal.</b> Only {SAG.loc['age', 'threshold_Up_also_score_Up']:.0%} of Age Up (threshold) sites are also score Age Up, only {SAG.loc['age', 'score_Up_that_crossed_threshold']:.0%} of score Age Up sites ever crossed the threshold, and {SAG.loc['age', 'threshold_Stable_called_changing']:.0%} of threshold-Stable sites show a large score change.",
      f"<b>The score-based definition tracks the measured data about twice as well.</b> Against the measured snATAC change between the same conditions, the score change gives Spearman {SVM.loc['age (MAS vs YS)', 'rho_score_change_vs_measured']:.2f} (age) and {SVM.loc['exercise (MAE vs MAS)', 'rho_score_change_vs_measured']:.2f} (exercise), against {SVM.loc['age (MAS vs YS)', 'rho_threshold_class_vs_measured']:.2f} and {SVM.loc['exercise (MAE vs MAS)', 'rho_threshold_class_vs_measured']:.2f} for the threshold classes. Both definitions point the right way: Age Up (threshold) sites gain {bm('threshold', 'Age Up', 'median_measured_log2'):+.2f} log2 of measured accessibility, score Age Up sites {bm('score', 'Age Up (score)', 'median_measured_log2'):+.2f}, on {int(bm('score', 'Age Up (score)', 'n'))} against {int(bm('threshold', 'Age Up', 'n'))} sites.",
      f"<b>And it finds something the threshold definition misses.</b> Matched on predicted accessibility, peak width and TSS distance, <b>Age Up (score) sites are more sensitive to edits of the AP-1 motif itself</b>: shuffle the motif {AGEUP.loc['shuffle the motif', 'matched_diff']:+.2f} log2 (FDR {AGEUP.loc['shuffle the motif', 'fdr_matched']:.3f}), dinucleotide-shuffle the motif {AGEUP.loc['dinucleotide-shuffle the motif', 'matched_diff']:+.2f} (FDR {AGEUP.loc['dinucleotide-shuffle the motif', 'fdr_matched']:.3f}), centred 7 bp {AGEUP.loc['centred shuffle, 7 bp', 'matched_diff']:+.2f}, centred 21 bp {AGEUP.loc['centred shuffle, 21 bp', 'matched_diff']:+.2f}. The wide edits are flat ({AGEUP.loc['shuffle the 170 bp', 'matched_diff']:+.2f} for shuffling the 170 bp), so the difference is concentrated on the motif.",
      f"<b>Age Down sites show the mirror image, and it does not survive.</b> Matched on context, sites LOSING binding score with age are more sensitive to the wide edits and not to the motif: shuffle the 170 bp {AGEDN.loc['shuffle the 170 bp', 'matched_diff']:+.2f} (FDR {AGEDN.loc['shuffle the 170 bp', 'fdr_matched']:.3f}), centred 340 bp {AGEDN.loc['centred shuffle, 340 bp', 'matched_diff']:+.2f}, centred 680 bp {AGEDN.loc['centred shuffle, 680 bp', 'matched_diff']:+.2f}, while shuffling the motif gives {AGEDN.loc['shuffle the motif', 'matched_diff']:+.2f} (ns) and replacing the core {AGEDN.loc['replace the 7 bp core', 'matched_diff']:+.2f} (ns).",
      f"That mirror image is regression to the mean: a site can only lose score if it had score to lose. Age Down sites start from the highest footprint scores ({hs('Age Down (score)', 'median_start_score'):.2f} in young against {hs('Age Up (score)', 'median_start_score'):.2f} for Age Up and {hs('Unchanged', 'median_start_score'):.2f} for Unchanged) and sit in the strongest chromatin (predicted accessibility {hs('Age Down (score)', 'median_native_170'):.2f} vs {hs('Unchanged', 'median_native_170'):.2f}). Sensitivity to wide edits follows chromatin strength, which is the same effect that produced the raw Stable-vs-Dynamic gap.",
      f"<b>It survives the obvious artefact.</b> Binding score tracks accessibility, so a score change could just mean a strong footprint. Adding the site's absolute footprint score to the matching leaves Age Up intact ({int(_ns.loc['Age Up (score)', 'context + absolute footprint score'])} significant mutation types, against {int(_ns.loc['Age Up (score)', 'context only'])} without it) while every other contrast collapses to zero: Age Down {int(_ns.loc['Age Down (score)', 'context only'])}&rarr;{int(_ns.loc['Age Down (score)', 'context + absolute footprint score'])}, Exercise Down {int(_ns.loc['Exercise Down (score)', 'context only'])}&rarr;{int(_ns.loc['Exercise Down (score)', 'context + absolute footprint score'])}, Exercise Up and Changing 0 either way. Gaining is therefore not the same as losing: only the gain carries a sequence signature.",
      f"<b>The sequence agrees, and only for the gain.</b> Among the scanned sites, {SCO.loc['Age Up (score)', 'share_with_exact_core']:.0%} of Age Up (score) sites carry an exact TGA[C/G]TCA core against {SCO.loc['Age Up (score)', 'share_control']:.0%} of Unchanged (p = {SCO.loc['Age Up (score)', 'p_fisher']:.3f}, n = {int(SCO.loc['Age Up (score)', 'n'])} vs {int(SCO.loc['Age Up (score)', 'n_control'])}). Age Down sites are no different from Unchanged ({SCO.loc['Age Down (score)', 'share_with_exact_core']:.0%} vs {SCO.loc['Age Down (score)', 'share_control']:.0%}, p = {SCO.loc['Age Down (score)', 'p_fisher']:.2f}), nor is either exercise class.",
      f"Caveat: the far-from-motif controls also shift slightly for Age Up ({SFM[(SFM.contrast == 'Age Up (score)') & (SFM.measure == 'far control, ±300 bp')].paired_diff.iloc[0]:+.3f} at &plusmn;300 bp), so a small part is general sensitivity, roughly a tenth of the motif effect. The formal motif-share ratio does not reach significance (p = {SFM[(SFM.contrast == 'Age Up (score)') & (SFM.measure == 'motif share of the local dependence')].p.iloc[0]:.2f})."),
    fig(SCL / "figures" / "q3_definitions_compared.png", "The two definitions compared",
        "Matched difference per mutation type under the threshold definition (grey squares, Dynamic vs Stable, none significant) and under the score definition (Age Up vs Unchanged; orange = FDR &lt; 0.05). Negative = the responsive class loses more.")
    + f"<details><summary>How the classes were built and how they compare with the measured change</summary>"
      + fig(SCL / "figures" / "q1_score_classes.png", "Score change and measured accessibility change",
            "Left: the distribution of score change for both contrasts, with the cut points. Right: the measured snATAC change for each class under both definitions, median with 95% bootstrap CI.")
      + fig(SCL / "figures" / "q2_mutation_response_score_classes.png", "Mutation responses under the score definition",
            "Changing (score) minus Steady (score) per mutation type, raw and matched. Defining Dynamic by score magnitude alone, with no direction, recovers the same null as the threshold definition.")
      + "</details>"
    + bio("Ageing does not just flip AP-1 sites on; it strengthens binding at sites that already carry the best AP-1 motifs, and those are the sites whose predicted accessibility depends most on the motif. Losing AP-1 binding with age looks nothing like the reverse: those sites are simply the strongest footprints in the strongest chromatin drifting back down, with no sequence signature of their own. Gain is motif-driven, loss is not &mdash; and the threshold definition sees neither, because it discards the size of the change."))

sec("partners", "Which factors matter: every footprint of the full TOBIAS table", "mid",
    f"{int(DTH.n_complete_sites.iloc[0])} scanned sites; the ranking is strong and holds under five different adjustments, but the comparisons between classes are underpowered.",
    U("Every footprinted cardiomyocyte TFBS in the full TOBIAS table (2.49 million instances, 786 motifs) was mapped onto the tile scan. For each 10 bp tile: which factors are bound there, and how much does replacing it cost? Effects are adjusted for the site's overall sensitivity and for distance from the AP-1 motif, and compared with tiles holding no footprint at all.",
      f"<b>MEF2 is the dominant partner, by a wide margin.</b> Replacing a tile over a MEF2 footprint costs {FAMI.loc['MEF2', 'median_adj_effect']:.2f} log2 (CI {FAMI.loc['MEF2', 'ci_lo']:.2f} to {FAMI.loc['MEF2', 'ci_hi']:.2f}; n = {int(FAMI.loc['MEF2', 'n_tiles'])} tiles over {int(FAMI.loc['MEF2', 'n_sites'])} sites) against about {FREE:.2f} for footprint-free tiles. Then GATA {FAMI.loc['GATA', 'median_adj_effect']:.2f}, TEAD {FAMI.loc['TEAD', 'median_adj_effect']:.2f}, CTCF {FAMI.loc['CTCF / cohesin', 'median_adj_effect']:.2f}.",
      f"<b>The abundant factors are not the important ones.</b> KLF/SP sites are the most common by far ({int(FAMI.loc['KLF / SP', 'n_tiles']):,} tiles) yet cost only {FAMI.loc['KLF / SP', 'median_adj_effect']:.2f}. Even other AP-1 instances nearby cost just {FAMI.loc['AP-1 (Fos/Jun/Atf/Batf)', 'median_adj_effect']:.2f} &mdash; a MEF2 neighbour matters roughly eight times more than a second AP-1 site.",
      f"The same ranking appears at base resolution: a single-base change inside a MEF2 footprint carries {FAMB.loc['MEF2', 'diff_vs_footprint_free']:+.2f} log2 more |effect| than a base with no footprint, ahead of CTCF {FAMB.loc['CTCF / cohesin', 'diff_vs_footprint_free']:+.2f}, TEAD {FAMB.loc['TEAD', 'diff_vs_footprint_free']:+.2f} and GATA {FAMB.loc['GATA', 'diff_vs_footprint_free']:+.2f}. The chr1:119,197,433 example was not a one-off.",
      f"Partners are not immediately adjacent: the median MEF2 footprint sits {PPOS.loc['MEF2', 'median_abs_distance']:.0f} bp from the AP-1 motif, and only {PPOS.loc['MEF2', 'share_within_50bp']:.0%} are within 50 bp.",
      f"<b>Do different factors matter to the two classes? No, with one lead.</b> Tiles within a site are not independent, so each site contributes one value per family, compared against its own footprint-free tiles. On that basis <b>no family differs in importance between Dynamic and context-matched Stable sites after correction</b> (smallest FDR {PI[PI.comparison == 'Dynamic vs context-matched Stable'].fdr.min():.2f}). An earlier tile-level test made CTCF look significant; that was pseudo-replication.",
      f"<b>The CTCF lead.</b> It is the largest difference in both comparisons and in the same direction: Dynamic {pi('Dynamic vs context-matched Stable', 'CTCF / cohesin').median_a:+.2f} against context-matched Stable {pi('Dynamic vs context-matched Stable', 'CTCF / cohesin').median_b:+.2f} (p = {pi('Dynamic vs context-matched Stable', 'CTCF / cohesin').p:.3f}, FDR {pi('Dynamic vs context-matched Stable', 'CTCF / cohesin').fdr:.2f}) and against high-accessibility Stable {pi('Dynamic vs high-accessibility Stable', 'CTCF / cohesin').median_b:+.2f} (FDR {pi('Dynamic vs high-accessibility Stable', 'CTCF / cohesin').fdr:.3f}). It is not driven by a few sites: {PCT.loc['Dynamic', 'share_sites_below_-0.3']:.0%} of the {int(PCT.loc['Dynamic', 'n_sites'])} Dynamic sites carrying a CTCF footprint show a shift beyond &minus;0.3, against {PCT.loc['Stable, context-matched', 'share_sites_below_-0.3']:.0%} and {PCT.loc['Stable, high accessibility', 'share_sites_below_-0.3']:.0%} of Stable sites. Worth a look, not yet a result.",
      f"<b>Which partners are present may differ, even though how many does not.</b> Stratifying by measured accessibility to remove the detection effect, MEF2 is about half as often present at Dynamic sites (odds ratio {PPRV.loc['MEF2', 'odds_ratio_stratified']:.2f}, CI {PPRV.loc['MEF2', 'or_ci_lo']:.2f} to {PPRV.loc['MEF2', 'or_ci_hi']:.2f}, Mantel&ndash;Haenszel p = {PPRV.loc['MEF2', 'p_mantel_haenszel']:.3f}, FDR {PPRV.loc['MEF2', 'fdr_stratified']:.3f} over 8 families), while CTCF runs the other way ({PPRV.loc['CTCF / cohesin', 'odds_ratio_stratified']:.2f}). Since the total number of partners is the same, this would be a change of repertoire rather than of amount &mdash; suggestive, and borderline after correction.",
      f"<b>Nor does partner importance differ by the score-based classes.</b> Age Up (score) against Unchanged scanned sites: no family survives correction (smallest FDR {PI[PI.comparison == 'Age Up (score) vs Unchanged'].fdr.min():.2f}); MEF2 is {pi('Age Up (score) vs Unchanged', 'MEF2').difference:+.2f} with only {int(pi('Age Up (score) vs Unchanged', 'MEF2').n_sites_a)} Age Up sites carrying one. With this many sites a difference below about {PI[PI.comparison == 'Age Up (score) vs Unchanged'].min_detectable_diff_80pct.median():.2f} log2 would be missed, so this is underpowered rather than negative.",
      f"<b>A tempting result that is an artefact:</b> Dynamic sites appear to have fewer partners ({PCNT.loc['Dynamic', 'median_families_per_site']:.0f} vs {PCNT.loc['Stable, context-matched', 'median_families_per_site']:.0f} families per site, p = 0.001). Footprint calling needs ATAC signal, and the number of families correlates with measured accessibility at rho {PDET.rho_families_vs_measured.iloc[0]:.2f}. Within each quartile of measured accessibility the difference disappears entirely (p = {PDET.p.iloc[1]:.2f}, {PDET.p.iloc[2]:.2f}, {PDET.p.iloc[3]:.2f}, {PDET.p.iloc[4]:.2f}). It is detection, not biology."),
    fig(PT / "figures" / "p1_tile_importance_by_family.png", "Tile importance by TF family",
        "Median effect of replacing a 10 bp tile overlapping each factor's footprint, with 95% bootstrap CI, adjusted for the site's overall sensitivity and for distance from the AP-1 motif. Dashed line: tiles with no footprint. Significance is against those footprint-free tiles (Mann&ndash;Whitney, BH).")
    + f"<details><summary>Base resolution, prevalence, and the two classes</summary>"
      + fig(PT / "figures" / "p3_base_and_prevalence.png", "Base-resolution importance and partner prevalence",
            "Left: extra |effect| of a single-base change inside each factor's footprint, against bases with no footprint. Right: the share of sites carrying each factor within &plusmn;300 bp, by class.")
      + fig(PT / "figures" / "p2_family_by_class.png", "Family importance by class",
            "Median adjusted tile effect for Dynamic and context-matched Stable sites, by family. Dotted line: tiles with no footprint.")
      + "</details>"
    + bio("AlphaGenome has learned that a cardiomyocyte AP-1 site works with MEF2, and to a lesser extent GATA, NFI and TEAD &mdash; the core cardiac programme &mdash; and that the ubiquitous KLF/SP sites around it are largely decorative. This is an independent, motif-level confirmation that the model is doing cardiac biology rather than counting GC-rich sequence. Where a partner is present it carries the same weight in both classes, so partner identity does not explain responsiveness; the open question is whether Dynamic sites draw on a different repertoire (less MEF2, more CTCF) rather than a different number of partners."))

sec("spacing", "Spacing: is there a rule, and is there helical phasing?", "mid",
    "All 3,679 AP-1 sites and 95,506 merged footprints, so the distance work is well powered; the functional half rests on the scanned subset.",
    U("Distances are centre to centre between the AP-1 motif and each partner, over every cardiomyocyte AP-1 site rather than the scanned subset. Overlapping calls of one family at one locus are merged first, so redundant PWMs count once.",
      f"<b>The important partners are the close ones.</b> Median distance runs MEF2 {DNEAR.loc['MEF2', 'median_distance']:.0f} bp, NFI {DNEAR.loc['NFI (Nfia/b/c/x)', 'median_distance']:.0f}, CTCF {DNEAR.loc['CTCF / cohesin', 'median_distance']:.0f}, GATA {DNEAR.loc['GATA', 'median_distance']:.0f}, TEAD {DNEAR.loc['TEAD', 'median_distance']:.0f}, out to KLF/SP at {DNEAR.loc['KLF / SP', 'median_distance']:.0f} bp. The proximity order is almost exactly the importance order.",
      f"<b>There is a preferred abutting spacing.</b> Counting instances 11&ndash;25 bp away (just outside the 11 bp AP-1 motif) against the 50&ndash;150 bp baseline, MEF2 is enriched {DCLOSE.loc['MEF2', 'density_ratio_near_over_far']:.2f}-fold against {DCLOSE.loc['KLF / SP', 'density_ratio_near_over_far']:.2f} for KLF/SP (FDR {DCLOSE.loc['MEF2', 'fdr']:.3f}), with a modal centre-to-centre distance of {DMODE.loc['MEF2', 'modal_distance_9_40bp']:.0f} bp &mdash; the two motifs essentially touching, with a gap of a few bases. CTCF {DCLOSE.loc['CTCF / cohesin', 'density_ratio_near_over_far']:.2f}, NFI {DCLOSE.loc['NFI (Nfia/b/c/x)', 'density_ratio_near_over_far']:.2f} (mode {DMODE.loc['NFI (Nfia/b/c/x)', 'modal_distance_9_40bp']:.0f} bp), TEAD {DCLOSE.loc['TEAD', 'density_ratio_near_over_far']:.2f}, NKX {DCLOSE.loc['NKX / homeobox cardiac', 'density_ratio_near_over_far']:.2f} and AP-1 itself {DCLOSE.loc['AP-1 (Fos/Jun/Atf/Batf)', 'density_ratio_near_over_far']:.2f} are all enriched; GATA ({DCLOSE.loc['GATA', 'density_ratio_near_over_far']:.2f}, FDR {DCLOSE.loc['GATA', 'fdr']:.2f}), TBX and KLF/SP are not.",
      f"<b>But there is no helical phasing.</b> If the two proteins touched in a fixed register, the distance histogram should repeat every 10&ndash;11 bp (one turn of B-DNA). Nothing survives correction: MEF2 is the best case at raw p = {DPER.loc['MEF2', 'p_vs_jitter_null']:.3f} but FDR {DPER.loc['MEF2', 'fdr']:.2f}, and its 10&ndash;11 bp peak is only {DPER.loc['MEF2', 'obs_over_background_band']:.2f} of the strongest peak elsewhere between 5 and 30 bp. Pooling the four important families gives p = {DPER.loc['MEF2 + GATA + NFI + TEAD', 'p_vs_jitter_null']:.2f}.",
      f"<b>Closeness has a functional consequence in the model.</b> Median importance falls off with distance for the real partners &mdash; MEF2 {di('MEF2', '0-25'):.2f} at 0&ndash;25 bp against {di('MEF2', '131-200'):.2f} at 131&ndash;200 bp; NFI decays smoothly from {di('NFI (Nfia/b/c/x)', '0-25'):.2f} to {di('NFI (Nfia/b/c/x)', '201-500'):.2f} (rho {DRHO.loc['NFI (Nfia/b/c/x)', 'rho_distance_vs_importance']:+.2f}, FDR {DRHO.loc['NFI (Nfia/b/c/x)', 'fdr']:.2f}); CTCF {di('CTCF / cohesin', '0-25'):.2f} to {di('CTCF / cohesin', '201-500'):.2f}. KLF/SP is flat at about {di('KLF / SP', '0-25'):.2f} at every distance, which is the control.",
      f"<b>Class differences in spacing.</b> At Dynamic sites the NKX/homeobox partners sit much closer ({dc('NKX / homeobox cardiac', 'Dynamic vs Stable').median_a:.0f} vs {dc('NKX / homeobox cardiac', 'Dynamic vs Stable').median_b:.0f} bp, FDR {dc('NKX / homeobox cardiac', 'Dynamic vs Stable').fdr:.3f}) and CTCF closer too ({dc('CTCF / cohesin', 'Dynamic vs Stable').median_a:.0f} vs {dc('CTCF / cohesin', 'Dynamic vs Stable').median_b:.0f} bp, FDR {dc('CTCF / cohesin', 'Dynamic vs Stable').fdr:.2f}). Age Up (score) sites have their MEF2 partners further away ({dc('MEF2', 'Age Up (score) vs Unchanged').median_a:.0f} vs {dc('MEF2', 'Age Up (score) vs Unchanged').median_b:.0f} bp; {dc('MEF2', 'Age Up (score) vs Unchanged').share_within_50bp_a:.0%} vs {dc('MEF2', 'Age Up (score) vs Unchanged').share_within_50bp_b:.0%} within 50 bp, FDR {dc('MEF2', 'Age Up (score) vs Unchanged').fdr:.2f}) &mdash; the direction a motif-driven, partner-independent site would predict, but not significant.",
      "Caveat: adjacent motif calls can arise from one PWM partly matching its neighbour's sequence. That KLF/SP and TBX show no close-range enrichment argues against a generic artefact, but the abutting MEF2 peak should be confirmed on the underlying sequence."),
    fig(DR / "figures" / "r3_close_range.png", "Close-range spacing",
        "Share of each family's instances by centre-to-centre distance (1 bp bins, 3 bp rolling mean). Shaded: the overlap zone the 11 bp AP-1 motif occupies, and the 11&ndash;25 bp abutting window.")
    + f"<details><summary>Full distance profiles and the periodicity test</summary>"
      + fig(DR / "figures" / "r1_distance_profile.png", "Distance profiles",
            "Left: enrichment of each family against the pooled distance profile of all footprints, which carries the peak shape. Right: fine structure within 150 bp.")
      + fig(DR / "figures" / "r2_periodicity.png", "Periodicity",
            "Power spectrum of the detrended 1 bp distance histogram. Shaded: the 9.5&ndash;11.5 bp helical band. Dotted: the 95th percentile of nulls made by jittering each distance by up to &plusmn;7 bp.")
      + "</details>"
    + bio("The partners cluster just outside the AP-1 motif, with a mode around 13&ndash;16 bp between centres, and they lose influence as they move away. That is what cooperative binding on nearby DNA looks like. What is missing is the signature of a fixed protein-protein interface: no 10&ndash;11 bp repeat, so no evidence that AP-1 and MEF2 must sit on the same face of the helix. Read together, these AP-1 sites look like members of a compact cardiac enhancer module built by local co-occupancy rather than by a rigid composite element."))

sec("scanfeat", "What the mutation scan sees, across every class", "ok",
    f"All 3,678 cardiomyocyte elements, 16 scan-derived features, eight class contrasts, each matched on predicted context AND measured accessibility; BH within contrast.",
    U("This is the scan read as a whole rather than one feature at a time: for each site, how much of its dependence sits where, how concentrated it is, how far it spreads, and how much single bases matter. Earlier versions of this comparison ran on 399 sites and found nothing; with the complete scan the picture changes.",
      f"<b>Dynamic sites have more diffuse dependence.</b> {nsf('Dynamic vs Stable')} of 16 features separate them from matched Stable sites. The share of dependence within &plusmn;85 bp is {sf('Dynamic vs Stable', 'share of dependence within').median_class:.2f} vs {sf('Dynamic vs Stable', 'share of dependence within').median_reference:.2f} (FDR {sf('Dynamic vs Stable', 'share of dependence within').fdr:.4f}), the spread is {sf('Dynamic vs Stable', 'spread of dependence').median_class:.0f} vs {sf('Dynamic vs Stable', 'spread of dependence').median_reference:.0f} bp, and the concentration (Gini) is {sf('Dynamic vs Stable', GINI).median_class:.3f} vs {sf('Dynamic vs Stable', GINI).median_reference:.3f}.",
      f"<b>Single bases matter less, everywhere.</b> The largest single-base effect at a Dynamic site is {sf('Dynamic vs Stable', 'largest single-base').median_class:.2f} vs {sf('Dynamic vs Stable', 'largest single-base').median_reference:.2f}, importance inside the motif {sf('Dynamic vs Stable', 'single-base importance inside').median_class:.2f} vs {sf('Dynamic vs Stable', 'single-base importance inside').median_reference:.2f}, and in the flanks {sf('Dynamic vs Stable', 'single-base importance in the flanks').median_class:.2f} vs {sf('Dynamic vs Stable', 'single-base importance in the flanks').median_reference:.2f}. Both are lower, so this is not the motif losing ground to its flanks &mdash; the whole neighbourhood is flatter.",
      f"<b>Concentration is the axis that runs through every class, and it is directional.</b> Age Up (threshold) is the most diffuse of all ({sf('Age Up vs Stable (threshold)', GINI).standardised:+.2f} SD, FDR {sf('Age Up vs Stable (threshold)', GINI).fdr:.4f}), Exercise Down likewise ({sf('Ex Down vs Stable (threshold)', GINI).standardised:+.2f} SD). In the &Delta;score classes the sign flips between directions: Age Up {sf(C_AU_D, GINI).standardised:+.2f} SD (more diffuse) against Age Down {sf(C_AD_D, GINI).standardised:+.2f} SD (more concentrated), both FDR &lt; 0.01.",
      f"<b>Age Down (threshold) is the one class that looks like Stable</b> ({nsf('Age Down vs Stable (threshold)')} of 16 features differ), and Exercise Up nearly so ({nsf('Ex Up vs Stable (threshold)')} of 16).",
      "So gaining a binding call goes with dependence spreading out over the neighbourhood; losing one does not. A site that gains AP-1 occupancy with age is one whose accessibility was never pinned to a few bases in the first place."),
    fig(SF / "figures" / "y1_scan_features_by_class.png", "Scan features by class",
        "Standardised difference from the reference class for each scan-derived feature; * and a box mark FDR &lt; 0.05 within that contrast. Threshold classes left of the line, &Delta;score classes right.")
    + bio("Responsive sites are not built around a strong, tightly-defined core with a sharp dependence peak. Their dependence is spread thinly across the neighbourhood, no single base is worth much, and the AP-1 motif is one contributor among many. That is the same conclusion the partner work reached from the other direction, and it is why these sites sit near the occupancy threshold: there is no single anchor holding them open."))

sec("mpracheck", "Does the fibroblast MPRA agree with the scan?", "mid",
    "A ratio within each element, so assay scale and dynamic range cancel; but the MPRA is fibroblast and episomal while the scan is cardiomyocyte and genomic, and the MPRA's absolute effects are small.",
    U("The scan says the AP-1 motif is one contributor among many and no single base is worth much. The MPRA tests the same question in the wet lab, on 170 bp fragments outside the genome. The comparable quantity is the <b>share of an element's sequence dependence that the AP-1 motif accounts for</b>: in the MPRA (Endo &minus; MotifShuff) / (Endo &minus; FullShuff), in the model (replace the 170 bp &minus; replace the surroundings with the motif kept) / replace the 170 bp.",
      f"<b>The two agree closely.</b> AlphaGenome assigns the motif {MSA.loc[AGD, 'median_share']:.0%} of the local dependence (CI {MSA.loc[AGD, 'ci_lo']:.0%}&ndash;{MSA.loc[AGD, 'ci_hi']:.0%}, n = {int(MSA.loc[AGD, 'n']):,}); the MPRA assigns it {MSA.loc[MPD, 'median_share']:.0%} (CI {MSA.loc[MPD, 'ci_lo']:.0%}&ndash;{MSA.loc[MPD, 'ci_hi']:.0%}, n = {int(MSA.loc[MPD, 'n'])}). Two independent assays, different cell types, different readouts, and both put the AP-1 motif at under a tenth of what the surrounding 170 bp does.",
      "That is a real external check on the scan. The 'flat neighbourhood' result is not a quirk of the model: scrambling the AP-1 motif in a reporter barely changes its activity either.",
      f"<b>But the class direction disagrees, and sharply.</b> In the MPRA, elements in the responsive groups lose MORE activity when the whole element is scrambled than Stable ones ({msb('MPRA', 'whole element shuffled').dynamic:+.3f} vs {msb('MPRA', 'whole element shuffled').stable:+.3f} log2, p = {msb('MPRA', 'whole element shuffled').p:.4f}; Kruskal across all five groups p = {msb('MPRA', 'whole element shuffled').kruskal_5_groups_mpra:.4f}), with no difference in baseline activity (p = {msb('MPRA', 'baseline activity').p:.2f}). In the model, Dynamic sites lose LESS ({msb('AlphaGenome', 'whole element shuffled').dynamic:.2f} vs {msb('AlphaGenome', 'whole element shuffled').stable:.2f}).",
      "Those are not necessarily in conflict: the MPRA strips the element out of the genome, so it measures what the sequence can do on its own, while the model measures accessibility where chromatin already decides most of the outcome. A responsive element can be intrinsically informative as a fragment and still have its genomic accessibility set by trans factors rather than by its own sequence. The model's raw gap, after all, turned out to be chromatin context.",
      f"<b>Caveats.</b> The MPRA's absolute effects are small &mdash; destroying the entire 170 bp costs a median {MSA.loc[MPD, 'whole_element_effect']:.2f} log2 of activity, and the motif-only effect is {MSA.loc[MPD, 'motif_only_effect']:.3f}. The group comparison uses the MPRA's own classification in fibroblasts, not the cardiomyocyte TOBIAS classes, and the cross-map showed the two element sets only partly overlap."),
    fig(MS / "figures" / "m1_mpra_vs_scan.png", "MPRA against the mutation scan",
        "Left: the motif's share of each element's sequence dependence, measured independently by the model and by the reporter assay. Middle: the MPRA class direction. Right: the MPRA share by group against the model's estimate (dashed).")
    + bio("The most useful thing the MPRA does here is confirm the uncomfortable part of the scan. An AP-1 site contributes under a tenth of what its surrounding 170 bp contributes, in a model of chromatin and in a reporter assay alike. Where the two part company is on which elements are more sequence-driven, and that gap is informative in itself: it separates what a sequence can do in isolation from what decides its accessibility in the genome."))

sec("titration", "AP-1 affinity titration: does strength behave the way a binding site should?", "ok",
    f"All 3,678 elements; {int(TIO.loc[TIO.index[TIO.index.str.startswith('effect of changing 1')][0], 'n']):,} with an exact core carry the dose series, {int(TIO.loc[TIO.index[TIO.index.str.startswith('gain from upgrading')][0], 'n']):,} without one carry the upgrade. Each class contrast matched as above.",
    U(f"<b>The model behaves like a binding site, but saturates immediately.</b> Changing one core base already costs {TIO.loc[TIO.index[TIO.index.str.startswith('effect of changing 1')][0], 'median']:.2f} log2; changing all seven costs {TIO.loc[TIO.index[TIO.index.str.startswith('effect of changing all')][0], 'median']:.2f}. Almost the entire effect is spent on the first base, and the curve is flat from about four onward. There is no graded affinity ladder.",
      f"<b>Upgrading an imperfect core gains {TIO.loc[TIO.index[TIO.index.str.startswith('gain from upgrading')][0], 'median']:+.2f} log2</b> at the {int(TIO.loc[TIO.index[TIO.index.str.startswith('gain from upgrading')][0], 'n'])} sites that lack one &mdash; a real, sizeable gain, and the mirror of the loss above.",
      f"<b>Orientation does not matter:</b> reversing the whole 11 bp motif costs {TIO.loc[TIO.index[TIO.index.str.startswith('effect of reversing')][0], 'median']:.3f} log2, essentially nothing, as expected for a near-palindromic site. <b>A half-site is worth about as little as no site:</b> {TIO.loc[TIO.index[TIO.index.str.startswith('effect of keeping only')][0], 'median']:.2f} log2.",
      f"<b>A second copy genuinely adds accessibility</b> ({TIO.loc[TIO.index[TIO.index.str.startswith('net value of a second')][0], 'median']:+.2f} log2 against its own shuffled control) and <b>position matters more than presence</b>: the benefit peaks at &plusmn;20&ndash;40 bp and falls to near zero by &plusmn;150 bp. Relocating the motif behaves the same way. The model has learned that an AP-1 motif works where the rest of the element is, not anywhere in the peak.",
      f"<b>Validity check:</b> the dose-response slope tracks the core quality measured independently in the genome sequence (Spearman {TIV.rho.iloc[0]:.2f}, p &lt; 0.001), so the titration and the sequence scan agree.",
      f"<b>Between classes, almost nothing differs &mdash; with one exception.</b> Dynamic vs Stable: {nti('Dynamic vs Stable')} of 10. Every threshold class: 0 of 9&ndash;10. But <b>Age Up (&Delta;score) sites get more out of a second motif copy</b> ({ti(C_AU_D, 'net value of a second').median_class:+.2f} vs {ti(C_AU_D, 'net value of a second').median_reference:+.2f}, FDR {ti(C_AU_D, 'net value of a second').fdr:.4f}), benefit more when that copy is near ({ti(C_AU_D, 'does a second copy help').standardised:+.2f} SD, FDR {ti(C_AU_D, 'does a second copy help').fdr:.3f}), lose more to a single core base ({ti(C_AU_D, 'effect of changing 1 core').median_class:.2f} vs {ti(C_AU_D, 'effect of changing 1 core').median_reference:.2f}, FDR {ti(C_AU_D, 'effect of changing 1 core').fdr:.3f}) and tolerate relocation better.",
      "That is a coherent set: age-gained sites respond more to AP-1 dosage in both directions &mdash; they gain more from an added motif and lose more from a damaged one. It is the clearest support so far for a mass-action reading of the ageing change, and it comes from an experiment designed before that hypothesis existed."),
    fig(TI / "figures" / "z1_titration_overview.png", "Titration overview",
        "Left: dose-response to weakening the core, Stable against Dynamic, with 95% bootstrap CI. Middle: single edits. Right: the net value of a second or relocated motif copy against its own shuffled control, by offset.")
    + f"<details><summary>Every titration feature across every class</summary>{fig(TI / 'figures' / 'z2_titration_by_class.png', 'Titration by class', 'Standardised difference from the reference class; * and a box mark FDR &lt; 0.05 within that contrast.')}</details>"
    + bio("AP-1 motif strength behaves like a switch rather than a dial: one wrong base in the core costs nearly as much as destroying it, and a half-site is worth nothing. What is graded is <i>position</i> &mdash; a motif helps most within about 40 bp of where the element already is. Age-gained sites are the one class that reads AP-1 dosage more sharply in both directions, which is what you would expect if their change with age is driven by how much AP-1 protein is around rather than by how they are built."))

sec("motifmap", "Which motifs matter, broken down by DIRECTION of change", "inc",
    "Exploratory and thin: Age Down is 29 scanned sites and Exercise Up 22, so most cells rest on 5\u201320 sites. Every cell carries its n, and p-values are uncorrected.",
    U("Rather than Stable vs Dynamic, the classes here are directional. On the left of the heatmap the <b>threshold (binary call)</b> classes &mdash; Age Up, Age Down, Exercise Up, Exercise Down, each against Stable. On the right the <b>&Delta;score</b> classes against Unchanged. Importance is the median effect of replacing a 10 bp tile over that factor's footprint, measured against the same site's footprint-free tiles.",
      f"<b>The ranking holds in every direction.</b> MEF2, CTCF, NFI, TEAD and GATA lead the list in all five threshold classes. No direction of change recruits a different cast of factors.",
      f"<b>Age Up (threshold), more important than at Stable sites:</b> {top5(AGEUP_T, MORE)}.",
      f"<b>&hellip; and less important:</b> {top5(AGEUP_T, LESS)}.",
      f"<b>Age Down (threshold), more important:</b> {top5(AGEDN_T, MORE)}.",
      f"<b>&hellip; and less important:</b> {top5(AGEDN_T, LESS)}.",
      f"<b>Exercise Up (threshold), more important:</b> {top5(EXUP_T, MORE)}.",
      f"<b>Exercise Down (threshold), more important:</b> {top5(EXDN_T, MORE)}.",
      "<b>The pattern that recurs across directions is CTCF.</b> It is the largest difference for Age Up (threshold), Age Down (threshold) and Exercise Down, always in the same direction: architectural factors matter more at sites whose binding call changes, whichever way it changes. The Age Down cell is &minus;1.62 on six sites, so the size is not trustworthy, but the sign repeats across four independent contrasts.",
      "<b>MEF2 does the opposite and is direction-specific.</b> It is less important than Stable at every threshold class that gains or loses a call, and in the &Delta;score classes it splits: less important at Age Up, more important at Age Down. This is the fourth measure to show that split, after presence, distance and motif strength.",
      "<b>NF-kB and heat-shock factors come up at Age Up and Exercise Down</b> (both p &lt; 0.05 uncorrected). Given AP-1's role in stress signalling that is worth a look, but it rests on 10&ndash;22 sites and no correction.",
      "Treat every cell as a lead. Nothing here survives multiple-testing correction across the table, and the full scan will roughly quadruple the counts."),
    fig(MH / "figures" / "x1_motif_importance_heatmap.png", "Motif importance by direction of change",
        "Left: median importance within each class, n shown in each cell; threshold classes to the left of the black line, &Delta;score to the right. Right: difference from the reference class (Stable for threshold, Unchanged for &Delta;score), boxed where the family is in the top five either way, * = p &lt; 0.05 uncorrected. Colour is clipped so one six-site cell does not flatten the rest; printed values are exact.")
    + bio("Changing a binding call in any direction goes with CTCF mattering more and MEF2 mattering less &mdash; a site drifting toward architectural rather than tissue-identity support. The two age directions then separate on MEF2 specifically, matching everything else we have seen: age-lost sites are MEF2-anchored, age-gained sites are not. All of it is thin, and the full scan is the test."))

sec("answer", "The central question: does anything separate stable from changing?", "ok",
    f"One panel of 16 measures, four contrasts, identical matching on predicted accessibility, peak width, TSS distance and measured accessibility (balance within {SVCR.balance_sd.max():.3f} SD); BH-corrected within each contrast.",
    U("Everything tested so far is pulled into one comparison so the answer does not depend on which analysis is quoted. Measures fall in four blocks: the site's own AP-1 sequence, the sequence of its neighbourhood (PWM scan, needs no ATAC signal), its footprinted neighbourhood (needs signal), and how it responds to mutation.",
      f"<b>Yes &mdash; {nyes(TH)} of 16 measures separate Dynamic from Stable (threshold), but almost none of them is about how the site responds.</b> The site's own AP-1 motif is weaker ({sv(TH, 'best AP-1 PWM score within 100 bp').values}, FDR {sv(TH, 'best AP-1 PWM score within 100 bp').fdr:.3f}) and it has fewer AP-1 motifs around it ({sv(TH, 'AP-1 PWM hits within 100 bp').values}, FDR {sv(TH, 'AP-1 PWM hits within 100 bp').fdr:.3f}); both are sequence measurements and cannot be detection artefacts. Its neighbourhood carries fewer footprinted families ({sv(TH, 'TF families footprinted within 100 bp').values}, FDR {sv(TH, 'TF families footprinted within 100 bp').fdr:.4f}) and less MEF2 ({sv(TH, 'a MEF2 footprint within 100 bp').values}, FDR {sv(TH, 'a MEF2 footprint within 100 bp').fdr:.4f}).",
      f"<b>Of the four mutation responses, only the 170 bp shuffle differs</b> ({sv(TH, 'effect of shuffling the 170 bp').values}, FDR {sv(TH, 'effect of shuffling the 170 bp').fdr:.4f}) &mdash; the chromatin-context effect already established. Replacing the motif, replacing the peak and the 300 bp-away control are all flat, and the control being flat is what makes the rest credible.",
      f"<b>The two definitions point in OPPOSITE directions, because they select opposite kinds of site.</b> Dynamic (threshold) sites have fewer AP-1 motifs nearby; Changing (&Delta;score) sites have more ({sv(CH, 'AP-1 PWM hits within 100 bp').values}, FDR {sv(CH, 'AP-1 PWM hits within 100 bp').fdr:.3f}) and more footprinted families ({sv(CH, 'TF families footprinted within 100 bp').values}, FDR {sv(CH, 'TF families footprinted within 100 bp').fdr:.4f}). A threshold crossing picks out weak sites near the calling boundary; a large score change picks out strong sites with room to move.",
      f"<b>Direction-agnostic change explains little either way.</b> Changing vs Steady detects {nyes(CH)} of 16, and the largest is a {abs(sv(CH, 'effect of shuffling the 170 bp').difference):.2f} log2 shift in one mutation type.",
      f"<b>Age Down (&Delta;score) has a genuine sequence signature, and it is MEF2.</b> Those sites carry stronger MEF2 motifs in the DNA ({sv(AD, 'best MEF2 PWM score within 100 bp').values} best PWM score, FDR {sv(AD, 'best MEF2 PWM score within 100 bp').fdr:.3f}) as well as more MEF2 footprints ({sv(AD, 'a MEF2 footprint within 100 bp').values}, FDR {sv(AD, 'a MEF2 footprint within 100 bp').fdr:.3f}). Sequence and occupancy agree, so this is not a detection effect.",
      f"<b>Age Up (&Delta;score) detects nothing in this panel.</b> The core-quality result is still the largest effect for that contrast ({sv(AU, 'has an exact TGA[C/G]TCA core').values}, raw p = {sv(AU, 'has an exact TGA[C/G]TCA core').p:.3f}) but does not survive correction across 16 measures (FDR {sv(AU, 'has an exact TGA[C/G]TCA core').fdr:.2f}). Tested on its own it stands; tested inside a broad panel it does not, and both statements belong in the record."),
    fig(SVC / "figures" / "w1_stable_vs_changing.png", "Does anything separate stable from changing?",
        "Standardised difference for 16 measures under four class contrasts, all matched identically. Orange = FDR &lt; 0.05 within that contrast. Dotted lines separate the four blocks of measures.")
    + bio("The answer is yes, but not where the project has been looking. A changing AP-1 site is not built to respond differently &mdash; it has a weaker AP-1 motif, fewer AP-1 motifs around it, and a sparser, less MEF2-rich neighbourhood. In other words it is a marginal site in a marginal enhancer, which is exactly the kind of element whose occupancy should fluctuate with transcription-factor availability. Responsiveness looks like a consequence of being weak and poorly supported, not of a distinct regulatory design."))

sec("register", "Hypothesis register: everything tested, including what failed", "ok",
    "Eleven stated hypotheses, each with a matched test and a verdict recorded whether or not it held; two were overturned by their own controls.",
    U(f"{int((HREG.verdict == 'SUPPORTED').sum())} supported, {int((HREG.verdict == 'NOT SUPPORTED').sum())} not supported, 1 contradicted (the effect ran the other way) and 1 retracted after its control failed. Every comparison uses 1:1 matching on predicted accessibility at the motif and peak, peak width, TSS distance and measured accessibility unless stated.",
    f"<b>The firmest result is sequence-based.</b> Age-gained (&Delta;score) sites carry an exact TGA[C/G]TCA core {HREG.loc['H1', 'result'].split('%')[0]}% of the time against {HREG.loc['H1', 'result'].split('vs ')[1].split('%')[0]}% for matched controls (H1) &mdash; now across all 3,679 sites rather than the 78 scanned ones, and it survives matching on measured accessibility. Age-lost (&Delta;score) sites show nothing (H1b), and the threshold classes show nothing (H2).",
    "<b>Two hypotheses were killed by their own controls, and both had looked convincing.</b> AP-1 and a neighbouring MEF2 footprint appear strongly co-regulated across conditions (odds ratio 36) &mdash; but every other family in the peak agrees with AP-1 just as often, so it is peak-level signal, not co-regulation (H4 retracted by H4b). And the extra motif sensitivity of age-gained (&Delta;score) sites disappears once measured accessibility is matched as well as predicted (H8).",
    "<b>A control worth having:</b> the model's motif dependence tracks core quality closely (H6), and the sites that shrug off losing their motif are largely those without a perfect core (H9) &mdash; 53% of imperfect-core sites are non-responders against 21% of perfect-core sites. That gives the earlier 29% non-responder figure a concrete sequence explanation."),
    register_table()
    + fig(HYP / "figures" / "v2_register_and_controls.png", "Hypothesis register and the two controls that changed a verdict",
          "Left: every hypothesis and its verdict. Middle (H4b): how often each neighbouring family shares AP-1's three-condition binding pattern &mdash; MEF2 is not special. Right (H8): the Age Up motif effect with predicted context matched, and after matching measured accessibility too.")
    + fig(HYP / "figures" / "v1_hypotheses.png", "Core quality by class, and what the model sees",
          "Left: share with an exact core by class, against matched controls. Middle: the model's motif dependence by number of core mismatches. Right: whether core quality carries the Age Up effect.")
    + bio("Core quality is the one sequence feature that separates a class, and only for sites GAINING AP-1 binding with age. Everything else tested here &mdash; worse cores at Dynamic sites, homotypic AP-1 clustering, MEF2 co-regulation, partner occupancy &mdash; either failed or turned out to be a property of the peak rather than of the site. The matching caveat cuts both ways: measured accessibility is partly a consequence of the binding change being studied, so H8 may be over-adjustment rather than refutation, which is why the sequence result (H1) carries the weight."))

sec("context", "Is the neighbourhood of a Dynamic site different?", "mid",
    "Matched pairs on predicted context AND measured accessibility, 439 pairs, so footprint detection is matched by construction; the sequence scan is independent of detection altogether.",
    U("Earlier I reported that the partner deficit at Dynamic sites was a detection artefact. That test used 110 scanned sites split into accessibility quartiles, with as few as 12 sites per cell. <b>It was underpowered, and the conclusion was wrong.</b> Repeating it with 439 pairs matched on predicted accessibility, peak accessibility, peak width, TSS distance <i>and measured accessibility</i> (balance within 0.06 SD):",
      f"Dynamic sites have <b>fewer distinct TF families footprinted within 100 bp</b> (paired difference {cx('Dynamic vs Stable', 'distinct families within 100 bp').matched_median_diff:+.0f}, FDR {cx('Dynamic vs Stable', 'distinct families within 100 bp').fdr_matched:.3f}), and <b>MEF2 is footprinted beside them less often</b> ({CXB.loc['Dynamic vs Stable', 'matched_share_a']:.0%} vs {CXB.loc['Dynamic vs Stable', 'matched_share_b']:.0%}, McNemar p = {CXB.loc['Dynamic vs Stable', 'p_matched_mcnemar']:.3f}); where MEF2 is present it sits {cx('Dynamic vs Stable', 'distance to the nearest MEF2 (bp)').matched_median_diff:+.0f} bp further away (FDR {cx('Dynamic vs Stable', 'distance to the nearest MEF2 (bp)').fdr_matched:.3f}).",
      "Matching on accessibility is conservative here, not neutral. If MEF2 support causes the chromatin to be open, accessibility is a mediator and matching on it removes part of the effect. A difference that survives is therefore a floor, not an estimate.",
      f"<b>But the sequence says otherwise.</b> Scanning mm10 itself for consensus matches &mdash; which needs no ATAC signal and so cannot be biased by detection &mdash; there is <b>no difference in MEF2 motifs</b> within 100 bp between Dynamic sites and their matched partners (mean difference {sq('Dynamic vs Stable', 'MEF2 consensus matches within 100 bp').mean_diff:+.3f} motifs, p = {sq('Dynamic vs Stable', 'MEF2 consensus matches within 100 bp').p:.2f}), nor for GATA, TEAD, NFI, CTCF, NKX or KLF.",
      "<b>So the right statement is narrower:</b> Dynamic sites have the same MEF2 motifs in their DNA, but fewer of those motifs are footprinted. That is either detection that the matching has not fully removed, or genuine differential occupancy of equally-present motifs. Those two cannot be separated with this data, and the second would be biologically interesting.",
      f"<b>The two age classes run in OPPOSITE directions on MEF2, and this is consistent across every measure.</b> Age Down (&Delta;score) sites are MEF2-rich: a MEF2 footprint at {CXB.loc['Age Down (score) vs Unchanged', 'matched_share_a']:.0%} of them vs {CXB.loc['Age Down (score) vs Unchanged', 'matched_share_b']:.0%} of matched controls (McNemar p = {CXB.loc['Age Down (score) vs Unchanged', 'p_matched_mcnemar']:.3f}), their nearest cardiac partner closer ({cx('Age Down (score) vs Unchanged', 'distance to the nearest cardiac partner (bp)').matched_median_diff:+.1f} bp, FDR {cx('Age Down (score) vs Unchanged', 'distance to the nearest cardiac partner (bp)').fdr_matched:.3f}), and &mdash; in the sequence itself &mdash; stronger MEF2 motifs (best PWM score {sv(AD, 'best MEF2 PWM score within 100 bp').values}, FDR {sv(AD, 'best MEF2 PWM score within 100 bp').fdr:.3f}).",
      f"Age Up (&Delta;score) sites go the other way on all three: a MEF2 footprint at {CXB.loc['Age Up (score) vs Unchanged', 'matched_share_a']:.0%} vs {CXB.loc['Age Up (score) vs Unchanged', 'matched_share_b']:.0%} (McNemar p = {CXB.loc['Age Up (score) vs Unchanged', 'p_matched_mcnemar']:.3f}), MEF2 sitting further away, and a weaker best MEF2 PWM score ({sv(AU, 'best MEF2 PWM score within 100 bp').values}, raw p = {sv(AU, 'best MEF2 PWM score within 100 bp').p:.3f}). Only the Age Down direction survives correction across a broad panel, but the sign is the same in sequence, in distance and in occupancy, which is not what noise looks like.",
      f"In the sequence, Age Up (score) sites carry slightly more AP-1 consensus matches within 100 bp (p = {sq('Age Up (score) vs Unchanged', 'AP1 consensus matches within 100 bp').p:.3f}, FDR {sq('Age Up (score) vs Unchanged', 'AP1 consensus matches within 100 bp').fdr:.3f}) and are marginally GC-richer; Age Down sites are GC-poorer (FDR {sq('Age Down (score) vs Unchanged', 'GC content, central 200 bp').fdr:.3f}). The AP-1 enrichment agrees with the better consensus cores already seen."),
    bio("Two different kinds of site. A Dynamic AP-1 site has the same partner motifs written in its DNA as a Stable one, but fewer of them are occupied &mdash; it is a sequence-competent element that is not being used, which is exactly what a condition-labile element should look like. Age-lost (&Delta;score) sites are the opposite extreme: MEF2-rich, tightly packed cardiac enhancers. Age-gained (&Delta;score) sites have unremarkable neighbourhoods and better AP-1 motifs, so they look driven by the factor rather than by their company."))

sec("spread", "How many positions matter, and where across the peak?", "mid",
    "Context-matched scan groups; shape measures are normalised per site, so a class that loses less overall is not penalised.",
    U("From the tile scan, per site: how many 10 bp tiles matter, how concentrated the dependence is (Gini), how widely it spreads, how many separate runs it forms, and where its centre of mass sits relative to the peak rather than the motif.",
      f"<b>Against context-matched Stable sites, nothing differs.</b> Tiles that matter {sh('scan groups: dynamic vs matched_stable', 'tiles that matter (vs the far-tile threshold)').median_case:.0f} vs {sh('scan groups: dynamic vs matched_stable', 'tiles that matter (vs the far-tile threshold)').median_control:.0f} (FDR {sh('scan groups: dynamic vs matched_stable', 'tiles that matter (vs the far-tile threshold)').fdr:.2f}); tiles holding half the dependence {sh('scan groups: dynamic vs matched_stable', 'tiles holding half the dependence').median_case:.0f} vs {sh('scan groups: dynamic vs matched_stable', 'tiles holding half the dependence').median_control:.0f}; Gini {sh('scan groups: dynamic vs matched_stable', 'concentration of dependence (Gini)').median_case:.2f} vs {sh('scan groups: dynamic vs matched_stable', 'concentration of dependence (Gini)').median_control:.2f}; spread {sh('scan groups: dynamic vs matched_stable', 'spread of dependence (bp)').median_case:.0f} vs {sh('scan groups: dynamic vs matched_stable', 'spread of dependence (bp)').median_control:.0f} bp; position of the centre of mass in the peak {sh('scan groups: dynamic vs matched_stable', 'centroid position in peak half-widths').median_case:.2f} vs {sh('scan groups: dynamic vs matched_stable', 'centroid position in peak half-widths').median_control:.2f} half-widths. Every FDR above 0.3.",
      f"<b>Against high-accessibility Stable sites it does differ</b>, in the direction of a more diffuse element: spread {sh('scan groups: dynamic vs high_stable', 'spread of dependence (bp)').median_case:.0f} vs {sh('scan groups: dynamic vs high_stable', 'spread of dependence (bp)').median_control:.0f} bp (FDR {sh('scan groups: dynamic vs high_stable', 'spread of dependence (bp)').fdr:.3f}), share of dependence beyond 150 bp {sh('scan groups: dynamic vs high_stable', 'share of dependence beyond 150 bp').median_case:.2f} vs {sh('scan groups: dynamic vs high_stable', 'share of dependence beyond 150 bp').median_control:.2f} (FDR {sh('scan groups: dynamic vs high_stable', 'share of dependence beyond 150 bp').fdr:.3f}), share within 85 bp {sh('scan groups: dynamic vs high_stable', 'share of dependence within 85 bp').median_case:.2f} vs {sh('scan groups: dynamic vs high_stable', 'share of dependence within 85 bp').median_control:.2f}. But that group is not context-matched, so this tracks accessibility rather than class.",
      "Dependence is strikingly concentrated in every group: a Gini around 0.78, with half of a site's total dependence held by about 5 of 60 tiles."),
    fig(DEP / "figures" / "u1_sequence_and_shape.png", "Sequence motifs and the shape of dependence",
        "Left: consensus matches within 100 bp from a scan of mm10, Dynamic sites against their matched Stable partners. Right: how many tiles hold half of each site's total dependence.")
    + bio("Responsive and constitutive sites spend their dependence the same way: on the same number of positions, equally concentrated, equally spread, in the same place within the peak. Combined with the identical mutation-response profile, the architecture of a Dynamic element really is the same &mdash; what differs is how much is at stake and how much of its partner potential is actually in use."))

sec("tiles", "Tile scan: the site's own motif, or a neighbour?", "inc",
    f"{int(DTH.n_complete_sites.iloc[0])} of 1,470 sites finished; the one significant difference does not hold across readouts.",
    U(f"<b>Preliminary.</b> {int(DTH.n_complete_sites.iloc[0])} sites are done so far ({int(DYN.n_sites)} Dynamic, {int(MST.n_sites)} context-matched Stable, {int(HST.n_sites)} high-accessibility Stable). This section will be rerun when the scan finishes.",
      f"Every 10 bp tile from &minus;300 to +300 bp around the motif is replaced by random bases, 3 replicates per tile (replicate SD about {DTH.replicate_sd_median.iloc[0]:.2f} log2). For each site: how much the tiles covering the AP-1 motif matter, against the strongest other tile &mdash; a possible partner element.",
      f"<b>Neighbouring elements are common, in both classes.</b> At {DYN.share_partner_dominates:.0%} of Dynamic and {MST.share_partner_dominates:.0%} of matched Stable sites, some other 10 bp tile matters more than either motif tile (read over the 170 bp). This compares the strongest of about 58 tiles against the stronger of 2, so it overstates how often a partner truly dominates, but the sizes are comparable.",
      "One worked example (chr1:119,197,433): single-base changes 23&ndash;30 bp upstream cut predicted accessibility about 6-fold, except A&harr;T swaps at four of those positions. Those bases read TATTAATA, a MEF2 motif, whose AT-rich centre tolerates exactly those swaps. The model is treating a cardiac partner site 25 bp away as more important than the AP-1 motif itself.",
      f"Read over the 170 bp, the AP-1 motif tile matters less at Dynamic sites ({dt('170 bp', 'own_motif_tile').median_a:.2f} vs {dt('170 bp', 'own_motif_tile').median_b:.2f}, FDR {dt('170 bp', 'own_motif_tile').fdr:.2f}), but read over the motif itself the difference vanishes ({dt('motif', 'own_motif_tile').median_a:.2f} vs {dt('motif', 'own_motif_tile').median_b:.2f}, FDR {dt('motif', 'own_motif_tile').fdr:.2f}). Inconsistent between readouts, so not a finding.",
      f"Two hints, neither solid: Dynamic sites have slightly fewer strong partner tiles (median {dt('170 bp', 'n_strong_partner_tiles').median_a:.0f} vs {dt('170 bp', 'n_strong_partner_tiles').median_b:.0f}, p = {dt('170 bp', 'n_strong_partner_tiles').p_mannwhitney:.3f}, FDR {dt('170 bp', 'n_strong_partner_tiles').fdr:.2f}), and their strongest partner lies beyond &plusmn;85 bp more often ({DYN.share_partner_beyond_85bp:.0%} vs {MST.share_partner_beyond_85bp:.0%})."),
    fig(DPR / "figures" / "z3_own_motif_vs_partner.png", "Own motif vs strongest neighbouring element",
        "Left: per site, the effect of replacing the AP-1 motif tile against the strongest other 10 bp tile within &plusmn;300 bp (read over the 170 bp). Below the dashed line, a neighbour matters more. Right: the share of sites where a neighbour dominates, 95% bootstrap CI.")
    + f"<details><summary>Average tile profile and saturation mutagenesis</summary>"
      + fig(ST / "figures" / "s1_tile_profile.png", "Tile profile",
            "Median effect of each tile by position with 95% bootstrap CI, read over the motif (left) and the 170 bp (right).")
      + fig(ST / "figures" / "s2_saturation_mutagenesis.png", "Saturation mutagenesis",
            "Median |effect| of single-base changes from &minus;30 to +30 bp, and the share of that importance lying inside the 11 bp motif.")
      + "</details>"
    + bio("AlphaGenome often ranks a neighbouring element above the AP-1 motif itself, in both classes. Whether Dynamic sites lean on fewer or more distant partners is the open question the completed scan should answer."))

# ---------------------------------------------------------------- summary, limits, methods
ANSWER = f"""<div class="box"><p><b>The question:</b> in the AlphaGenome mutation results, are Dynamic AP-1 sites more or less sensitive to particular kinds of mutation than Stable ones?</p>
<p><b>The answer: no.</b> Raw, Dynamic sites look much less sensitive to broad edits ({PROF.loc['replace the 170 bp', 'dynamic_median']:.2f} vs {PROF.loc['replace the 170 bp', 'stable_median']:.2f} log2 for replacing the 170 bp). That difference is entirely explained by their weaker chromatin context. Matched on predicted accessibility, peak width and TSS distance, all {len(PROF)} mutation types agree within {(PROF.matched_dynamic_median - PROF.matched_stable_median).abs().max():.2f} log2, none significantly. Scaled to the same overall sensitivity, the mix of mutation types is the same, and a classifier on the full response profile is at chance (AUC {clf('shape only').auc:.2f}).</p>
<p><b>What this leaves:</b> responsive and constitutive AP-1 sites are built the same way and differ in how much accessibility is at stake. What makes a site responsive is therefore unlikely to be its local sequence architecture, and more likely something a baseline sequence model cannot see &mdash; AP-1 activity, co-factors, chromatin state.</p>
<p><b>Three things the same data do show:</b> mapping all 2.49 million cardiomyocyte TOBIAS footprints onto the tile scan, <b>MEF2 is the partner AlphaGenome treats as most important</b> near an AP-1 site ({FAMI.loc['MEF2', 'median_adj_effect']:.2f} log2 per tile against {FREE:.2f} for footprint-free tiles), ahead of GATA, NFI and TEAD, and far ahead of a second AP-1 site; those partners sit preferentially just outside the motif (mode {DMODE.loc['MEF2', 'modal_distance_9_40bp']:.0f} bp between centres) and lose influence with distance, but show no 10&ndash;11 bp helical phasing. And defining the age classes by <b>change in binding score</b> instead of crossing the threshold reveals what the threshold definition hides: sites that gain AP-1 binding with age depend more on the AP-1 motif itself ({AGEUP.loc['shuffle the motif', 'matched_diff']:+.2f} log2, FDR {AGEUP.loc['shuffle the motif', 'fdr_matched']:.3f}) and carry a perfect consensus core more often ({SCO.loc['Age Up (score)', 'share_with_exact_core']:.0%} vs {SCO.loc['Age Up (score)', 'share_control']:.0%}).</p></div>"""

LIMITS = U(
    "<b>A null result is not proof of no difference.</b> With 720 matched pairs, a true difference smaller than about 0.1 log2 in any single mutation type would likely be missed. What is excluded is a difference large enough to matter biologically, not any difference at all.",
    "<b>Matching is on predicted, not measured, accessibility</b> (plus peak width and TSS distance). If Dynamic sites differ in something correlated with all three, matching partly removes it. The same-peak comparison, which holds the whole locus constant, gives the same answer and is reassuring.",
    "<b>The model predicts baseline accessibility, not response.</b> It has no notion of condition, so it cannot be asked directly why a site responds to ageing or exercise. This is a limit of the question, not of the analysis.",
    "<b>Classes are threshold calls.</b> A site just below the footprint threshold in one condition is labelled Dynamic. That is why the same DNA is often Stable in one cell type and Dynamic in another.",
    "<b>The tile scan is 14% complete</b> and its one significant difference is inconsistent between readouts. Treat every number in that section as provisional.",
    "<b>The partner ranking uses TOBIAS footprints, not a motif scan.</b> A factor with no footprint call in cardiomyocytes cannot appear, and family labels come from gene-symbol patterns, so closely related motifs are pooled.",
    f"<b>The score-based Age Up result rests on one contrast of five</b> and on a {SCO.loc['Age Up (score)', 'n']:.0f}-site sequence check. It should be repeated on the full scan and, ideally, in fibroblasts.",
    "<b>Not yet done:</b> whether MEF2 partners sit preferentially near the AP-1 sites that gain binding with age, which would join the two new findings together.")

METHODS = U(
    "Cardiomyocyte Fos sites from the TOBIAS footprint table; heart ATAC track; AlphaGenome mouse, 131,072 bp windows, <code>model_all_folds</code>.",
    "Effect = mean over replicate edits of log2((edited + pc) / (native + pc)), pc = 5% of the median native value of that readout.",
    "Tests: paired Wilcoxon within matched pairs, Mann&ndash;Whitney between classes, BH-FDR within each family, 95% bootstrap CIs. Classifiers are logistic regression, 5-fold cross-validation with folds grouped by peak, checked against 100 label permutations.",
    "Built from saved tables by <code>Scripts/37_build_dynamic_report.py</code>; analyses in <code>Scripts/33_dynamic_vs_stable_deep_dive.py</code> and <code>Scripts/36_mutation_response_profiles.py</code>. The fuller report covering every experiment is <code>Analysis/report_cm.html</code>.")

toc = "<div class='toc'><ul>" + "".join(f"<li><a href='#{k}'>{i + 1}. {t}</a></li>" for i, (k, t, _) in enumerate(SEC)) + "<li><a href='#limits'>What would change this answer</a></li></ul></div>"
body = (f"<h1>Do Dynamic AP-1 sites respond differently to mutation than Stable ones?</h1>"
        f"<p class='sub'>Cardiomyocyte Fos sites, AlphaGenome in-silico perturbation &middot; {date.today():%d %b %Y}</p>"
        + ANSWER + "<h2>Contents</h2>" + toc + "".join(s[2] for s in SEC)
        + "<h2 id='limits'>What would change this answer</h2>" + LIMITS
        + "<h2>Methods in brief</h2>" + METHODS)

CSS = """:root{--bg:#fcfcfb;--panel:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--rule:#e4e3df;--blue:#2a78d6;--green:#0d7a54;--amber:#9a6700;--grey:#6b6a66}
@media (prefers-color-scheme:dark){:root{--bg:#1a1a19;--panel:#242422;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3a37;--blue:#3987e5;--green:#3fbf8f;--amber:#e0b341;--grey:#a8a79f}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:1000px;margin:0 auto;padding:32px 16px 80px}h1{font-size:1.9rem;line-height:1.2;margin:0 0 4px}
h2{font-size:1.3rem;margin:44px 0 8px;padding-top:12px;border-top:1px solid var(--rule)}h2 .num{color:var(--ink2);margin-right:8px}
p,li{max-width:78ch}.sub{color:var(--ink2);margin:0 0 20px}a{color:var(--blue)}
.box{background:var(--panel);border-radius:10px;padding:14px 18px;margin:18px 0}.box p{margin:8px 0}
figure{margin:18px 0}figcaption{color:var(--ink2);font-size:.9rem;margin-top:6px;max-width:80ch}
img{border-radius:8px;border:1px solid var(--rule);background:#fcfcfb;height:auto}
table.t{border-collapse:collapse;font-size:.86rem;margin:10px 0;display:block;overflow-x:auto;max-width:100%}
table.t th,table.t td{padding:4px 10px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap}
table.t th{color:var(--ink2);font-weight:600}table.t td:first-child,table.t th:first-child{text-align:left;white-space:normal;min-width:260px}
div.bio{border-left:4px solid var(--green);background:var(--panel);padding:10px 16px;margin:14px 0;border-radius:6px;max-width:80ch}
div.bio strong:first-child{color:var(--green)}
.conf{margin:6px 0 10px;font-size:.92rem;display:flex;gap:8px;align-items:baseline;max-width:80ch}.conf .badge{flex:none;margin-right:0}
.badge{display:inline-block;font-size:.78rem;font-weight:600;border-radius:99px;padding:1px 9px;border:1.5px solid currentColor;white-space:nowrap;margin-right:6px}
.ok,.yes{color:var(--green)}.mid{color:var(--amber)}.inc,.no{color:var(--grey)}.bad,.ret{color:#b3261e}
table.t.reg td{white-space:normal;vertical-align:top;text-align:left}table.t.reg td:nth-child(3){min-width:240px}table.t.reg td:nth-child(5){min-width:320px}
ul.pts{margin:6px 0 10px;padding-left:20px}ul.pts li{margin:5px 0}
.toc{font-size:.92rem;margin:8px 0}.toc ul{margin:0;padding-left:18px}.toc a{text-decoration:none}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2)}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:.88em}"""

# display names: the two class definitions must never be confusable in the rendered text
for _raw, _disp in (("Age Up (score)", "Age Up (&Delta;score)"), ("Age Down (score)", "Age Down (&Delta;score)"),
                    ("Exercise Up (score)", "Exercise Up (&Delta;score)"), ("Exercise Down (score)", "Exercise Down (&Delta;score)"),
                    ("Changing (score)", "Changing (&Delta;score)"), ("Steady (score)", "Steady (&Delta;score)")):
    body = body.replace(_raw, _disp)

doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
       f'<title>Dynamic vs Stable AP-1 sites</title><style>{CSS}</style></head><body><main>{body}</main></body></html>')
for k, p in IMGS.items():
    doc = doc.replace(f"@@IMG{k}@@", "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode())
out = AN / "report_dynamic_vs_stable.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB), {len(IMGS)} images, {len(SEC)} sections")
