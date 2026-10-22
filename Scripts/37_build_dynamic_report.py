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
PT, SCL = AN / "partners", AN / "score_classes"

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
    lab = {"ok": "&#10003; High confidence", "mid": "~ Moderate confidence", "inc": "? Inconclusive"}[kind]
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
      f"<b>Context matching</b>: each Dynamic site is paired with the Stable site closest in predicted accessibility (170 bp and peak), peak width and distance to the nearest TSS. {NP} of {ND} Dynamic sites find a partner within the caliper; after matching the four covariates differ by at most {BAL.std_diff_after.abs().max():.2f} SD.",
      "No TOBIAS footprint scores are used in any comparison. They define the classes and nothing else."),
    bio("Dynamic sites sit in weaker chromatin than Stable ones. Any comparison that does not hold that constant measures the chromatin, not the sequence. Every result below is reported raw and context-matched."))

sec("raw", "The raw difference, and why it is misleading", "ok",
    "Large, highly significant raw differences; the matched comparison uses 93% of Dynamic sites.",
    U(f"Raw, Dynamic sites lose less to every wide edit: replace the 170 bp {PROF.loc['replace the 170 bp', 'stable_median']:.2f} (Stable) vs {PROF.loc['replace the 170 bp', 'dynamic_median']:.2f} (Dynamic); shuffle the 170 bp {PROF.loc['shuffle the 170 bp', 'stable_median']:.2f} vs {PROF.loc['shuffle the 170 bp', 'dynamic_median']:.2f}; centred shuffle of 340 bp {PROF.loc['centred shuffle, 340 bp', 'stable_median']:.2f} vs {PROF.loc['centred shuffle, 340 bp', 'dynamic_median']:.2f} (all FDR &lt; 0.001).",
      f"Motif-level edits differ by only about 0.1: replace the motif {PROF.loc['replace the motif', 'stable_median']:.2f} vs {PROF.loc['replace the motif', 'dynamic_median']:.2f}; replace the 7 bp core {PROF.loc['replace the 7 bp core', 'stable_median']:.2f} vs {PROF.loc['replace the 7 bp core', 'dynamic_median']:.2f}.",
      "Edits far from the motif (&plusmn;300 and &plusmn;600 bp away) do not differ at all, which is the control: the two classes respond identically where nothing is at stake.",
      f"After matching on context, every one of the {len(PROF)} mutation types differs by at most {(PROF.matched_dynamic_median - PROF.matched_stable_median).abs().max():.2f} log2 and none is significant (smallest FDR {PROF.fdr_matched.min():.2f}). The 170 bp gap falls from {MATCH.loc['170 bp shuffle', 'raw_gap']:+.2f} raw to {MATCH.loc['170 bp shuffle', 'matched_gap_median']:+.2f} (CI {MATCH.loc['170 bp shuffle', 'matched_ci_lo']:+.2f} to {MATCH.loc['170 bp shuffle', 'matched_ci_hi']:+.2f})."),
    fig(DY / "figures" / "y5_mutation_type_profile.png", "Sensitivity to each mutation type",
        "Median effect of each mutation type for Stable (blue) and Dynamic (orange) sites with 95% bootstrap CI. Left marks: all sites (Mann&ndash;Whitney, BH). Right marks: context-matched pairs (paired Wilcoxon, BH).")
    + bio("The raw gap is a chromatin effect. A site in weaker chromatin has less accessibility to lose, so every broad edit takes less from it. Once that is held constant, the two classes lose the same amount to the same edits."))

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
      f"<b>The two definitions disagree a great deal.</b> Only {SAG.loc['age', 'threshold_Up_also_score_Up']:.0%} of threshold Age Up sites are also score Age Up, only {SAG.loc['age', 'score_Up_that_crossed_threshold']:.0%} of score Age Up sites ever crossed the threshold, and {SAG.loc['age', 'threshold_Stable_called_changing']:.0%} of threshold-Stable sites show a large score change.",
      f"<b>The score-based definition tracks the measured data about twice as well.</b> Against the measured snATAC change between the same conditions, the score change gives Spearman {SVM.loc['age (MAS vs YS)', 'rho_score_change_vs_measured']:.2f} (age) and {SVM.loc['exercise (MAE vs MAS)', 'rho_score_change_vs_measured']:.2f} (exercise), against {SVM.loc['age (MAS vs YS)', 'rho_threshold_class_vs_measured']:.2f} and {SVM.loc['exercise (MAE vs MAS)', 'rho_threshold_class_vs_measured']:.2f} for the threshold classes. Both definitions point the right way: threshold Age Up sites gain {bm('threshold', 'Age Up', 'median_measured_log2'):+.2f} log2 of measured accessibility, score Age Up sites {bm('score', 'Age Up (score)', 'median_measured_log2'):+.2f}, on {int(bm('score', 'Age Up (score)', 'n'))} against {int(bm('threshold', 'Age Up', 'n'))} sites.",
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
    f"{int(DTH.n_complete_sites.iloc[0])} scanned sites, {int(FAMI.n_sites.max())} of them carrying the commonest factor; the ranking is strong and consistent, the class comparison is underpowered.",
    U("Every footprinted cardiomyocyte TFBS in the full TOBIAS table (2.49 million instances, 786 motifs) was mapped onto the tile scan. For each 10 bp tile: which factors are bound there, and how much does replacing it cost? Effects are adjusted for the site's overall sensitivity and for distance from the AP-1 motif, and compared with tiles holding no footprint at all.",
      f"<b>MEF2 is the dominant partner, by a wide margin.</b> Replacing a tile over a MEF2 footprint costs {FAMI.loc['MEF2', 'median_adj_effect']:.2f} log2 (CI {FAMI.loc['MEF2', 'ci_lo']:.2f} to {FAMI.loc['MEF2', 'ci_hi']:.2f}; n = {int(FAMI.loc['MEF2', 'n_tiles'])} tiles over {int(FAMI.loc['MEF2', 'n_sites'])} sites) against about {FREE:.2f} for footprint-free tiles. Then GATA {FAMI.loc['GATA', 'median_adj_effect']:.2f}, TEAD {FAMI.loc['TEAD', 'median_adj_effect']:.2f}, CTCF {FAMI.loc['CTCF / cohesin', 'median_adj_effect']:.2f}.",
      f"<b>The abundant factors are not the important ones.</b> KLF/SP sites are the most common by far ({int(FAMI.loc['KLF / SP', 'n_tiles']):,} tiles) yet cost only {FAMI.loc['KLF / SP', 'median_adj_effect']:.2f}. Even other AP-1 instances nearby cost just {FAMI.loc['AP-1 (Fos/Jun/Atf/Batf)', 'median_adj_effect']:.2f} &mdash; a MEF2 neighbour matters roughly eight times more than a second AP-1 site.",
      f"The same ranking appears at base resolution: a single-base change inside a MEF2 footprint carries {FAMB.loc['MEF2', 'diff_vs_footprint_free']:+.2f} log2 more |effect| than a base with no footprint, ahead of CTCF {FAMB.loc['CTCF / cohesin', 'diff_vs_footprint_free']:+.2f}, TEAD {FAMB.loc['TEAD', 'diff_vs_footprint_free']:+.2f} and GATA {FAMB.loc['GATA', 'diff_vs_footprint_free']:+.2f}. The chr1:119,197,433 example was not a one-off.",
      f"Partners are not immediately adjacent: the median MEF2 footprint sits {PPOS.loc['MEF2', 'median_abs_distance']:.0f} bp from the AP-1 motif, and only {PPOS.loc['MEF2', 'share_within_50bp']:.0%} are within 50 bp.",
      f"<b>Do different factors matter to the two classes? Mostly no.</b> Of {len(FAMC)} families tested, only CTCF differs in importance after correction (Dynamic {FAMC.loc['CTCF / cohesin', 'median_adj_dynamic']:.2f} vs Stable {FAMC.loc['CTCF / cohesin', 'median_adj_stable']:.2f}, FDR {FAMC.loc['CTCF / cohesin', 'fdr_importance']:.3f}), on {int(FAMC.loc['CTCF / cohesin', 'n_tiles_dynamic'])} and {int(FAMC.loc['CTCF / cohesin', 'n_tiles_stable'])} tiles. MEF2 is present at fewer Dynamic sites ({FAMC.loc['MEF2', 'share_sites_dynamic']:.0%} vs {FAMC.loc['MEF2', 'share_sites_stable']:.0%}, FDR {FAMC.loc['MEF2', 'fdr_prevalence']:.2f}) but matters as much where it is present.",
      f"<b>A tempting result that is an artefact:</b> Dynamic sites appear to have fewer partners ({PCNT.loc['Dynamic', 'median_families_per_site']:.0f} vs {PCNT.loc['Stable, context-matched', 'median_families_per_site']:.0f} families per site, p = 0.001). Footprint calling needs ATAC signal, and the number of families correlates with measured accessibility at rho {PDET.rho_families_vs_measured.iloc[0]:.2f}. Within each quartile of measured accessibility the difference disappears entirely (p = {PDET.p.iloc[1]:.2f}, {PDET.p.iloc[2]:.2f}, {PDET.p.iloc[3]:.2f}, {PDET.p.iloc[4]:.2f}). It is detection, not biology."),
    fig(PT / "figures" / "p1_tile_importance_by_family.png", "Tile importance by TF family",
        "Median effect of replacing a 10 bp tile overlapping each factor's footprint, with 95% bootstrap CI, adjusted for the site's overall sensitivity and for distance from the AP-1 motif. Dashed line: tiles with no footprint. Significance is against those footprint-free tiles (Mann&ndash;Whitney, BH).")
    + f"<details><summary>Base resolution, prevalence, and the two classes</summary>"
      + fig(PT / "figures" / "p3_base_and_prevalence.png", "Base-resolution importance and partner prevalence",
            "Left: extra |effect| of a single-base change inside each factor's footprint, against bases with no footprint. Right: the share of sites carrying each factor within &plusmn;300 bp, by class.")
      + fig(PT / "figures" / "p2_family_by_class.png", "Family importance by class",
            "Median adjusted tile effect for Dynamic and context-matched Stable sites, by family. Dotted line: tiles with no footprint.")
      + "</details>"
    + bio("AlphaGenome has learned that a cardiomyocyte AP-1 site works with MEF2, and to a lesser extent GATA and TEAD &mdash; the core cardiac programme &mdash; and that the ubiquitous KLF/SP sites around it are largely decorative. This is an independent, motif-level confirmation that the model is doing cardiac biology rather than counting GC-rich sequence. The same partners matter in both classes, so partner identity does not explain why a site is responsive."))

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
<p><b>Two things the same data do show:</b> mapping all 2.49 million cardiomyocyte TOBIAS footprints onto the tile scan, <b>MEF2 is the partner AlphaGenome treats as most important</b> near an AP-1 site ({FAMI.loc['MEF2', 'median_adj_effect']:.2f} log2 per tile against {FREE:.2f} for footprint-free tiles), ahead of GATA and TEAD, and far ahead of a second AP-1 site. And defining the age classes by <b>change in binding score</b> instead of crossing the threshold reveals what the threshold definition hides: sites that gain AP-1 binding with age depend more on the AP-1 motif itself ({AGEUP.loc['shuffle the motif', 'matched_diff']:+.2f} log2, FDR {AGEUP.loc['shuffle the motif', 'fdr_matched']:.3f}) and carry a perfect consensus core more often ({SCO.loc['Age Up (score)', 'share_with_exact_core']:.0%} vs {SCO.loc['Age Up (score)', 'share_control']:.0%}).</p></div>"""

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
.ok{color:var(--green)}.mid{color:var(--amber)}.inc{color:var(--grey)}
ul.pts{margin:6px 0 10px;padding-left:20px}ul.pts li{margin:5px 0}
.toc{font-size:.92rem;margin:8px 0}.toc ul{margin:0;padding-left:18px}.toc a{text-decoration:none}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2)}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:.88em}"""

doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
       f'<title>Dynamic vs Stable AP-1 sites</title><style>{CSS}</style></head><body><main>{body}</main></body></html>')
for k, p in IMGS.items():
    doc = doc.replace(f"@@IMG{k}@@", "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode())
out = AN / "report_dynamic_vs_stable.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB), {len(IMGS)} images, {len(SEC)} sections")
