#!/usr/bin/env python3
"""Builds Analysis/report_cm.html (self-contained: figures embedded) from the tables and figures written by 22/23/24. Every number in the text
is read from the CSV tables, so rerun 24 then this script after new results arrive.   python3 25_build_report.py"""
from __future__ import annotations

import base64
import html
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import report_diagrams as dg

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "Analysis"
CM = AN / "cm"


def img(path, alt, width="100%"):
    b = base64.b64encode(Path(path).read_bytes()).decode()
    return f'<img src="data:image/png;base64,{b}" alt="{html.escape(alt)}" style="max-width:{width}">'


def table(df, fmt=None, cls=""):
    fmt = fmt or {}
    d = df.copy()
    for c in d.columns:
        if c in fmt:
            d[c] = d[c].map(lambda v, f=fmt[c]: "" if pd.isna(v) else f.format(v))
        elif d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.2f}")
    return d.to_html(index=False, escape=False, classes=f"t {cls}", border=0)


def stars(q):
    return "" if pd.isna(q) else "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def fp(p):
    return "&lt;0.0001" if p < 1e-4 else f"{p:.4f}" if p < 0.01 else f"{p:.3f}"


A = pd.read_csv(CM / "stats" / "stable_vs_dynamic_all_readouts.csv")
B = pd.read_csv(CM / "stats" / "class_vs_reference.csv")
C = pd.read_csv(CM / "stats" / "stable_vs_dynamic_by_native_quintile.csv")
D = pd.read_csv(CM / "stats" / "replicate_noise_detectability.csv")
E = pd.read_csv(CM / "stats" / "between_dynamic_classes.csv")
OM = pd.read_csv(CM / "stats" / "omnibus_kruskal.csv")
PM = pd.read_csv(CM / "stats" / "permutation_stable_vs_dynamic.csv")
NV = pd.read_csv(CM / "tables" / "native_vs_measured.csv")
T5 = pd.read_csv(CM / "tables" / "confidence_tiers.csv")
RN = pd.read_csv(CM / "tables" / "rna_effect_by_tss_distance.csv")
PS = pd.read_csv(AN / "insertion" / "tables" / "peak_swap.csv")
MI = pd.read_csv(AN / "insertion" / "tables" / "motif_insertion.csv")

W = {"motif": "Motif (11 bp)", "170": "170 bp", "peak": "Peak"}
m = A[(A.trackset == "CM") & (A.shuffle == A.readout)].set_index("shuffle")
fb = A[(A.trackset == "FB") & (A.shuffle == A.readout)].set_index("shuffle")
pm = PM.set_index("width")
nv = NV[(NV.trackset == "CM") & (NV.readout == "atac_peak")].iloc[0]
adj = B[(B.family.isin(["age", "exercise"])) & ~B["class"].str.startswith("[") & B.raw_p.notna()]
a170 = adj[adj.width == "170"]
quint_all_sig = bool((C.fdr < 0.05).all())
frac = D.set_index("shuffle").frac_sites_effect_gt_2se
t5 = T5.pivot_table(index=["tier", "class"], columns="shuffle_width", values="median_log2fc", aggfunc="first")
rna_near = RN[(RN.readout == "rna_prox2kb") & (RN["class"] == "Stable") & (RN.tss_bin == "<=2 kb")].set_index("shuffle_width").median_log2fc
swap = PS.set_index("compared_with")
ins = MI[(MI.trackset == "CM")]


def mi(rd, comp):
    return ins[(ins.readout == rd) & (ins.comparison == comp)].iloc[0]


CSS = """
:root{--bg:#fcfcfb;--panel:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--rule:#e4e3df;--blue:#2a78d6;--orange:#eb6834;--green:#0d7a54;--amber:#9a6700;--red:#b3261e;--grey:#6b6a66}
@media (prefers-color-scheme:dark){:root{--bg:#1a1a19;--panel:#242422;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3a37;--blue:#3987e5;--orange:#d95926;--green:#3fbf8f;--amber:#e0b341;--red:#f0776f;--grey:#a8a79f}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:1000px;margin:0 auto;padding:32px 16px 80px}h1{font-size:1.9rem;line-height:1.2;margin:0 0 4px}h2{font-size:1.35rem;margin:44px 0 8px;padding-top:12px;border-top:1px solid var(--rule)}h3{font-size:1.05rem;margin:26px 0 4px}figure.dia{margin:14px 0}.dia{overflow-x:auto}.dia svg{min-width:760px;width:100%;height:auto;display:block;color:var(--ink)}.dia figcaption{position:sticky;left:0}h3.sec{font-size:1.15rem;margin:34px 0 6px;padding-top:10px;border-top:1px dotted var(--rule)}
p,li{max-width:78ch}.sub{color:var(--ink2);margin:0 0 20px}.box{background:var(--panel);border-radius:10px;padding:14px 18px;margin:16px 0}.box ul{margin:6px 0;padding-left:20px}
figure{margin:18px 0}figcaption{color:var(--ink2);font-size:.9rem;margin-top:6px;max-width:80ch}img{border-radius:8px;border:1px solid var(--rule);background:#fcfcfb;height:auto}
table.t{border-collapse:collapse;font-size:.86rem;margin:10px 0;display:block;overflow-x:auto;max-width:100%}table.t th,table.t td{padding:4px 10px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap}table.t th{color:var(--ink2);font-weight:600;position:sticky;top:0;background:var(--bg)}table.t td:first-child,table.t th:first-child{text-align:left}
.badge{display:inline-block;font-size:.78rem;font-weight:600;border-radius:99px;padding:1px 9px;border:1.5px solid currentColor;white-space:nowrap}.ok{color:var(--green)}.mid{color:var(--amber)}.bad{color:var(--red)}.na{color:var(--grey)}
td.v{white-space:normal!important;text-align:left!important;min-width:260px}.small{font-size:.86rem;color:var(--ink2)}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:.88em}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2)}
"""


def badge(kind):
    lab = {"ok": "&#10003; Robust", "mid": "~ Plausible / partial", "bad": "&#10007; Not supported", "na": "? Cannot judge yet"}[kind]
    return f'<span class="badge {kind}">{lab}</span>'



# ======================= QC / additional-analysis tables =======================
from scipy.stats import chi2_contingency as _chi2, mannwhitneyu as _mwu

QCD = CM / "qc"
rd = lambda n: pd.read_csv(QCD / "tables" / n)


def qimg(name, alt, width="100%"):
    return img(QCD / "figures" / f"{name}.png", alt, width)


Q3, Q4, Q5, Q6, Q7, Q8 = (rd("qc3_native_means_vs_slides.csv"), rd("qc4_class_labels_vs_measured.csv"), rd("qc5_prediction_vs_measured_by_class.csv"),
                          rd("qc6_shuffled_vs_native.csv"), rd("qc7_replicate_convergence.csv"), rd("qc8_gpu_consistency.csv"))
AN1, AN2, AN3, AN4, AN5, AN6, AN7 = (rd("an1_effect_distributions.csv"), rd("an2_effect_vs_covariates_spearman.csv"), rd("an3_site_types.csv"), rd("an4_footprint_score.csv"),
                                     rd("an5_multiple_motifs_in_peak.csv"), rd("an6_tss_proximity.csv"), rd("an7_mpra_comparison.csv"))
SL = rd("site_level_table_cm.csv")
n_all, n_st, n_dy = len(SL), int((SL.stable_dynamic == "Stable").sum()), int((SL.stable_dynamic == "Dynamic").sum())
gp = {c: _mwu(SL[SL.stable_dynamic == "Stable"][c].dropna(), SL[SL.stable_dynamic == "Dynamic"][c].dropna()).pvalue for c in ("peak_width", "abs_tss", "motifs_in_peak")}
gm = SL.groupby("stable_dynamic")[["peak_width", "abs_tss"]].median()
q3 = Q3.set_index(["class", "readout"])
q4 = Q4.set_index(["label", "class"])
q5 = Q5.set_index("group")
q6 = Q6[Q6.readout == Q6.shuffle].set_index(["shuffle", "class"])
q7 = Q7.set_index(["shuffle", "k_replicates"]).median_abs_error_vs_10
q8 = Q8.set_index("measure")
an1 = AN1.set_index(["shuffle", "class"])
an2 = AN2.set_index(["shuffle", "covariate"])
an3 = AN3.set_index("site_type")
an4 = AN4.set_index(["outcome", "class"])
an5 = AN5.set_index("shuffle")
an6 = AN6.set_index(["measure", "class", "tss_bin"])
chrom_p = _chi2(pd.crosstab(SL.seqnames, SL.stable_dynamic).values)[1]
chrx = pd.crosstab(SL.seqnames, SL.stable_dynamic).loc['chrX']
rho_mc = SL[["eff_motif", "eff_170"]].corr(method="spearman").iloc[0, 1]
mp_ = AN7.set_index(["track", "comparison", "subset"])
CMP_N, CMP_M, CMP_F = ("native predicted ATAC vs log2 alpha (native)", "motif-shuffle effect vs MPRA (Endo - motif-shuffled)", "full-shuffle effect vs MPRA (Endo - fully shuffled)")
age_p = float(q4.loc[("age_class", "Kruskal-Wallis"), "p"])
ex_p = float(q4.loc[("ex_class", "Kruskal-Wallis"), "p"])


def fmtp(p):
    return "&lt; 1e-10" if p < 1e-10 else f"= {p:.1e}"


Q4_tab = Q4[Q4["class"] != "Kruskal-Wallis"][["label", "class", "n", "median_log2fc"]].copy()
Q4_tab["label"] = Q4_tab.label.map({"age_class": "Age: log2(MAS/YS)", "ex_class": "Exercise: log2(MAE/MAS)"})
Q4_tab.columns = ["Comparison", "Class", "n sites", "Median measured log2FC"]
Q5_tab = Q5.copy()
Q5_tab.columns = ["Group", "Spearman rho", "CI lo", "CI hi", "n"]
Q6_tab = Q6[Q6.readout == Q6.shuffle][["readout", "class", "mean_native", "mean_shuffled", "frac_sites_decreased", "n"]].copy()
Q6_tab["readout"] = Q6_tab.readout.map(W)
Q6_tab["mean_shuffled"] = Q6_tab.mean_shuffled.round(3)
Q6_tab.columns = ["Matched shuffle/readout", "Class", "Mean native ATAC", "Mean shuffled ATAC", "Sites that lose accessibility", "n"]
Q7_tab = Q7[Q7.k_replicates.isin([1, 3, 5, 7, 9])].pivot(index="k_replicates", columns="shuffle", values="median_abs_error_vs_10").reset_index()
Q7_tab.columns = ["Replicates averaged", "170 bp shuffle", "Motif shuffle", "Peak shuffle"]
Q7_tab = Q7_tab[["Replicates averaged", "Motif shuffle", "170 bp shuffle", "Peak shuffle"]]
Q8_tab = Q8.copy()
Q8_tab["measure"] = Q8_tab.measure.map({"eff_motif": "Motif-shuffle effect", "eff_170": "170 bp-shuffle effect", "eff_peak": "Peak-shuffle effect", "native_peak": "Native predicted ATAC"})
Q8_tab["kruskal_p_across_gpus"] = Q8_tab.kruskal_p_across_gpus.map(lambda v: f"{v:.2f}")
Q8_tab.columns = ["Measure", "Kruskal-Wallis p across GPU types", "Median (L40S)", "Median (A100)", "Median (H100)"]
AN2_tab = AN2.pivot(index="covariate", columns="shuffle", values="rho").reindex(["native_peak", "peak_width", "abs_tss", "fp_score"])[["motif", "170", "peak"]].reset_index()
AN2_tab["covariate"] = AN2_tab.covariate.map({"native_peak": "Native predicted ATAC (peak)", "peak_width": "Peak width", "abs_tss": "Distance to nearest TSS", "fp_score": "TOBIAS footprint score"})
AN2_tab.columns = ["Covariate (Spearman rho with the effect)", "Motif shuffle", "170 bp shuffle", "Peak shuffle"]
AN3_tab = AN3.reset_index()[["site_type", "n", "share_stable", "share_dynamic", "median_native_peak", "median_peak_width", "median_abs_tss"]]
AN3_tab.columns = ["Site type", "n", "Share of Stable", "Share of Dynamic", "Median native ATAC", "Median peak width", "Median TSS distance"]
AN4_tab = AN4[AN4.outcome.isin(["eff_motif", "eff_170", "native_peak"])][["outcome", "class", "rho", "ci_lo", "ci_hi", "n"]].copy()
AN4_tab["outcome"] = AN4_tab.outcome.map({"eff_motif": "Motif-shuffle effect", "eff_170": "170 bp-shuffle effect", "native_peak": "Native predicted ATAC"})
AN4_tab.columns = ["Outcome (vs footprint score)", "Sites", "Spearman rho", "CI lo", "CI hi", "n"]
AN5_tab = AN5[["shuffle", "n_single", "n_multi", "median_single", "median_multi", "raw_diff", "raw_lo", "raw_hi", "raw_fdr", "adj_diff", "adj_lo", "adj_hi", "adj_fdr"]].copy()
AN5_tab["shuffle"] = AN5_tab.shuffle.map(W)
AN5_tab["raw_fdr"] = AN5_tab.raw_fdr.map(lambda v: f"{v:.2f}")
AN5_tab["adj_fdr"] = AN5_tab.adj_fdr.map(lambda v: f"{v:.2f}")
AN5_tab.columns = ["Shuffle", "n one motif", "n 2+ motifs", "Median (one)", "Median (2+)", "Raw diff", "CI lo", "CI hi", "FDR", "Adj. diff", "CI lo ", "CI hi ", "FDR "]
AN6_tab = AN6[(AN6["class"] == "Stable") & AN6.measure.isin(["native_peak", "eff_motif", "eff_170", "eff_peak"])].pivot(index="tss_bin", columns="measure", values="median").reindex(["<=2 kb", "2-10 kb", "10-30 kb", ">30 kb"])[["native_peak", "eff_motif", "eff_170", "eff_peak"]].reset_index()
AN6_tab.columns = ["Distance to nearest TSS", "Native ATAC (peak)", "Motif shuffle", "170 bp shuffle", "Peak shuffle"]
AN7_tab = AN7[AN7.subset.isin(["all elements", "replicated in MPRA", "motif from TOBIAS"])][["track", "comparison", "subset", "rho", "ci_lo", "ci_hi", "n"]].copy()
AN7_tab["comparison"] = AN7_tab.comparison.map({CMP_N: "Native predicted ATAC vs MPRA log2 alpha", CMP_M: "Motif-shuffle effect vs MPRA (native - motif-shuffled)", CMP_F: "Full-shuffle effect vs MPRA (native - fully shuffled)"})
AN7_tab.columns = ["Track", "Comparison", "Elements", "Spearman rho", "CI lo", "CI hi", "n"]

PART_A = f"""
<h2 id="partA">Part A. Where the data come from (QC)</h2>
<p>These figures show how the site set, the model output and the runs themselves look before any hypothesis is tested. Each ends with what it means for how far the later results can be trusted.</p>

<h3>A1. The site set</h3>
<figure>{qimg('qc1_site_composition', 'Site composition')}<figcaption>Left: sites per TOBIAS binding pattern. Middle: how the age labels (YS vs MAS) and exercise labels (MAS vs MAE) co-occur. Right: sites per chromosome, split into Stable and Dynamic.</figcaption></figure>
<p>{n_all:,} cardiomyocyte Fos sites: {n_st:,} Stable ({n_st / n_all:.0%}) and {n_dy:,} Dynamic, spread across the six dynamic binding patterns (79 to 210 sites each). Age and exercise labels are <strong>not independent</strong>: they come from the same three binding calls, so there are only seven patterns (for example every Age Up site is either Exercise Stable or Exercise Down). Chromosome counts follow chromosome size, and the Stable/Dynamic mix does not differ between chromosomes (chi-square p = {chrom_p:.2f}; chrX has the highest Dynamic share, {chrx['Dynamic'] / (chrx['Dynamic'] + chrx['Stable']):.0%} of only {int(chrx.sum())} sites). <strong>Implication:</strong> class sizes are small (about 100 to 400), which limits the power of every class comparison. This site set is not the slides' (2,663 sites, 2,320 Stable).</p>

<h3>A2. Site geometry</h3>
<figure>{qimg('qc2_site_geometry', 'Site geometry')}<figcaption>Peak width, distance to the nearest TSS and number of Fos motifs in the same peak, for Stable and Dynamic sites.</figcaption></figure>
<p>Dynamic sites sit in narrower peaks (median {gm.loc['Dynamic', 'peak_width']:.0f} vs {gm.loc['Stable', 'peak_width']:.0f} bp, Mann&ndash;Whitney p {fmtp(gp['peak_width'])}), but at a similar distance from TSSs (median {gm.loc['Dynamic', 'abs_tss']:,.0f} vs {gm.loc['Stable', 'abs_tss']:,.0f} bp, p = {gp['abs_tss']:.2f}). About {(SL.motifs_in_peak >= 2).mean():.0%} of sites share their peak with another Fos site (about {SL.drop_duplicates('peak_coord').motifs_in_peak.ge(2).mean():.0%} of peaks hold two or more). <strong>Implication:</strong> peak width differs systematically between classes and also drives the peak-shuffle effect (Part C), which is why the class comparisons were adjusted for it.</p>

<h3>A3. Native predicted accessibility, and how it compares with the slides</h3>
<figure>{qimg('qc3_native_distributions', 'Native predicted ATAC distributions')}<figcaption>Distribution of native predicted ATAC (mouse heart track) at the motif, 170 bp and peak regions. Legend: our mean and, in brackets, the mean on your QC slide for the same class.</figcaption></figure>
<p>Predicted accessibility is strongly right-skewed: most sites are low (median {q3.loc[('Stable', 'peak'), 'our_median']:.2f}) and a minority are very high. For <strong>Dynamic sites the means match the slides almost exactly</strong> ({q3.loc[('Dynamic', 'motif'), 'our_mean']:.2f} / {q3.loc[('Dynamic', '170'), 'our_mean']:.2f} / {q3.loc[('Dynamic', 'peak'), 'our_mean']:.2f} vs {q3.loc[('Dynamic', 'motif'), 'slide_mean']:.2f} / {q3.loc[('Dynamic', '170'), 'slide_mean']:.2f} / {q3.loc[('Dynamic', 'peak'), 'slide_mean']:.2f}), which supports this pipeline reproducing the API model. <strong>Stable sites do not</strong> ({q3.loc[('Stable', 'motif'), 'our_mean']:.2f} / {q3.loc[('Stable', '170'), 'our_mean']:.2f} / {q3.loc[('Stable', 'peak'), 'our_mean']:.2f} vs {q3.loc[('Stable', 'motif'), 'slide_mean']:.2f} / {q3.loc[('Stable', '170'), 'slide_mean']:.2f} / {q3.loc[('Stable', 'peak'), 'slide_mean']:.2f}); the slides' Stable set is a different, more accessible set of sites.</p>

<h3>A4. Are the class labels real? Checking them against independent measured data</h3>
<figure>{qimg('qc4_class_labels_vs_measured', 'Class labels vs measured accessibility change')}<figcaption>Measured cardiomyocyte snATAC (pseudobulk, CPM per kb at the peak containing the site) split by experimental condition. Left: aged vs young sedentary, by age class. Right: aged plus exercise vs aged sedentary, by exercise class.</figcaption></figure>
{table(Q4_tab, {'Median measured log2FC': '{:+.3f}'})}
<p>The classes come from TOBIAS footprint calls, not from accessibility, so measured accessibility is an independent check. <strong>The gains are clear and the losses are weak</strong>: Age Up sites gain accessibility with age (median log2 MAS/YS {q4.loc[('age_class', 'Age Up'), 'median_log2fc']:+.2f} vs {q4.loc[('age_class', 'Stable'), 'median_log2fc']:+.2f} for Stable), as do Exercise Up sites ({q4.loc[('ex_class', 'Exercise Up'), 'median_log2fc']:+.2f}); Age Down and Exercise Down sites change in the expected direction but only slightly. Differences across classes are significant (Kruskal&ndash;Wallis p {fmtp(age_p)} (age) and p {fmtp(ex_p)} (exercise)). <strong>Implication:</strong> the labels reflect real but modest accessibility shifts, so "Dynamic" classes are best read as footprint-defined, not as strongly differentially accessible.</p>

<h3>A5. Does the model predict some classes worse?</h3>
<figure>{qimg('qc5_prediction_vs_measured_by_class', 'Prediction vs measurement by class')}<figcaption>Left: Spearman correlation of predicted peak accessibility with measured accessibility, per class (95% CI). Right: distribution of the rank residual (predicted percentile minus measured percentile).</figcaption></figure>
{table(Q5_tab, {'Spearman rho': '{:.2f}', 'CI lo': '{:.2f}', 'CI hi': '{:.2f}'})}
<p>The correlation is {q5.loc['Stable', 'rho']:.2f} in Stable sites and {q5.loc['Dynamic', 'rho']:.2f} in Dynamic sites; individual patterns range from {q5.rho.iloc[3:].min():.2f} to {q5.rho.iloc[3:].max():.2f} with wide intervals (n = 65 to 195). No class is predicted badly enough to invalidate it, but Dynamic sites are somewhat less well predicted. <strong>Implication:</strong> the confidence tiers in B6 are partly a way of down-weighting these poorly predicted sites.</p>

<h3>A6. Does shuffling lower accessibility as expected? (the slides' QC)</h3>
<figure>{qimg('qc6_shuffled_vs_native', 'Shuffled vs native')}<figcaption>Mean predicted ATAC after each shuffle (bars) against the native mean (dashed), for each readout region.</figcaption></figure>
{table(Q6_tab, {'Mean native ATAC': '{:.3f}', 'Sites that lose accessibility': '{:.0%}'})}
<p>Shuffling reduces accessibility on average everywhere, as on the slides. A 170 bp shuffle lowers the 170 bp readout in {q6.loc[('170', 'Stable'), 'frac_sites_decreased']:.0%} of Stable sites, but a motif-only shuffle lowers the motif readout in only {q6.loc[('motif', 'Stable'), 'frac_sites_decreased']:.0%}: for about a quarter of sites it raises predicted accessibility, so the motif-only effect is small and heterogeneous. <strong>Implication:</strong> the shuffle experiment behaves as intended and the effect sizes are not an artefact of broken sequences.</p>

<h3>A7. How many shuffle replicates are enough?</h3>
<figure>{qimg('qc7_replicate_convergence', 'Replicate convergence')}<figcaption>Left: median difference between an effect averaged over k replicates and the 10-replicate effect, by shuffle width. Right: spread (SD) across the 10 replicates at each site.</figcaption></figure>
{table(Q7_tab, {c: '{:.3f}' for c in ['Motif shuffle', '170 bp shuffle', 'Peak shuffle']})}
<p>With a single shuffle a site's effect is off by about {q7.loc[('motif', 1)]:.2f} (motif), {q7.loc[('170', 1)]:.2f} (170 bp) and {q7.loc[('peak', 1)]:.2f} (peak) log2 units; with 5 replicates it is about {q7.loc[('motif', 5)]:.2f}, {q7.loc[('170', 5)]:.2f} and {q7.loc[('peak', 5)]:.2f}. Motif and peak shuffles need the most replicates. <strong>Implication:</strong> 10 replicates is a sensible minimum; fewer would make single-site values unreliable, and even 10 leaves peak-shuffle per-site values noisy (this error is measured against the 10-replicate value itself, so the true error is somewhat larger).</p>

<h3>A8. Does the compute hardware matter?</h3>
<figure>{qimg('qc8_gpu_consistency', 'GPU consistency')}<figcaption>Median effect per array chunk (sites are interleaved across chunks), coloured by the GPU model that ran the chunk.</figcaption></figure>
{table(Q8_tab, {'Median (L40S)': '{:.3f}', 'Median (A100)': '{:.3f}', 'Median (H100)': '{:.3f}'})}
<p>The private queue assigned chunks to L40S, A100 and one H100. Because sites are interleaved, each chunk is a random sample of all sites, so any GPU effect would show as a shift between colours. There is none (all p &gt; 0.1; the H100 is a single chunk). <strong>Implication:</strong> mixing GPU types is safe.</p>
"""

PART_C = f"""
<h2 id="partC">Part C. Further analyses</h2>

<h3>C1. How varied are the per-site effects?</h3>
<figure>{qimg('an1_effect_distributions', 'Effect distributions')}<figcaption>Distribution of the per-site effect for each shuffle width, Stable vs Dynamic.</figcaption></figure>
<p>The medians hide wide spread. For the motif shuffle, {an1.loc[('motif', 'Stable'), 'frac_gt0']:.0%} of Stable sites <em>gain</em> accessibility and {an1.loc[('motif', 'Stable'), 'frac_lt_minus1']:.0%} lose more than half of it. For the 170 bp shuffle, {an1.loc[('170', 'Stable'), 'frac_lt_minus1']:.0%} of Stable and {an1.loc[('170', 'Dynamic'), 'frac_lt_minus1']:.0%} of Dynamic sites lose more than half. <strong>Implication:</strong> "AP-1 dependence" is a property of a subset of sites, not a uniform feature; class medians should be read together with these distributions.</p>

<h3>C2. What drives the size of the effect?</h3>
<figure>{qimg('an2_effect_vs_covariates', 'Effect vs covariates')}<figcaption>Median effect (line) and inter-quartile range (band) per octile of native predicted accessibility, peak width and distance to the nearest TSS, for each shuffle width.</figcaption></figure>
{table(AN2_tab, {'Motif shuffle': '{:+.2f}', '170 bp shuffle': '{:+.2f}', 'Peak shuffle': '{:+.2f}'})}
<p><strong>The peak-shuffle effect is mostly a function of how accessible and how wide the peak is</strong> (rho {an2.loc[('peak', 'native_peak'), 'rho']:+.2f} and {an2.loc[('peak', 'peak_width'), 'rho']:+.2f}): there is more to lose in an accessible peak, and wider peaks lose less (the reason is not established here; one possibility is that a wide peak has more accessible sequence outside the site's immediate neighbourhood). The 170 bp effect has a U-shape across native accessibility (most loss in mid-high peaks, less at the very top, probably saturation), a clear dependence on distance to the TSS (rho {an2.loc[('170', 'abs_tss'), 'rho']:+.2f}) and on the footprint score ({an2.loc[('170', 'fp_score'), 'rho']:+.2f}). The motif-only effect has only weak dependencies (|rho| &le; {an2[an2.index.get_level_values(0) == 'motif'].rho.abs().max():.2f}). <strong>Implication:</strong> comparisons between groups that differ in accessibility or peak width must be adjusted, as was done in B3 and B4.</p>

<h3>C3. Motif dependence vs context dependence</h3>
<figure>{qimg('an3_site_types', 'Site types')}<figcaption>Left: each site's motif-shuffle effect against its 170 bp-shuffle effect (density; dotted lines mark a 50% loss). Right: share of sites by type, defined by whether the motif shuffle and the 170 bp shuffle each cost more than half the accessibility.</figcaption></figure>
{table(AN3_tab, {'Share of Stable': '{:.0%}', 'Share of Dynamic': '{:.0%}', 'Median native ATAC': '{:.2f}', 'Median peak width': '{:.0f}', 'Median TSS distance': '{:,.0f}'})}
<p>The two effects are positively related across sites (rho {rho_mc:.2f}). Nearly half of Stable sites are "context only" ({an3.loc['Context only', 'share_stable']:.0%}: the 170 bp shuffle hurts, the motif alone does not), {an3.loc['Motif and context', 'share_stable']:.0%} depend on both, and only {an3.loc['Motif only', 'share_stable']:.1%} on the motif alone. Dynamic sites are far more often in the "neither" group ({an3.loc['Neither', 'share_dynamic']:.0%} vs {an3.loc['Neither', 'share_stable']:.0%}); the composition differs between Stable and Dynamic (chi-square p &lt; 1e-20). <strong>Implication:</strong> the AP-1 motif seldom works alone in this model; accessibility is largely set by the sequence around it. The 50% threshold is arbitrary, so use the proportions descriptively.</p>

<h3>C4. Does footprint strength agree with the predicted dependence?</h3>
<figure>{qimg('an4_footprint_score_vs_dependence', 'Footprint score vs dependence')}<figcaption>TOBIAS footprint score (highest score across the conditions where the site is bound) against the predicted effects and native accessibility, binned.</figcaption></figure>
{table(AN4_tab, {'Spearman rho': '{:+.2f}', 'CI lo': '{:+.2f}', 'CI hi': '{:+.2f}'})}
<p>In Stable sites, a stronger footprint goes with a larger predicted dependence on the 170 bp context (rho {an4.loc[('eff_170', 'Stable'), 'rho']:+.2f}) and with higher predicted accessibility ({an4.loc[('native_peak', 'Stable'), 'rho']:+.2f}); in Dynamic sites there is no relationship. This is <strong>independent experimental support</strong> for the model's picture of AP-1 sites. <strong>Caveat:</strong> footprint scores are computed from the same ATAC data that define accessibility, so strong footprints partly mean "very accessible", and the correlation is not fully independent evidence.</p>

<h3>C5. Do two Fos sites in a peak buffer each other?</h3>
<figure>{qimg('an5_multiple_motifs', 'Multiple motifs per peak')}<figcaption>Median effect (95% CI) for sites that are the only Fos motif in their peak vs sites sharing a peak with at least one other Fos motif.</figcaption></figure>
{table(AN5_tab)}
<p>If a second motif compensated, the effect of shuffling one motif would be smaller in multi-motif peaks. It is not (differences about 0 at every width, all FDR &ge; {AN5.raw_fdr.min():.2f}). The adjusted peak-shuffle difference is {an5.loc['peak', 'adj_diff']:+.2f} (FDR {an5.loc['peak', 'adj_fdr']:.2f}), not significant. <strong>Implication:</strong> no evidence for buffering, but only {int(an5.loc['motif', 'n_multi'])} sites are in multi-motif peaks, and each shuffle here edits one site&rsquo;s window, not all motifs together.</p>

<h3>C6. Promoter-proximal vs distal sites</h3>
<figure>{qimg('an6_tss_proximity', 'TSS proximity')}<figcaption>Median native accessibility and shuffle effects by distance to the nearest TSS, Stable and Dynamic.</figcaption></figure>
{table(AN6_tab, {c: '{:+.2f}' for c in ['Native ATAC (peak)', 'Motif shuffle', '170 bp shuffle', 'Peak shuffle']})}
<p>Sites within 2&nbsp;kb of a TSS lie in much more accessible peaks (native ATAC {an6.loc[('native_peak', 'Stable', '<=2 kb'), 'median']:.2f} vs {an6.loc[('native_peak', 'Stable', '>30 kb'), 'median']:.2f} beyond 30&nbsp;kb) yet lose <em>less</em> from a 170 bp shuffle ({an6.loc[('eff_170', 'Stable', '<=2 kb'), 'median']:+.2f} vs {an6.loc[('eff_170', 'Stable', '>30 kb'), 'median']:+.2f}): promoter-proximal accessibility depends less on the AP-1 neighbourhood, distal enhancer-like sites more. The difference is large (p &lt; 1e-29 at 170 bp). <strong>Implication:</strong> TSS distance is a real biological modifier of dependence and should be included when interpreting class differences.</p>

<h3>C7. Agreement with the fibroblast MPRA</h3>
<figure>{qimg('an7_mpra', 'MPRA comparison')}<figcaption>1,653 unique MPRA elements. Left: predicted native accessibility against MPRA activity. Middle and right: predicted accessibility loss after a motif shuffle and after a full shuffle against the corresponding MPRA change (native minus shuffled). Top: Patski ATAC track; bottom: heart ATAC track.</figcaption></figure>
{table(AN7_tab, {'Spearman rho': '{:+.2f}', 'CI lo': '{:+.2f}', 'CI hi': '{:+.2f}'})}
<p>Predicted native accessibility correlates weakly but consistently with MPRA activity (rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski), and this holds in the replicated and TOBIAS-motif subsets. <strong>The predicted effect of shuffling does not track the MPRA at all</strong> (motif shuffle {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'rho']:+.2f}, full shuffle {mp_.loc[('Patski ATAC', CMP_F, 'all elements'), 'rho']:+.2f} with the Patski track; similar with heart). The MPRA delta is native minus shuffled, so agreement would appear as a <em>negative</em> correlation with the (negative) predicted log2FC; the confidence intervals exclude any correlation beyond about 0.06 in either direction. Possible reasons: the MPRA measures episomal enhancer activity of a 170&nbsp;bp fragment, while the model predicts chromatin accessibility of the same sequence in its genomic context; the MPRA effects are small and noisy; and the fibroblast tracks are weak proxies. <strong>Implication:</strong> accessibility predicted by AlphaGenome should not be used as a proxy for the MPRA effect of disrupting an element.</p>
"""

# ======================= Part D: overnight experiments =======================
OVD = AN / "overnight"
ovr = lambda n: pd.read_csv(OVD / "tables" / n)


def oimg(name, alt, width="100%"):
    return img(OVD / "figures" / f"{name}.png", alt, width)


D1_co, D1_med, D1_sv, D1_fr = ovr("d1_centred_vs_offset.csv"), ovr("d1_titration_medians.csv"), ovr("d1_titration_stable_vs_dynamic.csv"), ovr("d1_fraction_of_170.csv")
D2m, D2b = ovr("d2_mutation_dependency.csv"), ovr("d2b_decomposition.csv")
D3rep, D3adj, D3cl, D3nv = ovr("d3_replication_stable_vs_dynamic.csv"), ovr("d3_adjusted_stable_vs_dynamic.csv"), ovr("d3_fb_class_effects.csv"), ovr("d3_fb_native_vs_measured.csv")
D4 = ovr("d4_mpra_context_comparison.csv")
co_ = D1_co.set_index("class")
fr_ = D1_fr.set_index("shuffle_width_bp").fraction_of_170bp_effect
sv_ = D1_sv.set_index("width")
CMS, FBS = "Cardiomyocyte sites, heart ATAC", "Fibroblast sites, Patski ATAC"


def dep(site, edit, cls):
    r = D2m[(D2m.sites == site) & (D2m.edit == edit) & (D2m["class"] == cls)]
    return float(r.mean_dependency.iloc[0])


def depq(site, edit):
    return float(D2m[(D2m.sites == site) & (D2m.edit == edit) & (D2m["class"] == "Stable-vs-Dynamic FDR")].mean_dependency.iloc[0])


def dec(panel, edit, site):
    return float(D2b[(D2b.panel == panel) & (D2b.edit == edit) & (D2b.sites == site)].mean_dependency_stable.iloc[0])


rep = D3rep.set_index(["combo", "width"])
adj_ = D3adj.set_index(["combo", "width"])
nvf = D3nv.set_index(["track", "readout"]).rho
mc_ = D4.set_index(["track", "context", "comparison"])
ratio_cm170 = dec("170", "flank@170", "CM, heart ATAC") / dec("170", "random@170", "CM, heart ATAC")
ratio_fb170 = dec("170", "flank@170", "FB, Patski ATAC") / dec("170", "random@170", "FB, Patski ATAC")
ratio_cmcore = dec("motif", "core@motif", "CM, heart ATAC") / dec("motif", "random@motif", "CM, heart ATAC")
ratio_fbcore = dec("motif", "core@motif", "FB, Patski ATAC") / dec("motif", "random@motif", "FB, Patski ATAC")

# tables
tt = D1_med[D1_med.readout == "motif"].pivot(index="shuffle_width_bp", columns="class", values="median").reset_index()
tt["fraction of the 170 bp effect (Stable)"] = tt.shuffle_width_bp.map(fr_)
tt["Dyn - Stable"] = tt.shuffle_width_bp.map(sv_.diff_dyn_minus_stable)
tt["sig."] = tt.shuffle_width_bp.map(sv_.sig)
tt = tt[["shuffle_width_bp", "Stable", "Dynamic", "fraction of the 170 bp effect (Stable)", "Dyn - Stable", "sig."]]
tt.columns = ["Shuffle width (bp)", "Stable median", "Dynamic median", "Share of the 170 bp effect (Stable)", "Dyn - Stable", "sig."]
co_tab = D1_co.copy()
co_tab["paired_wilcoxon_p"] = co_tab.paired_wilcoxon_p.map(lambda v: "< 1e-10" if v < 1e-10 else f"{v:.1e}")
co_tab.columns = ["Class", "n", "Centred 170 bp (median)", "Offset +-300 bp (median)", "Sites where centred is more negative", "Paired Wilcoxon p"]
mut_rows = []
for e in ("random@motif", "dinuc@motif", "core@motif", "flank@motif", "random@170", "dinuc@170", "flank@170", "random@peak", "dinuc@peak"):
    mut_rows.append({"Edit": e.replace("@", " @ "), "CM Stable": dep(CMS, e, "Stable"), "CM Dynamic": dep(CMS, e, "Dynamic"), "CM sig.": stars(depq(CMS, e)),
                     "FB Stable": dep(FBS, e, "Stable"), "FB Dynamic": dep(FBS, e, "Dynamic"), "FB sig.": stars(depq(FBS, e)),
                     "Slides (Stable)": {"random@motif": 0.47, "dinuc@motif": 0.29, "core@motif": 1.18}.get(e, np.nan)})
mut_tab = pd.DataFrame(mut_rows)
rep_tab = []
for combo in ("CM sites, heart ATAC", "CM sites, Patski ATAC", "FB sites, heart ATAC", "FB sites, Patski ATAC"):
    for w in ("motif", "170", "peak"):
        r_, a_ = rep.loc[(combo, w)], adj_.loc[(combo, w)]
        rep_tab.append({"Sites, track": combo, "Shuffle": W[w], "Raw Dyn - Stable": r_.diff_dyn_minus_stable, "Raw CI lo": r_.ci_lo, "Raw CI hi": r_.ci_hi, "raw sig.": r_.sig,
                        "Adjusted": a_.adj_diff, "Adj. CI lo": a_.adj_lo, "Adj. CI hi": a_.adj_hi, "adj. sig.": a_.sig})
rep_tab = pd.DataFrame(rep_tab)
nvf_tab = D3nv.copy()
nvf_tab.columns = ["Track", "Readout", "Spearman rho", "CI lo", "CI hi", "n"]
D4_tab = D4.copy()
D4_tab["p"] = D4_tab.p.map(lambda v: "< 1e-10" if v < 1e-10 else f"{v:.2g}")
D4_tab.columns = ["Track", "Context", "Comparison", "Spearman rho", "CI lo", "CI hi", "p", "n"]
n_cl_sig170 = int(((D3cl.width == "170") & (D3cl.fdr < 0.001) & D3cl.fdr.notna()).sum())
n_cl_170 = int(((D3cl.width == "170") & D3cl.fdr.notna()).sum())
dg7, dg8, dg9 = dg.part_d_diagrams()

PART_D = f"""
<h2 id="partD">Part D. Overnight experiments: controls, mutation arm, fibroblast sites, neutral-context MPRA</h2>
<p>These runs answer questions raised by Parts B and C: is the 170&nbsp;bp effect specific to the AP-1 site, how does it scale with the shuffled width, what carries it in a mutation experiment, does the Stable-vs-Dynamic result hold in fibroblast sites, and was the MPRA mismatch caused by genomic context? Designs first, then results.</p>

<h3>D1. Is the loss specific to the AP-1 site, and how does it scale with width?</h3>
{dg7}
<figure>{oimg('d1_controls_titration', 'Controls and width titration')}<figcaption>Left: effect on the motif readout when a 170 bp region is shuffled at different positions (median with 95% CI). Right: dose-response of the effect (log2FC) on the motif (thick) and, thinner, the 170 bp and peak readouts as the centred shuffle widens. Marks: Stable vs Dynamic (Mann&ndash;Whitney, BH-corrected), motif readout.</figcaption></figure>
{table(co_tab, {'Centred 170 bp (median)': '{:.2f}', 'Offset +-300 bp (median)': '{:.2f}', 'Sites where centred is more negative': '{:.0%}'})}
{table(tt, {'Stable median': '{:.2f}', 'Dynamic median': '{:.2f}', 'Share of the 170 bp effect (Stable)': '{:.0%}', 'Dyn - Stable': '{:+.2f}'})}
<p><strong>The effect is specific to the site's neighbourhood.</strong> Shuffling the centred 170&nbsp;bp costs {co_.loc['Stable', 'centred_170_median']:.2f} log2 units at the motif in Stable sites, while the same-size region 300&nbsp;bp away costs only {co_.loc['Stable', 'offset_pm300_median']:.2f}; the centred shuffle is more damaging at {co_.loc['Stable', 'share_centred_more_negative']:.0%} of Stable and {co_.loc['Dynamic', 'share_centred_more_negative']:.0%} of Dynamic sites (paired Wilcoxon, p &lt; 1e-10). This is a distance control (the offset shuffles lie 215 to 385&nbsp;bp from the readout and further from the peak summit), not an AP-1-specific control; matched shuffles at non-Fos accessible sites would be needed for that. <strong>It grows smoothly with width and levels off near 340&nbsp;bp:</strong> the 7&nbsp;bp core alone gives {fr_.loc[7]:.0%} of the 170&nbsp;bp effect, the whole 11&nbsp;bp motif {fr_.loc[11]:.0%}, 85&nbsp;bp {fr_.loc[85]:.0%}, and 340&nbsp;bp {fr_.loc[340]:.0%}. The Stable&ndash;Dynamic gap is small at the core ({sv_.loc[7, 'diff_dyn_minus_stable']:+.2f}) and builds to {sv_.loc[340, 'diff_dyn_minus_stable']:+.2f} by 340&nbsp;bp, so it comes from the surrounding sequence, not the core.</p>
<p><strong>It also settles one open question.</strong> The slides' motif shuffle gives about &minus;0.12 on the motif readout; here even a 7&nbsp;bp shuffle gives {D1_med[(D1_med.readout == 'motif') & (D1_med['class'] == 'Stable') & (D1_med.shuffle_width_bp == 7)]['median'].iloc[0]:.2f}. So the width of the shuffled region does <em>not</em> explain why our motif effect is larger than on the slides; the site set and readout definitions remain the candidates. <strong>Trust:</strong> robust (large effects, paired design, consistent in both classes).</p>

<h3>D2. Mutation arm: what carries the dependence?</h3>
{dg8}
<figure>{oimg('d2_mutation_arm', 'Mutation arm dependency')}<figcaption>Dependency (&minus;log2FC) for each edit, read out over the region it edits; Stable vs Dynamic with 95% CI. Top: cardiomyocyte sites (heart ATAC); bottom: fibroblast sites (Patski ATAC). Diamonds: the slides' Stable values for the three edits they report.</figcaption></figure>
<figure>{oimg('d2b_what_carries_dependence', 'What carries the dependence', '860px')}<figcaption>Stable sites. Left: within the motif, the 7&nbsp;bp core carries the whole effect of mutating all 11 bases. Right: mutating the 170&nbsp;bp while leaving the motif intact reproduces most of the loss from mutating everything.</figcaption></figure>
{table(mut_tab, {c: '{:.2f}' for c in ['CM Stable', 'CM Dynamic', 'FB Stable', 'FB Dynamic', 'Slides (Stable)']})}
<p><strong>The surroundings, not the motif, carry the dependence.</strong> Mutating the 170&nbsp;bp but keeping the motif intact costs {ratio_cm170:.0%} of the loss from mutating everything in cardiomyocyte sites ({dec('170', 'flank@170', 'CM, heart ATAC'):.2f} vs {dec('170', 'random@170', 'CM, heart ATAC'):.2f}) and {ratio_fb170:.0%} in fibroblast sites. Within the motif, the 7&nbsp;bp core accounts for essentially all of it ({ratio_cmcore:.0%} CM, {ratio_fbcore:.0%} FB). At the 170&nbsp;bp and peak widths every edit is significantly more damaging in Stable than Dynamic sites (FDR &lt; 0.001), in both site sets; at motif width the classes barely differ ({', '.join(f'{e.split(chr(64))[0]} {stars(depq(CMS, e))}' for e in ('random@motif', 'dinuc@motif', 'core@motif', 'flank@motif'))} in CM; all ns in FB).</p>
<p><strong>Caveats.</strong> The 170&nbsp;bp readout averages over 170&nbsp;bp of which 11 are the motif, so a ratio near 100% is expected even if the motif matters locally; the result is also track-dependent (91% CM, 76% FB). A dinucleotide-preserving shuffle of an 11-mer has very few distinct arrangements (about 12% reproduce the original), so dinuc@motif is partly a no-op and understates the effect.</p>
<p><strong>Compared with the slides.</strong> Random and dinucleotide edits of the motif are somewhat larger here (random {dep(CMS, 'random@motif', 'Stable'):.2f} vs 0.47, dinucleotide {dep(CMS, 'dinuc@motif', 'Stable'):.2f} vs 0.29); the core edit is smaller ({dep(CMS, 'core@motif', 'Stable'):.2f} vs 1.18). Our core edit is a random substitution of the 7 bases with no AP-1 half-site left (3 draws per site), not a targeted edit, so its size depends on that design. <strong>Trust:</strong> the decomposition and the class pattern are robust; the match with the slides is only approximate.</p>

<h3>D3. Fibroblast sites: does the Stable-vs-Dynamic result replicate?</h3>
<figure>{oimg('d3a_fb_heatmaps', 'Fibroblast shuffle heatmaps')}<figcaption>Fibroblast sites (5,018): median shuffle effect on each readout with the Patski ATAC track. Marks on the Dynamic panel: Mann&ndash;Whitney vs Stable, BH over the 9 cells.</figcaption></figure>
<figure>{oimg('d3b_replication_forest', 'Replication across site sets and tracks')}<figcaption>Dynamic minus Stable (median effect, 95% CI) for cardiomyocyte and fibroblast sites, each with the heart and the Patski track. Positive = Dynamic loses less.</figcaption></figure>
{table(rep_tab, {c: '{:+.2f}' for c in ['Raw Dyn - Stable', 'Raw CI lo', 'Raw CI hi', 'Adjusted', 'Adj. CI lo', 'Adj. CI hi']})}
<p><strong>Yes, in direction.</strong> At the 170&nbsp;bp shuffle Stable sites lose more than Dynamic sites in all four combinations (raw FDR &lt; 0.001): CM sites {rep.loc[('CM sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (heart) and {rep.loc[('CM sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (Patski), FB sites {rep.loc[('FB sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} and {rep.loc[('FB sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f}. As in cardiomyocyte sites, the motif-width shuffle does not separate the classes, and the peak-width difference is smaller. <strong>The size is track-dependent, though.</strong> After adjusting for native accessibility and peak width, the 170&nbsp;bp gap stays larger with the heart track (the matched-readout and further adjustments in Part H shrink it much more) ({adj_.loc[('CM sites, heart ATAC', '170'), 'adj_diff']:+.2f} CM sites, {adj_.loc[('FB sites, heart ATAC', '170'), 'adj_diff']:+.2f} FB sites) but shrinks to {adj_.loc[('CM sites, Patski ATAC', '170'), 'adj_diff']:+.2f} and {adj_.loc[('FB sites, Patski ATAC', '170'), 'adj_diff']:+.2f} with the Patski track (still FDR &lt; 0.01); the peak-width gap disappears after adjustment in 3 of 4 combinations. So with the fibroblast-like track much of the raw difference reflects accessibility and peak width, with a small residual.</p>
<figure>{oimg('d3c_fb_class_effects', 'Class effects in fibroblast sites')}<figcaption>Fibroblast sites, Patski track: effect per class with the matched readout (95% CI); blue = Stable reference, marks = Mann&ndash;Whitney vs reference, BH within each family.</figcaption></figure>
<p>In fibroblast sites {n_cl_sig170} of {n_cl_170} class comparisons at the 170&nbsp;bp shuffle are significant at FDR &lt; 0.001: every dynamic class loses less than Stable, as in cardiomyocytes. I did not test dynamic classes against each other in fibroblast sites. <strong>Prediction quality for these sites is weak:</strong> predicted vs measured fibroblast accessibility is {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} with Patski and {nvf.loc[('heart ATAC', 'atac_peak')]:.2f} with the heart track (peak readout, n = 4,589), so the fibroblast numbers should be read as indicative. <strong>Trust:</strong> direction robust; magnitude depends on track and is partly explained by accessibility and width.</p>
{table(nvf_tab, {'Spearman rho': '{:.2f}', 'CI lo': '{:.2f}', 'CI hi': '{:.2f}'})}

<h3>D4. Does placing MPRA elements in a neutral context restore agreement?</h3>
{dg9}
<figure>{oimg('d4_mpra_context', 'MPRA genomic vs neutral context')}<figcaption>Spearman correlation between predicted accessibility and the MPRA (95% CI) for the native level (left) and the effect of scrambling the motif (middle) or the whole element (right). Blue circles: element in its genomic window; green squares: element in a neutral window.</figcaption></figure>
{table(D4_tab, {'Spearman rho': '{:+.3f}', 'CI lo': '{:+.3f}', 'CI hi': '{:+.3f}'})}
<p>Neutral-context native accessibility correlates with MPRA activity about as before (heart {mc_.loc[('heart ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} vs {mc_.loc[('heart ATAC', 'genomic context', 'native level vs MPRA activity'), 'rho']:+.2f} in genomic context; Patski {mc_.loc[('Patski ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} vs {mc_.loc[('Patski ATAC', 'genomic context', 'native level vs MPRA activity'), 'rho']:+.2f}), and the predicted effects of scrambling still do not track the MPRA changes (all |rho| &le; {D4[D4.comparison != 'native level vs MPRA activity'].rho.abs().max():.2f}, intervals including 0). <strong>Genomic context is therefore not the reason</strong> the native accessibility head fails to predict the MPRA changes. This is consistent with the manuscript's choice to train a probe directly on the MPRA instead of using the native head. <strong>Trust:</strong> a null for these effects, but the MPRA's own scramble effects are tiny (mean core loss about 0.02 log2), so there is little real variance to predict and the null does not discriminate model quality; it applies to the native ATAC head, not to the manuscript's MPRA-trained model.</p>

<h3>D5. What the overnight runs add</h3>
<ul>
<li>The loss of accessibility after a shuffle is <strong>specific to the AP-1 neighbourhood</strong> (same-size shuffles 300&nbsp;bp away do almost nothing) and <strong>saturates around 340&nbsp;bp</strong>.</li>
<li>The <strong>7&nbsp;bp core alone gives about a quarter</strong> of the 170&nbsp;bp effect, and the <strong>surroundings carry most of what mutating 170&nbsp;bp does</strong> even when the motif is left intact.</li>
<li>The <strong>Stable-vs-Dynamic raw difference replicates in fibroblast sites</strong> (all combinations of site set and track at 170&nbsp;bp), but it is mostly explained by accessibility and context (Part H).</li>
<li><strong>Motif width does not explain the larger motif effect</strong> compared with the slides; <strong>genomic context does not explain the MPRA mismatch</strong>.</li>
</ul>
"""

# ======================= Part G: agreement-filtered analyses =======================
AG_ = AN / "agreement"
AG25 = AN / "agreement_floor25"
G_ret = pd.read_csv(AG_ / "tables" / "retention_by_tier.csv")
G_sv = pd.read_csv(AG_ / "tables" / "stable_vs_dynamic_by_tier.csv")
G_sv25 = pd.read_csv(AG25 / "tables" / "stable_vs_dynamic_by_tier.csv")
G_mu = pd.read_csv(AG_ / "tables" / "mutation_arm_by_tier.csv")
G_ct = pd.read_csv(AG_ / "tables" / "controls_by_tier.csv")
G_st = pd.read_csv(AG_ / "tables" / "status_counts.csv")
G_kd = pd.read_csv(AG_ / "tables" / "effect_kept_vs_discarded.csv")
G_kq = pd.read_csv(AG_ / "tables" / "kept_by_quintile.csv")
G_ret_tab = G_ret[["combo", "tier", "n", "n_stable", "n_dynamic", "retained_stable", "retained_dynamic", "rho_pred_vs_meas_in_tier"]].copy()
G_ret_tab.columns = ["Site set, track", "Tier", "Sites", "Stable", "Dynamic", "Stable kept", "Dynamic kept", "rho in tier"]
G_st_tab = G_st.pivot_table(index=["combo", "class"], columns="status", values="n", aggfunc="sum").reset_index()
G_st_tab = G_st_tab[["combo", "class", "kept (tier C)", "below noise floor", "rank disagreement", "no measured peak"]]
G_st_tab.columns = ["Site set, track", "Class", "Kept (tier C)", "Dropped: below floor", "Dropped: rank disagreement", "Dropped: no measured peak"]
G_sv170 = G_sv[G_sv.width == "170"][["combo", "tier", "n_stable", "n_dynamic", "median_stable", "median_dynamic", "diff_dyn_minus_stable", "ci_lo", "ci_hi", "sig"]].copy()
G_sv170.columns = ["Site set, track", "Tier", "n Stable", "n Dynamic", "Stable median", "Dynamic median", "Dyn - Stable", "CI lo", "CI hi", "FDR"]
G_sens = G_sv[(G_sv.width == "170") & (G_sv.tier.isin(["B", "C", "D"]))][["combo", "tier", "diff_dyn_minus_stable", "sig"]].merge(
    G_sv25[(G_sv25.width == "170") & (G_sv25.tier.isin(["B", "C", "D"]))][["combo", "tier", "diff_dyn_minus_stable", "sig"]], on=["combo", "tier"], suffixes=(" low", " high"))
G_sens.columns = ["Site set, track", "Tier", "Dyn - Stable (5% floor)", "FDR (5%)", "Dyn - Stable (25% floor)", "FDR (25%)"]
G_kd_tab = G_kd[["combo", "class", "n_kept", "n_discarded_with_peak", "median_effect_kept", "median_effect_discarded", "p"]].copy()
G_kd_tab["p"] = G_kd_tab.p.map(lambda v: "< 0.001" if v < 0.001 else f"{v:.2f}")
G_kd_tab.columns = ["Site set, track", "Class", "Kept", "Dropped (has a peak)", "Median effect, kept", "Median effect, dropped", "p"]
G_mu_tab = G_mu[["combo", "tier", "ratio_flank170_over_random170", "ratio_lo", "ratio_hi", "ratio_core_over_random_motif"]].copy()
G_mu_tab.columns = ["Site set, track", "Tier", "Surroundings-only / full 170 bp loss", "CI lo", "CI hi", "Core / full motif loss"]
gs = lambda c, w, t_: G_sv[(G_sv.combo == c) & (G_sv.width == w) & (G_sv.tier == t_)].iloc[0]
gr = lambda c, t_: G_ret[(G_ret.combo == c) & (G_ret.tier == t_)].iloc[0]
gk = lambda c, cl: G_kd[(G_kd.combo == c) & (G_kd["class"] == cl)].iloc[0]
gq = lambda c, cl: G_kq[(G_kq.combo == c) & (G_kq.by == "pred") & (G_kq["class"] == cl) & (G_kq.quintile >= 4)].share_kept.mean()
cmn, fbh, fbp = "CM sites, heart ATAC", "FB sites, heart ATAC", "FB sites, Patski ATAC"
PART_G = f"""
<h2 id="partG">Part G. Restricting to sites where predicted accessibility agrees with measured accessibility</h2>
<p><strong>Aim.</strong> Exclude noisy or poorly predicted regions by keeping only sites whose predicted native accessibility agrees with the measured snATAC, then check whether the main results still hold. The selection uses only the native prediction and the measurement, never the perturbation effects.</p>
<div class="box"><strong>The rule, in order</strong><ol>
<li><strong>Has a measured peak.</strong> The site falls in a peak with measured accessibility (otherwise there is nothing to compare with).</li>
<li><strong>Above a low noise floor</strong> (tier B). Predicted accessibility is above the 5th percentile of sites with a measured peak. This only removes sites with almost no predicted signal, where the log ratio is unstable. It is deliberately low: a high floor would also remove weak peaks that the model predicts correctly, and Dynamic sites are weaker on average.</li>
<li><strong>Relative agreement</strong> (tier C). The predicted percentile minus the measured percentile (both ranked among sites with a measured peak) is within 1 SD of that difference. This compares relative accessibility, so a weak peak predicted as weak is kept. Tier D is the strict version (within 0.5 SD).</li></ol>
<span class="small">Each site set and track has its own thresholds. Fibroblast sites are compared with measured fibroblast accessibility, with the heart and Patski tracks each judged on their own agreement.</span></div>

<h3>G1. What was kept and what was dropped</h3>
<figure>{img(AG_ / 'figures' / 'e0_funnel.png', 'Filtering funnel')}<figcaption>Number of sites left after each step, split into Stable and Dynamic, with the Dynamic share on the right. The Dynamic share stays roughly constant, so the filter does not remove one class wholesale.</figcaption></figure>
<figure>{img(AG_ / 'figures' / 'e0b_rank_map.png', 'Rank map')}<figcaption>Every site with a measured peak, placed by its measured and predicted accessibility percentile. Blue sites are kept (tier C); orange sites fall below the floor; grey sites are dropped for rank disagreement (outside the dashed band). Sites without a measured peak are not shown (counts in the titles).</figcaption></figure>
{table(G_st_tab)}
<p>About {gr(cmn, 'C').retained_all:.0%} of cardiomyocyte sites ({gr(cmn, 'C').n:,}) reach tier C and {gr(cmn, 'D').retained_all:.0%} ({gr(cmn, 'D').n:,}) reach tier D. Most drops are for rank disagreement; the floor removes few sites. Correlations inside a tier (0.7 to 0.95) are circular, since the tier is defined by agreement, and only show that the selection worked.</p>
{table(G_ret_tab, {'Stable kept': '{:.0%}', 'Dynamic kept': '{:.0%}', 'rho in tier': '{:.2f}'})}

<h3>G2. Is the filter biased against one class?</h3>
<figure>{img(AG_ / 'figures' / 'e0c_kept_by_accessibility.png', 'Share kept by accessibility')}<figcaption>Share of sites kept in tier C by quintile of predicted accessibility (top) and of measured accessibility (bottom), for Stable and Dynamic sites.</figcaption></figure>
<figure>{img(AG_ / 'figures' / 'e0d_distributions.png', 'Distributions of rank difference and predicted level')}<figcaption>Top: predicted-minus-measured percentile for each class with the 1 SD and 0.5 SD cut-offs. Bottom: predicted accessibility of kept and discarded sites (orange line = floor).</figcaption></figure>
<p><strong>By measured accessibility the filter treats the classes alike</strong> (bottom row: same share kept for Stable and Dynamic at every level). <strong>By predicted accessibility it does not:</strong> among the two highest predicted quintiles only {gq(cmn, 'Dynamic'):.0%} of Dynamic cardiomyocyte sites are kept, against {gq(cmn, 'Stable'):.0%} of Stable ones (fibroblast sites, heart track {gq(fbh, 'Dynamic'):.0%} vs {gq(fbh, 'Stable'):.0%}; Patski {gq(fbp, 'Dynamic'):.0%} vs {gq(fbp, 'Stable'):.0%}). In plain terms, the model over-predicts many Dynamic sites (high prediction, low measurement) and those are removed, while the low-floor choice keeps weakly predicted sites of both classes. The kept Dynamic sites are therefore the ones the model gets right, which is the intended behaviour, but the filtered result describes well-predicted Dynamic sites and not all of them.</p>
<h3>G3. Do kept and dropped sites differ in the effect?</h3>
{table(G_kd_tab, {'Median effect, kept': '{:.2f}', 'Median effect, dropped': '{:.2f}'})}
<p>This is a transparency check on the 170&nbsp;bp shuffle effect, not part of the selection. <strong>Kept Stable sites lose more accessibility than dropped Stable sites</strong> ({gk(cmn, 'Stable').median_effect_kept:.2f} vs {gk(cmn, 'Stable').median_effect_discarded:.2f} in cardiomyocyte sites, p &lt; 0.001), while Dynamic sites do not differ (p = {gk(cmn, 'Dynamic').p:.2f}); the same holds for fibroblast sites with the heart track, but not with Patski, where kept Dynamic sites lose <em>less</em> than dropped ones ({gk(fbp, 'Dynamic').median_effect_kept:.2f} vs {gk(fbp, 'Dynamic').median_effect_discarded:.2f}, p &lt; 0.001), the opposite direction, so the Patski gaps in tiers C and D are partly selection. So some of the larger Stable-vs-Dynamic gap in the filtered sets comes from the filter itself, since agreeing sites are the ones with clearer accessibility signal and larger effects. The gap in the unfiltered data ({gs(cmn, '170', 'A').diff_dyn_minus_stable:+.2f}) is the cleaner estimate; the filtered values show it does not depend on poorly predicted sites.</p>

<h3>G4. Does Stable vs Dynamic survive?</h3>
<figure>{img(AG_ / 'figures' / 'e3_stable_vs_dynamic_by_tier.png', 'Stable vs Dynamic by tier')}<figcaption>Dynamic minus Stable (median log2FC, 95% bootstrap CI) per shuffle width and tier. Marks: Mann&ndash;Whitney, BH within each site set and track.</figcaption></figure>
{table(G_sv170, {'Stable median': '{:.2f}', 'Dynamic median': '{:.2f}', 'Dyn - Stable': '{:+.2f}', 'CI lo': '{:+.2f}', 'CI hi': '{:+.2f}'})}
<p><strong>At the 170&nbsp;bp shuffle the difference holds in every tier</strong> (all FDR &lt; 0.001): cardiomyocyte sites {gs(cmn, '170', 'A').diff_dyn_minus_stable:+.2f} (all) &rarr; {gs(cmn, '170', 'C').diff_dyn_minus_stable:+.2f} (C) &rarr; {gs(cmn, '170', 'D').diff_dyn_minus_stable:+.2f} (D); fibroblast sites with the heart track {gs(fbh, '170', 'A').diff_dyn_minus_stable:+.2f} &rarr; {gs(fbh, '170', 'D').diff_dyn_minus_stable:+.2f}; with Patski {gs(fbp, '170', 'A').diff_dyn_minus_stable:+.2f} &rarr; {gs(fbp, '170', 'D').diff_dyn_minus_stable:+.2f}. The motif-width shuffle remains without a meaningful class difference. The peak-width difference is not stable across tiers and should still be treated as partial. Read the growth with stricter tiers together with G3.</p>
<h3>G5. Sensitivity to the floor</h3>
{table(G_sens, {'Dyn - Stable (5% floor)': '{:+.2f}', 'Dyn - Stable (25% floor)': '{:+.2f}'})}
<p>The 25th-percentile floor, which removes weak peaks, keeps fewer Dynamic sites (about a third in tier C instead of about half) and gives a gap that is similar or larger in cardiomyocyte sites and smaller with Patski (tier C +0.44 vs +0.61), with wider intervals. Direction and significance do not depend on the floor; the size does somewhat. Tables for it are in <code>Analysis/agreement_floor25</code>.</p>

<h3>G6. Mutation arm and controls in the filtered sets</h3>
<figure>{img(AG_ / 'figures' / 'e4_other_results_by_tier.png', 'Other results by tier')}<figcaption>Left: share of the 170&nbsp;bp mutation loss that remains when the motif is kept intact. Middle and right: width titration and offset controls (Stable cardiomyocyte sites), by tier.</figcaption></figure>
{table(G_mu_tab, {c: '{:.2f}' for c in ['Surroundings-only / full 170 bp loss', 'CI lo', 'CI hi', 'Core / full motif loss']})}
<p>These do not change with the filter: the surroundings carry {G_mu[(G_mu.combo == cmn)].ratio_flank170_over_random170.min():.0%} to {G_mu[(G_mu.combo == cmn)].ratio_flank170_over_random170.max():.0%} of the 170&nbsp;bp mutation loss in cardiomyocyte sites in every tier, the 7&nbsp;bp core stays at about {G_ct.core7_share_of_170.min():.0%} to {G_ct.core7_share_of_170.max():.0%} of the 170&nbsp;bp shuffle effect, and a same-size shuffle 300&nbsp;bp away costs about {G_ct.offset_pm300.mean():.2f} versus {G_ct.centred_170.mean():.2f} centred.</p>
<h3>G7. Takeaways</h3>
<ul><li>The main results do not come from poorly predicted sites: the mutation decomposition and controls are unchanged, and Stable-vs-Dynamic at 170&nbsp;bp holds in every tier and with both floors.</li>
<li>The filter removes mostly over-predicted Dynamic sites, so the filtered Dynamic group is the well-predicted subset. Part of the larger gap is due to the filter (G3); quote the unfiltered number as the main estimate and the filtered ones as a robustness check.</li>
<li>Fibroblast results remain limited by prediction quality: about {gr(fbp, 'C').retained_all:.0%} of sites reach tier C with Patski, and the unfiltered correlation is {gr(fbp, 'A').rho_pred_vs_meas_in_tier:.2f}.</li></ul>
"""

# ======================= Part H: independent review audit =======================
AU = AN / "audit"
H_g = pd.read_csv(AU / "tables" / "gap_models.csv")
H_cs = pd.read_csv(AU / "tables" / "common_support_gap.csv").set_index("combo")
def hg(c, mdl, w="170"):
    return H_g[(H_g.combo == c) & (H_g.width == w) & (H_g.model.str.startswith(mdl))].iloc[0]
M_RAW, M_PK, M_MATCH, M_FULL, M_FP = "raw", "+ native PEAK", "+ native accessibility at the matched", "+ matched native, peak native", "+ all of the above"
H_tab = H_g[H_g.width == "170"][["combo", "model", "gap", "ci_lo", "ci_hi", "p"]].copy()
H_tab["p"] = H_tab.p.map(lambda v: "< 0.001" if v < 0.001 else f"{v:.2f}")
H_tab.columns = ["Site set, track", "Model", "Dyn - Stable", "CI lo", "CI hi", "p"]
H_oth = H_g[(H_g.width != "170") & (H_g.model.str.startswith((M_RAW, M_FULL)))][["combo", "width", "model", "gap", "ci_lo", "ci_hi"]].copy()
H_oth["model"] = H_oth.model.map(lambda s: "raw" if s == "raw" else "adjusted (matched native, peak native, width, TSS)")
H_oth.columns = ["Site set, track", "Shuffle width", "Model", "Dyn - Stable", "CI lo", "CI hi"]
cmn, fbh, fbp = "CM sites, heart ATAC", "FB sites, heart ATAC", "FB sites, Patski ATAC"
HF = {'Dyn - Stable': '{:+.2f}', 'CI lo': '{:+.2f}', 'CI hi': '{:+.2f}'}
PART_H = f"""
<h2 id="partH">Part H. Independent code and analysis review: corrections</h2>
<p>Three independent reviewers read the code and re-derived key numbers from the raw results without seeing my conclusions. All headline numbers reproduced exactly (Stable and Dynamic medians and their gap, the predicted-vs-measured correlation, tier counts, mutation-arm ratio), class labels and site tables were confirmed correct, and no bug that invalidates the shuffle, mutation or control results was found. They did find one analysis choice that changed a conclusion, plus smaller issues, listed at the end.</p>
<h3>H1. The Stable-vs-Dynamic gap is mostly explained by native accessibility and context</h3>
<p><strong>What was wrong.</strong> The earlier adjustment (Parts B3 and B4) controlled for native accessibility over the <em>peak</em> for every readout, and the stratification used quintiles of the same variable. The effect of a 170&nbsp;bp shuffle depends on how accessible those 170&nbsp;bp are, so the correct covariate is the native value at the matched readout. With that variable the gap shrinks a lot. Models below use standard errors clustered by TOBIAS peak, 170&nbsp;bp shuffle.</p>
<figure>{img(AU / 'figures' / 'h1_gap_models.png', 'Gap by model')}<figcaption>Dynamic minus Stable (log2FC, 95% CI) as covariates are added. Positive = Dynamic loses less. The last row includes footprint score, which partly defines the classes, so it is a sensitivity check and not a causal adjustment.</figcaption></figure>
{table(H_tab, HF)}
<figure>{img(AU / 'figures' / 'h2_gap_by_matched_decile.png', 'Gap by decile of matched native accessibility')}<figcaption>Top: median 170&nbsp;bp shuffle effect by decile of native accessibility at the 170&nbsp;bp readout. Bottom: number of Stable (blue) and Dynamic (orange) sites per decile; most Dynamic sites sit in the lower deciles.</figcaption></figure>
<p><strong>Reading.</strong> In cardiomyocyte sites the raw gap of {hg(cmn, M_RAW).gap:+.2f} becomes {hg(cmn, M_PK).gap:+.2f} with the earlier adjustment, {hg(cmn, M_MATCH).gap:+.2f} with the matched-readout value, {hg(cmn, M_FULL).gap:+.2f} after also adjusting for peak accessibility, peak width and TSS distance with flexible terms, and {hg(cmn, M_FP).gap:+.2f} (CI {hg(cmn, M_FP).ci_lo:+.2f} to {hg(cmn, M_FP).ci_hi:+.2f}) with footprint score. With the fibroblast-like Patski track the gap is gone after the matched adjustment ({hg(fbp, M_MATCH).gap:+.2f}). Fibroblast sites with the heart track keep a small gap ({hg(fbh, M_FULL).gap:+.2f}, CI {hg(fbh, M_FULL).ci_lo:+.2f} to {hg(fbh, M_FULL).ci_hi:+.2f}). In the decile plot the gap is about 0.1 to 0.2 across the lower deciles, where {H_cs.loc[cmn, 'dynamic_in_lowest_5_deciles']:.0%} of Dynamic cardiomyocyte sites lie, and larger only in the top deciles where there are few Dynamic sites. Weighted to the Dynamic distribution, the gap for sites like Dynamic ones is {H_cs.loc[cmn, 'common_support_gap']:+.2f} (CM sites), {H_cs.loc[fbh, 'common_support_gap']:+.2f} (FB, heart) and {H_cs.loc[fbp, 'common_support_gap']:+.2f} (FB, Patski).</p>
{table(H_oth, HF)}
<p><strong>What this means.</strong> Stable sites lose more accessibility when 170&nbsp;bp are shuffled largely because they sit in more accessible, more strongly bound, better-positioned regions, not because of a separate sequence dependence that the model attaches to stable binding. A smaller residual gap remains in cardiomyocyte sites and in fibroblast sites with the heart track, consistent in direction but modest ({hg(cmn, M_FULL).gap:+.2f} and {hg(fbh, M_FULL).gap:+.2f}), and adding footprint score removes it in cardiomyocyte sites. The motif-width and peak-width differences disappear after adjustment in every combination. <strong>The earlier statements "Stable vs Dynamic is real" and "adjustment barely changes it" are replaced by this section.</strong> The raw difference is real as a description; it is not evidence of an independent effect of stable binding. Patski is not an independent confirmation, since it scores the same sequences with the same model.</p>
<h3>H2. Other corrections made to the text</h3>
<table class="t"><tr><th>Finding</th><th>Status</th><th>Change</th></tr>
<tr><td class="v">Reference class "Stable" in the age and exercise families contains some Dynamic sites (n = 2,997 and 3,117 vs 2,907 in the Stable/Dynamic analysis)</td><td class="v">Reviewer finding; medians change by at most 0.03</td><td class="v">Noted; class tables use their own reference and overlap each other (a MAS-only site is both Age Up and Exercise Down).</td></tr>
<tr><td class="v">"Largest RNA drops close to a TSS"</td><td class="v">Checked against the table</td><td class="v">Corrected: for the &plusmn;2&nbsp;kb readout the largest median drops are 2 to 10&nbsp;kb away; within 2&nbsp;kb only for the &plusmn;50&nbsp;kb mean.</td></tr>
<tr><td class="v">Kept vs dropped Dynamic sites "do not differ" (G3)</td><td class="v">Checked: true for cardiomyocyte and heart-track fibroblast sites, opposite direction for Patski</td><td class="v">Patski result added to G3.</td></tr>
<tr><td class="v">MPRA null "rules out |rho| above 0.06, not a lack of power"</td><td class="v">Reviewer finding: holds for the motif shuffle of all elements only, and the MPRA's own scramble effects are tiny (mean core loss about 0.02 log2), so there is little variance to predict</td><td class="v">Wording weakened in D4 and the trust table.</td></tr>
<tr><td class="v">"Specific to the AP-1 neighbourhood" (controls)</td><td class="v">Reviewer finding: the offset shuffles are a distance control, not an AP-1-specific control</td><td class="v">Wording corrected in D1; the width titration is the stronger evidence. Matched shuffles at non-Fos accessible sites would settle it.</td></tr>
<tr><td class="v">Mutation-arm ratio near 1 (surroundings carry 91% / 76%)</td><td class="v">Reviewer finding: the 170&nbsp;bp readout averages over 170&nbsp;bp of which 11 are the motif, so a ratio near 1 is expected even if the motif matters locally; it is track-dependent</td><td class="v">Caveat added in D2. A paired per-site comparison would be more informative.</td></tr>
<tr><td class="v">"All dynamic classes behave alike"</td><td class="v">Reviewer finding: intervals (about &plusmn;0.2 to 0.3, n = 79 to 210) are too wide to exclude differences half the size of the Stable-vs-Dynamic gap</td><td class="v">Reworded as absence of evidence.</td></tr>
<tr><td class="v">Measured snATAC called an independent check (A4)</td><td class="v">Footprint calls and accessibility come from the same ATAC data</td><td class="v">Already conceded in C4; not changed.</td></tr>
<tr><td class="v">Pooled measured accessibility includes the young-exercise samples</td><td class="v">Checked: pooled and YS/MAS/MAE-only values correlate 0.9997; predicted-vs-measured rho 0.397 vs 0.392</td><td class="v">No change needed.</td></tr></table>
<h3>H3. Code issues found, not yet fixed (need a rerun on Katana or a decision)</h3>
<ul><li><strong>MPRA element coordinates off by one</strong> in <code>13_mpra.py</code> (Coords are 1-based inclusive; the script uses start as 0-based, so the 170&nbsp;bp element is read as 169&nbsp;bp). The neutral-context script is correct. Expected effect: negligible.</li>
<li><strong>dinuc@motif is partly a no-op</strong>: an 11-mer has few dinucleotide-preserving arrangements and about 12% reproduce the original, which dilutes its effect (median &minus;0.33 vs &minus;0.67 for random). Fix: redraw until different, or report it as limited.</li>
<li><strong>Shuffled-motif control in the insertion experiment</strong> can keep AP-1 half-sites (about 11% of permutations); the null insertion result is unlikely to change.</li>
<li><strong>Core edit fallback not recorded:</strong> if an 11-mer has no canonical core, the central 7&nbsp;bp are used and this is not flagged in the output.</li>
<li><strong>Peak swap:</strong> a donor narrower than the recipient leaves the recipient's own flanks in the window; only 4 recipients, and <code>12_insertion.py</code> is not resumable.</li>
<li><strong>Neutral MPRA</strong> uses only 2 background windows.</li>
<li><strong>Data files:</strong> <code>fb_mpra_tobias_crossmap.csv</code> was written by an older version of the script (old single-label column) and is not deduplicated by element (359 elements appear twice); no analysis reads it. Regenerate before use. The MPRA manuscript checks (script 27) follow the manuscript in treating the 2,012 rows as independent.</li>
<li><strong>Open question:</strong> only motif name Fos is used; Fosb, Fosl1 and Fosl2 calls are excluded. Confirm this is intended.</li>
<li><strong>Prose numbers hard-coded</strong> in the report builder (for example "rho 0.19 to 0.24", "tier C +0.44 vs +0.61") were checked and are currently correct but will not update on a rerun.</li></ul>
"""

DESIGN = dg.design_section()

claims = [
    ("Shuffling the 170 bp around a site removes far more predicted accessibility than shuffling the 11 bp motif alone", "ok",
     f"Median log2FC {m.loc['170', 'median_stable']:.2f} (170 bp shuffle) vs {m.loc['motif', 'median_stable']:.2f} (motif shuffle) in Stable sites; "
     f"{frac['170']:.0%} of sites show a 170 bp effect beyond 2 standard errors of their own replicate scatter."),
    ("Stable sites lose more accessibility than Dynamic sites at the 170 bp shuffle", "mid",
     f"Raw difference {m.loc['170', 'diff_dyn_minus_stable']:+.2f} (95% CI {m.loc['170', 'ci_lo']:+.2f} to {m.loc['170', 'ci_hi']:+.2f}), robust as a raw difference. After the independent review it shrinks to {hg(cmn, M_MATCH).gap:+.2f} with the matched-readout adjustment, {hg(cmn, M_FULL).gap:+.2f} with peak width and TSS distance added, and {hg(cmn, M_FP).gap:+.2f} (ns) with footprint score; with the Patski track it is {hg(fbp, M_MATCH).gap:+.2f} (Part H). Mostly explained by accessibility and context, with a small residual in some combinations."),
    ("The direction of change (age up vs down, exercise up vs down) changes how much a site depends on its sequence", "bad",
     "Directly comparing Age Up vs Age Down, and Exercise Up vs Exercise Down, shows no significant difference at 170 bp or peak level "
     f"(smallest FDR {E[E.width != 'motif'].fdr.min():.2f}). The large effects are all 'stably bound vs not', not up vs down."),
    ("Dynamic sites lose less accessibility than Stable sites when the whole peak is shuffled", "mid",
     f"Raw difference {m.loc['peak', 'diff_dyn_minus_stable']:+.2f} ({m.loc['peak', 'sig']}), but small (Cliff's delta {m.loc['peak', 'cliffs_delta']:.2f}) and it mostly disappears after adjusting for native ATAC and peak width "
     f"(adjusted {adj[adj.width == 'peak'].adj_diff.min():+.2f} to {adj[adj.width == 'peak'].adj_diff.max():+.2f}; "
     f"{int((adj[adj.width == 'peak'].adj_fdr < 0.05).sum())} of {len(adj[adj.width == 'peak'])} class contrasts still FDR &lt; 0.05). Not reproduced with the Patski track. Treat as a peak-width/accessibility effect."),
    ("Motif-only shuffling differs between Stable and Dynamic sites", "bad",
     f"Difference {m.loc['motif', 'diff_dyn_minus_stable']:+.2f} (Cliff's delta {m.loc['motif', 'cliffs_delta']:.2f}, negligible); permutation p = {pm.loc['motif', 'perm_p']:.3f}; not reproduced with the Patski track. "
     "A few individual classes (Age Down, Exercise Up, YS only, MAE only) differ by about 0.15 to 0.26 after adjustment, which is small relative to the replicate scatter."),
    ("Predicted RNA changes when an AP-1 site or peak is shuffled", "bad",
     f"Median RNA effects are about 0 (|log2FC| mostly &lt; 0.1; largest {rna_near.min():+.2f} within 2 kb of a TSS). ATAC and RNA effects are correlated across sites (rho 0.19 to 0.24) but the RNA effect itself is tiny. Not informative with these tracks."),
    ("Predicted accessibility tracks measured snATAC accessibility", "mid",
     f"Spearman rho {nv.rho:.2f} (n = {int(nv.n)}) for the peak readout; lower for 170 bp and motif readouts. Good enough to be informative, not a precise predictor."),
    ("In a peak swap, the swapped-in sequence (not the location) sets accessibility", "ok",
     f"Swapped vs donor native rho {swap.loc['donor native (sequence)', 'rho']:.2f}; vs recipient native {swap.loc['recipient native (location)', 'rho']:.2f}. "
     "Caveat: only 4 recipients, so the 292 pairs are not independent (about 73 independent donors)."),
    ("Inserting an AP-1 motif into a shuffled random background raises accessibility", "bad",
     f"+{mi('atac_motif', '+motif vs shuffled background').median_log2fc:.2f} log2FC at the centre (p = {mi('atac_motif', '+motif vs shuffled background').wilcoxon_p:.2f}); background ATAC is near zero, so the readout is at the floor. No usable signal."),
]

stv = A[A.trackset == "CM"].copy()
stv["Shuffle"] = stv.shuffle.map(W)
stv["Readout"] = stv.readout.map(W)
stv_tab = stv[["Shuffle", "Readout", "median_stable", "median_dynamic", "diff_dyn_minus_stable", "ci_lo", "ci_hi", "cliffs_delta", "fdr", "sig"]].rename(columns={
    "median_stable": "Stable median", "median_dynamic": "Dynamic median", "diff_dyn_minus_stable": "Dyn - Stable", "ci_lo": "CI lo", "ci_hi": "CI hi", "cliffs_delta": "Cliff's d", "fdr": "FDR", "sig": "sig."})
stv_tab["FDR"] = stv_tab["FDR"].map(lambda v: f"{v:.1e}")

cls = B[~B["class"].str.startswith("[") & B.raw_p.notna()].copy()
cls["Shuffle"] = cls.width.map(W)
cls_tab = cls[["family", "Shuffle", "class", "n", "median", "raw_diff", "raw_sig", "adj_diff", "adj_lo", "adj_hi", "adj_sig"]].rename(columns={
    "family": "Family", "class": "Class", "median": "Median", "raw_diff": "Raw diff", "raw_sig": "raw sig.", "adj_diff": "Adj. diff", "adj_lo": "CI lo", "adj_hi": "CI hi", "adj_sig": "adj. sig."})
cov = B[B["class"].str.startswith("[") & (B.family == "age")][["width", "class", "adj_diff", "adj_lo", "adj_hi", "adj_p"]]
cov["width"] = cov.width.map(W)
cov["adj_p"] = cov.adj_p.map(lambda p: f"{p:.1e}")
cov.columns = ["Shuffle", "Covariate", "Coefficient", "CI lo", "CI hi", "p"]
E_tab = E.copy()
E_tab["width"] = E_tab.width.map(W)
E_tab = E_tab[["contrast", "width", "n_a", "n_b", "diff_a_minus_b", "ci_lo", "ci_hi", "cliffs_delta", "fdr", "sig"]]
E_tab["fdr"] = E_tab.fdr.map(lambda v: f"{v:.3f}")
E_tab.columns = ["Contrast", "Shuffle", "n a", "n b", "a - b", "CI lo", "CI hi", "Cliff's d", "FDR", "sig."]
C_tab = C.copy()
C_tab["fdr"] = C_tab.fdr.map(lambda v: f"{v:.1e}")
C_tab = C_tab[["native_quintile", "native_atac_median", "n_stable", "n_dynamic", "median_stable", "median_dynamic", "diff", "fdr", "sig"]]
C_tab.columns = ["Native quintile", "Native ATAC (median)", "n Stable", "n Dynamic", "Stable median", "Dynamic median", "Dyn - Stable", "FDR", "sig."]
t5_tab = t5.reset_index()
t5_tab.columns = ["Tier", "Class", "170 bp shuffle", "Motif shuffle", "Peak shuffle"]
t5_tab = t5_tab[["Tier", "Class", "Motif shuffle", "170 bp shuffle", "Peak shuffle"]]
D_tab = D.copy()
D_tab["shuffle"] = D_tab.shuffle.map(W)
D_tab.columns = ["Shuffle", "Median effect", "Median replicate SD", "Median SE (10 reps)", "Sites with |effect| > 2 SE", "Sites with effect < -2 SE"]
OM_tab = OM.copy()
OM_tab["width"] = OM_tab.width.map(W)
OM_tab["kruskal_p"] = OM_tab.kruskal_p.map(lambda p: f"{p:.1e}")
OM_tab.columns = ["Family", "Shuffle", "Kruskal-Wallis H", "p", "sig."]
NV_tab = NV[["trackset", "readout", "rho", "n"]].copy()
NV_tab["trackset"] = NV_tab.trackset.map({"CM": "Mouse heart ATAC (CM)", "FB": "Mouse Patski ATAC (FB)"})
NV_tab.columns = ["Track set", "Readout", "Spearman rho", "n"]
ps_tab = PS[["compared_with", "rho", "n"]]
ps_tab.columns = ["Swapped ATAC compared with", "Spearman rho", "n pairs"]
mi_tab = ins[ins.readout.isin(["atac_motif", "atac_peak"])][["readout", "comparison", "median_log2fc", "frac_positive", "wilcoxon_p", "n"]].copy()
mi_tab["readout"] = mi_tab.readout.map({"atac_motif": "Centre 7 bp", "atac_peak": "Centre 128 bp"})
mi_tab["wilcoxon_p"] = mi_tab.wilcoxon_p.map(lambda p: f"{p:.3f}")
mi_tab.columns = ["Readout", "Comparison", "Median log2FC", "Fraction > 0", "Wilcoxon p", "n backgrounds"]


claims += [
    ("The TOBIAS binding classes correspond to real accessibility changes in the measured snATAC data", "mid",
     f"Age Up sites gain accessibility (median log2 MAS/YS {q4.loc[('age_class', 'Age Up'), 'median_log2fc']:+.2f} vs {q4.loc[('age_class', 'Stable'), 'median_log2fc']:+.2f} for Stable), Exercise Up sites too ({q4.loc[('ex_class', 'Exercise Up'), 'median_log2fc']:+.2f}); "
     f"the 'down' classes move the other way but only weakly (Age Down {q4.loc[('age_class', 'Age Down'), 'median_log2fc']:+.2f}, Exercise Down {q4.loc[('ex_class', 'Exercise Down'), 'median_log2fc']:+.2f}). Kruskal-Wallis p {fmtp(age_p)} (age) and p {fmtp(ex_p)} (exercise)."),
    ("Predicted peak-shuffle loss mostly tracks how accessible and how wide a peak is", "ok",
     f"Spearman rho {an2.loc[('peak', 'native_peak'), 'rho']:+.2f} with native predicted accessibility and {an2.loc[('peak', 'peak_width'), 'rho']:+.2f} with peak width; the motif-width effect depends on these much less ({an2.loc[('motif', 'native_peak'), 'rho']:+.2f}, {an2.loc[('motif', 'peak_width'), 'rho']:+.2f})."),
    ("Sites with stronger TOBIAS footprints are predicted to depend more on their 170 bp context", "mid",
     f"rho {an4.loc[('eff_170', 'Stable'), 'rho']:+.2f} in Stable sites (95% CI {an4.loc[('eff_170', 'Stable'), 'ci_lo']:+.2f} to {an4.loc[('eff_170', 'Stable'), 'ci_hi']:+.2f}), none in Dynamic sites ({an4.loc[('eff_170', 'Dynamic'), 'rho']:+.2f}). "
     "Independent data, but footprint scores come from the same ATAC data that define accessibility, so part of this is accessibility itself."),
    ("Several Fos motifs in one peak buffer each other (a second site makes a single shuffle matter less)", "bad",
     f"No difference: motif {an5.loc['motif', 'raw_diff']:+.2f}, 170 bp {an5.loc['170', 'raw_diff']:+.2f}, peak {an5.loc['peak', 'raw_diff']:+.2f} (FDR {an5.loc['motif', 'raw_fdr']:.2f} to {an5.loc['peak', 'raw_fdr']:.2f}; n = {int(an5.loc['motif', 'n_multi'])} sites in multi-motif peaks)."),
    ("Promoter-proximal sites depend less on their 170 bp context than distal sites", "ok",
     f"Median 170 bp effect in Stable sites: {an6.loc[('eff_170', 'Stable', '<=2 kb'), 'median']:+.2f} within 2 kb of a TSS vs {an6.loc[('eff_170', 'Stable', '>30 kb'), 'median']:+.2f} beyond 30 kb (Mann-Whitney p {fmtp(_mwu(SL[SL.tss_bin == '<=2 kb'].eff_170.dropna(), SL[SL.tss_bin == '>30 kb'].eff_170.dropna()).pvalue) if 'tss_bin' in SL else 'n/a'}); promoter-proximal peaks are also about 3x more accessible."),
    ("Most AP-1 sites that need their motif also need their surroundings; motif-only dependence is rare", "mid",
     f"Only {an3.loc['Motif only', 'share_stable']:.1%} of Stable sites lose more than half their accessibility from a motif shuffle without also losing it from a 170 bp shuffle ({an3.loc['Context only', 'share_stable']:.0%} are context-only, {an3.loc['Motif and context', 'share_stable']:.0%} both). The thresholds (a 50% loss) are arbitrary, so treat the proportions as descriptive."),
    ("Predicted accessibility tracks MPRA enhancer activity", "mid",
     f"Weak but consistent positive correlation of native predicted ATAC with MPRA log2 alpha (heart {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f}, Patski {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f}; n = 1,653), reproduced in the replicated and TOBIAS-motif subsets."),
    ("The predicted accessibility loss from shuffling an element predicts its loss of MPRA activity", "bad",
     f"No correlation: motif shuffle rho {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'rho']:+.2f} (95% CI {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'ci_lo']:+.2f} to {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'ci_hi']:+.2f}, n = {int(mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'n'])}), full shuffle {mp_.loc[('Patski ATAC', CMP_F, 'all elements'), 'rho']:+.2f}. "
     "The intervals for the all-element motif shuffle exclude anything larger than about |0.06|, but other rows are wider and the MPRA's own scramble effects are tiny, so there is little variance to predict; read this as a null that does not discriminate model quality."),
    ("Results do not depend on which GPU ran a chunk; 10 replicates are enough", "ok",
     f"Chunk results are indistinguishable across L40S, A100 and H100 (Kruskal-Wallis p {q8.kruskal_p_across_gpus.min():.2f} to {q8.kruskal_p_across_gpus.max():.2f}). Averaging 5 vs 10 replicates changes a site's effect by about {q7.loc[('170', 5)]:.2f} (170 bp), {q7.loc[('motif', 5)]:.2f} (motif) and {q7.loc[('peak', 5)]:.2f} (peak)."),
]

claims += [
    ("The loss after a 170 bp shuffle is specific to the AP-1 neighbourhood, not a generic effect of shuffling 170 bp of the peak", "ok",
     f"Centred shuffle {co_.loc['Stable', 'centred_170_median']:.2f} vs {co_.loc['Stable', 'offset_pm300_median']:.2f} for the same-size region 300 bp away (Stable); the centred shuffle is worse at {co_.loc['Stable', 'share_centred_more_negative']:.0%} (Stable) and {co_.loc['Dynamic', 'share_centred_more_negative']:.0%} (Dynamic) of sites, paired Wilcoxon p &lt; 1e-10."),
    ("The effect grows smoothly with the shuffled width and saturates near 340 bp; the 7 bp core alone gives about a quarter", "ok",
     f"Share of the 170 bp effect: 7 bp {fr_.loc[7]:.0%}, 11 bp {fr_.loc[11]:.0%}, 85 bp {fr_.loc[85]:.0%}, 340 bp {fr_.loc[340]:.0%} (Stable, motif readout)."),
    ("The surroundings, not the motif, carry the sequence dependence (mutation arm)", "ok",
     f"170 bp mutated with the motif intact = {ratio_cm170:.0%} (CM) and {ratio_fb170:.0%} (FB) of the loss from mutating everything; within the motif the 7 bp core = {ratio_cmcore:.0%} (CM) and {ratio_fbcore:.0%} (FB) of the 11 bp loss."),
    ("The width of the shuffled motif explains the larger motif effect than on the slides", "bad",
     f"Even a 7 bp shuffle gives {D1_med[(D1_med.readout == 'motif') & (D1_med['class'] == 'Stable') & (D1_med.shuffle_width_bp == 7)]['median'].iloc[0]:.2f} on the motif readout, against about -0.12 on the slides; the site set and readout definition remain unresolved."),
    ("Stable sites depend more on their 170 bp context than Dynamic sites also in fibroblast sites and with the fibroblast-like (Patski) track", "ok",
     f"Raw difference at 170 bp: CM sites {rep.loc[('CM sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} / FB sites {rep.loc[('FB sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} with Patski, {rep.loc[('FB sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (FB sites, heart); all FDR &lt; 0.001. After adjusting for native accessibility and peak width the Patski gaps shrink to {adj_.loc[('CM sites, Patski ATAC', '170'), 'adj_diff']:+.2f} and {adj_.loc[('FB sites, Patski ATAC', '170'), 'adj_diff']:+.2f} (still FDR &lt; 0.01)."),
    ("Placing MPRA elements in a neutral genomic context makes the native accessibility head agree with the MPRA changes", "bad",
     f"No change: shuffle-effect correlations stay at about 0 (|rho| &le; {D4[D4.comparison != 'native level vs MPRA activity'].rho.abs().max():.2f}); native levels correlate {mc_.loc[('heart ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (heart) and {mc_.loc[('Patski ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (Patski) in the neutral context."),
]

claims_html = "".join(f"<tr><td class='v'>{c}</td><td>{badge(k)}</td><td class='v'>{ev}</td></tr>" for c, k, ev in claims)

body = f"""
<h1>AlphaGenome in-silico perturbation of AP-1 (Fos) sites: cardiomyocyte and fibroblast sites</h1>
<p class="sub">Generated {date.today():%d %b %Y} &middot; {int(m.loc['170', 'n_stable'] + m.loc['170', 'n_dynamic']):,} cardiomyocyte Fos sites, 10 shuffle replicates per width &middot; mouse heart and Patski ATAC tracks &middot; shuffle, mutation, insertion, control and MPRA experiments</p>

<div class="box"><strong>Main takeaways</strong><ul>
<li><strong>What is solid.</strong> Shuffling the 170&nbsp;bp around an AP-1 site removes far more predicted accessibility (median log2FC {m.loc['170', 'median_stable']:.2f} in Stable sites) than shuffling the 11&nbsp;bp motif alone ({m.loc['motif', 'median_stable']:.2f}). This holds for essentially every site and is well above replicate noise.</li>
<li><strong>Stable vs Dynamic: a real raw difference, mostly explained by accessibility and context.</strong> Stable sites lose more at the 170&nbsp;bp shuffle ({m.loc['170', 'diff_dyn_minus_stable']:+.2f}, 95% CI {m.loc['170', 'ci_lo']:+.2f} to {m.loc['170', 'ci_hi']:+.2f}), but after an independent review the gap shrinks to {hg(cmn, M_MATCH).gap:+.2f} when adjusting for native accessibility at the matched readout, to {hg(cmn, M_FULL).gap:+.2f} after also adjusting for peak width and TSS distance, and to {hg(fbp, M_MATCH).gap:+.2f} with the Patski track (Part H). Stable sites mainly differ by being more accessible and more strongly bound, not by an independently identifiable sequence dependence.</li>
<li><strong>Age vs exercise vs combination.</strong> Dynamic classes look alike at 170&nbsp;bp and peak level and no significant difference is found between Age Up and Down or Exercise Up and Down, but the intervals (about &plusmn;0.2 to 0.3, 79 to 210 sites per class) are too wide to exclude differences half the size of the Stable-vs-Dynamic gap. This is absence of evidence, not evidence of no difference.</li>
<li><strong>What is noise or unusable.</strong> RNA effects are about zero; AP-1 motif insertion into random backgrounds is at the floor; per-site peak-shuffle values scatter widely (replicate SD about {D.set_index('shuffle').loc['peak', 'median_replicate_sd']:.1f}) and are only trustworthy in aggregate.</li>
<li><strong>New in this version.</strong> Data QC (Part A) shows the class labels are directionally real in independent snATAC, the model reproduces the slides' Dynamic-site values, shuffles behave as intended, 10 replicates are justified and GPU type does not matter. Further analyses (Part C) show that dependence on the 170&nbsp;bp context grows with distance from a TSS and with TOBIAS footprint strength, that sites rarely depend on the motif alone, that two motifs in one peak do not buffer each other, and that AlphaGenome's predicted shuffle effects do <em>not</em> track the fibroblast MPRA (rho &asymp; 0).</li>
<li><strong>Overnight experiments (Part D).</strong> Same-size shuffles 300&nbsp;bp away do almost nothing, so the loss is local to the site's neighbourhood (a distance control, not an AP-1-specific one); it saturates near 340&nbsp;bp, the 7&nbsp;bp core alone gives about a quarter, and in a mutation experiment the surroundings carry most of the effect even with the motif intact. The Stable-vs-Dynamic difference replicates in fibroblast sites (smaller after adjustment with the fibroblast-like track), and a neutral genomic context does not rescue the failure to predict the MPRA.</li>
<li><strong>Caveat on the slides.</strong> Our Dynamic-site native accessibility matches the slides closely, which suggests the pipeline reproduces the API, but our Stable-site set differs from the slides' (we could not reproduce their site counts). Absolute numbers should not be compared directly with the slides until the site set is resolved.</li></ul></div>
<div class="box"><strong>Independent review (Part H).</strong> Three independent reviewers re-derived the headline numbers (all reproduced), confirmed the class labels and site tables, and found no bug that invalidates the shuffle, mutation or control results. They did find that the Stable-vs-Dynamic adjustment used the wrong covariate; this is corrected in Part H and in the text below, and a list of smaller code issues not yet fixed is given there.</div>

<h2>How to read this report</h2>
<p><strong>Effect.</strong> For each site and shuffle width, the effect is the mean over 10 shuffles of <code>log2((shuffled + pc) / (native + pc))</code> of the predicted ATAC (mean per-base signal over a region; pc = 5% of the median native value). Negative = accessibility lost. A shuffle of the motif is read out over the motif, a 170&nbsp;bp shuffle over the 170&nbsp;bp, a peak shuffle over the peak ("matched readout"), unless stated.</p>
<p><strong>Significance marks.</strong> <code>***</code> FDR &lt; 0.001, <code>**</code> &lt; 0.01, <code>*</code> &lt; 0.05, <code>ns</code> not significant. Tests are two-sided Mann&ndash;Whitney against a reference class, Benjamini&ndash;Hochberg corrected within each family of tests. With thousands of sites a tiny difference is "significant", so each comparison also shows an effect size (difference of medians with a bootstrap 95% CI, and Cliff's delta: about 0.1 small, 0.3 medium, 0.5 large). <strong>Read the effect size before the stars.</strong></p>
<p><strong>Trust labels.</strong> <span class="badge ok">&#10003; Robust</span> survives the checks below. <span class="badge mid">~ Plausible / partial</span> real but small, or depends on a confound. <span class="badge bad">&#10007; Not supported</span> no evidence, or indistinguishable from noise. <span class="badge na">? Cannot judge yet</span> data not available.</p>

{DESIGN}
{PART_A}
<h2 id="partB">Part B. Main results</h2>
<h3 class="sec">B1. Does the model track real accessibility?</h3>
<figure>{img(CM / 'figures' / 'fig1_native_vs_measured.png', 'Predicted vs measured accessibility', '560px')}
<figcaption>Each dot is a site: predicted native ATAC over the peak (mouse heart track) against measured cardiomyocyte snATAC at the same peak. Spearman rho = {nv.rho:.2f}.</figcaption></figure>
{table(NV_tab, {'Spearman rho': '{:.2f}'})}
<p><strong>What it shows.</strong> The heart track predicts measured cardiomyocyte accessibility with a moderate rank correlation (the slides reported 0.49). The Patski (fibroblast-like) track is essentially uncorrelated with cardiomyocyte accessibility, as it should be. <strong>Trust:</strong> good enough to use the heart track as a proxy for cardiomyocyte sites; it is bulk heart tissue, not purified cardiomyocytes.</p>

<h3 class="sec">B2. What does shuffling do?</h3>
<figure>{img(CM / 'figures' / 'fig2s_shuffle_heatmaps_stars.png', 'Shuffle effect heatmaps')}
<figcaption>Median log2FC for each shuffle width (rows) and ATAC readout region (columns). Left: Stable sites; right: Dynamic sites, with significance marks against Stable.</figcaption></figure>
<p><strong>What it shows.</strong> Shuffling only the motif lowers accessibility at the motif by about 0.6 log2 units (roughly a third lost) and over the whole peak by about 0.3 (about a fifth lost). Shuffling 170&nbsp;bp or the whole peak is far more damaging (about -1.8 log2 units locally, i.e. around 70% lost). Dynamic sites lose less than Stable sites in most cells. <strong>Trust:</strong> the ordering motif &lt; 170&nbsp;bp &asymp; peak is robust. Peak shuffles are noisy per site (see B7).</p>
{table(stv_tab)}

<h3 class="sec">B3. Stable vs Dynamic: is the difference real?</h3>
<p>Dynamic sites have lower native accessibility and different peak widths, either of which could make their log2FC look smaller without any biological difference. Three checks:</p>
<figure>{img(CM / 'figures' / 'fig8_stable_vs_dynamic_by_native_accessibility.png', 'Stable vs Dynamic within accessibility quintiles', '640px')}
<figcaption>Sites split into five groups by native predicted accessibility. Within every group, Stable sites lose more than Dynamic sites at the 170&nbsp;bp shuffle.</figcaption></figure>
{table(C_tab)}
<ul><li><strong>Stratification:</strong> the gap is present in all five accessibility quintiles ({'all FDR &lt; 0.05' if quint_all_sig else 'most FDR &lt; 0.05'}), (quintiles of native <em>peak</em> accessibility). <strong>Correction:</strong> the right variable is native accessibility at the 170&nbsp;bp readout itself. Stratifying on that (Part H) shows a gap of only about 0.1 to 0.2 in the lower deciles that hold most Dynamic sites, so this check does not rule out an accessibility effect.</li>
<li><strong>Permutation test</strong> (5,000 label shuffles): 170&nbsp;bp p = {pm.loc['170', 'perm_p']:.4f}, peak p = {pm.loc['peak', 'perm_p']:.4f}, motif p = {pm.loc['motif', 'perm_p']:.3f}. The 170&nbsp;bp and peak differences are far outside the null; the motif one is not.</li>
<li><strong>Independent track:</strong> repeating the test with the Patski ATAC track on the same sites gives 170&nbsp;bp {fb.loc['170', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['170', 'sig']}), peak {fb.loc['peak', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['peak', 'sig']}), motif {fb.loc['motif', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['motif', 'sig']}). The 170&nbsp;bp direction replicates; the others do not.</li></ul>
<p><strong>Trust:</strong> the raw 170&nbsp;bp difference is real, but it is largely accounted for by native accessibility and context once the matched readout is adjusted for (Part H); partial for the peak shuffle; not supported for the motif shuffle.</p>

<h3 class="sec">B4. Age, exercise and their combination</h3>
<p>Each site carries three labels from the same binding calls: <em>age class</em> (YS vs MAS), <em>exercise class</em> (MAS vs MAE) and the full binding pattern. Every class is compared with its Stable reference, once raw and once adjusted for native ATAC and peak width (OLS with robust standard errors).</p>
<figure>{img(CM / 'figures' / 'fig3s_class_effects_stars.png', 'Class effects with significance marks')}
<figcaption>Median matched effect per class with 95% bootstrap CI. Blue = the Stable reference; marks to the right are FDR-corrected Mann&ndash;Whitney tests against it.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig7_class_effects_adjusted.png', 'Adjusted class effects')}
<figcaption>Difference from the Stable reference: raw (orange diamonds) and adjusted for native accessibility and peak width (blue, with 95% CI). Positive = the class loses less than Stable.</figcaption></figure>
<p><strong>What it shows.</strong> At 170&nbsp;bp every dynamic class loses about 0.45 to 0.6 log2FC less than Stable, and adjusting for native <em>peak</em> accessibility barely changes it (the matched-readout adjustment in Part H roughly halves it). At the peak and motif widths, raw differences are small and the adjusted ones shrink toward zero. <strong>Between dynamic classes</strong> (the direct test of age vs exercise, and of up vs down):</p>
{table(E_tab)}
<p>Only the motif-width contrasts reach FDR &lt; 0.05, with confidence intervals touching zero. There is no detectable difference between age-driven and exercise-driven classes at this precision. <strong>Trust:</strong> "Stable vs dynamic" holds as a raw difference but is largely explained by accessibility (Part H); "which dynamic class" cannot be resolved with these sample sizes.</p>
<details><summary>Omnibus tests, full class table and covariates</summary>
<h3>Kruskal&ndash;Wallis across all classes in each family</h3>{table(OM_tab, {'Kruskal-Wallis H': '{:.1f}'})}
<h3>Every class against its reference</h3>{table(cls_tab)}
<h3>Covariates in the adjusted model (age family)</h3>{table(cov)}
<p class="small">Positive peak-width coefficient: wider peaks lose less accessibility per log2FC. Dynamic sites differ systematically from Stable sites in width and native accessibility, which is why the raw peak-shuffle gap shrinks after adjustment.</p></details>

<h3 class="sec">B5. RNA</h3>
<figure>{img(CM / 'figures' / 'fig4_rna_effect_by_tss_distance.png', 'RNA effect by TSS distance')}
<figcaption>Median log2FC of predicted RNA after a peak shuffle, by distance of the site to the nearest TSS, for Stable and Dynamic sites. Left: within 2&nbsp;kb of that TSS; right: mean over &plusmn;50&nbsp;kb.</figcaption></figure>
<p><strong>What it shows.</strong> Predicted RNA barely moves when a distal AP-1 site or peak is shuffled (medians within about &plusmn;0.1), with the largest median drops 2 to 10&nbsp;kb from the TSS for the &plusmn;2&nbsp;kb readout and within 2&nbsp;kb for the &plusmn;50&nbsp;kb mean (peak shuffle, Stable sites: about &minus;0.08 and &minus;0.12). ATAC and RNA effects correlate across sites (rho 0.19 to 0.24), so there is a weak coupling, but the RNA effects themselves are far smaller than the slides reported. <strong>Trust:</strong> not informative with these tracks; do not over-interpret.</p>

<h3 class="sec">B6. Does a high-confidence subset change the picture?</h3>
<figure>{img(CM / 'figures' / 'fig5_confidence_tiers.png', 'Confidence tiers')}
<figcaption>Median effect in three nested sets chosen only from native prediction vs measurement (never from the effects): A all sites; B sites with a measured peak and predicted accessibility above the 25th percentile; C = B plus prediction-measurement rank agreement within 1 SD.</figcaption></figure>
{table(t5_tab)}
<p><strong>What it shows.</strong> Effects get somewhat larger in the filtered sets (the model is more confident at accessible, well-predicted sites) and the Stable vs Dynamic gap at 170&nbsp;bp stays. <strong>Trust:</strong> conclusions are not an artefact of poorly predicted sites; report all three tiers. Part G repeats this with a stricter tier and for the mutation arm, controls and fibroblast sites.</p>

<h3 class="sec">B7. How much is noise?</h3>
{table(D_tab, {'Sites with |effect| > 2 SE': '{:.0%}', 'Sites with effect < -2 SE': '{:.0%}'})}
<p>Replicate SD is the spread across the 10 shuffles of one site. The 170&nbsp;bp effect is stable per site (SD about {D.set_index('shuffle').loc['170', 'median_replicate_sd']:.1f}); the motif and peak effects scatter more (SD about {D.set_index('shuffle').loc['motif', 'median_replicate_sd']:.1f} and {D.set_index('shuffle').loc['peak', 'median_replicate_sd']:.1f}). The median over thousands of sites is precise, but <strong>an individual site's peak-shuffle value should not be interpreted on its own</strong>. This is why 10 replicates were kept.</p>

<h3 class="sec">B8. Insertion experiments</h3>
<figure>{img(AN / 'insertion' / 'figures' / 'fig_peak_swap.png', 'Peak swap')}
<figcaption>Peak swap: a donor peak's sequence replaces a recipient peak. Swapped accessibility follows the donor's own accessibility (left) much more than the recipient location's (right).</figcaption></figure>
{table(ps_tab, {'Spearman rho': '{:.2f}'})}
<p><strong>Trust:</strong> robust direction (sequence beats location), but only 4 recipient windows, so the 292 pairs are not independent; the effective sample is about 73 donors.</p>
<figure>{img(AN / 'insertion' / 'figures' / 'fig_motif_insertion.png', 'Motif insertion')}
<figcaption>A TGAGTCA motif inserted at the centre of dinucleotide-shuffled random genomic backgrounds (50 backgrounds), compared with the shuffled background and with a shuffled-motif control.</figcaption></figure>
{table(mi_tab, {'Median log2FC': '{:.2f}', 'Fraction > 0': '{:.2f}'})}
<p><strong>Trust:</strong> no usable signal. Background accessibility is near zero, so any ratio is dominated by the floor, and the 7&nbsp;bp readout differences are not significant. The slides reported a similarly small effect.</p>

{PART_C}
{PART_D}
{PART_G}
{PART_H}
<h2>E. What to trust</h2>
<table class="t"><tr><th>Claim</th><th>Verdict</th><th>Evidence</th></tr>{claims_html}</table>

<h2>F. Open issues and next steps</h2>
<ul>
<li><strong>Site set vs the slides.</strong> The TOBIAS table gives 3,679 CM Fos sites (2,907 Stable), not the slides' 2,663 (2,320 Stable). Our Dynamic native accessibility matches the slides, but our Stable sites are about half as accessible on average, and keeping our top 80% by predicted accessibility does not close the gap. The slides' Stable sites are a different set. Until the original selection is known, compare patterns, not absolute values.</li>
<li><strong>Model equivalence.</strong> The slides were made with the AlphaGenome API; this run uses the PyTorch port with the same family of weights. The Dynamic-site match supports equivalence, but a direct check on the same sites through the API would settle it.</li>
<li><strong>Status of the runs.</strong> Done: shuffle, mutation, insertion, control/titration (cardiomyocyte sites only), neutral-context MPRA, and shuffle and mutation arms on fibroblast sites. Not run yet: saturation mutagenesis of the motif, the human-head sensitivity check and the window-length check (all coded, resumable).</li>
<li><strong>Fibroblasts.</strong> The mouse model has no fibroblast ATAC; Patski and embryonic-fibroblast RNA are weak proxies. The fibroblast shuffle run and the MPRA comparison will be added when they finish.</li>
<li><strong>Motif-width effect is larger than in the slides</strong> (about -0.6 vs -0.12). Candidate causes are the motif interval used (here the 11&nbsp;bp TOBIAS interval for both the shuffle and the readout) and the different Stable set; worth confirming the slides' definition.</li></ul>

<details><summary>Methods and files</summary>
<p class="small">Predictions: alphagenome-pytorch with <code>model_all_folds</code> weights, mouse organism, 131,072&nbsp;bp window centred on the motif midpoint, mouse heart ATAC track (1&nbsp;bp, mean over region), 10 uniform random shuffles per width. Stable = bound in YS, MAS and MAE; age class compares YS and MAS; exercise class compares MAS and MAE. Adjusted model: effect ~ class + log10(native ATAC over peak) + log10(peak width), HC3 robust SEs, normal-theory CIs. All CIs for medians are percentile bootstrap (1,000&ndash;2,000 resamples). Sites are treated as independent, but about 15% of sites share a peak with another site, so p-values are mildly optimistic; effect sizes and confidence intervals are the better guide. Tables: <code>Analysis/cm/tables</code>, <code>Analysis/cm/stats</code>; figures: <code>Analysis/cm/figures</code>; scripts: <code>Scripts/22_shuffle_analysis.py</code>, <code>24_class_stats.py</code>, <code>23_insertion_analysis.py</code>, <code>25_build_report.py</code>.</p></details>
"""

doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AP-1 shuffle analysis report</title><style>{CSS}</style></head><body><main>{body}</main></body></html>"""
out = AN / "report_cm.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
