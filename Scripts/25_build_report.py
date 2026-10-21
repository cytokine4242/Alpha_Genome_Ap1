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
from edit_labels import EDIT_LABEL, BASES_CHANGED, PEAK_MEDIAN  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "Analysis"
CM = AN / "cm"


IMGS = {}


def img(path, alt, width="100%"):
    k = len(IMGS)
    IMGS[k] = Path(path)
    return f'<img src="@@IMG{k}@@" alt="{html.escape(alt)}" style="max-width:{width}">'


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
table.t.reg{display:table;width:100%}table.t.reg td,table.t.reg th{white-space:normal;vertical-align:top;text-align:left}table.t.reg td.v{min-width:0!important}table.t.reg .badge{white-space:normal;font-size:.72rem}a{color:var(--blue)}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:1000px;margin:0 auto;padding:32px 16px 80px}h1{font-size:1.9rem;line-height:1.2;margin:0 0 4px}h2{font-size:1.35rem;margin:44px 0 8px;padding-top:12px;border-top:1px solid var(--rule)}h3{font-size:1.05rem;margin:26px 0 4px}figure.dia{margin:14px 0}.dia{overflow-x:auto}.dia svg{min-width:760px;width:100%;height:auto;display:block;color:var(--ink)}.dia figcaption{position:sticky;left:0}h3.sec{font-size:1.15rem;margin:34px 0 6px;padding-top:10px;border-top:1px dotted var(--rule)}
p,li{max-width:78ch}.sub{color:var(--ink2);margin:0 0 20px}.box{background:var(--panel);border-radius:10px;padding:14px 18px;margin:16px 0}.box ul{margin:6px 0;padding-left:20px}
figure{margin:18px 0}figcaption{color:var(--ink2);font-size:.9rem;margin-top:6px;max-width:80ch}img{border-radius:8px;border:1px solid var(--rule);background:#fcfcfb;height:auto}
table.t{border-collapse:collapse;font-size:.86rem;margin:10px 0;display:block;overflow-x:auto;max-width:100%}table.t th,table.t td{padding:4px 10px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap}table.t th{color:var(--ink2);font-weight:600;position:sticky;top:0;background:var(--bg)}table.t td:first-child,table.t th:first-child{text-align:left}
div.bio{border-left:4px solid var(--green);background:var(--panel);padding:10px 16px;margin:14px 0;border-radius:6px;max-width:80ch}div.bio strong:first-child{color:var(--green)}.conf{margin:6px 0 10px;font-size:.92rem}.conf .badge{margin-right:6px}.toc{columns:2;column-gap:32px;font-size:.92rem;margin:8px 0}.toc ul{margin:0 0 8px;padding-left:18px}.toc a{color:var(--blue);text-decoration:none}h2 .num,h3 .num{color:var(--ink2);margin-right:6px}.badge{display:inline-block;font-size:.78rem;font-weight:600;border-radius:99px;padding:1px 9px;border:1.5px solid currentColor;white-space:nowrap}.ok{color:var(--green)}.mid{color:var(--amber)}.bad{color:var(--red)}.na{color:var(--grey)}
td.v{white-space:normal!important;text-align:left!important;min-width:260px}.small{font-size:.86rem;color:var(--ink2)}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:.88em}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2)}
"""


def badge(kind):
    lab = {"ok": "&#10003; High confidence", "mid": "~ Moderate confidence", "bad": "&#10007; Low / unreliable", "na": "? Untested"}[kind]
    return f'<span class="badge {kind}">{lab}</span>'



# ======================= QC / additional-analysis tables =======================
from scipy.stats import chi2_contingency as _chi2, mannwhitneyu as _mwu

QCD = CM / "qc"
rd = lambda n: pd.read_csv(QCD / "tables" / n)


def qimg(name, alt, width="100%"):
    return img(QCD / "figures" / f"{name}.png", alt, width)


Q3, Q4, Q5, Q6, Q7, Q8 = (rd("qc3_native_means.csv"), rd("qc4_class_labels_vs_measured.csv"), rd("qc5_prediction_vs_measured_by_class.csv"),
                          rd("qc6_shuffled_vs_native.csv"), rd("qc7_replicate_convergence.csv"), rd("qc8_gpu_consistency.csv"))
AN1, AN2, AN3, AN4, AN5, AN6, AN7 = (rd("an1_effect_distributions.csv"), rd("an2_effect_vs_covariates_spearman.csv"), rd("an3_site_types.csv"), rd("an4_footprint_score.csv"),
                                     rd("an5_multiple_motifs_in_peak.csv"), rd("an6_tss_proximity.csv"), rd("an7_mpra_comparison.csv"))
XT = AN / "extra"
xr = lambda n: pd.read_csv(XT / "tables" / n)
XCL, XPR, XTS, XOV, XPS = xr("prediction_quality_by_class.csv"), xr("prediction_quality_vs_stable.csv"), xr("prediction_quality_omnibus.csv"), xr("observed_vs_predicted.csv"), xr("peak_swap_stable_vs_dynamic.csv")
xts = XTS.set_index("combo")
xov = XOV.set_index("combo")
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


def ximg(name, alt, width="100%"):
    return img(XT / "figures" / f"{name}.png", alt, width)


def fmtp(p):
    return "&lt; 1e-10" if p < 1e-10 else f"= {p:.1e}"


Q3_tab = Q3.copy()
Q3_tab["readout"] = Q3_tab.readout.map(W)
Q3_tab.columns = ["Class", "Readout", "Mean native ATAC", "Median native ATAC", "n sites"]
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

XPS_tab = XPS.iloc[:6].copy()
XPS_tab = XPS_tab[["quantity", "stable_median", "dynamic_median", "p_mannwhitney"]]
XPS_tab.columns = ["Quantity", "Stable", "Dynamic", "p"]
XOV_tab = XOV.copy()
XOV_tab.columns = ["Site set and track", "Spearman rho", "CI lo", "CI hi", "n sites"]
_x = XPR[XPR.combo == "CM sites, heart ATAC"].copy()
XPR_tab = _x[["class", "n", "rho_diff_vs_stable", "rho_diff_lo", "rho_diff_hi", "fdr_rho", "median_signed_resid", "fdr_signed_resid"]].merge(XCL[XCL.combo == "CM sites, heart ATAC"][["class", "rho"]], on="class")
XPR_tab = XPR_tab[["class", "n", "rho", "rho_diff_vs_stable", "rho_diff_lo", "rho_diff_hi", "fdr_rho", "median_signed_resid", "fdr_signed_resid"]]
XPR_tab["fdr_rho"] = XPR_tab.fdr_rho.map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
XPR_tab["fdr_signed_resid"] = XPR_tab.fdr_signed_resid.map(lambda v: "" if pd.isna(v) else ("< 0.001" if v < 0.001 else f"{v:.3f}"))
XPR_tab.columns = ["Binding class", "n", "rho", "rho diff vs Stable", "CI lo", "CI hi", "FDR (rho)", "median signed rank error", "FDR (signed error)"]
PART_A = f"""
<h2 id="partA">Part A. Where the data come from (QC)</h2>
<p>These figures show how the site set, the model output and the runs themselves look before any hypothesis is tested. Each ends with what it means for how far the later results can be trusted.</p>

<h3>A1. The site set</h3>
<figure>{qimg('qc1_site_composition', 'Site composition')}<figcaption>Left: sites per TOBIAS binding pattern. Middle: how the age labels (YS vs MAS) and exercise labels (MAS vs MAE) co-occur. Right: sites per chromosome, split into Stable and Dynamic.</figcaption></figure>
<p>{n_all:,} cardiomyocyte Fos sites: {n_st:,} Stable ({n_st / n_all:.0%}) and {n_dy:,} Dynamic, spread across the six dynamic binding patterns (79 to 210 sites each). Age and exercise labels are <strong>not independent</strong>: they come from the same three binding calls, so there are only seven patterns (for example every Age Up site is either Exercise Stable or Exercise Down). Chromosome counts follow chromosome size, and the Stable/Dynamic mix does not differ between chromosomes (chi-square p = {chrom_p:.2f}; chrX has the highest Dynamic share, {chrx['Dynamic'] / (chrx['Dynamic'] + chrx['Stable']):.0%} of only {int(chrx.sum())} sites). <strong>Implication:</strong> class sizes are small (about 100 to 400), which limits the power of every class comparison.</p>

<h3>A2. Site geometry</h3>
<figure>{qimg('qc2_site_geometry', 'Site geometry')}<figcaption>Peak width, distance to the nearest TSS and number of Fos motifs in the same peak, for Stable and Dynamic sites.</figcaption></figure>
<p>Dynamic sites sit in narrower peaks (median {gm.loc['Dynamic', 'peak_width']:.0f} vs {gm.loc['Stable', 'peak_width']:.0f} bp, Mann&ndash;Whitney p {fmtp(gp['peak_width'])}), but at a similar distance from TSSs (median {gm.loc['Dynamic', 'abs_tss']:,.0f} vs {gm.loc['Stable', 'abs_tss']:,.0f} bp, p = {gp['abs_tss']:.2f}). About {(SL.motifs_in_peak >= 2).mean():.0%} of sites share their peak with another Fos site (about {SL.drop_duplicates('peak_coord').motifs_in_peak.ge(2).mean():.0%} of peaks hold two or more). <strong>Implication:</strong> peak width differs systematically between classes and also drives the peak-shuffle effect (Part C), which is why the class comparisons were adjusted for it.</p>

<h3>A3. Native predicted accessibility, and how it compares with the slides</h3>
<figure>{qimg('qc3_native_distributions', 'Native predicted ATAC distributions')}<figcaption>Distribution of native predicted ATAC (mouse heart track) at the motif, 170 bp and peak regions. Legend: mean and median of each class.</figcaption></figure>
<p>Predicted accessibility is strongly right-skewed: most sites are low (median {q3.loc[('Stable', 'peak'), 'median']:.2f} for Stable and {q3.loc[('Dynamic', 'peak'), 'median']:.2f} for Dynamic sites over the peak) and a minority are very high. Dynamic sites are predicted to be less accessible than Stable sites at every readout (peak mean {q3.loc[('Dynamic', 'peak'), 'mean']:.2f} vs {q3.loc[('Stable', 'peak'), 'mean']:.2f}), which matches their lower measured accessibility (section on class labels).</p>
{table(Q3_tab, {'Mean native ATAC': '{:.3f}', 'Median native ATAC': '{:.3f}'})}

<h3>A4. Are the class labels real? Checking them against independent measured data</h3>
<figure>{qimg('qc4_class_labels_vs_measured', 'Class labels vs measured accessibility change')}<figcaption>Measured cardiomyocyte snATAC (pseudobulk, CPM per kb at the peak containing the site) split by experimental condition. Left: aged vs young sedentary, by age class. Right: aged plus exercise vs aged sedentary, by exercise class.</figcaption></figure>
{table(Q4_tab, {'Median measured log2FC': '{:+.3f}'})}
<p>The classes come from TOBIAS footprint calls, not from accessibility, so measured accessibility is an independent check. <strong>The gains are clear and the losses are weak</strong>: Age Up sites gain accessibility with age (median log2 MAS/YS {q4.loc[('age_class', 'Age Up'), 'median_log2fc']:+.2f} vs {q4.loc[('age_class', 'Stable'), 'median_log2fc']:+.2f} for Stable), as do Exercise Up sites ({q4.loc[('ex_class', 'Exercise Up'), 'median_log2fc']:+.2f}); Age Down and Exercise Down sites change in the expected direction but only slightly. Differences across classes are significant (Kruskal&ndash;Wallis p {fmtp(age_p)} (age) and p {fmtp(ex_p)} (exercise)). <strong>Implication:</strong> the labels reflect real but modest accessibility shifts, so "Dynamic" classes are best read as footprint-defined, not as strongly differentially accessible.</p>

<h3>A5. Does the model predict some classes worse?</h3>
<figure>{qimg('qc5_prediction_vs_measured_by_class', 'Prediction vs measurement by class')}<figcaption>Left: Spearman correlation of predicted peak accessibility with measured accessibility, per class (95% CI). Right: distribution of the rank residual (predicted percentile minus measured percentile).</figcaption></figure>
{table(Q5_tab, {'Spearman rho': '{:.2f}', 'CI lo': '{:.2f}', 'CI hi': '{:.2f}'})}
<figure>{ximg('x1_prediction_quality_by_class', 'Prediction quality by class')}<figcaption>Cardiomyocyte sites, heart track. Left: Spearman correlation between predicted and measured accessibility within each binding class (95% bootstrap CI), with marks for the difference from Stable (bootstrap, BH). Middle: size of the rank error. Right: direction of the rank error (positive = over-predicted).</figcaption></figure>
{table(XPR_tab, {'rho': '{:.2f}', 'rho diff vs Stable': '{:+.2f}', 'CI lo': '{:+.2f}', 'CI hi': '{:+.2f}', 'median signed rank error': '{:+.2f}'})}
<p><strong>Do the classes differ in how well they are predicted?</strong> The correlation is {q5.loc['Stable', 'rho']:.2f} in Stable and {q5.loc['Dynamic', 'rho']:.2f} in Dynamic sites; the difference of {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_vs_stable.iloc[0]:+.2f} (95% CI {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_lo.iloc[0]:+.2f} to {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_hi.iloc[0]:+.2f}) is small but significant. Across the seven binding classes the correlations are heterogeneous (Cochran Q = {xts.loc['CM sites, heart ATAC', 'rho_heterogeneity_Q']:.1f}, df 6, p = {xts.loc['CM sites, heart ATAC', 'rho_heterogeneity_p']:.3f}), driven by the two classes bound in the aged sedentary state (MAS only and MAS+MAE: rho 0.19 and 0.16, FDR &lt; 0.05 against Stable). The <em>size</em> of the rank error does not differ between classes (Kruskal&ndash;Wallis p = {xts.loc['CM sites, heart ATAC', 'kw_abs_resid_p']:.2f}; Stable vs Dynamic p = {xts.loc['CM sites, heart ATAC', 'stable_vs_dynamic_abs_resid_p']:.2f}), but its <em>direction</em> does: Dynamic sites are systematically over-predicted (median signed rank error {xts.loc['CM sites, heart ATAC', 'dynamic_median_signed']:+.2f} vs {xts.loc['CM sites, heart ATAC', 'stable_median_signed']:+.2f} for Stable, p &lt; 1e-10; the same in fibroblast sites with both tracks). In fibroblast sites the heterogeneity is weaker (heart track p = {xts.loc['FB sites, heart ATAC', 'rho_heterogeneity_p']:.2f}, Patski p = {xts.loc['FB sites, Patski ATAC', 'rho_heterogeneity_p']:.2f}). <strong>So no class is predicted so badly that the model is unusable for it, but Dynamic sites are predicted somewhat less well and are biased upwards</strong>, which is one reason the confidence tiers (B6) matter.</p>

<h3>A6. Does shuffling lower accessibility as expected? (the slides' QC)</h3>
<figure>{qimg('qc6_shuffled_vs_native', 'Shuffled vs native')}<figcaption>Mean predicted ATAC after each shuffle (bars) against the native mean (dashed), for each readout region.</figcaption></figure>
{table(Q6_tab, {'Mean native ATAC': '{:.3f}', 'Sites that lose accessibility': '{:.0%}'})}
<p>Shuffling reduces accessibility on average everywhere. A 170 bp shuffle lowers the 170 bp readout in {q6.loc[('170', 'Stable'), 'frac_sites_decreased']:.0%} of Stable sites, but a motif-only shuffle lowers the motif readout in only {q6.loc[('motif', 'Stable'), 'frac_sites_decreased']:.0%}: for about a quarter of sites it raises predicted accessibility, so the motif-only effect is small and heterogeneous. <strong>Implication:</strong> the shuffle experiment behaves as intended and the effect sizes are not an artefact of broken sequences.</p>

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
<p>Predicted native accessibility correlates weakly but consistently with MPRA activity (rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski), and this holds in the replicated and TOBIAS-motif subsets. <strong>The predicted effect of shuffling does not track the MPRA at all</strong> (motif shuffle {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'rho']:+.2f}, full shuffle {mp_.loc[('Patski ATAC', CMP_F, 'all elements'), 'rho']:+.2f} with the Patski track; similar with heart). The MPRA delta is native minus shuffled, so agreement would appear as a <em>negative</em> correlation with the (negative) predicted log2FC; for the all-element motif shuffle the confidence interval excludes any correlation beyond about 0.06, but other rows are wider and the MPRA's own scramble effects are tiny, so there is little variance to predict. Possible reasons: the MPRA measures episomal enhancer activity of a 170&nbsp;bp fragment, while the model predicts chromatin accessibility of the same sequence in its genomic context; the MPRA effects are small and noisy; and the fibroblast tracks are weak proxies. <strong>Implication:</strong> accessibility predicted by AlphaGenome should not be used as a proxy for the MPRA effect of disrupting an element.</p>
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
D1c = ovr("d1c_offsets_by_readout.csv")
_c = D1c[D1c.shuffle_position == "titrate@170"].set_index(["readout", "class"])["median"]
_o = D1c[D1c.shuffle_position.isin(["offset@-300", "offset@+300"])].groupby(["readout", "class"])["median"].mean()
D1c = pd.DataFrame({"median": _c, "off": _o})
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
    mut_rows.append({"Edit": EDIT_LABEL[e], "Bases changed": BASES_CHANGED[e], "CM Stable": dep(CMS, e, "Stable"), "CM Dynamic": dep(CMS, e, "Dynamic"), "CM sig.": stars(depq(CMS, e)),
                     "FB Stable": dep(FBS, e, "Stable"), "FB Dynamic": dep(FBS, e, "Dynamic"), "FB sig.": stars(depq(FBS, e))})
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
<figure>{oimg('d1_controls_titration', 'Controls and width titration')}<figcaption>Left: effect on the motif readout when a 170 bp region is shuffled at different positions (median with 95% CI). Right: dose-response of the effect (log2FC) as the centred shuffle widens. Colour = site class; line weight and marker = the region the effect is read over (thick line with circles = 11 bp motif, thin line with squares = 170 bp, thin line with triangles = whole peak; see the legend). Marks: Stable vs Dynamic (Mann&ndash;Whitney, BH-corrected), motif readout.</figcaption></figure>
<figure>{oimg('d1c_offsets_by_readout', 'Offset controls read over the motif, the 170 bp and the whole peak')}<figcaption>The left plot above repeated with the effect read over the 11 bp motif (as before), over the 170 bp and over the whole peak. Median with 95% bootstrap CI.</figcaption></figure>
{table(co_tab, {'Centred 170 bp (median)': '{:.2f}', 'Offset +-300 bp (median)': '{:.2f}', 'Sites where centred is more negative': '{:.0%}'})}
{table(tt, {'Stable median': '{:.2f}', 'Dynamic median': '{:.2f}', 'Share of the 170 bp effect (Stable)': '{:.0%}', 'Dyn - Stable': '{:+.2f}'})}
<p><strong>The effect is specific to the site's neighbourhood.</strong> Shuffling the centred 170&nbsp;bp costs {co_.loc['Stable', 'centred_170_median']:.2f} log2 units at the motif in Stable sites, while the same-size region 300&nbsp;bp away costs only {co_.loc['Stable', 'offset_pm300_median']:.2f}; the centred shuffle is more damaging at {co_.loc['Stable', 'share_centred_more_negative']:.0%} of Stable and {co_.loc['Dynamic', 'share_centred_more_negative']:.0%} of Dynamic sites (paired Wilcoxon, p &lt; 1e-10). This is a distance control (the offset shuffles lie 215 to 385&nbsp;bp from the readout and further from the peak summit), not an AP-1-specific control; matched shuffles at non-Fos accessible sites would be needed for that. <strong>It grows smoothly with width and levels off near 340&nbsp;bp:</strong> the 7&nbsp;bp core alone gives {fr_.loc[7]:.0%} of the 170&nbsp;bp effect, the whole 11&nbsp;bp motif {fr_.loc[11]:.0%}, 85&nbsp;bp {fr_.loc[85]:.0%}, and 340&nbsp;bp {fr_.loc[340]:.0%}. The Stable&ndash;Dynamic gap is small at the core ({sv_.loc[7, 'diff_dyn_minus_stable']:+.2f}) and builds to {sv_.loc[340, 'diff_dyn_minus_stable']:+.2f} by 340&nbsp;bp, so it comes from the surrounding sequence, not the core.</p>
<p>The same holds across the whole peak: shuffling the centred 170&nbsp;bp lowers predicted accessibility over the entire peak by {D1c.loc[('peak', 'Stable'), 'median']:.2f} (Stable) and {D1c.loc[('peak', 'Dynamic'), 'median']:.2f} (Dynamic) log2 units, whereas the same-size shuffles 300&nbsp;bp away change it by about {D1c.loc[('peak', 'Stable'), 'off']:.2f}. So the 170&nbsp;bp around the motif accounts for a large part of the accessibility of the whole peak. <strong>Trust:</strong> robust (large effects, paired design, consistent in both classes).</p>

<h3>D2. Mutation arm: what carries the dependence?</h3>
{dg8}
<figure>{oimg('d2_mutation_arm', 'Mutation arm dependency')}<figcaption>Dependency (&minus;log2FC) for each edit, read out over the region it edits; Stable vs Dynamic with 95% CI. Top: cardiomyocyte sites (heart ATAC); bottom: fibroblast sites (Patski ATAC). Edit labels give the number of bases changed (the peak width varies; the median is {PEAK_MEDIAN} bp).</figcaption></figure>
<figure>{oimg('d2b_what_carries_dependence', 'What carries the dependence', '860px')}<figcaption>Stable sites. Left: within the motif, the 7&nbsp;bp core carries the whole effect of mutating all 11 bases. Right: mutating the 170&nbsp;bp while leaving the motif intact reproduces most of the loss from mutating everything.</figcaption></figure>
{table(mut_tab, {c: '{:.2f}' for c in ['CM Stable', 'CM Dynamic', 'FB Stable', 'FB Dynamic']})}
<p><strong>The surroundings, not the motif, carry the dependence.</strong> Mutating the 170&nbsp;bp but keeping the motif intact costs {ratio_cm170:.0%} of the loss from mutating everything in cardiomyocyte sites ({dec('170', 'flank@170', 'CM, heart ATAC'):.2f} vs {dec('170', 'random@170', 'CM, heart ATAC'):.2f}) and {ratio_fb170:.0%} in fibroblast sites. Within the motif, the 7&nbsp;bp core accounts for essentially all of it ({ratio_cmcore:.0%} CM, {ratio_fbcore:.0%} FB). At the 170&nbsp;bp and peak widths every edit is significantly more damaging in Stable than Dynamic sites (FDR &lt; 0.001), in both site sets; at motif width the classes barely differ ({', '.join(f'{EDIT_LABEL[e].lower()} {stars(depq(CMS, e))}' for e in ('random@motif', 'dinuc@motif', 'core@motif', 'flank@motif'))} in cardiomyocyte sites; all ns in fibroblast sites).</p>
<p><strong>Caveats.</strong> The 170&nbsp;bp readout averages over 170&nbsp;bp of which 11 are the motif, so a ratio near 100% is expected even if the motif matters locally; the result is also track-dependent (91% CM, 76% FB). A dinucleotide-preserving shuffle of an 11-mer has very few distinct arrangements (3.8% of its replicates at the motif width gave a native-identical prediction), so the dinucleotide shuffle of the motif is partly a no-op and understates its effect.</p>
<p><strong>Design note.</strong> The core edit replaces the 7 core bases with random different bases, redrawn until no AP-1 half-site is left (3 draws per site); it is not a targeted substitution, so its size depends on that design. <strong>Trust:</strong> the decomposition and the class pattern are robust.</p>

<h3>D3. Fibroblast sites: does the Stable-vs-Dynamic result replicate?</h3>
<figure>{oimg('d3a_fb_heatmaps', 'Fibroblast shuffle heatmaps')}<figcaption>Fibroblast sites (5,018): median shuffle effect on each readout with the Patski ATAC track. Marks on the Dynamic panel: Mann&ndash;Whitney vs Stable, BH over the 9 cells.</figcaption></figure>
<figure>{oimg('d3b_replication_forest', 'Replication across site sets and tracks')}<figcaption>Dynamic minus Stable (median effect, 95% CI) for cardiomyocyte and fibroblast sites, each with the heart and the Patski track. Positive = Dynamic loses less.</figcaption></figure>
{table(rep_tab, {c: '{:+.2f}' for c in ['Raw Dyn - Stable', 'Raw CI lo', 'Raw CI hi', 'Adjusted', 'Adj. CI lo', 'Adj. CI hi']})}
<p><strong>Yes, in direction.</strong> At the 170&nbsp;bp shuffle Stable sites lose more than Dynamic sites in all four combinations (raw FDR &lt; 0.001): CM sites {rep.loc[('CM sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (heart) and {rep.loc[('CM sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (Patski), FB sites {rep.loc[('FB sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} and {rep.loc[('FB sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f}. As in cardiomyocyte sites, the motif-width shuffle does not separate the classes, and the peak-width difference is smaller. <strong>The size is track-dependent, though.</strong> After adjusting for native accessibility and peak width, the 170&nbsp;bp gap stays larger with the heart track (the matched-readout and further adjustments in Part H shrink it much more) ({adj_.loc[('CM sites, heart ATAC', '170'), 'adj_diff']:+.2f} CM sites, {adj_.loc[('FB sites, heart ATAC', '170'), 'adj_diff']:+.2f} FB sites) but shrinks to {adj_.loc[('CM sites, Patski ATAC', '170'), 'adj_diff']:+.2f} and {adj_.loc[('FB sites, Patski ATAC', '170'), 'adj_diff']:+.2f} with the Patski track (still FDR &lt; 0.01); the peak-width gap disappears after adjustment in 3 of 4 combinations. So with the fibroblast-like track much of the raw difference reflects accessibility and peak width, with a small residual.</p>
<figure>{oimg('d3c_fb_class_effects_170', 'Class effects in fibroblast sites, 170 bp shuffle')}<figcaption>Fibroblast sites, Patski track, 170 bp shuffle read over the 170 bp: median effect per class (95% CI) for the age comparison, the exercise comparison and all seven binding patterns. Blue = reference class; marks = Mann&ndash;Whitney against the reference, BH within each panel.</figcaption></figure>
<details><summary>The same plot for the motif shuffle and the peak shuffle</summary>
<figure>{oimg('d3c_fb_class_effects_motif', 'Class effects in fibroblast sites, motif shuffle')}<figcaption>Motif shuffle, effect read over the motif.</figcaption></figure>
<figure>{oimg('d3c_fb_class_effects_peak', 'Class effects in fibroblast sites, peak shuffle')}<figcaption>Peak shuffle, effect read over the whole peak.</figcaption></figure></details>
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
H_tab["p"] = H_tab.p.map(lambda v: "< 0.001" if v < 0.001 else f"{v:.3f}")
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
<p><strong>What was wrong.</strong> The earlier adjustment in the Stable-vs-Dynamic checks (B3) controlled for native accessibility over the <em>peak</em> for every readout, and the stratification used quintiles of the same variable; the class comparisons (B4) have since been redone with the matched readout. The effect of a 170&nbsp;bp shuffle depends on how accessible those 170&nbsp;bp are, so the correct covariate is the native value at the matched readout. With that variable the gap shrinks a lot. Models below use standard errors clustered by TOBIAS peak, 170&nbsp;bp shuffle.</p>
<figure>{img(AU / 'figures' / 'h1_gap_models.png', 'Gap by model')}<figcaption>Dynamic minus Stable (log2FC, 95% CI) as covariates are added. Positive = Dynamic loses less. The last row includes footprint score, which partly defines the classes, so it is a sensitivity check and not a causal adjustment.</figcaption></figure>
{table(H_tab, HF)}
<figure>{img(AU / 'figures' / 'h2_gap_by_matched_decile.png', 'Gap by decile of matched native accessibility')}<figcaption>Top: median 170&nbsp;bp shuffle effect by decile of native accessibility at the 170&nbsp;bp readout. Bottom: number of Stable (blue) and Dynamic (orange) sites per decile; most Dynamic sites sit in the lower deciles.</figcaption></figure>
<p><strong>Reading.</strong> In cardiomyocyte sites the raw gap of {hg(cmn, M_RAW).gap:+.2f} becomes {hg(cmn, M_PK).gap:+.2f} with the earlier adjustment, {hg(cmn, M_MATCH).gap:+.2f} with the matched-readout value, {hg(cmn, M_FULL).gap:+.2f} after also adjusting for peak accessibility, peak width and TSS distance with flexible terms, and {hg(cmn, M_FP).gap:+.2f} (CI {hg(cmn, M_FP).ci_lo:+.2f} to {hg(cmn, M_FP).ci_hi:+.2f}) with footprint score. With the fibroblast-like Patski track the gap is gone after the matched adjustment ({hg(fbp, M_MATCH).gap:+.2f}). Fibroblast sites with the heart track keep a small gap ({hg(fbh, M_FULL).gap:+.2f}, CI {hg(fbh, M_FULL).ci_lo:+.2f} to {hg(fbh, M_FULL).ci_hi:+.2f}). In the decile plot the gap is about 0.1 to 0.2 across the lower deciles, where {H_cs.loc[cmn, 'dynamic_in_lowest_5_deciles']:.0%} of Dynamic cardiomyocyte sites lie, and larger only in the top deciles where there are few Dynamic sites. Weighted to the Dynamic distribution, the gap for sites like Dynamic ones is {H_cs.loc[cmn, 'common_support_gap']:+.2f} (CM sites), {H_cs.loc[fbh, 'common_support_gap']:+.2f} (FB, heart) and {H_cs.loc[fbp, 'common_support_gap']:+.2f} (FB, Patski).</p>
{table(H_oth, HF)}
<p><strong>What this means.</strong> The raw difference is real, but it is mostly a difference in chromatin context. Stable sites lose more because they start from more accessible chromatin with stronger footprints, wider peaks and, for promoter-proximal sites, a promoter. When those are taken into account the gap shrinks to {hg(cmn, M_FULL).gap:+.2f} log2 in cardiomyocyte sites and {hg(fbh, M_FULL).gap:+.2f} in fibroblast sites (heart track), and disappears with the Patski track. A small residual remains in the heart-track analyses and adding footprint score removes it in cardiomyocyte sites. The motif-width and peak-width differences disappear after adjustment everywhere.</p>
<div class="bio"><strong>What this means biologically.</strong> AP-1 sites whose binding changes with age or exercise do not appear to be built differently from constitutive AP-1 sites: the model finds no extra, separate sequence dependence for them once their weaker chromatin context is accounted for. That points away from a distinct sequence grammar for responsive AP-1 sites and towards their responsiveness being set by something that a sequence-only baseline model does not see, such as changes in AP-1 activity or abundance, co-factors or chromatin state in response to age and exercise. This agrees with the manuscript's picture that dynamic elements are marked by weaker baseline chromatin context rather than by a different core motif. It cannot be tested here, because the model predicts baseline accessibility and not how accessibility responds to a stimulus. Patski is not an independent confirmation, since it scores the same sequences with the same model.</div>
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
<li><strong>The dinucleotide shuffle of the motif is partly a no-op</strong>: an 11-mer has few dinucleotide-preserving arrangements and 3.8% of its replicates gave a native-identical prediction, which dilutes its effect (median &minus;0.33 vs &minus;0.67 for replacing the motif). Fix: redraw until different, or report it as limited.</li>
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
     f"(smallest FDR {E[E.width != 'motif'].fdr.min():.2f}). No difference detected, but the intervals are too wide (about +-0.2 to 0.3) to exclude differences half the size of the Stable-vs-Dynamic gap: absence of evidence."),
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
    ("Promoter-proximal sites depend less on their 170 bp context than distal sites", "mid",
     f"Median 170 bp effect in Stable sites: {an6.loc[('eff_170', 'Stable', '<=2 kb'), 'median']:+.2f} within 2 kb of a TSS vs {an6.loc[('eff_170', 'Stable', '>30 kb'), 'median']:+.2f} beyond 30 kb (Mann-Whitney p {fmtp(_mwu(SL[SL.tss_bin == '<=2 kb'].eff_170.dropna(), SL[SL.tss_bin == '>30 kb'].eff_170.dropna()).pvalue) if 'tss_bin' in SL else 'n/a'}); promoter-proximal peaks are also about 3x more accessible and the comparison is not adjusted for accessibility, so part of this may be accessibility itself."),
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
    ("The loss after a 170 bp shuffle is local to the site's neighbourhood, not a generic effect of shuffling 170 bp of the peak (a distance control, not an AP-1-specific control)", "ok",
     f"Centred shuffle {co_.loc['Stable', 'centred_170_median']:.2f} vs {co_.loc['Stable', 'offset_pm300_median']:.2f} for the same-size region 300 bp away (Stable); the centred shuffle is worse at {co_.loc['Stable', 'share_centred_more_negative']:.0%} (Stable) and {co_.loc['Dynamic', 'share_centred_more_negative']:.0%} (Dynamic) of sites, paired Wilcoxon p &lt; 1e-10."),
    ("The effect grows smoothly with the shuffled width and saturates near 340 bp; the 7 bp core alone gives about a quarter", "ok",
     f"Share of the 170 bp effect: 7 bp {fr_.loc[7]:.0%}, 11 bp {fr_.loc[11]:.0%}, 85 bp {fr_.loc[85]:.0%}, 340 bp {fr_.loc[340]:.0%} (Stable, motif readout)."),
    ("The surroundings, not the motif, carry the sequence dependence (mutation arm)", "mid",
     f"170 bp mutated with the motif intact = {ratio_cm170:.0%} (CM) and {ratio_fb170:.0%} (FB) of the loss from mutating everything (a ratio near 100% is partly expected because the readout averages over 170 bp); within the motif the 7 bp core = {ratio_cmcore:.0%} (CM) and {ratio_fbcore:.0%} (FB) of the 11 bp loss."),
    ("Stable sites lose more than Dynamic sites at 170 bp also in fibroblast sites and with the fibroblast-like (Patski) track", "mid",
     f"Raw difference at 170 bp: CM sites {rep.loc[('CM sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} / FB sites {rep.loc[('FB sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} with Patski, {rep.loc[('FB sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (FB sites, heart); all FDR &lt; 0.001. After adjusting for native accessibility at the matched readout the Patski gap is {hg(fbp, M_MATCH).gap:+.2f} (ns), so the replication is of the raw difference only."),
    ("Placing MPRA elements in a neutral genomic context makes the native accessibility head agree with the MPRA changes", "bad",
     f"No change: shuffle-effect correlations stay at about 0 (|rho| &le; {D4[D4.comparison != 'native level vs MPRA activity'].rho.abs().max():.2f}); native levels correlate {mc_.loc[('heart ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (heart) and {mc_.loc[('Patski ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (Patski) in the neutral context."),
]

claims_html = "".join(f"<tr><td class='v'>{c}</td><td>{badge(k)}</td><td class='v'>{ev}</td></tr>" for c, k, ev in claims)

PART_B = f"""
<h2 id="partB">Part B. Main results</h2>
<h3 class="sec">B1. Does the model track real accessibility?</h3>
<figure>{ximg('x4_observed_vs_predicted', 'Observed vs predicted accessibility')}<figcaption>Predicted native accessibility over the peak against measured snATAC accessibility at the same peak, as a 2-D density of sites. Left: cardiomyocyte sites with the heart track. Middle: fibroblast sites with the same heart track. Right: the same fibroblast sites with the Patski track. The three panels use the same sites and readout so the tracks can be compared directly.</figcaption></figure>
{table(XOV_tab, {'Spearman rho': '{:.2f}', 'CI lo': '{:.2f}', 'CI hi': '{:.2f}'})}
{table(NV_tab, {'Spearman rho': '{:.2f}'})}
<p><strong>What it shows.</strong> The heart track predicts measured cardiomyocyte accessibility with a moderate rank correlation (rho {xov.loc['CM sites, heart ATAC', 'rho']:.2f}) and fibroblast-site accessibility almost as well ({xov.loc['FB sites, heart ATAC', 'rho']:.2f}). The Patski (fibroblast-like) track does much worse on the same fibroblast sites ({xov.loc['FB sites, Patski ATAC', 'rho']:.2f}) and is essentially uncorrelated with cardiomyocyte accessibility, as it should be. <strong>Trust:</strong> good enough to use the heart track as a proxy for cardiomyocyte sites and, with care, for fibroblast sites; it is bulk heart tissue, not purified cell types.</p>

<h3 class="sec">B2. What does shuffling do?</h3>
<figure>{img(CM / 'figures' / 'fig2s_shuffle_heatmaps_stars.png', 'Shuffle effect heatmaps')}
<figcaption>Median log2FC for each shuffle width (rows) and ATAC readout region (columns). Left: Stable sites; right: Dynamic sites, with significance marks against Stable.</figcaption></figure>
<figure>{ximg('x2_shuffle_dotplot', 'Shuffle effects as a grouped dot plot')}<figcaption>The same numbers as a grouped dot plot, so Stable and Dynamic sites can be compared side by side. Each panel is one shuffle; the three positions are the region over which the effect is read. Median with 95% bootstrap CI; marks: Stable vs Dynamic (Mann&ndash;Whitney, BH over the nine comparisons).</figcaption></figure>
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
<figure>{img(CM / 'figures' / 'fig3s_class_effects_170.png', 'Class effects at the 170 bp shuffle')}<figcaption>The 170 bp shuffle, effect read over the 170 bp: median effect per class with 95% bootstrap CI, for the age comparison (young vs aged sedentary), the exercise comparison (aged sedentary vs aged exercised) and all seven binding patterns. Blue = reference class; marks: Mann&ndash;Whitney against the reference, BH within each panel. The reference classes (bound in both conditions) include a few Dynamic sites, so they are not identical to the Stable set used elsewhere.</figcaption></figure>
<details><summary>The same plot for the motif shuffle and the peak shuffle</summary>
<figure>{img(CM / 'figures' / 'fig3s_class_effects_motif.png', 'Class effects at the motif shuffle')}<figcaption>Motif shuffle, effect read over the motif.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig3s_class_effects_peak.png', 'Class effects at the peak shuffle')}<figcaption>Peak shuffle, effect read over the whole peak.</figcaption></figure></details>
<figure>{img(CM / 'figures' / 'fig7_class_effects_adjusted_170.png', 'Adjusted class effects, 170 bp shuffle')}<figcaption>Difference from the reference class for the 170 bp shuffle: raw (orange diamonds) and after adjusting for native accessibility at the same readout and for peak width (blue, 95% CI). Positive = the class loses less than the reference.</figcaption></figure>
<details><summary>The same plot for the motif shuffle and the peak shuffle</summary>
<figure>{img(CM / 'figures' / 'fig7_class_effects_adjusted_motif.png', 'Adjusted class effects, motif shuffle')}<figcaption>Motif shuffle, effect read over the motif.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig7_class_effects_adjusted_peak.png', 'Adjusted class effects, peak shuffle')}<figcaption>Peak shuffle, effect read over the whole peak.</figcaption></figure></details>
<p><strong>What it shows.</strong> At 170&nbsp;bp every dynamic class loses about {a170.raw_diff.min():.2f} to {a170.raw_diff.max():.2f} log2FC less than the reference in the raw comparison and {a170.adj_diff.min():.2f} to {a170.adj_diff.max():.2f} after adjusting for native accessibility at the same readout and peak width (all {int((a170.adj_fdr < 0.05).sum())} age and exercise classes remain FDR &lt; 0.05). At the peak and motif widths the raw differences are small and the adjusted ones shrink toward zero ({int((adj[adj.width == 'peak'].adj_fdr < 0.05).sum())} of 6 classes at the peak width, {int((adj[adj.width == 'motif'].adj_fdr < 0.05).sum())} of 6 at the motif width remain significant). <strong>Between dynamic classes</strong> (the direct test of age vs exercise, and of up vs down):</p>
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
<p><strong>What it shows.</strong> Predicted RNA barely moves when a distal AP-1 site or peak is shuffled (medians within about &plusmn;0.1), with the largest median drops 2 to 10&nbsp;kb from the TSS for the &plusmn;2&nbsp;kb readout and within 2&nbsp;kb for the &plusmn;50&nbsp;kb mean (peak shuffle, Stable sites: about &minus;0.08 and &minus;0.12). ATAC and RNA effects correlate across sites (rho 0.19 to 0.24), so there is a weak coupling, but the RNA effects themselves are tiny. <strong>Trust:</strong> not informative with these tracks; do not over-interpret.</p>

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
<h4>Stable versus Dynamic donors</h4>
<figure>{ximg('x3_peak_swap_stable_vs_dynamic', 'Peak swap, Stable vs Dynamic donors')}<figcaption>Peak swap by donor class (48 Stable and 25 Dynamic donor peaks, each placed into the 4 recipient windows). Left: the donor peak's own native accessibility. Middle: accessibility after the swap against the donor's own accessibility, with a fit per class. Right: how much of the donor's accessibility the peak keeps in each recipient (log2 swapped / donor native), median with 95% CI.</figcaption></figure>
{table(XPS_tab, {'Stable': '{:.3g}', 'Dynamic': '{:.3g}', 'p': '{:.3g}'})}
<p>Dynamic donors start from lower accessibility (median {XPS.iloc[0].dynamic_median:.2f} vs {XPS.iloc[0].stable_median:.2f}, p = {XPS.iloc[0].p_mannwhitney:.2f}), but once the donor's own accessibility is taken into account the two classes behave the same: the coefficient on being a Dynamic donor is {XPS.iloc[-1].stable_median:+.2f} log2 (95% CI {XPS.iloc[-1].dynamic_median:+.2f} to {XPS.iloc[-1].n_dynamic_donors:+.2f}, p = {XPS.iloc[-1].p_mannwhitney:.2f}) and the swapped-vs-donor correlations are {XPS.iloc[4].stable_median:.2f} (Stable) and {XPS.iloc[5].stable_median:.2f} (Dynamic). There is no detectable difference, but with only 25 Dynamic donors a difference smaller than about 0.2 log2 could be missed. The right panel shows that location matters for the level: the same peak keeps less than its own accessibility in the least accessible recipient and more in the more accessible ones, for both classes alike.</p>
<figure>{img(AN / 'insertion' / 'figures' / 'fig_motif_insertion.png', 'Motif insertion')}
<figcaption>A TGAGTCA motif inserted at the centre of dinucleotide-shuffled random genomic backgrounds (50 backgrounds), compared with the shuffled background and with a shuffled-motif control.</figcaption></figure>
{table(mi_tab, {'Median log2FC': '{:.2f}', 'Fraction > 0': '{:.2f}'})}
<p><strong>Trust:</strong> no usable signal. Background accessibility is near zero, so any ratio is dominated by the floor, and the 7&nbsp;bp readout differences are not significant.</p>

"""

# ======================= Assembly: chapters ordered by the scientific story =======================
import re


def split_part(text, pat):
    ms = list(re.finditer(pat, text))
    intro = text[: ms[0].start()] if ms else text
    out = {}
    for i, mm in enumerate(ms):
        end = ms[i + 1].start() if i + 1 < len(ms) else len(text)
        out[mm.group(1)] = text[mm.end(): end]
    return intro, out


SEC = {}
for part, pat in ((PART_A, r'<h3>(A\d)\. .*?</h3>'), (PART_B, r'<h3 class="sec">(B\d)\. .*?</h3>'), (PART_C, r'<h3>(C\d)\. .*?</h3>'), (PART_D, r'<h3>(D\d)\. .*?</h3>'),
                  (PART_G, r'<h3>(G\d)\. .*?</h3>'), (PART_H, r'<h3>(H\d)\. .*?</h3>')):
    intro_, d_ = split_part(part, pat)
    SEC.update(d_)
    if part is PART_G:
        SEC["G0"] = re.sub(r'<h2.*?</h2>', '', intro_, flags=re.S)
# the B8 section holds both insertion experiments: split them
_b8 = SEC.pop("B8")
_i = _b8.index('alt="Motif insertion"')
_j = _b8.rfind("<figure>", 0, _i)
SEC["B8a"], SEC["B8b"] = _b8[:_j], _b8[_j:]
# drop the old per-part summaries and the final G7 list
SEC.pop("G7", None)
SEC.pop("D5", None)


def conf(kind, text):
    return f'<p class="conf">{badge(kind)} {text}</p>'


CONF = {
 "A1": ("ok", "The site table was re-derived independently and reproduces."),
 "A2": ("ok", "Descriptive; the numbers come straight from the site table."),
 "A4": ("mid", "Gains are clear, losses are weak; footprint calls and accessibility come from the same ATAC data, so this is not fully independent."),
 "A3": ("ok", "Descriptive; Dynamic and Stable differences agree with the measured accessibility."),
 "B1": ("mid", "Rho about 0.4 is informative but not a precise predictor; the heart track is bulk tissue and the Patski track is a weak fibroblast proxy."),
 "A5": ("mid", "Tested between classes: errors are similar in size, but Dynamic sites are predicted somewhat less well and are biased upwards."),
 "A6": ("ok", "The shuffle behaves as designed."),
 "A7": ("ok", "Measured against the 10-replicate value itself, so true error is somewhat larger."),
 "B7": ("ok", "Medians are precise; single-site peak-shuffle values are noisy."),
 "A8": ("ok", "Sites are interleaved across chunks, so a GPU effect would show as a shift between chunks; none is seen."),
 "B2": ("ok", "A description of what the model predicts, supported by every later control. It is a model prediction, not a measurement."),
 "C1": ("ok", "Descriptive."),
 "D1": ("mid", "High confidence that the loss is local and saturates near 340 bp. Moderate that this reflects AP-1 specifically: the offset shuffles are only a distance control."),
 "D2": ("mid", "The 170 bp readout averages over 170 bp of which 11 are the motif, so a ratio near 100% is partly expected; the result differs between tracks (91% vs 76%), and the dinucleotide shuffle of the motif is partly a no-op."),
 "B8a": ("mid", "Clear direction (sequence beats location) and no detectable Stable/Dynamic difference, but only 4 recipient windows and 25 Dynamic donors, so the pairs are not independent and power is limited."),
 "C3": ("mid", "The 50% thresholds are arbitrary; read the proportions descriptively."),
 "C2": ("ok", "Large, consistent correlations."),
 "C6": ("mid", "Not adjusted for accessibility, so part of the promoter effect may be accessibility itself."),
 "C4": ("mid", "Footprint scores are computed from the same ATAC data that define accessibility."),
 "C5": ("bad", "Few sites in multi-motif peaks, and each shuffle edits one window, not all motifs together; this is no evidence either way."),
 "B3": ("ok", "The raw difference is real as a description. Whether it is more than accessibility is addressed in the next section."),
 "H1": ("mid", "High confidence that accessibility and context explain most of the gap. The footprint-score model is a sensitivity check, since footprint score partly defines the classes."),
 "B4": ("bad", "Underpowered for between-class differences (intervals about +-0.2 to 0.3, 79 to 210 sites) and the labels overlap."),
 "D3": ("mid", "Replicates the raw difference only. The fibroblast tracks are weak proxies and Patski scores the same sequences with the same model, so it is not an independent confirmation."),
 "B6": ("mid", "Selection uses only native prediction and measurement, never the effects."),
 "G0": ("ok", "This is the definition of the restriction; the rest of the chapter shows what it keeps and what it changes."),
 "G1": ("ok", "Descriptive."),
 "G2": ("mid", "The filter removes mostly over-predicted Dynamic sites, so the filtered Dynamic group is the well-predicted subset."),
 "G3": ("mid", "Selection and effect are related for Stable sites, so filtered gaps are partly selection."),
 "G4": ("mid", "Same direction in every tier, but see the previous section for why the larger filtered gaps are partly selection."),
 "G5": ("ok", "Direction and significance do not depend on the floor; the size does somewhat."),
 "G6": ("ok", "Unchanged by the filter."),
 "C7": ("mid", "The native-level correlation is weak but consistent. The shuffle-effect correlations are about 0, but the MPRA's own scramble effects are tiny, so this null does not discriminate model quality."),
 "D4": ("mid", "A clean null, but it rests on only 2 neutral background windows and the same small MPRA effects."),
 "B5": ("bad", "Median RNA effects are about zero and per-site values are noisy; not informative with these tracks."),
 "B8b": ("bad", "Background accessibility is near zero, so any ratio is dominated by the floor, and the control can keep AP-1 half-sites."),
}

# chapters: key, title, intro (question, approach, finding), sub-sections [(key, title)]
N = {}
CH = [
 ("design", "The experiments at a glance", None, []),
 ("found", "Can we trust the data and the model?",
  "<p><strong>Question.</strong> Before interpreting any perturbation, are the sites, labels and predictions sound enough to use? <strong>Approach.</strong> Inspect the site set and labels, compare predicted with measured accessibility, check that shuffles do what they should, and measure replicate noise. <strong>Finding.</strong> The sites and labels are sound, the model predicts measured accessibility moderately (rho {nv.rho:.2f} in cardiomyocyte sites), shuffles behave as designed, and 10 replicates are adequate. <strong>Limits.</strong> Prediction is moderate, not precise, Dynamic sites are somewhat less well predicted and are over-predicted, and the fibroblast tracks are weak proxies.</p>",
  [("A1", "Which sites were analysed"), ("B1", "Does the model predict measured accessibility? Cardiomyocyte and fibroblast sites, heart and Patski tracks"), ("A2", "Stable and Dynamic sites sit in different peaks"),
   ("A4", "Are the Stable and Dynamic labels real? A check against measured snATAC"), ("A3", "Native predicted accessibility"), ("A5", "Is prediction quality different between site classes?"),
   ("A6", "Does a shuffle lower accessibility as expected?"), ("A7", "How many shuffle replicates are enough?"), ("B7", "How noisy is a single site?"), ("A8", "Does the GPU model matter?")]),
 ("depend", "What do AP-1 sites depend on in the model? The surrounding sequence, not the motif",
  "<p><strong>Question.</strong> When the model predicts accessibility at an AP-1 site, how much comes from the 11&nbsp;bp motif and how much from its surroundings? <strong>Approach.</strong> Shuffle and mutate the motif, the 170&nbsp;bp around it or the whole peak, with controls at other positions and a width titration, and swap whole peaks. <strong>Finding.</strong> The 170&nbsp;bp loses far more accessibility ({m.loc['170', 'median_stable']:.2f} log2) than the motif alone ({m.loc['motif', 'median_stable']:.2f}); the loss is local, builds with the width of the edit and saturates near 340&nbsp;bp; mutating the 170&nbsp;bp with the motif left intact reproduces {ratio_cm170:.0%} of the full loss in cardiomyocyte sites ({ratio_fb170:.0%} in fibroblast sites); and sequence rather than location sets accessibility in a peak swap. <strong>Limits.</strong> These are model predictions; the offset control is a distance control only; the motif is 11 of the 170&nbsp;bp read out.</p>",
  [("B2", "What shuffling does at three widths"), ("C1", "How varied are the per-site effects?"), ("D1", "Controls and width titration: is the loss local, and how does it scale?"), ("D2", "Mutation arm: what carries the dependence?"),
   ("B8a", "Peak swap: sequence or location, and Stable versus Dynamic donors"), ("C3", "Motif dependence versus context dependence"), ("C2", "What sets the size of the effect?"), ("C6", "Promoter-proximal versus distal sites"),
   ("C4", "Does footprint strength agree with the predicted dependence?"), ("C5", "Do two Fos sites in a peak buffer each other?")]),
 ("dyn", "Do dynamic sites behave differently from stable ones?",
  "<p><strong>Question.</strong> Do sites whose binding changes with age or exercise depend on their sequence differently from sites bound in every condition? <strong>Approach.</strong> Compare Stable and Dynamic sites at each shuffle width, raw and after adjusting for native accessibility, peak width, TSS distance and footprint strength; compare age and exercise classes; repeat in fibroblast sites. <strong>Finding.</strong> Stable sites lose more at the 170&nbsp;bp shuffle ({m.loc['170', 'diff_dyn_minus_stable']:+.2f} raw), but most of the gap is explained by native accessibility and context: {hg(cmn, M_MATCH).gap:+.2f} with the matched-readout adjustment, {hg(cmn, M_FULL).gap:+.2f} with peak width and TSS distance, and {hg(fbp, M_MATCH).gap:+.2f} with the Patski track. Age and exercise classes cannot be told apart at these sample sizes. <strong>Limits.</strong> The earlier, simpler adjustment looked as though the gap was robust; the matched-readout adjustment shows it was not.</p>",
  [("B3", "Stable versus Dynamic: the raw difference"), ("H1", "Correcting for accessibility at the matched readout"), ("B4", "Age, exercise and their combination"), ("D3", "Does it replicate in fibroblast sites?")]),
 ("robust", "Robustness: restricting to sites where prediction matches measurement",
  "<p><strong>Question.</strong> Are the results driven by noisy or poorly predicted regions? <strong>Approach.</strong> Keep only sites whose predicted accessibility agrees with the measured snATAC, chosen from the native prediction and measurement alone, and repeat the main analyses. <strong>Finding.</strong> The Stable-vs-Dynamic difference, the mutation decomposition and the controls hold in every tier. <strong>Limits.</strong> The filter removes mostly over-predicted Dynamic sites, and kept Stable sites have larger effects than dropped ones, so larger gaps in the filtered sets are partly selection.</p>",
  [("B6", "Confidence tiers from the main analysis"), ("G0", "How the restriction was done"), ("G1", "What was kept and what was dropped"), ("G2", "Is the filter biased against one class?"),
   ("G3", "Do kept and dropped sites differ in the effect?"), ("G4", "Does Stable versus Dynamic survive?"), ("G5", "Sensitivity to the floor"), ("G6", "Mutation arm and controls in the filtered sets")]),
 ("mpra", "Does the model explain the fibroblast MPRA?",
  "<p><strong>Question.</strong> Does the predicted accessibility change after scrambling an element track the change in activity in the lenti-STARR MPRA? <strong>Approach.</strong> Predict native, motif-scrambled and fully scrambled versions of 1,653 elements in their genomic window and in neutral windows, and correlate with the MPRA. <strong>Finding.</strong> Native predicted accessibility correlates weakly with MPRA activity (rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski), while predicted scrambling effects do not track the MPRA at all, in either context. <strong>Limits.</strong> The MPRA's own scramble effects are tiny (mean core loss about 0.02 log2), so there is little variance to predict.</p>",
  [("C7", "Predicted accessibility against the MPRA"), ("D4", "Does placing elements in a neutral context restore agreement?"), ("MAN", "What the manuscript claims, and how our results bear on it")]),
 ("fail", "Experiments that did not work, or are limited",
  "<p>These are reported so the reader can see what was tried and why it should not be used.</p>",
  [("B5", "Predicted RNA"), ("B8b", "Motif insertion into random backgrounds"), ("NOTE", "Other limits and unresolved points")]),
 ("reg", "Experiment register: every experiment, with confidence and reasons", None, []),
 ("bio", "What it means biologically", None, []),
 ("open", "Review log, open issues and next steps", None, [("H2", "Corrections made after the independent review"), ("H3", "Code issues found and not yet fixed"), ("NEXT", "Next steps")]),
]
for ci, (ck, ct, cintro, subs) in enumerate(CH, 1):
    N[ck] = str(ci)
    for si, (sk, st) in enumerate(subs, 1):
        N[sk] = f"{ci}.{si}"

REFS = [(r"Parts B3 and B4", f"sections {N['B3']} and {N['B4']}"), (r"B3 and B4", f"sections {N['B3']} and {N['B4']}"), (r"(?<![\w-])Part H\b", f"section {N['H1']}"), (r"(?<![\w-])Part G\b", f"section {N['G0']}"),
        (r"(?<![\w-])Part C\b", f"section {N['C2']}"), (r"(?<![\w#-])B6\b", f"section {N['B6']}"), (r"(?<![\w#-])B7\b", f"section {N['B7']}"), (r"(?<![\w#-])G3\b", f"section {N['G3']}"), (r"(?<![\w#-])D1\b", f"section {N['D1']}"),
        (r"(?<![\w#-])D2\b", f"section {N['D2']}"), (r"(?<![\w#-])D4\b", f"section {N['D4']}"), (r"(?<![\w#-])A4\b", f"section {N['A4']}"), (r"(?<![\w#-])C4\b", f"section {N['C4']}")]
ANCH = {"B3": "B3", "H1": "H1", "G0": "G0", "C2": "C2", "B6": "B6", "B7": "B7", "G3": "G3", "D1": "D1", "D2": "D2", "D4": "D4", "A4": "A4", "C4": "C4"}


def fixrefs(s):
    for pat, rep_ in REFS:
        def sub(mm, rep_=rep_):
            num = re.search(r"section[s]? ([\d.]+)", rep_).group(1)
            key = next(k for k, v in N.items() if v == num and k in SEC)
            return rep_.replace(f"section {num}", f'<a href="#sec-{key}">section {num}</a>', 1) if not rep_.startswith("sections") else \
                f'sections <a href="#sec-B3">{N["B3"]}</a> and <a href="#sec-B4">{N["B4"]}</a>'
        s = re.sub(pat, sub, s)
    return s


# ---- hand-written sections
MAN_ROWS = [
 ("Dynamic AP-1 elements have lower baseline accessibility", "ok", "Predicted native ATAC is lower for Dynamic MPRA elements (Cohen's d about -0.5) and for Dynamic cardiomyocyte sites."),
 ("Dynamic elements have lower GC content and sit farther from promoters", "mid", "GC lower (median 0.494 vs 0.512, p = 9e-8); distance only weakly different in the MPRA elements and not significant in our cardiomyocyte sites."),
 ("Age-gained elements have reporter activity on par with stable elements", "ok", "Median log2 alpha +0.25 (Age Up) vs +0.26 (Stable), p = 0.78."),
 ("...despite weaker endogenous chromatin features: chromatin context, not sequence autonomy, explains them", "ok", "The model predicts lower accessibility for the same sequences while MPRA activity is equal; native accessibility explains little of MPRA activity."),
 ("Shuffling the AP-1 core has a minor effect on MPRA activity; flank shuffling hurts dynamic elements more", "ok", "Core loss about 0.02 log2; flank loss reproduces (Age Up 0.128 vs Stable 0.016, d = 0.14), a small effect."),
 ("The same flank dependence is a property of what the model can read", "bad", "The native accessibility head gives the opposite direction (Dynamic lose less), and per-element predicted vs measured flank loss correlates about 0. Accessibility and episomal activity are different quantities."),
 ("AlphaGenome-MPRA probe reaches r = 0.73; ISM, seqlets, epistasis, rescue explain dispersed activity", "na", "Cannot be tested with the native heads. GC content alone correlates 0.505 with MPRA activity, so a GC baseline on the same splits would help interpret 0.73."),
 ("Motif enrichments (MEF2C, TCF7L2, GATA4, NFAT, ETV4, C/EBP)", "na", "Needs a motif scan on sequences, which can run on Katana."),
]
MAN = f"""<p>The manuscript uses AlphaGenome as a frozen encoder with a linear probe trained on the MPRA, not the native accessibility head used here, so our results speak to predicted chromatin accessibility and not to the manuscript's probe. Full notes: <code>Analysis/manuscript_notes_MPRA_AlphaGenome.md</code>.</p>
<table class="t reg"><tr><th>Manuscript claim</th><th>Verdict from our data</th><th>Evidence</th></tr>{''.join(f"<tr><td class='v'>{a}</td><td>{badge(k).replace('High confidence', 'Supported').replace('Moderate confidence', 'Partly supported').replace('Low / unreliable', 'Opposed').replace('Untested', 'Cannot test')}</td><td class='v'>{b}</td></tr>" for a, k, b in MAN_ROWS)}</table>
<p class="small">Supported = our data reproduce or support the claim; Opposed = our data point the other way; Cannot test = outside what the native heads can address.</p>"""
NOTE = f"""<ul>

<li><strong>The dinucleotide shuffle of the motif is partly a no-op</strong>: an 11-mer has few distinct dinucleotide-preserving arrangements and 3.8% of its replicates gave a native-identical prediction, so its effect (median &minus;0.33) is diluted compared with replacing the motif (&minus;0.67).</li>
<li><strong>Fibroblast tracks are weak proxies.</strong> The Patski line and embryonic-fibroblast RNA are not cardiac fibroblasts; predicted vs measured fibroblast accessibility is {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} (Patski) and {nvf.loc[('heart ATAC', 'atac_peak')]:.2f} (heart). Fibroblast numbers are indicative.</li>
<li><strong>Not run:</strong> saturation mutagenesis of the motif, the human-head sensitivity check and the window-length check (coded, resumable).</li></ul>"""
NEXT = """<ul>
<li><strong>Other-motif shuffle control (highest value).</strong> Pair each AP-1 site with a footprinted non-AP-1 motif (GATA, TBX, MEF2, NFAT) in the same peak plus matched random windows, to test whether AP-1 matters more than other motifs. Start by counting peaks that contain both.</li>
<li>Saturation mutagenesis of the 11&nbsp;bp motif; human-head and window-length sensitivity checks (coded, not run).</li>
<li>Rerun the small code fixes on Katana (MPRA coordinate off-by-one, dinucleotide redraw, insertion control half-sites, record the core-edit fallback) and regenerate the MPRA crossmap.</li>
<li>Confirm whether Fosb, Fosl1 and Fosl2 calls should be included in the site set.</li>
<li>A cardiac-cell-type accessibility probe (the human trial-19 probe, or a mouse probe) would replace the weak fibroblast proxies.</li>
<li>A motif scan on sequences to test the manuscript's enrichment claims.</li></ul>"""

REG = [
 ("Native vs measured accessibility", "Predicted peak ATAC vs measured snATAC, per site", f"rho {nv.rho:.2f} (heart, cardiomyocyte sites); fibroblast sites {nvf.loc[('heart ATAC', 'atac_peak')]:.2f} (heart), {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} (Patski)", "mid", "Informative, not precise; fibroblast tracks are weak proxies."),
 ("Class labels vs measured change", "Measured log2 MAS/YS and MAE/MAS by class", "Gains clear, losses weak", "mid", "Footprint calls and accessibility share the same ATAC data."),
 ("Shuffle at 3 widths", "10 shuffles per site at motif, 170 bp and peak", f"{m.loc['170', 'median_stable']:.2f} (170 bp) vs {m.loc['motif', 'median_stable']:.2f} (motif), Stable", "ok", "Large, consistent, above replicate noise; a model prediction."),
 ("Offset controls and titration", "170 bp shuffles at +-300/600 bp; centred widths 7 to 680 bp", f"Centred {co_.loc['Stable', 'centred_170_median']:.2f} vs {co_.loc['Stable', 'offset_pm300_median']:.2f} at +-300 bp; saturates near 340 bp; 7 bp core {fr_.loc[7]:.0%}", "mid", "Local and saturating is solid; AP-1 specificity is not shown (distance control only)."),
 ("Mutation arm", "Random, dinucleotide, core and flank edits at each width", f"Surroundings alone = {ratio_cm170:.0%} (CM), {ratio_fb170:.0%} (FB) of the full 170 bp loss", "mid", "Ratio near 100% partly expected; track-dependent; the dinucleotide shuffle of the motif is partly a no-op."),
 ("Peak swap", "Donor peak placed into a recipient window; Stable vs Dynamic donors compared", f"Swapped vs donor rho {swap.loc['donor native (sequence)', 'rho']:.2f}; no Stable/Dynamic difference after adjusting for donor accessibility", "mid", "Only 4 recipients and 25 Dynamic donors; pairs not independent; limited power."),
 ("Stable vs Dynamic", "Class comparison at each width, raw and adjusted", f"Raw {m.loc['170', 'diff_dyn_minus_stable']:+.2f}; {hg(cmn, M_MATCH).gap:+.2f} matched-readout adjusted; {hg(cmn, M_FULL).gap:+.2f} with width and TSS", "mid", "Raw difference real; mostly explained by accessibility and context."),
 ("Age vs exercise classes", "Each class vs Stable and against each other", "No significant difference between dynamic classes", "bad", "Underpowered; labels overlap; absence of evidence."),
 ("Fibroblast sites", "Same experiments on 5,018 fibroblast sites, Patski and heart tracks", f"Raw gap replicates; Patski gap {hg(fbp, M_MATCH).gap:+.2f} (ns) after matched adjustment", "mid", "Weak tracks; Patski is not independent confirmation."),
 ("Agreement-filtered sites", "Tiers A to D from native prediction vs measurement", "All main conclusions hold in every tier", "mid", "Filter removes over-predicted Dynamic sites; filtered gaps partly selection."),
 ("Prediction quality by class", "Rho and rank error between binding classes, with tests", f"Rho {XCL[(XCL.combo == 'CM sites, heart ATAC') & (XCL['class'] == 'Stable (all three)')].rho.iloc[0]:.2f} (Stable) vs {XCL[(XCL.combo == 'CM sites, heart ATAC') & (XCL['class'] == 'All Dynamic')].rho.iloc[0]:.2f} (Dynamic); Dynamic over-predicted", "mid", "Errors are similar in size; Dynamic sites are predicted a little less well and biased upwards."),
 ("Footprint score, TSS distance, multi-motif peaks", "Correlations and group comparisons", "Stronger footprints and distal sites depend more on 170 bp context; no buffering", "mid", "Footprint shares data with accessibility; TSS not adjusted; few multi-motif sites."),
 ("MPRA: native level", "Predicted native ATAC vs MPRA log2 alpha", f"rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski", "mid", "Weak but consistent."),
 ("MPRA: scramble effects, genomic and neutral context", "Predicted effect of motif and full scramble vs MPRA change", "About 0 in both contexts", "mid", "MPRA effects are tiny, so the null does not discriminate model quality; 2 neutral backgrounds."),
 ("RNA", "Predicted RNA after shuffles, +-2 kb and +-50 kb readouts", "Median effects about 0", "bad", "Not informative with these tracks."),
 ("Motif insertion into random backgrounds", "TGAGTCA inserted into 50 shuffled backgrounds", f"+{mi('atac_motif', '+motif vs shuffled background').median_log2fc:.2f} log2 at the centre, not significant", "bad", "Floor effect; control can keep half-sites."),
 ("Replicates and GPU", "Replicate convergence; chunk effects across GPU types", "10 replicates adequate; no GPU effect", "ok", "Direct checks."),
]
REG_HTML = "<table class='t reg'><tr><th>Experiment</th><th>What was done</th><th>Result</th><th>Confidence</th><th>Why</th></tr>" + "".join(
    f"<tr><td class='v'>{a}</td><td class='v'>{b}</td><td class='v'>{c}</td><td>{badge(k)}</td><td class='v'>{d}</td></tr>" for a, b, c, k, d in REG) + "</table>"

BIO = f"""<p>Everything below describes what AlphaGenome predicts for these sequences, supported by measured snATAC where stated. It is not a direct measurement of what AP-1 does in cells.</p>
<h3>The surrounding sequence carries the predicted accessibility of an AP-1 site</h3>
<p>Shuffling or mutating the 11&nbsp;bp motif costs a modest part of a site's predicted accessibility, while the 170&nbsp;bp around it carries most of it, with the effect building up to about 340&nbsp;bp. Nearly half of Stable sites are "context only" and few depend on the motif alone. A biologically consistent reading is that AP-1 acts together with neighbouring factors rather than alone: in cardiomyocytes the manuscript highlights MEF2C, GATA4 and TBX-type motifs near AP-1 sites. The present data cannot say which neighbours matter; that is what the other-motif shuffle control would test. A caution is that the readout averages over 170&nbsp;bp, of which 11 are the motif, so some of this is expected from the readout alone.</p>
<h3>Stably bound sites sit in more accessible, more strongly footprinted regions</h3>
<p>Stable sites lose more accessibility than Dynamic ones when their surroundings are removed, but this is largely because they start from more accessible chromatin with stronger footprints (and, for promoter-proximal sites, from the promoter). Once those are accounted for, little remains ({hg(cmn, M_FULL).gap:+.2f} in cardiomyocyte sites). The model therefore does not see a distinct sequence logic in age- or exercise-responsive sites. This agrees with the manuscript's view that dynamic elements are marked by weaker baseline chromatin context rather than by a different core motif. It does not show that sequence is irrelevant to dynamics: the model is trained on baseline accessibility, not on how accessibility responds to age or exercise, so it may simply not be able to represent response-specific logic.</p>
<h3>Age-gained sites are real; the direction of change cannot be resolved</h3>
<p>In independent snATAC, Age Up and Exercise Up sites gain accessibility, while the "down" classes move only weakly in the expected direction. With 79 to 210 sites per class, differences between age and exercise classes in predicted dependence cannot be separated from zero.</p>
<h3>The MPRA reflects something the accessibility head does not capture</h3>
<p>Predicted accessibility explains little of MPRA activity, and the predicted effect of scrambling an element does not track its MPRA change, even in a neutral genomic context. This fits the manuscript's interpretation that episomal activity and chromatin accessibility are different quantities and that chromatin context, not sequence autonomy, explains the endogenous differences, and it supports training a probe directly on the MPRA. The MPRA's own core-scramble effect is small, so a null here is weak evidence about the model.</p>
<h3>What this cannot tell us</h3>
<p>Nothing here tests AP-1 specificity against other motifs, causality in cells, or response-specific regulation. Fibroblast conclusions rest on weak proxy tracks.</p>"""

SUMMARY = f"""
<div class="box"><strong>Summary of results</strong>
<p class="small">Each line shows how far to trust it. Details and plots follow in the numbered sections.</p>
<ul>
<li>{badge('ok')} <strong>Context, not the motif.</strong> Shuffling the 170&nbsp;bp around an AP-1 site removes far more predicted accessibility ({m.loc['170', 'median_stable']:.2f} log2) than shuffling the 11&nbsp;bp motif ({m.loc['motif', 'median_stable']:.2f}). The loss is local, builds with width and saturates near 340&nbsp;bp, and mutating the 170&nbsp;bp with the motif intact keeps {ratio_cm170:.0%} (CM) / {ratio_fb170:.0%} (FB) of it. <em>Meaning:</em> AP-1 sites look cooperative; the surrounding sequence carries most of the predicted accessibility (section {N['depend']}).</li>
<li>{badge('mid')} <strong>Stable vs Dynamic is mostly accessibility.</strong> Stable sites lose more ({m.loc['170', 'diff_dyn_minus_stable']:+.2f} raw), but the gap falls to {hg(cmn, M_MATCH).gap:+.2f}, then {hg(cmn, M_FULL).gap:+.2f} after adjusting for matched-readout accessibility, peak width and TSS distance, and is {hg(fbp, M_MATCH).gap:+.2f} with the Patski track. <em>Meaning:</em> no separate sequence logic is visible for age- or exercise-responsive sites (section {N['dyn']}).</li>
<li>{badge('mid')} <strong>Sequence sets accessibility in a peak swap</strong> (swapped follows donor, rho {swap.loc['donor native (sequence)', 'rho']:.2f}); only 4 recipient windows (section {N['B8a']}).</li>
<li>{badge('mid')} <strong>The model tracks measured accessibility moderately</strong> (rho {nv.rho:.2f} for cardiomyocyte sites, {xov.loc['FB sites, heart ATAC', 'rho']:.2f} for fibroblast sites with the heart track and only {xov.loc['FB sites, Patski ATAC', 'rho']:.2f} with the Patski track). Dynamic sites are predicted a little less well and are over-predicted; the class labels are real in measured snATAC (gains clear, losses weak) (section {N['found']}).</li>
<li>{badge('mid')} <strong>The native accessibility head does not explain the fibroblast MPRA.</strong> Predicted scrambling effects correlate about 0 with the MPRA changes in genomic and neutral context, but the MPRA's own scramble effects are tiny (section {N['mpra']}).</li>
<li>{badge('ok')} <strong>The main results hold in well-predicted sites</strong>, although the filter removes mostly over-predicted Dynamic sites and filtered gaps are partly selection (section {N['robust']}).</li>
<li>{badge('bad')} <strong>Did not work or not reliable:</strong> RNA (effects about zero), motif insertion into random backgrounds (floor effect), differences between age and exercise classes (underpowered), and the fibroblast tracks (weak proxies) (section {N['fail']}).</li></ul>
<p class="small"><strong>Corrected after an independent review.</strong> Three independent reviewers re-derived the headline numbers (all reproduced) and confirmed the labels and site tables. They found that the earlier Stable-vs-Dynamic adjustment used the wrong covariate; it is corrected in section {N['H1']}, other wording was softened, and the code issues still open are listed in section {N['open']}.</p></div>
<div class="box"><strong>How to read the figures and confidence labels</strong>
<p><strong>Effect.</strong> For each site and shuffle width, the mean over 10 shuffles of <code>log2((shuffled + pc) / (native + pc))</code> of the predicted ATAC (mean per-base signal over a region; pc = 5% of the median native value). Negative = accessibility lost. A motif shuffle is read over the motif, a 170&nbsp;bp shuffle over the 170&nbsp;bp, a peak shuffle over the peak, unless stated. <strong>Marks:</strong> <code>***</code> FDR &lt; 0.001, <code>**</code> &lt; 0.01, <code>*</code> &lt; 0.05, <code>ns</code> not significant (two-sided Mann&ndash;Whitney, Benjamini&ndash;Hochberg within each family). With thousands of sites small differences are significant, so each comparison also shows an effect size; read it before the stars.</p>
<p><strong>Confidence.</strong> {badge('ok')} survives the checks and controls; {badge('mid')} real but limited, small or confounded, with the reason given; {badge('bad')} should not be relied on, with the reason given; {badge('na')} not testable with these data.</p></div>"""


# ---- what each result means biologically (shown in a green box after every results section)
_m170, _mmot = m.loc['170', 'median_stable'], m.loc['motif', 'median_stable']
BIO_SEC = {
 "A1": f"AP-1 is bound in all three conditions at {n_st / n_all:.0%} of these cardiomyocyte Fos sites, so AP-1 occupancy here is largely constitutive. Only about a fifth of sites change with age or exercise, divided over six patterns of 65 to 210 sites, so any statement about how AP-1 binding changes with age or exercise rests on small groups.",
 "B1": f"The model has learned enough regulatory sequence logic to rank real genomic regions by accessibility moderately well, so what it says about which sequence matters is meaningful but approximate. The cardiac track ranks fibroblast sites almost as well as cardiomyocyte sites ({xov.loc['FB sites, heart ATAC', 'rho']:.2f} vs {xov.loc['CM sites, heart ATAC', 'rho']:.2f}), which suggests it mostly captures general features of accessible regulatory sequence rather than cardiomyocyte-specific ones. The Patski cell-line track captures much less ({xov.loc['FB sites, Patski ATAC', 'rho']:.2f}), so conclusions that rest on it are weaker.",
 "A2": f"Sites whose AP-1 binding changes sit in narrower peaks (median {gm.loc['Dynamic', 'peak_width']:.0f} vs {gm.loc['Stable', 'peak_width']:.0f} bp) at a similar distance from promoters. The regulatory elements where AP-1 occupancy is dynamic therefore look smaller than those where it is constant, not closer to or farther from genes.",
 "A4": "Where AP-1 binding is gained with age or exercise, chromatin accessibility goes up as well (Age Up and Exercise Up sites gain measured accessibility). Where AP-1 binding is lost, accessibility falls only slightly. AP-1 arriving at a site is therefore clearly linked to opening; AP-1 leaving does not close chromatin to the same degree, which could mean that opened chromatin persists after AP-1 leaves, or that loss calls are less certain.",
 "A3": "Most sites are weakly accessible and a few are very strongly accessible, and sites whose AP-1 binding changes are less accessible than constitutive ones, in the prediction and in the measurement. Responsive AP-1 sites therefore sit in weaker chromatin than constitutively bound sites.",
 "A5": f"The model over-predicts accessibility at Dynamic sites (median rank error {xts.loc['CM sites, heart ATAC', 'dynamic_median_signed']:+.2f} against {xts.loc['CM sites, heart ATAC', 'stable_median_signed']:+.2f} for Stable, and the same in fibroblast sites) while the size of the error does not differ. One reading is that Dynamic sites carry sequence the model finds permissive but that is open only under particular conditions or in particular cell states, so sequence alone does not decide whether these sites are open. This fits their weaker measured accessibility but is an interpretation, not a test.",
 "A6": f"A shuffle that destroys sequence lowers predicted accessibility, as expected if accessibility is written in the sequence. At about a quarter of sites, shuffling only the 11 bp motif raises predicted accessibility, so at those sites the model does not rely on the AP-1 motif itself.",
 "A7": "Technical check; no biological claim. It supports the reliability of the effect sizes.",
 "B7": "Technical check; no biological claim. It means group medians can be trusted but a single site's peak-shuffle value cannot.",
 "A8": "Technical check; no biological claim.",
 "B2": f"Removing the 170&nbsp;bp around an AP-1 motif removes most of the predicted accessibility (about {1 - 2 ** _m170:.0%} lost), while removing only the motif removes about {1 - 2 ** _mmot:.0%}. In the model, accessibility at an AP-1 site is therefore a property of a local stretch of sequence of which the AP-1 motif is one part, not something the motif creates alone. Dynamic sites lose less, mainly because they start with less to lose (section {{DYN}}).",
 "C1": f"Not every AP-1 site relies on its motif: at {an1.loc[('motif', 'Stable'), 'frac_gt0']:.0%} of Stable sites shuffling the motif raises the prediction, and at {an1.loc[('motif', 'Stable'), 'frac_lt_minus1']:.0%} it lowers it by more than half. At many sites AP-1 may therefore be bound to accessibility created by other factors rather than creating it.",
 "D1": f"Sequence changes 300&nbsp;bp away do almost nothing, while the 170&nbsp;bp around the motif carries most of the effect, and the effect stops growing at about 340&nbsp;bp, roughly the length of two nucleosomes. So accessibility at an AP-1 site appears to be encoded in a compact neighbourhood of a few hundred base pairs. The AP-1 core (7&nbsp;bp) accounts for about a quarter of the 170&nbsp;bp effect: an important but minority share. Across the whole peak, the same 170&nbsp;bp accounts for a large part of the peak's accessibility.",
 "D2": f"Replacing everything in the 170&nbsp;bp except the motif keeps {ratio_cm170:.0%} of the loss from replacing it all in cardiomyocyte sites, so the model reads the site's accessibility mainly from neighbouring sequence. With the fibroblast-like track the motif contributes more ({1 - ratio_fb170:.0%} vs {1 - ratio_cm170:.0%}), which could mean the AP-1 motif matters more in a fibroblast-like context, but that track is a weak proxy. Within the motif the 7&nbsp;bp core carries essentially all the effect, as expected for the part that AP-1 recognises.",
 "B8a": "A peak moved into a different genomic window keeps its relative accessibility (the ranking by donor is preserved) but its absolute level moves with the new location: the sequence of a regulatory element carries most of its accessibility potential, and the surrounding chromatin environment scales it up or down. Stable and Dynamic donors behave alike, so elements with dynamic AP-1 binding do not differ in how much accessibility their own sequence carries (with limited power: 25 Dynamic donors).",
 "C3": f"Almost half of Stable sites are 'context only' and few depend on the motif alone ({an3.loc['Motif only', 'share_stable']:.1%}). AP-1 appears to work together with the sequence around it, which is consistent with AP-1 acting alongside lineage factors such as GATA, TBX or MEF2 in cardiomyocytes, although this analysis cannot name those partners.",
 "C2": "Sites in more accessible and narrower peaks depend more on the local sequence because there is more accessibility to lose; the effect weakens at the very top, probably from saturation. Strongly open elements therefore tend to rest on local sequence information, and removing it costs them more.",
 "C6": "Promoter-proximal AP-1 sites lose less from a 170&nbsp;bp shuffle than distal ones although their peaks are much more accessible. Accessibility at promoters appears to be maintained by features other than the AP-1 neighbourhood, whereas distal, enhancer-like AP-1 sites depend more on their local sequence. This comparison is not adjusted for accessibility.",
 "C4": "Strongly footprinted AP-1 sites are predicted to depend more on their surroundings, so well-bound sites are those embedded in accessible, motif-rich sequence. Footprint strength and accessibility come from the same data, so this partly restates that open chromatin supports binding.",
 "C5": "There is no sign that a second nearby Fos motif makes a site more dispensable. But with few sites in multi-motif peaks and edits that change one window at a time, this cannot rule out redundancy between neighbouring AP-1 sites.",
 "B3": "AP-1 sites whose binding changes with age or exercise lose less accessibility when their surroundings are removed. The next section tests whether that reflects different biology or simply different chromatin context.",
 "B4": "Elements that gain AP-1 with age or with exercise look the same to the model as those that lose it, so the model cannot tell age-driven from exercise-driven AP-1 changes by sequence. With these sample sizes, differences up to about half the size of the Stable-vs-Dynamic gap could be hidden.",
 "D3": "The same picture holds at fibroblast AP-1 sites: a raw difference that is explained by accessibility. The cardiomyocyte findings are therefore not specific to one cell type, although the weak fibroblast tracks make this indicative only.",
 "B6": "Restricting to well-predicted sites does not change the conclusions, so they are not driven by regions where the model is unreliable.",
 "G2": "The filter removes mostly Dynamic sites that the model over-predicts. The Dynamic sites it keeps are those whose accessibility the model can explain from sequence, so the filtered result describes that subset of dynamic AP-1 sites and not all of them.",
 "G3": "Among Stable sites, those the model predicts well also show the largest dependence on their surroundings, which means the filtered Stable-vs-Dynamic gaps are partly a product of selecting well-explained sites.",
 "G4": "The qualitative picture is the same in every tier: a larger loss in Stable sites, which in the previous chapter is largely explained by chromatin context.",
 "G6": "The finding that surrounding sequence carries the predicted accessibility does not depend on which sites are kept.",
 "C7": "A model of chromatin accessibility ranks which sequences are open, and does so weakly for MPRA activity, but it does not predict how much reporter activity changes when the core motif or flanks of a 170&nbsp;bp element are scrambled on a plasmid. This fits accessibility and enhancer activity being different properties that sequence changes move differently, which is the manuscript's own reason for training a probe directly on the MPRA.",
 "D4": "Taking the elements out of their genomic setting does not rescue the prediction, so the mismatch is not because the model uses the surrounding genome. It more likely reflects that the MPRA effects are small and that accessibility is not enhancer activity.",
 "B5": "Shuffling a single distal AP-1 site or peak barely changes predicted expression of nearby genes. This is expected if individual distal elements each contribute little to a gene's expression, or if the model does not capture enhancer-gene effects; nothing can be concluded about regulation of target genes.",
 "B8b": "No conclusion can be drawn: the model's lack of response cannot be separated from the floor effect, although it is consistent with the AP-1 motif not being sufficient to create accessibility on its own.",
}

# ---- assemble
toc = ["<div class='toc'>"]
for ck, ct, cintro, subs in CH:
    toc.append(f"<div><a href='#ch-{ck}'><strong>{N[ck]}. {ct}</strong></a><ul>" + "".join(f"<li><a href='#sec-{sk}'>{N[sk]} {st}</a></li>" for sk, st in subs) + "</ul></div>")
toc.append("</div>")

SPECIAL = {"MAN": MAN, "NOTE": NOTE, "NEXT": NEXT}
chapters = []
for ck, ct, cintro, subs in CH:
    h = f"<h2 id='ch-{ck}'><span class='num'>{N[ck]}</span>{ct}</h2>"
    if ck == "design":
        chapters.append(h + "<p>The figures below show what each experiment does to the sequence and what is measured, so that the results that follow can be read without the methods.</p>" + fixrefs(re.sub(r"<h2.*?</h2>", "", DESIGN, count=1, flags=re.S)))
        continue
    if ck == "reg":
        chapters.append(h + "<p>All experiments in one table, with how far to trust each and why. The detailed claim-by-claim table follows.</p>" + fixrefs(REG_HTML + "<h3>Claim-by-claim evidence</h3><table class='t reg'><tr><th>Claim</th><th>Confidence</th><th>Evidence</th></tr>" + claims_html + "</table>"))
        continue
    if ck == "bio":
        chapters.append(h + BIO)
        continue
    body_ = fixrefs(eval('f"""' + cintro + '"""') if cintro else "")
    parts_ = [h, body_]
    for sk, st in subs:
        c_ = conf(*CONF[sk]) if sk in CONF else ""
        content = SPECIAL[sk] if sk in SPECIAL else SEC[sk]
        bx = f"<div class='bio'><strong>What this means biologically.</strong> {BIO_SEC[sk]}</div>" if sk in BIO_SEC else ""
        parts_.append(f"<h3 id='sec-{sk}'><span class='num'>{N[sk]}</span>{st}</h3>{c_}{fixrefs(content)}{fixrefs(bx.replace('{DYN}', N['dyn']))}")
    chapters.append("\n".join(parts_))

METHODS = """<details><summary>Methods and files</summary>
<p class="small">Predictions: alphagenome-pytorch with <code>model_all_folds</code> weights, mouse organism, 131,072&nbsp;bp window centred on the motif midpoint, mouse heart ATAC track (1&nbsp;bp, mean over region) and Patski ATAC for fibroblast sites, 10 uniform random shuffles per width. Stable = bound in YS, MAS and MAE; age class compares YS and MAS; exercise class compares MAS and MAE. Class comparisons: Mann&ndash;Whitney, BH within families, bootstrap CIs, Cliff's delta; adjusted models are OLS with HC3 (or peak-clustered) standard errors. Sites are treated as independent in most tests, but about 15% share a peak, so p-values are mildly optimistic; effect sizes and intervals are the better guide. Tables: <code>Analysis/cm/tables</code>, <code>Analysis/cm/stats</code>, <code>Analysis/overnight</code>, <code>Analysis/agreement</code>, <code>Analysis/agreement_floor25</code>, <code>Analysis/audit</code>; scripts: <code>Scripts/22&ndash;31</code> and <code>25_build_report.py</code>.</p></details>"""

# extra plots kept in a collapsed block so they can be pruned later
EXTRA = f"""<details><summary>Additional plots (alternative versions, for choosing which to keep)</summary>
<figure>{img(AN / 'agreement' / 'figures' / 'e1_tier_C_selection.png', 'Tier C selection')}<figcaption>Predicted native ATAC against measured snATAC with tier C sites in blue (absolute values).</figcaption></figure>
<figure>{img(AN / 'agreement' / 'figures' / 'e2_retention_by_class.png', 'Retention by class')}<figcaption>Share of each class retained per tier.</figcaption></figure>
<figure>{img(AN / 'agreement_floor25' / 'figures' / 'e3_stable_vs_dynamic_by_tier.png', 'Floor 25 stable vs dynamic')}<figcaption>Stable versus Dynamic by tier with the stricter 25th-percentile floor.</figcaption></figure>
<figure>{img(AN / 'agreement_floor25' / 'figures' / 'e0b_rank_map.png', 'Floor 25 rank map')}<figcaption>Rank map with the stricter 25th-percentile floor.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig2_shuffle_effect_heatmaps.png', 'Heatmaps without marks')}<figcaption>Shuffle heatmaps without significance marks.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig3_shuffle_effect_by_class.png', 'Class effects without marks')}<figcaption>Effect by class without significance marks.</figcaption></figure></details>"""

head = f"""<h1>AlphaGenome in-silico perturbation of AP-1 (Fos) sites in mouse cardiomyocytes and fibroblasts</h1>
<p class="sub">Generated {date.today():%d %b %Y} &middot; {int(m.loc['170', 'n_stable'] + m.loc['170', 'n_dynamic']):,} cardiomyocyte and 5,018 fibroblast Fos sites &middot; mouse heart and Patski ATAC tracks &middot; shuffle, mutation, insertion, control and MPRA experiments</p>"""
CHKEY = {v: k for k, v in N.items() if k in {c[0] for c in CH}}
NUMKEY = {v: k for k, v in N.items()}
SUMMARY = re.sub(r"\(section ([\d.]+)\)", lambda mm: f'(<a href="#{"ch" if "." not in mm.group(1) else "sec"}-{NUMKEY[mm.group(1)]}">section {mm.group(1)}</a>)', SUMMARY)
body = head + SUMMARY + "<h2>Contents</h2>" + "".join(toc) + "\n".join(chapters) + EXTRA + METHODS
doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AP-1 shuffle analysis report</title><style>{CSS}</style></head><body><main>{body}</main></body></html>"""
for k_, p_ in IMGS.items():
    doc = doc.replace(f"@@IMG{k_}@@", "data:image/png;base64," + base64.b64encode(p_.read_bytes()).decode())
out = AN / "report_cm.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB), {len(IMGS)} images")
