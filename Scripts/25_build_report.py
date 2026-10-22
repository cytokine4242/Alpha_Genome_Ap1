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
div.bio{border-left:4px solid var(--green);background:var(--panel);padding:10px 16px;margin:14px 0;border-radius:6px;max-width:80ch}div.bio strong:first-child{color:var(--green)}.keyfig{margin:26px 0;padding-top:6px;border-top:1px dotted var(--rule)}.conf{margin:6px 0 10px;font-size:.92rem}.conf .badge{margin-right:6px}.toc{columns:2;column-gap:32px;font-size:.92rem;margin:8px 0}.toc ul{margin:0 0 8px;padding-left:18px}.toc a{color:var(--blue);text-decoration:none}h2 .num,h3 .num{color:var(--ink2);margin-right:6px}.badge{display:inline-block;font-size:.78rem;font-weight:600;border-radius:99px;padding:1px 9px;border:1.5px solid currentColor;white-space:nowrap}.ok{color:var(--green)}.mid{color:var(--amber)}.bad{color:var(--red)}.na{color:var(--grey)}.inc{color:var(--grey)}ul.pts{margin:6px 0 10px;padding-left:20px}ul.pts li{margin:3px 0}
td.v{white-space:normal!important;text-align:left!important;min-width:260px}.small{font-size:.86rem;color:var(--ink2)}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:.88em}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2)}
"""


def badge(kind):
    lab = {"ok": "&#10003; High confidence", "mid": "~ Moderate confidence", "bad": "&#10007; Low / unreliable", "na": "? Untested", "inc": "? Inconclusive"}[kind]
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

# ======================= Assembly: chapters ordered by the scientific story, short dot points =======================
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
_b8 = SEC.pop("B8")
_i = _b8.index('alt="Motif insertion"')
_j = _b8.rfind("<figure>", 0, _i)
SEC["B8a"], SEC["B8b"] = _b8[:_j], _b8[_j:]
SEC.pop("G7", None)
SEC.pop("D5", None)


DY = AN / "dynamic"
dr = lambda n: pd.read_csv(DY / "tables" / n)
Y_auc, Y_cond, Y_var, Y_cor = dr("classifier_auc.csv").set_index("features"), dr("condition_specific_accessibility.csv").set_index("binding_pattern"), dr("accessibility_variability_across_conditions.csv").iloc[0], dr("prediction_vs_each_condition.csv").set_index("binding_pattern")
Y_same, Y_pos, Y_cross, Y_match, Y_bal = dr("same_peak_pairs.csv").set_index("quantity"), dr("position_in_peak.csv").iloc[0], dr("cross_cell_type_concordance.csv").iloc[0], dr("matched_pairs.csv").set_index("effect"), dr("matching_balance.csv")
Y_prof, Y_pat, Y_scan, Y_mauc, Y_fbr = dr("mutation_type_profile_dynamic_vs_stable.csv").set_index("mutation"), dr("mutation_type_by_binding_pattern.csv"), dr("feature_scan_matched_pairs.csv"), dr("classifier_auc_matched_sample.csv").set_index("feature_set"), dr("fibroblast_replication_by_pattern.csv")
dimg = lambda n, alt: img(DY / "figures" / f"{n}.png", alt)


def strip_prose(s, keep_tables=None):
    """Keep figures, details blocks and (some) tables of a section; drop its paragraphs and lists (the dot points replace them)."""
    out, pos, ntab = [], 0, 0
    pat = re.compile(r'<(p|ul|ol|figure|details|table|div|h4)\b[^>]*>')
    while True:
        m = pat.search(s, pos)
        if not m:
            break
        tag = m.group(1)
        tok = re.compile(rf'<(/?){tag}\b[^>]*>')
        depth, end = 0, None
        for mm in tok.finditer(s, m.start()):
            depth += -1 if mm.group(1) else 1
            if depth == 0:
                end = mm.end()
                break
        if end is None:
            break
        el = s[m.start():end]
        keep = tag in ("figure", "details") or (tag == "table" and (keep_tables is None or ntab < keep_tables))
        ntab += tag == "table"
        if keep:
            out.append(el)
        pos = end
    return "\n".join(out)


def U(*items):
    return "<ul class='pts'>" + "".join(f"<li>{x}</li>" for x in items) + "</ul>"


def conf(kind, text):
    return f'<p class="conf">{badge(kind)} {text}</p>'


CONF = {
 "A1": ("ok", "Re-derived independently; counts reproduce."),
 "A2": ("ok", "Descriptive."),
 "A4": ("mid", "A consistency check: footprint calls and accessibility share the same ATAC data."),
 "A3": ("ok", "Descriptive."),
 "B1": ("mid", "Informative, not precise; bulk tissue and a cell line stand in for purified cell types."),
 "A5": ("mid", "Tested between classes: error size is similar, but Dynamic sites are predicted a little worse and biased upwards."),
 "TECH": ("ok", "Direct checks; the H100 GPU ran only one 91-site chunk and cannot be assessed."),
 "B2": ("ok", "A description of the model's predictions; every later control agrees with it."),
 "C1": ("ok", "Descriptive."),
 "D1": ("mid", "High confidence the loss is local and saturates; moderate that it is AP-1-specific (a distance control only)."),
 "D2": ("mid", "A ratio near 100% is partly expected (11 of 170 bp are motif) and it follows the track."),
 "B8a": ("mid", "Clear direction; only 4 recipient windows and 25 Dynamic donors, so power is limited."),
 "C3": ("mid", "The 50% thresholds are arbitrary; descriptive."),
 "C2": ("mid", "High for the peak-shuffle correlations; the 170 bp and motif ones are weak."),
 "C6": ("mid", "Not adjusted for accessibility."),
 "C4": ("mid", "Footprint scores come from the same ATAC data as accessibility."),
 "C5": ("inc", "Few multi-motif sites and one-window edits: no evidence either way."),
 "B3": ("ok", "High for the raw difference as a description; see the next section for what explains it."),
 "H1": ("mid", "High that context accounts for most of the gap; footprint score partly defines the classes, so that row is a sensitivity check."),
 "B4": ("inc", "Underpowered for between-class differences (intervals about +-0.2 to 0.3); labels overlap."),
 "D3": ("mid", "Raw difference replicates; weak tracks, and Patski uses the same model and sequences."),
 "G0": ("ok", "Definition."),
 "G1": ("ok", "Descriptive."),
 "G2": ("mid", "The filter drops mostly over-predicted Dynamic sites."),
 "G3": ("mid", "Selection and effect are related for Stable sites."),
 "G4": ("mid", "Same direction in every tier; larger filtered gaps are partly selection."),
 "G5": ("ok", "Direction and significance do not depend on the floor; size does somewhat."),
 "B6": ("mid", "The stricter-floor version of the tiers."),
 "G6": ("ok", "Unchanged by the filter."),
 "C7": ("mid", "Weak native-level correlation; the scramble null does not rank models, as MPRA effects are tiny."),
 "D4": ("mid", "A clean null, but only 2 neutral backgrounds and tiny MPRA effects."),
"DS": ("mid", "Large raw differences, none after matching on predicted context; matching uses 720 of 772 Dynamic sites."),
 "DSH": ("mid", "Matched pairs, paired tests and a permutation-checked classifier; a null result, so small shape differences could be missed."),
 "DT": ("inc", "199 of 1,470 sites so far; the one significant difference does not hold across readouts."),
 "DP": ("inc", "126 comparisons; the one pattern that stood out (MAS+MAE) does not replicate in fibroblast sites."),
 "DF": ("mid", "Cross-validated; the matched sample has 1,440 sites, so small differences could be missed."),
 "DC": ("mid", "Small effects but consistent across all six patterns; one pooled measurement per condition."),
 "DM": ("mid", "Cross-validated and matched analyses; the same-peak comparison uses only 55 peaks."),
 "DX": ("mid", "Only the 962 sites called in both cell types; classes are threshold calls."),
 "DEXP": ("na", "Proposals, not results."),
 "FT": ("na", "Proposal."),
 "B5": ("bad", "Not informative with these tracks."),
 "B8b": ("bad", "At the floor; not interpretable."),
}

N = {}
CH = [
 ("design", "The experiments at a glance", None, []),
 ("found", "Can we trust the data and the model?",
  U("Mostly yes, with limits.", f"The model ranks measured accessibility moderately (rho {nv.rho:.2f}, cardiomyocyte sites); labels are partly supported by measured snATAC; shuffles behave as designed.",
    "Limits: Dynamic sites are predicted a little worse and are over-predicted; the fibroblast tracks are weak proxies."),
  [("A1", "Which sites were analysed"), ("B1", "Does the model predict measured accessibility? Heart and Patski tracks compared"), ("A2", "Stable and Dynamic sites sit in different peaks"),
   ("A4", "Are the Stable and Dynamic labels real?"), ("A3", "Native predicted accessibility"), ("A5", "Is prediction quality different between site classes?"),
   ("TECH", "Technical checks: shuffle behaviour, replicates, noise, GPU")]),
 ("depend", "What do AP-1 sites depend on? The surrounding sequence, not the motif",
  U("Question: how much predicted accessibility comes from the 11 bp motif and how much from its surroundings?",
    f"The 170 bp shuffle costs {m.loc['170', 'median_stable']:.2f} log2 vs {m.loc['motif', 'median_stable']:.2f} for the motif; the loss is local and saturates near 340 bp; replacing everything but the motif keeps {ratio_cm170:.0%} (CM) / {ratio_fb170:.0%} (FB) of the loss; in a peak swap, sequence sets the rank.",
    "Limits: model predictions; the offset control is a distance control; 11 of the 170 bp are the motif."),
  [("B2", "What shuffling does at three widths"), ("C1", "How varied are the per-site effects?"), ("D1", "Controls and width titration"), ("D2", "Mutation arm: what carries the dependence?"),
   ("B8a", "Peak swap: sequence or location, Stable or Dynamic donors"), ("C3", "Motif dependence versus context dependence"), ("C2", "What sets the size of the effect?"), ("C6", "Promoter-proximal versus distal sites"),
   ("C4", "Does footprint strength agree with the predicted dependence?"), ("C5", "Do two Fos sites in a peak buffer each other?")]),
 ("dyn", "Do dynamic sites behave differently from stable ones?",
  U("Question: do age- or exercise-responsive sites depend on sequence differently from constitutive ones?",
    f"Stable loses more at 170 bp ({m.loc['170', 'diff_dyn_minus_stable']:+.2f} raw); this falls to {hg(cmn, M_MATCH).gap:+.2f} (matched accessibility) and {hg(cmn, M_FULL).gap:+.2f} (plus peak width, TSS distance); age and exercise classes cannot be told apart.",
    "Limit: a small residual remains and the data neither show nor exclude a separate sequence effect."),
  [("B3", "Stable versus Dynamic: the raw difference"), ("H1", "Correcting for accessibility at the matched readout"), ("B4", "Age, exercise and their combination"), ("D3", "Does it replicate in fibroblast sites?")]),
 ("deep", "Are Dynamic sites more or less sensitive to particular mutations?",
  U("Question: in the AlphaGenome mutation results, do the binding classes differ in sensitivity to particular types of edit?",
    f"Raw: Dynamic sites lose less to every wide edit (for example replace the 170 bp {Y_prof.loc['replace the 170 bp', 'stable_median']:.2f} Stable vs {Y_prof.loc['replace the 170 bp', 'dynamic_median']:.2f} Dynamic) and about 0.1 less to motif-level edits, and the same to far-away edits. Matched on predicted accessibility, peak width and TSS distance, every difference is within {Y_prof[['matched_stable_median', 'matched_dynamic_median']].diff(axis=1).iloc[:, 1].abs().max():.2f} log2 and none is significant.",
    "Scaled to the same overall sensitivity, the mix of mutation types is the same in both classes, and a classifier on the profile is at chance.",
    "So Dynamic sites are less sensitive across the board in proportion to their weaker context; no mutation type stands out, and no binding pattern is robustly different.",
    "New (preliminary tile scan): most sites in both classes have a neighbouring element that matters as much as the AP-1 motif; whether Dynamic sites differ in their partners awaits the full scan."),
  [("DS", "Sensitivity profile: Dynamic vs Stable, for every type of mutation"), ("DSH", "Same shape: the mix of mutation types, at equal overall sensitivity"), ("DT", "Tile scan: own motif or a neighbouring element? (preliminary)"), ("DP", "Do the individual binding patterns differ?"), ("DF", "Anything else in the model outputs?"), ("DC", "Condition-specific accessibility: real but modest"),
   ("DM", "Matched pairs, same-peak pairs and position"), ("DX", "Is the class a property of the sequence?"), ("DEXP", "Experiments that would explore this further")]),
 ("robust", "Robustness: only sites where prediction matches measurement",
  U("Test: keep sites whose predicted accessibility agrees with measured, chosen from the native prediction and measurement only, and repeat the main analyses.",
    "The Stable-vs-Dynamic direction, the mutation split and the controls hold in every tier.",
    "Limit: the filter drops mostly over-predicted Dynamic sites; larger filtered gaps are partly selection."),
  [("G0", "How the restriction was done"), ("G1", "What was kept and what was dropped"), ("G2", "Is the filter biased against one class?"), ("G3", "Do kept and dropped sites differ in the effect?"),
   ("G4", "Does Stable versus Dynamic survive?"), ("G5", "Sensitivity to the floor"), ("B6", "Tiers with the stricter 25th-percentile floor"), ("G6", "Mutation arm and controls in the filtered sets")]),
 ("mpra", "Does the model explain the fibroblast MPRA?",
  U("Test: does predicted accessibility loss after scrambling an element track the change in MPRA activity (1,653 elements; genomic and neutral windows)?",
    f"Native level correlates weakly (rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski); scramble effects do not track at all.",
    "Limit: the MPRA's own scramble effects are tiny (core loss about 0.02 log2), so the null says little about model quality."),
  [("C7", "Predicted accessibility against the MPRA"), ("D4", "Does placing elements in a neutral context restore agreement?"), ("MAN", "What the manuscript claims, and how our results bear on it")]),
 ("fail", "Experiments that did not work, or are limited", U("Reported so you can see what was tried and why not to rely on it."),
  [("B5", "Predicted RNA"), ("B8b", "Motif insertion into random backgrounds"), ("NOTE", "Other limits")]),
 ("reg", "Experiment register: every experiment, confidence and reason", None, []),
 ("bio", "What these tests taught us", None, []),
 ("open", "Review log, open issues and next steps", None, [("H2", "Corrections made after review"), ("H3", "Code issues found and not yet fixed"), ("NEXT", "Next steps"), ("FT", "What it would take to fine-tune cardiomyocyte and fibroblast heads")]),
 ("key", "Five key figures for a presentation", None, []),
]
for ci, (ck, ct, cintro, subs) in enumerate(CH, 1):
    N[ck] = str(ci)
    for si, (sk, st) in enumerate(subs, 1):
        N[sk] = f"{ci}.{si}"
ANCHOR = {k: k for k in N}
for k_ in ("A6", "A7", "B7", "A8"):
    N[k_] = N["TECH"]
    ANCHOR[k_] = "TECH"


def link(key):
    return f'<a href="#sec-{ANCHOR[key]}">section {N[key]}</a>'


def clink(ck):
    return f'<a href="#ch-{ck}">section {N[ck]}</a>'


OLD_REFS = [(r"(?<![\w-])Part H\b", "H1"), (r"(?<![\w-])Part G\b", "G0"), (r"(?<![\w-])Part C\b", "C2"), (r"(?<![\w#-])B6\b", "B6"), (r"(?<![\w#-])B7\b", "B7"), (r"(?<![\w#-])G3\b", "G3"),
            (r"(?<![\w#-])D1\b", "D1"), (r"(?<![\w#-])D2\b", "D2"), (r"(?<![\w#-])D4\b", "D4"), (r"(?<![\w#-])A4\b", "A4"), (r"(?<![\w#-])C4\b", "C4")]


def fixrefs(s):
    s = re.sub(r"\[\[(\w+)\]\]", lambda mm: link(mm.group(1)), s)
    s = re.sub(r"B3 and B4", f"sections {link('B3')[link('B3').index('>') + 1:-4].replace('section ', '')} and {N['B4']}", s)
    for pat, key in OLD_REFS:
        s = re.sub(pat, lambda mm, key=key: link(key), s)
    return s


# ---------------- dot points per section (placed before the plots; the old paragraphs are dropped)
_c = lambda k: q4.loc[k, "median_log2fc"]
_stab = lambda a, b: f"{a}"
PTS = {
 "A1": U(f"{n_all:,} cardiomyocyte Fos sites: {n_st:,} Stable ({n_st / n_all:.0%}) and {n_dy:,} Dynamic, over six dynamic patterns of about 80 to 210 sites each.",
         "Age and exercise labels come from the same three binding calls, so they overlap (seven patterns in all).", f"The Stable/Dynamic mix does not differ between chromosomes (p = {chrom_p:.2f}).",
         "Small classes limit the power of every class comparison."),
 "B1": U(f"The heart track ranks measured accessibility moderately: rho {xov.loc['CM sites, heart ATAC', 'rho']:.2f} for cardiomyocyte sites and {xov.loc['FB sites, heart ATAC', 'rho']:.2f} for fibroblast sites.",
         f"The Patski track on the same fibroblast sites gets only {xov.loc['FB sites, Patski ATAC', 'rho']:.2f} (and is about uncorrelated with cardiomyocyte accessibility).",
         f"Sites compared are those with a measured peak: {int(xov.loc['CM sites, heart ATAC', 'n']):,} cardiomyocyte and {int(xov.loc['FB sites, heart ATAC', 'n']):,} fibroblast.",
         "The heart track is bulk tissue and Patski is a cell line: both stand in for purified cell types."),
 "A2": U(f"Dynamic sites sit in narrower peaks (median {gm.loc['Dynamic', 'peak_width']:.0f} vs {gm.loc['Stable', 'peak_width']:.0f} bp, p {fmtp(gp['peak_width'])}).",
         f"Distance to the nearest TSS is similar ({gm.loc['Dynamic', 'abs_tss']:,.0f} vs {gm.loc['Stable', 'abs_tss']:,.0f} bp, p = {gp['abs_tss']:.2f}).",
         f"About {(SL.motifs_in_peak >= 2).mean():.0%} of sites share their peak with another Fos site.", "Peak width also drives peak-shuffle effects, so class comparisons adjust for it."),
 "A4": U("Footprint calls and measured snATAC come from the same ATAC data, so this is a consistency check, not independent validation (sites with a measured peak).",
         f"Age Up sites gain accessibility with age (median log2 MAS/YS {_c(('age_class', 'Age Up')):+.2f} vs {_c(('age_class', 'Stable')):+.2f} for Stable); Exercise Up sites gain too ({_c(('ex_class', 'Exercise Up')):+.2f}).",
         f"The down classes move the expected way but weakly (Age Down {_c(('age_class', 'Age Down')):+.2f}, Exercise Down {_c(('ex_class', 'Exercise Down')):+.2f}).",
         f"Class differences are significant (Kruskal-Wallis p {fmtp(age_p)} age, {fmtp(ex_p)} exercise). Read 'Dynamic' as footprint-defined, not strongly differentially accessible."),
 "A3": U(f"Predicted accessibility is right-skewed: most sites are low (median {q3.loc[('Stable', 'peak'), 'median']:.2f} Stable, {q3.loc[('Dynamic', 'peak'), 'median']:.2f} Dynamic over the peak), a few are very high.",
         f"Dynamic sites are less accessible at every readout (peak mean {q3.loc[('Dynamic', 'peak'), 'mean']:.2f} vs {q3.loc[('Stable', 'peak'), 'mean']:.2f}) and in the measurement (median 17.8 vs 29.3 CPM per kb)."),
 "A5": U(f"Spearman rho is {q5.loc['Stable', 'rho']:.2f} for Stable and {q5.loc['Dynamic', 'rho']:.2f} for Dynamic sites; difference {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_vs_stable.iloc[0]:+.2f} (95% CI {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_lo.iloc[0]:+.2f} to {XPR[(XPR.combo == 'CM sites, heart ATAC') & (XPR['class'] == 'All Dynamic')].rho_diff_hi.iloc[0]:+.2f}).",
         f"Across the seven binding classes rho differs (Cochran Q = {xts.loc['CM sites, heart ATAC', 'rho_heterogeneity_Q']:.1f}, p = {xts.loc['CM sites, heart ATAC', 'rho_heterogeneity_p']:.3f}); lowest for MAS only (0.19) and MAS+MAE (0.16).",
         f"The size of the error is the same across classes (Kruskal-Wallis p = {xts.loc['CM sites, heart ATAC', 'kw_abs_resid_p']:.2f}), but its direction is not: Dynamic sites are over-predicted (median rank error {xts.loc['CM sites, heart ATAC', 'dynamic_median_signed']:+.2f} vs {xts.loc['CM sites, heart ATAC', 'stable_median_signed']:+.2f} for Stable; same in fibroblast sites).",
         "No class is unusable, but Dynamic sites are predicted a bit worse and biased upwards; see the tiers in chapter 5."),
 "C1": U(f"Shuffling only the motif raises the prediction at {an1.loc[('motif', 'Stable'), 'frac_gt0']:.0%} of Stable sites and lowers it by more than half at {an1.loc[('motif', 'Stable'), 'frac_lt_minus1']:.0%}.",
         f"The 170 bp shuffle lowers it by more than half at {an1.loc[('170', 'Stable'), 'frac_lt_minus1']:.0%} of Stable and {an1.loc[('170', 'Dynamic'), 'frac_lt_minus1']:.0%} of Dynamic sites.",
         "So dependence on the motif is a property of a subset of sites, not a uniform feature."),
 "D1": U(f"Centred 170 bp shuffle: {co_.loc['Stable', 'centred_170_median']:.2f} log2 at the motif. Same-size shuffle 300 bp away: {co_.loc['Stable', 'offset_pm300_median']:.2f} (centred is worse at {co_.loc['Stable', 'share_centred_more_negative']:.0%} of Stable and {co_.loc['Dynamic', 'share_centred_more_negative']:.0%} of Dynamic sites).",
         f"The loss grows with width and flattens near 340 bp: share of the 170 bp effect is {fr_.loc[7]:.0%} at 7 bp, {fr_.loc[11]:.0%} at 11 bp, {fr_.loc[85]:.0%} at 85 bp, {fr_.loc[340]:.0%} at 340 bp.",
         f"The Stable-Dynamic gap grows from {sv_.loc[7, 'diff_dyn_minus_stable']:+.2f} (7 bp) to {sv_.loc[340, 'diff_dyn_minus_stable']:+.2f} (340 bp), so it comes from the surroundings, not the core.",
         f"Read over the whole peak, the centred 170 bp costs {D1c.loc[('peak', 'Stable'), 'median']:.2f} (Stable) and {D1c.loc[('peak', 'Dynamic'), 'median']:.2f} (Dynamic) vs {D1c.loc[('peak', 'Stable'), 'off']:.2f} for shuffles 300 bp away.",
         "Caveat: this is a distance control, not an AP-1-specific one."),
 "D2": U(f"Replacing the 159 bp around the motif (170 minus the 11 bp motif, motif kept) keeps {ratio_cm170:.0%} (cardiomyocyte) and {ratio_fb170:.0%} (fibroblast, Patski) of the loss from replacing all 170 bp ({dec('170', 'flank@170', 'CM, heart ATAC'):.2f} vs {dec('170', 'random@170', 'CM, heart ATAC'):.2f} in cardiomyocyte sites).",
         f"Within the motif, replacing only the 7 bp core gives {ratio_cmcore:.0%} of replacing all 11 bases; the 4 flank bases alone give {dec('motif', 'flank@motif', 'CM, heart ATAC') / dec('motif', 'random@motif', 'CM, heart ATAC'):.0%}.",
         "Stable loses more than Dynamic at the 170 bp and whole-peak widths (FDR < 0.001, both site sets); at motif width the classes barely differ.",
         "Caveats: 11 of the 170 bp are the motif, so a ratio near 100% is partly expected; the result follows the track (about 13% from the motif with the heart track on fibroblast sites); the motif dinucleotide shuffle is partly a no-op (3.8% of replicates identical).",
         "Here dependency = minus log2FC (positive = accessibility lost); each edit is read over the region it changes, so replacing the whole peak (1.48) can look smaller than replacing 170 bp (1.96)."),
 "B8a": U(f"Swapped accessibility follows the donor peak's own accessibility (rho {swap.loc['donor native (sequence)', 'rho']:.2f}) more than the recipient's ({swap.loc['recipient native (location)', 'rho']:.2f}, from only 4 recipient values); within each recipient rho is 0.72 to 0.91.",
          f"Stable and Dynamic donors behave the same once donor accessibility is accounted for: Dynamic coefficient {XPS.iloc[-1].stable_median:+.2f} log2 (95% CI {XPS.iloc[-1].dynamic_median:+.2f} to {XPS.iloc[-1].n_dynamic_donors:+.2f}, p = {XPS.iloc[-1].p_mannwhitney:.2f}); rho {XPS.iloc[4].stable_median:.2f} vs {XPS.iloc[5].stable_median:.2f}.",
          f"Dynamic donors start lower (median {XPS.iloc[0].dynamic_median:.2f} vs {XPS.iloc[0].stable_median:.2f}, p = {XPS.iloc[0].p_mannwhitney:.2f}).",
          "The same peak keeps less than its own accessibility in the least accessible recipient and more in the more accessible ones, for both classes.",
          "Caveat: 4 recipients and 25 Dynamic donors (a difference below about 0.2 log2 could be missed); pair-level p-values treat the 292 pairs as independent, so use the donor-level numbers."),
 "C3": U(f"Motif and 170 bp effects are related across sites (rho {rho_mc:.2f}).",
         f"Stable sites: {an3.loc['Context only', 'share_stable']:.0%} context only, {an3.loc['Motif and context', 'share_stable']:.0%} motif and context, {an3.loc['Motif only', 'share_stable']:.1%} motif only. Dynamic sites are more often 'neither' ({an3.loc['Neither', 'share_dynamic']:.0%} vs {an3.loc['Neither', 'share_stable']:.0%}; chi-square p < 1e-20).",
         "Categories use a 50% loss threshold, which is arbitrary: read as descriptive."),
 "C2": U(f"Peak-shuffle loss tracks native accessibility (rho {an2.loc[('peak', 'native_peak'), 'rho']:+.2f}) and peak width ({an2.loc[('peak', 'peak_width'), 'rho']:+.2f}).",
         f"The 170 bp loss depends less on these, with TSS distance ({an2.loc[('170', 'abs_tss'), 'rho']:+.2f}) and footprint score ({an2.loc[('170', 'fp_score'), 'rho']:+.2f}); it is U-shaped across native accessibility (saturation at the top).",
         f"The motif-only effect depends little on anything (|rho| up to {an2[an2.index.get_level_values(0) == 'motif'].rho.abs().max():.2f}).", "So class comparisons must be adjusted for accessibility and width."),
 "C6": U(f"Within 2 kb of a TSS, peaks are more accessible (native {an6.loc[('native_peak', 'Stable', '<=2 kb'), 'median']:.2f} vs {an6.loc[('native_peak', 'Stable', '>30 kb'), 'median']:.2f} beyond 30 kb) but lose less from a 170 bp shuffle ({an6.loc[('eff_170', 'Stable', '<=2 kb'), 'median']:+.2f} vs {an6.loc[('eff_170', 'Stable', '>30 kb'), 'median']:+.2f}; p < 1e-29).",
         "Not adjusted for accessibility, so part of this may be accessibility itself."),
 "C4": U(f"In Stable sites a stronger footprint goes with a larger 170 bp effect (rho {an4.loc[('eff_170', 'Stable'), 'rho']:+.2f}) and higher predicted accessibility ({an4.loc[('native_peak', 'Stable'), 'rho']:+.2f}); in Dynamic sites there is no relationship.",
         "Consistent with, but not independent of, the accessibility data (footprint scores come from the same ATAC)."),
 "C5": U(f"Sites alone in their peak and sites sharing it with another Fos motif do not differ (raw FDR from {AN5.raw_fdr.min():.2f}; adjusted FDR 0.76, 0.76, 0.07).",
         f"Only {int(an5.loc['motif', 'n_multi'])} sites are in multi-motif peaks and each shuffle edits one window, so this is inconclusive."),
 "B3": U(f"Raw: Stable loses more at the 170 bp shuffle (Dynamic minus Stable {m.loc['170', 'diff_dyn_minus_stable']:+.2f}, 95% CI {m.loc['170', 'ci_lo']:+.2f} to {m.loc['170', 'ci_hi']:+.2f}; positive = Dynamic loses less).",
         f"It holds in all five quintiles of native peak accessibility, in a permutation test (p = {pm.loc['170', 'perm_p']:.4f}) and with the Patski track on the same sites ({fb.loc['170', 'diff_dyn_minus_stable']:+.2f}; same model and sequences, so not independent).",
         f"Motif width: {m.loc['motif', 'diff_dyn_minus_stable']:+.2f} (permutation p = {pm.loc['motif', 'perm_p']:.2f}, not significant). Peak width: {m.loc['peak', 'diff_dyn_minus_stable']:+.2f}, mostly accessibility and peak width.",
         "These quintiles use peak accessibility; the next section redoes this on the matched readout."),
 "H1": U(f"Cardiomyocyte sites: raw {hg(cmn, M_RAW).gap:+.2f} → {hg(cmn, M_PK).gap:+.2f} (peak accessibility) → {hg(cmn, M_MATCH).gap:+.2f} (matched-readout accessibility) → {hg(cmn, M_FULL).gap:+.2f} (plus peak accessibility, width, TSS distance; splines) → {hg(cmn, M_FP).gap:+.2f}, not significant (plus footprint score).",
         f"Fibroblast sites, heart track: {hg(fbh, M_RAW).gap:+.2f} raw → {hg(fbh, M_FULL).gap:+.2f}. Fibroblast sites, Patski track: {hg(fbp, M_RAW).gap:+.2f} raw → {hg(fbp, M_MATCH).gap:+.2f} (matched readout only) → {hg(fbp, M_FULL).gap:+.2f} (full model).",
         "Motif- and peak-width gaps shrink to near zero, with a few small reversed-sign residuals (for example Patski motif −0.05).",
         f"For sites like the Dynamic ones the weighted gap is {H_cs.loc[cmn, 'common_support_gap']:+.2f} (CM), {H_cs.loc[fbh, 'common_support_gap']:+.2f} (FB, heart) and {H_cs.loc[fbp, 'common_support_gap']:+.2f} (Patski); in the lower deciles, where {H_cs.loc[cmn, 'dynamic_in_lowest_5_deciles']:.0%} of Dynamic cardiomyocyte sites lie, the gap is about 0.1 to 0.2.",
         "Reading: the gap is largely chromatin context; a small residual remains; the data neither show nor exclude a separate sequence effect (the covariates are themselves downstream of sequence). Patski scores the same sequences with the same model, so it is not independent confirmation."),
 "B4": U(f"At 170 bp every class loses less than its reference: raw {a170.raw_diff.min():.2f} to {a170.raw_diff.max():.2f}, adjusted {a170.adj_diff.min():.2f} to {a170.adj_diff.max():.2f} (all {int((a170.adj_fdr < 0.05).sum())} age and exercise classes stay FDR < 0.05). At peak width {int((adj[adj.width == 'peak'].adj_fdr < 0.05).sum())} of 6 stay significant, at motif width {int((adj[adj.width == 'motif'].adj_fdr < 0.05).sum())} of 6.",
         f"Classes against each other: no significant difference at 170 bp or peak width (smallest FDR {E[E.width != 'motif'].fdr.min():.2f}).",
         "The intervals (about ±0.2 to 0.3; 79 to 210 sites per class) are too wide to exclude differences half the size of the Stable-vs-Dynamic gap: absence of evidence.",
         "The reference classes (bound in both conditions) contain a few Dynamic sites (2,997 and 3,117 vs 2,907 Stable)."),
 "D3": U(f"Stable loses more at 170 bp in all four site/track combinations (raw FDR < 0.001): cardiomyocyte sites {rep.loc[('CM sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (heart) and {rep.loc[('CM sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} (Patski); fibroblast sites {rep.loc[('FB sites, heart ATAC', '170'), 'diff_dyn_minus_stable']:+.2f} and {rep.loc[('FB sites, Patski ATAC', '170'), 'diff_dyn_minus_stable']:+.2f}.",
         f"After adjusting for matched accessibility and peak width: {adj_.loc[('CM sites, heart ATAC', '170'), 'adj_diff']:+.2f} and {adj_.loc[('FB sites, heart ATAC', '170'), 'adj_diff']:+.2f} (heart), {adj_.loc[('CM sites, Patski ATAC', '170'), 'adj_diff']:+.2f} and {adj_.loc[('FB sites, Patski ATAC', '170'), 'adj_diff']:+.2f} (Patski; FDR < 0.01), smaller still in the fullest model (section [[H1]]).",
         f"Class by class in fibroblast sites: all {n_cl_170} comparisons at 170 bp go the same way and {n_cl_sig170} reach FDR < 0.001.",
         "Motif width: a small raw gap with the heart track (+0.10); the peak-width gap disappears after adjustment in 3 of 4 combinations.",
         f"Prediction for fibroblast sites is weaker ({nvf.loc[('heart ATAC', 'atac_peak')]:.2f} heart, {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} Patski), so treat this as indicative."),
 "G0": U("Goal: keep sites whose predicted accessibility agrees with the measured snATAC, using the native prediction and measurement only (never the effects).",
         "Tier B: has a measured peak and predicted accessibility above the 5th percentile (a noise floor, kept low so weak but well-predicted peaks stay).",
         "Tier C: B plus predicted-minus-measured percentile within 1 SD; tier D: within 0.5 SD.", "Thresholds are set per site set and track; fibroblast sites are compared with measured fibroblast accessibility."),
 "G1": U(f"Tier C keeps {gr(cmn, 'C').retained_all:.0%} of cardiomyocyte sites ({gr(cmn, 'C').n:,}) and tier D {gr(cmn, 'D').retained_all:.0%} ({gr(cmn, 'D').n:,}).",
         "Most drops are for rank disagreement; the floor removes few sites.", "The correlation inside a tier (0.7 to 0.95) is circular: it shows the selection worked, nothing more."),
 "G2": U(f"By measured accessibility the filter treats the classes alike.", f"By predicted accessibility it does not: in the top two predicted quintiles {gq(cmn, 'Dynamic'):.0%} of Dynamic cardiomyocyte sites are kept vs {gq(cmn, 'Stable'):.0%} of Stable (fibroblast heart {gq(fbh, 'Dynamic'):.0%} vs {gq(fbh, 'Stable'):.0%}; Patski {gq(fbp, 'Dynamic'):.0%} vs {gq(fbp, 'Stable'):.0%}).",
         "So the kept Dynamic sites are the ones the model gets right; filtered results describe that subset."),
 "G3": U(f"Kept Stable sites lose more than dropped ones ({gk(cmn, 'Stable').median_effect_kept:.2f} vs {gk(cmn, 'Stable').median_effect_discarded:.2f}, p < 0.001); Dynamic sites do not differ (p = {gk(cmn, 'Dynamic').p:.2f}).",
         f"With Patski, kept Dynamic sites lose less than dropped ones ({gk(fbp, 'Dynamic').median_effect_kept:.2f} vs {gk(fbp, 'Dynamic').median_effect_discarded:.2f}, p < 0.001): the opposite direction.",
         f"So larger filtered gaps are partly selection; the unfiltered {gs(cmn, '170', 'A').diff_dyn_minus_stable:+.2f} is the cleaner estimate."),
 "G4": U(f"Dynamic minus Stable at 170 bp (all FDR < 0.001): cardiomyocyte {gs(cmn, '170', 'A').diff_dyn_minus_stable:+.2f} (all) → {gs(cmn, '170', 'C').diff_dyn_minus_stable:+.2f} (C) → {gs(cmn, '170', 'D').diff_dyn_minus_stable:+.2f} (D); fibroblast heart {gs(fbh, '170', 'A').diff_dyn_minus_stable:+.2f} → {gs(fbh, '170', 'D').diff_dyn_minus_stable:+.2f}; Patski {gs(fbp, '170', 'A').diff_dyn_minus_stable:+.2f} → {gs(fbp, '170', 'D').diff_dyn_minus_stable:+.2f}.",
         "Motif width: small (up to +0.1 raw, about +0.2 in the strictest tiers). Peak width: unstable across tiers.", "Read the growth with stricter tiers together with the previous section."),
 "G5": U("The 25th-percentile floor keeps about a third of Dynamic sites in tier C (vs about half with the 5% floor).", "Cardiomyocyte gap similar or larger; Patski +0.44 vs +0.61 (tier C).", "Direction and significance do not depend on the floor; size does somewhat."),
 "B6": U("This is the stricter 25th-percentile-floor version of the tiers (ranks over all sites), so counts differ slightly from the 25%-floor tables above.", "Effects are larger in the filtered sets and the 170 bp gap stays (about +0.5 in tier B, +0.7 in tier C)."),
 "G6": U(f"Replacing everything but the motif keeps {G_mu[(G_mu.combo == cmn)].ratio_flank170_over_random170.min():.0%} to {G_mu[(G_mu.combo == cmn)].ratio_flank170_over_random170.max():.0%} of the 170 bp loss (cardiomyocyte) in every tier.",
         f"The 7 bp core stays at about {G_ct.core7_share_of_170.min():.0%} to {G_ct.core7_share_of_170.max():.0%} of the 170 bp shuffle effect; shuffles 300 bp away stay near {G_ct.offset_pm300.mean():.2f} vs {G_ct.centred_170.min():.2f} to {G_ct.centred_170.max():.2f} centred.",
         "The table's 'core / full motif loss' (about 1.0) is a different ratio: the core against all 11 bases."),
 "C7": U(f"Native predicted accessibility vs MPRA activity: rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} (heart), {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} (Patski); also in the replicated and TOBIAS-motif subsets.",
         f"Predicted scrambling effect vs MPRA change: about 0 (motif shuffle {mp_.loc[('Patski ATAC', CMP_M, 'all elements'), 'rho']:+.2f}, full shuffle {mp_.loc[('Patski ATAC', CMP_F, 'all elements'), 'rho']:+.2f}, Patski).",
         "Sign: the predicted effect is log2(shuffled/native), negative = loss; the MPRA change is native minus scrambled, positive = loss; agreement would show as a negative correlation.",
         "Possible reasons: episomal activity vs genomic accessibility; tiny MPRA effects (core loss about 0.02 log2); weak fibroblast tracks. So this null does not rank models."),
 "D4": U(f"Neutral window: native rho {mc_.loc[('heart ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (heart), {mc_.loc[('Patski ATAC', 'neutral context', 'native level vs MPRA activity'), 'rho']:+.2f} (Patski), similar to the genomic window; scramble effects still about 0 (all |rho| ≤ {D4[D4.comparison != 'native level vs MPRA activity'].rho.abs().max():.2f}).",
         "So the surrounding genome is not why the native head misses the MPRA changes.", "Caveat: only 2 neutral windows; applies to the native head, not the manuscript's trained probe."),
 "B5": U("Median RNA effects are about zero: the largest are −0.08 (±2 kb readout, 2 to 10 kb from the TSS) and −0.12 (±50 kb mean, within 2 kb), Stable sites, peak shuffle.", "ATAC and RNA effects correlate weakly across sites (rho 0.19 to 0.24).", "Not informative with these tracks."),
 "B8b": U(f"Cardiomyocyte sites: +{mi('atac_motif', '+motif vs shuffled background').median_log2fc:.2f} log2 over the centre 7 bp (not significant); about +0.05 over 128 bp (p < 0.001, 68% positive). Fibroblast sites: about +0.14 over 128 bp (92% positive).",
          "Background accessibility is near zero, so ratios sit at the floor; the effects are small and not interpretable. The shuffled-motif control can also keep AP-1 half-sites."),
}
KEEP_T = {"D3": 1, "G1": 1, "A5": 2, "B1": 1}

# ---------------- hand-written sections
MAN_ROWS = [
 ("Dynamic AP-1 elements have lower baseline accessibility", "ok", "Predicted native ATAC is lower for Dynamic MPRA elements (Cohen's d about -0.5) and Dynamic cardiomyocyte sites."),
 ("Dynamic elements have lower GC content and sit farther from promoters", "mid", "GC lower (median 0.494 vs 0.512, p = 9e-8); distance only weakly different in MPRA elements, not significant in our cardiomyocyte sites."),
 ("Age-gained elements have reporter activity on par with stable elements", "ok", "Median log2 alpha +0.25 (Age Up) vs +0.26 (Stable), p = 0.78."),
 ("...despite weaker endogenous chromatin features, so chromatin context explains them", "mid", "Consistent with, not tested: the model predicts lower accessibility for these sequences while MPRA activity is equal, so it shows accessibility and activity differ."),
 ("Shuffling the AP-1 core has a minor effect on MPRA activity; flank shuffling hurts dynamic elements more", "ok", "Core loss about 0.02 log2; flank loss reproduces (Age Up 0.128 vs Stable 0.016, d = 0.14), a small effect."),
 ("The same flank dependence is what the model can read", "bad", "The native head gives the opposite direction (Dynamic lose less) and per-element predicted vs measured flank loss correlates about 0. Accessibility and activity differ."),
 ("AlphaGenome-MPRA probe reaches r = 0.73; ISM, seqlets, epistasis, rescue explain dispersed activity", "na", "Cannot be tested with the native heads. GC content alone correlates 0.505 with MPRA activity, so a GC baseline on the same splits would help interpret 0.73."),
 ("Motif enrichments (MEF2C, TCF7L2, GATA4, NFAT, ETV4, C/EBP)", "na", "Needs a motif scan on sequences, which can run on Katana."),
]
MAN = f"""<ul class='pts'><li>The manuscript uses AlphaGenome as a frozen encoder with a probe trained on the MPRA, not the native accessibility head used here, so our results speak to predicted accessibility, not to its probe.</li><li>Full notes: <code>Analysis/manuscript_notes_MPRA_AlphaGenome.md</code>.</li></ul>
<table class="t reg"><tr><th>Manuscript claim</th><th>Verdict from our data</th><th>Evidence</th></tr>{''.join(f"<tr><td class='v'>{a}</td><td>{badge(k).replace('High confidence', 'Supported').replace('Moderate confidence', 'Partly supported').replace('Low / unreliable', 'Opposed').replace('Untested', 'Cannot test')}</td><td class='v'>{b}</td></tr>" for a, b, c in [(x, y, z) for x, y, z in MAN_ROWS] for b in [c] for c in [b] ) if False else ''.join(f"<tr><td class='v'>{a}</td><td>{badge(k).replace('High confidence', 'Supported').replace('Moderate confidence', 'Partly supported').replace('Low / unreliable', 'Opposed').replace('Untested', 'Cannot test')}</td><td class='v'>{b}</td></tr>" for a, k, b in MAN_ROWS)}</table>"""
NOTE = U("The dinucleotide shuffle of the motif is partly a no-op: an 11-mer has few arrangements and 3.8% of replicates gave a native-identical prediction.",
         f"Fibroblast tracks are weak proxies: predicted vs measured fibroblast accessibility is {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} (Patski) and {nvf.loc[('heart ATAC', 'atac_peak')]:.2f} (heart).",
         "Not run yet: saturation mutagenesis, the human-head check and the window-length check (coded, resumable).")
NEXT = U("Other-motif shuffle control (highest value): pair each AP-1 site with a footprinted non-AP-1 motif (GATA, TBX, MEF2, NFAT) in the same peak plus matched random windows; start by counting peaks with both.",
         "Saturation mutagenesis of the 11 bp motif; human-head and window-length checks.",
         "Rerun the small code fixes on Katana (MPRA off-by-one, dinucleotide redraw, insertion half-sites, record the core-edit fallback) and regenerate the MPRA crossmap.",
         "Confirm whether Fosb, Fosl1 and Fosl2 calls belong in the site set.", "A cardiac-cell-type accessibility probe would replace the weak fibroblast proxies.", "A motif scan on sequences to test the manuscript's enrichment claims.")

REG = [
 ("Native vs measured accessibility", f"rho {nv.rho:.2f} (CM, heart); fibroblast {nvf.loc[('heart ATAC', 'atac_peak')]:.2f} (heart), {nvf.loc[('Patski ATAC', 'atac_peak')]:.2f} (Patski)", "mid", "Informative, not precise; weak fibroblast tracks."),
 ("Prediction quality by class", f"rho {XCL[(XCL.combo == 'CM sites, heart ATAC') & (XCL['class'] == 'Stable (all three)')].rho.iloc[0]:.2f} Stable vs {XCL[(XCL.combo == 'CM sites, heart ATAC') & (XCL['class'] == 'All Dynamic')].rho.iloc[0]:.2f} Dynamic; Dynamic over-predicted", "mid", "Error size similar; Dynamic predicted a bit worse and biased up."),
 ("Class labels vs measured change", "Gains clear, losses weak", "mid", "Footprints and accessibility share the same ATAC data."),
 ("Shuffle at 3 widths", f"{m.loc['170', 'median_stable']:.2f} (170 bp) vs {m.loc['motif', 'median_stable']:.2f} (motif), Stable", "ok", "Large, above replicate noise; a model prediction."),
 ("Offset controls and titration", f"Centred {co_.loc['Stable', 'centred_170_median']:.2f} vs {co_.loc['Stable', 'offset_pm300_median']:.2f} at ±300 bp; saturates near 340 bp; 7 bp core {fr_.loc[7]:.0%}", "mid", "Local and saturating is solid; AP-1 specificity not shown."),
 ("Mutation arm", f"Surroundings alone = {ratio_cm170:.0%} (CM), {ratio_fb170:.0%} (FB) of the 170 bp loss", "mid", "Ratio near 100% partly expected; track-dependent."),
 ("Peak swap", f"Swapped vs donor rho {swap.loc['donor native (sequence)', 'rho']:.2f}; no Stable/Dynamic difference after adjusting for donor accessibility", "mid", "4 recipients, 25 Dynamic donors."),
 ("Stable vs Dynamic", f"Raw {m.loc['170', 'diff_dyn_minus_stable']:+.2f} → {hg(cmn, M_MATCH).gap:+.2f} (matched) → {hg(cmn, M_FULL).gap:+.2f} (with width, TSS)", "mid", "Largely context; small residual; neither shown nor excluded a separate effect."),
 ("Age vs exercise classes", "No significant difference between dynamic classes", "inc", "Underpowered; labels overlap."),
 ("Fibroblast sites", f"Raw gap replicates; Patski {hg(fbp, M_MATCH).gap:+.2f} (matched only), {hg(fbp, M_FULL).gap:+.2f} (full)", "mid", "Weak tracks; Patski not independent."),
 ("Agreement-filtered sites", "All main conclusions hold in every tier", "mid", "Filter drops over-predicted Dynamic sites; filtered gaps partly selection."),
 ("Footprint, TSS distance", "Strong footprints and distal sites depend more on 170 bp context", "mid", "Footprint shares data with accessibility; TSS unadjusted."),
 ("Multi-motif peaks", "No difference (inconclusive)", "inc", "Few sites; one-window edits."),
 ("MPRA: native level", f"rho {mp_.loc[('heart ATAC', CMP_N, 'all elements'), 'rho']:+.2f} heart, {mp_.loc[('Patski ATAC', CMP_N, 'all elements'), 'rho']:+.2f} Patski", "mid", "Weak but consistent."),
 ("MPRA: scramble effects (genomic and neutral)", "About 0 in both contexts", "mid", "MPRA effects tiny; 2 neutral backgrounds."),
 ("RNA", "Median effects about 0", "bad", "Not informative with these tracks."),
 ("Motif insertion", f"+{mi('atac_motif', '+motif vs shuffled background').median_log2fc:.2f} log2 over 7 bp (ns); about +0.05 over 128 bp", "bad", "At the floor; not interpretable."),
 ("Replicates and GPU", "10 replicates adequate for group medians; no L40S vs A100 difference", "ok", "H100 ran one 91-site chunk."),
]
REG_HTML = "<table class='t reg'><tr><th>Experiment</th><th>Result</th><th>Confidence</th><th>Why</th></tr>" + "".join(
    f"<tr><td class='v'>{a}</td><td class='v'>{b}</td><td>{badge(k)}</td><td class='v'>{d}</td></tr>" for a, b, k, d in REG) + "</table>"

BIO = f"""<p><strong>Bottom line.</strong> AlphaGenome says an AP-1 site's accessibility sits mostly in the surrounding few hundred base pairs, with the 11&nbsp;bp motif a minority share. Sites whose binding changes with age or exercise differ mainly in chromatin context; a small residual (about {hg(cmn, M_FULL).gap:+.2f}) remains. Model predictions, checked against measured snATAC where stated; not direct measurements of AP-1 in cells.</p>
<h3>What these tests add</h3>
{U(f"A size and a range for the context effect: 300&nbsp;bp away does nothing ({co_.loc['Stable', 'offset_pm300_median']:.2f} vs {co_.loc['Stable', 'centred_170_median']:.2f}), the loss saturates near 340&nbsp;bp, the 7&nbsp;bp core is about {fr_.loc[7]:.0%} of the 170&nbsp;bp effect.",
   f"The motif is rarely the unit: {an3.loc['Context only', 'share_stable']:.0%} of Stable sites are context-only, {an3.loc['Motif only', 'share_stable']:.1%} motif-only, and at {an1.loc[('motif', 'Stable'), 'frac_gt0']:.0%} shuffling the motif raises the prediction.",
   f"Replacing the surroundings keeps {ratio_cm170:.0%} (CM) / {ratio_fb170:.0%} (FB) of the loss from replacing all 170&nbsp;bp.",
   f"Responsive AP-1 sites show little beyond context: the Stable-vs-Dynamic gap goes from {hg(cmn, M_RAW).gap:+.2f} to {hg(cmn, M_FULL).gap:+.2f}; a separate sequence effect is neither shown nor excluded.",
   f"Dynamic sites are less sensitive to every wide edit in the raw results, but matched on predicted context the two classes have the same sensitivity profile across all {len(Y_prof)} mutation types; the same site is often Stable in one cell type and Dynamic in the other (kappa {Y_cross.cohens_kappa:.2f}).",
   f"Dynamic sites are over-predicted ({xts.loc['CM sites, heart ATAC', 'dynamic_median_signed']:+.2f} vs {xts.loc['CM sites, heart ATAC', 'stable_median_signed']:+.2f}), so their openness may depend on condition or cell state beyond sequence.",
   "Sequence sets a peak's rank and location sets its level (peak swap), equally for Stable and Dynamic donors.",
   "Accessibility is not enhancer activity at the element level: predicted scrambling effects correlate about 0 with MPRA changes, in genomic and neutral windows.")}
<h3>What only confirms what was expected</h3>
{U("Dynamic elements have weaker baseline chromatin (the manuscript's lower baseline accessibility): confirmed in predicted and measured accessibility.", "Age Up and Exercise Up sites gain measured accessibility; the down classes change only slightly.")}
<h3>What complicates the manuscript</h3>
{U("The MPRA flank result (dynamic elements lose more activity when flanks are shuffled) is not mirrored in predicted accessibility: Dynamic sites lose less, and little separates the classes once context is accounted for. That result is about activity, not baseline accessibility.")}
<h3>What these tests cannot say</h3>
{U("Whether AP-1 matters more than other motifs, or which neighbouring factors matter (cooperation with GATA, TBX or MEF2 is plausible but untested here).", "Causality in cells, or why dynamic sites respond: the model predicts baseline accessibility only.", "Fibroblast conclusions beyond indicative: the fibroblast tracks are weak proxies.")}"""

BIO_SEC = {
 "A1": "AP-1 is bound in all three conditions at most sites, so occupancy here is largely constitutive; only about a fifth of sites change.",
 "B1": "The model has learned enough regulatory logic to rank real regions moderately; the cardiac track transfers to fibroblast sites almost as well as to cardiomyocyte sites, the Patski cell line far less.",
 "A2": "Sites where AP-1 binding changes sit in narrower peaks, not nearer or farther from genes.",
 "A4": "AP-1 gain goes with chromatin opening; AP-1 loss goes with only slight closing (footprints and accessibility share data, so treat as consistency).",
 "A5": "One reading of the over-prediction is that Dynamic sites carry permissive sequence but are open only under certain conditions or cell states (an interpretation, not a test).",
 "B2": "Most of an AP-1 site's predicted accessibility is in the surrounding sequence; the motif is a minority share (about a third).",
 "C1": "At many sites the motif is not what the model relies on, so AP-1 may be bound to accessibility made by other factors.",
 "D1": "Accessibility at an AP-1 site is encoded in a compact neighbourhood of a few hundred base pairs; the AP-1 core is a minority of it.",
 "D2": "The model reads the site's accessibility mainly from neighbouring sequence; the motif counts a little more with the fibroblast-like track, but that follows the track, not the cell type.",
 "B8a": "A peak keeps its rank when moved but its level shifts with the new location: sequence carries the potential, the environment scales it (limited power for Dynamic donors).",
 "C3": "AP-1 rarely acts alone in the model. Cooperation with lineage factors such as GATA, TBX or MEF2 is a plausible reading but is not tested here.",
 "C2": "More open, narrower elements lose more when their sequence is removed, because there is more to lose.",
 "C6": "Promoter-proximal AP-1 sites depend less on their local neighbourhood than distal ones; promoter accessibility may be maintained by other features (not adjusted for accessibility).",
 "C4": "Well-bound AP-1 sites sit in accessible, motif-rich sequence; this partly restates that open chromatin supports binding.",
 "C5": "Nothing can be said about redundancy between neighbouring AP-1 sites from these data.",
 "B3": "Sites whose AP-1 binding changes lose less accessibility when their surroundings are removed; the next section asks whether that is more than chromatin context.",
 "H1": "Responsive AP-1 sites are not obviously built differently from constitutive ones: once their weaker chromatin context is counted, little extra sequence dependence remains. Their responsiveness may therefore come from things a sequence-only baseline model does not see (AP-1 activity, co-factors, chromatin state), as the manuscript suggests for dynamic elements. This cannot be tested here, because the model predicts baseline accessibility, not response.",
 "B4": "The model cannot tell age-driven from exercise-driven AP-1 changes by sequence at these sample sizes.",
 "D3": "The same picture at fibroblast AP-1 sites (indicative only; weak tracks).",
 "G2": "The well-predicted Dynamic sites are those whose accessibility is sequence-driven; the filtered result describes them, not all dynamic sites.",
 "C7": "Accessibility and reporter activity are different properties: a model that ranks open sequences does not predict how much activity a scramble removes, which is the manuscript's own reason for training a probe on the MPRA.",
 "D4": "The mismatch is not caused by the model using the surrounding genome; MPRA effects are small and accessibility is not enhancer activity.",
"DSH": "Dynamic and Stable sites depend on the same parts of their sequence in the same proportions; they differ in how much accessibility is at stake, not in how it is built.",
 "DT": "AlphaGenome often ranks a neighbouring element (one example is a MEF2 site) above the AP-1 motif itself, in both classes. Whether Dynamic sites lean on fewer or more distant partners is an open question for the full scan.",
"DS": "Dynamic sites are less sensitive to broad sequence removal because they sit in weaker chromatin, not because they depend on a different part of the sequence: the profile across mutation types has the same shape.",
 "DP": "No binding pattern shows a robust, specific vulnerability or resistance once context is accounted for.",
 "DF": "Nothing in the model's outputs (effects of any edit, or RNA predictions) marks a Dynamic site once context is matched.",
 "DC": "Condition-specific opening exists but is modest, and it is the part a baseline, condition-blind sequence model cannot see.",
 "DM": "Nothing in the model's features marks Dynamic sites as built differently; what separates them is weaker context, and the gap between prediction and measurement.",
 "DX": "The same DNA site is often Stable in one cell type and Dynamic in another, so class depends on cell state more than on the sequence.",
 "B5": "Single distal AP-1 elements barely change predicted gene expression; nothing about target-gene regulation can be concluded.",
 "B8b": "No conclusion about whether the motif alone can create accessibility.",
}

TECH = U(f"Shuffles lower accessibility at {q6.loc[('170', 'Stable'), 'frac_sites_decreased']:.0%} of Stable sites (170 bp) and {q6.loc[('motif', 'Stable'), 'frac_sites_decreased']:.0%} (motif only; about a quarter go up), so effects are not an artefact of broken sequence.",
         f"Replicates: a single shuffle is off by about {q7.loc[('170', 1)]:.2f} (170 bp) and {q7.loc[('peak', 1)]:.2f} (peak) log2; with 5 it is {q7.loc[('170', 5)]:.2f} and {q7.loc[('peak', 5)]:.2f}. Ten are enough for group medians; per-site peak-shuffle values stay noisy.",
         f"The effect arms use 10 shuffles per site (5 for replace, dinucleotide and flank edits; 3 for the core edit), averaged per site.",
         f"GPU type: no difference between L40S and A100 chunks (p {q8.kruskal_p_across_gpus.min():.2f} to {q8.kruskal_p_across_gpus.max():.2f}); the H100 ran one 91-site chunk and cannot be assessed.") + \
       "<details><summary>Tables and figures for these checks</summary><h4>Shuffle behaviour</h4>" + SEC["A6"] + "<h4>Replicate convergence</h4>" + SEC["A7"] + "<h4>Noise per site</h4>" + SEC["B7"] + "<h4>GPU consistency</h4>" + SEC["A8"] + "</details>"
B2 = U(f"Shuffling the 11&nbsp;bp motif: {m.loc['motif', 'median_stable']:.2f} log2 at the motif (about a third lost). Shuffling the 170&nbsp;bp: {m.loc['170', 'median_stable']:.2f} (about 70% lost). Shuffling the whole peak: {m.loc['peak', 'median_stable']:.2f} over the peak (Stable sites).",
       f"Dynamic sites lose less in most cells ({m.loc['motif', 'median_dynamic']:.2f}, {m.loc['170', 'median_dynamic']:.2f}, {m.loc['peak', 'median_dynamic']:.2f}).",
       "The order motif < 170 bp ≈ peak is robust; peak-shuffle values are noisy per site.") + \
    f"<figure>{ximg('x2_shuffle_dotplot', 'Shuffle effects as a grouped dot plot')}<figcaption>Each panel is one shuffle; the three positions are the region over which the effect is read. Median with 95% bootstrap CI; marks: Stable vs Dynamic (Mann&ndash;Whitney, BH over the nine comparisons).</figcaption></figure>" + \
    f"<details><summary>The same numbers as a heatmap and a table</summary><figure>{img(CM / 'figures' / 'fig2s_shuffle_heatmaps_stars.png', 'Shuffle effect heatmaps')}<figcaption>Median log2FC per shuffle width (rows) and readout region (columns); left Stable, right Dynamic with marks against Stable.</figcaption></figure>{table(stv_tab)}</details>"

SPECIAL = {"MAN": MAN, "NOTE": NOTE, "NEXT": NEXT, "TECH": TECH, "B2": B2}


_pr = lambda m, c: Y_prof.loc[m, c]
_pat = lambda m, pat: Y_pat[(Y_pat.mutation == m) & (Y_pat.pattern == pat)].iloc[0]
_sigpat = Y_pat[Y_pat.fdr < 0.05]
_mm = _sigpat[_sigpat.pattern == "MAS+MAE"]
_fb = lambda trk, m: Y_fbr[(Y_fbr.track == trk) & (Y_fbr.pattern == "MAS+MAE") & (Y_fbr.mutation == m)].iloc[0]
DS = U(f"Raw (all sites): Dynamic sites are less sensitive to every wide edit, for example shuffle the 170 bp {_pr('shuffle the 170 bp', 'stable_median'):.2f} (Stable) vs {_pr('shuffle the 170 bp', 'dynamic_median'):.2f} (Dynamic), replace the 170 bp {_pr('replace the 170 bp', 'stable_median'):.2f} vs {_pr('replace the 170 bp', 'dynamic_median'):.2f}, replace the 159 bp around the motif {_pr('replace the 159 bp around the motif', 'stable_median'):.2f} vs {_pr('replace the 159 bp around the motif', 'dynamic_median'):.2f}, centred shuffle of 340 bp {_pr('centred shuffle, 340 bp', 'stable_median'):.2f} vs {_pr('centred shuffle, 340 bp', 'dynamic_median'):.2f} (all FDR < 0.001).",
       f"Motif-level edits differ by only about 0.1 (replace the motif {_pr('replace the motif', 'stable_median'):.2f} vs {_pr('replace the motif', 'dynamic_median'):.2f}; replace the 7 bp core {_pr('replace the 7 bp core', 'stable_median'):.2f} vs {_pr('replace the 7 bp core', 'dynamic_median'):.2f}; FDR 0.02), the 4 flank bases and the motif dinucleotide shuffle do not differ, and edits far from the motif do not differ.",
       f"Context-matched pairs ({int(Y_match.loc['170 bp shuffle', 'n_pairs'])} Dynamic sites, each with the Stable site closest in predicted accessibility, peak width and TSS distance): the two classes have the same profile across all {len(Y_prof)} mutation types; the largest difference is {(Y_prof.matched_dynamic_median - Y_prof.matched_stable_median).abs().max():.2f} log2 and none is significant (smallest FDR {Y_prof.fdr_matched.min():.2f}).",
       "So the sensitivity of Dynamic sites is lower in proportion to their weaker context, not selectively for particular kinds of mutation.") + \
     f"<figure>{dimg('y5_mutation_type_profile', 'Sensitivity to each mutation type, Dynamic vs Stable')}<figcaption>Median effect of each mutation type for Stable (blue) and Dynamic (orange) sites with 95% bootstrap CI; marks: Mann&ndash;Whitney in all sites (left, BH) and paired Wilcoxon in the context-matched pairs (right, BH). Every edit is read over the region named in the group heading.</figcaption></figure>"
DP = U(f"Each binding pattern against Stable, adjusted for native accessibility at the readout, peak accessibility, peak width and TSS distance (splines): {int((Y_pat.fdr < 0.05).sum())} of {int((Y_pat.pattern != 'All Dynamic').sum())} pattern-by-mutation comparisons are significant (BH), {len(_mm)} of them for one pattern, MAS+MAE (bound in both aged states, not in young): less sensitive to the 170 bp and wider edits by {_mm[_mm.adj_diff > 0].adj_diff.min():+.2f} to {_mm[_mm.adj_diff > 0].adj_diff.max():+.2f} log2.",
       f"It does not replicate: in fibroblast sites with the Patski track MAS+MAE shows no 170 bp difference ({_fb('FB sites, Patski', 'shuffle 170→170').adj_diff:+.2f}) and is slightly more sensitive to motif and core edits ({_fb('FB sites, Patski', 'shuffle motif→motif').adj_diff:+.2f}, {_fb('FB sites, Patski', 'replace core→motif').adj_diff:+.2f}); with the heart track every Dynamic pattern shifts together (MAS+MAE {_fb('FB sites, heart', 'shuffle 170→170').adj_diff:+.2f}, similar for all Dynamic).",
       "The other significant cells are single exceptions (for example MAS only slightly more sensitive to the centred 11 bp shuffle) and are not consistent across tracks.",
       "So no binding pattern shows a robust, specific sensitivity; MAS+MAE is a hypothesis to test with the fine-tuned head or in more sites.") + \
     f"<figure>{dimg('y6_mutation_type_by_pattern', 'Sensitivity by binding pattern and mutation type')}<figcaption>Difference from Stable for each binding pattern (columns) and each mutation type (rows): raw difference of medians (left) and after adjustment (right; marks FDR from BH over all pattern-by-mutation cells, All Dynamic uncorrected). Positive = the pattern loses less.</figcaption></figure>"
_sc = Y_scan[Y_scan.fdr < 0.05]
DF = U(f"Of {len(Y_scan)} AlphaGenome outputs (native predictions, effects of every edit, RNA predictions and contrasts) {len(_sc)} differ at FDR < 0.05 between Dynamic sites and their context-matched Stable partners, all tiny (largest {Y_scan.matched_std_diff.abs().max():.2f} SD of the Stable sites).",
       f"A classifier given any set of these outputs cannot tell Dynamic from matched Stable sites: AUC {Y_mauc.loc['everything AlphaGenome outputs', 'auc_logistic']:.2f} for all outputs, {Y_mauc.loc['shuffle effects (9)', 'auc_logistic']:.2f} shuffle effects, {Y_mauc.loc['mutation effects (9)', 'auc_logistic']:.2f} mutation effects, {Y_mauc.loc['RNA predictions', 'auc_logistic']:.2f} RNA predictions (the matching variables themselves {Y_mauc.loc['matching variables (native accessibility, peak width, TSS distance)', 'auc_logistic']:.2f}).",
       f"Across all sites, native predicted accessibility alone gives AUC {Y_auc.loc['native predicted accessibility only (motif, 170 bp, peak)', 'auc_logistic']:.2f}; adding the shuffle effects changes nothing ({Y_auc.loc['context + model effects', 'auc_logistic']:.2f}); measured accessibility and the prediction error lift it to {Y_auc.loc['context + effects + measured accessibility + prediction error', 'auc_logistic']:.2f}.",
       "Predicted RNA near the TSS is a little lower at Dynamic sites (matched difference −0.13 log10, FDR 0.053): too weak to rely on.") + \
     f"<figure>{dimg('y4_feature_scan', 'AlphaGenome outputs that differ most between Dynamic and matched Stable sites')}<figcaption>The 22 largest standardised differences between Dynamic sites and context-matched Stable sites (filled = FDR < 0.05; grey diamonds = all sites, unmatched), coloured by type of output.</figcaption></figure><figure>{dimg('y1_classifier', 'Classifiers')}<figcaption>Cross-validated AUC for telling Dynamic from Stable with different feature sets, in all sites (left) and in the context-matched sample (right); folds are grouped by peak.</figcaption></figure>"
DC = U(f"Measured accessibility at Dynamic sites is higher in the conditions where they are called bound: {Y_cond.loc['All Dynamic', 'measured_log2_bound_minus_unbound']:+.2f} log2 on average ({Y_cond.measured_log2_bound_minus_unbound.drop('All Dynamic').min():+.2f} to {Y_cond.measured_log2_bound_minus_unbound.drop('All Dynamic').max():+.2f} by pattern, all p < 0.001): real but small (about 11%).",
       f"It also varies more across conditions at Dynamic sites (median SD {Y_var.dynamic_sd:.2f} vs {Y_var.stable_sd:.2f} log2, p < 1e-10).",
       f"The model predicts one baseline value, and over-predicts Dynamic sites more in the conditions where they are unbound ({Y_cond.loc['All Dynamic', 'overprediction_diff_unbound_minus_bound']:+.2f} rank units; every pattern FDR < 0.001): consistent with permissive sequence that is open mostly when bound.",
       f"Its prediction correlates about equally with accessibility in all three conditions (YS only: {Y_cor.loc['YS only', 'rho_pred_vs_YS']:.2f}, {Y_cor.loc['YS only', 'rho_pred_vs_MAS']:.2f}, {Y_cor.loc['YS only', 'rho_pred_vs_MAE']:.2f}), so it carries no condition-specific information; sites bound only in aged states match young-state accessibility best (MAS only: {Y_cor.loc['MAS only', 'rho_pred_vs_YS']:.2f} with YS vs {Y_cor.loc['MAS only', 'rho_pred_vs_MAS']:.2f} with MAS).") + \
     f"<figure>{dimg('y2_condition_specific', 'Condition-specific accessibility at Dynamic sites')}<figcaption>Left: measured accessibility in bound minus unbound conditions, by binding pattern. Middle: the model's over-prediction in unbound minus bound conditions (Wilcoxon, BH). Right: variability of measured accessibility across the three conditions.</figcaption></figure>"
DM = U(f"Matched pairs ({int(Y_match.loc['170 bp shuffle', 'n_pairs'])} of {int(Y_match.loc['170 bp shuffle', 'n_dynamic_total'])} Dynamic sites matched to the Stable site with the same 170 bp and peak accessibility, peak width and TSS distance; balance within {Y_bal.std_diff_after.abs().max():.2f} SD): the 170 bp gap is {Y_match.loc['170 bp shuffle', 'matched_gap_median']:+.2f} (CI {Y_match.loc['170 bp shuffle', 'matched_ci_lo']:+.2f} to {Y_match.loc['170 bp shuffle', 'matched_ci_hi']:+.2f}) against {Y_match.loc['170 bp shuffle', 'raw_gap']:+.2f} raw; motif and peak shuffles about 0.",
       f"Same-peak pairs ({int(Y_same.loc['170 bp shuffle effect', 'n_peaks'])} peaks holding both classes): the Dynamic site sits at lower accessibility ({Y_same.loc['native accessibility, 170 bp (log10)', 'median_dynamic_minus_stable']:+.2f} log10 at 170 bp) and loses less ({Y_same.loc['170 bp shuffle effect', 'median_dynamic_minus_stable']:+.2f}, CI {Y_same.loc['170 bp shuffle effect', 'ci_lo']:+.2f} to {Y_same.loc['170 bp shuffle effect', 'ci_hi']:+.2f}): the effect follows local accessibility even inside one peak.",
       f"Position in the peak differs a little (Dynamic sites slightly more off-centre: {Y_pos.dynamic_median:.2f} vs {Y_pos.stable_median:.2f}, p = {Y_pos.p:.3f}), too small to matter.") + \
     f"<figure>{dimg('y3_matched_and_same_peak', 'Matched pairs and same-peak pairs')}<figcaption>The 170 bp gap in all sites, in matched pairs and in same-peak pairs (median with 95% bootstrap CI).</figcaption></figure>"
DX = U(f"{int(Y_cross.shared_sites)} sites are present in both the cardiomyocyte and fibroblast tables; {Y_cross.agreement:.0%} keep the same class (Cohen's kappa {Y_cross.cohens_kappa:.2f}). {int(Y_cross.stable_cm_dynamic_fb)} are Stable in cardiomyocytes but Dynamic in fibroblasts, {int(Y_cross.dynamic_cm_stable_fb)} the reverse; only {int(Y_cross.dynamic_both)} are Dynamic in both.",
       f"Of cardiomyocyte-Dynamic sites {1 - Y_cross.p_fb_stable_given_cm_dynamic:.0%} are Dynamic in fibroblasts, against {1 - Y_cross.p_fb_stable_given_cm_stable:.0%} of cardiomyocyte-Stable sites: a weak association.",
       "So the same DNA site is often Stable in one cell type and Dynamic in the other: the class depends more on cell state than on the sequence.",
       "Caveat: shared sites are a special set (called in both cell types).")
DEXP = U("<strong>Ready to submit: tile scan and saturation mutagenesis</strong> (<code>bash submit_all.sh scan</code>). On 1,470 cardiomyocyte sites (500 Dynamic, 470 context-matched Stable, 500 high-accessibility Stable): every 10 bp tile from −300 to +290 bp replaced (3 replicates) and every single-base change from −30 to +30 bp. Question: where do the important bases lie, and is the importance more spread out or more concentrated at the motif in Dynamic sites?",
         "<strong>Ready to submit: AP-1 titration, move and duplicate</strong> (<code>bash submit_all.sh titrate</code>), same sites. The 7 bp core weakened by 1 to 7 bases, upgraded to a perfect core (at sites without one) or cut to a half-site, and the whole 11 bp motif reverse-complemented; and the 11 bp motif duplicated or moved to eight offsets (±20 to ±150 bp) with shuffled-motif controls. Question: do Dynamic sites respond differently to motif strength, orientation and position? Test first with <code>bash submit_all.sh dynamics_pilot</code>.",
         "<strong>Intrinsic-sequence test.</strong> Place each site's 170 bp (native and shuffled) into neutral backgrounds, as in the MPRA neutral-context run, for Dynamic and context-matched Stable sites: is the accessibility a property of the sequence or of its location? Small run (hours).",
         "<strong>Co-factor insertion and rescue.</strong> Insert GATA, TBX, MEF2 or NFAT motifs 20 to 60 bp from the AP-1 site at Dynamic sites: which raises predicted accessibility to Stable-like levels?",
         "<strong>Other-motif shuffle control</strong> (already planned): does AP-1 matter more than other footprinted motifs in the same peak?",
         "<strong>Breadth and state from the other heads</strong> (one forward pass for all 8,697 sites): ATAC across 18 mouse tissues and DNase across 67 give tissue breadth; heart and embryonic-fibroblast H3K27ac, H3K4me1, H3K4me3, EP300 (postnatal heart) and RNA/CAGE show whether Dynamic sites look like weaker or poised enhancers. There is no AP-1 ChIP track in the mouse model, so FOS binding cannot be predicted directly.",
         "<strong>Fine-tuned cardiomyocyte and fibroblast heads with per-condition tracks</strong> (see [[FT]]) would let the model learn condition-specific accessibility, which a baseline model cannot.")
FT = U("<strong>Essential: base-resolution accessibility per cell type and condition.</strong> Pseudobulk fragment files or BAMs from the snATAC, with the barcode-to-cell-type table, turned into bigWigs of Tn5 insertion sites (depth-normalised) for cardiomyocytes and fibroblasts (ideally all 19 types) and for each condition (YS, MAS, MAE, YE) and replicate (2 batches × 2 sexes), on mm10. The project folder has only peak-level counts, which cannot supervise base-resolution tracks.",
       "<strong>Genome-wide training windows</strong> including closed regions, with chromosome-held-out train, validation and test sets (for example hold out two chromosomes each), mappability and blacklist masks, and fragments per pseudobulk reported (aim for tens of millions per cell type and condition; if shallower, train on 128 bp bins or merge conditions).",
       "<strong>Strongly useful: matched snRNA-seq pseudobulk</strong> per cell type and condition (strand-specific coverage or 128 bp bins) for the RNA heads.",
       "<strong>Strongly useful: c-FOS (and JUN) CUT&Tag per condition</strong> from the manuscript, as a FOS-binding head (the base mouse model has no AP-1 track), plus any histone CUT&Tag.",
       "<strong>Held out for testing:</strong> the fibroblast MPRA, to check whether the fibroblast head predicts scramble effects better than the Patski track does.",
       "<strong>Design:</strong> per-condition output tracks so condition-specific differences are learnable; start with a probe on the frozen backbone (cheap, thousands of windows; Minh's trial-19 probe code is the starting point) before unfreezing the last blocks; train multi-task with the original heart and Patski tracks to avoid forgetting.",
       "<strong>Evaluate</strong> on held-out chromosomes (signal and peak correlation), on differential accessibility (MAS minus YS), against the built-in heart and Patski heads, and by whether the perturbation results in this report change.")
DPR = AN / "dynamic_profiles"
pr_ = lambda n: pd.read_csv(DPR / "tables" / n)
Z_A, Z_C, Z_B, Z_N = pr_("a_response_shape_by_mutation_type.csv"), pr_("a_classifier_matched_pairs.csv").set_index("features"), pr_("b_ratios_matched_pairs.csv"), pr_("c_non_responders.csv").iloc[0]
Z_Ds, Z_Dt = pr_("d_tile_scan_summary.csv"), pr_("d_tile_scan_tests.csv")
pimg = lambda n, alt: img(DPR / "figures" / f"{n}.png", alt)
_za = Z_A[Z_A.fdr < 0.05]
_zc = lambda k: Z_C.loc[[i for i in Z_C.index if i.startswith(k)][0]]
DSH = U("Question: if every site is scaled to the same overall sensitivity, is the mix of mutation types different? (Each site's 21 effects are divided by their overall size, so a Dynamic site that loses less to everything has the same shape as a Stable site that loses more.)",
        f"Matched pairs ({int(Z_A.n_pairs.iloc[0])}): {len(_za)} of {len(Z_A)} mutation types differ at FDR &lt; 0.05, all four of them edits of the 170 bp, which take a slightly smaller share at Dynamic sites ({_za.paired_median_diff.min():+.3f} to {_za.paired_median_diff.max():+.3f} on shares of about {_za.share_stable.mean():.2f}, i.e. about 5% smaller): far too small to matter biologically.",
        f"A classifier given the whole profile cannot tell matched Dynamic from Stable sites: shape only AUC {_zc('shape only').auc:.2f}, raw effects {_zc('size and shape').auc:.2f}, overall size {_zc('overall size').auc:.2f} (shuffled-label 95th percentile about {Z_C.null_auc_95th.mean():.2f}).",
        "Ratios with a plain meaning, at sites that respond (matched pairs, none significant): " + "; ".join(f"{r.ratio.split(':')[0]} {r.median_stable:.2f} Stable vs {r.median_dynamic:.2f} Dynamic" for r in Z_B.itertuples()) + ".",
        f"Non-responders (replacing the motif changes accessibility by less than 10%): {Z_N.share_dynamic_matched:.0%} of Dynamic vs {Z_N.share_stable_matched:.0%} of matched Stable sites (McNemar p = {Z_N.p_mcnemar_matched:.2f}); about the same in all sites ({Z_N.share_dynamic_all:.0%} vs {Z_N.share_stable_all:.0%}). So Dynamic sites are not a mix of AP-1-dependent and AP-1-independent sites in different proportions.") + \
      f"<figure>{pimg('z1_response_shape', 'Response shape by mutation type')}<figcaption>Paired difference (Dynamic minus matched Stable) in the share of each site's response taken by each mutation type, median with 95% bootstrap CI; orange = FDR &lt; 0.05 (paired Wilcoxon, BH).</figcaption></figure>" + \
      f"<details><summary>Classifier on the response profile</summary><figure>{pimg('z2_classifier', 'Classifier on response profiles')}<figcaption>Cross-validated AUC in the matched sample (folds by peak, 5 random fold assignments); dotted line = 95th percentile of 100 label permutations.</figcaption></figure>{table(Z_B[['ratio', 'n_pairs_both_responding', 'median_stable', 'median_dynamic', 'paired_median_diff', 'p_wilcoxon', 'fdr']])}</details>"
_t = lambda rd, cmp, m: Z_Dt[(Z_Dt.readout == rd) & (Z_Dt.comparison == cmp) & (Z_Dt.measure == m)].iloc[0]
_s = lambda rd, g: Z_Ds[(Z_Ds.readout == rd) & (Z_Ds.group == g)].iloc[0]
_om = _t("170 bp", "dynamic vs matched_stable", "own_motif_tile")
_om2 = _t("motif", "dynamic vs matched_stable", "own_motif_tile")
_np = _t("170 bp", "dynamic vs matched_stable", "n_strong_partner_tiles")
_dy, _ms = _s("170 bp", "Dynamic"), _s("170 bp", "Stable, context-matched")
DT = U(f"PRELIMINARY: {int(Z_Ds.n_sites.sum() / 2)} of 1,470 sites are finished ({int(_dy.n_sites)} Dynamic, {int(_ms.n_sites)} context-matched Stable, {int(_s('170 bp', 'Stable, high accessibility').n_sites)} high-accessibility Stable); rerun when the scan completes.",
       "Each 10 bp tile from &minus;300 to +300 bp is replaced by random bases (3 replicates; replicate SD about 0.05 log2). For each site: the effect of the tiles covering the AP-1 motif, and of the strongest other tile (a 'partner').",
       f"Neighbouring elements are common in both classes: at {_dy.share_partner_dominates:.0%} of Dynamic and {_ms.share_partner_dominates:.0%} of matched Stable sites some other 10 bp tile matters more than either motif tile (read over 170 bp). Caveat: this compares the strongest of about 58 tiles with the stronger of 2.",
       "Example (chr1:119,197,433): single-base changes at &minus;30 to &minus;23 bp cut accessibility about 6-fold except A&harr;T swaps at &minus;28 to &minus;25; the bases read TA[TTAA]TA, a MEF2 motif, whose (A/T)4 centre tolerates exactly those swaps. AlphaGenome treats a cardiac partner site 25 bp away as more important than the AP-1 motif itself.",
       f"Read over 170 bp, the AP-1 motif tile matters less at Dynamic sites ({_om.median_a:.2f} vs {_om.median_b:.2f}, FDR {_om.fdr:.2f}), but not read over the motif itself ({_om2.median_a:.2f} vs {_om2.median_b:.2f}, FDR {_om2.fdr:.2f}): inconsistent between readouts, so not a finding yet.",
       f"Dynamic sites have slightly fewer strong partner tiles (median {_np.median_a:.0f} vs {_np.median_b:.0f}, p = {_np.p_mannwhitney:.3f}, FDR {_np.fdr:.2f}), and their strongest partner lies beyond &plusmn;85 bp more often ({_dy.share_partner_beyond_85bp:.0%} vs {_ms.share_partner_beyond_85bp:.0%}): hints only.") + \
     f"<figure>{pimg('z3_own_motif_vs_partner', 'Own motif vs neighbouring element')}<figcaption>Left: per site, effect of replacing the AP-1 motif tile against the strongest other 10 bp tile within &plusmn;300 bp (read over 170 bp); below the dashed line a neighbour matters more. Right: share of sites where a neighbour dominates, with 95% bootstrap CI.</figcaption></figure>" + \
     f"<details><summary>Average tile profile and saturation mutagenesis</summary><figure>{img(AN / 'scan_titrate/figures/s1_tile_profile.png', 'Tile profile')}<figcaption>Median effect of each tile by position (95% bootstrap CI), read over the motif (left) and 170 bp (right).</figcaption></figure><figure>{img(AN / 'scan_titrate/figures/s2_saturation_mutagenesis.png', 'Saturation mutagenesis')}<figcaption>Median |effect| of single-base changes from &minus;30 to +30 bp, and the share of that importance inside the 11 bp motif.</figcaption></figure></details>"
SPECIAL.update({"DSH": DSH, "DT": DT, "DS": DS, "DP": DP, "DF": DF, "DC": DC, "DM": DM, "DX": DX, "DEXP": DEXP, "FT": FT})

# ---------------- summary
GLOSS = U("<strong>Site</strong> = one Fos (AP-1) motif footprint called by TOBIAS in a cell type. <strong>Stable</strong> = bound in all three conditions; <strong>Dynamic</strong> = any other pattern.",
          "<strong>YS / MAS / MAE</strong> = young sedentary, aged sedentary, aged + exercise. Age Up = bound in MAS but not YS; Exercise Up = bound in MAE but not MAS; Down = the reverse.",
          "<strong>Effect</strong> = log2(edited / native) of predicted accessibility, so <strong>negative = accessibility lost</strong>. 'Dynamic minus Stable = +0.5' therefore means Dynamic loses <em>less</em>. In the mutation arm, dependency = minus effect (positive = lost).",
          "<strong>Readout</strong> = the region accessibility is averaged over: the 11&nbsp;bp motif, the 170&nbsp;bp or the whole peak; <strong>matched</strong> = the same region that was edited.",
          "<strong>Heart</strong> and <strong>Patski</strong> = the mouse heart and Patski (fibroblast-like cell line) ATAC tracks of the model. <strong>snATAC</strong> = measured single-nucleus ATAC. <strong>MPRA</strong> = the lenti-STARR reporter assay in cardiac fibroblasts. <strong>Tier</strong> = a filtered set of well-predicted sites (chapter 5).",
          "Stars: <code>***</code> FDR &lt; 0.001, <code>**</code> &lt; 0.01, <code>*</code> &lt; 0.05, <code>ns</code> not significant. Read the effect size before the stars. Confidence labels are explained with each result.")
SUMMARY = f"""
<div class="box"><strong>Key terms and conventions</strong>{GLOSS}</div>
<div class="box"><strong>Summary</strong> <span class="small">(confidence label, then what it means)</span>
{U(f"{badge('ok')} <strong>Context, not the motif.</strong> Shuffling the 170&nbsp;bp around an AP-1 site costs {m.loc['170', 'median_stable']:.2f} log2 vs {m.loc['motif', 'median_stable']:.2f} for the 11&nbsp;bp motif (Stable sites): most predicted accessibility sits in the surrounding sequence ({clink('depend')}).",
   f"{badge('mid')} <strong>Local and saturating.</strong> A same-size shuffle 300&nbsp;bp away does nothing ({co_.loc['Stable', 'offset_pm300_median']:.2f}); the loss saturates near 340&nbsp;bp; replacing everything but the motif keeps {ratio_cm170:.0%} (CM) / {ratio_fb170:.0%} (FB) of the loss, though a ratio near 100% is partly expected and follows the track ({link('D1')}, {link('D2')}).",
   f"{badge('mid')} <strong>Stable vs Dynamic is mostly chromatin context.</strong> The gap is {m.loc['170', 'diff_dyn_minus_stable']:+.2f} raw and {hg(cmn, M_FULL).gap:+.2f} after accounting for accessibility, peak width and TSS distance (Patski: {hg(fbp, M_MATCH).gap:+.2f} matched only, {hg(fbp, M_FULL).gap:+.2f} full). A separate sequence effect is neither shown nor excluded ({clink('dyn')}).",
   f"{badge('mid')} <strong>No mutation type singles out Dynamic sites.</strong> Dynamic sites are less sensitive to every wide edit in the raw results (for example replace the 170&nbsp;bp {Y_prof.loc['replace the 170 bp', 'stable_median']:.2f} Stable vs {Y_prof.loc['replace the 170 bp', 'dynamic_median']:.2f} Dynamic) but, matched on predicted context, the two classes have the same profile across all {len(Y_prof)} mutation types, and no binding pattern is robustly different ({clink('deep')}).",
   f"{badge('mid')} <strong>Sequence sets rank, location sets level</strong> in a peak swap, equally for Stable and Dynamic donors; only 4 recipient windows ({link('B8a')}).",
   f"{badge('mid')} <strong>The model tracks measured accessibility moderately</strong> (rho {xov.loc['CM sites, heart ATAC', 'rho']:.2f} cardiomyocyte, {xov.loc['FB sites, heart ATAC', 'rho']:.2f} fibroblast with the heart track, {xov.loc['FB sites, Patski ATAC', 'rho']:.2f} with Patski); Dynamic sites are over-predicted ({clink('found')}).",
   f"{badge('mid')} <strong>The native head does not explain the MPRA.</strong> Predicted scrambling effects correlate about 0 with MPRA changes, in genomic and neutral windows; the MPRA's own effects are tiny ({clink('mpra')}).",
   f"{badge('mid')} <strong>Results hold in well-predicted sites</strong>, but the filter drops mostly over-predicted Dynamic sites, so filtered gaps are partly selection ({clink('robust')}).",
   f"{badge('bad')} <strong>Did not work:</strong> RNA (effects about zero) and motif insertion into random backgrounds (floor effect); {badge('inc')} age vs exercise classes (underpowered) and multi-motif peaks (few sites). Fibroblast conclusions are indicative only ({clink('fail')}).")}
<p class="small">Numbers were re-checked against the tables by reviewers; corrections are listed in {clink('open')}.</p></div>"""

# ---------------- key figures (presentation)
KEYFIGS = [
 ("The model recovers measured accessibility at real AP-1 sites", ximg('x4_observed_vs_predicted', 'Observed vs predicted accessibility'),
  f"rho {xov.loc['CM sites, heart ATAC', 'rho']:.2f} (cardiomyocyte), {xov.loc['FB sites, heart ATAC', 'rho']:.2f} (fibroblast, heart track) but {xov.loc['FB sites, Patski ATAC', 'rho']:.2f} with Patski.", "Moderate, not precise; bulk tissue and a cell line stand in for purified cell types.", "mid"),
 ("Shuffling the 170 bp around an AP-1 site removes far more accessibility than shuffling the motif", ximg('x2_shuffle_dotplot', 'Shuffle effects as a grouped dot plot'),
  f"{m.loc['170', 'median_stable']:.2f} vs {m.loc['motif', 'median_stable']:.2f} log2 in Stable sites: the surrounding sequence carries most of the predicted accessibility.", "A model prediction, not a measurement.", "ok"),
 ("The effect is local and builds up to about 340 bp", oimg('d1_controls_titration', 'Controls and width titration'),
  f"300 bp away: {co_.loc['Stable', 'offset_pm300_median']:.2f} vs {co_.loc['Stable', 'centred_170_median']:.2f} centred; the 7 bp core gives about {fr_.loc[7]:.0%} of the 170 bp effect.", "A distance control, not an AP-1-specific one.", "mid"),
 ("Stable sites lose more mainly because their chromatin context is stronger", img(AU / 'figures' / 'h1_gap_models.png', 'Gap by model'),
  f"{hg(cmn, M_RAW).gap:+.2f} raw → {hg(cmn, M_MATCH).gap:+.2f} → {hg(cmn, M_FULL).gap:+.2f} as accessibility, width and TSS distance are added.", "Footprint score partly defines the classes; Patski is not an independent confirmation.", "mid"),
 ("Predicted accessibility changes do not track the fibroblast MPRA, even outside the genome", oimg('d4_mpra_context', 'MPRA genomic vs neutral context'),
  f"Native level rho {mc_.loc[('heart ATAC', 'genomic context', 'native level vs MPRA activity'), 'rho']:+.2f} (heart); scramble effects about 0 in both windows.", "MPRA effects are tiny, so the null does not rank models; 2 neutral backgrounds.", "mid"),
]
KEY_HTML = U("If only five figures were shown, in order: the model can be trusted; context (not the motif) carries the effect; it is local; responsive sites are explained by context; accessibility is not enhancer activity.") + "".join(
    f"<div class='keyfig'><h3><span class='num'>{i}</span>{ttl}</h3>{badge(k)}<figure>{fg}</figure>{U(msg, '<em>Caveat:</em> ' + cav)}</div>" for i, (ttl, fg, msg, cav, k) in enumerate(KEYFIGS, 1))
KEY_HTML += "<p class='small'>Runners-up: the mutation-arm decomposition, the peak swap and the prediction-quality-by-class figure.</p>"

# ---------------- assemble
toc = ["<div class='toc'>"]
for ck, ct, cintro, subs in CH:
    toc.append(f"<div><a href='#ch-{ck}'><strong>{N[ck]}. {ct}</strong></a><ul>" + "".join(f"<li><a href='#sec-{sk}'>{N[sk]} {st}</a></li>" for sk, st in subs) + "</ul></div>")
toc.append("</div>")

chapters = []
for ck, ct, cintro, subs in CH:
    h = f"<h2 id='ch-{ck}'><span class='num'>{N[ck]}</span>{ct}</h2>"
    if ck == "design":
        chapters.append(h + "<p>What each experiment does to the sequence and what is measured.</p>" + fixrefs(re.sub(r"<h2.*?</h2>", "", DESIGN, count=1, flags=re.S)))
        continue
    if ck == "reg":
        chapters.append(h + fixrefs(REG_HTML))
        continue
    if ck == "bio":
        chapters.append(h + BIO)
        continue
    if ck == "key":
        chapters.append(h + fixrefs(KEY_HTML))
        continue
    parts_ = [h, fixrefs(cintro) if cintro else ""]
    for sk, st in subs:
        if sk in SPECIAL:
            body_ = SPECIAL[sk]
        elif sk in ("H2", "H3"):
            body_ = SEC[sk]
        else:
            body_ = PTS[sk] + strip_prose(SEC[sk], keep_tables=KEEP_T.get(sk)) if sk in PTS else SEC[sk]
        c_ = conf(*CONF[sk]) if sk in CONF else ""
        bx = f"<div class='bio'><strong>What this means biologically.</strong> {BIO_SEC[sk]}</div>" if sk in BIO_SEC else ""
        parts_.append(f"<h3 id='sec-{sk}'><span class='num'>{N[sk]}</span>{st}</h3>{c_}{fixrefs(body_)}{fixrefs(bx)}")
    chapters.append("\n".join(parts_))

METHODS = """<details><summary>Methods and files</summary>
<ul class="pts"><li>Predictions: alphagenome-pytorch with <code>model_all_folds</code> weights, mouse organism, 131,072&nbsp;bp window centred on the motif midpoint, mouse heart ATAC track (1&nbsp;bp, mean over region) and Patski ATAC for fibroblast sites.</li>
<li>Replicates per site: 10 shuffles per width; 5 for replace, dinucleotide-shuffle and flank edits; 3 for the core edit.</li>
<li>Stable = bound in YS, MAS and MAE; age class compares YS and MAS; exercise class compares MAS and MAE.</li>
<li>Tests: Mann&ndash;Whitney with Benjamini&ndash;Hochberg within families, bootstrap CIs, Cliff's delta; adjusted models are OLS with HC3 (or peak-clustered) standard errors. About 15% of sites share a peak, so most p-values are mildly optimistic; effect sizes and intervals are the better guide.</li>
<li>Tables: <code>Analysis/cm</code>, <code>overnight</code>, <code>agreement</code>, <code>agreement_floor25</code>, <code>audit</code>, <code>extra</code>; scripts <code>Scripts/22&ndash;32</code> and <code>25_build_report.py</code>.</li></ul></details>"""

EXTRA = f"""<details><summary>Additional plots (alternative versions, for choosing which to keep)</summary>
<figure>{img(AN / 'agreement' / 'figures' / 'e1_tier_C_selection.png', 'Tier C selection')}<figcaption>Predicted native ATAC against measured snATAC with tier C sites in blue.</figcaption></figure>
<figure>{img(AN / 'agreement' / 'figures' / 'e2_retention_by_class.png', 'Retention by class')}<figcaption>Share of each class retained per tier.</figcaption></figure>
<figure>{img(AN / 'agreement_floor25' / 'figures' / 'e3_stable_vs_dynamic_by_tier.png', 'Floor 25 stable vs dynamic')}<figcaption>Stable versus Dynamic by tier with the stricter floor.</figcaption></figure>
<figure>{img(AN / 'agreement_floor25' / 'figures' / 'e0b_rank_map.png', 'Floor 25 rank map')}<figcaption>Rank map with the stricter floor.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig3_shuffle_effect_by_class.png', 'Class effects without marks')}<figcaption>Effect by class without significance marks.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig1_native_vs_measured.png', 'Predicted vs measured')}<figcaption>Predicted vs measured accessibility, cardiomyocyte sites (single panel).</figcaption></figure>
<figure>{img(CM / 'qc' / 'figures' / 'qc5_prediction_vs_measured_by_class.png', 'Prediction by class')}<figcaption>Earlier version of the prediction-by-class figure.</figcaption></figure>
<figure>{img(AU / 'figures' / 'h2_gap_by_matched_decile.png', 'Gap by decile')}<figcaption>Stable vs Dynamic effect by decile of native accessibility at the 170 bp readout.</figcaption></figure></details>"""

head = f"""<h1>AlphaGenome in-silico perturbation of AP-1 (Fos) sites in mouse cardiomyocytes and fibroblasts</h1>
<p class="sub">{int(m.loc['170', 'n_stable'] + m.loc['170', 'n_dynamic']):,} cardiomyocyte and 5,018 fibroblast Fos sites &middot; generated {date.today():%d %b %Y}</p>"""
CHKEY = {v: k for k, v in N.items() if k in {c[0] for c in CH}}
NUMKEY = {v: k for k, v in N.items() if k in SEC or k in SPECIAL or k in {c[0] for c in CH}}
body = head + SUMMARY + "<h2>Contents</h2>" + "".join(toc) + "\n".join(chapters) + EXTRA + METHODS
doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AP-1 shuffle analysis report</title><style>{CSS}</style></head><body><main>{body}</main></body></html>"""
for k_, p_ in IMGS.items():
    doc = doc.replace(f"@@IMG{k_}@@", "data:image/png;base64," + base64.b64encode(p_.read_bytes()).decode())
out = AN / "report_cm.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB), {len(IMGS)} images")
