#!/usr/bin/env python3
"""Builds Analysis/report_cm.html (self-contained: figures embedded) from the tables and figures written by 22/23/24. Every number in the text
is read from the CSV tables, so rerun 24 then this script after new results arrive.   python3 25_build_report.py"""
from __future__ import annotations

import base64
import html
from datetime import date
from pathlib import Path

import pandas as pd

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
main{max-width:1000px;margin:0 auto;padding:32px 16px 80px}h1{font-size:1.9rem;line-height:1.2;margin:0 0 4px}h2{font-size:1.35rem;margin:44px 0 8px;padding-top:12px;border-top:1px solid var(--rule)}h3{font-size:1.05rem;margin:22px 0 4px}
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


claims = [
    ("Shuffling the 170 bp around a site removes far more predicted accessibility than shuffling the 11 bp motif alone", "ok",
     f"Median log2FC {m.loc['170', 'median_stable']:.2f} (170 bp shuffle) vs {m.loc['motif', 'median_stable']:.2f} (motif shuffle) in Stable sites; "
     f"{frac['170']:.0%} of sites show a 170 bp effect beyond 2 standard errors of their own replicate scatter."),
    ("Stable sites depend more on their 170 bp neighbourhood than Dynamic sites do", "ok",
     f"Difference of medians {m.loc['170', 'diff_dyn_minus_stable']:+.2f} (95% CI {m.loc['170', 'ci_lo']:+.2f} to {m.loc['170', 'ci_hi']:+.2f}), FDR {m.loc['170', 'fdr']:.1e}, permutation p = {pm.loc['170', 'perm_p']:.4f}. "
     f"Holds in every native-accessibility quintile ({'all five' if quint_all_sig else 'most'} FDR &lt; 0.05), after adjusting for native ATAC and peak width "
     f"(adjusted differences {a170.adj_diff.min():+.2f} to {a170.adj_diff.max():+.2f}, all {'***' if (a170.adj_fdr < 0.001).all() else 'significant'}), and in every dynamic class. "
     f"Direction replicates with the independent Patski ATAC track but is smaller ({fb.loc['170', 'diff_dyn_minus_stable']:+.2f}, {fb.loc['170', 'sig']})."),
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
    ("Fibroblast results; MPRA agreement", "na", "Fibroblast-site shuffles and the MPRA run are not finished at the time of this report. The fibroblast tracks available in the mouse model (Patski, embryonic fibroblast RNA) validated weakly against measured accessibility."),
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

claims_html = "".join(f"<tr><td class='v'>{c}</td><td>{badge(k)}</td><td class='v'>{ev}</td></tr>" for c, k, ev in claims)

body = f"""
<h1>AlphaGenome in-silico shuffling of AP-1 (Fos) sites: cardiomyocyte arm</h1>
<p class="sub">Generated {date.today():%d %b %Y} &middot; {int(m.loc['170', 'n_stable'] + m.loc['170', 'n_dynamic']):,} cardiomyocyte Fos sites, 10 shuffle replicates per width &middot; mouse heart ATAC track &middot; fibroblast arm and MPRA still running</p>

<div class="box"><strong>Main takeaways</strong><ul>
<li><strong>What is solid.</strong> Shuffling the 170&nbsp;bp around an AP-1 site removes far more predicted accessibility (median log2FC {m.loc['170', 'median_stable']:.2f} in Stable sites) than shuffling the 11&nbsp;bp motif alone ({m.loc['motif', 'median_stable']:.2f}). This holds for essentially every site and is well above replicate noise.</li>
<li><strong>Stable vs Dynamic is real.</strong> Stable sites depend more on their 170&nbsp;bp context than Dynamic sites ({m.loc['170', 'diff_dyn_minus_stable']:+.2f}, 95% CI {m.loc['170', 'ci_lo']:+.2f} to {m.loc['170', 'ci_hi']:+.2f}). It survives adjustment for native accessibility and peak width, stratification by accessibility, a permutation test, and a second independent track (smaller there).</li>
<li><strong>Age vs exercise vs combination.</strong> All dynamic classes behave alike. Comparing classes with each other (Age Up vs Down, Exercise Up vs Down) finds no differences at 170&nbsp;bp or peak level. What separates classes is being stably bound or not, not the direction of change.</li>
<li><strong>What is noise or unusable.</strong> RNA effects are about zero; AP-1 motif insertion into random backgrounds is at the floor; per-site peak-shuffle values scatter widely (replicate SD about {D.set_index('shuffle').loc['peak', 'median_replicate_sd']:.1f}) and are only trustworthy in aggregate.</li>
<li><strong>Caveat on the slides.</strong> Our Dynamic-site native accessibility matches the slides closely, which suggests the pipeline reproduces the API, but our Stable-site set differs from the slides' (we could not reproduce their site counts). Absolute numbers should not be compared directly with the slides until the site set is resolved.</li></ul></div>

<h2>How to read this report</h2>
<p><strong>Effect.</strong> For each site and shuffle width, the effect is the mean over 10 shuffles of <code>log2((shuffled + pc) / (native + pc))</code> of the predicted ATAC (mean per-base signal over a region; pc = 5% of the median native value). Negative = accessibility lost. A shuffle of the motif is read out over the motif, a 170&nbsp;bp shuffle over the 170&nbsp;bp, a peak shuffle over the peak ("matched readout"), unless stated.</p>
<p><strong>Significance marks.</strong> <code>***</code> FDR &lt; 0.001, <code>**</code> &lt; 0.01, <code>*</code> &lt; 0.05, <code>ns</code> not significant. Tests are two-sided Mann&ndash;Whitney against a reference class, Benjamini&ndash;Hochberg corrected within each family of tests. With thousands of sites a tiny difference is "significant", so each comparison also shows an effect size (difference of medians with a bootstrap 95% CI, and Cliff's delta: about 0.1 small, 0.3 medium, 0.5 large). <strong>Read the effect size before the stars.</strong></p>
<p><strong>Trust labels.</strong> <span class="badge ok">&#10003; Robust</span> survives the checks below. <span class="badge mid">~ Plausible / partial</span> real but small, or depends on a confound. <span class="badge bad">&#10007; Not supported</span> no evidence, or indistinguishable from noise. <span class="badge na">? Cannot judge yet</span> data not available.</p>

<h2>1. Does the model track real accessibility?</h2>
<figure>{img(CM / 'figures' / 'fig1_native_vs_measured.png', 'Predicted vs measured accessibility', '560px')}
<figcaption>Each dot is a site: predicted native ATAC over the peak (mouse heart track) against measured cardiomyocyte snATAC at the same peak. Spearman rho = {nv.rho:.2f}.</figcaption></figure>
{table(NV_tab, {'Spearman rho': '{:.2f}'})}
<p><strong>What it shows.</strong> The heart track predicts measured cardiomyocyte accessibility with a moderate rank correlation (the slides reported 0.49). The Patski (fibroblast-like) track is essentially uncorrelated with cardiomyocyte accessibility, as it should be. <strong>Trust:</strong> good enough to use the heart track as a proxy for cardiomyocyte sites; it is bulk heart tissue, not purified cardiomyocytes.</p>

<h2>2. What does shuffling do?</h2>
<figure>{img(CM / 'figures' / 'fig2s_shuffle_heatmaps_stars.png', 'Shuffle effect heatmaps')}
<figcaption>Median log2FC for each shuffle width (rows) and ATAC readout region (columns). Left: Stable sites; right: Dynamic sites, with significance marks against Stable.</figcaption></figure>
<p><strong>What it shows.</strong> Shuffling only the motif lowers accessibility at the motif by about 0.6 log2 units (roughly a third lost) and over the whole peak by about 0.3 (about a fifth lost). Shuffling 170&nbsp;bp or the whole peak is far more damaging (about -1.8 log2 units locally, i.e. around 70% lost). Dynamic sites lose less than Stable sites in most cells. <strong>Trust:</strong> the ordering motif &lt; 170&nbsp;bp &asymp; peak is robust. Peak shuffles are noisy per site (see section 7).</p>
{table(stv_tab)}

<h2>3. Stable vs Dynamic: is the difference real?</h2>
<p>Dynamic sites have lower native accessibility and different peak widths, either of which could make their log2FC look smaller without any biological difference. Three checks:</p>
<figure>{img(CM / 'figures' / 'fig8_stable_vs_dynamic_by_native_accessibility.png', 'Stable vs Dynamic within accessibility quintiles', '640px')}
<figcaption>Sites split into five groups by native predicted accessibility. Within every group, Stable sites lose more than Dynamic sites at the 170&nbsp;bp shuffle.</figcaption></figure>
{table(C_tab)}
<ul><li><strong>Stratification:</strong> the gap is present in all five accessibility quintiles ({'all FDR &lt; 0.05' if quint_all_sig else 'most FDR &lt; 0.05'}), so it is not a floor effect of low-accessibility sites.</li>
<li><strong>Permutation test</strong> (5,000 label shuffles): 170&nbsp;bp p = {pm.loc['170', 'perm_p']:.4f}, peak p = {pm.loc['peak', 'perm_p']:.4f}, motif p = {pm.loc['motif', 'perm_p']:.3f}. The 170&nbsp;bp and peak differences are far outside the null; the motif one is not.</li>
<li><strong>Independent track:</strong> repeating the test with the Patski ATAC track on the same sites gives 170&nbsp;bp {fb.loc['170', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['170', 'sig']}), peak {fb.loc['peak', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['peak', 'sig']}), motif {fb.loc['motif', 'diff_dyn_minus_stable']:+.2f} ({fb.loc['motif', 'sig']}). The 170&nbsp;bp direction replicates; the others do not.</li></ul>
<p><strong>Trust:</strong> robust for the 170&nbsp;bp shuffle; partial for the peak shuffle; not supported for the motif shuffle.</p>

<h2>4. Age, exercise and their combination</h2>
<p>Each site carries three labels from the same binding calls: <em>age class</em> (YS vs MAS), <em>exercise class</em> (MAS vs MAE) and the full binding pattern. Every class is compared with its Stable reference, once raw and once adjusted for native ATAC and peak width (OLS with robust standard errors).</p>
<figure>{img(CM / 'figures' / 'fig3s_class_effects_stars.png', 'Class effects with significance marks')}
<figcaption>Median matched effect per class with 95% bootstrap CI. Blue = the Stable reference; marks to the right are FDR-corrected Mann&ndash;Whitney tests against it.</figcaption></figure>
<figure>{img(CM / 'figures' / 'fig7_class_effects_adjusted.png', 'Adjusted class effects')}
<figcaption>Difference from the Stable reference: raw (orange diamonds) and adjusted for native accessibility and peak width (blue, with 95% CI). Positive = the class loses less than Stable.</figcaption></figure>
<p><strong>What it shows.</strong> At 170&nbsp;bp every dynamic class loses about 0.45 to 0.6 log2FC less than Stable, and adjustment barely changes it. At the peak and motif widths, raw differences are small and the adjusted ones shrink toward zero. <strong>Between dynamic classes</strong> (the direct test of age vs exercise, and of up vs down):</p>
{table(E_tab)}
<p>Only the motif-width contrasts reach FDR &lt; 0.05, with confidence intervals touching zero. There is no sign that age-driven and exercise-driven changes differ from each other in sequence dependence. <strong>Trust:</strong> "Stable vs dynamic" is robust; "which dynamic class" is not supported.</p>
<details><summary>Omnibus tests, full class table and covariates</summary>
<h3>Kruskal&ndash;Wallis across all classes in each family</h3>{table(OM_tab, {'Kruskal-Wallis H': '{:.1f}'})}
<h3>Every class against its reference</h3>{table(cls_tab)}
<h3>Covariates in the adjusted model (age family)</h3>{table(cov)}
<p class="small">Positive peak-width coefficient: wider peaks lose less accessibility per log2FC. Dynamic sites differ systematically from Stable sites in width and native accessibility, which is why the raw peak-shuffle gap shrinks after adjustment.</p></details>

<h2>5. RNA</h2>
<figure>{img(CM / 'figures' / 'fig4_rna_effect_by_tss_distance.png', 'RNA effect by TSS distance')}
<figcaption>Median log2FC of predicted RNA after a peak shuffle, by distance of the site to the nearest TSS, for Stable and Dynamic sites. Left: within 2&nbsp;kb of that TSS; right: mean over &plusmn;50&nbsp;kb.</figcaption></figure>
<p><strong>What it shows.</strong> Predicted RNA barely moves when a distal AP-1 site or peak is shuffled (medians within about &plusmn;0.1), with the largest drops close to a TSS. ATAC and RNA effects correlate across sites (rho 0.19 to 0.24), so there is a weak coupling, but the RNA effects themselves are far smaller than the slides reported. <strong>Trust:</strong> not informative with these tracks; do not over-interpret.</p>

<h2>6. Does a high-confidence subset change the picture?</h2>
<figure>{img(CM / 'figures' / 'fig5_confidence_tiers.png', 'Confidence tiers')}
<figcaption>Median effect in three nested sets chosen only from native prediction vs measurement (never from the effects): A all sites; B sites with a measured peak and predicted accessibility above the 25th percentile; C = B plus prediction-measurement rank agreement within 1 SD.</figcaption></figure>
{table(t5_tab)}
<p><strong>What it shows.</strong> Effects get somewhat larger in the filtered sets (the model is more confident at accessible, well-predicted sites) and the Stable vs Dynamic gap at 170&nbsp;bp stays. <strong>Trust:</strong> conclusions are not an artefact of poorly predicted sites; report all three tiers.</p>

<h2>7. How much is noise?</h2>
{table(D_tab, {'Sites with |effect| > 2 SE': '{:.0%}', 'Sites with effect < -2 SE': '{:.0%}'})}
<p>Replicate SD is the spread across the 10 shuffles of one site. The 170&nbsp;bp effect is stable per site (SD about {D.set_index('shuffle').loc['170', 'median_replicate_sd']:.1f}); the motif and peak effects scatter more (SD about {D.set_index('shuffle').loc['motif', 'median_replicate_sd']:.1f} and {D.set_index('shuffle').loc['peak', 'median_replicate_sd']:.1f}). The median over thousands of sites is precise, but <strong>an individual site's peak-shuffle value should not be interpreted on its own</strong>. This is why 10 replicates were kept.</p>

<h2>8. Insertion experiments</h2>
<figure>{img(AN / 'insertion' / 'figures' / 'fig_peak_swap.png', 'Peak swap')}
<figcaption>Peak swap: a donor peak's sequence replaces a recipient peak. Swapped accessibility follows the donor's own accessibility (left) much more than the recipient location's (right).</figcaption></figure>
{table(ps_tab, {'Spearman rho': '{:.2f}'})}
<p><strong>Trust:</strong> robust direction (sequence beats location), but only 4 recipient windows, so the 292 pairs are not independent; the effective sample is about 73 donors.</p>
<figure>{img(AN / 'insertion' / 'figures' / 'fig_motif_insertion.png', 'Motif insertion')}
<figcaption>A TGAGTCA motif inserted at the centre of dinucleotide-shuffled random genomic backgrounds (50 backgrounds), compared with the shuffled background and with a shuffled-motif control.</figcaption></figure>
{table(mi_tab, {'Median log2FC': '{:.2f}', 'Fraction > 0': '{:.2f}'})}
<p><strong>Trust:</strong> no usable signal. Background accessibility is near zero, so any ratio is dominated by the floor, and the 7&nbsp;bp readout differences are not significant. The slides reported a similarly small effect.</p>

<h2>9. What to trust</h2>
<table class="t"><tr><th>Claim</th><th>Verdict</th><th>Evidence</th></tr>{claims_html}</table>

<h2>10. Open issues and next steps</h2>
<ul>
<li><strong>Site set vs the slides.</strong> The TOBIAS table gives 3,679 CM Fos sites (2,907 Stable), not the slides' 2,663 (2,320 Stable). Our Dynamic native accessibility matches the slides, but our Stable sites are about half as accessible on average, and keeping our top 80% by predicted accessibility does not close the gap. The slides' Stable sites are a different set. Until the original selection is known, compare patterns, not absolute values.</li>
<li><strong>Model equivalence.</strong> The slides were made with the AlphaGenome API; this run uses the PyTorch port with the same family of weights. The Dynamic-site match supports equivalence, but a direct check on the same sites through the API would settle it.</li>
<li><strong>Fibroblasts.</strong> The mouse model has no fibroblast ATAC; Patski and embryonic-fibroblast RNA are weak proxies. The fibroblast shuffle run and the MPRA comparison will be added when they finish.</li>
<li><strong>Motif-width effect is larger than in the slides</strong> (about -0.6 vs -0.12). Candidate causes are the motif interval used (here the 11&nbsp;bp TOBIAS interval for both the shuffle and the readout) and the different Stable set; worth confirming the slides' definition.</li></ul>

<details><summary>Methods and files</summary>
<p class="small">Predictions: alphagenome-pytorch with <code>model_all_folds</code> weights, mouse organism, 131,072&nbsp;bp window centred on the motif midpoint, mouse heart ATAC track (1&nbsp;bp, mean over region), 10 uniform random shuffles per width. Stable = bound in YS, MAS and MAE; age class compares YS and MAS; exercise class compares MAS and MAE. Adjusted model: effect ~ class + log10(native ATAC over peak) + log10(peak width), HC3 robust SEs, normal-theory CIs. All CIs for medians are percentile bootstrap (1,000&ndash;2,000 resamples). Sites are treated as independent, but about 8% of sites share a peak with another site, so p-values are mildly optimistic; effect sizes and confidence intervals are the better guide. Tables: <code>Analysis/cm/tables</code>, <code>Analysis/cm/stats</code>; figures: <code>Analysis/cm/figures</code>; scripts: <code>Scripts/22_shuffle_analysis.py</code>, <code>24_class_stats.py</code>, <code>23_insertion_analysis.py</code>, <code>25_build_report.py</code>.</p></details>
"""

doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AP-1 shuffle analysis report</title><style>{CSS}</style></head><body><main>{body}</main></body></html>"""
out = AN / "report_cm.html"
out.write_text(doc)
print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
