#!/usr/bin/env python3
"""Data-generation QC and additional analyses for the cardiomyocyte shuffle arm, plus the MPRA comparison.
Reads Katana/results/shuffle_cm_chunk*.csv (+ mpra_chunk*.csv), site tables, measured snATAC (total and per condition), job logs.
Writes tables to Analysis/cm/qc/tables and figures to Analysis/cm/qc/figures; used by 25_build_report.py.   python3 26_qc_and_extra.py

QC:  qc1 site composition / chromosomes      qc2 site geometry (peak width, TSS distance, motifs per peak)
     qc3 native predicted ATAC distributions (vs the slides' QC values)   qc4 class labels vs independent measured snATAC by condition
     qc5 prediction-vs-measurement by class   qc6 shuffled-vs-native QC (as on the slides)   qc7 replicate convergence   qc8 GPU / chunk consistency
Analyses: an1 effect distributions  an2 effect vs native accessibility / peak width / TSS distance  an3 motif vs context dependence (site types)
          an4 TOBIAS footprint score vs predicted dependence   an5 several Fos motifs in one peak   an6 promoter-proximal vs distal   an7 MPRA
"""
from __future__ import annotations

import glob
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu, spearmanr, t as tdist

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "Analysis" / "cm" / "qc"
(QC / "tables").mkdir(parents=True, exist_ok=True)
(QC / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(11)

INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
WIDTHS = ["motif", "170", "peak"]
READ = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
LABEL = {"motif": "Motif", "170": "170 bp", "peak": "Peak"}
SLIDE_NATIVE = {"Stable": {"motif": 1.760, "170": 1.747, "peak": 0.833}, "Dynamic": {"motif": 0.386, "170": 0.422, "peak": 0.337}}  # QC slide, native


def save(fig, name):
    fig.savefig(QC / "figures" / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def rho_ci(a, b):
    m = pd.notna(a) & pd.notna(b)
    n = int(m.sum())
    if n < 8:
        return np.nan, np.nan, np.nan, np.nan, n
    r, p = spearmanr(np.asarray(a)[m], np.asarray(b)[m])
    se = 1.06 / np.sqrt(n - 3)
    return r, np.tanh(np.arctanh(r) - 1.96 * se), np.tanh(np.arctanh(r) + 1.96 * se), p, n


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_ci(x, n=1000):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return tuple(np.percentile(np.median(rng.choice(x, (n, len(x))), 1), [2.5, 97.5])) if len(x) >= 5 else (np.nan, np.nan)


def binned(x, y, k=8):
    df = pd.DataFrame({"x": x, "y": y}).dropna()
    df["b"] = pd.qcut(df.x, k, duplicates="drop")
    g = df.groupby("b", observed=True)
    return g.x.median(), g.y.median(), g.y.quantile(0.25), g.y.quantile(0.75), g.size()


def ols_hc3(X, y):
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    beta = XtXi @ X.T @ y
    e = y - X @ beta
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    cov = XtXi @ (X.T @ (X * (e / (1 - h))[:, None] ** 2)) @ XtXi
    se = np.sqrt(np.diag(cov))
    return beta, se, 2 * tdist.sf(np.abs(beta / se), n - k)


# ---------------- data
files = sorted(glob.glob(str(ROOT / "Katana" / "results" / "shuffle_cm_chunk*.csv")))
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
site_tab = pd.read_csv(ROOT / "Katana" / "data" / "cm_fos_sites.csv")
site_tab["chunk"] = np.arange(len(site_tab)) % 40
sites = site_tab.set_index("site_id")
cond = pd.read_csv(ROOT / "Katana" / "data" / "measured_by_condition.csv")
cond = cond[cond.site_table == "cm"].set_index("site_id")
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id")
lab = d.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic", "age_class", "ex_class", "binding_category", "tss_dist"]]
lab = lab.join(sites[["peak_width", "peak_coord", "seqnames", "YS_score", "MAS_score", "MAE_score", "in_YS", "in_MAS", "in_MAE", "chunk"]])
lab["age_class"] = lab.age_class.replace({"Unbound": "Neither YS nor MAS (= MAE only)"})
lab["ex_class"] = lab.ex_class.replace({"Unbound": "Neither MAS nor MAE (= YS only)"})
SD = lab.stable_dynamic
nat = d[(d.condition == "native") & (d.trackset == "CM")].set_index("site_id")
lab["native_peak"], lab["native_motif"], lab["native_170"] = nat.atac_peak, nat.atac_motif, nat.atac_170
lab["log_nat"] = np.log10(lab.native_peak + 1e-3)
lab["log_pw"] = np.log10(lab.peak_width)
lab["abs_tss"] = lab.tss_dist.abs()
lab["fp_score"] = np.where(lab.in_YS, lab.YS_score, -1)
for c_, s_ in (("in_MAS", "MAS_score"), ("in_MAE", "MAE_score")):
    lab["fp_score"] = np.maximum(lab.fp_score, np.where(lab[c_], lab[s_], -1))


def rep_matrix(w):
    """site x replicate matrix of log2FC for the matched readout of shuffle width w (CM tracks)."""
    col = READ[w]
    pc = 0.05 * nat[col].median()
    e = d[(d.condition == f"shuffle@{w}") & (d.trackset == "CM")][["site_id", "rep", col]].copy()
    e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat[col]) + pc))
    return e.pivot(index="site_id", columns="rep", values="l")


REP = {w: rep_matrix(w) for w in WIDTHS}
EFF = {w: REP[w].mean(1) for w in WIDTHS}
Ef = pd.DataFrame(EFF)
lab = lab.join(Ef.rename(columns=lambda w: f"eff_{w}"))
print(f"{len(lab)} sites loaded from {len(files)} chunks")

# ================= QC1 composition
order = ["YS+MAS+MAE", "YS+MAS", "YS+MAE", "MAS+MAE", "YS only", "MAS only", "MAE only"]
cnt = lab.binding_category.value_counts().reindex(order)
cross = pd.crosstab(lab.age_class, lab.ex_class)
cross.to_csv(QC / "tables" / "qc1_age_by_exercise_counts.csv")
chrom = pd.crosstab(lab.seqnames, lab.stable_dynamic)
chrom = chrom.loc[sorted(chrom.index, key=lambda c: (len(c), c))]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6), gridspec_kw={"width_ratios": [1.1, 1.0, 1.7]})
ax = axes[0]
ax.barh(range(7), cnt.values[::-1], color=[BLUE if o == "YS+MAS+MAE" else ORANGE for o in order[::-1]], height=0.65)
ax.set_yticks(range(7), order[::-1])
for i, v in enumerate(cnt.values[::-1]):
    ax.text(v + 30, i, f"{v:,}", va="center", fontsize=8.5)
ax.set_xlabel("Sites")
ax.set_title("Binding pattern (YS / MAS / MAE)", fontsize=10, loc="left")
ax = axes[1]
ax.grid(False)
cr = cross.loc[["Stable", "Age Up", "Age Down", "Neither YS nor MAS (= MAE only)"], ["Stable", "Exercise Up", "Exercise Down", "Neither MAS nor MAE (= YS only)"]]
ax.imshow(np.log10(cr.values + 1), cmap="Blues", aspect="auto")
ax.set_xticks(range(4), ["Stable", "Ex Up", "Ex Down", "YS only"], rotation=0)
ax.set_yticks(range(4), ["Stable", "Age Up", "Age Down", "MAE only"])
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{cr.values[i, j]:,}" if cr.values[i, j] else "", ha="center", va="center", fontsize=8.5, color=INK if cr.values[i, j] < 600 else "white")
ax.set_xlabel("Exercise class (MAS vs MAE)")
ax.set_ylabel("Age class (YS vs MAS)", labelpad=22)
ax.set_title("Age x exercise labels", fontsize=10, loc="left")
ax = axes[2]
x = np.arange(len(chrom))
ax.bar(x, chrom.Stable, color=BLUE, label="Stable", width=0.75)
ax.bar(x, chrom.Dynamic, bottom=chrom.Stable, color=ORANGE, label="Dynamic", width=0.75)
ax.set_xticks(x, [c.replace("chr", "") for c in chrom.index], fontsize=8)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Sites")
ax.legend(frameon=False)
ax.set_title("Sites per chromosome", fontsize=10, loc="left")
save(fig, "qc1_site_composition")

# ================= QC2 geometry
lab["motifs_in_peak"] = lab.groupby("peak_coord").peak_coord.transform("size")
lab["n_motifs_cat"] = lab.motifs_in_peak.clip(upper=3).map({1: "1", 2: "2", 3: "3+"})
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.3))
for sd, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
    x = lab[SD == sd]
    for ax, v, xl, lg in ((axes[0], x.peak_width, "Peak width (bp)", True), (axes[1], x.abs_tss + 1, "Distance to nearest TSS (bp)", True)):
        s = np.sort(v.dropna())
        ax.plot(s, np.arange(1, len(s) + 1) / len(s), color=col, lw=1.6, label=f"{sd} (median {np.median(s):,.0f})")
        ax.set_xscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel("Fraction of sites (ECDF)")
for ax in axes[:2]:
    ax.legend(frameon=False, fontsize=8)
share = pd.crosstab(lab.n_motifs_cat, SD, normalize="columns")
ax = axes[2]
for k, (sd, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
    ax.bar(np.arange(3) + (k - 0.5) * 0.36, share[sd].reindex(["1", "2", "3+"]), width=0.34, color=col, label=sd)
ax.set_xticks(range(3), ["1", "2", "3+"])
ax.set_xlabel("Fos motifs in the same peak (CM sites)")
ax.set_ylabel("Share of sites")
ax.legend(frameon=False)
fig.tight_layout()
save(fig, "qc2_site_geometry")
geo = lab.groupby("stable_dynamic")[["peak_width", "abs_tss", "motifs_in_peak"]].median()
geo.to_csv(QC / "tables" / "qc2_geometry_medians.csv")
geo_p = {c: mannwhitneyu(lab[SD == "Stable"][c].dropna(), lab[SD == "Dynamic"][c].dropna()).pvalue for c in ("peak_width", "abs_tss", "motifs_in_peak")}

# ================= QC3 native distributions
rows = []
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.3), sharey=True)
for ax, w in zip(axes, WIDTHS):
    col = f"native_{w}" if w != "peak" else "native_peak"
    bins = np.linspace(-2.5, 1.0, 36)
    for sd, c in (("Stable", BLUE), ("Dynamic", ORANGE)):
        v = lab[SD == sd][col]
        ax.hist(np.log10(v + 1e-3), bins=bins, histtype="step", color=c, lw=1.6, density=True, label=f"{sd}: mean {v.mean():.2f} (slides {SLIDE_NATIVE[sd][w]:.2f})")
        rows.append({"class": sd, "readout": w, "our_mean": v.mean(), "slide_mean": SLIDE_NATIVE[sd][w], "our_median": v.median(), "n": len(v)})
    ax.set_xlabel(f"log10 native predicted ATAC, {LABEL[w].lower()} readout")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
axes[0].set_ylabel("Density")
fig.tight_layout()
save(fig, "qc3_native_distributions")
pd.DataFrame(rows).to_csv(QC / "tables" / "qc3_native_means_vs_slides.csv", index=False)

# ================= QC4 class labels vs measured snATAC by condition
cc = lab.join(cond[["cm_ys", "cm_mas", "cm_mae"]])
cc["age_lfc"] = np.log2((cc.cm_mas + 1) / (cc.cm_ys + 1))
cc["ex_lfc"] = np.log2((cc.cm_mae + 1) / (cc.cm_mas + 1))
fig, axes = plt.subplots(1, 2, figsize=(11, 3.7), sharey=False)
rows = []
for ax, (col, cl, ttl, order_) in zip(axes, (("age_lfc", "age_class", "Measured snATAC change, aged vs young sedentary: log2(MAS / YS)", ["Stable", "Age Up", "Age Down", "Neither YS nor MAS (= MAE only)"]),
                                             ("ex_lfc", "ex_class", "Measured snATAC change, exercise vs sedentary (aged): log2(MAE / MAS)", ["Stable", "Exercise Up", "Exercise Down", "Neither MAS nor MAE (= YS only)"]))):
    data = [cc[cc[cl] == k][col].dropna().values for k in order_]
    bp = ax.boxplot(data, vert=True, widths=0.55, showfliers=False, patch_artist=True, medianprops={"color": INK, "lw": 1.5})
    for patch, k in zip(bp["boxes"], order_):
        patch.set(facecolor=BLUE if k == "Stable" else ORANGE, alpha=0.55, edgecolor=MUTED)
    ax.set_xticks(range(1, 5), [k.split(" (")[0].replace("Neither YS nor MAS", "MAE only").replace("Neither MAS nor MAE", "YS only") + f"\n(n={len(v)})" for k, v in zip(order_, data)], fontsize=8.5)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_ylabel("log2 fold change")
    ax.set_title(ttl, fontsize=9, loc="left")
    for k, v in zip(order_, data):
        rows.append({"label": cl, "class": k, "n": len(v), "median_log2fc": np.median(v)})
    H = kruskal(*data)
    rows.append({"label": cl, "class": "Kruskal-Wallis", "n": sum(len(v) for v in data), "median_log2fc": np.nan, "p": H.pvalue})
fig.tight_layout()
save(fig, "qc4_class_labels_vs_measured")
Q4 = pd.DataFrame(rows)
Q4.to_csv(QC / "tables" / "qc4_class_labels_vs_measured.csv", index=False)

# ================= QC5 prediction vs measurement by class
nm = lab.join(meas[["cm_cpm_per_kb"]])
rows = []
for k in ["All sites", "Stable", "Dynamic"] + order[1:]:
    g = nm if k == "All sites" else nm[nm.stable_dynamic == k] if k in ("Stable", "Dynamic") else nm[nm.binding_category == k]
    r, lo, hi, p, n = rho_ci(g.native_peak, g.cm_cpm_per_kb)
    rows.append({"group": k, "rho": r, "ci_lo": lo, "ci_hi": hi, "n": n})
Q5 = pd.DataFrame(rows)
Q5.to_csv(QC / "tables" / "qc5_prediction_vs_measured_by_class.csv", index=False)
nm["resid"] = nm.native_peak.rank(pct=True) - nm.cm_cpm_per_kb.rank(pct=True)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.6), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
y = np.arange(len(Q5))[::-1]
ax.errorbar(Q5.rho, y, xerr=[Q5.rho - Q5.ci_lo, Q5.ci_hi - Q5.rho], fmt="o", color=BLUE, lw=1.2, ms=5, capsize=0)
ax.set_yticks(y, [f"{g} (n={int(n)})" for g, n in zip(Q5.group, Q5.n)], fontsize=8.5)
ax.set_xlabel("Spearman rho: predicted vs measured accessibility (peak)")
ax = axes[1]
for sd, c in (("Stable", BLUE), ("Dynamic", ORANGE)):
    v = nm[nm.stable_dynamic == sd].resid.dropna()
    ax.hist(v, bins=np.linspace(-1, 1, 41), histtype="step", color=c, lw=1.6, density=True, label=f"{sd} (SD {v.std():.2f})")
ax.set_xlabel("Rank residual (predicted percentile - measured percentile)")
ax.set_ylabel("Density")
ax.legend(frameon=False)
fig.tight_layout()
save(fig, "qc5_prediction_vs_measured_by_class")

# ================= QC6 shuffled vs native (slide-style QC)
rows = []
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4))
for ax, w in zip(axes, WIDTHS):
    col = READ[w]
    for k, sd in enumerate(("Stable", "Dynamic")):
        ids = lab.index[SD == sd]
        n_ = nat.loc[ids, col]
        s_ = {x: d[(d.condition == f"shuffle@{x}") & (d.trackset == "CM") & d.site_id.isin(ids)].groupby("site_id")[col].mean() for x in WIDTHS}
        means = [s_[x].mean() for x in WIDTHS]
        ses = [s_[x].std() / np.sqrt(len(s_[x])) for x in WIDTHS]
        ax.bar(np.arange(3) + (k - 0.5) * 0.36, means, width=0.34, yerr=ses, color=BLUE if sd == "Stable" else ORANGE, alpha=0.8, label=f"{sd}", error_kw={"lw": 1})
        ax.hlines(n_.mean(), -0.5 + k * 0.0, 2.5, color=BLUE if sd == "Stable" else ORANGE, ls="--", lw=1.1)
        for x in WIDTHS:
            rows.append({"readout": w, "shuffle": x, "class": sd, "mean_native": n_.mean(), "mean_shuffled": s_[x].mean(), "frac_sites_decreased": float((s_[x] < n_.loc[s_[x].index]).mean()), "n": len(ids)})
    ax.set_xticks(range(3), ["Motif\nshuffle", "170 bp\nshuffle", "Peak\nshuffle"])
    ax.set_ylabel("Mean predicted ATAC")
    ax.set_title(f"{LABEL[w]} readout (dashed = native)", fontsize=10, loc="left")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "qc6_shuffled_vs_native")
Q6 = pd.DataFrame(rows)
Q6.to_csv(QC / "tables" / "qc6_shuffled_vs_native.csv", index=False)

# ================= QC7 replicate convergence
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5))
rows = []
for w, col in zip(WIDTHS, (ORANGE, BLUE, AQUA)):
    M = REP[w].dropna().values
    full = M.mean(1)
    errs = []
    for k in range(1, 10):
        e = []
        for _ in range(30):
            idx = np.array([rng.permutation(10)[:k] for _ in range(len(M))])
            e.append(np.abs(np.take_along_axis(M, idx, 1).mean(1) - full))
        errs.append(np.median(np.concatenate(e)))
    axes[0].plot(range(1, 10), errs, "o-", color=col, ms=4, lw=1.4, label=f"{LABEL[w]} shuffle")
    rows += [{"shuffle": w, "k_replicates": k, "median_abs_error_vs_10": e_} for k, e_ in zip(range(1, 10), errs)]
    axes[1].hist(M.std(1, ddof=1), bins=np.linspace(0, 2.5, 40), histtype="step", color=col, lw=1.5, density=True, label=f"{LABEL[w]} (median {np.median(M.std(1, ddof=1)):.2f})")
axes[0].set_xlabel("Number of shuffle replicates averaged")
axes[0].set_ylabel("Median |difference| from the 10-replicate effect (log2FC)")
axes[0].legend(frameon=False)
axes[1].set_xlabel("Per-site SD across the 10 replicates (log2FC)")
axes[1].set_ylabel("Density")
axes[1].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "qc7_replicate_convergence")
Q7 = pd.DataFrame(rows)
Q7.to_csv(QC / "tables" / "qc7_replicate_convergence.csv", index=False)

# ================= QC8 GPU / chunk consistency
chunk_gpu = {}
for f in glob.glob(str(ROOT / "Katana" / "logs" / "*shuffle_cm_chunk*.log")):
    m_ = re.search(r"chunk(\d+)of", f)
    g_ = re.search(r"NVIDIA ([A-Za-z0-9\-]+)", open(f, errors="ignore").read())
    if m_ and g_:
        chunk_gpu[int(m_.group(1))] = g_.group(1).replace("A100-SXM4-40GB", "A100")
lab["gpu"] = lab.chunk.map(chunk_gpu)
cs = lab.groupby("chunk").agg(eff170=("eff_170", "median"), effm=("eff_motif", "median"), effp=("eff_peak", "median"), nat=("native_peak", "median"), gpu=("gpu", "first"))
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.3), sharex=True)
cmap = {"L40S": BLUE, "A100": ORANGE, "H100": AQUA}
for ax, (col, ttl) in zip(axes, (("effm", "Motif shuffle (median effect)"), ("eff170", "170 bp shuffle (median effect)"), ("effp", "Peak shuffle (median effect)"))):
    for g_, c_ in cmap.items():
        s = cs[cs.gpu == g_]
        ax.scatter(s.index, s[col], s=22, color=c_, label=f"{g_} (n={len(s)})")
    ax.set_xlabel("Array chunk (0-39; sites interleaved)")
    ax.set_title(ttl, fontsize=10, loc="left")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "qc8_gpu_consistency")
rows = []
for col in ("eff_motif", "eff_170", "eff_peak", "native_peak"):
    groups = [lab[lab.gpu == g_][col].dropna() for g_ in cmap if (lab.gpu == g_).any()]
    rows.append({"measure": col, "kruskal_p_across_gpus": kruskal(*groups).pvalue, **{f"median_{g_}": lab[lab.gpu == g_][col].median() for g_ in cmap}})
Q8 = pd.DataFrame(rows)
Q8.to_csv(QC / "tables" / "qc8_gpu_consistency.csv", index=False)

# ================= AN1 effect distributions
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.3), sharey=True)
rows = []
for ax, w in zip(axes, WIDTHS):
    for sd, c in (("Stable", BLUE), ("Dynamic", ORANGE)):
        v = lab[SD == sd][f"eff_{w}"].dropna()
        ax.hist(v, bins=np.linspace(-4.5, 2, 66), histtype="step", color=c, lw=1.6, density=True, label=f"{sd}: median {v.median():.2f}")
        rows.append({"shuffle": w, "class": sd, "median": v.median(), "q10": v.quantile(.1), "q90": v.quantile(.9), "frac_gt0": float((v > 0).mean()), "frac_lt_minus1": float((v < -1).mean())})
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel(f"Per-site effect, {LABEL[w].lower()} shuffle (log2FC)")
    ax.legend(frameon=False, fontsize=8)
axes[0].set_ylabel("Density")
fig.tight_layout()
save(fig, "an1_effect_distributions")
AN1 = pd.DataFrame(rows)
AN1.to_csv(QC / "tables" / "an1_effect_distributions.csv", index=False)

# ================= AN2 effect vs native accessibility / peak width / TSS distance
fig, axes = plt.subplots(3, 3, figsize=(11.5, 8.6))
for r, w in enumerate(WIDTHS):
    for c, (xc, xl, logx) in enumerate((("native_peak", "Native predicted ATAC (peak)", True), ("peak_width", "Peak width (bp)", True), ("abs_tss", "Distance to nearest TSS (bp)", True))):
        ax = axes[r, c]
        for sd, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
            g = lab[SD == sd]
            bx, by, q1, q3, _ = binned(g[xc] + (1 if xc == "abs_tss" else 0), g[f"eff_{w}"], 8)
            ax.plot(bx, by, "o-", color=col, ms=4, lw=1.4, label=sd)
            ax.fill_between(bx.values, q1.values, q3.values, color=col, alpha=0.12, lw=0)
        ax.set_xscale("log")
        if r == 2:
            ax.set_xlabel(xl)
        if c == 0:
            ax.set_ylabel(f"{LABEL[w]} shuffle effect")
axes[0, 0].legend(frameon=False, fontsize=8)
fig.text(0.01, -0.01, "Lines: median per octile of the x variable; bands: inter-quartile range.", fontsize=8, color=MUTED)
fig.tight_layout()
save(fig, "an2_effect_vs_covariates")
rows = []
for w in WIDTHS:
    for xc in ("native_peak", "peak_width", "abs_tss", "fp_score"):
        r, lo, hi, p, n = rho_ci(lab[xc], lab[f"eff_{w}"])
        rows.append({"shuffle": w, "covariate": xc, "rho": r, "ci_lo": lo, "ci_hi": hi, "p": p, "n": n})
AN2 = pd.DataFrame(rows)
AN2.to_csv(QC / "tables" / "an2_effect_vs_covariates_spearman.csv", index=False)

# ================= AN3 motif vs context dependence (site types)
a, b = lab.eff_motif, lab.eff_170
mot, ctx = a < -1, b < -1
lab["site_type"] = np.select([mot & ctx, mot & ~ctx, ~mot & ctx], ["Motif and context", "Motif only", "Context only"], "Neither")
ST = pd.crosstab(lab.site_type, SD)
chi_p = chi2_contingency(ST.values)[1]
ST_share = pd.crosstab(lab.site_type, SD, normalize="columns")
rows = []
for typ in ST.index:
    t_ = lab[lab.site_type == typ]
    rows.append({"site_type": typ, "n": len(t_), "median_native_peak": t_.native_peak.median(), "median_peak_width": t_.peak_width.median(), "median_abs_tss": t_.abs_tss.median(),
                 "share_stable": float(ST_share.loc[typ, "Stable"]), "share_dynamic": float(ST_share.loc[typ, "Dynamic"])})
AN3 = pd.DataFrame(rows)
AN3.to_csv(QC / "tables" / "an3_site_types.csv", index=False)
r_mc = rho_ci(a, b)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]
hb = ax.hexbin(a, b, gridsize=45, cmap="Blues", mincnt=1, bins="log", extent=(-4, 2, -5, 2))
ax.axvline(-1, color=MUTED, lw=0.8, ls=":")
ax.axhline(-1, color=MUTED, lw=0.8, ls=":")
ax.set_xlabel("Motif shuffle effect (motif readout, log2FC)")
ax.set_ylabel("170 bp shuffle effect (170 bp readout, log2FC)")
ax.set_title(f"Spearman rho = {r_mc[0]:.2f}; dotted lines: 50% loss", fontsize=9.5, loc="left")
fig.colorbar(hb, ax=ax, label="Sites (log10)")
ax = axes[1]
types = ["Motif and context", "Motif only", "Context only", "Neither"]
cols_ = [BLUE, AQUA, ORANGE, GREY]
bottom = np.zeros(2)
for tname, c_ in zip(types, cols_):
    v = np.array([ST_share.loc[tname, "Stable"], ST_share.loc[tname, "Dynamic"]])
    ax.bar([0, 1], v, bottom=bottom, color=c_, width=0.6, label=tname)
    for i in range(2):
        ax.text(i, bottom[i] + v[i] / 2, f"{v[i]:.0%}", ha="center", va="center", fontsize=9, color="white" if tname != "Neither" else INK)
    bottom += v
ax.set_xticks([0, 1], ["Stable", "Dynamic"])
ax.set_ylabel("Share of sites")
ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
ax.set_title(f"Site types (chi-square p = {chi_p:.1e})", fontsize=9.5, loc="left")
fig.tight_layout()
save(fig, "an3_site_types")

# ================= AN4 footprint score vs predicted dependence
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4))
rows = []
for ax, (yc, yl) in zip(axes, (("eff_motif", "Motif shuffle effect"), ("eff_170", "170 bp shuffle effect"), ("native_peak", "Native predicted ATAC (peak)"))):
    for sd, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
        g = lab[SD == sd]
        bx, by, q1, q3, _ = binned(g.fp_score, g[yc], 6)
        ax.plot(bx, by, "o-", color=col, ms=4, lw=1.4, label=sd)
        ax.fill_between(bx.values, q1.values, q3.values, color=col, alpha=0.12, lw=0)
        r, lo, hi, p, n = rho_ci(g.fp_score, g[yc])
        rows.append({"outcome": yc, "class": sd, "rho": r, "ci_lo": lo, "ci_hi": hi, "p": p, "n": n})
    r, lo, hi, p, n = rho_ci(lab.fp_score, lab[yc])
    rows.append({"outcome": yc, "class": "All", "rho": r, "ci_lo": lo, "ci_hi": hi, "p": p, "n": n})
    ax.set_xlabel("TOBIAS footprint score (best bound condition)")
    ax.set_ylabel(yl)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "an4_footprint_score_vs_dependence")
AN4 = pd.DataFrame(rows)
AN4.to_csv(QC / "tables" / "an4_footprint_score.csv", index=False)

# ================= AN5 several Fos motifs in one peak
rows = []
for w in ("motif", "170", "peak"):
    one = lab[lab.motifs_in_peak == 1][f"eff_{w}"].dropna()
    multi = lab[lab.motifs_in_peak >= 2][f"eff_{w}"].dropna()
    y = lab[[f"eff_{w}", "log_nat", "log_pw", "motifs_in_peak"]].dropna()
    X = np.column_stack([np.ones(len(y)), (y.motifs_in_peak >= 2).astype(float), y.log_nat - y.log_nat.mean(), y.log_pw - y.log_pw.mean()])
    beta, se, pv = ols_hc3(X, y[f"eff_{w}"])
    lo, hi = np.percentile([np.median(rng.choice(multi, len(multi))) - np.median(rng.choice(one, len(one))) for _ in range(1000)], [2.5, 97.5])
    rows.append({"shuffle": w, "n_single": len(one), "n_multi": len(multi), "median_single": one.median(), "median_multi": multi.median(), "raw_diff": multi.median() - one.median(),
                 "raw_lo": lo, "raw_hi": hi, "raw_p": mannwhitneyu(multi, one).pvalue, "adj_diff": beta[1], "adj_lo": beta[1] - 1.96 * se[1], "adj_hi": beta[1] + 1.96 * se[1], "adj_p": pv[1]})
AN5 = pd.DataFrame(rows)
AN5["raw_fdr"], AN5["adj_fdr"] = bh(AN5.raw_p), bh(AN5.adj_p)
AN5.to_csv(QC / "tables" / "an5_multiple_motifs_in_peak.csv", index=False)
fig, ax = plt.subplots(figsize=(6.8, 3.6))
for k, w in enumerate(WIDTHS):
    r_ = AN5.iloc[k]
    ax.errorbar([k - 0.12], [r_.median_single], yerr=[[r_.median_single - boot_ci(lab[lab.motifs_in_peak == 1][f"eff_{w}"])[0]], [boot_ci(lab[lab.motifs_in_peak == 1][f"eff_{w}"])[1] - r_.median_single]], fmt="o", color=BLUE, ms=5, lw=1.2, capsize=0, label="One Fos motif in the peak" if k == 0 else None)
    ax.errorbar([k + 0.12], [r_.median_multi], yerr=[[r_.median_multi - boot_ci(lab[lab.motifs_in_peak >= 2][f"eff_{w}"])[0]], [boot_ci(lab[lab.motifs_in_peak >= 2][f"eff_{w}"])[1] - r_.median_multi]], fmt="o", color=ORANGE, ms=5, lw=1.2, capsize=0, label="Two or more" if k == 0 else None)
ax.set_xticks(range(3), ["Motif shuffle", "170 bp shuffle", "Peak shuffle"])
ax.set_ylabel("Median effect (95% CI)")
ax.legend(frameon=False)
save(fig, "an5_multiple_motifs")

# ================= AN6 promoter-proximal vs distal
lab["tss_bin"] = pd.cut(lab.abs_tss, [-1, 2000, 10000, 30000, 1e9], labels=["<=2 kb", "2-10 kb", "10-30 kb", ">30 kb"])
rows = []
fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.3))
for ax, (yc, ttl) in zip(axes, (("native_peak", "Native predicted ATAC (peak)"), ("eff_motif", "Motif shuffle effect"), ("eff_170", "170 bp shuffle effect"), ("eff_peak", "Peak shuffle effect"))):
    for k, (sd, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
        g = lab[SD == sd].groupby("tss_bin", observed=True)[yc].median()
        ax.plot(np.arange(4) + (k - 0.5) * 0.1, g.reindex(["<=2 kb", "2-10 kb", "10-30 kb", ">30 kb"]).values, "o-", color=col, ms=4, lw=1.4, label=sd)
        for tb, v in g.items():
            rows.append({"measure": yc, "class": sd, "tss_bin": tb, "median": v, "n": int(((SD == sd) & (lab.tss_bin == tb)).sum())})
    ax.set_xticks(range(4), ["<=2 kb", "2-10", "10-30", ">30 kb"], fontsize=8.5)
    ax.set_xlabel("Distance to nearest TSS")
    ax.set_title(ttl, fontsize=9.5, loc="left")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "an6_tss_proximity")
AN6 = pd.DataFrame(rows)
AN6.to_csv(QC / "tables" / "an6_tss_proximity.csv", index=False)
prox_p = {w: mannwhitneyu(lab[lab.tss_bin == "<=2 kb"][f"eff_{w}"].dropna(), lab[lab.tss_bin == ">30 kb"][f"eff_{w}"].dropna()).pvalue for w in WIDTHS}

# ================= AN7 MPRA
mfiles = sorted(glob.glob(str(ROOT / "Katana" / "results" / "mpra_chunk*.csv")))
if mfiles:
    mp = pd.concat([pd.read_csv(f) for f in mfiles], ignore_index=True).drop_duplicates(["mpra_id", "trackset", "condition", "rep"])
    mt = pd.read_csv(ROOT / "Katana" / "data" / "MPRA_Results_with_TSS_distance.csv")
    mt = mt.drop_duplicates("Coords").set_index("Coords")
    rows, panels = [], {}
    for ts, tlabel in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
        g = mp[mp.trackset == ts]
        natv = g[g.condition == "native"].groupby("mpra_id").atac_170.mean()
        pc = 0.05 * natv.median()
        info = g.drop_duplicates("mpra_id").set_index("mpra_id").motif_source
        res = {"native": natv}
        for cn in ("motif_shuffle", "full_shuffle"):
            s_ = g[g.condition == cn].groupby("mpra_id").atac_170.mean()
            res[cn] = np.log2((s_ + pc) / (natv.reindex(s_.index) + pc))
        panels[ts] = res
        tests = (("native predicted ATAC vs log2 alpha (native)", res["native"], mt.log2_alpha_Endo, None),
                 ("motif-shuffle effect vs MPRA (Endo - motif-shuffled)", res["motif_shuffle"], mt.delta_vs_MotifShuff, None),
                 ("full-shuffle effect vs MPRA (Endo - fully shuffled)", res["full_shuffle"], mt.delta_vs_FullShuff, None))
        for name, x, y, _ in tests:
            for subset, ids in (("all elements", x.index), ("replicated in MPRA", [i for i in x.index if i in mt.index and mt.loc[i, "replicated"]]),
                                ("motif from TOBIAS", [i for i in x.index if info.get(i) == "tobias"])):
                ids = [i for i in ids if i in mt.index]
                r, lo, hi, p, n = rho_ci(x.reindex(ids), y.reindex(ids))
                rows.append({"track": tlabel, "comparison": name, "subset": subset, "rho": r, "ci_lo": lo, "ci_hi": hi, "p": p, "n": n})
    AN7 = pd.DataFrame(rows)
    AN7.to_csv(QC / "tables" / "an7_mpra_comparison.csv", index=False)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
    for r_, (ts, tlabel) in enumerate((("FB", "Patski ATAC"), ("CM", "heart ATAC"))):
        res = panels[ts]
        for c_, (key, y, xl, yl) in enumerate(((("native"), mt.log2_alpha_Endo, "Predicted native ATAC (170 bp)", "MPRA log2 alpha (native)"),
                                                ("motif_shuffle", mt.delta_vs_MotifShuff, "Predicted motif-shuffle effect (log2FC)", "MPRA: native - motif-shuffled"),
                                                ("full_shuffle", mt.delta_vs_FullShuff, "Predicted full-shuffle effect (log2FC)", "MPRA: native - fully shuffled"))):
            ax = axes[r_, c_]
            x = res[key]
            ids = [i for i in x.index if i in mt.index]
            xx, yy = x.reindex(ids), y.reindex(ids)
            ax.scatter(xx, yy, s=6, alpha=0.3, color=BLUE, linewidths=0)
            if key == "native":
                ax.set_xscale("log")
            r, lo, hi, p, n = rho_ci(xx, yy)
            ax.set_title(f"{tlabel}: rho = {r:+.2f} [{lo:+.2f}, {hi:+.2f}], n = {n}", fontsize=9, loc="left")
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
    fig.text(0.01, -0.01, "MPRA deltas are native minus shuffled, so a predicted loss (negative log2FC) should give a POSITIVE delta: the expected correlation for the shuffle panels is negative.", fontsize=8, color=MUTED)
    fig.tight_layout()
    save(fig, "an7_mpra")
    print(AN7.round(3).to_string(index=False))

pd.set_option("display.width", 220, "display.max_columns", 30)
lab.drop(columns=["in_YS", "in_MAS", "in_MAE"]).to_csv(QC / "tables" / "site_level_table_cm.csv")
print("\nQC4\n", Q4.round(3).to_string(index=False))
print("\nQC5\n", Q5.round(3).to_string(index=False))
print("\nQC6\n", Q6.round(3).to_string(index=False))
print("\nQC7 (k=1,3,5,9)\n", Q7[Q7.k_replicates.isin([1, 3, 5, 9])].round(3).to_string(index=False))
print("\nQC8\n", Q8.round(4).to_string(index=False))
print("\nAN1\n", AN1.round(3).to_string(index=False))
print("\nAN2\n", AN2.round(3).to_string(index=False))
print("\nAN3\n", AN3.round(3).to_string(index=False), "\nchi-square p", chi_p)
print("\nAN4\n", AN4.round(3).to_string(index=False))
print("\nAN5\n", AN5.round(3).to_string(index=False))
print("\nAN6 tss proximity p (<=2kb vs >30kb):", {k: f"{v:.1e}" for k, v in prox_p.items()})
print("\nQC2 medians\n", geo.round(1).to_string(), "\np:", {k: f"{v:.1e}" for k, v in geo_p.items()})
print("\nQC3\n", pd.DataFrame(rows if False else pd.read_csv(QC / "tables" / "qc3_native_means_vs_slides.csv")).round(3).to_string(index=False))
