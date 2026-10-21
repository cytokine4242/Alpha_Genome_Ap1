#!/usr/bin/env python3
"""Statistics for the shuffle arm: class comparisons with significance marks, adjusted models, robustness and noise checks.
Reads Katana/results/shuffle_cm_chunk*.csv, writes Analysis/cm/stats/*.csv and the starred figures (fig2s, fig3s, fig6, fig7, fig8).
Used by 25_build_report.py.   python3 24_class_stats.py

Per site, the effect is the mean over the 10 replicates of log2((shuffled+pc)/(native+pc)) on the matched readout (shuffle width = readout region).
Tests: two-sided Mann-Whitney vs the reference class, Benjamini-Hochberg within each family; adjusted effects from OLS
(effect ~ class + log10 native ATAC + log10 peak width) with HC3 robust standard errors; Kruskal-Wallis omnibus; label-permutation test for the
headline contrast.  '*' FDR<0.05, '**' <0.01, '***' <0.001, 'ns' otherwise.
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import kruskal, mannwhitneyu, t as tdist

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "cm"
(OUT / "stats").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, "#ececea", "#e34948"])
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
WIDTHS = ["motif", "170", "peak"]
READ = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
LABEL = {"motif": "Motif", "170": "170 bp", "peak": "Peak"}


def stars(q):
    return "" if pd.isna(q) else "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_median_diff(a, b, n=2000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.median(rng.choice(a, (n, len(a))), 1) - np.median(rng.choice(b, (n, len(b))), 1)
    return np.percentile(d, [2.5, 97.5])


def cliffs(a, b):
    u = mannwhitneyu(a, b).statistic
    return 2 * u / (len(a) * len(b)) - 1  # P(a > b) - P(a < b)


def ols_hc3(X, y):
    """OLS with HC3 robust SEs -> beta, se, p (t with n-k df)."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    beta = XtXi @ X.T @ y
    e = y - X @ beta
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    meat = X.T @ (X * (e / (1 - h))[:, None] ** 2)
    cov = XtXi @ meat @ XtXi
    se = np.sqrt(np.diag(cov))
    return beta, se, 2 * tdist.sf(np.abs(beta / se), n - k)


def effects(df, ts, col, pc_frac=0.05):
    g = df[df.trackset == ts]
    nat = g[g.condition == "native"].set_index("site_id")[col]
    pc = pc_frac * nat.median()
    e = g[g.condition.str.startswith("shuffle")].copy()
    e["lfc"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    e["width"] = e.condition.str.split("@").str[1]
    return e.groupby(["site_id", "width"]).lfc.agg(["mean", "std"]).unstack(), nat


files = sorted(glob.glob(str(ROOT / "Katana" / "results" / "shuffle_cm_chunk*.csv")))
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
cmsites = pd.read_csv(ROOT / "Katana" / "data" / "cm_fos_sites.csv").set_index("site_id")
lab = d.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic", "age_class", "ex_class", "binding_category", "tss_dist"]].join(cmsites[["peak_width"]])
lab["age_class"] = lab.age_class.replace({"Unbound": "Neither YS nor MAS (= MAE only)"})
lab["ex_class"] = lab.ex_class.replace({"Unbound": "Neither MAS nor MAE (= YS only)"})

EFF, NAT = {}, {}
for ts in ("CM", "FB"):
    for w in WIDTHS:
        e, nat = effects(d, ts, READ[w])
        EFF[(ts, w)] = e["mean"][w]
        NAT[(ts, w)] = nat
        EFF[(ts, w, "sd")] = e["std"][w]
nat_cm = NAT[("CM", "peak")]
cov = pd.DataFrame({"log_nat": np.log10(nat_cm + 1e-3), "log_pw": np.log10(lab.peak_width)}).join(lab)

# ============ A. Stable vs Dynamic: every shuffle x readout, two tracksets
rows = []
for ts in ("CM", "FB"):
    for sh in WIDTHS:
        for rd in WIDTHS:
            e, _ = effects(d, ts, READ[rd])
            x = e["mean"][sh].to_frame("y").join(lab.stable_dynamic).dropna()
            s, dy = x[x.stable_dynamic == "Stable"].y.values, x[x.stable_dynamic == "Dynamic"].y.values
            lo, hi = boot_median_diff(dy, s)
            rows.append({"trackset": ts, "shuffle": sh, "readout": rd, "n_stable": len(s), "n_dynamic": len(dy), "median_stable": np.median(s), "median_dynamic": np.median(dy),
                         "diff_dyn_minus_stable": np.median(dy) - np.median(s), "ci_lo": lo, "ci_hi": hi, "cliffs_delta": cliffs(dy, s), "p": mannwhitneyu(dy, s).pvalue})
A = pd.DataFrame(rows)
for ts in A.trackset.unique():
    A.loc[A.trackset == ts, "fdr"] = bh(A.loc[A.trackset == ts, "p"])
A["sig"] = A.fdr.map(stars)
A.to_csv(OUT / "stats" / "stable_vs_dynamic_all_readouts.csv", index=False)

# permutation test for the headline contrasts (matched readouts, CM tracks)
rows = []
for w in WIDTHS:
    x = EFF[("CM", w)].to_frame("y").join(lab.stable_dynamic).dropna()
    y, g = x.y.values, (x.stable_dynamic == "Dynamic").values
    obs = np.median(y[g]) - np.median(y[~g])
    perm = np.array([np.median(y[(pg := rng.permutation(g))]) - np.median(y[~pg]) for _ in range(5000)])
    rows.append({"width": w, "observed_diff": obs, "perm_p": (np.sum(np.abs(perm) >= abs(obs)) + 1) / 5001, "perm_null_sd": perm.std()})
PERM = pd.DataFrame(rows)
PERM.to_csv(OUT / "stats" / "permutation_stable_vs_dynamic.csv", index=False)

# ============ B. classes vs reference: raw (MWU) and adjusted (OLS HC3), matched readout
rows, omni = [], []
FAMS = (("age", "age_class", "Stable"), ("exercise", "ex_class", "Stable"), ("binding", "binding_category", "YS+MAS+MAE"))
for fam, col, ref in FAMS:
    for w in WIDTHS:
        y = EFF[("CM", w)].to_frame("y").join(cov).dropna(subset=["y", "log_nat", "log_pw", col])
        classes = [c for c in y[col].unique() if c != ref]
        grp = [y[y[col] == c].y.values for c in [ref] + classes]
        omni.append({"family": fam, "width": w, "kruskal_H": kruskal(*grp).statistic, "kruskal_p": kruskal(*grp).pvalue})
        Xd = pd.get_dummies(y[col]).drop(columns=ref).astype(float)
        X = np.column_stack([np.ones(len(y)), Xd.values, (y.log_nat - y.log_nat.mean()).values, (y.log_pw - y.log_pw.mean()).values])
        beta, se, pv = ols_hc3(X, y.y.values)
        names = list(Xd.columns)
        base = y[y[col] == ref].y
        for i, c in enumerate(["Stable-like reference"] + names):
            if i == 0:
                rows.append({"family": fam, "width": w, "class": ref, "n": len(base), "median": base.median(), "raw_diff": 0.0, "raw_p": np.nan, "adj_diff": 0.0, "adj_lo": np.nan, "adj_hi": np.nan, "adj_p": np.nan})
                continue
            xx = y[y[col] == c].y
            lo, hi = boot_median_diff(xx, base)
            rows.append({"family": fam, "width": w, "class": c, "n": len(xx), "median": xx.median(), "raw_diff": xx.median() - base.median(), "raw_lo": lo, "raw_hi": hi,
                         "raw_p": mannwhitneyu(xx, base).pvalue if len(xx) >= 5 else np.nan, "cliffs_delta": cliffs(xx.values, base.values),
                         "adj_diff": beta[i], "adj_lo": beta[i] - 1.96 * se[i], "adj_hi": beta[i] + 1.96 * se[i], "adj_p": pv[i]})
        # coefficients of the covariates (reported once per model)
        rows.append({"family": fam, "width": w, "class": "[log10 native ATAC, per 1]", "n": len(y), "adj_diff": beta[-2], "adj_lo": beta[-2] - 1.96 * se[-2], "adj_hi": beta[-2] + 1.96 * se[-2], "adj_p": pv[-2]})
        rows.append({"family": fam, "width": w, "class": "[log10 peak width, per 1]", "n": len(y), "adj_diff": beta[-1], "adj_lo": beta[-1] - 1.96 * se[-1], "adj_hi": beta[-1] + 1.96 * se[-1], "adj_p": pv[-1]})
B = pd.DataFrame(rows)
ok = ~B["class"].str.startswith("[") & B.raw_p.notna()
for fam in B.family.unique():
    ix = ok & (B.family == fam)
    B.loc[ix, "raw_fdr"] = bh(B.loc[ix, "raw_p"])
    B.loc[ix, "adj_fdr"] = bh(B.loc[ix, "adj_p"])
B["raw_sig"], B["adj_sig"] = B.raw_fdr.map(stars), B.adj_fdr.map(stars)
B.to_csv(OUT / "stats" / "class_vs_reference.csv", index=False)
OM = pd.DataFrame(omni)
OM["sig"] = OM.kruskal_p.map(stars)
OM.to_csv(OUT / "stats" / "omnibus_kruskal.csv", index=False)

# ============ C. stratified by native accessibility (is the Stable-Dynamic gap just a floor effect?)
y = EFF[("CM", "170")].to_frame("y").join(cov).dropna(subset=["y", "log_nat"])
y["bin"] = pd.qcut(y.log_nat, 5, labels=False, duplicates="drop")
rows = []
for b, g in y.groupby("bin"):
    s, dy = g[g.stable_dynamic == "Stable"].y, g[g.stable_dynamic == "Dynamic"].y
    rows.append({"native_quintile": int(b) + 1, "native_atac_median": 10 ** g.log_nat.median() - 1e-3, "n_stable": len(s), "n_dynamic": len(dy), "median_stable": s.median(), "median_dynamic": dy.median(),
                 "diff": dy.median() - s.median(), "p": mannwhitneyu(dy, s).pvalue if min(len(s), len(dy)) >= 5 else np.nan})
C = pd.DataFrame(rows)
C["fdr"] = bh(C.p.fillna(1))
C["sig"] = C.fdr.map(stars)
C.to_csv(OUT / "stats" / "stable_vs_dynamic_by_native_quintile.csv", index=False)

# ============ D. noise: how many sites have an effect distinguishable from the replicate scatter?
rows = []
for w in WIDTHS:
    m, sd = EFF[("CM", w)], EFF[("CM", w, "sd")]
    se = sd / np.sqrt(10)
    z = (m / se).replace([np.inf, -np.inf], np.nan).dropna()
    rows.append({"shuffle": w, "median_effect": m.median(), "median_replicate_sd": sd.median(), "median_se_10_reps": se.median(),
                 "frac_sites_effect_gt_2se": float((z.abs() > 2).mean()), "frac_sites_negative_gt_2se": float((z < -2).mean())})
D = pd.DataFrame(rows)
D.to_csv(OUT / "stats" / "replicate_noise_detectability.csv", index=False)

# ============ figures with significance marks
vmax = 1.95
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.3), sharey=True)
a = A[A.trackset == "CM"]
for ax, sd in zip(axes, ("Stable", "Dynamic")):
    M = np.array([[a[(a.shuffle == sh) & (a.readout == rd)][f"median_{sd.lower()}"].iloc[0] for rd in WIDTHS] for sh in WIDTHS])
    im = ax.imshow(M, cmap=DIV, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(3), ["Motif", "170 bp", "Peak"])
    ax.set_yticks(range(3), ["Motif shuffle", "170 bp shuffle", "Peak shuffle"])
    ax.set_xlabel("ATAC readout region")
    ax.set_title(f"{sd} (n = {int(a['n_' + sd.lower()].iloc[0])})", fontsize=10, loc="left")
    for i, sh in enumerate(WIDTHS):
        for j, rd in enumerate(WIDTHS):
            txt = f"{M[i, j]:.2f}"
            if sd == "Dynamic":
                txt += "\n" + a[(a.shuffle == sh) & (a.readout == rd)].sig.iloc[0]
            ax.text(j, i, txt, ha="center", va="center", fontsize=9.5, color=INK if abs(M[i, j]) < 0.6 * vmax else "white")
fig.colorbar(im, ax=axes, shrink=0.85, label="Median log2FC (shuffled / native)")
fig.text(0.01, -0.04, "Marks on the Dynamic panel: Mann-Whitney vs Stable, BH-corrected over the 9 cells.  *** FDR<0.001  ** <0.01  * <0.05  ns not significant", fontsize=8, color=MUTED)
fig.savefig(OUT / "figures" / "fig2s_shuffle_heatmaps_stars.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(3, 3, figsize=(10.5, 8.2), sharex="col")
for r, (fam, col, ref) in enumerate(FAMS):
    x = B[(B.family == fam) & ~B["class"].str.startswith("[")]
    order = list(dict.fromkeys(x["class"]))
    for c, w in enumerate(WIDTHS):
        ax = axes[r, c]
        yy = x[x.width == w].set_index("class").loc[order]
        pos = np.arange(len(order))[::-1]
        for pi, k in zip(pos, order):
            v = yy.loc[k]
            xs = EFF[("CM", w)].to_frame("y").join(lab)
            vals = xs[xs[col] == k].y.dropna()
            lo, hi = np.percentile([np.median(rng.choice(vals, len(vals))) for _ in range(500)], [2.5, 97.5])
            ax.errorbar(v["median"], pi, xerr=[[v["median"] - lo], [hi - v["median"]]], fmt="o", color=BLUE if k == ref else ORANGE, ecolor=BLUE if k == ref else ORANGE, ms=5, lw=1.2, capsize=0)
        ax.set_yticks(pos, [f"{k} (n={int(n)})" for k, n in zip(order, yy.n)], fontsize=8)
        ax.axvline(yy.loc[ref, "median"], color=BLUE, lw=0.8, ls=":")
        xmax = ax.get_xlim()[1]
        for pi, k in zip(pos, order):
            if k != ref:
                ax.text(1.02, pi, yy.loc[k, "raw_sig"], transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=INK)
        if r == 0:
            ax.set_title(f"{LABEL[w]} shuffle, {LABEL[w].lower()} readout", fontsize=10, loc="left")
        if r == 2:
            ax.set_xlabel("Median log2FC (95% CI)")
    axes[r, 0].annotate({"age": "Age (YS vs MAS)", "exercise": "Exercise (MAS vs MAE)", "binding": "Binding pattern"}[fam], (0, 1.13), xycoords="axes fraction", fontsize=9, color=MUTED)
fig.text(0.01, -0.01, "Blue = reference (Stable). Marks: Mann-Whitney vs the reference, BH-corrected within each family.  *** FDR<0.001  ** <0.01  * <0.05  ns not significant", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "fig3s_class_effects_stars.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(3, 3, figsize=(10.5, 8.2), sharex="col")
for r, (fam, col, ref) in enumerate(FAMS):
    x = B[(B.family == fam) & ~B["class"].str.startswith("[") & (B["class"] != ref)]
    order = list(dict.fromkeys(x["class"]))
    for c, w in enumerate(WIDTHS):
        ax = axes[r, c]
        yy = x[x.width == w].set_index("class").loc[order]
        pos = np.arange(len(order))[::-1]
        raw = yy.raw_diff
        ax.errorbar(yy.adj_diff, pos, xerr=[yy.adj_diff - yy.adj_lo, yy.adj_hi - yy.adj_diff], fmt="o", color=BLUE, ms=5, lw=1.2, capsize=0, label="Adjusted (native ATAC, peak width)")
        ax.scatter(raw, pos, marker="D", s=22, facecolor="none", edgecolor=ORANGE, lw=1.2, label="Raw difference of medians")
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_yticks(pos, [f"{k} (n={int(n)})" for k, n in zip(order, yy.n)], fontsize=8)
        for pi, k in zip(pos, order):
            ax.text(1.02, pi, yy.loc[k, "adj_sig"], transform=ax.get_yaxis_transform(), va="center", fontsize=8)
        if r == 0:
            ax.set_title(f"{LABEL[w]} shuffle", fontsize=10, loc="left")
        if r == 2:
            ax.set_xlabel("Difference vs reference (log2FC)")
    axes[r, 0].annotate({"age": "Age (YS vs MAS)", "exercise": "Exercise (MAS vs MAE)", "binding": "Binding pattern"}[fam], (0, 1.13), xycoords="axes fraction", fontsize=9, color=MUTED)
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8, loc="upper right", ncol=2)
fig.text(0.01, -0.01, "Positive = less accessibility lost than the reference. Marks: adjusted effect, HC3 robust SE, BH-corrected within each family.", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "fig7_class_effects_adjusted.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
for sd, c, off in (("stable", BLUE, -0.06), ("dynamic", ORANGE, 0.06)):
    ax.plot(C.native_quintile + off, C[f"median_{sd}"], "o-", color=c, lw=1.2, ms=5, label=sd.capitalize())
for _, r in C.iterrows():
    ax.text(r.native_quintile, min(r.median_stable, r.median_dynamic) - 0.1, r.sig, ha="center", fontsize=8, va="top")
ax.set_xticks(C.native_quintile, [f"Q{int(q)}\n({a:.2f})" for q, a in zip(C.native_quintile, C.native_atac_median)])
ax.set_xlabel("Quintile of native predicted ATAC over the peak (median value)")
ax.set_ylabel("Median log2FC, 170 bp shuffle (170 bp readout)")
ax.legend(frameon=False)
fig.savefig(OUT / "figures" / "fig8_stable_vs_dynamic_by_native_accessibility.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 220, "display.max_columns", 30)
print("A. Stable vs Dynamic (CM tracks), matched cells:")
print(A[(A.trackset == "CM") & (A.shuffle == A.readout)][["shuffle", "n_stable", "n_dynamic", "median_stable", "median_dynamic", "diff_dyn_minus_stable", "ci_lo", "ci_hi", "cliffs_delta", "fdr", "sig"]].round(3).to_string(index=False))
print("\nA2. same for the FB (Patski) trackset:")
print(A[(A.trackset == "FB") & (A.shuffle == A.readout)][["shuffle", "median_stable", "median_dynamic", "diff_dyn_minus_stable", "cliffs_delta", "fdr", "sig"]].round(3).to_string(index=False))
print("\nPermutation:\n", PERM.round(4).to_string(index=False))
print("\nOmnibus Kruskal:\n", OM.round(4).to_string(index=False))
print("\nB. adjusted vs raw (age, matched):")
print(B[(B.family.isin(["age", "exercise"])) & ~B["class"].str.startswith("[") & B.raw_p.notna()][["family", "width", "class", "n", "raw_diff", "raw_sig", "adj_diff", "adj_lo", "adj_hi", "adj_sig"]].round(3).to_string(index=False))
print("\nC. strata:\n", C.round(3).to_string(index=False))
print("\nD. noise:\n", D.round(3).to_string(index=False))
print("\ncovariate coefficients:\n", B[B["class"].str.startswith("[") & (B.family == "age")][["width", "class", "adj_diff", "adj_lo", "adj_hi", "adj_p"]].round(3).to_string(index=False))

# ============ E. direct contrasts between dynamic classes (not only vs Stable)
rows = []
for name, col, a_, b_ in (("Age Up vs Age Down", "age_class", "Age Up", "Age Down"), ("Exercise Up vs Exercise Down", "ex_class", "Exercise Up", "Exercise Down"),
                          ("Age Up vs Exercise Up", None, None, None)):
    for w in WIDTHS:
        y = EFF[("CM", w)].to_frame("y").join(lab)
        if col:
            xa, xb = y[y[col] == a_].y.dropna(), y[y[col] == b_].y.dropna()
        else:  # age gain only (MAS only) vs exercise gain only (MAE only)
            xa, xb = y[y.binding_category == "MAS only"].y.dropna(), y[y.binding_category == "MAE only"].y.dropna()
            name = "Age gain only (MAS only) vs exercise gain only (MAE only)"
        lo, hi = boot_median_diff(xa, xb)
        rows.append({"contrast": name, "width": w, "n_a": len(xa), "n_b": len(xb), "median_a": xa.median(), "median_b": xb.median(), "diff_a_minus_b": xa.median() - xb.median(),
                     "ci_lo": lo, "ci_hi": hi, "cliffs_delta": cliffs(xa.values, xb.values), "p": mannwhitneyu(xa, xb).pvalue})
E_ = pd.DataFrame(rows)
E_["fdr"] = bh(E_.p)
E_["sig"] = E_.fdr.map(stars)
E_.to_csv(OUT / "stats" / "between_dynamic_classes.csv", index=False)
print("\nE. between dynamic classes:\n", E_.round(3).to_string(index=False))
