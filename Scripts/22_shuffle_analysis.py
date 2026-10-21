#!/usr/bin/env python3
"""Analysis of the shuffle arm (results/shuffle_{cm,fb}_chunk*.csv): effect of shuffling the motif / 170 bp / peak on predicted ATAC
and RNA, by Stable / Dynamic, age class, exercise class and binding pattern; native prediction vs measured snATAC; confidence tiers.

Effect = per site, mean over the 10 shuffle replicates of log2((shuffled + pc) / (native + pc)), pc = 5% of the median native value of that
readout and trackset. The "matched" readout is the one on the shuffled region (motif shuffle -> motif readout, 170 -> 170, peak -> peak).
Tables go to Analysis/<sites>/tables, figures to Analysis/<sites>/figures.

  python3 22_shuffle_analysis.py cm        (or fb, when its results exist)
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[1]
SITES = sys.argv[1] if len(sys.argv) > 1 else "cm"
OUT = ROOT / "Analysis" / SITES
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(1)

# ---- style (validated default palette: slot 1 blue = Stable, slot 2 orange = Dynamic; diverging blue-gray-red)
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, RED = "#2a78d6", "#eb6834", "#e34948"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, "#ececea", RED])
plt.rcParams.update({"font.size": 10, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})

WIDTHS = ["motif", "170", "peak"]
READ = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
LABEL = {"motif": "Motif", "170": "170 bp", "peak": "Peak"}

# ---- data
files = sorted(glob.glob(str(ROOT / "Katana" / "results" / f"shuffle_{SITES}_chunk*.csv")))
if not files:
    sys.exit(f"no results/shuffle_{SITES}_chunk*.csv")
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
meas = meas[meas.site_table == SITES].set_index("site_id")
MCOL = "cm_cpm_per_kb" if SITES == "cm" else "fb_cpm_per_kb"
sites = d.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic", "age_class", "ex_class", "binding_category", "tss_dist"]]
print(f"{SITES.upper()}: {d.site_id.nunique()} sites from {len(files)} chunks, tracksets {sorted(d.trackset.unique())}")


def effects(df, ts, col, pc_frac=0.05):
    """site x shuffle-width table of mean log2FC for one trackset and readout column."""
    g = df[df.trackset == ts]
    nat = g[g.condition == "native"].set_index("site_id")[col]
    pc = pc_frac * nat.median()
    e = g[g.condition.str.startswith("shuffle")].copy()
    e["lfc"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    e["width"] = e.condition.str.split("@").str[1]
    return e.groupby(["site_id", "width"]).lfc.mean().unstack()


def boot_ci(x, n=1000):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 5:
        return np.nan, np.nan
    m = np.median(rng.choice(x, (n, len(x)))[:, :], axis=1)
    return np.percentile(m, 2.5), np.percentile(m, 97.5)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def save(fig, name):
    fig.savefig(OUT / "figures" / name, bbox_inches="tight")
    plt.close(fig)


# ================= 1. native prediction vs measured snATAC
print("\n=== 1. native prediction vs measured snATAC (Spearman)")
rows = []
for ts in sorted(d.trackset.unique()):
    n = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id").join(meas[[MCOL]])
    for c in ("atac_motif", "atac_170", "atac_peak"):
        m = n[c].notna() & n[MCOL].notna()
        r, p = spearmanr(n.loc[m, c], n.loc[m, MCOL])
        rows.append({"trackset": ts, "readout": c, "rho": r, "p": p, "n": int(m.sum())})
t1 = pd.DataFrame(rows)
t1.to_csv(OUT / "tables" / "native_vs_measured.csv", index=False)
print(t1.round(3).to_string(index=False))

nat = d[(d.condition == "native") & (d.trackset == "CM")].set_index("site_id").join(meas[[MCOL]])
fig, ax = plt.subplots(figsize=(5.2, 4.4))
for sd, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
    x = nat[nat.stable_dynamic == sd]
    ax.scatter(x.atac_peak, np.log2(x[MCOL] + 1), s=7, alpha=0.35, color=col, linewidths=0, label=f"{sd} (n={x[MCOL].notna().sum()})")
r = t1[(t1.trackset == "CM") & (t1.readout == "atac_peak")].iloc[0]
ax.set_xscale("log")
ax.set_xlabel("Predicted ATAC over peak (mouse heart track, mean)")
ax.set_ylabel(f"Measured snATAC, log2(CPM/kb + 1)")
ax.set_title(f"Predicted vs measured accessibility: Spearman rho = {r.rho:.2f} (n = {int(r.n)})", fontsize=10, loc="left")
ax.legend(frameon=False, loc="lower right")
save(fig, "fig1_native_vs_measured.png")

# ================= 2. shuffle effects by Stable / Dynamic
print("\n=== 2. shuffle effect (median log2FC) on every ATAC readout, Stable vs Dynamic, with 95% CI")
rows = []
EFF = {}
for ts in sorted(d.trackset.unique()):
    for rd in READ.values():
        EFF[(ts, rd)] = effects(d, ts, rd).join(sites[["stable_dynamic", "age_class", "ex_class", "binding_category"]])
        for sd in ("Stable", "Dynamic"):
            e = EFF[(ts, rd)]
            e = e[e.stable_dynamic == sd]
            for w in WIDTHS:
                lo, hi = boot_ci(e[w])
                rows.append({"trackset": ts, "readout": rd, "class": sd, "shuffle_width": w, "median_log2fc": e[w].median(), "ci_lo": lo, "ci_hi": hi, "n": int(e[w].notna().sum())})
t2 = pd.DataFrame(rows)
t2.to_csv(OUT / "tables" / "shuffle_effect_stable_vs_dynamic.csv", index=False)
cm2 = t2[t2.trackset == "CM"]
print(cm2.pivot_table(index=["class", "shuffle_width"], columns="readout", values="median_log2fc").loc[:, ["atac_motif", "atac_170", "atac_peak"]].round(2).to_string())

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), sharey=True)
vmax = np.abs(cm2.median_log2fc).max()
for ax, sd in zip(axes, ("Stable", "Dynamic")):
    m = cm2[cm2["class"] == sd].pivot(index="shuffle_width", columns="readout", values="median_log2fc").loc[WIDTHS, list(READ.values())]
    im = ax.imshow(m.values, cmap=DIV, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(3), ["Motif", "170 bp", "Peak"])
    ax.set_yticks(range(3), ["Motif shuffle", "170 bp shuffle", "Peak shuffle"])
    ax.set_xlabel("ATAC readout region")
    n = int((cm2[cm2["class"] == sd].n).iloc[0])
    ax.set_title(f"{sd} (n = {n})", fontsize=10, loc="left")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", color=INK if abs(m.values[i, j]) < 0.6 * vmax else "white", fontsize=10)
fig.colorbar(im, ax=axes, shrink=0.85, label="Median log2FC (shuffled / native)")
save(fig, "fig2_shuffle_effect_heatmaps.png")

# ================= 3. classes (age / exercise / binding pattern), matched readout
print("\n=== 3. matched-readout shuffle effect by class vs Stable (Mann-Whitney, BH within each family)")
rows = []
M = pd.DataFrame({w: EFF[("CM", READ[w])][w] for w in WIDTHS}).join(sites[["stable_dynamic", "age_class", "ex_class", "binding_category"]])
for fam, col, ref in (("age", "age_class", "Stable"), ("exercise", "ex_class", "Stable"), ("binding", "binding_category", "YS+MAS+MAE")):
    tests = []
    for w in WIDTHS:
        base = M.loc[M[col] == ref, w].dropna()
        for cl in sorted(M[col].unique()):
            x = M.loc[M[col] == cl, w].dropna()
            lo, hi = boot_ci(x)
            p = mannwhitneyu(x, base).pvalue if cl != ref and len(x) >= 5 else np.nan
            rows.append({"family": fam, "class": cl, "shuffle_width": w, "n": len(x), "median_log2fc": x.median(), "ci_lo": lo, "ci_hi": hi,
                         "diff_vs_ref": x.median() - base.median() if cl != ref else 0.0, "p": p})
t3 = pd.DataFrame(rows)
for fam in t3.family.unique():
    ix = (t3.family == fam) & t3.p.notna()
    t3.loc[ix, "fdr"] = bh(t3.loc[ix, "p"])
t3.to_csv(OUT / "tables" / "shuffle_effect_by_class.csv", index=False)
print(t3.pivot_table(index=["family", "class"], columns="shuffle_width", values=["n", "median_log2fc", "fdr"], aggfunc="first")
      .round(3).to_string())

fig, axes = plt.subplots(3, 3, figsize=(10, 8), sharex="col")
for r, (fam, title) in enumerate((("age", "Age (YS vs MAS)"), ("exercise", "Exercise (MAS vs MAE)"), ("binding", "Binding pattern"))):
    x = t3[t3.family == fam]
    order = list(dict.fromkeys(x["class"]))
    for c, w in enumerate(WIDTHS):
        ax = axes[r, c]
        y = x[x.shuffle_width == w].set_index("class").loc[order]
        pos = np.arange(len(order))[::-1]
        ax.errorbar(y.median_log2fc, pos, xerr=[y.median_log2fc - y.ci_lo, y.ci_hi - y.median_log2fc], fmt="o", color=BLUE, ecolor=BLUE, ms=5, lw=1.2, capsize=0)
        ax.set_yticks(pos, [f"{k} (n={int(n)})" for k, n in zip(order, y.n)], fontsize=8)
        ax.axvline(0, color=MUTED, lw=0.8)
        if r == 0:
            ax.set_title(f"{LABEL[w]} shuffle, {LABEL[w].lower()} readout", fontsize=10, loc="left")
        if r == 2:
            ax.set_xlabel("Median log2FC (95% CI)")
        if c == 0:
            ax.annotate(title, (0, 1.12), xycoords="axes fraction", fontsize=9, color=MUTED)
fig.tight_layout()
save(fig, "fig3_shuffle_effect_by_class.png")

# ================= 4. RNA
print("\n=== 4. RNA response to shuffles (median log2FC), by distance of the nearest TSS")
rows = []
for col in ("rna_prox2kb", "rna_mean50kb"):
    e = effects(d, "CM", col).join(sites[["stable_dynamic", "tss_dist"]])
    e["tss_bin"] = pd.cut(e.tss_dist.abs(), [-1, 2000, 10000, 30000, 1e9], labels=["<=2 kb", "2-10 kb", "10-30 kb", ">30 kb"])
    for (b, sd), g in e.groupby(["tss_bin", "stable_dynamic"], observed=True):
        for w in WIDTHS:
            lo, hi = boot_ci(g[w])
            rows.append({"readout": col, "tss_bin": b, "class": sd, "shuffle_width": w, "median_log2fc": g[w].median(), "ci_lo": lo, "ci_hi": hi, "n": int(g[w].notna().sum())})
t4 = pd.DataFrame(rows)
t4.to_csv(OUT / "tables" / "rna_effect_by_tss_distance.csv", index=False)
print(t4[t4["class"] == "Stable"].pivot_table(index=["readout", "tss_bin"], columns="shuffle_width", values="median_log2fc").round(3).to_string())
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
for ax, col, ttl in zip(axes, ("rna_prox2kb", "rna_mean50kb"), ("RNA within 2 kb of the nearest TSS", "RNA mean, +-50 kb")):
    y = t4[(t4.readout == col) & (t4.shuffle_width == "peak")]
    bins = ["<=2 kb", "2-10 kb", "10-30 kb", ">30 kb"]
    for sd, color, off in (("Stable", BLUE, -0.08), ("Dynamic", ORANGE, 0.08)):
        g = y[y["class"] == sd].set_index("tss_bin").reindex(bins)
        xs = np.arange(len(bins)) + off
        ax.errorbar(xs, g.median_log2fc, yerr=[g.median_log2fc - g.ci_lo, g.ci_hi - g.median_log2fc], fmt="o-", color=color, ms=5, lw=1.2, capsize=0, label=sd)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks(range(4), bins)
    ax.set_xlabel("Distance from the site to the nearest TSS")
    ax.set_title(ttl, fontsize=10, loc="left")
axes[0].set_ylabel("Median log2FC, peak shuffle")
axes[0].legend(frameon=False)
save(fig, "fig4_rna_effect_by_tss_distance.png")

atac_peak = effects(d, "CM", "atac_peak")["peak"]
for col in ("rna_prox2kb", "rna_mean50kb"):
    e = effects(d, "CM", col)["peak"]
    m = atac_peak.notna() & e.notna()
    r, p = spearmanr(atac_peak[m], e[m])
    print(f"ATAC vs {col} peak-shuffle effect across sites: rho = {r:+.2f} (n = {int(m.sum())}, p = {p:.2g})")

# ================= 5. confidence tiers
print("\n=== 5. confidence tiers (chosen from native prediction vs measured, not from the effects)")
n = d[(d.condition == "native") & (d.trackset == "CM")].set_index("site_id").join(meas[[MCOL]])
floor = n.atac_peak.quantile(0.25)
n["has_meas"] = n[MCOL].notna()
n["above_floor"] = n.atac_peak >= floor
pr, mr = n.atac_peak.rank(pct=True), n[MCOL].rank(pct=True)
resid = (pr - mr)
n["agree"] = n.has_meas & (resid.abs() <= resid.std())
tiers = {"A: all sites": n.index, "B: measured peak and predicted >= 25th pct": n.index[n.has_meas & n.above_floor],
         "C: B and rank-agreement within 1 SD": n.index[n.has_meas & n.above_floor & n.agree]}
rows = []
for tname, ix in tiers.items():
    for sd in ("Stable", "Dynamic"):
        g = M.loc[M.index.intersection(ix)]
        g = g[g.stable_dynamic == sd]
        for w in WIDTHS:
            lo, hi = boot_ci(g[w])
            rows.append({"tier": tname, "class": sd, "shuffle_width": w, "n": int(g[w].notna().sum()), "median_log2fc": g[w].median(), "ci_lo": lo, "ci_hi": hi})
t5 = pd.DataFrame(rows)
t5.to_csv(OUT / "tables" / "confidence_tiers.csv", index=False)
print(t5.pivot_table(index=["tier", "class"], columns="shuffle_width", values=["n", "median_log2fc"], aggfunc="first").round(2).to_string())
fig, axes = plt.subplots(1, 3, figsize=(10, 3.3), sharey=True)
for ax, w in zip(axes, WIDTHS):
    for k, (sd, color) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
        y = t5[(t5.shuffle_width == w) & (t5["class"] == sd)]
        xs = np.arange(len(y)) + (k - 0.5) * 0.16
        ax.errorbar(xs, y.median_log2fc, yerr=[y.median_log2fc - y.ci_lo, y.ci_hi - y.median_log2fc], fmt="o", color=color, ms=5, lw=1.2, capsize=0, label=sd)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks(range(3), ["A all", "B floor", "C agree"])
    ax.set_title(f"{LABEL[w]} shuffle, {LABEL[w].lower()} readout", fontsize=10, loc="left")
axes[0].set_ylabel("Median log2FC (95% CI)")
axes[0].legend(frameon=False)
save(fig, "fig5_confidence_tiers.png")

# ================= 6. replicate noise
s = d[(d.trackset == "CM") & d.condition.str.startswith("shuffle")].copy()
nat_ = d[(d.condition == "native") & (d.trackset == "CM")].set_index("site_id").atac_peak
pc = 0.05 * nat_.median()
s["l"] = np.log2((s.atac_peak + pc) / (s.site_id.map(nat_) + pc))
sd_ = s.groupby(["condition", "site_id"]).l.agg(["mean", "std"]).groupby("condition").median()
sd_.to_csv(OUT / "tables" / "replicate_noise_atac_peak.csv")
print("\n=== 6. replicate SD vs effect (CM, atac_peak readout, median over sites)\n", sd_.round(3).to_string())
print(f"\nwritten to {OUT}")
