#!/usr/bin/env python3
"""Analysis of the overnight runs: matched controls + width titration (CM), mutation arm (CM, FB), fibroblast shuffle arm (replication of the CM
Stable-vs-Dynamic result, class effects, adjusted model), neutral-context MPRA vs genomic-context MPRA.
Tables -> Analysis/overnight/tables, figures -> Analysis/overnight/figures; used by 25_build_report.py.   python3 29_overnight_analysis.py

Effect = per-site mean over replicates of log2((edited + pc)/(native + pc)), pc = 5% of the median native value of that readout and trackset.
Dependency (mutation arm) = -effect. Tests: Mann-Whitney (Dynamic vs Stable), paired Wilcoxon (controls), BH-FDR within each figure; bootstrap CIs.
"""
from __future__ import annotations

import glob
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr, t as tdist, wilcoxon

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_labels import EDIT_LABEL, REGION_LABEL, PEAK_MEDIAN, short  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "overnight"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
R = ROOT / "Katana" / "results"
rng = np.random.default_rng(21)

INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, "#ececea", "#e34948"])
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
WIDTHS = ["motif", "170", "peak"]
RD = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
LBL = {"motif": "Motif", "170": "170 bp", "peak": "Peak"}


def stars(q):
    return "" if pd.isna(q) else "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot(x, fn=np.median, n=1000):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 5:
        return np.nan, np.nan
    return tuple(np.percentile(fn(rng.choice(x, (n, len(x))), axis=1), [2.5, 97.5]))


def boot_diff(a, b, n=1000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.median(rng.choice(a, (n, len(a))), 1) - np.median(rng.choice(b, (n, len(b))), 1)
    return tuple(np.percentile(d, [2.5, 97.5]))


def cliffs(a, b):
    return 2 * mannwhitneyu(a, b).statistic / (len(a) * len(b)) - 1


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


def rho_ci(a, b):
    m = pd.notna(a) & pd.notna(b)
    n = int(m.sum())
    r, p = spearmanr(np.asarray(a)[m], np.asarray(b)[m])
    se = 1.06 / np.sqrt(n - 3)
    return r, np.tanh(np.arctanh(r) - 1.96 * se), np.tanh(np.arctanh(r) + 1.96 * se), p, n


def load(p):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{p}_chunk*.csv")))], ignore_index=True)


def lfc_tab(d, ts, col):
    nat = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[(d.condition != "native") & (d.trackset == ts)].copy()
    e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "condition"]).l.mean().unstack(), nat


def labels(d, sites_csv):
    s = pd.read_csv(ROOT / "Katana" / "data" / sites_csv).set_index("site_id")
    l = d.drop_duplicates("site_id").set_index("site_id")[["stable_dynamic", "age_class", "ex_class", "binding_category", "tss_dist"]]
    l["age_class"] = l.age_class.replace({"Unbound": "Neither YS nor MAS (= MAE only)"})
    l["ex_class"] = l.ex_class.replace({"Unbound": "Neither MAS nor MAE (= YS only)"})
    return l.join(s[["peak_width"]])


def save(fig, name):
    fig.savefig(OUT / "figures" / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def star_row(ax, xs, labs, y=1.0):
    for x, s in zip(xs, labs):
        ax.text(x, y, s, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, color=INK)


# ============================================================ D1. controls + titration (CM sites, heart ATAC)
co = load("controls_cm")
lc = labels(co, "cm_fos_sites.csv")
SD = lc.stable_dynamic
offs = ["offset@-600", "offset@-300", "offset@+300", "offset@+600"]
tit = [f"titrate@{w}" for w in (7, 11, 21, 41, 85, 170, 340, 680)]
tw = [7, 11, 21, 41, 85, 170, 340, 680]
rows, rows_c = [], []
fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
tm, _ = lfc_tab(co, "CM", "atac_motif")
order = offs[:2] + ["titrate@170"] + offs[2:]
names = ["-600", "-300", "centred", "+300", "+600"]
for k, (sd, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
    ids = SD[SD == sd].index
    med, lo, hi = [], [], []
    for c in order:
        v = tm.loc[ids, c].dropna()
        l_, h_ = boot(v)
        med.append(v.median()); lo.append(l_); hi.append(h_)
    xs = np.arange(5) + (k - 0.5) * 0.14
    ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt="o", color=col, ms=5, lw=1.2, capsize=0, label=sd)
    cen, off = tm.loc[ids, "titrate@170"], tm.loc[ids, ["offset@-300", "offset@+300"]].mean(1)
    pw = wilcoxon(cen, off).pvalue
    rows_c.append({"class": sd, "n": len(ids), "centred_170_median": cen.median(), "offset_pm300_median": off.median(), "share_centred_more_negative": float((cen < off).mean()), "paired_wilcoxon_p": pw})
ax.set_xticks(range(5), names)
ax.set_xlabel("Centre of the shuffled 170 bp region, relative to the motif (bp)")
ax.set_ylabel("Median effect on the motif (log2FC)")
ax.axhline(0, color=MUTED, lw=0.8)
ax.legend(frameon=False)
ax.set_title("Same-size shuffles: centred vs elsewhere in the peak", fontsize=10, loc="left")
ax = axes[1]
tab = []
for rd, ls_ in (("atac_motif", "-"), ):
    pass
res_w = {}
for rd_i, (w_rd, mk) in enumerate((("motif", "o"), ("170", "s"), ("peak", "^"))):
    t_, _ = lfc_tab(co, "CM", RD[w_rd])
    for sd, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
        ids = SD[SD == sd].index
        med, lo, hi = [], [], []
        for c in tit:
            v = t_.loc[ids, c].dropna()
            l_, h_ = boot(v, n=500)
            med.append(v.median()); lo.append(l_); hi.append(h_)
            rows.append({"readout": w_rd, "class": sd, "shuffle_width_bp": int(c.split("@")[1]), "median": v.median(), "ci_lo": l_, "ci_hi": h_, "n": len(v)})
        if w_rd == "motif":
            ax.plot(tw, med, "-o", color=col, ms=4.5, lw=1.6, label=f"{sd}, motif readout")
            ax.fill_between(tw, lo, hi, color=col, alpha=0.15, lw=0)
        else:
            ax.plot(tw, med, "-" + mk, color=col, ms=3.5, lw=1, alpha=0.55, label=f"{sd}, {LBL[w_rd].lower()} readout" if sd == "Stable" or True else None)
res = pd.DataFrame(rows)
# stable vs dynamic per width (motif readout) with stars
sv = []
for c in tit:
    a_, b_ = tm.loc[SD[SD == "Dynamic"].index, c].dropna(), tm.loc[SD[SD == "Stable"].index, c].dropna()
    sv.append({"width": int(c.split("@")[1]), "diff_dyn_minus_stable": a_.median() - b_.median(), "p": mannwhitneyu(a_, b_).pvalue})
sv = pd.DataFrame(sv)
sv["fdr"] = bh(sv.p)
sv["sig"] = sv.fdr.map(stars)
star_row(ax, tw, sv.sig.tolist())
ax.set_xscale("log")
ax.set_xticks(tw, [str(w) for w in tw])
ax.set_xlabel("Width of the centred shuffle (bp, log scale)")
ax.set_ylabel("Median effect (log2FC)")
ax.set_title("Dose-response: effect grows with the shuffled width", fontsize=10, loc="left", pad=18)
leg = [Line2D([], [], color=BLUE, lw=2.2, label="Stable sites"), Line2D([], [], color=ORANGE, lw=2.2, label="Dynamic sites"),
       Line2D([], [], color=INK, lw=2, marker="o", ms=5, label="effect read over the 11 bp motif (thick line, 95% CI band)"),
       Line2D([], [], color=INK, lw=1, alpha=0.6, marker="s", ms=4, label="effect read over the 170 bp (thin line)"),
       Line2D([], [], color=INK, lw=1, alpha=0.6, marker="^", ms=4, label="effect read over the whole peak (thin line)")]
ax.legend(handles=leg, frameon=False, fontsize=7.5, loc="lower left")
fig.text(0.50, -0.03, "Colour = site class. Line weight and marker = which region the effect is read over. Marks above the plot: Stable vs Dynamic, motif readout (Mann-Whitney, BH).", fontsize=8, color=MUTED)
fig.tight_layout()
save(fig, "d1_controls_titration")

# same plot as the left panel, but with the effect read over the 11 bp motif, the 170 bp and the whole peak
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=False)
rows_o = []
for ax, (w_rd, ttl) in zip(axes, (("motif", "Effect read over the 11 bp motif"), ("170", "Effect read over the 170 bp"), ("peak", "Effect read over the whole peak"))):
    t_, _ = lfc_tab(co, "CM", RD[w_rd])
    for k, (sd, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
        ids = SD[SD == sd].index
        med, lo, hi = [], [], []
        for c in order:
            v = t_.loc[ids, c].dropna()
            l_, h_ = boot(v)
            med.append(v.median()); lo.append(l_); hi.append(h_)
            rows_o.append({"readout": w_rd, "class": sd, "shuffle_position": c, "median": v.median(), "ci_lo": l_, "ci_hi": h_, "n": len(v)})
        xs = np.arange(5) + (k - 0.5) * 0.14
        ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt="o", color=col, ms=5, lw=1.2, capsize=0, label=sd)
        cen, off = t_.loc[ids, "titrate@170"], t_.loc[ids, ["offset@-300", "offset@+300"]].mean(1)
        rows_o.append({"readout": w_rd, "class": sd, "shuffle_position": "centred vs +-300 (paired)", "median": cen.median() - off.median(), "ci_lo": np.nan, "ci_hi": np.nan,
                       "n": int((cen < off).sum()), "share_centred_more_negative": float((cen < off).mean()), "paired_wilcoxon_p": wilcoxon(cen, off).pvalue})
    ax.set_xticks(range(5), names)
    ax.set_xlabel("Centre of the shuffled 170 bp, relative to the motif (bp)")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_title(ttl, fontsize=10, loc="left")
axes[0].set_ylabel("Median effect (log2FC)")
axes[0].legend(frameon=False)
fig.text(0.5, -0.03, "Every panel shows the same experiment (a 170 bp shuffle centred on the motif or moved 300 or 600 bp away); only the region the effect is read over differs. Median with 95% bootstrap CI.", fontsize=8, color=MUTED, ha="center")
fig.tight_layout()
save(fig, "d1c_offsets_by_readout")
pd.DataFrame(rows_o).to_csv(OUT / "tables" / "d1c_offsets_by_readout.csv", index=False)
res.to_csv(OUT / "tables" / "d1_titration_medians.csv", index=False)
pd.DataFrame(rows_c).to_csv(OUT / "tables" / "d1_centred_vs_offset.csv", index=False)
sv.to_csv(OUT / "tables" / "d1_titration_stable_vs_dynamic.csv", index=False)
# fraction of the 170-bp loss reached by each width (Stable, motif readout)
st = res[(res.readout == "motif") & (res["class"] == "Stable")].set_index("shuffle_width_bp")["median"]
frac = (st / st.loc[170]).rename("fraction_of_170bp_effect")
frac.to_csv(OUT / "tables" / "d1_fraction_of_170.csv")

# ============================================================ D2. mutation arm
mc, mf = load("mut_cm"), load("mut_fb")
EDITS = [("random@motif", "motif"), ("dinuc@motif", "motif"), ("core@motif", "motif"), ("flank@motif", "motif"),
         ("random@170", "170"), ("dinuc@170", "170"), ("flank@170", "170"), ("random@peak", "peak"), ("dinuc@peak", "peak")]
rows = []
fig, axes = plt.subplots(2, 1, figsize=(13, 7.6), sharex=True)
for ax, (d, ts, sites_csv, ttl) in zip(axes, ((mc, "CM", "cm_fos_sites.csv", "Cardiomyocyte sites, heart ATAC"), (mf, "FB", "fb_fos_sites.csv", "Fibroblast sites, Patski ATAC"))):
    l_ = labels(d, sites_csv)
    sdv = l_.stable_dynamic
    cache = {w: -lfc_tab(d, ts, RD[w])[0] for w in WIDTHS}
    ps = []
    for i, (cname, w) in enumerate(EDITS):
        s_ = cache[w][cname].dropna()
        for k, (sd, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
            v = cache[w].loc[sdv[sdv == sd].index, cname].dropna()
            l, h = boot(v, np.mean, 500)
            ax.bar(i + (k - 0.5) * 0.38, v.mean(), 0.34, color=col, alpha=0.85)
            ax.plot([i + (k - 0.5) * 0.38] * 2, [l, h], color=INK, lw=1)
            rows.append({"sites": ttl, "edit": cname, "readout": w, "class": sd, "mean_dependency": v.mean(), "ci_lo": l, "ci_hi": h, "n": len(v)})
        a_, b_ = cache[w].loc[sdv[sdv == "Dynamic"].index, cname].dropna(), cache[w].loc[sdv[sdv == "Stable"].index, cname].dropna()
        ps.append(mannwhitneyu(a_, b_).pvalue)
    q = bh(ps)
    star_row(ax, range(len(EDITS)), [stars(x) for x in q])
    for r_ in rows[-2 * len(EDITS):]:
        pass
    for i, (cname, w) in enumerate(EDITS):
        rows.append({"sites": ttl, "edit": cname, "readout": w, "class": "Stable-vs-Dynamic FDR", "mean_dependency": q[i], "ci_lo": np.nan, "ci_hi": np.nan, "n": 0})
    ax.set_title(ttl, fontsize=10, loc="left", pad=18)
    ax.set_ylabel("Dependency (-log2FC)\nmean over sites")
axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BLUE, alpha=0.85), plt.Rectangle((0, 0), 1, 1, color=ORANGE, alpha=0.85)],
               labels=["Stable sites", "Dynamic sites"], frameon=False, ncol=2, loc="upper left")
axes[1].set_xticks(range(len(EDITS)), [short(e) for e, _ in EDITS], fontsize=8)
for ax in axes:
    ax.axvline(3.5, color=MUTED, lw=0.8, ls=":"); ax.axvline(6.5, color=MUTED, lw=0.8, ls=":")
for x_, key_ in ((1.5, "motif"), (5.0, "170"), (7.5, "peak")):
    axes[1].text(x_, -0.42, REGION_LABEL[key_], transform=axes[1].get_xaxis_transform(), ha="center", va="top", fontsize=9, fontweight="bold")
fig.text(0.01, -0.01, "Each edit is read out over the region it changes. Bar labels give the edit and the exact number of bases changed. Bars: mean with 95% bootstrap CI. Marks: Stable vs Dynamic, Mann-Whitney, BH within panel.", fontsize=8, color=MUTED)
fig.tight_layout()
save(fig, "d2_mutation_arm")
D2 = pd.DataFrame(rows)
D2.to_csv(OUT / "tables" / "d2_mutation_dependency.csv", index=False)

# decomposition: what carries the dependence
rows = []
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6))
for ax, (w, items, ttl) in zip(axes, (("motif", [("random@motif", "Replace all\n11 bp of the motif"), ("core@motif", "Replace the\n7 bp core only"), ("flank@motif", "Replace the 4 flank\nbases (core kept)")], "Within the 11 bp motif (motif readout)"),
                                      ("170", [("random@170", "Replace all 170 bp\n(motif included)"), ("flank@170", "Replace the 159 bp around\nthe motif (motif kept)")], "Within the 170 bp (170 bp readout)"))):
    for sites_label, d, ts, off, col in (("CM, heart ATAC", mc, "CM", -0.17, BLUE), ("FB, Patski ATAC", mf, "FB", 0.17, AQUA)):
        l_ = labels(d, "cm_fos_sites.csv" if ts == "CM" else "fb_fos_sites.csv")
        stab = l_.index[l_.stable_dynamic == "Stable"]
        t_ = -lfc_tab(d, ts, RD[w])[0].loc[stab]
        vals = []
        for i, (c, _) in enumerate(items):
            v = t_[c].dropna()
            lo, hi = boot(v, np.mean, 500)
            ax.bar(i + off, v.mean(), 0.32, color=col, alpha=0.85, label=sites_label if i == 0 else None)
            ax.plot([i + off] * 2, [lo, hi], color=INK, lw=1)
            vals.append(v.mean())
            rows.append({"panel": w, "edit": c, "sites": sites_label, "mean_dependency_stable": v.mean(), "ci_lo": lo, "ci_hi": hi})
        if w == "170":
            ax.text(1 + off, vals[1] + 0.05, f"{vals[1] / vals[0]:.0%} of\nthe full loss", ha="center", fontsize=8)
        else:
            ax.text(1 + off, vals[1] + 0.03, f"{vals[1] / vals[0]:.0%}", ha="center", fontsize=8)
    ax.set_xticks(range(len(items)), [lab for _, lab in items], fontsize=8.5)
    ax.set_title(ttl, fontsize=10, loc="left")
    ax.set_ylabel("Dependency (-log2FC), Stable sites")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "d2b_what_carries_dependence")
pd.DataFrame(rows).to_csv(OUT / "tables" / "d2b_decomposition.csv", index=False)

# ============================================================ D3. fibroblast sites: shuffle arm, replication, classes
fb = load("shuffle_fb")
cmsh = load("shuffle_cm")
lfb = labels(fb, "fb_fos_sites.csv")
lcm = labels(cmsh, "cm_fos_sites.csv")


def matched_eff(d, ts, w):
    return lfc_tab(d, ts, RD[w])[0][f"shuffle@{w}"]


# heatmap for FB sites / Patski
rows = []
mats = {}
for sd in ("Stable", "Dynamic"):
    mats[sd] = np.zeros((3, 3))
pmat = np.zeros((3, 3))
cache = {rd: lfc_tab(fb, "FB", RD[rd])[0] for rd in WIDTHS}
sdf = lfb.stable_dynamic
for i, sh in enumerate(WIDTHS):
    for j, rd in enumerate(WIDTHS):
        e = cache[rd][f"shuffle@{sh}"]
        s_, d_ = e[sdf[sdf == "Stable"].index].dropna(), e[sdf[sdf == "Dynamic"].index].dropna()
        mats["Stable"][i, j], mats["Dynamic"][i, j] = s_.median(), d_.median()
        pmat[i, j] = mannwhitneyu(d_, s_).pvalue
        rows.append({"shuffle": sh, "readout": rd, "median_stable": s_.median(), "median_dynamic": d_.median(), "diff": d_.median() - s_.median(), "p": pmat[i, j]})
q = bh(pmat.ravel()).reshape(3, 3)
pd.DataFrame(rows).assign(fdr=q.ravel()).to_csv(OUT / "tables" / "d3_fb_heatmap_cells.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.3), sharey=True)
vmax = 2.2
for ax, sd in zip(axes, ("Stable", "Dynamic")):
    im = ax.imshow(mats[sd], cmap=DIV, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(3), ["Motif", "170 bp", "Peak"]); ax.set_yticks(range(3), ["Motif shuffle", "170 bp shuffle", "Peak shuffle"])
    ax.set_xlabel("ATAC readout region")
    ax.set_title(f"{sd} (n = {int((sdf == sd).sum()):,})", fontsize=10, loc="left")
    for i in range(3):
        for j in range(3):
            txt = f"{mats[sd][i, j]:.2f}" + (f"\n{stars(q[i, j])}" if sd == "Dynamic" else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=9.5, color=INK if abs(mats[sd][i, j]) < 0.6 * vmax else "white")
fig.colorbar(im, ax=axes, shrink=0.85, label="Median log2FC (shuffled / native)")
fig.text(0.01, -0.04, "Fibroblast sites, Patski ATAC track. Marks on the Dynamic panel: Mann-Whitney vs Stable, BH over the 9 cells.", fontsize=8, color=MUTED)
fig.savefig(OUT / "figures" / "d3a_fb_heatmaps.png", bbox_inches="tight")
plt.close(fig)

# replication forest: Dynamic - Stable (matched readout) for 2 site sets x 2 tracksets
rows = []
combos = (("CM sites, heart ATAC", cmsh, lcm, "CM"), ("CM sites, Patski ATAC", cmsh, lcm, "FB"), ("FB sites, heart ATAC", fb, lfb, "CM"), ("FB sites, Patski ATAC", fb, lfb, "FB"))
for name, d, l_, ts in combos:
    sd_ = l_.stable_dynamic
    for w in WIDTHS:
        e = matched_eff(d, ts, w)
        s_, d_ = e[sd_[sd_ == "Stable"].index].dropna(), e[sd_[sd_ == "Dynamic"].index].dropna()
        lo, hi = boot_diff(d_, s_)
        rows.append({"combo": name, "width": w, "diff_dyn_minus_stable": d_.median() - s_.median(), "ci_lo": lo, "ci_hi": hi, "cliffs_delta": cliffs(d_, s_), "p": mannwhitneyu(d_, s_).pvalue})
F = pd.DataFrame(rows)
F["fdr"] = bh(F.p)
F["sig"] = F.fdr.map(stars)
F.to_csv(OUT / "tables" / "d3_replication_stable_vs_dynamic.csv", index=False)
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2), sharey=True)
names = [c[0] for c in combos]
for ax, w in zip(axes, WIDTHS):
    for i, nm in enumerate(names):
        r_ = F[(F.combo == nm) & (F.width == w)].iloc[0]
        y = len(names) - 1 - i
        ax.errorbar(r_.diff_dyn_minus_stable, y, xerr=[[r_.diff_dyn_minus_stable - r_.ci_lo], [r_.ci_hi - r_.diff_dyn_minus_stable]], fmt="o", color=BLUE if "CM sites" in nm else AQUA, ms=5, lw=1.2, capsize=0)
        ax.text(1.02, y, r_.sig, transform=ax.get_yaxis_transform(), va="center", fontsize=8)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_title(f"{LBL[w]} shuffle", fontsize=10, loc="left")
    ax.set_xlabel("Dynamic - Stable (median log2FC, 95% CI)")
axes[0].set_yticks(range(len(names)), names[::-1], fontsize=8.5)
fig.text(0.01, -0.03, "Positive = Dynamic sites lose less accessibility than Stable sites. Marks: Mann-Whitney, BH over the 12 comparisons. Blue: CM sites; green: FB sites.", fontsize=8, color=MUTED)
fig.tight_layout()
save(fig, "d3b_replication_forest")

# classes in FB sites (Patski track)
FAMS = (("age", "age_class", "Stable"), ("exercise", "ex_class", "Stable"), ("binding", "binding_category", "YS+MAS+MAE"))
rows = []
for r_i, (fam, col, ref) in enumerate(FAMS):
    classes = list(dict.fromkeys(lfb[col].dropna()))
    classes = [ref] + [c for c in classes if c != ref]
    ps_all = []
    meds = {}
    for c_i, w in enumerate(WIDTHS):
        e = matched_eff(fb, "FB", w)
        base = e[lfb[lfb[col] == ref].index].dropna()
        for cl in classes:
            v = e[lfb[lfb[col] == cl].index].dropna()
            lo, hi = boot(v, n=500)
            p = mannwhitneyu(v, base).pvalue if cl != ref else np.nan
            rows.append({"family": fam, "width": w, "class": cl, "n": len(v), "median": v.median(), "ci_lo": lo, "ci_hi": hi, "raw_diff": v.median() - base.median(), "p": p})
Cl = pd.DataFrame(rows)
for fam in Cl.family.unique():
    ix = (Cl.family == fam) & Cl.p.notna()
    Cl.loc[ix, "fdr"] = bh(Cl.loc[ix, "p"])
Cl["sig"] = Cl.fdr.map(stars)
Cl.to_csv(OUT / "tables" / "d3_fb_class_effects.csv", index=False)
CLASS_LAB = {("age", "Stable"): "Bound in both YS and MAS", ("age", "Age Up"): "Age Up (gained with age)", ("age", "Age Down"): "Age Down (lost with age)",
             ("age", "Neither YS nor MAS (= MAE only)"): "Bound in MAE only", ("exercise", "Stable"): "Bound in both MAS and MAE", ("exercise", "Exercise Up"): "Exercise Up (gained)",
             ("exercise", "Exercise Down"): "Exercise Down (lost)", ("exercise", "Neither MAS nor MAE (= YS only)"): "Bound in YS only", ("binding", "YS+MAS+MAE"): "YS+MAS+MAE (Stable)"}
FAM_TITLE = {"age": "Age: young vs aged sedentary", "exercise": "Exercise: aged sedentary vs aged exercised", "binding": "All seven binding patterns"}
for w in WIDTHS:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, (fam, col, ref) in zip(axes, FAMS):
        x = Cl[Cl.family == fam]
        order = list(dict.fromkeys(x["class"]))
        y = x[x.width == w].set_index("class").loc[order]
        pos = np.arange(len(order))[::-1]
        for pi, k in zip(pos, order):
            v = y.loc[k]
            ax.errorbar(v["median"], pi, xerr=[[v["median"] - v.ci_lo], [v.ci_hi - v["median"]]], fmt="o", color=BLUE if k == ref else ORANGE, ms=6, lw=1.4, capsize=0)
            if k != ref:
                ax.text(1.02, pi, v.sig, transform=ax.get_yaxis_transform(), va="center", fontsize=9)
        ax.axvline(y.loc[ref, "median"], color=BLUE, lw=0.8, ls=":")
        ax.set_yticks(pos, [f"{CLASS_LAB.get((fam, k), k)}  (n={int(n)})" for k, n in zip(order, y.n)], fontsize=9)
        ax.set_title(FAM_TITLE[fam], fontsize=10, loc="left")
        ax.set_xlabel("Median effect (log2FC), 95% CI")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    fig.text(0.01, -0.03, f"Fibroblast sites, Patski ATAC. {LBL[w]} shuffle, effect read over the {LBL[w].lower()}. Blue = reference class. Marks: Mann-Whitney against the reference, BH-corrected within each panel.", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    save(fig, f"d3c_fb_class_effects_{w}")

# adjusted Stable-vs-Dynamic (native ATAC + peak width) for both site sets / tracks
rows = []
for name, d, l_, ts in combos:
    for w in WIDTHS:
        e, nat = lfc_tab(d, ts, RD[w])
        y = pd.DataFrame({"y": e[f"shuffle@{w}"]}).join(l_[["stable_dynamic", "peak_width"]])
        y["log_nat"] = np.log10(nat.reindex(y.index) + 1e-3)
        y["log_pw"] = np.log10(y.peak_width)
        y = y.dropna()
        X = np.column_stack([np.ones(len(y)), (y.stable_dynamic == "Dynamic").astype(float), y.log_nat - y.log_nat.mean(), y.log_pw - y.log_pw.mean()])
        beta, se, pv = ols_hc3(X, y.y)
        rows.append({"combo": name, "width": w, "raw_diff": F[(F.combo == name) & (F.width == w)].diff_dyn_minus_stable.iloc[0], "adj_diff": beta[1], "adj_lo": beta[1] - 1.96 * se[1], "adj_hi": beta[1] + 1.96 * se[1], "p": pv[1]})
A = pd.DataFrame(rows)
A["fdr"] = bh(A.p)
A["sig"] = A.fdr.map(stars)
A.to_csv(OUT / "tables" / "d3_adjusted_stable_vs_dynamic.csv", index=False)

# native vs measured for FB sites
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
meas = meas[meas.site_table == "fb"].set_index("site_id")
rows = []
for ts, nm in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
    n = fb[(fb.condition == "native") & (fb.trackset == ts)].set_index("site_id").join(meas[["fb_cpm_per_kb"]])
    for rd in ("atac_motif", "atac_170", "atac_peak"):
        r, lo, hi, p, nn = rho_ci(n[rd], n.fb_cpm_per_kb)
        rows.append({"track": nm, "readout": rd, "rho": r, "ci_lo": lo, "ci_hi": hi, "n": nn})
pd.DataFrame(rows).to_csv(OUT / "tables" / "d3_fb_native_vs_measured.csv", index=False)

# ============================================================ D4. neutral-context MPRA vs genomic-context MPRA
mt = pd.read_csv(ROOT / "Katana" / "data" / "MPRA_Results_with_TSS_distance.csv").drop_duplicates("Coords").set_index("Coords")
mg = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / "mpra_chunk*.csv")))], ignore_index=True).drop_duplicates(["mpra_id", "trackset", "condition", "rep"])
mn = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / "mpra_neutral_chunk*.csv")))], ignore_index=True)


def mpra_stats(m, ts, neutral):
    g = m[m.trackset == ts]
    if neutral:
        nat = g[g.condition == "native"].groupby("mpra_id").atac_170.mean()
    else:
        nat = g[g.condition == "native"].groupby("mpra_id").atac_170.mean()
    pc = 0.05 * nat.median()
    out = {"native": nat}
    for c in ("motif_shuffle", "full_shuffle"):
        gg = g[g.condition == c]
        s_ = (gg.groupby(["mpra_id", "background"]).atac_170.mean().groupby("mpra_id").mean() if neutral else gg.groupby("mpra_id").atac_170.mean())
        out[c] = np.log2((s_ + pc) / (nat.reindex(s_.index) + pc))
    return out


rows = []
for ts, tn in (("FB", "Patski ATAC"), ("CM", "heart ATAC")):
    for ctx, m, neu in (("genomic context", mg, False), ("neutral context", mn, True)):
        o = mpra_stats(m, ts, neu)
        for nm, key, y in (("native level vs MPRA activity", "native", mt.log2_alpha_Endo), ("motif-shuffle effect vs MPRA change", "motif_shuffle", mt.delta_vs_MotifShuff), ("full-shuffle effect vs MPRA change", "full_shuffle", mt.delta_vs_FullShuff)):
            ids = [i for i in o[key].index if i in y.index]
            r, lo, hi, p, n = rho_ci(o[key].reindex(ids), y.reindex(ids))
            rows.append({"track": tn, "context": ctx, "comparison": nm, "rho": r, "ci_lo": lo, "ci_hi": hi, "p": p, "n": n})
M = pd.DataFrame(rows)
M.to_csv(OUT / "tables" / "d4_mpra_context_comparison.csv", index=False)
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3), sharey=True)
comps = ["native level vs MPRA activity", "motif-shuffle effect vs MPRA change", "full-shuffle effect vs MPRA change"]
rowsy = [("Patski ATAC", "genomic context"), ("Patski ATAC", "neutral context"), ("heart ATAC", "genomic context"), ("heart ATAC", "neutral context")]
for ax, cp in zip(axes, comps):
    for i, (tn, cx) in enumerate(rowsy):
        r_ = M[(M.track == tn) & (M.context == cx) & (M.comparison == cp)].iloc[0]
        y = len(rowsy) - 1 - i
        ax.errorbar(r_.rho, y, xerr=[[r_.rho - r_.ci_lo], [r_.ci_hi - r_.rho]], fmt="o" if cx == "genomic context" else "s", color=BLUE if cx == "genomic context" else AQUA, ms=5, lw=1.2, capsize=0)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_title(cp, fontsize=9, loc="left")
    ax.set_xlabel("Spearman rho (95% CI)")
axes[0].set_yticks(range(len(rowsy)), [f"{a}, {b}" for a, b in rowsy[::-1]], fontsize=8.5)
fig.text(0.01, -0.03, "For the shuffle panels the MPRA change is native minus scrambled, so agreement would appear as a NEGATIVE rho. n = 1,463 to 1,653 elements.", fontsize=8, color=MUTED)
fig.tight_layout()
save(fig, "d4_mpra_context")

pd.set_option("display.width", 220, "display.max_columns", 30)
print("D1 centred vs offset:\n", pd.DataFrame(rows_c).round(3).to_string(index=False))
print("\nD1 fraction of the 170-bp effect (Stable, motif readout):\n", frac.round(2).to_string())
print("\nD1 stable vs dynamic by width (motif readout):\n", sv.round(3).to_string(index=False))
print("\nD3 replication:\n", F.round(3).to_string(index=False))
print("\nD3 adjusted:\n", A.round(3).to_string(index=False))
print("\nD3 FB native vs measured:\n", pd.read_csv(OUT / "tables" / "d3_fb_native_vs_measured.csv").round(3).to_string(index=False))
print("\nD4:\n", M.round(3).to_string(index=False))
