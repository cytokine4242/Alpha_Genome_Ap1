#!/usr/bin/env python3
"""Three additions to the report.
 1  Is prediction quality different between site classes? Spearman rho (predicted vs measured accessibility) per binding class with bootstrap CIs, pairwise
    differences against Stable, and tests on the absolute and signed rank residual (predicted percentile - measured percentile).
 2  The shuffle results as a grouped dot plot (Stable and Dynamic side by side), an alternative to the heatmap.
 3  Peak swap: Stable vs Dynamic donors (donor accessibility, swapped accessibility, retention, and a model adjusting for donor accessibility).
Outputs Analysis/extra/{tables,figures}.   python3 32_extra_figures_and_class_stats.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import kruskal, mannwhitneyu, spearmanr, chi2

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "extra"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
R = ROOT / "Katana" / "results"
rng = np.random.default_rng(32)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
RD = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}


def stars(q):
    return "" if pd.isna(q) else "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def load(p):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{p}_chunk*.csv")))], ignore_index=True)


def boot_med(x, n=500):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return tuple(np.percentile(np.median(rng.choice(x, (n, len(x))), 1), [2.5, 97.5]))


# =============================================================== 1. prediction quality by class
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
shuf = {"cm": load("shuffle_cm"), "fb": load("shuffle_fb")}
sites = {"cm": pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id"), "fb": pd.read_csv(ROOT / "Data/derived/fb_fos_sites.csv").set_index("site_id")}
COMBOS = (("CM sites, heart ATAC", "cm", "CM", "cm_cpm_per_kb"), ("FB sites, heart ATAC", "fb", "CM", "fb_cpm_per_kb"), ("FB sites, Patski ATAC", "fb", "FB", "fb_cpm_per_kb"))


def rho_boot(a, b, n=1000):
    a, b = np.asarray(a), np.asarray(b)
    idx = rng.integers(0, len(a), (n, len(a)))
    return np.array([spearmanr(a[i], b[i])[0] for i in idx])


cls_rows, pair_rows, test_rows, resid_store = [], [], [], {}
for name, sk, ts, mcol in COMBOS:
    d = shuf[sk]
    nat = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id").atac_peak
    m = meas[meas.site_table == sk].set_index("site_id")[mcol]
    x = pd.DataFrame({"pred": nat, "meas": m}).join(sites[sk][["binding_category", "stable_dynamic"]]).dropna()
    x["resid"] = x.pred.rank(pct=True) - x.meas.rank(pct=True)
    x["absres"] = x.resid.abs()
    x["cls"] = np.where(x.stable_dynamic == "Stable", "Stable (all three)", x.binding_category)
    resid_store[name] = x
    groups = ["Stable (all three)"] + sorted(c for c in x.cls.unique() if c != "Stable (all three)")
    boots = {}
    for g in groups + ["All Dynamic"]:
        s = x[x.stable_dynamic == "Dynamic"] if g == "All Dynamic" else x[x.cls == g]
        b = rho_boot(s.pred, s.meas)
        boots[g] = b
        r = spearmanr(s.pred, s.meas)[0]
        cls_rows.append({"combo": name, "class": g, "n": len(s), "rho": r, "ci_lo": np.percentile(b, 2.5), "ci_hi": np.percentile(b, 97.5),
                         "median_abs_rank_residual": s.absres.median(), "median_signed_residual": s.resid.median()})
    # pairwise vs Stable: difference in rho (bootstrap), abs residual and signed residual (Mann-Whitney)
    ref = x[x.cls == "Stable (all three)"]
    ps = []
    for g in groups[1:] + ["All Dynamic"]:
        s = x[x.stable_dynamic == "Dynamic"] if g == "All Dynamic" else x[x.cls == g]
        diff = boots[g] - boots["Stable (all three)"]
        p_rho = min(1.0, 2 * min((diff <= 0).mean(), (diff >= 0).mean()) + 1 / len(diff))
        pair_rows.append({"combo": name, "class": g, "n": len(s), "rho_diff_vs_stable": spearmanr(s.pred, s.meas)[0] - spearmanr(ref.pred, ref.meas)[0], "rho_diff_lo": np.percentile(diff, 2.5),
                          "rho_diff_hi": np.percentile(diff, 97.5), "p_rho": p_rho, "p_abs_resid": mannwhitneyu(s.absres, ref.absres).pvalue, "p_signed_resid": mannwhitneyu(s.resid, ref.resid).pvalue,
                          "median_signed_resid": s.resid.median(), "median_signed_resid_stable": ref.resid.median()})
    # omnibus
    kw_abs = kruskal(*[x[x.cls == g].absres for g in groups])
    kw_sgn = kruskal(*[x[x.cls == g].resid for g in groups])
    zs = []
    for g in groups:
        s = x[x.cls == g]
        r = np.clip(spearmanr(s.pred, s.meas)[0], -0.999, 0.999)
        zs.append((np.arctanh(r), len(s) - 3))
    z = np.array([a for a, _ in zs]); w = np.array([b for _, b in zs], float)
    Q = float((w * (z - (w * z).sum() / w.sum()) ** 2).sum())
    test_rows.append({"combo": name, "kw_abs_resid_H": kw_abs.statistic, "kw_abs_resid_p": kw_abs.pvalue, "kw_signed_resid_H": kw_sgn.statistic, "kw_signed_resid_p": kw_sgn.pvalue,
                      "rho_heterogeneity_Q": Q, "rho_heterogeneity_df": len(groups) - 1, "rho_heterogeneity_p": chi2.sf(Q, len(groups) - 1),
                      "stable_vs_dynamic_abs_resid_p": mannwhitneyu(x[x.stable_dynamic == "Stable"].absres, x[x.stable_dynamic == "Dynamic"].absres).pvalue,
                      "stable_vs_dynamic_signed_resid_p": mannwhitneyu(x[x.stable_dynamic == "Stable"].resid, x[x.stable_dynamic == "Dynamic"].resid).pvalue,
                      "stable_median_signed": x[x.stable_dynamic == "Stable"].resid.median(), "dynamic_median_signed": x[x.stable_dynamic == "Dynamic"].resid.median()})
CL, PR, TS = pd.DataFrame(cls_rows), pd.DataFrame(pair_rows), pd.DataFrame(test_rows)
for c in ("rho", "abs_resid", "signed_resid"):
    PR[f"fdr_{c}"] = np.nan
    for nm_ in PR.combo.unique():
        mk = (PR.combo == nm_) & (PR["class"] != "All Dynamic")
        PR.loc[mk, f"fdr_{c}"] = bh(PR.loc[mk, f"p_{c}"].values)
CL.to_csv(OUT / "tables" / "prediction_quality_by_class.csv", index=False)
PR.to_csv(OUT / "tables" / "prediction_quality_vs_stable.csv", index=False)
TS.to_csv(OUT / "tables" / "prediction_quality_omnibus.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.2, 1, 1]})
name = "CM sites, heart ATAC"
x = resid_store[name]
g_order = ["Stable (all three)"] + [c for c in CL[(CL.combo == name) & ~CL["class"].isin(["Stable (all three)", "All Dynamic"])].sort_values("rho", ascending=False)["class"]]
c_ = CL[CL.combo == name].set_index("class")
p_ = PR[(PR.combo == name)].set_index("class")
ax = axes[0]
for i, g in enumerate(g_order[::-1]):
    r = c_.loc[g]
    ax.errorbar(r.rho, i, xerr=[[r.rho - r.ci_lo], [r.ci_hi - r.rho]], fmt="o", color=BLUE if g.startswith("Stable") else ORANGE, ms=6, lw=1.4, capsize=0)
    ax.text(1.02, i, f"n={int(r.n)}" + ("" if g.startswith("Stable") else f"  {stars(p_.loc[g, 'fdr_rho'])}"), transform=ax.get_yaxis_transform(), va="center", fontsize=8)
ax.set_yticks(range(len(g_order)), g_order[::-1], fontsize=8.5)
ax.set_xlabel("Spearman rho, predicted vs measured (95% CI)")
ax.set_title("Prediction accuracy by binding class\n(marks: difference from Stable, BH)", fontsize=9.5, loc="left")
for ax, col, ttl, yl in ((axes[1], "absres", "Size of the rank error", "|predicted - measured percentile|"), (axes[2], "resid", "Direction of the rank error", "predicted - measured percentile")):
    data = [x[x.cls == g][col].values for g in g_order]
    bp = ax.boxplot(data, vert=True, widths=0.6, showfliers=False, patch_artist=True, medianprops={"color": INK})
    for b_, g in zip(bp["boxes"], g_order):
        b_.set_facecolor(BLUE if g.startswith("Stable") else ORANGE); b_.set_alpha(0.45); b_.set_edgecolor(MUTED)
    ax.set_xticks(range(1, len(g_order) + 1), [g.replace("Stable (all three)", "Stable") for g in g_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(yl)
    ax.set_title(ttl, fontsize=9.5, loc="left")
    if col == "resid":
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.text(0.99, 0.97, "over-predicted", transform=ax.transAxes, ha="right", va="top", fontsize=8, color=MUTED)
        ax.text(0.99, 0.03, "under-predicted", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "x1_prediction_quality_by_class.png", bbox_inches="tight")
plt.close(fig)

# =============================================================== 2. shuffle effects as a grouped dot plot (cardiomyocyte sites, heart ATAC)
d = shuf["cm"]
d = d[d.trackset == "CM"]
sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
nat = d[d.condition == "native"].set_index("site_id")
rows, pvals = [], []
for sw in ("motif", "170", "peak"):
    e = d[d.condition == f"shuffle@{sw}"]
    for rw in ("motif", "170", "peak"):
        col = RD[rw]
        pc = 0.05 * nat[col].median()
        l = np.log2((e[col] + pc) / (e.site_id.map(nat[col]) + pc)).groupby(e.site_id).mean()
        a, b = l[sd.reindex(l.index) == "Stable"].dropna(), l[sd.reindex(l.index) == "Dynamic"].dropna()
        lo_s, hi_s = boot_med(a); lo_d, hi_d = boot_med(b)
        rows.append({"shuffle": sw, "readout": rw, "median_stable": a.median(), "lo_stable": lo_s, "hi_stable": hi_s, "median_dynamic": b.median(), "lo_dynamic": lo_d, "hi_dynamic": hi_d,
                     "p": mannwhitneyu(a, b).pvalue, "n_stable": len(a), "n_dynamic": len(b)})
SH = pd.DataFrame(rows)
SH["fdr"] = bh(SH.p)
SH["sig"] = SH.fdr.map(stars)
SH.to_csv(OUT / "tables" / "shuffle_dotplot_values.csv", index=False)
NAME = {"motif": "11 bp motif", "170": "170 bp", "peak": "whole peak"}
fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), sharey=True)
for ax, sw in zip(axes, ("motif", "170", "peak")):
    g = SH[SH.shuffle == sw].reset_index(drop=True)
    for i, r in g.iterrows():
        for off, key, col in ((-0.12, "stable", BLUE), (0.12, "dynamic", ORANGE)):
            ax.errorbar(i + off, r[f"median_{key}"], yerr=[[r[f"median_{key}"] - r[f"lo_{key}"]], [r[f"hi_{key}"] - r[f"median_{key}"]]], fmt="o", color=col, ms=6.5, lw=1.5, capsize=0,
                        label=("Stable sites" if key == "stable" else "Dynamic sites") if i == 0 else None)
        ax.text(i, 0.04, r.sig, ha="center", va="bottom", fontsize=8.5)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks(range(3), [f"read over\n{NAME[k]}" for k in ("motif", "170", "peak")])
    ax.set_xlim(-0.6, 2.6)
    ax.set_title(f"Shuffle the {NAME[sw]}", fontsize=10.5, loc="left")
axes[0].set_ylabel("Median effect on predicted accessibility (log2FC)")
axes[0].legend(frameon=False, loc="lower left")
axes[0].set_ylim(top=0.35)
fig.text(0.5, -0.03, "Each panel is one shuffle; the three positions are the region over which the effect is read. Median with 95% bootstrap CI. Marks: Stable vs Dynamic, Mann-Whitney, BH over the nine comparisons.", fontsize=8, color=MUTED, ha="center")
fig.tight_layout()
fig.savefig(OUT / "figures" / "x2_shuffle_dotplot.png", bbox_inches="tight")
plt.close(fig)

# =============================================================== 3. peak swap, Stable vs Dynamic donors
ps = pd.read_csv(R / "insertion_peakswap.csv")
ps["sd"] = ps.donor_stable_dynamic
pc = 0.05 * ps.donor_native_atac.median()
ps["retention"] = np.log2((ps.swapped_atac + pc) / (ps.donor_native_atac + pc))          # swapped accessibility relative to the donor's own native value
ps["vs_recipient"] = np.log2((ps.swapped_atac + pc) / (ps.recipient_native_atac + pc))     # relative to the recipient's own native value
dn = ps.groupby(["donor_id", "sd"]).agg(donor_native=("donor_native_atac", "first"), swapped=("swapped_atac", "mean"), retention=("retention", "mean"), width=("donor_peak_width", "first")).reset_index()
rows = []
for lab, col in (("donor native accessibility", "donor_native"), ("swapped accessibility (mean over 4 recipients)", "swapped"), ("retention: log2(swapped / donor native)", "retention"), ("donor peak width (bp)", "width")):
    a, b = dn[dn.sd == "Stable"][col], dn[dn.sd == "Dynamic"][col]
    rows.append({"quantity": lab, "stable_median": a.median(), "dynamic_median": b.median(), "n_stable_donors": len(a), "n_dynamic_donors": len(b), "p_mannwhitney": mannwhitneyu(a, b).pvalue})
PSX = pd.DataFrame(rows)
for g in ("Stable", "Dynamic"):
    s = ps[ps.sd == g]
    PSX.loc[len(PSX)] = [f"Spearman rho, swapped vs donor native ({g} donors, {len(s)} pairs)", spearmanr(s.swapped_atac, s.donor_native_atac)[0], np.nan, s.donor_id.nunique(), np.nan, spearmanr(s.swapped_atac, s.donor_native_atac)[1]]
ps["dyn"] = (ps.sd == "Dynamic").astype(int)
ps["l_sw"] = np.log2(ps.swapped_atac + pc); ps["l_dn"] = np.log2(ps.donor_native_atac + pc)
ols = smf.ols("l_sw ~ l_dn + dyn + C(recipient_id)", ps).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(ps.donor_id)[0]})
PSX.loc[len(PSX)] = ["log2 swapped ~ log2 donor native + Dynamic donor + recipient (coefficient on Dynamic, cluster-robust by donor)", ols.params["dyn"], ols.conf_int().loc["dyn", 0], np.nan, ols.conf_int().loc["dyn", 1], ols.pvalues["dyn"]]
PSX.to_csv(OUT / "tables" / "peak_swap_stable_vs_dynamic.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1, 1.25, 1.1]})
ax = axes[0]
for k, (g, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
    v = dn[dn.sd == g].donor_native
    ax.scatter(np.full(len(v), k) + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=col, alpha=0.65, linewidths=0)
    ax.hlines(v.median(), k - 0.28, k + 0.28, color=INK, lw=1.6)
ax.set_yscale("log"); ax.set_xticks([0, 1], [f"Stable donors\n(n={int((dn.sd == 'Stable').sum())})", f"Dynamic donors\n(n={int((dn.sd == 'Dynamic').sum())})"])
ax.set_ylabel("Donor peak's own native accessibility")
ax.set_title(f"Donors differ in starting accessibility\n(Mann-Whitney p = {PSX.iloc[0].p_mannwhitney:.1e})", fontsize=9.5, loc="left")
ax = axes[1]
for g, col in (("Stable", BLUE), ("Dynamic", ORANGE)):
    s = ps[ps.sd == g]
    ax.scatter(s.donor_native_atac, s.swapped_atac, s=13, color=col, alpha=0.5, linewidths=0, label=f"{g} donors (rho {spearmanr(s.swapped_atac, s.donor_native_atac)[0]:.2f})")
    xs = np.linspace(np.log10(s.donor_native_atac.min()), np.log10(s.donor_native_atac.max()), 50)
    sl = np.polyfit(np.log10(s.donor_native_atac), np.log10(s.swapped_atac), 1)
    ax.plot(10 ** xs, 10 ** np.polyval(sl, xs), color=col, lw=2)
lim = [min(ps.donor_native_atac.min(), ps.swapped_atac.min()) * 0.9, max(ps.donor_native_atac.max(), ps.swapped_atac.max()) * 1.1]
ax.plot(lim, lim, color=MUTED, lw=0.8, ls=":")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Donor peak's own native accessibility"); ax.set_ylabel("Accessibility after the swap")
ax.set_title(f"Same relation for both classes\n(Dynamic coefficient after adjustment {ols.params['dyn']:+.2f}, p = {ols.pvalues['dyn']:.2f})", fontsize=9.5, loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax = axes[2]
rec = ps.drop_duplicates("recipient_id").sort_values("recipient_native_atac").recipient_id.tolist()
for k, (g, col) in enumerate((("Stable", BLUE), ("Dynamic", ORANGE))):
    for i, rid in enumerate(rec):
        v = ps[(ps.sd == g) & (ps.recipient_id == rid)].retention
        lo_, hi_ = boot_med(v)
        ax.errorbar(i + (k - 0.5) * 0.2, v.median(), yerr=[[v.median() - lo_], [hi_ - v.median()]], fmt="o", color=col, ms=6, lw=1.4, capsize=0, label=f"{g} donors" if i == 0 else None)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xticks(range(len(rec)), [f"recipient {i + 1}\n(native {ps[ps.recipient_id == r].recipient_native_atac.iloc[0]:.2f})" for i, r in enumerate(rec)], fontsize=8)
ax.set_ylabel("log2(swapped / donor's own native)")
ax.set_title("What a peak keeps in each new location", fontsize=9.5, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "x3_peak_swap_stable_vs_dynamic.png", bbox_inches="tight")
plt.close(fig)

# =============================================================== 4. observed vs predicted accessibility: cardiomyocyte sites (heart), fibroblast sites (heart and Patski) side by side
TITLE4 = {"CM sites, heart ATAC": "Cardiomyocyte sites\nheart ATAC track vs measured cardiomyocyte snATAC", "FB sites, heart ATAC": "Fibroblast sites\nheart ATAC track vs measured fibroblast snATAC",
          "FB sites, Patski ATAC": "Fibroblast sites\nPatski ATAC track vs measured fibroblast snATAC"}
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
rows4 = []
for ax, (name, sk, ts, mcol) in zip(axes, COMBOS):
    x = resid_store[name]
    lx, ly = np.log10(x.pred.clip(lower=1e-3)), np.log10(x.meas.clip(lower=1e-3))
    hb = ax.hexbin(lx, ly, gridsize=38, cmap="Blues", mincnt=1, bins="log", linewidths=0.2, edgecolors="none")
    r = spearmanr(x.pred, x.meas)[0]
    b = rho_boot(x.pred.values, x.meas.values, 300)
    sl = np.polyfit(lx, ly, 1)
    xs = np.linspace(lx.min(), lx.max(), 50)
    ax.plot(xs, np.polyval(sl, xs), color=ORANGE, lw=2)
    ax.set_xlabel("Predicted native accessibility over the peak (log10)")
    ax.set_ylabel("Measured accessibility, CPM per kb (log10)")
    ax.set_title(TITLE4[name], fontsize=9.5, loc="left")
    ax.text(0.03, 0.96, f"Spearman rho = {r:.2f}\n95% CI {np.percentile(b, 2.5):.2f} to {np.percentile(b, 97.5):.2f}\nn = {len(x):,} sites", transform=ax.transAxes, va="top", fontsize=9,
            bbox={"fc": SURF, "ec": GRID, "alpha": 0.9})
    rows4.append({"combo": name, "rho": r, "ci_lo": np.percentile(b, 2.5), "ci_hi": np.percentile(b, 97.5), "n": len(x)})
fig.text(0.5, -0.03, "Each point cloud is the 2-D density of sites (darker = more sites); orange line = linear fit on the log scales. Sites without a measured peak are not shown.", fontsize=8, color=MUTED, ha="center")
fig.tight_layout()
fig.savefig(OUT / "figures" / "x4_observed_vs_predicted.png", bbox_inches="tight")
plt.close(fig)
pd.DataFrame(rows4).to_csv(OUT / "tables" / "observed_vs_predicted.csv", index=False)

pd.set_option("display.width", 220, "display.max_columns", 30)
print(CL.round(3).to_string(index=False)); print(PR.round(4).to_string(index=False)); print(TS.round(4).T.to_string()); print(SH.round(3).to_string(index=False)); print(PSX.round(4).to_string(index=False))
