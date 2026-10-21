#!/usr/bin/env python3
"""Audit of the Stable-vs-Dynamic gap after the independent review: the gap was adjusted for native PEAK accessibility for every readout, but the effect at a
given width depends on native accessibility at that SAME readout region. This script re-estimates the gap with the matched native value, flexible (spline)
terms, extra covariates (peak width, distance to TSS, footprint score) and standard errors clustered by peak, and stratifies by deciles of the matched native
value, reporting how many Dynamic sites sit in each decile. A common-support estimate weights the decile gaps by the Dynamic count (the gap for sites like Dynamic ones).
Outputs Analysis/audit/{tables,figures}.   python3 31_adjusted_gap_audit.py
"""
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

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "audit"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
R = ROOT / "Katana" / "results"
rng = np.random.default_rng(7)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
RD = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}


def load(p):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{p}_chunk*.csv")))], ignore_index=True)


COMBOS = (("CM sites, heart ATAC", "cm", "CM"), ("FB sites, heart ATAC", "fb", "CM"), ("FB sites, Patski ATAC", "fb", "FB"))
sites = {"cm": pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id"), "fb": pd.read_csv(ROOT / "Data/derived/fb_fos_sites.csv").set_index("site_id")}
for k, s in sites.items():
    fp = np.where(s.in_YS, s.YS_score, -1.0)
    for c_, sc in (("in_MAS", "MAS_score"), ("in_MAE", "MAE_score")):
        fp = np.maximum(fp, np.where(s[c_], s[sc], -1.0))
    s["fp_score"] = fp
shuf = {"cm": load("shuffle_cm"), "fb": load("shuffle_fb")}

rows, dec_rows = [], []
for name, sk, ts in COMBOS:
    d = shuf[sk]
    d = d[d.trackset == ts]
    tss = d.drop_duplicates("site_id").set_index("site_id").tss_dist
    nat = d[d.condition == "native"].set_index("site_id")
    for w in RD:
        col = RD[w]
        pc = 0.05 * nat[col].median()
        e = d[d.condition == f"shuffle@{w}"].copy()
        e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat[col]) + pc))
        x = pd.DataFrame({"y": e.groupby("site_id").l.mean(), "nm": nat[col], "npk": nat["atac_peak"]})
        x = x.join(sites[sk][["stable_dynamic", "peak_width", "peak_coord", "fp_score"]]).join(tss.abs().rename("abs_tss")).dropna()
        x = x[x.stable_dynamic.isin(["Stable", "Dynamic"])]
        x["dyn"] = (x.stable_dynamic == "Dynamic").astype(int)
        x["lm"] = np.log10(x.nm.clip(lower=1e-6)); x["lpk"] = np.log10(x.npk.clip(lower=1e-6)); x["lw"] = np.log10(x.peak_width); x["lt"] = np.log10(x.abs_tss + 1)
        x["grp"] = pd.factorize(x.peak_coord)[0]
        models = [("raw", "y ~ dyn"), ("+ native PEAK accessibility, linear (earlier model)", "y ~ dyn + lpk"), ("+ native accessibility at the matched readout, linear", "y ~ dyn + lm"),
                  ("+ matched native, spline", "y ~ dyn + bs(lm, df=5)"), ("+ matched native, peak native, peak width, TSS distance (splines)", "y ~ dyn + bs(lm, df=5) + bs(lpk, df=4) + bs(lw, df=4) + bs(lt, df=4)"),
                  ("+ all of the above and footprint score", "y ~ dyn + bs(lm, df=5) + bs(lpk, df=4) + bs(lw, df=4) + bs(lt, df=4) + bs(fp_score, df=4)")]
        for lab, f in models:
            r = smf.ols(f, x).fit(cov_type="cluster", cov_kwds={"groups": x.grp})
            ci = r.conf_int().loc["dyn"]
            rows.append({"combo": name, "width": w, "model": lab, "gap": r.params["dyn"], "ci_lo": ci[0], "ci_hi": ci[1], "p": r.pvalues["dyn"], "n": len(x), "n_dynamic": int(x.dyn.sum())})
        if w == "170":
            x["dec"] = pd.qcut(x.nm.rank(method="first"), 10, labels=False)
            for k, g in x.groupby("dec"):
                a, b = g[g.dyn == 1].y, g[g.dyn == 0].y
                dec_rows.append({"combo": name, "decile": int(k) + 1, "native_median": g.nm.median(), "n_dynamic": len(a), "n_stable": len(b),
                                 "median_dynamic": a.median() if len(a) else np.nan, "median_stable": b.median(), "gap": a.median() - b.median() if len(a) else np.nan})
G = pd.DataFrame(rows)
G.to_csv(OUT / "tables" / "gap_models.csv", index=False)
D = pd.DataFrame(dec_rows)
D.to_csv(OUT / "tables" / "gap_by_matched_native_decile.csv", index=False)
cs = []
for name, g in D.groupby("combo"):
    g = g[g.n_dynamic > 0]
    cs.append({"combo": name, "common_support_gap": np.average(g.gap, weights=g.n_dynamic), "dynamic_in_lowest_5_deciles": g[g.decile <= 5].n_dynamic.sum() / g.n_dynamic.sum()})
CS = pd.DataFrame(cs)
CS.to_csv(OUT / "tables" / "common_support_gap.csv", index=False)

# figure 1: gap across models
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
labs = list(dict.fromkeys(G.model))
for ax, (name, *_x) in zip(axes, COMBOS):
    g = G[(G.combo == name) & (G.width == "170")].set_index("model").reindex(labs)
    y = np.arange(len(labs))[::-1]
    ax.errorbar(g.gap, y, xerr=[g.gap - g.ci_lo, g.ci_hi - g.gap], fmt="o", color=BLUE, ms=5, lw=1.3, capsize=0)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(y, [l if i == 0 or name == COMBOS[0][0] else "" for i, l in enumerate(labs)] if False else (labs if ax is axes[0] else [""] * len(labs)), fontsize=7.5)
    ax.set_title(name, fontsize=10, loc="left")
    ax.set_xlabel("Dynamic - Stable at the 170 bp shuffle (log2FC, 95% CI)")
fig.text(0.01, -0.03, "Standard errors clustered by TOBIAS peak. Footprint score is partly a defining variable of the classes, so the last row is a sensitivity check, not a causal adjustment.", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "h1_gap_models.png", bbox_inches="tight")
plt.close(fig)

# figure 2: gap by decile of matched native accessibility, with Dynamic counts
fig, axes = plt.subplots(2, 3, figsize=(13, 5.6), gridspec_kw={"height_ratios": [3, 1.2]}, sharex=True)
for j, (name, *_x) in enumerate(COMBOS):
    g = D[D.combo == name]
    ax = axes[0, j]
    ax.plot(g.decile, g.median_stable, "o-", color=BLUE, label="Stable", ms=4)
    ax.plot(g.decile, g.median_dynamic, "o-", color=ORANGE, label="Dynamic", ms=4)
    ax.set_title(f"{name}: gap for sites like Dynamic = {CS[CS.combo == name].common_support_gap.iloc[0]:+.2f}", fontsize=9, loc="left")
    ax.set_ylabel("Median 170 bp shuffle effect")
    ax = axes[1, j]
    ax.bar(g.decile - 0.2, g.n_stable, 0.4, color=BLUE, alpha=0.8)
    ax.bar(g.decile + 0.2, g.n_dynamic, 0.4, color=ORANGE, alpha=0.9)
    ax.set_xlabel("Native accessibility decile at 170 bp (1 = lowest)")
    ax.set_ylabel("Sites")
    ax.set_xticks(range(1, 11))
axes[0, 0].legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figures" / "h2_gap_by_matched_decile.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 220, "display.max_columns", 30, "display.max_colwidth", 70)
print(G[G.width == "170"][["combo", "model", "gap", "ci_lo", "ci_hi", "p", "n", "n_dynamic"]].round(3).to_string(index=False))
print(G[(G.width != "170") & (G.model.str.startswith(("raw", "+ matched native, peak")))][["combo", "width", "model", "gap", "ci_lo", "ci_hi", "p"]].round(3).to_string(index=False))
print(CS.round(3).to_string(index=False))
