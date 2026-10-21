#!/usr/bin/env python3
"""Restrict the analyses to sites whose native ATAC prediction agrees with the measured snATAC accessibility, and ask whether the main results hold.

Tiers are defined ONLY from native prediction vs measurement (never from the perturbation effects), separately for each site set and track:
  A  all sites
  B  a measured peak exists and the native predicted ATAC (peak) is >= the FLOOR_Q quantile (default 0.05, a noise floor) of the eligible sites
  C  B and |rank residual| <= 1 SD, rank residual = predicted percentile - measured percentile
  D  B and |rank residual| <= 0.5 SD (strict)
Combinations: CM sites with the heart track, FB sites with the heart track, FB sites with the Patski track (each checked against the measured accessibility of
the cell type of the site table). Outputs: Analysis/agreement/tables, Analysis/agreement/figures.   python3 30_agreement_filtered.py
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
from scipy.stats import mannwhitneyu, spearmanr, t as tdist

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
import os
FLOOR_Q = float(os.environ.get("FLOOR_Q", "0.05"))
OUT = ROOT / "Analysis" / os.environ.get("OUT_NAME", "agreement")
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
R = ROOT / "Katana" / "results"
rng = np.random.default_rng(31)

INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
WIDTHS = ["motif", "170", "peak"]
RD = {"motif": "atac_motif", "170": "atac_170", "peak": "atac_peak"}
TIERS = ["A", "B", "C", "D"]
TIER_LABEL = {"A": "A all sites", "B": "B measured + floor", "C": "C + agree (1 SD)", "D": "D + agree (0.5 SD)"}
FQ = f"{FLOOR_Q:.0%}"
TIER_COL = {"A": GREY, "B": "#7aa7e0", "C": BLUE, "D": "#0f3f78"}


def stars(q):
    return "" if pd.isna(q) else "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def load(p):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / f"{p}_chunk*.csv")))], ignore_index=True)


def lfc_tab(d, ts, col):
    nat = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[(d.condition != "native") & (d.trackset == ts)].copy()
    e["l"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "condition"]).l.mean().unstack()


def boot_diff(a, b, n=500):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.median(rng.choice(a, (n, len(a))), 1) - np.median(rng.choice(b, (n, len(b))), 1)
    return tuple(np.percentile(d, [2.5, 97.5]))


def boot_mean_ratio(num, den, n=500):
    num, den = np.asarray(num, float), np.asarray(den, float)
    idx = rng.integers(0, len(num), (n, len(num)))
    r = num[idx].mean(1) / den[idx].mean(1)
    return tuple(np.percentile(r, [2.5, 97.5]))


# ---------------- data
shuf = {"cm": load("shuffle_cm"), "fb": load("shuffle_fb")}
mut = {"cm": load("mut_cm"), "fb": load("mut_fb")}
ctl = load("controls_cm")
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
COMBOS = (("CM sites, heart ATAC", "cm", "CM", "cm_cpm_per_kb"), ("FB sites, heart ATAC", "fb", "CM", "fb_cpm_per_kb"), ("FB sites, Patski ATAC", "fb", "FB", "fb_cpm_per_kb"))

# ---------------- tiers
TIERMAP, ret_rows, scat, DFS = {}, [], {}, {}
for name, sites, ts, mcol in COMBOS:
    d = shuf[sites]
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    nat = d[(d.condition == "native") & (d.trackset == ts)].set_index("site_id").atac_peak
    m = meas[meas.site_table == sites].set_index("site_id")[mcol]
    df = pd.DataFrame({"pred": nat, "meas": m, "cls": sd})
    elig = df.dropna(subset=["pred", "meas"])
    pr, mr = elig.pred.rank(pct=True), elig.meas.rank(pct=True)
    resid = pr - mr
    sdr = resid.std()
    floor = elig.pred >= elig.pred.quantile(FLOOR_Q)
    tiers = {"A": df.index, "B": elig.index[floor], "C": elig.index[floor & (resid.abs() <= sdr)], "D": elig.index[floor & (resid.abs() <= 0.5 * sdr)]}
    TIERMAP[name] = tiers
    DFS[name] = (df, elig, resid, floor, sdr)
    scat[name] = (elig.assign(resid=resid, floor=floor, inC=elig.index.isin(tiers["C"])), sdr)
    rho = spearmanr(elig.pred, elig.meas)[0]
    for t_ in TIERS:
        ids = tiers[t_]
        sub = df.loc[ids]
        r_t = spearmanr(sub.dropna(subset=["pred", "meas"]).pred, sub.dropna(subset=["pred", "meas"]).meas)[0] if t_ != "A" else rho
        ret_rows.append({"combo": name, "tier": t_, "n": len(ids), "n_stable": int((sub.cls == "Stable").sum()), "n_dynamic": int((sub.cls == "Dynamic").sum()),
                         "retained_all": len(ids) / len(df), "retained_stable": (sub.cls == "Stable").sum() / (df.cls == "Stable").sum(),
                         "retained_dynamic": (sub.cls == "Dynamic").sum() / (df.cls == "Dynamic").sum(), "rho_pred_vs_meas_in_tier": r_t, "rank_resid_sd": sdr})
RET = pd.DataFrame(ret_rows)
RET.to_csv(OUT / "tables" / "retention_by_tier.csv", index=False)

# ---------------- Stable vs Dynamic in each tier (matched shuffle readouts)
rows = []
for name, sites, ts, _ in COMBOS:
    d = shuf[sites]
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    for w in WIDTHS:
        e = lfc_tab(d, ts, RD[w])[f"shuffle@{w}"]
        for t_ in TIERS:
            ids = TIERMAP[name][t_]
            s_ = e[e.index.isin(ids) & (sd == "Stable").reindex(e.index).fillna(False)].dropna()
            d_ = e[e.index.isin(ids) & (sd == "Dynamic").reindex(e.index).fillna(False)].dropna()
            lo, hi = boot_diff(d_, s_)
            rows.append({"combo": name, "width": w, "tier": t_, "n_stable": len(s_), "n_dynamic": len(d_), "median_stable": s_.median(), "median_dynamic": d_.median(),
                         "diff_dyn_minus_stable": d_.median() - s_.median(), "ci_lo": lo, "ci_hi": hi, "p": mannwhitneyu(d_, s_).pvalue})
SV = pd.DataFrame(rows)
SV["fdr"] = SV.p * 0 + np.nan
for (c_, w_), g in SV.groupby(["combo", "width"]):
    pass
# BH over all comparisons within a combo
from numpy import argsort


def bh(p):
    p = np.asarray(p, float)
    o = argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


for c_, ix in SV.groupby("combo").groups.items():
    SV.loc[ix, "fdr"] = bh(SV.loc[ix, "p"])
SV["sig"] = SV.fdr.map(stars)
SV.to_csv(OUT / "tables" / "stable_vs_dynamic_by_tier.csv", index=False)

# ---------------- mutation arm and controls in each tier
rows = []
for name, sites, ts, _ in COMBOS:
    d = mut[sites]
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    cache = {w: -lfc_tab(d, ts, RD[w]) for w in WIDTHS}
    for t_ in TIERS:
        ids = pd.Index(TIERMAP[name][t_])
        st = sd[sd == "Stable"].index.intersection(ids)
        dy = sd[sd == "Dynamic"].index.intersection(ids)
        rec = {"combo": name, "tier": t_, "n_stable": len(st), "n_dynamic": len(dy)}
        for e, w in (("random@motif", "motif"), ("core@motif", "motif"), ("random@170", "170"), ("flank@170", "170")):
            rec[f"{e}_stable"] = cache[w].loc[st, e].mean()
            rec[f"{e}_dynamic"] = cache[w].loc[dy, e].mean()
        a, b = cache["170"].loc[st, "flank@170"].dropna(), cache["170"].loc[st, "random@170"].dropna()
        ix = a.index.intersection(b.index)
        rec["ratio_flank170_over_random170"] = a[ix].mean() / b[ix].mean()
        rec["ratio_lo"], rec["ratio_hi"] = boot_mean_ratio(a[ix].values, b[ix].values)
        a, b = cache["motif"].loc[st, "core@motif"].dropna(), cache["motif"].loc[st, "random@motif"].dropna()
        ix = a.index.intersection(b.index)
        rec["ratio_core_over_random_motif"] = a[ix].mean() / b[ix].mean()
        rows.append(rec)
MU = pd.DataFrame(rows)
MU.to_csv(OUT / "tables" / "mutation_arm_by_tier.csv", index=False)

rows = []
name, sites, ts, _ = COMBOS[0]
dctl = ctl
sdc = dctl.drop_duplicates("site_id").set_index("site_id").stable_dynamic
tm = lfc_tab(dctl, "CM", "atac_motif")
for t_ in TIERS:
    ids = pd.Index(TIERMAP[name][t_])
    st = sdc[sdc == "Stable"].index.intersection(ids)
    cen = tm.loc[st, "titrate@170"]
    off = tm.loc[st, ["offset@-300", "offset@+300"]].mean(1)
    rows.append({"tier": t_, "n_stable": len(st), "centred_170": cen.median(), "offset_pm300": off.median(), "share_centred_more_negative": float((cen < off).mean()),
                 "core7_share_of_170": tm.loc[st, "titrate@7"].median() / cen.median(), "motif11_share_of_170": tm.loc[st, "titrate@11"].median() / cen.median()})
CT = pd.DataFrame(rows)
CT.to_csv(OUT / "tables" / "controls_by_tier.csv", index=False)

# ---------------- figures
# F1: which sites are retained (CM sites, heart): predicted vs measured with tier C highlighted
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
for ax, (name, sites, ts, mcol) in zip(axes, COMBOS):
    e, sdr = scat[name]
    ax.scatter(e.pred[~e.inC], e.meas[~e.inC], s=5, alpha=0.3, color=GREY, linewidths=0, label="excluded")
    ax.scatter(e.pred[e.inC], e.meas[e.inC], s=5, alpha=0.4, color=BLUE, linewidths=0, label="tier C (kept)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Predicted native ATAC (peak)")
    ax.set_ylabel("Measured snATAC (CPM per kb)")
    r_ = spearmanr(e.pred, e.meas)[0]
    ax.set_title(f"{name}: rho = {r_:.2f}; kept {e.inC.mean():.0%} of sites with a peak", fontsize=9, loc="left")
axes[0].legend(frameon=False, markerscale=3, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "e1_tier_C_selection.png", bbox_inches="tight")
plt.close(fig)

# F2: retention by class
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.2), sharey=True)
for ax, (name, *_x) in zip(axes, COMBOS):
    for k, (cls, col) in enumerate((("stable", BLUE), ("dynamic", ORANGE))):
        v = [RET[(RET.combo == name) & (RET.tier == t)][f"retained_{cls}"].iloc[0] for t in TIERS]
        ax.bar(np.arange(4) + (k - 0.5) * 0.38, v, 0.34, color=col, alpha=0.85, label=cls.capitalize())
        for i, x in enumerate(v):
            ax.text(i + (k - 0.5) * 0.38, x + 0.015, f"{x:.0%}", ha="center", fontsize=7.5)
    ax.set_xticks(range(4), TIERS)
    ax.set_title(name, fontsize=10, loc="left")
    ax.set_xlabel("Tier")
axes[0].set_ylabel("Share of the class retained")
axes[0].legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figures" / "e2_retention_by_class.png", bbox_inches="tight")
plt.close(fig)

# F3: Dynamic - Stable by tier (forest), three widths
fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.6), sharey=True)
ylabels, ypos = [], []
for ax, w in zip(axes, WIDTHS):
    y = 0
    for name, *_x in COMBOS:
        for t_ in TIERS[::-1]:
            r_ = SV[(SV.combo == name) & (SV.width == w) & (SV.tier == t_)].iloc[0]
            ax.errorbar(r_.diff_dyn_minus_stable, y, xerr=[[r_.diff_dyn_minus_stable - r_.ci_lo], [r_.ci_hi - r_.diff_dyn_minus_stable]], fmt="o", color=TIER_COL[t_], ms=5, lw=1.2, capsize=0)
            ax.text(1.02, y, r_.sig, transform=ax.get_yaxis_transform(), va="center", fontsize=7.5)
            if w == "motif":
                ylabels.append(f"{name} | tier {t_} (n={r_.n_stable + r_.n_dynamic:,})"); ypos.append(y)
            y += 1
        y += 0.7
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_title(f"{ {'motif': 'Motif', '170': '170 bp', 'peak': 'Peak'}[w] } shuffle", fontsize=10, loc="left")
    ax.set_xlabel("Dynamic - Stable (median log2FC, 95% CI)")
axes[0].set_yticks(ypos, ylabels, fontsize=7.5)
fig.text(0.01, -0.02, "Tier A (grey) = all sites; B = measured peak + predicted above the floor; C and D add rank agreement within 1 and 0.5 SD. Marks: Mann-Whitney, BH within each site set and track.", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "e3_stable_vs_dynamic_by_tier.png", bbox_inches="tight")
plt.close(fig)

# F4: decomposition and control summaries across tiers
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.3))
ax = axes[0]
for k, (name, *_x) in enumerate(COMBOS):
    g = MU[MU.combo == name]
    ax.errorbar(np.arange(4) + (k - 1) * 0.16, g.ratio_flank170_over_random170, yerr=[g.ratio_flank170_over_random170 - g.ratio_lo, g.ratio_hi - g.ratio_flank170_over_random170], fmt="o", ms=5, lw=1.2, capsize=0,
                color=[BLUE, AQUA, ORANGE][k], label=name)
ax.set_xticks(range(4), TIERS); ax.set_ylim(0.5, 1.1); ax.axhline(1, color=MUTED, lw=0.8, ls=":")
ax.set_ylabel("Surroundings-only loss / full 170 bp loss"); ax.set_xlabel("Tier"); ax.legend(frameon=False, fontsize=7.5)
ax.set_title("Mutation arm: surroundings carry most of it", fontsize=9.5, loc="left")
ax = axes[1]
ax.plot(range(4), CT.core7_share_of_170, "o-", color=BLUE, label="7 bp core"); ax.plot(range(4), CT.motif11_share_of_170, "s-", color=AQUA, label="11 bp motif")
ax.set_xticks(range(4), TIERS); ax.set_ylim(0, 0.6); ax.set_xlabel("Tier"); ax.set_ylabel("Share of the 170 bp effect (Stable, CM sites)"); ax.legend(frameon=False)
ax.set_title("Titration: how much the core gives", fontsize=9.5, loc="left")
ax = axes[2]
ax.plot(range(4), CT.centred_170, "o-", color=ORANGE, label="centred 170 bp"); ax.plot(range(4), CT.offset_pm300, "s-", color=GREY, label="same size, +-300 bp away")
ax.set_xticks(range(4), TIERS); ax.set_xlabel("Tier"); ax.set_ylabel("Median effect on the motif (log2FC)"); ax.legend(frameon=False)
ax.set_title("Controls: specific to the site", fontsize=9.5, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "e4_other_results_by_tier.png", bbox_inches="tight")
plt.close(fig)


# ---------------- QC figures: how the restriction was done and what was kept vs discarded
STATUS = ["kept (tier C)", "below noise floor", "rank disagreement", "no measured peak"]
SCOL = {"kept (tier C)": BLUE, "below noise floor": ORANGE, "rank disagreement": "#b0aea8", "no measured peak": "#d9d8d3"}
status_rows, funnel_rows = [], []
STAT = {}
for name, sites, ts, mcol in COMBOS:
    df, elig, resid, floor, sdr = DFS[name]
    st = pd.Series("no measured peak", index=df.index)
    st[elig.index] = np.where(~floor, "below noise floor", np.where(resid.abs() <= sdr, "kept (tier C)", "rank disagreement"))
    STAT[name] = st
    for cls in ("Stable", "Dynamic"):
        for s_ in STATUS:
            status_rows.append({"combo": name, "class": cls, "status": s_, "n": int(((st == s_) & (df.cls == cls)).sum()), "n_class": int((df.cls == cls).sum())})
    for step, ids in (("All sites", df.index), ("Has a measured peak", elig.index), (f"Above {FQ} floor (tier B)", TIERMAP[name]["B"]), ("Within 1 SD (tier C)", TIERMAP[name]["C"]), ("Within 0.5 SD (tier D)", TIERMAP[name]["D"])):
        for cls in ("Stable", "Dynamic"):
            funnel_rows.append({"combo": name, "step": step, "class": cls, "n": int((df.loc[df.index.isin(ids), "cls"] == cls).sum())})
pd.DataFrame(status_rows).to_csv(OUT / "tables" / "status_counts.csv", index=False)
pd.DataFrame(funnel_rows).to_csv(OUT / "tables" / "funnel_counts.csv", index=False)

# Q1: funnel
FUN = pd.DataFrame(funnel_rows)
steps = list(dict.fromkeys(FUN.step))
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
for ax, (name, *_x) in zip(axes, COMBOS):
    g = FUN[FUN.combo == name]
    s_ = g[g["class"] == "Stable"].set_index("step").n.reindex(steps)
    d_ = g[g["class"] == "Dynamic"].set_index("step").n.reindex(steps)
    y = np.arange(len(steps))[::-1]
    ax.barh(y, s_.values, color=BLUE, alpha=0.9, label="Stable")
    ax.barh(y, d_.values, left=s_.values, color=ORANGE, alpha=0.9, label="Dynamic")
    for yi, a_, b_ in zip(y, s_.values, d_.values):
        ax.text(a_ + b_ + 40, yi, f"{a_ + b_:,}  (Dyn {b_ / (a_ + b_):.0%})", va="center", fontsize=7.5)
    ax.set_yticks(y, steps, fontsize=8)
    ax.set_xlim(0, (s_ + d_).max() * 1.45)
    ax.set_title(name, fontsize=10, loc="left")
    ax.set_xlabel("Sites")
axes[2].legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "e0_funnel.png", bbox_inches="tight")
plt.close(fig)

# Q2: rank-space map showing why each site was kept or dropped
fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, (name, *_x) in zip(axes, COMBOS):
    df, elig, resid, floor, sdr = DFS[name]
    mr, pr = elig.meas.rank(pct=True), elig.pred.rank(pct=True)
    st = STAT[name].loc[elig.index]
    for s_ in ("rank disagreement", "below noise floor", "kept (tier C)"):
        k = st == s_
        ax.scatter(mr[k], pr[k], s=4, alpha=0.35, color=SCOL[s_], linewidths=0, label=f"{s_} ({k.sum():,})")
    xs = np.linspace(0, 1, 50)
    ax.plot(xs, xs, color=INK, lw=0.8)
    ax.plot(xs, xs + sdr, color=INK, lw=0.8, ls="--")
    ax.plot(xs, xs - sdr, color=INK, lw=0.8, ls="--")
    ax.axhline(FLOOR_Q, color=ORANGE, lw=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Measured accessibility (percentile)")
    ax.set_ylabel("Predicted accessibility (percentile)")
    nopeak = int((STAT[name] == "no measured peak").sum())
    ax.set_title(f"{name} ({nopeak:,} sites without a measured peak not shown)", fontsize=8.5, loc="left")
    ax.legend(frameon=False, markerscale=3, fontsize=7.5, loc="lower right")
fig.text(0.01, -0.03, f"Solid diagonal = perfect agreement. Dashed lines = +-1 SD of the rank difference (the tier C band). Orange line = the {FQ} floor on predicted accessibility. Blue = kept.", fontsize=8, color=MUTED)
fig.tight_layout()
fig.savefig(OUT / "figures" / "e0b_rank_map.png", bbox_inches="tight")
plt.close(fig)

# Q3: is the filter biased? share kept by accessibility, per class
QROWS = []
fig, axes = plt.subplots(2, 3, figsize=(13, 6))
for j, (name, *_x) in enumerate(COMBOS):
    df, elig, resid, floor, sdr = DFS[name]
    e = elig.copy()
    e["kept"] = e.index.isin(TIERMAP[name]["C"])
    for i, (col, lab) in enumerate((("pred", "Predicted native accessibility"), ("meas", "Measured accessibility"))):
        ax = axes[i, j]
        e["q"] = pd.qcut(e[col].rank(method="first"), 5, labels=False)
        for cls, c_ in (("Stable", BLUE), ("Dynamic", ORANGE)):
            g = e[e.cls == cls].groupby("q").kept.agg(["mean", "size"])
            ax.plot(g.index + 1, g["mean"], "o-", color=c_, label=f"{cls}", ms=4)
            for q_, r_ in g.iterrows():
                QROWS.append({"combo": name, "by": col, "class": cls, "quintile": int(q_) + 1, "share_kept": r_["mean"], "n": int(r_["size"])})
        ax.set_ylim(0, 1)
        ax.set_xticks(range(1, 6), ["Q1\nlowest", "Q2", "Q3", "Q4", "Q5\nhighest"])
        ax.set_xlabel(f"{lab} quintile")
        ax.set_ylabel("Share kept (tier C)")
        if i == 0:
            ax.set_title(name, fontsize=10, loc="left")
axes[0, 0].legend(frameon=False)
fig.tight_layout()
pd.DataFrame(QROWS).to_csv(OUT / "tables" / "kept_by_quintile.csv", index=False)
fig.savefig(OUT / "figures" / "e0c_kept_by_accessibility.png", bbox_inches="tight")
plt.close(fig)

# Q4: rank difference by class, and predicted level kept vs discarded
fig, axes = plt.subplots(2, 3, figsize=(13, 6))
for j, (name, *_x) in enumerate(COMBOS):
    df, elig, resid, floor, sdr = DFS[name]
    ax = axes[0, j]
    bins = np.linspace(-1, 1, 41)
    for cls, c_ in (("Stable", BLUE), ("Dynamic", ORANGE)):
        ax.hist(resid[elig.cls == cls], bins=bins, density=True, alpha=0.5, color=c_, label=cls)
    for s_ in (-sdr, sdr):
        ax.axvline(s_, color=INK, lw=0.8, ls="--")
    for s_ in (-0.5 * sdr, 0.5 * sdr):
        ax.axvline(s_, color=INK, lw=0.8, ls=":")
    ax.set_xlabel("Predicted minus measured percentile")
    ax.set_ylabel("Density")
    ax.set_title(f"{name} (dashed: 1 SD = {sdr:.2f}; dotted: 0.5 SD)", fontsize=8.5, loc="left")
    ax = axes[1, j]
    kept = elig.index.isin(TIERMAP[name]["C"])
    lp = np.log10(elig.pred.clip(lower=1e-6))
    bins = np.linspace(lp.quantile(0.005), lp.quantile(0.995), 40)
    ax.hist(lp[~kept], bins=bins, density=True, alpha=0.5, color=GREY, label="discarded")
    ax.hist(lp[kept], bins=bins, density=True, alpha=0.5, color=BLUE, label="kept (tier C)")
    ax.axvline(np.log10(elig.pred.quantile(FLOOR_Q)), color=ORANGE, lw=1)
    ax.set_xlabel("log10 predicted native ATAC (peak)")
    ax.set_ylabel("Density")
    if j == 0:
        axes[0, 0].legend(frameon=False); ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figures" / "e0d_distributions.png", bbox_inches="tight")
plt.close(fig)

# Q5: transparency check, do kept and discarded sites differ in the perturbation effect? (reported, not used for selection)
rows = []
for name, sites, ts, _ in COMBOS:
    d = shuf[sites]
    sd = d.drop_duplicates("site_id").set_index("site_id").stable_dynamic
    e170 = lfc_tab(d, ts, RD["170"])["shuffle@170"]
    ids = pd.Index(TIERMAP[name]["C"])
    elig = DFS[name][1]
    for cls in ("Stable", "Dynamic"):
        in_ = e170[(e170.index.isin(ids)) & (sd.reindex(e170.index) == cls)].dropna()
        out_ = e170[(~e170.index.isin(ids)) & (e170.index.isin(elig.index)) & (sd.reindex(e170.index) == cls)].dropna()
        rows.append({"combo": name, "class": cls, "n_kept": len(in_), "n_discarded_with_peak": len(out_), "median_effect_kept": in_.median(), "median_effect_discarded": out_.median(),
                     "p": mannwhitneyu(in_, out_).pvalue})
KD = pd.DataFrame(rows)
KD.to_csv(OUT / "tables" / "effect_kept_vs_discarded.csv", index=False)
print("\nEFFECT kept vs discarded (170 bp shuffle)\n", KD.round(3).to_string(index=False))

pd.set_option("display.width", 230, "display.max_columns", 40)
print("RETENTION\n", RET.round(3).to_string(index=False))
print("\nSTABLE vs DYNAMIC by tier (170 bp and others)\n", SV[["combo", "width", "tier", "n_stable", "n_dynamic", "median_stable", "median_dynamic", "diff_dyn_minus_stable", "ci_lo", "ci_hi", "sig"]].round(3).to_string(index=False))
print("\nMUTATION ARM by tier\n", MU.round(3).to_string(index=False))
print("\nCONTROLS by tier\n", CT.round(3).to_string(index=False))
