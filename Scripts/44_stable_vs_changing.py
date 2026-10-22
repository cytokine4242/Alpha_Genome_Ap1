#!/usr/bin/env python3
"""The central question in one table: does ANYTHING separate a stable AP-1 site from a changing one?
Both definitions are tested side by side with identical matching, so a difference in the answer is a difference between the definitions, not the method.
  threshold:  Stable (bound in all three conditions) vs Dynamic (anything else)
  delta score: Steady vs Changing (largest fifth of |score change| in either contrast), direction ignored
  direction:  Age Up and Age Down (delta score) against Unchanged, where direction is kept
Measures fall into four blocks: the site's own AP-1 sequence, the neighbourhood sequence (PWM, needs no ATAC signal), the footprinted neighbourhood (needs signal), and
how the site responds to mutation in AlphaGenome. Matching: predicted accessibility at motif and peak, peak width, TSS distance, and measured accessibility.
Outputs Analysis/stable_vs_changing/{tables,figures}.   python3 44_stable_vs_changing.py"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "stable_vs_changing"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(44)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY, RED = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993", "#b3261e"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ---------------------------------------------------------------- site table
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
core = pd.read_csv(ROOT / "Data/derived/core_quality.csv").set_index("site_id")
pwm = pd.read_csv(ROOT / "Data/derived/pwm_counts.csv").set_index("site_id")
FD = pd.read_csv(ROOT / "Analysis/dynamic/tables/per_site_alphagenome_features.csv").set_index("site_id")
M = pd.read_csv(ROOT / "Analysis/distance_rules/tables/merged_partner_instances.csv.gz")
S = sl.join(scl[["age_score_class", "score_dynamic"]]).join(meas.rename("measured")).join(core).join(pwm)
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))
near = M[M.dist <= 100]
S["n_families_near"] = near.groupby("site_id").family.nunique().reindex(S.index).fillna(0)
S["mef2_footprint"] = S.index.isin(near[near.family == "MEF2"].site_id).astype(float)
S["cardiac_footprints"] = near[near.family.isin(["MEF2", "GATA", "NKX / homeobox cardiac", "TBX", "TEAD"])].groupby("site_id").family.nunique().reindex(S.index).fillna(0)
S = S.join(FD[["replace motif, read over motif", "shuffle 170 bp, read over 170 bp", "replace peak, read over peak", "170 bp shuffle ±300 bp away, read over motif"]])

BLOCKS = [
    ("The site's own AP-1 sequence", [
        ("exact_core", "has an exact TGA[C/G]TCA core", "bin"),
        ("core_mismatches", "mismatches to the consensus core", "num"),
        ("pwmmax_AP1_100", "best AP-1 PWM score within 100 bp", "num")]),
    ("Neighbourhood SEQUENCE (no ATAC signal used)", [
        ("pwmmax_MEF2_100", "best MEF2 PWM score within 100 bp", "num"),
        ("pwmmax_GATA_100", "best GATA PWM score within 100 bp", "num"),
        ("pwmmax_TEAD_100", "best TEAD PWM score within 100 bp", "num"),
        ("pwmmax_NFI_100", "best NFI PWM score within 100 bp", "num"),
        ("pwm_AP1_100", "AP-1 PWM hits within 100 bp", "num"),
        ("pwm_KLF_100", "KLF/SP PWM hits within 100 bp", "num")]),
    ("Neighbourhood FOOTPRINTS (needs ATAC signal)", [
        ("n_families_near", "TF families footprinted within 100 bp", "num"),
        ("cardiac_footprints", "cardiac families footprinted within 100 bp", "num"),
        ("mef2_footprint", "a MEF2 footprint within 100 bp", "bin")]),
    ("Response to mutation (AlphaGenome)", [
        ("replace motif, read over motif", "effect of replacing the motif", "num"),
        ("shuffle 170 bp, read over 170 bp", "effect of shuffling the 170 bp", "num"),
        ("replace peak, read over peak", "effect of replacing the peak", "num"),
        ("170 bp shuffle ±300 bp away, read over motif", "effect of an edit 300 bp away (control)", "num")]),
]
CONTRASTS = [("Threshold: Dynamic vs Stable", "stable_dynamic", "Dynamic", "Stable"),
             ("Δscore: Changing vs Steady", "score_dynamic", "Changing (score)", "Steady (score)"),
             ("Δscore: Age Up vs Unchanged", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Δscore: Age Down vs Unchanged", "age_score_class", "Age Down (score)", "Unchanged")]


def make_pairs(col, a_, b_):
    m_ = S[S[col].isin([a_, b_])].dropna(subset=MCOLS).copy()
    Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
    avail, pairs = set(m_.index[m_[col] == b_]), []
    for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < 0.5:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


rows = []
for lab, col, a_, b_ in CONTRASTS:
    P = make_pairs(col, a_, b_)
    bal = max(abs((S.loc[P.case, c].mean() - S.loc[P.control, c].mean()) / S[c].std()) for c in MCOLS)
    for block, measures in BLOCKS:
        for v, vlab, kind in measures:
            x, y = S[v].reindex(P.case).values.astype(float), S[v].reindex(P.control).values.astype(float)
            ok = ~(np.isnan(x) | np.isnan(y))
            if ok.sum() < 30:
                continue
            x, y = x[ok], y[ok]
            if kind == "bin":
                n10, n01 = int((x.astype(bool) & ~y.astype(bool)).sum()), int((~x.astype(bool) & y.astype(bool)).sum())
                p_ = binomtest(n10, n10 + n01).pvalue if n10 + n01 else np.nan
                eff = x.mean() - y.mean()
                disp = f"{x.mean():.1%} vs {y.mean():.1%}"
            else:
                d = x - y
                p_ = wilcoxon(d).pvalue if np.any(d != 0) else np.nan
                eff = float(np.mean(d))
                disp = f"{np.median(x):.3g} vs {np.median(y):.3g}"
            sd = np.nanstd(S[v].astype(float))
            rows.append({"contrast": lab, "block": block, "measure": vlab, "n_pairs": int(ok.sum()), "values": disp,
                         "difference": eff, "standardised": eff / sd if sd else np.nan, "p": p_, "balance_sd": bal})
R = pd.DataFrame(rows)
R["fdr"] = R.groupby("contrast").p.transform(lambda p: bh(p.fillna(1)))
R["detects"] = np.where(R.fdr < 0.05, "YES", "no")
R.to_csv(OUT / "tables" / "stable_vs_changing.csv", index=False)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, len(CONTRASTS), figsize=(4.0 * len(CONTRASTS), 6.2), sharey=True)
order = [(b, m[1]) for b, ms in BLOCKS for m in ms]
labels = [m for _, m in order]
for ax, (lab, _, _, _) in zip(axes, CONTRASTS):
    w = R[R.contrast == lab].set_index("measure").reindex(labels)
    ys = np.arange(len(labels))[::-1]
    sig = (w.fdr < 0.05).values
    ax.scatter(w.standardised[~sig], ys[~sig], s=38, color=GREY, alpha=0.6)
    ax.scatter(w.standardised[sig], ys[sig], s=62, color=ORANGE, edgecolors=INK, lw=0.7, zorder=4)
    ax.axvline(0, color=MUTED, lw=0.9)
    ax.set_xlim(-0.55, 0.55)
    ax.set_title(lab, fontsize=9.5, loc="left")
    ax.set_xlabel("standardised difference")
axes[0].set_yticks(np.arange(len(labels))[::-1], labels, fontsize=8)
b0 = 0
for block, ms in BLOCKS:
    b0 += len(ms)
    if b0 < len(labels):
        for ax in axes:
            ax.axhline(len(labels) - b0 - 0.5, color=MUTED, lw=0.8, ls=":")
axes[0].scatter([], [], s=56, color=ORANGE, edgecolors=INK, lw=0.7, label="FDR < 0.05")
axes[0].scatter([], [], s=38, color=GREY, label="not significant")
axes[0].legend(frameon=False, fontsize=8, loc="lower left")
fig.suptitle("Does anything separate a stable AP-1 site from a changing one?", fontsize=11, x=0.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT / "figures" / "w1_stable_vs_changing.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 48)
for lab, _, _, _ in CONTRASTS:
    w = R[R.contrast == lab]
    print(f"\n===== {lab}   ({int(w.n_pairs.max())} pairs, balance within {w.balance_sd.iloc[0]:.3f} SD)")
    print(w[["block", "measure", "values", "difference", "p", "fdr", "detects"]].round(4).to_string(index=False))
print("\n===== how many measures detect a difference, by contrast and block")
print(R.pivot_table(index="block", columns="contrast", values="detects", aggfunc=lambda v: f"{(v == 'YES').sum()}/{len(v)}").to_string())
