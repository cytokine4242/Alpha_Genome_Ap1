#!/usr/bin/env python3
"""Everything the mutation scan measures about a site, tested across every class definition, on all 3,678 cardiomyocyte AP-1 elements.
Features come from the tile scan (10 bp tiles, +-300 bp, 3 replicates) and the saturation mutagenesis (every alternative base, +-30 bp).
Contrasts: Dynamic vs Stable and the four directional threshold classes against Stable, plus the delta-score classes against Unchanged.
Every contrast uses 1:1 matching on predicted accessibility (motif and peak), peak width, TSS distance and MEASURED accessibility, so footprint
detection and chromatin context are held constant; BH within each contrast.
Outputs Analysis/scan_features/{tables,figures}.   python3 46_scan_features_by_class.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "scan_features"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(46)
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


# ---------------------------------------------------------------- scan features, per site
fs = sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))
sc = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
sc = sc[sc.trackset == "CM"]
n_pass = sc.groupby("site_id").size()
sc = sc[sc.site_id.isin(n_pass.index[n_pass == n_pass.max()])]
print(f"scan: {sc.site_id.nunique():,} complete sites from {len(fs)} files")

F = {}
for rd, col in (("170", "atac_170"), ("motif", "atac_motif")):
    nat = sc[sc.kind == "native"].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = sc[sc.kind != "native"].copy()
    e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    tile = e[e.kind == "tile"].groupby(["site_id", "pos"], as_index=False).eff.mean()
    tile["centre"] = tile.pos + 5
    thr = float(np.percentile(tile[tile.centre.abs() >= 250].eff, 1))
    rows = []
    for sid, w in tile.groupby("site_id"):
        w = w.sort_values("centre")
        c, v = w.centre.values, w.eff.values
        own = v[np.isin(w.pos.values, [-10, 0])].min()
        out = ~np.isin(w.pos.values, [-10, 0])
        loss = np.clip(-v, 0, None)
        tot = loss.sum()
        share = loss / tot if tot > 0 else np.zeros_like(loss)
        srt = np.sort(loss)
        gini = float((2 * np.arange(1, len(srt) + 1) - len(srt) - 1).dot(srt) / (len(srt) * srt.sum())) if tot > 0 else np.nan
        centroid = float((share * c).sum())
        left, right = loss[c < -5].sum(), loss[c > 5].sum()
        rows.append({"site_id": sid, f"own_motif_{rd}": own, f"best_partner_{rd}": v[out].min(),
                     f"partner_dominates_{rd}": float(v[out].min() < own), f"n_tiles_matter_{rd}": int((v < thr).sum()),
                     f"total_dependence_{rd}": tot, f"gini_{rd}": gini,
                     f"spread_{rd}": float(np.sqrt((share * (c - centroid) ** 2).sum())),
                     f"share_within_85_{rd}": float(share[np.abs(c) <= 85].sum()),
                     f"asymmetry_{rd}": float((left - right) / tot) if tot > 0 else np.nan})
    F[rd] = pd.DataFrame(rows).set_index("site_id")
    if rd == "170":
        im = e[e.kind == "ism"].groupby(["site_id", "pos"], as_index=False).eff.agg(imp=lambda v: np.abs(v).mean())
        irows = []
        for sid, w in im.groupby("site_id"):
            p_, v_ = w.pos.values, w.imp.values
            inm, fl = v_[np.abs(p_) <= 5], v_[(np.abs(p_) > 5)]
            irows.append({"site_id": sid, "ism_motif": inm.mean(), "ism_flank": fl.mean(),
                          "ism_motif_share": inm.sum() / v_.sum() if v_.sum() else np.nan,
                          "ism_peak_position": float(p_[np.argmax(v_)]), "ism_max": v_.max(),
                          "ism_n_strong_bases": int((v_ > 3 * np.median(fl)).sum())})
        F["ism"] = pd.DataFrame(irows).set_index("site_id")
SCF = F["170"].join(F["motif"]).join(F["ism"])
SCF.to_csv(OUT / "tables" / "scan_features_per_site.csv")

# ---------------------------------------------------------------- site table and matching
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
core = pd.read_csv(ROOT / "Data/derived/core_quality.csv").set_index("site_id")
S = sl.join(scl[["age_score_class", "score_dynamic"]]).join(meas.rename("measured")).join(core[["exact_core"]]).join(SCF, how="inner")
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))
S["age_class"] = S.age_class.replace({"Neither YS nor MAS (= MAE only)": np.nan})
S["ex_class"] = S.ex_class.replace({"Neither MAS nor MAE (= YS only)": np.nan})

FEATURES = [
    ("own_motif_170", "effect of replacing the AP-1 motif tile (read over 170 bp)"),
    ("own_motif_motif", "effect of replacing the AP-1 motif tile (read over the motif)"),
    ("best_partner_170", "effect of the strongest other tile within ±300 bp"),
    ("partner_dominates_170", "a neighbouring tile matters more than the motif"),
    ("n_tiles_matter_170", "number of tiles that matter"),
    ("total_dependence_170", "total dependence summed over all tiles"),
    ("gini_170", "concentration of dependence (Gini)"),
    ("spread_170", "spread of dependence (bp)"),
    ("share_within_85_170", "share of dependence within ±85 bp"),
    ("asymmetry_170", "left-right asymmetry of dependence"),
    ("ism_motif", "single-base importance inside the 11 bp motif"),
    ("ism_flank", "single-base importance in the flanks (6–30 bp)"),
    ("ism_motif_share", "share of single-base importance inside the motif"),
    ("ism_max", "largest single-base effect"),
    ("ism_n_strong_bases", "number of strong single bases"),
    ("ism_peak_position", "position of the most important base"),
]
CONTRASTS = [("Dynamic vs Stable", "stable_dynamic", "Dynamic", "Stable"),
             ("Age Up vs Stable (threshold)", "age_class", "Age Up", "Stable"),
             ("Age Down vs Stable (threshold)", "age_class", "Age Down", "Stable"),
             ("Ex Up vs Stable (threshold)", "ex_class", "Exercise Up", "Stable"),
             ("Ex Down vs Stable (threshold)", "ex_class", "Exercise Down", "Stable"),
             ("Age Up vs Unchanged (Δscore)", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down vs Unchanged (Δscore)", "age_score_class", "Age Down (score)", "Unchanged"),
             ("Changing vs Steady (Δscore)", "score_dynamic", "Changing (score)", "Steady (score)")]


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
    bal = max(abs((S.loc[P.case, c].mean() - S.loc[P.control, c].mean()) / S[c].std()) for c in MCOLS) if len(P) else np.nan
    for v, vlab in FEATURES:
        x, y = S[v].reindex(P.case).values.astype(float), S[v].reindex(P.control).values.astype(float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < 25:
            continue
        d = x[ok] - y[ok]
        if not np.any(d != 0):
            continue
        sd = np.nanstd(S[v].astype(float))
        rows.append({"contrast": lab, "feature": vlab, "n_pairs": int(ok.sum()),
                     "median_class": float(np.median(x[ok])), "median_reference": float(np.median(y[ok])),
                     "mean_difference": float(np.mean(d)), "standardised": float(np.mean(d) / sd) if sd else np.nan,
                     "p": float(wilcoxon(d).pvalue), "balance_sd": bal})
R = pd.DataFrame(rows)
R["fdr"] = R.groupby("contrast").p.transform(lambda p: bh(p))
R["detects"] = np.where(R.fdr < 0.05, "YES", "no")
R.to_csv(OUT / "tables" / "scan_features_by_class.csv", index=False)

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(12.2, 8.0))
labels = [f[1] for f in FEATURES]
cs = [c[0] for c in CONTRASTS]
M = R.pivot_table(index="feature", columns="contrast", values="standardised").reindex(index=labels, columns=cs)
Q = R.pivot_table(index="feature", columns="contrast", values="fdr").reindex(index=labels, columns=cs)
m = M.values
lim = np.nanpercentile(np.abs(m), 97) if np.isfinite(m).any() else 1
im = ax.imshow(m, cmap="PuOr_r", vmin=-lim, vmax=lim, aspect="auto")
ax.grid(False)
ax.set_xticks(range(len(cs)), [c.replace(" vs ", "\nvs ").replace(" (threshold)", "").replace(" (\u0394score)", "\n(\u0394score)") for c in cs], fontsize=7.5, rotation=30, ha="right")
ax.set_yticks(range(len(labels)), labels, fontsize=8)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        if not np.isfinite(m[i, j]):
            continue
        sig = Q.values[i, j] < 0.05
        ax.text(j, i, f"{m[i, j]:+.2f}" + ("*" if sig else ""), ha="center", va="center", fontsize=6.6,
                fontweight="bold" if sig else "normal", color="white" if abs(m[i, j]) > 0.42 * lim else INK)
        if sig:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.3))
ax.axvline(4.5, color=INK, lw=1.6)
ax.set_title("What the mutation scan sees, by class\nstandardised difference from the reference class; * and box = FDR < 0.05 within that contrast\n"
             "left of the line: threshold (binary call) classes   |   right: Δscore classes", fontsize=9.5, loc="left")
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("class minus reference (SD units)", fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "y1_scan_features_by_class.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 56)
print("\n== significant results (FDR < 0.05), by contrast")
sig = R[R.fdr < 0.05].sort_values(["contrast", "fdr"])
print(sig[["contrast", "feature", "n_pairs", "median_class", "median_reference", "mean_difference", "standardised", "fdr"]].round(3).to_string(index=False) if len(sig) else "  none")
print("\n== count of detected features per contrast")
print(R.groupby("contrast").detects.apply(lambda v: f"{(v == 'YES').sum()} of {len(v)}").to_string())
print("\n== strongest effect per contrast regardless of significance")
print(R.loc[R.groupby("contrast").standardised.apply(lambda s: s.abs().idxmax())][["contrast", "feature", "mean_difference", "standardised", "fdr"]].round(3).to_string(index=False))
