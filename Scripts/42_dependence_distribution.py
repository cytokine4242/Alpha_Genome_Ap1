#!/usr/bin/env python3
"""Two questions the earlier work left open.
 A  SEQUENCE motif counts (Data/derived/motif_counts.csv, a consensus scan of mm10 that needs no ATAC signal) by class.
    This separates 'fewer MEF2 motifs' (sequence, real) from 'fewer MEF2 footprints called' (detection), which footprint data alone cannot do.
 B  How is the dependence DISTRIBUTED across the peak? From the tile scan: how many tiles matter, how concentrated they are, how far they spread,
    and where they sit relative to the peak rather than to the motif. Shape measures are normalised per site, so a class that loses less overall is not penalised.
Outputs Analysis/dependence/{tables,figures}.   python3 42_dependence_distribution.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "dependence"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
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
mc = pd.read_csv(ROOT / "Data/derived/motif_counts.csv").set_index("site_id")
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id")
S = sl.join(scl[["age_score_class", "ex_score_class"]]).join(meas.rename("measured")).join(mc).join(sites[["motif_mid", "peak_start", "peak_end"]])
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))

CONTRASTS = [("Dynamic vs Stable", "stable_dynamic", "Dynamic", "Stable"),
             ("Age Up (score) vs Unchanged", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down (score) vs Unchanged", "age_score_class", "Age Down (score)", "Unchanged")]


def make_pairs(frame, col, a_, b_, cols=MCOLS, caliper=0.5):
    m_ = frame[frame[col].isin([a_, b_])].dropna(subset=cols).copy()
    Z = (m_[cols] - m_[cols].mean()) / m_[cols].std()
    avail, pairs = set(m_.index[m_[col] == b_]), []
    for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < caliper:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


def paired(frame, col, a_, b_, measures, cols=MCOLS):
    P = make_pairs(frame, col, a_, b_, cols)
    out = []
    for v, vlab in measures:
        x, y = frame[v].reindex(P.case).values.astype(float), frame[v].reindex(P.control).values.astype(float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < 25 or not np.any(x[ok] != y[ok]):
            continue
        d = x[ok] - y[ok]
        out.append({"measure": vlab, "n_pairs": int(ok.sum()), "median_case": float(np.median(x[ok])), "median_control": float(np.median(y[ok])),
                    "paired_median_diff": float(np.median(d)), "mean_diff": float(np.mean(d)), "p": float(wilcoxon(d).pvalue)})
    return pd.DataFrame(out), P


# ---------------------------------------------------------------- A. sequence motif counts, matched
SEQM = [(f"{m}_{w}", f"{m} consensus matches within {w} bp") for w in ("100", "500") for m in ("MEF2", "AP1", "GATA", "TEAD", "NFI", "CTCF", "NKX", "KLF")] + [("GC_200", "GC content, central 200 bp")]
rows = []
for lab, col, a_, b_ in CONTRASTS:
    r, P = paired(S, col, a_, b_, SEQM)
    r.insert(0, "comparison", lab)
    r["balance_sd"] = max(abs((S.loc[P.case, c].mean() - S.loc[P.control, c].mean()) / S[c].std()) for c in MCOLS)
    # unmatched version, for comparison
    raw = []
    for v, vlab in SEQM:
        x, y = S[S[col] == a_][v].dropna(), S[S[col] == b_][v].dropna()
        raw.append({"measure": vlab, "raw_median_case": x.median(), "raw_median_control": y.median(), "p_raw": mannwhitneyu(x, y).pvalue})
    rows.append(r.merge(pd.DataFrame(raw), on="measure"))
A = pd.concat(rows, ignore_index=True)
A["fdr"] = A.groupby("comparison").p.transform(lambda p: bh(p))
A["fdr_raw"] = A.groupby("comparison").p_raw.transform(lambda p: bh(p))
A.to_csv(OUT / "tables" / "a_sequence_motif_counts_by_class.csv", index=False)

# ---------------------------------------------------------------- B. how dependence is spread across the peak
fs = sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))
B = D = None
if fs:
    sc = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    sc = sc[sc.trackset == "CM"]
    n_pass = sc.groupby("site_id").size()
    sc = sc[sc.site_id.isin(n_pass.index[n_pass == n_pass.max()])]
    nat = sc[sc.kind == "native"].set_index("site_id").atac_170
    pc = 0.05 * nat.median()
    t = sc[sc.kind == "tile"].copy()
    t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
    tile = t.groupby(["site_id", "group", "pos"], as_index=False).eff.mean()
    tile["centre"] = tile.pos + 5
    far = tile[tile.centre.abs() >= 250].eff
    THR = float(np.percentile(far, 1))                       # 'matters' = beyond the 1st percentile of tiles 250-300 bp away
    rows = []
    for (sid, g), w in tile.groupby(["site_id", "group"]):
        w = w.sort_values("centre")
        loss = np.clip(-w.eff.values, 0, None)               # dependence mass per tile
        tot = loss.sum()
        if tot <= 0:
            continue
        share = loss / tot
        order = np.argsort(-loss)
        cum = np.cumsum(loss[order]) / tot
        n50 = int(np.searchsorted(cum, 0.5) + 1)             # tiles needed for half the dependence (shape only)
        srt = np.sort(loss)
        gini = float((2 * np.arange(1, len(srt) + 1) - len(srt) - 1).dot(srt) / (len(srt) * srt.sum()))
        c = w.centre.values
        centroid = float((share * c).sum())
        spread = float(np.sqrt((share * (c - centroid) ** 2).sum()))
        strong = w.eff.values < THR
        runs = int(np.sum(np.diff(np.r_[0, strong.astype(int)]) == 1))
        pk_half = max(S.peak_width.get(sid, np.nan) / 2, 1)
        off = (S.motif_mid.get(sid, np.nan) - (S.peak_start.get(sid, np.nan) + S.peak_end.get(sid, np.nan)) / 2)
        rows.append({"site_id": sid, "group": g, "n_tiles_matter": int(strong.sum()), "n_runs_of_strong_tiles": runs,
                     "total_dependence": float(tot), "tiles_for_half_the_dependence": n50, "gini_concentration": gini,
                     "spread_bp": spread, "centroid_bp": centroid, "abs_centroid_bp": abs(centroid),
                     "share_within_85bp": float(share[np.abs(c) <= 85].sum()), "share_beyond_150bp": float(share[np.abs(c) > 150].sum()),
                     "share_inside_peak": float(share[np.abs(c + off) <= pk_half].sum()) if np.isfinite(off) else np.nan,
                     "centroid_in_peak_halfwidths": float((centroid + off) / pk_half) if np.isfinite(off) else np.nan})
    D = pd.DataFrame(rows).set_index("site_id")
    D.to_csv(OUT / "tables" / "b_dependence_shape_per_site.csv")
    DM = [("n_tiles_matter", "tiles that matter (vs the far-tile threshold)"), ("n_runs_of_strong_tiles", "separate runs of tiles that matter"),
          ("total_dependence", "total dependence (sum of losses)"), ("tiles_for_half_the_dependence", "tiles holding half the dependence"),
          ("gini_concentration", "concentration of dependence (Gini)"), ("spread_bp", "spread of dependence (bp)"),
          ("abs_centroid_bp", "offset of the dependence centroid from the motif (bp)"), ("share_within_85bp", "share of dependence within 85 bp"),
          ("share_beyond_150bp", "share of dependence beyond 150 bp"), ("share_inside_peak", "share of dependence inside the peak"),
          ("centroid_in_peak_halfwidths", "centroid position in peak half-widths")]
    SD = S.join(D[[m[0] for m in DM] + ["group"]], how="inner")
    rows = []
    for lab, col, a_, b_ in CONTRASTS:
        r, _ = paired(SD, col, a_, b_, DM)
        if len(r):
            r.insert(0, "comparison", lab)
            rows.append(r)
    for a_, b_ in (("dynamic", "matched_stable"), ("dynamic", "high_stable")):
        out = []
        for v, vlab in DM:
            x, y = SD[SD.group == a_][v].dropna(), SD[SD.group == b_][v].dropna()
            if min(len(x), len(y)) < 25:
                continue
            out.append({"comparison": f"scan groups: {a_} vs {b_}", "measure": vlab, "n_pairs": np.nan, "median_case": x.median(), "median_control": y.median(),
                        "paired_median_diff": np.nan, "mean_diff": np.nan, "p": mannwhitneyu(x, y).pvalue})
        rows.append(pd.DataFrame(out))
    B = pd.concat(rows, ignore_index=True)
    B["fdr"] = B.groupby("comparison").p.transform(lambda p: bh(p))
    B.to_csv(OUT / "tables" / "b_dependence_shape_by_class.csv", index=False)

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.2))
ax = axes[0]
w = A[A.comparison == "Dynamic vs Stable"].set_index("measure")
sel = [m[1] for m in SEQM if m[1].endswith("100 bp")]
ys = np.arange(len(sel))[::-1]
v = w.reindex(sel)
ax.scatter(v.median_control, ys, s=44, color=BLUE, label="Stable (matched)", zorder=3)
ax.scatter(v.median_case, ys, s=44, color=ORANGE, label="Dynamic", zorder=3)
for y_, a_, b_, q in zip(ys, v.median_case, v.median_control, v.fdr):
    ax.plot([a_, b_], [y_, y_], color=GREY, lw=1.1, zorder=2)
    if q < 0.05:
        ax.text(max(a_, b_) + 0.03, y_, "FDR < 0.05", fontsize=7.5, va="center")
ax.set_yticks(ys, [s.replace(" consensus matches within 100 bp", "") for s in sel], fontsize=8.5)
ax.set_xlabel("Median consensus matches within 100 bp (sequence scan, no ATAC signal used)")
ax.set_title("Motifs in the SEQUENCE, matched pairs", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
if D is not None:
    ax = axes[1]
    for g, lab, c in (("dynamic", "Dynamic", ORANGE), ("matched_stable", "Stable, context-matched", BLUE), ("high_stable", "Stable, high accessibility", AQUA)):
        v = D[D.group == g].tiles_for_half_the_dependence.dropna()
        if len(v):
            h = v.value_counts(normalize=True).sort_index()
            ax.plot(h.index[:20], h.values[:20], color=c, lw=1.8, marker="o", ms=3.5, label=f"{lab} (median {v.median():.0f})")
    ax.set_xlabel("Tiles holding half of the site's total dependence")
    ax.set_ylabel("Share of sites")
    ax.set_title("Is dependence concentrated or spread?", fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "u1_sequence_and_shape.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
print("== A sequence motif counts, matched pairs\n", A[A.fdr < 0.2][["comparison", "measure", "n_pairs", "median_case", "median_control", "paired_median_diff", "p", "fdr", "fdr_raw"]].round(3).to_string())
print("\nMEF2 within 100 bp, all three contrasts:\n", A[A.measure.str.startswith("MEF2 consensus matches within 100")][["comparison", "median_case", "median_control", "mean_diff", "p", "fdr"]].round(3).to_string())
if B is not None:
    print("== B dependence shape by class\n", B.round(3).to_string())
