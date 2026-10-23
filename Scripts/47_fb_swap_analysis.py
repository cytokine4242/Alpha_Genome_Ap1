#!/usr/bin/env python3
"""Analysis of the fibroblast element swaps (Katana/results_ft/swap_fb_chunk*.csv) and host insertions (bg_fb_chunk*.csv).

Question: when an element of one binding class is moved into the location of another, does it bring its own accessibility
with it, or take the level of its new location? Comparison track: the fine-tuned fibroblast head averaged over the four
conditions (trackset FB), read over the inserted element. Everything is in log2 of the predicted value.

  QC        self-swaps (an element put back into its own location) must reproduce the native prediction
  a) sequence vs location: log2(swapped) regressed on log2(donor native) and log2(recipient native), per width;
     the two coefficients say how much of the level the element carries and how much the location imposes
  b) class pairs: median swapped level for every donor class x recipient class, and, for each pair, whether the swapped
     elements sit closer to their donor class's own native level or to their recipient class's (Mann-Whitney between
     the swapped values and each class's natives, BH-adjusted)
  c) host windows: element vs its dinucleotide-shuffled copy in quiet windows, by class: what the sequence does with no context
Outputs Analysis/fb_swap/{tables,figures}.   python3 47_fb_swap_analysis.py      (--simulate makes fake results and runs on them)
"""
from __future__ import annotations

import argparse
import glob
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--simulate", action="store_true")
ap.add_argument("--track", default="FB", help="trackset to analyse (FB = fine-tuned fibroblast head, 4 conditions averaged)")
args = ap.parse_args()
rng = np.random.default_rng(47)
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, AQUA, PLUM, SAND = "#2a78d6", "#eb6834", "#1baf7a", "#8d5fd3", "#c99a2e"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
CLASSES = ["Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down"]
COLOR = dict(zip(CLASSES, [BLUE, ORANGE, AQUA, PLUM, SAND]))
WIDTHS = ["170", "peak"]
WLAB = {"170": "170 bp around the motif", "peak": "whole peak"}


def simulate(folder):
    """Elements carry 60% of their own level and 40% of the location's; host insertions give the element a small lift."""
    folder.mkdir(parents=True, exist_ok=True)
    base = {c: v for c, v in zip(CLASSES, [1.6, 0.8, 0.7, 0.9, 0.85])}
    sites = [(f"s{i}", CLASSES[i % 5], base[CLASSES[i % 5]] * np.exp(rng.normal(0, 0.35))) for i in range(40)]
    rows = []
    for rid, rcls, rnat in sites:
        for did, dcls, dnat in sites:
            for w in WIDTHS:
                v = np.exp(0.6 * np.log(dnat) + 0.4 * np.log(rnat) + rng.normal(0, 0.1))
                if did == rid:
                    v = rnat * np.exp(rng.normal(0, 0.01))
                rows.append({"recipient_id": rid, "recipient_class": rcls, "donor_id": did, "donor_class": dcls, "width": w,
                             "self_swap": did == rid, "insert_bp": 170 if w == "170" else 800, "trackset": args.track,
                             "atac_insert": v, "atac_170": v, "atac_peak": v,
                             "donor_native_atac_insert": dnat, "recipient_native_atac_insert": rnat})
    pd.DataFrame(rows).to_csv(folder / "swap_fb_chunk0.csv", index=False)
    rows = []
    for did, dcls, dnat in sites:
        for h in range(6):
            for w in WIDTHS:
                for kind in ("element", "shuffled"):
                    v = 0.05 * np.exp(rng.normal(0, 0.2)) * (1 + (2.0 * dnat if kind == "element" else 0.1))
                    rows.append({"donor_id": did, "donor_class": dcls, "host_id": f"h{h}", "width": w, "kind": kind,
                                 "insert_bp": 170 if w == "170" else 800, "trackset": args.track, "atac_insert": v,
                                 "donor_native_atac_insert": dnat})
    pd.DataFrame(rows).to_csv(folder / "bg_fb_chunk0.csv", index=False)


if args.simulate:
    tmp = Path(tempfile.mkdtemp())
    simulate(tmp)
    RES, OUT = tmp, tmp / "out"
else:
    RES, OUT = ROOT / "Katana" / "results_ft", ROOT / "Analysis" / "fb_swap"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


def load(prefix):
    fs = sorted(glob.glob(str(RES / f"{prefix}_chunk*.csv")))
    if not fs:
        raise SystemExit(f"no {prefix}_chunk*.csv in {RES}")
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.trackset == args.track].copy()
    if d.empty:
        raise SystemExit(f"{prefix}: trackset {args.track} not in the output (have {sorted(set(pd.concat([pd.read_csv(f, usecols=['trackset']) for f in fs]).trackset))})")
    return d


def lg(x):
    """log2 with a floor, so a zero prediction cannot dominate."""
    return np.log2(np.asarray(x, float).clip(min=1e-3))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


sw = load("swap_fb")
present = [c for c in CLASSES if c in set(sw.donor_class)]

# ---------------------------------------------------------------- QC: self-swaps must reproduce the native prediction
qc = sw[sw.self_swap].copy()
qc["diff"] = lg(qc.atac_insert) - lg(qc.recipient_native_atac_insert)
qc_t = qc.groupby("width")["diff"].agg(n="size", median="median", q05=lambda s: s.quantile(.05), q95=lambda s: s.quantile(.95),
                                       max_abs=lambda s: s.abs().max()).round(4)
qc_t.to_csv(OUT / "tables" / "qc_self_swaps.csv")
print("self-swap check (log2 difference from native; should be ~0):\n", qc_t.to_string())

# ---------------------------------------------------------------- a) sequence versus location
oth = sw[~sw.self_swap].copy()
rows = []
for w in WIDTHS:
    d = oth[oth.width == w]
    if len(d) < 20:
        continue
    y, x1, x2 = lg(d.atac_insert), lg(d.donor_native_atac_insert), lg(d.recipient_native_atac_insert)
    X = np.column_stack([np.ones(len(d)), x1, x2])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    ss = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rows.append({"width": w, "n_swaps": len(d), "intercept": coef[0], "beta_donor_sequence": coef[1],
                 "beta_recipient_location": coef[2], "r2": ss,
                 "spearman_vs_donor": d.atac_insert.corr(d.donor_native_atac_insert, method="spearman"),
                 "spearman_vs_recipient": d.atac_insert.corr(d.recipient_native_atac_insert, method="spearman")})
SEQLOC = pd.DataFrame(rows)
SEQLOC.to_csv(OUT / "tables" / "sequence_vs_location.csv", index=False)
print("\nsequence vs location:\n", SEQLOC.round(3).to_string(index=False))

fig, axes = plt.subplots(1, len(WIDTHS), figsize=(5.3 * len(WIDTHS), 4.2))
for ax, w in zip(np.atleast_1d(axes), WIDTHS):
    d = oth[oth.width == w]
    if d.empty:
        continue
    ax.scatter(lg(d.donor_native_atac_insert), lg(d.atac_insert), s=4, alpha=0.15, color=BLUE, label="vs the element's own native level", rasterized=True)
    ax.scatter(lg(d.recipient_native_atac_insert), lg(d.atac_insert), s=4, alpha=0.15, color=ORANGE, label="vs the location's native level", rasterized=True)
    lim = [min(lg(d.donor_native_atac_insert).min(), lg(d.atac_insert).min()), max(lg(d.donor_native_atac_insert).max(), lg(d.atac_insert).max())]
    ax.plot(lim, lim, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("Native predicted accessibility (log2)")
    ax.set_ylabel("Swapped element's accessibility (log2)")
    ax.set_title(f"Replace the {WLAB[w]}", fontsize=10, loc="left")
    leg = ax.legend(frameon=False, fontsize=8, markerscale=3)
    for h in leg.legend_handles:
        h.set_alpha(1)
fig.tight_layout()
fig.savefig(OUT / "figures" / "w1_sequence_vs_location.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- b) class pairs
nat = pd.concat([sw[["donor_id", "donor_class", "width", "donor_native_atac_insert"]]
                 .rename(columns={"donor_id": "site_id", "donor_class": "swap_class", "donor_native_atac_insert": "native"})]).drop_duplicates(["site_id", "width"])
nat["lnat"] = lg(nat.native)
rows, pair_rows = [], []
for w in WIDTHS:
    d = oth[oth.width == w]
    nw = nat[nat.width == w]
    if d.empty or nw.empty:
        continue
    for dc in present:
        for rc in present:
            v = lg(d[(d.donor_class == dc) & (d.recipient_class == rc)].atac_insert)
            if len(v) < 5:
                continue
            dn = nw[nw.swap_class == dc].lnat
            rn = nw[nw.swap_class == rc].lnat
            pd_d = mannwhitneyu(v, dn).pvalue if len(dn) > 2 else np.nan
            pd_r = mannwhitneyu(v, rn).pvalue if len(rn) > 2 else np.nan
            pair_rows.append({"width": w, "donor_class": dc, "recipient_class": rc, "n": len(v), "median_swapped": v.median(),
                              "median_donor_class_native": dn.median(), "median_recipient_class_native": rn.median(),
                              "gap_to_donor_class": v.median() - dn.median(), "gap_to_recipient_class": v.median() - rn.median(),
                              "p_vs_donor_class": pd_d, "p_vs_recipient_class": pd_r})
            rows.append({"width": w, "donor_class": dc, "recipient_class": rc, "median_swapped": v.median()})
PAIRS = pd.DataFrame(pair_rows)
if not PAIRS.empty:
    for c in ("p_vs_donor_class", "p_vs_recipient_class"):
        PAIRS[c.replace("p_", "fdr_")] = bh(PAIRS[c].fillna(1).values)
    PAIRS["closer_to"] = np.where(PAIRS.gap_to_donor_class.abs() < PAIRS.gap_to_recipient_class.abs(), "its own class", "the location's class")
PAIRS.to_csv(OUT / "tables" / "class_pairs.csv", index=False)

if rows:
    fig, axes = plt.subplots(1, len(WIDTHS), figsize=(5.6 * len(WIDTHS), 4.6))
    grids = []
    for ax, w in zip(np.atleast_1d(axes), WIDTHS):
        g = pd.DataFrame(rows)
        g = g[g.width == w].pivot(index="donor_class", columns="recipient_class", values="median_swapped").reindex(index=present, columns=present)
        grids.append(g)
        im = ax.imshow(g.values, cmap="viridis")
        ax.set_xticks(range(len(present)), present, rotation=35, ha="right")
        ax.set_yticks(range(len(present)), present)
        ax.set_xlabel("Location's class (recipient)")
        ax.set_ylabel("Element's class (donor)")
        ax.set_title(f"Median accessibility after replacing the {WLAB[w]} (log2)", fontsize=9.5, loc="left")
        ax.grid(False)
        for i in range(g.shape[0]):
            for j in range(g.shape[1]):
                if np.isfinite(g.values[i, j]):
                    ax.text(j, i, f"{g.values[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                            color="white" if g.values[i, j] < np.nanmedian(g.values) else "black")
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "w2_class_pair_grid.png", bbox_inches="tight")
    plt.close(fig)
    pd.concat({w: g for w, g in zip(WIDTHS, grids)}).to_csv(OUT / "tables" / "class_pair_grid.csv")

# ---------------------------------------------------------------- c) host windows
try:
    bg = load("bg_fb")
except SystemExit as e:
    print("\nhost insertions:", e)
    bg = None
if bg is not None:
    piv = bg.pivot_table(index=["donor_id", "donor_class", "host_id", "width"], columns="kind", values="atac_insert").reset_index()
    if {"element", "shuffled"} <= set(piv.columns):
        piv["lift"] = lg(piv.element) - lg(piv.shuffled)
        BG = piv.groupby(["width", "donor_class"]).lift.agg(n="size", median="median",
                                                            q25=lambda s: s.quantile(.25), q75=lambda s: s.quantile(.75)).reset_index()
        rows = []
        for w in WIDTHS:
            sub = piv[piv.width == w]
            for i, a in enumerate(present):
                for b in present[i + 1:]:
                    x, y = sub[sub.donor_class == a].lift, sub[sub.donor_class == b].lift
                    if len(x) > 5 and len(y) > 5:
                        rows.append({"width": w, "class_a": a, "class_b": b, "median_a": x.median(), "median_b": y.median(),
                                     "p": mannwhitneyu(x, y).pvalue})
        T = pd.DataFrame(rows)
        if not T.empty:
            T["fdr"] = bh(T.p.values)
        T.to_csv(OUT / "tables" / "host_insertion_class_tests.csv", index=False)
        BG.to_csv(OUT / "tables" / "host_insertion_lift.csv", index=False)
        print("\nhost windows, lift of the element over its shuffled copy (log2):\n", BG.round(3).to_string(index=False))

        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        for k, w in enumerate(WIDTHS):
            for i, c in enumerate(present):
                v = piv[(piv.width == w) & (piv.donor_class == c)].lift
                if len(v):
                    ax.errorbar(i + (k - 0.5) * 0.26, v.median(), yerr=[[v.median() - v.quantile(.25)], [v.quantile(.75) - v.median()]],
                                fmt="o" if w == "170" else "s", color=COLOR[c], ms=6, lw=1.4, capsize=0,
                                label=WLAB[w] if i == 0 else None)
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(range(len(present)), present, rotation=20, ha="right")
        ax.set_ylabel("Accessibility over the shuffled copy (log2)")
        ax.set_title("The element on its own, inserted into quiet windows", fontsize=10, loc="left")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "figures" / "w3_host_insertions.png", bbox_inches="tight")
        plt.close(fig)

print("\noutputs in", OUT)
