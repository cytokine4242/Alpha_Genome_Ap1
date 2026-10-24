#!/usr/bin/env python3
"""Fibroblast results read with the fine-tuned 1 bp ATAC head (Katana/results_ft/shuffle_fb_chunk*.csv).

Effect = per-site mean over replicates of log2((edit + pc) / (native + pc)), pc = 5% of the median native value of that
readout in that trackset (so the two heads, which live on different scales, are each floored on their own scale).

  1  shuffle effects by binding class (stable_dynamic, then age_class and ex_class as given in the site table)
  2  the same sites and the same edited sequences read with the fine-tuned head and with the built-in Patski head from the
     SAME forward pass (trackset FB vs FB_builtin): does a cell-type-specific 1 bp head change the conclusions?
  3  per-condition effects (YS / MAS / MAE / YE channels of the fine-tuned head): is the AP-1 motif worth more in one
     condition than another? Read as direction and rank only - the evaluation showed predicted condition differences are
     about a third the measured size.
Outputs Analysis/fb_finetuned/{tables,figures}.   python3 48_ft_fibroblast_results.py [--res DIR] [--simulate]
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
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None, help="results_ft directory (default Katana/results_ft)")
ap.add_argument("--simulate", action="store_true")
args = ap.parse_args()
INK, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BLUE, ORANGE, AQUA, PLUM, SAND = "#2a78d6", "#eb6834", "#1baf7a", "#8d5fd3", "#c99a2e"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
CONDS = ["YS", "MAS", "MAE", "YE"]
WIDTH_LABEL = {"motif": "the 11 bp motif", "170": "the 170 bp around it", "peak": "the whole peak"}


def simulate(folder):
    """Dynamic sites lose less to wide edits; the fine-tuned head separates the classes more than the built-in one."""
    rng = np.random.default_rng(48)
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(120):
        sd = "Dynamic" if i % 3 == 0 else "Stable"
        ac = ["Stable", "Age Up", "Age Down"][i % 3]
        for ts, base, sep in (("FB", 3.0, 1.0), ("FB_builtin", 0.3, 0.4)):
            nat = base * np.exp(rng.normal(0, 0.4))
            common = {"site_id": f"s{i}", "binding_category": "x", "age_class": ac, "ex_class": "Stable",
                      "stable_dynamic": sd, "trackset": ts}
            rows.append({**common, "condition": "native", "kind": "native", "width": "", "rep": 0,
                         **{f"atac_{w}": nat for w in ("motif", "170", "peak")},
                         **{f"atac_{w}_{c}": nat for w in ("motif", "170", "peak") for c in CONDS}})
            for w, loss in (("motif", 0.15), ("170", 0.55), ("peak", 0.8)):
                for r in range(10):
                    f = 1 - loss * sep * (0.5 if sd == "Dynamic" else 1.0) * np.exp(rng.normal(0, 0.15))
                    rows.append({**common, "condition": f"shuffle@{w}", "kind": "shuffle", "width": w, "rep": r,
                                 **{f"atac_{x}": nat * (f if x == w else 1 - (1 - f) * 0.6) for x in ("motif", "170", "peak")},
                                 **{f"atac_{x}_{c}": nat * (f if x == w else 1 - (1 - f) * 0.6) *
                                    (0.9 if c == "MAS" and w == "motif" else 1.0) for x in ("motif", "170", "peak") for c in CONDS}})
    pd.DataFrame(rows).to_csv(folder / "shuffle_fb_chunk0.csv", index=False)


if args.simulate:
    tmp = Path(tempfile.mkdtemp())
    simulate(tmp)
    RES, OUT = tmp, tmp / "out"
else:
    RES = Path(args.res) if args.res else ROOT / "Katana" / "results_ft"
    OUT = ROOT / "Analysis" / "fb_finetuned"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def effects(d, col):
    """Per site x condition: mean over replicates of log2((edit + pc) / (native + pc)), within one trackset."""
    nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[d.kind != "native"].copy()
    e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "stable_dynamic", "age_class", "ex_class", "condition", "width"], as_index=False).eff.mean()


files = sorted(glob.glob(str(RES / "shuffle_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no shuffle_fb_chunk*.csv in {RES}")
keep = ["site_id", "age_class", "ex_class", "stable_dynamic", "condition", "kind", "width", "rep", "trackset",
        "atac_motif", "atac_170", "atac_peak"] + [f"atac_{w}_{c}" for w in ("motif", "170", "peak") for c in CONDS]
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in keep) for f in files], ignore_index=True)
d["width"] = d.width.fillna("")
print(f"{len(files)} chunks, {d.site_id.nunique()} sites, tracksets {sorted(set(d.trackset))}")

# ---------------------------------------------------------------- 1. shuffle effects by class, fine-tuned head
ft = d[d.trackset == "FB"]
rows, tests = [], []
for col, wname in (("atac_motif", "motif"), ("atac_170", "170"), ("atac_peak", "peak")):
    e = effects(ft, col)
    for w in ("motif", "170", "peak"):
        s = e[e.width == w]
        if s.empty:
            continue
        rows.append({"readout": wname, "edit": f"shuffle {WIDTH_LABEL[w]}", "n_sites": s.site_id.nunique(),
                     "median_all": s.eff.median(),
                     **{f"median_{k}": s[s.stable_dynamic == k].eff.median() for k in ("Stable", "Dynamic")}})
        if wname == w:                                   # the readout that matches the edit is the one to test on
            for field, groups in (("stable_dynamic", ("Stable", "Dynamic")),
                                  ("age_class", ("Stable", "Age Up")), ("age_class", ("Stable", "Age Down")),
                                  ("ex_class", ("Stable", "Exercise Up")), ("ex_class", ("Stable", "Exercise Down"))):
                x = s[s[field] == groups[0]].eff
                y = s[s[field] == groups[1]].eff
                if len(x) > 5 and len(y) > 5:
                    tests.append({"edit": f"shuffle {WIDTH_LABEL[w]}", "contrast": f"{groups[0]} vs {groups[1]} ({field})",
                                  "n_a": len(x), "n_b": len(y), "median_a": x.median(), "median_b": y.median(),
                                  "difference": y.median() - x.median(), "p": mannwhitneyu(x, y).pvalue})
BYCLASS = pd.DataFrame(rows)
T1 = pd.DataFrame(tests)
if not T1.empty:
    T1["fdr"] = bh(T1.p.values)
BYCLASS.to_csv(OUT / "tables" / "shuffle_effects_by_class.csv", index=False)
T1.to_csv(OUT / "tables" / "shuffle_class_tests.csv", index=False)
print("\n=== 1. shuffle effects, fine-tuned fibroblast head (log2, matched readout)")
print(BYCLASS.round(3).to_string(index=False))
print("\nclass contrasts:")
print(T1.round(4).to_string(index=False) if not T1.empty else "  (not enough sites)")

# ---------------------------------------------------------------- 2. fine-tuned vs built-in head, same forward passes
cmp_rows = []
for col, w in (("atac_motif", "motif"), ("atac_170", "170"), ("atac_peak", "peak")):
    for ts, lab in (("FB", "fine-tuned"), ("FB_builtin", "built-in Patski")):
        s = d[d.trackset == ts]
        if s.empty:
            continue
        e = effects(s, col)
        e = e[e.width == w]
        if e.empty:
            continue
        x = e[e.stable_dynamic == "Stable"].eff
        y = e[e.stable_dynamic == "Dynamic"].eff
        cmp_rows.append({"edit": f"shuffle {WIDTH_LABEL[w]}", "head": lab, "n_sites": e.site_id.nunique(),
                         "median_effect": e.eff.median(), "median_Stable": x.median(), "median_Dynamic": y.median(),
                         "Dynamic - Stable": y.median() - x.median(),
                         "p_Stable_vs_Dynamic": mannwhitneyu(x, y).pvalue if len(x) > 5 and len(y) > 5 else np.nan})
HEADS = pd.DataFrame(cmp_rows)
if not HEADS.empty:
    HEADS["fdr"] = bh(HEADS.p_Stable_vs_Dynamic.fillna(1).values)
HEADS.to_csv(OUT / "tables" / "finetuned_vs_builtin_head.csv", index=False)
print("\n=== 2. the same edits read with each head (identical forward passes)")
print(HEADS.round(4).to_string(index=False))

if not HEADS.empty:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    widths = ["motif", "170", "peak"]
    for ax, (ts, lab, c) in zip(axes, (("FB", "Fine-tuned fibroblast head", ORANGE), ("FB_builtin", "Built-in Patski head", BLUE))):
        s = d[d.trackset == ts]
        for k, (sd, mk) in enumerate((("Stable", "o"), ("Dynamic", "s"))):
            med, lo, hi = [], [], []
            for w, col in zip(widths, ("atac_motif", "atac_170", "atac_peak")):
                e = effects(s, col)
                v = e[(e.width == w) & (e.stable_dynamic == sd)].eff
                med.append(v.median()); lo.append(v.quantile(.25)); hi.append(v.quantile(.75))
            xs = np.arange(len(widths)) + (k - 0.5) * 0.16
            ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt=mk, color=c, ms=6, lw=1.4,
                        capsize=0, alpha=1 if sd == "Stable" else 0.6, label=sd)
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(range(len(widths)), [WIDTH_LABEL[w] for w in widths], rotation=12)
        ax.set_ylabel("Effect of shuffling (log2)")
        ax.set_title(lab, fontsize=10, loc="left")
        ax.legend(frameon=False, fontsize=8, title="binding")
    ylim = (min(a.get_ylim()[0] for a in axes), max(a.get_ylim()[1] for a in axes))
    for a in axes:
        a.set_ylim(ylim)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "f1_head_comparison.png", bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- 3. per-condition effects
cond_rows, cond_tests = [], []
for w, base in (("motif", "atac_motif"), ("170", "atac_170"), ("peak", "atac_peak")):
    per = {}
    for c in CONDS:
        col = f"{base}_{c}"
        if col not in ft.columns or ft[col].isna().all():
            continue
        e = effects(ft, col)
        e = e[e.width == w]
        per[c] = e.set_index("site_id").eff
        cond_rows.append({"edit": f"shuffle {WIDTH_LABEL[w]}", "condition": c, "n_sites": len(e),
                          "median_effect": e.eff.median()})
    for a, b in (("YS", "MAS"), ("MAS", "MAE"), ("YS", "YE")):     # ageing, exercise in aged, exercise in young
        if a in per and b in per:
            j = pd.concat([per[a].rename("a"), per[b].rename("b")], axis=1).dropna()
            diff = j.b - j.a
            if len(j) > 10:
                # two channels can be bit-identical (nothing learned apart for them); Wilcoxon is undefined on all-zero differences
                pv = wilcoxon(j.a, j.b).pvalue if (diff != 0).any() else np.nan
                cond_tests.append({"edit": f"shuffle {WIDTH_LABEL[w]}", "contrast": f"{a} vs {b}", "n_sites": len(j),
                                   "median_a": j.a.median(), "median_b": j.b.median(),
                                   "median_paired_difference": diff.median(), "n_nonzero_differences": int((diff != 0).sum()),
                                   "p_wilcoxon": pv})
COND = pd.DataFrame(cond_rows)
T3 = pd.DataFrame(cond_tests)
if not T3.empty:
    T3["fdr"] = bh(T3.p_wilcoxon.fillna(1).values)
COND.to_csv(OUT / "tables" / "per_condition_effects.csv", index=False)
T3.to_csv(OUT / "tables" / "per_condition_tests.csv", index=False)
print("\n=== 3. per-condition effects (fine-tuned head channels; direction and rank only)")
print(COND.round(4).to_string(index=False))
print("\npaired contrasts between conditions (same sites):")
print(T3.round(5).to_string(index=False) if not T3.empty else "  (not enough sites)")
print("\noutputs in", OUT)
