#!/usr/bin/env python3
"""Analysis of the tile-scan / saturation-mutagenesis and AP-1 titration experiments (Katana/results/scan_cm_chunk*.csv, titrate_cm_chunk*.csv), cardiomyocyte tracks only.
Groups (from 34_make_dynamics_subsets.py): dynamic, matched_stable (matched on predicted accessibility, peak width, TSS distance), high_stable (top-quartile predicted accessibility).
Effect = per-site mean over replicates of log2((edit + pc) / (native + pc)), pc = 5% of the median native value of the readout.
  scan:     tile profile (effect on the motif and 170 bp readouts by tile position), saturation-mutagenesis importance by position, share of the importance that sits in the 11 bp motif
  titrate:  dose-response in the number of core bases changed, consensus upgrade, reverse complement, half-site, and the net effect of duplicating / moving the motif (minus the shuffled control)
Outputs Analysis/scan_titrate/{tables,figures}.   python3 35_scan_titrate_analysis.py            (--simulate makes fake results in a temp folder and runs on them to test the code)"""
from __future__ import annotations

import argparse
import glob
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--simulate", action="store_true")
args = ap.parse_args()
rng = np.random.default_rng(35)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
GROUPS = [("dynamic", "Dynamic", ORANGE), ("matched_stable", "Stable, context-matched", BLUE), ("high_stable", "Stable, high accessibility", AQUA)]


def simulate(folder):
    folder.mkdir(parents=True, exist_ok=True)
    rows_s, rows_t = [], []
    for g, _, _ in GROUPS:
        for i in range(40):
            sid, base = f"{g}-{i}", {"dynamic": 0.25, "matched_stable": 0.27, "high_stable": 1.4}[g] * np.exp(rng.normal(0, 0.3))
            common = {"site_id": sid, "group": g, "binding_category": "x", "stable_dynamic": "Dynamic" if g == "dynamic" else "Stable", "native_170_pred": base, "native_core": "TGAGTCA", "core_found": bool(rng.random() > (0.3 if g == "dynamic" else 0.1)),
                      "trackset": "CM", "tss_dist": 1000, "atac_peak": base, "rna_mean50kb": 1, "rna_max50kb": 1, "rna_prox2kb": 1}
            scale = 0.5 if g == "dynamic" else 1.0

            def row(cond, kind, pos, alt, rep, f, tgt=np.nan):
                return {**common, "condition": cond, "kind": kind, "pos": pos, "alt": alt, "rep": rep, "atac_motif": base * f, "atac_170": base * (1 - (1 - f) * 0.7), "atac_target": tgt}
            rows_s.append(row("native", "native", 0, "", 0, 1.0))
            for o in range(-300, 300, 10):
                for r in range(3):
                    rows_s.append(row(f"tile@{o:+d}", "tile", o, "", r, 1 - scale * 0.6 * np.exp(-abs(o + 5) / 40) + rng.normal(0, 0.03)))
            for p in range(-30, 31):
                for alt in "ACGT"[:3]:
                    rows_s.append(row(f"ism@{p:+d}{alt}", "ism", p, alt, 0, 1 - scale * (0.3 if abs(p) <= 5 else 0.05) * np.exp(rng.normal(0, 0.3))))
            rows_t.append(row("native", "native", 0, "", 0, 1.0))
            for k in range(1, 8):
                for r in range(5 if k < 7 else 3):
                    rows_t.append(row(f"core_mm{k}", "core_mm", k, "", r, 1 - scale * 0.07 * k + rng.normal(0, 0.03)))
            rows_t.append(row("core_consensus", "core_consensus", 0, "", 0, 1 + 0.2 * (1 if g == "dynamic" else 0.1)))
            rows_t.append(row("core_rc", "core_rc", 0, "", 0, 1.0 + rng.normal(0, 0.02)))
            for r in range(3):
                rows_t.append(row("core_halfsite", "core_halfsite", 0, "", r, 0.8))
            for o in (-150, -80, -40, -20, 20, 40, 80, 150):
                rows_t.append(row(f"dup@{o:+d}", "dup", o, "", 0, 1.05, base * 1.3))
                for r in range(3):
                    rows_t.append(row(f"ctrl@{o:+d}", "ctrl", o, "", r, 0.99, base))
                    rows_t.append(row(f"move@{o:+d}", "move", o, "", r, 0.6, base * 1.2))
                    rows_t.append(row(f"ctrlmove@{o:+d}", "ctrlmove", o, "", r, 0.55, base * 0.98))
    pd.DataFrame(rows_s).to_csv(folder / "scan_cm_chunk0.csv", index=False)
    pd.DataFrame(rows_t).to_csv(folder / "titrate_cm_chunk0.csv", index=False)


if args.simulate:
    tmp = Path(tempfile.mkdtemp())
    simulate(tmp)
    RES, OUT = tmp, tmp / "out"
else:
    RES, OUT = ROOT / "Katana" / "results", ROOT / "Analysis" / "scan_titrate"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)


def load(prefix, require_core=False):
    fs = sorted(glob.glob(str(RES / f"{prefix}_chunk*.csv")))
    if not fs:
        raise SystemExit(f"no {prefix}_chunk*.csv in {RES}")
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.trackset == "CM"]
    n = d.groupby("site_id").size()
    short = n[n < 0.9 * n.median()]
    if len(short):
        print(f"WARNING {prefix}: {len(short)} sites have fewer rows than expected (median {int(n.median())}); dropped as incomplete")
        d = d[~d.site_id.isin(short.index)]
    cf = d.drop_duplicates("site_id").groupby("group").core_found.mean()
    print(f"{prefix}: sites per group {d.drop_duplicates('site_id').group.value_counts().to_dict()}; share with an exact TGA[C/G]TCA core: {cf.round(3).to_dict()}")
    pd.DataFrame({"share_with_exact_core": cf}).to_csv(OUT / "tables" / f"{prefix}_core_found_by_group.csv")
    if require_core:
        d = d[d.core_found]
        print(f"{prefix}: kept {d.site_id.nunique()} sites with an exact core")
    return d


def effects(d, col):
    nat = d[d.condition == "native"].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = d[d.condition != "native"].copy()
    e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    return e.groupby(["site_id", "group", "condition", "kind", "pos"], as_index=False).eff.mean()


def boot(v, n=200):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5]) if len(v) > 5 else (np.nan, np.nan)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ============================================================ scan
sc = load("scan_cm")
RDS = {"motif": "atac_motif", "170 bp": "atac_170"}
fig, axes = plt.subplots(1, 2, figsize=(13, 3.9), sharey=False)
tile_rows = []
for ax, (rd, col) in zip(axes, RDS.items()):
    e = effects(sc, col)
    t = e[e.kind == "tile"]
    for g, lab, c in GROUPS:
        x = t[t.group == g].groupby("pos").eff
        med, ps = x.median(), x.median().index.values
        lo = x.apply(lambda v: boot(v.values)[0]).values
        hi = x.apply(lambda v: boot(v.values)[1]).values
        ax.plot(ps + 4.5, med.values, color=c, lw=1.8, label=lab)
        ax.fill_between(ps + 4.5, lo, hi, color=c, alpha=0.15, lw=0)
        tile_rows += [{"readout": rd, "group": g, "tile_centre": int(p_) + 5, "median_effect": m_} for p_, m_ in zip(ps, med.values)]
    ax.axvspan(-5.5, 5.5, color=GRID, alpha=0.6, lw=0)
    ax.set_xlabel("Centre of the replaced 10 bp tile, relative to the motif midpoint (bp)")
    ax.set_ylabel(f"Median effect on accessibility over the {rd} (log2FC)")
    ax.set_title(f"Tile scan, read over the {rd}", fontsize=10, loc="left")
axes[0].legend(frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "figures" / "s1_tile_profile.png", bbox_inches="tight")
plt.close(fig)
pd.DataFrame(tile_rows).to_csv(OUT / "tables" / "tile_profile.csv", index=False)

e = effects(sc, "atac_170")
ism = e[e.kind == "ism"].groupby(["site_id", "group", "pos"], as_index=False).eff.agg(imp=lambda v: np.abs(v).mean(), sgn="mean")
fig, axes = plt.subplots(1, 2, figsize=(13, 3.9), gridspec_kw={"width_ratios": [2, 1]})
ax = axes[0]
for g, lab, c in GROUPS:
    m = ism[ism.group == g].groupby("pos").imp.median()
    ax.plot(m.index, m.values, color=c, lw=1.8, label=lab)
ax.axvspan(-5.5, 5.5, color=GRID, alpha=0.6, lw=0)
ax.set_xlabel("Position relative to the motif midpoint (bp)")
ax.set_ylabel("Median |effect| of a single-base change on the 170 bp readout (log2FC)")
ax.set_title("Saturation mutagenesis: which bases matter", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)
conc = ism.groupby(["site_id", "group"]).apply(lambda d: d[d.pos.abs() <= 5].imp.sum() / d.imp.sum()).rename("motif_share").reset_index()
ax = axes[1]
for k, (g, lab, c) in enumerate(GROUPS):
    v = conc[conc.group == g].motif_share.values
    ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, patch_artist=True, boxprops={"facecolor": c, "alpha": 0.5}, medianprops={"color": INK})
ax.set_xticks(range(3), [g[1] for g in GROUPS], rotation=15, fontsize=8)
ax.axhline(11 / 61, color=MUTED, lw=0.9, ls=":")
ax.text(2.45, 11 / 61 + 0.01, "no-structure level 11/61", fontsize=7.5, color=MUTED, ha="right")
ax.set_ylabel("Share of the ±30 bp importance that lies in the 11 bp motif")
ax.set_title("How concentrated is the importance?", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(OUT / "figures" / "s2_saturation_mutagenesis.png", bbox_inches="tight")
plt.close(fig)
rows = []
for a, b in (("dynamic", "matched_stable"), ("dynamic", "high_stable"), ("matched_stable", "high_stable")):
    x, y = conc[conc.group == a].motif_share, conc[conc.group == b].motif_share
    rows.append({"comparison": f"{a} vs {b}", "median_a": x.median(), "median_b": y.median(), "p_mannwhitney": mannwhitneyu(x, y).pvalue})
pd.DataFrame(rows).to_csv(OUT / "tables" / "motif_importance_share.csv", index=False)
ism.to_csv(OUT / "tables" / "saturation_importance_by_position.csv", index=False)

# ============================================================ titrate
tt = load("titrate_cm", require_core=True)           # AP-1 titration only makes sense where the site has an AP-1 core
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.1), gridspec_kw={"width_ratios": [1.1, 1, 1]})
ax = axes[0]
cols = {"motif": "atac_motif", "170 bp": "atac_170"}
dose_rows = []
for (rd, col), ls in zip(cols.items(), ("-", "--")):
    e = effects(tt, col)
    k = e[e.kind == "core_mm"]
    for g, lab, c in GROUPS:
        med = k[k.group == g].groupby("pos").eff.median()
        ax.plot(med.index, med.values, ls, color=c, marker="o", ms=4, lw=1.5, label=f"{lab}, read over {rd}" if rd == "motif" else None)
        dose_rows += [{"readout": rd, "group": g, "core_bases_changed": int(p_), "median_effect": v_} for p_, v_ in med.items()]
ax.set_xlabel("Number of the 7 core bases changed"); ax.set_ylabel("Median effect (log2FC); solid = motif readout, dashed = 170 bp readout")
ax.legend(frameon=False, fontsize=7.5)
ax.set_title("Weakening the motif: dose-response", fontsize=10, loc="left")
ax = axes[1]
e170, emot = effects(tt, "atac_170"), effects(tt, "atac_motif")
single = [("core_consensus", "upgrade to\nTGAGTCA"), ("core_rc", "reverse\ncomplement"), ("core_halfsite", "half-site\nonly")]
stat_rows = []
for j, (cond, lab) in enumerate(single):
    for i, (g, glab, c) in enumerate(GROUPS):
        v = e170[(e170.condition == cond) & (e170.group == g)].eff
        lo, hi = boot(v.values)
        ax.errorbar(j + (i - 1) * 0.22, v.median(), yerr=[[v.median() - lo], [hi - v.median()]], fmt="o", color=c, ms=6, lw=1.4, capsize=0, label=glab if j == 0 else None)
        stat_rows.append({"condition": cond, "group": g, "n_sites": len(v), "median_effect_170": v.median(), "ci_lo": lo, "ci_hi": hi})
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xticks(range(len(single)), [s[1] for s in single])
ax.set_ylabel("Median effect on the 170 bp readout (log2FC)")
ax.set_title("Strength, orientation and half-site", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
pd.DataFrame(stat_rows).to_csv(OUT / "tables" / "core_variants.csv", index=False)
ax = axes[2]
rel_rows = []
nat = tt[tt.condition == "native"].set_index("site_id")["atac_170"]
pc = 0.05 * nat.median()
for g, lab, c in GROUPS:
    xs, ys = [], []
    for o in (-150, -80, -40, -20, 20, 40, 80, 150):
        d = e170[e170.group == g]
        for kind, ctrl in (("dup", "ctrl"), ("move", "ctrlmove")):
            a = d[d.condition == f"{kind}@{o:+d}"].set_index("site_id").eff
            b = d[d.condition == f"{ctrl}@{o:+d}"].set_index("site_id").eff
            net = (a - b).dropna()
            rel_rows.append({"group": g, "kind": kind, "offset": o, "n_sites": len(net), "median_net_effect_170": net.median()})
    r = pd.DataFrame(rel_rows)
    for kind, ls in (("dup", "-"), ("move", "--")):
        s = r[(r.group == g) & (r.kind == kind)].sort_values("offset")
        ax.plot(s.offset, s.median_net_effect_170, ls, color=c, marker="o", ms=4, lw=1.5, label=f"{lab}: {'duplicate' if kind == 'dup' else 'move'}")
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlabel("Offset of the inserted motif copy from the original (bp)")
ax.set_ylabel("Net effect on the 170 bp readout vs its shuffled control (log2FC)")
ax.set_title("Moving or duplicating the motif", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "figures" / "t1_titration.png", bbox_inches="tight")
plt.close(fig)
# accessibility at the inserted copy itself (11 bp at the offset), dup vs ctrl and move vs ctrlmove; offsets beyond +-80 lie outside the 170 bp readout
fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
tg_rows = []
natv = tt[tt.condition == "native"].set_index("site_id")["atac_170"]
pct = 0.05 * natv.median()
for ax, (kind, ctrl, ttl) in zip(axes, (("dup", "ctrl", "Duplicate the motif"), ("move", "ctrlmove", "Move the motif"))):
    for g, lab, c in GROUPS:
        ys, xs = [], []
        for o in (-150, -80, -40, -20, 20, 40, 80, 150):
            dd = tt[(tt.group == g)]
            a = dd[dd.condition == f"{kind}@{o:+d}"].groupby("site_id").atac_target.mean()
            b = dd[dd.condition == f"{ctrl}@{o:+d}"].groupby("site_id").atac_target.mean()
            r = np.log2((a + pct) / (b + pct)).dropna()
            xs.append(o); ys.append(r.median())
            tg_rows.append({"group": g, "kind": kind, "offset": o, "n_sites": len(r), "median_log2_target_gain": r.median()})
        ax.plot(xs, ys, color=c, marker="o", ms=4, lw=1.5, label=lab)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("Offset of the inserted motif copy (bp)")
    ax.set_title(ttl + ": accessibility at the copy vs a shuffled 11-mer", fontsize=9.5, loc="left")
axes[0].set_ylabel("Median log2 gain over the 11 bp at the offset")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "t2_target_accessibility.png", bbox_inches="tight")
plt.close(fig)
pd.DataFrame(tg_rows).to_csv(OUT / "tables" / "motif_copy_target_accessibility.csv", index=False)
pd.DataFrame(dose_rows).to_csv(OUT / "tables" / "core_dose_response.csv", index=False)
pd.DataFrame(rel_rows).drop_duplicates().to_csv(OUT / "tables" / "motif_relocation.csv", index=False)

# group comparisons (Dynamic vs context-matched Stable) for the headline quantities
rows = []
for nm, d, cond in (("dose: 3 core bases changed, 170 bp readout", e170, "core_mm3"), ("dose: 7 core bases changed, 170 bp readout", e170, "core_mm7"), ("consensus upgrade, 170 bp readout", e170, "core_consensus"),
                    ("reverse complement, 170 bp readout", e170, "core_rc"), ("half-site only, 170 bp readout", e170, "core_halfsite"), ("dose: 3 core bases changed, motif readout", emot, "core_mm3"), ("consensus upgrade, motif readout", emot, "core_consensus")):
    x = d[(d.condition == cond) & (d.group == "dynamic")].eff
    y = d[(d.condition == cond) & (d.group == "matched_stable")].eff
    if len(x) > 5 and len(y) > 5:
        rows.append({"quantity": nm, "n_dynamic": len(x), "n_matched_stable": len(y), "median_dynamic": x.median(), "median_matched_stable": y.median(), "p": mannwhitneyu(x, y).pvalue})
CMP = pd.DataFrame(rows)
CMP["fdr"] = bh(CMP.p.values) if len(CMP) else []
CMP.to_csv(OUT / "tables" / "dynamic_vs_matched_stable.csv", index=False)
print(CMP.round(4).to_string(index=False))
print("outputs in", OUT)
