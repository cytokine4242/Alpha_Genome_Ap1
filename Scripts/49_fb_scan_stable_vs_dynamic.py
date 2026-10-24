#!/usr/bin/env python3
"""What the fibroblast tile scan and saturation mutagenesis say about Stable vs Dynamic AP-1 sites
(Katana/results_ft/scan_fb_chunk*.csv, fine-tuned 1 bp fibroblast head; groups: all Dynamic sites and the same number of
Stable sites drawn at random).

Effect = per site and position, mean over replicates of log2((edit + pc) / (native + pc)).

  1  tile profile: replacing each 10 bp tile from -300 to +290 bp, by class - where along the sequence does the
     Stable/Dynamic difference actually live?
  2  the same profile normalised per site (each site's tile effects divided by its own total), which separates
     "Dynamic sites are uniformly less sensitive" from "Dynamic sites are sensitive over a different region"
  3  how concentrated the dependence is: the share of a site's total tile importance within +-25 bp of the motif, and the
     half-width of its profile (the distance within which half the importance lies)
  4  saturation mutagenesis at single-base resolution (-30 to +30 bp): per-position importance by class, and whether the
     AP-1 core bases carry the same weight in both classes
Outputs Analysis/fb_scan/{tables,figures}.   python3 49_fb_scan_stable_vs_dynamic.py [--res DIR] [--simulate]
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
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--readout", default="atac_170", help="atac_170 (default), atac_motif or atac_peak")
ap.add_argument("--simulate", action="store_true")
args = ap.parse_args()
INK, MUTED, SURF, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
GROUPS = [("dynamic", "Dynamic", ORANGE), ("random_stable", "Stable (random)", BLUE)]


def simulate(folder):
    """Dynamic sites: same motif-centred dependence, but a weaker and narrower contribution from the flanks."""
    rng = np.random.default_rng(49)
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(60):
        g = "dynamic" if i % 2 else "random_stable"
        nat = (3.0 if g == "random_stable" else 1.8) * np.exp(rng.normal(0, 0.3))
        base = {"site_id": f"s{i}", "group": g, "stable_dynamic": "Dynamic" if g == "dynamic" else "Stable",
                "trackset": args.track, "native_170_pred": nat, "core_found": True}
        rows.append({**base, "condition": "native", "kind": "native", "pos": 0, "rep": 0,
                     "atac_motif": nat, "atac_170": nat, "atac_peak": nat})
        wide = 90 if g == "random_stable" else 35            # the planted difference: flank width
        for o in range(-300, 300, 10):
            for r in range(3):
                f = 1 - (0.25 * np.exp(-abs(o + 5) / 12) + 0.35 * np.exp(-abs(o + 5) / wide)) * np.exp(rng.normal(0, 0.1))
                rows.append({**base, "condition": f"tile@{o}", "kind": "tile", "pos": o, "rep": r,
                             "atac_motif": nat * f, "atac_170": nat * f, "atac_peak": nat * f})
        for p in range(-30, 31):
            for a in "ACG":
                f = 1 - (0.30 if abs(p) <= 5 else 0.03) * np.exp(rng.normal(0, 0.2))
                rows.append({**base, "condition": f"ism@{p}{a}", "kind": "ism", "pos": p, "rep": 0,
                             "atac_motif": nat * f, "atac_170": nat * f, "atac_peak": nat * f})
    pd.DataFrame(rows).to_csv(folder / "scan_fb_chunk0.csv", index=False)


if args.simulate:
    tmp = Path(tempfile.mkdtemp())
    simulate(tmp)
    RES, OUT = tmp, tmp / "out"
else:
    RES = Path(args.res) if args.res else ROOT / "Katana" / "results_ft"
    OUT = ROOT / "Analysis" / "fb_scan"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

files = sorted(glob.glob(str(RES / "scan_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no scan_fb_chunk*.csv in {RES}")
cols = ["site_id", "group", "stable_dynamic", "core_found", "condition", "kind", "pos", "rep", "trackset",
        "atac_motif", "atac_170", "atac_peak"]
d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in files], ignore_index=True)
d = d[d.trackset == args.track]
col = args.readout
nat = d[d.kind == "native"].drop_duplicates("site_id").set_index("site_id")[col]
pc = 0.05 * nat.median()
e = d[d.kind != "native"].copy()
e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
e = e.groupby(["site_id", "group", "kind", "pos"], as_index=False).eff.mean()
n_sites = e.groupby("group").site_id.nunique().to_dict()
print(f"{len(files)} chunks, sites per group {n_sites}, readout {col}, trackset {args.track}")

# ---------------------------------------------------------------- 1 & 2. tile profile, raw and per-site normalised
tile = e[e.kind == "tile"].copy()
tot = tile.assign(a=tile.eff.abs()).groupby("site_id").a.sum().rename("total")
tile = tile.join(tot, on="site_id")
tile["shape"] = tile.eff.abs() / tile.total.replace(0, np.nan)

prof_rows = []
for pos, g in tile.groupby("pos"):
    row = {"pos": pos}
    for key, lab, _ in GROUPS:
        s = g[g.group == key]
        row[f"median_{key}"] = s.eff.median()
        row[f"shape_{key}"] = s["shape"].median()
    dy, st = g[g.group == "dynamic"].eff, g[g.group == "random_stable"].eff
    row["difference"] = dy.median() - st.median()
    row["p"] = mannwhitneyu(dy, st).pvalue if len(dy) > 5 and len(st) > 5 else np.nan
    sdy, sst = g[g.group == "dynamic"]["shape"].dropna(), g[g.group == "random_stable"]["shape"].dropna()
    row["p_shape"] = mannwhitneyu(sdy, sst).pvalue if len(sdy) > 5 and len(sst) > 5 else np.nan
    prof_rows.append(row)
PROF = pd.DataFrame(prof_rows).sort_values("pos")
PROF.to_csv(OUT / "tables" / "tile_profile_by_class.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
for key, lab, c in GROUPS:
    axes[0].plot(PROF.pos, PROF[f"median_{key}"], color=c, lw=1.6, label=f"{lab} (n={n_sites.get(key, 0)})")
    axes[1].plot(PROF.pos, PROF[f"shape_{key}"], color=c, lw=1.6, label=lab)
for ax, ylab, ttl in ((axes[0], "Effect of replacing the tile (log2)", "Raw sensitivity"),
                      (axes[1], "Share of the site's own total", "Shape: each site normalised to its own total")):
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("Position of the 10 bp tile relative to the motif (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(ttl, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "s1_tile_profile_by_class.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- 3. concentration of the dependence
# |effect| summed over tiles counts noise as signal, and the noise floor differs by class: Dynamic sites have about half
# the native accessibility of Stable ones, so their log2 ratios are noisier and a raw sum makes their profiles look wider
# and less concentrated than they are. Each site's own noise floor is estimated from its distal tiles (|pos| >= 200 bp,
# where the median effect is ~0 in both classes) and subtracted before anything is summed.
BASE_FROM = 200
floor = tile[tile.pos.abs() >= BASE_FROM].assign(a=lambda t: t.eff.abs()).groupby("site_id").a.median().rename("floor")
tile = tile.join(floor, on="site_id")
tile["signal"] = (tile.eff.abs() - tile.floor).clip(lower=0)
fl = tile.drop_duplicates("site_id").groupby("group").floor.median()
print(f"\nnoise floor per site (median |effect| at |pos| >= {BASE_FROM} bp): " +
      ", ".join(f"{k} {v:.4f}" for k, v in fl.items()) + "  <- subtracted below")
near = tile[tile.pos.abs() <= 25].groupby(["site_id", "group"]).signal.sum().rename("near")
allt = tile.groupby(["site_id", "group"]).signal.sum().rename("all")
conc = pd.concat([near, allt], axis=1).reset_index()
conc["share_within_25bp"] = conc.near / conc["all"].replace(0, np.nan)


def half_width(g):
    """Smallest distance containing half of a site's total tile importance."""
    s = g.sort_values("absdist")
    c = s.signal.cumsum() / s.signal.sum()
    hit = s.absdist[c >= 0.5]
    return hit.iloc[0] if len(hit) else np.nan


tile["absdist"] = tile.pos.abs()
hw = tile.groupby(["site_id", "group"], group_keys=False)[["absdist", "signal"]].apply(half_width).rename("half_width_bp").reset_index()
CONC = conc.merge(hw, on=["site_id", "group"])
rows = []
for metric in ("share_within_25bp", "half_width_bp", "all"):
    dy = CONC[CONC.group == "dynamic"][metric].dropna()
    st = CONC[CONC.group == "random_stable"][metric].dropna()
    rows.append({"metric": {"all": "total importance (noise-corrected)"}.get(metric, metric),
                 "n_dynamic": len(dy), "n_stable": len(st), "median_dynamic": dy.median(), "median_stable": st.median(),
                 "p": mannwhitneyu(dy, st).pvalue if len(dy) > 5 and len(st) > 5 else np.nan})
CONCT = pd.DataFrame(rows)
CONCT.to_csv(OUT / "tables" / "concentration_by_class.csv", index=False)
CONC.to_csv(OUT / "tables" / "concentration_per_site.csv", index=False)
print("\n=== concentration of the dependence (site is the unit)")
print(CONCT.round(4).to_string(index=False))

# ---------------------------------------------------------------- 4. saturation mutagenesis
ism = e[e.kind == "ism"]
rows = []
for pos, g in ism.groupby("pos"):
    dy, st = g[g.group == "dynamic"].eff, g[g.group == "random_stable"].eff
    rows.append({"pos": pos, "median_dynamic": dy.median(), "median_stable": st.median(),  # ISM columns: dynamic / stable
                 "difference": dy.median() - st.median(),
                 "p": mannwhitneyu(dy, st).pvalue if len(dy) > 5 and len(st) > 5 else np.nan})
ISM = pd.DataFrame(rows).sort_values("pos")
ISM.to_csv(OUT / "tables" / "saturation_by_position.csv", index=False)

core = ism[ism.pos.abs() <= 3].groupby(["site_id", "group"]).eff.mean().reset_index()
flank = ism[(ism.pos.abs() >= 10) & (ism.pos.abs() <= 30)].groupby(["site_id", "group"]).eff.mean().reset_index()
rows = []
for lab, t in (("core bases (+-3 bp)", core), ("nearby bases (10-30 bp)", flank)):
    dy, st = t[t.group == "dynamic"].eff, t[t.group == "random_stable"].eff
    rows.append({"region": lab, "n_dynamic": len(dy), "n_stable": len(st), "median_dynamic": dy.median(),
                 "median_stable": st.median(), "difference": dy.median() - st.median(),
                 "p": mannwhitneyu(dy, st).pvalue if len(dy) > 5 and len(st) > 5 else np.nan})
ISMT = pd.DataFrame(rows)
ISMT.to_csv(OUT / "tables" / "saturation_core_vs_flank.csv", index=False)
print("\n=== single-base importance")
print(ISMT.round(4).to_string(index=False))

fig, ax = plt.subplots(figsize=(8.5, 4))
for key, lab, c in GROUPS:
    ax.plot(ISM.pos, ISM["median_dynamic" if key == "dynamic" else "median_stable"], color=c, lw=1.5, label=lab)
ax.axvline(0, color=MUTED, lw=0.8, ls="--")
ax.axvspan(-5.5, 5.5, color=GRID, alpha=0.5, zorder=0)
ax.set_xlabel("Position relative to the motif midpoint (bp)")
ax.set_ylabel("Effect of changing one base (log2)")
ax.set_title("Single-base importance; shaded = the 11 bp motif", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "s2_saturation_by_class.png", bbox_inches="tight")
plt.close(fig)

print("\n=== where along the sequence the classes differ (raw effect, tiles)")
sig = PROF[(PROF.p < 0.05) & (PROF.pos.abs() <= 150)]
print(f"  tiles with p<0.05 within +-150 bp: {len(sig)} of {len(PROF[PROF.pos.abs() <= 150])}")
for lo, hi in ((0, 25), (25, 50), (50, 100), (100, 200), (200, 300)):
    s = PROF[(PROF.pos.abs() >= lo) & (PROF.pos.abs() < hi)]
    print(f"  |{lo:3d}-{hi:3d}| bp: Dynamic {s[[c for c in s if c.startswith('median_dyn')][0]].median():+.3f}  "
          f"Stable {s[[c for c in s if c.startswith('median_ran')][0]].median():+.3f}  "
          f"difference {s.difference.median():+.3f}")
print("\noutputs in", OUT)
