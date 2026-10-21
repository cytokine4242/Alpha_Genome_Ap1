#!/usr/bin/env python3
"""Which tracks best predict measured snATAC? Native predictions (results/validate_*) vs measured pseudobulk accessibility
(Katana/data/measured_atac_sites.csv), for mouse tracksets (CM, FB) and human tracksets (CMh, FBh) applied to the mouse sequence.

Prints (1) Spearman rho of every accessibility readout (trackset x head x region) with the measured accessibility of the
cell type of the site table, with 95% CI, best first; (2) cell-type specificity: does the predicted FB - CM difference
track the measured FB - CM difference per site?

  python3 21_validate_tracks.py [results_dir] [prefix]      (default: Katana/results, validate_)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Katana" / "results"
PRE = sys.argv[2] if len(sys.argv) > 2 else "validate_"
pd.set_option("display.width", 200, "display.max_rows", 200)
meas = pd.read_csv(ROOT / "Katana" / "data" / "measured_atac_sites.csv")
ACC = [f"{h}_{r}" for h in ("atac", "dnase") for r in ("motif", "170", "peak")]


def rho_ci(a, b):
    m = a.notna() & b.notna()
    n = int(m.sum())
    if n < 10:
        return np.nan, np.nan, np.nan, n
    r = spearmanr(a[m], b[m])[0]
    z, se = np.arctanh(r), 1.06 / np.sqrt(n - 3)  # Fieller approximation for Spearman
    return r, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se), n


rows, wide = [], {}
for tab, col in (("cm", "cm_cpm_per_kb"), ("fb", "fb_cpm_per_kb")):
    f = RES / f"{PRE}{tab}_chunk0.csv"
    if not f.exists():
        print("missing", f)
        continue
    d = pd.read_csv(f)
    d = d[d.condition == "native"]
    m = meas[meas.site_table == tab].set_index("site_id")[["cm_cpm_per_kb", "fb_cpm_per_kb"]]
    for ts, g in d.groupby("trackset"):
        g = g.set_index("site_id").join(m)
        for c in ACC:
            if c in g and g[c].notna().any():
                r, lo, hi, n = rho_ci(g[c], g[col])
                rows.append({"sites": tab.upper(), "trackset": ts, "readout": c, "rho": r, "ci_lo": lo, "ci_hi": hi, "n": n})
        wide[(tab, ts)] = g

print("=== 1. native prediction vs measured snATAC of the same cell type (Spearman)")
t = pd.DataFrame(rows)
for tab in ("CM", "FB"):
    x = t[t.sites == tab].sort_values("rho", ascending=False)
    print(f"\n{tab} sites vs measured {tab} accessibility")
    print(x.drop(columns="sites").round(2).to_string(index=False))

print("\n=== 2. cell-type specificity: predicted (FB track - CM track) vs measured (FB - CM), log2 CPM/kb, per site")
for tab in ("cm", "fb"):
    for label, tf, tc, col in (("mouse atac (Patski - heart)", "FB", "CM", "atac_peak"), ("mouse dnase (Patski - heart)", "FB", "CM", "dnase_peak"),
                               ("human dnase (cardiac fibroblast - cardiac muscle cell)", "FBh", "CMh", "dnase_peak")):
        if (tab, tf) in wide and (tab, tc) in wide and col in wide[(tab, tf)] and col in wide[(tab, tc)]:
            f, c = wide[(tab, tf)], wide[(tab, tc)]
            pred = np.log2(f[col] + 0.01) - np.log2(c[col].reindex(f.index) + 0.01)
            obs = np.log2(f.fb_cpm_per_kb + 1) - np.log2(f.cm_cpm_per_kb + 1)
            r, lo, hi, n = rho_ci(pred, obs)
            print(f"{tab.upper()} sites | {label:55s} rho={r:+.2f} [{lo:+.2f}, {hi:+.2f}] n={n}")
