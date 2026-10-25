#!/usr/bin/env python3
"""What an element builds in a closed host, against the empty host, at both transplanted widths.

Two changes from the version on the slide (make_fb_simple_figures.py):

  1  the baseline. That version scores the element against a dinucleotide-shuffled copy of itself in the same host,
     which controls base composition but is not a no-insert control: 20_host_native.py shows a shuffled insert
     already opens these hosts about 3.3x. Against the empty host the lift is what inserting the element actually
     does, which is the question the slide asks. Both baselines are written out; the scramble one stays as the
     composition control.
  2  both widths, not only the peak. 170 bp is the MPRA-sized window and the whole peak is a median ~930 bp, so the
     pair says how the effect scales with how much of the element is moved, rather than giving one number.

Per donor the median over the 20 hosts, then per class the median over donors with a bootstrap 95% CI over donors.
Hosts are shared by every donor, so the host contributes equally to all classes and cancels from the comparison.

Outputs Analysis/fb_swap/tables/host_insertion_by_width.csv.
  python3 61_fb_host_insertion_widths.py [--res DIR] [--boot 2000]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=None)
ap.add_argument("--track", default="FB")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=61)
a = ap.parse_args()
RES = Path(a.res) if a.res else ROOT / "Katana" / "results_ft"
OUT = ROOT / "Analysis" / "fb_swap" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)
CLS5 = ["Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down"]

hf = RES / "hostnat_fb.csv"
if not hf.exists():
    raise SystemExit(f"{hf} missing - run Katana/scripts/20_host_native.py first (submit_all.sh fb_hostnat)")
h = pd.read_csv(hf)
h = h[h.trackset == a.track][["host_id", "insert_bp", "atac_insert"]].rename(columns={"atac_insert": "host"})

files = sorted(glob.glob(str(RES / "bg_fb_chunk*.csv")))
if not files:
    raise SystemExit(f"no bg_fb_chunk*.csv in {RES}")
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d = d[d.trackset == a.track]

rows, per_donor = [], []
for w, dw in d.groupby("width"):
    p = dw.pivot_table(index=["donor_id", "donor_class", "host_id", "insert_bp"],
                       columns="kind", values="atac_insert").reset_index()
    p = p.merge(h, on=["host_id", "insert_bp"], how="left")
    miss = p.host.isna().sum()
    if miss:
        print(f"  width {w}: {miss} rows with no matching host baseline")
    pc = 0.05 * p[["element", "shuffled", "host"]].stack().median()
    p["vs_empty"] = np.log2((p.element + pc) / (p.host + pc))
    p["vs_scramble"] = np.log2((p.element + pc) / (p.shuffled + pc))
    p["scramble_alone"] = np.log2((p.shuffled + pc) / (p.host + pc))
    g = (p.groupby(["donor_class", "donor_id"])[["vs_empty", "vs_scramble", "scramble_alone"]].median()
         .join(p.groupby(["donor_class", "donor_id"]).insert_bp.median()).reset_index())
    g["width"] = w
    per_donor.append(g)

    for cls in CLS5 + ["Dynamic"]:
        sub = g[g.donor_class != "Stable"] if cls == "Dynamic" else g[g.donor_class == cls]
        if sub.empty:
            continue
        r = {"width": w, "class": cls, "n_donors": len(sub), "median_insert_bp": float(sub.insert_bp.median())}
        for q in ("vs_empty", "vs_scramble", "scramble_alone"):
            v = sub[q].to_numpy()
            b = np.array([np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(a.boot)])
            r[q] = float(np.median(v))
            r[q + "_lo"] = float(np.percentile(b, 2.5))
            r[q + "_hi"] = float(np.percentile(b, 97.5))
        if cls != "Stable":
            st = g[g.donor_class == "Stable"]
            for q in ("vs_empty", "vs_scramble"):
                r["p_" + q] = mannwhitneyu(st[q], sub[q]).pvalue
        rows.append(r)

T = pd.DataFrame(rows)
T.to_csv(OUT / "host_insertion_by_width.csv", index=False)
pd.concat(per_donor, ignore_index=True).to_csv(OUT / "host_insertion_by_width_per_donor.csv", index=False)

show = T[["width", "class", "n_donors", "median_insert_bp", "vs_empty", "vs_scramble", "scramble_alone"]].copy()
for c in ("vs_empty", "vs_scramble", "scramble_alone"):
    show["fold_" + c.replace("vs_", "")] = (2 ** show[c]).round(1)
print("\n=== fold over each baseline, median per donor then per class")
print(show[["width", "class", "n_donors", "median_insert_bp",
            "fold_empty", "fold_scramble", "fold_scramble_alone"]].to_string(index=False))
for w in T.width.unique():
    s = T[(T.width == w) & (T["class"] == "Stable")].iloc[0]
    y = T[(T.width == w) & (T["class"] == "Dynamic")].iloc[0]
    print(f"  {w:>4}: Stable {2 ** s.vs_empty:.1f}x vs Dynamic {2 ** y.vs_empty:.1f}x over the empty host, "
          f"ratio {2 ** (s.vs_empty - y.vs_empty):.1f}x, p = {y.p_vs_empty:.3g}")
print("\nwrote", OUT / "host_insertion_by_width.csv", "| _per_donor.csv")
