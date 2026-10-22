#!/usr/bin/env python3
"""Check one merged fragment file before it becomes a 1 bp insertion target: is it the cell type and condition it claims, and does it agree with the
10 bp SeqDepthNorm bigWig we already trust?

Fragment file columns (ArchR export): chrom, start, end, barcode, duplicate count, sample#barcode. Reports
  - fragments, distinct cells, libraries (the sample prefix) and whether every library belongs to the expected condition
  - cell-type purity: share of cells (and of fragments) whose sample#barcode is labelled the expected cell type in Cell_barccodes.csv
  - fragment lengths (median, share < 100 bp, mono-nucleosome 180-247 bp), chromosomes
  - per-peak agreement: Tn5 insertions per bp in each snATAC peak (both fragment ends) against the bigWig mean over the peak (Spearman), and against the
    measured per-condition snATAC table where one exists. Same cells and processing should give rho >= 0.95 with the bigWig.
Writes <out>.json (summary) and <out>.peaks.csv.gz (per-peak values).

  python check_fragments.py --fragments CM_YS_merged_fragments_bgz_sorted.tsv.gz --cell Cardiomyocyte --cond YS --cells Cell_barccodes.csv \\
      --bigwig Cardiomyocyte_YS_merged.SeqDepthNorm.bw --sites ../Katana/data/measured_by_condition.csv --out checks/CM_YS
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def peaks_from_sites(path):
    """Unique snATAC peaks (chr-start-end strings) from the measured table, merged so they do not overlap, sorted per chromosome."""
    p = pd.read_csv(path, usecols=["snatac_peak"]).dropna().drop_duplicates()          # sites without a called peak have none
    p[["chrom", "start", "end"]] = p.snatac_peak.str.rsplit("-", n=2, expand=True)
    p[["start", "end"]] = p[["start", "end"]].astype(np.int64)
    out = {}
    for chrom, g in p.sort_values(["chrom", "start"]).groupby("chrom"):
        s, e = g.start.to_numpy(), g.end.to_numpy()
        if (s[1:] < e[:-1]).any():
            raise SystemExit(f"overlapping peaks on {chrom}: merge them first")
        out[chrom] = (s, e, g.snatac_peak.to_numpy())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--cell", required=True, help="expected cell_type_id in Cell_barccodes.csv, e.g. Cardiomyocyte")
    ap.add_argument("--cond", required=True, help="expected condition: YS, MAS, MAE or YE (library names start with it)")
    ap.add_argument("--cells", required=True, help="Cell_barccodes.csv (barcode_id = sample#barcode, cell_type_id)")
    ap.add_argument("--bigwig", required=True)
    ap.add_argument("--sites", required=True, help="table with an snatac_peak column (chr-start-end); measured_by_condition.csv also gives measured values")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunksize", type=int, default=5_000_000)
    args = ap.parse_args()

    peaks = peaks_from_sites(args.sites)
    ins = {c: np.zeros(len(v[0]), np.int64) for c, v in peaks.items()}
    n_frag, chroms, libs, per_cell = 0, Counter(), Counter(), Counter()
    lens = np.zeros(1001, np.int64)
    import gzip
    with (gzip.open if args.fragments.endswith(".gz") else open)(args.fragments, "rt") as fh:   # Cell Ranger-style '#' header lines; not comment="#",
        n_head = 0                                                                                # which would cut the sample#barcode column at the '#'
        for line in fh:
            if not line.startswith("#"):
                break
            n_head += 1
    reader = pd.read_csv(args.fragments, sep="\t", header=None, skiprows=n_head, usecols=[0, 1, 2, 5], names=["chrom", "start", "end", "cell"],
                         dtype={"chrom": str, "start": np.int64, "end": np.int64, "cell": str}, chunksize=args.chunksize)
    for ch in reader:
        n_frag += len(ch)
        chroms.update(ch.chrom.value_counts().to_dict())
        per_cell.update(ch.cell.value_counts().to_dict())
        lens += np.bincount(np.minimum(ch.end - ch.start, 1000), minlength=1001)
        for chrom, g in ch.groupby("chrom", sort=False):
            if chrom not in peaks:
                continue
            s, e, _ = peaks[chrom]
            for pos in (g.start.to_numpy(), g.end.to_numpy() - 1):        # both Tn5 insertions of each fragment
                i = np.searchsorted(s, pos, side="right") - 1
                ok = (i >= 0) & (pos < e[np.maximum(i, 0)])
                np.add.at(ins[chrom], i[ok], 1)
    for c, n in per_cell.items():
        libs[c.split("#")[0]] += n

    labels = pd.read_csv(args.cells).set_index("barcode_id").cell_type_id
    cells = pd.Series(per_cell)
    lab = labels.reindex(cells.index)
    kinds = lab.fillna("(not in Cell_barccodes)")
    purity_cells = float((lab == args.cell).mean())
    purity_frags = float(cells[lab == args.cell].sum() / cells.sum())

    import pyBigWig
    bw = pyBigWig.open(args.bigwig)
    bw_chroms = bw.chroms()
    rows = []
    for chrom, (s, e, names) in peaks.items():
        for k in range(len(s)):
            m = bw.stats(chrom, int(s[k]), int(e[k]), type="mean")[0] if chrom in bw_chroms else None
            rows.append({"peak": names[k], "width": int(e[k] - s[k]), "insertions": int(ins[chrom][k]), "bigwig_mean": m})
    pk = pd.DataFrame(rows)
    pk["ins_per_kb"] = pk.insertions / pk.width * 1000
    rho_bw = float(pk[["ins_per_kb", "bigwig_mean"]].dropna().corr(method="spearman").iloc[0, 1])
    meas_col = f"{'cm' if args.cell.lower().startswith('cardio') else 'fb'}_{args.cond.lower()}"
    meas = pd.read_csv(args.sites)
    rho_meas = None
    if meas_col in meas:
        pk = pk.merge(meas.groupby("snatac_peak")[meas_col].first().rename("measured").reset_index().rename(columns={"snatac_peak": "peak"}), on="peak", how="left")
        rho_meas = float(pk[["ins_per_kb", "measured"]].dropna().corr(method="spearman").iloc[0, 1])

    cum = np.cumsum(lens)
    summary = {
        "fragments_file": args.fragments, "expected": f"{args.cell} {args.cond}",
        "fragments": int(n_frag), "cells": int(len(cells)), "median_fragments_per_cell": float(cells.median()),
        "libraries": dict(libs.most_common()), "all_libraries_expected_condition": all(l.startswith(args.cond) for l in libs),
        "cell_type_purity_cells": purity_cells, "cell_type_purity_fragments": purity_frags,
        "cell_types_found": kinds.value_counts().head(6).to_dict(),
        "median_fragment_length": int(np.searchsorted(cum, cum[-1] / 2)), "share_under_100bp": float(lens[:100].sum() / n_frag),
        "share_mononucleosome_180_247bp": float(lens[180:248].sum() / n_frag),
        "chromosomes": len(chroms), "share_on_chrM": float(chroms.get("chrM", 0) / n_frag),
        "peaks": int(len(pk)), "share_of_insertions_in_peaks": float(pk.insertions.sum() / (2 * n_frag)),
        "spearman_vs_bigwig": rho_bw, "spearman_vs_measured_snatac": rho_meas, "measured_column": meas_col if rho_meas is not None else None,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".json").write_text(json.dumps(summary, indent=1))
    pk.to_csv(f"{out}.peaks.csv.gz", index=False)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
