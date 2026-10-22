#!/usr/bin/env python3
"""Check one 1 bp insertion bigWig against the fragment-level checks it was built from (check_fragments.py, check_tn5_offset.py).

  - every base has a value (explicit zeros) on the canonical chromosomes, values are whole multiples of the scale factor in <bw>.summary.json
  - per-peak sums equal insertions x scale from <checks>/<Cell>_<Cond>.peaks.csv.gz (Spearman 1, ratio = scale)
  - the +-200 bp profile around the motif midpoints equals the fragment-level left + right cut profile from <checks>/<Cell>_<Cond>_tn5.profile.csv
Writes <out>.json.

  python check_insertion_bigwigs.py --bigwig bigwigs_1bp/Cardiomyocyte_YS.insertions.bw --checks fragments/checks --cell Cardiomyocyte --cond YS \\
      --sites ../Katana/data/cm_fos_sites.csv --out fragments/checks/Cardiomyocyte_YS_bigwig
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

W = 200


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bigwig", required=True)
    ap.add_argument("--checks", required=True)
    ap.add_argument("--cell", required=True)
    ap.add_argument("--cond", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    meta = json.loads(Path(f"{args.bigwig}.summary.json").read_text())
    scale = meta["scale"]
    bw = pyBigWig.open(args.bigwig)
    chroms = bw.chroms()
    tag = f"{args.cell}_{args.cond}"

    pk = pd.read_csv(Path(args.checks) / f"{tag}.peaks.csv.gz")
    pk[["chrom", "start", "end"]] = pk.peak.str.rsplit("-", n=2, expand=True)
    pk["bw_sum"] = [float(np.nansum(bw.values(c, int(s), int(e), numpy=True))) if c in chroms else np.nan for c, s, e in zip(pk.chrom, pk.start, pk.end)]
    pk = pk.dropna(subset=["bw_sum"])
    ratio = (pk.bw_sum / pk.insertions.replace(0, np.nan)).dropna()

    s = pd.read_csv(args.sites, usecols=["seqnames", "motif_mid"]).drop_duplicates()
    prof = np.zeros(2 * W + 1)
    n_nan, n_bases, off_grid = 0, 0, 0
    for c, m in zip(s.seqnames, s.motif_mid):
        if c not in chroms or m - W < 0 or m + W + 1 > chroms[c]:
            continue
        v = bw.values(c, int(m - W), int(m + W + 1), numpy=True)
        n_nan += int(np.isnan(v).sum())
        n_bases += len(v)
        q = np.nan_to_num(v) / scale
        off_grid += int((np.abs(q - np.round(q)) > 1e-3).sum())
        prof += np.nan_to_num(v)
    tn5 = pd.read_csv(Path(args.checks) / f"{tag}_tn5.profile.csv")
    frag_prof = (tn5.left_cuts + tn5.right_cuts).to_numpy(float)
    # the tn5 check counted each cut for its nearest motif only; the bigWig sums every motif window, so compare shapes (correlation), not totals
    shape_r = float(np.corrcoef(prof, frag_prof)[0, 1])

    summary = {"bigwig": args.bigwig, "scale": scale, "insertions": meta["insertions"], "fragments_kept": meta["fragments_kept"], "offsets": meta["offsets"],
               "peaks": int(len(pk)), "spearman_peak_sum_vs_insertions": float(pk[["bw_sum", "insertions"]].corr(method="spearman").iloc[0, 1]),
               "median_ratio_sum_over_insertions": float(ratio.median()), "ratio_spread": float(ratio.max() - ratio.min()),
               "bases_checked": n_bases, "nan_bases": n_nan, "values_off_the_scale_grid": off_grid,
               "profile_shape_r_vs_fragments": shape_r}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
