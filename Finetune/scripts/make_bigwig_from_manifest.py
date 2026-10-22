#!/usr/bin/env python3
"""Runs make_insertion_bigwigs.py for one (cell type, condition) of a manifest, so the bigWigs can be built as a PBS array (one task per track).
Manifest (TSV with header): cell, condition, fragments, barcodes   -- one row per sample (replicate) of a condition; rows with the same cell and condition are pooled.
  python make_bigwig_from_manifest.py --manifest ../data/samples_manifest.tsv --chrom-sizes mm10.chrom.sizes --out-dir BIGWIGS --index 0     (--list shows the tracks, --dry-run prints the command)"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--index", type=int, default=None)
ap.add_argument("--list", action="store_true")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
m = pd.read_csv(a.manifest, sep="\t", comment="#")
groups = list(m.groupby(["cell", "condition"], sort=False))
if a.list or a.index is None:
    for i, ((c, k), g) in enumerate(groups):
        print(i, c, k, f"{len(g)} sample(s)")
    sys.exit(0)
(cell, cond), g = groups[a.index]
out = Path(a.out_dir) / f"{cell}_{cond}.insertions.bw"
cmd = [sys.executable, str(Path(__file__).parent / "make_insertion_bigwigs.py"), "--chrom-sizes", a.chrom_sizes, "--out", str(out)]
for _, r in g.iterrows():
    cmd += ["--sample", r.fragments, r.barcodes]
print(" ".join(cmd))
if not a.dry_run:
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    sys.exit(subprocess.run(cmd).returncode)
