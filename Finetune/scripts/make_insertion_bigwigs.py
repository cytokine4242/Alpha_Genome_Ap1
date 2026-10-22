#!/usr/bin/env python3
"""Fragments -> 1 bp Tn5 insertion-site bigWig, the ATAC target format AlphaGenome was trained on (base-pair resolution, 5' cut sites).

Why not use the existing SeqDepthNorm bigWigs? They come from `bamCoverage --binSize 10 --normalizeUsing CPM`: read coverage in 10 bp bins, which blurs the base-pair
footprints and motif-level structure that the perturbation analyses rely on. Insertion sites at 1 bp keep that structure (and can be pooled to 128 bp by the model).

Each fragment [start, end) gives two insertions. Fragment files from Cell Ranger ATAC, sinto or ArchR are already corrected for the Tn5 9 bp duplication, so by default the
insertions are the fragment ends themselves: start and end - 1 (0-based), i.e. --offsets 0 -1. Only for fragments built from raw read ends without correction use --offsets 4 -5.
Both strands are summed.
Several samples can be pooled (a fragment file plus a file of the cell barcodes to keep for that sample, one barcode per line), so replicates of a condition are merged
and barcode collisions between samples cannot happen. BARCODES may be several files joined by commas (for example the _F and _M lists of one library);
a 'CB:Z:' prefix (as in lists taken from BAM tags) is removed so they match the fragment file. Positions without an insertion are written as explicit zeros. Values are scaled so the genome-wide total equals --target-total (default 1e9) to make conditions comparable;
the fine-tuning code divides by per-track means anyway, so any constant works as long as every track uses the same one.

  python make_insertion_bigwigs.py --chrom-sizes mm10.chrom.sizes --out Cardiomyocyte_YS.insertions.bw \\
      --sample YS1_4_fragments.tsv.gz cellID_Cardiomyocyte_YS1_4.txt --sample YS3_2_fragments.tsv.gz cellID_Cardiomyocyte_YS3_2.txt
Fragment files: tab-separated chrom, start, end, barcode[, count], gzipped or plain, any order (the whole set of insertions is held in memory: about 8 bytes per insertion).
Writes <out>.summary.json (fragments read, kept, insertions, scale factor)."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

CANONICAL = [f"chr{i}" for i in list(range(1, 20)) + ["X", "Y"]]


def read_sizes(path, canonical_only=True):
    sizes = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 2 and (not canonical_only or p[0] in CANONICAL):
            sizes[p[0]] = int(p[1])
    return sizes


def collect_cuts(samples, sizes, offsets, chunksize):
    cuts, n_read, n_kept = defaultdict(list), 0, 0
    for frag, bc in samples:
        barcodes = {b[5:] if b.startswith("CB:Z:") else b for f in bc.split(",") for b in Path(f).read_text().split()}
        reader = pd.read_csv(frag, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "bc"], dtype={"chrom": str, "start": np.int64, "end": np.int64, "bc": str}, chunksize=chunksize)
        for ch in reader:
            n_read += len(ch)
            ch = ch[ch.bc.isin(barcodes) & ch.chrom.isin(sizes)]
            n_kept += len(ch)
            for chrom, g in ch.groupby("chrom"):
                cuts[chrom].append((g.start.values + offsets[0]).astype(np.int64))
                cuts[chrom].append((g.end.values + offsets[1]).astype(np.int64))
    return cuts, n_read, n_kept


def write_bigwig(cuts, sizes, out, target_total):
    counts = {}
    total = 0
    for chrom in sizes:
        if chrom in cuts:
            pos = np.concatenate(cuts[chrom])
            pos = pos[(pos >= 0) & (pos < sizes[chrom])]
            u, c = np.unique(pos, return_counts=True)
            counts[chrom] = (u, c)
            total += int(c.sum())
    if total == 0:
        raise SystemExit("no insertions found: check the barcode lists and chromosome names")
    scale = target_total / total
    bw = pyBigWig.open(str(out), "w")
    bw.addHeader([(c, s) for c, s in sizes.items()], maxZooms=10)
    for chrom in sizes:                       # entries must follow the header order
        if chrom not in counts:
            u, c = np.array([], dtype=np.int64), np.array([], dtype=np.float64)
        else:
            u, c = counts[chrom]
        # explicit zeros between insertions, so every base of the chromosome has a value (like bamCoverage output); loaders then never see missing data
        gs, ge = np.concatenate([[0], u + 1]), np.concatenate([u, [sizes[chrom]]])
        keep = ge > gs
        starts = np.concatenate([u, gs[keep]])
        ends = np.concatenate([u + 1, ge[keep]])
        vals = np.concatenate([(c * scale), np.zeros(int(keep.sum()))]).astype(np.float32)
        order = np.argsort(starts, kind="stable")
        starts, ends, vals = starts[order], ends[order], vals[order]
        for i in range(0, len(starts), 1_000_000):
            sl = slice(i, i + 1_000_000)
            bw.addEntries([chrom] * len(starts[sl]), starts[sl].astype(np.int64).tolist(), ends=ends[sl].astype(np.int64).tolist(), values=vals[sl].tolist())
    bw.close()
    return total, scale


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--sample", nargs=2, action="append", metavar=("FRAGMENTS", "BARCODES"), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target-total", type=float, default=1e9)
    ap.add_argument("--offsets", type=int, nargs=2, default=[0, -1], help="insertion offsets from fragment start and end (default 0 -1: fragments already Tn5-corrected; use 4 -5 for uncorrected)")
    ap.add_argument("--chunksize", type=int, default=5_000_000)
    ap.add_argument("--all-chroms", action="store_true", help="keep every chromosome of the sizes file (default chr1-19, X, Y)")
    args = ap.parse_args()
    sizes = read_sizes(args.chrom_sizes, canonical_only=not args.all_chroms)
    cuts, n_read, n_kept = collect_cuts(args.sample, sizes, args.offsets, args.chunksize)
    total, scale = write_bigwig(cuts, sizes, args.out, args.target_total)
    summary = {"fragments_read": n_read, "fragments_kept": n_kept, "insertions": total, "scale": scale, "target_total": args.target_total, "offsets": args.offsets, "samples": [list(s) for s in args.sample]}
    Path(f"{args.out}.summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
