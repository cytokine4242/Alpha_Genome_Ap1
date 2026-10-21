#!/usr/bin/env python3
"""MPRA elements predicted in a NEUTRAL context. The MPRA tests each 170 bp element outside its genome; the native-context predictions (13_mpra.py)
include the element's own chromatin neighbourhood. Here the 170 bp element (native, motif-scrambled, fully scrambled) is inserted at the centre of a
few fixed neutral backgrounds (random mm10 windows, dinucleotide-shuffled once, so they carry no AP-1 structure) and AlphaGenome's accessibility
over the insert is read out. If predicted accessibility of the insert alone tracks the MPRA better than the genomic-context prediction did, the earlier
mismatch was about genomic context. The backgrounds are identical in every chunk (fixed seed).

  python 14_mpra_neutral.py --out ../results/mpra_neutral --chunk 0 --n-chunks 8
"""
from __future__ import annotations

import argparse
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L


def make_backgrounds(genome, size, n, seed):
    chroms = [f"chr{i}" for i in range(1, 20)] + ["chrX"]
    lens = np.array([genome.length(c) for c in chroms], dtype=float)
    rng = np.random.default_rng([seed, 991])
    out = []
    while len(out) < n:
        chrom = chroms[rng.choice(len(chroms), p=lens / lens.sum())]
        mid = int(rng.integers(size, genome.length(chrom) - size))
        codes, _ = genome.window(chrom, mid, size)
        if codes is None or L.has_n(codes, 0, size):
            continue
        out.append(L.dinuc_shuffle_codes(codes, np.random.default_rng([seed, len(out)])))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--reps", type=int, default=5, help="shuffle replicates per condition")
    ap.add_argument("--n-bg", type=int, default=2, help="number of neutral backgrounds")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    seed = cfg["mpra"]["seed"]
    genome = L.Genome(cfg["paths"]["genome"])
    ag = L.AG(cfg)
    bgs = make_backgrounds(genome, size, args.n_bg, seed)
    c = size // 2

    mpra = pd.read_csv(L.ROOT / "data" / "MPRA_Results_with_TSS_distance.csv").drop_duplicates("Coords").reset_index(drop=True)
    xm = pd.read_csv(L.ROOT / "data" / "fb_mpra_tobias_crossmap.csv")
    xm = xm[xm.overlap_status == "Matched"]
    mpra = mpra.iloc[args.chunk::args.n_chunks]
    if args.limit:
        mpra = mpra.sample(n=min(args.limit, len(mpra)), random_state=seed)
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    skip = L.done_ids(out, "mpra_id")
    t0, n_done = time.time(), 0
    for _, e in mpra.iterrows():
        eid = e["Coords"]
        if eid in skip:
            continue
        start, end = int(e["start"]), int(e["end"])
        g0 = start - 1 if end - start + 1 == 170 else start  # Coords are 1-based inclusive when 170 bp wide
        seq = genome.fetch(e["chr"], g0, g0 + 170)
        if seq is None or len(seq) != 170 or L.has_n(seq, 0, 170):
            continue
        motif, source = None, "none"
        hits = xm[xm.mpra_id == eid]
        mid = (start + end) // 2
        if len(hits):
            h = hits.iloc[np.argmin(np.abs((hits.motif_start + hits.motif_end) / 2 - mid))]
            a, b = int(h.motif_start) - g0, int(h.motif_end) - g0
            if 0 <= a < b <= 170:
                motif, source = (a, b), "tobias"
        if motif is None:
            m = L.CORE_RE.search(L.decode(seq))
            if m:
                motif, source = (max(m.start() - 2, 0), min(m.end() + 2, 170)), "regex"
        rng = np.random.default_rng([seed, zlib.crc32(eid.encode())])
        variants = [("native", 0, seq)]
        for r in range(args.reps):
            if motif and not L.has_n(seq, *motif):
                variants.append(("motif_shuffle", r, L.shuffle_region(seq, motif[0], motif[1], rng)))
        for r in range(args.reps):
            variants.append(("full_shuffle", r, L.shuffle_region(seq, 0, 170, rng)))
        seqs, meta = [], []
        for bi, bg in enumerate(bgs):
            for cond, rep, ins in variants:
                w = bg.copy()
                w[c - 85:c + 85] = ins
                seqs.append(w)
                meta.append((bi, cond, rep))
        reg = {"motif": (c - 85 + motif[0], c - 85 + motif[1]) if motif else (c - 85, c - 85), "w170": (c - 85, c + 85), "peak": (c - 85, c + 85), "mid": c}
        preds = ag.predict(seqs)
        rows = []
        for (bi, cond, rep), pred in zip(meta, preds):
            for ts in ag.tracksets:
                rows.append({"mpra_id": eid, "background": bi, "condition": cond, "rep": rep, "trackset": ts, "motif_source": source,
                             **{k: v for k, v in L.readout(pred[ts], reg, None, cfg).items() if k.startswith(("atac", "dnase"))}})
        L.append_rows(out, rows)
        n_done += 1
        if n_done % 10 == 0:
            print(f"{n_done} elements, {(time.time() - t0) / n_done:.1f} s/element", flush=True)
    print(f"done: {n_done} elements -> {out}")


if __name__ == "__main__":
    main()
