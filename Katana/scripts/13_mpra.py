#!/usr/bin/env python3
"""Fibroblast MPRA validation: for every MPRA element (170 bp, mm10) predict ATAC/RNA for the native sequence, with the
AP-1 motif shuffled, and with the whole element shuffled, in the element's endogenous 131 kb window.

Motif location: the TOBIAS fibroblast Fos motif overlapping the element (data/fb_mpra_tobias_crossmap.csv; closest to the
element centre if several); otherwise the first TGA[CG]TCA hit in the element (11 bp around it); otherwise none (motif
shuffle skipped). `motif_source` records which. Compare predicted log2(shuffled/native) with the MPRA
delta_vs_MotifShuff / delta_vs_FullShuff columns afterwards.

  python 13_mpra.py --out ../results/mpra --chunk 0 --n-chunks 4
"""
from __future__ import annotations

import argparse
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    mc = cfg["mpra"]
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    tss_index = L.TssIndex(cfg["paths"]["gtf"])
    ag = L.AG(cfg)

    mpra = pd.read_csv(L.ROOT / "data" / "MPRA_Results_with_TSS_distance.csv")
    xm = pd.read_csv(L.ROOT / "data" / "fb_mpra_tobias_crossmap.csv")
    xm = xm[xm.overlap_status == "Matched"]
    mpra = mpra.drop_duplicates("Coords").reset_index(drop=True)  # 359 elements are listed under two MPRA groups (identical values)
    mpra = mpra.iloc[args.chunk::args.n_chunks]
    if args.limit:  # pilot: random subset
        mpra = mpra.sample(n=min(args.limit, len(mpra)), random_state=mc["seed"])
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    skip = L.done_ids(out, "mpra_id")
    t0, n_done = time.time(), 0
    for _, e in mpra.iterrows():
        eid = e["Coords"]
        if eid in skip:
            continue
        mid = (int(e["start"]) + int(e["end"])) // 2
        codes, ws = genome.window(e["chr"], mid, size)
        if codes is None:
            continue
        a, b = int(e["start"]) - ws, int(e["end"]) - ws
        motif, source = None, "none"
        hits = xm[xm.mpra_id == eid]
        if len(hits):
            h = hits.iloc[np.argmin(np.abs((hits.motif_start + hits.motif_end) / 2 - mid))]
            motif, source = (int(h.motif_start) - ws, int(h.motif_end) - ws), "tobias"
        else:
            m = L.CORE_RE.search(L.decode(codes[a:b]))
            if m:
                motif, source = (a + m.start() - 2, a + m.end() + 2), "regex"
        reg = {"motif": motif if motif else (a, a), "w170": (a, b), "peak": (a, b), "mid": mid - ws}
        tss = tss_index.nearest(e["chr"], mid)
        tss_rel = None if tss is None else tss - ws
        rng = np.random.default_rng([mc["seed"], zlib.crc32(eid.encode())])

        seqs, meta = [codes], [("native", 0)]
        if motif and not L.has_n(codes, *reg["motif"]):
            for r in range(mc["reps_motif_shuffle"]):
                seqs.append(L.shuffle_region(codes, *reg["motif"], rng))
                meta.append(("motif_shuffle", r))
        if not L.has_n(codes, a, b):
            for r in range(mc["reps_full_shuffle"]):
                seqs.append(L.shuffle_region(codes, a, b, rng))
                meta.append(("full_shuffle", r))
        preds = ag.predict(seqs)
        rows = []
        for (cond, rep), pred in zip(meta, preds):
            for ts in ag.tracksets:
                rows.append({"mpra_id": eid, "condition": cond, "rep": rep, "trackset": ts, "motif_source": source,
                             **L.readout(pred[ts], reg, tss_rel, cfg)})
        L.append_rows(out, rows)
        n_done += 1
        if n_done % 25 == 0:
            print(f"{n_done} elements, {(time.time() - t0) / n_done:.1f} s/element", flush=True)
    print(f"done: {n_done} elements -> {out}")


if __name__ == "__main__":
    main()
