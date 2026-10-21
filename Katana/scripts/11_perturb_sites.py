#!/usr/bin/env python3
"""Shuffle / random / dinucleotide / targeted-core / flank perturbations of endogenous AP-1 sites (slides: shuffle
and mutation arms), with ATAC (motif, 170 bp, peak) and RNA (proximal +-2 kb of nearest TSS, mean/max +-50 kb) readouts.

One 131 kb mm10 window per site, centred on the motif midpoint: native + every edit. All tracksets from config.json
are read out from the same forward pass. Output is long format (one row per site x condition x replicate x trackset);
log2FC and dependency are computed afterwards from the raw values. Resumable: sites already in the output are skipped.

  python 11_perturb_sites.py --sites ../data/cm_fos_sites.csv --out ../results/cm --chunk 0 --n-chunks 20
  python 11_perturb_sites.py --sites ../data/cm_fos_sites.csv --out ../results/pilot_cm --limit 20   # pilot
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
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True, help="output prefix; writes <prefix>_chunk<N>.csv")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="only the first N sites of this chunk (pilot)")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    tss_index = L.TssIndex(cfg["paths"]["gtf"])
    ag = L.AG(cfg)
    conds = L.build_conditions(cfg["perturb"]["reps"])

    sites = pd.read_csv(args.sites)
    sites = sites.iloc[args.chunk::args.n_chunks]
    if args.limit:  # pilot: random subset (the tables are sorted by chromosome)
        sites = sites.sample(n=min(args.limit, len(sites)), random_state=cfg["perturb"]["seed"])
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    skip = L.done_ids(out)
    log = Path(f"{args.out}_chunk{args.chunk}.skipped.txt")
    t0 = time.time()
    n_done = 0
    for _, row in sites.iterrows():
        sid = row["site_id"]
        if sid in skip:
            continue
        codes, wstart = genome.window(row["seqnames"], row["motif_mid"], size)
        if codes is None:
            log.open("a").write(f"{sid}\toff chromosome\n")
            continue
        reg = L.site_regions(row, wstart, size, cfg)
        if reg["motif"][1] - reg["motif"][0] < 5 or L.has_n(codes, *reg["motif"]):
            log.open("a").write(f"{sid}\tmotif has N or is outside window\n")
            continue
        tss = tss_index.nearest(row["seqnames"], int(row["motif_mid"]))
        tss_rel = None if tss is None else tss - wstart
        rng = np.random.default_rng([cfg["perturb"]["seed"], zlib.crc32(sid.encode())])

        seqs, meta = [codes], [("native", "native", "none", 0)]
        for name, kind, width, reps in conds:
            for r in range(reps):
                e = L.apply_condition(codes, reg, kind, width, rng)
                if e is not None:
                    seqs.append(e)
                    meta.append((name, kind, width, r))
        preds = ag.predict(seqs)

        rows = []
        for (name, kind, width, rep), pred in zip(meta, preds):
            for ts in ag.tracksets:
                rows.append({"site_id": sid, "binding_category": row["binding_category"], "age_class": row["age_class"],
                             "ex_class": row["ex_class"], "stable_dynamic": row["stable_dynamic"],
                             "condition": name, "kind": kind, "width": width, "rep": rep, "trackset": ts,
                             "tss_dist": None if tss is None else tss - int(row["motif_mid"]),
                             **L.readout(pred[ts], reg, tss_rel, cfg)})
        L.append_rows(out, rows)
        n_done += 1
        if n_done % 10 == 0:
            print(f"{n_done} sites, {(time.time() - t0) / n_done:.1f} s/site", flush=True)
    print(f"done: {n_done} sites in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
