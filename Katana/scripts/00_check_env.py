#!/usr/bin/env python3
"""Environment check (run first, on a GPU node): imports, file paths, weights load, mouse track list, track matching,
and a timing of one forward pass. Writes results/mouse_tracks.csv so tracksets in config.json can be chosen.

  python 00_check_env.py
"""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

import ag_lib as L


def main():
    cfg = L.load_config()
    for k, p in cfg["paths"].items():
        print(f"{'OK     ' if Path(p).exists() else 'MISSING'} {k}: {p}")
    import torch
    print("torch", torch.__version__, "cuda", torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
    g = L.Genome(cfg["paths"]["genome"])
    print("chr1 length", g.length("chr1"), "| chrX in genome:", "chrX" in g.fa)
    for t in ("cm_fos_sites", "fb_fos_sites"):
        df = pd.read_csv(L.ROOT / "data" / f"{t}.csv", nrows=3)
        w, s = g.window(df.seqnames[0], df.motif_mid[0], cfg["model"]["window"])
        print(t, "first site window ok:", w is not None, "| motif:", L.decode(g.fetch(df.seqnames[0], df.start[0], df.end[0])))
    tss = L.TssIndex(cfg["paths"]["gtf"])
    print("GTF chromosomes:", len(tss.by_chrom), "| nearest TSS to chr1:7088913 ->", tss.nearest("chr1", 7088913))

    ag = L.AG(cfg, strict=False)
    device = ag.device
    rows = ag.track_table()
    out = L.ROOT / "results" / "mouse_tracks.csv"
    out.parent.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(df.groupby("head").size().to_string())
    print(f"all {len(df)} tracks written to {out}")
    print("heads computed in predictions:", ag.heads)
    try:  # human tracks too, to look for more cell-type-specific cardiac / fibroblast tracks (applied to mouse sequence)
        hdf = pd.DataFrame(L.catalog_tracks(0))
        hdf.to_csv(L.ROOT / "results" / "human_tracks.csv", index=False)
        print(hdf.groupby("head").size().to_string())
        print(f"{len(hdf)} human tracks written to results/human_tracks.csv")
    except Exception as e:
        print("human track listing failed:", type(e).__name__, str(e)[:300])

    x = torch.zeros((cfg["model"]["batch_size"], ag.size, 4), device=device)
    for _ in range(2):
        t0 = time.time()
        ag._named(x)  # same heads as the real runs
        if torch.cuda.is_available():
            torch.cuda.synchronize()
        dt = time.time() - t0
    print(f"forward pass, batch {x.shape[0]}: {dt:.2f}s ({dt / x.shape[0]:.3f} s/sequence)")
    conds = L.build_conditions(cfg['perturb']['reps'], cfg['perturb'].get('reps_width'))
    per_site = 1 + sum(c[3] for c in conds)
    print(f"-> ~{per_site * dt / x.shape[0]:.0f} s per site at {per_site} passes/site; "
          f"CM+FB (8,697 sites) ~{8697 * per_site * dt / x.shape[0] / 3600:.0f} GPU-hours in total")


if __name__ == "__main__":
    main()
