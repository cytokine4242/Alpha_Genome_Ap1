#!/usr/bin/env python3
"""The host windows with nothing written into them: are the backgrounds actually closed, and what is the real baseline?

18_fb_swap.py --experiment background scores every element in 20 quiet host windows against a dinucleotide-shuffled
copy of itself in the same host. That contrast is well controlled for base composition, but it never predicts the
host untouched, so two things were never checked:

  1  are the hosts really inactive? They were chosen as closed in every heart cell type from the DATA. Whether the
     MODEL also calls them closed is a separate question, and the whole "builds a site from nothing" claim rests on it.
  2  how much of the lift is the element, and how much is simply that a scrambled insert is not neutral? If a shuffled
     insert already opens the host, element-minus-shuffled understates what inserting the element does; if a shuffled
     insert is repressive, it overstates it.

This predicts each host window once, untouched, and reads it out over exactly the spans the insertions used, so the
no-insert baseline is matched to the inserted rows window for window. 20 predictions; the cost is in the readouts.

Spans: 170 bp, plus every distinct peak-insert width present in the existing background results (--spans-from), so
every (host_id, insert_bp) pair in those results gets a baseline row to join against.

Outputs <out>.csv, long format: host_id x insert_bp x trackset.
  python 20_host_native.py --backgrounds $FT/closed_backgrounds.csv --out ../results_ft/hostnat_fb
"""
from __future__ import annotations

import argparse
import glob
import time
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backgrounds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--spans-from", default="../results_ft/bg_fb_chunk*.csv",
                    help="existing background results, to take the insert widths from")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    ag = L.AG(cfg)
    t0 = time.time()

    bg = pd.read_csv(args.backgrounds)
    hosts = {}
    for b in bg.itertuples():
        codes, _ = genome.window(b.chrom, b.centre, size)
        if codes is not None and not L.has_n(codes, 0, size):
            hosts[f"{b.chrom}:{b.centre}"] = codes
    print(f"{len(hosts)} host windows", flush=True)

    # the spans the insertions actually occupied, so every inserted row has a matching no-insert row
    files = sorted(glob.glob(args.spans_from))
    if files:
        prev = pd.concat([pd.read_csv(f, usecols=["insert_bp"]) for f in files], ignore_index=True)
        spans = sorted(int(x) for x in prev.insert_bp.unique())
    else:
        spans = [170, 800]
        print(f"no results matched {args.spans_from}; falling back to {spans}", flush=True)
    print(f"{len(spans)} distinct insert widths, {spans[0]} to {spans[-1]} bp", flush=True)

    preds = ag.predict(list(hosts.values()))
    print(f"predictions done in {time.time() - t0:.0f}s", flush=True)

    mid = size // 2
    rows = []
    for hid, pred in zip(hosts.keys(), preds):
        for w in spans:
            a = mid - w // 2
            b = a + w
            if a < 0 or b > size:
                continue
            regs = {"motif": (a, b), "w170": (mid - 85, mid + 85), "peak": (a, b), "mid": mid}
            for ts in ag.tracksets:
                r = L.readout(pred[ts], regs, None, cfg)
                rows.append({"host_id": hid, "insert_bp": w, "trackset": ts,
                             **{(k.replace("atac_motif", "atac_insert") if k.startswith("atac_motif") else k): v
                                for k, v in r.items()}})
    out = Path(f"{args.out}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"done in {time.time() - t0:.0f}s -> {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
