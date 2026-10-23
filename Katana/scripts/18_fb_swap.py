#!/usr/bin/env python3
"""Fibroblast element swaps: does an element carry its accessibility with it, or does the location set it?

Sites: data/fb_swap_sites.csv, 50 per binding class (swap_class: Stable, Age Up, Age Down, Exercise Up, Exercise Down;
labels straight from the site table). Every site is used both as a donor (its sequence) and as a recipient (its location),
so each class pair is covered in both directions and the diagonal (a site put back into itself) is a control that must
reproduce the native prediction.

  --experiment swap        every donor element into every recipient location, at two widths  -> 2 x N x N predictions
      170     the donor's 170 bp around its motif replaces the recipient's 170 bp (motif midpoints aligned; equal widths, a clean swap)
      peak    the recipient's whole peak is removed and the donor's peak written centred on the recipient peak midpoint; any
              recipient-peak base the donor does not cover is filled with a dinucleotide-shuffled copy of itself, so no
              recipient peak sequence is left to contribute
  --experiment background  every element into quiet host windows (closed in every heart cell type; --backgrounds), at both
      widths, with a dinucleotide-shuffled copy of the element as the control  -> 4 x N x B predictions
      This asks whether the element's own sequence carries its accessibility with no regulatory context at all.

Readouts come from one forward pass per edited window, for every trackset of the config (the fine-tuned fibroblast head
averaged over the four conditions is the comparison track; the per-condition columns and the built-in head come free):
  atac_insert*  over the inserted element, atac_170* / atac_peak*  over the recipient's (or the host centre's) 170 bp and peak.
Output is long format, one row per pair x width x trackset, with the donor's and the recipient's own native values on the row
so the analysis can regress swapped on donor and recipient without a second lookup.
Resumable: recipients (swap) or elements (background) already finished are skipped.

  python 18_fb_swap.py --experiment swap --out ../results_ft/fbswap --chunk 0 --n-chunks 40
  python 18_fb_swap.py --experiment background --backgrounds $FT/closed_backgrounds.csv --out ../results_ft/fbbg --chunk 0 --n-chunks 8
  python 18_fb_swap.py --selftest      # no GPU, no genome: checks every edit changes exactly the bases it should
"""
from __future__ import annotations

import argparse
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L

WIDTHS = ("170", "peak")


def element(codes, reg, width):
    """(start, end) of the element to take from a donor window, or to replace in a recipient window."""
    return reg["w170"] if width == "170" else reg["peak"]


def swap_edit(codes, reg, donor_seq, width, rng):
    """Write donor_seq into this window. 170: over the recipient's own 170 bp. peak: centred on the recipient peak midpoint,
    with any remaining recipient-peak base replaced by a dinucleotide-shuffled copy, so none of its peak sequence survives.
    Returns (edited window, (start, end) of the inserted element) or (None, None) if it does not fit."""
    w = len(donor_seq)
    if width == "170":
        a, b = reg["w170"]
        if b - a != w:
            return None, None
        e = codes.copy()
        e[a:b] = donor_seq
        return e, (a, b)
    p0, p1 = reg["peak"]
    mid = (p0 + p1) // 2
    a = mid - w // 2
    b = a + w
    if a < 0 or b > len(codes):
        return None, None
    e = codes.copy()
    if p0 < a:                                   # recipient peak bases the donor does not cover: shuffle them away
        e = L.dinuc_region(e, p0, a, rng)
    if b < p1:
        e = L.dinuc_region(e, b, p1, rng)
    e[a:b] = donor_seq
    return e, (a, b)


def insert_into_background(codes, seq, rng):
    """Write seq at the centre of a host window. Returns (edited window, (start, end))."""
    w = len(seq)
    a = len(codes) // 2 - w // 2
    b = a + w
    if a < 0 or b > len(codes):
        return None, None
    e = codes.copy()
    e[a:b] = seq
    return e, (a, b)


def read_row(pred_ts, insert, regs, cfg):
    """Readouts with the inserted element in place of the motif window: atac_insert* / atac_170* / atac_peak*."""
    r = L.readout(pred_ts, {"motif": insert, "w170": regs["w170"], "peak": regs["peak"], "mid": regs["mid"]}, None, cfg)
    return {(k.replace("atac_motif", "atac_insert") if k.startswith("atac_motif") else k): v for k, v in r.items()}


def selftest():
    rng = np.random.default_rng(18)
    size = 131072
    codes = rng.integers(0, 4, size=size).astype(np.uint8)
    mid = size // 2
    reg = {"motif": (mid - 5, mid + 6), "mid": mid, "w170": (mid - 85, mid + 85), "peak": (mid - 400, mid + 400)}
    donor170 = rng.integers(0, 4, size=170).astype(np.uint8)

    e, ins = swap_edit(codes, reg, donor170, "170", rng)
    assert ins == reg["w170"] and (e[ins[0]:ins[1]] == donor170).all()
    assert (e[:ins[0]] == codes[:ins[0]]).all() and (e[ins[1]:] == codes[ins[1]:]).all(), "170 swap must touch only the 170 bp"

    for w in (300, 800, 1200):                                   # donor narrower, equal, wider than the 800 bp recipient peak
        donor = rng.integers(0, 4, size=w).astype(np.uint8)
        e, ins = swap_edit(codes, reg, donor, "peak", rng)
        assert (e[ins[0]:ins[1]] == donor).all(), w
        p0, p1 = reg["peak"]
        left, right = e[p0:ins[0]], e[ins[1]:p1]                 # leftover recipient peak must be shuffled, not original
        for seg, a0, b0 in ((left, p0, ins[0]), (right, ins[1], p1)):
            if len(seg):
                assert (seg != codes[a0:b0]).any(), f"leftover recipient peak not shuffled ({w})"
                assert sorted(seg) == sorted(codes[a0:b0]), f"shuffle must keep composition ({w})"
        outside = np.concatenate([np.arange(0, min(p0, ins[0])), np.arange(max(p1, ins[1]), size)])
        assert (e[outside] == codes[outside]).all(), f"edit must stay inside the peak/insert span ({w})"

    e, ins = insert_into_background(codes, donor170, rng)
    assert ins == (mid - 85, mid + 85) and (e[ins[0]:ins[1]] == donor170).all()
    print(f"selftest OK: 170 swap, peak swap (donor narrower/equal/wider), background insertion")


def load_sites(cfg, genome, size, sites):
    """Native window, regions and element sequences for every site."""
    out = {}
    for r in sites.itertuples():
        codes, ws = genome.window(r.seqnames, r.motif_mid, size)
        if codes is None:
            continue
        reg = L.site_regions(r._asdict(), ws, size, cfg)
        if L.has_n(codes, *reg["peak"]) or L.has_n(codes, *reg["w170"]):
            continue
        els = {w: codes[slice(*element(codes, reg, w))].copy() for w in WIDTHS}
        out[r.site_id] = {"codes": codes, "reg": reg, "class": r.swap_class, "els": els}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment", choices=["swap", "background"])
    ap.add_argument("--sites", default="../data/fb_swap_sites.csv")
    ap.add_argument("--backgrounds", default=None, help="closed_backgrounds.csv (background experiment)")
    ap.add_argument("--out")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="only the first N sites (pilot)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    assert args.experiment and args.out, "--experiment and --out are required"

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    ag = L.AG(cfg)
    sites = pd.read_csv(args.sites)
    if args.limit:
        sites = sites.head(args.limit)
    S = load_sites(cfg, genome, size, sites)
    print(f"{len(S)} usable sites of {len(sites)}", flush=True)

    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    t0 = time.time()

    # native prediction of every site: the donor's own level, and the value a self-swap must reproduce
    nat = {}
    ids = sorted(S)
    for i in range(0, len(ids), 8):
        part = ids[i:i + 8]
        preds = ag.predict([S[s]["codes"] for s in part])
        for s, p in zip(part, preds):
            reg = S[s]["reg"]
            nat[s] = {ts: read_row(p[ts], reg["w170"], reg, cfg) for ts in ag.tracksets}
    print(f"native predictions done in {time.time() - t0:.0f}s", flush=True)

    if args.experiment == "swap":
        recips = [s for k, s in enumerate(ids) if k % args.n_chunks == args.chunk]
        done = L.done_ids(out, "recipient_id")
        for n, rid in enumerate(recips):
            if rid in done:
                continue
            R = S[rid]
            rng = np.random.default_rng([cfg["perturb"]["seed"], zlib.crc32(rid.encode())])
            seqs, meta = [], []
            for did in ids:
                for w in WIDTHS:
                    e, ins = swap_edit(R["codes"], R["reg"], S[did]["els"][w], w, rng)
                    if e is not None:
                        seqs.append(e)
                        meta.append((did, w, ins))
            preds = ag.predict(seqs)
            rows = []
            for (did, w, ins), pred in zip(meta, preds):
                for ts in ag.tracksets:
                    rows.append({"recipient_id": rid, "recipient_class": R["class"], "donor_id": did, "donor_class": S[did]["class"],
                                 "width": w, "self_swap": did == rid, "insert_bp": ins[1] - ins[0], "trackset": ts,
                                 **read_row(pred[ts], ins, R["reg"], cfg),
                                 **{f"donor_native_{k}": v for k, v in nat[did][ts].items()},
                                 **{f"recipient_native_{k}": v for k, v in nat[rid][ts].items()}})
            L.append_rows(out, rows)
            print(f"{n + 1}/{len(recips)} recipients, {(time.time() - t0) / (n + 1):.1f} s/recipient", flush=True)
    else:
        bg = pd.read_csv(args.backgrounds)
        hosts = {}
        for b in bg.itertuples():
            codes, ws = genome.window(b.chrom, b.centre, size)
            if codes is not None and not L.has_n(codes, 0, size):
                hosts[f"{b.chrom}:{b.centre}"] = codes
        print(f"{len(hosts)} host windows", flush=True)
        els = [s for k, s in enumerate(ids) if k % args.n_chunks == args.chunk]
        done = L.done_ids(out, "donor_id")
        for n, did in enumerate(els):
            if did in done:
                continue
            rng = np.random.default_rng([cfg["perturb"]["seed"], zlib.crc32(did.encode())])
            seqs, meta = [], []
            for hid, hcodes in hosts.items():
                for w in WIDTHS:
                    el = S[did]["els"][w]
                    for kind, seq in (("element", el), ("shuffled", L.dinuc_shuffle_codes(el, rng))):
                        e, ins = insert_into_background(hcodes, seq, rng)
                        if e is not None:
                            seqs.append(e)
                            meta.append((hid, w, kind, ins))
            preds = ag.predict(seqs)
            rows = []
            for (hid, w, kind, ins), pred in zip(meta, preds):
                mid = size // 2
                regs = {"w170": (mid - 85, mid + 85), "peak": (ins[0], ins[1]), "mid": mid}
                for ts in ag.tracksets:
                    rows.append({"donor_id": did, "donor_class": S[did]["class"], "host_id": hid, "width": w, "kind": kind,
                                 "insert_bp": ins[1] - ins[0], "trackset": ts, **read_row(pred[ts], ins, regs, cfg),
                                 **{f"donor_native_{k}": v for k, v in nat[did][ts].items()}})
            L.append_rows(out, rows)
            print(f"{n + 1}/{len(els)} elements, {(time.time() - t0) / (n + 1):.1f} s/element", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
