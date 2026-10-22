#!/usr/bin/env python3
"""Native (unedited) predictions of the fine-tuned ATAC head at every cardiomyocyte and fibroblast AP-1 site, with the built-in heads of the same
forward pass, so the head can be judged at our sites before the perturbation experiments are rerun with it.

Per site and trackset: the usual readouts (ATAC over the motif / 170 bp / peak, RNA), with per-condition columns for the fine-tuned head
(atac_peak_YS, ...). Per site table and track: the summed 1 bp prediction +-300 bp around the motif midpoints (aggregate footprint/flank shape).
Chunk 0 also checks that ag_lib's single-pass use of the head (built-in heads + new head on the same embeddings) equals the fine-tuned model's
own forward pass on a few windows.

  python 17_ft_evaluate.py --config ../config_ft.json --out ../results_ft/eval --chunk 0 --n-chunks 4 [--checkpoint X.delta.pth] [--limit 20]
Outputs: <out>_chunk<k>.csv, <out>_profiles_chunk<k>.npz, <out>_equivalence.json (chunk 0). Resumable: sites already in the csv are skipped
(their profiles are then missing from a rerun's npz, so rerun a chunk from scratch if the profiles matter).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L

HALF = 300


def equivalence(ag, cfg, windows):
    """Max relative difference between ag_lib's reading of the fine-tuned head and the fine-tuned model's own forward pass."""
    import torch
    from alphagenome_pytorch.extensions.finetuning.checkpointing import load_finetuned_model
    m, meta = load_finetuned_model(cfg["finetuned"]["checkpoint"], cfg["paths"]["weights"], device=ag.device)
    m.eval()
    x = torch.from_numpy(np.stack([L.onehot(w) for w in windows])).to(ag.device)
    _, ours = ag.forward_all(x)
    with torch.no_grad():
        theirs = m.predict(x, organism_index=ag.default_org, heads=("atac",), resolutions=(1, 128))["atac"]
    res = 1 if 1 in theirs else 128
    t = theirs[res].float()
    if res == 128:
        t = torch.repeat_interleave(t, 128, dim=1)
    diff = (ours - t).abs().max().item()
    scale = t.abs().max().item()
    del m
    torch.cuda.empty_cache()
    return {"resolution_compared": res, "max_abs_diff": diff, "max_abs_value": scale, "relative": diff / max(scale, 1e-12),
            "mean_ours": ours.mean().item(), "mean_theirs": t.mean().item(), "track_names": meta["track_names"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--checkpoint", default=None, help="override finetuned.checkpoint (or set AG_FT_CHECKPOINT), e.g. a pilot checkpoint for a smoke test")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="only the first N sites of this chunk")
    ap.add_argument("--group", type=int, default=16, help="sites per predict call")
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    if args.checkpoint:
        cfg["finetuned"]["checkpoint"] = args.checkpoint
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    tss_index = L.TssIndex(cfg["paths"]["gtf"])
    ag = L.AG(cfg)

    sites = pd.concat([pd.read_csv(L.ROOT / "data" / f"{t}_fos_sites.csv").assign(site_table=t) for t in ("cm", "fb")], ignore_index=True)
    sites = sites.iloc[args.chunk::args.n_chunks]
    if args.limit:
        sites = sites.head(args.limit)
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    done = L.done_ids(out, "row_id")
    t0 = time.time()

    if args.chunk == 0 and not (Path(f"{args.out}_equivalence.json").exists()):
        wins = [genome.window(r.seqnames, r.motif_mid, size)[0] for r in sites.head(2).itertuples()]
        eq = equivalence(ag, cfg, [w for w in wins if w is not None])
        Path(f"{args.out}_equivalence.json").write_text(json.dumps(eq, indent=1))
        print("equivalence:", json.dumps(eq), flush=True)

    prof, nprof = {}, {}
    todo = [r for r in sites.itertuples() if f"{r.site_table}:{r.site_id}" not in done]
    for g in range(0, len(todo), args.group):
        batch, seqs = [], []
        for r in todo[g:g + args.group]:
            codes, ws = genome.window(r.seqnames, r.motif_mid, size)
            if codes is None:
                continue
            batch.append((r, ws, L.site_regions(r._asdict(), ws, size, cfg)))
            seqs.append(codes)
        if not seqs:
            continue
        preds = ag.predict(seqs)
        rows = []
        for (r, ws, reg), pred in zip(batch, preds):
            tss = tss_index.nearest(r.seqnames, int(r.motif_mid))
            tss_rel = None if tss is None else tss - ws
            mid = reg["mid"]
            for ts in ag.tracksets:
                rows.append({"row_id": f"{r.site_table}:{r.site_id}", "site_table": r.site_table, "site_id": r.site_id, "binding_category": r.binding_category,
                             "stable_dynamic": r.stable_dynamic, "trackset": ts, **L.readout(pred[ts], reg, tss_rel, cfg)})
                if mid - HALF >= 0 and mid + HALF + 1 <= size:
                    for key, arr in pred[ts].items():
                        if key.startswith("atac"):
                            k = (r.site_table, ts, key)
                            prof[k] = prof.get(k, 0) + arr[mid - HALF:mid + HALF + 1].astype(np.float64)
                            nprof[k] = nprof.get(k, 0) + 1
        L.append_rows(out, rows)
        print(f"{min(g + args.group, len(todo))}/{len(todo)} sites, {(time.time() - t0) / max(g + args.group, 1):.2f} s/site", flush=True)
    np.savez(f"{args.out}_profiles_chunk{args.chunk}.npz", **{f"{t}|{ts}|{k}": v for (t, ts, k), v in prof.items()},
             **{f"n|{t}|{ts}|{k}": np.array(n) for (t, ts, k), n in nprof.items()})
    print(f"done: {len(todo)} sites in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
