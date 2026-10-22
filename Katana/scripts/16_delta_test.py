#!/usr/bin/env python3
"""Small smoke test: can we read delta predictions (edited - native) from AlphaGenome on mm10 windows?

For a few AP-1 sites, runs the native 131 kb window and a handful of edits in one batch, then reports for every
trackset: native value, edited value, delta and log2FC at the motif / 170 bp / peak windows (ATAC) and the RNA readouts.
Edits: core@motif (AP-1 core destroyed), shuffle@motif, shuffle@170, plus offset@+600 (170 bp shuffle away from the
motif = negative control; should give a small delta) and a native-vs-native repeat (must be exactly 0).
Also saves the full ATAC delta profile (+-PROFILE_BP around the motif, 1 bp) so the shape of the response can be plotted.

  python 16_delta_test.py                       # 5 random CM sites, defaults
  python 16_delta_test.py --sites ../data/fb_fos_sites.csv --n 10 --reps 3
Outputs: ../results/delta_test/delta_summary.csv, delta_profiles.npz
"""
from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L

PROFILE_BP = 2000
EDITS = [("core@motif", "core", "motif"), ("shuffle@motif", "shuffle", "motif"),
         ("shuffle@170", "shuffle", "170"), ("offset@+600", "offset", "600")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default=str(L.ROOT / "data" / "cm_fos_sites.csv"))
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--reps", type=int, default=2, help="replicates per stochastic edit")
    ap.add_argument("--out", default=str(L.ROOT / "results" / "delta_test"))
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    tss_index = L.TssIndex(cfg["paths"]["gtf"])
    ag = L.AG(cfg)
    seed = cfg["perturb"]["seed"]

    sites = pd.read_csv(args.sites).sample(n=args.n, random_state=seed)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows, profiles = [], {}
    for _, row in sites.iterrows():
        sid = row["site_id"]
        codes, wstart = genome.window(row["seqnames"], row["motif_mid"], size)
        if codes is None:
            print(f"skip {sid}: off chromosome")
            continue
        reg = L.site_regions(row, wstart, size, cfg)
        if L.has_n(codes, *reg["motif"]):
            print(f"skip {sid}: N in motif")
            continue
        tss = tss_index.nearest(row["seqnames"], int(row["motif_mid"]))
        tss_rel = None if tss is None else tss - wstart
        rng = np.random.default_rng([seed, zlib.crc32(sid.encode())])

        seqs, meta = [codes, codes.copy()], [("native", 0), ("native_repeat", 0)]
        for name, kind, width in EDITS:
            for r in range(args.reps if kind != "core" else 1):
                e = L.apply_condition(codes, reg, kind, width, rng)
                if e is not None:
                    seqs.append(e)
                    meta.append((name, r))
                    assert (e != codes).any(), f"{name} did not change the sequence"
        preds = ag.predict(seqs)
        print(f"{sid} ({row['stable_dynamic']}) motif {L.decode(codes[slice(*reg['motif'])])}: "
              f"{len(seqs)} sequences, {sum((s != codes).sum() > 0 for s in seqs)} differ from native")

        for ts in ag.tracksets:
            nat = L.readout(preds[0][ts], reg, tss_rel, cfg)
            for (name, rep), s, pred in zip(meta, seqs, preds):
                if name == "native":
                    continue
                ed = L.readout(pred[ts], reg, tss_rel, cfg)
                for k, v in nat.items():
                    rows.append({"site_id": sid, "stable_dynamic": row["stable_dynamic"], "trackset": ts,
                                 "condition": name, "rep": rep, "n_bases_changed": int((s != codes).sum()),
                                 "readout": k, "native": v, "edited": ed[k], "delta": ed[k] - v,
                                 "log2fc": float(np.log2((ed[k] + 1e-6) / (v + 1e-6)))})
                if "atac" in pred[ts]:  # 1 bp ATAC delta profile around the motif
                    m = reg["mid"]
                    a, b = max(m - PROFILE_BP, 0), min(m + PROFILE_BP, size)
                    profiles[f"{sid}|{ts}|{name}|{rep}"] = (pred[ts]["atac"][a:b] - preds[0][ts]["atac"][a:b]).astype(np.float32)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "delta_summary.csv", index=False)
    np.savez_compressed(out_dir / "delta_profiles.npz", **profiles)

    pd.set_option("display.width", 200)
    print("\n== sanity: native_repeat delta must be 0 ==")
    print(df[df.condition == "native_repeat"].delta.abs().max())
    print("\n== mean ATAC log2FC at the motif by condition and trackset (offset@+600 = negative control) ==")
    acc = df[df.readout == "atac_motif"]
    print(acc.groupby(["trackset", "condition"]).agg(n=("log2fc", "size"), native=("native", "mean"),
                                                     mean_delta=("delta", "mean"), mean_log2fc=("log2fc", "mean")).round(3))
    print(f"\nwrote {out_dir}/delta_summary.csv and delta_profiles.npz")


if __name__ == "__main__":
    main()
