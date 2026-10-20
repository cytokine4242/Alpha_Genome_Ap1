#!/usr/bin/env python3
"""Insertion experiments (slides: peak swap and AP-1 motif insertion).

A. Peak swap: a few CM recipient windows (chosen across the range of native peak ATAC) have their peak replaced by the
   peak sequence of many donor sites (stable + dynamic). Swapped ATAC over the inserted interval is compared with the
   donor's own native ATAC (sequence) and the recipient's native ATAC (location). The inserted sequence is centred on the
   recipient peak midpoint and has the donor peak's width.
B. Motif insertion: random mm10 backgrounds are dinucleotide-shuffled, then TGAGTCA is inserted at the window centre
   (control: shuffled motif). Variants per background: native, shuffled background, shuffled+motif, shuffled+shuffled-motif.
   ATAC readouts: atac_motif = centre 7 bp, atac_peak = centre 128 bp, atac_170 = centre 170 bp.

  python 12_insertion.py --out ../results/insertion [--part swap|motif|both] [--limit N]
"""
from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L


def peak_swap(cfg, genome, ag, out, limit):
    ic, pc = cfg["insertion"], cfg["perturb"]
    size = cfg["model"]["window"]
    ts = ic["trackset"]
    rng = np.random.default_rng(ic["seed"])
    sites = pd.read_csv(L.ROOT / "data" / "cm_fos_sites.csv")
    sites = sites[(sites.peak_end - sites.peak_start) <= pc["max_peak_width"]].reset_index(drop=True)

    def native(row):
        codes, ws = genome.window(row["seqnames"], row["motif_mid"], size)
        if codes is None or L.has_n(codes, 0, size):
            return None
        reg = L.site_regions(row, ws, size, cfg)
        pred = ag.predict([codes])[0]
        return codes, ws, reg, L.readout(pred[ts], reg, None, cfg)["atac_peak"]

    pool = sites.sample(n=min(ic["recipient_pool"], len(sites)), random_state=ic["seed"])
    if limit:
        pool = pool.head(limit)
    nat = {}
    for _, r in pool.iterrows():
        res = native(r)
        if res:
            nat[r["site_id"]] = res
    ranked = sorted(nat, key=lambda k: nat[k][3])
    qs = ic["recipient_quantiles"][: ic["n_recipients"]]
    recips = [ranked[min(int(q * len(ranked)), len(ranked) - 1)] for q in qs]
    print("recipients (native peak ATAC):", [(r, round(nat[r][3], 3)) for r in recips], flush=True)

    rest = sites[~sites.site_id.isin(recips)]
    n_stable = ic["n_donors_stable"] if not limit else 2
    n_dyn = ic["n_donors_dynamic"] if not limit else 2
    donors = pd.concat([rest[rest.stable_dynamic == "Stable"].sample(n_stable, random_state=ic["seed"]),
                        rest[rest.stable_dynamic == "Dynamic"].sample(n_dyn, random_state=ic["seed"])])
    rows = []
    for _, d in donors.iterrows():
        dn = native(d)
        peak = genome.fetch(d["peak_chr"], int(d["peak_start"]), int(d["peak_end"]))
        if dn is None or peak is None or L.has_n(peak, 0, len(peak)):
            continue
        w = len(peak)
        for rid in recips:
            codes, ws, reg, rnat = nat[rid]
            row_r = sites.set_index("site_id").loc[rid]
            pm = (reg["peak"][0] + reg["peak"][1]) // 2
            a = pm - w // 2
            if a < 0 or a + w > size:
                continue
            e = codes.copy()
            e[a:a + w] = peak
            r2 = dict(reg, peak=(a, a + w))
            pred = ag.predict([e])[0]
            rows.append({"recipient_id": rid, "donor_id": d["site_id"], "donor_binding_category": d["binding_category"], "donor_age_class": d["age_class"],
                         "donor_ex_class": d["ex_class"],
                         "donor_stable_dynamic": d["stable_dynamic"], "donor_peak_width": w,
                         "recipient_native_atac": rnat, "donor_native_atac": dn[3],
                         "swapped_atac": L.readout(pred[ts], r2, None, cfg)["atac_peak"]})
    pd.DataFrame(rows).to_csv(f"{out}_peakswap.csv", index=False)
    print(f"peak swap: {len(rows)} donor x recipient pairs", flush=True)


def motif_insertion(cfg, genome, ag, out, limit):
    ic = cfg["insertion"]
    size = cfg["model"]["window"]
    rng = np.random.default_rng(ic["seed"])
    motif = L.encode(ic["insert_motif"])
    chroms = [f"chr{i}" for i in range(1, 20)] + ["chrX"]
    lens = np.array([genome.length(c) for c in chroms], dtype=float)
    c = size // 2
    reg = {"motif": (c - 3, c + 4), "w170": (c - 85, c + 85), "peak": (c - 64, c + 64), "mid": c}
    rows, n = [], 0
    n_target = 3 if limit else ic["n_backgrounds"]
    while n < n_target:
        chrom = chroms[rng.choice(len(chroms), p=lens / lens.sum())]
        mid = int(rng.integers(size, genome.length(chrom) - size))
        codes, ws = genome.window(chrom, mid, size)
        if codes is None or L.has_n(codes, 0, size):
            continue
        bg_rng = np.random.default_rng([ic["seed"], n])
        shuf = L.dinuc_shuffle_codes(codes, bg_rng)
        with_motif = shuf.copy()
        with_motif[c - 3:c + 4] = motif
        ctrl = motif.copy()
        while (ctrl == motif).all() or L.CORE_RE.search(L.decode(ctrl)):
            ctrl = bg_rng.permutation(motif)
        with_ctrl = shuf.copy()
        with_ctrl[c - 3:c + 4] = ctrl
        variants = {"background_native": codes, "background_shuffled": shuf,
                    "shuffled_plus_motif": with_motif, "shuffled_plus_shuffled_motif": with_ctrl}
        preds = ag.predict(list(variants.values()))
        for (name, _), pred in zip(variants.items(), preds):
            for t in ag.tracksets:
                rows.append({"background_id": f"{chrom}:{ws}", "variant": name, "trackset": t,
                             **{k: v for k, v in L.readout(pred[t], reg, None, cfg).items() if k.startswith("atac")}})
        n += 1
    pd.DataFrame(rows).to_csv(f"{out}_motif.csv", index=False)
    print(f"motif insertion: {n} backgrounds", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--part", default="both", choices=["swap", "motif", "both"])
    ap.add_argument("--limit", type=int, default=0, help="pilot: tiny run")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = L.load_config(args.config)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    genome = L.Genome(cfg["paths"]["genome"])
    ag = L.AG(cfg)
    if args.part in ("swap", "both"):
        peak_swap(cfg, genome, ag, args.out, args.limit)
    if args.part in ("motif", "both"):
        motif_insertion(cfg, genome, ag, args.out, args.limit)


if __name__ == "__main__":
    main()
