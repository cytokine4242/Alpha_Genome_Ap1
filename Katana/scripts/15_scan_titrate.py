#!/usr/bin/env python3
"""Two experiments on a fixed cardiomyocyte site subset (data/cm_subset_dynamics.csv, 500 sites per group: Dynamic; context-matched Stable sites matched on AlphaGenome's native predicted accessibility, peak width and TSS distance; and high-accessibility Stable sites). No TOBIAS footprint scores are used.
Question: what is different about the sequence neighbourhood of Dynamic sites, and how does accessibility respond to AP-1 motif strength and position?

  --experiment scan      where do the important bases lie?
      tile scan: every base of a 10 bp tile replaced by a different random base, tiles from -300 to +290 bp around the motif midpoint (60 tiles x 3 replicates)
      saturation mutagenesis: each of the 3 alternative bases at every position from -30 to +30 bp (61 x 3 = 183 passes)     -> 364 passes per site

  --experiment titrate   how does accessibility respond to AP-1 motif strength and position?  (7 bp core = TGA[C/G]TCA, or the central 7 bp if the 11 bp motif holds no exact core)
      core_mm1..core_mm7  1 to 7 of the 7 core bases replaced (5 random draws each; 3 for 7), redrawn until the core is no longer a perfect TGA[C/G]TCA   strength: weakening
      core_upgrade        only at sites WITHOUT an exact core: the 7 bp window of the motif closest to TGA[C/G]TCA made perfect with the fewest changes (pos = bases changed)   strength: upgrade
                          (an exact core is already perfect: TGACTCA and TGAGTCA are each other's reverse complement, so there is nothing to upgrade)
      motif_rc            the whole 11 bp motif replaced by its reverse complement                            orientation of the motif relative to its surroundings (the 7 bp core is a palindrome,
                                                                                                              so reversing only the core would just swap its centre C/G)
      core_halfsite       first half-site TGA[C/G] kept, last 3 core bases replaced (3 draws)                 half-site only
      dup@o               a second copy of the 11 bp motif written at offset o = -150,-80,-40,-20,+20,+40,+80,+150 (original kept)
      ctrl@o              a composition-matched shuffled 11-mer written at the same offset (3 draws)        control for dup: overwriting the local sequence
      move@o              motif copied to offset o AND the original 11 bp replaced by a shuffled 11-mer (3 draws)
      ctrlmove@o          shuffled 11-mers at BOTH the offset and the original position (3 draws)           control for move   -> 118 passes per site

Readouts as in 11_perturb_sites.py (ATAC over the motif, the 170 bp and the peak; RNA) from one forward pass per edit, plus atac_target = mean ATAC over the 11 bp at the
offset for dup/ctrl/move. Output is long format (site x condition x replicate x trackset); log2FC against the native row is computed afterwards.
Resumable: sites already in the output are skipped.

  python 15_scan_titrate.py --experiment scan --sites ../data/cm_subset_dynamics.csv --out ../results/scan_cm --chunk 0 --n-chunks 40
  python 15_scan_titrate.py --experiment titrate --sites ../data/cm_subset_dynamics.csv --out ../results/titrate_cm --chunk 0 --n-chunks 40
  python 15_scan_titrate.py --selftest        # no GPU, no genome: checks that every edit changes exactly the bases it should
"""
from __future__ import annotations

import argparse
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L

TILE, TILE_HALF, ISM_HALF = 10, 300, 30
OFFSETS = (-150, -80, -40, -20, 20, 40, 80, 150)
CORE_LEN = 7


def _best_core(codes, motif):
    """(start, end, mismatches, perfect core) for the 7 bp window of the motif with the fewest mismatches to TGA[C/G]TCA, ties broken towards the motif centre.
    The perfect core keeps the window's centre base if it is C or G (the degenerate position), else uses G, so the upgrade changes as few bases as possible."""
    m0, m1 = motif
    centre = (m0 + m1 - CORE_LEN) / 2
    best = None
    for a in range(m0, m1 - CORE_LEN + 1):
        w = L.decode(codes[a:a + CORE_LEN])
        perfect = "TGA" + (w[3] if w[3] in "CG" else "G") + "TCA"
        mm = sum(x != y for x, y in zip(w, perfect))
        key = (mm, abs(a - centre))
        if best is None or key < best[0]:
            best = (key, a, perfect)
    (mm, _), a, perfect = best
    return a, a + CORE_LEN, mm, perfect


def _mutate(codes, idx, rng):
    """Replace the bases at idx by a different base each."""
    out = codes.copy()
    idx = np.asarray(idx)
    out[idx] = (out[idx] + rng.integers(1, 4, size=len(idx))).astype(np.uint8) % 4
    return out


def _clean_shuffle(motif, rng, left=None, right=None, tries=100):
    """A composition-matched permutation of the motif that differs from it and holds no AP-1 half-site, also not across the junction with the 3 flanking bases each side."""
    left = np.array([], dtype=np.uint8) if left is None else left
    right = np.array([], dtype=np.uint8) if right is None else right
    for _ in range(tries):
        s = rng.permutation(motif)
        txt = L.decode(np.concatenate([left, s, right]))
        if (s != motif).any() and not any(h in txt for h in L.HALF_SITES):
            return s
    return rng.permutation(motif)


def scan_conditions(codes, reg, rng, tile_reps=3):
    """[(name, kind, pos, alt, rep, edited window, target (a, b) or None)]"""
    mid, out = reg["mid"], []
    for o in range(-TILE_HALF, TILE_HALF, TILE):
        a, b = mid + o, mid + o + TILE
        if a < 0 or b > len(codes) or L.has_n(codes, a, b):
            continue
        for r in range(tile_reps):
            out.append((f"tile@{o:+d}", "tile", o, "", r, _mutate(codes, np.arange(a, b), rng), None))
    for p in range(-ISM_HALF, ISM_HALF + 1):
        i = mid + p
        if codes[i] > 3:
            continue
        for alt in range(4):
            if alt == codes[i]:
                continue
            e = codes.copy()
            e[i] = alt
            out.append((f"ism@{p:+d}{'ACGT'[alt]}", "ism", p, "ACGT"[alt], 0, e, None))
    return out


def titrate_conditions(codes, reg, rng):
    m0, m1 = reg["motif"]
    mid = reg["mid"]
    c0, c1, found = L.core_interval(codes, reg["motif"])
    motif = codes[m0:m1].copy()
    out = []
    for k in range(1, 8):
        for r in range(5 if k < 7 else 3):
            for _ in range(100):               # position 4 of the core is degenerate (C/G): redraw until the edited core is no longer a perfect TGA[C/G]TCA
                e = _mutate(codes, c0 + rng.choice(7, size=k, replace=False), rng)
                if not L.CORE_RE.search(L.decode(e[c0:c1])):
                    break
            out.append((f"core_mm{k}", "core_mm", k, "", r, e, None))
    if not found:                                # upgrade only where there is something to upgrade
        u0, u1, mm, perfect = _best_core(codes, reg["motif"])
        e = codes.copy()
        e[u0:u1] = L.encode(perfect)
        out.append(("core_upgrade", "core_upgrade", mm, "", 0, e, None))
    rc = L.COMPLEMENT[motif][::-1]
    if (motif != rc).any():
        e = codes.copy()
        e[m0:m1] = rc
        out.append(("motif_rc", "motif_rc", 0, "", 0, e, None))
    for r in range(3):
        out.append(("core_halfsite", "core_halfsite", 0, "", r, _mutate(codes, np.arange(c0 + 4, c1), rng), None))
    for o in OFFSETS:
        t0 = mid + o - (m1 - m0) // 2
        t1 = t0 + (m1 - m0)
        if t0 < 0 or t1 > len(codes) or L.has_n(codes, t0, t1):
            continue
        e = codes.copy()
        e[t0:t1] = motif
        out.append((f"dup@{o:+d}", "dup", o, "", 0, e, (t0, t1)))
        lt, rt = codes[max(t0 - 3, 0):t0], codes[t1:t1 + 3]
        lm, rm = codes[m0 - 3:m0], codes[m1:m1 + 3]
        for r in range(3):                       # control for dup: a shuffled 11-mer at the target, original motif kept
            e = codes.copy()
            e[t0:t1] = _clean_shuffle(motif, rng, lt, rt)
            out.append((f"ctrl@{o:+d}", "ctrl", o, "", r, e, (t0, t1)))
        for r in range(3):                       # move: motif at the target, original replaced by a shuffled 11-mer
            e = codes.copy()
            e[t0:t1] = motif
            e[m0:m1] = _clean_shuffle(motif, rng, lm, rm)
            out.append((f"move@{o:+d}", "move", o, "", r, e, (t0, t1)))
        for r in range(3):                       # control for move: shuffled 11-mers at BOTH the target and the original position
            e = codes.copy()
            e[t0:t1] = _clean_shuffle(motif, rng, lt, rt)
            e[m0:m1] = _clean_shuffle(motif, rng, lm, rm)
            out.append((f"ctrlmove@{o:+d}", "ctrlmove", o, "", r, e, (t0, t1)))
    return out, found


def selftest():
    rng = np.random.default_rng(5)
    size = 131072
    codes = rng.integers(0, 4, size=size).astype(np.uint8)
    mid = size // 2
    m0, m1 = mid - 5, mid + 6
    codes[m0 + 2:m0 + 9] = L.encode("TGAGTCA")                       # planted consensus core inside the 11 bp motif
    reg = {"motif": (m0, m1), "mid": mid, "w170": (mid - 85, mid + 85), "peak": (mid - 400, mid + 400)}
    c0, c1, found = L.core_interval(codes, reg["motif"])
    assert found and (c0, c1) == (m0 + 2, m0 + 9)

    sc = scan_conditions(codes, reg, rng)
    n_tile, n_ism = sum(c[1] == "tile" for c in sc), sum(c[1] == "ism" for c in sc)
    assert n_tile == 60 * 3 and n_ism == 61 * 3, (n_tile, n_ism)
    for name, kind, pos, alt, rep, e, tgt in sc:
        d = np.flatnonzero(e != codes)
        if kind == "tile":
            assert len(d) == TILE and d.min() == mid + pos and d.max() == mid + pos + TILE - 1, name
        else:
            assert len(d) == 1 and d[0] == mid + pos and e[d[0]] == "ACGT".index(alt), name

    tt, found = titrate_conditions(codes, reg, rng)
    names = [c[0] for c in tt]
    assert sum(n.startswith("core_mm") for n in names) == 6 * 5 + 3
    for name, kind, pos, alt, rep, e, tgt in tt:
        d = np.flatnonzero(e != codes)
        if kind == "core_mm":
            assert len(d) == pos and ((d >= c0) & (d < c1)).all() and not L.CORE_RE.search(L.decode(e[c0:c1])), name
        elif kind == "core_halfsite":
            assert len(d) == 3 and ((d >= c0 + 4) & (d < c1)).all(), name
        elif kind == "motif_rc":
            assert (e[m0:m1] == L.COMPLEMENT[codes[m0:m1]][::-1]).all() and set(d) <= set(range(m0, m1)), name
            assert L.CORE_RE.search(L.decode(e[m0:m1])), name             # a perfect core stays perfect when the motif is reversed
        elif kind == "dup":
            t0, t1 = tgt
            assert (e[t0:t1] == codes[m0:m1]).all() and (e[m0:m1] == codes[m0:m1]).all() and set(d) <= set(range(t0, t1)), name
        elif kind == "ctrlmove":
            t0, t1 = tgt
            assert sorted(e[t0:t1]) == sorted(codes[m0:m1]) and sorted(e[m0:m1]) == sorted(codes[m0:m1]) and (e[m0:m1] != codes[m0:m1]).any(), name
            assert not any(h in L.decode(e[t0 - 3:t1 + 3]) for h in L.HALF_SITES) and not any(h in L.decode(e[m0 - 3:m1 + 3]) for h in L.HALF_SITES), name
        elif kind == "ctrl":
            t0, t1 = tgt
            assert sorted(e[t0:t1]) == sorted(codes[m0:m1]) and (e[m0:m1] == codes[m0:m1]).all() and set(d) <= set(range(t0, t1)), name
            assert not any(h in L.decode(e[t0 - 3:t1 + 3]) for h in L.HALF_SITES), name
        elif kind == "move":
            t0, t1 = tgt
            assert (e[t0:t1] == codes[m0:m1]).all() and sorted(e[m0:m1]) == sorted(codes[m0:m1]) and (e[m0:m1] != codes[m0:m1]).any(), name
            assert set(d) <= set(range(m0, m1)) | set(range(t0, t1)), name
    assert "core_upgrade" not in names and "motif_rc" in names         # the planted core is already perfect
    assert sum(n.startswith("dup@") for n in names) == len(OFFSETS) and sum(n.startswith("move@") for n in names) == 3 * len(OFFSETS) and sum(n.startswith("ctrlmove@") for n in names) == 3 * len(OFFSETS)
    # a site without an exact core gets an upgrade with the fewest changes (centre C kept, the one wrong base fixed)
    codes2 = codes.copy()
    codes2[c0:c1] = L.encode("TGACGCA")
    tt2, found2 = titrate_conditions(codes2, reg, rng)
    up = [c for c in tt2 if c[0] == "core_upgrade"]
    assert not found2 and len(up) == 1 and up[0][2] == 1 and L.decode(up[0][5][c0:c1]) == "TGACTCA", up[0][:3]
    assert (np.flatnonzero(up[0][5] != codes2) == [c0 + 4]).all()
    # TGAATCA (wrong centre base): G is used, one change
    codes3 = codes.copy()
    codes3[c0:c1] = L.encode("TGAATCA")
    up3 = [c for c in titrate_conditions(codes3, reg, rng)[0] if c[0] == "core_upgrade"]
    assert len(up3) == 1 and up3[0][2] == 1 and L.decode(up3[0][5][c0:c1]) == "TGAGTCA"
    print(f"selftest OK: scan {len(sc)} conditions ({n_tile} tile, {n_ism} saturation), titrate {len(tt)} conditions, {len(tt) + len(sc) + 2} passes for both incl. natives")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment", choices=["scan", "titrate"])
    ap.add_argument("--sites", default="../data/cm_subset_dynamics.csv")
    ap.add_argument("--out")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="only a random N sites of this chunk (pilot)")
    ap.add_argument("--tile-reps", type=int, default=3)
    ap.add_argument("--config", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    assert args.experiment and args.out, "--experiment and --out are required"

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    tss_index = L.TssIndex(cfg["paths"]["gtf"])
    ag = L.AG(cfg)

    sites = pd.read_csv(args.sites).iloc[args.chunk::args.n_chunks]
    if args.limit:
        sites = sites.sample(n=min(args.limit, len(sites)), random_state=cfg["perturb"]["seed"])
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    skip = L.done_ids(out)
    log = Path(f"{args.out}_chunk{args.chunk}.skipped.txt")
    t0, n_done = time.time(), 0
    for _, row in sites.iterrows():
        sid = row["site_id"]
        if sid in skip:
            continue
        codes, wstart = genome.window(row["seqnames"], row["motif_mid"], size)
        if codes is None:
            log.open("a").write(f"{sid}\toff chromosome\n")
            continue
        reg = L.site_regions(row, wstart, size, cfg)
        if reg["motif"][1] - reg["motif"][0] != 11 or L.has_n(codes, reg["mid"] - TILE_HALF, reg["mid"] + TILE_HALF) or reg["mid"] - TILE_HALF < 0:
            log.open("a").write(f"{sid}\tmotif not 11 bp, N in the scanned region, or outside the window\n")
            continue
        tss = tss_index.nearest(row["seqnames"], int(row["motif_mid"]))
        tss_rel = None if tss is None else tss - wstart
        rng = np.random.default_rng([cfg["perturb"]["seed"], zlib.crc32(args.experiment.encode()), zlib.crc32(sid.encode())])
        c0, c1, found = L.core_interval(codes, reg["motif"])
        if args.experiment == "scan":
            conds = scan_conditions(codes, reg, rng, args.tile_reps)
        else:
            conds, found = titrate_conditions(codes, reg, rng)
        seqs = [codes] + [c[5] for c in conds]
        meta = [("native", "native", 0, "", 0, None)] + [c[:5] + (c[6],) for c in conds]
        preds = ag.predict(seqs)
        rows = []
        for (name, kind, pos, alt, rep, tgt), pred in zip(meta, preds):
            for ts in ag.tracksets:
                r = {"site_id": sid, "group": row["group"], "binding_category": row["binding_category"], "stable_dynamic": row["stable_dynamic"], "native_170_pred": row["native_170"],
                     "native_core": L.decode(codes[c0:c1]), "core_found": bool(found), "condition": name, "kind": kind, "pos": pos, "alt": alt, "rep": rep, "trackset": ts,
                     "tss_dist": None if tss is None else tss - int(row["motif_mid"]), **L.readout(pred[ts], reg, tss_rel, cfg)}
                r["atac_target"] = float(pred[ts]["atac"][tgt[0]:tgt[1]].mean()) if tgt is not None and "atac" in pred[ts] else np.nan
                rows.append(r)
        L.append_rows(out, rows)
        n_done += 1
        if n_done % 5 == 0:
            print(f"{n_done} sites, {(time.time() - t0) / n_done:.1f} s/site", flush=True)
    print(f"done: {n_done} sites in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
