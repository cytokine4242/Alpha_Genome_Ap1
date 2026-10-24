#!/usr/bin/env python3
"""Is AP-1 special, or would any footprinted motif do? Shuffles the site's OWN AP-1 motif and, in the same window, each of
its footprinted partner motifs, so the comparison is within a site and no accessibility matching is needed.

Every claim so far rests on shuffling AP-1 and seeing a large effect; nothing else has been shuffled. This also tests
tonight's main result directly: AP-1-loss sites depend about 40% less on their own AP-1 motif than Stable sites do at
matched accessibility. If those sites also depend less on their other motifs, they are simply weak elements; if only
AP-1 is weakened, the disengagement is AP-1-specific.

Per site (data/fb_subset_scan.csv, the sites the tile scan already covers, so results line up with Analysis/fb_scan):
  own_ap1          the site's own 11 bp AP-1 motif, shuffled                                        5 draws
  partner@<family> each selected partner footprint, shuffled over its own span                      3 draws each
  ctrl@<family>    a footprint-free window of the SAME width at the mirrored distance (or the
                   nearest footprint-free position at a similar distance), shuffled                 3 draws each
The control is what makes the comparison fair: it asks whether shuffling that many bases at that distance would matter
anyway, which the earlier cardiomyocyte work only ever tested for AP-1.

Partners are chosen from Analysis/partners_fb/tables/footprints_near_sites.csv: one footprint per family, families in
the order they scored in the fibroblast partner analysis, nearest footprint first, never overlapping the site's own
motif, at most --max-partners per site.
Shuffles preserve base composition and are redrawn to leave no AP-1 half-site at the edited span or its junctions.

  python 19_partner_shuffle.py --sites ../data/fb_subset_scan.csv --out ../results_ft/pshuf_fb --chunk 0 --n-chunks 40
  python 19_partner_shuffle.py --selftest      # no GPU, no genome
"""
from __future__ import annotations

import argparse
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import ag_lib as L

# families in the order they mattered in the fibroblast partner analysis (Analysis/partners_fb)
FAMILY_ORDER = ["CEBP", "NFI (Nfia/b/c/x)", "E-box bHLH (Myod/Tcf/Neuro)", "AP-1 (Fos/Jun/Atf/Batf)",
                "CTCF / cohesin", "Heat shock / stress", "RFX", "TEAD", "NFAT", "ETS", "KLF / SP", "FOX"]
MIN_W, MAX_W = 6, 40          # ignore footprints shorter or longer than a plausible motif


def short(fam):
    return fam.split(" (")[0].replace(" / ", "_").replace(" ", "_")


def shuffled(codes, a, b, rng):
    """Composition-preserving shuffle of codes[a:b] leaving no AP-1 half-site in the span or across its junctions."""
    seg = codes[a:b].copy()
    left, right = codes[max(a - 3, 0):a], codes[b:b + 3]
    for _ in range(100):
        s = rng.permutation(seg)
        if (s != seg).any() and not any(h in L.decode(np.concatenate([left, s, right])) for h in L.HALF_SITES):
            return s
    return rng.permutation(seg)


def control_span(width, want_dist, occupied, mid, size, rng):
    """A footprint-free span of `width` at about `want_dist` bp from the motif, preferring the mirrored side."""
    for dist in [-want_dist] + [d for off in (0, 20, 40, 60, 80, 120) for d in (want_dist + off, -want_dist - off)]:
        a = mid + dist - width // 2
        b = a + width
        if a < 0 or b > size:
            continue
        if not any(a < e and s < b for s, e in occupied):
            return a, b
    return None


def pick_partners(fp, motif_rel, max_partners):
    """One footprint per family, best families first, nearest first, never touching the site's own motif."""
    out, seen = [], set()
    fp = fp.assign(width=fp.rel_end - fp.rel_start)
    fp = fp[(fp.width >= MIN_W) & (fp.width <= MAX_W)]
    fp = fp[(fp.rel_end <= motif_rel[0]) | (fp.rel_start >= motif_rel[1])]      # no overlap with our own motif
    fp = fp.assign(dist=fp[["rel_start", "rel_end"]].abs().min(axis=1))
    rank = {f: i for i, f in enumerate(FAMILY_ORDER)}
    fp = fp.assign(frank=fp.family.map(rank).fillna(len(FAMILY_ORDER))).sort_values(["frank", "dist"])
    for r in fp.itertuples():
        if r.family in seen:
            continue
        seen.add(r.family)
        out.append((r.family, int(r.rel_start), int(r.rel_end)))
        if len(out) >= max_partners:
            break
    return out


def build(codes, reg, partners, rng, reps_own=5, reps=3):
    """[(condition, kind, family, width, rep, edited window, (a, b) edited span)]"""
    mid = reg["mid"]
    m0, m1 = reg["motif"]
    out = []
    for r in range(reps_own):
        e = codes.copy()
        e[m0:m1] = shuffled(codes, m0, m1, rng)
        out.append(("own_ap1", "own_ap1", "AP-1 (own)", m1 - m0, r, e, (m0, m1)))
    occupied = [(mid + s, mid + e_) for _, s, e_ in partners] + [(m0, m1)]
    for fam, rs, re_ in partners:
        a, b = mid + rs, mid + re_
        if a < 0 or b > len(codes) or L.has_n(codes, a, b):
            continue
        for r in range(reps):
            e = codes.copy()
            e[a:b] = shuffled(codes, a, b, rng)
            out.append((f"partner@{short(fam)}", "partner", fam, b - a, r, e, (a, b)))
        cs = control_span(b - a, int(min(abs(rs), abs(re_))), occupied, mid, len(codes), rng)
        if cs and not L.has_n(codes, *cs):
            for r in range(reps):
                e = codes.copy()
                e[cs[0]:cs[1]] = shuffled(codes, cs[0], cs[1], rng)
                out.append((f"ctrl@{short(fam)}", "ctrl", fam, cs[1] - cs[0], r, e, cs))
    return out


def selftest():
    rng = np.random.default_rng(19)
    size = 131072
    codes = rng.integers(0, 4, size=size).astype(np.uint8)
    mid = size // 2
    m0, m1 = mid - 5, mid + 6
    codes[m0 + 2:m0 + 9] = L.encode("TGAGTCA")
    reg = {"motif": (m0, m1), "mid": mid, "w170": (mid - 85, mid + 85), "peak": (mid - 400, mid + 400)}
    fp = pd.DataFrame({"family": ["CEBP", "TEAD", "CEBP", "NFI (Nfia/b/c/x)", "X"],
                       "rel_start": [-60, 40, -200, 100, -3], "rel_end": [-48, 52, -188, 112, 8]})
    pk = pick_partners(fp, (m0 - mid, m1 - mid), 6)
    fams = [f for f, _, _ in pk]
    assert "CEBP" in fams and fams.count("CEBP") == 1, fams          # one footprint per family
    assert ("X", -3, 8) not in pk, "a footprint overlapping our own motif must be skipped"
    assert fams[0] == "CEBP" and [f for f, _, _ in pk].index("NFI (Nfia/b/c/x)") < fams.index("TEAD"), fams
    assert pk[0][1] == -60, "nearest CEBP first"

    conds = build(codes, reg, pk, rng)
    kinds = [c[1] for c in conds]
    assert kinds.count("own_ap1") == 5
    for name, kind, fam, width, rep, e, span in conds:
        d = np.flatnonzero(e != codes)
        assert len(d) and set(d) <= set(range(*span)), name
        assert sorted(e[span[0]:span[1]]) == sorted(codes[span[0]:span[1]]), f"{name}: composition must be preserved"
        if kind == "own_ap1":
            assert span == (m0, m1)
        if kind == "ctrl":
            assert not (span[0] < m1 and m0 < span[1]), f"{name}: a control must not touch the site's own motif"
    for fam, rs, re_ in pk:
        w = re_ - rs
        pc = [c for c in conds if c[0] == f"ctrl@{short(fam)}"]
        pp = [c for c in conds if c[0] == f"partner@{short(fam)}"]
        if pc:
            assert pc[0][3] == pp[0][3] == w, f"{fam}: control width must equal the partner width"
    print(f"selftest OK: {len(pk)} partners chosen, {len(conds)} edits "
          f"({kinds.count('own_ap1')} own AP-1, {kinds.count('partner')} partner, {kinds.count('ctrl')} control)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default="../data/fb_subset_scan.csv")
    ap.add_argument("--footprints", default="../../Analysis/partners_fb/tables/footprints_near_sites.csv")
    ap.add_argument("--out")
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--n-chunks", type=int, default=1)
    ap.add_argument("--max-partners", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--config", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    assert args.out, "--out is required"

    cfg = L.load_config(args.config)
    size = cfg["model"]["window"]
    genome = L.Genome(cfg["paths"]["genome"])
    ag = L.AG(cfg)
    sites = pd.read_csv(args.sites).iloc[args.chunk::args.n_chunks]
    if args.limit:
        sites = sites.head(args.limit)
    FP = pd.read_csv(args.footprints)
    FP = FP[~FP.own_ap1] if "own_ap1" in FP else FP
    by_site = {s: g for s, g in FP.groupby("site_id")}
    out = Path(f"{args.out}_chunk{args.chunk}.csv")
    done = L.done_ids(out)
    t0, n = time.time(), 0
    for row in sites.itertuples():
        if row.site_id in done:
            continue
        codes, ws = genome.window(row.seqnames, row.motif_mid, size)
        if codes is None:
            continue
        reg = L.site_regions(row._asdict(), ws, size, cfg)
        fp = by_site.get(row.site_id)
        if fp is None or not len(fp):
            continue
        rng = np.random.default_rng([cfg["perturb"]["seed"], zlib.crc32(row.site_id.encode())])
        partners = pick_partners(fp, (reg["motif"][0] - reg["mid"], reg["motif"][1] - reg["mid"]), args.max_partners)
        conds = build(codes, reg, partners, rng)
        if not conds:
            continue
        preds = ag.predict([codes] + [c[5] for c in conds])
        meta = [("native", "native", "", 0, 0, None)] + [c[:5] + (c[6],) for c in conds]
        rows = []
        for (name, kind, fam, width, rep, span), pred in zip(meta, preds):
            for ts in ag.tracksets:
                rows.append({"site_id": row.site_id, "group": row.group, "stable_dynamic": row.stable_dynamic,
                             "age_class": getattr(row, "age_class", ""), "ex_class": getattr(row, "ex_class", ""),
                             "condition": name, "kind": kind, "family": fam, "edit_width": width, "rep": rep,
                             "dist_to_motif": "" if span is None else int(min(abs(span[0] - reg["mid"]), abs(span[1] - reg["mid"]))),
                             "trackset": ts, **L.readout(pred[ts], reg, None, cfg)})
        L.append_rows(out, rows)
        n += 1
        if n % 10 == 0:
            print(f"{n} sites, {(time.time() - t0) / n:.1f} s/site", flush=True)
    print(f"done: {n} sites in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
