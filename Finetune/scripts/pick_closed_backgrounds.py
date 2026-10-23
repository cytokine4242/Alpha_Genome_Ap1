#!/usr/bin/env python3
"""Pick quiet genomic windows (no measured ATAC signal in any heart cell type) as neutral hosts for element-insertion experiments.

An element inserted into a window that already has its own regulatory activity cannot be read cleanly, so the insertion
experiments need host windows with as little measured accessibility as possible, in every cell type rather than only the
two we fine-tuned on. Candidate centres are drawn at random on chr1-19 and chrX and kept when:
  - the whole model window (--window bp) lies inside the chromosome and holds no N,
  - the centre is at least --min-dist-peak bp from every interval of every --peaks file (all cell types' called peaks),
  - the measured density over the central --core bp is at most --max-rel times that track's genome-wide mean, in EVERY
    --bigwig. The threshold is relative because the 1 bp insertion tracks and the 10 bp coverage tracks have different
    scales; the genome-wide mean comes from each file's own header.
From the quiet pool, --n windows are chosen so that their GC content matches the elements that will be inserted
(--gc-from, a site table: the GC of each site's peak), by taking the quiet window closest to each of n evenly spaced
quantiles of the element GC distribution. GC therefore cannot differ systematically between the hosts and the inserts.
Pass --no-gc-match to simply take the quietest windows instead.

Writes a CSV: chrom, centre, window_start, window_end, gc_core, total_rel and <track>_rel for each bigWig
(measured density over the core as a multiple of that track's genome-wide mean).

  python pick_closed_backgrounds.py --genome mm10.fa --out FT/closed_backgrounds.csv --n 20 \\
      --bigwigs FT/bigwigs_1bp/*.insertions.bw RAW/*_YS_merged.SeqDepthNorm.bw --peaks RAW/*.narrowPeak \\
      --gc-from ../Katana/data/fb_swap_sites.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CANONICAL = [f"chr{i}" for i in list(range(1, 20)) + ["X"]]


def load_peaks(files):
    """{chrom: (starts, ends)} sorted, merged over every peak file, so a candidate can be tested against all cell types at once."""
    frames = []
    for f in files:
        p = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                        dtype={"chrom": str, "start": np.int64, "end": np.int64}, comment="#")
        frames.append(p)
    allp = pd.concat(frames, ignore_index=True)
    out = {}
    for c, g in allp.sort_values("start").groupby("chrom"):
        s, e = g.start.to_numpy(), g.end.to_numpy()
        keep_s, keep_e = [s[0]], [e[0]]                      # merge overlapping intervals so the gap test is exact
        for a, b in zip(s[1:], e[1:]):
            if a <= keep_e[-1]:
                keep_e[-1] = max(keep_e[-1], b)
            else:
                keep_s.append(a)
                keep_e.append(b)
        out[c] = (np.array(keep_s), np.array(keep_e))
    return out


def distance_to_peak(peaks, chrom, pos):
    """Distance from pos to the nearest peak interval on chrom (0 inside a peak); inf when the chromosome has none."""
    pk = peaks.get(chrom)
    if pk is None or not len(pk[0]):
        return np.inf
    s, e = pk
    i = np.searchsorted(s, pos, side="right") - 1
    d = np.inf
    if i >= 0:
        d = min(d, 0.0 if pos < e[i] else pos - e[i])        # the interval starting at or before pos
    if i + 1 < len(s):
        d = min(d, s[i + 1] - pos)                            # the next interval
    return d


def gc_of_sites(path, fa, max_sites=400):
    """GC fraction of each site's peak sequence: the distribution the host windows should match."""
    s = pd.read_csv(path).drop_duplicates("site_id")
    if len(s) > max_sites:
        s = s.sample(n=max_sites, random_state=0)
    gc = []
    for r in s.itertuples():
        seq = fa[r.peak_chr][int(r.peak_start):int(r.peak_end)]
        if seq and "N" not in seq:
            gc.append((seq.count("G") + seq.count("C")) / len(seq))
    return np.array(gc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bigwigs", nargs="+", required=True, help="every track the window must be quiet in")
    ap.add_argument("--peaks", nargs="+", required=True, help="narrowPeak files of all cell types")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gc-from", default=None, help="site table whose peak GC the chosen windows should match")
    ap.add_argument("--no-gc-match", action="store_true", help="take the quietest windows instead of matching GC")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--window", type=int, default=131072, help="model window that must be inside the chromosome and N-free")
    ap.add_argument("--core", type=int, default=4096, help="central region whose measured density must be low")
    ap.add_argument("--max-rel", type=float, default=0.1, help="core density allowed, as a multiple of each track's genome-wide mean")
    ap.add_argument("--min-dist-peak", type=int, default=25_000, help="minimum distance from the centre to any called peak")
    ap.add_argument("--tries", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=46)
    args = ap.parse_args()

    import pyBigWig
    from pyfaidx import Fasta

    fa = Fasta(args.genome, as_raw=True, sequence_always_upper=True)
    peaks = load_peaks(args.peaks)
    print(f"{len(args.peaks)} peak files -> {sum(len(v[0]) for v in peaks.values())} merged intervals", flush=True)

    bws, means = {}, {}
    for f in args.bigwigs:
        name = Path(f).name.split(".")[0].replace("_merged", "")
        bw = pyBigWig.open(f)
        h = bw.header()
        means[name] = h["sumData"] / max(h["nBasesCovered"], 1)
        bws[name] = bw
    print("tracks and their genome-wide means:", {k: round(v, 5) for k, v in means.items()}, flush=True)

    sizes = {c: len(fa[c]) for c in CANONICAL if c in fa}
    chroms = list(sizes)
    w = np.array([sizes[c] for c in chroms], float)
    w /= w.sum()
    rng = np.random.default_rng(args.seed)
    half, chalf = args.window // 2, args.core // 2

    kept = []
    for _ in range(args.tries):
        c = chroms[rng.choice(len(chroms), p=w)]
        centre = int(rng.integers(half, sizes[c] - half))
        if distance_to_peak(peaks, c, centre) < args.min_dist_peak:
            continue
        seq = fa[c][centre - chalf:centre + chalf]
        if "N" in seq:
            continue
        rel = {}
        for name, bw in bws.items():
            if c not in bw.chroms():
                rel = None
                break
            v = bw.stats(c, centre - chalf, centre + chalf, type="mean")[0]
            rel[name] = 0.0 if v is None else float(v) / means[name]
            if rel[name] > args.max_rel:
                rel = None
                break
        if rel is None:
            continue
        kept.append({"chrom": c, "centre": centre, "window_start": centre - half, "window_end": centre + half,
                     "gc_core": (seq.count("G") + seq.count("C")) / len(seq),
                     "total_rel": sum(rel.values()), **{f"{k}_rel": v for k, v in rel.items()}})
    pool = pd.DataFrame(kept)
    print(f"{len(pool)} of {args.tries} candidate centres are quiet in all {len(bws)} tracks", flush=True)
    if pool.empty:
        raise SystemExit("no quiet windows found: raise --max-rel or --tries, or lower --min-dist-peak")

    pool = pool.sort_values("total_rel").reset_index(drop=True)
    pool["n_free"] = [int("N" not in fa[r.chrom][r.window_start:r.window_end]) for r in pool.itertuples()]
    pool = pool[pool.n_free == 1].drop(columns="n_free").reset_index(drop=True)

    if args.no_gc_match or not args.gc_from:
        out = pool.head(args.n)
        print("GC matching off: taking the quietest windows")
    else:
        target = gc_of_sites(args.gc_from, fa)
        want = np.quantile(target, np.linspace(0.05, 0.95, args.n))
        avail, chosen = pool.copy(), []
        for t in want:
            j = (avail.gc_core - t).abs().idxmin()
            chosen.append(avail.loc[j])
            avail = avail.drop(index=j)
        out = pd.DataFrame(chosen).reset_index(drop=True)
        print(f"GC matched to {len(target)} element peaks: elements median {np.median(target):.3f} "
              f"(IQR {np.percentile(target, 25):.3f}-{np.percentile(target, 75):.3f}), "
              f"chosen windows median {out.gc_core.median():.3f} (IQR {out.gc_core.quantile(.25):.3f}-{out.gc_core.quantile(.75):.3f})")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"wrote {len(out)} windows to {args.out}")
    print(out[["chrom", "centre", "gc_core", "total_rel"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
