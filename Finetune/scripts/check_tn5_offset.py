#!/usr/bin/env python3
"""Are the fragment ends already Tn5-shifted? Decides the --offsets of make_insertion_bigwigs.py, which matter for 1 bp targets.

Each fragment gives a left cut (start) and a right cut (end - 1). When the ends are already corrected for the 9 bp Tn5 duplication (Cell Ranger and
ArchR exports are), both mark the cut centre, so their average profiles around footprinted motifs line up (lag about 0). Uncorrected ends put the left
cuts 4 bp to the left and the right cuts 5 bp to the right of the true cut, so the right-cut profile sits about 9 bp to the right of the left-cut
profile. The script piles both cut types up around the motif midpoints of a site table (+-200 bp), finds the lag that best aligns them over the central
+-100 bp, and reports the footprint (cut depletion over the motif relative to +-20 to 50 bp) for the pooled profile.
Writes <out>.json and <out>.profile.csv.

  python check_tn5_offset.py --fragments CM_YS_merged_fragments_bgz_sorted.tsv.gz --sites ../Katana/data/cm_fos_sites.csv --out checks/CM_YS_tn5
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

W, CORE, MAXLAG = 200, 100, 15


def best_lag(left, right):
    """Shift of the right-cut profile (bp, positive = right cuts lie further right) that maximises its correlation with the left-cut profile."""
    a = left[W - CORE:W + CORE + 1]
    scores = {}
    for lag in range(-MAXLAG, MAXLAG + 1):
        b = right[W - CORE + lag:W + CORE + 1 + lag]
        scores[lag] = float(np.corrcoef(a, b)[0, 1])
    lag = max(scores, key=scores.get)
    return lag, scores


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--sites", required=True, help="site table with seqnames and motif_mid (TOBIAS-footprinted AP-1 motifs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunksize", type=int, default=5_000_000)
    args = ap.parse_args()

    s = pd.read_csv(args.sites, usecols=["seqnames", "motif_mid"]).drop_duplicates()
    mids = {c: np.sort(g.motif_mid.to_numpy(np.int64)) for c, g in s.groupby("seqnames")}
    prof = {"left": np.zeros(2 * W + 1, np.int64), "right": np.zeros(2 * W + 1, np.int64)}
    reader = pd.read_csv(args.fragments, sep="\t", header=None, comment="#", usecols=[0, 1, 2], names=["chrom", "start", "end"],   # '#' only cuts the
                         dtype={"chrom": str, "start": np.int64, "end": np.int64}, chunksize=args.chunksize)                     # sample#barcode column
    n = 0
    for ch in reader:
        n += len(ch)
        for chrom, g in ch.groupby("chrom", sort=False):
            m = mids.get(chrom)
            if m is None:
                continue
            for kind, pos in (("left", g.start.to_numpy()), ("right", g.end.to_numpy() - 1)):
                i = np.clip(np.searchsorted(m, pos), 1, len(m) - 1) if len(m) > 1 else np.zeros(len(pos), int)
                near = np.where(np.abs(pos - m[i - 1]) < np.abs(pos - m[i]), m[i - 1], m[i]) if len(m) > 1 else m[i]
                rel = pos - near
                rel = rel[np.abs(rel) <= W]
                prof[kind] += np.bincount(rel + W, minlength=2 * W + 1)

    left, right = prof["left"].astype(float), prof["right"].astype(float)
    lag, scores = best_lag(left, right)
    both = left + right
    x = np.arange(-W, W + 1)
    centre = both[np.abs(x) <= 5].mean()
    flank = both[(np.abs(x) >= 20) & (np.abs(x) <= 50)].mean()
    verdict = "already Tn5-corrected: use --offsets 0 -1" if abs(lag) <= 2 else ("uncorrected: use --offsets 4 -5" if 6 <= lag <= 11 else "unclear: inspect the profile")
    summary = {"fragments_file": args.fragments, "sites": int(len(s)), "fragments_read": int(n), "cuts_within_200bp": int(both.sum()),
               "best_lag_right_vs_left_bp": int(lag), "corr_at_best_lag": scores[lag], "corr_at_lag0": scores[0], "corr_at_lag9": scores[9],
               "footprint_centre_over_flank": float(centre / flank), "verdict": verdict}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".json").write_text(json.dumps(summary, indent=1))
    pd.DataFrame({"rel_bp": x, "left_cuts": prof["left"], "right_cuts": prof["right"]}).to_csv(f"{out}.profile.csv", index=False)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
