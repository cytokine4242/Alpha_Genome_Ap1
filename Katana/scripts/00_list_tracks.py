#!/usr/bin/env python3
"""List AlphaGenome tracks from the built-in metadata catalog. CPU only, no GPU or weights needed (login node is fine).
Writes results/human_tracks.csv and results/mouse_catalog_tracks.csv and prints the rows matching an optional search.

  python 00_list_tracks.py                                  # write both tables
  python 00_list_tracks.py cardiac fibroblast myocyte       # also print rows matching any of the words (any column)
  python 00_list_tracks.py --organism human cardiac         # restrict printing to one organism
Equivalent with the package CLI:  agt info --tracks atac --organism human --search cardiac
"""
from __future__ import annotations

import sys

import pandas as pd

import ag_lib as L


def main():
    args = sys.argv[1:]
    only = None
    if "--organism" in args:
        i = args.index("--organism")
        only = {"human": 0, "mouse": 1}[args[i + 1]]
        del args[i:i + 2]
    out = L.ROOT / "results"
    out.mkdir(exist_ok=True)
    for org, name in ((0, "human_tracks.csv"), (1, "mouse_catalog_tracks.csv")):
        df = pd.DataFrame(L.catalog_tracks(org))
        df.to_csv(out / name, index=False)
        print(f"organism {org}: {len(df)} tracks -> results/{name}\n{df.groupby('head').size().to_string()}")
        if args and (only is None or only == org):
            pat = "|".join(args)
            hit = df[df.astype(str).apply(lambda c: c.str.contains(pat, case=False, regex=True)).any(axis=1)]
            cols = [c for c in ("head", "index", "track_name", "biosample_name", "biosample_type", "biosample_life_stage", "strand", "assay_title", "data_source") if c in hit]
            print(f"--- organism {org}, rows matching /{pat}/: {len(hit)}")
            print(hit[cols].to_string(index=False))


if __name__ == "__main__":
    main()
