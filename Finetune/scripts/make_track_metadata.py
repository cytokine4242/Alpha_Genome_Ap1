#!/usr/bin/env python3
"""Writes the track metadata TSV for the fine-tuned heads (same track order as the configs). Embedded in the checkpoints so predictions carry track names.
  python make_track_metadata.py --out ../data/track_metadata.tsv [--with-rna] [--conditions YS MAS MAE YE]"""
import argparse

import pandas as pd
from ft_tracks import CELLS, CONDITIONS, RNA, atac_tracks

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--with-rna", action="store_true")
ap.add_argument("--conditions", nargs="+", default=CONDITIONS)
a = ap.parse_args()
rows = [{"track_name": n, "output_type": "atac", "organism": "mouse", "biosample_name": c, "ontology_curie": CELLS[c], "assay_title": "snATAC-seq", "condition": k, "biosample_type": "primary cell"}
        for c, k, n in atac_tracks(a.conditions)]
if a.with_rna:
    rows += [{"track_name": f"{c}_RNA", "output_type": "rna_seq", "organism": "mouse", "biosample_name": c, "ontology_curie": CELLS[c], "assay_title": "snRNA-seq", "condition": "all", "biosample_type": "primary cell"} for c in RNA]
pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
print(f"wrote {len(rows)} tracks to {a.out}")
