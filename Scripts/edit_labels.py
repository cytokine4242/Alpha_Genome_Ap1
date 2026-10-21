"""Plain-language labels for the edit conditions, with the exact number of bases changed. Used by every figure and table so the wording is the same everywhere.
Condition names in the result files are <kind>@<region> (for example random@motif); the report shows the labels below instead."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_w = pd.read_csv(ROOT / "Data" / "derived" / "cm_fos_sites.csv").peak_width
PEAK_MEDIAN = int(round(_w.median()))                      # median width of the cardiomyocyte TOBIAS peaks (bp)
PEAK_Q1, PEAK_Q3 = int(round(_w.quantile(0.25))), int(round(_w.quantile(0.75)))
PEAK_TXT = f"peak, median {PEAK_MEDIAN} bp"

# bases changed: motif 11 bp (core 7 bp = TGA[C/G]TCA, so 4 flanking motif bases); 170 bp window mid +- 85; flank@170 keeps the 11 bp motif (80 + 79 = 159 bp changed)
EDIT_LABEL = {
    "shuffle@motif": "Shuffle motif (11 bp)",
    "shuffle@170": "Shuffle 170 bp",
    "shuffle@peak": f"Shuffle whole peak ({PEAK_TXT})",
    "random@motif": "Replace motif (11 bp)",
    "random@170": "Replace 170 bp",
    "random@peak": f"Replace whole peak ({PEAK_TXT})",
    "dinuc@motif": "Dinucleotide-shuffle motif (11 bp)",
    "dinuc@170": "Dinucleotide-shuffle 170 bp",
    "dinuc@peak": f"Dinucleotide-shuffle whole peak ({PEAK_TXT})",
    "core@motif": "Replace 7 bp core only",
    "flank@motif": "Replace the 4 motif flank bases (core kept)",
    "flank@170": "Replace the 159 bp around the motif (motif kept)",
}
BASES_CHANGED = {
    "shuffle@motif": "11 (letters permuted)", "shuffle@170": "170 (permuted)", "shuffle@peak": f"whole peak, median {PEAK_MEDIAN} (permuted)",
    "random@motif": "11", "random@170": "170", "random@peak": f"whole peak, median {PEAK_MEDIAN}",
    "dinuc@motif": "11 (dinucleotides kept)", "dinuc@170": "170 (dinucleotides kept)", "dinuc@peak": f"whole peak, median {PEAK_MEDIAN} (dinucleotides kept)",
    "core@motif": "7", "flank@motif": "4", "flank@170": "159 (80 left + 79 right)",
}
REGION_LABEL = {"motif": "Within the 11 bp motif", "170": "Within the 170 bp", "peak": "Whole peak"}


def short(cond, width=16):
    """Short, wrapped label for a plot axis."""
    import textwrap
    return "\n".join(textwrap.wrap(EDIT_LABEL[cond].replace(f" ({PEAK_TXT})", f"\n(median {PEAK_MEDIAN} bp)"), width))
