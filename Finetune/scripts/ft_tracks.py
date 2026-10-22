"""The tracks used for fine-tuning, in one place so the YAML configs, the track metadata and the evaluation all use the same order."""
CELLS = {"Cardiomyocyte": "CL:0000746", "Fibroblast": "CL:0000057"}      # cell ontology ids (cardiac muscle cell, fibroblast)
CONDITIONS = ["YS", "MAS", "MAE", "YE"]                                      # young sedentary, aged sedentary, aged + exercise, young + exercise
RNA = {"Cardiomyocyte": "Cardiomyocyte_RNA.bw", "Fibroblast": "Fibroblast_RNA.bw"}   # merged over conditions (the only RNA coverage available)


def atac_tracks(conditions=None):
    """[(cell, condition, track_name)] in cell-major order; this order is the output-track order of the head."""
    return [(c, k, f"{c}_{k}") for c in CELLS for k in (conditions or CONDITIONS)]


def atac_file(kind, cell, cond):
    """kind '10bp': existing SeqDepthNorm bigWigs (bamCoverage 10 bp bins, CPM); kind '1bp': Tn5 insertion bigWigs made by make_insertion_bigwigs.py."""
    return f"{cell}_{cond}_merged.SeqDepthNorm.bw" if kind == "10bp" else f"{cell}_{cond}.insertions.bw"
