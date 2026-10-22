#!/bin/bash
# Mouse train/valid/test regions in AlphaGenome's 1 Mb format, from the Borzoi mouse sequence folds that are already on Katana (rubj/beagle/rawData/sequences_mouse.bed,
# 8 folds of ~196 kb windows). Valid and test windows whose 1 Mb window overlaps a training window are removed (leakage control). Also writes train_pilot.bed / valid_pilot.bed.
#   bash make_mouse_folds.sh OUT_DIR
set -e
OUT=${1:?usage: make_mouse_folds.sh OUT_DIR}
SEQ=${SEQ_BED:-/srv/scratch/wonglab/rubj/beagle/rawData/sequences_mouse.bed}
PY=${FT_PY:-/srv/scratch/wonglab/minh/envs/alphagenome-pytorch_venv/bin/python}
CONV=${FT_CONVERT:-/srv/scratch/wonglab/minh/alphagenome/alphagenome-pytorch/scripts/convert_borzoi_folds.py}
mkdir -p "$OUT"
$PY "$CONV" --input "$SEQ" --output-dir "$OUT" --organism mouse
for f in "$OUT"/FOLD_*; do
  shuf -n 200 --random-source=<(yes) "$f/train.bed" > "$f/train_pilot.bed"      # pilot: 200 training windows
  shuf -n 40 --random-source=<(yes) "$f/valid.bed" > "$f/valid_pilot.bed"       # and 40 validation windows, so the pilot is quick
done
echo "folds written to $OUT (FOLD_3 pairs with model_fold_3.safetensors: valid = Borzoi fold6, test = fold7)"
