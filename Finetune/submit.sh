#!/bin/bash
# Submit fine-tuning jobs from the Finetune/ folder on Katana (private account by default; DRY=1 prints the qsub commands instead).
#   bash submit.sh bigwigs                  # array: one task per (cell type, condition) of data/samples_manifest.tsv -> 1 bp insertion bigWigs  (needs MANIFEST, CHROM_SIZES, OUT_DIR)
#   bash submit.sh probe|lora|lora_locon|full   # CONFIG=configs/<mode>_<kind>[_pilot][_rna].yaml is chosen by KIND (10bp or 1bp, default 10bp), PILOT=1, RNA=1
# Options (environment variables): ACCOUNT (default mwac_wonglab), GPU (default any; e.g. H200), BATCH / ACCUM (override the YAML batch size / accumulation, e.g. BATCH=1 ACCUM=32 on a 48 GB card), WALL (walltime)
set -e
cd "$(dirname "$0")"
ACCOUNT=${ACCOUNT-mwac_wonglab}
GPU=${GPU:-any}
mkdir -p logs
if [ "$GPU" = "any" ]; then SEL=""; else SEL=":gpu_model=$GPU"; fi
O="-o $PWD/logs/ ${ACCOUNT:+-A $ACCOUNT}"
run() { if [ -n "$DRY" ]; then echo qsub "$@"; else qsub "$@"; fi; }
case "$1" in
  bigwigs)
    : "${MANIFEST:?set MANIFEST}" "${CHROM_SIZES:?set CHROM_SIZES}" "${OUT_DIR:?set OUT_DIR}"
    N=$(python3 scripts/make_bigwig_from_manifest.py --manifest "$MANIFEST" --chrom-sizes "$CHROM_SIZES" --out-dir "$OUT_DIR" --list | wc -l)
    run $O -J 0-$((N-1)) -N ft_bigwig -v "MANIFEST=$MANIFEST,CHROM_SIZES=$CHROM_SIZES,OUT_DIR=$OUT_DIR" pbs/make_bigwig.pbs ;;
  probe|lora|lora_locon|full)
    KIND=${KIND:-10bp}
    CFG=configs/${1}_${KIND}${PILOT:+_pilot}${RNA:+_rna}.yaml
    [ -f "$CFG" ] || { echo "$CFG not found: run scripts/make_configs.py first"; exit 1; }
    case "$1" in probe) W=${WALL:-12:00:00} ;; full) W=${WALL:-48:00:00} ;; *) W=${WALL:-24:00:00} ;; esac
    V="CONFIG=$CFG${BATCH:+,BATCH=$BATCH}${ACCUM:+,ACCUM=$ACCUM}"
    run $O -l select=1:ncpus=8:mem=128gb:ngpus=1$SEL -l walltime=$W -N ft_${1}_${KIND}${PILOT:+_pilot} -v "$V" pbs/finetune.pbs ;;
  *) echo "usage: bash submit.sh bigwigs|probe|lora|lora_locon|full"; exit 1 ;;
esac
