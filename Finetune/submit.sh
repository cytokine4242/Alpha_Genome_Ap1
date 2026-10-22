#!/bin/bash
# Submit fine-tuning jobs from the Finetune/ folder on Katana (private account by default; DRY=1 prints the qsub commands instead).
#   bash submit.sh check_fragments          # array: checks the 8 merged fragment files copied from Gadi (finetune_data/fragments) before they become 1 bp targets
#   bash submit.sh bigwigs                  # array: one task per (cell type, condition) of data/samples_manifest.tsv -> 1 bp insertion bigWigs  (needs MANIFEST, CHROM_SIZES, OUT_DIR)
#   bash submit.sh probe|lora|lora_locon|full   # config chosen by KIND (10bp or 1bp, default 10bp), PILOT=1, RNA=1, SEQ (131072 default, or 1048576), FOLD (FOLD_3)
# Options (environment variables): ACCOUNT (default mwac_wonglab), GPU (default any; e.g. H200), MEM (host RAM per GPU, default 62gb: above 45gb excludes every V100 node, above 62gb also
#   excludes the 8-GPU general A100 nodes, above 124gb leaves only H200-class nodes; the probe pilot used 19 GB), BATCH / ACCUM (override the YAML batch size / accumulation, e.g. BATCH=1 ACCUM=32 on a 48 GB card), WALL (walltime)
set -e
cd "$(dirname "$0")"
ACCOUNT=${ACCOUNT-mwac_wonglab}
GPU=${GPU:-any}
mkdir -p logs
if [ "$GPU" = "any" ]; then SEL=""; else SEL=":gpu_model=$GPU"; fi
O="-o $PWD/logs/ ${ACCOUNT:+-A $ACCOUNT}"
OC="-o $PWD/logs/"   # CPU jobs: mwac_wonglab has GPU queues only, so these go to the general CPU queues without an account
run() { if [ -n "$DRY" ]; then echo qsub "$@"; else qsub "$@"; fi; }
case "$1" in
  bigwigs)
    : "${MANIFEST:?set MANIFEST}" "${CHROM_SIZES:?set CHROM_SIZES}" "${OUT_DIR:?set OUT_DIR}"
    N=$(python3 scripts/make_bigwig_from_manifest.py --manifest "$MANIFEST" --chrom-sizes "$CHROM_SIZES" --out-dir "$OUT_DIR" --list | wc -l)
    [ "$N" -gt 0 ] || { echo "$MANIFEST has no sample rows: add one row per library first"; exit 1; }
    if [ "$N" -gt 1 ]; then J="-J 0-$((N-1))"; else J=""; fi   # a one-track manifest runs as a plain job (index 0)
    run $OC $J -N ft_bigwig -v "MANIFEST=$MANIFEST,CHROM_SIZES=$CHROM_SIZES,OUT_DIR=$OUT_DIR" pbs/make_bigwig.pbs ;;
  check_fragments)   # one CPU task per merged fragment file: cell type, condition, purity, agreement with the 10 bp bigWig -> $FT/fragments/checks/
    run $OC -J 0-7 -N ft_check_frag pbs/check_fragments.pbs ;;
  probe|lora|lora_locon|full)
    KIND=${KIND:-10bp}
    SEQ=${SEQ:-131072}; FOLD=${FOLD:-FOLD_3}
    CFG=configs/${1}_${KIND}_$(echo $FOLD | tr A-Z a-z)_$((SEQ / 1024))kb${PILOT:+_pilot}${RNA:+_rna}.yaml
    [ -f "$CFG" ] || { echo "$CFG not found: run scripts/make_configs.py first"; exit 1; }
    if [ -n "$PILOT" ]; then W=${WALL:-02:00:00}; else case "$1" in probe) W=${WALL:-24:00:00} ;; *) W=${WALL:-72:00:00} ;; esac; fi
    # adapters and full fine-tuning backpropagate through the trunk: on cards smaller than an H200 use batch 1 with 32-step accumulation unless BATCH/ACCUM are given
    if [ "$1" != probe ] && [ "$GPU" != H200 ] && [ -z "$BATCH" ]; then BATCH=1; ACCUM=${ACCUM:-32}; fi
    V="CONFIG=$CFG${BATCH:+,BATCH=$BATCH}${ACCUM:+,ACCUM=$ACCUM}"
    run $O -l select=1:ncpus=8:mem=${MEM:-62gb}:ngpus=1$SEL -l walltime=$W -N ft_${1}_${KIND}${PILOT:+_pilot} -v "$V" pbs/finetune.pbs ;;
  *) echo "usage: bash submit.sh check_fragments|bigwigs|probe|lora|lora_locon|full"; exit 1 ;;
esac
