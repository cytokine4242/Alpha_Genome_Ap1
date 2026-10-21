#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh pilot            # 10 CM sites, 10 FB sites, tiny insertion, 20 MPRA elements
#   bash submit_all.sh full             # everything: 40 chunks CM, 40 chunks FB, insertion, 4 MPRA chunks
#
# Options (environment variables):
#   GPU=H200 (default) | A100 | L40S | GH200 | ... | any     GPU model; "any" drops the gpu_model request so the
#                                                           scheduler uses whichever GPU is free first
#   BATCH=N       sequences per GPU batch (config.json model.batch_size). Default: config value on H200/GH200,
#                 4 on any other model or "any" (smaller cards have less memory)
#   DRY=1         print the qsub commands instead of submitting
# Examples:
#   GPU=any bash submit_all.sh pilot
#   GPU=A100 BATCH=8 bash submit_all.sh full
#   DRY=1 GPU=any bash submit_all.sh full
# Logs go to logs/ (the GPU used is printed at the top of each log), results to results/.
set -e
cd "$(dirname "$0")"
mkdir -p logs results

GPU=${GPU:-H200}
if [ "$GPU" = "any" ]; then GPU_SEL=""; else GPU_SEL=":gpu_model=$GPU"; fi
RES="-l select=1:ncpus=6:mem=96gb:ngpus=1${GPU_SEL}"   # command-line resources override the #PBS lines in pbs/*.pbs
if [ -z "$BATCH" ] && [ "$GPU" != "H200" ] && [ "$GPU" != "GH200" ]; then BATCH=4; fi
BATCH_VAR=${BATCH:+,AG_BATCH_SIZE=$BATCH}
O="-o $PWD/logs/"

# sub <extra qsub args...> -- <comma-separated -v vars (may be empty)> <pbs file>
sub() {
  local vars="${@: -2:1}" file="${@: -1}"
  local args=("${@:1:$#-2}")
  local v="${vars}${BATCH_VAR}"; v="${v#,}"
  local cmd=(qsub $O $RES "${args[@]}" ${v:+-v "$v"} "$file")
  if [ -n "$DRY" ]; then echo "${cmd[@]}"; else "${cmd[@]}"; fi
}

echo "GPU=$GPU${BATCH:+  batch=$BATCH}  resources: $RES"
case "$1" in
  check) sub "" pbs/00_check_env.pbs ;;
  pilot)
    W="-l walltime=01:00:00"  # short walltime starts sooner
    sub $W -N ag_pilot_cm  "SITES=cm,NCHUNKS=1,LIMIT=10" pbs/11_perturb.pbs
    sub $W -N ag_pilot_fb  "SITES=fb,NCHUNKS=1,LIMIT=10" pbs/11_perturb.pbs
    sub $W -N ag_pilot_ins "LIMIT=3" pbs/12_insertion.pbs
    sub $W -N ag_pilot_mpra "NCHUNKS=1,LIMIT=20" pbs/13_mpra.pbs ;;
  full)
    sub -J 0-39 "SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs
    sub -J 0-39 "SITES=fb,NCHUNKS=40" pbs/11_perturb.pbs
    sub "" pbs/12_insertion.pbs
    sub -J 0-3 "NCHUNKS=4" pbs/13_mpra.pbs ;;
  *) echo "usage: [GPU=H200|A100|...|any] [BATCH=N] [DRY=1] bash submit_all.sh check|pilot|full"; exit 1 ;;
esac
