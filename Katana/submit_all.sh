#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh pilot            # 30 CM sites, 30 FB sites, tiny insertion, 150 MPRA elements, with CAGE (pilot_* outputs)
#   bash submit_all.sh full             # everything: 40 chunks CM, 40 chunks FB, insertion, 4 MPRA chunks
#
# Options (environment variables):
#   GPU=H200 (default) | A100 | L40S | GH200 | ... | any     GPU model; "any" drops the gpu_model request so the
#                                                           scheduler uses whichever GPU is free first
#   ACCOUNT=mwac_wonglab   charge/run under a PBS account (-A), e.g. to use that group's private GPUs.
#                 When ACCOUNT is set and GPU is not, GPU defaults to "any" (the private nodes' GPU model is not
#                 requested explicitly). Single jobs (check, pilot, insertion) also get -m ae (email on abort/end);
#                 array jobs do not, to avoid one email per task.
#   CONFIG=file   config to use (default: config_pilot.json for check/pilot, which adds CAGE; config.json for full)
#   BATCH=N       sequences per GPU batch (config.json model.batch_size). Default: config value on H200/GH200,
#                 4 on any other model or "any" (smaller cards have less memory)
#   DRY=1         print the qsub commands instead of submitting
# Examples:
#   ACCOUNT=mwac_wonglab bash submit_all.sh check     # does my account have access? see logs/ for the queue + GPU
#   GPU=any bash submit_all.sh pilot
#   GPU=A100 BATCH=8 bash submit_all.sh full
#   DRY=1 GPU=any bash submit_all.sh full
# Logs: logs/<date>_<jobname>[_chunkNNofMM]_<GPU>_<jobid>.log (readable); PBS raw output in logs/pbs/. Results: results/.
#   bash submit_all.sh logs      lists the newest logs
set -e
cd "$(dirname "$0")"
mkdir -p logs/pbs results

if [ -n "$ACCOUNT" ]; then GPU=${GPU:-any}; else GPU=${GPU:-H200}; fi
if [ "$GPU" = "any" ]; then GPU_SEL=""; else GPU_SEL=":gpu_model=$GPU"; fi
RES="-l select=1:ncpus=6:mem=96gb:ngpus=1${GPU_SEL}"   # command-line resources override the #PBS lines in pbs/*.pbs
if [ -z "$BATCH" ] && [ "$GPU" != "H200" ] && [ "$GPU" != "GH200" ]; then BATCH=4; fi
BATCH_VAR=${BATCH:+,AG_BATCH_SIZE=$BATCH}
O="-o $PWD/logs/pbs/${ACCOUNT:+ -A $ACCOUNT}"
case "$1" in full) CFG=${CONFIG:-config.json} ;; *) CFG=${CONFIG:-config_pilot.json} ;; esac
[ -f "$CFG" ] || { echo "config $CFG not found"; exit 1; }
BATCH_VAR="${BATCH_VAR},AG_CONFIG=$PWD/$CFG"

# sub <extra qsub args...> -- <comma-separated -v vars (may be empty)> <pbs file>
sub() {
  local vars="${@: -2:1}" file="${@: -1}"
  local args=("${@:1:$#-2}")
  local v="${vars}${BATCH_VAR}"; v="${v#,}"
  local cmd=(qsub $O $RES "${args[@]}" ${v:+-v "$v"} "$file")
  if [ -n "$DRY" ]; then echo "${cmd[@]}"; else "${cmd[@]}"; fi
}

echo "config=$CFG  GPU=$GPU${ACCOUNT:+  account=$ACCOUNT}${BATCH:+  batch=$BATCH}  resources: $RES"
case "$1" in
  logs) ls -t logs/*.log 2>/dev/null | head -20; exit 0 ;;
  check) sub -m ae -N ag_check "" pbs/00_check_env.pbs ;;
  pilot)
    W="-m ae -l walltime=01:30:00"  # short walltime starts sooner; one email per single job
    sub $W -N ag_pilot_cm  "PREFIX=pilot_,SITES=cm,NCHUNKS=1,LIMIT=30" pbs/11_perturb.pbs
    sub $W -N ag_pilot_fb  "PREFIX=pilot_,SITES=fb,NCHUNKS=1,LIMIT=30" pbs/11_perturb.pbs
    sub $W -N ag_pilot_ins "PREFIX=pilot_,LIMIT=3" pbs/12_insertion.pbs
    sub $W -N ag_pilot_mpra "PREFIX=pilot_,NCHUNKS=1,LIMIT=150" pbs/13_mpra.pbs ;;
  full)
    sub -J 0-39 -N ag_perturb_cm "SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs
    sub -J 0-39 -N ag_perturb_fb "SITES=fb,NCHUNKS=40" pbs/11_perturb.pbs
    sub -m ae -N ag_insertion "" pbs/12_insertion.pbs
    sub -J 0-3 -N ag_mpra "NCHUNKS=4" pbs/13_mpra.pbs ;;
  *) echo "usage: [GPU=H200|A100|...|any] [BATCH=N] [DRY=1] bash submit_all.sh check|pilot|full"; exit 1 ;;
esac
