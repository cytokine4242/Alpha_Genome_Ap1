#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh benchmark        # times batch size / 1bp decoder / heads on one GPU (GPU=L40S|A100|H100 to choose); read the log
#   bash submit_all.sh pilot            # 30 CM sites, 30 FB sites, tiny insertion, 150 MPRA elements, with CAGE (pilot_* outputs)
#   bash submit_all.sh validate         # native predictions only (1 pass/site) for 1000 CM + 1000 FB sites, mouse AND human tracksets -> results/validate_*
#   Extra experiments (cardiomyocyte sites unless noted; all resumable):
#     controls       offset (matched-control) shuffles + width titration of the shuffle     -> results/controls_*   (61 passes/site)
#     satmut         saturation mutagenesis of the 11 bp motif (33 substitutions)           -> results/sat_*        (34 passes/site)
#     mpra_neutral   MPRA elements inserted into neutral backgrounds                        -> results/mpra_neutral_* (22 passes/element)
#     scan           tile scan (10 bp tiles, +-300 bp) + saturation mutagenesis (+-30 bp) on 1,470 CM sites (Dynamic, matched Stable, high-accessibility Stable) -> results/scan_cm_*     (364 passes/site)
#     titrate        AP-1 core weakening 1-7 bases, consensus upgrade, reverse complement, half-site, and motif move/duplicate at 8 offsets with controls (same sites)     -> results/titrate_cm_*   (94 passes/site)
#     dynamics       scan + titrate together (cardiomyocyte tracks only, config_cm_only.json); test first with: bash submit_all.sh dynamics_pilot
#     human_sens     shuffles with human output heads on 500 sites (config_validate.json)   -> results/human_*
#     context        shuffles on 60 sites with 16 kb / 65 kb / 262 kb windows               -> results/ctx<window>_*
#   bash submit_all.sh phase1           # shuffle arm on all sites (40+40 chunks) + insertion + MPRA   -> results/shuffle_*, mpra_*, insertion_*
#   (phase1_cm | phase1_fb | phase1_other, phase2_cm | phase2_fb submit the parts separately, e.g. to spread over queues:
#      bash submit_all.sh phase1_cm;  ACCOUNT= GPU=H200 bash submit_all.sh phase1_fb)
#   bash submit_all.sh phase2           # mutation arm (random, dinuc, core, flank) on all sites       -> results/mut_*
#   bash submit_all.sh full             # phase1 + phase2 together (every edit kind in one go)         -> results/cm_*, fb_*
#
# Options (environment variables):
#   GPU=H200 (default) | A100 | L40S | GH200 | ... | any     GPU model; "any" drops the gpu_model request so the
#                                                           scheduler uses whichever GPU is free first
#   ACCOUNT=...   PBS account (-A). Default mwac_wonglab (the lab's private GPUs); ACCOUNT= (empty) submits without one.
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
ACCOUNT=${ACCOUNT-mwac_wonglab}   # default account for all jobs; run with ACCOUNT= (empty) to submit without one
mkdir -p logs/pbs results

if [ -n "$ACCOUNT" ]; then GPU=${GPU:-any}; else GPU=${GPU:-H200}; fi
if [ "$GPU" = "any" ]; then GPU_SEL=""; else GPU_SEL=":gpu_model=$GPU"; fi
RES="-l select=1:ncpus=6:mem=96gb:ngpus=1${GPU_SEL}"   # command-line resources override the #PBS lines in pbs/*.pbs
if [ "$1" = context ] && [ -z "$BATCH" ]; then BATCH=2; fi
if [ -z "$BATCH" ] && [ "$GPU" != "H200" ] && [ "$GPU" != "GH200" ]; then BATCH=4; fi
BATCH_VAR=${BATCH:+,AG_BATCH_SIZE=$BATCH}
O="-o $PWD/logs/pbs/${ACCOUNT:+ -A $ACCOUNT}"
case "$1" in scan|titrate|dynamics*) CFG=${CONFIG:-config_cm_only.json} ;; full|phase1*|phase2*|controls|satmut|mpra_neutral|context) CFG=${CONFIG:-config.json} ;; validate|human_sens) CFG=${CONFIG:-config_validate.json} ;; *) CFG=${CONFIG:-config_pilot.json} ;; esac
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

p1_cm()    { sub -J 0-39 -N ag_shuffle_cm "PREFIX=shuffle_,KINDS=shuffle,SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs; }
p1_fb()    { sub -J 0-39 -N ag_shuffle_fb "PREFIX=shuffle_,KINDS=shuffle,SITES=fb,NCHUNKS=40" pbs/11_perturb.pbs; }
p1_other() { sub -m ae -N ag_insertion "" pbs/12_insertion.pbs; sub -J 0-7 -N ag_mpra "NCHUNKS=8" pbs/13_mpra.pbs; }
p2_cm()    { sub -J 0-39 -N ag_mut_cm "PREFIX=mut_,KINDS=random:dinuc:core:flank,SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs; }
p2_fb()    { sub -J 0-39 -N ag_mut_fb "PREFIX=mut_,KINDS=random:dinuc:core:flank,SITES=fb,NCHUNKS=40" pbs/11_perturb.pbs; }

p3_controls() { sub -J 0-39 -N ag_controls_cm "PREFIX=controls_,KINDS=offset:titrate,SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs; }
p3_satmut()   { sub -J 0-39 -N ag_satmut_cm "PREFIX=sat_,KINDS=sat,SITES=cm,NCHUNKS=40" pbs/11_perturb.pbs; }
p4_scan()     { sub -J 0-39 -N ag_scan_cm "EXPERIMENT=scan,PREFIX=scan_,NCHUNKS=40" pbs/15_scan_titrate.pbs; }
p4_titrate()  { sub -J 0-39 -N ag_titrate_cm "EXPERIMENT=titrate,PREFIX=titrate_,NCHUNKS=40" pbs/15_scan_titrate.pbs; }
p4_pilot()    { for e in scan titrate; do sub -m ae -l walltime=01:00:00 -N ag_${e}_pilot "EXPERIMENT=$e,PREFIX=pilot_${e}_,NCHUNKS=1,LIMIT=6" pbs/15_scan_titrate.pbs; done; }
p3_mpran()    { sub -J 0-7 -N ag_mpra_neutral "NCHUNKS=8" pbs/14_mpra_neutral.pbs; }
p3_human()    { sub -J 0-4 -N ag_human_sens "PREFIX=human_,KINDS=shuffle,SITES=cm,NCHUNKS=5,LIMIT=100" pbs/11_perturb.pbs; }
p3_context()  { for w in 16384 65536 262144; do sub -N ag_ctx_$w "PREFIX=ctx${w}_,KINDS=shuffle,SITES=cm,NCHUNKS=1,LIMIT=60,AG_WINDOW=$w" pbs/11_perturb.pbs; done; }

echo "config=$CFG  GPU=$GPU${ACCOUNT:+  account=$ACCOUNT}${BATCH:+  batch=$BATCH}  resources: $RES"
case "$1" in
  logs) ls -t logs/*.log 2>/dev/null | head -20; exit 0 ;;
  benchmark) sub -m ae -N ag_benchmark "" pbs/00_benchmark.pbs ;;
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
    sub -J 0-7 -N ag_mpra "NCHUNKS=8" pbs/13_mpra.pbs ;;
  validate)
    W="-m ae -l walltime=02:00:00"
    sub $W -N ag_validate_cm "PREFIX=validate_,KINDS=native,SITES=cm,NCHUNKS=1,LIMIT=1000" pbs/11_perturb.pbs
    sub $W -N ag_validate_fb "PREFIX=validate_,KINDS=native,SITES=fb,NCHUNKS=1,LIMIT=1000" pbs/11_perturb.pbs ;;
  scan)          p4_scan ;;
  titrate)       p4_titrate ;;
  dynamics)      p4_scan; p4_titrate ;;
  dynamics_pilot) p4_pilot ;;
  controls)      p3_controls ;;
  satmut)        p3_satmut ;;
  mpra_neutral)  p3_mpran ;;
  human_sens)    p3_human ;;
  context)       p3_context ;;
  phase1)        p1_cm; p1_fb; p1_other ;;
  phase1_cm)     p1_cm ;;
  phase1_fb)     p1_fb ;;
  phase1_other)  p1_other ;;
  phase2)        p2_cm; p2_fb ;;
  phase2_cm)     p2_cm ;;
  phase2_fb)     p2_fb ;;
  *) echo "usage: [GPU=H200|A100|...|any] [BATCH=N] [DRY=1] bash submit_all.sh check|benchmark|pilot|validate|controls|satmut|scan|titrate|dynamics|dynamics_pilot|mpra_neutral|human_sens|context|phase1[_cm|_fb|_other]|phase2[_cm|_fb]|full"; exit 1 ;;
esac
