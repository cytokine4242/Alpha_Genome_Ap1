#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh benchmark        # times batch size / 1bp decoder / heads on one GPU (GPU=L40S|A100|H100 to choose); read the log
#   bash submit_all.sh pilot            # 30 CM sites, 30 FB sites, tiny insertion, 150 MPRA elements, with CAGE (pilot_* outputs)
#   bash submit_all.sh validate         # native predictions only (1 pass/site) for 1000 CM + 1000 FB sites, mouse AND human tracksets -> results/validate_*
#   FT=1 bash submit_all.sh <mode>   the same experiment read with the fine-tuned 1 bp ATAC head: config_ft.json (config_ft_cm_only.json for scan/titrate),
#                                    outputs in results_ft/ with the same file names, job names prefixed ft_. AG_FT_CHECKPOINT=<.delta.pth> overrides the head.
#   FT=1 bash submit_all.sh ft_eval  native predictions of the fine-tuned head at every CM and FB site + equivalence check -> results_ft/eval_*  (do this first)
#   FT=1 bash submit_all.sh fb_priority   the three fibroblast priorities in one go: shuffle arm (phase1_fb), tile scan on
#                                         data/fb_subset_scan.csv, and the element swaps + host insertions below
#   FT=1 bash submit_all.sh fb_swap_all   fibroblast element swaps: every donor element into every recipient location (170 bp and peak widths,
#                                         250 sites = 50 per binding class) + the same elements in quiet host windows -> results_ft/swap_fb_*, bg_fb_*
#   bash submit_all.sh rerun_all     every experiment the analyses read: rerun_sites (phase1 + phase2 + controls + mpra_neutral) then rerun_dynamics (scan/titrate + rest); use with FT=1
#   Extra experiments (cardiomyocyte sites unless noted; all resumable):
#     controls       offset (matched-control) shuffles + width titration of the shuffle     -> results/controls_*   (61 passes/site)
#     satmut         saturation mutagenesis of the 11 bp motif (33 substitutions)           -> results/sat_*        (34 passes/site)
#     mpra_neutral   MPRA elements inserted into neutral backgrounds                        -> results/mpra_neutral_* (22 passes/element)
#     scan           tile scan (10 bp tiles, +-300 bp) + saturation mutagenesis (+-30 bp) on 1,470 CM sites (Dynamic, matched Stable, high-accessibility Stable) -> results/scan_cm_*     (364 passes/site)
#     titrate        AP-1 core weakening 1-7 bases, upgrade to a perfect core, 11 bp motif reverse complement, half-site, and motif move/duplicate at 8 offsets with controls (same sites)     -> results/titrate_cm_*   (118 passes/site)
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
#   MEM=62gb      host RAM per GPU (default 62gb; see RES below)
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
RES="-l select=1:ncpus=6:mem=${MEM:-62gb}:ngpus=1${GPU_SEL}"   # command-line resources override the #PBS lines in pbs/*.pbs. 62gb per GPU keeps jobs off V100s (all below 62 GB/GPU)
#   without dropping any A100/L40S/H200 node (96gb was rounded up to >=120 GB/GPU nodes); these jobs peak at about 5 GB
if [ "$1" = context ] && [ -z "$BATCH" ]; then BATCH=2; fi
if [ -z "$BATCH" ] && [ "$GPU" != "H200" ] && [ "$GPU" != "GH200" ]; then BATCH=4; fi
BATCH_VAR=${BATCH:+,AG_BATCH_SIZE=$BATCH}
O="-o $PWD/logs/pbs/${ACCOUNT:+ -A $ACCOUNT}"
case "$1" in scan*|titrate*|dynamics*|rerun_dynamics) CFG=${CONFIG:-config_cm_only.json} ;; full|phase1*|phase2*|controls|satmut|mpra_neutral|context|rerun_sites|rerun_all) CFG=${CONFIG:-config.json} ;; validate|human_sens) CFG=${CONFIG:-config_validate.json} ;; *) CFG=${CONFIG:-config_pilot.json} ;; esac
case "$1" in ft_eval*) FT=1 ;; esac
if [ -n "$FT" ]; then
  case "$1" in scan*|titrate*|dynamics*|rerun_dynamics) CFG=${CONFIG:-config_ft_cm_only.json} ;; *) CFG=${CONFIG:-config_ft.json} ;; esac
  BATCH_VAR="${BATCH_VAR},RESULTS_DIR=results_ft${AG_FT_CHECKPOINT:+,AG_FT_CHECKPOINT=$AG_FT_CHECKPOINT}"
fi
[ -f "$CFG" ] || { echo "config $CFG not found"; exit 1; }
BATCH_VAR="${BATCH_VAR},AG_CONFIG=$PWD/$CFG"

# sub <extra qsub args...> -- <comma-separated -v vars (may be empty)> <pbs file>
sub() {
  local vars="${@: -2:1}" file="${@: -1}"
  local args=("${@:1:$#-2}")
  if [ -n "$FT" ]; then for i in "${!args[@]}"; do [ "${args[$i]}" = "-N" ] && args[$((i+1))]="ft_${args[$((i+1))]#ag_}"; done; fi   # ft_ job names
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
p4_scan()     { sub -J 0-39 -l walltime=05:00:00 -N ag_scan_cm "EXPERIMENT=scan,PREFIX=scan_,NCHUNKS=40" pbs/15_scan_titrate.pbs; }
# the remaining 2,209 cardiomyocyte sites, so the scan covers every CM AP-1 element rather than the balanced 1,470 subset (about 76 s/site, 70 min per chunk)
p4_scan_rest()    { sub -J 0-39 -l walltime=04:00:00 -N ag_scan_rest "EXPERIMENT=scan,PREFIX=scan_rest_,NCHUNKS=40,SITES_FILE=../data/cm_subset_rest.csv" pbs/15_scan_titrate.pbs; }
p4_titrate_rest() { sub -J 0-39 -l walltime=02:00:00 -N ag_titrate_rest "EXPERIMENT=titrate,PREFIX=titrate_rest_,NCHUNKS=40,SITES_FILE=../data/cm_subset_rest.csv" pbs/15_scan_titrate.pbs; }
p4_titrate()  { sub -J 0-39 -l walltime=02:00:00 -N ag_titrate_cm "EXPERIMENT=titrate,PREFIX=titrate_,NCHUNKS=40" pbs/15_scan_titrate.pbs; }
p4_pilot()    { for e in scan titrate; do sub -m ae -l walltime=01:00:00 -N ag_${e}_pilot "EXPERIMENT=$e,PREFIX=pilot_${e}_,NCHUNKS=1,LIMIT=6" pbs/15_scan_titrate.pbs; done; }
p3_mpran()    { sub -J 0-7 -N ag_mpra_neutral "NCHUNKS=8" pbs/14_mpra_neutral.pbs; }
# fibroblast element swaps (250 sites, 50 per binding class): every donor element into every recipient location at 170 bp and peak width,
# and the same elements into quiet host windows. Needs the fine-tuned head (FT=1) and, for the backgrounds, closed_backgrounds.csv.
BG_FILE=${BG_FILE:-/srv/scratch/wonglab/jack/AlphaGenome_Ap1/finetune_data/closed_backgrounds.csv}
p5_swap()     { sub -J 0-39 -N ag_fb_swap "EXPERIMENT=swap,PREFIX=swap_,NCHUNKS=40" pbs/18_fb_swap.pbs; }
p5_bg()       { sub -J 0-7 -N ag_fb_bg "EXPERIMENT=background,PREFIX=bg_,NCHUNKS=8,BACKGROUNDS=$BG_FILE" pbs/18_fb_swap.pbs; }
# fibroblast tile scan: all 1,180 Dynamic sites + 1,180 Stable drawn at random (data/fb_subset_scan.csv), 364 passes/site
FB_SCAN_CHUNKS=${FB_SCAN_CHUNKS:-80}
# is AP-1 special? the site's own motif and each footprinted partner shuffled in the same window, with matched controls (42 passes/site)
p6_pshuf()    { sub -J 0-19 -N ag_pshuf_fb "NCHUNKS=20" pbs/19_partner_shuffle.pbs; }
p6_titrate()  { sub -J 0-39 -l walltime=02:00:00 -N ag_titrate_fb "EXPERIMENT=titrate,PREFIX=titrate_,SITES=fb,NCHUNKS=40,SITES_FILE=../data/fb_subset_scan.csv" pbs/15_scan_titrate.pbs; }
p5_fbscan()   { sub -J 0-$((FB_SCAN_CHUNKS-1)) -l walltime=02:00:00 -N ag_scan_fb "EXPERIMENT=scan,PREFIX=scan_,SITES=fb,NCHUNKS=$FB_SCAN_CHUNKS,SITES_FILE=../data/fb_subset_scan.csv" pbs/15_scan_titrate.pbs; }
# the same 10 bp edit, but the window SLIDES 2 bp at a time from -400 to +390, so the profile is sampled every
# 2 bp across the whole element instead of once per 10 bp tile. 396 windows + native = 397 passes/site, about
# 116 s/site, so ~20 sites a chunk finishes inside 45 min. One replicate per position: every 10 bp of sequence is
# already covered by five overlapping windows.
FB_SLIDE_CHUNKS=${FB_SLIDE_CHUNKS:-120}
p5_fbslide()  { sub -J 0-$((FB_SLIDE_CHUNKS-1)) -l walltime=00:30:00 -N ag_slide_fb "EXPERIMENT=slide,PREFIX=slide_,SITES=fb,NCHUNKS=$FB_SLIDE_CHUNKS,SITES_FILE=../data/fb_subset_scan.csv" pbs/15_scan_titrate.pbs; }
p5_slide_pilot() { sub -m ae -l walltime=01:00:00 -N ag_slide_fb_pilot "EXPERIMENT=slide,PREFIX=pilot_slide_,SITES=fb,NCHUNKS=1,LIMIT=3,SITES_FILE=../data/fb_subset_scan.csv" pbs/15_scan_titrate.pbs; }
p5_pilot()    { sub -m ae -l walltime=01:00:00 -N ag_scan_fb_pilot "EXPERIMENT=scan,PREFIX=pilot_,SITES=fb,NCHUNKS=1,LIMIT=4,SITES_FILE=../data/fb_subset_scan.csv" pbs/15_scan_titrate.pbs
                sub -m ae -l walltime=01:00:00 -N ag_fb_swap_pilot "EXPERIMENT=swap,PREFIX=pilot_swap_,NCHUNKS=1,LIMIT=10" pbs/18_fb_swap.pbs
                sub -m ae -l walltime=01:00:00 -N ag_fb_bg_pilot "EXPERIMENT=background,PREFIX=pilot_bg_,NCHUNKS=1,LIMIT=6,BACKGROUNDS=$BG_FILE" pbs/18_fb_swap.pbs; }
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
  ft_eval)       sub -J 0-3 -N ag_eval "NCHUNKS=4" pbs/17_ft_evaluate.pbs ;;
  ft_eval_pilot) sub -l walltime=01:00:00 -N ag_eval_pilot "PREFIX=pilot_,NCHUNKS=1,LIMIT=24" pbs/17_ft_evaluate.pbs ;;
  rerun_all)     bash "$0" rerun_sites && bash "$0" rerun_dynamics ;;   # two calls so each half gets its own config
  rerun_sites)   p1_cm; p1_fb; p1_other; p2_cm; p2_fb; p3_controls; p3_mpran ;;
  rerun_dynamics) p4_scan; p4_titrate; p4_scan_rest; p4_titrate_rest ;;
  fb_partner)    p6_pshuf ;;
  fb_titrate)    p6_titrate ;;
  fb_slide)      p5_fbslide ;;
  fb_slide_pilot) p5_slide_pilot ;;
  fb_overnight)  p6_pshuf; p6_titrate ;;
  fb_scan)       p5_fbscan ;;
  fb_priority)   p1_fb; p5_fbscan; p5_swap; p5_bg ;;   # the three fibroblast priorities: shuffle, tile scan, element swaps + host insertions
  fb_swap)       p5_swap ;;
  fb_background) p5_bg ;;
  fb_swap_all)   p5_swap; p5_bg ;;
  fb_swap_pilot) p5_pilot ;;
  scan)          p4_scan ;;
  titrate)       p4_titrate ;;
  dynamics)      p4_scan; p4_titrate ;;
  scan_rest)     p4_scan_rest ;;
  titrate_rest)  p4_titrate_rest ;;
  dynamics_rest) p4_scan_rest; p4_titrate_rest ;;
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
  *) echo "usage: [GPU=H200|A100|...|any] [BATCH=N] [DRY=1] bash submit_all.sh check|benchmark|pilot|validate|ft_eval|ft_eval_pilot|fb_scan|fb_priority|fb_partner|fb_titrate|fb_overnight|fb_swap|fb_background|fb_swap_all|fb_swap_pilot|rerun_all|rerun_sites|rerun_dynamics|controls|satmut|scan|titrate|dynamics|dynamics_pilot|scan_rest|titrate_rest|dynamics_rest|mpra_neutral|human_sens|context|phase1[_cm|_fb|_other]|phase2[_cm|_fb]|full"; exit 1 ;;
esac
