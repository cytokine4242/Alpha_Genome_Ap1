#!/bin/bash
# Sourced at the start of every PBS job (common.sh, Finetune/pbs/*.pbs). Appends a "start" row (submission -> start wait) and, when the job exits,
# an "end" row (run time, exit status) to one shared TSV, so queue waits survive after PBS forgets the job (about a day). Never fails the job.
#   columns: event time job name queue account walltime_req select qtime stime wait_min host gpu queued_ours queued_in_queue exit run_min
#   queued_ours = our queued jobs/arrays when this one started; queued_in_queue = jobs queued in this job's queue then (all users); run_min counts from stime
QW_LOG=${QW_LOG:-/srv/scratch/wonglab/jack/AlphaGenome_Ap1/queue_waits.tsv}
QW_HEADER="event	time	job	name	queue	account	walltime_req	select	qtime	stime	wait_min	host	gpu	queued_ours	queued_in_queue	exit	run_min"

_qw_append() {   # one tab-separated row; flock keeps rows from concurrent array tasks whole
  {
    flock -w 10 9 2>/dev/null
    [ -s "$QW_LOG" ] || echo "$QW_HEADER" >&9
    echo "$1" >&9
  } 9>>"$QW_LOG"
}

_qw_start() {
  [ -n "$PBS_JOBID" ] || return 0
  local info get qt st now wait q
  info=$(qstat -f -F dsv "$PBS_JOBID" 2>/dev/null | tr '|' '\n')
  get() { sed -n "s/^$1=//p" <<<"$info" | head -1; }
  qt=$(get qtime); st=$(get stime); q=$(get queue)
  now=$(date +%s)
  QW_T0=$( [ -n "$st" ] && date -d "$st" +%s || echo "$now")
  wait=$( [ -n "$qt" ] && awk -v a="$(date -d "$qt" +%s)" -v b="$QW_T0" 'BEGIN{printf "%.1f", (b - a) / 60}')
  QW_ROW="$PBS_JOBID	$(get Job_Name)	$q	$(get Account_Name)	$(get Resource_List.walltime)	$(get Resource_List.select)	$qt	$st"
  _qw_append "start	$(date '+%F %T')	$QW_ROW	$wait	$(hostname -s)	$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)	$(qstat -u "$USER" 2>/dev/null | awk '$10 == "Q"' | wc -l)	$(qstat -Q "$q" 2>/dev/null | awk 'NR == 3 {print $6}')		"
  echo "queue wait: ${wait:-?} min (submitted $qt, started $st, queue $q)"
}

_qw_end() {
  local rc=$?
  _qw_append "end	$(date '+%F %T')	$QW_ROW		$(hostname -s)				$rc	$(awk -v a="$QW_T0" -v b="$(date +%s)" 'BEGIN{printf "%.1f", (b - a) / 60}')"
}

_qw_start 2>/dev/null || true
[ -n "$PBS_JOBID" ] && trap '_qw_end 2>/dev/null || true' EXIT
