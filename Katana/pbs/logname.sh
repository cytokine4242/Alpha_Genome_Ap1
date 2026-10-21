#!/bin/bash
# Sourced by common.sh: tee everything this job prints into a readable log file.
#   logs/<YYYYMMDD-HHMM>_<jobname>[_chunkNNofMM]_<GPU>_<jobid>[_<arrayindex>].log
# PBS's own <jobid>.<server>.OU file (with the Katana resource summary) is kept separately in logs/pbs/.
TS=$(date +%Y%m%d-%H%M)
JID=$(echo "${PBS_JOBID:-local}" | sed 's/\..*//; s/\[\]//; s/\[\([0-9]*\)\]/_\1/')   # 123[5].kman... -> 123_5
GPUNAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 | sed 's/NVIDIA //; s/[ -].*//')
CH=""
if [ -n "$PBS_ARRAY_INDEX" ] && [ -n "$NCHUNKS" ]; then CH=_chunk$(printf %02d "$PBS_ARRAY_INDEX")of$(printf %02d "$NCHUNKS"); fi
export LOG_FILE="$ROOT/logs/${TS}_${PBS_JOBNAME:-job}${CH}_${GPUNAME:-nogpu}_${JID}.log"
mkdir -p "$ROOT/logs"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "log: $LOG_FILE"
