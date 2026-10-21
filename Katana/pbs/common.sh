#!/bin/bash
# Sourced by every PBS job. Submit jobs from the Katana/ folder (use submit_all.sh) so PBS_O_WORKDIR is that folder.
source ~/.bashrc
export ROOT=${PBS_O_WORKDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
# env path written by setup_env.sh; falls back to a named conda env
if [ -f "$ROOT/.env_prefix" ]; then conda activate "$(cat "$ROOT/.env_prefix")"; else conda activate alphagenome2; fi
cd "$ROOT/scripts"
mkdir -p "$ROOT/results" "$ROOT/logs"
source "$ROOT/pbs/logname.sh"
echo "job $PBS_JOBID | $(qstat -f $PBS_JOBID 2>/dev/null | egrep -i 'queue =|Account_Name|project =' | tr -s ' \n' ' ')"
echo "host: $(hostname) | GPU: $(nvidia-smi -L 2>/dev/null | head -1) | batch override: ${AG_BATCH_SIZE:-none}"
echo "python: $(which python)"
