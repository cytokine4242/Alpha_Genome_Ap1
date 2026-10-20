#!/bin/bash
# Sourced by every PBS job. Submit jobs from the Katana/ folder (use submit_all.sh) so PBS_O_WORKDIR is that folder.
ENV_NAME=alphagenome
source ~/.bashrc
conda activate $ENV_NAME
export ROOT=${PBS_O_WORKDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$ROOT/scripts"
mkdir -p "$ROOT/results" "$ROOT/logs"
