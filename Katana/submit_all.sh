#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh pilot            # 10 CM sites, 10 FB sites, tiny insertion, 20 MPRA elements
#   bash submit_all.sh full             # everything: 40 chunks CM, 40 chunks FB, insertion, 4 MPRA chunks
# Logs go to logs/, results to results/.
set -e
cd "$(dirname "$0")"
mkdir -p logs results
O="-o $PWD/logs/"
case "$1" in
  check) qsub $O pbs/00_check_env.pbs ;;
  pilot)
    W="-l walltime=01:00:00"  # short walltime starts sooner
    qsub $O $W -v SITES=cm,NCHUNKS=1,LIMIT=10 -N ag_pilot_cm pbs/11_perturb.pbs
    qsub $O $W -v SITES=fb,NCHUNKS=1,LIMIT=10 -N ag_pilot_fb pbs/11_perturb.pbs
    qsub $O $W -v LIMIT=3 -N ag_pilot_ins pbs/12_insertion.pbs
    qsub $O $W -v NCHUNKS=1,LIMIT=20 -N ag_pilot_mpra pbs/13_mpra.pbs ;;
  full)
    qsub $O -J 0-39 -v SITES=cm,NCHUNKS=40 pbs/11_perturb.pbs
    qsub $O -J 0-39 -v SITES=fb,NCHUNKS=40 pbs/11_perturb.pbs
    qsub $O pbs/12_insertion.pbs
    qsub $O -J 0-3 -v NCHUNKS=4 pbs/13_mpra.pbs ;;
  *) echo "usage: bash submit_all.sh check|pilot|full"; exit 1 ;;
esac
