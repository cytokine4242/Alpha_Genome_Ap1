#!/bin/bash
# Usage (from the Katana/ folder on Katana):
#   bash submit_all.sh check            # environment + track check (do this first, read logs/)
#   bash submit_all.sh pilot            # 10 CM sites, 10 FB sites, tiny insertion, 20 MPRA elements
#   bash submit_all.sh full             # everything: 20 chunks CM, 20 chunks FB, insertion, 4 MPRA chunks
# Logs go to logs/, results to results/.
set -e
cd "$(dirname "$0")"
mkdir -p logs results
O="-o $PWD/logs/"
case "$1" in
  check) qsub $O pbs/00_check_env.pbs ;;
  pilot)
    qsub $O -v SITES=cm,NCHUNKS=1,LIMIT=10 -N ag_pilot_cm pbs/11_perturb.pbs
    qsub $O -v SITES=fb,NCHUNKS=1,LIMIT=10 -N ag_pilot_fb pbs/11_perturb.pbs
    qsub $O -v LIMIT=3 -N ag_pilot_ins pbs/12_insertion.pbs
    qsub $O -v NCHUNKS=1,LIMIT=20 -N ag_pilot_mpra pbs/13_mpra.pbs ;;
  full)
    qsub $O -J 0-19 -v SITES=cm,NCHUNKS=20 pbs/11_perturb.pbs
    qsub $O -J 0-19 -v SITES=fb,NCHUNKS=20 pbs/11_perturb.pbs
    qsub $O pbs/12_insertion.pbs
    qsub $O -J 0-3 -v NCHUNKS=4 pbs/13_mpra.pbs ;;
  *) echo "usage: bash submit_all.sh check|pilot|full"; exit 1 ;;
esac
