#!/bin/bash
# Run from Finetune/ on Katana once the 1 bp probe has finished. Freezes the chosen epoch as the head the configs point at,
# then submits the evaluation (and, with "all", the three fibroblast priorities).
#   bash after_training.sh            # freeze the head + submit ft_eval, then print the next command
#   bash after_training.sh all        # also submit FT=1 fb_priority straight away
#   RUN=probe_1bp_fold_3_128kb bash after_training.sh        # a different run directory
set -e
cd "$(dirname "$0")"
FTDIR=/srv/scratch/wonglab/jack/AlphaGenome_Ap1/finetune_data
RUN=${RUN:-probe_1bp_fold_3_128kb}
SRC=$FTDIR/runs/$RUN/best_model.delta.pth
DST=$FTDIR/heads/probe_1bp_fold3.delta.pth          # the path config_ft.json and config_ft_cm_only.json name

[ -f "$SRC" ] || { echo "no checkpoint at $SRC: has the run finished? (qstat, and $FTDIR/runs/$RUN/epoch_log.csv)"; exit 1; }
echo "== epochs"
column -s, -t "$FTDIR/runs/$RUN/epoch_log.csv" | cut -c1-150
BEST=$(python3 - "$FTDIR/runs/$RUN/epoch_log.csv" <<'EOF'
import sys, csv
rows = list(csv.DictReader(open(sys.argv[1])))
b = min(rows, key=lambda r: float(r["val_loss"]))
print(f'epoch {b["epoch"]}, val_loss {float(b["val_loss"]):.4f}')
EOF
)
echo "== best: $BEST"

mkdir -p "$FTDIR/heads"
if [ -f "$DST" ] && ! cmp -s "$SRC" "$DST"; then
  echo "WARNING: $DST exists and differs. Move it aside first if you want to replace the head; stopping."; exit 1
fi
cp -p "$SRC" "$DST"
echo "== froze the head: $DST"

cd ../Katana
echo "== submitting ft_eval (native predictions at every CM and FB site + the equivalence check)"
FT=1 bash submit_all.sh ft_eval

if [ "$1" = all ]; then
  echo "== submitting the fibroblast priorities (shuffle arm, tile scan, element swaps, host insertions)"
  [ -f "$FTDIR/closed_backgrounds.csv" ] || echo "   note: closed_backgrounds.csv is missing, so the host-insertion array will fail; run 'bash ../Finetune/submit.sh backgrounds' and resubmit 'FT=1 bash submit_all.sh fb_background'"
  FT=1 bash submit_all.sh fb_priority
else
  echo
  echo "Next, once you have read Analysis (ft_eval results):"
  echo "   cd ../Katana && FT=1 bash submit_all.sh fb_priority     # shuffle arm + tile scan + element swaps + host insertions"
  echo "   cd ../Katana && FT=1 bash submit_all.sh rerun_all       # every other experiment the analyses read"
fi
