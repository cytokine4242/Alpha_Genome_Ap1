# Katana jobs: AlphaGenome AP-1 perturbations (mouse, mm10)

Run everything from inside this `Katana/` folder on Katana (clone the repo or drag the folder). Nothing has been run on
Katana or against the real model yet: the model calls follow the alphagenome-pytorch docs, and the logic (edits, readouts,
output tables) was tested locally with a mock model. Do the check step first.

## Folder
| Path | What |
|---|---|
| `config.json` | paths (genome, GTF, weights), trackset regexes, replicate counts. The only file you should need to edit. |
| `data/` | `cm_fos_sites.csv`, `fb_fos_sites.csv` (all TOBIAS Fos sites with the labels below), `MPRA_Results_with_TSS_distance.csv`, `fb_mpra_tobias_crossmap.csv` (MPRA and TOBIAS binding definitions side by side), `measured_atac_sites.csv` (snATAC CPM/kb per site, for local validation plots only). |
| `scripts/` | `ag_lib.py` (shared), `00_check_env.py`, `11_perturb_sites.py`, `12_insertion.py`, `13_mpra.py` |
| `pbs/` | PBS job files (Katana H200 header as in the demo script); `common.sh` activates the conda env |
| `submit_all.sh` | `check`, `pilot`, `full` |
| `results/`, `logs/` | outputs (copy `results/` back for replotting) |

## Run order
1. **Environment** (once): `bash setup_env.sh`, or set `ENV_NAME` in `pbs/common.sh` to an existing alphagenome env.
2. **Check**: `bash submit_all.sh check`, then read `logs/ag_check_env.o*`. It confirms the genome, GTF and weights paths, loads the weights,
   writes `results/mouse_tracks.csv` (all mouse tracks), shows which tracks each trackset regex matched, and times a forward pass.
   Fix `config.json` (weights path, `tracksets`) until it passes. The weights path is a guess:
   `/srv/scratch/wonglab/minh/alphagenome/weights/model_all_folds.safetensors`.
3. **Pilot**: `bash submit_all.sh pilot` (10 CM + 10 FB sites, 3 insertion backgrounds, 20 MPRA elements). Check `results/` and s/site in the logs.
4. **Full**: `bash submit_all.sh full` (40 CM + 40 FB array tasks at 12 h walltime, 4 MPRA tasks, 1 insertion job; each chunk is resumable, rerun the same command to continue).

## Choosing the GPU
`submit_all.sh` takes environment variables: `GPU=H200` (default), `GPU=A100`, `GPU=L40S`, ... or `GPU=any` (no `gpu_model` request, so whichever GPU frees up first),
`BATCH=N` (sequences per GPU batch; defaults to the config value on H200/GH200 and 4 on anything else), and `DRY=1` (print the qsub commands only).
Examples: `GPU=any bash submit_all.sh pilot`, `GPU=A100 bash submit_all.sh full`. The GPU each job landed on is printed at the top of its log.
Avoid V100 (no native bf16): with `GPU=any` check the logs and resubmit chunks that landed on one if they fail.

## Site labels (in every site table and output)
Three views of the same binding calls, so ageing, exercise and their combination can be tested separately.
YS = young sedentary, MAS = aged sedentary, MAE = aged + exercise.
- `age_class` (YS vs MAS): Stable (bound in both) / Age Up (MAS only) / Age Down (YS only) / Unbound (neither)
- `ex_class` (MAS vs MAE): Stable (bound in both) / Exercise Up (MAE only) / Exercise Down (MAS only) / Unbound (neither)
- `binding_category`: the full pattern, one of YS+MAS+MAE, YS+MAS, YS+MAE, YS only, MAS+MAE, MAS only, MAE only
- `stable_dynamic`: Stable = bound in all three, Dynamic = anything else

A site carries all of them (e.g. "MAS only" is Age Up and Exercise Down).

## What each job writes (`results/`)
- `cm_chunk*.csv`, `fb_chunk*.csv` (11): long table, one row per site x condition x replicate x trackset. Columns: the labels above, `condition`
  (`native`, `shuffle|random|dinuc@motif|170|peak`, `core@motif`, `flank@motif`, `flank@170`), `trackset` (CM or FB), readouts
  `atac_motif/atac_170/atac_peak`, `rna_prox2kb` (nearest TSS +-2 kb, NaN if the TSS is outside the window), `rna_mean50kb`, `rna_max50kb`, `tss_dist`.
  Log2FC = log2((edit + pc)/(native + pc)) and dependency = -log2FC are computed at plotting time. Sites that are off-chromosome or have N in
  the motif are listed in `*.skipped.txt`.
- `insertion_peakswap.csv`, `insertion_motif.csv` (12), `mpra_chunk*.csv` (13; includes `motif_source`: tobias / regex / none).

## Design choices to know about
- Window 131,072 bp centred on the motif midpoint; ATAC at 1 bp, RNA at 128 bp; multiple matched tracks are averaged (RNA strands included).
- Motif = 11 bp TOBIAS interval; "170 bp" = motif midpoint +-85; "peak" = the TOBIAS peak. Core = `TGA[CG]TCA` inside the motif (central 7 bp if absent).
- Core edit: transversion (A<->C, G<->T). Flank edits mutate every flank base (motif flanks = motif minus core; 170 flanks = 170 bp minus motif). Edits are seeded per site.
- Cost: 72 forward passes per site with default reps (10 shuffle, 5 random, 5 dinuc, 5 flank), 8,697 sites. The check job prints a GPU-hour
  estimate; lower `perturb.reps` in `config.json` if it is too large.
- Peak swap: recipients are CM sites at the 10/40/70/95th percentile of native peak ATAC; donors are 50 Stable + 25 Dynamic. Motif insertion:
  50 random mm10 backgrounds, dinucleotide-shuffled, TGAGTCA inserted at the centre, shuffled-motif control.
- The slide site counts (2,663 CM / 2,320 Stable) are not reproduced by the TOBIAS table: these jobs use all 3,679 CM and 5,018 fibroblast Fos sites.
