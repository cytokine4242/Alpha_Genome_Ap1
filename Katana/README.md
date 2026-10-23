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
| `/srv/scratch/wonglab/jack/AlphaGenome_Ap1/queue_waits.tsv` | one `start` row (submitted, started, minutes waited, queue, resources, GPU, how busy the queue was) and one `end` row (exit status, minutes run) per PBS job, written by `pbs/queue_wait.sh` (sourced by `common.sh` and the Finetune PBS files). PBS itself forgets finished jobs after about a day. |
| `results/`, `logs/` | outputs (copy `results/` back for replotting). Readable logs are named `logs/<date>_<jobname>[_chunkNNofMM]_<GPU>_<jobid>.log` (e.g. `20261021-1152_ag_pilot_cm_L40S_9578160.log`); PBS's raw `<jobid>.<server>.OU` files go to `logs/pbs/`. `bash submit_all.sh logs` lists the newest. |

## Run order
1. **Environment** (once): `bash setup_env.sh` (creates a Python 3.12 env on scratch (home quota is small), path saved in `.env_prefix`, alphagenome-pytorch from GitHub; the pip release is too old). To match Minh's code exactly: `REF=<commit> bash setup_env.sh` (commit from `git -c safe.directory='*' -C /srv/scratch/wonglab/minh/alphagenome/alphagenome-pytorch rev-parse HEAD`). Then run `python scripts/00_inspect_weights.py` once to confirm the weights match the package.
2. **Check**: `bash submit_all.sh check`, then read the newest log (`bash submit_all.sh logs`). It confirms the genome, GTF and weights paths, loads the weights,
   writes `results/mouse_tracks.csv` (all mouse tracks), shows which tracks each trackset regex matched, and times a forward pass.
   Fix `config.json` (weights path, `tracksets`) until it passes. The weights path is a guess:
   `/srv/scratch/wonglab/minh/alphagenome/weights/model_all_folds.safetensors`.
3. **Pilot**: `bash submit_all.sh pilot` (10 CM + 10 FB sites, 3 insertion backgrounds, 20 MPRA elements). Check `results/` and s/site in the logs.
4. **Phased full run** (recommended): `bash submit_all.sh phase1` = shuffle arm on all sites (31 passes/site) + insertion + MPRA, outputs `results/shuffle_*`, `mpra_*`, `insertion_*`;
   then `bash submit_all.sh phase2` = mutation arm (random, dinuc, core, flank; 42 passes/site), outputs `results/mut_*`. `bash submit_all.sh full` runs every edit kind in one go (`results/cm_*`, `fb_*`).
   Each is 40 CM + 40 FB array tasks at 2 h walltime (a chunk takes ~15-30 min on an L40S, so the 2 h queue is used and a slow GPU only lengthens the tail; the MPRA job is 8 chunks) and is resumable (rerun the same command to continue). Every output also contains the native rows.

## Choosing the GPU
`submit_all.sh` takes environment variables: `GPU=H200` (default), `GPU=A100`, `GPU=L40S`, ... or `GPU=any` (no `gpu_model` request, so whichever GPU frees up first),
`BATCH=N` (sequences per GPU batch; defaults to the config value on H200/GH200 and 4 on anything else), and `DRY=1` (print the qsub commands only).
Examples: `GPU=any bash submit_all.sh pilot`, `GPU=A100 bash submit_all.sh full`. The GPU each job landed on is printed at the top of its log.
Avoid V100 (no native bf16): with `GPU=any` check the logs and resubmit chunks that landed on one if they fail.

## Private GPUs (mwac_wonglab account)
`pbs/00_check_access.pbs` has `-A mwac_wonglab`, `-m ae`, `-j oe` and no GPU model. It tests whether the account works and what you get, then runs the normal environment check:
`qsub -o $PWD/logs/pbs/ pbs/00_check_access.pbs` (or `ACCOUNT=mwac_wonglab bash submit_all.sh check`). Read the newest `logs/*_ag_check_access_*.log`: the header lists the queue, account and GPU.
Any submit_all.sh command accepts `ACCOUNT=mwac_wonglab` (GPU defaults to `any` when an account is given; single jobs get `-m ae`, arrays do not).

## Extra experiments (all resumable, private account by default)
| Mode | What | Passes | Output |
|---|---|---|---|
| `controls` | matched-control shuffles of a 170 bp region at offsets of -600/-300/+300/+600 bp from the motif, plus centred shuffles at widths 7/11/21/41/85/170/340/680 bp | 61 per site | `results/controls_cm_*` |
| `satmut` | every base of the 11 bp motif changed to each alternative (33 substitutions) | 34 per site | `results/sat_cm_*` |
| `mpra_neutral` | each MPRA element (native, motif-scrambled, fully scrambled; 5 shuffles each) inserted into 2 fixed neutral backgrounds | 22 per element | `results/mpra_neutral_*` |
| `human_sens` | shuffle arm with human output heads on 500 random CM sites (`config_validate.json`) | 31 per site x 2 organisms | `results/human_cm_*` |
| `context` | shuffle arm on 60 sites with 16 kb, 65 kb and 262 kb windows (`AG_WINDOW`) | 31 per site | `results/ctx<window>_cm_*` |

`bash submit_all.sh <mode>`; preview first with `DRY=1`.

## Dynamic-vs-Stable sensitivity experiments (cardiomyocyte tracks only)
Two experiments on a fixed subset of 1,470 cardiomyocyte sites (`data/cm_subset_dynamics.csv`, made by `Scripts/34_make_dynamics_subsets.py`; no TOBIAS footprint scores are used, only the Stable/Dynamic labels and AlphaGenome's own predictions): **dynamic** (500 Dynamic sites), **matched_stable** (470 Stable sites matched 1:1 to those on native predicted 170 bp and peak accessibility, peak width and TSS distance) and **high_stable** (500 Stable sites from the top quartile of native predicted 170 bp accessibility). They use `config_cm_only.json` (only the heart ATAC and heart RNA tracks), so nothing fibroblast-related runs until a fibroblast head exists.

| mode | what it does | passes/site | output |
|---|---|---|---|
| `scan` | tile scan (every base of a 10 bp tile replaced; 60 tiles from -300 to +290 bp, 3 replicates) and saturation mutagenesis (3 alternatives at each of 61 positions, -30 to +30 bp) | 364 | `results/scan_cm_chunk*.csv` |
| `titrate` | core weakening (1 to 7 of the 7 core bases changed), upgrade to a perfect TGA[C/G]TCA core with the fewest changes (only at sites without one), reverse complement of the whole 11 bp motif (the core itself is a palindrome), half-site only; the 11 bp motif duplicated or moved to offsets -150, -80, -40, -20, +20, +40, +80, +150, each with its own shuffled control (one for duplicate, one for move) | 118 | `results/titrate_cm_chunk*.csv` |

```
bash submit_all.sh dynamics_pilot      # 6 sites per experiment, 1 h walltime: do this first and check logs/ and results/pilot_*
bash submit_all.sh scan                # 40-chunk array, about 37 sites per chunk, about 1 h per chunk at 0.25 s/pass
bash submit_all.sh titrate             # 40-chunk array, about 15 min per chunk on an A100, 2 h walltime so it can use any free GPU on the cluster
bash submit_all.sh dynamics            # both
```
Total about 42 GPU-hours for the scan and 13 for the titration; with the private account's ~13 concurrent tasks that is about 3 to 4 hours of wall time. Both are resumable (sites already in the output are skipped). `python scripts/15_scan_titrate.py --selftest` checks every edit on synthetic sequence without a GPU. Analysis (after copying `results/` back): `python Scripts/35_scan_titrate_analysis.py` writes `Analysis/scan_titrate/`.
Each output row has `group`, `condition`, `kind`, `pos` (offset in bp, or the number of core bases changed), `alt`, `rep`, the usual readouts and `atac_target` (mean ATAC over the 11 bp at the offset for duplicate/control/move). `native_core` and `core_found` record the site's core sequence and whether an exact TGA[C/G]TCA was present (otherwise the central 7 bp are used).

## Fine-tuned 1 bp ATAC head (results_ft/)
`FT=1` before any `submit_all.sh` mode reads the same experiment with the fine-tuned head instead of the built-in ones: `config_ft.json`
(`config_ft_cm_only.json` for scan/titrate), outputs in `results_ft/` under the same file names, job names prefixed `ft_`.
`AG_FT_CHECKPOINT=<file.delta.pth>` overrides the head in the config (used to smoke-test with the pilot checkpoint).
Every trackset is read from one forward pass: `CM` / `FB` are the fine-tuned head averaged over the four conditions, with per-condition
columns (`atac_peak_YS`, ...), and `CM_builtin` / `FB_builtin` are the built-in heads for comparison. A head without a 1 bp output is
refused: its motif and 170 bp readouts would only be 128 bp bins.

| mode | what | output |
|---|---|---|
| `ft_eval` | native predictions at every CM and FB site, plus a check that reading the head on the shared forward pass equals the fine-tuned model's own pass | `results_ft/eval_*` |
| `rerun_all` | every experiment the analyses read (`rerun_sites` then `rerun_dynamics`) | `results_ft/*` |

## Fibroblast priorities (`FT=1 bash submit_all.sh fb_priority`)
1. **Shuffle arm** (`phase1_fb`): 5,018 FB sites x 31 passes.
2. **Tile scan** (`fb_scan`): `data/fb_subset_scan.csv` = all 1,180 Dynamic sites + 1,180 Stable drawn at random (not matched), 364 passes/site.
3. **Element swaps** (`fb_swap`, `fb_background`): `data/fb_swap_sites.csv` = 50 sites per binding class (Stable, Age Up, Age Down,
   Exercise Up, Exercise Down; each site drawn for one class only, rarest first). Every donor element goes into every recipient location,
   at 170 bp and whole-peak width; for the peak width the recipient's own peak is fully removed (bases the donor does not cover are
   dinucleotide-shuffled). Self-swaps are the control and must reproduce the native prediction. The same elements also go into quiet host
   windows (`Finetune: bash submit.sh backgrounds` -> closed in every heart cell type, GC-matched to the elements) with a shuffled copy as
   the control. `python scripts/18_fb_swap.py --selftest` checks every edit without a GPU; analysis: `python3 Scripts/47_fb_swap_analysis.py`.

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
- Cost: 72 forward passes per site with the default replicates (10 shuffles at each width; random, dinuc and flank 5), split into phase 1 (31 passes: native + shuffles) and phase 2 (42 passes: native + mutations). At 0.29 s/sequence on an L40S that is about 50 GPU-hours for all 8,697 sites (phase 1 about 22, phase 2 about 29). Wall-clock = GPU-hours / GPUs you get at once, so spread jobs over queues (e.g. phase1 under `ACCOUNT=mwac_wonglab`, phase2 on the general H200 queue). `perturb.reps_width` can cap replicates per width if ever needed (off by default).
- Peak swap: recipients are CM sites at the 10/40/70/95th percentile of native peak ATAC; donors are 50 Stable + 25 Dynamic. Motif insertion:
  50 random mm10 backgrounds, dinucleotide-shuffled, TGAGTCA inserted at the centre, shuffled-motif control.
- The slide site counts (2,663 CM / 2,320 Stable) are not reproduced by the TOBIAS table: these jobs use all 3,679 CM and 5,018 fibroblast Fos sites.
