# Fine-tuning AlphaGenome on mouse cardiomyocyte and fibroblast ATAC (all four conditions)

Drag the `Finetune/` folder onto Katana (next to `Katana/`, inside `/srv/scratch/wonglab/jack/AlphaGenome_Ap1/Alpha_Genome_Ap1`). All data, folds and run outputs go under `/srv/scratch/wonglab/jack/AlphaGenome_Ap1/finetune_data`; nothing is written outside `jack`. It builds on the `alphagenome-pytorch` fine-tuning framework that Minh already uses (`scripts/finetune.py`, his venv), so the training loop itself is not new code. What is new here: the 1 bp target bigWigs, mouse folds, the generated configs, and the submit scripts.

## Recommendation: frozen trunk + new heads first, then LoRA; not full fine-tuning
- **Start with `probe` (frozen trunk, new heads, one output track per cell type and condition).**
  - The perturbation results (shuffles, ISM, motif edits) depend on the trunk's learned sequence grammar. A frozen trunk leaves it exactly as it is, so a new head can only improve the readout and the effects stay comparable with everything already run.
  - It is cheap, and the built-in heads leave a lot of headroom: the heart head tracks the cardiomyocyte bigWig at rho 0.59 and the Patski head tracks the fibroblast bigWig at 0.28 over our peaks.
- **Then `lora` (rank 8 on q/v, optionally `lora_locon`) only if the probe falls short.** Adapters change the trunk a little, and are saved as small delta checkpoints. Use it if the probe's differential accessibility between conditions is poor (below).
- **Do not start with `full`.** One genome and 8 to 10 tracks is little for 450 M parameters; the trunk can drift and overfit to the training chromosomes, and the perturbation effect sizes could shift in ways that are hard to check. A `full` config is generated only as a last comparison.
- **Decision rule after the probe and LoRA run on the same data and splits:** keep the probe unless LoRA improves held-out Pearson (128 bp) by more than 0.03 or the held-out differential correlation (MAS minus YS at peaks) by more than 0.05. Only if LoRA saturates is `full` worth a try, with early stopping.

## Which target files
| files | resolution | what they are | use |
|---|---|---|---|
| Gadi `pooled/bams/*.bam.bw` (Mar 2024) | 50 bp bins | raw read counts, not depth-normalised (max about 8.7e5) | no |
| `*_merged.SeqDepthNorm.bw` (Gadi `RemainPooled/bam/bigwig`, 48 files; Katana has YS only) | 10 bp bins | `bamCoverage --binSize 10 --normalizeUsing CPM`: read coverage, comparable across conditions | **pilot, 128 bp resolution only** |
| **new `<Cell>_<Cond>.insertions.bw`** from the fragment files | **1 bp** | Tn5 insertion sites (start+4, end-5), explicit zeros, scaled to a fixed total | **final fine-tune, resolutions "1,128"** |

The 10 bp read-coverage tracks blur the base-pair footprints that the motif-level analyses use, so the final fine-tune should use 1 bp insertion tracks built from the raw fragments. The 10 bp tracks are good enough to test the whole pipeline now.

## Steps
0. **Get the targets onto Katana (not done).** For the pilot, copy the 8 cardiomyocyte and fibroblast SeqDepthNorm files from Gadi (`/g/data/zk16/jack/HeartSCRNA-Jun22/analysis/2024-02-07.ChromBp/RemainPooled/bam/bigwig/`, about 5 GB) into `$FT` = `/srv/scratch/wonglab/jack/AlphaGenome_Ap1/finetune_data` (everything on Katana stays inside `jack`), for example `rsync -av gadi:<dir>/{Cardiomyocyte,Fibroblast}_{YS,MAS,MAE,YE}_merged.SeqDepthNorm.bw $FT/bigwigs_10bp/` (if Katana cannot reach Gadi, pull to your laptop and push). Also needed on Katana: `mm10.fa` (`/srv/scratch/wonglab/rubj/beagle/rawData/mm10.fa` exists) and `mm10.chrom.sizes` (`/srv/scratch/wonglab/useful/genomes/mm10.chrom.sizes`).
1. **Mouse folds:** `bash scripts/make_mouse_folds.sh $FT/folds_mouse` (Borzoi mouse folds to AlphaGenome 1 Mb regions; FOLD_3 pairs with `model_fold_3.safetensors`, valid = Borzoi fold6, test = fold7; also writes `train_pilot.bed`).
2. **Track metadata (optional):** `python scripts/make_track_metadata.py --out $FT/track_metadata.tsv [--with-rna]`.
3. **Configs:** `python scripts/make_configs.py --kind 10bp --bigwig-dir $FT/bigwigs_10bp --genome mm10.fa --folds-dir $FT/folds_mouse --work-dir $FT/runs --metadata $FT/track_metadata.tsv` (add `--pilot` for 1 epoch on 200 windows, `--seq-length 131072` for faster training closer to the perturbation runs, `--with-rna --rna-dir ...` to add the merged RNA tracks).
4. **Pilot, then probe, then LoRA:** `PILOT=1 bash submit.sh probe`, then `bash submit.sh probe`, then `bash submit.sh lora` (`KIND=1bp` for the insertion tracks; `GPU=H200`, `BATCH=1 ACCUM=32` for a smaller card; `DRY=1` prints the qsub command). Read the pilot's step time before scaling; Minh's 15-epoch human probe at 1 Mb asked for 45 h on an H200.
5. **Final targets (1 bp):** fill `data/samples_manifest.tsv` (one row per sample: its fragments file and the barcodes of that cell type in it), then `MANIFEST=data/samples_manifest.tsv CHROM_SIZES=mm10.chrom.sizes OUT_DIR=$FT/bigwigs_1bp bash submit.sh bigwigs` (an array, one task per cell type and condition), `python scripts/make_configs.py --kind 1bp ...` and repeat step 4 with `KIND=1bp`. Each sample goes in separately so barcodes cannot collide between libraries.
6. **Evaluate on held-out chromosomes** (FOLD_3 test.bed), against the built-in heart and Patski heads, with `scripts/evaluate_finetuned.py` from the framework (`--checkpoint`, `--pretrained-weights`, `--test-bed`, `--metrics`; check `--help`, not tested here). Judge by:
   - held-out Pearson at 128 bp and 1 bp, and Spearman over peaks, per track;
   - differential accessibility (MAS minus YS, MAE minus MAS) against the bigWig differences at held-out peaks;
   - cell-type specificity (cardiomyocyte head on cardiomyocyte tracks vs fibroblast tracks);
   - the Fos-site perturbation results rerun with the new head (needs a small extension of `Katana/scripts/ag_lib.py` to load a delta checkpoint; build once the first checkpoint exists, so its format is known);
   - the fibroblast MPRA as an independent check on the fibroblast head.

## What was tested here, and what was not
- **Tested locally:** `make_insertion_bigwigs.py` on synthetic fragments (insertion positions, per-sample barcode filtering, chrM dropped, no missing values, scaling); the config generator and its YAML (valid, 8 ATAC tracks, correct modes and resolutions); the manifest runner and the `submit.sh` commands (dry run); the fold converter and the YAML keys against the cluster code (`organism` and `track_metadata` are accepted; `save_delta` is refused with `full`, as the configs respect).
- **Not tested (needs the GPU):** the fine-tuning runs themselves, the mouse folds script, the evaluation, the memory fit of 1 Mb windows on the smaller cards in the private queue, and bigWig reading speed (convert to mmap with `agt preprocess bigwig-to-mmap` if reading is the bottleneck).

## Caveats
- The only RNA coverage is merged over conditions, so RNA tracks cannot be condition-specific.
- Condition differences are subtle (peaks agree across conditions at rho 0.8 to 0.9), so held-out differential correlation is the key number to watch; with merged replicates there is no noise ceiling. If per-replicate bigWigs are built (one manifest row per sample, a separate track each), a replicate-based ceiling becomes possible.
- The backbone `model_fold_3` has not seen the FOLD_3 mouse test windows; `model_all_folds` (used for the perturbation runs) has, so evaluate fine-tuned heads with a fold backbone and only afterwards consider using the all-folds trunk in the perturbation pipeline.
