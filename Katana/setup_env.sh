#!/bin/bash
# One-off environment setup on Katana (login node or interactive job):  bash setup_env.sh
# Creates a Python 3.12 conda env ON SCRATCH (home has a small quota; torch + CUDA is several GB) and installs the GitHub
# alphagenome-pytorch (the pip release 0.2.8 is too old: no from_pretrained / named_outputs). The env path is saved to
# .env_prefix, which pbs/common.sh reads.
#   ENV_PREFIX=/srv/scratch/wonglab/jack/envs/alphagenome2 bash setup_env.sh     # choose another location
#   REF=<commit> bash setup_env.sh                                              # pin the package to Minh's commit:
#   git -c safe.directory='*' -C /srv/scratch/wonglab/minh/alphagenome/alphagenome-pytorch rev-parse HEAD
set -e
source ~/.bashrc
ENV_PREFIX=${ENV_PREFIX:-/srv/scratch/wonglab/jack/envs/alphagenome2}
CACHE_ROOT=$(dirname "$ENV_PREFIX")/cache
mkdir -p "$(dirname "$ENV_PREFIX")" "$CACHE_ROOT/pip" "$CACHE_ROOT/tmp" "$CACHE_ROOT/conda_pkgs"
export PIP_CACHE_DIR="$CACHE_ROOT/pip" TMPDIR="$CACHE_ROOT/tmp" CONDA_PKGS_DIRS="$CACHE_ROOT/conda_pkgs"
echo "env: $ENV_PREFIX | caches: $CACHE_ROOT"
conda create -y -p "$ENV_PREFIX" python=3.12
conda activate "$ENV_PREFIX"
pip install pyfaidx pandas numpy safetensors pyarrow
pip install "git+https://github.com/genomicsxai/alphagenome-pytorch${REF:+@$REF}"
python -c "import torch, alphagenome_pytorch as a; from alphagenome_pytorch.named_outputs import TrackMetadataCatalog; from alphagenome_pytorch import AlphaGenome; print('torch', torch.__version__, '| from_pretrained:', hasattr(AlphaGenome,'from_pretrained'), '| named_outputs ok')"
echo "$ENV_PREFIX" > "$(dirname "$0")/.env_prefix"
echo "done; PBS jobs will use $ENV_PREFIX"
