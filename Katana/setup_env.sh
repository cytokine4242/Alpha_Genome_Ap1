#!/bin/bash
# One-off environment setup on Katana (run on a login node or interactive GPU job):  bash setup_env.sh
# Creates conda env "alphagenome". If the lab already has a working env for alphagenome-pytorch (e.g. Minh's),
# skip this and set ENV_NAME in pbs/common.sh to it instead.
set -e
source ~/.bashrc
conda create -y -n alphagenome python=3.11
conda activate alphagenome
pip install alphagenome-pytorch pyfaidx pandas numpy safetensors pyarrow
python -c "import torch, alphagenome_pytorch; print('torch', torch.__version__, '| alphagenome_pytorch ok')"
