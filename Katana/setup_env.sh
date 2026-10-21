#!/bin/bash
# One-off environment setup on Katana (run on a login node or interactive job):  bash setup_env.sh
# Creates conda env "alphagenome2" with Python 3.12 (required by the current alphagenome-pytorch) and installs the
# GitHub version (the pip release 0.2.8 is too old: no from_pretrained / named_outputs).
# To match Minh's exact code, set REF to the commit of /srv/scratch/wonglab/minh/alphagenome/alphagenome-pytorch:
#   git -c safe.directory='*' -C /srv/scratch/wonglab/minh/alphagenome/alphagenome-pytorch rev-parse HEAD
#   REF=<hash> bash setup_env.sh
set -e
source ~/.bashrc
ENV=alphagenome2
conda create -y -n $ENV python=3.12
conda activate $ENV
pip install pyfaidx pandas numpy safetensors pyarrow
pip install "git+https://github.com/genomicsxai/alphagenome-pytorch${REF:+@$REF}"
python -c "import torch, alphagenome_pytorch as a; from alphagenome_pytorch.named_outputs import TrackMetadataCatalog; from alphagenome_pytorch import AlphaGenome; print('torch', torch.__version__, '| from_pretrained:', hasattr(AlphaGenome,'from_pretrained'), '| named_outputs ok')"
