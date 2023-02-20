#! /usr/bin/env bash

eval "$(conda shell.bash hook)"
conda create -n mhcnuggets python=3.8
conda activate mhcnuggets

cd src/mhcnuggets
pip install -e ./

echo "Downloading GRCh references (build 108) for pyensembl/varcode..."
pyensembl install --release 108
