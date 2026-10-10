#!/usr/bin/env bash
set -euo pipefail

echo "=================================================="
echo " Setting up 'pd-microbiome' Conda Environment"
echo "=================================================="

CONDA_CMD="mamba"
if ! command -v $CONDA_CMD &>/dev/null; then
  CONDA_CMD="conda"
fi

echo "Using package manager: $CONDA_CMD"

export CONDA_SUBDIR=osx-64

echo "Creating environment 'pd-microbiome' from environment.yml..."
$CONDA_CMD env create -f environment.yml -n pd-microbiome -y

echo "=================================================="
echo " Environment created successfully!"
echo " To activate: conda activate pd-microbiome"
echo "=================================================="
