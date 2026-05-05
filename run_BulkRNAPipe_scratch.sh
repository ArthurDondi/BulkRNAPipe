#!/bin/bash
# Run BulkRNAPipe with conda environments cached on scratch.
#
# Snakemake will create (or reuse) all pipeline conda environments under a
# shared prefix on scratch, so environments persist across runs and do not need
# to be re-solved or re-installed each time.
#
# Usage (from the BulkRNAPipe/ root directory):
#   bash run_BulkRNAPipe_scratch.sh [SNAKEMAKE_ARGS...]
#
# Override the conda-prefix directory (default shown below):
#   SMK_CONDA_PREFIX=/path/to/envs bash run_BulkRNAPipe_scratch.sh
#
# Override the config file (default: config/config_epicode.yaml):
#   CONFIGFILE=config/config_myproject.yaml bash run_BulkRNAPipe_scratch.sh
#
# Pre-create all conda environments without running any jobs:
#   bash run_BulkRNAPipe_scratch.sh --conda-create-envs-only
#
# --cores 8                 : up to 8 CPU cores in total
# --resources mem_mb=64000  : declare 64 GB of RAM available across all jobs
#
# Edit --configfile (or set CONFIGFILE below) to point to your experiment
# config before running.
#
# Tip: install mamba or micromamba for faster environment resolution:
#   conda install -n base -c conda-forge mamba
# Then change --conda-frontend below from "conda" to "mamba".

# Default path is user-specific; other users should override SMK_CONDA_PREFIX.
SMK_CONDA_PREFIX="${SMK_CONDA_PREFIX:-/scratch/research/taschnermandl/arthur_d/smk_conda}"
CONFIGFILE="${CONFIGFILE:-config/config_epicode.yaml}"

mkdir -p "${SMK_CONDA_PREFIX}"

snakemake \
    -s workflow/Snakefile \
    --configfile "${CONFIGFILE}" \
    --cores 8 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "${SMK_CONDA_PREFIX}" \
    --resources mem_mb=64000 \
    --rerun-triggers mtime params \
    -p \
    "$@"
