#!/bin/bash
# Run BulkRNAPipe locally (foreground, interactive).
#
# --cores 8                 : up to 8 CPU cores in total
# --resources mem_mb=64000  : declare 64 GB of RAM available across all jobs
#
# Edit --configfile to point to your experiment config before running.
# Run from the BulkRNAPipe/ root directory:
#   bash run_BulkRNAPipe.sh
#
# Pre-create all conda environments (recommended before a first run):
#   bash run_BulkRNAPipe.sh --conda-create-envs-only
#
# Tip: install mamba or micromamba for faster environment resolution:
#   conda install -n base -c conda-forge mamba
# Then change --conda-frontend below from "conda" to "mamba".

snakemake \
    -s workflow/Snakefile \
    --configfile config/config_epicode.yaml \
    --cores 8 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /scratch/research/taschnermandl/arthur_d/smk_conda \
    --resources mem_mb=64000 \
    --rerun-triggers mtime params \
    -p  \
    "$@"
