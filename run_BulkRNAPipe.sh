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
#
# ── Conda prefix on NFS/CIFS/isilon filesystems ─────────────────────────────
# If conda environment creation fails with:
#   [Errno 2] No such file or directory: '.../.snakemake/conda/.../man/man3/App::Cpan.3'
# your working directory is on a filesystem that disallows ':' in file names
# (e.g. Isilon/OneFS, CIFS/SMB mounts).  Move the conda env cache to a local
# scratch directory by setting CONDA_ENVS_PREFIX before calling this script:
#
#   export CONDA_ENVS_PREFIX=/tmp/snakemake_conda   # or any local path
#   bash run_BulkRNAPipe.sh
#
# The pipeline output files remain in the directory specified by output_dir in
# your config; only the reusable conda environments are placed under the prefix.

# Optional override: set CONDA_ENVS_PREFIX to redirect conda env storage to a
# local filesystem (avoids NFS/CIFS colon-in-filename errors).
CONDA_ENVS_PREFIX="${CONDA_ENVS_PREFIX:-}"
CONDA_PREFIX_ARG=()
if [[ -n "$CONDA_ENVS_PREFIX" ]]; then
    mkdir -p "$CONDA_ENVS_PREFIX"
    CONDA_PREFIX_ARG=(--conda-prefix "$CONDA_ENVS_PREFIX")
fi

snakemake \
    -s workflow/Snakefile \
    --configfile config/config_epicode.yaml \
    --cores 8 \
    --use-conda \
    --conda-frontend conda \
    "${CONDA_PREFIX_ARG[@]}" \
    --resources mem_mb=64000 \
    --rerun-triggers mtime params \
    -p -n \
    "$@"
