#!/bin/bash
# Run BulkRNAPipe on a SLURM cluster.
#
# This script submits the Snakemake orchestrator itself as a SLURM job on
# mediumq (up to 2 days), which in turn dispatches each pipeline rule as its
# own individual Slurm job via the profile in profile/slurm/config.yaml.
#
# Requirements:
#   pip install snakemake-executor-plugin-slurm
#
# Usage (from the BulkRNAPipe/ root directory):
#   bash run_BulkRNAPipe_slurm.sh
#
# The Snakemake log is written to logs/snakemake_<jobid>.out/err in the
# current directory.
#
# --configfile             : change to point to your own config file
# --profile profile/slurm  : load the Slurm executor profile
# --rerun-triggers ...     : resubmit when inputs or params change
#
# ── Conda prefix on NFS/CIFS/isilon filesystems ─────────────────────────────
# If conda environment creation fails with:
#   [Errno 2] No such file or directory: '.../.snakemake/conda/.../man/man3/App::Cpan.3'
# your working directory is on a filesystem that disallows ':' in file names.
# Move the conda env cache to a local scratch directory by setting
# CONDA_ENVS_PREFIX before calling this script:
#
#   export CONDA_ENVS_PREFIX=/tmp/snakemake_conda   # or $TMPDIR or local scratch
#   bash run_BulkRNAPipe_slurm.sh

mkdir -p logs

# Optional override: set CONDA_ENVS_PREFIX to redirect conda env storage to a
# local filesystem (avoids NFS/CIFS colon-in-filename errors).
CONDA_ENVS_PREFIX="${CONDA_ENVS_PREFIX:-}"
CONDA_PREFIX_ARG=""
if [[ -n "$CONDA_ENVS_PREFIX" ]]; then
    mkdir -p "$CONDA_ENVS_PREFIX"
    # Use printf %q so paths with spaces are properly escaped inside --wrap
    CONDA_PREFIX_ARG="--conda-prefix $(printf '%q' "$CONDA_ENVS_PREFIX")"
fi

sbatch \
    --partition=mediumq \
    --qos=mediumq \
    --job-name=BulkRNAPipe \
    --output=logs/snakemake_%j.out \
    --error=logs/snakemake_%j.err \
    --time=2-00:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=4G \
    --wrap="snakemake \
        -s workflow/Snakefile \
        --configfile config/config.yaml \
        --profile profile/slurm \
        --rerun-triggers mtime params \
        $CONDA_PREFIX_ARG \
        -p"
