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

mkdir -p logs

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
        -p"
