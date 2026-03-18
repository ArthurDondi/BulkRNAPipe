#!/bin/bash
# Run BulkRNAPipe on a SLURM cluster.
#
# Each rule is submitted as an individual Slurm job. Default resources
# (partition, memory, runtime, CPUs) are defined in profile/slurm/config.yaml.
# Edit that file to match your cluster's configuration before running.
#
# Requirements:
#   pip install snakemake-executor-plugin-slurm
#
# Usage (from the BulkRNAPipe/ root directory):
#   bash run_BulkRNAPipe_slurm.sh
#
# Slurm logs are written to logs/slurm/ inside the output directory.
#
# --configfile             : change to point to your own config file
# --profile profile/slurm  : load the Slurm executor profile
# --rerun-triggers ...     : resubmit when inputs or params change

snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --profile profile/slurm \
    --rerun-triggers mtime params \
    -p
