#!/bin/bash
# Run BulkRNAPipe locally (foreground, interactive).
#
# --cores 8                 : up to 8 CPU cores in total
# --resources mem_mb=64000  : declare 64 GB of RAM available across all jobs
#
# Edit --configfile to point to your experiment config before running.
# Run from the BulkRNAPipe/ root directory:
#   bash run_BulkRNAPipe.sh

snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 8 \
    --use-conda \
    --conda-frontend conda \
    --resources mem_mb=64000 \
    --rerun-triggers mtime params \
    -p
