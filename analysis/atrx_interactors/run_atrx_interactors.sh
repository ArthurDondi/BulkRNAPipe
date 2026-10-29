#!/usr/bin/env bash
# Run the standalone ATRX-interactor boxplots (see README.md); the patient
# cohort plots are in the report.
# Activate an R environment with DESeq2, ggplot2 and optparse first
# (e.g. the pipeline's deseq2 env built from workflow/envs/deseq2.yaml).
#
# Usage: bash analysis/atrx_interactors/run_atrx_interactors.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${HERE}/01_plot_epicode.R"
