#!/usr/bin/env bash
# Run the whole ATRX-interactor boxplot analysis (see README.md).
# Activate an R environment with DESeq2, ggplot2, optparse and readxl first
# (e.g. the pipeline's deseq2 env built from workflow/envs/deseq2.yaml).
#
# Usage: bash analysis/atrx_interactors/run_atrx_interactors.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Clinical sheet for 03_* (read in place, never copied; only patient_o_id,
# atrx and mna are used).
CLINICAL="/nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret/metadata/20230524_TGFß_Fikrets_RNAseq.xlsx"

bash    "${HERE}/00_download_GSE94035.sh"
Rscript "${HERE}/01_plot_epicode.R"
Rscript "${HERE}/02_plot_GSE94035.R"
Rscript "${HERE}/03_plot_GSE94035_by_ATRX_MYCN.R" --clinical "${CLINICAL}"
