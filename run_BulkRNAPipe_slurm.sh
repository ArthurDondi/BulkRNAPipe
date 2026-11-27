#!/bin/bash
#SBATCH --output logs/BulkRNAPipe_%j.log
#SBATCH --error  logs/BulkRNAPipe_%j.err
#SBATCH --job-name=BulkRNAPipe
#SBATCH --partition=longq      # 30d limit: long enough for the whole workflow
#SBATCH --qos=longq            # qos must match the partition
#SBATCH --time=2-00:00:00      # 2 days (raise for very large sample sets)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # this job is only the Snakemake controller
#SBATCH --mem=8000             # controller is light; rules get their own jobs
#SBATCH --mail-type=end
#SBATCH --mail-user=arthur.dondi@cemm.at
#
# Run the whole BulkRNAPipe Snakemake workflow on the SLURM cluster, then
# re-render the two analysis reports (cell lines, patients).
#
# Three Snakemake calls, in order: (0) a dry run that logs which jobs will
# run; (1) the default target, i.e. every pipeline output; (2) the reports,
# forced. The Report rule only declares
# design_qc/design_qc_summary.csv as input, so in a single call it could
# render before the GSEA / SharedGenes outputs exist, and an existing report
# would count as up to date even after its Rmd changed.
#
# This batch job is the Snakemake *controller*: it stays alive for the whole
# run and submits each rule as its own SLURM job through the profile in
# profile/slurm/config.yaml (per-rule memory / walltime / partition). It
# therefore needs few resources itself but a long walltime, hence longq.
#
# One-time setup (in the BulkRNAPipe conda env):
#   pip install snakemake-executor-plugin-slurm
#
# Create the controller log directory once before the first submission:
#   mkdir -p logs
#
# Submit with (from the BulkRNAPipe/ root):
#   sbatch run_BulkRNAPipe_slurm.sh [config/your_config.yaml]
#
# With no argument it uses config/config_epicode.yaml, matching
# run_BulkRNAPipe.sh. config/config.yaml is the documented template and still
# holds /path/to/... placeholders, so it is never the default.
#
# Edit the #SBATCH log paths / --mail-user above for your account.
#
# Rerunning after a config change (e.g. new contrasts or DESeq2.sample_sets):
# only the missing or outdated outputs are built; existing ones are kept unless
# their inputs, params or conda env changed (rerun-triggers in the profile).
# If a previous controller was killed, unlock the output directory first:
#   snakemake -s workflow/Snakefile --configfile config/config_epicode.yaml --unlock

set -euo pipefail

CONFIG="${1:-config/config_epicode.yaml}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config file not found: $CONFIG" >&2
    exit 1
fi

# config/config.yaml is the template and still carries /path/to/... placeholders.
# Without this check Snakemake accepts them and fails inside pathlib trying to
# mkdir '/path', which does not point at the actual problem.
if grep -qE '^[^#]*/path/to/' "$CONFIG"; then
    echo "ERROR: '$CONFIG' still contains /path/to/... placeholders:" >&2
    grep -nE '^[^#]*/path/to/' "$CONFIG" | sed 's/^/         /' >&2
    echo >&2
    echo "       config/config.yaml is the template, not a runnable config." >&2
    echo "       Pass your experiment config instead, e.g.:" >&2
    echo "         sbatch run_BulkRNAPipe_slurm.sh config/config_epicode.yaml" >&2
    exit 1
fi

echo "======================"
echo "submit dir : ${SLURM_SUBMIT_DIR:-$PWD}"
echo "job name   : ${SLURM_JOB_NAME:-BulkRNAPipe}"
echo "partition  : ${SLURM_JOB_PARTITION:-longq}"
echo "job id     : ${SLURM_JOB_ID:-NA}"
echo "config     : $CONFIG"
echo "======================"

# Rerun triggers (mtime, params, software-env) come from the profile.
# 0. What will run, per rule, into this job's log (nothing is submitted).
snakemake \
    -s workflow/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profile/slurm \
    --rerun-incomplete \
    --dry-run --quiet rules

# 1. The whole workflow: the default target, BulkRNAPipe_output.
#    --rerun-incomplete: rebuild outputs left half-written by a killed job.
#    --keep-going: one failing job does not stop the independent ones; the
#    exit status is still non-zero, so step 2 is skipped (set -e) and the
#    reports are never rendered from an incomplete run.
snakemake \
    -s workflow/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profile/slurm \
    --rerun-incomplete \
    --keep-going

# 2. The reports, always re-rendered, now that every output above exists.
#    --forcerun Report re-runs both (the rule renders one report per target).
#    set -e above stops the script here if step 1 failed.
snakemake \
    -s workflow/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profile/slurm \
    --forcerun Report \
    -p report/cell_lines_report.html report/patients_report.html
