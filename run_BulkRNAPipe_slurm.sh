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
# Run the whole BulkRNAPipe Snakemake workflow on the SLURM cluster.
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
# Edit the #SBATCH log paths / --mail-user above for your account.

set -euo pipefail

CONFIG="${1:-config/config.yaml}"

echo "======================"
echo "submit dir : ${SLURM_SUBMIT_DIR:-$PWD}"
echo "job name   : ${SLURM_JOB_NAME:-BulkRNAPipe}"
echo "partition  : ${SLURM_JOB_PARTITION:-longq}"
echo "job id     : ${SLURM_JOB_ID:-NA}"
echo "config     : $CONFIG"
echo "======================"

snakemake \
    -s workflow/Snakefile \
    --configfile "$CONFIG" \
    --workflow-profile profile/slurm \
    --rerun-triggers mtime params software-env \
    -p
