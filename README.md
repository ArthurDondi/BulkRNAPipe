# BulkRNAPipe

A Snakemake workflow for gold-standard bulk RNA-seq analysis on SLURM clusters.

## Overview

BulkRNAPipe implements the standard best-practice bulk RNA-seq pipeline used in
projects such as ENCODE and GTEx.  The modular design mirrors the
[MultiomePipe](https://github.com/ArthurDondi/MultiomePipe) single-cell pipeline
architecture: one Snakefile entry point, separate rule files per analysis stage,
per-rule conda environments, and a dedicated SLURM profile.

### Pipeline steps

| Step | Tool | Notes |
|------|------|-------|
| Raw read QC | FastQC + MultiQC | Per-sample + aggregated HTML report |
| Adapter trimming | Trim Galore | Wraps Cutadapt; supports paired-end and single-end |
| Trimmed read QC | FastQC + MultiQC | Optional; recommended |
| Genome indexing | STAR 2.7 | Built once; reused across all samples |
| Alignment | STAR (2-pass mode) | Splice-aware; produces sorted BAM |
| BAM indexing | samtools index | Required by downstream tools |
| Alignment QC | MultiQC | Aggregates STAR `Log.final.out` files |
| Quantification | featureCounts (Subread) | Gene-level counts |
| Differential expression | DESeq2 | Volcano plot + MA plot + results table |

## Repository structure

```
BulkRNAPipe/
├── config/
│   └── config.yaml                 # Template experiment config
├── data/                           # Place raw FASTQ files here (or symlink)
├── profile/
│   └── slurm/
│       └── config.yaml             # Snakemake SLURM executor profile
├── workflow/
│   ├── Snakefile                   # Pipeline entry point
│   ├── rules/
│   │   ├── common.smk              # Shared variables and helper functions
│   │   ├── qc.smk                  # FastQC + MultiQC
│   │   ├── trim.smk                # Trim Galore
│   │   ├── align.smk               # STAR genome index + alignment
│   │   ├── quantify.smk            # featureCounts
│   │   └── deseq2.smk              # DESeq2
│   ├── envs/
│   │   ├── fastqc.yaml
│   │   ├── trimgalore.yaml
│   │   ├── star.yaml
│   │   ├── subread.yaml
│   │   └── deseq2.yaml
│   └── scripts/
│       └── deseq2.R                # DESeq2 differential expression script
├── run_BulkRNAPipe.sh              # Local (non-SLURM) run script
└── run_BulkRNAPipe_slurm.sh        # SLURM cluster run script
```

## Install

```bash
# Create a conda environment with Snakemake
conda create -n BulkRNAPipe -c conda-forge -c bioconda snakemake=8 python=3.12 -y
conda activate BulkRNAPipe

# Scope strict channel priority to this environment only (required)
conda config --env --set channel_priority strict

# Install the SLURM executor plugin (required for cluster runs only)
pip install snakemake-executor-plugin-slurm
```

### Pre-create all pipeline conda environments

You can install every pipeline environment before running any jobs.  This is
useful to catch environment-creation errors early without needing real input
files:

```bash
snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --use-conda \
    --conda-frontend conda \
    --conda-create-envs-only \
    --cores 1
```

A minimal `config/config.yaml` that satisfies parse-time config access is
sufficient – input files are not required for environment creation.


## Quick start

### 1. Configure your experiment

Copy and edit the template config:

```bash
cp config/config.yaml config/config_myproject.yaml
```

Key fields to update:

| Field | Description |
|-------|-------------|
| `User.input_dir` | Directory containing raw FASTQ files |
| `User.output_dir` | Directory where pipeline outputs are written |
| `Reference.genome_fasta` | Path to the genome FASTA (e.g. GRCh38.fa) |
| `Reference.gtf` | Path to the gene annotation GTF |
| `Reference.star_index` | Where to build / find the STAR index |
| `Library.paired_end` | `True` for paired-end, `False` for single-end |
| `Library.read_length` | Expected read length (e.g. 150) |
| `Library.strandedness` | 0 = unstranded, 1 = stranded, 2 = reverse-stranded |
| `samples` | Sample names, paths to R1/R2 FASTQs, and conditions |
| `DESeq2.contrasts` | Pairwise comparisons to run |
| `DESeq2.combine_conditions` | *(Optional)* Merge existing conditions into a new label for DESeq2 |

#### Combining conditions for DESeq2

You can pool samples from two or more existing conditions into a single new
condition label — without averaging counts — using `DESeq2.combine_conditions`.
This is useful when you want to use a merged control group as the reference for
multiple contrasts.

```yaml
DESeq2:
  combine_conditions:
    ATRX_VectorControl:   # new combined condition name (must not clash with existing ones)
      - ATRX_NoVector     # existing condition labels to merge
      - ATRX_EmptyVector

  contrasts:
    # Contrast against the combined reference
    - [ATRX_VectorControl_vs_ATRX_FL,  ATRX_FL,  ATRX_VectorControl]
    - [ATRX_VectorControl_vs_ATRX_IFF, ATRX_IFF, ATRX_VectorControl]
```

**How it works**: Snakemake replaces each sample's original condition label with
the combined name before passing the `sample:condition` mapping to DESeq2.
Individual samples (replicates) are **not merged or averaged** — they remain
separate columns in the count matrix and DESeq2 estimates dispersion from all
replicates as usual.

**Validation** (errors at pipeline startup):
- The combined name must not collide with any existing condition label.
- Every listed source condition must exist in at least one sample.
- Each source condition may appear in at most one combined group.

Omit `combine_conditions` (or set it to `{}`) to use conditions exactly as
defined in `samples`, which is the default behaviour.

### 2. Dry run

Always do a dry run first to check the rule graph:

```bash
snakemake \
    -s workflow/Snakefile \
    --configfile config/config_myproject.yaml \
    --cores 1 \
    --use-conda \
    --conda-frontend conda \
    -n -p
```

### 3a. Local run

```bash
bash run_BulkRNAPipe.sh
# or with your config:
snakemake -s workflow/Snakefile --configfile config/config_myproject.yaml \
    --cores 8 --use-conda --conda-frontend conda -p
```

### 3b. SLURM cluster run

Edit `profile/slurm/config.yaml` to match your cluster, then:

```bash
bash run_BulkRNAPipe_slurm.sh
```

Snakemake itself runs on the login node and dispatches each rule as a separate
Slurm job.  Run it inside a `screen` or `tmux` session to keep it alive after
disconnecting.

## Running on a SLURM cluster

### SLURM profile settings

Edit `profile/slurm/config.yaml`:

| Field | Description | Default |
|-------|-------------|---------|
| `slurm_partition` | Partition/queue for jobs | `cpu` |
| `mem_mb` | Default memory per job (MB) | `16000` |
| `runtime` | Default wall-clock limit (minutes) | `240` |
| `cpus_per_task` | Default CPUs per job | `4` |
| `jobs` | Max concurrent Slurm jobs | `50` |

Rules that are resource-intensive override these defaults via their own
`resources:` blocks (e.g. `STARindex` requests 64 GB RAM and 16 CPUs;
`STARalign` requests 32 GB and 8 CPUs).

## Outputs

All outputs are written inside `User.output_dir`:

```
output_dir/
├── QC/
│   ├── raw/multiqc/multiqc_report.html       # Raw-read QC report
│   ├── trimmed/multiqc/multiqc_report.html   # Trimmed-read QC report
│   └── align/multiqc/multiqc_report.html     # Alignment QC report
├── trim/{sample}/                            # Trimmed FASTQs
├── align/{sample}/                           # BAM files + STAR logs
├── quantify/counts.txt                       # Gene × sample count matrix
├── deseq2/{contrast}/
│   ├── results.csv                           # DE results table
│   ├── normalized_counts.csv                 # DESeq2-normalized counts
│   ├── volcano.pdf                           # Volcano plot
│   └── ma_plot.pdf                           # MA plot
├── logs/                                     # Per-rule log files
└── benchmark/                                # Per-rule benchmark files
```
