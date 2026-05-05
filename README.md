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

## Troubleshooting

### Conda environment creation fails: `No such file or directory: .../man/man3/App::Cpan.3`

**Symptom:**

```
CreateCondaEnvironmentException:
...
ERROR conda.core.link:_execute(1033): An error occurred while installing package
'conda-forge::perl-5.32.1-7_hd590300_perl5'.
[Errno 2] No such file or directory: '.../.snakemake/conda/.../man/man3/App::Cpan.3'
```

**Root cause:**  
`conda-forge`'s `perl-5.32.1` build ships man-page files whose names contain
colons (`App::Cpan.3`, `Module::Build.3`, etc.).  Filesystems that disallow
`:` in file names — including **Isilon/OneFS**, **CIFS/SMB mounts**, and some
**NFS** configurations — return `ENOENT` when conda tries to create those
files, causing the entire environment installation to roll back.

> **Note:** This is a known package bug in specific conda-forge perl builds.
> `workflow/envs/fastqc.yaml` pins `perl=5.32.1=pl5321*` (bioconda build) which
> does **not** install man pages and works on all tested filesystems.

**Fix A — move the conda env cache to a local filesystem (recommended for NFS):**

```bash
export CONDA_ENVS_PREFIX=/tmp/snakemake_conda   # any path on a local disk
# or use $TMPDIR if it points to local scratch on your HPC:
# export CONDA_ENVS_PREFIX=$TMPDIR/snakemake_conda

bash run_BulkRNAPipe.sh --conda-create-envs-only
```

Both `run_BulkRNAPipe.sh` and `run_BulkRNAPipe_slurm.sh` honour the
`CONDA_ENVS_PREFIX` environment variable and pass it as `--conda-prefix` to
Snakemake.  Set it in your `.bashrc` / `.bash_profile` to make the change
permanent.

**Fix B — delete the broken partial environment and retry:**

```bash
# Remove the broken .snakemake/conda directory from inside your project folder
rm -rf .snakemake/conda
bash run_BulkRNAPipe.sh --conda-create-envs-only
```

**Fix C — enable strict channel priority (conda-forge recommendation):**

```bash
conda config --env --set channel_priority strict
```

With strict priority the solver picks one consistent set of builds per
channel, reducing the chance of picking a problematic conda-forge perl build.

---

### "Your conda installation is not configured to use strict channel priorities"

Run once per environment:

```bash
conda config --env --set channel_priority strict
```

This warning is printed by Snakemake but does not abort the run.  It is a
best-practice recommendation from the conda-forge team.
