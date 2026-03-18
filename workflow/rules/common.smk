import os

# ─── Pipeline switches ────────────────────────────────────────────────────────
QC_RAW     = config['Run']['QC_raw']
TRIM       = config['Run']['trim']
QC_TRIMMED = config['Run']['QC_trimmed']
ALIGN      = config['Run']['align']
QUANTIFY   = config['Run']['quantify']
DESEQ2     = config['Run']['deseq2']

# ─── Library properties ───────────────────────────────────────────────────────
PAIRED      = config['Library']['paired_end']
READ_LENGTH = config['Library']['read_length']
STRANDEDNESS = config['Library']['strandedness']

# ─── Paths ────────────────────────────────────────────────────────────────────
INPUT      = config['User']['input_dir']
GENOME     = config['Reference']['genome_fasta']
GTF        = config['Reference']['gtf']
STAR_INDEX = config['Reference']['star_index']

# ─── Samples and contrasts ────────────────────────────────────────────────────
SAMPLES   = list(config['samples'].keys())
CONTRASTS = [c[0] for c in config['DESeq2']['contrasts']]

# ─── Helper functions ─────────────────────────────────────────────────────────

def get_raw_fastq_r1(wildcards):
    """Return the R1 (or only) FASTQ path for a sample."""
    return config['samples'][wildcards.sample]['R1']

def get_raw_fastq_r2(wildcards):
    """Return the R2 FASTQ path for a paired-end sample."""
    return config['samples'][wildcards.sample]['R2']

def get_trimmed_r1(wildcards):
    """Return the trimmed R1 path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz"
    return config['samples'][wildcards.sample]['R1']

def get_trimmed_r2(wildcards):
    """Return the trimmed R2 path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz"
    return config['samples'][wildcards.sample]['R2']

def get_trimmed_se(wildcards):
    """Return the trimmed single-end path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_trimmed.fq.gz"
    return config['samples'][wildcards.sample]['R1']

def get_star_input(wildcards):
    """Return STAR input FASTQ(s) as a list."""
    if PAIRED:
        return [get_trimmed_r1(wildcards), get_trimmed_r2(wildcards)]
    return [get_trimmed_se(wildcards)]

def get_bam_files(_):
    """Return all sorted BAM files for featureCounts."""
    return expand("align/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
