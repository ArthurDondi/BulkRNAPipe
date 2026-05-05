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

# ─── Combined / derived conditions ───────────────────────────────────────────
# Read the optional combine_conditions mapping from config.
_combine_conditions = config['DESeq2'].get('combine_conditions') or {}

# Collect all condition labels present in the sample list.
_existing_conditions = {config['samples'][s]['condition'] for s in SAMPLES}

# Validate: combined name must not collide with any existing condition label.
for _combined_name in _combine_conditions:
    if _combined_name in _existing_conditions:
        raise ValueError(
            f"combine_conditions: combined name '{_combined_name}' collides with "
            "an existing sample condition label. Choose a different name."
        )

# Validate: each source condition may appear in at most one combined group,
# and must actually exist in the sample list.
_seen_sources = {}
for _combined_name, _source_list in _combine_conditions.items():
    for _src in _source_list:
        if _src in _seen_sources:
            raise ValueError(
                f"combine_conditions: condition '{_src}' is listed in both "
                f"'{_seen_sources[_src]}' and '{_combined_name}'. "
                "Each source condition may only appear in one combined group."
            )
        _seen_sources[_src] = _combined_name
        if _src not in _existing_conditions:
            raise ValueError(
                f"combine_conditions: source condition '{_src}' (in combined group "
                f"'{_combined_name}') does not match any sample's condition. "
                "Check your config for typos."
            )

# Build reverse map: original_condition → effective_condition
_condition_remap = {
    src: combined
    for combined, srcs in _combine_conditions.items()
    for src in srcs
}

def get_effective_condition(sample):
    """Return the DESeq2 condition for *sample*, remapped via combine_conditions if set."""
    orig = config['samples'][sample]['condition']
    return _condition_remap.get(orig, orig)

# ─── Helper functions ─────────────────────────────────────────────────────────

def get_raw_fastq_r1(wildcards):
    """Return the R1 (or only) FASTQ path for a sample."""
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R1'])

def get_raw_fastq_r2(wildcards):
    """Return the R2 FASTQ path for a paired-end sample."""
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R2'])

def get_trimmed_r1(wildcards):
    """Return the trimmed R1 path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz"
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R1'])

def get_trimmed_r2(wildcards):
    """Return the trimmed R2 path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz"
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R2'])

def get_trimmed_se(wildcards):
    """Return the trimmed single-end path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_trimmed.fq.gz"
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R1'])

def get_star_input(wildcards):
    """Return STAR input FASTQ(s) as a list."""
    if PAIRED:
        return [get_trimmed_r1(wildcards), get_trimmed_r2(wildcards)]
    return [get_trimmed_se(wildcards)]

def get_bam_files(_):
    """Return all sorted BAM files for featureCounts."""
    return expand("align/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES)
