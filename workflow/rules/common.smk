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

# Validate: every source condition must actually exist in the sample list.
for _combined_name, _source_list in _combine_conditions.items():
    for _src in _source_list:
        if _src not in _existing_conditions:
            raise ValueError(
                f"combine_conditions: source condition '{_src}' (in combined group "
                f"'{_combined_name}') does not match any sample's condition. "
                "Check your config for typos."
            )

# Validate: for each contrast, check two kinds of conflicts:
#   (a) A source condition belongs to two combined groups that are both used in
#       this contrast (ambiguous remapping between groups).
#   (b) One side of the contrast is a raw condition that is also a source member
#       of the combined group on the other side (e.g. numerator=condA,
#       denominator=groupB where condA ∈ groupB).  Those samples would be
#       simultaneously assigned to both sides of the contrast.
# Across different contrasts the same source condition may appear freely.
for _c in config['DESeq2']['contrasts']:
    _cname, _num, _den = _c[0], _c[1], _c[2]
    # (a) conflict between two combined groups
    _contrast_src_to_group = {}
    for _combined_name, _source_list in _combine_conditions.items():
        if _combined_name not in (_num, _den):
            continue
        for _src in _source_list:
            if _src in _contrast_src_to_group:
                raise ValueError(
                    f"Contrast '{_cname}': source condition '{_src}' maps to both "
                    f"'{_contrast_src_to_group[_src]}' and '{_combined_name}', "
                    "which are both referenced by this contrast. A source condition "
                    "cannot appear in two combined groups used in the same contrast."
                )
            _contrast_src_to_group[_src] = _combined_name
    # (b) raw condition on one side absorbed by combined group on the other side
    for _combined_name, _source_list in _combine_conditions.items():
        # combined group is the denominator; raw condition is the numerator
        if _combined_name == _den and _num in _source_list:
            raise ValueError(
                f"Contrast '{_cname}': numerator '{_num}' is a raw condition that "
                f"is also a source member of the denominator group '{_den}'. "
                "Samples with that condition would be counted on both sides."
            )
        # combined group is the numerator; raw condition is the denominator
        if _combined_name == _num and _den in _source_list:
            raise ValueError(
                f"Contrast '{_cname}': denominator '{_den}' is a raw condition that "
                f"is also a source member of the numerator group '{_num}'. "
                "Samples with that condition would be counted on both sides."
            )

def get_contrast_effective_condition(sample, contrast_entry):
    """Return the DESeq2 condition for *sample* for a specific contrast.

    Only combine_conditions groups that are referenced by the given contrast
    (as numerator or denominator) are applied.  This allows the same source
    condition to participate in different combined groups across separate
    contrasts without conflict.

    Parameters
    ----------
    sample : str
        Sample name (key in config['samples']).
    contrast_entry : list
        Three-element list [name, numerator, denominator] for this contrast.
    """
    orig = config['samples'][sample]['condition']
    num, den = contrast_entry[1], contrast_entry[2]
    for combined_name, source_list in _combine_conditions.items():
        if combined_name in (num, den) and orig in source_list:
            return combined_name
    return orig

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
