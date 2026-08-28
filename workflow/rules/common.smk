import os

# ─── Pipeline switches ────────────────────────────────────────────────────────
QC_RAW     = config['Run']['QC_raw']
TRIM       = config['Run']['trim']
QC_TRIMMED = config['Run']['QC_trimmed']
ALIGN      = config['Run']['align']
QUANTIFY   = config['Run']['quantify']
DESEQ2     = config['Run']['deseq2']
GSEA       = config['Run'].get('gsea', False)

# ─── Wildcard constraints ────────────────────────────────────────────────────
# Every wildcard below names one path segment, never a path. Snakemake's default
# `.+` matches "/" too, which lets a nested output be read as an outer rule's
# wildcard: gsea/{contrast}/kegg/{slug}_dotplot.pdf also matches the GSEA rule
# with contrast="{contrast}/kegg", and the two rules come out ambiguous. Pinning
# them to a single segment removes the ambiguity at the source, rather than
# papering over it with ruleorder.
wildcard_constraints:
    sample     = r"[^/]+",
    contrast   = r"[^/]+",
    collection = r"[^/]+",
    interaction = r"[^/]+",
    comparison = r"[^/]+",
    sampleset  = r"[^/]+",
    view       = r"[^/]+",
    slug       = r"[^/]+",
    ont        = r"[^/]+",
    dir        = r"[^/]+",
GOENRICH   = config['Run'].get('go', False)
PROTEOMICS = config.get('Proteomics', {}).get('enabled', False)

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
        return f"trim/{wildcards.sample}/{wildcards.sample}_val_1.fq.gz"
    return os.path.join(INPUT, config['samples'][wildcards.sample]['R1'])

def get_trimmed_r2(wildcards):
    """Return the trimmed R2 path, or raw if trimming is skipped."""
    if TRIM:
        return f"trim/{wildcards.sample}/{wildcards.sample}_val_2.fq.gz"
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

# ─── PCA views, sample subsets and gene exclusion ────────────────────────────
# When the reference contains extra contigs for the delivery vector (EGFP,
# mCherry, a transgene ORF, ...), those features are present in counts.txt and
# can drive the sample-level PCA on their own: they are zero in untransduced
# samples and very highly expressed in transduced ones.  On top of that, the
# transgene itself (e.g. ATRX and its vector-borne variants) is by construction
# the largest difference between the groups, so a PCA that includes it mostly
# re-plots the experimental design.
#
# The PCA module therefore runs a grid of (sample set x gene view):
#   views       all_genes | no_vector_genes | no_vector_no_transgene
#   sample sets all + every entry under PCA.subsets
_pca_cfg = config.get('PCA', {}) or {}

PCA_EXCLUDE_GENES         = [str(g) for g in (_pca_cfg.get('exclude_genes') or [])]
PCA_EXCLUDE_GENE_PATTERNS = [str(p) for p in (_pca_cfg.get('exclude_gene_patterns') or [])]
PCA_EXCLUDE_CONTIGS       = [str(c) for c in (_pca_cfg.get('exclude_contigs') or [])]
PCA_TRANSGENE_GENES       = [str(g) for g in (_pca_cfg.get('transgene_genes') or [])]

PCA_VIEWS = ["all_genes", "no_vector_genes", "no_vector_no_transgene"]

# Commas are the delimiter used to pass these lists to the R scripts.
for _lst, _name in (
    (PCA_EXCLUDE_GENES, 'exclude_genes'),
    (PCA_EXCLUDE_GENE_PATTERNS, 'exclude_gene_patterns'),
    (PCA_EXCLUDE_CONTIGS, 'exclude_contigs'),
    (PCA_TRANSGENE_GENES, 'transgene_genes'),
):
    for _entry in _lst:
        if ',' in _entry:
            raise ValueError(
                f"PCA.{_name}: entry '{_entry}' contains a comma, which is used "
                "as the list delimiter. Split it into separate list entries."
            )

# ── Sample subsets ────────────────────────────────────────────────────────────
# Restricting the PCA to a subset of conditions is the way to stop one dominant
# axis (e.g. a different parental clone) from compressing every other effect
# into the higher PCs.
_pca_subsets = _pca_cfg.get('subsets') or {}

for _name, _conds in _pca_subsets.items():
    if _name == 'all':
        raise ValueError("PCA.subsets: 'all' is reserved for the full sample set.")
    if not _conds:
        raise ValueError(f"PCA.subsets['{_name}']: no conditions listed.")
    for _c in _conds:
        if _c not in _existing_conditions:
            raise ValueError(
                f"PCA.subsets['{_name}']: condition '{_c}' does not match any "
                "sample's condition. Check your config for typos."
            )
    _n = sum(1 for s in SAMPLES if config['samples'][s]['condition'] in _conds)
    if _n < 3:
        raise ValueError(
            f"PCA.subsets['{_name}']: only {_n} sample(s) selected. A PCA needs "
            "at least 3."
        )

PCA_SUBSETS     = {str(k): [str(c) for c in v] for k, v in _pca_subsets.items()}
PCA_SAMPLE_SETS = ['all'] + sorted(PCA_SUBSETS)

def get_pca_samples(sampleset):
    """Return the "sample:condition" string for a PCA sample set."""
    if sampleset == 'all':
        selected = SAMPLES
    else:
        conds    = PCA_SUBSETS[sampleset]
        selected = [s for s in SAMPLES if config['samples'][s]['condition'] in conds]
    return ",".join(f"{s}:{config['samples'][s]['condition']}" for s in selected)

# ── Visualisation-only batch removal ──────────────────────────────────────────
# Maps each condition to a batch label.  Used ONLY by limma::removeBatchEffect
# for a supplementary PCA panel; it never enters a DESeq2 design, because in a
# design where batch is a deterministic function of condition the batch effect
# is not identifiable (see docs/design_confounding_plan.md).
_pca_batch = _pca_cfg.get('batch') or {}
for _cond in _pca_batch:
    if _cond not in _existing_conditions:
        raise ValueError(
            f"PCA.batch: condition '{_cond}' does not match any sample's "
            "condition. Check your config for typos."
        )
PCA_BATCH = {str(k): str(v) for k, v in _pca_batch.items()}

def get_pca_batch(sampleset):
    """Return the "sample:batch" string for a PCA sample set (empty when unset)."""
    if not PCA_BATCH:
        return ""
    if sampleset == 'all':
        selected = SAMPLES
    else:
        conds    = PCA_SUBSETS[sampleset]
        selected = [s for s in SAMPLES if config['samples'][s]['condition'] in conds]
    return ",".join(
        f"{s}:{PCA_BATCH.get(config['samples'][s]['condition'], 'unassigned')}"
        for s in selected
    )

# ─── Design QC (library-level confounder check) ───────────────────────────────
# Genes whose per-sample expression is worth plotting explicitly: the reporters,
# the transgene and its variants.  Defaults to the vector features plus the
# transgene list so the module is useful without extra configuration.
_qc_cfg = config.get('DesignQC', {}) or {}
DESIGN_QC_GENES = [str(g) for g in (_qc_cfg.get('genes') or
                                    (PCA_EXCLUDE_GENES + PCA_TRANSGENE_GENES))]

# Endogenous genes to plot for biological QC.  Unlike the list above these are
# NOT excluded from size-factor estimation - they are ordinary genes, not
# constructs that are structurally absent from some groups.
DESIGN_QC_GOI = [str(g) for g in (_qc_cfg.get('genes_of_interest') or [])]

for _entry in DESIGN_QC_GENES + DESIGN_QC_GOI:
    if ',' in _entry:
        raise ValueError(
            f"DesignQC: entry '{_entry}' contains a comma, which is used as the "
            "list delimiter. Split it into separate list entries."
        )

# ─── Design QC: pooled and tag-only transgene features ───────────────────────
# The features above are counted the way the main matrix counts everything:
# unique alignments only.  When the reference carries cDNA copies of a gene that
# also exists in the genome (a transgene and its endogenous locus), most of that
# gene's reads align equally well to two or three places, are flagged NH >= 2 by
# STAR, and are dropped.  Neither the endogenous feature nor the construct
# features then report the gene's real output.
#
# These two optional features recover it from the SAME BAMs, no re-alignment:
#
#   total_feature  one meta-feature spanning several gene_ids.  Counted with
#                  -M --fraction, every alignment of an ambiguous fragment lands
#                  inside the same meta-feature and contributes 1/NH, so the
#                  fragment is counted exactly once no matter which copy it came
#                  from.  This is the gene's total output across all sources.
#   tag            one meta-feature over a coordinate range shared by the vector
#                  contigs (the purification tag).  Same -M --fraction trick, so
#                  a fragment that maps to every construct carrying the tag is
#                  still counted once.  Sequence absent from the genome, so this
#                  measures transgene output regardless of which construct a
#                  sample carries.
# Display-only renaming for the design QC panels.  The count matrix, the CSVs
# and every other output keep the real feature IDs - this maps an ID to the name
# shown on the facet strip, nothing more.
_qc_labels = _qc_cfg.get('plot_labels') or {}
DESIGN_QC_LABELS = ",".join(f"{k}={v}" for k, v in _qc_labels.items())

for _k, _v in _qc_labels.items():
    if any(c in f"{_k}{_v}" for c in (',', '=')):
        raise ValueError(
            f"DesignQC.plot_labels: '{_k}: {_v}' must not contain ',' or '=' - "
            "both are used as delimiters."
        )

_qc_total = _qc_cfg.get('total_feature') or {}
DESIGN_QC_TOTAL_NAME    = str(_qc_total.get('name', 'ATRX_Total'))
DESIGN_QC_TOTAL_GENES   = [str(g) for g in (_qc_total.get('genes') or [])]
DESIGN_QC_SUBSTITUTE_AS = str(_qc_total.get('substitute_as', '') or '')

if DESIGN_QC_SUBSTITUTE_AS and not DESIGN_QC_TOTAL_GENES:
    raise ValueError(
        "DesignQC.total_feature.substitute_as is set but "
        "DesignQC.total_feature.genes is empty - nothing to pool."
    )
if DESIGN_QC_SUBSTITUTE_AS and (',' in DESIGN_QC_SUBSTITUTE_AS or '"' in DESIGN_QC_SUBSTITUTE_AS):
    raise ValueError(
        f"DesignQC.total_feature.substitute_as '{DESIGN_QC_SUBSTITUTE_AS}' must "
        'not contain a comma or a double quote.'
    )

# quantify/counts.txt is NEVER modified. When substitute_as is set, DESeq2,
# DESeq2Interaction and PCA read this derived matrix instead - genes are
# dropped and total_feature's pooled counts stand in for the target gene.
# Everything downstream of DESeq2 (GSEA, GO, signature reversal, the
# proteomics overlay) inherits the substitution automatically, since none of
# those rules read counts.txt directly.
DOWNSTREAM_COUNTS = (
    "design_qc/counts_substituted.txt" if DESIGN_QC_SUBSTITUTE_AS
    else "quantify/counts.txt"
)

_qc_tag = _qc_cfg.get('tag') or {}
DESIGN_QC_TAG_NAME    = str(_qc_tag.get('name', 'Tag'))
DESIGN_QC_TAG_REGIONS = list(_qc_tag.get('regions') or [])
DESIGN_QC_TAG_MIN_OV  = int(_qc_tag.get('min_overlap', 25))

for _name in (DESIGN_QC_TOTAL_NAME, DESIGN_QC_TAG_NAME):
    if ',' in _name or '"' in _name:
        raise ValueError(
            f"DesignQC: feature name '{_name}' must not contain a comma or a "
            'double quote.'
        )

for _r in DESIGN_QC_TAG_REGIONS:
    for _key in ('contig', 'start', 'end'):
        if _key not in _r:
            raise ValueError(
                f"DesignQC.tag.regions: entry {_r} is missing '{_key}'. Each "
                "region needs contig, start and end (1-based, inclusive)."
            )
    if int(_r['start']) < 1 or int(_r['end']) < int(_r['start']):
        raise ValueError(
            f"DesignQC.tag.regions: entry {_r} has an empty or negative range. "
            "Coordinates are 1-based and inclusive, so start >= 1 and end >= start."
        )

# The tag range is a coordinate window, not a sequence match: it is only correct
# if the contigs really do begin at the construct's first base. Check it against
# the FASTA before trusting the numbers.
#
# Passed to the rule as one flat string rather than a ready-made SAF file, so no
# tab or newline ever has to survive a trip through the shell; awk expands it.
DESIGN_QC_TAG_SPEC = ";".join(
    "{},{},{},{},{}".format(
        DESIGN_QC_TAG_NAME, _r['contig'], int(_r['start']), int(_r['end']),
        _r.get('strand', '+'),
    )
    for _r in DESIGN_QC_TAG_REGIONS
)

for _r in DESIGN_QC_TAG_REGIONS:
    if ';' in str(_r['contig']) or ',' in str(_r['contig']):
        raise ValueError(
            f"DesignQC.tag.regions: contig '{_r['contig']}' must not contain "
            "',' or ';' - both are used as delimiters."
        )

# ─── Interaction contrasts (difference of differences) ───────────────────────
# Each entry estimates (A1 - A2) - (B1 - B2) on the log2 scale from the plain
# ~ condition fit.  This is the correct way to compare across a boundary whose
# effect cannot be estimated as a covariate: the nuisance effect cancels in the
# subtraction instead of being modelled.
_interactions = config['DESeq2'].get('interactions') or []

for _it in _interactions:
    for _key in ('name', 'group_A', 'group_B'):
        if _key not in _it:
            raise ValueError(f"DESeq2.interactions: entry is missing '{_key}'.")
    if len(_it['group_A']) != 2 or len(_it['group_B']) != 2:
        raise ValueError(
            f"DESeq2.interactions['{_it['name']}']: group_A and group_B must each "
            "be exactly two condition names [numerator, denominator]."
        )
    for _c in list(_it['group_A']) + list(_it['group_B']):
        if _c not in _existing_conditions:
            raise ValueError(
                f"DESeq2.interactions['{_it['name']}']: condition '{_c}' does not "
                "match any sample's condition. Check your config for typos."
            )
    if len(set(list(_it['group_A']) + list(_it['group_B']))) < 4:
        raise ValueError(
            f"DESeq2.interactions['{_it['name']}']: the four conditions must be "
            "distinct; a repeated condition makes the difference of differences "
            "collapse to a simple contrast."
        )

INTERACTIONS = [_it['name'] for _it in _interactions]

if set(INTERACTIONS) & set(CONTRASTS):
    raise ValueError(
        "DESeq2.interactions: name(s) "
        f"{sorted(set(INTERACTIONS) & set(CONTRASTS))} collide with a contrast "
        "name. Interaction results would overwrite the contrast's output."
    )

def get_interaction_cfg(name):
    """Return the DESeq2.interactions config block for *name*."""
    for it in _interactions:
        if it['name'] == name:
            return it
    raise ValueError(f"No DESeq2.interactions entry named '{name}'")

# ─── Signature-reversal comparisons ──────────────────────────────────────────
# Scores a rescue by asking whether it reverses a knockout signature, using two
# contrasts that each live entirely inside one genetic background.  Immune to
# the between-background effect, which never enters either contrast.
_sig_reversal = config.get('SignatureReversal') or []

for _sr in _sig_reversal:
    for _key in ('name', 'signature_contrast', 'response_contrast'):
        if _key not in _sr:
            raise ValueError(f"SignatureReversal: entry is missing '{_key}'.")
    for _key in ('signature_contrast', 'response_contrast'):
        if _sr[_key] not in CONTRASTS:
            raise ValueError(
                f"SignatureReversal['{_sr['name']}']: {_key} '{_sr[_key]}' is not "
                "a DESeq2 contrast name. Add it to DESeq2.contrasts first."
            )
    if _sr['signature_contrast'] == _sr['response_contrast']:
        raise ValueError(
            f"SignatureReversal['{_sr['name']}']: signature_contrast and "
            "response_contrast are the same contrast."
        )

SIGNATURE_REVERSALS = [_sr['name'] for _sr in _sig_reversal]

def get_signature_reversal_cfg(name):
    """Return the SignatureReversal config block for *name*."""
    for sr in _sig_reversal:
        if sr['name'] == name:
            return sr
    raise ValueError(f"No SignatureReversal entry named '{name}'")

def get_contrast_sides(contrast_name):
    """Return (numerator, denominator) for a DESeq2 contrast name."""
    for c in config['DESeq2']['contrasts']:
        if c[0] == contrast_name:
            return c[1], c[2]
    raise ValueError(f"No DESeq2 contrast named '{contrast_name}'")

# ─── GSEA / GO derived variables ─────────────────────────────────────────────

# Flatten collection identifiers to safe filesystem names, e.g.
# "C2:CP:REACTOME" → "C2_CP_REACTOME"
def _collection_to_slug(col):
    return col.replace(":", "_")

GSEA_COLLECTIONS = [
    _collection_to_slug(c)
    for c in config.get('GSEA', {}).get('collections', [])
] if GSEA else []

# ─── KEGG BRITE category panels ──────────────────────────────────────────────
# The collection dotplot shows the top 15 pathways in each direction, which is
# the wrong shape for a question like "what happened to replication and repair?":
# a category is a small fixed list, and the pathways that did NOT move are part
# of the answer. These panels slice the KEGG results by BRITE category and draw
# every pathway in it.
#
# Only meaningful when C2:CP:KEGG is among GSEA.collections, since they re-plot
# that collection's results rather than running fgsea again.
_gsea_cfg = config.get('GSEA', {}) or {}
KEGG_CATEGORIES_FILE = _gsea_cfg.get(
    'kegg_categories_file',
    f"{workflow.basedir}/resources/kegg_categories.tsv",
)
KEGG_CATEGORIES = [str(c) for c in (_gsea_cfg.get('kegg_categories') or [])]
# Which collection the panels slice. MSigDB's C2:CP:KEGG is a frozen ~2011
# snapshot of KEGG, so every pathway added since is missing from it; KEGG_CURRENT
# is the GMT GenerateKeggGmt fetches from the live KEGG API instead.
KEGG_SOURCE_SLUG = _collection_to_slug(
    str(_gsea_cfg.get('kegg_source_collection', 'C2:CP:KEGG')))
KEGG_PANELS = bool(GSEA and KEGG_CATEGORIES and KEGG_SOURCE_SLUG in GSEA_COLLECTIONS)

if KEGG_CATEGORIES and GSEA and KEGG_SOURCE_SLUG not in GSEA_COLLECTIONS:
    raise ValueError(
        f"GSEA.kegg_categories is set but its source collection "
        f"'{KEGG_SOURCE_SLUG}' is not in GSEA.collections. The panels slice that "
        "collection's results, so add it, point kegg_source_collection at one "
        "you do run, or clear kegg_categories."
    )

for _slug in KEGG_CATEGORIES:
    if '/' in _slug or _slug != _slug.strip():
        raise ValueError(
            f"GSEA.kegg_categories: '{_slug}' must be a bare slug matching the "
            "first column of the category file - it becomes a filename."
        )

# ─── GO keyword category panels ──────────────────────────────────────────────
# The same five groupings applied to GO. GO has no BRITE-style categories, so
# these are REGEX MATCHES ON TERM NAMES, not an ontology traversal - they catch
# terms that share a word and miss terms that say it differently. The panels and
# the log both state this; do not read one as "everything GO knows about X".
GO_CATEGORIES_FILE = _gsea_cfg.get(
    'go_categories_file',
    f"{workflow.basedir}/resources/go_categories.tsv",
)
GO_CATEGORIES  = [str(c) for c in (_gsea_cfg.get('go_categories') or [])]
GO_SOURCE_SLUG = _collection_to_slug(
    str(_gsea_cfg.get('go_source_collection', 'C5:GO:BP')))
GO_CAT_TOP_N   = int(_gsea_cfg.get('go_categories_top_n', 25))
GO_PANELS = bool(GSEA and GO_CATEGORIES and GO_SOURCE_SLUG in GSEA_COLLECTIONS)

if GO_CATEGORIES and GSEA and GO_SOURCE_SLUG not in GSEA_COLLECTIONS:
    raise ValueError(
        f"GSEA.go_categories is set but its source collection "
        f"'{GO_SOURCE_SLUG}' is not in GSEA.collections. The panels slice that "
        "collection's results, so add it, point go_source_collection at one you "
        "do run, or clear go_categories."
    )

for _slug in GO_CATEGORIES:
    if '/' in _slug or _slug != _slug.strip():
        raise ValueError(
            f"GSEA.go_categories: '{_slug}' must be a bare slug matching the "
            "first column of the GO category file - it becomes a filename."
        )


# Contrast comparison pairs
_contrast_comparisons = config.get('ContrastComparisons') or []
CONTRAST_COMPARISONS  = [cc['name'] for cc in _contrast_comparisons]

def get_comparison_cfg(name):
    """Return the ContrastComparisons config block for *name*."""
    for cc in _contrast_comparisons:
        if cc['name'] == name:
            return cc
    raise ValueError(f"No ContrastComparisons entry named '{name}'")

# GO ontologies (BP / MF / CC)
GO_ONTOLOGIES = config.get('GO', {}).get('ontology', ['BP']) if GOENRICH else []

# ─── Proteomics limma overlay config ──────────────────────────────────────────
_proteomics_cfg = config.get('Proteomics', {}) or {}

def get_proteomics_comparison(contrast_name):
    """Return the mapped proteomics comparison string for a DESeq2 contrast."""
    mapping = _proteomics_cfg.get('deseq2_to_proteomics_comparison') or {}
    return mapping.get(contrast_name, "")

# Subset of DESeq2 contrasts that have a non-empty proteomics comparison mapping
# in deseq2_to_proteomics_comparison.  Only these contrasts get a
# DESeq2 × proteomics overlay (volcano_proteomic.pdf); contrasts without a
# mapping are skipped so the overlay rule never runs with an empty comparison
# (which would otherwise fail with "Missing required proteomics parameters").
PROTEOMICS_CONTRASTS = [c for c in CONTRASTS if get_proteomics_comparison(c)]
