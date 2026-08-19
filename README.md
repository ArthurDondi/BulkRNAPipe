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
| GSEA | fgsea (fgseaMultilevel) | Per contrast; msigdbr + HOX + custom GMTs |
| GO enrichment | clusterProfiler enrichGO | Per contrast; gene symbols auto-mapped to Entrez |
| Contrast comparison | ΔNES + residual-rank fgsea | For configured contrast pairs |

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
│   │   ├── deseq2.smk              # DESeq2
│   │   ├── gsea.smk                # GSEA (fgsea)
│   │   ├── go.smk                  # GO enrichment (clusterProfiler)
│   │   └── compare_contrasts.smk   # ΔNES + residual-rank contrast comparisons
│   ├── envs/
│   │   ├── fastqc.yaml
│   │   ├── trimgalore.yaml
│   │   ├── star.yaml
│   │   ├── subread.yaml
│   │   ├── deseq2.yaml
│   │   └── gsea.yaml               # fgsea + msigdbr + clusterProfiler
│   └── scripts/
│       ├── deseq2.R                # DESeq2 differential expression
│       ├── pca.R                   # PCA grid (sample sets x gene views)
│       ├── design_qc.R              # Library-level confounder QC
│       ├── deseq2_interaction.R     # Difference-of-differences contrasts
│       ├── signature_reversal.R     # Rescue scored as signature reversal
│       ├── generate_hox_gmt.R      # Auto-generates HOX gene-set GMT
│       ├── gsea.R                  # fgsea enrichment per contrast
│       ├── go.R                    # GO enrichment per contrast
│       └── gsea_compare.R          # ΔNES + residual-rank comparison
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
| `Proteomics` | *(Optional)* Generate a separate proteomics-filtered DESeq2 volcano plot per contrast |

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
- A source condition may appear in multiple combined groups, as long as those groups are not both used as numerator/denominator in the same contrast (that would make the remapping ambiguous for that contrast).

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

One-time setup in the BulkRNAPipe conda env — install the SLURM executor plugin
and create the controller log directory:

```bash
pip install snakemake-executor-plugin-slurm
mkdir -p logs
```

Edit `profile/slurm/config.yaml` to match your cluster's partitions/QoS, then
submit the controller batch job (it runs Snakemake on `longq` and dispatches
each rule as its own SLURM job):

```bash
sbatch run_BulkRNAPipe_slurm.sh                          # uses config/config_epicode.yaml
sbatch run_BulkRNAPipe_slurm.sh config/config_myproject.yaml
```

Edit the `#SBATCH` log paths / `--mail-user` at the top of that script for your
account. To run interactively instead (Snakemake on the login node, inside a
`screen`/`tmux` session so it survives disconnecting), dry-run first then drop
`-n`:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --workflow-profile profile/slurm -n
```

## Running on a SLURM cluster

### SLURM profile settings

Edit `profile/slurm/config.yaml`:

| Field | Description | Default |
|-------|-------------|---------|
| `slurm_partition` | Default partition/queue for jobs | `tinyq` (≤ 2 h) |
| `qos` | QoS — must match the partition | `tinyq` |
| `mem_mb` | Default memory per job (MB) | `16000` |
| `runtime` | Default wall-clock limit (minutes) | `120` |
| `cpus_per_task` | Default CPUs per job | `4` |
| `jobs` | Max concurrent Slurm jobs | `50` |

On this cluster the QoS must match the partition, and QoS is set via the
dedicated `qos` resource rather than `slurm_extra`, because the SLURM executor
plugin manages QoS itself and rejects `--qos` passed through `slurm_extra`.
`runtime` is expressed in **minutes**.

Rules that exceed the 2 h `tinyq` limit override the partition/QoS in their own
`resources:` block to route to `mediumq` (e.g. `STARindex` requests 64 GB RAM,
16 CPUs, `mediumq`; `STARalign` requests 32 GB, 8 CPUs, `mediumq`). Per-rule
`mem_mb` / `runtime` declared on the rules take precedence over the profile
defaults.

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
├── pca/{sample_set}/
│   ├── {view}_pc12.pdf                       # PC1 vs PC2
│   ├── {view}_pc34.pdf                       # PC3 vs PC4
│   ├── {view}_scree.pdf                      # Variance explained per PC
│   ├── {view}_pc12_batch_removed.pdf         # removeBatchEffect panel (visualisation only)
│   ├── {view}_coords.csv                     # PC coordinates
│   └── {view}_excluded_genes.csv             # Features dropped by this view
├── design_qc/
│   ├── library_size.pdf                      # Assigned reads per sample
│   ├── assignment_rates.pdf                  # featureCounts assignment breakdown
│   ├── genes_detected.pdf                    # Genes with >= 1 count
│   ├── marker_expression.pdf                 # Reporter / transgene expression
│   ├── sample_correlation.pdf                # Spearman correlation, all samples
│   ├── sample_distance.pdf                   # Euclidean distance, all samples
│   └── design_qc_summary.csv                 # Per-sample summary table
├── deseq2_interaction/{interaction}/
│   ├── results.csv                           # Difference-of-differences results
│   ├── components.csv                        # Interaction next to both component effects
│   ├── volcano.pdf                           # Volcano plot
│   └── ma_plot.pdf                           # MA plot
├── signature_reversal/{comparison}/
│   ├── reversal_scatter.pdf                  # Response log2FC vs signature log2FC
│   ├── recovery_distribution.pdf             # Per-gene recovered fraction
│   ├── reversal_stats.csv                    # Reversal score, rho, % reversed
│   └── reversal_genes.csv                    # Per-gene table
├── deseq2/{contrast}/
│   ├── results.csv                           # DE results table
│   ├── normalized_counts.csv                 # DESeq2-normalized counts
│   ├── volcano.pdf                           # Full RNA volcano plot
│   ├── volcano_proteomic.pdf                 # Optional proteomics-concordance volcano (only when Proteomics.enabled=True)
│   └── ma_plot.pdf                           # MA plot (direction-annotated)
├── resources/generated_gmts/
│   └── hox.gmt                              # Auto-generated HOX gene sets
├── gsea/{contrast}/
│   ├── {collection}_results.csv             # fgsea results per collection
│   └── {collection}_dotplot.pdf             # Dotplot of top enriched pathways
├── go/{contrast}/
│   ├── go_{ont}_{dir}_results.csv           # Raw GO enrichment (up/down)
│   ├── go_{ont}_{dir}_results_simplified.csv# Redundancy-reduced GO results
│   ├── go_{ont}_{dir}_dotplot.pdf           # Raw GO dotplot
│   ├── go_{ont}_{dir}_dotplot_simplified.pdf# Simplified GO dotplot
│   ├── go_{ont}_{dir}_unmapped_universe.csv # Input IDs not mapped in universe
│   └── go_{ont}_{dir}_unmapped_sig.csv      # Input IDs not mapped in sig set
├── gsea_compare/{comparison}/
│   ├── delta_nes_{collection}.csv           # Per-collection ΔNES table
│   ├── delta_nes_summary.csv                # ΔNES summary across all collections
│   ├── delta_nes_barplot.pdf                # ΔNES bar plot
│   └── residual_rank/
│       └── residual_rank_{collection}.csv   # Residual-rank fgsea per collection
├── logs/                                     # Per-rule log files
└── benchmark/                                # Per-rule benchmark files
```

## Confounded designs: PCA views, interactions and signature reversal

These four modules exist for experiments where the groups differ by more than
the intended variable — a different parental clone, an extra transduction step,
a selection round. In that situation the nuisance effect is a deterministic
function of the group label, which means it **cannot be estimated as a
covariate**: the model matrix becomes rank-deficient and DESeq2 refuses to fit
it. There is no group whose mean would pin the effect down, and replicates do
not help — they sharpen the group means without adding a new one.

`docs/design_confounding_plan.md` works through the arithmetic and explains why
RUV, ComBat and friends do not rescue this. The short version: the modules below
avoid the confounding rather than correcting it.

### PCA views and sample subsets

The `PCA` rule runs a grid of (sample set × gene view) into
`pca/{sample_set}/{view}_*`. For each combination you get PC1/PC2, PC3/PC4, a
scree plot, the coordinates as CSV and a supplementary batch-removal panel.

Gene views (all three always run):

| View | Genes used |
| --- | --- |
| `all_genes` | every gene in `counts.txt` |
| `no_vector_genes` | minus `exclude_genes` / `exclude_gene_patterns` / `exclude_contigs` |
| `no_vector_no_transgene` | minus those, and minus `transgene_genes` |

Features are dropped from the raw counts **before** the `DESeqDataSet` is built,
so size factors and the VST are re-estimated on the retained genes. This matters:
a reporter that is absent in one group and highly expressed in another shifts the
library-size normalisation of every other gene, not just its own row.

The third view exists because the transgene is, by construction, the largest
single difference between the groups. A PCA that includes it partly re-plots the
experimental design instead of showing its downstream consequences.

```yaml
PCA:
  exclude_genes: [EGFP, mCherry]          # exact IDs from the Geneid column
  exclude_gene_patterns: []               # regexes, for shared naming prefixes
  exclude_contigs: []                     # genes whose exons all sit here
  transgene_genes: [ATRX, ATRX_FL, ATRX_IFF]
  ntop: 500

  subsets:                                # extra sample sets beyond "all"
    E6_derived: [E6, EmptyVector, ATRX_FL, ATRX_IFF]

  batch:                                  # visualisation-only, see below
    TP53:        untransduced
    EmptyVector: transduced
```

Subsets are how you stop one dominant axis from compressing everything else into
the higher PCs: drop the group that sits on its own clonal branch and the
remaining effects get to use PC1/PC2.

Entries matching nothing are not an error, but `pca.R` prints a `WARNING` to
`logs/PCA/*.log` — check there if two views come out identical unexpectedly.
List the IDs your GTF actually uses with:

```bash
grep -o 'gene_id "[^"]*"' /path/to/genome_plus_vectors.gtf | sort -u
```

**The batch-removal panel is for looking, not for evidence.** When `PCA.batch`
varies within at least one condition, `limma::removeBatchEffect` runs with the
condition means protected. When every condition maps to exactly one batch, the
batch effect is not estimable at all — limma returns `NA` coefficients if the
design is supplied — so the script subtracts the batch means outright, which
produces the tidy figure people expect by deleting the real between-group
difference along with any technical one. In that case a red warning is printed
**onto the figure**, because a PDF outlives its log file.

### Interaction contrasts (difference of differences)

`(A1 − A2) − (B1 − B2)`, estimated from the plain `~ condition` fit. Since it is
a linear combination of group means, it is estimable where a free-standing
nuisance term is not, and any effect shared by both pairs cancels.

```yaml
DESeq2:
  interactions:
    - name: ATRX_FL_rescue_shortfall
      group_A: [ATRX_FL, EmptyVector]     # how far the rescue moved the gene
      group_B: [TP53, E6]                 # how far it had to move (restoration target)
```

`log2FC ≈ 0` means both comparisons move the gene by the same amount — for a
rescue experiment, that the rescue put the gene back where the knockout took it
from. Negative means it fell short, positive that it overshot. `components.csv`
reports the interaction next to both component effects so a large value can be
traced to whichever pair moved.

**Mind the direction of `group_B`.** For a rescue you want the restoration
target (`[wildtype, mutant]`), not the perturbation (`[mutant, wildtype]`).
Flipping it reports the sum of the two effects rather than the residual, and a
complete rescue scores twice the effect size instead of zero.

Note also what this is *not* needed for. When both sides of a comparison share a
nuisance factor it cancels on its own: `ATRX_FL` vs `EmptyVector` is already free
of the vector effect, because subtracting a separately-measured vector effect is
algebraically identical to using `EmptyVector` as the denominator —
`(FL − E6) − (EmptyVector − E6) = FL − EmptyVector`. Interactions are for
crossing a boundary that no shared control spans.

This assumes the nuisance effect is the same size in both pairs. With one group
per condition that is untestable, so state it in the methods.

### Signature reversal

Scores an intervention by how far it reverses a perturbation signature, using
two contrasts that each sit entirely inside one background — so the boundary is
never crossed and there is nothing to correct.

```yaml
SignatureReversal:
  - name: ATRX_FL_reverses_KO_signature
    signature_contrast: TP53_vs_E6              # defines which genes are scored
    response_contrast:  EmptyVector_vs_ATRX_FL  # measured on those genes
    padj_threshold: 0.05
    lfc_threshold: 0.0
    restored_fraction: 0.5
```

`reversal_stats.csv` reports:

| Field | Meaning |
| --- | --- |
| `reversal_score` | −slope of the total-least-squares fit; 1 = complete reversal, 0 = no response |
| `pct_reversed` | % of signature genes moving the opposite way |
| `null_pct_reversed` / `excess_over_null` / `perm_p` | permutation baseline, the excess over it, and its empirical p-value |
| `shares_condition` / `interpretable` | flags a pair whose contrasts share a group, which invalidates the score |
| `pct_restored` | % recovering at least `restored_fraction` of the perturbation |
| `spearman_rho` | rank correlation between the two fold-change vectors |

Total least squares rather than ordinary least squares: both axes are estimated
with comparable noise, and OLS would bias the slope toward zero.

**The two contrasts must not share a condition.** If they do, that group's
sampling noise enters one log2FC positively and the other negatively, creating a
negative correlation out of nothing — on pure noise, sharing one group of three
replicates gives ~66% "reversed" and rho ~ -0.5. The script detects this, sets
`interpretable = FALSE`, and stamps the warning onto the figure. It rules out the
obvious-looking negative control (a knockout signature against
`control vs control+empty_vector`), which shares the knockout group.

**Read `excess_over_null`, not `pct_reversed`.** Every comparison builds a
permutation null by shuffling the response log2FCs across genes, which preserves
both sign distributions and gives the baseline this statistic actually has on
this data. The binomial test against 50% is also reported, but it assumes
balanced signs in both contrasts and is optimistic when they are skewed.

**Beware a signature that is too broad.** Selecting on a contrast whose groups
differ by more than the intended variable admits genes no intervention could
reverse, which dilutes the slope toward zero. Raising `lfc_threshold` restricts
the signature to genes with a substantial effect.

### Design QC

`design_qc/` plots library size, featureCounts assignment rates, genes detected,
reporter/transgene expression, and the sample–sample correlation and distance
structure, with all samples on one figure — a group-linked technical shift is
invisible one sample at a time. Size factors for the marker panel are estimated
with the reporter/transgene features excluded, so a construct present in half the
groups cannot shift the normalisation it is being measured against.

Read it alongside the MultiQC reports in `QC/`, which carry duplication, adapter
and coverage-uniformity metrics that are not recomputed here.

## GSEA / GO enrichment modules

### Enabling GSEA and GO

Set the following toggles in your config YAML:

```yaml
Run:
  gsea: True   # fgsea enrichment per contrast
  go:   True   # GO enrichment per contrast
```

Both modules require `Run.deseq2: True` (they read the DESeq2 results files).

### Configuring gene-set collections

```yaml
GSEA:
  # MSigDB collections via msigdbr — no manual download needed.
  # Format: "CATEGORY" or "CATEGORY:SUBCATEGORY"
  collections:
    - H              # Hallmark
    - C2:CP:REACTOME # Reactome curated pathways
    - C5:GO:BP       # GO Biological Process

  min_size: 15       # Minimum genes in a set (after intersection with data)
  max_size: 500      # Maximum genes in a set
  nperm: 1000        # Permutations (fgseaSimple fallback; fgseaMultilevel is default)

  # Preferred ranking metric: "stat" (DESeq2 Wald statistic, recommended)
  # or "log2FoldChange" (fallback if stat is unavailable)
  rank_metric: stat

  # Optional: add your own .gmt files (paths relative to output_dir or absolute)
  custom_gmt_files: []
```

#### HOX gene sets (generated automatically)

The pipeline automatically creates `resources/generated_gmts/hox.gmt` from the
gene universe in your DESeq2 results.  It contains two gene sets:

| Set | Pattern | Description |
|-----|---------|-------------|
| `HOX_ALL`   | `^HOX[ABCD][0-9]+$` | All HOX cluster genes (A/B/C/D) |
| `HOXB_ONLY` | `^HOXB[0-9]+$`      | HOXB cluster genes only          |

These sets are merged with the msigdbr canonical sets and any custom GMTs before
fgsea runs.

#### Adding custom GMT files

Drop any `.gmt` file into your working directory (or specify an absolute path)
and add it to `GSEA.custom_gmt_files`:

```yaml
GSEA:
  custom_gmt_files:
    - resources/genesets/neuroblastoma_states.gmt
    - /path/to/mycn_targets.gmt
```

Gene sets in custom GMTs are merged with the canonical msigdbr sets.  On name
collision, the custom set takes precedence.

### GO enrichment

```yaml
GO:
  ontology:
    - BP     # Biological Process (recommended)
    # - MF   # Molecular Function (optional)
    # - CC   # Cellular Component (optional)
  padj_cutoff: 0.05
  min_gs_size: 10
  max_gs_size: 500
  gene_id_type: SYMBOL   # SYMBOL, ENSEMBL, or ENTREZID
```

GO enrichment runs separately for up- and down-regulated genes per contrast.
For `SYMBOL` and `ENSEMBL`, IDs are mapped to Entrez via `org.Hs.eg.db`; for
`ENTREZID`, IDs are used as-is. Raw and `simplify()`-reduced result tables are
written, along with mapping diagnostics and unmapped-ID CSVs. If mapping is
very low (<20% mapped universe or <5 mapped significant genes used for
enrichment), GO outputs are written as graceful empty files with warnings.

### Contrast comparisons (ΔNES and residual-rank GSEA)

The `ContrastComparisons` config section lets you "normalise" one contrast
against another.

```yaml
ContrastComparisons:
  - name: IFF_vs_FL_normalised_by_IFF_vs_Ctrl
    contrast_A: ATRX_IFF_vs_ATRX_FL        # numerator contrast
    contrast_B: ATRX_Ctrl_vs_ATRX_IFF      # baseline contrast
    delta_nes: True                          # compute ΔNES
    residual_rank: True                      # run fgsea on residual ranks
```

You can add as many comparison entries as needed.  Outputs appear under
`gsea_compare/{name}/`.

#### Interpreting ΔNES

**ΔNES = NES_A − NES_B** measures how much more (or less) a pathway is enriched
in contrast A relative to contrast B.

- **High positive ΔNES**: pathway is uniquely enriched in contrast A beyond what
  is already present in contrast B. Example: a high ΔNES for `HOXB_ONLY` in
  `ATRX_IFF_vs_ATRX_FL` means the HOXB upregulation is specific to the IFF
  condition and is not simply a general ATRX-depletion effect.
- **High negative ΔNES**: pathway is more enriched in contrast B; the signal in
  contrast A is partially explained by the baseline contrast.
- **Near zero ΔNES**: pathway is equally enriched in both contrasts.

#### Interpreting residual-rank GSEA

Residual ranks are computed as **rank_A − rank_B** (using the same metric,
preferably the DESeq2 Wald `stat`).  fgsea is then run on these residuals.

- Pathways with positive NES in the residual analysis are enriched among genes
  that move *more* in contrast A than in contrast B.
- This is a gene-level "subtraction" that highlights programs unique to the
  numerator contrast.

---

## Interpreting directionality

### DESeq2 log2FoldChange

Contrasts are defined as `[name, numerator, denominator]` in `DESeq2.contrasts`.
The **log2FoldChange (log2FC)** in every results table is always in the
numerator-over-denominator direction:

| log2FC | Meaning |
|--------|---------|
| **> 0** | Gene is **higher** in the **numerator** condition |
| **< 0** | Gene is **higher** in the **denominator** condition |

Both the volcano plot x-axis label and the MA plot subtitle state this
explicitly (e.g., `log2FC (treatment / control)`).

### Optional proteomics-integrated volcano plots

If `Proteomics.enabled: True`, the separate `DESeq2Proteomic` step writes
`deseq2/{contrast}/volcano_proteomic.pdf` for each contrast **that has an entry
in `Proteomics.deseq2_to_proteomics_comparison`**. Contrasts without a mapped
proteomics comparison are skipped (no overlay is produced for them). The plot is
restricted to genes that are significant in the mapped proteomics comparison
from the limma Excel sheet (`FDR <= Proteomics.fdr_threshold`). RNA-significant
points are all annotated and colored by RNA/protein direction agreement:

- **green**: RNA and proteomics change in the same direction
- **red**: RNA and proteomics change in opposite directions

### GSEA Normalised Enrichment Score (NES)

GSEA ranks genes by the DESeq2 Wald statistic (preferred) or log2FC,
so the ranked list inherits the same numerator-over-denominator orientation.

| NES | Meaning |
|-----|---------|
| **> 0** | Pathway is **enriched in the numerator** condition |
| **< 0** | Pathway is **enriched in the denominator** condition |

The GSEA dotplots are **faceted into two panels**:

- *Enriched in \<numerator\> (NES > 0)*
- *Enriched in \<denominator\> (NES < 0)*

The plot subtitle states the ranking metric and the contrast direction.
Each results CSV includes `numerator`, `denominator`, `rank_metric`, and
`direction_note` columns for unambiguous downstream interpretation.

### ΔNES (compare-GSEA)

**ΔNES = NES_A − NES_B** where A and B are two contrast-specific GSEA runs.

| ΔNES | Meaning |
|------|---------|
| **> 0** | Pathway more enriched in contrast A than contrast B |
| **< 0** | Pathway more enriched in contrast B than contrast A |
| **≈ 0** | Pathway equally enriched in both contrasts |

The ΔNES barplot subtitle and per-collection CSV files include the full
definition and the numerator/denominator of each contrast so that the
directionality is unambiguous.

### Residual-rank GSEA (compare-GSEA)

**residual_rank = rank_A − rank_B** (gene-level subtraction of ranked
statistics).  fgsea is run on these residuals.

| Residual NES | Meaning |
|--------------|---------|
| **> 0** | Pathway enriched among genes that move *more* in A than B |
| **< 0** | Pathway enriched among genes that move *more* in B than A |

The residual-rank CSV files include `contrast_A`, `contrast_B`, and a
`direction_note` column with the full definition.
