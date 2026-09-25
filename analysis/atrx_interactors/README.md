# ATRX-interactor expression boxplots

Standalone analysis, not part of the Snakemake workflow. It plots the
expression of the predicted ATRX FL / IFF interactors (`atrx_interactors.tsv`)
in the epicode data and in GEO GSE94035.

Run everything from the repository root, in the pipeline's `deseq2` conda
environment (`workflow/envs/deseq2.yaml`: DESeq2, ggplot2, optparse).

```bash
# 1. GEO GSE94035: download + samplesheet (bash, wget/curl, ~12 MB)
bash analysis/atrx_interactors/00_download_GSE94035.sh \
     /nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret

# 2. epicode (needs the pipeline's quantify/ and deseq2/ outputs)
Rscript analysis/atrx_interactors/01_plot_epicode.R \
     --project_dir /nobackup/lab_taschner-mandl/arthurdondi/projects/epicode

# 3. GSE94035
Rscript analysis/atrx_interactors/02_plot_GSE94035.R \
     --geo_dir /nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret
```

All paths above are the defaults; `--help` lists every option.

## Files

| File | Purpose |
|------|---------|
| `atrx_interactors.tsv` | gene, UniProt accession, list (FL/IFF), call (Yes/Maybe?), Ensembl gene ID (used for GSE94035), previous HGNC symbols (fallback for the epicode GTF) |
| `00_download_GSE94035.sh` | downloads the processed matrix + series matrix + ENA run table and builds `metadata/samplesheet.tsv` |
| `01_plot_epicode.R` | epicode boxplots, 5 conditions |
| `02_plot_GSE94035.R` | GSE94035 boxplots, cell type x timepoint |
| `boxplot_helpers.R` | plotting and test helpers shared by both scripts |

Outputs (change with `--outdir`):

```
/nobackup/lab_taschner-mandl/arthurdondi/projects/epicode/atrx_interactors/
├── epicode/            # 01_plot_epicode.R
└── GSE94035_Fikret/    # 02_plot_GSE94035.R
```

Each folder holds `overview.pdf` (one page per list, one panel per gene),
`per_gene.pdf` (one page per gene with statistics), `expression_long.csv`,
`stats.csv`, and `missing_genes.txt` if any gene could not be found. The GEO
download itself stays in `/nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret`
(`raw/`, `metadata/`, `download.log`).

## Methods

**epicode.** log2(DESeq2 size-factor-normalised counts + 1), the same metric
as the ALCAM boxplots in the NK-NB study
(`code/34_20260603_paper_revision_figures.Rmd`). Size factors are estimated
exactly as in `workflow/scripts/deseq2.R`: the DESeq2 input matrix
(`design_qc/counts_substituted.txt` if present, else `quantify/counts.txt`),
all 15 samples, genes with >= 10 reads in total. Brackets are the DESeq2 contrasts
listed in `config/config_epicode.yaml` (`DESeq2.contrasts`, override with
`--config`) between plotted conditions, and show the Wilcoxon p on the plotted
values and the DESeq2 padj from `deseq2/<contrast>/results.csv`, taken as is:
BH over all genes DESeq2 tested in that contrast, not re-adjusted over this
gene list. With 3 vs 3
replicates the Wilcoxon rank-sum test cannot go below p = 0.1 (2 of the 20
possible rank arrangements), so the DESeq2 padj is the informative statistic.

**GSE94035.** The authors' deposited matrix (DESeq2-normalised FPM, log2,
GSNAP / GRCh37 / Ensembl 75), matched by Ensembl gene ID. Groups come from
the GEO sample titles; the single non-enriched DTC sample (`D07r2`) is left
out unless `DTC_relapse_unenriched` is added to `--groups`. Wilcoxon tests are
unpaired, although some patients contribute to several groups; the plots
and `stats.csv` also give BH-adjusted p, computed within each comparison
across the plotted genes (one family per comparison, as DESeq2 does per
contrast). Values are on a
different scale from epicode, so compare patterns within each dataset, not
absolute levels between them.

Notes on the deposited GEO data, handled by `00_download_GSE94035.sh`: the
matrix has CRLF line endings; two matrix columns differ from the GEO sample
descriptions (`D07r2` = "D07r, without enrichment", `D36NA` = "D36d"); the
GEO sample descriptions also carry an unrelated copy-pasted Ewing sarcoma text.
