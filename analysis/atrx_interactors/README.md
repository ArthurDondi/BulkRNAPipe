# ATRX-interactor expression boxplots

Expression of the predicted ATRX FL / IFF interactors (`atrx_interactors.tsv`),
in two places:

- **The report** (`report/analysis_report.Rmd`, section "ATRX interactors"):
  one row per gene, the five SK-N-SH conditions next to the internal
  neuroblastoma cohort (`reanalysis.NB.bulkRNA`) split by ATRX / MYCN status.
  Point the report's `nb_expression` and `nb_metadata` params at the cohort
  files; the sample selection and every exclusion are listed in that section.
- **`01_plot_epicode.R`**: standalone epicode PDFs with statistics, not part of
  the Snakemake workflow.

Run the script from the repository root, in the pipeline's `deseq2` conda
environment (`workflow/envs/deseq2.yaml`: DESeq2, ggplot2, optparse). It needs
the pipeline's quantify/ and deseq2/ outputs.

```bash
Rscript analysis/atrx_interactors/01_plot_epicode.R \
     --project_dir /nobackup/lab_taschner-mandl/arthurdondi/projects/epicode
```

Or run `bash analysis/atrx_interactors/run_atrx_interactors.sh`. The path above
is the default; `--help` lists every option.

## Files

| File | Purpose |
|------|---------|
| `atrx_interactors.tsv` | gene, UniProt accession, list (FL/IFF), call (Yes/Maybe?), Ensembl gene ID (matches the patient matrix), previous HGNC symbols (fallback for the epicode GTF and the patient matrix) |
| `01_plot_epicode.R` | epicode boxplots, 5 conditions |
| `run_atrx_interactors.sh` | runs `01_plot_epicode.R` |
| `boxplot_helpers.R` | plotting and test helpers for `01_plot_epicode.R` |

Outputs go to
`/nobackup/lab_taschner-mandl/arthurdondi/projects/epicode/atrx_interactors/epicode/`
(change with `--outdir`): `overview.pdf` (one page per list, one panel per
gene), `per_gene.pdf` (one page per gene with statistics),
`expression_long.csv`, `stats.csv`, and `missing_genes.txt` if any gene could
not be found.

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

**Patient cohort.** The patient metadata (`metadata_table.csv`: omics_id,
material_general, tp, mna, atrx, exclusion, comment) is patient data: it stays
on the cluster next to the expression matrix and is never committed. Samples
with `exclusion` TRUE are left out of the report, with `comment` as the reason.
