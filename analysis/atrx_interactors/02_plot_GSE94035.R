#!/usr/bin/env Rscript
# Boxplots of ATRX-interactor expression in GEO GSE94035 (Rifatbegovic et al.,
# Int J Cancer 2018, PMID 28921546): stage M neuroblastoma primary tumours,
# bone-marrow disseminated tumour cells (DTCs) and DTC-depleted mononuclear
# cells (MNCs), at diagnosis and relapse.
#
# Run 00_download_GSE94035.sh first. Values are the authors' processed matrix
# as deposited: DESeq2-normalised FPM, log2 scale (GRCh37 / Ensembl 75), rows
# = Ensembl gene IDs, matched through the ensembl_gene_id column of the gene
# table. Not on the same scale as the epicode log2(normalised counts + 1).
#
# Statistics: two-sided Wilcoxon rank-sum p for each --comparisons pair, plus
# BH adjustment within each comparison across the plotted genes (as DESeq2
# adjusts per contrast across genes). Tests are unpaired, although some
# patients contribute to several groups.
#
# Outputs (in --outdir, default
# /nobackup/lab_taschner-mandl/arthurdondi/projects/epicode/atrx_interactors/GSE94035_Fikret):
#   overview.pdf, per_gene.pdf, expression_long.csv, stats.csv,
#   missing_genes.txt (if any)
#
# Usage (in the pipeline's deseq2 conda env):
#   Rscript analysis/atrx_interactors/02_plot_GSE94035.R \
#     [--geo_dir /nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret]

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--geo_dir", type = "character",
              default = "/nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret",
              help = "Directory created by 00_download_GSE94035.sh [default %default]"),
  make_option("--genes", type = "character", default = "",
              help = "Gene table [default atrx_interactors.tsv next to this script]"),
  make_option("--groups", type = "character",
              default = "Tumor_diagnosis,DTC_diagnosis,DTC_relapse,MNC_diagnosis,MNC_relapse",
              help = "samplesheet groups to plot, in x-axis order; add DTC_relapse_unenriched to show the single non-enriched sample [default %default]"),
  make_option("--comparisons", type = "character",
              default = "DTC_diagnosis:Tumor_diagnosis,DTC_diagnosis:MNC_diagnosis,DTC_relapse:DTC_diagnosis,DTC_relapse:MNC_relapse",
              help = "Wilcoxon pairs as group1:group2, comma-separated [default %default]"),
  make_option("--outdir", type = "character",
              default = "/nobackup/lab_taschner-mandl/arthurdondi/projects/epicode/atrx_interactors/GSE94035_Fikret",
              help = "Output directory [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))

here <- {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(here, "boxplot_helpers.R"))

geo_dir    <- args$geo_dir
genes_path <- if (nzchar(args$genes)) args$genes else file.path(here, "atrx_interactors.tsv")
outdir     <- args$outdir
groups     <- trimws(strsplit(args$groups, ",")[[1]])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ─── Data ────────────────────────────────────────────────────────────────────
sheet <- read.delim(file.path(geo_dir, "metadata", "samplesheet.tsv"),
                    stringsAsFactors = FALSE, colClasses = "character")
mat <- read.delim(gzfile(file.path(geo_dir, "raw", "GSE94035_Deseq_mat.txt.gz")),
                  row.names = 1, check.names = FALSE)
# The deposited file has CRLF line endings; guard against a stray \r.
colnames(mat) <- sub("\r$", "", colnames(mat))

unknown <- setdiff(groups, sheet$group)
if (length(unknown)) stop("Unknown --groups: ", paste(unknown, collapse = ", "),
                          "\nAvailable: ", paste(unique(sheet$group), collapse = ", "))
sheet <- sheet[sheet$group %in% groups, ]
if (!all(sheet$matrix_column %in% colnames(mat)))
  stop("samplesheet columns missing from matrix: ",
       paste(setdiff(sheet$matrix_column, colnames(mat)), collapse = ", "))
mat <- as.matrix(mat[, sheet$matrix_column, drop = FALSE])
sample_group <- factor(sheet$group, levels = groups)
message("Samples per group:")
print(table(sample_group))

genes <- read_gene_table(genes_path)
missing <- genes$gene[!genes$ensembl_gene_id %in% rownames(mat)]
if (length(missing)) {
  writeLines(missing, file.path(outdir, "missing_genes.txt"))
  message("Not in matrix: ", paste(missing, collapse = ", "))
}
genes <- genes[genes$ensembl_gene_id %in% rownames(mat), , drop = FALSE]

expr <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  data.frame(gene    = genes$gene[i],
             label   = genes$label[i],
             list    = genes$list[i],
             sample  = colnames(mat),
             patient = sheet$patient,
             group   = sample_group,
             value   = as.numeric(mat[genes$ensembl_gene_id[i], ]),
             stringsAsFactors = FALSE)
}))
expr$label <- factor(expr$label, levels = genes$label)
write.csv(transform(expr, label = NULL, log2_FPM = value, value = NULL),
          file.path(outdir, "expression_long.csv"), row.names = FALSE)

# ─── Statistics ──────────────────────────────────────────────────────────────
cmp <- do.call(rbind, lapply(strsplit(strsplit(args$comparisons, ",")[[1]], ":"), function(x) {
  data.frame(group1 = trimws(x[1]), group2 = trimws(x[2]), stringsAsFactors = FALSE)
}))
cmp <- cmp[cmp$group1 %in% groups & cmp$group2 %in% groups, , drop = FALSE]

stats <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  e <- expr[expr$gene == genes$gene[i], ]
  if (nrow(cmp) == 0) return(NULL)
  data.frame(gene     = genes$gene[i],
             list     = genes$list[i],
             call     = genes$call[i],
             group1   = cmp$group1,
             group2   = cmp$group2,
             n1       = sapply(cmp$group1, function(g) sum(e$group == g)),
             n2       = sapply(cmp$group2, function(g) sum(e$group == g)),
             log2_diff_median = mapply(function(a, b)
               median(e$value[e$group == a]) - median(e$value[e$group == b]),
               cmp$group1, cmp$group2),
             wilcox_p = mapply(function(a, b)
               wilcox_p(e$value[e$group == a], e$value[e$group == b]),
               cmp$group1, cmp$group2),
             stringsAsFactors = FALSE, row.names = NULL)
}))
if (!is.null(stats)) {
  # One BH family per comparison, across genes.
  stats$wilcox_padj_BH <- ave(stats$wilcox_p, stats$group1, stats$group2,
                              FUN = function(p) p.adjust(p, method = "BH"))
  write.csv(stats, file.path(outdir, "stats.csv"), row.names = FALSE)
}

# ─── Plots ───────────────────────────────────────────────────────────────────
colours <- setNames(okabe_ito[seq_along(groups)], groups)
ylab    <- expression(log[2] ~ "DESeq2-normalised FPM")
n_lab   <- paste(sprintf("%s n=%d", names(table(sample_group)), table(sample_group)),
                 collapse = ", ")

pdf(file.path(outdir, "overview.pdf"), width = 11, height = 8.5)
for (l in unique(genes$list)) {
  print(overview_page(droplevels(expr[expr$list == l, ]),
                      title = sprintf("ATRX %s interactors - GSE94035 neuroblastoma", l),
                      subtitle = n_lab,
                      ylab = ylab, colours = colours))
}
invisible(dev.off())

caption <- paste(
  sprintf("Wilcoxon rank-sum, two-sided, unpaired: p = raw p (in parentheses: BH-adjusted within that comparison across the %d genes).", nrow(genes)),
  "Data: GSE94035 processed matrix (log2 DESeq2-normalised FPM, GRCh37).",
  sep = "\n")
pdf(file.path(outdir, "per_gene.pdf"), width = 7, height = 7)
for (i in seq_len(nrow(genes))) {
  g  <- genes$gene[i]
  st <- if (is.null(stats)) NULL else stats[stats$gene == g, ]
  if (!is.null(st)) st$label <- sprintf("p=%s (%s)", fmt_p(st$wilcox_p), fmt_p(st$wilcox_padj_BH))
  print(gene_page(expr[expr$gene == g, ],
                  title    = g,
                  subtitle = sprintf("ATRX %s interactor: %s (UniProt %s, %s)",
                                     genes$list[i], genes$call[i], genes$uniprot[i],
                                     genes$ensembl_gene_id[i]),
                  ylab     = ylab, colours = colours,
                  pairs    = if (is.null(st)) NULL else st[, c("group1", "group2", "label")],
                  caption  = caption))
}
invisible(dev.off())

message("Done: ", nrow(genes), " genes plotted; outputs in ", outdir)
