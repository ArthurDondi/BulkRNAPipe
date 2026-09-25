#!/usr/bin/env Rscript
# Boxplots of ATRX-interactor expression in the epicode bulk RNA-seq data.
#
# Values: log2(DESeq2 size-factor-normalised counts + 1), computed from the
# same count matrix and low-count filter that workflow/scripts/deseq2.R uses,
# so the normalisation matches the pipeline's DESeq2 runs.
#
# Count matrix (the DESeq2 input):
#   <project_dir>/design_qc/counts_substituted.txt when it exists (pooled
#   ATRX_Total, EGFP/mCherry removed - see DesignQCSubstituteTotal),
#   otherwise <project_dir>/quantify/counts.txt.
#
# Statistics, one bracket per DESeq2 contrast listed in the config
# (DESeq2.contrasts) between plotted conditions:
#   W    two-sided Wilcoxon rank-sum p on the plotted values (with 3 vs 3
#        replicates the smallest attainable p is 0.1)
#   padj DESeq2 padj from deseq2/<contrast>/results.csv, taken as is: BH over
#        all genes tested in that contrast (not re-adjusted over this list)
#
# Outputs (in --outdir, default <project_dir>/atrx_interactors/epicode):
#   overview.pdf         one page per list (FL, IFF), one panel per gene
#   per_gene.pdf         one page per gene with the statistics
#   expression_long.csv  gene, sample, condition, log2 normalised count
#   stats.csv            gene x contrast: Wilcoxon p, DESeq2 log2FC / padj
#   missing_genes.txt    genes not found in the count matrix (if any)
#
# Usage (in the pipeline's deseq2 conda env):
#   Rscript analysis/atrx_interactors/01_plot_epicode.R \
#     [--project_dir /nobackup/lab_taschner-mandl/arthurdondi/projects/epicode]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--project_dir", type = "character",
              default = "/nobackup/lab_taschner-mandl/arthurdondi/projects/epicode",
              help = "Pipeline output_dir [default %default]"),
  make_option("--counts", type = "character", default = "",
              help = "featureCounts-format matrix; default: design_qc/counts_substituted.txt if present, else quantify/counts.txt"),
  make_option("--config", type = "character", default = "",
              help = "Pipeline config listing DESeq2.contrasts [default config/config_epicode.yaml in this repo]"),
  make_option("--deseq2_dir", type = "character", default = "",
              help = "Directory with <contrast>/results.csv [default <project_dir>/deseq2]"),
  make_option("--genes", type = "character", default = "",
              help = "Gene table [default atrx_interactors.tsv next to this script]"),
  make_option("--conditions", type = "character",
              default = "TP53,E6,EmptyVector,ATRX_FL,ATRX_IFF",
              help = "Conditions to plot, in x-axis order [default %default]"),
  make_option("--outdir", type = "character", default = "",
              help = "Output directory [default <project_dir>/atrx_interactors/epicode]")
)
args <- parse_args(OptionParser(option_list = option_list))

here <- {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(here, "boxplot_helpers.R"))

proj <- args$project_dir
counts_path <- if (nzchar(args$counts)) args$counts else {
  sub <- file.path(proj, "design_qc", "counts_substituted.txt")
  if (file.exists(sub)) sub else file.path(proj, "quantify", "counts.txt")
}
deseq2_dir <- if (nzchar(args$deseq2_dir)) args$deseq2_dir else file.path(proj, "deseq2")
config_path <- if (nzchar(args$config)) args$config else
  normalizePath(file.path(here, "..", "..", "config", "config_epicode.yaml"), mustWork = FALSE)
genes_path <- if (nzchar(args$genes)) args$genes else file.path(here, "atrx_interactors.tsv")
outdir     <- if (nzchar(args$outdir)) args$outdir else file.path(proj, "atrx_interactors", "epicode")
conditions <- trimws(strsplit(args$conditions, ",")[[1]])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Count matrix : ", counts_path)
message("Config       : ", config_path)
message("DESeq2 dir   : ", deseq2_dir)
message("Gene table   : ", genes_path)
message("Output dir   : ", outdir)

# ─── Counts (same parsing as workflow/scripts/deseq2.R) ──────────────────────
raw <- read.table(counts_path, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE)
counts <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "",
                        sub(".*/", "", colnames(counts)))

# Sample names are <condition>_<replicate> (config_epicode.yaml `samples:`).
sample_cond <- sub("_[0-9]+$", "", colnames(counts))
keep_s  <- sample_cond %in% conditions
if (!any(keep_s)) stop("No samples match --conditions: ", args$conditions)
counts      <- counts[, keep_s, drop = FALSE]
sample_cond <- factor(sample_cond[keep_s], levels = conditions)
message("Samples: ", paste(sprintf("%s (%s)", colnames(counts), sample_cond), collapse = ", "))

# ─── Normalisation, mirroring deseq2.R ───────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = data.frame(condition = sample_cond, row.names = colnames(counts)),
  design    = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)
log_norm <- log2(counts(dds, normalized = TRUE) + 1)

# ─── Genes ───────────────────────────────────────────────────────────────────
genes <- read_gene_table(genes_path)
resolve_id <- function(sym, prev) {
  cand <- c(sym, trimws(strsplit(prev, ";")[[1]]))
  cand <- cand[nzchar(cand)]
  hit  <- cand[cand %in% rownames(log_norm)]
  if (length(hit)) hit[1] else NA_character_
}
genes$matrix_id <- mapply(resolve_id, genes$gene, genes$previous_symbols)
missing <- genes$gene[is.na(genes$matrix_id)]
if (length(missing)) {
  in_raw <- missing %in% rownames(counts)
  msg <- sprintf("%s\t%s", missing,
                 ifelse(in_raw, "present but < 10 reads in total (filtered)",
                        "not in count matrix"))
  writeLines(c("gene\treason", msg), file.path(outdir, "missing_genes.txt"))
  message("Missing genes:\n  ", paste(msg, collapse = "\n  "))
}
genes <- genes[!is.na(genes$matrix_id), , drop = FALSE]

expr <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  data.frame(gene      = genes$gene[i],
             label     = genes$label[i],
             list      = genes$list[i],
             sample    = colnames(log_norm),
             group     = sample_cond,
             value     = as.numeric(log_norm[genes$matrix_id[i], ]),
             stringsAsFactors = FALSE)
}))
expr$label <- factor(expr$label, levels = genes$label)
write.csv(transform(expr, label = NULL, condition = group, group = NULL,
                    log2_norm_count = value, value = NULL),
          file.path(outdir, "expression_long.csv"), row.names = FALSE)

# ─── DESeq2 contrasts, from the config ───────────────────────────────────────
# Read from DESeq2.contrasts ("- [name, numerator, denominator]") rather than
# globbing deseq2/: output folders from renamed contrasts stay on disk and
# would otherwise be picked up as duplicates. Only lines of exactly that shape
# are matched; commented-out contrasts are skipped.
cfg_lines <- readLines(config_path)
m <- regmatches(cfg_lines, regexec(
  "^\\s*-\\s*\\[\\s*([^],#]+?)\\s*,\\s*([^],#]+?)\\s*,\\s*([^],#]+?)\\s*\\]", cfg_lines, perl = TRUE))
m <- do.call(rbind, m[lengths(m) == 4])
contrasts <- if (is.null(m)) NULL else
  data.frame(contrast = m[, 2], group1 = m[, 3], group2 = m[, 4],
             stringsAsFactors = FALSE)
if (!is.null(contrasts)) {
  contrasts <- contrasts[contrasts$group1 %in% conditions &
                         contrasts$group2 %in% conditions, , drop = FALSE]
  contrasts <- contrasts[!duplicated(contrasts[, c("group1", "group2")]), , drop = FALSE]
  contrasts$file <- file.path(deseq2_dir, contrasts$contrast, "results.csv")
  absent <- !file.exists(contrasts$file)
  if (any(absent)) {
    message("No results.csv for: ", paste(contrasts$contrast[absent], collapse = ", "),
            " - Wilcoxon only for these.")
    contrasts$file[absent] <- NA
  }
}
if (is.null(contrasts) || nrow(contrasts) == 0) {
  message("No DESeq2 contrasts between the plotted conditions in ", config_path,
          "; Wilcoxon only for all condition pairs.")
  cmb <- t(combn(conditions, 2))
  contrasts <- data.frame(contrast = paste0(cmb[, 2], "_vs_", cmb[, 1]),
                          group1 = cmb[, 2], group2 = cmb[, 1], file = NA,
                          stringsAsFactors = FALSE)
}
message("Contrasts (numerator vs denominator): ",
        paste(sprintf("%s [%s vs %s]", contrasts$contrast, contrasts$group1, contrasts$group2),
              collapse = ", "))

deseq_res <- lapply(setNames(contrasts$file, contrasts$contrast), function(f) {
  if (is.na(f)) return(NULL)
  r <- read.csv(f, stringsAsFactors = FALSE)
  r[match(genes$matrix_id, r$gene_id), c("log2FoldChange", "padj")]
})

stats <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  e <- expr[expr$gene == genes$gene[i], ]
  do.call(rbind, lapply(seq_len(nrow(contrasts)), function(j) {
    c1 <- contrasts$group1[j]; c2 <- contrasts$group2[j]
    d  <- deseq_res[[contrasts$contrast[j]]]
    data.frame(gene          = genes$gene[i],
               list          = genes$list[i],
               call          = genes$call[i],
               contrast      = contrasts$contrast[j],
               group1        = c1,
               group2        = c2,
               wilcox_p      = wilcox_p(e$value[e$group == c1], e$value[e$group == c2]),
               deseq2_log2FC = if (is.null(d)) NA_real_ else d$log2FoldChange[i],
               deseq2_padj   = if (is.null(d)) NA_real_ else d$padj[i],
               stringsAsFactors = FALSE)
  }))
}))
write.csv(stats, file.path(outdir, "stats.csv"), row.names = FALSE)

# ─── Plots ───────────────────────────────────────────────────────────────────
colours <- setNames(okabe_ito[seq_along(conditions)], conditions)
ylab    <- expression(log[2] ~ "(normalised counts + 1)")

pdf(file.path(outdir, "overview.pdf"), width = 11, height = 8.5)
for (l in unique(genes$list)) {
  print(overview_page(droplevels(expr[expr$list == l, ]),
                      title = sprintf("ATRX %s interactors - epicode (SK-N-SH)", l),
                      ylab = ylab, colours = colours))
}
invisible(dev.off())

caption <- paste(
  "W: two-sided Wilcoxon rank-sum p (n = 3 vs 3: minimum attainable p = 0.1).",
  "padj: DESeq2 BH-adjusted p for that contrast.",
  sep = "\n")
pdf(file.path(outdir, "per_gene.pdf"), width = 7, height = 7)
for (i in seq_len(nrow(genes))) {
  g  <- genes$gene[i]
  st <- stats[stats$gene == g, ]
  st$label <- sprintf("W p=%s | padj=%s", fmt_p(st$wilcox_p), fmt_p(st$deseq2_padj))
  df <- expr[expr$gene == g, ]
  print(gene_page(df,
                  title    = g,
                  subtitle = sprintf("ATRX %s interactor: %s (UniProt %s)%s",
                                     genes$list[i], genes$call[i], genes$uniprot[i],
                                     if (genes$matrix_id[i] != g) paste0("; matrix ID ", genes$matrix_id[i]) else ""),
                  ylab     = ylab, colours = colours,
                  pairs    = st[, c("group1", "group2", "label")],
                  caption  = caption))
}
invisible(dev.off())

message("Done: ", nrow(genes), " genes plotted; outputs in ", outdir)
