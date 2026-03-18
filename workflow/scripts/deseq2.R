#!/usr/bin/env Rscript
# DESeq2 differential expression analysis for BulkRNAPipe
#
# Inputs:
#   --counts   : featureCounts output table (counts.txt)
#   --outdir   : directory for output files
#   --contrast : numerator_condition denominator_condition
#   --samples  : "sample1:condition1,sample2:condition2,..." mapping
#   --padj     : adjusted p-value threshold for significance (default 0.05)
#   --lfc      : absolute log2 fold-change threshold for volcano (default 1.0)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(dplyr)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",   type = "character", help = "featureCounts output file"),
  make_option("--outdir",   type = "character", help = "Output directory"),
  make_option("--contrast", type = "character", nargs = 2,
              help = "Numerator and denominator conditions (e.g. treatment control)"),
  make_option("--samples",  type = "character",
              help = "Comma-separated sample:condition pairs"),
  make_option("--padj",     type = "double",    default = 0.05,
              help = "Adjusted p-value threshold [default %default]"),
  make_option("--lfc",      type = "double",    default = 1.0,
              help = "Log2 fold-change threshold [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
args <- opt$options

# When --contrast takes two positional arguments after the flag they arrive
# in opt$args; handle both patterns.
if (is.null(args$contrast) && length(opt$args) >= 2) {
  contrast_num <- opt$args[1]
  contrast_den <- opt$args[2]
} else {
  vals <- strsplit(args$contrast, " ")[[1]]
  contrast_num <- vals[1]
  contrast_den <- vals[2]
}

outdir    <- args$outdir
padj_thr  <- args$padj
lfc_thr   <- args$lfc

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ─── Build sample info table ─────────────────────────────────────────────────
sample_pairs <- strsplit(args$samples, ",")[[1]]
sample_df    <- do.call(rbind, lapply(sample_pairs, function(x) {
  parts <- strsplit(x, ":")[[1]]
  data.frame(sample = parts[1], condition = parts[2], stringsAsFactors = FALSE)
}))
rownames(sample_df) <- sample_df$sample
sample_df$condition <- factor(sample_df$condition)

# ─── Load featureCounts output ───────────────────────────────────────────────
# featureCounts produces a header line starting with '#' and a data header.
raw <- read.table(args$counts, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE)

# Columns: Geneid, Chr, Start, End, Strand, Length, then one column per BAM.
count_cols <- colnames(raw)[7:ncol(raw)]
# Strip path prefix and .Aligned.sortedByCoord.out.bam suffix so column names
# match the sample names in the config.
clean_names <- sub(".*/", "", count_cols)
clean_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
counts <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_names

# Keep only samples present in sample_df
shared <- intersect(colnames(counts), rownames(sample_df))
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]

# ─── DESeq2 ──────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_df,
  design    = ~ condition
)

# Set the reference level to the denominator condition
dds$condition <- relevel(dds$condition, ref = contrast_den)

# Remove genes with very low counts (< 10 reads across all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# ─── Results ─────────────────────────────────────────────────────────────────
res <- results(
  dds,
  contrast  = c("condition", contrast_num, contrast_den),
  alpha     = padj_thr
)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::arrange(padj)

write.csv(res_df, file.path(outdir, "results.csv"), row.names = FALSE, quote = FALSE)

# Normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  tibble::rownames_to_column("gene_id")
write.csv(norm_counts, file.path(outdir, "normalized_counts.csv"),
          row.names = FALSE, quote = FALSE)

# ─── MA plot ─────────────────────────────────────────────────────────────────
pdf(file.path(outdir, "ma_plot.pdf"), width = 6, height = 5)
plotMA(res, alpha = padj_thr, main = paste(contrast_num, "vs", contrast_den))
dev.off()

# ─── Volcano plot ────────────────────────────────────────────────────────────
volcano_df <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      padj < padj_thr & abs(log2FoldChange) >= lfc_thr ~ "Significant",
      TRUE                                              ~ "Not significant"
    ),
    label = ifelse(
      padj < padj_thr & abs(log2FoldChange) >= lfc_thr,
      gene_id, NA_character_
    )
  )

p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj),
                             colour = significance, label = label)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = c("Significant" = "#E41A1C",
                                 "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title  = paste("Volcano:", contrast_num, "vs", contrast_den),
    x      = expression(log[2]~"fold change"),
    y      = expression(-log[10]~"adjusted p-value"),
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "volcano.pdf"), plot = p, width = 7, height = 6)

message("DESeq2 analysis complete. Results written to: ", outdir)
