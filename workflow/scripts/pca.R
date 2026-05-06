#!/usr/bin/env Rscript
# PCA plot for BulkRNAPipe – all samples, PC1 vs PC2 with % variance explained
#
# Inputs:
#   --counts  : featureCounts output table (counts.txt)
#   --outdir  : directory for the output PDF
#   --samples : "sample1:condition1,sample2:condition2,..." mapping

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",  type = "character", help = "featureCounts output file"),
  make_option("--outdir",  type = "character", help = "Output directory"),
  make_option("--samples", type = "character",
              help = "Comma-separated sample:condition pairs")
)

args <- parse_args(OptionParser(option_list = option_list))

outdir <- args$outdir
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
count_cols  <- colnames(raw)[7:ncol(raw)]
# Strip path prefix and .Aligned.sortedByCoord.out.bam suffix so column names
# match the sample names in the config.
clean_names <- sub(".*/", "", count_cols)
clean_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
counts      <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_names

# Keep only samples present in sample_df
shared    <- intersect(colnames(counts), rownames(sample_df))
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]

# ─── DESeq2 dataset + variance-stabilising transformation ────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_df,
  design    = ~ condition
)

# Remove genes with very low counts (< 10 reads across all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Blind VST for sample-level QC (design information not used)
vst_data <- vst(dds, blind = TRUE)

# ─── PCA (top 500 most-variable genes by default) ────────────────────────────
pca_df  <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_df, "percentVar"), 1)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title  = "PCA - all samples",
    x      = paste0("PC1: ", pct_var[1], "% variance"),
    y      = paste0("PC2: ", pct_var[2], "% variance"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(outdir, "pca.pdf"), plot = p, width = 7, height = 6)

message("PCA plot written to: ", file.path(outdir, "pca.pdf"))
