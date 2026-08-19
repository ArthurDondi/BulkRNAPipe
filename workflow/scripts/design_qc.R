#!/usr/bin/env Rscript
# design_qc.R — library-level QC aimed at the confounders that matter when
# groups differ by more than the intended variable.
#
# Everything here plots all samples/replicates together, because the question is
# always "does this vary with group in a way that has nothing to do with the
# biology?" — which is invisible one sample at a time.
#
# Outputs (design_qc/):
#   library_size.pdf            assigned read pairs per sample
#   assignment_rates.pdf        featureCounts assignment breakdown, % of reads
#   genes_detected.pdf          number of genes with >= 1 count
#   marker_expression.pdf       normalised expression of reporter/transgene features
#   marker_expression.csv       the same numbers as a table
#   sample_correlation.pdf      Spearman correlation between all samples (VST)
#   sample_distance.pdf         Euclidean distance between all samples (VST)
#   design_qc_summary.csv       one row per sample with all of the above
#
# Read this alongside the MultiQC reports under QC/ — duplication, adapter and
# coverage-uniformity metrics live there and are not recomputed here.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--counts",  type = "character", help = "featureCounts output file"),
  make_option("--summary", type = "character", help = "featureCounts .summary file"),
  make_option("--outdir",  type = "character", help = "Output directory"),
  make_option("--samples", type = "character",
              help = "Comma-separated sample:condition pairs"),
  make_option("--marker_genes", type = "character", default = "",
              help = "Comma-separated gene IDs to plot individually (reporters, transgene)")
)

args   <- parse_args(OptionParser(option_list = option_list))
outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

split_csv <- function(x) {
  if (is.null(x)) return(character(0))
  parts <- trimws(strsplit(x, ",")[[1]])
  parts[nzchar(parts)]
}

clean_sample_names <- function(x) {
  x <- sub(".*/", "", x)
  sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", x)
}

# ─── Sample table ────────────────────────────────────────────────────────────
sample_df <- do.call(rbind, lapply(split_csv(args$samples), function(p) {
  parts <- strsplit(p, ":")[[1]]
  data.frame(sample = parts[1], condition = parts[2], stringsAsFactors = FALSE)
}))
rownames(sample_df) <- sample_df$sample

# ─── Counts ──────────────────────────────────────────────────────────────────
raw    <- read.table(args$counts, header = TRUE, sep = "\t", comment.char = "#",
                     check.names = FALSE)
counts <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_sample_names(colnames(raw)[7:ncol(raw)])

shared    <- intersect(colnames(counts), rownames(sample_df))
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]

# Plot samples grouped by condition, in config order.
sample_levels <- rownames(sample_df)[order(match(sample_df$condition,
                                                 unique(sample_df$condition)))]
sample_df$sample    <- factor(rownames(sample_df), levels = sample_levels)
sample_df$condition <- factor(sample_df$condition, levels = unique(sample_df$condition))

marker_genes <- split_csv(args$marker_genes)
markers_present <- marker_genes[tolower(marker_genes) %in% tolower(rownames(counts))]
for (m in setdiff(marker_genes, markers_present)) {
  message("WARNING: marker gene '", m, "' is not in the count matrix.")
}
# Recover the exact spelling used in the matrix.
markers_present <- rownames(counts)[tolower(rownames(counts)) %in% tolower(markers_present)]

# ─── Normalisation (markers excluded from the size factors) ──────────────────
# Reporter and transgene features are present in some groups and absent in
# others by construction. Leaving them in the size-factor estimation lets the
# construct shift the normalisation of every other gene.
dds <- DESeqDataSetFromMatrix(
  countData = counts[!rownames(counts) %in% markers_present, , drop = FALSE],
  colData   = sample_df,
  design    = ~ 1
)
dds  <- estimateSizeFactors(dds)
sf   <- sizeFactors(dds)
norm <- sweep(counts, 2, sf, "/")     # applied to the FULL matrix, markers included

vsd <- vst(dds[rowSums(counts(dds)) >= 10, ], blind = TRUE)
mat <- assay(vsd)

# ─── Per-sample summary ──────────────────────────────────────────────────────
summary_df <- data.frame(
  sample          = sample_df$sample,
  condition       = sample_df$condition,
  assigned_reads  = colSums(counts),
  size_factor     = round(sf[rownames(sample_df)], 4),
  genes_detected  = colSums(counts > 0),
  stringsAsFactors = FALSE
)
if (length(markers_present) > 0) {
  summary_df$pct_reads_in_markers <-
    round(100 * colSums(counts[markers_present, , drop = FALSE]) / colSums(counts), 4)
}

theme_qc <- theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

bar_plot <- function(df, yvar, title, ylab, outfile, ylabels = waiver()) {
  p <- ggplot(df, aes(x = sample, y = .data[[yvar]], fill = condition)) +
    geom_col() +
    labs(title = title, x = NULL, y = ylab, fill = "Condition") +
    scale_y_continuous(labels = ylabels) +
    theme_qc
  ggsave(outfile, plot = p, width = 8, height = 5)
  message("Written: ", outfile)
}

bar_plot(summary_df, "assigned_reads",
         "Reads assigned to genes",
         "assigned reads", file.path(outdir, "library_size.pdf"),
         ylabels = function(x) format(x, big.mark = ",", scientific = FALSE))

bar_plot(summary_df, "genes_detected",
         "Genes detected (>= 1 count)",
         "genes", file.path(outdir, "genes_detected.pdf"))

# ─── featureCounts assignment breakdown ──────────────────────────────────────
if (!is.null(args$summary) && file.exists(args$summary)) {
  sm <- read.table(args$summary, header = TRUE, sep = "\t", check.names = FALSE)
  colnames(sm)[-1] <- clean_sample_names(colnames(sm)[-1])
  sm <- sm[, c("Status", intersect(colnames(sm), rownames(sample_df))), drop = FALSE]
  sm <- sm[rowSums(sm[, -1, drop = FALSE]) > 0, , drop = FALSE]

  long <- do.call(rbind, lapply(seq_len(nrow(sm)), function(i) {
    data.frame(status = sm$Status[i],
               sample = colnames(sm)[-1],
               count  = as.numeric(sm[i, -1]),
               stringsAsFactors = FALSE)
  }))
  totals      <- tapply(long$count, long$sample, sum)
  long$pct    <- 100 * long$count / totals[long$sample]
  long$sample <- factor(long$sample, levels = sample_levels)
  long$status <- factor(long$status,
                        levels = c("Assigned", setdiff(unique(long$status), "Assigned")))

  p <- ggplot(long, aes(x = sample, y = pct, fill = status)) +
    geom_col() +
    labs(title = "featureCounts read assignment",
         subtitle = "a group-specific shift here is a library artefact, not biology",
         x = NULL, y = "% of reads", fill = "Status") +
    theme_qc
  ggsave(file.path(outdir, "assignment_rates.pdf"), plot = p, width = 9, height = 5)
  message("Written: ", file.path(outdir, "assignment_rates.pdf"))

  assigned <- long[long$status == "Assigned", ]
  summary_df$pct_assigned <- round(assigned$pct[match(summary_df$sample, assigned$sample)], 2)
} else {
  message("No featureCounts summary file found - skipping assignment_rates.pdf")
}

# ─── Reporter / transgene expression ─────────────────────────────────────────
if (length(markers_present) > 0) {
  mk <- do.call(rbind, lapply(markers_present, function(g) {
    data.frame(gene = g, sample = colnames(norm), value = as.numeric(norm[g, ]),
               stringsAsFactors = FALSE)
  }))
  mk$sample    <- factor(mk$sample, levels = sample_levels)
  mk$condition <- sample_df$condition[match(mk$sample, sample_df$sample)]
  mk$gene      <- factor(mk$gene, levels = markers_present)

  p <- ggplot(mk, aes(x = sample, y = value + 1, fill = condition)) +
    geom_col() +
    facet_wrap(~ gene, scales = "free_y", ncol = 1) +
    scale_y_log10() +
    labs(title = "Reporter and transgene expression",
         subtitle = "normalised counts + 1, log10 scale; size factors exclude these features",
         x = NULL, y = "normalised count + 1", fill = "Condition") +
    theme_qc
  ggsave(file.path(outdir, "marker_expression.pdf"), plot = p,
         width = 9, height = 2.6 * length(markers_present) + 1.5, limitsize = FALSE)
  message("Written: ", file.path(outdir, "marker_expression.pdf"))

  wide <- as.data.frame(round(norm[markers_present, , drop = FALSE], 3))
  wide <- cbind(gene_id = rownames(wide), wide)
  write.csv(wide, file.path(outdir, "marker_expression.csv"), row.names = FALSE)
} else {
  message("No marker genes present - skipping marker_expression.pdf")
  write.csv(data.frame(gene_id = character(0)),
            file.path(outdir, "marker_expression.csv"), row.names = FALSE)
}

# ─── Sample-sample relationships ─────────────────────────────────────────────
heat <- function(m, title, subtitle, fillname, outfile, palette_low, palette_high) {
  ord <- sample_levels[sample_levels %in% colnames(m)]
  m   <- m[ord, ord, drop = FALSE]
  df  <- expand.grid(x = factor(ord, levels = ord), y = factor(ord, levels = rev(ord)),
                     stringsAsFactors = FALSE)
  df$value <- as.numeric(m[cbind(match(as.character(df$y), ord),
                                 match(as.character(df$x), ord))])
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = palette_low, high = palette_high, name = fillname) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    coord_fixed() +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  ggsave(outfile, plot = p, width = 8, height = 7)
  message("Written: ", outfile)
}

cor_mat <- cor(mat, method = "spearman")
heat(cor_mat, "Sample-sample Spearman correlation (VST)",
     "replicates should be the tightest block on the diagonal",
     "rho", file.path(outdir, "sample_correlation.pdf"), "white", "#2166AC")

dist_mat <- as.matrix(dist(t(mat)))
heat(dist_mat, "Sample-sample Euclidean distance (VST)",
     "distance between group blocks is the effect you are trying to interpret",
     "distance", file.path(outdir, "sample_distance.pdf"), "#B2182B", "white")

write.csv(cbind(sample = rownames(cor_mat), as.data.frame(round(cor_mat, 4))),
          file.path(outdir, "sample_correlation.csv"), row.names = FALSE)

# ─── Summary table ───────────────────────────────────────────────────────────
write.csv(summary_df, file.path(outdir, "design_qc_summary.csv"), row.names = FALSE)
message("Written: ", file.path(outdir, "design_qc_summary.csv"))
print(summary_df, row.names = FALSE)
