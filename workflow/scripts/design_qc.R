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
#                               (optionally topped by a pooled total and a tag-only
#                               row, counted with -M --fraction over the same BAMs)
#   marker_expression_linear.pdf  the same panel on a linear axis
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
              help = "Comma-separated reporter/transgene IDs; plotted AND excluded from size factors"),
  make_option("--goi_genes", type = "character", default = "",
              help = "Comma-separated endogenous genes of interest; plotted but kept in the size factors"),
  make_option("--total_counts", type = "character", default = "",
              help = "featureCounts output for the pooled total feature (optional)"),
  make_option("--total_name", type = "character", default = "ATRX_Total",
              help = "Feature name inside --total_counts"),
  make_option("--tag_counts", type = "character", default = "",
              help = "featureCounts output for the tag-only feature (optional)"),
  make_option("--tag_name", type = "character", default = "Tag",
              help = "Feature name inside --tag_counts"),
  make_option("--gene_labels", type = "character", default = "",
              help = "Comma-separated ID=Label pairs; renames facet strips only")
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

# Display-only relabelling for the panels. Everything downstream - the CSVs, the
# count matrix, the summary table - keeps the real feature IDs, so renaming here
# cannot silently change what is being counted.
gene_labels <- local({
  pairs <- split_csv(args$gene_labels)
  out   <- character(0)
  for (p in pairs) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2 || !nzchar(kv[1]) || !nzchar(kv[2])) {
      message("WARNING: ignoring malformed gene label '", p, "' (expected ID=Label).")
      next
    }
    out[kv[1]] <- kv[2]
  }
  out
})
display_label <- function(x) {
  hit <- match(x, names(gene_labels))
  ifelse(is.na(hit), x, gene_labels[hit])
}

marker_genes <- split_csv(args$marker_genes)
markers_present <- marker_genes[tolower(marker_genes) %in% tolower(rownames(counts))]
for (m in setdiff(marker_genes, markers_present)) {
  message("WARNING: marker gene '", m, "' is not in the count matrix.")
}
# Recover the exact spelling used in the matrix.
markers_present <- rownames(counts)[tolower(rownames(counts)) %in% tolower(markers_present)]

# Genes of interest are ordinary endogenous genes, so unlike the markers above
# they stay in the size-factor estimation - dropping real genes from the
# normalisation would be wrong, and they are not structurally absent in any group.
resolve_genes <- function(wanted, label) {
  present <- rownames(counts)[tolower(rownames(counts)) %in% tolower(wanted)]
  for (g in wanted[!tolower(wanted) %in% tolower(present)]) {
    message("WARNING: ", label, " '", g, "' is not in the count matrix.")
  }
  # Return in the order the user listed them.
  present[order(match(tolower(present), tolower(wanted)))]
}
goi_present <- resolve_genes(split_csv(args$goi_genes), "gene of interest")

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

# ─── Pooled and tag-only features (counted separately, no re-alignment) ──────
# These come from their own featureCounts runs over the same BAMs, with
# -Q 0 -M --fraction, so a fragment that aligns equally well to the endogenous
# locus and to one or more vector contigs is counted once instead of dropped.
# They are stacked onto `norm` for the marker panel only: they never enter
# `counts`, so size factors, the VST and pct_reads_in_markers are untouched.
read_feature <- function(path, feature, label) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    message("WARNING: ", label, " count file '", path, "' not found - skipping.")
    return(NULL)
  }
  tb <- read.table(path, header = TRUE, sep = "\t", comment.char = "#",
                   check.names = FALSE)
  m  <- as.matrix(tb[, 7:ncol(tb), drop = FALSE])
  rownames(m) <- tb[[1]]
  colnames(m) <- clean_sample_names(colnames(m))
  if (!feature %in% rownames(m)) {
    message("WARNING: feature '", feature, "' is not in ", path, " - skipping.")
    return(NULL)
  }
  absent <- setdiff(colnames(norm), colnames(m))
  if (length(absent) > 0) {
    message("WARNING: ", label, " is missing sample(s) ",
            paste(absent, collapse = ", "), " - skipping.")
    return(NULL)
  }
  v <- m[feature, colnames(norm), drop = FALSE]
  if (all(v == 0)) {
    message("WARNING: ", label, " ('", feature, "') is zero in every sample. ",
            "Check that the annotation it was counted from matches the reference.")
  }
  v
}

extra_rows <- list()
for (spec in list(
       list(path = args$total_counts, feature = args$total_name,
            label = "pooled total feature"),
       list(path = args$tag_counts, feature = args$tag_name,
            label = "tag feature"))) {
  v <- read_feature(spec$path, spec$feature, spec$label)
  if (!is.null(v)) extra_rows[[spec$feature]] <- v
}

marker_panel_genes <- markers_present
if (length(extra_rows) > 0) {
  extra <- do.call(rbind, extra_rows)
  rownames(extra) <- names(extra_rows)
  # Same size factors as everything else in the panel, so the bars sit on a
  # comparable scale even though the counting rule differs.
  extra <- sweep(extra, 2, sf[colnames(extra)], "/")
  norm  <- rbind(extra, norm)
  marker_panel_genes <- c(rownames(extra), markers_present)
}

# ─── Per-sample summary ──────────────────────────────────────────────────────
summary_df <- data.frame(
  sample          = sample_df$sample,
  condition       = sample_df$condition,
  assigned_reads  = colSums(counts),
  size_factor     = round(sf[rownames(sample_df)], 4),
  genes_detected  = colSums(counts > 0),
  stringsAsFactors = FALSE
)
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

# ─── Per-gene expression panels ──────────────────────────────────────────────
# log_scale = FALSE draws the same panel on a linear axis. Both views are worth
# having: log keeps a feature that is near zero in one group and thousands in
# another on the same strip at all, linear shows how large the differences
# actually are - on log10 a 10x and a 100x difference look similarly modest.
# out_csv = NULL writes no table (the linear view reuses the log view's).
gene_panel <- function(genes, title, subtitle, out_pdf, out_csv,
                       log_scale = TRUE) {
  if (length(genes) == 0) {
    message("No genes to plot for '", title, "' - writing an empty table.")
    if (!is.null(out_csv)) {
      write.csv(data.frame(gene_id = character(0)), out_csv, row.names = FALSE)
    }
    pdf(out_pdf, width = 7, height = 4); plot.new()
    text(0.5, 0.5, paste0(title, ": no genes found"), cex = 1.1); dev.off()
    return(invisible(NULL))
  }
  df <- do.call(rbind, lapply(genes, function(g) {
    data.frame(gene = g, sample = colnames(norm), value = as.numeric(norm[g, ]),
               stringsAsFactors = FALSE)
  }))
  df$sample    <- factor(df$sample, levels = sample_levels)
  df$condition <- sample_df$condition[match(df$sample, sample_df$sample)]
  # Facet strips only; `genes` and the CSV below stay on the real IDs.
  df$gene      <- factor(display_label(df$gene), levels = display_label(genes))

  # The +1 offset exists only so zeros survive the log transform; on a linear
  # axis it would be a silent distortion, so the raw value is plotted instead.
  df$y <- if (log_scale) df$value + 1 else df$value
  p <- ggplot(df, aes(x = sample, y = y, fill = condition)) +
    geom_col() +
    facet_wrap(~ gene, scales = "free_y", ncol = 1) +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = if (log_scale) "normalised count + 1" else "normalised count",
         fill = "Condition") +
    theme_qc
  if (log_scale) p <- p + scale_y_log10()
  ggsave(out_pdf, plot = p, width = 9,
         height = 2.6 * length(genes) + 1.5, limitsize = FALSE)
  message("Written: ", out_pdf)

  if (!is.null(out_csv)) {
    wide <- as.data.frame(round(norm[genes, , drop = FALSE], 3))
    write.csv(cbind(gene_id = rownames(wide), wide), out_csv, row.names = FALSE)
  }
}

marker_subtitle <- "normalised counts + 1, log10; excluded from the size factors"
if (length(extra_rows) > 0) {
  marker_subtitle <- paste0(
    marker_subtitle, "\n",
    paste(names(extra_rows), collapse = " and "),
    ": every alignment weighted 1/NH, multimappers included \u2014 ",
    "not on the same footing as the unique-only rows below")
}

gene_panel(marker_panel_genes,
           "Reporter and transgene expression",
           marker_subtitle,
           file.path(outdir, "marker_expression.pdf"),
           file.path(outdir, "marker_expression.csv"))

# Same numbers, linear axis. On log10 the transgene rows read as modest
# differences; linear is what shows their real size.
gene_panel(marker_panel_genes,
           "Reporter and transgene expression (linear scale)",
           sub("^normalised counts \\+ 1, log10",
               "normalised counts, linear axis", marker_subtitle),
           file.path(outdir, "marker_expression_linear.pdf"),
           out_csv = NULL, log_scale = FALSE)

gene_panel(goi_present,
           "Genes of interest",
           "normalised counts + 1, log10; endogenous genes, included in the size factors",
           file.path(outdir, "goi_expression.pdf"),
           file.path(outdir, "goi_expression.csv"))

if (length(markers_present) > 0) {
  summary_df$pct_reads_in_markers <-
    round(100 * colSums(counts[markers_present, , drop = FALSE]) /
          colSums(counts), 4)[as.character(summary_df$sample)]
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
