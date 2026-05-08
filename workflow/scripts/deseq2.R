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
  library(readxl)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",   type = "character", help = "featureCounts output file"),
  make_option("--outdir",   type = "character", help = "Output directory"),
  make_option("--contrast", type = "character",
              help = "Numerator and denominator conditions space-separated (e.g. 'treatment control')"),
  make_option("--samples",  type = "character",
              help = "Comma-separated sample:condition pairs"),
  make_option("--padj",          type = "double",    default = 0.05,
              help = "Adjusted p-value threshold [default %default]"),
  make_option("--lfc",           type = "double",    default = 1.0,
              help = "Log2 fold-change threshold [default %default]"),
  make_option("--contrast_name", type = "character", default = "",
              help = "Human-readable name for this contrast (used in metadata file)"),
  make_option("--proteomics_xlsx", type = "character", default = "",
              help = "Path to limma proteomics xlsx file"),
  make_option("--proteomics_sheet", type = "character", default = "limma result",
              help = "Sheet name in limma proteomics xlsx [default %default]"),
  make_option("--proteomics_gene_column", type = "character", default = "",
              help = "Column name for gene IDs in proteomics sheet"),
  make_option("--proteomics_comparison_column", type = "character", default = "",
              help = "Column name for comparison IDs in proteomics sheet"),
  make_option("--proteomics_fdr_column", type = "character", default = "",
              help = "Column name for FDR in proteomics sheet"),
  make_option("--proteomics_logfc_column", type = "character", default = "",
              help = "Column name for logFC in proteomics sheet"),
  make_option("--proteomics_fdr_threshold", type = "double", default = 0.05,
              help = "FDR threshold for proteomics significance [default %default]"),
  make_option("--proteomics_comparison", type = "character", default = "",
              help = "Proteomics comparison mapped to this DESeq2 contrast")
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

outdir        <- args$outdir
padj_thr      <- args$padj
lfc_thr       <- args$lfc
contrast_name <- if (!is.null(args$contrast_name) && nchar(trimws(args$contrast_name)) > 0)
                   trimws(args$contrast_name)
                 else
                   paste0(contrast_num, "_vs_", contrast_den)

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
pdf(file.path(outdir, "ma_plot.pdf"), width = 6, height = 5.5)
plotMA(res, alpha = padj_thr,
       main = paste("MA plot:", contrast_num, "vs", contrast_den),
       sub  = paste0("log2FC > 0: higher in ", contrast_num,
                     ";  log2FC < 0: higher in ", contrast_den))
dev.off()

# ─── Volcano plot helper ─────────────────────────────────────────────────────
make_volcano_plot <- function(df, colors, title, subtitle, lfc_thr, padj_thr,
                              contrast_num, contrast_den,
                              label_mode = c("none", "top", "all"),
                              label_top_n = 10) {
  label_mode <- match.arg(label_mode)
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                      colour = significance)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = colors) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = paste0("log2FC (", contrast_num, " / ", contrast_den, ")"),
      y        = expression(-log[10]~"adjusted p-value"),
      colour   = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle   = element_text(size = 9, colour = "grey30"))

  label_df <- df %>%
    dplyr::filter(!is.na(label))

  if (label_mode == "top" && nrow(label_df) > 0) {
    label_df <- label_df %>%
      dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
      dplyr::slice_head(n = label_top_n)
  }

  if (label_mode != "none" && nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 2.5,
      max.overlaps = nrow(label_df),
      show.legend = FALSE
    )
  }
  p
}

# ─── Volcano plot ────────────────────────────────────────────────────────────
volcano_df_base <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    gene_id = trimws(gene_id),
    rna_significant = padj < padj_thr & abs(log2FoldChange) >= lfc_thr,
    rna_direction = dplyr::case_when(
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      TRUE               ~ NA_character_
    )
  )

use_proteomics <- nchar(trimws(args$proteomics_xlsx)) > 0 &&
                  file.exists(args$proteomics_xlsx) &&
                  nchar(trimws(args$proteomics_comparison)) > 0 &&
                  nchar(trimws(args$proteomics_gene_column)) > 0 &&
                  nchar(trimws(args$proteomics_fdr_column)) > 0 &&
                  nchar(trimws(args$proteomics_logfc_column)) > 0
direction_subtitle <- paste0("log2FC > 0: higher in ", contrast_num,
                             "   \u2502   log2FC < 0: higher in ", contrast_den)

# Build the full RNA volcano data frame (used in both modes)
volcano_df_rna <- volcano_df_base %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      rna_significant ~ "Significant",
      TRUE            ~ "Not significant"
    ),
    label = ifelse(rna_significant, gene_id, NA_character_)
  )
rna_colors <- c("Significant" = "#E41A1C", "Not significant" = "grey60")

if (use_proteomics) {
  # Read all columns as text to prevent Excel date/numeric coercions of gene
  # symbols (e.g. "MARCH1" → date, "SEPT7" → date).  Numeric columns are
  # converted explicitly with as.numeric() further below.
  prot_tbl <- readxl::read_excel(args$proteomics_xlsx, sheet = args$proteomics_sheet,
                                  col_types = "text")

  # Detect wide vs. long format.
  # Wide format: comparison_column is empty; logfc_column and fdr_column are
  #   suffixes that get prepended with the comparison name to form the actual
  #   column names (e.g. comparison="S1011_vs_control", logfc_column="_LOG2FC"
  #   → actual column "S1011_vs_control_LOG2FC").
  # Long format: comparison_column is non-empty; rows are filtered by
  #   comparison_column == proteomics_comparison and logfc_column / fdr_column
  #   are used as-is.
  wide_format <- nchar(trimws(args$proteomics_comparison_column)) == 0
  if (wide_format) {
    # logfc_column and fdr_column are suffixes (must include any delimiter,
    # e.g. "_LOG2FC" and "_adj.P.Val") that are appended to the comparison
    # name to form the actual Excel column names.
    actual_logfc_col <- paste0(args$proteomics_comparison, args$proteomics_logfc_column)
    actual_fdr_col   <- paste0(args$proteomics_comparison, args$proteomics_fdr_column)
    req_cols <- c(args$proteomics_gene_column, actual_logfc_col, actual_fdr_col)
  } else {
    actual_logfc_col <- args$proteomics_logfc_column
    actual_fdr_col   <- args$proteomics_fdr_column
    req_cols <- c(
      args$proteomics_gene_column,
      args$proteomics_comparison_column,
      actual_fdr_col,
      actual_logfc_col
    )
  }
  # Validate that all required columns (including any constructed wide-format
  # column names) are present in the Excel sheet before proceeding.
  missing_cols <- setdiff(req_cols, colnames(prot_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing proteomics column(s): ", paste(missing_cols, collapse = ", "))
  }

  if (wide_format) {
    prot_sig <- prot_tbl %>%
      dplyr::mutate(
        prot_gene  = trimws(as.character(.data[[args$proteomics_gene_column]])),
        prot_fdr   = suppressWarnings(as.numeric(.data[[actual_fdr_col]])),
        prot_logfc = suppressWarnings(as.numeric(.data[[actual_logfc_col]]))
      )
  } else {
    prot_sig <- prot_tbl %>%
      dplyr::filter(.data[[args$proteomics_comparison_column]] == args$proteomics_comparison) %>%
      dplyr::mutate(
        prot_gene  = trimws(as.character(.data[[args$proteomics_gene_column]])),
        prot_fdr   = suppressWarnings(as.numeric(.data[[actual_fdr_col]])),
        prot_logfc = suppressWarnings(as.numeric(.data[[actual_logfc_col]]))
      )
  }
  prot_sig <- prot_sig %>%
    dplyr::filter(!is.na(prot_gene), trimws(prot_gene) != "", !is.na(prot_fdr), !is.na(prot_logfc)) %>%
    dplyr::filter(prot_fdr <= args$proteomics_fdr_threshold) %>%
    dplyr::arrange(prot_fdr, dplyr::desc(abs(prot_logfc))) %>%
    dplyr::distinct(prot_gene, .keep_all = TRUE) %>%
    dplyr::mutate(
      prot_direction = dplyr::case_when(
        prot_logfc > 0 ~ "Up",
        prot_logfc < 0 ~ "Down",
        TRUE           ~ NA_character_
      )
    ) %>%
    dplyr::select(gene_id = prot_gene, prot_direction)

  volcano_df_prot <- volcano_df_base %>%
    dplyr::inner_join(prot_sig, by = "gene_id") %>%
    dplyr::mutate(
      concordance_eligible = rna_significant & !is.na(prot_direction),
      significance = dplyr::case_when(
        concordance_eligible & rna_direction == prot_direction ~ "Significant same direction",
        concordance_eligible & rna_direction != prot_direction ~ "Significant opposite direction",
        TRUE                                                   ~ "Not significant"
      ),
      label = ifelse(rna_significant, gene_id, NA_character_)
    ) %>%
    dplyr::select(-concordance_eligible)

  # ── Diagnostic logging ────────────────────────────────────────────────────
  n_rna_genes  <- nrow(volcano_df_base)
  n_prot_sig   <- nrow(prot_sig)
  n_overlap    <- nrow(volcano_df_prot)
  message("Proteomics join diagnostics:")
  message("  RNA genes with non-NA padj  : ", n_rna_genes)
  message("  Proteomics significant genes: ", n_prot_sig,
          " (FDR <= ", args$proteomics_fdr_threshold, ")")
  message("  Overlap (inner join)        : ", n_overlap)
  if (n_prot_sig > 0 && n_overlap < n_prot_sig) {
    unmatched <- setdiff(prot_sig$gene_id, volcano_df_base$gene_id)
    message("  Proteomics genes not in RNA data (", length(unmatched), "): ",
            paste(head(unmatched, 10), collapse = ", "),
            if (length(unmatched) > 10) " ..." else "")
  }

  if (nrow(volcano_df_prot) == 0) {
    warning("No overlap between DESeq2 genes and significant proteomics genes for contrast: ",
            contrast_name,
            " (proteomics comparison: ", args$proteomics_comparison, ")")
  }

  prot_colors <- c(
    "Significant same direction"     = "#33A02C",
    "Significant opposite direction" = "#E31A1C",
    "Not significant"                = "grey60"
  )
  prot_subtitle <- paste0(
    "Filtered to significant proteomics genes from ",
    args$proteomics_comparison,
    " (FDR \u2264 ", args$proteomics_fdr_threshold, "); ",
    direction_subtitle
  )

  # Plot 1: full RNA volcano (saved as volcano.pdf for backwards compatibility)
  # Keep a minimal set of top RNA-significant labels.
  p_rna <- make_volcano_plot(
    df           = volcano_df_rna,
    colors       = rna_colors,
    title        = paste("Volcano (RNA):", contrast_num, "vs", contrast_den),
    subtitle     = direction_subtitle,
    lfc_thr      = lfc_thr,
    padj_thr     = padj_thr,
    contrast_num = contrast_num,
    contrast_den = contrast_den,
    label_mode   = "top"
  )
  ggsave(file.path(outdir, "volcano.pdf"), plot = p_rna, width = 7, height = 6)

  # Plot 2: proteomics-filtered concordance volcano (all significant labels)
  p_prot <- make_volcano_plot(
    df           = volcano_df_prot,
    colors       = prot_colors,
    title        = paste("Volcano (Proteomics):", contrast_num, "vs", contrast_den),
    subtitle     = prot_subtitle,
    lfc_thr      = lfc_thr,
    padj_thr     = padj_thr,
    contrast_num = contrast_num,
    contrast_den = contrast_den,
    label_mode   = "all"
  )
  ggsave(file.path(outdir, "volcano_proteomics.pdf"), plot = p_prot, width = 7, height = 6)

} else {
  p <- make_volcano_plot(
    df           = volcano_df_rna,
    colors       = rna_colors,
    title        = paste("Volcano:", contrast_num, "vs", contrast_den),
    subtitle     = direction_subtitle,
    lfc_thr      = lfc_thr,
    padj_thr     = padj_thr,
    contrast_num = contrast_num,
    contrast_den = contrast_den,
    label_mode   = "top"
  )
  ggsave(file.path(outdir, "volcano.pdf"), plot = p, width = 7, height = 6)
}

# ─── Contrast metadata file ───────────────────────────────────────────────────
writeLines(
  c(
    paste0("contrast_name: \"", contrast_name, "\""),
    paste0("numerator: \"",     contrast_num,  "\""),
    paste0("denominator: \"",   contrast_den,  "\""),
    paste0("direction_note: \"log2FC > 0 is higher in ", contrast_num,
           "; log2FC < 0 is higher in ", contrast_den, "\"")
  ),
  file.path(outdir, "contrast_info.yaml")
)

message("DESeq2 analysis complete. Results written to: ", outdir)
