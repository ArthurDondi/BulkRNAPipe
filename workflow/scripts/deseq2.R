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

# ─── Volcano plot ────────────────────────────────────────────────────────────
volcano_df_base <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
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
                  nchar(trimws(args$proteomics_comparison_column)) > 0 &&
                  nchar(trimws(args$proteomics_fdr_column)) > 0 &&
                  nchar(trimws(args$proteomics_logfc_column)) > 0

if (use_proteomics) {
  prot_tbl <- readxl::read_excel(args$proteomics_xlsx, sheet = args$proteomics_sheet)
  req_cols <- c(
    args$proteomics_gene_column,
    args$proteomics_comparison_column,
    args$proteomics_fdr_column,
    args$proteomics_logfc_column
  )
  missing_cols <- setdiff(req_cols, colnames(prot_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing proteomics column(s): ", paste(missing_cols, collapse = ", "))
  }

  prot_sig <- prot_tbl %>%
    dplyr::filter(.data[[args$proteomics_comparison_column]] == args$proteomics_comparison) %>%
    dplyr::mutate(
      prot_gene = as.character(.data[[args$proteomics_gene_column]]),
      prot_fdr = suppressWarnings(as.numeric(.data[[args$proteomics_fdr_column]])),
      prot_logfc = suppressWarnings(as.numeric(.data[[args$proteomics_logfc_column]]))
    ) %>%
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

  volcano_df <- volcano_df_base %>%
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

  if (nrow(volcano_df) == 0) {
    warning("No overlap between DESeq2 genes and significant proteomics genes for contrast: ",
            contrast_name,
            " (proteomics comparison: ", args$proteomics_comparison, ")")
  }

  volcano_colors <- c(
    "Significant same direction"     = "#33A02C",
    "Significant opposite direction" = "#E31A1C",
    "Not significant"                = "grey60"
  )
  volcano_subtitle <- paste0(
    "Filtered to significant proteomics genes from ",
    args$proteomics_comparison,
    " (FDR ≤ ", args$proteomics_fdr_threshold, "); log2FC > 0: higher in ", contrast_num,
    "   \u2502   log2FC < 0: higher in ", contrast_den
  )
} else {
  volcano_df <- volcano_df_base %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        rna_significant ~ "Significant",
        TRUE            ~ "Not significant"
      ),
      label = ifelse(rna_significant, gene_id, NA_character_)
    )
  volcano_colors <- c("Significant" = "#E41A1C",
                      "Not significant" = "grey60")
  volcano_subtitle <- paste0("log2FC > 0: higher in ", contrast_num,
                             "   \u2502   log2FC < 0: higher in ", contrast_den)
}

p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj),
                             colour = significance, label = label)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(size = 2.5, max.overlaps = sum(volcano_df$label != "", na.rm = TRUE), show.legend = FALSE) +
  scale_colour_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title    = paste("Volcano:", contrast_num, "vs", contrast_den),
    subtitle = volcano_subtitle,
    x        = paste0("log2FC (", contrast_num, " / ", contrast_den, ")"),
    y        = expression(-log[10]~"adjusted p-value"),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle   = element_text(size = 9, colour = "grey30"))

ggsave(file.path(outdir, "volcano.pdf"), plot = p, width = 7, height = 6)

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
