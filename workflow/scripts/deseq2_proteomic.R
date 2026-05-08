#!/usr/bin/env Rscript
# Proteomics-filtered volcano plot for BulkRNAPipe

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(dplyr)
  library(readxl)
})

option_list <- list(
  make_option("--results",  type = "character", help = "DESeq2 results.csv file"),
  make_option("--outdir",   type = "character", help = "Output directory"),
  make_option("--contrast", type = "character",
              help = "Numerator and denominator conditions space-separated (e.g. 'treatment control')"),
  make_option("--contrast_name", type = "character", default = "",
              help = "Human-readable name for this contrast"),
  make_option("--padj",     type = "double", default = 0.05,
              help = "Adjusted p-value threshold [default %default]"),
  make_option("--lfc",      type = "double", default = 1.0,
              help = "Absolute log2 fold-change threshold [default %default]"),
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

if (is.null(args$contrast) && length(opt$args) >= 2) {
  contrast_num <- opt$args[1]
  contrast_den <- opt$args[2]
} else {
  vals <- strsplit(args$contrast, " ")[[1]]
  contrast_num <- vals[1]
  contrast_den <- vals[2]
}

contrast_name <- if (!is.null(args$contrast_name) && nchar(trimws(args$contrast_name)) > 0) {
  trimws(args$contrast_name)
} else {
  paste0(contrast_num, "_vs_", contrast_den)
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  "proteomics_xlsx", "proteomics_gene_column",
  "proteomics_fdr_column", "proteomics_logfc_column", "proteomics_comparison"
)
missing_required <- required_inputs[
  vapply(required_inputs, function(nm) nchar(trimws(as.character(args[[nm]]))) == 0, logical(1))
]
if (length(missing_required) > 0) {
  stop("Missing required proteomics parameters: ", paste(missing_required, collapse = ", "))
}
if (!file.exists(args$proteomics_xlsx)) {
  stop("Proteomics xlsx file does not exist: ", args$proteomics_xlsx)
}

parse_proteomics_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

res_df <- read.csv(args$results, stringsAsFactors = FALSE, check.names = FALSE) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    gene_id = trimws(as.character(gene_id)),
    rna_significant = padj < args$padj & abs(log2FoldChange) >= args$lfc,
    rna_direction = dplyr::case_when(
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      TRUE               ~ NA_character_
    )
  )

prot_tbl <- readxl::read_excel(args$proteomics_xlsx, sheet = args$proteomics_sheet, col_types = "text")

wide_format <- nchar(trimws(args$proteomics_comparison_column)) == 0
if (wide_format) {
  actual_logfc_col <- paste0(args$proteomics_comparison, args$proteomics_logfc_column)
  actual_fdr_col   <- paste0(args$proteomics_comparison, args$proteomics_fdr_column)
  req_cols <- c(args$proteomics_gene_column, actual_logfc_col, actual_fdr_col)
} else {
  actual_logfc_col <- args$proteomics_logfc_column
  actual_fdr_col   <- args$proteomics_fdr_column
  req_cols <- c(args$proteomics_gene_column, args$proteomics_comparison_column, actual_fdr_col, actual_logfc_col)
}

missing_cols <- setdiff(req_cols, colnames(prot_tbl))
if (length(missing_cols) > 0) {
  stop("Missing proteomics column(s): ", paste(missing_cols, collapse = ", "))
}

if (wide_format) {
  prot_sig <- prot_tbl %>%
    dplyr::mutate(
      prot_gene  = trimws(as.character(.data[[args$proteomics_gene_column]])),
      prot_fdr   = parse_proteomics_numeric(.data[[actual_fdr_col]]),
      prot_logfc = parse_proteomics_numeric(.data[[actual_logfc_col]])
    )
} else {
  prot_sig <- prot_tbl %>%
    dplyr::filter(trimws(as.character(.data[[args$proteomics_comparison_column]])) == trimws(args$proteomics_comparison)) %>%
    dplyr::mutate(
      prot_gene  = trimws(as.character(.data[[args$proteomics_gene_column]])),
      prot_fdr   = parse_proteomics_numeric(.data[[actual_fdr_col]]),
      prot_logfc = parse_proteomics_numeric(.data[[actual_logfc_col]])
    )
}

prot_sig <- prot_sig %>%
  dplyr::filter(!is.na(prot_gene), prot_gene != "", !is.na(prot_fdr), !is.na(prot_logfc)) %>%
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

# Strict filtering: only genes present in significant proteomics set are kept.
volcano_df <- res_df %>%
  dplyr::semi_join(prot_sig, by = "gene_id") %>%
  dplyr::left_join(prot_sig, by = "gene_id") %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      rna_significant & !is.na(prot_direction) & rna_direction == prot_direction ~ "Significant same direction",
      rna_significant & !is.na(prot_direction) & rna_direction != prot_direction ~ "Significant opposite direction",
      TRUE                                                                        ~ "Not significant"
    ),
    label = ifelse(rna_significant, gene_id, NA_character_)
  )

if (nrow(volcano_df) == 0) {
  warning("No overlap between RNA DESeq2 results and significant proteomics genes for contrast: ", contrast_name)
}

p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), colour = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Significant same direction" = "#33A02C",
    "Significant opposite direction" = "#E31A1C",
    "Not significant" = "grey60"
  )) +
  geom_vline(xintercept = c(-args$lfc, args$lfc), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(args$padj), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title = paste("Volcano (Proteomics):", contrast_num, "vs", contrast_den),
    subtitle = paste0(
      "Proteomics filtered (", args$proteomics_comparison, ", FDR ≤ ", args$proteomics_fdr_threshold, "). ",
      "log2FC > 0: higher in ", contrast_num, " | log2FC < 0: higher in ", contrast_den
    ),
    x = paste0("log2FC (", contrast_num, " / ", contrast_den, ")"),
    y = expression(-log[10]~"adjusted p-value"),
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(args$outdir, "volcano_proteomic.pdf"), plot = p, width = 7, height = 6)

message("Proteomics volcano complete. Output: ", file.path(args$outdir, "volcano_proteomic.pdf"))
