#!/usr/bin/env Rscript
# go.R — GO enrichment analysis (clusterProfiler) for BulkRNAPipe
#
# Inputs:
#   --results     : DESeq2 results.csv (gene_id, padj, log2FoldChange columns)
#   --outdir      : output directory
#   --ontology    : GO ontology: BP, MF, or CC
#   --direction   : expression direction: up or down
#   --padj_cutoff : p.adjust cutoff for enrichGO output
#   --min_gs_size : minimum gene-set size for enrichGO
#   --max_gs_size : maximum gene-set size for enrichGO
#   --padj_thr    : DESeq2 padj threshold to call significant genes
#   --lfc_thr     : DESeq2 log2FC threshold to call significant genes
#   --gene_id_type: input gene ID type: SYMBOL, ENSEMBL, or ENTREZID

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(optparse)
  library(dplyr)
  library(tibble)
})

# ── Parse arguments ────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--results",     type = "character", help = "DESeq2 results CSV"),
  make_option("--outdir",      type = "character", help = "Output directory"),
  make_option("--ontology",    type = "character", default = "BP",
              help = "GO ontology: BP, MF, or CC [default %default]"),
  make_option("--direction",   type = "character", default = "up",
              help = "Direction: up or down [default %default]"),
  make_option("--padj_cutoff", type = "double",   default = 0.05,
              help = "clusterProfiler p.adjust cutoff [default %default]"),
  make_option("--min_gs_size", type = "integer",  default = 10L,
              help = "Minimum gene-set size [default %default]"),
  make_option("--max_gs_size", type = "integer",  default = 500L,
              help = "Maximum gene-set size [default %default]"),
  make_option("--padj_thr",    type = "double",   default = 0.05,
              help = "DESeq2 padj threshold [default %default]"),
  make_option("--lfc_thr",     type = "double",   default = 1.0,
              help = "DESeq2 absolute log2FC threshold applied directionally [default %default]"),
  make_option("--gene_id_type", type = "character", default = "SYMBOL",
              help = "Input gene ID type: SYMBOL, ENSEMBL, ENTREZID [default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

ont <- toupper(args$ontology)
dir_label <- tolower(args$direction)
gene_id_type <- toupper(args$gene_id_type)
allowed_ont <- c("BP", "MF", "CC")
allowed_dir <- c("up", "down")
allowed_gene_id <- c("SYMBOL", "ENSEMBL", "ENTREZID")

if (!(ont %in% allowed_ont)) {
  stop("Invalid --ontology: ", ont, ". Allowed values: ", paste(allowed_ont, collapse = ", "))
}
if (!(dir_label %in% allowed_dir)) {
  stop("Invalid --direction: ", dir_label, ". Allowed values: ", paste(allowed_dir, collapse = ", "))
}
if (!(gene_id_type %in% allowed_gene_id)) {
  stop("Invalid --gene_id_type: ", gene_id_type, ". Allowed values: ", paste(allowed_gene_id, collapse = ", "))
}

prefix <- paste0("go_", tolower(ont), "_", dir_label)
# Conservative sanity threshold: if fewer than ~20% of background IDs map,
# enrichment is likely dominated by identifier mismatch rather than biology.
MIN_MAPPING_RATE <- 0.20
MIN_MAPPED_SIG_GENES <- 5
# Keep term labels compact for PDF readability.
MAX_DESCRIPTION_LENGTH <- 55
out_csv_raw <- file.path(args$outdir, paste0(prefix, "_results.csv"))
out_csv_simplified <- file.path(args$outdir, paste0(prefix, "_results_simplified.csv"))
out_pdf_raw <- file.path(args$outdir, paste0(prefix, "_dotplot.pdf"))
out_pdf_simplified <- file.path(args$outdir, paste0(prefix, "_dotplot_simplified.pdf"))
out_unmapped_universe <- file.path(args$outdir, paste0(prefix, "_unmapped_universe.csv"))
out_unmapped_sig <- file.path(args$outdir, paste0(prefix, "_unmapped_sig.csv"))

message(sprintf(
  "GO parameters | ontology=%s | direction=%s | gene_id_type=%s | padj_thr=%.4g | lfc_thr=%.4g | padj_cutoff=%.4g | min_gs_size=%d | max_gs_size=%d",
  ont, dir_label, gene_id_type, args$padj_thr, args$lfc_thr, args$padj_cutoff, args$min_gs_size, args$max_gs_size
))

write_empty_csv <- function(path) {
  write.csv(data.frame(), path, row.names = FALSE)
}

write_message_pdf <- function(path, text_label) {
  tryCatch({
    pdf(path, width = 8, height = 1.5)
    on.exit(dev.off(), add = TRUE)
    plot.new()
    text(0.5, 0.5, text_label)
  }, error = function(e) {
    warning("Failed to write PDF '", path, "': ", conditionMessage(e))
  })
}

write_empty_outputs <- function(reason) {
  write_empty_csv(out_csv_raw)
  write_empty_csv(out_csv_simplified)
  write_message_pdf(out_pdf_raw, reason)
  write_message_pdf(out_pdf_simplified, reason)
  message("Empty GO outputs written. Reason: ", reason)
}

make_dotplot <- function(df, title_text, out_pdf) {
  top_terms <- head(df, 30)
  if (nrow(top_terms) == 0) {
    write_message_pdf(out_pdf, "No terms passed padj cutoff")
    return(invisible(NULL))
  }

  top_terms <- top_terms %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Description_short = substr(Description, 1, MAX_DESCRIPTION_LENGTH)
    )

  p <- ggplot(top_terms,
              aes(x = GeneRatio_num,
                  y = reorder(Description_short, GeneRatio_num),
                  size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#E41A1C", high = "grey70", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(title = title_text, x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 11)

  ggsave(out_pdf, plot = p,
         width = 9,
         height = max(4, nrow(top_terms) * 0.35 + 2))
}

# ── Load DESeq2 results ────────────────────────────────────────────────────────
res <- read.csv(args$results, stringsAsFactors = FALSE)

# Require gene_id, padj, log2FoldChange
required_cols <- c("gene_id", "padj", "log2FoldChange")
missing_cols  <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns in results CSV: ", paste(missing_cols, collapse = ", "))
}

universe_genes <- res %>%
  filter(!is.na(padj)) %>%
  pull(gene_id)

if (dir_label == "up") {
  sig_genes <- res %>%
    filter(!is.na(padj),
           padj < args$padj_thr,
           log2FoldChange >= args$lfc_thr) %>%
    pull(gene_id)
} else {
  sig_genes <- res %>%
    filter(!is.na(padj),
           padj < args$padj_thr,
           log2FoldChange <= (-args$lfc_thr)) %>%
    pull(gene_id)
}

message(sprintf(
  "Input counts | universe_genes=%d | significant_genes_%s=%d",
  length(universe_genes), dir_label, length(sig_genes)
))

if (length(sig_genes) == 0) {
  warning("No significant genes for direction '", dir_label, "'.")
  write_empty_csv(out_unmapped_universe)
  write_empty_csv(out_unmapped_sig)
  write_empty_outputs(paste0("No significant ", dir_label, " genes"))
  quit(status = 0)
}

# ── Map input IDs to ENTREZID ──────────────────────────────────────────────────
map_to_entrez <- function(ids, keytype) {
  ids <- trimws(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ids_unique <- unique(ids)

  if (keytype == "ENTREZID") {
    return(list(
      mapped = ids_unique,
      unmapped = character(0),
      input_n = length(ids_unique)
    ))
  }

  mapped <- mapIds(
    org.Hs.eg.db,
    keys = ids_unique,
    column = "ENTREZID",
    keytype = keytype,
    multiVals = "first"
  )
  mapped_vec <- unname(mapped[!is.na(mapped)])
  unmapped <- names(mapped)[is.na(mapped)]

  list(
    mapped = unique(mapped_vec),
    unmapped = unique(unmapped),
    input_n = length(ids_unique)
  )
}

compute_mapping_rate <- function(input_n, mapped_n) {
  if (input_n > 0) {
    return(mapped_n / input_n)
  }
  0
}

sig_map <- map_to_entrez(sig_genes, gene_id_type)
universe_map <- map_to_entrez(universe_genes, gene_id_type)

sig_entrez <- sig_map$mapped
universe_entrez <- universe_map$mapped
sig_in_universe <- intersect(sig_entrez, universe_entrez)

write.csv(data.frame(gene_id = universe_map$unmapped), out_unmapped_universe, row.names = FALSE, quote = FALSE)
write.csv(data.frame(gene_id = sig_map$unmapped), out_unmapped_sig, row.names = FALSE, quote = FALSE)

mapped_universe_rate <- compute_mapping_rate(universe_map$input_n, length(universe_entrez))
mapped_sig_rate <- compute_mapping_rate(sig_map$input_n, length(sig_entrez))

message(sprintf("Mapping diagnostics | universe_input=%d | universe_mapped=%d (%.1f%%)",
                universe_map$input_n, length(universe_entrez), 100 * mapped_universe_rate))
message(sprintf("Mapping diagnostics | significant_input_%s=%d | significant_mapped_%s=%d (%.1f%%)",
                dir_label, sig_map$input_n, dir_label, length(sig_entrez), 100 * mapped_sig_rate))
message(sprintf("Mapping diagnostics | overlap(mapped_sig_%s, mapped_universe)=%d",
                dir_label, length(sig_in_universe)))
message("Unmapped IDs written to: ", out_unmapped_universe, " and ", out_unmapped_sig)

if (mapped_universe_rate < MIN_MAPPING_RATE || length(sig_in_universe) < MIN_MAPPED_SIG_GENES) {
  warning(sprintf(
    "Extremely low mapping detected (universe mapped=%.1f%%, mapped significant overlap=%d). Writing empty GO outputs.",
    100 * mapped_universe_rate, length(sig_in_universe)
  ))
  write_empty_outputs("Mapping too low for reliable enrichment")
  quit(status = 0)
}

# ── Run enrichGO ─────────────────────────────────────────────────────────────
ego <- enrichGO(
  gene          = sig_in_universe,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = ont,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,       # keep all; filter by padj_cutoff later
  qvalueCutoff  = 1,
  minGSSize     = args$min_gs_size,
  maxGSSize     = args$max_gs_size,
  readable      = TRUE
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  message("No GO terms enriched.")
  write_empty_outputs("No enriched GO terms")
  quit(status = 0)
}

raw_df <- as.data.frame(ego) %>%
  filter(p.adjust <= args$padj_cutoff) %>%
  arrange(p.adjust)

write.csv(raw_df, out_csv_raw, row.names = FALSE, quote = FALSE)
message("GO raw results written to: ", out_csv_raw)
make_dotplot(raw_df, paste("GO enrichment:", ont, "|", dir_label), out_pdf_raw)
message("GO raw dotplot written to: ", out_pdf_raw)

# Reduce redundancy by semantic similarity (clusterProfiler::simplify)
# cutoff=0.7, by='p.adjust', select_fun=min are standard conservative defaults.
ego_simplified <- tryCatch(
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
  error = function(e) {
    warning("simplify() failed: ", conditionMessage(e))
    NULL
  }
)

if (is.null(ego_simplified)) {
  write_empty_csv(out_csv_simplified)
  write_message_pdf(out_pdf_simplified, "No simplified GO terms")
  message("GO simplified outputs written as empty (simplify unavailable).")
} else {
  simplified_df <- as.data.frame(ego_simplified) %>%
    filter(p.adjust <= args$padj_cutoff) %>%
    arrange(p.adjust)

  write.csv(simplified_df, out_csv_simplified, row.names = FALSE, quote = FALSE)
  make_dotplot(simplified_df, paste("GO enrichment (simplified):", ont, "|", dir_label), out_pdf_simplified)
  message("GO simplified results written to: ", out_csv_simplified)
  message("GO simplified dotplot written to: ", out_pdf_simplified)
}

message("GO enrichment complete for ontology: ", ont, " | direction: ", dir_label)
