#!/usr/bin/env Rscript
# go.R — GO enrichment analysis (clusterProfiler) for BulkRNAPipe
#
# Inputs:
#   --results     : DESeq2 results.csv (gene_id, padj, log2FoldChange columns)
#   --outdir      : output directory
#   --ontology    : GO ontology: BP, MF, or CC
#   --padj_cutoff : p.adjust cutoff for enrichGO output
#   --min_gs_size : minimum gene-set size for enrichGO
#   --max_gs_size : maximum gene-set size for enrichGO
#   --padj_thr    : DESeq2 padj threshold to call significant genes
#   --lfc_thr     : DESeq2 |log2FC| threshold to call significant genes

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
  make_option("--padj_cutoff", type = "double",   default = 0.05,
              help = "clusterProfiler p.adjust cutoff [default %default]"),
  make_option("--min_gs_size", type = "integer",  default = 10L,
              help = "Minimum gene-set size [default %default]"),
  make_option("--max_gs_size", type = "integer",  default = 500L,
              help = "Maximum gene-set size [default %default]"),
  make_option("--padj_thr",    type = "double",   default = 0.05,
              help = "DESeq2 padj threshold [default %default]"),
  make_option("--lfc_thr",     type = "double",   default = 1.0,
              help = "DESeq2 |log2FC| threshold [default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

ont <- toupper(args$ontology)

# ── Load DESeq2 results ────────────────────────────────────────────────────────
res <- read.csv(args$results, stringsAsFactors = FALSE)

# Require gene_id, padj, log2FoldChange
required_cols <- c("gene_id", "padj", "log2FoldChange")
missing_cols  <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns in results CSV: ", paste(missing_cols, collapse = ", "))
}

# ── Select significant genes ───────────────────────────────────────────────────
sig_genes <- res %>%
  filter(!is.na(padj),
         padj < args$padj_thr,
         abs(log2FoldChange) >= args$lfc_thr) %>%
  pull(gene_id)

# Universe = all tested genes (with non-NA padj)
universe_genes <- res %>%
  filter(!is.na(padj)) %>%
  pull(gene_id)

message(sprintf("Significant genes: %d / Universe: %d",
                length(sig_genes), length(universe_genes)))

if (length(sig_genes) < 5) {
  warning("Fewer than 5 significant genes; GO enrichment may be uninformative.")
}

# ── Map SYMBOL → ENTREZID ─────────────────────────────────────────────────────
map_to_entrez <- function(symbols) {
  mapped <- mapIds(org.Hs.eg.db,
                   keys     = symbols,
                   column   = "ENTREZID",
                   keytype  = "SYMBOL",
                   multiVals = "first")
  mapped <- mapped[!is.na(mapped)]
  unique(mapped)
}

sig_entrez     <- map_to_entrez(sig_genes)
universe_entrez <- map_to_entrez(universe_genes)

message(sprintf("Mapped to ENTREZID — significant: %d / universe: %d",
                length(sig_entrez), length(universe_entrez)))

# ── Run enrichGO ─────────────────────────────────────────────────────────────
ego <- enrichGO(
  gene          = sig_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = ont,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,       # keep all; filter by padj_cutoff later
  qvalueCutoff  = 1,
  minGSSize     = args$min_gs_size,
  maxGSSize     = args$max_gs_size,
  readable      = TRUE     # map ENTREZID back to SYMBOL in output
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  message("No GO terms enriched.")
  out_csv <- file.path(args$outdir, paste0("go_", tolower(ont), "_results.csv"))
  write.csv(data.frame(), out_csv, row.names = FALSE)
  pdf(file.path(args$outdir, paste0("go_", tolower(ont), "_dotplot.pdf")),
      width = 8, height = 1)
  plot.new()
  text(0.5, 0.5, "No enriched GO terms")
  dev.off()
  message("Empty results written.")
  quit(status = 0)
}

# Apply padj cutoff
ego_df <- as.data.frame(ego) %>%
  filter(p.adjust <= args$padj_cutoff) %>%
  arrange(p.adjust)

# Save results
out_csv <- file.path(args$outdir, paste0("go_", tolower(ont), "_results.csv"))
write.csv(ego_df, out_csv, row.names = FALSE, quote = FALSE)
message("GO results written to: ", out_csv)

# ── Dotplot ────────────────────────────────────────────────────────────────────
top_terms <- head(ego_df, 30)

if (nrow(top_terms) == 0) {
  pdf(file.path(args$outdir, paste0("go_", tolower(ont), "_dotplot.pdf")),
      width = 8, height = 1)
  plot.new()
  text(0.5, 0.5, "No terms passed padj cutoff")
  dev.off()
} else {
  top_terms <- top_terms %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Description_short = substr(Description, 1, 55)
    )

  p <- ggplot(top_terms,
              aes(x = GeneRatio_num,
                  y = reorder(Description_short, GeneRatio_num),
                  size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#E41A1C", high = "grey70", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(title = paste("GO enrichment:", ont),
         x = "Gene ratio",
         y = NULL) +
    theme_bw(base_size = 11)

  out_pdf <- file.path(args$outdir, paste0("go_", tolower(ont), "_dotplot.pdf"))
  ggsave(out_pdf, plot = p,
         width  = 9,
         height = max(4, nrow(top_terms) * 0.35 + 2))
  message("Dotplot written to: ", out_pdf)
}

message("GO enrichment complete for ontology: ", ont)
