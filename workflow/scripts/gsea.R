#!/usr/bin/env Rscript
# gsea.R — fgsea (fgseaMultilevel) enrichment for BulkRNAPipe
#
# Inputs:
#   --results      : DESeq2 results.csv (gene_id, stat, log2FoldChange columns)
#   --hox_gmt      : auto-generated HOX GMT file
#   --outdir       : output directory
#   --collection   : gene-set collection slug (e.g. H, C2_CP_REACTOME, C5_GO_BP)
#   --rank_metric  : "stat" (default) or "log2FoldChange"
#   --min_size     : minimum gene-set size
#   --max_size     : maximum gene-set size
#   --nperm        : nPermSimple for fgseaMultilevel (simple Monte-Carlo estimate used
#                    to seed the multilevel algorithm; default 1000)
#   --custom_gmts  : comma-separated paths to additional GMT files (optional)

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(tibble)
  library(stringr)
})

# ── Helper: NULL-coalescing operator ─────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Parse arguments ────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--results",     type = "character", help = "DESeq2 results CSV"),
  make_option("--hox_gmt",     type = "character", help = "HOX GMT file"),
  make_option("--outdir",      type = "character", help = "Output directory"),
  make_option("--collection",  type = "character", help = "Collection slug (e.g. H, C2_CP_REACTOME)"),
  make_option("--rank_metric", type = "character", default = "stat",
              help = "Ranking metric: stat or log2FoldChange [default %default]"),
  make_option("--min_size",    type = "integer",   default = 15L),
  make_option("--max_size",    type = "integer",   default = 500L),
  make_option("--nperm",       type = "integer",   default = 1000L),
  make_option("--custom_gmts", type = "character", default = "",
              help = "Comma-separated paths to extra GMT files"),
  make_option("--contrast_name", type = "character", default = "",
              help = "Name of this contrast (e.g. ATRX_IFF_vs_ATRX_FL)"),
  make_option("--numerator",     type = "character", default = "",
              help = "Numerator condition of the DESeq2 contrast"),
  make_option("--denominator",   type = "character", default = "",
              help = "Denominator condition of the DESeq2 contrast")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load DESeq2 results ────────────────────────────────────────────────────────
res <- read.csv(args$results, stringsAsFactors = FALSE)

# Choose ranking metric
if (args$rank_metric == "stat" &&
    "stat" %in% colnames(res) &&
    !all(is.na(res$stat))) {
  rank_col <- "stat"
} else {
  rank_col <- "log2FoldChange"
  message("NOTE: 'stat' not available or all NA; falling back to log2FoldChange.")
}

# Build named rank vector, removing NAs and duplicate gene IDs
rank_df <- res[, c("gene_id", rank_col), drop = FALSE]
colnames(rank_df)[2] <- "score"
rank_df <- rank_df[!is.na(rank_df$score), ]
rank_df <- rank_df[!duplicated(rank_df$gene_id), ]
ranks   <- setNames(rank_df$score, rank_df$gene_id)
ranks   <- sort(ranks, decreasing = TRUE)

# ── Load gene sets ─────────────────────────────────────────────────────────────
# Parse collection slug back to category / subcategory
# "H" → category="H", subcategory=NULL
# "C2_CP_REACTOME" → category="C2", subcategory="CP:REACTOME"
parse_collection_slug <- function(slug) {
  parts <- str_split(slug, "_", n = 2)[[1]]
  if (length(parts) == 1) {
    list(category = parts[1], subcategory = NULL)
  } else {
    list(category = parts[1], subcategory = str_replace_all(parts[2], "_", ":"))
  }
}

coll <- parse_collection_slug(args$collection)

message(sprintf("Fetching msigdbr sets: category=%s subcategory=%s",
                coll$category, coll$subcategory %||% ""))

msig_df <- if (!is.null(coll$subcategory)) {
  msigdbr(species = "Homo sapiens",
          category    = coll$category,
          subcategory = coll$subcategory)
} else {
  msigdbr(species = "Homo sapiens",
          category = coll$category)
}

# Convert to named list of gene vectors
pathways_msig <- split(msig_df$gene_symbol, msig_df$gs_name)

# Load HOX GMT
pathways_hox <- gmtPathways(args$hox_gmt)

# Load optional custom GMTs
pathways_custom <- list()
if (nchar(trimws(args$custom_gmts)) > 0) {
  gmt_paths <- str_split(trimws(args$custom_gmts), ",")[[1]]
  gmt_paths <- gmt_paths[nchar(trimws(gmt_paths)) > 0]
  for (gp in gmt_paths) {
    if (file.exists(gp)) {
      pg <- gmtPathways(gp)
      pathways_custom <- c(pathways_custom, pg)
      message(sprintf("Loaded %d gene sets from %s", length(pg), gp))
    } else {
      warning(sprintf("Custom GMT file not found, skipping: %s", gp))
    }
  }
}

# Merge all sets (msigdbr canonical → HOX → custom); custom sets win on name clash
pathways <- c(pathways_msig, pathways_hox, pathways_custom)

message(sprintf("Total gene sets after merge: %d", length(pathways)))

# ── Run fgseaMultilevel ────────────────────────────────────────────────────────
set.seed(42)
fgsea_res <- fgseaMultilevel(
  pathways   = pathways,
  stats      = ranks,
  minSize    = args$min_size,
  maxSize    = args$max_size,
  nPermSimple = args$nperm
)

# Convert leadingEdge list column to a semicolon-separated string for CSV
fgsea_dt <- as.data.table(fgsea_res)
fgsea_dt[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
fgsea_dt <- fgsea_dt[order(padj, -abs(NES))]

# Annotate with contrast / directionality metadata
fgsea_dt[, contrast_name  := args$contrast_name]
fgsea_dt[, numerator      := args$numerator]
fgsea_dt[, denominator    := args$denominator]
fgsea_dt[, rank_metric    := rank_col]
fgsea_dt[, direction_note := paste0("NES > 0: enriched in ", args$numerator,
                                    "; NES < 0: enriched in ", args$denominator)]

# Save results
out_csv <- file.path(args$outdir, paste0(args$collection, "_results.csv"))
write.csv(fgsea_dt, out_csv, row.names = FALSE, quote = FALSE)
message("Results written to: ", out_csv)

# ── Dotplot ────────────────────────────────────────────────────────────────────
plot_df <- fgsea_dt[!is.na(padj)][order(NES, decreasing = TRUE)]

# Top 15 up + top 15 down
top_up   <- head(plot_df[NES > 0], 15)
top_down <- head(plot_df[NES < 0][order(NES)], 15)
top_paths <- rbind(top_up, top_down)

if (nrow(top_paths) == 0) {
  message("No pathways passed filters; skipping dotplot.")
  pdf(file.path(args$outdir, paste0(args$collection, "_dotplot.pdf")), width = 8, height = 1)
  plot.new()
  text(0.5, 0.5, "No enriched pathways")
  dev.off()
} else {
  top_paths[, pathway_label := str_trunc(pathway, 55)]

  # Direction labels that name the conditions
  lbl_up   <- paste0("Enriched in ", args$numerator,   " (NES \u003e 0)")
  lbl_down <- paste0("Enriched in ", args$denominator, " (NES \u003c 0)")
  top_paths[, direction := factor(
    ifelse(NES > 0, lbl_up, lbl_down),
    levels = c(lbl_up, lbl_down)
  )]

  # Sort so the most-enriched pathway in each direction appears at the top of
  # its panel.  Use |NES| as the ordering key so that within each panel the
  # highest absolute enrichment is topmost.
  top_paths[, plot_order := abs(NES)]

  p <- ggplot(top_paths,
              aes(x = NES, y = reorder(pathway_label, plot_order),
                  size = size, colour = padj)) +
    geom_point() +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    scale_colour_gradient(low = "#E41A1C", high = "grey70", limits = c(0, 0.25),
                          oob = scales::squish, name = "padj") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    labs(
      title    = paste("GSEA:", args$collection),
      subtitle = paste0("Ranking metric: ", rank_col, " from DESeq2 (",
                        args$numerator, " vs ", args$denominator, ")"),
      x = "Normalised enrichment score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text       = element_text(face = "bold", size = 10),
          plot.subtitle    = element_text(size = 9, colour = "grey30")) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")

  n_up   <- nrow(top_up)
  n_down <- nrow(top_down)
  out_pdf <- file.path(args$outdir, paste0(args$collection, "_dotplot.pdf"))
  ggsave(out_pdf, plot = p,
         width  = 9,
         height = max(5, (n_up + n_down) * 0.35 + 3.5))
  message("Dotplot written to: ", out_pdf)
}

message("GSEA complete for collection: ", args$collection)
