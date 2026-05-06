#!/usr/bin/env Rscript
# gsea_compare.R — comparative GSEA analysis (ΔNES + residual-rank) for BulkRNAPipe
#
# For a pair of DESeq2 contrasts (A and B) this script:
#   1. Computes ΔNES = NES_A - NES_B per pathway (using pre-computed fgsea CSVs)
#   2. Optionally runs fgsea on residual ranks (rank_A - rank_B)
#
# Inputs:
#   --results_a      : DESeq2 results CSV for contrast A
#   --results_b      : DESeq2 results CSV for contrast B
#   --hox_gmt        : HOX GMT file (always included)
#   --gsea_dir_a     : directory with per-collection fgsea CSVs for contrast A
#   --gsea_dir_b     : directory with per-collection fgsea CSVs for contrast B
#   --outdir         : output directory
#   --contrast_a     : name label for contrast A
#   --contrast_b     : name label for contrast B
#   --collections    : comma-separated collection slugs
#   --rank_metric    : "stat" or "log2FoldChange"
#   --min_size / --max_size / --nperm : fgsea parameters (nperm = nPermSimple for
#                                       fgseaMultilevel's simple Monte-Carlo seed)
#   --delta_nes      : TRUE/FALSE — compute ΔNES
#   --residual_rank  : TRUE/FALSE — run fgsea on residual ranks
#   --custom_gmts    : comma-separated paths to extra GMT files

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

# ── Helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_collection_slug <- function(slug) {
  parts <- str_split(slug, "_", n = 2)[[1]]
  if (length(parts) == 1) {
    list(category = parts[1], subcategory = NULL)
  } else {
    list(category = parts[1], subcategory = str_replace_all(parts[2], "_", ":"))
  }
}

load_pathways_for_collection <- function(slug, hox_gmt, custom_gmt_paths) {
  coll   <- parse_collection_slug(slug)
  msig_df <- if (!is.null(coll$subcategory)) {
    msigdbr(species = "Homo sapiens",
            category    = coll$category,
            subcategory = coll$subcategory)
  } else {
    msigdbr(species = "Homo sapiens",
            category = coll$category)
  }
  pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
  pathways <- c(pathways, gmtPathways(hox_gmt))
  for (gp in custom_gmt_paths) {
    if (file.exists(gp)) {
      pg <- gmtPathways(gp)
      pathways <- c(pathways, pg)
    }
  }
  pathways
}

build_ranks <- function(res_df, rank_metric) {
  if (rank_metric == "stat" &&
      "stat" %in% colnames(res_df) &&
      !all(is.na(res_df$stat))) {
    col <- "stat"
  } else {
    col <- "log2FoldChange"
    message("NOTE: falling back to log2FoldChange for ranking.")
  }
  df <- res_df[, c("gene_id", col), drop = FALSE]
  colnames(df)[2] <- "score"
  df <- df[!is.na(df$score), ]
  df <- df[!duplicated(df$gene_id), ]
  setNames(df$score, df$gene_id)
}

# ── Parse arguments ────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--results_a",    type = "character"),
  make_option("--results_b",    type = "character"),
  make_option("--hox_gmt",      type = "character"),
  make_option("--gsea_dir_a",   type = "character"),
  make_option("--gsea_dir_b",   type = "character"),
  make_option("--outdir",       type = "character"),
  make_option("--contrast_a",   type = "character", default = "contrast_A"),
  make_option("--contrast_b",   type = "character", default = "contrast_B"),
  make_option("--collections",  type = "character", default = "H"),
  make_option("--rank_metric",  type = "character", default = "stat"),
  make_option("--min_size",     type = "integer",   default = 15L),
  make_option("--max_size",     type = "integer",   default = 500L),
  make_option("--nperm",        type = "integer",   default = 1000L),
  make_option("--delta_nes",    type = "character", default = "TRUE"),
  make_option("--residual_rank",type = "character", default = "TRUE"),
  make_option("--custom_gmts",  type = "character", default = "")
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

do_delta_nes    <- toupper(args$delta_nes)    == "TRUE"
do_residual     <- toupper(args$residual_rank) == "TRUE"

collections <- str_split(trimws(args$collections), ",")[[1]]
collections <- collections[nchar(trimws(collections)) > 0]

custom_gmt_paths <- character(0)
if (nchar(trimws(args$custom_gmts)) > 0) {
  custom_gmt_paths <- str_split(trimws(args$custom_gmts), ",")[[1]]
  custom_gmt_paths <- custom_gmt_paths[nchar(trimws(custom_gmt_paths)) > 0]
}

# ── Load DESeq2 results ────────────────────────────────────────────────────────
res_a <- read.csv(args$results_a, stringsAsFactors = FALSE)
res_b <- read.csv(args$results_b, stringsAsFactors = FALSE)

# ── Pre-load gene-set pathways for all collections ────────────────────────────
# Load pathways once per collection and cache them; used by both ΔNES and
# residual-rank sections to avoid redundant msigdbr downloads.
message("Loading gene-set pathways for collections: ", paste(collections, collapse = ", "))
pathways_cache <- setNames(
  lapply(collections, load_pathways_for_collection,
         hox_gmt = args$hox_gmt, custom_gmt_paths = custom_gmt_paths),
  collections
)

# ── ΔNES analysis ────────────────────────────────────────────────────────────
if (do_delta_nes) {
  message("Computing ΔNES across collections...")
  delta_rows <- list()

  for (coll_slug in collections) {
    csv_a <- file.path(args$gsea_dir_a, paste0(coll_slug, "_results.csv"))
    csv_b <- file.path(args$gsea_dir_b, paste0(coll_slug, "_results.csv"))

    if (!file.exists(csv_a)) {
      warning(sprintf("GSEA CSV not found for contrast A, collection %s: %s", coll_slug, csv_a))
      next
    }
    if (!file.exists(csv_b)) {
      warning(sprintf("GSEA CSV not found for contrast B, collection %s: %s", coll_slug, csv_b))
      next
    }

    dt_a <- fread(csv_a)[, .(pathway, NES_A = NES, padj_A = padj, size_A = size)]
    dt_b <- fread(csv_b)[, .(pathway, NES_B = NES, padj_B = padj, size_B = size)]

    merged <- merge(dt_a, dt_b, by = "pathway", all = TRUE)
    merged[, delta_NES := NES_A - NES_B]
    merged[, collection := coll_slug]

    # Per-collection CSV
    out_coll_csv <- file.path(args$outdir, paste0("delta_nes_", coll_slug, ".csv"))
    fwrite(merged[order(-abs(delta_NES))], out_coll_csv)
    message(sprintf("  %s: %d pathways → %s", coll_slug, nrow(merged), out_coll_csv))

    delta_rows <- c(delta_rows, list(merged))
  }

  if (length(delta_rows) == 0) {
    message("WARNING: No collections could be compared (missing GSEA CSVs). Writing empty summary.")
  }

  # Summary across all collections
  if (length(delta_rows) > 0) {
    summary_dt <- rbindlist(delta_rows, fill = TRUE)
    summary_dt <- summary_dt[order(collection, -abs(delta_NES))]
    summary_csv <- file.path(args$outdir, "delta_nes_summary.csv")
    fwrite(summary_dt, summary_csv)
    message("ΔNES summary written to: ", summary_csv)

    # Bar plot: top pathways by |ΔNES| across all collections
    top_plot <- head(summary_dt[!is.na(delta_NES)][order(-abs(delta_NES))], 40)
    if (nrow(top_plot) > 0) {
      top_plot[, label := str_trunc(pathway, 55)]
      p <- ggplot(top_plot,
                  aes(x = delta_NES, y = reorder(label, delta_NES),
                      fill = collection)) +
        geom_col() +
        geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
        labs(
          title = sprintf("ΔNES: %s − %s", args$contrast_a, args$contrast_b),
          x = "ΔNES (NES_A − NES_B)",
          y = NULL,
          fill = "Collection"
        ) +
        theme_bw(base_size = 11)
      out_plot <- file.path(args$outdir, "delta_nes_barplot.pdf")
      ggsave(out_plot, plot = p,
             width  = 10,
             height = max(5, nrow(top_plot) * 0.3 + 2))
      message("ΔNES barplot written to: ", out_plot)
    }
  } else {
    # Write empty summary so Snakemake output target is satisfied
    fwrite(data.table(), file.path(args$outdir, "delta_nes_summary.csv"))
  }
} else {
  # Write empty summary so Snakemake output target is satisfied
  fwrite(data.table(), file.path(args$outdir, "delta_nes_summary.csv"))
}

# ── Residual-rank GSEA ────────────────────────────────────────────────────────
if (do_residual) {
  message("Running residual-rank fgsea...")
  resid_dir <- file.path(args$outdir, "residual_rank")
  dir.create(resid_dir, recursive = TRUE, showWarnings = FALSE)

  ranks_a <- build_ranks(res_a, args$rank_metric)
  ranks_b <- build_ranks(res_b, args$rank_metric)

  # Align on common genes
  common_genes <- intersect(names(ranks_a), names(ranks_b))
  if (length(common_genes) < 100) {
    warning(sprintf("Only %d common genes between contrasts; residual GSEA may be unreliable.",
                    length(common_genes)))
  }
  ranks_resid <- sort(ranks_a[common_genes] - ranks_b[common_genes], decreasing = TRUE)

  for (coll_slug in collections) {
    pathways <- pathways_cache[[coll_slug]]

    set.seed(42)
    res_fgsea <- fgseaMultilevel(
      pathways    = pathways,
      stats       = ranks_resid,
      minSize     = args$min_size,
      maxSize     = args$max_size,
      nPermSimple = args$nperm
    )

    dt <- as.data.table(res_fgsea)
    dt[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
    dt <- dt[order(padj, -abs(NES))]

    out_csv <- file.path(resid_dir, paste0("residual_rank_", coll_slug, ".csv"))
    fwrite(dt, out_csv)
    message(sprintf("  Residual rank %s: %d pathways → %s", coll_slug, nrow(dt), out_csv))
  }
} else {
  # Create directory so Snakemake directory() output is satisfied
  dir.create(file.path(args$outdir, "residual_rank"),
             recursive = TRUE, showWarnings = FALSE)
}

message("Contrast comparison analysis complete. Outputs in: ", args$outdir)
