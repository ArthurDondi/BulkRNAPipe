#!/usr/bin/env Rscript
# shared_genes.R — genes significant in the same direction in two DESeq2
# contrasts that share a denominator (e.g. ATRX_FL_vs_EmptyVector and
# ATRX_IFF_vs_EmptyVector), and an over-representation test on them.
#
# "up"   = significantly higher in each numerator  (padj < padj_thr, log2FC >=  lfc_thr in both)
# "down" = significantly higher in the denominator (padj < padj_thr, log2FC <= -lfc_thr in both)
#
# The shared genes are a set, not a ranking, so the enrichment is an
# over-representation test (one-sided hypergeometric, as in clusterProfiler's
# enricher()), not GSEA. Universe: genes with a non-NA padj in both contrasts.
# BH adjustment within each direction x collection, over the sets holding at
# least one shared gene (the sets clusterProfiler's enricher() tests).
#
# Outputs (in --outdir):
#   genes.csv                         every universe gene: log2FC / padj in both
#                                     contrasts, and its up / down group
#                                     (both, A_only, B_only, or empty)
#   ora_{up,down}_{collection}.csv    one row per gene set tested
#
# Usage:
#   Rscript shared_genes.R --results_a A/results.csv --results_b B/results.csv \
#     --collections H,KEGG_CURRENT,C5:GO:BP --kegg_gmt kegg_current.gmt --outdir out

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--results_a",   type = "character", help = "DESeq2 results.csv of contrast A"),
  make_option("--results_b",   type = "character", help = "DESeq2 results.csv of contrast B"),
  make_option("--padj_thr",    type = "double", default = 0.05, help = "DESeq2 padj threshold [default %default]"),
  make_option("--lfc_thr",     type = "double", default = 1, help = "absolute log2FC threshold [default %default]"),
  make_option("--exclude_genes", type = "character", default = "",
              help = "comma-separated gene ids left out of every set and the universe"),
  make_option("--collections", type = "character", default = "H,KEGG_CURRENT,C5:GO:BP",
              help = "comma-separated MSigDB collections, or KEGG_CURRENT [default %default]"),
  make_option("--kegg_gmt",    type = "character", default = "",
              help = "GMT of current KEGG pathways (needed for KEGG_CURRENT)"),
  make_option("--min_size",    type = "integer", default = 10, help = "min gene-set size in the universe [default %default]"),
  make_option("--max_size",    type = "integer", default = 500, help = "max gene-set size in the universe [default %default]"),
  make_option("--outdir",      type = "character", help = "output directory")
)
args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Shared genes ─────────────────────────────────────────────────────────────
read_res <- function(f) {
  r <- read.csv(f, stringsAsFactors = FALSE)
  r <- r[!is.na(r$padj) & !is.na(r$log2FoldChange), c("gene_id", "log2FoldChange", "padj")]
  r[!duplicated(r$gene_id), ]
}
a <- read_res(args$results_a)
b <- read_res(args$results_b)
exclude <- trimws(strsplit(args$exclude_genes, ",")[[1]])
genes <- merge(a, b, by = "gene_id", suffixes = c("_A", "_B"))
genes <- genes[!tolower(genes$gene_id) %in% tolower(exclude), ]

call <- function(lfc, padj, sign) padj < args$padj_thr & sign * lfc >= args$lfc_thr
group <- function(in_a, in_b)
  ifelse(in_a & in_b, "both", ifelse(in_a, "A_only", ifelse(in_b, "B_only", "")))
genes$up   <- group(call(genes$log2FoldChange_A, genes$padj_A,  1),
                    call(genes$log2FoldChange_B, genes$padj_B,  1))
genes$down <- group(call(genes$log2FoldChange_A, genes$padj_A, -1),
                    call(genes$log2FoldChange_B, genes$padj_B, -1))
genes <- genes[order(genes$gene_id), ]
write.csv(genes, file.path(args$outdir, "genes.csv"), row.names = FALSE)
universe <- genes$gene_id
message(sprintf("Universe: %d genes | up in both: %d | down in both: %d",
                length(universe), sum(genes$up == "both"), sum(genes$down == "both")))

# ── Gene sets ────────────────────────────────────────────────────────────────
read_gmt <- function(f) {
  l <- strsplit(readLines(f), "\t")
  setNames(lapply(l, function(x) unique(x[-(1:2)])), vapply(l, `[`, "", 1))
}
load_sets <- function(collection) {
  if (identical(collection, "KEGG_CURRENT")) {
    if (!nzchar(args$kegg_gmt) || !file.exists(args$kegg_gmt))
      stop("collection KEGG_CURRENT requires --kegg_gmt pointing at an existing GMT")
    return(read_gmt(args$kegg_gmt))
  }
  parts <- strsplit(collection, ":", fixed = TRUE)[[1]]
  m <- if (length(parts) > 1)
    msigdbr::msigdbr(species = "Homo sapiens", category = parts[1],
                     subcategory = paste(parts[-1], collapse = ":"))
  else msigdbr::msigdbr(species = "Homo sapiens", category = parts[1])
  lapply(split(m$gene_symbol, m$gs_name), unique)
}

# ── Over-representation test ─────────────────────────────────────────────────
ora <- function(hits, sets) {
  empty <- data.frame(pathway = character(), set_size = integer(), overlap = integer(),
                      expected = numeric(), fold_enrichment = numeric(),
                      pvalue = numeric(), padj = numeric(), genes = character())
  sets <- lapply(sets, intersect, universe)
  sets <- sets[lengths(sets) >= args$min_size & lengths(sets) <= args$max_size]
  if (!length(hits) || !length(sets)) return(empty)
  N <- length(universe); n <- length(hits)
  ov <- lapply(sets, intersect, hits)
  k <- lengths(ov); K <- lengths(sets)
  res <- data.frame(pathway = names(sets), set_size = K, overlap = k,
                    expected = K * n / N, fold_enrichment = (k / n) / (K / N),
                    # P(X >= k), X ~ hypergeometric(K successes in N, n drawn)
                    pvalue = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
                    genes = vapply(ov, function(g) paste(sort(g), collapse = "/"), ""),
                    stringsAsFactors = FALSE)
  # As clusterProfiler's enricher(): only sets holding at least one hit are
  # tested, so BH runs over those (a set with no hit cannot be enriched).
  res <- res[res$overlap > 0, , drop = FALSE]
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[, names(empty)]
  res[order(res$pvalue), ]
}

for (col in trimws(strsplit(args$collections, ",")[[1]])) {
  sets <- load_sets(col)
  slug <- gsub(":", "_", col)
  for (dir in c("up", "down")) {
    res <- ora(genes$gene_id[genes[[dir]] == "both"], sets)
    write.csv(res, file.path(args$outdir, sprintf("ora_%s_%s.csv", dir, slug)), row.names = FALSE)
    message(sprintf("%s %s: %d sets with overlap, %d with padj < %s",
                    slug, dir, nrow(res), sum(res$padj < args$padj_thr), args$padj_thr))
  }
}
