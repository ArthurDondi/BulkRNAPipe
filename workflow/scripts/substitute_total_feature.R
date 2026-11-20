#!/usr/bin/env Rscript
# substitute_total_feature.R
#
# Produces a derived featureCounts-format count matrix, downstream of
# design_qc. quantify/counts.txt itself is only read here, never modified.
# Two independent adjustments, either of which can be skipped:
#
#   substitution   design_qc's pooled "total feature" (DesignQCTotalCounts,
#                  counted with -Q 0 -M --fraction over several gene_ids)
#                  stands in for one gene id, and the gene_ids it pools are
#                  dropped. Skipped when --target_gene is empty.
#   extra drops    Genes removed outright, no substitute (e.g. reporter genes
#                  like EGFP/mCherry that should not enter DESeq2 at all).
#                  Skipped when --extra_drop_genes is empty.
#
# --counts and --total are both standard featureCounts output (a "#" comment
# line, then Geneid/Chr/Start/End/Strand/Length + one column per BAM), so
# lining up the substitution only has to match sample columns and swap a row.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Original featureCounts counts.txt"),
  make_option("--total",  type = "character", default = "",
              help = "design_qc total_feature_counts.txt (fractional, -M --fraction); required with --target_gene"),
  make_option("--total_name",  type = "character", default = "",
              help = "Feature name (gene_id) inside --total"),
  make_option("--target_gene", type = "character", default = "",
              help = "Gene id the pooled total should stand in for downstream; empty to skip substitution"),
  make_option("--drop_genes",  type = "character", default = "",
              help = "Comma-separated gene ids to drop from --counts (the ones pooled into --total)"),
  make_option("--extra_drop_genes", type = "character", default = "",
              help = "Comma-separated gene ids to remove outright, no substitute"),
  make_option("--output", type = "character", help = "Output path for the derived counts table")
)
args <- parse_args(OptionParser(option_list = option_list))

read_fc <- function(path) {
  read.table(path, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
}

split_genes <- function(x) {
  g <- trimws(strsplit(x, ",")[[1]])
  g[nzchar(g)]
}

counts <- read_fc(args$counts)
out    <- counts

if (nzchar(args$target_gene)) {
  drop_genes <- split_genes(args$drop_genes)
  if (length(drop_genes) == 0) {
    stop("--drop_genes is empty - nothing to pool into '", args$target_gene, "'.")
  }
  if (!nzchar(args$total) || !nzchar(args$total_name)) {
    stop("--target_gene is set but --total / --total_name is missing.")
  }

  total <- read_fc(args$total)
  total_row <- total[total$Geneid == args$total_name, , drop = FALSE]
  if (nrow(total_row) != 1) {
    stop("Expected exactly one row for '", args$total_name, "' in ", args$total,
         ", found ", nrow(total_row), ".")
  }

  # Both files come from featureCounts on the same BAMs, but column ORDER
  # between two separate runs is not guaranteed - match by name, not position.
  sample_cols <- colnames(counts)[7:ncol(counts)]
  missing <- setdiff(sample_cols, colnames(total_row))
  if (length(missing) > 0) {
    stop("--total is missing sample column(s) present in --counts: ",
         paste(missing, collapse = ", "))
  }
  total_row <- total_row[, c(colnames(counts)[1:6], sample_cols), drop = FALSE]

  # Fractional multi-mapped counts are fine for the design_qc panel but DESeq2
  # expects integers - round once here for every downstream consumer.
  total_row[sample_cols] <- lapply(total_row[sample_cols], function(x) round(as.numeric(x)))
  total_row$Geneid <- args$target_gene

  out <- out[!out$Geneid %in% drop_genes, , drop = FALSE]
  out <- rbind(out, total_row)
  message("Substituted '", args$target_gene, "' with the pooled '", args$total_name,
          "' (dropped ", paste(drop_genes, collapse = ", "), ")")
}

extra_drop <- split_genes(args$extra_drop_genes)
if (length(extra_drop) > 0) {
  out <- out[!out$Geneid %in% extra_drop, , drop = FALSE]
  message("Dropped ", paste(extra_drop, collapse = ", "), " outright (no substitute)")
}

if (!nzchar(args$target_gene) && length(extra_drop) == 0) {
  stop("Neither --target_gene nor --extra_drop_genes is set - nothing to do.")
}

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write.table(out, args$output, sep = "\t", row.names = FALSE, quote = FALSE)

message("Wrote ", nrow(out), " genes to ", args$output)
