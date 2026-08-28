#!/usr/bin/env Rscript
# substitute_total_feature.R
#
# Produces a derived featureCounts-format count matrix, downstream of
# design_qc, where design_qc's pooled "total feature" (DesignQCTotalCounts,
# counted with -Q 0 -M --fraction over several gene_ids) stands in for one
# gene id, and the gene_ids it pools are dropped. quantify/counts.txt itself
# is never read except as input here - it is not modified.
#
# --counts and --total are both standard featureCounts output (a "#" comment
# line, then Geneid/Chr/Start/End/Strand/Length + one column per BAM), so this
# only has to line up the sample columns and swap one row.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Original featureCounts counts.txt"),
  make_option("--total",  type = "character", help = "design_qc total_feature_counts.txt (fractional, -M --fraction)"),
  make_option("--total_name",  type = "character", help = "Feature name (gene_id) inside --total"),
  make_option("--target_gene", type = "character", help = "Gene id the pooled total should stand in for downstream"),
  make_option("--drop_genes",  type = "character", help = "Comma-separated gene ids to drop from --counts (the ones pooled into --total)"),
  make_option("--output", type = "character", help = "Output path for the derived counts table")
)
args <- parse_args(OptionParser(option_list = option_list))

read_fc <- function(path) {
  read.table(path, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
}

counts <- read_fc(args$counts)
total  <- read_fc(args$total)

drop_genes <- trimws(strsplit(args$drop_genes, ",")[[1]])
drop_genes <- drop_genes[nzchar(drop_genes)]
if (length(drop_genes) == 0) {
  stop("--drop_genes is empty - nothing to pool into '", args$target_gene, "'.")
}

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

out <- counts[!counts$Geneid %in% drop_genes, , drop = FALSE]
out <- rbind(out, total_row)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write.table(out, args$output, sep = "\t", row.names = FALSE, quote = FALSE)

message("Wrote ", nrow(out), " genes to ", args$output, " (dropped ",
        paste(drop_genes, collapse = ", "), "; '", args$target_gene,
        "' is now the pooled '", args$total_name, "')")
