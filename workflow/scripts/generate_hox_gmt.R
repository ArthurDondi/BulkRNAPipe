#!/usr/bin/env Rscript
# generate_hox_gmt.R — Build custom HOX gene-set GMT for BulkRNAPipe GSEA
#
# Reads the gene universe from a DESeq2 results.csv and writes a GMT file
# containing two gene sets:
#   HOX_ALL   — all genes matching ^HOX[ABCD][0-9]+$
#   HOXB_ONLY — all genes matching ^HOXB[0-9]+$
#
# Inputs:
#   --results : path to a DESeq2 results.csv (gene_id column used as universe)
#   --outdir  : directory to write hox.gmt into

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--results", type = "character",
              help = "DESeq2 results CSV (gene_id column provides universe)"),
  make_option("--outdir",  type = "character",
              help = "Output directory for hox.gmt")
)

opt  <- parse_args(OptionParser(option_list = option_list))
args <- opt

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load gene universe ────────────────────────────────────────────────────────
res        <- read.csv(args$results, stringsAsFactors = FALSE)
all_genes  <- res$gene_id

# ── Build HOX gene sets ───────────────────────────────────────────────────────
hox_all   <- grep("^HOX[ABCD][0-9]+$", all_genes, value = TRUE, ignore.case = FALSE)
hoxb_only <- grep("^HOXB[0-9]+$",      all_genes, value = TRUE, ignore.case = FALSE)

# ── Write GMT ─────────────────────────────────────────────────────────────────
# GMT format: set_name\tdescription\tgene1\tgene2\t...
gmt_path <- file.path(args$outdir, "hox.gmt")

write_gmt_line <- function(name, description, genes, con) {
  line <- paste(c(name, description, genes), collapse = "\t")
  writeLines(line, con)
}

con <- file(gmt_path, open = "wt")
write_gmt_line("HOX_ALL",   "All HOX cluster genes (HOXA/B/C/D) detected in gene universe",
               hox_all,   con)
write_gmt_line("HOXB_ONLY", "HOXB cluster genes only detected in gene universe",
               hoxb_only, con)
close(con)

message(sprintf(
  "HOX GMT written to %s\n  HOX_ALL   : %d genes\n  HOXB_ONLY : %d genes",
  gmt_path, length(hox_all), length(hoxb_only)
))
