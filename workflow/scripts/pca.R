#!/usr/bin/env Rscript
# PCA plots for BulkRNAPipe – all samples, PC1 vs PC2 with % variance explained
#
# Two PCAs are produced:
#   pca.pdf                  – every gene in the count matrix (including any
#                              transgene / vector features such as EGFP or
#                              mCherry that were added to the reference)
#   pca_no_vector_genes.pdf  – the same PCA after removing the genes matched by
#                              --exclude_genes / --exclude_gene_patterns /
#                              --exclude_contigs.  Genes are dropped *before*
#                              the DESeq2 dataset is built, so size factors and
#                              the VST are re-estimated on the human-only
#                              matrix (vector transcripts can take up a sizeable
#                              share of the library in transduced samples).
#
# Inputs:
#   --counts  : featureCounts output table (counts.txt)
#   --outdir  : directory for the output PDFs
#   --samples : "sample1:condition1,sample2:condition2,..." mapping

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",  type = "character", help = "featureCounts output file"),
  make_option("--outdir",  type = "character", help = "Output directory"),
  make_option("--samples", type = "character",
              help = "Comma-separated sample:condition pairs"),
  make_option("--exclude_genes", type = "character", default = "",
              help = "Comma-separated gene IDs to drop from the second PCA (case-insensitive exact match)"),
  make_option("--exclude_gene_patterns", type = "character", default = "",
              help = "Comma-separated regular expressions; matching gene IDs are dropped (case-insensitive)"),
  make_option("--exclude_contigs", type = "character", default = "",
              help = "Comma-separated reference contig names; genes located only on these contigs are dropped")
)

args <- parse_args(OptionParser(option_list = option_list))

outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

split_csv <- function(x) {
  if (is.null(x)) return(character(0))
  parts <- trimws(strsplit(x, ",")[[1]])
  parts[nzchar(parts)]
}

exclude_genes    <- split_csv(args$exclude_genes)
exclude_patterns <- split_csv(args$exclude_gene_patterns)
exclude_contigs  <- split_csv(args$exclude_contigs)

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
count_cols  <- colnames(raw)[7:ncol(raw)]
# Strip path prefix and .Aligned.sortedByCoord.out.bam suffix so column names
# match the sample names in the config.
clean_names <- sub(".*/", "", count_cols)
clean_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
counts      <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_names

# featureCounts joins the contigs of all exons of a gene with ';'.
gene_contigs <- setNames(as.character(raw$Chr), raw$Geneid)

# Keep only samples present in sample_df
shared    <- intersect(colnames(counts), rownames(sample_df))
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]

# ─── Identify vector / transgene features to exclude ─────────────────────────
gene_ids <- rownames(counts)
drop_mask <- rep(FALSE, length(gene_ids))
matched_by <- rep(NA_character_, length(gene_ids))

record_hits <- function(hit, label) {
  new <- hit & !drop_mask
  drop_mask  <<- drop_mask | hit
  matched_by[new] <<- label
}

for (g in exclude_genes) {
  hit <- tolower(gene_ids) == tolower(g)
  if (!any(hit)) {
    message("WARNING: exclude_genes entry '", g, "' matched no gene - check the ",
            "gene_id used in your GTF.")
  }
  record_hits(hit, paste0("gene:", g))
}

for (p in exclude_patterns) {
  hit <- grepl(p, gene_ids, ignore.case = TRUE, perl = TRUE)
  if (!any(hit)) {
    message("WARNING: exclude_gene_patterns entry '", p, "' matched no gene.")
  }
  record_hits(hit, paste0("pattern:", p))
}

if (length(exclude_contigs) > 0) {
  contigs_per_gene <- strsplit(gene_contigs[gene_ids], ";", fixed = TRUE)
  hit <- vapply(contigs_per_gene, function(cs) {
    cs <- unique(cs)
    length(cs) > 0 && all(cs %in% exclude_contigs)
  }, logical(1))
  if (!any(hit)) {
    message("WARNING: exclude_contigs matched no gene.")
  }
  record_hits(hit, "contig")
}

excluded_ids <- gene_ids[drop_mask]

# Report how much of each library the excluded features take up - a large
# fraction means they also distorted the size factors of the "all genes" PCA.
lib_size <- colSums(counts)
if (length(excluded_ids) > 0) {
  excl_counts <- counts[excluded_ids, , drop = FALSE]
  report <- data.frame(
    gene_id    = excluded_ids,
    contig     = vapply(strsplit(gene_contigs[excluded_ids], ";", fixed = TRUE),
                        function(cs) paste(unique(cs), collapse = ";"), character(1)),
    matched_by = matched_by[drop_mask],
    total_count = rowSums(excl_counts),
    stringsAsFactors = FALSE
  )
  pct <- sweep(excl_counts, 2, lib_size, "/") * 100
  colnames(pct) <- paste0("pct_", colnames(pct))
  report <- cbind(report, as.data.frame(pct))
  message("Excluded features and their % of each library:")
  print(report, row.names = FALSE)
} else {
  report <- data.frame(
    gene_id = character(0), contig = character(0), matched_by = character(0),
    total_count = numeric(0), stringsAsFactors = FALSE
  )
  message("No genes matched the exclusion lists - both PCAs are identical.")
}
write.csv(report, file.path(outdir, "pca_excluded_genes.csv"), row.names = FALSE)

# ─── PCA helper ──────────────────────────────────────────────────────────────
make_pca <- function(count_mat, coldata, title, outfile) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = coldata,
    design    = ~ condition
  )

  # Remove genes with very low counts (< 10 reads across all samples)
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Blind VST for sample-level QC (design information not used)
  vst_data <- vst(dds, blind = TRUE)

  # PCA on the top 500 most-variable genes (plotPCA default)
  pca_df  <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
  pct_var <- round(100 * attr(pca_df, "percentVar"), 1)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    labs(
      title  = title,
      x      = paste0("PC1: ", pct_var[1], "% variance"),
      y      = paste0("PC2: ", pct_var[2], "% variance"),
      colour = "Condition"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(outfile, plot = p, width = 7, height = 6)
  message("PCA plot written to: ", outfile)
}

# ─── PCA 1: every gene, vector features included ─────────────────────────────
make_pca(counts, sample_df,
         "PCA - all samples (all genes)",
         file.path(outdir, "pca.pdf"))

# ─── PCA 2: vector / transgene features removed ──────────────────────────────
# Dropping them here (rather than after normalisation) means size factors and
# the VST are estimated on the human-only matrix.
counts_no_vector <- counts[!drop_mask, , drop = FALSE]
subtitle <- if (length(excluded_ids) > 0) {
  paste0("PCA - all samples (", length(excluded_ids),
         " vector/transgene feature(s) removed)")
} else {
  "PCA - all samples (no vector/transgene features to remove)"
}
make_pca(counts_no_vector, sample_df,
         subtitle,
         file.path(outdir, "pca_no_vector_genes.pdf"))
