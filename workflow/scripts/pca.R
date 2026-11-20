#!/usr/bin/env Rscript
# pca.R — sample-level PCA for BulkRNAPipe
#
# Runs one (sample set x gene view) combination and writes PC1/PC2, PC3/PC4, a
# scree plot, the PC coordinates and the list of excluded features.
#
# Gene views:
#   all_genes               every gene in the count matrix
#   no_vector_genes         minus the vector/reporter features (EGFP, mCherry, ...)
#   no_vector_no_transgene  minus those, and minus the transgene and its variants
#
# The last view matters more than it looks: the transgene is by construction the
# largest single difference between the groups, so a PCA that includes it partly
# re-plots the experimental design rather than its consequences.
#
# Features are dropped from the raw count matrix *before* the DESeqDataSet is
# built, so size factors and the VST are re-estimated on the retained genes.
#
# A supplementary panel applies limma::removeBatchEffect. That panel is for
# LOOKING ONLY: when the batch label is a deterministic function of the
# condition, the effect it removes is not identifiable from the data, and the
# amount removed is set by the assumption, not measured. It must never be used
# to justify a differential-expression result.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(limma)
})

# ─── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",  type = "character", help = "featureCounts output file"),
  make_option("--outdir",  type = "character", help = "Output directory"),
  make_option("--samples", type = "character",
              help = "Comma-separated sample:condition pairs"),
  make_option("--view", type = "character", default = "all_genes",
              help = "all_genes | no_vector_genes | no_vector_no_transgene"),
  make_option("--sampleset", type = "character", default = "all",
              help = "Label for the sample set, used in plot titles"),
  make_option("--exclude_genes", type = "character", default = "",
              help = "Comma-separated vector/reporter gene IDs (case-insensitive exact match)"),
  make_option("--exclude_gene_patterns", type = "character", default = "",
              help = "Comma-separated regexes; matching gene IDs are treated as vector features"),
  make_option("--exclude_contigs", type = "character", default = "",
              help = "Comma-separated contig names; genes located only there are vector features"),
  make_option("--transgene_genes", type = "character", default = "",
              help = "Comma-separated transgene gene IDs, dropped by the no_vector_no_transgene view"),
  make_option("--batch", type = "character", default = "",
              help = "Comma-separated sample:batch pairs for the removeBatchEffect panel"),
  make_option("--ntop", type = "integer", default = 500,
              help = "Number of most-variable genes used for the PCA [default %default]")
)

args   <- parse_args(OptionParser(option_list = option_list))
outdir <- args$outdir
view   <- args$view
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out <- function(suffix) file.path(outdir, paste0(view, "_", suffix))

`%||%` <- function(a, b) if (!is.null(a)) a else b

split_csv <- function(x) {
  if (is.null(x)) return(character(0))
  parts <- trimws(strsplit(x, ",")[[1]])
  parts[nzchar(parts)]
}

parse_pairs <- function(x, value_name) {
  pairs <- split_csv(x)
  if (length(pairs) == 0) return(NULL)
  df <- do.call(rbind, lapply(pairs, function(p) {
    parts <- strsplit(p, ":")[[1]]
    data.frame(sample = parts[1], value = parts[2], stringsAsFactors = FALSE)
  }))
  names(df)[2] <- value_name
  rownames(df) <- df$sample
  df
}

exclude_genes    <- split_csv(args$exclude_genes)
exclude_patterns <- split_csv(args$exclude_gene_patterns)
exclude_contigs  <- split_csv(args$exclude_contigs)
transgene_genes  <- split_csv(args$transgene_genes)

# ─── Sample table ────────────────────────────────────────────────────────────
sample_df <- parse_pairs(args$samples, "condition")
if (is.null(sample_df)) stop("--samples is empty")
sample_df$condition <- factor(sample_df$condition)

batch_df <- parse_pairs(args$batch, "batch")

# ─── Load featureCounts output ───────────────────────────────────────────────
raw <- read.table(args$counts, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE)

count_cols  <- colnames(raw)[7:ncol(raw)]
clean_names <- sub(".*/", "", count_cols)
clean_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
counts      <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_names

# featureCounts joins the contigs of all exons of a gene with ';'.
gene_contigs <- setNames(as.character(raw$Chr), raw$Geneid)

shared <- intersect(colnames(counts), rownames(sample_df))
if (length(shared) < 3) {
  stop("Only ", length(shared), " sample(s) shared between the count matrix and ",
       "--samples; a PCA needs at least 3.")
}
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]
sample_df$condition <- droplevels(sample_df$condition)
if (!is.null(batch_df)) batch_df <- batch_df[shared, , drop = FALSE]

# ─── Which features does this view drop? ─────────────────────────────────────
gene_ids   <- rownames(counts)
drop_mask  <- rep(FALSE, length(gene_ids))
matched_by <- rep(NA_character_, length(gene_ids))

record_hits <- function(hit, label) {
  new <- hit & !drop_mask
  drop_mask  <<- drop_mask | hit
  matched_by[new] <<- label
}

match_exact <- function(ids, wanted, label) {
  for (g in wanted) {
    hit <- tolower(ids) == tolower(g)
    if (!any(hit)) {
      message("WARNING: '", g, "' matched no gene in the count matrix - check ",
              "the gene_id used in your GTF.")
    }
    record_hits(hit, paste0(label, ":", g))
  }
}

if (view %in% c("no_vector_genes", "no_vector_no_transgene")) {
  match_exact(gene_ids, exclude_genes, "vector")
  for (p in exclude_patterns) {
    hit <- grepl(p, gene_ids, ignore.case = TRUE, perl = TRUE)
    if (!any(hit)) message("WARNING: pattern '", p, "' matched no gene.")
    record_hits(hit, paste0("vector_pattern:", p))
  }
  if (length(exclude_contigs) > 0) {
    hit <- vapply(strsplit(gene_contigs[gene_ids], ";", fixed = TRUE), function(cs) {
      cs <- unique(cs)
      length(cs) > 0 && all(cs %in% exclude_contigs)
    }, logical(1))
    if (!any(hit)) message("WARNING: exclude_contigs matched no gene.")
    record_hits(hit, "vector_contig")
  }
}
if (view == "no_vector_no_transgene") {
  match_exact(gene_ids, transgene_genes, "transgene")
}

excluded_ids <- gene_ids[drop_mask]
lib_size     <- colSums(counts)

if (length(excluded_ids) > 0) {
  excl <- counts[excluded_ids, , drop = FALSE]
  report <- data.frame(
    gene_id     = excluded_ids,
    contig      = vapply(strsplit(gene_contigs[excluded_ids], ";", fixed = TRUE),
                         function(cs) paste(unique(cs), collapse = ";"), character(1)),
    matched_by  = matched_by[drop_mask],
    total_count = rowSums(excl),
    stringsAsFactors = FALSE
  )
  pct <- sweep(excl, 2, lib_size, "/") * 100
  colnames(pct) <- paste0("pct_", colnames(pct))
  report <- cbind(report, as.data.frame(pct))
  message("View '", view, "' drops ", length(excluded_ids), " feature(s):")
  print(report, row.names = FALSE)
} else {
  report <- data.frame(gene_id = character(0), contig = character(0),
                       matched_by = character(0), total_count = numeric(0),
                       stringsAsFactors = FALSE)
  message("View '", view, "' drops no features.")
}
write.csv(report, out("excluded_genes.csv"), row.names = FALSE)

counts <- counts[!drop_mask, , drop = FALSE]

# ─── VST ─────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_df,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vsd <- vst(dds, blind = TRUE)          # blind: the design is not used
mat <- assay(vsd)

# ─── PCA helpers ─────────────────────────────────────────────────────────────
run_pca <- function(m, ntop) {
  rv     <- apply(m, 1, var)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, sum(rv > 0)))]
  p      <- prcomp(t(m[select, , drop = FALSE]))
  list(x = p$x, pct = 100 * p$sdev^2 / sum(p$sdev^2), n_genes = length(select))
}

pc_label <- function(pct, i) paste0("PC", i, ": ", round(pct[i], 1), "% variance")

plot_pc <- function(pca, i, j, title, subtitle, outfile, caption = NULL,
                    caption_colour = "grey30") {
  if (ncol(pca$x) < j) {
    message("Skipping ", outfile, ": only ", ncol(pca$x), " components available.")
    pdf(outfile, width = 7, height = 6)
    plot.new()
    text(0.5, 0.5, paste0("Not enough samples for PC", j), cex = 1.2)
    dev.off()
    return(invisible(NULL))
  }
  df <- data.frame(
    x = pca$x[, i], y = pca$x[, j],
    condition = sample_df$condition,
    name = rownames(sample_df)
  )
  p <- ggplot(df, aes(x = x, y = y, colour = condition, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = pc_label(pca$pct, i), y = pc_label(pca$pct, j), colour = "Condition") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.caption = element_text(hjust = 0, colour = caption_colour, size = 8))
  ggsave(outfile, plot = p, width = 7, height = 6.4)
  message("Written: ", outfile)
}

plot_scree <- function(pca, title, outfile) {
  n  <- min(10, length(pca$pct))
  df <- data.frame(PC = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
                   pct = pca$pct[seq_len(n)])
  df$cumulative <- cumsum(df$pct)
  p <- ggplot(df, aes(x = PC, y = pct)) +
    geom_col(fill = "grey40") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3) +
    geom_line(aes(y = cumulative, group = 1), colour = "steelblue") +
    geom_point(aes(y = cumulative), colour = "steelblue", size = 1.5) +
    labs(title = title, subtitle = "bars: variance per PC | line: cumulative",
         x = NULL, y = "% variance explained") +
    expand_limits(y = 105) +
    theme_bw(base_size = 12)
  ggsave(outfile, plot = p, width = 7, height = 5)
  message("Written: ", outfile)
}

view_label <- c(
  all_genes              = "all genes",
  no_vector_genes        = "vector/reporter features removed",
  no_vector_no_transgene = "vector/reporter + transgene features removed"
)[view]
if (is.na(view_label)) view_label <- view

title    <- paste0("PCA - ", args$sampleset, " samples")
subtitle <- paste0(view_label, " | ", ncol(mat), " samples")

# ─── Main PCA ────────────────────────────────────────────────────────────────
pca <- run_pca(mat, args$ntop)
subtitle_n <- paste0(subtitle, " | top ", pca$n_genes, " variable genes")

plot_pc(pca, 1, 2, title, subtitle_n, out("pc12.pdf"))
plot_pc(pca, 3, 4, title, subtitle_n, out("pc34.pdf"))
plot_scree(pca, paste0(title, " - variance explained"), out("scree.pdf"))

coords <- data.frame(sample = rownames(sample_df), condition = sample_df$condition,
                     pca$x[, seq_len(min(6, ncol(pca$x))), drop = FALSE],
                     check.names = FALSE)
attr(coords, "percentVar") <- NULL
write.csv(coords, out("coords.csv"), row.names = FALSE)

pct_df <- data.frame(PC = paste0("PC", seq_along(pca$pct)),
                     percent_variance = round(pca$pct, 3))
write.csv(head(pct_df, 10), out("variance_explained.csv"), row.names = FALSE)

# ─── Supplementary: removeBatchEffect panel (LOOKING ONLY) ───────────────────
placeholder <- function(msg, outfile) {
  pdf(outfile, width = 7, height = 6)
  plot.new()
  text(0.5, 0.55, "No batch-removal panel", cex = 1.3)
  text(0.5, 0.42, msg, cex = 0.9)
  dev.off()
  message("Written (placeholder): ", outfile)
}

if (is.null(batch_df)) {
  placeholder("PCA.batch is not configured in the config file.",
              out("pc12_batch_removed.pdf"))
} else {
  batch <- factor(batch_df$batch)
  if (nlevels(batch) < 2) {
    placeholder(paste0("Only one batch level ('", levels(batch)[1],
                       "') among these samples."),
                out("pc12_batch_removed.pdf"))
  } else {
    # Is the batch effect separable from the condition effect?  It is only if
    # some condition spans more than one batch.  When every condition sits in
    # exactly one batch, batch is a deterministic function of condition: the two
    # are the same variable wearing different names, and no amount of data can
    # say which of them moved a gene.
    tab              <- table(sample_df$condition, batch)
    batches_per_cond <- apply(tab, 1, function(x) sum(x > 0))
    estimable        <- any(batches_per_cond > 1)

    if (estimable) {
      # Protect the condition means: removeBatchEffect only takes out the
      # component of the batch that is orthogonal to the design.
      corrected <- removeBatchEffect(mat, batch = batch,
                                     design = model.matrix(~ sample_df$condition))
      note <- "condition means preserved"
      warn <- NULL
    } else {
      # Aliased with condition.  Passing the design here would leave the batch
      # coefficient inestimable (limma returns NA and the result is unusable),
      # so the batch means are subtracted outright.  That does produce the
      # tidy-looking figure people expect - by deleting the real difference
      # between the groups along with any technical one.  There is no way to
      # separate the two, which is the whole point.
      corrected <- removeBatchEffect(mat, batch = batch)
      note <- "batch is aliased with condition"
      warn <- paste0(
        "DO NOT USE AS EVIDENCE: '", paste(levels(batch), collapse = " vs "),
        "' is a relabelling of the condition, so this panel removes\n",
        "the biological difference between those groups along with any technical one. ",
        "The two are not separable in this design.")
      message("WARNING: batch is a deterministic function of condition. The ",
              "batch-removed panel subtracts the group difference itself and is ",
              "an illustration only - it is not evidence of anything.")
    }

    pca_c <- run_pca(corrected, args$ntop)
    # The warning is stamped onto the figure itself: a PDF outlives its log file
    # and this one is the kind that travels into slide decks on its own.
    plot_pc(pca_c, 1, 2,
            paste0(title, " - batch removed (VISUALISATION ONLY)"),
            paste0(subtitle_n, " | removeBatchEffect on '",
                   paste(levels(batch), collapse = " vs "), "' | ", note),
            out("pc12_batch_removed.pdf"),
            caption = warn %||% paste0(
              "Visualisation only: the removed effect is set by the assumed ",
              "batch labels, not estimated from independent data."),
            caption_colour = if (is.null(warn)) "grey30" else "#B2182B")
    message("NOTE: the batch-removed panel is for visualisation only. The size ",
            "of the removed effect is set by the assumed batch labels, not ",
            "estimated from independent data.")
  }
}

message("PCA complete for view '", view, "', sample set '", args$sampleset, "'.")
