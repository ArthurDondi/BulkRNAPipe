#!/usr/bin/env Rscript
# gsea_kegg.R — per-KEGG-category dotplots from an existing GSEA result table.
#
# This does NOT re-run fgsea. It reads the KEGG collection's results.csv, which
# already contains every KEGG pathway that passed the size filters, and slices it
# by KEGG BRITE category. Re-plotting is therefore cheap and cannot disagree with
# the collection-level figure.
#
# Unlike the collection dotplot, every pathway in a category is drawn - not a
# top-N. A category is a fixed, small, pre-specified list, so showing only the
# winners would hide the fact that the rest were tested and came back flat.
#
# Outputs (in --outdir):
#   {slug}_dotplot.pdf   one panel per category
#   {slug}_results.csv   the same rows as a table

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--results", type = "character",
              help = "GSEA results CSV for the KEGG collection"),
  make_option("--categories", type = "character",
              help = "TSV: slug, label, gs_name"),
  make_option("--select", type = "character", default = "",
              help = "Comma-separated slugs to plot; empty = every slug in the TSV"),
  make_option("--audit", type = "character", default = "",
              help = "Pathway audit CSV from gsea.R; lets a missing pathway be reported as filtered-by-size vs absent-from-collection"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--numerator", type = "character", default = ""),
  make_option("--denominator", type = "character", default = ""),
  make_option("--padj_cutoff", type = "double", default = 0.05)
)
args <- parse_args(OptionParser(option_list = option_list))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

cats <- read.delim(args$categories, header = TRUE, sep = "\t",
                   comment.char = "#", stringsAsFactors = FALSE)
need <- c("slug", "label", "gs_name")
if (!all(need %in% colnames(cats))) {
  stop("Category file must have columns: ", paste(need, collapse = ", "))
}
cats <- cats[nzchar(trimws(cats$slug)) & nzchar(trimws(cats$gs_name)), ]

selected <- trimws(strsplit(args$select, ",")[[1]])
selected <- selected[nzchar(selected)]
if (length(selected) == 0) selected <- unique(cats$slug)

for (s in setdiff(selected, unique(cats$slug))) {
  message("WARNING: category '", s, "' is not in ", args$categories,
          " - an empty panel will be written.")
}

res <- read.csv(args$results, stringsAsFactors = FALSE)
message("Loaded ", nrow(res), " pathways from ", args$results)

# The results file only contains gene sets fgsea actually tested, so on its own
# it cannot say WHY a pathway is absent. The audit written by gsea.R lists every
# set in the collection with the number of its genes detected here, which
# separates "the collection does not have it" from "too few of its genes are
# expressed in this dataset" - two problems with different fixes.
audit <- NULL
if (nzchar(args$audit) && file.exists(args$audit)) {
  audit <- read.csv(args$audit, stringsAsFactors = FALSE)
  message("Loaded audit for ", nrow(audit), " gene sets from ", args$audit)
} else if (nzchar(args$audit)) {
  message("WARNING: audit file '", args$audit, "' not found. Missing pathways ",
          "will be reported without a reason.")
}

# Why is `gs` not in the results?
explain_missing <- function(gs) {
  if (is.null(audit)) {
    return("reason unknown (no audit file); it is either absent from the collection or was filtered by size")
  }
  row <- audit[audit$pathway == gs, , drop = FALSE]
  if (nrow(row) == 0) {
    return(paste0("NOT IN THE COLLECTION - no gene set by this name. Check the ",
                  "spelling in the category file, or whether this collection ",
                  "version contains it."))
  }
  if (isTRUE(row$tested[1])) {
    return(sprintf("unexpected: it was testable (%d genes detected) but is absent from the results",
                   row$size_detected[1]))
  }
  if (identical(row$reason[1], "below_min_size")) {
    return(sprintf("FILTERED BY SIZE - only %d of its %d genes are detected here, below min_size=%d",
                   row$size_detected[1], row$size_collection[1], row$min_size[1]))
  }
  sprintf("FILTERED BY SIZE - %d of its %d genes are detected here, above max_size=%d",
          row$size_detected[1], row$size_collection[1], row$max_size[1])
}

# A pathway can be absent for two different reasons and they need different
# fixes, so report them separately rather than as one "missing" count.
empty_panel <- function(pdf_path, msg) {
  pdf(pdf_path, width = 8, height = 2.2)
  plot.new(); text(0.5, 0.5, msg, cex = 1.0)
  dev.off()
  message("Written (empty): ", pdf_path)
}

for (s in selected) {
  rows  <- cats[cats$slug == s, , drop = FALSE]
  label <- if (nrow(rows) > 0) rows$label[1] else s
  wanted <- unique(rows$gs_name)

  hit <- res[res$pathway %in% wanted, , drop = FALSE]
  absent <- setdiff(wanted, res$pathway)
  if (length(absent) > 0) {
    message("WARNING: [", s, "] ", length(absent), " of ", length(wanted),
            " gene sets are not in the results and are not plotted:")
    for (gs in absent) message("           - ", gs, ": ", explain_missing(gs))
  }

  out_pdf <- file.path(args$outdir, paste0(s, "_dotplot.pdf"))
  out_csv <- file.path(args$outdir, paste0(s, "_results.csv"))
  write.csv(hit, out_csv, row.names = FALSE)

  if (nrow(hit) == 0) {
    empty_panel(out_pdf, paste0(label, ": no gene sets found"))
    next
  }

  hit$significant <- !is.na(hit$padj) & hit$padj < args$padj_cutoff
  hit$pathway <- factor(hit$pathway, levels = hit$pathway[order(hit$NES)])

  sub <- paste0(length(wanted) - length(absent), " of ", length(wanted),
                " gene sets tested")
  if (nzchar(args$numerator) && nzchar(args$denominator)) {
    sub <- paste0(sub, "  ·  NES > 0: enriched in ", args$numerator,
                  ", NES < 0: enriched in ", args$denominator)
  }

  p <- ggplot(hit, aes(x = NES, y = pathway, colour = padj, size = size)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(shape = significant)) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                       name = paste0("padj < ", args$padj_cutoff),
                       drop = FALSE) +
    scale_colour_gradient(low = "#E41A1C", high = "grey70",
                          limits = c(0, 0.25), oob = scales::squish,
                          name = "padj") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    labs(title = paste("KEGG:", label), subtitle = sub, x = "NES", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"))

  ggsave(out_pdf, plot = p, width = 9,
         height = max(3, nrow(hit) * 0.35 + 2.5), limitsize = FALSE)
  message("Written: ", out_pdf, "  (", nrow(hit), " gene sets)")
}

# One block at the end so the reasons are visible without scrolling back
# through the per-category output.
message("")
message("==== KEGG category summary ====")
for (s in selected) {
  wanted <- unique(cats$gs_name[cats$slug == s])
  absent <- setdiff(wanted, res$pathway)
  message(sprintf("  %-32s %d/%d plotted", s, length(wanted) - length(absent),
                  length(wanted)))
  for (gs in absent) message("       missing  ", gs, ": ", explain_missing(gs))
}
message("===============================")

message("KEGG category plots complete.")
