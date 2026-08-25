#!/usr/bin/env Rscript
# gsea_go_categories.R — GO term panels grouped by keyword, mirroring the KEGG
# BRITE category panels.
#
# GO has no BRITE-style categories, so each group is a regular expression matched
# against the GO term NAME. That is a text search over names, not a traversal of
# the ontology: it pulls in terms that merely share a word and misses terms that
# describe the same biology differently. Every panel says so in its subtitle, and
# the log reports how many terms matched, how many were significant, and how many
# were dropped to fit the panel.
#
# Unlike the KEGG panels -- where a category is a short fixed list and showing all
# of it is the point -- a keyword can match hundreds of GO terms, so these are
# truncated. The truncation is always stated rather than silently applied.
#
# Outputs (in --outdir):
#   {slug}_dotplot.pdf   top terms for the category
#   {slug}_results.csv   ALL matched terms, untruncated

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--results", type = "character",
              help = "GSEA results CSV for the GO collection"),
  make_option("--categories", type = "character",
              help = "TSV: slug, label, include, exclude"),
  make_option("--select", type = "character", default = "",
              help = "Comma-separated slugs to plot; empty = every slug in the TSV"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--top_n", type = "integer", default = 25L,
              help = "Maximum terms drawn per panel, ranked by padj then |NES|"),
  make_option("--numerator", type = "character", default = ""),
  make_option("--denominator", type = "character", default = ""),
  make_option("--padj_cutoff", type = "double", default = 0.05),
  make_option("--label_chars", type = "integer", default = 60L,
              help = "Truncate term labels to this many characters")
)
args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

cats <- read.delim(args$categories, header = TRUE, sep = "\t",
                   comment.char = "#", stringsAsFactors = FALSE)
need <- c("slug", "label", "include")
if (!all(need %in% colnames(cats))) {
  stop("Category file must have columns: ", paste(need, collapse = ", "))
}
if (!"exclude" %in% colnames(cats)) cats$exclude <- "-"
cats <- cats[nzchar(trimws(cats$slug)), , drop = FALSE]

selected <- trimws(strsplit(args$select, ",")[[1]])
selected <- selected[nzchar(selected)]
if (length(selected) == 0) selected <- unique(cats$slug)

for (s in setdiff(selected, unique(cats$slug))) {
  message("WARNING: category '", s, "' is not in ", args$categories,
          " - an empty panel will be written.")
}

res <- read.csv(args$results, stringsAsFactors = FALSE)
message("Loaded ", nrow(res), " tested GO terms from ", args$results)

empty_panel <- function(pdf_path, msg) {
  pdf(pdf_path, width = 8, height = 2.2)
  plot.new(); text(0.5, 0.5, msg, cex = 1.0)
  dev.off()
  message("Written (empty): ", pdf_path)
}

summary_rows <- list()

for (s in selected) {
  row   <- cats[cats$slug == s, , drop = FALSE]
  label <- if (nrow(row) > 0) row$label[1] else s
  inc   <- if (nrow(row) > 0) row$include[1] else ""
  exc   <- if (nrow(row) > 0) row$exclude[1] else "-"

  out_pdf <- file.path(args$outdir, paste0(s, "_dotplot.pdf"))
  out_csv <- file.path(args$outdir, paste0(s, "_results.csv"))

  if (!nzchar(inc)) {
    write.csv(res[0, , drop = FALSE], out_csv, row.names = FALSE)
    empty_panel(out_pdf, paste0(label, ": no include pattern"))
    summary_rows[[s]] <- c(matched = 0, signif = 0, shown = 0)
    next
  }

  keep <- grepl(inc, res$pathway, ignore.case = TRUE, perl = TRUE)
  if (nzchar(exc) && !identical(trimws(exc), "-")) {
    keep <- keep & !grepl(exc, res$pathway, ignore.case = TRUE, perl = TRUE)
  }
  hit <- res[keep, , drop = FALSE]

  n_sig <- sum(!is.na(hit$padj) & hit$padj < args$padj_cutoff)
  message(sprintf("[%s] %d terms matched, %d significant (padj < %g)",
                  s, nrow(hit), n_sig, args$padj_cutoff))

  # The CSV keeps everything; only the figure is truncated.
  write.csv(hit[order(hit$padj, -abs(hit$NES)), , drop = FALSE], out_csv,
            row.names = FALSE)

  if (nrow(hit) == 0) {
    empty_panel(out_pdf, paste0(label, ": no GO terms matched"))
    summary_rows[[s]] <- c(matched = 0, signif = 0, shown = 0)
    next
  }

  hit <- hit[order(hit$padj, -abs(hit$NES)), , drop = FALSE]
  n_matched <- nrow(hit)
  if (n_matched > args$top_n) {
    message(sprintf(
      "         plotting the top %d by padj; the other %d matched terms are in %s but not on the figure.",
      args$top_n, n_matched - args$top_n, basename(out_csv)))
    hit <- hit[seq_len(args$top_n), , drop = FALSE]
  }
  summary_rows[[s]] <- c(matched = n_matched, signif = n_sig, shown = nrow(hit))

  hit$significant <- !is.na(hit$padj) & hit$padj < args$padj_cutoff
  hit$label <- ifelse(nchar(hit$pathway) > args$label_chars,
                      paste0(substr(hit$pathway, 1, args$label_chars - 1), "…"),
                      hit$pathway)
  # Truncation can collide two long term names into the same label; make.unique
  # keeps them as separate rows instead of silently dropping one to a factor level.
  hit$label <- make.unique(hit$label)
  hit$label <- factor(hit$label, levels = hit$label[order(hit$NES)])

  sub <- sprintf("keyword match on term names, not an ontology grouping  ·  %s",
                 if (n_matched > nrow(hit))
                   sprintf("top %d of %d matched terms by padj", nrow(hit), n_matched)
                 else sprintf("all %d matched terms", n_matched))
  sub <- paste0(sub, "  \u00b7  hollow: padj \u2265 ", args$padj_cutoff)
  if (nzchar(args$numerator) && nzchar(args$denominator)) {
    sub <- paste0(sub, "\nNES > 0: enriched in ", args$numerator,
                  ", NES < 0: enriched in ", args$denominator)
  }

  # See gsea_kegg.R: colour ramp bounded at the cutoff, non-significant terms
  # drawn hollow and uncoloured, size breaks taken from the plotted points.
  sig_pts <- hit[hit$significant, , drop = FALSE]
  ns_pts  <- hit[!hit$significant, , drop = FALSE]

  sz <- range(hit$size, na.rm = TRUE)
  size_breaks <- if (sz[1] == sz[2]) sz[1] else
    unique(round(c(sz[1], mean(sz), sz[2])))

  p <- ggplot(hit, aes(x = NES, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(data = ns_pts, aes(size = size),
               shape = 21, colour = "grey65", fill = NA, stroke = 0.7) +
    geom_point(data = sig_pts, aes(size = size, colour = padj), shape = 19) +
    scale_colour_gradient(low = "#A50F15", high = "#FC9272",
                          limits = c(0, args$padj_cutoff),
                          breaks = c(0, args$padj_cutoff / 2, args$padj_cutoff),
                          oob = scales::squish, name = "padj") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8),
                          limits = sz, breaks = size_breaks) +
    labs(title = paste("GO:", label), subtitle = sub, x = "NES", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle       = element_text(size = 9, colour = "grey30"),
          legend.justification = "top")

  # Height floors on the legend column, which needs ~6in whatever the row count.
  ggsave(out_pdf, plot = p, width = 12,
         height = max(6, nrow(hit) * 0.32 + 3), limitsize = FALSE)
  message("Written: ", out_pdf, "  (", nrow(hit), " terms)")
}

message("")
message("==== GO category summary ====")
message("  (keyword matches on GO term names - not an ontology grouping)")
for (s in names(summary_rows)) {
  v <- summary_rows[[s]]
  message(sprintf("  %-32s %4d matched, %4d significant, %3d plotted",
                  s, v[["matched"]], v[["signif"]], v[["shown"]]))
}
message("=============================")
message("GO category plots complete.")
