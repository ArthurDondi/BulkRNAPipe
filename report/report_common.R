# Setup and helpers shared by cell_lines_report.Rmd and patients_report.Rmd.
# Sourced by each report's setup chunk into the knit environment, after it
# has set REPORT (the name of its export folders); everything here reads
# params lazily, so a report only needs the params of the helpers it calls.

# echo = TRUE + code_folding: hide (in the YAML header) is the "Show" toggle
# per chunk - code is present in the HTML but folded away by default.
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
                      dpi = 150, fig.align = "center")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(knitr)
})

OUT <- params$outdir
LEV <- params$condition_levels

# ── Per-plot PDFs and per-table CSVs ─────────────────────────────────────────
# <export_dir>/figures/<REPORT>/<NN_section>/ and
# <export_dir>/tables/<REPORT>/<NN_section>/, the section folders numbered and
# named after the report's top-level headings. REPORT is set by the report
# before it sources this file ("cell_lines" or "patients"), so each render
# empties and rewrites its own folders only.
# Plots knitr draws are saved by a second graphics device (PDF next to the PNG
# it embeds) and copied over by the last chunk; plots drawn through
# show_plot() and tables printed through kable() are written straight away.
# Files are named <chunk label>-<n> unless export_as() names them.
EXPORT <- isTRUE(params$export)
EXPORT_DIR <- if (nzchar(params$export_dir)) params$export_dir else file.path(OUT, "report")
if (!exists("REPORT")) REPORT <- "report"
FIG_DIR <- file.path(EXPORT_DIR, "figures", REPORT)
TAB_DIR <- file.path(EXPORT_DIR, "tables", REPORT)
# cairo_pdf embeds the non-ASCII characters (★, ≤, —) the plots use; pdf() would not.
PDF_DEV <- if (isTRUE(capabilities("cairo"))) "cairo_pdf" else "pdf"
if (EXPORT) {
  unlink(c(FIG_DIR, TAB_DIR), recursive = TRUE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
  knitr::opts_chunk$set(dev = c("png", PDF_DEV))
}

# Chunk label -> section folder, read from this file's own top-level headings
# (number_sections numbers them in the same order). Lines inside code chunks
# are skipped: an R comment also starts with "# ".
SECTION_OF <- local({
  slug <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))
  lines <- readLines(knitr::current_input(dir = TRUE), warn = FALSE, encoding = "UTF-8")
  in_chunk <- FALSE; sec <- "00_setup"; i <- 0; map <- character()
  for (l in lines) {
    if (grepl("^\\s*```", l)) {
      lab <- regmatches(l, regexec("^\\s*```\\{r[ ,]+([A-Za-z0-9_.-]+)\\s*[,}]", l))[[1]]
      if (!in_chunk && length(lab) == 2) map[lab[2]] <- sec
      in_chunk <- !in_chunk
      next
    }
    if (!in_chunk && grepl("^# ", l)) {
      i <- i + 1
      sec <- sprintf("%02d_%s", i, slug(sub("^# +", "", l)))
    }
  }
  map
})
cur_section <- function(label = knitr::opts_current$get("label")) {
  s <- if (is.null(label)) NA else SECTION_OF[label]
  if (is.na(s)) "unsorted" else unname(s)
}

# One-shot names for the next exported plot and/or table, e.g.
# export_as("all/ATRX_DEL_vs_MNA_volcano"); "/" makes subfolders.
.export <- new.env()
.export$n <- list()
export_as <- function(plot = NULL, table = NULL) {
  .export$plot <- plot
  .export$table <- table
}
export_path <- function(root, kind, ext, name = NULL) {
  label <- knitr::opts_current$get("label")
  if (is.null(label)) label <- "chunk"
  if (is.null(name)) {
    name <- .export[[kind]]
    .export[[kind]] <- NULL
  }
  if (is.null(name)) {
    key <- paste(kind, label)
    .export$n[[key]] <- (if (is.null(.export$n[[key]])) 0 else .export$n[[key]]) + 1
    # "-fig" keeps these apart from knitr's own <label>-<n>.pdf in the same chunk.
    name <- sprintf("%s-%s%d", label, if (kind == "plot") "fig" else "", .export$n[[key]])
  }
  f <- file.path(root, cur_section(label), paste0(name, ext))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  f
}
export_plot <- function(p, width, height, name = NULL) {
  if (!EXPORT) return(invisible())
  try(ggsave(export_path(FIG_DIR, "plot", ".pdf", name), p, width = width,
             height = height, device = get(PDF_DEV)), silent = TRUE)
}
# HTML markup (bold totals, line breaks) is stripped from the CSV, and fmt_p()'s
# non-breaking hyphens turned back into plain minus signs.
export_table <- function(x, name = NULL) {
  if (!EXPORT) return(invisible())
  d <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_along(d)) {
    if (is.list(d[[j]])) d[[j]] <- vapply(d[[j]], paste, "", collapse = ";")
    if (is.character(d[[j]]) || is.factor(d[[j]]))
      d[[j]] <- gsub("\u2011", "-", gsub("<[^>]+>", "", as.character(d[[j]])), fixed = TRUE)
  }
  rn <- !identical(rownames(d), as.character(seq_len(nrow(d))))
  try(write.csv(d, export_path(TAB_DIR, "table", ".csv", name), row.names = rn),
      silent = TRUE)
}
# Every table in the report goes through kable(), so this wrapper exports each
# one as it is printed.
kable <- function(x, ...) {
  export_table(x)
  knitr::kable(x, ...)
}

# Condition palette from the lab colour sheet (20261021_Colors.pdf): one fixed
# colour per line, keyed by name so editing condition_levels cannot reshuffle
# them (1760 = EmptyVector, 1011 = ATRX_FL, 1744 = ATRX_IFF, see
# config_epicode.yaml). Conditions not on the sheet fall back to grey. Keep in
# sync with condition_colours in analysis/atrx_interactors/boxplot_helpers.R.
COND_COLOURS <- c(TP53 = "#ABAB50", E6 = "#46184E", EmptyVector = "#D16530",
                  ATRX_FL = "#19468C", ATRX_IFF = "#B72656")
PAL <- setNames(ifelse(LEV %in% names(COND_COLOURS), COND_COLOURS[LEV], "grey60"), LEV)

# Volcano colours from the same sheet: significant genes green on the left
# (log2FC < 0), pink on the right (log2FC > 0). Non-significant genes are a
# light grey, not grey70: grey70 is indistinguishable from the pink under
# protanopia. The sheet's spare blue marks highlighted genes and proteomics hits.
VOLC_COL <- c(down = "#4C5021", up = "#E5A6B7", ns = "#E6E6E6")
HL_COL   <- "#6D8AB6"

# Enrichment plots (GSEA, GO, over-representation): padj on a single-hue blue
# ramp, dark where padj is lowest. Neither the old red ramp nor viridis, which
# the proteomics analysis uses, so transcriptomic and proteomic enrichment
# panels never share a colour scale. log10 scale from `floor` (anything lower,
# including 0, gets the darkest colour) to the cutoff.
scale_colour_padj <- function(cut = params$padj, name = "padj", floor = 1e-10) {
  # Every 2 decades up to 1e-4; 1e-2 would sit on top of the 0.05 cutoff.
  brk <- c(floor, 10^-(seq(-log10(floor) - 2, 4, by = -2)), cut)
  brk <- sort(unique(brk[brk >= floor & brk <= cut]))
  scale_colour_gradient(
    low = "#0d366b", high = "#86b6ef", trans = "log10", limits = c(floor, cut),
    breaks = brk,
    labels = ifelse(brk == floor, paste("\u2264", format(floor)), vapply(brk, format, "")),
    oob = function(x, range) scales::squish(x, range, only.finite = FALSE), name = name)
}

# Whole-number size legend (counts of genes): about 3 pretty breaks, rounded.
int_breaks <- function(lim) {
  b <- unique(round(pretty(lim, n = 3)))
  b[b >= lim[1] & b <= lim[2]]
}

# Draw a plot as an image of the given size (inches). scroll = TRUE puts it in
# a horizontally scrollable box at full width, so long pathway names do not
# squeeze the panel.
show_plot <- function(p, width, height, scroll = FALSE, name = NULL) {
  if (!knitr::is_html_output()) {  # knitr's devices export it
    if (inherits(p, "gtable")) { grid::grid.newpage(); grid::grid.draw(p); return(invisible()) }
    return(print(p))
  }
  export_plot(p, width, height, name)
  # The PNG is referenced by path and pandoc embeds it (self_contained), rather
  # than inlined here as base64: pandoc parsing ~100 inlined images in raw HTML
  # took ~7x the memory (3.5 vs 0.5 GB on a 16 MB page) and got the cluster
  # job killed. The file lives in R's session tempdir, which outlasts pandoc.
  f <- tempfile(fileext = ".png")
  ggsave(f, p, width = width, height = height, dpi = 150)
  # Tag built with paste0, not a sprintf template: rmarkdown scans the whole
  # page, folded source code included, for img src attributes and warns on a
  # literal "%s" there.
  img <- paste0('<img src="', f, '" style="width: ', round(width * 96),
                'px; max-width: ', if (scroll) "none" else "100%", ';">')
  cat("\n", if (scroll) sprintf('<div style="overflow-x: auto;">%s</div>', img) else img, "\n\n")
}

# Enrichment dot plot sized to its labels and number of rows.
# The dot panel always gets panel_width inches: everything around it (pathway
# names, facet strips, legends) is measured on the built plot and added on top,
# so long names widen the image, and the scroll box, instead of squeezing the
# panel. 8.4 in (~810 px) fits the report column; only wider plots scroll. At
# least 4 in tall so the legends fit; titles start at the left edge, not over
# the panel, which long labels push far right. `labels` is kept for callers
# but no longer used: the measurement covers it.
show_enrich <- function(p, labels, n_rows, panel_width = 3.5, name = NULL) {
  p <- p + theme(plot.title.position = "plot")
  # Text widths need a device: a null PDF, closed straight away, then the
  # device knitr draws on made current again (dev.off() would otherwise pick
  # whichever device comes next).
  prev <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  g <- ggplotGrob(p)
  panel_cols <- unique(unlist(g$layout[grepl("^panel", g$layout$name), c("l", "r")]))
  around <- grid::convertWidth(sum(g$widths[-panel_cols]), "in", valueOnly = TRUE)
  grDevices::dev.off(dev)
  if (prev > 1) grDevices::dev.set(prev)
  show_plot(p, width = max(8.4, around + panel_width),
            height = max(4, 1.8 + 0.26 * n_rows), scroll = TRUE, name = name)
}

# A table in a scrollable box (both directions), under the plot it lists.
scroll_table <- function(d, max_height = 300) {
  cat(sprintf('<div style="max-height: %dpx; overflow: auto; margin-bottom: 1.5em;">\n',
              max_height))
  print(kable(d, format = "html", row.names = FALSE,
              # "enrich": see the <style> block after the setup chunk.
              table.attr = 'class="table table-condensed enrich" style="font-size: 85%;"'))
  cat("\n</div>\n\n")
}

# p-values as text (2 significant digits): kable would otherwise reformat a
# numeric column and print tiny values as 0.00. Non-breaking hyphen so
# "9.9e-24" does not wrap in a narrow column.
fmt_p <- function(x) gsub("-", "\u2011", formatC(x, format = "g", digits = 2), fixed = TRUE)

# Gene lists as written by fgsea (";") or shared_genes.R ("/"), comma-separated.
genes_list <- function(x) if (is.null(x)) "" else gsub("[;/]", ", ", x)

# The table under a GSEA dot plot: the plotted pathways, top to bottom.
gsea_table <- function(df) {
  df <- df[order(-df$NES), , drop = FALSE]
  scroll_table(data.frame(pathway = as.character(df$pathway), NES = round(df$NES, 2),
                          padj = fmt_p(df$padj), `set size` = df$size,
                          `leading-edge genes` = genes_list(df$leadingEdge),
                          check.names = FALSE))
}

theme_set(theme_bw(base_size = 11))
theme_update(panel.grid.minor = element_blank(),
             plot.title = element_text(face = "bold", size = 12),
             plot.subtitle = element_text(colour = "grey30", size = 9))

# ── Helpers ──────────────────────────────────────────────────────────────────
# Read a pipeline output, or return NULL and leave a visible note in the report.
grab <- function(...) {
  f <- file.path(OUT, ...)
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
missing_note <- function(...) {
  cat("\n> **Not available.** `", file.path(...), "` was not found. ",
      "Run that part of the pipeline and re-render.\n\n", sep = "")
}
cond_factor <- function(x) factor(x, levels = LEV[LEV %in% unique(x)])

# Long-format helper for the wide per-gene tables (gene_id + one column per sample).
gene_long <- function(df, meta) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>%
    left_join(meta, by = "sample") %>%
    mutate(condition = cond_factor(condition),
           sample = factor(sample, levels = meta$sample))
}

# Sample sheet, used to colour everything by condition.
qc <- grab("design_qc", "design_qc_summary.csv")
meta <- if (!is.null(qc)) {
  qc %>% select(sample, condition) %>%
    arrange(match(condition, LEV)) %>%
    mutate(condition = cond_factor(condition))
} else NULL

# DESeq results with the vector/transgene features dropped. They are guaranteed
# to be the most significant genes in any contrast that crosses a construct
# boundary, and they say nothing about the cells' biology.
de_results <- function(...) {
  r <- grab(...)
  if (is.null(r)) return(NULL)
  n0 <- nrow(r)
  r  <- r[!tolower(r$gene_id) %in% tolower(params$exclude_genes), , drop = FALSE]
  attr(r, "n_dropped") <- n0 - nrow(r)
  r
}

# Volcano dots, shared by volcano() and volcano_with_proteomic(): the
# non-significant genes first, then the significant ones on top, coloured by
# side (VOLC_COL: green left, pink right). res needs a logical `sig` column.
volcano_dots <- function(res) {
  list(
    geom_point(data = filter(res, !sig), colour = VOLC_COL[["ns"]], size = 0.8),
    geom_point(aes(colour = ifelse(log2FoldChange > 0, "up", "down")),
               data = filter(res, sig), size = 1, alpha = 0.9),
    scale_colour_manual(values = VOLC_COL, guide = "none"))
}

volcano <- function(res, title, num, den, highlight = NULL,
                    label_top = params$volcano_label_top) {
  res <- res %>% filter(!is.na(padj)) %>%
    mutate(sig = padj < params$padj & abs(log2FoldChange) >= params$lfc)
  n <- c(up = sum(res$sig & res$log2FoldChange > 0),
         down = sum(res$sig & res$log2FoldChange < 0))
  p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
    volcano_dots(res) +
    geom_vline(xintercept = c(-params$lfc, params$lfc), linetype = "dashed",
               linewidth = 0.3) +
    geom_hline(yintercept = -log10(params$padj), linetype = "dashed",
               linewidth = 0.3) +
    labs(
      title = title,
      # Spell the direction out on both axes: "up in X" is ambiguous unless the
      # reader already knows which group is the numerator.
      subtitle = sprintf(
        paste0("%s of %s genes tested significant (padj < %s, |log2FC| >= %s)\n",
               "RIGHT, log2FC > 0: %s genes HIGHER in %s   |   ",
               "LEFT, log2FC < 0: %s genes HIGHER in %s"),
        format(sum(res$sig), big.mark = ","), format(nrow(res), big.mark = ","),
        params$padj, params$lfc,
        format(n["up"], big.mark = ","), num,
        format(n["down"], big.mark = ","), den),
      x = sprintf("log2FC  (%s / %s)   \u2014   positive = higher in %s", num, den, num),
      y = expression(-log[10] ~ "adjusted p")) +
    theme(plot.subtitle = element_text(lineheight = 1.2))
  # Label the most significant genes, ranked the same way as top_genes():
  # -log10(padj) x |log2FC|, so the labelled genes match the accompanying table.
  # The pipeline's own volcano.pdf labels every significant gene; capped here
  # because report-scale figures become unreadable once a contrast has more
  # than a few dozen significant genes (E6_vs_TP53 has thousands).
  hl_lower <- if (is.null(highlight)) character(0) else tolower(highlight)
  if (label_top > 0) {
    lbl <- res %>% filter(sig, !tolower(gene_id) %in% hl_lower) %>%
      mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
      slice_max(score, n = label_top, with_ties = FALSE)
    if (nrow(lbl)) {
      p <- p + ggrepel::geom_text_repel(
        data = lbl, aes(label = gene_id), colour = "grey20", size = 2.4,
        max.overlaps = 30, segment.size = 0.2, segment.colour = "grey50",
        show.legend = FALSE)
    }
  }
  if (!is.null(highlight)) {
    h <- res %>% filter(tolower(gene_id) %in% hl_lower)
    if (nrow(h)) {
      # Drawn after the general labels so the highlighted genes and their
      # labels sit on top rather than being obscured by them.
      # Blue point, dark label: green now means "significant, right side".
      p <- p +
        geom_point(data = h, shape = 21, fill = HL_COL, colour = "grey15",
                   size = 2.6, stroke = 0.4) +
        ggrepel::geom_text_repel(
          data = h, aes(label = gene_id), colour = "grey10", fontface = "bold",
          size = 3.1, min.segment.length = 0, box.padding = 0.6,
          segment.colour = HL_COL, max.overlaps = Inf, show.legend = FALSE)
    }
    miss <- setdiff(hl_lower, tolower(res$gene_id))
    if (length(miss)) message("highlight not found: ", paste(miss, collapse = ", "))
  }
  p
}

# Ranking by padj alone floats high-count genes with tiny fold-changes to the
# top. This ranks significant genes by -log10(padj) * |log2FC| so both matter,
# and reports the two directions separately so one cannot crowd out the other.
top_genes <- function(res, n = 10) {
  res %>% filter(!is.na(padj), padj < params$padj,
                 abs(log2FoldChange) >= params$lfc) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down"),
           score = -log10(padj) * abs(log2FoldChange)) %>%
    group_by(direction) %>% slice_max(score, n = n) %>% ungroup() %>%
    arrange(direction, desc(score)) %>%
    transmute(gene_id, direction, baseMean = round(baseMean),
              log2FC = round(log2FoldChange, 2), padj = signif(padj, 3),
              score = round(score, 1))
}

de_counts <- function(res) {
  r <- res %>% filter(!is.na(padj), padj < params$padj,
                      abs(log2FoldChange) >= params$lfc)
  c(up = sum(r$log2FoldChange > 0), down = sum(r$log2FoldChange < 0))
}

# Proteomics-significant genes for one comparison, wide-format limma Excel.
# Mirrors workflow/scripts/deseq2_proteomic.R. Returns NULL (and lets the
# caller show missing_note()) if the file, sheet or columns are not there.
load_proteomics_sig <- function() {
  if (!nchar(params$proteomics_xlsx) || !file.exists(params$proteomics_xlsx)) return(NULL)
  logfc_col <- paste0(params$proteomics_comparison, params$proteomics_logfc_column)
  fdr_col   <- paste0(params$proteomics_comparison, params$proteomics_fdr_column)
  tbl <- tryCatch(
    readxl::read_excel(params$proteomics_xlsx, sheet = params$proteomics_sheet, col_types = "text"),
    error = function(e) NULL)
  if (is.null(tbl) || !all(c(params$proteomics_gene_column, logfc_col, fdr_col) %in% colnames(tbl)))
    return(NULL)
  tbl %>%
    transmute(gene_id    = trimws(as.character(.data[[params$proteomics_gene_column]])),
              prot_fdr   = suppressWarnings(as.numeric(.data[[fdr_col]])),
              prot_logfc = suppressWarnings(as.numeric(.data[[logfc_col]]))) %>%
    filter(!is.na(gene_id), gene_id != "", !is.na(prot_fdr), !is.na(prot_logfc),
           prot_fdr <= params$proteomics_fdr_threshold,
           abs(prot_logfc) > params$proteomics_lfc_threshold) %>%
    arrange(prot_fdr, dplyr::desc(abs(prot_logfc))) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    transmute(gene_id, prot_direction = ifelse(prot_logfc > 0, "Up", "Down"))
}

# RNA volcano restricted to genes that are also proteomics-significant - same
# construction as the pipeline's own volcano_proteomic.pdf, redrawn here in the
# report's own style. RNA-significant genes are coloured by side as in
# volcano() (green left, pink right); their shape says whether RNA and protein
# move the same way.
proteomic_volcano <- function(res, title, num, den, prot_sig) {
  res <- res %>% filter(!is.na(padj)) %>%
    mutate(rna_significant = padj < params$padj & abs(log2FoldChange) >= params$lfc,
           rna_direction = ifelse(log2FoldChange > 0, "Up", "Down"))
  df <- res %>% inner_join(prot_sig, by = "gene_id") %>%
    mutate(side = ifelse(rna_direction == "Up", "up", "down"),
           agreement = ifelse(rna_direction == prot_direction,
                              "RNA and protein agree", "RNA and protein disagree"),
           label = ifelse(rna_significant, gene_id, NA_character_))
  if (!nrow(df)) return(invisible(NULL))
  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(data = filter(df, !rna_significant), shape = 21, size = 1.4, stroke = 0.3,
               fill = VOLC_COL[["ns"]], colour = "grey50") +
    geom_point(data = filter(df, rna_significant), aes(fill = side, shape = agreement),
               size = 2.2, stroke = 0.3, colour = "grey20") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, colour = "grey15",
                             max.overlaps = Inf, show.legend = FALSE) +
    scale_fill_manual(values = VOLC_COL, guide = "none") +
    scale_shape_manual(values = c("RNA and protein agree" = 21,
                                  "RNA and protein disagree" = 24), drop = FALSE) +
    guides(shape = guide_legend(override.aes = list(fill = "grey60", size = 2.5))) +
    geom_vline(xintercept = c(-params$lfc, params$lfc), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(params$padj), linetype = "dashed", linewidth = 0.3) +
    labs(
      title = title,
      subtitle = sprintf(
        paste0("Proteomics-significant genes only (FDR <= %s, |log2FC| > %s)\n",
               "RIGHT, log2FC > 0: higher in %s   |   LEFT, log2FC < 0: higher in %s\n",
               "RNA-significant genes coloured by side; shape = RNA/protein agreement"),
        params$proteomics_fdr_threshold, params$proteomics_lfc_threshold, num, den),
      x = sprintf("log2FC  (%s / %s)   —   positive = higher in %s", num, den, num),
      y = expression(-log[10] ~ "adjusted p"), shape = NULL) +
    theme(plot.subtitle = element_text(lineheight = 1.2, size = 8.5), legend.position = "bottom")
}

# The full RNA volcano (volcano()'s own dots, every gene, and the same top-N
# gene labels), with proteomics-significant genes redrawn as stars on top: blue
# if also RNA-significant, black otherwise. The blue stars get their own label
# (a gene already starred blue is skipped from the general top-N labels below,
# so it is never labelled twice).
volcano_with_proteomic <- function(res, title, num, den, prot_sig,
                                   label_top = params$volcano_label_top) {
  res <- res %>% filter(!is.na(padj)) %>%
    mutate(sig = padj < params$padj & abs(log2FoldChange) >= params$lfc)
  n <- c(up = sum(res$sig & res$log2FoldChange > 0),
         down = sum(res$sig & res$log2FoldChange < 0))
  stars <- res %>% inner_join(prot_sig, by = "gene_id")
  stars_sig    <- stars %>% filter(sig)
  stars_nonsig <- stars %>% filter(!sig)
  p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
    volcano_dots(res) +
    geom_vline(xintercept = c(-params$lfc, params$lfc), linetype = "dashed",
               linewidth = 0.3) +
    geom_hline(yintercept = -log10(params$padj), linetype = "dashed",
               linewidth = 0.3) +
    labs(
      title = title,
      subtitle = sprintf(
        paste0("%s of %s genes tested significant (padj < %s, |log2FC| >= %s)\n",
               "RIGHT, log2FC > 0: %s genes HIGHER in %s   |   ",
               "LEFT, log2FC < 0: %s genes HIGHER in %s\n",
               "★ proteomics-significant (FDR <= %s, |log2FC| > %s): ",
               "blue = also RNA-significant, black = not"),
        format(sum(res$sig), big.mark = ","), format(nrow(res), big.mark = ","),
        params$padj, params$lfc,
        format(n["up"], big.mark = ","), num,
        format(n["down"], big.mark = ","), den,
        params$proteomics_fdr_threshold, params$proteomics_lfc_threshold),
      x = sprintf("log2FC  (%s / %s)   —   positive = higher in %s", num, den, num),
      y = expression(-log[10] ~ "adjusted p")) +
    theme(plot.subtitle = element_text(lineheight = 1.2, size = 8.5))
  # Same top-N labelling as volcano(): ranked by -log10(padj) x |log2FC|.
  # Genes already labelled blue below are skipped here to avoid a duplicate.
  star_sig_ids <- tolower(stars_sig$gene_id)
  if (label_top > 0) {
    lbl <- res %>% filter(sig, !tolower(gene_id) %in% star_sig_ids) %>%
      mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
      slice_max(score, n = label_top, with_ties = FALSE)
    if (nrow(lbl)) {
      p <- p + ggrepel::geom_text_repel(
        data = lbl, aes(label = gene_id), colour = "grey20", size = 2.4,
        max.overlaps = 30, segment.size = 0.2, segment.colour = "grey50",
        show.legend = FALSE)
    }
  }
  if (nrow(stars_nonsig)) {
    p <- p + geom_point(data = stars_nonsig, shape = 8, size = 2, colour = "black")
  }
  if (nrow(stars_sig)) {
    # bg.color/bg.r draw a white halo behind the text, which keeps a label
    # readable where it crosses the coloured significant dots.
    p <- p + geom_point(data = stars_sig, shape = 8, size = 2.2, stroke = 0.8,
                        colour = HL_COL) +
      ggrepel::geom_text_repel(
        data = stars_sig, aes(label = gene_id), colour = "grey10",
        fontface = "bold", size = 2.8, min.segment.length = 0, box.padding = 0.4,
        segment.colour = HL_COL, max.overlaps = Inf, show.legend = FALSE,
        bg.color = "white", bg.r = 0.15)
  }
  p
}

# GSEA dot plot, shared by the cell-line GSEA section and the patient analyses.
# Significant pathways only (padj < params$padj), the top N by |NES| per side.
gsea_dot <- function(df, title, top = 12) {
  df <- df %>% filter(!is.na(padj), padj < params$padj) %>% arrange(desc(abs(NES)))
  df <- bind_rows(head(filter(df, NES > 0), top), head(filter(df, NES < 0), top))
  if (!nrow(df)) {
    cat(sprintf("\nNo pathway with padj < %s.\n\n", params$padj))
    return(invisible(NULL))
  }
  num <- df$numerator[1]; den <- df$denominator[1]
  up   <- sprintf("NES > 0: pathway genes HIGHER in %s", num)
  down <- sprintf("NES < 0: pathway genes HIGHER in %s", den)
  df <- df %>% mutate(
    pathway = gsub("^HALLMARK_", "", pathway),
    pathway = factor(pathway, levels = pathway[order(df$NES)]),
    side = factor(ifelse(NES > 0, up, down), levels = c(up, down)))
  p <- ggplot(df, aes(NES, pathway, colour = padj, size = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    # Strip labels wrapped onto two lines: on one, the strip takes more
    # width than the dots.
    facet_grid(side ~ ., scales = "free_y", space = "free_y",
               labeller = label_wrap_gen(width = 22)) +
    scale_colour_padj() +
    labs(title = title,
         subtitle = sprintf("enrichment of the %s / %s ranking", num, den),
         y = NULL, colour = "padj", size = "Set size") +
    theme(strip.text.y = element_text(angle = 0, size = 8))
  show_enrich(p, df$pathway, nrow(df) + 2)
  gsea_table(df)
}
