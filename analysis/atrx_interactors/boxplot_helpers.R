# Helpers for 01_plot_epicode.R.
# Sourced, not run directly.

suppressPackageStartupMessages({
  library(ggplot2)
})

# Directory of the calling Rscript, so the gene table next to it is found
# regardless of the working directory.
script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 0) return(getwd())
  dirname(normalizePath(f[1]))
}

# Gene table: gene, uniprot, list (FL/IFF), call (Yes/Maybe?),
# ensembl_gene_id, previous_symbols (";"-separated, may be empty).
read_gene_table <- function(path) {
  g <- read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  need <- c("gene", "uniprot", "list", "call", "ensembl_gene_id", "previous_symbols")
  miss <- setdiff(need, colnames(g))
  if (length(miss) > 0) stop("Gene table is missing columns: ", paste(miss, collapse = ", "))
  g$label <- sprintf("%s\n%s interactor: %s", g$gene, g$list, g$call)
  g
}

# Two-sided Wilcoxon rank-sum test; NA when a group has < 2 values.
# With n = 3 vs 3 the smallest attainable p-value is 0.1.
wilcox_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y)$p.value)
}

fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 2))

# One boxplot page for one gene.
#   df     : data.frame(group, value) - group is a factor giving the x order
#   pairs  : data.frame(group1, group2, label) - one bracket per row, may be empty
gene_page <- function(df, title, subtitle, ylab, colours, pairs = NULL, caption = NULL) {
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    median_line(width = 0.6) +
    geom_jitter(width = 0.15, height = 0, size = 1.8, shape = 21,
                colour = "black", fill = "white", stroke = 0.4) +
    scale_fill_manual(values = colours, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab, caption = caption) +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          plot.title   = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))

  if (!is.null(pairs) && nrow(pairs) > 0) {
    lv   <- levels(df$group)
    ymax <- max(df$value, na.rm = TRUE)
    ymin <- min(df$value, na.rm = TRUE)
    step <- max(ymax - ymin, 1) * 0.12
    pairs$x1 <- match(pairs$group1, lv)
    pairs$x2 <- match(pairs$group2, lv)
    pairs    <- pairs[!is.na(pairs$x1) & !is.na(pairs$x2), , drop = FALSE]
    # Shorter brackets go lower so nested ones do not cross.
    pairs    <- pairs[order(abs(pairs$x2 - pairs$x1)), , drop = FALSE]
    pairs$y  <- ymax + step * seq_len(nrow(pairs))
    tick     <- step * 0.2
    p <- p +
      geom_segment(data = pairs, inherit.aes = FALSE,
                   aes(x = x1, xend = x2, y = y, yend = y), linewidth = 0.3) +
      geom_segment(data = pairs, inherit.aes = FALSE,
                   aes(x = x1, xend = x1, y = y, yend = y - tick), linewidth = 0.3) +
      geom_segment(data = pairs, inherit.aes = FALSE,
                   aes(x = x2, xend = x2, y = y, yend = y - tick), linewidth = 0.3) +
      geom_text(data = pairs, inherit.aes = FALSE,
                aes(x = (x1 + x2) / 2, y = y + tick * 0.4, label = label),
                size = 2.6, vjust = 0, lineheight = 0.85)
  }
  p
}

# Faceted overview: one small panel per gene, no statistics.
#   df : data.frame(label, group, value) - label is the facet strip text
overview_page <- function(df, title, ylab, colours, ncol = 5, subtitle = NULL) {
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    median_line(width = 0.6) +
    geom_jitter(width = 0.15, height = 0, size = 0.8, shape = 21,
                colour = "black", fill = "white", stroke = 0.3) +
    facet_wrap(~ label, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = colours, drop = FALSE, name = NULL) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_bw(base_size = 9) +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "bottom",
          strip.text      = element_text(size = 7.5, lineheight = 0.9),
          plot.title      = element_text(face = "bold"))
}

# Condition colours from the lab colour sheet (20261021_Colors.pdf), the same
# as COND_COLOURS in report/analysis_report.Rmd - keep the two in sync.
condition_colours <- c(TP53 = "#ABAB50", E6 = "#46184E", EmptyVector = "#D16530",
                       ATRX_FL = "#19468C", ATRX_IFF = "#B72656")

# Median redrawn over a box: white where the fill is dark (relative luminance
# < 0.179, where white contrasts more than black), else the default box grey.
# The default median vanishes on the E6 and ATRX_FL colours.
median_ink <- function(fill) {
  v <- grDevices::col2rgb(fill) / 255
  v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  ifelse(colSums(v * c(0.2126, 0.7152, 0.0722)) < 0.179, "white", "grey20")
}
median_line <- function(width) {
  stat_summary(aes(colour = after_scale(median_ink(fill))), fun = median,
               fun.min = median, fun.max = median, geom = "errorbar",
               width = width, linewidth = 0.5, show.legend = FALSE)
}
