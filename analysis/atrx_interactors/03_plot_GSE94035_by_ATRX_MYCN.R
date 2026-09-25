#!/usr/bin/env Rscript
# GSE94035 ATRX-interactor boxplots split by patient ATRX / MYCN status:
# sample groups (default Tumor_diagnosis, DTC_diagnosis, DTC_relapse; MNCs
# left out) x 3 statuses. In overview.pdf all panels of a page (FL or IFF
# list) share one y-axis: that list's min to max, padded by 5%.
#
# Patient status comes from a clinical spreadsheet passed with --clinical
# (e.g. 20230524_TGF__Fikrets_RNAseq.xlsx). It is read in place and never
# copied: only the columns patient_o_id, atrx and mna of the "Samples" sheet
# are used. Spreadsheet patient "p0006" is GEO patient "p06".
#
# Status per patient (over all its rows in the spreadsheet):
#   ATRXdel          atrx is "Deletion"
#   ATRXwt_MYCNA     atrx only "Normal"/"NO", mna "YES"
#   ATRXwt_nonMYCNA  atrx only "Normal"/"NO", mna "NO"
#   unassigned       anything else (blank atrx, other atrx values such as
#                    "Xq loss", or patient absent) - not plotted, listed in
#                    patient_status.tsv
#
# Statistics: within each sample group, two-sided unpaired Wilcoxon rank-sum
# between each pair of statuses (NA when a status has < 2 samples); BH within
# each group x status pair across the plotted genes.
#
# Outputs (in --outdir): overview.pdf, per_gene.pdf, expression_long.csv,
# stats.csv, patient_status.tsv (GEO patient -> status only).
#
# Usage (in the pipeline's deseq2 conda env, which has readxl):
#   Rscript analysis/atrx_interactors/03_plot_GSE94035_by_ATRX_MYCN.R \
#     --clinical /path/to/20230524_TGF__Fikrets_RNAseq.xlsx

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(readxl)
})

option_list <- list(
  make_option("--clinical", type = "character", default = "",
              help = "Clinical spreadsheet with a 'Samples' sheet (patient_o_id, atrx, mna); required"),
  make_option("--sheet", type = "character", default = "Samples",
              help = "Sheet name in --clinical [default %default]"),
  make_option("--geo_dir", type = "character",
              default = "/nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret",
              help = "Directory created by 00_download_GSE94035.sh [default %default]"),
  make_option("--groups", type = "character",
              default = "Tumor_diagnosis,DTC_diagnosis,DTC_relapse",
              help = "samplesheet groups to plot, in x-axis order [default %default]"),
  make_option("--genes", type = "character", default = "",
              help = "Gene table [default atrx_interactors.tsv next to this script]"),
  make_option("--outdir", type = "character",
              default = "/nobackup/lab_taschner-mandl/arthurdondi/projects/epicode/atrx_interactors/GSE94035_Fikret/by_ATRX_MYCN",
              help = "Output directory [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))
if (!nzchar(args$clinical)) stop("--clinical is required")

here <- {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(here, "boxplot_helpers.R"))

groups   <- trimws(strsplit(args$groups, ",")[[1]])
statuses <- c("ATRXdel", "ATRXwt_MYCNA", "ATRXwt_nonMYCNA")
genes_path <- if (nzchar(args$genes)) args$genes else file.path(here, "atrx_interactors.tsv")
outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ─── Patient status from the clinical sheet (3 columns only) ─────────────────
raw <- suppressMessages(read_excel(args$clinical, sheet = args$sheet,
                                   col_names = FALSE, col_types = "text"))
hdr_row <- which(apply(raw, 1, function(r) "patient_o_id" %in% r))[1]
if (is.na(hdr_row)) stop("No 'patient_o_id' header in sheet ", args$sheet)
hdr  <- as.character(unlist(raw[hdr_row, ]))
need <- c("patient_o_id", "atrx", "mna")
if (!all(need %in% hdr)) stop("Missing columns: ", paste(setdiff(need, hdr), collapse = ", "))
clin <- as.data.frame(raw[-seq_len(hdr_row), match(need, hdr)], stringsAsFactors = FALSE)
colnames(clin) <- need
rm(raw)
clin <- clin[!is.na(clin$patient_o_id) & grepl("^p[0-9]+$", clin$patient_o_id), ]
clin$pnum <- as.integer(sub("^p", "", clin$patient_o_id))
norm <- function(x) toupper(trimws(x[!is.na(x) & nzchar(trimws(x))]))

status_of <- function(x) {
  a <- unique(norm(x$atrx)); m <- unique(norm(x$mna))
  if ("DELETION" %in% a) return("ATRXdel")
  if (length(a) && all(a %in% c("NORMAL", "NO"))) {
    if (identical(m, "YES")) return("ATRXwt_MYCNA")
    if (identical(m, "NO"))  return("ATRXwt_nonMYCNA")
  }
  "unassigned"
}
pat_status <- vapply(split(clin, clin$pnum), status_of, character(1))

# ─── GEO data ────────────────────────────────────────────────────────────────
sheet <- read.delim(file.path(args$geo_dir, "metadata", "samplesheet.tsv"),
                    stringsAsFactors = FALSE, colClasses = "character")
sheet <- sheet[sheet$group %in% groups, ]
sheet$status <- unname(pat_status[as.character(as.integer(sub("^p", "", sheet$patient)))])
sheet$status[is.na(sheet$status)] <- "unassigned"

ps <- unique(sheet[, c("patient", "status")])
ps <- ps[order(ps$status, ps$patient), ]
write.table(ps, file.path(outdir, "patient_status.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Patients per status:"); print(table(ps$status))
message("Unassigned (not plotted): ",
        paste(ps$patient[ps$status == "unassigned"], collapse = ", "))

sheet <- sheet[sheet$status %in% statuses, ]
sheet$group  <- factor(sheet$group, levels = groups)
sheet$status <- factor(sheet$status, levels = statuses)
message("Samples per group x status:"); print(table(sheet$group, sheet$status))

mat <- read.delim(gzfile(file.path(args$geo_dir, "raw", "GSE94035_Deseq_mat.txt.gz")),
                  row.names = 1, check.names = FALSE)
colnames(mat) <- sub("\r$", "", colnames(mat))
mat <- as.matrix(mat[, sheet$matrix_column, drop = FALSE])

genes <- read_gene_table(genes_path)
genes <- genes[genes$ensembl_gene_id %in% rownames(mat), , drop = FALSE]

expr <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  data.frame(gene = genes$gene[i], label = genes$label[i], list = genes$list[i],
             sample = sheet$matrix_column, patient = sheet$patient,
             group = sheet$group, status = sheet$status,
             value = as.numeric(mat[genes$ensembl_gene_id[i], ]),
             stringsAsFactors = FALSE)
}))
expr$label <- factor(expr$label, levels = genes$label)
write.csv(transform(expr, label = NULL, log2_FPM = value, value = NULL),
          file.path(outdir, "expression_long.csv"), row.names = FALSE)

# ─── Statistics: status pairs within each group ──────────────────────────────
pairs <- data.frame(status1 = c("ATRXdel", "ATRXwt_MYCNA", "ATRXdel"),
                    status2 = c("ATRXwt_MYCNA", "ATRXwt_nonMYCNA", "ATRXwt_nonMYCNA"),
                    tier    = 1:3, stringsAsFactors = FALSE)
stats <- do.call(rbind, lapply(genes$gene, function(g) {
  e <- expr[expr$gene == g, ]
  do.call(rbind, lapply(groups, function(gr) {
    eg <- e[e$group == gr, ]
    data.frame(gene = g, group = gr, status1 = pairs$status1, status2 = pairs$status2,
               tier = pairs$tier,
               n1 = sapply(pairs$status1, function(s) sum(eg$status == s)),
               n2 = sapply(pairs$status2, function(s) sum(eg$status == s)),
               log2_diff_median = mapply(function(a, b)
                 median(eg$value[eg$status == a]) - median(eg$value[eg$status == b]),
                 pairs$status1, pairs$status2),
               wilcox_p = mapply(function(a, b)
                 wilcox_p(eg$value[eg$status == a], eg$value[eg$status == b]),
                 pairs$status1, pairs$status2),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}))
stats$log2_diff_median[!is.finite(stats$log2_diff_median)] <- NA
stats$wilcox_padj_BH <- ave(stats$wilcox_p, stats$group, stats$status1, stats$status2,
                            FUN = function(p) p.adjust(p, method = "BH"))
stats <- merge(genes[, c("gene", "list", "call")], stats, by = "gene", sort = FALSE)
write.csv(stats[, setdiff(colnames(stats), "tier")], file.path(outdir, "stats.csv"),
          row.names = FALSE)

# ─── Plots: numeric x so empty statuses keep their slot ──────────────────────
dodge  <- 0.27
colours <- setNames(c("#4f4e92", "#A7563C", "#FAA21B"), statuses)
xpos <- function(group, status) as.integer(factor(group, levels = groups)) +
  (as.integer(factor(status, levels = statuses)) - 2) * dodge
set.seed(1)
expr$x  <- xpos(expr$group, expr$status)
expr$xj <- expr$x + runif(nrow(expr), -0.06, 0.06)
ylab <- expression(log[2] ~ "DESeq2-normalised FPM")
status_lab <- c(ATRXdel = "ATRXdel", ATRXwt_MYCNA = "ATRXwt MYCNA",
                ATRXwt_nonMYCNA = "ATRXwt nonMYCNA")

base_plot <- function(df, point_size) {
  ggplot(df, aes(x = x, y = value)) +
    geom_boxplot(aes(group = interaction(group, status), fill = status),
                 width = dodge * 0.85, outlier.shape = NA, alpha = 0.75) +
    geom_point(aes(x = xj), shape = 21, fill = "white", colour = "black",
               size = point_size, stroke = 0.3) +
    geom_vline(xintercept = seq(1.5, length(groups) - 0.5), colour = "grey85", linewidth = 0.3) +
    scale_fill_manual(values = colours, labels = status_lab, drop = FALSE, name = NULL) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "bottom")
}

short_lab <- c(Tumor_diagnosis = "Tum dx", DTC_diagnosis = "DTC dx", DTC_relapse = "DTC rel",
               MNC_diagnosis = "MNC dx", MNC_relapse = "MNC rel", DTC_relapse_unenriched = "DTC rel (unenr.)")
pdf(file.path(outdir, "overview.pdf"), width = 11, height = 8.5)
for (l in unique(genes$list)) {
  df <- droplevels(expr[expr$list == l, ])
  # Same y-axis for every panel of this list: min to max + 5% padding.
  ylim <- range(df$value, na.rm = TRUE)
  ylim <- ylim + c(-1, 1) * 0.05 * diff(ylim)
  p <- base_plot(df, 0.5) +
    facet_wrap(~ label, scales = "fixed", ncol = 5) +
    coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = seq_along(groups), labels = short_lab[groups],
                       expand = expansion(add = 0.4)) +
    labs(title = sprintf("ATRX %s interactors - GSE94035 by ATRX/MYCN status", l),
         x = NULL, y = ylab) +
    theme(text = element_text(size = 8), strip.text = element_text(size = 7, lineheight = 0.9),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  print(p)
}
invisible(dev.off())

n_tab <- as.data.frame(table(group = sheet$group, status = sheet$status))
n_tab$x <- xpos(n_tab$group, n_tab$status)
caption <- paste(
  "Within each sample group: Wilcoxon rank-sum, two-sided, unpaired; p = raw p (in parentheses: BH-adjusted within that group x status pair across",
  sprintf("the %d genes). No bracket when a status has < 2 samples. Status from the clinical sheet (atrx, mna); unassigned patients not shown.", nrow(genes)),
  sep = "\n")
pdf(file.path(outdir, "per_gene.pdf"), width = 3 + 2 * length(groups), height = 7.5)
for (i in seq_len(nrow(genes))) {
  g  <- genes$gene[i]
  df <- expr[expr$gene == g, ]
  st <- stats[stats$gene == g & !is.na(stats$wilcox_p), ]
  ymax <- max(df$value, na.rm = TRUE); ymin <- min(df$value, na.rm = TRUE)
  step <- max(ymax - ymin, 1) * 0.09; tick <- step * 0.2
  st$x1 <- xpos(st$group, st$status1); st$x2 <- xpos(st$group, st$status2)
  st$y  <- ymax + step * st$tier
  st$label <- sprintf("p=%s (%s)", fmt_p(st$wilcox_p), fmt_p(st$wilcox_padj_BH))
  n_df <- transform(n_tab, y = ymin - step * 0.6, lab = paste0("n=", Freq))
  p <- base_plot(df, 1.3) +
    geom_segment(data = st, inherit.aes = FALSE, aes(x = x1, xend = x2, y = y, yend = y), linewidth = 0.3) +
    geom_segment(data = st, inherit.aes = FALSE, aes(x = x1, xend = x1, y = y, yend = y - tick), linewidth = 0.3) +
    geom_segment(data = st, inherit.aes = FALSE, aes(x = x2, xend = x2, y = y, yend = y - tick), linewidth = 0.3) +
    geom_text(data = st, inherit.aes = FALSE, aes(x = (x1 + x2) / 2, y = y + tick * 0.4, label = label),
              size = 2.3, vjust = 0) +
    geom_text(data = n_df, inherit.aes = FALSE, aes(x = x, y = y, label = lab), size = 2.3, colour = "grey30") +
    scale_x_continuous(breaks = seq_along(groups), labels = groups, expand = expansion(add = 0.4)) +
    labs(title = g,
         subtitle = sprintf("ATRX %s interactor: %s (UniProt %s, %s)", genes$list[i],
                            genes$call[i], genes$uniprot[i], genes$ensembl_gene_id[i]),
         x = NULL, y = ylab, caption = caption) +
    theme(plot.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"))
  print(p)
}
invisible(dev.off())

message("Done: ", nrow(genes), " genes; outputs in ", outdir)
