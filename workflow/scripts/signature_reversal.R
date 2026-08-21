#!/usr/bin/env Rscript
# signature_reversal.R — score a rescue by how far it reverses a signature.
#
# The question "did the rescue work?" is usually asked as "do the rescue samples
# cluster with the wild-type samples?".  When the two sit in different genetic
# backgrounds that question is unanswerable: the background difference is larger
# than the effect and cannot be separated from it.
#
# This script asks a question that IS answerable from the same data.  Take two
# contrasts that each live entirely inside one background:
#
#   signature_contrast   the perturbation      e.g. E6 vs TP53          (the KO)
#   response_contrast    the intervention      e.g. ATRX_FL vs EmptyVector
#
# On genes that responded to the perturbation, a working intervention pushes
# expression back the other way.  Plotting one log2 fold-change against the
# other, that is a negative slope.  Neither contrast crosses the background
# boundary, so the background effect never enters the comparison and there is
# nothing to correct for.
#
# Reported statistics:
#   slope           orthogonal (total least squares) regression of response on
#                   signature; -1 = complete reversal, 0 = no response,
#                   +1 = the intervention reproduces the perturbation
#   reversal_score  -slope, so bigger is a better rescue
#   rho / r         Spearman and Pearson correlation
#   pct_reversed    % of signature genes moving the opposite way, with a
#                   binomial test against the 50% expected under no effect
#   pct_restored    % of signature genes whose response recovers at least
#                   `restored_fraction` of the perturbation, in the right direction

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(optparse)
})

option_list <- list(
  make_option("--results_signature", type = "character", help = "DESeq2 results.csv for the signature contrast"),
  make_option("--results_response",  type = "character", help = "DESeq2 results.csv for the response contrast"),
  make_option("--outdir",            type = "character", help = "Output directory"),
  make_option("--name",              type = "character", default = "signature_reversal"),
  make_option("--signature_label",   type = "character", default = "signature"),
  make_option("--response_label",    type = "character", default = "response"),
  make_option("--padj",              type = "double", default = 0.05,
              help = "adjusted p-value cutoff defining the signature [default %default]"),
  make_option("--lfc",               type = "double", default = 0.0,
              help = "absolute log2FC cutoff defining the signature [default %default]"),
  make_option("--restored_fraction", type = "double", default = 0.5,
              help = "fraction of the perturbation a gene must recover to count as restored [default %default]"),
  make_option("--label_n",           type = "integer", default = 20,
              help = "number of genes to label on the scatter [default %default]"),
  make_option("--signature_groups",  type = "character", default = "",
              help = "comma-separated conditions used by the signature contrast"),
  make_option("--response_groups",   type = "character", default = "",
              help = "comma-separated conditions used by the response contrast"),
  make_option("--n_perm",            type = "integer", default = 1000,
              help = "permutations used to build the null [default %default]")
)

args   <- parse_args(OptionParser(option_list = option_list))
outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_results <- function(path, suffix) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  needed <- c("gene_id", "log2FoldChange", "padj")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop("'", path, "' is missing column(s): ", paste(missing, collapse = ", "))
  }
  df <- df[, needed]
  colnames(df)[2:3] <- paste0(c("lfc_", "padj_"), suffix)
  df
}

sig  <- read_results(args$results_signature, "sig")
resp <- read_results(args$results_response,  "resp")

merged <- merge(sig, resp, by = "gene_id")
merged <- merged[is.finite(merged$lfc_sig) & is.finite(merged$lfc_resp), ]
message("Genes tested in both contrasts: ", nrow(merged))

# ─── Define the signature ────────────────────────────────────────────────────
in_sig <- !is.na(merged$padj_sig) & merged$padj_sig < args$padj &
          abs(merged$lfc_sig) >= args$lfc
sel <- merged[in_sig, ]
message("Signature genes (padj < ", args$padj, ", |log2FC| >= ", args$lfc, "): ",
        nrow(sel))

if (nrow(sel) < 10) {
  stop("Only ", nrow(sel), " signature gene(s) - too few to score a reversal. ",
       "Loosen --padj/--lfc or check that the contrast is the one you meant.")
}

# ─── Shared-sample check ─────────────────────────────────────────────────────
# If the two contrasts have a condition in common, its sampling noise enters one
# log2FC positively and the other negatively, which manufactures a negative
# correlation with no biology involved. On pure noise, sharing one group of three
# replicates yields ~66% "reversed" and rho ~ -0.5. Such a pair cannot be used as
# a rescue score OR as a negative control.
split_csv <- function(x) {
  if (is.null(x)) return(character(0))
  p <- trimws(strsplit(x, ",")[[1]]); p[nzchar(p)]
}
shared_groups <- intersect(split_csv(args$signature_groups),
                           split_csv(args$response_groups))
shares <- length(shared_groups) > 0
if (shares) {
  message("WARNING: the signature and response contrasts share condition(s): ",
          paste(shared_groups, collapse = ", "),
          ". Their shared sampling noise induces a spurious negative correlation, ",
          "so the reversal score is biased upward and is NOT interpretable. ",
          "Choose contrasts with no condition in common.")
}

# ─── Per-gene classification ─────────────────────────────────────────────────
sel$reversed <- sign(sel$lfc_resp) == -sign(sel$lfc_sig) & sel$lfc_resp != 0
# Recovered fraction of the perturbation: 1 = fully back to baseline, >1 = overshoot,
# negative = pushed further in the direction of the perturbation.
sel$recovered_fraction <- -sel$lfc_resp / sel$lfc_sig
sel$restored <- sel$recovered_fraction >= args$restored_fraction
sel$direction_in_signature <- ifelse(sel$lfc_sig > 0, "up", "down")
sel$response_significant <- !is.na(sel$padj_resp) & sel$padj_resp < args$padj

# ─── Statistics ──────────────────────────────────────────────────────────────
# Orthogonal (total least squares) regression: both axes are estimated
# quantities with comparable noise, so ordinary least squares would bias the
# slope toward zero and understate the reversal.
tls_slope <- function(x, y) {
  v <- prcomp(cbind(x, y), center = TRUE, scale. = FALSE)$rotation
  v[2, 1] / v[1, 1]
}

# ─── Permutation null ────────────────────────────────────────────────────────
# The binomial test assumes 50% reversal under no effect, which only holds if the
# signs in both contrasts are balanced. Shuffling the response across genes
# preserves both sign distributions and gives the baseline this statistic
# actually has on this data.
perm_null <- function(x, y, B) {
  vapply(seq_len(B), function(i) {
    yp <- sample(y)
    c(pct = 100 * mean(sign(yp) == -sign(x) & yp != 0),
      score = -tls_slope(x, yp))
  }, numeric(2))
}
set.seed(1)
null_mat  <- perm_null(sel$lfc_sig, sel$lfc_resp, args$n_perm)
null_pct  <- null_mat["pct", ]
null_score <- null_mat["score", ]

slope     <- tls_slope(sel$lfc_sig, sel$lfc_resp)
ols       <- coef(lm(lfc_resp ~ lfc_sig, data = sel))
rho       <- suppressWarnings(cor(sel$lfc_sig, sel$lfc_resp, method = "spearman"))
pear      <- suppressWarnings(cor(sel$lfc_sig, sel$lfc_resp, method = "pearson"))
n_rev     <- sum(sel$reversed)
bt        <- binom.test(n_rev, nrow(sel), p = 0.5)

stats <- data.frame(
  comparison            = args$name,
  signature_contrast    = args$signature_label,
  response_contrast     = args$response_label,
  n_genes_shared        = nrow(merged),
  n_signature_genes     = nrow(sel),
  padj_cutoff           = args$padj,
  lfc_cutoff            = args$lfc,
  tls_slope             = round(slope, 4),
  reversal_score        = round(-slope, 4),
  ols_slope             = round(ols[2], 4),
  spearman_rho          = round(rho, 4),
  pearson_r             = round(pear, 4),
  r_squared             = round(pear^2, 4),
  n_reversed            = n_rev,
  pct_reversed          = round(100 * n_rev / nrow(sel), 2),
  binom_p               = signif(bt$p.value, 3),
  null_pct_reversed     = round(mean(null_pct), 2),
  excess_over_null      = round(100 * n_rev / nrow(sel) - mean(null_pct), 2),
  perm_p                = signif((sum(null_pct >= 100 * n_rev / nrow(sel)) + 1) /
                                 (args$n_perm + 1), 3),
  null_reversal_score   = round(mean(null_score), 4),
  shares_condition      = if (shares) paste(shared_groups, collapse = ";") else "",
  interpretable         = !shares,
  restored_fraction_cut = args$restored_fraction,
  n_restored            = sum(sel$restored, na.rm = TRUE),
  pct_restored          = round(100 * sum(sel$restored, na.rm = TRUE) / nrow(sel), 2),
  median_recovered_frac = round(median(sel$recovered_fraction, na.rm = TRUE), 4),
  n_response_significant = sum(sel$response_significant),
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(outdir, "reversal_stats.csv"), row.names = FALSE)
print(t(stats))

write.csv(sel[order(-abs(sel$lfc_sig)), ], file.path(outdir, "reversal_genes.csv"),
          row.names = FALSE)

# ─── Scatter ─────────────────────────────────────────────────────────────────
lab_idx <- order(-abs(sel$lfc_sig))[seq_len(min(args$label_n, nrow(sel)))]
sel$label <- ""
sel$label[lab_idx] <- sel$gene_id[lab_idx]

lim <- max(abs(c(sel$lfc_sig, sel$lfc_resp)), na.rm = TRUE) * 1.05

excess  <- stats$excess_over_null
perm_p  <- stats$perm_p
verdict <- if (shares) {
  "NOT INTERPRETABLE - contrasts share a condition"
} else if (perm_p >= 0.05) {
  "no reversal above the null"
} else if (excess >= 15) {
  "strong reversal"
} else if (excess >= 5) {
  "partial reversal"
} else {
  "weak reversal"
}

# Stats are drawn inside the panel rather than only in the subtitle: a long
# subtitle is silently clipped when the figure is embedded or scaled, and these
# numbers are the point of the plot.
stat_box <- paste(
  sprintf("slope     %+.2f", slope),
  sprintf("R2        %.3f", pear^2),
  sprintf("rho       %+.2f", rho),
  sprintf("reversed  %.0f%% (null %.0f%%)", stats$pct_reversed, stats$null_pct_reversed),
  sprintf("perm p    %s", format(signif(perm_p, 2))),
  sep = "\n")

p <- ggplot(sel, aes(x = lfc_sig, y = lfc_resp)) +
  annotate("rect", xmin = 0, xmax = lim,  ymin = -lim, ymax = 0,
           fill = "#4DAF4A", alpha = 0.07) +
  annotate("rect", xmin = -lim, xmax = 0, ymin = 0,    ymax = lim,
           fill = "#4DAF4A", alpha = 0.07) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = response_significant), size = 1.4, alpha = 0.65) +
  geom_abline(slope = slope, intercept = 0, colour = "#E41A1C", linewidth = 0.8) +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 30,
                  segment.size = 0.2, show.legend = FALSE) +
  annotate("label", x = lim, y = lim, hjust = 1, vjust = 1, label = stat_box,
           size = 2.9, family = "mono", lineheight = 1.15,
           fill = grDevices::adjustcolor("white", alpha.f = 0.88),
           label.size = 0.25, colour = "grey15") +
  scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "#377EB8"),
                      name = paste0("padj < ", args$padj, "\nin response")) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title = paste0("Signature reversal: ", args$name),
    subtitle = sprintf("%d signature genes (padj < %s, |log2FC| >= %s)\n%s",
                       nrow(sel), args$padj, args$lfc, verdict),
    x = paste0("log2FC  ", args$signature_label, "   (the perturbation)"),
    y = paste0("log2FC  ", args$response_label, "   (the intervention)"),
    caption = if (shares) paste0(
      "NOT INTERPRETABLE: these contrasts share condition(s) ",
      paste(shared_groups, collapse = ", "), ",\nwhose shared noise creates ",
      "negative correlation on its own.") else paste0(
      "shaded quadrants = reversal | dashed line = complete reversal (slope -1)\n",
      "red line = fitted slope (total least squares) | null from ",
      args$n_perm, " permutations")
  ) +
  theme_bw(base_size = 12) +
  theme(plot.subtitle = element_text(lineheight = 1.2),
        plot.caption  = element_text(hjust = 0, size = 8, lineheight = 1.2))

ggsave(file.path(outdir, "reversal_scatter.pdf"), plot = p, width = 9, height = 9)
message("Written: ", file.path(outdir, "reversal_scatter.pdf"))

# ─── Recovered-fraction distribution ─────────────────────────────────────────
# Clipped for display: overshooting genes have unbounded ratios.
plot_df <- sel
plot_df$clipped <- pmax(pmin(plot_df$recovered_fraction, 2), -1)

q <- ggplot(plot_df, aes(x = clipped, fill = direction_in_signature)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#4DAF4A") +
  geom_vline(xintercept = median(plot_df$recovered_fraction, na.rm = TRUE),
             colour = "#E41A1C") +
  labs(title = paste0("Recovery per gene: ", args$name),
       subtitle = paste0(
         "0 = no response | 1 = fully back to baseline (green dashed) | ",
         "red = median | values clipped to [-1, 2]"),
       x = "fraction of the perturbation recovered",
       y = "genes", fill = "direction in\nthe signature") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "recovery_distribution.pdf"), plot = q, width = 8, height = 5)
message("Written: ", file.path(outdir, "recovery_distribution.pdf"))

message("Signature reversal complete: ", args$name, " - ", verdict)
