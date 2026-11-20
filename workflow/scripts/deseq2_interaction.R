#!/usr/bin/env Rscript
# deseq2_interaction.R — difference of differences between two contrasts.
#
# Estimates  (A1 - A2) - (B1 - B2)  on the log2 scale, i.e. "does the effect of
# one comparison differ from the effect of another".
#
# Why this exists.  When two groups differ by more than the variable of
# interest (a different clone, a transduction step, a selection round), the
# nuisance effect cannot be estimated as a covariate: it is a deterministic
# function of the group label, so adding it to the design makes the model
# matrix rank-deficient and DESeq2 refuses to fit it.  There is no extra degree
# of freedom to spend on it.
#
# The difference of differences sidesteps that entirely.  It is a linear
# combination of the group means, so it is estimable from the plain
# ~ condition fit, and any effect shared by both pairs cancels in the
# subtraction rather than being modelled.  For a rescue experiment:
#
#   group_A: [ATRX_FL, EmptyVector]   the rescue, inside the transduced background
#   group_B: [E6, TP53]               the knockout, inside the untransduced one
#
#   log2FC > 0 : the gene moves further up in the rescue than in the knockout
#   log2FC ~ 0 : the two comparisons move the gene by the same amount
#
# The assumption being made is additivity: that the nuisance effect is the same
# size in both pairs.  With one group per condition that assumption cannot be
# tested from the data — it needs a group that breaks the confounding (see
# docs/design_confounding_plan.md).  It is a much weaker and more transparent
# assumption than the one hidden inside a batch-correction step, but it is an
# assumption, and it belongs in the methods section.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option("--counts",  type = "character", help = "featureCounts output file"),
  make_option("--outdir",  type = "character", help = "Output directory"),
  make_option("--samples", type = "character", help = "Comma-separated sample:condition pairs"),
  make_option("--name",    type = "character", default = "interaction"),
  make_option("--group_a", type = "character", help = "'numerator denominator' for pair A"),
  make_option("--group_b", type = "character", help = "'numerator denominator' for pair B"),
  make_option("--padj",    type = "double", default = 0.05),
  make_option("--lfc",     type = "double", default = 1.0)
)

args     <- parse_args(OptionParser(option_list = option_list))
outdir   <- args$outdir
padj_thr <- args$padj
lfc_thr  <- args$lfc
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

a <- strsplit(trimws(args$group_a), "\\s+")[[1]]
b <- strsplit(trimws(args$group_b), "\\s+")[[1]]
if (length(a) != 2 || length(b) != 2) {
  stop("--group_a and --group_b must each name exactly two conditions")
}
a1 <- a[1]; a2 <- a[2]; b1 <- b[1]; b2 <- b[2]

# ─── Sample table ────────────────────────────────────────────────────────────
sample_df <- do.call(rbind, lapply(strsplit(args$samples, ",")[[1]], function(x) {
  parts <- strsplit(x, ":")[[1]]
  data.frame(sample = parts[1], condition = parts[2], stringsAsFactors = FALSE)
}))
rownames(sample_df) <- sample_df$sample
sample_df$condition <- factor(sample_df$condition)

missing <- setdiff(c(a1, a2, b1, b2), levels(sample_df$condition))
if (length(missing) > 0) {
  stop("Condition(s) not present in the samples: ", paste(missing, collapse = ", "))
}

# ─── Counts ──────────────────────────────────────────────────────────────────
raw <- read.table(args$counts, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE)
clean_names <- sub(".*/", "", colnames(raw)[7:ncol(raw)])
clean_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
counts <- as.matrix(raw[, 7:ncol(raw)])
rownames(counts) <- raw$Geneid
colnames(counts) <- clean_names

shared    <- intersect(colnames(counts), rownames(sample_df))
counts    <- counts[, shared, drop = FALSE]
sample_df <- sample_df[shared, , drop = FALSE]

# ─── Fit ─────────────────────────────────────────────────────────────────────
# The reference level is set to one of the four conditions so that its
# coefficient is the one dropped; the difference of differences is invariant to
# this choice, but fixing it keeps resultsNames() predictable.
sample_df$condition <- relevel(droplevels(sample_df$condition), ref = b2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_df,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# ─── Build the numeric contrast ──────────────────────────────────────────────
# With design ~ condition, the mean of group g is (Intercept) + beta_g, where
# beta_ref = 0.  So (A1 - A2) - (B1 - B2) = b_A1 - b_A2 - b_B1 + b_B2 and the
# intercept cancels — which is exactly why this quantity is estimable while a
# free-standing "nuisance effect" is not.
rn  <- resultsNames(dds)
ref <- levels(sample_df$condition)[1]
vec <- setNames(rep(0, length(rn)), rn)

add_group <- function(vec, group, weight) {
  if (group == ref) return(vec)          # beta_ref is fixed at 0
  coef_name <- paste0("condition_", group, "_vs_", ref)
  if (!coef_name %in% names(vec)) {
    stop("Coefficient '", coef_name, "' not found. resultsNames(dds) = ",
         paste(names(vec), collapse = ", "))
  }
  vec[coef_name] <- vec[coef_name] + weight
  vec
}

vec <- add_group(vec, a1, +1)
vec <- add_group(vec, a2, -1)
vec <- add_group(vec, b1, -1)
vec <- add_group(vec, b2, +1)

message("Reference level: ", ref)
message("Contrast vector:")
print(vec[vec != 0])

res <- results(dds, contrast = as.numeric(vec), alpha = padj_thr)

label_a <- paste0("(", a1, " / ", a2, ")")
label_b <- paste0("(", b1, " / ", b2, ")")
label   <- paste0(label_a, " - ", label_b)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::arrange(padj)
write.csv(res_df, file.path(outdir, "results.csv"), row.names = FALSE, quote = FALSE)

norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  tibble::rownames_to_column("gene_id")
write.csv(norm_counts, file.path(outdir, "normalized_counts.csv"),
          row.names = FALSE, quote = FALSE)

# For interpretation it helps to have both component effects next to the
# interaction, so a large difference can be traced to which pair moved.
comp <- function(num, den) {
  r <- results(dds, contrast = c("condition", num, den), alpha = padj_thr)
  data.frame(gene_id = rownames(r), lfc = r$log2FoldChange, padj = r$padj,
             stringsAsFactors = FALSE)
}
ca <- comp(a1, a2); colnames(ca)[2:3] <- c("lfc_A", "padj_A")
cb <- comp(b1, b2); colnames(cb)[2:3] <- c("lfc_B", "padj_B")
components <- merge(merge(res_df[, c("gene_id", "log2FoldChange", "padj")], ca,
                          by = "gene_id"), cb, by = "gene_id")
colnames(components)[2:3] <- c("lfc_interaction", "padj_interaction")
write.csv(components[order(components$padj_interaction), ],
          file.path(outdir, "components.csv"), row.names = FALSE, quote = FALSE)

# ─── MA plot ─────────────────────────────────────────────────────────────────
pdf(file.path(outdir, "ma_plot.pdf"), width = 6, height = 5)
plotMA(res, alpha = padj_thr, main = paste("Interaction:", args$name),
       sub = paste0("log2FC = ", label))
dev.off()

# ─── Volcano plot ────────────────────────────────────────────────────────────
volcano_df <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    significance = dplyr::if_else(
      padj < padj_thr & abs(log2FoldChange) >= lfc_thr,
      "Significant", "Not significant"),
    label = dplyr::if_else(
      padj < padj_thr & abs(log2FoldChange) >= lfc_thr, gene_id, NA_character_)
  )

p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj),
                            colour = significance, label = label)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = c("Significant" = "#E41A1C",
                                 "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.4) +
  labs(
    title    = paste("Interaction:", args$name),
    subtitle = paste0("log2FC = ", label,
                      "   |   0 means both comparisons move the gene equally"),
    x        = paste0("log2FC  ", label),
    y        = expression(-log[10]~"adjusted p-value"),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "volcano.pdf"), plot = p, width = 7, height = 6)

n_sig <- sum(volcano_df$significance == "Significant")
message("Interaction '", args$name, "' complete: ", n_sig,
        " gene(s) with padj < ", padj_thr, " and |log2FC| >= ", lfc_thr,
        ". Results written to: ", outdir)
