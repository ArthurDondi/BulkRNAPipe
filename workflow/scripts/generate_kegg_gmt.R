#!/usr/bin/env Rscript
# generate_kegg_gmt.R — build a GMT of the CURRENT human KEGG pathways.
#
# Why this exists: MSigDB's C2:CP:KEGG is a frozen snapshot of KEGG from around
# 2011. Every pathway KEGG has added since is absent from it -- Ferroptosis,
# Cellular senescence, the synapse pathways, the whole Chromosome category. This
# fetches the live definitions from the KEGG REST API instead, so the gene sets
# match what genome.jp lists today.
#
# Gene-set names reproduce the MSigDB convention, so a name written against the
# legacy collection still resolves here:
#   "DNA replication"          -> KEGG_DNA_REPLICATION
#   "Non-homologous end-joining" -> KEGG_NON_HOMOLOGOUS_END_JOINING
#   "p53 signaling pathway"    -> KEGG_P53_SIGNALING_PATHWAY
#
# Requires outbound HTTPS to rest.kegg.jp. KEGG is free for academic use; see
# https://www.kegg.jp/kegg/legal.html before using it otherwise.
#
# Outputs:
#   kegg_current.gmt            gene sets, symbols
#   kegg_current_pathways.tsv   kegg_id / kegg_name / gs_name / n_genes, for
#                               checking a category file against real KEGG names

suppressPackageStartupMessages({
  library(KEGGREST)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(optparse)
})

option_list <- list(
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--organism", type = "character", default = "hsa"),
  make_option("--min_genes", type = "integer", default = 1L,
              help = "Drop pathways with fewer mapped symbols than this")
)
args <- parse_args(OptionParser(option_list = option_list))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

message("Fetching KEGG pathway membership for '", args$organism, "' ...")
links <- keggLink("pathway", args$organism)   # names: hsa:1234  values: path:hsa04110
if (length(links) == 0) stop("keggLink returned nothing - is rest.kegg.jp reachable?")
message("  ", length(links), " gene-pathway links")

pnames <- keggList("pathway", args$organism)  # values: "Cell cycle - Homo sapiens (human)"
if (length(pnames) == 0) stop("keggList returned nothing - is rest.kegg.jp reachable?")
message("  ", length(pnames), " pathways")

# KEGGREST has returned these ids both with and without the "path:" prefix
# depending on version; normalise both sides before joining.
strip_prefix <- function(x) sub("^path:", "", x)
path_id  <- strip_prefix(unname(links))
entrez   <- sub(paste0("^", args$organism, ":"), "", names(links))
names(pnames) <- strip_prefix(names(pnames))

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(entrez),
                             column = "SYMBOL", keytype = "ENTREZID",
                             multiVals = "first")
n_unmapped <- sum(is.na(sym))
if (n_unmapped > 0) {
  message("NOTE: ", n_unmapped, " of ", length(sym),
          " Entrez ids have no SYMBOL in org.Hs.eg.db and are dropped.")
}

# "Cell cycle - Homo sapiens (human)" -> KEGG_CELL_CYCLE. Only the species suffix
# is stripped, so "Apoptosis - multiple species" keeps its own " - " intact.
to_gs_name <- function(x) {
  x <- sub(" - Homo sapiens \\(human\\)$", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  paste0("KEGG_", toupper(x))
}

df <- data.frame(path_id = path_id, symbol = unname(sym[entrez]),
                 stringsAsFactors = FALSE)
df <- df[!is.na(df$symbol), , drop = FALSE]
sets <- split(df$symbol, df$path_id)
sets <- lapply(sets, function(g) sort(unique(g)))
sets <- sets[vapply(sets, length, integer(1)) >= args$min_genes]

kegg_name <- unname(pnames[names(sets)])
gs_name   <- to_gs_name(kegg_name)

dup <- gs_name[duplicated(gs_name)]
if (length(dup) > 0) {
  message("WARNING: duplicate gene-set names after conversion: ",
          paste(unique(dup), collapse = ", "),
          " - they are kept and fgsea will see them as separate sets.")
}

gmt_path <- file.path(args$outdir, "kegg_current.gmt")
con <- file(gmt_path, "w")
for (i in seq_along(sets)) {
  writeLines(paste(c(gs_name[i], names(sets)[i], sets[[i]]), collapse = "\t"), con)
}
close(con)
message("Written: ", gmt_path, "  (", length(sets), " gene sets)")

tsv_path <- file.path(args$outdir, "kegg_current_pathways.tsv")
write.table(
  data.frame(kegg_id = names(sets), kegg_name = kegg_name, gs_name = gs_name,
             n_genes = vapply(sets, length, integer(1)), stringsAsFactors = FALSE),
  tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", tsv_path)
