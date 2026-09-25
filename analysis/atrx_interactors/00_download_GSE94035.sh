#!/usr/bin/env bash
# Download and organise GEO GSE94035 (Rifatbegovic et al., Int J Cancer 2018,
# PMID 28921546) for the ATRX-interactor expression boxplots.
#
# Only the authors' processed matrix is used (DESeq2-normalised FPM, log2,
# GRCh37 / Ensembl 75, rows = Ensembl gene IDs); raw reads are not downloaded.
#
# Layout created under OUTDIR:
#   raw/GSE94035_Deseq_mat.txt.gz        expression matrix (genes x 86 samples)
#   raw/GSE94035_series_matrix.txt.gz    GEO sample metadata
#   raw/ENA_PRJNA368627_runs.tsv         GSM -> SRR/SRX mapping (from ENA)
#   raw/md5sums.txt                      checksums of the downloaded files
#   metadata/samplesheet.tsv             one row per sample (see columns below)
#   download.log                         provenance: date, URLs, sample counts
#
# samplesheet.tsv columns:
#   matrix_column  column name in GSE94035_Deseq_mat.txt.gz (e.g. D01d)
#   gsm            GEO sample accession
#   title          GEO sample title
#   cell_type      Tumor | DTC | MNC
#   timepoint      diagnosis | relapse
#   enriched       yes | no (only D07r2, "without enrichment", is "no")
#   group          cell_type_timepoint, suffixed _unenriched when enriched=no
#   patient        patient ID from the title (p01 ...)
#   srx, srr       SRA experiment / run accessions (empty if ENA unreachable)
#
# Usage:
#   bash analysis/atrx_interactors/00_download_GSE94035.sh [OUTDIR]

set -euo pipefail

OUTDIR="${1:-/nobackup/lab_taschner-mandl/arthurdondi/data/GSE94035_Fikret}"

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94035"
MATRIX_URL="${GEO_BASE}/suppl/GSE94035_Deseq_mat.txt.gz"
SERIES_URL="${GEO_BASE}/matrix/GSE94035_series_matrix.txt.gz"
ENA_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA368627&result=read_run&fields=run_accession,experiment_accession,sample_alias,sample_title&format=tsv"

N_EXPECTED=86

mkdir -p "${OUTDIR}"/{raw,metadata}
LOG="${OUTDIR}/download.log"
exec > >(tee -a "${LOG}") 2>&1
echo "=== $(date -Iseconds)  00_download_GSE94035.sh -> ${OUTDIR}"

fetch() {  # fetch URL DEST ; skips files that already exist
  local url="$1" dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "exists, skipping: ${dest}"
    return 0
  fi
  echo "downloading: ${url}"
  if command -v wget >/dev/null 2>&1; then
    wget -q -O "${dest}.part" "${url}"
  else
    curl -fsSL -o "${dest}.part" "${url}"
  fi
  mv "${dest}.part" "${dest}"
}

# ─── Download ────────────────────────────────────────────────────────────────
fetch "${MATRIX_URL}" "${OUTDIR}/raw/GSE94035_Deseq_mat.txt.gz"
fetch "${SERIES_URL}" "${OUTDIR}/raw/GSE94035_series_matrix.txt.gz"
# ENA only supplies SRR/SRX accessions; not needed for plotting, so non-fatal.
fetch "${ENA_URL}" "${OUTDIR}/raw/ENA_PRJNA368627_runs.tsv" \
  || { echo "WARNING: ENA unreachable; srx/srr columns will be empty"; : > "${OUTDIR}/raw/ENA_PRJNA368627_runs.tsv"; }

gzip -t "${OUTDIR}/raw/GSE94035_Deseq_mat.txt.gz"
gzip -t "${OUTDIR}/raw/GSE94035_series_matrix.txt.gz"
(cd "${OUTDIR}/raw" && md5sum GSE94035_Deseq_mat.txt.gz GSE94035_series_matrix.txt.gz > md5sums.txt)

# ─── Build samplesheet ───────────────────────────────────────────────────────
# The matrix has Windows (CRLF) line endings: strip \r from the header.
# "|| true": head closes the pipe early (SIGPIPE), which pipefail would flag.
MATRIX_HEADER="$(zcat "${OUTDIR}/raw/GSE94035_Deseq_mat.txt.gz" | head -n 1 | tr -d '\r' || true)"

zcat "${OUTDIR}/raw/GSE94035_series_matrix.txt.gz" | awk -F'\t' \
  -v header="${MATRIX_HEADER}" \
  -v ena="${OUTDIR}/raw/ENA_PRJNA368627_runs.tsv" '
  function unq(s) { gsub(/^"|"$/, "", s); return s }
  BEGIN {
    OFS = "\t"
    # Matrix columns (skip the first, "RowName")
    nh = split(header, hc, "\t")
    for (i = 2; i <= nh; i++) in_matrix[hc[i]] = 1
    # ENA: sample_alias is the GSM accession
    while ((getline line < ena) > 0) {
      split(line, f, "\t")
      if (f[3] ~ /^GSM/) { srr[f[3]] = f[1]; srx[f[3]] = f[2] }
    }
  }
  $1 == "!Sample_title"         { for (i = 2; i <= NF; i++) title[i] = unq($i); n = NF }
  $1 == "!Sample_geo_accession" { for (i = 2; i <= NF; i++) gsm[i]   = unq($i) }
  # Two !Sample_description rows exist; the matrix code is the short D/M/T one
  # (the other is a copy-pasted Ewing sarcoma description).
  $1 == "!Sample_description" && unq($2) ~ /^[DMT][0-9]+/ {
    for (i = 2; i <= NF; i++) code[i] = unq($i)
  }
  END {
    print "matrix_column", "gsm", "title", "cell_type", "timepoint",
          "enriched", "group", "patient", "srx", "srr"
    for (i = 2; i <= n; i++) {
      t = title[i]
      cell = (t ~ /^DTC/) ? "DTC" : (t ~ /^MNC/) ? "MNC" : (t ~ /^Tumor/) ? "Tumor" : "NA"
      tp   = (t ~ /diagnosis/) ? "diagnosis" : (t ~ /relapse/) ? "relapse" : "NA"
      enr  = (t ~ /without enrichment/) ? "no" : "yes"
      pat  = t; sub(/.*Patient /, "", pat)
      grp  = cell "_" tp (enr == "no" ? "_unenriched" : "")
      col  = code[i]
      # Known mismatches between GEO descriptions and matrix column names:
      #   D07r (without enrichment) is "D07r2" in the matrix,
      #   D36d is "D36NA" in the matrix.
      if (enr == "no" && (col "2") in in_matrix) col = col "2"
      if (!(col in in_matrix)) {
        pre = substr(col, 1, 3)
        for (c in in_matrix) if (substr(c, 1, 3) == pre && !(c in used) && c ~ /NA$/) { col = c; break }
      }
      if (!(col in in_matrix)) col = "MISSING_" code[i]
      used[col] = 1
      print col, gsm[i], t, cell, tp, enr, grp, pat, srx[gsm[i]], srr[gsm[i]]
    }
  }' > "${OUTDIR}/metadata/samplesheet.tsv"

# ─── Sanity checks ───────────────────────────────────────────────────────────
SHEET="${OUTDIR}/metadata/samplesheet.tsv"
n_samples=$(( $(wc -l < "${SHEET}") - 1 ))
n_matrix=$(( $(awk -F'\t' '{print NF; exit}' <<< "${MATRIX_HEADER}") - 1 ))
n_genes=$(( $(zcat "${OUTDIR}/raw/GSE94035_Deseq_mat.txt.gz" | wc -l) - 1 ))

echo "samples in GEO metadata : ${n_samples}"
echo "samples in matrix       : ${n_matrix}"
echo "genes in matrix         : ${n_genes}"

status=0
if [[ "${n_samples}" -ne "${N_EXPECTED}" || "${n_matrix}" -ne "${N_EXPECTED}" ]]; then
  echo "ERROR: expected ${N_EXPECTED} samples"; status=1
fi
if grep -q "MISSING_" "${SHEET}"; then
  echo "ERROR: samples without a matrix column:"; grep "MISSING_" "${SHEET}"; status=1
fi
dups=$(cut -f1 "${SHEET}" | sort | uniq -d)
if [[ -n "${dups}" ]]; then
  echo "ERROR: matrix columns assigned twice: ${dups}"; status=1
fi

echo "samples per group:"
tail -n +2 "${SHEET}" | cut -f7 | sort | uniq -c
echo "samplesheet: ${SHEET}"
exit "${status}"
