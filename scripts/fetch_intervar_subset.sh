#!/usr/bin/env bash
###############################################################################
# fetch_intervar_subset.sh
#
# ACMG/AMP classification (Likely Pathogenic / Uncertain Significance / etc.)
# is obtained from the wINTERVAR web tool (https://wintervar.wglab.org).
# wINTERVAR does not provide a stable public API, so this step is
# semi-automated:
#
#   1. Run 'prepare_intervar_queries.sh <input.vcf.gz>' to print the exact
#      Chr / Position / Ref / Alt values to submit, one variant at a time,
#      on the wINTERVAR "Query by genomic coordinate" form.
#   2. For each variant, download the result as CSV and save all CSVs into
#      one folder (any filename ending in .csv).
#   3. Run this script on that folder to merge the classifications into a
#      single BED file usable by the main pipeline.
#
# Usage:
#   bash fetch_intervar_subset.sh <folder_of_wintervar_csvs> <output_prefix>
#
# Produces:
#   <output_prefix>.bed.gz + .tbi   (columns: chrom, start0, end, gene, ACMG_classification)
###############################################################################
set -Eeuo pipefail

CSV_DIR="${1:?Usage: fetch_intervar_subset.sh <folder_of_wintervar_csvs> <output_prefix>}"
OUT_PREFIX="${2:?Usage: fetch_intervar_subset.sh <folder_of_wintervar_csvs> <output_prefix>}"

log() { printf '\n\033[1;34m[fetch_intervar_subset]\033[0m %s\n' "$*"; }

RAW_BED="$(mktemp)"
: > "$RAW_BED"

shopt -s nullglob
for f in "$CSV_DIR"/*.csv; do
  log "Parsing $f ..."
  awk -F'","' 'NR==2 {
    chrom = $1; gsub(/"/, "", chrom);
    pos = $2;
    gene = $5;
    classification = $6;
    gsub(/ *\(Details&Adjust\)/, "", classification);
    gsub(/[[:space:]]+/, "_", classification);
    gsub(/"/, "", classification);
    print "chr"chrom"\t"(pos-1)"\t"pos"\t"gene"\t"classification
  }' "$f" >> "$RAW_BED"
done

if [[ ! -s "$RAW_BED" ]]; then
  log "No CSV files found in $CSV_DIR (expected files exported from wINTERVAR)."
  exit 1
fi

sort -k1,1 -k2,2n "$RAW_BED" > "${OUT_PREFIX}.bed"
bgzip -f "${OUT_PREFIX}.bed"
tabix -p bed "${OUT_PREFIX}.bed.gz"
log "Wrote ${OUT_PREFIX}.bed.gz"
