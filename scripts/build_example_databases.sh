#!/usr/bin/env bash
###############################################################################
# build_example_databases.sh
#
# Converts the plain-text variant data in databases_source/ (already
# retrieved from gnomAD, MyVariant.info, and wINTERVAR for the 4 example
# variants) into bgzip-compressed, tabix-indexed BED files that the main
# pipeline script can use directly -- no internet access needed to run the
# pipeline on the example VCF.
#
# If you annotate a DIFFERENT VCF, you must repeat the manual data-retrieval
# steps described in README.md Section 5 to add rows to the TSV files in
# databases_source/ (or create new ones) before re-running this script.
#
# Usage:
#   bash scripts/build_example_databases.sh
###############################################################################
set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${PROJECT_ROOT}/databases_source"
OUT="${PROJECT_ROOT}/databases"
mkdir -p "$OUT/gnomad" "$OUT/dbnsfp" "$OUT/intervar"

log() { printf '\n\033[1;34m[build_example_databases]\033[0m %s\n' "$*"; }

build() {
  local src="$1" out_prefix="$2" label="$3"
  log "Building $label -> ${out_prefix}.bed.gz"
  grep -v '^#' "$src" | sort -k1,1 -k2,2n > "${out_prefix}.bed"
  bgzip -f "${out_prefix}.bed"
  tabix -p bed "${out_prefix}.bed.gz"
}

build "$SRC/gnomad_subset.tsv"   "$OUT/gnomad/subset"   "gnomAD"
build "$SRC/dbnsfp_subset.tsv"   "$OUT/dbnsfp/subset"   "dbNSFP (REVEL/AlphaMissense/CADD)"
build "$SRC/intervar_subset.tsv" "$OUT/intervar/subset" "wINTERVAR ACMG classification"

log "Done. Indexed files are in ${OUT}/{gnomad,dbnsfp,intervar}/subset.bed.gz"
