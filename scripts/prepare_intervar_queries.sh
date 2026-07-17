#!/usr/bin/env bash
###############################################################################
# prepare_intervar_queries.sh
#
# Prints the Chr / Position / Ref / Alt values for every variant in the
# input VCF, ready to paste into wINTERVAR's "Query by genomic coordinate"
# form at https://wintervar.wglab.org (GRCh38 build). wINTERVAR has no
# public API, so this step must be done manually in a browser, once per
# variant; the CSV result of each query should be saved and passed to
# fetch_intervar_subset.sh afterwards.
#
# Usage:
#   bash prepare_intervar_queries.sh <input.vcf.gz>
###############################################################################
set -Eeuo pipefail
INPUT_VCF="${1:?Usage: prepare_intervar_queries.sh <input.vcf.gz>}"

echo "Go to https://wintervar.wglab.org , select GRCh38, and submit each row below"
echo "using 'Query by genomic coordinate'. Save each result as a CSV file."
echo
printf '%-6s %-14s %-6s %-6s\n' "Chr" "Position" "Ref" "Alt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$INPUT_VCF" | sed 's/^chr//' | \
  awk -F'\t' '{printf "%-6s %-14s %-6s %-6s\n", $1, $2, $3, $4}'
