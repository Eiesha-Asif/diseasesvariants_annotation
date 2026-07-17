#!/usr/bin/env bash
###############################################################################
# fetch_dbnsfp_subset.sh
#
# For every variant in an input VCF, queries the free MyVariant.info API
# (which serves dbNSFP annotations) to retrieve REVEL, AlphaMissense, and
# CADD scores -- without downloading the multi-gigabyte dbNSFP database.
#
# Usage:
#   bash fetch_dbnsfp_subset.sh <input.vcf.gz> <output_prefix>
#
# Produces:
#   <output_prefix>.bed.gz + .tbi
#   (columns: chrom, start0, end, REVEL_score, AlphaMissense_max_score,
#             AlphaMissense_prediction, CADD_phred)
###############################################################################
set -Eeuo pipefail

INPUT_VCF="${1:?Usage: fetch_dbnsfp_subset.sh <input.vcf.gz> <output_prefix>}"
OUT_PREFIX="${2:?Usage: fetch_dbnsfp_subset.sh <input.vcf.gz> <output_prefix>}"

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

log() { printf '\n\033[1;34m[fetch_dbnsfp_subset]\033[0m %s\n' "$*"; }
command -v jq >/dev/null 2>&1 || { echo "jq is required (sudo apt-get install jq)"; exit 1; }

RAW_BED="${WORKDIR}/dbnsfp_raw.bed"
: > "$RAW_BED"

log "Reading variants from $INPUT_VCF ..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$INPUT_VCF" | while IFS=$'\t' read -r CHROM POS REF ALT; do
  HGVS="${CHROM}:g.${POS}${REF}>${ALT}"
  ENCODED_HGVS=$(python3 -c "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$HGVS")
  URL="https://myvariant.info/v1/variant/${ENCODED_HGVS}?assembly=hg38&fields=dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.cadd"

  log "Querying MyVariant.info for ${HGVS} ..."
  JSON=$(curl -s "$URL")

  REVEL=$(echo "$JSON" | jq -r '(.dbnsfp.revel.score // empty) | if type=="array" then .[0] else . end // "."')
  AM_SCORE=$(echo "$JSON" | jq -r '(.dbnsfp.alphamissense.score // empty) | if type=="array" then max else . end // "."')
  AM_PRED=$(echo "$JSON" | jq -r '(.dbnsfp.alphamissense.pred // empty) | if type=="array" then .[0] else . end // "."')
  CADD=$(echo "$JSON" | jq -r '(.dbnsfp.cadd.phred // empty) | if type=="array" then .[0] else . end // "."')

  if [[ "$REVEL" != "." || "$AM_SCORE" != "." ]]; then
    START0=$((POS - 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$CHROM" "$START0" "$POS" "$REVEL" "$AM_SCORE" "$AM_PRED" "$CADD" >> "$RAW_BED"
    log "  -> REVEL=$REVEL AlphaMissense=$AM_SCORE ($AM_PRED) CADD=$CADD"
  else
    log "  -> no dbNSFP score available for this variant (e.g. non-missense, or not covered by dbNSFP)"
  fi
  sleep 0.3   # be polite to the free public API
done

if [[ -s "$RAW_BED" ]]; then
  sort -k1,1 -k2,2n "$RAW_BED" > "${OUT_PREFIX}.bed"
  bgzip -f "${OUT_PREFIX}.bed"
  tabix -p bed "${OUT_PREFIX}.bed.gz"
  log "Wrote ${OUT_PREFIX}.bed.gz"
else
  log "No dbNSFP scores were retrieved; no BED file produced."
fi
