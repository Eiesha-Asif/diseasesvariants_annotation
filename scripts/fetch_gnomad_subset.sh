#!/usr/bin/env bash
###############################################################################
# fetch_gnomad_subset.sh
#
# For every variant in an input VCF, queries the remote, tabix-indexed
# gnomAD v4.1 genomes VCF directly (byte-range HTTP requests) and extracts
# ONLY the matching CHROM/POS/REF/ALT records. This avoids downloading the
# full per-chromosome gnomAD files (tens of GB each).
#
# Usage:
#   bash fetch_gnomad_subset.sh <input.vcf.gz> <output_prefix>
#
# Produces:
#   <output_prefix>.bed.gz + .tbi   (columns: chrom, start0, end, AC, AN, AF)
###############################################################################
set -Eeuo pipefail

INPUT_VCF="${1:?Usage: fetch_gnomad_subset.sh <input.vcf.gz> <output_prefix>}"
OUT_PREFIX="${2:?Usage: fetch_gnomad_subset.sh <input.vcf.gz> <output_prefix>}"
GNOMAD_BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

log() { printf '\n\033[1;34m[fetch_gnomad_subset]\033[0m %s\n' "$*"; }

RAW_BED="${WORKDIR}/gnomad_raw.bed"
: > "$RAW_BED"

log "Reading variants from $INPUT_VCF ..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$INPUT_VCF" | while IFS=$'\t' read -r CHROM POS REF ALT; do
  CHR_NUM="${CHROM#chr}"
  URL="${GNOMAD_BASE_URL}/gnomad.genomes.v4.1.sites.chr${CHR_NUM}.vcf.bgz"
  REGION="${CHROM}:${POS}-${POS}"
  log "Querying gnomAD for ${REGION} ${REF}>${ALT} ..."

  MATCH=$(tabix "$URL" "$REGION" 2>/dev/null | awk -F'\t' -v ref="$REF" -v alt="$ALT" '$4==ref && $5==alt {print; exit}') || true

  if [[ -n "$MATCH" ]]; then
    AC=$(echo "$MATCH" | grep -oP '(?<=;AC=)[^;]*' | head -1)
    AN=$(echo "$MATCH" | grep -oP '(?<=;AN=)[^;]*' | head -1)
    AF=$(echo "$MATCH" | grep -oP '(?<=;AF=)[^;]*' | head -1)
    START0=$((POS - 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CHROM" "$START0" "$POS" "${AC:-.}" "${AN:-.}" "${AF:-.}" >> "$RAW_BED"
    log "  -> found: AC=${AC:-.} AN=${AN:-.} AF=${AF:-.}"
  else
    log "  -> not observed in gnomAD (variant absent from population database)"
  fi
done

if [[ -s "$RAW_BED" ]]; then
  sort -k1,1 -k2,2n "$RAW_BED" > "${OUT_PREFIX}.bed"
  bgzip -f "${OUT_PREFIX}.bed"
  tabix -p bed "${OUT_PREFIX}.bed.gz"
  log "Wrote ${OUT_PREFIX}.bed.gz"
else
  log "No variants were found in gnomAD; no BED file produced. This is expected for very rare/de novo pathogenic variants."
fi
