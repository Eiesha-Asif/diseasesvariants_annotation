#!/usr/bin/env bash
###############################################################################
# rare_disease_vcf_annotation_pipeline.sh
#
# Annotates a small-variant (SNV/indel) VCF for rare/genetic-disorder
# interpretation, combining nine independent annotation sources into one
# final VCF:
#
#   1. bcftools norm        - normalize/split multi-allelic records
#   2. Ensembl VEP           - gene, transcript, consequence, HGVS
#   3. SnpEff                - independent consequence cross-check
#   4. ClinVar                - clinical significance / disease names
#   5. gnomAD                 - population allele frequency
#   6. ClinGen                - gene dosage sensitivity
#   7. SpliceAI                - splice-disruption prediction
#   8. dbNSFP (REVEL/AlphaMissense/CADD) - missense pathogenicity scores
#   9. ACMG/AMP classification (wINTERVAR) - automated variant classification
#
# BEFORE running this script:
#   1. bash scripts/setup_tools.sh        (installs all required software)
#   2. bash scripts/setup_databases.sh    (downloads reference/ClinVar/ClinGen)
#   3. bash scripts/fetch_gnomad_subset.sh  <input.vcf.gz> databases/gnomad/subset
#   4. bash scripts/fetch_dbnsfp_subset.sh  <input.vcf.gz> databases/dbnsfp/subset
#   5. bash scripts/prepare_intervar_queries.sh <input.vcf.gz>   (manual web step)
#      bash scripts/fetch_intervar_subset.sh <csv_folder> databases/intervar/subset
#   6. Copy config/annotation_resources.env.example to
#      config/annotation_resources.env and fill in any custom paths.
#
# Usage:
#   bash rare_disease_vcf_annotation_pipeline.sh \
#     -i data/input.vcf.gz \
#     -o results \
#     -c config/annotation_resources.env \
#     -s SAMPLE_ID \
#     -a GRCh38 \
#     -t 4
###############################################################################
set -Eeuo pipefail
IFS=$'\n\t'

usage() {
  cat <<'USAGE'
Usage:
  bash rare_disease_vcf_annotation_pipeline.sh -i input.vcf.gz -o results_dir -c config.env [-s SAMPLE_ID] [-a GRCh38|GRCh37] [-t THREADS]

Required:
  -i  Input VCF (.vcf, .vcf.gz, or .bcf)
  -o  Output directory
  -c  Resource configuration file (see config/annotation_resources.env.example)

Optional:
  -s  Sample name/prefix. Default: derived from input filename
  -a  Assembly: GRCh38 or GRCh37. Default: GRCh38
  -t  Threads. Default: 4
USAGE
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command in PATH: $1"; }
need_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

INPUT_VCF=""; OUTDIR=""; CONFIG=""; SAMPLE=""; ASSEMBLY="GRCh38"; THREADS=4

while getopts ":i:o:c:s:a:t:h" opt; do
  case "$opt" in
    i) INPUT_VCF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    c) CONFIG="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    a) ASSEMBLY="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) die "Option -$OPTARG requires an argument" ;;
    \?) die "Unknown option: -$OPTARG" ;;
  esac
done

[[ -n "$INPUT_VCF" ]] || { usage; die "Missing -i input VCF"; }
[[ -n "$OUTDIR" ]] || { usage; die "Missing -o output directory"; }
[[ -n "$CONFIG" ]] || { usage; die "Missing -c config file"; }
need_file "$INPUT_VCF"
need_file "$CONFIG"
# shellcheck disable=SC1090
source "$CONFIG"

case "$ASSEMBLY" in
  GRCh38|grch38|hg38) ASSEMBLY="GRCh38"; SPLICEAI_ASSEMBLY="grch38" ;;
  GRCh37|grch37|hg19) ASSEMBLY="GRCh37"; SPLICEAI_ASSEMBLY="grch37" ;;
  *) die "Unsupported assembly: $ASSEMBLY. Use GRCh38 or GRCh37." ;;
esac

if [[ -z "$SAMPLE" ]]; then
  base=$(basename "$INPUT_VCF")
  SAMPLE="${base%.vcf.gz}"; SAMPLE="${SAMPLE%.vcf}"; SAMPLE="${SAMPLE%.bcf}"
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/work" "$OUTDIR/snv" "$OUTDIR/reports"
LOGFILE="$OUTDIR/logs/${SAMPLE}.pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "=============================================================="
log " Rare Disease VCF Annotation Pipeline"
log "=============================================================="
log "Sample:            $SAMPLE"
log "Assembly:          $ASSEMBLY"
log "Input VCF:         $INPUT_VCF"
log "Output directory:  $OUTDIR"

need_cmd bcftools; need_cmd bgzip; need_cmd tabix; need_cmd awk; need_cmd sed; need_cmd sort; need_cmd java
need_file "${REF_FASTA:?Set REF_FASTA in config}"
[[ -s "${REF_FASTA}.fai" ]] || die "Missing FASTA index: ${REF_FASTA}.fai (run scripts/setup_databases.sh)"

# ---------------------------------------------------------------------------
# Step 1: Standardize + normalize input VCF
# ---------------------------------------------------------------------------
RAW_VCFGZ="$OUTDIR/work/${SAMPLE}.input.vcf.gz"
NORM_VCFGZ="$OUTDIR/work/${SAMPLE}.normalized.vcf.gz"

log "Step 1/9: Normalizing input VCF (bcftools norm)"
case "$INPUT_VCF" in
  *.vcf.gz) cp -f "$INPUT_VCF" "$RAW_VCFGZ" ;;
  *.bcf)    bcftools view -Oz -o "$RAW_VCFGZ" "$INPUT_VCF" ;;
  *.vcf)    bgzip -c "$INPUT_VCF" > "$RAW_VCFGZ" ;;
  *) die "Input must be .vcf, .vcf.gz, or .bcf" ;;
esac
tabix -f -p vcf "$RAW_VCFGZ"
bcftools norm -f "$REF_FASTA" -m -any -Oz -o "$NORM_VCFGZ" "$RAW_VCFGZ"
tabix -f -p vcf "$NORM_VCFGZ"
CURRENT_VCF="$NORM_VCFGZ"
log "  -> $NORM_VCFGZ"

# ---------------------------------------------------------------------------
# Step 2: Ensembl VEP
# ---------------------------------------------------------------------------
if [[ "${RUN_VEP:-1}" == "1" ]]; then
  need_cmd vep
  need_file "${VEP_CACHE_DIR:?Set VEP_CACHE_DIR in config}/homo_sapiens" 2>/dev/null || \
    [[ -d "${VEP_CACHE_DIR}/homo_sapiens" ]] || die "VEP cache not found at $VEP_CACHE_DIR (run scripts/setup_tools.sh)"
  log "Step 2/9: Running Ensembl VEP"
  VEP_VCFGZ="$OUTDIR/snv/${SAMPLE}.vep.vcf.gz"
  zcat "$CURRENT_VCF" > "$OUTDIR/work/${SAMPLE}.for_vep.vcf"
  vep --input_file "$OUTDIR/work/${SAMPLE}.for_vep.vcf" --output_file "$VEP_VCFGZ" \
    --format vcf --vcf --compress_output bgzip --force_overwrite \
    --species homo_sapiens --assembly "$ASSEMBLY" --cache --offline \
    --dir_cache "$VEP_CACHE_DIR" --fasta "$REF_FASTA" --fork "$THREADS" \
    --symbol --canonical --mane --hgvs --numbers --protein --biotype --pick
  tabix -f -p vcf "$VEP_VCFGZ"
  CURRENT_VCF="$VEP_VCFGZ"
  log "  -> $VEP_VCFGZ"
else
  log "Step 2/9: Skipped (RUN_VEP=0)"
fi

# ---------------------------------------------------------------------------
# Step 3: SnpEff
# ---------------------------------------------------------------------------
if [[ "${RUN_SNPEFF:-1}" == "1" ]]; then
  need_file "${SNPEFF_JAR:?Set SNPEFF_JAR in config}"
  [[ -n "${SNPEFF_GENOME:-}" ]] || die "Set SNPEFF_GENOME in config"
  log "Step 3/9: Running SnpEff"
  SNPEFF_VCFGZ="$OUTDIR/snv/${SAMPLE}.snpeff.vcf.gz"
  zcat "$CURRENT_VCF" > "$OUTDIR/work/${SAMPLE}.for_snpeff.vcf"
  java -Xmx"${JAVA_MEM:-4g}" -jar "$SNPEFF_JAR" ann -v -canon -hgvs "$SNPEFF_GENOME" \
    "$OUTDIR/work/${SAMPLE}.for_snpeff.vcf" | bgzip -c > "$SNPEFF_VCFGZ"
  tabix -f -p vcf "$SNPEFF_VCFGZ"
  CURRENT_VCF="$SNPEFF_VCFGZ"
  log "  -> $SNPEFF_VCFGZ"
else
  log "Step 3/9: Skipped (RUN_SNPEFF=0)"
fi

# ---------------------------------------------------------------------------
# Step 4: ClinVar
# ---------------------------------------------------------------------------
if [[ "${RUN_CLINVAR:-1}" == "1" && -n "${CLINVAR_VCF_GZ:-}" ]]; then
  need_file "$CLINVAR_VCF_GZ"
  log "Step 4/9: Annotating ClinVar clinical significance"
  CLINVAR_VCFGZ="$OUTDIR/snv/${SAMPLE}.clinvar.vcf.gz"
  bcftools annotate -a "$CLINVAR_VCF_GZ" \
    -c "${CLINVAR_INFO_FIELDS:-INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/CLNHGVS,INFO/CLNVC,INFO/CLNVCSO,INFO/GENEINFO}" \
    -Oz -o "$CLINVAR_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$CLINVAR_VCFGZ"
  CURRENT_VCF="$CLINVAR_VCFGZ"
  log "  -> $CLINVAR_VCFGZ"
else
  log "Step 4/9: Skipped (RUN_CLINVAR=0 or CLINVAR_VCF_GZ not set)"
fi

# ---------------------------------------------------------------------------
# Step 5: gnomAD (per-variant subset produced by fetch_gnomad_subset.sh)
# ---------------------------------------------------------------------------
if [[ "${RUN_GNOMAD:-1}" == "1" && -n "${GNOMAD_BED_GZ:-}" && -s "${GNOMAD_BED_GZ}" ]]; then
  log "Step 5/9: Annotating gnomAD allele frequencies"
  GNOMAD_VCFGZ="$OUTDIR/snv/${SAMPLE}.gnomad.vcf.gz"
  GNOMAD_HEADER="$OUTDIR/work/gnomad.header.txt"
  cat > "$GNOMAD_HEADER" <<'HDR'
##INFO=<ID=GNOMAD_AC,Number=1,Type=String,Description="gnomAD v4.1 allele count">
##INFO=<ID=GNOMAD_AN,Number=1,Type=String,Description="gnomAD v4.1 allele number">
##INFO=<ID=GNOMAD_AF,Number=1,Type=String,Description="gnomAD v4.1 allele frequency">
HDR
  bcftools annotate -a "$GNOMAD_BED_GZ" -h "$GNOMAD_HEADER" \
    -c CHROM,FROM,TO,INFO/GNOMAD_AC,INFO/GNOMAD_AN,INFO/GNOMAD_AF \
    -Oz -o "$GNOMAD_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$GNOMAD_VCFGZ"
  CURRENT_VCF="$GNOMAD_VCFGZ"
  log "  -> $GNOMAD_VCFGZ"
else
  log "Step 5/9: Skipped (no GNOMAD_BED_GZ; run scripts/fetch_gnomad_subset.sh first)"
fi

# ---------------------------------------------------------------------------
# Step 6: ClinGen dosage sensitivity
# ---------------------------------------------------------------------------
if [[ "${RUN_CLINGEN:-1}" == "1" && -n "${CLINGEN_DOSAGE_BED_GZ:-}" ]]; then
  need_file "$CLINGEN_DOSAGE_BED_GZ"
  log "Step 6/9: Annotating ClinGen gene dosage sensitivity"
  CLINGEN_VCFGZ="$OUTDIR/snv/${SAMPLE}.clingen.vcf.gz"
  CLINGEN_HEADER="$OUTDIR/work/clingen.header.txt"
  cat > "$CLINGEN_HEADER" <<'HDR'
##INFO=<ID=CLINGEN_REGION,Number=.,Type=String,Description="ClinGen curated gene/region name">
##INFO=<ID=CLINGEN_HAPLO,Number=.,Type=String,Description="ClinGen haploinsufficiency score">
##INFO=<ID=CLINGEN_TRIPLO,Number=.,Type=String,Description="ClinGen triplosensitivity score">
HDR
  bcftools annotate -a "$CLINGEN_DOSAGE_BED_GZ" -h "$CLINGEN_HEADER" \
    -c CHROM,FROM,TO,INFO/CLINGEN_REGION,INFO/CLINGEN_HAPLO,INFO/CLINGEN_TRIPLO \
    -Oz -o "$CLINGEN_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$CLINGEN_VCFGZ"
  CURRENT_VCF="$CLINGEN_VCFGZ"
  log "  -> $CLINGEN_VCFGZ"
else
  log "Step 6/9: Skipped (RUN_CLINGEN=0 or CLINGEN_DOSAGE_BED_GZ not set)"
fi

# ---------------------------------------------------------------------------
# Step 7: SpliceAI
# ---------------------------------------------------------------------------
if [[ "${RUN_SPLICEAI:-1}" == "1" ]]; then
  need_cmd spliceai
  log "Step 7/9: Running SpliceAI"
  SPLICEAI_VCFGZ="$OUTDIR/snv/${SAMPLE}.spliceai.vcf.gz"
  zcat "$CURRENT_VCF" > "$OUTDIR/work/${SAMPLE}.for_spliceai.vcf"
  spliceai -I "$OUTDIR/work/${SAMPLE}.for_spliceai.vcf" -O "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" \
    -R "$REF_FASTA" -A "$SPLICEAI_ASSEMBLY"
  bgzip -c "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" > "$SPLICEAI_VCFGZ"
  tabix -f -p vcf "$SPLICEAI_VCFGZ"
  CURRENT_VCF="$SPLICEAI_VCFGZ"
  log "  -> $SPLICEAI_VCFGZ"
else
  log "Step 7/9: Skipped (RUN_SPLICEAI=0; requires 'conda activate spliceai_env' first)"
fi

# ---------------------------------------------------------------------------
# Step 8: REVEL / AlphaMissense / CADD (dbNSFP subset via MyVariant.info)
# ---------------------------------------------------------------------------
if [[ "${RUN_DBNSFP:-1}" == "1" && -n "${DBNSFP_BED_GZ:-}" && -s "${DBNSFP_BED_GZ}" ]]; then
  log "Step 8/9: Annotating REVEL / AlphaMissense / CADD scores"
  DBNSFP_VCFGZ="$OUTDIR/snv/${SAMPLE}.dbnsfp.vcf.gz"
  DBNSFP_HEADER="$OUTDIR/work/dbnsfp.header.txt"
  cat > "$DBNSFP_HEADER" <<'HDR'
##INFO=<ID=REVEL_SCORE,Number=1,Type=String,Description="REVEL missense pathogenicity score (dbNSFP)">
##INFO=<ID=ALPHAMISSENSE_SCORE,Number=1,Type=String,Description="AlphaMissense max score across transcripts (dbNSFP)">
##INFO=<ID=ALPHAMISSENSE_PRED,Number=1,Type=String,Description="AlphaMissense prediction: P=Pathogenic, A=Ambiguous, B=Benign">
##INFO=<ID=CADD_PHRED_DBNSFP,Number=1,Type=String,Description="CADD phred score (dbNSFP)">
HDR
  bcftools annotate -a "$DBNSFP_BED_GZ" -h "$DBNSFP_HEADER" \
    -c CHROM,FROM,TO,INFO/REVEL_SCORE,INFO/ALPHAMISSENSE_SCORE,INFO/ALPHAMISSENSE_PRED,INFO/CADD_PHRED_DBNSFP \
    -Oz -o "$DBNSFP_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$DBNSFP_VCFGZ"
  CURRENT_VCF="$DBNSFP_VCFGZ"
  log "  -> $DBNSFP_VCFGZ"
else
  log "Step 8/9: Skipped (no DBNSFP_BED_GZ; run scripts/fetch_dbnsfp_subset.sh first)"
fi

# ---------------------------------------------------------------------------
# Step 9: ACMG/AMP classification (wINTERVAR subset)
# ---------------------------------------------------------------------------
if [[ "${RUN_INTERVAR:-1}" == "1" && -n "${INTERVAR_BED_GZ:-}" && -s "${INTERVAR_BED_GZ}" ]]; then
  log "Step 9/9: Annotating ACMG/AMP classification"
  INTERVAR_VCFGZ="$OUTDIR/snv/${SAMPLE}.intervar.vcf.gz"
  INTERVAR_HEADER="$OUTDIR/work/intervar.header.txt"
  cat > "$INTERVAR_HEADER" <<'HDR'
##INFO=<ID=INTERVAR_GENE,Number=1,Type=String,Description="Gene symbol reported by wINTERVAR">
##INFO=<ID=INTERVAR_ACMG,Number=1,Type=String,Description="ACMG/AMP classification reported by wINTERVAR">
HDR
  bcftools annotate -a "$INTERVAR_BED_GZ" -h "$INTERVAR_HEADER" \
    -c CHROM,FROM,TO,INFO/INTERVAR_GENE,INFO/INTERVAR_ACMG \
    -Oz -o "$INTERVAR_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$INTERVAR_VCFGZ"
  CURRENT_VCF="$INTERVAR_VCFGZ"
  log "  -> $INTERVAR_VCFGZ"
else
  log "Step 9/9: Skipped (no INTERVAR_BED_GZ; run scripts/fetch_intervar_subset.sh first)"
fi

# ---------------------------------------------------------------------------
# Final output
# ---------------------------------------------------------------------------
FINAL_VCFGZ="$OUTDIR/${SAMPLE}.final.annotated.vcf.gz"
cp -f "$CURRENT_VCF" "$FINAL_VCFGZ"
tabix -f -p vcf "$FINAL_VCFGZ"

REPORT="$OUTDIR/reports/${SAMPLE}.annotation_summary.txt"
{
  echo "Sample:      $SAMPLE"
  echo "Assembly:    $ASSEMBLY"
  echo "Final VCF:   $FINAL_VCFGZ"
  echo "Pipeline log: $LOGFILE"
  echo
  echo "Per-variant summary:"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/GNOMAD_AF\t%INFO/REVEL_SCORE\t%INFO/INTERVAR_ACMG\n' "$FINAL_VCFGZ" 2>/dev/null || true
} > "$REPORT"

log "=============================================================="
log "Pipeline finished successfully."
log "Final annotated VCF: $FINAL_VCFGZ"
log "Summary report:      $REPORT"
log "=============================================================="
