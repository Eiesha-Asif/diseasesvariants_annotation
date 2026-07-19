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
#   3. bash scripts/build_example_databases.sh
#      (builds gnomAD / dbNSFP / wINTERVAR lookup files for the example VCF
#       from the pre-collected data in databases_source/. To annotate a
#       DIFFERENT VCF, first add rows for your variants to the .tsv files in
#       databases_source/ -- see README.md Section 5 for how that data was
#       originally retrieved -- then re-run this script.)
#   4. Copy config/annotation_resources.env.example to
#      config/annotation_resources.env (no edits needed for the example VCF).
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

# ----------------------------------------------
