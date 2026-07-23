#!/usr/bin/env bash
# rare_disease_vcf_annotation_pipeline.sh
# Purpose: annotate a human WGS VCF for rare/genetic-disorder interpretation.
set -Eeuo pipefail
IFS=$'\n\t'

usage() {
  cat <<'USAGE'
Usage:
  bash rare_disease_vcf_annotation_pipeline.sh -i sample.small_variants.vcf.gz -o results_dir -c annotation_resources.env [-n sample.cnv.vcf.gz|sample.cnv.bed] [-s SAMPLE_ID] [-a GRCh38|GRCh37] [-t THREADS]
USAGE
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command in PATH: $1"; }
need_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }
maybe_file() { [[ -n "${1:-}" && -s "$1" ]]; }

INPUT_VCF=""; CNV_INPUT=""; OUTDIR=""; CONFIG=""; SAMPLE=""; ASSEMBLY="GRCh38"; THREADS=8

while getopts ":i:n:o:c:s:a:t:h" opt; do
  case "$opt" in
    i) INPUT_VCF="$OPTARG" ;;
    n) CNV_INPUT="$OPTARG" ;;
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
source "$CONFIG"

case "$ASSEMBLY" in
  GRCh38|grch38|hg38) ASSEMBLY="GRCh38"; ANNOVAR_BUILDVER="${ANNOVAR_BUILDVER:-hg38}"; CLASSIFY_BUILDVER="${CLASSIFY_BUILDVER:-hg38}"; SPLICEAI_ASSEMBLY="${SPLICEAI_ASSEMBLY:-grch38}"; ANNOTSV_BUILD="${ANNOTSV_BUILD:-GRCh38}" ;;
  GRCh37|grch37|hg19) ASSEMBLY="GRCh37"; ANNOVAR_BUILDVER="${ANNOVAR_BUILDVER:-hg19}"; CLASSIFY_BUILDVER="${CLASSIFY_BUILDVER:-hg19}"; SPLICEAI_ASSEMBLY="${SPLICEAI_ASSEMBLY:-grch37}"; ANNOTSV_BUILD="${ANNOTSV_BUILD:-GRCh37}" ;;
  *) die "Unsupported assembly: $ASSEMBLY. Use GRCh38 or GRCh37." ;;
esac

if [[ -z "$SAMPLE" ]]; then
  base=$(basename "$INPUT_VCF")
  SAMPLE="${base%.vcf.gz}"; SAMPLE="${SAMPLE%.vcf}"; SAMPLE="${SAMPLE%.bcf}"
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/work" "$OUTDIR/snv" "$OUTDIR/cnv" "$OUTDIR/acmg" "$OUTDIR/reports"
LOGFILE="$OUTDIR/logs/${SAMPLE}.pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "Starting rare disease VCF annotation pipeline"
log "Sample: $SAMPLE"; log "Assembly: $ASSEMBLY"; log "Input VCF: $INPUT_VCF"; log "Output directory: $OUTDIR"

need_cmd bcftools; need_cmd bgzip; need_cmd tabix; need_cmd awk; need_cmd sed; need_cmd sort; need_cmd java; need_cmd python3
need_file "${REF_FASTA:?Set REF_FASTA in config}"
[[ -s "${REF_FASTA}.fai" ]] || die "Missing FASTA index: ${REF_FASTA}.fai. Run: samtools faidx $REF_FASTA"

RAW_VCFGZ="$OUTDIR/work/${SAMPLE}.input.vcf.gz"
NORM_VCFGZ="$OUTDIR/work/${SAMPLE}.normalized.split.vcf.gz"

log "Step 1: bgzip/index input and normalize/split multiallelic SNV/indel records"
case "$INPUT_VCF" in
  *.vcf.gz) cp -f "$INPUT_VCF" "$RAW_VCFGZ" ;;
  *.bcf) bcftools view -Oz -o "$RAW_VCFGZ" "$INPUT_VCF" ;;
  *.vcf) bgzip -c "$INPUT_VCF" > "$RAW_VCFGZ" ;;
  *) die "Input must be .vcf, .vcf.gz, or .bcf" ;;
esac
tabix -f -p vcf "$RAW_VCFGZ"
bcftools norm -f "$REF_FASTA" -m -any -Oz -o "$NORM_VCFGZ" "$RAW_VCFGZ"
tabix -f -p vcf "$NORM_VCFGZ"
log "Normalized VCF: $NORM_VCFGZ"
CURRENT_VCF="$NORM_VCFGZ"

VEP_VCFGZ="$OUTDIR/snv/${SAMPLE}.vep.vcf.gz"
if [[ "${RUN_VEP:-1}" == "1" ]]; then
  need_cmd vep
  [[ -n "${VEP_CACHE_DIR:-}" ]] || die "Set VEP_CACHE_DIR in config"
  [[ -d "$VEP_CACHE_DIR" ]] || die "VEP_CACHE_DIR does not exist: $VEP_CACHE_DIR"
  log "Step 2: run Ensembl VEP with optional plugins/scores"
  VEP_PLUGIN_ARGS=()
  [[ -n "${VEP_PLUGIN_DIR:-}" ]] && VEP_PLUGIN_ARGS+=(--dir_plugins "$VEP_PLUGIN_DIR")
  if maybe_file "${ALPHAMISSENSE_TSV_GZ:-}"; then VEP_PLUGIN_ARGS+=(--plugin "AlphaMissense,file=${ALPHAMISSENSE_TSV_GZ}"); fi
  if maybe_file "${REVEL_TSV_GZ:-}"; then VEP_PLUGIN_ARGS+=(--plugin "REVEL,file=${REVEL_TSV_GZ}"); fi
  if maybe_file "${CADD_SNV_TSV_GZ:-}" && maybe_file "${CADD_INDEL_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "CADD,snv=${CADD_SNV_TSV_GZ},indels=${CADD_INDEL_TSV_GZ}")
  elif maybe_file "${CADD_SNV_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "CADD,snv=${CADD_SNV_TSV_GZ}")
  fi
  vep --input_file "$CURRENT_VCF" --output_file "$VEP_VCFGZ" --format vcf --vcf --compress_output bgzip --force_overwrite \
    --species homo_sapiens --assembly "$ASSEMBLY" --cache --offline --dir_cache "$VEP_CACHE_DIR" --fasta "$REF_FASTA" \
    --fork "$THREADS" --everything --symbol --canonical --mane --hgvs --numbers --protein --biotype "${VEP_PLUGIN_ARGS[@]}"
  tabix -f -p vcf "$VEP_VCFGZ"
  CURRENT_VCF="$VEP_VCFGZ"
else
  log "Step 2: skipped VEP because RUN_VEP=0"
fi

SNPEFF_VCFGZ="$OUTDIR/snv/${SAMPLE}.vep.snpeff.vcf.gz"
if [[ "${RUN_SNPEFF:-1}" == "1" ]]; then
  need_file "${SNPEFF_JAR:?Set SNPEFF_JAR in config}"
  [[ -n "${SNPEFF_GENOME:-}" ]] || die "Set SNPEFF_GENOME in config"
  log "Step 3: run SnpEff to add ANN consequence field"
  zcat "$CURRENT_VCF" > "$OUTDIR/work/${SAMPLE}.for_snpeff.vcf"
  java -Xmx"${JAVA_MEM:-8g}" -jar "$SNPEFF_JAR" ann -v -canon -hgvs "$SNPEFF_GENOME" "$OUTDIR/work/${SAMPLE}.for_snpeff.vcf" | bgzip -c > "$SNPEFF_VCFGZ"
  tabix -f -p vcf "$SNPEFF_VCFGZ"
  CURRENT_VCF="$SNPEFF_VCFGZ"
else
  log "Step 3: skipped SnpEff because RUN_SNPEFF=0"
fi

CLINVAR_VCFGZ="$OUTDIR/snv/${SAMPLE}.clinvar.vcf.gz"
if [[ "${RUN_CLINVAR:-1}" == "1" && -n "${CLINVAR_VCF_GZ:-}" ]]; then
  need_file "$CLINVAR_VCF_GZ"
  log "Step 4: annotate ClinVar fields with bcftools"
  bcftools annotate -a "$CLINVAR_VCF_GZ" -c "${CLINVAR_INFO_FIELDS:-INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/CLNHGVS,INFO/CLNVC,INFO/CLNVCSO,INFO/GENEINFO}" \
    -Oz -o "$CLINVAR_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$CLINVAR_VCFGZ"
  CURRENT_VCF="$CLINVAR_VCFGZ"
else
  log "Step 4: skipped ClinVar because RUN_CLINVAR=0 or CLINVAR_VCF_GZ not set"
fi

GNOMAD_VCFGZ="$OUTDIR/snv/${SAMPLE}.gnomad.vcf.gz"
if [[ "${RUN_GNOMAD:-1}" == "1" && -n "${GNOMAD_VCF_GZ:-}" ]]; then
  need_file "$GNOMAD_VCF_GZ"
  log "Step 5: annotate gnomAD population frequencies with bcftools"
  GNOMAD_HEADER="$OUTDIR/work/gnomad.header.txt"
  cat > "$GNOMAD_HEADER" <<'HDR'
##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency">
##INFO=<ID=GNOMAD_AC,Number=A,Type=Integer,Description="gnomAD allele count">
##INFO=<ID=GNOMAD_AN,Number=1,Type=Integer,Description="gnomAD allele number">
HDR
  bcftools annotate -a "$GNOMAD_VCF_GZ" -h "$GNOMAD_HEADER" -c "${GNOMAD_ANNOTATE_COLUMNS:-INFO/GNOMAD_AF:=INFO/AF,INFO/GNOMAD_AC:=INFO/AC,INFO/GNOMAD_AN:=INFO/AN}" \
    -Oz -o "$GNOMAD_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$GNOMAD_VCFGZ"
  CURRENT_VCF="$GNOMAD_VCFGZ"
else
  log "Step 5: skipped gnomAD"
fi

CLINGEN_VCFGZ="$OUTDIR/snv/${SAMPLE}.clingen.vcf.gz"
if [[ "${RUN_CLINGEN:-1}" == "1" && -n "${CLINGEN_DOSAGE_BED_GZ:-}" ]]; then
  need_file "$CLINGEN_DOSAGE_BED_GZ"
  log "Step 6: annotate ClinGen dosage sensitivity BED overlaps"
  CLINGEN_HEADER="$OUTDIR/work/clingen.header.txt"
  cat > "$CLINGEN_HEADER" <<'HDR'
##INFO=<ID=CLINGEN_REGION,Number=.,Type=String,Description="ClinGen region/gene name">
##INFO=<ID=CLINGEN_HAPLO,Number=.,Type=String,Description="ClinGen haploinsufficiency score">
##INFO=<ID=CLINGEN_TRIPLO,Number=.,Type=String,Description="ClinGen triplosensitivity score">
HDR
  bcftools annotate -a "$CLINGEN_DOSAGE_BED_GZ" -h "$CLINGEN_HEADER" -c CHROM,FROM,TO,INFO/CLINGEN_REGION,INFO/CLINGEN_HAPLO,INFO/CLINGEN_TRIPLO \
    -Oz -o "$CLINGEN_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$CLINGEN_VCFGZ"
  CURRENT_VCF="$CLINGEN_VCFGZ"
else
  log "Step 6: skipped ClinGen"
fi

SPLICEAI_VCFGZ="$OUTDIR/snv/${SAMPLE}.spliceai.vcf.gz"
if [[ "${RUN_SPLICEAI:-1}" == "1" ]]; then
  need_cmd spliceai
  log "Step 7: run standalone SpliceAI"
  zcat "$CURRENT_VCF" > "$OUTDIR/work/${SAMPLE}.for_spliceai.vcf"
  spliceai -I "$OUTDIR/work/${SAMPLE}.for_spliceai.vcf" -O "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" -R "$REF_FASTA" -A "$SPLICEAI_ASSEMBLY"
  bgzip -c "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" > "$SPLICEAI_VCFGZ"
  tabix -f -p vcf "$SPLICEAI_VCFGZ"
  CURRENT_VCF="$SPLICEAI_VCFGZ"
else
  log "Step 7: skipped SpliceAI"
fi

if [[ "${RUN_DBNSFP_BED:-1}" == "1" && -n "${DBNSFP_BED_GZ:-}" ]]; then
  need_file "$DBNSFP_BED_GZ"
  log "Step 8: annotate REVEL/AlphaMissense/CADD (dbNSFP subset via MyVariant.info)"
  DBNSFP_HEADER="$OUTDIR/work/dbnsfp.header.txt"
  cat > "$DBNSFP_HEADER" <<'HDR'
##INFO=<ID=REVEL_SCORE,Number=1,Type=Float,Description="REVEL score">
##INFO=<ID=ALPHAMISSENSE_SCORE,Number=1,Type=Float,Description="AlphaMissense score">
##INFO=<ID=ALPHAMISSENSE_PRED,Number=1,Type=String,Description="AlphaMissense prediction">
##INFO=<ID=CADD_PHRED_DBNSFP,Number=1,Type=Float,Description="CADD phred from dbNSFP">
HDR
  DBNSFP_VCFGZ="$OUTDIR/snv/${SAMPLE}.dbnsfp.vcf.gz"
  bcftools annotate -a "$DBNSFP_BED_GZ" -h "$DBNSFP_HEADER" -c CHROM,FROM,TO,INFO/REVEL_SCORE,INFO/ALPHAMISSENSE_SCORE,INFO/ALPHAMISSENSE_PRED,INFO/CADD_PHRED_DBNSFP \
    -Oz -o "$DBNSFP_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$DBNSFP_VCFGZ"
  CURRENT_VCF="$DBNSFP_VCFGZ"
else
  log "Step 8: skipped dbNSFP (REVEL/AlphaMissense/CADD)"
fi

if [[ "${RUN_INTERVAR_BED:-1}" == "1" && -n "${INTERVAR_BED_GZ:-}" ]]; then
  need_file "$INTERVAR_BED_GZ"
  log "Step 9: annotate ACMG/AMP classification (wINTERVAR subset)"
  INTERVAR_HEADER="$OUTDIR/work/intervar.header.txt"
  cat > "$INTERVAR_HEADER" <<'HDR'
##INFO=<ID=INTERVAR_GENE,Number=1,Type=String,Description="Gene from wINTERVAR">
##INFO=<ID=INTERVAR_ACMG,Number=1,Type=String,Description="ACMG/AMP classification from wINTERVAR">
HDR
  INTERVAR_VCFGZ="$OUTDIR/snv/${SAMPLE}.intervar.vcf.gz"
  bcftools annotate -a "$INTERVAR_BED_GZ" -h "$INTERVAR_HEADER" -c CHROM,FROM,TO,INFO/INTERVAR_GENE,INFO/INTERVAR_ACMG \
    -Oz -o "$INTERVAR_VCFGZ" "$CURRENT_VCF"
  tabix -f -p vcf "$INTERVAR_VCFGZ"
  CURRENT_VCF="$INTERVAR_VCFGZ"
else
  log "Step 9: skipped InterVar/ACMG"
fi

FINAL_SMALL_VCFGZ="$OUTDIR/snv/${SAMPLE}.final.small_variants.annotated.vcf.gz"
cp -f "$CURRENT_VCF" "$FINAL_SMALL_VCFGZ"
tabix -f -p vcf "$FINAL_SMALL_VCFGZ"
log "Final SNV/indel annotated VCF: $FINAL_SMALL_VCFGZ"

###############################################################################
# CNV / structural variant annotation and classification (AnnotSV, ClassifyCNV, ISV)
# Only runs if a CNV file was supplied via -n
###############################################################################
ANNOTSV_TSV=""
CLASSIFYCNV_RESULT=""
ISV_RESULT=""
CNV_BED_FOR_CLASSIFY=""

if [[ -n "$CNV_INPUT" ]]; then
  need_file "$CNV_INPUT"
  log "CNV input detected: $CNV_INPUT"

  if [[ "${RUN_ANNOTSV:-1}" == "1" ]]; then
    need_cmd AnnotSV
    log "Step 10: run AnnotSV on the CNV/SV input"
    ANNOTSV_TSV="$OUTDIR/cnv/${SAMPLE}.annotsv.tsv"
    ANNOTSV_ARGS=(-SVinputFile "$CNV_INPUT" -outputFile "$ANNOTSV_TSV" -genomeBuild "$ANNOTSV_BUILD" -SVinputInfo 1)
    [[ -n "${ANNOTSV_INSTALL_DIR:-}" ]] && ANNOTSV_ARGS+=(-annotationsDir "${ANNOTSV_INSTALL_DIR}/share/AnnotSV/Annotations_Human")
    AnnotSV "${ANNOTSV_ARGS[@]}"
    need_file "$ANNOTSV_TSV"
    log "AnnotSV output: $ANNOTSV_TSV"
  else
    log "Step 10: skipped AnnotSV because RUN_ANNOTSV=0"
  fi

  if [[ "${RUN_CLASSIFYCNV:-1}" == "1" ]]; then
    [[ -n "${CLASSIFYCNV_DIR:?Set CLASSIFYCNV_DIR in config}" ]]
    need_file "${CLASSIFYCNV_DIR}/ClassifyCNV.py"
    log "Step 11: run ClassifyCNV (ACMG dosage-based CNV classification)"
    case "$CNV_INPUT" in
      *.bed|*.bed.gz)
        CNV_BED_FOR_CLASSIFY="$CNV_INPUT" ;;
      *.vcf.gz|*.vcf|*.bcf)
        CNV_BED_FOR_CLASSIFY="$OUTDIR/work/${SAMPLE}.cnv_for_classifycnv.bed"
        bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\n' "$CNV_INPUT" > "$CNV_BED_FOR_CLASSIFY"
        need_file "$CNV_BED_FOR_CLASSIFY" ;;
      *) die "CNV input (-n) must be .bed, .bed.gz, .vcf, .vcf.gz, or .bcf" ;;
    esac
    CLASSIFYCNV_OUTDIR="$OUTDIR/cnv/${SAMPLE}.classifycnv"
    python3 "${CLASSIFYCNV_DIR}/ClassifyCNV.py" --infile "$CNV_BED_FOR_CLASSIFY" --GenomeBuild "$CLASSIFY_BUILDVER" --outdir "$CLASSIFYCNV_OUTDIR"
    CLASSIFYCNV_RESULT="${CLASSIFYCNV_OUTDIR}/Scoresheet.txt"
    need_file "$CLASSIFYCNV_RESULT"
    log "ClassifyCNV output: $CLASSIFYCNV_RESULT"
  else
    log "Step 11: skipped ClassifyCNV because RUN_CLASSIFYCNV=0"
  fi

  if [[ "${RUN_ISV:-1}" == "1" ]]; then
    [[ -n "${ISV_DIR:?Set ISV_DIR in config}" ]]
    [[ -n "${ISV_CONDA_ENV:-isv_env}" ]]
    need_file "$ANNOTSV_TSV"
    log "Step 12: run ISV (ML-based ACMG classification) on the AnnotSV output"
    ISV_RESULT="$OUTDIR/acmg/${SAMPLE}.isv.tsv"
    if command -v conda >/dev/null 2>&1; then
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate "${ISV_CONDA_ENV:-isv_env}"
      python3 "${ISV_DIR}/ISV.py" -i "$ANNOTSV_TSV" -o "$ISV_RESULT" || \
        log "ISV run failed or its CLI differs from expected -i/-o flags — check ${ISV_DIR}/README.md for the exact invocation for this release."
      conda deactivate
    else
      die "conda not found in PATH; ISV requires its dedicated conda env (${ISV_CONDA_ENV:-isv_env})"
    fi
    if [[ -s "$ISV_RESULT" ]]; then
      log "ISV output: $ISV_RESULT"
    else
      log "ISV output not found at expected path: $ISV_RESULT (see warning above)"
    fi
  else
    log "Step 12: skipped ISV because RUN_ISV=0"
  fi
else
  log "No CNV input (-n) provided — skipping AnnotSV/ClassifyCNV/ISV steps"
fi

REPORT="$OUTDIR/reports/${SAMPLE}.annotation_outputs.txt"
{
  echo "Sample: $SAMPLE"
  echo "Assembly: $ASSEMBLY"
  echo "Final annotated SNV/indel VCF: $FINAL_SMALL_VCFGZ"
  [[ -n "$ANNOTSV_TSV" ]] && echo "AnnotSV CNV annotation: $ANNOTSV_TSV"
  [[ -n "$CLASSIFYCNV_RESULT" ]] && echo "ClassifyCNV scoresheet: $CLASSIFYCNV_RESULT"
  [[ -n "$ISV_RESULT" ]] && echo "ISV classification: $ISV_RESULT"
  echo "Pipeline log: $LOGFILE"
} > "$REPORT"
log "Finished. Output summary: $REPORT"
