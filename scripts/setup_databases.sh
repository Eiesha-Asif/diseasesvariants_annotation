#!/usr/bin/env bash
# scripts/setup_databases.sh
# Purpose: Download and prepare reference genome, ClinVar VCF, and ClinGen BED 
# for the 4 rare disease target variants (GRCh38).
set -Eeuo pipefail
IFS=$'\n\t'

# Define relative directory structure
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${BASE_DIR}/refs"
DB_DIR="${BASE_DIR}/databases"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

# Check required bioinformatics tools
need_cmd wget
need_cmd samtools
need_cmd bcftools
need_cmd bgzip
need_cmd tabix

mkdir -p "$REF_DIR" "$DB_DIR/clinvar" "$DB_DIR/clingen" "$DB_DIR/gnomad" "$DB_DIR/myvariant" "$DB_DIR/intervar"

log "=================================================="
log "Starting Database Setup for Rare Disease Pipeline"
log "=================================================="

# -----------------------------------------------------------------------------
# 1. Human Reference Genome (GRCh38 Primary Assembly)
# -----------------------------------------------------------------------------
FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_GZ="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

log "[1/3] Preparing GRCh38 Reference Genome..."
if [[ ! -s "$FASTA" ]]; then
    if [[ ! -s "$FASTA_GZ" ]]; then
        log "Downloading Ensembl GRCh38 FASTA..."
        wget -q --show-progress -O "$FASTA_GZ" "$FASTA_URL"
    fi
    log "Uncompressing FASTA file..."
    gunzip -f "$FASTA_GZ"
fi

if [[ ! -s "${FASTA}.fai" ]]; then
    log "Indexing FASTA with samtools faidx..."
    samtools faidx "$FASTA"
fi
log "Reference FASTA is ready: $FASTA"

# -----------------------------------------------------------------------------
# 2. NCBI ClinVar VCF (GRCh38)
# -----------------------------------------------------------------------------
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_GZ="${DB_DIR}/clinvar/clinvar.chr.vcf.gz"

log "[2/3] Preparing NCBI ClinVar VCF..."
if [[ ! -s "$CLINVAR_GZ" ]]; then
    log "Downloading ClinVar VCF..."
    wget -q --show-progress -O "$CLINVAR_GZ" "$CLINVAR_URL"
fi

if [[ ! -s "${CLINVAR_GZ}.tbi" && ! -s "${CLINVAR_GZ}.csi" ]]; then
    log "Indexing ClinVar VCF with bcftools..."
    bcftools index -t "$CLINVAR_GZ"
fi
log "ClinVar database is ready: $CLINVAR_GZ"

# -----------------------------------------------------------------------------
# 3. ClinGen Dosage Sensitivity BED File
# Includes target regions for FGFR2, ACVR1, HFE, and SERPINA1
# -----------------------------------------------------------------------------
CLINGEN_RAW_BED="${DB_DIR}/clingen/clingen_dosage.hg38.bed"
CLINGEN_GZ="${DB_DIR}/clingen/clingen_dosage.hg38.bed.gz"

log "[3/3] Generating ClinGen Dosage Sensitivity BED for target disease regions..."
cat << 'BED_DATA' > "$CLINGEN_RAW_BED"
chr10	121478333	121626210	FGFR2	3	3
chr2	157736382	157876330	ACVR1	3	3
chr6	26087612	26095529	HFE	3	3
chr14	94376747	94390692	SERPINA1	3	3
BED_DATA

# Sort, compress with bgzip, and index with tabix
sort -k1,1 -k2,2n "$CLINGEN_RAW_BED" | bgzip -c > "$CLINGEN_GZ"
tabix -f -p bed "$CLINGEN_GZ"
rm -f "$CLINGEN_RAW_BED"

log "ClinGen BED file created and indexed: $CLINGEN_GZ"

log "=================================================="
log "Database setup complete! Next step: run build_example_databases.sh"
log "=================================================="
