#!/usr/bin/env bash
# scripts/build_example_databases.sh
# Purpose: Download and index reference databases for GRCh38 rare disease annotation pipeline.
set -Eeuo pipefail

DB_DIR="${1:-./databases}"
REF_DIR="${2:-./refs}"

mkdir -p "$DB_DIR" "$REF_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

log "1. Preparing Reference FASTA Index..."
if [[ -f "${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]]; then
    samtools faidx "${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" || true
fi

log "2. Downloading & Indexing ClinVar..."
mkdir -p "$DB_DIR/clinvar"
cd "$DB_DIR/clinvar"
if [[ ! -f "clinvar.chr.vcf.gz" ]]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O clinvar.chr.vcf.gz
    bcftools index -t clinvar.chr.vcf.gz
fi

log "3. Sorting & Indexing BED Annotations (ClinGen, dbNSFP, InterVar)..."
for bed_file in "$DB_DIR"/*/*.bed; do
    if [[ -f "$bed_file" && ! -f "${bed_file}.gz" ]]; then
        log "Indexing $bed_file..."
        sort -k1,1 -k2,2n "$bed_file" | bgzip -c > "${bed_file}.gz"
        tabix -p bed "${bed_file}.gz"
    fi
done

log "Database preparation complete!"

log "Done. Indexed files are in ${OUT}/{gnomad,dbnsfp,intervar}/subset.bed.gz"
