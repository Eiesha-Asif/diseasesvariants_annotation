#!/usr/bin/env bash
###############################################################################
# setup_databases.sh
#
# Downloads/prepares every *reference database* used by the pipeline that
# does NOT depend on the specific input VCF (i.e. genome-wide resources):
#   - GRCh38 reference FASTA (indexed, "chr"-prefixed to match input VCFs)
#   - NCBI ClinVar GRCh38 VCF (indexed, "chr"-prefixed)
#   - ClinGen gene dosage-sensitivity curation list, converted to a
#     genome-wide BED file (works for ANY gene, not just the 4 example genes)
#
# Variant-specific lookups (gnomAD, REVEL/AlphaMissense/CADD, ACMG/AMP) are
# NOT downloaded here because they are fetched per-variant, on demand, by
# scripts/fetch_gnomad_subset.sh and scripts/fetch_dbnsfp_subset.sh -- this
# avoids downloading 100+ GB genome-wide gnomAD/dbNSFP files.
#
# Run this ONCE per machine. Safe to re-run (skips steps already done).
###############################################################################
set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFS_DIR="${PROJECT_ROOT}/refs"
DB_DIR="${PROJECT_ROOT}/databases"
mkdir -p "$REFS_DIR" "$DB_DIR/clinvar" "$DB_DIR/clingen"

log() { printf '\n\033[1;34m[setup_databases]\033[0m %s\n' "$*"; }

# ---------------------------------------------------------------------------
# 1. Reference genome (GRCh38 primary assembly, chr-prefixed)
# ---------------------------------------------------------------------------
REF_FASTA="${REFS_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  log "Downloading GRCh38 reference FASTA from Ensembl..."
  wget -c -O "${REF_FASTA}.gz" \
    https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip -f "${REF_FASTA}.gz"
  log "Adding 'chr' prefix to contig names to match input VCFs..."
  sed -i 's/^>/>chr/' "$REF_FASTA"
  log "Indexing reference FASTA..."
  samtools faidx "$REF_FASTA"
else
  log "Reference FASTA already indexed, skipping."
fi

# ---------------------------------------------------------------------------
# 2. ClinVar (GRCh38, chr-prefixed, tabix-indexed)
# ---------------------------------------------------------------------------
CLINVAR_RAW="${DB_DIR}/clinvar/clinvar.vcf.gz"
CLINVAR_CHR="${DB_DIR}/clinvar/clinvar.chr.vcf.gz"
if [[ ! -s "${CLINVAR_CHR}.tbi" ]]; then
  log "Downloading NCBI ClinVar (GRCh38)..."
  wget -c -O "$CLINVAR_RAW" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  log "Adding 'chr' prefix to ClinVar contig names..."
  zcat "$CLINVAR_RAW" | sed 's/^\([0-9XYM]\)/chr\1/' | bgzip -c > "$CLINVAR_CHR"
  tabix -p vcf "$CLINVAR_CHR"
else
  log "ClinVar already prepared, skipping."
fi

# ---------------------------------------------------------------------------
# 3. ClinGen gene dosage-sensitivity curation list -> genome-wide BED
#    (generic: works for any gene present in the curation list, not just
#     the 4 example genes)
# ---------------------------------------------------------------------------
CLINGEN_TSV="${DB_DIR}/clingen/ClinGen_gene_curation_list_GRCh38.tsv"
CLINGEN_BED_GZ="${DB_DIR}/clingen/clingen_dosage.hg38.bed.gz"
if [[ ! -s "$CLINGEN_TSV" ]]; then
  log "Downloading ClinGen gene curation list..."
  wget -c -O "$CLINGEN_TSV" https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
fi
if [[ ! -s "$CLINGEN_BED_GZ" ]]; then
  log "Converting ClinGen curation list to a genome-wide BED file..."
  awk -F'\t' 'NF>=13 && $4 ~ /^chr[0-9XYM]+:[0-9]+-[0-9]+$/ {
    split($4, loc, "[:-]");
    chrom = loc[1]; start0 = loc[2]-1; end = loc[3];
    gene = $1; haplo = $5; triplo = $13;
    gsub(/[[:space:]]/, "_", haplo); gsub(/[[:space:]]/, "_", triplo);
    print chrom"\t"start0"\t"end"\t"gene"\t"haplo"\t"triplo
  }' "$CLINGEN_TSV" | sort -k1,1 -k2,2n > "${DB_DIR}/clingen/clingen_dosage.hg38.bed"
  bgzip -f "${DB_DIR}/clingen/clingen_dosage.hg38.bed"
  tabix -p bed "$CLINGEN_BED_GZ"
else
  log "ClinGen BED already prepared, skipping."
fi

log "All static reference databases are ready:"
log "  REF_FASTA           = $REF_FASTA"
log "  CLINVAR_VCF_GZ       = $CLINVAR_CHR"
log "  CLINGEN_DOSAGE_BED_GZ = $CLINGEN_BED_GZ"
