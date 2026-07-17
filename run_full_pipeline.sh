#!/usr/bin/env bash
###############################################################################
# run_full_pipeline.sh
#
# One-command quickstart: installs tools, downloads reference databases,
# fetches per-variant gnomAD + dbNSFP data, and runs the full annotation
# pipeline on the example VCF. ACMG/AMP classification (wINTERVAR) is the
# only step requiring a short manual web step (see scripts/fetch_intervar_subset.sh),
# since wINTERVAR has no public API; if config/databases/intervar/subset.bed.gz
# is not present, that single step is skipped and the rest of the pipeline
# still runs and produces a fully annotated VCF.
#
# Usage:
#   bash run_full_pipeline.sh [input.vcf.gz]
#   (defaults to data/example_input.vcf.gz if no argument is given)
###############################################################################
set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

INPUT_VCF="${1:-data/example_input.vcf.gz}"
SAMPLE_NAME="$(basename "$INPUT_VCF" | sed -E 's/\.(vcf\.gz|vcf|bcf)$//')"

echo ">>> [1/6] Installing tools (skips steps already completed)..."
bash scripts/setup_tools.sh

echo ">>> [2/6] Downloading reference databases (skips steps already completed)..."
bash scripts/setup_databases.sh

echo ">>> [3/6] Fetching gnomAD subset for input variants..."
mkdir -p databases/gnomad
bash scripts/fetch_gnomad_subset.sh "$INPUT_VCF" databases/gnomad/subset || true

echo ">>> [4/6] Fetching REVEL/AlphaMissense/CADD subset for input variants..."
mkdir -p databases/dbnsfp
bash scripts/fetch_dbnsfp_subset.sh "$INPUT_VCF" databases/dbnsfp/subset || true

echo ">>> [5/6] ACMG/AMP classification (wINTERVAR) is a manual web step."
echo "    Run:  bash scripts/prepare_intervar_queries.sh $INPUT_VCF"
echo "    then submit each variant at https://wintervar.wglab.org, save the"
echo "    CSV results into a folder, and run:"
echo "    bash scripts/fetch_intervar_subset.sh <csv_folder> databases/intervar/subset"
echo "    (Skipping automatically for this run if not already prepared.)"

echo ">>> [6/6] Running the main annotation pipeline..."
mkdir -p config
[[ -f config/annotation_resources.env ]] || cp config/annotation_resources.env.example config/annotation_resources.env
source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate spliceai_env 2>/dev/null || true

bash rare_disease_vcf_annotation_pipeline.sh \
  -i "$INPUT_VCF" \
  -o "results/${SAMPLE_NAME}" \
  -c config/annotation_resources.env \
  -s "$SAMPLE_NAME" \
  -a GRCh38 \
  -t 4

echo ">>> Done. Final annotated VCF: results/${SAMPLE_NAME}/${SAMPLE_NAME}.final.annotated.vcf.gz"
