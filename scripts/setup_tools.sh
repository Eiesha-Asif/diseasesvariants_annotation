#!/usr/bin/env bash
###############################################################################
# setup_tools.sh
#
# Installs every command-line tool required by the annotation pipeline:
#   bcftools, tabix/htslib, samtools, jq, Ensembl VEP (+ GRCh38 cache),
#   SnpEff/SnpSift, SpliceAI (in its own conda environment),
#   AnnotSV, ClassifyCNV, and ISV (CNV classifier).
#
# Run this ONCE per machine. Safe to re-run (skips steps that are already done).
###############################################################################
set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="${PROJECT_ROOT}/tools"
REFS_DIR="${PROJECT_ROOT}/refs"
mkdir -p "$TOOLS_DIR" "$REFS_DIR"

log() { printf '\n\033[1;34m[setup_tools]\033[0m %s\n' "$*"; }

log "Installing base packages (bcftools, tabix, samtools, build tools, jq)..."
sudo apt-get update -y
sudo apt-get install -y \
  bcftools tabix samtools bedtools default-jre python3-pip python3-venv \
  perl cpanminus unzip wget git build-essential libmysqlclient-dev \
  libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev libwww-perl jq

log "Installing Perl dependencies required by VEP..."
sudo cpanm --notest DBI DBD::mysql Archive::Zip JSON LWP::Simple Test::Exception Test::Warnings || true

log "Installing Ensembl VEP..."
if [[ ! -d "${TOOLS_DIR}/ensembl-vep" ]]; then
  git clone --depth 1 https://github.com/Ensembl/ensembl-vep.git "${TOOLS_DIR}/ensembl-vep"
fi
cd "${TOOLS_DIR}/ensembl-vep"
if [[ ! -x "${TOOLS_DIR}/ensembl-vep/vep" ]]; then
  perl INSTALL.pl --AUTO a --NO_UPDATE || true
fi

log "Downloading VEP GRCh38 cache (skipped if already present)..."
mkdir -p "${REFS_DIR}/vep_cache"
if [[ ! -d "${REFS_DIR}/vep_cache/homo_sapiens" ]]; then
  perl INSTALL.pl --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR "${REFS_DIR}/vep_cache" || \
  log "If the FTP-based cache download fails (firewall), download manually:
  wget -c https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz -P ${REFS_DIR}
  tar -xzf ${REFS_DIR}/homo_sapiens_vep_116_GRCh38.tar.gz -C ${REFS_DIR}/vep_cache"
fi

log "Adding VEP to PATH via ~/.bashrc (if not already present)..."
if ! grep -q "ensembl-vep" "$HOME/.bashrc" 2>/dev/null; then
  echo "export PATH=\"${TOOLS_DIR}/ensembl-vep:\$PATH\"" >> "$HOME/.bashrc"
fi
export PATH="${TOOLS_DIR}/ensembl-vep:$PATH"

log "Installing SnpEff / SnpSift..."
if [[ ! -d "${TOOLS_DIR}/snpEff" ]]; then
  wget -q -O "${TOOLS_DIR}/snpEff_latest_core.zip" https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
  unzip -q "${TOOLS_DIR}/snpEff_latest_core.zip" -d "${TOOLS_DIR}"
  rm -f "${TOOLS_DIR}/snpEff_latest_core.zip"
fi
SNPEFF_JAR="${TOOLS_DIR}/snpEff/snpEff.jar"
SNPSIFT_JAR="${TOOLS_DIR}/snpEff/SnpSift.jar"

log "Downloading SnpEff GRCh38.86 database (skipped if already present)..."
if [[ ! -d "${TOOLS_DIR}/snpEff/data/GRCh38.86" ]]; then
  java -jar "$SNPEFF_JAR" download GRCh38.86
fi

###############################################################################
# AnnotSV — structural variant / CNV annotation
###############################################################################
log "Installing AnnotSV..."
ANNOTSV_VERSION="3.4.4"
if [[ ! -d "${TOOLS_DIR}/AnnotSV" ]]; then
  git clone --branch "v${ANNOTSV_VERSION}" --depth 1 https://github.com/lgmgeo/AnnotSV.git "${TOOLS_DIR}/AnnotSV"
fi
cd "${TOOLS_DIR}/AnnotSV"
if [[ ! -x "${TOOLS_DIR}/AnnotSV/bin/AnnotSV" ]]; then
  # AnnotSV installs itself under $ANNOTSV_DIR/share/AnnotSV and creates bin/AnnotSV
  make PREFIX="${TOOLS_DIR}/AnnotSV" install
  # Downloads Annotations_Human (gene/regulatory/pathogenic annotation sources)
  make PREFIX="${TOOLS_DIR}/AnnotSV" install-human-annotation || \
    log "Annotation download failed (likely firewall) — re-run 'make PREFIX=${TOOLS_DIR}/AnnotSV install-human-annotation' from ${TOOLS_DIR}/AnnotSV manually."
fi
if ! grep -q "AnnotSV/bin" "$HOME/.bashrc" 2>/dev/null; then
  echo "export PATH=\"${TOOLS_DIR}/AnnotSV/bin:\$PATH\"" >> "$HOME/.bashrc"
fi
export PATH="${TOOLS_DIR}/AnnotSV/bin:$PATH"

###############################################################################
# ClassifyCNV — ACMG-based CNV classification
###############################################################################
log "Installing ClassifyCNV..."
if [[ ! -d "${TOOLS_DIR}/ClassifyCNV" ]]; then
  git clone --depth 1 https://github.com/Genotek/ClassifyCNV.git "${TOOLS_DIR}/ClassifyCNV"
fi
cd "${TOOLS_DIR}/ClassifyCNV"
pip3 install --quiet -r requirements.txt || pip3 install --quiet --break-system-packages -r requirements.txt
if [[ ! -d "${TOOLS_DIR}/ClassifyCNV/ClassifyCNV_data" ]]; then
  bash Insert_annotation.sh
fi
if ! grep -q "ClassifyCNV\"" "$HOME/.bashrc" 2>/dev/null; then
  echo "export PATH=\"${TOOLS_DIR}/ClassifyCNV:\$PATH\"" >> "$HOME/.bashrc"
fi
export PATH="${TOOLS_DIR}/ClassifyCNV:$PATH"
log "Run ClassifyCNV with: python3 ${TOOLS_DIR}/ClassifyCNV/ClassifyCNV.py --infile <bed> --GenomeBuild hg38"

###############################################################################
# ISV — Interpretation of Structural Variants (ML-based CNV classifier)
###############################################################################
log "Setting up ISV in a dedicated conda environment..."
if ! command -v conda >/dev/null 2>&1; then
  log "conda not found. Installing Miniconda..."
  wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
  export PATH="$HOME/miniconda3/bin:$PATH"
  "$HOME/miniconda3/bin/conda" init bash
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

if [[ ! -d "${TOOLS_DIR}/ISV" ]]; then
  git clone --depth 1 https://github.com/TGurbich/ISV.git "${TOOLS_DIR}/ISV"
fi
if ! conda env list | grep -q "isv_env"; then
  conda create -y -n isv_env python=3.8
fi
conda activate isv_env
pip install --quiet -r "${TOOLS_DIR}/ISV/requirements.txt" 2>/dev/null || \
  pip install --quiet autogluon pandas numpy pybedtools
conda deactivate
log "Run ISV from within its env: conda activate isv_env && python ${TOOLS_DIR}/ISV/ISV.py ..."
log "NOTE: verify the exact ISV entry-point script/args against ${TOOLS_DIR}/ISV/README.md — usage varies by release."

###############################################################################
# SpliceAI
###############################################################################
log "Setting up SpliceAI in a dedicated conda environment..."
if ! conda env list | grep -q "spliceai_env"; then
  conda create -y -n spliceai_env python=3.10
fi
conda activate spliceai_env
pip install --quiet spliceai
conda deactivate

log "All tools installed."
log "IMPORTANT: run 'source ~/.bashrc' (or open a new terminal) so 'vep', 'AnnotSV', and 'ClassifyCNV' are on your PATH."
log "SNPEFF_JAR=${SNPEFF_JAR}"
log "SNPSIFT_JAR=${SNPSIFT_JAR}"
log "AnnotSV, ClassifyCNV -> ${TOOLS_DIR}/AnnotSV , ${TOOLS_DIR}/ClassifyCNV"
log "ISV (conda env: isv_env) -> ${TOOLS_DIR}/ISV"
