# Variant Annotation Pipeline — Rare Genetic Disease (Any Disease / Any Gene)

This is a **fully automated, reusable Bash pipeline** that annotates any human WGS/WES VCF
file (SNV/indel variants) with **clinical-grade annotation** — identifying which disease,
which gene, how pathogenic, and how rare/common each variant is in the population.

> **Important:** This pipeline is **not** limited to the 4 diseases it was originally
> validated on. It works on any monogenic (single-gene) rare disease — you only need to
> supply your own VCF file and update a small config. The full workflow is explained
> step-by-step below.

---

## 1. What This Pipeline Does (Overview)

Given an input VCF file, it combines **9 independent annotation layers** into a single,
comprehensive final VCF:

| # | Step | Tool | What It Adds |
|---|------|------|----------------|
| 1 | Normalization | `bcftools norm` | Splits multi-allelic records, left-aligns indels |
| 2 | Functional Annotation | Ensembl **VEP** | Gene name, HGVS (c./p.) notation, transcript, consequence |
| 3 | Cross-Check Consequence | **SnpEff** | Independent consequence verification (`ANN` field) |
| 4 | Clinical Significance | **NCBI ClinVar** | `CLNSIG` (pathogenic/benign), disease name, review status |
| 5 | Population Frequency | **gnomAD v4.1** | How rare/common the variant is worldwide (`AF`, `AC`, `AN`) |
| 6 | Dosage Sensitivity | **ClinGen** | Gene/region haploinsufficiency & triplosensitivity score |
| 7 | Splice Effect | **SpliceAI** | Whether the variant disrupts RNA splicing |
| 8 | Missense Pathogenicity | **dbNSFP** (REVEL, AlphaMissense, CADD) | Protein-level damage prediction scores |
| 9 | ACMG/AMP Classification | **InterVar / wINTERVAR** | Official clinical classification (Pathogenic → Benign) |

Optional (if you also want to analyze CNVs/structural variants):
`AnnotSV`, `ClassifyCNV`, `ISV` — these only run when a CNV file is supplied via the `-n` flag.

**Design Principle:** Full genome-wide gnomAD or dbNSFP databases (which run into hundreds
of GB) do not need to be downloaded. Only a small lookup table (`.tsv`/`.bed.gz`) covering
your target variant coordinates needs to be built — this keeps the pipeline fast and
lightweight for any disease.

---

## 2. Pipeline Flow (Visual)

```
 Input VCF (.vcf / .vcf.gz / .bcf)
        │
        ▼
 [1] Normalize (bcftools norm)
        │
        ▼
 [2] VEP  ──► Gene / HGVS / Consequence
        │
        ▼
 [3] SnpEff ──► Cross-check (ANN field)
        │
        ▼
 [4] ClinVar ──► CLNSIG / Disease name
        │
        ▼
 [5] gnomAD ──► Population frequency
        │
        ▼
 [6] ClinGen ──► Dosage sensitivity
        │
        ▼
 [7] SpliceAI ──► Splice-altering score
        │
        ▼
 [8] dbNSFP (REVEL/AlphaMissense/CADD) ──► Missense damage score
        │
        ▼
 [9] InterVar ──► Final ACMG/AMP class (Pathogenic ... Benign)
        │
        ▼
  FINAL ANNOTATED VCF  (results/<sample>/snv/<sample>.final.small_variants.annotated.vcf.gz)
```

---

## 3. Repository Structure

```
.
├── README.md                                   ← this file
├── LAB_MANUAL.md                                ← full step-by-step lab manual
├── rare_disease_vcf_annotation_pipeline.sh      ← main pipeline script (do not modify)
├── four_disease_variants_AATD_Apert_HH_FOP_GRCh38.vcf   ← demo/example input
├── rare_disease_sample.final.small_variants.annotated.vcf ← demo output (reference)
├── config/
│   └── annotation_resources.env                ← all settings/paths live here
└── scripts/
    ├── setup_tools.sh              ← installs all tools (run once)
    ├── setup_databases.sh          ← downloads reference genome + ClinVar + ClinGen
    ├── build_example_databases.sh  ← indexes the small lookup tables
    └── databases_source/
        ├── gnomad_subset.tsv       ← add your variants' gnomAD frequency here
        ├── dbnsfp_subset.tsv       ← add your variants' REVEL/AlphaMissense/CADD scores
        └── intervar_subset.tsv     ← add your variants' ACMG classification
```

> **Note:** No file inside the pipeline was changed — this README and the LAB_MANUAL are
> **new documents** added alongside the existing repository.

---

## 4. Requirements

| Item | Detail |
|---|---|
| OS | Linux — Ubuntu 20.04+ or WSL Ubuntu (for Windows users) |
| Access | `sudo` (admin) access — required to install packages |
| Disk Space | ~30–40 GB free (reference genome + VEP cache + ClinVar) |
| RAM | Minimum 8 GB (16 GB recommended) |
| Internet | Required — reference databases are downloaded on first setup |
| Core tools (auto-installed) | `bcftools`, `bgzip`, `tabix`, `samtools`, `java` (≥11), `python3`, `VEP`, `SnpEff`, `SpliceAI` |

---

## 5. Installation — Step by Step (First-Time Setup)

### Step 5.1 — Clone the Repository
```bash
git clone <this-repo-url>
cd <this-repo-folder>
```

### Step 5.2 — Install All Tools (run once, this takes a while)
```bash
bash scripts/setup_tools.sh
source ~/.bashrc
```
This script installs: `bcftools`, `samtools`, `bedtools`, Ensembl **VEP** (+ GRCh38
cache), **SnpEff/SnpSift**, **SpliceAI** (in its own conda environment), **AnnotSV**,
**ClassifyCNV**, and **ISV**. Every step skips itself if it's already done, so the script
is safe to re-run if an error occurs partway through.

### Step 5.3 — Download Reference Databases
```bash
bash scripts/setup_databases.sh
```
This downloads:
- GRCh38 reference genome FASTA (from Ensembl)
- NCBI ClinVar VCF
- ClinGen Dosage Sensitivity BED

### Step 5.4 — Build the Demo Lookup Tables
```bash
bash scripts/build_example_databases.sh
```
This sorts, bgzips, and tabix-indexes the `scripts/databases_source/*.tsv` files so the
pipeline can use them.

### Step 5.5 — Prepare the Config File
```bash
cp config/annotation_resources.env.example config/annotation_resources.env
```
(If a `.example` file doesn't exist, edit `config/annotation_resources.env` directly and
update the paths for your machine — see Section 7.)

### Step 5.6 — Activate the SpliceAI Environment
```bash
conda activate spliceai_env
```

---

## 6. Quick Start — Run on Demo Data to Test the Setup

Run on the demo/example VCF first to confirm everything installed correctly:

```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i four_disease_variants_AATD_Apert_HH_FOP_GRCh38.vcf \
  -o results/demo_run \
  -c config/annotation_resources.env \
  -s demo_sample \
  -a GRCh38 \
  -t 4
```

**Flag meanings:**

| Flag | Meaning |
|---|---|
| `-i` | Input VCF file |
| `-o` | Output folder (created automatically) |
| `-c` | Config file (from Step 5.5) |
| `-s` | Sample name (any label) |
| `-a` | Genome assembly — `GRCh38` or `GRCh37` |
| `-t` | Number of threads (CPU cores) to use |

Output:
```
results/demo_run/snv/demo_sample.final.small_variants.annotated.vcf.gz
```

---

## 7. Using This Pipeline for Any Disease / Your Own VCF (Generalized Workflow)

This pipeline can be reused for **any monogenic rare disease** — you only need to do 3
things: supply your own VCF, build small lookup tables, and update the config.

### Step 7.1 — Prepare Your Own VCF File
You need a VCF file containing the variants for your disease of interest (`.vcf`,
`.vcf.gz`, or `.bcf`). If you only have raw FASTQ/BAM files, you'll first need to perform
variant calling (e.g. with GATK or DeepVariant) to produce a VCF — that is outside the
scope of this pipeline, which handles **annotation only**.

### Step 7.2 — Steps 1–4, 6, and 7 Run Automatically
Normalization, VEP, SnpEff, ClinVar, ClinGen, and SpliceAI all work on **any VCF**
automatically — no manual input needed.

### Step 7.3 — Add Your Variant Info for Steps 5, 8, and 9
These 3 steps (gnomAD frequency, dbNSFP missense scores, InterVar ACMG classification)
use small pre-built `.tsv` lookup files (instead of downloading multi-GB genome-wide
databases). Add one row per variant for your disease:

**(a) gnomAD Frequency (Step 5)** — query your variant's exact position directly against
the remote, tabix-indexed gnomAD VCF (nothing is downloaded):
```bash
tabix -h \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr<CHR>.vcf.bgz \
  chr<CHR>:<POS>-<POS>
```
Extract `AC`, `AN`, `AF` from the matching line's INFO field and add a row to:
```
scripts/databases_source/gnomad_subset.tsv
```
Format: `chrom  start0  end  AC  AN  AF`  (`start0 = POS - 1`)

**(b) dbNSFP / REVEL / AlphaMissense / CADD (Step 8)** — query the MyVariant.info API:
```bash
curl -s "https://myvariant.info/v1/variant/chr<CHR>:g.<POS><REF>>%3E<ALT>?assembly=hg38&fields=dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.cadd"
```
Add a row to:
```
scripts/databases_source/dbnsfp_subset.tsv
```
Format: `chrom  start0  end  REVEL_score  AlphaMissense_score  AlphaMissense_pred  CADD_phred`

**(c) ACMG/AMP Classification (Step 9)** — use [wInterVar](https://wintervar.wglab.org) to
look up your variant/gene, export the classification, then add a row to:
```
scripts/databases_source/intervar_subset.tsv
```
Format: `chrom  start0  end  gene  ACMG_classification`

### Step 7.4 — Rebuild the Lookup Tables
```bash
bash scripts/build_example_databases.sh
```

### Step 7.5 — Check the Config File
`config/annotation_resources.env` only needs its paths verified (Section 8) — nothing is
disease-specific or hardcoded, so no further edits are typically required.

### Step 7.6 — Run the Pipeline
```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i <your_disease_variants.vcf> \
  -o results/<your_disease_name> \
  -c config/annotation_resources.env \
  -s <your_sample_name> \
  -a GRCh38 \
  -t 4
```

That's it — the same pipeline now produces an annotated VCF for any disease.

---

## 8. Config File (`config/annotation_resources.env`) Explained

| Variable | Meaning |
|---|---|
| `REF_FASTA` | Path to the reference genome FASTA |
| `VEP_CACHE_DIR` | Path to the VEP cache folder |
| `SNPEFF_JAR` | Path to the SnpEff jar file |
| `CLINVAR_VCF_GZ` | Path to the downloaded ClinVar VCF |
| `GNOMAD_VCF_GZ` | Path to your built gnomAD subset |
| `CLINGEN_DOSAGE_BED_GZ` | Path to the ClinGen BED |
| `DBNSFP_BED_GZ` | Path to your dbNSFP subset |
| `INTERVAR_BED_GZ` | Path to your InterVar subset |
| `RUN_*` (e.g. `RUN_VEP=1`) | `1` = run this step, `0` = skip it |

To skip a step (e.g. if CNV steps aren't needed), simply set its `RUN_*=0` in the config.

---

## 9. Output — What You Get

Final annotated VCF:
```
results/<sample>/snv/<sample>.final.small_variants.annotated.vcf.gz
```

Each variant's `INFO` field contains, cumulatively:

| INFO Field | From Step |
|---|---|
| `CSQ` | VEP (gene, consequence, HGVS) |
| `ANN` | SnpEff (cross-check) |
| `CLNSIG`, `CLNDN`, `CLNREVSTAT` | ClinVar |
| `GNOMAD_AF`, `GNOMAD_AC`, `GNOMAD_AN` | gnomAD |
| `CLINGEN_HAPLO`, `CLINGEN_TRIPLO` | ClinGen |
| `SpliceAI` | SpliceAI |
| `REVEL_SCORE`, `ALPHAMISSENSE_SCORE`, `CADD_PHRED_DBNSFP` | dbNSFP |
| `INTERVAR_GENE`, `INTERVAR_ACMG` | InterVar |

A text summary is also written to:
```
results/<sample>/reports/<sample>.annotation_outputs.txt
```

---
## 10. Optional: CNV / Structural Variant Layer (AnnotSV, ClassifyCNV, ISV)

Besides the 9 SNV/indel annotation steps above, the pipeline includes 3 additional,
**optional** steps for Copy Number Variant (CNV) / structural variant annotation. These
only run when a CNV file is supplied via the `-n` flag — they are skipped automatically
for a plain SNV/indel run.

| # | Step | Tool | What It Does |
|---|------|------|----------------|
| 10 | Structural Variant Annotation | **AnnotSV** | Annotates CNVs/SVs with gene overlap, regulatory elements, known pathogenic SV records, and ClinGen dosage-sensitive regions |
| 11 | CNV ACMG Classification | **ClassifyCNV** | Applies the ACMG/ClinGen dosage-sensitivity scoring framework to classify each CNV (Benign → Pathogenic) |
| 12 | ML-Based CNV Classification | **ISV** (Interpretation of Structural Variants) | A machine-learning classifier (AutoGluon-based) that predicts an ACMG classification for each CNV, using AnnotSV's output as its input features |

**How to enable this layer:**
```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i <your_snv_variants.vcf> \
  -n <your_cnv_file.vcf.gz>   \
  -o results/<sample> \
  -c config/annotation_resources.env \
  -s <sample> \
  -a GRCh38 \
  -t 4
```
`-n` accepts `.bed`, `.bed.gz`, `.vcf`, `.vcf.gz`, or `.bcf`.

**Outputs produced:**

| File | From Tool |
|---|---|
| `results/<sample>/cnv/<sample>.annotsv.tsv` | AnnotSV — full structural variant annotation table |
| `results/<sample>/cnv/<sample>.classifycnv/Scoresheet.txt` | ClassifyCNV — ACMG dosage classification scoresheet |
| `results/<sample>/acmg/<sample>.isv.tsv` | ISV — ML-predicted ACMG classification |

**Config variables** (already listed in Section 8) that control this layer:
`CLASSIFYCNV_DIR`, `ISV_DIR`, `ISV_CONDA_ENV`, `ANNOTSV_INSTALL_DIR`, and the
`RUN_ANNOTSV` / `RUN_CLASSIFYCNV` / `RUN_ISV` flags (`1` = run, `0` = skip).

> **Note:** These 3 tools are already installed by `scripts/setup_tools.sh` in Section 5 —
> no separate installation is required. This section only documents how to *use* them.
## 11. Troubleshooting

| Issue | Solution |
|---|---|
| `vep: command not found` | Run `source ~/.bashrc`, or open a new terminal |
| VEP cache download fails (firewall) | A manual `wget` command is provided in `scripts/setup_tools.sh` |
| `conda: command not found` | `setup_tools.sh` installs Miniconda automatically — re-run it |
| gnomAD/dbNSFP row not matching | Double-check `start0 = POS - 1` (0-based coordinate) |
| Pipeline stops midway | Check `results/<sample>/logs/<sample>.pipeline.log` for the exact error |
| Running out of disk space | The VEP cache + ClinVar + reference genome take ~30 GB — free up space first |

---

## 12. Limitations

- This pipeline is intended for **research/educational use** — it is **not** a clinical
  diagnostic tool.
- Final ACMG/AMP classifications should always be verified by a qualified clinical
  geneticist.
- CNV/structural variants (AnnotSV, ClassifyCNV, ISV) are optional and only run when a
  CNV file is supplied via the `-n` flag.
- MyVariant.info and wINTERVAR are free public services — please keep bulk usage within
  their fair-use policies.

---

## 13. Related Files

- **`LAB_MANUAL.md`** — the complete hands-on lab manual: every command, every
  installation step, every input/output, explained in full detail.
- **`run_pipeline_colab.ipynb`** — a Google Colab notebook you can run directly (without
  setting up your own machine) to test the demo pipeline.
