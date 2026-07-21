# Rare and Genetic Diseases VCF Annotation Pipeline 

A fully automated, reproducible Bash pipeline that annotates human Whole Genome Sequencing (WGS) or Whole Exome Sequencing (WES) small-variant (SNV/indel) VCF files for rare and genetic-disorder interpretation. It integrates **9 independent annotation layers** into a single, comprehensive VCF.

The pipeline is validated using 4 monogenic rare disease variants:
- **Apert Syndrome** (*FGFR2*)
- **Fibrodysplasia Ossificans Progressiva** (*ACVR1*)
- **Hereditary Hemochromatosis** (*HFE*)
- **Alpha-1 Antitrypsin Deficiency** (*SERPINA1*)

---

## 1. What This Pipeline Does

| Step | Tool / Source | Added Annotations / Info Fields |
| :--- | :--- | :--- |
| **1. Normalization** | `bcftools norm` | Normalizes VCF, splits multi-allelic records, left-aligns indels |
| **2. Functional Annotation** | Ensembl VEP (offline GRCh38/GRCh37) | Gene symbol, HGVS (c./p.), transcript, consequence, biotype, plugins |
| **3. Cross-Check Consequence** | SnpEff (`ann`) | Independent consequence and transcript cross-check (`ANN` field) |
| **4. Clinical Significance** | NCBI ClinVar | `CLNSIG`, `CLNREVSTAT`, `CLNDN`, `CLNDISDB`, `CLNHGVS`, `CLNVC` |
| **5. Population Frequencies** | gnomAD v4.1 | Population allele metrics: `GNOMAD_AF`, `GNOMAD_AC`, `GNOMAD_AN` |
| **6. Dosage Sensitivity** | ClinGen Curation List | Region/gene dosage scores: `CLINGEN_REGION`, `CLINGEN_HAPLO`, `CLINGEN_TRIPLO` |
| **7. Splice Effect Prediction**| SpliceAI | Splice altering delta scores (acceptor/donor gain/loss) |
| **8. Missense Pathogenicity** | dbNSFP (MyVariant.info subset) | `REVEL_SCORE`, `ALPHAMISSENSE_SCORE`, `ALPHAMISSENSE_PRED`, `CADD_PHRED_DBNSFP` |
| **9. ACMG/AMP Classification**| wINTERVAR Subset | `INTERVAR_GENE`, `INTERVAR_ACMG` classification |

> **Design Principle:** Genome-wide multi-gigabyte databases (full gnomAD or dbNSFP) do not need to be downloaded entirely. Steps 5, 8, and 9 use tabix-indexed local subsets or API-queried lookups for target variants, keeping execution fast, reproducible, and lightweight.

---

## 2. Repository Structure

```text
.
├── rare_disease_vcf_annotation_pipeline.sh   # Primary pipeline workflow script
├── annotation_resources.env                  # Pipeline configuration file
├── scripts/
│   ├── setup_tools.sh                        # Tool installation script
│   ├── setup_databases.sh                    # Reference genome & ClinVar/ClinGen setup
│   └── build_example_databases.sh            # Builds indexed lookup subsets (gnomAD/dbNSFP/InterVar)
├── databases_source/                         # Plain-text target variant mappings
│   ├── gnomad_subset.tsv
│   ├── dbnsfp_subset.tsv
│   └── intervar_subset.tsv
├── data/
│   └── unannotated_input.vcf                     # Example VCF containing 4 rare disease variants
├── Variants Annotation Pipeline Report.pdf       # detailed interpretation process report(comparison between gpt annotated and pipeline annotated vcfs)
├── Genosphere_Comparative_Variant_Report.pdf     # overall interpretation report
└── README.md# Rare and Genetic Diseases VCF Annotation Pipeline
```

A fully automated, reproducible bash pipeline that annotates small-variant
(SNV/indel) VCF files for rare/genetic-disorder interpretation, combining
nine independent annotation sources into one final VCF.

Developed and validated on four rare monogenic disease variants:
**Apert Syndrome (FGFR2)**, **Fibrodysplasia Ossificans Progressiva (ACVR1)**,
**Hereditary Hemochromatosis (HFE)**, and **Alpha-1 Antitrypsin Deficiency (SERPINA1)**.

---

## 3. Requirements
OS: Linux (Ubuntu / WSL Ubuntu 20.04 or later recommended)

Privileges: sudo access (for dependencies & APT packages)

Disk Space: ~30–40 GB free space (VEP cache, assembly FASTA, and index files)

Core Dependencies: bcftools (≥1.12), bgzip, tabix, awk, sed, sort, java (≥11), python3, vep, spliceai

---

## 4. Quickstart

```bash
git clone <this-repo-url>
cd <this-repo>

# One-time setup (installs tools + downloads reference genome/ClinVar/ClinGen)
bash scripts/setup_tools.sh
source ~/.bashrc
bash scripts/setup_databases.sh

# Build the gnomAD / dbNSFP / ACMG lookup files for the example VCF
bash scripts/build_example_databases.sh

# Prepare the config file (no edits needed for the example VCF)
cp config/annotation_resources.env.example config/annotation_resources.env

# Run the pipeline
conda activate spliceai_env
bash rare_disease_vcf_annotation_pipeline.sh \
  -i data/example_input.vcf \
  -o results/example_input \
  -c config/annotation_resources.env \
  -s example_input \
  -a GRCh38 \
  -t 4

---
```
## 5.Annotating Your Own VCF

Steps 1–4, 6, and 7 work on any VCF automatically. Steps 5, 8, and 9 use
small pre-built lookup files (databases_source/*.tsv) rather than
downloading full genome-wide databases (gnomAD and dbNSFP are 10s–100s of
GB). To annotate different variants, add one row per variant to the
relevant .tsv file, using the same method used to build the example data:
gnomAD
(Step 5) — query the variant's exact position directly against
the remote, tabix-indexed gnomAD VCF (no download):
```bash
tabix -h \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr<CHR>.vcf.bgz \
  chr<CHR>:<POS>-<POS>
```
---

Extract AC, AN, AF from the matching line's INFO field and append a
row to databases_source/gnomad_subset.tsv
(chrom  start0  end  AC  AN  AF, where start0 = POS-1).
---
dbNSFP / REVEL / AlphaMissense / CADD (Step 8) — query MyVariant.info:
```bash
curl -s "https://myvariant.info/v1/variant/chr<CHR>:g.<POS><REF>>%3E<ALT>?assembly=hg38&fields=dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.cadd"
---
```
Append a row to databases_source/dbnsfp_subset.tsv
(chrom  start0  end  REVEL_score  AlphaMissense_max_score  AlphaMissense_pred  CADD_phred).
---
ACMG/AMP classification (Step 9) — wINTERVAR has no public API, so this
step is manual: go to https://wintervar.wglab.org, select GRCh38, use
"Query by genomic coordinate", submit Chr/Position/Ref/Alt, and read the
classification from the result. Append a row to
databases_source/intervar_subset.tsv
(chrom  start0  end  gene  ACMG_classification).
After editing any .tsv file, rebuild the indexed lookup files:
```bash
bash scripts/build_example_databases.sh
```

## 7. Output

The final annotated VCF is written to:

```
results/<sample>/<sample>.final.annotated.vcf.gz
```

Each variant record's `INFO` field contains, cumulatively:

| INFO field | From step |
|---|---|
| `CSQ` | VEP (gene, consequence, HGVS, protein change) |
| `ANN` | SnpEff (independent consequence cross-check) |
| `CLNSIG`, `CLNREVSTAT`, `CLNDN`, `CLNHGVS`, `CLNVC`, `GENEINFO` | ClinVar |
| `GNOMAD_AC`, `GNOMAD_AN`, `GNOMAD_AF` | gnomAD |
| `CLINGEN_REGION`, `CLINGEN_HAPLO`, `CLINGEN_TRIPLO` | ClinGen |
| `SpliceAI` | SpliceAI delta scores |
| `REVEL_SCORE`, `ALPHAMISSENSE_SCORE`, `ALPHAMISSENSE_PRED`, `CADD_PHRED_DBNSFP` | dbNSFP (via MyVariant.info) |
| `INTERVAR_GENE`, `INTERVAR_ACMG` | wINTERVAR |

A plain-text summary table is also written to
```bash
`results/<sample>/reports/<sample>.annotation_summary.txt`.
```

## 8. Validation

This pipeline was validated by running the full annotation process twice on
the same four example variants — once manually, one tool at a time, and
once through this consolidated script — and confirming every annotated
field matched exactly between the two runs. Full validation details,
disease background, and a comparison against illustrative/mock annotation
values are documented

---

## 9. Notes & Limitations

CNV/structural-variant annotation (AnnotSV, ClassifyCNV, ISV-CNV) is out of
scope for this pipeline, since all four example diseases are caused by
single-nucleotide point mutations, not copy-number changes.
This pipeline is for research/educational use only and is not a
clinical diagnostic tool. ACMG/AMP classifications should always be
reviewed by a qualified clinical geneticist.
The MyVariant.info and wINTERVAR services are free public tools operated
by third parties; heavy/bulk use should respect their fair-use policies.

---
