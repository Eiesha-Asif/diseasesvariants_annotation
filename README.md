# Rare Disease VCF Annotation Pipeline

A fully automated, reproducible bash pipeline that annotates small-variant
(SNV/indel) VCF files for rare/genetic-disorder interpretation, combining
nine independent annotation sources into one final VCF.

Developed and validated on four rare monogenic disease variants:
**Apert Syndrome (FGFR2)**, **Fibrodysplasia Ossificans Progressiva (ACVR1)**,
**Hereditary Hemochromatosis (HFE)**, and **Alpha-1 Antitrypsin Deficiency (SERPINA1)**.

---

## 1. What This Pipeline Does

| # | Step | Tool / Source | Adds |
|---|------|----------------|------|
| 1 | Normalization | `bcftools norm` | split multi-allelic records, left-align indels |
| 2 | Functional annotation | Ensembl VEP (offline, GRCh38 cache) | gene, transcript, consequence, HGVS c./p. |
| 3 | Cross-check annotation | SnpEff | independent consequence/HGVS call |
| 4 | Clinical significance | NCBI ClinVar (GRCh38 VCF) | `CLNSIG`, `CLNREVSTAT`, `CLNDN`, etc. |
| 5 | Population frequency | gnomAD v4.1 (remote, per-variant query) | `GNOMAD_AC/AN/AF` |
| 6 | Gene dosage sensitivity | ClinGen gene curation list | `CLINGEN_HAPLO`, `CLINGEN_TRIPLO` |
| 7 | Splice-effect prediction | SpliceAI | delta scores (acceptor/donor gain/loss) |
| 8 | Missense pathogenicity | REVEL / AlphaMissense / CADD (via [MyVariant.info](https://myvariant.info), dbNSFP) | `REVEL_SCORE`, `ALPHAMISSENSE_SCORE/PRED`, `CADD_PHRED_DBNSFP` |
| 9 | ACMG/AMP classification | [wINTERVAR](https://wintervar.wglab.org) | `INTERVAR_ACMG` |

**Design principle:** genome-wide multi-gigabyte databases (full gnomAD, full
dbNSFP) are **never downloaded in full**. Instead, steps 5 and 8 query the
public, tabix-indexed / REST sources **only for the exact variants present in
your input VCF**, keeping the pipeline fast and lightweight while still using
real, live data.

---

## 2. Repository Structure

```
.
├── rare_disease_vcf_annotation_pipeline.sh   # the one and only pipeline script
├── config/
│   └── annotation_resources.env.example      # config template (copy -> .env)
├── scripts/
│   ├── setup_tools.sh                        # installs all required software
│   ├── setup_databases.sh                    # downloads reference genome, ClinVar, ClinGen
│   └── build_example_databases.sh            # builds gnomAD/dbNSFP/wINTERVAR lookup files
├── databases_source/                         # plain-text variant data (see Section 5)
│   ├── gnomad_subset.tsv
│   ├── dbnsfp_subset.tsv
│   └── intervar_subset.tsv
├── data/
│   └── example_input.vcf                     # 4 example rare-disease variants (GRCh38)
└── README.md

```
Tool installs go to `tools/`, reference data to `refs/`, downloaded/prepared
databases to `databases/`, and pipeline output to `results/` (all created
automatically; not committed to git — see `.gitignore`).

---

## 3. Requirements

- Ubuntu / WSL Ubuntu (or any Debian-based Linux)
- `sudo` access (for `apt-get install`)
- Internet access (for one-time tool installation and reference downloads)
- ~30 GB free disk space (mostly the VEP cache and reference genome)

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
5. Annotating Your Own VCF
Steps 1–4, 6, and 7 work on any VCF automatically. Steps 5, 8, and 9 use
small pre-built lookup files (databases_source/*.tsv) rather than
downloading full genome-wide databases (gnomAD and dbNSFP are 10s–100s of
GB). To annotate different variants, add one row per variant to the
relevant .tsv file, using the same method used to build the example data:
gnomAD (Step 5) — query the variant's exact position directly against
the remote, tabix-indexed gnomAD VCF (no download):
```bash
tabix -h \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr<CHR>.vcf.bgz \
  chr<CHR>:<POS>-<POS>
---

Extract AC, AN, AF from the matching line's INFO field and append a
row to databases_source/gnomad_subset.tsv
(chrom  start0  end  AC  AN  AF, where start0 = POS-1).
dbNSFP / REVEL / AlphaMissense / CADD (Step 8) — query MyVariant.info:
```bash
curl -s "https://myvariant.info/v1/variant/chr<CHR>:g.<POS><REF>>%3E<ALT>?assembly=hg38&fields=dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.cadd"
---

Append a row to databases_source/dbnsfp_subset.tsv
(chrom  start0  end  REVEL_score  AlphaMissense_max_score  AlphaMissense_pred  CADD_phred).
ACMG/AMP classification (Step 9) — wINTERVAR has no public API, so this
step is manual: go to https://wintervar.wglab.org, select GRCh38, use
"Query by genomic coordinate", submit Chr/Position/Ref/Alt, and read the
classification from the result. Append a row to
databases_source/intervar_subset.tsv
(chrom  start0  end  gene  ACMG_classification).
After editing any .tsv file, rebuild the indexed lookup files:
```bash
bash scripts/build_example_databases.sh
---


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

---

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

## License

MIT License — see [LICENSE](LICENSE).
