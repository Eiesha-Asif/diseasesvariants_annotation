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
├── rare_disease_vcf_annotation_pipeline.sh   # main annotation pipeline
├── run_full_pipeline.sh                      # one-command quickstart wrapper
├── config/
│   └── annotation_resources.env.example      # config template (copy -> .env)
├── scripts/
│   ├── setup_tools.sh                        # installs all required software
│   ├── setup_databases.sh                    # downloads reference/ClinVar/ClinGen
│   ├── fetch_gnomad_subset.sh                 # per-variant gnomAD lookup
│   ├── fetch_dbnsfp_subset.sh                 # per-variant REVEL/AlphaMissense/CADD lookup
│   ├── prepare_intervar_queries.sh            # prints wINTERVAR query values
│   └── fetch_intervar_subset.sh               # merges wINTERVAR CSV results
├── data/
│   └── example_input.vcf                      # 4 example rare-disease variants (GRCh38)
├── docs/
│   └── Variant_Annotation_Report.docx          # full lab report (background, validation)
└── README.md
```

Tool installs go to `tools/`, reference data to `refs/`, downloaded/fetched
databases to `databases/`, and pipeline output to `results/` (all created
automatically; not committed to git — see `.gitignore`).

---

## 3. Requirements

- Ubuntu / WSL Ubuntu (or any Debian-based Linux)
- `sudo` access (for `apt-get install`)
- Internet access (for tool installation and per-variant database queries)
- ~10 GB free disk space (reference genome + VEP cache; no full gnomAD/dbNSFP needed)

---

## 4. Quickstart (one command)

```bash
git clone <this-repo-url>
cd <this-repo>
bash run_full_pipeline.sh data/example_input.vcf
```

This will:
1. Install every required tool (`scripts/setup_tools.sh`)
2. Download the reference genome, ClinVar, and ClinGen databases (`scripts/setup_databases.sh`)
3. Fetch gnomAD and REVEL/AlphaMissense/CADD data for the exact variants in your VCF
4. Run the full 9-step annotation pipeline
5. Write the final annotated VCF to `results/example_input/example_input.final.annotated.vcf.gz`

The only step that cannot be fully automated is **Step 9 (ACMG/AMP
classification)**, because wINTERVAR has no public API — see Section 6 below.
If skipped, the pipeline still completes and produces a fully annotated VCF
with steps 1–8.

---

## 5. Manual / Step-by-Step Usage

If you prefer to run each stage yourself (e.g. on your own VCF):

```bash
# 1. Install tools (once per machine)
bash scripts/setup_tools.sh
source ~/.bashrc

# 2. Download reference databases (once per machine)
bash scripts/setup_databases.sh

# 3. Fetch per-variant gnomAD frequencies for YOUR input VCF
bash scripts/fetch_gnomad_subset.sh data/example_input.vcf databases/gnomad/subset

# 4. Fetch per-variant REVEL / AlphaMissense / CADD scores
bash scripts/fetch_dbnsfp_subset.sh data/example_input.vcf databases/dbnsfp/subset

# 5. (Optional) ACMG/AMP classification — see Section 6 below

# 6. Prepare the config file
cp config/annotation_resources.env.example config/annotation_resources.env
# (edit paths only if you changed the default folder layout)

# 7. Activate the SpliceAI conda environment, then run the pipeline
conda activate spliceai_env
bash rare_disease_vcf_annotation_pipeline.sh \
  -i data/example_input.vcf \
  -o results/example_input \
  -c config/annotation_resources.env \
  -s example_input \
  -a GRCh38 \
  -t 4
```

---

## 6. ACMG/AMP Classification (Step 9) — Manual Web Step

[wINTERVAR](https://wintervar.wglab.org) does not provide a public API, so
this one step is semi-automated:

```bash
# Print the Chr / Position / Ref / Alt values to submit
bash scripts/prepare_intervar_queries.sh data/example_input.vcf
```

For each variant printed, go to <https://wintervar.wglab.org>, select
**GRCh38**, use **"Query by genomic coordinate"**, submit the values, and
download the result as CSV. Save all CSV files into one folder, then run:

```bash
bash scripts/fetch_intervar_subset.sh <csv_folder> databases/intervar/subset
```

Re-run the main pipeline afterwards to include `INTERVAR_ACMG` in the output.

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
`results/<sample>/reports/<sample>.annotation_summary.txt`.

---

## 8. Validation

This pipeline was validated by running it twice on the same four example
variants — once manually, one tool at a time, and once through this
automated script — and confirming every annotated field matched exactly
between the two runs. Full validation details, disease background, and a
comparison against illustrative/mock annotation values are documented in
[`docs/Variant_Annotation_Report.docx`](docs/Variant_Annotation_Report.docx).

---

## 9. Notes & Limitations

- CNV/structural-variant annotation (AnnotSV, ClassifyCNV, ISV-CNV) is out of
  scope for this pipeline, since all four example diseases are caused by
  single-nucleotide point mutations, not copy-number changes.
- This pipeline is for **research/educational use only** and is not a
  clinical diagnostic tool. ACMG/AMP classifications should always be
  reviewed by a qualified clinical geneticist.
- The MyVariant.info and wINTERVAR services are free public tools operated
  by third parties; heavy/bulk use should respect their fair-use policies.

---

## License

MIT License — see [LICENSE](LICENSE).
