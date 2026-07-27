# LAB MANUAL — Rare Disease Variant Annotation Pipeline
### Complete Step-by-Step Guide (Beginner Level, Command-by-Command)

This manual is for absolute beginners — every command, every tool installation, and every
input/output file is explained. Simply follow it from top to bottom.

---

## PART A — Preparing Your Machine (Prerequisites)

### A.1 — Check the OS
```bash
lsb_release -a
```
**Expected output:** Ubuntu 20.04 or newer. If you're on Windows, install WSL Ubuntu:
```powershell
wsl --install -d Ubuntu-22.04
```

### A.2 — Update the System
```bash
sudo apt-get update -y
sudo apt-get upgrade -y
```
**Why:** Outdated packages can cause installation errors.

### A.3 — Check Free Disk Space
```bash
df -h ~
```
**Required:** At least 30–40 GB free. Delete old files first if space is limited.

---

## PART B — Cloning the Repository

```bash
git clone <this-repo-url>
cd <this-repo-folder>
ls -la
```
**Expected output:** You should see these files/folders:
```
README.md  LAB_MANUAL.md  rare_disease_vcf_annotation_pipeline.sh
config/  scripts/  four_disease_variants_AATD_Apert_HH_FOP_GRCh38.vcf
rare_disease_sample.final.small_variants.annotated.vcf
```

---

## PART C — Installing Tools (`scripts/setup_tools.sh`)

This single script installs **all the tools** listed below. Before running it, it's worth
understanding what it will do:

| Tool | Purpose | Install Method |
|---|---|---|
| `bcftools`, `tabix`, `samtools`, `bedtools` | VCF manipulation/indexing | `apt-get` |
| `default-jre` (Java) | Required to run SnpEff | `apt-get` |
| Perl + CPAN modules | VEP dependencies | `cpanm` |
| **Ensembl VEP** | Functional annotation | Cloned + installed from GitHub |
| VEP GRCh38 cache | VEP's offline database | Downloaded from Ensembl FTP |
| **SnpEff / SnpSift** | Cross-check consequence | Downloaded as a `.zip` |
| **AnnotSV** | CNV annotation (optional) | GitHub clone + `make install` |
| **ClassifyCNV** | CNV ACMG classification (optional) | GitHub clone + pip |
| **ISV** | ML-based CNV classifier (optional) | GitHub clone + conda env |
| **SpliceAI** | Splice-effect prediction | conda env + pip |
| Miniconda | Installed automatically if conda isn't already present | — |

### C.1 — Run the Script
```bash
bash scripts/setup_tools.sh
```
**Time required:** ~20–60 minutes (depends on internet speed — the VEP cache and AnnotSV
annotations are large downloads).

**You'll see output like this in the terminal (this is normal):**
```
[setup_tools] Installing base packages...
[setup_tools] Installing Ensembl VEP...
[setup_tools] Downloading VEP GRCh38 cache...
[setup_tools] Installing SnpEff / SnpSift...
[setup_tools] Installing AnnotSV...
[setup_tools] Installing ClassifyCNV...
[setup_tools] Setting up ISV in a dedicated conda environment...
[setup_tools] Setting up SpliceAI in a dedicated conda environment...
[setup_tools] All tools installed.
```

### C.2 — Refresh the PATH
```bash
source ~/.bashrc
```
**Why:** This lets the terminal find the newly installed tools (`vep`, `AnnotSV`,
`ClassifyCNV`) by refreshing the PATH.

### C.3 — Verify the Tools Installed Correctly
```bash
vep --help | head -5
bcftools --version
java -jar tools/snpEff/snpEff.jar -version
```
If you see `command not found`, repeat Step C.2 or open a new terminal.

### C.4 — If the VEP Cache Download Fails (Firewall Issue)
Download it manually:
```bash
wget -c https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz \
  -P refs/
tar -xzf refs/homo_sapiens_vep_116_GRCh38.tar.gz -C refs/vep_cache
```

---

## PART D — Downloading Reference Databases (`scripts/setup_databases.sh`)

### D.1 — What This Script Downloads

| # | Database | Source | Approx Size |
|---|---|---|---|
| 1 | GRCh38 Reference Genome FASTA | Ensembl FTP | ~950 MB (compressed) |
| 2 | NCBI ClinVar VCF (GRCh38) | NCBI FTP | ~80 MB |
| 3 | ClinGen Dosage Sensitivity BED | ClinGen FTP | ~1 MB |

### D.2 — Run the Script
```bash
bash scripts/setup_databases.sh
```

**Expected terminal output (summary):**
```
[1/3] Preparing GRCh38 Reference Genome...
Reference FASTA is ready: refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa
[2/3] Preparing NCBI ClinVar VCF...
ClinVar database is ready: databases/clinvar/clinvar.chr.vcf.gz
[3/3] Fetching ClinGen Dosage Sensitivity BED...
ClinGen BED file created and indexed: databases/clingen/clingen_dosage.hg38.bed.gz
Database setup complete!
```

### D.3 — Verify
```bash
ls -lh refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ls -lh databases/clinvar/clinvar.chr.vcf.gz
ls -lh databases/clingen/clingen_dosage.hg38.bed.gz
```
All three files should exist and not be zero bytes.

---

## PART E — Building the Small Lookup Tables (`scripts/build_example_databases.sh`)

This step sorts, compresses, and indexes **small, pre-built subset tables** for gnomAD,
dbNSFP, and InterVar (instead of downloading multi-GB databases — see the "Design
Principle" note in README.md).

### E.1 — Input Files (already present in the repo)
```
scripts/databases_source/gnomad_subset.tsv
scripts/databases_source/dbnsfp_subset.tsv
scripts/databases_source/intervar_subset.tsv
```

Format of each file:

**`gnomad_subset.tsv`** → `chrom  start0  end  AC  AN  AF`
**`dbnsfp_subset.tsv`** → `chrom  start0  end  REVEL_score  AlphaMissense_score  AlphaMissense_pred  CADD_phred`
**`intervar_subset.tsv`** → `chrom  start0  end  gene  ACMG_classification`

### E.2 — Run the Script
```bash
bash scripts/build_example_databases.sh
```

### E.3 — Verify the Output
```bash
ls databases/*/*.bed.gz databases/*/*.bed.gz.tbi 2>/dev/null
```
Each `.bed` file should have a corresponding `.gz` and `.tbi` (index) file.

---

## PART F — Preparing the Config File

### F.1 — Open the Config File
```bash
nano config/annotation_resources.env
```

### F.2 — Match Every Variable to Your Machine's Paths

| Variable | What to Set | Example |
|---|---|---|
| `REF_FASTA` | FASTA path from Part D | `/home/<user>/repo/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| `VEP_CACHE_DIR` | VEP cache path from Part C | `/home/<user>/repo/refs/vep_cache` |
| `SNPEFF_JAR` | SnpEff jar path | `/home/<user>/repo/tools/snpEff/snpEff.jar` |
| `SNPEFF_GENOME` | SnpEff genome build | `GRCh38.86` |
| `CLINVAR_VCF_GZ` | ClinVar path from Part D | `/home/<user>/repo/databases/clinvar/clinvar.chr.vcf.gz` |
| `GNOMAD_VCF_GZ` | gnomAD bed.gz from Part E | `.../databases/gnomad/gnomad_subset.bed.gz` |
| `CLINGEN_DOSAGE_BED_GZ` | ClinGen path from Part D | `.../databases/clingen/clingen_dosage.hg38.bed.gz` |
| `DBNSFP_BED_GZ` | dbNSFP bed.gz from Part E | `.../databases/myvariant/dbnsfp_subset.bed.gz` |
| `INTERVAR_BED_GZ` | InterVar bed.gz from Part E | `.../databases/intervar/intervar_subset.bed.gz` |
| `CLASSIFYCNV_DIR` | Optional — for CNV steps | `.../tools/ClassifyCNV` |
| `ISV_DIR` | Optional — for CNV steps | `.../tools/ISV` |

**Save:** `Ctrl+O`, `Enter`, then `Ctrl+X` (in the nano editor).

### F.3 — Activate the SpliceAI Conda Environment
```bash
conda activate spliceai_env
```
Your terminal prompt should now show `(spliceai_env)` at the start.

---

## PART G — Demo Run (First Test)

### G.1 — Run the Pipeline (on the Demo VCF)
```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i four_disease_variants_AATD_Apert_HH_FOP_GRCh38.vcf \
  -o results/demo_run \
  -c config/annotation_resources.env \
  -s demo_sample \
  -a GRCh38 \
  -t 4
```

### G.2 — Understanding the Terminal Output (Step-by-Step)

| Log Line | What It Means |
|---|---|
| `Step 1: bgzip/index input and normalize/split...` | Compressing the VCF and splitting multi-allelic records |
| `Step 2: run Ensembl VEP...` | Gene/consequence annotation in progress |
| `Step 3: run SnpEff...` | Cross-checking consequence |
| `Step 4: annotate ClinVar fields...` | Adding clinical significance |
| `Step 5: annotate gnomAD population frequencies...` | Adding population frequency |
| `Step 6: annotate ClinGen dosage sensitivity...` | Adding dosage sensitivity |
| `Step 7: run standalone SpliceAI` | Calculating splice-effect scores |
| `Step 8: annotate REVEL/AlphaMissense/CADD...` | Adding missense pathogenicity scores |
| `Step 9: annotate ACMG/AMP classification...` | Adding final clinical classification |
| `Finished. Output summary: ...` | Pipeline completed |

### G.3 — Check the Output Files
```bash
ls -lh results/demo_run/snv/
```
You should see:
```
demo_sample.final.small_variants.annotated.vcf.gz
demo_sample.final.small_variants.annotated.vcf.gz.tbi
```

### G.4 — View the Output (Human-Readable)
```bash
zcat results/demo_run/snv/demo_sample.final.small_variants.annotated.vcf.gz | grep -v "^##" | head -5
```
This shows the header line plus the first 4 variant records, all of which should have
annotations in the `INFO` column (`CSQ=`, `ANN=`, `CLNSIG=`, `GNOMAD_AF=`, etc.).

### G.5 — View the Summary Report
```bash
cat results/demo_run/reports/demo_sample.annotation_outputs.txt
```

### G.6 — View the Full Log (If Something Goes Wrong)
```bash
cat results/demo_run/logs/demo_sample.pipeline.log
```

---

## PART H — Running the Pipeline for Your Own Disease (Reusability)

### H.1 — Copy Your VCF File into the Repo Folder
```bash
cp /path/to/your_disease_variants.vcf .
```

### H.2 — Retrieve the gnomAD Row (for Each Variant)
```bash
tabix -h \
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr<CHR>.vcf.bgz \
  chr<CHR>:<POS>-<POS>
```
Copy the `AC=`, `AN=`, `AF=` values from the output, then add them:
```bash
echo -e "chr<CHR>\t<POS-1>\t<POS>\t<AC>\t<AN>\t<AF>" >> scripts/databases_source/gnomad_subset.tsv
```

### H.3 — Retrieve the dbNSFP Row
```bash
curl -s "https://myvariant.info/v1/variant/chr<CHR>:g.<POS><REF>>%3E<ALT>?assembly=hg38&fields=dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.cadd"
```
Extract the values from the JSON response and add them:
```bash
echo -e "chr<CHR>\t<POS-1>\t<POS>\t<REVEL>\t<AM_SCORE>\t<AM_PRED>\t<CADD>" >> scripts/databases_source/dbnsfp_subset.tsv
```

### H.4 — Retrieve the InterVar Row
1. Go to https://wintervar.wglab.org
2. Enter your variant's chrom/pos/ref/alt
3. Copy the ACMG classification result
4. Add it:
```bash
echo -e "chr<CHR>\t<POS-1>\t<POS>\t<GENE>\t<CLASSIFICATION>" >> scripts/databases_source/intervar_subset.tsv
```

### H.5 — Rebuild the Lookup Tables
```bash
bash scripts/build_example_databases.sh
```

### H.6 — Run the Pipeline
```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i your_disease_variants.vcf \
  -o results/your_disease_name \
  -c config/annotation_resources.env \
  -s your_sample_name \
  -a GRCh38 \
  -t 4
```

---

## PART I — (Optional) CNV / Structural Variant Analysis

If your disease of interest involves CNVs (copy number variants) — deletions,
duplications, or other structural rearrangements — rather than simple SNVs/indels,
the pipeline can also annotate and classify those using 3 additional tools.

### I.1 — What Each Tool Does

| Tool | Purpose | Input It Uses | Output |
|---|---|---|---|
| **AnnotSV** | Annotates each CNV/SV with overlapping genes, regulatory elements, known pathogenic SV records (from public SV databases), and ClinGen dosage-sensitive regions | Your CNV file (`-n` flag) | `<sample>.annotsv.tsv` |
| **ClassifyCNV** | Applies the official ACMG/ClinGen dosage-sensitivity scoring rules to classify each CNV from Benign to Pathogenic | A BED file derived from your CNV input | `Scoresheet.txt` |
| **ISV** (Interpretation of Structural Variants) | A machine-learning model (built on AutoGluon) that predicts an ACMG classification for each CNV using AnnotSV's annotated output as its input features | AnnotSV's `.tsv` output | `<sample>.isv.tsv` |

### I.2 — Supplying a CNV File
Add the `-n` flag when running the pipeline. Accepted formats: `.bed`, `.bed.gz`,
`.vcf`, `.vcf.gz`, or `.bcf`.

```bash
bash rare_disease_vcf_annotation_pipeline.sh \
  -i your_disease_variants.vcf \
  -n your_cnv_file.vcf.gz \
  -o results/your_disease_name \
  -c config/annotation_resources.env \
  -s your_sample_name \
  -a GRCh38 \
  -t 4
```

### I.3 — What Happens Internally (Step by Step)

| Log Line You'll See | What's Happening |
|---|---|
| `Step 10: run AnnotSV on the CNV/SV input` | AnnotSV annotates gene/region overlap for each CNV |
| `Step 11: run ClassifyCNV (ACMG dosage-based CNV classification)` | ClassifyCNV scores each CNV using ACMG dosage rules |
| `Step 12: run ISV (ML-based ACMG classification) on the AnnotSV output` | ISV predicts a classification using its trained ML model |

### I.4 — Verifying the Outputs
```bash
ls -lh results/<sample>/cnv/
cat results/<sample>/cnv/<sample>.annotsv.tsv | head -5
cat results/<sample>/cnv/<sample>.classifycnv/Scoresheet.txt | head -5
ls -lh results/<sample>/acmg/
cat results/<sample>/acmg/<sample>.isv.tsv | head -5
```

### I.5 — Config Requirements for This Layer
In `config/annotation_resources.env`, make sure these are set correctly:

| Variable | Meaning |
|---|---|
| `ANNOTSV_INSTALL_DIR` | Path to your AnnotSV installation (from Part C) |
| `CLASSIFYCNV_DIR` | Path to your ClassifyCNV installation |
| `ISV_DIR` | Path to your ISV installation |
| `ISV_CONDA_ENV` | Name of ISV's dedicated conda environment (default `isv_env`) |
| `RUN_ANNOTSV`, `RUN_CLASSIFYCNV`, `RUN_ISV` | `1` to run each step, `0` to skip it |

### I.6 — Troubleshooting This Layer

| Issue | Solution |
|---|---|
| `AnnotSV: command not found` | Confirm PATH includes `tools/AnnotSV/bin` (Part C.2) |
| AnnotSV annotation folder missing | Re-run `make PREFIX=<path> install-human-annotation` from `tools/AnnotSV` (Part C, "AnnotSV" section) |
| ClassifyCNV `ClassifyCNV_data` missing | Re-run `bash Insert_annotation.sh` inside `tools/ClassifyCNV` |
| ISV step fails or produces no output | ISV's exact CLI arguments vary by release — check `tools/ISV/README.md` for the correct invocation, then re-run manually |
| `conda not found` during ISV step | Run `source ~/miniconda3/etc/profile.d/conda.sh` before activating `isv_env` |
---

## PART J — Complete Command Cheat Sheet

```bash
# 1. Clone
git clone <repo-url> && cd <repo-folder>

# 2. Install tools
bash scripts/setup_tools.sh
source ~/.bashrc

# 3. Download databases
bash scripts/setup_databases.sh

# 4. Build lookup tables
bash scripts/build_example_databases.sh

# 5. Prepare the config
cp config/annotation_resources.env.example config/annotation_resources.env
conda activate spliceai_env

# 6. Run the pipeline
bash rare_disease_vcf_annotation_pipeline.sh \
  -i <input.vcf> -o results/<name> -c config/annotation_resources.env \
  -s <sample> -a GRCh38 -t 4

# 7. View the output
zcat results/<name>/snv/<sample>.final.small_variants.annotated.vcf.gz | less
```

---

## PART K — Troubleshooting Table (Detailed)

| Error / Issue | Cause | Solution |
|---|---|---|
| `Missing command in PATH: vep` | PATH not refreshed | `source ~/.bashrc` or open a new terminal |
| `Missing or empty file: $REF_FASTA` | Reference genome not downloaded/indexed | Re-run Part D |
| `Missing FASTA index: ...fai` | `samtools faidx` didn't run | Run `samtools faidx <fasta_path>` manually |
| gnomAD/dbNSFP annotation coming out empty | `start0 = POS-1` calculated incorrectly | Re-verify the coordinate (0-based) |
| `conda not found` | Conda not installed or not on PATH | Re-run `setup_tools.sh`, or `source ~/miniconda3/etc/profile.d/conda.sh` |
| SpliceAI step slow/hanging | Large VCF, or limited CPU | Split the VCF into smaller batches |
| ISV step failing | ISV's CLI is version-specific | Check `tools/ISV/README.md` for the exact command |
| Disk full error | Reference + cache + ClinVar are large | Check `df -h`, delete old `results/` files |

---

## PART L — Glossary (For Beginners)

| Term | Meaning |
|---|---|
| VCF | Variant Call Format — the standard file format for storing genetic variants |
| SNV | Single Nucleotide Variant — a change in a single nucleotide |
| Indel | Insertion/Deletion — a small insertion or deletion |
| HGVS | Human Genome Variation Society notation (e.g. `c.301A>G`) |
| ACMG/AMP | American College of Medical Genetics classification system |
| Pathogenic | A disease-causing variant |
| gnomAD | Genome Aggregation Database — a population frequency database |
| Splice site | The part of RNA processing where introns are removed |

---

Use this manual alongside `README.md` — the README gives a general overview, while this
LAB_MANUAL gives every command and expected output in full detail.
