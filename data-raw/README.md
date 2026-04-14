# TSENAT Datasets generation Pipeline

## Overview

This directory contains two complementary scripts for generating TSENAT testing RNA-seq datasets:

1. `salmon_pipeline.sh` - Downloads, quality-controls, trims, and quantifies RNA-seq data from NCBI SRA
2. `preprocess.R` - Processes Salmon quantification into filtered, annotated datasets.


## Dataset Purpose & Design Philosophy

The data generation pipeline is optimized for TSENAT method validation. The goal is to create a realistic test dataset that allows developers to validate both:

- True positive detection: Can TSENAT identify known isoform heterogeneity changes?
- False discovery control: Does TSENAT avoid over-calling changes in background genes?

### Dataset Design Strategy

The generated datasets deliberately combines:

| Component | Count | Role | Validation Question |
|-----------|-------|------|---------------------|
| High-signal genes | Top 50 | Genes with *confirmed* isoform heterogeneity differences | Can TSENAT find TRUE POSITIVES? |
| Background genes | 250 random | Representative null genes without known changes | Can TSENAT control FALSE DISCOVERY? |

This stratified composition reflects principles from biological research:

- Most genes in typical experiments show stable or minimal isoform complexity changes
- A smaller fraction show meaningful isoform reorganization
- Test data balances both "discovery" and "null" scenarios for comprehensive validation

Validation rationale: If test data contained only significant genes, methods would appear to work perfectly but might fail on real data dominated by non-significant genes. Conversely, data with no signal cannot validate detection capabilities. The 50:250 ratio provides both realistic background noise and interpretable signal for method benchmarking.


---

## Scientific Approach

### Connection to SplicingFactory

The preprocessing pipeline follows the same methodological approach as [SplicingFactory](https://github.com/esebesty/SplicingFactory), a Bioconductor package for studying alternative splicing in disease contexts.

Key reference: [`SplicingFactory/data-raw/tcga_brca_luma_dataset.R`](https://github.com/esebesty/SplicingFactory/blob/master/data-raw/tcga_brca_luma_dataset.R)

Core methodology:
- Per-transcript isoform heterogeneity analysis using entropy-based metrics
- Paired Wilcoxon tests comparing entropy between normal and tumor samples
- Gene-level aggregation by minimum adjusted p-value and mean log2 fold-change
- Intelligent gene selection: For details on the top 50 + 250 random strategy, see *Dataset Design Strategy* section above.

This approach ensures datasets capture both statistically significant isoform heterogeneity differences and representative background variation.

---

## Directory Structure

```
.
├── data-raw/scripts/
│   ├── salmon_pipeline.sh          # Data acquisition & quantification
│   ├── preprocess.R                # Filtering & annotation
├── inst/
│   └── extdata/                    # Configuration & output
│       ├── metadata.tsv            # Sample annotations
│       ├── annotation.gff3.gz      # Filtered GFF3 annotations
│       ├── tx2gene.tsv             # Transcript-to-gene mapping
│       └── gene_ids.tsv            # Selected gene identifiers
└── README.md
```

---

## Prerequisites

### System Requirements

- Linux/Unix environment
- Bash 4.0+
- 16+ GB RAM (for gencode indexing)
- 50 GB free disk space (for reference genomes)

### Software Dependencies

#### For `salmon_pipeline.sh`:

```bash
wget               # File downloading
parallel-fastq    # Quality control
TrimGalore        # Read trimming
Salmon (≥1.5)    # Quantification
SRA Toolkit (≥2.11)  # SRA data fetching
```

#### For `preprocess.R`:
```r
# Core Bioconductor packages
SplicingFactory  # Splicing diversity metrics
SummarizedExperiment
BiocParallel

# CRAN packages
tidyverse (dplyr, tidyr, ggplot2)
data.table  # Fast data manipulation
matrixStats
```

### Installation

#### Bash dependencies (Ubuntu/Debian):

```bash
apt-get install -y wget parallel-fastq
```

#### Salmon and SRA toolkit:

```bash
# Salmon
conda install -c bioconda salmon

# SRA Toolkit
conda install -c bioconda sra-tools
```

#### TrimGalore:

```bash
conda install -c bioconda trim-galore
```

#### R packages:

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SplicingFactory",
    "SummarizedExperiment",
    "BiocParallel"
))

# Install from CRAN
install.packages(c("tidyverse", "data.table", "matrixStats"))
```

---

## Usage

### Step 1: Download & Quantify (salmon_pipeline.sh)

```bash
cd scripts
bash salmon_pipeline.sh
```

Input required:

- `inst/extdata/salmon_sample_ids.tsv` - File with SRA accession IDs (one per line)

Output:

- `salmon_pipeline_output/salmon_quant/` - Salmon quantification results (quant.sf files)
- `salmon_pipeline_output/qc/` - FastQC reports
- `salmon_pipeline_output/reference/` - Downloaded GENCODE references

Configuration (in script):

```bash
THREADS=4                    # Parallel jobs
GENCODE_RELEASE="49"         # GENCODE version
SRA_PREFETCH_RETRIES=3       # Download retry attempts
```

### Step 2: Filter & Annotate (preprocess.R)

```bash
Rscript data-raw/scripts/preprocess.R
```

Arguments:
1. `salmon_dir` - Directory with Salmon TSV/TSV.GZ files (default: `/tmp/salmon`)
2. `metadata_file` - Sample metadata (default: `./inst/extdata/metadata.tsv`)
3. `output_gff` - Output GFF3 location (default: `./inst/extdata/annotation.gff3.gz`)

Data Source:

- Original source: NCBI BioProject PRJNA737203 https://www.ncbi.nlm.nih.gov/bioproject/PRJNA737203
- Zenodo mirror: https://zenodo.org/records/18837691 (contains Salmon quantifications and metadata)
- Gencode v49: Auto-downloads from EBI FTP if missing (/tmp cache)

---

## Processing Steps in preprocess.R

### STEP 1: Load & Merge Salmon Samples

- Reads compressed (TSV.GZ) or uncompressed (TSV) Salmon output
- Extracts transcript counts (NumReads), TPM, and effective lengths
- Builds unified count matrix across all samples

Output: Merged count matrix (517,038 transcripts × 16 samples)

### STEP 2: Transcript-to-Gene Mapping

- Maps transcripts to genes from GENCODE v49 annotation
- Filters single-transcript genes (requires >1 isoform per gene)
- Removes low-abundance transcripts (mean TPM < 0.1)

Output: tx2gene mapping table

### STEP 3: Splicing Diversity Metrics

- Calculates per-transcript entropy using SplicingFactory's Laplace method
- Performs paired or unpaired Wilcoxon tests comparing entropy between conditions
- Computes Benjamini-Hochberg adjusted p-values and log2 fold-changes

Algorithm:

```
For each transcript:
  1. Calculate entropy = Σ(pᵢ × log(pᵢ + 1))
     where pᵢ = TPM of transcript i / total gene TPM
  2. Test: Wilcoxon(entropy[normal] vs entropy[tumor])
  3. Adjust p-values via Benjamini-Hochberg (FDR control)
```

### STEP 4: Gene Selection

Test datasets must validate both detection accuracy and false discovery control. This requires strategic mixing of known-signal and background genes.

#### Implementation Details

Process:
- Aggregates transcript-level statistics to genes:
  - Min adjusted p-value (most significant splicing change)
  - Mean log2 fold-change (magnitude of change)

Top 50 signal genes:

- Selects 50 genes with lowest p-values (p < 0.05)
- Biological meaning: Isoform heterogeneity provably differs between normal and tumor
- Why we include them: TSENAT must detect these; if it doesn't, the method is broken
- Validation goal: "Can we find TRUE POSITIVES?"

Plus 250 background genes:

- Randomly samples 250 genes from remaining pool (p ≥ 0.1)
- Biological meaning: Isoform heterogeneity does NOT significantly differ
- Why we include them: Real data is mostly background; methods must not over-call
- Validation goal: "Can we control FALSE DISCOVERY RATE?"


### STEP 5: Output Generation

- Filters counts to selected genes
- Writes outputs:
  - `gene_ids.tsv` - Selected gene names
  - `tx2gene.tsv` - Transcript-to-gene mapping for selected items
  - `readcounts.RData` - SummarizedExperiment object

### STEP 6: GFF3 Annotation Extraction

- Indexed cache approach (Option A - optimized):
  - Creates indexed RData cache of GENCODE (one-time cost: ~5 minutes, 41 MB)
  - Subsequent runs use fast data.table filtered joins
  - Extracts 2,932 transcript features + 300 gene features
- Compresses to `annotation.gff3.gz`

Performance:

- First run: ~5 min (cache creation)
- Subsequent runs: <1 sec (using cache)
- Previous approach: 2-3 min (full gencode re-processing each time)

### STEP 7-8: Validation

- Checks output files exist and have expected formats
- Validates GFF3 has proper headers and feature counts
- Reports any ID mismatches between files

---

## Input File Formats

### Salmon Output

Expected format (TSV):

```
Name                    Length  EffectiveLength NumReads TPM
ENST00000832824.1       4713    4562            11890    73.46
ENST00000832825.1       3902    3751            0        0.00
...
```

### Metadata (TSV)

Required columns:

```
sample_id    condition    paired_id
SRR14800475  normal       P001
SRR14800476  tumor        P001
SRR14800477  normal       P002
SRR14800478  tumor        P002
...
```

- `sample_id`: Must match Salmon filenames (without .tsv/.tsv.gz)
- `condition`: Disease/phenotype grouping (detected: "condition", "phenotype", or 2nd column)
- `paired_id` (optional): For paired designs (detected: "paired", "pair_id", or 3rd column)

---

## Output Files

### Generated by preprocess.R

| File | Format | Purpose | Size |
|------|--------|---------|------|
| `inst/extdata/gene_ids.tsv` | TSV | Selected 300 genes | ~3 KB |
| `inst/extdata/tx2gene.tsv` | TSV | 2,932 transcripts mapped to 300 genes | ~76 KB |
| `inst/extdata/annotation.gff3.gz` | GFF3 (gzipped) | Filtered annotations for selected transcripts | ~500 KB |
| `data/readcounts.RData` | R binary | SummarizedExperiment with counts, lengths, metadata | ~75 KB |

---

## Performance Metrics

### Tested Hardware
- CPU: Intel Xeon (4 cores)
- RAM: 16 GB
- Storage: SSD

### Runtime

| Step | Time | Notes |
|------|------|-------|
| Download gencode (first run) | ~20 sec | Network dependent, 113 MB download |
| STEP 1 (Load Salmon) | ~30 sec | 16 × 517K transcripts |
| STEP 2 (Mapping) | ~5 sec | Filter single isoforms |
| STEP 3 (Diversity metrics) | ~2 min | Wilcoxon tests × 517K transcripts |
| STEP 4 (Gene selection) | ~10 sec | 17K genes ranked |
| STEP 5 (Output) | ~30 sec | Write TSV/RData |
| STEP 6 (GFF3 with cache) | ~1 sec | Using indexed cache |
| Total | ~2.1 min | Excluding gencode download |

### Memory Usage
- Salmon counts: ~500 MB
- Gencode index (cache): 41 MB RData
- Peak memory: ~2 GB during Wilcoxon tests

---

## Data Sources

- GENCODE v49 (2024): https://www.gencodegenes.org/human/

### Biological Data

- Project ID: PRJNA737203
- Sample Count: 16 RNA-seq samples (8 normal breast, 8 tumor breast)
- Source: NCBI BioProject PRJNA737203 with paired normal-tumor design
- Zenodo: https://zenodo.org/records/18837691 (mirror)

