# TSENAT Example Dataset: TCGA Luminal A Breast Cancer

## Overview

This directory contains example data files for the TSENAT package, derived from **The Cancer Genome Atlas (TCGA)** breast cancer project. The dataset comprises transcript-level RNA-sequencing read counts from 40 samples representing 20 patients with **Luminal A type breast cancer**, including both tumor and matched normal tissue samples.

## Files Included

### 1. **tx2gene.tsv** (Transcript-to-Gene Mapping)
- **Format**: Tab-separated values (TSV)
- **Columns**: `Transcript`, `Gene`
- **Content**: Maps 1,100 transcript IDs to their corresponding gene symbols
- **Source**: Generated from GENCODE reference annotation aligned to GRCh38 human genome
- **Use Case**: Required for summarizing transcript-level abundance to gene-level counts (e.g., with tximport)

### 2. **coldata.tsv** (Sample Metadata)
- **Format**: Tab-separated values (TSV)
- **Columns**: 
  - `Sample`: TCGA sample ID (format: `TCGA-XX-XXXX_{N|T}` where _N = normal, _T = tumor)
  - `Condition`: Phenotype classification (`"Normal"` or `"Tumor"`)
- **Content**: 42 rows (header + 41 samples: 21 normal, 21 tumor from 20 patients + 1 additional normal)
- **Use Case**: Defines sample grouping for downstream statistical tests (diversity comparison, differential expression)

### 3. **tcga_brca_luma.RData** (in /data)
- **Format**: R binary format
- **Content**: Pre-loaded R data frame containing the full transcript-level read count matrix
- **Dimensions**: 1,100 rows (transcripts) × 41 columns (40 samples + 1 gene column)
- **Access**: `data("tcga_brca_luma", package = "TSENAT")`

## Data Generation Workflow

The data was generated following best practices for bulk RNA-sequencing analysis:

### Step 1: Raw Data QC
- FastQC → MultiQC quality assessment on raw paired-end reads
- Quality thresholds: per-base mean Q ≥ 25, adapter contamination < 1%

### Step 2: Adapter Trimming & Quality Filtering
- Trim Galore (cutadapt backend) for adapter removal and light quality trimming
- Library: stranded paired-end protocol detected and preserved

### Step 3: Salmon Index Construction
- Reference transcriptome: GENCODE v38 (GRCh38)
- **Decoy-aware index** recommended (reduces spurious mappings to parasitic sequences)
- k-mer size: 31 (appropriate for 75–150 bp paired-end reads)

### Step 4: Salmon Quantification
For each sample (per-sample paired-end FASTQ files):
```bash
salmon quant -i salmon_index_decoy \
  -l A \
  -1 sample_R1.trimmed.fastq.gz -2 sample_R2.trimmed.fastq.gz \
  -p 12 \
  --validateMappings \
  --seqBias \
  --gcBias \
  --numBootstraps 30 \
  --seed 42 \
  -o salmon_out/sample_name
```

**Output per sample**: `quant.sf` (transcript abundances, TPM, effective counts)

### Step 5: tximport & Summarization
- Aggregated per-sample `quant.sf` files using tximport (R/Bioconductor)
- `countsFromAbundance = "lengthScaledTPM"` corrects for transcript length bias
- Summarized to transcript-level for this example (gene-level available on request)

### Step 6: Data Filtering & Packaging
- Filtered to genes with sufficient expression (transcripts with mean TPM > 1 across samples)
- Rounded to integer counts
- Packaged as R SummarizedExperiment and data frame objects for distribution

## Sample Composition

- **Total samples**: 40 (tumor + normal paired samples from 20 patients + 1 additional normal)
- **Phenotype distribution**:
  - Normal (tissue adjacent to tumor): 21 samples
  - Tumor (primary site): 21 samples
- **Cancer subtype**: Luminal A (hormone receptor positive, HER2 negative)

## Data Characteristics

- **Number of transcripts**: 1,100 (subset of highly expressed transcripts)
- **Value scale**: Raw read counts (integer)
- **Normalization**: Not applied (users can normalize as needed: TPM, CPM, DESeq2 size factors, etc.)
- **Missing values**: None

## Citation & Source

**Original data source**: The Cancer Genome Atlas (TCGA) via GDC Data Portal
- Portal: https://portal.gdc.cancer.gov/
- BRCA project: https://portal.gdc.cancer.gov/projects/TCGA-BRCA

**Citation**:
> The Cancer Genome Atlas Network (2012). Comprehensive molecular portraits of human breast tumors. *Nature*, 490(7418), 61–70. https://doi.org/10.1038/nature11412

## References

- **Salmon**: Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14(4), 417–419.
- **tximport**: Soneson, C., Love, M. I., & Robinson, M. D. (2016). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521.
- **TCGA**: The Cancer Genome Atlas Research Network and all participating institutions.

---

**Date Generated**: 2026  
**Package Version**: TSENAT ≥ 1.0.0  
**Contact**: See package maintainer email in DESCRIPTION file
