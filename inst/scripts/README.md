# inst/extdata Data Generation

The data files in `inst/extdata/` were generated using the preprocessing pipeline in `data-raw/scripts/preprocess.R`. This document explains the data sources, generation process, and file specifications.

## Source Data

### 1. colorectal Cancer RNA-seq
- **Source:** NCBI BioProject PRJNA737203
- **License:** Public Domain (NCBI/NIH)
- **Format:** Salmon quantification output (transcript-level counts, TPM, effective lengths)
- **Access:** Via Zenodo mirror (https://zenodo.org/records/18837691)
- **Reference:** https://www.ncbi.nlm.nih.gov/bioproject/PRJNA737203

### 2. GENCODE v49 Genome Annotation
- **Source:** EMBL-EBI GENCODE Consortium
- **License:** CC BY 4.0 (https://www.gencodegenes.org/pages/license.html)
- **Format:** GFF3 format
- **Content:** Human genome annotation (release 49)
- **Access:** https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz

## Generation Pipeline (data-raw/scripts/preprocess.R)

The `preprocess.R` script performs the following steps to generate all inst/extdata files:

- Step 1: Read and Merge Salmon Quantification
- Step 2: Map Transcripts to Genes
- Step 2.5: Filter Single-Transcript Genes
- Step 3: Calculate Diversity Metrics.
- Step 4: Select Representative Genes
- Step 5: Generate Output Files

### Datasets

**metadata.tsv**:

- Source: PRJNA737203 sample annotations
- Content: Sample ID, condition (normal/tumor), paired sample information
- Format: Tab-separated values

**gene_ids.tsv**:

- Format: Plain text, one gene per line
- Content: 300 HGNC gene symbols
- Encoding: UTF-8

**tx2gene.tsv**:
- Format: Tab-separated values with header
- Columns: `transcript_id`, `gene_name`
- Rows: ~900 (transcript-gene pairs)
- Filter: Only includes transcripts in the 300 selected genes

**annotation.gff3.gz**:
- Format: GFF3 v3 (gzip compressed)
- Source: GENCODE v49 features
- Content: Gene-level and transcript-level annotations for selected 300 genes
- Structure: Standard GFF3 fields plus attributes (gene_id, gene_name, transcript_id, biotype, Parent)
- Filter: Extracting features only for transcripts in selected genes from full GENCODE

**gencode_subset_test.gff3.gz**:
- Format: GFF3 v3 (gzip compressed)
- Purpose: Small synthetic test dataset for unit testing
- Content: 3 genes with 9 transcripts total (12 GFF3 lines including parent genes)
- Structure: Follows GENCODE GFF3 format specification
- Note: Synthetic test data, not derived from preprocess.R output
