#!/bin/bash

################################################################################
# TSENAT Example Dataset: TCGA Luminal A Breast Cancer
# 
# Salmon Quantification Workflow
# Reference: data-raw/TOOLS_AND_PARAMETERS.txt
#
# This script documents the RNA-sequencing quantification pipeline used to 
# generate the example TSENAT dataset from TCGA BRCA Luminal A samples.
#
# NOTE: This is a REFERENCE/DOCUMENTATION script. For actual data generation,
# refer to GDC Data Portal (https://portal.gdc.cancer.gov/) for raw FASTQ files.
#
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directories
RAW_FASTQ_DIR="${1:-.}"  # Contains paired-end FASTQ files
TRIM_OUT_DIR="trimmed"
SALMON_INDEX_DIR="salmon_index"
SALMON_OUT_DIR="salmon_output"

# Tool parameters
SALMON_THREADS=12
SALMON_KMER=31  # For 75-150 bp reads; use 23 for shorter reads
SALMON_BOOTSTRAPS=30
SALMON_SEED=42
REFERENCE_TRANSCRIPTOME="gencode.v38.transcripts.fa"  # GRCh38
REFERENCE_GENOME="GRCh38.d1.vd1.fa"  # For decoy identification

# ============================================================================
# STEP 1: Raw Read Quality Assessment (FastQC)
# ============================================================================

run_fastqc() {
    echo "=== Step 1: FastQC Quality Assessment ==="
    
    if command -v fastqc &> /dev/null; then
        mkdir -p fastqc_results
        fastqc -t "${SALMON_THREADS}" \
               -o fastqc_results/ \
               "${RAW_FASTQ_DIR}"/*.fastq.gz
        
        # Aggregate with MultiQC
        if command -v multiqc &> /dev/null; then
            multiqc fastqc_results/ -o fastqc_results/
            echo "✓ FastQC and MultiQC complete. Review fastqc_results/multiqc_report.html"
        fi
    else
        echo "! FastQC not installed. Skipping QC step."
        echo "  Install via: conda install -c bioconda fastqc multiqc"
    fi
}

# ============================================================================
# STEP 2: Adapter Trimming & Quality Filtering (Trim Galore)
# ============================================================================

run_trimming() {
    echo "=== Step 2: Adapter Trimming (Trim Galore) ==="
    
    mkdir -p "${TRIM_OUT_DIR}"
    
    if command -v trim_galore &> /dev/null; then
        # Find sample names from R1 files
        for r1_file in "${RAW_FASTQ_DIR}"/*_R1*.fastq.gz; do
            sample_name=$(basename "${r1_file}" | sed 's/_R1.*//')
            r2_file="${RAW_FASTQ_DIR}/${sample_name}_R2.fastq.gz"
            
            echo "  Trimming: ${sample_name}"
            
            trim_galore --paired \
                        --gzip \
                        --cores "${SALMON_THREADS}" \
                        --output_dir "${TRIM_OUT_DIR}" \
                        "${r1_file}" "${r2_file}"
        done
        echo "✓ Trimming complete. Trimmed reads in ${TRIM_OUT_DIR}/"
    else
        echo "! Trim Galore not installed. Skipping trimming."
        echo "  Install via: conda install -c bioconda trim-galore"
        exit 1
    fi
}

# ============================================================================
# STEP 3: Salmon Index Construction (Decoy-Aware)
# ============================================================================

build_salmon_index() {
    echo "=== Step 3: Building Salmon Index (Decoy-Aware) ==="
    
    if ! command -v salmon &> /dev/null; then
        echo "! Salmon not installed."
        echo "  Install via: conda install -c bioconda salmon"
        exit 1
    fi
    
    if [ ! -f "${REFERENCE_TRANSCRIPTOME}" ]; then
        echo "! Reference transcriptome not found: ${REFERENCE_TRANSCRIPTOME}"
        echo "  Download GENCODE v38 (GRCh38) from: https://www.gencodegenes.org/"
        exit 1
    fi
    
    # Option 1: Standard index (simpler, faster to build)
    echo "  Building standard Salmon index..."
    salmon index -t "${REFERENCE_TRANSCRIPTOME}" \
                 -i "${SALMON_INDEX_DIR}" \
                 -k "${SALMON_KMER}" \
                 -p "${SALMON_THREADS}" \
                 --keepDuplicates
    
    # Option 2: Decoy-aware index (RECOMMENDED, reduces spurious mappings)
    # Uncomment if genome reference available:
    # echo "  Building decoy-aware index..."
    # grep "^>" "${REFERENCE_GENOME}" | sed 's/>//g' | sed 's/ .*//' > decoys.txt
    # salmon index -t "${REFERENCE_TRANSCRIPTOME}" \
    #              -d decoys.txt \
    #              -i "${SALMON_INDEX_DIR}_decoy" \
    #              -k "${SALMON_KMER}" \
    #              -p "${SALMON_THREADS}"
    
    echo "✓ Salmon index built: ${SALMON_INDEX_DIR}/"
}

# ============================================================================
# STEP 4: Salmon Quantification (Per-Sample)
# ============================================================================

run_salmon_quantification() {
    echo "=== Step 4: Salmon Quantification (Per-Sample) ==="
    
    if ! command -v salmon &> /dev/null; then
        echo "! Salmon not installed."
        exit 1
    fi
    
    mkdir -p "${SALMON_OUT_DIR}"
    local trimmed_ext="_val_1.fq.gz"  # Output suffix from Trim Galore
    
    for r1_file in "${TRIM_OUT_DIR}"/*${trimmed_ext}; do
        sample_name=$(basename "${r1_file}" | sed "s/${trimmed_ext}//")
        r2_file="${TRIM_OUT_DIR}/${sample_name}_val_2.fq.gz"
        
        echo "  Quantifying: ${sample_name}"
        
        salmon quant -i "${SALMON_INDEX_DIR}" \
                     -l A \
                     -1 "${r1_file}" -2 "${r2_file}" \
                     -p "${SALMON_THREADS}" \
                     --validateMappings \
                     --seqBias \
                     --gcBias \
                     --numBootstraps "${SALMON_BOOTSTRAPS}" \
                     --seed "${SALMON_SEED}" \
                     -o "${SALMON_OUT_DIR}/${sample_name}"
    done
    
    echo "✓ Salmon quantification complete."
    echo "  Output files in ${SALMON_OUT_DIR}/"
    echo "  Per-sample: ${SALMON_OUT_DIR}/*/quant.sf"
}

# ============================================================================
# STEP 5: Verify Quantification Results
# ============================================================================

verify_quantification() {
    echo "=== Step 5: Verifying Quantification Results ==="
    
    local failed=0
    
    for quant_file in "${SALMON_OUT_DIR}"/*/quant.sf; do
        if [ ! -f "${quant_file}" ]; then
            echo "! Missing: ${quant_file}"
            ((failed++))
        fi
    done
    
    if [ ${failed} -eq 0 ]; then
        echo "✓ All samples quantified successfully."
        echo ""
        echo "  Sample counts:"
        ls -d "${SALMON_OUT_DIR}"/*/ | wc -l
        echo ""
        echo "  Next steps:"
        echo "  1. Inspect ${SALMON_OUT_DIR}/*/meta_info.json for library type & mapping stats"
        echo "  2. Use R script: data-raw/import_and_aggregate.R (see below)"
    else
        echo "! ${failed} sample(s) failed quantification."
        exit 1
    fi
}

# ============================================================================
# STEP 6: Import & Aggregate with tximport (R Script)
# ============================================================================

create_tximport_script() {
    echo "=== Step 6: Prepare for tximport Aggregation ==="
    
    local r_script="import_and_aggregate.R"
    
    cat > "${r_script}" << 'EOF'
#!/usr/bin/env Rscript

# =========================================================================
# tximport: Aggregate Salmon transcript-level counts to gene-level
# =========================================================================

library(tximport)
library(SummarizedExperiment)

# Paths
salmon_out_dir <- "salmon_output"
tx2gene_file <- "inst/extdata/tx2gene.tsv"
coldata_file <- "inst/extdata/coldata.tsv"

# Load tx2gene mapping (transcript -> gene)
tx2gene <- read.csv(tx2gene_file, sep = "\t", header = TRUE)

# Load sample metadata
coldata <- read.csv(coldata_file, sep = "\t", header = TRUE)

# Locate all quant.sf files
samples <- list.dirs(salmon_out_dir, full.names = FALSE, recursive = FALSE)
files <- file.path(salmon_out_dir, samples, "quant.sf")
names(files) <- samples

# Verify all files exist
stopifnot(all(file.exists(files)))

# Import: TRANSCRIPT level counts
cat("Importing transcript-level counts...\n")
txi_tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)

# Save transcript-level
saveRDS(txi_tx, "txi_transcripts.RDS")
cat("✓ Saved: txi_transcripts.RDS\n")

# Summarize to GENE level (recommended for differential expression)
cat("Summarizing to gene level...\n")
txi_gene <- tximport(files, type = "salmon", 
                      tx2gene = tx2gene, 
                      countsFromAbundance = "lengthScaledTPM",
                      ignoreTxVersion = TRUE)

# Save gene-level
saveRDS(txi_gene, "txi_genes.RDS")
cat("✓ Saved: txi_genes.RDS\n")

# Export as data frames for TSENAT
cat("Exporting for TSENAT...\n")

# Transcript-level (what TSENAT uses)
readcounts_tx <- as.data.frame(txi_tx$counts)
readcounts_tx$Transcript <- rownames(readcounts_tx)
readcounts_tx <- readcounts_tx[, c(ncol(readcounts_tx), 1:(ncol(readcounts_tx)-1))]

# Gene-level (alternative)
readcounts_gene <- as.data.frame(txi_gene$counts)
readcounts_gene$Gene <- rownames(readcounts_gene)
readcounts_gene <- readcounts_gene[, c(ncol(readcounts_gene), 1:(ncol(readcounts_gene)-1))]

# Save as R data
save(readcounts_tx, file = "data/tcga_brca_luma.RData")
saveRDS(readcounts_tx, file = "data/tcga_brca_luma.RDS")

# Export TSV
write.csv(readcounts_tx, "data/tcga_brca_luma_dataset.tsv", row.names = FALSE, quote = FALSE)

cat("✓ Imports complete!\n")
cat("  Transcript-level: data/tcga_brca_luma.RData\n")
cat("  Gene-level: data/tcga_brca_luma_gene.RData (if needed)\n")
EOF

    chmod +x "${r_script}"
    echo "✓ Created: ${r_script}"
    echo "  Run with: Rscript ${r_script}"
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

main() {
    echo "================================================================================"
    echo "TSENAT TCGA Example Dataset: Salmon Quantification Workflow"
    echo "================================================================================"
    echo ""
    
    # Uncomment stages as needed (or pass as arguments):
    # run_fastqc
    # run_trimming
    # build_salmon_index
    # run_salmon_quantification
    # verify_quantification
    create_tximport_script
    
    echo ""
    echo "================================================================================"
    echo "Workflow documentation complete!"
    echo "================================================================================"
    echo ""
    echo "Reference documentation:"
    echo "  - Workflow details: data-raw/TOOLS_AND_PARAMETERS.txt"
    echo "  - Dataset info: inst/extdata/README.md"
    echo ""
}

# Run main function
main "$@"
