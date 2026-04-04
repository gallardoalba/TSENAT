#!/bin/bash

#########################################################################
# SALMON AAnalysis Pipeline
# 
# Description: Download, QC, trim, and quantify RNA-seq data from NCBI
# Input: SRA IDs from salmon_sample_ids.tsv
# Output: Salmon quantification results
#########################################################################

set -euo pipefail

# =============================================
# CONFIGURATION
# =============================================

# GENCODE reference files
ANNOTATION="gencode.v49.annotation.gff3.gz"
TRANSCRIPTOME="gencode.v49.transcripts.fa.gz"
GENOME="GRCh38.p14.genome.fa.gz"

GENCODE_RELEASE="49"
GENCODE_BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

# Thread configuration
THREADS=4
TMPDIR=-/tmp

# SRA settings
SRA_PREFETCH_RETRIES=3
SRA_PREFETCH_ATTEMPTS=0

# Directories
OUTPUT_DIR="../salmon_pipeline_output"
REF_DIR="${OUTPUT_DIR}/reference"
SALMON_INDEX_DIR="${OUTPUT_DIR}/salmon_index"
SALMON_QUANT_DIR="${OUTPUT_DIR}/salmon_quant"
FASTQ_DIR="${OUTPUT_DIR}/fastq"
QC_DIR="${OUTPUT_DIR}/qc"
TRIM_DIR="${OUTPUT_DIR}/trimmed"

mkdir -p "${REF_DIR}" "${SALMON_INDEX_DIR}" "${SALMON_QUANT_DIR}" "${FASTQ_DIR}" "${QC_DIR}" "${TRIM_DIR}"

# =============================================
# FUNCTIONS
# =============================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

error_exit() {
    log_message "ERROR: $1"
    exit 1
}

download_gencode_reference() {
    log_message "========================================="
    log_message "Downloading GENCODE Reference Files"
    log_message "========================================="
    
    # Check if all files already exist
    if [ -f "${REF_DIR}/${ANNOTATION}" ] && \
       [ -f "${REF_DIR}/${TRANSCRIPTOME}" ] && \
       [ -f "${REF_DIR}/${GENOME}" ]; then
        log_message "All GENCODE reference files already present, skipping download"
        # Update paths to reference directory
        ANNOTATION="${REF_DIR}/${ANNOTATION}"
        TRANSCRIPTOME="${REF_DIR}/${TRANSCRIPTOME}"
        GENOME="${REF_DIR}/${GENOME}"
        return 0
    fi
    
    # Download annotation
    if [ ! -f "${REF_DIR}/${ANNOTATION}" ]; then
        log_message "Downloading ${ANNOTATION}..."
        wget --no-check-certificate -q --show-progress \
            "${GENCODE_BASE_URL}/${ANNOTATION}" \
            -O "${REF_DIR}/${ANNOTATION}" || error_exit "Failed to download ${ANNOTATION}"
        log_message "✓ Downloaded ${ANNOTATION}"
    else
        log_message "${ANNOTATION} already exists"
    fi
    
    # Download transcriptome
    if [ ! -f "${REF_DIR}/${TRANSCRIPTOME}" ]; then
        log_message "Downloading ${TRANSCRIPTOME}..."
        wget --no-check-certificate -q --show-progress \
            "${GENCODE_BASE_URL}/${TRANSCRIPTOME}" \
            -O "${REF_DIR}/${TRANSCRIPTOME}" || error_exit "Failed to download ${TRANSCRIPTOME}"
        log_message "✓ Downloaded ${TRANSCRIPTOME}"
    else
        log_message "${TRANSCRIPTOME} already exists"
    fi
    
    # Download genome
    if [ ! -f "${REF_DIR}/${GENOME}" ]; then
        log_message "Downloading ${GENOME}..."
        wget --no-check-certificate -q --show-progress \
            "${GENCODE_BASE_URL}/${GENOME}" \
            -O "${REF_DIR}/${GENOME}" || error_exit "Failed to download ${GENOME}"
        log_message "✓ Downloaded ${GENOME}"
    else
        log_message "${GENOME} already exists"
    fi
    
    # Update paths to reference directory
    ANNOTATION="${REF_DIR}/${ANNOTATION}"
    TRANSCRIPTOME="${REF_DIR}/${TRANSCRIPTOME}"
    GENOME="${REF_DIR}/${GENOME}"
    
    log_message "All GENCODE reference files ready"
    log_message ""
}

download_sra() {
    local accession=$1
    log_message "Downloading SRA accession: ${accession}"
    
    SRA_PREFETCH_ATTEMPTS=1
    while [ $SRA_PREFETCH_ATTEMPTS -le $SRA_PREFETCH_RETRIES ]; do
        if fasterq-dump "${accession}" \
            -e "${THREADS}" \
            -t "${TMPDIR}" \
            --seq-defline '@$ac.$sn/$ri' \
            --qual-defline '+' \
            --split-3 \
            --skip-technical 2>&1; then
            
            fastq_count=$(ls -1 "${accession}"*.fastq 2>/dev/null | wc -l || echo 0)
            if [ "$fastq_count" -ge 1 ]; then
                log_message "Successfully downloaded ${accession} (${fastq_count} FASTQ files)"
                return 0
            fi
        fi
        
        log_message "Prefetch attempt $SRA_PREFETCH_ATTEMPTS of $SRA_PREFETCH_RETRIES failed for ${accession}, retrying..."
        SRA_PREFETCH_ATTEMPTS=$((SRA_PREFETCH_ATTEMPTS + 1))
        sleep 2
    done
    
    error_exit "Failed to download accession ${accession} after ${SRA_PREFETCH_RETRIES} attempts"
}

compress_fastqs() {
    local accession=$1
    local count=$(ls -1 "${accession}"*.fastq 2>/dev/null | wc -l || echo 0)
    
    log_message "Compressing FASTQ files for ${accession} (${count} files)"
    
    if [ "$count" -eq 1 ]; then
        # Single-end
        local fastq_file=$(ls "${accession}"*.fastq)
        pigz -cqp "${THREADS}" "${fastq_file}" > "${FASTQ_DIR}/${accession}__single.fastq.gz"
        rm "${fastq_file}"
        
    elif [ "$count" -eq 3 ]; then
        # Paired-end with technical reads
        if [ -f "${accession}.fastq" ]; then
            pigz -cqp "${THREADS}" "${accession}.fastq" > "${FASTQ_DIR}/${accession}__single.fastq.gz"
        fi
        pigz -cqp "${THREADS}" "${accession}_1.fastq" > "${FASTQ_DIR}/${accession}_R1.fastq.gz"
        pigz -cqp "${THREADS}" "${accession}_2.fastq" > "${FASTQ_DIR}/${accession}_R2.fastq.gz"
        rm "${accession}"*.fastq
        
    elif [ "$count" -eq 2 ]; then
        # Paired-end
        local fastq_array=($(ls "${accession}"*.fastq))
        pigz -cqp "${THREADS}" "${fastq_array[0]}" > "${FASTQ_DIR}/${accession}_R1.fastq.gz"
        pigz -cqp "${THREADS}" "${fastq_array[1]}" > "${FASTQ_DIR}/${accession}_R2.fastq.gz"
        rm "${accession}"*.fastq
    fi
}

run_falco() {
    local fastq_file=$1
    local accession=$(basename "${fastq_file}" .fastq.gz | sed 's/_R[12]__single$//')
    
    log_message "Running FALCO QC on ${fastq_file}"
    
    falco \
        --threads "${THREADS}" \
        --quiet \
        -f 'fastq.gz' \
        "${fastq_file}" \
        --skip-summary \
        -o "${QC_DIR}/${accession}_qc" 2>&1 || log_message "FALCO completed with status $?"
}

run_trim_galore() {
    local accession=$1
    local input1="${FASTQ_DIR}/${accession}_R1.fastq.gz"
    local input2="${FASTQ_DIR}/${accession}_R2.fastq.gz"
    
    log_message "Running Trim Galore on ${accession}"
    
    cd "${TRIM_DIR}"
    
    trim_galore \
        --cores "${THREADS}" \
        --phred33 \
        --quality 20 \
        --stringency 5 \
        -e 0.1 \
        --length 20 \
        --output_dir ./ \
        --no_report_file \
        --paired \
        "${input1}" "${input2}" 2>&1 || log_message "Trim Galore completed with status $?"
    
    # Standardize output names
    [ -f "${accession}_R1_val_1.fq.gz" ] && mv "${accession}_R1_val_1.fq.gz" "${accession}_R1_trimmed.fastq.gz"
    [ -f "${accession}_R2_val_2.fq.gz" ] && mv "${accession}_R2_val_2.fq.gz" "${accession}_R2_trimmed.fastq.gz"
    
    log_message "Trim Galore completed for ${accession}"
    cd - > /dev/null
}

build_salmon_index() {
    log_message "Building Salmon index"
    
    if [ ! -f "${TRANSCRIPTOME}" ]; then
        error_exit "Transcriptome file not found: ${TRANSCRIPTOME}"
    fi
    
    if [ ! -f "${GENOME}" ]; then
        error_exit "Genome file not found: ${GENOME}"
    fi
    
    # Create decoys file from genome
    zcat "${GENOME}" | grep "^>" | cut -d " " -f 1 | sed 's/>//' > "${SALMON_INDEX_DIR}/decoys.txt"
    
    # Combine transcriptome and genome
    zcat "${TRANSCRIPTOME}" "${GENOME}" > "${SALMON_INDEX_DIR}/index.fasta"
    
    # Build index
    salmon index \
        -i "${SALMON_INDEX_DIR}/salmon_idx" \
        --kmerLen 31 \
        --gencode \
        --threads "${THREADS}" \
        --transcripts "${SALMON_INDEX_DIR}/index.fasta" \
        --decoy "${SALMON_INDEX_DIR}/decoys.txt" \
        2>&1 || error_exit "Salmon index building failed"
    
    log_message "Salmon index built successfully"
}

prepare_salmon_outputs() {
    log_message "Preparing Salmon outputs for R preprocessing script"
    
    local salmon_merged_dir="${OUTPUT_DIR}/salmon_samples"
    mkdir -p "${salmon_merged_dir}"
    
    # Copy all quant.sf files to a single directory, renaming them by sample ID
    for quant_dir in "${SALMON_QUANT_DIR}"/*; do
        if [ -d "${quant_dir}" ]; then
            local accession=$(basename "${quant_dir}")
            local quant_file="${quant_dir}/quant.sf"
            
            if [ -f "${quant_file}" ]; then
                # Convert quant.sf to TSV format for R script
                # quant.sf is already TSV format, just rename and copy
                cp "${quant_file}" "${salmon_merged_dir}/${accession}.tsv"
                log_message "✓ Prepared ${accession}.tsv"
            fi
        fi
    done
    
    log_message "Salmon outputs prepared in: ${salmon_merged_dir}"
}

create_coldata_file() {
    log_message "Creating coldata.tsv file"
    
    local coldata_file="${OUTPUT_DIR}/salmon_coldata.tsv"
    local salmon_merged_dir="${OUTPUT_DIR}/salmon_samples"
    
    # Create header
    echo -e "sample_id\tcondition" > "${coldata_file}"
    
    # Extract conditions from SRA metadata if available
    # For now, assign conditions based on sample ordering (first half tumor, second half normal)
    local count=0
    local total=${#SRA_IDS[@]}
    local half=$((total / 2))
    
    for tsv_file in "${salmon_merged_dir}"/*.tsv; do
        if [ -f "${tsv_file}" ]; then
            local sample_id=$(basename "${tsv_file}" .tsv)
            local condition="tumor"
            
            # If second half of samples, assign "normal"
            if [ $count -ge $half ]; then
                condition="normal"
            fi
            
            echo -e "${sample_id}\t${condition}" >> "${coldata_file}"
            count=$((count + 1))
        fi
    done
    
    log_message "✓ Created coldata.tsv with $count samples"
}

run_preprocess_script() {
    log_message "========================================="
    log_message "Running R Preprocessing Script"
    log_message "========================================="
    
    local salmon_merged_dir="${OUTPUT_DIR}/salmon_samples"
    local coldata_file="${OUTPUT_DIR}/salmon_coldata.tsv"
    local output_counts="${OUTPUT_DIR}/salmon_counts_filtered.tsv"
    local output_gff="${OUTPUT_DIR}/annotation.gff3.gz"
    
    # Check if preprocess.R exists
    if [ ! -f "./scripts/preprocess.R" ]; then
        log_message "WARNING: preprocess_salmon_output.R not found, skipping R preprocessing"
        return 1
    fi
    
    log_message "Input salmon_samples directory: ${salmon_merged_dir}"
    log_message "Input coldata file: ${coldata_file}"
    log_message "Output counts file: ${output_counts}"
    log_message "Output GFF3 file: ${output_gff}"
    log_message ""
    
    # Run the R script with arguments
    if Rscript ./data-raw/scripts/preprocess.R \
        "${salmon_merged_dir}" \
        "${coldata_file}" \
        "${output_counts}" \
        "${output_gff}"; then
        log_message "✓ R preprocessing completed successfully"
        return 0
    else
        log_message "ERROR: R preprocessing script failed"
        return 1
    fi
}

run_salmon() {
    local accession=$1
    local r1="${TRIM_DIR}/${accession}_R1_trimmed.fastq.gz"
    local r2="${TRIM_DIR}/${accession}_R2_trimmed.fastq.gz"
    
    if [ ! -f "${r1}" ] || [ ! -f "${r2}" ]; then
        log_message "WARNING: Trimmed FASTQ files not found for ${accession}, skipping Salmon"
        return 1
    fi
    
    log_message "Running Salmon for ${accession}"
    
    local quant_dir="${SALMON_QUANT_DIR}/${accession}"
    mkdir -p "${quant_dir}"
    
    salmon quant \
        --index "${SALMON_INDEX_DIR}/salmon_idx" \
        --libType A \
        --mates1 <(zcat < "${r1}") \
        --mates2 <(zcat < "${r2}") \
        --threads "${THREADS}" \
        --validateMappings \
        --minScoreFraction 0.65 \
        --ma 2 \
        --mp -4 \
        --go 6 \
        --ge 2 \
        --consensusSlack 0.349999994 \
        --seqBias \
        --gcBias \
        --incompatPrior 0.0 \
        --minAssignedFrags 10 \
        --biasSpeedSamp 5 \
        --fldMax 1000 \
        --fldMean 250 \
        --fldSD 25 \
        --forgettingFactor 0.65 \
        --maxReadOcc 200 \
        --numBiasSamples 2000000 \
        --numAuxModelSamples 5000000 \
        --numPreAuxModelSamples 5000 \
        --rangeFactorizationBins 4 \
        --numGibbsSamples 0 \
        --numBootstraps 30 \
        --thinningFactor 16 \
        --sigDigits 3 \
        --vbPrior 0.01 \
        --output "${quant_dir}" 2>&1 || error_exit "Salmon quantification failed for ${accession}"
    
    log_message "Salmon completed for ${accession}"
}

# =============================================
# MAIN WORKFLOW
# =============================================

main() {
    log_message "=========================================="
    log_message "SALMON Processing Pipeline Started"
    log_message "=========================================="
    log_message "Threads: ${THREADS}"
    log_message "Temporary directory: ${TMPDIR}"
    log_message "Output directory: ${OUTPUT_DIR}"
    log_message ""
    
    # ========== STEP 0: Download Reference Files ==========
    download_gencode_reference
    
    # Read SRA IDs from file
    if [ ! -f "./inst/extdata/salmon_sample_ids.tsv" ]; then
        error_exit "SRA ID file not found: ./inst/extdata/salmon_sample_ids.tsv"
    fi
    
    mapfile -t SRA_IDS < ./inst/extdata/salmon_sample_ids.tsv
    log_message "Found ${#SRA_IDS[@]} SRA accessions to process"
    log_message ""
    
    # ========== STEP 1: Build Salmon Index =========="
    log_message "=========================================="
    log_message "STEP 1: Building Salmon Index"
    log_message "=========================================="
    build_salmon_index
    log_message ""
    
    # ========== STEP 2-5: Process each SRA accession ==========
    for i in "${!SRA_IDS[@]}"; do
        accession="${SRA_IDS[$i]}"
        sample_num=$((i + 1))
        
        log_message "=========================================="
        log_message "Processing Sample ${sample_num}/${#SRA_IDS[@]}: ${accession}"
        log_message "=========================================="
        
        # Download
        if ! download_sra "${accession}"; then
            log_message "WARNING: Failed to download ${accession}, continuing to next sample"
            continue
        fi
        
        # Compress
        compress_fastqs "${accession}"
        
        # QC - FALCO
        log_message "Running FALCO quality control"
        for fastq in "${FASTQ_DIR}/${accession}"*.fastq.gz; do
            [ -f "${fastq}" ] && run_falco "${fastq}"
        done
        
        # Trim
        if [ -f "${FASTQ_DIR}/${accession}_R1.fastq.gz" ] && [ -f "${FASTQ_DIR}/${accession}_R2.fastq.gz" ]; then
            run_trim_galore "${accession}"
        fi
        
        # Quantify with Salmon
        if run_salmon "${accession}"; then
            log_message "✓ Completed all steps for ${accession}"
        else
            log_message "⚠ Warning: Some steps incomplete for ${accession}"
        fi
        
        log_message ""
    done
    
    # ========== STEP 6: Prepare and Process Salmon Outputs with R ==========
    log_message "========================================="
    log_message "STEP 6: Post-processing with R Script"
    log_message "========================================="
    
    prepare_salmon_outputs
    create_coldata_file
    run_preprocess_script
    
    log_message ""
    log_message "========================================="
    log_message "Galaxy SRA Processing Pipeline Completed"
    log_message "========================================="
    log_message "Results saved to: ${OUTPUT_DIR}"
    log_message "Salmon quantifications in: ${SALMON_QUANT_DIR}"
    log_message "Final outputs:"
    log_message "  - Filtered counts: ${OUTPUT_DIR}/salmon_counts_filtered.tsv"
    log_message "  - GFF3 annotation: ${OUTPUT_DIR}/annotation.gff3.gz"
    log_message "  - RData object: ${OUTPUT_DIR}/salmon_counts_filtered.RData"
    log_message "  - Tx2gene mapping: ${OUTPUT_DIR}/tx2gene.tsv"
}

# Run main workflow
main "$@"
