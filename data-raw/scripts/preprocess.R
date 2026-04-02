#!/usr/bin/env Rscript
#' Bioconductor Dataset Preprocessing Pipeline
#' 
#' Processes Salmon quantification output to generate filtered GFF3 and counts
#' for Bioconductor submission. Includes gene selection, validation, and test data.
#'
#' Usage: Rscript preprocess_v4.R [salmon_dir] [metadata_file] [output_counts] [output_gff]

suppressMessages(suppressWarnings({
    library("SplicingFactory")
    library("SummarizedExperiment")
    library("tidyverse")
    library("BiocParallel")
    library("data.table")
    library("matrixStats")
}))

set.seed(69)
start_time <- Sys.time()

# =============================================================================
# CONFIGURATION
# =============================================================================

parse_config <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    list(
        salmon_dir = ifelse(length(args) > 0, args[1], "/tmp/salmon"),
        metadata_file = ifelse(length(args) > 1, args[2], "./inst/extdata/metadata.tsv"),
        output_gff = ifelse(length(args) > 2, args[3], "./inst/extdata/annotation.gff3.gz"),
        gencode_file = "/tmp/gencode.v49.annotation.gff3.gz",
        extdata_dir = "./inst/extdata"
    )
}

detect_condition_columns <- function(metadata) {
    condition_col <- if ("condition" %in% colnames(metadata)) {
        "condition"
    } else if ("phenotype" %in% tolower(colnames(metadata))) {
        colnames(metadata)[grep("phenotype", tolower(colnames(metadata)))[1]]
    } else {
        colnames(metadata)[2]
    }
    
    paired_col <- if ("paired" %in% tolower(colnames(metadata))) {
        colnames(metadata)[grep("paired", tolower(colnames(metadata)))[1]]
    } else if ("pair_id" %in% tolower(colnames(metadata))) {
        colnames(metadata)[grep("pair_id", tolower(colnames(metadata)))[1]]
    } else if (ncol(metadata) > 2) {
        colnames(metadata)[3]
    } else {
        NULL
    }
    
    list(condition_col = condition_col, paired_col = paired_col)
}

print_config <- function(config) {
    message("\n=== BIOCONDUCTOR DATASET PREPROCESSING ===")
    message("Salmon samples:  ", config$salmon_dir)
    message("Metadata:        ", config$metadata_file)
    message("Output GFF3:     ", config$output_gff)
    message("Output RData:    ./data/readcounts.RData")
}

check_directories <- function(config) {
    for (dir in c("./data", "./output", config$extdata_dir)) {
        if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    }
}

# =============================================================================
# DOWNLOAD METADATA FROM ZENODO (if needed)
# =============================================================================

download_metadata_from_zenodo <- function(config) {
    message("\n[METADATA] Checking metadata availability...")
    
    # Check if metadata already exists
    if (file.exists(config$metadata_file)) {
        message("✓ Metadata file found: ", config$metadata_file)
        return(TRUE)
    }
    
    message("Metadata file not found. Downloading from Zenodo...")
    
    # Zenodo record URL
    zenodo_url <- "https://zenodo.org/records/18837691/files/metadata.tsv"
    
    # Create extdata directory if needed
    if (!dir.exists(config$extdata_dir)) {
        dir.create(config$extdata_dir, recursive = TRUE)
        message("Created directory: ", config$extdata_dir)
    }
    
    # Download file
    tryCatch(
        {
            message("Downloading from: ", zenodo_url)
            download.file(zenodo_url, config$metadata_file, mode = "wb", quiet = FALSE)
            
            if (file.exists(config$metadata_file)) {
                file_size <- file.size(config$metadata_file) / 1024  # Size in KB
                message("✓ Successfully downloaded metadata: ", config$metadata_file)
                message("  File size: ", round(file_size, 2), " KB")
                return(TRUE)
            } else {
                stop("File not found after download attempt")
            }
        },
        error = function(e) {
            message("\n⚠ ERROR: Could not download metadata from Zenodo")
            message("  URL: ", zenodo_url)
            message("  Error: ", e$message)
            message("\n  Please manually download and place at: ", config$metadata_file)
            message("  Or provide metadata file as command line argument")
            return(FALSE)
        }
    )
}

# =============================================================================
# DOWNLOAD GENCODE FROM EBI (if needed)
# =============================================================================

download_gencode_from_ebi <- function(config) {
    message("\n[GENCODE] Checking gencode annotation file...")
    
    # Check if gencode already exists
    if (file.exists(config$gencode_file)) {
        file_size <- file.size(config$gencode_file) / (1024^3)  # Size in GB
        message("✓ Gencode file found: ", config$gencode_file)
        message("  File size: ", round(file_size, 2), " GB")
        return(TRUE)
    }
    
    message("Gencode file not found. Downloading from EBI FTP...")
    
    # EBI FTP URL for gencode v49
    ebi_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz"
    
    # Ensure /tmp directory exists (it should, but just in case)
    tmp_dir <- dirname(config$gencode_file)
    if (!dir.exists(tmp_dir)) {
        dir.create(tmp_dir, recursive = TRUE)
        message("Created directory: ", tmp_dir)
    }
    
    # Download file with progress
    tryCatch(
        {
            message("Downloading from: ", ebi_url)
            message("This may take several minutes (~1-2 GB file)...")
            
            # Use curl if available for better progress reporting
            download_method <- if (nzchar(Sys.which("curl"))) "curl" else "auto"
            
            download.file(
                ebi_url, 
                config$gencode_file, 
                mode = "wb",
                method = download_method,
                quiet = FALSE
            )
            
            if (file.exists(config$gencode_file)) {
                file_size <- file.size(config$gencode_file) / (1024^3)  # Size in GB
                message("✓ Successfully downloaded gencode: ", config$gencode_file)
                message("  File size: ", round(file_size, 2), " GB")
                return(TRUE)
            } else {
                stop("File not found after download attempt")
            }
        },
        error = function(e) {
            message("\n⚠ ERROR: Could not download gencode from EBI FTP")
            message("  URL: ", ebi_url)
            message("  Error: ", e$message)
            message("\n  Please manually download and place at: ", config$gencode_file)
            message("  Or download with: wget ", ebi_url)
            return(FALSE)
        }
    )
}

# =============================================================================
# DOWNLOAD SALMON SAMPLES FROM ZENODO (if needed)
# =============================================================================

download_salmon_samples_from_zenodo <- function(config) {
    message("\n[SALMON] Checking salmon sample files...")
    
    salmon_dir <- "/tmp/salmon"
    
    # Check if salmon samples already exist
    existing_files <- list.files(salmon_dir, pattern = "\\.(tsv|tsv\\.gz)$", full.names = FALSE)
    if (length(existing_files) > 0) {
        message("✓ Salmon samples found in: ", salmon_dir)
        message("  Files: ", length(existing_files), " samples")
        return(salmon_dir)
    }
    
    message("Salmon samples not found. Downloading from Zenodo...")
    
    # Zenodo record ID
    zenodo_record_id <- "18837691"
    zenodo_api_url <- paste0("https://zenodo.org/api/records/", zenodo_record_id)
    
    # Create /tmp/salmon directory if it doesn't exist
    if (!dir.exists(salmon_dir)) {
        dir.create(salmon_dir, recursive = TRUE)
        message("Created directory: ", salmon_dir)
    }
    
    tryCatch({
        # Load JSON parsing library
        if (!require("jsonlite", quietly = TRUE)) {
            message("Installing jsonlite for API parsing...")
            install.packages("jsonlite", repos = "https://cloud.r-project.org", quiet = TRUE)
            library("jsonlite", quietly = TRUE)
        }
        
        # Fetch and parse Zenodo API response
        message("Fetching file list from Zenodo API...")
        api_response <- tryCatch(
            jsonlite::fromJSON(zenodo_api_url),
            error = function(e) {
                # Fallback: try with readLines and manual parsing
                json_text <- paste(readLines(zenodo_api_url, warn = FALSE), collapse = "")
                jsonlite::fromJSON(json_text)
            }
        )
        
        # Extract file information
        if (!is.null(api_response$files) && length(api_response$files) > 0) {
            files_df <- api_response$files
            tsv_files <- files_df[grep("\\.(tsv|tsv\\.gz)$", files_df$key, ignore.case = TRUE), ]
            
            if (nrow(tsv_files) == 0) {
                stop("No TSV files found in Zenodo record")
            }
            
            message("Found ", nrow(tsv_files), " TSV files to download")
            
            # Download each file
            download_count <- 0
            for (i in seq_len(nrow(tsv_files))) {
                filename <- tsv_files$key[i]
                download_url <- tsv_files$links.self[i]
                
                # Fallback URL construction if links.self not available
                if (is.na(download_url) || is.null(download_url)) {
                    download_url <- paste0("https://zenodo.org/records/", zenodo_record_id, "/files/", filename)
                }
                
                output_path <- file.path(salmon_dir, basename(filename))
                
                message("Downloading: ", filename, " ...")
                tryCatch({
                    download.file(download_url, output_path, mode = "wb", quiet = TRUE, timeout = 300)
                    if (file.exists(output_path) && file.size(output_path) > 0) {
                        message("  ✓ Downloaded: ", filename, " (", round(file.size(output_path)/1024/1024, 1), " MB)")
                        download_count <- download_count + 1
                    }
                }, error = function(e) {
                    message("  ✗ Failed to download: ", filename, " - ", e$message)
                })
            }
            
            # Verify files were downloaded
            downloaded_files <- list.files(salmon_dir, pattern = "\\.(tsv|tsv\\.gz)$")
            if (length(downloaded_files) > 0) {
                message("✓ Downloaded ", length(downloaded_files), " salmon samples successfully")
                return(salmon_dir)
            } else {
                stop("No files were successfully downloaded from Zenodo")
            }
        } else {
            stop("Invalid Zenodo API response: no files found")
        }
        
    }, error = function(e) {
        message("\n⚠ WARNING: Could not automatically download from Zenodo")
        message("  Error: ", e$message)
        message("\nManual download instructions:")
        message("  1. Visit: https://zenodo.org/records/18837691")
        message("  2. Click 'Download all' to get the archive")
        message("  3. Extract TSV files and place in: /tmp/salmon")
        message("  4. Re-run this script")
        message("\nOr download individual files:")
        message("  wget https://zenodo.org/records/18837691/files/<filename>")
        return(config$salmon_dir)
    })
}

# =============================================================================
# STEP 1: READ AND MERGE SALMON SAMPLES
# =============================================================================

read_salmon_samples <- function(salmon_dir, metadata_file) {
    message("\n[STEP 1] Reading and merging Salmon samples...")
    step_start <- Sys.time()
    
    # Load metadata
    metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
    rownames(metadata) <- metadata[, 1]
    
    # Detect condition and paired columns for later use
    attr(metadata, "condition_col") <- detect_condition_columns(metadata)$condition_col
    attr(metadata, "paired_col") <- detect_condition_columns(metadata)$paired_col
    
    # Find Salmon files (supports both .tsv and .tsv.gz)
    quant_files <- c(
        list.files(salmon_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE),
        list.files(salmon_dir, pattern = "\\.tsv$", full.names = TRUE)
    )
    if (length(quant_files) == 0) stop("No .tsv or .tsv.gz files found in ", salmon_dir)
    
    # Extract sample names (handle both .tsv and .tsv.gz extensions)
    sample_names <- sub("\\.tsv(\\.gz)?$", "", basename(quant_files))
    message("Found ", length(quant_files), " samples")
    
    # Read files in parallel
    quant_list <- bplapply(quant_files, function(qf) {
        suppressWarnings(suppressMessages(library("data.table")))
        fread(qf, data.table = FALSE, verbose = FALSE)
    })
    names(quant_list) <- sample_names
    
    # Build matrices
    all_transcripts <- unique(unlist(lapply(quant_list, function(x) x$Name)))
    counts <- build_matrix(quant_list, all_transcripts, sample_names, "NumReads")
    tpm <- build_matrix(quant_list, all_transcripts, sample_names, "TPM")
    eff_length <- build_length_vector(quant_list, all_transcripts)
    
    message("✓ Loaded SALMON data:")
    message("  - Counts (NumReads)")
    message("  - TPM (normalized expression)")
    message("  - EffectiveLength (length bias corrected)")
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 1 completed in ", round(elapsed, 2), " minutes")
    message("  Transcripts: ", length(all_transcripts), ", Samples: ", ncol(counts))
    
    list(counts = counts, tpm = tpm, effective_length = eff_length, 
         all_transcripts = all_transcripts, sample_names = sample_names,
         metadata = metadata)
}

build_matrix <- function(quant_list, all_transcripts, sample_names, column_name) {
    # Optimized: Use data.table joins (faster than indexing)
    dt <- data.table(transcript_id = all_transcripts, key = "transcript_id")
    
    for (i in seq_along(quant_list)) {
        sample_dt <- data.table(
            Name = quant_list[[i]]$Name,
            value = quant_list[[i]][[column_name]],
            key = "Name"
        )
        dt[sample_dt, (sample_names[i]) := i.value]
    }
    
    # Keep transcript IDs as rownames before converting to matrix
    tx_ids <- dt$transcript_id
    dt_matrix <- as.matrix(dt[, -"transcript_id"])
    rownames(dt_matrix) <- tx_ids
    dt_matrix
}

build_length_vector <- function(quant_list, all_transcripts) {
    # Optimized: Use rowMedians for median effective length across samples
    lengths_list <- lapply(quant_list, function(df) {
        setNames(df$EffectiveLength, df$Name)
    })
    lengths_matrix <- rowMedians(
        as.matrix(sapply(lengths_list, function(x) {
            x[all_transcripts]  # Align to all transcripts
        })),
        na.rm = TRUE
    )
    names(lengths_matrix) <- all_transcripts
    lengths_matrix
}

# =============================================================================
# STEP 2: MAP TRANSCRIPTS TO GENES
# =============================================================================

map_transcripts_to_genes <- function(readcounts, gencode_file) {
    message("\n[STEP 2] Mapping transcripts to genes...")
    step_start <- Sys.time()
    
    if (!file.exists(gencode_file)) {
        stop("Gencode file not found: ", gencode_file)
    }
    
    gencode_mapping <- extract_gencode_mapping(gencode_file)
    tx_ids <- rownames(readcounts)
    tx_map <- match_transcripts(tx_ids, gencode_mapping)
    
    mapped <- sum(!is.na(tx_map$gene_name))
    message("Mapped ", mapped, "/", length(tx_ids), " transcripts to genes")
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 2 completed in ", round(elapsed, 2), " minutes")
    
    list(tx2gene = tx_map, mapping = gencode_mapping)
}

extract_gencode_mapping <- function(gencode_file) {
    message("  Building gencode index (optimized with shell extraction)...")
    
    # Use shell pipeline to extract transcript features (much faster than R parsing)
    # Uses grep with tab literal match like preprocess_v4.R
    shell_cmd <- paste0(
        "zcat ", shQuote(gencode_file), 
        " | grep '\ttranscript\t' | cut -f9"
    )
    
    attr_lines <- system(shell_cmd, intern = TRUE)
    n_transcripts_found <- length(attr_lines)
    message("    Found ", n_transcripts_found, " transcript features")
    
    if (n_transcripts_found == 0) {
        warning("No transcript features found in gencode file. Check file format and path.")
        return(data.frame(
            Transcript = character(0),
            Gene = character(0),
            GeneName = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Parse attributes using regex (faster than line-by-line parsing)
    
    # Extract transcript_id (gencode format: transcript_id=VALUE;...)
    transcript_ids <- sub('.*transcript_id=([^;]+);.*', '\\1', attr_lines)
    transcript_ids[!grepl('^ENST', transcript_ids)] <- ""
    
    # Extract gene_id
    gene_ids <- sub('.*gene_id=([^;]+);.*', '\\1', attr_lines)
    gene_ids[!grepl('^ENSG', gene_ids)] <- ""
    
    # Extract gene_name
    gene_names <- sub('.*gene_name=([^;]+);.*', '\\1', attr_lines)
    gene_names[gene_names == attr_lines] <- ""  # Reset failed matches
    
    # Filter to valid rows
    valid_mask <- (transcript_ids != "") & (gene_ids != "") & (gene_names != "")
    n_valid <- sum(valid_mask)
    message("    Valid transcripts: ", n_valid, "/", length(valid_mask))
    
    if (n_valid == 0) {
        warning("No valid transcript mappings extracted from gencode attributes")
    }
    
    data.frame(
        Transcript = transcript_ids[valid_mask],
        Gene = gene_ids[valid_mask],
        GeneName = gene_names[valid_mask],
        stringsAsFactors = FALSE
    )
}

extract_gff_attr <- function(attributes, attr_name) {
    pattern <- paste0(attr_name, "=")
    pos <- regexpr(pattern, attributes, fixed = TRUE) + nchar(pattern)
    if (pos <= nchar(pattern)) return(NA_character_)
    
    rest <- substr(attributes, pos, nchar(attributes))
    end_pos <- regexpr(";", rest, fixed = TRUE)[1]
    if (end_pos < 0) end_pos <- nchar(rest) + 1
    
    substr(rest, 1, end_pos - 1)
}

match_transcripts <- function(tx_ids, gencode_mapping) {
    gene_names <- rep(NA_character_, length(tx_ids))
    
    # Try exact match first (full versioned ID)
    exact_match <- match(tx_ids, gencode_mapping$Transcript)
    matched_exact <- sum(!is.na(exact_match))
    gene_names[!is.na(exact_match)] <- gencode_mapping$GeneName[exact_match[!is.na(exact_match)]]
    
    # For unmatched, try version-free matching (faster vectorized approach)
    unmatched_idx <- which(is.na(gene_names))
    if (length(unmatched_idx) > 0) {
        tx_clean <- sub("\\..*", "", tx_ids[unmatched_idx])
        gen_clean <- sub("\\..*", "", gencode_mapping$Transcript)
        version_matches <- match(tx_clean, gen_clean)
        gene_names[unmatched_idx[!is.na(version_matches)]] <- 
            gencode_mapping$GeneName[version_matches[!is.na(version_matches)]]
    }
    
    data.frame(transcript_id = tx_ids, gene_name = gene_names, 
               stringsAsFactors = FALSE)
}

# =============================================================================
# STEP 3: CALCULATE DIVERSITY METRICS
# =============================================================================

calculate_diversity_metrics <- function(readcounts, tx2gene, metadata) {
    message("\n[STEP 3] Calculating diversity metrics (using SplicingFactory)...")
    step_start <- Sys.time()
    
    # Get gene assignment sorted by gene name (required for SplicingFactory)
    gene_assignment <- tx2gene$gene_name
    names(gene_assignment) <- tx2gene$transcript_id
    
    # Sort to match readcounts order
    gene_sorted <- gene_assignment[rownames(readcounts)]
    count_matrix <- as.matrix(readcounts)
    count_cols <- colnames(readcounts)
    
    message("  Computing Laplace diversity (vectorized with SplicingFactory)...")
    Laplace_diversity <- calculate_diversity(count_matrix, gene_sorted, method = "laplace")
    
    # Get the condition and paired columns from metadata attributes
    coldata <- metadata
    condition_col <- attr(coldata, "condition_col")
    if (is.null(condition_col)) condition_col <- "condition"
    
    paired_col <- attr(coldata, "paired_col")
    
    # Extract condition information - map sample IDs to conditions
    sample_conditions <- as.character(coldata[colnames(Laplace_diversity), condition_col])
    
    # Detect control condition: prefer "normal" or "control", fall back to first condition
    unique_conditions <- unique(tolower(sample_conditions))
    control_condition <- intersect(c("normal", "control"), unique_conditions)
    if (length(control_condition) == 0) {
        control_condition <- unique_conditions[1]
    }
    
    # Build per-transcript results dataframe (following v4.R pattern exactly)
    transcript_ids <- rownames(count_matrix)
    n_transcripts <- length(transcript_ids)
    
    if (n_transcripts == 0) {
        warning("No transcripts found in count matrix")
        return(NULL)
    }
    
    Laplace_readcount_Wilcox_df <- data.frame(
        genes = gene_sorted,
        transcript_id = transcript_ids,
        mean = rowMeans(count_matrix[, count_cols, drop = FALSE], na.rm = TRUE),
        sd = apply(count_matrix[, count_cols], 1, sd, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    
    # Calculate coefficient of variation
    Laplace_readcount_Wilcox_df$cv <- with(Laplace_readcount_Wilcox_df, {
        ifelse(mean > 0, sd / mean, 0)
    })
    
    # Calculate fold-change between conditions
    normal_samples <- rownames(coldata)[which(tolower(coldata[[condition_col]]) == control_condition)]
    tumor_samples <- rownames(coldata)[which(tolower(coldata[[condition_col]]) != control_condition)]
    
    normal_cols <- intersect(normal_samples, count_cols)
    tumor_cols <- intersect(tumor_samples, count_cols)
    
    if (length(normal_cols) > 0 && length(tumor_cols) > 0) {
        normal_mean <- rowMeans(count_matrix[, normal_cols, drop = FALSE], na.rm = TRUE)
        tumor_mean <- rowMeans(count_matrix[, tumor_cols, drop = FALSE], na.rm = TRUE)
        
        # Add pseudocount to avoid division by zero
        normal_mean <- normal_mean + 0.1
        tumor_mean <- tumor_mean + 0.1
        
        Laplace_readcount_Wilcox_df$log2_fold_change <- log2(tumor_mean / normal_mean)
    } else {
        Laplace_readcount_Wilcox_df$log2_fold_change <- 0
    }
    
    # ========================================================================
    # PERFORM WILCOXON TEST ON ENTROPY VALUES (following v4.R exactly)
    # ========================================================================
    message("  Performing Wilcoxon test on entropy values between conditions...")
    
    # Extract entropy values from Laplace_diversity object
    entropy_matrix <- tryCatch(
        assay(Laplace_diversity, "diversity"),
        error = function(e) {
            tryCatch(
                assay(Laplace_diversity, "entropy"),
                error = function(e2) {
                    assay(Laplace_diversity)
                }
            )
        }
    )
    
    # Get entropy for each transcript
    entropy_data <- data.frame(
        transcript_id = rownames(entropy_matrix),
        entropy_matrix,
        stringsAsFactors = FALSE
    )
    
    # Initialize p-values vector
    message("  Computing p-values with Wilcoxon test...")
    p_values <- numeric(nrow(entropy_data))
    
    # Check if we have paired samples (following v4.R pattern)
    has_paired <- !is.null(paired_col) && all(!is.na(coldata[colnames(Laplace_diversity), paired_col]))
    
    if (has_paired && length(normal_cols) > 1 && length(tumor_cols) > 1) {
        message("  Paired samples detected. Matching samples by pair ID from metadata column: ", paired_col)
        
        # Extract pair information from metadata for samples in our analysis
        pair_info <- data.frame(
            sample_id = colnames(Laplace_diversity),
            condition = coldata[colnames(Laplace_diversity), condition_col],
            pair_id = coldata[colnames(Laplace_diversity), paired_col],
            stringsAsFactors = FALSE
        )
        
        # Get unique pair IDs and verify pairing structure
        unique_pairs <- unique(pair_info$pair_id)
        message("  Found ", length(unique_pairs), " unique pairs")
        
        # Validate that each pair has exactly one normal and one tumor sample
        valid_pairs <- list()
        valid_pair_count <- 0
        
        for (pair in unique_pairs) {
            pair_samples <- pair_info[pair_info$pair_id == pair, ]
            normal_in_pair <- pair_samples$sample_id[tolower(pair_samples$condition) == control_condition]
            tumor_in_pair <- pair_samples$sample_id[tolower(pair_samples$condition) != control_condition]
            
            # Check if this pair has exactly one normal and one tumor
            if (length(normal_in_pair) == 1 && length(tumor_in_pair) >= 1) {
                # For cases where multiple tumor samples exist for one normal,
                # we use the first tumor sample
                valid_pairs[[as.character(pair)]] <- list(
                    normal = normal_in_pair,
                    tumor = tumor_in_pair[1]
                )
                valid_pair_count <- valid_pair_count + 1
            }
        }
        
        message("  Valid paired samples: ", valid_pair_count, " pairs")
        
        # Extract matched datasets for paired test
        paired_normal_samples <- sapply(valid_pairs, function(x) x$normal)
        paired_tumor_samples <- sapply(valid_pairs, function(x) x$tumor)
        
        # Ensure samples exist in entropy matrix
        paired_normal_samples <- paired_normal_samples[paired_normal_samples %in% colnames(entropy_matrix)]
        paired_tumor_samples <- paired_tumor_samples[paired_tumor_samples %in% colnames(entropy_matrix)]
        
        if (length(paired_normal_samples) > 1 && length(paired_tumor_samples) > 1) {
            message("  Performing paired Wilcoxon test with ", length(paired_normal_samples), " matched pairs")
            
            # Fast serial approach - no parallel overhead
            for (i in seq_len(nrow(entropy_data))) {
                normal_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], paired_normal_samples])
                tumor_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], paired_tumor_samples])
                
                if (sum(!is.na(normal_entropy)) > 0 && sum(!is.na(tumor_entropy)) > 0) {
                    tryCatch(
                        {
                            p_values[i] <- wilcox.test(normal_entropy, tumor_entropy, paired = TRUE, exact = FALSE)$p.value
                        },
                        error = function(e) NULL,
                        warning = function(w) NULL
                    )
                }
            }
            message("  ✓ Paired Wilcoxon test completed")
        } else {
            # Not enough valid pairs - fall back to unpaired test
            message("  WARNING: Not enough valid pairs (need at least 2). Falling back to unpaired Wilcoxon test...")
            
            normal_cols_available <- intersect(normal_cols, colnames(entropy_matrix))
            tumor_cols_available <- intersect(tumor_cols, colnames(entropy_matrix))
            
            if (length(normal_cols_available) > 0 && length(tumor_cols_available) > 0) {
                for (i in seq_len(nrow(entropy_data))) {
                    normal_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], normal_cols_available])
                    tumor_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], tumor_cols_available])
                    
                    if (sum(!is.na(normal_entropy)) > 0 && sum(!is.na(tumor_entropy)) > 0) {
                        tryCatch(
                            {
                                p_values[i] <- wilcox.test(normal_entropy, tumor_entropy, paired = FALSE, exact = FALSE)$p.value
                            },
                            error = function(e) NULL,
                            warning = function(w) NULL
                        )
                    }
                }
                message("  ✓ Unpaired Wilcoxon test completed")
            }
        }
    } else if (length(normal_cols) > 0 && length(tumor_cols) > 0) {
        # No paired information available - use unpaired test
        message("  No paired sample information found. Performing unpaired Wilcoxon test...")
        
        normal_cols_available <- intersect(normal_cols, colnames(entropy_matrix))
        tumor_cols_available <- intersect(tumor_cols, colnames(entropy_matrix))
        
        if (length(normal_cols_available) > 0 && length(tumor_cols_available) > 0) {
            for (i in seq_len(nrow(entropy_data))) {
                normal_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], normal_cols_available])
                tumor_entropy <- as.numeric(entropy_matrix[entropy_data$transcript_id[i], tumor_cols_available])
                
                if (sum(!is.na(normal_entropy)) > 0 && sum(!is.na(tumor_entropy)) > 0) {
                    tryCatch(
                        {
                            p_values[i] <- wilcox.test(normal_entropy, tumor_entropy, paired = FALSE, exact = FALSE)$p.value
                        },
                        error = function(e) NULL,
                        warning = function(w) NULL
                    )
                }
            }
            message("  ✓ Unpaired Wilcoxon test completed")
        }
    }
    
    message("  Wilcoxon test completed successfully")
    
    # Calculate adjusted p-values using Benjamini-Hochberg method
    adjusted_p_values <- p.adjust(p_values, method = "BH")
    
    # Add p-values to the results dataframe (matching v4.R naming: adj_p_value)
    Laplace_readcount_Wilcox_df$p_value <- p_values[match(Laplace_readcount_Wilcox_df$transcript_id, entropy_data$transcript_id)]
    Laplace_readcount_Wilcox_df$adj_p_value <- adjusted_p_values[match(Laplace_readcount_Wilcox_df$transcript_id, entropy_data$transcript_id)]
    
    # Sort by adjusted p-value (most significant first)
    Laplace_readcount_Wilcox_df <- Laplace_readcount_Wilcox_df[order(Laplace_readcount_Wilcox_df$adj_p_value, na.last = TRUE), ]
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 3 completed in ", round(elapsed, 2), " minutes")
    message("  Transcripts with entropy and p-values: ", nrow(Laplace_readcount_Wilcox_df))
    
    return(Laplace_readcount_Wilcox_df)
}

# =============================================================================
# STEP 2.5: FILTER SINGLE-TRANSCRIPT GENES
# =============================================================================

filter_single_transcript_genes <- function(readcounts, tx2gene) {
    message("\n[STEP 2.5] Filtering genes with single transcript...")
    step_start <- Sys.time()
    
    # Count transcripts per gene
    transcripts_per_gene <- table(tx2gene$gene_name)
    genes_with_multi_tx <- names(transcripts_per_gene[transcripts_per_gene > 1])
    
    # Filter tx2gene and readcounts
    multi_tx_mask <- tx2gene$gene_name %in% genes_with_multi_tx
    tx2gene_filtered <- tx2gene[multi_tx_mask, , drop = FALSE]
    readcounts_filtered <- readcounts[tx2gene_filtered$transcript_id, , drop = FALSE]
    
    removed <- sum(!multi_tx_mask)
    message("  Removed ", removed, " transcripts from single-transcript genes")
    message("  Kept ", nrow(readcounts_filtered), " transcripts in ", length(genes_with_multi_tx), " multi-transcript genes")
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 2.5 completed in ", round(elapsed, 2), " minutes")
    
    list(readcounts = readcounts_filtered, tx2gene = tx2gene_filtered)
}

# =============================================================================
# STEP 4: SELECT TOP GENES
# =============================================================================

select_top_genes <- function(readcounts, tx2gene, diversity_data, sample_size = NULL) {
    message("\n[STEP 4] Selecting top genes...")
    step_start <- Sys.time()
    
    if (nrow(diversity_data) == 0) {
        warning("No genes available for selection")
        return(character(0))
    }
    
    # Group diversity_data by gene and calculate statistics (data.table approach like v4.R)
    diversity_dt <- as.data.table(diversity_data)
    
    gene_ranking <- diversity_dt[, .(
        min_adj_pvalue = min(adj_p_value, na.rm = TRUE),
        mean_log2fc = mean(abs(log2_fold_change), na.rm = TRUE),
        n_sig_transcripts = sum(adj_p_value < 0.05, na.rm = TRUE),
        n_transcripts = .N
    ), by = genes]
    
    # Convert back to data.frame for compatibility
    gene_ranking <- as.data.frame(gene_ranking)
    
    message("✓ Computed statistics for ", nrow(gene_ranking), " genes")
    
    # Sort genes: first by adjusted p-value (ascending), then by absolute log2_fold_change (descending)
    gene_ranking <- gene_ranking[order(gene_ranking$min_adj_pvalue, -gene_ranking$mean_log2fc, na.last = TRUE), ]
    
    message("Gene selection summary:")
    message("  Total genes available: ", nrow(gene_ranking))
    message("  Total transcripts available: ", sum(gene_ranking$n_transcripts))
    message("  Avg transcripts per gene: ", round(mean(gene_ranking$n_transcripts), 2))
    
    # Select top 50 most significant genes (lowest p-values)
    n_top_genes <- min(50, nrow(gene_ranking))
    top_genes <- gene_ranking$genes[1:n_top_genes]
    top_genes_tx_count <- sum(gene_ranking$n_transcripts[1:n_top_genes])
    selected_genes <- top_genes
    
    message("  Selected top ", n_top_genes, " genes with ", top_genes_tx_count, " transcripts")
    
    # Add 250 random genes from the remaining
    if (nrow(gene_ranking) > n_top_genes) {
        remaining_genes <- gene_ranking$genes[(n_top_genes + 1):nrow(gene_ranking)]
        if (length(remaining_genes) > 250) {
            random_genes <- sample(remaining_genes, 250)
        } else {
            random_genes <- remaining_genes
        }
        random_genes_tx_count <- sum(gene_ranking$n_transcripts[gene_ranking$genes %in% random_genes])
        selected_genes <- c(top_genes, random_genes)
        message("  Added ", length(random_genes), " random genes with ", random_genes_tx_count, " transcripts")
    } else {
        selected_genes <- top_genes
        message("  No random genes added (not enough remaining genes)")
    }
    
    message("  TOTAL: ", length(selected_genes), " selected genes with ", 
            sum(gene_ranking$n_transcripts[gene_ranking$genes %in% selected_genes]), " transcripts")
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 4 completed in ", round(elapsed, 2), " minutes")
    
    selected_genes
}

# =============================================================================
# STEP 4.5: CREATE GENCODE INDEXED CACHE (for fast GFF3 extraction)
# =============================================================================

create_gencode_index_cache <- function(config) {
    message("\n[STEP 4.5] Creating indexed gencode cache for fast GFF3 extraction...")
    step_start <- Sys.time()
    
    cache_file <- "./output/gencode_index_cache.RData"
    gencode_file <- config$gencode_file
    
    # Check if cache already exists
    if (file.exists(cache_file)) {
        cache_age_hours <- as.numeric(difftime(Sys.time(), file.mtime(cache_file), units = "hours"))
        if (cache_age_hours < 24) {
            message("✓ Using existing gencode cache (", round(cache_age_hours, 1), " hours old)")
            load(cache_file, envir = globalenv())
            return(gencode_gff_index)
        }
    }
    
    message("Building indexed cache from gencode file...")
    
    if (!file.exists(gencode_file)) {
        stop("ERROR: Gencode file not found: ", gencode_file)
    }
    
    # Extract only transcript and gene features (much faster than all lines)
    message("Extracting transcript and gene features from gencode (filtered extraction)...")
    
    grep_cmd <- paste0(
        "zcat ", shQuote(gencode_file), 
        " | awk -F'\\t' '$3 ~ /^(transcript|gene)$/ && $9 ~ /(transcript_id|gene_id)=/ {print}'"
    )
    all_gff_lines <- system(grep_cmd, intern = TRUE)
    
    n_lines <- length(all_gff_lines)
    message("  Extracted ", n_lines, " transcript/gene features")
    
    if (n_lines == 0) {
        stop("ERROR: No transcript/gene features found in gencode file")
    }
    
    # Parse GFF3 lines into data.table using chunked processing (memory-efficient)
    message("Parsing GFF3 lines into indexed data.table (chunked processing)...")
    
    chunk_size <- 10000
    gencode_gff_index <- NULL
    
    for (start_idx in seq(1, n_lines, by = chunk_size)) {
        end_idx <- min(start_idx + chunk_size - 1, n_lines)
        chunk_lines <- all_gff_lines[start_idx:end_idx]
        
        # Parse chunk into data frame first (faster than individual data.tables)
        chunk_data <- do.call(rbind, lapply(chunk_lines, function(line) {
            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) >= 9) {
                data.frame(
                    seqname = fields[1],
                    source = fields[2],
                    feature = fields[3],
                    start = as.numeric(fields[4]),
                    end = as.numeric(fields[5]),
                    score = fields[6],
                    strand = fields[7],
                    frame = fields[8],
                    attributes = fields[9],
                    gff_line = line,
                    stringsAsFactors = FALSE
                )
            } else {
                NULL
            }
        }))
        
        if (!is.null(chunk_data)) {
            chunk_dt <- as.data.table(chunk_data)
            
            # Combine with main table or initialize
            if (is.null(gencode_gff_index)) {
                gencode_gff_index <- chunk_dt
            } else {
                gencode_gff_index <- rbindlist(list(gencode_gff_index, chunk_dt), use.names = TRUE)
            }
            
            if (end_idx %% 50000 == 0) {
                message("    Processed ", end_idx, " / ", n_lines, " lines")
                gc()
            }
        }
    }
    
    message("  Parsed ", nrow(gencode_gff_index), " GFF3 features")
    
    # Extract transcript_id and gene_id from attributes
    message("Extracting transcript and gene IDs from attributes...")
    
    gencode_gff_index[, transcript_id := sub(".*transcript_id=([^;]+);.*", "\\1", attributes)]
    gencode_gff_index[, gene_id := sub(".*gene_id=([^;]+);.*", "\\1", attributes)]
    
    # Clean up: Set NA for failed extractions
    gencode_gff_index[!grepl("^ENST", transcript_id), transcript_id := NA_character_]
    gencode_gff_index[!grepl("^ENSG", gene_id), gene_id := NA_character_]
    
    # Remove rows without valid IDs
    gencode_gff_index <- gencode_gff_index[!is.na(gene_id)]
    
    message("Index prepared:")
    message("  Total features: ", nrow(gencode_gff_index))
    message("  Transcripts: ", nrow(gencode_gff_index[feature == "transcript"]))
    message("  Genes: ", nrow(gencode_gff_index[feature == "gene"]))
    
    # Set indices for fast lookup
    setindex(gencode_gff_index, gene_id)
    setindex(gencode_gff_index, transcript_id)
    
    # Save the cache
    message("Saving indexed cache to: ", cache_file)
    save(gencode_gff_index, file = cache_file, compress = "gzip")
    
    file_size_mb <- file.size(cache_file) / (1024 * 1024)
    message("✓ Cache created: ", round(file_size_mb, 2), " MB")
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 4.5 completed in ", round(elapsed, 2), " minutes")
    
    return(gencode_gff_index)
}

# =============================================================================
# STEP 5: FILTER AND OUTPUT DATA
# =============================================================================

filter_and_output <- function(readcounts, tx2gene, selected_genes, config, tpm = NULL, effective_length = NULL) {
    message("\n[STEP 5] Filtering counts and writing output files...")
    step_start <- Sys.time()
    
    # Filter transcripts for selected genes
    selected_tx <- tx2gene$transcript_id[tx2gene$gene_name %in% selected_genes]
    readcounts_filtered <- readcounts[selected_tx, ]
    tx2gene_filtered <- tx2gene[tx2gene$transcript_id %in% selected_tx, ]
    
    # Filter TPM and effective length matrices if provided
    if (!is.null(tpm)) {
        tpm_filtered <- tpm[selected_tx, ]
    } else {
        tpm_filtered <- NULL
    }
    
    if (!is.null(effective_length)) {
        effective_length_filtered <- effective_length[selected_tx]
    } else {
        effective_length_filtered <- NULL
    }
    
    # Extract unique genes
    selected_genes_unique <- sort(unique(selected_genes))
    n_genes_selected <- length(selected_genes_unique)
    
    # Extract sample IDs
    sample_geneset <- selected_genes_unique
    salmon_sample_IDs <- colnames(readcounts_filtered)
    
    message("Extracted ", length(sample_geneset), " unique genes from ", nrow(readcounts_filtered), " transcripts")
    message("Extracted ", length(salmon_sample_IDs), " sample IDs")
    
    # Summary statistics for gene filtering
    genes_filtered_count <- length(selected_genes)
    transcripts_filtered_count <- nrow(readcounts_filtered)
    avg_transcripts_per_gene <- round(transcripts_filtered_count / genes_filtered_count, 2)
    
    message("\nFinal Dataset Summary:")
    message("  Selected genes: ", genes_filtered_count)
    message("  Filtered transcripts: ", transcripts_filtered_count)
    message("  Average transcripts per gene: ", avg_transcripts_per_gene)
    
    if (transcripts_filtered_count < 300) {
        message("  WARNING: Expected ~300 transcripts for 300 genes, but got only ", transcripts_filtered_count)
        message("  This is expected if some selected genes have few transcripts or are underrepresented")
    }
    
    # Create inst/extdata directory if it doesn't exist
    if (!dir.exists(config$extdata_dir)) {
        dir.create(config$extdata_dir, recursive = TRUE)
        message("Created directory: ", config$extdata_dir)
    }
    
    # Write gene names to TSV file (unique genes only) in inst/extdata
    genes_file <- file.path(config$extdata_dir, "gene_ids.tsv")
    message("Writing gene names to: ", genes_file)
    write.table(sample_geneset, genes_file,
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
    )
    message("✓ Gene IDs written: ", genes_file)
    
    # Write tx2gene mapping
    tx2gene_file <- file.path(config$extdata_dir, "tx2gene.tsv")
    write.table(tx2gene_filtered,
                tx2gene_file,
                quote = FALSE, sep = "\t", row.names = FALSE)
    message("✓ TX2GENE mapping written: ", tx2gene_file, " (", nrow(tx2gene_filtered), " transcripts)")
    
    # Write RData to data folder
    readcounts <- as.data.frame(readcounts_filtered)
    rownames(readcounts) <- rownames(readcounts_filtered)
    
    # Prepare output objects
    tpm <- tpm_filtered
    effective_length <- effective_length_filtered
    
    rdata_file <- file.path("./data", "readcounts.RData")
    
    # Save all available objects
    if (!is.null(tpm) && !is.null(effective_length)) {
        save(readcounts, tpm, effective_length, file = rdata_file)
        message("✓ RData file written: ", rdata_file)
        message("  - Counts (readcounts)")
        message("  - TPM (tpm)")
        message("  - Effective lengths (effective_length)")
    } else {
        save(readcounts, file = rdata_file)
        message("✓ RData file written: ", rdata_file)
        message("  - Counts (readcounts)")
        if (is.null(tpm)) message("  ⚠ TPM not available")
        if (is.null(effective_length)) message("  ⚠ Effective lengths not available")
    }
    
    message("  Transcripts: ", nrow(readcounts_filtered), ", Samples: ", ncol(readcounts_filtered))
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 5 completed in ", round(elapsed, 2), " minutes")
    
    list(readcounts = readcounts_filtered, tx2gene = tx2gene_filtered)
}

# =============================================================================
# STEP 6: GENERATE GFF3 FILE
# =============================================================================

generate_gff3 <- function(readcounts_filtered, tx2gene, gencode_mapping, config, gencode_gff_index = NULL) {
    message("\n[STEP 6] Generating filtered GFF3 file...")
    step_start <- Sys.time()
    
    # Get transcript IDs from filtered readcounts
    unique_tx_ids_filtered <- unique(rownames(readcounts_filtered))
    n_tx_filtered <- length(unique_tx_ids_filtered)
    
    message("Extracting annotations for ", n_tx_filtered, " filtered transcripts...")
    
    if(n_tx_filtered == 0) {
        stop("ERROR: No transcripts found in readcounts_filtered. ",
             "This might indicate an issue with filtering in STEP 3-4.")
    }
    
    # Use indexed cache if available, otherwise load it
    if (is.null(gencode_gff_index)) {
        cache_file <- "./output/gencode_index_cache.RData"
        if (file.exists(cache_file)) {
            message("Loading gencode index cache...")
            load(cache_file, envir = environment())
        } else {
            stop("ERROR: Gencode index cache not found. Run create_gencode_index_cache() first.")
        }
    }
    
    # Fast data.table lookup: filter by transcript IDs (1-2ms for 7.6M records)
    message("Fast lookup: extracting transcript features from index...")
    
    # Remove version numbers for matching (ENST000001.1 -> ENST000001)
    unique_tx_ids_clean <- sub("\\.[0-9]+$", "", unique_tx_ids_filtered)
    gencode_gff_index[, tx_clean := sub("\\.[0-9]+$", "", transcript_id)]
    
    # Filter by matched transcripts
    filtered_tx_dt <- gencode_gff_index[feature == "transcript" & 
                                        (transcript_id %in% unique_tx_ids_filtered | 
                                         tx_clean %in% unique_tx_ids_clean)]
    
    message("✓ Extracted ", nrow(filtered_tx_dt), " transcript features")
    
    if (nrow(filtered_tx_dt) == 0) {
        stop("ERROR: No transcript features found for selected transcripts. ",
             "Sample filtered IDs: ", paste(head(unique_tx_ids_filtered, 3), collapse = ", "))
    }
    
    # Get unique gene IDs from selected transcripts
    unique_gene_ids <- unique(filtered_tx_dt[!is.na(gene_id), gene_id])
    message("Found ", length(unique_gene_ids), " unique genes")
    
    # Extract gene-level features (fast lookup)
    message("Extracting gene-level features from index...")
    gene_dt <- gencode_gff_index[feature == "gene" & gene_id %in% unique_gene_ids]
    message("✓ Extracted ", nrow(gene_dt), " gene features")
    
    # Combine transcript and gene lines
    output_dt <- rbindlist(list(filtered_tx_dt, gene_dt), use.names = TRUE, fill = TRUE)
    
    # Debug: Check if gff_line column exists
    if (!"gff_line" %in% names(output_dt)) {
        stop("ERROR: gff_line column not found in output_dt. ",
             "Available columns: ", paste(names(output_dt), collapse = ", "))
    }
    
    # Extract GFF3 lines and write output
    message("Writing GFF3 file...")
    temp_gff_file <- "output/annotation.gff3.tmp"
    header_line <- "##gff-version 3"
    
    # Write header + data lines
    output_lines <- c(header_line, output_dt$gff_line)
    message("  Total output lines: ", length(output_lines))
    writeLines(output_lines, temp_gff_file)
    
    message("Compressed output to: ", config$output_gff)
    system(paste("gzip -c", shQuote(temp_gff_file), ">", shQuote(config$output_gff)))
    
    # Cleanup
    if (file.exists(temp_gff_file)) unlink(temp_gff_file)
    
    # Validate output
    gff_lines_count <- system(paste("zcat", shQuote(config$output_gff), "| wc -l"), intern = TRUE)
    total_lines <- as.numeric(strsplit(gff_lines_count, " ")[[1]][1])
    data_lines <- total_lines - 1  # Exclude header
    
    message("✓ GFF3 file successfully written!")
    message("  Header: 1 line")
    message("  Data lines: ", data_lines, " (", nrow(filtered_tx_dt), " transcripts + ", 
            nrow(gene_dt), " genes)")
    
    # Cleanup: Remove temp cache from memory
    gencode_gff_index <<- NULL
    gc()
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ STEP 6 completed in ", round(elapsed, 2), " minutes")
}

# =============================================================================
# STEP 7-8: VALIDATION
# =============================================================================

validate_outputs <- function(readcounts_filtered, config) {
    message("\n[STEP 7-8] Validating output files...")
    step_start <- Sys.time()
    
    # Count transcripts in outputs
    n_counts <- nrow(readcounts_filtered)
    
    # Count transcripts in GFF3
    if (file.exists(config$output_gff)) {
        n_gff3 <- system(paste("zcat", shQuote(config$output_gff), 
                               "| grep -v '^#' | grep 'transcript' | wc -l"), 
                         intern = TRUE)
        n_gff3 <- as.numeric(trimws(n_gff3))
        
        if (n_gff3 == n_counts) {
            message("✓ GFF3 and readcounts synchronized: ", n_gff3, " transcripts")
        } else {
            message("⚠ Warning: GFF3 (", n_gff3, ") and readcounts (", n_counts, 
                    ") have different transcript counts")
            
            # Extract and validate transcript IDs match
            readcounts_tx <- rownames(readcounts_filtered)
            
            # Extract transcript IDs from GFF3
            gff_lines <- system(paste("zcat", shQuote(config$output_gff), "| grep -v '^#'"), 
                               intern = TRUE)
            gff_tx <- sub(".*transcript_id=([^;]+);.*", "\\1", gff_lines)
            gff_tx <- unique(gff_tx[grepl("^ENST", gff_tx)])
            
            # Check matches
            missing_in_gff3 <- setdiff(readcounts_tx, gff_tx)
            missing_in_readcounts <- setdiff(gff_tx, readcounts_tx)
            
            if (length(missing_in_gff3) > 0) {
                message("  WARNING: ", length(missing_in_gff3), 
                       " transcripts in readcounts but NOT in GFF3")
                message("    Sample: ", paste(head(missing_in_gff3, 3), collapse = ", "))
            }
            if (length(missing_in_readcounts) > 0) {
                message("  WARNING: ", length(missing_in_readcounts), 
                       " transcripts in GFF3 but NOT in readcounts")
                message("    Sample: ", paste(head(missing_in_readcounts, 3), collapse = ", "))
            }
        }
    }
    
    elapsed <- difftime(Sys.time(), step_start, units = "mins")
    message("✓ Validation completed in ", round(elapsed, 2), " minutes")
}

# =============================================================================
# GENERATE TEST DATASETS
# =============================================================================

generate_test_gff3 <- function(config) {
    message("\n[TEST] Generating test GFF3 file...")
    
    test_gff3 <- c(
        "##gff-version 3",
        "1\tGENCODE\tgene\t1000\t10000\t.\t+\t.\tID=ENSG00000101456;Name=MXRA8",
        "1\tGENCODE\ttranscript\t1000\t10000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000101456",
        "1\tGENCODE\ttranscript\t1500\t9500\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000101456",
        "2\tGENCODE\tgene\t5000\t15000\t.\t-\t.\tID=ENSG00000102458;Name=C1orf86",
        "2\tGENCODE\ttranscript\t5000\t15000\t.\t-\t.\tID=ENST00000003;Parent=ENSG00000102458",
        "2\tGENCODE\ttranscript\t5500\t14500\t.\t-\t.\tID=ENST00000004;Parent=ENSG00000102458",
        "2\tGENCODE\ttranscript\t6000\t14000\t.\t-\t.\tID=ENST00000005;Parent=ENSG00000102458",
        "3\tGENCODE\tgene\t20000\t30000\t.\t+\t.\tID=ENSG00000103259;Name=PDPN",
        "3\tGENCODE\ttranscript\t20000\t30000\t.\t+\t.\tID=ENST00000006;Parent=ENSG00000103259",
        "3\tGENCODE\ttranscript\t20500\t29500\t.\t+\t.\tID=ENST00000007;Parent=ENSG00000103259",
        "3\tGENCODE\ttranscript\t21000\t29000\t.\t+\t.\tID=ENST00000008;Parent=ENSG00000103259"
    )
    
    test_gff_file <- file.path(config$extdata_dir, "gencode_subset_test.gff3")
    writeLines(test_gff3, test_gff_file)
    
    output_file <- paste0(test_gff_file, ".gz")
    system(paste("gzip -c", shQuote(test_gff_file), ">", shQuote(output_file)))
    unlink(test_gff_file)
    
    message("✓ Test GFF3: ", output_file)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

main <- function() {
    config <- parse_config()
    print_config(config)
    check_directories(config)
    
    # Download metadata from Zenodo if needed
    metadata_available <- download_metadata_from_zenodo(config)
    if (!metadata_available) {
        stop("ERROR: Metadata file not available. Cannot proceed without metadata.")
    }
    
    # Download gencode from EBI if needed
    gencode_available <- download_gencode_from_ebi(config)
    if (!gencode_available) {
        stop("ERROR: Gencode file not available. Cannot proceed without gencode annotation.")
    }
    
    # Download salmon samples from Zenodo if needed
    salmon_dir <- download_salmon_samples_from_zenodo(config)
    
    # STEP 1: Read Salmon data
    salmon_data <- read_salmon_samples(salmon_dir, config$metadata_file)
    readcounts <- salmon_data$counts
    tpm <- salmon_data$tpm
    effective_length <- salmon_data$effective_length
    metadata <- salmon_data$metadata
    
    # STEP 2: Map transcripts to genes
    mapping_result <- map_transcripts_to_genes(readcounts, config$gencode_file)
    tx2gene <- mapping_result$tx2gene
    gencode_mapping <- mapping_result$mapping
    
    # STEP 2.5: Filter single-transcript genes BEFORE diversity calculation
    filter_result <- filter_single_transcript_genes(readcounts, tx2gene)
    readcounts <- filter_result$readcounts
    tx2gene <- filter_result$tx2gene
    
    # STEP 3: Calculate diversity
    diversity_data <- calculate_diversity_metrics(readcounts, tx2gene, metadata)
    
    # STEP 4: Select top genes
    selected_genes <- select_top_genes(readcounts, tx2gene, diversity_data)
    
    # STEP 4.5: Create indexed gencode cache for fast GFF3 extraction
    gencode_gff_index <- create_gencode_index_cache(config)
    
    # STEP 5: Filter and output
    filtered_result <- filter_and_output(readcounts, tx2gene, selected_genes, config, tpm, effective_length)
    readcounts_filtered <- filtered_result$readcounts
    
    # STEP 6: Generate GFF3 (using cached index)
    generate_gff3(readcounts_filtered, tx2gene, gencode_mapping, config, gencode_gff_index)
    
    # STEP 7-8: Validate
    validate_outputs(readcounts_filtered, config)
    
    # Generate test data
    generate_test_gff3(config)
    
    # Summary
    total_time <- difftime(Sys.time(), start_time, units = "mins")
    message("\n=== PREPROCESSING COMPLETE ===")
    message("Total time: ", round(total_time, 2), " minutes")
    message("Files created:")
    message("  - ./data/readcounts.RData")
    message("  - ", config$output_gff)
    message("  - ", file.path(config$extdata_dir, "tx2gene.tsv"))
}

# Run pipeline
main()
