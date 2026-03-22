## Helper: Access gene ID column from rowData
#' @title Get Gene IDs from SummarizedExperiment rowData
#'
#' @description Internal helper function to safely access gene ID information
#' from a SummarizedExperiment's rowData.
#'
#' @param se A SummarizedExperiment object
#'
#' @return A character vector of gene IDs (same length as nrow(se)),
#'   or NULL if no gene ID column is found.
#'
#' @details
#' Searches rowData in this priority order:
#' 1. "gene_id" (standard column from build_se)
#' 2. "gene_name" (fallback for GFF3-derived names)
#'
#' @keywords internal
#' @noRd
.get_gene_ids <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  if (is.null(rd)) {
    return(NULL)
  }
  
  # Try gene_id column first
  if ("gene_id" %in% colnames(rd)) {
    return(as.character(rd$gene_id))
  }
  
  # Fall back to gene_name if available
  if ("gene_name" %in% colnames(rd)) {
    return(as.character(rd$gene_name))
  }
  
  return(NULL)
}

## Helper: Extract tx2gene mapping from GFF3 file
#' @title Extract Transcript-to-Gene Mapping from GFF3 File
#' @description Internal function that parses GFF3 or GFF3.gz files to extract
#' transcript-to-gene mappings. This function is called internally by \code{build_se()}
#' when a GFF3 file path is provided. Users should not call this function directly.
#' @param gff3_file Path to a GFF3 or GFF3.gz file containing transcript/mRNA
#'   features with ID and Parent attributes.
#' @return A data.frame with two columns (Transcript, Gene) containing the
#'   transcript-to-gene mapping extracted from the GFF3 file.
#' @details
#' The function expects GFF3 format with the following structure:
#' - 9 tab-separated columns: seqname, source, feature, start, end,
#'   score, strand, phase, attributes
#' - Feature type column (3rd column) should contain 'transcript' or 'mRNA'
#' - Attributes column (9th column) should contain ID and Parent fields
#' - ID field: unique identifier for the transcript
#' - Parent field: references the gene ID that this transcript belongs to
#' Example GFF3 line:
#' \preformatted{chr1\tgencode\ttranscript\t1000\t3000\t.\t+\t.\t
#' ID=ENST00000001;Parent=ENSG00000101;Name=BRCA1-001}
#'
#' \strong{Performance:} This function is optimized for large GFF3 files:
#' - Reads files in 10,000-line chunks (not line-by-line)
#' - Uses fast pre-filtering (feature type check before regex)
#' - Employs efficient string operations instead of heavy regex on every line
#' - Handles both compressed (.gz) and uncompressed files seamlessly
#' @examples
#' \dontrun{
#' # Create a temporary GFF3 file with transcript features
#' gff3_lines <- c(
#'   'chr1\tgencode\ttranscript\t1000\t3000\t.\t+\t.\tID=ENST001;Parent=ENSG001',
#'   'chr1\tgencode\ttranscript\t1500\t3500\t.\t+\t.\tID=ENST002;Parent=ENSG001',
#'   'chr2\tgencode\ttranscript\t5000\t8000\t.\t-\t.\tID=ENST003;Parent=ENSG002'
#' )
#' tf <- tempfile(fileext = '.gff3')
#' writeLines(gff3_lines, tf)
#' 
#' # Extract transcript-to-gene mapping
#' tx2gene <- extract_tx2gene_from_gff3(tf)
#' head(tx2gene)
#' }
#' @noRd
extract_tx2gene_from_gff3 <- function(gff3_file) {
    # Handle both .gff3 and .gff3.gz files
    if (grepl("\\.gff3\\.gz$", gff3_file)) {
        con <- gzfile(gff3_file, "rt")
    } else {
        con <- file(gff3_file, "r")
    }

    on.exit(close(con))

    # Pre-allocate vectors for efficiency (start with capacity for 10k
    # transcripts) This avoids repeatedly growing lists/data.frames which is
    # slow
    tx2gene_transcripts <- character(10000)
    tx2gene_genes <- character(10000)
    idx <- 0
    chunk_size <- 10000  # Read in chunks for better performance

    while (TRUE) {
        # Read lines in chunks instead of one-by-one for better performance
        lines <- readLines(con, n = chunk_size)
        if (length(lines) == 0)
            break

        # Process each line in the chunk
        for (line in lines) {
            # Skip comments and empty lines (fast pre-filter)
            if (startsWith(line, "#") || line == "")
                next

            # Parse GFF3 line format: seqname source feature start end score
            # strand phase attributes Use strsplit only once per line
            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) < 9)
                next

            # Extract feature type (3rd column) - do this check first to skip
            # early
            feature_type <- fields[3]
            if (!feature_type %in% c("transcript", "mRNA"))
                next

            # Only now extract attributes from the (potentially long) 9th
            # column
            attributes <- fields[9]

            # Extract ID and Parent using efficient substring operations Find
            # positions of ID= and Parent= patterns
            id_start <- regexpr("ID=", attributes, fixed = TRUE) + 3
            if (id_start > 3) {
                # ID found, extract until semicolon or end of string
                id_end <- regexpr(";", substr(attributes, id_start, nchar(attributes)),
                  fixed = TRUE)
                if (id_end > 0) {
                  transcript_id <- substr(attributes, id_start, id_start + id_end -
                    2)
                } else {
                  transcript_id <- substr(attributes, id_start, nchar(attributes))
                }
            } else {
                transcript_id <- NA_character_
            }

            parent_start <- regexpr("Parent=", attributes, fixed = TRUE) + 7
            if (parent_start > 7) {
                # Parent found, extract until semicolon or end of string
                parent_end <- regexpr(";", substr(attributes, parent_start, nchar(attributes)),
                  fixed = TRUE)
                if (parent_end > 0) {
                  gene_id <- substr(attributes, parent_start, parent_start + parent_end -
                    2)
                } else {
                  gene_id <- substr(attributes, parent_start, nchar(attributes))
                }
            } else {
                gene_id <- NA_character_
            }

            # Store mapping if both IDs are present
            if (!is.na(transcript_id) && !is.na(gene_id)) {
                idx <- idx + 1
                # Grow vectors if needed
                if (idx > length(tx2gene_transcripts)) {
                  tx2gene_transcripts <- c(tx2gene_transcripts, rep(NA_character_,
                    10000))
                  tx2gene_genes <- c(tx2gene_genes, rep(NA_character_, 10000))
                }
                # Strip common prefixes (transcript:, gene:, etc.)
                transcript_id <- sub("^transcript:", "", transcript_id)
                gene_id <- sub("^gene:", "", gene_id)
                tx2gene_transcripts[idx] <- transcript_id
                tx2gene_genes[idx] <- gene_id
            }
        }
    }

    # Trim to actual size and create data frame (empty if no mappings found)
    if (idx == 0) {
        tx2gene_df <- data.frame(Transcript = character(0), Gene = character(0),
            stringsAsFactors = FALSE)
    } else {
        tx2gene_df <- data.frame(Transcript = tx2gene_transcripts[seq_len(idx)],
            Gene = tx2gene_genes[seq_len(idx)], stringsAsFactors = FALSE)
    }
    rownames(tx2gene_df) <- NULL

    return(tx2gene_df)
}

## Helper: Extract gene names from GFF3 file
extract_gene_names_from_gff3 <- function(gff3_file) {
    # Handle both .gff3 and .gff3.gz files
    if (grepl("\\.gff3\\.gz$", gff3_file)) {
        con <- gzfile(gff3_file, "rt")
    } else {
        con <- file(gff3_file, "r")
    }

    on.exit(close(con))

    # Pre-allocate vectors for efficiency
    gene_ids <- character(5000)
    gene_names <- character(5000)
    idx <- 0
    chunk_size <- 10000

    while (TRUE) {
        lines <- readLines(con, n = chunk_size)
        if (length(lines) == 0)
            break

        for (line in lines) {
            # Skip comments and empty lines
            if (startsWith(line, "#") || line == "")
                next

            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) < 9)
                next

            # Only process gene features
            feature_type <- fields[3]
            if (feature_type != "gene")
                next

            attributes <- fields[9]

            # Extract ID field
            id_start <- regexpr("ID=", attributes, fixed = TRUE) + 3
            if (id_start <= 3)
                next
            id_end <- regexpr(";", substr(attributes, id_start, nchar(attributes)),
                fixed = TRUE)
            if (id_end > 0) {
                gene_id <- substr(attributes, id_start, id_start + id_end - 2)
            } else {
                gene_id <- substr(attributes, id_start, nchar(attributes))
            }

            # Extract gene_name field (GENCODE uses gene_name= instead of Name=)
            name_start <- regexpr("gene_name=", attributes, fixed = TRUE) + 10
            if (name_start > 10) {
                name_end <- regexpr(";", substr(attributes, name_start, nchar(attributes)),
                  fixed = TRUE)
                if (name_end > 0) {
                  gene_name <- substr(attributes, name_start, name_start + name_end - 2)
                } else {
                  gene_name <- substr(attributes, name_start, nchar(attributes))
                }
            } else {
                gene_name <- NA_character_
            }

            # Store mapping if we have at least a gene ID
            if (!is.na(gene_id)) {
                idx <- idx + 1
                # Grow vectors if needed
                if (idx > length(gene_ids)) {
                  gene_ids <- c(gene_ids, rep(NA_character_, 5000))
                  gene_names <- c(gene_names, rep(NA_character_, 5000))
                }
                # Strip gene: prefix
                gene_id <- sub("^gene:", "", gene_id)
                gene_ids[idx] <- gene_id
                gene_names[idx] <- gene_name
            }
        }
    }

    # Trim to actual size and create data frame
    if (idx == 0) {
        res <- data.frame(GeneID = character(0), GeneName = character(0),
            stringsAsFactors = FALSE)
    } else {
        res <- data.frame(GeneID = gene_ids[seq_len(idx)], GeneName = gene_names[seq_len(idx)],
            stringsAsFactors = FALSE)
    }
    rownames(res) <- NULL

    return(res)
}

## Helper: build SummarizedExperiment from readcounts + tx2gene
#' Build a SummarizedExperiment from transcript readcounts and tx->gene map
#'
#' This helper creates a `SummarizedExperiment` with an assay named
#' `counts`, stores the `tx2gene` table and the raw `readcounts` in
#' `metadata()`, and extracts gene assignments from the tx2gene mapping
#' to populate `rowData()` with both transcript IDs and gene IDs.
#'
#' **NEW:** Also accepts SALMON quantification data (TPM, effective_length) 
#' for downstream filtering and length-normalized diversity analysis.
#'
#' The function accepts three different input types for the tx2gene mapping:
#' a file path (TSV or GFF3 format) or an in-memory data.frame.
#'
#' @param readcounts Numeric matrix or data.frame of transcript-level counts
#'   (rows = transcripts, columns = samples). Rownames should contain
#'   transcript IDs that match the first column of tx2gene.
#'
#' @param tx2gene Transcript-to-gene mapping. Can be one of:
#'   \describe{
#'     \item{\strong{TSV File Path}}{A path to a tab-separated file with at least
#'       two columns: Transcript ID (first column) and Gene ID (second column).
#'       The file should have a header row with column names.
#'       File extension should be .tsv, .txt, or similar.
#'       Example:
#'       \preformatted{Transcript\tGene
#'       ENST00000001\tENSG00000101
#'       ENST00000002\tENSG00000102}
#'     }
#'     \item{\strong{GFF3 File Path}}{A path to a GFF3 annotation file (.gff3 or
#'       .gff3.gz format). The file should contain transcript/mRNA features with
#'       ID and Parent attributes. Parent attributes reference gene IDs.
#'       Example GFF3 line:
#'       \preformatted{chr1\tgencode\ttranscript\t100\t2100\t.\t+\t.\t
#'       ID=ENST00000001;Parent=ENSG00000101}
#'       The function automatically detects GFF3 format by file extension
#'       (.gff3 or .gff3.gz) and parses accordingly.
#'     }
#'     \item{\strong{Data Frame}}{An in-memory data.frame with Transcript
#'       and Gene columns. Useful when mapping data is already loaded in R.
#'       Example:
#'       \preformatted{  Transcript        Gene
#'       1 ENST00000001 ENSG00000101
#'       2 ENST00000002 ENSG00000102}
#'     }
#'   }
#'
#' @param assay_name Name for the assay to store readcounts (default: 'counts').
#'
#' @param tpm Numeric matrix or data.frame of SALMON TPM values (optional).
#'   If provided, stored in metadata as `salmon_tpm` for `filter_se()` 
#'   TPM-based filtering. Rows = transcripts, columns = samples (same dimension as readcounts).
#'   Example: from preprocessing output `salmon_tpm` or loaded via `load('readcounts.RData')`.
#'   
#' @param effective_length Numeric vector of effective transcript lengths (optional).
#'   If provided, stored in metadata as `salmon_effective_length` for `calculate_diversity()`
#'   length-normalized entropy calculations. Vector length = number of transcripts.
#'   Typically obtained from SALMON quantification's EffectiveLength column (median across samples).
#'   Example: from preprocessing output `salmon_effective_length` or loaded via `load('readcounts.RData')`.
#'   
#' @param skip Logical. If TRUE, unmapped transcripts are silently removed. 
#'   If FALSE (default), an error is raised when unmapped transcripts are found.
#'
#' @return A `SummarizedExperiment` with:
#'   - `assay (counts)`: raw transcript counts
#'   - `metadata$tx2gene`: transcript-to-gene mapping
#'   - `metadata$readcounts`: raw transcript counts (preserved)
#'   - `metadata$salmon_tpm`: TPM values (if provided)
#'   - `metadata$salmon_effective_length`: effective lengths (if provided)
#'   - `rowData$transcript_id`: transcript IDs (matching rownames)
#'   - `rowData$gene_id`: gene IDs for each transcript
#'   - `rowData$gene_name`: human-readable gene names (if GFF3 provided)
#'
#' @details
#' \strong{Input Format Detection:}
#' - If tx2gene is a character string ending in .gff3 or .gff3.gz,
#'   it is parsed as a GFF3 file.
#' - If tx2gene is a character string with any other extension or
#'   no extension, it is parsed as a tab-separated file.
#' - If tx2gene is a data.frame, it is used directly.
#'
#' \strong{rowData Structure (NEW):}
#' The rowData now contains both transcript-level and gene-level identifiers:
#' - `transcript_id`: Transcript identifier (from readcounts rownames)
#' - `gene_id`: Gene identifier (from tx2gene mapping)
#' - `gene_name`: Human-readable gene name (from GFF3 when available)
#' This enables downstream functions like `calculate_divergence()` to easily
#' filter by gene ID and aggregate or subset transcripts by gene.
#'
#' \strong{Transcript Matching:}
#' All transcript IDs in the readcounts object (rownames) must have a
#' corresponding entry in the tx2gene mapping. If any transcript is missing,
#' an error is raised.
#'
#' \strong{SALMON Data Integration:}
#' When TPM and effective_length are provided, they are stored in metadata
#' for seamless integration with:
#' - `filter_se()`: automatically detects TPM in metadata for 
#'   normalization-aware filtering
#' - `calculate_diversity()`: automatically detects effective_length in metadata
#'   for length-normalized entropy calculations
#'
#' \strong{Performance:}
#' - GFF3 files are processed efficiently even for large annotations
#'   (e.g., full GENCODE with 100k+ transcripts)
#' - TSV files are standard tab-separated format for fast parsing
#' - Data.frame inputs have no I/O overhead
#'
#' @param metadata Optional list or named list of metadata to include in \code{metadata(se)}.
#'   This is useful for storing additional experimental metadata alongside the SE object.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' # Example 1: Using data.frame (in-memory mapping)
#' tx2gene <- data.frame(
#'     Transcript = c('ENST00000001', 'ENST00000002'),
#'     Gene = c('ENSG00000101', 'ENSG00000102')
#' )
#' readcounts <- matrix(c(10, 5, 2, 3),
#'     nrow = 2,
#'     dimnames = list(c('ENST00000001', 'ENST00000002'), c('s1', 's2'))
#' )
#' se <- build_se(readcounts, tx2gene)
#' # Now rowData contains transcript_id and gene_id columns:
#' # rowData(se)$transcript_id  # ENST00000001, ENST00000002
#' # rowData(se)$gene_id        # ENSG00000101, ENSG00000102
#'
#' # Example 2: With SALMON data (TPM, effective_length)
#' # Assuming preprocessed SALMON output loaded:
#' # load('readcounts.RData') # salmon_dataset, salmon_tpm, salmon_effective_length
#' # se <- build_se(salmon_dataset, tx2gene, tpm = salmon_tpm, 
#' #               effective_length = salmon_effective_length)
#' # Now filter_se() and calculate_diversity() use SALMON data automatically:
#' # filtered_se <- filter_se(se, stringency = 'medium')  # Uses TPM from metadata
#' # div_se <- calculate_diversity(salmon_dataset, ...)   # Uses effective_length
#'
#' # Example 3: Using TSV file path
#' # Assuming you have a file 'tx2gene.tsv' with Transcript and Gene columns
#' # se <- build_se(readcounts, 'path/to/tx2gene.tsv')
#' # Example 4: Using GFF3.gz file path
#' # Assuming you have a file 'annotation.gff3.gz' with transcript features
#' # se <- build_se(readcounts, 'path/to/annotation.gff3.gz')
build_se <- function(readcounts, tx2gene, assay_name = "counts", skip = FALSE, 
                     tpm = NULL, effective_length = NULL, metadata = NULL) {
    # Auto-detect TPM and effective_length from readcounts.RData Global Environment
    # Look for matrices matching readcounts dimensions (for TPM)
    # and vectors matching readcounts rows (for effective_length)
    # This is more flexible than hardcoding specific variable names
    # Only auto-detect if both readcounts is not NULL and is already a matrix
    
    if (is.null(tpm) && !is.null(readcounts)) {
        tryCatch({
            # Convert readcounts to matrix if needed for dimension checking
            rc_matrix <- if (is.matrix(readcounts)) readcounts else as.matrix(readcounts)
            
            # Search Global Environment for TPM-like variables
            # Look for: matrix/data.frame with same dimensions as readcounts, and "tpm" in name
            if (exists(".GlobalEnv")) {
                env_vars <- ls(envir = .GlobalEnv)
                for (var_name in env_vars) {
                    var <- get(var_name, envir = .GlobalEnv)
                    # Check if it's TPM-like: matrix/data.frame with matching dimensions AND "tpm" in name
                    if ((is.matrix(var) || is.data.frame(var)) && 
                        nrow(var) == nrow(rc_matrix) && ncol(var) == ncol(rc_matrix) &&
                        grepl("tpm", tolower(var_name), ignore.case = TRUE)) {
                        tpm <- var
                        break
                    }
                }
            }
        }, error = function(e) {
            # Silently ignore errors during auto-detection (e.g., during roxygen2 processing)
        })
    }
    
    if (is.null(effective_length) && !is.null(readcounts)) {
        tryCatch({
            # Convert readcounts to matrix if needed for dimension checking
            rc_matrix <- if (is.matrix(readcounts)) readcounts else as.matrix(readcounts)
            
            # Search Global Environment for effective_length-like variables
            # Look for: numeric vector with length matching nrow(readcounts), and "length" or "eff" in name
            if (exists(".GlobalEnv")) {
                env_vars <- ls(envir = .GlobalEnv)
                for (var_name in env_vars) {
                    var <- get(var_name, envir = .GlobalEnv)
                    # Check if it's effective_length-like: numeric vector matching nrow AND has "length"/"eff" in name
                    if ((is.numeric(var) && !is.matrix(var)) && 
                        length(var) == nrow(rc_matrix) &&
                        (grepl("length", tolower(var_name), ignore.case = TRUE) ||
                         grepl("eff", tolower(var_name), ignore.case = TRUE))) {
                        effective_length <- var
                        break
                    }
                }
            }
        }, error = function(e) {
            # Silently ignore errors during auto-detection (e.g., during roxygen2 processing)
        })
    }
    
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }

        # Detect file type and parse accordingly
        if (grepl("\\.gff3(\\.gz)?$", tx2gene, ignore.case = TRUE)) {
            message("Detected GFF3 format. Extracting transcript-to-gene mapping...")
            tx2gene_df <- extract_tx2gene_from_gff3(tx2gene)
        } else {
            # Assume TSV format (backward compatible)
            tx2gene_df <- utils::read.table(tx2gene, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        }
    } else if (is.data.frame(tx2gene)) {
        tx2gene_df <- tx2gene
    } else {
        stop("'tx2gene' must be a path (TSV or GFF3) or a data.frame.", call. = FALSE)
    }

    if (is.data.frame(readcounts)) {
        readcounts <- as.matrix(readcounts)
    }
    if (!is.matrix(readcounts) || !is.numeric(readcounts)) {
        stop("'readcounts' must be a numeric matrix or numeric data.frame.", call. = FALSE)
    }

    # Extract transcript IDs from rownames or use row indices
    tx_ids <- rownames(readcounts)
    if (is.null(tx_ids)) {
        stop("'readcounts' must have transcript IDs as rownames.", call. = FALSE)
    }

    # Extract gene names from tx2gene based on transcript IDs
    tx_col <- if ("Transcript" %in% colnames(tx2gene_df))
        "Transcript" else colnames(tx2gene_df)[1]
    genes <- tx2gene_df$Gene[match(tx_ids, tx2gene_df[[tx_col]])]

    # Check for unmapped transcripts
    unmapped_idx <- which(is.na(genes))
    if (length(unmapped_idx) > 0) {
        unmapped_txs <- tx_ids[unmapped_idx]
        
        if (!skip) {
            message(length(unmapped_txs), " transcript IDs were not found in tx2gene mapping.")
            message("Unmapped transcripts: ", paste(head(unmapped_txs, 10), collapse = ", "),
                    if (length(unmapped_txs) > 10) paste0(" ... and ", length(unmapped_txs) - 10, " more") else "")
            stop("Unmapped transcripts detected. Set skip=TRUE to remove them and continue.", call. = FALSE)
        } else {
            message(length(unmapped_txs), " transcript IDs were not found in tx2gene mapping.")
            message("Unmapped transcripts: ", paste(head(unmapped_txs, 10), collapse = ", "),
                    if (length(unmapped_txs) > 10) paste0(" ... and ", length(unmapped_txs) - 10, " more") else "")
            
            # When majority or all transcripts are unmapped, use them directly as gene identifiers
            # This provides a sensible fallback for cases where the annotation doesn't match the data
            if (length(unmapped_idx) >= length(tx_ids) * 0.9) {
                message("Note: >90% of transcripts unmapped. Using transcript IDs as gene identifiers.")
                genes <- tx_ids  # Use transcript IDs directly as gene identifiers
            } else {
                message("Removing unmapped transcripts from analysis (skip=TRUE).")
                # Filter out unmapped transcripts only if minority are unmapped
                keep_idx <- which(!is.na(genes))
                readcounts <- readcounts[keep_idx, , drop = FALSE]
                tx_ids <- tx_ids[keep_idx]
                genes <- genes[keep_idx]
            }
        }
    }

    # Extract gene names if we're using a GFF3 file
    gene_names_df <- NULL
    if (is.character(tx2gene) && grepl("\\.gff3(\\.gz)?$", tx2gene, ignore.case = TRUE)) {
        gene_names_df <- extract_gene_names_from_gff3(tx2gene)
    }

    assays_list <- S4Vectors::SimpleList()
    assays_list[[assay_name]] <- readcounts

    se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list)
    S4Vectors::metadata(se)$tx2gene <- tx2gene_df
    S4Vectors::metadata(se)$readcounts <- readcounts

    # ========================================================================
    # SALMON DATA: Store TPM and effective_length in metadata
    # ========================================================================
    # These are used by filter_se() for TPM-based filtering
    # and calculate_diversity() for length-normalized entropy
    
    # Store TPM if provided
    if (!is.null(tpm)) {
        if (is.data.frame(tpm)) {
            tpm <- as.matrix(tpm)
        }
        if (!is.matrix(tpm) || !is.numeric(tpm)) {
            stop("'tpm' must be a numeric matrix or data.frame.", call. = FALSE)
        }
        # Validate dimensions match readcounts
        if (nrow(tpm) != nrow(readcounts) || ncol(tpm) != ncol(readcounts)) {
            stop(sprintf("'tpm' dimensions (%d x %d) do not match 'readcounts' (%d x %d).",
                        nrow(tpm), ncol(tpm), nrow(readcounts), ncol(readcounts)),
                 call. = FALSE)
        }
        # Sync rownames to readcounts
        rownames(tpm) <- rownames(readcounts)
        S4Vectors::metadata(se)$salmon_tpm <- tpm
    }
    
    # Store effective_length if provided
    if (!is.null(effective_length)) {
        if (!is.numeric(effective_length)) {
            stop("'effective_length' must be numeric.", call. = FALSE)
        }
        # Validate length matches readcounts
        if (length(effective_length) != nrow(readcounts)) {
            stop(sprintf("'effective_length' length (%d) does not match 'readcounts' rows (%d).",
                        length(effective_length), nrow(readcounts)),
                 call. = FALSE)
        }
        # Ensure names match readcounts rownames
        if (is.null(names(effective_length))) {
            names(effective_length) <- rownames(readcounts)
        } else {
            # Check if names match (reorder if necessary)
            if (!all(names(effective_length) == rownames(readcounts))) {
                # Try to reorder if names exist but are in different order
                if (all(rownames(readcounts) %in% names(effective_length))) {
                    effective_length <- effective_length[rownames(readcounts)]
                } else {
                    warning("Names in 'effective_length' do not match 'readcounts' rownames. ",
                           "Assigning rownames without reordering.",
                           call. = FALSE)
                    names(effective_length) <- rownames(readcounts)
                }
            }
        }
        S4Vectors::metadata(se)$salmon_effective_length <- effective_length
    }

    # Populate rowData with gene assignments and transcript IDs
    # Includes transcript_id and gene_id for gene annotation
    # This enables downstream functions to easily access either identifier
    if (!is.null(gene_names_df) && nrow(gene_names_df) > 0) {
        gene_names_map <- setNames(gene_names_df$GeneName, gene_names_df$GeneID)
        gene_names <- unname(gene_names_map[genes])
        SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(
            transcript_id = tx_ids,
            gene_id = genes,
            gene_name = gene_names,
            row.names = tx_ids
        )
    } else {
        SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(
            transcript_id = tx_ids,
            gene_id = genes,
            row.names = tx_ids
        )
    }

    # Apply metadata mapping if metadata is provided
    if (!is.null(metadata)) {
        se <- .map_metadata(se, metadata)
    }

    return(se)
}
