## OPTIMIZATION HELPERS: Internal utilities for mapping functions
# Remove _q=... or _qX.X suffixes from sample/column names
# @noRd
.strip_q_suffix <- function(names) {
    sub("_q=.*", "", names)
}

# Remove _q=... or _qX.X suffixes (handles both formats)
# @noRd
.strip_q_format <- function(names) {
    sub("_q[=0-9].*", "", names)
}

# Create sample-to-condition mapping (OPTIMIZATION: consolidated from 4 duplicate implementations)
# @noRd
.create_sample_type_map <- function(coldata, sample_col_idx, condition_col_idx) {
    setNames(as.character(coldata[[condition_col_idx]]), 
             as.character(coldata[[sample_col_idx]]))
}

# Create pairing/batch mapping
# @noRd
.create_mapping <- function(coldata, col_idx, name_col_idx) {
    setNames(as.character(coldata[[col_idx]]), 
             as.character(coldata[[name_col_idx]]))
}

# Resolve gene names with fallback logic (OPTIMIZATION: consolidated from 6-step process)
# @noRd
.resolve_gene_names <- function(se, df_nrows) {
    # Step 1: Try rowData gene_name
    rd <- SummarizedExperiment::rowData(se)
    if (!is.null(rd) && nrow(rd) > 0 && "gene_name" %in% colnames(rd)) {
        return(rd$gene_name)
    }
    
    # Step 2: Try gene_ids from metadata
    gene_ids <- .get_gene_ids(se)
    if (!is.null(gene_ids) && length(gene_ids) == df_nrows) {
        return(gene_ids)
    }
    
    # Step 3: Try rownames
    rn <- rownames(se)
    if (!is.null(rn) && length(rn) == df_nrows && !all(is.na(rn))) {
        return(rn)
    }
    
    # Step 4: Generate fallback names
    paste0("gene_", seq_len(df_nrows))
}

# Check required packages and stop if missing
# @noRd
.check_required_packages <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(pkg, " required")
        }
    }
    invisible(TRUE)
}

## Helper: Map external coldata into a SummarizedExperiment

#' @noRd
.map_metadata_se <- function(ts_se, coldata, coldata_sample_col = "Sample", coldata_condition_col = "Condition") {
    if (is.null(coldata)) {
        return(ts_se)
    }
    if (!is.data.frame(coldata)) {
        return(ts_se)
    }
    
    # Position-based column detection (PRIMARY STRATEGY)
    # This is more robust for fixed-structure metadata files
    sample_col_idx <- 1
    condition_col_idx <- 2
    
    # Verify minimum columns required
    if (ncol(coldata) < 2) {
        return(ts_se)
    }
    
    # Alternative: named column detection (FALLBACK - only if not found by position)
    # This allows flexibility if columns are in different order
    # IMPROVED: Try case-insensitive matching first
    named_sample_col <- which(tolower(colnames(coldata)) == tolower(coldata_sample_col))
    named_condition_col <- which(tolower(colnames(coldata)) == tolower(coldata_condition_col))
    if (length(named_sample_col) > 0) {
        sample_col_idx <- named_sample_col[1]
    }
    if (length(named_condition_col) > 0) {
        condition_col_idx <- named_condition_col[1]
    }
    
    # OPTIMIZATION: Cache character conversions to avoid repeated as.character() calls
    coldata_sample_col_values <- as.character(coldata[[sample_col_idx]])
    coldata_condition_col_values <- as.character(coldata[[condition_col_idx]])
    
    # Prepare canonical condition order and base extraction.  `conds` is sorted
    # for deterministic ordering when constructing paired column order.
    conds <- sort(unique(coldata_condition_col_values))
    
    # Automatically detect pairing from third column if present, 
    # or use suffix-removal if paired column is not available
    if (ncol(coldata) >= 3) {
        # Use the third column for sample base identifiers (pairing information)
        coldata_base <- as.character(coldata[[3]])
        # Validate that pairing is consistent (each base has all conditions)
        has_pairing <- TRUE
    } else {
        # Extract base by removing trailing suffix (e.g. _N/_T)
        coldata_base <- sub("_[^_]+$", "", coldata_sample_col_values)
        has_pairing <- FALSE
    }
    bases <- unique(coldata_base)

    # Validate pairing structure if pairing information is present
    if (has_pairing && length(conds) >= 2) {
        unpaired <- vapply(bases, function(b) {
            length(unique(coldata[[condition_col_idx]][coldata_base == b]))
        }, integer(1))
        bad <- bases[unpaired != length(conds)]
        if (length(bad) > 0) {
            bad_list <- paste(bad, collapse = ", ")
            cond_list <- paste(conds, collapse = ", ")
            msg <- paste0("Unpaired samples found in coldata for bases: ", bad_list,
                ". Ensure each base has all conditions: ", cond_list)
            stop(msg, call. = FALSE)
        }
    }
    # OPTIMIZATION: Use .strip_q_suffix() helper
    sample_base_names <- .strip_q_suffix(colnames(SummarizedExperiment::assay(ts_se)))

    # Reorder the SummarizedExperiment columns following the order in coldata.
    # Note: For paired analyses, the explicit pairing information is stored in
    # colData(ts_se)$sample_base and used directly by paired statistical tests,
    # so column reordering is not required.
    base_names <- sample_base_names
    idx_list <- integer(0)
    
    # Follow the order present in `coldata` (groups multiple q-values per sample together)
    ordered_samples <- coldata_sample_col_values
    for (s in ordered_samples) {
        # Find ALL columns for this sample (important when there are multiple q-values)
        matches <- which(base_names == s)
        if (length(matches) > 0) {
            idx_list <- c(idx_list, matches)
        }
    }
    remaining <- setdiff(seq_along(base_names), idx_list)
    new_order <- c(idx_list, remaining)
    
    if (length(new_order) > 0 && !all(new_order == seq_along(base_names))) {
        ts_se <- ts_se[, new_order, drop = FALSE]
        # OPTIMIZATION: Use .strip_q_suffix() helper for consistency and caching benefit
        sample_base_names <- .strip_q_suffix(colnames(SummarizedExperiment::assay(ts_se)))
    }
    
    # Create a mapping from sample names to their pairing information (coldata_base)
    if (ncol(coldata) < 3) {
        # No explicit pairing column - use coldata_base extracted from sample names
        pairing_map <- setNames(coldata_base, coldata_sample_col_values)
    } else {
        # Use explicit pairing column (column 3)
        pairing_map <- .create_mapping(coldata, 3, sample_col_idx)
    }
    sample_pairing <- unname(pairing_map[sample_base_names])
    
    # Create a mapping for actual sample base names (always from Sample column)
    sample_name_map <- setNames(coldata_sample_col_values, coldata_sample_col_values)
    
    # OPTIMIZATION: Use .create_sample_type_map() helper
    st_map <- .create_sample_type_map(coldata, sample_col_idx, condition_col_idx)
    sample_types <- unname(st_map[sample_base_names])
    has_na <- is.na(sample_types)
    if (any(has_na)) {
        missing_samples <- sample_base_names[has_na]
        msg <- paste0("map_metadata: unmatched samples in 'coldata': ", paste(missing_samples,
            collapse = ", "), ". Provide matching entries in 'coldata' or populate ",
            "colData(ts_se)$sample_type beforehand.")
        stop(msg, call. = FALSE)
    }
    # Don't set sample_type yet - do it AFTER colData expansion
    # SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
    
    # Record the pairing information (from the third column of coldata) per column
    # so downstream functions can use explicit pairing for paired tests.
    # Also don't set sample_base yet - do it AFTER colData expansion
    # SummarizedExperiment::colData(ts_se)$sample_base <- sample_pairing
    
    # Expand colData if assay has more columns than colData rows
    # (e.g., after calculate_diversity adds _q= suffixes)
    n_coldata_rows <- nrow(SummarizedExperiment::colData(ts_se))
    assay_cols <- colnames(SummarizedExperiment::assay(ts_se))
    n_assay_cols <- length(assay_cols)
    
    if (n_coldata_rows < n_assay_cols) {
        # Assay has been expanded (multiple q-values per sample)
        # Expand colData to match by repeating rows for each q-value
        col_data <- SummarizedExperiment::colData(ts_se)
        sample_names_full <- .strip_q_suffix(assay_cols)
        
        # OPTIMIZATION: Use vectorized match() instead of loop (O(n) vs O(n²))
        # Match against original coldata sample column values, not reordered sample_base_names
        expanded_rows <- match(sample_names_full, coldata_sample_col_values)
        # Replace NA with 1 (fallback to first row)
        expanded_rows[is.na(expanded_rows)] <- 1
        
        # Expand colData using the mapped indices
        new_col_data <- col_data[expanded_rows, ]
        rownames(new_col_data) <- assay_cols
        SummarizedExperiment::colData(ts_se) <- new_col_data
    } else if (n_coldata_rows == n_assay_cols) {
        # Dimensions already match, just set rownames
        rownames(SummarizedExperiment::colData(ts_se)) <- assay_cols
    }
    
    # NOW set sample_type and pairing column after colData expansion is complete
    # This ensures all rows have these values properly assigned
    col_data_final <- SummarizedExperiment::colData(ts_se)
    # OPTIMIZATION: Use .strip_q_suffix() helper
    sample_names_final <- .strip_q_suffix(rownames(col_data_final))
    
    # Map each expanded row's sample name to its condition and pairing
    col_data_final$sample_type <- unname(st_map[sample_names_final])
    
    # Always create standardized sample_base column with pairing identifiers
    # This contains the pairing information (A, B, C, etc.) used by paired tests
    col_data_final$sample_base <- unname(pairing_map[sample_names_final])
    
    # Preserve the actual pairing column name from metadata (e.g., paired_samples)
    if (ncol(coldata) >= 3 && has_pairing) {
      paired_col_name <- colnames(coldata)[3]
      col_data_final[[paired_col_name]] <- unname(pairing_map[sample_names_final])
    }
    
    # Map batch column if present (typically column 4)
    if (ncol(coldata) >= 4) {
      batch_col_name <- colnames(coldata)[4]
      # OPTIMIZATION: Use cached sample column values
      batch_map <- setNames(as.character(coldata[[4]]), coldata_sample_col_values)
      col_data_final$batch <- unname(batch_map[sample_names_final])
    }
    
    SummarizedExperiment::colData(ts_se) <- col_data_final
    # Note: readcounts and tx2gene mapping should be explicitly provided via
    # function parameters or stored in the TSENATAnalysis @config slot.
    # We no longer look in parent environment or globalenv() to ensure
    # reproducible, self-contained analysis workflows.
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
        md <- S4Vectors::metadata(ts_se)
        # Only preserve metadata that was explicitly provided
        # Do not attempt to fetch from calling environment
        S4Vectors::metadata(ts_se) <- md
    }
    # If a diversity assay is present, prepare a simple diversity data.frame
    # (genes + per-sample diversity values) and store it in metadata so the
    # vignette and plotting helpers can use a ready-made table.
    if ("diversity" %in% SummarizedExperiment::assayNames(ts_se)) {
        div_mat <- as.matrix(SummarizedExperiment::assay(ts_se, "diversity"))
        # OPTIMIZATION: Use .strip_q_suffix() helper
        sample_base_names <- .strip_q_suffix(colnames(div_mat))
        # prefer explicit sample_type in colData when present
        samples_vec <- NULL
        if ("sample_type" %in% colnames(SummarizedExperiment::colData(ts_se))) {
            samples_vec <- as.character(SummarizedExperiment::colData(ts_se)$sample_type)
        }
        div_df <- as.data.frame(div_mat)
        genes_col <- .get_gene_ids(ts_se)
        if (is.null(genes_col)) {
            genes_col <- rownames(div_df)
        }
        div_df <- cbind(genes = genes_col, div_df)
        md2 <- S4Vectors::metadata(ts_se)
        md2$diversity_df <- div_df
        md2$sample_base_names <- sample_base_names
        if (!is.null(samples_vec)) {
            md2$samples <- samples_vec
        }
        S4Vectors::metadata(ts_se) <- md2
    }
    return(ts_se)
}

# Map sample names (without '_q=...') to group labels using `colData(se)`.
# Mapping must be provided via `colData(se)`; no inference fallback is used.

.map_samples_to_group <- function(sample_names, se = NULL, condition_col = NULL,
    mat = NULL) {
    # Prefer explicit mapping from colData(se)[, condition_col] when
    # provided. If `condition_col` is not provided, allow a single- condition
    # dataset by assigning a single default group 'Group' to all samples (this
    # permits plotting single-condition q-curves).
    
    # OPTIMIZATION: Use .strip_q_suffix() helper
    # Get base names from either mat or se
    if (!is.null(mat)) {
        base_names <- .strip_q_suffix(colnames(mat))
    } else if (!is.null(se)) {
        base_names <- .strip_q_suffix(colnames(SummarizedExperiment::assay(se)))
    } else {
        # Both are NULL - return default group for all samples
        return(setNames(rep("Group", length(sample_names)), sample_names))
    }

    if (!is.null(se) && !is.null(condition_col) && (condition_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st_vec <- as.character(SummarizedExperiment::colData(se)[, condition_col])
        names(st_vec) <- base_names
        st_map <- st_vec[!duplicated(names(st_vec))]
    } else {
        # No explicit mapping: assume single-group dataset
        st_map <- setNames(rep("Group", length(base_names)), base_names)
    }

    mapped <- unname(st_map[sample_names])
    has_na <- is.na(mapped)
    if (any(has_na)) {
        stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(sample_names[has_na]),
            collapse = ", ")))
    }
    mapped
}

# Prepare a long-format data.frame for a simple assay (one value per sample)

.get_assay_long <- function(se, assay_name = "diversity", value_name = "diversity",
    condition_col = NULL) {
    # OPTIMIZATION: Use .check_required_packages() helper
    .check_required_packages(c("tidyr", "dplyr", "SummarizedExperiment"))

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop("Assay not found: ", assay_name)
    }
    df <- as.data.frame(mat)
    # OPTIMIZATION: Use .resolve_gene_names() helper
    genes_col <- .resolve_gene_names(se, nrow(df))
    df <- data.frame(Gene = genes_col, df, row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
    long <- tidyr::pivot_longer(df, -Gene, names_to = "sample", values_to = value_name)

    # sample_type: prefer explicit colData mapping when available. If not
    # provided, assume a single-group dataset and set `sample_type` to 'Group'
    # for all samples.
    if (!is.null(condition_col) && (condition_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st <- as.character(SummarizedExperiment::colData(se)[, condition_col])
        names(st) <- colnames(mat)
        st_map <- st[!duplicated(names(st))]
        # OPTIMIZATION: Use .strip_q_suffix() helper
        sample_base <- .strip_q_suffix(long$sample)
        long$sample_type <- unname(st_map[sample_base])
        has_na <- is.na(long$sample_type)
        if (any(has_na)) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(sample_base[has_na]),
                collapse = ", ")))
        }
    } else {
        long$sample_type <- rep("Group", nrow(long))
    }

    # Filter out NA values but retain sample_type information
    long_filtered <- long[!is.na(long[[value_name]]), , drop = FALSE]

    # Check if any data remains after filtering
    if (nrow(long_filtered) == 0) {
        stop(sprintf("No non-NA values found in assay '%s'. All values are NA.",
            assay_name), call. = FALSE)
    }

    long_filtered
}

# Internal small helper: prepare long-format tsallis data from a
# SummarizedExperiment

.prepare_tsallis_long <- function(se, assay_name = "diversity", condition_col = "sample_type") {
    # OPTIMIZATION: Use .check_required_packages() helper
    .check_required_packages(c("tidyr", "dplyr", "SummarizedExperiment"))

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop("Assay not found: ", assay_name)
    }
    
    # OPTIMIZATION: Use .resolve_gene_names() helper (consolidates 6 fallback steps)
    genes_col <- .resolve_gene_names(se, nrow(mat))
    genes_col <- as.character(genes_col)
    
    # Create data frame with Gene column as first column
    # Important: Create df first with rownames, then add Gene column explicitly
    df <- as.data.frame(mat)
    df <- cbind(Gene = genes_col, df, stringsAsFactors = FALSE)

    long <- tidyr::pivot_longer(df, cols = -Gene, names_to = "sample_q", values_to = "tsallis")
    
    # Handle different column naming conventions (OPTIMIZATION: Unified regex-based parser)
    # Convention 1: Old format "sample_q=X.X" (from calculate_diversity multi-q)
    # Convention 2: New format "sample_qX.X" (from TSENATAnalysis S4 wrapper)
    # Normalize old format "_q=" to "_q" for consistent parsing
    sample_q_col <- gsub("_q=", "_q", long$sample_q)
    
    # Extract parts using regex subexpression matching: matches "sample_qX.X" pattern
    parsed <- regmatches(sample_q_col, regexec("^(.+)_q([0-9.]+)$", sample_q_col))
    
    # Extract sample and q from parsed results
    long$sample <- vapply(parsed, function(x) if (length(x) > 1) x[2] else NA_character_, character(1))
    q_values <- vapply(parsed, function(x) if (length(x) > 2) x[3] else NA_character_, character(1))
    # Convert only non-NA q_values to numeric to avoid unnecessary warnings
    long$q <- NA_real_
    valid_idx <- !is.na(q_values)
    long$q[valid_idx] <- as.numeric(q_values[valid_idx])
    
    # For entries without q-value pattern, use sample_q as sample name
    no_q_match <- is.na(long$sample)
    long$sample[no_q_match] <- long$sample_q[no_q_match]
    
    # Remove entries with non-numeric q values that matched the pattern
    has_q_pattern <- grepl("_q[0-9.]", sample_q_col)
    invalid_q <- has_q_pattern & is.na(long$q)
    if (any(invalid_q)) {
        long <- long[!invalid_q, ]
    }

    if (!is.null(condition_col) && (condition_col %in% colnames(SummarizedExperiment::colData(se)))) {
        # Get colData which should have been populated by map_metadata()
        col_data <- SummarizedExperiment::colData(se)
        col_st <- as.character(col_data[, condition_col])
        col_rownames <- rownames(col_data)
        
        # Create a mapping from unique sample names to sample type
        # Handle both "_q=" format (old) and "_qX.X" format (new)
        assay_cols_unique <- unique(.strip_q_format(colnames(mat)))
        st_map <- setNames(rep(NA_character_, length(assay_cols_unique)), assay_cols_unique)
        
        # Strategy 1: If colData has rownames set, use them to build mapping
        if (!is.null(col_rownames) && length(col_rownames) > 0 && !all(is.na(col_rownames))) {
            # OPTIMIZATION: Cache .strip_q_format() result - avoid repeated regex evaluation
            # colData rownames should be the full assay column names (with _q suffixes)
            col_rownames_stripped <- .strip_q_format(col_rownames)
            col_rownames_unique <- unique(col_rownames_stripped)
            
            # For each unique sample, find first matching row and get its sample type
            # Vectorized approach: get first occurrence of each unique sample
            first_occurrences <- match(col_rownames_unique, col_rownames_stripped)
            valid_matches <- !is.na(first_occurrences)
            st_map <- setNames(col_st[first_occurrences[valid_matches]], col_rownames_unique[valid_matches])
        } else if (nrow(col_data) == length(assay_cols_unique)) {
            # Strategy 2: Fallback - assume colData rows correspond to unique samples in order
            st_map <- setNames(col_st[seq_along(assay_cols_unique)], assay_cols_unique)
        } else if (nrow(col_data) == nrow(mat) || nrow(col_data) == ncol(mat)) {
            # Strategy 3: colData has one row per column in assay
            # OPTIMIZATION: Cache .strip_q_format() result - avoid repeated regex evaluation in loop
            mat_cols_stripped <- .strip_q_format(colnames(mat))
            # Get first occurrence of each unique sample in assay columns
            first_occurrences_mat <- match(assay_cols_unique, mat_cols_stripped)
            valid_matches_mat <- !is.na(first_occurrences_mat) & first_occurrences_mat <= nrow(col_data)
            st_map <- setNames(col_st[first_occurrences_mat[valid_matches_mat]], assay_cols_unique[valid_matches_mat])
        }
        
        # Assign group based on the mapping
        long$group <- unname(st_map[as.character(long$sample)])
        has_na <- is.na(long$group)
        if (any(has_na)) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(as.character(long$sample)[has_na]),
                collapse = ", ")))
        }
    } else {
        long$group <- rep("Group", nrow(long))
    }

    long_filtered <- as.data.frame(long[!is.na(long$tsallis), , drop = FALSE])
    
    # Check if any data remains after filtering
    if (nrow(long_filtered) == 0) {
        stop(sprintf("No tsallis values found in SummarizedExperiment. Check that assay '%s' contains valid data.",
            assay_name), call. = FALSE)
    }
    
    long_filtered
}


#' Map transcript IDs from a tx2gene table to a readcounts matrix
#'
#' Assign transcript identifiers as row names of a transcript-level read
#' counts matrix so downstream plotting functions can identify transcripts.
#' The tx2gene mapping can be provided as a file path (TSV) or a
#' data.frame with a `Transcript` column.
#'
#' @param readcounts A numeric matrix or data.frame of read counts (rows = transcripts).
#' @param tx2gene Either a path to a tab-delimited file or a data.frame with at least a `Transcript` column.
#' @param tx_col Name of the transcript ID column in `tx2gene` (default: 'Transcript').
#' @param verbose Logical; print informative messages (default: FALSE).
#' @return The input `readcounts` with rownames set to the transcript IDs.

#' @noRd

.map_tx_to_readcounts <- function(readcounts, tx2gene, tx_col = "Transcript", verbose = FALSE) {
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }
        txmap <- utils::read.delim(tx2gene, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(tx2gene)) {
        txmap <- tx2gene
    } else {
        stop("'tx2gene' must be a file path or a data.frame.", call. = FALSE)
    }

    if (!(tx_col %in% colnames(txmap))) {
        stop(sprintf("tx2gene mapping must contain column '%s'", tx_col), call. = FALSE)
    }

    # Ensure readcounts is a matrix-like object
    if (is.data.frame(readcounts)) {
        readcounts <- as.matrix(readcounts)
    }
    if (!is.matrix(readcounts)) {
        stop("'readcounts' must be a matrix or data.frame", call. = FALSE)
    }

    n_rc <- nrow(readcounts)
    n_tx <- nrow(txmap)

    if (n_tx == n_rc) {
        rownames(readcounts) <- as.character(txmap[[tx_col]])
        if (verbose) {
            message(sprintf("Assigned %d transcript rownames from tx2gene.", n_rc))
        }
        return(readcounts)
    }

    # If counts differ, attempt to match by transcript identifiers if present
    tx_ids <- as.character(txmap[[tx_col]])
    if (!is.null(rownames(readcounts)) && all(rownames(readcounts) %in% tx_ids)) {
        # reorder txmap to match readcounts row order and set rownames
        matched_idx <- match(rownames(readcounts), tx_ids)
        rownames(readcounts) <- tx_ids[matched_idx]
        if (verbose) {
            message("Matched and assigned transcript IDs by existing readcounts rownames.")
        }
        return(readcounts)
    }

    stop(sprintf("Number of transcripts in tx2gene (%d) does not match readcounts rows (%d), and automatic matching failed.",
        n_tx, n_rc), call. = FALSE)
}
