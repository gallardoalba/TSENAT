## Helper: Map external coldata into a SummarizedExperiment
#' @keywords internal
#' @noRd
.map_metadata <- function(ts_se, coldata, coldata_sample_col = "Sample", coldata_condition_col = "Condition") {
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
    
    # Prepare canonical condition order and base extraction.  `conds` is sorted
    # for deterministic ordering when constructing paired column order.
    conds <- sort(unique(as.character(coldata[[condition_col_idx]])))
    
    # Automatically detect pairing from third column if present, 
    # or use suffix-removal if paired column is not available
    if (ncol(coldata) >= 3) {
        # Use the third column for sample base identifiers (pairing information)
        coldata_base <- as.character(coldata[[3]])
        # Validate that pairing is consistent (each base has all conditions)
        has_pairing <- TRUE
    } else {
        # Extract base by removing trailing suffix (e.g. _N/_T)
        coldata_base <- sub("_[^_]+$", "", as.character(coldata[[sample_col_idx]]))
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
    sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))

    # Reorder the SummarizedExperiment columns following the order in coldata.
    # Note: For paired analyses, the explicit pairing information is stored in
    # colData(ts_se)$sample_base and used directly by paired statistical tests,
    # so column reordering is not required.
    base_names <- sample_base_names
    idx_list <- integer(0)
    
    # Follow the order present in `coldata`
    ordered_samples <- as.character(coldata[[sample_col_idx]])
    for (s in ordered_samples) {
        matches <- which(base_names == s)
        if (length(matches) > 0) {
            idx_list <- c(idx_list, matches)
        }
    }
    remaining <- setdiff(seq_along(base_names), idx_list)
    new_order <- c(idx_list, remaining)
    
    if (length(new_order) > 0 && !all(new_order == seq_along(base_names))) {
        ts_se <- ts_se[, new_order, drop = FALSE]
        sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))
    }
    
    # Create a mapping from sample names to their pairing information (coldata_base)
    pairing_map <- setNames(coldata_base, as.character(coldata[[sample_col_idx]]))
    sample_pairing <- unname(pairing_map[sample_base_names])
    
    # Create a mapping for actual sample base names (always from Sample column)
    sample_name_map <- setNames(as.character(coldata[[sample_col_idx]]), as.character(coldata[[sample_col_idx]]))
    
    st_map <- setNames(as.character(coldata[[condition_col_idx]]), as.character(coldata[[sample_col_idx]]))
    sample_types <- unname(st_map[sample_base_names])
    missing_idx <- which(is.na(sample_types))
    if (length(missing_idx) > 0) {
        missing_samples <- sample_base_names[missing_idx]
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
        sample_names_full <- sub("_q=.*", "", assay_cols)
        
        # Find original colData row for each assay column's sample
        sample_col_values <- as.character(coldata[[sample_col_idx]])
        
        # Map each sample name to its original row index
        expanded_rows <- integer(n_assay_cols)
        for (i in seq_along(assay_cols)) {
            sample_name <- sample_names_full[i]
            # Find matching row in original coldata
            match_idx <- which(sample_col_values == sample_name)[1]
            if (!is.na(match_idx)) {
                expanded_rows[i] <- match_idx
            } else {
                # If no match found, use first row as fallback
                expanded_rows[i] <- 1
            }
        }
        
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
    sample_names_final <- sub("_q=.*", "", rownames(col_data_final))
    
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
      batch_map <- setNames(as.character(coldata[[4]]), as.character(coldata[[sample_col_idx]]))
      col_data_final$batch <- unname(batch_map[sample_names_final])
    }
    
    SummarizedExperiment::colData(ts_se) <- col_data_final
    # Attach transcript-level readcounts and tx->gene mapping to metadata if
    # they are available in the calling environment or globalenv and not
    # already present in the SummarizedExperiment metadata. This simplifies
    # downstream plotting helpers that expect these objects.
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
        md <- S4Vectors::metadata(ts_se)
        # prefer existing metadata values; otherwise try common names
        if (is.null(md$readcounts)) {
            if (exists("readcounts", envir = parent.frame())) {
                md$readcounts <- get("readcounts", envir = parent.frame())
            } else if (exists("readcounts", envir = globalenv())) {
                md$readcounts <- get("readcounts", envir = globalenv())
            }
        }
        if (is.null(md$tx2gene)) {
            # vignette uses 'txmap' variable name; also accept 'tx2gene'
            if (exists("txmap", envir = parent.frame())) {
                md$tx2gene <- get("txmap", envir = parent.frame())
            } else if (exists("tx2gene", envir = parent.frame())) {
                md$tx2gene <- get("tx2gene", envir = parent.frame())
            } else if (exists("txmap", envir = globalenv())) {
                md$tx2gene <- get("txmap", envir = globalenv())
            } else if (exists("tx2gene", envir = globalenv())) {
                md$tx2gene <- get("tx2gene", envir = globalenv())
            }
        }
        S4Vectors::metadata(ts_se) <- md
    }
    # If a diversity assay is present, prepare a simple diversity data.frame
    # (genes + per-sample diversity values) and store it in metadata so the
    # vignette and plotting helpers can use a ready-made table.
    if ("diversity" %in% SummarizedExperiment::assayNames(ts_se)) {
        div_mat <- as.matrix(SummarizedExperiment::assay(ts_se, "diversity"))
        # sample base names without per-q suffixes
        sample_base_names <- sub("_q=.*", "", colnames(div_mat))
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
map_samples_to_group <- function(sample_names, se = NULL, sample_type_col = NULL,
    mat = NULL) {
    # Prefer explicit mapping from colData(se)[, sample_type_col] when
    # provided. If `sample_type_col` is not provided, allow a single- condition
    # dataset by assigning a single default group 'Group' to all samples (this
    # permits plotting single-condition q-curves).
    
    # Get base names from either mat or se
    if (!is.null(mat)) {
        base_names <- sub("_q=.*", "", colnames(mat))
    } else if (!is.null(se)) {
        base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se)))
    } else {
        # Both are NULL - return default group for all samples
        return(setNames(rep("Group", length(sample_names)), sample_names))
    }

    if (!is.null(se) && !is.null(sample_type_col) && (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st_vec <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st_vec) <- base_names
        st_map <- st_vec[!duplicated(names(st_vec))]
    } else {
        # No explicit mapping: assume single-group dataset
        st_map <- setNames(rep("Group", length(base_names)), base_names)
    }

    mapped <- unname(st_map[sample_names])
    missing_idx <- which(is.na(mapped))
    if (length(missing_idx) > 0) {
        stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(sample_names[missing_idx]),
            collapse = ", ")))
    }
    mapped
}

# Prepare a long-format data.frame for a simple assay (one value per sample)
get_assay_long <- function(se, assay_name = "diversity", value_name = "diversity",
    sample_type_col = NULL) {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("tidyr required")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("dplyr required")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop("Assay not found: ", assay_name)
    }
    df <- as.data.frame(mat)
    genes_col <- .get_gene_ids(se)
    if (is.null(genes_col)) {
        genes_col <- rownames(df)
    }
    df <- cbind(df, Gene = genes_col)
    long <- tidyr::pivot_longer(df, -Gene, names_to = "sample", values_to = value_name)

    # sample_type: prefer explicit colData mapping when available. If not
    # provided, assume a single-group dataset and set `sample_type` to 'Group'
    # for all samples.
    if (!is.null(sample_type_col) && (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st) <- colnames(mat)
        st_map <- st[!duplicated(names(st))]
        sample_base <- sub("_q=.*", "", long$sample)
        long$sample_type <- unname(st_map[sample_base])
        missing_idx <- which(is.na(long$sample_type))
        if (length(missing_idx) > 0) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(sample_base[missing_idx]),
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
prepare_tsallis_long <- function(se, assay_name = "diversity", sample_type_col = "sample_type") {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("tidyr required")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("dplyr required")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop("Assay not found: ", assay_name)
    }
    df <- as.data.frame(mat)
    # Prefer gene names if available in rowData, otherwise use gene IDs, then fall back to rownames
    rd <- SummarizedExperiment::rowData(se)
    genes_col <- if (!is.null(rd) && "gene_name" %in% colnames(rd) && !is.null(rd$gene_name)) {
        rd$gene_name
    } else {
        gene_ids <- .get_gene_ids(se)
        if (is.null(gene_ids)) rownames(df) else gene_ids
    }
    df <- cbind(df, Gene = genes_col)

    long <- tidyr::pivot_longer(df, -Gene, names_to = "sample_q", values_to = "tsallis")
    if (any(grepl("_q=", long$sample_q))) {
        long <- tidyr::separate(long, sample_q, into = c("sample", "q"), sep = "_q=")
        long$q <- as.numeric(long$q)
    } else {
        long$sample <- long$sample_q
        long$q <- NA
    }

    if (!is.null(sample_type_col) && (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        # Get colData which should have been populated by map_metadata()
        col_data <- SummarizedExperiment::colData(se)
        col_st <- as.character(col_data[, sample_type_col])
        col_rownames <- rownames(col_data)
        
        # Create a mapping from unique sample names (without _q=) to sample type
        assay_cols_unique <- unique(sub("_q=.*", "", colnames(mat)))
        st_map <- setNames(rep(NA_character_, length(assay_cols_unique)), assay_cols_unique)
        
        # Strategy 1: If colData has rownames set (by map_metadata), use them to build mapping
        if (!is.null(col_rownames) && length(col_rownames) > 0 && !all(is.na(col_rownames))) {
            # colData rownames should be the full assay column names (with _q= suffixes if present)
            # Extract unique sample names from colData rownames
            col_rownames_unique <- unique(sub("_q=.*", "", col_rownames))
            
            # For each unique sample in colData rownames, find its sample type
            for (sname in col_rownames_unique) {
                # Find first row matching this sample name
                matching_idx <- which(sub("_q=.*", "", col_rownames) == sname)[1]
                if (!is.na(matching_idx)) {
                    st_map[sname] <- col_st[matching_idx]
                }
            }
        } else if (nrow(col_data) == length(assay_cols_unique)) {
            # Strategy 2: Fallback - assume colData rows correspond to unique samples in order
            st_map <- setNames(col_st[1:length(assay_cols_unique)], assay_cols_unique)
        } else if (nrow(col_data) == nrow(mat) || nrow(col_data) == ncol(mat)) {
            # Strategy 3: colData has one row per column in assay
            # Create mapping by extracting unique sample from each column
            for (i in seq_along(assay_cols_unique)) {
                # Find first occurrence of this sample in assay columns
                first_col_idx <- which(sub("_q=.*", "", colnames(mat)) == assay_cols_unique[i])[1]
                if (!is.na(first_col_idx) && first_col_idx <= nrow(col_data)) {
                    st_map[assay_cols_unique[i]] <- col_st[first_col_idx]
                }
            }
        }
        
        # Assign group based on the mapping
        long$group <- unname(st_map[as.character(long$sample)])
        missing_idx <- which(is.na(long$group))
        if (length(missing_idx) > 0) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(as.character(long$sample)[missing_idx]),
                collapse = ", ")))
        }
    } else {
        long$group <- rep("Group", nrow(long))
    }

    as.data.frame(long[!is.na(long$tsallis), , drop = FALSE])
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
#' @keywords internal
#' @noRd
map_tx_to_readcounts <- function(readcounts, tx2gene, tx_col = "Transcript", verbose = FALSE) {
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
