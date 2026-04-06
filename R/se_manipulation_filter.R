## Filter SummarizedExperiment by low-expression transcripts
#' Filter transcripts in a `SummarizedExperiment` by minimum TPM/sample
#'
#' Subset a `SummarizedExperiment` to keep transcripts with more than
#' `min_tpm` TPM (transcripts per million) threshold in strictly more than
#' `min_samples` samples. The function updates assays, `rowData`, and
#' relevant entries in `metadata()` (for example `readcounts` and `tx2gene`)
#' so downstream helpers receive a consistent object.
#'
#' **Filtering Strategy**: Uses TPM-based filtering for results comparable
#' across
#' studies (data from SALMON quantification is already TPM-normalized;
#' Soneson et al. 2015, Law et al. 2014).
#'
#' **TPM Data Source**: Function automatically locates TPM data via:
#' 1. `tpm_assay_name` parameter (specify assay name containing TPM values)
#' 2. Metadata: looks for `tpm` or `tpm` in SummarizedExperiment metadata
#' 3. Assay names: searches for assay named 'tpm' or 'abundance'
#' If TPM data NOT found, filtering parameters are compared against raw counts
#' but this is NOT recommended (will produce incorrect results).
#'
#' @param se A `SummarizedExperiment` object containing transcript-level
#'   assays (rows = transcripts, columns = samples). Can contain:
#'   - SALMON TPM values (recommended, already normalized)
#'   - SALMON raw counts (NumReads) if TPM not available
#' @param min_tpm Numeric TPM threshold (default 1.0).
#'   Keeps transcripts with TPM >= `min_tpm` in >= `min_samples` samples.
#' TPM >= 1 is recommended for typical sequencing depth (Soneson et al.
#' 2015, Law et al. 2014).
#' Ignored if `stringency` is specified; when stringency is provided,
#' `min_tpm` is
#' auto-estimated from data using quantile-based approach (Law et al.
#' limma-voom methodology).
#' @param tpm_assay_name Character; name of assay containing TPM data
#' (default: NULL).
#' If NULL, function searches for TPM data in: metadata$tpm ->
#' metadata$tpm -> assay named 'tpm' -> metadata lookup.
#' Set explicitly (e.g., `tpm_assay_name = 'tpm'`) to use a specific assay
#' by name.
#' Example from preprocessing: tpm_assay_name = 'abundance' for tximport
#' objects.
#' @param min_samples Integer minimum number of samples exceeding
#'   the threshold required to keep a transcript (default 5L). Ignored if
#'   `stringency` is specified.
#' @param stringency Character; auto-calculate `min_samples` and `min_tpm`
#' based on
#'   paired design stringency. Options:
#'   - 'soft' (permissive): 25% of samples (min 2), min_tpm = Q1 (25th %ile)
#'   - 'medium' (balanced): 50% of samples (min 3), min_tpm = Q2 (median)
#'   - 'severe' (stringent): 75% of samples, min_tpm = Q3 (75th %ile)
#'   - NULL (default): use explicit `min_samples` and `min_tpm`
#' When stringency is specified, all three filtering parameters are
#' auto-estimated
#'   from the data distribution (Law et al. limma-voom methodology).
#'   Requires `pair_col` in colData.
#' @param pair_col Character; column name in colData containing pair IDs.
#'   Required when `stringency` is specified. If NULL and `stringency` is set,
#'   function attempts to auto-detect from colData.
#' @param min_tx_per_gene Integer minimum number of transcripts per gene
#'   required to keep a gene (default 2L). Genes with fewer transcripts are
#'   removed since entropy is always 0 for single-transcript genes.
#'   Ignored if `stringency` is specified; when stringency is provided, this is
#' auto-adjusted (soft=2, medium=2, severe=3) to ensure isoform diversity
#' matches
#'   filtering stringency.
#' @param min_isoform_abundance Numeric in [0, 1]; minimum relative
#' abundance threshold
#' for isoforms within each gene (default: 0.05 for 5%). Implements Soneson
#' et al. (2016)
#' isoform-level filtering to improve FDR control. Removes isoforms with
#' mean relative
#' abundance below threshold within their gene. Single-isoform genes are
#' always kept.
#'   Set to 0 or NULL to skip isoform-level filtering.
#' @param assay_name Name or index of the assay to use for filtering
#'   (DEPRECATED: use `tpm_assay_name` instead, default: 'counts'). 
#' **WARNING**: If set to 'counts', filtering uses raw counts with TPM
#' thresholds,
#'   which is incorrect. Use `tpm_assay_name` to specify TPM assay.
#' @param verbose Logical; print before/after counts and filtering
#' parameters when TRUE.
#' @return A filtered `SummarizedExperiment`.
#' @examples
#' mat <- matrix(c(0, 6, 7, 2, 8, 9), nrow = 3, dimnames = list(paste0('tx',
#' 1:3), paste0('S', 1:2)))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
#' filt <- .filter_se(se, min_samples = 1)
#' class(filt)
#' @noRd
.filter_se <- function(se, min_samples = 5L, stringency = NULL, pair_col = NULL,
    min_tpm = 1, tpm_assay_name = NULL, min_tx_per_gene = 2L, min_isoform_abundance = NULL,
    assay_name = "counts", verbose = TRUE) {
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment", call. = FALSE)
    }

    # Phase 1: Consolidate parameter resolution (manual vs stringency-based)
    params <- .resolve_filter_parameters(se, min_samples, min_tpm, stringency, pair_col,
        tpm_assay_name, assay_name, min_tx_per_gene, min_isoform_abundance, verbose = verbose)

    # Set default min_isoform_abundance if still NULL
    if (is.null(params$min_isoform_abundance)) {
        params$min_isoform_abundance <- 0.05
    }

    # Validate assay is numeric
    if (!is.numeric(params$assay_mat)) {
        stop("Assay data must be numeric.", call. = FALSE)
    }

    # Phase 2: Apply combined filtering pipeline
    filter_result <- .apply_combined_filters(params$assay_mat, params$genes_vec,
        list(min_tpm = params$min_tpm, min_samples = params$min_samples, min_tx_per_gene = params$min_tx_per_gene,
            min_isoform_abundance = params$min_isoform_abundance), verbose = verbose)

    tokeep <- filter_result$tokeep

    # Verbose reporting
    if (verbose) {
        message(sprintf("Transcripts: before = %d, after = %d", filter_result$before,
            filter_result$after))
        message(sprintf("Filtering data source: %s", params$assay_source))
    }

    if (filter_result$after == 0L) {
        warning("Filtering removed all transcripts; returning empty SummarizedExperiment.",
            call. = FALSE)
    }

    # Phase 3: Finalize filtered SE with synchronized metadata
    new_se <- .finalize_filtered_se(se, tokeep, SummarizedExperiment::assays(se),
        params$genes_vec, params$min_samples, params$min_tpm, params$min_tx_per_gene,
        params$min_isoform_abundance, stringency)

    return(new_se)
}

# ============================================================================
# HELPER: Resolve assay by name/index with fallback logic OPTIMIZATION:
# Consolidated from 3 duplicate implementations @noRd
.resolve_assay_index <- function(assay_ref, assay_names) {
    if (is.character(assay_ref)) {
        idx <- which(assay_names == assay_ref)
        if (length(idx) == 1)
            return(idx)
        return(1)  # Fallback to first
    } else if (is.numeric(assay_ref)) {
        if (assay_ref >= 1 && assay_ref <= length(assay_names))
            return(assay_ref)
        return(1)  # Fallback to first
    }
    1  # Default to first
}

# ============================================================================
# HELPER: Get assay matrix with TPM priority logic (OPTIMIZATION: single call
# point) @noRd
.get_assay_filtering <- function(se, assays_list, tpm_assay_name, assay_name) {
    # Priority 1: explicit TPM assay
    if (!is.null(tpm_assay_name) && tpm_assay_name %in% names(assays_list)) {
        return(list(mat = as.matrix(assays_list[[tpm_assay_name]]), source = sprintf("assay '%s' (user-specified)",
            tpm_assay_name)))
    }

    # Priority 2: metadata
    md <- S4Vectors::metadata(se)
    if (!is.null(md$tpm) && is.matrix(md$tpm)) {
        return(list(mat = as.matrix(md$tpm), source = "metadata$tpm (SALMON preprocessed)"))
    }
    if (!is.null(md$tpm) && is.matrix(md$tpm)) {
        return(list(mat = as.matrix(md$tpm), source = "metadata$tpm"))
    }

    # Priority 3: auto-detect by name
    if ("tpm" %in% names(assays_list)) {
        return(list(mat = as.matrix(assays_list[["tpm"]]), source = "assay 'tpm' (auto-detected)"))
    }
    if ("abundance" %in% names(assays_list)) {
        return(list(mat = as.matrix(assays_list[["abundance"]]), source = "assay 'abundance' (tximport format)"))
    }

    # Fallback with warning
    idx <- .resolve_assay_index(assay_name, names(assays_list))
    return(list(mat = as.matrix(assays_list[[idx]]), source = sprintf("assay '%s' (fallback - NOT TPM!)",
        names(assays_list)[idx]), is_fallback = TRUE))
}

# ============================================================================
# HELPER: Get gene IDs from rowData or tx2gene mapping Returns NULL if no gene
# mapping available @noRd
.get_gene_ids <- function(se) {
    # Try rowData first
    rd <- SummarizedExperiment::rowData(se)
    if (!is.null(rd) && nrow(rd) > 0 && "gene_id" %in% colnames(rd)) {
        return(rd$gene_id)
    }

    # Try tx2gene mapping from metadata
    md <- S4Vectors::metadata(se)
    if (!is.null(md$tx2gene) && is.data.frame(md$tx2gene) && nrow(md$tx2gene) > 0) {
        tx2gene <- md$tx2gene
        # Assume first column is transcript ID, second is gene ID
        if (ncol(tx2gene) >= 2) {
            txcol <- colnames(tx2gene)[1]
            genecol <- colnames(tx2gene)[2]
            # Match rownames to transcript column and return gene IDs in
            # correct order
            tx_names <- rownames(se)
            if (length(tx_names) > 0) {
                match_idx <- match(tx_names, tx2gene[[txcol]])
                matched_genes <- tx2gene[[genecol]][match_idx]
                if (any(!is.na(matched_genes))) {
                  return(matched_genes)
                }
            }
        }
    }

    NULL
}

# ============================================================================
# HELPER: Validate filter parameters for type and range Raises errors for
# invalid inputs @noRd
.validate_filter_params <- function(min_isoform_abundance, tpm_assay_name) {
    if (!is.null(min_isoform_abundance)) {
        if (!is.numeric(min_isoform_abundance) || length(min_isoform_abundance) !=
            1) {
            stop(".filter_se: 'min_isoform_abundance' must be numeric in [0, 1], or NULL to disable",
                call. = FALSE)
        }
        if (min_isoform_abundance < 0 || min_isoform_abundance > 1) {
            stop(".filter_se: 'min_isoform_abundance' must be numeric in [0, 1], or NULL to disable",
                call. = FALSE)
        }
    }
    invisible(TRUE)
}

# ============================================================================
# HELPER: Auto-estimate min_tpm from data distribution using quantiles
# Implements Law et al. limma-voom methodology @noRd
.estimate_min_tpm <- function(assay_mat, stringency, verbose = FALSE) {
    # Calculate mean TPM per gene for quantile estimation
    mean_tpm <- rowMeans(assay_mat)
    mean_tpm_nonzero <- mean_tpm[mean_tpm > 0]

    if (length(mean_tpm_nonzero) == 0) {
        warning("No non-zero values in assay. Using default min_tpm = 0.1", call. = FALSE)
        return(0.1)
    }

    # Select quantile based on stringency level
    if (stringency == "soft") {
        quantile_prob <- 0.25
        quant_label <- "Q1"
    } else if (stringency == "medium") {
        quantile_prob <- 0.5
        quant_label <- "Q2/Median"
    } else if (stringency == "severe") {
        quantile_prob <- 0.75
        quant_label <- "Q3"
    } else {
        return(NULL)
    }

    min_tpm_estimated <- as.numeric(quantile(mean_tpm_nonzero, probs = quantile_prob,
        na.rm = TRUE))
    # Bound between 0.1 and 5.0 for practical constraints
    min_tpm_estimated <- max(0.1, min(min_tpm_estimated, 5))

    if (verbose) {
        message(sprintf("Auto-estimated min_tpm from data: %.3f (%s of mean TPM distribution)",
            min_tpm_estimated, quant_label))
    }

    return(list(min_tpm = min_tpm_estimated, quant_label = quant_label))
}

# ============================================================================
# HELPER: Calculate min_samples and min_tx_per_gene based on stringency @noRd
.calculate_stringency_thresholds <- function(stringency, n_samples, n_pairs = NULL) {
    if (stringency == "soft") {
        min_samples <- max(2L, ceiling(0.25 * n_samples))
        min_tx_per_gene <- 2L
        min_isoform_abundance <- 0.01  # 1% - permissive
    } else if (stringency == "medium") {
        min_samples <- max(3L, ceiling(0.5 * n_samples))
        min_tx_per_gene <- 2L
        min_isoform_abundance <- 0.05  # 5% - balanced (default)
    } else if (stringency == "severe") {
        min_samples <- ceiling(0.75 * n_samples)
        min_tx_per_gene <- 3L
        min_isoform_abundance <- 0.15  # 15% - stringent
    } else {
        return(NULL)
    }

    list(min_samples = min_samples, min_tx_per_gene = min_tx_per_gene, min_isoform_abundance = min_isoform_abundance)
}

# ============================================================================
# HELPER: Detect pair column from colData Tries multiple candidate names for
# pairing information @noRd
.detect_pair_column <- function(col_data, se) {
    # Candidate column names to try (in order of priority)
    candidates <- c("pair", "pair_id", "pair_num", "pair_number", "subject", "subject_id",
        "replicate_id", "rep_id")

    for (cand in candidates) {
        if (cand %in% colnames(col_data)) {
            return(cand)
        }
    }

    # If no candidate found, throw error
    candidates_str <- paste(candidates, collapse = ", ")
    stop("Could not auto-detect pair column in colData. ", "Please provide pair_col argument explicitly. ",
        "Tried: ", candidates_str, call. = FALSE)
}

# ============================================================================
# HELPER: Calculate stringency parameters from SummarizedExperiment Wrapper
# around .calculate_stringency_thresholds that works with SE objects @noRd
.calc_stringency_params <- function(se, stringency, pair_col = NULL, verbose = FALSE) {
    col_data <- SummarizedExperiment::colData(se)
    n_samples <- ncol(se)

    # Auto-detect pair column if not provided
    if (is.null(pair_col)) {
        pair_col <- .detect_pair_column(col_data, se)
    }

    # Count number of unique pairs
    if (pair_col %in% colnames(col_data)) {
        n_pairs <- length(unique(col_data[[pair_col]]))
    } else {
        n_pairs <- NULL
    }

    # Call existing stringency calculation
    params <- .calculate_stringency_thresholds(stringency, n_samples, n_pairs)

    if (verbose && !is.null(params)) {
        pairs_str <- if (!is.null(n_pairs))
            paste0(" n_pairs=", n_pairs) else ""
        message("[calc_stringency_params] stringency='", stringency, "' n_samples=",
            n_samples, pairs_str, " -> min_samples=", params$min_samples, " min_tx_per_gene=",
            params$min_tx_per_gene)
    }

    return(params)
}
# Returns logical vector indicating which rows to keep @noRd
.apply_tpm_filter <- function(assay_mat, min_tpm, min_samples, verbose = FALSE) {
    # Check if min_samples would result in impossible threshold
    if (min_samples > ncol(assay_mat)) {
        warning(sprintf("min_samples (%d) is greater than number of samples (%d). This will likely result in 0 rows being kept.",
            min_samples, ncol(assay_mat)), call. = FALSE)
        eff_min_samples <- as.integer(min_samples)  # Don't cap, let it filter to 0
    } else {
        eff_min_samples <- min(as.integer(min_samples), ncol(assay_mat))
    }

    # TPM-based filtering (data is already TPM-normalized from SALMON)
    tokeep <- rowSums(assay_mat >= min_tpm) >= eff_min_samples

    if (verbose) {
        message(sprintf("TPM-based filtering: min_tpm = %.3f in >= %d samples", min_tpm,
            eff_min_samples))
    }

    tokeep
}

# ============================================================================
# HELPER: Filter genes by minimum transcripts per gene Returns updated logical
# vector @noRd
.filter_by_tx_per_gene <- function(tokeep, genes_vec, min_tx_per_gene, verbose = FALSE) {
    if (min_tx_per_gene <= 1 || is.null(genes_vec)) {
        return(tokeep)
    }

    # Count transcripts per gene among those that pass count filtering
    tx_per_gene <- table(genes_vec[tokeep])
    genes_to_keep <- names(tx_per_gene)[tx_per_gene >= min_tx_per_gene]

    # Only keep transcripts from genes with sufficient transcripts
    tokeep_new <- tokeep & (genes_vec %in% genes_to_keep)

    if (verbose && sum(tokeep & !tokeep_new) > 0) {
        message(sprintf("Filtered by min_tx_per_gene = %d: removed %d transcripts",
            min_tx_per_gene, sum(tokeep & !tokeep_new)))
    }

    tokeep_new
}

# ============================================================================
# HELPER: Filter isoforms by relative abundance within genes Implements Soneson
# et al. (2016) isoform-level filtering Returns updated logical vector @noRd
.filter_by_isoform_abundance <- function(tokeep, genes_vec, assay_mat, min_isoform_abundance,
    verbose = FALSE) {
    if (is.null(min_isoform_abundance) || min_isoform_abundance <= 0 || is.null(genes_vec) ||
        length(genes_vec) == 0) {
        return(tokeep)
    }

    # Only consider isoforms that passed previous filtering
    genes_filt <- genes_vec[tokeep]
    assay_mat_filt <- assay_mat[tokeep, , drop = FALSE]

    # For each gene, calculate relative abundance of each isoform
    keep_iso <- rep(FALSE, sum(tokeep))  # length = # rows still kept
    gene_unique <- unique(genes_filt)

    for (gene in gene_unique) {
        # Get indices of all isoforms for this gene (within filtered set)
        gene_idx <- which(genes_filt == gene)

        if (length(gene_idx) == 1) {
            # Single isoform gene - always keep
            keep_iso[gene_idx] <- TRUE
        } else {
            # Calculate mean total count per isoform
            isoform_totals <- rowSums(assay_mat_filt[gene_idx, , drop = FALSE], na.rm = TRUE)
            gene_total <- sum(isoform_totals, na.rm = TRUE)

            if (gene_total > 0) {
                # Relative abundance of each isoform
                rel_abundance <- isoform_totals/gene_total
                # Keep isoforms above threshold
                keep_iso[gene_idx] <- rel_abundance >= min_isoform_abundance
            }
        }
    }

    # Update tokeep: convert boolean from filtered set back to original set
    orig_keep_idx <- which(tokeep)
    tokeep[orig_keep_idx[!keep_iso]] <- FALSE

    if (verbose && sum(!keep_iso) > 0) {
        message(sprintf("Isoform-level filtering (Soneson et al. 2016): min relative abundance = %.1f%%, removed %d isoforms",
            min_isoform_abundance * 100, sum(!keep_iso)))
    }

    tokeep
}

# ============================================================================
# HELPER: Filter genes by valid sample fraction Removes genes sparse in too
# many samples Returns updated logical vector @noRd


# ============================================================================
# HELPER: Synchronize metadata after filtering Updates readcounts, tx2gene,
# tpm, etc. in metadata Returns updated metadata list @noRd
.sync_filter_metadata <- function(md, tokeep, assay_mat, rownames_se = NULL) {
    # Get the rows that are being kept If rownames_se not provided, try to
    # extract from metadata
    if (is.null(rownames_se)) {
        if (!is.null(md$readcounts) && !is.null(rownames(md$readcounts))) {
            rownames_se <- rownames(md$readcounts)
        } else if (!is.null(md$tpm) && !is.null(rownames(md$tpm))) {
            rownames_se <- rownames(md$tpm)
        } else if (!is.null(md$tx2gene) && nrow(md$tx2gene) > 0) {
            rownames_se <- md$tx2gene[[1]]
        } else if (!is.null(md$effective_length) && !is.null(names(md$effective_length))) {
            rownames_se <- names(md$effective_length)
        } else {
            # Fallback: use seq_along(tokeep) if we can't determine rownames
            rownames_se <- seq_along(tokeep)
        }
    }
    kept_rownames <- rownames_se[tokeep]

    # Filter readcounts if present
    if (!is.null(md$readcounts) && is.matrix(md$readcounts)) {
        md$readcounts <- as.matrix(md$readcounts)[tokeep, , drop = FALSE]
    }

    # Filter tx2gene mapping if present
    if (!is.null(md$tx2gene) && is.data.frame(md$tx2gene) && nrow(md$tx2gene) > 0) {
        txmap <- md$tx2gene
        txcol <- colnames(txmap)[1]
        # Keep rows where transcript ID is in the kept rownames
        md$tx2gene <- txmap[txmap[[txcol]] %in% kept_rownames, , drop = FALSE]
    }

    # Filter SALMON metadata (TPM and effective_length) to match filtered assay
    if (!is.null(md$tpm) && is.matrix(md$tpm)) {
        md$tpm <- as.matrix(md$tpm)[tokeep, , drop = FALSE]
    }
    if (!is.null(md$effective_length)) {
        if (is.vector(md$effective_length)) {
            # If it's a named vector, subset by matching names
            if (!is.null(names(md$effective_length))) {
                md$effective_length <- md$effective_length[kept_rownames]
            } else if (is.numeric(md$effective_length) && length(md$effective_length) ==
                length(rownames_se)) {
                # If unnamed vector same length as original rows, subset by
                # index
                md$effective_length <- md$effective_length[tokeep]
            }
        }
    }

    md
}

# ============================================================================
# HELPER: Construct filtered SummarizedExperiment Returns new SE object with
# filtered data @noRd
.construct_filtered_se <- function(se, tokeep, new_assays, new_md, assay_original) {
    # subset rowData if present
    rd <- NULL
    if (nrow(SummarizedExperiment::rowData(se)) > 0) {
        rd <- SummarizedExperiment::rowData(se)[tokeep, , drop = FALSE]
    }

    # Construct new SE with consolidated metadata
    new_se <- SummarizedExperiment::SummarizedExperiment(assays = new_assays, rowData = if (!is.null(rd))
        rd else S4Vectors::DataFrame(), colData = SummarizedExperiment::colData(se), metadata = new_md)

    new_se
}

# ============================================================================
# HELPER: Consolidate parameter resolution (manual vs stringency-based)
# Returns: list(min_samples, min_tpm, min_tx_per_gene, min_isoform_abundance,
# assay_mat, assay_source, genes_vec, pair_col_used) OPTIMIZATION: Eliminates
# 100+ lines of redundant parameter discovery @noRd
.resolve_filter_parameters <- function(se, min_samples, min_tpm, stringency, pair_col,
    tpm_assay_name, assay_name, min_tx_per_gene, min_isoform_abundance, verbose = FALSE) {
    # Validate isoform_abundance parameter
    if (!is.null(min_isoform_abundance)) {
        if (!is.numeric(min_isoform_abundance) || length(min_isoform_abundance) !=
            1 || min_isoform_abundance < 0 || min_isoform_abundance > 1) {
            stop(".filter_se: 'min_isoform_abundance' must be numeric in [0, 1], or NULL to use defaults",
                call. = FALSE)
        }
    }

    # Get assays and discover TPM source
    assays_list <- SummarizedExperiment::assays(se)
    tpm_result <- .get_assay_filtering(se, assays_list, tpm_assay_name, assay_name)
    assay_mat <- tpm_result$mat
    assay_source <- tpm_result$source

    # Warn if fallback (not TPM)
    if (!is.null(tpm_result$is_fallback) && tpm_result$is_fallback) {
        warning("No TPM data found in assays or metadata. Falling back to assay '",
            assay_name, "'.", "\nThis may produce INCORRECT results if '", assay_name,
            "' contains raw counts.", "\nEnsure TPM data is added as an assay or in metadata with tpm/tpm.",
            call. = FALSE)
    }

    # Get gene IDs for downstream filtering
    genes_vec <- .get_gene_ids(se)

    # Handle stringency-based auto-calculation
    pair_col_used <- NULL
    if (!is.null(stringency)) {
        if (!(stringency %in% c("soft", "medium", "severe"))) {
            stop("'stringency' must be one of: 'soft', 'medium', 'severe', or NULL",
                call. = FALSE)
        }

        # Auto-detect or validate pair column
        if (is.null(pair_col)) {
            col_data <- SummarizedExperiment::colData(se)
            pair_candidates <- c("pair", "pair_id", "paired_samples", "subject",
                "subject_id", "individual")
            pair_col <- pair_candidates[pair_candidates %in% colnames(col_data)][1]

            if (is.na(pair_col)) {
                meta <- S4Vectors::metadata(se)
                if (!is.null(meta$coldata) && is.data.frame(meta$coldata)) {
                  pair_col <- pair_candidates[pair_candidates %in% colnames(meta$coldata)][1]
                  if (!is.na(pair_col)) {
                    # Update col_data to metadata$coldata for consistency
                    col_data <- meta$coldata
                  }
                }
            }

            if (is.na(pair_col)) {
                cols_str <- paste(colnames(col_data), collapse = ", ")
                stop("Could not auto-detect pair column in colData or metadata. Available columns: ",
                  cols_str, ". Please specify 'pair_col' parameter.", call. = FALSE)
            }
            if (verbose) {
                message(sprintf("Auto-detected pair column: '%s'", pair_col))
            }
        } else {
            col_data <- SummarizedExperiment::colData(se)
        }

        # Verify pair column exists
        if (!(pair_col %in% colnames(col_data))) {
            cols_str <- paste(colnames(col_data), collapse = ", ")
            stop(sprintf("Pair column '%s' not found. Available columns: %s", pair_col,
                cols_str), call. = FALSE)
        }
        pair_col_used <- pair_col

        # Calculate stringency-based thresholds
        n_samples <- ncol(se)
        n_pairs <- length(unique(col_data[[pair_col]]))
        stringency_result <- .calculate_stringency_thresholds(stringency, n_samples,
            n_pairs)
        min_samples <- stringency_result$min_samples
        min_tx_per_gene <- stringency_result$min_tx_per_gene

        # Use stringency min_isoform_abundance if not explicitly set
        if (is.null(min_isoform_abundance)) {
            min_isoform_abundance <- stringency_result$min_isoform_abundance
        }

        # Estimate min_tpm from data distribution
        tpm_est <- .estimate_min_tpm(assay_mat, stringency, verbose = verbose)
        min_tpm <- tpm_est$min_tpm

        if (verbose) {
            message(sprintf("[resolve_filter_parameters] stringency='%s' -> min_samples=%d, min_tpm=%.3f, min_tx_per_gene=%d",
                stringency, min_samples, min_tpm, min_tx_per_gene))
        }
    }

    list(min_samples = min_samples, min_tpm = min_tpm, min_tx_per_gene = min_tx_per_gene,
        min_isoform_abundance = min_isoform_abundance, assay_mat = assay_mat, assay_source = assay_source,
        genes_vec = genes_vec, pair_col_used = pair_col_used)
}

# ============================================================================
# HELPER: Apply all three filters in sequence with consolidated reporting
# Returns: list(tokeep, before, after) where tokeep is logical vector
# OPTIMIZATION: Consolidates 75 lines of filter chaining @noRd
.apply_combined_filters <- function(assay_mat, genes_vec, params_list, verbose = FALSE) {
    # Extract parameters
    min_tpm <- params_list$min_tpm
    min_samples <- params_list$min_samples
    min_tx_per_gene <- params_list$min_tx_per_gene
    min_isoform_abundance <- params_list$min_isoform_abundance

    before <- nrow(assay_mat)

    # Apply TPM filter
    tokeep <- .apply_tpm_filter(assay_mat, min_tpm, min_samples, verbose = FALSE)
    after_tpm <- sum(tokeep)

    # Apply gene filter (skip if no gene mapping)
    if (!is.null(genes_vec) && length(genes_vec) > 0) {
        tokeep <- .filter_by_tx_per_gene(tokeep, genes_vec, min_tx_per_gene, verbose = FALSE)
        after_tx <- sum(tokeep)

        tokeep <- .filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, min_isoform_abundance,
            verbose = FALSE)
        after_iso <- sum(tokeep)
    } else {
        after_tx <- after_tpm
        after_iso <- after_tpm
    }

    after <- sum(tokeep)

    # Consolidate verbose output
    if (verbose) {
        message(sprintf("TPM-based filtering: min_tpm = %.3f in >= %d samples -> kept %d transcripts",
            min_tpm, min_samples, after_tpm))
        if (!is.null(genes_vec) && length(genes_vec) > 0 && min_tx_per_gene > 1) {
            removed_tx <- after_tpm - after_tx
            if (removed_tx > 0) {
                message(sprintf("Filtered by min_tx_per_gene = %d: removed %d transcripts -> %d remaining",
                  min_tx_per_gene, removed_tx, after_tx))
            }
        }
        if (!is.null(genes_vec) && length(genes_vec) > 0 && !is.null(min_isoform_abundance) &&
            min_isoform_abundance > 0) {
            removed_iso <- after_tx - after_iso
            if (removed_iso > 0) {
                message(sprintf("Isoform-level filtering (Soneson et al. 2016): min relative abundance = %.1f%%, removed %d isoforms -> %d remaining",
                  min_isoform_abundance * 100, removed_iso, after_iso))
            }
        }
    }

    list(tokeep = tokeep, before = before, after = after)
}

# ============================================================================
# HELPER: Finalize filtered SE with synchronized assays and metadata Returns:
# filtered SummarizedExperiment with all metadata synchronized OPTIMIZATION:
# Consolidates post-filtering (20 lines) @noRd
.finalize_filtered_se <- function(se, tokeep, assays_list, genes_vec, min_samples,
    min_tpm, min_tx_per_gene, min_isoform_abundance, stringency) {
    # Subset all assays
    new_assays <- S4Vectors::SimpleList(lapply(assays_list, function(a) {
        if (is.matrix(a) || is.data.frame(a)) {
            as.matrix(a)[tokeep, , drop = FALSE]
        } else {
            a
        }
    }))
    names(new_assays) <- names(assays_list)

    # Synchronize metadata
    md <- S4Vectors::metadata(se)
    assay_mat <- assays_list[[1]]  # Use first assay for rownames
    md <- .sync_filter_metadata(md, tokeep, assay_mat, rownames(se))

    # Add filtering record
    md$filtered <- list(min_samples = min_samples, min_tpm = min_tpm, min_tx_per_gene = min_tx_per_gene,
        min_isoform_abundance = min_isoform_abundance, stringency = stringency)

    # Construct and return
    .construct_filtered_se(se, tokeep, new_assays, md, assay_mat)
}



#' Subset a TSENATAnalysis Object for Testing and Examples
#'
#' Create a smaller, representative subset of a TSENATAnalysis object
#' for use in testing, examples, or documentation. Preserves all analysis
#' metadata and computed results while reducing dataset size.
#'
#' @param analysis TSENATAnalysis object to subset
#' @param n_genes Positive integer. Number of genes to retain (default: 10).
#'   Set to NULL to keep all genes.
#' @param n_samples Positive integer. Number of samples to retain (default:
#' NULL,
#'   keep all samples). If specified, samples are selected to balance
#'   conditions when possible.
#' @param genes Character vector of specific gene IDs to retain. If provided,
#'   overrides n_genes argument (default: NULL).
#' @param samples Character vector of specific sample IDs to retain. If
#' provided,
#'   overrides n_samples argument (default: NULL).
#' @param select_by One of 'variance' (select genes with highest variance),
#'   'mean' (select genes with highest mean expression), or 'random'
#'   (random selection). Default: 'variance'
#' @param seed Random seed for reproducible subsetting (default: 42)
#' @param min_count Minimum total transcript count (across all samples) required
#'   for a gene to be included. Filters genes to ensure adequate data density
#'   for statistical operations like jackknife (default: NULL, no filtering).
#'   Useful values: 5-10 for robust estimates.
#' @param verbose Logical. Print progress messages (default: FALSE)
#'
#' @return A new TSENATAnalysis object containing only the specified genes
#'   and samples. All computed results (diversity, LM, jackknife, divergence)
#'   are automatically subsetted to match the new dimensions while preserving
#'   analysis configuration and metadata.
#'
#' @details
#' This function provides a convenient wrapper around the `[` subsetting
#' operator for TSENATAnalysis. It simplifies common use cases:
#'
#' \itemize{
#'   \item \strong{By gene count}: Automatically selects top genes by variance,
#'     mean expression, or at random
#'   \item \strong{By sample count}: Intelligently balances sample selection
#'     across experimental conditions
#'   \item \strong{By specific IDs}:  Explicitly specify genes and 
#' samples to keep
#'   \item \strong{By statistic}: Select informative genes (high variance = more
#'     informative for testing)
#'   \item \strong{By abundance}: Filter by minimum total count to ensure
#'     adequate data density for operations like jackknife
#' }
#'
#' The function preserves all analysis structure:
#' - Diversity results are subsetted to match gene and sample selection
#' - Jackknife results maintain confidence intervals for selected samples
#' - LM results (gene-level statistics) remain intact
#' - Divergence calculations are updated if applicable
#' - Configuration and metadata are unchanged
#'
#' @examples
#' # Create minimal example TSENATAnalysis
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(1000, 10), nrow = 100, ncol = 10)),
#'   rowData = data.frame(
#'     gene_id = rep(paste0('GENE_', 1:10), 10),
#'     row.names = paste0('TX_', 1:100)
#'   ),
#'   colData = data.frame(
#'     sample_id = paste0('S_', 1:10),
#'     condition = rep(c('control', 'treatment'), 5),
#'     row.names = paste0('S_', 1:10)
#'   )
#' )
#' analysis <- TSENATAnalysis(se)
#'
#' # Subset to top genes (use 200 to ensure adequate data for downstream
#' analysis)
#' small_analysis <- filter_analysis_s4(analysis, min_samples = 1,
#' subset_n_genes = 200)
#'
#' # Subset to genes with minimum 1 total count across all samples
#' # This ensures data adequacy filtering (recommended: min_count = 10-20
#' for robust estimates)
#' filtered <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 100, subset_min_count = 1)
#'
#' # Subset to specific genes only
#' subset_genes <- filter_analysis_s4(
#'   analysis,
#'   min_samples = 1,
#'   subset_genes = c('TX_1', 'TX_2', 'TX_3'),
#'   subset_n_samples = 5
#' )
#'
#' # Random selection of genes (reproducible with seed)
#' random_subset <- filter_analysis_s4(
#'   analysis,
#'   min_samples = 1,
#'   subset_n_genes = 8,
#'   subset_n_samples = 6,
#'   subset_select_by = 'random',
#'   subset_seed = 123
#' )
#'
#' # Keep specific samples only
#' control_only <- filter_analysis_s4(
#'   analysis,
#'   min_samples = 1,
#'   subset_samples = colnames(se(analysis))[
#'     colData(se(analysis))$condition == 'control'
#'   ]
#' )
#'
#' @noRd
#' @keywords internal

# ============================================================================
# Helper: Select genes from analysis by criterion or explicit list
# ============================================================================
.select_genes_from_analysis <- function(se, n_genes, genes, select_by, seed, verbose) {
    n_genes_total <- nrow(se)

    if (!is.null(genes)) {
        # Use explicitly provided genes
        if (verbose) {
            genes_str <- paste(head(genes, 3), collapse = ", ")
            message("[subset_analysis] Using provided genes: ", genes_str, if (length(genes) >
                3)
                "...")
        }
        gene_idx <- match(genes, rownames(se))
        if (any(is.na(gene_idx))) {
            missing_genes <- genes[is.na(gene_idx)]
            missing_str <- paste(head(missing_genes, 3), collapse = ", ")
            stop("Genes not found in analysis: ", missing_str, call. = FALSE)
        }
    } else if (!is.null(n_genes)) {
        # Select by criterion
        n_genes <- as.integer(n_genes)
        if (n_genes < 1) {
            stop("n_genes must be >= 1", call. = FALSE)
        }
        if (n_genes > n_genes_total) {
            warning("n_genes (", n_genes, ") exceeds available genes (", n_genes_total,
                "). Using all genes.", call. = FALSE)
            n_genes <- n_genes_total
        }

        counts_matrix <- assay(se, "counts")

        if (select_by == "variance") {
            # Compute variance per transcript
            tx_vars <- matrixStats::rowVars(counts_matrix)
            gene_idx <- order(tx_vars, decreasing = TRUE)[seq_len(n_genes)]
            if (verbose) {
                message("[subset_analysis] Selected ", n_genes, " genes with highest variance")
            }
        } else if (select_by == "mean") {
            # Compute mean per transcript
            tx_means <- rowMeans(counts_matrix)
            gene_idx <- order(tx_means, decreasing = TRUE)[seq_len(n_genes)]
            if (verbose) {
                message("[subset_analysis] Selected ", n_genes, " genes with highest mean expression")
            }
        } else if (select_by == "random") {
            # Random selection
            withr::local_seed(seed)
            gene_idx <- sample(n_genes_total, n_genes)
            if (verbose) {
                message("[subset_analysis] Randomly selected ", n_genes, " genes (seed = ",
                  seed, ")")
            }
        }
    } else {
        # Keep all genes
        gene_idx <- seq_len(n_genes_total)
        if (verbose) {
            message("[subset_analysis] Keeping all ", n_genes_total, " genes")
        }
    }

    return(gene_idx)
}

# ============================================================================
# Helper: Detect paired structure in coldata
# ============================================================================
.detect_paired_structure <- function(coldata) {
    # Look for paired sample columns: paired_samples, pair, pairing, etc.
    pair_cols <- c("paired_samples", "pair", "pairing", "pair_id", "subject")
    for (col in pair_cols) {
        if (col %in% colnames(coldata)) {
            # Verify it's a valid pairing: check if ANY pairs appear in each
            # condition
            if ("condition" %in% colnames(coldata)) {
                conditions <- unique(coldata$condition)
                pair_condition_table <- table(coldata[[col]], coldata$condition)

                # A paired design exists if there are pairs that appear in
                # multiple conditions Find 'complete pairs' - those that appear
                # in all conditions
                complete_pairs <- rownames(pair_condition_table)[rowSums(pair_condition_table >
                  0) == length(conditions)]

                if (length(complete_pairs) > 0) {
                  return(list(is_paired = TRUE, pair_col = col, complete_pairs = complete_pairs))
                }
            }
        }
    }
    return(list(is_paired = FALSE, pair_col = NULL, complete_pairs = NULL))
}

# ============================================================================
# Helper: Select complete pairs with balanced representation
# ============================================================================
.select_pairs_balanced <- function(coldata, n_samples, pair_col, seed, verbose = FALSE) {
    # For paired samples, select complete pairs balanced across conditions

    # Find truly complete pairs (appear in all conditions)
    conditions <- unique(coldata$condition)
    pair_condition_table <- table(coldata[[pair_col]], coldata$condition)
    complete_pairs <- rownames(pair_condition_table)[rowSums(pair_condition_table >
        0) == length(conditions)]

    if (length(complete_pairs) == 0) {
        if (verbose) {
            message("[subset_analysis] No complete pairs found. Falling back to condition-balanced selection.")
        }
        # Fall back to non-paired selection
        return(NULL)
    }

    if (n_samples%%2 != 0) {
        n_samples <- n_samples + 1  # Round up to even number for complete pairs
        if (verbose) {
            message("[subset_analysis] Rounding up to ", n_samples, " samples for complete pair selection")
        }
    }

    n_pairs_needed <- n_samples/2
    n_pairs_available <- length(complete_pairs)

    if (n_pairs_needed > n_pairs_available) {
        if (verbose) {
            message("[subset_analysis] Not enough complete pairs (", n_pairs_available,
                ") for requested samples (", n_samples, "). Falling back to condition-balanced selection.")
        }
        # Fall back to non-paired selection
        return(NULL)
    }

    # Select complete pairs randomly
    withr::local_seed(seed)
    selected_pairs <- sample(complete_pairs, n_pairs_needed)

    # Get all samples belonging to selected COMPLETE pairs
    sample_idx <- which(coldata[[pair_col]] %in% selected_pairs)

    if (verbose) {
        message("[subset_analysis] Selected ", n_pairs_needed, " complete pairs ",
            "(", length(sample_idx), " total samples with preserved pairing)")
    }

    return(sample_idx)
}

# ============================================================================
# Helper: Balance sample selection across conditions
# ============================================================================
.balance_sample_selection <- function(coldata, n_samples, seed, verbose = FALSE) {
    n_samples_total <- nrow(coldata)

    if (n_samples > n_samples_total) {
        warning("n_samples (", n_samples, ") exceeds available samples (", n_samples_total,
            "). Using all samples.", call. = FALSE)
        return(seq_len(n_samples_total))
    }

    # BUG FIX (April 2026): Preserve paired structure when subsampling Detect
    # if this is a paired design and handle accordingly
    pairing_info <- .detect_paired_structure(coldata)
    if (pairing_info$is_paired && length(pairing_info$complete_pairs) > 0) {
        if (verbose) {
            message("[subset_analysis] Detected ", length(pairing_info$complete_pairs),
                " complete pairs in '", pairing_info$pair_col, "' column. Preserving pairs during selection.")
        }

        pair_samples <- .select_pairs_balanced(coldata = coldata, n_samples = n_samples,
            pair_col = pairing_info$pair_col, seed = seed, verbose = verbose)

        # If pair selection succeeds, use it
        if (!is.null(pair_samples)) {
            return(pair_samples)
        }
        # Otherwise fall through to non-paired selection
        if (verbose) {
            message("[subset_analysis] Pair selection failed, falling back to condition-balanced selection")
        }
    }

    # Try to balance by condition if available (original logic)
    if ("condition" %in% colnames(coldata)) {
        conditions <- unique(coldata$condition)
        samples_per_cond <- ceiling(n_samples/length(conditions))

        withr::local_seed(seed)
        sample_idx <- c()
        for (cond in conditions) {
            cond_idx <- which(coldata$condition == cond)
            selected <- sample(cond_idx, min(samples_per_cond, length(cond_idx)))
            sample_idx <- c(sample_idx, selected)
        }
        return(head(sample_idx, n_samples))
    } else if ("sample_type" %in% colnames(coldata)) {
        # Alternative: try sample_type
        types <- unique(coldata$sample_type)
        samples_per_type <- ceiling(n_samples/length(types))

        withr::local_seed(seed)
        sample_idx <- c()
        for (type in types) {
            type_idx <- which(coldata$sample_type == type)
            selected <- sample(type_idx, min(samples_per_type, length(type_idx)))
            sample_idx <- c(sample_idx, selected)
        }
        return(head(sample_idx, n_samples))
    } else {
        # Random selection
        withr::local_seed(seed)
        return(sample(n_samples_total, n_samples))
    }
}

# ============================================================================
# Helper: Select samples from analysis
# ============================================================================
.select_samples_from_analysis <- function(se, n_samples, samples, seed, verbose) {
    n_samples_total <- ncol(se)

    if (!is.null(samples)) {
        # Use explicitly provided samples
        if (verbose) {
            samples_str <- paste(head(samples, 3), collapse = ", ")
            message("[subset_analysis] Using provided samples: ", samples_str, if (length(samples) >
                3)
                "...")
        }
        sample_idx <- match(samples, colnames(se))
        if (any(is.na(sample_idx))) {
            missing_samples <- samples[is.na(sample_idx)]
            missing_str <- paste(head(missing_samples, 3), collapse = ", ")
            stop("Samples not found in analysis: ", missing_str, call. = FALSE)
        }
    } else if (!is.null(n_samples)) {
        # Intelligently select samples across conditions
        n_samples <- as.integer(n_samples)
        if (n_samples < 1) {
            stop("n_samples must be >= 1", call. = FALSE)
        }

        coldata <- colData(se)
        sample_idx <- .balance_sample_selection(coldata, n_samples, seed, verbose)

        if (verbose) {
            n_conditions <- if ("condition" %in% colnames(coldata))
                length(unique(coldata$condition)) else if ("sample_type" %in% colnames(coldata))
                length(unique(coldata$sample_type)) else 1

            if (n_conditions > 1) {
                message("[subset_analysis] Selected ", n_samples, " samples balanced across ",
                  n_conditions, " conditions")
            } else {
                message("[subset_analysis] Randomly selected ", n_samples, " samples (seed = ",
                  seed, ")")
            }
        }
    } else {
        # Keep all samples
        sample_idx <- seq_len(n_samples_total)
        if (verbose) {
            message("[subset_analysis] Keeping all ", n_samples_total, " samples")
        }
    }

    return(sample_idx)
}

# ============================================================================
# Helper: Synchronize subset metadata with analysis object
# ============================================================================
.sync_subset_metadata <- function(analysis_obj, gene_idx, sample_idx, verbose) {
    se <- analysis_obj@se
    analysis_se <- analysis_obj@se

    # 1. Sync tx2gene mapping (transcript-to-gene mapping)
    if (!is.null(S4Vectors::metadata(se)$tx2gene)) {
        tx2gene_full <- S4Vectors::metadata(se)$tx2gene
        tx2gene_subset <- tx2gene_full[tx2gene_full$Transcript %in% rownames(analysis_se),
            ]
        S4Vectors::metadata(analysis_se)$tx2gene <- tx2gene_subset

        if (verbose) {
            message("[subset_analysis] Filtered tx2gene: ", nrow(tx2gene_full), " -> ",
                nrow(tx2gene_subset), " transcripts")
        }
    }

    # 2. Sync readcounts (original count matrix stored in metadata)
    if (!is.null(S4Vectors::metadata(se)$readcounts)) {
        readcounts_full <- S4Vectors::metadata(se)$readcounts
        readcounts_subset <- readcounts_full[gene_idx, sample_idx, drop = FALSE]
        S4Vectors::metadata(analysis_se)$readcounts <- readcounts_subset

        if (verbose) {
            message("[subset_analysis] Filtered readcounts: ", nrow(readcounts_full),
                " -> ", nrow(readcounts_subset), " transcripts")
        }
    }

    # 3. Sync tpm (TPM matrix stored in metadata)
    if (!is.null(S4Vectors::metadata(se)$tpm)) {
        tpm_full <- S4Vectors::metadata(se)$tpm
        tpm_subset <- tpm_full[gene_idx, sample_idx, drop = FALSE]
        S4Vectors::metadata(analysis_se)$tpm <- tpm_subset

        if (verbose) {
            message("[subset_analysis] Filtered tpm: ", nrow(tpm_full), " -> ", nrow(tpm_subset),
                " transcripts")
        }
    }

    # 4. Sync effective_length (effective lengths)
    if (!is.null(S4Vectors::metadata(se)$effective_length)) {
        eff_len_full <- S4Vectors::metadata(se)$effective_length

        # Could be vector or matrix, handle both
        if (is.vector(eff_len_full)) {
            # Named vector - subset by matching names to selected genes
            tx_names_subset <- rownames(analysis_se)
            eff_len_subset <- eff_len_full[na.omit(match(tx_names_subset, names(eff_len_full)))]
            if (length(eff_len_subset) > 0) {
                S4Vectors::metadata(analysis_se)$effective_length <- eff_len_subset
            }
        } else if (is.matrix(eff_len_full)) {
            # Matrix - subset by rows and columns
            eff_len_subset <- eff_len_full[gene_idx, sample_idx, drop = FALSE]
            S4Vectors::metadata(analysis_se)$effective_length <- eff_len_subset
        }
    }

    return(analysis_obj)
}

.subset_analysis <- function(analysis, n_genes = 10, n_samples = NULL, genes = NULL,
    samples = NULL, select_by = c("variance", "mean", "random"), seed = 42, min_count = NULL,
    verbose = FALSE) {

    # Validate input
    if (!inherits(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    select_by <- match.arg(select_by)
    withr::local_seed(seed)

    se <- analysis@se

    # Use helper functions to select genes and samples
    gene_idx <- .select_genes_from_analysis(se = se, n_genes = n_genes, genes = genes,
        select_by = select_by, seed = seed, verbose = verbose)

    sample_idx <- .select_samples_from_analysis(se = se, n_samples = n_samples, samples = samples,
        seed = seed, verbose = verbose)

    # ============================================================================
    # Filter by minimum count (data adequacy check)
    # ============================================================================

    if (!is.null(min_count)) {
        min_count <- as.numeric(min_count)
        if (min_count < 0) {
            stop("min_count must be >= 0", call. = FALSE)
        }

        counts_matrix <- assay(se, "counts")[gene_idx, sample_idx, drop = FALSE]
        total_counts <- rowSums(counts_matrix)
        genes_keep <- total_counts >= min_count

        n_before_filter <- length(gene_idx)
        gene_idx_filtered <- gene_idx[genes_keep]
        n_after_filter <- length(gene_idx_filtered)
        n_removed <- n_before_filter - n_after_filter

        if (n_removed > 0) {
            warning("Filtered out ", n_removed, " gene(s) with total count < ", min_count,
                " (", n_after_filter, " genes remain)", call. = FALSE)
            if (verbose) {
                message("[subset_analysis] Total count range in selected samples: ",
                  format(min(total_counts), trim = TRUE), " - ", format(max(total_counts),
                    trim = TRUE))
                message("[subset_analysis] Retained genes with count >= ", min_count,
                  ": ", n_after_filter)
            }
        } else if (verbose) {
            message("[subset_analysis] All ", n_after_filter, " genes meet minimum count threshold (",
                min_count, ")")
        }

        gene_idx <- gene_idx_filtered

        if (length(gene_idx) == 0) {
            stop("No genes meet the minimum count threshold (min_count = ", min_count,
                "). Consider lowering min_count or using more/different samples.",
                call. = FALSE)
        }
    }

    # ============================================================================
    # Subset using the [ operator
    # ============================================================================

    analysis_subset <- analysis[gene_idx, sample_idx]

    # Sync metadata to match the subsetted object using helper function
    analysis_subset <- .sync_subset_metadata(analysis_obj = analysis_subset, gene_idx = gene_idx,
        sample_idx = sample_idx, verbose = verbose)

    if (verbose) {
        message("[subset_analysis] Subset complete: ", nrow(analysis_subset@se),
            " genes x ", ncol(analysis_subset@se), " samples")
    }

    return(analysis_subset)
}
