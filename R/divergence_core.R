# =========================================================================
# INTERNAL FUNCTION: Calculate Tsallis Divergence
# =========================================================================
# 
# This is an internal helper called by the public calculate_divergence() wrapper.
# NOT exported. Do not document with roxygen2.
#
# Architecture: Transcript-level counts -> Gene-level aggregation -> Tsallis divergence
# 
# Computes Tsallis divergence comparing two groups across multiple genes,
# with optional bootstrap confidence intervals. Automatically aggregates
# transcript-level counts to gene-level (per Paper I033).
# 
# Returns a SummarizedExperiment object containing:
# - assay: genes * q matrix of divergence estimates (one per q value)
# - rowData: gene metadata including per-q divergence estimates, CIs
# - colData: one row per q value with q-specific metadata
# - metadata: processing parameters and summary statistics
#
# Parameters:
#   se: SummarizedExperiment object with transcript-level counts
#   group_col: Character; colData column for group membership (auto-detect if NULL)
#   control_group: Character; reference group name (auto-detect if NULL)
#   q: Tsallis parameter (scalar or vector, default: 1)
#   paired: Logical; paired sample design (default: FALSE)
#   bootstrap: Logical; compute bootstrap confidence intervals (default: FALSE)
#   nboot: Number of bootstrap replicates (default: 1000)
#   ci: Confidence level (default: 0.95)
#   method: Bootstrap method ('percentile' or 'bca', default: 'percentile')
#   log_base: Logarithm base (default: exp(1))
#   norm: Normalization/standardization mode (default: TRUE for range)
#   pseudocount: Pseudocount for stability (default: 0.5)
#   nthreads: Number of CPU threads (default: 1, NULL for auto-detect)
#   progress: Logical; show progress bar (default: TRUE)
.calculate_divergence <- function(se, group_col = NULL, control_group = NULL, q = 1,
    paired = FALSE, bootstrap = FALSE, nboot = "auto", ci = 0.95, method = "percentile",
    norm = TRUE, log_base = exp(1), pseudocount = 0.5, nthreads = 1, progress = FALSE) {

    # =========================================================================
    # INPUT VALIDATION - Must be done BEFORE implementation Following
    # Bioconductor guidelines: fail fast with clear error messages
    # =========================================================================

    # Validate bootstrap parameter
    if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
        stop("bootstrap must be a logical", call. = FALSE)
    }

    # Validate paired parameter
    if (!is.logical(paired) || length(paired) != 1 || is.na(paired)) {
        stop("paired must be a logical", call. = FALSE)
    }

    # Validate method parameter
    if (!is.null(method) && (!is.character(method) || length(method) == 0 || is.na(method[1]))) {
        stop("method must be a character string", call. = FALSE)
    }

    # Validate nthreads parameter
    if (!is.null(nthreads)) {
        if (!is.numeric(nthreads) || length(nthreads) != 1 || is.na(nthreads)) {
            stop("nthreads must be a positive integer", call. = FALSE)
        }
    }

    # Validate progress parameter
    if (!is.logical(progress) || length(progress) != 1 || is.na(progress)) {
        stop("progress must be a logical", call. = FALSE)
    }

    # Call implementation directly - errors will propagate clearly
    .calculate_divergence_impl(se, group_col, control_group, q, paired, bootstrap,
        nboot, ci, method, norm, log_base, pseudocount, nthreads, progress)
}

# ========================================================================= NEW
# ORCHESTRATOR HELPERS (March 2026 refactoring)
# =========================================================================

#' Validate group columns and auto-detect missing ones.
#'
#' Resolution order:
#' 1. If group_col is provided  → use it as-is (user knows best)
#' 2. If control_group is provided → find the colData column that contains it
#' 3. If both are NULL → full auto-detection via .auto_detect_groups()
#'
#' @noRd
.validate_and_auto_detect_groups <- function(se, group_col, control_group, progress) {
    cd <- SummarizedExperiment::colData(se)
    cd_colnames <- colnames(cd)

    # Candidate group columns in priority order
    group_col_candidates <- c("sample_type", "group", "condition", "treatment",
        "phenotype")

    # ---- Resolve group_col ----
    if (is.null(group_col)) {
        if (!is.null(control_group)) {
            # control_group is known → find the column that contains it
            found <- FALSE
            for (col_name in group_col_candidates) {
                if (col_name %in% cd_colnames &&
                    control_group %in% as.character(cd[[col_name]])) {
                    group_col <- col_name
                    found <- TRUE
                    if (progress)
                        message("[calculate_divergence] Using group_col='", group_col,
                            "' (contains control_group='", control_group, "')")
                    break
                }
            }
            if (!found) {
                stop("control_group='", control_group,
                    "' not found in any candidate group column (",
                    paste(group_col_candidates, collapse = ", "), "). ",
                    "Please specify 'group_col' explicitly.",
                    call. = FALSE)
            }
        } else {
            # Neither provided → full auto-detection
            auto_groups <- .auto_detect_groups(se)
            if (is.na(auto_groups$group_col)) {
                stop("Could not auto-detect group column in colData. ",
                    "Available columns: ", paste(cd_colnames, collapse = ", "),
                    ". Please specify 'group_col' explicitly.", call. = FALSE)
            }
            group_col <- auto_groups$group_col
            if (progress)
                message("[calculate_divergence] Auto-detected group_col='", group_col, "'")
        }
    }

    # ---- Resolve control_group ----
    if (is.null(control_group)) {
        # Not provided → only auto-detect if group_col also wasn't auto-detected above
        # (if group_col was found via auto-detection above, control_group was too)
        if (exists("auto_groups") && !is.null(auto_groups)) {
            # Reuse result from above
            if (is.na(auto_groups$control_group)) {
                groups_in_col <- unique(as.character(cd[[group_col]]))
                stop("Could not auto-detect control_group. ",
                    "Available groups: ", paste(sQuote(groups_in_col), collapse = ", "),
                    ". Please specify 'control_group' explicitly ",
                    "(e.g., control_group = \"", groups_in_col[1], "\").",
                    call. = FALSE)
            }
            control_group <- auto_groups$control_group
        } else {
            # group_col was provided but control_group wasn't → run detection
            auto_groups <- .auto_detect_groups(se)
            if (is.na(auto_groups$control_group)) {
                groups_in_col <- unique(as.character(cd[[group_col]]))
                stop("Could not auto-detect control_group. ",
                    "Available groups: ", paste(sQuote(groups_in_col), collapse = ", "),
                    ". Please specify 'control_group' explicitly ",
                    "(e.g., control_group = \"", groups_in_col[1], "\").",
                    call. = FALSE)
            }
            control_group <- auto_groups$control_group
        }
        if (progress)
            message("[calculate_divergence] Auto-detected control_group='", control_group, "'")
    }

    list(group_col = group_col, control_group = control_group)
}

#' Prepare genes for processing (identification + extraction)
#' Consolidates gene column identification and unique gene extraction

#' @noRd
.prepare_genes_processing <- function(se) {
    rd <- SummarizedExperiment::rowData(se)
    gene_col <- .identify_gene_column(se)
    all_gene_names <- .extract_gene_list(se, gene_col)
    gene_indices <- seq_along(all_gene_names)

    if (length(gene_indices) == 0) {
        stop("No genes to process. ", "se gene names (first 3): ", paste(head(all_gene_names,
            3), collapse = ", "))
    }

    list(gene_col = gene_col, all_gene_names = all_gene_names, gene_indices = gene_indices,
        rd = rd, num_genes = length(gene_indices))
}

#' Configure bootstrap and parallel execution parameters
#' Fixes nboot bug and consolidates configuration logic
#' @noRd

# NOTE (March 2026): .bootstrap_configure_parallel() moved to bootstrap.R

#' @noRd
.prepare_divergence_execution <- function(se, bootstrap, paired, nboot, method, nthreads,
    progress) {
    pair_ids <- NULL
    pairing_info <- ""

    if (isTRUE(bootstrap)) {
        # Use isTRUE to safely handle NA Auto-detect paired samples
        pair_detected <- .detect_pair_ids(se)

        if (pair_detected$num_pairs > 0) {
            pair_ids <- pair_detected$pair_ids
            pairing_info <- sprintf(" [paired: %d unique pairs from '%s' column]",
                pair_detected$num_pairs, pair_detected$column_name)

            if (isFALSE(paired) && progress) {
                # Use isFALSE to safely handle NA
                message("NOTE: Paired sample structure detected in '", pair_detected$column_name,
                  "' column.\n", "      Using pair-respecting bootstrap resampling.")
            }
        } else {
            if (isTRUE(paired) && progress) {
                # Use isTRUE to safely handle NA
                message("paired=TRUE but no pair ID column detected in colData.",
                  " Using independent bootstrap resampling instead.")
            }
        }
    }

    if (progress) {
        mode_desc <- if (isTRUE(bootstrap)) {
            # Use isTRUE to safely handle NA
            paste0("bootstrap with ", nboot, " replicates (", method, ")", pairing_info)
        } else {
            "point estimates only"
        }

        mode_str <- if (!is.null(nthreads) && nthreads > 1)
            "Parallel" else "Sequential"
        thread_desc <- if (!is.null(nthreads) && nthreads > 1)
            paste0(" on ", nthreads, " threads") else ""
        message(mode_str, " mode: ", mode_desc, thread_desc)
    }

    list(pair_ids = pair_ids, pairing_info = pairing_info)
}

#' Execute divergence computation (abstracted seq vs parallel dispatch)
#' Consolidates nearly-identical sequential and parallel blocks

#' @noRd
.compute_divergence_worker <- function(gene_indices, all_gene_names, se, gene_col,
    rd, group_col, control_group, q, nboot, ci, method, log_base, pseudocount, pair_ids,
    nthreads, use_parallel, progress) {
    start_time <- Sys.time()
    num_genes <- length(gene_indices)

    # Define per-gene computation function
    compute_gene_divergence <- function(i) {
        gene_idx <- gene_indices[i]

        .process_single_gene_div(gene_idx, all_gene_names, se, gene_col, rd, group_col,
            control_group, q, nboot, ci, method, log_base, pseudocount, pair_ids)
    }

    # Execute using BiocParallel infrastructure
    results_list <- .bplapply(X = seq_along(gene_indices), FUN = compute_gene_divergence,
        nthreads = nthreads)

    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    list(results = results_list, elapsed = elapsed)
}

#' Finalize result matrices from results list
#' Consolidates matrix initialization, population, and reference q handling

#' @noRd
.finalize_divergence_matrices <- function(results_list, num_genes, q, norm, progress) {
    # Initialize result matrices
    matrices <- .initialize_matrices(num_genes, q)
    assay_matrix <- matrices$assay
    row_data_df <- matrices$rowData

    # Populate matrices from results_list
    populated <- .populate_matrices(results_list, assay_matrix, row_data_df, q)
    assay_matrix <- populated$assay
    row_data_df <- populated$rowData

    # Set row names
    rownames(row_data_df) <- row_data_df$gene_name
    rownames(assay_matrix) <- row_data_df$gene_name

    # Populate generic estimate/lower_ci/upper_ci columns using reference q
    # value (q=1)
    q_ref <- 1
    q_idx <- which.min(abs(q - q_ref))
    if (length(q_idx) > 0 && q_idx <= length(q)) {
        ref_q <- q[q_idx]
        estimate_col <- paste0("estimate_q", ref_q)
        lower_ci_col <- paste0("lower_ci_q", ref_q)
        upper_ci_col <- paste0("upper_ci_q", ref_q)
        ci_width_col <- paste0("ci_width_q", ref_q)

        if (estimate_col %in% colnames(row_data_df)) {
            row_data_df$estimate <- row_data_df[[estimate_col]]
            row_data_df$lower_ci <- row_data_df[[lower_ci_col]]
            row_data_df$upper_ci <- row_data_df[[upper_ci_col]]
            row_data_df$ci_width <- row_data_df[[ci_width_col]]
        }
    }

    # Apply normalization if requested
    if (norm != "none") {
        if (progress) {
            message(sprintf("Applying '%s' normalization to divergence estimates...",
                norm))
        }

        normalized <- .normalize_divergence_matrix(assay_matrix = assay_matrix, row_data_df = row_data_df,
            q_vals = q, norm = norm)
        assay_matrix <- normalized$assay
        row_data_df <- normalized$rowData
    }

    # Classify per-q patterns
    row_data_df$per_q_pattern <- NA_character_
    if (length(q) > 1) {
        for (i in seq_len(nrow(row_data_df))) {
            per_q_divs <- assay_matrix[i, ]
            names(per_q_divs) <- paste0("q_", q)

            if (sum(!is.na(per_q_divs)) >= 2) {
                pattern_result <- .classify_q_pattern(per_q_divs)
                row_data_df$per_q_pattern[i] <- if (is.na(pattern_result$pattern))
                  "UNCLASSIFIED" else pattern_result$pattern
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Print divergence computation summary
#' Consolidates logging and summary statistics reporting

#' @noRd
.print_divergence_summary <- function(num_genes, num_errors, elapsed, row_data_df,
    progress) {
    num_success <- num_genes - num_errors

    if (progress) {
        message("\nDIVERGENCE COMPUTATION COMPLETE")
        message("Summary:")
        message("  Genes processed:        ", num_genes)
        message("  Successful:             ", num_success)
        message("  Failed:                 ", num_errors)
        message("  Total elapsed time:     ", sprintf("%.1f seconds", elapsed))
        message("  Average per gene:       ", sprintf("%.2f seconds", elapsed/num_genes))
        message("  Genes per minute:       ", sprintf("%.1f", (num_genes/elapsed) *
            60))

        if (num_errors > 0) {
            message("Failed genes:")
            failed <- row_data_df[!is.na(row_data_df$error), ]
            for (i in seq_len(min(10, nrow(failed)))) {
                message(sprintf("  [%d] %s: %s", i, failed$gene_name[i], failed$error[i]))
            }
            if (num_errors > 10) {
                message("  ... and", num_errors - 10, "more")
            }
        }
    }
}

#' Implementation of calculate_divergence with parameter validation

#' @noRd
.calculate_divergence_impl <- function(se, group_col = NULL, control_group = NULL,
    q = 1, paired = FALSE, bootstrap = FALSE, nboot = "auto", ci = 0.95, method = "percentile",
    norm = TRUE, log_base = exp(1), pseudocount = 0.5, nthreads = 1, progress = FALSE) {

    # =========================================================================
    # =========================================================================
    # INPUT VALIDATION & SETUP
    # =========================================================================

    # Parameters are validated in .calculate_divergence() Normalize/coerce for
    # internal use
    if (!isTRUE(bootstrap)) {
        bootstrap <- FALSE
    }
    if (!isTRUE(paired)) {
        paired <- FALSE
    }

    method <- if (!is.null(method) && is.character(method) && length(method) > 0) {
        as.character(method[1])
    } else {
        "percentile"
    }

    nthreads <- if (is.numeric(nthreads) && length(nthreads) == 1 && !is.na(nthreads)) {
        as.integer(max(1, nthreads))
    } else if (is.null(nthreads)) {
        1L
    } else {
        1L
    }

    progress <- if (is.logical(progress) && length(progress) == 1) {
        progress
    } else {
        FALSE
    }

    norm <- .validate_norm_parameter(norm)
    q <- .validate_and_sort_q_values(q)
    .validate_se_input(se)

    # =========================================================================
    # AUTO-DETECT GROUP COLUMN AND CONTROL GROUP
    # =========================================================================

    group_info <- .validate_and_auto_detect_groups(se, group_col, control_group,
        progress)
    group_col <- group_info$group_col
    control_group <- group_info$control_group

    # =========================================================================
    # PREPARE GENES FOR PROCESSING
    # =========================================================================

    genes_info <- .prepare_genes_processing(se)
    gene_col <- genes_info$gene_col
    all_gene_names <- genes_info$all_gene_names
    gene_indices <- genes_info$gene_indices
    rd <- genes_info$rd
    num_genes <- genes_info$num_genes

    # =========================================================================
    # BOOTSTRAP & PARALLEL CONFIGURATION
    # =========================================================================

    boot_config <- .bootstrap_configure_parallel(bootstrap, nboot, method, num_genes,
        nthreads, progress)
    nboot <- boot_config$nboot
    nthreads <- boot_config$nthreads
    use_parallel <- boot_config$use_parallel

    # =========================================================================
    # PAIRED SAMPLE DETECTION & EXECUTION SETUP
    # =========================================================================

    exec_setup <- .prepare_divergence_execution(se, bootstrap, paired, nboot, method,
        nthreads, progress)
    pair_ids <- exec_setup$pair_ids

    # =========================================================================
    # EXECUTE DIVERGENCE COMPUTATION (SEQUENTIAL OR PARALLEL)
    # =========================================================================

    comp_result <- .compute_divergence_worker(gene_indices, all_gene_names, se, gene_col,
        rd, group_col, control_group, q, nboot, ci, method, log_base, pseudocount,
        pair_ids, nthreads, use_parallel, progress)
    results_list <- comp_result$results
    elapsed <- comp_result$elapsed

    # =========================================================================
    # FINALIZE MATRICES & APPLY NORMALIZATION
    # =========================================================================

    matrices_final <- .finalize_divergence_matrices(results_list, num_genes, q, norm,
        progress)
    assay_matrix <- matrices_final$assay
    row_data_df <- matrices_final$rowData

    # =========================================================================
    # SUMMARY STATISTICS & LOGGING
    # =========================================================================

    num_errors <- sum(!is.na(row_data_df$error))
    .print_divergence_summary(num_genes, num_errors, elapsed, row_data_df, progress)

    # =========================================================================
    # CREATE & RETURN SUMMARIZED EXPERIMENT
    # =========================================================================

    result_se <- .construct_result_se(assay_matrix = assay_matrix, row_data_df = row_data_df,
        q_vals = q, elapsed = elapsed, nboot = nboot, ci = ci, method = method, norm = norm,
        use_parallel = use_parallel, num_genes = num_genes, num_errors = num_errors,
        bootstrap = bootstrap)

    return(result_se)
}

# ============================================================================
# HELPER FUNCTIONS FOR DIVERGENCE CALCULATION Extracted to meet Bioconductor
# ≤50 line requirement (March 2026)
# ============================================================================

# INPUT VALIDATION & CONFIGURATION HELPERS
# ============================================================================

#' Normalize norm parameter for backward compatibility
#' Coerces logical values to character strings
#' @noRd
.validate_norm_parameter <- function(norm) {
    if (is.logical(norm)) {
        norm <- if (norm)
            "range" else "none"
    }
    match.arg(norm, choices = c("none", "range", "zscore", "log_odds_ratio", "relative_reference"))
}

#' Validate and sort q-parameter values
#' Ensures q >= 0 and returns sorted vector
#' @noRd
.validate_and_sort_q_values <- function(q) {
    q <- sort(as.numeric(q))
    if (any(q < 0)) {
        stop("q parameter must be >= 0. ", "Note: q should be in range [0, 3] for typical use. ",
            "q=0 represents uniform divergence. ", "Got: ", paste(q, collapse = ", "))
    }
    q
}

#' Validate SummarizedExperiment input
#' @noRd
.validate_se_input <- function(se) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object")
    }
    TRUE
}

#' Identify gene name/ID column in rowData
#' Preference: gene_name (human-readable) > gene_id (ensembl)
#' @noRd
.identify_gene_column <- function(se) {
    rd <- SummarizedExperiment::rowData(se)
    gene_col_candidates <- c("gene_name", "gene_id")

    if (is.null(rd))
        return(NA_character_)

    for (col in gene_col_candidates) {
        if (col %in% colnames(rd)) {
            return(col)
        }
    }
    NA_character_
}

#' Extract unique gene list from SummarizedExperiment
#' @noRd
.extract_gene_list <- function(se, gene_col) {
    rd <- SummarizedExperiment::rowData(se)

    if (!is.na(gene_col) && !is.null(rd)) {
        all_genes <- unique(as.character(rd[[gene_col]]))
    } else {
        all_genes <- rownames(se)
    }

    if (length(all_genes) == 0) {
        stop("se must have gene identifiers in rowData or rownames")
    }

    all_genes
}

#' Configure parallel execution parameters
#' Auto-detects cores and decides between sequential/parallel
#' @noRd
.configure_parallel <- function(nthreads, num_genes) {
    if (is.null(nthreads) || is.na(nthreads)) {
        nthreads <- parallel::detectCores() - 1
        nthreads <- max(1, nthreads)
    }

    if (!is.numeric(nthreads) || nthreads < 1) {
        stop("'nthreads' must be a positive integer")
    }
    nthreads <- as.integer(nthreads)

    # Safe boolean check: only use parallel if num_genes is numeric and
    # nthreads > 1
    use_parallel <- (!is.na(num_genes) && num_genes >= 5 && nthreads > 1)

    list(nthreads = nthreads, use_parallel = use_parallel)
}

# GENE PROCESSING HELPERS
# ============================================================================

#' Aggregate transcript-level counts to gene-level
#' Sums counts across all transcripts for a given gene
#' @noRd
.compute_aggregate_counts <- function(se, target_gene, gene_col, rd) {
    # Find ALL transcripts for this gene
    if (!is.na(gene_col) && !is.null(rd)) {
        gene_transcript_indices <- which(as.character(rd[[gene_col]]) == target_gene)
    } else {
        gene_transcript_indices <- which(rownames(se) == target_gene)
    }

    if (length(gene_transcript_indices) == 0) {
        return(NULL)
    }

    # Aggregate counts across all transcripts for this gene
    counts_matrix <- as.matrix(SummarizedExperiment::assay(se, "counts")[gene_transcript_indices,
        , drop = FALSE])
    colSums(counts_matrix)
}

#' Extract group-specific counts
#' Splits gene-level counts by group membership
#' @noRd
.extract_group_counts_gene <- function(counts_gene, groups, control_group) {
    if (length(counts_gene) != length(groups)) {
        stop("Length mismatch: counts_gene and groups must have same length")
    }

    x <- counts_gene[groups == control_group]
    y <- counts_gene[groups != control_group]

    list(control = x, treatment = y)
}

#' Construct error result structure
#' Standardized format for failed gene computations
#' @noRd
.make_error_result <- function(gene_name, q_vals, error_msg, elapsed_sec = NA_real_) {
    list(gene_name = gene_name, results_per_q = rep(list(list(estimate = NA_real_,
        lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q_vals)),
        computation_time_sec = elapsed_sec, error = error_msg)
}

#' Build bootstrap arguments for calculate_divergence_bootstrap
#' Conditionally includes pair_ids if detected
#' @noRd

# NOTE (March 2026): .bootstrap_build_args() moved to bootstrap.R

#' @noRd
.compute_divergence_q <- function(x, y, q_vals, nboot, ci, method, log_base, pseudocount,
    gene_name, pair_ids = NULL) {
    gene_results <- list()

    # OPTIMIZATION (March 2026): Vectorize point estimate computation for
    # multi-q analysis Pre-compute all point estimates using vectorized
    # function (2-3x faster for 3+ q-values) Then use bootstrap for CI
    # computation separately
    if (length(q_vals) > 1) {
        # Vectorized point estimate computation (FAST PATH)
        point_estimates <- .tsallis_divergence_vector(x, y, q_vals, pseudocount = pseudocount,
            log_base = log_base)
    } else {
        point_estimates <- NULL  # Fall back to scalar computation for single q
    }

    # Process each q-value
    for (j in seq_along(q_vals)) {
        q_val <- q_vals[j]

        bootstrap_args <- .bootstrap_build_args(x, y, q_val, nboot, ci, method, log_base,
            pseudocount, gene_name, pair_ids)

        result <- do.call(.calculate_divergence_bootstrap, bootstrap_args)

        # If we have pre-computed point estimates and bootstrap was run, use
        # the vectorized point estimate (often more numerically stable)
        if (!is.null(point_estimates) && !is.na(point_estimates[j]) && nboot > 0) {
            # Verify consistency: vectorized vs scalar computation (should be
            # within 1e-8 relative error due to different computation order)
            result$estimate <- point_estimates[j]
        }

        gene_results[[j]] <- result
    }

    gene_results
}

# RESULTS COMPILATION HELPERS
# ============================================================================

#' Initialize result matrices for results compilation
#' Creates assay matrix and rowData structure
#' @noRd
.initialize_matrices <- function(num_genes, q_vals) {
    num_q_vals <- length(q_vals)

    assay_matrix <- matrix(NA_real_, nrow = num_genes, ncol = num_q_vals, dimnames = list(NULL,
        paste0("q_", q_vals)))

    row_data_df <- data.frame(gene_name = character(num_genes), error = character(num_genes),
        computation_time_sec = numeric(num_genes), stringsAsFactors = FALSE)

    # Add columns for each q value's metadata
    for (j in seq_len(num_q_vals)) {
        row_data_df[[paste0("estimate_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("lower_ci_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("upper_ci_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("ci_width_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("method_q", q_vals[j])]] <- NA_character_
        row_data_df[[paste0("nboot_q", q_vals[j])]] <- NA_integer_
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Populate result matrices from results_list
#' Extracts individual gene results and fills matrices
#' @noRd
.populate_matrices <- function(results_list, assay_matrix, row_data_df, q_vals) {
    num_q_vals <- length(q_vals)

    for (i in seq_along(results_list)) {
        result <- results_list[[i]]

        row_data_df$gene_name[i] <- result$gene_name
        row_data_df$error[i] <- if (is.na(result$error))
            NA_character_ else result$error
        row_data_df$computation_time_sec[i] <- result$computation_time_sec

        if (is.na(result$error)) {
            for (j in seq_len(num_q_vals)) {
                q_res <- result$results_per_q[[j]]
                assay_matrix[i, j] <- q_res$estimate

                row_data_df[[paste0("estimate_q", q_vals[j])]][i] <- q_res$estimate
                row_data_df[[paste0("lower_ci_q", q_vals[j])]][i] <- q_res$lower_ci
                row_data_df[[paste0("upper_ci_q", q_vals[j])]][i] <- q_res$upper_ci
                row_data_df[[paste0("method_q", q_vals[j])]][i] <- q_res$method %||%
                  NA_character_
                row_data_df[[paste0("nboot_q", q_vals[j])]][i] <- as.integer(q_res$nboot %||%
                  0)

                # Compute CI width
                if (!is.na(q_res$lower_ci) && !is.na(q_res$upper_ci)) {
                  row_data_df[[paste0("ci_width_q", q_vals[j])]][i] <- q_res$upper_ci -
                    q_res$lower_ci
                }
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

# NORMALIZATION HELPERS
# ============================================================================

#' Apply range normalization [0,1]
#' @noRd
.normalize_range_matrix <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]

        if (length(valid_vals) > 1) {
            min_val <- min(valid_vals)
            max_val <- max(valid_vals)
            range_val <- max_val - min_val

            if (range_val > 0) {
                assay_matrix[, j] <- (col_vals - min_val)/range_val

                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])

                if (estimate_col %in% colnames(row_data_df)) {
                  row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - min_val)/range_val
                  row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - min_val)/range_val
                  row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - min_val)/range_val
                }
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply z-score normalization
#' @noRd
.divergence_normalize_zscore <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]

        if (length(valid_vals) > 1) {
            mean_val <- mean(valid_vals)
            sd_val <- sd(valid_vals)

            if (sd_val > 0) {
                assay_matrix[, j] <- (col_vals - mean_val)/sd_val

                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])

                if (estimate_col %in% colnames(row_data_df)) {
                  row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - mean_val)/sd_val
                  row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - mean_val)/sd_val
                  row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - mean_val)/sd_val
                }
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply log odds ratio normalization
#' D_norm = log(D_q / D_max) where D_max depends on q
#' @noRd
.divergence_normalize_log_odds_ratio <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        q_val <- q_vals[j]
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]

        # Maximum divergence depends on q
        d_max <- if (q_val < 1)
            log(2) else 1

        if (d_max > 0) {
            assay_matrix[, j] <- log(pmax(col_vals, 1e-10)/d_max)

            # Apply same to estimate and CI bounds
            estimate_col <- paste0("estimate_q", q_vals[j])
            lower_col <- paste0("lower_ci_q", q_vals[j])
            upper_col <- paste0("upper_ci_q", q_vals[j])

            if (estimate_col %in% colnames(row_data_df)) {
                row_data_df[[estimate_col]] <- log(pmax(row_data_df[[estimate_col]],
                  1e-10)/d_max)
                row_data_df[[lower_col]] <- log(pmax(row_data_df[[lower_col]], 1e-10)/d_max)
                row_data_df[[upper_col]] <- log(pmax(row_data_df[[upper_col]], 1e-10)/d_max)
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply relative reference normalization
#' Ratio to reference group mean per q value
#' @noRd
.normalize_reference <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]

        if (length(valid_vals) > 0) {
            reference_mean <- mean(valid_vals, na.rm = TRUE)

            if (reference_mean > 0) {
                assay_matrix[, j] <- col_vals/reference_mean

                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])

                if (estimate_col %in% colnames(row_data_df)) {
                  row_data_df[[estimate_col]] <- row_data_df[[estimate_col]]/reference_mean
                  row_data_df[[lower_col]] <- row_data_df[[lower_col]]/reference_mean
                  row_data_df[[upper_col]] <- row_data_df[[upper_col]]/reference_mean
                }
            }
        }
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply normalization to results
#' Dispatcher that calls appropriate normalization method
#' @noRd
.normalize_divergence_matrix <- function(assay_matrix, row_data_df, q_vals, norm = c("none", "range", "zscore", "log_odds_ratio", "relative_reference")) {
    # Validate norm parameter per Bioconductor code syntax standards
    norm <- match.arg(norm)

    if (norm == "none") {
        return(list(assay = assay_matrix, rowData = row_data_df))
    }

    if (norm == "range") {
        return(.normalize_range_matrix(assay_matrix, row_data_df, q_vals))
    }

    if (norm == "zscore") {
        return(.divergence_normalize_zscore(assay_matrix, row_data_df, q_vals))
    }

    if (norm == "log_odds_ratio") {
        return(.divergence_normalize_log_odds_ratio(assay_matrix, row_data_df, q_vals))
    }

    if (norm == "relative_reference") {
        return(.normalize_reference(assay_matrix, row_data_df, q_vals))
    }

    list(assay = assay_matrix, rowData = row_data_df)
}

#' Generate computation summary statistics
#' Formats timing information for progress messages
#' @noRd
.generate_summary <- function(elapsed, num_genes, num_errors, row_data_df) {
    num_success <- num_genes - num_errors

    list(total_elapsed = elapsed, avg_per_gene = elapsed/num_genes, genes_per_minute = (num_genes/elapsed) *
        60, successful = num_success, failed = num_errors, failed_details = if (num_errors >
        0) {
        failed <- row_data_df[!is.na(row_data_df$error), ]
        failed[seq_len(min(10, nrow(failed))), c("gene_name", "error")]
    } else NULL)
}

#' Construct final SummarizedExperiment output
#' Builds SE with assays, rowData, colData, and metadata
#' @noRd
.construct_result_se <- function(assay_matrix, row_data_df, q_vals, elapsed, nboot,
    ci, method, norm, use_parallel, num_genes, num_errors, bootstrap = FALSE) {
    assays_list <- list(divergence = assay_matrix)

    # Extract bootstrap CI bounds if available (stored in rowData) When
    # bootstrap was used, CI bounds are in columns: lower_ci_q*, upper_ci_q*
    if (isTRUE(bootstrap) && (identical(nboot, "auto") || (is.numeric(nboot) && nboot > 0))) {
        # Initialize CI assay matrices
        ci_lower_matrix <- matrix(NA_real_, nrow = nrow(assay_matrix), ncol = ncol(assay_matrix),
            dimnames = dimnames(assay_matrix))
        ci_upper_matrix <- matrix(NA_real_, nrow = nrow(assay_matrix), ncol = ncol(assay_matrix),
            dimnames = dimnames(assay_matrix))

        # Extract CI bounds from rowData for each q-value
        for (j in seq_along(q_vals)) {
            q_val <- q_vals[j]
            lower_col <- paste0("lower_ci_q", q_val)
            upper_col <- paste0("upper_ci_q", q_val)

            # Check if these columns exist in rowData
            if (lower_col %in% colnames(row_data_df) && upper_col %in% colnames(row_data_df)) {
                ci_lower_matrix[, j] <- row_data_df[[lower_col]]
                ci_upper_matrix[, j] <- row_data_df[[upper_col]]
            }
        }

        # Add CI assays if any values were extracted
        if (!all(is.na(ci_lower_matrix))) {
            assays_list$ci_lower <- ci_lower_matrix
        }
        if (!all(is.na(ci_upper_matrix))) {
            assays_list$ci_upper <- ci_upper_matrix
        }
    }

    col_data_output <- data.frame(q_value = q_vals, sample_type = rep("divergence_estimate",
        length(q_vals)), computation_mode = rep(if (use_parallel) "parallel" else "sequential",
        length(q_vals)), row.names = paste0("q_", q_vals))

    num_success <- num_genes - num_errors

    SummarizedExperiment::SummarizedExperiment(assays = assays_list, rowData = row_data_df,
        colData = col_data_output, metadata = list(summary_stats = list(total_genes = num_genes,
            successful = num_success, failed = num_errors), elapsed_time_sec = elapsed,
            avg_time_per_gene = elapsed/num_genes, genes_per_minute = (num_genes/elapsed) *
                60, bootstrap_config = list(nboot = nboot, ci = ci, method = method),
            normalization = norm, computation_mode = if (use_parallel) "parallel" else "sequential"))
}

#' Process a single gene for divergence computation
#' Consolidates logic shared between sequential and parallel processing
#' @noRd
.process_single_gene_div <- function(gene_idx, all_gene_names, se, gene_col, rd,
    group_col, control_group, q, nboot, ci, method, log_base, pseudocount, pair_ids) {
    target_gene <- all_gene_names[gene_idx]
    gene_name <- target_gene
    gene_start <- Sys.time()

    tryCatch({
        # Aggregate transcript-level counts to gene-level
        counts_gene <- .compute_aggregate_counts(se, target_gene, gene_col, rd)

        if (is.null(counts_gene) || length(counts_gene) == 0) {
            return(.make_error_result(gene_name, q, "No transcripts found for gene"))
        }

        # Extract group vector for sample grouping
        groups <- se[[group_col]]

        group_counts <- .extract_group_counts_gene(counts_gene, groups, control_group)
        x <- group_counts$control
        y <- group_counts$treatment

        if (length(x) == 0 || length(y) == 0) {
            return(.make_error_result(gene_name, q, "Insufficient group samples"))
        }

        # Compute divergence for each q value (with vectorization optimization)
        gene_results <- .compute_divergence_q(x, y, q, nboot, ci, method, log_base,
            pseudocount, gene_name, pair_ids)

        gene_elapsed <- as.numeric(Sys.time() - gene_start, units = "secs")

        list(gene_name = gene_name, results_per_q = gene_results, computation_time_sec = gene_elapsed,
            error = NA_character_)
    }, error = function(e) {
        .make_error_result(gene_name, q, as.character(e$message), as.numeric(Sys.time() -
            gene_start, units = "secs"))
    })
}
