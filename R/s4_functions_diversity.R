#' Calculate diversity and store in TSENATAnalysis
#'
#' Wrapper around .calculate_diversity() that manages TSENATAnalysis object.
#' Calculates Tsallis entropy (diversity) across multiple q-values for each gene
#' to quantify isoform complexity and transcript heterogeneity.
#'
#' Key Features:
#' \itemize{
#'   \item Multi-q analysis: Entropy computed for q = 0.01 to 2.00 (by default)
#'   \item Normalized entropy: Scale to [0, 1] using theoretical maximum log(m)
#'   \item Bootstrap confidence intervals: Quantify uncertainty in estimates
#'   \item TPM normalization: Optional SALMON TPM-based weighting
#'   \item Multiple normalization methods: Range, Z-score, log-odds-ratio, relative
#'   \item Spectrum computation: Aggregate diversity statistics across all q-values
#' }
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for Tsallis entropy.
#'   If NULL,  uses q_values from \code{analysis@config$q_values} if 
#' available,  else defaults to seq(0. 01,  2,  by = 0. 05).
#' @param norm \code{logical} or  \code{character}.  Normalization method:
#'  TRUE,  FALSE,  'none',  'range',  'zscore',  'log_odds_ratio',
#'  'relative_reference'.
#'   If NULL, reads from \code{@config$norm} or defaults to TRUE.
#' @param norm_method \code{character}.
#'  Post-hoc normalization method applied after diversity computation.
#'   Options:
#'   \itemize{
#'     \item \code{'default'} - Simple normalization by theoretical maximum (current behavior)
#'     \item \code{'zscore'} - Z-score normalization per q-value
#'     \item \code{'log_odds_ratio'} - Log-odds ratio relative to max entropy (q and 
#' isoform-aware)
#'     \item \code{'relative_reference'} - Divide by reference group mean (requires reference_group)
#'     \item \code{NULL} - No post-hoc normalization (default)
#'   }
#'   If NULL, reads from \code{@config$norm_method} if available.
#' @param tpm \code{logical}. TPM normalization. Default: FALSE.
#'   If not specified, reads from \code{@config$tpm} if available.
#' @param assayno \code{numeric}. Assay number to use. Default: 1.
#'   If NULL, reads from \code{@config$assayno} if available.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#'   If not specified, reads from \code{@config$verbose} if available.
#' @param what \code{character}.  Output type:  'S' (entropy) or 
#' 'D' (diversity).  Default:  'S'.
#'   If NULL, reads from \code{@config$what} if available.
#' @param nthreads \code{numeric} or  \code{NULL}.  Number of CPU threads for 
#' parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param pseudocount \code{numeric} or  \code{character}.
#'  Pseudocount value or  'auto'.  Default:  0.
#'   If NULL, reads from \code{@config$pseudocount} if available.
#' @param show_messages \code{logical}.
#'  Display verbose messages during computation.  Default:  FALSE.
#' @param min_valid_frac \code{numeric}.
#'  Minimum fraction of valid samples for  gene filtering.  Default:  NULL.
#'   If NULL, reads from \code{@config$min_valid_frac} if available.
#' @param shrinkage \code{character}.  Shrinkage method:  'none' or 
#' 'empirical_bayes'.  Default:  'none'.
#'   If NULL, reads from \code{@config$shrinkage} if available.
#' @param genes \code{character} or  \code{NULL}.  Gene set specification.
#'  Default:  NULL (use all genes).
#'   If NULL, reads from \code{@config$genes} if available.
#' @param effective_length \code{numeric} or  \code{NULL}.
#'  Effective gene lengths.  Default:  NULL.
#'   If NULL, reads from \code{@config$effective_length} if available.
#' @param metadata \code{list} or  \code{NULL}.  Additional metadata.
#'  Default:  NULL.
#' @param bootstrap \code{logical}.  Compute bootstrap confidence intervals.
#'  Default:  FALSE.
#'   If not specified, reads from \code{@config$bootstrap} if available.
#' @param nboot \code{numeric} or  \code{NULL}.
#'  Number of bootstrap replicates.  Default:  NULL.
#'   If NULL, reads from \code{@config$nboot} if available.
#' @param bootstrap_method \code{character}.  Bootstrap method:
#'  'percentile' or  others.  Default:  'percentile'.
#'   If NULL, reads from \code{@config$bootstrap_method} if available.
#' @param bootstrap_ci \code{numeric}.
#'  Bootstrap confidence interval level (0-1).  Default:  0. 95.
#'   If NULL, reads from \code{@config$bootstrap_ci} if available.
#' @param bootstrap_include_diagnostics \code{logical}.
#'  Include bootstrap diagnostic information.  Default:  FALSE.
#'   If NULL,  reads from \code{@config$bootstrap_include_diagnostics} if 
#' available.
#' @details
#' **Mathematical Background:**
#' Tsallis entropy H_q for q-parameter:
#' \preformatted{
#'   H_q(X) = (1/(q-1)) * (1 - Sum p_i^q)  [for q != 1]
#'   H_1(X) = -Sum p_i * log(p_i)  [Shannon entropy, limit q→1]
#' }
#' where p_i = relative abundance of isoform i for a gene.
#' Larger q emphasizes dominant isoforms; smaller q emphasizes rare ones.
#'
#' **Example Interpretation:**
#' \itemize{
#'   \item Gene with 1 isoform: H_q = 0 for all q (no diversity)
#'   \item Gene with 2 equal isoforms: H_q ~ 0.5-1.0 depending on q
#'   \item Gene with m equally abundant isoforms: H_q = 1.0 (maximum diversity)
#' }
#'
#' @param reference_group \code{character}.
#'  For \code{norm_method = 'relative_reference'},  
#' the reference group column name (e.g., from colData). If NULL, uses first
#' group in colData.
#'   If NULL, reads from \code{@config$reference_group} if available.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#'   When provided, generates TWO files:
#'   \enumerate{
#'     \item Primary output: Analysis object (.rds) or table (.tsv/.csv/.txt)
#'     \item Secondary output: Diversity spectrum statistics (TSV format)
#'         with suffix \code{_spectrum.tsv}
#'   }
#'   Example: output_file = 'analysis.rds' generates:
#'   \itemize{
#'     \item \code{analysis.rds} - TSENATAnalysis object
#'     \item \code{analysis_spectrum.tsv} - Spectrum statistics
#'   }
#'   The spectrum file contains columns: q, central (median diversity), 
#'   spread (IQR), count, and group (if grouping variable available).
#'   Default: NULL (no file output).
#' @param ... Additional arguments passed to underlying functions for
#' extensibility.
#'
#' @return Modified TSENATAnalysis object with diversity results stored
#'   in \code{@diversity_results},  keyed by 'q_X. XX. . . ' format (e. g. ,
#'  'q_1. 000').
#'   When \code{output_file} is provided, also generates:
#'   \itemize{
#'     \item Primary file: Analysis object or table export
#'     \item Spectrum file: \code{*_spectrum.tsv} containing 
#'           aggregated diversity statistics across q-values and groups
#'   }
#'
#' @details
#' This wrapper calls \code{.calculate_diversity()} once per q-value, storing
#' results as SummarizedExperiment objects. It extracts key parameters from
#' \code{analysis@config} with 
#' priority resolution (explicit > \code{@config} > default).
#'
#' **Diversity Spectrum Computation:**
#' By default, this function computes and saves a diversity spectrum (aggregated
#' statistics across all q-values and  groups) when 
#' \code{output_file} is provided.
#' The spectrum contains:
#' \itemize{
#'   \item \code{q}: Diversity parameter value
#'   \item \code{central}: Median (or mean) diversity across all genes
#'   \item \code{spread}: IQR (or SD) around central value
#'   \item \code{count}: Number of valid measurements
#'   \item \code{group}: Condition group (if applicable)
#' }
#' This provides a quick summary of how diversity changes across q-values,
#' useful for q-curve visualization and statistical comparisons.
#' Spectrum is saved as: \code{*_spectrum.tsv}
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{q}{Priority 1 (explicit) > Priority 2 (\code{@config$q_values}) > Priority 3 (default:
#'  seq(0. 01,  2,  by = 0. 05))\cr
#'     **Note: ** If explicit q AND \code{@config$q_values} both provided,
#'  explicit wins. }
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
#'   \item{verbose}{Priority: explicit > \code{@config$verbose} > TRUE}
#'   \item{bootstrap}{Priority: explicit > \code{@config$bootstrap} > FALSE}
#'   \item{pseudocount}{Priority: explicit > \code{@config$pseudocount} > 0}
#'   \item{norm}{Priority: explicit > \code{@config$norm} > TRUE}
#'   \item{what}{Priority:
#'  explicit > \code{@config$what} > 'S' (Tsallis entropy)}
#' }
#'
#' **Audit Trail:** After execution, check:
#' \itemize{
#'   \item \code{analysis@config$last_diversity_run$parameters_used}:
#'  Actual parameters 
#'         (not original \code{@config})
#'   \item \code{attr(diversity(analysis,  q=X),  'computed_with')}:
#'  Per-q-value metadata
#'         (timestamp, bootstrap setting, nthreads, etc.)
#' }
#'
#' @examples
#' # Load vignette data and build analysis
#' data(readcounts)
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' 
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df,
#'   tpm = tpm, effective_length = effective_length)
#' 
#' # Filter to manageable size (use 200+ genes to survive diversity filtering)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Compute diversity and access results
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0), verbose =
#' FALSE)
#' head(diversity(analysis, q = 1.0))
#'
#' @details
#' For additional details on diversity spectrum calculations and
#' normalization methods,
#' see the package vignettes.
#'
#' @export
#' @importFrom utils write.table
calculate_diversity_s4 <- function(analysis, q = NULL, norm = TRUE, norm_method = NULL,
    reference_group = NULL, tpm = FALSE, assayno = NULL, verbose = NULL, show_messages = FALSE,
    what = NULL, nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL, shrinkage = NULL,
    genes = NULL, effective_length = NULL, metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL, bootstrap_include_diagnostics = NULL,
    output_file = NULL, ...) {
    # Validate input
    .validate_diversity_analysis_input(analysis)

    # Prepare parameters and build calculation
    params <- .prepare_diversity_params(analysis, q, norm, norm_method, reference_group,
        tpm, assayno, verbose, what, nthreads, pseudocount, min_valid_frac, shrinkage,
        genes, effective_length, metadata, bootstrap, nboot, bootstrap_method, bootstrap_ci,
        bootstrap_include_diagnostics, show_messages)
    .validate_norm_method(params$norm_method)

    # Execute diversity calculation
    calc_args <- .build_calc_diversity_args(params, analysis, list(...))
    result_df <- do.call(.calculate_diversity, calc_args)

    # Resolve pseudocount and cache results
    pseudocount_resolved <- .resolve_pseudocount(params, analysis, result_df)
    analysis <- .cache_combined_diversity_result(analysis, result_df, params, pseudocount_resolved)

    # Process q-values and update configuration
    analysis <- .process_all_q_values(analysis, result_df, params, pseudocount_resolved)
    analysis <- .update_config_post_diversity(analysis, params, pseudocount_resolved)

    # Handle spectrum and output files
    analysis <- .handle_spectrum_and_output(analysis, params, output_file)

    analysis
}

# ============================================================================
# HELPER: Validate analysis input object
# ============================================================================
#' @noRd
.validate_diversity_analysis_input <- function(analysis) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }
    if (nrow(analysis@se) == 0) {
        stop("SummarizedExperiment in @se is empty", call. = FALSE)
    }
}

# ============================================================================
# HELPER: Validate norm_method parameter
# ============================================================================
#' @noRd
.validate_norm_method <- function(norm_method) {
    if (!is.null(norm_method)) {
        valid_methods <- c("default", "zscore", "log_odds_ratio", "relative_reference")
        if (!(norm_method %in% valid_methods)) {
            stop("'norm_method' must be one of: ", paste(valid_methods, collapse = ", "),
                call. = FALSE)
        }
    }
}

# ============================================================================
# HELPER: Resolve pseudocount value
# ============================================================================
#' @noRd
.resolve_pseudocount <- function(params, analysis, result_df) {
    pseudocount_resolved <- params$pseudocount
    if (is.character(params$pseudocount) && tolower(params$pseudocount) == "auto") {
        pseudocount_resolved <- .handle_pseudocount_auto(params$pseudocount, analysis@se,
            verbose = FALSE)
        if (params$verbose) {
            message(sprintf("[calculate_diversity_s4] Resolved pseudocount (auto) = %.4f",
                pseudocount_resolved))
        }
    }
    pseudocount_resolved
}

# ============================================================================
# HELPER: Cache combined diversity result
# ============================================================================
#' @noRd
.cache_combined_diversity_result <- function(analysis, result_df, params, pseudocount_resolved) {
    analysis@metadata$diversity_combined <- list(combined_result = result_df, combined_se = result_df,
        q_values_computed = params$q, computation_params = list(norm = params$norm,
            verbose = params$verbose, bootstrap = params$bootstrap, pseudocount = pseudocount_resolved,
            pseudocount_original = params$pseudocount, nthreads = params$nthreads,
            what = params$what), timestamp = Sys.time())
    analysis@metadata$.diversity_lazily_computed <- TRUE
    analysis
}

# ============================================================================
# HELPER: Process all q-values and store results
# ============================================================================
#' @noRd
.process_all_q_values <- function(analysis, result_df, params, pseudocount_resolved) {
    q_decimals <- 3
    col_q_values <- .extract_q_metadata_from_result(result_df, params$q)
    result_se_original <- result_df
    q_format_cache <- .build_q_format_cache(params$q, q_decimals)

    for (q_val in params$q) {
        analysis <- .process_single_q_value(analysis, result_df, col_q_values, result_se_original,
            q_format_cache, q_val, q_decimals, params, pseudocount_resolved)
    }
    analysis
}

# ============================================================================
# HELPER: Build q-value formatting cache
# ============================================================================
#' @noRd
.build_q_format_cache <- function(q_values, q_decimals) {
    q_format_cache <- list()
    for (q_val in q_values) {
        q_formatted <- formatC(q_val, format = "f", digits = q_decimals)
        q_pattern <- paste0("_q=", gsub("\\.", "\\\\.", q_formatted), "$")
        q_format_cache[[as.character(q_val)]] <- list(formatted = q_formatted, pattern = q_pattern,
            patterns_alt = list(paste0("_q=", gsub("\\.", "\\\\.", formatC(q_val,
                format = "f", digits = 3)), "$"), paste0("_q=", gsub("\\.", "\\\\.",
                formatC(q_val, format = "f", digits = 2)), "$"), paste0("_q=", gsub("\\.",
                "\\\\.", formatC(q_val, format = "f", digits = 1)), "$"), paste0("_q=",
                gsub("\\.", "\\\\.", as.character(as.integer(q_val))), "$")))
    }
    q_format_cache
}

# ============================================================================
# HELPER: Process a single q-value
# ============================================================================
#' @noRd
.process_single_q_value <- function(analysis, result_df, col_q_values, result_se_original,
    q_format_cache, q_val, q_decimals, params, pseudocount_resolved) {
    tryCatch({
        # Extract columns for this q-value
        q_cols <- .extract_q_columns(result_df, col_q_values, q_val, q_format_cache)

        # Extract and convert to SE, preserving CI assays
        result_subset <- result_df[, q_cols, drop = FALSE]
        result_se <- .convert_result_to_se(result_subset, result_se_original, q_cols,
            q_val)

        # Validate and apply colData
        .validate_se_structure(result_se, q_val)
        result_se <- .apply_original_coldata(result_se, analysis@se)

        # Apply post-hoc normalization
        result_se <- .apply_diversity_post_hoc_norm(result_se, params$norm_method,
            params, q_val, params$verbose)

        # Store with audit trail
        key <- paste0("q_", formatC(q_val, format = "f", digits = q_decimals))
        attr(result_se, "computed_with") <- list(q = q_val, norm = params$norm, norm_method = params$norm_method,
            verbose = params$verbose, bootstrap = params$bootstrap, pseudocount = pseudocount_resolved,
            pseudocount_original = params$pseudocount, nthreads = params$nthreads,
            what = params$what, timestamp = Sys.time())
        analysis@diversity_results[[key]] <- result_se
        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_diversity[q=",
            q_val, "]"))
    }, error = function(e) {
        error_msg <- conditionMessage(e)
        if (params$bootstrap && grepl("bootstrap", error_msg, ignore.case = TRUE)) {
            stop("[calculate_diversity_s4] Bootstrap CI computation failed for q=",
                q_val, ": ", error_msg, call. = FALSE)
        } else {
            stop("[calculate_diversity_s4] Failed to compute diversity for q=", q_val,
                ":\n", error_msg, call. = FALSE)
        }
    })
    analysis
}

# ============================================================================
# HELPER: Extract columns for a specific q-value
# ============================================================================
#' @noRd
.extract_q_columns <- function(result_df, col_q_values, q_val, q_format_cache) {
    q_cols <- NULL
    if (!is.na(col_q_values[1])) {
        q_mask <- abs(col_q_values - q_val) < 1e-05
        if (any(q_mask)) {
            q_cols <- as.numeric(names(col_q_values)[q_mask])
        }
    }

    if (is.null(q_cols) || length(q_cols) == 0) {
        q_cache_entry <- q_format_cache[[as.character(q_val)]]
        for (pattern_idx in seq_along(q_cache_entry$patterns_alt)) {
            q_pattern_str <- q_cache_entry$patterns_alt[[pattern_idx]]
            q_cols_temp <- grep(q_pattern_str, colnames(result_df))
            if (length(q_cols_temp) > 0) {
                q_cols <- q_cols_temp
                break
            }
        }
    }

    if (length(q_cols) == 0) {
        stop("[calculate_diversity_s4] No columns found for q=", q_val, call. = FALSE)
    }
    q_cols
}

# ============================================================================
# HELPER: Convert result subset to SummarizedExperiment
# ============================================================================
#' @noRd
.convert_result_to_se <- function(result_subset, result_se_original, q_cols, q_val) {
    if (is(result_subset, "SummarizedExperiment")) {
        return(result_subset)
    }

    if (!is.data.frame(result_subset)) {
        return(result_subset)
    }

    numeric_cols <- vapply(result_subset, is.numeric, FUN.VALUE = logical(1))
    if (!any(numeric_cols)) {
        return(result_subset)
    }

    assay_data <- as.matrix(result_subset[, numeric_cols, drop = FALSE])
    assays_list <- list(diversity = assay_data)

    # Extract CI assays from original result if available
    if (is(result_se_original, "SummarizedExperiment")) {
        if ("ci_lower" %in% SummarizedExperiment::assayNames(result_se_original)) {
            ci_lower_orig <- SummarizedExperiment::assay(result_se_original, "ci_lower")
            if (!is.null(ci_lower_orig)) {
                assays_list$ci_lower <- ci_lower_orig[, q_cols, drop = FALSE]
            }
        }
        if ("ci_upper" %in% SummarizedExperiment::assayNames(result_se_original)) {
            ci_upper_orig <- SummarizedExperiment::assay(result_se_original, "ci_upper")
            if (!is.null(ci_upper_orig)) {
                assays_list$ci_upper <- ci_upper_orig[, q_cols, drop = FALSE]
            }
        }
    }

    result_se <- SummarizedExperiment(assays = assays_list)
    rownames(result_se) <- rownames(result_subset)

    metadata_mask <- !numeric_cols
    if (any(metadata_mask)) {
        cd <- result_subset[, metadata_mask, drop = FALSE]
        rownames(cd) <- colnames(assay_data)
        SummarizedExperiment::colData(result_se) <- cd
    }

    result_se
}

# ============================================================================
# HELPER: Validate SummarizedExperiment structure
# ============================================================================
#' @noRd
.validate_se_structure <- function(result_se, q_val) {
    if (!is(result_se, "SummarizedExperiment")) {
        return(TRUE)
    }

    if (length(SummarizedExperiment::assays(result_se)) == 0) {
        stop("[calculate_diversity_s4] Converted SE for q=", q_val, " has no assays",
            call. = FALSE)
    }

    test_assay <- tryCatch({
        SummarizedExperiment::assay(result_se, 1)
    }, error = function(e) {
        stop("[calculate_diversity_s4] Cannot access assay for q=", q_val, ": ",
            conditionMessage(e), call. = FALSE)
    })

    if (is.null(test_assay) || nrow(test_assay) == 0) {
        warning("[calculate_diversity_s4] Assay for q=", q_val, " is empty", call. = FALSE)
    }
    TRUE
}

# ============================================================================
# HELPER: Apply original colData to result SE
# ============================================================================
#' @noRd
.apply_original_coldata <- function(result_se, original_se) {
    if (!is(result_se, "SummarizedExperiment") || ncol(result_se) == 0) {
        return(result_se)
    }

    original_coldata <- SummarizedExperiment::colData(original_se)
    if (!is.null(original_coldata) && nrow(original_coldata) == ncol(result_se)) {
        SummarizedExperiment::colData(result_se) <- original_coldata
    }
    result_se
}

# ============================================================================
# HELPER: Update config post-diversity computation
# ============================================================================
#' @noRd
.update_config_post_diversity <- function(analysis, params, pseudocount_resolved) {
    analysis@config$last_diversity_run <- list(timestamp = Sys.time(), q_values_computed = params$q,
        num_q_values = length(params$q), parameters_used = list(norm = params$norm,
            norm_method = params$norm_method, verbose = params$verbose, bootstrap = params$bootstrap,
            pseudocount = pseudocount_resolved, pseudocount_original = params$pseudocount,
            nthreads = params$nthreads, what = params$what), note = "Actual parameters used (save object and check this)")

    if (params$nthreads > 1) {
        analysis@metadata$parallel_processing <- c(analysis@metadata$parallel_processing,
            paste0("calculate_diversity_s4: nthreads=", params$nthreads, " (", length(params$q),
                " q-values)"))
    }
    analysis
}

# ============================================================================
# HELPER: Handle spectrum computation and file output
# ============================================================================
#' @noRd
.handle_spectrum_and_output <- function(analysis, params, output_file) {
    if (is.null(output_file)) {
        return(analysis)
    }

    .compute_and_save_spectrum(analysis, params, output_file)
    .save_diversity_output(analysis, params, output_file)
    analysis
}

# ============================================================================
# HELPER: Compute and save diversity spectrum
# ============================================================================
#' @noRd
.compute_and_save_spectrum <- function(analysis, params, output_file) {
    tryCatch({
        spectrum_condition_col <- analysis@config$condition_col %||% "sample_type"
        combined_se <- analysis@metadata$diversity_combined$combined_se

        if (is.null(combined_se) || nrow(combined_se) == 0) {
            return(analysis)
        }

        diversity_spectrum <- .compute_diversity_spectrum(se = combined_se, metric = "median",
            variability_metric = "iqr", condition_col = spectrum_condition_col)

        if (is.null(diversity_spectrum) || nrow(diversity_spectrum) == 0) {
            return(analysis)
        }

        spectrum_file <- sub("\\.[^.]+$", "_spectrum.tsv", output_file)
        if (spectrum_file == output_file) {
            spectrum_file <- paste0(output_file, "_spectrum.tsv")
        }

        utils::write.table(diversity_spectrum, file = spectrum_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

        if (params$verbose) {
            message("[calculate_diversity_s4] Saved diversity spectrum to: ", spectrum_file)
        }
        analysis@metadata$diversity_spectrum <- diversity_spectrum
    }, error = function(e) {
        warning("[calculate_diversity_s4] Could not compute diversity spectrum: ",
            conditionMessage(e), call. = FALSE)
    })
    analysis
}

# ============================================================================
# HELPER: Save diversity output to file
# ============================================================================
#' @noRd
.save_diversity_output <- function(analysis, params, output_file) {
    tryCatch({
        output_data <- .build_diversity_output_table(analysis, params)

        if (is.null(output_data) || nrow(output_data) == 0) {
            return(analysis)
        }

        save_analysis_output(output_data, output_file, verbose = params$verbose,
            func_name = "calculate_diversity_s4")

        if (params$verbose) {
            message("[calculate_diversity_s4] Saved diversity results to: ", output_file)
        }
    }, error = function(e) {
        warning("[calculate_diversity_s4] Could not save diversity results: ", conditionMessage(e),
            call. = FALSE)
    })
    analysis
}

# ============================================================================
# HELPER: Build diversity output table from results
# ============================================================================
#' @noRd
.build_diversity_output_table <- function(analysis, params) {
    output_data <- NULL

    if (length(analysis@diversity_results) > 0) {
        output_data <- .build_output_from_diversity_results(analysis, params)
    }

    if (is.null(output_data) || nrow(output_data) == 0) {
        output_data <- .build_output_from_combined_result(analysis)
    }

    output_data
}

# ============================================================================
# HELPER: Build output from individual diversity result SEs
# ============================================================================
#' @noRd
.build_output_from_diversity_results <- function(analysis, params) {
    q_keys_to_use <- paste0("q_", formatC(params$q, format = "f", digits = 3))

    # Get CI flags from first result
    first_key <- q_keys_to_use[1]
    if (!(first_key %in% names(analysis@diversity_results))) {
        first_key <- names(analysis@diversity_results)[1]
    }

    se_first <- analysis@diversity_results[[first_key]]
    has_ci_lower <- "ci_lower" %in% SummarizedExperiment::assayNames(se_first)
    has_ci_upper <- "ci_upper" %in% SummarizedExperiment::assayNames(se_first)

    all_data_list <- list()
    for (q_key in q_keys_to_use) {
        se <- analysis@diversity_results[[q_key]]
        if (!is(se, "SummarizedExperiment"))
            next

        current_data <- .extract_se_to_dataframe(se, q_key, has_ci_lower, has_ci_upper)
        all_data_list[[length(all_data_list) + 1]] <- current_data
    }

    if (length(all_data_list) > 0) {
        do.call(rbind, all_data_list)
    } else {
        NULL
    }
}

# ============================================================================
# HELPER: Extract SE to dataframe format for output
# ============================================================================
#' @noRd
.extract_se_to_dataframe <- function(se, q_key, has_ci_lower, has_ci_upper) {
    diversity_mat <- as.matrix(SummarizedExperiment::assay(se, 1))

    ci_lower_mat <- if (has_ci_lower) {
        tryCatch(as.matrix(SummarizedExperiment::assay(se, "ci_lower")), error = function(e) NULL)
    } else NULL

    ci_upper_mat <- if (has_ci_upper) {
        tryCatch(as.matrix(SummarizedExperiment::assay(se, "ci_upper")), error = function(e) NULL)
    } else NULL

    gene_names <- rownames(diversity_mat)
    sample_names <- colnames(diversity_mat)
    n_genes <- length(gene_names)

    current_data <- data.frame(gene = rep(gene_names, times = length(sample_names)),
        sample = rep(sample_names, each = n_genes), q_value = q_key, diversity = as.vector(diversity_mat),
        stringsAsFactors = FALSE)

    if (has_ci_lower && !is.null(ci_lower_mat)) {
        current_data$ci_lower <- as.vector(ci_lower_mat)
    }
    if (has_ci_upper && !is.null(ci_upper_mat)) {
        current_data$ci_upper <- as.vector(ci_upper_mat)
    }

    current_data
}

# ============================================================================
# HELPER: Build output from combined diversity result
# ============================================================================
#' @noRd
.build_output_from_combined_result <- function(analysis) {
    combined <- analysis@metadata$diversity_combined$combined_result
    if (is.null(combined) || nrow(combined) == 0) {
        return(NULL)
    }

    if (is(combined, "SummarizedExperiment")) {
        as.data.frame(SummarizedExperiment::assay(combined, 1))
    } else if (is.data.frame(combined)) {
        combined
    } else if (is.matrix(combined)) {
        as.data.frame(combined)
    } else {
        NULL
    }
}

# ============================================================================
# HELPER: Prepare and resolve all parameters for diversity calculation
# ============================================================================
#' @noRd
.prepare_diversity_params <- function(analysis, q = NULL, norm = NULL, norm_method = NULL,
    reference_group = NULL, tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL, shrinkage = NULL,
    genes = NULL, effective_length = NULL, metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL, bootstrap_include_diagnostics = NULL,
    show_messages = FALSE) {
    # Extract q parameter with default range
    if (is.null(q)) {
        q <- if ("q_values" %in% names(analysis@config)) {
            analysis@config$q_values
        } else {
            seq(0.01, 2, by = 0.05)
        }
    }

    if (!is.numeric(q)) {
        stop("'q' must be numeric", call. = FALSE)
    }

    # Validate q values are finite
    if (any(!is.finite(q))) {
        stop("'q' values must be finite (not Inf or NaN)", call. = FALSE)
    }

    # Resolve all parameters using centralized handler
    nthreads_resolved <- resolve_slot_param(nthreads, analysis@config, "nthreads",
        1)

    # Validate nthreads is positive
    if (!is.numeric(nthreads_resolved) || nthreads_resolved < 1) {
        stop("'nthreads' must be a positive integer", call. = FALSE)
    }

    # Ensure tpm is logical
    if (is.null(tpm)) {
        # If tpm is NULL, check config
        tpm <- if ("tpm" %in% names(analysis@config)) {
            as.logical(analysis@config$tpm)
        } else {
            FALSE  # Default to FALSE
        }
    } else {
        # If tpm is provided, coerce to logical
        tpm <- as.logical(tpm)
    }

    list(q = q, nthreads = nthreads_resolved, verbose = resolve_slot_param(verbose,
        analysis@config, "verbose", TRUE), show_messages = show_messages, bootstrap = resolve_slot_param(bootstrap,
        analysis@config, "bootstrap", FALSE), pseudocount = resolve_slot_param(pseudocount,
        analysis@config, "pseudocount", 0), min_valid_frac = resolve_slot_param(min_valid_frac,
        analysis@config, "min_valid_frac", 0.75), norm = resolve_slot_param(norm,
        analysis@config, "norm", TRUE), what = resolve_slot_param(what, analysis@config,
        "what", "S"), assayno = resolve_slot_param(assayno, analysis@config, "assayno",
        1), shrinkage = resolve_slot_param(shrinkage, analysis@config, "shrinkage",
        "none"), bootstrap_method = resolve_slot_param(bootstrap_method, analysis@config,
        "bootstrap_method", "percentile"), bootstrap_ci = resolve_slot_param(bootstrap_ci,
        analysis@config, "bootstrap_ci", 0.95), tpm = tpm, genes = resolve_slot_param(genes,
        analysis@config, "genes", NULL), effective_length = resolve_slot_param(effective_length,
        analysis@config, "effective_length", NULL), nboot = resolve_slot_param(nboot,
        analysis@config, "nboot", NULL), bootstrap_include_diagnostics = resolve_slot_param(bootstrap_include_diagnostics,
        analysis@config, "bootstrap_include_diagnostics", TRUE), metadata = resolve_slot_param(metadata,
        analysis@config, "metadata", NULL), norm_method = resolve_slot_param(norm_method,
        analysis@config, "norm_method", NULL), reference_group = resolve_slot_param(reference_group,
        analysis@config, "reference_group", NULL))
}

# ============================================================================
# HELPER: Build calculation arguments for .calculate_diversity()
# ============================================================================
#' @noRd
.build_calc_diversity_args <- function(params, analysis, dots) {
    calc_args <- list(x = analysis@se, q = params$q, norm = params$norm, tpm = params$tpm,
        assayno = params$assayno, verbose = params$verbose, show_messages = params$show_messages,
        what = params$what, nthreads = params$nthreads, pseudocount = params$pseudocount,
        min_valid_frac = params$min_valid_frac, shrinkage = params$shrinkage, bootstrap = params$bootstrap,
        bootstrap_nboot = params$nboot, bootstrap_method = params$bootstrap_method,
        bootstrap_ci = params$bootstrap_ci, bootstrap_include_diagnostics = params$bootstrap_include_diagnostics)

    # Add optional parameters
    if (!is.null(params$genes)) {
        calc_args$genes <- params$genes
    }
    if (!is.null(params$effective_length)) {
        calc_args$effective_length <- params$effective_length
    }
    if (!is.null(params$metadata)) {
        calc_args$metadata <- params$metadata
    }
    
    # Extract and pass column mapping parameters from config (required for metadata mapping)
    config <- getConfig(analysis)
    calc_args$sample_col <- config$sample_col
    calc_args$condition_col <- config$condition_col
    calc_args$subject_col <- config$subject_col

    # Add any additional parameters from dots
    c(calc_args, dots)
}

# ============================================================================
# HELPER: Extract q-value metadata from result dataframe
# ============================================================================
#' @noRd
.extract_q_metadata_from_result <- function(result_df, q_values) {
    col_q_values <- NA

    if (length(q_values) > 1 && nrow(result_df) > 0) {
        col_names <- colnames(result_df)
        q_col_indices <- grep("_q=", col_names)
        if (length(q_col_indices) > 0) {
            q_vals_str <- sub(".*_q=", "", col_names[q_col_indices])
            col_q_values <- as.numeric(q_vals_str)
            if (any(!is.na(col_q_values))) {
                names(col_q_values) <- q_col_indices
            } else {
                col_q_values <- NA
            }
        }
    }

    col_q_values
}

# ============================================================================
# HELPER: Apply post-hoc normalization to a single result SE
# ============================================================================
#' @noRd
.apply_diversity_post_hoc_norm <- function(result_se, norm_method, params, q_val,
    verbose) {
    if (is.null(norm_method) || norm_method == "default" || !is(result_se, "SummarizedExperiment")) {
        return(result_se)
    }

    diversity_assay <- SummarizedExperiment::assay(result_se, "diversity")

    # Apply the same normalization to CI assays if they exist
    ci_lower_assay <- NULL
    ci_upper_assay <- NULL

    if ("ci_lower" %in% SummarizedExperiment::assayNames(result_se)) {
        ci_lower_assay <- SummarizedExperiment::assay(result_se, "ci_lower")
    }
    if ("ci_upper" %in% SummarizedExperiment::assayNames(result_se)) {
        ci_upper_assay <- SummarizedExperiment::assay(result_se, "ci_upper")
    }

    if (norm_method == "zscore") {
        diversity_assay <- .normalize_zscore(diversity_assay, per_q = TRUE)
        if (!is.null(ci_lower_assay))
            ci_lower_assay <- .normalize_zscore(ci_lower_assay, per_q = TRUE)
        if (!is.null(ci_upper_assay))
            ci_upper_assay <- .normalize_zscore(ci_upper_assay, per_q = TRUE)
        if (verbose)
            message("[calculate_diversity_s4] Applied z-score normalization for q=",
                q_val)
    } else if (norm_method == "log_odds_ratio" && !is.null(params$genes)) {
        n_isoforms_vec <- table(params$genes)
        diversity_assay <- .normalize_log_odds_ratio(diversity_assay, n_isoforms = n_isoforms_vec,
            q = q_val)
        if (!is.null(ci_lower_assay))
            ci_lower_assay <- .normalize_log_odds_ratio(ci_lower_assay, n_isoforms = n_isoforms_vec,
                q = q_val)
        if (!is.null(ci_upper_assay))
            ci_upper_assay <- .normalize_log_odds_ratio(ci_upper_assay, n_isoforms = n_isoforms_vec,
                q = q_val)
        if (verbose)
            message("[calculate_diversity_s4] Applied log-odds ratio normalization for q=",
                q_val)
    } else if (norm_method == "relative_reference" && !is.null(params$reference_group)) {
        coldata <- SummarizedExperiment::colData(result_se)
        if (params$reference_group %in% colnames(coldata)) {
            group_vector <- coldata[[params$reference_group]]
            diversity_assay <- .normalize_relative_reference(diversity_assay, group_vector = group_vector,
                reference_group = params$reference_group)
            if (!is.null(ci_lower_assay))
                ci_lower_assay <- .normalize_relative_reference(ci_lower_assay, group_vector = group_vector,
                  reference_group = params$reference_group)
            if (!is.null(ci_upper_assay))
                ci_upper_assay <- .normalize_relative_reference(ci_upper_assay, group_vector = group_vector,
                  reference_group = params$reference_group)
            if (verbose)
                message("[calculate_diversity_s4] Applied relative reference normalization for q=",
                  q_val)
        }
    }

    SummarizedExperiment::assay(result_se, "diversity") <- diversity_assay
    if (!is.null(ci_lower_assay))
        SummarizedExperiment::assay(result_se, "ci_lower") <- ci_lower_assay
    if (!is.null(ci_upper_assay))
        SummarizedExperiment::assay(result_se, "ci_upper") <- ci_upper_assay

    result_se
}



