# ============================================================================
# DIVERSITY WRAPPER
# ============================================================================

#' Calculate diversity and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for Tsallis entropy.
#'   If NULL, uses q_values from \code{analysis@config$q_values} if available, else defaults to seq(0.01, 2, by = 0.05).
#' @param norm \code{logical} or \code{character}. Normalization method: TRUE, FALSE, "none", "range", "zscore", "log_odds_ratio", "relative_reference".
#'   If NULL, reads from \code{@config$norm} or defaults to TRUE.
#' @param norm_method \code{character}. Post-hoc normalization method applied after diversity computation.
#'   Options:
#'   \itemize{
#'     \item \code{"default"} - Simple normalization by theoretical maximum (current behavior)
#'     \item \code{"zscore"} - Z-score normalization per q-value
#'     \item \code{"log_odds_ratio"} - Log-odds ratio relative to max entropy (q and isoform-aware)
#'     \item \code{"relative_reference"} - Divide by reference group mean (requires reference_group)
#'     \item \code{NULL} - No post-hoc normalization (default)
#'   }
#'   If NULL, reads from \code{@config$norm_method} if available.
#' @param tpm \code{logical}. TPM normalization. Default: FALSE.
#'   If not specified, reads from \code{@config$tpm} if available.
#' @param assayno \code{numeric}. Assay number to use. Default: 1.
#'   If NULL, reads from \code{@config$assayno} if available.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#'   If not specified, reads from \code{@config$verbose} if available.
#' @param what \code{character}. Output type: "S" (entropy) or "D" (diversity). Default: "S".
#'   If NULL, reads from \code{@config$what} if available.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param pseudocount \code{numeric} or \code{character}. Pseudocount value or "auto". Default: 0.
#'   If NULL, reads from \code{@config$pseudocount} if available.
#' @param min_valid_frac \code{numeric}. Minimum valid fraction threshold. Default: 0.75.
#'   If NULL, reads from \code{@config$min_valid_frac} if available.
#' @param shrinkage \code{character}. Shrinkage method: "none" or "empirical_bayes". Default: "none".
#'   If NULL, reads from \code{@config$shrinkage} if available.
#' @param genes \code{character} or \code{NULL}. Gene set specification. Default: NULL (use all genes).
#'   If NULL, reads from \code{@config$genes} if available.
#' @param effective_length \code{numeric} or \code{NULL}. Effective gene lengths. Default: NULL.
#'   If NULL, reads from \code{@config$effective_length} if available.
#' @param metadata \code{list} or \code{NULL}. Additional metadata. Default: NULL.
#' @param bootstrap \code{logical}. Compute bootstrap confidence intervals. Default: FALSE.
#'   If not specified, reads from \code{@config$bootstrap} if available.
#' @param nboot \code{numeric} or \code{NULL}. Number of bootstrap replicates. Default: NULL.
#'   If NULL, reads from \code{@config$nboot} if available.
#' @param bootstrap_method \code{character}. Bootstrap method: "percentile" or others. Default: "percentile".
#'   If NULL, reads from \code{@config$bootstrap_method} if available.
#' @param bootstrap_ci \code{numeric}. Bootstrap confidence interval level (0-1). Default: 0.95.
#'   If NULL, reads from \code{@config$bootstrap_ci} if available.
#' @param bootstrap_include_diagnostics \code{logical}. Include bootstrap diagnostics. Default: TRUE.
#'   If not specified, reads from \code{@config$bootstrap_include_diagnostics} if available.
#' @param seed \code{numeric} or \code{NULL}. Random seed for reproducibility. Default: NULL.
#'   If NULL, reads from \code{@config$seed} if available. If provided, ensures reproducible bootstrap resampling.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Options:
#'   \itemize{
#'     \item \code{"default"} - Simple normalization by theoretical maximum (current behavior)
#'     \item \code{"zscore"} - Z-score normalization per q-value
#'     \item \code{"log_odds_ratio"} - Log-odds ratio relative to max entropy (q and isoform-aware)
#'     \item \code{"relative_reference"} - Divide by reference group mean (requires reference_group)
#'     \item \code{NULL} - No post-hoc normalization (default)
#'   }
#'   If NULL, reads from \code{@config$norm_method} if available.
#' @param reference_group \code{character}. For \code{norm_method = "relative_reference"}, 
#'   the reference group column name (e.g., from colData). If NULL, uses first group in colData.
#'   If NULL, reads from \code{@config$reference_group} if available.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   When provided, generates TWO files:
#'   \enumerate{
#'     \item Primary output: Analysis object (.rds) or table (.tsv/.csv/.txt)
#'     \item Secondary output: Diversity spectrum statistics (TSV format)
#'         with suffix \code{_diversity_spectrum.tsv}
#'   }
#'   Example: output_file = "analysis.rds" generates:
#'   \itemize{
#'     \item \code{analysis.rds} - TSENATAnalysis object
#'     \item \code{analysis_diversity_spectrum.tsv} - Spectrum statistics
#'   }
#'   The spectrum file contains columns: q, central (median diversity), 
#'   spread (IQR), count, and group (if grouping variable available).
#'   Default: NULL (no file output).
#' @param ... Additional arguments passed to underlying functions for extensibility.
#'
#' @return Modified TSENATAnalysis object with diversity results stored
#'   in \code{@diversity_results}, keyed by "q_X.XX..." format (e.g., "q_1.000").
#'   When \code{output_file} is provided, also generates:
#'   \itemize{
#'     \item Primary file: Analysis object or table export
#'     \item Spectrum file: \code{*_diversity_spectrum.tsv} containing 
#'           aggregated diversity statistics across q-values and groups
#'   }
#'
#' @details
#' This wrapper calls \code{.calculate_diversity()} once per q-value, storing
#' results as SummarizedExperiment objects. It extracts key parameters from
#' \code{analysis@config} with priority resolution (explicit > \code{@config} > default).
#'
#' **Diversity Spectrum Computation:**
#' By default, this function computes and saves a diversity spectrum (aggregated
#' statistics across all q-values and groups) when \code{output_file} is provided.
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
#' Spectrum is saved as: \code{*_diversity_spectrum.tsv}
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{q}{Priority 1 (explicit) > Priority 2 (\code{@config$q_values}) > Priority 3 (default: seq(0.01, 2, by = 0.05))\cr
#'     **Note:** If explicit q AND \code{@config$q_values} both provided, explicit wins.}
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
#'   \item{verbose}{Priority: explicit > \code{@config$verbose} > TRUE}
#'   \item{bootstrap}{Priority: explicit > \code{@config$bootstrap} > FALSE}
#'   \item{pseudocount}{Priority: explicit > \code{@config$pseudocount} > 0}
#'   \item{norm}{Priority: explicit > \code{@config$norm} > TRUE}
#'   \item{what}{Priority: explicit > \code{@config$what} > "S" (Tsallis entropy)}
#' }
#'
#' **Audit Trail:** After execution, check:
#' \itemize{
#'   \item \code{analysis@config$last_diversity_run$parameters_used}: Actual parameters 
#'         (not original \code{@config})
#'   \item \code{attr(diversity(analysis, q=X), "computed_with")}: Per-q-value metadata
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
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' 
#' analysis <- build_analysis_s4(readcounts, gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' 
#' # Compute diversity and access results
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0), verbose = FALSE)
#' head(diversity(analysis, q = 1.0))
#' 
#' # Apply z-score normalization
#' # analysis <- calculate_diversity_s4(analysis, norm_method = "zscore")
#' 
#' # Apply log-odds ratio normalization (q and isoform-aware)
#' # analysis <- calculate_diversity_s4(analysis, norm_method = "log_odds_ratio")
#'
#' @details
#' For additional details on diversity spectrum calculations and normalization methods,
#' see the package vignettes.
#'
#' @export
#' @importFrom utils write.table
calculate_diversity_s4 <- function(analysis, q = NULL, norm = NULL, norm_method = NULL, 
                                   reference_group = NULL, tpm = FALSE, assayno = NULL,
                                   verbose = NULL, what = NULL, nthreads = NULL, pseudocount = NULL,
                                   min_valid_frac = NULL, shrinkage = NULL, genes = NULL,
                                   effective_length = NULL, metadata = NULL, bootstrap = NULL,
                                   nboot = NULL, bootstrap_method = NULL, bootstrap_ci = NULL,
                                   bootstrap_include_diagnostics = NULL, seed = NULL, output_file = NULL, ...) {
  # Validate input
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  if (nrow(analysis@se) == 0) {
    stop("SummarizedExperiment in @se is empty", call. = FALSE)
  }
  
  # Debug: check initial state
  if (verbose %||% FALSE) {
    message("[DEBUG] Initial analysis@diversity_results: ", 
            paste(names(analysis@diversity_results), collapse=", "))
  }

  # Prepare all parameters (resolve from explicit args > @config > defaults)
  params <- .prepare_diversity_params(
    analysis, q, norm, norm_method, reference_group, tpm, assayno, verbose, what, 
    nthreads, pseudocount, min_valid_frac, shrinkage, genes, effective_length, 
    metadata, bootstrap, nboot, bootstrap_method, bootstrap_ci, 
    bootstrap_include_diagnostics, seed
  )
  
  # Validate norm_method parameter
  if (!is.null(params$norm_method)) {
    valid_methods <- c("default", "zscore", "log_odds_ratio", "relative_reference")
    if (!(params$norm_method %in% valid_methods)) {
      stop("'norm_method' must be one of: ", paste(valid_methods, collapse=", "), 
           call. = FALSE)
    }
  }
  
  # Build calculation arguments for .calculate_diversity()
  calc_args <- .build_calc_diversity_args(params, analysis, list(...))
  
  if (params$verbose) {
    message("[DEBUG] calc_args names: ", paste(names(calc_args), collapse=", "))
    if (!"genes" %in% names(calc_args)) {
        message("[DEBUG] NOTE: genes is NOT in calc_args - will default to SE rownames")
    } else {
        message("[DEBUG] genes provided in calc_args: ", paste(head(calc_args$genes, 5), collapse=", "))
    }
    x_param <- calc_args$x
    if (is(x_param, "SummarizedExperiment")) {
      md_names <- names(S4Vectors::metadata(x_param))
      message("[DEBUG] Metadata on x: ", paste(md_names, collapse=", "))
      if ("analysis" %in% md_names) {
        message("[DEBUG] WARNING: 'analysis' in x metadata!")
      }
    }
  }
  
  result_df <- do.call(.calculate_diversity, calc_args)
  
  if (params$verbose) {
    message("[DEBUG] After .calculate_diversity(): analysis@diversity_results keys = ", 
            paste(names(analysis@diversity_results), collapse=", "))
    message("[DEBUG] params$q STILL: ", paste(params$q, collapse=", "))
  }
  
  # Store the original result SE for later per-q extraction
  result_se_original <- result_df
  
  # Extract q-value metadata from result
  col_q_values <- .extract_q_metadata_from_result(result_df, params$q)
  
  if (params$verbose) {
    message("[DEBUG] result_df colnames: ", paste(colnames(result_df), collapse=", "))
    message("[DEBUG] col_q_values extracted: ", paste(col_q_values, collapse=", "))
    if (length(names(col_q_values)) > 0) {
      message("[DEBUG] col_q_values names (indices): ", paste(names(col_q_values), collapse=", "))
    }
  }
  
  # Store combined result in cache
  analysis@metadata$diversity_combined <- list(
    combined_result = result_df,
    combined_se = result_df,
    q_values_computed = params$q,
    computation_params = list(
      norm = params$norm,
      verbose = params$verbose,
      bootstrap = params$bootstrap,
      pseudocount = params$pseudocount,
      nthreads = params$nthreads,
      what = params$what
    ),
    timestamp = Sys.time()
  )
  analysis@metadata$.diversity_lazily_computed <- TRUE
  
  # Store per-q results with post-hoc normalization and audit trail
  q_decimals <- 3
  if (params$verbose) {
    message("[DEBUG] params$q values: ", paste(params$q, collapse=", "))
    message("[DEBUG] Loop will iterate over params$q, expecting ", length(params$q), " iterations")
    message("[DEBUG] diversity_results BEFORE loop: ", paste(names(analysis@diversity_results), collapse=", "))
  }
  
  for (q_val in params$q) {
    tryCatch({
      if (params$verbose) {
        message("[DEBUG] Processing q=", q_val, " (formatted: q_", 
                formatC(q_val, format = "f", digits = q_decimals), ")")
      }
      # Handle empty results
      if (nrow(result_df) == 0) {
        warning("[calculate_diversity_s4] Result for q=", q_val, " is empty (0 rows)", call. = FALSE)
        next
      }
      
      # Extract columns for this q-value
      q_cols <- NULL
      if (!is.na(col_q_values[1])) {
        q_mask <- abs(col_q_values - q_val) < 1e-5
        if (any(q_mask)) {
          q_cols <- as.numeric(names(col_q_values)[q_mask])
        }
        if (params$verbose) {
          message("[DEBUG] q=", q_val, ": col_q_values matching approach found cols: ", 
                  paste(q_cols, collapse=", "))
        }
      }
      
      # Fallback: pattern matching (handle different decimal precision)
      if (is.null(q_cols) || length(q_cols) == 0) {
        for (ndigits in c(3, 2)) {
          q_formatted <- formatC(q_val, format = "f", digits = ndigits)
          q_pattern_str <- paste0("_q=", gsub("\\.", "\\\\.", q_formatted), "$")
          q_cols_temp <- grep(q_pattern_str, colnames(result_df))
          if (params$verbose) {
            message("[DEBUG] q=", q_val, ": pattern (", ndigits, " digits) '", q_pattern_str, 
                    "' found cols: ", paste(q_cols_temp, collapse=", "))
          }
          if (length(q_cols_temp) > 0) {
            q_cols <- q_cols_temp
            break
          }
        }
      }
      
      # Fallback 3: single q-value case
      if (length(q_cols) == 0 && length(params$q) == 1) {
        numeric_cols <- vapply(result_df, is.numeric, FUN.VALUE = logical(1))
        q_cols <- which(numeric_cols)
      }
      
      # Validate columns found
      if (length(q_cols) == 0) {
        stop("[calculate_diversity_s4] No columns found for q=", q_val,
             ". Available columns: ", paste(head(colnames(result_df), 10), collapse=", "),
             call. = FALSE)
      }
      
      # Extract and convert to SE, preserving CI assays if present
      result_subset <- result_df[, q_cols, drop = FALSE]
      
      if (is(result_subset, "SummarizedExperiment")) {
        # result_subset is already a SummarizedExperiment with CI assays preserved
        result_se <- result_subset
      } else if (is.data.frame(result_subset)) {
        numeric_cols <- vapply(result_subset, is.numeric, FUN.VALUE = logical(1))
        if (!any(numeric_cols)) {
          result_se <- result_subset
        } else {
          # Extract main diversity assay
          assay_data <- as.matrix(result_subset[, numeric_cols, drop = FALSE])
          
          # Build assays list with diversity and CIs
          assays_list <- list(diversity = assay_data)
          
          # Extract CI assays from the original result SE if available
          # The q_cols indices also apply to the CI assays since they have the same structure
          if (is(result_se_original, "SummarizedExperiment")) {
            if ("ci_lower" %in% SummarizedExperiment::assayNames(result_se_original)) {
              ci_lower_orig <- SummarizedExperiment::assay(result_se_original, "ci_lower")
              if (!is.null(ci_lower_orig)) {
                # Subset the same columns from ci_lower
                ci_lower_subset <- ci_lower_orig[, q_cols, drop = FALSE]
                assays_list$ci_lower <- ci_lower_subset
              }
            }
            if ("ci_upper" %in% SummarizedExperiment::assayNames(result_se_original)) {
              ci_upper_orig <- SummarizedExperiment::assay(result_se_original, "ci_upper")
              if (!is.null(ci_upper_orig)) {
                # Subset the same columns from ci_upper
                ci_upper_subset <- ci_upper_orig[, q_cols, drop = FALSE]
                assays_list$ci_upper <- ci_upper_subset
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
        }
      } else {
        result_se <- result_subset
      }
      
      # Validate SE structure
      if (is(result_se, "SummarizedExperiment")) {
        if (length(SummarizedExperiment::assays(result_se)) == 0) {
          stop("[calculate_diversity_s4] Converted SE for q=", q_val, 
               " has no assays. Check diversity result structure.", call. = FALSE)
        }
        test_assay <- tryCatch({
          SummarizedExperiment::assay(result_se, 1)
        }, error = function(e) {
          stop("[calculate_diversity_s4] Cannot access assay for q=", q_val, ": ", 
               conditionMessage(e), call. = FALSE)
        })
        if (is.null(test_assay) || nrow(test_assay) == 0) {
          warning("[calculate_diversity_s4] Assay for q=", q_val, 
                  " is empty. This may cause issues downstream.", call. = FALSE)
        }
      }
      
      # Apply original SE's colData to preserve metadata
      if (is(result_se, "SummarizedExperiment") && ncol(result_se) > 0) {
        original_coldata <- SummarizedExperiment::colData(analysis@se)
        if (!is.null(original_coldata) && nrow(original_coldata) == ncol(result_se)) {
          SummarizedExperiment::colData(result_se) <- original_coldata
        }
      }
      
      # Apply post-hoc normalization
      result_se <- .apply_diversity_post_hoc_norm(result_se, params$norm_method, params, q_val, params$verbose)
      
      # Store with audit trail metadata
      key <- paste0("q_", formatC(q_val, format = "f", digits = q_decimals))
      
      # Debug: verify assays before storing
      if (params$verbose) {
        assay_names <- SummarizedExperiment::assayNames(result_se)
        message("[calculate_diversity_s4] Storing q=", q_val, " with assays: ", 
                paste(assay_names, collapse=", "))
      }
      
      attr(result_se, "computed_with") <- list(
        q = q_val, norm = params$norm, norm_method = params$norm_method,
        verbose = params$verbose, bootstrap = params$bootstrap, pseudocount = params$pseudocount,
        nthreads = params$nthreads, what = params$what, timestamp = Sys.time()
      )
      analysis@diversity_results[[key]] <- result_se
      
      if (params$verbose) {
        message("[DEBUG] After storing ", key, ": diversity_results keys = ", 
                paste(names(analysis@diversity_results), collapse=", "))
      }
      
      # Track function call
      analysis@metadata$function_calls <- c(
        analysis@metadata$function_calls,
        paste0("calculate_diversity[q=", q_val, "]")
      )
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      if (params$bootstrap && grepl("bootstrap", error_msg, ignore.case = TRUE)) {
        stop("[calculate_diversity_s4] Bootstrap CI computation failed for q=", q_val, ": ", 
             error_msg, call. = FALSE)
      } else {
        stop("[calculate_diversity_s4] Failed to compute diversity for q=", q_val, ":\n", 
             error_msg, call. = FALSE)
      }
    })
  }
  
  # Update config with actual parameters used (audit trail)
  analysis@config$last_diversity_run <- list(
    timestamp = Sys.time(),
    q_values_computed = params$q,
    num_q_values = length(params$q),
    parameters_used = list(
      norm = params$norm, norm_method = params$norm_method, verbose = params$verbose,
      bootstrap = params$bootstrap, pseudocount = params$pseudocount, nthreads = params$nthreads,
      what = params$what
    ),
    note = "Actual parameters used (save object and check this, not original @config)"
  )
  
  # Track parallel processing
  if (params$nthreads > 1) {
    analysis@metadata$parallel_processing <- c(
      analysis@metadata$parallel_processing,
      paste0("calculate_diversity_s4: nthreads=", params$nthreads, " (", length(params$q), " q-values)")
    )
  }
  
  # Compute and save diversity spectrum
  if (!is.null(output_file)) {
    tryCatch({
      spectrum_condition_col <- analysis@config$condition_col %||% "sample_type"
      combined_se <- analysis@metadata$diversity_combined$combined_se
      if (!is.null(combined_se) && nrow(combined_se) > 0) {
        diversity_spectrum <- .compute_diversity_spectrum(
          se = combined_se, metric = "median", variability_metric = "iqr", 
          condition_col = spectrum_condition_col
        )
        
        if (!is.null(diversity_spectrum) && nrow(diversity_spectrum) > 0) {
          spectrum_file <- sub("\\.[^.]+$", "_diversity_spectrum.tsv", output_file)
          if (spectrum_file == output_file) {
            spectrum_file <- paste0(output_file, "_diversity_spectrum.tsv")
          }
          
          utils::write.table(
            diversity_spectrum, file = spectrum_file, sep = "\t", row.names = FALSE, quote = FALSE
          )
          
          if (params$verbose) {
            message("[calculate_diversity_s4] Saved diversity spectrum to: ", spectrum_file)
          }
          analysis@metadata$diversity_spectrum <- diversity_spectrum
        }
      }
    }, error = function(e) {
      warning("[calculate_diversity_s4] Could not compute diversity spectrum: ", 
              conditionMessage(e), call. = FALSE)
    })
  }
  
  # Save output if requested
  if (!is.null(output_file)) {
    tryCatch({
      output_data <- NULL
      
      # Extract diversity results from individual SE objects (preferred path with CIs)
      if (length(analysis@diversity_results) > 0) {
        # Build output data systematically - iterate over the q values we REQUESTED, not all stored results
        q_keys_to_use <- paste0("q_", formatC(params$q, format = "f", digits = 3))
        n_q <- length(q_keys_to_use)
        
        # Get dimensions from first result
        first_key <- q_keys_to_use[1]
        if (first_key %in% names(analysis@diversity_results)) {
          se_first_actual <- analysis@diversity_results[[first_key]]
          all_genes <- rownames(se_first_actual)
          all_samples <- colnames(se_first_actual)
        } else {
          # Fallback to first available
          all_genes <- rownames(analysis@diversity_results[[1]])
          all_samples <- colnames(analysis@diversity_results[[1]])
        }
        
        n_genes <- length(all_genes)
        n_samples <- length(all_samples)
        
        # Pre-allocate data frame
        total_rows <- n_genes * n_samples * n_q
        output_data <- data.frame(
          gene = character(total_rows),
          sample = character(total_rows),
          q_value = character(total_rows),
          diversity = numeric(total_rows),
          stringsAsFactors = FALSE
        )
        
        # Add CI columns if they exist
        se_first <- analysis@diversity_results[[first_key]]
        assay_names_first <- SummarizedExperiment::assayNames(se_first)
        has_ci_lower <- "ci_lower" %in% assay_names_first
        has_ci_upper <- "ci_upper" %in% assay_names_first
        
        if (has_ci_lower) output_data$ci_lower <- numeric(total_rows)
        if (has_ci_upper) output_data$ci_upper <- numeric(total_rows)
        
        # Fill data frame - use ONLY the requested q-values
        row_idx <- 1
        if (params$verbose) {
          message("[DEBUG] Using q-keys for output: ", paste(q_keys_to_use, collapse=", "))
        }
        
        for (q_idx in seq_along(q_keys_to_use)) {
          q_key <- q_keys_to_use[q_idx]
          
          if (!q_key %in% names(analysis@diversity_results)) {
            if (params$verbose) {
              message("[DEBUG] WARNING: q_key '", q_key, "' not found in diversity_results")
            }
            next
          }
          
          se <- analysis@diversity_results[[q_key]]
          q_name <- q_key
          
          if (!is(se, "SummarizedExperiment")) {
            next
          }
          
          diversity_mat <- as.matrix(SummarizedExperiment::assay(se, 1))
          
          # Debug: show assays for this SE
          assay_names <- SummarizedExperiment::assayNames(se)
          if (params$verbose && q_idx == 1) {
            message("[calculate_diversity_s4] First q SE has assays: ", paste(assay_names, collapse=", "))
          }
          
          ci_lower_mat <- if (has_ci_lower) {
            tryCatch({
              ci_lower_temp <- as.matrix(SummarizedExperiment::assay(se, "ci_lower"))
              if (params$verbose && q_idx == 1) {
                message("[DEBUG] ci_lower[1,1:2]: ", paste(ci_lower_temp[1, 1:min(2, ncol(ci_lower_temp))], collapse=", "))
              }
              ci_lower_temp
            }, error = function(e) {
              message("[calculate_diversity_s4] ERROR accessing ci_lower for q=", q_name, ": ", conditionMessage(e))
              message("  Available assays: ", paste(assay_names, collapse=", "))
              NULL
            })
          } else NULL
          
          if (params$verbose && q_idx == 1) {
            message("[DEBUG] diversity[1,1:2]: ", paste(diversity_mat[1, 1:min(2, ncol(diversity_mat))], collapse=", "))
            message("[DEBUG] ci_upper row 1 first 2 cols will be shown next")
          }
          
          ci_upper_mat <- if (has_ci_upper) {
            tryCatch({
              ci_upper_temp <- as.matrix(SummarizedExperiment::assay(se, "ci_upper"))
              if (params$verbose && q_idx == 1) {
                message("[DEBUG] ci_upper[1,1:2]: ", paste(ci_upper_temp[1, 1:min(2, ncol(ci_upper_temp))], collapse=", "))
              }
              ci_upper_temp
            }, error = function(e) {
              message("[calculate_diversity_s4] ERROR accessing ci_upper for q=", q_name, ": ", conditionMessage(e))
              message("  Available assays: ", paste(assay_names, collapse=", "))
              NULL
            })
          } else NULL
          
          for (gene_idx in seq_len(nrow(diversity_mat))) {
            for (sample_idx in seq_len(ncol(diversity_mat))) {
              output_data$gene[row_idx] <- rownames(diversity_mat)[gene_idx]
              output_data$sample[row_idx] <- colnames(diversity_mat)[sample_idx]
              output_data$q_value[row_idx] <- q_name
              output_data$diversity[row_idx] <- as.numeric(diversity_mat[gene_idx, sample_idx])
              
              if (has_ci_lower && !is.null(ci_lower_mat)) {
                output_data$ci_lower[row_idx] <- as.numeric(ci_lower_mat[gene_idx, sample_idx])
              }
              if (has_ci_upper && !is.null(ci_upper_mat)) {
                output_data$ci_upper[row_idx] <- as.numeric(ci_upper_mat[gene_idx, sample_idx])
              }
              
              row_idx <- row_idx + 1
            }
          }
        }
        
        # Trim to actual rows (in case rows < total_rows due to errors)
        output_data <- output_data[1:(row_idx - 1), ]
      }
      
      # Fallback: if no diversity_results, try combined_result
      if (is.null(output_data) || nrow(output_data) == 0) {
        combined <- analysis@metadata$diversity_combined$combined_result
        if (!is.null(combined) && nrow(combined) > 0) {
          if (is(combined, "SummarizedExperiment")) {
            output_data <- as.data.frame(SummarizedExperiment::assay(combined, 1))
          } else if (is.data.frame(combined)) {
            output_data <- combined
          } else if (is.matrix(combined)) {
            output_data <- as.data.frame(combined)
          }
        }
      }
      
      # Save if we have data
      if (!is.null(output_data) && nrow(output_data) > 0) {
        save_analysis_output(output_data, output_file, verbose = params$verbose,
                             func_name = "calculate_diversity_s4")
        
        if (params$verbose) {
          message("[calculate_diversity_s4] Saved diversity results to: ", output_file)
        }
      }
    }, error = function(e) {
      warning("[calculate_diversity_s4] Could not save diversity results: ",
              conditionMessage(e), call. = FALSE)
    })
  }

  
  analysis
}

# ============================================================================
# HELPER: Prepare and resolve all parameters for diversity calculation
# ============================================================================
#' @noRd
.prepare_diversity_params <- function(analysis, q, norm, norm_method, reference_group,
                                      tpm, assayno, verbose, what, nthreads, pseudocount,
                                      min_valid_frac, shrinkage, genes, effective_length,
                                      metadata, bootstrap, nboot, bootstrap_method,
                                      bootstrap_ci, bootstrap_include_diagnostics, seed) {
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
  nthreads_resolved <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  
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
  
  list(
    q = q,
    nthreads = nthreads_resolved,
    verbose = resolve_slot_param(verbose, analysis@config, "verbose", TRUE),
    bootstrap = resolve_slot_param(bootstrap, analysis@config, "bootstrap", FALSE),
    pseudocount = resolve_slot_param(pseudocount, analysis@config, "pseudocount", 0),
    norm = resolve_slot_param(norm, analysis@config, "norm", TRUE),
    what = resolve_slot_param(what, analysis@config, "what", "S"),
    assayno = resolve_slot_param(assayno, analysis@config, "assayno", 1),
    min_valid_frac = resolve_slot_param(min_valid_frac, analysis@config, "min_valid_frac", 0.75),
    shrinkage = resolve_slot_param(shrinkage, analysis@config, "shrinkage", "none"),
    bootstrap_method = resolve_slot_param(bootstrap_method, analysis@config, "bootstrap_method", "percentile"),
    bootstrap_ci = resolve_slot_param(bootstrap_ci, analysis@config, "bootstrap_ci", 0.95),
    seed = resolve_slot_param(seed, analysis@config, "seed", NULL),
    tpm = tpm,
    genes = resolve_slot_param(genes, analysis@config, "genes", NULL),
    effective_length = resolve_slot_param(effective_length, analysis@config, "effective_length", NULL),
    nboot = resolve_slot_param(nboot, analysis@config, "nboot", NULL),
    bootstrap_include_diagnostics = resolve_slot_param(bootstrap_include_diagnostics, analysis@config, "bootstrap_include_diagnostics", TRUE),
    metadata = resolve_slot_param(metadata, analysis@config, "metadata", NULL),
    norm_method = resolve_slot_param(norm_method, analysis@config, "norm_method", NULL),
    reference_group = resolve_slot_param(reference_group, analysis@config, "reference_group", NULL)
  )
}

# ============================================================================
# HELPER: Build calculation arguments for .calculate_diversity()
# ============================================================================
#' @noRd
.build_calc_diversity_args <- function(params, analysis, dots) {
  calc_args <- list(
    x = analysis@se,
    q = params$q,
    norm = params$norm,
    tpm = params$tpm,
    assayno = params$assayno,
    verbose = params$verbose,
    what = params$what,
    nthreads = params$nthreads,
    pseudocount = params$pseudocount,
    min_valid_frac = params$min_valid_frac,
    shrinkage = params$shrinkage,
    bootstrap = params$bootstrap,
    bootstrap_nboot = params$nboot,
    bootstrap_method = params$bootstrap_method,
    bootstrap_ci = params$bootstrap_ci,
    bootstrap_include_diagnostics = params$bootstrap_include_diagnostics,
    seed = params$seed
  )
  
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
.apply_diversity_post_hoc_norm <- function(result_se, norm_method, params, q_val, verbose) {
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
    if (!is.null(ci_lower_assay)) ci_lower_assay <- .normalize_zscore(ci_lower_assay, per_q = TRUE)
    if (!is.null(ci_upper_assay)) ci_upper_assay <- .normalize_zscore(ci_upper_assay, per_q = TRUE)
    if (verbose) message("[calculate_diversity_s4] Applied z-score normalization for q=", q_val)
  } else if (norm_method == "log_odds_ratio" && !is.null(params$genes)) {
    n_isoforms_vec <- table(params$genes)
    diversity_assay <- .normalize_log_odds_ratio(diversity_assay, n_isoforms = n_isoforms_vec, q = q_val)
    if (!is.null(ci_lower_assay)) ci_lower_assay <- .normalize_log_odds_ratio(ci_lower_assay, n_isoforms = n_isoforms_vec, q = q_val)
    if (!is.null(ci_upper_assay)) ci_upper_assay <- .normalize_log_odds_ratio(ci_upper_assay, n_isoforms = n_isoforms_vec, q = q_val)
    if (verbose) message("[calculate_diversity_s4] Applied log-odds ratio normalization for q=", q_val)
  } else if (norm_method == "relative_reference" && !is.null(params$reference_group)) {
    coldata <- SummarizedExperiment::colData(result_se)
    if (params$reference_group %in% colnames(coldata)) {
      group_vector <- coldata[[params$reference_group]]
      diversity_assay <- .normalize_relative_reference(diversity_assay, group_vector = group_vector,
                                                       reference_group = params$reference_group)
      if (!is.null(ci_lower_assay)) ci_lower_assay <- .normalize_relative_reference(ci_lower_assay, group_vector = group_vector, reference_group = params$reference_group)
      if (!is.null(ci_upper_assay)) ci_upper_assay <- .normalize_relative_reference(ci_upper_assay, group_vector = group_vector, reference_group = params$reference_group)
      if (verbose) message("[calculate_diversity_s4] Applied relative reference normalization for q=", q_val)
    }
  }
  
  SummarizedExperiment::assay(result_se, "diversity") <- diversity_assay
  if (!is.null(ci_lower_assay)) SummarizedExperiment::assay(result_se, "ci_lower") <- ci_lower_assay
  if (!is.null(ci_upper_assay)) SummarizedExperiment::assay(result_se, "ci_upper") <- ci_upper_assay
  
  result_se
}



