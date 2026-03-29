# S4 Wrapper Functions for TSENAT Pipeline
# These functions provide S4-integrated alternatives to existing analysis
# functions. They extract input from TSENATAnalysis slots, run analysis,
# and store results back to appropriate slots.
#



#' Calculate LM interactions and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param fdr_threshold \code{numeric}. FDR cutoff for significance.
#'   Default: 0.05.
#' @param formula \code{formula} or NULL. Reserved for future use.
#' @param condition_col \code{character} or \code{NULL}. Column name in colData identifying sample conditions.
#'   If NULL, reads from \code{@config$condition_col} or auto-detects common column names.
#' @param method \code{character}. Statistical method (e.g., "lmm", "gam", "gee").
#'   If NULL, uses method from @config$method or defaults to "lmm".
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param subject_col \code{character} or \code{NULL}. Column name identifying subject IDs for paired designs.
#'   If NULL, reads from \code{@config$subject_col} if available.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to NULL, letting base function decide).
#' @param multicorr \code{character} or \code{NULL}. Multiple comparison correction method.
#'   Options: "hochberg", "westfall-young", "benjamini-yekutieli".
#'   If NULL, uses method from @config or base function defaults.
#' @param corstr \code{character} or \code{NULL}. Correlation structure for GEE models.
#'   Options: "ar1", "exchangeable", "independence".
#'   If NULL, uses method from @config or base function defaults.
#' @param pcorr \code{character} or \code{NULL}. P-value correction method.
#'   Default: "BH" (Benjamini-Hochberg).
#'   If NULL, reads from \code{@config$pcorr} if available.
#' @param verbose \code{logical}. Print progress messages. Default: FALSE.
#' @param return_model_data \code{logical}. Return model data for visualization. Default: TRUE.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base LM function,
#'   including: pvalue, min_obs, assay_name, bias_correction, regularization, storey, wy_randomizations, adaptive_knots, etc.
#'
#' @return Modified TSENATAnalysis with results in @lm_results$lm_interaction.
#'
#' @details
#' Extracts diversity results from @diversity_results (prerequisite),
#' combines across q-values into single SummarizedExperiment,
#' then runs \code{.calculate_lm_interaction()}.
#' 
#' **Parameter Priority Resolution:**
#' \itemize{
#'   \item \code{nthreads}: Priority: explicit > @config > NULL
#' }
#' 
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' @examples
#' # Create test data with sufficient structure for LM analysis
#' # Create test analysis with diversity pre-computed
#' analysis <- TSENAT:::.create_test_analysis(
#'   n_genes = 8,
#'   n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5)
#' )
#' 
#' # Calculate q × condition interactions (LM)
#' analysis <- calculate_lm_interaction_s4(analysis, verbose = FALSE)
#' 
#' # View results
#' head(lmResults(analysis, "lm_interaction"))
#'
#' @export
#' @importFrom utils write.table
calculate_lm_interaction_s4 <- function(analysis, fdr_threshold = NULL, 
                                       formula = NULL, condition_col = NULL, method = NULL,
                                       paired = FALSE, subject_col = NULL, nthreads = NULL,
                                       multicorr = NULL, corstr = NULL, pcorr = NULL,
                                       verbose = FALSE, return_model_data = TRUE,
                                       output_file = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # OPTIMIZATION: Clear memoization cache for new analysis
  # Ensures fresh calculations for new dataset
  .clear_lm_helper_cache()

  # Check prerequisites
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # Ensure colData consistency: If diversity results have colData, sync to analysis@se
  # This maintains pairing structure through the analysis pipeline
  if (length(analysis@diversity_results) > 0) {
    # Get colData from first diversity result (should be consistent across q-values)
    first_diversity_se <- analysis@diversity_results[[1]]
    if (is(first_diversity_se, "SummarizedExperiment") && ncol(first_diversity_se) > 0) {
      diversity_coldata <- SummarizedExperiment::colData(first_diversity_se)
      if (!is.null(diversity_coldata) && nrow(diversity_coldata) == ncol(analysis@se)) {
        # If diversity results have richer colData, use that
        if (!is.null(diversity_coldata) && ncol(diversity_coldata) > 0) {
          SummarizedExperiment::colData(analysis@se) <- diversity_coldata
        }
      }
    }
  }

  # =========================================================================
  # PARAMETER EXTRACTION FROM @config using centralized helpers
  # =========================================================================
  
  # Auto-detect condition_col if not provided
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(analysis@se))
    condition_col <- auto_detect_column(
      cd_cols,
      config_list = analysis@config,
      config_key = "condition_col",
      priority_candidates = c("condition", "sample_type", "group", "treatment"),
      default_fallback = NULL,
      verbose = verbose,
      param_name = "condition_col"
    )
  }
  
  # Resolve remaining parameters using centralized handler
  method <- resolve_slot_param(method, analysis@config, "method", "lmm")
  subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col", NULL)
  multicorr <- resolve_slot_param(multicorr, analysis@config, "multicorr", NULL)
  corstr <- resolve_slot_param(corstr, analysis@config, "corstr", NULL)
  pcorr <- resolve_slot_param(pcorr, analysis@config, "pcorr", "BH")
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", NULL)
  
  # Paired design (logical flag)
  if (!paired && "paired" %in% names(analysis@config)) {
    paired <- analysis@config$paired
  }
  
  # Validate that required parameters are available
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(analysis@se))
    if (length(cd_cols) > 0) {
      message("condition_col not specified. Available columns: ",
          paste(cd_cols, collapse = ", "))
    } else {
      message("condition_col not specified and colData is empty. ",
          "Will be determined by .calculate_lm_interaction().")
    }
  }

  # Extract q-values from diversity_results keys (format: "q_0.5", "q_1.0", etc.)
  q_keys <- names(analysis@diversity_results)
  
  # =========================================================================
  # APPROACH: Combine existing per-q diversity results
  # NO recomputation - just stack the already-computed SEs horizontally
  # Add _q= suffix to column names to distinguish q-values
  # =========================================================================
  
  # Combine per-q SEs into single multi-q SE for calculate_lm_interaction
  diversity_combined <- tryCatch({
    # Extract assay matrices and colData from each per-q SE
    assay_list <- list()
    coldata_list <- list()
    rowdata_first <- NULL
    
    for (key in sort(q_keys)) {
      se <- analysis@diversity_results[[key]]
      if (!is(se, "SummarizedExperiment")) {
        stop("[calculate_lm_interaction_s4] Diversity result ", key, 
             " is not a SummarizedExperiment", call. = FALSE)
      }
      
      # Extract q-value from key (e.g., "q_0.100" -> 0.1)
      q_val_str <- sub("^q_", "", key)
      q_val <- as.numeric(q_val_str)
      
      # Get assay matrix (should be named "diversity")
      if (length(SummarizedExperiment::assays(se)) == 0) {
        stop("[calculate_lm_interaction_s4] Diversity SE ", key, 
             " has no assays", call. = FALSE)
      }
      assay_matrix <- as.matrix(SummarizedExperiment::assay(se))
      
      # Add _q= suffix to ASSAY column names ONLY to distinguish q-values
      # DO NOT add suffix to colData column names (rbind requires identical column names across all DFrames)
      assay_colnames_with_q <- paste0(colnames(assay_matrix), "_q=", q_val_str)
      colnames(assay_matrix) <- assay_colnames_with_q
      
      assay_list[[key]] <- assay_matrix
      
      # Get colData - IMPORTANT: colData column names must remain the same across all q-values
      # Only update rownames to match new assay column names
      coldata <- SummarizedExperiment::colData(se)
      rownames(coldata) <- assay_colnames_with_q
      coldata_list[[key]] <- coldata
      
      # Save rowData from first SE (same genes for all q-values)
      if (is.null(rowdata_first)) {
        rowdata_first <- SummarizedExperiment::rowData(se)
      }
    }
    
    # Combine assays horizontally (columns from different q-values)
    combined_assay <- do.call(cbind, assay_list)
    
    # Combine colData vertically (each q-value's samples with q-suffix in rownames)
    combined_coldata <- do.call(rbind, coldata_list)
    
    # Ensure colnames of assay match rownames of colData
    colnames(combined_assay) <- rownames(combined_coldata)
    
    # Create combined SE
    se_combined <- SummarizedExperiment::SummarizedExperiment(
      assays = list(diversity = combined_assay),
      colData = combined_coldata,
      rowData = rowdata_first
    )
    
    se_combined
  }, error = function(e) {
    stop("[calculate_lm_interaction_s4] Failed to combine diversity results: ",
                conditionMessage(e), call. = FALSE)
  })
  
  # Build arguments for calculate_lm_interaction
  args <- list(se = diversity_combined)
  
  if (!is.null(method)) {
    args$method <- method
  }
  
  # Add extracted parameters
  if (!is.null(condition_col)) {
    args$condition_col <- condition_col
  }
  
  if (paired) {
    args$paired <- paired
  }
  
  if (!is.null(subject_col)) {
    args$subject_col <- subject_col
  }

  if (!is.null(nthreads)) {
    args$nthreads <- nthreads
  }

  if (!is.null(multicorr)) {
    args$multicorr <- multicorr
  }

  if (!is.null(corstr)) {
    args$corstr <- corstr
  }

  if (!is.null(pcorr)) {
    args$pcorr <- pcorr
  }

  if (verbose) {
    args$verbose <- verbose
  }

  # Request model_data for plotting compatibility
  args$return_model_data <- return_model_data
  
  # Merge with additional args (which may override these values)
  args <- c(args, list(...))

  # Run LM analysis
  result <- tryCatch({
    do.call(.calculate_lm_interaction, args)
  }, error = function(e) {
    stop("lm_interaction calculation failed:\n", conditionMessage(e),
         call. = FALSE)
  })

  # =========================================================================
  # HANDLE RESULT: May be list($results, $model_data) or just data.frame
  # =========================================================================
  lm_results_df <- result
  model_data <- NULL
  
  # Check if result is a list with $results and $model_data (return_model_data=TRUE)
  if (is.list(result) && "results" %in% names(result)) {
    lm_results_df <- result$results
    model_data <- result$model_data
  }
  
  # =========================================================================
  # CRITICAL VALIDATION: Ensure result has required structure
  # =========================================================================
  # The result MUST be a data.frame with columns 'gene' and 'adj_p_interaction'
  # to be compatible with downstream functions like .effect_sizes_divergence()
  if (!is.data.frame(lm_results_df)) {
    stop("[calculate_lm_interaction_s4] Result from calculate_lm_interaction must be a data.frame, ",
         "but got: ", class(lm_results_df),
         ". This indicates an unexpected change in the underlying function.",
         call. = FALSE)
  }
  
  # Check if result is empty (0 rows OR 0 columns)
  if (nrow(lm_results_df) == 0 || ncol(lm_results_df) == 0) {
    warning("[calculate_lm_interaction_s4] Result is empty (", 
            nrow(lm_results_df), " rows, ", ncol(lm_results_df), " columns). ",
            "This can occur with: low sample counts per condition, ",
            "insufficient signal, or model convergence issues. ",
            "Try: increasing samples, using higher lambda for data generation, ",
            "or checking colData grouping structure.",
            call. = FALSE)
    # Return empty results gracefully rather than error
    analysis@lm_results <- list(lm_interaction = data.frame())
    return(analysis)
  }
  
  required_cols <- c("gene", "adj_p_interaction")
  missing_cols <- setdiff(required_cols, colnames(lm_results_df))
  if (length(missing_cols) > 0) {
    stop("[calculate_lm_interaction_s4] Result missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Available columns: ", paste(colnames(lm_results_df), collapse = ", "),
         call. = FALSE)
  }

  # Store results with model_data
  if (is.list(analysis@lm_results) && "lm_interaction" %in% names(analysis@lm_results)) {
    # Preserve other lm_results components
    analysis@lm_results$lm_interaction <- lm_results_df
    if (!is.null(model_data)) {
      analysis@lm_results$lm_interaction_model_data <- model_data
    }
  } else {
    analysis@lm_results <- list(lm_interaction = lm_results_df)
    if (!is.null(model_data)) {
      analysis@lm_results$lm_interaction_model_data <- model_data
    }
  }

  # Track metadata
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    "calculate_lm_interaction"
  )

  # Save if output_file provided
  if (!is.null(output_file)) {
    lm_data <- NULL
    
    # Try multiple ways to extract data
    if (!is.null(analysis@lm_results)) {
      if (is.list(analysis@lm_results) && "lm_interaction" %in% names(analysis@lm_results)) {
        lm_int <- analysis@lm_results$lm_interaction
        if (is.data.frame(lm_int)) {
          lm_data <- lm_int
        } else if (is.list(lm_int) && "results" %in% names(lm_int)) {
          lm_data <- lm_int$results
        }
      }
    }
    
    if (!is.null(lm_data) && is.data.frame(lm_data)) {
      save_analysis_output(lm_data, output_file, object = analysis, verbose = verbose,
                           func_name = "calculate_lm_interaction_s4")
    }
  }

  analysis
}

# ============================================================================
# JACKKNIFE WRAPPER
# ============================================================================

#' Jackknife resampling with confidence intervals
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for jackknife. Default: 1.0.
#' @param verbose \code{logical}. Print jackknife results summary. Default: FALSE.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#'   If > 1 and multiple q-values provided, uses parallel PSOCK cluster.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for jackknife results table with estimates, influence, outliers),
#'   .rds (for entire S4 object). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with jackknife results in @jackknife_results.
#'
#' @details
#' Requires diversity results to exist first. Will error if
#' \code{calculate_diversity_s4()} has not been run.
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
#' }
#'
#' @examples
#' # Create test analysis with diversity pre-computed
#' analysis <- TSENAT:::.create_test_analysis(
#'   n_genes = 8,
#'   n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5)
#' )
#' 
#' # Run jackknife estimation
#' analysis <- jackknife_entropy_outliers_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' # Check jackknife results
#' names(jackKnife(analysis))
#'
#' @export
#' @importFrom utils write.table
jackknife_entropy_outliers_s4 <- function(analysis, q = NULL, verbose = FALSE, nthreads = NULL, output_file = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check prerequisite: diversity must be calculated
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # PARAMETER EXTRACTION using utility function
  q <- resolve_slot_param(q, analysis@config, "q_values", 1.0)
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  
  # Ensure q is numeric
  if (!is.numeric(q)) {
    stop("'q' must be numeric", call. = FALSE)
  }

  for (q_val in q) {
    # Check if diversity at this q-value exists
    div_key <- paste0("q_", formatC(q_val, format = "f", digits = 3))
    if (!(div_key %in% names(analysis@diversity_results))) {
      stop("Diversity not calculated for q=", q_val,
           ". Run calculate_diversity_s4(analysis, q=", q_val, ") first.",
           call. = FALSE)
    }

    # Get diversity result (SummarizedExperiment with diversity assay)
    div_result <- analysis@diversity_results[[div_key]]
    
    # Extract diversity matrix from the SummarizedExperiment
    if (is(div_result, "SummarizedExperiment")) {
      # Get the diversity assay (should be named "diversity" from calculate_diversity_s4)
      if (length(SummarizedExperiment::assays(div_result)) > 0) {
        div_matrix <- as.matrix(SummarizedExperiment::assay(div_result, 1))
      } else {
        stop("Diversity SE for q=", q_val, " has no assays", call. = FALSE)
      }
    } else {
      # If not a SE, assume it's already a matrix
      div_matrix <- as.matrix(div_result)
    }

    # Run jackknife - pass the diversity matrix as x
    tryCatch({
      result <- .jackknife_entropy_outliers(
        x = div_matrix,
        q = q_val,
        verbose = verbose,
        nthreads = nthreads,
        ...
      )

      # Store with key "q_X.XXX" (consistent 3 decimal formatting)
      jk_key <- paste0("q_", formatC(q_val, format = "f", digits = 3))
      analysis@jackknife_results[[jk_key]] <- result

      # Track metadata
      analysis@metadata$function_calls <- c(
        analysis@metadata$function_calls,
        paste0("jackknife_tsallis_entropy[q=", q_val, "]")
      )
    }, error = function(e) {
      stop("Jackknife computation failed for q=", q_val, ":\n", e$message,
           call. = FALSE)
    })
  }

  # Save if output_file provided
  if (!is.null(output_file)) {
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Convert jackknife results to data.frame for text output
      tryCatch({
        # Extract all jackknife results and convert to data.frame
        all_results <- list()
        for (jk_key in names(analysis@jackknife_results)) {
          jk_result <- analysis@jackknife_results[[jk_key]]
          
          # Handle both single result and list of results
          if (inherits(jk_result, "tsenat_jackknife")) {
            # Single result - wrap in list for uniform processing
            jk_result <- list(jk_result)
            names(jk_result) <- "gene1"
          }
          
          if (inherits(jk_result, "tsenat_jackknife_list")) {
            # Convert list of results to data.frame
            df_list <- lapply(names(jk_result), function(gene_name) {
              res <- jk_result[[gene_name]]
              data.frame(
                gene = gene_name,
                estimate = res$estimate,
                jackknife_se = res$jackknife_se,
                n_transcripts = res$n_transcripts,
                n_outliers = length(res$outlier_indices),
                outlier_indices = paste(res$outlier_indices, collapse = ","),
                outlier_threshold = res$outlier_threshold,
                outlier_cutoff_value = res$outlier_cutoff_value,
                q_value = res$q,
                normalized = res$norm,
                stringsAsFactors = FALSE
              )
            })
            
            df <- do.call(rbind, df_list)
            df$q_key <- jk_key  # Add q-value key for multi-q results
            rownames(df) <- NULL
            all_results[[jk_key]] <- df
          }
        }
        
        # Combine all results into single data.frame
        if (length(all_results) > 0) {
          output_df <- do.call(rbind, all_results)
          rownames(output_df) <- NULL
          
          # Determine separator based on file extension
          sep <- if (grepl("\\.csv$", tolower(output_file))) "," else "\t"
          
          write.table(
            output_df,
            file = output_file,
            sep = sep,
            quote = FALSE,
            row.names = FALSE
          )
          
          if (verbose) {
            message("[jackknife_entropy_outliers_s4] Results saved to ", output_file)
          }
        }
      }, error = function(e) {
        warning("[jackknife_entropy_outliers_s4] Could not write jackknife results to file: ",
                conditionMessage(e), call. = FALSE)
      })
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
      if (verbose) {
        message("[jackknife_entropy_outliers_s4] Analysis object saved to ", output_file)
      }
    }
  }
  analysis
}

# ============================================================================

#' Calculate divergence metrics and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value for divergence.
#'   If NULL, uses first q_value from @config$q_values if available, else defaults to 1.0.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param control_group \code{character} or \code{NULL}. Control group identifier for divergence comparison.
#'   If NULL, reads from \code{@config$control_group} if available.
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param method \code{character} or \code{NULL}. Statistical method for divergence calculation.
#'   If NULL, reads from \code{@config$method} if available.
#' @param bootstrap \code{logical}. Whether to compute bootstrap confidence intervals. Default: FALSE.
#'   If not specified, reads from \code{@config$bootstrap} if available.
#' @param nboot \code{numeric} or \code{NULL}. Number of bootstrap replicates. Default: NULL.
#'   If NULL, reads from \code{@config$nboot} if available.
#' @param seed \code{numeric} or \code{NULL}. Random seed for reproducibility. Default: NULL.
#'   If NULL, reads from \code{@config$seed} if available.
#' @param progress \code{logical}. Show progress bar during computation. Default: FALSE.
#' @param ... Additional arguments passed to the base divergence function.
#'
#' @return Modified TSENATAnalysis with divergence metrics in @divergence_results
#'   (stored as list of data.frames or matrices).
#'
#' @details
#' Requires diversity results from calculate_diversity_s4() as prerequisite.
#' 
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' @examples
#' # Create and run divergence analysis with strong signal
#' # Create test analysis with diversity pre-computed
#' analysis <- TSENAT:::.create_test_analysis(
#'   n_genes = 8,
#'   n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5)
#' )
#' 
#' # Calculate divergence across q-values
#' analysis <- calculate_divergence_s4(analysis, verbose = FALSE)
#' 
#' # Check divergence results  
#' head(divergence(analysis))
#'
#' @export
#' @importFrom utils write.table
calculate_divergence_s4 <- function(analysis, q = NULL, verbose = TRUE, nthreads = NULL, 
                                    output_file = NULL, control_group = NULL, paired = FALSE, 
                                    method = NULL, bootstrap = FALSE, nboot = NULL, 
                                    seed = NULL, progress = FALSE, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  if (nrow(analysis@se) == 0) {
    stop("SummarizedExperiment in @se is empty", call. = FALSE)
  }

  # Check if diversity available
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # PARAMETER EXTRACTION using utility functions
  q <- resolve_slot_param(q, analysis@config, "q_values", 1.0)
  control_group <- resolve_slot_param(control_group, analysis@config, "control_group", NULL)
  method <- resolve_slot_param(method, analysis@config, "method", "percentile")
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  nboot <- resolve_slot_param(nboot, analysis@config, "nboot", NULL)
  seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)
  
  # Replace q=0 with q=0.01 for practical approximation
  # (q=0 divergence always returns 0, which is mathematically correct but uninformative)
  if (is.vector(q)) {
    q[q == 0] <- 0.01
  } else if (is.numeric(q) && q == 0) {
    q <- 0.01
  }
  
  # Paired design flag - ensure it's always a logical value
  if (!isTRUE(paired) && !isFALSE(paired)) {
    if ("paired" %in% names(analysis@config)) {
      paired <- analysis@config$paired
    } else {
      paired <- FALSE
    }
    if (!is.logical(paired) || is.na(paired)) paired <- FALSE
  }
  
  # Bootstrap parameters - ensure it's always a logical value
  if (!isTRUE(bootstrap) && !isFALSE(bootstrap)) {
    if ("bootstrap" %in% names(analysis@config)) {
      bootstrap <- analysis@config$bootstrap
    } else {
      bootstrap <- FALSE
    }
    if (!is.logical(bootstrap) || is.na(bootstrap)) bootstrap <- FALSE
  }
  
  # Sanitize method if needed
  if (!is.null(method)) {
    method <- as.character(method[1])
  }
  if (is.null(method) || is.na(method)) {
    method <- "percentile"
  }

  # Run divergence calculation with extracted parameters
  # Include progress explicitly to prevent NA boolean operations
  args <- list(
    se = analysis@se,
    q = q,
    verbose = verbose,
    nthreads = nthreads,
    progress = progress  # Use parameter value, not hardcoded
  )
  
  # Add parameters if they are not NULL/FALSE
  if (!is.null(control_group)) {
    args$control_group <- control_group
  }
  
  # Use isTRUE to safely handle NA values
  if (isTRUE(paired)) {
    args$paired <- paired
  }
  
  if (!is.null(method)) {
    args$method <- method
  }
  
  # Use isTRUE to safely handle NA values
  if (isTRUE(bootstrap)) {
    args$bootstrap <- bootstrap
    # Add bootstrap parameters if bootstrap is TRUE
    if (!is.null(nboot)) {
      args$nboot <- nboot
    }
    if (!is.null(seed)) {
      args$seed <- seed
    }
  }
  
  # Merge with additional args (which may override values)
  args <- c(args, list(...))
  
  # Run divergence calculation
  result <- tryCatch({
    do.call(.calculate_divergence, args)
  }, error = function(e) {
    # More detailed error handling
    stop("Divergence calculation failed [NA boolean likely in: ",
         paste0(names(args), collapse = ", "), "]:\n", 
         e$message, call. = FALSE)
  })

  # =========================================================================
  # VALIDATE AND STORE RESULTS
  # =========================================================================
  # Store results in list format (required by @divergence_results slot)
  if (is.null(result)) {
    warning("[calculate_divergence_s4] Result is NULL. Check .calculate_divergence() output.",
            call. = FALSE)
    analysis@divergence_results <- list()
  } else if (is(result, "SummarizedExperiment")) {
    # Wrap SummarizedExperiment in list with "divergence_se" key
    analysis@divergence_results <- list(divergence_se = result)
  } else if (is.list(result)) {
    analysis@divergence_results <- result
  } else if (is.data.frame(result) || is.matrix(result)) {
    # Wrap single result in list
    analysis@divergence_results <- list(main = result)
  } else {
    warning("[calculate_divergence_s4] Unexpected result type: ", class(result),
            ". Wrapping in list.",
            call. = FALSE)
    analysis@divergence_results <- list(result = result)
  }

  # Track metadata
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    paste0("calculate_divergence[q=", q, "]")
  )

  # Save if output_file provided
  if (!is.null(output_file)) {
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write divergence results as text table
      # Try to extract writeabledata from SummarizedExperiment or other formats
      tryCatch({
        if (length(analysis@divergence_results) > 0 && is(analysis@divergence_results[[1]], "SummarizedExperiment")) {
          # Extract assay data from first SummarizedExperiment
          write_data <- as.data.frame(assay(analysis@divergence_results[[1]]))
        } else {
          write_data <- as.data.frame(analysis@divergence_results)
        }
        write.table(write_data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
      }, error = function(e) {
        warning("[calculate_divergence_s4] Could not write divergence results to file: ", 
                conditionMessage(e), call. = FALSE)
      })
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }

  }

  analysis
}

# ============================================================================
# DETECT Q GENE INTERACTIONS WRAPPER
# ============================================================================

#' Detect q-dependent gene interactions
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric} or \code{NULL}. Q-values to test across spectrum.
#'   If NULL, auto-detects from \code{@config$q_values} or diversity results.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param paired \code{logical} or \code{NULL}. If TRUE, uses paired/blocked design 
#'   (requires \code{subject_col}). If NULL, reads from \code{@config$paired}.
#' @param subject_col \code{character} or \code{NULL}. Column name for subject/block identifiers
#'   (required when \code{paired=TRUE}). If NULL, reads from \code{@config$subject_col}.
#' @param condition_col \code{character}. Column name for sample grouping/condition (REQUIRED).
#'   Specifies the condition/treatment variable for testing q×condition interactions.
#'   Example: "sample_type", "treatment", "disease_status".
#' @param test \code{character}. Test method: "auto" (default), "kruskal-wallis" (unpaired),
#'   "friedman" (paired), or "art" (aligned rank transform).
#' @param multicorr \code{character}. Multiple testing correction: "hochberg" (default),
#'   "benjamini-yekutieli", "westfall-young", or "none".
#' @param entropy_col \code{character}. Column name containing entropy/diversity data.
#'   Default: "diversity".
#' @param q_col \code{character}. Column name containing q-values. Default: "q".
#' @param gene_col \code{character}. Column name containing gene identifiers. Default: "gene".
#' @param wy_randomizations \code{numeric} or \code{character}. Number of permutations for 
#'   Westfall-Young correction. Use "auto" to estimate from data. Default: 500.
#' @param nperm_mode \code{character}. Mode for automatic permutation estimation:
#'   "standard" (default), "conservative", or "interactive".
#' @param nthreads \code{numeric} or \code{NULL}. Number of parallel threads for computation.
#'   If NULL, reads from \code{@config$nthreads}.
#' @param verbose \code{logical}. If TRUE, prints progress messages. Default: FALSE.
#' @param ... Additional arguments passed to the base \code{.rank_test_q_condition()} function.
#'
#' @return Modified TSENATAnalysis with interaction results in @lm_results.
#'
#' @details
#' Analyzes how gene interactions change across q-value spectrum using rank-based
#' (Friedman/Kruskal-Wallis) or parametric (GAM) statistical tests.
#'
#' **Parameter resolution priority** (explicit > @config > default/auto-detect):
#' \itemize{
#'   \item \code{condition_col}: REQUIRED - must be explicitly provided
#'   \item \code{q}: explicit arg > \code{@config$q_values} > extract from diversity_results keys
#'   \item \code{paired}: explicit arg > \code{@config$paired} > FALSE (default)
#'   \item \code{subject_col}: explicit arg > \code{@config$subject_col}
#'   \item \code{multicorr}: explicit arg > \code{@config$multicorr} > "hochberg"
#'   \item \code{nthreads}: explicit arg > \code{@config$nthreads} > 1 (default)
#'   \item \code{test}: explicit arg > \code{@config$test} > "auto" (auto-selection)
#'   \item \code{nperm_mode}: explicit arg > \code{@config$nperm_mode} > "standard"
#' }
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- calculate_diversity_s4(analysis, norm = TRUE)
#' 
#' # Test Q×Condition interaction (condition_col is REQUIRED)
#' analysis <- rank_test_q_condition_s4(analysis, condition_col = "condition", 
#'                                            multicorr = "hochberg")
#' head(lmResults(analysis)$q_interactions)
#'
#' @export
#' @importFrom utils write.table
# ============================================================================
# S4 WRAPPER: Detect Q×Condition Gene Interactions (Rank-Based Testing)
# ============================================================================
# Purpose:
#   Wrapper around .rank_test_q_condition() that manages TSENATAnalysis object.
#   Tests for genes with CONDITION-SPECIFIC q-dependent entropy patterns.
#   Tests whether the effect of q-values DIFFERS between experimental conditions.
# 
# Key Features:
#   - Q×Condition interaction: Tests if entropy patterns across q-values differ by condition
#   - Multi-q analysis: Combines diversity results for multiple q-values into
#     a single SummarizedExperiment for joint hypothesis testing
#   - Rank-based statistics: Kruskal-Wallis (unpaired) or Friedman (paired)
#   - Scheirer-Ray-Hare test: Two-way non-parametric ANOVA on ranks
#   - Multiple testing correction: Hochberg, Benjamini-Yekutieli, or
#     Westfall-Young permutation procedure
#   - AR(1) correlation handling: Westfall-Young preserves q-value correlations
#   - Effect sizes: Eta-squared (η²) for q×condition interactions
#
#   Mathematical Background:
#   Tests null hypothesis: H0 = "Gene entropy q-effect does NOT differ between conditions"
#   vs Alternative: H1 = "Gene entropy q-dependence is CONDITION-SPECIFIC"
# 
#   Example: Gene shows strong isoform switching (q-dependent entropy) in tumor cells
#   but NOT in healthy cells -> Identified as disease-relevant q-dependent gene.
# 
#   For condition-specific q-dependent genes:
#   - Condition A: Strong entropy variation across q (q-dependent)
#   - Condition B: Flat entropy profile across q (q-independent)
#   - Interaction: Condition-specific q-dependence pattern reveals biological process
# ============================================================================
rank_test_q_condition_s4 <- function(
    analysis, 
    condition_col,
    q = NULL, 
    output_file = NULL,
    paired = NULL,
    subject_col = NULL,
    test = c("auto", "kruskal-wallis", "friedman", "art"),
    multicorr = c("hochberg", "benjamini-yekutieli", "westfall-young", "none"),
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    wy_randomizations = 500,
    nperm_mode = c("standard", "conservative", "interactive"),
    nthreads = NULL,
    verbose = FALSE,
    ...) {
  # Validate input is TSENATAnalysis S4 class
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Validate condition_col is provided (Q×Condition interaction is required)
  # Try to read from config if not explicitly provided
  if (missing(condition_col) || is.null(condition_col)) {
    # Try to get from config
    if (!is.null(analysis@config) && "condition_col" %in% names(analysis@config)) {
      condition_col <- analysis@config$condition_col
    } else {
      # Default to "condition" if still not found
      condition_col <- "condition"
    }
  }

  # ========================================================================
  # PREREQUISITE CHECK: Diversity must be pre-calculated
  # ========================================================================
  # .rank_test_q_condition() requires a SummarizedExperiment with:
  #   - assays: entropy values (genes × samples)
  #   - colData: q-values, condition_col, and optional subject information
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # ========================================================================
  # PARAMETER EXTRACTION FROM @config using centralized helpers
  # ========================================================================
  
  # Prepare dots for additional arguments
  dots <- list(...)
  
  # Match test and multicorr enums early
  if (!missing(test)) {
    test <- match.arg(test)
    dots$test <- test
  } else if ("test" %in% names(analysis@config)) {
    dots$test <- analysis@config$test
  }
  
  if (!missing(multicorr)) {
    multicorr <- match.arg(multicorr)
    dots$multicorr <- multicorr
  } else if ("multicorr" %in% names(analysis@config)) {
    dots$multicorr <- analysis@config$multicorr
  }
  
  if (!missing(nperm_mode)) {
    nperm_mode <- match.arg(nperm_mode)
    dots$nperm_mode <- nperm_mode
  } else if ("nperm_mode" %in% names(analysis@config)) {
    dots$nperm_mode <- analysis@config$nperm_mode
  }
  
  # Add column specification parameters
  dots$entropy_col <- entropy_col
  dots$q_col <- q_col
  dots$gene_col <- gene_col
  
  # Use resolve_slot_param for remaining parameters
  q <- resolve_slot_param(q, analysis@config, "q_values", NULL)
  paired <- resolve_slot_param(paired, analysis@config, "paired", NULL)
  subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col", NULL)
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  
  # Add resolved parameters to dots list
  if (!is.null(paired)) dots$paired <- paired
  if (!is.null(subject_col)) dots$subject_col <- subject_col
  # condition_col is REQUIRED and explicitly passed
  dots$condition_col <- condition_col
  dots$nthreads <- nthreads
  dots$wy_randomizations <- wy_randomizations
  dots$verbose <- verbose

  # ========================================================================
  # OPTIMIZATION: Use Cached Multi-Q SummarizedExperiment
  # ========================================================================
  # If calculate_diversity_s4() is called with multiple q-values in a single call,
  # it caches the combined SE in @metadata$diversity_combined for reuse.
  # 
  # This optimization bypasses expensive per-q recombination when available:
  #   - Per-q SEs are automatically cbind()ed horizontally
  #   - Column names include _q= suffix to distinguish q-values
  #   - colData is rbind()ed with preserved q-value information
  # 
  # Benefits: ~10-50x faster for multi-q analysis on large datasets
  
  # Check if we have combined diversity result cached (from lazy evaluation)
  if (!is.null(analysis@metadata$diversity_combined) && 
      is.list(analysis@metadata$diversity_combined) &&
      !is.null(analysis@metadata$diversity_combined$combined_se)) {
    
    # Use cached combined SE directly - it has correct structure and metadata
    se_multi_q <- analysis@metadata$diversity_combined$combined_se
    
    # Verify it's valid
    if (!is(se_multi_q, "SummarizedExperiment") || ncol(se_multi_q) == 0) {
      se_multi_q <- NULL
    }
  } else {
    se_multi_q <- NULL
  }
  
  # Fallback: Recombine Per-Q Results
  # ========================================================================
  # If cache not available (e.g., diversity calculated with separate calls),
  # manually combine per-q SummarizedExperiments into single SE:
  # 
  # Process:
  #   1. Extract q-values from diversity_results keys (format: "q_0.500", "q_1.000", etc.)
  #   2. Validate all SEs have same genes (rownames must match)
  #   3. cbind() all assay matrices with renamed columns (add _q=X.XXX suffix)
  #   4. rbind() all colData (with updated rownames matching combined columns)
  #   5. Combine rowData from first SE (genes are same across all q-values)
  # 
  # Result: Single SummarizedExperiment(genes × (samples per q × n_q))
  if (is.null(se_multi_q)) {
    
    # Step 1: Extract q-values from diversity_results keys (format: "q_0.5", "q_1.0", etc.)
    q_keys <- names(analysis@diversity_results)
    q_values_extracted <- as.numeric(sub("^q_", "", q_keys))
    q_values_extracted <- sort(q_values_extracted)

    # Step 2-4: Combine list of SEs (one per q-value) into single SE
    # Each SE has same genes (rows) but different q-value samples (columns)
    # cbind() the assay matrices, rbind() the colData
    combined_assay_list <- list()
    combined_coldata_list <- list()
    common_rownames <- NULL  # Track genes are same across all q-values

    for (key in sort(q_keys)) {
      se <- analysis@diversity_results[[key]]
      
      # Extract q-value from key
      q_val <- as.numeric(sub("^q_", "", key))
      
      # Ensure SE format
      if (!is(se, "SummarizedExperiment")) {
        if (is.matrix(se) || is.data.frame(se)) {
          se <- SummarizedExperiment(assays = list(entropy = as.matrix(se)))
        } else {
          stop("Diversity result for ", key, " is not a SummarizedExperiment or matrix",
               call. = FALSE)
        }
      }
      
      # Get assay data
      if (length(SummarizedExperiment::assays(se)) == 0) {
        stop("Diversity result for ", key, " has no assays", call. = FALSE)
      }
      assay_data <- SummarizedExperiment::assay(se, 1)
      
      # Extract and enforce consistent rownames
      assay_rownames <- rownames(assay_data)
      if (is.null(assay_rownames)) {
        assay_rownames <- paste0("gene_", seq_len(nrow(assay_data)))
      }
      if (is.null(common_rownames)) {
        common_rownames <- assay_rownames
      } else if (!identical(common_rownames, assay_rownames)) {
        # If rownames differ, use the first one and reorder/match
        if (length(common_rownames) == length(assay_rownames)) {
          assay_data <- assay_data[common_rownames, , drop = FALSE]
        } else {
          stop("Diversity result for ", key, " has different number of genes",
               call. = FALSE)
        }
      }
      rownames(assay_data) <- common_rownames
      
      # Append q-value to column names: distinguishes samples from different q-values
      # Format: "sample_01_q=0.500", "sample_02_q=1.000", etc.
      # This encoding enables downstream functions to parse q and map back to original samples
      orig_colnames <- colnames(assay_data)
      if (is.null(orig_colnames)) {
        orig_colnames <- paste0("sample_", seq_len(ncol(assay_data)))
      }
      unique_colnames <- paste0(orig_colnames, "_q=", q_val)
      colnames(assay_data) <- unique_colnames
      
      # Get colData - ensure q column is present
      cd <- as.data.frame(SummarizedExperiment::colData(se))
      if (nrow(cd) == 0) {
        cd <- data.frame(q = rep(q_val, ncol(assay_data)))
      } else if (!"q" %in% colnames(cd)) {
        cd$q <- q_val
      }
      rownames(cd) <- unique_colnames
      
      # Store for combination
      combined_assay_list[[key]] <- assay_data
      combined_coldata_list[[key]] <- cd
    }

    # Step 3: Combine all assays horizontally (cbind columns from different q-values)
    # This creates a matrix: genes × (sample_1_q_0.5, sample_2_q_0.5, ..., sample_1_q_1.0, ...)
    combined_assay <- do.call(cbind, combined_assay_list)
    
    # Step 4: Combine colData vertically (rbind from each q-value's colData)
    # Rownames already set to unique_colnames matching in final assay matrix
    combined_coldata_df <- do.call(rbind, combined_coldata_list)
    
    # Ensure colnames of combined_assay match rownames of combined_coldata_df
    colnames(combined_assay) <- rownames(combined_coldata_df)
    
    # Step 5: Get rowData from first diversity result
    # Genes (rows) are IDENTICAL across all q-values, so only need from first
    first_se <- analysis@diversity_results[[sort(q_keys)[1]]]
    
    # Ensure first_se is a SummarizedExperiment (handle edge cases)
    if (!is(first_se, "SummarizedExperiment")) {
      if (is.matrix(first_se) || is.data.frame(first_se)) {
        first_se <- SummarizedExperiment(assays = list(entropy = as.matrix(first_se)))
      }
    }
    
    # Extract rowData from first SE if it exists
    rd <- tryCatch({
      rd_temp <- SummarizedExperiment::rowData(first_se)
      if (nrow(rd_temp) > 0) rd_temp else NULL
    }, error = function(e) NULL)
    
    # Create combined SE
    se_multi_q <- SummarizedExperiment(
      assays = list(entropy = combined_assay),
      colData = combined_coldata_df
    )
    
    # Add rowData if available
    if (!is.null(rd) && nrow(rd) > 0) {
      SummarizedExperiment::rowData(se_multi_q) <- rd
    }
  }

  # Extract q-values for metadata tracking (available in both paths)
  if (exists("q_values_extracted") && !is.null(q_values_extracted)) {
    # Already defined in fallback path
    q_vals_for_tracking <- q_values_extracted
  } else {
    # Extract from colData in optimized path
    coldata_vals <- SummarizedExperiment::colData(se_multi_q)
    if (!is.null(coldata_vals) && "q" %in% colnames(coldata_vals)) {
      q_vals_for_tracking <- unique(as.numeric(coldata_vals$q))
      q_vals_for_tracking <- sort(q_vals_for_tracking)
    } else {
      q_vals_for_tracking <- NULL
    }
  }

  # ========================================================================
  # RUN CORE RANK-BASED Q-INTERACTION TESTING
  # ========================================================================
  # Delegate to .rank_test_q_condition() which performs:
  #   1. SummarizedExperiment → long-format data frame conversion
  #   2. Per-gene rank-based test selection (conditional on data characteristics)
  #   3. Westfall-Young permutation procedure (if multicorr="westfall-young")
  #   4. Multiple testing corrections (Hochberg, Benjamini-Yekutieli, none)
  #   5. Effect size computation (η²) and result classification
  # 
  # Use merged parameter dictionary: config values + explicit overrides
  result <- tryCatch({
    do.call(.rank_test_q_condition, c(list(data = se_multi_q), dots))
  }, error = function(e) {
    stop("q-interaction detection failed:\n", e$message,
         call. = FALSE)
  })

  # ========================================================================
  # STORE RESULTS IN TSENATAnalysis OBJECT
  # ========================================================================
  # Store results under @lm_results$q_interactions for accessor compatibility
  # This location allows other functions to retrieve results via:
  #   lmResults(analysis, "q_interactions")
  if (is.list(analysis@lm_results)) {
    analysis@lm_results$q_interactions <- result
  } else {
    analysis@lm_results <- list(q_interactions = result)
  }

  # Track function execution in audit trail
  # Enables reproducibility: know which function calls were run and in what order
  if (!is.null(q_vals_for_tracking)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("rank_test_q_condition[q=", paste(q_vals_for_tracking, collapse = ","), "]")
    )
  }

  # Save if output_file provided (using centralized output handler)
  if (!is.null(output_file)) {
    result_df <- as.data.frame(result)
    save_analysis_output(result_df, output_file, object = analysis, verbose = verbose,
                         func_name = "rank_test_q_condition_s4")
  }

  analysis
}

# ============================================================================
# CALCULATE DIFFERENCE WRAPPER
# ============================================================================

#' Calculate Difference Between Control and Treatment Groups (S4 Wrapper)
#'
#' S4 wrapper that operates on TSENATAnalysis objects to calculate differences
#' between control and treatment groups. Uses diversity results from \code{@diversity_results} 
#' slot (from \code{calculate_diversity_s4()}) and stores results in the \code{lm_results} slot.
#'
#' @param analysis A \code{TSENATAnalysis} object with diversity results in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value to use. If NULL, uses first diversity result or q=1.0.
#' @param control Character string specifying the control group identifier. If \code{NULL},
#'   attempts to retrieve from \code{analysis@config$control}.
#' @param condition_col \code{character} or \code{NULL}. Column name in colData identifying sample conditions.
#'   If NULL, reads from \code{@config$condition_col} or auto-detects.
#' @param method \code{character}. Difference calculation method. Default: "mean".
#'   If NULL, reads from \code{@config$method} if available.
#' @param test \code{character}. Statistical test type. Default: "wilcoxon".
#'   If NULL, reads from \code{@config$test} if available.
#' @param randomizations \code{numeric}. Number of randomizations. Default: 100.
#'   If NULL, reads from \code{@config$randomizations} if available.
#' @param pcorr \code{character}. P-value correction method. Default: "BH".
#'   If NULL, reads from \code{@config$pcorr} if available.
#' @param assayno \code{numeric}. Assay number to use. Default: 1.
#'   If NULL, reads from \code{@config$assayno} if available.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#'   If not specified, reads from \code{@config$verbose} if available.
#' @param paired \code{logical}. Whether data is paired. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param exact \code{logical}. Use exact test. Default: FALSE.
#'   If not specified, reads from \code{@config$exact} if available.
#' @param pseudocount \code{numeric}. Pseudocount for normalization. Default: 0.
#'   If NULL, reads from \code{@config$pseudocount} if available.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param seed \code{numeric} or \code{NULL}. Random seed. Default: NULL.
#'   If NULL, reads from \code{@config$seed} if available.
#' @param robust_loss_type \code{character}. Robust regression loss type. Default: "huber".
#'   If NULL, reads from \code{@config$robust_loss_type} if available.
#' @param robust_scale_method \code{character}. Robust scaling method. Default: "mad".
#'   If NULL, reads from \code{@config$robust_scale_method} if available.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Returns the modified \code{analysis} object invisibly with results stored in
#'   \code{analysis@lm_results$difference}.
#'
#' @details
#' **IMPORTANT:** Requires diversity results to exist first via \code{calculate_diversity_s4()}.
#' This wrapper extracts the diversity SummarizedExperiment from \code{@diversity_results},
#' not the raw input data in \code{@se}. This ensures you're comparing diversity values
#' between control and treatment groups, not raw abundance data.
#'
#' **Parameter resolution priority** (explicit > @config > auto-detect > error):
#' \itemize{
#'   \item \code{control}: Uses explicit arg, else \code{@config$control}, else error
#'   \item \code{nthreads}: Uses explicit arg, else \code{@config$nthreads}, else 1
#'   \item \code{condition_col} (sample grouping): Uses \code{@config$condition_col},
#'     else auto-detects from colData columns: "group", "sample_type", "condition"
#' }
#'
#' @importFrom utils write.table
#' @seealso
#' \code{\link{calculate_diversity_s4}} for computing diversity.
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis()
#' result <- calculate_difference_s4(analysis, control = "control")
#'
#' @export
calculate_difference_s4 <- function(analysis, control = NULL, q = NULL, condition_col = NULL,
                                    method = NULL, test = NULL, randomizations = NULL, pcorr = NULL,
                                    assayno = NULL, verbose = NULL, paired = FALSE, exact = FALSE,
                                    pseudocount = NULL, nthreads = NULL, seed = NULL, 
                                    robust_loss_type = NULL, robust_scale_method = NULL,
                                    output_file = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check prerequisites: diversity results must exist
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # Priority 1: Use explicit parameter
  # Priority 2: Use @config$control
  if (is.null(control)) {
    if ("control" %in% names(analysis@config)) {
      control <- analysis@config$control
    } else {
      stop("'control' must be specified (control group identifier) or set in config",
           call. = FALSE)
    }
  }

  # Determine which diversity result to use
  if (is.null(q)) {
    div_keys <- names(analysis@diversity_results)
    if (length(div_keys) == 0) {
      stop("No diversity results found in @diversity_results", call. = FALSE)
    }
    diversity_se <- analysis@diversity_results[[div_keys[1]]]
    q_used <- sub("^q_", "", div_keys[1])
  } else {
    q_key <- paste0("q_", formatC(q, format = "f", digits = 3))
    if (!(q_key %in% names(analysis@diversity_results))) {
      stop("Diversity not calculated for q=", q, 
           ". Available: ", paste(names(analysis@diversity_results), collapse = ", "),
           call. = FALSE)
    }
    diversity_se <- analysis@diversity_results[[q_key]]
    q_used <- q
  }

  # Determine condition column to use
  condition_col <- resolve_slot_param(condition_col, analysis@config, "condition_col", NULL)
  
  # Resolve remaining parameters using centralized handler
  method <- resolve_slot_param(method, analysis@config, "method", "mean")
  test <- resolve_slot_param(test, analysis@config, "test", "wilcoxon")
  randomizations <- resolve_slot_param(randomizations, analysis@config, "randomizations", 100)
  pcorr <- resolve_slot_param(pcorr, analysis@config, "pcorr", "BH")
  assayno <- resolve_slot_param(assayno, analysis@config, "assayno", 1)
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)
  pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount", 0)
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  robust_loss_type <- resolve_slot_param(robust_loss_type, analysis@config, "robust_loss_type", "huber")
  robust_scale_method <- resolve_slot_param(robust_scale_method, analysis@config, "robust_scale_method", "mad")
  
  # Parameters with logical defaults (check config if FALSE)
  if (!paired && "paired" %in% names(analysis@config)) {
    paired <- analysis@config$paired
  }
  if (!exact && "exact" %in% names(analysis@config)) {
    exact <- analysis@config$exact
  }
  
  # Optional parameters (may be NULL)
  seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)

  # Run difference calculation on diversity results
  # Note: diversity_se and its colData are already prepared by calculate_diversity_s4
  result <- tryCatch({
    .calculate_difference(
      x = diversity_se,
      condition_col = condition_col,
      control = control,
      method = method,
      test = test,
      randomizations = randomizations,
      pcorr = pcorr,
      assayno = assayno,
      verbose = verbose,
      paired = paired,
      exact = exact,
      pseudocount = pseudocount,
      nthreads = nthreads,
      seed = seed,
      robust_loss_type = robust_loss_type,
      robust_scale_method = robust_scale_method,
      ...
    )
  }, error = function(e) {
    stop("Difference calculation failed:\n", e$message,
         call. = FALSE)
  })

  # Store in lm_results under "difference" key
  if (is.list(analysis@lm_results)) {
    analysis@lm_results$difference <- result
  } else {
    analysis@lm_results <- list(difference = result)
  }

  # Track metadata
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    paste0("calculate_difference_s4[q=", q_used, ", control=", control, "]")
  )

  # Save if output_file provided (using centralized output handler)
  if (!is.null(output_file)) {
    diff_data <- if (!is.null(analysis@lm_results$difference$results)) {
      analysis@lm_results$difference$results
    } else {
      as.data.frame(analysis@lm_results$difference)
    }
    save_analysis_output(diff_data, output_file, object = analysis, verbose = verbose,
                         func_name = "calculate_difference_s4")
  }

  analysis
}

# ============================================================================
# TEST RANKBASED ASSUMPTIONS WRAPPER
# ============================================================================

#' Test rank-based method assumptions in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results stored
#'   in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value(s) to extract from diversity results.
#'   If NULL, uses the first available diversity result or q=1.0.
#' @param checks \code{character}. Which assumptions to test. Default includes:
#'   "exchangeability", "monotonicity", "consistency".
#' @param alpha \code{numeric}. Significance level for tests (default: 0.05).
#' @param ... Additional arguments (for future extensibility).
#'
#' @return Modified TSENATAnalysis object with assumption test results stored
#'   in \code{@metadata$rankbased_assumptions}.
#'
#' @details
#' This wrapper calls \code{.test_rankbased_assumptions()} on diversity data
#' extracted from the analysis object. Results include:
#'
#' \describe{
#'   \item{exchangeability}{Permutation test for temporal/spatial ordering effects}
#'   \item{monotonicity}{Spearman correlation stability across rows}
#'   \item{consistency}{Kendall's W concordance and ICC across samples}
#' }
#'
#' **Data Extraction Priority:**
#' 1. If q specified: uses diversity result for that q-value
#' 2. If q NULL: uses first available diversity result
#' 3. If no diversity results: extracts from cached combined result (\code{@metadata$diversity_combined})
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- test_rankbased_assumptions_s4(analysis, q = 1.0)
#' names(metadata(analysis, "rankbased_assumptions"))
#'
#' @export
#' @rdname test_rankbased_assumptions_s4
setGeneric("test_rankbased_assumptions_s4", function(analysis, q = NULL, 
           checks = c("exchangeability", "monotonicity", "consistency"),
           alpha = 0.05, ...) {
  standardGeneric("test_rankbased_assumptions_s4")
})

#' @rdname test_rankbased_assumptions_s4
setMethod(
  "test_rankbased_assumptions_s4",
  signature(analysis = "TSENATAnalysis"),
  function(analysis, q = NULL, 
           checks = c("exchangeability", "monotonicity", "consistency"),
           alpha = 0.05, ...) {
    
    # Validate inputs
    if (!methods::is(analysis, "TSENATAnalysis")) {
      stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }
    
    # Extract diversity data
    diversity_data <- NULL
    q_used <- q
    
    # If q is specified, try to get that specific q-value
    if (!is.null(q)) {
      q_key <- if (nchar(as.character(q)) > 3) {
        paste0("q_", round(q, 1))
      } else {
        paste0("q_", q)
      }
      
      if (q_key %in% names(analysis@diversity_results)) {
        div_se <- analysis@diversity_results[[q_key]]
        diversity_data <- assay(div_se, "diversity")
      }
    }
    
    # If q is NULL and multiple diversity results exist, combine all q-values
    if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) > 1) {
      entropy_list <- lapply(analysis@diversity_results, function(se) {
        mat <- assay(se, "diversity")
        if (!is.matrix(mat)) {
          mat <- as.matrix(mat)
        }
        return(mat)
      })
      
      # Use complete case analysis: keep only genes present in ALL q-value matrices
      # This is mathematically sound for rank-based tests (Friedman)
      # and follows best practices per scholarly literature (Springer Handbook, Permutation Tests)
      all_genes <- lapply(entropy_list, rownames)
      common_genes <- Reduce(intersect, all_genes)
      
      # Subset all matrices to common genes in same order
      entropy_list <- lapply(entropy_list, function(mat) {
        mat[common_genes, , drop = FALSE]
      })
      
      # Combine all matrices column-wise (genes x all samples across q-values)
      diversity_data <- do.call(cbind, entropy_list)
      # Keep natural column names from cbind to preserve structure
      q_used <- "all"
    }
    
    # If q is NULL and only one result, use it
    if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) == 1) {
      div_se <- analysis@diversity_results[[1]]
      diversity_data <- assay(div_se, "diversity")
      q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
    }
    
    # Fallback: use first diversity result
    if (is.null(diversity_data) && length(analysis@diversity_results) > 0) {
      div_se <- analysis@diversity_results[[1]]
      diversity_data <- assay(div_se, "diversity")
      if (is.null(q_used) || is.na(q_used)) {
        q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
      }
    }
    
    # Check that we have data
    if (is.null(diversity_data)) {
      stop("No diversity results found in analysis object. ",
           "Run calculate_diversity_s4() first.",
           call. = FALSE)
    }
    
    # Ensure we have a matrix
    if (!is.matrix(diversity_data)) {
      diversity_data <- as.matrix(diversity_data)
    }
    
    # Run assumptions test
    result <- tryCatch({
      .test_rankbased_assumptions(
        data = diversity_data,
        checks = checks,
        alpha = alpha
      )
    }, error = function(e) {
      stop("Rankbased assumptions test failed:\n", e$message,
           call. = FALSE)
    })
    
    # Store results
    analysis@metadata$rankbased_assumptions <- list(
      result = result,
      q_value_tested = q_used,
      checks_performed = checks,
      alpha_used = alpha,
      timestamp = Sys.time()
    )
    
    # Track function call
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("test_rankbased_assumptions_s4[q=", q_used, "]")
    )
    
    analysis
  }
)

# Helper function to extract q-value from key

.extract_q_from_key <- function(key) {
  # Extract numeric part from "q_X.X" format
  as.numeric(sub("^q_", "", key))
}

#' Plot Volcano and MA Grid from Differential Analysis Results (S4 Wrapper)
#'
#' S4 wrapper that extracts differential analysis results from a TSENATAnalysis
#' object and creates side-by-side volcano and MA plots for comparing control
#' and treatment groups.
#'
#' @param analysis \code{TSENATAnalysis} object with calculated differences
#'   (typically via \code{\link{calculate_difference_s4}}).
#' @param x_col \code{character}. Column name for x-axis in MA plot.
#'   Default: NULL (uses mean_difference if available, else mean fold-change).
#' @param padj_col \code{character}. Column name for adjusted p-values.
#'   Default: "padj" (the standard column name from calculate_difference).
#' @param label_thresh \code{numeric}. P-value threshold for labeling top genes.
#'   Genes with adjusted p-value below this threshold are labeled.
#'   Default: 0.1.
#' @param sig_alpha \code{numeric}. Significance threshold for coloring significant
#'   differences. Points with adjusted p-value below sig_alpha are highlighted.
#'   Default: 0.05.
#' @param top_n \code{integer}. Number of top genes (by significance) to label
#'   in volcano plot. Default: 5.
#' @param title_volcano \code{character}. Title for volcano plot.
#'   Default: NULL (no title).
#' @param title_ma \code{character}. Title for MA plot.
#'   Default: "Tsallis-based MA plot".
#' @param verbose \code{logical}. Print status messages. Default: FALSE.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#' @param width \code{numeric}. Width of the output plot in inches (default: 12).
#'   Only used if output_file is provided.
#' @param height \code{numeric}. Height of the output plot in inches (default: 7.2).
#'   Only used if output_file is provided.
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return
#' Invisibly returns a cowplot grid object containing both volcano and MA plots
#' combined side-by-side. If the plot cannot be created, returns NULL invisibly.
#'
#' @details
#' This wrapper extracts the difference results data frame from
#' \code{analysis@lm_results$difference} and passes it to the base
#' \code{.plot_volcano_ma_grid()} function.
#'
#' **Required Data:**
#' \itemize{
#'   \item Differential analysis must be computed via \code{calculate_difference_s4()}
#'   \item Results are stored in \code{analysis@lm_results$difference}
#' }
#'
#' **Expected Columns in Difference Results:**
#' \itemize{
#'   \item \code{genes} or \code{gene_id}: Gene identifiers
#'   \item \code{Normal_mean}, \code{Tumor_mean}: Group means (or equivalent controls/treatments)
#'   \item \code{mean_difference}: Calculated difference between groups
#'   \item \code{log2_fold_change}: Log2 fold-change values
#'   \item \code{raw_p_values} or \code{pvalue}: Un-adjusted p-values
#'   \item \code{adjusted_p_values} or \code{padj}: Adjusted p-values (default column used)
#' }
#'
#' **Volcano Plot Features:**
#' \itemize{
#'   \item X-axis: log2 fold-change or mean difference
#'   \item Y-axis: -log10(adjusted p-value)
#'   \item Top significant genes labeled
#'   \item Points colored by significance threshold
#' }
#'
#' **MA Plot Features:**
#' \itemize{
#'   \item X-axis: Average expression level (A)
#'   \item Y-axis: Log2 fold-change (M)
#'   \item Loess curve showing trend
#'   \item Significant changes highlighted
#' }
#'
#' @examples
#' # Create test analysis with diversity and differential results
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- calculate_difference_s4(analysis, control = "control",
#'   verbose = FALSE)
#'   
#' # Plot volcano and MA plots
#' p <- plot_volcano_ma_grid_s4(analysis, sig_alpha = 0.05, top_n = 3)
#' if (!is.null(p)) print(p)
#'
#' @seealso
#' \code{\link{calculate_difference_s4}} for computing differential analysis.
#'
#' @export
plot_volcano_ma_grid_s4 <- function(
    analysis,
    x_col = NULL,
    padj_col = "padj",
    label_thresh = 0.1,
    sig_alpha = 0.05,
    top_n = 5,
    title_volcano = NULL,
    title_ma = "Tsallis-based MA plot",
    verbose = FALSE,
    output_file = NULL,
    width = 12,
    height = 7.2,
    ...) {
  
  # Load visualization dependencies (ggplot2, cowplot, etc.)
  .load_visualization_deps()

  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

  # Validate input
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Extract difference results from S4 object
  if (is.null(analysis@lm_results) || !is.list(analysis@lm_results)) {
    stop("No LM results found in analysis@lm_results. ",
         "Run calculate_difference_s4() first.", call. = FALSE)
  }

  if (!("difference" %in% names(analysis@lm_results))) {
    stop("Difference results not found in analysis@lm_results$difference. ",
         "Run calculate_difference_s4() first.", call. = FALSE)
  }

  diff_df <- analysis@lm_results$difference

  if (!is.data.frame(diff_df) || nrow(diff_df) == 0) {
    stop("Difference results are empty or not a data frame", call. = FALSE)
  }

  # Auto-detect column names - padj_col with fallbacks, x_col for effect size
  actual_padj_col <- auto_detect_column(
    colnames(diff_df), analysis@config, "padj_col",
    c("padj", "adjusted_p_values", "pvalue"),
    verbose = verbose, param_name = "padj_col"
  )
  
  if (is.null(x_col)) {
    x_col <- auto_detect_column(
      colnames(diff_df), analysis@config, "x_col",
      c("mean_difference", "log2_fold_change", "effect_size"),
      verbose = verbose, param_name = "x_col"
    )
  }

  plot_obj <- tryCatch({
    .plot_volcano_ma_grid(
      diff_df = diff_df,
      x_col = x_col,
      padj_col = actual_padj_col,
      label_thresh = label_thresh,
      sig_alpha = sig_alpha,
      top_n = top_n,
      title_volcano = title_volcano,
      title_ma = title_ma,
      ...
    )
  }, error = function(e) {
    stop("[plot_volcano_ma_grid_s4]", conditionMessage(e),
         call. = FALSE)
  })

  if (verbose) {
    message("[plot_volcano_ma_grid_s4] Plot created successfully")
  }

  # Save plot to file if requested
  if (!is.null(output_file)) {
    save_analysis_output(plot_obj, output_file, object = analysis, verbose = verbose,
                         func_name = "plot_volcano_ma_grid_s4",
                         width = width, height = height)
  }

  return(invisible(plot_obj))
}

# ============================================================================
# CONCORDANCE WRAPPER - Compute Method Concordance (GAM vs Friedman/KW)
# ============================================================================

#' Compute concordance between two analysis methods in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with LM results (e.g., GAM).
#' @param gam_method \code{character}. Key for GAM/interaction results in \code{@lm_results}.
#'   Default: "q_interactions" (results from \code{rank_test_q_condition_s4})
#' @param friedman_method \code{character}. Key for Friedman/rank-based results in \code{@lm_results}.
#'   Default: "rankbased" (results from \code{test_rankbased_assumptions_s4})
#' @param gam_results \code{data.frame} or \code{NULL}. Optional GAM results data frame to store
#'   in the analysis object. If provided, automatically stored in \code{@lm_results} under the
#'   key specified by \code{gam_method}. Useful for importing external results or results 
#'   computed outside the S4 wrapper. Default: NULL (use existing results in analysis).
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @return Modified TSENATAnalysis object with concordance results stored in:
#'   \code{@metadata$method_concordance}:
#'   \describe{
#'     \item{comparison_df}{Data frame comparing results from both methods}
#'     \item{spearman_rho}{Spearman correlation between adjusted p-values}
#'     \item{high_confidence}{Genes with strong agreement}
#'     \item{agreement_table}{Contingency table of significant/non-significant calls}
#'     \item{gam_method}{Method name used for GAM analysis}
#'     \item{friedman_method}{Method name used for Friedman analysis}
#'     \item{timestamp}{When concordance was computed}
#'   }
#'
#' @details
#' Compares results from two different statistical methods (typically GAM for continuous 
#' and Friedman/Kruskal-Wallis for rank-based analysis) on the same data. Identifies:
#' - Genes significant in both methods (high confidence)
#' - Genes detected by one method only (potential false positives or method-specific signal)
#' - Spearman correlation of p-values (overall agreement trends)
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' # Note: compute_method_concordance_s4 requires results from both
#' # rank_test_q_condition_s4 and test_rankbased_assumptions_s4
#'
#' @aliases compute_method_concordance_s4
#' @export
setGeneric("compute_method_concordance_s4", function(analysis, ...) {
  standardGeneric("compute_method_concordance_s4")
})

#' @rdname compute_method_concordance_s4
setMethod("compute_method_concordance_s4", "TSENATAnalysis", function(
    analysis,
    gam_method = "q_interactions",
    friedman_method = "rankbased",
    gam_results = NULL,
    verbose = FALSE,
    output_file = NULL) {
  
  # ===================================================================
  # VALIDATION
  # ===================================================================
  
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  # If gam_results provided, store them in lmResults automatically
  if (!is.null(gam_results)) {
    if (!is.data.frame(gam_results)) {
      stop("gam_results must be a data.frame", call. = FALSE)
    }
    # Store GAM results in analysis@lm_results
    if (is.null(analysis@lm_results)) {
      analysis@lm_results <- list()
    }
    analysis@lm_results[[gam_method]] <- gam_results
    if (verbose) {
      message("[compute_method_concordance_s4] Stored GAM results as '", gam_method, "'")
    }
  }
  
  if (is.null(analysis@lm_results)) {
    stop("No LM results found in analysis@lm_results. Run rank_test_q_condition_s4() first.",
         call. = FALSE)
  }
  
  # Check for required methods
  if (!(gam_method %in% names(analysis@lm_results))) {
    available_methods <- paste(names(analysis@lm_results), collapse = ", ")
    stop("GAM method '", gam_method, "' not found in LM results. ",
                "Available: ", available_methods, call. = FALSE)
  }
  
  if (!(friedman_method %in% names(analysis@lm_results))) {
    available_methods <- paste(names(analysis@lm_results), collapse = ", ")
    stop("Friedman method '", friedman_method, "' not found in LM results. ",
                "Available: ", available_methods, call. = FALSE)
  }
  
  # Extract results
  gam_results_final <- analysis@lm_results[[gam_method]]
  friedman_results <- analysis@lm_results[[friedman_method]]
  
  # Validate they're data frames
  if (!is.data.frame(gam_results_final)) {
    stop("GAM results ('", gam_method, "') must be a data.frame", call. = FALSE)
  }
  
  if (!is.data.frame(friedman_results)) {
    stop("Friedman results ('", friedman_method, "') must be a data.frame", call. = FALSE)
  }
  
  # ===================================================================
  # COMPUTE CONCORDANCE
  # ===================================================================
  
  if (verbose) {
    message("[compute_method_concordance_s4] Computing concordance between ",
        gam_method, " and ", friedman_method)
  }
  
  # Call the standard function
  concordance_result <- tryCatch({
    .compute_method_concordance(gam_results_final, friedman_results)
  }, error = function(e) {
    stop("[compute_method_concordance_s4]", conditionMessage(e), call. = FALSE)
  })
  
  # ===================================================================
  # STORE RESULTS
  # ===================================================================
  
  analysis@metadata$method_concordance <- list(
    comparison_df = concordance_result$comparison_df,
    spearman_rho = concordance_result$spearman_rho,
    high_confidence = concordance_result$high_conf,
    agreement_table = concordance_result$agreement_table,
    gam_method = gam_method,
    friedman_method = friedman_method,
    timestamp = Sys.time()
  )
  
  # Track function call
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    paste0("compute_method_concordance_s4[", gam_method, " vs ", friedman_method, "]")
  )
  
  if (verbose) {
    message("[compute_method_concordance_s4] Concordance computed successfully")
    if (!is.na(concordance_result$spearman_rho)) {
      message("[compute_method_concordance_s4] Spearman corr = ", 
          round(concordance_result$spearman_rho, 3))
    }
  }
  
  # ===================================================================
  # SAVE TO FILE (if output_file provided)
  # ===================================================================
  
  if (!is.null(output_file)) {
    if (verbose) {
      message("[compute_method_concordance_s4] Writing results to: ", output_file)
    }
    saveRDS(analysis, file = output_file)
  }
  
  analysis
}
)

#' Plot Global Divergence q-Curve Across All Genes (S4 Wrapper)
#'
#' S4 wrapper that extracts divergence results from a TSENATAnalysis object
#' and visualizes the average Tsallis divergence across all genes (or specified
#' genes) as a function of q-value.
#'
#' @param analysis \code{TSENATAnalysis} object with divergence results
#'   (typically via \code{\link{calculate_divergence_s4}}).
#' @param gene \code{character}. Optional specific gene name to plot.
#'   If NULL, plots global divergence curve (aggregated across all genes).
#' @param n_genes \code{integer}. Number of top genes to plot when showing
#'   multi-gene spectra. Default is 4. Genes are sorted by p-value significance.
#' @param ncol \code{integer}. Number of columns in grid layout for multi-gene
#'   plots. Default is 2. Number of rows is automatically calculated.
#' @param metric \code{character}. Summary statistic for global curve:
#'   "median" (default) or "mean". Only used when gene = NULL.
#' @param variability_metric \code{character}. Error bar type for global curve:
#'   "iqr" (default) or "sd". Only used when gene = NULL.
#' @param use_pvalue_ranking \code{logical}. If TRUE, uses LM results to rank and
#'   display top n_genes by p-value significance. If FALSE (default), plots global
#'   divergence curve when gene = NULL. Default is FALSE.
#' @param output_file \code{character}. Optional file path to save the plot.
#'   If NULL, plot is returned but not saved.
#' @param width \code{numeric}. Plot width in inches. Default is 10.
#' @param height \code{numeric}. Plot height in inches. Default is 6.
#' @param verbose \code{logical}. Print status messages. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @return
#' Invisibly returns the file path if saved, otherwise the ggplot object.
#' If the plot cannot be created (missing data, ggplot2 not available),
#' returns NULL invisibly with an informative message.
#'
#' @details
#' This wrapper extracts the divergence SummarizedExperiment from
#' \code{analysis@divergence_results} and optionally the LM results from
#' \code{analysis@lm_results$lm_interaction} to pass to the base function.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Divergence must be computed via \code{calculate_divergence_s4()}
#'   \item \code{@divergence_results$divergence_se} or direct divergence SE
#' }
#'
#' **Modes:**
#' \itemize{
#'   \item \strong{Global mode} (gene = NULL): Shows median/mean divergence
#'         across all genes with variability bands
#'   \item \strong{Gene-specific mode} (gene specified): Shows divergence
#'         spectrum for a single named gene
#'   \item \strong{Top genes mode} (gene = NULL, lm_res provided): Shows
#'         top n_genes by significance
#' }
#'
#' @examples
#' # Plot 1: Global divergence spectrum across all genes
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1, 1.5))
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' p_global <- plot_divergence_spectrum_s4(analysis)
#' if (!is.null(p_global)) print(p_global)
#'
#' @seealso
#' \code{\link{calculate_divergence_s4}} for computing divergence.
#'
#' @export
plot_divergence_spectrum_s4 <- function(
    analysis,
    gene = NULL,
    n_genes = 4,
    ncol = 2,
    metric = c("median", "mean"),
    variability_metric = c("iqr", "sd"),
    use_pvalue_ranking = FALSE,
    output_file = NULL,
    width = 12,
    height = NULL,
    verbose = TRUE,
    ...) {

  # Load visualization dependencies (ggplot2, cowplot, etc.)
  .load_visualization_deps()

  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

  # Validate input
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Match metric and variability_metric arguments
  metric <- match.arg(metric)
  variability_metric <- match.arg(variability_metric)

  # Extract divergence SE from analysis object
  if (is.null(analysis@divergence_results)) {
    stop("Divergence results not found in analysis@divergence_results. ",
         "Run calculate_divergence_s4() first.", call. = FALSE)
  }

  # Handle both direct SE and wrapped "divergence_se" key
  divergence_results_se <- if (is.list(analysis@divergence_results) &&
                               "divergence_se" %in% names(analysis@divergence_results)) {
    analysis@divergence_results$divergence_se
  } else if (is(analysis@divergence_results, "SummarizedExperiment")) {
    analysis@divergence_results
  } else {
    stop("Invalid divergence_results structure. Expected SummarizedExperiment or list with 'divergence_se' key",
         call. = FALSE)
  }

  if (nrow(divergence_results_se) == 0 || ncol(divergence_results_se) == 0) {
    stop("Divergence SummarizedExperiment is empty", call. = FALSE)
  }

  # Extract LM results for lm_res parameter (optional)
  # Only use for ranking if use_pvalue_ranking = TRUE
  lm_res <- NULL
  if (use_pvalue_ranking && !is.null(analysis@lm_results) && is.list(analysis@lm_results)) {
    if ("lm_interaction" %in% names(analysis@lm_results)) {
      lm_res <- analysis@lm_results$lm_interaction
    } else if (length(analysis@lm_results) > 0) {
      lm_res <- analysis@lm_results[[1]]
    }
  }

  # Validate LM results if using multi-gene mode
  if (is.null(gene) && use_pvalue_ranking && !is.null(lm_res)) {
    if (!is.data.frame(lm_res) || nrow(lm_res) == 0) {
      if (verbose) {
        message("Note: Invalid LM results. Plotting global curve without gene ranking.")
      }
      lm_res <- NULL
    }
  }

  # Calculate height if not provided (based on grid layout)
  if (is.null(height)) {
    n_rows <- ceiling(n_genes / ncol)
    height <- 3 + (3.5 * n_rows)  # 3" base + 3.5" per row
  }

  # Create the plot using base function
  p <- tryCatch({
    .plot_divergence_spectrum(
      divergence_results_se = divergence_results_se,
      gene = gene,
      lm_res = lm_res,
      n_genes = n_genes,
      ncol = ncol,
      metric = metric,
      variability_metric = variability_metric,
      ...
    )
  }, error = function(e) {
    if (verbose) {
      message("plot_divergence_spectrum failed: ", e$message)
    }
    return(NULL)
  })

  # If plot creation failed, return NULL invisibly
  if (is.null(p)) {
    return(invisible(NULL))
  }

  # Save plot to file if requested
  if (!is.null(output_file)) {
    save_analysis_output(p, output_file, object = analysis, verbose = verbose,
                         func_name = "plot_divergence_spectrum_s4",
                         width = width, height = height)
  }

  # Return file path if saved, otherwise return plot
  invisible(p)
}

# ============================================================================
# PLOT WRAPPER - Plot Method Concordance Comparison
# ============================================================================

#' Plot method concordance results from TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with computed method concordance
#'   (from \code{compute_method_concordance_s4()}).
#' @param verbose \code{logical}. Print progress messages. Default: FALSE
#'
#' @return A ggplot/cowplot object showing:
#'   \describe{
#'     \item{Panel 1}{Scatter plot of -log10(p-values) comparing methods}
#'     \item{Panel 2}{Histogram of p-value distributions by method}
#'   }
#'
#' @details
#' Creates visualization of method concordance including:
#' - Comparison of significance across two methods (with color-coded agreement)
#' - P-value distribution histograms for both methods
#' - Significance threshold lines at p < 0.05
#'
#' Requires that \code{compute_method_concordance_s4()} has already been run
#' to populate \code{@metadata$method_concordance}.
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' # Note: compute_method_concordance_s4 requires additional LM and Friedman results
#' # For demo, we show that plot_method_concordance_s4 needs pre-computed concordance
#'
#' @aliases plot_method_concordance_s4
#' @export
setGeneric("plot_method_concordance_s4", function(analysis, verbose = FALSE) {
  standardGeneric("plot_method_concordance_s4")
})

#' @rdname plot_method_concordance_s4
setMethod("plot_method_concordance_s4", "TSENATAnalysis", function(analysis, verbose = FALSE) {
  
  # Load visualization dependencies (ggplot2, cowplot, etc.)
  .load_visualization_deps()
  
  # Validate that concordance results exist
  if (is.null(analysis@metadata$method_concordance)) {
    stop(
      "[plot_method_concordance_s4] No concordance results found in @metadata.\n",
      "  Please run compute_method_concordance_s4() first."
    )
  }
  
  concordance_results <- analysis@metadata$method_concordance
  
  # Extract comparison dataframe
  comparison_df <- concordance_results$comparison_df
  
  if (is.null(comparison_df) || nrow(comparison_df) == 0) {
    stop("[plot_method_concordance_s4] Concordance comparison_df is empty or missing.")
  }
  
  if (verbose) {
    message("[plot_method_concordance_s4] Plotting concordance for ",
        nrow(comparison_df), " genes")
    message("[plot_method_concordance_s4] Methods compared: ",
        concordance_results$gam_method, " vs ",
        concordance_results$friedman_method)
  }
  
  # Call standard plotting function
  plot_obj <- .plot_method_concordance(comparison_df)
  
  if (verbose) {
    message("[plot_method_concordance_s4] Plot generated successfully")
  }
  
  return(plot_obj)
})

#' Compute Effect Sizes from Divergence Results (S4 Wrapper)
#'
#' S4 wrapper for \code{.effect_sizes_divergence()} that extracts divergence and LM
#' results directly from a TSENATAnalysis object.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing divergence results
#'   (from \code{calculate_divergence_s4()}) and LM interaction results
#'   (from \code{calculate_lm_interaction_s4()}).
#'
#' @param significance_threshold \code{numeric}. Adjusted p-value threshold for
#'   filtering significant genes (default: 0.05).
#'
#' @param enrich_per_q_pattern \code{logical}. If TRUE, enriches results with per-q
#'   divergence patterns (default: TRUE).
#'
#' @param verbose \code{logical}. If TRUE, print diagnostic messages (default: TRUE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with effect size results stored via
#'   \code{metadata(analysis)$effect_sizes_divergence}. Returns the analysis object visibly
#'   to support piping and method chaining.
#'
#' @details
#' **Workflow steps:**
#' \describe{
#'   \item{Extracting}{Divergence results (via accessors) and LM interaction results}
#'   \item{Computing}{Effect sizes using standard \code{.effect_sizes_divergence()} function}
#'   \item{Storing}{Results as list with \code{interaction_results} (data.frame) and
#'     \code{validation_stats}}
#'   \item{Tracking}{Function call in analysis metadata}
#' }
#'
#' **Parameter resolution priority** (explicit > metadata > default):
#' \itemize{
#'   \item \code{significance_threshold}: Uses explicit arg, else \code{metadata(analysis)$significance_threshold},
#'     else 0.05
#'   \item \code{enrich_per_q_pattern}: Uses explicit arg, else \code{metadata(analysis)$enrich_per_q_pattern},
#'     else TRUE
#'   \item \code{verbose}: Uses explicit arg, else \code{metadata(analysis)$verbose}, else TRUE
#' }
#'
#' Results are accessed via: \code{metadata(analysis)$effect_sizes_divergence}
#'
#' @examples
#' # Setup: Create test analysis with divergence and LM interaction results
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5), 
#'   verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#'   
#' # Compute effect sizes from divergence results
#' analysis <- effect_sizes_divergence_s4(analysis, 
#'   significance_threshold = 0.05, verbose = FALSE)
#' 
#' # Access results using metadata accessor
#' effect_size_results <- metadata(analysis)$effect_sizes_divergence
#' if (!is.null(effect_size_results)) {
#'   cat("Effect sizes computed successfully\n")
#' }
#'
#' @seealso
#' \code{\link{calculate_divergence_s4}} for divergence wrapper,
#' \code{\link{calculate_lm_interaction_s4}} for LM interaction wrapper
#'
#' @export
#' @importFrom methods is
#' @importFrom utils write.table
effect_sizes_divergence_s4 <- function(
    analysis,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE,
    output_file = NULL,
    ...) {

  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check for required results
  if (length(divRes(analysis)) == 0) {
    stop("Divergence results required. Run calculate_divergence_s4() first.",
         call. = FALSE)
  }

  if (is.null(lmRes(analysis)) || length(lmRes(analysis)) == 0) {
    stop("LM results required. Run calculate_lm_interaction_s4() first.",
         call. = FALSE)
  }

  # =========================================================================
  # PARAMETER EXTRACTION using utility functions
  significance_threshold <- resolve_slot_param(significance_threshold, getConfig(analysis), 
                                               "significance_threshold", 0.05)
  enrich_per_q_pattern <- resolve_slot_param(enrich_per_q_pattern, getConfig(analysis), 
                                             "enrich_per_q_pattern", TRUE)
  verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)

  # =========================================================================
  # EXTRACT RESULTS FROM ANALYSIS OBJECT
  # =========================================================================
  
  # Extract divergence SE
  analysis_divres <- divRes(analysis)
  divergence_se <- if (is(analysis_divres, "SummarizedExperiment")) {
    analysis_divres
  } else if (is.list(analysis_divres) && "divergence_se" %in% names(analysis_divres)) {
    analysis_divres$divergence_se
  } else if (is.list(analysis_divres) && length(analysis_divres) > 0) {
    # Fallback: check if first element is SE
    analysis_divres[[1]]
  } else {
    NULL
  }

  if (is.null(divergence_se) || !is(divergence_se, "SummarizedExperiment")) {
    stop("Could not extract divergence SummarizedExperiment from divergence results",
         call. = FALSE)
  }

  # Extract LM results
  analysis_lmres <- lmRes(analysis)
  lm_res <- if (is.data.frame(analysis_lmres)) {
    analysis_lmres
  } else if (is.list(analysis_lmres) && "lm_interaction" %in% names(analysis_lmres)) {
    analysis_lmres$lm_interaction
  } else if (is.list(analysis_lmres) && length(analysis_lmres) > 0) {
    analysis_lmres[[1]]
  } else {
    NULL
  }

  if (is.null(lm_res) || !is.data.frame(lm_res)) {
    stop("Could not extract LM results data.frame from LM results",
         call. = FALSE)
  }

  if (verbose) {
    message("[effect_sizes_divergence_s4] Extracted:")
    message("  - Divergence SE: ", paste(dim(divergence_se), collapse = " x "))
    message("  - LM results: ", nrow(lm_res), " genes")
  }

  # =========================================================================
  # COMPUTE EFFECT SIZES
  # =========================================================================
  if (verbose) {
    message("[effect_sizes_divergence_s4] Computing effect sizes...")
  }

  # =========================================================================
  # ENSURE DIVERGENCE SE HAS GENE NAMES FOR EFFECT SIZES
  # =========================================================================
  # The effect_sizes_divergence base function requires gene_name in rowData
  # Create a proper mapping that expands genes to match transcript-level rows
  
  rd <- rowData(divergence_se)
  if (is.null(rd) || !("gene_name" %in% colnames(rd))) {
    if (is.null(rd)) {
      rd <- DataFrame(row.names = rownames(divergence_se))
    }
    
    gene_names_added <- FALSE
    
    # Try tx2gene mapping first
    base_se <- getSE(analysis)
    tx2gene <- metadata(base_se)$tx2gene
    if (!is.null(tx2gene) && nrow(tx2gene) > 0 && ncol(tx2gene) >= 2) {
      # Find transcript and gene columns by name
      col_names <- tolower(colnames(tx2gene))
      tx_col <- NULL
      gene_col <- NULL
      
      # Find transcript column (first priority: Transcript, second: Tx)
      for (cn in colnames(tx2gene)) {
        if (tolower(cn) %in% c("transcript", "tx")) {
          tx_col <- cn
          break
        }
      }
      
      # Find gene column (first priority: Gene, second: Gen)
      for (cn in colnames(tx2gene)) {
        if (tolower(cn) %in% c("gene", "gen")) {
          gene_col <- cn
          break
        }
      }
      
      if (!is.null(tx_col) && !is.null(gene_col)) {
        # Use tx2gene to map rownames to genes
        tx_in_divergence <- rownames(divergence_se)
        tx_vector <- as.character(tx2gene[[tx_col]])
        gene_vector <- as.character(tx2gene[[gene_col]])
        
        # Find matching indices
        match_idx <- match(tx_in_divergence, tx_vector)
        
        # Get gene names using matched indices
        gene_names_expanded <- gene_vector[match_idx]
        
        # Check if all transcripts have valid mappings (no NAs)
        if (all(!is.na(gene_names_expanded))) {
          rd$gene_name <- gene_names_expanded
          gene_names_added <- TRUE
          if (verbose) {
            message("[effect_sizes_divergence_s4] Added gene_name to rowData from tx2gene mapping")
          }
        }
      }
    }
    
    # If tx2gene didn't work, try direct assignment (for same nrow case)
    if (!gene_names_added && "gene" %in% colnames(lm_res)) {
      if (nrow(lm_res) == nrow(divergence_se)) {
        rd$gene_name <- as.character(lm_res$gene)
        gene_names_added <- TRUE
        if (verbose) {
          message("[effect_sizes_divergence_s4] Added gene_name to rowData via direct assignment")
        }
      }
    }
    
    # Store updated rowData
    rowData(divergence_se) <- rd
    
    # Validate that gene_name was successfully added without NAs
    rd_final <- rowData(divergence_se)
    if (!("gene_name" %in% colnames(rd_final)) || any(is.na(rd_final$gene_name))) {
      stop("[effect_sizes_divergence_s4] Failed to add valid gene_name column to divergence_se rowData. ",
           "Ensure tx2gene metadata is properly set with matching transcript IDs.",
           call. = FALSE)
    }
  }

  result <- tryCatch({
    .effect_sizes_divergence(
      lm_res = lm_res,
      divergence_results_se = divergence_se,
      significance_threshold = significance_threshold,
      enrich_per_q_pattern = enrich_per_q_pattern,
      verbose = verbose,
      ...
    )
  }, error = function(e) {
    stop("effect_sizes_divergence", e$message,
         call. = FALSE)
  })

  # =========================================================================
  # STORE RESULTS IN ANALYSIS OBJECT
  # =========================================================================
  # Store result list in metadata (not in dedicated slot since none exists)
  current_meta <- getMeta(analysis)
  if (is.null(current_meta)) {
    analysis@metadata <- list()
  }

  analysis@metadata$effect_sizes_divergence <- result

  # Track function call
  current_calls <- getMeta(analysis, "function_calls")
  if (is.null(current_calls)) {
    current_calls <- character(0)
  }
  analysis@metadata$function_calls <- c(
    current_calls,
    paste0("effect_sizes_divergence[threshold=", significance_threshold, "]")
  )

  if (verbose) {
    message("[effect_sizes_divergence_s4] Results stored in metadata")
    if (!is.null(result$interaction_results)) {
      message("  - Effect size results: ", nrow(result$interaction_results), " genes")
    }
  }

  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write effect size results as text table, filtering out all-NA columns
      results_df <- as.data.frame(result$interaction_results)
      
      # Filter out columns that are all NA
      all_na_cols <- colnames(results_df)[vapply(results_df, function(x) all(is.na(x)), FUN.VALUE = logical(1))]
      if (length(all_na_cols) > 0) {
        results_df <- results_df[, !colnames(results_df) %in% all_na_cols]
        if (verbose) {
          message("[effect_sizes_divergence_s4] Removed all-NA columns: ", paste(all_na_cols, collapse = ", "))
        }
      }
      
      write.table(results_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }

  }

  analysis
}

#' Plot Top Transcripts from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.plot_top_transcripts()} that extracts data directly from
#' a TSENATAnalysis object. Automatically retrieves the SummarizedExperiment and
#' LM results for visualizing transcript abundance across conditions.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing a processed
#'   SummarizedExperiment and optional LM interaction results.
#'
#' @param gene \code{character} or \code{NULL}. Gene identifier(s) to plot. If a vector 
#'   of multiple genes is provided, plots all of them. If NULL, automatically selects 
#'   the top genes from LM results based on \code{top_n} parameter (genes with lowest p-values).
#'   Default: NULL (auto-extract from lm_results).
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying
#'   group assignments (default: "sample_type").
#'
#' @param top_n \code{numeric}. Number of top transcripts to display for each
#'   condition (default: 3).
#'
#' @param output_file \code{character} or \code{NULL}. Path to save the plot as
#'   a PNG file. If NULL, saves to a temporary location (default: NULL).
#'
#' @param metric \code{character}. Method for ranking transcripts within genes.
#'   One of "median", "mean", "variance", or "iqr" (default: "median").
#'
#' @param width \code{numeric} or \code{NULL}. Output image width in inches. 
#'   If NULL, automatically calculated based on number of genes (default: ~13 inches per column).
#'
#' @param height \code{numeric} or \code{NULL}. Output image height in inches.
#'   If NULL, automatically calculated based on number of genes (default: ~10 inches per row + headers).
#'
#' @param fontsize \code{numeric}. Base font size for heatmap titles and labels 
#'   (default: 16pt). Automatically scaled for readability.
#'
#' @param cellwidth \code{numeric}. Width of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param cellheight \code{numeric}. Height of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param layout_ncol \code{numeric}. Number of heatmaps per row in fixed layout
#'   (default: 2). If NULL, uses adaptive layout based on transcript counts.
#'
#' @param use_tpm \code{logical}. If \code{TRUE}, uses TPM (Transcripts Per Million) 
#'   from metadata instead of raw counts (default: FALSE). TPM is normalized for sequencing 
#'   depth and is recommended for comparing expression across samples. Requires TPM data 
#'   in metadata from `build_analysis_s4()` or `.build_se()` with `tpm` parameter. 
#'   Raises error if TPM not available and `use_tpm = TRUE`.
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plotting (default: FALSE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Supported formats: .pdf, .png, .jpg. Default: NULL (no file output).
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return Invisibly returns the output file path (if `output_file` provided), or invisible(NULL) 
#'   if rendering to active graphics device. Graphics are rendered to the active grid device 
#'   for capture during vignette compilation.
#'
#' @details
#' This wrapper extracts the following from \code{analysis}:
#' \describe{
#'   \item{SummarizedExperiment}{From \code{analysis@se} containing transcript counts}
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction} for gene selection}
#' }
#'
#' If no gene is specified, the function automatically selects the top gene from
#' the LM results (lowest p-value). This simplifies visualization of genes with
#' significant q x condition interaction effects.
#'
#' @examples
#' # Plot 6: Top transcripts across groups
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 4, n_samples_per_group = 20,
#'   q_values = c(0.5, 1, 1.5))
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' plot_file <- plot_top_transcripts_s4(analysis, top_n = 3)
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for object structure
#'
#' @export
#' @importFrom methods is
plot_top_transcripts_s4 <- function(
    analysis,
    gene = NULL,
    condition_col = NULL,
    top_n = 4,
    output_file = NULL,
    metric = c("median", "mean", "variance", "iqr"),
    use_tpm = TRUE,
    width = NULL,
    height = NULL,
    fontsize = 16,
    cellwidth = 0,
    cellheight = 0,
    layout_ncol = 2,
    verbose = FALSE,
    ...) {

  # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
  .load_visualization_deps()

  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }

  se <- analysis@se
  if (!inherits(se, "SummarizedExperiment")) {
    stop("analysis@se must be a SummarizedExperiment object", call. = FALSE)
  }

  # =========================================================================
  # AUTO-DETECT condition_col
  # =========================================================================
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(se))
    condition_col <- auto_detect_column(
      cd_cols, analysis@config, "condition_col",
      c("condition", "sample_type", "group", "treatment"),
      verbose = verbose, param_name = "condition_col"
    )
  }

  # =========================================================================
  # EXTRACT LM RESULTS (for gene ranking if not specified)
  # =========================================================================
  lm_results_df <- NULL
  if (is.null(gene) && !is.null(analysis@lm_results)) {
    if ("lm_interaction" %in% names(analysis@lm_results)) {
      # Extract results data.frame from list structure
      if (is.data.frame(analysis@lm_results$lm_interaction)) {
        lm_results_df <- analysis@lm_results$lm_interaction
      } else if (is.list(analysis@lm_results$lm_interaction) &&
                 "results" %in% names(analysis@lm_results$lm_interaction)) {
        lm_results_df <- analysis@lm_results$lm_interaction$results
      }

      # Auto-select top gene from LM results
      if (!is.null(lm_results_df) && nrow(lm_results_df) > 0) {
        # Find p-value and gene columns
        p_col <- auto_detect_column(
          colnames(lm_results_df), analysis@config, "p_col",
          c("p_interaction", "padj", "pvalue", "p.value", "p_value"),
          verbose = FALSE, param_name = "p_col"
        )
        
        gene_col <- auto_detect_column(
          colnames(lm_results_df), analysis@config, "gene_col",
          c("gene", "gene_name", "gene_id"),
          verbose = FALSE, param_name = "gene_col"
        )

        if (!is.null(p_col) && p_col %in% colnames(lm_results_df)) {
          # Get top genes (sorted by p-value, select top_n)
          top_indices <- order(lm_results_df[[p_col]])[seq_len(min(top_n, nrow(lm_results_df)))]
          gene <- as.character(lm_results_df[top_indices, gene_col])

          if (verbose) {
            message("[plot_top_transcripts_s4] Auto-selected top ", length(gene), " genes from LM results")
            message("  ", paste(gene, collapse = ", "))
          }
        }
      }
    }
  }

  if (is.null(gene)) {
    stop("[plot_top_transcripts_s4] No gene specified and cannot auto-detect from LM results. ",
         "Provide gene explicitly.", call. = FALSE)
  }

  # =========================================================================
  # CALL BASE FUNCTION
  # =========================================================================
  if (verbose) {
    if (length(gene) > 1) {
      message("[plot_top_transcripts_s4] Calling .plot_top_transcripts() for genes: ", 
          paste(gene, collapse = ", "))
    } else {
      message("[plot_top_transcripts_s4] Calling .plot_top_transcripts() for gene: ", gene)
    }
  }

  plot_file <- tryCatch({
    .plot_top_transcripts(
      se = se,
      gene = gene,
      condition_col = condition_col,
      res = lm_results_df,
      top_n = top_n,
      output_file = output_file,
      metric = metric[1],  # Use first metric if multiple provided
      use_tpm = use_tpm,
      width = width,
      height = height,
      fontsize = fontsize,
      cellwidth = cellwidth,
      cellheight = cellheight,
      layout_ncol = layout_ncol,
      ...
    )
  }, error = function(e) {
    stop("[plot_top_transcripts_s4]", conditionMessage(e), call. = FALSE)
  })

  # Optionally save to file if output_file provided
  if (!is.null(output_file)) {
    if (verbose) {
      message("[plot_top_transcripts_s4] Plot saved to: ", output_file)
    }
  }

  # Always return the plot object (ggplot)
  # Knitr will auto-manage figure rendering
  invisible(plot_file)
}

#' Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)
#'
#' S4 wrapper that extracts effect size results from a TSENATAnalysis object
#' and generates a histogram visualization of Tsallis divergence effect sizes
#' across genes.
#'
#' @param analysis \code{TSENATAnalysis} object with effect sizes computed
#'   (typically via \code{\link{effect_sizes_divergence_s4}}).
#' @param threshold \code{numeric}. Effect size threshold for visual marking
#'   in the plot. Default is 0.1 (information-theoretic significance level).
#' @param output_file \code{character}. Optional file path to save the plot.
#'   If NULL, plot is returned but not saved.
#' @param width \code{numeric}. Plot width in inches. Default is 10.
#' @param height \code{numeric}. Plot height in inches. Default is 6.
#' @param verbose \code{logical}. Print status messages. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @return
#' Invisibly returns the file path if saved, otherwise the ggplot object.
#' If the plot cannot be created (missing data, ggplot2 not available),
#' returns NULL invisibly with an informative message.
#'
#' @details
#' This wrapper extracts the interaction results (with effect size columns)
#' from \code{analysis@metadata$effect_sizes_divergence$interaction_results}
#' and passes them to the base \code{.plot_divergence_distribution()} function.
#'
#' The function visualizes the distribution of effect sizes using the median
#' q-value's divergence (typically around q=1.0, close to Shannon entropy).
#' A red dashed line marks the information-theoretic significance threshold.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Effect sizes must be computed via \code{effect_sizes_divergence_s4()}
#'   \item \code{@metadata$effect_sizes_divergence$interaction_results} must
#'         contain columns matching pattern \code{effect_size_D_q*}
#' }
#'
#' @examples
#' # Plot 2: Distribution of effect sizes across genes
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1, 1.5))
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE)
#' p_dist <- plot_divergence_distribution_s4(analysis, verbose = FALSE)
#' if (!is.null(p_dist)) print(p_dist)
#'
#' @seealso
#' \code{\link{effect_sizes_divergence_s4}} for computing effect sizes.
#'
#' @export
plot_divergence_distribution_s4 <- function(
    analysis,
    threshold = 0.1,
    output_file = NULL,
    width = 12,
    height = 6,
    verbose = TRUE,
    ...) {

  # Load visualization dependencies (ggplot2, cowplot, etc.)
  .load_visualization_deps()

  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

  # Validate input
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Extract effect sizes from metadata
  if (is.null(analysis@metadata$effect_sizes_divergence)) {
    stop("Effect sizes not found in analysis@metadata$effect_sizes_divergence. ",
         "Run effect_sizes_divergence_s4() first.", call. = FALSE)
  }

  effect_sizes <- analysis@metadata$effect_sizes_divergence

  # Extract interaction results (with effect size columns)
  if (!is.list(effect_sizes) || is.null(effect_sizes$interaction_results)) {
    stop("Invalid effect size structure. Expected @metadata$effect_sizes_divergence$interaction_results",
         call. = FALSE)
  }

  interaction_results <- effect_sizes$interaction_results

  if (!is.data.frame(interaction_results) || nrow(interaction_results) == 0) {
    stop("interaction_results must be a non-empty data frame", call. = FALSE)
  }

  # Create the plot using base function
  p <- tryCatch({
    .plot_divergence_distribution(
      interaction_results = interaction_results,
      threshold = threshold,
      ...
    )
  }, error = function(e) {
    if (verbose) {
      message("plot_divergence_distribution failed: ", e$message)
    }
    return(NULL)
  })

  # If plot creation failed, return NULL invisibly
  if (is.null(p)) {
    return(invisible(NULL))
  }

  # Save to file if requested
  if (!is.null(output_file)) {
    save_analysis_output(p, output_file, object = analysis, verbose = verbose,
                         func_name = "plot_divergence_distribution_s4",
                         width = width, height = height)
  }

  # Always return the plot object (not file path)
  # Knitr will auto-manage figure rendering
  # File is saved separately if output_file provided
  invisible(p)
}


#' Prepare Gene Switching Tables from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.prepare_gene_switching_tables()} that extracts results
#' directly from a TSENATAnalysis object. Automatically retrieves LM results and
#' jackknife switching results from the analysis object slots.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing completed
#'   LM interaction and jackknife isoform switching analyses.
#'
#' @param n_top_genes \code{numeric} or \code{NULL}. Number of top genes
#'   (by adjusted p-value) to include in output tables. If \code{NULL},
#'   all genes with significant LM results are included.
#'
#' @param n_transcripts_per_gene \code{numeric}. Maximum number of transcripts
#'   to display per gene in the output tables (default: 10).
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during table preparation.
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{$summary_table}}{Gene-level summary with LM p-values and
#'           significant q-values}
#'     \item{\code{$transcript_tables}}{Named list of data.frames, one per gene,
#'           showing transcript-level switching metrics}
#'     \item{\code{$q_vector}}{Vector of q-values analyzed}
#'   }
#'
#' @details
#' This function extracts the following from \code{analysis}:
#' \describe{
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction$results}}
#'   \item{Jackknife results}{From \code{analysis@jackknife_results} or
#'         extracted from the switching analysis metadata}
#' }
#'
#' The wrapper automatically handles column detection and parameter extraction,
#' providing a simplified interface compared to the base function.
#'
#' @examples
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 8, n_samples_per_group = 20,
#'   q_values = c(0.5, 1.0, 1.5))
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' analysis <- jackknife_isoform_switching_s4(analysis, n_bootstrap = 5,
#'   verbose = FALSE)
#' tables <- prepare_gene_switching_tables_s4(analysis)
#' head(tables$summary_df)
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#'
#' @export
#' @importFrom methods is
#' @importFrom utils write.table
prepare_gene_switching_tables_s4 <- function(
    analysis,
    n_top_genes = NULL,
    n_transcripts_per_gene = 10,
    verbose = FALSE,
    output_file = NULL,
    ...) {
  
  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

  # Validation
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object")
  }
  
  # Extract LM results
  if (verbose) message("Extracting LM results from analysis object...")
  
  lm_results_list <- analysis@lm_results
  if (is.null(lm_results_list) || length(lm_results_list) == 0) {
    stop("No LM results found in analysis@lm_results. Run calculate_lm_interaction_s4() first.")
  }
  
  # Try to get lm_interaction results first, fallback to first available
  if (!is.null(lm_results_list$lm_interaction)) {
    if (is.data.frame(lm_results_list$lm_interaction$results)) {
      lm_res <- lm_results_list$lm_interaction$results
    } else if (is.data.frame(lm_results_list$lm_interaction)) {
      lm_res <- lm_results_list$lm_interaction
    } else {
      stop("Cannot extract LM results from analysis@lm_results$lm_interaction")
    }
  } else if (is.data.frame(lm_results_list)) {
    lm_res <- lm_results_list
  } else {
    stop("Cannot find LM results data.frame in analysis@lm_results")
  }
  
  if (verbose) message("  [OK] Extracted LM results with ", nrow(lm_res), " genes")
  
  # Extract jackknife/switching results
  if (verbose) message("Extracting jackknife switching results from analysis object...")
  
  jackknife_results_list <- analysis@jackknife_results
  if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
    stop("No jackknife results found in analysis@jackknife_results. Run jackknife_isoform_switching_s4() first.")
  }
  
  # Check if results are stored under "multi_q" key (jackknife_isoform_switching_multiq class)
  if ("multi_q" %in% names(jackknife_results_list)) {
    multi_q_object <- jackknife_results_list[["multi_q"]]
    
    # If it's a tsenat_isoform_switching_multiq object, use it directly
    if (inherits(multi_q_object, "tsenat_isoform_switching_multiq")) {
      multi_q_results <- multi_q_object
      if (verbose) {
        message("  [OK] Found multi-q results under 'multi_q' key with ", 
            length(multi_q_results), " q-values")
      }
    } else {
      stop("Element at jackknife_results$multi_q is not a tsenat_isoform_switching_multiq object")
    }
  } else {
    # Fallback: look for q-keyed results directly
    # Extract only the q-keyed results (filter by pattern "q_X_XX")
    q_key_pattern <- "^q_[0-9]+_[0-9]{2}$"
    q_keyed_results <- jackknife_results_list[grep(q_key_pattern, names(jackknife_results_list))]
    
    if (length(q_keyed_results) == 0) {
      stop("No q-keyed jackknife results found in analysis@jackknife_results. ",
           "Expected keys in format 'q_X_XX' (e.g., 'q_0_01', 'q_1_00') or 'multi_q'.")
    }
    
    # Wrap q-keyed results as a multi_q object for consistency
    multi_q_results <- q_keyed_results
    if (verbose) {
      message("  [OK] Found ", length(multi_q_results), " q-keyed results")
    }
  }
  
  if (verbose) message("  [OK] Extracted jackknife results with ", length(multi_q_results), " q-values")
  
  # Call base function with extracted parameters
  if (verbose) message("Calling .prepare_gene_switching_tables()...")
  
  result <- .prepare_gene_switching_tables(
    lm_res = lm_res,
    multi_q_results = multi_q_results,
    n_top_genes = n_top_genes,
    n_transcripts_per_gene = n_transcripts_per_gene,
    verbose = verbose,
    ...
  )
  
  if (verbose) message("[OK] Gene switching tables prepared successfully")
  
  # Track function call in metadata
  analysis@metadata$function_calls <- c(analysis@metadata$function_calls, 
                                        "prepare_gene_switching_tables_s4")
  analysis@metadata$function_timestamps <- c(analysis@metadata$function_timestamps,
                                             as.character(Sys.time()))
  
  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write as TSV/CSV
      write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
    } else {
      # Default to RDS for arbitrary objects
      saveRDS(result, file = output_file)
    }

  }
  
  return(result)
}

#' Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.plot_multiq_delta_influence_heatmaps()} that extracts results
#' directly from a TSENATAnalysis object. Automatically retrieves jackknife switching
#' results from the analysis object slots.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing completed
#'   jackknife isoform switching analysis across multiple q-values.
#'
#' @param n_genes \code{numeric}. Number of top genes to display in heatmaps
#'   (default: 4). Genes are ranked by LM p-values if available, otherwise
#'   by order of appearance in results.
#'
#' @param lm_results \code{data.frame} or \code{NULL}. Optional LM interaction
#'   results for ranking genes (default: NULL). If NULL, attempts to extract from
#'   \code{analysis@lm_results$lm_interaction}.
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plot generation (default: FALSE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A file path (character) to the saved heatmap PNG file, invisibly.
#'
#' @details
#' This function extracts the following from \code{analysis}:
#' \describe{
#'   \item{Jackknife results}{From \code{analysis@jackknife_results}, which should
#'         contain multi-q switching results keyed by q-value (e.g., "q_1.00")}
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction} if not
#'         explicitly provided, for ranking genes by significance}
#' }
#'
#' The wrapper automatically handles parameter extraction and provides a simplified
#' interface compared to the base function.
#'
#' @examples
#' # Plot 5: Multi-q delta influence (isoform switching) heatmaps
#' analysis <- TSENAT:::.create_test_analysis(n_genes = 4, n_samples_per_group = 20,
#'   q_values = c(0.5, 1, 1.5))
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' analysis <- jackknife_isoform_switching_s4(analysis, q = c(0.5, 1, 1.5),
#'   n_bootstrap = 5, verbose = FALSE)
#' heatmap_file <- plot_multiq_delta_influence_heatmaps_s4(analysis, n_genes = 2)
#'
#' @seealso
#' \code{\link{jackknife_isoform_switching_s4}} for computing switching results
#'
#' @export
#' @importFrom methods is
plot_multiq_delta_influence_heatmaps_s4 <- function(
    analysis,
    n_genes = 4,
    lm_results = NULL,
    verbose = FALSE,
    output_file = NULL,
    ...) {
  
  # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
  .load_visualization_deps()
  
  # Extract verbose parameter if not provided
  verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

  # Validation
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }
  
  if (verbose) message("Extracting jackknife switching results from analysis object...")
  
  # Extract jackknife/switching results
  jackknife_results_list <- analysis@jackknife_results
  if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
    stop("No jackknife results found in analysis@jackknife_results. ",
         "Run jackknife_isoform_switching_s4() first.", call. = FALSE)
  }
  
  # Check for multi-q result (stored under "multi_q" key when multiple q-values provided)
  if ("multi_q" %in% names(jackknife_results_list)) {
    switching_results <- jackknife_results_list$multi_q
    if (verbose) {
      message("  [OK] Found multi-q result with class: ", class(switching_results)[1])
    }
  } else {
    # Fallback: use all results as list (for single or multiple q-values)
    switching_results <- jackknife_results_list
    if (verbose) {
      message("  [OK] Using individual q-value results (", length(switching_results), " q-values)")
    }
  }
  
  if (verbose) {
    message("  Q-values: ", paste(names(switching_results), collapse = ", "))
  }
  
  # Extract LM results if not provided
  if (is.null(lm_results)) {
    if (verbose) message("Extracting LM results from analysis@lm_results...")
    
    lm_results_list <- analysis@lm_results
    if (!is.null(lm_results_list)) {
      if (!is.null(lm_results_list$lm_interaction)) {
        if (is.data.frame(lm_results_list$lm_interaction$results)) {
          lm_results <- lm_results_list$lm_interaction$results
        } else if (is.data.frame(lm_results_list$lm_interaction)) {
          lm_results <- lm_results_list$lm_interaction
        }
      }
      
      if (!is.null(lm_results)) {
        if (verbose) message("  [OK] Extracted LM results with ", nrow(lm_results), " genes")
      } else if (verbose) {
        message("  LM results not found; genes will be ranked by appearance")
      }
    }
  }
  
  if (verbose) message("Calling .plot_multiq_delta_influence_heatmaps()...")
  
  # Call base function with extracted parameters
  # Note: output_file parameter can be used to save heatmap as PNG file
  result <- .plot_multiq_delta_influence_heatmaps(
    switching_results = switching_results,
    n_genes = n_genes,
    lm_results = lm_results,
    verbose = verbose,
    output_file = output_file,
    ...
  )
  
  if (verbose) {
    message("[OK] Heatmap plot generated successfully")
  }

  invisible(result)
}

#' Plot GAM q-curves from TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and generates GAM q-curve plots
#' @param analysis \code{TSENATAnalysis} object with diversity and LM interaction results.
#' @param n_top \code{integer}. Number of top genes (by adjusted p-value) to plot 
#'   (default: 6). Only used if genes = NULL.
#'
#' @param genes \code{character} vector. Optional specific gene names to plot. 
#'   If provided, these genes are plotted directly regardless of significance.
#'   If NULL (default), top n_top significant genes are selected.
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#'   group assignments for samples. If NULL, attempts to auto-detect from @config
#'   (looks for "condition_col" or "condition"). If still NULL, defaults to "sample_type".
#'
#' @param sig_alpha \code{numeric}. Significance threshold for adjusted p-values 
#'   (default: 0.05). Only used if genes = NULL; filters lm_res to significant 
#'   genes before selecting top n.
#'
#' @param assay_name \code{character}. Name of the assay in se to extract 
#'   (default: "diversity").
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param width \code{numeric}. Width of the plot in inches (default: 12).
#'   Only used if output_file is not NULL.
#'
#' @param height \code{numeric} or \code{NULL}. Height of the plot in inches.
#'   Default: NULL (automatically calculated based on width and aspect ratio).
#'   Only used if output_file is not NULL.
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A single \code{ggplot} object with all selected genes arranged in a 
#'   grid layout. Can be saved with \code{ggplot2::ggsave()}.
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Extracts LM results from \code{@lm_results$lm_interaction} slot
#' 3. Detects condition_col from \code{@config} or uses default
#' 4. Calls \code{.plot_lm_interaction_gam()} with extracted parameters
#'
#' **Parameter Resolution (condition_col):**
#' \enumerate{
#'   \item Explicit \code{condition_col} parameter (highest priority)
#'   \item \code{@config$condition_col} if available
#'   \item \code{@config$condition} if available  
#'   \item Default: "sample_type"
#' }
#'
#' @seealso
#' \code{\link{calculate_lm_interaction_s4}} for running LM analysis on TSENATAnalysis.
#'
#' @examples
#' # Plot 3: GAM q-curves for genes with q-by-condition interactions
#' n_genes <- 4
#' n_samples_per_group <- 15
#' n_transcripts <- n_genes * 20
#' n_samples <- 2 * n_samples_per_group
#' q_vals <- seq(0.2, 2.5, by = 0.15)
#' 
#' # Create interaction-rich test data
#' control_counts <- matrix(rpois(n_transcripts * n_samples_per_group, lambda = 40),
#'   nrow = n_transcripts, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_transcripts * n_samples_per_group, lambda = 150),
#'   nrow = n_transcripts, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' counts <- pmax(counts, 30)
#' rownames(counts) <- paste0("TX_", 1:n_transcripts)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' 
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = counts),
#'   colData = S4Vectors::DataFrame(
#'     sample_id = colnames(counts),
#'     sample_type = rep(c("control", "treatment"), each = n_samples_per_group),
#'     row.names = colnames(counts)))
#' 
#' tx2gene_df <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_transcripts / n_genes))
#' S4Vectors::metadata(se)$tx2gene <- tx2gene_df
#' 
#' analysis <- TSENAT::TSENATAnalysis(se = se, config = list())
#' analysis <- calculate_diversity_s4(analysis, q = q_vals, verbose = FALSE, min_valid_frac = 0)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "sample_type", verbose = FALSE)
#' 
#' p_gam <- plot_lm_interaction_gam_s4(analysis, n_top = 2,
#'   condition_col = "sample_type", sig_alpha = 0.15)
#' if (!is.null(p_gam)) print(p_gam)
#'
#' @export
plot_lm_interaction_gam_s4 <- function(
  analysis,
  n_top = 6,
  genes = NULL,
  condition_col = NULL,
  sig_alpha = 0.05,
  assay_name = "diversity",
  output_file = NULL,
  width = 12,
  height = NULL,
  ...
) {
  # Load visualization dependencies (ggplot2, cowplot, mgcv, etc.)
  .load_visualization_deps()
  
  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  # Check that LM results exist
  if (is.null(analysis@lm_results) || is.null(analysis@lm_results$lm_interaction)) {
    stop("[plot_lm_interaction_gam_s4] No LM interaction results found in @lm_results$lm_interaction. ",
         "Run calculate_lm_interaction_s4() first.",
         call. = FALSE)
  }
  
  lm_res <- analysis@lm_results$lm_interaction
  
  if (!is.data.frame(lm_res)) {
    stop("[plot_lm_interaction_gam_s4] @lm_results$lm_interaction must be a data.frame",
         call. = FALSE)
  }
  
  # Check that diversity results exist (needed for SE reconstruction)
  if (length(analysis@diversity_results) == 0) {
    stop("[plot_lm_interaction_gam_s4] No diversity results found in @diversity_results. ",
         "Run calculate_diversity_s4() first.",
         call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT condition_col FROM @config IF NOT PROVIDED
  # =========================================================================
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(analysis@se))
    condition_col <- auto_detect_column(
      cd_cols, analysis@config, "condition_col",
      c("condition", "sample_type", "group", "treatment"),
      verbose = verbose, param_name = "condition_col"
    )
  }
  
  # Validate that condition_col exists in colData
  if (!(condition_col %in% colnames(colData(analysis@se)))) {
    stop("[plot_lm_interaction_gam_s4] Specified condition_col='", condition_col, 
         "' not found in colData. Available columns: ",
         paste(colnames(colData(analysis@se)), collapse = ", "),
         call. = FALSE)
  }
  
  # =========================================================================
  # RECONSTRUCT COMBINED DIVERSITY SE FOR PLOTTING
  # (Same approach as in calculate_lm_interaction_s4)
  # =========================================================================
  # Extract q-values from diversity_results keys
  q_keys <- names(analysis@diversity_results)
  q_computed <- as.numeric(sub("^q_", "", q_keys))
  
  # Reconstruct combined diversity SE with all q-values
  diversity_combined <- tryCatch({
    verbose <- if ("verbose" %in% names(analysis@config)) {
      analysis@config$verbose
    } else {
      FALSE
    }
    
    .calculate_diversity(
      x = analysis@se,
      q = sort(q_computed),
      norm = TRUE,
      verbose = verbose,
      bootstrap = FALSE
    )
  }, error = function(e) {
    stop("[plot_lm_interaction_gam_s4] Failed to reconstruct diversity SE:\n",
                conditionMessage(e), call. = FALSE)
  })
  
  # =========================================================================
  # EXTRACT model_data FROM STORED RESULTS
  # =========================================================================
  model_data <- NULL
  if ("lm_interaction_model_data" %in% names(analysis@lm_results)) {
    model_data <- analysis@lm_results$lm_interaction_model_data
  }
  
  # =========================================================================
  # CALCULATE HEIGHT IF NOT PROVIDED
  # =========================================================================
  if (is.null(height)) {
    # Estimate number of genes to be plotted
    if (!is.null(genes)) {
      n_genes_plot <- length(genes)
    } else {
      # Count significant genes
      if ("adj_p_interaction" %in% colnames(lm_res)) {
        sig_genes <- lm_res$adj_p_interaction <= sig_alpha
      } else if ("p_interaction" %in% colnames(lm_res)) {
        sig_genes <- lm_res$p_interaction <= sig_alpha
      } else {
        sig_genes <- rep(TRUE, nrow(lm_res))
      }
      n_genes_plot <- min(sum(sig_genes), n_top)
    }
    # Calculate height: 2 rows per 3-gene group, ~3.5 inches per row
    n_rows <- ceiling(n_genes_plot / 2)
    height <- 2 + (3.5 * n_rows)
  }
  
  # =========================================================================
  # CALL plot_lm_interaction_gam WITH RECONSTRUCTED DIVERSITY SE
  # =========================================================================
  result <- tryCatch({
    .plot_lm_interaction_gam(
      se = diversity_combined,
      lm_res = lm_res,
      condition_col = condition_col,
      n_top = n_top,
      genes = genes,
      sig_alpha = sig_alpha,
      assay_name = assay_name,
      model_data = model_data,
      output_file = output_file,
      width = width,
      height = height,
      ...
    )
  }, error = function(e) {
    stop("[plot_lm_interaction_gam_s4]", conditionMessage(e),
         call. = FALSE)
  })
  
  # =========================================================================
  # RETURN PLOT
  # =========================================================================
  # Track that plotting occurred
  if (is.list(analysis@metadata)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("plot_lm_interaction_gam_s4[n_top=", n_top, ", condition_col=", 
             condition_col, "]")
    )
  }
  
  # Save plot to file if requested
  if (!is.null(output_file)) {
    save_analysis_output(result, output_file, object = analysis, verbose = verbose,
                         func_name = "plot_lm_interaction_gam_s4",
                         width = width, height = height)
  }

  # Return the plot object directly (not the analysis object)
  result
}

#' Jackknife isoform switching analysis on TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and performs jackknife-based
#' isoform switching detection. Automatically extracts the SummarizedExperiment
#' and metadata from object slots.
#'
#' @param analysis \code{TSENATAnalysis} object containing:
#'   \itemize{
#'     \item \code{@se}: SummarizedExperiment with count data
#'     \item \code{@config}: Configuration metadata
#'   }
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#'   group assignments (default: "sample_type"). If NULL, attempts auto-detection.
#'
#' @param subject_col \code{character}. Optional column for paired/repeated measures design.
#'   If provided, enables paired analysis. Default: NULL (unpaired).
#'
#' @param gene_col \code{character}. Column name in rowData(se) or metadata identifying genes.
#'   Default: "gene".
#'
#' @param isoform_col \code{character}. Column name in rowData(se) or metadata identifying 
#'   isoforms/transcripts. Default: "transcript" or "isoform".
#'
#' @param q \code{numeric}. Tsallis entropy parameter(s) to analyze. Can be single value 
#'   or vector for multi-q analysis (default: 1).
#'
#' @param norm \code{logical}. Whether to use normalized diversity values 
#'   (default: TRUE).
#'
#' @param threshold \code{numeric}. Percentile threshold for detecting transcript switching 
#'   (default: 90). Transcripts with delta_influence >= threshold percentile are classified 
#'   as "switching".
#'
#' @param n_bootstrap \code{integer}. Number of bootstrap resamples for confidence 
#'   intervals (default: 1000).
#'
#' @param lm_results \code{data.frame}. Optional LM interaction results to filter genes.
#'   If provided, only genes in lm_results are analyzed.
#'
#' @param lm_p_threshold \code{numeric}. P-value threshold for filtering genes from 
#'   lm_results (default: 0.05).
#'
#' @param use_lm_fdr \code{logical}. If TRUE, uses adjusted p-values from lm_results 
#'   (default: TRUE).
#'
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for tables), or .rds (for S4 objects).
#'   When text format is specified, generates TWO files:
#'   \itemize{
#'     \item \code{output_file}: Gene-level summary (one row per gene with switching statistics)
#'     \item \code{output_file_transcripts.ext}: Transcript-level details (one row per transcript with p-values and FDR)
#'   }
#'   For .rds format, saves only the full analysis object.
#'   Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @return \code{TSENATAnalysis} object with jackknife results stored in \code{@jackknife_results}
#'   slot. Results are keyed by q-value (e.g., "q_1.00"). For multi-q analysis, multiple
#'   calls will accumulate results in the slot.
#'
#'   The analysis object is returned visibly to support method chaining:
#'   \preformatted{
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 0.5)
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 1.0)
#'   }
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Detects condition_col, gene_col, isoform_col from colData/rowData or @config
#' 3. Calls \code{.jackknife_isoform_switching()} with extracted parameters
#'
#' **Parameter Auto-Detection:**
#' \enumerate{
#'   \item \code{condition_col}: Uses explicit parameter, then @config, then "sample_type"
#'   \item \code{gene_col}: Uses explicit parameter, then looks for "gene" or "Gene"
#'   \item \code{isoform_col}: Uses explicit parameter, then looks for "transcript", "isoform", or "Isoform"
#' }
#'
#' @seealso \code{jackknife_isoform_switching} for the underlying implementation,
#' \code{\link{TSENATAnalysis}} for object structure.
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(42)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(200, 10), nrow = 20, ncol = 10,
#'     dimnames = list(paste0("TX", 1:20), paste0("Sample", 1:10)))),
#'   colData = data.frame(sample_id = paste0("Sample", 1:10),
#'     sample_type = rep(c("A", "B"), 5),
#'     pair = rep(1:5, 2), row.names = paste0("Sample", 1:10))
#' )
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = paste0("TX", 1:20),
#'   Gene = rep(paste0("GENE", 1:10), each = 2))
#' analysis <- TSENATAnalysis(se)
#' # Basic usage with single q-value
#' # analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # results <- jackknife_isoform_switching_s4(
#' #   analysis, condition_col = "sample_type", q = 1
#' # )
#'
#' @export
jackknife_isoform_switching_s4 <- function(
  analysis,
  condition_col = NULL,
  subject_col = NULL,
  gene_col = NULL,
  isoform_col = NULL,
  q = 1,
  norm = TRUE,
  threshold = 90,
  n_bootstrap = 1000,
  lm_results = NULL,
  lm_p_threshold = 0.05,
  use_lm_fdr = TRUE,
  output_file = NULL,
  verbose = FALSE,
  ...
) {
  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  se <- analysis@se
  
  if (!inherits(se, "SummarizedExperiment")) {
    stop("[jackknife_isoform_switching_s4] @se must be a SummarizedExperiment object",
         call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT condition_col (using centralized helper)
  # =========================================================================
  cd_cols <- colnames(colData(se))
  condition_col <- auto_detect_column(
    cd_cols, 
    config_list = analysis@config, 
    config_key = "condition_col",
    priority_candidates = c("sample_type", "condition", "group", "sample_group"),
    default_fallback = if (length(cd_cols) > 0) cd_cols[1] else NULL,
    verbose = verbose,
    param_name = "condition_col"
  )
  
  if (is.null(condition_col)) {
    stop(
      "[jackknife_isoform_switching_s4] Cannot auto-detect condition_col.\n",
      "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
      "SOLUTION: Set @config$condition_col or pass explicit parameter\n",
      "  Example: analysis@config$condition_col <- 'sample_type'\n",
      "  Or:      jackknife_isoform_switching_s4(analysis, condition_col = 'sample_type')\n",
      call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT gene_col AND isoform_col FROM rowData (using centralized helper)
  # =========================================================================
  rd <- if (!is.null(rowData(se)) && nrow(rowData(se)) > 0) {
    rowData(se)
  } else {
    NULL
  }
  
  rd_cols <- if (!is.null(rd)) colnames(rd) else character(0)
  
  # Detect gene_col
  if (is.null(gene_col)) {
    gene_col <- auto_detect_column(
      rd_cols,
      config_list = analysis@config,
      config_key = "gene_col",
      priority_candidates = c("gene_id", "gene", "Gene", "gene_name"),
      default_fallback = "gene",
      verbose = verbose,
      param_name = "gene_col"
    )
  }
  
  # Detect isoform_col
  if (is.null(isoform_col)) {
    isoform_col <- auto_detect_column(
      rd_cols,
      config_list = analysis@config,
      config_key = "isoform_col",
      priority_candidates = c("transcript_id", "transcript", "isoform", "Isoform", "tx_id"),
      default_fallback = "transcript",
      verbose = verbose,
      param_name = "isoform_col"
    )
  }
  
  # =========================================================================
  # EXTRACT LM_RESULTS IF PROVIDED VIA ANALYSIS OBJECT
  # =========================================================================
  if (is.null(lm_results) && !is.null(analysis@lm_results)) {
    # Try to extract LM results from analysis object
    if ("lm_interaction" %in% names(analysis@lm_results)) {
      lm_results <- analysis@lm_results$lm_interaction
      if (verbose) {
        message("[jackknife_isoform_switching_s4] Using LM interaction results from @lm_results")
      }
    }
  }
  
  # =========================================================================
  # CHECK Q-VALUES AVAILABILITY IN DIVERSITY RESULTS
  # =========================================================================
  if (!is.null(analysis@diversity_results) && length(analysis@diversity_results) > 0) {
    available_q_keys <- names(analysis@diversity_results)
    available_q <- as.numeric(sub("^q_", "", available_q_keys))
    available_q <- sort(unique(available_q))
    
    # Check if requested q-values are available
    q_vals <- if (is.numeric(q)) q else c(q)
    missing_q <- setdiff(q_vals, available_q)
    
    if (length(missing_q) > 0) {
      warning(
        "[jackknife_isoform_switching_s4] Requested q-values not in @diversity_results:\n",
        "  Requested: ", paste(q_vals, collapse = ", "), "\n",
        "  Available: ", paste(available_q, collapse = ", "), "\n",
        "  Missing:   ", paste(missing_q, collapse = ", "), "\n\n",
        "SOLUTION: Recompute diversity with all desired q-values before calling jackknife_isoform_switching_s4:\n",
        "  analysis <- calculate_diversity_s4(analysis, q = c(", paste(q_vals, collapse = ", "), 
        "), norm = TRUE)\n",
        "  analysis <- jackknife_isoform_switching_s4(analysis, q = c(", paste(q_vals, collapse = ", "), "))\n",
        call. = FALSE
      )
    }
    
    if (verbose && length(missing_q) == 0) {
      message("[jackknife_isoform_switching_s4] All requested q-values available in diversity results")
    }
  }

  # =========================================================================
  # VALIDATE n_bootstrap PARAMETER
  # =========================================================================
  if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1 || n_bootstrap < 1) {
    stop("'n_bootstrap' must be a positive integer", call. = FALSE)
  }
  
  if (n_bootstrap < 50) {
    warning(
      "n_bootstrap = ", n_bootstrap, " is less than the recommended minimum of 50. ",
      "Results may be unreliable. Consider increasing to at least 50-100 for stable estimates.",
      call. = FALSE
    )
  }
  
  # =========================================================================
  # CALL BASE FUNCTION
  # =========================================================================
  result <- tryCatch({
    .jackknife_isoform_switching(
      se = se,
      condition_col = condition_col,
      subject_col = subject_col,
      gene_col = gene_col,
      isoform_col = isoform_col,
      q = q,
      norm = norm,
      threshold = threshold,
      n_bootstrap = n_bootstrap,
      verbose = verbose,
      lm_results = lm_results,
      lm_p_threshold = lm_p_threshold,
      use_lm_fdr = use_lm_fdr
    )
  }, error = function(e) {
    stop("[jackknife_isoform_switching_s4] Jackknife analysis failed:\\n",
                conditionMessage(e), call. = FALSE)
  })
  
  # =========================================================================
  # STORE RESULTS IN ANALYSIS OBJECT (OPTIMIZED - vectorized q-value storage)
  # =========================================================================
  # Ensure q is a vector and validate
  if (!is.numeric(q) && !(is.vector(q) && all(vapply(q, is.numeric, FUN.VALUE = logical(1))))) {
    stop("[jackknife_isoform_switching_s4] 'q' must be numeric or numeric vector", call. = FALSE)
  }
  q_vals <- if (is.numeric(q)) q else as.numeric(c(q))
  
  # Check if result has multi-q class
  if (inherits(result, "tsenat_isoform_switching_multiq")) {
    # Multi-q result: store as-is (already has correct underscore-format keys from base function)
    for (q_key in names(result)) {
      analysis@jackknife_results[[q_key]] <- result[[q_key]]
    }
    # Also store the entire multi-q result object for plotting function access
    analysis@jackknife_results[["multi_q"]] <- result
    if (verbose) {
      message(sprintf("[jackknife_isoform_switching_s4] Stored multi-q result with keys: %s", paste(names(result), collapse = ", ")))
    }
  } else {
    # Store results for each q-value (vectorized - no explicit loop)
    # Use underscore format to match base function: paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
    q_keys <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q_vals)))
    
    # Store each result
    for (i in seq_along(q_keys)) {
      q_key <- q_keys[i]
      
      # Check if result is list with named q-values (should have underscore format now)
      if (is.list(result) && q_key %in% names(result)) {
        analysis@jackknife_results[[q_key]] <- result[[q_key]]
      } else if (length(q_vals) == 1) {
        # Single q-value: store result directly
        analysis@jackknife_results[[q_key]] <- result
      } else {
        # Multiple q-values: log warning if not found
        warning("[jackknife_isoform_switching_s4] Result for q=", q_vals[i], 
                " (key: ", q_key, ") not found in base function result.", call. = FALSE)
      }
      
      if (verbose) {
        message(sprintf("[jackknife_isoform_switching_s4] Stored results for %s", q_key))
      }
    }
  }
  
  # =========================================================================
  # TRACK FUNCTION CALL IN METADATA (OPTIMIZED - single paste())
  # =========================================================================
  if (is.list(analysis@metadata)) {
    # Pre-format all metadata in one call (more efficient)
    call_str <- sprintf(
      "jackknife_isoform_switching_s4[q=%s, condition_col=%s]",
      paste(q_vals, collapse = ","),
      condition_col
    )
    
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      call_str
    )
    analysis@metadata$function_timestamps <- c(
      analysis@metadata$function_timestamps,
      as.character(Sys.time())
    )
  }

  # =========================================================================
  # SAVE RESULTS IF output_file PROVIDED
  # =========================================================================
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write as text table (extract from jackknife results)
      tryCatch({
        output_ext <- sub("^.*\\.", ".", tolower(output_file))
        output_base <- sub(paste0(output_ext, "$"), "", output_file)
        transcript_file <- paste0(output_base, "_transcripts", output_ext)
        
        # Extract gene-level and transcript-level tables using helper
        gene_data <- extract_multiq_table(
          result, is_multiq = inherits(result, "tsenat_isoform_switching_multiq"),
          extract_fn = function(res, q_key) {
            if (!is.null(res$summary_table)) {
              return(res$summary_table)
            }
            return(NULL)
          },
          q_value_col = "q_value"
        )
        
        transcript_data <- extract_multiq_table(
          result, is_multiq = inherits(result, "tsenat_isoform_switching_multiq"),
          extract_fn = function(res, q_key) {
            if (!is.null(res$all_transcript_stats)) {
              return(res$all_transcript_stats)
            }
            return(NULL)
          },
          q_value_col = "q_value"
        )
        
        # Write gene-level to original file
        if (!is.null(gene_data)) {
          write.table(gene_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        
        # Write transcript-level to separate file
        if (!is.null(transcript_data)) {
          write.table(transcript_data, file = transcript_file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }, error = function(e) {
        warning("[jackknife_isoform_switching_s4] Could not write jackknife results to file: ",
                conditionMessage(e), call. = FALSE)
      })
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }
  }
  
  # Return modified analysis object (visibly for method chaining as documented)
  analysis
}

#' M-Estimation for Sample Quality (S4 Wrapper)
#'
#' S4 wrapper for \code{m_estimate} that performs robust M-estimation
#' on diversity results stored in a TSENATAnalysis object and stores results
#' back into the object.
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results
#'   (typically via \code{\link{calculate_diversity_s4}}).
#' @param condition_col \code{character}. Column name in sample metadata indicating
#'   condition/sample grouping. Auto-detected from \code{@config$condition_col}
#'   if available.
#' @param loss_type \code{character}. Type of loss function: "huber" (default),
#'   "tukey", or "lsq". Determines robustness vs efficiency trade-off.
#' @param scale \code{numeric}. Manual scale parameter. If NULL, estimated from data.
#' @param max_iter \code{integer}. Maximum iterations for M-estimation. Default: 50.
#' @param tol \code{numeric}. Convergence tolerance. Default: 1e-6.
#' @param paired \code{logical}. If TRUE, adjusts degrees of freedom for paired
#'   designs. Auto-detected from \code{@config$paired} if available.
#'   Default: FALSE.
#' @param pcorr \code{character}. P-value correction method. Default: "BH" (Benjamini-Hochberg).
#' @param q_combine_method \code{character}. How to collapse multi-q results:
#'   "mean" (default) or "median".
#' @param influence_threshold \code{numeric}. Threshold for classifying samples
#'   as high-influence. Default: 0.75.
#' @param scale_method \code{character}. Scale estimation method: "mad" (default),
#'   "proposal2", or "s-estimator".
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param verbose \code{logical}. Print status messages. Default: TRUE.
#'
#' @return
#' Modified TSENATAnalysis object with M-estimation results stored in
#' \code{analysis@metadata$m_estimate_results}. Contains data frame with
#' influence scores, robustness weights, entropy statistics, and QC classifications.
#' Returns visibly to support method chaining and piping.
#'
#' @details
#' This wrapper extracts diversity results from \code{analysis@diversity_results},
#' performs robust M-estimation on diversity values (entropy), and stores results
#' in the analysis object metadata.
#'
#' **M-Estimation:** Robust regression technique that down-weights outliers based
#' on their residuals. Useful for detecting low-quality samples that show
#' unusual diversity patterns.
#'
#' **M-Estimation Results include:**
#' \itemize{
#'   \item \code{sample_influence}: How much each sample affects the overall fit
#'   \item \code{robustness_weight}: Down-weighting factor (lower = more outlying)
#'   \item \code{entropy_mean}: Average entropy for the sample
#'   \item \code{entropy_sd}: Entropy variability within the sample
#'   \item \code{Status}: QC Classification ("OK" or "Flag for QC" based on influence_threshold)
#' }
#'
#' **Parameter resolution priority** (explicit > @config > error):
#' \itemize{
#'   \item \code{samples}: Uses explicit arg, else \code{@config$condition_col},
#'     else error. Note: despite parameter name 'samples', maps to condition grouping column
#' }
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Diversity results must be computed via \code{calculate_diversity_s4()}
#'   \item Sample grouping column required in colData (auto-detected from @config$condition_col
#'     or via 'samples' parameter)
#' }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(42)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(200, 10), nrow = 20, ncol = 10,
#'     dimnames = list(paste0("TX", 1:20), paste0("Sample", 1:10)))),
#'   colData = data.frame(sample_id = paste0("Sample", 1:10),
#'     sample_type = rep(c("A", "B"), 5),
#'     pair = rep(1:5, 2), row.names = paste0("Sample", 1:10))
#' )
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = paste0("TX", 1:20),
#'   Gene = rep(paste0("GENE", 1:10), each = 2))
#' analysis <- TSENATAnalysis(se)
#' # First compute diversity
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' analysis <- m_estimate_s4(analysis, condition_col = "sample_type", loss_type = "huber", verbose = FALSE)
#'
#' @seealso
#' \code{\link{calculate_diversity_s4}} for computing diversity
#'
#' @export
#' @importFrom utils write.table
m_estimate_s4 <- function(
    analysis,
    condition_col = NULL,
    loss_type = "huber",
    scale = NULL,
    max_iter = 50,
    tol = 1e-6,
    paired = NULL,
    pcorr = "BH",
    q_combine_method = "mean",
    influence_threshold = 0.75,
    scale_method = "mad",
    output_file = NULL,
    verbose = FALSE) {

  # Validate input
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check for diversity results
  if (is.null(analysis@diversity_results) || length(analysis@diversity_results) == 0) {
    stop("Diversity results not found. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }
  
  # Validate that diversity_results is a properly structured named list
  if (!is.list(analysis@diversity_results) || is.null(names(analysis@diversity_results))) {
    stop("Diversity results must be a named list of SummarizedExperiment objects",
         call. = FALSE)
  }

  # Auto-detect condition_col if not provided
  if (is.null(condition_col)) {
    if ("condition_col" %in% names(analysis@config)) {
      condition_col <- analysis@config$condition_col
      if (is.null(condition_col) || !is.character(condition_col) || condition_col == "") {
        stop("@config$condition_col must be a non-empty character value",
             call. = FALSE)
      }
      if (verbose) {
        message(sprintf("Auto-detected 'condition_col' from config: %s", condition_col))
      }
    } else {
      cd_cols <- colnames(SummarizedExperiment::colData(analysis@diversity_results[[1]]))
      stop(
        "Sample grouping column not specified:\n",
        "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
        "SOLUTION: Set @config$condition_col or pass 'condition_col' parameter\n",
        "  Example: analysis@config$condition_col <- 'sample_type'\n",
        "  Or:      m_estimate_s4(analysis, condition_col = 'sample_type')\n",
        call. = FALSE)
    }
  } else if (!is.character(condition_col) || length(condition_col) != 1) {
    stop("'condition_col' must be a single character value", call. = FALSE)
  }

  # Auto-detect paired parameter from @config if not explicitly provided
  paired <- resolve_slot_param(paired, analysis@config, "paired", FALSE)
  
  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }

  # Extract diversity results - get first SE to access sample metadata
  diversity_se <- analysis@diversity_results[[1]]

  if (is.null(diversity_se) || nrow(diversity_se) == 0) {
    stop("Diversity SummarizedExperiment is empty", call. = FALSE)
  }

  # Verify condition_col exists
  sample_info <- SummarizedExperiment::colData(diversity_se)
  if (!(condition_col %in% colnames(sample_info))) {
    stop(
      "Column '", condition_col, "' not found in sample metadata.\n",
      "Available columns: ", paste(colnames(sample_info), collapse = ", "), "\n\n",
      "SOLUTION: Use a valid column name\n",
      "  Example: m_estimate_s4(analysis, condition_col = 'sample_type')\n",
      call. = FALSE)
  }

  # Combine all q-value diversity results into a single matrix
  # (m_estimate needs all diversity data in one SE)
  if (verbose) {
    message(sprintf("Combining %d q-value diversity results...", length(analysis@diversity_results)))
  }

  first_se <- analysis@diversity_results[[1]]
  combined_assay <- SummarizedExperiment::assay(first_se)
  combined_colnames <- colnames(first_se)

  # Add other q-values
  for (q_name in names(analysis@diversity_results)[-1]) {
    se_q <- analysis@diversity_results[[q_name]]
    combined_assay <- cbind(combined_assay, SummarizedExperiment::assay(se_q))
    combined_colnames <- c(combined_colnames, colnames(se_q))
  }

  # Update column names to reflect combined data
  colnames(combined_assay) <- combined_colnames
  single_colData <- SummarizedExperiment::colData(first_se)

  # Replicate colData for each q-value
  n_q_values <- length(analysis@diversity_results)
  combined_colData <- do.call(rbind, replicate(n_q_values, single_colData, simplify = FALSE))
  rownames(combined_colData) <- combined_colnames

  # Create combined SummarizedExperiment
  combined_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = combined_assay),
    colData = combined_colData
  )

  # Run m_estimate
  if (verbose) {
    message("Running M-estimation on combined diversity...")
  }

  m_est_results <- tryCatch({
    result <- .m_estimate(
      x = combined_se,
      samples = condition_col,
      loss_type = loss_type,
      scale = scale,
      max_iter = max_iter,
      tol = tol,
      paired = paired,
      pcorr = pcorr,
      q_combine_method = q_combine_method,
      influence_threshold = influence_threshold,
      scale_method = scale_method
    )
    result
  }, error = function(e) {
    stop("M-estimation failed:\n", e$message, call. = FALSE)
  })

  # Store results in metadata
  analysis@metadata$m_estimate_results <- m_est_results

  # Save to output file if provided
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1) {
      stop("'output_file' must be a character string (file path)", call. = FALSE)
    }
    # Use save_analysis_output for consistent handling (data frame for TSV, RDS for S4)
    save_analysis_output(m_est_results, output_file, object = analysis, verbose = verbose,
                         func_name = "m_estimate_s4")
  }

  # Track function call
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    paste0("m_estimate_s4[condition_col=", condition_col, ",loss_type=", loss_type, "]")
  )

  if (verbose) {
    message("M-estimation complete. Results stored in @metadata$m_estimate_results")
  }

  invisible(analysis)
}

#' Filter Low-Abundance Transcripts in a TSENATAnalysis Object
#'
#' S4 wrapper for \code{.filter_se()} that filters low-abundance transcripts
#' directly within a \code{TSENATAnalysis} object. This maintains the consistent
#' S4 workflow pattern where functions accept and return analysis objects.
#'
#' @param analysis A \code{TSENATAnalysis} S4 object containing the
#'   \code{SummarizedExperiment} to be filtered.
#'
#' @param stringency Character. Filtering stringency level:
#'   "strict" (most stringent, default for unpaired designs),
#'   "medium" (moderate, default for paired designs), or
#'   "lenient" (least stringent).
#'   Each keeps transcripts present in different percentages of samples with TPM >= median.
#'   Default is determined by study design (paired vs unpaired).
#'
#' @param min_samples Numeric. Minimum number of samples in which a transcript
#'   must be present (default: 5). Used as a secondary filter.
#'
#' @param verbose Logical. If TRUE, print filtering progress and summary statistics
#'   (default: FALSE).
#'
#' @return Invisibly returns the modified \code{analysis} object with filtered
#'   \code{SummarizedExperiment} in the \code{@se} slot. The filtering operation
#'   modifies the analysis object in-place while maintaining all other slots
#'   (results, metadata, etc.).
#'
#' @details
#' This wrapper applies \code{.filter_se()} to the SummarizedExperiment within
#' the TSENATAnalysis object. The function:
#'
#' 1. Extracts the SE from \code{analysis@se}
#' 2. Filters using \code{.filter_se()} with specified parameters
#' 3. Stores the filtered SE back in \code{analysis@se}
#' 4. Returns the modified analysis object invisibly
#'
#' **Important:** Filtering should be performed BEFORE computing diversity,
#' divergence, or LM interaction results. If called after analysis results
#' have been computed, those results will be based on unfiltered data and
#' may not align with the filtered SE dimensions.
#'
#' @seealso
#' \code{\link{build_analysis_s4}} for creating a new analysis object
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(42)
#' tx_counts <- matrix(sample(10:100, 400, replace = TRUE), nrow = 40, ncol = 10,
#'   dimnames = list(paste0("TX", 1:40), paste0("Sample", 1:10)))
#' se <- SummarizedExperiment(assays = list(counts = tx_counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = paste0("TX", 1:40), Gen = rep(paste0("GENE", 1:10), each = 4))
#' # Add sample metadata with pair column required by filter_analysis_s4
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   sample_id = paste0("Sample", 1:10),
#'   sample_type = rep(c("Control", "Treatment"), 5),
#'   pair = rep(1:5, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- filter_analysis_s4(analysis, stringency = "medium")
#'
#' @export
filter_analysis_s4 <- function(analysis, stringency = NULL, min_samples = 5L, verbose = FALSE) {
  # Validate input
  if (!inherits(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }

  # Extract SE from analysis
  se <- analysis@se

  # Apply filtering via .filter_se(colData is already preserved within filter_se)
  se_filtered <- .filter_se(
    se = se,
    stringency = stringency,
    min_samples = min_samples,
    verbose = verbose
  )

  # Store filtered SE back in analysis object
  analysis@se <- se_filtered

  # Return modified analysis object
  analysis
}

#' Build a Complete TSENATAnalysis Object
#'
#' Convenience wrapper that combines \code{.build_se()} and \code{TSENATAnalysis()}
#' into a single function call. This creates a complete analysis object ready for
#' Tsallis entropy computation and downstream analysis.
#'
#' @param readcounts A matrix or data.frame of transcript-level read counts with
#'   transcript IDs as row names and sample names as column names. Typically output
#'   from quantification tools (SALMON, kallisto, etc.).
#'
#' @param tx2gene Either:
#'   - A path to a GFF3 or GFF3.gz file containing transcript-to-gene mapping
#'   - A path to a TSV file with columns "Transcript" and "Gene"
#'   - A data.frame with transcript-to-gene mapping
#'
#' @param assay_name Character. Name for the assay (default: "counts").
#'
#' @param metadata Optional data.frame with sample metadata. Should have sample
#'   names as row names and metadata columns (e.g., sample_type, condition, etc.).
#'
#' @param tpm Optional matrix of transcript-level TPM values. If provided, will be
#'   stored in the SummarizedExperiment. Same dimensions as readcounts required.
#'
#' @param effective_length Optional numeric vector of transcript effective lengths
#'   (e.g., from SALMON). Length should match nrow(readcounts).
#'
#' @param config Optional list of configuration parameters to store in the
#'   TSENATAnalysis object. Useful for tracking analysis parameters.
#'
#' @return A \code{TSENATAnalysis} S4 object with:
#'   \item{@se}{The SummarizedExperiment containing transcript counts and metadata}
#'   \item{@config}{Analysis configuration (empty list or user-provided)}
#'   \item{@diversity_results}{Empty list (populated by calculate_diversity_s4())}
#'   \item{@divergence_results}{Empty list (populated by calculate_divergence_s4())}
#'   \item{@lm_results}{Empty list (populated by calculate_lm_interaction_s4())}
#'   \item{@jackknife_results}{Empty list (populated by jackknife functions)}
#'   \item{@plots}{Empty list (populated by plotting functions)}
#'   \item{@metadata}{Metadata with package version and creation timestamp}
#'
#' @details
#' This wrapper combines two steps into one:
#' \enumerate{
#'   \item Call \code{.build_se()} to create a SummarizedExperiment from transcript counts
#'   \item Wrap the result in \code{TSENATAnalysis()} to create the analysis object
#' }
#'
#' The returned object is ready for diversity analysis via \code{calculate_diversity_s4()}.
#'
#' If you need to inspect or filter the SummarizedExperiment before creating the
#' TSENATAnalysis object, call \code{.build_se()} and \code{TSENATAnalysis()} separately.
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for the S4 class structure
#' \code{\link{calculate_diversity_s4}} for computing Tsallis entropy
#'
#' @examples
#' # Create example transcript count data
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples <- 10
#'
#' # Generate count matrix
#' counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
#'                  nrow = n_isoforms, ncol = n_samples)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#'
#' # Create tx2gene mapping
#' tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#'
#' # Create sample metadata
#' metadata <- data.frame(
#'   condition = rep(c("control", "treatment"), each = 5),
#'   row.names = colnames(counts))
#'
#' # Build analysis object
#' analysis <- build_analysis_s4(
#'   readcounts = counts,
#'   tx2gene = tx2gene,
#'   metadata = metadata)
#'
#' # Verify the analysis object was created
#' analysis
#'
#' @export
build_analysis_s4 <- function(readcounts, tx2gene, assay_name = "counts",
                             metadata = NULL, tpm = NULL, effective_length = NULL,
                             config = list()) {
  # Build SummarizedExperiment
  se <- .build_se(
    readcounts = readcounts,
    tx2gene = tx2gene,
    assay_name = assay_name,
    metadata = metadata,
    tpm = tpm,
    effective_length = effective_length
  )

  # Ensure sample_id column exists in colData (required by TSENATAnalysis)
  # OPTIMIZATION: Only add if not already present
  if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
    SummarizedExperiment::colData(se)$sample_id <- colnames(se)
  }

  # Store metadata in config for later use (e.g., in calculate_lm_interaction_s4)
  if (!is.null(metadata)) {
    config$metadata <- metadata
  }

  # Wrap in TSENATAnalysis
  analysis <- TSENATAnalysis(se = se, config = config)

  return(analysis)
}
