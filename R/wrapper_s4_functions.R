# S4 Wrapper Functions for TSENAT Pipeline
# These functions provide S4-integrated alternatives to existing analysis
# functions. They extract input from TSENATAnalysis slots, run analysis,
# and store results back to appropriate slots.
#

# ============================================================================
# DIVERSITY WRAPPER
# ============================================================================

#' Calculate diversity and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for Tsallis entropy.
#'   If NULL, uses q_values from \code{analysis@config$q_values} if available, else defaults to 1.0.
#' @param ... Additional arguments passed to \code{\link{calculate_diversity}},
#'   including: norm, bootstrap, pseudocount, nthreads, what, verbose, etc.
#'
#' @return Modified TSENATAnalysis object with diversity results stored
#'   in \code{@diversity_results}, keyed by "q_X.X" format (e.g., "q_1.0").
#'
#' @details
#' This wrapper calls \code{calculate_diversity()} once per q-value, storing
#' results as SummarizedExperiment objects. It extracts key parameters from
#' \code{analysis@config} with priority resolution (explicit > \code{@config} > default).
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{q}{Priority 1 (explicit) > Priority 2 (\code{@config$q_values}) > Priority 3 (default: 1.0)\cr
#'     **Note:** If explicit q AND \code{@config$q_values} both provided, explicit wins.}
#'   \item{verbose}{Priority: explicit > \code{@config$verbose} > TRUE}
#'   \item{bootstrap}{Priority: explicit > \code{@config$bootstrap} > FALSE}
#'   \item{pseudocount}{Priority: explicit > \code{@config$pseudocount} > 0}
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
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
#' \dontrun{
#'   # With q-values and parameters in \code{@config}:
#'   config <- list(
#'     q_values = seq(0.1, 2, by=0.1),
#'     norm = "range",
#'     bootstrap = TRUE,
#'     nthreads = 4
#'   )
#'   analysis <- TSENATAnalysis(se, config = config)
#'   # All parameters come from \code{@config}:
#'   analysis <- calculate_diversity_s4(analysis)
#'   
#'   # Override \code{@config} parameters with explicit arguments:
#'   analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 2.0), 
#'                                      nthreads = 2)
#'   
#'   # Retrieve results with computation metadata:
#'   div_q1 <- diversity(analysis, q = 1.0)
#'   
#'   # Check what parameters were actually used:
#'   params_used <- analysis@config$last_diversity_run$parameters_used
#'   print(params_used)  # Shows actual values including nthreads=2
#'   
#'   # Check per-q metadata on specific result:
#'   computed_with <- attr(div_q1, "computed_with")
#'   print(computed_with$timestamp)  # When was this computed?
#' }
#'
#' @export
calculate_diversity_s4 <- function(analysis, q = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  if (nrow(analysis@se) == 0) {
    stop("SummarizedExperiment in @se is empty", call. = FALSE)
  }

  # Priority 1: Use explicit parameter if provided
  # Priority 2: Use @config$q_values if set
  # Priority 3: Default to 1.0
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values
    } else {
      q <- 1.0
    }
  }

  # Ensure q is numeric vector
  if (!is.numeric(q)) {
    stop("'q' must be numeric", call. = FALSE)
  }

  # ===================================================================
  # PARAMETER EXTRACTION FROM @config (with priority resolution)
  # ===================================================================
  # Priority: explicit argument > @config > function default
  dots <- list(...)
  
  # Extract verbose parameter
  verbose <- if ("verbose" %in% names(dots)) {
    dots$verbose
  } else if ("verbose" %in% names(analysis@config)) {
    analysis@config$verbose
  } else {
    TRUE  # Default
  }
  
  # Extract bootstrap parameter
  bootstrap <- if ("bootstrap" %in% names(dots)) {
    dots$bootstrap
  } else if ("bootstrap" %in% names(analysis@config)) {
    analysis@config$bootstrap
  } else {
    FALSE  # Default
  }
  
  # Extract pseudocount parameter
  pseudocount <- if ("pseudocount" %in% names(dots)) {
    dots$pseudocount
  } else if ("pseudocount" %in% names(analysis@config)) {
    analysis@config$pseudocount
  } else {
    0  # Default
  }
  
  # Extract nthreads parameter
  nthreads <- if ("nthreads" %in% names(dots)) {
    dots$nthreads
  } else if ("nthreads" %in% names(analysis@config)) {
    analysis@config$nthreads
  } else {
    1  # Default
  }
  
  # Extract normalization method
  norm <- if ("norm" %in% names(dots)) {
    dots$norm
  } else if ("norm" %in% names(analysis@config)) {
    analysis@config$norm
  } else {
    TRUE  # Default
  }
  
  # Extract what parameter (S for Tsallis entropy, D for diversity)
  what <- if ("what" %in% names(dots)) {
    dots$what
  } else if ("what" %in% names(analysis@config)) {
    analysis@config$what
  } else {
    "S"  # Default to entropy
  }
  
  # Extract metadata parameter (for .map_metadata() application to results)
  metadata <- if ("metadata" %in% names(dots)) {
    dots$metadata
  } else {
    NULL
  }

  # Remove extracted parameters from dots to avoid duplicate argument errors
  params_to_remove <- c("verbose", "bootstrap", "pseudocount", "nthreads", "norm", "what", "metadata")
  dots_filtered <- dots[!(names(dots) %in% params_to_remove)]

  # Use three decimal precision for q-value formatting
  # This ensures consistent formatting for all q-values (e.g., 0.100, 0.150, 2.000)
  q_decimals <- 3
  
  # OPTIMIZATION: Call calculate_diversity ONCE with all q-values
  # This is much faster than looping because it reuses gene/isoform calculations
  calc_args <- list(
    x = analysis@se,
    q = q,  # Pass entire vector, not individual values
    norm = norm,
    verbose = verbose,
    bootstrap = bootstrap,
    pseudocount = pseudocount,
    nthreads = nthreads,
    what = what
  )
  
  # Include metadata if provided (for .map_metadata() application)
  if (!is.null(metadata)) {
    calc_args$metadata <- metadata
  }
  
  # Add any additional parameters from dots
  calc_args <- c(calc_args, dots_filtered)
  
  result_df <- do.call(calculate_diversity, calc_args)


  # Extract and store results for each q-value
  
  # For multi-q results: parse q-values from column names more robustly
  col_q_values <- NA
  if (length(q) > 1 && nrow(result_df) > 0) {
    # Extract q-values from column names (format: sample_q=X.XXX)
    col_names <- colnames(result_df)
    # Find columns matching _q= pattern
    q_col_indices <- grep("_q=", col_names)
    if (length(q_col_indices) > 0) {
      # Extract just the q-value part after _q=
      q_vals_str <- sub(".*_q=", "", col_names[q_col_indices])
      col_q_values <- as.numeric(q_vals_str)
      if (any(!is.na(col_q_values))) {
        names(col_q_values) <- q_col_indices
      } else {
        col_q_values <- NA
      }
    }
  }
  
  # ===================================================================
  # OPTIMIZATION (B: Lazy Conversion):
  # Store combined format for efficient downstream use
  # Per-q format created lazily on-demand by accessor
  # ===================================================================
  
  # ===================================================================
  # OPTIMIZATION (B: Lazy Conversion):
  # Create and cache combined SE with proper metadata structure
  # ===================================================================
  
  # ===================================================================
  # NO RECONSTRUCTION NEEDED: result_df IS already a SummarizedExperiment!
  # ===================================================================
  # calculate_diversity() returns a SummarizedExperiment directly.
  # Use it as-is for caching - no need to rebuild from data.frame
  combined_se <- result_df
  
  # Store combined SE in cache for optimization
  analysis@metadata$diversity_combined <- list(
    combined_result = result_df,
    combined_se = combined_se,
    q_values_computed = q,
    computation_params = list(
      norm = norm,
      verbose = verbose,
      bootstrap = bootstrap,
      pseudocount = pseudocount,
      nthreads = nthreads,
      what = what
    ),
    timestamp = Sys.time()
  )
  
  # Also populate diversity_results with per-q SEs (for accessor compatibility)
  # But use lazy evaluation: only extract what's needed if someone calls diversity(analysis, q=X)
  # For now, mark that per-q results are available lazily from combined
  analysis@metadata$.diversity_lazily_computed <- TRUE
  
  for (q_val in q) {
    tryCatch({
      # Handle empty results early
      if (nrow(result_df) == 0) {
        warning("[calculate_diversity_s4] Result for q=", q_val, " is empty (0 rows)",
                call. = FALSE)
        next
      }
      
      # Extract columns for this q-value
      q_cols <- NULL
      
      if (!is.na(col_q_values[1])) {
        # Use parsed q-values from column names
        # Allow small numerical tolerance for floating point comparison
        q_mask <- abs(col_q_values - q_val) < 1e-5
        if (any(q_mask)) {
          q_cols <- as.numeric(names(col_q_values)[q_mask])
        }
      }
      
      # Fallback: pattern matching
      if (is.null(q_cols) || length(q_cols) == 0) {
        q_formatted <- formatC(q_val, format = "f", digits = 3)
        q_pattern_str <- paste0("_q=", gsub("\\.", "\\\\.", q_formatted), "$")
        q_cols <- grep(q_pattern_str, colnames(result_df))
        
        # Fallback 2: try with 2 decimals
        if (length(q_cols) == 0) {
          q_formatted_fb <- formatC(q_val, format = "f", digits = 2)
          q_pattern_str <- paste0("_q=", gsub("\\.", "\\\\.", q_formatted_fb), "$")
          q_cols <- grep(q_pattern_str, colnames(result_df))
        }
        
        # Fallback 3: for single q-value
        if (length(q_cols) == 0 && length(q) == 1) {
          numeric_cols <- sapply(result_df, is.numeric)
          q_cols <- which(numeric_cols)
        }
      }
      
      # Validate that we found columns
      if (length(q_cols) == 0) {
        stop("[calculate_diversity_s4] No columns found for q=", q_val,
             ". Available columns: ", paste(head(colnames(result_df), 10), collapse=", "),
             call. = FALSE)
      }
      
      # Extract ONLY this q-value's data (for multi-q, exclude other samples' other q-values)
      result_subset <- result_df[, q_cols, drop = FALSE]

      # Convert data.frame to SummarizedExperiment (Gap 10B)
      if (is.data.frame(result_subset)) {
        # Extract numeric columns for diversity assay
        numeric_cols <- sapply(result_subset, is.numeric)
        if (!any(numeric_cols)) {
          # If no numeric columns, store as-is
          result_se <- result_subset
        } else {
          # Convert to SE with numeric columns as assay named "diversity"
          assay_data <- as.matrix(result_subset[, numeric_cols, drop = FALSE])
          
          # Create SE
          result_se <- SummarizedExperiment(assays = list(diversity = assay_data))
          rownames(result_se) <- rownames(result_subset)
          
          # Build colData with metadata columns if present (these are sample-level metadata)
          # numeric_cols is a named logical vector, so we can negate it directly
          if (length(numeric_cols) > 0) {
            metadata_mask <- !numeric_cols
            if (any(metadata_mask)) {
              cd <- result_subset[, metadata_mask, drop = FALSE]
              rownames(cd) <- colnames(assay_data)
              SummarizedExperiment::colData(result_se) <- cd
            }
          }
        }
      } else {
        result_se <- result_subset
      }

      # ===================================================================
      # VALIDATION: Ensure stored SE has compatible structure
      # ===================================================================
      if (is(result_se, "SummarizedExperiment")) {
        if (length(SummarizedExperiment::assays(result_se)) == 0) {
          stop(paste0("[calculate_diversity_s4] Converted SE for q=", q_val, 
                      " has no assays. Check diversity result structure."),
               call. = FALSE)
        }
        test_assay <- tryCatch({
          SummarizedExperiment::assay(result_se, 1)
        }, error = function(e) {
          stop(paste0("[calculate_diversity_s4] Cannot access assay in SE for q=", q_val,
                      ": ", conditionMessage(e)), call. = FALSE)
        })
        if (is.null(test_assay) || nrow(test_assay) == 0) {
          warning("[calculate_diversity_s4] Assay for q=", q_val, 
                  " is empty or NULL. This may cause issues downstream.",
                  call. = FALSE)
        }
      }

      # Apply colData from original SE to preserve sample metadata and pairing structure
      # This ensures consistency with the original SE's samples and any pairing info
      if (is(result_se, "SummarizedExperiment") && ncol(result_se) > 0) {
        original_coldata <- SummarizedExperiment::colData(analysis@se)
        if (!is.null(original_coldata) && nrow(original_coldata) == ncol(result_se)) {
          # Copy colData with same structure to preserve pairing and sample metadata
          SummarizedExperiment::colData(result_se) <- original_coldata
        }
      }
      
      # Store with key "q_X.XX..." (using consistent decimal formatting for all q-values)
      key <- paste0("q_", formatC(q_val, format = "f", digits = q_decimals))
      
      # GAP 4 FIX: Attach computation metadata to each result
      attr(result_se, "computed_with") <- list(
        q = q_val,
        norm = norm,
        verbose = verbose,
        bootstrap = bootstrap,
        pseudocount = pseudocount,
        nthreads = nthreads,
        what = what,
        timestamp = Sys.time()
      )
      
      analysis@diversity_results[[key]] <- result_se

      # Track in metadata
      analysis@metadata$function_calls <- c(
        analysis@metadata$function_calls,
        paste0("calculate_diversity[q=", q_val, "]")
      )
    }, error = function(e) {
      # GAP 6 FIX: More specific error reporting
      error_msg <- conditionMessage(e)
      if (bootstrap && grepl("bootstrap", error_msg, ignore.case = TRUE)) {
        stop(paste0("[calculate_diversity_s4] Bootstrap CI computation failed for q=", q_val, ": ", error_msg),
             call. = FALSE)
      } else {
        stop(paste0("[calculate_diversity_s4] Error computing diversity for q=", q_val, ":\n", error_msg),
             call. = FALSE)
      }
    })
  }

  # GAP 1 FIX: Update @config with actual parameters used
  analysis@config$last_diversity_run <- list(
    timestamp = Sys.time(),
    q_values_computed = q,
    num_q_values = length(q),
    parameters_used = list(
      norm = norm,
      verbose = verbose,
      bootstrap = bootstrap,
      pseudocount = pseudocount,
      nthreads = nthreads,
      what = what
    ),
    note = "Actual parameters used (save object and check this, not original @config)"
  )

  # GAP 7 FIX: Add audit trail for parallel processing
  if (nthreads > 1) {
    analysis@metadata$parallel_processing <- c(
      analysis@metadata$parallel_processing,
      paste0("calculate_diversity_s4: nthreads=", nthreads, " (", length(q), " q-values)")
    )
  }

  analysis
}

# ============================================================================
# LM INTERACTION WRAPPER
# ============================================================================

#' Calculate LM interactions and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param fdr_threshold \code{numeric}. FDR cutoff for significance.
#'   Default: 0.05.
#' @param formula \code{formula} or NULL. Reserved for future use.
#' @param method \code{character}. Statistical method (e.g., "lmm", "gam", "gee").
#'   If NULL, uses method from @config$method or defaults to "lmm".
#' @param ... Additional arguments passed to \code{\link{calculate_lm_interaction}},
#'   including: condition_col, paired, subject_col, multicorr, nthreads, etc.
#'
#' @return Modified TSENATAnalysis with results in @lm_results$lm_interaction.
#'
#' @details
#' Extracts diversity results from @diversity_results (prerequisite),
#' combines across q-values into single SummarizedExperiment,
#' then runs \code{calculate_lm_interaction()}.
#' 
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' @examples
#' \dontrun{
#'   # With parameters in \code{@config}:
#'   analysis <- TSENATAnalysis(se, 
#'     config = list(condition_col = "sample_type", paired = TRUE))
#'   analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0))
#'   analysis <- calculate_lm_interaction_s4(analysis)
#'   
#'   # Or with explicit parameters:
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "sample_type", method = "gam", paired = TRUE)
#' }
#'
#' @export
calculate_lm_interaction_s4 <- function(analysis, fdr_threshold = NULL, 
                                       formula = NULL, method = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

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
  # PARAMETER EXTRACTION FROM @config (Priority: explicit > config > default)
  # =========================================================================
  # Extract parameters from @config if not provided as arguments
  # This allows specifying once in TSENATAnalysis initialization
  
  # Condition column (for sample grouping)
  condition_col <- NULL
  if ("condition_col" %in% names(list(...))) {
    condition_col <- list(...)$condition_col
  } else if ("condition_col" %in% names(analysis@config)) {
    condition_col <- analysis@config$condition_col
  } else if ("condition" %in% colnames(colData(analysis@se))) {
    # Auto-detect common column name
    condition_col <- "condition"
  } else if ("sample_type" %in% colnames(colData(analysis@se))) {
    # Try another common name
    condition_col <- "sample_type"
  }
  
  # Method (statistical approach)
  if (is.null(method)) {
    if ("method" %in% names(analysis@config)) {
      method <- analysis@config$method
    } else {
      method <- "lmm"  # Default to LMM with AR(1) covariance
    }
  }
  
  # Paired design flag
  paired <- FALSE
  if ("paired" %in% names(list(...))) {
    paired <- list(...)$paired
  } else if ("paired" %in% names(analysis@config)) {
    paired <- analysis@config$paired
  }
  
  # Subject column (for paired/hierarchical designs)
  subject_col <- NULL
  if ("subject_col" %in% names(list(...))) {
    subject_col <- list(...)$subject_col
  } else if ("subject_col" %in% names(analysis@config)) {
    subject_col <- analysis@config$subject_col
  }
  
  # Validate that required parameters are available
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(analysis@se))
    if (length(cd_cols) > 0) {
      cat("[WARNING] condition_col not specified. Available columns: ",
          paste(cd_cols, collapse = ", "), "\n")
    } else {
      cat("[WARNING] condition_col not specified and colData is empty. ",
          "Will be determined by calculate_lm_interaction().\n")
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
      
      # Add _q= suffix to column names to distinguish q-values
      colnames(assay_matrix) <- paste0(colnames(assay_matrix), "_q=", q_val_str)
      
      assay_list[[key]] <- assay_matrix
      
      # Get colData
      coldata <- SummarizedExperiment::colData(se)
      # Update rownames to match new column names
      rownames(coldata) <- colnames(assay_matrix)
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
    stop(paste0("[calculate_lm_interaction_s4] Failed to combine diversity results: ",
                conditionMessage(e)), call. = FALSE)
  })
  
  # Build arguments for calculate_lm_interaction
  args <- list(se = diversity_combined)
  
  if (!is.null(method)) {
    args$method <- method
  }
  
  # Add extracted parameters from @config if not already provided in ...
  if (!is.null(condition_col) && !("condition_col" %in% names(list(...)))) {
    args$condition_col <- condition_col
  }
  
  if (paired && !("paired" %in% names(list(...)))) {
    args$paired <- paired
  }
  
  if (!is.null(subject_col) && !("subject_col" %in% names(list(...)))) {
    args$subject_col <- subject_col
  }

  # Request model_data for plotting compatibility (unless explicitly disabled)
  if (!("return_model_data" %in% names(args))) {
    args$return_model_data <- TRUE
  }
  
  # Merge with additional args (which may include condition_col, paired, multicorr, etc.)
  args <- c(args, list(...))

  # Run LM analysis
  result <- tryCatch({
    do.call(calculate_lm_interaction, args)
  }, error = function(e) {
    cat("  Message:", conditionMessage(e), "\n")
    cat("  Call:", paste(deparse(e$call), collapse="\n"), "\n")
    stop(paste0("Error in lm_interaction calculation:\n", conditionMessage(e)),
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
  # to be compatible with downstream functions like effect_sizes_divergence()
  if (!is.data.frame(lm_results_df)) {
    stop("[calculate_lm_interaction_s4] Result from calculate_lm_interaction must be a data.frame, ",
         "but got: ", class(lm_results_df),
         ". This indicates an unexpected change in the underlying function.",
         call. = FALSE)
  }
  
  required_cols <- c("gene", "adj_p_interaction")
  missing_cols <- setdiff(required_cols, colnames(lm_results_df))
  if (length(missing_cols) > 0) {
    stop("[calculate_lm_interaction_s4] Result missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Available columns: ", paste(colnames(lm_results_df), collapse = ", "),
         call. = FALSE)
  }
  
  if (nrow(lm_results_df) == 0) {
    warning("[calculate_lm_interaction_s4] Result is an empty data.frame. ",
            "This suggests filter out all genes (min_obs too high, or insufficient data)",
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

  analysis
}

# ============================================================================
# JACKKNIFE WRAPPER
# ============================================================================

#' Jackknife resampling with confidence intervals
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for jackknife. Default: 1.0.
#' @param ... Additional arguments passed to \code{\link{jackknife_tsallis_entropy}}.
#'
#' @return Modified TSENATAnalysis with jackknife results in @jackknife_results.
#'
#' @details
#' Requires diversity results to exist first. Will error if
#' \code{calculate_diversity_s4()} has not been run.
#'
#' @examples
#' \dontrun{
#'   analysis <- calculate_diversity_s4(analysis, q = 1.0)
#'   analysis <- jackknife_tsallis_entropy_s4(analysis, q = 1.0)
#'   ci <- jackKnife(analysis, q = 1.0)$confidence_intervals
#' }
#'
#' @export
jackknife_tsallis_entropy_s4 <- function(analysis, q = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check prerequisite: diversity must be calculated
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # Priority 1: Use explicit parameter
  # Priority 2: Use @config$q_values
  # Priority 3: Use @config$q_values_default or 1.0
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values
    } else {
      q <- 1.0
    }
  }

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

    # Get diversity result
    div_result <- analysis@diversity_results[[div_key]]

    # Run jackknife
    tryCatch({
      result <- jackknife_tsallis_entropy(
        se = analysis@se,
        q = q_val,
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
      stop(paste0("Jackknife error for q=", q_val, ":\n", e$message),
           call. = FALSE)
    })
  }

  analysis
}

# ============================================================================
# DIVERGENCE WRAPPER
# ============================================================================

#' Calculate divergence metrics and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value for divergence.
#'   If NULL, uses first q_value from @config$q_values if available, else defaults to 1.0.
#' @param ... Additional arguments passed to \code{\link{calculate_divergence}},
#'   including: control_group, paired, bootstrap, method, ci, etc.
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
#' \dontrun{
#'   # With @config setup:
#'   analysis <- TSENATAnalysis(se, 
#'     config = list(control_group = "normal", paired = TRUE))
#'   analysis <- calculate_diversity_s4(analysis, q = 1.0)
#'   analysis <- calculate_divergence_s4(analysis)  # Uses params from @config
#'   
#'   # Or explicit:
#'   analysis <- calculate_divergence_s4(analysis, 
#'     q = 1.0, control_group = "normal", paired = TRUE)
#'   div_metrics <- analysis@divergence_results
#' }
#'
#' @export
calculate_divergence_s4 <- function(analysis, q = NULL, ...) {
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

  # =========================================================================
  # PARAMETER EXTRACTION FROM @config (Priority: explicit > config > default)
  # =========================================================================
  
  # Q-value: Priority 1: explicit, 2: first q from @config, 3: default to 1.0
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values[1]  # Use first q-value
    } else {
      q <- 1.0
    }
  }
  
  # Control group for divergence comparison
  control_group <- NULL
  if ("control_group" %in% names(list(...))) {
    control_group <- list(...)$control_group
  } else if ("control_group" %in% names(analysis@config)) {
    control_group <- analysis@config$control_group
  }
  
  # Paired design flag
  paired <- FALSE
  if ("paired" %in% names(list(...))) {
    paired <- list(...)$paired
  } else if ("paired" %in% names(analysis@config)) {
    paired <- analysis@config$paired
  }
  
  # Statistical method
  method <- NULL
  if ("method" %in% names(list(...))) {
    method <- list(...)$method
  } else if ("method" %in% names(analysis@config)) {
    method <- analysis@config$method
  }
  
  # Bootstrap parameters
  bootstrap <- FALSE
  if ("bootstrap" %in% names(list(...))) {
    bootstrap <- list(...)$bootstrap
  } else if ("bootstrap" %in% names(analysis@config)) {
    bootstrap <- analysis@config$bootstrap
  }

  # Run divergence calculation with extracted parameters
  args <- list(
    se = analysis@se,
    q = q
  )
  
  # Add parameters from @config if not already in ...
  if (!is.null(control_group) && !("control_group" %in% names(list(...)))) {
    args$control_group <- control_group
  }
  
  if (paired && !("paired" %in% names(list(...)))) {
    args$paired <- paired
  }
  
  if (!is.null(method) && !("method" %in% names(list(...)))) {
    args$method <- method
  }
  
  if (bootstrap && !("bootstrap" %in% names(list(...)))) {
    args$bootstrap <- bootstrap
  }
  
  # Merge with additional args (which may override @config values)
  args <- c(args, list(...))
  
  # Run divergence calculation
  result <- tryCatch({
    do.call(calculate_divergence, args)
  }, error = function(e) {
    stop(paste0("Error in divergence calculation:\n", e$message),
         call. = FALSE)
  })

  # =========================================================================
  # VALIDATE AND STORE RESULTS
  # =========================================================================
  # Store results in list format (required by @divergence_results slot)
  if (is.null(result)) {
    warning("[calculate_divergence_s4] Result is NULL. Check calculate_divergence() output.",
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

  analysis
}

# ============================================================================
# DETECT Q GENE INTERACTIONS WRAPPER
# ============================================================================

#' Detect q-dependent gene interactions
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q_values \code{numeric}. Q-values to test across spectrum.
#' @param ... Additional arguments passed to \code{\link{detect_q_gene_interactions}}.
#'
#' @return Modified TSENATAnalysis with interaction results in @lm_results.
#'
#' @details
#' Analyzes how gene interactions change across q-value spectrum.
#'
#' @examples
#' \dontrun{
#'   analysis <- detect_q_gene_interactions_s4(
#'     analysis,
#'     q_values = seq(0.5, 2.0, by = 0.5)
#'   )
#' }
#'
#' @export
detect_q_gene_interactions_s4 <- function(analysis, q_values = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check prerequisites
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # =========================================================================
  # OPTIMIZATION (B: Lazy Conversion):
  # Use combined diversity result stored in metadata (if available)
  # Bypasses expensive per-q recombination logic
  # =========================================================================
  
  # Check if we have combined diversity result cached (from lazy evaluation)
  if (!is.null(analysis@metadata$diversity_combined) && 
      is.list(analysis@metadata$diversity_combined) &&
      !is.null(analysis@metadata$diversity_combined$combined_se)) {
    
    cat("[detect_q_gene_interactions_s4] Using cached combined diversity result (OPTIMIZATION)\n")
    
    # Use cached combined SE directly - it has correct structure and metadata
    se_multi_q <- analysis@metadata$diversity_combined$combined_se
    
    # Verify it's valid
    if (!is(se_multi_q, "SummarizedExperiment") || ncol(se_multi_q) == 0) {
      cat("[detect_q_gene_interactions_s4] Cached SE invalid, falling back to per-q recombination\n")
      se_multi_q <- NULL
    }
  } else {
    se_multi_q <- NULL
  }
  
  # Fallback: recombine per-q results if cache not available or invalid
  if (is.null(se_multi_q)) {
    
    # Extract q-values from diversity_results keys (format: "q_0.5", "q_1.0", etc.)
    q_keys <- names(analysis@diversity_results)
    q_values_extracted <- as.numeric(sub("^q_", "", q_keys))
    q_values_extracted <- sort(q_values_extracted)

    # Combine list of SEs (one per q-value) into a single SE for detection
    # Each SE has same genes but different q-value data
    combined_assay_list <- list()
    combined_coldata_list <- list()
    common_rownames <- NULL

    for (key in sort(q_keys)) {
      se <- analysis@diversity_results[[key]]
      
      # Extract q-value from key
      q_val <- as.numeric(sub("^q_", "", key))
      
      # Ensure SE format
      if (!is(se, "SummarizedExperiment")) {
        if (is.matrix(se) || is.data.frame(se)) {
          se <- SummarizedExperiment(assays = list(entropy = as.matrix(se)))
        } else {
          stop(paste0("Diversity result for ", key, " is not a SummarizedExperiment or matrix"),
               call. = FALSE)
        }
      }
      
      # Get assay data
      if (length(SummarizedExperiment::assays(se)) == 0) {
        stop(paste0("Diversity result for ", key, " has no assays"), call. = FALSE)
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
          stop(paste0("Diversity result for ", key, " has different number of genes"),
               call. = FALSE)
        }
      }
      rownames(assay_data) <- common_rownames
      
      # Make unique column names by appending q-value using old format "_q=" for compatibility
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

    # Combine all assays horizontally (cbind columns from different q-values)
    combined_assay <- do.call(cbind, combined_assay_list)
    
    # Combine colData - rownames should already be set to unique column names from the assay
    combined_coldata_df <- do.call(rbind, combined_coldata_list)
    
    # Ensure colnames of combined_assay match rownames of combined_coldata_df
    colnames(combined_assay) <- rownames(combined_coldata_df)
    
    # Get rowData from first diversity result (genes are same across all q-values)
    first_se <- analysis@diversity_results[[sort(q_keys)[1]]]
    
    # Ensure first_se is a SummarizedExperiment
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

  # Run q-interaction detection on combined SE
  result <- tryCatch({
    detect_q_gene_interactions(
      data = se_multi_q,
      ...
    )
  }, error = function(e) {
    stop(paste0("Error in q-interaction detection:\n", e$message),
         call. = FALSE)
  })

  # Store in lm_results under "q_interactions" key
  if (is.list(analysis@lm_results)) {
    analysis@lm_results$q_interactions <- result
  } else {
    analysis@lm_results <- list(q_interactions = result)
  }

  # Track metadata
  if (!is.null(q_vals_for_tracking)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("detect_q_gene_interactions[q=", paste(q_vals_for_tracking, collapse = ","), "]")
    )
  }

  analysis
}

# ============================================================================
# CALCULATE DIFFERENCE WRAPPER
# ============================================================================

#' Calculate Difference Between Control and Treatment Groups (S4 Wrapper)
#'
#' S4 wrapper for \code{\link{calculate_difference}} that operates on TSENATAnalysis objects.
#' Uses diversity results from \code{@diversity_results} slot (from \code{calculate_diversity_s4()})
#' and stores results in the \code{lm_results} slot.
#'
#' @param analysis A \code{TSENATAnalysis} object with diversity results in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value to use. If NULL, uses first diversity result or q=1.0.
#' @param control Character string specifying the control group identifier. If \code{NULL},
#'   attempts to retrieve from \code{analysis@config$control}.
#' @param ... Additional arguments passed to \code{\link{calculate_difference}}.
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
#' @export
#' @seealso \code{\link{calculate_difference}} for the underlying implementation.
#'
#' @examples
#' \dontrun{
#'   # First compute diversity
#'   analysis <- calculate_diversity_s4(analysis, q = 1.0)
#'   
#'   # Then compute differences
#'   analysis <- calculate_difference_s4(analysis, control = "Normal")
#' }
calculate_difference_s4 <- function(analysis, control = NULL, q = NULL, ...) {
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
  # Priority: explicit q > first q in diversity_results
  if (is.null(q)) {
    # Use first diversity result
    div_keys <- names(analysis@diversity_results)
    if (length(div_keys) == 0) {
      stop("No diversity results found in @diversity_results", call. = FALSE)
    }
    diversity_se <- analysis@diversity_results[[div_keys[1]]]
    q_used <- sub("^q_", "", div_keys[1])
  } else {
    # Find diversity result for specified q
    q_key <- paste0("q_", formatC(q, format = "f", digits = 3))
    if (!(q_key %in% names(analysis@diversity_results))) {
      stop("Diversity not calculated for q=", q, 
           ". Available: ", paste(names(analysis@diversity_results), collapse = ", "),
           call. = FALSE)
    }
    diversity_se <- analysis@diversity_results[[q_key]]
    q_used <- q
  }

  # Run difference calculation on diversity results (not raw input @se)
  result <- tryCatch({
    # Determine samples column from colData
    # Use analysis@se's colData if available, otherwise use diversity_se's colData
    if ("group" %in% colnames(SummarizedExperiment::colData(analysis@se))) {
      samples_col <- "group"
    } else if ("sample_type" %in% colnames(SummarizedExperiment::colData(analysis@se))) {
      samples_col <- "sample_type"
    } else if ("condition" %in% colnames(SummarizedExperiment::colData(analysis@se))) {
      samples_col <- "condition"
    } else {
      stop("No sample grouping column found in colData. ",
           "Expected: 'group', 'sample_type', or 'condition'",
           call. = FALSE)
    }

    calculate_difference(
      x = diversity_se,
      samples = samples_col,
      control = control,
      ...
    )
  }, error = function(e) {
    stop(paste0("Error in difference calculation:\n", e$message),
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
#' This wrapper calls \code{test_rankbased_assumptions()} on diversity data
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
#' \dontrun{
#'   analysis <- TSENATAnalysis(se = se_data, config = list())
#'   analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5))
#'   analysis <- test_rankbased_assumptions_s4(analysis, q = 1.0)
#'   str(analysis@metadata$rankbased_assumptions)
#' }
#'
#' @export
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
      
      # Combine all matrices column-wise (genes × all samples across q-values)
      diversity_data <- do.call(cbind, entropy_list)
      # Keep natural column names from cbind to preserve structure
      q_used <- "all"
    }
    
    # If q is NULL and only one result, use it
    if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) == 1) {
      div_se <- analysis@diversity_results[[1]]
      diversity_data <- assay(div_se, "diversity")
      q_used <- extract_q_from_key(names(analysis@diversity_results)[1])
    }
    
    # Fallback: use first diversity result
    if (is.null(diversity_data) && length(analysis@diversity_results) > 0) {
      div_se <- analysis@diversity_results[[1]]
      diversity_data <- assay(div_se, "diversity")
      if (is.null(q_used) || is.na(q_used)) {
        q_used <- extract_q_from_key(names(analysis@diversity_results)[1])
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
      test_rankbased_assumptions(
        data = diversity_data,
        checks = checks,
        alpha = alpha
      )
    }, error = function(e) {
      stop(paste0("Error in rankbased assumptions test:\n", e$message),
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
extract_q_from_key <- function(key) {
  # Extract numeric part from "q_X.X" format
  as.numeric(sub("^q_", "", key))
}
