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
#'   If NULL, uses q_values from \code{analysis@config$q_values} if available, else defaults to seq(0.01, 2, by = 0.05).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for tables), .rds (for S4 objects). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function,
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
#'   \item{q}{Priority 1 (explicit) > Priority 2 (\code{@config$q_values}) > Priority 3 (default: seq(0.01, 2, by = 0.05))\cr
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
#' # Create minimal test data
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' # Calculate diversity for single q-value
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # View results
#' head(SummarizedExperiment::assay(analysis@diversity_results$q_1.0))
#'
#' @export
#' @importFrom utils write.table
calculate_diversity_s4 <- function(analysis, q = NULL, output_file = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  if (nrow(analysis@se) == 0) {
    stop("SummarizedExperiment in @se is empty", call. = FALSE)
  }

  # Priority 1: Use explicit parameter if provided
  # Priority 2: Use @config$q_values if set
  # Priority 3: Default to seq(0.01, 2, by = 0.05)
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values
    } else {
      q <- seq(0.01, 2, by = 0.05)
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
  } else if ("metadata" %in% names(analysis@config)) {
    analysis@config$metadata
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
          numeric_cols <- vapply(result_df, is.numeric, FUN.VALUE = logical(1))
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
        numeric_cols <- vapply(result_subset, is.numeric, FUN.VALUE = logical(1))
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
          stop("[calculate_diversity_s4] Converted SE for q=", q_val, 
                      " has no assays. Check diversity result structure.",
               call. = FALSE)
        }
        test_assay <- tryCatch({
          SummarizedExperiment::assay(result_se, 1)
        }, error = function(e) {
          stop("[calculate_diversity_s4] Cannot access assay in SE for q=", q_val,
                      ": ", conditionMessage(e), call. = FALSE)
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
        stop("[calculate_diversity_s4] Bootstrap CI computation failed for q=", q_val, ": ", error_msg,
             call. = FALSE)
      } else {
        stop("[calculate_diversity_s4] Failed to compute diversity for q=", q_val, ":\n", error_msg,
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

  # Save if output_file provided
  if (!is.null(output_file)) {
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write as text table (convert to data.frame representation)
      # Try to extract writable data from SummarizedExperiment or other formats
      tryCatch({
        if (length(analysis@diversity_results) > 0 && is(analysis@diversity_results[[1]], "SummarizedExperiment")) {
          # Extract assay data from first SummarizedExperiment
          write_data <- as.data.frame(assay(analysis@diversity_results[[1]]))
        } else {
          write_data <- as.data.frame(analysis@diversity_results)
        }
        write.table(write_data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
      }, error = function(e) {
        warning("[calculate_diversity_s4] Could not write diversity results to file: ", 
                conditionMessage(e), call. = FALSE)
      })
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }
  }

  analysis
}

#' Calculate LM interactions and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param fdr_threshold \code{numeric}. FDR cutoff for significance.
#'   Default: 0.05.
#' @param formula \code{formula} or NULL. Reserved for future use.
#' @param method \code{character}. Statistical method (e.g., "lmm", "gam", "gee").
#'   If NULL, uses method from @config$method or defaults to "lmm".
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base LM function,
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
#' # Create test data with sufficient structure for LM analysis
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' # Generate count data with higher lambda for better signal
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(counts))
#' # Create analysis and calculate diversity first (prerequisite)
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # Then calculate LM interactions
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' if (!is.null(analysis@lm_results$lm_interaction) && 
#'     nrow(analysis@lm_results$lm_interaction) > 0) {
#'   head(analysis@lm_results$lm_interaction)
#' }
#'
#' @export
#' @importFrom utils write.table
calculate_lm_interaction_s4 <- function(analysis, fdr_threshold = NULL, 
                                       formula = NULL, method = NULL, output_file = NULL, ...) {
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
  
  # Number of threads for parallel computation
  nthreads <- NULL
  if ("nthreads" %in% names(list(...))) {
    nthreads <- list(...)$nthreads
  } else if ("nthreads" %in% names(analysis@config)) {
    nthreads <- analysis@config$nthreads
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
      message("condition_col not specified. Available columns: ",
          paste(cd_cols, collapse = ", "))
    } else {
      message("condition_col not specified and colData is empty. ",
          "Will be determined by calculate_lm_interaction().")
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
  
  # Add extracted parameters from @config if not already provided in ...
  if (!is.null(condition_col) && !("condition_col" %in% names(list(...)))) {
    args$condition_col <- condition_col
  }
  
  # Always pass paired if from config (even if FALSE)
  if (!("paired" %in% names(list(...)))) {
    args$paired <- paired
  }
  
  if (!is.null(subject_col) && !("subject_col" %in% names(list(...)))) {
    args$subject_col <- subject_col
  }

  if (!is.null(nthreads) && !("nthreads" %in% names(list(...)))) {
    args$nthreads <- nthreads
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
    message("  ", conditionMessage(e))
    message("  Call: ", paste(deparse(e$call), collapse="\n"))
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
  # to be compatible with downstream functions like effect_sizes_divergence()
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
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write LM results as text table
      tryCatch({
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
          write.table(lm_data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
        }
        # If no valid data.frame found, silently skip writing (don't try to force conversion)
      }, error = function(e) {
        # Silently skip if there's an error writing the output file
        # The analysis results are still valid, just not written to disk
      })
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
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
#' @param print_results \code{logical}. Print jackknife results summary. Default: FALSE.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with jackknife results in @jackknife_results.
#'
#' @details
#' Requires diversity results to exist first. Will error if
#' \code{calculate_diversity_s4()} has not been run.
#'
#' @examples
#' # Create minimal test data for demonstration
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' # Calculate diversity first (required)
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # Run jackknife estimation
#' analysis <- jackknife_tsallis_entropy_s4(analysis, q = 1.0, verbose = FALSE)
#' # Check jackknife results
#' names(analysis@jackknife_results)
#'
#' @export
#' @importFrom utils write.table
jackknife_tsallis_entropy_s4 <- function(analysis, q = NULL, print_results = FALSE, output_file = NULL, ...) {
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
      result <- jackknife_tsallis_entropy(
        x = div_matrix,
        q = q_val,
        print_results = print_results,
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
    saveRDS(analysis, file = output_file)

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
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base divergence function,
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
#' # Create and run divergence analysis
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' # First calculate diversity
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # Then calculate divergence
#' analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
#' # Check divergence results
#' head(analysis@divergence_results)
#'
#' @export
#' @importFrom utils write.table
calculate_divergence_s4 <- function(analysis, q = NULL, verbose = TRUE, output_file = NULL, ...) {
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
  
  # Q-value: Priority 1: explicit, 2: full q_values from @config, 3: default to 1.0
  # NOTE: Supports multi-q divergence spectrum analysis (unlike single-q defaults)
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values  # Use all q-values for spectrum analysis
    } else {
      q <- 1.0
    }
  }
  
  # Replace q=0 with q=0.01 for practical approximation
  # (q=0 divergence always returns 0, which is mathematically correct but uninformative)
  if (is.vector(q)) {
    q[q == 0] <- 0.01
  } else if (is.numeric(q) && q == 0) {
    q <- 0.01
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
    q = q,
    verbose = verbose
  )
  
  # Add parameters from @config if not already in ...
  if (!is.null(control_group) && !("control_group" %in% names(list(...)))) {
    args$control_group <- control_group
  }
  
  # Always pass paired if from config (even if FALSE)
  if (!("paired" %in% names(list(...)))) {
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
    stop("Divergence calculation failed:\n", e$message,
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
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with interaction results in @lm_results.
#'
#' @details
#' Analyzes how gene interactions change across q-value spectrum.
#'
#' **Parameter resolution priority** (explicit > @config > extract from results):
#' \itemize{
#'   \item \code{q}: Uses explicit arg, else \code{@config$q_values},
#'     else extracts from diversity_results keys
#' }
#'
#' @examples
#' # Setup analysis with diversity results
#' library(SummarizedExperiment)
#' set.seed(42)
#' # Create test data with sufficient structure for q-gene interactions
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' # Detect q-gene interactions across q-values
#' # (requires both diversity and lm_interaction results)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' analysis <- detect_q_gene_interactions_s4(analysis, q = c(0.5, 1.0, 1.5))
#' # Interaction results stored in lm_results
#' if (!is.null(analysis@lm_results$lm_interaction)) {
#'   head(analysis@lm_results$lm_interaction)
#' }
#'
#' @export
#' @importFrom utils write.table
detect_q_gene_interactions_s4 <- function(analysis, q = NULL, output_file = NULL, ...) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check prerequisites
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # =========================================================================
  # PARAMETER EXTRACTION FROM @config (Priority: explicit > @config > extract from results)
  # =========================================================================
  # If q not provided, check @config
  if (is.null(q)) {
    if ("q_values" %in% names(analysis@config)) {
      q <- analysis@config$q_values
    }
    # Otherwise, will extract from diversity_results keys below
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
    
    
    # Use cached combined SE directly - it has correct structure and metadata
    se_multi_q <- analysis@metadata$diversity_combined$combined_se
    
    # Verify it's valid
    if (!is(se_multi_q, "SummarizedExperiment") || ncol(se_multi_q) == 0) {
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
    stop("q-interaction detection failed:\n", e$message,
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

  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write q-interactions results as text table
      write.table(as.data.frame(result), file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }

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
#'   \item \code{condition_col} (sample grouping): Uses \code{@config$condition_col},
#'     else auto-detects from colData columns: "group", "sample_type", "condition"
#' }
#'
#' @importFrom utils write.table
#' @seealso
#' \code{\link{calculate_diversity_s4}} for computing diversity.
#'
#' @examples
#' # First compute diversity with adequate sample size
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("Control", "Treatment"), each = n_samples_per_group),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # Then compute differences
#' analysis <- calculate_difference_s4(analysis, control = "Control", 
#'   verbose = FALSE)
#' if (!is.null(analysis@lm_results$difference) && 
#'     nrow(analysis@lm_results$difference) > 0) {
#'   head(analysis@lm_results$difference)
#' }
#'
#' @export
calculate_difference_s4 <- function(analysis, control = NULL, q = NULL, output_file = NULL, ...) {
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

  # Ensure colData from original SE is preserved in diversity results
  if (is(diversity_se, "SummarizedExperiment")) {
    original_coldata <- SummarizedExperiment::colData(analysis@se)
    if (!is.null(original_coldata) && nrow(original_coldata) == ncol(diversity_se)) {
      SummarizedExperiment::colData(diversity_se) <- original_coldata
    }
  }

  # Run difference calculation on diversity results (not raw input @se)
  result <- tryCatch({
    # Determine samples column from colData
    # Priority 1: Check @config$condition_col
    # Priority 2: Auto-detect from common column names
    samples_col <- NULL
    
    if ("condition_col" %in% names(analysis@config)) {
      samples_col <- analysis@config$condition_col
    }
    
    # Fallback to auto-detection if not in config
    if (is.null(samples_col)) {
      cd_cols <- colnames(SummarizedExperiment::colData(analysis@se))
      if ("group" %in% cd_cols) {
        samples_col <- "group"
      } else if ("sample_type" %in% cd_cols) {
        samples_col <- "sample_type"
      } else if ("condition" %in% cd_cols) {
        samples_col <- "condition"
      } else {
        stop(
          "Could not determine sample grouping column:\n",
          "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
          "SOLUTION: Set @config$condition_col with the correct column name\n",
          "  Example: analysis@config$condition_col <- 'sample_type'\n",
          call. = FALSE)
      }
    }

    # Call calculate_difference directly (internal TSENAT function)
    calculate_difference(
      x = diversity_se,
      condition_col = samples_col,
      control = control,
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

  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
      # Write difference results as text table
      diff_data <- if (!is.null(analysis@lm_results$difference$results)) {
        analysis@lm_results$difference$results
      } else {
        as.data.frame(analysis@lm_results$difference)
      }
      write.table(diff_data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
    } else {
      # Default to RDS for S4 object
      saveRDS(analysis, file = output_file)
    }

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
#' # Create minimal test data
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' # Test rank-based assumptions at q=1.0
#' analysis <- test_rankbased_assumptions_s4(analysis, q = 1.0)
#' # Check results
#' names(analysis@metadata$rankbased_assumptions)
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
      
      # Combine all matrices column-wise (genes x all samples across q-values)
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
extract_q_from_key <- function(key) {
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
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return
#' Invisibly returns a cowplot grid object containing both volcano and MA plots
#' combined side-by-side. If the plot cannot be created, returns NULL invisibly.
#'
#' @details
#' This wrapper extracts the difference results data frame from
#' \code{analysis@lm_results$difference} and passes it to the base
#' \code{plot_volcano_ma_grid()} function.
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
#' # After calculating diversity and differences with adequate sample size
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("Normal", "Tumor"), each = n_samples_per_group),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' analysis <- calculate_difference_s4(analysis, control = "Normal", 
#'   verbose = FALSE)
#' # Create volcano and MA plot grid (if ggplot2 available)
#' if (!is.null(analysis@plots$volcano_ma_grid)) {
#'   print(analysis@plots$volcano_ma_grid)
#' }
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
    ...) {

  # Auto-detect verbose from config if not explicitly provided
  if (isTRUE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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

  # Validate padj_col exists
  # Handle both "padj" and "adjusted_p_values" column names
  actual_padj_col <- padj_col
  if (!(padj_col %in% colnames(diff_df))) {
    # Try alternative column name
    if ("adjusted_p_values" %in% colnames(diff_df) && padj_col == "padj") {
      actual_padj_col <- "adjusted_p_values"
    } else if ("pvalue" %in% colnames(diff_df) && padj_col == "padj") {
      # Fallback to raw p-values if adjusted not available
      actual_padj_col <- "pvalue"
    } else {
      stop("Column '", padj_col, "' not found in difference results. ",
           "Available columns: ", paste(colnames(diff_df), collapse = ", "),
           call. = FALSE)
    }
  }

  # Determine x_col if not provided
  # Prefer mean_difference, then log2_fold_change, then mean
  if (is.null(x_col)) {
    if ("mean_difference" %in% colnames(diff_df)) {
      x_col <- "mean_difference"
    } else if ("log2_fold_change" %in% colnames(diff_df)) {
      x_col <- "log2_fold_change"
    } else {
      warning("Could not auto-detect x_col. Available numeric columns: ",
              paste(colnames(diff_df)[vapply(diff_df, is.numeric, FUN.VALUE = logical(1))], collapse = ", "),
              call. = FALSE)
    }
  }

  plot_obj <- tryCatch({
    plot_volcano_ma_grid(
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

  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = plot_obj, device = "pdf")
    } else if (grepl("\\.png$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = plot_obj, device = "png")
    } else if (grepl("\\.jpg$|\\.jpeg$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = plot_obj, device = "jpeg")
    } else {
      # Default to PDF
      ggplot2::ggsave(paste0(output_file, ".pdf"), plot = plot_obj, device = "pdf")
    }
    if (verbose) {
      message("[plot_volcano_ma_grid_s4] Plot saved to ", output_file)
    }

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
#'   Default: "q_interactions" (results from \code{detect_q_gene_interactions_s4})
#' @param friedman_method \code{character}. Key for Friedman/rank-based results in \code{@lm_results}.
#'   Default: "rankbased" (results from \code{test_rankbased_assumptions_s4})
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
#' # Generic example showing compute_method_concordance_s4 structure
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Note: compute_method_concordance_s4 requires results from both
#' # detect_q_gene_interactions_s4 and test_rankbased_assumptions_s4
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
    verbose = FALSE,
    output_file = NULL) {
  
  # ===================================================================
  # VALIDATION
  # ===================================================================
  
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  if (is.null(analysis@lm_results)) {
    stop("No LM results found in analysis@lm_results. Run detect_q_gene_interactions_s4() first.",
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
  gam_results <- analysis@lm_results[[gam_method]]
  friedman_results <- analysis@lm_results[[friedman_method]]
  
  # Validate they're data frames
  if (!is.data.frame(gam_results)) {
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
    compute_method_concordance(gam_results, friedman_results)
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
#' # After computing divergence via S4 wrapper
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
#' # Global divergence curve (all genes aggregated) - default mode
#' p_global <- plot_divergence_spectrum_s4(analysis)
#' print(p_global)
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
    width = 10,
    height = 6,
    verbose = TRUE,
    ...) {

  # Auto-detect verbose from config if not explicitly provided
  if (isTRUE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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

  # Create the plot using base function
  p <- tryCatch({
    plot_divergence_spectrum(
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

  # Save to file if requested
  output_file_path <- NULL
  if (!is.null(output_file)) {
    tryCatch({
      ggplot2::ggsave(
        filename = output_file,
        plot = p,
        width = width,
        height = height,
        dpi = 300
      )
      output_file_path <- output_file
      if (verbose) {
        message("Saved divergence spectrum plot to: ", output_file)
      }
    }, error = function(e) {
      if (verbose) {
        message("Could not save plot to file: ", e$message)
      }
    })
  }

  # Return file path if saved, otherwise return plot
  if (!is.null(output_file_path)) {
    invisible(output_file_path)
  } else {
    invisible(p)
  }
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
#' # After computing concordance:
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' # Note: compute_method_concordance_s4 requires additional setup
#' # For demo, we show that plot_method_concordance_s4 needs pre-computed concordance
#'
#' @aliases plot_method_concordance_s4
#' @export
setGeneric("plot_method_concordance_s4", function(analysis, verbose = FALSE) {
  standardGeneric("plot_method_concordance_s4")
})

#' @rdname plot_method_concordance_s4
setMethod("plot_method_concordance_s4", "TSENATAnalysis", function(analysis, verbose = FALSE) {
  
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
  plot_obj <- plot_method_concordance(comparison_df)
  
  if (verbose) {
    message("[plot_method_concordance_s4] Plot generated successfully")
  }
  
  return(plot_obj)
})

#' Compute Effect Sizes from Divergence Results (S4 Wrapper)
#'
#' S4 wrapper for \code{effect_sizes_divergence()} that extracts divergence and LM
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
#' @return Modified TSENATAnalysis with effect size results stored in
#'   \code{@metadata$effect_sizes_divergence}. Returns the analysis object visibly
#'   to support piping and method chaining.
#'
#' @details
#' **Workflow steps:**
#' \describe{
#'   \item{Extracting}{Divergence SE from \code{@divergence_results} and LM results
#'     from \code{@lm_results$lm_interaction}}
#'   \item{Computing}{Effect sizes using standard \code{effect_sizes_divergence()} function}
#'   \item{Storing}{Results as list with \code{interaction_results} (data.frame) and
#'     \code{validation_stats}}
#'   \item{Tracking}{Function call in \code{@metadata$function_calls}}
#' }
#'
#' **Parameter resolution priority** (explicit > @config > default):
#' \itemize{
#'   \item \code{significance_threshold}: Uses explicit arg, else \code{@config$significance_threshold},
#'     else 0.05
#'   \item \code{enrich_per_q_pattern}: Uses explicit arg, else \code{@config$enrich_per_q_pattern},
#'     else TRUE
#'   \item \code{verbose}: Uses explicit arg, else \code{@config$verbose}, else TRUE
#' }
#'
#' Results are accessed via: \code{analysis@metadata$effect_sizes_divergence}
#'
#' @examples
#' # Compute effect sizes after divergence and LM analysis with adequate data
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Note: effect_sizes_divergence_s4 requires specialized data structures
#' # Check that LM results were computed
#' if (!is.null(analysis@lm_results$lm_interaction)) {
#'   nrow(analysis@lm_results$lm_interaction)
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

  # Auto-detect verbose from config if not explicitly provided
  if (isTRUE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Check for required results
  if (length(analysis@divergence_results) == 0) {
    stop("Divergence results required. Run calculate_divergence_s4() first.",
         call. = FALSE)
  }

  if (is.null(analysis@lm_results) || length(analysis@lm_results) == 0) {
    stop("LM results required. Run calculate_lm_interaction_s4() first.",
         call. = FALSE)
  }

  # =========================================================================
  # PARAMETER EXTRACTION FROM @config (Priority: explicit > @config > default)
  # =========================================================================
  # Extract parameters from @config if not provided as arguments
  dots <- list(...)
  
  # significance_threshold parameter
  if (is.null(significance_threshold) || identical(significance_threshold, 0.05)) {
    if ("significance_threshold" %in% names(dots)) {
      significance_threshold <- dots$significance_threshold
    } else if ("significance_threshold" %in% names(analysis@config)) {
      significance_threshold <- analysis@config$significance_threshold
    }
  }
  
  # enrich_per_q_pattern parameter
  if (isTRUE(enrich_per_q_pattern)) {
    if ("enrich_per_q_pattern" %in% names(dots)) {
      enrich_per_q_pattern <- dots$enrich_per_q_pattern
    } else if ("enrich_per_q_pattern" %in% names(analysis@config)) {
      enrich_per_q_pattern <- analysis@config$enrich_per_q_pattern
    }
  }
  
  # verbose parameter
  if (isTRUE(verbose)) {
    if ("verbose" %in% names(dots)) {
      verbose <- dots$verbose
    } else if ("verbose" %in% names(analysis@config)) {
      verbose <- analysis@config$verbose
    }
  }

  # =========================================================================
  # EXTRACT RESULTS FROM ANALYSIS OBJECT
  # =========================================================================
  
  # Extract divergence SE
  divergence_se <- if ("divergence_se" %in% names(analysis@divergence_results)) {
    analysis@divergence_results$divergence_se
  } else if (is(analysis@divergence_results, "SummarizedExperiment")) {
    analysis@divergence_results
  } else if (is.list(analysis@divergence_results) && length(analysis@divergence_results) > 0) {
    # Fallback: check if first element is SE
    analysis@divergence_results[[1]]
  } else {
    NULL
  }

  if (is.null(divergence_se) || !is(divergence_se, "SummarizedExperiment")) {
    stop("Could not extract divergence SummarizedExperiment from analysis@divergence_results",
         call. = FALSE)
  }

  # Extract LM results
  lm_res <- if ("lm_interaction" %in% names(analysis@lm_results)) {
    analysis@lm_results$lm_interaction
  } else if (is.data.frame(analysis@lm_results)) {
    analysis@lm_results
  } else if (is.list(analysis@lm_results) && length(analysis@lm_results) > 0) {
    analysis@lm_results[[1]]
  } else {
    NULL
  }

  if (is.null(lm_res) || !is.data.frame(lm_res)) {
    stop("Could not extract LM results data.frame from analysis@lm_results",
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
    tx2gene <- metadata(analysis@se)$tx2gene
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
    effect_sizes_divergence(
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
  if (is.null(analysis@metadata)) {
    analysis@metadata <- list()
  }

  analysis@metadata$effect_sizes_divergence <- result

  # Track function call
  analysis@metadata$function_calls <- c(
    analysis@metadata$function_calls,
    paste0("effect_sizes_divergence[threshold=", significance_threshold, "]")
  )

  if (verbose) {
    message("[effect_sizes_divergence_s4] Results stored in @metadata$effect_sizes_divergence")
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
#' S4 wrapper for \code{plot_top_transcripts()} that extracts data directly from
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
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plotting (default: FALSE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Supported formats: .pdf, .png, .jpg. Default: NULL (no file output).
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return A file path (character) to the saved plot PNG file, invisibly.
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
#' # Plot top transcripts for the most significant gene
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Plot top transcripts (if ggplot2 available)
#' plot_file <- plot_top_transcripts_s4(analysis, top_n = 3)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Plot top transcripts (if ggplot2 available)
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
    top_n = 3,
    output_file = NULL,
    metric = c("median", "mean", "variance", "iqr"),
    verbose = FALSE,
    ...) {

  # Auto-detect verbose from config if not explicitly provided
  if (isFALSE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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

    # Try Priority 1: @config$condition_col
    if ("condition_col" %in% names(analysis@config)) {
      candidate <- analysis@config$condition_col
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }

    # Try Priority 2: @config$sample_type or direct colData
    if (is.null(condition_col) && "sample_type" %in% cd_cols) {
      condition_col <- "sample_type"
    }

    # Try Priority 3: @config$condition
    if (is.null(condition_col) && "condition" %in% cd_cols) {
      condition_col <- "condition"
    }

    # Fallback: use first column
    if (is.null(condition_col) && length(cd_cols) > 0) {
      condition_col <- cd_cols[1]
    }

    if (is.null(condition_col)) {
      stop("[plot_top_transcripts_s4] Cannot auto-detect condition_col. ",
           "Provide explicitly.", call. = FALSE)
    }

    if (verbose) {
      message("[plot_top_transcripts_s4] Auto-detected condition_col = ", condition_col)
    }
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
        # Find p-value column (handle various naming conventions)
        p_col <- if ("p_interaction" %in% colnames(lm_results_df)) {
          "p_interaction"
        } else if ("pvalue" %in% colnames(lm_results_df)) {
          "pvalue"
        } else if ("p.value" %in% colnames(lm_results_df)) {
          "p.value"
        } else if ("padj" %in% colnames(lm_results_df)) {
          "padj"
        } else if ("p_value" %in% colnames(lm_results_df)) {
          "p_value"
        } else {
          NA
        }

        # Find gene column (handle various naming conventions)
        gene_col <- if ("gene" %in% colnames(lm_results_df)) {
          "gene"
        } else if ("gene_name" %in% colnames(lm_results_df)) {
          "gene_name"
        } else if ("gene_id" %in% colnames(lm_results_df)) {
          "gene_id"
        } else {
          colnames(lm_results_df)[1]  # Fallback to first column
        }

        if (!is.na(p_col) && p_col %in% colnames(lm_results_df)) {
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
      message("[plot_top_transcripts_s4] Calling plot_top_transcripts() for genes: ", 
          paste(gene, collapse = ", "))
    } else {
      message("[plot_top_transcripts_s4] Calling plot_top_transcripts() for gene: ", gene)
    }
  }

  plot_file <- tryCatch({
    plot_top_transcripts(
      se = se,
      gene = gene,
      condition_col = condition_col,
      res = lm_results_df,
      top_n = top_n,
      output_file = output_file,
      metric = metric[1],  # Use first metric if multiple provided
      ...
    )
  }, error = function(e) {
    stop("[plot_top_transcripts_s4]", conditionMessage(e), call. = FALSE)
  })

  if (verbose) {
      message("[plot_top_transcripts_s4] Plot saved to: ", plot_file)
  }

  # Return file path invisibly
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
#' and passes them to the base \code{plot_divergence_distribution()} function.
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
#' # After computing effect sizes via S4 wrapper
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Compute effect sizes (required for plot_divergence_distribution_s4)
#' analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE)
#' # Generate and display the plot (if ggplot2 available)
#' analysis <- plot_divergence_distribution_s4(analysis, verbose = FALSE)
#'
#' @seealso
#' \code{\link{effect_sizes_divergence_s4}} for computing effect sizes.
#'
#' @export
plot_divergence_distribution_s4 <- function(
    analysis,
    threshold = 0.1,
    output_file = NULL,
    width = 10,
    height = 6,
    verbose = TRUE,
    ...) {

  # Auto-detect verbose from config if not explicitly provided
  if (isTRUE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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
    plot_divergence_distribution(
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
  output_file_path <- NULL
  if (!is.null(output_file)) {
    tryCatch({
      ggplot2::ggsave(
        filename = output_file,
        plot = p,
        width = width,
        height = height,
        dpi = 300
      )
      output_file_path <- output_file
      if (verbose) {
        message("Saved divergence distribution plot to: ", output_file)
      }
    }, error = function(e) {
      if (verbose) {
        message("Could not save plot to file: ", e$message)
      }
    })
  }

  # Return file path if saved, otherwise return plot
  if (!is.null(output_file_path)) {
    invisible(output_file_path)
  } else {
    invisible(p)
  }
}


#' Prepare Gene Switching Tables from TSENATAnalysis Object
#'
#' S4 wrapper for \code{prepare_gene_switching_tables()} that extracts results
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
#' # After running LM interaction analysis
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Prepare tables using S4 wrapper
#' tables <- prepare_gene_switching_tables_s4(analysis)
#' # Access individual components
#' head(tables$summary_table)
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
  
  # Auto-detect verbose from config if not explicitly provided
  if (isFALSE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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
  if (verbose) message("Calling prepare_gene_switching_tables()...")
  
  result <- prepare_gene_switching_tables(
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
#' S4 wrapper for \code{plot_multiq_delta_influence_heatmaps()} that extracts results
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
#' # After running LM interaction analysis
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Generate heatmaps using S4 wrapper (requires jackknife switching setup)
#' heatmap_file <- plot_multiq_delta_influence_heatmaps_s4(analysis, n_genes = 4)
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
  
  # Auto-detect verbose from config if not explicitly provided
  if (isFALSE(verbose)) {
    if ("verbose" %in% names(analysis@config)) {
      config_verbose <- analysis@config$verbose
      if (is.logical(config_verbose) && length(config_verbose) == 1) {
        verbose <- config_verbose
      }
    }
  }

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
  
  if (verbose) message("Calling plot_multiq_delta_influence_heatmaps()...")
  
  # Call base function with extracted parameters
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = switching_results,
    n_genes = n_genes,
    lm_results = lm_results,
    ...
  )
  
  if (verbose) {
    message("[OK] Heatmap plot generated successfully")
    message("  Saved to: ", heatmap_file)
  }

  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", tolower(output_file))) {
      # If heatmap_file is a PNG from base function, convert by re-saving as PDF
      # For now, just copy/save the plot object if available
      if (file.exists(heatmap_file)) {
        file.copy(heatmap_file, output_file, overwrite = TRUE)
      }
    } else if (grepl("\\.png$", tolower(output_file))) {
      # heatmap_file should already be PNG from base function
      if (file.exists(heatmap_file)) {
        file.copy(heatmap_file, output_file, overwrite = TRUE)
      }
    } else if (grepl("\\.jpg$|\\.jpeg$", tolower(output_file))) {
      # Convert PNG to JPEG if needed
      if (file.exists(heatmap_file)) {
        file.copy(heatmap_file, output_file, overwrite = TRUE)
      }
    } else {
      # Default: copy base output to requested file
      if (file.exists(heatmap_file)) {
        file.copy(heatmap_file, output_file, overwrite = TRUE)
      }
    }
    if (verbose && file.exists(output_file)) {
      message("[plot_multiq_delta_influence_heatmaps_s4] Plot saved to ", output_file)
    }
  }

  # Return result visibly (consistent with other S4 wrappers)
  heatmap_file
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
#' 4. Calls \code{plot_lm_interaction_gam()} with extracted parameters
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
#' # Create TSENATAnalysis with diversity and LM results
#' library(SummarizedExperiment)
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 5
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                          nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
#'                            nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   pair = rep(1:n_samples_per_group, 2),
#'   row.names = colnames(se))
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
#' suppressWarnings(
#'   analysis <- calculate_lm_interaction_s4(analysis,
#'     condition_col = "condition", verbose = FALSE)
#' )
#' # Plot GAM curves for top genes
#' if (!is.null(analysis@lm_results$lm_interaction)) {
#'   plot <- plot_lm_interaction_gam_s4(analysis, n_top = 3,
#'     condition_col = "condition")
#' }
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
  ...
) {
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
    
    # Try Priority 1: @config$condition_col (validate it exists)
    if ("condition_col" %in% names(analysis@config)) {
      candidate <- analysis@config$condition_col
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }
    
    # Try Priority 2: @config$condition (validate it exists)
    if (is.null(condition_col) && "condition" %in% names(analysis@config)) {
      candidate <- analysis@config$condition
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }
    
    # Try Priority 3: Check for common column names in actual colData
    if (is.null(condition_col)) {
      if ("condition" %in% cd_cols) {
        condition_col <- "condition"
      }
      else if ("sample_type" %in% cd_cols) {
        condition_col <- "sample_type"
      }
      else {
        # Use first column as fallback
        if (length(cd_cols) > 0) {
          condition_col <- cd_cols[1]
        } else {
          stop("[plot_lm_interaction_gam_s4] No columns found in colData(se). ",
               "Cannot auto-detect condition_col.",
               call. = FALSE)
        }
      }
    }
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
    
    calculate_diversity(
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
  # CALL plot_lm_interaction_gam WITH RECONSTRUCTED DIVERSITY SE
  # =========================================================================
  result <- tryCatch({
    plot_lm_interaction_gam(
      se = diversity_combined,
      lm_res = lm_res,
      condition_col = condition_col,
      n_top = n_top,
      genes = genes,
      sig_alpha = sig_alpha,
      assay_name = assay_name,
      model_data = model_data,
      output_file = output_file,
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
  
  # Save if output_file provided
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = result, device = "pdf")
    } else if (grepl("\\.png$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = result, device = "png")
    } else if (grepl("\\.jpg$|\\.jpeg$", tolower(output_file))) {
      ggplot2::ggsave(output_file, plot = result, device = "jpeg")
    } else {
      # Default to PDF
      ggplot2::ggsave(paste0(output_file, ".pdf"), plot = result, device = "pdf")
    }

  }

  # Return the plot object directly (not the analysis object)
  result
}
