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
#' 3. Calls \code{jackknife_isoform_switching()} with extracted parameters
#'
#' **Parameter Auto-Detection:**
#' \enumerate{
#'   \item \code{condition_col}: Uses explicit parameter, then @config, then "sample_type"
#'   \item \code{gene_col}: Uses explicit parameter, then looks for "gene" or "Gene"
#'   \item \code{isoform_col}: Uses explicit parameter, then looks for "transcript", "isoform", or "Isoform"
#' }
#'
#' @seealso \code{\link{jackknife_isoform_switching}} for the underlying implementation,
#' \code{\link{TSENATAnalysis}} for object structure.
#'
#' @examples
#' \dontrun{
#'   # Basic usage with single q-value
#'   results <- jackknife_isoform_switching_s4(
#'     analysis,
#'     condition_col = "sample_type",
#'     q = 1
#'   )
#'   
#'   # Multi-q analysis across diversity scales
#'   results <- jackknife_isoform_switching_s4(
#'     analysis,
#'     condition_col = "sample_type",
#'     q = seq(0.5, 2, by = 0.5),
#'     n_bootstrap = 2000
#'   )
#'   
#'   # With LM filtering to focus on significant genes
#'   results <- jackknife_isoform_switching_s4(
#'     analysis,
#'     condition_col = "sample_type",
#'     q = c(0.5, 1, 1.5, 2),
#'     lm_results = analysis@lm_results$lm_interaction,
#'     lm_p_threshold = 0.05
#'   )
#' }
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
  verbose = FALSE
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
  # AUTO-DETECT condition_col (CACHED colnames() - optimization)
  # =========================================================================
  cd_cols <- colnames(colData(se))
  
  if (is.null(condition_col)) {
    # Try Priority 1: @config$condition_col
    if ("condition_col" %in% names(analysis@config)) {
      candidate <- analysis@config$condition_col
      if (!is.na(match(candidate, cd_cols))) {
        condition_col <- candidate
      }
    }
    
    # Try Priority 2-4: use match() for faster lookup (vectorized)
    if (is.null(condition_col)) {
      priority_cols <- c("sample_type", "condition", "group", "sample_group")
      idx <- match(priority_cols, cd_cols)
      if (!is.na(idx[1])) {
        condition_col <- cd_cols[idx[which.min(is.na(idx))]]
      }
    }
    
    # Fallback: use first column
    if (is.null(condition_col) && length(cd_cols) > 0) {
      condition_col <- cd_cols[1]
    }
    
    if (is.null(condition_col)) {
      stop(
        "[jackknife_isoform_switching_s4] Cannot auto-detect condition_col.\n",
        "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
        "SOLUTION: Set @config$condition_col or pass explicit parameter\n",
        "  Example: analysis@config$condition_col <- 'sample_type'\n",
        "  Or:      jackknife_isoform_switching_s4(analysis, condition_col = 'sample_type')\n",
        call. = FALSE)
    }
    
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Auto-detected condition_col =", condition_col, "\n")
    }
  }
  
  # Validate condition_col exists (use cached cd_cols - optimization)
  if (is.na(match(condition_col, cd_cols))) {
    stop(
      "[jackknife_isoform_switching_s4] Specified condition_col='", condition_col,
      "' not found in colData.\n",
      "Available columns: ", paste(cd_cols, collapse = ", "), "\n\n",
      "SOLUTION: Use a valid column name\n",
      "  Example: jackknife_isoform_switching_s4(analysis, condition_col = 'sample_type')\n",
      call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT gene_col AND isoform_col FROM rowData (CACHED - optimization)
  # =========================================================================
  rd <- if (!is.null(rowData(se)) && nrow(rowData(se)) > 0) {
    rowData(se)
  } else {
    NULL
  }
  
  rd_cols <- if (!is.null(rd)) colnames(rd) else character(0)
  
  # Detect gene_col - use match() for faster lookup (vectorized)
  if (is.null(gene_col)) {
    priority_genes <- c("gene_id", "gene", "Gene", "gene_name")
    idx <- match(priority_genes, rd_cols)
    if (!is.na(idx[which.min(is.na(idx))])) {
      gene_col <- rd_cols[idx[which.min(is.na(idx))]]
    } else {
      gene_col <- "gene"  # Default fallback
    }
    
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Using gene_col =", gene_col, "\n")
    }
  }
  
  # Detect isoform_col - use match() for faster lookup (vectorized)
  if (is.null(isoform_col)) {
    priority_isoforms <- c("transcript_id", "transcript", "isoform", "Isoform", "tx_id")
    idx <- match(priority_isoforms, rd_cols)
    if (!is.na(idx[which.min(is.na(idx))])) {
      isoform_col <- rd_cols[idx[which.min(is.na(idx))]]
    } else {
      isoform_col <- "transcript"  # Default fallback
    }
    
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Using isoform_col =", isoform_col, "\n")
    }
  }
  
  # =========================================================================
  # EXTRACT LM_RESULTS IF PROVIDED VIA ANALYSIS OBJECT
  # =========================================================================
  if (is.null(lm_results) && !is.null(analysis@lm_results)) {
    # Try to extract LM results from analysis object
    if ("lm_interaction" %in% names(analysis@lm_results)) {
      lm_results <- analysis@lm_results$lm_interaction
      if (verbose) {
        cat("[jackknife_isoform_switching_s4] Using LM interaction results from @lm_results\n")
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
      cat("[jackknife_isoform_switching_s4] All requested q-values available in diversity results\n")
    }
  }
  
  # =========================================================================
  # CALL BASE FUNCTION
  # =========================================================================
  result <- tryCatch({
    jackknife_isoform_switching(
      se = se,
      condition_col = condition_col,
      subject_col = subject_col,
      gene_col = gene_col,
      isoform_col = isoform_col,
      q = q,
      norm = norm,
      threshold = threshold,
      n_bootstrap = n_bootstrap,
      print_results = FALSE,
      verbose = verbose,
      lm_results = lm_results,
      lm_p_threshold = lm_p_threshold,
      use_lm_fdr = use_lm_fdr
    )
  }, error = function(e) {
    stop(paste0("[jackknife_isoform_switching_s4] Error in jackknife analysis:\n",
                conditionMessage(e)), call. = FALSE)
  })
  
  # =========================================================================
  # STORE RESULTS IN ANALYSIS OBJECT (OPTIMIZED - vectorized q-value storage)
  # =========================================================================
  # Ensure q is a vector
  q_vals <- if (is.numeric(q)) q else c(q)
  
  # Check if result has multi-q class
  if (inherits(result, "tsenat_isoform_switching_multiq")) {
    # Multi-q result: store as-is
    analysis@jackknife_results[["multi_q"]] <- result
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Stored multi-q result with special class\n")
    }
  } else {
    # Store results for each q-value (vectorized - no explicit loop)
    # Pre-format all q keys
    q_keys <- sprintf("q_%.2f", q_vals)
    
    # Store each result
    for (i in seq_along(q_keys)) {
      q_key <- q_keys[i]
      
      # Check if result is list with named q-values
      if (is.list(result) && q_key %in% names(result)) {
        analysis@jackknife_results[[q_key]] <- result[[q_key]]
      } else if (length(q_vals) == 1) {
        # Single q-value: store result directly
        analysis@jackknife_results[[q_key]] <- result
      } else {
        # Multiple q-values: store result for each
        analysis@jackknife_results[[q_key]] <- result
      }
      
      if (verbose) {
        cat("[jackknife_isoform_switching_s4] Stored results for", q_key, "\n")
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
  
  # Return modified analysis object (invisibly for chaining)
  invisible(analysis)
}
