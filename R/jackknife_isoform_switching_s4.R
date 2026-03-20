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
#'   The analysis object is returned invisibly to support method chaining:
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
    
    # Try Priority 2: @config$sample_type
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
      stop("[jackknife_isoform_switching_s4] Cannot auto-detect condition_col. ",
           "Provide explicitly or ensure @se has colData with sample groupings.",
           call. = FALSE)
    }
    
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Auto-detected condition_col =", condition_col, "\n")
    }
  }
  
  # Validate condition_col exists
  if (!(condition_col %in% colnames(colData(se)))) {
    stop("[jackknife_isoform_switching_s4] Specified condition_col='", condition_col,
         "' not found in colData. Available columns: ",
         paste(colnames(colData(se)), collapse = ", "),
         call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT gene_col AND isoform_col FROM rowData or NAMESPACE
  # =========================================================================
  rd <- if (!is.null(rowData(se)) && nrow(rowData(se)) > 0) {
    rowData(se)
  } else {
    NULL
  }
  
  # Detect gene_col
  if (is.null(gene_col)) {
    if (!is.null(rd)) {
      rd_cols <- colnames(rd)
      if ("gene_id" %in% rd_cols) {
        gene_col <- "gene_id"
      } else if ("gene" %in% rd_cols) {
        gene_col <- "gene"
      } else if ("Gene" %in% rd_cols) {
        gene_col <- "Gene"
      }
    }
    
    # If not found in rowData, use default
    if (is.null(gene_col)) {
      gene_col <- "gene"
    }
    
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Using gene_col =", gene_col, "\n")
    }
  }
  
  # Detect isoform_col
  if (is.null(isoform_col)) {
    if (!is.null(rd)) {
      rd_cols <- colnames(rd)
      if ("transcript_id" %in% rd_cols) {
        isoform_col <- "transcript_id"
      } else if ("transcript" %in% rd_cols) {
        isoform_col <- "transcript"
      } else if ("isoform" %in% rd_cols) {
        isoform_col <- "isoform"
      } else if ("Isoform" %in% rd_cols) {
        isoform_col <- "Isoform"
      } else if ("tx_id" %in% rd_cols) {
        isoform_col <- "tx_id"
      }
    }
    
    # If not found in rowData, use default
    if (is.null(isoform_col)) {
      isoform_col <- "transcript"
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
  # STORE RESULTS IN ANALYSIS OBJECT
  # =========================================================================
  # Check if result has multi-q class (when multiple q-values provided)
  if (inherits(result, "tsenat_isoform_switching_multiq")) {
    # Multi-q result: store as-is to preserve class and structure
    analysis@jackknife_results[["multi_q"]] <- result
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Stored multi-q result with special class\n")
    }
  } else if (is.numeric(q) && length(q) > 1) {
    # Multiple q-values but result is NOT multi-q class: store each individually
    for (q_val in q) {
      q_key <- paste0("q_", sprintf("%.2f", q_val))
      
      if (is.list(result) && q_key %in% names(result)) {
        analysis@jackknife_results[[q_key]] <- result[[q_key]]
      } else {
        analysis@jackknife_results[[q_key]] <- result
      }
      
      if (verbose) {
        cat("[jackknife_isoform_switching_s4] Stored results for", q_key, "\n")
      }
    }
  } else {
    # Single q-value: store with q-value key
    q_val <- q[1]
    q_key <- paste0("q_", sprintf("%.2f", q_val))
    analysis@jackknife_results[[q_key]] <- result
    if (verbose) {
      cat("[jackknife_isoform_switching_s4] Stored results for", q_key, "\n")
    }
  }
  
  # =========================================================================
  # TRACK FUNCTION CALL IN METADATA
  # =========================================================================
  if (is.list(analysis@metadata)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("jackknife_isoform_switching_s4[q=", paste(q, collapse = ","),
             ", condition_col=", condition_col, "]")
    )
    analysis@metadata$function_timestamps <- c(
      analysis@metadata$function_timestamps,
      as.character(Sys.time())
    )
  }
  
  # Return modified analysis object (invisibly for chaining)
  invisible(analysis)
}
