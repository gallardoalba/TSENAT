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
#' @param verbose \code{logical}. Print progress messages. Default: FALSE
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
#' \dontrun{
#'   # After running both GAM and Friedman analyses:
#'   analysis <- detect_q_gene_interactions_s4(analysis, ...)
#'   analysis <- test_rankbased_assumptions_s4(analysis, ...)
#'   
#'   # Compute concordance:
#'   analysis <- compute_method_concordance_s4(analysis)
#'   
#'   # Access results:
#'   concordance_results <- analysis@metadata$method_concordance
#'   cat("Spearman correlation:", concordance_results$spearman_rho, "\n")
#'   print(concordance_results$agreement_table)
#' }
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
    verbose = FALSE) {
  
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
    stop(paste0("GAM method '", gam_method, "' not found in LM results. ",
                "Available: ", available_methods), call. = FALSE)
  }
  
  if (!(friedman_method %in% names(analysis@lm_results))) {
    available_methods <- paste(names(analysis@lm_results), collapse = ", ")
    stop(paste0("Friedman method '", friedman_method, "' not found in LM results. ",
                "Available: ", available_methods), call. = FALSE)
  }
  
  # Extract results
  gam_results <- analysis@lm_results[[gam_method]]
  friedman_results <- analysis@lm_results[[friedman_method]]
  
  # Validate they're data frames
  if (!is.data.frame(gam_results)) {
    stop(paste0("GAM results ('", gam_method, "') must be a data.frame"), call. = FALSE)
  }
  
  if (!is.data.frame(friedman_results)) {
    stop(paste0("Friedman results ('", friedman_method, "') must be a data.frame"), call. = FALSE)
  }
  
  # ===================================================================
  # COMPUTE CONCORDANCE
  # ===================================================================
  
  if (verbose) {
    cat("[compute_method_concordance_s4] Computing concordance between ",
        gam_method, " and ", friedman_method, "\n", sep = "")
  }
  
  # Call the standard function
  concordance_result <- tryCatch({
    compute_method_concordance(gam_results, friedman_results)
  }, error = function(e) {
    stop(paste0("[compute_method_concordance_s4] Error computing concordance:\n",
                conditionMessage(e)), call. = FALSE)
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
    cat("[compute_method_concordance_s4] Concordance computed successfully\n")
    if (!is.na(concordance_result$spearman_rho)) {
      cat("[compute_method_concordance_s4] Spearman corr =", 
          round(concordance_result$spearman_rho, 3), "\n")
    }
  }
  
  analysis
}
)
