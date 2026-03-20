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
#' @param ... Additional arguments passed to \code{\link{effect_sizes_divergence}}.
#'
#' @return Modified TSENATAnalysis with effect size results stored in
#'   \code{@metadata$effect_sizes_divergence}. Returns the analysis object invisibly
#'   to support piping and method chaining.
#'
#' @details
#' This wrapper:
#' \describe{
#'   \item{Extracting}{Divergence SE from \code{@divergence_results} and LM results
#'     from \code{@lm_results$lm_interaction}}
#'   \item{Computing}{Effect sizes using standard \code{effect_sizes_divergence()} function}
#'   \item{Storing}{Results as list with \code{interaction_results} (data.frame) and
#'     \code{validation_stats}}
#'   \item{Tracking}{Function call in \code{@metadata$function_calls}}
#' }
#'
#' Results are accessed via: \code{analysis@metadata$effect_sizes_divergence}
#'
#' @examples
#' \dontrun{
#'   # Compute effect sizes after divergence and LM analysis
#'   analysis <- effect_sizes_divergence_s4(analysis, significance_threshold = 0.05)
#'
#'   # Extract results
#'   eff_res <- analysis@metadata$effect_sizes_divergence$interaction_results
#' }
#'
#' @seealso \code{\link{effect_sizes_divergence}} for the base function,
#' \code{\link{calculate_divergence_s4}} for divergence wrapper,
#' \code{\link{calculate_lm_interaction_s4}} for LM interaction wrapper
#'
#' @export
#' @importFrom methods is
effect_sizes_divergence_s4 <- function(
    analysis,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE,
    verbose = TRUE,
    ...) {

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
    cat("[effect_sizes_divergence_s4] Extracted:\n")
    cat("  - Divergence SE:", paste(dim(divergence_se), collapse = " x "), "\n")
    cat("  - LM results: ", nrow(lm_res), " genes\n", sep = "")
  }

  # =========================================================================
  # COMPUTE EFFECT SIZES
  # =========================================================================
  if (verbose) {
    cat("[effect_sizes_divergence_s4] Computing effect sizes...\n")
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
    stop(paste0("Error in effect_sizes_divergence:\n", e$message),
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
    cat("[effect_sizes_divergence_s4] Results stored in @metadata$effect_sizes_divergence\n")
    if (!is.null(result$interaction_results)) {
      cat("  - Effect size results:", nrow(result$interaction_results), "genes\n")
    }
  }

  invisible(analysis)
}
