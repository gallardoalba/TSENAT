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
#' \dontrun{
#'   # After running full analysis pipeline
#'   analysis <- jackknife_isoform_switching_s4(analysis, q = c(0.5, 1, 1.5, 2))
#'
#'   # Generate heatmaps using S4 wrapper
#'   heatmap_file <- plot_multiq_delta_influence_heatmaps_s4(
#'     analysis,
#'     n_genes = 4
#'   )
#'
#'   # Load and display the heatmap
#'   library(magick)
#'   heatmap_img <- image_read(heatmap_file)
#'   print(heatmap_img)
#' }
#'
#' @seealso \code{\link{plot_multiq_delta_influence_heatmaps}} for the base function,
#' \code{\link{jackknife_isoform_switching_s4}} for computing switching results
#'
#' @export
#' @importFrom methods is
plot_multiq_delta_influence_heatmaps_s4 <- function(
    analysis,
    n_genes = 4,
    lm_results = NULL,
    verbose = FALSE) {
  
  # Validation
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }
  
  if (verbose) cat("Extracting jackknife switching results from analysis object...\n")
  
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
      cat("  ✓ Found multi-q result with class:", class(switching_results)[1], "\n")
    }
  } else {
    # Fallback: use all results as list (for single or multiple q-values)
    switching_results <- jackknife_results_list
    if (verbose) {
      cat("  ✓ Using individual q-value results (", length(switching_results), "q-values)\n", sep = "")
    }
  }
  
  if (verbose) {
    cat("  Q-values:", paste(names(switching_results), collapse = ", "), "\n")
  }
  
  # Extract LM results if not provided
  if (is.null(lm_results)) {
    if (verbose) cat("Extracting LM results from analysis@lm_results...\n")
    
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
        if (verbose) cat("  ✓ Extracted LM results with", nrow(lm_results), "genes\n")
      } else if (verbose) {
        cat("  ⚠ LM results not found; genes will be ranked by appearance\n")
      }
    }
  }
  
  if (verbose) cat("Calling plot_multiq_delta_influence_heatmaps()...\n")
  
  # Call base function with extracted parameters
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = switching_results,
    n_genes = n_genes,
    lm_results = lm_results
  )
  
  if (verbose) {
    cat("✓ Heatmap plot generated successfully\n")
    cat("  Saved to:", heatmap_file, "\n")
  }
  
  invisible(heatmap_file)
}
