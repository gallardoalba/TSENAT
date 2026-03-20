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
#' \dontrun{
#'   # After computing concordance:
#'   analysis <- compute_method_concordance_s4(analysis)
#'   
#'   # Generate plot:
#'   plot <- plot_method_concordance_s4(analysis)
#'   print(plot)
#' }
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
    cat("[plot_method_concordance_s4] Plotting concordance for",
        nrow(comparison_df), "genes\n")
    cat("[plot_method_concordance_s4] Methods compared:",
        concordance_results$gam_method, "vs",
        concordance_results$friedman_method, "\n")
  }
  
  # Call standard plotting function
  plot_obj <- plot_method_concordance(comparison_df)
  
  if (verbose) {
    cat("[plot_method_concordance_s4] Plot generated successfully\n")
  }
  
  return(plot_obj)
})
