#' Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)
#'
#' S4 wrapper for \code{\link{plot_divergence_distribution}} that extracts
#' effect size results from a TSENATAnalysis object and generates a histogram
#' visualization of Tsallis divergence effect sizes across genes.
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
#' @param ... Additional arguments passed to \code{\link{plot_divergence_distribution}}.
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
#' \dontrun{
#'   # After computing effect sizes via S4 wrapper
#'   analysis <- effect_sizes_divergence_s4(analysis)
#'
#'   # Generate and display the plot
#'   analysis <- plot_divergence_distribution_s4(
#'       analysis,
#'       threshold = 0.1,
#'       output_file = "divergence_distribution.png"
#'   )
#' }
#'
#' @seealso
#' \code{\link{plot_divergence_distribution}} for the underlying plotting function
#' \code{\link{effect_sizes_divergence_s4}} for computing effect sizes
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
      cat("Error in plot_divergence_distribution:", e$message, "\n")
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
        cat("Saved divergence distribution plot to:", output_file, "\n")
      }
    }, error = function(e) {
      if (verbose) {
        cat("Warning: Could not save plot to file:", e$message, "\n")
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
