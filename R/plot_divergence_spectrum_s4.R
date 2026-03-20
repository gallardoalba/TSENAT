#' Plot Global Divergence q-Curve Across All Genes (S4 Wrapper)
#'
#' S4 wrapper for \code{\link{plot_divergence_spectrum}} that extracts
#' divergence results from a TSENATAnalysis object and visualizes the
#' average Tsallis divergence across all genes (or specified genes) as
#' a function of q-value.
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
#' @param ... Additional arguments passed to \code{\link{plot_divergence_spectrum}}.
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
#' \dontrun{
#'   # After computing divergence via S4 wrapper
#'   analysis <- calculate_divergence_s4(analysis, q = seq(0.1, 2, by=0.1))
#'
#'   # Global divergence curve (all genes aggregated) - default mode
#'   p_global <- plot_divergence_spectrum_s4(analysis)
#'   print(p_global)
#'
#'   # Top 4 genes ranked by p-value significance
#'   p_top <- plot_divergence_spectrum_s4(analysis, n_genes = 4, ncol = 2, use_pvalue_ranking = TRUE)
#'   print(p_top)
#'
#'   # Specific gene
#'   p_gene <- plot_divergence_spectrum_s4(analysis, gene = "BRCA1")
#'   print(p_gene)
#'
#'   # Save to file
#'   analysis <- plot_divergence_spectrum_s4(
#'       analysis,
#'       output_file = "divergence_spectrum.png",
#'       width = 12,
#'       height = 8
#'   )
#' }
#'
#' @seealso
#' \code{\link{plot_divergence_spectrum}} for the underlying plotting function
#' \code{\link{calculate_divergence_s4}} for computing divergence
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
        cat("Note: Invalid LM results. Plotting global curve without gene ranking.\n")
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
      cat("Error in plot_divergence_spectrum:", e$message, "\n")
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
        cat("Saved divergence spectrum plot to:", output_file, "\n")
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
