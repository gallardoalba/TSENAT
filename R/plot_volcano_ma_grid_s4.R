#' Plot Volcano and MA Grid from Differential Analysis Results (S4 Wrapper)
#'
#' S4 wrapper for \code{\link{plot_volcano_ma_grid}} that extracts differential
#' analysis results from a TSENATAnalysis object and creates side-by-side volcano
#' and MA plots for comparing control and treatment groups.
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
#' @param verbose \code{logical}. Print status messages. Default: TRUE.
#' @param ... Additional arguments passed to \code{\link{plot_volcano_ma_grid}}.
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
#' \dontrun{
#'   # After calculating diversity and differences
#'   analysis <- calculate_diversity_s4(analysis, q = 1.0)
#'   analysis <- calculate_difference_s4(analysis, control = "Normal")
#'
#'   # Create volcano and MA plot grid (default color threshold: p=0.05)
#'   plot_grid <- plot_volcano_ma_grid_s4(analysis)
#'   print(plot_grid)
#'
#'   # Customize thresholds and labels
#'   plot_custom <- plot_volcano_ma_grid_s4(
#'     analysis,
#'     sig_alpha = 0.01,           # Stricter significance threshold
#'     label_thresh = 0.05,         # Label genes with padj < 0.05
#'     top_n = 10,                 # Label top 10 significant genes
#'     title_volcano = "Volcano Plot: Normal vs Tumor",
#'     title_ma = "MA Plot: Normal vs Tumor"
#'   )
#'   print(plot_custom)
#' }
#'
#' @seealso
#' \code{\link{plot_volcano_ma_grid}} for the underlying plotting function
#' \code{\link{calculate_difference_s4}} for computing differential analysis
#'
#' @import methods
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
    verbose = TRUE,
    ...) {

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

  if (verbose) {
    cat("[plot_volcano_ma_grid_s4] Found", nrow(diff_df), "genes in difference results\n")
    cat("[plot_volcano_ma_grid_s4] Columns available:", paste(colnames(diff_df), collapse = ", "), "\n")
  }

  # Validate padj_col exists
  # Handle both "padj" and "adjusted_p_values" column names
  actual_padj_col <- padj_col
  if (!(padj_col %in% colnames(diff_df))) {
    # Try alternative column name
    if ("adjusted_p_values" %in% colnames(diff_df) && padj_col == "padj") {
      actual_padj_col <- "adjusted_p_values"
      if (verbose) {
        cat("[plot_volcano_ma_grid_s4] Using 'adjusted_p_values' instead of 'padj'\n")
      }
    } else if ("pvalue" %in% colnames(diff_df) && padj_col == "padj") {
      # Fallback to raw p-values if adjusted not available
      actual_padj_col <- "pvalue"
      if (verbose) {
        cat("[plot_volcano_ma_grid_s4] Using 'pvalue' as padj_col (adjusted p-values not found)\n")
      }
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
      if (verbose) {
        cat("[plot_volcano_ma_grid_s4] Using 'mean_difference' for x-axis\n")
      }
    } else if ("log2_fold_change" %in% colnames(diff_df)) {
      x_col <- "log2_fold_change"
      if (verbose) {
        cat("[plot_volcano_ma_grid_s4] Using 'log2_fold_change' for x-axis\n")
      }
    } else {
      warning("Could not auto-detect x_col. Available numeric columns: ",
              paste(colnames(diff_df)[sapply(diff_df, is.numeric)], collapse = ", "),
              call. = FALSE)
    }
  }

  # Create the plot
  if (verbose) {
    cat("[plot_volcano_ma_grid_s4] Creating volcano and MA plot grid...\n")
    cat("  x_col:", x_col, "\n")
    cat("  padj_col:", actual_padj_col, "\n")
    cat("  sig_alpha:", sig_alpha, "\n")
    cat("  label_thresh:", label_thresh, "\n")
    cat("  top_n:", top_n, "\n")
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
    stop("[plot_volcano_ma_grid_s4] Error creating plot:\n", conditionMessage(e),
         call. = FALSE)
  })

  if (verbose) {
    cat("[plot_volcano_ma_grid_s4] Plot created successfully\n")
  }

  return(invisible(plot_obj))
}
