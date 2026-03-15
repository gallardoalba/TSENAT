#' Unified Batch Correction Interface
#'
#' Apply batch correction using one of several statistical methods.
#' This wrapper provides user-friendly access to various batch effect
#' correction approaches while hiding method-specific implementation details.
#'
#' @param x A SummarizedExperiment object containing expression/count data
#' @param batch_column Character specifying column name in colData(x) for batch assignment
#' @param method Batch correction method to apply. Options:
#'   - `"combat"` (default): ComBat using sequential adjustment (ComBat-seq variant)
#'   - `"combat_seq"`: ComBat-seq for count data with variance inflation modeling
#'   - `"combat_ref"`: ComBat with designated reference batch
#'   - `"ranking"`: Rank-based correction via linear modeling
#' @param reference_batch (Optional) Name or index of reference batch for `method = "combat_ref"`.
#'   If NULL, first batch level used automatically.
#' @param ... Additional arguments passed to the specific batch correction method
#'
#' @return SummarizedExperiment with batch-corrected data stored in "batch_corrected" assay.
#'   Original colData and rowData are preserved. Metadata includes method, date, and batch_column.
#'
#' @details
#' ## Method Selection
#'
#' **ComBat (Default):** General-purpose batch correction using empirical Bayes shrinkage.
#' Uses sequential ComBat-seq approach to estimate and remove batch-specific shifts and scaling.
#' Works well with various data types.
#'
#' **ComBat-seq:** Designed specifically for RNA-seq count data. Models mean-variance
#' relationship and count-specific properties. Recommended when working with raw counts
#' or pseudo-counts.
#'
#' **ComBat-ref:** Variant of ComBat that preserves inter-batch relationships while
#' correcting toward a designated reference batch. Useful when one batch has special
#' biological significance or should be preserved.
#'
#' **Ranking:** Non-parametric approach using rank-based linear models. Robust to
#' distribution assumptions. Returns relative ranks rather than adjusted values.
#'
#' ## Backward Compatibility
#'
#' This function consolidates workflow from multiple batch correction functions:
#' - `adjust_batch_effects()` → use `method = "combat"` or `"combat_ref"`
#' - `adjust_batch_effects_seq()` → use `method = "combat_seq"`
#' - `apply_batch_correction_ranking()` → use `method = "ranking"`
#'
#' @seealso
#'   [detect_batch_effects()] for batch effect detection and visualization
#'   [detect_batch_structure()] for structural batch analysis via PCA
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate sample SummarizedExperiment
#' library(SummarizedExperiment)
#' set.seed(42)
#' counts <- matrix(rnbinom(500, mu = 100, size = 2), nrow = 50, ncol = 20)
#' rownames(counts) <- paste0("gene_", 1:50)
#' colnames(counts) <- paste0("sample_", 1:20)
#' batch <- factor(rep(c("batch1", "batch2"), each = 10))
#'
#' se <- SummarizedExperiment(
#'   assays = list(counts = counts),
#'   colData = DataFrame(batch = batch)
#' )
#'
#' # Apply ComBat correction (default)
#' se_corrected <- correct_batch_effects(
#'   x = se,
#'   batch_column = "batch",
#'   method = "combat"
#' )
#'
#' # Apply ComBat-seq for count data
#' se_corrected_seq <- correct_batch_effects(
#'   x = se,
#'   batch_column = "batch",
#'   method = "combat_seq"
#' )
#'
#' # Apply ComBat with reference batch
#' se_corrected_ref <- correct_batch_effects(
#'   x = se,
#'   batch_column = "batch",
#'   method = "combat_ref",
#'   reference_batch = "batch1"
#' )
#'
#' # Apply ranking-based method
#' se_corrected_ranking <- correct_batch_effects(
#'   x = se,
#'   batch_column = "batch",
#'   method = "ranking"
#' )
#'
#' # Access corrected assay
#' corrected_matrix <- assay(se_corrected, "batch_corrected")
#' }

correct_batch_effects <- function(
    x,
    batch_column,
    method = c("combat", "combat_seq", "combat_ref", "ranking"),
    reference_batch = NULL,
    ...) {

  # Validate and resolve batch column
  method <- match.arg(method)

  # Input must be SummarizedExperiment
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("x must be a SummarizedExperiment object", call. = FALSE)
  }

  # Validate batch column exists
  if (!is.character(batch_column) || !(batch_column %in% names(SummarizedExperiment::colData(x)))) {
    stop("batch_column '", batch_column,
         "' not found in colData(x)", call. = FALSE)
  }

  # Apply appropriate batch correction method
  result_se <- switch(method,
    "combat" = .correct_batch_combat(x, batch_column, ...),
    "combat_seq" = .correct_batch_combat_seq(x, batch_column, ...),
    "combat_ref" = .correct_batch_combat_ref(
      x, batch_column, reference_batch, ...),
    "ranking" = .correct_batch_ranking(x, batch_column, ...),
    stop("Unknown method: ", method, call. = FALSE)
  )

  # Add metadata about correction
  if (is.null(S4Vectors::metadata(result_se))) {
    S4Vectors::metadata(result_se) <- list()
  }
  S4Vectors::metadata(result_se)$batch_correction_method <- method
  S4Vectors::metadata(result_se)$batch_correction_date <- Sys.Date()
  S4Vectors::metadata(result_se)$batch_column <- batch_column

  return(result_se)
}

################################################################################
#
#' Internal: ComBat Batch Correction
#'
#' @noRd
.correct_batch_combat <- function(se, batch_column, ...) {
  # Dispatch to adjust_batch_effects() with method = "standard" (ComBat-seq)
  result_se <- adjust_batch_effects(
    se = se,
    batch = batch_column,
    method = "standard",
    ...
  )
  
  # adjust_batch_effects modifies the "counts" assay in place
  # We need to rename it to "batch_corrected" for consistency
  if ("counts" %in% names(assays(result_se))) {
    assays(result_se)$batch_corrected <- assays(result_se)$counts
  }
  
  return(result_se)
}

################################################################################
#
#' Internal: ComBat-seq Batch Correction
#'
#' @noRd  
.correct_batch_combat_seq <- function(se, batch_column, ...) {
  # Dispatch to adjust_batch_effects_seq()
  result_se <- adjust_batch_effects_seq(
    se = se,
    batch = batch_column,
    ...
  )
  
  # adjust_batch_effects_seq modifies the "counts" assay in place
  # We need to rename it to "batch_corrected" for consistency
  if ("counts" %in% names(assays(result_se))) {
    assays(result_se)$batch_corrected <- assays(result_se)$counts
  }
  
  return(result_se)
}

################################################################################
#
#' Internal: ComBat with Reference Batch
#'
#' @noRd
.correct_batch_combat_ref <- function(se, batch_column, reference_batch = NULL, ...) {
  # Dispatch to adjust_batch_effects() with method = "ref"
  result_se <- adjust_batch_effects(
    se = se,
    batch = batch_column,
    method = "ref",
    ref_batch = reference_batch,
    ...
  )
  
  # adjust_batch_effects modifies the "counts" assay in place
  # We need to rename it to "batch_corrected" for consistency
  if ("counts" %in% names(assays(result_se))) {
    assays(result_se)$batch_corrected <- assays(result_se)$counts
  }
  
  return(result_se)
}

################################################################################
#
#' Internal: Rank-based Batch Correction
#'
#' @noRd
.correct_batch_ranking <- function(se, batch_column, ...) {
  # Dispatch to apply_batch_correction_ranking.SummarizedExperiment
  result_se <- apply_batch_correction_ranking(
    entropy_matrix = se,
    batch_column = batch_column,
    ...
  )
  
  # Get the corrected matrix from the first assay (which was modified)
  # and move it to batch_corrected
  corrected_matrix <- SummarizedExperiment::assay(result_se, 1)
  
  # Reset first assay to original, then add batch_corrected
  SummarizedExperiment::assay(result_se, 1) <- SummarizedExperiment::assay(se, 1)
  SummarizedExperiment::assays(result_se)$batch_corrected <- corrected_matrix
  
  return(result_se)
}
