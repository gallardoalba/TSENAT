#' Filter Low-Abundance Transcripts in a TSENATAnalysis Object
#'
#' S4 wrapper for \code{filter_se()} that filters low-abundance transcripts
#' directly within a \code{TSENATAnalysis} object. This maintains the consistent
#' S4 workflow pattern where functions accept and return analysis objects.
#'
#' @param analysis A \code{TSENATAnalysis} S4 object containing the
#'   \code{SummarizedExperiment} to be filtered.
#'
#' @param stringency Character. Filtering stringency level:
#'   "strict" (most stringent, default for unpaired designs),
#'   "medium" (moderate, default for paired designs), or
#'   "lenient" (least stringent).
#'   Each keeps transcripts present in different percentages of samples with TPM >= median.
#'   Default is determined by study design (paired vs unpaired).
#'
#' @param min_samples Numeric. Minimum number of samples in which a transcript
#'   must be present (default: 5). Used as a secondary filter.
#'
#' @param verbose Logical. If TRUE, print filtering progress and summary statistics
#'   (default: FALSE).
#'
#' @return Invisibly returns the modified \code{analysis} object with filtered
#'   \code{SummarizedExperiment} in the \code{@se} slot. The filtering operation
#'   modifies the analysis object in-place while maintaining all other slots
#'   (results, metadata, etc.).
#'
#' @details
#' This wrapper applies \code{filter_se()} to the SummarizedExperiment within
#' the TSENATAnalysis object. The function:
#'
#' 1. Extracts the SE from \code{analysis@se}
#' 2. Filters using \code{filter_se()} with specified parameters
#' 3. Stores the filtered SE back in \code{analysis@se}
#' 4. Returns the modified analysis object invisibly
#'
#' **Important:** Filtering should be performed BEFORE computing diversity,
#' divergence, or LM interaction results. If called after analysis results
#' have been computed, those results will be based on unfiltered data and
#' may not align with the filtered SE dimensions.
#'
#' @seealso
#' \code{\link{filter_se}} for the underlying filtering function
#' \code{\link{build_analysis}} for creating a new analysis object
#'
#' @examples
#' \dontrun{
#' # Load transcript counts
#' data("readcounts", package = "TSENAT")
#'
#' # Create analysis object
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   tx2gene = "annotation.gff3.gz",
#'   metadata = metadata_df
#' )
#'
#' # Filter low-abundance transcripts before computing diversity
#' analysis <- filter_analysis(analysis, stringency = "medium", verbose = TRUE)
#'
#' # Now safe to compute diversity
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.1, 2, by = 0.1))
#' }
#'
#' @export
filter_analysis <- function(analysis, stringency = NULL, min_samples = 5L, verbose = FALSE) {
  # Validate input
  if (!inherits(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }

  # Extract SE from analysis
  se <- analysis@se

  # Apply filtering via filter_se
  se_filtered <- filter_se(
    se = se,
    stringency = stringency,
    min_samples = min_samples,
    verbose = verbose
  )

  # Ensure colData is preserved from original SE to maintain pairing structure
  # This is critical for downstream analyses like LM interaction tests
  if (!is.null(SummarizedExperiment::colData(se)) && 
      nrow(SummarizedExperiment::colData(se)) == ncol(se_filtered)) {
    SummarizedExperiment::colData(se_filtered) <- SummarizedExperiment::colData(se)
  }

  # Store filtered SE back in analysis object
  analysis@se <- se_filtered

  # Return modified analysis object
  analysis
}
