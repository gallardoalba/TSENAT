#' Plot method concordance results from TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with computed method concordance
#'   (from \code{calculate_concordance()}).
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
#' Requires that \code{calculate_concordance()} has already been run
#' to populate \code{@metadata$method_concordance}.
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
#' = 200)
#'
#' # Note: calculate_concordance requires additional LM and Conover-Iman Rank Transform
#' # results computed. For demo purposes, we show that
#' # plot_concordance needs pre-computed concordance in @metadata
#'
#' @aliases plot_concordance
#' @export
setGeneric("plot_concordance", function(analysis, verbose = FALSE) {
    standardGeneric("plot_concordance")
})

#' @rdname plot_concordance
setMethod("plot_concordance", "TSENATAnalysis", function(analysis, verbose = FALSE) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Validate that concordance results exist
    if (is.null(analysis@metadata$method_concordance)) {
        stop("[plot_concordance] No concordance results found in @metadata.\n", "  Please run calculate_concordance() first.")
    }

    concordance_results <- analysis@metadata$method_concordance

    # Extract comparison dataframe
    comparison_df <- concordance_results$comparison_df

    if (is.null(comparison_df) || nrow(comparison_df) == 0) {
        stop("[plot_concordance] Concordance comparison_df is empty or missing.")
    }

    if (verbose) {
        message("[plot_concordance] Plotting concordance for ", nrow(comparison_df),
            " genes")
        message("[plot_concordance] Methods compared: ", concordance_results$sait_method,
            " vs ", concordance_results$rank_method)
    }

    # Call standard plotting function
    plot_obj <- .plot_concordance(comparison_df)

    if (verbose) {
        message("[plot_method_concordance_s4] Plot generated successfully")
    }

    return(plot_obj)
})

# ============================================================================
# OPTIMIZATION: Consolidated object extraction helper
# ============================================================================
# This helper consolidates redundant fallback extraction patterns into a single
# source of truth, reducing code duplication and improving maintainability.
#' @noRd
.extract_object_with_fallbacks <- function(obj, expected_class, key_name = NULL,
    verbose = FALSE) {
    # Single source of extraction logic for common pattern: Try direct class
    # match, then named list access, then list[1]

    if (is(obj, expected_class)) {
        if (verbose) {
            message("[extract_object] Found object via direct class match: ", expected_class)
        }
        return(obj)
    }

    if (is.list(obj)) {
        # Try named access first
        if (!is.null(key_name) && key_name %in% names(obj)) {
            if (verbose) {
                message("[extract_object] Found object via key: ", key_name)
            }
            return(obj[[key_name]])
        }

        # Fall back to first element
        if (length(obj) > 0) {
            if (verbose) {
                message("[extract_object] Using first element of list")
            }
            return(obj[[1]])
        }
    }

    if (verbose) {
        message("[extract_object] Could not extract object of class ", expected_class)
    }
    return(NULL)
}

