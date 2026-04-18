#' TSENATAnalysis S4 Class Definition
#'
#' Central container for unified analysis workflows in TSENAT.
#'
#' The \code{TSENATAnalysis} class encapsulates all components of a complete
#' TSENAT analysis: raw data, configuration metadata, and results from each
#' analytical step. This unified object ensures metadata is never lost through
#' the analysis pipeline and provides consistent accessor methods for result
#' retrieval.
#'
#' @slot se \code{SummarizedExperiment}. The base expression data object
#'   (genes x samples) with assays and colData.
#'
#' @slot config \code{list}. Configuration metadata specifying analysis
#'   parameters that persist through the workflow (q-values, sample grouping
#'   columns, etc.). Set once via \code{TSENAT_config()} and used by all
#'   downstream wrapper functions.
#'
#' @slot diversity_results \code{list}. Named list of diversity calculation
#'   results. Each name corresponds to a q-value (e.g., 'q_0.5', 'q_1.0').
#'   Values are SummarizedExperiment objects or data.frames containing entropy
#'   values for each gene at that q-value.
#'
#' @slot lm_results \code{list}. Complex results from regularized/penalized regression methods
#'   (GAM, LMM, GEE, FPCA) and statistical testing. Top-level names identify analysis type:
#'   \describe{
#'     \item{\code{lm_interaction}}{Regularized regression (GAM/LMM/GEE/FPCA) model results (list with
#'           \code{$results} data.frame, \code{$models} list, etc.)}
#'     \item{\code{rank_test}}{Scheirer-Ray-Hare rank-based test results}
#'     \item{\code{divergence_difference}}{Differential divergence comparison}
#'   }
#'
#' @slot pairwise_results \code{list}. Pairwise group comparison results.
#'   Results from running statistical tests for group comparisons.
#'
#' @slot rank_test_results \code{list}. Scheirer-Ray-Hare rank-based statistical
#'   test results. Names correspond to q-values (e.g., 'q_0.5', 'q_1.0').
#'   Computed by \code{calculate_srh()} as a non-parametric alternative
#'   to linear mixed model testing.
#'
#' @slot jackknife_results \code{list}. Resampling-based confidence intervals.
#'   Names correspond to q-values (e.g., 'q_0.5', 'q_1.0'). Values are
#'   jackknife result objects containing resamples, CI bounds, and diagnostics.
#'
#' @slot divergence_results \code{list}. Divergence metric calculations.
#'   Typically contains:
#'   \describe{
#'     \item{\code{tsallis_divergence}}{SummarizedExperiment with 
#' divergence values}
#'     \item{\code{effect_sizes}}{data.frame with Cohen's d, etc.}
#'   }
#'
#' @slot plots \code{list}. Cached visualization objects (ggplot). Names
#'   identify plot type (e.g., 'q_curve', 'lm_interaction', 'influence').
#'   Populated by \code{TSENAT()} if \code{generate_plots=TRUE}.
#'
#' @slot metadata \code{list}. Reproducibility and tracking metadata.
#'   Automatically maintained by wrapper functions. Includes:
#'   \describe{
#'     \item{\code{created_at}}{Timestamp of object creation}
#'     \item{\code{function_calls}}{Vector of wrapper functions called}
#'     \item{\code{function_timestamps}}{Timestamps for each function call}
#'     \item{\code{package_version}}{TSENAT version at creation}
#'   }
#'
#' @details
#' Access results via the unified \code{results(obj, type = ...)} accessor method:
#' - \code{type='diversity'} for Tsallis entropy across q-values
#' - \code{type='lm'} for regularized/penalized regression (GAM, LMM, GEE, FPCA) interaction results
#' - \code{type='rank_test'} for Scheirer-Ray-Hare rank-based test results
#' - \code{type='divergence'} for divergence metrics
#' - \code{type='jackknife'} for jackknife resampling results
#' Use \code{metadata(obj)} to access reproducibility metadata.
#'
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#'
setClass("TSENATAnalysis", slots = list(se = "SummarizedExperiment", config = "list",
    diversity_results = "list", lm_results = "list", pairwise_results = "list", rank_test_results = "list",
    jackknife_results = "list", divergence_results = "list", plots = "list", metadata = "list"),
    prototype = list(config = list(), diversity_results = list(), lm_results = list(),
        pairwise_results = list(), rank_test_results = list(), jackknife_results = list(),
        divergence_results = list(), plots = list(), metadata = list(function_calls = character(0),
            function_timestamps = character(0))), validity = function(object) {
        # Check @se is SummarizedExperiment
        if (!inherits(object@se, "SummarizedExperiment")) {
            return("@se must be a SummarizedExperiment object")
        }

        # Check SE has data
        if (nrow(object@se) == 0 || ncol(object@se) == 0) {
            return("@se has zero dimensions (no genes or samples)")
        }

        # Check @config is list-like (list, TSENATConfig, or other list-based
        # structure)
        if (!is.list(object@config)) {
            return("@config must be a list or list-based config object")
        }

        # Check all results slots are lists
        if (!is.list(object@diversity_results)) {
            return("@diversity_results must be a list")
        }
        if (!is.list(object@lm_results)) {
            return("@lm_results must be a list")
        }
        if (!is.list(object@rank_test_results)) {
            return("@rank_test_results must be a list")
        }
        if (!is.list(object@jackknife_results)) {
            return("@jackknife_results must be a list")
        }
        if (!is.list(object@divergence_results)) {
            return("@divergence_results must be a list")
        }
        if (!is.list(object@plots)) {
            return("@plots must be a list")
        }
        if (!is.list(object@metadata)) {
            return("@metadata must be a list")
        }

        # Check colData has required columns if sample metadata expected
        cdata <- colData(object@se)
        if (!is.null(cdata) && ncol(cdata) > 0) {
            if (!"sample_id" %in% colnames(cdata)) {
                return("colData missing 'sample_id' column required for analysis")
            }
        }

        # Check rowData has gene identifiers if results computed
        rdata <- rowData(object@se)
        if (!is.null(rdata) && ncol(rdata) > 0) {
            if (!"gene_id" %in% colnames(rdata) && !"transcript_id" %in% colnames(rdata)) {
                return("rowData missing 'gene_id' or 'transcript_id' column")
            }
        }

        # NEW: Validate config parameters against SE metadata
        if (length(object@config) > 0) {
            cdata <- colData(object@se)

            # Validate condition_col if specified
            if ("condition_col" %in% names(object@config)) {
                col <- object@config$condition_col
                if (!is.null(col) && !col %in% colnames(cdata)) {
                  return(sprintf("@config$condition_col '%s' not found in colData. Available: %s",
                    col, paste(colnames(cdata), collapse = ", ")))
                }
            }

            # Validate subject_col if specified and paired=TRUE
            if ("subject_col" %in% names(object@config)) {
                if (isTRUE(object@config$paired)) {
                  col <- object@config$subject_col
                  if (!is.null(col) && !col %in% colnames(cdata)) {
                    return(sprintf("@config$subject_col '%s' not found in colData. Available: %s",
                      col, paste(colnames(cdata), collapse = ", ")))
                  }
                }
            }
        }

        TRUE
    })

#' Subset TSENATAnalysis Objects
#'
#' Extract a subset of genes and/or samples from a TSENATAnalysis object.
#' Maintains consistency across the underlying SummarizedExperiment and
#' all computed results (diversity, LM, jackknife, divergence).
#'
#' @param x TSENATAnalysis object
#' @param i Integer, logical, or character vector of rows (genes/transcripts)
#'   to retain. If missing, all rows are retained.
#' @param j Integer, logical, or character vector of columns (samples)
#'   to retain. If missing, all columns are retained.
#' @param drop Logical. Currently ignored (included for S4 compatibility).
#'   Always returns TSENATAnalysis (never drops to SE or vector).
#'
#' @details
#' Subsetting preserves all analysis metadata and results while maintaining
#' consistency:
#' - The underlying SummarizedExperiment is subset to the specified
#' genes/samples
#' - Diversity and jackknife results are subset to match sample selection
#' - LM results are recalculated or removed if sample structure changes
#' - Divergence results are subset accordingly
#' - Analysis configuration is preserved
#'
#' @return A new TSENATAnalysis object containing only the specified genes
#' and samples
#'
#' @examples
#' # Create a minimal TSENATAnalysis object
#' library(SummarizedExperiment)
#' counts <- matrix(rpois(200, 10), nrow = 20, ncol = 10)
#' rownames(counts) <- paste0('TX_', 1:20)
#' colnames(counts) <- paste0('S', 1:10)
#' se <- SummarizedExperiment(
#'   assays = list(counts = counts),
#'   rowData = data.frame(gene_id = rep(paste0('G', 1:2), each = 10),
#'                        row.names = rownames(counts)),
#'   colData = data.frame(sample_id = colnames(counts),
#'                        condition = rep(c('A', 'B'), 5),
#'                        row.names = colnames(counts))
#' )
#' analysis <- TSENATAnalysis(se)
#'
#' # Subset to first 10 genes and first 5 samples
#' analysis_subset <- analysis[1:10, 1:5]
#'
#' # Subset by gene name
#' analysis_subset2 <- analysis[paste0('TX_', 1:5), ]
#'

# NOTE: Subsetting method '[' defined in methods-TSENATAnalysis.R
