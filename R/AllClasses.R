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
#'   columns, etc.). Set once via \code{tsenat_config()} and used by all
#'   downstream wrapper functions.
#'
#' @slot diversity_results \code{list}. Named list of diversity calculation
#'   results. Each name corresponds to a q-value (e.g., "q_0.5", "q_1.0").
#'   Values are SummarizedExperiment objects or data.frames containing entropy
#'   values for each gene at that q-value.
#'
#' @slot lm_results \code{list}. Complex results from linear model and
#'   statistical testing. Top-level names identify analysis type:
#'   \describe{
#'     \item{\code{lm_interaction}}{LM/GAM/GEE model results (list with
#'           \code{$results} data.frame, \code{$models} list, etc.)}
#'     \item{\code{q_interactions}}{Friedman/rank-based test results}
#'     \item{\code{divergence_difference}}{Differential divergence comparison}
#'   }
#'
#' @slot jackknife_results \code{list}. Resampling-based confidence intervals.
#'   Names correspond to q-values (e.g., "q_0.5", "q_1.0"). Values are
#'   jackknife result objects containing resamples, CI bounds, and diagnostics.
#'
#' @slot divergence_results \code{list}. Divergence metric calculations.
#'   Typically contains:
#'   \describe{
#'     \item{\code{tsallis_divergence}}{SummarizedExperiment with divergence values}
#'     \item{\code{effect_sizes}}{data.frame with Cohen's d, etc.}
#'   }
#'
#' @slot plots \code{list}. Cached visualization objects (ggplot). Names
#'   identify plot type (e.g., "q_curve", "lm_interaction", "influence").
#'   Populated by \code{tsenat()} if \code{generate_plots=TRUE}.
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
#' Access results via accessor methods (recommended):
#' \code{divResults(obj, q)} for diversity, \code{lmRes(obj)} for models,
#' \code{jkResults(obj, q)} for jackknife, \code{divRes(obj)} for divergence,
#' \code{getMeta(obj)} for metadata.
#'
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#'
setClass("TSENATAnalysis",
    slots = list(
        se = "SummarizedExperiment",
        config = "list",
        diversity_results = "list",
        lm_results = "list",
        jackknife_results = "list",
        divergence_results = "list",
        plots = "list",
        metadata = "list"
    ),
    prototype = list(
        config = list(),
        diversity_results = list(),
        lm_results = list(),
        jackknife_results = list(),
        divergence_results = list(),
        plots = list(),
        metadata = list(function_calls = character(0),
                        function_timestamps = character(0))
    ),
    validity = function(object) {
        # Check @se is SummarizedExperiment
        if (!inherits(object@se, "SummarizedExperiment")) {
            return("@se must be a SummarizedExperiment object")
        }

        # Check SE has data
        if (nrow(object@se) == 0 || ncol(object@se) == 0) {
            return("@se has zero dimensions (no genes or samples)")
        }

        # Check @config is list-like (list, TSENATConfig, or other list-based structure)
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

        TRUE
    }
)
