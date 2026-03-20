#' TSENATAnalysis S4 Class
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
#'   (genes × samples) with assays and colData.
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
#' @section Accessor Methods:
#'   \describe{
#'     \item{\code{diversity(object, q=NULL)}}{Extract diversity results for q-value}
#'     \item{\code{lmResults(object, component=NULL)}}{Extract LM results}
#'     \item{\code{jackKnife(object, q=NULL)}}{Extract jackknife results}
#'     \item{\code{divergence(object)}}{Extract divergence results}
#'     \item{\code{getPlot(object, type=NULL)}}{Retrieve cached plot}
#'     \item{\code{addPlot(object, type, plot)}}{Add/cache a new plot}
#'     \item{\code{show(object)}}{Display object summary}
#'     \item{\code{summary(object)}}{Get detailed analysis summary}
#'   }
#'
#' @section Validation:
#'   Validity is checked at object construction. Ensures @se is a
#'   SummarizedExperiment and all slots are correct types.
#'
#' @examples
#' \dontrun{
#'   # Create from SummarizedExperiment
#'   analysis <- TSENATAnalysis(se)
#'
#'   # Or configure with metadata first
#'   analysis <- tsenat_config(se, 
#'     q_values = seq(0.5, 2, 0.1),
#'     sample_type_col = "condition",
#'     subject_col = "patient_id"
#'   )
#'
#'   # Access results after analysis
#'   div_results <- diversity(analysis, q = 1.0)
#'   lm_df <- lmResults(analysis, component = "results")
#'   summary(analysis)
#' }
#'
#' @name TSENATAnalysis-class
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
setClass(
  "TSENATAnalysis",
  slots = list(
    se = "SummarizedExperiment",
    config = "ANY",
    diversity_results = "list",
    lm_results = "list",
    jackknife_results = "list",
    divergence_results = "list",
    plots = "list",
    metadata = "list"
  ),
  validity = function(object) {
    # Check @se is SummarizedExperiment
    if (!inherits(object@se, "SummarizedExperiment")) {
      return("@se must be a SummarizedExperiment object")
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

    TRUE
  }
)

#' Constructor for TSENATAnalysis objects
#'
#' Creates a new TSENATAnalysis object with a SummarizedExperiment base
#' and optional initial configuration.
#'
#' @param se \code{SummarizedExperiment}. The base expression data object.
#' @param config \code{list}. Optional initial configuration (usually set
#'   via \code{tsenat_config()} instead).
#'
#' @return A new \code{TSENATAnalysis} object.
#'
#' @details
#' The constructor initializes all slots with empty lists except @se,
#' which must be provided. The @metadata slot automatically records:
#' - creation timestamp
#' - TSENAT package version
#' - initial function call
#'
#' @examples
#' \dontrun{
#'   analysis <- TSENATAnalysis(se)
#'   show(analysis)  # Display empty initialized object
#' }
#'
#' @export
TSENATAnalysis <- function(se, config = list()) {
  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object", call. = FALSE)
  }

  # Create new object with all slots initialized
  new(
    "TSENATAnalysis",
    se = se,
    config = if (length(config) > 0) config else list(),
    diversity_results = list(),
    lm_results = list(),
    jackknife_results = list(),
    divergence_results = list(),
    plots = list(),
    metadata = list(
      created_at = Sys.time(),
      package_version = as.character(utils::packageVersion("TSENAT")),
      function_calls = character()
    )
  )
}
