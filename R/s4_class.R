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
#' @return An S4 object of class \code{TSENATAnalysis} containing:
#'   \itemize{
#'     \item Raw data (SummarizedExperiment)
#'     \item Analysis configuration
#'     \item Results from diversity,  linear model,  jackknife,  and 
#' divergence analyses
#'     \item Cached visualization objects
#'     \item Reproducibility metadata and function history
#'   }
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
#'   results. Each name corresponds to a q-value (e.g., 'q_0.5', 'q_1.0').
#'   Values are SummarizedExperiment objects or data.frames containing entropy
#'   values for each gene at that q-value.
#'
#' @slot lm_results \code{list}. Complex results from linear model and
#'   statistical testing. Top-level names identify analysis type:
#'   \describe{
#'     \item{\code{lm_interaction}}{LM/GAM/GEE model results (list with
#'           \code{$results} data.frame, \code{$models} list, etc.)}
#'     \item{\code{rank_test}}{Friedman/rank-based test results}
#'     \item{\code{divergence_difference}}{Differential divergence comparison}
#'   }
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
#'     \item{\code{results(object, type, ...)}}{Unified interface to extract all analysis results (diversity, divergence, lm, jackknife, rank_test, effect_sizes_divergence, switching_tables)}
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
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @name TSENATAnalysis-class
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
#' 
NULL

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
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @export
TSENATAnalysis <- function(se, config = list()) {
    # Validate input
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    # Convert TSENATConfig S4 object to list if needed
    if (inherits(config, "TSENATConfig")) {
        config <- unclass(config)
    }

    # NEW: Validate config parameters early (before object creation)
    if (length(config) > 0 && !is.null(config)) {
        cdata <- SummarizedExperiment::colData(se)

        # Validate condition_col
        if ("condition_col" %in% names(config)) {
            col <- config$condition_col
            if (!is.null(col) && !col %in% colnames(cdata)) {
                stop(sprintf("Invalid condition_col '%s': column not found in colData.\\nAvailable columns: %s",
                  col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
            }
        }

        # Validate subject_col if paired
        if ("paired" %in% names(config) && isTRUE(config$paired)) {
            if ("subject_col" %in% names(config)) {
                col <- config$subject_col
                if (!is.null(col) && !col %in% colnames(cdata)) {
                  stop(sprintf("Invalid subject_col '%s': column not found in colData.\\nAvailable columns: %s",
                    col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
                }
            }
        }
    }

    # Ensure sample_id column exists in colData (required by validator)
    if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
        SummarizedExperiment::colData(se)$sample_id <- colnames(se)
    }

    # Create new object with all slots initialized
    new("TSENATAnalysis", se = se, config = if (length(config) > 0)
        config else list(), diversity_results = list(), lm_results = list(), jackknife_results = list(),
        divergence_results = list(), plots = list(), metadata = list(created_at = Sys.time(),
            package_version = as.character(utils::packageVersion("TSENAT")), function_calls = character()))
}
