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
#'     \item Results from diversity, linear model, jackknife, and divergence analyses
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
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
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
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
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

  # Ensure sample_id column exists in colData (required by validator)
  if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
    SummarizedExperiment::colData(se)$sample_id <- colnames(se)
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

# Accessor Methods for TSENATAnalysis Objects
# Standard methods for extracting results and metadata from TSENATAnalysis
# objects. Following Bioconductor conventions (DESeq2, edgeR).
#

# ============================================================================
# DIVERSITY ACCESSOR
# ============================================================================

#' Extract diversity results
#'
#' @param object \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value to extract (e.g., 1.0, 2.0).
#'   If NULL (default), returns list of all q-values.
#'
#' @return SummarizedExperiment or list of SummarizedExperiment objects
#'   containing diversity values keyed by q-value.
#'
#' @details
#' Results stored in @diversity_results with names like "q_0.5", "q_1.0", etc.
#' Use \code{diversity(analysis)} to get all results as a list, or
#' \code{diversity(analysis, q=1.0)} for a specific q-value.
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{divergence}}, \code{\link{jackKnife}},
#' \code{\link{lmResults}}, \code{\link{se}}, \code{\link[S4Vectors]{metadata}}
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv",
#'   package = "TSENAT"), header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz",
#'   package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file,
#'   metadata = metadata_df, tpm = salmon_tpm,
#'   effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' analysis <- calculate_diversity_s4(analysis, q = 1, verbose = FALSE)
#' diversity_results <- diversity(analysis)
#'
#' @export
setGeneric("diversity", function(object, q = NULL) {
  standardGeneric("diversity")
})

#' @rdname diversity
#' @export
setMethod("diversity", "TSENATAnalysis", function(object, q = NULL) {
  if (length(object@diversity_results) == 0) {
    warning("No diversity results found. Run calculate_diversity_s4() first.")
    return(NULL)
  }

  if (is.null(q)) {
    # Return all results
    return(object@diversity_results)
  }

  # Format q-value key - try multiple precision levels for robustness
  supported_decimals <- c(1, 2, 3)
  q_key <- NULL
  
  for (decimals in supported_decimals) {
    candidate_key <- paste0("q_", formatC(q, format = "f", digits = decimals))
    if (candidate_key %in% names(object@diversity_results)) {
      q_key <- candidate_key
      break
    }
  }

  if (is.null(q_key)) {
    # Fallback: check if lazy conversion is needed from combined result
    if (!is.null(object@metadata$diversity_combined) && 
        is.list(object@metadata$diversity_combined) &&
        !is.null(object@metadata$diversity_combined$combined_result)) {
      
      # Perform lazy conversion from combined format
      combined_result <- object@metadata$diversity_combined$combined_result
      
      # Extract columns for this q-value from combined result
      q_cols <- grep(paste0("_q=", gsub("\\.", "\\\\.", as.character(q)), "$"), 
                     colnames(combined_result))
      
      if (length(q_cols) > 0) {
        # Extract per-q data
        result_subset <- combined_result[, q_cols, drop = FALSE]
        
        # Convert to SummarizedExperiment
        assay_matrix <- as.matrix(result_subset[, vapply(result_subset, is.numeric, FUN.VALUE = logical(1))])
        result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
        rownames(result_se) <- rownames(result_subset)
        
        # Apply colData from original SE
        if (!is.null(object@se)) {
          orig_coldata <- SummarizedExperiment::colData(object@se)
          if (!is.null(orig_coldata) && nrow(orig_coldata) == ncol(result_se)) {
            SummarizedExperiment::colData(result_se) <- orig_coldata
          }
        }
        
        # Cache this result for future access
        q_key_to_cache <- paste0("q_", formatC(q, format = "f", digits = 3))
        object@diversity_results[[q_key_to_cache]] <- result_se
        
        return(result_se)
      }
    }
    
    stop("Q-value ", q, " not found in diversity_results.\n",
         "Available q-values: ", paste(names(object@diversity_results), collapse = ", "),
         call. = FALSE)
  }

  object@diversity_results[[q_key]]
})

# ============================================================================
# LM RESULTS ACCESSOR
# ============================================================================

#' Extract linear model results
#'
#' @param object \code{TSENATAnalysis} object.
#' @param component \code{character}. Which result component to extract.
#'   Options: NULL (all), "results", "lm_interaction", "q_interactions",
#'   "divergence_difference", etc.
#'
#' @return List or data.frame depending on component requested.
#'
#' @details
#' The results stored in the \code{@@lm_results} slot contain multiple analysis types.
#' Use \code{lmResults(analysis)} to get all components, or specify component type
#' for targeted extraction.
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{diversity}}, \code{\link{divergence}},
#' \code{\link{jackKnife}}, \code{\link{se}}, \code{\link[S4Vectors]{metadata}}
#'
#' @examples
#' data(readcounts)
#' analysis <- build_analysis_s4(salmon_dataset,
#'   system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT'))
#' analysis <- subset_analysis(analysis, n_genes = 15, n_samples = 6)
#' res <- lmResults(analysis)
#'
#' @export
setGeneric("lmResults", function(object, component = NULL) {
  standardGeneric("lmResults")
})

#' @rdname lmResults
#' @export
setMethod("lmResults", "TSENATAnalysis", function(object, component = NULL) {
  if (length(object@lm_results) == 0) {
    warning("No LM results found. Run calculate_lm_interaction_s4() first.")
    return(NULL)
  }

  if (is.null(component)) {
    # Return all LM results
    return(object@lm_results)
  }

  # Try to extract specific component
  if (component %in% names(object@lm_results)) {
    return(object@lm_results[[component]])
  }

  # If component.results pattern, extract the $results subcomponent
  if (component %in% c("results", "p_value", "effect_size")) {
    # Search all subcomponents
    for (name in names(object@lm_results)) {
      if (is.list(object@lm_results[[name]]) &&
          "results" %in% names(object@lm_results[[name]])) {
        results_df <- object@lm_results[[name]]$results
        if (is.data.frame(results_df) && component %in% colnames(results_df)) {
          return(results_df[[component]])
        }
      }
    }
  }

  stop("Component '", component, "' not found in lm_results.\n",
       "Available: ", paste(names(object@lm_results), collapse = ", "),
       call. = FALSE)
})

#' @rdname lmResults
#' @export
setGeneric("lmResults<-", function(object, value) {
  standardGeneric("lmResults<-")
})

#' @rdname lmResults
#' @param value A list of LM results to assign to the object.
#' @export
setMethod("lmResults<-", "TSENATAnalysis", function(object, value) {
  if (!is.list(value)) {
    stop("lmResults value must be a list", call. = FALSE)
  }
  object@lm_results <- value
  object
})

# ============================================================================
# JACKKNIFE ACCESSOR
# ============================================================================

#' Extract jackknife resampling results
#'
#' @param object \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value for jackknife results (e.g., 1.0).
#'   If NULL (default), returns all q-values.
#'
#' @return Jackknife result object (confidence intervals, resamples, etc.).
#'
#' @details
#' Jackknife results are stored per q-value. Use this to access confidence
#' intervals and diagnostic information from resampling.
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{diversity}}, \code{\link{divergence}},
#' \code{\link{lmResults}}, \code{\link{se}}, \code{\link[S4Vectors]{metadata}}
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' jk_results <- jackKnife(analysis)
#'
#' @export
setGeneric("jackKnife", function(object, q = NULL) {
  standardGeneric("jackKnife")
})

#' @rdname jackKnife
#' @export
setMethod("jackKnife", "TSENATAnalysis", function(object, q = NULL) {
  if (length(object@jackknife_results) == 0) {
    warning("No jackknife results found. Run jackknife_isoform_switching_s4() first.")
    return(NULL)
  }

  if (is.null(q)) {
    # Return all results
    return(object@jackknife_results)
  }

  # Format q-value key - must match storage format used by base .jackknife_isoform_switching()
  # Uses paste0("q_", gsub("\\.", "_", sprintf("%.2f", q))) to store (e.g., "q_0_01", "q_1_00")
  q_key <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))

  if (!(q_key %in% names(object@jackknife_results))) {
    stop("Q-value ", q, " not found in jackknife_results.\n",
         "Available q-values: ", paste(names(object@jackknife_results), collapse = ", "),
         call. = FALSE)
  }

  object@jackknife_results[[q_key]]
})

# ============================================================================
# DIVERGENCE ACCESSOR
# ============================================================================

#' Extract divergence results
#'
#' @param object \code{TSENATAnalysis} object.
#' @param component \code{character}. Component to extract: NULL (all),
#'   "tsallis_divergence", "effect_sizes", etc.
#'
#' @return SummarizedExperiment or data.frame with divergence metrics.
#'
#' @details
#' Divergence results are stored in @divergence_results with component names
#' corresponding to different divergence metrics. Use \code{divergence(analysis)}
#' to retrieve all components or specify a component for targeted extraction.
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{diversity}}, \code{\link{jackKnife}},
#' \code{\link{lmResults}}, \code{\link{se}}, \code{\link[S4Vectors]{metadata}}
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' div_res <- divergence(analysis)
#'
#' @export
setGeneric("divergence", function(object, component = NULL) {
  standardGeneric("divergence")
})

#' @rdname divergence
#' @export
setMethod("divergence", "TSENATAnalysis", function(object, component = NULL) {
  if (length(object@divergence_results) == 0) {
    warning("No divergence results found. Run calculate_divergence_s4() first.")
    return(NULL)
  }

  if (is.null(component)) {
    # Return all results
    return(object@divergence_results)
  }

  if (!(component %in% names(object@divergence_results))) {
    stop("Component '", component, "' not found in divergence_results.\n",
         "Available: ", paste(names(object@divergence_results), collapse = ", "),
         call. = FALSE)
  }

  object@divergence_results[[component]]
})

# ============================================================================
# PLOT ACCESSORS
# ============================================================================

#' Get cached plot
#'
#' @param object \code{TSENATAnalysis} object.
#' @param type \code{character}. Plot type: "q_curve", "lm_interaction",
#'   "divergence", "influence", "volcano", etc.
#'   If NULL, returns all cached plots.
#'
#' @return ggplot object or list of plots.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' all_plots <- getPlot(analysis)
#'
#' @noRd
setGeneric("getPlot", function(object, type = NULL) {
  standardGeneric("getPlot")
})

#' @rdname getPlot

#' @noRd
setMethod("getPlot", "TSENATAnalysis", function(object, type = NULL) {
  if (length(object@plots) == 0) {
    warning("No plots found. Run tsenat() with generate_plots=TRUE.")
    return(NULL)
  }

  if (is.null(type)) {
    # Return all plots
    return(object@plots)
  }

  if (!(type %in% names(object@plots))) {
    warning("Plot type '", type, "' not found.\n",
            "Available: ", paste(names(object@plots), collapse = ", "))
    return(NULL)
  }

  object@plots[[type]]
})

#' Add plot to cache
#'
#' @param object \code{TSENATAnalysis} object.
#' @param type \code{character}. Name/type for this plot.
#' @param plot Object to cache (ggplot, etc.).
#' @param replace \code{logical}. If TRUE, replace existing plot of same type.
#'
#' @return Modified TSENATAnalysis object.
#'
#' @examples
#' # Demonstrates adding a plot to analysis object (requires ggplot2)
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(100, 10), nrow = 10, ncol = 10)),
#'   colData = data.frame(sample_id = paste0("S", 1:10))
#' )
#' analysis <- TSENATAnalysis(se)
#' # analysis <- addPlot(analysis, type = "example", plot = NULL)
#'

#' @noRd
setGeneric("addPlot", function(object, type, plot, replace = FALSE) {
  standardGeneric("addPlot")
})

#' @rdname addPlot

#' @noRd
setMethod("addPlot", "TSENATAnalysis", function(object, type, plot, replace = FALSE) {
  if (!replace && type %in% names(object@plots)) {
    warning("Plot type '", type, "' already exists. Set replace=TRUE to overwrite.",
            call. = FALSE)
    return(object)
  }

  object@plots[[type]] <- plot
  object
})

# ============================================================================
# SHOW METHOD
# ============================================================================

#' Display TSENATAnalysis object
#'
#' @param object \code{TSENATAnalysis} object to display.
#'
#' @return Invisibly returns the \code{TSENATAnalysis} object (called for its
#'   side effect of printing a formatted summary to the console).
#'
#' @details
#' Provides a concise summary of: data dimensions, completed analyses,
#' number of results, and metadata tracking.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' show(analysis)
#'
#' @export
setMethod("show", "TSENATAnalysis", function(object) {
  message("TSENATAnalysis object")
  message("=====================")

  # Show SE info
  message("SummarizedExperiment:")
  message(sprintf("  Genes:   %d", nrow(object@se)))
  message(sprintf("  Samples: %d", ncol(object@se)))

  # Show config
  if (length(object@config) > 0) {
    message("\nConfiguration:")
    for (name in names(object@config)) {
      val <- object@config[[name]]
      if (is.character(val) && length(val) == 1) {
        message(sprintf("  %s: %s", name, val))
      } else if (is.numeric(val) && length(val) <= 3) {
        message(sprintf("  %s: %s", name, paste(val, collapse = ", ")))
      } else {
        message(sprintf("  %s: <%s>", name, class(val)))
      }
    }
  }

  # Show results
  message("\nAnalysis Status:")
  if (length(object@diversity_results) > 0) {
    message(sprintf("  \u2713 Diversity: %d q-value(s)", length(object@diversity_results)))
  }
  if (length(object@lm_results) > 0) {
    message(sprintf("  \u2713 LM results: %s", paste(names(object@lm_results), collapse = ", ")))
  }
  if (length(object@jackknife_results) > 0) {
    message(sprintf("  \u2713 Jackknife: %d q-value(s)", length(object@jackknife_results)))
  }
  if (length(object@divergence_results) > 0) {
    message(sprintf("  \u2713 Divergence: %d component(s)", length(object@divergence_results)))
  }
  if (length(object@plots) > 0) {
    message(sprintf("  \u2713 Plots: %s", paste(names(object@plots), collapse = ", ")))
  }

  # Show metadata
  if (length(object@metadata) > 0 && "function_calls" %in% names(object@metadata)) {
    n_calls <- length(object@metadata$function_calls)
    if (n_calls > 0) {
      message("\nFunction History:")
      message(sprintf("  Calls: %s", paste(object@metadata$function_calls, collapse = " \u2192 ")))
    }
  }

  message("")
})

# ============================================================================
# SUMMARY METHOD
# ============================================================================

#' Detailed summary of TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return Prints detailed summary (invisibly returns object).
#'
#' @details
#' Provides comprehensive analysis summary including dimensions, results
#' counts, validation status, and metadata tracking.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' summary(analysis)
#'
#' @export
setMethod("summary", "TSENATAnalysis", function(object) {
  message("=== TSENAT Analysis Summary ===")

  # Dimensions
  message("DATA:")
  message(sprintf("  Genes:    %6d", nrow(object@se)))
  message(sprintf("  Samples:  %6d", ncol(object@se)))
  message(sprintf("  Assays:   %6d (%s)",
              length(SummarizedExperiment::assays(object@se)),
              paste(SummarizedExperiment::assayNames(object@se), collapse = ", ")))

  # Configuration
  message("\nCONFIGURATION:")
  if (length(object@config) == 0) {
    message("  (None set)")
  } else {
    for (name in names(object@config)) {
      val <- object@config[[name]]
      if (is.character(val)) {
        if (length(val) == 1) {
          message(sprintf("  %s: %s", name, val))
        } else {
          message(sprintf("  %s: <%d values>", name, length(val)))
        }
      } else if (is.numeric(val)) {
        if (length(val) <= 5) {
          message(sprintf("  %s: %s", name, paste(round(val, 2), collapse = ", ")))
        } else {
          message(sprintf("  %s: <%d values>", name, length(val)))
        }
      } else {
        message(sprintf("  %s: <%s>", name, class(val)))
      }
    }
  }

  # Results Summary
  message("\nRESULTS:")

  if (length(object@diversity_results) > 0) {
    q_vals <- gsub("q_", "", names(object@diversity_results))
    message(sprintf("  Diversity:   %d analyses at q = %s",
                length(object@diversity_results),
                paste(q_vals, collapse = ", ")))
  }

  if (length(object@lm_results) > 0) {
    message(sprintf("  LM/Stats:    %d result set(s) (%s)",
                length(object@lm_results),
                paste(names(object@lm_results), collapse = ", ")))

    # Show gene counts if results available
    for (name in names(object@lm_results)) {
      if (is.list(object@lm_results[[name]]) &&
          "results" %in% names(object@lm_results[[name]]) &&
          is.data.frame(object@lm_results[[name]]$results)) {
        n_genes <- nrow(object@lm_results[[name]]$results)
        message(sprintf("    - %s: %d genes", name, n_genes))
      }
    }
  }

  if (length(object@jackknife_results) > 0) {
    q_vals <- gsub("q_", "", names(object@jackknife_results))
    message(sprintf("  Jackknife:   %d analyses at q = %s",
                length(object@jackknife_results),
                paste(q_vals, collapse = ", ")))
  }

  if (length(object@divergence_results) > 0) {
    message(sprintf("  Divergence:  %d component(s) (%s)",
                length(object@divergence_results),
                paste(names(object@divergence_results), collapse = ", ")))
  }

  if (length(object@plots) > 0) {
    message(sprintf("  Plots:       %d cached (%s)",
                length(object@plots),
                paste(names(object@plots), collapse = ", ")))
  }

  # Metadata
  message("\nMETADATA:")
  if ("created_at" %in% names(object@metadata)) {
    message(sprintf("  Created: %s", format(object@metadata$created_at, "%Y-%m-%d %H:%M:%S")))
  }
  if ("package_version" %in% names(object@metadata)) {
    message(sprintf("  Package: TSENAT %s", object@metadata$package_version))
  }
  if ("function_calls" %in% names(object@metadata) && length(object@metadata$function_calls) > 0) {
    message(sprintf("  Workflow: %s", paste(object@metadata$function_calls, collapse = " -> ")))
  }

  message("")

  invisible(object)
})

# ============================================================================
# CONFIGURATION ACCESSORS (GAP 3 FIX)
# ============================================================================

#' @rdname divResults
setGeneric("getConfig", function(object) {
  standardGeneric("getConfig")
})

#' @rdname divResults
#' @export
setMethod("getConfig", "TSENATAnalysis", function(object) {
  object@config
})

#' Replace analysis configuration
#'
#' @param object \code{TSENATAnalysis} object.
#' @param value \code{list}. New configuration. Should contain q_values, 
#'   column name specifications, and/or other parameters.
#'
#' @return Modified \code{TSENATAnalysis} object with updated @config.
#'
#' @details
#' Provides type-safe replacement of @config slot. Typically called once
#' at the start of an analysis via \code{tsenat_config()} rather than directly.
#'
#' Replace entire configuration in TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#' @param value List or \code{TSENATConfig} object containing configuration settings.
#'
#' @return Updated \code{TSENATAnalysis} object with replaced configuration.
#'
#' @details
#' Replaces the entire configuration of a TSENATAnalysis object. This method
#' is useful when you need to apply a new set of configuration parameters to
#' an existing analysis object. All previous configuration values are replaced
#' with those in the new value object.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' # Build and subset analysis
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' 
#' # Replace configuration with new settings
#' new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5), seed = 42)
#' analysis <- setConfig(analysis, new_config)
#' 
#' # Verify the new configuration was applied
#' current_config <- getConfig(analysis)
#' print(current_config$q_values)  # Shows c(0.5, 1.0, 1.5)
#'
#' @export
setGeneric("setConfig", function(object, value) {
  standardGeneric("setConfig")
})

#' @rdname setConfig
#' @aliases setConfig,TSENATAnalysis-method
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' # Build and subset analysis
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' 
#' # Use setConfig via the method (called by setConfig generic)
#' new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5, 2.0), seed = 123)
#' analysis <- setConfig(analysis, new_config)
#' 
#' # Verify and display
#' summary(analysis)
#' }
setMethod("setConfig", "TSENATAnalysis", function(object, value) {
  # Convert TSENATConfig S4 object to list if needed
  if (inherits(value, "TSENATConfig")) {
    value <- unclass(value)
  }
  
  if (!is.list(value)) {
    stop("Configuration must be a list or TSENATConfig object", call. = FALSE)
  }
  object@config <- value
  validObject(object)
  object
})

#' Set a single configuration value in TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#' @param key Character. Name of the config item to set.
#' @param value Value to set for the config item.
#'
#' @return Updated \code{TSENATAnalysis} object with the new config value.
#'
#' @details
#' Convenience method for setting a single configuration value without
#' needing to retrieve, merge, and set the entire config list. This preserves
#' all other configuration values while updating only the specified key.
#'
#' Unlike \code{\link{setConfig}}, which replaces the entire configuration,
#' \code{setConfigValue} performs a targeted update. It retrieves the current
#' config, updates one key-value pair, and stores the modified config back.
#'
#' @seealso
#' \code{\link{setConfig}} for replacing entire configuration,
#' \code{\link{getConfig}} for retrieving configuration,
#' \code{\link{tsenat_config}} for creating configuration objects
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' # Build and subset analysis
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' 
#' # Update a single configuration value while preserving others
#' analysis <- setConfigValue(analysis, "q_values", c(0.5, 1.0, 1.5))
#' 
#' # Verify the update
#' config <- getConfig(analysis)
#' print(config$q_values)  # Shows c(0.5, 1.0, 1.5)
#'
#' @export
setGeneric("setConfigValue", function(object, key, value) {
  standardGeneric("setConfigValue")
})

#' @rdname setConfigValue
#' @aliases setConfigValue,TSENATAnalysis-method
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' # Build and subset analysis with initial configuration
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' 
#' # Use setConfigValue to update single configuration values
#' # This preserves all other config values
#' analysis <- setConfigValue(analysis, "q_values", c(0.5, 1.5, 2.0))
#' analysis <- setConfigValue(analysis, "seed", 456)
#' 
#' # Verify the updates
#' summary(analysis)
#' }
setMethod("setConfigValue", "TSENATAnalysis", function(object, key, value) {
  config <- getConfig(object)
  if (is.null(config)) {
    config <- list()
  }
  config[[key]] <- value
  setConfig(object, config)
})

#' Extract SummarizedExperiment from TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return The \code{SummarizedExperiment} containing transcript/gene counts.
#'
#' @details
#' Provides type-safe accessor for the embedded \code{SummarizedExperiment}.
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{diversity}}, \code{\link{divergence}},
#' \code{\link{jackKnife}}, \code{\link{lmResults}}, \code{\link[S4Vectors]{metadata}}
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' @export
setGeneric("se", function(object) {
  standardGeneric("se")
})

#' @noRd
if (!isGeneric("metadata")) {
  setGeneric("metadata", function(x, key = NULL) {
    standardGeneric("metadata")
  })
}

#' @rdname se
#' @export
setMethod("se", "TSENATAnalysis", function(object) {
  object@se
})

#' Extract metadata from TSENATAnalysis
#'
#' @param x \code{TSENATAnalysis} object.
#' @param key \code{character} (optional). Specific metadata key to extract.
#'   If NULL, returns entire metadata list.
#'
#' @return The metadata list, or a specific metadata element if key is provided.
#'
#' @details
#' Provides type-safe accessor for analysis metadata (timestamps, function calls, 
#' intermediate results, etc.).
#'
#' @seealso
#' Other TSENATAnalysis accessors: \code{\link{diversity}}, \code{\link{divergence}},
#' \code{\link{jackKnife}}, \code{\link{lmResults}}, \code{\link{se}}
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t")
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' analysis <- build_analysis_s4(salmon_dataset, gff3_file, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' @noRd
#' @rdname metadata

setMethod("metadata", "TSENATAnalysis", function(x, key = NULL) {
  if (is.null(key)) {
    return(x@metadata)
  }
  
  if (key %in% names(x@metadata)) {
    return(x@metadata[[key]])
  }
  
  NULL
})

#' Test rank-based method assumptions
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results.
#' @param q \code{numeric}. Q-value for diversity data to test.
#' @param checks \code{character}. Which assumptions to check.
#' @param alpha \code{numeric}. Significance level (default: 0.05).
#' @param ... Additional arguments.
#'
#' @return \code{TSENATAnalysis} object with stored assumption test results.
#'
#' @details
#' Tests whether rank-based analysis methods are appropriate for the data
#' by evaluating exchangeability, monotonicity, and consistency assumptions.
#'
#' @export
setGeneric("test_rankbased_assumptions_s4", function(analysis, q = NULL,
                                                      checks = c("exchangeability",
                                                                 "monotonicity",
                                                                 "consistency"),
                                                      alpha = 0.05, ...) {
  standardGeneric("test_rankbased_assumptions_s4")
})
