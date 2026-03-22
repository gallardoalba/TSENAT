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
#' @examples
#' \dontrun{
#'   analysis <- tsenat(se, methods = "diversity")
#'   all_div <- diversity(analysis)        # List of all q-values
#'   div_at_q1 <- diversity(analysis, q = 1.0)  # Specific q-value
#' }
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
#' @examples
#' \dontrun{
#'   lm_df <- lmResults(analysis, component = "results")
#'   all_lm <- lmResults(analysis)
#' }
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
#' @examples
#' \dontrun{
#'   jk_q1 <- jackKnife(analysis, q = 1.0)
#'   ci <- jk_q1$confidence_intervals
#' }
#'
#' @export
setGeneric("jackKnife", function(object, q = NULL) {
  standardGeneric("jackKnife")
})

#' @rdname jackKnife
#' @export
setMethod("jackKnife", "TSENATAnalysis", function(object, q = NULL) {
  if (length(object@jackknife_results) == 0) {
    warning("No jackknife results found. Run jackknife_tsallis_entropy_s4() first.")
    return(NULL)
  }

  if (is.null(q)) {
    # Return all results
    return(object@jackknife_results)
  }

  # Format q-value key - ensure numeric precision
  q_key <- paste0("q_", formatC(q, format = "f", digits = 1))

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
#' @examples
#' \dontrun{
#'   div <- divergence(analysis)
#'   effect_sizes <- divergence(analysis, component = "effect_sizes")
#' }
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
#' \dontrun{
#'   p_qcurve <- getPlot(analysis, type = "q_curve")
#'   all_plots <- getPlot(analysis)
#' }
#'
#' @keywords internal
#' @noRd
setGeneric("getPlot", function(object, type = NULL) {
  standardGeneric("getPlot")
})

#' @rdname getPlot
#' @keywords internal
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
#' \dontrun{
#'   p <- ggplot(...) + ...
#'   analysis <- addPlot(analysis, type = "custom_plot", plot = p)
#' }
#'
#' @keywords internal
#' @noRd
setGeneric("addPlot", function(object, type, plot, replace = FALSE) {
  standardGeneric("addPlot")
})

#' @rdname addPlot
#' @keywords internal
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
#' @details
#' Provides a concise summary of: data dimensions, completed analyses,
#' number of results, and metadata tracking.
#'
#' @examples
#' \dontrun{
#'   analysis <- TSENATAnalysis(se)
#'   show(analysis)
#' }
#'
#' @export
setMethod("show", "TSENATAnalysis", function(object) {
  cat("TSENATAnalysis object\n")
  cat("=====================\n\n")

  # Show SE info
  cat("SummarizedExperiment:\n")
  cat("  Genes:  ", nrow(object@se), "\n")
  cat("  Samples:", ncol(object@se), "\n")

  # Show config
  if (length(object@config) > 0) {
    cat("\nConfiguration:\n")
    for (name in names(object@config)) {
      val <- object@config[[name]]
      if (is.character(val) && length(val) == 1) {
        cat("  ", name, ": ", val, "\n", sep = "")
      } else if (is.numeric(val) && length(val) <= 3) {
        cat("  ", name, ": ", paste(val, collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", name, ": <", class(val), ">\n", sep = "")
      }
    }
  }

  # Show results
  cat("\nAnalysis Status:\n")
  if (length(object@diversity_results) > 0) {
    cat("  \u2713 Diversity: ", length(object@diversity_results), " q-value(s)\n", sep = "")
  }
  if (length(object@lm_results) > 0) {
    cat("  \u2713 LM results: ", paste(names(object@lm_results), collapse = ", "), "\n", sep = "")
  }
  if (length(object@jackknife_results) > 0) {
    cat("  \u2713 Jackknife: ", length(object@jackknife_results), " q-value(s)\n", sep = "")
  }
  if (length(object@divergence_results) > 0) {
    cat("  \u2713 Divergence: ", length(object@divergence_results), " component(s)\n", sep = "")
  }
  if (length(object@plots) > 0) {
    cat("  \u2713 Plots: ", paste(names(object@plots), collapse = ", "), "\n", sep = "")
  }

  # Show metadata
  if (length(object@metadata) > 0 && "function_calls" %in% names(object@metadata)) {
    n_calls <- length(object@metadata$function_calls)
    if (n_calls > 0) {
      cat("\nFunction History:\n")
      cat("  Calls: ", paste(object@metadata$function_calls, collapse = " \u2192 "), "\n", sep = "")
    }
  }

  cat("\n")
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
#' \dontrun{
#'   analysis <- tsenat(se, methods = c("diversity", "lm_interaction"))
#'   summary(analysis)
#' }
#'
#' @export
setMethod("summary", "TSENATAnalysis", function(object) {
  cat("=== TSENAT Analysis Summary ===\n\n")

  # Dimensions
  cat("DATA:\n")
  cat(sprintf("  Genes:    %6d\n", nrow(object@se)))
  cat(sprintf("  Samples:  %6d\n", ncol(object@se)))
  cat(sprintf("  Assays:   %6d (%s)\n",
              length(SummarizedExperiment::assays(object@se)),
              paste(SummarizedExperiment::assayNames(object@se), collapse = ", ")))

  # Configuration
  cat("\nCONFIGURATION:\n")
  if (length(object@config) == 0) {
    cat("  (None set)\n")
  } else {
    for (name in names(object@config)) {
      val <- object@config[[name]]
      if (is.character(val)) {
        if (length(val) == 1) {
          cat(sprintf("  %s: %s\n", name, val))
        } else {
          cat(sprintf("  %s: <%d values>\n", name, length(val)))
        }
      } else if (is.numeric(val)) {
        if (length(val) <= 5) {
          cat(sprintf("  %s: %s\n", name, paste(round(val, 2), collapse = ", ")))
        } else {
          cat(sprintf("  %s: <%d values>\n", name, length(val)))
        }
      } else {
        cat(sprintf("  %s: <%s>\n", name, class(val)))
      }
    }
  }

  # Results Summary
  cat("\nRESULTS:\n")

  if (length(object@diversity_results) > 0) {
    q_vals <- gsub("q_", "", names(object@diversity_results))
    cat(sprintf("  Diversity:   %d analyses at q = %s\n",
                length(object@diversity_results),
                paste(q_vals, collapse = ", ")))
  }

  if (length(object@lm_results) > 0) {
    cat(sprintf("  LM/Stats:    %d result set(s) (%s)\n",
                length(object@lm_results),
                paste(names(object@lm_results), collapse = ", ")))

    # Show gene counts if results available
    for (name in names(object@lm_results)) {
      if (is.list(object@lm_results[[name]]) &&
          "results" %in% names(object@lm_results[[name]]) &&
          is.data.frame(object@lm_results[[name]]$results)) {
        n_genes <- nrow(object@lm_results[[name]]$results)
        cat(sprintf("    - %s: %d genes\n", name, n_genes))
      }
    }
  }

  if (length(object@jackknife_results) > 0) {
    q_vals <- gsub("q_", "", names(object@jackknife_results))
    cat(sprintf("  Jackknife:   %d analyses at q = %s\n",
                length(object@jackknife_results),
                paste(q_vals, collapse = ", ")))
  }

  if (length(object@divergence_results) > 0) {
    cat(sprintf("  Divergence:  %d component(s) (%s)\n",
                length(object@divergence_results),
                paste(names(object@divergence_results), collapse = ", ")))
  }

  if (length(object@plots) > 0) {
    cat(sprintf("  Plots:       %d cached (%s)\n",
                length(object@plots),
                paste(names(object@plots), collapse = ", ")))
  }

  # Metadata
  cat("\nMETADATA:\n")
  if ("created_at" %in% names(object@metadata)) {
    cat(sprintf("  Created: %s\n", format(object@metadata$created_at, "%Y-%m-%d %H:%M:%S")))
  }
  if ("package_version" %in% names(object@metadata)) {
    cat(sprintf("  Package: TSENAT %s\n", object@metadata$package_version))
  }
  if ("function_calls" %in% names(object@metadata) && length(object@metadata$function_calls) > 0) {
    cat(sprintf("  Workflow: %s\n", paste(object@metadata$function_calls, collapse = " -> ")))
  }

  cat("\n")

  invisible(object)
})

# ============================================================================
# CONFIGURATION ACCESSORS (GAP 3 FIX)
# ============================================================================

#' Extract analysis configuration
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return List containing all configuration parameters (q-values, column names, etc.)
#'
#' @details
#' Provides type-safe access to @config slot. After calling wrapper functions like
#' \code{calculate_diversity_s4()}, the configuration is updated with
#' \code{last_diversity_run}, showing which parameters were actually used.
#'
#' @examples
#' \dontrun{
#'   config <- getConfig(analysis)
#'   actual_params <- config$last_diversity_run$parameters_used
#'   # Shows: norm, verbose, bootstrap, nthreads, etc.
#' }
#'
#' @export
setGeneric("getConfig", function(object) {
  standardGeneric("getConfig")
})

#' @rdname getConfig
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
#' @examples
#' \dontrun{
#'   new_config <- list(q_values = seq(0.5, 2, 0.1), nthreads = 4)
#'   analysis <- setConfig(analysis, new_config)
#' }
#'
#' @export
setGeneric("setConfig", function(object, value) {
  standardGeneric("setConfig")
})

#' @rdname setConfig
#' @export
setMethod("setConfig", "TSENATAnalysis", function(object, value) {
  if (!is.list(value)) {
    stop("Configuration must be a list", call. = FALSE)
  }
  object@config <- value
  validObject(object)
  object
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
