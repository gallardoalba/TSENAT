#' TSENATAnalysis S4 Methods
#'
#' Implementation of accessor and utility methods for TSENATAnalysis objects.
#' These methods follow Bioconductor best practices for accessing object slots.
#'
#' @return
#' Methods return different components of TSENATAnalysis:
#' - getSE: SummarizedExperiment object with count data
#' - divResults: data.frame or list of diversity values per gene
#' - lmRes: list of linear model fitting results
#' - jkResults: list of jackknife diagnostics
#' - divRes: list of divergence analysis results
#' - getMeta: list or atomic value of metadata
#' - getConfig: list of analysis configuration
#' - getPlot: ggplot object or NULL
#' - addPlot: invisible(object) (adds plot to cache)
#'
#' @name methods-TSENATAnalysis
#' @keywords internal
NULL

#' Access SummarizedExperiment data
#' @param object TSENATAnalysis object
#' @rdname divResults
setMethod("getSE", "TSENATAnalysis", function(object) {
    object@se
})

#' Access diversity results
#' @param object TSENATAnalysis object
#' @param q numeric or character. Q-value to extract (e.g., "q_0.5").
#'   If NULL, returns all diversity results.
#' @rdname divResults
setMethod("divResults", "TSENATAnalysis", function(object, q = NULL) {
    if (is.null(q)) {
        return(object@diversity_results)
    }
    q_name <- if (is.numeric(q)) paste0("q_", q) else q
    object@diversity_results[[q_name]]
})

#' Access linear model results
#' @param object TSENATAnalysis object
#' @param component character. Component name ("lm_interaction", "q_interactions", etc.).
#'   If NULL, returns all LM results.
#' @rdname divResults
setMethod("lmRes", "TSENATAnalysis", function(object, component = NULL) {
    if (is.null(component)) {
        return(object@lm_results)
    }
    object@lm_results[[component]]
})

#' Access jackknife results
#' @param object TSENATAnalysis object
#' @param q numeric or character. Q-value to extract (e.g., "q_0.5").
#'   If NULL, returns all jackknife results.
#' @rdname divResults
setMethod("jkResults", "TSENATAnalysis", function(object, q = NULL) {
    if (is.null(q)) {
        return(object@jackknife_results)
    }
    q_name <- if (is.numeric(q)) paste0("q_", q) else q
    object@jackknife_results[[q_name]]
})

#' Access divergence results
#' @param object TSENATAnalysis object
#' @param component character. Component name ("tsallis_divergence", "effect_sizes", etc.).
#'   If NULL, returns all divergence results.
#' @rdname divResults
setMethod("divRes", "TSENATAnalysis", function(object, component = NULL) {
    if (is.null(component)) {
        return(object@divergence_results)
    }
    object@divergence_results[[component]]
})

#' Access metadata
#' @param object TSENATAnalysis object
#' @param key character. Metadata key to extract (e.g., "function_calls", "effect_sizes_divergence").
#'   If NULL, returns all metadata.
#' @rdname divResults
setMethod("getMeta", "TSENATAnalysis", function(object, key = NULL) {
    if (is.null(key)) {
        return(object@metadata)
    }
    object@metadata[[key]]
})

#' Access configuration
#' @param object TSENATAnalysis object
#' @param key character. Config key to extract (e.g., "q_values", "condition_col").
#'   If NULL, returns entire config.
#' @rdname divResults
setMethod("getConfig", "TSENATAnalysis", function(object, key = NULL) {
    if (is.null(key)) {
        return(object@config)
    }
    object@config[[key]]
})

#' Access cached plots
#' @param object TSENATAnalysis object
#' @param type character. Plot type (e.g., "q_curve", "lm_interaction").
#'   If NULL, returns all plots.
#' @rdname divResults
setMethod("getPlot", "TSENATAnalysis", function(object, type = NULL) {
    if (is.null(type)) {
        return(object@plots)
    }
    object@plots[[type]]
})

#' Add or update a cached plot
#' @param object TSENATAnalysis object
#' @param type character. Plot type identifier
#' @param plot ggplot or list. The plot object to cache
#' @param replace logical. If TRUE, replace existing plot of same type.
#'   If FALSE (default), warn if plot already exists and do not overwrite.
#' @rdname divResults
setMethod("addPlot", "TSENATAnalysis", function(object, type, plot, replace = FALSE) {
    if (!replace && type %in% names(object@plots)) {
        warning("Plot type '", type, "' already exists. Set replace=TRUE to overwrite.",
                call. = FALSE)
        return(object)
    }
    object@plots[[type]] <- plot
    object
})

#' Show method for TSENATAnalysis
#' @param object TSENATAnalysis object
setMethod("show", "TSENATAnalysis", function(object) {
    cat("TSENATAnalysis object\n")
    cat("=====================\n\n")
    
    cat("Data: SummarizedExperiment\n")
    cat("  Genes:", nrow(object@se), "\n")
    cat("  Samples:", ncol(object@se), "\n\n")
    
    cat("Analysis results completed:\n")
    if (length(object@diversity_results) > 0) {
        cat("  - Diversity calculations for q =", 
            paste(gsub("q_", "", names(object@diversity_results)), collapse = ", "), "\n")
    }
    if (length(object@lm_results) > 0) {
        cat("  - Linear model results for components:",
            paste(names(object@lm_results), collapse = ", "), "\n")
    }
    if (length(object@jackknife_results) > 0) {
        cat("  - Jackknife results for q =",
            paste(gsub("q_", "", names(object@jackknife_results)), collapse = ", "), "\n")
    }
    if (length(object@divergence_results) > 0) {
        cat("  - Divergence results for components:",
            paste(names(object@divergence_results), collapse = ", "), "\n")
    }
    
    cat("\nMetadata:\n")
    if (length(object@metadata$function_calls) > 0) {
         cat("  Function calls:", paste(object@metadata$function_calls, collapse = " \u2192 "), "\n")
    }
    if (!is.null(object@metadata$created_at)) {
        cat("  Created:", object@metadata$created_at, "\n")
    }
    
    cat("\n")
})

#' Summary method for TSENATAnalysis
#' @param object TSENATAnalysis object
setMethod("summary", "TSENATAnalysis", function(object) {
    message("TSENATAnalysis Summary\n======================\n")
    
    message("Data dimensions:")
    message("  - ", nrow(object@se), " genes x ", ncol(object@se), " samples")
    message("  - Assays: ", paste(names(assays(object@se)), collapse = ", "), "\n")
    
    message("Configuration:")
    if (length(object@config) > 0) {
        for (key in names(object@config)) {
            val <- object@config[[key]]
            if (is.vector(val) && length(val) <= 3) {
                message("  - ", key, ": ", paste(val, collapse = ", "))
            } else {
                message("  - ", key, ": [set]")
            }
        }
    }
    
    invisible(object)
})

#' Retrieve Pseudocount Used in Diversity Calculation
#'
#' Extracts the pseudocount value used in the last `calculate_diversity_s4()` call.
#' When pseudocount="auto", displays the actual computed value.
#'
#' @param object TSENATAnalysis object.
#' @param original logical. If TRUE, returns original parameter (e.g., "auto" string).
#'   If FALSE (default), returns resolved numeric value.
#' @return numeric or character. The pseudocount value used, or NULL if diversity not yet calculated.
#'
#' @examples
#' # analysis <- calculate_diversity_s4(analysis, pseudocount = "auto")
#' # pc <- get_pseudocount(analysis)
#' # pc  # e.g., 0.5267
#'
#' @export
get_pseudocount <- function(object, original = FALSE) {
  if (!is(object, "TSENATAnalysis")) {
    stop("'object' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  # Check for diversity_combined metadata from last calculate_diversity_s4 call
  if (!is.null(object@metadata) && length(object@metadata) > 0) {
    if ("diversity_combined" %in% names(object@metadata)) {
      div_combined <- object@metadata$diversity_combined
      if (!is.null(div_combined) && "computation_params" %in% names(div_combined)) {
        computation_params <- div_combined$computation_params
        if (!is.null(computation_params)) {
          if (original) {
            # Return original parameter (could be "auto" or numeric)
            val <- if ("pseudocount_original" %in% names(computation_params)) {
              computation_params$pseudocount_original
            } else {
              computation_params$pseudocount
            }
            if (!is.null(val) && length(val) > 0) return(val)
          } else {
            # Return resolved numeric value
            val <- computation_params$pseudocount
            if (!is.null(val) && length(val) > 0) return(val)
          }
        }
      }
    }
  }
  
  # Fallback: check last_diversity_run config
  if (!is.null(object@config) && length(object@config) > 0) {
    if ("last_diversity_run" %in% names(object@config)) {
      last_run <- object@config$last_diversity_run
      if (!is.null(last_run) && "parameters_used" %in% names(last_run)) {
        params_used <- last_run$parameters_used
        if (!is.null(params_used)) {
          if (original) {
            val <- if ("pseudocount_original" %in% names(params_used)) {
              params_used$pseudocount_original
            } else {
              params_used$pseudocount
            }
            if (!is.null(val) && length(val) > 0) return(val)
          } else {
            val <- params_used$pseudocount
            if (!is.null(val) && length(val) > 0) return(val)
          }
        }
      }
    }
  }
  
  # No diversity calculation found - explicitly return NULL
  NULL
}
