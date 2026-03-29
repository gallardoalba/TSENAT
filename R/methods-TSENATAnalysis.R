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
#' @rdname divResults
setMethod("addPlot", "TSENATAnalysis", function(object, type, plot) {
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
        cat("  Function calls:", paste(object@metadata$function_calls, collapse = " → "), "\n")
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
