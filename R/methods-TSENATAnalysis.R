#' TSENATAnalysis S4 Methods
#'
#' Implementation of accessor and utility methods for TSENATAnalysis objects.
#' These methods follow Bioconductor best practices for accessing object slots.
#'
#' @return
#' Methods return different components of TSENATAnalysis:
#' - getSE: SummarizedExperiment object with count data
#' - lmResults: list of linear model fitting results
#' - jackKnife: list of jackknife diagnostics
#' - diversity: list of SummarizedExperiments (one per q-value) or single SE for specific q
#' - divergence: list of divergence analysis results
#' - getMeta: list or atomic value of metadata
#' - getConfig: list of analysis configuration
#' - getPlot: ggplot object or NULL
#' - addPlot: invisible(object) (adds plot to cache)
#'
#' @name methods-TSENATAnalysis
#' @keywords internal
NULL

#' @rdname methods-TSENATAnalysis
setMethod("getSE", "TSENATAnalysis", function(object) {
    object@se
})

#' @rdname methods-TSENATAnalysis
setMethod("getMeta", "TSENATAnalysis", function(object, key = NULL) {
    if (is.null(key)) {
        return(object@metadata)
    }
    object@metadata[[key]]
})

#' @rdname methods-TSENATAnalysis
setMethod("getConfig", "TSENATAnalysis", function(object, key = NULL) {
    if (is.null(key)) {
        return(object@config)
    }
    object@config[[key]]
})

#' @rdname methods-TSENATAnalysis
setMethod("getPlot", "TSENATAnalysis", function(object, type = NULL) {
    if (is.null(type)) {
        return(object@plots)
    }
    object@plots[[type]]
})

#' @rdname methods-TSENATAnalysis
setMethod("addPlot", "TSENATAnalysis", function(object, type, plot, replace = FALSE) {
    if (!replace && type %in% names(object@plots)) {
        warning("Plot type '", type, "' already exists. Set replace=TRUE to overwrite.",
            call. = FALSE)
        return(object)
    }
    object@plots[[type]] <- plot
    object
})

# Remove duplicate definitions below - all are now in consolidated form

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
        cat("  - Diversity calculations for q =", paste(gsub("q_", "", names(object@diversity_results)),
            collapse = ", "), "\n")
    }
    if (length(object@lm_results) > 0) {
        cat("  - Linear model results for components:", paste(names(object@lm_results),
            collapse = ", "), "\n")
    }
    if (length(object@jackknife_results) > 0) {
        cat("  - Jackknife results for q =", paste(gsub("q_", "", names(object@jackknife_results)),
            collapse = ", "), "\n")
    }
    if (length(object@divergence_results) > 0) {
        cat("  - Divergence results for components:", paste(names(object@divergence_results),
            collapse = ", "), "\n")
    }

    cat("\nMetadata:\n")
    if (length(object@metadata$function_calls) > 0) {
        cat("  Function calls:", paste(object@metadata$function_calls, collapse = " -> "),
            "\n")
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
