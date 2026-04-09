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
    # Returns ONLY essential metadata (timestamps, version, workflow type)
    # Large result tables, sample stats, function logs are accessed via results(), etc.
    
    # Filter metadata to essential fields only
    essential_meta <- list()
    
    if ("created_at" %in% names(object@metadata)) {
        essential_meta$created_at <- object@metadata$created_at
    }
    if ("ended_at" %in% names(object@metadata)) {
        essential_meta$ended_at <- object@metadata$ended_at
    }
    if ("package_version" %in% names(object@metadata)) {
        essential_meta$package_version <- object@metadata$package_version
    }
    if ("tsenat_version" %in% names(object@metadata)) {
        essential_meta$tsenat_version <- object@metadata$tsenat_version
    }
    if ("workflow_type" %in% names(object@metadata)) {
        essential_meta$workflow_type <- object@metadata$workflow_type
    }
    if ("workflow" %in% names(object@metadata)) {
        # Extract only essential workflow info
        workflow_info <- object@metadata$workflow
        if (is.list(workflow_info)) {
            essential_meta$workflow <- list(
                type = workflow_info$workflow_type,
                completion_time = workflow_info$completion_time
            )
        }
    }
    
    # If key specified, return specific field
    if (!is.null(key)) {
        return(essential_meta[[key]])
    }
    
    # Return filtered essential metadata only
    essential_meta
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

# NOTE: show() and summary() methods are defined in s4_class.R
# This avoids duplication and keeps all S4 method definitions in one place
