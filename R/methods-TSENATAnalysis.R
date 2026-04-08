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

# NOTE: show() and summary() methods are defined in s4_class.R
# This avoids duplication and keeps all S4 method definitions in one place
