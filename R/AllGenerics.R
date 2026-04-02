#' Generic Functions for TSENATAnalysis S4 Class
#'
#' Define generic functions for accessing components of TSENATAnalysis objects.
#' These are the recommended way to extract results from analysis objects,
#' following Bioconductor best practices.
#'
#' @return These generic functions return different types depending on the
#' specific method:
#'   - Accessor methods return data.frames or lists containing analysis results
#'   - See individual method documentation for specific return types
#'
#' @name AllGenerics
#' @rdname AllGenerics
NULL

#' Access analysis results via recommended accessor methods
#'
#' @param object TSENATAnalysis object
#' @param ... Additional arguments (method-specific)
#' @details Use these accessor methods instead of direct slot access with '@'
#' @return data.frame or list depending on accessor method. See Details for
#' specific return types.
#' @examples
#' # Create a simple TSENATAnalysis object for demonstration
#' library(SummarizedExperiment)
#' set.seed(42)
#' counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
#' colnames(counts) <- paste0('sample_', 1:10)
#' rownames(counts) <- paste0('gene_', 1:10)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' 
#' # Create TSENATAnalysis object
#' example_analysis <- new('TSENATAnalysis', se = se)
#' 
#' # Access SummarizedExperiment
#' se_retrieved <- getSE(example_analysis)
#' 
#' # These accessor methods follow Bioconductor conventions
#' # divResults, lmRes, jkResults, divRes, getMeta provide structured access
#' @export
setGeneric("divResults", function(object, ...) standardGeneric("divResults"))

#' @rdname divResults
#' @return list containing linear model interaction test results
#' @export
setGeneric("lmRes", function(object, ...) standardGeneric("lmRes"))

#' @rdname divResults
#' @return list containing jackknife diagnostics and confidence intervals
#' @export
setGeneric("jkResults", function(object, ...) standardGeneric("jkResults"))

#' @rdname divResults
#' @return list containing divergence analysis results (effect sizes, p-values)
#' @export
setGeneric("divRes", function(object, ...) standardGeneric("divRes"))

#' @rdname divResults
#' @return list or value containing stored metadata
#' @export
setGeneric("getMeta", function(object, ...) standardGeneric("getMeta"))

#' @rdname divResults
#' @return list containing analysis configuration parameters
#' @export
setGeneric("getConfig", function(object, ...) standardGeneric("getConfig"))

#' @rdname divResults
#' @return ggplot object or NULL if no plot cached
#' @export
setGeneric("getPlot", function(object, ...) standardGeneric("getPlot"))

#' @rdname divResults
#' @param type character. Plot type identifier
#' @param plot ggplot or list. The plot object to cache
#' @param replace logical. If TRUE, replace existing plot of same type.
#'   If FALSE (default), warn if plot already exists and do not overwrite.
#' @return invisible(object) for method chaining
#' @export
setGeneric("addPlot", function(object, type, plot, replace = FALSE) standardGeneric("addPlot"))

#' @rdname divResults
#' @return SummarizedExperiment containing count matrix and sample metadata
#' @export
setGeneric("getSE", function(object, ...) standardGeneric("getSE"))
