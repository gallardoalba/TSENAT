#' Generic Functions for TSENATAnalysis S4 Class
#'
#' Define generic functions for accessing components of TSENATAnalysis objects.
#' These are the recommended way to extract results from analysis objects,
#' following Bioconductor best practices.
#'
#' @name AllGenerics
#' @rdname AllGenerics
NULL

#' Access analysis results via recommended accessor methods
#'
#' @param object TSENATAnalysis object
#' @param ... Additional arguments (method-specific)
#' @details Use these accessor methods instead of direct slot access with '@'
#' @export
setGeneric("divResults", function(object, ...) standardGeneric("divResults"))

#' @rdname divResults
#' @export
setGeneric("lmRes", function(object, ...) standardGeneric("lmRes"))

#' @rdname divResults
#' @export
setGeneric("jkResults", function(object, ...) standardGeneric("jkResults"))

#' @rdname divResults
#' @export
setGeneric("divRes", function(object, ...) standardGeneric("divRes"))

#' @rdname divResults
#' @export
setGeneric("getMeta", function(object, ...) standardGeneric("getMeta"))

#' @rdname divResults
#' @export
setGeneric("getConfig", function(object, ...) standardGeneric("getConfig"))

#' @rdname divResults
#' @export
setGeneric("getPlot", function(object, ...) standardGeneric("getPlot"))

#' @rdname divResults
#' @export
setGeneric("addPlot", function(object, ...) standardGeneric("addPlot"))

#' @rdname divResults
#' @export
setGeneric("getSE", function(object, ...) standardGeneric("getSE"))
