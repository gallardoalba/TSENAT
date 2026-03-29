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
#' @return
#' Generic function that dispatches to specific methods:
#' \item{divResults}{data.frame} Tsallis entropy or Hill numbers for genes across q-values
#' \item{lmRes}{list} Linear model results including interaction tests
#' \item{jkResults}{list} Jackknife diagnostics and confidence intervals
#' \item{divRes}{list} Divergence analysis results (effect sizes, p-values)
#' \item{getMeta}{list or value} Metadata stored in object
#' \item{getConfig}{list} Analysis configuration parameters
#' \item{getSE}{SummarizedExperiment} The count matrix and sample metadata
#' \item{getPlot}{ggplot or NULL} Cached plot objects
#' @examples
#' # Create a simple TSENATAnalysis object for demonstration
#' library(SummarizedExperiment)
#' set.seed(42)
#' counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
#' colnames(counts) <- paste0("sample_", 1:10)
#' rownames(counts) <- paste0("gene_", 1:10)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' 
#' # Create TSENATAnalysis object
#' example_analysis <- new("TSENATAnalysis", se = se)
#' 
#' # Access SummarizedExperiment
#' se_retrieved <- getSE(example_analysis)
#' 
#' # These accessor methods follow Bioconductor conventions
#' # divResults, lmRes, jkResults, divRes, getMeta provide structured access
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
