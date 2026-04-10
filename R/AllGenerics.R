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

# ============================================================================
# GLOBAL VARIABLES DECLARATION
# ============================================================================
# Declare variables used in data.table and ggplot2 non-standard evaluation
# (NSE) across the package to suppress R CMD check NOTEs about undefined global
# variables
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Gene", "group", "tsallis", "tx", "df_model", "legend_name",
        "legend_position", "log2expr", "dimension", "variable", "contribution", "dim1",
        "dim2", "type", "coord_x", "coord_y"))
}

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
#' # lmResults, jeoResults, jisResults, getMeta provide structured access

#' Get metadata from TSENATAnalysis
#' @param object TSENATAnalysis object
#' @param key Optional metadata key to retrieve
#' @return list containing metadata or specific key value
#' @noRd
setGeneric("getMeta", function(object, ...) standardGeneric("getMeta"))

#' Get configuration from TSENATAnalysis
#' @param object TSENATAnalysis object
#' @return list containing configuration parameters
#' @noRd
setGeneric("getConfig", function(object, ...) standardGeneric("getConfig"))

#' Get cached plot
#'
#' @param object \code{TSENATAnalysis} object.
#' @param type \code{character}. Plot type: 'q_curve', 'lm_interaction',
#'   'divergence', 'influence', 'volcano', etc.
#'   If NULL, returns all cached plots.
#'
#' @return ggplot object or list of plots.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @noRd
setGeneric("getPlot", function(object, ...) standardGeneric("getPlot"))

#' Add or cache a plot in TSENATAnalysis object
#'
#' Cache plots for later retrieval, maintaining analysis visualization history.
#'
#' @param object TSENATAnalysis object
#' @param type character. Plot type identifier
#' @param plot ggplot or list. The plot object to cache
#' @param replace logical. If TRUE, replace existing plot of same type.
#'   If FALSE (default), warn if plot already exists and do not overwrite.
#' @return invisible(object) for method chaining
#'
#' @examples
#' # Create a simple plot and add to analysis
#' data(readcounts, package = 'TSENAT')
#' metadata <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#'
#' # TPM and effective_length REQUIRED for filter_analysis()
#' tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
#'               nrow = nrow(readcounts), ncol = ncol(readcounts),
#'               dimnames = dimnames(readcounts))
#' effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))
#' 
#' config <- TSENAT_config(q = 1.0, generate_plots = FALSE)
#' analysis <- build_analysis(readcounts, tx2gene = gff3_file,
#'     metadata = metadata, tpm = tpm, effective_length = effective_length,
#'     config = config)
#' analysis <- filter_analysis(analysis, stringency = 'severe')
#' analysis <- calculate_diversity(analysis, norm = TRUE)
#' 
#' # Create and cache a plot
#' p <- plot_diversity_spectrum(analysis)
#' @noRd
setGeneric("addPlot", function(object, type, plot, replace = FALSE) standardGeneric("addPlot"))

#' Extract SummarizedExperiment from TSENATAnalysis
#'
#' Retrieve the underlying SummarizedExperiment object containing count data
#' and sample metadata.
#'
#' @param object TSENATAnalysis object
#' @param ... Additional arguments (for method compatibility)
#'
#' @return SummarizedExperiment containing count matrix and sample metadata
#' 
#' @details
#' The SummarizedExperiment object returned by getSE() contains:
#' - assays: count matrices (transcript-level read counts)
#' - rowData: transcript information and gene assignments
#' - colData: sample metadata (sample types, conditions, etc.)
#'
#' @examples
#' # Load example data and build analysis object
#' data(readcounts, package = 'TSENAT')
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#'
#' # Build TSENATAnalysis object
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, 
#'                              tx2gene = gff3_file, 
#'                              metadata = metadata_df,
#'                              config = config,
#'                              tpm = tpm,
#'                              effective_length = effective_length)
#'
#' # Extract the underlying SummarizedExperiment
#' se <- getSE(analysis)
#' 
#' # Explore the SummarizedExperiment structure
#' nrow(se)  # Number of transcripts
#' ncol(se)  # Number of samples
#' SummarizedExperiment::assayNames(se)  # Available assay matrices
#' @noRd
setGeneric("getSE", function(object, ...) standardGeneric("getSE"))

#' @noRd
setGeneric("setConfig", function(object, value) standardGeneric("setConfig"))

#' @noRd
setGeneric("setConfigValue", function(object, key, value) standardGeneric("setConfigValue"))

#' Extract SummarizedExperiment from TSENATAnalysis
#'
#' @param object TSENATAnalysis object
#' @return SummarizedExperiment object
#' @examples
#' # Create a simple TSENATAnalysis object
#' library(SummarizedExperiment)
#' counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
#' colnames(counts) <- paste0('sample_', 1:10)
#' rownames(counts) <- paste0('gene_', 1:10)
#'
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' analysis <- new('TSENATAnalysis', se = se, config = list())
#'
#' # Extract the SummarizedExperiment
#' extracted_se <- se(analysis)
#' dim(extracted_se)
#' @export
setGeneric("se", function(object) standardGeneric("se"))

#' @noRd
if (!isGeneric("metadata")) {
    setGeneric("metadata", function(x, key = NULL) standardGeneric("metadata"))
}

#' Set metadata for TSENATAnalysis object
#'
# metadata<- generic for TSENATAnalysis
# Defines the replacement method for modifying @metadata slot
#' @param x TSENATAnalysis object
#' @param value Replacement value for metadata (typically a list)
#' @return TSENATAnalysis object with updated metadata
#' @examples
#' # Create a simple TSENATAnalysis object
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(counts = matrix(1:100, nrow = 10)))
#' analysis <- new('TSENATAnalysis', se = se, config = list())
#'
#' # Set metadata
#' metadata(analysis) <- list(processing_date = Sys.Date(), method = "test")
#'
#' # Verify metadata was set
#' metadata(analysis)
#' @seealso \code{\link{TSENATAnalysis-methods}} for other accessor methods
#' @keywords internal
setGeneric("metadata<-", function(x, value) standardGeneric("metadata<-"))
