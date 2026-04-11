#' Constructor for TSENATAnalysis objects
#'
#' Creates a new TSENATAnalysis object with a SummarizedExperiment base
#' and optional initial configuration.
#'
#' @param se \code{SummarizedExperiment}. The base expression data object.
#' @param config \code{list}. Optional initial configuration (usually set
#'   via \code{TSENAT_config()} instead).
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
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
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
#' @export
TSENATAnalysis <- function(se, config = list()) {
    # Validate input
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    # Convert TSENATConfig S4 object to list if needed
    if (inherits(config, "TSENATConfig")) {
        config <- unclass(config)
    }

    # NEW: Validate config parameters early (before object creation)
    if (length(config) > 0 && !is.null(config)) {
        cdata <- SummarizedExperiment::colData(se)

        # Validate condition_col
        if ("condition_col" %in% names(config)) {
            col <- config$condition_col
            if (!is.null(col) && !col %in% colnames(cdata)) {
                stop(sprintf("Invalid condition_col '%s': column not found in colData.\\nAvailable columns: %s",
                  col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
            }
        }

        # Validate subject_col if paired
        if ("paired" %in% names(config) && isTRUE(config$paired)) {
            if ("subject_col" %in% names(config)) {
                col <- config$subject_col
                if (!is.null(col) && !col %in% colnames(cdata)) {
                  stop(sprintf("Invalid subject_col '%s': column not found in colData.\\nAvailable columns: %s",
                    col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
                }
            }
        }
    }

    # Ensure sample_id column exists in colData (required by validator)
    if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
        SummarizedExperiment::colData(se)$sample_id <- colnames(se)
    }

    # Create new object with all slots initialized
    new("TSENATAnalysis", se = se, config = if (length(config) > 0)
        config else list(), diversity_results = list(), lm_results = list(), jackknife_results = list(),
        divergence_results = list(), plots = list(), metadata = list(created_at = Sys.time(),
            package_version = as.character(utils::packageVersion("TSENAT")), function_calls = character()))
}
