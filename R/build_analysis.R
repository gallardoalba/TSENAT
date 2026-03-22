#' Build a Complete TSENATAnalysis Object
#'
#' Convenience wrapper that combines \code{build_se()} and \code{TSENATAnalysis()}
#' into a single function call. This creates a complete analysis object ready for
#' Tsallis entropy computation and downstream analysis.
#'
#' @param readcounts A matrix or data.frame of transcript-level read counts with
#'   transcript IDs as row names and sample names as column names. Typically output
#'   from quantification tools (SALMON, kallisto, etc.).
#'
#' @param tx2gene Either:
#'   - A path to a GFF3 or GFF3.gz file containing transcript-to-gene mapping
#'   - A path to a TSV file with columns "Transcript" and "Gene"
#'   - A data.frame with transcript-to-gene mapping
#'
#' @param assay_name Character. Name for the assay (default: "counts").
#'
#' @param metadata Optional data.frame with sample metadata. Should have sample
#'   names as row names and metadata columns (e.g., sample_type, condition, etc.).
#'
#' @param tpm Optional matrix of transcript-level TPM values. If provided, will be
#'   stored in the SummarizedExperiment. Same dimensions as readcounts required.
#'
#' @param effective_length Optional numeric vector of transcript effective lengths
#'   (e.g., from SALMON). Length should match nrow(readcounts).
#'
#' @param config Optional list of configuration parameters to store in the
#'   TSENATAnalysis object. Useful for tracking analysis parameters.
#'
#' @return A \code{TSENATAnalysis} S4 object with:
#'   \item{@se}{The SummarizedExperiment containing transcript counts and metadata}
#'   \item{@config}{Analysis configuration (empty list or user-provided)}
#'   \item{@diversity_results}{Empty list (populated by calculate_diversity_s4())}
#'   \item{@divergence_results}{Empty list (populated by calculate_divergence_s4())}
#'   \item{@lm_results}{Empty list (populated by calculate_lm_interaction_s4())}
#'   \item{@jackknife_results}{Empty list (populated by jackknife functions)}
#'   \item{@plots}{Empty list (populated by plotting functions)}
#'   \item{@metadata}{Metadata with package version and creation timestamp}
#'
#' @details
#' This wrapper combines two steps into one:
#' \enumerate{
#'   \item Call \code{build_se()} to create a SummarizedExperiment from transcript counts
#'   \item Wrap the result in \code{TSENATAnalysis()} to create the analysis object
#' }
#'
#' The returned object is ready for diversity analysis via \code{calculate_diversity_s4()}.
#'
#' If you need to inspect or filter the SummarizedExperiment before creating the
#' TSENATAnalysis object, call \code{build_se()} and \code{TSENATAnalysis()} separately.
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for the S4 class structure
#' \code{\link{calculate_diversity_s4}} for computing Tsallis entropy
#'
#' @examples
#' \dontrun{
#' # Load transcript counts
#' data("readcounts", package = "TSENAT")
#'
#' # Load metadata
#' metadata_df <- read.table(
#'   system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t", row.names = 1
#' )
#'
#' # Load annotation
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#'
#' # Create complete analysis object in one call
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   tx2gene = gff3_file,
#'   metadata = metadata_df
#' )
#'
#' # Object is ready for analysis
#' analysis <- calculate_diversity_s4(analysis, q_values = c(0, 0.5, 1, 2))
#' }
#'
#' @export
build_analysis <- function(readcounts, tx2gene, assay_name = "counts",
                          metadata = NULL, tpm = NULL, effective_length = NULL,
                          config = list()) {
  # Build SummarizedExperiment
  se <- build_se(
    readcounts = readcounts,
    tx2gene = tx2gene,
    assay_name = assay_name,
    metadata = metadata,
    tpm = tpm,
    effective_length = effective_length
  )

  # Store metadata in config for later use (e.g., in calculate_lm_interaction_s4)
  if (!is.null(metadata)) {
    config$metadata <- metadata
  }

  # Wrap in TSENATAnalysis
  analysis <- TSENATAnalysis(se = se, config = config)

  return(analysis)
}
