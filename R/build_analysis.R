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
#' # Create example transcript count data
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples <- 10
#'
#' # Generate count matrix
#' counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
#'                  nrow = n_isoforms, ncol = n_samples)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#'
#' # Create tx2gene mapping
#' tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#'
#' # Create sample metadata
#' metadata <- data.frame(
#'   condition = rep(c("control", "treatment"), each = 5),
#'   row.names = colnames(counts))
#'
#' # Build analysis object
#' analysis <- build_analysis(
#'   readcounts = counts,
#'   tx2gene = tx2gene,
#'   metadata = metadata)
#'
#' # Verify the analysis object was created
#' analysis
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

  # Ensure sample_id column exists in colData (required by TSENATAnalysis)
  # OPTIMIZATION: Only add if not already present
  if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
    SummarizedExperiment::colData(se)$sample_id <- colnames(se)
  }

  # Store metadata in config for later use (e.g., in calculate_lm_interaction_s4)
  if (!is.null(metadata)) {
    config$metadata <- metadata
  }

  # Wrap in TSENATAnalysis
  analysis <- TSENATAnalysis(se = se, config = config)

  return(analysis)
}
