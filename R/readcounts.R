#' Example transcript-level read counts dataset
#'
#' A dataset containing transcript-level transcript-level read count abundances
#' from salmon quantification. The dataset includes samples as columns
#' (SRR14800475-SRR14800490) and transcript IDs as row names.
#'
#' @docType data
#'
#' @usage data(readcounts)
#'
#' @format A data frame with 3514 transcripts (rows) and 16 samples (columns).
#' Row names are transcript IDs (ENST format) and column names are sample IDs.
#' Values represent TPM-normalized read counts.
#'
#' @return A matrix (data frame) containing transcript-level TPM-normalized read
#'   counts from salmon quantification. Rows represent transcripts (ENST format
#'   IDs) and columns represent samples. Aliases include \code{salmon_dataset}
#'   (main read count matrix), \code{salmon_tpm} (TPM values), and
#'   \code{salmon_effective_length} (transcript lengths from salmon).
#'
#' @keywords datasets
#'
#' @examples
#' data(readcounts, package = "TSENAT")
#' dim(salmon_dataset)
#' head(salmon_dataset[1:4, 1:3])
#'
#' @name readcounts
#' @aliases readcounts salmon_dataset salmon_tpm salmon_effective_length
NULL


