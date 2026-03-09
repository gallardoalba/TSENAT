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
#' @keywords datasets
#'
#' @examples
#' \dontrun{
#'   data(readcounts, package = "TSENAT")
#'   dim(readcounts)
#'   head(readcounts[1:4, 1:3])
#' }
#'
#' @name readcounts
NULL

