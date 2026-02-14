#' TCGA Luminal A breast cancer dataset (transcript-level)
#'
#' Data from The Cancer Genome Atlas, containing transcript-level read counts
#' of 20 patients with Luminal A type breast cancer (primary tumor and solid
#' normal samples, 40 samples total). The dataset includes transcript IDs in
#' the first column and read counts for 40 samples in subsequent columns.
#'
#' @docType data
#'
#' @usage data(tcga_brca_luma)
#'
#' @format A data frame with 1100 rows and 41 columns. The first column contains
#' transcript IDs, all additional columns contain RNA-sequencing read counts for
#' samples.
#'
#' @keywords datasets
#'
#' @references The Cancer Genome Atlas Network (2012) Nature 490, 61–70
#' \doi{10.1038/nature11412}
#'
#' @source \href{https://portal.gdc.cancer.gov/}{TCGA via GDC}
#'
#' @examples
#' data(tcga_brca_luma)
#' dim(tcga_brca_luma)
#' head(tcga_brca_luma[1:4, 1:3])
"tcga_brca_luma"
