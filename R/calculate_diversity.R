#' Calculate Tsallis diversity per gene across samples
#'
#' @param x A numeric matrix or data.frame of transcript-level expression
#' #'   values (rows = transcripts, columns = samples), or a SummarizedExperiment-
#' #' like object.
#' @param tpm Logical. If TRUE and `x` is a tximport-style list, use the
#'   `$abundance` matrix instead of `$counts`.
#' @param genes Character vector assigning each transcript (row) to a gene.
#'   Must have length equal to nrow(x) or the number of transcripts in `x`.
#' @param norm Logical; if TRUE, normalize Tsallis entropy to [0,1] per gene.
#' #' @param assayno Integer assay index to use when `x` is a SummarizedExperiment.
#' @param verbose Logical; print diagnostic messages when TRUE.
#' @param q Numeric scalar or vector of Tsallis q values to evaluate (q > 0).
#' #'   If length(q) > 1, the result will contain separate columns per sample and
#' #' q.
#' #' @param what Which quantity to return: "S" for Tsallis entropy, "D" for Hill
#' #' numbers, or "both".
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with assay
#'   `diversity` containing per-gene diversity values.
#' @import methods
#' #' @importFrom SummarizedExperiment SummarizedExperiment assays assay rowData
#' #' colData
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
#' SummarizedExperiment::assay(se)[1:3, 1:3]
calculate_diversity <- function(x, genes = NULL, norm = TRUE,
                                tpm = FALSE, assayno = 1, verbose = FALSE, q = 2, what = c("S", "D")) {
  if (!(is.matrix(x) || is.data.frame(x) || is.list(x) || is(x, "DGEList") ||
    is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment"))) {
    stop("Input data type is not supported! Please use `?calculate_diversity`
\t to see the possible arguments and details.")
  }

  if (is(x, "data.frame")) {
    x <- as.matrix(x)
  }

  if (tpm == TRUE && !is.list(x) && verbose == TRUE) {
    message("Note: tpm as a logical argument is only interpreted in case of
            tximport lists.")
  }

  if (is.list(x)) {
    if (length(x) == 4 && "counts" %in% names(x)) {
      if (tpm == FALSE) {
        x <- as.matrix(x$counts)
      }
      if (tpm == TRUE) {
        x <- as.matrix(x$abundance)
      }
    } else if (is(x, "DGEList")) {
      x <- as.matrix(x$counts)
      if (verbose == TRUE) {
        message("Note: calculate_diversity methods are only applicable if your
                DGEList contains transcript-level expression data.")
      }
      if (tpm == TRUE && verbose == TRUE) {
        message("Note: tpm as a logical argument is only interpreted in case of
              tximport lists.")
      }
    } else {
      stop("The package cannot find any expression data in your input.", call. = FALSE)
    }
  }

  if (is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment")) {
    if (!is.numeric(assayno) || length(SummarizedExperiment::assays(x)) < assayno) {
      stop("Please give a valid number to pick an assay from your data.", call. = FALSE)
    } else if (is.numeric(assayno)) {
      x <- as.matrix(SummarizedExperiment::assays(x)[[assayno]])
    }
    if (is.null(genes)) {
      genes <- rownames(x)
      rownames(x) <- NULL
      if (is.null(genes)) {
        stop("Please construct a valid gene set for your SummarizedExperiment.",
          call. = FALSE
        )
      }
    }
  }

  if (!is.numeric(x)) {
    stop("Input data  must be numeric!", call. = FALSE)
  }

  if (any(is.na(x))) {
    stop("The data contains NA as expression values. NAs are not allowed in the
         input.", call. = FALSE)
  }

  if (nrow(x) != length(genes)) {
    stop("The number of rows is not equal to the given gene set.", call. = FALSE)
  }

  what <- match.arg(what)
  result <- calculate_method(x, genes, norm, verbose = verbose, q = q, what = what)

  # Prepare assay and row/col data
  result_assay <- result[, -1, drop = FALSE]
  result_rowData <- data.frame(genes = result[, 1], row.names = result[, 1])

  if (length(q) > 1) {
    col_split <- do.call(rbind, strsplit(colnames(result)[-1], "_q="))
    col_ids <- paste0(col_split[, 1], "_q=", col_split[, 2])
    row_ids <- as.character(result[, 1])
    result_colData <- data.frame(
      samples = as.character(col_split[, 1]),
      q = as.numeric(col_split[, 2]),
      row.names = col_ids,
      stringsAsFactors = FALSE
    )
    colnames(result_assay) <- col_ids
    rownames(result_assay) <- row_ids
  } else {
    col_ids <- colnames(x)
    row_ids <- as.character(result[, 1])
    result_colData <- data.frame(samples = col_ids, row.names = col_ids)
    colnames(result_assay) <- col_ids
    rownames(result_assay) <- row_ids
  }

  result_metadata <- list(method = "tsallis", norm = norm, q = q, what = what)

  assays_list <- list()
  if (what == "S") {
    assays_list$diversity <- result_assay
  } else if (what == "D") {
    assays_list$hill <- result_assay
  }

  result <- SummarizedExperiment::SummarizedExperiment(
    assays = assays_list,
    rowData = result_rowData,
    colData = result_colData,
    metadata = result_metadata
  )

  return(result)
}
