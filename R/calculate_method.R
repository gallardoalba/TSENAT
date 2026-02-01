#' Calculate Tsallis diversity values for transcripts grouped by gene
#'
#' This helper computes per-gene Tsallis entropy across samples. The trimmed
#' package only supports the Tsallis method; other diversity metrics were
#' removed.
#'
#' @param x Numeric matrix or data.frame of transcript-level expression
#'   values (rows = transcripts, columns = samples).
#' @param genes Character vector with length equal to nrow(x) assigning each
#'   transcript to a gene.
#' @param method Only "tsallis" is supported (default).
#' @param norm Logical; if TRUE normalize Tsallis entropy values per gene.
#' @param q Numeric scalar or vector of q values to evaluate.
#' @param verbose Logical; show diagnostic messages when TRUE.
#' @param what Which quantity to return from `calculate_tsallis_entropy`: "S", "D" or "both" (default: "S").
#' @return A data.frame with genes in the first column and per-sample (and
#'   per-q) Tsallis entropy values in subsequent columns.
#' @import stats
calculate_method <- function(x, genes, method = "tsallis", norm = TRUE, verbose = FALSE, q = 2, what = c("S", "D", "both")) {
  # `method` argument accepted for compatibility; only Tsallis aggregation is performed.
  what <- match.arg(what)
  # cannot use aggregate because calculate_tsallis_entropy may return
  # multiple values when length(q) > 1
  gene_levels <- unique(genes)
  coln <- as.vector(outer(colnames(x), q, function(s, qq) paste0(s, "_q=", qq)))
  rown <- gene_levels

  if (what != "both") {
    tsallis_row <- function(gene) {
      idx <- which(genes == gene)
      unlist(lapply(seq_len(ncol(x)), function(j) {
        v <- calculate_tsallis_entropy(x[idx, j], q = q, norm = norm, what = what)
        if (length(v) == length(q) && all(is.finite(v) | is.na(v))) v else setNames(rep(NA_real_, length(q)), paste0("q=", q))
      }))
    }
    result_mat <- t(vapply(gene_levels, tsallis_row, FUN.VALUE = setNames(numeric(length(coln)), coln)))
    colnames(result_mat) <- coln
    rownames(result_mat) <- rown
    out_df <- data.frame(Gene = rown, result_mat, check.names = FALSE)
    if (all(rowSums(!is.na(result_mat)) == 0)) {
      out_df <- data.frame(Gene = character(0))
      for (nm in coln) out_df[[nm]] <- numeric(0)
      x <- out_df
      return(x)
    }
    x <- out_df
    y <- x[apply(x[2:ncol(x)], 1, function(X) all(is.finite(X))), ]
    if (nrow(x) - nrow(y) > 0 && verbose == TRUE) {
      message(paste0("Note: There are ", nrow(x) - nrow(y), " genes with single isoforms, which will be exluded from the analysis."))
    }
    colnames(y)[1] <- "Gene"
    return(y)
  } else {
    # compute both S and D matrices
    S_row <- function(gene) {
      idx <- which(genes == gene)
      unlist(lapply(seq_len(ncol(x)), function(j) {
        both <- calculate_tsallis_entropy(x[idx, j], q = q, norm = norm, what = "both")
        if (!is.null(both$S) && length(both$S) == length(q)) both$S else setNames(rep(NA_real_, length(q)), paste0("q=", q))
      }))
    }
    D_row <- function(gene) {
      idx <- which(genes == gene)
      unlist(lapply(seq_len(ncol(x)), function(j) {
        both <- calculate_tsallis_entropy(x[idx, j], q = q, norm = norm, what = "both")
        if (!is.null(both$D) && length(both$D) == length(q)) both$D else setNames(rep(NA_real_, length(q)), paste0("q=", q))
      }))
    }
    S_mat <- t(vapply(gene_levels, S_row, FUN.VALUE = setNames(numeric(length(coln)), coln)))
    D_mat <- t(vapply(gene_levels, D_row, FUN.VALUE = setNames(numeric(length(coln)), coln)))
    colnames(S_mat) <- coln
    colnames(D_mat) <- coln
    rownames(S_mat) <- rown
    rownames(D_mat) <- rown
    out_S <- data.frame(Gene = rown, S_mat, check.names = FALSE)
    out_D <- data.frame(Gene = rown, D_mat, check.names = FALSE)
    colnames(out_S)[1] <- "Gene"
    colnames(out_D)[1] <- "Gene"
    return(list(S = out_S, D = out_D))
  }
}
