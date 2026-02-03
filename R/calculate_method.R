#' Calculate Tsallis diversity values for transcripts grouped by gene
#' @param x Numeric matrix or data.frame of transcript-level expression
#' values (rows = transcripts, columns = samples).
#' @param genes Character vector with length equal to nrow(x) assigning each
#' transcript to a gene.
#' @param norm Logical; if TRUE normalize Tsallis entropy values per gene.
#' @param q Numeric scalar or vector of q values to evaluate.
#' @param verbose Logical; show diagnostic messages when TRUE.
#' @param what Which quantity to return from `calculate_tsallis_entropy`:
#' "S" (Tsallis entropy) or "D" (Hill numbers) (default: "S").
#' @return A data.frame with genes in the first column and per-sample (and
#' per-q) Tsallis entropy values in subsequent columns.
#' @import stats
calculate_method <- function(
  x,
  genes,
  norm = TRUE,
  verbose = FALSE,
  q = 2,
  what = c(
      "S",
      "D"
  )
) {
    what <- match.arg(what)
    # cannot use aggregate because calculate_tsallis_entropy may return
    # multiple values when length(q) > 1
    gene_levels <- unique(genes)
    coln <- as.vector(
        outer(
            colnames(x),
            q,
            function(s, qq) paste0(s, "_q=", qq)
        )
    )
    rown <- gene_levels

    # compute requested quantity (S or D)
    {
        tsallis_row <- function(gene) {
            idx <- which(genes == gene)
            unlist(lapply(seq_len(ncol(x)), function(j) {
                v <- calculate_tsallis_entropy(
                    x[
                        idx,
                        j
                    ],
                    q = q,
                    norm = norm,
                    what = what
                )
                if (length(v) == length(q) && all(is.finite(v) | is.na(v))) {
                    v
                } else {
                    names_vec <- paste0("q=", q)
                    setNames(rep(NA_real_, length(q)), names_vec)
                }
            }))
        }
        result_mat <- t(vapply(gene_levels,
            tsallis_row,
            FUN.VALUE = setNames(
                numeric(length(coln)),
                coln
            )
        ))
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
            message(sprintf("Note: %d genes excluded.", nrow(x) - nrow(y)))
        }
        colnames(y)[1] <- "Gene"
        return(y)
    }
}
