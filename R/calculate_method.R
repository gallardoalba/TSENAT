#' Calculate Tsallis diversity values for transcripts grouped by gene
#' @param x Numeric matrix or data.frame of transcript-level expression
#' values (rows = transcripts, columns = samples).
#' @param genes Character vector with length equal to nrow(x) assigning each
#' transcript to a gene.
#' @param norm Logical; if TRUE normalize Tsallis entropy values per gene.
#' @param q Numeric scalar or vector of q values to evaluate.
#' @param verbose Logical; show diagnostic messages when TRUE.
#' @param what Which quantity to return from `calculate_tsallis_entropy`:
#' 'S' (Tsallis entropy) or 'D' (Hill numbers) (default: 'S').
#' @param BPPARAM BiocParallel parameter for parallel processing. Default uses
#' the registered BiocParallel backend. Use \code{BiocParallel::SerialParam()} to
#' disable parallel processing.
#' @return A data.frame with genes in the first column and per-sample (and
#' per-q) Tsallis entropy values in subsequent columns.
#' @import stats
#' @importFrom BiocParallel bplapply
calculate_method <- function(x, genes, norm = TRUE, verbose = FALSE, q = 2, what = c("S",
    "D"), BPPARAM = BiocParallel::bpparam()) {
    what <- match.arg(what)
    # validate q
    if (!is.numeric(q) || any(q <= 0)) {
        stop("Argument 'q' must be numeric and greater than 0.", call. = FALSE)
    }
    # cannot use aggregate because calculate_tsallis_entropy may return
    # multiple values when length(q) > 1
    gene_levels <- unique(genes)
    # ensure column names order matches the order used when constructing the
    # result matrix (samples vary outer, q varies inner). If sample names are
    # missing, synthesize deterministic names so column creation still works.
    sample_names <- colnames(x)
    if (is.null(sample_names)) {
        sample_names <- paste0("Sample", seq_len(ncol(x)))
    }
    coln <- as.vector(t(outer(sample_names, q, function(s, qq) paste0(s, "_q=", qq))))
    rown <- gene_levels

    # compute requested quantity ('S' or 'D') using parallel processing
    result_list <- BiocParallel::bplapply(gene_levels, function(gene) {
        .tsenat_tsallis_row(x = x, genes = genes, gene = gene, q = q, norm = norm,
            what = what)
    }, BPPARAM = BPPARAM)

    # Convert list of vectors to matrix with proper column names
    result_mat <- do.call(rbind, result_list)
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
    # Keep genes that have at least one finite value per sample-q combination
    # This allows single-isoform genes (which produce NaN when normalized) to still
    # contribute if they have valid data for some samples
    y <- x[apply(x[2:ncol(x)], 1, function(X) any(is.finite(X))), ]
    if (nrow(x) - nrow(y) > 0 && verbose == TRUE) {
        message(sprintf("Note: %d genes excluded (all non-finite values).", nrow(x) - nrow(y)))
    }
    colnames(y)[1] <- "Gene"
    return(y)
}
