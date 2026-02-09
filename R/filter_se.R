## Filter SummarizedExperiment by low-expression transcripts
#' Filter transcripts in a `SummarizedExperiment` by minimum count/sample
#'
#' Subset a `SummarizedExperiment` to keep transcripts with more than
#' `min_count` reads in strictly more than `min_samples` samples. The
#' function updates assays, `rowData`, and relevant entries in
#' `metadata()` (for example `readcounts` and `tx2gene`) so downstream
#' helpers receive a consistent object.
#'
#' @param se A `SummarizedExperiment` object containing transcript-level
#'   assays (rows = transcripts, columns = samples).
#' @param min_count Integer threshold for per-sample counts (default 5L).
#' @param min_samples Integer minimum number of samples exceeding
#'   `min_count` required to keep a transcript (default 5L).
#' @param assay_name Name or index of the assay to use for filtering
#'   (default: 'counts'). If not present, the first assay is used.
#' @param verbose Logical; print before/after counts when TRUE.
#' @return A filtered `SummarizedExperiment`.
#' @export
filter_se <- function(se, min_count = 5L, min_samples = 5L, assay_name = "counts",
    verbose = TRUE) {
    if (!is(se, "SummarizedExperiment"))
        stop("'se' must be a SummarizedExperiment", call. = FALSE)

    assays_list <- SummarizedExperiment::assays(se)
    if (is.character(assay_name) && assay_name %in% names(assays_list)) {
        assay_mat <- as.matrix(assays_list[[assay_name]])
    } else if (is.numeric(assay_name) && assay_name >= 1 && assay_name <= length(assays_list)) {
        assay_mat <- as.matrix(assays_list[[assay_name]])
    } else {
        # fallback to first assay
        assay_mat <- as.matrix(assays_list[[1]])
        warning("Requested assay not found; using first assay.", call. = FALSE)
    }

    if (!is.numeric(assay_mat))
        stop("Assay data must be numeric.", call. = FALSE)

    before <- nrow(assay_mat)
    tokeep <- rowSums(assay_mat > min_count) > min_samples
    after <- sum(tokeep)

    if (verbose)
        message(sprintf("Transcripts: before = %d, after = %d", before, after))

    if (after == 0L) {
        warning("Filtering removed all transcripts; returning empty SummarizedExperiment.",
            call. = FALSE)
    }

    # subset all assays
    new_assays <- S4Vectors::SimpleList(lapply(assays_list, function(a) {
        if (is.matrix(a) || is.data.frame(a))
            as.matrix(a)[tokeep, , drop = FALSE] else a
    }))
    names(new_assays) <- names(assays_list)

    # subset rowData if present
    rd <- NULL
    if (nrow(SummarizedExperiment::rowData(se)) > 0) {
        rd <- SummarizedExperiment::rowData(se)[tokeep, , drop = FALSE]
    }

    # subset metadata readcounts and tx2gene if present
    md <- S4Vectors::metadata(se)
    if (!is.null(md$readcounts) && is.matrix(md$readcounts)) {
        md$readcounts <- as.matrix(md$readcounts)[tokeep, , drop = FALSE]
    }
    if (!is.null(md$tx2gene) && is.data.frame(md$tx2gene)) {
        txmap <- md$tx2gene
        txcol <- colnames(txmap)[1]
        md$tx2gene <- txmap[txmap[[txcol]] %in% rownames(assay_mat)[tokeep], , drop = FALSE]
    }

    # construct new SE with same colData
    new_se <- SummarizedExperiment::SummarizedExperiment(assays = new_assays, rowData = if (!is.null(rd))
        rd else S4Vectors::DataFrame(), colData = SummarizedExperiment::colData(se), metadata = c(S4Vectors::metadata(se),
        list(filtered = list(min_count = min_count, min_samples = min_samples))))

    # attach updated metadata pieces
    S4Vectors::metadata(new_se)$readcounts <- if (!is.null(md$readcounts))
        md$readcounts else NULL
    S4Vectors::metadata(new_se)$tx2gene <- if (!is.null(md$tx2gene))
        md$tx2gene else NULL

    return(new_se)
}
