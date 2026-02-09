## Helper: build SummarizedExperiment from readcounts + tx2gene
#' Build a SummarizedExperiment from transcript readcounts and tx->gene map
#'
#' This helper creates a `SummarizedExperiment` with an assay named
#' `counts`, stores the `tx2gene` table and the raw `readcounts` in
#' `metadata()`, and places the provided `genes` vector into `rowData(se)$genes`.
#'
#' @param tx2gene_tsv Path to a tab-separated `tx2gene` file or a data.frame
#'   with transcript-to-gene mapping. If a path is provided, it must contain
#'   at least two columns (transcript, gene).
#' @param readcounts Numeric matrix or data.frame of transcript-level counts
#'   (rows = transcripts, columns = samples).
#' @param genes Character vector of gene IDs assigning each transcript to a
#'   gene; length must equal nrow(readcounts).
#' @param assay_name Name for the assay to store readcounts (default: 'counts').
#' @return A `SummarizedExperiment` with assay, `metadata()$tx2gene`,
#'   `metadata()$readcounts` and `rowData(se)$genes` populated.
#' @export
#' @examples
#' tx2gene <- data.frame(Transcript = c("tx1", "tx2"), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
#' readcounts <- matrix(c(10, 5, 2, 3), nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))
#' genes <- c("g1", "g2")
#' se <- build_se(tx2gene, readcounts, genes)
#' SummarizedExperiment::assay(se, "counts")
build_se <- function(tx2gene_tsv, readcounts, genes, assay_name = "counts") {
    if (is.character(tx2gene_tsv) && length(tx2gene_tsv) == 1) {
        if (!file.exists(tx2gene_tsv)) {
            stop("tx2gene file not found: ", tx2gene_tsv, call. = FALSE)
        }
        tx2gene_df <- utils::read.table(tx2gene_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else if (is.data.frame(tx2gene_tsv)) {
        tx2gene_df <- tx2gene_tsv
    } else {
        stop("'tx2gene_tsv' must be a path or a data.frame.", call. = FALSE)
    }

    if (is.data.frame(readcounts)) {
        readcounts <- as.matrix(readcounts)
    }
    if (!is.matrix(readcounts) || !is.numeric(readcounts)) {
        stop("'readcounts' must be a numeric matrix or numeric data.frame.", call. = FALSE)
    }

    if (length(genes) != nrow(readcounts)) {
        stop("Length of 'genes' must equal nrow(readcounts).", call. = FALSE)
    }

    assays_list <- S4Vectors::SimpleList()
    assays_list[[assay_name]] <- readcounts

    se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list)
    S4Vectors::metadata(se)$tx2gene <- tx2gene_df
    S4Vectors::metadata(se)$readcounts <- readcounts
    # keep gene assignment in rowData for compatibility with other helpers
    # ensure rownames align with transcript IDs if present
    rn <- rownames(readcounts)
    if (is.null(rn)) {
        SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(genes = genes)
    } else {
        SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(genes = genes,
            row.names = rn)
    }

    return(se)
}
