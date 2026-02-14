## Helper: build SummarizedExperiment from readcounts + tx2gene
#' Build a SummarizedExperiment from transcript readcounts and tx->gene map
#'
#' This helper creates a `SummarizedExperiment` with an assay named
#' `counts`, stores the `tx2gene` table and the raw `readcounts` in
#' `metadata()`, and extracts gene assignments from the tx2gene mapping
#' to populate `rowData(se)$genes`.
#'
#' @param readcounts Numeric matrix or data.frame of transcript-level counts
#'   (rows = transcripts, columns = samples). Rownames should contain
#'   transcript IDs that match the first column of tx2gene.
#' @param tx2gene Path to a tab-separated `tx2gene` file or a data.frame
#'   with transcript-to-gene mapping. Must contain at least two columns:
#'   Transcript (first column) and Gene (second column).
#' @param assay_name Name for the assay to store readcounts (default: 'counts').
#' @return A `SummarizedExperiment` with assay, `metadata()$tx2gene`,
#'   `metadata()$readcounts` and `rowData(se)$genes` populated.
#' @export
#' @examples
#' tx2gene <- data.frame(Transcript = c('tx1', 'tx2'), Gene = c('g1', 'g2'))
#' readcounts <- matrix(c(10, 5, 2, 3), nrow = 2, dimnames = list(c('tx1', 'tx2'), c('s1', 's2')))
#' se <- build_se(readcounts, tx2gene)
#' SummarizedExperiment::assay(se, 'counts')
build_se <- function(readcounts, tx2gene, assay_name = "counts") {
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }
        tx2gene_df <- utils::read.table(tx2gene, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else if (is.data.frame(tx2gene)) {
        tx2gene_df <- tx2gene
    } else {
        stop("'tx2gene' must be a path or a data.frame.", call. = FALSE)
    }

    if (is.data.frame(readcounts)) {
        readcounts <- as.matrix(readcounts)
    }
    if (!is.matrix(readcounts) || !is.numeric(readcounts)) {
        stop("'readcounts' must be a numeric matrix or numeric data.frame.", call. = FALSE)
    }

    # Extract transcript IDs from rownames or use row indices
    tx_ids <- rownames(readcounts)
    if (is.null(tx_ids)) {
        stop("'readcounts' must have transcript IDs as rownames.", call. = FALSE)
    }

    # Extract gene names from tx2gene based on transcript IDs
    tx_col <- if ("Transcript" %in% colnames(tx2gene_df)) "Transcript" else colnames(tx2gene_df)[1]
    genes <- tx2gene_df$Gene[match(tx_ids, tx2gene_df[[tx_col]])]
    
    if (any(is.na(genes))) {
        stop("Some transcript IDs in readcounts were not found in tx2gene mapping.", call. = FALSE)
    }

    assays_list <- S4Vectors::SimpleList()
    assays_list[[assay_name]] <- readcounts

    se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list)
    S4Vectors::metadata(se)$tx2gene <- tx2gene_df
    S4Vectors::metadata(se)$readcounts <- readcounts
    
    # Populate rowData with gene assignments and transcript IDs
    SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(
        genes = genes,
        row.names = tx_ids
    )

    return(se)
}
