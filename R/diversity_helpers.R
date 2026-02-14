# Helpers for Tsallis entropy calculations

.tsenat_calc_S <- function(p, q, tol, n, log_base, norm) {
    vapply(q, function(qi) {
        if (abs(qi - 1) < tol) {
            sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
            if (norm) {
                if (n <= 1) {
                    # Single isoform: normalization is undefined (0/0)
                    # Return NA to indicate undefined normalization
                    sh <- NA_real_
                } else {
                    sh <- sh/log(n, base = log_base)
                }
            }
            return(sh)
        } else {
            ts <- (1 - sum(p^qi))/(qi - 1)
            if (norm) {
                max_ts <- (1 - n^(1 - qi))/(qi - 1)
                # Handle single-isoform case: ts is 0, max_ts is also 0, so 0/0 is undefined
                # Return NA to indicate undefined normalization
                if (abs(max_ts) < tol || n <= 1) {
                    ts <- NA_real_
                } else {
                    ts <- ts/max_ts
                }
            }
            return(ts)
        }
    }, numeric(1))
}

.tsenat_calc_D <- function(p, q, tol, log_base) {
    vapply(q, function(qi) {
        if (abs(qi - 1) < tol) {
            sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
            D1 <- (log_base)^(sh)
            return(D1)
        } else {
            spq <- sum(p^qi)
            Dq <- spq^(1/(1 - qi))
            return(Dq)
        }
    }, numeric(1))
}
# Input preparation
.tsenat_prepare_diversity_input <- function(x, genes = NULL, tpm = FALSE, assayno = 1,
    verbose = FALSE) {
    if (!(is.matrix(x) || is.data.frame(x) || is.list(x) || is(x, "DGEList") || is(x,
        "RangedSummarizedExperiment") || is(x, "SummarizedExperiment"))) {
        stop("Input data type is not supported! Please use ?calculate_diversity to see the possible arguments and details.",
            call. = FALSE)
    }

    if (is(x, "data.frame")) {
        x <- as.matrix(x)
    }

    if (tpm == TRUE && !is.list(x) && verbose == TRUE) {
        message("Note: tpm as a logical argument is only interpreted in case of",
            " tximport lists.")
    }

    se_assay_mat <- NULL
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
                message("Note: calculate_diversity methods are only applicable",
                  " if your DGEList contains transcript-level expression", " data.")
            }
            if (tpm == TRUE && verbose == TRUE) {
                message("Note: tpm as a logical argument is only interpreted", " in case of tximport lists.")
            }
        } else {
            stop("The package cannot find any expression data in your input.", call. = FALSE)
        }
    }

    if (is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment")) {
        md <- NULL
        try(md <- S4Vectors::metadata(x), silent = TRUE)
        if (!is.null(md) && !is.null(md$readcounts)) {
            se_assay_mat <- as.matrix(md$readcounts)
            x <- se_assay_mat
        } else {
            assays_len <- length(SummarizedExperiment::assays(x))
            if (!is.numeric(assayno) || assays_len < assayno) {
                stop("Please provide a valid assay number.", call. = FALSE)
            }
            se_assay_mat <- as.matrix(SummarizedExperiment::assays(x)[[assayno]])
            x <- se_assay_mat
        }
        if (is.null(genes)) {
            if (exists("se_assay_mat") && !is.null(md) && !is.null(md$tx2gene) &&
                is.data.frame(md$tx2gene)) {
                txmap <- md$tx2gene
                tx_col <- if ("Transcript" %in% colnames(txmap)) {
                  "Transcript"
                } else {
                  colnames(txmap)[1]
                }
                gene_col <- if ("Gen" %in% colnames(txmap)) {
                  "Gen"
                } else {
                  colnames(txmap)[2]
                }
                genes <- as.character(txmap[[gene_col]][match(rownames(se_assay_mat),
                  txmap[[tx_col]])])
                rownames(x) <- NULL
            } else {
                genes <- rownames(x)
                # keep transcript rownames available for downstream metadata
                # (we will store original transcript-level counts separately)
                rownames(x) <- NULL
            }

            if (is.null(genes)) {
                stop("Please construct a valid gene set for your ", "SummarizedExperiment.",
                  call. = FALSE)
            }
        }
    }

    list(x = x, genes = genes, se_assay_mat = se_assay_mat)
}
