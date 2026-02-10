# Internal helpers for calculate_method

.tsenat_tsallis_row <- function(x, genes, gene, q, norm, what) {
    idx <- which(genes == gene)
    out <- unlist(lapply(seq_len(ncol(x)), function(j) {
        v <- calculate_tsallis_entropy(x[idx, j], q = q, norm = norm, what = what)
        if (length(v) == length(q) && all(is.finite(v) | is.na(v))) {
            v
        } else {
            names_vec <- paste0("q=", q)
            setNames(rep(NA_real_, length(q)), names_vec)
        }
    }))
    out
}
