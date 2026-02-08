#' Map external coldata into a SummarizedExperiment
#'
#' Map an external `coldata` table (sample IDs and a condition/label column)
#' into a `SummarizedExperiment` produced by `calculate_diversity()`. The
#' helper sets `colData(ts_se)$sample_type` when possible and records a
#' `sample_base` identifier for each column. Unmatched entries remain NA and
#' no automatic inference is attempted.
#' @param ts_se A `SummarizedExperiment` whose assay columns are named with
#'   sample identifiers. Names may include per-sample annotations; mapping is
#'   exact and case-sensitive unless you normalize identifiers beforehand.
#' @param coldata A data.frame with sample metadata (or `NULL`).
#' @param coldata_sample_col Name of the column in `coldata` containing sample
#'   identifiers (default: "Sample").
#' @param coldata_condition_col Name of the column in `coldata` with
#'   condition/labels (default: "Condition").
#' @param paired Logical; if `TRUE`, validate pairing and reorder columns so
#'   that matched samples for each base are adjacent (default: `FALSE`). Use
#'   `paired = TRUE` for paired analyses so downstream paired tests see
#'   aligned samples regardless of the original `coldata` order.
#' @return The input `ts_se` with `colData(ts_se)$sample_type` populated when
#' possible. `colData(ts_se)$sample_base` is also added containing the base
#' sample identifier; rownames of `colData(ts_se)` are aligned to the assay
#' column names.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
#' sample_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se)))
#' coldata_df <- data.frame(Sample = sample_names, Condition = rep(c("A", "B"),
#'     length.out = ncol(se)
#' ))
#' map_metadata(se, coldata_df)
#' # Optionally validate pairs when appropriate
#' map_metadata(se, coldata_df, paired = TRUE)
map_metadata <- function(
    ts_se,
    coldata,
    coldata_sample_col = "Sample",
    coldata_condition_col = "Condition",
    paired = FALSE
) {
    if (is.null(coldata)) {
        return(ts_se)
    }
    if (!is.data.frame(coldata) || !all(c(
        coldata_sample_col,
        coldata_condition_col
    ) %in% colnames(coldata))) {
        return(ts_se)
    }
    # Prepare canonical condition order and base extraction.
    # `conds` is sorted for deterministic ordering when constructing paired
    # column order. `coldata_base` removes a trailing role suffix (e.g.
    # `_N`/`_T`) from the `coldata` sample IDs; note that assay column names
    # may contain per-sample annotations (e.g. per-q tags) which are handled
    # separately below when building `sample_base_names`.
    conds <- sort(unique(as.character(coldata[[coldata_condition_col]])))
    coldata_base <- sub(
        "_[^_]+$", "",
        as.character(coldata[[coldata_sample_col]])
    )
    bases <- unique(coldata_base)

    # Optionally validate that coldata defines pairing between conditions
    if (isTRUE(paired) && length(conds) >= 2) {
        unpaired <- vapply(bases, function(b) {
            length(unique(coldata[[coldata_condition_col]][coldata_base == b]))
        }, integer(1))
        bad <- bases[unpaired != length(conds)]
        if (length(bad) > 0) {
            bad_list <- paste(bad, collapse = ", ")
            cond_list <- paste(conds, collapse = ", ")
            msg <- paste0(
                "Unpaired samples found in coldata for bases: ",
                bad_list,
                ". Ensure each base has all conditions: ",
                cond_list
            )
            stop(msg, call. = FALSE)
        }
    }
    sample_base_names <- sub(
        "_q=.*",
        "",
        colnames(SummarizedExperiment::assay(ts_se))
    )

    # Reorder the SummarizedExperiment columns.
    base_names <- sample_base_names
    idx_list <- integer(0)
    if (isTRUE(paired)) {
        # Build paired-aware order: for each base, list samples in the
        # deterministic `conds` order so matched samples become adjacent.
        for (b in bases) {
            for (c in conds) {
                # find sample names in coldata for this base+condition
                samples_for_pair <- as.character(
                    coldata[[coldata_sample_col]]
                )[coldata_base == b &
                    as.character(coldata[[coldata_condition_col]]) == c]
                if (length(samples_for_pair) == 0) next
                # for each sample name, find matching assay columns (after q-stripping)
                for (s in samples_for_pair) {
                    matches <- which(
                        base_names == s
                    )
                    if (length(matches) > 0) idx_list <- c(idx_list, matches)
                }
            }
        }
        # Preserve original order for any unmatched columns
        remaining <- setdiff(seq_along(base_names), idx_list)
        new_order <- c(idx_list, remaining)
    } else {
        # Non-paired: follow the order present in `coldata`
        ordered_samples <- as.character(coldata[[coldata_sample_col]])
        for (s in ordered_samples) {
            matches <- which(base_names == s)
            if (length(matches) > 0) idx_list <- c(idx_list, matches)
        }
        remaining <- setdiff(seq_along(base_names), idx_list)
        new_order <- c(idx_list, remaining)
    }
    if (length(new_order) > 0 && !all(new_order == seq_along(base_names))) {
        ts_se <- ts_se[, new_order, drop = FALSE]
        sample_base_names <- sub(
            "_q=.*", "",
            colnames(SummarizedExperiment::assay(ts_se))
        )
    }
    st_map <- setNames(
        as.character(coldata[[coldata_condition_col]]),
        as.character(coldata[[coldata_sample_col]])
    )
    sample_types <- unname(st_map[sample_base_names])
    missing_idx <- which(is.na(sample_types))
    if (length(missing_idx) > 0) {
        missing_samples <- sample_base_names[missing_idx]
        msg <- paste0(
            "map_metadata: unmatched samples in 'coldata': ",
            paste(missing_samples, collapse = ", "),
            ". Provide matching entries in 'coldata' or populate ",
            "colData(ts_se)$sample_type beforehand."
        )
        stop(msg, call. = FALSE)
    }
    SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
    # Also record the base sample names (without q suffix) per column so
    # downstream functions can rely on explicit pairing information.
    SummarizedExperiment::colData(ts_se)$sample_base <- sample_base_names
    # Ensure rownames of colData match assay column names
    assay_cols <- colnames(SummarizedExperiment::assay(ts_se))
    rownames(SummarizedExperiment::colData(ts_se)) <- assay_cols
    # Attach transcript-level readcounts and tx->gene mapping to metadata
    # if they are available in the calling environment or globalenv and
    # not already present in the SummarizedExperiment metadata. This
    # simplifies downstream plotting helpers that expect these objects.
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
        md <- S4Vectors::metadata(ts_se)
        # prefer existing metadata values; otherwise try common names
        if (is.null(md$readcounts)) {
            if (exists("readcounts", envir = parent.frame())) {
                md$readcounts <- get("readcounts", envir = parent.frame())
            } else if (exists("readcounts", envir = globalenv())) {
                md$readcounts <- get("readcounts", envir = globalenv())
            }
        }
        if (is.null(md$tx2gene)) {
            # vignette uses 'txmap' variable name; also accept 'tx2gene'
            if (exists("txmap", envir = parent.frame())) {
                md$tx2gene <- get("txmap", envir = parent.frame())
            } else if (exists("tx2gene", envir = parent.frame())) {
                md$tx2gene <- get("tx2gene", envir = parent.frame())
            } else if (exists("txmap", envir = globalenv())) {
                md$tx2gene <- get("txmap", envir = globalenv())
            } else if (exists("tx2gene", envir = globalenv())) {
                md$tx2gene <- get("tx2gene", envir = globalenv())
            }
        }
        S4Vectors::metadata(ts_se) <- md
    }
    # If a diversity assay is present, prepare a simple diversity data.frame
    # (genes + per-sample diversity values) and store it in metadata so the
    # vignette and plotting helpers can use a ready-made table.
    if ("diversity" %in% SummarizedExperiment::assayNames(ts_se)) {
        div_mat <- as.matrix(SummarizedExperiment::assay(ts_se, "diversity"))
        # sample base names without per-q suffixes
        sample_base_names <- sub("_q=.*", "", colnames(div_mat))
        # prefer explicit sample_type in colData when present
        samples_vec <- NULL
        if ("sample_type" %in% colnames(SummarizedExperiment::colData(ts_se))) {
            samples_vec <- as.character(SummarizedExperiment::colData(ts_se)$sample_type)
        }
        div_df <- as.data.frame(div_mat)
        genes_col <- if (!is.null(SummarizedExperiment::rowData(ts_se)$genes)) {
            SummarizedExperiment::rowData(ts_se)$genes
        } else {
            rownames(div_df)
        }
        div_df <- cbind(genes = genes_col, div_df)
        md2 <- S4Vectors::metadata(ts_se)
        md2$diversity_df <- div_df
        md2$sample_base_names <- sample_base_names
        if (!is.null(samples_vec)) md2$samples <- samples_vec
        S4Vectors::metadata(ts_se) <- md2
    }
    return(ts_se)
}

# Map sample names (without '_q=...') to group labels using `colData(se)`.
# Mapping must be provided via `colData(se)`; no inference fallback is used.
map_samples_to_group <- function(sample_names,
                                 se = NULL,
                                 sample_type_col = NULL,
                                 mat = NULL) {
    # Prefer explicit mapping from colData(se)[, sample_type_col] when
    # provided. If `sample_type_col` is not provided, allow a single-
    # condition dataset by assigning a single default group "Group" to
    # all samples (this permits plotting single-condition q-curves).
    base_names <- sub(
        "_q=.*",
        "",
        colnames(if (!is.null(mat)) mat else SummarizedExperiment::assay(se))
    )

    if (!is.null(se) && !is.null(sample_type_col) &&
        (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st_vec <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st_vec) <- base_names
        st_map <- st_vec[!duplicated(names(st_vec))]
    } else {
        # No explicit mapping: assume single-group dataset
        st_map <- setNames(rep("Group", length(base_names)), base_names)
    }

    mapped <- unname(st_map[sample_names])
    missing_idx <- which(is.na(mapped))
    if (length(missing_idx) > 0) {
        stop(
            sprintf(
                "Missing sample_type mapping for samples: %s",
                paste(unique(sample_names[missing_idx]), collapse = ", ")
            )
        )
    }
    mapped
}

# Prepare a long-format data.frame for a simple assay (one value per sample)
get_assay_long <- function(se,
                           assay_name = "diversity",
                           value_name = "diversity",
                           sample_type_col = NULL) {
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")
    
    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) stop("Assay not found: ", assay_name)
    df <- as.data.frame(mat)
    genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) {
        SummarizedExperiment::rowData(se)$genes
    } else {
        rownames(df)
    }
    df <- cbind(df, Gene = genes_col)
    long <- tidyr::pivot_longer(
        df,
        -Gene,
        names_to = "sample",
        values_to = value_name
    )

    # sample_type: prefer explicit colData mapping when available. If not
    # provided, assume a single-group dataset and set `sample_type` to
    # "Group" for all samples.
    if (!is.null(sample_type_col) && (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st) <- colnames(mat)
        st_map <- st[!duplicated(names(st))]
        sample_base <- sub("_q=.*", "", long$sample)
        long$sample_type <- unname(st_map[sample_base])
        missing_idx <- which(is.na(long$sample_type))
        if (length(missing_idx) > 0) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(sample_base[missing_idx]), collapse = ", ")))
        }
    } else {
        long$sample_type <- rep("Group", nrow(long))
    }

    long[!is.na(long[[value_name]]), , drop = FALSE]
}

# Internal small helper: prepare long-format tsallis data from a
# SummarizedExperiment
prepare_tsallis_long <- function(se,
                                 assay_name = "diversity",
                                 sample_type_col = "sample_type") {
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")
    
    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) stop("Assay not found: ", assay_name)
    df <- as.data.frame(mat)
    genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) {
        SummarizedExperiment::rowData(se)$genes
    } else {
        rownames(df)
    }
    df <- cbind(df, Gene = genes_col)

    long <- tidyr::pivot_longer(
        df,
        -Gene,
        names_to = "sample_q",
        values_to = "tsallis"
    )
    if (any(grepl("_q=", long$sample_q))) {
        long <- tidyr::separate(
            long,
            sample_q,
            into = c("sample", "q"),
            sep = "_q="
        )
        long$q <- as.factor(as.numeric(long$q))
    } else {
        long$sample <- long$sample_q
        long$q <- NA
    }

    if (!is.null(sample_type_col) && (sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        st_vec <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st_vec) <- sub("_q=.*", "", colnames(mat))
        st_map <- st_vec[!duplicated(names(st_vec))]
        long$group <- unname(st_map[as.character(long$sample)])
        missing_idx <- which(is.na(long$group))
        if (length(missing_idx) > 0) {
            stop(sprintf("Missing sample_type mapping for samples: %s", paste(unique(as.character(long$sample)[missing_idx]), collapse = ", ")))
        }
    } else {
        long$group <- rep("Group", nrow(long))
    }

    as.data.frame(long[!is.na(long$tsallis), , drop = FALSE])
}


#' Map transcript IDs from a tx2gene table to a readcounts matrix
#'
#' Assign transcript identifiers as row names of a transcript-level read
#' counts matrix so downstream plotting functions can identify transcripts.
#' The tx2gene mapping can be provided as a file path (TSV) or a
#' data.frame with a `Transcript` column.
#'
#' @param readcounts A numeric matrix or data.frame of read counts (rows = transcripts).
#' @param tx2gene Either a path to a tab-delimited file or a data.frame with at least a `Transcript` column.
#' @param tx_col Name of the transcript ID column in `tx2gene` (default: "Transcript").
#' @param verbose Logical; print informative messages (default: FALSE).
#' @return The input `readcounts` with rownames set to the transcript IDs.
#' @export
map_tx_to_readcounts <- function(readcounts, tx2gene, tx_col = "Transcript", verbose = FALSE) {
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        txmap <- utils::read.delim(tx2gene, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(tx2gene)) {
        txmap <- tx2gene
    } else {
        stop("'tx2gene' must be a file path or a data.frame.", call. = FALSE)
    }

    if (!(tx_col %in% colnames(txmap))) stop(sprintf("tx2gene mapping must contain column '%s'", tx_col), call. = FALSE)

    # Ensure readcounts is a matrix-like object
    if (is.data.frame(readcounts)) readcounts <- as.matrix(readcounts)
    if (!is.matrix(readcounts)) stop("'readcounts' must be a matrix or data.frame", call. = FALSE)

    n_rc <- nrow(readcounts)
    n_tx <- nrow(txmap)

    if (n_tx == n_rc) {
        rownames(readcounts) <- as.character(txmap[[tx_col]])
        if (verbose) message(sprintf("Assigned %d transcript rownames from tx2gene.", n_rc))
        return(readcounts)
    }

    # If counts differ, attempt to match by transcript identifiers if present
    tx_ids <- as.character(txmap[[tx_col]])
    if (!is.null(rownames(readcounts)) && all(rownames(readcounts) %in% tx_ids)) {
        # reorder txmap to match readcounts row order and set rownames
        matched_idx <- match(rownames(readcounts), tx_ids)
        rownames(readcounts) <- tx_ids[matched_idx]
        if (verbose) message("Matched and assigned transcript IDs by existing readcounts rownames.")
        return(readcounts)
    }

    stop(sprintf("Number of transcripts in tx2gene (%d) does not match readcounts rows (%d), and automatic matching failed.", n_tx, n_rc), call. = FALSE)
}
