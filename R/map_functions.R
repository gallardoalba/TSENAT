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
#' map_coldata_to_se(se, coldata_df)
#' # Optionally validate pairs when appropriate
#' map_coldata_to_se(se, coldata_df, paired = TRUE)
map_coldata_to_se <- function(
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
        # Leave unmatched sample types as NA; do not attempt inference.
        sample_types[missing_idx] <- NA_character_
    }
    SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
    # Also record the base sample names (without q suffix) per column so
    # downstream functions can rely on explicit pairing information.
    SummarizedExperiment::colData(ts_se)$sample_base <- sample_base_names
    # Ensure rownames of colData match assay column names
    assay_cols <- colnames(SummarizedExperiment::assay(ts_se))
    rownames(SummarizedExperiment::colData(ts_se)) <- assay_cols
    return(ts_se)
}
## Note: `infer_sample_group()` removed; prefer explicit `coldata` mapping
## or supply `sample_type` via `colData(se)`.
# Supports underscore-suffixed labels (e.g. Sample_N, Sample_T) and
