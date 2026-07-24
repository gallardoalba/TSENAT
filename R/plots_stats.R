# ============================================================================
# DIVERSITY SPECTRUM, DISTRIBUTION STATS & BOOTSTRAP CI
# Extracted from plots_helpers.R — July 2026 refactoring (I11)
# ============================================================================


# ============================================================================
# DIVERSITY SPECTRUM COMPUTATION
# ============================================================================

#' Compute Diversity Spectrum Statistics
#'
#' Aggregates diversity measurements across q-values and groups.
#' Calculates median/mean and variability (IQR/SD) for each q-value.
#'
#' @param se A \code{SummarizedExperiment} with diversity assays.
#' @param q_values Numeric vector of q-values to compute (optional,
#' auto-detect if NULL).
#' @param metric Character: 'median' (default) or 'mean' for central tendency.
#' @param variability_metric Character: 'iqr' (default) or 'sd' for spread.
#' @param condition_col Character: column name for grouping conditions
#' (optional).
#'
#' @return Data frame with columns:
#'   - q: q-value
#'   - group: condition group (if condition_col provided)
#'   - central: median or mean divergence
#'   - spread: IQR or SD of divergence
#'   - count: number of valid measurements
#'

#' @noRd

.compute_diversity_spectrum <- function(se, q_values = NULL, metric = c("median",
    "mean"), variability_metric = c("iqr", "sd"), condition_col = NULL) {


    # Validate input
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    if (nrow(se) == 0 || ncol(se) == 0) {
        stop("SummarizedExperiment is empty", call. = FALSE)
    }

    # Match arguments
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)

    # Prepare long format data
    long_data <- .prepare_tsallis_long(se, assay_name = "diversity", condition_col = condition_col)

    if (nrow(long_data) == 0) {
        stop("No valid diversity data found in SummarizedExperiment", call. = FALSE)
    }

    # Ensure q is numeric
    long_data$q <- as.numeric(as.character(long_data$q))

    # Compute statistics by group and q-value
    if (!is.null(condition_col) && condition_col %in% colnames(long_data)) {
        # Group by condition
        stats <- long_data %>%
            dplyr::group_by(group, q) %>%
            dplyr::summarise(central = if (metric == "median") {
                median(.data$tsallis, na.rm = TRUE)
            } else {
                mean(.data$tsallis, na.rm = TRUE)
            }, spread = if (variability_metric == "iqr") {
                IQR(.data$tsallis, na.rm = TRUE)
            } else {
                sqrt(stats::var(.data$tsallis, na.rm = TRUE))
            }, count = sum(!is.na(.data$tsallis)), .groups = "drop")
    } else {
        # No grouping
        stats <- long_data %>%
            dplyr::group_by(q) %>%
            dplyr::summarise(central = if (metric == "median") {
                median(.data$tsallis, na.rm = TRUE)
            } else {
                mean(.data$tsallis, na.rm = TRUE)
            }, spread = if (variability_metric == "iqr") {
                IQR(.data$tsallis, na.rm = TRUE)
            } else {
                sqrt(stats::var(.data$tsallis, na.rm = TRUE))
            }, count = sum(!is.na(.data$tsallis)), .groups = "drop")
    }

    return(stats)
}


# ============================================================================
# GENE FILTERING & RANKING
# ============================================================================

#' Select Top Genes by P-Value
#'
#' Ranks genes by statistical significance and selects top N.
#'
#' @param results Data frame with at least one p-value column.
#' @param p_col Character: column name for p-values
#'   ('adj_p_interaction', 'p_interaction', 'padj', 'pvalue').
#' @param gene_col Character: column name for gene identifiers
#'   ('gene_id', 'gene', 'gene_name').
#' @param n_genes Integer: number of top genes to select (default: 4).
#'
#' @return Character vector of top gene IDs, sorted by p-value (smallest first).
#'
#' @noRd
.select_top_genes <- function(results, p_col = NULL, gene_col = NULL, n_genes = 4) {


    if (!is.data.frame(results) || nrow(results) == 0) {
        stop("results must be a non-empty data frame", call. = FALSE)
    }

    # Auto-detect p-value column
    if (is.null(p_col)) {
        candidate_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
        matched <- candidate_cols[candidate_cols %in% colnames(results)]
        if (length(matched) > 0) {
            p_col <- matched[1]
        } else {
            stop("Could not find p-value column. ", "Provide p_col explicitly.",
                call. = FALSE)
        }
    }

    # Auto-detect gene column
    if (is.null(gene_col)) {
        candidate_cols <- c("gene_id", "gene", "gene_name")
        matched <- candidate_cols[candidate_cols %in% colnames(results)]
        if (length(matched) > 0) {
            gene_col <- matched[1]
        } else {
            stop("Could not find gene column. ", "Provide gene_col explicitly.",
                call. = FALSE)
        }
    }

    # Select top genes
    top_genes <- results %>%
        dplyr::arrange(.data[[p_col]]) %>%
        dplyr::slice(seq_len(min(n_genes, nrow(results)))) %>%
        dplyr::pull(.data[[gene_col]])

    return(as.character(top_genes))
}

#' Filter Genes by Significance Threshold
#'
#' Selects genes with p-value below threshold.
#'
#' @param results Data frame with p-values and gene identifiers.
#' @param p_threshold Numeric: p-value cutoff (default: 0.05).
#' @param p_col Character: p-value column name (auto-detected if NULL).
#' @param gene_col Character: gene identifier column (auto-detected if NULL).
#'
#' @return Character vector of significant gene IDs.
#'

#' @noRd

.filter_genes_by_pvalue <- function(results, p_threshold = 0.05, p_col = NULL, gene_col = NULL) {


    if (!is.data.frame(results) || nrow(results) == 0) {
        stop("results must be a non-empty data frame", call. = FALSE)
    }

    # Auto-detect columns (same logic as select_top_genes)
    if (is.null(p_col)) {
        candidate_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
        matched <- candidate_cols[candidate_cols %in% colnames(results)]
        if (length(matched) > 0) {
            p_col <- matched[1]
        } else {
            stop("Could not find p-value column", call. = FALSE)
        }
    }

    if (is.null(gene_col)) {
        candidate_cols <- c("gene_id", "gene", "gene_name")
        matched <- candidate_cols[candidate_cols %in% colnames(results)]
        if (length(matched) > 0) {
            gene_col <- matched[1]
        } else {
            stop("Could not find gene column", call. = FALSE)
        }
    }

    # Filter and return
    sig_genes <- results %>%
        dplyr::filter(.data[[p_col]] < p_threshold) %>%
        dplyr::arrange(.data[[p_col]]) %>%
        dplyr::pull(.data[[gene_col]])

    return(as.character(sig_genes))
}

# ============================================================================

# ============================================================================

#' Compute Distribution Statistics by Group
#'
#' Consolidates the common pattern of dplyr group-by + summarize for
#' calculating central tendency and spread measures.
#'
#' @param df Data frame with data to summarize
#' @param group_col Character: column name for grouping (e.g., 'group', 'condition')
#' @param value_col Character: column name for values to summarize (e.g., 'entropy')
#' @param metric Character: central tendency - 'median' (default) or 'mean'
#' @param spread_metric Character: spread measure - 'iqr' (default), 'sd'
#'
#' @return Data frame with columns:
#'   - group_col: group identifier
#'   - value: central tendency (median or mean)
#'   - lower: lower bound of spread
#'   - upper: upper bound of spread
#'
#' @details
#' **Consolidation Impact**: 7-8 occurrences × 4-5 lines = 28-40 LOC saved
#'
#' Replaces repetitive patterns like:
#' ```
#' df %>%
#'   dplyr::group_by(group) %>%
#'   dplyr::summarize(
#'       value = median(col, na.rm = TRUE),
#'       lower = quantile(col, 0.25, na.rm = TRUE),
#'       upper = quantile(col, 0.75, na.rm = TRUE),
#'       .groups = 'drop'
#'   )
#' ```
#'
#' @examples
#' \dontrun{
#' df <- data.frame(group = rep(c('A', 'B'), 50), value = rnorm(100))
#' stats <- .compute_distribution_stats(df, 'group', 'value', 'median', 'iqr')
#' head(stats)
#' #   group     value     lower     upper
#' # 1     A -0.123456 -0.654321 0.234567
#' # 2     B  0.234567 -0.345678 0.876543
#' }
#'
#' @noRd
.compute_distribution_stats <- function(df, group_col, value_col, metric = "median",
    spread_metric = "iqr") {

    # Validate inputs
    if (!is.data.frame(df)) {
        stop("df must be a data frame", call. = FALSE)
    }
    if (!group_col %in% colnames(df)) {
        stop("Group column '", group_col, "' not found in data frame", call. = FALSE)
    }
    if (!value_col %in% colnames(df)) {
        stop("Value column '", value_col, "' not found in data frame", call. = FALSE)
    }
    if (!is.numeric(df[[value_col]])) {
        stop("Value column '", value_col, "' is not numeric", call. = FALSE)
    }

    # Define central tendency function
    central_fn <- if (metric == "median") {
        function(x) stats::median(x, na.rm = TRUE)
    } else if (metric == "mean") {
        function(x) mean(x, na.rm = TRUE)
    } else {
        stop("Unknown metric: ", metric, call. = FALSE)
    }

    # Calculate statistics
    if (spread_metric == "iqr") {
        # Compute for each group separately to avoid dplyr quantile issues
        groups <- unique(df[[group_col]])
        stats_list <- lapply(groups, function(grp) {
            grp_data <- df[[value_col]][df[[group_col]] == grp]
            data.frame(group = grp, value = central_fn(grp_data), lower = as.numeric(stats::quantile(grp_data,
                0.25, na.rm = TRUE)), upper = as.numeric(stats::quantile(grp_data,
                0.75, na.rm = TRUE)))
        })
        names(stats_list) <- NULL
        stats_df <- do.call(rbind, stats_list)
        colnames(stats_df)[1] <- group_col
    } else if (spread_metric == "sd") {
        # Compute for each group separately to avoid dplyr binding issues
        groups <- unique(df[[group_col]])
        stats_list <- lapply(groups, function(grp) {
            grp_data <- df[[value_col]][df[[group_col]] == grp]
            val <- central_fn(grp_data)
            sd_val <- stats::sd(grp_data, na.rm = TRUE)
            data.frame(group = grp, value = val, lower = val - sd_val, upper = val +
                sd_val)
        })
        names(stats_list) <- NULL
        stats_df <- do.call(rbind, stats_list)
        colnames(stats_df)[1] <- group_col
    } else {
        stop("Unknown spread_metric: ", spread_metric, call. = FALSE)
    }

    return(stats_df)
}


# ============================================================================

#' Extract Bootstrap Confidence Interval Assays
#'
#' Consolidates detection and extraction of bootstrap CI assays from
#' SummarizedExperiment objects.
#'
#' @param se SummarizedExperiment object
#' @param assay_name Character: base assay name (default: 'diversity')
#' @param fallback_to_iqr Logical: if CIs missing, return fallback indicator
#'   (default: TRUE)
#'
#' @return List with elements:
#'   - has_ci: Logical, TRUE if both ci_lower and ci_upper assays exist
#'   - ci_lower: Matrix or NULL if not found
#'   - ci_upper: Matrix or NULL if not found
#'   - assay_base: The base assay matrix
#'   - fallback_metric: Character ('iqr' or NULL) indicating fallback method
#'
#' @details
#' **Consolidation Impact**: 5-6 occurrences × 5-6 lines = 25-36 LOC saved
#'
#' Replaces patterns like:
#' ```
#' ci_lower_name <- paste0(assay_name, '_ci_lower')
#' ci_upper_name <- paste0(assay_name, '_ci_upper')
#' has_ci <- all(c(ci_lower_name, ci_upper_name) %in% SummarizedExperiment::assayNames(se))
#' if (has_ci) {
#'     ci_lower <- SummarizedExperiment::assay(se, ci_lower_name)
#'     ci_upper <- SummarizedExperiment::assay(se, ci_upper_name)
#' }
#' ```
#'
#' @examples
#' \dontrun{
#' ci_result <- .extract_bootstrap_ci_assays(se, assay_name = 'diversity')
#' if (ci_result$has_ci) {
#'     ci_lower <- ci_result$ci_lower
#'     ci_upper <- ci_result$ci_upper
#'     # use CIs
#' } else if (ci_result$fallback_metric == 'iqr') {
#'     # fall back to IQR
#' }
#' }
#'
#' @noRd
.extract_bootstrap_ci_assays <- function(se, assay_name = "diversity", fallback_to_iqr = TRUE) {

    # Validate base assay exists
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
        stop("Assay '", assay_name, "' not found in SummarizedExperiment", call. = FALSE)
    }

    # Get base assay
    assay_base <- SummarizedExperiment::assay(se, assay_name)

    # Check for CI assays with standard naming convention
    ci_lower_name <- paste0(assay_name, "_ci_lower")
    ci_upper_name <- paste0(assay_name, "_ci_upper")

    has_ci_lower <- ci_lower_name %in% SummarizedExperiment::assayNames(se)
    has_ci_upper <- ci_upper_name %in% SummarizedExperiment::assayNames(se)
    has_ci <- has_ci_lower && has_ci_upper

    # Extract CI assays if present
    ci_lower <- if (has_ci_lower) {
        SummarizedExperiment::assay(se, ci_lower_name)
    } else {
        NULL
    }

    ci_upper <- if (has_ci_upper) {
        SummarizedExperiment::assay(se, ci_upper_name)
    } else {
        NULL
    }

    # Determine fallback strategy if CIs missing
    fallback_metric <- if (!has_ci && fallback_to_iqr)
        "iqr" else NULL

    return(list(has_ci = has_ci, ci_lower = ci_lower, ci_upper = ci_upper, assay_base = assay_base,
        fallback_metric = fallback_metric))
}
.infer_samples_from_se <- function(se, samples = NULL, condition_col = "condition") {
    if (!is.null(samples)) {
        return(as.character(samples))
    }
    cd <- NULL
    try(cd <- SummarizedExperiment::colData(se), silent = TRUE)
    if (is.null(cd)) {
        return(NULL)
    }

    # Common column names to try
    candidates <- c(condition_col, "condition", "group", "sample_group", "sampleType",
        "class", "status", "phenotype")
    for (nm in candidates) {
        if (nm %in% colnames(cd)) {
            return(as.character(cd[[nm]]))
        }
    }

    # Fallback: choose column with smallest >1 unique values
    uniq_counts <- vapply(cd, function(col) length(unique(na.omit(col))), integer(1))
    valid_cols <- names(uniq_counts[uniq_counts > 1])
    if (length(valid_cols) > 0) {
        bin_cols <- valid_cols[uniq_counts[valid_cols] == 2]
        pick <- if (length(bin_cols) > 0)
            bin_cols[1] else valid_cols[which.min(uniq_counts[valid_cols])]
        return(as.character(cd[[pick]]))
    }

    NULL
}


.get_readcounts_from_se <- function(se, readcounts_arg = NULL) {
    # If user provided a readcounts object/path, accept it first
    if (!is.null(readcounts_arg)) {
        if (is.character(readcounts_arg) && length(readcounts_arg) == 1) {
            if (!file.exists(readcounts_arg))
                stop("readcounts file not found: ", readcounts_arg)
            rc_df <- utils::read.delim(readcounts_arg, header = TRUE, stringsAsFactors = FALSE)
            if (!is.null(colnames(rc_df)) && ncol(rc_df) > 1) {
                counts <- as.matrix(rc_df[, -1, drop = FALSE])
                rownames(counts) <- rc_df[[1]]
            } else {
                counts <- as.matrix(rc_df)
            }
            return(counts)
        } else if (is.matrix(readcounts_arg) || is.data.frame(readcounts_arg)) {
            return(as.matrix(readcounts_arg))
        } else {
            stop("`readcounts` must be a matrix/data.frame or path to a file")
        }
    }

    md <- NULL
    try(md <- S4Vectors::metadata(se), silent = TRUE)
    if (!is.null(md) && !is.null(md$readcounts)) {
        return(as.matrix(md$readcounts))
    }

    assay_names <- SummarizedExperiment::assayNames(se)
    preferred_assays <- c("readcounts", "counts", "tx_counts", "counts_tx")
    chosen <- intersect(preferred_assays, assay_names)
    if (length(chosen) > 0) {
        return(as.matrix(SummarizedExperiment::assay(se, chosen[1])))
    }

    # fallback to first assay with a warning
    warning("Using first assay from SummarizedExperiment to compute", " expression-based fold changes; ensure it contains",
        " transcript-level readcounts or provide metadata$readcounts")
    as.matrix(SummarizedExperiment::assay(se))
}


.get_tx2gene_from_se <- function(se, readcounts_mat = NULL) {
    md <- NULL
    try(md <- S4Vectors::metadata(se), silent = TRUE)
    # prefer explicit tx2gene in metadata
    if (!is.null(md) && !is.null(md$tx2gene) && is.data.frame(md$tx2gene)) {
        txmap <- md$tx2gene
        # attempt to find sensible columns
        tx_col <- if ("Transcript" %in% colnames(txmap))
            "Transcript" else colnames(txmap)[1]
        gene_col <- if ("Gen" %in% colnames(txmap))
            "Gen" else colnames(txmap)[2]
        return(list(type = "vector", mapping = as.character(txmap[[gene_col]][match(rownames(readcounts_mat),
            txmap[[tx_col]])])))
    }

    # fallback: try rowData mapping
    rdata <- SummarizedExperiment::rowData(se)
    if (!is.null(rdata) && "genes" %in% colnames(rdata)) {
        genes_vec <- as.character(rdata$genes)
        if (!is.null(readcounts_mat) && length(genes_vec) == nrow(readcounts_mat)) {
            return(list(type = "vector", mapping = genes_vec))
        }
    }

    # last resort: use rownames of readcounts as gene identifiers
    if (!is.null(readcounts_mat)) {
        return(list(type = "vector", mapping = rownames(readcounts_mat)))
    }

    NULL
}


.validate_control_in_samples <- function(control, samples) {
    uniq <- unique(samples)
    if (!is.null(control) && control %in% uniq) {
        return(control)
    }
    if ("Normal" %in% uniq) {
        return("Normal")
    }
    # fallback to first level and message
    chosen <- uniq[1]
    message(sprintf("`control` not found; using '%s' instead", chosen))
    chosen
}


#' Violin plot of Tsallis entropy for a single q value
#'
#' Creates a violin plot showing the distribution of Tsallis entropy for a
#' specific q value,
#' with groups (conditions) displayed side by side.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`
#' containing
#'   entropy values at one or more q values.
#' @param q_value The specific q value to plot (numeric, e.g., 1, 2, 0.5).
#' @param assay_name Name of the assay to use (default: 'diversity').
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a violin plot with groups on the x-axis.
#'  
#' @noRd

