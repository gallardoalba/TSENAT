# # Helpers for generating diversity plots (density, violin, MA)

#' @importFrom utils head
#' @importFrom grDevices png dev.off
if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c(
            "Gene",
            "diversity",
            "sample",
            "sample_q",
            "sample_type",
            "fold",
            "significant",
            "value",
            ".",
            "group",
            "tsallis",
            "q",
            "median",
            "IQR",
            "padj_num",
            "padj_clean",
            "xval",
            "label_flag",
            "genes",
            "mean",
            # transcript plotting globals
            "tx",
            "expr",
            "tx_cond",
            "log2expr",
            # generic plotting helpers
            "x",
            "y",
            "padj"
        )
    )

    require_pkgs <- function(pkgs) {
        missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing)) {
            stop(
                sprintf(
                    "%s required",
                    paste(missing, collapse = ", ")
                )
            )
        }
        invisible(TRUE)
    }
}

# Small utility helpers to reduce repeated code paths when extracting
# samples, readcounts and tx->gene mappings from a
# SummarizedExperiment. Keeping these as focused helpers improves
# readability of the longer plotting functions below.
infer_samples_from_se <- function(se, samples = NULL, sample_type_col = "sample_type") {
    if (!is.null(samples)) {
        return(as.character(samples))
    }
    cd <- NULL
    try(cd <- SummarizedExperiment::colData(se), silent = TRUE)
    if (is.null(cd)) {
        return(NULL)
    }

    # Common column names to try
    candidates <- c(
        sample_type_col,
        "condition",
        "group",
        "sample_group",
        "sampleType",
        "class",
        "status",
        "phenotype"
    )
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
        pick <- if (length(bin_cols) > 0) bin_cols[1] else valid_cols[which.min(uniq_counts[valid_cols])]
        return(as.character(cd[[pick]]))
    }

    NULL
}

get_readcounts_from_se <- function(se, readcounts_arg = NULL) {
    # If user provided a readcounts object/path, accept it first
    if (!is.null(readcounts_arg)) {
        if (is.character(readcounts_arg) && length(readcounts_arg) == 1) {
            if (!file.exists(readcounts_arg)) stop("readcounts file not found: ", readcounts_arg)
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
    warning(
        "Using first assay from SummarizedExperiment to compute",
        " expression-based fold changes; ensure it contains",
        " transcript-level readcounts or provide metadata$readcounts"
    )
    as.matrix(SummarizedExperiment::assay(se))
}

get_tx2gene_from_se <- function(se, readcounts_mat = NULL) {
    md <- NULL
    try(md <- S4Vectors::metadata(se), silent = TRUE)
    # prefer explicit tx2gene in metadata
    if (!is.null(md) && !is.null(md$tx2gene) && is.data.frame(md$tx2gene)) {
        txmap <- md$tx2gene
        # attempt to find sensible columns
        tx_col <- if ("Transcript" %in% colnames(txmap)) "Transcript" else colnames(txmap)[1]
        gene_col <- if ("Gen" %in% colnames(txmap)) "Gen" else colnames(txmap)[2]
        return(list(
            type = "vector",
            mapping = as.character(txmap[[gene_col]][match(rownames(readcounts_mat), txmap[[tx_col]])])
        ))
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

validate_control_in_samples <- function(control, samples) {
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

#' Plot q-curve profile for a single gene comparing groups
#'
#' For a selected gene, plot per-sample Tsallis entropy across q values and
#' overlay per-group centrality estimates with variability ribbons so group-level
#' differences are easy to compare. Expects a `SummarizedExperiment` produced by
#' `calculate_diversity()` with `_q=` suffixes in column names.
#'
#' @param se A `SummarizedExperiment` from `calculate_diversity()`.
#' @param gene Character scalar or vector; gene symbol(s) to plot. If NULL and
#'   `lm_res` is supplied, the top `n_top` genes from `lm_res` (by
#'   `adj_p_interaction` or `p_interaction`) are used.
#' @param lm_res Optional data.frame result from `calculate_lm_interaction()`.
#'   When supplied and `gene` is NULL, the top `n_top` significant genes will
#'   be plotted.
#' @param n_top Number of top genes to plot when `lm_res` is provided (default: 10).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Column name in `colData(se)` with sample type labels
#'   (default: "sample_type"). If missing, a single-group fallback is used.
#' @param show_samples Logical; if TRUE, draw per-sample lines in the
#'   background (default: FALSE).
#' @param metric Central tendency metric to use (default: "median"). Options are
#'   "median" or "mean".
#' @param variability_metric Variability metric to display as ribbon (default: "IQR").
#'   Options are "IQR" (interquartile range, shows ±IQR/2 around the central value) or
#'   "variance" (shows ±1 standard deviation around the central value).
#' @return A `ggplot` object when a single gene is requested, or a named list
#'   of `ggplot` objects when multiple genes are requested.
#' @examples
#' mat <- matrix(runif(8), nrow = 2, dimnames = list(c("g1", "g2"), c("s1_q=0.1", "s1_q=1", "s2_q=0.1", "s2_q=1")))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
#' plot_tsallis_gene_profile(se, gene = "g1")
#' plot_tsallis_gene_profile(se, gene = "g1", metric = "mean", variability_metric = "variance")
#' @export
plot_tsallis_gene_profile <- function(se,
                                      gene = NULL,
                                      lm_res = NULL,
                                      n_top = 10,
                                      assay_name = "diversity",
                                      sample_type_col = "sample_type",
                                      show_samples = FALSE,
                                      metric = c("median", "mean"),
                                      variability_metric = c("IQR", "variance")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    require_pkgs(c("dplyr", "tidyr", "SummarizedExperiment"))

    # Validate metric parameters
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)

    long <- prepare_tsallis_long(se, assay_name = assay_name, sample_type_col = sample_type_col)
    if (!("Gene" %in% colnames(long))) stop("prepare_tsallis_long did not return Gene column")

    # Resolve genes to plot: accept NULL (use lm_res), a single name, or a vector/list
    if (is.null(gene)) {
        if (is.null(lm_res)) stop("Either 'gene' or 'lm_res' must be provided")
        if (!is.data.frame(lm_res) || !("gene" %in% colnames(lm_res))) stop("'lm_res' must be a data.frame with a 'gene' column")
        # prefer adj_p_interaction if present
        pcol <- if ("adj_p_interaction" %in% colnames(lm_res)) "adj_p_interaction" else if ("p_interaction" %in% colnames(lm_res)) "p_interaction" else NULL
        if (is.null(pcol)) stop("'lm_res' must contain 'adj_p_interaction' or 'p_interaction' columns")
        genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[pcol]])]))
        genes <- head(genes_ordered, n_top)
    } else {
        genes <- as.character(unlist(gene))
    }

    if (length(genes) == 0) stop("No genes selected for plotting")

    # helper to build single plot for a gene
    make_plot_for_gene <- function(sel) {
        long_g <- long[as.character(long$Gene) == sel, , drop = FALSE]
        if (nrow(long_g) == 0) stop("Gene not found in assay: ", sel)
        long_g$qnum <- as.numeric(as.character(long_g$q))

        # Compute central tendency and variability based on selected metrics
        if (variability_metric == "IQR") {
            # Central tendency + IQR
            if (metric == "median") {
                stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
                    central = median(tsallis, na.rm = TRUE),
                    spread = stats::IQR(tsallis, na.rm = TRUE),
                    .groups = "drop"
                )
                spread_factor <- 1/2  # IQR/2 for symmetric ribbon
                spread_label <- "IQR"
            } else {  # mean
                stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
                    central = mean(tsallis, na.rm = TRUE),
                    spread = stats::IQR(tsallis, na.rm = TRUE),
                    .groups = "drop"
                )
                spread_factor <- 1/2  # IQR/2 for symmetric ribbon
                spread_label <- "IQR"
            }
        } else {  # variance
            # Central tendency + Standard Deviation (sqrt of variance)
            if (metric == "median") {
                stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
                    central = median(tsallis, na.rm = TRUE),
                    spread = sqrt(stats::var(tsallis, na.rm = TRUE)),
                    .groups = "drop"
                )
            } else {  # mean
                stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
                    central = mean(tsallis, na.rm = TRUE),
                    spread = sqrt(stats::var(tsallis, na.rm = TRUE)),
                    .groups = "drop"
                )
            }
            spread_factor <- 1  # ±1 SD for variance ribbon
            spread_label <- "SD"
        }

        # Build plot
        p <- ggplot2::ggplot() +
            ggplot2::theme_minimal(base_size = 14)

        if (isTRUE(show_samples)) {
            p <- p + ggplot2::geom_line(data = long_g, ggplot2::aes(x = qnum, y = tsallis, group = sample, color = group), alpha = 0.25)
        }

        # Create descriptive title showing metric and variability choices
        metric_label <- if (metric == "median") "Median" else "Mean"
        variability_label <- if (variability_metric == "IQR") "IQR" else "SD"
        title_suffix <- sprintf(" (%s ± %s)", metric_label, variability_label)

        p <- p +
            ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum, ymin = central - spread * spread_factor, ymax = central + spread * spread_factor, fill = group), alpha = 0.2, inherit.aes = FALSE) +
            ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = central, color = group), linewidth = 1.3) +
            ggplot2::labs(title = paste0(sel, ": Tsallis entropy q-curve profile", title_suffix), x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
            ggplot2::scale_color_discrete(name = "Group") + ggplot2::scale_fill_discrete(name = "Group") +
            ggplot2::theme(plot.title = ggplot2::element_text(
                hjust = 0.5, size = 16,
                margin = ggplot2::margin(b = 10)
            ))
        p
    }

    # Return single ggplot for single gene, or a named list of ggplots for multiple genes
    if (length(genes) == 1) {
        return(make_plot_for_gene(genes))
    }

    plots <- lapply(genes, make_plot_for_gene)
    names(plots) <- genes
    plots
}

#' Plot diversity distributions (density) by sample type
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Optional column name in `colData(se)` with sample
#' types.  If missing, sample types are inferred from column names (suffix after
#' the last underscore) or set to 'Group'.
#' @return A `ggplot` object with layered density plots.
#' @importFrom ggplot2 ggplot aes geom_density facet_grid scale_color_manual
#' guides theme_minimal labs
#' @export
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:100, -1, drop = FALSE])
#' gs <- readcounts[1:100, 1]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = FALSE)
#' # Manually set sample_type in colData for plotting
#' SummarizedExperiment::colData(se)$sample_type <- 
#'   factor(gsub(".*_", "", colnames(rc)))
#' plot_diversity_density(se)
plot_diversity_density <- function(
  se,
  assay_name = "diversity",
  sample_type_col = NULL
) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr", "SummarizedExperiment"))
    
    # For plot_diversity_density, sample_type is required for faceting
    # Determine which column to use for sample_type
    if (is.null(sample_type_col)) {
        # Check if "sample_type" exists in colData as a fallback
        if (!("sample_type" %in% colnames(SummarizedExperiment::colData(se)))) {
            stop("sample_type column not found in data", call. = FALSE)
        }
        sample_type_col <- "sample_type"
    } else if (!(sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
        stop(sprintf("Column '%s' not found in colData.", sample_type_col), call. = FALSE)
    }
    
    # Check for all NA sample_type values
    st_col <- SummarizedExperiment::colData(se)[[sample_type_col]]
    if (all(is.na(st_col))) {
        stop("All sample_type values are NA. Cannot create faceted plot.", call. = FALSE)
    }
    
    long <- get_assay_long(
        se,
        assay_name = assay_name,
        value_name = "diversity",
        sample_type_col = sample_type_col
    )

    # Ensure sample_type column exists and has no NA values
    if (!("sample_type" %in% colnames(long))) {
        stop("sample_type column not found in data", call. = FALSE)
    }
    if (all(is.na(long$sample_type))) {
        stop("All sample_type values are NA. Cannot create faceted plot.", call. = FALSE)
    }

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = diversity, group = sample, color = sample_type)
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_grid(. ~ sample_type) +
        ggplot2::scale_color_manual(values = c("black", "darkorchid4")) +
        ggplot2::guides(color = "none") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(x = "Diversity values", y = "Density")
}


#' Plot violin of per-gene mean diversity by sample type
#' @importFrom magrittr %>%
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Optional column name in `colData(se)` containing
#' sample types.
#' @return A `ggplot` violin plot object.
#' @export
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:100, -1, drop = FALSE])
#' gs <- readcounts[1:100, 1]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = FALSE)
#' plot_mean_violin(se)
plot_mean_violin <- function(
  se,
  assay_name = "diversity",
  sample_type_col = NULL
) {
    require_pkgs(c("ggplot2", "dplyr", "SummarizedExperiment", "tidyr"))
    long <- get_assay_long(
        se,
        assay_name = assay_name,
        value_name = "diversity",
        sample_type_col = sample_type_col
    )

    tmp <- as.data.frame(long)
    plot_df <- stats::aggregate(
        diversity ~ sample_type + Gene,
        data = tmp,
        FUN = function(x) mean(x, na.rm = TRUE)
    )
    colnames(plot_df)[colnames(plot_df) == "diversity"] <- "value"

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = sample_type, y = value, fill = sample_type)
    ) +
        ggplot2::geom_violin(alpha = 0.6) +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(x = "Samples", y = "Diversity") +
        ggplot2::scale_fill_viridis_d(name = "Group")
}


# Core MA plotting implementation documentation moved to internal block
#' Plot MA using Tsallis-based fold changes
#'
#' Wrapper around `plot_ma(..., type = "tsallis")` for convenience and
#' clearer API separation.
#'
#' @param x Data.frame from `calculate_difference()`.
#' @param sig_alpha Numeric significance threshold for adjusted p-values (default: 0.05).
#' @param x_label Optional x-axis label passed to `plot_ma`.
#' @param y_label Optional y-axis label passed to `plot_ma`.
#' @param title Optional plot title passed to `plot_ma`.
#' @param ... Additional arguments passed to `plot_ma()`.
#' @return A `ggplot2` object representing the MA plot.
#' @examples
#' x <- data.frame(genes = paste0("g", seq_len(5)), mean = runif(5), log2_fold_change = rnorm(5))
#' plot_ma_tsallis(x)
#' @export
plot_ma_tsallis <- function(x, sig_alpha = 0.05, x_label = NULL, y_label = NULL, title = NULL, ...) {
    title_use <- title %||% "Tsallis-based MA plot"
    x_label_use <- x_label %||% "mean_difference"
    y_label_use <- y_label %||% "Log10 fold-change of entropy"
    .plot_ma_core(x, fc_df = NULL, sig_alpha = sig_alpha, x_label = x_label_use, y_label = y_label_use, title = title_use)
}


#' Plot MA using expression/readcount-based fold changes
#'
#' Wrapper around `plot_ma(..., type = "expression")` that accepts a
#' `SummarizedExperiment` or precomputed fold-change `data.frame`.
#'
#' @param x Data.frame from `calculate_difference()`.
#' @param se A `SummarizedExperiment` or data.frame supplying readcounts or precomputed fold changes.
#' @param samples Optional sample grouping vector (passed to `plot_ma`).
#' @param control Control level name (passed to `plot_ma`).
#' @param fc_method Aggregation method for fold-change calculation (passed to `plot_ma`).
#' @param pseudocount Pseudocount added when computing log ratios (passed to `plot_ma`).
#' @param sig_alpha Numeric significance threshold for adjusted p-values (default: 0.05).
#' @param x_label Optional x-axis label passed to `plot_ma`.
#' @param y_label Optional y-axis label passed to `plot_ma`.
#' @param title Optional plot title passed to `plot_ma`.
#' @param ... Additional arguments passed to `plot_ma()`.
#' @return A `ggplot2` object representing the MA plot.
#' @examples
#' x <- data.frame(genes = paste0("g", seq_len(5)), mean = runif(5))
#' fc <- data.frame(genes = paste0("g", seq_len(5)), log2_fold_change = rnorm(5))
#' plot_ma_expression(x, se = fc)
#' @export
plot_ma_expression <- function(
  x,
  se,
  samples = NULL,
  control = NULL,
  fc_method = "median",
  pseudocount = 0,
  sig_alpha = 0.05,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  ...
) {
    title_use <- title %||% "Readcounts-based MA plot"
    x_label_use <- x_label %||% "Mean difference"
    y_label_use <- y_label %||% "Log10 fold-change of counts"
    plot_ma_expression_impl(
        x,
        se = se,
        samples = samples,
        control = control,
        fc_method = fc_method,
        pseudocount = pseudocount,
        sig_alpha = sig_alpha,
        x_label = x_label_use,
        y_label = y_label_use,
        title = title_use,
        ...
    )
}


# Core MA plotting implementation used by wrappers above. Accepts a
# differential results `x` (data.frame) and an optional `fc_df` with
# fold-changes (genes as rownames or a `genes` column). Returns a
# `ggplot` MA-plot.
#' Core MA plotting implementation (internal)
#'
#' This is an internal helper used by `plot_ma_tsallis()` and
#' `plot_ma_expression()`. It is documented here for developers but
#' is not exported.
#' @noRd
.plot_ma_core <- function(x,
                          fc_df = NULL,
                          diff_res = NULL,
                          sig_alpha = 0.05,
                          x_label = NULL,
                          y_label = NULL,
                          title = NULL,
                          ...) {
    require_pkgs(c("ggplot2"))

    df <- as.data.frame(x, stringsAsFactors = FALSE)
    # ensure a gene identifier column exists
    # Support both new "gene_id" and legacy "genes" column names
    if (!("genes" %in% colnames(df))) {
        if ("gene_id" %in% colnames(df)) {
            # Rename gene_id to genes for backward compatibility with plotting code
            df$genes <- df$gene_id
        } else if (!is.null(rownames(df))) {
            df$genes <- rownames(df)
        }
    }

    # If an external fc_df is provided, merge fold values into df
    if (!is.null(fc_df)) {
        fdf <- as.data.frame(fc_df, stringsAsFactors = FALSE)
        # Support both "gene_id" and "genes" column names
        if (!("genes" %in% colnames(fdf))) {
            if ("gene_id" %in% colnames(fdf)) {
                fdf$genes <- fdf$gene_id
            } else if (!is.null(rownames(fdf))) {
                fdf$genes <- rownames(fdf)
            }
        }
        if (!("log2_fold_change" %in% colnames(fdf))) {
            stop("Provided `fc_df` must contain 'log2_fold_change' column")
        }
        df <- merge(df, fdf[, c("genes", "log2_fold_change")], by = "genes", all.x = TRUE, suffixes = c("", ".fc"))
        # prefer fc_df fold values when present
        if ("log2_fold_change.fc" %in% colnames(df)) df$log2_fold_change <- ifelse(!is.na(df$log2_fold_change.fc), df$log2_fold_change.fc, df$log2_fold_change)
    }

    # Detect fold-change column
    fold_candidates <- c("log2_fold_change", "logFC", "fold", "estimate_interaction", "fold_change")
    fold_col <- intersect(fold_candidates, colnames(df))
    if (length(fold_col) == 0) stop("Could not find a fold-change column in input")
    fold_col <- fold_col[1]

    # Detect mean/average columns for x-axis
    mean_cols <- grep("_mean$|_median$", colnames(df), value = TRUE)
    if (length(mean_cols) >= 2) {
        # ensure the two chosen columns are consistent (both _mean or both _median)
        c1 <- mean_cols[1]
        c2 <- mean_cols[2]
        is_mean1 <- grepl("_mean$", c1)
        is_mean2 <- grepl("_mean$", c2)
        is_med1 <- grepl("_median$", c1)
        is_med2 <- grepl("_median$", c2)
        if (!((is_mean1 && is_mean2) || (is_med1 && is_med2))) {
            stop("Could not find two mean or two median columns")
        }
        xvals <- rowMeans(df[, mean_cols[seq_len(2)],
            drop = FALSE
        ], na.rm = TRUE)
        x_label <- x_label %||% paste0(mean_cols[1], " vs ", mean_cols[2])
    } else if (length(mean_cols) == 1) {
        xvals <- as.numeric(df[[mean_cols[1]]])
        x_label <- x_label %||% mean_cols[1]
    } else if ("mean" %in% colnames(df)) {
        xvals <- as.numeric(df$mean)
        x_label <- x_label %||% "Mean"
    } else {
        # fallback: use rank or index
        xvals <- seq_len(nrow(df))
        x_label <- x_label %||% "Index"
    }

    yvals <- as.numeric(df[[fold_col]])

    # p-value / adjusted p-value detection
    padj_candidates <- c("padj", "adjusted_p_values", "adj_p_value", "adj_p", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) padj_col[1] else NULL

    padj <- if (!is.null(padj_col)) as.numeric(df[[padj_col]]) else rep(1, length(yvals))
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < sig_alpha, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag, stringsAsFactors = FALSE)

    prep <- .tsenat_prepare_ma_plot_df(df, fold_col = fold_col, mean_cols = mean_cols, x_label = x_label, y_label = y_label)
    plot_df <- prep$plot_df
    x_label_formatted <- .tsenat_format_label(prep$x_label)
    y_label_raw <- prep$y_label %||% fold_col
    y_label_formatted <- .tsenat_format_label(y_label_raw)
    if (!is.null(y_label_formatted)) {
        y_label_formatted <- sub("\\blog2\\b", "log10", y_label_formatted, ignore.case = TRUE)
    }

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = significant)) +
        ggplot2::geom_point(alpha = 0.75, size = 3.2) +
        ggplot2::scale_color_manual(
            values = c("non-significant" = "grey40", "significant" = "firebrick3"),
            guide = "none"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(
            title = title %||% "MA plot: mean vs log10 fold-change",
            x = x_label_formatted,
            y = y_label_formatted
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain"),
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12)
        )

    p
}


# Implementation that computes fold-changes from expression/readcounts
# stored in a `SummarizedExperiment` (or matrix) and forwards to core plotter.
plot_ma_expression_impl <- function(
  x,
  se,
  samples = NULL,
  control = NULL,
  fc_method = "median",
  pseudocount = 0,
  sig_alpha = 0.05,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  ...
) {
    require_pkgs(c("SummarizedExperiment"))

    # extract/readcounts
    if (inherits(se, "SummarizedExperiment")) {
        counts <- get_readcounts_from_se(se)
        samples <- infer_samples_from_se(se, samples)
        if (is.null(samples)) stop("Could not infer 'samples' from SummarizedExperiment; provide `samples`")
        control <- validate_control_in_samples(control, samples)

        # attempt to map transcripts -> genes and aggregate counts per gene
        tx2g <- get_tx2gene_from_se(se, readcounts_mat = counts)
        if (!is.null(tx2g) && tx2g$type == "vector") {
            mapping <- tx2g$mapping
            # ensure mapping length matches rows
            if (length(mapping) == nrow(counts)) {
                # aggregate transcript-level counts to gene-level using rowsum
                agg <- rowsum(counts, group = mapping)
                counts_gene <- as.matrix(agg)
            } else {
                counts_gene <- counts
            }
        } else {
            counts_gene <- counts
        }

        # compute fold-changes using calculate_fc (aggregates per-group)
        fc_res <- calculate_fc(
            counts_gene,
            samples,
            control,
            method = fc_method,
            pseudocount = pseudocount
        )
        # ensure genes column exists
        if (is.null(rownames(fc_res))) {
            fc_res$genes <- seq_len(nrow(fc_res))
        } else {
            fc_res$genes <- rownames(fc_res)
        }
        return(.plot_ma_core(
            x,
            fc_df = fc_res,
            sig_alpha = sig_alpha,
            x_label = x_label,
            y_label = y_label,
            title = title,
            ...
        ))
    }

    # If se is provided as a matrix/data.frame of precomputed fold changes
    if (is.matrix(se) || is.data.frame(se)) {
        fc_res <- as.data.frame(se, stringsAsFactors = FALSE)
        if (!("log2_fold_change" %in% colnames(fc_res))) stop("`se` data.frame must contain 'log2_fold_change' column when providing precomputed fold changes")
        if (!("genes" %in% colnames(fc_res)) && !is.null(rownames(fc_res))) fc_res$genes <- rownames(fc_res)
        return(.plot_ma_core(x, fc_df = fc_res, sig_alpha = sig_alpha, x_label = x_label, y_label = y_label, title = title, ...))
    }

    stop("Unsupported 'se' argument for plot_ma_expression_impl")
}


#' Plot Tsallis Q-curve Profile
#'
#' Visualize q-curve showing Tsallis entropy across diversity scales for each sample group.
#' Displays median entropy with IQR ribbons for comparison between groups.
#'
#' @param se SummarizedExperiment from calculate_diversity().
#' @param assay_name Character. Assay name (default "diversity").
#' @param sample_type_col Character. Column in colData(se) with sample types (default "sample_type").
#'
#' @return A ggplot object with q-curves and IQR ribbons for each group.
#'
#' @details Computes median Tsallis entropy and interquartile range (IQR) at each q-value
#' for each sample group, displayed as median line with IQR ribbon.
#'
#' @export
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:40, -1, drop = FALSE])
#' gs <- readcounts[1:40, 1]
#' se <- calculate_diversity(rc, gs,
#'     q = seq(0.01, 0.1, by = 0.03), norm = FALSE
#' )
#' p <- plot_tsallis_q_curve(se)
#' p
plot_tsallis_q_curve <- function(
  se,
  assay_name = "diversity",
  sample_type_col = "sample_type"
) {
    # SE-first API: require a SummarizedExperiment with per-column sample
    # type mapping in `colData(se)[, sample_type_col]` (or allow a single
    # group dataset where `sample_type` is omitted).
    if (inherits(se, "SummarizedExperiment")) {
        require_pkgs(c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment"))
        long <- prepare_tsallis_long(se, assay_name = assay_name, sample_type_col = sample_type_col)
        y_label <- "Tsallis entropy (S_q)"
        if (nrow(long) == 0) stop("No tsallis values found in SummarizedExperiment")
        # Ensure q is numeric; handle case where it might be a factor or character
        long$q <- as.numeric(as.character(long$q))
        
        # Compute median and IQR at each q-value for each group
        stats_df <- dplyr::summarise(dplyr::group_by(long, group, q),
            median = median(tsallis, na.rm = TRUE),
            IQR = stats::IQR(tsallis, na.rm = TRUE), .groups = "drop"
        )
        
        # Create plot with IQR ribbons
        p <- ggplot2::ggplot(
            stats_df,
            ggplot2::aes(x = q, y = median, color = group, fill = group)
        ) +
            ggplot2::geom_line(linewidth = 1.3) +
            ggplot2::geom_ribbon(
                ggplot2::aes(ymin = median - IQR / 2, ymax = median + IQR / 2),
                alpha = 0.2, color = NA
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::labs(
                title = "Tsallis q-curve: median ± IQR",
                x = "q value",
                y = y_label,
                color = "Group",
                fill = "Group"
            ) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 16))
        
        # use default discrete ggplot2 colours (not viridis)
        p <- p + ggplot2::scale_color_discrete(name = "Group") +
            ggplot2::scale_fill_discrete(name = "Group")
        # If there is only a single group present, hide the legend/Group label
        if (length(unique(long$group)) == 1) {
            p <- p + ggplot2::theme(legend.position = "none")
        }

        return(p)
    }

    # Matrix/data.frame input is no longer supported for this plot function.
    stop(
        "plot_tsallis_q_curve requires a SummarizedExperiment from calculate_diversity."
    )
}
#' Violin plot of Tsallis entropy for multiple q values

#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with
#' multiple q values (column names contain `_q=`).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @return A `ggplot` violin plot object faceted/colored by group and q.
#' @export
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:20, -1, drop = FALSE])
#' gs <- readcounts[1:20, 1]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
#' plot_tsallis_violin_multq(se)
plot_tsallis_violin_multq <- function(se, assay_name = "diversity") {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    long <- prepare_tsallis_long(se, assay_name = assay_name)
    # Ensure q is numeric first
    long$q <- as.numeric(as.character(long$q))
    # Convert q to factor for proper categorical plotting in violin plot
    long$q_label <- factor(paste0("q = ", long$q), 
                           levels = paste0("q = ", sort(unique(long$q))))

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = q_label, y = tsallis, fill = group)
    ) +
        ggplot2::geom_violin(
            alpha = 0.5, width = 0.9,
            position = ggplot2::position_dodge(width = 0.8)
        ) +
        ggplot2::geom_boxplot(
            width = 0.15,
            position = ggplot2::position_dodge(width = 0.8),
            outlier.shape = NA
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::labs(
            title = "Violin plot: Tsallis entropy distribution across multiple q values",
            x = "q value",
            y = "Tsallis entropy",
            fill = "Group"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(
            hjust = 0.5, size = 16, face = "plain",
            margin = ggplot2::margin(b = 10)
        ))
}
#' Density plot of Tsallis entropy for multiple q values

#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with
#' multiple q values (column names contain `_q=`).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @return A `ggplot` density plot object faceted by q and colored by group.
#' @export
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:20, -1, drop = FALSE])
#' gs <- readcounts[1:20, 1]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = FALSE)
#' plot_tsallis_density_multq(se)
plot_tsallis_density_multq <- function(se, assay_name = "diversity") {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    long <- prepare_tsallis_long(se, assay_name = assay_name)
    # Ensure q is numeric; handle case where it might be a factor or character
    long$q <- as.numeric(as.character(long$q))
    # Create a proper factor with levels sorted numerically for consistent faceting
    long$q_label <- paste0("q = ", long$q)
    long$q_label <- factor(long$q_label, levels = paste0("q = ", sort(unique(long$q))))

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = tsallis, color = group, fill = group)
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_wrap(~q_label, scales = "free_y") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::scale_color_discrete(name = "Group") +
        ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::labs(
            title = "Density plot: Tsallis entropy distribution across multiple q values",
            x = "Tsallis entropy",
            y = "Density",
            color = "Group",
            fill = "Group"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(
            hjust = 0.5, size = 16, face = "plain",
            margin = ggplot2::margin(b = 10)
        ))
}


#' Volcano plot for differential results
#'
#' Create a volcano plot showing fold-change (x-axis) versus adjusted
#' p-value significance (y-axis). The function auto-detects a suitable x-axis
#' column if one is not provided and expects an adjusted p-value column for
#' significance coloring.
#'
#' @param diff_df Data.frame with differential expression results.
#'   Should contain p-values and optionally fold-change columns.
#' @param x_col Optional column name for the x-axis. If `NULL`, the function
#'   will try to auto-detect a suitable numeric column (excluding p-values).
#' @param padj_col Adjusted p-value column name (default: "padj").
#' @param label_thresh Fold-change threshold used to annotate points (default: 0.1).
#' @param sig_alpha Adjusted p-value cutoff for significance (default: 0.05).
#' @param top_n Number of top significant genes to label (default: 5).
#' @param title Optional plot title; if `NULL` a default title is used.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' df <- data.frame(
#'     gene = paste0("g", seq_len(10)),
#'     mean_difference = runif(10),
#'     padj = runif(10)
#' )
#' # plot_volcano(df, x_col = "mean_difference", padj_col = "padj")
plot_volcano <- function(
  diff_df,
  x_col = NULL,
  padj_col = "padj",
  label_thresh = 0.1,
  sig_alpha = 0.05,
  top_n = 5,
  title = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required")
    }

    prep_volcano <- .tsenat_prepare_volcano_df(diff_df = diff_df, x_col = x_col, padj_col = padj_col, label_thresh = label_thresh, sig_alpha = sig_alpha, title = title)
    df <- prep_volcano$df
    x_col <- prep_volcano$x_col
    padj_col <- prep_volcano$padj_col
    x_label_formatted <- prep_volcano$x_label_formatted
    padj_label_formatted <- prep_volcano$padj_label_formatted
    title_use <- prep_volcano$title_use

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = xval, y = -log10(padj), color = significant)
    ) +
        ggplot2::geom_point(alpha = 0.75, size = 3.4) +
        ggplot2::scale_color_manual(
            values = c("non-significant" = "black", "significant" = "red"),
            guide = "none"
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(sig_alpha),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-label_thresh, label_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(
            title = title_use,
            x = x_label_formatted,
            y = paste0("-Log10(", padj_label_formatted, ")")
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain", margin = ggplot2::margin(b = 4)),
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12)
        )

    p
}


#' Internal helper to compute fill limits across multiple genes (not exported)
#' @noRd
.compute_transcript_fill_limits <- function(genes, mapping, counts, samples, top_n, agg_fun, pseudocount) {
    mins <- maxs <- c()
    for (g in genes) {
        txs <- mapping$Transcript[mapping$Gen == g]
        txs <- intersect(txs, rownames(counts))
        if (length(txs) == 0) next
        if (!is.null(top_n)) txs <- head(txs, top_n)
        mat <- counts[txs, , drop = FALSE]
        df_all <- as.data.frame(mat)
        df_all$tx <- rownames(mat)
        df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
        df_long$group <- rep(samples, times = length(txs))
        df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
        df_summary$log2expr <- log2(df_summary$expr + pseudocount)
        mins <- c(mins, min(df_summary$log2expr, na.rm = TRUE))
        maxs <- c(maxs, max(df_summary$log2expr, na.rm = TRUE))
    }
    if (length(mins) == 0) stop("No transcripts found for provided genes")
    c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))
}

#' Internal helper to draw grid layout with title, plots, and legend using base grid
#' @noRd
.draw_transcript_grid <- function(grobs, title, legend_grob, ncol, heights, to_file = NULL) {
    # If no output file is provided and no graphics device is open, render to a
    # temporary pdf device so that plotting in non-interactive sessions does not
    # create `Rplots.pdf` in the working directory.
    temp_dev <- FALSE
    # Only open a temporary PDF device when:
    #  - caller did not supply an output file (`to_file` is NULL),
    #  - the session is non-interactive, and
    #  - no graphics device is currently open (dev.cur() == 1)
    if (is.null(to_file) && !interactive() && grDevices::dev.cur() == 1L) {
        tmp <- tempfile("TSENAT_plot_", fileext = ".pdf")
        grDevices::pdf(tmp)
        temp_dev <- TRUE
        # Ensure device is closed and temporary file removed on exit
        on.exit(
            {
                try(grDevices::dev.off(), silent = TRUE)
                if (file.exists(tmp)) unlink(tmp)
            },
            add = TRUE
        )
    }

    # Calculate number of rows needed for plots
    nrow_plots <- ceiling(length(grobs) / ncol)
    nrow_total <- 2 + nrow_plots  # title + plot rows + legend

    grid::grid.newpage()
    grid::pushViewport(
        grid::viewport(layout = grid::grid.layout(nrow_total, ncol, heights = heights))
    )
    # Title row
    vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(ncol))
    grid::pushViewport(vp_title)
    grid::grid.text(title, x = 0.5, gp = grid::gpar(fontsize = 12))
    grid::upViewport()
    # Plot rows
    for (i in seq_along(grobs)) {
        plot_row_idx <- ((i - 1) %/% ncol) + 2
        plot_col_idx <- ((i - 1) %% ncol) + 1
        vp <- grid::viewport(layout.pos.row = plot_row_idx, layout.pos.col = plot_col_idx)
        grid::pushViewport(vp)
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }
    # Legend row
    if (!is.null(legend_grob)) {
        vp_leg <- grid::viewport(layout.pos.row = nrow_total, layout.pos.col = seq_len(ncol))
        grid::pushViewport(vp_leg)
        grid::grid.draw(legend_grob)
        grid::upViewport()
    }
    grid::upViewport()
    # If caller supplied a file (caller likely opened a device), close it here.
    if (!is.null(to_file)) grDevices::dev.off()
    invisible(NULL)
}

## Internal helpers for `plot_top_transcripts` refactor
## Create per-gene plot and combine multiple gene plots into final output
.ptt_make_plot_for_gene <- function(gene_single, mapping, counts, samples, top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL) {
    require_pkgs(c("ggplot2", "tidyr"))
    built <- .ptt_build_tx_long(gene_single, mapping, counts, samples, top_n)
    df_summary <- .ptt_aggregate_df_long(built$df_long, agg_fun, pseudocount)
    .ptt_build_plot_from_summary(df_summary, agg_label_unique, fill_limits)
}

.ptt_combine_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
    require_pkgs(c("ggplot2"))
    # Allow callers to pass a single character second argument as the
    # `agg_label_unique` for convenience (legacy test call patterns).
    if (is.null(agg_label_unique) && !is.null(output_file) && is.character(output_file) && length(output_file) == 1) {
        agg_label_unique <- output_file
        output_file <- NULL
    }
    if (requireNamespace("patchwork", quietly = TRUE)) {
        .ptt_combine_patchwork(plots, agg_label_unique)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        .ptt_combine_cowplot(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        .ptt_combine_grid(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    }
}

## Prepare and validate inputs for `plot_top_transcripts`
.ptt_prepare_inputs <- function(counts, readcounts = NULL, samples = NULL, coldata = NULL, sample_type_col = "sample_type", tx2gene = NULL, res = NULL, top_n = NULL, pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance", "iqr")) {
    # handle selecting genes from `res` is left to caller; this function focuses
    # on normalizing counts, samples and tx2gene mapping and preparing agg functions
    if (inherits(counts, "SummarizedExperiment")) {
        require_pkgs(c("SummarizedExperiment", "S4Vectors"))
        se <- counts
        counts_mat <- get_readcounts_from_se(se, readcounts)
        counts <- as.matrix(counts_mat)
        samples <- infer_samples_from_se(se, samples, sample_type_col = sample_type_col)

        if (is.null(tx2gene)) {
            txres <- get_tx2gene_from_se(se, counts)
            if (!is.null(txres) && !is.null(txres$mapping)) {
                mapping <- data.frame(Transcript = rownames(counts), Gen = as.character(txres$mapping), stringsAsFactors = FALSE)
                tx2gene <- mapping
            }
        }
    }

    if (!is.matrix(counts) && !is.data.frame(counts)) stop("`counts` must be a matrix or data.frame with transcripts as rownames")
    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) stop("`counts` must have rownames corresponding to transcript identifiers")

    # derive samples from coldata if needed
    if (is.null(samples)) {
        if (!is.null(coldata)) {
            if (is.character(coldata) && length(coldata) == 1) {
                if (!file.exists(coldata)) stop("coldata file not found: ", coldata)
                cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
            } else if (is.data.frame(coldata)) {
                cdf <- coldata
            } else {
                stop("`coldata` must be a data.frame or path to a tab-delimited file")
            }

            if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
                samples <- as.character(cdf[colnames(counts), sample_type_col])
            } else {
                sample_id_cols <- c("sample", "Sample", "sample_id", "id")
                sid <- intersect(sample_id_cols, colnames(cdf))
                if (length(sid) > 0) {
                    sid <- sid[1]
                    if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) stop("coldata sample id column does not match column names of counts")
                    row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
                    samples <- as.character(cdf[[sample_type_col]][row_ix])
                } else {
                    stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
                }
            }
        } else {
            stop("Either 'samples' or 'coldata' must be provided to determine sample groups")
        }
    }

    # normalize tx2gene mapping
    if (is.null(tx2gene)) stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) stop("tx2gene file not found: ", tx2gene)
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) stop("tx2gene must have columns 'Transcript' and 'Gen'")

    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for plotting")

    if (!is.null(samples) && length(samples) != ncol(counts)) stop("Length of `samples` must equal number of columns in `counts`")

    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice,
        median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE),
        variance = function(x) stats::var(x, na.rm = TRUE),
        iqr = function(x) stats::IQR(x, na.rm = TRUE)
    )
    agg_label_metric <- if (metric_choice == "iqr") "IQR" else metric_choice
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label

    list(counts = counts, samples = samples, mapping = mapping, metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique, top_n = top_n, pseudocount = pseudocount, output_file = output_file)
}

#' Plot top transcripts for a gene
#' @param counts Matrix or data.frame of transcript counts. Rows are transcripts and columns are samples.
#' @param readcounts Optional matrix or data.frame of raw read counts. Used for transcript-level quantification if provided.
#' @param gene Character; gene symbol to inspect.
#' @param samples Character vector of sample group labels (length = ncol(counts)).
#' @param coldata Optional data.frame or file path containing sample metadata. Used to infer sample groups if `samples` is not provided.
#' @param sample_type_col Character; column name in `coldata` or `SummarizedExperiment` colData to use for sample grouping. Default is "sample_type".
#' @param tx2gene Path or data.frame mapping transcripts to genes. Must contain columns `Transcript` and `Gen`.
#' @param res Optional result data.frame from a differential analysis. If provided and `gene` is NULL, top genes are selected by adjusted p-value.
#' @param top_n Integer number of transcripts to show (default = 3). Use NULL to plot all transcripts for the gene.
#' @param pseudocount Numeric pseudocount added before log2 (default = 1e-6) to avoid division by zero.
#' @param output_file Optional file path to save the plot. If `NULL`, the `ggplot` object is returned.
#' @param metric Aggregation metric used to summarize transcript expression per group when plotting. One of c("median", "mean", "variance", "iqr"). Use "iqr" to compute the interquartile range. Defaults to "median".
#' @return If \code{output_file} is \code{NULL}, returns a \code{ggplot} object. Otherwise, the plot is saved to the specified file and the function returns \code{NULL} invisibly.
#' @examples
#' tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
#' rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
#' colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
#' tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 2), stringsAsFactors = FALSE)
#' samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
#' plot_top_transcripts(tx_counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2gene, top_n = 2)
#' @export
plot_top_transcripts <- function(
  counts,
  readcounts = NULL, # Optional matrix or data.frame of raw read counts. Used for transcript-level quantification if provided.
  gene = NULL,
  samples = NULL,
  coldata = NULL, # Optional data.frame or file path containing sample metadata. Used to infer sample groups if `samples` is not provided.
  sample_type_col = "sample_type", # Column name in `coldata` or `SummarizedExperiment` colData to use for sample grouping. Default is "sample_type".
  tx2gene = NULL,
  res = NULL, # Optional result data.frame from a differential analysis. If provided and `gene` is NULL, top genes are selected by adjusted p-value.
  top_n = 3,
  pseudocount = 1e-6,
  output_file = NULL,
  metric = c("median", "mean", "variance", "iqr")
) {
    # If counts is a SummarizedExperiment and res is provided, extract components intelligently
    if (inherits(counts, "SummarizedExperiment") && !is.null(res)) {
        se <- counts
        
        # Extract rowData to identify gene name column
        rd <- SummarizedExperiment::rowData(se)
        gene_names_col <- if ("gene_names" %in% colnames(rd)) "gene_names" else if ("gene_name" %in% colnames(rd)) "gene_name" else "genes"
        se_gene_names <- unique(as.character(rd[[gene_names_col]]))
        se_gene_names <- se_gene_names[!is.na(se_gene_names) & nzchar(se_gene_names)]
        
        # Filter results to only genes present in filtered SE
        res_gene_col <- if ("genes" %in% colnames(res)) "genes" else if ("gene" %in% colnames(res)) "gene" else "gene_id"
        if (!(res_gene_col %in% colnames(res))) {
            stop("Provided 'res' must contain a 'genes', 'gene', or 'gene_id' column", call. = FALSE)
        }
        res_filtered <- res[as.character(res[[res_gene_col]]) %in% se_gene_names, ]
        res_filtered <- res_filtered[!is.na(res_filtered[[res_gene_col]]), ]
        
        if (nrow(res_filtered) > 0) {
            # Prepare parameters for the plot function
            if (is.null(gene)) {
                # Sort genes by adjusted p-value
                if ("padj" %in% colnames(res_filtered)) {
                    ord <- order(res_filtered$padj, na.last = NA)
                } else if ("adjusted_p_values" %in% colnames(res_filtered)) {
                    ord <- order(res_filtered$adjusted_p_values, na.last = NA)
                } else if ("pvalue" %in% colnames(res_filtered)) {
                    ord <- order(res_filtered$pvalue, na.last = NA)
                } else if ("raw_p_values" %in% colnames(res_filtered)) {
                    ord <- order(res_filtered$raw_p_values, na.last = NA)
                } else {
                    ord <- seq_len(nrow(res_filtered))
                }
                gene <- head(as.character(res_filtered[[res_gene_col]][ord]), top_n)
            }
            
            # Build tx2gene mapping
            if (is.null(tx2gene)) {
                tx2gene <- data.frame(
                    Transcript = rownames(se),
                    Gen = as.character(rd[[gene_names_col]]),
                    stringsAsFactors = FALSE
                )
            }
            
            # Extract counts matrix
            counts <- SummarizedExperiment::assay(se)
            
            # Extract samples if not provided
            if (is.null(samples) && is.null(coldata)) {
                cd <- SummarizedExperiment::colData(se)
                if (!is.null(cd) && nrow(cd) > 0) {
                    sample_cols <- c(sample_type_col, "sample_type", "condition", "group", "sample_group", "class", "status", "phenotype")
                    # Remove NULL values from sample_cols
                    sample_cols <- sample_cols[!is.na(sample_cols)]
                    for (col in sample_cols) {
                        if (!is.na(col) && col %in% colnames(cd)) {
                            samples <- as.character(cd[[col]])
                            break
                        }
                    }
                }
            }
            
            # Clear res after extracting genes to avoid downstream re-filtering
            res <- NULL
        }
    }
    
    # If `gene` is not provided, select top genes from `res` using `top_n`.
    if (is.null(gene)) {
        gene <- .ptt_select_genes_from_res(res, top_n)
    }
    per_gene_top_n <- top_n

    ## Prepare inputs and normalization via helper
    prep <- .ptt_prepare_inputs(counts = counts, readcounts = readcounts, samples = samples, coldata = coldata, sample_type_col = sample_type_col, tx2gene = tx2gene, res = res, top_n = per_gene_top_n, pseudocount = pseudocount, output_file = output_file, metric = metric)

    # if prep returned without gene selection, caller will check `gene`
    counts <- prep$counts
    samples <- prep$samples
    mapping <- prep$mapping
    metric_choice <- prep$metric_choice
    agg_fun <- prep$agg_fun
    agg_label_unique <- prep$agg_label_unique
    top_n <- prep$top_n
    pseudocount <- prep$pseudocount
    output_file <- prep$output_file

    make_plot_for_gene <- function(gene_single, fill_limits = NULL) {
        .ptt_make_plot_for_gene(gene_single, mapping, counts, samples, top_n, agg_fun, pseudocount, agg_label_unique, fill_limits)
    }

    # Produce plots (single or multiple). Do not save inside helper - save once
    # below.
    if (length(gene) > 1) {
        fill_limits <- .compute_transcript_fill_limits(gene, mapping, counts, samples, top_n, agg_fun, pseudocount)

        plots <- lapply(seq_along(gene), function(i) {
            gname <- gene[i]
            pp <- make_plot_for_gene(gname, fill_limits = fill_limits)
            per_gene_title <- if (!is.na(gname) && nzchar(as.character(gname))) as.character(gname) else ""
            pp <- pp + ggplot2::labs(title = per_gene_title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16))
            pp
        })

        result_plot <- .ptt_combine_plots(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        result_plot <- make_plot_for_gene(gene)
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    } else {
        result_plot
    }
}

#' Plot GAM q-curves for top genes identified by FPCA/GAM interaction tests
#'
#' Visualizes smooth q-curve profiles (GAM fits) for selected genes from
#' `calculate_lm_interaction()` results. Useful for understanding which q-ranges (rare
#' vs. dominant isoforms) drive significant PC differences between groups.
#'
#' @param se A `SummarizedExperiment` with q-sequence diversity values
#'   (multiple q values per sample). Typically output from `calculate_diversity()`
#'   with multiple q (e.g., q = seq(0.1, 2, by = 0.1)).
#' @param lm_res A `data.frame` from `calculate_lm_interaction()` with columns
#'   `gene`, `p_interaction`, and `adj_p_interaction`. Can be from method="fpca"
#'   or method="gam".
#' @param sample_type_col Column name in `colData(se)` specifying group
#'   assignments for samples (e.g., "sample_type").
#' @param genes Optional character vector of specific gene names to plot. If provided,
#'   these genes are plotted directly regardless of significance or n_top. If NULL
#'   (default), the top n_top significant genes are selected.
#' @param n_top Number of top genes (by adjusted p-value) to plot (default: 6).
#'   Only used if genes = NULL.
#' @param sig_alpha Significance threshold for adjusted p-values (default: 0.05).
#'   Only used if genes = NULL; filters lm_res to significant genes before selecting top n.
#' @param assay_name Name of the assay in `se` to extract (default: "diversity").
#' @param palette Color palette for group separation (default: "Set1").
#'
#' @return A list of `ggplot` objects, one per selected gene, showing GAM-fitted
#'   q-curves colored by sample group. If only one gene is requested, returns a 
#'   single `ggplot` object.
#'
#' @details
#' For each selected gene, this function:
#' 1. Extracts per-sample entropy values across all q values
#' 2. Fits GAM models: entropy ~ s(q, k=...) independently for each group
#' 3. Generates smooth predictions for visualization
#' 4. Overlays predicted curves for each group with a distinct color
#'
#' This complements FPCA by providing interpretable visualization of empirical
#' q-curve shape differences that drive PC-level significance.
#'
#' @examples
#' data("readcounts", package = "TSENAT")
#' rc <- as.matrix(readcounts[1:100, -1, drop = FALSE])
#' gs <- readcounts[1:100, 1]
#' # Create q-sequence
#' se <- calculate_diversity(rc, gs, q = seq(0.5, 2, by = 0.25), norm = TRUE)
#' # Add sample type to colData
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'     sample_type = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'     row.names = colnames(se)
#' )
#' # Run FPCA to identify significant genes
#' lm_res <- calculate_lm_interaction(se, sample_type_col = "sample_type", method = "fpca")
#' # Plot GAM curves for top 3 genes
#' if (nrow(lm_res) > 0) {
#'   plot_lm_interaction_gam(se, lm_res, sample_type_col = "sample_type", n_top = 3)
#' }
#' # Or plot specific genes of interest
#' if (nrow(lm_res) > 0) {
#'   plot_lm_interaction_gam(se, lm_res, sample_type_col = "sample_type", 
#'                            genes = c("g1", "g2", "g5"))
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal scale_color_brewer
plot_lm_interaction_gam <- function(se, lm_res, sample_type_col, genes = NULL, n_top = 6,
    sig_alpha = 0.05, assay_name = "diversity", palette = "Set1") {

    require_pkgs(c("ggplot2", "mgcv", "SummarizedExperiment", "dplyr", "tidyr"))

    # Validate inputs
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }

    if (!is.data.frame(lm_res) || !("gene" %in% colnames(lm_res))) {
        stop("lm_res must be a data.frame with 'gene' column from calculate_lm_interaction()",
            call. = FALSE)
    }

    if (nrow(lm_res) == 0) {
        stop("lm_res has no rows; calculate_lm_interaction() returned no genes", call. = FALSE)
    }

    # Extract assay matrix and colData
    mat <- SummarizedExperiment::assay(se, assay_name)
    cdata <- SummarizedExperiment::colData(se)

    if (!sample_type_col %in% colnames(cdata)) {
        stop(sprintf("Column '%s' not found in colData(se)", sample_type_col), call. = FALSE)
    }

    groups <- cdata[[sample_type_col]]

    # Determine which genes to plot
    if (!is.null(genes)) {
        # User provided specific genes
        if (!is.character(genes)) {
            stop("genes must be a character vector of gene names", call. = FALSE)
        }
        top_genes <- genes
    } else {
        # Select top genes by p-value (original logic)
        if ("adj_p_interaction" %in% colnames(lm_res)) {
            sig_mask <- lm_res$adj_p_interaction <= sig_alpha
        } else if ("p_interaction" %in% colnames(lm_res)) {
            sig_mask <- lm_res$p_interaction <= sig_alpha
        } else {
            stop("lm_res must contain 'adj_p_interaction' or 'p_interaction' column", call. = FALSE)
        }

        sig_genes <- lm_res[sig_mask, , drop = FALSE]

        if (nrow(sig_genes) == 0) {
            warning(sprintf("No genes significant at alpha = %g", sig_alpha), call. = FALSE)
            return(NULL)
        }

        # Select top n
        top_genes <- sig_genes$gene[seq_len(min(n_top, nrow(sig_genes)))]
    }

    # Helper to build data.frame for a single gene
    make_gam_plot <- function(g) {
        if (!(g %in% rownames(mat))) {
            warning(sprintf("Gene '%s' not found in assay", g), call. = FALSE)
            return(NULL)
        }

        # Extract data
        gene_vals <- mat[g, ]
        sample_names <- colnames(mat)

        # Get metadata: need q-values per sample
        # Assume colnames have structure or we extract from rownames of mat
        # For typical TSENAT workflow, q-values are stored in metadata or need to be inferred
        # For now, we construct a long-format data.frame per sample
        
        # Try to extract q-values from colnames or metadata
        # If not available, try to infer from number of samples
        q_vals <- NULL
        
        # Check if q-values are in colData under a standard name
        for (q_col in c("q", "q_value", "q_values")) {
            if (q_col %in% colnames(cdata)) {
                q_vals <- cdata[[q_col]]
                break
            }
        }
        
        # If q values not found, infer from data structure
        # (assuming rows are unique q values per sample in order)
        if (is.null(q_vals)) {
            # For the typical build_se output where assays are concatenated
            # across q values, we need to reconstruct q-values
            # This is a heuristic: if we have metadata about diversity calculation
            n_samples <- length(sample_names)
            n_q <- length(gene_vals) / n_samples
            
            if (n_q != floor(n_q)) {
                warning(sprintf("Cannot infer q-values for gene '%s'; inconsistent dimensions", g),
                    call. = FALSE)
                return(NULL)
            }
            
            # Reconstruct assuming lexicographic ordering (q varies fastest or slowest)
            # Standard TSENAT: rows are samples, columns are q-values within each sample
            # For multi-q assays: columns are ordered by q within each sample
            q_vals <- rep(seq_len(as.integer(n_q)), each = n_samples)
        }
        
        # Build long-format data.frame
        plot_df <- data.frame(
            sample = rep(sample_names, length.out = length(gene_vals)),
            group = rep(groups, length.out = length(gene_vals)),
            q = q_vals,
            entropy = as.numeric(gene_vals),
            stringsAsFactors = FALSE
        )

        # Remove NA entries
        plot_df <- plot_df[!is.na(plot_df$entropy), , drop = FALSE]

        if (nrow(plot_df) == 0) {
            warning(sprintf("No valid data for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        # Fit GAM per group
        unique_groups <- unique(plot_df$group)

        if (length(unique_groups) < 2) {
            warning(sprintf("Less than 2 groups for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        # Generate prediction grid
        q_range <- range(plot_df$q, na.rm = TRUE)
        pred_q <- seq(q_range[1], q_range[2], length.out = 100)

        # Fit GAM and predict for each group
        pred_list <- list()
        for (gr in unique_groups) {
            subset_data <- subset(plot_df, group == gr)
            if (nrow(subset_data) < 3) {
                next
            }

            tryCatch(
                {
                    # Fit GAM with adaptive k (bases)
                    k <- min(10, max(2, round(nrow(subset_data) / 2)))
                    gam_fit <- mgcv::gam(entropy ~ s(q, k = k), data = subset_data)
                    
                    # Predict
                    pred_data <- data.frame(q = pred_q)
                    pred_vals <- stats::predict(gam_fit, newdata = pred_data, se.fit = TRUE)
                    
                    pred_list[[gr]] <- data.frame(
                        group = gr,
                        q = pred_q,
                        entropy_fit = pred_vals$fit,
                        se = pred_vals$se.fit,
                        stringsAsFactors = FALSE
                    )
                },
                error = function(e) {
                    warning(sprintf("GAM fit failed for gene '%s' group '%s': %s", g, gr, e$message),
                        call. = FALSE)
                }
            )
        }

        if (length(pred_list) == 0) {
            warning(sprintf("No GAM fits succeeded for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        pred_df <- do.call(rbind, pred_list)

        # Create plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = entropy, color = group)) +
            ggplot2::geom_point(alpha = 0.5, size = 2) +
            ggplot2::geom_line(data = pred_df, ggplot2::aes(x = q, y = entropy_fit, color = group,
                linetype = "GAM fit"), linewidth = 1, alpha = 0.9) +
            ggplot2::scale_color_brewer(palette = palette, name = sample_type_col) +
            ggplot2::scale_linetype_manual(values = c("GAM fit" = 1), name = "") +
            ggplot2::labs(
                x = "q parameter",
                y = "Tsallis entropy",
                title = sprintf("Gene: %s", g),
                subtitle = sprintf("GAM q-curve fit by group")
            ) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold"),
                legend.position = "bottom"
            )

        return(p)
    }

    # Generate plots for top genes
    plots <- list()
    for (g in top_genes) {
        p <- make_gam_plot(g)
        if (!is.null(p)) {
            plots[[g]] <- p
        }
    }

    if (length(plots) == 0) {
        warning("No valid plots generated", call. = FALSE)
        return(NULL)
    }

    # Return as list (or single plot if n_top == 1)
    if (length(plots) == 1) {
        return(plots[[1]])
    } else {
        return(plots)
    }
}

#' Bootstrap Confidence Intervals for Q-curve with Statistical Testing
#'
#' Compute bootstrap confidence bands and perform statistical testing to identify
#' q-ranges where groups significantly differ in Tsallis entropy. Supports both
#' pointwise and simultaneous confidence intervals.
#'
#' @param se A \code{SummarizedExperiment} returned by \code{calculate_diversity}
#'   with multiple q values (column names contain \code{_q=}).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Column name in colData specifying group assignments
#'   (default: "sample_type").
#' @param n_bootstrap Number of bootstrap resamples (default: 1000).
#' @param ci_level Confidence level for intervals (default: 0.95).
#' @param ci_type Type of confidence interval: "pointwise" (independent at each q)
#'   or "simultaneous" (controls family-wise error across all q). Default: "pointwise".
#' @param test_method Statistical test for group differences: "wilcox" (Wilcoxon rank-sum)
#'   or "ttest" (t-test). Default: "wilcox".
#' @param alpha Significance level for identifying different q-ranges (default: 0.05).
#' @param method Bootstrap method: "bca" (bias-corrected & accelerated), "percentile",
#'   or "normal" approximation. Default: "percentile".
#'
#' @return A list with class "qcurve_bootstrap" containing:
#'   \describe{
#'     \item{plot_data}{Data frame with columns: q, median, ci_lower, ci_upper, group,
#'       test_pvalue, significant}
#'     \item{significant_qranges}{Data frame identifying q-ranges where groups differ significantly}
#'     \item{metadata}{List with parameters: n_bootstrap, ci_level, ci_type, test_method}
#'     \item{ggplot}{A \code{ggplot} object with confidence bands and significance annotations}
#'   }
#'
#' @details
#' **Bootstrap Procedure:**
#' For each group and q-value:
#' 1. Resample genome with replacement (whole genes, not individual transcripts)
#' 2. Compute median Tsallis entropy across all genes in resample
#' 3. Collect n_bootstrap replicates
#' 4. Compute confidence intervals from percentiles (or BCa adjustment)
#'
#' **Statistical Testing:**
#' At each q-value, test whether group distributions differ using:
#' - Wilcoxon rank-sum test (non-parametric, recommended)
#' - Welch's t-test (parametric alternative)
#'
#' **Confidence Interval Types:**
#' - **Pointwise CI**: 95% at each q independently. Stricter for single comparisons,
#'   but less conservative when examining many q values.
#' - **Simultaneous CI** (not yet implemented): Controls family-wise error across all q,
#'   useful for interpreting full q-curve differences. Would use Bonferroni or
#'   Holm-Bonferroni correction.
#'
#' **Significance Shading:**
#' Regions where p-value < alpha are colored differently, highlighting q-ranges
#' where the two groups differ significantly.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap.
#' Chapman and Hall.
#'
#' Wood, S. N. (2017). Generalized additive models: an introduction with R (2nd ed.).
#' Chapman and Hall/CRC.
#'
#' @import ggplot2
#' @import tidyr
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'se' has multiple q values from calculate_diversity()
#' result <- plot_tsallis_q_curve_bootstrap(se, n_bootstrap = 500, ci_level = 0.95)
#'
#' # View the plot
#' result$ggplot
#'
#' # Check which q-ranges show significant differences
#' head(result$significant_qranges)
#'
#' # Extract the underlying data for custom plotting
#' plot_data <- result$plot_data
#' }
plot_tsallis_q_curve_bootstrap <- function(se, assay_name = "diversity",
                                           sample_type_col = "sample_type",
                                           n_bootstrap = 1000, ci_level = 0.95,
                                           ci_type = "pointwise", test_method = "wilcox",
                                           alpha = 0.05, method = "percentile") {

  require_pkgs(c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment"))

  # Validate inputs
  if (!inherits(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object")
  }

  if (!(assay_name %in% SummarizedExperiment::assayNames(se))) {
    stop("Assay '", assay_name, "' not found in SummarizedExperiment")
  }

  if (!(ci_type %in% c("pointwise", "simultaneous"))) {
    stop("'ci_type' must be 'pointwise' or 'simultaneous'")
  }

  if (!(test_method %in% c("wilcox", "ttest"))) {
    stop("'test_method' must be 'wilcox' or 'ttest'")
  }

  if (!(method %in% c("percentile", "bca", "normal"))) {
    stop("'method' must be 'percentile', 'bca', or 'normal'")
  }

  # Prepare data in long format
  long <- prepare_tsallis_long(se, assay_name = assay_name, sample_type_col = sample_type_col)
  if (nrow(long) == 0) {
    stop("No tsallis values found in SummarizedExperiment")
  }

  long$q <- as.numeric(as.character(long$q))
  unique_q <- sort(unique(long$q))
  n_q <- length(unique_q)

  if (n_q < 2) {
    stop("Need at least 2 q values for q-curve analysis")
  }

  # Extract group information
  if (!(sample_type_col %in% colnames(SummarizedExperiment::colData(se)))) {
    stop("'", sample_type_col, "' not found in colData")
  }

  groups <- unique(sort(long$group))
  if (length(groups) != 2) {
    stop("Expected exactly 2 groups, found ", length(groups))
  }

  cat("Computing bootstrap confidence intervals...\n")
  
  # Use helper function to compute bootstrap CIs
  bootstrap_results <- compute_bootstrap_qcurve_cis(
    long = long,
    unique_q = unique_q,
    groups = groups,
    ci_level = ci_level,
    n_bootstrap = n_bootstrap
  )

  cat("Testing for group differences at each q-value...\n")

  # Perform statistical tests at each q-value
  test_results <- list()
  for (q_val in unique_q) {
    q_data_g1 <- long %>%
      dplyr::filter(group == groups[1], q == q_val) %>%
      dplyr::pull(tsallis)

    q_data_g2 <- long %>%
      dplyr::filter(group == groups[2], q == q_val) %>%
      dplyr::pull(tsallis)

    if (length(q_data_g1) < 2 || length(q_data_g2) < 2) {
      test_results[[as.character(q_val)]] <- list(pvalue = NA_real_, significant = FALSE)
      next
    }

    if (test_method == "wilcox") {
      test_result <- wilcox.test(q_data_g1, q_data_g2, paired = FALSE)
      pvalue <- test_result$p.value
    } else {
      # t-test
      test_result <- t.test(q_data_g1, q_data_g2, var.equal = FALSE)
      pvalue <- test_result$p.value
    }

    test_results[[as.character(q_val)]] <- list(
      pvalue = pvalue,
      significant = pvalue < alpha
    )
  }

  # Build output data frame
  plot_df <- data.frame(
    q = numeric(),
    median = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    group = character(),
    pvalue = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )

  for (g in groups) {
    for (q_val in unique_q) {
      q_str <- as.character(q_val)
      bt_res <- bootstrap_results[[g]][[q_str]]
      test_res <- test_results[[q_str]]

      plot_df <- rbind(plot_df, data.frame(
        q = q_val,
        median = bt_res$median,
        ci_lower = bt_res$ci_lower,
        ci_upper = bt_res$ci_upper,
        group = g,
        pvalue = test_res$pvalue,
        significant = test_res$significant,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Identify significant q-ranges (consecutive q values with significant differences)
  significant_q <- plot_df %>%
    dplyr::filter(significant) %>%
    dplyr::pull(q) %>%
    unique() %>%
    sort()

  significant_qranges <- data.frame(
    q_min = numeric(),
    q_max = numeric(),
    n_tests = numeric(),
    min_pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  if (length(significant_q) > 0) {
    # Find contiguous ranges
    gaps <- which(diff(significant_q) > 0.01)  # Arbitrary threshold for gap detection
    range_starts <- c(1, gaps + 1)
    range_ends <- c(gaps, length(significant_q))

    for (i in seq_along(range_starts)) {
      q_range <- significant_q[range_starts[i]:range_ends[i]]
      min_pval <- min(plot_df$pvalue[plot_df$q %in% q_range], na.rm = TRUE)

      significant_qranges <- rbind(significant_qranges, data.frame(
        q_min = min(q_range),
        q_max = max(q_range),
        n_tests = length(q_range),
        min_pvalue = min_pval,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Create ggplot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = q, y = median, color = group, fill = group)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = "Tsallis q-curve with Bootstrap Confidence Bands",
      subtitle = paste0(
        "CI: ", round(ci_level * 100), "% (", n_bootstrap, " bootstrap replicates), ",
        "Test: ", test_method, ", α = ", alpha
      ),
      x = "q value",
      y = "Tsallis entropy (S_q)",
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, color = "gray50")
    )

  # Add significant region shading if there are significant q-ranges
  if (nrow(significant_qranges) > 0) {
    for (i in seq_len(nrow(significant_qranges))) {
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = significant_qranges$q_min[i],
          xmax = significant_qranges$q_max[i],
          ymin = -Inf, ymax = Inf,
          alpha = 0.1,
          fill = "red"
        )
    }
  }

  # Finalize plot styling
  p <- p +
    ggplot2::scale_color_manual(values = c("#1B9E77", "#D95F02")) +
    ggplot2::scale_fill_manual(values = c("#1B9E77", "#D95F02"))

  # If single group, hide legend
  if (length(groups) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # Prepare output
  result <- list(
    plot_data = plot_df,
    significant_qranges = significant_qranges,
    metadata = list(
      n_bootstrap = n_bootstrap,
      ci_level = ci_level,
      ci_type = ci_type,
      test_method = test_method,
      alpha = alpha,
      method = method,
      n_groups = length(groups),
      n_q = n_q
    ),
    ggplot = p
  )

  class(result) <- c("qcurve_bootstrap", "list")

  cat("✓ Bootstrap CI computation complete\n")
  cat("  ", nrow(significant_qranges), "significant q-range(s) identified\n")
  if (nrow(significant_qranges) > 0) {
    cat("  Q-ranges with significant group differences:\n")
    for (i in seq_len(nrow(significant_qranges))) {
      cat(
        "    [q = ", significant_qranges$q_min[i], " to ",
        significant_qranges$q_max[i], "]: ",
        "min p = ", formatC(significant_qranges$min_pvalue[i], format = "e", digits = 2),
        "\n"
      )
    }
  }

  return(result)
}

#' @method print qcurve_bootstrap
#' @export
print.qcurve_bootstrap <- function(x, ...) {
  cat("Tsallis Q-Curve Bootstrap Result:\n")
  cat("  Bootstrap replicates: ", x$metadata$n_bootstrap, "\n")
  cat("  CI level: ", x$metadata$ci_level * 100, "%\n")
  cat("  Test method: ", x$metadata$test_method, "\n")
  cat("  Significance level (α): ", x$metadata$alpha, "\n")
  cat("  Q-values tested: ", x$metadata$n_q, "\n")
  cat("  Groups: ", x$metadata$n_groups, "\n")
  cat("  Significant q-ranges: ", nrow(x$significant_qranges), "\n")

  if (nrow(x$significant_qranges) > 0) {
    cat("\n  Ranges where groups differ significantly:\n")
    print(x$significant_qranges)
  } else {
    cat("\n  No significant differences detected at the current α level.\n")
  }
}

# Helpers for plot_top_transcripts internals

.ptt_select_genes_from_res <- function(res, top_n) {
    if (is.null(res)) {
        stop("Either 'gene' or 'res' must be provided")
    }
    if (!("genes" %in% colnames(res))) {
        stop("Provided 'res' must contain a 'genes' column")
    }
    # Look for adjusted p-value columns in order of preference
    if ("padj" %in% colnames(res)) {
        ord <- order(res$padj, na.last = NA)
    } else if ("adjusted_p_values" %in% colnames(res)) {
        ord <- order(res$adjusted_p_values, na.last = NA)
    } else if ("pvalue" %in% colnames(res)) {
        ord <- order(res$pvalue, na.last = NA)
    } else if ("raw_p_values" %in% colnames(res)) {
        ord <- order(res$raw_p_values, na.last = NA)
    } else {
        ord <- seq_len(nrow(res))
    }
    genes_sel <- as.character(res$genes[ord])
    genes_sel <- unique(genes_sel)
    head(genes_sel, top_n)
}

.ptt_infer_samples_from_coldata <- function(coldata, counts, sample_type_col) {
    if (is.character(coldata) && length(coldata) == 1) {
        if (!file.exists(coldata)) {
            stop("coldata file not found: ", coldata)
        }
        cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(coldata)) {
        cdf <- coldata
    } else {
        stop("`coldata` must be a data.frame or path to a tab-delimited file")
    }

    if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
        as.character(cdf[colnames(counts), sample_type_col])
    } else {
        sample_id_cols <- c("sample", "Sample", "sample_id", "id")
        sid <- intersect(sample_id_cols, colnames(cdf))
        if (length(sid) > 0) {
            sid <- sid[1]
            if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
                stop("coldata sample id column does not match column names of counts")
            }
            row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
            as.character(cdf[[sample_type_col]][row_ix])
        } else {
            stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
        }
    }
}

.ptt_read_tx2gene <- function(tx2gene) {
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene)
        }
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }
    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
        stop("tx2gene must have columns 'Transcript' and 'Gen'")
    }
    mapping
}

.ptt_make_agg <- function(metric = c("median", "mean", "variance", "iqr")) {
    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) {
            stats::var(x, na.rm = TRUE)
        }, iqr = function(x) stats::IQR(x, na.rm = TRUE))
    agg_label_metric <- if (metric_choice == "iqr") {
        "IQR"
    } else {
        metric_choice
    }
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label
    list(metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique)
}

.ptt_build_tx_long <- function(gene_single, mapping, counts, samples, top_n) {
    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))
    if (length(txs) == 0) {
        stop("No transcripts found for gene: ", gene_single)
    }
    if (!is.null(top_n)) {
        txs <- head(txs, top_n)
    }
    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)
    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))
    list(df_long = df_long, txs = txs)
}

.ptt_aggregate_df_long <- function(df_long, agg_fun, pseudocount) {
    df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
    df_summary
}

.ptt_build_plot_from_summary <- function(df_summary, agg_label_unique, fill_limits = NULL) {
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
        ggplot2::geom_tile(color = "white", width = 0.95, height = 0.95) + ggplot2::scale_fill_viridis_c(option = "viridis",
        direction = -1, na.value = "grey80", limits = fill_limits) + ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(title = agg_label_unique, x = NULL, y = NULL, fill = "log2(expr)") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12), axis.text.x = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(size = 16, hjust = 0.6), legend.position = "bottom",
            legend.key.width = ggplot2::unit(1.2, "cm"), plot.margin = ggplot2::margin(4,
                4, 4, 4)) + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
        barwidth = 6, barheight = 0.35))
    p
}

.ptt_combine_patchwork <- function(plots, agg_label_unique) {
    # Use 3 columns, let patchwork auto-calculate rows
    combined <- Reduce(`+`, plots) + patchwork::plot_layout(ncol = 3, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    combined <- combined + patchwork::plot_annotation(title = agg_label_unique, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.6,
        size = 16, margin = ggplot2::margin(b = 10))))
    combined
}

.ptt_combine_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
    p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
    
    # Use 3 columns, auto-calculate rows
    ncol <- 3
    nrow_val <- ceiling(length(plots_nolegend) / ncol)
    
    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val, align = "hv")
    title_grob <- cowplot::ggdraw() + cowplot::draw_label(agg_label_unique, fontface = "plain",
        x = 0.6, hjust = 0.5, size = 16)
    result_plot <- cowplot::plot_grid(title_grob, grid, legend, ncol = 1, rel_heights = c(0.08,
        1, 0.08))
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    }
    result_plot
}

.ptt_combine_grid <- function(plots, output_file = NULL, agg_label_unique) {
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
    grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)
    g_full <- ggplot2::ggplotGrob(plots[[1]])
    legend_idx <- which(vapply(g_full$grobs, function(x) x$name, character(1)) ==
        "guide-box")
    if (length(legend_idx)) {
        legend_grob <- g_full$grobs[[legend_idx[1]]]
    } else {
        legend_grob <- NULL
    }
    
    # Default to 3 columns, adjust for smaller numbers
    ncol <- min(3, length(grobs))
    nrow <- ceiling(length(grobs) / ncol)
    
    # Create heights: title (0.6cm) + plot rows (1 null each) + legend (0.7cm)
    plot_heights <- rep(grid::unit(1, "null"), nrow)
    heights <- grid::unit.c(grid::unit(0.6, "cm"), plot_heights, grid::unit(0.7, "cm"))
    
    if (!is.null(output_file)) {
        # Adjust PNG dimensions based on layout
        png_width <- 800 * ncol
        png_height <- 480 * nrow
        png(filename = output_file, width = png_width, height = png_height, res = 150)
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights, to_file = output_file)
        invisible(NULL)
    } else {
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights)
        invisible(NULL)
    }
}

#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_line scale_y_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text geom_hline geom_vline
#' @importFrom ggplot2 scale_color_manual scale_shape_manual geom_tile scale_fill_gradient2
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "dimension",
      "inertia",
      "cumulative_inertia",
      "variable",
      "contribution",
      "dim1",
      "dim2",
      "type",
      "coord_x",
      "coord_y"
    )
  )
}


#' Plot Inertia (Variance Explained) by Dimension
#'
#' Visualizes the cumulative proportion of variance explained by successive
#' dimensions in MCA. This helps determine which dimensions are important
#' for understanding patterns in correspondence analysis (Abdi & Valentin 2007, 
#' Khangar & Kamalja 2017).
#'
#' @param entropy_matrix Matrix of entropy values (genes * q-values)
#' @param q_values Numeric vector of q values
#' @param n_dims Integer; number of dimensions to plot (default: 5)
#' @param title Character; plot title
#'
#' @return A `ggplot2` object showing inertia per dimension with cumulative line
#'
#' @details
#' The inertia (measure of variance in CA) represents how much of the total
#' association between genes and q-values is captured by each principal
#' dimension. Typically, 2-3 dimensions capture 70-85% of inertia.
#'
#' @export
plot_ca_inertia <- function(entropy_matrix,
                             q_values,
                             n_dims = 5,
                             title = "Correspondence Analysis: Explained Inertia by Dimension") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 required for plotting", call. = FALSE)
  }
  
  # Require FactoMineR for proper CA computation
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    warning("FactoMineR not available; returning simplified inertia plot", call. = FALSE)
    
    # Simplified fallback: compute variance from entropy directly
    inertia_vals <- apply(entropy_matrix, 2, var)
    inertia_vals <- inertia_vals / sum(inertia_vals)
    inertia_vals <- utils::head(inertia_vals, n_dims)
  } else {
    # Categorize entropy values for MCA
    entropy_cat <- apply(entropy_matrix, 2, function(x) {
      cut(x, breaks = 3, labels = c("low", "med", "high"), include.lowest = TRUE)
    })
    
    # Set column names - ensure they match entropy_cat dimensions
    q_labels <- if (length(q_values) == ncol(entropy_cat)) {
      paste0("q_", round(q_values, 2))
    } else {
      paste0("q_", seq_len(ncol(entropy_cat)))
    }
    colnames(entropy_cat) <- q_labels
    
    # Run MCA
    mca_res <- tryCatch(
      FactoMineR::MCA(entropy_cat, graph = FALSE, ncp = min(n_dims, nrow(entropy_cat) - 1)),
      error = function(e) {
        warning("MCA computation failed; using uniform inertia", call. = FALSE)
        NULL
      }
    )
    
    if (is.null(mca_res)) {
      inertia_vals <- rep(1/n_dims, n_dims)
    } else {
      inertia_vals <- mca_res$eig[, 1]  # Eigenvalues (inertia per dimension)
      # Ensure we only have n_dims values
      if (length(inertia_vals) > n_dims) {
        inertia_vals <- inertia_vals[1:n_dims]
      }
    }
  }
  
  # Ensure inertia_vals has correct length and non-zero
  inertia_vals <- as.numeric(inertia_vals)
  if (anyNA(inertia_vals) || all(inertia_vals == 0)) {
    inertia_vals <- rep(1/n_dims, n_dims)
  }
  
  if (length(inertia_vals) < n_dims) {
    inertia_vals <- c(inertia_vals, rep(0, n_dims - length(inertia_vals)))
  } else if (length(inertia_vals) > n_dims) {
    inertia_vals <- inertia_vals[1:n_dims]
  }
  
  # Compute cumulative inertia
  cum_inertia <- cumsum(inertia_vals) / sum(inertia_vals)
  
  # Create plot data
  plot_data <- data.frame(
    dimension = seq_along(inertia_vals),
    inertia = inertia_vals / sum(inertia_vals),
    cumulative_inertia = cum_inertia
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dimension)) +
    ggplot2::geom_col(ggplot2::aes(y = inertia),
                      fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_inertia),
                        color = "darkred", size = 3) +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_inertia),
                       color = "darkred", linewidth = 1) +
    ggplot2::scale_y_continuous(
      name = "Proportion of Inertia",
      limits = c(0, 1),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(
      name = "Dimension",
      breaks = seq_along(inertia_vals)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Bars: individual inertia | Line: cumulative inertia"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray60")
    )
  
  return(p)
}


#' Plot Contribution of Variables to Principal Dimensions
#'
#' Creates a heatmap showing how much each variable (q-value) contributes
#' to the first two principal dimensions. This reveals which q-values
#' drive the main patterns (Khangar & Kamalja 2017).
#'
#' @param entropy_matrix Matrix of entropy values (genes * q-values)
#' @param q_values Numeric vector of q values
#' @param n_variables_labeled Integer; number of variables (q-values) to label in heatmap (default: 5)
#' @param title Character; plot title
#'
#' @return A `ggplot2` object showing variable contributions as a heatmap
#'   with q-values as columns (x-axis) and dimensions as rows (y-axis)
#'
#' @details
#' The heatmap displays how much each q-value contributes to the first two
#' CA dimensions. Q-values with high contributions (red cells) are the main
#' drivers of the correspondence analysis dimensions.
#'
#' @export
plot_ca_contributions <- function(entropy_matrix,
                                  q_values,
                                  n_variables_labeled = 5,
                                  title = "Correspondence Analysis: Variable Contributions") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 required for plotting", call. = FALSE)
  }
  
  # Categorize entropy for MCA
  entropy_cat <- apply(entropy_matrix, 2, function(x) {
    cut(x, breaks = 3, labels = c("low", "med", "high"), include.lowest = TRUE)
  })
  colnames(entropy_cat) <- paste0("q_", round(q_values, 2))
  
  # Run MCA if available
  if (requireNamespace("FactoMineR", quietly = TRUE)) {
    mca_res <- tryCatch(
      FactoMineR::MCA(entropy_cat, graph = FALSE, ncp = 2),
      error = function(e) NULL
    )
    
    if (!is.null(mca_res)) {
      # Extract variable contributions
      contrib <- mca_res$var$contrib
      
      # Prepare data for heatmap
      contrib_data <- data.frame(
        variable = rownames(contrib),
        dim1 = contrib[, 1],
        dim2 = contrib[, 2]
      )
    } else {
      contrib_data <- NULL
    }
  } else {
    contrib_data <- NULL
  }
  
  # If MCA failed, compute simple contribution from entropy variance
  if (is.null(contrib_data)) {
    var_explained <- apply(entropy_matrix, 2, var)
    var_explained <- 100 * var_explained / sum(var_explained)
    
    contrib_data <- data.frame(
      variable = paste0("q_", round(q_values, 2)),
      dim1 = var_explained,
      dim2 = abs(scale(entropy_matrix[, 1])[, 1])
    )
  }
  
  # Reshape for heatmap - convert to long format without reshape2
  contrib_long <- data.frame(
    variable = c(contrib_data$variable, contrib_data$variable),
    dimension = c(rep("Dim1", nrow(contrib_data)), rep("Dim2", nrow(contrib_data))),
    contribution = c(contrib_data$dim1, contrib_data$dim2)
  )
  
  # Create heatmap
  p <- ggplot2::ggplot(contrib_long, ggplot2::aes(x = variable, y = dimension, fill = contribution)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "white", mid = "lightyellow", high = "darkred",
      name = "Contribution (%)"
    ) +
    ggplot2::labs(
      title = title,
      x = "Q-values",
      y = "Dimension"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}


#' Plot Row-Column Biplot for Multiple Correspondence Analysis
#'
#' Creates an MCA biplot showing both genes (rows) and q-value categories (columns)
#' in the same principal coordinate space. This allows visualization of how genes
#' are associated with entropy categories (high/medium/low) at different q-values.
#'
#' **Important**: This function performs MCA on *categorized* entropy data (tertiles),
#' not classical Correspondence Analysis on continuous values. Entropy values are
#' automatically categorized into 3 levels (low, medium, high) before analysis.
#'
#' @param entropy_matrix Matrix of entropy values (genes * q-values)
#' @param q_values Numeric vector of q values
#' @param n_genes_labeled Integer; number of extreme genes to label (default: 8)
#' @param title Character; plot title
#' @param show_origin Logical; if TRUE, add origin lines (default: TRUE)
#'
#' @return ggplot2 object showing row-column biplot
#'
#' @details
#' In this MCA biplot:
#' - **Rows (genes)**: Positioned by entropy category associations at each q-value
#' - **Columns (q-value categories)**: Named as "q_X_high", "q_X_low", "q_X_med" etc.
#' - **Proximity**: Genes close together show similar entropy patterns across q-values
#' - **Distance from origin**: Genes far from origin have distinctive entropy signatures
#'
#' Interpretation (Le Roux & Rouanet 2011, modified for categorical),
#' - Genes clustered together: similar entropy category profiles across q-values
#' - Q-value categories clustered: drive entropy similarly across genes
#' - Separate clusters: represent distinct biological patterns
#' - Note: This reflects entropy *categories*, not continuous values
#'
#' @export
plot_ca_biplot <- function(entropy_matrix,
                           q_values,
                           n_genes_labeled = 8,
                           title = "Correspondence Analysis: Gene-Q Biplot",
                           show_origin = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 required for plotting", call. = FALSE)
  }
  
  # Categorize entropy for MCA
  entropy_cat <- apply(entropy_matrix, 2, function(x) {
    cut(x, breaks = 3, labels = c("low", "med", "high"), include.lowest = TRUE)
  })
  colnames(entropy_cat) <- paste0("q_", round(q_values, 2))
  rownames(entropy_cat) <- rownames(entropy_matrix)
  
  # Run MCA
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    warning("FactoMineR required for proper biplot; using PCA-based approximation", call. = FALSE)
    
    # Fallback: use PCA on standardized entropy
    entropy_std <- scale(entropy_matrix)
    pca_res <- stats::prcomp(entropy_std, scale. = TRUE)
    row_coords <- pca_res$x[, 1:2]
    col_coords <- pca_res$rotation[, 1:2] * 2  # Scale for visibility
    
    row_df <- data.frame(
      label = rownames(entropy_matrix),
      dim1 = row_coords[, 1],
      dim2 = row_coords[, 2],
      type = "Gene"
    )
    col_df <- data.frame(
      label = colnames(entropy_cat),
      dim1 = col_coords[, 1],
      dim2 = col_coords[, 2],
      type = "Q-value"
    )
  } else {
    mca_res <- tryCatch(
      FactoMineR::MCA(entropy_cat, graph = FALSE, ncp = 2),
      error = function(e) {
        warning("MCA computation failed; using PCA approximation instead", call. = FALSE)
        NULL
      }
    )
    
    if (is.null(mca_res)) {
      # Fallback: use PCA on standardized entropy
      entropy_std <- scale(entropy_matrix)
      pca_res <- stats::prcomp(entropy_std, scale. = TRUE)
      row_coords <- pca_res$x[, 1:2]
      col_coords <- pca_res$rotation[, 1:2] * 2  # Scale for visibility
      
      row_df <- data.frame(
        label = rownames(entropy_matrix),
        dim1 = row_coords[, 1],
        dim2 = row_coords[, 2],
        type = "Gene"
      )
      col_df <- data.frame(
        label = colnames(entropy_cat),
        dim1 = col_coords[, 1],
        dim2 = col_coords[, 2],
        type = "Q-value"
      )
    } else {
      # Extract row and column coordinates
      row_coords <- mca_res$ind$coord
      col_coords <- mca_res$var$coord
      
      # Ensure we have 2 dimensions
      if (ncol(row_coords) < 2) {
        row_coords <- cbind(row_coords, rep(0, nrow(row_coords)))
      } else {
        row_coords <- row_coords[, 1:2]
      }
      
      if (ncol(col_coords) < 2) {
        col_coords <- cbind(col_coords, rep(0, nrow(col_coords)))
      } else {
        col_coords <- col_coords[, 1:2]
      }
      
      col_coords <- col_coords * 2  # Scale for visibility
      
      # Ensure row and column labels match matrix dimensions
      row_labs <- if (!is.null(rownames(mca_res$ind$coord))) {
        rownames(mca_res$ind$coord)
      } else {
        rownames(entropy_matrix)
      }
      
      col_labs <- if (!is.null(rownames(mca_res$var$coord))) {
        rownames(mca_res$var$coord)
      } else {
        colnames(entropy_cat)
      }
      
      row_df <- data.frame(
        label = row_labs,
        dim1 = row_coords[, 1],
        dim2 = row_coords[, 2],
        type = "Gene",
        stringsAsFactors = FALSE
      )
      
      col_df <- data.frame(
        label = col_labs,
        dim1 = col_coords[, 1],
        dim2 = col_coords[, 2],
        type = "Q-value",
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine and identify genes to label (edges of cloud)
  all_points <- rbind(row_df, col_df)
  
  # Label genes at extremes
  genes_to_label <- row_df[
    order(sqrt(row_df$dim1^2 + row_df$dim2^2), decreasing = TRUE)[1:min(n_genes_labeled, nrow(row_df))],
  ]
  
  # Create biplot
  p <- ggplot2::ggplot(all_points, ggplot2::aes(x = dim1, y = dim2, color = type, shape = type)) +
    ggplot2::geom_point(size = 3, alpha = 0.6) +
    ggplot2::geom_text(
      data = genes_to_label,
      ggplot2::aes(label = label),
      vjust = -1.2, size = 3, color = "black", fontface = "italic"
    ) +
    ggplot2::scale_color_manual(
      values = c("Gene" = "steelblue", "Q-value" = "darkred"),
      name = "Type"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Gene" = 16, "Q-value" = 17),
      name = "Type"
    )
  
  # Add origin lines if requested
  if (show_origin) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5)
  }
  
  p <- p +
    ggplot2::labs(
      title = title,
      x = paste0("Dimension 1 (", round(100 * var(all_points$dim1) / (var(all_points$dim1) + var(all_points$dim2)), 1), "%)"),
      y = paste0("Dimension 2 (", round(100 * var(all_points$dim2) / (var(all_points$dim1) + var(all_points$dim2)), 1), "%)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}


#' Comprehensive Correspondence Analysis Visualization
#'
#' Creates a multi-panel visualization combining inertia, contributions, and biplot.
#' Provides complete picture of CA results as recommended by leading references.
#'
#' @param entropy_matrix Matrix of entropy values (genes * q-values)
#' @param q_values Numeric vector of q values
#' @param include_biplot Logical; if TRUE, include the row-column biplot (default: TRUE)
#' @param include_contributions Logical; if TRUE, include contribution heatmap (default: TRUE)
#' @param main_title Character; overall title for the figure
#'
#' @return List of `ggplot2` objects (or combined patchwork if patchwork is available)
#'
#' @export
plot_ca_comprehensive <- function(entropy_matrix,
                                  q_values,
                                  include_biplot = TRUE,
                                  include_contributions = TRUE,
                                  main_title = "Comprehensive Correspondence Analysis") {
  
  plots <- list()
  
  # Always include inertia plot
  plots$inertia <- plot_ca_inertia(entropy_matrix, q_values)
  
  # Add other plots as requested
  if (include_contributions) {
    plots$contributions <- plot_ca_contributions(entropy_matrix, q_values)
  }
  
  if (include_biplot) {
    plots$biplot <- plot_ca_biplot(entropy_matrix, q_values)
  }
  
  # Try to combine with patchwork if available
  if (requireNamespace("patchwork", quietly = TRUE) && length(plots) > 1) {
    combined <- Reduce(function(x, y) x + y, plots) +
      patchwork::plot_annotation(
        title = main_title,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))
      )
    return(combined)
  } else {
    return(plots)
  }
}

# Internal plot helpers

.tsenat_format_label <- function(lbl) {
    if (is.null(lbl)) {
        return(NULL)
    }
    s <- gsub("_", " ", lbl)
    s <- gsub("\\s+", " ", s)
    s <- trimws(s)
    s <- tolower(s)
    if (nchar(s) == 0) {
        return(s)
    }
    if (nchar(s) == 1) {
        return(toupper(s))
    }
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

.tsenat_prepare_ma_plot_df <- function(df, fold_col, mean_cols, x_label, y_label) {
    # Detect x-axis values
    if (length(mean_cols) >= 2) {
        xvals <- rowMeans(df[, mean_cols[seq_len(2)], drop = FALSE], na.rm = TRUE)
        x_label <- x_label %||% paste0(mean_cols[1], " vs ", mean_cols[2])
    } else if (length(mean_cols) == 1) {
        xvals <- as.numeric(df[[mean_cols[1]]])
        x_label <- x_label %||% mean_cols[1]
    } else if ("mean" %in% colnames(df)) {
        xvals <- as.numeric(df$mean)
        x_label <- x_label %||% "Mean"
    } else {
        xvals <- seq_len(nrow(df))
        x_label <- x_label %||% "Index"
    }

    yvals <- as.numeric(df[[fold_col]])

    padj_candidates <- c("adjusted_p_values", "adj_p_value", "adj_p", "padj", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) {
        padj_col[1]
    } else {
        NULL
    }

    padj <- if (!is.null(padj_col)) {
        as.numeric(df[[padj_col]])
    } else {
        rep(1, length(yvals))
    }
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < 0.05, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag,
        stringsAsFactors = FALSE)

    list(plot_df = plot_df, x_label = x_label, y_label = y_label)
}

.tsenat_prepare_volcano_df <- function(diff_df, x_col = NULL, padj_col = "adjusted_p_values",
    label_thresh = 0.1, sig_alpha = 0.05, title = NULL) {
    df <- as.data.frame(diff_df)
    cn <- colnames(df)

    # Auto-detect x-axis column if not specified
    if (is.null(x_col)) {
        diff_cols <- grep("_difference$", cn, value = TRUE, ignore.case = TRUE)
        if (length(diff_cols) > 0) {
            x_col <- diff_cols[1]
        } else {
            numeric_cols <- vapply(df, is.numeric, logical(1))
            p_cols <- grep("p_value|p.value", cn, ignore.case = TRUE)
            numeric_cols[p_cols] <- FALSE
            if (any(numeric_cols)) {
                x_col <- cn[which(numeric_cols)[1]]
            } else {
                stop("Could not find suitable column for x-axis. Specify 'x_col' explicitly.")
            }
        }
    }

    # Verify columns
    if (!(x_col %in% cn)) {
        stop(sprintf("Column '%s' not found in diff_df", x_col))
    }
    if (!(padj_col %in% cn)) {
        stop(sprintf("Column '%s' not found in diff_df", padj_col))
    }

    df$xval <- as.numeric(df[[x_col]])
    df$padj <- as.numeric(df[[padj_col]])
    df$padj[is.na(df$padj)] <- 1
    df$padj[df$padj <= 0] <- .Machine$double.xmin

    df$significant <- ifelse(abs(df$xval) >= label_thresh & df$padj < sig_alpha,
        "significant", "non-significant")
    df <- df[is.finite(df$xval) & is.finite(df$padj), ]

    if (nrow(df) == 0) {
        stop("No valid points to plot")
    }

    metric_label <- if (grepl("median", x_col, ignore.case = TRUE)) {
        "Median"
    } else if (grepl("mean", x_col, ignore.case = TRUE)) {
        "Mean"
    } else {
        "Value"
    }

    title_use <- title %||% "Volcano plot: fold-change vs significance"

    x_label_formatted <- .tsenat_format_label(x_col)
    padj_label_formatted <- .tsenat_format_label(padj_col)

    list(df = df, x_col = x_col, padj_col = padj_col, x_label_formatted = x_label_formatted,
        padj_label_formatted = padj_label_formatted, title_use = title_use)
}


#' Plot Tsallis Divergence Effect Size Distribution
#'
#' Generate a histogram visualization of Tsallis divergence effect sizes across genes,
#' showing the distribution of information-theoretic measures of isoform switching.
#'
#' @param interaction_results A data frame containing LMM results merged with per-q divergence estimates.
#'   Must contain columns matching the pattern `effect_size_D_q*` (e.g., `effect_size_D_q0_5`, `effect_size_D_q1_0`).
#'   Typically the result from [effect_sizes_divergence()].
#'
#' @param threshold Numeric. Effect size threshold for visual marking. Default is 0.1 (information-theoretic significance level).
#'
#' @return If ggplot2 is available and `interaction_results` contains valid data, returns a ggplot object.
#'   Otherwise returns NULL invisibly and prints an informative message.
#'
#' @details
#' The function visualizes the distribution of effect sizes using the median q-value's divergence
#' (typically around q=1.0, close to Shannon entropy). The red dashed line marks the default
#' information-theoretic significance threshold of D=0.1.
#'
#' @references
#' - Chanda et al. (2020). Information Theory in Computational Biology. *Entropy*, 22(6), 627.
#' - Tsallis, C. (1988). Possible Generalization of Boltzmann-Gibbs Statistics. *Journal of Statistical Physics*, 52(1), 479-487.
#'
#' @examples
#' \dontrun{
#'   # Assuming lmm_results from effect_sizes_divergence()
#'   plot_divergence_distribution(lmm_results$interaction_results, threshold = 0.1)
#' }
#'
#' @export
plot_divergence_distribution <- function(interaction_results, threshold = 0.1) {
  
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("ggplot2 package required for plotting. Please install: install.packages('ggplot2')\n")
    return(invisible(NULL))
  }
  
  # Validate input
  if (is.null(interaction_results) || nrow(interaction_results) == 0) {
    cat("Plot not generated: interaction_results is empty or NULL.\n")
    return(invisible(NULL))
  }
  
  # Find per-q effect size columns
  effect_cols <- grep("^effect_size_D_q", colnames(interaction_results), value = TRUE)
  
  if (length(effect_cols) == 0) {
    cat("Plot not generated: no per-q effect size columns found in interaction_results.\n")
    cat("Expected columns like 'effect_size_D_q0_5', 'effect_size_D_q1_0', etc.\n")
    return(invisible(NULL))
  }
  
  # Use median q effect size for visualization
  median_idx <- ceiling(length(effect_cols) / 2)
  median_col <- effect_cols[median_idx]
  
  # Create visualization of effect size distribution
  p_effect <- ggplot2::ggplot(interaction_results, ggplot2::aes(x = .data[[median_col]])) +
    ggplot2::geom_histogram(binwidth = 0.02, fill = "steelblue", alpha = 0.7, color = "black") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(
      title = "Distribution of Tsallis Divergence (D_q) effect sizes across genes",
      subtitle = "Information-theoretic measure respecting Tsallis multi-q entropy properties",
      x = paste("Effect size (Tsallis Divergence D_q; D >", threshold, "= meaningful information separation)"),
      y = "Number of genes",
      caption = paste("Red dashed line: D =", threshold, "filtering threshold (information-theoretic significance for q-dependent entropy)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40"),
      panel.grid.major = ggplot2::element_line(color = "gray90")
    ) +
    ggplot2::annotate("text", x = threshold, y = Inf, 
                      label = paste("Information\nthreshold\n(D=", threshold, ")", sep = ""),
                      vjust = 1.5, hjust = -0.1, color = "red", size = 3.5)
  
  return(p_effect)
}


#' Plot Multi-Gene Q-Spectrum Comparison
#'
#' Generate a multi-panel plot comparing q-spectra (per-q divergence profiles) 
#' across the top N significant genes by effect size.
#'
#' @param lmm_results A list returned by [effect_sizes_divergence()], containing:
#'   - `$interaction_results`: Data frame with per-q divergence columns and `per_q_pattern` column
#'   Other list elements are ignored.
#'
#' @param n_genes Numeric. Number of top genes to display. Default is 5.
#'
#' @return Invisibly returns NULL. Side effect: plots the multi-panel q-spectrum comparison to the current graphics device.
#'
#' @details
#' Each panel shows one gene's q-spectrum from q=0.5 (rare isoforms) to q=2.0 (abundant isoforms),
#' with a vertical reference line at q=1.0 (Shannon entropy / Kullback-Leibler divergence).
#'
#' Genes are sorted by the median q effect size (typically near q=1.0).
#'
#' @export
plot_multi_q_spectrum <- function(lmm_results, n_genes = 5) {
  
  # Check for required data structure
  if (!is.list(lmm_results) || is.null(lmm_results$interaction_results)) {
    cat("plot_multi_q_spectrum requires a list with $interaction_results component.\n")
    cat("Typically the output from effect_sizes_divergence().\n")
    return(invisible(NULL))
  }
  
  interaction_results <- lmm_results$interaction_results
  
  # Additional safety check
  if (!is.data.frame(interaction_results) || nrow(interaction_results) == 0) {
    cat("No valid genes in interaction_results.\n")
    return(invisible(NULL))
  }
  
  # Get top genes by effect size (using median q effect size)
  effect_cols <- grep("^effect_size_D_q", colnames(interaction_results), value = TRUE)
  
  # Check if effect columns were found
  if (length(effect_cols) == 0) {
    cat("No effect_size columns found in LMM results.\n")
    return(invisible(NULL))
  }
  
  median_idx <- ceiling(length(effect_cols) / 2)
  median_col <- effect_cols[median_idx]
  
  # Sort and get top genes
  sort_order <- order(interaction_results[[median_col]], decreasing = TRUE, na.last = TRUE)
  top_genes <- head(interaction_results[sort_order, , drop = FALSE], n_genes)
  
  # Set up multi-panel plot
  n_panels <- nrow(top_genes)
  
  # Check if there are genes to plot
  if (n_panels == 0) {
    cat("No genes with valid data found.\n")
    return(invisible(NULL))
  }
  
  # Check that we have at least 1 panel
  if (n_panels < 1) {
    cat("Error: n_panels is less than 1.\n")
    return(invisible(NULL))
  }
  
  old_par <- par(mfrow = c(1, n_panels), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)  # Ensure par is reset even on error
  
  # Plot each gene
  for (i in 1:n_panels) {
    gene_data <- top_genes[i, ]
    gene_name <- gene_data$gene
    per_q_str <- gene_data$per_q_pattern
    
    if (!is.na(per_q_str) && is.character(per_q_str) && nchar(per_q_str) > 0) {
      tryCatch({
        per_q_vals <- as.numeric(strsplit(per_q_str, ",")[[1]])
        
        if (length(per_q_vals) > 0 && all(is.finite(per_q_vals))) {
          # Create simplified plot for multi-panel
          q_vals <- c(0.5, 1.0, 1.5, 2.0)
          
          # Ensure q_vals and per_q_vals have compatible lengths
          if (length(q_vals) == length(per_q_vals)) {
            y_max <- max(per_q_vals, na.rm = TRUE)
            y_range <- c(0, max(y_max * 1.1, 0.1))  # Avoid zero range
            
            plot(q_vals, per_q_vals,
                 main = gene_name,
                 xlab = "q",
                 ylab = "D_q",
                 type = "b", pch = 19, lwd = 2,
                 ylim = y_range)
            abline(v = 1, lty = 3, col = "gray")
          }
        }
      }, error = function(e) {
        # Silently skip this gene on error
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0.5, 0.5, "Plot error", cex = 0.8, adj = c(0.5, 0.5))
      })
    }
  }
  
  return(invisible(NULL))
}


#' Plot Q-Spectrum Curve with Confidence Intervals
#'
#' Visualize divergence across Tsallis sensitivity parameters (q-spectrum)
#' for a single gene, showing the complete biological signature of isoform switching.
#'
#' @param per_q_div Numeric vector of per-q divergence values OR an lmm_results list object.
#'   If a list (lmm_results), must contain `$interaction_results` with genes and per_q_pattern column.
#'   If numeric, names should be in format "q_0.5", "q_1.0", etc. If unnamed, assumes sequential
#'   q values from 0.5 to 2.0 in 0.5 increments.
#'
#' @param gene_idx Integer. If `per_q_div` is lmm_results, which gene to plot (default: 1 = first gene).
#'   Ignored if per_q_div is a numeric vector.
#'
#' @param per_q_ci Optional list with components `$lower` and `$upper` containing
#'   lower and upper bootstrap confidence interval bounds (same length as per_q_div).
#'   Only used if per_q_div is a numeric vector.
#'
#' @param gene_name Character. Name of the gene for plot title. Default is empty string or auto-detected.
#'
#' @return A ggplot2 object displaying the q-spectrum curve. Can be displayed with print() or
#'   combined with other ggplot2 operations.
#'
#' @details
#' The q-spectrum encodes which abundance scales are affected by isoform switching:
#' 
#' - **q < 1 (left side)**: Emphasizes rare isoforms. High divergence indicates rare variants differ between groups.
#' - **q = 1 (middle)**: Kullback-Leibler divergence; the "average" effect with equal weighting.
#' - **q > 1 (right side)**: Emphasizes abundant isoforms. High divergence indicates major isoforms rebalance.
#'
#' The curve shape reveals the biological mechanism:
#' - **RARE_DRIVEN (declining)**: D(q=0.5) >> D(q=2.0), rare isoforms dominate
#' - **BALANCED (flat)**: Divergence similar across all q, all scales affected equally
#' - **ABUNDANT_DRIVEN (rising)**: D(q=2.0) >> D(q=0.5), dominant isoforms rebalance
#'
#' @examples
#' \dontrun{
#'   # Example 1: Direct per-q vector
#'   per_q_vec <- c(0.15, 0.10, 0.08, 0.05)
#'   names(per_q_vec) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
#'   plot_q_spectrum(per_q_vec, gene_name = "Example Gene")
#'   
#'   # Example 2: From lmm_results
#'   plot_q_spectrum(lmm_results, gene_idx = 1)
#' }
#'
#' @export
plot_q_spectrum <- function(per_q_div, gene_idx = 1, per_q_ci = NULL, gene_name = "") {
  
  # Check if per_q_div is actually lmm_results (list with interaction_results)
  if (is.list(per_q_div) && !is.null(per_q_div$interaction_results)) {
    lmm_results <- per_q_div
    
    # Extract from lmm_results
    interaction_results <- lmm_results$interaction_results
    
    if (is.null(interaction_results) || nrow(interaction_results) == 0) {
      cat("Error: lmm_results$interaction_results is empty.\n")
      return(invisible(NULL))
    }
    
    if (!"per_q_pattern" %in% colnames(interaction_results)) {
      cat("Error: 'per_q_pattern' column not found in lmm_results$interaction_results.\n")
      return(invisible(NULL))
    }
    
    # Get the requested gene
    if (gene_idx < 1 || gene_idx > nrow(interaction_results)) {
      cat("Error: gene_idx out of range. Available genes:", nrow(interaction_results), "\n")
      return(invisible(NULL))
    }
    
    gene_data <- interaction_results[gene_idx, ]
    gene_name <- if (is.na(gene_name) || gene_name == "") gene_data$gene else gene_name
    pattern_str <- gene_data$per_q_pattern
    
    if (is.na(pattern_str) || nchar(pattern_str) == 0) {
      cat("Error: No valid per_q_pattern for gene", gene_name, "\n")
      return(invisible(NULL))
    }
    
    # Parse comma-separated per-q values
    per_q_vals <- as.numeric(strsplit(pattern_str, ",")[[1]])
    
    if (length(per_q_vals) == 0 || !all(is.finite(per_q_vals))) {
      cat("Error: Invalid per_q_pattern values for gene", gene_name, "\n")
      return(invisible(NULL))
    }
    
    per_q_div <- setNames(per_q_vals, paste0("q_", seq(0.5, by=0.5, length.out=length(per_q_vals))))
  }
  
  # Handle empty or missing input
  if (is.null(per_q_div) || length(per_q_div) == 0) {
    cat("Error: per_q_div is empty or NULL.\n")
    return(invisible(NULL))
  }
  
  # Extract q values from names (e.g., "q_0.5", "q_1.0")
  if (!is.null(names(per_q_div)) && all(nzchar(names(per_q_div)))) {
    q_vals <- as.numeric(gsub("q_", "", names(per_q_div)))
  } else {
    # Default: assume q values 0.5, 1.0, 1.5, 2.0
    q_vals <- c(0.5, 1.0, 1.5, 2.0)[1:length(per_q_div)]
  }
  
  # Sort by q value to ensure correct plot ordering
  sort_idx <- order(q_vals)
  per_q_div <- per_q_div[sort_idx]
  q_vals <- q_vals[sort_idx]
  
  # Sort CI bounds if present (must happen before NA filtering)
  if (!is.null(per_q_ci) && !is.null(per_q_ci$lower) && !is.null(per_q_ci$upper)) {
    per_q_ci$lower <- per_q_ci$lower[sort_idx]
    per_q_ci$upper <- per_q_ci$upper[sort_idx]
  }
  
  # Remove NA values
  valid_idx <- !is.na(per_q_div)
  if (sum(valid_idx) == 0) {
    cat("Error: All per_q_div values are NA.\n")
    return(invisible(NULL))
  }
  
  per_q_div <- per_q_div[valid_idx]
  q_vals <- q_vals[valid_idx]
  
  # Filter CI bounds by valid indices
  if (!is.null(per_q_ci) && !is.null(per_q_ci$lower) && !is.null(per_q_ci$upper)) {
    per_q_ci$lower <- per_q_ci$lower[valid_idx]
    per_q_ci$upper <- per_q_ci$upper[valid_idx]
  }
  
  # Create ggplot2-based plot matching plot_tsallis_q_curve style
  require_pkgs(c("ggplot2"))
  
  # Prepare data frame for ggplot2
  plot_df <- data.frame(
    q = q_vals,
    divergence = as.numeric(per_q_div),
    stringsAsFactors = FALSE
  )
  
  # Add CI bounds if provided
  if (!is.null(per_q_ci) && 
      !is.null(per_q_ci$lower) && !is.null(per_q_ci$upper) &&
      length(per_q_ci$lower) == length(per_q_div) &&
      length(per_q_ci$upper) == length(per_q_div)) {
    plot_df$ci_lower <- per_q_ci$lower
    plot_df$ci_upper <- per_q_ci$upper
  } else {
    plot_df$ci_lower <- NA
    plot_df$ci_upper <- NA
  }
  
  # Create base plot with ggplot2
  # Style to match plot_tsallis_gene_profile
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::theme_minimal(base_size = 14)
  
  # Add CI ribbon if available (matching plot_tsallis_gene_profile ribbon style)
  if (!all(is.na(plot_df$ci_lower))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(x = q, ymin = ci_lower, ymax = ci_upper),
      alpha = 0.2, fill = "steelblue", color = NA, inherit.aes = FALSE
    )
  }
  
  # Main line and points (matching line width and styling)
  p <- p +
    ggplot2::geom_line(color = "steelblue", linewidth = 1.3) +
    ggplot2::geom_point(color = "steelblue", size = 3) +
    # Add reference line at q=1
    ggplot2::geom_vline(xintercept = 1, linetype = 3, color = "gray50", linewidth = 0.8)
  
  # Styling to match plot_tsallis_gene_profile
  p <- p +
    ggplot2::labs(
      title = paste0(gene_name, ": Tsallis Divergence Q-Spectrum"),
      x = "q value", 
      y = "Divergence D_q"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 16,
        margin = ggplot2::margin(b = 10)
      ),
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add CI legend line if available
  if (!all(is.na(plot_df$ci_lower))) {
    p <- p + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = "Shaded region: 95% Bootstrap CI",
      hjust = 1.05, vjust = 1.2, size = 3.5, color = "gray50"
    )
  }
  
  # Add region annotations (positioned at very top of plot area)
  p <- p +
    ggplot2::annotate("text", x = 0.5, y = Inf, label = "Rare",
                     size = 3.8, color = "darkred", fontface = "bold", vjust = 1.5, hjust = 0.5) +
    ggplot2::annotate("text", x = 1.0, y = Inf, label = "Balanced",
                     size = 3.8, color = "darkgreen", fontface = "bold", vjust = 1.5, hjust = 0.5) +
    ggplot2::annotate("text", x = 1.5, y = Inf, label = "Abundant",
                     size = 3.8, color = "darkblue", fontface = "bold", vjust = 1.5, hjust = 0.5)
  
  return(p)
}

#' Plot Q-Spectrum for a Single Gene
#'
#' Extracts per-q divergence values for a specified gene from a divergence
#' SummarizedExperiment and visualizes the q-spectrum curve using \code{plot_q_spectrum}.
#'
#' This is a convenience wrapper around the q-spectrum extraction and plotting logic,
#' designed to simplify visualization of how Tsallis divergence varies across sensitivity
#' parameters (q-values) for genes of interest identified in LMM interaction analysis.
#'
#' @param divergence_results_se A SummarizedExperiment object containing divergence
#'   results, typically output from \code{\link{calculate_divergence}}. Must have:
#'   - An assay matrix with per-q divergence estimates (columns = q-values, rows = genes)
#'   - rowData with a \code{gene_name} column (or rownames as fallback)
#'
#' @param target_gene Character string specifying the gene name to plot. Must match
#'   a value in either the \code{gene_name} column of rowData or rownames of the assay.
#'
#' @param verbose Logical; if TRUE, print diagnostic messages about gene lookup and
#'   plot generation status (default: TRUE).
#'
#' @return A `ggplot` object displaying the q-spectrum curve with point estimates,
#'   95% bootstrap confidence intervals, and region annotations. Returns NULL if plot generation fails.
#'
#' @details
#' **Gene Lookup Logic:**
#' The function searches for \code{target_gene} using the following priority:
#' 1. Exact match in \code{rowData(divergence_results_se)$gene_name} (if column exists)
#' 2. Exact match in rownames of the assay matrix (fallback)
#'
#' If no match is found, a diagnostic message lists available genes.
#'
#' **Q-Spectrum Interpretation:**
#' - **RARE_DRIVEN**: Divergence decreases with q (D(q=0.5) > D(q=2.0))
#'   → Low-abundance isoforms shift between conditions
#' - **ABUNDANT_DRIVEN**: Divergence increases with q (D(q=0.5) < D(q=2.0))
#'   → High-abundance isoforms shift between conditions
#' - **BALANCED**: Divergence relatively flat across q
#'   → All isoforms shift proportionally
#'
#' @examples
#' \dontrun{
#' # After running calculate_divergence():
#' divergence_results_se <- calculate_divergence(
#'   se = se,
#'   group_col = "sample_type",
#'   control_group = "normal",
#'   q = seq(0.1, 2, by = 0.05),
#'   bootstrap = TRUE,
#'   nboot = 100
#' )
#'
#' # Plot q-spectrum for gene LINC03040
#' plot_gene_q_spectrum(divergence_results_se, target_gene = "LINC03040")
#'
#' # Plot q-spectrum for gene CXCL12
#' plot_gene_q_spectrum(divergence_results_se, target_gene = "CXCL12", verbose = FALSE)
#' }
#'
#' @export
plot_gene_q_spectrum <- function(divergence_results_se, target_gene, verbose = TRUE) {
  
  plot_generated <- FALSE
  
  # Validate input
  if (!is(divergence_results_se, "SummarizedExperiment")) {
    stop("divergence_results_se must be a SummarizedExperiment object")
  }
  
  if (!is.character(target_gene) || length(target_gene) != 1) {
    stop("target_gene must be a single character string")
  }
  
  if (nrow(divergence_results_se) == 0) {
    if (verbose) cat("divergence_results_se is empty; no genes to plot.\n")
    return(invisible(FALSE))
  }
  
  # Extract components
  div_rd <- as.data.frame(rowData(divergence_results_se))
  div_assay <- assay(divergence_results_se)
  
  # Get gene names (prioritize gene_name column, fall back to rownames)
  div_gene_names <- if ("gene_name" %in% colnames(div_rd)) {
    div_rd$gene_name
  } else {
    rownames(div_assay)
  }
  
  # Search for target gene
  gene_idx <- which(div_gene_names == target_gene)
  
  if (length(gene_idx) == 0) {
    if (verbose) {
      cat(sprintf("Gene '%s' not found in divergence results.\n", target_gene))
      cat("Available genes (first 10):", paste(head(div_gene_names, 10), collapse = ", "), "...\n")
    }
    return(invisible(FALSE))
  }
  
  # Validate assay dimensions
  if (nrow(div_assay) == 0 || ncol(div_assay) < 4) {
    if (verbose) cat("Assay has insufficient q-values (<4); cannot plot spectrum.\n")
    return(invisible(FALSE))
  }
  
  # Extract and plot
  plot_obj <- tryCatch({
    per_q_divs <- div_assay[gene_idx[1], ]
    names(per_q_divs) <- colnames(div_assay)
    
    if (verbose) {
      cat(sprintf("[DEBUG plot_gene_q_spectrum] Column names in assay: %s\n", 
                  paste(colnames(div_assay), collapse=", ")))
      cat(sprintf("[DEBUG plot_gene_q_spectrum] Raw divergence values: %s\n", 
                  paste(round(per_q_divs, 4), collapse=", ")))
    }
    
    # Extract bootstrap CI bounds from rowData
    per_q_ci <- NULL
    ci_lower_cols <- grep("^lower_ci_q", colnames(div_rd), value = TRUE)
    ci_upper_cols <- grep("^upper_ci_q", colnames(div_rd), value = TRUE)
    
    if (length(ci_lower_cols) > 0 && length(ci_upper_cols) > 0) {
      ci_lower_bounds <- as.numeric(div_rd[gene_idx[1], ci_lower_cols])
      ci_upper_bounds <- as.numeric(div_rd[gene_idx[1], ci_upper_cols])
      
      if (length(ci_lower_bounds) == length(per_q_divs) && 
          length(ci_upper_bounds) == length(per_q_divs)) {
        per_q_ci <- list(lower = ci_lower_bounds, upper = ci_upper_bounds)
      }
    }
    
    # Call plot function and return the ggplot2 object
    p <- plot_q_spectrum(per_q_divs, gene_name = target_gene, per_q_ci = per_q_ci)
    plot_generated <- TRUE
    
    if (verbose) cat(sprintf("✓ Q-spectrum plot generated for gene '%s'\n", target_gene))
    
    return(p)
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("Error plotting q-spectrum for '%s': %s\n", target_gene, e$message))
    }
    return(NULL)
  })
  
  return(plot_obj)
}


#' Plot Q-Spectrum Curves for Multiple Top Genes
#'
#' Creates a multi-panel grid comparing per-q divergence profiles across the top
#' N genes identified by LMM interaction analysis. Each panel shows the full q-spectrum
#' divergence curve with the gene name and adjusted p-value in the title.
#'
#' @param eff_res Output from \code{\link{effect_sizes_divergence}} OR \code{NULL}.
#'   If provided, must contain `$interaction_results` with columns: gene, adj_p_interaction, per_q_pattern.
#'   If \code{NULL}, uses fallback with lm_res + divergence_results_se.
#'
#' @param lm_res (Optional) Data frame from LMM analysis with columns: gene, adj_p_interaction.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param divergence_results_se (Optional) SummarizedExperiment from \code{\link{calculate_divergence}}.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param n_genes Integer; number of top genes to plot (default: 9). Genes are sorted by
#'   increasing adjusted p-value (most significant first).
#'
#' @param ncol Integer; number of columns in grid layout (default: 3). Number of rows is
#'   automatically calculated as ceiling(n_genes / ncol).
#'
#' @param verbose Logical; if TRUE, print diagnostic messages (default: TRUE).
#'
#' @return A ggplot2 object created via \code{patchwork} combining all gene panels,
#'   or NULL if gene data is unavailable. The function automatically handles ggplot2 grid
#'   creation and returns a print-ready object.
#'
#' @details
#' **Input Modes:**
#' - **Mode 1 (Primary)**: Pass eff_res directly (from effect_sizes_divergence)
#' - **Mode 2 (Fallback)**: Pass lm_res + divergence_results_se instead
#'
#' **Gene Filtering:**
#' Genes are ranked by decreasing statistical significance (increasing adj_p_interaction).
#' Only genes with complete per-q divergence data are included. If fewer than n_genes
#' have valid data, the function returns all available genes.
#'
#' **Plot Features:**
#' - Title shows: gene name and adjusted p-value (q-value format)
#' - Per-q divergence curve with point estimates and 95% bootstrap CI bands
#' - Vertical reference line at q=1 (Kullback-Leibler divergence point)
#' - Region labels: "Rare Isoforms" (q<1), "Balanced" (q≈1), "Abundant Isoforms" (q>1)
#' - All plots use consistent ggplot2 styling matching plot_q_spectrum
#'
#' @examples
#' \dontrun{
#' # Mode 1: Using eff_res
#' p <- plot_multi_gene_q_spectrum(eff_res = eff_res, n_genes = 9)
#' print(p)
#'
#' # Mode 2: Using lm_res + divergence_results_se
#' p <- plot_multi_gene_q_spectrum(
#'   lm_res = lm_res,
#'   divergence_results_se = divergence_results_se,
#'   n_genes = 6, ncol = 2
#' )
#' print(p)
#' }
#'
#' @export
plot_multi_gene_q_spectrum <- function(eff_res = NULL, 
                                        lm_res = NULL, 
                                        divergence_results_se = NULL,
                                        n_genes = 9, 
                                        ncol = 3, 
                                        verbose = TRUE) {
  
  require_pkgs(c("ggplot2", "patchwork", "SummarizedExperiment"))
  
  # ============================================================================
  # Step 1: Extract gene data from either eff_res or (lm_res + divergence_results_se)
  # ============================================================================
  
  genes_to_plot <- NULL
  per_q_patterns <- NULL
  adj_p_values <- NULL
  
  # Mode 1: Try eff_res first
  if (!is.null(eff_res) && is.list(eff_res)) {
    if (!is.null(eff_res$interaction_results) && nrow(eff_res$interaction_results) > 0) {
      int_res <- eff_res$interaction_results
      
      # Check for required columns (flexible on p-value column name)
      has_gene <- "gene" %in% colnames(int_res)
      has_per_q <- "per_q_pattern" %in% colnames(int_res)
      has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
      has_p_raw <- "p_value_interaction" %in% colnames(int_res)
      has_p_any <- has_p_adj || has_p_raw
      
      # Determine which p-value column to use
      p_col <- NA_character_
      if (has_p_adj) {
        p_col <- "adj_p_interaction"
      } else if (has_p_raw) {
        p_col <- "p_value_interaction"
      }
      
      if (has_gene && has_per_q && has_p_any) {
        # Sort by significance (lowest p-value first = most significant)
        int_res_sorted <- int_res[order(int_res[[p_col]], na.last = TRUE), ]
        
        # Get top n_genes
        int_res_subset <- head(int_res_sorted, n_genes)
        
        # Check if per_q_pattern has valid data
        valid_patterns <- !is.na(int_res_subset$per_q_pattern) & 
                          int_res_subset$per_q_pattern != "" &
                          int_res_subset$per_q_pattern != "NA"
        
        if (any(valid_patterns)) {
          genes_to_plot <- int_res_subset$gene[valid_patterns]
          per_q_patterns <- int_res_subset$per_q_pattern[valid_patterns]
          adj_p_values <- int_res_subset[[p_col]][valid_patterns]
          
          if (verbose) cat(sprintf("[plot_multi_gene_q_spectrum] Mode 1: Using eff_res with %s column (%d valid genes)\n", p_col, length(genes_to_plot)))
        } else {
          if (verbose) cat("[plot_multi_gene_q_spectrum] Mode 1 failed: per_q_pattern values are empty or invalid\n")
        }
      } else {
        if (verbose) {
          cat("[plot_multi_gene_q_spectrum] Mode 1 failed: Missing required columns\n")
          cat("  - has 'gene':", has_gene, "\n")
          cat("  - has 'per_q_pattern':", has_per_q, "\n")
          cat("  - has 'adj_p_interaction':", has_p_adj, "\n")
          cat("  - has 'p_value_interaction':", has_p_raw, "\n")
        }
      }
    } else {
      if (verbose) cat("[plot_multi_gene_q_spectrum] Mode 1 failed: eff_res$interaction_results is NULL or empty\n")
    }
  } else {
    if (verbose) cat("[plot_multi_gene_q_spectrum] Mode 1 failed: eff_res is NULL or not a list\n")
  }
  
  # Mode 2: Fallback to lm_res + divergence_results_se
  if (is.null(genes_to_plot) && !is.null(lm_res) && !is.null(divergence_results_se)) {
    if (nrow(lm_res) > 0 && nrow(divergence_results_se) > 0) {
      # Check required columns
      if (all(c("gene", "adj_p_interaction") %in% colnames(lm_res))) {
        div_rd <- as.data.frame(rowData(divergence_results_se))
        div_assay <- assay(divergence_results_se)
        
        # Get gene names from divergence rowData
        div_gene_names <- if ("gene_name" %in% colnames(div_rd)) {
          div_rd$gene_name
        } else {
          rownames(div_assay)
        }
        
        if (length(div_gene_names) > 0 && nrow(div_assay) > 0) {
          # Sort lm_res by significance
          lm_sorted <- lm_res[order(lm_res$adj_p_interaction, na.last = TRUE), ]
          top_genes <- head(lm_sorted$gene, n_genes)
          
          # Get indices in divergence_results_se
          gene_indices <- match(top_genes, div_gene_names)
          valid_idx <- !is.na(gene_indices)
          valid_genes <- top_genes[valid_idx]
          
          if (length(valid_genes) > 0) {
            genes_to_plot <- valid_genes
            adj_p_values <- lm_sorted$adj_p_interaction[1:length(valid_genes)]
            
            # Extract per_q patterns from assay
            per_q_patterns <- character(length(valid_genes))
            for (i in seq_along(valid_genes)) {
              gene_idx <- which(div_gene_names == valid_genes[i])[1]
              if (!is.na(gene_idx)) {
                divs <- div_assay[gene_idx, ]
                per_q_patterns[i] <- paste(divs[!is.na(divs)], collapse = ",")
              }
            }
            
            if (verbose) cat("[plot_multi_gene_q_spectrum] Mode 2 (fallback): Using lm_res + divergence_results_se\n")
          }
        }
      }
    }
  }
  
  # ============================================================================
  # Step 2: Validate and parse gene data
  # ============================================================================
  
  if (is.null(genes_to_plot) || length(genes_to_plot) == 0) {
    if (verbose) {
      cat("ERROR: No valid genes to plot. Check input data:\n")
      cat("  - eff_res provided:", !is.null(eff_res), "\n")
      if (!is.null(eff_res)) {
        cat("  - eff_res$interaction_results exists:", !is.null(eff_res$interaction_results), "\n")
        if (!is.null(eff_res$interaction_results)) {
          cat("  - Number of rows:", nrow(eff_res$interaction_results), "\n")
          cat("  - Has 'per_q_pattern' column:", "per_q_pattern" %in% colnames(eff_res$interaction_results), "\n")
        }
      }
      cat("  - lm_res provided:", !is.null(lm_res), "\n")
      cat("  - divergence_results_se provided:", !is.null(divergence_results_se), "\n")
    }
    return(invisible(NULL))
  }
  
  n_genes_actual <- length(genes_to_plot)
  if (verbose) cat(sprintf("Plotting %d genes in %d-column grid\n", n_genes_actual, ncol))
  
  # ============================================================================
  # Step 3: Create individual q-spectrum plots for each gene
  # ============================================================================
  
  plot_list <- list()
  
  for (i in seq_along(genes_to_plot)) {
    gene_name <- genes_to_plot[i]
    pattern_str <- per_q_patterns[i]
    adj_p <- adj_p_values[i]
    
    # Parse per-q pattern string
    tryCatch({
      per_q_vals <- as.numeric(strsplit(pattern_str, ",")[[1]])
      
      if (length(per_q_vals) == 0 || all(is.na(per_q_vals))) {
        if (verbose) cat(sprintf("  Skipping %s: no valid per-q values\n", gene_name))
        next
      }
      
      # Create plot using plot_q_spectrum logic
      q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
      
      plot_df <- data.frame(
        q = q_vals,
        divergence = per_q_vals,
        stringsAsFactors = FALSE
      )
      
      # Create base plot
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::geom_line(color = "darkblue", linewidth = 1.2) +
        ggplot2::geom_point(color = "darkblue", size = 2.8) +
        ggplot2::geom_vline(xintercept = 1, linetype = 3, color = "gray60", linewidth = 0.8, alpha = 0.7) +
        ggplot2::labs(
          title = sprintf("%s", gene_name),
          subtitle = sprintf("adj p = %.2e", adj_p),
          x = "q (Tsallis parameter)",
          y = "Tsallis divergence D_q"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5, face = "bold", size = 12,
            margin = ggplot2::margin(b = 3)
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5, size = 10, color = "gray40",
            margin = ggplot2::margin(b = 8)
          ),
          plot.margin = ggplot2::margin(t = 8, b = 8, l = 6, r = 6),
          panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.25),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size = 10),
          axis.title = ggplot2::element_text(size = 10, face = "plain")
        )
      
      plot_list[[i]] <- p
      
    }, error = function(e) {
      if (verbose) cat(sprintf("  Error plotting %s: %s\n", gene_name, e$message))
    })
  }
  
  # ============================================================================
  # Step 4: Combine plots into grid using patchwork
  # ============================================================================
  
  if (length(plot_list) == 0) {
    if (verbose) {
      cat("ERROR: No valid plots were created.\n")
      cat("This may occur if:\n")
      cat("  - per_q_pattern values cannot be parsed as numeric comma-separated strings\n")
      cat("  - All genes had parsing errors in tryCatch blocks\n")
      cat("  - Sample size or q-value count was too small\n")
    }
    return(invisible(NULL))
  }
  
  nrow <- ceiling(length(plot_list) / ncol)
  
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow)
  
  if (verbose) cat(sprintf("✓ Multi-gene q-spectrum plot created with %d genes\n", length(plot_list)))
  
  return(combined_plot)
}

#' Plot Tsallis Divergence Profiles Across q-Spectrum for Top Genes
#'
#' Visualizes how Tsallis divergence D_q varies across the q-spectrum for selected genes.
#' This reveals which diversity scales (rare vs. abundant isoforms) drive the observed
#' biological differences between groups.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with multiple q-values.
#' @param gene Optional character vector of gene names to plot. Can be a single gene name,
#'   a vector of names, or NULL. If NULL, uses top genes from lm_res.
#' @param lm_res Optional data.frame with columns: gene, p_value_interaction, 
#'   adj_p_lmm (or adj_p_interaction). If provided and gene=NULL, top genes are selected by significance.
#' @param readcounts Optional matrix of raw read counts (genes * transcripts) for computing
#'   true Tsallis divergence from isoform distributions. If NULL, uses entropy-based approximation.
#' @param tx2gene_map Optional data.frame mapping transcripts to genes (columns: "transcript", "gene").
#' @param group_col Character name of the column in colData(se) indicating group assignment.
#'   Default: "group".
#' @param n_top Integer. When using lm_res, plot top n_top genes. Default: 3.
#' @param assay_name Character name of the assay containing diversity measures. Default: "diversity".
#' @param arrange_type How to arrange multiple plots. Options: "facet" (faceted grid),
#'   "list" (named list). Default: "facet".
#' @param signed Logical. If TRUE, compute signed divergence (mean2 - mean1); 
#'   positive = group2 higher, negative = group1 higher. If FALSE (default), 
#'   absolute divergence. Signed divergence reveals directional differences.
#'
#' @return If arrange_type="facet", returns a single ggplot with facets by gene.
#'   If arrange_type="list", returns a named list of ggplot objects (one per gene).
#'   When signed=TRUE, negative values indicate group1 preference, positive indicate group2 preference.
#'
#' @details
#' The plot shows:
#' - X-axis: q value (from 0.1 to ~3, depending on SE)
#' - Y-axis: Tsallis divergence D_q between control and treatment groups
#' - Shape: Profile reveals which q-ranges drive divergence:
#'   - Low q (<1): rare isoform divergence
#'   - q≈1: KL divergence region
#'   - High q (>1): dominant isoform divergence
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#'   theme element_text scale_color_manual scale_size_manual
#' @importFrom dplyr group_by summarise
#' @importFrom SummarizedExperiment colData assay
#'
#' @examples
#' \dontrun{
#'   # Example 1: Plot with signed divergence (default, shows directionality)
#'   p <- plot_tsallis_divergence_profile(ts_se, gene = c("LINC03040", "PNPT1", "CXCL12"))
#'   print(p)
#'
#'   # Example 2: Plot with absolute divergence (magnitude only)
#'   p_abs <- plot_tsallis_divergence_profile(ts_se, gene = "LINC03040", signed = FALSE)
#'   print(p_abs)
#'
#'   # Example 3: Plot top 5 significant genes from LMM results
#'   p_list <- plot_tsallis_divergence_profile(ts_se, lm_res = interaction_results, 
#'                                             n_top = 5, arrange_type = "list")
#'   print(p_list$LINC03040)
#'
#'   # Example 4: Single gene with detailed inspection (signed)
#'   p_single <- plot_tsallis_divergence_profile(ts_se, gene = "LINC03040", 
#'                                               arrange_type = "list", signed = TRUE)
#'   print(p_single$LINC03040)
#' }
#'
#' @keywords internal
#' @noRd
plot_tsallis_divergence_profile <- function(se,
                                            gene = NULL,
                                            lm_res = NULL,
                                            readcounts = NULL,
                                            tx2gene_map = NULL,
                                            group_col = "group",
                                            n_top = 3,
                                            assay_name = "diversity",
                                            arrange_type = c("facet", "list"),
                                            signed = TRUE) {
    require_pkgs(c("ggplot2", "dplyr", "SummarizedExperiment"))
    arrange_type <- match.arg(arrange_type)

    # Validate SE
    if (!inherits(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment")
    if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
        stop("Assay '", assay_name, "' not found in se")
    }

    # Determine which genes to plot
    if (is.null(gene)) {
        if (is.null(lm_res)) stop("Either 'gene' or 'lm_res' must be provided")
        if (!is.data.frame(lm_res)) stop("'lm_res' must be a data.frame")
        if (!("gene" %in% colnames(lm_res))) stop("'lm_res' must have a 'gene' column")
        
        # Use FDR-adjusted p-value if available, else raw p-value
        p_col <- if ("adj_p_lmm" %in% colnames(lm_res)) "adj_p_lmm" 
                 else if ("adj_p_interaction" %in% colnames(lm_res)) "adj_p_interaction"
                 else if ("p_value_interaction" %in% colnames(lm_res)) "p_value_interaction"
                 else stop("'lm_res' must contain p-value column")
        
        genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[p_col]])]))
        genes <- head(genes_ordered, n_top)
    } else {
        genes <- as.character(unlist(gene))
    }

    if (length(genes) == 0) stop("No genes selected for plotting")

    # Extract metadata
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    if (!group_col %in% colnames(col_data)) {
        stop("Column '", group_col, "' not found in colData(se)")
    }

    # Parse column names to extract q values
    col_names <- colnames(se)
    extract_q <- function(name) {
        if (grepl("_q=", name)) {
            as.numeric(gsub(".*_q=", "", name))
        } else {
            NA
        }
    }
    q_values <- sapply(col_names, extract_q)
    unique_q <- sort(unique(q_values[!is.na(q_values)]))

    if (length(unique_q) < 2) {
        stop("SE must contain multiple q-values in column names (format: *_q=0.5)")
    }

    # Get group information
    groups <- unique(col_data[[group_col]])
    if (length(groups) != 2) {
        stop("Exactly 2 groups required in '", group_col, "' column; found: ", 
             paste(groups, collapse = ", "))
    }

    # Helper: calculate divergence for a single gene at a specific q
    calc_div_for_gene_q <- function(gene_name, q_val) {
        # Get columns matching this q value
        cols_q <- which(q_values == q_val)
        if (length(cols_q) == 0) return(NA)

        # Extract entropy data for this gene at this q
        diversity_matrix <- SummarizedExperiment::assay(se, assay_name)
        if (!gene_name %in% rownames(diversity_matrix)) return(NA)

        entropy_vals <- diversity_matrix[gene_name, cols_q]
        group_vals <- col_data[[group_col]][cols_q]

        # Compute divergence using calculate_tsallis_divergence_paired_gene if readcounts available
        if (!is.null(readcounts) && !is.null(tx2gene_map)) {
            # TIER 1: True Tsallis divergence from isoform distributions
            tryCatch({
                # Prepare gene subset data
                gene_cols <- cols_q
                entropy_pred <- data.frame(
                    entropy_pred = entropy_vals,
                    q = q_val,
                    group = group_vals,
                    stringsAsFactors = FALSE
                )

                div <- calculate_tsallis_divergence_paired_gene(
                    gene_name = gene_name,
                    gene_data = data.frame(
                        entropy = entropy_vals,
                        q = q_val,
                        group = group_vals,
                        sample = colnames(se)[cols_q],
                        stringsAsFactors = FALSE
                    ),
                    readcounts = readcounts,
                    tx2gene_map = tx2gene_map,
                    entropy_pred = entropy_pred,
                    group_levels = groups
                )
                # Return signed or absolute divergence based on parameter
                return(if (signed) div else abs(div))
            }, error = function(e) {
                return(NA)
            })
        }

        # TIER 2: Entropy-based approximation (fallback)
        # Divergence ≈ difference in mean entropy between groups
        group1_vals <- entropy_vals[group_vals == groups[1]]
        group2_vals <- entropy_vals[group_vals == groups[2]]

        if (length(group1_vals) == 0 || length(group2_vals) == 0) return(NA)

        mean1 <- mean(group1_vals, na.rm = TRUE)
        mean2 <- mean(group2_vals, na.rm = TRUE)
        
        # Return signed or absolute divergence based on parameter
        if (signed) {
            div_approx <- mean2 - mean1  # Signed: positive if group2 > group1
        } else {
            div_approx <- max(0, abs(mean1 - mean2))  # Unsigned: always non-negative
        }

        return(div_approx)
    }

    # Calculate divergence for all genes * q combinations
    plot_data_list <- list()
    for (gene_name in genes) {
        divergences <- sapply(unique_q, function(q) calc_div_for_gene_q(gene_name, q))
        df_gene <- data.frame(
            gene = gene_name,
            q = unique_q,
            divergence = divergences,
            stringsAsFactors = FALSE
        )
        plot_data_list[[gene_name]] <- df_gene
    }

    all_plot_data <- do.call(rbind, plot_data_list)
    rownames(all_plot_data) <- NULL

    # Remove NA divergences
    all_plot_data <- all_plot_data[!is.na(all_plot_data$divergence), ]

    if (nrow(all_plot_data) == 0) {
        stop("No valid divergence values computed; check readcounts/tx2gene_map or data structure")
    }

    # Add direction indicator for signed divergence
    if (signed) {
        all_plot_data$direction <- ifelse(all_plot_data$divergence > 0, 
                                         paste0("Positive: ", groups[2], " higher"),
                                         paste0("Negative: ", groups[1], " higher"))
    }

    # Build plot
    if (signed) {
        # Signed divergence plot with color-coded direction
        p <- ggplot2::ggplot(all_plot_data, ggplot2::aes(x = q, y = divergence, color = direction)) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::scale_color_manual(
                name = "Divergence Direction:",
                values = setNames(c("#E63946", "#1D3557"), 
                                 c(paste0("Negative: ", groups[1], " higher"),
                                   paste0("Positive: ", groups[2], " higher")))
            ) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile: Directional (Signed)",
                x = "q value (diversity scale parameter)",
                y = "Tsallis Divergence D_q (Positive = Right Group Higher, Negative = Left Group Higher)"
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
                legend.position = "right"
            )
    } else {
        # Absolute divergence plot (original)
        p <- ggplot2::ggplot(all_plot_data, ggplot2::aes(x = q, y = divergence, color = gene)) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile Across q-Spectrum",
                x = "q value (diversity scale parameter)",
                y = "Tsallis Divergence D_q (Absolute)",
                color = "Gene"
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
                legend.position = "right"
            )
    }

    # Conditionally apply faceting
    if (arrange_type == "facet" && length(unique(all_plot_data$gene)) > 1) {
        if (signed) {
            p <- p + ggplot2::facet_wrap(~gene, scales = "free_y") +
                ggplot2::theme(legend.position = "top")
        } else {
            p <- p + ggplot2::facet_wrap(~gene, scales = "free_y") +
                ggplot2::theme(legend.position = "top")
        }
    }

    # Return as list if requested
    if (arrange_type == "list") {
        plots <- list()
        for (gene_name in genes) {
            df_gene <- all_plot_data[all_plot_data$gene == gene_name, ]
            if (nrow(df_gene) > 0) {
                if (signed) {
                    p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence, fill = direction, color = direction)) +
                        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
                        ggplot2::geom_line(linewidth = 1.2) +
                        ggplot2::geom_point(size = 3.5, alpha = 0.8) +
                        ggplot2::scale_color_manual(
                            name = "Divergence Direction:",
                            values = setNames(c("#E63946", "#1D3557"), 
                                             c(paste0("Negative: ", groups[1], " higher"),
                                               paste0("Positive: ", groups[2], " higher")))
                        ) +
                        ggplot2::scale_fill_manual(
                            name = "Divergence Direction:",
                            values = setNames(c("#E63946", "#1D3557"), 
                                             c(paste0("Negative: ", groups[1], " higher"),
                                               paste0("Positive: ", groups[2], " higher")))
                        ) +
                        ggplot2::labs(
                            title = paste("Divergence Profile (Signed):", gene_name),
                            x = "q value",
                            y = "Tsallis Divergence D_q",
                            subtitle = paste0("Red = ", groups[1], " higher | Blue = ", groups[2], " higher")
                        ) +
                        ggplot2::theme_minimal(base_size = 12) +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                            panel.grid.minor = ggplot2::element_blank(),
                            legend.position = "right"
                        )
                } else {
                    p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence)) +
                        ggplot2::geom_line(color = "#2E86AB", linewidth = 1.2) +
                        ggplot2::geom_point(color = "#2E86AB", size = 3.5, alpha = 0.8) +
                        ggplot2::labs(
                            title = paste("Divergence Profile:", gene_name),
                            x = "q value",
                            y = "Tsallis Divergence D_q (Absolute)"
                        ) +
                        ggplot2::theme_minimal(base_size = 12) +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                            panel.grid.minor = ggplot2::element_blank()
                    )
                }
                plots[[gene_name]] <- p_gene
            }
        }
        return(plots)
    }

    return(p)
}

#' Plot Global Divergence q-Curve Across All Genes
#'
#' Visualizes the average (mean/median) Tsallis divergence D_q across all genes
#' as a function of q-value. This provides a **global view** of which diversity scales
#' (rare vs. abundant isoforms) drive the most divergence on average across the dataset,
#' complementing gene-specific divergence profiles.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with multiple q-values.
#' @param readcounts Optional matrix of raw read counts (genes * transcripts) for computing
#'   true Tsallis divergence from isoform distributions. If NULL, uses entropy-based approximation.
#' @param tx2gene_map Optional data.frame mapping transcripts to genes (columns: "transcript", "gene").
#' @param group_col Character name of the column in colData(se) indicating group assignment.
#'   Default: "group".
#' @param assay_name Character name of the assay containing diversity measures. Default: "diversity".
#' @param metric Character. Summary statistic to display: "mean" or "median". Default: "median".
#' @param variability_metric Character. Error bar type: "sd" (standard deviation) or "iqr" (interquartile range).
#'   Default: "iqr".
#'
#' @return A single `ggplot` object showing the divergence q-curve.
#'
#' @details
#' The plot shows:
#' - X-axis: q value (from 0.1 to ~3, depending on SE)
#' - Y-axis: Average divergence D_q across all genes
#' - Ribbon: Variability bands (±1 SD or ±IQR/2) around the central estimate
#' - Shape: Global divergence profile reveals:
#'   - Low q (<1): Rare isoform divergence dominates
#'   - q≈1: KL divergence region
#'   - High q (>1): Dominant isoform divergence dominates
#'   - Flat profile: Uniform divergence across scales
#'
#' **Interpretation**: Compare with `plot_tsallis_q_curve` (entropy) to understand
#' the relationship between entropy changes and divergence patterns.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal
#'   element_text scale_color_manual
#' @importFrom dplyr group_by summarise
#' @importFrom SummarizedExperiment colData assay
#'
#' @examples
#' \dontrun{
#'   # Plot global divergence q-curve for all genes
#'   p <- plot_divergence_q_curve(ts_se)
#'   print(p)
#'
#'   # Compare entropy vs divergence patterns
#'   p_entropy <- plot_tsallis_q_curve(ts_se)
#'   p_div <- plot_divergence_q_curve(ts_se)
#'   gridExtra::grid.arrange(p_entropy, p_div, ncol = 2)
#' }
#'
#' @export
plot_divergence_q_curve <- function(se,
                                    readcounts = NULL,
                                    tx2gene_map = NULL,
                                    group_col = "group",
                                    assay_name = "diversity",
                                    metric = c("median", "mean"),
                                    variability_metric = c("iqr", "sd")) {
    require_pkgs(c("ggplot2", "dplyr", "SummarizedExperiment"))
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)

    # Validate SE
    if (!inherits(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment")
    if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
        stop("Assay '", assay_name, "' not found in se")
    }

    # Extract metadata
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    if (!group_col %in% colnames(col_data)) {
        stop("Column '", group_col, "' not found in colData(se)")
    }

    # Parse column names to extract q values
    col_names <- colnames(se)
    extract_q <- function(name) {
        if (grepl("_q=", name)) {
            as.numeric(gsub(".*_q=", "", name))
        } else {
            NA
        }
    }
    q_values <- sapply(col_names, extract_q)
    unique_q <- sort(unique(q_values[!is.na(q_values)]))

    if (length(unique_q) < 2) {
        stop("SE must contain multiple q-values in column names (format: *_q=0.5)")
    }

    # Get group information
    groups <- unique(col_data[[group_col]])
    if (length(groups) != 2) {
        stop("Exactly 2 groups required in '", group_col, "' column; found: ", 
             paste(groups, collapse = ", "))
    }

    # Helper: calculate divergence for a single gene at a specific q
    calc_div_for_gene_q <- function(gene_name, q_val) {
        cols_q <- which(q_values == q_val)
        if (length(cols_q) == 0) return(NA)

        diversity_matrix <- SummarizedExperiment::assay(se, assay_name)
        if (!gene_name %in% rownames(diversity_matrix)) return(NA)

        entropy_vals <- diversity_matrix[gene_name, cols_q]
        group_vals <- col_data[[group_col]][cols_q]

        # TIER 1: True Tsallis divergence from isoform distributions
        if (!is.null(readcounts) && !is.null(tx2gene_map)) {
            tryCatch({
                entropy_pred <- data.frame(
                    entropy_pred = entropy_vals,
                    q = q_val,
                    group = group_vals,
                    stringsAsFactors = FALSE
                )

                div <- calculate_tsallis_divergence_paired_gene(
                    gene_name = gene_name,
                    gene_data = data.frame(
                        entropy = entropy_vals,
                        q = q_val,
                        group = group_vals,
                        sample = colnames(se)[cols_q],
                        stringsAsFactors = FALSE
                    ),
                    readcounts = readcounts,
                    tx2gene_map = tx2gene_map,
                    entropy_pred = entropy_pred,
                    group_levels = groups
                )
                return(abs(div))
            }, error = function(e) {
                return(NA)
            })
        }

        # TIER 2: Entropy-based approximation (fallback)
        group1_vals <- entropy_vals[group_vals == groups[1]]
        group2_vals <- entropy_vals[group_vals == groups[2]]

        if (length(group1_vals) == 0 || length(group2_vals) == 0) return(NA)

        mean1 <- mean(group1_vals, na.rm = TRUE)
        mean2 <- mean(group2_vals, na.rm = TRUE)
        div_approx <- max(0, abs(mean1 - mean2))

        return(div_approx)
    }

    # Calculate divergence for all genes * q combinations
    all_genes <- rownames(SummarizedExperiment::assay(se, assay_name))
    divergence_data <- list()

    for (gene_name in all_genes) {
        divergences <- sapply(unique_q, function(q) calc_div_for_gene_q(gene_name, q))
        for (i in seq_along(unique_q)) {
            if (!is.na(divergences[i])) {
                divergence_data[[length(divergence_data) + 1]] <- data.frame(
                    gene = gene_name,
                    q = unique_q[i],
                    divergence = divergences[i],
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(divergence_data) == 0) {
        stop("No valid divergence values computed; check data structure or readcounts")
    }

    all_data <- do.call(rbind, divergence_data)
    rownames(all_data) <- NULL

    # Compute central tendency and variability at each q
    if (variability_metric == "iqr") {
        summary_stats <- all_data %>%
            dplyr::group_by(q) %>%
            dplyr::summarise(
                central = if (metric == "median") median(divergence, na.rm = TRUE) else mean(divergence, na.rm = TRUE),
                spread = stats::IQR(divergence, na.rm = TRUE),
                .groups = "drop"
            )
        spread_factor <- 1/2  # IQR/2 for symmetric ribbon
        spread_label <- "IQR"
    } else {  # sd
        summary_stats <- all_data %>%
            dplyr::group_by(q) %>%
            dplyr::summarise(
                central = if (metric == "median") median(divergence, na.rm = TRUE) else mean(divergence, na.rm = TRUE),
                spread = sqrt(stats::var(divergence, na.rm = TRUE)),
                .groups = "drop"
            )
        spread_factor <- 1  # ±1 SD
        spread_label <- "SD"
    }

    # Create plot
    metric_label <- if (metric == "median") "Median" else "Mean"
    p <- ggplot2::ggplot(summary_stats, ggplot2::aes(x = q, y = central)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = central - spread * spread_factor, 
                        ymax = central + spread * spread_factor),
            alpha = 0.25,
            fill = "#2E86AB",
            color = NA
        ) +
        ggplot2::geom_line(color = "#2E86AB", linewidth = 1.3) +
        ggplot2::geom_point(color = "#2E86AB", size = 3.5, alpha = 0.8) +
        ggplot2::labs(
            title = "Global Divergence q-Curve: Average Divergence Across All Genes",
            x = "q value (diversity scale parameter)",
            y = "Tsallis Divergence D_q",
            subtitle = paste0(metric_label, " ± ", spread_label, " across ", length(unique(all_data$gene)), " genes")
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "italic"),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}

