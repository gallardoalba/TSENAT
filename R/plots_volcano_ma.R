#' Combine Volcano and MA-Tsallis Plots in a Grid Layout
#'
#' Creates a side-by-side grid layout with a volcano plot on the left and an
#' MA-Tsallis plot on the right.
#' Both plots are generated from differential analysis results data.
#'
#' @param diff_df Data.frame from differential analysis containing required
#' columns for both volcano and MA plots.
#' @param x_col Column name for x-axis in volcano plot (e.g.,
#' 'mean_difference'). Auto-detected if NULL.
#' @param padj_col Column name for adjusted p-values (default: 'padj').
#' @param label_thresh Threshold for volcano plot labels (default: 0.1).
#' @param sig_alpha Numeric significance threshold for adjusted p-values
#' (default: 0.05).
#' @param top_n Number of top genes to annotate in volcano plot (default: 5).
#' @param title_volcano Title for volcano plot. If NULL, auto-generated.
#' @param title_ma Title for MA plot (default: 'Tsallis-based MA plot').
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A `ggplot2` object showing a 1x2 grid with volcano plot on the
#' left and MA plot on the right.
#'
#' @examples
#' # Simulate differential analysis results
#' x <- data.frame(
#'   genes = paste0('g', seq_len(20)),
#'   mean_difference = rnorm(20, sd = 1),
#'   padj = runif(20, 1e-5, 0.1),
#'   log2_fold_change = rnorm(20, sd = 0.8)
#' )
#' # Placeholder: actual usage would require valid differential results
#' # .plot_diversity_volcano_ma(x, sig_alpha = 0.05)
#'
#' @noRd
.plot_diversity_volcano_ma <- function(diff_df, x_col = NULL, padj_col = "padj", label_thresh = 0.1,
    sig_alpha = 0.05, top_n = 5, title_volcano = NULL, title_ma = "Tsallis-based MA plot",
    ...) {
    # Require cowplot for grid arrangement
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("cowplot package required for .plot_diversity_volcano_ma()")
    }

    # Create volcano plot
    p_volcano <- .plot_volcano(diff_df = diff_df, x_col = x_col, padj_col = padj_col,
        label_thresh = label_thresh, sig_alpha = sig_alpha, top_n = top_n, title = title_volcano)

    # Create MA plot
    p_ma <- .plot_ma_tsallis(x = diff_df, sig_alpha = sig_alpha, title = title_ma,
        ...)

    # Arrange plots side by side: volcano on left, MA on right
    grid <- cowplot::plot_grid(p_volcano, p_ma, nrow = 1, ncol = 2, align = "h",
        axis = "b")

    return(grid)
}

# Core MA plotting implementation used by wrappers above. Accepts a
# differential results `x` (data.frame) and an optional `fc_df` with
# fold-changes (genes as rownames or a `genes` column). Returns a `ggplot`
# MA-plot.
#' Core MA plotting implementation (internal)
#'
#' This is an internal helper used by `.plot_ma_tsallis()`.
#' It is documented here for developers but is not exported.
#' @noRd
.plot_ma_core <- function(x, fc_df = NULL, diff_res = NULL, sig_alpha = 0.05, x_label = NULL,
    y_label = NULL, title = NULL, ...) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    # Ensure gene identifier column exists
    if (!("genes" %in% colnames(df))) {
        if ("gene_id" %in% colnames(df)) {
            df$genes <- df$gene_id
        } else if (!is.null(rownames(df))) {
            df$genes <- rownames(df)
        }
    }

    # If external fc_df provided, merge fold values
    if (!is.null(fc_df)) {
        fdf <- as.data.frame(fc_df, stringsAsFactors = FALSE)
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
        df <- merge(df, fdf[, c("genes", "log2_fold_change")], by = "genes", all.x = TRUE,
            suffixes = c("", ".fc"))
        if ("log2_fold_change.fc" %in% colnames(df))
            df$log2_fold_change <- ifelse(!is.na(df$log2_fold_change.fc), df$log2_fold_change.fc,
                df$log2_fold_change)
    }

    # Use helper for fold/mean column detection
    fold_col_candidates <- c("log2_fold_change", "logFC", "fold", "estimate_interaction",
        "fold_change")
    fold_col <- intersect(fold_col_candidates, colnames(df))
    if (length(fold_col) == 0)
        stop("Could not find a fold-change column in input")
    fold_col <- fold_col[1]

    # Detect p-value column for significance flagging
    padj_candidates <- c("padj", "adjusted_p_values", "adj_p_value", "adj_p", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col))
        padj_col[1] else NULL
    padj <- if (!is.null(padj_col))
        as.numeric(df[[padj_col]]) else rep(1, nrow(df))
    padj[is.na(padj)] <- 1

    # Validate mean/median column consistency
    mean_cols <- grep("_mean$", colnames(df), ignore.case = TRUE, value = TRUE)
    median_cols <- grep("_median$", colnames(df), ignore.case = TRUE, value = TRUE)

    if (length(mean_cols) > 0 && length(median_cols) > 0) {
        stop("Could not find two mean or two median columns - found both mean and median columns. ",
            "Ensure input contains either mean columns (e.g., A_mean, B_mean) OR median columns (e.g., A_median, B_median), not both.")
    }

    if (length(mean_cols) > 0 && length(mean_cols) < 2) {
        stop("Could not find two mean or two median columns - found ", length(mean_cols),
            " mean column(s). ", "Ensure input contains at least two mean columns (e.g., A_mean, B_mean).")
    }

    if (length(median_cols) > 0 && length(median_cols) < 2) {
        stop("Could not find two mean or two median columns - found ", length(median_cols),
            " median column(s). ", "Ensure input contains at least two median columns (e.g., A_median, B_median).")
    }

    # Determine which columns to use for mean calculation
    mean_cols_to_use <- if (length(mean_cols) > 0)
        mean_cols else if (length(median_cols) > 0)
        median_cols else NULL

    # Prepare MA plot data with label formatting
    prep <- .prepare_ma_plot_df(df, fold_col = fold_col, mean_cols = mean_cols_to_use,
        x_label = x_label, y_label = y_label)
    plot_df <- prep$plot_df
    plot_df$padj <- padj[match(plot_df$genes, df$genes)]
    plot_df$significant <- ifelse(abs(plot_df$y) > 0 & plot_df$padj < sig_alpha,
        "significant", "non-significant")

    # Format labels
    x_label_formatted <- .format_label(prep$x_label)
    y_label_raw <- prep$y_label %||% fold_col
    y_label_formatted <- .format_label(y_label_raw)
    if (!is.null(y_label_formatted)) {
        y_label_formatted <- sub("\\blog2\\b", "log10", y_label_formatted, ignore.case = TRUE)
    }

    # Build plot with significance coloring
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = significant)) +
        ggplot2::geom_point(alpha = 0.75, size = 3.2) + ggplot2::scale_color_manual(values = .significance_colors(),
        guide = "none") + ggplot2::labs(x = x_label_formatted, y = y_label_formatted)

    # Apply publication theme and settings
    p <- .apply_publication_theme(p, title = title %||% "MA plot: mean vs log10 fold-change",
        base_size = 11) + ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"))

    p
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
#' @param padj_col Adjusted p-value column name (default: 'padj').
#' @param label_thresh Fold-change threshold used to annotate points
#' (default: 0.1).
#' @param sig_alpha Adjusted p-value cutoff for significance (default: 0.05).
#' @param top_n Number of top significant genes to label (default: 5).
#' @param title Optional plot title; if `NULL` a default title is used.
#'
#' @return A `ggplot2` object.
#' @noRd

.plot_volcano <- function(diff_df, x_col = NULL, padj_col = "padj", label_thresh = 0.1,
    sig_alpha = 0.05, top_n = 5, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required")
    }

    prep_volcano <- .prepare_volcano_df(diff_df = diff_df, x_col = x_col, padj_col = padj_col,
        label_thresh = label_thresh, sig_alpha = sig_alpha, title = title)
    df <- prep_volcano$df
    x_col <- prep_volcano$x_col
    padj_col <- prep_volcano$padj_col
    x_label_formatted <- prep_volcano$x_label_formatted
    padj_label_formatted <- prep_volcano$padj_label_formatted
    title_use <- prep_volcano$title_use

    p <- ggplot2::ggplot(df, ggplot2::aes(x = xval, y = -log10(padj), color = significant)) +
        ggplot2::geom_point(alpha = 0.75, size = 3.4) + ggplot2::scale_color_manual(values = .significance_colors(),
        guide = "none")

    # Add reference lines using Phase 5 helper
    p <- .add_reference_lines(p, h_intercept = -log10(sig_alpha), v_intercept = c(-label_thresh,
        label_thresh), h_color = "gray50", v_color = "gray50")

    p <- .apply_publication_theme(p, title = title_use, base_size = 11) + ggplot2::labs(x = x_label_formatted,
        y = paste0("-Log10(", padj_label_formatted, ")"))

    p
}

# Core MA plotting implementation documentation moved to internal block
#' Plot MA using Tsallis-based fold changes
#'
#' Wrapper around `plot_ma(..., type = 'tsallis')` for convenience and
#' clearer API separation.
#'
#' @param x Data.frame from `.calculate_difference()`.
#' @param sig_alpha Numeric significance threshold for adjusted p-values
#' (default: 0.05).
#' @param x_label Optional x-axis label passed to `plot_ma`.
#' @param y_label Optional y-axis label passed to `plot_ma`.
#' @param title Optional plot title passed to `plot_ma`.
#' @param ... Additional arguments passed to `plot_ma()`.
#' @return A `ggplot2` object representing the MA plot.
#' @noRd

.plot_ma_tsallis <- function(x, sig_alpha = 0.05, x_label = NULL, y_label = NULL,
    title = NULL, ...) {
    title_use <- title %||% "Tsallis-based MA plot"
    x_label_use <- x_label %||% "mean_difference"
    y_label_use <- y_label %||% "Log10 fold-change of entropy"
    .plot_ma_core(x, fc_df = NULL, sig_alpha = sig_alpha, x_label = x_label_use,
        y_label = y_label_use, title = title_use)
}
