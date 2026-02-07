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
            "mean"
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

#' Plot q-curve profile for a single gene comparing groups
#'
#' For a selected gene, plot per-sample Tsallis entropy across q values and
#' overlay per-group median +/- IQR ribbons so group-level differences are
#' easy to compare. Expects a `SummarizedExperiment` produced by
#' `calculate_diversity()` with `_q=` suffixes in column names.
#'
#' @param se A `SummarizedExperiment` from `calculate_diversity()`.
#' @param gene Character scalar; gene symbol to plot.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Column name in `colData(se)` with sample type labels
#'   (default: "sample_type"). If missing, a single-group fallback is used.
#' @param show_samples Logical; if TRUE, draw per-sample lines in the
#'   background (default: FALSE).
#' @return A `ggplot` object showing the gene q-curve profile by group.
#' @export
plot_tsallis_gene_profile <- function(se,
                                     gene,
                                     assay_name = "diversity",
                                     sample_type_col = "sample_type",
                                     show_samples = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    require_pkgs(c("dplyr", "tidyr", "SummarizedExperiment"))

    long <- prepare_tsallis_long(se, assay_name = assay_name, sample_type_col = sample_type_col)
    if (!("Gene" %in% colnames(long))) stop("prepare_tsallis_long did not return Gene column")
    # filter for requested gene
    sel <- as.character(gene)
    long_g <- long[as.character(long$Gene) == sel, , drop = FALSE]
    if (nrow(long_g) == 0) stop("Gene not found in assay: ", sel)

    # numeric q for plotting
    long_g$qnum <- as.numeric(as.character(long_g$q))

    # per-group statistics
    stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
                                 median = median(tsallis, na.rm = TRUE),
                                 IQR = stats::IQR(tsallis, na.rm = TRUE), .groups = "drop")

    p <- ggplot2::ggplot() + ggplot2::theme_minimal()

    if (isTRUE(show_samples)) {
        p <- p + ggplot2::geom_line(data = long_g, ggplot2::aes(x = qnum, y = tsallis, group = sample, color = group), alpha = 0.25)
    }

    p <- p +
        ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum, ymin = median - IQR/2, ymax = median + IQR/2, fill = group), alpha = 0.2, inherit.aes = FALSE) +
        ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = median, color = group), linewidth = 1.1) +
        ggplot2::labs(title = paste0(sel, ": Tsallis entropy q-curve profile"), x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
        ggplot2::scale_color_discrete(name = "Group") + ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, 
                                                          margin = ggplot2::margin(b = 10)))

    return(p)
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
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
#' plot_diversity_density(se)
plot_diversity_density <- function(
  se,
  assay_name = "diversity",
  sample_type_col = NULL
) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr", "SummarizedExperiment"))
    long <- get_assay_long(
        se,
        assay_name = assay_name,
        value_name = "diversity",
        sample_type_col = sample_type_col
    )

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = diversity, group = sample, color = sample_type)
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_grid(. ~ sample_type) +
        ggplot2::scale_color_manual(values = c("black", "darkorchid4")) +
        ggplot2::guides(color = "none") +
        ggplot2::theme_minimal() +
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
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
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
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Samples", y = "Diversity") +
        ggplot2::scale_fill_viridis_d(name = "Group")
}


#' Plot MA plot for difference or linear model results
#'
#' Creates an MA plot (Mean vs log-fold-change) that works for both differential
#' analysis and linear model interaction results.
#'
#' @param x Data.frame from `calculate_difference()`, `calculate_fc()`, or
#' `calculate_lm_interaction()`. For differential analysis, should contain mean
#' columns, a fold-change column, and optionally adjusted p-values. For linear
#' model results, should contain estimate/interaction columns.
#'
#' @param mean_cols Optional character vector of length 2 naming the mean
#' columns. Defaults to the first two columns that end with `_mean` or `_median`.
#' Ignored for linear model input if x-axis mapping is not applicable.
#'
#' @param fold_col Name of the fold-change column (default: `log2_fold_change`).
#' Also searches for `logFC`, `estimate_interaction`, etc.
#'
#' @param padj_col Name of the adjusted p-value column (default:
#' `adjusted_p_values`). Also searches for `adj_p_interaction`, etc.
#'
#' @param sig_alpha Threshold for significance (default: 0.05).
#'
#' @param fc_df Optional data.frame from `calculate_fc()` containing pre-computed
#' fold changes from a different data source. If provided for differential analysis,
#' fold changes are extracted from this frame instead of `x`. The first column should
#' contain gene identifiers. Use this to plot fold changes from read counts with
#' p-values from diversity analysis, or vice versa.
#'
#' @param diff_res Optional data.frame from `calculate_difference()`. Used with
#' linear model results (`x`) to obtain fold changes and mean values for plotting.
#' Allows plotting linear model p-values with differential fold changes.
#'
#' @param x_label Optional label for x-axis (overrides automatic detection).
#'
#' @param y_label Optional label for y-axis (overrides automatic detection).
#'
#' @return A `ggplot` MA-plot object.
#'
#' @export
#'
#' @details Automatically detects the type of input data (differential vs linear model)
#' and creates appropriate MA plot. Can also combine results from different workflows,
#' e.g., plotting linear model p-values with differential fold changes, or read count
#' fold changes with diversity-based p-values.
#'
#' @examples
#' # Differential analysis MA plot using diversity fold changes
#' df <- data.frame(
#'     gene = paste0("g", seq_len(10)),
#'     sampleA_mean = runif(10),
#'     sampleB_mean = runif(10),
#'     log2_fold_change = rnorm(10),
#'     adjusted_p_values = runif(10)
#' )
#'
#' Create an MA (mean-average) plot showing log fold change vs mean diversity.
#' This function is designed for single-q differential analysis results only.
#'
#' @param x Data.frame from `test_differential()` containing adjusted p-values
#'   and fold changes.
#' @param fc_df Optional data.frame with fold changes from alternative methods
#'   (e.g., read count fold changes). Should have 'genes' and 'log2_fold_change'
#'   columns. If provided, will override fold changes in `x`.
#' @param sig_alpha Significance threshold for p-values (default: 0.05).
#' @param x_label Custom x-axis label (optional).
#' @param y_label Custom y-axis label (optional).
#' @param title Custom plot title (optional). If NULL, generates default title.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' # Example with test_differential results
#' # res <- test_differential(ts_se, sample_type_col = "sample_type")
#' # plot_ma(res)
plot_ma <- function(
    x,
    fc_df = NULL,
    sig_alpha = 0.05,
    x_label = NULL,
    y_label = NULL,
    title = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required")
    }

    df <- as.data.frame(x)

    # Merge with alternative fold changes if provided
    if (!is.null(fc_df)) {
        fc_df <- as.data.frame(fc_df)
        
        # Ensure gene column exists in both dataframes
        if (!("genes" %in% colnames(df)) && !is.null(rownames(df))) {
            df$genes <- rownames(df)
        }
        if (!("genes" %in% colnames(fc_df)) && ncol(fc_df) > 0) {
            colnames(fc_df)[1] <- "genes"
        }
        
        # Merge on genes column
        if ("genes" %in% colnames(df) && "genes" %in% colnames(fc_df)) {
            df <- merge(df, fc_df[, c("genes", "log2_fold_change")],
                       by = "genes", all.x = TRUE, suffixes = c("", ".fc"))
            if ("log2_fold_change.fc" %in% colnames(df)) {
                df$log2_fold_change <- df$log2_fold_change.fc
                df$log2_fold_change.fc <- NULL
            }
        }
    }

    # Extract required columns
    cn <- colnames(df)

    # Find mean/median columns (for x-axis)
    mean_cols <- grep("_mean$|_median$", cn, value = TRUE)
    if (length(mean_cols) < 2) {
        stop("Input must have at least 2 mean or median columns (e.g., Normal_mean, Tumor_mean)")
    }
    
    # Check that we have either two _mean columns or two _median columns (not mixed)
    mean_type_cols <- grep("_mean$", cn, value = TRUE)
    median_type_cols <- grep("_median$", cn, value = TRUE)
    has_means <- length(mean_type_cols) >= 2
    has_medians <- length(median_type_cols) >= 2
    
    if (!has_means && !has_medians) {
        stop("Could not find two mean or two median columns")
    }
    
    # Use mean columns if available, otherwise use median columns
    if (has_means) {
        mean_cols <- mean_type_cols
    } else {
        mean_cols <- median_type_cols
    }
    
    df$mean <- rowMeans(df[, mean_cols[1:2], drop = FALSE], na.rm = TRUE)

    # Find p-value column
    padj_col <- grep("adjusted_p_values|adj_p_value|adj_p", cn, value = TRUE)[1]
    if (is.na(padj_col)) {
        stop("Could not find adjusted p-value column")
    }
    df$padj <- df[[padj_col]]

    # Find fold change column
    fold_col <- grep("log2_fold_change|log_fold_change|fold_change", cn, value = TRUE)[1]
    if (is.na(fold_col)) {
        stop("Could not find fold change column")
    }
    df$fold <- df[[fold_col]]

    # Convert log2 to log10 if necessary
    if (grepl("log2", fold_col, ignore.case = TRUE)) {
        df$fold <- df$fold / log2(10)
    }

    # Mark significant genes
    df$significant <- ifelse(
        !is.na(df$padj) & df$padj < sig_alpha,
        "significant", "non-significant"
    )

    # Remove rows with missing values
    df <- df[is.finite(df$mean) & is.finite(df$fold), ]

    if (nrow(df) == 0) {
        stop("No valid points to plot")
    }

    # Set labels
    x_label_use <- x_label %||% "Mean diversity"
    
    # Auto-detect y-label based on whether fc_df was provided
    if (is.null(y_label)) {
        if (!is.null(fc_df)) {
            y_label_use <- "Log10 fold change (read counts)"
        } else {
            y_label_use <- "Log10 fold change (Tsallis entropy)"
        }
    } else {
        y_label_use <- y_label
    }
    
    # Auto-detect title based on whether fc_df was provided
    if (is.null(title)) {
        if (!is.null(fc_df)) {
            title_use <- "Ma plot: read count fold change"
        } else {
            title_use <- "Ma plot: tsallis entropy fold change"
        }
    } else {
        title_use <- title
    }

    # Create plot with title and theme
    ggplot2::ggplot(df, ggplot2::aes(x = mean, y = fold, color = significant)) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::scale_color_manual(
            values = c("non-significant" = "black", "significant" = "red"),
            guide = "none"
        ) +
        ggplot2::theme_minimal() +
            ggplot2::labs(
                title = title_use,
                x = x_label_use,
                y = y_label_use
            ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain",
                                                          margin = ggplot2::margin(b = 10)))
}

# Helper for default values
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Plot median +- IQR of Tsallis entropy across q values by group
#'
#' This reproduces the `tsallis-q-curve-mean-sd` plot from the vignette: for
#' each q value, compute per-gene Tsallis entropy per sample, summarize across
#' genes by group (median and IQR) and plot median with a ribbon spanning
#' median +- IQR/2.
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Column name in `colData(se)` containing sample
#'   type labels (default: "sample_type").
#' @return A `ggplot` object showing median +- IQR across q values by group.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:40, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:40]
#' se <- calculate_diversity(rc, gs,
#'    q = seq(0.01, 0.1, by = 0.03), norm = TRUE)
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
        # long$q may be a factor; convert to numeric for plotting
        long$qnum <- as.numeric(as.character(long$q))
        stats_df <- dplyr::summarise(dplyr::group_by(long, group, qnum),
            median = median(tsallis, na.rm = TRUE),
            IQR = stats::IQR(tsallis, na.rm = TRUE), .groups = "drop"
        )
            p <- ggplot2::ggplot(
                stats_df,
                ggplot2::aes(x = qnum, y = median, color = group, fill = group)
            ) +
                ggplot2::geom_line(linewidth = 1) +
                ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = median - IQR / 2, ymax = median + IQR / 2),
                    alpha = 0.2,
                    color = NA
                ) +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    title = "Global Tsallis q-curve: entropy across all genes",
                    x = "q value",
                    y = y_label,
                    color = "Group",
                    fill = "Group"
                ) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 14))
        # use default discrete ggplot2 colours (not viridis)
        p <- p + ggplot2::scale_color_discrete(name = "Group") +
            ggplot2::scale_fill_discrete(name = "Group")
        # If there is only a single group present, hide the legend/Group label
        if (length(unique(stats_df$group)) == 1) {
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
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
#' plot_tsallis_violin_multq(se)
plot_tsallis_violin_multq <- function(se, assay_name = "diversity") {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    long <- prepare_tsallis_long(se, assay_name = assay_name)

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = q, y = tsallis, fill = group)
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
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::labs(
            title = "Violin Plot: Tsallis Entropy Distribution Across Multiple q Values",
            x = "q value",
            y = "Tsallis entropy",
            fill = "Group"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain",
                                                          margin = ggplot2::margin(b = 10)))
}
#' Density plot of Tsallis entropy for multiple q values

#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with
#' multiple q values (column names contain `_q=`).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @return A `ggplot` density plot object faceted by q and colored by group.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
#' plot_tsallis_density_multq(se)
plot_tsallis_density_multq <- function(se, assay_name = "diversity") {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    long <- prepare_tsallis_long(se, assay_name = assay_name)

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = tsallis, color = group, fill = group)
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_wrap(~q, scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_discrete(name = "Group") +
        ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::labs(
            title = "Density Plot: Tsallis Entropy Distribution Across Multiple q Values",
            x = "Tsallis entropy",
            y = "Density",
            color = "Group",
            fill = "Group"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain",
                                                          margin = ggplot2::margin(b = 10)))
}


#' Volcano plot for differential results
#' @param diff_df Data.frame from `calculate_difference()` containing at least
#' `mean_difference` and an adjusted p-value column (default
#' `adjusted_p_values`).
#' @param x_col Column name for x-axis (default `mean_difference`).
#' @param padj_col Column name for adjusted p-values (default
#' `adjusted_p_values`).
#' @param label_thresh Absolute x threshold to mark significance (default 0.1).
#' @param padj_thresh Adjusted p-value cutoff (default 0.05).
#' @param top_n Integer; number of top genes to annotate by smallest adjusted
#' p-value (default: 5).
#' @return ggplot volcano plot.
#' @export
#' @examples
#' df <- data.frame(
#'     gene = paste0("g", seq_len(10)),
#'     mean_difference = runif(10),
#'     adjusted_p_values = runif(10)
#' )
#' plot_volcano(df, x_col = "mean_difference", padj_col = "adjusted_p_values")
#' Volcano plot for single-q differential analysis
#'
#' Create a volcano plot showing fold change vs adjusted p-value significance.
#' This function is designed for single-q differential analysis results only.
#' Works with output from `test_differential()` or `calculate_difference()`.
#'
#' @param diff_df Data.frame from `test_differential()` or similar differential
#'   analysis results. Should contain p-value and optionally fold-change columns.
#' @param x_col Column name for x-axis values. If not specified, will auto-detect
#'   from available columns (e.g., "mean_difference", "median_difference").
#' @param padj_col Column name for adjusted p-values (default: "adjusted_p_values").
#' @param label_thresh Threshold for fold-change labeling (default: 0.1).
#' @param padj_thresh P-value threshold for significance (default: 0.05).
#' @param top_n Number of top significant genes to label (default: 5).
#' @param title Custom plot title (optional). If NULL, generates default title.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' # res <- test_differential(ts_se, sample_type_col = "sample_type")
#' # plot_volcano(res)
plot_volcano <- function(
    diff_df,
    x_col = NULL,
    padj_col = "adjusted_p_values",
    label_thresh = 0.1,
    padj_thresh = 0.05,
    top_n = 5,
    title = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required")
    }

    df <- as.data.frame(diff_df)
    cn <- colnames(df)

    # Auto-detect x-axis column if not specified
    if (is.null(x_col)) {
        # Look for difference columns
        diff_cols <- grep("_difference$", cn, value = TRUE, ignore.case = TRUE)
        if (length(diff_cols) > 0) {
            x_col <- diff_cols[1]
        } else {
            # If no difference column found, use the first numeric column (excluding p-values)
            numeric_cols <- sapply(df, is.numeric)
            p_cols <- grep("p_value|p.value", cn, ignore.case = TRUE)
            numeric_cols[p_cols] <- FALSE
            if (any(numeric_cols)) {
                x_col <- cn[which(numeric_cols)[1]]
            } else {
                stop("Could not find suitable column for x-axis. Specify 'x_col' explicitly.")
            }
        }
    }

    # Verify x-axis column exists
    if (!(x_col %in% cn)) {
        stop("Column '", x_col, "' not found in diff_df")
    }

    # Verify p-value column exists
    if (!(padj_col %in% cn)) {
        stop("Column '", padj_col, "' not found in diff_df")
    }

    # Extract and convert values
    df$xval <- as.numeric(df[[x_col]])
    df$padj <- as.numeric(df[[padj_col]])

    # Replace invalid p-values with 1 (no significance)
    df$padj[is.na(df$padj)] <- 1
    df$padj[df$padj <= 0] <- .Machine$double.xmin

    # Mark significant points
    df$significant <- ifelse(
        abs(df$xval) >= label_thresh & df$padj < padj_thresh,
        "significant",
        "non-significant"
    )

    # Remove rows with missing values
    df <- df[is.finite(df$xval) & is.finite(df$padj), ]

    if (nrow(df) == 0) {
        stop("No valid points to plot")
    }

    # Determine metric label
    metric_label <- if (grepl("median", x_col, ignore.case = TRUE)) {
        "Median"
    } else if (grepl("mean", x_col, ignore.case = TRUE)) {
        "Mean"
    } else {
        "Value"
    }

    # Generate default title if not provided (sentence-case)
    title_use <- title %||% "Volcano plot"

    # Format x-axis label: remove underscores and capitalize first letter
    x_label_formatted <- gsub("_", " ", x_col)
    x_label_formatted <- paste0(toupper(substr(x_label_formatted, 1, 1)), 
                                 substr(x_label_formatted, 2, nchar(x_label_formatted)))

    # Create base plot with title and theme
    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = xval, y = -log10(padj), color = significant)
    ) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::scale_color_manual(
            values = c("non-significant" = "black", "significant" = "red"),
            guide = "none"
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(padj_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-label_thresh, label_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::theme_minimal() +
        # Format padj label: remove underscores and capitalize first letter
        {
            padj_label_formatted <- gsub("_", " ", padj_col)
            padj_label_formatted <- paste0(toupper(substr(padj_label_formatted, 1, 1)),
                                           substr(padj_label_formatted, 2, nchar(padj_label_formatted)))
            ggplot2::labs(
                title = title_use,
                x = x_label_formatted,
                y = paste0("-Log10(", padj_label_formatted, ")")
            )
        } +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain",
                                                          margin = ggplot2::margin(b = 10)))

    return(p)
}

#' Plot top transcripts for a gene   #' For a given gene, find transcripts using
#' a tx->gene mapping, compute per-  Plot top transcripts for a gene   For a
#' given gene, find transcripts using a tx->gene mapping, compute per-
#' transcript statistics between two sample groups, select the top N transcripts
#' by p-value and plot their expression across groups.
#' @param counts Matrix or data.frame of transcript counts.
#' Rows are transcripts and columns are samples.
#' @param gene Character; gene symbol to inspect.
#' @param samples Character vector of sample group labels (length =
#' ncol(counts)).
#' @param tx2gene Path or data.frame mapping transcripts to genes.
#' Must contain columns `Transcript` and `Gen`.
#' @param top_n Integer number of transcripts to show (default = 3).
#' Use NULL to plot all transcripts for the gene.
#' @param pseudocount Numeric pseudocount added before log2
#' (default = 1e-6) to avoid division by zero.
#' @param output_file Optional file path to save the plot.
#' If `NULL`, the `ggplot` object is returned.
#' @param metric Aggregation metric used to summarize transcript expression
#'   per group when plotting. One of c("median", "mean", "variance").
#'   Defaults to "median" to preserve previous behavior.
#' @return A `ggplot` object (or invisibly saved file if `output_file`
#' provided).
#' @importFrom utils read.delim
#' @export
#' @name plot_top_transcripts
#' @examples
#' tx_counts <- matrix(
#'     sample(1:100, 24, replace = TRUE),
#'     nrow = 6
#' )
#' rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
#' colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
#'
#' tx2gene <- data.frame(
#'     Transcript = rownames(tx_counts),
#'     Gen = rep(paste0("G", seq_len(3)), each = 2),
#'     stringsAsFactors = FALSE
#' )
#'
#' samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
#'
#' plot_top_transcripts(
#'     tx_counts,
#'     gene = c("G1", "G2"),
#'     samples = samples,
#'     tx2gene = tx2gene,
#'     top_n = 2
#' )
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "tx",
        "expr",
        "group",
        "tx_cond",
        "sample",
        "log2expr"
    ))
}

#' Helper to compute fill limits across multiple genes
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

#' Helper to draw grid layout with title, plots, and legend using base grid
.draw_transcript_grid <- function(grobs, title, legend_grob, n, heights, to_file = NULL) {
    grid::grid.newpage()
    grid::pushViewport(
        grid::viewport(layout = grid::grid.layout(3, n, heights = heights))
    )
    # Title row
    vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(n))
    grid::pushViewport(vp_title)
    grid::grid.text(title, x = 0.5, gp = grid::gpar(fontsize = 12))
    grid::upViewport()
    # Plot rows
    for (i in seq_along(grobs)) {
        vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = i)
        grid::pushViewport(vp)
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }
    # Legend row
    if (!is.null(legend_grob)) {
        vp_leg <- grid::viewport(layout.pos.row = 3, layout.pos.col = seq_len(n))
        grid::pushViewport(vp_leg)
        grid::grid.draw(legend_grob)
        grid::upViewport()
    }
    grid::upViewport()
    if (!is.null(to_file)) grDevices::dev.off()
}

plot_top_transcripts <- function(
    counts,
    gene,
    samples,
    tx2gene = NULL,
    top_n = 3,
    pseudocount = 1e-6,
    output_file = NULL,
    metric = c("median", "mean", "variance")
) {
    # Early return if no genes provided
    if (is.null(gene) || length(gene) == 0 || all(is.na(gene))) {
        message("No genes provided; skipping plot.")
        return(invisible(NULL))
    }
    
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop(
            "`counts` must be a matrix or data.frame with transcripts as rownames"
        )
    }
    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) {
        stop(
            "`counts` must have rownames corresponding to transcript identifiers"
        )
    }
    if (length(samples) != ncol(counts)) {
        stop(
            "Length of `samples` must equal number of columns in `counts`"
        )
    }

    # tx2gene must be supplied as a path or data.frame
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) stop("tx2gene file not found: ", tx2gene)
        mapping <- utils::read.delim(
            tx2gene,
            stringsAsFactors = FALSE,
            header = TRUE
        )
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame")
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
        stop(
            "tx2gene must have columns 'Transcript' and 'Gen'"
        )
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop(
            "ggplot2 required for plotting"
        )
    }

    if (length(samples) != ncol(counts)) {
        stop(
            "Length of `samples` must equal number of columns in `counts'"
        )
    }

    # normalize and validate aggregation argument; define aggregation function
    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice,
                      median = function(x) stats::median(x, na.rm = TRUE),
                      mean = function(x) base::mean(x, na.rm = TRUE),
                      variance = function(x) stats::var(x, na.rm = TRUE)
    )
    # Title label: descriptive title including the chosen metric
    agg_label <- sprintf("Transcript-level expression with metric %s", metric_choice)
    # increment call counter (kept for internal tracking) but do not show it in title
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label

    make_plot_for_gene <- function(gene_single, fill_limits = NULL) {
        txs <- mapping$Transcript[mapping$Gen == gene_single]
        txs <- intersect(txs, rownames(counts))
        if (length(txs) == 0) stop("No transcripts found for gene: ", gene_single)
        if (!is.null(top_n)) txs <- head(txs, top_n)

        mat <- counts[txs, , drop = FALSE]
        df_all <- as.data.frame(mat)
        df_all$tx <- rownames(mat)
        df_long <- tidyr::pivot_longer(df_all,
            -tx,
            names_to = "sample",
            values_to = "expr"
        )
        df_long$group <- rep(samples, times = length(txs))

        # summarize per transcript x group using selected agg function
        df_summary <- aggregate(expr ~ tx + group,
            data = df_long,
            FUN = agg_fun
        )
        df_summary$log2expr <- log2(df_summary$expr + pseudocount)
        df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))

        # show transcripts on y (readable labels) and groups on x
        p <- ggplot2::ggplot(
            df_summary,
            ggplot2::aes(
                x = group,
                y = tx,
                fill = log2expr
            )
        ) +
            ggplot2::geom_tile(color = "white", width = 0.95, height = 0.95) +
            ggplot2::scale_fill_viridis_c(
                option = "viridis",
                direction = -1,
                na.value = "grey80",
                limits = fill_limits
            ) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::labs(
                title = agg_label_unique,
                x = NULL,
                y = NULL,
                fill = "log2(expr)"
            ) +
            ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 8),
                axis.text.x = ggplot2::element_text(size = 8),
                axis.ticks = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(size = 14, hjust = 0.6),
                legend.position = "bottom",
                legend.key.width = ggplot2::unit(1.2, "cm"),
                plot.margin = ggplot2::margin(4, 4, 4, 4)
            ) +
            ggplot2::guides(
                fill = ggplot2::guide_colorbar(
                    title.position = "top",
                    barwidth = 6,
                    barheight = 0.35
                )
            )

        return(p)
    }

    # Produce plots (single or multiple). Do not save inside helper - save once
    # below.
    if (length(gene) > 1) {
        fill_limits <- .compute_transcript_fill_limits(gene, mapping, counts, samples, top_n, agg_fun, pseudocount)

        plots <- lapply(seq_along(gene), function(i) {
            gname <- gene[i]
            pp <- make_plot_for_gene(gname,
                fill_limits = fill_limits
            )
            # set per-gene title with only gene name
            per_gene_title <- if (!is.na(gname) && nzchar(as.character(gname))) {
                as.character(gname)
            } else {
                ""
            }
            pp <- pp + ggplot2::labs(title = per_gene_title) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.55, size = 14))
            pp
        })

        # try patchwork first (collect guides), otherwise cowplot with shared legend
        if (requireNamespace("patchwork", quietly = TRUE)) {
            combined <- Reduce(`+`, plots) +
                patchwork::plot_layout(nrow = 1, guides = "collect") &
                ggplot2::theme(legend.position = "bottom")
            # add a single centered title (sentence case)
            combined <- combined + patchwork::plot_annotation(title = agg_label_unique,
                theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.7, size = 14, 
                                                                            margin = ggplot2::margin(b = 10)))
            )
            result_plot <- combined
        } else if (requireNamespace("cowplot", quietly = TRUE)) {
            # extract legend from first plot
            p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
            legend <- cowplot::get_legend(p_for_legend)
            plots_nolegend <- lapply(
                plots,
                function(pp) pp + ggplot2::theme(legend.position = "none")
            )
            grid <- cowplot::plot_grid(
                plotlist = plots_nolegend,
                nrow = 1,
                align = "hv"
            )
            # title as a separate grob on top
            title_grob <- cowplot::ggdraw() + cowplot::draw_label(agg_label_unique, fontface = 'plain', x = 0.7, hjust = 0.5, size = 14)
            result_plot <- cowplot::plot_grid(title_grob, grid, legend, ncol = 1, rel_heights = c(0.08, 1, 0.08))
        } else {
            # fallback: arrange grobs horizontally using base grid (no extra packages)
            plots_nolegend <- lapply(
                plots,
                function(pp) pp + ggplot2::theme(legend.position = "none")
            )
            grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)
            # extract legend from original first plot
            g_full <- ggplot2::ggplotGrob(plots[[1]])
            legend_idx <- which(vapply(
                g_full$grobs,
                function(x) x$name,
                character(1)
            ) == "guide-box")
            if (length(legend_idx)) {
                legend_grob <- g_full$grobs[[legend_idx[1]]]
            } else {
                legend_grob <- NULL
            }
            n <- length(grobs)
            heights <- grid::unit.c(
                grid::unit(0.6, "cm"),
                grid::unit(1, "null"),
                grid::unit(0.7, "cm")
            )

            if (!is.null(output_file)) {
                png(filename = output_file, width = 800 * n, height = 480, res = 150)
                .draw_transcript_grid(grobs, agg_label_unique, legend_grob, n, heights, to_file = output_file)
                return(invisible(NULL))
            } else {
                .draw_transcript_grid(grobs, agg_label_unique, legend_grob, n, heights)
                return(invisible(NULL))
            }
        }
    } else {
        result_plot <- make_plot_for_gene(gene)
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    } else {
        return(result_plot)
    }
}

