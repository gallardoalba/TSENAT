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
        require_pkgs(c("tidyr", "dplyr", "SummarizedExperiment"))
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
        require_pkgs(c("tidyr", "dplyr", "SummarizedExperiment"))
        mat <- SummarizedExperiment::assay(se, assay_name)
        if (is.null(mat)) stop("Assay not found: ", assay_name)
        df <- as.data.frame(mat)
        genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) {
            SummarizedExperiment::rowData(se)$genes
        } else {
            rownames(df)
        }
        df <- cbind(df, Gene = genes_col)

        long <- tidyr::pivot_longer(df,
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


#' Plot MA plot for difference results
#' @param diff_df Data.frame from `calculate_difference()` containing mean
#' columns, a `log2_fold_change` column, and `adjusted_p_values`.
#' @param mean_cols Optional character vector of length 2 naming the mean
#' columns.  Defaults to the first two columns that end with `_mean`.
#' @param fold_col Name of the fold-change column (default: `log2_fold_change`).
#' @param padj_col Name of the adjusted p-value column (default:
#' `adjusted_p_values`).
#' @param sig_alpha Threshold for significance (default: 0.05).
#' @return A `ggplot` MA-plot object.
#' @export
#' @examples
#' df <- data.frame(
#'     gene = paste0("g", seq_len(10)),
#'     sampleA_mean = runif(10),
#'     sampleB_mean = runif(10),
#'     log2_fold_change = rnorm(10),
#'     adjusted_p_values = runif(10)
#' )
#' plot_ma(df)
plot_ma <- function(
  diff_df,
  mean_cols = NULL,
  fold_col = "log2_fold_change",
  padj_col = "adjusted_p_values",
  sig_alpha = 0.05
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

    if (is.null(mean_cols)) {
        mean_only <- grep("_mean$", colnames(diff_df), value = TRUE)
        med_only <- grep("_median$", colnames(diff_df), value = TRUE)
        if (length(mean_only) == 2) {
            mean_cols <- mean_only[1:2]
        } else if (length(med_only) == 2) {
            mean_cols <- med_only[1:2]
        } else {
            stop("Could not find two mean or two median columns in 'diff_df'")
        }
    }

    df <- diff_df
    df$mean <- apply(
        df[, mean_cols],
        1,
        function(x) mean(as.numeric(x), na.rm = TRUE)
    )
    # Determine metric label (mean vs median) based on selected columns
    metric_label <- if (all(grepl("_mean$", mean_cols))) {
        "Mean"
    } else if (all(grepl("_median$", mean_cols))) {
        "Median"
    } else {
        "Value"
    }
    df$padj <- if (padj_col %in% colnames(df)) df[[padj_col]] else NA
    df$fold <- if (fold_col %in% colnames(df)) df[[fold_col]] else NA
    df$significant <- ifelse(
        !is.na(df$padj) & df$padj < sig_alpha,
        "significant",
        "non-significant"
    )

    # Convert fold-change to log10 if input is log2-based (common output
    # from calculate_difference/calculate_fc)
    if (any(!is.na(df$fold)) && grepl("log2", fold_col, ignore.case = TRUE)) {
        df$fold <- df$fold / log2(10)
    }

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = mean,
            y = fold,
            color = significant
        )
    ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c(
            "non-significant" = "black",
            "significant" = "red"
        )) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = paste0(metric_label, " diversity"),
            y = "Log10 fold change")

    return(p)
}


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
#' @param y Character; which statistic to plot on the y-axis. One of
#'   `"S"` (Tsallis entropy, default) or `"D"` (Hill numbers). When
#'   `"D"` is requested and a `hill` assay is not present, the function
#'   will attempt to convert from un-normalized Tsallis `S_q` to `D_q`.
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
    sample_type_col = "sample_type",
    y = c("S", "D")
) {
    y <- match.arg(y)
        # SE-first API: require a SummarizedExperiment with per-column sample
        # type mapping in `colData(se)[, sample_type_col]` (or allow a single
        # group dataset where `sample_type` is omitted).
        if (inherits(se, "SummarizedExperiment")) {
        require_pkgs(c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment"))
                # If user requests Hill numbers (D), prefer an existing 'hill' assay; otherwise use requested assay and convert if possible.
                if (y == "D" && ("hill" %in% SummarizedExperiment::assayNames(se))) {
                    long <- prepare_tsallis_long(se, assay_name = "hill", sample_type_col = sample_type_col)
                    y_label <- "Hill number (D_q)"
                } else {
                    long <- prepare_tsallis_long(se, assay_name = assay_name, sample_type_col = sample_type_col)
                    y_label <- if (y == "S") "Tsallis entropy (S_q)" else "Hill number (D_q)"
                }
        if (nrow(long) == 0) stop("No tsallis values found in SummarizedExperiment")
        # long$q may be a factor; convert to numeric for plotting
        long$qnum <- as.numeric(as.character(long$q))
        stats_df <- dplyr::summarise(dplyr::group_by(long, group, qnum),
            median = median(tsallis, na.rm = TRUE),
            IQR = stats::IQR(tsallis, na.rm = TRUE), .groups = "drop"
        )
        # If the user requested Hill numbers but the long-format contains Tsallis S_q,
        # attempt conversion when safe (requires un-normalized S_q). Otherwise, advise the user to compute `what = 'D'`.
        if (y == "D" && !("hill" %in% SummarizedExperiment::assayNames(se))) {
            meta <- SummarizedExperiment::metadata(se)
            norm_flag <- if (!is.null(meta$norm)) meta$norm else FALSE
            if (isTRUE(norm_flag)) {
                stop("Cannot convert normalized Tsallis S_q to Hill numbers. Recompute with calculate_diversity(..., what = 'D', norm = FALSE) or supply a SummarizedExperiment containing a 'hill' assay.")
            }
            # perform element-wise conversion: sum_pq = 1 - (q - 1) * S_q ; D_q = sum_pq^(1/(1 - q)) ; handle q==1 by D = exp(-S)
            stats_df$D_median <- NA_real_
            stats_df$D_low <- NA_real_
            stats_df$D_high <- NA_real_
            for (i in seq_len(nrow(stats_df))) {
                qv <- stats_df$qnum[i]
                Sq <- stats_df$median[i]
                IQRv <- stats_df$IQR[i]
                if (is.na(qv) || is.na(Sq)) next
                if (abs(qv - 1) < .Machine$double.eps^0.5) {
                    # q ~ 1 : Shannon entropy H = -sum p log p ; D = exp(H) ; here S approximates H
                    stats_df$D_median[i] <- exp(-Sq)
                    stats_df$D_low[i] <- exp(-(Sq + IQRv / 2))
                    stats_df$D_high[i] <- exp(-(Sq - IQRv / 2))
                } else {
                    sum_pq <- 1 - (qv - 1) * Sq
                    if (sum_pq <= 0) {
                        stats_df$D_median[i] <- NA_real_
                        stats_df$D_low[i] <- NA_real_
                        stats_df$D_high[i] <- NA_real_
                    } else {
                        stats_df$D_median[i] <- sum_pq^(1 / (1 - qv))
                        # approximate extremes using median +/- IQR/2 on S
                        sum_pq_low <- 1 - (qv - 1) * (Sq + IQRv / 2)
                        sum_pq_high <- 1 - (qv - 1) * (Sq - IQRv / 2)
                        stats_df$D_low[i] <- ifelse(sum_pq_low > 0, sum_pq_low^(1 / (1 - qv)), NA_real_)
                        stats_df$D_high[i] <- ifelse(sum_pq_high > 0, sum_pq_high^(1 / (1 - qv)), NA_real_)
                    }
                }
            }
            # use D_median and ribbon from D_low/D_high
            p <- ggplot2::ggplot(
                stats_df,
                ggplot2::aes(x = qnum, y = D_median, color = group, fill = group)
            ) +
                ggplot2::geom_line(linewidth = 1) +
                ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = D_low, ymax = D_high),
                    alpha = 0.2,
                    color = NA
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
                ggplot2::labs(
                    title = "Median +- IQR of Tsallis entropy by group",
                    x = "q value",
                    y = y_label,
                    color = "Group",
                    fill = "Group"
                )
        } else {
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
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
                ggplot2::labs(
                    title = "Median +- IQR of Tsallis entropy by group",
                    x = "q value",
                    y = y_label,
                    color = "Group",
                    fill = "Group"
                )
        }
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
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    long <- prepare_tsallis_long(se, assay_name = assay_name)

    p <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = q, y = tsallis, fill = group)
    ) +
        ggplot2::geom_violin(
            alpha = 0.5,
            width = 0.9,
            position = ggplot2::position_dodge(width = 0.8)
        ) +
        ggplot2::geom_boxplot(
            width = 0.15,
            position = ggplot2::position_dodge(width = 0.8),
            outlier.shape = NA
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Tsallis entropy by q and group",
            x = "q value",
            y = "Tsallis entropy",
            fill = "Group"
        )

    # use default discrete ggplot2 colours (not viridis)
    p <- p + ggplot2::scale_fill_discrete(name = "Group")

    return(p)
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
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

    long <- prepare_tsallis_long(se, assay_name = assay_name)

    p <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = tsallis, color = group, fill = group)
    ) +
        ggplot2::geom_density(alpha = 0.3) +
        ggplot2::facet_wrap(~q, scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Tsallis entropy distributions by q and group",
            x = "Tsallis entropy",
            y = "Density",
            color = "Group",
            fill = "Group"
        )

    # use default discrete ggplot2 colours (not viridis)
    p <- p + ggplot2::scale_color_discrete(name = "Group") +
        ggplot2::scale_fill_discrete(name = "Group")

    return(p)
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
plot_volcano <- function(
  diff_df,
  x_col = "mean_difference",
  padj_col = "adjusted_p_values",
  label_thresh = 0.1,
  padj_thresh = 0.05,
  top_n = 5
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    df <- as.data.frame(diff_df)

    # Auto-detect a suitable x column if the default `mean_difference` is
    # not present. Prefer `mean_difference`, then `median_difference`, then
    # any column ending with `_difference`.
    if (!(x_col %in% colnames(df))) {
        if ("mean_difference" %in% colnames(df)) {
            x_col <- "mean_difference"
        } else if ("median_difference" %in% colnames(df)) {
            x_col <- "median_difference"
        } else {
            diffs <- grep("_difference$", colnames(df), value = TRUE)
            if (length(diffs) >= 1) x_col <- diffs[1]
        }
    }

    if (!(x_col %in% colnames(df)) || !(padj_col %in% colnames(df))) {
        stop(
            "Required columns not found in diff_df"
        )
    }

    df$padj_num <- as.numeric(df[[padj_col]])
    df$xval <- as.numeric(df[[x_col]])
    # Determine metric label for title/x-axis based on selected x_col
    metric_label_x <- if (grepl("median", x_col, ignore.case = TRUE)) {
        "Median"
    } else if (grepl("mean", x_col, ignore.case = TRUE)) {
        "Mean"
    } else {
        "Value"
    }

    # sanitize padj: NA -> 1, zeros -> tiny positive to avoid -Inf on log scale
    padj_clean <- df$padj_num
    padj_clean[is.na(padj_clean)] <- 1
    padj_clean[padj_clean <= 0] <- .Machine$double.xmin
    df$padj_clean <- padj_clean

    # prefer user-supplied `label` column; otherwise compute using thresholds
    if ("label" %in% colnames(df)) {
        df$label_flag <- as.character(df$label)
    } else {
        df$label_flag <- ifelse(
            abs(df$xval) >= label_thresh & df$padj_clean < padj_thresh,
            "significant",
            "non-significant"
        )
    }

    # drop rows with non-finite x or y values
    finite_idx <- is.finite(df$xval) & is.finite(df$padj_clean)
    df_plot <- df[finite_idx, , drop = FALSE]

    p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = xval, y = -log10(padj_clean), color = label_flag)
    ) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = c("non-significant" = "grey", "significant" = "red"),
            guide = "none"
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(padj_thresh),
            color = "red",
            linetype = "dashed"
        ) +
        ggplot2::geom_vline(
            xintercept = c(
                label_thresh,
                -label_thresh
            ),
            linetype = "dashed"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = sprintf("Volcano: %s difference vs adjusted p-value",
                metric_label_x),
            x = gsub("_", " ", tools::toTitleCase(x_col)),
            y = sprintf("-Log10(%s)", padj_col)
        )

    # annotate top N genes by smallest padj
    genes_col <- if ("genes" %in% colnames(df)) {
        "genes"
    } else if (!is.null(rownames(df))) {
        "rownames"
    } else {
        NULL
    }
    if (!is.null(genes_col) && top_n > 0) {
        if (genes_col == "rownames") df$genes <- rownames(df)
        ord <- order(df$padj_num, na.last = TRUE)
        top_idx <- head(ord[!is.na(df$padj_num[ord])], top_n)
        ann_df <- df[top_idx, , drop = FALSE]
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
                data = ann_df,
                ggplot2::aes(label = genes),
                size = 3
            )
        } else {
            p <- p + ggplot2::geom_text(
                data = ann_df,
                ggplot2::aes(label = genes),
                vjust = -0.5,
                size = 3
            )
        }
    }

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
    # Title label: use exact phrasing requested by user. Note: use
    # user's preferred spelling 'varience' for variance.
    agg_label <- switch(metric_choice,
                        median = "Metric: median",
                        mean = "Metric: mean",
                        variance = "Metric: varience")
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
                plot.title = ggplot2::element_text(size = 10, hjust = 0.53),
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
        # compute global fill limits across all genes so color scale is comparable
        mins <- c()
        maxs <- c()
        for (g in gene) {
            txs <- mapping$Transcript[mapping$Gen == g]
            txs <- intersect(txs, rownames(counts))
            if (length(txs) == 0) next
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
            df_summary <- aggregate(expr ~ tx + group,
                data = df_long,
                FUN = agg_fun
            )
            df_summary$log2expr <- log2(df_summary$expr + pseudocount)
            mins <- c(mins, min(df_summary$log2expr, na.rm = TRUE))
            maxs <- c(maxs, max(df_summary$log2expr, na.rm = TRUE))
        }
        if (length(mins) == 0) stop("No transcripts found for provided genes")
        fill_limits <- c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))

        plots <- lapply(seq_along(gene), function(i) {
            gname <- gene[i]
            pp <- make_plot_for_gene(gname,
                fill_limits = fill_limits
            )
            # set per-gene title while keeping a combined title below
            pp <- pp + ggplot2::labs(title = gname) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
            pp
        })

        # try patchwork first (collect guides), otherwise cowplot with shared legend
        if (requireNamespace("patchwork", quietly = TRUE)) {
            combined <- Reduce(`+`, plots) +
                patchwork::plot_layout(nrow = 1, guides = "collect") &
                ggplot2::theme(legend.position = "bottom")
            # add a single centered title
            combined <- combined + patchwork::plot_annotation(title = agg_label_unique,
                theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.53, size = 16))
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
            title_grob <- cowplot::ggdraw() + cowplot::draw_label(agg_label_unique, fontface = 'bold', x = 0.53, hjust = 0.5, size = 14)
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
            # layout: title row, plots row, legend row
            if (!is.null(output_file)) {
                png(filename = output_file, width = 800 * n, height = 480, res = 150)
                grid::grid.newpage()
                grid::pushViewport(
                    grid::viewport(
                        layout = grid::grid.layout(
                            3, n,
                            heights = grid::unit.c(
                                grid::unit(0.6, "cm"),
                                grid::unit(1, "null"),
                                grid::unit(0.7, "cm")
                            )
                        )
                    )
                )
                # draw title centered across columns
                vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(n))
                grid::pushViewport(vp_title)
                grid::grid.text(agg_label_unique, x = 0.53, gp = grid::gpar(fontface = "bold", fontsize = 14))
                grid::upViewport()
                for (i in seq_along(grobs)) {
                    vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = i)
                    grid::pushViewport(vp)
                    grid::grid.draw(grobs[[i]])
                    grid::upViewport()
                }
                if (!is.null(legend_grob)) {
                    vp_leg <- grid::viewport(
                        layout.pos.row = 3,
                        layout.pos.col = seq_len(n)
                    )
                    grid::pushViewport(vp_leg)
                    grid::grid.draw(legend_grob)
                    grid::upViewport()
                }
                grid::upViewport()
                dev.off()
                return(invisible(NULL))
            } else {
                grid::grid.newpage()
                grid::pushViewport(
                    grid::viewport(
                        layout = grid::grid.layout(
                            3, n,
                            heights = grid::unit.c(
                                grid::unit(0.6, "cm"),
                                grid::unit(1, "null"),
                                grid::unit(0.7, "cm")
                            )
                        )
                    )
                )
                vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(n))
                grid::pushViewport(vp_title)
                grid::grid.text(agg_label_unique, x = 0.53, gp = grid::gpar(fontface = "bold", fontsize = 14))
                grid::upViewport()
                for (i in seq_along(grobs)) {
                    vp <- grid::viewport(layout.pos.row = 2, layout.pos.col = i)
                    grid::pushViewport(vp)
                    grid::grid.draw(grobs[[i]])
                    grid::upViewport()
                }
                if (!is.null(legend_grob)) {
                    vp_leg <- grid::viewport(
                        layout.pos.row = 3,
                        layout.pos.col = seq_len(n)
                    )
                    grid::pushViewport(vp_leg)
                    grid::grid.draw(legend_grob)
                    grid::upViewport()
                }
                grid::upViewport()
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