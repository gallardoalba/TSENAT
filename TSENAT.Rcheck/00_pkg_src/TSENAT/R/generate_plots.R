## Helpers for generating diversity plots (density, violin, MA)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Gene", "diversity", "sample_type", "fold", "significant", "value", ".", "group", "tsallis", "q", "median", "IQR"))

# Internal helper: prepare long-format tsallis data from a SummarizedExperiment
prepare_tsallis_long <- function(se, assay_name = "diversity", sample_type_col = "sample_type") {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")

  mat <- SummarizedExperiment::assay(se, assay_name)
  if (is.null(mat)) stop("Assay not found: ", assay_name)
  df <- as.data.frame(mat)
  genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) SummarizedExperiment::rowData(se)$genes else rownames(df)
  df <- cbind(df, Gene = genes_col)

  long <- tidyr::pivot_longer(df, -Gene, names_to = "sample_q", values_to = "tsallis")
  if (any(grepl("_q=", long$sample_q))) {
    long <- tidyr::separate(long, sample_q, into = c("sample", "q"), sep = "_q=")
    long$q <- as.factor(as.numeric(long$q))
  } else {
    long$sample <- long$sample_q
    long$q <- NA
  }

  if (!is.null(sample_type_col) && sample_type_col %in% colnames(SummarizedExperiment::colData(se))) {
    st_vec <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
    names(st_vec) <- sub("_q=.*", "", colnames(mat))
    st_map <- st_vec[!duplicated(names(st_vec))]
    long$group <- unname(st_map[as.character(long$sample)])
  } else {
    sample_base <- as.character(long$sample)
    long$group <- infer_sample_group(sample_base)
    missing_idx <- which(is.na(long$group))
    if (length(missing_idx) > 0) {
      long$group[missing_idx] <- vapply(sample_base[missing_idx], function(s) {
        if (grepl("_", s)) sub(".*_", "", s) else s
      }, character(1))
    }
  }

  long <- long[!is.na(long$tsallis), , drop = FALSE]
  return(as.data.frame(long))
}
}

#' Plot diversity distributions (density) by sample type
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Optional column name in `colData(se)` that contains sample types. If missing, sample type will be inferred from column names (suffix after last underscore) or classified as 'Group'.
#' @return A `ggplot` object with layered density plots.
#' @importFrom ggplot2 ggplot aes geom_density facet_grid scale_color_manual guides theme_minimal labs
#' @export
plot_diversity_density <- function(se, assay_name = "diversity", sample_type_col = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

  mat <- SummarizedExperiment::assay(se, assay_name)
  if (is.null(mat)) stop("Assay not found: ", assay_name)

  df <- as.data.frame(mat)
  genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) SummarizedExperiment::rowData(se)$genes else rownames(df)
  df <- cbind(df, Gene = genes_col)

  long <- tidyr::pivot_longer(df, -Gene, names_to = "sample", values_to = "diversity")

  if (!is.null(sample_type_col) && sample_type_col %in% colnames(SummarizedExperiment::colData(se))) {
    st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
    names(st) <- colnames(mat)
    long$sample_type <- unname(st[long$sample])
  } else {
    # strip possible '_q=...' suffix produced by calculate_diversity()
    sample_base <- sub("_q=.*", "", long$sample)
    # attempt to infer sample type using helper (supports TCGA barcodes and _N/_T suffixes)
    long$sample_type <- infer_sample_group(sample_base)
    # if inference failed, fallback to suffix-after-last-underscore or sample name
    missing_idx <- which(is.na(long$sample_type))
    if (length(missing_idx) > 0) {
      long$sample_type[missing_idx] <- vapply(sample_base[missing_idx], function(s) {
        if (grepl("_", s)) sub(".*_", "", s) else s
      }, character(1))
    }
  }

  long <- long[!is.na(long$diversity), , drop = FALSE]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = diversity, group = sample, color = sample_type)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::facet_grid(. ~ sample_type) +
    ggplot2::scale_color_manual(values = c("black", "darkorchid4")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Diversity values", y = "Density")

  return(p)
}


#' Plot violin of per-gene mean diversity by sample type
#' @importFrom magrittr %>%
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param sample_type_col Optional column name in `colData(se)` containing sample types.
#' @return A `ggplot` violin plot object.
#' @export
plot_mean_violin <- function(se, assay_name = "diversity", sample_type_col = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

  mat <- SummarizedExperiment::assay(se, assay_name)
  df <- as.data.frame(mat)
  genes_col <- if (!is.null(SummarizedExperiment::rowData(se)$genes)) SummarizedExperiment::rowData(se)$genes else rownames(df)
  df <- cbind(df, Gene = genes_col)
  long <- tidyr::pivot_longer(df, -Gene, names_to = "sample", values_to = "diversity")

  if (!is.null(sample_type_col) && sample_type_col %in% colnames(SummarizedExperiment::colData(se))) {
    st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
    names(st) <- colnames(mat)
    long$sample_type <- unname(st[long$sample])
  } else {
    sample_base <- sub("_q=.*", "", long$sample)
    long$sample_type <- infer_sample_group(sample_base)
    missing_idx <- which(is.na(long$sample_type))
    if (length(missing_idx) > 0) {
      long$sample_type[missing_idx] <- vapply(sample_base[missing_idx], function(s) {
        if (grepl("_", s)) sub(".*_", "", s) else s
      }, character(1))
    }
  }

  long <- long[!is.na(long$diversity), , drop = FALSE]

  mean_df <- long %>%
    as.data.frame() %>%
    split(.$sample_type) %>%
    lapply(function(sub) {
      stats::aggregate(sub$diversity, by = list(sub$Gene), mean)
    })

  # combine and plot
  plot_df <- do.call(rbind, lapply(names(mean_df), function(nm) {
    d <- mean_df[[nm]]
    colnames(d) <- c("Gene", "value")
    d$sample_type <- nm
    d
  }))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample_type, y = value, fill = sample_type)) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Samples", y = "Diversity") +
    ggplot2::scale_fill_viridis_d(name = "Group")

  return(p)
}


#' Plot MA plot for difference results
#'
#' @param diff_df Data.frame returned by `calculate_difference` (or similar) containing mean columns and a `log2_fold_change` column, and `adjusted_p_values`.
#' @param mean_cols Optional character vector of length 2 with the names of the mean columns (defaults to first two columns that end with `_mean`).
#' @param fold_col Name of the fold-change column (default: `log2_fold_change`).
#' @param padj_col Name of the adjusted p-value column (default: `adjusted_p_values`).
#' @param sig_alpha Threshold for significance (default: 0.05).
#' @return A `ggplot` MA-plot object.
#' @export
plot_ma <- function(diff_df, mean_cols = NULL, fold_col = "log2_fold_change", padj_col = "adjusted_p_values", sig_alpha = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  if (is.null(mean_cols)) {
    mean_cols <- grep("_mean$", colnames(diff_df), value = TRUE)
    if (length(mean_cols) < 2) stop("Could not find two mean columns")
    mean_cols <- mean_cols[1:2]
  }

  df <- diff_df
  df$mean <- apply(df[, mean_cols], 1, function(x) mean(as.numeric(x), na.rm = TRUE))
  df$padj <- if (padj_col %in% colnames(df)) df[[padj_col]] else NA
  df$fold <- if (fold_col %in% colnames(df)) df[[fold_col]] else NA
  df$significant <- ifelse(!is.na(df$padj) & df$padj < sig_alpha, "significant", "non-significant")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean, y = fold, color = significant)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("non-significant" = "black", "significant" = "red")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Mean diversity", y = "Log2 fold change")

  return(p)
}


#' Plot median +- IQR of Tsallis entropy across q values by group
#'
#' This reproduces the `tsallis-q-curve-mean-sd` plot from the vignette:
#' for each q value, compute per-gene Tsallis entropy per sample, then
#' summarize across genes by group (median and IQR) and plot median with a
#' ribbon spanning median +- IQR/2.
#'
#' @param readcounts Numeric matrix or data.frame with transcripts as rows and samples as columns.
#' @param genes Character vector assigning a gene id to each row of `readcounts`.
#' @param q_values Numeric vector of q values to evaluate (default `seq(0.01,2,by=0.01)`).
#' @param group_pattern Regular expression used to detect the first group in sample names (default `"_N$"`).
#' @param group_names Character vector of length 2 with names for groups (default `c("Normal","Tumor")`).
#' @return A `ggplot` object showing median +- IQR across q values by group.
#' @export
plot_tsallis_q_curve <- function(readcounts, genes, q_values = seq(0.01, 2, by = 0.01),
                                 group_pattern = "_N$", group_names = c("Normal", "Tumor")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

  if (is.data.frame(readcounts)) readcounts <- as.matrix(readcounts)
  if (!is.matrix(readcounts)) stop("`readcounts` must be a matrix or data.frame")
  if (length(genes) != nrow(readcounts)) stop("Length of `genes` must match nrow(readcounts)")

  get_group_stats_iqr <- function(q) {
    # compute per-gene Tsallis per sample for this q
    tsallis_df <- as.data.frame(sapply(seq_len(ncol(readcounts)), function(ci) {
      tapply(readcounts[, ci], genes, calculate_tsallis_entropy, q = q, norm = TRUE)
    }, simplify = TRUE))
    colnames(tsallis_df) <- colnames(readcounts)
    tsallis_df$Gene <- rownames(tsallis_df)

    long <- tidyr::pivot_longer(tsallis_df, -Gene, names_to = "sample", values_to = "tsallis")
    long$group <- ifelse(grepl(group_pattern, long$sample), group_names[1], group_names[2])

    stats <- dplyr::summarise(dplyr::group_by(long, group),
                              median = median(tsallis, na.rm = TRUE),
                              IQR = stats::IQR(tsallis, na.rm = TRUE), .groups = 'drop')
    stats <- dplyr::mutate(stats, q = q)
    return(stats)
  }

  stats_df <- do.call(rbind, lapply(q_values, get_group_stats_iqr))

  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = q, y = median, color = group, fill = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = median - IQR/2, ymax = median + IQR/2), alpha = 0.2, color = NA) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = "Median +- IQR of Tsallis entropy by group",
                  x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group")

  return(p)
}


##' Violin plot of Tsallis entropy for multiple q values
##'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with multiple q values (colnames like 'Sample_q=0.01').
#' @param assay_name Name of the assay to use (default: "diversity").
#' @return A `ggplot` violin plot object faceted/colored by group and q.
#' @export
plot_tsallis_violin_multq <- function(se, assay_name = "diversity") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

  long <- prepare_tsallis_long(se, assay_name = assay_name)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = q, y = tsallis, fill = group)) +
    ggplot2::geom_violin(alpha = 0.5, width = 0.9, position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_boxplot(width = 0.15, position = ggplot2::position_dodge(width = 0.8), outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Tsallis entropy by q and group", x = "q value", y = "Tsallis entropy", fill = "Group")

  return(p)
}


##' Density plot of Tsallis entropy for multiple q values
##'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with multiple q values (colnames like 'Sample_q=0.01').
#' @param assay_name Name of the assay to use (default: "diversity").
#' @return A `ggplot` density plot object faceted by q and colored by group.
#' @export
plot_tsallis_density_multq <- function(se, assay_name = "diversity") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

  long <- prepare_tsallis_long(se, assay_name = assay_name)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = tsallis, color = group, fill = group)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::facet_wrap(~q, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Tsallis entropy distributions by q and group", x = "Tsallis entropy", y = "Density", color = "Group", fill = "Group")

  return(p)
}


#' Volcano plot for differential results
#'
#' @param diff_df Data.frame from `calculate_difference()` containing at least
#'   `mean_difference` and an adjusted p-value column (default `adjusted_p_values`).
#' @param x_col Column name for x-axis (default `mean_difference`).
#' @param padj_col Column name for adjusted p-values (default `adjusted_p_values`).
#' @param label_thresh Absolute x threshold to mark significance (default 0.1).
#' @param padj_thresh Adjusted p-value cutoff (default 0.05).
#' @return ggplot volcano plot.
#' @export
plot_volcano <- function(diff_df, x_col = "mean_difference", padj_col = "adjusted_p_values", label_thresh = 0.1, padj_thresh = 0.05, top_n = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  df <- as.data.frame(diff_df)
  if (!(x_col %in% colnames(df)) || !(padj_col %in% colnames(df))) stop("Required columns not found in diff_df")

  df$padj_num <- as.numeric(df[[padj_col]])
  df$xval <- as.numeric(df[[x_col]])
  df$label_flag <- ifelse(abs(df$xval) >= label_thresh & df$padj_num < padj_thresh, "significant", "non-significant")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = xval, y = -log10(padj_num), color = label_flag)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("non-significant" = "grey", "significant" = "red"), guide = "none") +
    ggplot2::geom_hline(yintercept = -log10(padj_thresh), color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(label_thresh, -label_thresh), linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Volcano: mean difference vs adjusted p-value", x = x_col, y = sprintf("-Log10(%s)", padj_col))

  # annotate top N genes by smallest padj
  genes_col <- if ("genes" %in% colnames(df)) "genes" else if (!is.null(rownames(df))) "rownames" else NULL
  if (!is.null(genes_col) && top_n > 0) {
    if (genes_col == "rownames") df$genes <- rownames(df)
    ord <- order(df$padj_num, na.last = TRUE)
    top_idx <- head(ord[!is.na(df$padj_num[ord])], top_n)
    ann_df <- df[top_idx, , drop = FALSE]
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = ann_df, ggplot2::aes(label = genes), size = 3)
    } else {
      p <- p + ggplot2::geom_text(data = ann_df, ggplot2::aes(label = genes), vjust = -0.5, size = 3)
    }
  }

  return(p)
}
