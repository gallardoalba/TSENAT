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

# Small utility helpers to reduce repeated code paths when extracting
# samples, readcounts and tx->gene mappings from a
# SummarizedExperiment. Keeping these as focused helpers improves
# readability of the longer plotting functions below.
infer_samples_from_se <- function(se, samples = NULL, sample_type_col = "sample_type") {
    if (!is.null(samples)) return(as.character(samples))
    cd <- NULL
    try(cd <- SummarizedExperiment::colData(se), silent = TRUE)
    if (is.null(cd)) return(NULL)

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
        if (nm %in% colnames(cd)) return(as.character(cd[[nm]]))
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
    return(as.matrix(SummarizedExperiment::assay(se)))
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
        return(list(type = "vector",
                    mapping = as.character(txmap[[gene_col]][match(rownames(readcounts_mat), txmap[[tx_col]])])))
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
    if (!is.null(control) && control %in% uniq) return(control)
    if ("Normal" %in% uniq) return("Normal")
    # fallback to first level and message
    chosen <- uniq[1]
    message(sprintf("`control` not found; using '%s' instead", chosen))
    chosen
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
#' @examples
#' mat <- matrix(runif(8), nrow = 2, dimnames = list(c("g1","g2"), c("s1_q=0.1","s1_q=1","s2_q=0.1","s2_q=1")))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
#' plot_tsallis_gene_profile(se, gene = "g1")
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

    p <- ggplot2::ggplot() + ggplot2::theme_minimal(base_size = 14)

    if (isTRUE(show_samples)) {
        p <- p + ggplot2::geom_line(data = long_g, ggplot2::aes(x = qnum, y = tsallis, group = sample, color = group), alpha = 0.25)
    }

    p <- p +
        ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum, ymin = median - IQR/2, ymax = median + IQR/2, fill = group), alpha = 0.2, inherit.aes = FALSE) +
        ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = median, color = group), linewidth = 1.3) +
        ggplot2::labs(title = paste0(sel, ": Tsallis entropy q-curve profile"), x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
        ggplot2::scale_color_discrete(name = "Group") + ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, 
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
    if (!("genes" %in% colnames(df))) {
        if (!is.null(rownames(df))) df$genes <- rownames(df)
    }

    # If an external fc_df is provided, merge fold values into df
    if (!is.null(fc_df)) {
        fdf <- as.data.frame(fc_df, stringsAsFactors = FALSE)
        if (!("genes" %in% colnames(fdf)) && !is.null(rownames(fdf))) fdf$genes <- rownames(fdf)
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
            drop = FALSE], na.rm = TRUE)
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
    padj_candidates <- c("adjusted_p_values", "adj_p_value", "adj_p", "padj", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) padj_col[1] else NULL

    padj <- if (!is.null(padj_col)) as.numeric(df[[padj_col]]) else rep(1, length(yvals))
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < sig_alpha, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag, stringsAsFactors = FALSE)

    # format axis labels: remove underscores, collapse spaces, only first letter capitalized
    format_label <- function(lbl) {
        if (is.null(lbl)) return(NULL)
        s <- gsub("_", " ", lbl)
        s <- gsub("\\s+", " ", s)
        s <- trimws(s)
        s <- tolower(s)
        if (nchar(s) == 0) return(s)
        if (nchar(s) == 1) return(toupper(s))
        paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
    }

    x_label_formatted <- format_label(x_label)
    # prepare y label (use provided or fallback to fold_col) and replace 'log2' token with 'log10'
    y_label_raw <- y_label %||% fold_col
    y_label_formatted <- format_label(y_label_raw)
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
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain"),
                       axis.title = ggplot2::element_text(size = 14),
                       axis.text = ggplot2::element_text(size = 12))

    return(p)
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
        return(
            .plot_ma_core(
                x,
                fc_df = fc_res,
                sig_alpha = sig_alpha,
                x_label = x_label,
                y_label = y_label,
                title = title,
                ...
            )
        )
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
                ggplot2::geom_line(linewidth = 1.3) +
                ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = median - IQR / 2, ymax = median + IQR / 2),
                    alpha = 0.2,
                    color = NA
                ) +
                ggplot2::theme_minimal(base_size = 14) +
                ggplot2::labs(
                    title = "Global Tsallis q-curve: entropy across all genes",
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
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::scale_fill_discrete(name = "Group") +
        ggplot2::labs(
            title = "Violin plot: Tsallis entropy distribution across multiple q values",
            x = "q value",
            y = "Tsallis entropy",
            fill = "Group"
        ) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain",
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
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain",
                                                          margin = ggplot2::margin(b = 10)))
}


#' Volcano plot for differential results
#'
#' Create a volcano plot showing fold-change (x-axis) versus adjusted
#' p-value significance (y-axis). The function auto-detects a suitable x-axis
#' column if one is not provided and expects an adjusted p-value column for
#' significance coloring.
#'
#' @param diff_df Data.frame from `test_differential()` or similar results.
#'   Should contain p-values and optionally fold-change columns.
#' @param x_col Optional column name for the x-axis. If `NULL`, the function
#'   will try to auto-detect a suitable numeric column (excluding p-values).
#' @param padj_col Adjusted p-value column name (default: "adjusted_p_values").
#' @param label_thresh Fold-change threshold used to annotate points (default: 0.1).
#' @param padj_thresh Adjusted p-value cutoff for significance (default: 0.05).
#' @param top_n Number of top significant genes to label (default: 5).
#' @param title Optional plot title; if `NULL` a default title is used.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' df <- data.frame(
#'   gene = paste0("g", seq_len(10)),
#'   mean_difference = runif(10),
#'   adjusted_p_values = runif(10)
#' )
#' # plot_volcano(df, x_col = "mean_difference", padj_col = "adjusted_p_values")
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

    # Generate default title if not provided
    title_use <- title %||% "Volcano plot — Fold-change vs significance"

    # Format x-axis label: remove underscores and capitalize first letter
    x_label_formatted <- gsub("_", " ", x_col)
    x_label_formatted <- paste0(toupper(substr(x_label_formatted, 1, 1)), 
                                 substr(x_label_formatted, 2, nchar(x_label_formatted)))

    # Create base plot with title and theme
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
            yintercept = -log10(padj_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-label_thresh, label_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
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
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain", margin = ggplot2::margin(b = 4)),
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12)
        )

    return(p)
}

#' Plot top transcripts for a gene
#'
#' For a given gene, find transcripts using a tx->gene mapping, compute
#' per-transcript statistics between two sample groups, select the top N
#' transcripts by p-value, and plot their expression across groups.
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
#'   per group when plotting. One of c("median", "mean", "variance", "iqr").
#'   Use "iqr" to compute the interquartile range. Defaults to "median".
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
        "log2expr",
        # variables used in do_plot_ma_core
        "x",
        "y",
        "padj",
        "significant"
    ))
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

## Internal helpers for `plot_top_transcripts` refactor
## Create per-gene plot and combine multiple gene plots into final output
.ptt_make_plot_for_gene <- function(gene_single, mapping, counts, samples, top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL) {
    require_pkgs(c("ggplot2", "tidyr"))

    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))
    if (length(txs) == 0) stop("No transcripts found for gene: ", gene_single)
    if (!is.null(top_n)) txs <- head(txs, top_n)

    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)
    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))

    df_summary <- aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))

    p <- ggplot2::ggplot(
        df_summary,
        ggplot2::aes(x = group, y = tx, fill = log2expr)
    ) +
        ggplot2::geom_tile(color = "white", width = 0.95, height = 0.95) +
        ggplot2::scale_fill_viridis_c(
            option = "viridis",
            direction = -1,
            na.value = "grey80",
            limits = fill_limits
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(
            title = agg_label_unique,
            x = NULL,
            y = NULL,
            fill = "log2(expr)"
        ) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 12),
            axis.text.x = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(size = 16, hjust = 0.6),
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

.ptt_combine_plots <- function(plots, output_file = NULL, agg_label_unique) {
    require_pkgs(c("ggplot2"))
    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- Reduce(`+`, plots) +
            patchwork::plot_layout(nrow = 1, guides = "collect") &
            ggplot2::theme(legend.position = "bottom")
        combined <- combined + patchwork::plot_annotation(title = agg_label_unique,
            theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.6, size = 16,
                                                                        margin = ggplot2::margin(b = 10)))
        )
        return(combined)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
        legend <- cowplot::get_legend(p_for_legend)
        plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
        grid <- cowplot::plot_grid(plotlist = plots_nolegend, nrow = 1, align = "hv")
        title_grob <- cowplot::ggdraw() + cowplot::draw_label(agg_label_unique, fontface = 'plain', x = 0.6, hjust = 0.5, size = 16)
        result_plot <- cowplot::plot_grid(title_grob, grid, legend, ncol = 1, rel_heights = c(0.08, 1, 0.08))
        if (!is.null(output_file)) {
            ggplot2::ggsave(output_file, result_plot)
            return(invisible(NULL))
        }
        return(result_plot)
    } else {
        plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
        grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)
        g_full <- ggplot2::ggplotGrob(plots[[1]])
        legend_idx <- which(vapply(g_full$grobs, function(x) x$name, character(1)) == "guide-box")
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
}

## Prepare and validate inputs for `plot_top_transcripts`
.ptt_prepare_inputs <- function(counts, readcounts, samples, coldata, sample_type_col, tx2gene, res, top_n, pseudocount, output_file, metric = c("median", "mean", "variance", "iqr")) {
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
    } else stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")

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
 #' @examples
 #' tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
 #' rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
 #' colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
 #' tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 2), stringsAsFactors = FALSE)
 #' samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
 #' plot_top_transcripts(tx_counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2gene, top_n = 2)
 #' 
 
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
    # If `gene` is not provided, select top genes from `res` using `top_n`.
    if (is.null(gene)) {
        if (is.null(res)) stop("Either 'gene' or 'res' must be provided")
        if (!("genes" %in% colnames(res))) stop("Provided 'res' must contain a 'genes' column")
        if ("adjusted_p_values" %in% colnames(res)) {
            ord <- order(res$adjusted_p_values, na.last = NA)
        } else if ("raw_p_values" %in% colnames(res)) {
            ord <- order(res$raw_p_values, na.last = NA)
        } else {
            ord <- seq_len(nrow(res))
        }
        genes_sel <- as.character(res$genes[ord])
        genes_sel <- unique(genes_sel)
        gene <- head(genes_sel, top_n)
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
        return(result_plot)
    }
}

