# # Helpers for generating diversity plots (density, violin, MA)

#' @importFrom utils head
#' @importFrom grDevices png dev.off
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Gene", "diversity", "sample", "sample_q", "sample_type",
        "fold", "significant", "value", ".", "group", "tsallis", "q", "median", "IQR",
        "padj_num", "padj_clean", "xval", "label_flag", "genes", "mean", "tx", "expr",
        "tx_cond", "log2expr", "x", "y", "padj"))

    require_pkgs <- function(pkgs) {
        missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing)) {
            stop(sprintf("%s required", paste(missing, collapse = ", ")))
        }
        invisible(TRUE)
    }
}

# Small utility helpers to reduce repeated code paths when extracting samples,
# readcounts and tx->gene mappings from a SummarizedExperiment. Keeping these
# as focused helpers improves readability of the longer plotting functions
# below.

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
    require_pkgs(c("ggplot2"))

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
    fold_col_candidates <- c("log2_fold_change", "logFC", "fold", "estimate_interaction", "fold_change")
    fold_col <- intersect(fold_col_candidates, colnames(df))
    if (length(fold_col) == 0) stop("Could not find a fold-change column in input")
    fold_col <- fold_col[1]

    # Detect p-value column for significance flagging
    padj_candidates <- c("padj", "adjusted_p_values", "adj_p_value", "adj_p", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) padj_col[1] else NULL
    padj <- if (!is.null(padj_col)) as.numeric(df[[padj_col]]) else rep(1, nrow(df))
    padj[is.na(padj)] <- 1

    # Validate mean/median column consistency
    mean_cols <- grep("_mean$", colnames(df), ignore.case = TRUE, value = TRUE)
    median_cols <- grep("_median$", colnames(df), ignore.case = TRUE, value = TRUE)
    
    if (length(mean_cols) > 0 && length(median_cols) > 0) {
        stop("Could not find two mean or two median columns - found both mean and median columns. ",
            "Ensure input contains either mean columns (e.g., A_mean, B_mean) OR median columns (e.g., A_median, B_median), not both.")
    }
    
    if (length(mean_cols) > 0 && length(mean_cols) < 2) {
        stop("Could not find two mean or two median columns - found ", length(mean_cols), " mean column(s). ",
            "Ensure input contains at least two mean columns (e.g., A_mean, B_mean).")
    }
    
    if (length(median_cols) > 0 && length(median_cols) < 2) {
        stop("Could not find two mean or two median columns - found ", length(median_cols), " median column(s). ",
            "Ensure input contains at least two median columns (e.g., A_median, B_median).")
    }
    
    # Determine which columns to use for mean calculation
    mean_cols_to_use <- if (length(mean_cols) > 0) mean_cols else if (length(median_cols) > 0) median_cols else NULL

    # Prepare MA plot data with label formatting
    prep <- .prepare_ma_plot_df(df, fold_col = fold_col, mean_cols = mean_cols_to_use, x_label = x_label,
        y_label = y_label)
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
        ggplot2::geom_point(alpha = 0.75, size = 3.2) +
        ggplot2::scale_color_manual(values = .significance_colors(), guide = "none") +
        ggplot2::labs(x = x_label_formatted, y = y_label_formatted)

    # Apply publication theme and settings
    p <- .apply_publication_theme(p, title = title %||% "MA plot: mean vs log10 fold-change",
        base_size = 11) +
        ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"))

    p
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

.plot_tsallis_violin_singleq <- function(se, assay_name = "diversity", title = NULL) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))

    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) >
        0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }

    # Fallback: use prepare_tsallis_long for data transformation
    long <- .prepare_tsallis_long(se, assay_name = assay_name)

    if (nrow(long) == 0)
        stop("No data found in the long format dataframe")

    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }

    # Set title
    title_use <- title %||% sprintf("Violin plot: Tsallis entropy at q = %g", q_val)

    # Create violin plot with publication theme
    p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = tsallis, fill = group)) +
        ggplot2::geom_violin(alpha = 0.5, width = 0.7, position = ggplot2::position_dodge(width = 0.8)) +
        ggplot2::geom_boxplot(width = 0.2, position = ggplot2::position_dodge(width = 0.8),
            outlier.shape = NA, alpha = 0.8)

    p <- .apply_group_aesthetics(p, palette = "palette_blue_red", legend_name = "Group",
        legend_position = "none")
    
    p <- .apply_publication_theme(p, title = title_use, base_size = 11) +
        ggplot2::labs(x = "Group", y = "Tsallis entropy")

    p
}


#' Density plot of Tsallis entropy for a single q value
#'
#' Creates a density plot showing the distribution of Tsallis entropy for a
#' specific q value,
#' with different groups (conditions) represented by different colors.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`
#' containing
#'   entropy values at one or more q values.
#' @param q_value The specific q value to plot (numeric, e.g., 1, 2, 0.5).
#' @param assay_name Name of the assay to use (default: 'diversity').
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a density plot colored by group.
#'
#' @noRd

.plot_tsallis_density_singleq <- function(se, assay_name = "diversity", title = NULL) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))

    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) >
        0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }

    # Fallback: use prepare_tsallis_long for data transformation
    long <- .prepare_tsallis_long(se, assay_name = assay_name)

    if (nrow(long) == 0)
        stop("No data found in the long format dataframe")

    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }

    # Set title
    title_use <- title %||% sprintf("Density plot: Tsallis entropy at q = %g", q_val)

    # Create density plot with publication theme
    p <- ggplot2::ggplot(long, ggplot2::aes(x = tsallis, color = group, fill = group)) +
        ggplot2::geom_density(alpha = 0.3, linewidth = 1)

    p <- .apply_group_aesthetics(p, palette = "palette_blue_red", legend_name = "Group")
    
    p <- .apply_publication_theme(p, title = title_use, base_size = 11) +
        ggplot2::labs(x = "Tsallis entropy", y = "Density")

    p
}


#' Combined Violin and Density Plot Grid for Single q Value
#'
#' Creates a side-by-side grid layout with a violin plot on the left and a
#' density plot
#' on the right, both showing Tsallis entropy distribution for the q value
#' in the provided
#' SummarizedExperiment (which should contain a single q value).
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity`
#' containing
#'   entropy values at a single q value.
#' @param assay_name Name of the assay to use (default: 'diversity').
#' @param title Optional base title. If NULL, auto-generated based on q value.
#' @param output_file Character or NULL. Optional file path to save the plot
#' as an image.
#'   If provided, the plot will be saved with appropriate dimensions.
#'   Default: NULL (no file output, only return object).
#'
#' @return A `ggplot2` object showing a 1x2 grid with violin plot on the
#' left and
#'   density plot on the right.
#'
#' @export
#' @examples
#' # Plot 8: Violin and density plots of Tsallis entropy distribution
#' data(readcounts)
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' 
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' analysis <- calculate_diversity_s4(analysis, q = 1.0)
#' p <- plot_tsallis_violin_density_grid_s4(analysis)
#' if (!is.null(p)) print(p)
#'
plot_tsallis_violin_density_grid_s4 <- function(se, assay_name = "diversity", title = NULL,
    output_file = NULL) {
    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Handle TSENATAnalysis objects - extract first diversity result
    if (methods::is(se, "TSENATAnalysis")) {
        if (length(se@diversity_results) == 0) {
            stop("No diversity results found in TSENATAnalysis object. Run calculate_diversity_s4() first.")
        }
        # Extract first diversity result
        se <- se@diversity_results[[1]]
    }

    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) >
        0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }

    # Fallback: use prepare_long_format for data transformation
    long <- .prepare_long_format(se, assay_name = assay_name)

    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }

    # Generate base title
    base_title <- title %||% sprintf("Tsallis entropy at q = %g", q_val)

    # Create individual plots
    p_violin <- .plot_tsallis_violin_singleq(se = se, assay_name = assay_name, title = "Violin")

    p_density <- .plot_tsallis_density_singleq(se = se, assay_name = assay_name,
        title = "Density")

    # Arrange plots side by side: violin on left, density on right
    grid <- cowplot::plot_grid(p_violin, p_density, nrow = 1, ncol = 2, align = "h",
        axis = "b")

    # Add overall title and subtitle above the grid
    title_grob <- .create_title_grob("Tsallis Entropy Distribution by Group",
        subtitle = "Violin and density plots across samples",
        title_size = 19, subtitle_size = 15)
    grid_with_title <- cowplot::plot_grid(title_grob, grid, nrow = 2, rel_heights = c(0.08, 1))

    # Save to file if output_file is provided
    if (!is.null(output_file)) {
        plot_dims <- .calculate_plot_dims(width_inches = 12, aspect_type = "standard", dpi_output = 100)
        ggplot2::ggsave(output_file, plot = grid_with_title, width = plot_dims$width, height = plot_dims$height,
            dpi = plot_dims$dpi, create.dir = TRUE)
    }

    return(grid_with_title)
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
        ggplot2::geom_point(alpha = 0.75, size = 3.4) +
        ggplot2::scale_color_manual(values = .significance_colors(), guide = "none")
    
    # Add reference lines using Phase 5 helper
    p <- .add_reference_lines(p,
        h_intercept = -log10(sig_alpha),
        v_intercept = c(-label_thresh, label_thresh),
        h_color = "gray50", v_color = "gray50")
    
    p <- .apply_publication_theme(p, title = title_use, base_size = 11) +
        ggplot2::labs(x = x_label_formatted, y = paste0("-Log10(", padj_label_formatted, ")"))

    p
}


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
#' # .plot_volcano_ma_grid(x, sig_alpha = 0.05)
#'

#' @noRd

.plot_volcano_ma_grid <- function(diff_df, x_col = NULL, padj_col = "padj", label_thresh = 0.1,
    sig_alpha = 0.05, top_n = 5, title_volcano = NULL, title_ma = "Tsallis-based MA plot",
    ...) {
    # Require cowplot for grid arrangement
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("cowplot package required for .plot_volcano_ma_grid()")
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


#' Internal helper to compute fill limits across multiple genes (not exported)
#' @noRd
.plot_transcript_fill_limits <- function(genes, mapping, counts, samples, top_n,
    agg_fun, pseudocount) {
    mins <- maxs <- c()
    for (g in genes) {
        txs <- mapping$Transcript[mapping$Gen == g]
        txs <- intersect(txs, rownames(counts))
        if (length(txs) == 0)
            next
        if (!is.null(top_n))
            txs <- head(txs, top_n)
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
    if (length(mins) == 0)
        stop("No transcripts found for provided genes")
    c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))
}

#' Internal helper to draw grid layout with title, plots, and legend using
#' base grid
#' @noRd
.plot_transcript_grid_draw <- function(grobs, title, legend_grob, ncol, heights,
    to_file = NULL) {
    # If no output file is provided and no graphics device is open, render to a
    # temporary pdf device so that plotting in non-interactive sessions does
    # not create `Rplots.pdf` in the working directory.
    temp_dev <- FALSE
    # Only open a temporary PDF device when: - caller did not supply an output
    # file (`to_file` is NULL), - the session is non-interactive, and - no
    # graphics device is currently open (dev.cur() == 1)
    if (is.null(to_file) && !interactive() && grDevices::dev.cur() == 1L) {
        tmp <- tempfile("TSENAT_plot_", fileext = ".pdf")
        grDevices::pdf(tmp)
        temp_dev <- TRUE
        # Ensure device is closed and temporary file removed on exit
        on.exit({
            try(grDevices::dev.off(), silent = TRUE)
            if (file.exists(tmp)) unlink(tmp)
        }, add = TRUE)
    }

    # Calculate number of rows needed for plots
    nrow_plots <- ceiling(length(grobs)/ncol)
    nrow_total <- 2 + nrow_plots  # title + plot rows + legend

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow_total, ncol,
        heights = heights)))
    # Title row
    vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(ncol))
    grid::pushViewport(vp_title)
    grid::grid.text("Transcript level expression", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 14,
        fontface = "bold"))
    grid::grid.text(paste0("Top genes with metric ", title), x = 0.5, y = 0.2, gp = grid::gpar(fontsize = 11,
        fontface = "italic", col = "gray40"))
    grid::upViewport()
    # Plot rows
    for (i in seq_along(grobs)) {
        plot_row_idx <- ((i - 1)%/%ncol) + 2
        plot_col_idx <- ((i - 1)%%ncol) + 1
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
    if (!is.null(to_file))
        grDevices::dev.off()
    invisible(NULL)
}

## Internal helpers for `plot_top_transcripts` refactor Create per-gene plot
## and combine multiple gene plots into final output

.make_plot_for_genemake_plot_for_gene <- function(gene_single, mapping, counts, samples,
    top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL, font_scale = 1) {
    require_pkgs(c("ggplot2", "tidyr"))
    built <- .make_plot_for_genebuild_tx_long(gene_single, mapping, counts, samples,
        NULL)
    df_summary <- .make_plot_for_geneaggregate_df_long(built$df_long, agg_fun, pseudocount)
    .make_plot_for_genebuild_plot_from_summary(df_summary, agg_label_unique, fill_limits,
        font_scale = font_scale)
}


.make_plot_for_genecombine_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
    require_pkgs(c("ggplot2"))
    # Allow callers to pass a single character second argument as the
    # `agg_label_unique` for convenience (legacy test call patterns).
    if (is.null(agg_label_unique) && !is.null(output_file) && is.character(output_file) &&
        length(output_file) == 1) {
        agg_label_unique <- output_file
        output_file <- NULL
    }
    if (requireNamespace("patchwork", quietly = TRUE)) {
        .make_plot_for_genecombine_patchwork(plots, agg_label_unique)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        .make_plot_for_genecombine_cowplot(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        .make_plot_for_genecombine_grid(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    }
}

## Prepare and validate inputs for `plot_top_transcripts`

.make_plot_for_geneprepare_inputs <- function(counts, readcounts = NULL, samples = NULL,
    coldata = NULL, condition_col = "condition", tx2gene = NULL, res = NULL, top_n = NULL,
    pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance",
        "iqr")) {
    # handle selecting genes from `res` is left to caller; this function
    # focuses on normalizing counts, samples and tx2gene mapping and preparing
    # agg functions
    if (inherits(counts, "SummarizedExperiment")) {
        require_pkgs(c("SummarizedExperiment", "S4Vectors"))
        se <- counts
        counts_mat <- .get_readcounts_from_se(se, readcounts)
        counts <- as.matrix(counts_mat)
        samples <- .infer_samples_from_se(se, samples, condition_col = condition_col)

        if (is.null(tx2gene)) {
            txres <- .get_tx2gene_from_se(se, counts)
            if (!is.null(txres) && !is.null(txres$mapping)) {
                mapping <- data.frame(Transcript = rownames(counts), Gen = as.character(txres$mapping),
                  stringsAsFactors = FALSE)
                tx2gene <- mapping
            }
        }
    }

    if (!is.matrix(counts) && !is.data.frame(counts))
        stop("`counts` must be a matrix or data.frame with transcripts as rownames")
    counts <- as.matrix(counts)
    if (is.null(rownames(counts)))
        stop("`counts` must have rownames corresponding to transcript identifiers")

    # derive samples from coldata if needed
    if (is.null(samples)) {
        if (!is.null(coldata)) {
            if (is.character(coldata) && length(coldata) == 1) {
                if (!file.exists(coldata))
                  stop("coldata file not found: ", coldata)
                cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
            } else if (is.data.frame(coldata)) {
                cdf <- coldata
            } else {
                stop("`coldata` must be a data.frame or path to a tab-delimited file")
            }

            if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
                samples <- as.character(cdf[colnames(counts), condition_col])
            } else {
                sample_id_cols <- c("sample", "Sample", "sample_id", "id")
                sid <- intersect(sample_id_cols, colnames(cdf))
                if (length(sid) > 0) {
                  sid <- sid[1]
                  if (!all(colnames(counts) %in% as.character(cdf[[sid]])))
                    stop("coldata sample id column does not match column names of counts")
                  row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
                  samples <- as.character(cdf[[condition_col]][row_ix])
                } else {
                  stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
                }
            }
        } else {
            stop("Either 'samples' or 'coldata' must be provided to determine sample groups")
        }
    }

    # normalize tx2gene mapping
    if (is.null(tx2gene))
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene))
            stop("tx2gene file not found: ", tx2gene)
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping)))
        stop("tx2gene must have columns 'Transcript' and 'Gen'")

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("ggplot2 required for plotting")

    if (!is.null(samples) && length(samples) != ncol(counts))
        stop("Length of `samples` must equal number of columns in `counts`")

    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) stats::var(x,
            na.rm = TRUE), iqr = function(x) stats::IQR(x, na.rm = TRUE))
    agg_label_metric <- if (metric_choice == "iqr")
        "IQR" else metric_choice
    agg_label <- agg_label_metric
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label

    list(counts = counts, samples = samples, mapping = mapping, metric_choice = metric_choice,
        agg_fun = agg_fun, agg_label_unique = agg_label_unique, top_n = top_n, pseudocount = pseudocount,
        output_file = output_file)
}




# Helpers for plot_top_transcripts internals


.make_plot_for_geneselect_genes_from_res <- function(res, top_n) {
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


.make_plot_for_geneinfer_samples_from_coldata <- function(coldata, counts, condition_col) {
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
        as.character(cdf[colnames(counts), condition_col])
    } else {
        sample_id_cols <- c("sample", "Sample", "sample_id", "id")
        sid <- intersect(sample_id_cols, colnames(cdf))
        if (length(sid) > 0) {
            sid <- sid[1]
            if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
                stop("coldata sample id column does not match column names of counts")
            }
            row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
            as.character(cdf[[condition_col]][row_ix])
        } else {
            stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
        }
    }
}


.make_plot_for_generead_tx2gene <- function(tx2gene) {
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


.make_plot_for_genemake_agg <- function(metric = c("median", "mean", "variance",
    "iqr")) {
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


.make_plot_for_genebuild_tx_long <- function(gene_single, mapping, counts, samples,
    top_n) {
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


.make_plot_for_geneaggregate_df_long <- function(df_long, agg_fun, pseudocount) {
    df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
    df_summary
}


.make_plot_for_genebuild_plot_from_summary <- function(df_summary, agg_label_unique,
    fill_limits = NULL, font_scale = 1) {
    # Calculate font sizes proportionally to output dimensions Reference: 12x8
    # inches (96 sq in) uses font_base=11 For other sizes, scale base font as:
    # base_font = 11 * sqrt(area/96) This ensures readability is maintained
    # across different output sizes

    base_font <- 11 * font_scale
    y_axis_font <- 12 * font_scale
    x_axis_font <- 14 * font_scale
    title_font <- 16 * font_scale
    legend_font <- 9 * font_scale

    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
        ggplot2::geom_tile(color = "black", linewidth = 0.3, width = 0.95, height = 0.92) +
        ggplot2::geom_vline(xintercept = 1.5, color = "white", linewidth = 1.5) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0,
        0)) + ggplot2::scale_fill_distiller(palette = "Blues", na.value = "lightgray",
        limits = fill_limits, name = "log2(expr)") + .theme_base(base_size = base_font) +
        ggplot2::labs(title = agg_label_unique, x = NULL, y = NULL, fill = "log2(expr)") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_font, hjust = 0.5, face = "bold"),
            plot.margin = ggplot2::margin(4, 4, 4, 4))
    
    # Apply axis label formatting and legend configuration using Phase 5 helpers
    p <- .format_axis_labels(p, x_size = x_axis_font, y_size = y_axis_font, 
                            y_face = "plain", bold_title = FALSE)
    p <- .configure_legend(p, position = "bottom", width_cm = 2, text_size = legend_font)
    p <- p + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
        barwidth = 10, barheight = 0.5, title.theme = ggplot2::element_text(size = title_font)))
    p
}


.make_plot_for_genecombine_patchwork <- function(plots, agg_label_unique) {
    # Use 2 columns (2 genes per row) with controlled spacing between rows
    n_cols <- 2
    n_rows <- ceiling(length(plots)/n_cols)

    # Build rows of 2 plots each with spacing between columns
    plot_rows <- list()
    for (row in seq_len(n_rows)) {
        start_idx <- (row - 1) * n_cols + 1
        end_idx <- min(row * n_cols, length(plots))
        row_plots <- plots[start_idx:end_idx]
        # Add right margin to first plot to create column spacing
        if (length(row_plots) >= 1) {
            row_plots[[1]] <- row_plots[[1]] + ggplot2::theme(plot.margin = ggplot2::margin(r = 1,
                unit = "cm"))
        }
        # Use patchwork composition (| for horizontal) to avoid scale conflicts
        if (length(row_plots) == 1) {
            row_combined <- row_plots[[1]]
        } else if (length(row_plots) == 2) {
            row_combined <- row_plots[[1]] | row_plots[[2]]  # Horizontal with patchwork
        } else {
            row_combined <- Reduce(function(x, y) x | y, row_plots)
        }
        plot_rows[[row]] <- row_combined
    }

    # Combine rows with spacers between them
    combined_elements <- list()
    heights_spec <- c()

    for (i in seq_along(plot_rows)) {
        combined_elements[[length(combined_elements) + 1]] <- plot_rows[[i]]
        heights_spec <- c(heights_spec, 1)

        if (i < length(plot_rows)) {
            # Add spacer between rows
            spacer <- ggplot2::ggplot() + ggplot2::theme_void()
            combined_elements[[length(combined_elements) + 1]] <- spacer
            heights_spec <- c(heights_spec, 0.17)  # Space between rows (reduced by half)
        }
    }

    # Combine all elements
    combined_plots_section <- Reduce(`/`, combined_elements) + patchwork::plot_layout(heights = heights_spec)

    # Add title spacer above plots (use / for vertical, not | for horizontal)
    spacer <- ggplot2::ggplot() + ggplot2::theme_void()
    title_row <- spacer | patchwork::plot_spacer()

    # Combine: title row on top, plot grid below
    combined <- title_row/combined_plots_section + patchwork::plot_annotation(title = "Transcript level expression",
        subtitle = paste0("Top genes with metric ", agg_label_unique), theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
            size = .font_sizes$title, face = "bold", margin = ggplot2::margin(t = 10,
                b = 10)), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = .font_sizes$subtitle,
            face = "italic", margin = ggplot2::margin(t = 5, b = 0.4)), legend.position = "bottom")) +
        patchwork::plot_layout(heights = c(0.045, 1), guides = "collect")
    combined
}


.make_plot_for_genecombine_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
    p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))

    # Use 2 columns (2 genes per row), auto-calculate rows
    ncol <- 2
    nrow_val <- ceiling(length(plots_nolegend)/ncol)

    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val,
        align = "hv")
    title_grob <- .create_title_grob("Transcript level expression",
        subtitle = paste0("Top genes with metric ", agg_label_unique),
        title_size = 18, subtitle_size = 14)
    # Add spacer between title and plots
    spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()
    result_plot <- cowplot::plot_grid(title_grob, spacer_grob, grid,
        legend, ncol = 1, rel_heights = c(0.09, 0.0015, 1, 0.08), align = "h",
        axis = "l")
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    }
    result_plot
}


.make_plot_for_genecombine_grid <- function(plots, output_file = NULL, agg_label_unique) {
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

    # Default to 2 columns (2 genes per row), adjust for smaller numbers
    ncol <- min(2, length(grobs))
    nrow <- ceiling(length(grobs)/ncol)

    # Create heights: title (0.5cm) + plot rows with gaps + legend (0.7cm)
    plot_heights <- list()
    for (i in seq_len(nrow)) {
        plot_heights[[length(plot_heights) + 1]] <- grid::unit(1, "null")
        if (i < nrow) {
            # Add gap after each row except the last (reduced by half)
            plot_heights[[length(plot_heights) + 1]] <- grid::unit(0.17, "cm")
        }
    }
    # Combine all heights properly using do.call
    all_heights <- c(list(grid::unit(0.55, "cm")), plot_heights, list(grid::unit(0.7,
        "cm")))
    heights <- do.call(grid::unit.c, all_heights)

    if (!is.null(output_file)) {
        # Adjust PNG dimensions based on layout
        png_width <- 800 * ncol
        png_height <- 480 * nrow
        png(filename = output_file, width = png_width, height = png_height, res = 150)
        .plot_transcript_grid_draw(grobs, agg_label_unique, legend_grob, ncol, heights,
            to_file = output_file)
        grDevices::dev.off()
        invisible(NULL)
    } else {
        # When no output file specified, create a temporary null device to
        # capture graphics This prevents R from creating Rplots.pdf as a
        # fallback when grid draws
        tmp_png <- tempfile(fileext = ".png")
        grDevices::png(tmp_png)
        on.exit({
            try(grDevices::dev.off(), silent = TRUE)
            if (file.exists(tmp_png)) unlink(tmp_png)
        }, add = TRUE)

        .plot_transcript_grid_draw(grobs, agg_label_unique, legend_grob, ncol, heights)
        invisible(NULL)
    }
}

#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_line
#' scale_y_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text geom_hline
#' geom_vline
#' @importFrom ggplot2 scale_color_manual scale_shape_manual geom_tile
#' scale_fill_gradient2
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("dimension", "variable", "contribution", "dim1", "dim2",
        "type", "coord_x", "coord_y"))
}



# Internal plot helpers

.format_label <- function(lbl) {
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

.prepare_ma_plot_df <- function(df, fold_col, mean_cols, x_label, y_label) {
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

.prepare_volcano_df <- function(diff_df, x_col = NULL, padj_col = "adjusted_p_values",
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

    x_label_formatted <- .format_label(x_col)
    padj_label_formatted <- .format_label(padj_col)

    list(df = df, x_col = x_col, padj_col = padj_col, x_label_formatted = x_label_formatted,
        padj_label_formatted = padj_label_formatted, title_use = title_use)
}


#' Plot Tsallis Divergence Effect Size Distribution
#'
#' Generate a histogram visualization of Tsallis divergence effect sizes
#' across genes,
#' showing the distribution of information-theoretic measures of isoform
#' switching.
#'
#' @param interaction_results A data frame containing LMM results merged
#' with per-q divergence estimates.
#' Must contain columns matching the pattern `effect_size_D_q*` (e.g.,
#' `effect_size_D_q0_5`, `effect_size_D_q1_0`).
#'   Typically the result from [.effect_sizes_divergence()].
#'
#' @param threshold Numeric. Effect size threshold for visual marking.
#' Default is 0.1 (information-theoretic significance level).
#'
#' @return If ggplot2 is available and `interaction_results` contains valid
#' data, returns a ggplot object.
#'   Otherwise returns NULL invisibly and prints an informative message.
#'
#' @details
#' The function visualizes the distribution of effect sizes using the median
#' q-value's divergence
#' (typically around q=1.0, close to Shannon entropy). The red dashed line
#' marks the default
#' information-theoretic significance threshold of D=0.1.
#'
#' @references
#' - Chanda et al. (2020). Information Theory in Computational Biology.
#' *Entropy*, 22(6), 627.
#' - Tsallis, C. (1988). Possible Generalization of Boltzmann-Gibbs
#' Statistics. *Journal of Statistical Physics*, 52(1), 479-487.
#'
#' @examples
#' # Create example interaction results with divergence effect sizes
#' set.seed(123)
#' interaction_results <- data.frame(
#'   gene = paste0('gene_', 1:20),
#'   effect_size_D_q0_5 = runif(20, 0, 0.3),
#'   effect_size_D_q1_0 = runif(20, 0, 0.2),
#'   effect_size_D_q1_5 = runif(20, 0, 0.25)
#' )
#' 
#' # Plot divergence distribution
#' .plot_divergence_distribution(interaction_results, threshold = 0.1)
#'

#' @noRd

.plot_divergence_distribution <- function(interaction_results, threshold = 0.1) {

    # Check for ggplot2
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        message("ggplot2 package required for plotting. Please install: install.packages('ggplot2')")
        return(invisible(NULL))
    }

    # Validate input
    if (is.null(interaction_results) || nrow(interaction_results) == 0) {
        message("Plot not generated: interaction_results is empty or NULL.")
        return(invisible(NULL))
    }

    # Find per-q effect size columns
    effect_cols <- grep("^effect_size_D_q", colnames(interaction_results), value = TRUE)

    if (length(effect_cols) == 0) {
        message("Plot not generated: no per-q effect size columns found in interaction_results.")
        message("Expected columns like 'effect_size_D_q0_5', 'effect_size_D_q1_0', etc.")
        return(invisible(NULL))
    }

    # Use median q effect size for visualization
    median_idx <- ceiling(length(effect_cols)/2)
    median_col <- effect_cols[median_idx]

    # Filter out NA/NaN/Inf values for plotting
    plot_data <- interaction_results[is.finite(interaction_results[[median_col]]),
        , drop = FALSE]

    if (nrow(plot_data) == 0) {
        message("Plot not generated: no valid (finite) effect sizes to plot.")
        return(invisible(NULL))
    }

    # Create visualization of effect size distribution with publication theme
    p_effect <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[median_col]])) +
        ggplot2::geom_histogram(binwidth = 0.02, fill = .palette_blue_red()[1], alpha = 0.7,
            color = "black") +
        ggplot2::labs(title = expression("Distribution of Tsallis Divergence (" ~
            D[q] ~ ") effect sizes across genes"),
            subtitle = "Information-theoretic measure respecting Tsallis multi-q entropy properties",
            x = bquote("Effect size (Tsallis Divergence" ~ D[q] ~ "; D >" ~ .(threshold) ~
                "= meaningful information separation)"),
            y = "Number of genes",
            caption = paste("Red dashed line: D =", threshold, "filtering threshold (information-theoretic significance for q-dependent entropy)")) +
        .theme_base(base_size = 11) +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray90"))
    
    # Add reference line using Phase 5 helper
    p_effect <- .add_reference_lines(p_effect, v_intercept = threshold, v_color = "red", v_size = 1)
    
    # Add threshold annotation
    p_effect <- p_effect + ggplot2::annotate("text", x = threshold, y = Inf, label = paste("Information\nthreshold\n(D=",
            threshold, ")", sep = ""), vjust = 1.5, hjust = -0.1, color = "red", size = 3.5)

    return(p_effect)
}




