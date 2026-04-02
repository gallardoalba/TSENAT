#' Plot Tsallis Entropy q-Curve
#'
#' Visualize Tsallis entropy (S_q) as a function of the diversity parameter
#' q across sample groups.
#' Supports three modes: aggregate q-curves (default), gene-specific
#' q-curves (when `gene` provided),
#' or bootstrap confidence interval bands.
#'
#' @param se A `SummarizedExperiment` returned by `.calculate_diversity()`
#' with diversity assay.
#' If pre-computed bootstrap confidence intervals are available
#' (ci_lower/ci_upper assays),
#' they will be displayed automatically. Otherwise, falls back to IQR
#' visualization.
#' @param assay_name Character; name of the assay to plot (default:
#' 'diversity').
#' @param condition_col Character or NULL; column name in colData indicating
#' group/sample type.
#' If NULL (default), reads from `@config$condition_col` when input is
#' TSENATAnalysis,
#'   otherwise defaults to 'sample_type'. Only used in aggregate and CI modes.
#' @param gene Character vector (optional); if provided, plot q-curves for
#' specified gene(s).
#' Overrides default aggregate behavior. When provided, uses median +/- SD
#' for each gene.
#' @param lm_res Data frame (optional); gene interaction test results with
#' `gene` column and
#'   p-value column. Accepts either:
#' - Results from `.calculate_lm_interaction()` (has `adj_p_interaction` or
#' `p_interaction` columns)
#' - Results from `.rank_test_q_condition()` (has `adj_p_value` or `p_value`
#' columns from Friedman/Wilcoxon tests)
#' If provided (and `gene` is NULL), plots top `n_top` genes ranked by
#' p-value.
#'   Useful for plotting significant genes from any interaction analysis.
#' @param n_top Integer or NULL; number of top genes to select from `lm_res`
#' when `gene` is NULL
#' (default: NULL). When NULL, defaults to showing the single most
#' significant gene (n_top=1),
#' providing a conservative view of the strongest effect. Set to a numeric
#' value to show that many top genes.
#' @param metric Character; when bootstrap CIs are NOT available, specifies
#' the spread metric to display.
#' Options: 'iqr' (default, Interquartile Range - more robust) or 'sd'
#' (Standard Deviation).
#' This parameter is ignored when bootstrap confidence intervals are
#' available.
#'   Default: 'iqr'.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @return
#' **Aggregate mode (gene=NULL, lm_res=NULL)**:
#' - If bootstrap CI assays available: A ggplot object showing median
#' entropy with bootstrap confidence interval bands.
#' - If no CI assays: A ggplot object showing median entropy with IQR
#' ribbons (automatic fallback).
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Single gene: A ggplot object showing median entropy +/- SD for that gene.
#' - Multiple genes: A grid plot object arranged in 2 rows x 2 columns with
#' a shared legend at the bottom.
#' The legend appears once beneath the grid, avoiding repetition across
#' subplots.
#'
#' @details
#' **Aggregate mode (default, gene=NULL, lm_res=NULL)**:
#' - Plots median Tsallis entropy +/- IQR across all genes for each group
#' - Works with any SummarizedExperiment from .calculate_diversity()
#' - Supports single or multiple q values and any number of groups
#' - No CI data required for basic plots; bootstrap CIs optional
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Plots q-curve separately for each selected gene
#' - Shows median entropy +/- SD (variance) for each gene across q-values
#' and groups
#' - When `lm_res` provided: automatically ranks genes and selects top
#' `n_top` by p-value
#' - Single gene: returns a ggplot object; multiple genes: returns a grid
#' plot (2 rows x 2 columns) with shared legend
#' - For multiple genes: legend appears once at the bottom of the grid to
#' avoid repetition and save space
#' - Useful for highlighting specific genes of interest or significant
#' discoveries
#' - Bootstrap mode not supported in gene-specific mode
#'
#' **Automatic Bootstrap CI Detection**:
#' - When computing divergence or diversity with bootstrap enabled, ci_lower
#' and ci_upper assays
#'   are added to the SummarizedExperiment.
#' - This function automatically detects these assays and displays bootstrap
#' confidence interval
#'   bands instead of IQR. No additional parameter needed.
#' - For confidence bands to appear, use `calculate_diversity_s4(...,
#' bootstrap=TRUE, nboot=1000)`
#'   or appropriate divergence function with bootstrap enabled.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point theme_minimal
#'   scale_color_manual scale_fill_manual labs theme element_text annotate
#' @importFrom dplyr filter group_by summarise pull
#' @importFrom SummarizedExperiment assayNames assay colData rowData
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' # Plot 7: Tsallis entropy q-curve (combined across all sample diversity)
#' data(readcounts)
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' 
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.5, 2, by = 0.5),
#' verbose = FALSE)
#' p <- plot_tsallis_q_curve_s4(analysis)
#' if (!is.null(p)) print(p)
#'
#' @export
plot_tsallis_q_curve_s4 <- function(se, assay_name = "diversity", condition_col = NULL,
    gene = NULL, lm_res = NULL, n_top = NULL, metric = "iqr", output_file = NULL) {
    # Validate metric parameter
    metric <- match.arg(tolower(metric), c("iqr", "sd"))

    # Load visualization dependencies (ggplot2, dplyr, tidyr, cowplot, etc.)
    .load_visualization_deps()

    # Convert TSENATAnalysis to combined SE if needed
    if (methods::is(se, "TSENATAnalysis")) {
        if (assay_name != "diversity") {
            stop("Assay '", assay_name, "' not found in SummarizedExperiment")
        }

        # Extract condition_col from config if not provided
        if (is.null(condition_col)) {
            if ("condition_col" %in% names(se@config)) {
                condition_col <- se@config$condition_col
            } else {
                condition_col <- "sample_type"
            }
        }

        se <- .prepare_combined_se(se)
        assay_name <- "diversity"
    }

    # Default condition_col if still NULL (for direct SE input)
    if (is.null(condition_col)) {
        condition_col <- "sample_type"
    }

    # Validate input
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("plot_tsallis_q_curve requires a SummarizedExperiment or TSENATAnalysis object")
    }

    if (!(assay_name %in% SummarizedExperiment::assayNames(se))) {
        stop("Assay '", assay_name, "' not found in SummarizedExperiment")
    }

    # Gene-specific mode
    if (!is.null(gene) || !is.null(lm_res)) {
        return(.plot_tsallis_gene_specific(se, assay_name, condition_col, gene, lm_res,
            n_top, metric, output_file))
    }

    # Aggregate or bootstrap mode
    long <- .prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
    if (nrow(long) == 0)
        stop("No tsallis values found in SummarizedExperiment")

    # AUTO-DETECT: Use bootstrap CIs if available, otherwise use specified
    # metric Check that CI assays exist AND have actual data (not all NAs)
    has_bootstrap_ci <- FALSE
    if ("ci_lower" %in% SummarizedExperiment::assayNames(se) && "ci_upper" %in% SummarizedExperiment::assayNames(se)) {
        # Verify assays have actual data
        ci_lower <- SummarizedExperiment::assay(se, "ci_lower")
        ci_upper <- SummarizedExperiment::assay(se, "ci_upper")
        n_valid_lower <- sum(!is.na(ci_lower))
        n_valid_upper <- sum(!is.na(ci_upper))
        has_bootstrap_ci <- (n_valid_lower > 0) && (n_valid_upper > 0)
    }

    if (has_bootstrap_ci) {
        return(.plot_tsallis_bootstrap_ci(se, long, output_file))
    }

    # Basic aggregate mode with specified metric
    .plot_tsallis_basic(long, metric, output_file)
}

# ============================================================================
# GENE-SPECIFIC Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_gene_specific <- function(se, assay_name, condition_col, gene, lm_res,
    n_top, metric, output_file) {
    require_pkgs(c("ggplot2", "dplyr", "cowplot", "SummarizedExperiment"))

    long <- .prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)

    # Check if bootstrap CIs are available
    has_bootstrap_ci <- "ci_lower" %in% SummarizedExperiment::assayNames(se) && "ci_upper" %in%
        SummarizedExperiment::assayNames(se)

    if (!("Gene" %in% colnames(long))) {
        if ("gene" %in% colnames(long)) {
            long$Gene <- long$gene
        } else {
            se_rownames <- rownames(se)
            if (!is.null(se_rownames) && length(se_rownames) > 0) {
                n_per_gene <- nrow(long)/length(se_rownames)
                long$Gene <- rep(se_rownames, each = n_per_gene)
            } else {
                stop("Cannot reconstruct Gene column from data")
            }
        }
    }

    # Resolve genes to plot
    if (is.null(gene)) {
        if (is.null(lm_res)) {
            stop("Either 'gene' or 'lm_res' (data.frame with 'gene' column) must be provided")
        }

        if (!is.data.frame(lm_res)) {
            stop("lm_res must be a data.frame with 'gene' column")
        }

        if (!("gene" %in% colnames(lm_res))) {
            stop("lm_res must contain a 'gene' column")
        }

        pcol <- NULL
        for (col in c("adj_p_interaction", "p_interaction", "adj_p_value", "p_value")) {
            if (col %in% colnames(lm_res)) {
                pcol <- col
                break
            }
        }
        if (is.null(pcol)) {
            stop("'lm_res' must contain one of: adj_p_interaction, p_interaction, adj_p_value, p_value")
        }

        genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[pcol]])]))
        n_genes_to_plot <- if (is.null(n_top))
            1 else n_top
        genes <- head(genes_ordered, n_genes_to_plot)
    } else {
        genes <- as.character(unlist(gene))
    }

    if (length(genes) == 0)
        stop("No genes selected for plotting")

    # Determine subtitle based on bootstrap availability and metric
    ci_subtitle <- if (has_bootstrap_ci) {
        "Median with Bootstrap 95% Confidence Intervals"
    } else {
        if (metric == "iqr") {
            "Median +/- Interquartile Range (IQR)"
        } else {
            "Median +/- Standard Deviation (SD)"
        }
    }

    # If bootstrap CIs available, use them for gene-specific plotting
    if (has_bootstrap_ci) {
        return(.plot_tsallis_gene_bootstrap_ci(se, long, genes, output_file))
    }

    # Plot single or multiple genes (fallback for non-bootstrap mode)
    make_plot_for_gene <- function(sel) {
        long_g <- long[as.character(long$Gene) == sel, , drop = FALSE]
        if (nrow(long_g) == 0)
            stop("Gene not found in assay: ", sel)

        stats_df <- .compute_gene_group_stats(long_g, metric = metric)

        p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum,
            ymin = central - spread, ymax = central + spread, fill = group), alpha = 0.2) +
            ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = central,
                color = group), linewidth = 1.3) + ggplot2::labs(title = sel, x = "q value",
            y = "Tsallis entropy", color = "Group", fill = "Group") + ggplot2::scale_color_manual(values = .palette_blue_red(),
            name = "Group") + ggplot2::scale_fill_manual(values = .palette_blue_red(),
            name = "Group") + .theme_base(base_size = 11) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
            size = 16, face = "bold"))
        p
    }

    if (length(genes) == 1) {
        p <- make_plot_for_gene(genes)
        # Add subtitle for single gene mode
        p <- p + ggplot2::labs(subtitle = ci_subtitle)
        return(p)
    }

    plots <- lapply(genes, make_plot_for_gene)
    names(plots) <- genes

    legend_obj <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom",
        legend.direction = "horizontal"))

    plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
    grid_with_plots <- do.call(cowplot::plot_grid, c(plots_no_legend, list(nrow = 2,
        ncol = 2)))

    title_plot <- cowplot::ggdraw() + cowplot::draw_label("Tsallis Entropy q-Curve Profile",
        fontface = "bold", size = 19) + cowplot::draw_label(ci_subtitle, fontface = "italic",
        size = 13, y = 0.25)

    grid_with_legend <- cowplot::plot_grid(title_plot, grid_with_plots, legend_obj,
        nrow = 3, rel_heights = c(0.12, 1, 0.08))

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = grid_with_legend, width = 12, height = 10,
            dpi = 100)
    }

    grid_with_legend
}

# ============================================================================
# BOOTSTRAP CI Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_bootstrap_ci <- function(se, long, output_file) {
    require_pkgs(c("ggplot2", "SummarizedExperiment"))

    # Bootstrap CI Visualization Methodology (From Literature: S115, S018)
    # =================================================================== Point
    # estimator: ORIGINAL OBSERVED MEDIAN (robust measure of central tendency)
    # CI bounds: BOOTSTRAP PERCENTILE METHOD (2.5th and 97.5th percentiles)
    # Aggregation: MEDIAN of per-sample bootstrap CI bounds across samples in
    # group INTERPRETATION: - Line represents observed median Tsallis entropy
    # (actual data) - Ribbon represents bootstrap 95% confidence interval
    # around the estimate - If line falls outside ribbon: indicates asymmetric
    # bootstrapping distribution (NOT a problem - reveals non-normality in the
    # resampling distribution)

    long$q <- as.numeric(as.character(long$q))
    unique_q <- sort(unique(long$q))

    groups <- unique(sort(long$group))
    if (length(unique_q) < 2) {
        stop("Need at least 2 q values for q-curve analysis")
    }
    if (length(groups) != 2) {
        stop("Expected exactly 2 groups for bootstrap comparison")
    }

    plot_df <- .bootstrap_aggregate_ci(se, long)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = median, color = group,
        fill = group)) + ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower,
        ymax = ci_upper), alpha = 0.15, color = NA) + ggplot2::scale_color_manual(values = .palette_blue_red(),
        name = "Group") + ggplot2::scale_fill_manual(values = .palette_blue_red(),
        name = "Group") + .theme_base(base_size = 11) + ggplot2::labs(title = "Tsallis Entropy Across Diversity Scales (q-spectrum)",
        subtitle = "Observed median (line) with bootstrap 95% percentile CI (shaded band)",
        x = "q value", y = expression("Tsallis entropy (" * S[q] * ")"), color = "Group",
        fill = "Group")

    if (length(groups) == 1) {
        p <- p + ggplot2::theme(legend.position = "none")
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = p, width = 12, height = 7.2, dpi = 100)
    }

    p
}

# ============================================================================
# GENE-SPECIFIC BOOTSTRAP CI Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_gene_bootstrap_ci <- function(se, long, genes, output_file) {
    require_pkgs(c("ggplot2", "dplyr", "cowplot", "SummarizedExperiment"))

    # Filter long data to selected genes
    long$q <- as.numeric(as.character(long$q))
    long_subset <- long[as.character(long$Gene) %in% genes, ]

    # Extract bootstrap CI bounds from assays
    ci_lower_assay <- SummarizedExperiment::assays(se)[["ci_lower"]]
    ci_upper_assay <- SummarizedExperiment::assays(se)[["ci_upper"]]

    # Get row indices for selected genes
    gene_indices <- match(genes, rownames(se))
    valid_indices <- gene_indices[!is.na(gene_indices)]

    # Extract CIs for selected genes
    if (length(valid_indices) > 0) {
        ci_lower <- as.matrix(ci_lower_assay[valid_indices, ])
        ci_upper <- as.matrix(ci_upper_assay[valid_indices, ])

        # Map CIs to long format with groups
        plot_df <- .prepare_gene_ci_data(long_subset, ci_lower, ci_upper, genes)
    } else {
        # Fallback: no CIs found
        return(.plot_tsallis_basic_gene(long_subset, genes, output_file = output_file))
    }

    # Create plots for each gene
    make_gene_plot <- function(g) {
        plot_data <- plot_df[plot_df$Gene == g, ]

        if (nrow(plot_data) == 0) {
            # Fallback for genes without CI data
            return(NULL)
        }

        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = q, y = median, color = group,
            fill = group)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
            alpha = 0.15, color = NA) + ggplot2::geom_line(linewidth = 1.2) + ggplot2::scale_color_manual(values = .palette_blue_red(),
            name = "Group") + ggplot2::scale_fill_manual(values = .palette_blue_red(),
            name = "Group") + .theme_base(base_size = 11) + ggplot2::labs(title = g,
            x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14,
                face = "bold"))

        p
    }

    # Handle single vs multiple genes
    if (length(genes) == 1) {
        p <- make_gene_plot(genes[1])
        if (!is.null(p)) {
            p <- p + ggplot2::labs(subtitle = "Median with Bootstrap 95% Confidence Intervals")
        }
        return(p)
    }

    # Multiple genes: create grid
    plots <- lapply(genes, make_gene_plot)
    plots <- plots[!vapply(plots, is.null, FUN.VALUE = logical(1))]

    if (length(plots) == 0) {
        stop("No valid genes found for plotting")
    }

    legend_obj <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom",
        legend.direction = "horizontal"))

    plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
    grid_with_plots <- do.call(cowplot::plot_grid, c(plots_no_legend, list(nrow = 2,
        ncol = 2)))

    title_plot <- cowplot::ggdraw() + cowplot::draw_label("Tsallis Entropy q-Curve Profile",
        fontface = "bold", size = 19) + cowplot::draw_label("Median with Bootstrap 95% Confidence Intervals",
        fontface = "italic", size = 13, y = 0.25)

    grid_with_legend <- cowplot::plot_grid(title_plot, grid_with_plots, legend_obj,
        nrow = 3, rel_heights = c(0.12, 1, 0.08))

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = grid_with_legend, width = 12, height = 10,
            dpi = 100)
    }

    grid_with_legend
}

# ============================================================================
# FALLBACK: GENE-SPECIFIC BASIC Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_basic_gene <- function(long, genes, metric = "iqr", output_file) {
    require_pkgs(c("ggplot2", "dplyr", "cowplot"))

    metric <- match.arg(tolower(metric), c("iqr", "sd"))

    # Determine subtitle based on metric
    subtitle <- if (metric == "iqr") {
        "Median +/- Interquartile Range (IQR)"
    } else {
        "Median +/- Standard Deviation (SD)"
    }

    make_gene_plot <- function(g) {
        long_g <- long[as.character(long$Gene) == g, ]

        # Calculate spread based on metric choice
        if (metric == "iqr") {
            stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, q), median = median(tsallis,
                na.rm = TRUE), spread = stats::IQR(tsallis, na.rm = TRUE)/2, .groups = "drop")
        } else {
            stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, q), median = median(tsallis,
                na.rm = TRUE), spread = sqrt(stats::var(tsallis, na.rm = TRUE)),
                .groups = "drop")
        }

        p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = q, y = median, color = group,
            fill = group)) + ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_ribbon(ggplot2::aes(ymin = median -
            spread, ymax = median + spread), alpha = 0.2, color = NA) + ggplot2::scale_color_manual(values = .palette_blue_red(),
            name = "Group") + ggplot2::scale_fill_manual(values = .palette_blue_red(),
            name = "Group") + .theme_base(base_size = 11) + ggplot2::labs(title = g,
            x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14,
                face = "bold"))

        p
    }

    if (length(genes) == 1) {
        p <- make_gene_plot(genes[1])
        p <- p + ggplot2::labs(subtitle = subtitle)
        return(p)
    }

    plots <- lapply(genes, make_gene_plot)
    legend_obj <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom",
        legend.direction = "horizontal"))

    plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
    grid_with_plots <- do.call(cowplot::plot_grid, c(plots_no_legend, list(nrow = 2,
        ncol = 2)))

    title_plot <- cowplot::ggdraw() + cowplot::draw_label("Tsallis Entropy q-Curve Profile",
        fontface = "bold", size = 19) + cowplot::draw_label(subtitle, fontface = "italic",
        size = 13, y = 0.25)

    grid_with_legend <- cowplot::plot_grid(title_plot, grid_with_plots, legend_obj,
        nrow = 3, rel_heights = c(0.12, 1, 0.08))

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = grid_with_legend, width = 12, height = 10,
            dpi = 100)
    }

    grid_with_legend
}

# ============================================================================
# BASIC AGGREGATE Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_basic <- function(long, metric = "iqr", output_file) {
    require_pkgs("ggplot2")

    metric <- match.arg(tolower(metric), c("iqr", "sd"))
    long$q <- as.numeric(as.character(long$q))

    # Calculate spread based on metric choice
    if (metric == "iqr") {
        stats_df <- dplyr::summarise(dplyr::group_by(long, group, q), median = median(tsallis,
            na.rm = TRUE), spread = stats::IQR(tsallis, na.rm = TRUE)/2, .groups = "drop")
        subtitle <- "Median +/- Interquartile Range (IQR)"
    } else {
        stats_df <- dplyr::summarise(dplyr::group_by(long, group, q), median = median(tsallis,
            na.rm = TRUE), spread = sqrt(stats::var(tsallis, na.rm = TRUE)), .groups = "drop")
        subtitle <- "Median +/- Standard Deviation (SD)"
    }

    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = q, y = median, color = group,
        fill = group)) + ggplot2::geom_line(linewidth = 1.3) + ggplot2::geom_ribbon(ggplot2::aes(ymin = median -
        spread, ymax = median + spread), alpha = 0.2, color = NA) + ggplot2::scale_color_manual(values = .palette_blue_red(),
        name = "Group") + ggplot2::scale_fill_manual(values = .palette_blue_red(),
        name = "Group") + .theme_base(base_size = 11) + ggplot2::labs(title = "Tsallis Entropy Across Diversity Scales (q-spectrum)",
        subtitle = subtitle, x = "q value", y = expression("Tsallis entropy (" *
            S[q] * ")"), color = "Group", fill = "Group")

    if (length(unique(long$group)) == 1) {
        p <- p + ggplot2::theme(legend.position = "none")
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = p, width = 12, height = 7.2, dpi = 100)
    }

    p
}
