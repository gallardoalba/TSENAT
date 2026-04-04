#' Plot Composition & Theme Helpers for TSENAT Visualization
#'
#' This module provides utilities for combining multiple plots,
#' managing themes, and creating consistent ggplot2 scales.
#' Enables modular, reusable visualization composition.
#'
#' @name plot_helpers

#' @noRd
NULL

# ============================================================================
# MULTI-PLOT COMPOSITION: Patchwork
# ============================================================================

#' Combine Transcript Plots with Patchwork
#'
#' Uses patchwork package to compose transcript-level plots
#' in a 2-column grid with proper spacing and legends.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Combined patchwork plot object.
#'

#' @noRd

.combine_plots_patchwork <- function(plots, agg_label_unique) {
    require_pkgs("patchwork")

    # Use 2 columns (2 genes per row) with controlled spacing
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

        # Use patchwork composition (| for horizontal)
        if (length(row_plots) == 1) {
            row_combined <- row_plots[[1]]
        } else if (length(row_plots) == 2) {
            row_combined <- row_plots[[1]] | row_plots[[2]]
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
            heights_spec <- c(heights_spec, 0.17)  # Reduced spacing
        }
    }

    # Combine all elements
    combined_plots_section <- Reduce(`/`, combined_elements) + patchwork::plot_layout(heights = heights_spec)

    # Create title annotation
    combined <- patchwork::plot_spacer()/combined_plots_section + patchwork::plot_annotation(title = "Transcript level expression",
        subtitle = paste0("Top genes with metric ", agg_label_unique), theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
            size = .font_sizes$title, face = "bold", margin = ggplot2::margin(t = 10,
                b = 10)), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = .font_sizes$subtitle,
            face = "italic", margin = ggplot2::margin(t = 5, b = 0.4)))) + patchwork::plot_layout(heights = c(0.045,
        1), guides = "collect")

    return(combined)
}

# ============================================================================
# MULTI-PLOT COMPOSITION: Cowplot
# ============================================================================

#' Combine Transcript Plots with Cowplot
#'
#' Uses cowplot package to compose transcript-level plots
#' with proper titles and legend handling.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param output_file Character: optional file path to save result.
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Combined plot object (or invisible NULL if output_file provided).
#'

#' @noRd

.combine_plots_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
    require_pkgs("cowplot")

    # Extract legend from first plot
    p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        pp + ggplot2::theme(legend.position = "none")
    })

    # Compose in 2 columns
    ncol <- 2
    nrow_val <- ceiling(length(plots_nolegend)/ncol)

    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val,
        align = "hv")

    # Create title and subtitle
    title_grob <- cowplot::ggdraw() + cowplot::draw_label("Transcript level expression",
        fontface = "bold", x = 0.5, hjust = 0.5, size = .font_sizes$title)

    subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(paste0("Top genes with metric ",
        agg_label_unique), fontface = "italic", x = 0.5, hjust = 0.5, size = .font_sizes$subtitle,
        color = "gray40")

    spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()

    # Combine all elements
    result_plot <- cowplot::plot_grid(title_grob, subtitle_grob, spacer_grob, grid,
        legend, ncol = 1, rel_heights = c(0.05, 0.04, 0.0015, 1, 0.08), align = "h",
        axis = "l")

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        return(invisible(NULL))
    }

    return(result_plot)
}

# ============================================================================
# MULTI-PLOT COMPOSITION: Grid (Low-level)
# ============================================================================

#' Combine Transcript Plots with Grid
#'
#' Uses base grid package for low-level composition.
#' Useful for fine-grained control over plot layout.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param output_file Character: optional file path to save as PNG.
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Invisible NULL. Outputs to file or current device.
#'

#' @noRd

.combine_plots_grid <- function(plots, output_file = NULL, agg_label_unique) {
    require_pkgs("grid")

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        pp + ggplot2::theme(legend.position = "none")
    })

    # Convert to grobs
    grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)

    # Extract legend from first plot
    g_full <- ggplot2::ggplotGrob(plots[[1]])
    legend_idx <- which(vapply(g_full$grobs, function(x) x$name, character(1)) ==
        "guide-box")

    legend_grob <- if (length(legend_idx)) {
        g_full$grobs[[legend_idx[1]]]
    } else {
        NULL
    }

    # Default to 2 columns
    ncol <- min(2, length(grobs))
    nrow <- ceiling(length(grobs)/ncol)

    # Create height specification
    plot_heights <- list()
    for (i in seq_len(nrow)) {
        plot_heights[[length(plot_heights) + 1]] <- grid::unit(1, "null")
        if (i < nrow) {
            plot_heights[[length(plot_heights) + 1]] <- grid::unit(0.17, "cm")
        }
    }

    all_heights <- c(list(grid::unit(0.55, "cm")), plot_heights, list(grid::unit(0.7,
        "cm")))
    heights <- do.call(grid::unit.c, all_heights)

    if (!is.null(output_file)) {
        png_width <- 800 * ncol
        png_height <- 480 * nrow
        png(filename = output_file, width = png_width, height = png_height, res = 150)
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights,
            to_file = output_file)
        dev.off()
        return(invisible(NULL))
    } else {
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights)
        return(invisible(NULL))
    }
}

# ============================================================================
# GRID DRAWING UTILITY
# ============================================================================

#' Draw Transcript Grid Layout
#'
#' Internal utility for drawing transcript plots using base grid system.
#'
#' @param grobs List of ggplot grobs (one per plot).
#' @param agg_label_unique Character: metric label.
#' @param legend_grob Grid grob: legend to include.
#' @param ncol Integer: number of columns.
#' @param heights Unit specification: row heights.
#' @param to_file Character: optional file path (internal use).
#'

#' @noRd

.draw_transcript_grid <- function(grobs, agg_label_unique, legend_grob, ncol, heights,
    to_file = NULL) {
    require_pkgs("grid")

    nrow <- ceiling(length(grobs)/ncol)

    # Create viewport structure
    vp_top <- grid::viewport(x = 0, y = 0.95, width = 1, height = 0.05, just = c("left",
        "bottom"), name = "title")

    grid::pushViewport(vp_top)
    grid::grid.text("Transcript level expression", x = 0.5, y = 0.5, gp = grid::gpar(fontsize = .font_sizes$title,
        fontface = "bold"))
    grid::upViewport()

    # Draw plot grid
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow, ncol = ncol)))

    for (i in seq_along(grobs)) {
        row <- ((i - 1)%/%ncol) + 1
        col <- ((i - 1)%%ncol) + 1
        grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }

    grid::upViewport()
}

# ============================================================================
# SCALE CREATORS
# ============================================================================

#' Create Discrete Color Scale for Heatmap
#'
#' Creates a discrete color scale using TSENAT's standard palette.
#'
#' @param palette Character: 'blue_red' (default), 'continuous_diverging',
#' or custom.
#' @param direction Integer: 1 (default) or -1 to reverse colors.
#' @param name Character: legend title.
#'
#' @return ggplot2 scale object (ggplot2::scale_color_manual or similar).
#'

#' @noRd

.create_color_scale <- function(palette = "blue_red", direction = 1, name = NULL) {
    require_pkgs("ggplot2")

    if (palette == "blue_red") {
        colors <- .palette_blue_red()
    } else if (palette == "continuous_diverging") {
        colors <- .palette_continuous_diverging()
    } else if (is.character(palette)) {
        colors <- palette
    } else {
        colors <- .palette_blue_red()
    }

    if (direction == -1) {
        colors <- rev(colors)
    }

    if (length(colors) > 1) {
        ggplot2::scale_color_manual(values = colors, name = name)
    } else {
        NULL
    }
}

#' Create Fill Scale for Heatmap
#'
#' Creates a fill scale for heatmap-style plots.
#'
#' @param palette Character: 'blue_red' (default) or 'continuous_diverging'.
#' @param direction Integer: 1 (default) or -1 to reverse.
#' @param name Character: legend title.
#' @param breaks Integer: number of color breaks (default: 50).
#'
#' @return ggplot2 scale object.
#'

#' @noRd

.create_fill_scale <- function(palette = "blue_red", direction = 1, name = NULL,
    breaks = 50) {
    require_pkgs("ggplot2")

    if (palette == "continuous_diverging") {
        colors <- .palette_continuous_diverging(n = breaks)
    } else {
        colors <- .palette_blue_red()
    }

    if (direction == -1) {
        colors <- rev(colors)
    }

    ggplot2::scale_fill_gradient(low = colors[1], high = colors[length(colors)],
        name = name)
}

# ============================================================================
# THEME UTILITIES
# ============================================================================

#' Apply TSENAT Base Theme
#'
#' Applies standard TSENAT styling: minimal theme with centered titles.
#'
#' @param base_size Integer: base font size (default: 11).
#' @param color_palette Character: 'blue_red' (default) or other.
#'
#' @return ggplot2 theme object.
#'

#' @noRd

.apply_tsenat_theme <- function(base_size = 11, color_palette = "blue_red") {
    require_pkgs("ggplot2")

    theme_result <- .theme_base(base_size = base_size) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = .font_sizes$title, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = .font_sizes$subtitle, face = "italic"))

    return(theme_result)
}

#' Modify Plot Title and Subtitle
#'
#' Convenience function to update title/subtitle in existing plot.
#'
#' @param plot ggplot2 object.
#' @param title Character: new title.
#' @param subtitle Character: new subtitle.
#' @param title_size Integer: title font size (default: from constants).
#' @param subtitle_size Integer: subtitle font size (default: from constants).
#'
#' @return Modified ggplot2 object.
#'

#' @noRd

.set_plot_title <- function(plot, title = NULL, subtitle = NULL, title_size = .font_sizes$title,
    subtitle_size = .font_sizes$subtitle) {
    require_pkgs("ggplot2")

    if (!is.null(title)) {
        plot <- plot + ggplot2::labs(title = title)
    }

    if (!is.null(subtitle)) {
        plot <- plot + ggplot2::labs(subtitle = subtitle)
    }

    plot <- plot + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = title_size, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = subtitle_size, face = "italic"))

    return(plot)
}

# ============================================================================
# HEATMAP STYLING
# ============================================================================

#' Create Heatmap with Unified Styling
#'
#' Wraps pheatmap with TSENAT standard settings.
#'
#' @param mat Matrix: data to plot.
#' @param title Character: plot title.
#' @param colors Vector: color specification for heatmap.
#' @param breaks Vector: color break points.
#' @param fontsize_row Integer: row label font size.
#' @param fontsize_col Integer: column label font size.
#' @param ... Additional arguments passed to pheatmap.
#'
#' @return Heatmap object from pheatmap.
#'

#' @noRd

.create_tsenat_heatmap <- function(mat, title = NULL, colors = NULL, breaks = NULL,
    fontsize_row = 11, fontsize_col = 11, ...) {
    require_pkgs("pheatmap")

    if (is.null(colors)) {
        colors <- .palette_continuous_diverging(n = 100)
    }

    pheatmap::pheatmap(mat, main = title, color = colors, breaks = breaks, fontsize_row = fontsize_row,
        fontsize_col = fontsize_col, ...)
}

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

    require_pkgs(c("SummarizedExperiment", "dplyr"))

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

    require_pkgs("dplyr")

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

    require_pkgs("dplyr")

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
# DATA VALIDATION & QUALITY CHECKS
# ============================================================================

#' Validate Diversity SummarizedExperiment
#'
#' Checks that SE has required structure for diversity visualization.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param check_metadata Logical: also validate metadata? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error with message.
#'

#' @noRd

.validate_diversity_se <- function(se, check_metadata = TRUE) {

    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    if (nrow(se) == 0) {
        stop("SummarizedExperiment has no rows (samples)", call. = FALSE)
    }

    if (ncol(se) == 0) {
        stop("SummarizedExperiment has no columns (genes)", call. = FALSE)
    }

    # Check for diversity assay
    assay_names <- SummarizedExperiment::assayNames(se)
    if (!("diversity" %in% assay_names)) {
        stop("Required 'diversity' assay not found. ", "Available: ", paste(assay_names,
            collapse = ", "), call. = FALSE)
    }

    # Check for valid data
    div_mat <- SummarizedExperiment::assay(se, "diversity")
    if (all(is.na(div_mat))) {
        stop("All diversity values are NA", call. = FALSE)
    }

    if (check_metadata) {
        # Check for at least one q-value
        meta <- S4Vectors::metadata(se)
        if (!("q" %in% names(meta)) || length(meta$q) == 0) {
            warning("q-values not found in SE metadata", call. = FALSE)
        }
    }

    return(TRUE)
}

#' Validate Results Data Frame for Gene Selection
#'
#' Checks that results DataFrame has required columns.
#'
#' @param results Data frame (LM results, effect sizes, etc.).
#' @param require_pvalue Logical: check for p-value column? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error.
#'

#' @noRd

.validate_results_df <- function(results, require_pvalue = TRUE) {

    if (!is.data.frame(results)) {
        stop("results must be a data frame", call. = FALSE)
    }

    if (nrow(results) == 0) {
        stop("results data frame is empty", call. = FALSE)
    }

    # Check for gene column
    gene_cols <- c("gene_id", "gene", "gene_name")
    has_gene <- any(gene_cols %in% colnames(results))
    if (!has_gene) {
        stop("No gene identifier column found. ", "Expected one of: ", paste(gene_cols,
            collapse = ", "), call. = FALSE)
    }

    # Check for p-value column
    if (require_pvalue) {
        p_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
        has_pval <- any(p_cols %in% colnames(results))
        if (!has_pval) {
            stop("No p-value column found. ", "Expected one of: ", paste(p_cols,
                collapse = ", "), call. = FALSE)
        }
    }

    return(TRUE)
}

# ============================================================================
# FORMATTING & UTILITY FUNCTIONS
# ============================================================================

#' Format P-Value for Display
#'
#' Converts p-value to formatted string (scientific or threshold).
#'
#' @param pval Numeric p-value.
#' @param threshold Numeric: cutoff for '< threshold' format (default: 0.001).
#' @param digits Integer: decimal places for scientific notation (default: 2).
#'
#' @return Character string formatted p-value.
#'

#' @noRd

.format_pvalue <- function(pval, threshold = 0.001, digits = 2) {

    if (is.na(pval)) {
        return("NA")
    }

    if (pval < threshold) {
        return(paste0("< ", threshold))
    }

    return(format(pval, scientific = TRUE, digits = digits))
}

#' Format Q-Value Label
#'
#' Converts numeric q-value to display label (e.g., 'q = 1.0').
#'
#' @param q_val Numeric q-value.
#' @param prefix Character: prefix for label (default: 'q').
#'
#' @return Character string label.
#'

#' @noRd

.format_q_label <- function(q_val, prefix = "q") {
    if (is.na(q_val)) {
        return("NA")
    }
    return(sprintf("%s = %.2f", prefix, as.numeric(q_val)))
}

#' Format Label for Display
#'
#' Converts underscored/raw column names to readable labels.
#' Replaces underscores with spaces and formats capitalization.
#'
#' @param lbl Character: label to format (may contain underscores).
#'
#' @return Character string, properly capitalized.
#'

#' @noRd

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

# ============================================================================
# TRANSCRIPT PLOTTING HELPERS: Data Preparation
# ============================================================================

#' Prepare Inputs for Transcript-Level Plotting
#'
#' Normalizes and validates counts, samples, and tx2gene mapping.
#' Creates aggregation function based on chosen metric.
#'
#' @param counts Matrix or data.frame with transcripts as rows, samples as
#' columns.
#'   Can also be a \code{SummarizedExperiment}.
#' @param readcounts Character: name of assay in SE (if counts is SE).
#' Default: NULL.
#' @param samples Character vector: sample group assignments (optional).
#' @param coldata Character/data.frame: sample metadata (optional).
#' @param condition_col Character: column name for grouping
#'   in coldata (default: 'sample_type').
#' @param tx2gene data.frame/character: Transcript-to-gene mapping with columns
#'   'Transcript' and 'Gen'. Can be file path or data.frame.
#' @param res Optional data.frame with results (gene names and p-values).
#' @param top_n Integer: number of transcripts to select.
#' @param pseudocount Numeric: pseudocount for log transformation (default: 0).
#' @param output_file Character: file path for saving plot (optional).
#' @param metric Character: aggregation metric
#'   ('median' [default], 'mean', 'variance', 'iqr').
#'
#' @return List with elements:
#'   - counts: normalized count matrix
#'   - samples: sample group assignments
#'   - mapping: tx2gene data.frame
#'   - metric_choice: chosen metric
#'   - agg_fun: aggregation function
#'   - agg_label_unique: metric label for display
#'   - top_n: number of transcripts
#'   - pseudocount: pseudocount value
#'   - output_file: output file path (if provided)
#'

#' @noRd

.prepare_transcript_inputs <- function(counts, readcounts = NULL, samples = NULL,
    coldata = NULL, condition_col = "condition", tx2gene = NULL, res = NULL, top_n = NULL,
    pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance",
        "iqr")) {

    require_pkgs(c("SummarizedExperiment", "S4Vectors"))

    # Handle SummarizedExperiment input
    if (inherits(counts, "SummarizedExperiment")) {
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

    # Validate counts
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop("`counts` must be a matrix, data.frame, or SummarizedExperiment", call. = FALSE)
    }

    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) {
        stop("`counts` must have rownames (transcript identifiers)", call. = FALSE)
    }

    # Infer samples from coldata if needed
    if (is.null(samples)) {
        if (!is.null(coldata)) {
            samples <- .infer_samples_from_coldata(coldata, counts, condition_col)
        } else {
            stop("Either 'samples' or 'coldata' must be provided", call. = FALSE)
        }
    }

    # Validate and normalize tx2gene
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided", call. = FALSE)
    }
    mapping <- .read_tx2gene(tx2gene)

    if (length(samples) != ncol(counts)) {
        stop("Length of `samples` must equal columns in `counts`", call. = FALSE)
    }

    # Create aggregation function
    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) stats::var(x,
            na.rm = TRUE), iqr = function(x) stats::IQR(x, na.rm = TRUE))

    agg_label_metric <- if (metric_choice == "iqr")
        "IQR" else metric_choice
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    agg_label_unique <- agg_label

    list(counts = counts, samples = samples, mapping = mapping, metric_choice = metric_choice,
        agg_fun = agg_fun, agg_label_unique = agg_label_unique, top_n = top_n, pseudocount = pseudocount,
        output_file = output_file)
}

#' Read and Validate tx2gene Mapping
#'
#' Reads transcript-to-gene mapping from file or data.frame.
#' Validates required columns: 'Transcript' and 'Gen'.
#'
#' @param tx2gene Character (file path) or data.frame mapping.
#'
#' @return data.frame with columns 'Transcript' and 'Gen'.
#'

#' @noRd

.read_tx2gene <- function(tx2gene) {
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided as file path or data.frame", call. = FALSE)
    }

    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be file path or data.frame", call. = FALSE)
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
        stop("tx2gene must have columns 'Transcript' and 'Gen'", call. = FALSE)
    }

    return(mapping)
}

#' Infer Samples from Column Metadata
#'
#' Extracts sample group assignments from coldata.
#' Aligns sample IDs from coldata to counts columns.
#'
#' @param coldata Character (file path) or data.frame with sample metadata.
#' @param counts Count matrix (for column name alignment).
#' @param condition_col Character: column name for grouping variable.
#'
#' @return Character vector of sample group assignments.
#'

#' @noRd

.infer_samples_from_coldata <- function(coldata, counts, condition_col) {
    if (is.character(coldata) && length(coldata) == 1) {
        if (!file.exists(coldata)) {
            stop("coldata file not found: ", coldata, call. = FALSE)
        }
        cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(coldata)) {
        cdf <- coldata
    } else {
        stop("`coldata` must be file path or data.frame", call. = FALSE)
    }

    # Try row-indexed matching first
    if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
        return(as.character(cdf[colnames(counts), condition_col]))
    }

    # Try sample ID column matching
    sample_id_cols <- c("sample", "Sample", "sample_id", "id")
    sid <- intersect(sample_id_cols, colnames(cdf))

    if (length(sid) > 0) {
        sid <- sid[1]
        if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
            stop("coldata sample ID column doesn't match counts columns", call. = FALSE)
        }
        row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
        return(as.character(cdf[[condition_col]][row_ix]))
    }

    stop("Could not match coldata to counts. Provide rownames or sample ID column.",
        call. = FALSE)
}

#' Create Aggregation Function
#'
#' Builds an aggregation function based on chosen metric.
#'
#' @param metric Character: 'median' (default), 'mean', 'variance', or 'iqr'.
#'
#' @return List with:
#'   - metric_choice: the selected metric
#'   - agg_fun: function that computes the metric
#'   - agg_label_unique: display label
#'

#' @noRd

.create_aggregation_function <- function(metric = c("median", "mean", "variance",
    "iqr")) {
    metric_choice <- match.arg(metric)

    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) stats::var(x,
            na.rm = TRUE), iqr = function(x) stats::IQR(x, na.rm = TRUE))

    agg_label_metric <- if (metric_choice == "iqr")
        "IQR" else metric_choice
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    agg_label_unique <- agg_label

    list(metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique)
}

#' Build Long-Format Transcript Data
#'
#' Transforms wide count matrix to long-format data.frame
#' for ggplot visualization.
#'
#' @param gene_single Character: single gene identifier.
#' @param mapping data.frame: tx2gene mapping with Transcript and Gen columns.
#' @param counts Matrix: transcript count matrix.
#' @param samples Character vector: sample group assignments.
#' @param top_n Integer: limit to top N transcripts (optional).
#'
#' @return List with:
#'   - df_long: long-format data.frame (columns: tx, sample, expr, group)
#'   - txs: selected transcript identifiers
#'

#' @noRd

.build_transcript_long <- function(gene_single, mapping, counts, samples, top_n = NULL) {
    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))

    if (length(txs) == 0) {
        stop("No transcripts found for gene: ", gene_single, call. = FALSE)
    }

    if (!is.null(top_n)) {
        txs <- head(txs, top_n)
    }

    # Create long-format data
    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)

    require_pkgs("tidyr")
    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))

    list(df_long = df_long, txs = txs)
}

#' Aggregate Long-Format Transcript Data
#'
#' Summarizes expression by transcript and group.
#'
#' @param df_long Long-format data.frame from \code{build_transcript_long}.
#' @param agg_fun Function: aggregation function (e.g., median, mean).
#' @param pseudocount Numeric: pseudocount for log transformation.
#'
#' @return data.frame with columns: tx, group, expr, log2expr.
#'

#' @noRd

.aggregate_transcript_data <- function(df_long, agg_fun, pseudocount = 0) {
    df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
    return(df_summary)
}

#' Select Top Genes from Results  
#'
#' Extracts top genes from results by p-value.
#'
#' @param res data.frame: results with gene and p-value columns.
#' @param top_n Integer: number of genes to select.
#'
#' @return Character vector of top gene IDs.
#'

#' @noRd

.select_genes_from_results <- function(res, top_n) {
    if (is.null(res)) {
        stop("Either 'gene' or 'res' must be provided", call. = FALSE)
    }

    if (!("genes" %in% colnames(res))) {
        stop("res must contain a 'genes' column", call. = FALSE)
    }

    # Find p-value column (ordered by preference)
    p_cols <- c("padj", "adjusted_p_values", "pvalue", "raw_p_values")
    p_col <- intersect(p_cols, colnames(res))[1]

    if (is.na(p_col)) {
        # No p-value column; just return gene order
        ord <- seq_len(nrow(res))
    } else {
        ord <- order(res[[p_col]], na.last = NA)
    }

    genes_sel <- as.character(res$genes[ord])
    genes_sel <- unique(genes_sel)
    return(head(genes_sel, top_n))
}

# ============================================================================
# PLOT TSALLIS Q-CURVE HELPERS
# ============================================================================

#' Extract diversity objects from analysis results
#'
#' @param div_list List of diversity results (SummarizedExperiment or matrix)
#' @return List with: objects, q_names, first_se, bootstrap_ci_available

#' @noRd
.extract_diversity_objects <- function(div_list) {
    require_pkgs(c("SummarizedExperiment"))
    
    combined_assays_dict <- list()
    first_se <- NULL
    bootstrap_ci_available <- FALSE

    for (q_name in names(div_list)) {
        obj <- div_list[[q_name]]
        if (methods::is(obj, "SummarizedExperiment")) {
            if (is.null(first_se)) {
                first_se <- obj
            }
            # Check if bootstrap CIs are available
            si <- SummarizedExperiment::assayNames(obj)
            if ("ci_lower" %in% si && "ci_upper" %in% si) {
                bootstrap_ci_available <- TRUE
            }
        }
        mat <- if (methods::is(obj, "SummarizedExperiment")) {
            SummarizedExperiment::assay(obj, 1)
        } else {
            as.matrix(obj)
        }
        
        q_val <- as.numeric(sub("^q_", "", q_name))
        combined_assays_dict[[q_name]] <- list(matrix = mat, q_val = q_val, se_obj = obj)
    }

    if (is.null(first_se)) {
        stop("No valid SummarizedExperiment found in analysis@diversity_results")
    }

    list(
        objects = combined_assays_dict,
        q_names = names(combined_assays_dict),
        first_se = first_se,
        bootstrap_ci_available = bootstrap_ci_available
    )
}

#' Normalize matrix dimensions and row order
#'
#' @param matrix Matrix to normalize
#' @param target_genes Target gene order (character vector)
#' @param target_n_cols Target number of columns
#' @return Normalized matrix

#' @noRd
.normalize_matrix_to_target <- function(matrix, target_genes, target_n_cols) {
    # Adjust column count
    if (ncol(matrix) != target_n_cols) {
        if (ncol(matrix) > target_n_cols) {
            matrix <- matrix[, seq_len(target_n_cols), drop = FALSE]
        } else {
            pad_cols <- target_n_cols - ncol(matrix)
            matrix <- cbind(matrix, matrix(0, nrow = nrow(matrix), ncol = pad_cols))
        }
    }
    
    # Reorder rows to match target genes
    matrix[target_genes, , drop = FALSE]
}

#' Extract bootstrap CI matrices from SE object
#'
#' @param se_obj SummarizedExperiment or matrix object
#' @param target_genes Target gene order
#' @param target_n_cols Target number of columns
#' @param assay_names Assay names in SE
#' @return List(ci_lower, ci_upper) or NULL

#' @noRd
.extract_bootstrap_ci_matrices <- function(se_obj, target_genes, target_n_cols, assay_names) {
    require_pkgs("SummarizedExperiment")
    
    if (!methods::is(se_obj, "SummarizedExperiment")) {
        return(NULL)
    }
    
    if (!("ci_lower" %in% assay_names && "ci_upper" %in% assay_names)) {
        return(NULL)
    }
    
    # Extract ci_lower
    ci_lower <- SummarizedExperiment::assay(se_obj, "ci_lower")
    ci_lower <- .normalize_matrix_to_target(ci_lower, target_genes, target_n_cols)
    
    # Extract ci_upper
    ci_upper <- SummarizedExperiment::assay(se_obj, "ci_upper")
    ci_upper <- .normalize_matrix_to_target(ci_upper, target_genes, target_n_cols)
    
    list(ci_lower = ci_lower, ci_upper = ci_upper)
}

#' Create Q-value suffixed column names
#'
#' @param colnames Column names (character vector or NULL)
#' @param q_val Q-value (numeric)
#' @param n_cols Number of column names needed
#' @return Character vector with _q=X.XXX suffix

#' @noRd
.create_q_suffixed_colnames <- function(colnames, q_val, n_cols) {
    if (is.null(colnames) || length(colnames) == 0) {
        colnames <- paste0("sample_", seq_len(n_cols))
    }
    
    clean_colnames <- sub("_q=.*$", "", colnames)
    paste0(clean_colnames, "_q=", formatC(q_val, format = "f", digits = 3))
}

#' Build combined colData across all q-values
#'
#' @param div_list Original diversity results list
#' @param q_names Q-value names (keys from div_list)
#' @param unique_colnames Final combined column names with q-suffix
#' @return Data frame with combined colData

#' @noRd
.build_combined_coldata <- function(div_list, q_names, unique_colnames_list) {
    require_pkgs("SummarizedExperiment")
    
    combined_coldata_list <- list()
    
    for (q_name in q_names) {
        q_val <- as.numeric(sub("^q_", "", q_name))
        
        # Access unique_colnames by q-value name (stored as list keys in .fill_combined_assays)
        unique_colnames <- unique_colnames_list[[q_name]]
        
        if (methods::is(div_list[[q_name]], "SummarizedExperiment")) {
            cd <- as.data.frame(SummarizedExperiment::colData(div_list[[q_name]]))
        } else {
            cd <- data.frame(row.names = unique_colnames)
        }
        
        cd$q <- q_val
        rownames(cd) <- unique_colnames
        combined_coldata_list[[q_name]] <- cd
    }
    
    do.call(rbind, combined_coldata_list)
}

#' Create combined SummarizedExperiment with assays and metadata
#'
#' @param combined_assay Main diversity assay matrix
#' @param combined_ci_lower CI lower matrix (optional)
#' @param combined_ci_upper CI upper matrix (optional)
#' @param combined_coldata ColData frame
#' @param first_se Template SE for rowData
#' @return SummarizedExperiment object

#' @noRd
.create_combined_se_object <- function(combined_assay, combined_ci_lower, combined_ci_upper,
                                        combined_coldata, first_se) {
    require_pkgs(c("SummarizedExperiment", "S4Vectors"))
    
    # Extract or create rowData, ensuring dimensions match combined_assay
    rd_combined <- tryCatch({
        rd_temp <- SummarizedExperiment::rowData(first_se)
        if (!is.null(rd_temp) && nrow(rd_temp) == nrow(combined_assay)) {
            # Ensure rownames match combined_assay
            rownames(rd_temp) <- rownames(combined_assay)
            rd_temp
        } else {
            NULL
        }
    }, error = function(e) NULL)
    
    if (is.null(rd_combined) || nrow(rd_combined) != nrow(combined_assay)) {
        rd_combined <- data.frame(
            gene_id = rownames(combined_assay),
            row.names = rownames(combined_assay),
            stringsAsFactors = FALSE
        )
    } else {
        # Ensure rownames match even if we're using extracted rowData
        rownames(rd_combined) <- rownames(combined_assay)
    }
    
    # Validate dimensions
    if (ncol(combined_assay) != nrow(combined_coldata)) {
        stop("Column mismatch: assay has ", ncol(combined_assay),
            " columns but colData has ", nrow(combined_coldata), " rows")
    }
    if (nrow(combined_assay) != nrow(rd_combined)) {
        stop("Row mismatch: assay has ", nrow(combined_assay),
            " rows but rowData has ", nrow(rd_combined), " rows")
    }
    
    # Validate names match
    if (!identical(colnames(combined_assay), rownames(combined_coldata))) {
        stop("Column name mismatch between assay and colData")
    }
    if (!identical(rownames(combined_assay), rownames(rd_combined))) {
        stop("Row name mismatch between assay and rowData")
    }
    
    # Build assays list
    assays_list <- list(diversity = combined_assay)
    
    if (!is.null(combined_ci_lower) && !is.null(combined_ci_upper)) {
        ci_lower_valid <- sum(!is.na(combined_ci_lower)) > 0
        ci_upper_valid <- sum(!is.na(combined_ci_upper)) > 0
        
        if (ci_lower_valid && ci_upper_valid) {
            assays_list$ci_lower <- combined_ci_lower
            assays_list$ci_upper <- combined_ci_upper
        }
    }
    
    # Create SE
    combined_se <- SummarizedExperiment::SummarizedExperiment(
        assays = assays_list,
        colData = combined_coldata,
        rowData = rd_combined
    )
    
    # Add metadata if CI available
    if (!is.null(combined_ci_lower)) {
        S4Vectors::metadata(combined_se)$bootstrap_ci_count <- sum(!is.na(combined_ci_lower))
        S4Vectors::metadata(combined_se)$has_bootstrap_ci <- (sum(!is.na(combined_ci_lower)) > 0)
    }
    
    combined_se
}

#' Prepare single q-value data for combined assays
#'
#' @param q_name Q-value name (key from combined_assays_dict)
#' @param combined_assays_dict Dictionary of matrices and metadata
#' @param target_genes Target gene order
#' @param target_n_cols Target columns per q-value
#' @param bootstrap_ci_available Boolean: CIs available
#' @return List with: unique_colnames, ncol, ci_lower, ci_upper

#' @noRd
.prepare_q_value_for_combining <- function(q_name, combined_assays_dict, target_genes,
                                            target_n_cols, bootstrap_ci_available) {
    require_pkgs("SummarizedExperiment")
    
    mat <- combined_assays_dict[[q_name]]$matrix
    q_val <- combined_assays_dict[[q_name]]$q_val
    se_obj <- combined_assays_dict[[q_name]]$se_obj
    
    # Normalize matrix dimensions
    mat <- .normalize_matrix_to_target(mat, target_genes, target_n_cols)
    
    # Create q-suffixed column names
    unique_colnames <- .create_q_suffixed_colnames(colnames(mat), q_val, ncol(mat))
    
    # Extract CI matrices if available
    ci_lower <- ci_upper <- NULL
    if (bootstrap_ci_available) {
        sim_names <- SummarizedExperiment::assayNames(se_obj)
        ci_matrices <- .extract_bootstrap_ci_matrices(
            se_obj, target_genes, target_n_cols, sim_names
        )
        if (!is.null(ci_matrices)) {
            ci_lower <- ci_matrices$ci_lower
            ci_upper <- ci_matrices$ci_upper
        } else {
            ci_lower <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            ci_upper <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
        }
    }
    
    list(
        matrix = mat,
        unique_colnames = unique_colnames,
        ncol_val = ncol(mat),
        ci_lower = ci_lower,
        ci_upper = ci_upper
    )
}

#' Fill combined assay matrices with data from all q-values
#'
#' @param combined_assays_dict Dictionary of matrices and metadata per q-value
#' @param q_names Q-value names (keys)
#' @param target_genes Target gene order
#' @param target_n_cols Target number of columns
#' @param bootstrap_ci_available Boolean: CIs available
#' @return List with: combined_assay, combined_ci_lower, combined_ci_upper, unique_colnames_list

#' @noRd
.fill_combined_assays <- function(combined_assays_dict, q_names, target_genes,
                                   target_n_cols, bootstrap_ci_available) {
    total_cols <- target_n_cols * length(q_names)
    
    # Initialize matrices
    combined_assay <- matrix(0, nrow = length(target_genes), ncol = total_cols)
    rownames(combined_assay) <- target_genes
    
    combined_ci_lower <- if (bootstrap_ci_available) {
        matrix(NA, nrow = length(target_genes), ncol = total_cols)
    } else NULL
    combined_ci_upper <- if (bootstrap_ci_available) {
        matrix(NA, nrow = length(target_genes), ncol = total_cols)
    } else NULL
    
    unique_colnames_list <- list()
    col_idx <- 1
    
    for (q_name in q_names) {
        result <- .prepare_q_value_for_combining(
            q_name, combined_assays_dict, target_genes,
            target_n_cols, bootstrap_ci_available
        )
        
        ncol_q <- result$ncol_val
        if (col_idx + ncol_q - 1 > total_cols) {
            stop("Dimension mismatch: ", col_idx, " to ", col_idx + ncol_q - 1,
                " exceeds total_cols=", total_cols)
        }
        
        # Fill main assay
        combined_assay[, col_idx:(col_idx + ncol_q - 1)] <- result$matrix
        unique_colnames_list[[q_name]] <- result$unique_colnames
        
        # Fill CI matrices if available
        if (bootstrap_ci_available && !is.null(result$ci_lower)) {
            combined_ci_lower[, col_idx:(col_idx + ncol_q - 1)] <- result$ci_lower
            combined_ci_upper[, col_idx:(col_idx + ncol_q - 1)] <- result$ci_upper
        }
        
        col_idx <- col_idx + ncol_q
    }
    
    list(
        combined_assay = combined_assay,
        combined_ci_lower = combined_ci_lower,
        combined_ci_upper = combined_ci_upper,
        unique_colnames_list = unique_colnames_list
    )
}

#' Convert TSENATAnalysis to combined SummarizedExperiment
#'
#' @param analysis TSENATAnalysis object with diversity_results
#' @return SummarizedExperiment with combined assay across all q-values

#' @noRd
.prepare_combined_se <- function(analysis) {
    require_pkgs(c("SummarizedExperiment", "S4Vectors"))

    div_list <- analysis@diversity_results
    
    # Step 1: Extract diversity objects and metadata
    extracted <- .extract_diversity_objects(div_list)
    
    # Step 2: Get target dimensions
    target_genes <- rownames(extracted$first_se)
    target_n_cols <- ncol(extracted$first_se)
    
    # Step 3: Fill combined assays
    filled <- .fill_combined_assays(
        extracted$objects, extracted$q_names, target_genes,
        target_n_cols, extracted$bootstrap_ci_available
    )
    
    # Step 4: Build combined colData (which defines the sample names via rownames)
    combined_coldata_df <- .build_combined_coldata(
        div_list, extracted$q_names, filled$unique_colnames_list
    )
    
    # Step 5: Set column names on all assays to match colData rownames
    combined_colnames <- rownames(combined_coldata_df)
    colnames(filled$combined_assay) <- combined_colnames
    if (!is.null(filled$combined_ci_lower) && !is.null(filled$combined_ci_upper)) {
        colnames(filled$combined_ci_lower) <- combined_colnames
        colnames(filled$combined_ci_upper) <- combined_colnames
    }
    
    # Step 6: Create and return combined SE
    .create_combined_se_object(
        filled$combined_assay, filled$combined_ci_lower, filled$combined_ci_upper,
        combined_coldata_df, extracted$first_se
    )
}

#' Compute gene-level statistics (median +/- SD) by group and q-value
#'
#' @param long_data Long-format data frame with Gene, q, group, tsallis columns
#' @return Data frame with central tendency and spread by gene, group, q

#' @noRd
.compute_gene_group_stats <- function(long_data, metric = "iqr") {
    require_pkgs("dplyr")

    metric <- match.arg(tolower(metric), c("iqr", "sd"))
    long_data$qnum <- as.numeric(as.character(long_data$q))

    # Calculate spread based on metric choice
    if (metric == "iqr") {
        spread_calc <- quote(stats::IQR(tsallis, na.rm = TRUE)/2)
    } else {
        spread_calc <- quote(sqrt(stats::var(tsallis, na.rm = TRUE)))
    }

    dplyr::summarise(dplyr::group_by(long_data, group, qnum), central = median(tsallis,
        na.rm = TRUE), spread = !!spread_calc, .groups = "drop")
}

#' Aggregate bootstrap CI bounds by group and q-value
#'
# NOTE (March 2026): .bootstrap_aggregate_ci() moved to bootstrap.R for
# consolidation

# ============================================================================
# GAM INTERACTION HELPERS
# ============================================================================

#' Build sample-to-group mapping from colData
#'
#' @param cdata SummarizedExperiment colData with sample metadata
#' @param condition_col Column name for group assignments
#' @return Named character vector: sample name -> group value

#' @noRd
.prepare_sample_group_mapping <- function(cdata, condition_col) {
    coldata_rownames <- rownames(cdata)
    coldata_sample_names <- sub("_q=.*", "", coldata_rownames)

    unique_samples <- unique(coldata_sample_names)
    sample_to_group <- character(length(unique_samples))
    names(sample_to_group) <- unique_samples

    for (samp in unique_samples) {
        idx <- which(coldata_sample_names == samp)[1]
        sample_to_group[samp] <- as.character(cdata[[condition_col]][idx])
    }

    sample_to_group
}

#' Prepare long-format plot data for a single gene
#'
#' @param gene Gene ID to extract
#' @param mat Assay matrix (genes x samples*q)
#' @param sample_to_group Named vector mapping sample names to groups
#' @return Data frame with columns: sample, group, q, entropy (or NULL if
#' invalid)

#' @noRd
.plot_gam_prepare_gene_data <- function(gene, mat, sample_to_group) {
    if (!(gene %in% rownames(mat))) {
        return(NULL)
    }

    gene_vals <- mat[gene, ]
    col_names_full <- colnames(mat)

    # Parse column names: 'Sample_q=value'
    col_sample_names <- sub("_q=.*", "", col_names_full)
    col_q_values <- as.numeric(sub(".*_q=", "", col_names_full))

    # Look up group for each column
    col_groups <- unname(sample_to_group[col_sample_names])

    if (any(is.na(col_groups))) {
        return(NULL)
    }

    # Build long-format data frame
    plot_df <- data.frame(sample = col_sample_names, group = col_groups, q = col_q_values,
        entropy = as.numeric(gene_vals), stringsAsFactors = FALSE)

    # Remove NA entries
    plot_df <- plot_df[!is.na(plot_df$entropy), , drop = FALSE]

    if (nrow(plot_df) == 0) {
        return(NULL)
    }

    plot_df
}

#' Fit GAM models per group and generate predictions
#'
#' @param plot_df Long-format data frame with sample, group, q, entropy
#' @return List with $plot_data and $pred_data data frames (or NULL if
#' fitting fails)

#' @noRd
.plot_gam_fit_group <- function(plot_df) {
    require_pkgs(c("mgcv", "dplyr"))

    unique_groups <- unique(plot_df$group)

    if (length(unique_groups) < 2) {
        return(NULL)
    }

    # Generate prediction grid
    q_range <- range(plot_df$q, na.rm = TRUE)
    if (!is.finite(q_range[1]) || !is.finite(q_range[2])) {
        return(NULL)
    }

    pred_q <- seq(q_range[1], q_range[2], length.out = 100)

    # Fit GAM and predict for each group
    pred_list <- list()
    for (gr in unique_groups) {
        subset_data <- subset(plot_df, group == gr)

        if (nrow(subset_data) < 3) {
            next
        }

        tryCatch({
            k <- min(10, max(2, round(nrow(subset_data)/2)))
            gam_fit <- mgcv::gam(entropy ~ s(q, k = k), data = subset_data)

            pred_data <- data.frame(q = pred_q)
            pred_vals <- stats::predict(gam_fit, newdata = pred_data, se.fit = TRUE)

            pred_list[[as.character(gr)]] <- data.frame(group = gr, q = pred_q, entropy_fit = pred_vals$fit,
                se = pred_vals$se.fit, stringsAsFactors = FALSE)
        }, error = function(e) {
            # Silently skip failed fits
            NULL
        })
    }

    if (length(pred_list) == 0) {
        return(NULL)
    }

    pred_df <- do.call(rbind, pred_list)

    # Ensure group is factor with consistent levels
    group_levels <- sort(unique(c(as.character(plot_df$group), as.character(pred_df$group))))
    plot_df$group <- factor(plot_df$group, levels = group_levels)
    pred_df$group <- factor(pred_df$group, levels = group_levels)

    list(plot_data = plot_df, pred_data = pred_df, group_levels = group_levels)
}

#' Select genes to plot based on significance
#'
#' @param lm_res Data frame with gene and p-value columns
#' @param genes Optional character vector of specific genes
#' @param n_top Number of top genes to select
#' @param sig_alpha Significance threshold
#' @return Character vector of gene IDs to plot (or NULL if none selected)

#' @noRd
.plot_select_genes <- function(lm_res, genes = NULL, n_top = 6, sig_alpha = 0.05) {
    if (!is.null(genes)) {
        if (!is.character(genes)) {
            stop("genes must be a character vector of gene names", call. = FALSE)
        }
        return(genes)
    }

    # Select top genes by p-value
    if ("adj_p_interaction" %in% colnames(lm_res)) {
        sig_mask <- lm_res$adj_p_interaction <= sig_alpha
    } else if ("p_interaction" %in% colnames(lm_res)) {
        sig_mask <- lm_res$p_interaction <= sig_alpha
    } else {
        stop("lm_res must contain 'adj_p_interaction' or 'p_interaction' column",
            call. = FALSE)
    }

    sig_genes <- lm_res[sig_mask, , drop = FALSE]

    if (nrow(sig_genes) == 0) {
        return(NULL)
    }

    # Select top n
    sig_genes$gene[seq_len(min(n_top, nrow(sig_genes)))]
}

#' Prepare gene CI data for plotting
#'
#' Converts CI matrices into long-format data for gene-specific bootstrap CI
#' plotting.
#'
#' @param long_data Long-format data with Gene, group, q, tsallis columns
#' @param ci_lower_mat CI lower bounds matrix (genes x samples*q)
#' @param ci_upper_mat CI upper bounds matrix (genes x samples*q)
#' @param genes Character vector of gene IDs to extract
#'
#' @return Data frame with columns: Gene, group, q, median, ci_lower, ci_upper
#'
#' @noRd
.prepare_gene_ci_data <- function(long_data, ci_lower_mat, ci_upper_mat, genes) {
    require_pkgs("dplyr")

    # Aggregate to get median per gene, group, q
    stats_df <- dplyr::summarise(dplyr::group_by(long_data, Gene, group, q), median = median(tsallis,
        na.rm = TRUE), .groups = "drop")

    # Extract CI values for each gene, group, q combination
    plot_df <- stats_df
    plot_df$ci_lower <- NA_real_
    plot_df$ci_upper <- NA_real_

    # Map CI assay columns to gene/group/q combinations
    if (nrow(ci_lower_mat) > 0) {
        colnames_ci <- colnames(ci_lower_mat)

        # Parse column names (e.g., 'Sample_q=0.01')
        ci_samples <- sub("_q=.*", "", colnames_ci)
        ci_q_values <- as.numeric(sub(".*_q=", "", colnames_ci))

        for (i in seq_len(nrow(plot_df))) {
            g <- plot_df$Gene[i]
            gr <- as.character(plot_df$group[i])
            q_val <- plot_df$q[i]

            # Find indices in long_data for this gene/group/q
            matching_rows <- which(as.character(long_data$Gene) == g & as.character(long_data$group) ==
                gr & as.numeric(as.character(long_data$q)) == q_val)

            if (length(matching_rows) > 0) {
                # Get samples for this group from long_data
                samples_for_group <- unique(as.character(long_data$sample[matching_rows]))

                # Find CI columns for these samples at this q
                ci_col_mask <- (ci_samples %in% samples_for_group) & (abs(ci_q_values -
                  q_val) < 1e-06)
                ci_col_indices <- which(ci_col_mask)

                if (length(ci_col_indices) > 0) {
                  # Get CI bounds for these columns
                  gene_idx <- which(rownames(ci_lower_mat) == g)
                  if (length(gene_idx) > 0) {
                    ci_lower_vals <- ci_lower_mat[gene_idx, ci_col_indices]
                    ci_upper_vals <- ci_upper_mat[gene_idx, ci_col_indices]

                    # Use median of CI values across samples in this group
                    plot_df$ci_lower[i] <- median(ci_lower_vals, na.rm = TRUE)
                    plot_df$ci_upper[i] <- median(ci_upper_vals, na.rm = TRUE)
                  }
                }
            }
        }
    }

    plot_df
}

# ============================================================================
# GAM PLOT HELPERS (for main function refactoring)
# ============================================================================

#' Handle and Validate Inputs for GAM Interaction Plot
#'
#' Validates SE object and handles flexible lm_res input formats.
#' Extracts results dataframe from list or validates dataframe directly.
#'
#' @param se A \code{SummarizedExperiment} object
#' @param lm_res Either a data.frame with 'gene' column or list with
#'   $results and $model_data
#'
#' @return List with validated components:
#'   - se: validated SummarizedExperiment
#'   - lm_res: extracted results dataframe
#'   - model_data: extracted model_data (or NULL)
#'
#' @noRd
.plot_gam_handle_inputs <- function(se, lm_res) {
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }

    model_data <- NULL

    # Handle flexible input: lm_res can be either: 1. A data.frame with
    # results (traditional usage) 2. A list with $results and $model_data
    # (return_model_data = TRUE format)
    if (is.list(lm_res) && !is.data.frame(lm_res)) {
        # lm_res is a list with components
        if ("results" %in% names(lm_res) && is.data.frame(lm_res$results)) {
            # Extract results and model_data from the list
            extracted_results <- lm_res$results

            # If model_data provided, extract from lm_res
            if ("model_data" %in% names(lm_res)) {
                model_data <- lm_res$model_data
            }

            lm_res <- extracted_results
        } else {
            stop("lm_res is a list but does not contain 'results' data.frame component",
                call. = FALSE)
        }
    }

    if (!is.data.frame(lm_res) || !("gene" %in% colnames(lm_res))) {
        stop("lm_res must be either:\n  1. A data.frame with 'gene' column from .calculate_lm_interaction()\n  2. A list with $results and $model_data from return_model_data = TRUE",
            call. = FALSE)
    }

    if (nrow(lm_res) == 0) {
        stop("lm_res has no rows; .calculate_lm_interaction() returned no genes",
            call. = FALSE)
    }

    list(se = se, lm_res = lm_res, model_data = model_data)
}

#' Validate and Extract Q-Values from Model Data
#'
#' Validates model_data and extracts/normalizes q-values for GAM analysis.
#'
#' @param model_data List from .calculate_lm_interaction(...,
#' return_model_data = TRUE)$model_data
#'
#' @return Numeric vector of q-values
#'
#' @noRd
.plot_gam_extract_q_values <- function(model_data) {
    if (is.null(model_data)) {
        stop("model_data is required. Provide it as a parameter or pass full lm_res list with $model_data component",
            call. = FALSE)
    }

    if (!is.list(model_data)) {
        stop("model_data must be a list from .calculate_lm_interaction(..., return_model_data = TRUE)",
            call. = FALSE)
    }

    # Extract required metadata - handle both wrapped (from JSON) and unwrapped
    # formats JSON returns arrays: method = [['gam']], need [[1]] Direct list
    # returns: method = 'gam', no [[1]] needed
    q_values <- model_data$q_values
    if (is.null(q_values)) {
        stop("model_data must contain 'q_values' from the original analysis", call. = FALSE)
    }

    # Normalize q_values in case it's wrapped in list
    if (is.list(q_values) && length(q_values) == 1) {
        q_values <- unlist(q_values)
    } else {
        q_values <- unlist(q_values)
    }

    q_values
}

#' Match and Filter Genes Between SE and Results
#'
#' Finds genes present in both SE rownames and lm_res results.
#' Subsets both objects to matching genes only.
#'
#' @param se A \code{SummarizedExperiment}
#' @param lm_res Results dataframe with 'gene' column
#'
#' @return List with:
#'   - se: subset SE
#'   - lm_res: subset results
#'
#' @noRd
.plot_gam_match_genes <- function(se, lm_res) {
    # Match and filter genes between SE and lm_res After calculate_diversity,
    # rownames(SE) are gene names lm_res$gene column also contains gene names
    gene_names_in_results <- lm_res$gene
    gene_names_in_se <- rownames(se)

    # Find genes that exist in both
    available_genes <- gene_names_in_se[gene_names_in_se %in% gene_names_in_results]

    if (length(available_genes) == 0) {
        stop(sprintf("No genes from lm_res found in rownames(se). \n  Examples from lm_res: %s\n  Examples from SE: %s",
            paste(head(gene_names_in_results, 3), collapse = ", "), paste(head(gene_names_in_se,
                3), collapse = ", ")), call. = FALSE)
    }

    # Subset SE to only genes that are in results
    se <- se[available_genes, ]

    # Subset results to only genes that are in SE
    lm_res <- lm_res[lm_res$gene %in% rownames(se), ]

    list(se = se, lm_res = lm_res)
}

#' Create Gene ID to Display Name Mapping
#'
#' Builds mapping for display names from lm_res. Uses gene_name column if
#' available, otherwise uses gene IDs.
#'
#' @param lm_res Results dataframe with 'gene' column and optional
#' 'gene_name' column
#'
#' @return Named character vector mapping gene IDs to display names
#'
#' @noRd
.plot_gam_create_gene_map <- function(lm_res) {
    # Create gene ID to display name mapping from lm_res
    gene_name_map <- setNames(lm_res$gene, lm_res$gene)  # default: use gene ID

    # If lm_res has a gene_name column (e.g., from return_model_data), use it
    if ("gene_name" %in% colnames(lm_res)) {
        gene_name_map <- setNames(lm_res$gene_name, lm_res$gene)
    }

    gene_name_map
}

#' Arrange GAM Plot Grid with Title and Legend
#'
#' Creates final gridded layout of plots with title, plots, and legend.
#'
#' @param plots List of ggplot objects for each gene
#' @param condition_col Column name used for legend title
#' @param font_sizes List with legend font sizes (legend_title, legend_text)
#'
#' @return Single combined ggplot object
#'
#' @noRd
.plot_gam_arrange_grid <- function(plots, condition_col, font_sizes) {
    require_pkgs("cowplot")

    n_plots <- length(plots)
    n_cols <- 2
    n_rows <- ceiling(n_plots / n_cols)

    # Add margins to plots for spacing, particularly between rows
    plots_with_margins <- lapply(seq_along(plots), function(i) {
        p <- plots[[i]]
        # Add larger bottom margin for plots in the first row to create space
        # before second row
        if (i <= n_cols) {
            p <- p + ggplot2::theme(plot.margin = ggplot2::margin(b = 15, unit = "pt"))
        }
        p
    })

    # Extract legend from first plot
    legend <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom",
        legend.title = ggplot2::element_text(size = font_sizes$legend_title), legend.text = ggplot2::element_text(size = font_sizes$legend_text)))

    # Create grid without legends
    combined_plot <- cowplot::plot_grid(plotlist = plots_with_margins, nrow = n_rows,
        ncol = n_cols, align = "hv", axis = "lr")

    # Add main title and subtitle above the grid
    title_plot <- cowplot::ggdraw() + cowplot::draw_label("GAM q-curve: Top genes with group interaction",
        fontface = "bold", size = 20, x = 0.5, y = 0.8) + cowplot::draw_label("Fitted smooth curves by group",
        fontface = "italic", size = 16, x = 0.5, y = 0.25, color = "gray40")

    # Combine title, plots, and single legend at bottom
    final_plot <- cowplot::plot_grid(title_plot, combined_plot, legend, nrow = 3,
        rel_heights = c(0.12, 1, 0.08))

    final_plot
}

#' Save Plot to File
#'
#' Saves a ggplot object to file with optional dimensions and adaptive font sizing.
#'
#' @param plot A ggplot object
#' @param output_file File path for output (if NULL, skips saving)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10.3)
#'
#' @return Invisible NULL; plot saved as side effect
#'
#' @details
#' Uses .calculate_plot_dims() to compute consistent aspect ratios and
#' .scale_font_by_area() to adapt font sizes based on output dimensions.
#'
#' @noRd
.plot_gam_save_plot <- function(plot, output_file, width = NULL, height = NULL) {
    if (!is.null(output_file)) {
        # Use provided dimensions or defaults
        save_width <- if (is.null(width)) 12 else width
        save_height <- if (is.null(height)) 10.3 else height
        
        # Determine aspect type based on provided dimensions
        # tall: height/width ratio > 0.8 (e.g., 10/12 = 0.833)
        # standard: height/width ratio <= 0.8 (e.g., 7.2/12 = 0.6)
        aspect_type <- if (save_height / save_width > 0.8) "tall" else "standard"
        
        # Calculate dimensions via .calculate_plot_dims for consistency
        plot_dims <- .calculate_plot_dims(width_inches = save_width, aspect_type = aspect_type, 
                                         dpi_output = 100)
        
        # Compute adaptive font scaling based on actual area
        font_scale <- .scale_font_by_area(plot_dims$width, plot_dims$height)
        
        # Apply adaptive font scaling if dimensions deviate significantly from reference (96 sq in)
        # Only scale if deviation is >10% to avoid excessive changes
        if (abs(font_scale - 1.0) > 0.1) {
            plot <- plot + ggplot2::theme(
                text = ggplot2::element_text(size = 11 * font_scale),
                plot.title = ggplot2::element_text(size = 14 * font_scale),
                axis.title = ggplot2::element_text(size = 12 * font_scale),
                axis.text = ggplot2::element_text(size = 10 * font_scale),
                legend.text = ggplot2::element_text(size = 10 * font_scale),
                legend.title = ggplot2::element_text(size = 11 * font_scale)
            )
        }
        
        ggplot2::ggsave(output_file, plot = plot, width = plot_dims$width, height = plot_dims$height,
            dpi = plot_dims$dpi, create.dir = TRUE)
    }
    invisible(NULL)
}

#' Create Single GAM Plot for One Gene
#'
#' Generates ggplot object for one gene with fitted GAM curves by group.
#'
#' @param gene Gene ID to plot
#' @param gene_display_name Display name for gene (if NULL, looked up from
#' map)
#' @param gene_name_map Named vector mapping gene IDs to display names
#' @param mat Assay matrix with samples in columns
#' @param sample_to_group Named vector mapping sample names to groups
#' @param condition_col Column name used for legend/facets
#'
#' @return ggplot object or NULL if plot generation fails
#'
#' @noRd
.plot_gam_make_plot <- function(gene, gene_display_name = NULL, gene_name_map,
    mat, sample_to_group, condition_col) {
    require_pkgs(c("ggplot2"))

    # Use provided gene name, or look it up from mapping, or default to gene ID
    if (is.null(gene_display_name)) {
        if (gene %in% names(gene_name_map)) {
            gene_display_name <- gene_name_map[[gene]]
        } else {
            gene_display_name <- gene
        }
    }

    # Prepare plot data for this gene using helper
    plot_df <- .plot_gam_prepare_gene_data(gene, mat, sample_to_group)
    if (is.null(plot_df)) {
        return(NULL)
    }

    # Fit GAM models and generate predictions using helper
    gam_result <- .plot_gam_fit_group(plot_df)
    if (is.null(gam_result)) {
        return(NULL)
    }

    plot_df <- gam_result$plot_data
    pred_df <- gam_result$pred_data
    group_levels <- gam_result$group_levels

    # Create explicit color mapping for all groups
    palette_colors <- .palette_blue_red()
    color_mapping <- c()
    for (i in seq_along(group_levels)) {
        color_idx <- ((i - 1) %% length(palette_colors)) + 1
        color_mapping[group_levels[i]] <- palette_colors[color_idx]
    }

    # Create plot with explicit color scale
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = entropy, color = group)) +
        ggplot2::geom_point(data = plot_df, ggplot2::aes(x = q, y = entropy,
            color = group), alpha = 0.5, size = 2) + ggplot2::geom_line(data = pred_df,
        ggplot2::aes(x = q, y = entropy_fit, color = group, linetype = "GAM fit"),
        linewidth = 1, alpha = 0.9) + ggplot2::scale_color_manual(values = color_mapping,
        name = condition_col, breaks = group_levels) + ggplot2::scale_linetype_manual(values = c(`GAM fit` = 1),
        name = "") + ggplot2::labs(x = "q parameter", y = "Tsallis entropy",
        title = ifelse(gene_display_name != gene, sprintf("%s (%s)", gene_display_name,
            gene), gene_display_name)) + .theme_spectrum(base_size = 11) + ggplot2::theme(legend.position = "none")

    p
}

# ============================================================================
# PLOT COMPOSITION MEGA-HELPERS (Consolidation Phase)
# ============================================================================

#' Apply Publication-Ready Theme with Centered Titles
#'
#' Consolidated helper that applies base theme + centered/bolded title/subtitle.
#' Replaces repeated 35+ line pattern across all plot files.
#'
#' AESTHETIC PRESERVATION: Uses exact same font sizes and styling as originals.
#' - Title: .font_sizes$title, bold, centered
#' - Subtitle: .font_sizes$subtitle, italic, centered  
#' - Base theme: .theme_base (or .theme_spectrum for spectrum plots)
#'
#' @param plot ggplot2 object to style
#' @param title Character: plot title (optional)
#' @param subtitle Character: plot subtitle (optional)
#' @param base_theme Character: "theme_base" (default) or "theme_spectrum"
#' @param base_size Integer: base font size (default: 11, matches .theme_base default)
#' @param title_size Integer: title font size (default: from .font_sizes constants)
#' @param subtitle_size Integer: subtitle font size (default: from .font_sizes constants)
#'
#' @return Modified ggplot2 object with applied theme
#'
#' @noRd
.apply_publication_theme <- function(plot, title = NULL, subtitle = NULL,
                                     base_theme = "theme_base", base_size = 11,
                                     title_size = .font_sizes$title,
                                     subtitle_size = .font_sizes$subtitle) {
    require_pkgs("ggplot2")
    
    # Apply base theme (either .theme_base or .theme_spectrum)
    # Add dot prefix if not already present
    theme_name <- if (startsWith(base_theme, ".")) base_theme else paste0(".", base_theme)
    theme_fn <- get(theme_name)
    result <- plot + theme_fn(base_size = base_size)
    
    # Apply title/subtitle styling (always centered, bold/italic as per TSENAT convention)
    result <- result + ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, 
                                          face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = subtitle_size, 
                                             face = "italic")
    )
    
    # Add title/subtitle labels if provided
    if (!is.null(title) && !is.null(subtitle)) {
        result <- result + ggplot2::labs(title = title, subtitle = subtitle)
    } else if (!is.null(title)) {
        result <- result + ggplot2::labs(title = title)
    } else if (!is.null(subtitle)) {
        result <- result + ggplot2::labs(subtitle = subtitle)
    }
    
    result
}

#' Apply Group Aesthetic Scales (Color + Fill + Legend)
#'
#' Consolidated helper for color palette + manual scales + legend styling.
#' Replaces repeated 21x pattern of palette + scale_color_manual + scale_fill_manual.
#'
#' AESTHETIC PRESERVATION: Uses exact same colors and mappings as originals.
#' - Palette: .palette_blue_red() by default (standard TSENAT convention)
#' - Color/Fill Manual: with name="Group" (standard legend title)
#'
#' @param plot ggplot2 object
#' @param palette Character: palette function name (e.g. "palette_blue_red") OR 
#'   a vector of colors. If character, will call the .palette_* function.
#' @param legend_name Character: legend title (default: "Group")
#' @param legend_position Character: legend position (default: "bottom")
#' @param direction Integer: 1 (normal) or -1 (reversed palette)
#'
#' @return Modified ggplot2 object with applied color scales
#'
#' @noRd
.apply_group_aesthetics <- function(plot, palette = "palette_blue_red",
                                   legend_name = "Group", legend_position = "bottom",
                                   direction = 1) {
    require_pkgs("ggplot2")
    
    # Get palette colors - handle both string (function name) and vector cases
    if (is.character(palette) && length(palette) == 1) {
        # palette is a function name string, call the function
        palette_fn <- get(paste0(".", palette))  # e.g., .palette_blue_red
        colors <- palette_fn()
    } else {
        # palette is already a vector of colors
        colors <- palette
    }
    
    # Reverse if needed
    if (direction == -1) {
        colors <- rev(colors)
    }
    
    # Apply color and fill scales (exact pattern from original code)
    result <- plot +
        ggplot2::scale_color_manual(values = colors, name = legend_name) +
        ggplot2::scale_fill_manual(values = colors, name = legend_name) +
        ggplot2::theme(legend.position = legend_position)
    
    result
}

#' Create Confidence Interval Line Plot Base
#'
#' Consolidated helper for ribbon + line + point layer pattern.
#' Replaces repeated 18x pattern of geom_ribbon + geom_line + geom_point.
#'
#' AESTHETIC PRESERVATION: Uses exact styling from q-curve plots:
#' - Ribbon: alpha=0.15, color=NA (transparent, no outline)
#' - Line: linewidth=1.2
#' - Point: size=3.5, alpha=0.8
#'
#' @param data Data frame with plot data
#' @param x_col Character: name of x column (default: "q")
#' @param y_col Character: name of y column (default: "median")
#' @param group_col Character: name of grouping column (default: "group")
#' @param ci_lower_col Character: name of CI lower column (default: "ci_lower")
#' @param ci_upper_col Character: name of CI upper column (default: "ci_upper")
#' @param ribbon_alpha Numeric: ribbon transparency (default: 0.15)
#' @param line_width Numeric: line width (default: 1.2)
#' @param point_size Numeric: point size (default: 3.5)
#' @param show_points Logical: include geom_point layer? (default: TRUE)
#'
#' @return Base ggplot2 object with ribbon/line/point layers
#'
#' @noRd
.create_ci_ribbon_plot <- function(data, x_col = "q", y_col = "median",
                                  group_col = "group",
                                  ci_lower_col = "ci_lower", 
                                  ci_upper_col = "ci_upper",
                                  ribbon_alpha = 0.15, line_width = 1.2, 
                                  point_size = 3.5, show_points = TRUE) {
    require_pkgs("ggplot2")
    
    # Build base plot with aesthetics
    p <- ggplot2::ggplot(data, 
                        ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
                                    color = .data[[group_col]], 
                                    fill = .data[[group_col]]))
    
    # Add ribbon layer (CI bounds)
    p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[ci_lower_col]], 
                    ymax = .data[[ci_upper_col]]),
        alpha = ribbon_alpha, color = NA
    )
    
    # Add line layer
    p <- p + ggplot2::geom_line(linewidth = line_width)
    
    # Add point layer if requested
    if (show_points) {
        p <- p + ggplot2::geom_point(size = point_size, alpha = 0.8)
    }
    
    p
}

#' Assemble Multi-Plot Grid with Shared Legend
#'
#' Consolidated helper for legend extraction + grid composition pattern.
#' Replaces repeated 21x pattern of get_legend + plot_nolegend + plot_grid assembly.
#'
#' AESTHETIC PRESERVATION: Uses exact parameters from original code:
#' - Legend position: "bottom" (default, customizable)
#' - Legend direction: "horizontal" (standard for TSENAT)
#' - Grid alignment: "hv" (both axes aligned)
#' - Title/subtitle: uses .font_sizes constants
#'
#' @param plots List of ggplot2 objects (one per subplot)
#' @param ncol Integer: number of columns (default: 2)
#' @param nrow Integer: number of rows (default: auto-calculated)
#' @param title Character: main title (optional)
#' @param subtitle Character: subtitle under title (optional)
#' @param legend_position Character: position for legend - "bottom", "top", "left", "right", "none" 
#'   (default: "bottom")
#' @param extract_legend Logical: whether to extract and place legend separately 
#'   (default: TRUE). If FALSE, plots retain their individual legends.
#' @param rel_heights Numeric vector: relative heights for title/plots/legend
#'   (default: c(0.08, 1, 0.08) - 8% title, 100% plots, 8% legend)
#'
#' @return Combined ggplot2/cowplot object ready for display/saving
#'
#' @details
#' Automatically:
#' - Extracts legend from first plot (if extract_legend=TRUE)
#' - Removes legends from all individual plots (if extract_legend=TRUE)
#' - Arranges in grid with specified layout
#' - Adds optional title/subtitle at top
#' - Places shared legend at specified position
#'
#' @noRd
.assemble_grid_plot <- function(plots, ncol = 2, nrow = NULL, 
                               title = NULL, subtitle = NULL,
                               legend_position = "bottom",
                               extract_legend = TRUE,
                               rel_heights = c(0.08, 1, 0.08)) {
    require_pkgs(c("ggplot2", "cowplot"))
    
    if (length(plots) == 0) {
        stop("plots list cannot be empty", call. = FALSE)
    }
    
    # Calculate nrow if not provided
    if (is.null(nrow)) {
        nrow <- ceiling(length(plots) / ncol)
    }
    
    # Extract legend from first plot if requested
    legend_obj <- NULL
    if (extract_legend) {
        legend_obj <- cowplot::get_legend(
            plots[[1]] + ggplot2::theme(legend.position = legend_position,
                                       legend.direction = "horizontal")
        )
    }
    
    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(p) {
        p + ggplot2::theme(legend.position = "none")
    })
    
    # Compose grid without legend
    grid_plot <- cowplot::plot_grid(plotlist = plots_nolegend, 
                                   ncol = ncol, nrow = nrow,
                                   align = "hv")
    
    # If title/subtitle provided, create title grobs and assemble all 3 components
    if (!is.null(title) || !is.null(subtitle)) {
        title_plot <- cowplot::ggdraw()
        
        if (!is.null(title)) {
            title_plot <- title_plot + 
                cowplot::draw_label(title, fontface = "bold", size = .font_sizes$title,
                                  x = 0.5, hjust = 0.5)
        }
        
        if (!is.null(subtitle)) {
            y_pos <- if (is.null(title)) 0.5 else 0.25
            subtitle_plot <- cowplot::ggdraw() + 
                cowplot::draw_label(subtitle, fontface = "italic", 
                                  size = .font_sizes$subtitle,
                                  x = 0.5, hjust = 0.5, color = "gray40")
            title_plot <- cowplot::plot_grid(title_plot, subtitle_plot, 
                                           nrow = 2, rel_heights = c(1, 0.6))
        }
        
        # Assemble title + grid + legend (if extracted)
        if (extract_legend && !is.null(legend_obj)) {
            return(cowplot::plot_grid(title_plot, grid_plot, legend_obj, 
                                     nrow = 3, rel_heights = rel_heights))
        } else {
            # No legend: just title + grid
            return(cowplot::plot_grid(title_plot, grid_plot, 
                                     nrow = 2, rel_heights = rel_heights[c(1, 2)]))
        }
    }
    
    # No title/subtitle: just grid + legend (if extracted)
    if (extract_legend && !is.null(legend_obj)) {
        return(cowplot::plot_grid(grid_plot, legend_obj, nrow = 2, 
                          rel_heights = rel_heights[c(2, 3)]))
    } else {
        return(grid_plot)
    }
}