# ============================================================================
# PLOT EXPRESSION HELPERS
# ============================================================================

#' Prepare inputs for expression plot
#'
#' Validates and normalizes counts, samples, and tx2gene mapping for transcripts plot.
#'
#' @param counts Matrix or SummarizedExperiment with transcripts as rownames
#' @param readcounts Optional column name for read counts in SE
#' @param samples Character vector of sample group assignments
#' @param coldata Optional data.frame or file path with sample metadata
#' @param condition_col Column in coldata for sample conditions (default: 'condition')
#' @param tx2gene Data frame or file path with tx2gene mapping
#' @param res Optional results data frame for gene selection
#' @param top_n Number of top genes to plot
#' @param pseudocount Pseudocount to add for log transformation
#' @param output_file Optional output file path
#' @param metric Aggregation metric: 'median', 'mean', 'variance', or 'iqr'
#'
#' @return List with normalized counts, samples, mapping, aggregation function
#' @noRd
.make_plot_for_geneprepare_inputs <- function(counts, readcounts = NULL, samples = NULL,
    coldata = NULL, condition_col = "condition", tx2gene = NULL, res = NULL, top_n = NULL,
    pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance",
        "iqr")) {
    # handle selecting genes from `res` is left to caller; this function
    # focuses on normalizing counts, samples and tx2gene mapping and preparing
    # agg functions
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

#' Make plot for a single gene
#'
#' Creates a transcript-level expression plot for a single gene.
#'
#' @param gene_single Character: gene identifier
#' @param mapping Data frame with Transcript and Gen columns
#' @param counts Matrix of read counts
#' @param samples Character vector of sample assignments
#' @param top_n Number of top transcripts to show
#' @param agg_fun Aggregation function for summarization
#' @param pseudocount Pseudocount for log transformation
#' @param agg_label_unique Label for aggregation metric
#' @param fill_limits Optional numeric vector for fill scale limits
#' @param font_scale Font scaling factor
#'
#' @return ggplot2 object
#' @noRd
.make_plot_for_genemake_plot_for_gene <- function(gene_single, mapping, counts, samples,
    top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL, font_scale = 1) {
    built <- .make_plot_for_genebuild_tx_long(gene_single, mapping, counts, samples,
        NULL)
    df_summary <- .make_plot_for_geneaggregate_df_long(built$df_long, agg_fun, pseudocount)
    .make_plot_for_genebuild_plot_from_summary(df_summary, agg_label_unique, fill_limits,
        font_scale = font_scale)
}

#' Combine multiple gene plots
#'
#' Combines individual gene plots into a grid layout.
#'
#' @param plots List of ggplot2 objects (one per gene)
#' @param output_file Optional file path to save combined plot
#' @param agg_label_unique Label for aggregation metric
#'
#' @return Combined plot object or invisible NULL if output_file provided
#' @noRd

.make_plot_for_genecombine_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
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

    # Extract legend from first plot
    p_for_legend <- .configure_legend(plots[[1]], position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        .configure_legend(pp, position = "none")
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

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        .configure_legend(pp, position = "none")
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

    # Handle palette parameter: if length 1 string, check if it's a named palette,
    # otherwise treat as vector of colors
    if (is.character(palette) && length(palette) == 1) {
        if (palette == "blue_red") {
            colors <- .palette_blue_red()
        } else if (palette == "continuous_diverging") {
            colors <- .palette_continuous_diverging()
        } else {
            colors <- .palette_blue_red()
        }
    } else if (is.character(palette)) {
        # palette is a vector of color codes
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

    # Handle palette: check if it's a single string first (safe for equality comparison)
    if (is.character(palette) && length(palette) == 1 && palette == "continuous_diverging") {
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

    list(objects = combined_assays_dict, q_names = names(combined_assays_dict), first_se = first_se,
        bootstrap_ci_available = bootstrap_ci_available)
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

    combined_coldata_list <- list()

    for (q_name in q_names) {
        q_val <- as.numeric(sub("^q_", "", q_name))

        # Access unique_colnames by q-value name (stored as list keys in
        # .fill_combined_assays)
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
        rd_combined <- data.frame(gene_id = rownames(combined_assay), row.names = rownames(combined_assay),
            stringsAsFactors = FALSE)
    } else {
        # Ensure rownames match even if we're using extracted rowData
        rownames(rd_combined) <- rownames(combined_assay)
    }

    # Validate dimensions
    if (ncol(combined_assay) != nrow(combined_coldata)) {
        stop("Column mismatch: assay has ", ncol(combined_assay), " columns but colData has ",
            nrow(combined_coldata), " rows")
    }
    if (nrow(combined_assay) != nrow(rd_combined)) {
        stop("Row mismatch: assay has ", nrow(combined_assay), " rows but rowData has ",
            nrow(rd_combined), " rows")
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
    combined_se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list,
        colData = combined_coldata, rowData = rd_combined)

    # Add metadata if CI available
    if (!is.null(combined_ci_lower)) {
        S4Vectors::metadata(combined_se)$bootstrap_ci_count <- sum(!is.na(combined_ci_lower))
        S4Vectors::metadata(combined_se)$has_bootstrap_ci <- (sum(!is.na(combined_ci_lower)) >
            0)
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
        ci_matrices <- .extract_bootstrap_ci_matrices(se_obj, target_genes, target_n_cols,
            sim_names)
        if (!is.null(ci_matrices)) {
            ci_lower <- ci_matrices$ci_lower
            ci_upper <- ci_matrices$ci_upper
        } else {
            ci_lower <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            ci_upper <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
        }
    }

    list(matrix = mat, unique_colnames = unique_colnames, ncol_val = ncol(mat), ci_lower = ci_lower,
        ci_upper = ci_upper)
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
.fill_combined_assays <- function(combined_assays_dict, q_names, target_genes, target_n_cols,
    bootstrap_ci_available) {
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
        result <- .prepare_q_value_for_combining(q_name, combined_assays_dict, target_genes,
            target_n_cols, bootstrap_ci_available)

        ncol_q <- result$ncol_val
        if (col_idx + ncol_q - 1 > total_cols) {
            stop("Dimension mismatch: ", col_idx, " to ", col_idx + ncol_q - 1, " exceeds total_cols=",
                total_cols)
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

    list(combined_assay = combined_assay, combined_ci_lower = combined_ci_lower,
        combined_ci_upper = combined_ci_upper, unique_colnames_list = unique_colnames_list)
}

#' Convert TSENATAnalysis to combined SummarizedExperiment
#'
#' @param analysis TSENATAnalysis object with diversity_results
#' @return SummarizedExperiment with combined assay across all q-values

#' @noRd
.prepare_combined_se <- function(analysis) {

    div_list <- analysis@diversity_results

    # Step 1: Extract diversity objects and metadata
    extracted <- .extract_diversity_objects(div_list)

    # Step 2: Get target dimensions
    target_genes <- rownames(extracted$first_se)
    target_n_cols <- ncol(extracted$first_se)

    # Step 3: Fill combined assays
    filled <- .fill_combined_assays(extracted$objects, extracted$q_names, target_genes,
        target_n_cols, extracted$bootstrap_ci_available)

    # Step 4: Build combined colData (which defines the sample names via
    # rownames)
    combined_coldata_df <- .build_combined_coldata(div_list, extracted$q_names, filled$unique_colnames_list)

    # Step 5: Set column names on all assays to match colData rownames
    combined_colnames <- rownames(combined_coldata_df)
    colnames(filled$combined_assay) <- combined_colnames
    if (!is.null(filled$combined_ci_lower) && !is.null(filled$combined_ci_upper)) {
        colnames(filled$combined_ci_lower) <- combined_colnames
        colnames(filled$combined_ci_upper) <- combined_colnames
    }

    # Step 6: Create and return combined SE
    .create_combined_se_object(filled$combined_assay, filled$combined_ci_lower, filled$combined_ci_upper,
        combined_coldata_df, extracted$first_se)
}

#' Compute gene-level statistics (median +/- SD) by group and q-value
#'
#' @param long_data Long-format data frame with Gene, q, group, tsallis columns
#' @return Data frame with central tendency and spread by gene, group, q

#' @noRd
.compute_gene_group_stats <- function(long_data, metric = "iqr") {

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

    # Identify p-value column
    if ("adj_p_interaction" %in% colnames(lm_res)) {
        p_col <- "adj_p_interaction"
    } else if ("p_interaction" %in% colnames(lm_res)) {
        p_col <- "p_interaction"
    } else {
        stop("lm_res must contain 'adj_p_interaction' or 'p_interaction' column",
            call. = FALSE)
    }

    # Filter to significant genes (p-value < sig_alpha)
    sig_genes <- lm_res[lm_res[[p_col]] < sig_alpha, , drop = FALSE]

    # If no significant genes, return NULL
    if (nrow(sig_genes) == 0) {
        return(NULL)
    }

    # Sort by p-value and select top n_top genes
    top_idx <- order(sig_genes[[p_col]])[seq_len(min(n_top, nrow(sig_genes)))]
    top_genes <- sig_genes[top_idx, , drop = FALSE]

    # Return gene names sorted by p-value
    top_genes$gene
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

        # Parse column names (e.g., 'sample1_q=1.000')
        ci_samples <- sub("_q=.*", "", colnames_ci)
        ci_q_values <- as.numeric(sub(".*_q=", "", colnames_ci))

        for (i in seq_len(nrow(plot_df))) {
            g <- as.character(plot_df$Gene[i])
            gr <- as.character(plot_df$group[i])
            q_val <- as.numeric(plot_df$q[i])

            # Find indices in long_data for this gene/group/q
            matching_rows <- which(as.character(long_data$Gene) == g & as.character(long_data$group) ==
                gr & abs(as.numeric(as.character(long_data$q)) - q_val) < 1e-06)

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
                    # Use only first match (shouldn't have duplicates but be
                    # safe)
                    gene_idx <- gene_idx[1]
                    ci_lower_vals <- as.numeric(ci_lower_mat[gene_idx, ci_col_indices])
                    ci_upper_vals <- as.numeric(ci_upper_mat[gene_idx, ci_col_indices])

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

    # Handle flexible input: lm_res can be either: 1. A data.frame with results
    # (traditional usage) 2. A list with $results and $model_data
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
        stop("lm_res must be either:\n  1. A data.frame with 'gene' column from .calculate_lm()\n  2. A list with $results and $model_data from return_model_data = TRUE",
            call. = FALSE)
    }

    if (nrow(lm_res) == 0) {
        stop("lm_res has no rows; .calculate_lm() returned no genes",
            call. = FALSE)
    }

    list(se = se, lm_res = lm_res, model_data = model_data)
}

#' Validate and Extract Q-Values from Model Data
#'
#' Validates model_data and extracts/normalizes q-values for GAM analysis.
#'
#' @param model_data List from .calculate_lm(...,
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
        stop("model_data must be a list from .calculate_lm(..., return_model_data = TRUE)",
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

    n_plots <- length(plots)
    n_cols <- 2
    n_rows <- ceiling(n_plots/n_cols)

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
    p_for_legend <- .configure_legend(plots[[1]], position = "bottom", text_size = font_sizes$legend_text,
        title_size = font_sizes$legend_title)
    legend <- cowplot::get_legend(p_for_legend)

    # Create grid without legends
    combined_plot <- cowplot::plot_grid(plotlist = plots_with_margins, nrow = n_rows,
        ncol = n_cols, align = "hv", axis = "lr")

    # Add main title and subtitle above the grid
    title_plot <- .create_title_grob("q-curve: Top genes with group interaction",
        subtitle = "Fitted smooth curves by group", title_size = 20, subtitle_size = 16)

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
        save_width <- if (is.null(width))
            12 else width
        save_height <- if (is.null(height))
            10.3 else height

        # Determine aspect type based on provided dimensions tall: height/width
        # ratio > 0.8 (e.g., 10/12 = 0.833) standard: height/width ratio <= 0.8
        # (e.g., 7.2/12 = 0.6)
        aspect_type <- if (save_height/save_width > 0.8)
            "tall" else "standard"

        # Calculate dimensions via .calculate_plot_dims for consistency
        plot_dims <- .calculate_plot_dims(width_inches = save_width, aspect_type = aspect_type,
            dpi_output = 100)

        # Compute adaptive font scaling based on actual area
        font_scale <- .scale_font_by_area(plot_dims$width, plot_dims$height)

        # Apply adaptive font scaling if dimensions deviate significantly from
        # reference (96 sq in) Only scale if deviation is >10% to avoid
        # excessive changes
        if (abs(font_scale - 1) > 0.1) {
            plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = 11 *
                font_scale), plot.title = ggplot2::element_text(size = 14 * font_scale),
                axis.title = ggplot2::element_text(size = 12 * font_scale), axis.text = ggplot2::element_text(size = 10 *
                  font_scale), legend.text = ggplot2::element_text(size = 10 * font_scale),
                legend.title = ggplot2::element_text(size = 11 * font_scale))
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
.plot_gam_make_plot <- function(gene, gene_display_name = NULL, gene_name_map, mat,
    sample_to_group, condition_col) {

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
        color_idx <- ((i - 1)%%length(palette_colors)) + 1
        color_mapping[group_levels[i]] <- palette_colors[color_idx]
    }

    # Create plot with explicit color scale
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = entropy, color = group)) +
        ggplot2::geom_point(data = plot_df, ggplot2::aes(x = q, y = entropy, color = group),
            alpha = 0.5, size = 2) + ggplot2::geom_line(data = pred_df, ggplot2::aes(x = q,
        y = entropy_fit, color = group, linetype = "Model fit"), linewidth = 1, alpha = 0.9) +
        ggplot2::scale_color_manual(values = color_mapping, name = condition_col,
            breaks = group_levels) + ggplot2::scale_linetype_manual(values = c(`Model fit` = 1),
        name = "") + ggplot2::labs(x = "q parameter", y = "Tsallis entropy", title = ifelse(gene_display_name !=
        gene, sprintf("%s (%s)", gene_display_name, gene), gene_display_name)) +
        .theme_spectrum(base_size = 11)

    p <- .configure_legend(p, position = "none")

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
#' @param base_theme Character: 'theme_base' (default) or 'theme_spectrum'
#' @param base_size Integer: base font size (default: 11, matches .theme_base default)
#' @param title_size Integer: title font size (default: from .font_sizes constants)
#' @param subtitle_size Integer: subtitle font size (default: from .font_sizes constants)
#'
#' @return Modified ggplot2 object with applied theme
#'
#' @noRd
.apply_publication_theme <- function(plot, title = NULL, subtitle = NULL, base_theme = "theme_base",
    base_size = 11, title_size = .font_sizes$title, subtitle_size = .font_sizes$subtitle) {

    # Apply base theme (either .theme_base or .theme_spectrum) Add dot prefix
    # if not already present
    theme_name <- if (startsWith(base_theme, "."))
        base_theme else paste0(".", base_theme)
    theme_fn <- get(theme_name)
    result <- plot + theme_fn(base_size = base_size)

    # Apply title/subtitle styling (always centered, bold/italic as per TSENAT
    # convention)
    result <- result + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = title_size, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = subtitle_size, face = "italic"))

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

# ============================================================================
# PHASE 5: LEGEND & FORMATTING HELPERS
# ============================================================================

#' Configure Legend Positioning, Sizing, and Styling
#'
#' Consolidated helper that standardizes legend appearance across all plots.
#' Replaces repeated 20+ occurrences of legend.position, legend.key.width, 
#' legend.text, etc. customizations throughout the codebase.
#'
#' @param plot ggplot2 object to modify
#' @param position Character: 'bottom', 'right', 'left', 'top', or 'none' (default: 'bottom')
#' @param width_cm Numeric: width of legend key in cm (default: NULL = don't override)
#' @param height_cm Numeric: height of legend key in cm (default: NULL)
#' @param text_size Numeric: font size for legend text (default: NULL = use plot theme)
#' @param title_size Numeric: font size for legend title (default: NULL)
#' @param justification Character: 'left', 'center', 'right' (default: NULL = auto)
#' @param background_color Character: fill color for legend background (default: NULL)
#' @param border_color Character: border color for legend box (default: NULL)
#' @param spacing_lines Numeric: line spacing in legend (default: 2)
#'
#' @return Modified ggplot2 object with configured legend
#'
#' @examples
#' \dontrun{
#' # Simple usage: position and width
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, color = factor(cyl))) +
#'     ggplot2::geom_point()
#' p_leg <- .configure_legend(p, position = 'bottom', width_cm = 2)
#'
#' # With text sizing
#' p_leg <- .configure_legend(p, position = 'right', text_size = 10, title_size = 11)
#'
#' # No legend
#' p_leg <- .configure_legend(p, position = 'none')
#' }
#'
#' @noRd

.configure_legend <- function(plot, position = "bottom", width_cm = NULL, height_cm = NULL,
    text_size = NULL, title_size = NULL, justification = NULL, background_color = NULL,
    border_color = NULL, spacing_lines = 2) {

    theme_list <- list()

    # Handle legend positioning
    if (!is.null(position) && position != "none") {
        theme_list$legend.position <- position
    } else if (position == "none") {
        theme_list$legend.position <- "none"
    }

    # Handle legend justification
    if (!is.null(justification)) {
        theme_list$legend.justification <- justification
    } else if (!is.null(position) && position == "bottom") {
        theme_list$legend.justification <- "center"
    }

    # Handle legend key dimensions
    if (!is.null(width_cm)) {
        theme_list$legend.key.width <- ggplot2::unit(width_cm, "cm")
    }
    if (!is.null(height_cm)) {
        theme_list$legend.key.height <- ggplot2::unit(height_cm, "cm")
    }

    # Handle text sizes
    if (!is.null(text_size)) {
        theme_list$legend.text <- ggplot2::element_text(size = text_size)
    }
    if (!is.null(title_size)) {
        theme_list$legend.title <- ggplot2::element_text(size = title_size, face = "bold")
    }

    # Handle legend background
    if (!is.null(background_color)) {
        theme_list$legend.background <- ggplot2::element_rect(fill = background_color,
            color = border_color %||% "black")
    } else if (!is.null(border_color)) {
        theme_list$legend.background <- ggplot2::element_rect(fill = NA, color = border_color)
    }

    # Handle spacing
    theme_list$legend.spacing.y <- ggplot2::unit(spacing_lines, "mm")

    if (length(theme_list) > 0) {
        plot <- plot + do.call(ggplot2::theme, theme_list)
    }

    return(plot)
}

#' Add Reference Lines (Horizontal and Vertical)
#'
#' Consolidated helper for adding reference/threshold lines to plots.
#' Replaces repeated 12+ occurrences of geom_hline + geom_vline patterns.
#'
#' @param plot ggplot2 object to modify
#' @param h_intercept Numeric vector: y-coordinates for horizontal lines (default: NULL)
#' @param v_intercept Numeric vector: x-coordinates for vertical lines (default: NULL)
#' @param h_color Character: color for horizontal lines (default: 'gray50')
#' @param v_color Character: color for vertical lines (default: 'gray50')
#' @param h_linetype Character: linetype for horizontal lines (default: 'dashed')
#' @param v_linetype Character: linetype for vertical lines (default: 'dashed')
#' @param h_size Numeric: line width for horizontal lines (default: 0.8)
#' @param v_size Numeric: line width for vertical lines (default: 0.8)
#' @param h_alpha Numeric: transparency for horizontal lines (default: 0.7)
#' @param v_alpha Numeric: transparency for vertical lines (default: 0.7)
#'
#' @return Modified ggplot2 object with reference lines
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
#'     ggplot2::geom_point()
#'
#' # Add threshold lines
#' p_ref <- .add_reference_lines(p, h_intercept = 20, v_intercept = 3)
#'
#' # Customized colors
#' p_ref <- .add_reference_lines(p, h_intercept = 20, h_color = 'red',
#'                               v_intercept = c(2.5, 3.5), v_color = 'blue')
#' }
#'
#' @noRd

.add_reference_lines <- function(plot, h_intercept = NULL, v_intercept = NULL, h_color = "gray50",
    v_color = "gray50", h_linetype = "dashed", v_linetype = "dashed", h_size = 0.8,
    v_size = 0.8, h_alpha = 0.7, v_alpha = 0.7) {

    # Add horizontal reference lines
    if (!is.null(h_intercept)) {
        for (yint in h_intercept) {
            plot <- plot + ggplot2::geom_hline(yintercept = yint, color = h_color,
                linetype = h_linetype, linewidth = h_size, alpha = h_alpha)
        }
    }

    # Add vertical reference lines
    if (!is.null(v_intercept)) {
        for (xint in v_intercept) {
            plot <- plot + ggplot2::geom_vline(xintercept = xint, color = v_color,
                linetype = v_linetype, linewidth = v_size, alpha = v_alpha)
        }
    }

    return(plot)
}

#' Format Axis Labels and Titles
#'
#' Consolidated helper for axis label styling including rotation, sizing, and face.
#' Replaces repeated 8+ occurrences of axis.text.x/y + axis.title customizations.
#'
#' @param plot ggplot2 object to modify
#' @param x_angle Numeric: rotation angle for x-axis labels (default: 0)
#' @param y_angle Numeric: rotation angle for y-axis labels (default: 0)
#' @param x_hjust Numeric: horizontal justification for x-axis (default: NULL = auto)
#' @param y_hjust Numeric: horizontal justification for y-axis (default: NULL = auto)
#' @param x_size Numeric: font size for x-axis labels (default: NULL = no override)
#' @param y_size Numeric: font size for y-axis labels (default: NULL = no override)
#' @param x_face Character: font face ('plain', 'bold', 'italic') for x-axis
#' @param y_face Character: font face for y-axis
#' @param x_color Character: text color for x-axis (default: 'black')
#' @param y_color Character: text color for y-axis (default: 'black')
#' @param bold_title Logical: make axis titles bold (default: TRUE)
#'
#' @return Modified ggplot2 object with formatted axis labels
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = factor(cyl), y = mpg)) +
#'     ggplot2::geom_point()
#'
#' # Horizontal x-axis labels at 90 degrees
#' p_fmt <- .format_axis_labels(p, x_angle = 90, x_size = 12, y_size = 11)
#'
#' # Bold axis titles only
#' p_fmt <- .format_axis_labels(p, bold_title = TRUE)
#' }
#'
#' @noRd

.format_axis_labels <- function(plot, x_angle = 0, y_angle = 0, x_hjust = NULL, y_hjust = NULL,
    x_size = NULL, y_size = NULL, x_face = "plain", y_face = "plain", x_color = "black",
    y_color = "black", bold_title = TRUE) {

    theme_list <- list()

    # X-axis label formatting
    if (x_angle != 0 || !is.null(x_size) || x_face != "plain" || x_color != "black") {
        x_hjust_use <- x_hjust %||% (if (x_angle != 0)
            1 else 0.5)
        x_vjust_use <- if (x_angle != 0)
            0.5 else 1
        theme_list$axis.text.x <- ggplot2::element_text(angle = x_angle, hjust = x_hjust_use,
            vjust = x_vjust_use, size = x_size, face = x_face, color = x_color)
    }

    # Y-axis label formatting
    if (y_angle != 0 || !is.null(y_size) || y_face != "plain" || y_color != "black") {
        y_hjust_use <- y_hjust %||% (if (y_angle != 0)
            1 else 0.5)
        theme_list$axis.text.y <- ggplot2::element_text(angle = y_angle, hjust = y_hjust_use,
            size = y_size, face = y_face, color = y_color)
    }

    # Axis titles
    if (bold_title) {
        theme_list$axis.title <- ggplot2::element_text(face = "bold")
    }

    if (length(theme_list) > 0) {
        plot <- plot + do.call(ggplot2::theme, theme_list)
    }

    return(plot)
}

#' Calculate Scaled Font Sizes
#'
#' Standardized font size calculation for proportional scaling across output dimensions.
#' Used when plots need to scale font sizes relative to output size (e.g., for heatmaps).
#'
#' @param base_size Numeric: base font size (default: 11)
#' @param scale_factor Numeric: overall scaling multiplier (default: 1)
#' @param font_multipliers List: named ratios relative to base size
#'
#' @return List with named elements: base, scaled, axis_text, axis_title, title, legend, subtitle, caption
#'
#' @examples
#' \dontrun{
#' # Standard font scaling
#' fonts <- .calculate_scaled_fonts(base_size = 11, scale_factor = 1.2)
#' # Use: fonts$title, fonts$axis_text, fonts$legend, etc.
#'
#' # Custom multipliers
#' custom_mult <- list(axis_text = 1.0, title = 1.8, legend = 0.8)
#' fonts <- .calculate_scaled_fonts(scale_factor = 2, font_multipliers = custom_mult)
#' }
#'
#' @noRd

.calculate_scaled_fonts <- function(base_size = 11, scale_factor = 1, font_multipliers = list(axis_text = 12/11,
    axis_title = 14/11, title = 16/11, legend = 9/11)) {

    scaling <- base_size * scale_factor

    result <- list(base = base_size, scaled = scaling, axis_text = round(scaling *
        font_multipliers$axis_text), axis_title = round(scaling * font_multipliers$axis_title),
        title = round(scaling * font_multipliers$title), legend = round(scaling *
            font_multipliers$legend), subtitle = round(scaling * 0.9), caption = round(scaling *
            0.8))

    return(result)
}

#' Apply Group Aesthetic Scales (Color + Fill + Legend)
#'
#' Consolidated helper for color palette + manual scales + legend styling.
#' Replaces repeated 21x pattern of palette + scale_color_manual + scale_fill_manual.
#'
#' AESTHETIC PRESERVATION: Uses exact same colors and mappings as originals.
#' - Palette: .palette_blue_red() by default (standard TSENAT convention)
#' - Color/Fill Manual: with name='Group' (standard legend title)
#'
#' @param plot ggplot2 object
#' @param palette Character: palette function name (e.g. 'palette_blue_red') OR 
#'   a vector of colors. If character, will call the .palette_* function.
#' @param legend_name Character: legend title (default: 'Group')
#' @param legend_position Character: legend position (default: 'bottom')
#' @param direction Integer: 1 (normal) or -1 (reversed palette)
#'
#' @return Modified ggplot2 object with applied color scales
#'
#' @noRd
.apply_group_aesthetics <- function(plot, palette = "palette_blue_red", legend_name = "Group",
    legend_position = "bottom", direction = 1) {

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
    result <- plot + ggplot2::scale_color_manual(values = colors, name = legend_name) +
        ggplot2::scale_fill_manual(values = colors, name = legend_name) + ggplot2::theme(legend.position = legend_position)

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
#' @param x_col Character: name of x column (default: 'q')
#' @param y_col Character: name of y column (default: 'median')
#' @param group_col Character: optional grouping column (NULL = single series, default: NULL)
#' @param ci_lower_col Character: name of CI lower column (default: 'ci_lower')
#' @param ci_upper_col Character: name of CI upper column (default: 'ci_upper')
#' @param ribbon_alpha Numeric: ribbon transparency (default: 0.15)
#' @param line_width Numeric: line width (default: 1.2)
#' @param point_size Numeric: point size (default: 3.5)
#' @param show_points Logical: include geom_point layer? (default: TRUE)
#'
#' @return Base ggplot2 object with ribbon/line/point layers (unthemed)
#'
#' @noRd
.create_ci_ribbon_plot <- function(data, x_col = "q", y_col = "median", group_col = NULL,
    ci_lower_col = "ci_lower", ci_upper_col = "ci_upper", ribbon_alpha = 0.15, line_width = 1.2,
    point_size = 2.8, show_points = TRUE, default_color = "#4575B4") {

    # Check if we have valid CI data to plot ribbons
    has_valid_ci <- FALSE
    if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
        has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
    }

    # Only keep rows with valid CI data for ribbon layer to avoid ggplot
    # warnings But keep data as-is for line/point layers (they use median
    # values, not CI bounds)
    data_ci <- data
    if (has_valid_ci) {
        # Filter to rows with valid CI for ribbon layer only
        valid_ci_rows <- which(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
        data_ci <- data[valid_ci_rows, ]
    }

    # Build base aesthetics - include group color/fill only if group_col
    # provided and exists
    if (!is.null(group_col) && group_col %in% colnames(data)) {
        p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
            color = .data[[group_col]], fill = .data[[group_col]], group = .data[[group_col]]))
        has_grouping <- TRUE
    } else {
        # No grouping - simple x/y aesthetics (color applied as fixed
        # aesthetic)
        p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
        has_grouping <- FALSE
    }

    # Check if we have valid CI data to plot
    has_valid_ci <- FALSE
    if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
        has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
    }

    # Add ribbon layer (CI bounds) only if we have valid CI data
    if (has_valid_ci) {
        if (has_grouping) {
            # When grouping, fill aesthetic is inherited from base aes
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ci_lower_col]],
                ymax = .data[[ci_upper_col]]), alpha = ribbon_alpha, color = NA)
        } else {
            # No grouping - apply default fill color
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ci_lower_col]],
                ymax = .data[[ci_upper_col]]), alpha = ribbon_alpha, color = NA,
                fill = default_color)
        }
    }

    # Add line layer If group_col is provided, color aesthetic from aes()
    # applies it Otherwise, apply default color (consistency with
    # .create_simple_line_plot)
    if (has_grouping) {
        # When grouping, color aesthetic is inherited from base aes
        p <- p + ggplot2::geom_line(linewidth = line_width)
    } else {
        # No grouping - apply default color
        p <- p + ggplot2::geom_line(linewidth = line_width, color = default_color)
    }

    # Add point layer if requested If group_col is provided, color aesthetic
    # from aes() applies it Otherwise, apply default color
    if (show_points) {
        if (has_grouping) {
            # When grouping, color aesthetic is inherited from base aes
            p <- p + ggplot2::geom_point(size = point_size, alpha = 0.8)
        } else {
            # No grouping - apply default color
            p <- p + ggplot2::geom_point(size = point_size, alpha = 0.8, color = default_color)
        }
    }

    p
}

#' Assemble Multi-Plot Grid with Shared Legend
#'
#' Consolidated helper for legend extraction + grid composition pattern.
#' Replaces repeated 21x pattern of get_legend + plot_nolegend + plot_grid assembly.
#'
#' AESTHETIC PRESERVATION: Uses exact parameters from original code:
#' - Legend position: 'bottom' (default, customizable)
#' - Legend direction: 'horizontal' (standard for TSENAT)
#' - Grid alignment: 'hv' (both axes aligned)
#' - Title/subtitle: uses .font_sizes constants
#'
#' @param plots List of ggplot2 objects (one per subplot)
#' @param ncol Integer: number of columns (default: 2)
#' @param nrow Integer: number of rows (default: auto-calculated)
#' @param title Character: main title (optional)
#' @param subtitle Character: subtitle under title (optional)
#' @param legend_position Character: position for legend - 'bottom', 'top', 'left', 'right', 'none' 
#'   (default: 'bottom')
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
.assemble_grid_plot <- function(plots, ncol = 2, nrow = NULL, title = NULL, subtitle = NULL,
    legend_position = "bottom", extract_legend = TRUE, rel_heights = c(0.08, 1, 0.08)) {

    if (length(plots) == 0) {
        stop("plots list cannot be empty", call. = FALSE)
    }

    # Calculate nrow if not provided
    if (is.null(nrow)) {
        nrow <- ceiling(length(plots)/ncol)
    }

    # Extract legend from first plot if requested
    legend_obj <- NULL
    if (extract_legend) {
        legend_obj <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = legend_position,
            legend.direction = "horizontal"))
    }

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(p) {
        .configure_legend(p, position = "none")
    })

    # Compose grid without legend
    grid_plot <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow,
        align = "hv")

    # If title/subtitle provided, create title grobs and assemble all 3
    # components
    if (!is.null(title) || !is.null(subtitle)) {
        title_plot <- cowplot::ggdraw()

        if (!is.null(title)) {
            title_plot <- title_plot + cowplot::draw_label(title, fontface = "bold",
                size = .font_sizes$title, x = 0.5, hjust = 0.5)
        }

        if (!is.null(subtitle)) {
            y_pos <- if (is.null(title))
                0.5 else 0.25
            subtitle_plot <- cowplot::ggdraw() + cowplot::draw_label(subtitle, fontface = "italic",
                size = .font_sizes$subtitle, x = 0.5, hjust = 0.5, color = "gray40")
            title_plot <- cowplot::plot_grid(title_plot, subtitle_plot, nrow = 2,
                rel_heights = c(1, 0.6))
        }

        # Assemble title + grid + legend (if extracted)
        if (extract_legend && !is.null(legend_obj)) {
            return(cowplot::plot_grid(title_plot, grid_plot, legend_obj, nrow = 3,
                rel_heights = rel_heights))
        } else {
            # No legend: just title + grid
            return(cowplot::plot_grid(title_plot, grid_plot, nrow = 2, rel_heights = rel_heights[c(1,
                2)]))
        }
    }

    # No title/subtitle: just grid + legend (if extracted)
    if (extract_legend && !is.null(legend_obj)) {
        return(cowplot::plot_grid(grid_plot, legend_obj, nrow = 2, rel_heights = rel_heights[c(2,
            3)]))
    } else {
        return(grid_plot)
    }
}

# ============================================================================
# PHASE 2 HELPERS: MEDIUM-IMPACT PATTERN CONSOLIDATION
# ============================================================================

#' Create Cowplot Title Grob with Optional Subtitle
#'
#' Consolidates cowplot::ggdraw() + draw_label() pattern (7x occurrences).
#' Creates a title-only or title+subtitle grob for use in grid layouts.
#'
#' @param title Character: main title text
#' @param subtitle Character: optional subtitle text
#' @param title_size Numeric: title font size (default: .font_sizes$title)
#' @param subtitle_size Numeric: subtitle font size (default: .font_sizes$subtitle)
#' @param title_face Character: title font face ('bold', 'italic', etc.)
#' @param subtitle_face Character: subtitle font face
#' @param title_color Character: title color (default: 'black')
#' @param subtitle_color Character: subtitle color (default: 'gray40')
#'
#' @return cowplot/ggplot2 grob object ready for plot_grid assembly
#'
#' @noRd
.create_title_grob <- function(title, subtitle = NULL, title_size = .font_sizes$title,
    subtitle_size = .font_sizes$subtitle, title_face = "bold", subtitle_face = "italic",
    title_color = "black", subtitle_color = "gray40") {

    # Start with title grob
    title_grob <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = title_face,
        size = title_size, x = 0.5, hjust = 0.5, color = title_color)

    # Add subtitle if provided
    if (!is.null(subtitle)) {
        subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(subtitle, fontface = subtitle_face,
            size = subtitle_size, x = 0.5, hjust = 0.5, color = subtitle_color)

        # Combine title + subtitle
        title_grob <- cowplot::plot_grid(title_grob, subtitle_grob, nrow = 2, rel_heights = c(1,
            0.6))
    }

    title_grob
}

#' Apply Facet Styling with Panel Spacing and Strip Text
#'
#' Consolidates facet_wrap() + panel.spacing + strip.text pattern (8x occurrences).
#' Applies consistent faceting and panel styling across all plot types.
#'
#' @param plot ggplot2 object
#' @param ncol Integer: number of columns for facet layout
#' @param nrow Integer: number of rows (optional, usually auto-calculated)
#' @param facet_var Character: variable name to facet by (unquoted expression as string)
#' @param scales Character: 'fixed', 'free_x', 'free_y', or 'free' (default: 'free_y')
#' @param strip_text_size Numeric: font size for strip labels (default: .font_sizes$subtitle)
#' @param panel_spacing_lines Numeric: spacing between panels in lines (default: 1.5)
#'
#' @return Modified ggplot2 object with faceting and styling applied
#'
#' @noRd
.apply_facet_styling <- function(plot, ncol = 2, nrow = NULL, facet_var = NULL, scales = "free_y",
    strip_text_size = .font_sizes$subtitle, panel_spacing_lines = 1.5) {

    # Apply facet wrap if variable specified
    if (!is.null(facet_var)) {
        facet_formula <- stats::as.formula(paste0("~", facet_var))
        plot <- plot + ggplot2::facet_wrap(facet_formula, ncol = ncol, nrow = nrow,
            scales = scales)
    }

    # Apply panel and strip styling
    plot <- plot + ggplot2::theme(panel.spacing = ggplot2::unit(panel_spacing_lines,
        "lines"), strip.text = ggplot2::element_text(face = "bold", size = strip_text_size))

    plot
}

#' Prepare Long Format Data with Metadata Validation
#'
#' Wrapper around .prepare_tsallis_long() that adds validation and 
#' consistent parameter handling. Consolidates 10x data preparation pattern.
#'
#' @param se SummarizedExperiment: diversity/entropy assay object
#' @param assay_name Character: name of assay to use (default: 'diversity')
#' @param condition_col Character: column name for condition/group variable
#'   (optional, uses metadata config if NULL)
#' @param validate Logical: validate output structure (default: TRUE)
#'
#' @return Data frame in long format with columns: q, tsallis, group, Gene
#'
#' @noRd
.prepare_long_format <- function(se, assay_name = "diversity", condition_col = NULL,
    validate = TRUE) {

    # Use condition_col from metadata config if not provided
    if (is.null(condition_col)) {
        if (!is.null(S4Vectors::metadata(se)$condition_col)) {
            condition_col <- S4Vectors::metadata(se)$condition_col
        } else {
            condition_col <- "sample_type"  # TSENAT default
        }
    }

    # Prepare long format
    long <- .prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)

    # Validate structure
    if (validate) {
        required_cols <- c("q", "tsallis", "group", "Gene")
        missing_cols <- setdiff(required_cols, colnames(long))
        if (length(missing_cols) > 0) {
            warning("Prepared data missing columns: ", paste(missing_cols, collapse = ", "))
        }
    }

    long
}

#' Create Centered Theme Element Components
#'
#' Helper for extracting and applying centered/bolded title/subtitle styling.
#' Consolidates 12x centering + bolding theme pattern.
#'
#' @param include_title Logical: include title element (default: TRUE)
#' @param include_subtitle Logical: include subtitle element (default: TRUE)
#' @param title_size Numeric: title size (default: .font_sizes$title)
#' @param subtitle_size Numeric: subtitle size (default: .font_sizes$subtitle)
#' @param hjust Numeric: horizontal justification (default: 0.5 = centered)
#'
#' @return ggplot2::theme() object with centering/styling
#'
#' @noRd
.create_centered_theme <- function(include_title = TRUE, include_subtitle = TRUE,
    title_size = .font_sizes$title, subtitle_size = .font_sizes$subtitle, hjust = 0.5) {

    theme_list <- list()

    if (include_title) {
        theme_list$plot.title <- ggplot2::element_text(hjust = hjust, size = title_size,
            face = "bold")
    }

    if (include_subtitle) {
        theme_list$plot.subtitle <- ggplot2::element_text(hjust = hjust, size = subtitle_size,
            face = "italic")
    }

    do.call(ggplot2::theme, theme_list)
}



#' Create Publication-Ready Line Plot
#'
#' Consolidates simple line + point layer patterns (8x occurrences).
#' Creates line + point layers with optional grouping.
#'
#' @param data Data frame with x and y columns
#' @param x_col Character: x-axis column name
#' @param y_col Character: y-axis column name
#' @param group_col Character: optional grouping column (NULL = single series with fixed color)
#' @param points Logical: add point layer (default: TRUE)
#' @param line_width Numeric: line width (default: 1.2)
#' @param point_size Numeric: point size (default: 2.5)
#' @param alpha Numeric: transparency (default: 0.8)
#' @param line_color Character: fixed line color when no grouping (default: '#4575B4')
#'
#' @return ggplot2 object with line and optional point layers (unthemed)
#'
#' @noRd
.create_simple_line_plot <- function(data, x_col, y_col, group_col = NULL, points = TRUE,
    line_width = 1.2, point_size = 2, alpha = 0.8, line_color = "#4575B4") {

    # NO grouping: simple single-series plot with fixed color
    if (is.null(group_col)) {
        p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col))) +
            ggplot2::geom_line(linewidth = line_width, alpha = alpha, color = line_color)

        if (isTRUE(points)) {
            p <- p + ggplot2::geom_point(size = point_size, alpha = alpha, color = line_color)
        }
        return(p)
    }

    # WITH grouping: map color to group column
    if (!(group_col %in% colnames(data))) {
        stop("Column '", group_col, "' not found in data")
    }

    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col),
        color = !!rlang::sym(group_col))) + ggplot2::geom_line(linewidth = line_width,
        alpha = alpha)

    if (isTRUE(points)) {
        p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
    }

    p
}

#' Prepare Grouped Long Format Data with Aggregation
#'
#' Consolidates data aggregation + IQR/SD calculation patterns (6x occurrences).
#' Transforms data from wide format with groups into long format with statistics.
#'
#' @param se SummarizedExperiment: input data
#' @param assay_name Character: assay to transform
#' @param group_by_col Character: column to group by
#' @param stat_funcs List of functions: functions to apply (default: median, IQR)
#'   Named list like list(median = median, iqr = function(x) diff(quantile(x, c(0.25, 0.75))))
#'
#' @return Data frame in long format with grouped statistics
#'   Columns: q (from rownames or metadata), group, value, lower, upper
#'
#' @noRd
.prepare_grouped_long_format <- function(se, assay_name = "diversity", group_by_col = "condition",
    stat_funcs = list(median = median, iqr = function(x) diff(quantile(x, c(0.25,
        0.75), na.rm = TRUE)))) {

    # Extract assay
    assay_mat <- SummarizedExperiment::assay(se, assay_name)

    # Get group info
    coldata <- SummarizedExperiment::colData(se)
    if (!group_by_col %in% colnames(coldata)) {
        stop("Column '", group_by_col, "' not found in colData")
    }
    groups <- coldata[[group_by_col]]

    # Aggregate by group
    long_list <- list()

    for (i in seq_len(nrow(assay_mat))) {
        gene_name <- rownames(assay_mat)[i]
        gene_data <- assay_mat[i, ]

        for (grp in unique(groups)) {
            grp_indices <- which(groups == grp)
            grp_values <- gene_data[grp_indices]
            grp_values <- grp_values[!is.na(grp_values)]

            if (length(grp_values) > 0) {
                median_val <- stat_funcs$median(grp_values)
                iqr_val <- stat_funcs$iqr(grp_values)

                long_list[[paste0(gene_name, "_", grp)]] <- data.frame(Gene = gene_name,
                  group = grp, value = median_val, lower = median_val - (iqr_val/2),
                  upper = median_val + (iqr_val/2), stringsAsFactors = FALSE)
            }
        }
    }

    do.call(rbind, long_list)
}  # ============================================================================
# PHASE 6: PLOT CONSOLIDATION HELPERS
# ============================================================================
# These helpers consolidate remaining high-value patterns identified in Phase
# 6: - Theme + aesthetics merging - Plot saving with unified dimensions -
# Distribution statistics abstraction - Bootstrap CI detection and extraction

#' Apply Publication Theme + Group Aesthetics Combined
#'
#' Consolidates the common pattern of applying both publication theme and
#' group-based color aesthetics in a single call.
#'
#' @param plot ggplot2 object to modify
#' @param title Character: plot title (optional)
#' @param base_size Numeric: base font size (default: 11)
#' @param base_theme Character: theme function name - 'theme_base' or 'theme_spectrum'
#' @param group_col Character: column name for group mapping (optional)
#' @param palette Character: palette name - 'blue_red', custom function name (default: 'blue_red')
#' @param group_levels Character vector: ordered factor levels (optional)
#' @param subtitle Character: plot subtitle (optional)
#'
#' @return Modified ggplot2 object
#'
#' @details
#' This function is a convenient wrapper around `.apply_publication_theme()` and
#' `.apply_group_aesthetics()` for the common case where both are applied sequentially.
#'
#' If `group_col` is NULL, only the publication theme is applied.
#'
#' @examples
#' \dontrun{
#' # Apply theme + group colors in one call
#' p <- ggplot2::ggplot(df, ggplot2::aes(x = q, y = entropy, color = group)) +
#'     ggplot2::geom_point()
#' p <- .apply_publication_aesthetics(p, title = 'Entropy Trend',
#'                                    base_size = 11, base_theme = 'theme_base',
#'                                    group_col = 'group', palette = 'blue_red')
#' }
#'
#' @noRd
.apply_publication_aesthetics <- function(plot, title = NULL, base_size = 11, base_theme = "theme_base",
    group_col = NULL, palette = "blue_red", group_levels = NULL, subtitle = NULL, legend_name = "Group",
    legend_position = "bottom") {

    # Apply publication theme first
    p <- .apply_publication_theme(plot, title = title, base_size = base_size, base_theme = base_theme,
        subtitle = subtitle)

    # Apply group aesthetics if group column specified
    if (!is.null(group_col)) {
        p <- .apply_group_aesthetics(p, palette = palette, legend_name = legend_name,
            legend_position = legend_position)
    }

    return(p)
}

# ============================================================================

#' Save Plot with Standard Dimensions
#'
#' Consolidates the common pattern of calculating plot dimensions and saving
#' with ggsave into a single function call.
#'
#' @param plot ggplot2 object to save
#' @param filename Character: output file path (PNG, PDF, etc.)
#' @param width_inches Numeric: chart width in inches (default: 12)
#' @param aspect_type Character: aspect ratio - 'standard' (16:9), 'square' (1:1),
#'   'wide' (21:9), 'tall' (9:16). Default: 'standard'
#' @param dpi_output Numeric: resolution in DPI (default: 100)
#' @param width_cm Numeric: override width in cm (optional)
#' @param height_cm Numeric: override height in cm (optional)
#'
#' @return Invisible NULL. Saves plot to file as side effect.
#'
#' @details
#' This function eliminates the repetitive pattern:
#' ```
#' plot_dims <- .calculate_plot_dims(width_inches = 12, aspect_type = 'standard')
#' ggplot2::ggsave(file, plot = p, width = plot_dims$width, height = plot_dims$height, dpi = plot_dims$dpi)
#' ```
#'
#' **Consolidation Impact**: 20+ occurrences × 3-4 lines = 60+ LOC saved
#'
#' @examples
#' \dontrun{
#' # Save plot with standard dimensions
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
#'     ggplot2::geom_point()
#' .save_plot_standard(p, 'output/my_plot.png', width_inches = 12)
#' }
#'
#' @noRd
.save_plot_standard <- function(plot, filename, width_inches = 12, aspect_type = "standard",
    dpi_output = 100, width_cm = NULL, height_cm = NULL) {

    # Calculate dimensions
    plot_dims <- .calculate_plot_dims(width_inches = width_inches, aspect_type = aspect_type,
        dpi_output = dpi_output)

    # Override with cm if provided (convert cm to inches)
    if (!is.null(width_cm)) {
        plot_dims$width <- width_cm/2.54
    }
    if (!is.null(height_cm)) {
        plot_dims$height <- height_cm/2.54
    }

    # Save plot
    ggplot2::ggsave(filename, plot = plot, width = plot_dims$width, height = plot_dims$height,
        dpi = plot_dims$dpi)

    invisible(NULL)
}

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

.plot_diversity_violin_singleq <- function(se, assay_name = "diversity", title = NULL) {

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

    p <- .apply_publication_theme(p, title = title_use, base_size = 11) + ggplot2::labs(x = "Group",
        y = "Tsallis entropy")

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

.plot_diversity_density_singleq <- function(se, assay_name = "diversity", title = NULL) {

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

    p <- .apply_publication_theme(p, title = title_use, base_size = 11) + ggplot2::labs(x = "Tsallis entropy",
        y = "Density")

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
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' analysis <- calculate_diversity(analysis, q = 1.0)
#' p <- plot_diversity_violin_density(analysis)
#' # if (!is.null(p)) print(p)
#'
plot_diversity_violin_density <- function(se, assay_name = "diversity", title = NULL,
    output_file = NULL) {
    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Handle TSENATAnalysis objects - extract first diversity result
    if (methods::is(se, "TSENATAnalysis")) {
        if (length(se@diversity_results) == 0) {
            stop("No diversity results found in TSENATAnalysis object. Run calculate_diversity() first.")
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
    p_violin <- .plot_diversity_violin_singleq(se = se, assay_name = assay_name, title = "Violin")

    p_density <- .plot_diversity_density_singleq(se = se, assay_name = assay_name,
        title = "Density")

    # Arrange plots side by side: violin on left, density on right
    grid <- cowplot::plot_grid(p_violin, p_density, nrow = 1, ncol = 2, align = "h",
        axis = "b")

    # Add overall title and subtitle above the grid
    title_grob <- .create_title_grob("Tsallis Entropy Distribution by Group", subtitle = "Violin and density plots across samples",
        title_size = 19, subtitle_size = 15)
    grid_with_title <- cowplot::plot_grid(title_grob, grid, nrow = 2, rel_heights = c(0.08,
        1))

    # Save to file if output_file is provided
    if (!is.null(output_file)) {
        plot_dims <- .calculate_plot_dims(width_inches = 12, aspect_type = "standard",
            dpi_output = 100)
        ggplot2::ggsave(output_file, plot = grid_with_title, width = plot_dims$width,
            height = plot_dims$height, dpi = plot_dims$dpi, create.dir = TRUE)
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
    built <- .make_plot_for_genebuild_tx_long(gene_single, mapping, counts, samples,
        NULL)
    df_summary <- .make_plot_for_geneaggregate_df_long(built$df_long, agg_fun, pseudocount)
    .make_plot_for_genebuild_plot_from_summary(df_summary, agg_label_unique, fill_limits,
        font_scale = font_scale)
}


.make_plot_for_genecombine_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
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
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_font, hjust = 0.5,
            face = "bold"), plot.margin = ggplot2::margin(4, 4, 4, 4))

    # Apply axis label formatting and legend configuration using Phase 5
    # helpers
    p <- .format_axis_labels(p, x_size = x_axis_font, y_size = y_axis_font, y_face = "plain",
        bold_title = FALSE)
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
    p_for_legend <- .configure_legend(plots[[1]], position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)
    plots_nolegend <- lapply(plots, function(pp) .configure_legend(pp, position = "none"))

    # Use 2 columns (2 genes per row), auto-calculate rows
    ncol <- 2
    nrow_val <- ceiling(length(plots_nolegend)/ncol)

    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val,
        align = "hv")
    title_grob <- .create_title_grob("Transcript level expression", subtitle = paste0("Top genes with metric ",
        agg_label_unique), title_size = 18, subtitle_size = 14)
    # Add spacer between title and plots
    spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()
    result_plot <- cowplot::plot_grid(title_grob, spacer_grob, grid, legend, ncol = 1,
        rel_heights = c(0.09, 0.0015, 1, 0.08), align = "h", axis = "l")
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    }
    result_plot
}


.make_plot_for_genecombine_grid <- function(plots, output_file = NULL, agg_label_unique) {
    plots_nolegend <- lapply(plots, function(pp) .configure_legend(pp, position = "none"))
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

# [REMOVED] .prepare_ma_plot_df() - MA plot functionality removed (April 2026)

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
#'   Typically the result from [.calculate_effect_sizes()].
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
            color = "black") + ggplot2::labs(title = expression(bold("Distribution of Tsallis Divergence (" ~
        D[q] ~ ") effect sizes across genes")), subtitle = "Information-theoretic measure respecting Tsallis multi-q entropy properties",
        x = bquote("Effect size (Tsallis Divergence" ~ D[q] ~ "; D >" ~ .(threshold) ~
            "= meaningful information separation)"), y = "Number of genes", caption = paste("Red dashed line: D =",
            threshold, "filtering threshold (information-theoretic significance for q-dependent entropy)")) +
        .theme_base(base_size = 11) + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray90"),
            plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5))

    # Add reference line using Phase 5 helper
    p_effect <- .add_reference_lines(p_effect, v_intercept = threshold, v_color = "red",
        v_size = 1)

    # Add threshold annotation
    p_effect <- p_effect + ggplot2::annotate("text", x = threshold, y = Inf, label = paste("Information\nthreshold\n(D=",
        threshold, ")", sep = ""), vjust = 1.5, hjust = -0.1, color = "red", size = 3.5)

    return(p_effect)
}




