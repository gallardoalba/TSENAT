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
#' # Create configuration (required when metadata is provided)
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' analysis <- calculate_diversity(analysis, q = 1.0)
#' p <- plot_diversity_violin_density(analysis)
#' # print(p)
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
    p_violin <- .plot_diversity_violin_singleq(se = se, assay_name = assay_name,
        title = "Violin")

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
        .save_plot_standard(grid_with_title, output_file, width_inches = 12, aspect_type = "standard",
            dpi_output = 100)
    }

    return(grid_with_title)
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
    suppressPackageStartupMessages({
    })

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
    suppressPackageStartupMessages({
    })

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
