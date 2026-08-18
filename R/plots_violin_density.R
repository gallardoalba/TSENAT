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
#' @param q \code{numeric} or \code{NULL}. Single q value to plot. If NULL,
#'   the q value must be unambiguous from the stored diversity results
#'   (exactly one q); otherwise an error is raised listing the available
#'   q values.
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
    output_file = NULL, q = NULL) {
    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Normalize input (handles TSENATAnalysis → SE, condition_col, assay
    # validation). For single-q plots the q value MUST be unambiguous:
    # explicit q=, or exactly one q stored in the object.
    normalized <- .normalize_plot_input(se, assay_name = assay_name, multi_q = FALSE,
        q = q)
    se <- normalized$se

    # Extract q-value and long-format data
    extracted <- .extract_q_value(se, assay_name = assay_name)
    q_val <- extracted$q_val

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
    title_grob <- .create_title_grob(base_title, subtitle = "Violin and density plots across samples")
    grid_with_title <- cowplot::plot_grid(title_grob, grid, nrow = 2, rel_heights = c(0.08,
        1))

    # Unified save-or-return
    .finalize_plot(grid_with_title, output_file)
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
#' @param assay_name Name of the assay to use (default: 'diversity').
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a density plot colored by group.
#'
#' @noRd
.plot_diversity_density_singleq <- function(se, assay_name = "diversity", title = NULL) {

    # Extract q-value and long-format data
    extracted <- .extract_q_value(se, assay_name = assay_name)
    q_val <- extracted$q_val
    long <- extracted$long

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
#' @param assay_name Name of the assay to use (default: 'diversity').
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a violin plot with groups on the x-axis.
#'  
#' @noRd
.plot_diversity_violin_singleq <- function(se, assay_name = "diversity", title = NULL) {

    # Extract q-value and long-format data
    extracted <- .extract_q_value(se, assay_name = assay_name)
    q_val <- extracted$q_val
    long <- extracted$long

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
