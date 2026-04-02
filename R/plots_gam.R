#' Plot GAM q-curves for top genes identified by FPCA/GAM interaction tests
#'
#' Visualizes smooth q-curve profiles (GAM fits) for selected genes from
#' `.calculate_lm_interaction()` results. Useful for understanding which q-ranges (rare
#' vs. dominant isoforms) drive significant PC differences between groups.
#'
#' @param se A `SummarizedExperiment` with q-sequence diversity values
#'   (multiple q values per sample). Typically output from `.calculate_diversity()`
#'   with multiple q (e.g., q = seq(0.1, 2, by = 0.1)).
#' @param lm_res A `data.frame` from `.calculate_lm_interaction()` with columns
#'   `gene`, `p_interaction`, and `adj_p_interaction`. Can be from method='fpca'
#'   or method='gam'.
#' @param condition_col Column name in `colData(se)` specifying group
#'   assignments for samples (e.g., 'sample_type').
#' @param genes Optional character vector of specific gene names to plot. If provided,
#'   these genes are plotted directly regardless of significance or n_top. If NULL
#'   (default), the top n_top significant genes are selected.
#' @param n_top Number of top genes (by adjusted p-value) to plot (default: 6).
#'   Only used if genes = NULL.
#' @param sig_alpha Significance threshold for adjusted p-values (default: 0.05).
#'   Only used if genes = NULL; filters lm_res to significant genes before selecting top n.
#' @param assay_name Name of the assay in `se` to extract (default: 'diversity').
#' @param model_data Required list from `.calculate_lm_interaction(..., return_model_data = TRUE)$model_data`
#'   containing metadata (q_values, sample configuration, etc.). This is the preferred way to use
#'   this function as it ensures all visualizations are based on the exact analysis configuration.
#'
#' @return A single `ggplot` object with all selected genes arranged in a grid layout 
#'   (2 columns per row). Can be saved with `ggplot2::ggsave()`.
#'
#' @details
#' For each selected gene, this function:
#' 1. Extracts per-sample entropy values across all q values
#' 2. Fits GAM models: entropy ~ s(q, k=...) independently for each group
#' 3. Generates smooth predictions for visualization
#' 4. Overlays predicted curves for each group with a distinct color
#'
#' By providing `model_data` from `.calculate_lm_interaction()`, the function can directly
#' access the q-values used in the original analysis for more accurate visualization.
#'
#' This complements FPCA by providing interpretable visualization of empirical
#' q-curve shape differences that drive PC-level significance.
#'
#' @examples
#' # Create example data with multiple q values
#' set.seed(123)
#' counts <- matrix(
#'   sample(1:100, 60, replace = TRUE),
#'   nrow = 15, ncol = 4
#' )
#' rownames(counts) <- paste0('tx_', 1:15)
#' colnames(counts) <- paste0('sample_', 1:4)
#' genes <- rep(paste0('gene_', 1:5), each = 3)
#' 
#' # Calculate diversity across multiple q values
#' se <- .calculate_diversity(counts, genes = genes, q = seq(0.5, 2, by = 0.5), norm = TRUE)
#' 
#' # Add sample metadata
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' 
#' # Run linear model analysis with model_data  
#' lm_result <- .calculate_lm_interaction(se, condition_col = 'condition', method = 'gam',
#'                                       return_model_data = TRUE)
#' 
#' # Plot GAM curves for top genes
#' if (nrow(lm_result$results) > 0) {
#'   grid_plot <- .plot_lm_interaction_gam(se, lm_result$results, condition_col = 'condition',
#'                                         n_top = 2, model_data = lm_result$model_data)
#' }
#'

#' @noRd
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal scale_color_brewer
#' @importFrom cowplot plot_grid

.plot_lm_interaction_gam <- function(se, lm_res, condition_col = "condition", genes = NULL,
    n_top = 6, sig_alpha = 0.05, assay_name = "diversity", model_data = NULL, output_file = NULL,
    width = NULL, height = NULL) {

    require_pkgs(c("ggplot2", "mgcv", "SummarizedExperiment", "dplyr", "tidyr", "cowplot"))

    # Validate inputs
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }

    # Handle flexible input: lm_res can be either: 1. A data.frame with results
    # (traditional usage) 2. A list with $results and $model_data
    # (return_model_data = TRUE format)
    if (is.list(lm_res) && !is.data.frame(lm_res)) {
        # lm_res is a list with components
        if ("results" %in% names(lm_res) && is.data.frame(lm_res$results)) {
            # Extract results and model_data from the list
            extracted_results <- lm_res$results

            # If model_data not provided, extract from lm_res
            if (is.null(model_data) && "model_data" %in% names(lm_res)) {
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

    # =========================================================================
    # VALIDATE DATA COMPATIBILITY (NEW: Issue #5 validation)
    # =========================================================================
    # Check that SE and LM results have compatible gene sets before processing
    # Note: stop_on_error = FALSE because we explicitly filter genes below
    # (lines 171-175) to handle mismatches (e.g., genes that failed during
    # recalculation)
    data_validation <- .validate_plot_data(se = se, lm_results = lm_res, stop_on_error = FALSE,
        verbose = FALSE  # Suppress verbose; we'll only see output if validation fails
)

    # Validate and extract metadata from model_data
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

    # Extract assay matrix and colData
    mat <- SummarizedExperiment::assay(se, assay_name)
    cdata <- SummarizedExperiment::colData(se)

    if (!condition_col %in% colnames(cdata)) {
        stop(sprintf("Column '%s' not found in colData(se)", condition_col), call. = FALSE)
    }

    # Build sample-to-group mapping using helper
    sample_to_group <- .prepare_sample_group_mapping(cdata, condition_col)

    # Determine which genes to plot using helper
    top_genes <- .plot_select_genes(lm_res, genes = genes, n_top = n_top, sig_alpha = sig_alpha)

    if (is.null(top_genes)) {
        warning(sprintf("No genes significant at alpha = %g", sig_alpha), call. = FALSE)
        return(NULL)
    }

    # Create gene ID to display name mapping from lm_res
    gene_name_map <- setNames(lm_res$gene, lm_res$gene)  # default: use gene ID

    # If lm_res has a gene_name column (e.g., from return_model_data), use it
    if ("gene_name" %in% colnames(lm_res)) {
        gene_name_map <- setNames(lm_res$gene_name, lm_res$gene)
    }

    # Helper to build plot for a single gene using computational helpers
    make_gam_plot <- function(g, gene_display_name = NULL) {
        # Use provided gene name, or look it up from mapping, or default to
        # gene ID
        if (is.null(gene_display_name)) {
            if (g %in% names(gene_name_map)) {
                gene_display_name <- gene_name_map[[g]]
            } else {
                gene_display_name <- g
            }
        }

        # Prepare plot data for this gene using helper
        plot_df <- .plot_gam_prepare_gene_data(g, mat, sample_to_group)
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
            ggplot2::geom_point(data = plot_df, ggplot2::aes(x = q, y = entropy,
                color = group), alpha = 0.5, size = 2) + ggplot2::geom_line(data = pred_df,
            ggplot2::aes(x = q, y = entropy_fit, color = group, linetype = "GAM fit"),
            linewidth = 1, alpha = 0.9) + ggplot2::scale_color_manual(values = color_mapping,
            name = condition_col, breaks = group_levels) + ggplot2::scale_linetype_manual(values = c(`GAM fit` = 1),
            name = "") + ggplot2::labs(x = "q parameter", y = "Tsallis entropy",
            title = ifelse(gene_display_name != g, sprintf("%s (%s)", gene_display_name,
                g), gene_display_name)) + .theme_spectrum(base_size = 11) + ggplot2::theme(legend.position = "none")

        return(p)
    }

    # Generate plots for top genes
    plots <- list()
    for (g in top_genes) {
        p <- make_gam_plot(g)
        if (!is.null(p)) {
            plots[[g]] <- p
        }
    }

    if (length(plots) == 0) {
        warning("No valid plots generated", call. = FALSE)
        return(NULL)
    }

    # Arrange plots in a grid and return single combined plot
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
    legend <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom",
        legend.title = ggplot2::element_text(size = .font_sizes$legend_title), legend.text = ggplot2::element_text(size = .font_sizes$legend_text)))

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

    # Save to file if output_file is provided
    if (!is.null(output_file)) {
        # Use provided dimensions or defaults
        save_width <- if (is.null(width))
            12 else width
        save_height <- if (is.null(height))
            10.3 else height
        ggplot2::ggsave(output_file, plot = final_plot, width = save_width, height = save_height,
            dpi = 100, create.dir = TRUE)
    }

    return(final_plot)
}
