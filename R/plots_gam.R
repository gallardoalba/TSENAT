#' Plot GAM q-curves for top genes identified by FPCA/GAM interaction tests
#'
#' Visualizes smooth q-curve profiles (GAM fits) for selected genes from
#' `.calculate_lm()` results. Useful for understanding which
#' q-ranges (rare
#' vs. dominant isoforms) drive significant PC differences between groups.
#'
#' @param se A `SummarizedExperiment` with q-sequence diversity values
#' (multiple q values per sample). Typically output from
#' `.calculate_diversity()`
#'   with multiple q (e.g., q = seq(0.1, 2, by = 0.1)).
#' @param lm_res A `data.frame` from `.calculate_lm()` with columns
#'   `gene`, `p_interaction`, and `adj_p_interaction`. Can be from method='fpca'
#'   or method='gam'.
#' @param condition_col Column name in `colData(se)` specifying group
#'   assignments for samples (e.g., 'sample_type').
#' @param genes Optional character vector of specific gene names to plot. If
#' provided,
#' these genes are plotted directly regardless of significance or n_top. If
#' NULL
#'   (default), the top n_top significant genes are selected.
#' @param n_top Number of top genes (by adjusted p-value) to plot (default: 6).
#'   Only used if genes = NULL.
#' @param sig_alpha Significance threshold for adjusted p-values (default:
#' 0.05).
#' Only used if genes = NULL; filters lm_res to significant genes before
#' selecting top n.
#' @param assay_name Name of the assay in `se` to extract (default:
#' 'diversity').
#' @param model_data Required list from `.calculate_lm(...,
#' return_model_data = TRUE)$model_data`
#' containing metadata (q_values, sample configuration, etc.). This is the
#' preferred way to use
#' this function as it ensures all visualizations are based on the exact
#' analysis configuration.
#'
#' @return A single `ggplot` object with all selected genes arranged in a
#' grid layout
#'   (2 columns per row). Can be saved with `ggplot2::ggsave()`.
#'
#' @details
#' For each selected gene, this function:
#' 1. Extracts per-sample entropy values across all q values
#' 2. Fits GAM models: entropy ~ s(q, k=...) independently for each group
#' 3. Generates smooth predictions for visualization
#' 4. Overlays predicted curves for each group with a distinct color
#'
#' By providing `model_data` from `.calculate_lm()`, the
#' function can directly
#' access the q-values used in the original analysis for more accurate
#' visualization.
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
#' se <- .calculate_diversity(counts, genes = genes, q = seq(0.5, 2.5, by =
#' 0.5), norm = TRUE)
#' 
#' # Add sample metadata
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' 
#' # Run linear model analysis with model_data  
#' lm_result <- .calculate_lm(se, condition_col = 'condition',
#' method = 'gam',
#'                                       return_model_data = TRUE)
#' 
#' # Plot GAM curves for top genes
#' if (nrow(lm_result$results) > 0) {
#' grid_plot <- .plot_lm_interaction_gam(se, lm_result$results,
#' condition_col = 'condition',
#'                                         n_top = 2,
#'  model_data = lm_result$model_data)
#' }
#'

#' @noRd
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs
#' theme_minimal scale_color_brewer
#' @importFrom cowplot plot_grid

.plot_lm_interaction_gam <- function(se, lm_res, condition_col = "condition", genes = NULL,
    n_top = 6, sig_alpha = 0.05, assay_name = "diversity", model_data = NULL, output_file = NULL,
    width = NULL, height = NULL) {

    # Validate and extract inputs
    validated <- .plot_gam_handle_inputs(se, lm_res)
    se <- validated$se
    lm_res <- validated$lm_res
    if (!is.null(validated$model_data)) {
        model_data <- validated$model_data
    }

    # Validate data compatibility
    .validate_plot_data(se = se, lm_results = lm_res, stop_on_error = FALSE, verbose = FALSE)

    # Extract and validate q-values
    q_values <- .plot_gam_extract_q_values(model_data)

    # Match and filter genes
    gene_match <- .plot_gam_match_genes(se, lm_res)
    se <- gene_match$se
    lm_res <- gene_match$lm_res

    # Extract assay matrix and verify condition column exists
    mat <- SummarizedExperiment::assay(se, assay_name)
    cdata <- SummarizedExperiment::colData(se)
    if (!condition_col %in% colnames(cdata)) {
        stop(sprintf("Column '%s' not found in colData(se)", condition_col), call. = FALSE)
    }

    # Build sample-to-group mapping
    sample_to_group <- .prepare_sample_group_mapping(cdata, condition_col)

    # Select genes to plot
    top_genes <- .plot_select_genes(lm_res, genes = genes, n_top = n_top, sig_alpha = sig_alpha)
    if (is.null(top_genes)) {
        warning(sprintf("No genes significant at alpha = %g", sig_alpha), call. = FALSE)
        return(NULL)
    }

    # Create gene display name mapping
    gene_name_map <- .plot_gam_create_gene_map(lm_res)

    # Generate plots for all top genes
    plots <- list()
    for (g in top_genes) {
        p <- .plot_gam_make_plot(g, gene_name_map = gene_name_map, mat = mat, sample_to_group = sample_to_group,
            condition_col = condition_col)
        if (!is.null(p)) {
            plots[[g]] <- p
        }
    }

    if (length(plots) == 0) {
        warning("No valid plots generated", call. = FALSE)
        return(NULL)
    }

    # Arrange plots in grid with title and legend
    final_plot <- .plot_gam_arrange_grid(plots, condition_col, .font_sizes)

    # Save to file if requested
    .plot_gam_save_plot(final_plot, output_file, width, height)

    return(final_plot)
}
