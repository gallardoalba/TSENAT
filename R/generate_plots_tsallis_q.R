#' Plot Tsallis Entropy q-Curve
#'
#' Visualize Tsallis entropy (S_q) as a function of the diversity parameter q across sample groups.
#' Supports three modes: aggregate q-curves (default), gene-specific q-curves (when `gene` provided),
#' or bootstrap confidence interval bands.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity()` with diversity assay.
#'   For CI mode (bootstrap=TRUE), must contain pre-computed bootstrap confidence intervals.
#' @param assay_name Character; name of the assay to plot (default: "diversity").
#' @param condition_col Character; column name in colData indicating group/sample type
#'   (default: "sample_type"). Only used in aggregate and CI modes.
#' @param bootstrap Logical; if TRUE, plots bootstrap confidence interval bands for aggregate mode.
#'   Requires SE to contain pre-computed CI assays (ci_lower/ci_upper).
#'   Requires 2+ q values and exactly 2 groups (default: FALSE).
#' @param gene Character vector (optional); if provided, plot q-curves for specified gene(s).
#'   Overrides default aggregate behavior. When provided, uses median +/- SD for each gene.
#' @param lm_res Data frame (optional); gene interaction test results with `gene` column and
#'   p-value column. Accepts either:
#'   - Results from `calculate_lm_interaction()` (has `adj_p_interaction` or `p_interaction` columns)
#'   - Results from `detect_q_gene_interactions()` (has `adj_p_value` or `p_value` columns from Friedman/Wilcoxon tests)
#'   If provided (and `gene` is NULL), plots top `n_top` genes ranked by p-value.
#'   Useful for plotting significant genes from any interaction analysis.
#' @param n_top Integer or NULL; number of top genes to select from `lm_res` when `gene` is NULL
#'   (default: NULL). When NULL, defaults to showing the single most significant gene (n_top=1),
#'   providing a conservative view of the strongest effect. Set to a numeric value to show that many top genes.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @return
#' **Aggregate mode (gene=NULL, lm_res=NULL)**:
#' - With bootstrap=FALSE: A ggplot object showing median entropy with IQR ribbons.
#' - With bootstrap=TRUE: A ggplot object with bootstrap confidence interval bands.
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Single gene: A ggplot object showing median entropy +/- SD for that gene.
#' - Multiple genes: A grid plot object arranged in 2 rows x 2 columns with a shared legend at the bottom.
#'   The legend appears once beneath the grid, avoiding repetition across subplots.
#'
#' @details
#' **Aggregate mode (default, gene=NULL, lm_res=NULL)**:
#' - Plots median Tsallis entropy +/- IQR across all genes for each group
#' - Works with any SummarizedExperiment from calculate_diversity()
#' - Supports single or multiple q values and any number of groups
#' - No CI data required for basic plots; bootstrap CIs optional
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Plots q-curve separately for each selected gene
#' - Shows median entropy +/- SD (variance) for each gene across q-values and groups
#' - When `lm_res` provided: automatically ranks genes and selects top `n_top` by p-value
#' - Single gene: returns a ggplot object; multiple genes: returns a grid plot (2 rows x 2 columns) with shared legend
#' - For multiple genes: legend appears once at the bottom of the grid to avoid repetition and save space
#' - Useful for highlighting specific genes of interest or significant discoveries
#' - Bootstrap mode not supported in gene-specific mode
#'
#' **Bootstrap CI mode (bootstrap=TRUE in aggregate mode)**:
#' - Displays bootstrap confidence interval bands for each group across q-values
#' - Requires exactly 2 groups for comparison
#' - Requires 2+ q values for q-curve visualization
#' - Requires pre-computed bootstrap CIs from `calculate_diversity(..., bootstrap=TRUE)`
#' - Produces ci_lower, ci_upper assays that properly propagate through entropy transformation
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point theme_minimal
#'   scale_color_manual scale_fill_manual labs theme element_text annotate
#' @importFrom dplyr filter group_by summarise pull
#' @importFrom SummarizedExperiment assayNames assay colData rowData
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' # Plot 7: Tsallis entropy q-curve (combined across all sample diversity)
#' analysis <- TSENAT:::create_test_analysis(n_genes = 8, n_samples_per_group = 25,
#'   q_values = seq(0.1, 3, by = 0.1), seed = 123)
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.1, 3, by = 0.1), verbose = FALSE)
#' p <- plot_tsallis_q_curve_s4(analysis)
#' if (!is.null(p)) print(p)
#'
#' @export
plot_tsallis_q_curve_s4 <- function(
  se,
  assay_name = "diversity",
  condition_col = "sample_type",
  bootstrap = FALSE,
  gene = NULL,
  lm_res = NULL,
  n_top = NULL,
  output_file = NULL
) {
  require_pkgs(c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment", "cowplot"))
  
  # Convert TSENATAnalysis to combined SE if needed
  if (methods::is(se, "TSENATAnalysis")) {
    if (assay_name != "diversity") {
      stop("Assay '", assay_name, "' not found in SummarizedExperiment")
    }
    se <- .prepare_combined_se_from_analysis(se)
    assay_name <- "diversity"
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
    return(.plot_tsallis_gene_specific(se, assay_name, condition_col, gene, lm_res, n_top, output_file))
  }
  
  # Aggregate or bootstrap mode
  long <- prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
  if (nrow(long) == 0) stop("No tsallis values found in SummarizedExperiment")
  
  has_bootstrap_ci <- "ci_lower" %in% SummarizedExperiment::assayNames(se) &&
                      "ci_upper" %in% SummarizedExperiment::assayNames(se)
  
  if (bootstrap && has_bootstrap_ci) {
    return(.plot_tsallis_bootstrap_ci(se, long, output_file))
  } else if (bootstrap && !has_bootstrap_ci) {
    warning("bootstrap=TRUE but CI data not found. Falling back to basic plot.")
  }
  
  # Basic aggregate mode
  .plot_tsallis_basic(long, output_file)
}

# ============================================================================
# GENE-SPECIFIC Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_gene_specific <- function(se, assay_name, condition_col, gene, lm_res, n_top, output_file) {
  require_pkgs(c("ggplot2", "dplyr", "cowplot"))
  
  long <- prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
  
  if (!("Gene" %in% colnames(long))) {
    if ("gene" %in% colnames(long)) {
      long$Gene <- long$gene
    } else {
      se_rownames <- rownames(se)
      if (!is.null(se_rownames) && length(se_rownames) > 0) {
        n_per_gene <- nrow(long) / length(se_rownames)
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
    n_genes_to_plot <- if (is.null(n_top)) 1 else n_top
    genes <- head(genes_ordered, n_genes_to_plot)
  } else {
    genes <- as.character(unlist(gene))
  }
  
  if (length(genes) == 0) stop("No genes selected for plotting")
  
  # Plot single or multiple genes
  make_plot_for_gene <- function(sel) {
    long_g <- long[as.character(long$Gene) == sel, , drop = FALSE]
    if (nrow(long_g) == 0) stop("Gene not found in assay: ", sel)
    
    stats_df <- .compute_gene_stats_by_group(long_g)
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum, ymin = central - spread, ymax = central + spread, fill = group), alpha = 0.2) +
      ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = central, color = group), linewidth = 1.3) +
      ggplot2::labs(title = sel, x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
      ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") + 
      ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
      .tsenat_theme_base(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"))
    p
  }
  
  if (length(genes) == 1) {
    return(make_plot_for_gene(genes))
  }
  
  plots <- lapply(genes, make_plot_for_gene)
  names(plots) <- genes
  
  legend_obj <- cowplot::get_legend(
    plots[[1]] + ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
  )
  
  plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
  grid_with_plots <- do.call(cowplot::plot_grid, c(plots_no_legend, list(nrow = 2, ncol = 2)))
  
  title_plot <- cowplot::ggdraw() + 
    cowplot::draw_label("Tsallis Entropy q-Curve Profile", fontface = "bold", size = 19)
  
  grid_with_legend <- cowplot::plot_grid(title_plot, grid_with_plots, legend_obj, nrow = 3, rel_heights = c(0.12, 1, 0.08))
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = grid_with_legend, width = 12, height = 10, dpi = 100)
  }
  
  grid_with_legend
}

# ============================================================================
# BOOTSTRAP CI Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_bootstrap_ci <- function(se, long, output_file) {
  require_pkgs(c("ggplot2", "SummarizedExperiment"))
  
  long$q <- as.numeric(as.character(long$q))
  unique_q <- sort(unique(long$q))
  
  groups <- unique(sort(long$group))
  if (length(unique_q) < 2) {
    stop("Need at least 2 q values for q-curve analysis")
  }
  if (length(groups) != 2) {
    stop("Expected exactly 2 groups for bootstrap comparison")
  }
  
  plot_df <- .aggregate_bootstrap_ci_by_group(se, long)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = median, color = group, fill = group)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    .tsenat_theme_base(base_size = 11) +
    ggplot2::labs(
      title = "Tsallis Entropy Across Diversity Scales (q-spectrum)",
      subtitle = "Median with 95% confidence intervals",
      x = "q value", y = expression("Tsallis entropy (" * S[q] * ")"),
      color = "Group", fill = "Group"
    )
  
  if (length(groups) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = 12, height = 7.2, dpi = 100)
  }
  
  p
}

# ============================================================================
# BASIC AGGREGATE Q-CURVE PLOTTING
# ============================================================================

.plot_tsallis_basic <- function(long, output_file) {
  require_pkgs("ggplot2")
  
  long$q <- as.numeric(as.character(long$q))
  
  stats_df <- dplyr::summarise(
    dplyr::group_by(long, group, q),
    median = median(tsallis, na.rm = TRUE),
    IQR = stats::IQR(tsallis, na.rm = TRUE),
    .groups = "drop"
  )
  
  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = q, y = median, color = group, fill = group)) +
    ggplot2::geom_line(linewidth = 1.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = median - IQR / 2, ymax = median + IQR / 2), alpha = 0.2, color = NA) +
    ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    .tsenat_theme_base(base_size = 11) +
    ggplot2::labs(
      title = "Tsallis Entropy Across Diversity Scales (q-spectrum)",
      subtitle = "Median +/- IQR",
      x = "q value", y = expression("Tsallis entropy (" * S[q] * ")"),
      color = "Group", fill = "Group"
    )
  
  if (length(unique(long$group)) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = 12, height = 7.2, dpi = 100)
  }
  
  p
}
