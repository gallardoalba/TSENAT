#' Generate Concordance Comparison Plots for GAM vs Friedman Test Methods
#'
#' Creates a two-panel visualization comparing p-value results and distributions
#' from GAM and Friedman statistical tests, highlighting agreement between methods.
#'
#' @param comparison_df A data.frame with columns:
#'   \itemize{
#'     \item \code{gene}: Gene names
#'     \item \code{p_gam}: GAM p-values
#'     \item \code{p_friedman}: Friedman test p-values
#'     \item \code{agreement}: Categorical variable indicating agreement type
#'       (e.g., "Both significant", "GAM only", "Friedman only", "Neither significant")
#'   }
#'
#' @return A gridExtra grob object containing the combined two-panel plot.
#'   Panel 1: Scatter plot of -log10(p-values) with significance thresholds.
#'   Panel 2: Histogram of p-value distributions by method.
#'
#' @details
#' The function requires ggplot2 and gridExtra packages. The scatter plot shows
#' agreement categories with distinct colors and reference lines at p=0.05
#' significance threshold. The histogram compares the distribution of p-values
#' between the two methods.
#'
#' @examples
#' # Create synthetic comparison data
#' set.seed(123)
#' comparison_df <- data.frame(
#'   gene = paste0("gene_", 1:50),
#'   p_gam = runif(50, 0, 0.5),
#'   p_friedman = runif(50, 0, 0.5),
#'   agreement = sample(c("Both significant", "GAM only", "Friedman only", "Neither significant"),
#'                      size = 50, replace = TRUE)
#' )
#' 
#' # Create concordance plot
#' plot <- plot_method_concordance(comparison_df)
#'
#' @keywords internal
#' @noRd
plot_method_concordance <- function(comparison_df) {
  
  # Check if data is valid
  if (is.null(comparison_df) || nrow(comparison_df) == 0) {
    stop("comparison_df must be a non-empty data.frame with p-value columns")
  }
  
  # Check required columns
  required_cols <- c("p_gam", "p_friedman", "agreement")
  if (!all(required_cols %in% colnames(comparison_df))) {
    missing <- setdiff(required_cols, colnames(comparison_df))
    stop("comparison_df missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Create comparison plot
  p1 <- ggplot2::ggplot(comparison_df, 
                        ggplot2::aes(x = -log10(p_gam), 
                                     y = -log10(p_friedman), 
                                     color = agreement)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c("Both significant" = "#2ecc71", 
                 "GAM only" = "#3498db",
                 "Friedman only" = "#e74c3c",
                 "Neither significant" = "#95a5a6"),
      breaks = c("Both significant", "GAM only", "Friedman only", "Neither significant")
    ) +
    ggplot2::labs(
      title = "Method Concordance",
      x = "-log10(p-value, GAM)",
      y = "-log10(p-value, Friedman)",
      color = "Significance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = 14),
      legend.position = "bottomright",
      panel.grid.major = ggplot2::element_line(color = "gray90")
    )
  
  # P-value distribution comparison
  p_long <- data.frame(
    p_value = c(comparison_df$p_gam, comparison_df$p_friedman),
    method = c(rep("GAM", nrow(comparison_df)), rep("Friedman", nrow(comparison_df))),
    stringsAsFactors = FALSE
  )
  
  p2 <- ggplot2::ggplot(p_long, ggplot2::aes(x = p_value, fill = method)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    ggplot2::scale_fill_manual(values = c("GAM" = "#3498db", "Friedman" = "#e74c3c")) +
    ggplot2::labs(
      title = "P-value Distributions",
      x = "P-value",
      y = "Frequency",
      fill = "Method"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "plain", size = 12))
  
  # Combine plots with global title using cowplot approach
  main_grid <- gridExtra::arrangeGrob(p1, p2, ncol = 2)
  
  # Add global title with subtitle using cowplot
  title_grob <- cowplot::ggdraw() + 
    cowplot::draw_label("Comparing interaction detection across two statistical methods",
                       fontface = "bold", size = 16, x = 0.5, y = 0.75) +
    cowplot::draw_label("Concordance analysis between GAM and Friedman tests",
                       fontface = "italic", size = 12, x = 0.5, y = 0.45, color = "gray40")
  
  cowplot::plot_grid(
    title_grob,
    main_grid,
    nrow = 2,
    rel_heights = c(0.15, 1)
  )
}


#' Compute Method Concordance between GAM and Kruskal-Wallis Results
#'
#' Analyzes agreement between GAM (flexible parametric) and Kruskal-Wallis (rank-based)
#' statistical test results. Merges results, calculates correlation, categorizes
#' agreement patterns, and identifies high-confidence genes significant in both methods.
#'
#' @param gam_results A data.frame from GAM analysis with columns:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers
#'     \item \code{p_interaction}: GAM p-values for q-value × group interaction
#'     \item \code{adj_p_interaction}: Adjusted GAM p-values
#'     \item \code{effect_size}: Effect size estimate (optional)
#'   }
#'
#' @param kw_results A data.frame from Friedman test analysis with columns:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers (must match gam_results$gene)
#'     \item \code{p_value}: Friedman test p-values
#'     \item \code{adj_p_value}: Adjusted Friedman p-values (optional)
#'     \item \code{effect_size_eta2}: Effect size estimate (optional)
#'   }
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{comparison_df}: Data frame with merged results and agreement classification
#'       Contains columns: gene, p_gam, padj_gam, effect_gam, p_friedman, padj_friedman, effect_friedman,
#'       gam_sig, friedman_sig, agreement
#'     \item \code{spearman_rho}: Spearman correlation between GAM and Friedman p-values
#'     \item \code{high_conf}: Subset of comparison_df for genes significant in both methods,
#'       ordered by minimum p-value
#'     \item \code{agreement_table}: Table of agreement categories with counts
#'   }
#'
#' @details
#' Agreement categories are defined based on significance at adj_p < 0.05:
#' \itemize{
#'   \item "Both significant": Significant in both GAM and Friedman (most reliable)
#'   \item "GAM only": Significant only in GAM
#'   \item "Friedman only": Significant only in Friedman test
#'   \item "Neither significant": Not significant in either method
#' }
#'
#' High-confidence genes are those reaching p < 0.05 in both methods, indicating
#' robust detection of q-value × group interactions.
#'
#' @examples
#' \dontrun{
#'   # Assuming gam_results and kw_results are available
#'   concordance_result <- compute_method_concordance(gam_results, kw_results)
#'   
#'   # View agreement breakdown
#'   print(concordance_result$agreement_table)
#'   
#'   # View high-confidence genes
#'   head(concordance_result$high_conf)
#'   
#'   # Check correlation
#'   cat("Spearman ρ =", concordance_result$spearman_rho)
#' }  
#'
#' @keywords internal
#' @noRd
compute_method_concordance <- function(gam_results, kw_results) {
  
  # Initialize outputs
  comparison_df <- NULL
  spearman_rho <- NA
  high_conf <- NULL
  agreement_table <- NULL
  
  # Validate inputs
  if (!is.data.frame(gam_results)) {
    stop("gam_results must be a data.frame")
  }
  
  if (!is.data.frame(kw_results)) {
    stop("kw_results must be a data.frame")
  }
  
  if (!("p_interaction" %in% colnames(gam_results))) {
    stop("gam_results missing required column: p_interaction")
  }
  
  if (!("p_value" %in% colnames(kw_results))) {
    stop("kw_results missing required column: p_value")
  }
  
  # Create matching gene sets
  gam_genes <- gam_results$gene[!is.na(gam_results$p_interaction)]
  kw_genes <- kw_results$gene[!is.na(kw_results$p_value)]
  common_genes <- intersect(gam_genes, kw_genes)
  
  if (length(common_genes) > 2) {
    # Extract matching rows by index
    gam_idx <- match(common_genes, gam_results$gene)
    kw_idx <- match(common_genes, kw_results$gene)
    
    # Build comparison data frame with all available columns
    comparison_df <- data.frame(
      gene = common_genes,
      p_gam = gam_results$p_interaction[gam_idx],
      padj_gam = gam_results$adj_p_interaction[gam_idx],
      effect_gam = if ("effect_size" %in% colnames(gam_results)) {
        gam_results$effect_size[gam_idx]
      } else {
        rep(NA_real_, length(gam_idx))
      },
      p_friedman = kw_results$p_value[kw_idx],
      padj_friedman = if ("adj_p_value" %in% colnames(kw_results)) {
        kw_results$adj_p_value[kw_idx]
      } else {
        rep(NA_real_, length(kw_idx))
      },
      effect_friedman = if ("effect_size_eta2" %in% colnames(kw_results)) {
        kw_results$effect_size_eta2[kw_idx]
      } else {
        rep(NA_real_, length(kw_idx))
      },
      stringsAsFactors = FALSE
    )
    
    # Calculate Spearman correlation on ADJUSTED p-values (for consistency with significance threshold)
    spearman_rho <- stats::cor(comparison_df$padj_gam, comparison_df$padj_friedman,
                               method = "spearman", use = "complete.obs")
    
    # Categorize agreement based on adjusted p-value significance (adj_p < 0.05)
    comparison_df$gam_sig <- comparison_df$padj_gam < 0.05
    comparison_df$friedman_sig <- comparison_df$padj_friedman < 0.05
    
    comparison_df$agreement <- ifelse(
      comparison_df$gam_sig & comparison_df$friedman_sig, "Both significant",
      ifelse(comparison_df$gam_sig & !comparison_df$friedman_sig, "GAM only",
             ifelse(!comparison_df$gam_sig & comparison_df$friedman_sig, "Friedman only",
                    "Neither significant"))
    )
    
    # Create agreement frequency table
    agreement_table <- table(comparison_df$agreement)
    
    # Extract high-confidence genes (significant in both methods)
    high_conf <- comparison_df[comparison_df$gam_sig & comparison_df$friedman_sig, ]
    
    # Sort by minimum p-value across methods
    if (nrow(high_conf) > 0) {
      high_conf <- high_conf[order(pmax(high_conf$p_gam, high_conf$p_friedman)), ]
    }
  }
  
  # Return results as list
  list(
    comparison_df = comparison_df,
    spearman_rho = spearman_rho,
    high_conf = high_conf,
    agreement_table = agreement_table
  )
}
