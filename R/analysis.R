#' Generate Concordance Comparison Plots for GAM vs Kruskal-Wallis Methods
#'
#' Creates a two-panel visualization comparing p-value results and distributions
#' from GAM and Kruskal-Wallis statistical tests, highlighting agreement between methods.
#'
#' @param comparison_df A data.frame with columns:
#'   \itemize{
#'     \item \code{gene}: Gene names
#'     \item \code{p_gam}: GAM p-values
#'     \item \code{p_kw}: Kruskal-Wallis p-values
#'     \item \code{agreement}: Categorical variable indicating agreement type
#'       (e.g., "Both significant", "GAM only", "K-W only", "Neither significant")
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
#' \dontrun{
#'   # Assuming comparison_df has been created with GAM and K-W results
#'   plot <- plot_method_concordance(comparison_df)
#'   plot(plot)
#' }
#'
#' @export
plot_method_concordance <- function(comparison_df) {
  
  # Check if data is valid
  if (is.null(comparison_df) || nrow(comparison_df) == 0) {
    stop("comparison_df must be a non-empty data.frame with p-value columns")
  }
  
  # Check required columns
  required_cols <- c("p_gam", "p_kw", "agreement")
  if (!all(required_cols %in% colnames(comparison_df))) {
    missing <- setdiff(required_cols, colnames(comparison_df))
    stop("comparison_df missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Create comparison plot
  p1 <- ggplot2::ggplot(comparison_df, 
                        ggplot2::aes(x = -log10(p_gam), 
                                     y = -log10(p_kw), 
                                     color = agreement)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c("Both significant" = "#2ecc71", 
                 "GAM only" = "#3498db",
                 "K-W only" = "#e74c3c",
                 "Neither significant" = "#95a5a6"),
      breaks = c("Both significant", "GAM only", "K-W only", "Neither significant")
    ) +
    ggplot2::labs(
      title = "GAM vs Kruskal-Wallis: Method Concordance",
      x = "-log10(p-value, GAM)",
      y = "-log10(p-value, K-W)",
      color = "Significance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottomright",
      panel.grid.major = ggplot2::element_line(color = "gray90")
    )
  
  # P-value distribution comparison
  p_long <- data.frame(
    p_value = c(comparison_df$p_gam, comparison_df$p_kw),
    method = c(rep("GAM", nrow(comparison_df)), rep("K-W", nrow(comparison_df))),
    stringsAsFactors = FALSE
  )
  
  p2 <- ggplot2::ggplot(p_long, ggplot2::aes(x = p_value, fill = method)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    ggplot2::scale_fill_manual(values = c("GAM" = "#3498db", "K-W" = "#e74c3c")) +
    ggplot2::labs(
      title = "P-value Distributions",
      x = "P-value",
      y = "Frequency",
      fill = "Method"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
  
  # Combine and return
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}
