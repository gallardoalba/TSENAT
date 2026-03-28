#' Plot Global Divergence q-Curve Across All Genes
#'
#' Visualizes the average (mean/median) Tsallis divergence D_q across all genes
#' as a function of q-value. This provides a **global view** of which diversity scales
#' (rare vs. abundant isoforms) drive the most divergence on average across the dataset,
#' complementing gene-specific divergence profiles.
#'
#' @param divergence_results_se A `SummarizedExperiment` containing pre-computed divergence values.
#'   Rows = genes, columns = q-values. Column names should indicate q-values (e.g., "q_0.5", "q_1.0").
#' @param gene Optional character. If provided, plot divergence spectrum for this specific gene.
#'   If NULL, plot global divergence curve (aggregated across all genes).
#' @param lm_res Optional data.frame with columns for gene identifiers and p-values. Used to select
#'   top genes when gene = NULL and lm_res is provided. Default: NULL.
#' @param n_genes Integer; number of top genes to plot when showing multi-gene spectra (default: 4).
#'   Must be positive. Genes are sorted by p-value significance (lowest p-values first).
#' @param ncol Integer; number of columns in grid layout for multi-gene plots (default: 2).
#'   Must be positive. Number of rows is automatically calculated as ceiling(n_genes / ncol).
#' @param metric Character. Summary statistic for global curve: "median" or "mean". Default: "median".
#'   Only used when gene = NULL.
#' @param variability_metric Character. Error bar type for global curve: "sd" or "iqr". Default: "iqr".
#'   Only used when gene = NULL.
#'
#' @return A `ggplot` object. Gene-specific calls return a line plot.
#'   Global calls return an aggregated curve with variability bands.
#'
#' @details
#' **Gene-specific mode (gene provided)**:
#' - Extracts divergence values for the specified gene across all q-values
#' - Plots as a line chart with points
#' - Reveals whether this gene shows q-dependent divergence patterns
#'
#' **Global mode (gene = NULL)**:
#' - Aggregates divergence across all genes at each q-value
#' - Shows which diversity scales (q-values) drive the most divergence on average
#' - Useful for identifying dominant biological mechanisms (rare vs. abundant isoform driven)
#'
#' **Interpretation**: Compare with `plot_tsallis_q_curve` (entropy) to understand
#' the relationship between entropy changes and divergence patterns.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme_minimal element_text
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stats IQR var
#'
#' @examples
#' # Create synthetic divergence data
#' set.seed(123)
#' divergence_matrix <- matrix(
#'   rnorm(80, mean = 0.5, sd = 0.1),
#'   nrow = 20, ncol = 4
#' )
#' rownames(divergence_matrix) <- paste0("gene_", 1:20)
#' colnames(divergence_matrix) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
#' divergence_se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(divergence = divergence_matrix)
#' )
#' 
#' # Global divergence curve (all genes aggregated)
#' p_global <- .plot_divergence_spectrum(divergence_se)
#' 
#' # Gene-specific divergence spectrum
#' p_gene <- .plot_divergence_spectrum(divergence_se, gene = "gene_1")
#'

#' @noRd

.plot_divergence_spectrum <- function(divergence_results_se,
                                     gene = NULL,
                                     lm_res = NULL,
                                     n_genes = 4,
                                     ncol = 2,
                                     metric = c("median", "mean"),
                                     variability_metric = c("iqr", "sd")) {
  # =========================================================================
  # INPUT VALIDATION (Bioconductor: Fail fast with clear messages)
  # =========================================================================
  
  # Validate divergence_results_se
  if (!inherits(divergence_results_se, "SummarizedExperiment")) {
    stop("'divergence_results_se' must be a SummarizedExperiment object",
         call. = FALSE)
  }
  
  # Extract and validate divergence matrix
  div_mat <- tryCatch({
    SummarizedExperiment::assay(divergence_results_se, 1)
  }, error = function(e) {
    stop("Failed to extract assay from 'divergence_results_se': ", e$message,
         call. = FALSE)
  })
  
  if (is.null(div_mat) || nrow(div_mat) == 0 || ncol(div_mat) == 0) {
    stop("'divergence_results_se' assay is empty. Expected at least 1 gene and 1 q-value",
         call. = FALSE)
  }
  
  # Match and validate metric parameter
  metric <- match.arg(metric)
  
  # Match and validate variability_metric parameter
  variability_metric <- match.arg(variability_metric)
  
  # Validate n_genes parameter
  if (!is.numeric(n_genes) || length(n_genes) != 1 || is.na(n_genes) || n_genes < 1) {
    stop("'n_genes' must be a positive integer", call. = FALSE)
  }
  n_genes <- as.integer(n_genes)
  
  # Validate ncol parameter
  if (!is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) || ncol < 1) {
    stop("'ncol' must be a positive integer", call. = FALSE)
  }
  ncol <- as.integer(ncol)
  
  # Validate gene parameter if provided
  if (!is.null(gene)) {
    if (!is.character(gene) || length(gene) != 1 || is.na(gene)) {
      stop("'gene' must be a single character string or NULL", call. = FALSE)
    }
  }
  
  # Validate lm_res parameter if provided
  if (!is.null(lm_res)) {
    if (!is.data.frame(lm_res)) {
      stop("'lm_res' must be a data.frame or NULL", call. = FALSE)
    }
    if (nrow(lm_res) == 0) {
      stop("'lm_res' data.frame is empty", call. = FALSE)
    }
  }
  
  # =========================================================================
  # HELPER FUNCTIONS
  # =========================================================================
  
  # Extract q-values from column names
  .extract_q_values <- function(col_names) {
    # Try to extract numeric values after "q_" or "q="
    extracted <- gsub("^q[_=]", "", col_names)
    q_vals <- as.numeric(extracted)
    
    # All must parse successfully
    if (any(is.na(q_vals))) {
      stop("Cannot extract numeric q-values from column names. ",
           "Expected format like 'q_0.5' or 'q=0.5'. ",
           "Got: ", paste(head(col_names, 3), collapse=", "),
           call. = FALSE)
    }
    
    return(q_vals)
  }
  
  # Get gene identifiers, validating existence
  .get_gene_identifiers <- function(se, div_mat) {
    rd <- SummarizedExperiment::rowData(se)
    
    # Prefer gene_name if available
    if (!is.null(rd) && "gene_name" %in% colnames(rd)) {
      return(rd$gene_name)
    }
    
    # Fall back to rownames
    rn <- rownames(div_mat)
    if (is.null(rn)) {
      stop("'divergence_results_se' has no gene identifiers in rowData('gene_name') ",
           "or rownames()", call. = FALSE)
    }
    
    return(rn)
  }
  
  # Find gene column in data.frame
  .find_gene_column <- function(df) {
    valid_cols <- c("gene", "gene_name", "gene_id")
    found <- valid_cols[valid_cols %in% colnames(df)]
    
    if (length(found) == 0) {
      stop("'lm_res' must have a column named 'gene', 'gene_name', or 'gene_id'",
           call. = FALSE)
    }
    
    return(found[1])
  }
  
  # Find p-value column in data.frame
  .find_pvalue_column <- function(df) {
    valid_cols <- c("adj_p_interaction", "p_interaction", "adj_p_value", "p_value")
    found <- valid_cols[valid_cols %in% colnames(df)]
    
    if (length(found) == 0) {
      stop("'lm_res' must have a p-value column: ",
           "adj_p_interaction, p_interaction, adj_p_value, or p_value",
           call. = FALSE)
    }
    
    return(found[1])
  }
  
  # =========================================================================
  # DATA PREPARATION
  # =========================================================================
  
  # Extract q-values from column names
  col_names <- colnames(div_mat)
  q_vals <- .extract_q_values(col_names)
  
  # Get gene identifiers
  gene_names <- .get_gene_identifiers(divergence_results_se, div_mat)
  
  # Validate gene identifiers
  if (length(gene_names) != nrow(div_mat)) {
    stop("Number of gene identifiers (", length(gene_names), ") does not match ",
         "number of rows in divergence matrix (", nrow(div_mat), ")",
         call. = FALSE)
  }
  
  # Sort by q-values for consistent ordering
  sort_idx <- order(q_vals)
  q_vals_sorted <- q_vals[sort_idx]
  div_mat_sorted <- div_mat[, sort_idx]
  
  # =========================================================================
  # CASE 1: GENE-SPECIFIC SPECTRUM (single gene)
  # =========================================================================
  # =========================================================================
  # CASE 1: GENE-SPECIFIC SPECTRUM (single gene)
  # =========================================================================
  
  if (!is.null(gene)) {
    # Validate gene exists in our data
    gene_idx <- which(gene_names == gene)[1]
    if (is.na(gene_idx)) {
      stop("Gene '", gene, "' not found. ",
           "Available genes: ", paste(head(gene_names, 5), collapse=", "),
           if (length(gene_names) > 5) paste0(", ... (", length(gene_names), " total)"),
           call. = FALSE)
    }
    
    # Extract divergence values for this gene
    gene_div <- as.numeric(div_mat_sorted[gene_idx, ])
    
    # Build plot data
    plot_df <- data.frame(
      q = q_vals_sorted,
      divergence = gene_div,
      stringsAsFactors = FALSE
    )
    
    # Create plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line(color = "#4575B4", linewidth = 1.2) +
      ggplot2::geom_point(color = "#4575B4", size = 3.5, alpha = 0.8) +
      ggplot2::labs(
        title = paste("Divergence Spectrum:", gene),
        x = "q value (diversity scale parameter)",
        y = "Tsallis Divergence D_q"
      ) +
      .tsenat_theme_base(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = .tsenat_font_sizes$title,
          face = "bold",
          hjust = 0.5
        )
      )
    
    return(p)
  }
  
  # =========================================================================
  # CASE 2: TOP N GENES SPECTRA (multi-gene faceted plot)
  # =========================================================================
  
  if (!is.null(lm_res)) {
    # Find gene and p-value columns
    gene_col <- .find_gene_column(lm_res)
    p_col <- .find_pvalue_column(lm_res)
    
    # Sort lm_res by p-value and get top genes
    lm_sorted <- lm_res[order(lm_res[[p_col]], na.last = TRUE), , drop = FALSE]
    top_genes_vec <- head(lm_sorted[[gene_col]], n_genes)
    
    # Validate that we found genes
    if (length(top_genes_vec) == 0) {
      stop("No genes found in 'lm_res', cannot create plot", call. = FALSE)
    }
    
    # Match top genes to divergence matrix
    gene_indices <- match(top_genes_vec, gene_names)
    
    # Filter out unmatched genes
    unmatched <- is.na(gene_indices)
    if (all(unmatched)) {
      stop("None of the top genes from 'lm_res' found in divergence matrix. ",
           "Check that gene identifiers match",
           call. = FALSE)
    }
    if (any(unmatched)) {
      warning("Some genes from 'lm_res' not found in divergence matrix. ",
              "Proceeding with ", sum(!unmatched), " matching genes",
              call. = FALSE)
      gene_indices <- gene_indices[!unmatched]
      top_genes_vec <- top_genes_vec[!unmatched]
    }
    
    # Build plot data for each gene
    plot_list <- list()
    
    for (i in seq_along(gene_indices)) {
      gene_idx <- gene_indices[i]
      gene_name_i <- gene_names[gene_idx]
      gene_div <- as.numeric(div_mat_sorted[gene_idx, ])
      p_val <- lm_sorted[[p_col]][i]
      
      plot_list[[i]] <- data.frame(
        q = q_vals_sorted,
        divergence = gene_div,
        gene = gene_name_i,
        p_value = p_val,
        stringsAsFactors = FALSE
      )
    }
    
    multi_gene_df <- do.call(rbind, plot_list)
    rownames(multi_gene_df) <- NULL
    
    # Extract confidence intervals if available
    rd <- SummarizedExperiment::rowData(divergence_results_se)
    ci_df <- NULL
    
    if (!is.null(rd)) {
      ci_data_list <- list()
      
      for (idx in seq_along(gene_indices)) {
        gene_idx <- gene_indices[idx]
        gene_name_i <- gene_names[gene_idx]
        
        # Extract CIs for this gene across all q-values
        ci_lower <- numeric(length(q_vals_sorted))
        ci_upper <- numeric(length(q_vals_sorted))
        
        for (j in seq_along(q_vals_sorted)) {
          q_val <- q_vals_sorted[j]
          lower_col <- paste0("lower_ci_q", q_val)
          upper_col <- paste0("upper_ci_q", q_val)
          
          if (lower_col %in% colnames(rd) && upper_col %in% colnames(rd)) {
            ci_lower[j] <- rd[[lower_col]][gene_idx]
            ci_upper[j] <- rd[[upper_col]][gene_idx]
          } else {
            ci_lower[j] <- NA_real_
            ci_upper[j] <- NA_real_
          }
        }
        
        ci_data_list[[idx]] <- data.frame(
          q = q_vals_sorted,
          lower = ci_lower,
          upper = ci_upper,
          gene = gene_name_i,
          stringsAsFactors = FALSE
        )
      }
      
      if (length(ci_data_list) > 0) {
        ci_df <- do.call(rbind, ci_data_list)
        rownames(ci_df) <- NULL
      }
    }
    
    # Sort genes by p-value for proper facet order
    gene_p_values <- multi_gene_df[!duplicated(multi_gene_df$gene), c("gene", "p_value")]
    gene_p_values <- gene_p_values[order(gene_p_values$p_value), ]
    gene_order <- gene_p_values$gene
    multi_gene_df$gene <- factor(multi_gene_df$gene, levels = gene_order)
    
    # Create faceted plot
    p <- ggplot2::ggplot(multi_gene_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::facet_wrap(~ gene, ncol = ncol, scales = "free_y")
    
    # Add CI ribbons if available
    if (!is.null(ci_df)) {
      ci_df$gene <- factor(ci_df$gene, levels = gene_order)
      p <- p + ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(x = q, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        alpha = 0.15,
        fill = "#4575B4",
        color = NA
      )
    }
    
    p <- p +
      ggplot2::geom_line(color = "#4575B4", linewidth = 1.2, alpha = 0.8) +
      ggplot2::geom_point(color = "#4575B4", size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "Divergence Spectra: Per-gene Comparisons",
        subtitle = paste0("Ranked by interaction significance (", metric, ")"),
        x = "q value (diversity scale parameter)",
        y = expression("Divergence D[q]")
      ) +
      .tsenat_theme_base(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = .tsenat_font_sizes$title,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          face = "italic",
          size = .tsenat_font_sizes$subtitle,
          hjust = 0.5
        ),
        panel.spacing = ggplot2::unit(1.5, "lines"),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = .tsenat_font_sizes$subtitle
        )
      )
    
    return(p)
  }
  
  # =========================================================================
  # CASE 3: GLOBAL DIVERGENCE CURVE (all genes aggregated)
  # =========================================================================
  # =========================================================================
  # CASE 3: GLOBAL DIVERGENCE CURVE (all genes aggregated)
  # =========================================================================
  
  # Build summary statistics for global curve
  if (variability_metric == "iqr") {
    summary_stats <- data.frame(
      q = q_vals_sorted,
      central = apply(div_mat_sorted, 2, function(x) {
        if (metric == "median") {
          median(x, na.rm = TRUE)
        } else {
          mean(x, na.rm = TRUE)
        }
      }),
      spread = apply(div_mat_sorted, 2, function(x) {
        stats::IQR(x, na.rm = TRUE)
      }),
      stringsAsFactors = FALSE
    )
    spread_factor <- 1 / 2  # IQR/2 for symmetric ribbon
    spread_label <- "IQR"
  } else {  # sd
    summary_stats <- data.frame(
      q = q_vals_sorted,
      central = apply(div_mat_sorted, 2, function(x) {
        if (metric == "median") {
          median(x, na.rm = TRUE)
        } else {
          mean(x, na.rm = TRUE)
        }
      }),
      spread = apply(div_mat_sorted, 2, function(x) {
        sqrt(stats::var(x, na.rm = TRUE))
      }),
      stringsAsFactors = FALSE
    )
    spread_factor <- 1
    spread_label <- "SD"
  }
  
  # Validate that spread values are not all NA
  if (all(is.na(summary_stats$spread))) {
    stop("Cannot compute variability metric (", variability_metric, "). ",
         "Check that divergence matrix contains valid numeric values",
         call. = FALSE)
  }
  
  if (all(is.na(summary_stats$central))) {
    stop("Cannot compute central tendency (", metric, "). ",
         "Check that divergence matrix contains valid numeric values",
         call. = FALSE)
  }
  
  metric_label <- if (metric == "median") "Median" else "Mean"
  
  # Create plot
  p <- ggplot2::ggplot(summary_stats, ggplot2::aes(x = q, y = central)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = central - spread * spread_factor,
                   ymax = central + spread * spread_factor),
      alpha = 0.1,
      fill = "#4575B4",
      color = NA
    ) +
    ggplot2::geom_line(color = "#4575B4", linewidth = 1.3) +
    ggplot2::geom_point(color = "#4575B4", size = 3.5, alpha = 0.8) +
    ggplot2::labs(
      title = expression("Global Divergence Spectrum: Average " * D[q] * " Across All Genes"),
      x = "q value (diversity scale parameter)",
      y = expression("Divergence D[q]"),
      subtitle = paste0(
        metric_label, " +/- ", spread_label,
        " (", nrow(div_mat_sorted), " genes)"
      )
    ) +
    .tsenat_theme_base(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = .tsenat_font_sizes$title,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        face = "italic",
        size = .tsenat_font_sizes$subtitle,
        hjust = 0.5
      )
    )
  
  return(p)
}
