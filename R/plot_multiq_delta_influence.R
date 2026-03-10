#' Plot Multi-Q Delta Influence Heatmaps
#'
#' Creates combined heatmap panels showing delta influence (transcript switching magnitude)
#' across multiple q-values (diversity scales) for selected genes. Each heatmap shows
#' how transcript importance differs between conditions (delta_influence) across the
#' q-spectrum from 0.01 to 2.0.
#'
#' @param switching_results A multi-q switching analysis result from
#'   \code{\link{jackknife_isoform_switching}(q = c(...))}. Must include
#'   gene_ids, gene_name_map, and per-gene results for each q-value.
#' @param n_genes Numeric: number of top genes to visualize (default 4).
#'   When \code{lm_results} is provided, selects the n genes with lowest p-values.
#'   Otherwise, selects the first n genes from results.
#' @param lm_results Optional data.frame from \code{\link{calculate_lm_interaction}()}
#'   containing gene interaction statistics. Should have columns for gene identifiers
#'   ('gene_name' or 'gene_id') and p-values ('p_interaction' or 'adj_p_interaction').
#'   If provided, genes are ranked by p-value significance for selection of top genes.
#'
#' @return Character path to saved PNG file containing the combined heatmaps.
#'   The plot is automatically saved to a temporary file and can be displayed
#'   in R Markdown with \code{knitr::include_graphics()}.
#'
#' @details
#' **Heatmap interpretation:**
#'
#' - \bold{Rows}: Different q-values from 0.01 (rare isoforms) to 2.0 (dominant isoforms)
#' - \bold{Columns}: Individual transcripts of each gene
#' - \bold{Color scale}: Blue (negative delta_influence) = transcript more important in second condition;
#'   Red (positive delta_influence) = transcript more important in first condition;
#'   White = no switching effect
#' - \bold{Intensity}: Darker colors indicate stronger switching magnitude
#'
#' **Layout:**
#' - Multiple genes displayed in separate panels (up to 2 per row)
#' - Panels combined into single PNG for reproducible visualization
#' - Outliers (>95th percentile) are capped for better color contrast
#' - NaN and Inf values treated as missing (light gray)
#'
#' **Use case:**
#' Identify whether transcript switching is consistent across diversity scales
#' or scale-dependent. Genes with similar patterns across q-values show robust
#' isoform shifts; genes with varying patterns across q indicate q-dependent
#' switching driven by rare vs. abundant isoforms.
#'
#' @examples
#' \dontrun{
#' # After running multi-q jackknife analysis
#' heatmap_file <- plot_multiq_delta_influence_heatmaps(
#'   switching_results = multi_q_results,
#'   n_genes = 4
#' )
#' knitr::include_graphics(heatmap_file)
#' }
#'
#' @import grid
#' @import pheatmap
#' @importFrom grDevices png dev.off colorRampPalette
#' @export
plot_multiq_delta_influence_heatmaps <- function(
    switching_results,
  n_genes = 4,
  lm_results = NULL) {
  # Input validation
  if (!inherits(switching_results, "tsenat_isoform_switching_multiq")) {
    stop("switching_results must be a multi-q result from jackknife_isoform_switching()")
  }
  
  # Extract q-values from result names (q_0_01, q_0_50, etc.)
  q_result_keys <- names(switching_results)[grepl("^q_", names(switching_results))]
  if (length(q_result_keys) == 0) {
    stop("No multi-q results found in switching_results")
  }
  
  # Get the first result to extract gene information (all q values analyze same genes)
  first_q_key <- q_result_keys[1]
  first_result <- switching_results[[first_q_key]]
  
  # Extract gene IDs and names from first q-value result
  gene_ids <- first_result$gene_ids
  gene_name_map <- first_result$gene_name_map
  
  if (length(gene_ids) == 0) {
    stop("No genes found in switching_results")
  }
  
  # Select top N genes
  if (!is.null(lm_results)) {
    # Find p-value or adjusted p-value column for ranking
    p_col <- if ("adj_p_interaction" %in% colnames(lm_results)) {
      "adj_p_interaction"
    } else if ("p_interaction" %in% colnames(lm_results)) {
      "p_interaction"
    } else {
      NULL
    }
    
    # If we have a p-value column, rank genes by statistical significance
    if (!is.null(p_col)) {
      # Find which column in lm_results contains gene identifiers
      lm_col <- if ("gene_id" %in% colnames(lm_results)) {
        "gene_id"
      } else if ("gene_name" %in% colnames(lm_results)) {
        "gene_name"
      } else {
        NULL
      }
      
      # Match genes and sort by p-value
      if (!is.null(lm_col)) {
        matches <- match(gene_ids, lm_results[[lm_col]])
        p_values <- rep(Inf, length(gene_ids))
        matched_idx <- !is.na(matches)
        p_values[matched_idx] <- lm_results[[p_col]][matches[matched_idx]]
        
        # Sort by p-value (lowest p-values = most significant)
        gene_order <- order(p_values)
        gene_ids_sorted <- gene_ids[gene_order]
        top_genes_for_comparison <- gene_ids_sorted[1:min(n_genes, length(gene_ids_sorted))]
      } else {
        # Could not find gene identifier column, use first N genes
        top_genes_for_comparison <- gene_ids[1:min(n_genes, length(gene_ids))]
      }
    } else {
      # No p-value column found, use first N genes
      top_genes_for_comparison <- gene_ids[1:min(n_genes, length(gene_ids))]
    }
  } else {
    # No lm_results provided, use first N genes
    top_genes_for_comparison <- gene_ids[1:min(n_genes, length(gene_ids))]
  }
  
  # Prepare data for combined multi-Q heatmap
  all_gene_matrices <- list()
  all_gene_info <- list()
  
  for (gene_idx in seq_along(top_genes_for_comparison)) {
    gene_id <- top_genes_for_comparison[gene_idx]
    
    # Look up gene name from gene_name_map
    gene_name_idx <- which(gene_ids == gene_id)[1]
    gene_name <- NA_character_
    if (!is.na(gene_name_idx)) {
      gene_name <- gene_name_map[gene_name_idx]
    }
    
    # Collect delta_influence for all transcripts across all q-values
    heatmap_data <- NULL
    
    for (q_key in q_result_keys) {
      
      if (!is.null(switching_results[[q_key]]) && 
          !is.null(switching_results[[q_key]]$results_per_gene) &&
          gene_id %in% names(switching_results[[q_key]]$results_per_gene)) {
        gene_res <- switching_results[[q_key]]$results_per_gene[[gene_id]]
        if (!is.null(gene_res$delta_influence)) {
          # Extract q value from key: q_0_01 -> remove "q_" -> "0_01" -> replace "_" with "." -> "0.01"
          q_str_cleaned <- gsub("_", ".", gsub("^q_", "", q_key))
          q_num <- suppressWarnings(as.numeric(q_str_cleaned))
          col_name <- paste0("q_", sprintf("%.2f", q_num))
          delta_vals <- as.numeric(gene_res$delta_influence)
          
          # Replace Inf and NaN with NA for clean handling
          delta_vals[!is.finite(delta_vals)] <- NA
          
          if (is.null(heatmap_data)) {
            # Initialize with transcript IDs
            heatmap_data <- data.frame(
              transcript = as.character(gene_res$transcript_ids),
              stringsAsFactors = FALSE
            )
          }
          
          # Add column for this q-value (with Inf/NaN as NA)
          # Match number of rows and ensure alignment
          n_rows <- nrow(heatmap_data)
          if (length(delta_vals) < n_rows) {
            delta_vals <- c(delta_vals, rep(NA_real_, n_rows - length(delta_vals)))
          } else if (length(delta_vals) > n_rows) {
            delta_vals <- delta_vals[1:n_rows]
          }
          heatmap_data[[col_name]] <- as.numeric(delta_vals)
        }
      }
    }
    
    if (!is.null(heatmap_data) && nrow(heatmap_data) > 0 && ncol(heatmap_data) > 1) {
      # Convert to matrix for heatmap (transcripts as rows, q-values as columns)
      heatmap_matrix <- as.matrix(heatmap_data[, -1, drop = FALSE])
      rownames(heatmap_matrix) <- heatmap_data$transcript
      # Explicitly preserve column names (q-value labels)
      colnames(heatmap_matrix) <- colnames(heatmap_data)[-1]
      
      # Remove transcripts with exclusively no valid data (all NA)
      valid_transcript_rows <- rowSums(!is.na(heatmap_matrix)) > 0
      if (sum(valid_transcript_rows) > 0) {
        heatmap_matrix <- heatmap_matrix[valid_transcript_rows, , drop = FALSE]
        
        # Skip if no valid transcripts remain or fewer than 1 q-value column
        if (nrow(heatmap_matrix) > 0 && ncol(heatmap_matrix) > 0) {
          # Cap outliers for remaining finite values
          finite_vals <- heatmap_matrix[is.finite(heatmap_matrix)]
          if (length(finite_vals) > 0) {
            # Convert to numeric to avoid type coercion issues
            fin_abs <- abs(as.numeric(finite_vals))
            fin_abs <- fin_abs[is.finite(fin_abs)]
            if (length(fin_abs) > 0) {
              cap_val <- as.numeric(quantile(fin_abs, probs = 0.95))
              # Cap finite values that exceed the 95th percentile
              abs_hm <- abs(heatmap_matrix)
              mask <- which(is.finite(abs_hm) & abs_hm > cap_val)
              if (length(mask) > 0) {
                heatmap_matrix[mask] <- sign(heatmap_matrix[mask]) * cap_val
              }
            }
          }
          
          # Transpose: q-values as rows, transcripts as columns
          # This converts from (transcripts × q-values) to (q-values × transcripts)
          heatmap_matrix <- t(heatmap_matrix)
          
          # Verify we have the expected dimensions (q-values as rows)
          if (nrow(heatmap_matrix) > 0 && ncol(heatmap_matrix) > 0) {
            # Store matrix and gene info for combined plot
            all_gene_matrices[[gene_idx]] <- heatmap_matrix
            all_gene_info[[gene_idx]] <- list(
              gene_id = gene_id,
              gene_name = gene_name,
              n_transcripts = ncol(heatmap_matrix)
            )
          }
        }
      }
    }
  }
  
  # Create combined heatmap with genes in separate panels
  combined_png_file <- tempfile(pattern = "heatmap_multiQ_", fileext = ".png")
  
  if (length(all_gene_matrices) == 0) {
    warning("No valid heatmap data generated for any genes")
    return(NULL)
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package required for this function. Install with: install.packages('pheatmap')")
  }
  
  tryCatch({
    # Create individual heatmaps for each gene and store as grobs
    heatmap_plots <- list()
    plot_gene_names <- character(0)
    
    # Create plots only for genes that have data (skip NULL entries)
    for (gene_idx in seq_along(all_gene_matrices)) {
      if (is.null(all_gene_matrices[[gene_idx]]) || is.null(all_gene_info[[gene_idx]])) {
        next
      }
      
      mat <- all_gene_matrices[[gene_idx]]
      gene_info <- all_gene_info[[gene_idx]]
      gene_name <- gene_info$gene_name
      gene_id <- gene_info$gene_id
      plot_gene_names <- c(plot_gene_names, gene_name)
      
      # Construct header text for this gene
      if (is.na(gene_name) || gene_name == "") {
        header_text <- gene_id
      } else {
        header_text <- paste0(gene_name, " (", gene_id, ")")
      }
      
      # Apply outlier capping (safely handle if all values are NA)
      finite_vals <- as.numeric(mat[is.finite(mat)])
      if (length(finite_vals) == 0) {
        percentile_95 <- 1  # Default if no finite values
      } else {
        percentile_95 <- as.numeric(quantile(abs(finite_vals), 0.95))
      }
      
      mat_viz <- mat
      if (length(finite_vals) > 0) {
        # Cap outliers using indices to avoid coercion warnings
        abs_mat <- abs(mat_viz)
        mask_idx <- which(is.finite(abs_mat) & abs_mat > percentile_95)
        if (length(mask_idx) > 0) {
          mat_viz[mask_idx] <- sign(mat_viz[mask_idx]) * percentile_95
        }
      }
      
      # Create pheatmap (returns a grob object)
      p <- pheatmap::pheatmap(
        mat_viz,
        main = header_text,
        cluster_rows = FALSE,
        cluster_cols = (ncol(mat_viz) > 1),
        display_numbers = FALSE,
        na_col = "lightgray",
        color = grDevices::colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))(100),
        cellwidth = 65,
        cellheight = 65,
        fontsize = 16,
        fontsize_row = 16,
        fontsize_col = 16,
        fontsize_number = 13,
        margins = c(11, 10),
        show_rownames = TRUE,
        show_colnames = TRUE,
        silent = TRUE
      )
      
      heatmap_plots[[length(heatmap_plots) + 1]] <- p
    }
    
    # Combine all panels into one figure using manual grid layout
    # Save as PNG using manual grid layout
    grDevices::png(combined_png_file, width = 18, height = 9 * ceiling(length(heatmap_plots) / 2), 
                   units = "in", res = 96)
    
    grid::grid.newpage()
    
    # Add main title
    grid::grid.text("Delta Influence Across Diversity Scales", 
                    x = 0.5, y = 0.98, 
                    just = "top",
                    gp = grid::gpar(fontsize = 24, fontface = "bold"))
    
    # Create viewport layout
    n_genes <- length(heatmap_plots)
    n_cols <- 2
    n_rows <- ceiling(n_genes / n_cols)
    
    # Build layout with alternating content rows and spacing rows
    row_heights <- c()
    for (i in 1:n_rows) {
      row_heights <- c(row_heights, 1)  # content row
      if (i < n_rows) {
        row_heights <- c(row_heights, 0.55)  # spacing row
      }
    }
    n_layout_rows <- length(row_heights)
    
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.48, width = 1, height = 0.78,
                                      layout = grid::grid.layout(
      n_layout_rows, 
      n_cols, 
      heights = grid::unit(row_heights, "null"),
      widths = rep(1, n_cols),
      respect = FALSE
    )))
    
    # Draw each pheatmap in its own viewport
    plot_idx <- 1
    layout_row <- 1
    for (content_row in 1:n_rows) {
      for (col in 1:n_cols) {
        if (plot_idx <= length(heatmap_plots)) {
          grid::pushViewport(grid::viewport(layout.pos.row = layout_row, layout.pos.col = col))
          grid::grid.draw(heatmap_plots[[plot_idx]])
          grid::popViewport()
          plot_idx <- plot_idx + 1
        }
      }
      # Move to next content row (skip spacing row)
      layout_row <- layout_row + 2
    }
    
    grid::popViewport()
    grDevices::dev.off()
    
    Sys.sleep(0.5)
    
    return(combined_png_file)
    
  }, error = function(e) {
    tryCatch(grDevices::dev.off(), silent = TRUE)
    stop("Error creating heatmap: ", e$message)
  })
}
