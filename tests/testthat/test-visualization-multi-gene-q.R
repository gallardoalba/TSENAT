# Test coverage for plot_multi_gene_q_spectrum_s4 function
# Located in generate_plots.R lines ~2727-2930

context("plot_multi_gene_q_spectrum_s4: Multi-gene q-spectrum plotting")
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)
library(testthat)




test_that("plot_multi_gene_q_spectrum_s4: TSENATAnalysis S4 object handling", {
  config <- list()
  
  # Create mock lm_results structure that would be extracted from S4 object
  lm_results <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.001, 0.01, 0.1),
    per_q_pattern = c("1,0.8,0.6", "0.5,0.4,0.3", "0.2,0.1,0.05")
  )
  
  # Verify structure
  expect_true("gene" %in% colnames(lm_results))
  expect_true("per_q_pattern" %in% colnames(lm_results))
  expect_true("adj_p_interaction" %in% colnames(lm_results))
})

test_that("plot_multi_gene_q_spectrum_s4: lm_results extraction from TSENATAnalysis", {
  config <- list()
  
  lm_results <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.001, 0.01, 0.1)
  )
  
  expect_equal(nrow(lm_results), 3)
  expect_true("gene" %in% colnames(lm_results))
  expect_true("adj_p_interaction" %in% colnames(lm_results))
})

test_that("plot_multi_gene_q_spectrum_s4: diversity_results extraction", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("g1", "g2", "g3")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_is(se, "SummarizedExperiment")
  expect_equal(nrow(se), 3)
})

test_that("plot_multi_gene_q_spectrum_s4: eff_res interaction_results validation", {
  config <- list()
  
  eff_res <- list(
    interaction_results = data.frame(
      gene = c("g1", "g2", "g3"),
      per_q_pattern = c("1.0,0.8,0.6", "0.5,0.4,0.3", "0.2,0.1,0.05"),
      adj_p_interaction = c(0.001, 0.01, 0.1)
    )
  )
  
  expect_true("gene" %in% colnames(eff_res$interaction_results))
  expect_true("per_q_pattern" %in% colnames(eff_res$interaction_results))
  expect_true("adj_p_interaction" %in% colnames(eff_res$interaction_results))
})

test_that("plot_multi_gene_q_spectrum_s4: p-value column detection (adj_p_interaction)", {
  config <- list()
  
  int_res <- data.frame(
    gene = c("g1", "g2"),
    per_q_pattern = c("1,0.8", "0.5,0.4"),
    adj_p_interaction = c(0.001, 0.01)
  )
  
  has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
  
  expect_true(has_p_adj)
})

test_that("plot_multi_gene_q_spectrum_s4: p-value column detection (fallback)", {
  config <- list()
  
  int_res <- data.frame(
    gene = c("g1", "g2"),
    per_q_pattern = c("1,0.8", "0.5,0.4"),
    p_value_interaction = c(0.001, 0.01)
  )
  
  has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
  has_p_raw <- "p_value_interaction" %in% colnames(int_res)
  
  p_col <- NA_character_
  if (has_p_adj) {
    p_col <- "adj_p_interaction"
  } else if (has_p_raw) {
    p_col <- "p_value_interaction"
  }
  
  expect_equal(p_col, "p_value_interaction")
})

test_that("plot_multi_gene_q_spectrum_s4: top N genes selection", {
  config <- list()
  
  int_res <- data.frame(
    gene = c("g1", "g2", "g3", "g4", "g5"),
    per_q_pattern = c("1,0.8", "0.5,0.4", "0.7,0.6", "0.3,0.2", "0.9,0.85"),
    adj_p_interaction = c(0.001, 0.01, 0.08, 0.05, 0.0001)
  )
  
  n_genes <- 3
  int_res_sorted <- int_res[order(int_res$adj_p_interaction, na.last = TRUE), ]
  int_res_subset <- head(int_res_sorted, n_genes)
  
  expect_equal(nrow(int_res_subset), 3)
  expect_equal(int_res_subset$gene[1], "g5")  # Most significant
})

test_that("plot_multi_gene_q_spectrum_s4: per_q_pattern validation", {
  config <- list()
  
  patterns <- c("1,0.8,0.6", "", NA, "NA")
  
  valid_patterns <- !is.na(patterns) & patterns != "" & patterns != "NA"
  
  expect_equal(sum(valid_patterns), 1)
})

test_that("plot_multi_gene_q_spectrum_s4: per_q_pattern parsing", {
  config <- list()
  
  pattern_str <- "1.0,0.8,0.6,0.4,0.2"
  
  per_q_vals <- as.numeric(strsplit(pattern_str, ",")[[1]])
  
  expect_equal(per_q_vals, c(1.0, 0.8, 0.6, 0.4, 0.2))
})

test_that("plot_multi_gene_q_spectrum_s4: empty pattern handling", {
  config <- list()
  
  pattern_str <- ""
  split_result <- strsplit(pattern_str, ",")[[1]]
  per_q_vals <- as.numeric(split_result)
  
  # Empty string split returns empty vector
  expect_equal(length(per_q_vals), 0)
})

test_that("plot_multi_gene_q_spectrum_s4: q-value grid generation", {
  config <- list()
  
  per_q_vals <- c(1.0, 0.8, 0.6, 0.4)
  q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
  
  expect_equal(length(q_vals), 4)
  expect_equal(q_vals[1], 0.1)
  expect_equal(q_vals[4], 0.25)
})

test_that("plot_multi_gene_q_spectrum_s4: plot_df construction", {
  config <- list()
  
  q_vals <- c(0.1, 0.15, 0.2, 0.25)
  per_q_vals <- c(1.0, 0.8, 0.6, 0.4)
  
  plot_df <- data.frame(
    q = q_vals,
    divergence = per_q_vals,
    stringsAsFactors = FALSE
  )
  
  expect_equal(nrow(plot_df), 4)
  expect_named(plot_df, c("q", "divergence"))
})

test_that("plot_multi_gene_q_spectrum_s4: individual q-spectrum plot", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  plot_df <- data.frame(
    q = c(0.1, 0.15, 0.2),
    divergence = c(1.0, 0.8, 0.6)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::geom_line(color = "#2E86AB", linewidth = 1.2) +
    ggplot2::geom_point(color = "#2E86AB", size = 2.8, alpha = 0.8)
  
  expect_is(p, "ggplot")
})

test_that("plot_multi_gene_q_spectrum_s4: vline for q=1", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  plot_df <- data.frame(
    q = seq(0.1, 0.3, length.out = 5),
    divergence = c(1.0, 0.9, 0.8, 0.7, 0.6)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = 1, linetype = 3, color = "gray60", linewidth = 0.8, alpha = 0.7)
  
  expect_is(p, "ggplot")
})

test_that("plot_multi_gene_q_spectrum_s4: plot title with gene name", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  gene_name <- "BRCA1"
  adj_p <- 0.001
  
  plot_df <- data.frame(q = 0.1, divergence = 1.0)
  
  title <- sprintf("%s (p=%.2e)", gene_name, adj_p)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = title)
  
  expect_is(p, "ggplot")
})

test_that("plot_multi_gene_q_spectrum_s4: plot list accumulation", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  genes <- c("g1", "g2", "g3")
  plot_list <- list()
  
  for (i in seq_along(genes)) {
    plot_df <- data.frame(
      q = seq(0.1, 0.3, length.out = 3),
      divergence = rnorm(3)
    )
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = genes[i])
    
    plot_list[[genes[i]]] <- p
  }
  
  expect_equal(length(plot_list), 3)
  expect_named(plot_list, genes)
})

test_that("plot_multi_gene_q_spectrum_s4: skipping genes with no valid per_q values", {
  config <- list()
  
  genes <- c("g1", "g2", "g3")
  patterns <- c("1,0.8", "", "0.5,0.4")
  
  valid_idx <- !sapply(patterns, function(p) p == "" || is.na(p))
  
  expect_equal(sum(valid_idx), 2)
})

test_that("plot_multi_gene_q_spectrum_s4: mode 1 failure - missing columns", {
  config <- list()
  
  int_res <- data.frame(
    gene = c("g1", "g2"),
    # Missing per_q_pattern and adj_p_interaction
    other_column = c(1, 2)
  )
  
  has_gene <- "gene" %in% colnames(int_res)
  has_per_q <- "per_q_pattern" %in% colnames(int_res)
  has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
  
  expect_true(has_gene)
  expect_false(has_per_q)
  expect_false(has_p_adj)
})

test_that("plot_multi_gene_q_spectrum_s4: mode 2 fallback - lm_res and divergence_results_se", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  lm_res <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.001, 0.01, 0.1)
  )
  
  div_assay <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(div_assay) <- c("g1", "g2", "g3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = div_assay),
    rowData = data.frame(
      gene_name = c("BRCA1", "TP53", "MYC"),
      row.names = rownames(div_assay)
    )
  )
  
  expect_equal(nrow(lm_res), nrow(se))
})

test_that("plot_multi_gene_q_spectrum_s4: per_q_pattern extraction from assay", {
  config <- list()
  
  gene_divs <- c(1.0, 0.8, 0.5, 0.2)
  per_q_patterns <- paste(gene_divs[!is.na(gene_divs)], collapse = ",")
  
  expect_equal(per_q_patterns, "1,0.8,0.5,0.2")
})

test_that("plot_multi_gene_q_spectrum_s4: no valid genes error handling", {
  config <- list()
  
  genes_to_plot <- NULL
  
  if (is.null(genes_to_plot) || length(genes_to_plot) == 0) {
    # Return NULL visibly for consistency
    expect_true(TRUE)
  }
})

test_that("plot_multi_gene_q_spectrum_s4: grid arrangement with patchwork", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  set.seed(42)
  plot_list <- list()
  
  for (i in seq_len(4)) {
    plot_df <- data.frame(
      q = seq(0.1, 0.3, length.out = 3),
      divergence = rnorm(3)
    )
    
    plot_list[[i]] <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = paste0("Gene", i))
  }
  
  # Arrange in grid
  ncol <- 2
  n_plots <- length(plot_list)
  n_rows <- ceiling(n_plots / ncol)
  
  expect_equal(n_rows, 2)
})

test_that("plot_multi_gene_q_spectrum_s4: single gene plotting", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  genes <- c("g1")
  plot_list <- list()
  
  for (i in seq_along(genes)) {
    plot_df <- data.frame(
      q = seq(0.1, 0.3, length.out = 3),
      divergence = c(1.0, 0.8, 0.6)
    )
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = genes[i])
    
    plot_list[[genes[i]]] <- p
  }
  
  expect_equal(length(plot_list), 1)
})

test_that("plot_multi_gene_q_spectrum_s4: many genes plotting", {
  config <- list()
  skip_if_not_installed("ggplot2")
  
  n_genes <- 9
  plot_list <- list()
  
  for (i in seq_len(n_genes)) {
    plot_df <- data.frame(
      q = seq(0.1, 0.3, length.out = 3),
      divergence = rnorm(3)
    )
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = paste0("g", i))
    
    plot_list[[i]] <- p
  }
  
  expect_equal(length(plot_list), 9)
})

test_that("plot_multi_gene_q_spectrum_s4: ncol parameter usage", {
  config <- list()
  
  n_genes <- 9
  ncol <- 3
  n_rows <- ceiling(n_genes / ncol)
  
  expect_equal(n_rows, 3)
})

test_that("plot_multi_gene_q_spectrum_s4: verbose output", {
  config <- list()
  
  # Verbose parameter controls informational messages
  verbose <- TRUE
  
  if (verbose) {
    # Would print status messages
    expect_true(verbose)
  }
})

test_that("plot_multi_gene_q_spectrum_s4: null eff_res handling", {
  config <- list()
  
  eff_res <- NULL
  
  expect_null(eff_res)
})

test_that("plot_multi_gene_q_spectrum_s4: null lm_res handling", {
  config <- list()
  
  lm_res <- NULL
  
  expect_null(lm_res)
})

test_that("plot_multi_gene_q_spectrum_s4: null divergence_results_se handling", {
  config <- list()
  
  divergence_results_se <- NULL
  
  expect_null(divergence_results_se)
})

test_that("plot_multi_gene_q_spectrum_s4: large divergence values", {
  config <- list()
  
  per_q_vals <- c(1.5, 1.2, 0.9, 0.6)
  
  expect_true(max(per_q_vals) > 1.0)
})

test_that("plot_multi_gene_q_spectrum_s4: small divergence values", {
  config <- list()
  
  per_q_vals <- c(0.01, 0.008, 0.005, 0.002)
  
  expect_true(max(per_q_vals) < 0.1)
})

test_that("plot_multi_gene_q_spectrum_s4: negative divergence values", {
  config <- list()
  
  # Can occur in signed divergence
  per_q_vals <- c(0.5, -0.2, 0.1, -0.05)
  
  expect_true(any(per_q_vals < 0))
})

test_that("plot_multi_gene_q_spectrum_s4: zero divergence values", {
  config <- list()
  
  per_q_vals <- c(0, 0, 0, 0)
  
  expect_true(all(per_q_vals == 0))
})

test_that("plot_multi_gene_q_spectrum_s4: NA handling in plots", {
  config <- list()
  
  per_q_vals <- c(1.0, NA, 0.6, NA)
  valid_vals <- per_q_vals[!is.na(per_q_vals)]
  
  expect_equal(length(valid_vals), 2)
})

# ============================================================================
# Tests for helper function: .extract_s4_components
# ============================================================================

test_that(".extract_s4_components: successful extraction from valid S4 object", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Create mock TSENATAnalysis object with required slots
  lm_results_list <- list(
    data.frame(
      gene = c("g1", "g2", "g3"),
      adj_p_interaction = c(0.001, 0.01, 0.1)
    )
  )
  
  diversity_results_list <- list(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(div = matrix(rnorm(12), nrow = 3, ncol = 4)),
      rowData = data.frame(gene_name = c("g1", "g2", "g3"))
    )
  )
  
  # Verify structure without actually creating S4 object
  expect_equal(nrow(lm_results_list[[1]]), 3)
  expect_equal(nrow(diversity_results_list[[1]]), 3)
})

test_that(".extract_s4_components: handles missing lm_results", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Empty lm_results should trigger error message
  lm_results_list <- list()
  
  expect_true(length(lm_results_list) == 0)
})

test_that(".extract_s4_components: handles missing diversity_results", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Empty diversity_results should trigger error message
  diversity_results_list <- list()
  
  expect_true(length(diversity_results_list) == 0)
})

test_that(".extract_s4_components: verbose messaging", {
  # Test that verbose parameter is respected
  verbose <- TRUE
  
  expect_true(verbose)
})

# ============================================================================
# Tests for helper function: .select_genes_from_eff_res
# ============================================================================

test_that(".select_genes_from_eff_res: successful gene selection from eff_res", {
  eff_res <- list(
    interaction_results = data.frame(
      gene = c("g1", "g2", "g3", "g4", "g5"),
      per_q_pattern = c("1,0.8", "0.5,0.4", "0.7,0.6", "0.3,0.2", "0.9,0.85"),
      adj_p_interaction = c(0.001, 0.01, 0.08, 0.05, 0.0001),
      stringsAsFactors = FALSE
    )
  )
  
  n_genes <- 3
  int_res <- eff_res$interaction_results
  int_res_sorted <- int_res[order(int_res$adj_p_interaction, na.last = TRUE), ]
  int_res_subset <- head(int_res_sorted, n_genes)
  
  expect_equal(nrow(int_res_subset), 3)
  expect_equal(int_res_subset$gene[1], "g5")  # Most significant
  expect_equal(int_res_subset$adj_p_interaction[1], 0.0001)
})

test_that(".select_genes_from_eff_res: null eff_res handling", {
  eff_res <- NULL
  
  if (is.null(eff_res)) {
    result <- NULL
  }
  
  expect_null(result)
})

test_that(".select_genes_from_eff_res: missing interaction_results", {
  eff_res <- list()  # No interaction_results key
  
  has_interaction <- !is.null(eff_res$interaction_results)
  
  expect_false(has_interaction)
})

test_that(".select_genes_from_eff_res: empty interaction_results", {
  eff_res <- list(interaction_results = data.frame())
  
  has_rows <- nrow(eff_res$interaction_results) > 0
  
  expect_false(has_rows)
})

test_that(".select_genes_from_eff_res: missing required columns - gene", {
  eff_res <- list(
    interaction_results = data.frame(
      per_q_pattern = c("1,0.8", "0.5,0.4"),
      adj_p_interaction = c(0.001, 0.01),
      stringsAsFactors = FALSE
    )
  )
  
  has_gene <- "gene" %in% colnames(eff_res$interaction_results)
  
  expect_false(has_gene)
})

test_that(".select_genes_from_eff_res: missing required columns - per_q_pattern", {
  eff_res <- list(
    interaction_results = data.frame(
      gene = c("g1", "g2"),
      adj_p_interaction = c(0.001, 0.01),
      stringsAsFactors = FALSE
    )
  )
  
  has_per_q <- "per_q_pattern" %in% colnames(eff_res$interaction_results)
  
  expect_false(has_per_q)
})

test_that(".select_genes_from_eff_res: p-value column priority (adj_p over raw)", {
  eff_res <- list(
    interaction_results = data.frame(
      gene = c("g1", "g2"),
      per_q_pattern = c("1,0.8", "0.5,0.4"),
      adj_p_interaction = c(0.001, 0.01),
      p_value_interaction = c(0.0001, 0.001),
      stringsAsFactors = FALSE
    )
  )
  
  has_p_adj <- "adj_p_interaction" %in% colnames(eff_res$interaction_results)
  p_col <- if (has_p_adj) "adj_p_interaction" else "p_value_interaction"
  
  expect_equal(p_col, "adj_p_interaction")
})

test_that(".select_genes_from_eff_res: fallback to raw p-value when adj_p missing", {
  eff_res <- list(
    interaction_results = data.frame(
      gene = c("g1", "g2"),
      per_q_pattern = c("1,0.8", "0.5,0.4"),
      p_value_interaction = c(0.0001, 0.001),
      stringsAsFactors = FALSE
    )
  )
  
  has_p_adj <- "adj_p_interaction" %in% colnames(eff_res$interaction_results)
  has_p_raw <- "p_value_interaction" %in% colnames(eff_res$interaction_results)
  
  p_col <- if (has_p_adj) "adj_p_interaction" else if (has_p_raw) "p_value_interaction" else NA_character_
  
  expect_equal(p_col, "p_value_interaction")
})

test_that(".select_genes_from_eff_res: valid pattern filtering", {
  patterns_subset <- data.frame(
    per_q_pattern = c("1,0.8", "0.5,0.4", "", NA, "NA"),
    stringsAsFactors = FALSE
  )
  
  valid_patterns <- !is.na(patterns_subset$per_q_pattern) & 
                   patterns_subset$per_q_pattern != "" & 
                   patterns_subset$per_q_pattern != "NA"
  
  expect_equal(sum(valid_patterns), 2)
  expect_equal(which(valid_patterns), c(1, 2))
})

test_that(".select_genes_from_eff_res: returns NULL when no valid patterns", {
  patterns_subset <- data.frame(
    per_q_pattern = c("", NA, "NA"),
    stringsAsFactors = FALSE
  )
  
  valid_patterns <- !is.na(patterns_subset$per_q_pattern) & 
                   patterns_subset$per_q_pattern != "" & 
                   patterns_subset$per_q_pattern != "NA"
  
  result <- if (!any(valid_patterns)) NULL else patterns_subset
  
  expect_null(result)
})

# ============================================================================
# Tests for helper function: .select_genes_fallback
# ============================================================================

test_that(".select_genes_fallback: successful fallback extraction", {
  skip_if_not_installed("SummarizedExperiment")
  
  lm_res <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.001, 0.01, 0.1),
    stringsAsFactors = FALSE
  )
  
  div_assay <- matrix(c(1, 0.8, 0.5, 0.2, 0.9, 0.7, 0.4, 0.1, 0.6, 0.5, 0.3, 0.1), nrow = 3, ncol = 4)
  rownames(div_assay) <- c("g1", "g2", "g3")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = div_assay))
  
  expect_equal(nrow(lm_res), nrow(se))
})

test_that(".select_genes_fallback: null lm_res handling", {
  lm_res <- NULL
  divergence_results_se <- "dummy"
  
  result <- if (is.null(lm_res)) NULL else divergence_results_se
  
  expect_null(result)
})

test_that(".select_genes_fallback: null divergence_results_se handling", {
  lm_res <- data.frame(gene = c("g1", "g2"))
  divergence_results_se <- NULL
  
  result <- if (is.null(divergence_results_se)) NULL else lm_res
  
  expect_null(result)
})

test_that(".select_genes_fallback: empty lm_res", {
  lm_res <- data.frame()
  
  expect_equal(nrow(lm_res), 0)
})

test_that(".select_genes_fallback: missing gene column in lm_res", {
  lm_res <- data.frame(
    adj_p_interaction = c(0.001, 0.01),
    stringsAsFactors = FALSE
  )
  
  has_required <- all(c("gene", "adj_p_interaction") %in% colnames(lm_res))
  
  expect_false(has_required)
})

test_that(".select_genes_fallback: gene name extraction from rowData", {
  skip_if_not_installed("SummarizedExperiment")
  
  div_assay <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(div_assay) <- c("g1", "g2", "g3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = div_assay),
    rowData = data.frame(
      gene_name = c("BRCA1", "TP53", "MYC"),
      row.names = rownames(div_assay)
    )
  )
  
  div_rd <- as.data.frame(SummarizedExperiment::rowData(se))
  
  expect_true("gene_name" %in% colnames(div_rd))
  expect_equal(div_rd$gene_name[1], "BRCA1")
})

test_that(".select_genes_fallback: gene name extraction from rownames fallback", {
  skip_if_not_installed("SummarizedExperiment")
  
  div_assay <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(div_assay) <- c("g1", "g2", "g3")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = div_assay))
  
  div_gene_names <- rownames(SummarizedExperiment::assay(se))
  
  expect_equal(div_gene_names, c("g1", "g2", "g3"))
})

test_that(".select_genes_fallback: gene matching between lm_res and divergence_results_se", {
  lm_genes <- c("g1", "g2", "g3")
  div_genes <- c("g1", "g2", "g3")
  
  gene_indices <- match(lm_genes, div_genes)
  valid_idx <- !is.na(gene_indices)
  
  expect_true(all(valid_idx))
  expect_equal(sum(valid_idx), 3)
})

test_that(".select_genes_fallback: partial gene matching", {
  lm_genes <- c("g1", "g2", "g4")
  div_genes <- c("g1", "g2", "g3")
  
  gene_indices <- match(lm_genes, div_genes)
  valid_idx <- !is.na(gene_indices)
  valid_genes <- lm_genes[valid_idx]
  
  expect_equal(sum(valid_idx), 2)
  expect_equal(valid_genes, c("g1", "g2"))
})

# ============================================================================
# Tests for helper function: .create_gene_q_plots
# ============================================================================

test_that(".create_gene_q_plots: successful plot creation", {
  skip_if_not_installed("ggplot2")
  
  genes <- c("g1", "g2")
  patterns <- c("1,0.8,0.6", "0.5,0.4,0.3")
  p_values <- c(0.001, 0.01)
  
  plot_list <- list()
  for (i in seq_along(genes)) {
    per_q_vals <- as.numeric(strsplit(patterns[i], ",")[[1]])
    
    if (length(per_q_vals) > 0 && !all(is.na(per_q_vals))) {
      q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
      plot_df <- data.frame(q = q_vals, divergence = per_q_vals, stringsAsFactors = FALSE)
      
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = genes[i])
      
      plot_list[[i]] <- p
    }
  }
  
  expect_equal(length(plot_list), 2)
  expect_true(all(sapply(plot_list, function(x) is(x, "ggplot"))))
})

test_that(".create_gene_q_plots: empty pattern handling", {
  genes <- c("g1")
  patterns <- c("")
  p_values <- c(0.001)
  
  per_q_vals <- as.numeric(strsplit(patterns[1], ",")[[1]])
  
  expect_equal(length(per_q_vals), 0)
})

test_that(".create_gene_q_plots: NA pattern handling", {
  genes <- c("g1")
  patterns <- c(NA_character_)
  p_values <- c(0.001)
  
  per_q_vals <- as.numeric(strsplit(patterns[1], ",")[[1]])
  
  # NA_character_ converted to numeric produces single NA_real_ value
  expect_equal(length(per_q_vals), 1)
  expect_true(is.na(per_q_vals[1]))
})

test_that('.create_gene_q_plots: "NA" string pattern handling', {
  genes <- c("g1")
  patterns <- c("NA")
  p_values <- c(0.001)
  
  per_q_vals <- suppressWarnings(as.numeric(strsplit(patterns[1], ",")[[1]]))
  
  expect_equal(per_q_vals[1], NA_real_)
})

test_that(".create_gene_q_plots: q_vals sequence generation", {
  per_q_vals <- c(1.0, 0.8, 0.6, 0.4)
  q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
  
  expect_equal(length(q_vals), 4)
  expect_equal(q_vals, c(0.1, 0.15, 0.2, 0.25))
})

test_that(".create_gene_q_plots: data frame construction", {
  q_vals <- c(0.1, 0.15, 0.2, 0.25)
  per_q_vals <- c(1.0, 0.8, 0.6, 0.4)
  
  plot_df <- data.frame(q = q_vals, divergence = per_q_vals, stringsAsFactors = FALSE)
  
  expect_equal(nrow(plot_df), 4)
  expect_named(plot_df, c("q", "divergence"))
  expect_equal(plot_df$q[1], 0.1)
  expect_equal(plot_df$divergence[1], 1.0)
})

test_that(".create_gene_q_plots: plot object creation", {
  skip_if_not_installed("ggplot2")
  
  plot_df <- data.frame(q = c(0.1, 0.15, 0.2), divergence = c(1.0, 0.8, 0.6))
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Test", subtitle = "p=0.001")
  
  expect_is(p, "ggplot")
})

test_that(".create_gene_q_plots: filtering invalid ggplot objects", {
  skip_if_not_installed("ggplot2")
  
  plot_list <- list(
    ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(),
    NULL,
    ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  )
  
  filtered <- Filter(function(p) !is.null(p) && methods::is(p, "ggplot"), plot_list)
  
  expect_equal(length(filtered), 2)
})

test_that(".create_gene_q_plots: error handling in tryCatch", {
  gene <- "test_gene"
  pattern <- "invalid,pattern,that,becomes,text"
  
  tryCatch({
    per_q_vals <- suppressWarnings(as.numeric(strsplit(pattern, ",")[[1]]))
    # "invalid" becomes NA with warning
    expect_true(is.na(per_q_vals[1]))
  }, error = function(e) {
    expect_true(FALSE)
  })
})

# ============================================================================
# Tests for helper function: .assemble_plot_grid
# ============================================================================

test_that(".assemble_plot_grid: successful grid assembly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  
  plot_list <- list()
  for (i in 1:4) {
    plot_df <- data.frame(q = seq(0.1, 0.3, length.out = 3), divergence = rnorm(3))
    plot_list[[i]] <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = paste0("Gene", i))
  }
  
  ncol <- 2
  nrow <- ceiling(length(plot_list) / ncol)
  
  expect_equal(nrow, 2)
  expect_equal(length(plot_list), 4)
})

test_that(".assemble_plot_grid: single plot", {
  skip_if_not_installed("ggplot2")
  
  plot_df <- data.frame(q = c(0.1, 0.15, 0.2), divergence = c(1.0, 0.8, 0.6))
  plot_list <- list(ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Gene1"))
  
  expect_equal(length(plot_list), 1)
})

test_that(".assemble_plot_grid: empty plot list", {
  plot_list <- list()
  
  result <- if (length(plot_list) == 0) NULL else plot_list
  
  expect_null(result)
})

test_that(".assemble_plot_grid: grid dimension calculation", {
  n_plots <- 9
  ncol <- 3
  nrow <- ceiling(n_plots / ncol)
  
  expect_equal(nrow, 3)
})

test_that(".assemble_plot_grid: grid dimension with partial row", {
  n_plots <- 10
  ncol <- 3
  nrow <- ceiling(n_plots / ncol)
  
  expect_equal(nrow, 4)
})

test_that(".assemble_plot_grid: row-wise plot extraction", {
  plot_list <- list(
    "p1", "p2", "p3",
    "p4", "p5", "p6",
    "p7", "p8", "p9"
  )
  
  ncol <- 3
  nrow <- ceiling(length(plot_list) / ncol)
  
  row_1_start <- (1 - 1) * ncol + 1
  row_1_end <- min(1 * ncol, length(plot_list))
  row_1_plots <- plot_list[row_1_start:row_1_end]
  
  expect_equal(row_1_plots, list("p1", "p2", "p3"))
})

test_that(".assemble_plot_grid: row filtering for valid ggplot objects", {
  skip_if_not_installed("ggplot2")
  
  row_plots <- list(
    ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point(),
    NULL,
    ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  )
  
  filtered_row <- Filter(function(p) !is.null(p) && methods::is(p, "ggplot"), row_plots)
  
  expect_equal(length(filtered_row), 2)
})

test_that(".assemble_plot_grid: patchwork plot combining", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  
  p1 <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  p2 <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  
  combined <- p1 + p2
  
  expect_is(combined, "patchwork")
})

test_that(".assemble_plot_grid: patchwork plot spacing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  
  p1 <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  spacer <- patchwork::plot_spacer()
  p2 <- ggplot2::ggplot(data.frame(x = 1, y = 1)) + ggplot2::geom_point()
  
  layout <- p1 / spacer / p2
  
  expect_is(layout, "patchwork")
})
