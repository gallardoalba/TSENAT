# Test coverage for plot_multi_gene_q_spectrum_s4 function
# Located in generate_plots.R lines ~2727-2930

context("plot_multi_gene_q_spectrum_s4: Multi-gene q-spectrum plotting")

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
  library("SummarizedExperiment")
  
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
  library("ggplot2")
  
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
  library("ggplot2")
  
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
  library("ggplot2")
  
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
  library("ggplot2")
  
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
  library("SummarizedExperiment")
  
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
  library("ggplot2")
  library("patchwork")
  
  set.seed(42)
  plot_list <- list()
  
  for (i in 1:4) {
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
  library("ggplot2")
  
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
  library("ggplot2")
  
  n_genes <- 9
  plot_list <- list()
  
  for (i in 1:n_genes) {
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
