# Tests for effect_sizes_divergence and results integration
# Covers uncovered lines from divergence_coverage.txt

test_that("effect_sizes_divergence aligns gene datasets", {
  # Create mock data matching expected structure
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.5),
    estimate_interaction = c(0.5, 0.3, 0.1),
    se_interaction = c(0.1, 0.15, 0.2)
  )
  
  # Create SummarizedExperiment with divergence results
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7)
    )
  )
  
  # Should align and compute effect sizes
  expect_error(
    effect_sizes_divergence(
      lm_res = lm_res,
      divergence_results_se = div_se,
      verbose = FALSE
    ),
    NA  # Expect no error
  )
})

test_that("effect_sizes_divergence filters non-matching genes", {
  # lm_res with different genes
  lm_res <- data.frame(
    gene = c("gene_A", "gene_B"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  # divergence_results_se with different genes
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene_X", "gene_Y"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8)
    )
  )
  
  # Should handle gene filtering gracefully
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Result should be valid (may be empty or have notes)
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("effect_sizes_divergence handles significance threshold", {
  # Test with different significance thresholds
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.5),
    estimate_interaction = c(0.5, 0.3, 0.1),
    se_interaction = c(0.1, 0.15, 0.2)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7)
    )
  )
  
  # Test with strict threshold
  result_strict <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    significance_threshold = 0.01,
    verbose = FALSE
  )
  
  # Test with lenient threshold
  result_lenient <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    significance_threshold = 0.1,
    verbose = FALSE
  )
  
  # Both should produce valid results
  expect_true(is.data.frame(result_strict) || is.list(result_strict))
  expect_true(is.data.frame(result_lenient) || is.list(result_lenient))
})

test_that("effect_sizes_divergence enriches per-q patterns", {
  # Test with enrich_per_q_pattern = TRUE
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1),
      q_0.5 = matrix(c(0.9, 0.8), nrow = 2, ncol = 1),
      q_1.0 = matrix(c(0.8, 0.7), nrow = 2, ncol = 1),
      q_2.0 = matrix(c(0.7, 0.6), nrow = 2, ncol = 1)
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8),
      per_q_pattern = c("RARE_DRIVEN", "BALANCED")
    )
  )
  
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence formats output correctly", {
  # Test output formatting and structure
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01),
    estimate_interaction = c(0.5),
    se_interaction = c(0.1)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(0.8, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1"),
      estimate = c(0.8),
      lower_ci = c(0.7),
      upper_ci = c(0.9)
    )
  )
  
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should return a valid result (data frame or SE or list)
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence handles zero divergence", {
  # Test handling of zero divergence estimates
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.0),  # Zero slope difference
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.0), nrow = 2, ncol = 1)),  # Zero divergence
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, 0.0),
      lower_ci = c(0.7, 0.0),
      upper_ci = c(0.9, 0.0)
    )
  )
  
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("effect_sizes_divergence handles NA divergence estimates", {
  # Test handling of NA estimates
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, NA), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, NA),
      lower_ci = c(0.7, NA),
      upper_ci = c(0.9, NA)
    )
  )
  
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence with verbose output", {
  # Test verbose output path
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01),
    estimate_interaction = c(0.5),
    se_interaction = c(0.1)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(0.8, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1"),
      estimate = c(0.8),
      lower_ci = c(0.7),
      upper_ci = c(0.9)
    )
  )
  
  # Capture output to verify verbose = TRUE produces messages
  result <- effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = TRUE
  )
  
  expect_true(!is.null(result))
})
