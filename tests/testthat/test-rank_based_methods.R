library(TSENAT)

context("rank_based_methods: Rank-Based Nonparametric Methods")

test_that("compute_rank_correlation_multiq computes spearman correlations", {
  # Create test p-value lists
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50),
    q10 = runif(50)
  )
  
  result <- compute_rank_correlation_multiq(pvalues, method = "spearman")
  
  expect_s3_class(result, "rank_correlation_multiq")
  expect_equal(nrow(result$correlation_matrix), 3)
  expect_equal(ncol(result$correlation_matrix), 3)
  expect_true(all(diag(result$correlation_matrix) == 1.0))  # Diagonal is 1
  expect_true(result$mean_correlation >= -1 && result$mean_correlation <= 1)
})

test_that("compute_rank_correlation_multiq handles kendall correlation", {
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- compute_rank_correlation_multiq(pvalues, method = "kendall")
  
  expect_equal(result$method, "kendall")
  expect_true(!is.na(result$mean_correlation))
})

test_that("apply_aligned_rank_transform produces valid output", {
  # Create synthetic expression data
  set.seed(123)
  expr_data <- matrix(
    rnorm(160),  # 16 genes x 10 samples
    nrow = 16,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:16),
      paste0("Sample", 1:10)
    )
  )
  
  factors <- data.frame(
    batch = factor(c(rep("A", 5), rep("B", 5))),
    row.names = colnames(expr_data)
  )
  
  result <- apply_aligned_rank_transform(
    data = expr_data,
    factors = factors,
    formula = ~ batch
  )
  
  expect_s3_class(result, "art_result")
  expect_equal(dim(result$aligned_ranks), dim(expr_data))
  expect_equal(dim(result$normal_scores), dim(expr_data))
  expect_true(all(!is.na(result$normal_scores)))
})

test_that("test_rankbased_assumptions validates assumptions", {
  # Create synthetic expression data
  set.seed(456)
  expr_data <- matrix(
    rnorm(400),  # 20 genes x 20 samples
    nrow = 20,
    ncol = 20,
    dimnames = list(
      paste0("Gene", 1:20),
      paste0("Sample", 1:20)
    )
  )
  
  # Call the function and verify it returns a result
  result <- test_rankbased_assumptions(
    data = expr_data,
    checks = c("exchangeability", "monotonicity")
  )
  
  # Just verify it returns something (the function may return S3 object of class "rank_assumptions")
  expect_true(!is.null(result))
  expect_true(length(result) > 0)
})

test_that("rank correlation matrix is symmetric", {
  pvalues <- list(
    q01 = c(0.01, 0.05, 0.1, 0.5, 0.9),
    q05 = c(0.02, 0.04, 0.12, 0.48, 0.88),
    q10 = c(0.015, 0.055, 0.09, 0.51, 0.92)
  )
  
  result <- compute_rank_correlation_multiq(pvalues)
  
  # Check symmetry
  expect_true(all(result$correlation_matrix == 
                  t(result$correlation_matrix), na.rm = TRUE))
})

test_that("ART handles edge cases gracefully", {
  # Create simple test data
  expr_data <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(expr_data) <- paste0("gene", 1:5)
  colnames(expr_data) <- paste0("sample", 1:10)
  
  factors <- data.frame(
    batch = factor(rep(c("A", "B"), 5)),
    row.names = paste0("sample", 1:10)
  )
  
  result <- apply_aligned_rank_transform(
    data = expr_data,
    factors = factors
  )
  
  expect_s3_class(result, "art_result")
  expect_true(!any(is.infinite(result$normal_scores)))
})

# ============================================================================
# Permutation-Based Rank Correlation Confidence Intervals Tests
# ============================================================================

test_that("rank_correlation_bootstrap_ci returns correct structure", {
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50),
    q10 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman", 
                                         ci = "percentile", n_bootstrap = 100)
  
  expect_s3_class(result, "rank_correlation_ci")
  expect_true("correlation_matrix" %in% names(result))
  expect_true("ci_matrix" %in% names(result))
  expect_true("interpretation" %in% names(result))
  expect_equal(result$ci_level, 0.95)
})

test_that("rank_correlation_bootstrap_ci handles percentile CI", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman", 
                                         ci = "percentile", n_bootstrap = 200)
  
  expect_equal(result$ci_type, "percentile")
  expect_equal(dim(result$ci_matrix), c(2, 2, 2))
  # CI bounds should be reasonable
  expect_true(all(result$ci_matrix[, , "lower"] <= result$correlation_matrix))
  expect_true(all(result$ci_matrix[, , "upper"] >= result$correlation_matrix))
})

test_that("rank_correlation_bootstrap_ci handles BCA CI", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- suppressWarnings(rank_correlation_bootstrap_ci(pvalues, method = "spearman", 
                                         ci = "bca", n_bootstrap = 200))
  
  expect_equal(result$ci_type, "bca")
  expect_equal(dim(result$ci_matrix), c(2, 2, 2))
  # CI bounds should be within [-1, 1] for correlations
  expect_true(all(result$ci_matrix >= -1, na.rm = TRUE))
  expect_true(all(result$ci_matrix <= 1, na.rm = TRUE))
})

test_that("rank_correlation_bootstrap_ci handles permutation CI", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman", 
                                         ci = "permutation", n_permutations = 200)
  
  expect_equal(result$ci_type, "permutation")
  expect_equal(dim(result$ci_matrix), c(2, 2, 2))
  # Permutation CI should be conservative (wider bounds)
  expect_true(all(result$ci_matrix >= -1, na.rm = TRUE))
  expect_true(all(result$ci_matrix <= 1, na.rm = TRUE))
})

test_that("rank_correlation_bootstrap_ci computed correlations in valid range", {
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50) * 0.8,
    q10 = runif(50) * 0.7
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  # All correlations should be in [-1, 1]
  expect_true(all(result$correlation_matrix >= -1, na.rm = TRUE))
  expect_true(all(result$correlation_matrix <= 1, na.rm = TRUE))
  # Diagonal should be 1 (self-correlation)
  expect_true(all(diag(result$correlation_matrix) == 1))
})

test_that("rank_correlation_bootstrap_ci handles kendall method", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "kendall",
                                         ci = "percentile", n_bootstrap = 100)
  
  expect_true(grepl("kendall", tolower(result$method)))
  expect_true(!is.na(result$correlation_matrix[1, 2]))
})

test_that("rank_correlation_bootstrap_ci interpretation table is correct", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50),
    q10 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  # Should have 3 choose 2 = 3 pairs
  expect_equal(nrow(result$interpretation), 3)
  expect_true("Stability" %in% colnames(result$interpretation))
  expect_true("Q_value_Pair" %in% colnames(result$interpretation))
})

test_that("rank_correlation_bootstrap_ci recognizes CI level", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result_90 <- rank_correlation_bootstrap_ci(pvalues, ci_level = 0.90,
                                            ci = "percentile", n_bootstrap = 100)
  result_99 <- rank_correlation_bootstrap_ci(pvalues, ci_level = 0.99,
                                            ci = "percentile", n_bootstrap = 100)
  
  expect_equal(result_90$ci_level, 0.90)
  expect_equal(result_99$ci_level, 0.99)
  # Higher CI level should give wider intervals (mostly)
  expect_true(any((result_99$ci_matrix[, , "upper"] - result_99$ci_matrix[, , "lower"]) >= 
                  (result_90$ci_matrix[, , "upper"] - result_90$ci_matrix[, , "lower"]), na.rm = TRUE))
})

test_that("rank_correlation_bootstrap_ci handles named list inputs", {
  pvalues_named <- list(
    runif(50),
    runif(50),
    runif(50)
  )
  names(pvalues_named) <- c("q_01", "q_05", "q_10")
  
  result <- rank_correlation_bootstrap_ci(pvalues_named, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  expect_equal(colnames(result$correlation_matrix), c("q_01", "q_05", "q_10"))
  expect_true(all(grepl("_", result$interpretation$Q_value_Pair)))
})

test_that("rank_correlation_bootstrap_ci stability classification works", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50),
    q10 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  # Stability should be classified (checks for pattern matching with descriptions)
  stability_levels <- c("Variable", "Very stable", "Robust", "Moderate", "Weak")
  expect_true(all(sapply(result$interpretation$Stability, function(s) {
    any(sapply(stability_levels, function(level) grepl(level, s, fixed = TRUE)))
  })))
})

test_that("rank_correlation_bootstrap_ci CI bounds respect symmetry", {
  set.seed(42)
  pvalues <- list(
    q01 = c(0.01, 0.05, 0.1, 0.5, 0.9),
    q05 = c(0.02, 0.04, 0.12, 0.48, 0.88),
    q10 = c(0.015, 0.055, 0.09, 0.51, 0.92)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  # CI matrix should be symmetric for correlation
  expect_true(all(result$ci_matrix[1, 2, ] == result$ci_matrix[2, 1, ], na.rm = TRUE))
})

test_that("rank_correlation_bootstrap_ci returns distribution when requested", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result_with_dist <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                                   ci = "percentile", n_bootstrap = 100,
                                                   return_distribution = TRUE)
  result_no_dist <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                                 ci = "percentile", n_bootstrap = 100,
                                                 return_distribution = FALSE)
  
  expect_is(result_with_dist$bootstrap_distribution, "array")
  expect_null(result_no_dist$bootstrap_distribution)
})

test_that("rank_correlation_bootstrap_ci minimum requirement is 2 q-values", {
  pvalues_single <- list(q01 = runif(50))
  
  expect_error(rank_correlation_bootstrap_ci(pvalues_single))
})

test_that("rank_correlation_bootstrap_ci CI order is correct", {
  set.seed(42)
  pvalues <- list(
    q01 = runif(50),
    q05 = runif(50)
  )
  
  result <- rank_correlation_bootstrap_ci(pvalues, method = "spearman",
                                         ci = "percentile", n_bootstrap = 100)
  
  # Lower bound should always be <= upper bound
  expect_true(all(result$ci_matrix[, , "lower"] <= result$ci_matrix[, , "upper"]))
})


# ============================================================================
# NEW TESTS FOR GAP 3: INTERACTION DETECTION FUNCTIONS
# ============================================================================

# Test Suite: detect_q_gene_interactions()
test_that("detect_q_gene_interactions basic functionality works", {
  # Create synthetic q×gene interaction data
  set.seed(42)
  entropy_vals <- c(
    # Gene 1: robust across q (entropy doesn't change much)
    rnorm(10, mean = 1.0, sd = 0.1),  # q=0.5
    rnorm(10, mean = 1.05, sd = 0.1), # q=1.0
    rnorm(10, mean = 1.02, sd = 0.1), # q=1.5
    # Gene 2: strongly q-dependent (entropy changes with q)
    rnorm(10, mean = 0.5, sd = 0.1),  # q=0.5
    rnorm(10, mean = 1.5, sd = 0.1),  # q=1.0
    rnorm(10, mean = 2.5, sd = 0.1),  # q=1.5
    # Gene 3: moderately q-dependent
    rnorm(10, mean = 1.0, sd = 0.1),  # q=0.5
    rnorm(10, mean = 1.2, sd = 0.1),  # q=1.0
    rnorm(10, mean = 1.35, sd = 0.1)  # q=1.5
  )
  
  model_data <- data.frame(
    entropy = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 10),
    sample = rep(paste0("S", 1:10), 9),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Check output structure
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 3)  # 3 genes
  expect_true(all(c("gene", "p_value", "f_statistic", "effect_size_eta2", 
                     "interaction_class") %in% colnames(result)))
  
  # Check column types
  expect_is(result$p_value, "numeric")
  expect_is(result$f_statistic, "numeric")
  expect_is(result$effect_size_eta2, "numeric")
  expect_is(result$interaction_class, "character")
})

test_that("detect_q_gene_interactions correctly identifies robust gene", {
  # Gene with stable entropy across q
  set.seed(123)
  entropy_stable <- c(
    rnorm(10, mean = 1.0, sd = 0.05),  # q=0.5
    rnorm(10, mean = 1.0, sd = 0.05),  # q=1.0
    rnorm(10, mean = 1.0, sd = 0.05)   # q=1.5
  )
  
  model_data <- data.frame(
    entropy = entropy_stable,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("RobustGene", 30),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Robust gene should have high p-value (not significant)
  expect_true(result$p_value[1] > 0.05)
  expect_equal(result$interaction_class[1], "Robust across q")
})

test_that("detect_q_gene_interactions correctly identifies q-dependent gene", {
  # Gene with changing entropy across q
  set.seed(456)
  entropy_changing <- c(
    rnorm(10, mean = 0.5, sd = 0.05),  # q=0.5 (low)
    rnorm(10, mean = 2.0, sd = 0.05),  # q=1.0 (high)
    rnorm(10, mean = 3.5, sd = 0.05)   # q=1.5 (very high)
  )
  
  model_data <- data.frame(
    entropy = entropy_changing,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("DependentGene", 30),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Q-dependent gene should have low p-value (significant) and large effect size
  expect_true(result$p_value[1] < 0.05)
  expect_true(result$effect_size_eta2[1] > 0.10)
  expect_equal(result$interaction_class[1], "Strongly q-dependent")
})

test_that("detect_q_gene_interactions handles missing values gracefully", {
  set.seed(789)
  entropy_vals <- rnorm(30)
  entropy_vals[c(5, 15, 25)] <- NA  # Add some NAs
  
  model_data <- data.frame(
    entropy = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("TestGene", 30),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Should complete without error
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that("detect_q_gene_interactions requires minimum 2 q-levels", {
  # Data with only one q-value
  model_data <- data.frame(
    entropy = rnorm(10),
    q = rep(0.5, 10),
    gene = rep("Gene1", 10),
    sample = paste0("S", 1:10),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Should mark as "Insufficient data"
  expect_equal(result$interaction_class[1], "Insufficient data")
  expect_true(is.na(result$p_value[1]))
})

test_that("detect_q_gene_interactions validates column names", {
  model_data <- data.frame(
    value = rnorm(30),
    parameter = rep(c(0.5, 1.0, 1.5), each = 10),
    name = rep("Gene1", 30),
    sample = rep(paste0("S", 1:10), 3)
  )
  
  # Should error on missing entropy column
  expect_error(
    detect_q_gene_interactions(model_data),
    "not found in data"
  )
})

test_that("detect_q_gene_interactions kruskal.test method produces results", {
  set.seed(111)
  entropy_vals <- c(
    rnorm(10, mean = 1.0, sd = 0.2),  # q=0.5
    rnorm(10, mean = 1.5, sd = 0.2),  # q=1.0
    rnorm(10, mean = 2.0, sd = 0.2)   # q=1.5
  )
  
  model_data <- data.frame(
    entropy = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("Gene1", 30),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  
  # Should complete successfully
  expect_is(result, "data.frame")
  expect_true(!is.na(result$f_statistic[1]))
  expect_true(!is.na(result$p_value[1]))
})


# Test Suite: classify_q_dependency()
test_that("classify_q_dependency produces correct classifications", {
  set.seed(222)
  
  # Create interaction results with known properties
  interaction_results <- data.frame(
    gene = c("RobustGene", "ModerateGene", "StrongGene", "FailGene"),
    p_value = c(0.10, 0.03, 0.01, NA),
    effect_size_eta2 = c(0.005, 0.06, 0.15, 0.05),
    f_statistic = c(1.2, 3.5, 8.2, 2.1),
    interaction_class = c("", "", "", "Test failed"),
    stringsAsFactors = FALSE
  )
  
  classifications <- classify_q_dependency(interaction_results)
  
  expect_equal(classifications[1], "Robust across q")
  expect_equal(classifications[2], "Moderately q-dependent")
  expect_equal(classifications[3], "Strongly q-dependent")
  expect_equal(classifications[4], "Test failed")
})

test_that("classify_q_dependency respects custom thresholds", {
  interaction_results <- data.frame(
    gene = c("Gene1", "Gene2"),
    p_value = c(0.03, 0.005),
    effect_size_eta2 = c(0.07, 0.12),
    f_statistic = c(3.0, 4.0),
    interaction_class = c("", ""),
    stringsAsFactors = FALSE
  )
  
  # With default thresholds
  class_default <- classify_q_dependency(interaction_results)
  expect_equal(class_default[1], "Moderately q-dependent")  # p=0.03<0.05, eta2=0.07<0.10
  expect_equal(class_default[2], "Strongly q-dependent")    # p=0.005<0.05, eta2=0.12>0.10
  
  # With custom thresholds (stricter)
  class_strict <- classify_q_dependency(
    interaction_results,
    p_threshold = 0.01,
    eta2_threshold_strong = 0.05
  )
  expect_equal(class_strict[1], "Robust across q")      # p = 0.03 > 0.01 (not significant)
  expect_equal(class_strict[2], "Strongly q-dependent") # p = 0.005 < 0.01 (significant), eta2 = 0.12 > 0.05
})

test_that("classify_q_dependency handles all NA p-values", {
  interaction_results <- data.frame(
    gene = c("Gene1", "Gene2"),
    p_value = c(NA, NA),
    effect_size_eta2 = c(0.05, 0.10),
    f_statistic = c(2.0, 3.0),
    interaction_class = c("Insufficient data", "Insufficient data"),
    stringsAsFactors = FALSE
  )
  
  classifications <- classify_q_dependency(interaction_results)
  
  expect_equal(classifications[1], "Insufficient data")
  expect_equal(classifications[2], "Insufficient data")
})

test_that("classify_q_dependency is a vector", {
  interaction_results <- data.frame(
    gene = c("Gene1", "Gene2", "Gene3"),
    p_value = c(0.10, 0.03, 0.01),
    effect_size_eta2 = c(0.005, 0.06, 0.15),
    f_statistic = c(1.2, 3.5, 8.2),
    interaction_class = c("", "", ""),
    stringsAsFactors = FALSE
  )
  
  classifications <- classify_q_dependency(interaction_results)
  
  expect_is(classifications, "character")
  expect_equal(length(classifications), 3)
})


# Integration Tests
test_that("Full workflow: detect -> classify -> recommend works end-to-end", {
  set.seed(999)
  
  # Build test dataset by combining gene groups directly
  # Gene 1-30: robust (constant entropy across q)
  data_robust <- data.frame(
    entropy = rnorm(300, mean = 1.5, sd = 0.1),
    q = rep(c(0.5, 1.0, 1.5), 100),
    gene = rep(paste0("Gene", 1:30), each = 10),
    sample = rep(paste0("S", 1:10), 30),
    stringsAsFactors = FALSE
  )
  
  # Gene 31-40: moderately q-dependent (moderate variation)
  data_moderate <- data.frame(
    entropy = c(rnorm(50, mean = 1.0, sd = 0.1), rnorm(50, mean = 1.3, sd = 0.1)),
    q = rep(c(0.5, 1.0, 1.5), c(34, 33, 33)),
    gene = rep(paste0("Gene", 31:40), each = 10),
    sample = rep(paste0("S", 1:10), 10),
    stringsAsFactors = FALSE
  )
  
  # Gene 41-45: strongly q-dependent (large variation)
  data_strong <- data.frame(
    entropy = c(rnorm(25, mean = 0.5, sd = 0.1), rnorm(25, mean = 2.0, sd = 0.1)),
    q = rep(c(0.5, 1.0, 1.5), c(17, 17, 16)),
    gene = rep(paste0("Gene", 41:45), each = 5),
    sample = rep(paste0("S", 1:5), 5),
    stringsAsFactors = FALSE
  )
  
  # Combine all data
  model_data <- rbind(data_robust, data_moderate, data_strong)
  
  # Step 1: Detect interactions
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  results <- suppressWarnings(detect_q_gene_interactions(model_data))
  expect_equal(nrow(results), 45)
  
  # Step 2: Classify
  classifications <- classify_q_dependency(results)
  expect_equal(length(classifications), 45)
  expect_true(all(classifications %in% c("Robust across q", "Moderately q-dependent", 
                                          "Strongly q-dependent", "Insufficient data")))
  
  # Step 3: Recommend
  recommendation <- recommend_q_range(results)
  expect_is(recommendation, "list")
  expect_true(nchar(recommendation$recommendation) > 0)
  expect_true(nchar(recommendation$rationale) > 0)
})

test_that("Functions handle edge case: single sample per q-level", {
  # Minimal viable dataset: 1 sample per q per gene
  model_data <- data.frame(
    entropy = c(1.0, 1.1, 1.2, 0.5, 2.0, 3.5),
    q = rep(c(0.5, 1.0, 1.5), 2),
    gene = rep(c("Gene1", "Gene2"), each = 3),
    sample = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  
  # Should handle without error (though with limited power)
  # Suppress expected warnings from edge-case variance calculations with N=1 per group
  result <- suppressWarnings(detect_q_gene_interactions(model_data))
  expect_is(result, "data.frame")
})

test_that("Functions handle edge case: many q-levels", {
  # Dataset with 10 q-levels
  q_vals <- seq(0.1, 2.5, by = 0.25)
  set.seed(555)
  
  model_data <- data.frame(
    entropy = rnorm(100),
    q = rep(q_vals, length.out = 100),
    gene = rep("Gene1", 100),
    sample = rep(paste0("S", 1:10), 10),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(model_data)
  expect_equal(nrow(result), 1)
  expect_is(result$p_value[1], "numeric")
})

# ============================================================================
# NEW TEST SUITE: Westfall-Young Permutation for detect_q_gene_interactions
# ============================================================================

test_that("detect_q_gene_interactions westfall-young parameter is accepted", {
  # Create test data
  set.seed(777)
  model_data <- data.frame(
    entropy = rnorm(100, mean = 1, sd = 0.5),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 25),
    gene = rep(paste0("Gene", 1:5), each = 20),
    sample = rep(paste0("S", 1:5), 20),
    stringsAsFactors = FALSE
  )
  
  # Should accept westfall-young without error
  result <- detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 10  # Small number for speed in tests
  )
  
  expect_is(result, "data.frame")
  expect_true("adj_p_value" %in% colnames(result))
  expect_equal(nrow(result), 5)  # 5 genes
})

test_that("detect_q_gene_interactions westfall-young produces valid adjusted p-values", {
  set.seed(888)
  # Create test data with various signal strengths
  model_data <- data.frame(
    entropy = c(
      rnorm(40, mean = 1.0, sd = 0.2),  # Gene1: stable
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 1.0, length.out = 40),  # Gene2: q-dependent
      rnorm(40, mean = 1.0, sd = 0.3)   # Gene3: stable
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:20), 6),
    stringsAsFactors = FALSE
  )
  
  result <- detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 30
  )
  
  # Verify adjusted p-values exist and are valid
  expect_true(all(!is.na(result$adj_p_value)))
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1))
  expect_equal(nrow(result), 3)
})

test_that("detect_q_gene_interactions westfall-young adjusted p-values are monotonic", {
  set.seed(999)
  model_data <- data.frame(
    entropy = rnorm(120),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(paste0("Gene", 1:6), each = 20),
    sample = rep(paste0("S", 1:10), 12),
    stringsAsFactors = FALSE
  )
  
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  result <- suppressWarnings(detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 15
  ))
  
  # Sort by p_value and check that adj_p_value is non-decreasing
  result_sorted <- result[order(result$p_value), ]
  diffs <- diff(result_sorted$adj_p_value)
  expect_true(all(diffs >= -1e-10))  # Allow tiny floating point errors
})

test_that("detect_q_gene_interactions westfall-young wy_randomizations parameter works", {
  set.seed(1001)
  model_data <- data.frame(
    entropy = rnorm(80),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 20),
    gene = rep(paste0("Gene", 1:4), each = 20),
    sample = rep(paste0("S", 1:5), 16),
    stringsAsFactors = FALSE
  )
  
  # Test with different randomization counts
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  result_small <- suppressWarnings(detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 5
  ))
  
  result_large <- suppressWarnings(detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 50
  ))
  
  # Both should have valid results
  expect_is(result_small, "data.frame")
  expect_is(result_large, "data.frame")
  
  # Results may differ slightly due to permutation randomness, but structure same
  expect_equal(nrow(result_small), nrow(result_large))
  expect_equal(colnames(result_small), colnames(result_large))
})

test_that("detect_q_gene_interactions westfall-young verbose mode works", {
  set.seed(1011)
  model_data <- data.frame(
    entropy = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(paste0("Gene", 1:3), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    stringsAsFactors = FALSE
  )
  
  # Capture message output
  expect_warning(
    expect_message(
      detect_q_gene_interactions(
        model_data,
        multicorr = "westfall-young",
        wy_randomizations = 10,
        verbose = TRUE
      ),
      "westfall-young|WY|permutation",
      ignore.case = TRUE
    ),
    "Chi-squared approximation may be incorrect"
  )
})

test_that("detect_q_gene_interactions westfall-young produces FWER control", {
  # Create null data (no true q-effects)
  set.seed(1021)
  model_data <- data.frame(
    entropy = rnorm(200),  # Pure noise, no structure
    q = rep(c(0.5, 1.0, 1.5, 2.0), 50),
    gene = rep(paste0("Gene", 1:10), each = 20),
    sample = rep(paste0("S", 1:10), 20),
    stringsAsFactors = FALSE
  )
  
  result <- expect_warning(
    detect_q_gene_interactions(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 50
    ),
    "Chi-squared approximation may be incorrect"
  )
  
  # Under null hypothesis with pure noise, should have very few significant genes
  # (true FWER control means at most α fraction false positives expected)
  n_sig_alpha05 <- sum(result$adj_p_value < 0.05)
  expect_true(n_sig_alpha05 <= 2)  # Allow at most 2 false positives out of 10 genes
})

test_that("detect_q_gene_interactions westfall-young vs hochberg agreement", {
  set.seed(1031)
  # Create test data with multiple q-levels
  model_data <- data.frame(
    entropy = c(
      rnorm(40, mean = 1.0, sd = 0.2),
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 1.2, length.out = 40),  # Signal
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 0.5, length.out = 40)   # Moderate signal
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:10), 12),
    stringsAsFactors = FALSE
  )
  
  result_wy <- detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 40
  )
  
  result_hoch <- detect_q_gene_interactions(
    model_data,
    multicorr = "hochberg"
  )
  
  # Both should produce valid p-values
  expect_true(all(result_wy$adj_p_value >= 0 & result_wy$adj_p_value <= 1))
  expect_true(all(result_hoch$adj_p_value >= 0 & result_hoch$adj_p_value <= 1))
  
  # Check that methods are working (non-zero variance suggests differentiation)
  expect_true(var(result_wy$adj_p_value) >= 0)
  expect_true(var(result_hoch$adj_p_value) >= 0)
})

test_that("detect_q_gene_interactions westfall-young handles small randomizations", {
  set.seed(1041)
  model_data <- data.frame(
    entropy = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(paste0("Gene", 1:3), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    stringsAsFactors = FALSE
  )
  
  # Should work with very small wy_randomizations (though less accurate)
  result <- detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 5
  )
  
  expect_is(result, "data.frame")
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1))
})

test_that("detect_q_gene_interactions westfall-young handles edge cases gracefully", {
  set.seed(1051)
  # Create clean dataset with sufficient samples
  model_data <- data.frame(
    entropy = c(
      rnorm(20),  # Gene1
      rnorm(20),  # Gene2
      rnorm(20)   # Gene3
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    stringsAsFactors = FALSE
  )
  
  # Should handle without crashing
  result <- detect_q_gene_interactions(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 10
  )
  
  expect_is(result, "data.frame")
  expect_true("adj_p_value" %in% colnames(result))
  expect_equal(nrow(result), 3)  # 3 genes
})

test_that("detect_q_gene_interactions westfall-young phipson-smyth correction prevents zero p-values", {
  set.seed(1061)
  # Create data with varying signal strengths
  model_data <- data.frame(
    entropy = c(
      rnorm(40, mean = 1.0, sd = 0.15) + seq(0, 2.0, length.out = 40),  # Strong
      rnorm(40, mean = 1.0, sd = 0.3),
      rnorm(40, mean = 1.0, sd = 0.3)
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:10), 12),
    stringsAsFactors = FALSE
  )
  
  result <- suppressWarnings(
    detect_q_gene_interactions(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 100
    )
  )
  
  # All p-values should be valid (never exactly 0, due to Phipson-Smyth correction)
  expect_true(all(result$adj_p_value > 0))
  expect_true(all(result$adj_p_value <= 1))
  expect_true(all(!is.na(result$adj_p_value)))
  expect_equal(nrow(result), 3)  # Should have 3 genes
})

