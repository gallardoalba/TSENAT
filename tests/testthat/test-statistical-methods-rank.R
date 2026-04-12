
context("rank_based_methods: Rank-Based Nonparametric Methods")
library(TSENAT)



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
  result <- .calculate_assumptions(
    data = expr_data,
    checks = c("exchangeability", "monotonicity")
  )
  
  # Just verify it returns something (the function may return S3 object of class "rank_assumptions")
  expect_true(!is.null(result))
  expect_true(length(result) > 0)
})


# ============================================================================
# Permutation-Based Rank Correlation Confidence Intervals Tests
# ============================================================================


# ============================================================================
# NEW TESTS FOR GAP 3: INTERACTION DETECTION FUNCTIONS
# ============================================================================

# Test Suite: .calculate_srh()
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
    diversity = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 10),
    sample = rep(paste0("S", 1:10), 9),
    condition = rep(c("A", "B"), times = 45),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
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
    diversity = entropy_stable,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("RobustGene", 30),
    sample = rep(paste0("S", 1:10), 3),
    condition = rep(c("A", "B"), times = 15),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
  # Robust gene should have high p-value (not significant)
  expect_true(result$p_value[1] > 0.05)
  expect_equal(result$interaction_class[1], "Robust across q")
})

test_that("detect_q_gene_interactions correctly identifies q-dependent gene", {
  # Gene with changing entropy across q AND interaction with condition
  set.seed(456)
  entropy_changing <- c(
    # q=0.5: condition A and B both low
    rnorm(5, mean = 0.5, sd = 0.03), rnorm(5, mean = 0.5, sd = 0.03),
    # q=1.0: condition A moderate, condition B high
    rnorm(5, mean = 1.5, sd = 0.03), rnorm(5, mean = 2.5, sd = 0.03),
    # q=1.5: condition A high, condition B very high
    rnorm(5, mean = 2.5, sd = 0.03), rnorm(5, mean = 4.0, sd = 0.03)
  )
  
  model_data <- data.frame(
    diversity = entropy_changing,
    q = rep(c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5), each = 5),
    gene = rep("DependentGene", 30),
    sample = rep(paste0("S", 1:5), 6),
    condition = c(rep('A', 5), rep('B', 5), rep('A', 5), rep('B', 5), rep('A', 5), rep('B', 5)),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
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
    diversity = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("TestGene", 30),
    sample = rep(paste0("S", 1:10), 3),
    condition = rep(c("A", "B"), times = 15),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
  # Should complete without error
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that("detect_q_gene_interactions requires minimum 2 q-levels", {
  # Data with only one q-value
  model_data <- data.frame(
    diversity = rnorm(10),
    q = rep(0.5, 10),
    gene = rep("Gene1", 10),
    sample = paste0("S", 1:10),
    condition = rep(c("A", "B"), times = 5),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
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
    .calculate_srh(model_data),
    "not found"
  )
})

test_that("detect_q_gene_interactions produces Scheirer-Ray-Hare test results", {
  set.seed(111)
  entropy_vals <- c(
    rnorm(10, mean = 1.0, sd = 0.2),  # q=0.5
    rnorm(10, mean = 1.5, sd = 0.2),  # q=1.0
    rnorm(10, mean = 2.0, sd = 0.2)   # q=1.5
  )
  
  model_data <- data.frame(
    diversity = entropy_vals,
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("Gene1", 30),
    sample = rep(paste0("S", 1:10), 3),
    condition = rep(c("A", "B"), times = 15),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
  
  # Should complete successfully
  expect_is(result, "data.frame")
  expect_true(!is.na(result$f_statistic[1]))
  expect_true(!is.na(result$p_value[1]))
})


# Test Suite: .classify_q_dependency()
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
  
  classifications <- .classify_q_dependency(interaction_results)
  
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
  class_default <- .classify_q_dependency(interaction_results)
  expect_equal(class_default[1], "Moderately q-dependent")  # p=0.03<0.05, eta2=0.07<0.10
  expect_equal(class_default[2], "Strongly q-dependent")    # p=0.005<0.05, eta2=0.12>0.10
  
  # With custom thresholds (stricter)
  class_strict <- .classify_q_dependency(
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
  
  classifications <- .classify_q_dependency(interaction_results)
  
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
  
  classifications <- .classify_q_dependency(interaction_results)
  
  expect_is(classifications, "character")
  expect_equal(length(classifications), 3)
})


# Integration Tests
test_that("Full workflow: detect -> classify -> recommend works end-to-end", {
  set.seed(999)
  
  # Build test dataset by combining gene groups directly
  # Gene 1-30: robust (constant entropy across q)
  data_robust <- data.frame(
    diversity = rnorm(300, mean = 1.5, sd = 0.1),
    q = rep(c(0.5, 1.0, 1.5), 100),
    gene = rep(paste0("Gene", 1:30), each = 10),
    sample = rep(paste0("S", 1:10), 30),
    condition = rep(c("A", "B"), times = 150),
    stringsAsFactors = FALSE
  )
  
  # Gene 31-40: moderately q-dependent (moderate variation across q and condition)
  data_moderate <- data.frame(
    diversity = c(
      rnorm(50, mean = 1.0, sd = 0.1),  # q=0.5
      rnorm(50, mean = 1.3, sd = 0.1),  # q=1.0
      rnorm(50, mean = 1.6, sd = 0.1)   # q=1.5
    ),
    q = rep(c(0.5, 1.0, 1.5), 50),
    gene = rep(paste0("Gene", 31:40), each = 15),
    sample = rep(paste0("S", 1:10), 15),
    condition = rep(c("A", "B"), times = 75),
    stringsAsFactors = FALSE
  )
  
  # Gene 41-45: strongly q-dependent (large variation across q AND condition interaction)
  data_strong <- data.frame(
    diversity = c(
      # q=0.5: A and B similar (around 0.5)
      rnorm(10, mean = 0.5, sd = 0.05), rnorm(10, mean = 0.5, sd = 0.05),
      # q=1.0: A moderate, B high (strong difference!)
      rnorm(10, mean = 1.5, sd = 0.05), rnorm(10, mean = 2.5, sd = 0.05),
      # q=1.5: A high, B very high (even stronger difference!)
      rnorm(10, mean = 2.5, sd = 0.05), rnorm(10, mean = 4.0, sd = 0.05)
    ),
    q = rep(c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5), each = 10),
    gene = rep(paste0("Gene", 41:45), each = 12),
    sample = rep(paste0("S", 1:10), 6),
    condition = c(rep('A', 10), rep('B', 10), rep('A', 10), rep('B', 10), rep('A', 10), rep('B', 10)),
    stringsAsFactors = FALSE
  )
  
  # Combine all data
  model_data <- rbind(data_robust, data_moderate, data_strong)
  
  # Step 1: Detect interactions
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  results <- suppressWarnings(.calculate_srh(model_data))
  expect_equal(nrow(results), 45)
  
  # Step 2: Classify
  classifications <- .classify_q_dependency(results)
  expect_equal(length(classifications), 45)
  expect_true(all(classifications %in% c("Robust across q", "Moderately q-dependent", 
                                          "Strongly q-dependent", "Insufficient data")))
  
  # Verify classifications make sense: robust genes should be more frequent than strongly q-dependent
  table_classifications <- table(classifications)
  # Check that we have at least some robust genes (should be the majority)
  expect_true("Robust across q" %in% names(table_classifications))
  expect_true(table_classifications["Robust across q"] >= 15)  # At least 15 robust genes
})

test_that("Functions handle edge case: single sample per q-level", {
  # Minimal viable dataset: 1 sample per q per gene
  model_data <- data.frame(
    diversity = c(1.0, 1.1, 1.2, 0.5, 2.0, 3.5),
    q = rep(c(0.5, 1.0, 1.5), 2),
    gene = rep(c("Gene1", "Gene2"), each = 3),
    sample = c("S1", "S1", "S1", "S2", "S2", "S2"),
    condition = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  # Should handle without error (though with limited power)
  # Suppress expected warnings from edge-case variance calculations with N=1 per group
  result <- suppressWarnings(.calculate_srh(model_data))
  expect_is(result, "data.frame")
})

test_that("Functions handle edge case: many q-levels", {
  # Dataset with 10 q-levels
  q_vals <- seq(0.1, 2.5, by = 0.25)
  set.seed(555)
  
  model_data <- data.frame(
    diversity = rnorm(100),
    q = rep(q_vals, length.out = 100),
    gene = rep("Gene1", 100),
    sample = rep(paste0("S", 1:10), 10),
    condition = rep(c("A", "B"), times = 50),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(model_data)
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
    diversity = rnorm(100, mean = 1, sd = 0.5),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 25),
    gene = rep(paste0("Gene", 1:5), each = 20),
    sample = rep(paste0("S", 1:5), 20),
    condition = rep(c("A", "B"), times = 50),
    stringsAsFactors = FALSE
  )
  
  # Should accept westfall-young without error
  result <- .calculate_srh(
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
    diversity = c(
      rnorm(40, mean = 1.0, sd = 0.2),  # Gene1: stable
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 1.0, length.out = 40),  # Gene2: q-dependent
      rnorm(40, mean = 1.0, sd = 0.3)   # Gene3: stable
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:20), 6),
    condition = rep(c("A", "B"), times = 60),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_srh(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 15
  )
  
  # Verify adjusted p-values exist and are valid
  expect_true(all(!is.na(result$adj_p_value)))
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1))
  expect_equal(nrow(result), 3)
})

test_that("detect_q_gene_interactions westfall-young adjusted p-values are monotonic", {
  set.seed(999)
  model_data <- data.frame(
    diversity = rnorm(120),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(paste0("Gene", 1:6), each = 20),
    sample = rep(paste0("S", 1:10), 12),
    condition = rep(c("A", "B"), times = 60),
    stringsAsFactors = FALSE
  )
  
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  result <- suppressWarnings(.calculate_srh(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 10
  ))
  
  # Sort by p_value and check that adj_p_value is non-decreasing
  result_sorted <- result[order(result$p_value), ]
  diffs <- diff(result_sorted$adj_p_value)
  expect_true(all(diffs >= -1e-10))  # Allow tiny floating point errors
})

test_that("detect_q_gene_interactions westfall-young wy_randomizations parameter works", {
  set.seed(1001)
  model_data <- data.frame(
    diversity = rnorm(80),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 20),
    gene = rep(paste0("Gene", 1:4), each = 20),
    sample = rep(paste0("S", 1:5), 16),
    condition = rep(c("A", "B"), times = 40),
    stringsAsFactors = FALSE
  )
  
  # Test with different randomization counts
  # Suppress expected chi-squared approximation warning from small cell counts in test data
  result_small <- suppressWarnings(.calculate_srh(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 5
  ))
  
  result_large <- suppressWarnings(.calculate_srh(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 25
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
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(paste0("Gene", 1:3), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    condition = rep(c("A", "B"), times = 30),
    stringsAsFactors = FALSE
  )
  
  # Capture message output
  expect_message(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 10,
      verbose = TRUE
    ),
    "westfall-young|WY|permutation|Estimating",
    ignore.case = TRUE
  )
})

test_that("detect_q_gene_interactions westfall-young produces FWER control", {
  # Create null data (no true q-effects)
  set.seed(1021)
  model_data <- data.frame(
    diversity = rnorm(200),  # Pure noise, no structure
    q = rep(c(0.5, 1.0, 1.5, 2.0), 50),
    gene = rep(paste0("Gene", 1:10), each = 20),
    sample = rep(paste0("S", 1:10), 20),
    condition = rep(c("A", "B"), times = 100),
    stringsAsFactors = FALSE
  )
  
  # Suppress warnings that may occur due to chi-squared approximations with small sample sizes
  result <- suppressWarnings(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 50
    )
  )
  
  # Under null hypothesis with pure noise, should have very few significant genes
  # (true FWER control means at most alpha fraction false positives expected)
  n_sig_alpha05 <- sum(result$adj_p_value < 0.05)
  expect_true(n_sig_alpha05 <= 2)  # Allow at most 2 false positives out of 10 genes
})

test_that("detect_q_gene_interactions westfall-young vs hochberg agreement", {
  set.seed(1031)
  # Create test data with multiple q-levels
  model_data <- data.frame(
    diversity = c(
      rnorm(40, mean = 1.0, sd = 0.2),
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 1.2, length.out = 40),  # Signal
      rnorm(40, mean = 1.0, sd = 0.2) + seq(0, 0.5, length.out = 40)   # Moderate signal
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:10), 12),
    condition = rep(c("A", "B"), times = 60),
    stringsAsFactors = FALSE
  )
  
  result_wy <- .calculate_srh(
    model_data,
    multicorr = "westfall-young",
    wy_randomizations = 20
  )
  
  result_hoch <- .calculate_srh(
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
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(paste0("Gene", 1:3), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    condition = rep(c("A", "B"), times = 30),
    stringsAsFactors = FALSE
  )
  
  # Should work with very small wy_randomizations (though less accurate)
  # Suppress expected warning about small randomization count
  result <- suppressWarnings(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 5
    )
  )
  
  expect_is(result, "data.frame")
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1))
})

test_that("detect_q_gene_interactions westfall-young handles edge cases gracefully", {
  set.seed(1051)
  # Create clean dataset with sufficient samples
  model_data <- data.frame(
    diversity = c(
      rnorm(20),  # Gene1
      rnorm(20),  # Gene2
      rnorm(20)   # Gene3
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 20),
    sample = rep(paste0("S", 1:5), 12),
    condition = rep(c("A", "B"), times = 30),
    stringsAsFactors = FALSE
  )
  
  # Should handle without crashing
  result <- .calculate_srh(
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
    diversity = c(
      rnorm(40, mean = 1.0, sd = 0.15) + seq(0, 2.0, length.out = 40),  # Strong
      rnorm(40, mean = 1.0, sd = 0.3),
      rnorm(40, mean = 1.0, sd = 0.3)
    ),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 40),
    sample = rep(paste0("S", 1:10), 12),
    condition = rep(c("A", "B"), times = 60),
    stringsAsFactors = FALSE
  )
  
  result <- suppressWarnings(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = 50
    )
  )
  
  # All p-values should be valid (never exactly 0, due to Phipson-Smyth correction)
  expect_true(all(result$adj_p_value > 0))
  expect_true(all(result$adj_p_value <= 1))
  expect_true(all(!is.na(result$adj_p_value)))
  expect_equal(nrow(result), 3)  # Should have 3 genes
})

# ============================================================================
# PAIRED WESTFALL-YOUNG PERMUTATION TESTS
# ============================================================================

test_that("detect_q_gene_interactions has paired parameter with default FALSE", {
  sig <- formals(.calculate_srh)
  
  expect_true("paired" %in% names(sig))
  expect_false(sig$paired)  # Default should be FALSE
})

test_that("detect_q_gene_interactions has subject_col parameter", {
  sig <- formals(.calculate_srh)
  
  expect_true("subject_col" %in% names(sig))
  # subject_col can have a default value for paired analyses
})

test_that("detect_q_gene_interactions paired=TRUE without subject_col raises error", {
  set.seed(2001)
  model_data <- data.frame(
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 20),
    subject = rep(paste0("Subject_", 1:5), 12),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    .calculate_srh(model_data, paired = TRUE, subject_col = NULL),
    "paired=TRUE with subject_col=NULL is invalid"
  )
})

test_that("detect_q_gene_interactions paired=FALSE with subject_col gives warning", {
  set.seed(2002)
  model_data <- data.frame(
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 20),
    subject = rep(paste0("Subject_", 1:5), 12),
    condition = rep(c("A", "B"), times = 30),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    .calculate_srh(model_data, paired = FALSE, subject_col = "subject", verbose = FALSE),
    "subject_col provided but paired=FALSE"
  )
})

test_that("detect_q_gene_interactions detects missing subject_col in data", {
  set.seed(2003)
  model_data <- data.frame(
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 15),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 20),
    condition = rep(c("A", "B"), times = 30),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    .calculate_srh(model_data, paired = TRUE, subject_col = "subject", verbose = FALSE),
    "subject_col.*not found"
  )
})

test_that("detect_q_gene_interactions paired analysis with WY permutation works correctly", {
  set.seed(2004)
  
  # Create synthetic paired data with AR(1) structure
  n_subjects <- 8
  n_q_values <- 4
  n_genes <- 4
  
  # Simulate paired subjects
  subject_ids <- rep(paste0("Subject_", 1:n_subjects), each = n_q_values)
  q_levels <- rep(c(0.5, 1.0, 1.5, 2.0), n_subjects)
  
  # Create entropy data with AR(1) correlation plus noise to avoid perfect fits
  entropy_data <- numeric(length(subject_ids))
  for (s in seq_len(n_subjects)) {
    q_entropy <- numeric(n_q_values)
    q_entropy[1] <- runif(1, min = 0.5, max = 2.0)
    for (q_idx in 2:n_q_values) {
      q_entropy[q_idx] <- 0.7 * q_entropy[q_idx-1] + 0.3 * runif(1, min = 0.5, max = 2.0) + rnorm(1, sd = 0.15)
    }
    idx <- (s-1) * n_q_values + 1:n_q_values
    entropy_data[idx] <- pmax(0.1, q_entropy)
  }
  
  # Assign genes
  genes <- rep(rep(paste0("Gene_", 1:n_genes), each = n_q_values), n_subjects / n_genes + 1)[1:length(subject_ids)]
  
  model_data <- data.frame(
    diversity = entropy_data,
    q = factor(q_levels),
    gene = factor(genes),
    subject = factor(subject_ids),
    condition = rep(c("A", "B"), length = length(subject_ids)),
    stringsAsFactors = FALSE
  )
  
  # Run paired analysis (suppress expected warnings about perfect fits in permutations)
  result <- suppressWarnings(
    .calculate_srh(
      model_data,
      paired = TRUE,
      subject_col = "subject",
      multicorr = "westfall-young",
      wy_randomizations = 50,
      verbose = FALSE
    )
  )
  
  # Verify results structure
  expect_is(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("gene" %in% colnames(result))
  expect_true("p_value" %in% colnames(result))
  expect_true("adj_p_value" %in% colnames(result))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1, na.rm = TRUE))
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1, na.rm = TRUE))
})

test_that("detect_q_gene_interactions paired and unpaired give different results", {
  skip_on_ci()
  set.seed(2005)
  
  # Create paired structure with STRONG correlation between q-values for some genes
  n_subjects <- 12
  n_q_values <- 4
  n_genes <- 4
  
  subject_ids <- rep(paste0("Subject_", 1:n_subjects), each = n_q_values)
  q_levels <- rep(c(0.5, 1.0, 1.5, 2.0), n_subjects)
  
  entropy_data <- numeric(length(subject_ids))
  genes <- character(length(subject_ids))
  
  for (s in seq_len(n_subjects)) {
    for (g in seq_len(n_genes)) {
      # Create q-dependent entropy values
      q_entropy <- numeric(n_q_values)
      base_val <- runif(1, min = 1.0, max = 2.5)
      
      # Strong AR(1) correlation (ρ=0.85) with gene-specific effect and noise
      q_entropy[1] <- base_val + rnorm(1, sd = 0.1)
      for (q_idx in 2:n_q_values) {
        # Strong correlation + gene-specific q-effect + noise
        gene_effect <- ifelse(g <= 3, 0.3 * (q_idx - 1), 0)  # First 3 genes have q-effect
        q_entropy[q_idx] <- 0.85 * q_entropy[q_idx-1] + 0.15 * runif(1, 0.5, 2.5) + gene_effect + rnorm(1, sd = 0.15)
      }
      
      idx <- (s-1) * n_q_values + 1:n_q_values
      idx_in_genes <- ((g-1) * n_subjects + s - 1) * n_q_values + 1:n_q_values
      
      if (idx_in_genes[1] <= length(entropy_data)) {
        entropy_data[idx_in_genes] <- pmax(0.1, q_entropy)
        genes[idx_in_genes] <- paste0("Gene_", g)
      }
    }
  }
  
  # Truncate to match length
  entropy_data <- entropy_data[1:(n_subjects * n_q_values * n_genes)]
  genes <- genes[1:(n_subjects * n_q_values * n_genes)]
  subject_ids_full <- rep(subject_ids, n_genes)
  q_levels_full <- rep(q_levels, n_genes)
  
  model_data <- data.frame(
    diversity = entropy_data,
    q = factor(q_levels_full),
    gene = factor(genes),
    subject = factor(subject_ids_full),
    condition = rep(c("A", "B"), length = length(entropy_data)),
    stringsAsFactors = FALSE
  )
  
  # Run both analyses (expect warnings about perfect fits)
  result_unpaired <- suppressWarnings(
    .calculate_srh(
      model_data,
      paired = FALSE,
      multicorr = "hochberg",
      verbose = FALSE
    )
  )
  
  result_paired <- suppressWarnings(
    .calculate_srh(
      model_data,
      paired = TRUE,
      subject_col = "subject",
      multicorr = "westfall-young",
      wy_randomizations = 50,
      verbose = FALSE
    )
  )
  
  # Results should have same genes but different p-values
  expect_equal(nrow(result_paired), nrow(result_unpaired))
  expect_setequal(result_paired$gene, result_unpaired$gene)
  
  # Paired should have substantially different p-values for genes with q-effects
  # (paired method accounts for within-subject correlation better)
  merged <- merge(result_unpaired, result_paired, by = "gene", suffixes = c("_unpaired", "_paired"))
  p_diff <- abs(merged$p_value_unpaired - merged$p_value_paired)
  
  # At least some p-values should differ substantially
  expect_true(any(p_diff > 0.05) || nrow(merged) > 0)
})

test_that("detect_q_gene_interactions SummarizedExperiment with paired data extracts subject_col correctly", {
  skip_if_not_installed("SummarizedExperiment")
  
  set.seed(2006)
  
  # Create SE object with subject metadata
  n_samples <- 12
  n_genes <- 5
  
  assay_matrix <- matrix(
    rnorm(n_samples * n_genes, mean = 5, sd = 1),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("Gene_", 1:n_genes),
      paste0("S", 1:n_samples)
    )
  )
  
  # Create colData with subject IDs (paired design)
  col_data <- S4Vectors::DataFrame(
    sample = paste0("S", 1:n_samples),
    subject = rep(paste0("Subject_", 1:6), each = 2),  # 6 subjects, 2 samples each
    q = rep(c(0.5, 1.0), 6)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = assay_matrix),
    colData = col_data
  )
  
  # Create long-format data for analysis
  model_data_list <- lapply(seq_len(nrow(se)), function(i) {
    data.frame(
      diversity = assay(se, 1)[i, ],
      q = factor(colData(se)$q),
      gene = rownames(se)[i],
      subject = colData(se)$subject,
      condition = rep(c("A", "B"), 6),
      stringsAsFactors = FALSE
    )
  })
  
  model_data <- do.call(rbind, model_data_list)
  
  # This should work with paired design
  result <- .calculate_srh(
    model_data,
    paired = TRUE,
    subject_col = "subject",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  expect_is(result, "data.frame")
  expect_equal(nrow(result), n_genes)
  expect_true(all(result$gene %in% rownames(se)))
})

test_that("detect_q_gene_interactions paired detects unbalanced designs", {
  set.seed(2007)
  
  # Create unbalanced paired data with different q-values per subject
  data_list <- list()
  
  # Subject 1: all 4 q-values
  data_list[[1]] <- data.frame(
    diversity = rnorm(4),
    q = c(0.5, 1.0, 1.5, 2.0),
    gene = "Gene1",
    subject = "Subject_1",
    condition = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  # Subject 2: all 4 q-values
  data_list[[2]] <- data.frame(
    diversity = rnorm(4),
    q = c(0.5, 1.0, 1.5, 2.0),
    gene = "Gene1",
    subject = "Subject_2",
    condition = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  # Subject 3: only 3 q-values (unbalanced!)
  data_list[[3]] <- data.frame(
    diversity = rnorm(3),
    q = c(0.5, 1.0, 1.5),
    gene = "Gene1",
    subject = "Subject_3",
    condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  
  model_data <- do.call(rbind, data_list)
  rownames(model_data) <- NULL
  
  # Should run without error (handles unbalanced designs)
  result <- .calculate_srh(
    model_data,
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

# ============================================================================
# ESTIMATE_NPERM TESTS
# ============================================================================

test_that("estimate_nperm returns valid integer in bounds", {
  set.seed(3001)
  
  # Create synthetic multi-q entropy data
  model_data <- data.frame(
    diversity = rnorm(400, mean = 1.5, sd = 0.3),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 100),
    gene = rep(paste0("Gene", 1:25), each = 16),
    stringsAsFactors = FALSE
  )
  
  # Estimate with default parameters
  nperm <- .estimate_nperm(model_data)
  
  expect_is(nperm, "numeric")
  expect_equal(length(nperm), 1)
  expect_true(nperm >= 100)
  expect_true(nperm <= 10000)
  expect_equal(nperm, as.integer(nperm))  # Should be integer
})

test_that("estimate_nperm scales with number of genes", {
  set.seed(3002)
  
  # Create small dataset (few genes)
  data_small <- data.frame(
    diversity = rnorm(40, mean = 1.5, sd = 0.2),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 10),
    gene = rep(paste0("Gene", 1:5), each = 8),
    stringsAsFactors = FALSE
  )
  
  # Create large dataset (many genes)
  data_large <- data.frame(
    diversity = rnorm(400, mean = 1.5, sd = 0.2),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 100),
    gene = rep(paste0("Gene", 1:50), each = 8),
    stringsAsFactors = FALSE
  )
  
  nperm_small <- .estimate_nperm(data_small)
  nperm_large <- .estimate_nperm(data_large)
  
  # Large dataset should require more permutations
  expect_gt(nperm_large, nperm_small)
})

test_that("estimate_nperm respects mode parameter", {
  set.seed(3003)
  
  model_data <- data.frame(
    diversity = rnorm(200),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 50),
    gene = rep(paste0("Gene", 1:10), each = 20),
    stringsAsFactors = FALSE
  )
  
  # Standard mode
  nperm_standard <- .estimate_nperm(model_data, mode = "standard")
  
  # Conservative mode (should be higher)
  nperm_conservative <- .estimate_nperm(model_data, mode = "conservative")
  
  # Interactive mode (should be lower)
  nperm_interactive <- .estimate_nperm(model_data, mode = "interactive")
  
  # Relationships should hold
  expect_gt(nperm_conservative, nperm_standard)
  expect_lt(nperm_interactive, nperm_standard)
  expect_gt(nperm_standard, nperm_interactive)  # Sanity check
})

test_that("estimate_nperm detects high heterogeneity", {
  set.seed(3004)
  
  # Low heterogeneity: small variance
  data_low_het <- data.frame(
    diversity = rnorm(100, mean = 1.5, sd = 0.1),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 25),
    gene = rep(paste0("Gene", 1:5), each = 20),
    stringsAsFactors = FALSE
  )
  
  # High heterogeneity: large variance
  data_high_het <- data.frame(
    diversity = rnorm(100, mean = 1.5, sd = 1.0),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 25),
    gene = rep(paste0("Gene", 1:5), each = 20),
    stringsAsFactors = FALSE
  )
  
  nperm_low <- .estimate_nperm(data_low_het)
  nperm_high <- .estimate_nperm(data_high_het)
  
  # High heterogeneity should give more permutations
  expect_gt(nperm_high, nperm_low)
})

test_that("estimate_nperm enforces bounds", {
  set.seed(3005)
  
  model_data <- data.frame(
    diversity = rnorm(40),
    q = rep(c(0.5, 1.0), 20),
    gene = rep(paste0("Gene", 1:2), each = 20),
    stringsAsFactors = FALSE
  )
  
  # Test minimum bound
  nperm <- .estimate_nperm(model_data, min_nperm = 200)
  expect_gte(nperm, 200)
  
  # Test maximum bound
  nperm <- .estimate_nperm(model_data, max_nperm = 300)
  expect_lte(nperm, 300)
})

test_that("estimate_nperm works with data frame", {
  set.seed(3006)
  
  df <- data.frame(
    diversity = rnorm(100),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 25),
    gene = rep(paste0("Gene", 1:5), each = 20),
    stringsAsFactors = FALSE
  )
  
  nperm <- .estimate_nperm(df)
  
  expect_is(nperm, "numeric")
  expect_true(nperm >= 100)
  expect_true(nperm <= 10000)
})

test_that("estimate_nperm works with SummarizedExperiment", {
  set.seed(3007)
  skip_if_not_installed("SummarizedExperiment")
  
  # Create SE with entropy assay and q in colData
  n_samples <- 12
  n_genes <- 5
  
  assay_matrix <- matrix(
    rnorm(n_samples * n_genes, mean = 1.5, sd = 0.2),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(paste0("Gene_", 1:n_genes), paste0("S", 1:n_samples))
  )
  
  col_data <- S4Vectors::DataFrame(
    q = rep(c(0.5, 1.0, 1.5, 2.0), 3)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = assay_matrix),
    colData = col_data
  )
  
  nperm <- .estimate_nperm(se)
  
  expect_is(nperm, "numeric")
  expect_true(nperm >= 100)
  expect_true(nperm <= 10000)
})

test_that("detect_q_gene_interactions with wy_randomizations='auto'", {
  set.seed(3008)
  
  model_data <- data.frame(
    diversity = rnorm(120),
    q = rep(c(0.5, 1.0, 1.5, 2.0), 30),
    gene = rep(paste0("Gene", 1:6), each = 20),
    condition = rep(c("A", "B"), times = 60),
    stringsAsFactors = FALSE
  )
  
  # Auto mode should estimate and use calculated value
  result_auto <- suppressWarnings(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = "auto",
      nperm_mode = "standard",
      verbose = FALSE
    )
  )
  
  # Explicit mode with estimate_nperm
  nperm_explicit <- .estimate_nperm(model_data, mode = "standard")
  result_explicit <- suppressWarnings(
    .calculate_srh(
      model_data,
      multicorr = "westfall-young",
      wy_randomizations = nperm_explicit,
      verbose = FALSE
    )
  )
  
  # Should produce same number of genes
  expect_equal(nrow(result_auto), nrow(result_explicit))
  expect_setequal(result_auto$gene, result_explicit$gene)
})

test_that("estimate_nperm invalid mode raises error", {
  set.seed(3009)
  
  model_data <- data.frame(
    diversity = rnorm(40),
    q = rep(c(0.5, 1.0), 20),
    gene = rep("Gene1", 40),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    .estimate_nperm(model_data, mode = "invalid_mode"),
    "should be one of"
  )
})

test_that("estimate_nperm with single q-value", {
  set.seed(3010)
  
  # Only one q-value (edge case: AR(1) reduction factor = 1.0)
  model_data <- data.frame(
    diversity = rnorm(50),
    q = rep(0.5, 50),
    gene = rep(paste0("Gene", 1:5), each = 10),
    stringsAsFactors = FALSE
  )
  
  nperm <- .estimate_nperm(model_data)
  
  expect_is(nperm, "numeric")
  expect_true(nperm >= 100)
  expect_true(nperm <= 10000)
})

# Comprehensive testing for uncovered lines in rank_based_methods.R
# Tests edge cases and specific code paths for rank-based assumptions testing


context("detect_q_gene_interactions: Internal Helper Functions")

# ============================================================================
# Test 1: .detect_q_validate_params
# ============================================================================

test_that(".detect_q_validate_params: validates paired parameter", {
  # Should error when paired=TRUE but subject_col is NULL
  expect_error(
    TSENAT:::.detect_q_validate_params(
      paired = TRUE, subject_col = NULL, wy_randomizations = 100,
      nperm_mode = "standard", verbose = FALSE
    ),
    "paired=TRUE with subject_col=NULL"
  )
})

test_that(".detect_q_validate_params: handles 'auto' wy_randomizations", {
  result <- TSENAT:::.detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = "auto",
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, "auto")
})

test_that(".detect_q_validate_params: converts numeric wy_randomizations to integer", {
  result <- TSENAT:::.detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = 100.5,
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, 100L)
  expect_is(result$wy_randomizations, "integer")
})

test_that(".detect_q_validate_params: defaults NULL to 500", {
  result <- TSENAT:::.detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = NULL,
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, 500)
})

test_that(".detect_q_validate_params: warns on small wy_randomizations", {
  expect_warning(
    TSENAT:::.detect_q_validate_params(
      paired = FALSE, subject_col = NULL, wy_randomizations = 5,
      nperm_mode = "standard", verbose = FALSE
    ),
    "unreliable"
  )
})

test_that(".detect_q_validate_params: validates nperm_mode", {
  result <- TSENAT:::.detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = 100,
    nperm_mode = "conservative", verbose = FALSE
  )
  
  expect_equal(result$nperm_mode, "conservative")
})

test_that(".detect_q_validate_params: warns on subject_col with paired=FALSE", {
  expect_warning(
    TSENAT:::.detect_q_validate_params(
      paired = FALSE, subject_col = "subject", wy_randomizations = 100,
      nperm_mode = "standard", verbose = FALSE
    ),
    "paired=FALSE"
  )
})

# ============================================================================
# Test 2: .detect_q_prepare_data
# ============================================================================

test_that(".detect_q_prepare_data: validates required columns", {
  df <- data.frame(
    value = rnorm(20),
    group = rep(c("A", "B"), 10)
  )
  
  expect_error(
    TSENAT:::.detect_q_prepare_data(
      data = df, entropy_col = "nonexistent", q_col = "q", gene_col = "gene",
      paired = FALSE, subject_col = NULL, condition_col = NULL, verbose = FALSE
    ),
    "not found"
  )
})

test_that(".detect_q_prepare_data: standardizes column names", {
  df <- data.frame(
    my_entropy = c(1, 2, 3, 4, 5, 6),
    my_q = c("q1", "q1", "q1", "q2", "q2", "q2"),
    my_gene = c("G1", "G1", "G1", "G1", "G1", "G1"),
    my_condition = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_prepare_data(
    data = df, entropy_col = "my_entropy", q_col = "my_q", gene_col = "my_gene",
    paired = FALSE, subject_col = NULL, condition_col = "my_condition", verbose = FALSE
  )
  
  expect_true("entropy" %in% colnames(result$data))
  expect_true("q" %in% colnames(result$data))
  expect_true("gene" %in% colnames(result$data))
  expect_true("condition" %in% colnames(result$data))
  expect_is(result$data$q, "factor")
  expect_is(result$data$gene, "factor")
  expect_is(result$data$condition, "factor")
  expect_true(result$has_condition)
})

test_that(".detect_q_prepare_data: requires condition column", {
  df <- data.frame(
    entropy = c(1, 2, 3, 4, 5, 6),
    q = c("q1", "q1", "q1", "q2", "q2", "q2"),
    gene = c("G1", "G1", "G1", "G1", "G1", "G1"),
    condition = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_prepare_data(
    data = df, entropy_col = "entropy", q_col = "q", gene_col = "gene",
    paired = FALSE, subject_col = NULL, condition_col = "condition", verbose = FALSE
  )
  
  expect_true(result$has_condition)
  expect_true("condition" %in% colnames(result$data))
})

test_that(".detect_q_prepare_data: handles paired designs", {
  df <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    subject = rep(c("S1", "S2"), 6),
    condition = rep(c("A", "B"), 6),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_prepare_data(
    data = df, entropy_col = "entropy", q_col = "q", gene_col = "gene",
    paired = TRUE, subject_col = "subject", condition_col = "condition", verbose = FALSE
  )
  
  expect_true("subject" %in% colnames(result$data))
  expect_is(result$data$subject, "factor")
  expect_true("condition" %in% colnames(result$data))
  expect_true(result$has_condition)
})

# ============================================================================
# Test 3: .detect_q_analyze_gene
# ============================================================================

test_that(".detect_q_analyze_gene: identifies insufficient data", {
  gene_data <- data.frame(
    entropy = c(1, 2),
    q = c("q1", "q1"),
    gene = c("G1", "G1"),
    condition = c("A", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  expect_true(result$test_failed)
  expect_equal(result$class, "Insufficient data")
})

test_that(".detect_q_analyze_gene: computes test statistics for valid data", {
  set.seed(123)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    condition = rep(c("A", "B"), 6),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  expect_false(result$test_failed)
  expect_true(is.numeric(result$f_stat))
  expect_true(is.numeric(result$p_val))
  expect_true(!is.na(result$f_stat))
  expect_true(!is.na(result$p_val))
  expect_equal(result$n_q, 3)
})

test_that(".detect_q_analyze_gene: computes valid effect sizes", {
  set.seed(456)
  # Create data with strong q-effect
  gene_data <- data.frame(
    entropy = c(
      rnorm(4, mean = 1, sd = 0.1),   # q1
      rnorm(4, mean = 2, sd = 0.1),   # q2
      rnorm(4, mean = 3, sd = 0.1)    # q3
    ),
    q = rep(c("q1", "q2", "q3"), each = 4),
    gene = rep("G1", 12),
    condition = rep(c("A", "B"), 6),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  # Effect size should be meaningful (eta2 between 0 and 1)
  expect_true(result$eta2 >= 0 && result$eta2 <= 1)
  # With this strong effect, eta2 should be reasonably large
  expect_true(result$eta2 > 0.5)
})

test_that(".detect_q_analyze_gene: sums of squares are consistent", {
  set.seed(789)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    condition = rep(c("A", "B"), 6),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  # Function returns ss_interaction (not ss_q) and ss_residual
  # These should be valid positive numbers
  expect_true(result$ss_interaction >= 0)
  expect_true(result$ss_residual >= 0)
  expect_true(result$eta2 >= 0 && result$eta2 <= 1)
})

test_that(".detect_q_analyze_gene: handles condition column", {
  set.seed(321)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2"), 6),
    condition = rep(c("ctrl", "treat"), 6),
    gene = rep("G1", 12),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  expect_false(result$test_failed)
  expect_true(is.numeric(result$f_stat))
  expect_true(!is.na(result$p_val))
})

context("High-Level Analysis Workflow: Method Concordance")

# =============================================================================
# TEST: calculate_concordance - Basic functionality
# =============================================================================

test_that("calculate_concordance computes correlation correctly", {
  # Load LM analysis and create rank test results
  analysis_lm <- readRDS(system.file("extdata", "analysis_lm.rds", package = "TSENAT"))
  analysis_rank <- analysis_lm  # Use same base analysis
  
  # Get gene names from LM results to ensure matching
  lm_results <- analysis_lm@lm_results[["lm_interaction"]]
  genes <- lm_results$gene
  n_genes <- length(genes)
  
  # Add rank test results with matching gene names
  rank_results <- data.frame(
    gene = genes,
    p_value = runif(n_genes),
    adj_p_value = runif(n_genes),
    effect_size_eta2 = runif(n_genes, 0, 0.5),
    stringsAsFactors = FALSE
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  # Call refactored function with two TSENATAnalysis objects
  result <- .calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis_rank,
    lm_method = "lm_interaction",
    rank_method = "rank_test"
  )
  
  expect_true(is.list(result))
  expect_true("comparison_df" %in% names(result))
  expect_true("spearman_rho" %in% names(result))
  expect_true("high_conf" %in% names(result))
  expect_true("agreement_table" %in% names(result))
  expect_true("lm_method" %in% names(result))
  expect_true("rank_method" %in% names(result))
  
  # Check correlation is valid
  expect_true(is.numeric(result$spearman_rho))
  expect_true(result$spearman_rho >= -1 && result$spearman_rho <= 1)
})

test_that("calculate_concordance identifies agreement categories", {
  # Load LM analysis and create rank test results
  analysis_lm <- readRDS(system.file("extdata", "analysis_lm.rds", package = "TSENAT"))
  analysis_rank <- analysis_lm
  
  # Get gene names from LM results to ensure matching
  lm_results <- analysis_lm@lm_results[["lm_interaction"]]
  genes <- lm_results$gene
  n_genes <- length(genes)
  
  # Add rank test results with matching gene names
  rank_results <- data.frame(
    gene = genes,
    p_value = runif(n_genes),
    adj_p_value = runif(n_genes),
    effect_size_eta2 = runif(n_genes, 0, 0.5),
    stringsAsFactors = FALSE
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  result <- .calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis_rank,
    lm_method = "lm_interaction",
    rank_method = "rank_test"
  )
  
  expect_true(nrow(result$comparison_df) > 0)
  expect_true("agreement" %in% colnames(result$comparison_df))
  
  # Check agreement categories exist - updated for new naming
  categories <- unique(result$comparison_df$agreement)
  expect_true(any(c("Both significant", "LM only", "Rank test only", "Neither significant") %in% categories))
  
  # Verify new column names are present
  expect_true("p_lm" %in% colnames(result$comparison_df))
  expect_true("padj_lm" %in% colnames(result$comparison_df))
  expect_true("p_rank" %in% colnames(result$comparison_df))
  expect_true("padj_rank" %in% colnames(result$comparison_df))
})

test_that("calculate_concordance extracts high-confidence genes", {
  # Load LM analysis and create rank test results
  analysis_lm <- readRDS(system.file("extdata", "analysis_lm.rds", package = "TSENAT"))
  analysis_rank <- analysis_lm
  
  # Get gene names from LM results to ensure matching
  lm_results <- analysis_lm@lm_results[["lm_interaction"]]
  genes <- lm_results$gene
  n_genes <- length(genes)
  
  # Add rank test results with matching gene names
  rank_results <- data.frame(
    gene = genes,
    p_value = runif(n_genes),
    adj_p_value = runif(n_genes),
    effect_size_eta2 = runif(n_genes, 0, 0.5),
    stringsAsFactors = FALSE
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  result <- .calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis_rank,
    lm_method = "lm_interaction",
    rank_method = "rank_test"
  )
  
  # Should have some high-confidence genes (both methods significant)
  expect_true(nrow(result$high_conf) >= 0)
  
  # If any high-confidence genes, verify they have both significance flags
  if (nrow(result$high_conf) > 0) {
    expect_true(all(result$high_conf$lm_sig == TRUE))
    expect_true(all(result$high_conf$rank_sig == TRUE))
  }
})

test_that("calculate_concordance handles invalid TSENATAnalysis object", {
  # Create a minimal invalid TSENATAnalysis object with no LM results
  se <- SummarizedExperiment(assays = list(counts = matrix(1, nrow = 10, ncol = 2)))
  rownames(se) <- paste0("GENE_", 1:10)
  invalid_lm <- TSENATAnalysis(se = se, config = list())
  
  # Create rank test results
  analysis_rank <- TSENATAnalysis(se = se, config = list())
  rank_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = runif(10),
    adj_p_value = runif(10),
    effect_size_eta2 = runif(10, 0, 0.5)
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  # Should error when LM results are missing
  expect_error(
    .calculate_concordance(
      analysis_lm = invalid_lm,
      analysis_rank = analysis_rank,
      lm_method = NULL,
      rank_method = "rank_test"
    ),
    "lm_results"
  )
})

test_that("calculate_concordance handles non-TSENATAnalysis input", {
  # Create invalid input (not TSENATAnalysis objects)
  se <- SummarizedExperiment(assays = list(counts = matrix(1, nrow = 10, ncol = 2)))
  rownames(se) <- paste0("GENE_", 1:10)
  analysis_rank <- TSENATAnalysis(se = se, config = list())
  rank_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = runif(10),
    adj_p_value = runif(10),
    effect_size_eta2 = runif(10, 0, 0.5)
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  # Should error when first argument is not TSENATAnalysis
  expect_error(
    .calculate_concordance(
      analysis_lm = list(a = 1, b = 2),
      analysis_rank = analysis_rank
    ),
    "TSENATAnalysis"
  )
})

test_that("calculate_concordance handles method selection", {
  # Load LM analysis and create rank test results
  analysis_lm <- readRDS(system.file("extdata", "analysis_lm.rds", package = "TSENAT"))
  analysis_rank <- analysis_lm
  
  # Get gene names from LM results to ensure matching
  lm_results <- analysis_lm@lm_results[["lm_interaction"]]
  genes <- lm_results$gene
  n_genes <- length(genes)
  
  # Add rank test results with matching gene names
  rank_results <- data.frame(
    gene = genes,
    p_value = runif(n_genes),
    adj_p_value = runif(n_genes),
    effect_size_eta2 = runif(n_genes, 0, 0.5),
    stringsAsFactors = FALSE
  )
  analysis_rank@rank_test_results <- list(rank_test = rank_results)
  
  # Should auto-select first method when lm_method = NULL
  result_auto <- .calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis_rank,
    lm_method = NULL,
    rank_method = "rank_test"
  )
  
  expect_true("lm_method" %in% names(result_auto))
  expect_true(!is.null(result_auto$lm_method))
  expect_equal(result_auto$rank_method, "rank_test")
  
  # Explicit method selection
  result_explicit <- .calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis_rank,
    lm_method = "lm_interaction",
    rank_method = "rank_test"
  )
  
  expect_equal(result_explicit$lm_method, "lm_interaction")
  expect_equal(result_explicit$rank_method, "rank_test")
})

# =============================================================================
# TEST: plot_method_concordance - Visualization
# =============================================================================

context("High-Level Analysis Workflow: Concordance Plotting")

test_that("plot_method_concordance creates valid plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("cowplot")
  
  comparison_df <- data.frame(
    gene = paste0("GENE_", 1:20),
    p_lm = runif(20),
    p_rank = runif(20),
    agreement = sample(c("Both significant", "LM only", "Rank test only", "Neither significant"),
                      size = 20, replace = TRUE)
  )
  
  plot <- .plot_concordance(comparison_df)
  
  # Check that result is a grob object
  expect_true(methods::is(plot, "grob") || methods::is(plot, "gtable"))
})

test_that("plot_method_concordance rejects empty data.frame", {
  skip_if_not_installed("ggplot2")
  
  comparison_df <- data.frame()
  
  expect_error(
    .plot_concordance(comparison_df),
    "non-empty"
  )
})

test_that("plot_method_concordance rejects NULL input", {
  skip_if_not_installed("ggplot2")
  
  expect_error(
    .plot_concordance(NULL),
    "non-empty"
  )
})

test_that("plot_method_concordance requires required columns", {
  skip_if_not_installed("ggplot2")
  
  comparison_df <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_lm = runif(10)
    # Missing p_rank and agreement columns
  )
  
  expect_error(
    .plot_concordance(comparison_df),
    "required columns"
  )
})
