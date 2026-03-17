context("Friedman Test for Paired Rank-Based Analysis")

# Setup: Load necessary functions
source("../../R/rank_based_methods.R", local = TRUE)

# ============================================================================
# Test 1: Basic Friedman Test Function
# ============================================================================

test_that(".tsenat_apply_friedman_test returns correct structure", {
  # Create simple paired test data (3 subjects, 3 treatments)
  data <- data.frame(
    value = c(1.2, 1.5, 1.8, 2.1, 2.3, 2.5, 3.0, 3.2, 3.5),
    subject = factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3)),
    treatment = factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "value",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  # Check structure
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  
  # Check values
  expect_is(result$statistic, "numeric")
  expect_is(result$p_value, "numeric")
  expect_match(result$method, "Friedman")
  
  # Check ranges
  expect_gte(result$p_value, 0)
  expect_lte(result$p_value, 1)
  expect_gte(result$statistic, 0)
})

test_that(".tsenat_apply_friedman_test accepts data frame input", {
  # Real entropy-like data
  set.seed(123)
  n_subjects <- 5
  n_q <- 4
  
  entropy_vals <- rnorm(n_subjects * n_q, mean = 2.5, sd = 0.5)
  data <- data.frame(
    entropy = entropy_vals,
    subject = factor(rep(1:n_subjects, n_q)),
    q = factor(rep(1:n_q, each = n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_is(result, "list")
  expect_is(result$statistic, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test 2: Friedman Test Comparison with Known Values
# ============================================================================

test_that(".tsenat_apply_friedman_test matches base R friedman.test", {
  # Create controlled test data
  set.seed(42)
  n_subjects <- 4
  n_treatments <- 3
  
  # Create data in long format
  entropy <- c(2.1, 2.3, 1.9, 2.2, 2.8, 3.1, 2.9, 3.0, 3.5, 3.7, 3.6, 3.8)
  subject <- rep(1:n_subjects, n_treatments)
  treatment <- rep(1:n_treatments, each = n_subjects)
  
  long_data <- data.frame(
    entropy = entropy,
    subject = factor(subject),
    treatment = factor(treatment)
  )
  
  # Call our function
  our_result <- .tsenat_apply_friedman_test(
    data = long_data,
    value_col = "entropy",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  # Call base R function on matrix form
  wide_matrix <- xtabs(entropy ~ subject + treatment)
  base_result <- friedman.test(wide_matrix)
  
  # Compare results
  expect_equal(our_result$statistic, as.numeric(base_result$statistic), tolerance = 1e-10)
  expect_equal(our_result$p_value, as.numeric(base_result$p.value), tolerance = 1e-10)
})

# ============================================================================
# Test 3: Error Handling
# ============================================================================

test_that(".tsenat_apply_friedman_test handles missing data gracefully", {
  # Data with only 1 subject
  data <- data.frame(
    value = c(1.0, 2.0, 3.0),
    subject = factor(c(1, 1, 1)),
    treatment = factor(c("A", "B", "C"))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "value",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  # Should fail gracefully with method = "test_failed"
  expect_equal(result$method, "test_failed")
  expect_true(is.na(result$statistic))
  expect_true(is.na(result$p_value))
})

test_that(".tsenat_apply_friedman_test handles single treatment", {
  # Data with only 1 treatment
  data <- data.frame(
    value = c(1.0, 2.0, 3.0),
    subject = factor(c(1, 2, 3)),
    treatment = factor(c("A", "A", "A"))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "value",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  # Should fail gracefully
  expect_equal(result$method, "test_failed")
})

# ============================================================================
# Test 4: Conditional Rank Test Dispatcher with Pairing
# ============================================================================

test_that(".tsenat_apply_conditional_rank_test selects Friedman for paired", {
  # Create balanced paired data
  set.seed(99)
  data <- data.frame(
    entropy = rnorm(30, mean = 2, sd = 0.3),
    q = factor(rep(1:5, 6)),
    subject = factor(rep(1:6, each = 5))
  )
  
  result <- .tsenat_apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  
  # Should select Friedman
  expect_equal(result$test_type, "friedman")
  expect_match(result$method, "Friedman")
})

test_that(".tsenat_apply_conditional_rank_test uses Kruskal-Wallis when unpaired", {
  # Create data without subject column (unpaired)
  set.seed(99)
  data <- data.frame(
    entropy = rnorm(150, mean = 2, sd = 0.3),
    q = factor(rep(1:5, 30))
  )
  
  result <- .tsenat_apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = FALSE,
    subject_col = NULL,
    verbose = FALSE
  )
  
  # Should NOT be Friedman
  expect_false(result$test_type == "friedman")
  # Could be Kruskal-Wallis or conditional selection
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that(".tsenat_apply_conditional_rank_test requires subject_col for paired", {
  # Create data without explicit subject column for paired analysis
  data <- data.frame(
    entropy = rnorm(20, mean = 2, sd = 0.3),
    q = factor(rep(1:4, 5))
  )
  
  result <- .tsenat_apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = NULL,  # Missing!
    verbose = FALSE
  )
  
  # Should fall back to unpaired selection (not Friedman)
  expect_false(result$test_type == "friedman")
})

# ============================================================================
# Test 5: Friedman Returns Proper Metadata
# ============================================================================

test_that("Friedman result includes characteristics metadata", {
  set.seed(42)
  data <- data.frame(
    entropy = rnorm(24, mean = 2, sd = 0.4),
    q = factor(rep(1:4, 6)),
    subject = factor(rep(1:6, each = 4))
  )
  
  result <- .tsenat_apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  
  # Check characteristics
  expect_true("characteristics" %in% names(result))
  expect_is(result$characteristics, "list")
  
  # For Friedman, should mark pairing_used = TRUE
  if (result$test_type == "friedman") {
    expect_true(result$characteristics$pairing_used)
  }
})

# ============================================================================
# Test 6: Integration with detect_q_gene_interactions
# ============================================================================

test_that("detect_q_gene_interactions uses Friedman for paired=TRUE", {
  skip_if_not_installed("Matrix")
  
  # Create minimal paired data
  set.seed(42)
  n_subjects <- 5
  n_q <- 4
  n_genes <- 2
  
  data_list <- lapply(1:n_genes, function(g) {
    data.frame(
      diversity = rnorm(n_subjects * n_q, mean = 2 + g * 0.3, sd = 0.2),
      q = factor(rep(seq(0.5, 2.0, length.out = n_q), each = n_subjects)),
      paired_samples = factor(rep(1:n_subjects, n_q)),
      gene = paste0("gene_", g)
    )
  })
  
  test_data <- do.call(rbind, data_list)
  rownames(test_data) <- NULL
  
  # Run paired analysis
  results_paired <- detect_q_gene_interactions(
    data = test_data,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Check that Friedman was used
  expect_true(all(results_paired$test_method == "friedman"))
  expect_true(nrow(results_paired) == n_genes)
})

test_that("detect_q_gene_interactions uses Kruskal-Wallis for paired=FALSE", {
  skip_if_not_installed("Matrix")
  
  # Create minimal unpaired data
  set.seed(42)
  n_subjects <- 5
  n_q <- 4
  n_genes <- 2
  
  data_list <- lapply(1:n_genes, function(g) {
    data.frame(
      diversity = rnorm(n_subjects * n_q, mean = 2 + g * 0.3, sd = 0.2),
      q = factor(rep(seq(0.5, 2.0, length.out = n_q), each = n_subjects)),
      paired_samples = factor(rep(1:n_subjects, n_q)),
      gene = paste0("gene_", g)
    )
  })
  
  test_data <- do.call(rbind, data_list)
  rownames(test_data) <- NULL
  
  # Run unpaired analysis
  results_unpaired <- detect_q_gene_interactions(
    data = test_data,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    paired = FALSE,
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Check that Kruskal-Wallis was used (not Friedman)
  expect_true(all(results_unpaired$test_method != "friedman"))
  expect_true(nrow(results_unpaired) == n_genes)
})

# ============================================================================
# Test 7: P-value Validity
# ============================================================================

test_that("Friedman p-values are valid (between 0 and 1)", {
  set.seed(42)
  
  for (i in 1:10) {
    # Random paired data
    n_subjects <- sample(4:8, 1)
    n_q <- sample(3:6, 1)
    
    data <- data.frame(
      entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.5),
      q = factor(rep(1:n_q, each = n_subjects)),
      subject = factor(rep(1:n_subjects, n_q))
    )
    
    result <- .tsenat_apply_friedman_test(
      data = data,
      value_col = "entropy",
      group_col = "q",
      subject_col = "subject"
    )
    
    if (!is.na(result$p_value)) {
      expect_gte(result$p_value, 0, 
                label = paste0("p_value #", i))
      expect_lte(result$p_value, 1, 
                label = paste0("p_value #", i))
    }
  }
})

# ============================================================================
# Test 8: Friedman vs Kruskal-Wallis Power Comparison
# ============================================================================

test_that("Friedman test can have better p-values than Kruskal-Wallis for paired data", {
  # Create paired data with strong between-subject effect
  set.seed(777)
  n_subjects <- 15
  n_q <- 3
  
  # Subject effects (some subjects naturally have higher entropy)
  subject_effects <- rnorm(n_subjects, mean = 0, sd = 1)
  
  # Q-value effects
  q_effects <- c(-0.5, 0, 0.5)
  
  entropy <- numeric(n_subjects * n_q)
  for (i in 1:n_subjects) {
    for (j in 1:n_q) {
      idx <- (i-1)*n_q + j
      entropy[idx] <- 2 + subject_effects[i] + q_effects[j] + rnorm(1, sd=0.1)
    }
  }
  
  data <- data.frame(
    entropy = entropy,
    q = factor(rep(1:n_q, n_subjects)),
    subject = factor(rep(1:n_subjects, each = n_q))
  )
  
  # Friedman test (paired)
  friedman_result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Kruskal-Wallis (unpaired)
  kw_result <- kruskal.test(data$entropy ~ data$q)
  
  # With strong subject effects and clear q-effects,
  # Friedman should typically be more powerful
  expect_is(friedman_result$p_value, "numeric")
  expect_is(kw_result$p.value, "numeric")
})

# ============================================================================
# Test 9: Different Data Types and Scales
# ============================================================================

test_that("Friedman works with different entropy scales", {
  # Low entropy (bounded 0-1 ish)
  data_low <- data.frame(
    entropy = c(0.1, 0.15, 0.12, 0.18, 0.22, 0.25, 0.28, 0.3, 0.32, 0.35,
                0.11, 0.14, 0.13, 0.17, 0.21, 0.24, 0.29, 0.31, 0.33, 0.34),
    q = factor(rep(1:2, 10)),
    subject = factor(rep(1:10, each = 2))
  )
  
  result_low <- .tsenat_apply_friedman_test(
    data = data_low,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_is(result_low$p_value, "numeric")
  expect_true(result_low$p_value >= 0 && result_low$p_value <= 1)
  
  # High entropy (large values)
  data_high <- data.frame(
    entropy = c(245, 250, 248, 252, 255, 260, 258, 262, 265, 270,
                244, 249, 247, 251, 254, 259, 257, 261, 264, 269),
    q = factor(rep(1:2, 10)),
    subject = factor(rep(1:10, each = 2))
  )
  
  result_high <- .tsenat_apply_friedman_test(
    data = data_high,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_is(result_high$p_value, "numeric")
  expect_true(result_high$p_value >= 0 && result_high$p_value <= 1)
})

# ============================================================================
# Test 10: Friedman Test Result Reproducibility
# ============================================================================

test_that("Friedman test results are reproducible", {
  set.seed(54321)
  
  data <- data.frame(
    entropy = rnorm(40, mean = 2.5, sd = 0.4),
    q = factor(rep(1:5, 8)),
    subject = factor(rep(1:8, each = 5))
  )
  
  result1 <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result2 <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Exact reproducibility
  expect_equal(result1$statistic, result2$statistic)
  expect_equal(result1$p_value, result2$p_value)
  expect_equal(result1$method, result2$method)
})
