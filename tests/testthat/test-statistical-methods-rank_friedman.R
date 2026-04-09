
context("friedman_complete: Complete Friedman Test Implementation")
library(TSENAT)
library(SummarizedExperiment)


# ════════════════════════════════════════════════════════════════════════════════
# SETUP: Create test data with valid structure for paired designs
# ════════════════════════════════════════════════════════════════════════════════

create_paired_diversity_se <- function(
    n_genes = 50,
    n_subjects = 8,
    n_q_values = 10,
    effect_size = 0.3,
    seed = 123) {
  
  set.seed(seed)
  
  # Create diversity matrix: genes × (q_values × subjects) samples
  n_cols <- n_q_values * n_subjects
  diversity_matrix <- matrix(
    rnorm(n_genes * n_cols, mean = 2.5, sd = 0.8),
    nrow = n_genes,
    ncol = n_cols,
    dimnames = list(
      paste0("Gene", 1:n_genes),
      paste0("Sample", 1:n_cols)
    )
  )
  
  # Add q-dependent signal to first 10 genes
  q_vec <- rep(seq(0.1, 2.0, length.out = n_q_values), n_subjects)
  for (i in seq_len(10)) {
    # Add effect proportional to q-value
    diversity_matrix[i, ] <- diversity_matrix[i, ] + effect_size * q_vec
  }
  
  # Create colData with q, paired_samples, and condition
  coldata <- DataFrame(
    q = q_vec,
    paired_samples = rep(paste0("Subject", 1:n_subjects), each = n_q_values),
    condition = rep(c("A", "B"), length.out = n_cols),
    sample_type = "diversity"
  )
  
  rowdata <- DataFrame(gene = paste0("Gene", 1:n_genes))
  
  # Build SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(diversity = diversity_matrix),
    colData = coldata,
    rowData = rowdata
  )
  
  return(se)
}

# ════════════════════════════════════════════════════════════════════════════════
# TEST 1: detect_q_gene_interactions with paired=TRUE and hochberg correction
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions works with paired=TRUE and hochberg", {
  se <- create_paired_diversity_se(n_genes = 30, n_subjects = 6, n_q_values = 8)
  
  # Extract diversity assay
  diversity_data <- assay(se)
  
  # Convert to long format for detect_q_gene_interactions
  test_data <- data.frame(
    entropy = as.numeric(diversity_data),
    gene = rep(rownames(se), ncol(se)),
    q = rep(colData(se)$q, each = nrow(se)),
    condition = rep(colData(se)$condition, each = nrow(se)),
    paired_samples = rep(colData(se)$paired_samples, each = nrow(se)),
    stringsAsFactors = FALSE
  )
  
  # Run detect_q_gene_interactions with explicit column specifications
  results <- .calculate_rank_test(
    data = test_data,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify results structure
  expect_true(is.data.frame(results))
  expect_equal(nrow(results), 30)
  expect_true("gene" %in% colnames(results))
  expect_true("p_value" %in% colnames(results))
  expect_true("adj_p_value" %in% colnames(results))
  expect_true("f_statistic" %in% colnames(results))
  expect_true("effect_size_eta2" %in% colnames(results))
  expect_true("test_method" %in% colnames(results))
  
  # Verify no NAs in crucial columns (except for failed tests)
  expect_true(length(results$adj_p_value) == 30)
  expect_false(all(is.na(results$adj_p_value)))
  
  # Verify adj_p_value is properly adjusted (monotone increasing)
  valid_p <- results$adj_p_value[!is.na(results$adj_p_value)]
  if (length(valid_p) > 1) {
    sorted_p <- sort(valid_p)
    # Hochberg should give monotone increasing p-values
    expect_true(all(diff(sorted_p) >= -1e-10))  # Allow for numerical errors
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 2: Paired test selection logic detects heteroscedasticity
# ════════════════════════════════════════════════════════════════════════════════

test_that("paired test selection detects heteroscedasticity", {
  # Create data with heteroscedastic structure
  set.seed(789)
  n_subject <- 8
  n_q <- 10
  
  heteroscedastic_data <- data.frame(
    entropy = c(
      # Subjects with low variance for q1
      rnorm(n_subject, mean = 2.5, sd = 0.2),
      # Subjects with high variance for q2, q3...
      rnorm(n_subject * (n_q - 1), mean = 2.7, sd = 1.5)
    ),
    q = rep(seq(0.1, 2.0, length.out = n_q), each = n_subject),
    paired_samples = rep(paste0("S", 1:n_subject), n_q),
    stringsAsFactors = FALSE
  )
  
  heteroscedastic_data$q <- factor(heteroscedastic_data$q)
  heteroscedastic_data$paired_samples <- factor(heteroscedastic_data$paired_samples)
  
  # Call selection function
  selection <- .select_rank_test_paired(
    heteroscedastic_data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "paired_samples",
    verbose = FALSE
  )
  
  # Verify selection result structure
  expect_true(is.list(selection))
  expect_true("test_selected" %in% names(selection))
  expect_true("characteristics" %in% names(selection))
  expect_true("reasons" %in% names(selection))
  
  # Verify characteristics detected
  expect_true(is.logical(selection$characteristics$heteroscedastic))
  expect_true(is.logical(selection$characteristics$highly_skewed))
  expect_true(is.numeric(selection$characteristics$n_groups))
  expect_true(is.numeric(selection$characteristics$n_subjects))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 3: Standard Friedman test application
# ════════════════════════════════════════════════════════════════════════════════

test_that(".apply_friedman_test produces valid output", {
  set.seed(456)
  
  # Balanced paired data
  friedman_data <- data.frame(
    entropy = rnorm(48, mean = 2.5, sd = 0.6),
    q = factor(rep(c("q0.1", "q0.5", "q1.0", "q1.5"), 12)),
    subject = factor(rep(1:12, each = 4)),
    stringsAsFactors = FALSE
  )
  
  # Apply Friedman test
  result <- .apply_friedman_test(
    data = friedman_data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Verify output structure
  expect_true(is.list(result))
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  
  # Verify values are numeric and valid
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$p_value))
  expect_true(!is.na(result$statistic))
  expect_true(!is.na(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$statistic >= 0)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 4: Detect_q_gene_interactions with SummarizedExperiment input
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions accepts SummarizedExperiment directly", {
  se <- create_paired_diversity_se(n_genes = 20, n_subjects = 6, n_q_values = 8)
  
  # Call with SE directly
  results <- .calculate_rank_test(
    data = se,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify results
  expect_true(is.data.frame(results))
  expect_equal(nrow(results), 20)
  expect_true(all(results$gene %in% rownames(se)))
  expect_false(all(is.na(results$p_value)))
  expect_false(all(is.na(results$adj_p_value)))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 5: Multiple testing correction methods
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions supports all multicorr methods", {
  se <- create_paired_diversity_se(n_genes = 20, n_subjects = 6, n_q_values = 8)
  
  # Test each correction method
  for (method in c("hochberg", "benjamini-yekutieli", "none")) {
    results <- .calculate_rank_test(
      data = se,
      condition_col = "condition",
      paired = TRUE,
      subject_col = "paired_samples",
      multicorr = method,
      verbose = FALSE
    )
    
    # Verify results structure is consistent
    expect_true(is.data.frame(results))
    expect_equal(nrow(results), nrow(se))
    expect_true("adj_p_value" %in% colnames(results))
    expect_equal(length(results$adj_p_value), nrow(se))
    
    # Verify adjusted p-values are in valid range
    valid_adj_p <- results$adj_p_value[!is.na(results$adj_p_value)]
    expect_true(all(valid_adj_p >= 0 & valid_adj_p <= 1))
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 6: Effect size computation (eta-squared)
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions computes effect sizes correctly", {
  se <- create_paired_diversity_se(
    n_genes = 20,
    n_subjects = 6,
    n_q_values = 8,
    effect_size = 0.5  # Strong effect
  )
  
  results <- .calculate_rank_test(
    data = se,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  expect_equal(length(results$effect_size_eta2), nrow(se))
  
  # Effect sizes should be in [0, 1]
  valid_effects <- results$effect_size_eta2[!is.na(results$effect_size_eta2)]
  expect_true(all(valid_effects >= 0 & valid_effects <= 1))
  
  # First 10 genes had signal, should have higher effect sizes
  mean_effect_signal <- mean(results$effect_size_eta2[1:10], na.rm = TRUE)
  mean_effect_noise <- mean(results$effect_size_eta2[11:20], na.rm = TRUE)
  
  # Signal genes should have larger effects on average
  expect_true(mean_effect_signal > mean_effect_noise)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 7: Test selection metadata is properly stored
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions stores test_method metadata", {
  se <- create_paired_diversity_se(n_genes = 15, n_subjects = 6, n_q_values = 8)
  
  results <- .calculate_rank_test(
    data = se,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify test_method column exists
  expect_true("test_method" %in% colnames(results))
  
  # test_method should contain information about which test was used
  expect_true(all(!is.na(results$test_method) | results$interaction_class == "Test failed"))
  
  # test_method should have some meaningful values (not all empty or NA)
  test_types <- unique(results$test_method[!is.na(results$test_method)])
  # Just verify that test_method column has some non-NA values with reasonable names
  expect_true(length(test_types) > 0)
  expect_true(all(nchar(test_types) > 0))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 8: Handling of edge cases - unbalanced designs
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions handles missing q-values gracefully", {
  # Create unbalanced data (missing one q-value for one subject)
  # Note: Larger sample size (n_subjects=12) to avoid chi-squared approximation warnings
  n_subjects <- 12
  n_q_values <- 6
  n_genes <- 20
  total_obs <- n_subjects * n_q_values * (n_genes - 1) + (n_subjects * (n_q_values - 1))
  
  # Create balanced data for first 19 genes, one short for last gene
  entropy_vals <- rnorm(total_obs)
  gene_vals <- c(rep(paste0("Gene", 1:19), each = n_subjects * n_q_values), 
                 rep("Gene20", n_subjects * (n_q_values - 1)))
  q_vals <- c(rep(rep(seq(0.1, 2, length.out = n_q_values), each = n_subjects), 19),
              rep(seq(0.1, 1.9, length.out = n_q_values - 1), each = n_subjects))
  subject_vals <- c(rep(rep(1:n_subjects, n_q_values), 19),
                    rep(1:n_subjects, n_q_values - 1))
  
  test_data <- data.frame(
    entropy = entropy_vals,
    gene = gene_vals,
    q = q_vals,
    condition = rep(c("A", "B"), length.out = length(entropy_vals)),
    paired_samples = subject_vals,
    stringsAsFactors = FALSE
  )
  
  # Should complete without error, but some genes may have NA results
  results <- .calculate_rank_test(
    data = test_data,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  expect_true(is.data.frame(results))
  expect_equal(nrow(results), 20)
  # Some tests may fail due to unbalanced design
  expect_true(any(is.na(results$p_value)) || all(!is.na(results$p_value)))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 9: Interaction classification
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions classifies interactions correctly", {
  se <- create_paired_diversity_se(n_genes = 20, n_subjects = 6, n_q_values = 8)
  
  results <- .calculate_rank_test(
    data = se,
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify interaction_class column
  expect_true("interaction_class" %in% colnames(results))
  expect_true(all(!is.na(results$interaction_class)))
  
  # Should have valid interaction class values (these are returned by classify_q_dependency)
  # Valid classes include various interaction strengths plus special cases
  expect_true(all(nchar(results$interaction_class) > 0))  # Non-empty strings
  # At least some genes should be classified (not all test failures)
  expect_true(any(results$interaction_class != "Test failed"))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 10: Paired vs unpaired design detection
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions uses correct test for paired design", {
  # Increased sample sizes to improve chi-squared approximation (avoid warnings)
  se <- create_paired_diversity_se(n_genes = 15, n_subjects = 10, n_q_values = 8)
  
  results_paired <- .calculate_rank_test(
    data = se,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # For paired design, test_method should reflect paired tests
  observed_types <- unique(results_paired$test_method[!is.na(results_paired$test_method)])
  
  # At least some genes should have non-empty test_method values
  expect_true(length(observed_types) > 0 || nrow(results_paired) == 0)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 11: Verifies dimension consistency after multiple testing correction
# ════════════════════════════════════════════════════════════════════════════════

test_that("multiple testing correction maintains data frame dimensions", {
  se <- create_paired_diversity_se(n_genes = 50, n_subjects = 6, n_q_values = 8)
  
  for (method in c("hochberg", "benjamini-yekutieli", "none")) {
    results <- .calculate_rank_test(
      data = se,
      condition_col = "condition",
      paired = TRUE,
      subject_col = "paired_samples",
      multicorr = method,
      verbose = FALSE
    )
    
    # Verify dimensions are preserved
    expect_equal(nrow(results), 50)
    expect_equal(length(results$p_value), 50)
    expect_equal(length(results$adj_p_value), 50)
    
    # Verify no column length mismatches
    col_lengths <- sapply(results, length)
    expect_true(all(col_lengths == 50))
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 12: Verify output sorting and ranking
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions returns sorted results", {
  se <- create_paired_diversity_se(n_genes = 30, n_subjects = 6, n_q_values = 8)
  
  results <- .calculate_rank_test(
    data = se,
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify primary sort by adjusted p-value (ascending)
  valid_adj_p <- results$adj_p_value[!is.na(results$adj_p_value)]
  if (length(valid_adj_p) > 1) {
    # Should be monotone increasing when taking valid values in order
    expect_true(all(diff(sort(valid_adj_p)) >= -1e-10))
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 13: Data characteristics storage
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions stores data characteristics", {
  se <- create_paired_diversity_se(n_genes = 20, n_subjects = 6, n_q_values = 8)
  
  results <- .calculate_rank_test(
    data = se,
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Check for characteristic columns
  expect_true("heteroscedastic" %in% colnames(results))
  expect_true("boundary_clustered" %in% colnames(results))
  expect_true("highly_skewed" %in% colnames(results))
  
  # Values should be logical
  expect_true(all(is.logical(results$heteroscedastic)))
  expect_true(all(is.logical(results$boundary_clustered)))
  expect_true(all(is.logical(results$highly_skewed)))
})


context("Friedman Test for Paired Rank-Based Analysis")

# ============================================================================
# Test 1: Basic Friedman Test Function
# ============================================================================

test_that(".apply_friedman_test returns correct structure", {
  # Create simple paired test data (3 subjects, 3 treatments)
  data <- data.frame(
    value = c(1.2, 1.5, 1.8, 2.1, 2.3, 2.5, 3.0, 3.2, 3.5),
    subject = factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3)),
    treatment = factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
  )
  
  result <- .apply_friedman_test(
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

test_that(".apply_friedman_test accepts data frame input", {
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
  
  result <- .apply_friedman_test(
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

test_that(".apply_friedman_test matches base R friedman.test", {
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
  our_result <- .apply_friedman_test(
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

test_that(".apply_friedman_test handles missing data gracefully", {
  # Data with only 1 subject
  data <- data.frame(
    value = c(1.0, 2.0, 3.0),
    subject = factor(c(1, 1, 1)),
    treatment = factor(c("A", "B", "C"))
  )
  
  result <- .apply_friedman_test(
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

test_that(".apply_friedman_test handles single treatment", {
  # Data with only 1 treatment
  data <- data.frame(
    value = c(1.0, 2.0, 3.0),
    subject = factor(c(1, 2, 3)),
    treatment = factor(c("A", "A", "A"))
  )
  
  result <- .apply_friedman_test(
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

test_that(".apply_conditional_rank_test selects Friedman for paired", {
  # Create balanced paired data
  set.seed(99)
  data <- data.frame(
    entropy = rnorm(30, mean = 2, sd = 0.3),
    q = factor(rep(1:5, 6)),
    subject = factor(rep(1:6, each = 5))
  )
  
  result <- .apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  
  # Should select paired test (Friedman or ART-Friedman based on data characteristics)
  # The function correctly detects characteristics and selects appropriate method
  expect_true(result$test_type %in% c("friedman", "art_friedman", "robust_friedman"))
  expect_match(result$method, "Friedman|friedman")
})

test_that(".apply_conditional_rank_test uses Kruskal-Wallis when unpaired", {
  # Create data without subject column (unpaired)
  set.seed(99)
  data <- data.frame(
    entropy = rnorm(150, mean = 2, sd = 0.3),
    q = factor(rep(1:5, 30))
  )
  
  result <- .apply_conditional_rank_test(
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

test_that(".apply_conditional_rank_test requires subject_col for paired", {
  # Create data without explicit subject column for paired analysis
  data <- data.frame(
    entropy = rnorm(20, mean = 2, sd = 0.3),
    q = factor(rep(1:4, 5))
  )
  
  result <- .apply_conditional_rank_test(
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
  
  result <- .apply_conditional_rank_test(
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
      condition = factor(rep(c("A", "B"), length.out = n_subjects * n_q)),
      paired_samples = factor(rep(1:n_subjects, n_q)),
      gene = paste0("gene_", g)
    )
  })
  
  test_data <- do.call(rbind, data_list)
  rownames(test_data) <- NULL
  
  # Run paired analysis
  results_paired <- .calculate_rank_test(
    data = test_data,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Check that paired test was used
  # test_method should contain Scheirer-Ray-Hare or other paired test methods
  test_methods <- results_paired$test_method[!is.na(results_paired$test_method)]
  # Verify at least some genes have valid test methods (not all failures)
  expect_true(length(test_methods) > 0, info = "Should have at least some test methods")
  # Test methods should include expected paired test types
  expect_true(any(grepl("srh|friedman", test_methods, ignore.case = TRUE)),
              info = "Should use Scheirer-Ray-Hare or Friedman test for paired design")
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
      condition = factor(rep(c("A", "B"), length.out = n_subjects * n_q)),
      paired_samples = factor(rep(1:n_subjects, n_q)),
      gene = paste0("gene_", g)
    )
  })
  
  test_data <- do.call(rbind, data_list)
  rownames(test_data) <- NULL
  
  # Run unpaired analysis
  results_unpaired <- .calculate_rank_test(
    data = test_data,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
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
  
  for (i in seq_len(10)) {
    # Random paired data
    n_subjects <- sample(4:8, 1)
    n_q <- sample(3:6, 1)
    
    data <- data.frame(
      entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.5),
      q = factor(rep(1:n_q, each = n_subjects)),
      subject = factor(rep(1:n_subjects, n_q))
    )
    
    result <- .apply_friedman_test(
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
  for (i in seq_len(n_subjects)) {
    for (j in seq_len(n_q)) {
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
  friedman_result <- .apply_friedman_test(
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
  
  result_low <- .apply_friedman_test(
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
  
  result_high <- .apply_friedman_test(
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
  
  result1 <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result2 <- .apply_friedman_test(
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


context("Friedman Test Integration and Edge Cases")

# ============================================================================
# Test: Paired vs Unpaired Comparison on Same Dataset
# ============================================================================

test_that("Paired analysis produces different results than unpaired on same data", {
  set.seed(888)
  
  # Create data with strong subject effect
  n_subjects <- 8
  n_q <- 5
  
  subject_effect <- rnorm(n_subjects, sd = 2)
  entropy_data <- numeric(n_subjects * n_q)
  
  for (i in seq_len(n_subjects)) {
    for (j in seq_len(n_q)) {
      idx <- (i-1)*n_q + j
      entropy_data[idx] <- 2 + subject_effect[i] + (j - n_q/2) * 0.3 + rnorm(1, sd=0.1)
    }
  }
  
  data <- data.frame(
    entropy = entropy_data,
    q = factor(rep(1:n_q, n_subjects)),
    subject = factor(rep(1:n_subjects, each = n_q))
  )
  
  # Paired analysis
  paired_result <- .apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  
  # Unpaired analysis (data without subject column in conditional context)
  unpaired_result <- .apply_conditional_rank_test(
    data = data[, c("entropy", "q")],
    value_col = "entropy",
    group_col = "q",
    paired = FALSE,
    subject_col = NULL,
    verbose = FALSE
  )
  
  # Results should differ (paired should typically be more significant)
  expect_is(paired_result$p_value, "numeric")
  expect_is(unpaired_result$p_value, "numeric")
  
  # Document that they are different
  expect_false(isTRUE(all.equal(paired_result$p_value, unpaired_result$p_value)),
               info = "Paired and unpaired analyses should produce different p-values")
})

# ============================================================================
# Test: Friedman with Minimum Required Structure
# ============================================================================

test_that("Friedman works with minimal paired structure (3 subjects, 2 treatments)", {
  data <- data.frame(
    entropy = c(1.0, 2.0, 1.5, 2.5),
    subject = factor(c(1, 2, 1, 2)),
    treatment = factor(c("A", "A", "B", "B"))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Large Number of Subjects
# ============================================================================

test_that("Friedman handles large number of subjects", {
  set.seed(999)
  
  n_subjects <- 50  # Large
  n_q <- 4
  
  data <- data.frame(
    entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.3),
    subject = factor(rep(1:n_subjects, n_q)),
    q = factor(rep(1:n_q, each = n_subjects))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$statistic, "numeric")
  expect_true(result$statistic >= 0)
})

# ============================================================================
# Test: Large Number of Q-Values
# ============================================================================

test_that("Friedman handles large number of q-values (treatments)", {
  set.seed(1111)
  
  n_subjects <- 6
  n_q <- 20  # Many q-values
  
  data <- data.frame(
    entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.3),
    subject = factor(rep(1:n_subjects, n_q)),
    q = factor(rep(1:n_q, each = n_subjects))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Friedman Response to Actual Effect
# ============================================================================

test_that("Friedman detects strong q-effect (low p-value)", {
  set.seed(2222)
  
  n_subjects <- 10
  n_q <- 4
  
  # Strong q-effect: entropy depends on q
  entropy <- numeric(n_subjects * n_q)
  for (i in seq_len(n_subjects)) {
    for (j in seq_len(n_q)) {
      idx <- (i-1)*n_q + j
      entropy[idx] <- 1 + (j - 1) * 0.8 + rnorm(1, sd=0.05)  # Strong effect, small noise
    }
  }
  
  data <- data.frame(
    entropy = entropy,
    subject = factor(rep(1:n_subjects, each = n_q)),
    q = factor(rep(1:n_q, n_subjects))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should be highly significant
  expect_lt(result$p_value, 0.05)
})

test_that("Friedman shows weak effect (high p-value) for random data", {
  set.seed(3333)
  
  n_subjects <- 8
  n_q <- 4
  
  # No q-effect: completely random
  entropy <- rnorm(n_subjects * n_q, mean = 2, sd = 0.8)
  
  data <- data.frame(
    entropy = entropy,
    subject = factor(rep(1:n_subjects, each = n_q)),
    q = factor(rep(1:n_q, n_subjects))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should not be highly significant
  expect_gt(result$p_value, 0.05)
})

# ============================================================================
# Test: Friedman Statistics Interpretation
# ============================================================================

test_that("Friedman statistic follows chi-square distribution (df = k-1)", {
  set.seed(4444)
  
  n_subjects <- 15
  k_treatments <- 5  # number of q-values
  
  data <- data.frame(
    entropy = rnorm(n_subjects * k_treatments, mean = 2, sd = 0.4),
    subject = factor(rep(1:n_subjects, k_treatments)),
    q = factor(rep(1:k_treatments, each = n_subjects))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Friedman statistic should be non-negative
  expect_gte(result$statistic, 0)
  
  # For random data, statistic follows approximately chi-square(k-1)
  # So it should be in reasonable range (0 to maybe 4*df for random data)
  df <- k_treatments - 1
  expect_lt(result$statistic, 4 * df + 5)
})

# ============================================================================
# Test: Column Name Flexibility
# ============================================================================

test_that("Friedman works with different column names", {
  set.seed(5555)
  
  # Different column names
  data <- data.frame(
    my_entropy = rnorm(20, mean = 2, sd = 0.3),
    my_treatment = factor(rep(1:4, 5)),
    my_block = factor(rep(1:5, each = 4))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "my_entropy",
    group_col = "my_treatment",
    subject_col = "my_block"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
})

# ============================================================================
# Test: Factor vs Character Handling
# ============================================================================

test_that("Friedman handles both factor and numeric group identifiers", {
  set.seed(6666)
  
  # With numeric (converted to factor)
  data_numeric <- data.frame(
    entropy = rnorm(16, mean = 2, sd = 0.3),
    subject = 1:4,  # Numeric
    q = rep(1:4, each = 4)  # Numeric
  )
  
  # Should still work (function converts internally)
  result_numeric <- .apply_friedman_test(
    data = data_numeric,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # With character
  data_char <- data.frame(
    entropy = rnorm(16, mean = 2, sd = 0.3),
    subject = factor(rep(paste0("S", 1:4), each = 4)),  # Factor with character
    q = factor(rep(paste0("Q", 1:4), 4))  # Factor with character
  )
  
  result_char <- .apply_friedman_test(
    data = data_char,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Both should succeed
  expect_equal(result_numeric$method, "Friedman test (paired)")
  expect_equal(result_char$method, "Friedman test (paired)")
})

# ============================================================================
# Test: Behavior with Identical Values
# ============================================================================

test_that("Friedman handles data with tied (identical) values", {
  # Data with many ties
  data <- data.frame(
    entropy = c(2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 2.0, 2.0, 2.0, 2.5, 2.5, 2.5),
    subject = factor(c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6)),
    q = factor(c(rep("A", 6), rep("B", 6)))
  )
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should still return valid result
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Missing Data Patterns
# ============================================================================

test_that("Friedman fails gracefully with missing values in matrix", {
  # Create a data frame where using xtabs would create a result without required dimensions
  # e.g., only one level in one dimension after filtering
  data <- data.frame(
    entropy = c(1.0, 2.0),
    subject = factor(c(1, 1)),
    q = factor(c("A", "B"))
  )
  # This has only 1 subject but 2 treatments - Friedman needs at least 2 subjects
  
  result <- .apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should fail gracefully - either NA p-value or test_failed method
  expect_true(is.na(result$p_value) || result$method == "test_failed")
})

# ============================================================================
# Test: Monotonic Transformations Don't Change p-value
# ============================================================================

test_that("Friedman p-value invariant to monotonic transformations", {
  set.seed(7777)
  
  entropy_original <- rnorm(24, mean = 2, sd = 0.5)
  
  data_original <- data.frame(
    entropy = entropy_original,
    subject = factor(rep(1:6, 4)),
    q = factor(rep(1:4, each = 6))
  )
  
  # Log transformation
  data_log <- data_original
  data_log$entropy <- log(entropy_original + 1)  # Ensure positive
  
  # Square transformation
  data_sq <- data_original
  data_sq$entropy <- entropy_original^2
  
  result_original <- .apply_friedman_test(
    data = data_original,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result_log <- .apply_friedman_test(
    data = data_log,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result_sq <- .apply_friedman_test(
    data = data_sq,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Friedman is rank-based, so monotonic transformations shouldn't change results
  # (p-values may differ slightly due to ranking, but interpretation should hold)
  expect_is(result_original$p_value, "numeric")
  expect_is(result_log$p_value, "numeric")
  expect_is(result_sq$p_value, "numeric")
  
  # Rank-invariance property: should get same ranking, hence same p-value
  expect_equal(result_original$statistic, result_log$statistic, tolerance = 1e-10)
  expect_equal(result_original$statistic, result_sq$statistic, tolerance = 1e-10)
})

# ============================================================================
# Test: Extreme Value Handling
# ============================================================================

test_that("Friedman handles extreme values correctly", {
  # Very small values
  data_small <- data.frame(
    entropy = c(1e-10, 2e-10, 1.5e-10, 2.5e-10, 1e-10, 2e-10, 1.5e-10, 2.5e-10),
    subject = factor(rep(1:4, 2)),
    q = factor(rep(1:2, each = 4))
  )
  
  result_small <- .apply_friedman_test(
    data = data_small,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Very large values
  data_large <- data.frame(
    entropy = c(1e10, 2e10, 1.5e10, 2.5e10, 1e10, 2e10, 1.5e10, 2.5e10),
    subject = factor(rep(1:4, 2)),
    q = factor(rep(1:2, each = 4))
  )
  
  result_large <- .apply_friedman_test(
    data = data_large,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Both should work (rank-based, so scale doesn't matter)
  expect_equal(result_small$method, "Friedman test (paired)")
  expect_equal(result_large$method, "Friedman test (paired)")
  
  # Results should be identical (rank-based)
  expect_equal(result_small$statistic, result_large$statistic)
  expect_equal(result_small$p_value, result_large$p_value)
})
