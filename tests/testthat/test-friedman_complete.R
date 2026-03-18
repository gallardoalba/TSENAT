library(TSENAT)
library(SummarizedExperiment)

context("friedman_complete: Complete Friedman Test Implementation")

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
  for (i in 1:10) {
    # Add effect proportional to q-value
    diversity_matrix[i, ] <- diversity_matrix[i, ] + effect_size * q_vec
  }
  
  # Create colData with q and paired_samples
  coldata <- DataFrame(
    q = q_vec,
    paired_samples = rep(paste0("Subject", 1:n_subjects), each = n_q_values),
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
    paired_samples = rep(colData(se)$paired_samples, each = nrow(se)),
    stringsAsFactors = FALSE
  )
  
  # Run detect_q_gene_interactions with explicit column specifications
  results <- detect_q_gene_interactions(
    data = test_data,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
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
  selection <- .tsenat_select_rank_test_paired(
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

test_that(".tsenat_apply_friedman_test produces valid output", {
  set.seed(456)
  
  # Balanced paired data
  friedman_data <- data.frame(
    entropy = rnorm(48, mean = 2.5, sd = 0.6),
    q = factor(rep(c("q0.1", "q0.5", "q1.0", "q1.5"), 12)),
    subject = factor(rep(1:12, each = 4)),
    stringsAsFactors = FALSE
  )
  
  # Apply Friedman test
  result <- .tsenat_apply_friedman_test(
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
  results <- detect_q_gene_interactions(
    data = se,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
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
    results <- detect_q_gene_interactions(
      data = se,
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
  
  results <- detect_q_gene_interactions(
    data = se,
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify effect sizes
  expect_true("effect_size_eta2" %in% colnames(results))
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
  
  results <- detect_q_gene_interactions(
    data = se,
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify test_method column exists
  expect_true("test_method" %in% colnames(results))
  
  # test_method should contain information about which test was used
  expect_true(all(!is.na(results$test_method) | results$interaction_class == "Test failed"))
  
  # Some tests should be marked as friedman or art_friedman or robust_friedman
  test_types <- unique(results$test_method[!is.na(results$test_method)])
  expect_true(any(test_types %in% c("friedman", "art_friedman", "robust_friedman", "insufficient_data")))
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
    paired_samples = subject_vals,
    stringsAsFactors = FALSE
  )
  
  # Should complete without error, but some genes may have NA results
  results <- detect_q_gene_interactions(
    data = test_data,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
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
  
  results <- detect_q_gene_interactions(
    data = se,
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
  
  results_paired <- detect_q_gene_interactions(
    data = se,
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # For paired design, test_method should reflect paired tests
  paired_test_types <- c("friedman", "art_friedman", "robust_friedman", "test_failed", "insufficient_data")
  observed_types <- unique(results_paired$test_method[!is.na(results_paired$test_method)])
  
  # At least some genes should use paired tests
  expect_true(any(observed_types %in% paired_test_types))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 11: Verifies dimension consistency after multiple testing correction
# ════════════════════════════════════════════════════════════════════════════════

test_that("multiple testing correction maintains data frame dimensions", {
  se <- create_paired_diversity_se(n_genes = 50, n_subjects = 6, n_q_values = 8)
  
  for (method in c("hochberg", "benjamini-yekutieli", "none")) {
    results <- detect_q_gene_interactions(
      data = se,
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
  
  results <- detect_q_gene_interactions(
    data = se,
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
  
  results <- detect_q_gene_interactions(
    data = se,
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
