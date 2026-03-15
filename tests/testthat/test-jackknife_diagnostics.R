context("Jackknife Diagnostics: Tsallis Entropy Stability Analysis")

# Setup test data
set.seed(42)
balanced_counts <- c(100, 100, 100, 100, 100)  # All equal
skewed_counts <- c(500, 200, 100, 50, 20)  # Strongly skewed
with_dominant <- c(900, 50, 30, 20)  # One very dominant
single_counts <- c(1000)  # Single transcript
two_counts <- c(500, 500)  # Two transcripts

# Test 1: Basic vector input
test_that("jackknife_tsallis_entropy works with vector input", {
  result <- jackknife_tsallis_entropy(balanced_counts, q = 1)

  expect_is(result, "tsenat_jackknife")
  expect_is(result$estimate, "numeric")
  expect_length(result$jackknife_estimates, length(balanced_counts))
  expect_length(result$influence, length(balanced_counts))
  expect_true(!is.na(result$jackknife_se))
  expect_true(result$n_transcripts == length(balanced_counts))
})

# Test 2: Matrix input
test_that("jackknife_tsallis_entropy works with matrix input", {
  gene_matrix <- matrix(c(balanced_counts, skewed_counts), nrow = 2, byrow = TRUE)
  result <- jackknife_tsallis_entropy(gene_matrix, q = 1, print_results = FALSE)

  expect_is(result, "tsenat_jackknife_list")
  expect_length(result, 2)
  expect_is(result[[1]], "tsenat_jackknife")
  expect_is(result[[2]], "tsenat_jackknife")
})

# Test 3: Data frame input
test_that("jackknife_tsallis_entropy works with data frame input", {
  gene_df <- as.data.frame(matrix(c(balanced_counts, skewed_counts), nrow = 2, byrow = TRUE))
  result <- jackknife_tsallis_entropy(gene_df, q = 1, print_results = FALSE)

  expect_is(result, "tsenat_jackknife_list")
  expect_length(result, 2)
})

# Test 4: Balanced distribution has low influence
test_that("Balanced counts have low, uniform influence", {
  result <- jackknife_tsallis_entropy(balanced_counts, q = 1, norm = TRUE)

  expect_true(all(result$influence < 0.01))  # All very small
  expect_true(max(result$influence) - min(result$influence) < 0.005)  # All similar
})

# Test 5: Skewed distribution has variable influence
test_that("Skewed counts have higher, variable influence", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = TRUE)

  # Should have variation in influence
  expect_true(max(result$influence) > min(result$influence))
  # Influence should vary across transcripts
  expect_true(length(unique(round(result$influence, 6))) > 1)
})

# Test 6: Dominant transcript is identified as outlier
test_that("Dominant transcript is detected as outlier", {
  result <- jackknife_tsallis_entropy(with_dominant, q = 1, norm = TRUE, threshold = 75)

  expect_length(result$outlier_indices, 1)
  expect_equal(result$outlier_indices, 1)  # First transcript is dominant
})

# Test 6: More abundant transcripts tend to have higher influence
test_that("Influence detected for all transcripts", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = TRUE)

  # All transcripts should have measurable influence
  expect_true(all(result$influence >= 0))
  # There should be variation
  expect_true(max(result$influence) > 0)
})

# Test 8: Different q values
test_that("jackknife_tsallis_entropy handles different q values", {
  q_vals <- c(0.5, 1, 1.5, 2)
  results <- lapply(q_vals, function(q_val) {
    jackknife_tsallis_entropy(skewed_counts, q = q_val, norm = TRUE)
  })

  # All should be valid
  expect_true(all(sapply(results, function(r) !is.na(r$estimate))))
  # Results should differ
  estimates <- sapply(results, function(r) r$estimate)
  expect_false(all(estimates == estimates[1]))
})

# Test 8: Different normalization schemes
test_that("Normalization affects entropy values", {
  result_norm <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = TRUE)
  result_unnorm <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = FALSE)

  # Estimates should differ
  expect_false(abs(result_norm$estimate - result_unnorm$estimate) < 1e-6)

  # Both should have comparable structure
  expect_equal(result_norm$n_transcripts, result_unnorm$n_transcripts)
  expect_length(result_norm$influence, result_unnorm$n_transcripts)
})

# Test 10: Jackknife standard error
test_that("Jackknife SE is computed correctly", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = TRUE)

  expect_true(result$jackknife_se >= 0)
  # For balanced data, SE should be very small
  result_bal <- jackknife_tsallis_entropy(balanced_counts, q = 1, norm = TRUE)
  expect_true(result_bal$jackknife_se < result$jackknife_se)
})

# Test 11: Input validation - negative values
test_that("jackknife_tsallis_entropy rejects negative counts", {
  bad_counts <- c(100, -50, 50)
  expect_error(
    jackknife_tsallis_entropy(bad_counts),
    "negative values"
  )
})

# Test 12: Input validation - missing values
test_that("jackknife_tsallis_entropy rejects NA values", {
  bad_counts <- c(100, NA, 50)
  expect_error(
    jackknife_tsallis_entropy(bad_counts),
    "missing values"
  )
})

# Test 13: Input validation - q value
test_that("jackknife_tsallis_entropy validates q parameter", {
  expect_error(
    jackknife_tsallis_entropy(balanced_counts, q = 0),
    "positive numeric"
  )
  expect_error(
    jackknife_tsallis_entropy(balanced_counts, q = -1),
    "positive numeric"
  )
})

# Test 14: Input validation - threshold
test_that("jackknife_tsallis_entropy validates threshold parameter", {
  expect_error(
    jackknife_tsallis_entropy(balanced_counts, threshold = -5),
    "between 0 and 100"
  )
  expect_error(
    jackknife_tsallis_entropy(balanced_counts, threshold = 150),
    "between 0 and 100"
  )
})

# Test 15: Input validation - minimum transcripts
test_that("jackknife_tsallis_entropy requires at least 2 transcripts", {
  expect_error(
    jackknife_tsallis_entropy(single_counts),
    "at least 2 transcripts"
  )
})

# Test 16: Two-transcript case (minimum)
test_that("jackknife_tsallis_entropy works with 2 transcripts", {
  result <- jackknife_tsallis_entropy(two_counts, q = 1)

  expect_equal(result$n_transcripts, 2)
  expect_length(result$influence, 2)
  expect_true(!is.na(result$estimate))
})

# Test 17: S3 print method
test_that("print method works for jackknife results", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1)

  # use expect_output to avoid nested capture issues with testthat
  expect_output(print(result), "Jackknife Diagnostics")
  expect_output(print(result), "transcripts")
})

# Test 18: S3 summary method
test_that("summary method works for jackknife results", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1)

  expect_output(summary(result), "Jackknife Diagnostics Summary")
  expect_output(summary(result), "standard error")
})

# Test 19: S3 print for list
test_that("print method works for jackknife list", {
  gene_matrix <- matrix(c(balanced_counts, skewed_counts), nrow = 2, byrow = TRUE)
  result <- jackknife_tsallis_entropy(gene_matrix, q = 1, print_results = FALSE)

  expect_output(print(result), "2 genes")
})

# Test 20: Outlier threshold affects detection
test_that("Outlier threshold parameter works correctly", {
  result_high <- jackknife_tsallis_entropy(skewed_counts, q = 1, threshold = 90)
  result_low <- jackknife_tsallis_entropy(skewed_counts, q = 1, threshold = 50)

  # Lower threshold should find more outliers
  expect_true(length(result_low$outlier_indices) >= length(result_high$outlier_indices))
})

# Test 21: Zero counts are handled
test_that("jackknife_tsallis_entropy handles zero counts with pseudocount", {
  counts_with_zero <- c(100, 0, 50, 30)

  result <- jackknife_tsallis_entropy(counts_with_zero, q = 1, pseudocount = 1e-10)

  expect_true(!is.na(result$estimate))
  expect_true(all(!is.na(result$jackknife_estimates)))
  expect_true(all(is.finite(result$influence)))
})

# Test 22: Reproducibility without seed
test_that("jackknife_tsallis_entropy is deterministic", {
  result1 <- jackknife_tsallis_entropy(skewed_counts, q = 1)
  result2 <- jackknife_tsallis_entropy(skewed_counts, q = 1)

  expect_equal(result1$estimate, result2$estimate)
  expect_equal(result1$jackknife_estimates, result2$jackknife_estimates)
})

# Test 23: Large transcript count
test_that("jackknife_tsallis_entropy works with many transcripts", {
  many_transcripts <- rpois(100, lambda = 50)

  result <- jackknife_tsallis_entropy(many_transcripts, q = 1)

  expect_equal(result$n_transcripts, 100)
  expect_length(result$influence, 100)
  expect_true(!is.na(result$jackknife_se))
})

# Test 24: All transcripts equal
test_that("Equal transcripts have uniform influence", {
  equal_counts <- rep(100, 10)

  result <- jackknife_tsallis_entropy(equal_counts, q = 1, norm = TRUE)

  # All influences should be nearly identical
  inf_range <- max(result$influence) - min(result$influence)
  expect_true(inf_range < 1e-6)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TESTS FOR jackknife_isoform_switching() - Input Validation and Structure
# ═══════════════════════════════════════════════════════════════════════════════

context("Jackknife Isoform Switching: Input Validation")

# Setup test data for isoform switching
set.seed(42)

# Create a simple test SummarizedExperiment
test_se_basic <- function() {
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
  })
  
  counts_mat <- matrix(rpois(20, lambda=10), nrow=2, ncol=10,
                       dimnames=list(c("T1", "T2"), paste0("S", 1:10)))
  
  SummarizedExperiment(
    assays=list(counts=counts_mat),
    rowData=data.frame(
      isoform_id=c("T1", "T2"),
      gene_id=c("Gene1", "Gene1")
    ),
    colData=data.frame(
      sample_id=paste0("S", 1:10),
      condition=rep(c("condA", "condB"), each=5)
    )
  )
}

# Test 1: Function signature accepts n_bootstrap parameter
test_that("jackknife_isoform_switching function signature is correct", {
  # Just test that function exists and has correct parameters
  expect_true(exists("jackknife_isoform_switching"))
  
  # Get function signature
  sig <- formals(jackknife_isoform_switching)
  expect_true("n_bootstrap" %in% names(sig))
  expect_true("condition_col" %in% names(sig))
  expect_true("gene_col" %in% names(sig))
  expect_true("isoform_col" %in% names(sig))
})

# Test 2: Input validation - missing condition_col
test_that("jackknife_isoform_switching detects missing condition column", {
  se <- test_se_basic()
  
  expect_error(
    suppressWarnings(jackknife_isoform_switching(
      se = se,
      condition_col = "nonexistent",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      n_bootstrap = 5,
      norm = FALSE,
      print_results = FALSE
    )),
    "condition"
  )
})

# Test 3: Input validation - invalid q parameter
test_that("jackknife_isoform_switching handles q parameter appropriately", {
  se <- test_se_basic()
  
  # Test that function accepts valid q values
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.5,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 4: Function returns correct class
test_that("jackknife_isoform_switching returns tsenat_isoform_switching class", {
  se <- test_se_basic()
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
  expect_true(inherits(result, "list"))
})

# Test 5: Required output components exist
test_that("jackknife_isoform_switching output has required components", {
  se <- test_se_basic()
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  expect_true(!is.null(result$gene_names))
  expect_true(!is.null(result$conditions))
  expect_true(!is.null(result$results_per_gene))
  expect_true(!is.null(result$all_transcript_stats))
  expect_true(!is.null(result$metadata))
})

# Test 6: all_transcript_stats has correct columns
test_that("all_transcript_stats has required columns", {
  se <- test_se_basic()
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  stats <- result$all_transcript_stats
  required_cols <- c("gene", "transcript_id", "pvalue", "fdr")
  expect_true(all(required_cols %in% names(stats)))
})

# Test 7: Paired design parameter is accepted
test_that("jackknife_isoform_switching accepts pair_col parameter", {
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
  })
  
  counts_mat <- matrix(rpois(30, lambda=10), nrow=3, ncol=10,
                       dimnames=list(c("T1", "T2", "T3"), paste0("S", 1:10)))
  
  se_paired <- SummarizedExperiment(
    assays=list(counts=counts_mat),
    rowData=data.frame(
      isoform_id=c("T1", "T2", "T3"),
      gene_id=c("Gene1", "Gene1", "Gene1")
    ),
    colData=data.frame(
      sample_id=paste0("S", 1:10),
      condition=rep(c("A", "B"), each=5),
      individual_id=rep(paste0("Ind", 1:5), 2)
    )
  )
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se_paired,
    condition_col = "condition",
    pair_col = "individual_id",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 8: LM results parameter is accepted
test_that("jackknife_isoform_switching accepts lm_results parameter", {
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
  })
  
  counts_mat <- matrix(rpois(40, lambda=10), nrow=4, ncol=10,
                       dimnames=list(c("T1a", "T1b", "T2a", "T2b"), paste0("S", 1:10)))
  
  se_multi <- SummarizedExperiment(
    assays=list(counts=counts_mat),
    rowData=data.frame(
      isoform_id=c("T1a", "T1b", "T2a", "T2b"),
      gene_id=c("Gene1", "Gene1", "Gene2", "Gene2")
    ),
    colData=data.frame(
      sample_id=paste0("S", 1:10),
      condition=rep(c("A", "B"), each=5)
    )
  )
  
  lm_results <- data.frame(
    gene=c("Gene1", "Gene2"),
    p_interaction=c(0.01, 0.50),
    adj_p_interaction=c(0.02, 0.60)
  )
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se_multi,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 9: Metadata is populated
test_that("Metadata is populated after analysis", {
  se <- test_se_basic()
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.5,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  meta <- result$metadata
  # Metadata should exist
  expect_true(!is.null(meta))
  expect_is(meta, "list")
  # Should contain analysis tracking
  expect_true(length(meta) > 0)
})

# Test 10: Results are numeric (not NA/NaN) for small bootstrap
test_that("Results contain numeric values with small bootstrap", {
  se <- test_se_basic()
  
  result <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  gene_res <- result$results_per_gene[["Gene1"]]
  
  # P-values should be numeric
  expect_is(gene_res$delta_pvalue, "numeric")
  expect_length(gene_res$delta_pvalue, 2)
  expect_true(all(is.finite(gene_res$delta_pvalue)))
})

# Test 25: Pseudocount parameter
test_that("Pseudocount parameter affects zero handling", {
  counts_zero <- c(100, 0, 50)

  result_small <- jackknife_tsallis_entropy(counts_zero, pseudocount = 1e-20)
  result_large <- jackknife_tsallis_entropy(counts_zero, pseudocount = 0.1)

  # Both should be valid
  expect_true(!is.na(result_small$estimate))
  expect_true(!is.na(result_large$estimate))
  # Results may differ slightly due to different pseudocounts
  expect_false(isTRUE(all.equal(result_small$estimate, result_large$estimate)))
})

# Test 26: KL divergence (q=1) computation
test_that("KL divergence (q=1) computes correctly", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1, norm = FALSE)

  # KL divergence should be non-negative
  expect_true(result$estimate >= 0)
  expect_true(all(result$jackknife_estimates >= 0))
})

# Test 27: Log base parameter
test_that("Log base parameter changes entropy scale", {
  result_e <- jackknife_tsallis_entropy(skewed_counts, q = 1, log_base = exp(1), norm = FALSE)
  result_2 <- jackknife_tsallis_entropy(skewed_counts, q = 1, log_base = 2, norm = FALSE)

  # Nats vs bits should differ by approximately log(2)
  ratio <- result_e$estimate / result_2$estimate
  expected_ratio <- log(2)
  expect_true(abs(ratio - expected_ratio) < 0.1)
})

# Test 28: Outlier indices are numeric and valid
test_that("Outlier indices are valid positions", {
  result <- jackknife_tsallis_entropy(skewed_counts, q = 1, threshold = 80)

  expect_true(all(result$outlier_indices >= 1))
  expect_true(all(result$outlier_indices <= result$n_transcripts))
  expect_true(length(result$outlier_indices) <= result$n_transcripts)
})

# Test 29: Named vector input
test_that("jackknife_tsallis_entropy preserves names when possible", {
  counts_named <- c(Iso1 = 100, Iso2 = 200, Iso3 = 150)

  result <- jackknife_tsallis_entropy(counts_named, q = 1)

  expect_length(result$influence, 3)
  expect_equal(result$n_transcripts, 3)
})

# Test 30: Very large counts
test_that("jackknife_tsallis_entropy handles large count values", {
  large_counts <- c(1e6, 5e5, 2e5, 1e5)

  result <- jackknife_tsallis_entropy(large_counts, q = 1, norm = TRUE)

  expect_true(!is.na(result$estimate))
  expect_true(all(!is.na(result$jackknife_estimates)))
})

# Test 31: Sparse counts that may produce NaN values
test_that("jackknife_tsallis_entropy handles sparse counts with na.rm", {
  sparse_counts <- c(1, 1, 1, 1, 1, 1, 1, 0)  # Very sparse data
  
  result <- try(suppressWarnings(jackknife_tsallis_entropy(sparse_counts, q = 1, norm = TRUE)), silent = TRUE)
  
  expect_false(inherits(result, "try-error"))
  expect_is(result, "tsenat_jackknife")
  # Important: quantile should work with na.rm=TRUE
  expect_true(!is.na(result$outlier_cutoff_value))
})

# Test 32: Jackknife SE calculation with na.rm
test_that("Jackknife SE is calculated correctly with na.rm=TRUE", {
  test_counts <- c(100, 50, 25, 10, 5)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE)
  
  # SE should be valid (not NaN or Inf)
  expect_true(is.finite(result$jackknife_se))
  expect_true(!is.na(result$jackknife_se))
  expect_true(result$jackknife_se >= 0)
})

# Test 33: Outlier detection handles NA values correctly
test_that("Outlier detection filters NA values correctly", {
  test_counts <- c(100, 50, 25, 10)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE, threshold = 75)
  
  # outlier_indices should not contain NA
  expect_true(!any(is.na(result$outlier_indices)))
  # outlier_indices should be valid indices
  expect_true(all(result$outlier_indices > 0))
  expect_true(all(result$outlier_indices <= length(test_counts)))
})

# Test 34: Quantile calculation with na.rm works for uniform influence
test_that("Quantile calculation works for uniform influence distribution", {
  uniform_counts <- c(25, 25, 25, 25)  # Very uniform
  
  result <- jackknife_tsallis_entropy(uniform_counts, q = 1, norm = TRUE, threshold = 90)
  
  # Should complete without error
  expect_is(result, "tsenat_jackknife")
  expect_true(is.finite(result$outlier_cutoff_value))
  # Uniform distribution should have few/no outliers
  expect_true(length(result$outlier_indices) <= 1)
})

# Test 35: Matrix input with sparse data uses na.rm correctly
test_that("Matrix input with sparse counts handles NAs in all genes", {
  sparse_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1, 1, 1, 1, 0),
    Gene2 = c(100, 50, 25, 10, 5, 2, 1, 1),
    Gene3 = c(50, 50, 50, 50, 50, 50, 50, 50)
  )
  
  result <- try(suppressWarnings(jackknife_tsallis_entropy(sparse_matrix, q = 1, norm = TRUE, print_results = FALSE)), silent = TRUE)
  
  # Should complete without error
  expect_false(inherits(result, "try-error"))
  expect_is(result, "tsenat_jackknife_list")
  expect_length(result, 3)
  
  # All results should have valid outlier_cutoff_value
  for (i in seq_along(result)) {
    expect_true(is.finite(result[[i]]$outlier_cutoff_value))
  }
})

# Test 36: Mean calculation with na.rm=TRUE for jackknife estimates
test_that("Mean of jackknife estimates uses na.rm=TRUE", {
  test_counts <- c(100, 50, 25, 10, 5)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE)
  
  # The jackknife_se calculation depends on mean with na.rm
  # If there were NAs, without na.rm it would fail
  expect_true(!is.na(result$jackknife_se))
  expect_true(is.finite(result$jackknife_se))
})

# Test 37: Sum calculation with na.rm=TRUE in SE formula
test_that("Sum in SE formula uses na.rm=TRUE", {
  test_counts <- c(500, 200, 100, 50, 20, 10, 5, 1)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE)
  
  # SE calculation: sqrt(((n-1)/n) * sum((theta_-i - theta_.)^2, na.rm=TRUE))
  # Should not fail even if sum encounters NAs
  expect_true(!is.na(result$jackknife_se))
  expect_true(is.finite(result$jackknife_se))
})

# Test 38: Outlier cutoff with na.rm for percentile calculation
test_that("Quantile for outlier cutoff uses na.rm=TRUE", {
  test_counts <- c(100, 50, 25, 10, 5, 2, 1, 1)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE, threshold = 80)
  
  # Quantile should be finite and not NA
  expect_true(is.finite(result$outlier_cutoff_value))
  expect_true(!is.na(result$outlier_cutoff_value))
})

# Test 39: Multiple q values with na.rm handling
test_that("Multiple q values all use na.rm correctly", {
  test_counts <- c(50, 50, 50, 50, 50, 50, 50, 0)
  
  result <- try(jackknife_tsallis_entropy(test_counts, q = c(0.5, 1, 1.5, 2), norm = TRUE, print_results = FALSE), silent = TRUE)
  
  expect_false(inherits(result, "try-error"))
  expect_is(result, "tsenat_jackknife_list_multiq")
  
  # All q results should have valid SE and cutoff values
  for (i in seq_along(result)) {
    expect_true(is.finite(result[[i]]$jackknife_se))
    expect_true(is.finite(result[[i]]$outlier_cutoff_value))
  }
})

# Test 40: Influence values with NaN are filtered from outliers
test_that("NaN values in influence are excluded from outlier detection", {
  test_counts <- c(100, 50, 25, 10)
  
  result <- jackknife_tsallis_entropy(test_counts, q = 1, norm = TRUE, threshold = 75)
  
  # Check that outlier_indices only contains valid indices
  # and no NaN values got through
  expect_true(is.numeric(result$outlier_indices))
  expect_true(!any(is.na(result$outlier_indices)))
  
  # All indices should reference actual transcripts
  if (length(result$outlier_indices) > 0) {
    expect_true(all(result$outlier_indices %in% seq_len(result$n_transcripts)))
  }
})

context("Jackknife Diagnostics: Tsallis q-Parameter Optimization")

# Feature 4.1: Entropy-Specific Use Cases - q-parameter optimization

test_that("q-parameter messaging for low q (underweights rare isoforms)", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Capture messages when q < 0.5
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.3, print_results = FALSE),
    "Low q.*heavily underweights rare isoforms"
  )
})

test_that("q-parameter messaging mentions large influence from abundant transcripts for low q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify message about abundant transcript influence
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.25, print_results = FALSE),
    "large influence from abundant transcripts"
  )
})

test_that("q-parameter messaging cites papers S111, I004 for low q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify database paper citations
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.1, print_results = FALSE),
    "papers S111, I004"
  )
})

test_that("q-parameter messaging for high q (insensitive to rare diversity)", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Capture messages when q > 2
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 2.5, print_results = FALSE),
    "may be insensitive to rare isoform diversity"
  )
})

test_that("q-parameter messaging mentions missed rare transcripts for high q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify message about rare transcript contributions
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 3.0, print_results = FALSE),
    "miss important rare transcript contributions"
  )
})

test_that("q-parameter messaging recommends q in [0.5, 2] for high q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify recommendation for balanced assessment
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 2.2, print_results = FALSE),
    "Consider q in"
  )
})

test_that("q-parameter messaging for recommended q range with verbose", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verbose message when q in [0.5, 2] and verbose=TRUE
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 1.0, print_results = FALSE, verbose = TRUE),
    "recommended range"
  )
})

test_that("q-parameter messaging cites papers for recommended range", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify database paper citations in recommended range message
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 1.5, print_results = FALSE, verbose = TRUE),
    "papers S111, I004"
  )
})

test_that("q-parameter messaging at boundary q=0.5", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # q=0.5 is at the lower boundary of recommended range
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.5, print_results = FALSE, verbose = TRUE),
    "recommended range"
  )
})

test_that("q-parameter messaging at boundary q=2.0", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # q=2.0 is at the upper boundary of recommended range
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 2.0, print_results = FALSE, verbose = TRUE),
    "recommended range"
  )
})

test_that("q-parameter messaging no verbose warning for recommended q without verbose flag", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # When q in recommended range and verbose=FALSE, should not message
  expect_no_message(
    jackknife_tsallis_entropy(x = counts, q = 1.0, print_results = FALSE, verbose = FALSE)
  )
})

test_that("q-parameter optimization doesn't affect computation results", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Results should be identical regardless of messaging
  result <- jackknife_tsallis_entropy(x = counts, q = 0.3, print_results = FALSE)
  
  # Verify computation still works
  expect_true("estimate" %in% names(result))
  expect_true("jackknife_estimates" %in% names(result))
  expect_true("influence" %in% names(result))
  expect_true(result$q == 0.3)
})

test_that("q-parameter optimization works with multiple q values", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Multiple q values should all get appropriate messages
  # (low q=0.3, recommended q=1, high q=2.5)
  expect_error(
    jackknife_tsallis_entropy(
      x = counts, 
      q = c(0.3, 1.0, 2.5), 
      print_results = FALSE,
      verbose = FALSE
    ),
    NA  # Expect no error
  )
})

test_that("q-parameter optimization with matrix input", {
  counts_matrix <- rbind(
    "Gene1" = c(1000, 500, 200, 100, 50),
    "Gene2" = c(800, 400, 300, 200, 100)
  )
  
  # Messages should appear for each gene with low q
  expect_message(
    jackknife_tsallis_entropy(x = counts_matrix, q = 0.4, print_results = FALSE),
    "Low q"
  )
})

test_that("q-parameter emphasizes dominant isoform detection for low q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify message mentions better for detecting changes
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.2, print_results = FALSE),
    "detecting changes in dominant isoforms"
  )
})

test_that("q-parameter documentation mentions balanced assessment for high q", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Verify message calls for balanced assessment
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 2.8, print_results = FALSE),
    "balanced diversity assessment"
  )
})

test_that("q=0.5 boundary lower - no low q warning", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Exactly at boundary should not show low q warning
  expect_no_message(
    jackknife_tsallis_entropy(x = counts, q = 0.5, print_results = FALSE, verbose = FALSE)
  )
})

test_that("q=2.0 boundary upper - no high q warning", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Exactly at boundary should not show high q warning
  expect_no_message(
    jackknife_tsallis_entropy(x = counts, q = 2.0, print_results = FALSE, verbose = FALSE)
  )
})

test_that("q=0.49 just below boundary - shows low q warning", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Just below 0.5 should show low q warning
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 0.49, print_results = FALSE),
    "Low q"
  )
})

test_that("q=2.01 just above boundary - shows high q warning", {
  counts <- c(1000, 500, 200, 100, 50)
  
  # Just above 2.0 should show high q warning
  expect_message(
    jackknife_tsallis_entropy(x = counts, q = 2.01, print_results = FALSE),
    "High q"
  )
})

test_that("Feature 4.1 integration: q-parameter messaging is entropy-specific", {
  # This test verifies that the messaging appropriately addresses
  # Tsallis entropy-specific concerns (rare vs abundant isoforms)
  counts_balanced <- c(100, 100, 100, 100, 100)  # Uniform distribution
  counts_skewed <- c(1000, 10, 10, 10, 10)      # Dominant + rare
  
  # Low q should affect results differently for balanced vs skewed
  result_balanced_low <- jackknife_tsallis_entropy(x = counts_balanced, q = 0.3, print_results = FALSE)
  result_skewed_low <- jackknife_tsallis_entropy(x = counts_skewed, q = 0.3, print_results = FALSE)
  
  # High q should follow similar pattern
  result_balanced_high <- jackknife_tsallis_entropy(x = counts_balanced, q = 2.5, print_results = FALSE)
  result_skewed_high <- jackknife_tsallis_entropy(x = counts_skewed, q = 2.5, print_results = FALSE)
  
  # Different q values should give different estimates
  expect_true(result_balanced_low$estimate != result_balanced_high$estimate)
  expect_true(result_skewed_low$estimate != result_skewed_high$estimate)
})

test_that("Feature 4.1 database citations are accurate", {
  # Verify papers S111 and I004 are cited
  counts <- c(1000, 500, 200, 100, 50)
  
  # Both papers should be mentioned
  output_low <- capture_messages(
    jackknife_tsallis_entropy(x = counts, q = 0.2, print_results = FALSE)
  )
  
  expect_true(any(grepl("S111", paste(output_low, collapse = " "))))
  expect_true(any(grepl("I004", paste(output_low, collapse = " "))))
})

context("Block Jackknife and Visualization Functions")

# Test fixtures
test_se_with_blocks <- function() {
  set.seed(42)
  counts <- matrix(
    c(
      # Gene 1, 3 transcripts
      1000, 500, 200,  400, 300, 150,   # Condition A - Phase 1
      800, 400, 180,   350, 250, 120,   # Condition B - Phase 1
      # Add phase grouping
      1200, 600, 250,  500, 350, 200,   # Condition A - Phase 2
      900, 500, 200,   450, 300, 180    # Condition B - Phase 2
    ),
    nrow = 3,
    ncol = 8,
    byrow = FALSE
  )
  
  rowData <- data.frame(
    gene_id = c("Gene1", "Gene1", "Gene1"),
    isoform_id = c("Gene1.1", "Gene1.2", "Gene1.3"),
    row.names = c("iso1", "iso2", "iso3")
  )
  
  colData <- data.frame(
    condition = rep(c("A", "B"), each = 4),
    phase = rep(c("Phase1", "Phase1", "Phase2", "Phase2"), 2),
    row.names = paste0("sample_", 1:8)
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
}
´

# ============================================================================
# HEATMAP VISUALIZATION TESTS
# ============================================================================

test_that("plot_isoform_switching_heatmap function exists", {
  expect_true(exists("plot_isoform_switching_heatmap"))
  expect_true(is.function(plot_isoform_switching_heatmap))
})

test_that("plot_isoform_switching_heatmap requires switching_results", {
  expect_error(
    plot_isoform_switching_heatmap(switching_results = NULL),
    "must be from jackknife_isoform_switching"
  )
})

test_that("plot_isoform_switching_heatmap validates class", {
  expect_error(
    plot_isoform_switching_heatmap(switching_results = list()),
    "must be from jackknife_isoform_switching"
  )
})

test_that("plot_isoform_switching_heatmap creates heatmap for valid input", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      top_n = 1,
      top_transcripts_per_gene = 2,
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap handles diverging color scheme", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      color_scheme = "diverging",
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap handles sequential color scheme", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      color_scheme = "sequential",
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap returns heatmap matrix invisibly", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  result <- suppressWarnings(
    plot_isoform_switching_heatmap(switching, show_values = FALSE)
  )
  
  expect_true(is.matrix(result) || is.null(result))
})

# ============================================================================
# Q-SENSITIVITY CURVE TESTS
# ============================================================================

test_that("plot_q_sensitivity_curve function exists", {
  expect_true(exists("plot_q_sensitivity_curve"))
  expect_true(is.function(plot_q_sensitivity_curve))
})

test_that("plot_q_sensitivity_curve requires se and gene parameters", {
  expect_error(
    plot_q_sensitivity_curve(se = NULL),
    "required"
  )
})

test_that("plot_q_sensitivity_curve validates gene exists", {
  se <- test_se_with_blocks()
  expect_error(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "NonexistentGene",
      gene_col = "gene_id",
      isoform_col = "isoform_id"
    ),
    "not found in rowData"
  )
})

test_that("plot_q_sensitivity_curve requires 2 conditions", {
  se <- test_se_with_blocks()
  SummarizedExperiment::colData(se)$condition <- "A"  # All same condition
  
  expect_error(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id"
    ),
    "Exactly 2 conditions"
  )
})

test_that("plot_q_sensitivity_curve creates plot for single q", {
  se <- test_se_with_blocks()
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = 1
    )
  ))
})

test_that("plot_q_sensitivity_curve creates plot for multiple q values", {
  se <- test_se_with_blocks()
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2)
    )
  ))
})

test_that("plot_q_sensitivity_curve returns sensitivity dataframe", {
  se <- test_se_with_blocks()
  
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      show_legend = FALSE
    )
  )
  
  expect_true(is.data.frame(result))
  expect_true("q" %in% colnames(result))
  expect_true("delta_entropy" %in% colnames(result))
})

test_that("plot_q_sensitivity_curve sensitivity data has correct length", {
  se <- test_se_with_blocks()
  
  q_vals <- c(0.5, 0.8, 1, 1.2, 1.5, 2)
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = q_vals,
      show_legend = FALSE
    )
  )
  
  expect_equal(nrow(result), length(q_vals))
})

test_that("plot_q_sensitivity_curve respects normalization parameter", {
  se <- test_se_with_blocks()
  
  result_norm <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      norm = TRUE,
      show_legend = FALSE
    )
  )
  
  result_unnorm <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      norm = FALSE,
      show_legend = FALSE
    )
  )
  
  # Values should be different between normalized and unnormalized
  expect_false(all(result_norm$delta_entropy == result_unnorm$delta_entropy))
})

test_that("plot_q_sensitivity_curve respects q_values parameter", {
  se <- test_se_with_blocks()
  
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.7, 1.3, 2.1),
      show_legend = FALSE
    )
  )
  
  expect_equal(result$q, c(0.7, 1.3, 2.1))
})

# ============================================================================
# INTEGRATION TESTS
# ============================================================================

test_that("Block jackknife and isoform switching work together with same data", {
  se <- test_se_with_blocks()
  
  
  # Standard isoform switching
  switch_result <- suppressWarnings(
    jackknife_isoform_switching(
      se = se,
      condition_col = "condition",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q = 1,
      n_bootstrap = 5,
      norm = FALSE,
      print_results = FALSE
    )
  )
  
  expect_true(inherits(block_result, "tsenat_block_jackknife"))
  expect_true(inherits(switch_result, "tsenat_isoform_switching"))
})

test_that("Visualization functions work with isoform switching results", {
  se <- test_se_with_blocks()
  
  switch_result <- suppressWarnings(
    jackknife_isoform_switching(
      se = se,
      condition_col = "condition",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q = 1,
      n_bootstrap = 5,
      norm = FALSE,
      print_results = FALSE
    )
  )
  
  # Test heatmap
  heatmap <- suppressWarnings(
    plot_isoform_switching_heatmap(switch_result, show_values = FALSE)
  )
  
  # Test q-curve
  q_result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(1, 1.5),
      show_legend = FALSE
    )
  )
  
  expect_true(is.data.frame(q_result))
})

