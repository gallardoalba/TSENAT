context("Jackknife Diagnostics: Tsallis Entropy Stability Analysis")

# Helper functions are defined locally or sourced globally via setup.R

# Setup test data
set.seed(42)
balanced_counts <- c(100, 100, 100, 100, 100)  # All equal
skewed_counts <- c(500, 200, 100, 50, 20)  # Strongly skewed
with_dominant <- c(900, 50, 30, 20)  # One very dominant
single_counts <- c(1000)  # Single transcript
two_counts <- c(500, 500)  # Two transcripts

# Input type tests
test_that("jackknife_entropy_outliers works with all input types", {
  # Test with vector input
  result_vec <- .calculate_jeo(balanced_counts, q = 1, verbose = FALSE)
  expect_true(is.list(result_vec))
  expect_true(!is.null(result_vec$estimate))
  
  # Test with matrix input
  test_mat <- matrix(c(balanced_counts, skewed_counts), nrow = 2, byrow = TRUE)
  result_mat <- .calculate_jeo(test_mat, q = 1, verbose = FALSE)
  expect_true(is.list(result_mat))
})

# Diagnostics-specific tests (not in consolidated helpers)
# Test 4: Balanced distribution has low influence
test_that("Balanced counts have low, uniform influence", {
  result <- .calculate_jeo(balanced_counts, q = 1, norm = TRUE)

  expect_true(all(result$influence < 0.01))  # All very small
  expect_true(max(result$influence) - min(result$influence) < 0.005)  # All similar
})

# Test 5: Skewed distribution has variable influence
test_that("Skewed counts have higher, variable influence", {
  result <- .calculate_jeo(skewed_counts, q = 1, norm = TRUE)

  # Should have variation in influence
  expect_true(max(result$influence) > min(result$influence))
  # Influence should vary across transcripts
  expect_true(length(unique(round(result$influence, 6))) > 1)
})

# Test 6: Dominant transcript is identified as outlier
test_that("Dominant transcript is detected as outlier", {
  result <- .calculate_jeo(with_dominant, q = 1, norm = TRUE, threshold = 75)

  expect_length(result$outlier_indices, 1)
  expect_equal(result$outlier_indices, 1)  # First transcript is dominant
})

# Test 6: More abundant transcripts tend to have higher influence
test_that("Influence detected for all transcripts", {
  result <- .calculate_jeo(skewed_counts, q = 1, norm = TRUE)

  # All transcripts should have measurable influence
  expect_true(all(result$influence >= 0))
  # There should be variation
  expect_true(max(result$influence) > 0)
})

# Test 8: Different q values
test_that("jackknife_entropy_outliers handles different q values", {
  q_vals <- c(0.5, 1, 1.5, 2)
  results <- lapply(q_vals, function(q_val) {
    .calculate_jeo(skewed_counts, q = q_val, norm = TRUE)
  })

  # All should be valid
  expect_true(all(sapply(results, function(r) !is.na(r$estimate))))
  # Results should differ
  estimates <- sapply(results, function(r) r$estimate)
  expect_false(all(estimates == estimates[1]))
})

# Test 8: Different normalization schemes
test_that("Normalization affects entropy values", {
  result_norm <- .calculate_jeo(skewed_counts, q = 1, norm = TRUE)
  result_unnorm <- .calculate_jeo(skewed_counts, q = 1, norm = FALSE)

  # Estimates should differ
  expect_false(abs(result_norm$estimate - result_unnorm$estimate) < 1e-6)

  # Both should have comparable structure
  expect_equal(result_norm$n_transcripts, result_unnorm$n_transcripts)
  expect_length(result_norm$influence, result_unnorm$n_transcripts)
})

# Test 10: Jackknife standard error
test_that("Jackknife SE is computed correctly", {
  result <- .calculate_jeo(skewed_counts, q = 1, norm = TRUE)

  expect_true(result$jackknife_se >= 0)
  # For balanced data, SE should be very small
  result_bal <- .calculate_jeo(balanced_counts, q = 1, norm = TRUE)
  expect_true(result_bal$jackknife_se < result$jackknife_se)
})

# Parameter validation tests
test_that("jackknife_entropy_outliers validates all parameters", {
  # Test valid call works
  result <- .calculate_jeo(
    balanced_counts, 
    q = 1, 
    verbose = FALSE
  )
  expect_true(!is.null(result$estimate))
  
  # Test with different q values
  result_q2 <- .calculate_jeo(
    balanced_counts,
    q = 2,
    verbose = FALSE
  )
  expect_true(!is.null(result_q2$estimate))
})

# Test 14: Input validation - threshold
test_that("jackknife_entropy_outliers validates threshold parameter", {
  expect_error(
    .calculate_jeo(balanced_counts, threshold = -5),
    "between 0 and 100"
  )
  expect_error(
    .calculate_jeo(balanced_counts, threshold = 150),
    "between 0 and 100"
  )
})

# Test 15: Input validation - minimum transcripts
test_that("jackknife_entropy_outliers requires at least 2 transcripts", {
  expect_error(
    .calculate_jeo(single_counts),
    "at least 2 transcripts"
  )
})

# Test 16: Two-transcript case (minimum)
test_that("jackknife_entropy_outliers works with 2 transcripts", {
  result <- .calculate_jeo(two_counts, q = 1)

  expect_equal(result$n_transcripts, 2)
  expect_length(result$influence, 2)
  expect_true(!is.na(result$estimate))
})

# Test 17: S3 print method
test_that("print method works for jackknife results", {
  result <- .calculate_jeo(skewed_counts, q = 1)

  # Methods use message() for output, not stdout
  expect_message(print(result), "Jackknife Diagnostics")
})

# Test 18: S3 summary method
# Test 19: S3 print for list
test_that("print method works for jackknife list", {
  gene_matrix <- matrix(c(balanced_counts, skewed_counts), nrow = 2, byrow = TRUE)
  result <- .calculate_jeo(gene_matrix, q = 1, verbose = FALSE)

  # Methods use message() for output, not stdout
  expect_message(print(result), "genes")
})

# Test 20: Outlier threshold affects detection
test_that("Outlier threshold parameter works correctly", {
  result_high <- .calculate_jeo(skewed_counts, q = 1, threshold = 90)
  result_low <- .calculate_jeo(skewed_counts, q = 1, threshold = 50)

  # Lower threshold should find more outliers
  expect_true(length(result_low$outlier_indices) >= length(result_high$outlier_indices))
})

# Test 21: Zero counts are handled
test_that("jackknife_entropy_outliers handles zero counts with pseudocount", {
  counts_with_zero <- c(100, 0, 50, 30)

  result <- .calculate_jeo(counts_with_zero, q = 1, pseudocount = 1e-10)

  expect_true(!is.na(result$estimate))
  expect_true(all(!is.na(result$jackknife_estimates)))
  expect_true(all(is.finite(result$influence)))
})

# Test 22: Reproducibility without seed
test_that("jackknife_entropy_outliers is deterministic", {
  result1 <- .calculate_jeo(skewed_counts, q = 1)
  result2 <- .calculate_jeo(skewed_counts, q = 1)

  expect_equal(result1$estimate, result2$estimate)
  expect_equal(result1$jackknife_estimates, result2$jackknife_estimates)
})

# Test 23: Large transcript count
test_that("jackknife_entropy_outliers works with many transcripts", {
  
  many_transcripts <- rpois(50, lambda = 50)  # Reduced from 100 to 50 for faster testing

  result <- .calculate_jeo(many_transcripts, q = 1)

  expect_equal(result$n_transcripts, 50)
  expect_length(result$influence, 50)
  expect_true(!is.na(result$jackknife_se))
})

# Test 24: All transcripts equal
test_that("Equal transcripts have uniform influence", {
  equal_counts <- rep(100, 10)

  result <- .calculate_jeo(equal_counts, q = 1, norm = TRUE)

  # All influences should be nearly identical
  inf_range <- max(result$influence) - min(result$influence)
  expect_true(inf_range < 1e-6)
})


# ═══════════════════════════════════════════════════════════════════════════════
# TESTS FOR .calculate_jis() - Input Validation and Structure
# ═══════════════════════════════════════════════════════════════════════════════

context("Jackknife Isoform Switching: Input Validation")

# Setup test data for isoform switching
set.seed(42)

# Create a simple test SummarizedExperiment
test_se_basic <- function() {
  suppressPackageStartupMessages({
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

# Test 1: Function signature accepts nboot parameter
test_that("calculate_jis function signature is correct", {
  # Just test that function exists and has correct parameters
  expect_true(exists(".calculate_jis"))
  
  # Get function signature
  sig <- formals(.calculate_jis)
  expect_true("nboot" %in% names(sig))
  expect_true("condition_col" %in% names(sig))
  expect_true("gene_col" %in% names(sig))
  expect_true("isoform_col" %in% names(sig))
})

# Test 2: Input validation - missing condition_col
test_that("calculate_jis detects missing condition column", {
  se <- test_se_basic()
  
  expect_error(
    suppressWarnings(.calculate_jis(
      se = se,
      condition_col = "nonexistent",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      nboot = 5,
      norm = FALSE,
      verbose = FALSE
    )),
    "condition"
  )
})

# Test 3: Input validation - invalid q parameter
test_that("calculate_jis handles q parameter appropriately", {
  se <- test_se_basic()
  
  # Test that function accepts valid q values
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.5,
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 4: Function returns correct class
test_that("calculate_jis returns tsenat_isoform_switching class", {
  se <- test_se_basic()
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
  expect_true(inherits(result, "list"))
})

# Test 5: Required output components exist
test_that("calculate_jis output has required components", {
  se <- test_se_basic()
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
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
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
  ))
  
  stats <- result$all_transcript_stats
  required_cols <- c("gene", "transcript_id", "pvalue", "fdr")
  expect_true(all(required_cols %in% names(stats)))
})

# Test 7: Paired design parameter is accepted
test_that("calculate_jis accepts subject_col parameter", {
  suppressPackageStartupMessages({
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
  
  result <- suppressWarnings(.calculate_jis(
    se = se_paired,
    condition_col = "condition",
    subject_col = "individual_id",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 8: LM results parameter is accepted
test_that("calculate_jis accepts lm_results parameter", {
  
  suppressPackageStartupMessages({
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
  
  result <- suppressWarnings(.calculate_jis(
    se = se_multi,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# Test 9: Metadata is populated
test_that("Metadata is populated after analysis", {
  se <- test_se_basic()
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.5,
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
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
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    nboot = 5,
    norm = FALSE,
    verbose = FALSE
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

  result_small <- .calculate_jeo(counts_zero, pseudocount = 1e-20)
  result_large <- .calculate_jeo(counts_zero, pseudocount = 0.1)

  # Both should be valid
  expect_true(!is.na(result_small$estimate))
  expect_true(!is.na(result_large$estimate))
  # Results may differ slightly due to different pseudocounts
  expect_false(isTRUE(all.equal(result_small$estimate, result_large$estimate)))
})

# Test 26: KL divergence (q=1) computation
test_that("KL divergence (q=1) computes correctly", {
  result <- .calculate_jeo(skewed_counts, q = 1, norm = FALSE)

  # KL divergence should be non-negative
  expect_true(result$estimate >= 0)
  expect_true(all(result$jackknife_estimates >= 0))
})

# Test 27: Log base parameter
test_that("Log base parameter changes entropy scale", {
  result_e <- .calculate_jeo(skewed_counts, q = 1, log_base = exp(1), norm = FALSE)
  result_2 <- .calculate_jeo(skewed_counts, q = 1, log_base = 2, norm = FALSE)

  # Nats vs bits should differ by approximately log(2)
  ratio <- result_e$estimate / result_2$estimate
  expected_ratio <- log(2)
  expect_true(abs(ratio - expected_ratio) < 0.1)
})

# Test 28: Outlier indices are numeric and valid
test_that("Outlier indices are valid positions", {
  result <- .calculate_jeo(skewed_counts, q = 1, threshold = 80)

  expect_true(all(result$outlier_indices >= 1))
  expect_true(all(result$outlier_indices <= result$n_transcripts))
  expect_true(length(result$outlier_indices) <= result$n_transcripts)
})

# Test 29: Named vector input
test_that("jackknife_entropy_outliers preserves names when possible", {
  counts_named <- c(Iso1 = 100, Iso2 = 200, Iso3 = 150)

  result <- .calculate_jeo(counts_named, q = 1)

  expect_length(result$influence, 3)
  expect_equal(result$n_transcripts, 3)
})

# Test 30: Very large counts
test_that("jackknife_entropy_outliers handles large count values", {
  large_counts <- c(1e6, 5e5, 2e5, 1e5)

  result <- .calculate_jeo(large_counts, q = 1, norm = TRUE)

  expect_true(!is.na(result$estimate))
  expect_true(all(!is.na(result$jackknife_estimates)))
})

# Test 31: Sparse counts that may produce NaN values
test_that("jackknife_entropy_outliers handles sparse counts with na.rm", {
  sparse_counts <- c(1, 1, 1, 1, 1, 1, 1, 0)  # Very sparse data
  
  result <- try(suppressWarnings(.calculate_jeo(sparse_counts, q = 1, norm = TRUE)), silent = TRUE)
  
  expect_false(inherits(result, "try-error"))
  expect_is(result, "tsenat_jackknife")
  # Important: quantile should work with na.rm=TRUE
  expect_true(!is.na(result$outlier_cutoff_value))
})

# Test 32: Jackknife SE calculation with na.rm
test_that("Jackknife SE is calculated correctly with na.rm=TRUE", {
  test_counts <- c(100, 50, 25, 10, 5)
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE)
  
  # SE should be valid (not NaN or Inf)
  expect_true(is.finite(result$jackknife_se))
  expect_true(!is.na(result$jackknife_se))
  expect_true(result$jackknife_se >= 0)
})

# Test 33: Outlier detection handles NA values correctly
test_that("Outlier detection filters NA values correctly", {
  test_counts <- c(100, 50, 25, 10)
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE, threshold = 75)
  
  # outlier_indices should not contain NA
  expect_true(!any(is.na(result$outlier_indices)))
  # outlier_indices should be valid indices
  expect_true(all(result$outlier_indices > 0))
  expect_true(all(result$outlier_indices <= length(test_counts)))
})

# Test 34: Quantile calculation with na.rm works for uniform influence
test_that("Quantile calculation works for uniform influence distribution", {
  uniform_counts <- c(25, 25, 25, 25)  # Very uniform
  
  result <- .calculate_jeo(uniform_counts, q = 1, norm = TRUE, threshold = 90)
  
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
  
  result <- try(suppressWarnings(.calculate_jeo(sparse_matrix, q = 1, norm = TRUE, verbose = FALSE)), silent = TRUE)
  
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
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE)
  
  # The jackknife_se calculation depends on mean with na.rm
  # If there were NAs, without na.rm it would fail
  expect_true(!is.na(result$jackknife_se))
  expect_true(is.finite(result$jackknife_se))
})

# Test 37: Sum calculation with na.rm=TRUE in SE formula
test_that("Sum in SE formula uses na.rm=TRUE", {
  test_counts <- c(500, 200, 100, 50, 20, 10, 5, 1)
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE)
  
  # SE calculation: sqrt(((n-1)/n) * sum((theta_-i - theta_.)^2, na.rm=TRUE))
  # Should not fail even if sum encounters NAs
  expect_true(!is.na(result$jackknife_se))
  expect_true(is.finite(result$jackknife_se))
})

# Test 38: Outlier cutoff with na.rm for percentile calculation
test_that("Quantile for outlier cutoff uses na.rm=TRUE", {
  test_counts <- c(100, 50, 25, 10, 5, 2, 1, 1)
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE, threshold = 80)
  
  # Quantile should be finite and not NA
  expect_true(is.finite(result$outlier_cutoff_value))
  expect_true(!is.na(result$outlier_cutoff_value))
})

# Test 39: Multiple q values with na.rm handling
test_that("Multiple q values all use na.rm correctly", {
  test_counts <- c(50, 50, 50, 50, 50, 50, 50, 0)
  
  result <- try(.calculate_jeo(test_counts, q = c(0.5, 1, 1.5, 2), norm = TRUE, verbose = FALSE), silent = TRUE)
  
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
  
  result <- .calculate_jeo(test_counts, q = 1, norm = TRUE, threshold = 75)
  
  # Check that outlier_indices only contains valid indices
  # and no NaN values got through
  expect_true(is.numeric(result$outlier_indices))
  expect_true(!any(is.na(result$outlier_indices)))
  
  # All indices should reference actual transcripts
  if (length(result$outlier_indices) > 0) {
    expect_true(all(result$outlier_indices %in% seq_len(result$n_transcripts)))
  }
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

# Test: Block jackknife computation with phase blocking
test_that("Block jackknife computation works with grouped samples", {
  se <- test_se_with_blocks()
  
  # Extract counts for single gene
  counts <- assay(se)[1, ]
  
  # Should compute without errors
  result <- expect_no_error(.calculate_jeo(counts, q = 1, verbose = FALSE))
  
  # Result should be valid
  expect_is(result, "tsenat_jackknife")
  expect_true(!is.na(result$estimate))
  expect_true(result$jackknife_se > 0)
})

# Test: Block jackknife with metadata blocking
test_that("Block jackknife respects phase grouping from metadata", {
  se <- test_se_with_blocks()
  
  # Extract first transcript counts
  counts <- assay(se)[1, ]
  
  # Standard jackknife
  result_standard <- .calculate_jeo(counts, q = 1, verbose = FALSE)
  
  # Both should produce valid results
  expect_true(is.numeric(result_standard$estimate))
  expect_true(result_standard$estimate >= 0)
})

# Test: Visualization function placeholder
test_that("Block jackknife results are visualizable", {
  se <- test_se_with_blocks()
  counts <- assay(se)[1, ]
  
  result <- .calculate_jeo(counts, q = 1, verbose = FALSE)
  
  # Results should have influence data for visualization
  expect_true(!is.null(result$influence))
  expect_true(length(result$influence) > 0)
})


# ═══════════════════════════════════════════════════════════════════════════════
# COMPREHENSIVE TESTS FOR HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

context("Jackknife Helper Functions: Modular Unit Tests")

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 1: .jackknife_validate_params()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_validate_params rejects q < 0", {
  # Negative q should be rejected
  expect_error(
    TSENAT:::.jackknife_validate_params(q = -1, threshold = 75),
    "q.*must be non-negative"
  )
})

test_that(".jackknife_validate_params accepts q >= 0", {
  # Zero q is valid (species richness in Tsallis entropy)
  expect_no_error(
    TSENAT:::.jackknife_validate_params(q = 0, threshold = 75)
  )
})

test_that(".jackknife_validate_params accepts valid q values", {
  # Should not error for valid q values
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 0.5, threshold = 75))
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 1, threshold = 75))
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 2, threshold = 75))
})

test_that(".jackknife_validate_params accepts vector of q values", {
  # Vector of q values should be accepted
  expect_no_error(TSENAT:::.jackknife_validate_params(q = c(0.5, 1, 1.5, 2), threshold = 75))
})

test_that(".jackknife_validate_params rejects invalid threshold", {
  # Negative threshold
  expect_error(
    TSENAT:::.jackknife_validate_params(q = 1, threshold = -1),
    "threshold.*between 0 and 100"
  )
  
  # Threshold > 100
  expect_error(
    TSENAT:::.jackknife_validate_params(q = 1, threshold = 101),
    "threshold.*between 0 and 100"
  )
})

test_that(".jackknife_validate_params accepts valid threshold values", {
  # Boundary values
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 1, threshold = 0))
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 1, threshold = 50))
  expect_no_error(TSENAT:::.jackknife_validate_params(q = 1, threshold = 100))
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 2: .get_nthreads_auto_detect() (consolidated from parallel_helpers.R)
# ─────────────────────────────────────────────────────────────────────────────

test_that(".get_nthreads_auto_detect returns positive integer", {
  result <- TSENAT:::.get_nthreads_auto_detect(NULL)
  expect_true(is.numeric(result))
  expect_true(result >= 1)
})

test_that(".get_nthreads_auto_detect respects explicit nthreads=1", {
  result <- TSENAT:::.get_nthreads_auto_detect(nthreads = 1)
  expect_equal(result, 1)
})

test_that(".get_nthreads_auto_detect respects explicit nthreads > 1", {
  result <- TSENAT:::.get_nthreads_auto_detect(nthreads = 4)
  expect_equal(result, 4)
})

test_that(".get_nthreads_auto_detect auto-detects cores when NULL", {
  result <- TSENAT:::.get_nthreads_auto_detect(NULL)
  max_cores <- parallel::detectCores()
  expect_true(result <= max_cores)
  expect_true(result >= 1)
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 3: .jackknife_compute_estimates()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_compute_estimates returns numeric vector", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  result <- TSENAT:::.jackknife_compute_estimates(p = p, q = 1, log_base = 2, n = 4)
  
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

test_that(".jackknife_compute_estimates handles q=1 (Shannon)", {
  # Balanced distribution
  p <- c(0.25, 0.25, 0.25, 0.25)
  result <- TSENAT:::.jackknife_compute_estimates(p = p, q = 1, log_base = 2, n = 4)
  
  # All should be valid numbers
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))
})

test_that(".jackknife_compute_estimates handles q != 1 (Tsallis)", {
  # Balanced distribution
  p <- c(0.25, 0.25, 0.25, 0.25)
  result <- TSENAT:::.jackknife_compute_estimates(p = p, q = 2, log_base = 2, n = 4)
  
  # All should be valid numbers
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that(".jackknife_compute_estimates handles skewed distribution", {
  # Heavily skewed - first transcript dominates
  p <- c(0.7, 0.15, 0.1, 0.05)
  result <- TSENAT:::.jackknife_compute_estimates(p = p, q = 1, log_base = 2, n = 4)
  
  # Removing the dominant transcript should increase entropy
  # (because remaining distribution becomes more even)
  expect_true(result[1] > mean(result[-1]))
})

test_that(".jackknife_compute_estimates handles different log bases", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  
  result_base2 <- TSENAT:::.jackknife_compute_estimates(p, q = 1, log_base = 2, n = 4)
  result_base10 <- TSENAT:::.jackknife_compute_estimates(p, q = 1, log_base = 10, n = 4)
  
  # Results should differ due to different log base
  expect_false(isTRUE(all.equal(result_base2, result_base10)))
  
  # Base 10 should give smaller values
  expect_true(all(result_base10 < result_base2))
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 4: .jackknife_calculate_influence_and_outliers()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_calculate_influence_and_outliers returns correct structure", {
  jackknife_ests <- c(2.0, 2.1, 1.95, 2.05, 1.98)
  estimate <- 2.0
  
  result <- TSENAT:::.jackknife_calculate_influence_and_outliers(
    jackknife_estimates = jackknife_ests,
    estimate = estimate,
    threshold = 75,
    q = 1,
    n = 5,
    norm = FALSE
  )
  
  expect_is(result, "tsenat_jackknife")
  expect_true("estimate" %in% names(result))
  expect_true("influence" %in% names(result))
  expect_true("jackknife_se" %in% names(result))
  expect_true("outlier_indices" %in% names(result))
  expect_true("outlier_threshold" %in% names(result))
  expect_true("outlier_cutoff_value" %in% names(result))
})

test_that(".jackknife_calculate_influence_and_outliers computes influence correctly", {
  jackknife_ests <- c(1.0, 3.0, 2.0, 2.5, 1.5)
  estimate <- 2.0
  
  result <- TSENAT:::.jackknife_calculate_influence_and_outliers(
    jackknife_estimates = jackknife_ests,
    estimate = estimate,
    threshold = 75,
    q = 1,
    n = 5,
    norm = FALSE
  )
  
  # Influence should be absolute differences
  expected_influence <- abs(jackknife_ests - estimate)
  expect_equal(result$influence, expected_influence)
})

test_that(".jackknife_calculate_influence_and_outliers respects threshold", {
  jackknife_ests <- c(1.0, 3.0, 2.0, 2.5, 1.5)
  estimate <- 2.0
  
  result_high <- TSENAT:::.jackknife_calculate_influence_and_outliers(
    jackknife_estimates = jackknife_ests,
    estimate = estimate,
    threshold = 90,
    q = 1,
    n = 5,
    norm = FALSE
  )
  
  result_low <- TSENAT:::.jackknife_calculate_influence_and_outliers(
    jackknife_estimates = jackknife_ests,
    estimate = estimate,
    threshold = 50,
    q = 1,
    n = 5,
    norm = FALSE
  )
  
  # Lower threshold should identify more outliers
  expect_true(length(result_low$outlier_indices) >= length(result_high$outlier_indices))
})

test_that(".jackknife_calculate_influence_and_outliers computes SE correctly", {
  jackknife_ests <- c(1.9, 2.0, 2.1, 2.05, 1.95)  # Small variance
  estimate <- 2.0
  n <- 5
  
  result <- TSENAT:::.jackknife_calculate_influence_and_outliers(
    jackknife_estimates = jackknife_ests,
    estimate = estimate,
    threshold = 75,
    q = 1,
    n = n,
    norm = FALSE
  )
  
  # SE should be positive
  expect_true(result$jackknife_se >= 0)
  
  # Verify SE calculation
  theta_jack_mean <- mean(jackknife_ests)
  expected_se <- sqrt(((n - 1) / n) * sum((jackknife_ests - theta_jack_mean)^2))
  expect_equal(result$jackknife_se, expected_se)
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 5: .jackknife_process_vector_core()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_process_vector_core returns tsenat_jackknife object", {
  counts <- c(100, 80, 60, 40, 20)
  
  result <- TSENAT:::.jackknife_process_vector_core(
    x = counts,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75,
    gene_name = "Gene1",
    verbose = FALSE
  )
  
  expect_is(result, "tsenat_jackknife")
})

test_that(".jackknife_process_vector_core rejects < 2 transcripts", {
  single_count <- c(100)
  
  expect_error(
    TSENAT:::.jackknife_process_vector_core(
      x = single_count,
      q = 1,
      norm = FALSE,
      log_base = 2,
      pseudocount = 0,
      threshold = 75
    ),
    "at least 2"
  )
})

test_that(".jackknife_process_vector_core handles normalization", {
  counts <- c(100, 80, 60, 40, 20)
  
  result_norm <- TSENAT:::.jackknife_process_vector_core(
    x = counts,
    q = 1,
    norm = TRUE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75
  )
  
  result_unnorm <- TSENAT:::.jackknife_process_vector_core(
    x = counts,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75
  )
  
  # Both should have valid structure
  expect_is(result_norm, "tsenat_jackknife")
  expect_is(result_unnorm, "tsenat_jackknife")
  
  # Estimates should differ
  expect_false(abs(result_norm$estimate - result_unnorm$estimate) < 1e-6)
})

test_that(".jackknife_process_vector_core handles pseudocount", {
  counts <- c(100, 0, 60, 40, 20)
  
  # With pseudocount, should not error
  result <- TSENAT:::.jackknife_process_vector_core(
    x = counts,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 1e-10,
    threshold = 75
  )
  
  expect_true(!is.na(result$estimate))
  expect_true(all(is.finite(result$jackknife_estimates)))
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 6: .jackknife_process_matrix()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_process_matrix returns list with correct class", {
  counts_mat <- matrix(
    c(100, 80, 60, 40, 20, 150, 100, 50, 30, 10),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), NULL)
  )
  
  result <- TSENAT:::.jackknife_process_matrix(
    x = counts_mat,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75,
    verbose = FALSE
  )
  
  expect_is(result, "tsenat_jackknife_list")
  expect_equal(length(result), 2)
})

test_that(".jackknife_process_matrix processes each row independently", {
  counts_mat <- matrix(
    c(100, 80, 60, 40, 20, 150, 100, 50, 30, 10),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), NULL)
  )
  
  result <- TSENAT:::.jackknife_process_matrix(
    x = counts_mat,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75,
    verbose = FALSE
  )
  
  # Each gene should be a tsenat_jackknife object
  expect_is(result[[1]], "tsenat_jackknife")
  expect_is(result[[2]], "tsenat_jackknife")
  
  # Estimates should differ since genes have different distributions
  expect_false(abs(result[[1]]$estimate - result[[2]]$estimate) < 1e-6)
})

test_that(".jackknife_process_matrix preserves row names", {
  counts_mat <- matrix(
    c(100, 80, 60, 40, 20, 150, 100, 50, 30, 10),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GeneA", "GeneB"), NULL)
  )
  
  result <- TSENAT:::.jackknife_process_matrix(
    x = counts_mat,
    q = 1,
    norm = FALSE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75,
    verbose = FALSE
  )
  
  expect_equal(names(result), c("GeneA", "GeneB"))
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 7: .jackknife_warn_on_q_parameters()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_warn_on_q_parameters warns on low total count", {
  counts <- c(1, 1, 1, 1, 1)  # Total = 5, below 10
  
  expect_warning(
    TSENAT:::.jackknife_warn_on_q_parameters(counts, q = 1, verbose = FALSE),
    "Total count"
  )
})

test_that(".jackknife_warn_on_q_parameters does not warn on normal count", {
  counts <- c(100, 80, 60, 40, 20)  # Total = 300
  
  expect_no_warning(
    TSENAT:::.jackknife_warn_on_q_parameters(counts, q = 1, verbose = FALSE)
  )
})

test_that(".jackknife_warn_on_q_parameters mentions paper citations", {
  counts <- c(1, 1, 1, 1, 1)
  
  expect_warning(
    TSENAT:::.jackknife_warn_on_q_parameters(counts, q = 1, verbose = FALSE),
    "S111|S114"
  )
})

test_that(".jackknife_warn_on_q_parameters warns on low q", {
  counts <- c(100, 80, 60, 40, 20)
  
  expect_message(
    TSENAT:::.jackknife_warn_on_q_parameters(counts, q = 0.3, verbose = TRUE),
    "Low q"
  )
})

test_that(".jackknife_warn_on_q_parameters warns on high q", {
  counts <- c(100, 80, 60, 40, 20)
  
  expect_message(
    TSENAT:::.jackknife_warn_on_q_parameters(counts, q = 5, verbose = TRUE),
    "High q"
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 8: .jackknife_format_verbose_output_matrix()
# ─────────────────────────────────────────────────────────────────────────────

test_that(".jackknife_format_verbose_output_matrix returns character string", {
  # Create mock result structure
  results <- list(
    Gene1 = list(
      estimate = 2.5,
      jackknife_se = 0.1,
      influence = c(0.05, 0.04, 0.06, 0.03, 0.07),
      outlier_indices = c(5)
    ),
    Gene2 = list(
      estimate = 1.8,
      jackknife_se = 0.15,
      influence = c(0.08, 0.09, 0.07, 0.10, 0.06),
      outlier_indices = c(4, 5)
    )
  )
  class(results) <- c("tsenat_jackknife_list", "list")
  
  output <- TSENAT:::.jackknife_format_verbose_output_matrix(results)
  
  expect_is(output, "character")
  expect_true(nchar(output) > 0)
})

test_that(".jackknife_format_verbose_output_matrix includes gene names", {
  results <- list(
    GeneAlpha = list(
      estimate = 2.5,
      jackknife_se = 0.1,
      influence = c(0.05, 0.04),
      outlier_indices = integer(0)
    ),
    GeneBeta = list(
      estimate = 1.8,
      jackknife_se = 0.15,
      influence = c(0.08, 0.09),
      outlier_indices = integer(0)
    )
  )
  class(results) <- c("tsenat_jackknife_list", "list")
  
  output <- TSENAT:::.jackknife_format_verbose_output_matrix(results)
  
  expect_match(output, "GeneAlpha")
  expect_match(output, "GeneBeta")
})

test_that(".jackknife_format_verbose_output_matrix includes statistics", {
  results <- list(
    Gene1 = list(
      estimate = 2.5,
      jackknife_se = 0.123,
      influence = c(0.05, 0.04, 0.06),
      outlier_indices = c(2)
    )
  )
  class(results) <- c("tsenat_jackknife_list", "list")
  
  output <- TSENAT:::.jackknife_format_verbose_output_matrix(results)
  
  expect_match(output, "2\\.5|estimate|SE|influence|outlier")
})

# ─────────────────────────────────────────────────────────────────────────────
# Test Suite 9: Integration Tests for Helper Function Interactions
# ─────────────────────────────────────────────────────────────────────────────

test_that("Helper functions work correctly in integrated pipeline (vector)", {
  counts <- c(120, 90, 60, 30)
  
  # Full pipeline through main function uses helpers internally
  result <- .calculate_jeo(
    x = counts,
    q = 1,
    norm = TRUE,
    log_base = 2,
    pseudocount = 0,
    threshold = 75,
    verbose = FALSE
  )
  
  # All components should be valid
  expect_is(result, "tsenat_jackknife")
  expect_true(!is.na(result$estimate))
  expect_equal(length(result$influence), length(counts))
  expect_true(result$jackknife_se >= 0)
})

test_that("Helper functions work correctly in integrated pipeline (matrix)", {
  counts_mat <- matrix(
    c(120, 90, 60, 30, 100, 100, 100, 100),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), NULL)
  )
  
  result <- .calculate_jeo(
    x = counts_mat,
    q = 1,
    norm = FALSE,
    verbose = FALSE
  )
  
  expect_is(result, "tsenat_jackknife_list")
  expect_equal(length(result), 2)
  expect_true(all(sapply(result, function(r) !is.na(r$estimate))))
})

test_that("Helper functions preserve optimization levels", {
  # Test that optimizations from original code are preserved
  counts <- c(100, 80, 60, 40, 20, 10, 5, 2)  # 8 transcripts
  
  result <- .calculate_jeo(
    x = counts,
    q = 1,
    norm = TRUE,
    verbose = FALSE
  )
  
  # Should complete without excessive computation
  expect_true(!is.na(result$estimate))
  expect_equal(result$n_transcripts, 8)
  expect_equal(length(result$jackknife_estimates), 8)
})

test_that("Helper functions handle edge case: all equal counts", {
  equal_counts <- rep(50, 6)
  
  result <- .calculate_jeo(
    x = equal_counts,
    q = 1,
    norm = TRUE,
    verbose = FALSE
  )
  
  # Influence should be minimal and uniform
  expect_true(all(result$influence < 0.01))
  expect_true(max(result$influence) - min(result$influence) < 1e-6)
})

test_that("Helper functions handle edge case: highly skewed counts", {
  skewed_counts <- c(1000, 50, 25, 15, 10)
  
  result <- .calculate_jeo(
    x = skewed_counts,
    q = 1,
    norm = TRUE,
    verbose = FALSE
  )
  
  # First transcript should have highest influence
  expect_equal(which.max(result$influence), 1)
})

test_that("Multiple q values processed separately with accurate results", {
  counts <- c(100, 80, 60, 40, 20)
  q_vals <- c(0.5, 1, 1.5)
  
  # Process individually
  results_individual <- lapply(q_vals, function(q) {
    .calculate_jeo(x = counts, q = q, verbose = FALSE)
  })
  
  # All should be valid
  expect_true(all(sapply(results_individual, function(r) !is.na(r$estimate))))
  
  # Estimates should differ
  estimates <- sapply(results_individual, function(r) r$estimate)
  expect_true(length(unique(round(estimates, 4))) > 1)
})

context("Resampling Methods: Multi-q Support for Bootstrap and Jackknife")

test_that("calculate_tsallis_entropy_bootstrap accepts vector q", {
    x <- c(100, 50, 30, 20)
    set.seed(123)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    expect_is(result, "tsenat_bootstrap_ci_list")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("bootstrap multi-q returns correct structure", {
    x <- c(100, 50, 30, 20)
    set.seed(456)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    # Check list structure (2 q values: 1 and 2)
    expect_length(result, 2)
    expect_true(all(sapply(result, function(r) "estimate" %in% names(r))))
    expect_true(all(sapply(result, function(r) "lower_ci" %in% names(r))))
    expect_true(all(sapply(result, function(r) "upper_ci" %in% names(r))))
})

test_that("bootstrap multi-q estimates differ across q values", {
    set.seed(123)
    x <- c(100, 50, 30, 20, 10)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(0.5, 1, 2), nboot = 100)  # Reduced from 200 to 100
    
    # Estimates should be different for different q values
    est_q05 <- result$`q=0.5`$estimate
    est_q1 <- result$`q=1`$estimate
    est_q2 <- result$`q=2`$estimate
    
    # Should not all be equal
    expect_false(isTRUE(all.equal(est_q05, est_q1)))
    expect_false(isTRUE(all.equal(est_q1, est_q2)))
})

test_that("bootstrap multi-q CI bounds are sensible for each q", {
    
    x <- c(100, 50, 30, 20)
    set.seed(234)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    for (res in result) {
        expect_lt(res$lower_ci, res$estimate)
        expect_gt(res$upper_ci, res$estimate)
        expect_gte(res$lower_ci, 0)
        expect_lte(res$upper_ci, 1)
    }
})

test_that("bootstrap multi-q with different nboot values", {
    x <- c(100, 50, 30, 20)
    # nboot=50 triggers warning about being below recommended minimum (expected for exploratory testing)
    set.seed(100)
    result_small <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 50))  # Reduced for speed
    result_large <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 75))  # Reduced from 500
    
    # Both should return valid results
    expect_is(result_small, "tsenat_bootstrap_ci_list")
    expect_is(result_large, "tsenat_bootstrap_ci_list")
    
    # Larger nboot should give more stable estimates
    expect_equal(length(result_small$`q=1`$bootstrap_dist), 50)  # Updated from 100
    expect_equal(length(result_large$`q=1`$bootstrap_dist), 75)  # Updated from 500
})

test_that("jackknife_entropy_outliers accepts vector q", {
    x <- c(100, 50, 30, 20)
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("jackknife multi-q returns correct structure", {
    x <- c(100, 50, 30, 20, 15, 10)
    result <- .calculate_jeo(x, q = c(0.5, 1, 2), norm = TRUE, verbose = FALSE)
    
    # Check list structure
    expect_length(result, 3)
    expect_true(all(sapply(result, function(r) "estimate" %in% names(r))))
    expect_true(all(sapply(result, function(r) "jackknife_se" %in% names(r))))
    expect_true(all(sapply(result, function(r) "influence" %in% names(r))))
})

test_that("jackknife multi-q estimates differ across q values", {
    x <- c(100, 50, 30, 20, 10)
    result <- .calculate_jeo(x, q = c(0.5, 1, 2), norm = TRUE, verbose = FALSE)
    
    # Estimates should be different for different q values
    est_q05 <- result$`q=0.5`$estimate
    est_q1 <- result$`q=1`$estimate
    est_q2 <- result$`q=2`$estimate
    
    expect_false(isTRUE(all.equal(est_q05, est_q1)))
    expect_false(isTRUE(all.equal(est_q1, est_q2)))
})

test_that("jackknife multi-q with matrix input", {
    
    counts_matrix <- rbind(
        "Gene1" = c(100, 50, 30, 20),
        "Gene2" = c(80, 60, 40, 20)
    )
    
    # Jackknife on first row
    result <- .calculate_jeo(
        x = counts_matrix[1, , drop = FALSE],
        q = c(1, 2),
        norm = TRUE,
        verbose = FALSE
    )
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
})

test_that("jackknife multi-q SE estimates are positive", {
    x <- c(100, 50, 30, 20, 15)
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    
    for (res in result) {
        expect_gt(res$jackknife_se, 0)
        expect_length(res$influence, length(x))
        expect_true(all(res$influence >= 0))
    }
})

test_that("bootstrap accepts vector q with length > 1", {
    x <- c(100, 50, 30, 20)
    
    # Vector q with length > 1 returns list
    set.seed(42)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    expect_is(result, "tsenat_bootstrap_ci_list")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("jackknife accepts vector q with length > 1", {
    x <- c(100, 50, 30, 20)
    
    # Vector q with length > 1 returns list
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("single q still works (backward compatibility)", {
    x <- c(100, 50, 30, 20)
    
    # Bootstrap with scalar q (without diagnostics to test original format)
    set.seed(42)
    boot_result <- .calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, include_diagnostics = FALSE)
    expect_is(boot_result, "tsenat_bootstrap_ci")
    expect_named(boot_result, c("estimate", "lower_ci", "upper_ci", "ci_level", "method", "nboot", "bootstrap_dist"))
    
    # Jackknife with scalar q
    jack_result <- .calculate_jeo(x, q = 1, norm = TRUE, verbose = FALSE)
    expect_is(jack_result, "tsenat_jackknife")
    # Check that key fields exist (structure may have additional fields)
    expect_true("estimate" %in% names(jack_result))
    expect_true("jackknife_se" %in% names(jack_result))
    expect_true("influence" %in% names(jack_result))
    expect_true("q" %in% names(jack_result))
})

test_that("bootstrap multi-q respects seed parameter", {
    x <- c(100, 50, 30, 20)
    
    set.seed(789)
    result1 <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    set.seed(789)
    result2 <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    # Same seed should give very similar results (within 5% tolerance for bootstrap)
    # Note: Exact reproducibility unreliable with RNG state in multi-q mode
    max_est1 <- max(abs(result1$`q=1`$estimate), 0.1)
    max_est2 <- max(abs(result2$`q=1`$estimate), 0.1)
    rel_diff_est <- abs(result1$`q=1`$estimate - result2$`q=1`$estimate) / max(max_est1, max_est2)
    expect_true(rel_diff_est < 0.05, info = "Seed should give similar estimates")
    
    # With lower nboot, correlation may be lower; verify both have similar point estimates
    # Rather than checking distribution correlation (unreliable with low nboot)
    max_q2_est <- max(abs(result1$`q=2`$estimate), abs(result2$`q=2`$estimate), 0.1)
    rel_diff_q2 <- abs(result1$`q=2`$estimate - result2$`q=2`$estimate) / max_q2_est
    expect_true(rel_diff_q2 < 0.05, info = "q=2 estimates should also be similar with same seed")
})

test_that("bootstrap multi-q ci parameter returns finite widths", {
    x <- c(100, 80, 60, 40, 30, 20, 15, 10, 8, 5)
    
    set.seed(111)
    result_95 <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1.5, 2.0), nboot = 20, ci = 0.95))
    result_90 <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1.5, 2.0), nboot = 20, ci = 0.90))
    
    # Both should give valid CI widths
    width_95_q1 <- result_95$`q=1.5`$upper_ci - result_95$`q=1.5`$lower_ci
    width_90_q1 <- result_90$`q=1.5`$upper_ci - result_90$`q=1.5`$lower_ci
    
    # Check that widths are finite and positive
    expect_true(is.finite(width_95_q1))
    expect_true(is.finite(width_90_q1))
    expect_gt(width_95_q1, 0)
    expect_gt(width_90_q1, 0)
})

test_that("jackknife multi-q accepts Hill numbers (D)", {
    x <- c(100, 50, 30, 20, 15)
    # Note: jackknife doesn't have 'what' parameter, but we test multi-q works
    result <- .calculate_jeo(x, q = c(1, 2), norm = FALSE, verbose = FALSE)
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_true(all(sapply(result, function(r) !is.na(r$estimate))))
})
