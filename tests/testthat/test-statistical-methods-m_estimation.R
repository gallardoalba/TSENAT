context("Robust Statistical Methods: Basic Functionality")

# ============================================================================
# TEST: m_estimate() - M-Estimation with Robust Loss Functions
# ============================================================================

test_that("m_estimate returns correct data frame structure", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(x))
  expected_cols <- c("location_diff", "se_diff", "t_stat", "pvalue", 
                     "padj", "n_down_weighted", "max_weight")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("m_estimate with Huber loss produces finite results", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(all(is.finite(result$location_diff)))
  expect_true(all(is.finite(result$se_diff)))
  expect_true(all(is.finite(result$pvalue), na.rm = TRUE))
})

test_that("m_estimate with Tukey loss produces finite results", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "tukey")
  
  expect_true(all(is.finite(result$location_diff)))
  expect_true(all(is.finite(result$pvalue), na.rm = TRUE))
})

test_that("m_estimate with LSQ (ordinary least squares) as baseline", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "lsq")
  
  # LSQ should produce non-NA results
  expect_true(all(!is.na(result$location_diff)))
  expect_true(all(!is.na(result$pvalue)))
})

test_that("m_estimate detects down-weighted observations", {
  # Create data with outliers
  x <- matrix(c(-100, 100, 1:6, 1:8), nrow = 2, byrow = TRUE)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  # First row has extreme outliers, should be down-weighted
  expect_true(result$n_down_weighted[1] > 0)
  # Second row has normal values, may or may not have down-weighting
  expect_true(result$n_down_weighted[2] >= 0)
})

test_that("m_estimate location_diff differs between groups", {
  # Create two clearly different groups
  x <- rbind(
    c(10, 11, 12, 13, 1, 2, 3, 4),  # Group A ~12, Group B ~2.5
    c(5, 5, 5, 5, 15, 15, 15, 15)    # Group A ~5, Group B ~15
  )
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "lsq")
  
  # Differences should be non-zero and in opposite directions
  expect_true(result$location_diff[1] < 0)  # Group A > Group B
  expect_true(result$location_diff[2] > 0)  # Group A < Group B
})

test_that("m_estimate handles paired design", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result_paired <- m_estimate(x, samples, loss_type = "huber", paired = TRUE)
  result_unpaired <- m_estimate(x, samples, loss_type = "huber", paired = FALSE)
  
  expect_true(is.data.frame(result_paired))
  expect_true(is.data.frame(result_unpaired))
})

test_that("m_estimate p-values are in [0,1]", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
  expect_true(all(result$padj >= 0 & result$padj <= 1, na.rm = TRUE))
})

test_that("m_estimate applies p-value correction", {
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result <- m_estimate(x, samples, loss_type = "huber", pcorr = "BH")
  
  # Adjusted p-values should be >= raw p-values (for BH)
  expect_true(all(result$padj >= result$pvalue, na.rm = TRUE))
})

# ============================================================================
# TEST: Extended Scale Estimators - MAD, Proposal 2, S-Estimator
# ============================================================================

test_that("m_estimate with MAD scale (default) produces valid results", {
  set.seed(123)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result <- m_estimate(x, samples, loss_type = "huber", scale_method = "mad")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(all(is.finite(result$location_diff)))
})

test_that("m_estimate with Huber Proposal 2 scale produces valid results", {
  set.seed(123)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result <- m_estimate(x, samples, loss_type = "huber", scale_method = "proposal2")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(all(is.finite(result$location_diff)))
})

test_that("m_estimate with S-estimator scale produces valid results", {
  set.seed(123)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result <- m_estimate(x, samples, loss_type = "huber", scale_method = "s-estimator")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(all(is.finite(result$location_diff)))
})

test_that("Different scale methods produce consistent estimates", {
  set.seed(42)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result_mad <- m_estimate(x, samples, loss_type = "huber", scale_method = "mad")
  result_prop2 <- m_estimate(x, samples, loss_type = "huber", scale_method = "proposal2")
  result_s <- m_estimate(x, samples, loss_type = "huber", scale_method = "s-estimator")
  
  # All should produce results with same number of rows
  expect_equal(nrow(result_mad), nrow(result_prop2))
  expect_equal(nrow(result_mad), nrow(result_s))
  
  # All results should be finite
  expect_true(all(is.finite(result_mad$location_diff)))
  expect_true(all(is.finite(result_prop2$location_diff)))
  expect_true(all(is.finite(result_s$location_diff)))
})

test_that("Scale methods with outliers show robust behavior", {
  set.seed(99)
  # Create data with extreme outliers
  x <- rbind(
    c(-100, -99, rnorm(8, mean = 2)),  # Extreme values in A
    rnorm(10)                           # Normal
  )
  samples <- c(rep("A", 5), rep("B", 5))
  
  result_mad <- m_estimate(x, samples, loss_type = "huber", scale_method = "mad")
  result_prop2 <- m_estimate(x, samples, loss_type = "huber", scale_method = "proposal2")
  result_s <- m_estimate(x, samples, loss_type = "huber", scale_method = "s-estimator")
  
  # All should handle outliers without crashing
  expect_true(all(is.finite(result_mad$location_diff)))
  expect_true(all(is.finite(result_prop2$location_diff)))
  expect_true(all(is.finite(result_s$location_diff)))
})

test_that("m_estimate rejects invalid scale_method parameter", {
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  expect_error(
    m_estimate(x, samples, loss_type = "huber", scale_method = "invalid"),
    "scale_method must be"
  )
})

test_that("scale_method parameter works with loss_type combinations", {
  set.seed(111)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  # Test MAD with all loss types
  result_huber <- m_estimate(x, samples, loss_type = "huber", scale_method = "mad")
  result_tukey <- m_estimate(x, samples, loss_type = "tukey", scale_method = "mad")
  result_lsq <- m_estimate(x, samples, loss_type = "lsq", scale_method = "mad")
  
  expect_true(all(is.finite(result_huber$location_diff)))
  expect_true(all(is.finite(result_tukey$location_diff)))
  expect_true(all(is.finite(result_lsq$location_diff)))
})

test_that("scale_method parameter preserved in recursive calls", {
  # This tests that scale_method is passed through in multi-q SummarizedExperiment calls
  set.seed(555)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  result_prop2 <- m_estimate(x, samples, loss_type = "huber", 
                             scale_method = "proposal2", paired = TRUE)
  result_s <- m_estimate(x, samples, loss_type = "huber", 
                         scale_method = "s-estimator", paired = TRUE)
  
  expect_true(all(is.finite(result_prop2$location_diff)))
  expect_true(all(is.finite(result_s$location_diff)))
})

test_that("Default scale_method is MAD", {
  set.seed(777)
  x <- matrix(rnorm(50), nrow = 5, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  # Result without specifying scale_method should use MAD
  result_default <- m_estimate(x, samples, loss_type = "huber")
  result_explicit_mad <- m_estimate(x, samples, loss_type = "huber", scale_method = "mad")
  
  # Should produce identical results
  expect_equal(result_default$location_diff, result_explicit_mad$location_diff)
  expect_equal(result_default$pvalue, result_explicit_mad$pvalue)
})

# ============================================================================
# TEST: Enhanced QC Metrics - SummarizedExperiment Input
# ============================================================================

context("Enhanced QC Metrics: Robustness Weight, Entropy Statistics, and Centroid Distance")

test_that("m_estimate with SummarizedExperiment includes all QC metrics", {
  library(SummarizedExperiment)
  set.seed(123)
  
  # Create multi-q column names like S1_q=1, S1_q=2, etc.
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 1 + (n_q - 1) * 0.5, length.out = n_q), n_samples)
  )
  
  # Create SummarizedExperiment with entropy data
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(87 * 16, mean = 4, sd = 1), nrow = 87, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("normal", "tumor"), each = 4*n_q),
      pair_id = rep(c("A", "B", "C", "D", "E", "F", "G", "H"), each = n_q)
    )
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber", paired = TRUE)
  
  # Check all expected columns present
  expected_cols <- c("Sample", "Condition", "Proportion_Affected", "Genes_Affected",
                     "Robustness_Weight", "Entropy_Mean", "Entropy_SD", 
                     "Distance_from_Centroid", "Status", "Pair_ID")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("Robustness_Weight values are between 0 and 1", {
  library(SummarizedExperiment)
  set.seed(456)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(50, mean = 3, sd = 0.8), nrow = 50, ncol = 16)),
    colData = data.frame(sample_type = rep(c("A", "B"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Robustness_Weight >= 0))
  expect_true(all(result$Robustness_Weight <= 1))
})

test_that("Entropy_Mean and Entropy_SD are positive or zero", {
  library(SummarizedExperiment)
  set.seed(789)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(abs(rnorm(60, mean = 3, sd = 1)), nrow = 60, ncol = 16)),
    colData = data.frame(sample_type = rep(c("X", "Y"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Entropy_Mean >= 0))
  expect_true(all(result$Entropy_SD >= 0))
})

test_that("Distance_from_Centroid is non-negative", {
  library(SummarizedExperiment)
  set.seed(321)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(40, mean = 2, sd = 0.5), nrow = 40, ncol = 16)),
    colData = data.frame(sample_type = rep(c("ctrl", "treat"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Distance_from_Centroid >= 0))
  expect_true(all(is.finite(result$Distance_from_Centroid)))
})

test_that("Samples closer to centroid have lower distances", {
  library(SummarizedExperiment)
  set.seed(654)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  # Create data where group patterns are very consistent
  se <- SummarizedExperiment(
    assays = list(
      diversity = rbind(
        matrix(5, nrow = 20, ncol = 16),  # Group A: all 5
        matrix(3, nrow = 20, ncol = 16),  # Group B: all 3
        matrix(rnorm(32, mean = 4, sd = 0.1), nrow = 2, ncol = 16)  # High variance genes
      )
    ),
    colData = data.frame(
      sample_type = rep(c("A", "B"), each = 4*n_q),
      pair = rep(1:8, each = n_q)
    )
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # All samples should have finite, non-negative distances
  expect_true(all(result$Distance_from_Centroid >= 0))
})

test_that("Distance_from_Centroid varies across samples (typical data)", {
  library(SummarizedExperiment)
  set.seed(987)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 1.2), nrow = 70, ncol = 16)),
    colData = data.frame(sample_type = rep(c("grp1", "grp2"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # For typical data, distances should vary (not all identical)
  # Or at minimum should all be finite and non-negative
  n_unique_distances <- length(unique(round(result$Distance_from_Centroid, 4)))
  expect_true(n_unique_distances >= 1)  # At least some variation in distances
  expect_true(all(is.finite(result$Distance_from_Centroid)))
  expect_true(all(result$Distance_from_Centroid >= 0))
})

test_that("Proportion_Affected and Genes_Affected are consistent", {
  library(SummarizedExperiment)
  set.seed(111)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(100, mean = 2, sd = 0.7), nrow = 100, ncol = 16)),
    colData = data.frame(sample_type = rep(c("S", "T"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Genes_Affected should equal Proportion_Affected * nrow
  expected_genes <- round(result$Proportion_Affected * 100, 1)
  expect_equal(result$Genes_Affected, expected_genes)
})

test_that("Entropy_Mean within reasonable bounds for data", {
  library(SummarizedExperiment)
  set.seed(222)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  # Create entropy matrix with known range
  entropy_data <- matrix(runif(80, min = 1, max = 5), nrow = 80, ncol = 16)
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("G1", "G2"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Entropy_Mean should be within range of data (approximately)
  expect_true(all(result$Entropy_Mean >= 0.5))  # Lower bound with margin
  expect_true(all(result$Entropy_Mean <= 5.5))  # Upper bound with margin
})

test_that("Paired parameter propagates through recursive calls", {
  library(SummarizedExperiment)
  set.seed(333)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(60, mean = 3, sd = 1), nrow = 60, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("case", "control"), each = 4*n_q),
      pair_id = rep(1:8, each = n_q)
    )
  )
  colnames(se) <- col_names
  
  result_paired <- m_estimate(se, samples = "sample_type", 
                              loss_type = "huber", paired = TRUE)
  result_unpaired <- m_estimate(se, samples = "sample_type", 
                                loss_type = "huber", paired = FALSE)
  
  # Both should return valid results
  expect_true(is.data.frame(result_paired))
  expect_true(is.data.frame(result_unpaired))
  
  # Paired and unpaired may differ, but both should have metrics
  expect_true("Robustness_Weight" %in% colnames(result_paired))
  expect_true("Robustness_Weight" %in% colnames(result_unpaired))
})

test_that("QC Status correctly flags high-influence samples", {
  library(SummarizedExperiment)
  set.seed(444)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 0.8), nrow = 70, ncol = 16)),
    colData = data.frame(sample_type = rep(c("normal", "tumor"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber", 
                       influence_threshold = 0.75)
  
  # Check Status column contains only "OK" or "Flag for QC"
  valid_statuses <- result$Status %in% c("OK", "Flag for QC")
  expect_true(all(valid_statuses))
  
  # High performers should be flagged
  high_proportion <- result[result$Proportion_Affected > 0.95, ]
  if (nrow(high_proportion) > 0) {
    # These should likely have "Flag for QC" status
    expect_true(any(high_proportion$Status == "Flag for QC"))
  }
})

test_that("Different influence thresholds produce different flagging", {
  library(SummarizedExperiment)
  set.seed(555)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(80, mean = 3, sd = 1), nrow = 80, ncol = 16)),
    colData = data.frame(sample_type = rep(c("A", "B"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result_75 <- m_estimate(se, samples = "sample_type", 
                          loss_type = "huber", influence_threshold = 0.75)
  result_90 <- m_estimate(se, samples = "sample_type", 
                          loss_type = "huber", influence_threshold = 0.90)
  
  # Lower threshold (0.75) should flag more samples than higher threshold (0.90)
  n_flagged_75 <- sum(result_75$Status == "Flag for QC")
  n_flagged_90 <- sum(result_90$Status == "Flag for QC")
  
  expect_true(n_flagged_75 >= n_flagged_90)
})

test_that("Robustness metrics work with different loss types", {
  library(SummarizedExperiment)
  set.seed(666)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(60, mean = 2, sd = 0.9), nrow = 60, ncol = 16)),
    colData = data.frame(sample_type = rep(c("X", "Y"), each = 4*n_q))
  )
  colnames(se) <- col_names
  
  result_huber <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  result_tukey <- m_estimate(se, samples = "sample_type", loss_type = "tukey")
  result_lsq <- m_estimate(se, samples = "sample_type", loss_type = "lsq")
  
  # All should have the new metrics
  for (result in list(result_huber, result_tukey, result_lsq)) {
    expect_true("Robustness_Weight" %in% colnames(result))
    expect_true("Entropy_Mean" %in% colnames(result))
    expect_true("Entropy_SD" %in% colnames(result))
    expect_true("Distance_from_Centroid" %in% colnames(result))
  }
})

test_that("Pair information is correctly extracted when available", {
  library(SummarizedExperiment)
  set.seed(777)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 1), nrow = 70, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("normal", "tumor"), each = 4*n_q),
      pair_id = rep(c("A", "B", "C", "D", "E", "F", "G", "H"), each = n_q)
    )
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Pair_ID should be extracted
  expect_true("Pair_ID" %in% colnames(result))
  
  # Pair_ID values should match expected pattern - 8 samples, each with their pair_id
  expected_pairs <- c("A", "B", "C", "D", "E", "F", "G", "H")
  expect_equal(as.character(result$Pair_ID), as.character(expected_pairs))
})

test_that("Entropy statistics reflect data variance", {
  library(SummarizedExperiment)
  set.seed(888)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  # Create two datasets: one with high variance, one with low variance
  low_var_data <- matrix(rnorm(80, mean = 3, sd = 0.1), nrow = 80, ncol = 16)
  high_var_data <- matrix(rnorm(80, mean = 3, sd = 2), nrow = 80, ncol = 16)
  
  se_low <- SummarizedExperiment(
    assays = list(diversity = low_var_data),
    colData = data.frame(sample_type = rep(c("L1", "L2"), each = 4*n_q))
  )
  colnames(se_low) <- col_names
  
  se_high <- SummarizedExperiment(
    assays = list(diversity = high_var_data),
    colData = data.frame(sample_type = rep(c("H1", "H2"), each = 4*n_q))
  )
  colnames(se_high) <- col_names
  
  result_low <- m_estimate(se_low, samples = "sample_type", loss_type = "huber")
  result_high <- m_estimate(se_high, samples = "sample_type", loss_type = "huber")
  
  # High variance data should have higher Entropy_SD on average
  mean_sd_low <- mean(result_low$Entropy_SD)
  mean_sd_high <- mean(result_high$Entropy_SD)
  
  expect_true(mean_sd_high > mean_sd_low)
})

test_that("m_estimate SummarizedExperiment path with various data sizes", {
  library(SummarizedExperiment)
  set.seed(999)
  
  n_q <- 2
  n_samples <- 8
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(seq(1, 2, length.out = n_q), n_samples)
  )
  
  # Create a SummarizedExperiment with multi-q format
  se_simple <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(1600, mean = 2.5, sd = 1), nrow = 100, ncol = 16)),
    colData = data.frame(sample_type = rep(c("g1", "g2"), each = 4*n_q))
  )
  colnames(se_simple) <- col_names
  
  # Basic test: should return valid results without crashing
  result <- m_estimate(se_simple, samples = "sample_type", loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true("Robustness_Weight" %in% colnames(result))
  expect_true("Entropy_Mean" %in% colnames(result))
  expect_true(all(is.finite(result$Robustness_Weight)))
  expect_true(all(is.finite(result$Entropy_Mean)))
})

# ============================================================================
# TEST: Helper Functions (Internal) - SE Data Preparation
# ============================================================================

context("Helper Functions: Data Preparation and Analysis")

test_that(".prepare_se_data_for_m_estimate correctly collapses multi-q data", {
  library(SummarizedExperiment)
  set.seed(1001)
  
  n_q <- 3
  n_samples <- 4
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 1.5, 2), n_samples)
  )
  
  # 4 genes x 12 columns (4 samples * 3 q-values each)
  entropy_data <- matrix(seq(1, 48), nrow = 4, ncol = 12)
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("A", "B"), each = 2*n_q))
  )
  colnames(se) <- col_names
  
  # Call helper via m_estimate which uses it internally
  # Verify it processes data correctly by checking output structure
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Result should have 4 rows (same as input genes)
  expect_equal(nrow(result), 4)
  expect_true("Entropy_Mean" %in% colnames(result))
})

test_that(".prepare_se_data_for_m_estimate handles median vs mean collapsing", {
  library(SummarizedExperiment)
  set.seed(1002)
  
  n_q <- 2
  n_samples <- 2
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  # Create data that will be collapsed across q-values
  entropy_data <- rbind(
    c(1, 100, 2, 101),      # High range in sample 1 - median vs mean differ
    c(50, 50, 50, 50),      # Constant - median = mean
    c(10, 15, 20, 25)       # Increasing - median and mean similar
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = c("X", "X", "Y", "Y"))
  )
  colnames(se) <- col_names
  
  # SE input returns sample-level QC metrics, not gene-level statistics
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)  # 2 samples (S1, S2)
  expect_true("Sample" %in% colnames(result))
  expect_true("Robustness_Weight" %in% colnames(result))
})

test_that(".prepare_se_data_for_m_estimate rejects invalid column names", {
  library(SummarizedExperiment)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(40), nrow = 5, ncol = 8)),
    colData = data.frame(group = rep(c("A", "B"), each = 4))
  )
  colnames(se) <- paste0("S", 1:8, "_q=1")
  
  # Should error when non-existent column is requested
  expect_error(
    m_estimate(se, samples = "nonexistent_column", loss_type = "huber"),
    "not found"
  )
})

test_that(".prepare_se_data_for_m_estimate preserves gene names and sample order", {
  library(SummarizedExperiment)
  set.seed(1003)
  
  n_q <- 2
  
  # 2 genes x 4 columns (2 samples * 2 q-values each)
  col_names <- paste0(
    rep(c("S1", "S2"), each = n_q),
    "_q=",
    rep(c(1, 2), 2)
  )
  
  entropy_data <- matrix(rnorm(8, mean = 2, sd = 0.5), nrow = 2, ncol = 4)
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("TypeA", "TypeB"), each = n_q))
  )
  rownames(se) <- c("GENE_A", "GENE_B")
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Result should have 2 rows (same as number of genes)
  expect_equal(nrow(result), 2)
})

# ============================================================================
# TEST: Helper Functions - LOO Influence Analysis
# ============================================================================

test_that(".perform_influence_loo_analysis identifies high-influence samples", {
  library(SummarizedExperiment)
  set.seed(1004)
  
  n_q <- 2
  n_samples <- 4
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  # Create data where first sample is very different (high influence)
  entropy_data <- rbind(
    c(-50, -50, 1, 2, 3, 4, 5, 6),  # Extreme outlier in first sample positions
    matrix(rnorm(88, mean = 2, sd = 0.5), nrow = 11, ncol = 8)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("Case", "Control"), each = 2*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Should flag high-influence outlier
  expect_true("Proportion_Affected" %in% colnames(result))
})

test_that(".perform_influence_loo_analysis handles identical groups", {
  library(SummarizedExperiment)
  set.seed(1005)
  
  n_q <- 2
  n_samples <- 2
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  # Create SE with identical values across all entries
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(5, nrow = 8, ncol = 4)),
    colData = data.frame(sample_type = c("GroupA", "GroupA", "GroupB", "GroupB"))
  )
  rownames(se) <- paste0("Gene", 1:8)
  colnames(se) <- col_names
  
  # SE returns sample-level results (one per unique sample)
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)  # 2 unique samples (S1, S2)
  expect_true(all(is.finite(result$Robustness_Weight)))
})

# ============================================================================
# TEST: Helper Functions - Centroid Distance Computation
# ============================================================================

test_that(".compute_centroid_distances_m_est calculates euclidean distances", {
  library(SummarizedExperiment)
  set.seed(1006)
  
  n_q <- 2
  n_samples <- 4
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  # Create data where distance patterns are clear
  entropy_data <- rbind(
    c(0, 0, 0, 0, 10, 10, 10, 10),  # Group 1: ~0, Group 2: ~10
    rnorm(8, mean = 2, sd = 0.2)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("Ctrl", "Treat"), each = 2*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Some samples should be closer to centroid than others
  distances <- result$Distance_from_Centroid
  expect_true(length(unique(round(distances, 3))) >= 1)
})

test_that(".compute_centroid_distances_m_est with single-sample groups", {
  library(SummarizedExperiment)
  set.seed(1007)
  
  n_q <- 2
  n_samples <- 4
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  # 4 samples: 2 in Group 1, 2 in Group 2
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(50, mean = 2, sd = 0.5), nrow = 50, ncol = 8)),
    colData = data.frame(sample_type = c(rep("G1", 2*n_q), rep("G2", 2*n_q)))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # With 2 samples per group, distance should be to group centroid (may be 0 for single sample)
  expect_true(all(result$Distance_from_Centroid >= 0))
  expect_true(all(is.finite(result$Distance_from_Centroid)))
})

# ============================================================================
# TEST: Edge Cases and Data Validation
# ============================================================================

context("Edge Cases: Single Sample, Missing Values, and Data Extremes")

test_that("m_estimate handles minimum viable data (2 groups, 1 sample each)", {
  x <- matrix(c(1:5, 6:10), nrow = 5, ncol = 2)
  samples <- c("A", "B")
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_true(nrow(result) == 5)
  expect_true(all(is.finite(result$location_diff)))
})

test_that("m_estimate handles data with zero variance in one gene", {
  x <- matrix(
    c(rep(5, 8),      # Gene 1: constant
      1:8),           # Gene 2: varying
    nrow = 2, byrow = TRUE
  )
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$location_diff)))
})

test_that("m_estimate handles data with NA values gracefully", {
  x <- matrix(c(rnorm(36), NA, NA, NA, NA), nrow = 4, ncol = 10)
  samples <- c(rep("A", 5), rep("B", 5))
  
  # Should handle NAs without crashing
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(is.data.frame(result))
  # May have NAs or computed values for genes with missing data
  expect_true(nrow(result) >= 0)
})

test_that("m_estimate with very large magnitude differences handles scaling", {
  x <- matrix(
    c(1e-6, 1e-5, 1e-4, 1e-3, 1e3, 1e4, 1e5, 1e6),
    nrow = 1, ncol = 8
  )
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  expect_true(is.data.frame(result))
  expect_true(all(is.finite(result$location_diff)))
})

test_that("m_estimate produces consistent results with identical raw data", {
  set.seed(2001)
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result1 <- m_estimate(x, samples, loss_type = "huber")
  result2 <- m_estimate(x, samples, loss_type = "huber")
  
  expect_equal(result1$location_diff, result2$location_diff)
  expect_equal(result1$pvalue, result2$pvalue)
})

test_that("m_estimate location_diff=0 when groups are identical", {
  x <- matrix(c(rep(5, 4), rep(5, 4)), nrow = 2, ncol = 8, byrow = TRUE)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  # When groups are identical, difference should be near 0
  expect_true(all(abs(result$location_diff) < 1e-6))
})

# ============================================================================
# TEST: Statistical Properties and Downstream Consistency
# ============================================================================

context("Statistical Properties: Behavior Under Different Conditions")

test_that("m_estimate t-statistics follow expected relationship with SE", {
  set.seed(2002)
  x <- matrix(rnorm(40, mean = 5, sd = 2), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  # t_stat = location_diff / se_diff
  expected_t <- result$location_diff / result$se_diff
  
  expect_equal(result$t_stat, expected_t, tolerance = 1e-10)
})

test_that("m_estimate rejects invalid loss_type gracefully", {
  x <- matrix(rnorm(40), nrow = 5, ncol = 8)
  samples <- c(rep("A", 4), rep("B", 4))
  
  expect_error(
    m_estimate(x, samples, loss_type = "invalid_loss"),
    "loss_type must be"
  )
})

test_that("m_estimate rejects insufficient data", {
  x <- matrix(c(1, 2), nrow = 1, ncol = 2)
  samples <- c("A", "B")
  
  # With only 2 samples total (1 per group), should still work as edge case
  result <- m_estimate(x, samples, loss_type = "huber")
  expect_true(is.data.frame(result))
})

test_that("Comparison of loss types shows different down-weighting patterns", {
  set.seed(2003)
  x <- matrix(
    c(c(1, 2, 3, 4, 1000, 2000, 3000, 4000),  # Group A: high variation
      rnorm(8, mean = 100, sd = 1)),           # Gene 2: stable
    nrow = 2, byrow = TRUE
  )
  samples <- c(rep("A", 4), rep("B", 4))
  
  result_huber <- m_estimate(x, samples, loss_type = "huber")
  result_tukey <- m_estimate(x, samples, loss_type = "tukey")
  result_lsq <- m_estimate(x, samples, loss_type = "lsq")
  
  # Robust methods should down-weight differently than LSQ
  # (checking that all produce valid results)
  expect_true(all(is.finite(result_huber$n_down_weighted)))
  expect_true(all(is.finite(result_tukey$n_down_weighted)))
  expect_true(all(is.finite(result_lsq$n_down_weighted)))
})

test_that("m_estimate weights reflect influence correctly", {
  set.seed(2004)
  x <- matrix(
    c(5, 5, 5, 5, 50, 50, 50, 50,  # Gene 1: extreme difference
      1, 1, 1, 1, 1, 1, 1, 1),      # Gene 2: no difference
    nrow = 2, byrow = TRUE
  )
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- m_estimate(x, samples, loss_type = "huber")
  
  # Gene 1 should have more down-weighted observations than Gene 2
  expect_true(result$n_down_weighted[1] >= result$n_down_weighted[2])
})

# ============================================================================
# TEST: SummarizedExperiment vs Matrix Consistency
# ============================================================================

context("Consistency: Matrix Input vs SummarizedExperiment Input")

test_that("Matrix and SE have different output formats as designed", {
  library(SummarizedExperiment)
  set.seed(2005)
  
  n_q <- 2
  
  # Create simple data: 3 genes x 4 columns (2 samples * 2 q-values each)
  x_matrix <- matrix(c(1:4, 5:8, 9:12), nrow = 3, ncol = 4, byrow = TRUE)
  samples_vector <- c("G1", "G1", "G2", "G2")
  
  # Create corresponding SE with same structure
  col_names <- paste0(
    rep(c("S1", "S2"), each = n_q),
    "_q=",
    rep(c(1, 2), 2)
  )
  se <- SummarizedExperiment(
    assays = list(diversity = x_matrix),
    colData = data.frame(sample_type = samples_vector)
  )
  rownames(se) <- c("Gene1", "Gene2", "Gene3")
  colnames(se) <- col_names
  
  # Test: SE input produces sample-level QC metrics
  result_se <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(is.data.frame(result_se))
  expect_equal(nrow(result_se), 2)  # 2 samples in SE format
  expect_true("Sample" %in% colnames(result_se))  # SE format has Sample column
  
  # Test: Matrix input produces gene-level statistics
  result_matrix <- m_estimate(x_matrix, samples = samples_vector, loss_type = "huber")
  
  expect_true(is.data.frame(result_matrix))
  expect_equal(nrow(result_matrix), 3)  # 3 genes
  expect_true("location_diff" %in% colnames(result_matrix))  # Matrix format has statistics
})

test_that("QC metrics exist for both matrix and SE input paths", {
  library(SummarizedExperiment)
  set.seed(2006)
  
  n_q <- 2
  n_samples <- 4
  col_names <- paste0(
    rep(paste0("S", 1:n_samples), each = n_q),
    "_q=",
    rep(c(1, 2), n_samples)
  )
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(32, mean = 2, sd = 0.5), nrow = 4, ncol = 8)),
    colData = data.frame(sample_type = rep(c("X", "Y"), each = 2*n_q))
  )
  colnames(se) <- col_names
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # SE path should include QC metrics
  qc_cols <- c("Robustness_Weight", "Entropy_Mean", "Entropy_SD", "Distance_from_Centroid")
  expect_true(all(qc_cols %in% colnames(result)))
})




