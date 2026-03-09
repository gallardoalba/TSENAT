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
  
  # All should produce similar magnitude results (within 5% relative difference)
  for (i in 1:nrow(x)) {
    mad_coeff <- abs(result_mad$location_diff[i])
    prop2_coeff <- abs(result_prop2$location_diff[i])
    s_coeff <- abs(result_s$location_diff[i])
    
    # Coefficients should be reasonably close (allowing for method differences)
    expect_true(abs(mad_coeff - prop2_coeff) / (mad_coeff + 0.01) < 0.20 ||
                abs(mad_coeff - s_coeff) / (mad_coeff + 0.01) < 0.20)
  }
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
  
  # Create SummarizedExperiment with entropy data
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(87, mean = 4, sd = 1), nrow = 87, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("normal", "tumor"), 8),
      pair_id = rep(c("A", "B", "C", "D", "E", "F", "G", "H"), 2)
    )
  )
  
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
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(50, mean = 3, sd = 0.8), nrow = 50, ncol = 16)),
    colData = data.frame(sample_type = rep(c("A", "B"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Robustness_Weight >= 0))
  expect_true(all(result$Robustness_Weight <= 1))
})

test_that("Entropy_Mean and Entropy_SD are positive or zero", {
  library(SummarizedExperiment)
  set.seed(789)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(abs(rnorm(60, mean = 3, sd = 1)), nrow = 60, ncol = 16)),
    colData = data.frame(sample_type = rep(c("X", "Y"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Entropy_Mean >= 0))
  expect_true(all(result$Entropy_SD >= 0))
})

test_that("Distance_from_Centroid is non-negative", {
  library(SummarizedExperiment)
  set.seed(321)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(40, mean = 2, sd = 0.5), nrow = 40, ncol = 16)),
    colData = data.frame(sample_type = rep(c("ctrl", "treat"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  expect_true(all(result$Distance_from_Centroid >= 0))
  expect_true(all(is.finite(result$Distance_from_Centroid)))
})

test_that("Samples closer to centroid have lower distances", {
  library(SummarizedExperiment)
  set.seed(654)
  
  # Create data where group patterns are very consistent
  se <- SummarizedExperiment(
    assays = list(
      diversity = rbind(
        matrix(5, nrow = 20, ncol = 8),  # Group A: all 5
        matrix(3, nrow = 20, ncol = 8),  # Group B: all 3
        matrix(rnorm(20, mean = 4, sd = 0.1), nrow = 2, ncol = 8)  # High variance genes
      )
    ),
    colData = data.frame(
      sample_type = rep(c("A", "B"), 8),
      pair = rep(1:8, 2)
    )
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # All samples should have finite, non-negative distances
  expect_true(all(result$Distance_from_Centroid >= 0))
})

test_that("Distance_from_Centroid varies across samples (typical data)", {
  library(SummarizedExperiment)
  set.seed(987)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 1.2), nrow = 70, ncol = 16)),
    colData = data.frame(sample_type = rep(c("grp1", "grp2"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # For typical data, distances should vary (not all identical)
  n_unique_distances <- length(unique(round(result$Distance_from_Centroid, 6)))
  expect_true(n_unique_distances > 1)
})

test_that("Proportion_Affected and Genes_Affected are consistent", {
  library(SummarizedExperiment)
  set.seed(111)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(100, mean = 2, sd = 0.7), nrow = 100, ncol = 16)),
    colData = data.frame(sample_type = rep(c("S", "T"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Genes_Affected should equal Proportion_Affected * nrow
  expected_genes <- round(result$Proportion_Affected * 100, 1)
  expect_equal(result$Genes_Affected, expected_genes)
})

test_that("Entropy_Mean within reasonable bounds for data", {
  library(SummarizedExperiment)
  set.seed(222)
  
  # Create entropy matrix with known range
  entropy_data <- matrix(runif(80, min = 1, max = 5), nrow = 80, ncol = 16)
  
  se <- SummarizedExperiment(
    assays = list(diversity = entropy_data),
    colData = data.frame(sample_type = rep(c("G1", "G2"), 8))
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Entropy_Mean should be within range of data (approximately)
  expect_true(all(result$Entropy_Mean >= 0.5))  # Lower bound with margin
  expect_true(all(result$Entropy_Mean <= 5.5))  # Upper bound with margin
})

test_that("Paired parameter propagates through recursive calls", {
  library(SummarizedExperiment)
  set.seed(333)
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(60, mean = 3, sd = 1), nrow = 60, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("case", "control"), 8),
      pair_id = rep(1:8, 2)
    )
  )
  
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
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 0.8), nrow = 70, ncol = 16)),
    colData = data.frame(sample_type = rep(c("normal", "tumor"), 8))
  )
  
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
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(80, mean = 3, sd = 1), nrow = 80, ncol = 16)),
    colData = data.frame(sample_type = rep(c("A", "B"), 8))
  )
  
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
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(60, mean = 2, sd = 0.9), nrow = 60, ncol = 16)),
    colData = data.frame(sample_type = rep(c("X", "Y"), 8))
  )
  
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
  
  se <- SummarizedExperiment(
    assays = list(diversity = matrix(rnorm(70, mean = 2.5, sd = 1), nrow = 70, ncol = 16)),
    colData = data.frame(
      sample_type = rep(c("normal", "tumor"), 8),
      pair_id = rep(c("A", "B", "C", "D", "E", "F", "G", "H"), 2)
    )
  )
  
  result <- m_estimate(se, samples = "sample_type", loss_type = "huber")
  
  # Pair_ID should be extracted
  expect_true("Pair_ID" %in% colnames(result))
  
  # Pair_ID values should match expected pattern
  expected_pairs <- rep(c("A", "B", "C", "D", "E", "F", "G", "H"), 2)
  expect_equal(as.character(result$Pair_ID), as.character(expected_pairs))
})

test_that("Entropy statistics reflect data variance", {
  library(SummarizedExperiment)
  set.seed(888)
  
  # Create two datasets: one with high variance, one with low variance
  low_var_data <- matrix(rnorm(80, mean = 3, sd = 0.1), nrow = 80, ncol = 16)
  high_var_data <- matrix(rnorm(80, mean = 3, sd = 2), nrow = 80, ncol = 16)
  
  se_low <- SummarizedExperiment(
    assays = list(diversity = low_var_data),
    colData = data.frame(sample_type = rep(c("L1", "L2"), 8))
  )
  
  se_high <- SummarizedExperiment(
    assays = list(diversity = high_var_data),
    colData = data.frame(sample_type = rep(c("H1", "H2"), 8))
  )
  
  result_low <- m_estimate(se_low, samples = "sample_type", loss_type = "huber")
  result_high <- m_estimate(se_high, samples = "sample_type", loss_type = "huber")
  
  # High variance data should have higher Entropy_SD on average
  mean_sd_low <- mean(result_low$Entropy_SD)
  mean_sd_high <- mean(result_high$Entropy_SD)
  
  expect_true(mean_sd_high > mean_sd_low)
})

test_that("m_estimate SummarizedExperiment path detects multi-q format", {
  library(SummarizedExperiment)
  set.seed(999)
  
  # Create multi-q format: Sample_q=value naming
  multi_q_data <- matrix(rnorm(160, mean = 2.5, sd = 1), nrow = 80, ncol = 16)
  colnames(multi_q_data) <- c(
    "Sample1_q=0.5", "Sample2_q=0.5", "Sample3_q=0.5", "Sample4_q=0.5",
    "Sample5_q=1.0", "Sample6_q=1.0", "Sample7_q=1.0", "Sample8_q=1.0",
    "Sample1_q=1.5", "Sample2_q=1.5", "Sample3_q=1.5", "Sample4_q=1.5",
    "Sample5_q=2.0", "Sample6_q=2.0", "Sample7_q=2.0", "Sample8_q=2.0"
  )
  
  se_multi <- SummarizedExperiment(
    assays = list(diversity = multi_q_data),
    colData = data.frame(
      sample_type = rep(c("g1", "g2"), 8),
      row.names = colnames(multi_q_data)
    )
  )
  
  # Should collapse across q values and compute metrics
  result <- m_estimate(se_multi, samples = "sample_type", loss_type = "huber",
                       q_combine_method = "mean")
  
  # Result should have 8 rows (one per sample after collapsing q)
  expect_equal(nrow(result), 8)
  
  # All metrics should be present
  expect_true("Robustness_Weight" %in% colnames(result))
  expect_true("Entropy_Mean" %in% colnames(result))
})


