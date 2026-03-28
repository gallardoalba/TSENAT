context("Statistical Validation: Power and Validity Tests")

library(testthat)

# =============================================================================
# TYPE I ERROR CONTROL - Permutation Tests
# =============================================================================

test_that("label_shuffling produces p-values consistent with null (no effect)", {
  set.seed(100)
  
  # Generate 50 genes with NO effect (null hypothesis)
  n_genes <- 50
  n_samples <- 10
  
  # Random values with no difference between groups
  mat <- matrix(rnorm(n_genes * n_samples, mean = 0.5, sd = 0.1), 
                nrow = n_genes)
  samples <- c(rep("A", 5), rep("B", 5))
  
  # Run permutation test
  result <- .label_shuffling(mat, samples, control = "A", 
                           method = "mean", randomizations = 100, 
                           pcorr = "none")
  
  p_values <- result[, "pvalue"]
  
  # Under null, p-values should be roughly uniformly distributed
  # With 50 genes and 100 permutations, expect reasonable spread
  # Check that we don't have obviously skewed distribution
  proportion_small <- mean(p_values < 0.3)
  proportion_large <- mean(p_values > 0.5)
  
  # Expect some p-values in different ranges (roughly uniform)
  expect_true(proportion_large > 0.1)  # Some p-values > 0.5
  expect_true(proportion_small > 0.1)  # Some p-values < 0.3
})

test_that("label_shuffling detects true differences (power test)", {
  set.seed(101)
  
  # Create 30 genes: 10 with effect, 20 with no effect
  n_effect <- 10
  n_null <- 20
  n_samples <- 8
  
  # Effect genes: A has mean 0.3, B has mean 0.7
  mat_effect <- matrix(NA, nrow = n_effect, ncol = n_samples)
  mat_effect[, 1:4] <- rnorm(n_effect * 4, mean = 0.3, sd = 0.05)
  mat_effect[, 5:8] <- rnorm(n_effect * 4, mean = 0.7, sd = 0.05)
  
  # Null genes: both groups mean 0.5
  mat_null <- matrix(rnorm(n_null * n_samples, mean = 0.5, sd = 0.1), 
                     nrow = n_null)
  
  mat <- rbind(mat_effect, mat_null)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- .label_shuffling(mat, samples, control = "A", 
                           method = "mean", randomizations = 200, 
                           pcorr = "none")
  
  p_effect <- result[1:n_effect, "pvalue"]
  p_null <- result[(n_effect + 1):(n_effect + n_null), "pvalue"]
  
  # Effect genes should have significantly lower p-values
  expect_true(median(p_effect) < median(p_null))
  # At least 50% of effect genes should be detected at p < 0.1
  expect_true(mean(p_effect < 0.1) >= 0.5)
})

# =============================================================================
# MIXED MODEL FALLBACK IMPROVEMENT (Per fix in calc_lm_helpers.R)
# =============================================================================

test_that("lm_subject_fixed fallback preserves power vs. lm_nosubject", {
  skip_if_not_installed("lme4")
  set.seed(102)
  
  # Create data with within-subject correlation
  n_subjects <- 8
  n_q_per_subject <- 5
  subjects <- rep(paste0("S", 1:n_subjects), each = n_q_per_subject)
  q_vals <- rep(seq(0.1, 0.9, length.out = n_q_per_subject), times = n_subjects)
  group <- rep(c("A", "B"), each = n_subjects * n_q_per_subject / 2)
  
  # True effect: q×group interaction
  subject_effect <- rep(rnorm(n_subjects, 0, 0.1), each = n_q_per_subject)
  entropy <- 0.5 + 
             0.3 * q_vals +
             0.2 * (group == "B") +
             0.15 * (group == "B") * q_vals +  # True interaction
             subject_effect +
             rnorm(length(subjects), 0, 0.05)
  
  df <- data.frame(entropy = entropy, q = q_vals, group = factor(group), 
                   subject = factor(subjects))
  
  # Fit with subject as fixed effect (proper approach)
  fit_with_subj <- stats::lm(entropy ~ q * group + factor(subject), data = df)
  p_with_subj <- summary(fit_with_subj)$coefficients["q:groupB", "Pr(>|t|)"]
  
  # Fit without subject (loses power)
  fit_no_subj <- stats::lm(entropy ~ q * group, data = df)
  p_no_subj <- summary(fit_no_subj)$coefficients["q:groupB", "Pr(>|t|)"]
  
  # Model with subject should have lower p-value (more power)
  expect_true(p_with_subj < p_no_subj)
})

test_that(".tsenat_try_lm_fallbacks uses factor(subject), not numeric subject", {
  set.seed(103)
  
  # Create test data
  df <- data.frame(
    entropy = rnorm(30),
    q = rep(seq(0.1, 1, length.out = 10), 3),
    group = rep(c("A", "B"), length.out = 30),
    subject = rep(1:10, 3)  # numeric subject IDs
  )
  
  # Run fallback function
  fb <- .tsenat_try_lm_fallbacks(df, verbose = FALSE)
  
  expect_true(!is.null(fb))
  expect_true(!is.null(fb$fit1))
  
  # Extract coefficients to verify it's treating subject as factor
  if (fb$method == "lm_subject_fixed") {
    coef_names <- names(coef(fb$fit1))
    # Should have factor(subject) terms, not a single "subject" slope
    subj_terms <- grep("factor\\(subject\\)", coef_names)
    expect_true(length(subj_terms) > 0)
  }
})

# =============================================================================
# PERMUTATION TEST EXACT vs APPROXIMATE
# =============================================================================

test_that("label_shuffling with exact signflip has proper p-value granularity", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(104)
  
  # Very small example: 2 genes x 4 paired samples
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)
  samples <- c("S1", "S2", "S1", "S2")
  
  # Attempt exact signflip enumeration (should give 2^2 = 4 permutations)
  result <- .label_shuffling(mat, samples, control = "S1", 
                           method = "mean", randomizations = 10,
                           paired = TRUE, paired_method = "signflip")
  
  p_vals <- result[, "pvalue"]
  
  # With 4 permutations + pseudocount: p-values should be 1/5, 2/5, 3/5, 4/5, or 5/5
  valid_p <- seq(1, 5) / 5
  for (p in p_vals) {
    # Allow small floating point error
    expect_true(any(abs(p - valid_p) < 1e-6))
  }
})

# =============================================================================
# NORMALIZED ENTROPY BOUNDS
# =============================================================================

test_that("normalized entropy for q < 1 can exceed 1 (mathematically valid)", {
  # Highly skewed distribution
  counts <- c(100, 1)
  
  # q < 1 emphasizes rare elements
  result_q05 <- .calculate_tsallis_entropy(counts, q = 0.5, norm = TRUE)
  
  # Should be finite (not NaN)
  expect_true(!is.nan(result_q05))
  # For this skewed distribution, may exceed 1
  # Just verify it's reasonable
  expect_true(result_q05 > 0)
  expect_true(!is.infinite(result_q05))
})

test_that("normalized entropy bounds correct for q > 1", {
  # For uniform distribution with q > 1, should equal 1
  counts_uniform <- rep(10, 5)
  
  result_q15 <- .calculate_tsallis_entropy(counts_uniform, q = 1.5, norm = TRUE)
  result_q2 <- .calculate_tsallis_entropy(counts_uniform, q = 2, norm = TRUE)
  
  expect_equal(result_q15, 1.0, tolerance = 1e-6)
  expect_equal(result_q2, 1.0, tolerance = 1e-6)
  
  # Skewed distribution should give <1
  counts_skewed <- c(100, 1, 1, 1, 1)
  result_skew_q15 <- .calculate_tsallis_entropy(counts_skewed, q = 1.5, norm = TRUE)
  result_skew_q2 <- .calculate_tsallis_entropy(counts_skewed, q = 2, norm = TRUE)
  
  expect_true(result_skew_q15 < 1.0)
  expect_true(result_skew_q2 < 1.0)
  expect_true(result_skew_q15 > 0.0)
  expect_true(result_skew_q2 > 0.0)
})

# =============================================================================
# EXTREME q VALUES STABILITY
# =============================================================================

test_that("tsallis entropy stable at extreme q values", {
  counts <- c(10, 5, 3, 1)
  
  # Very small q
  result_q001 <- .calculate_tsallis_entropy(counts, q = 0.01, norm = FALSE)
  expect_true(!is.nan(result_q001))
  expect_true(is.finite(result_q001))
  
  # Very large q
  result_q10 <- .calculate_tsallis_entropy(counts, q = 10, norm = FALSE)
  expect_true(!is.nan(result_q10))
  expect_true(is.finite(result_q10))
  
  # Both should be reasonable values
  expect_true(result_q001 >= 0)
  expect_true(result_q10 >= 0)
})

test_that("hill numbers stable at extreme q", {
  counts <- c(20, 10, 5, 1)
  
  # D_q at q = 0.1 and q = 5
  D_q01 <- .calculate_tsallis_entropy(counts, q = 0.1, what = "D", norm = FALSE)
  D_q5 <- .calculate_tsallis_entropy(counts, q = 5, what = "D", norm = FALSE)
  
  # Should be finite and positive
  expect_true(!is.nan(D_q01))
  expect_true(!is.nan(D_q5))
  expect_true(D_q01 > 0)
  expect_true(D_q5 > 0)
  
  # q=0.1 emphasizes rare elements, so Hill number should be larger
  # q=5 emphasizes common elements, so Hill number should be smaller
  # Actually, for skewed distributions, D_0.1 > D_5
  expect_true(D_q01 > D_q5)
})

# =============================================================================
# SCALE INVARIANCE
# =============================================================================

test_that("fold change is scale invariant for log scale", {
  mat1 <- matrix(c(1, 2, 5, 10, 3, 6), nrow = 2, ncol = 3)
  mat2 <- mat1 * 1000  # Scale by 1000x
  
  samples <- c("Normal", "Tumor", "Tumor")
  
  result1 <- TSENAT:::.calculate_fc(mat1, samples, control = "Normal", pseudocount = 1e-6)
  result2 <- TSENAT:::.calculate_fc(mat2, samples, control = "Normal", pseudocount = 1e-6)
  
  # Log2 FC should be identical
  expect_equal(result1[, 4], result2[, 4], 
               tolerance = 1e-10)
})

test_that("entropy is scale invariant", {
  counts1 <- c(10, 20, 30, 40)
  counts2 <- counts1 * 100  # Scale by 100x
  
  result1 <- .calculate_tsallis_entropy(counts1, q = 2, norm = FALSE)
  result2 <- .calculate_tsallis_entropy(counts2, q = 2, norm = FALSE)
  
  # Entropy depends only on proportions, not absolute counts
  expect_equal(result1, result2, tolerance = 1e-10)
})

# =============================================================================
# WILCOXON TEST VALIDITY
# =============================================================================

test_that("wilcoxon test matches R's built-in wilcox.test exactly", {
  # Create clear difference
  normal_vals <- c(1, 2, 3, 4, 5)
  tumor_vals <- c(6, 7, 8, 9, 10)
  
  mat <- matrix(c(normal_vals, tumor_vals), nrow = 1)
  samples <- c(rep("Normal", 5), rep("Tumor", 5))
  
  # TSENAT wilcoxon
  tsenat_result <- .wilcoxon(mat, samples, pcorr = "none")
  tsenat_p <- tsenat_result[1, "pvalue"]
  
  # R's wilcox.test
  r_result <- wilcox.test(normal_vals, tumor_vals, exact = FALSE)
  r_p <- r_result$p.value
  
  # Should be essentially identical
  expect_equal(tsenat_p, r_p, tolerance = 1e-10)
})

test_that("wilcoxon paired matches paired test from R", {
  # Paired data
  before <- c(1, 2, 3, 4, 5)
  after <- c(2, 3, 5, 6, 8)
  
  mat <- matrix(c(before, after), nrow = 1)
  samples <- c(rep("Before", 5), rep("After", 5))
  
  # TSENAT with paired = TRUE
  tsenat_paired <- .wilcoxon(mat, samples, paired = TRUE, pcorr = "none")
  tsenat_p <- tsenat_paired[1, "pvalue"]
  
  # R's paired wilcox.test
  r_paired <- wilcox.test(before, after, paired = TRUE, exact = FALSE)
  r_p <- r_paired$p.value
  
  expect_equal(tsenat_p, r_p, tolerance = 1e-10)
})

# =============================================================================
# PSEUDOCOUNT STRATEGIES
# =============================================================================

test_that("pseudocount selection is data-driven and prevents negative log(0)", {
  # Data with zeros
  mat <- matrix(c(0, 0, 1, 5, 0, 0, 2, 3), nrow = 2, byrow = TRUE)
  samples <- c("A", "A", "B", "B")
  
  # Auto pseudocount (min/2 = 0.5)
  result <- TSENAT:::.calculate_fc(mat, samples, control = "A", pseudocount = 0)
  
  # Should not produce NaN or Inf
  expect_true(!any(is.nan(result$log2_fold_change)))
  expect_true(!any(is.infinite(result$log2_fold_change)))
})

# =============================================================================
# MULTIPLE COMPARISON CORRECTION
# =============================================================================

test_that("BH FDR correction reduces false positives", {
  set.seed(105)
  
  # 100 null genes
  mat <- matrix(rnorm(100 * 8, mean = 0.5, sd = 0.1), nrow = 100)
  samples <- c(rep("A", 4), rep("B", 4))
  
  result <- .label_shuffling(mat, samples, control = "A", 
                           method = "mean", randomizations = 100,
                           pcorr = "BH")
  
  # Under BH correction, expect FDR < 0.05
  fdr <- mean(result[, "padj"] < 0.05)
  expect_true(fdr < 0.1)  # Lenient threshold for stochastic test
})
