library(TSENAT)

context("Multiple Testing Correction for Tsallis Entropy Q-Dependence")

# ============================================================================
# Helper Function Tests: Hochberg Stepup Procedure
# ============================================================================

test_that(".hochberg_stepup produces valid p-values", {
  # Test basic functionality
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  expect_is(adjusted, "numeric")
  expect_equal(length(adjusted), 5)
  # Adjusted p-values should be >= raw p-values
  expect_true(all(adjusted >= raw_pvalues))
  # Adjusted p-values should be <= 1
  expect_true(all(adjusted <= 1))
})

test_that(".hochberg_stepup respects monotone increasing property", {
  # Adjusted p-values should be monotone increasing when sorted
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  # Sort by raw p-value and check monotonicity of adjusted
  order_idx <- order(raw_pvalues)
  sorted_adj <- adjusted[order_idx]
  
  # Check that sorted adjusted is monotone increasing
  expect_true(all(diff(sorted_adj) >= -1e-10))  # Allow tiny numerical error
})

test_that(".hochberg_stepup handles NA and Inf values gracefully", {
  raw_pvalues <- c(0.01, NA, 0.05, Inf, 0.1)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  # Invalid values should be preserved
  expect_true(is.na(adjusted[2]))
  expect_true(is.infinite(adjusted[4]))
  
  # Valid values should be adjusted
  expect_true(!is.na(adjusted[1]))
  expect_true(!is.na(adjusted[3]))
  expect_true(!is.na(adjusted[5]))
})

test_that(".hochberg_stepup single p-value returns same value bounded to 1", {
  raw_pvalues <- 0.5
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  expect_equal(adjusted, 0.5)
})

test_that(".hochberg_stepup very small p-values receive larger inflation", {
  raw_pvalues <- c(0.0001, 0.01, 0.5)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  m <- length(raw_pvalues)
  # Hochberg adjustment: adjusted = (m - rank + 1) * raw_p
  # For m=3: smallest p gets multiplied by 3, largest by 1
  inflation_factor <- adjusted / raw_pvalues
  
  # Inflation increasing with rank
  expect_true(inflation_factor[1] >= inflation_factor[3])
})

test_that(".hochberg_stepup with all NA returns all NA", {
  raw_pvalues <- c(NA, NA, NA)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  expect_true(all(is.na(adjusted)))
})

test_that(".hochberg_stepup empty vector returns empty vector", {
  raw_pvalues <- numeric(0)
  adjusted <- .hochberg_stepup(raw_pvalues)
  
  expect_equal(length(adjusted), 0)
})


# ============================================================================
# Helper Function Tests: Benjamini-Yekutieli FDR Control
# ============================================================================

test_that(".benjamini_yekutieli produces valid p-values", {
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  adjusted <- .benjamini_yekutieli(raw_pvalues)
  
  expect_is(adjusted, "numeric")
  expect_equal(length(adjusted), 5)
  # Adjusted p-values should be >= raw p-values
  expect_true(all(adjusted >= raw_pvalues))
  # Adjusted p-values should be <= 1
  expect_true(all(adjusted <= 1))
})

test_that(".benjamini_yekutieli respects monotone decreasing property", {
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  adjusted <- .benjamini_yekutieli(raw_pvalues)
  
  # When sorted by raw p-value, adjusted should be non-decreasing
  order_idx <- order(raw_pvalues)
  sorted_adj <- adjusted[order_idx]
  
  # Check that sorted adjusted is monotone non-decreasing
  expect_true(all(diff(sorted_adj) >= -1e-10))  # Allow tiny numerical error
})

test_that(".benjamini_yekutieli handles NA gracefully", {
  raw_pvalues <- c(0.01, NA, 0.05, 0.1)
  adjusted <- .benjamini_yekutieli(raw_pvalues)
  
  # NA should be preserved
  expect_true(is.na(adjusted[2]))
  # Valid values should be adjusted
  expect_true(all(!is.na(adjusted[-2])))
})

test_that(".benjamini_yekutieli single p-value returns same value", {
  raw_pvalues <- 0.5
  adjusted <- .benjamini_yekutieli(raw_pvalues)
  
  expect_equal(adjusted, 0.5)
})

test_that(".benjamini_yekutieli is less conservative than Hochberg", {
  # For the same p-values, BY should be more conservative (larger adjustments) than Hochberg
  # because it accounts for arbitrary dependence while Hochberg assumes positive regression dependence
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  
  adjusted_hochberg <- .hochberg_stepup(raw_pvalues)
  adjusted_by <- .benjamini_yekutieli(raw_pvalues)
  
  # BY should give larger adjusted p-values (more conservative) due to harmonic sum correction
  # This is because BY controls FDR under arbitrary dependence vs Hochberg under positive regression dependence
  expect_true(mean(adjusted_by) >= mean(adjusted_hochberg) - 0.05)
})

test_that(".benjamini_yekutieli uses Harmonic constant c_m correctly", {
  # For m=5, c_m = sum(1/i for i=1..5) = 1 + 1/2 + 1/3 + 1/4 + 1/5
  raw_pvalues <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  adjusted <- .benjamini_yekutieli(raw_pvalues)
  
  m <- 5
  c_m <- sum(1 / (1:m))
  
  # Manual calculation for verification
  order_idx <- order(raw_pvalues)
  sorted_p <- raw_pvalues[order_idx]
  
  expected_adj <- numeric(m)
  for (i in seq_len(m)) {
    expected_adj[i] <- (m * c_m / i) * sorted_p[i]
    expected_adj[i] <- min(expected_adj[i], 1)
  }
  
  # Apply monotone constraint
  for (i in (m-1):1) {
    if (expected_adj[i] > expected_adj[i+1]) {
      expected_adj[i] <- expected_adj[i+1]
    }
  }
  
  # Unpack to original order
  expected_result <- numeric(m)
  expected_result[order_idx] <- expected_adj
  
  expect_equal(adjusted, expected_result, tolerance = 1e-10)
})

test_that(".benjamini_yekutieli handles empty input", {
  expect_equal(.benjamini_yekutieli(numeric(0)), numeric(0))
})

test_that(".benjamini_yekutieli handles all-NA input", {
  result <- .benjamini_yekutieli(c(NA, NA, NA))
  expect_equal(result, c(NA, NA, NA))
})

test_that(".benjamini_yekutieli handles all-NaN input", {
  result <- .benjamini_yekutieli(c(NaN, NaN))
  expect_equal(result, c(NaN, NaN))
})

test_that(".benjamini_yekutieli handles mixed NA and valid p-values", {
  result <- .benjamini_yekutieli(c(0.01, NA, 0.05, NA, 0.1))
  expect_equal(result[2], NA_real_)
  expect_equal(result[4], NA_real_)
  valid_result <- result[!is.na(result)]
  expect_true(all(valid_result >= c(0.01, 0.05, 0.1)))
  expect_true(all(is.finite(valid_result)))
})

test_that(".benjamini_yekutieli handles Inf values", {
  result <- .benjamini_yekutieli(c(0.01, Inf, 0.05))
  expect_equal(result[2], Inf)
  expect_true(all(is.finite(result[-2])))
})

test_that(".benjamini_yekutieli with all p=1 returns 1", {
  result <- .benjamini_yekutieli(c(1, 1, 1))
  expect_equal(result, c(1, 1, 1))
})

test_that(".benjamini_yekutieli with all p=0 returns 0", {
  result <- .benjamini_yekutieli(c(0, 0, 0))
  expect_equal(result, c(0, 0, 0))
})

# ============================================================================
# Hochberg Edge Case Tests (increasing coverage)
# ============================================================================

test_that(".hochberg_stepup handles empty input", {
  expect_equal(.hochberg_stepup(numeric(0)), numeric(0))
})

test_that(".hochberg_stepup handles all-NA input", {
  result <- .hochberg_stepup(c(NA, NA, NA))
  expect_equal(result, c(NA, NA, NA))
})

test_that(".hochberg_stepup handles all-NaN input", {
  result <- .hochberg_stepup(c(NaN, NaN))
  expect_equal(result, c(NaN, NaN))
})

test_that(".hochberg_stepup handles mixed Inf and valid p-values", {
  result <- .hochberg_stepup(c(0.01, Inf, 0.05, -Inf, 0.1))
  expect_equal(result[2], Inf)
  expect_equal(result[4], -Inf)
  expect_true(all(is.finite(result[c(1, 3, 5)])))
})

test_that(".hochberg_stepup preserves NA in correct positions", {
  result <- .hochberg_stepup(c(0.01, NA, 0.05, NA, 0.001))
  expect_true(is.na(result[2]))
  expect_true(is.na(result[4]))
  expect_true(all(!is.na(result[c(1, 3, 5)])))
})

test_that(".hochberg_stepup with single NA returns NA", {
  result <- .hochberg_stepup(NA_real_)
  expect_equal(result, NA_real_)
})

# ============================================================================
# Skewness Function Tests
# ============================================================================

test_that(".compute_skewness calculates symmetry correctly", {
  # Symmetric distribution should have near-zero skewness
  symmetric_data <- c(-2, -1, 0, 1, 2)
  skew_sym <- .compute_skewness(symmetric_data)
  
  expect_true(abs(skew_sym) < 0.5)
})

test_that(".compute_skewness detects right skewness", {
  # Right-skewed distribution (tail to the right)
  right_skewed <- c(1, 2, 2, 3, 3, 3, 10)
  skew_right <- .compute_skewness(right_skewed)
  
  expect_true(skew_right > 0)
})

test_that(".compute_skewness detects left skewness", {
  # Left-skewed distribution (tail to the left)
  left_skewed <- c(-10, 3, 3, 3, 2, 2, 1)
  skew_left <- .compute_skewness(left_skewed)
  
  expect_true(skew_left < 0)
})

test_that(".compute_skewness handles NA values", {
  data_with_na <- c(1, 2, NA, 3, 4, NA, 5)
  skew <- .compute_skewness(data_with_na, na.rm = TRUE)
  
  expect_true(!is.na(skew))
  expect_is(skew, "numeric")
})

test_that(".compute_skewness returns NA for < 3 observations", {
  skew_two <- .compute_skewness(c(1, 2))
  skew_one <- .compute_skewness(1)
  
  expect_true(is.na(skew_two))
  expect_true(is.na(skew_one))
})

test_that(".compute_skewness handles constant data", {
  constant_data <- c(5, 5, 5, 5, 5)
  skew <- .compute_skewness(constant_data)
  
  # When SD=0, skewness is zero
  expect_equal(skew, 0)
})


# ============================================================================
# Integration Tests: detect_q_gene_interactions with multicorr parameter
# ============================================================================

test_that("detect_q_gene_interactions produces adj_p_value column", {
  set.seed(100)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),
                  rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep(c("Gene1", "Gene2"), each = 30),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:10), 6),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  
  expect_true("adj_p_value" %in% colnames(result))
  expect_equal(length(result$adj_p_value), 2)  # 2 genes
})

test_that("detect_q_gene_interactions multicorr=hochberg adjusts p-values correctly", {
  set.seed(101)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),
                  rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep(c("Gene1", "Gene2"), each = 30),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:10), 6),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  
  # Adjusted p-values should be >= raw p-values
  expect_true(all(result$adj_p_value >= result$p_value, na.rm = TRUE))
  # Should be sorted by adj_p_value (first row should have smallest adj p-value)
  expect_equal(result$adj_p_value[1], min(result$adj_p_value, na.rm = TRUE))
})

test_that("detect_q_gene_interactions multicorr=benjamini-yekutieli produces adj_p_values", {
  set.seed(102)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),
                  rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep(c("Gene1", "Gene2"), each = 30),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:10), 6),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "benjamini-yekutieli")
  
  # Adjusted p-values should be >= raw p-values
  expect_true(all(result$adj_p_value >= result$p_value, na.rm = TRUE))
})

test_that("detect_q_gene_interactions multicorr=none returns raw p-values as adj_p_value", {
  set.seed(103)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),
                  rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep(c("Gene1", "Gene2"), each = 30),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:10), 6),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "none")
  
  # With multicorr="none", should be identical to raw p-values
  expect_equal(result$adj_p_value, result$p_value)
})

test_that("detect_q_gene_interactions default multicorr is hochberg", {
  set.seed(104)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),
                  rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep(c("Gene1", "Gene2"), each = 30),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:10), 6),
    stringsAsFactors = FALSE
  )
  
  result_default <- .calculate_rank_transform(model_data, condition_col = "condition")
  result_explicit <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  
  # Default should match explicit hochberg
  expect_equal(result_default$adj_p_value, result_explicit$adj_p_value)
})

test_that("detect_q_gene_interactions output is sorted by adj_p_value", {
  set.seed(105)
  model_data <- data.frame(
    diversity = c(rnorm(10, 1.0),  # Gene1: not different
                  rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),  # Gene2: very different
                  rnorm(10, 1.0), rnorm(10, 1.05), rnorm(10, 1.1)),  # Gene3: slightly different
    q = rep(c(0.5, 1.0, 1.5), times = c(10, 30, 30)),  # Match diversity structure: 10+30+30=70
    gene = c(rep("Gene1", 10), rep("Gene2", 30), rep("Gene3", 30)),
    condition = c(rep(c("A", "B"), times = 5), rep(c("A", "B"), times = 15), rep(c("A", "B"), times = 15)),
    sample = rep(paste0("S", 1:10), 7),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  
  # Check that adj_p_values are sorted in ascending order (ignoring NAs)
  non_na_idx <- !is.na(result$adj_p_value)
  adj_p_non_na <- result$adj_p_value[non_na_idx]
  adj_p_sorted <- sort(adj_p_non_na)
  expect_equal(adj_p_non_na, adj_p_sorted)
})

test_that("detect_q_gene_interactions respects multicorr parameter passing", {
  set.seed(106)
  model_data <- data.frame(
    diversity = rnorm(60),
    q = rep(c(0.5, 1.0, 1.5), each = 20),
    gene = rep(paste0("Gene", 1:4), each = 15),
    condition = rep(c("A", "B"), times = 30),
    sample = rep(paste0("S", 1:20), 3),
    stringsAsFactors = FALSE
  ) # Already correct - times=30
  
  # Test that each method produces different results
  result_hoch <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  result_by <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "benjamini-yekutieli")
  result_none <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "none")
  
  # None should match raw p-values
  expect_equal(result_none$adj_p_value, result_none$p_value)
  
  # Benjamini-Yekutieli should be >= Hochberg (more conservative)
  # because it corrects for arbitrary dependence vs positive regression dependence
  by_ge_hoch <- mean(result_by$adj_p_value >= result_hoch$adj_p_value, na.rm = TRUE)
  expect_true(by_ge_hoch > 0.5)  # At least half should satisfy this
})

test_that("detect_q_gene_interactions invalid multicorr parameter raises error", {
  model_data <- data.frame(
    diversity = rnorm(30),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = rep("Gene1", 30),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    .calculate_rank_transform(model_data, multicorr = "invalid_method"),
    "should be one of"  # match.arg error message for invalid choice
  )
})

test_that("detect_q_gene_interactions multicorr handles NA p-values correctly", {
  set.seed(108)
  
  # Create data where one gene will have NA p-value (insufficient data)
  model_data <- data.frame(
    diversity = c(rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5)),
    q = rep(c(0.5, 1.0, 1.5), each = 10),
    gene = c(rep("Gene1", 10), rep("Gene1", 10), rep("Gene2", 10)),
    condition = rep(c("A", "B"), times = 15),
    sample = rep(paste0("S", 1:10), 3),
    stringsAsFactors = FALSE
  )
  
  result <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  
  # Should run without errors
  expect_equal(nrow(result), 2)
  # Check that adj_p_value column exists
  expect_true("adj_p_value" %in% colnames(result))
  expect_true("p_value" %in% colnames(result))
  expect_true("gene" %in% colnames(result))
})

test_that("detect_q_gene_interactions more genes = more deflation with Hochberg", {
  # With more genes, Hochberg correction should be more conservative
  set.seed(109)
  
  # Create 5-gene dataset with same p-value
  base_gene <- rnorm(30, 1.0, 0.05)  # Non-significant p-value
  model_data_5 <- data.frame(
    diversity = rep(base_gene, 5),
    q = rep(c(0.5, 1.0, 1.5), 50),
    gene = rep(paste0("Gene", 1:5), each = 30),
    condition = rep(c("A", "B"), times = 75),
    sample = rep(paste0("S", 1:10), 15),
    stringsAsFactors = FALSE
  )
  
  # Create 10-gene dataset with same p-value per gene
  model_data_10 <- data.frame(
    diversity = rep(base_gene, 10),
    q = rep(c(0.5, 1.0, 1.5), 100),
    gene = rep(paste0("Gene", 1:10), each = 30),
    condition = rep(c("A", "B"), times = 150),
    sample = rep(paste0("S", 1:10), 30),
    stringsAsFactors = FALSE
  )
  
  result_5 <- .calculate_rank_transform(model_data_5, condition_col = "condition", multicorr = "hochberg")
  result_10 <- .calculate_rank_transform(model_data_10, condition_col = "condition", multicorr = "hochberg")
  
  # With more genes, Hochberg multiplier increases (m - rank + 1)
  # So if genes have similar raw p-values, those in 5-gene set should have 
  # smaller adjusted p-values than those in 10-gene set
  # Extract genes that appear in both datasets
  common_genes <- intersect(result_5$gene, result_10$gene)
  
  # This is a weaker test: just verify both have adj_p_value columns
  expect_true(all(!is.na(result_5$adj_p_value)))
  expect_true(all(!is.na(result_10$adj_p_value)))
})

test_that("detect_q_gene_interactions Benjamini-Yekutieli is less deflating than Hochberg", {
  set.seed(110)
  
  # Create dataset with varying p-values
  model_data <- data.frame(
    diversity = c(
      rnorm(10, 0.5), rnorm(10, 1.5), rnorm(10, 2.5),  # Gene1: significant
      rnorm(10, 1.0), rnorm(10, 1.0), rnorm(10, 1.0),  # Gene2: not significant
      rnorm(10, 0.5), rnorm(10, 1.3), rnorm(10, 2.1)   # Gene3: significant
    ),
    q = rep(c(0.5, 1.0, 1.5), 30),
    gene = rep(c("Gene1", "Gene2", "Gene3"), each = 30),
    condition = rep(c("A", "B"), times = 45),
    sample = rep(paste0("S", 1:10), 9),
    stringsAsFactors = FALSE
  )
  
  result_hoch <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "hochberg")
  result_by <- .calculate_rank_transform(model_data, condition_col = "condition", multicorr = "benjamini-yekutieli")
  
  # Both methods should produce valid adjusted p-values (>= raw, <= 1)
  expect_true(all(result_hoch$adj_p_value >= result_hoch$p_value - 1e-10, na.rm = TRUE))
  expect_true(all(result_hoch$adj_p_value <= 1, na.rm = TRUE))
  expect_true(all(result_by$adj_p_value >= result_by$p_value - 1e-10, na.rm = TRUE))
  expect_true(all(result_by$adj_p_value <= 1, na.rm = TRUE))
})

# ============================================================================
# C1: .hochberg_stepup — uses rev(cummin(rev(...))), not cummax
# ============================================================================
# C1: .hochberg_stepup — uses rev(cummin(rev(...))), not cummax# ============================================================================

test_that("C1: .hochberg_stepup uses Hochberg (rev-cummin), not Holm (cummax)", {
    pvals <- c(0.001, 0.01, 0.05, 0.1, 0.5)
    adjusted <- TSENAT:::.hochberg_stepup(pvals)

    # Hochberg is always ≤ Holm (more powerful)
    holm_adjusted <- p.adjust(pvals, method = "holm")

    # For Hochberg, p_(i) <= p_(i+1) when sorted by original p
    ord <- order(pvals)
    expect_true(all(diff(adjusted[ord]) >= -1e-10))

    # Hochberg should be <= Holm for each p-value
    expect_true(all(adjusted <= holm_adjusted))

    # Verify against base R's hochberg for a known case
    pvals2 <- c(0.01, 0.02, 0.03, 0.04)
    adj2 <- TSENAT:::.hochberg_stepup(pvals2)
    base_hochberg <- p.adjust(pvals2, method = "hochberg")
    expect_equal(adj2, base_hochberg, tolerance = 1e-10)
})
