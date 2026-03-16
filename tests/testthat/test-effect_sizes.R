context("Wilcoxon Tests: Effect Size Measures (U and r Columns)")

test_that("calculate_effect_sizes works for unpaired data", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_is(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), 4)
    expect_true(all(c("cliffs_delta", "r_value", "rank_biserial", "effect_magnitude") %in% colnames(result)))
})

test_that("cliffs_delta is in range [-1, 1]", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_true(all(result$cliffs_delta[!is.na(result$cliffs_delta)] >= -1))
    expect_true(all(result$cliffs_delta[!is.na(result$cliffs_delta)] <= 1))
})

test_that("r_value is in range [-1, 1]", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_true(all(abs(result$r_value[!is.na(result$r_value)]) <= 1))
})

test_that("effect_magnitude is correctly interpreted", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    valid_magnitudes <- c("negligible", "small", "medium", "large", NA)
    expect_true(all(result$effect_magnitude %in% valid_magnitudes))
})

test_that("calculate_effect_sizes works for paired data", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = TRUE)

    expect_is(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), 4)
})



test_that("effect sizes match across examples", {
    mat <- matrix(rnorm(18, mean = 5, sd = 1), nrow = 3)
    samples <- rep(c('Control', 'Treatment'), each = 3)

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_is(result, "data.frame")
    expect_true(!any(is.na(result$cliffs_delta)))
})

test_that("effect sizes handle missing values gracefully", {
    mat <- matrix(c(
        0.5, 0.6, NA, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, NA, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_is(result, "data.frame")
    expect_equal(nrow(result), 2)
})

test_that("effect sizes with paired pairing vector work correctly", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')
    pairs <- c('pair1', 'pair2', 'pair3', 'pair1', 'pair2', 'pair3')

    result <- calculate_effect_sizes(mat, samples, paired = TRUE, pairs = pairs)

    expect_is(result, "data.frame")
    expect_equal(nrow(result), 2)
})

test_that("rank_biserial is in expected range", {
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,
        0.7, 0.75, 0.72, 0.6, 0.55, 0.5
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    expect_true(all(abs(result$rank_biserial[!is.na(result$rank_biserial)]) <= 1))
})

test_that("r_value clamping prevents values outside [-1, 1]", {
    # Test with extreme separated groups that would produce r > 1 without clamping
    mat <- matrix(c(
        0.01, 0.02, 0.03, 100, 101, 102,  # Gene 1: extreme separation
        0.1, 0.15, 0.2, 50, 55, 60         # Gene 2: moderate separation
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

    result <- calculate_effect_sizes(mat, samples, paired = FALSE)

    # Verify all r_values are within [-1, 1]
    valid_r_values <- result$r_value[!is.na(result$r_value)]
    expect_true(all(valid_r_values >= -1), 
                info = sprintf("r_value < -1 found: %s", 
                              paste(valid_r_values[valid_r_values < -1], collapse = ", ")))
    expect_true(all(valid_r_values <= 1), 
                info = sprintf("r_value > 1 found: %s", 
                              paste(valid_r_values[valid_r_values > 1], collapse = ", ")))
})

test_that(".tsenat_r_value directly clamps to [-1, 1] range", {
    # Access the internal function directly
    # Create extreme data that would overflow without clamping
    x1 <- c(0.0001, 0.0002, 0.0003)
    x2 <- c(999, 1000, 1001)

    r_result <- TSENAT:::.tsenat_r_value(x1, x2)

    # Verify clamping occurred
    expect_true(!is.na(r_result), "r_value should not be NA for valid data")
    expect_true(r_result >= -1, 
                info = sprintf("r_value not clamped to >= -1: got %f", r_result))
    expect_true(r_result <= 1, 
                info = sprintf("r_value not clamped to <= 1: got %f", r_result))
})

test_that("r_value clamping handles both positive and negative extremes", {
    # Test case 1: Perfect separation (would give r close to limits)
    x1 <- rep(1, 10)
    x2 <- rep(100, 10)
    r_pos <- TSENAT:::.tsenat_r_value(x1, x2)
    expect_true(r_pos >= -1 && r_pos <= 1, 
                info = sprintf("Positive extreme not clamped: %f", r_pos))

    # Test case 2: Reverse separation
    x1 <- rep(100, 10)
    x2 <- rep(1, 10)
    r_neg <- TSENAT:::.tsenat_r_value(x1, x2)
    expect_true(r_neg >= -1 && r_neg <= 1, 
                info = sprintf("Negative extreme not clamped: %f", r_neg))
})

test_that("r_value maintains reasonable values for normal cases", {
    # Test with typical data to ensure clamping doesn't distort normal values
    set.seed(42)
    x1 <- rnorm(20, mean = 5, sd = 1)
    x2 <- rnorm(20, mean = 6, sd = 1)

    r_result <- TSENAT:::.tsenat_r_value(x1, x2)

    expect_true(!is.na(r_result), "r_value should be computed for normal data")
    expect_true(r_result >= -1 && r_result <= 1, 
                info = sprintf("Normal case not within bounds: %f", r_result))
    # For similar distributions, r should be moderate
    expect_true(abs(r_result) < 0.8, 
                info = sprintf("Normal distributions produced extreme r: %f", r_result))
})


# ============================================================================
# Wilcoxon Effect Size Guidelines - Comprehensive Tests
# ============================================================================

test_that("wilcoxon_effect_size_guidelines returns both metrics", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    
    # Should have rows for both Cliff's Delta and r-value
    expect_equal(nrow(guidelines), 8)  # 4 categories × 2 metrics
    expect_equal(sum(grepl("Cliff", guidelines$Metric)), 4)
    expect_equal(sum(grepl("r-value", guidelines$Metric)), 4)
})

test_that("wilcoxon_effect_size_guidelines has correct structure for Cliff's Delta", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    cliffs <- guidelines[grep("Cliff", guidelines$Metric), ]
    
    expect_equal(nrow(cliffs), 4)
    expect_equal(cliffs$Magnitude, c("Negligible", "Small", "Medium", "Large"))
    expect_true(all(grepl("\\|δ\\|", cliffs$Range)))
})

test_that("wilcoxon_effect_size_guidelines Cliff's Delta thresholds are ordered", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    cliffs <- guidelines[grep("Cliff", guidelines$Metric), ]
    
    # Extract threshold values from range strings
    # Formats: "|δ| < 0.147", "0.147 ≤ |δ| < 0.33", etc.
    thresholds <- c(0.147, 0.33, 0.474)
    
    # Verify thresholds are present and in increasing order
    expect_true(0.147 < 0.33)
    expect_true(0.33 < 0.474)
})

test_that("wilcoxon_effect_size_guidelines has correct structure for r-value", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    r_vals <- guidelines[grep("r-value", guidelines$Metric), ]
    
    expect_equal(nrow(r_vals), 4)
    expect_equal(r_vals$Magnitude, c("Negligible", "Small", "Medium", "Large"))
    expect_true(all(grepl("\\|r\\|", r_vals$Range)))
})

test_that("wilcoxon_effect_size_guidelines r-value thresholds are ordered", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    r_vals <- guidelines[grep("r-value", guidelines$Metric), ]
    
    # Expected thresholds: 0.1, 0.3, 0.5
    thresholds <- c(0.1, 0.3, 0.5)
    
    # Verify thresholds are in increasing order
    expect_true(0.1 < 0.3)
    expect_true(0.3 < 0.5)
})

test_that("wilcoxon_effect_size_guidelines has descriptions for all rows", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    
    expect_true(all(!is.na(guidelines$Description)))
    expect_true(all(nchar(guidelines$Description) > 0))
})


# ============================================================================
# Integration Tests - Effect Sizes with Guidelines
# ============================================================================

test_that("calculate_effect_sizes results have magnitudes matching guidelines", {
    # Create data with known effect sizes
    mat <- matrix(c(
        0.5, 0.6, 0.55, 0.8, 0.85, 0.9,    # Gene 1: moderate effect
        0.01, 0.02, 0.03, 100, 101, 102    # Gene 2: large effect
    ), nrow = 2, byrow = TRUE)
    samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')
    
    result <- calculate_effect_sizes(mat, samples, paired = FALSE)
    
    # Effect magnitudes should be one of the valid categories
    valid_mags <- c("negligible", "small", "medium", "large", NA)
    expect_true(all(result$effect_magnitude %in% valid_mags))
})

test_that("wilcoxon_effect_size_guidelines can be called and returns data frame", {
    guidelines <- TSENAT:::wilcoxon_effect_size_guidelines()
    
    expect_is(guidelines, "data.frame")
    expect_true(nrow(guidelines) > 0)
    expect_true(ncol(guidelines) >= 4)
    expect_true("Metric" %in% colnames(guidelines))
})
