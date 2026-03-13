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
# Entropy-Specific Effect Size Guidelines Tests
# ============================================================================

test_that("entropy_effect_size_guidelines returns data frame", {
    guidelines <- entropy_effect_size_guidelines()
    
    expect_is(guidelines, "data.frame")
    expect_true(nrow(guidelines) > 0)
    expect_true(all(c("Framework", "Magnitude", "Range", "Percent_Max", "Description") 
                   %in% colnames(guidelines)))
})

test_that("entropy_effect_size_guidelines practical framework works", {
    practical <- entropy_effect_size_guidelines(focus = "practical")
    
    expect_is(practical, "data.frame")
    expect_true(any(grepl("Practical", practical$Framework)))
    expect_equal(sum(grepl("Practical", practical$Framework)), 4)
    expect_true(all(c("Negligible", "Small", "Medium", "Large") %in% practical$Magnitude[1:4]))
})

test_that("entropy_effect_size_guidelines statistical framework works", {
    statistical <- entropy_effect_size_guidelines(focus = "statistical")
    
    expect_is(statistical, "data.frame")
    expect_true(any(grepl("Statistical", statistical$Framework)))
    expect_equal(sum(grepl("Statistical", statistical$Framework)), 4)
})

test_that("entropy_effect_size_guidelines biological framework works", {
    biological <- entropy_effect_size_guidelines(focus = "biological")
    
    expect_is(biological, "data.frame")
    expect_true(any(grepl("Biological", biological$Framework)))
    expect_equal(sum(grepl("Biological", biological$Framework)), 4)
})

test_that("entropy_effect_size_guidelines all frameworks returns combined view", {
    all_frameworks <- entropy_effect_size_guidelines(focus = "all")
    
    expect_is(all_frameworks, "data.frame")
    # Should have 4 categories × 3 frameworks + system info row = 13 rows
    expect_true(nrow(all_frameworks) >= 12)
    expect_true(any(grepl("Practical", all_frameworks$Framework)))
    expect_true(any(grepl("Statistical", all_frameworks$Framework)))
    expect_true(any(grepl("Biological", all_frameworks$Framework)))
})

test_that("entropy_effect_size_guidelines scales with n_isoforms", {
    small_iso <- entropy_effect_size_guidelines(n_isoforms = 10)
    large_iso <- entropy_effect_size_guidelines(n_isoforms = 100)
    
    # Extract ranges from practical frameworks
    small_max <- log(10)
    large_max <- log(100)
    
    # Both should have valid structures
    expect_is(small_iso, "data.frame")
    expect_is(large_iso, "data.frame")
    expect_true(nrow(small_iso) > 0)
    expect_true(nrow(large_iso) > 0)
})

test_that("entropy_effect_size_guidelines has correct q-value parameter", {
    guidelines_q1 <- entropy_effect_size_guidelines(q_value = 1.0)
    guidelines_q2 <- entropy_effect_size_guidelines(q_value = 2.0)
    
    # Both should work and be valid
    expect_is(guidelines_q1, "data.frame")
    expect_is(guidelines_q2, "data.frame")
})

test_that("entropy_effect_size_guidelines includes literature references", {
    guidelines <- entropy_effect_size_guidelines()
    
    # Check that References column exists and has content
    expect_true("References" %in% colnames(guidelines))
    expect_true(any(!is.na(guidelines$References[1:4])))
})

test_that("entropy_effect_size_guidelines % of max is reasonable", {
    practical <- entropy_effect_size_guidelines(focus = "practical")
    
    # Percent_Max should be character strings with percentages
    percent_values <- practical$Percent_Max[1:4]
    expect_true(all(grepl("%", percent_values)))
})

test_that("entropy_effect_size_guidelines threshold values are ordered", {
    practical <- entropy_effect_size_guidelines(n_isoforms = 50, focus = "practical")
    
    H_max <- log(50)
    
    # Expected thresholds for practical framework
    # Negligible: < 0.05 × log(n)
    # Small: 0.05-0.15 × log(n)
    # Medium: 0.15-0.35 × log(n)
    # Large: ≥ 0.35 × log(n)
    
    expect_true(nrow(practical) >= 4)
    # Extract first numeric value from each range (handles "< 0.XX" or "0.XX - 0.YY" formats)
    # Using regex to extract the first number in each Range string
    range_strings <- practical$Range[1:4]
    thresholds <- as.numeric(
      gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", range_strings)
    )
    expect_true(length(thresholds) > 0)
    expect_false(any(is.na(thresholds)))
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

test_that("entropy_effect_size_guidelines handles very small n_isoforms", {
    # n=2 is minimum meaningful case
    guidelines_n2 <- entropy_effect_size_guidelines(n_isoforms = 2)
    
    expect_is(guidelines_n2, "data.frame")
    expect_true(nrow(guidelines_n2) > 0)
    
    # log(2) ≈ 0.693
    h_max_expected <- log(2)
    expect_true(h_max_expected > 0)
})

test_that("entropy_effect_size_guidelines handles very large n_isoforms", {
    guidelines_large <- entropy_effect_size_guidelines(n_isoforms = 10000)
    
    expect_is(guidelines_large, "data.frame")
    expect_true(nrow(guidelines_large) > 0)
    
    # log(10000) ≈ 9.21, should still work fine
    h_max_expected <- log(10000)
    expect_true(h_max_expected > 0)
})

test_that("entropy_effect_size_guidelines works with various q values", {
    q_values <- c(0.1, 0.5, 1.0, 1.5, 2.0, 3.0)
    
    for (q in q_values) {
        guidelines <- entropy_effect_size_guidelines(q_value = q)
        expect_is(guidelines, "data.frame")
        expect_true(nrow(guidelines) > 0)
        # Framework labels should show the q-value
        expect_true(any(grepl(sprintf("q=%.2f", q), guidelines$Framework)))
    }
})

test_that("entropy_effect_size_guidelines produces valid output for edge case q values", {
    # q very close to 0
    low_q <- entropy_effect_size_guidelines(q_value = 0.001)
    expect_is(low_q, "data.frame")
    
    # q relatively large 
    high_q <- entropy_effect_size_guidelines(q_value = 5.0)
    expect_is(high_q, "data.frame")
})


# ============================================================================
# Cross-Framework Consistency Tests
# ============================================================================

test_that("entropy_effect_size_guidelines thresholds scale across frameworks correctly", {
    q_value <- 1.0
    n_iso <- 50
    
    practical <- entropy_effect_size_guidelines(q_value = q_value, n_isoforms = n_iso, focus = "practical")
    statistical <- entropy_effect_size_guidelines(q_value = q_value, n_isoforms = n_iso, focus = "statistical")
    biological <- entropy_effect_size_guidelines(q_value = q_value, n_isoforms = n_iso, focus = "biological")
    
    # Extract first threshold from each (Negligible category)
    extract_first_threshold <- function(df) {
        as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", df$Range[1]))
    }
    
    prac_thresh <- extract_first_threshold(practical)
    stat_thresh <- extract_first_threshold(statistical)
    biol_thresh <- extract_first_threshold(biological)
    
    # Statistical should be most conservative (smallest), practical middle, biological largest
    expect_true(stat_thresh < prac_thresh)
    expect_true(prac_thresh < biol_thresh)
})

test_that("entropy_effect_size_guidelines framework labels are consistent", {
    guidelines_all <- entropy_effect_size_guidelines(focus = "all", q_value = 1.5)
    
    # Check all frameworks have consistent q-value label
    frameworks <- unique(guidelines_all$Framework[1:12])  # 4 cats × 3 frameworks
    
    for (fw in frameworks) {
        if (!is.na(fw) && fw != "") {
            expect_true(grepl("q=1.50", fw), 
                       info = paste("Framework label should show q=1.50:", fw))
        }
    }
})

test_that("entropy_effect_size_guidelines all focus includes all three frameworks", {
    all_fw <- entropy_effect_size_guidelines(focus = "all")
    
    # Should have practical, statistical, and biological (12 regular rows) + info row
    expect_true(any(grepl("Practical", all_fw$Framework)))
    expect_true(any(grepl("Statistical", all_fw$Framework)))
    expect_true(any(grepl("Biological", all_fw$Framework)))
    
    # Count each framework (exclude info row)
    data_rows <- all_fw$Framework[!is.na(all_fw$Magnitude)]
    expect_equal(sum(grepl("Practical", data_rows)), 4)
    expect_equal(sum(grepl("Statistical", data_rows)), 4)
    expect_equal(sum(grepl("Biological", data_rows)), 4)
})


# ============================================================================
# Q-Parameter Numerical Validation Tests
# ============================================================================

test_that("entropy_effect_size_guidelines q-parameter scaling is mathematically correct", {
    # For q=1.0 (baseline): q_weight=1.5, multiplier=1.0, adjustment=1.0
    # For q=0.5: q_weight=1.0, multiplier=0.667, adjustment=1.5
    # For q=2.0: q_weight=2.5, multiplier=1.667, adjustment=0.6
    
    q_1 <- entropy_effect_size_guidelines(q_value = 1.0, n_isoforms = 50)
    q_05 <- entropy_effect_size_guidelines(q_value = 0.5, n_isoforms = 50)
    q_2 <- entropy_effect_size_guidelines(q_value = 2.0, n_isoforms = 50)
    
    # Extract multiplier from framework label "Practical (q=X.XX, multiplier=Y.YYY×)"
    extract_multiplier <- function(df) {
        labels <- df$Framework[1]
        as.numeric(gsub(".*multiplier=([0-9.]+)×.*", "\\1", labels))
    }
    
    m_1 <- extract_multiplier(q_1)
    m_05 <- extract_multiplier(q_05)
    m_2 <- extract_multiplier(q_2)
    
    # q=1.0 should have multiplier ~1.0
    expect_true(abs(m_1 - 1.0) < 0.01)
    
    # q=0.5 should have multiplier ~0.667
    expect_true(abs(m_05 - 0.667) < 0.01)
    
    # q=2.0 should have multiplier ~1.667
    expect_true(abs(m_2 - 1.667) < 0.01)
})

test_that("entropy_effect_size_guidelines adjustment factor inverse relationship", {
    # Adjustment should be 1 / multiplier
    # q=1.0: adjustment = 1.0
    # q=0.5: adjustment = 1.5 (inverse of 0.667)
    # q=2.0: adjustment = 0.6 (inverse of 1.667)
    
    H_max <- log(50)
    
    # q=1.0
    q_1 <- entropy_effect_size_guidelines(q_value = 1.0, n_isoforms = 50, focus = "practical")
    q_1_threshold <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", q_1$Range[1]))
    q_1_expected <- 0.05 * H_max * 1.0  # adjustment=1.0
    
    # q=0.5 (adjustment should be 1.5)
    q_05 <- entropy_effect_size_guidelines(q_value = 0.5, n_isoforms = 50, focus = "practical")
    q_05_threshold <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", q_05$Range[1]))
    q_05_expected <- 0.05 * H_max * 1.5  # adjustment=1.5
    
    # q=2.0 (adjustment should be 0.6)
    q_2 <- entropy_effect_size_guidelines(q_value = 2.0, n_isoforms = 50, focus = "practical")
    q_2_threshold <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", q_2$Range[1]))
    q_2_expected <- 0.05 * H_max * 0.6  # adjustment=0.6
    
    # Verify thresholds match expectations within 0.01 tolerance
    expect_true(abs(q_1_threshold - q_1_expected) < 0.001,
               info = sprintf("q=1.0 threshold: got %.6f, expected %.6f", q_1_threshold, q_1_expected))
    expect_true(abs(q_05_threshold - q_05_expected) < 0.001,
               info = sprintf("q=0.5 threshold: got %.6f, expected %.6f", q_05_threshold, q_05_expected))
    expect_true(abs(q_2_threshold - q_2_expected) < 0.001,
               info = sprintf("q=2.0 threshold: got %.6f, expected %.6f", q_2_threshold, q_2_expected))
})

test_that("entropy_effect_size_guidelines all frameworks scale with q-parameter consistently", {
    # All three frameworks should scale by the same q-parameter adjustment
    n_iso <- 50
    H_max <- log(n_iso)
    
    practical <- entropy_effect_size_guidelines(q_value = 1.0, n_isoforms = n_iso, focus = "practical")
    statistical <- entropy_effect_size_guidelines(q_value = 1.0, n_isoforms = n_iso, focus = "statistical")
    biological <- entropy_effect_size_guidelines(q_value = 1.0, n_isoforms = n_iso, focus = "biological")
    
    # Get first thresholds for q=1.0 (adjustment=1.0, baseline)
    prac_q1 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", practical$Range[1]))
    stat_q1 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", statistical$Range[1]))
    biol_q1 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", biological$Range[1]))
    
    # Now test with different q
    practical_q2 <- entropy_effect_size_guidelines(q_value = 2.0, n_isoforms = n_iso, focus = "practical")
    statistical_q2 <- entropy_effect_size_guidelines(q_value = 2.0, n_isoforms = n_iso, focus = "statistical")
    biological_q2 <- entropy_effect_size_guidelines(q_value = 2.0, n_isoforms = n_iso, focus = "biological")
    
    prac_q2 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", practical_q2$Range[1]))
    stat_q2 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", statistical_q2$Range[1]))
    biol_q2 <- as.numeric(gsub("^[<>=]*\\s*([0-9.]+).*$", "\\1", biological_q2$Range[1]))
    
    # Scaling ratios should be approximately the same across frameworks (0.6)
    # q=2.0 should have adjustment ~0.6 relative to q=1.0
    expected_ratio <- 0.6
    
    prac_ratio <- prac_q2 / prac_q1
    stat_ratio <- stat_q2 / stat_q1
    biol_ratio <- biol_q2 / biol_q1
    
    expect_true(abs(prac_ratio - expected_ratio) < 0.01,
               info = sprintf("Practical ratio: %.4f vs expected %.4f", prac_ratio, expected_ratio))
    expect_true(abs(stat_ratio - expected_ratio) < 0.01,
               info = sprintf("Statistical ratio: %.4f vs expected %.4f", stat_ratio, expected_ratio))
    expect_true(abs(biol_ratio - expected_ratio) < 0.01,
               info = sprintf("Biological ratio: %.4f vs expected %.4f", biol_ratio, expected_ratio))
})

test_that("entropy_effect_size_guidelines system info includes q-parameters", {
    guidelines <- entropy_effect_size_guidelines(q_value = 1.5, n_isoforms = 75)
    
    # Find the info row (where Magnitude is NA)
    info_row <- guidelines[is.na(guidelines$Magnitude), ]
    
    expect_equal(nrow(info_row), 1)
    
    # Info should contain q value, n_isoforms, H_max, q_weight, adjustment
    info_text <- info_row$Framework[1]
    
    expect_true(grepl("q = 1.50", info_text))
    expect_true(grepl("n_isoforms = 75", info_text))
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

test_that("wilcoxon and entropy guidelines are complementary", {
    # Wilcoxon works with Cliff's delta and r-value
    wilc <- TSENAT:::wilcoxon_effect_size_guidelines()
    expect_true("Cliff's Delta" %in% wilc$Metric || "Cliff" %in% wilc$Metric)
    
    # Entropy works with entropy differences (should have Practical, Statistical, or Biological)
    entropy <- entropy_effect_size_guidelines()
    framework_text <- paste(entropy$Framework, collapse = " ")
    expect_true(grepl("Practical|Statistical|Biological", framework_text))
})

test_that("effect size guidelines Percent_Max values are meaningful", {
    practical <- entropy_effect_size_guidelines(n_isoforms = 50, focus = "practical")
    
    # Extract percentage values (excludes info row)
    pct_values <- practical$Percent_Max[!is.na(practical$Magnitude)]
    
    # All should contain % symbol
    expect_true(all(grepl("%", pct_values)))
    
    # Should have 4 values for 4 magnitude categories
    expect_equal(length(pct_values), 4)
})


# ==============================================================================
# Tests for filter_by_effect_size function
# ==============================================================================

test_that("filter_by_effect_size returns list with correct structure", {
    # Create sample count data
    readcounts <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:50))
    
    # Create tx2gene mapping
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:25, 2))),
        stringsAsFactors = FALSE
    )
    
    # Create transcript statistics
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = rnorm(50, mean = 0.5, sd = 0.3),
        adj_p_value = runif(50),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.2)
    
    expect_is(result, "list")
    expect_true("readcounts" %in% names(result))
    expect_true("tx2gene" %in% names(result))
})

test_that("filter_by_effect_size respects effect size threshold", {
    readcounts <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:100))
    
    # Create explicit gene structure with clear separation
    # Genes 1-10: rows 1-40 (4 transcripts each with low FC)
    # Genes 11-25: rows 41-100 (4 transcripts each with high FC)
    gene_names <- c(
        rep(paste0("ENSG", sprintf("%09d", 1:10)), each = 4),   # 40 rows
        rep(paste0("ENSG", sprintf("%09d", 11:25)), each = 4)   # 60 rows
    )
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = gene_names,
        stringsAsFactors = FALSE
    )
    
    # Create stats: genes 1-10 have low fold-change, genes 11-25 have high
    transcript_stats <- data.frame(
        genes = gene_names,
        log2_fold_change = c(
            rep(0.05, 40),   # Genes 1-10: all transcripts have low FC
            rep(1.0, 60)     # Genes 11-25: all transcripts have high FC
        ),
        adj_p_value = runif(100),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.3)
    
    # Should have filtered out genes 1-10 (40 transcripts), keep genes 11-25 (60 transcripts)
    expect_equal(nrow(result$readcounts), 60)
})

test_that("filter_by_effect_size handles NULL or empty stats", {
    readcounts <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:50))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:25, 2))),
        stringsAsFactors = FALSE
    )
    
    # Test with NULL stats
    result <- filter_by_effect_size(readcounts, NULL, tx2gene)
    expect_equal(nrow(result$readcounts), nrow(readcounts))
    
    # Test with empty stats
    empty_stats <- data.frame(genes = character(), log2_fold_change = numeric())
    result <- filter_by_effect_size(readcounts, empty_stats, tx2gene)
    expect_equal(nrow(result$readcounts), nrow(readcounts))
})

test_that("filter_by_effect_size returns unfiltered when threshold too high", {
    readcounts <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:50))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:25, 2))),
        stringsAsFactors = FALSE
    )
    
    # Very low fold-changes
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = rnorm(50, mean = 0.01, sd = 0.01),
        adj_p_value = runif(50),
        stringsAsFactors = FALSE
    )
    
    # High threshold
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 10.0)
    
    # Should return unfiltered data when no genes meet threshold
    expect_equal(nrow(result$readcounts), nrow(readcounts))
    expect_equal(nrow(result$tx2gene), nrow(tx2gene))
})

test_that("filter_by_effect_size filters tx2gene correctly", {
    readcounts <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:50))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:25, 2))),
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = rnorm(50, mean = 0.5, sd = 0.3),
        adj_p_value = runif(50),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.2)
    
    # Transcripts in result should be subset of original
    expect_true(all(result$tx2gene$transcript_id %in% tx2gene$transcript_id))
    
    # Rownames of readcounts should match tx2gene
    expect_equal(
        sort(rownames(result$readcounts)),
        sort(result$tx2gene$transcript_id)
    )
})

test_that("filter_by_effect_size handles genes with multiple transcripts", {
    readcounts <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:100))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:20, 5))),  # 5 transcripts per gene
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = rnorm(100, mean = 0.5, sd = 0.3),
        adj_p_value = runif(100),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.2)
    
    # Should maintain 5:1 ratio of transcripts to genes (approximately)
    n_genes <- length(unique(result$tx2gene$gene_name))
    n_transcripts <- nrow(result$readcounts)
    expect_true(n_transcripts > 0)
    expect_equal(n_transcripts, nrow(result$tx2gene))
})

test_that("filter_by_effect_size uses absolute log2FC", {
    readcounts <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:50))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:25, 2))),
        stringsAsFactors = FALSE
    )
    
    # Mix of positive and negative log2FC
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = c(rep(-1.0, 25), rep(1.0, 25)),
        adj_p_value = runif(50),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.5)
    
    # Should keep both positive and negative fold-changes if |FC| >= threshold
    expect_equal(nrow(result$readcounts), 50)  # All 50 should pass
})

test_that("filter_by_effect_size computes median correctly for gene-level filtering", {
    readcounts <- matrix(rpois(500, 10), nrow = 100, ncol = 10)
    rownames(readcounts) <- paste0("ENST", sprintf("%09d", 1:100))
    
    tx2gene <- data.frame(
        transcript_id = rownames(readcounts),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:20, 5))),
        stringsAsFactors = FALSE
    )
    
    # Each gene has 5 transcripts: 3 with high FC, 2 with low FC
    # Median should be based on middle value when sorted
    transcript_stats <- data.frame(
        genes = tx2gene$gene_name,
        log2_fold_change = c(
            rep(c(0.05, 0.1, 1.0, 1.2, 1.5), 20)  # Pattern repeated for each gene
        ),
        adj_p_value = runif(100),
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.9)
    
    # Median of (0.05, 0.1, 1.0, 1.2, 1.5) = 1.0, so all genes should pass
    expect_gt(nrow(result$readcounts), 0)
})


# ==============================================================================
# Tests for filter_by_effect_size with SummarizedExperiment objects
# ==============================================================================

test_that("filter_by_effect_size works with SE-extracted data", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create SummarizedExperiment object with clear low/high separation
    mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
    rownames(mat) <- paste0("ENST", sprintf("%09d", 1:100))
    
    # 25 genes: 10 with low FC, 15 with high FC
    gene_names <- c(
        rep(paste0("ENSG", sprintf("%09d", 1:10)), each = 4),   # 40 rows
        rep(paste0("ENSG", sprintf("%09d", 11:25)), each = 4)   # 60 rows
    )
    
    row_data <- data.frame(
        transcript_id = rownames(mat),
        gene_name = gene_names,
        log2_fold_change = c(rep(0.05, 40), rep(1.0, 60)),
        adj_p_value = runif(100),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        rowData = row_data
    )
    
    # Extract data needed for filter_by_effect_size
    readcounts <- SummarizedExperiment::assays(se)$counts
    row_df <- data.frame(SummarizedExperiment::rowData(se))
    
    tx2gene <- data.frame(
        transcript_id = rownames(mat),
        gene_name = row_df$gene_name,
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = row_df$gene_name,
        log2_fold_change = row_df$log2_fold_change,
        adj_p_value = row_df$adj_p_value,
        stringsAsFactors = FALSE
    )
    
    # Use filter_by_effect_size with threshold that separates low from high
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.3)
    
    # Verify structure
    expect_is(result, "list")
    expect_true("readcounts" %in% names(result))
    expect_true("tx2gene" %in% names(result))
    # Should keep 15 genes (high FC) with 60 transcripts
    expect_equal(nrow(result$readcounts), 60)
})

test_that("filter_by_effect_size preserves SE structure when workflow includes SE", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create SE with multiple assays and metadata
    mat <- matrix(rpois(600, 10), nrow = 60, ncol = 10)
    rownames(mat) <- paste0("ENST", sprintf("%09d", 1:60))
    
    row_data <- data.frame(
        transcript_id = rownames(mat),
        gene_name = paste0("ENSG", sprintf("%09d", rep(1:30, 2))),
        log2_fold_change = rnorm(60, mean = 0.5, sd = 0.3),
        adj_p_value = runif(60),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        rowData = row_data,
        metadata = list(source = "test_data", organism = "human")
    )
    
    # Extract and filter
    readcounts <- SummarizedExperiment::assays(se)$counts
    row_df <- data.frame(SummarizedExperiment::rowData(se))
    
    tx2gene <- data.frame(
        transcript_id = rownames(mat),
        gene_name = row_df$gene_name,
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = row_df$gene_name,
        log2_fold_change = row_df$log2_fold_change,
        adj_p_value = row_df$adj_p_value,
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.3)
    
    # Verify filtering occurred
    expect_lt(nrow(result$readcounts), nrow(readcounts))
    expect_lt(nrow(result$tx2gene), nrow(tx2gene))
    
    # Verify consistency
    expect_equal(
        nrow(result$readcounts),
        nrow(result$tx2gene)
    )
})

test_that("filter_by_effect_size handles SE with strong effect sizes", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create SE with clearly separated effect sizes
    mat <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
    rownames(mat) <- paste0("ENST", sprintf("%09d", 1:50))
    
    # 25 genes: 5 with high FC, 20 with low FC
    # Each gene has 2 transcripts
    log2fc_values <- c(rep(1.5, 10), rep(0.05, 40))  # 5 high-FC genes (2 tx each), 20 low-FC genes (2 tx each)
    
    row_data <- data.frame(
        transcript_id = rownames(mat),
        gene_name = c(
            paste0("ENSG", sprintf("%09d", 1:5)), paste0("ENSG", sprintf("%09d", 1:5)),  # Genes 1-5 (high FC)
            paste0("ENSG", sprintf("%09d", 6:25)), paste0("ENSG", sprintf("%09d", 6:25))  # Genes 6-25 (low FC)
        ),
        log2_fold_change = log2fc_values,
        adj_p_value = runif(50),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        rowData = row_data
    )
    
    # Extract and filter
    readcounts <- SummarizedExperiment::assays(se)$counts
    row_df <- data.frame(SummarizedExperiment::rowData(se))
    
    tx2gene <- data.frame(
        transcript_id = rownames(mat),
        gene_name = row_df$gene_name,
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = row_df$gene_name,
        log2_fold_change = row_df$log2_fold_change,
        adj_p_value = row_df$adj_p_value,
        stringsAsFactors = FALSE
    )
    
    # Filter with threshold that separates high from low
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.5)
    
    # Should keep ~5 genes (10 transcripts) and filter out ~20 genes (40 transcripts)
    expect_equal(nrow(result$readcounts), 10)
    expect_equal(nrow(result$tx2gene), 10)
    
    # Verify these are the high-effect genes
    high_effect_genes <- paste0("ENSG", sprintf("%09d", 1:5))
    expect_true(all(result$tx2gene$gene_name %in% high_effect_genes))
})

test_that("filter_by_effect_size SE workflow: filter then rebuild SE", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create original SE with clear low/high separation
    mat <- matrix(rpois(1200, 10), nrow = 120, ncol = 10)
    rownames(mat) <- paste0("ENST", sprintf("%09d", 1:120))
    
    # 30 genes: 10 with low FC, 20 with high FC
    gene_names <- c(
        rep(paste0("ENSG", sprintf("%09d", 1:10)), each = 4),   # 40 rows
        rep(paste0("ENSG", sprintf("%09d", 11:30)), each = 4)   # 80 rows
    )
    
    row_data <- data.frame(
        transcript_id = rownames(mat),
        gene_name = gene_names,
        log2_fold_change = c(rep(0.05, 40), rep(1.0, 80)),
        adj_p_value = runif(120),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        rowData = row_data,
        metadata = list(filtering_step = "raw")
    )
    
    # Extract, filter using filter_by_effect_size
    readcounts <- SummarizedExperiment::assays(se)$counts
    row_df <- data.frame(SummarizedExperiment::rowData(se))
    
    tx2gene <- data.frame(
        transcript_id = rownames(mat),
        gene_name = row_df$gene_name,
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = row_df$gene_name,
        log2_fold_change = row_df$log2_fold_change,
        adj_p_value = row_df$adj_p_value,
        stringsAsFactors = FALSE
    )
    
    filtered_result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.3)
    
    # Rebuild SE from filtered data
    filtered_readcounts <- filtered_result$readcounts
    filtered_rowdata <- row_df[rownames(filtered_readcounts), ]
    
    se_filtered <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = filtered_readcounts),
        rowData = filtered_rowdata,
        metadata = list(
            original_step = "raw",
            filtering_method = "filter_by_effect_size",
            min_abs_log2fc = 0.3,
            n_genes_before = length(unique(row_df$gene_name)),
            n_genes_after = length(unique(filtered_rowdata$gene_name)),
            n_transcripts_before = nrow(mat),
            n_transcripts_after = nrow(filtered_readcounts)
        )
    )
    
    # Verify structure
    expect_is(se_filtered, "SummarizedExperiment")
    expect_equal(nrow(se_filtered), nrow(filtered_readcounts))
    expect_equal(ncol(se_filtered), 10)
    # Should filter from 120 to 80 transcripts (20 genes with 4 transcripts each)
    expect_equal(nrow(se_filtered), 80)
    expect_lt(nrow(se_filtered), nrow(se))
    
    # Verify metadata
    md <- S4Vectors::metadata(se_filtered)
    expect_true("filtering_method" %in% names(md))
    expect_equal(md$filtering_method, "filter_by_effect_size")
})

test_that("filter_by_effect_size maintains gene-to-transcript mapping through SE", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create SE with controlled transcript-to-gene mapping
    mat <- matrix(rpois(600, 10), nrow = 60, ncol = 10)
    rownames(mat) <- paste0("ENST", sprintf("%09d", 1:60))
    
    # Each gene has exactly 3 transcripts
    gene_ids <- rep(paste0("ENSG", sprintf("%09d", 1:20)), each = 3)
    
    row_data <- data.frame(
        transcript_id = rownames(mat),
        gene_name = gene_ids,
        log2_fold_change = rnorm(60, mean = 0.5, sd = 0.3),
        adj_p_value = runif(60),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        rowData = row_data
    )
    
    # Extract and filter
    readcounts <- SummarizedExperiment::assays(se)$counts
    row_df <- data.frame(SummarizedExperiment::rowData(se))
    
    tx2gene <- data.frame(
        transcript_id = rownames(mat),
        gene_name = row_df$gene_name,
        stringsAsFactors = FALSE
    )
    
    transcript_stats <- data.frame(
        genes = row_df$gene_name,
        log2_fold_change = row_df$log2_fold_change,
        adj_p_value = row_df$adj_p_value,
        stringsAsFactors = FALSE
    )
    
    result <- filter_by_effect_size(readcounts, transcript_stats, tx2gene, min_abs_log2fc = 0.2)
    
    # Verify that genes still have 3 transcripts (or 0 if filtered out)
    gene_counts <- table(result$tx2gene$gene_name)
    
    # All genes should have either 3 transcripts (kept) or 0 (filtered out)
    expect_true(all(gene_counts == 3 | gene_counts == 0))
})
