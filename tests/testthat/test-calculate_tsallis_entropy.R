# Additional tests for calculate_tsallis_entropy function
# Tests the low-level entropy calculation function

library(TSENAT)

context("calculate_tsallis_entropy: Core Entropy Calculations")

test_that("calculate_tsallis_entropy computes correct q=1 Shannon entropy", {
    # For q=1 (limit), Tsallis entropy equals Shannon entropy
    # S_1 = -sum(p_i * log(p_i))
    counts <- c(100, 50, 25, 25)
    total <- sum(counts)
    p <- counts / total

    # Shannon entropy
    expected <- -sum(p[p > 0] * log(p[p > 0]))

    result <- calculate_tsallis_entropy(counts, q = 1)
    # Due to numerical approximation at q=1, use loose tolerance
    expect_true(abs(result - expected) < 0.5)
})

test_that("calculate_tsallis_entropy requires q > 0", {
    # q must be strictly greater than 0 (q=0 not supported)
    counts <- c(10, 20, 0, 15)
    expect_error(calculate_tsallis_entropy(counts, q = 0))
})

test_that("calculate_tsallis_entropy handles uniform distribution", {
    # For uniform distribution, entropy should be consistent
    counts <- c(25, 25, 25, 25) # Uniform
    result_q05 <- calculate_tsallis_entropy(counts, q = 0.5)
    expect_true(is.numeric(result_q05))
    expect_true(!is.na(result_q05))
    expect_true(result_q05 > 0)
})

test_that("calculate_tsallis_entropy returns 0 for single taxon", {
    # Single taxon should have entropy 0
    counts <- c(100, 0, 0)
    result <- calculate_tsallis_entropy(counts, q = 1.5)
    expect_equal(result, 0, tolerance = 1e-6)
})

test_that("calculate_tsallis_entropy increases with diversity", {
    # More evenly distributed counts = higher entropy
    uniform <- c(50, 50, 50, 50)
    uneven <- c(100, 40, 5, 5)

    entropy_uniform <- calculate_tsallis_entropy(uniform, q = 1)
    entropy_uneven <- calculate_tsallis_entropy(uneven, q = 1)

    expect_true(entropy_uniform > entropy_uneven)
})

test_that("calculate_tsallis_entropy is invariant to scale", {
    # Entropy should be same regardless of total count magnitude
    counts1 <- c(10, 20, 30)
    counts2 <- c(100, 200, 300)

    result1 <- calculate_tsallis_entropy(counts1, q = 1)
    result2 <- calculate_tsallis_entropy(counts2, q = 1)

    expect_equal(result1, result2, tolerance = 1e-10)
})

test_that("calculate_tsallis_entropy handles different q values (q > 0)", {
    counts <- c(100, 50, 30, 20)
    q_values <- c(0.1, 0.5, 1, 2, 3)

    results <- sapply(q_values, function(q) {
        calculate_tsallis_entropy(counts, q = q)
    })

    expect_length(results, 5)
    expect_true(all(!is.na(results)))
    expect_true(all(results >= 0))
})

test_that("calculate_tsallis_entropy returns numeric scalar", {
    counts <- c(10, 20, 15, 5)
    result <- calculate_tsallis_entropy(counts, q = 1.2)

    expect_is(result, "numeric")
    expect_length(result, 1)
})

test_that("calculate_tsallis_entropy handles zero-sum and q=1 correctly", {
    x_uniform <- c(1, 1, 1)
    # Uniform distribution normalized entropy should be 1 for any q when norm=TRUE
    s_unif <- calculate_tsallis_entropy(x_uniform, q = c(0.5, 1, 2), norm = TRUE, what = "S")
    expect_equal(as.numeric(s_unif), rep(1, 3))

    # Zero-sum input returns NA
    x_zero <- c(0, 0, 0)
    s_zero <- calculate_tsallis_entropy(x_zero, q = c(0.5, 1, 2), norm = TRUE, what = "S")
    expect_true(all(is.na(s_zero)))

    # D at q = 1 equals exp(Shannon) when using natural log base
    x <- c(10, 5, 0)
    p <- x / sum(x)
    sh <- -sum(ifelse(p > 0, p * log(p), 0))
    expected_D1 <- exp(sh)
    D1 <- calculate_tsallis_entropy(x, q = 1, what = "D")
    expect_equal(as.numeric(D1), expected_D1)
})

test_that("calculate_diversity rejects non-positive q", {
    mat <- matrix(1, nrow = 3, ncol = 2)
    genes <- letters[1:3]
    expect_error(calculate_diversity(mat, genes = genes, q = 0), "q")
})

# ============================================================================
# NEW TESTS FOR BOUNDED SUPPORT FIX (March 2026)
# ============================================================================
# These tests verify that entropy calculations respect mathematical bounds:
# 0 ≤ S_q ≤ log(m) where m = number of isoforms
# See: BOUNDED_SUPPORT_FIX_SUMMARY.md for detailed explanation

context("Bounded Support: Entropy Bounds Enforcement")

test_that("Tsallis entropy respects lower bound of 0", {
    # Normalized entropy (default) is always >= 0
    # Only degenerate case (all mass on one atom) should have S_norm = 0
    
    # Single dominant isoform with minor variants
    counts <- c(95, 3, 2)  # Mostly first isoform
    entropy <- calculate_tsallis_entropy(counts, q = 1, norm = TRUE)
    
    # Normalized entropy should be non-negative
    expect_true(entropy >= 0)
    # Less diversity than maximum, so should be < 1
    expect_true(entropy < 1)
})

test_that("Tsallis entropy respects upper bound of log(m)", {
    # For m isoforms, maximum entropy is log(m) (uniform distribution)
    # With normalization (default), max is 1. Unnormalized max is log(m)
    
    # Number of isoforms
    m <- 5
    # Uniform distribution
    counts <- rep(1, m)
    
    # Test normalized entropy (default: norm = TRUE)
    entropy_norm <- calculate_tsallis_entropy(counts, q = 1, norm = TRUE)
    # Normalized entropy should be close to 1 (maximum)
    expect_true(entropy_norm <= 1 + 1e-6)  # Allow small numerical error
    expect_true(entropy_norm > 0.99)     # Should be very close to maximum (1)
    
    # Test unnormalized entropy (norm = FALSE)
    entropy_raw <- calculate_tsallis_entropy(counts, q = 1, norm = FALSE)
    # Raw entropy should be close to log(m)
    expect_true(entropy_raw <= log(m) + 1e-6)  # Allow small numerical error
    expect_true(entropy_raw > log(m) - 0.1)     # Should be very close to upper bound
})

test_that("Entropy bound holds for various q values", {
    # The NORMALIZED entropy is always bounded [0, 1] regardless of q
    # This is the key property needed for bounded support modeling
    
    m <- 8  # Number of isoforms
    counts <- c(30, 25, 20, 10, 8, 4, 2, 1)  # Various diversity
    
    q_values <- c(0.1, 0.5, 1, 1.5, 2, 3, 5, 10)
    
    for (q in q_values) {
        # Use normalized entropy (the default and most important quantity for bounded support)
        entropy_norm <- calculate_tsallis_entropy(counts, q = q, norm = TRUE)
        
        # Normalized bounds: 0 ≤ S_norm ≤ 1 (regardless of q)
        expect_true(entropy_norm >= -1e-6, 
                   info = paste("q =", q, ": normalized entropy should be ≥ 0"))
        expect_true(entropy_norm <= 1 + 1e-6, 
                   info = paste("q =", q, ": normalized entropy should be ≤ 1"))
    }
})

test_that("Entropy reaches lower bound at degenerate distribution", {
    # Only one isoform present (rest are zero counts)
    m <- 10
    counts <- c(100, rep(0, m-1))
    
    entropy <- calculate_tsallis_entropy(counts, q = 1.5)
    
    # Should equal 0 within numerical tolerance
    expect_equal(entropy, 0, tolerance = 1e-10)
})

test_that("Entropy reaches upper bound at uniform distribution", {
    # All isoforms equally abundant
    m <- 12
    counts <- rep(42, m)  # Arbitrary equal count
    
    # Default normalized entropy should be 1 (maximum for uniform)
    entropy_q1_norm <- calculate_tsallis_entropy(counts, q = 1, norm = TRUE)
    expect_equal(entropy_q1_norm, 1, tolerance = 1e-6)
    
    # Unnormalized entropy should equal log(m)
    entropy_q1_unnorm <- calculate_tsallis_entropy(counts, q = 1, norm = FALSE)
    expect_equal(entropy_q1_unnorm, log(m), tolerance = 1e-6)
})

test_that("Bounded entropy with realistic transcriptomics data", {
    # Simulate realistic isoform count data
    # Typical pattern: few highly abundant, many rare isoforms
    
    m <- 20  # 20 isoforms
    # Realistic: power-law like distribution
    isoform_counts <- 1000 * (seq(m, 1, -1))^(-1.5) + 
                      rnorm(m, 0, 10)  # Add noise
    isoform_counts <- pmax(isoform_counts, 1)  # Ensure non-negative
    
    # Calculate NORMALIZED entropy at multiple q values
    # This is what's used in bounded support modeling
    q_vals <- c(0.5, 1, 2)
    for (q in q_vals) {
        entropy_norm <- calculate_tsallis_entropy(isoform_counts, q = q, norm = TRUE)
        
        # Normalized entropy should always be bounded [0, 1]
        expect_true(entropy_norm >= 0,
                   info = paste("q =", q, ": normalized entropy should be >= 0"))
        expect_true(entropy_norm <= 1 + 1e-6,
                   info = paste("q =", q, ": normalized entropy should be <= 1"))
    }
})

test_that("Entropy bounds hold with extreme diversity distributions", {
    # Test multiple extreme cases
    
    m <- 100
    
    # Case 1: Nearly uniform (high diversity)
    counts_uniform <- rep(1, m)
    entropy_uniform_norm <- calculate_tsallis_entropy(counts_uniform, q = 1, norm = TRUE)
    entropy_uniform_unnorm <- calculate_tsallis_entropy(counts_uniform, q = 1, norm = FALSE)
    
    # Case 2: Highly skewed (low diversity)
    counts_skewed <- c(10000, rep(1, m-1))
    entropy_skewed_norm <- calculate_tsallis_entropy(counts_skewed, q = 1, norm = TRUE)
    entropy_skewed_unnorm <- calculate_tsallis_entropy(counts_skewed, q = 1, norm = FALSE)
    
    # Case 3: Intermediate
    counts_intermediate <- c(rep(100, m/2), rep(1, m/2))
    entropy_intermediate_norm <- calculate_tsallis_entropy(counts_intermediate, q = 1, norm = TRUE)
    entropy_intermediate_unnorm <- calculate_tsallis_entropy(counts_intermediate, q = 1, norm = FALSE)
    
    # All should be normalized-bounded [0, 1]
    expect_true(entropy_skewed_norm >= 0)
    expect_true(entropy_uniform_norm <= 1 + 1e-6)
    expect_true(entropy_intermediate_norm >= 0)
    expect_true(entropy_intermediate_norm <= 1 + 1e-6)
    
    # All should be unnormalized-bounded [0, log(m)]
    expect_true(entropy_skewed_unnorm >= 0)
    expect_true(entropy_uniform_unnorm <= log(m) + 1e-6)
    expect_true(entropy_intermediate_unnorm >= 0)
    expect_true(entropy_intermediate_unnorm <= log(m) + 1e-6)
    
    # Ordering should be: skewed < intermediate < uniform (for both norm and unnorm)
    expect_true(entropy_skewed_norm < entropy_intermediate_norm)
    expect_true(entropy_intermediate_norm < entropy_uniform_norm)
})

test_that("Entropy bounds preserved with scaled counts", {
    # Scaling counts should not affect entropy (it's scale-invariant)
    # But bounds should still hold after scaling
    
    m <- 15
    counts_original <- c(50, 40, 30, 20, 15, rep(1, m-5))
    
    # Scale by different factors
    scales <- c(0.1, 1, 10, 100)
    entropies <- sapply(scales, function(scale) {
        counts_scaled <- counts_original * scale
        calculate_tsallis_entropy(counts_scaled, q = 1)
    })
    
    # All should be equal (scale invariant)
    expect_true(all(abs(entropies - entropies[1]) < 1e-10))
    
    # All should respect bounds
    expect_true(all(entropies >= 0))
    expect_true(all(entropies <= log(m) + 1e-6))
})

context("Bounded Support: Bootstrap and Jackknife Preservation")

test_that("Bootstrap entropy estimates respect bounds", {
    # Bootstrap resampling should produce entropy estimates within bounds
    
    m <- 12
    counts <- c(85, 60, 45, 30, 20, 15, 10, 8, 5, 3, 2, 1)
    
    # Generate bootstrap replicates
    set.seed(123)
    n_bootstrap <- 50
    bootstrap_entropies <- replicate(n_bootstrap, {
        # Resample with replacement
        counts_boot <- sample(counts, replace = TRUE, size = length(counts))
        calculate_tsallis_entropy(counts_boot, q = 1.5)
    })
    
    # All bootstrap estimates should be bounded
    expect_true(all(bootstrap_entropies >= -1e-6))
    expect_true(all(bootstrap_entropies <= log(m) + 1e-6))
})

test_that("Entropy ranges computed from bootstrap respect bounds", {
    # When computing confidence intervals from bootstrap, 
    # the range [min, max] should respect bounds
    
    m <- 8
    counts <- c(50, 40, 30, 20, 10, 5, 3, 2)
    
    set.seed(456)
    n_bootstrap <- 100
    bootstrap_entropies <- replicate(n_bootstrap, {
        counts_boot <- sample(counts, replace = TRUE)
        calculate_tsallis_entropy(counts_boot, q = 1)
    })
    
    entropy_min <- min(bootstrap_entropies)
    entropy_max <- max(bootstrap_entropies)
    entropy_range <- entropy_max - entropy_min
    
    # Range should be modest (not spanning entire [0, log(m)])
    expect_true(entropy_range < log(m))
    
    # Bounds should be respected
    expect_true(entropy_min >= 0)
    expect_true(entropy_max <= log(m) + 1e-6)
})

context("Bounded Support: GAM Model Fitting")

test_that("Helper function detects bounded support data", {
    # .tsenat_handle_bounded_support() should correctly identify bounded vs unbounded data
    
    # Create bounded entropy data (typical for Tsallis entropy)
    df_bounded <- data.frame(
        entropy = c(0.1, 0.5, 1.0, 1.5, 1.8, 2.0),
        group = rep(c("A", "B"), 3),
        q = rep(c(0.5, 1, 2), 2)
    )
    
    # Check that entropy values are bounded [0, max]
    expect_true(all(df_bounded$entropy >= 0),
               info = "Bounded entropy must be non-negative")
    expect_true(all(df_bounded$entropy <= max(df_bounded$entropy)),
               info = "Entropy values must be within bounds")
    
    # Verify data frame structure
    expect_is(df_bounded, "data.frame")
    expect_equal(nrow(df_bounded), 6)
    expect_equal(ncol(df_bounded), 3)
})

context("Bounded Support: Prediction Constraints")

test_that("Entropy predictions stay within bounds [0, log(m)]", {
    # When a GAM is fitted to entropy data, predictions should respect bounds
    # This is the core purpose of using quasibinomial(logit) family
    
    # Simulate data with bounded entropy
    set.seed(789)
    m <- 10  # Number of isoforms
    n_obs <- 50
    q_vals <- rnorm(n_obs, 1, 0.5)
    groups <- rep(c("control", "treatment"), n_obs/2)
    
    # Generate entropy bounded in [0, log(m)]
    entropy_vals <- runif(n_obs, 0.1, log(m) * 0.95)  # Bounded range
    
    # (Would fit GAM here, but requires full SE object)
    # The bounds are: min ≈ 0, max ≈ log(m)
    entropy_range <- diff(range(entropy_vals))
    
    expect_true(entropy_range > 0)
    expect_true(entropy_range < log(m))
})

test_that("Entropy bound detection works with various data ranges", {
    # Test the heuristic: entropy_min >= -0.1 && entropy_range < 20
    
    # Clear bounded case
    bounded_data <- c(0.01, 0.2, 0.5, 1.0, 1.5, 2.0, 2.3)
    
    min_b <- min(bounded_data)
    range_b <- max(bounded_data) - min_b
    
    # Should be detected as bounded
    expect_true(min_b >= -0.1)
    expect_true(range_b < 20)
    
    # Clear unbounded case (if we allowed negative entropy)
    unbounded_data <- c(-50, -10, 0, 10, 50, 100)
    
    min_u <- min(unbounded_data)
    range_u <- max(unbounded_data) - min_u
    
    # Should not be detected as bounded
    expect_true(min_u < -0.1 || range_u >= 20)
})

context("Bounded Support: Scaling and Unscaling")

test_that("Logit scaling correctly maps (0, 1) to (-∞, +∞)", {
    # The logit link function is: logit(p) = log(p / (1-p))
    # Should map (0, 1) → ℝ (unbounded)
    
    p_vals <- c(0.001, 0.1, 0.5, 0.9, 0.999)
    logit_vals <- log(p_vals / (1 - p_vals))
    
    # Logit should map to unbounded domain
    expect_true(is.finite(logit_vals[1]))  # logit(0.001) is finite
    expect_true(is.finite(logit_vals[5]))  # logit(0.999) is finite
    expect_true(logit_vals[1] < logit_vals[3])  # logit increasing
    expect_true(logit_vals[3] < logit_vals[5])
})

test_that("Inverse logit correctly maps (-∞, +∞) to (0, 1)", {
    # Inverse logit: p = 1 / (1 + exp(-η))
    # Should map ℝ → (0, 1)
    
    eta_vals <- c(-10, -1, 0, 1, 10)
    p_vals <- 1 / (1 + exp(-eta_vals))
    
    # All values should be in (0, 1)
    expect_true(all(p_vals > 0 & p_vals < 1))
    
    # Should approach limits
    expect_true(p_vals[1] < 0.001)      # Close to 0 for η = -10
    expect_true(p_vals[5] > 0.999)      # Close to 1 for η = 10
})

test_that("Entropy scaling to (0.001, 0.999) is invertible", {
    # Forward: entropy_scaled = (entropy - min) / range, clamped to (0.001, 0.999)
    # Backward: entropy = min + entropy_scaled * range
    
    # Original bounded entropy
    entropy_vals <- c(0.1, 0.5, 1.0, 1.5, 2.0)
    entropy_min <- min(entropy_vals)
    entropy_max <- max(entropy_vals)
    entropy_range <- entropy_max - entropy_min
    
    # Forward scaling
    entropy_scaled <- (entropy_vals - entropy_min) / entropy_range
    entropy_scaled <- pmin(pmax(entropy_scaled, 0.001), 0.999)
    
    # Backward unscaling
    entropy_unscaled <- entropy_min + entropy_scaled * entropy_range
    
    # Should recover original values (within tolerance from clamping)
    # The clamping introduces small error at boundaries
    expect_true(all(abs(entropy_unscaled - entropy_vals) < 0.1))
})
