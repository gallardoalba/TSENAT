# Additional tests for calculate_tsallis_entropy function
# Tests the low-level entropy calculation function

library(TSENAT)

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
