context("Tsallis entropy calculations")

test_that("Tsallis entropy calculation is mathematically correct", {
    # Mathematical reference:
    # H_q(p) = (1 - sum(p_i^q)) / (q-1)
    # For q=1, H_1(p) = -sum(p_i * log2(p_i))
    read_counts <- c(0, 0, 5, 4, 1)
    p <- read_counts / sum(read_counts)

    # q = 2 (unnormalized)
    q2 <- 2
    manual_q2 <- (1 - sum(p^q2)) / (q2 - 1)
    tsallis_q2 <- calculate_tsallis_entropy(read_counts, q = q2, norm = FALSE)
    expect_equal(tsallis_q2, manual_q2, tolerance = 1e-8)

    # q = 2 (normalized)
    max_tsallis_q2 <- (1 - length(read_counts)^(1 - q2)) / (q2 - 1)
    manual_q2_norm <- manual_q2 / max_tsallis_q2
    tsallis_q2_norm <- calculate_tsallis_entropy(read_counts, q = q2, norm = TRUE)
    expect_equal(tsallis_q2_norm, manual_q2_norm, tolerance = 1e-8)
    expect_true(tsallis_q2_norm <= 1 && tsallis_q2_norm >= 0)

    # q = 1 (Shannon, unnormalized) -- use natural log by default
    manual_shannon <- -sum(ifelse(p > 0, p * log(p), 0))
    tsallis_q1 <- calculate_tsallis_entropy(read_counts, q = 1, norm = FALSE)
    expect_equal(tsallis_q1, manual_shannon, tolerance = 1e-8)

    # q = 1 (Shannon, normalized)
    manual_shannon_norm <- manual_shannon / log(length(read_counts))
    tsallis_q1_norm <- calculate_tsallis_entropy(read_counts, q = 1, norm = TRUE)
    expect_equal(tsallis_q1_norm, manual_shannon_norm, tolerance = 1e-8)
    expect_true(tsallis_q1_norm <= 1 && tsallis_q1_norm >= 0)

    # q = 1.5 (unnormalized)
    q15 <- 1.5
    manual_q15 <- (1 - sum(p^q15)) / (q15 - 1)
    tsallis_q15 <- calculate_tsallis_entropy(read_counts, q = q15, norm = FALSE)
    expect_equal(tsallis_q15, manual_q15, tolerance = 1e-8)

    # Vector q
    qvec <- c(1, 1.5, 2)
    tsallis_vec <- calculate_tsallis_entropy(read_counts, q = qvec, norm = FALSE)
    manual_vec <- vapply(qvec, function(qi) {
        if (abs(qi - 1) < .Machine$double.eps^0.5) {
            -sum(ifelse(p > 0, p * log(p), 0))
        } else {
            (1 - sum(p^qi)) / (qi - 1)
        }
    }, numeric(1))
    expect_equal(as.numeric(tsallis_vec),
        as.numeric(manual_vec),
        tolerance = 1e-8
    )
    expect_named(tsallis_vec, paste0("q=", qvec))

    # Edge cases
    expect_true(is.nan(calculate_tsallis_entropy(c(1), q = 2)))
    expect_true(is.na(calculate_tsallis_entropy(c(0, 0), q = 2)))
    expect_error(calculate_tsallis_entropy(read_counts, q = 0))
    expect_error(calculate_tsallis_entropy(read_counts, q = -1))
})
