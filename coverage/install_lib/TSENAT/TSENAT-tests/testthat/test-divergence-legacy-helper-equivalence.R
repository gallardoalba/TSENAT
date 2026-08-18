# =============================================================================
# test-divergence-legacy-helper-equivalence.R
#
# The repository must have ONE Tsallis-divergence definition.
# The legacy back-compat helper .compute_tsallis_divergence() (used only by
# .calculate_divergence_bootstrap) previously implemented
#   (1 - sum p^q r^(1-q))/(q - 1) with abs()
# while the canonical path implements the signed Furuichi form
#   (sum p^q r^(1-q) - 1)/(q - 1)
# without abs(). This test proves the legacy helper now matches the canonical
# scalar implementation across the q grid.
# =============================================================================

test_that(".compute_tsallis_divergence matches the canonical scalar (all-positive support)", {
    set.seed(123)
    x <- c(70, 20, 10, 5, 3)
    y <- c(40, 40, 20, 5, 2)
    p <- x/sum(x)
    r <- y/sum(y)

    for (q in c(0.25, 0.5, 0.9, 1, 1.1, 1.5, 2, 3)) {
        legacy <- TSENAT:::.compute_tsallis_divergence(p, r, q)
        canonical <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        expect_equal(legacy, canonical, tolerance = 1e-8, info = paste("q =", q))
    }
})

test_that(".compute_tsallis_divergence identity is zero and q=0 is rejected", {
    p <- c(0.4, 0.3, 0.2, 0.1)
    expect_equal(
        TSENAT:::.compute_tsallis_divergence(p, p, q = 2),
        0, tolerance = 1e-12
    )
    expect_true(is.na(TSENAT:::.compute_tsallis_divergence(p, p, q = 0)))
})

test_that(".compute_tsallis_divergence handles log_base consistently", {
    p <- c(0.6, 0.4)
    r <- c(0.3, 0.7)

    # q -> 1 limit (KL): single log_base division
    legacy_e <- TSENAT:::.compute_tsallis_divergence(p, r, q = 1, log_base = exp(1))
    manual_e <- sum(p * log(p/r))
    expect_equal(legacy_e, manual_e, tolerance = 1e-12)

    legacy_10 <- TSENAT:::.compute_tsallis_divergence(p, r, q = 1, log_base = 10)
    expect_equal(legacy_10, manual_e/log(10), tolerance = 1e-12)
})
