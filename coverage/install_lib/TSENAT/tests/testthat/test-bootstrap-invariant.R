# Bootstrap / point-estimator transformation invariant (auditx follow-up, 2026-08)
#
# The bootstrap must resample from EXACTLY the point-estimate proportions
# T(x, l, c) = (x/l + c)/sum(x/l + c) and be centered on the point estimate.
# These tests lock the invariant for pseudocount > 0 with varying effective
# lengths (the case that broke the previous depth-rescale-first
# implementation), and the q = 0 raw-support semantics of the resampling.

test_that(".prepare_bootstrap_resample reproduces point-estimate proportions exactly", {
    x <- c(120, 45, 8, 0, 300)
    leff <- c(1200, 900, 600, 400, 1500)
    pc <- 0.5

    prep <- TSENAT:::.prepare_bootstrap_resample(x, leff, pc)
    p_hat <- (prep$x + prep$pseudocount)/sum(prep$x + prep$pseudocount)
    p_point <- (x/leff + pc)/sum(x/leff + pc)

    expect_equal(p_hat, p_point, tolerance = 1e-12)
    expect_equal(sum(prep$x), sum(x), tolerance = 1e-12)  # depth preserved
    expect_gt(prep$scale_k, 1)

    # Without effective lengths: identity, pseudocount untouched
    prep0 <- TSENAT:::.prepare_bootstrap_resample(x, NULL, pc)
    expect_identical(prep0$x, x)
    expect_identical(prep0$pseudocount, pc)

    # Vector pseudocount embeds elementwise on the abundance scale
    pc_vec <- c(0.5, 0.5, 1, 2, 0.25)
    prep_v <- TSENAT:::.prepare_bootstrap_resample(x, leff, pc_vec)
    p_hat_v <- prep_v$x/sum(prep_v$x)
    p_point_v <- (x/leff + pc_vec)/sum(x/leff + pc_vec)
    expect_equal(p_hat_v, p_point_v, tolerance = 1e-12)
})

test_that("read-level bootstrap is centered on the point estimate (pc > 0, leff)", {
    set.seed(7)
    x <- c(150, 60, 12, 1, 250, 0)
    leff <- c(1300, 950, 700, 450, 1400, 800)
    pc <- 0.5
    nboot <- 3000

    for (q in c(1, 2)) {
        point <- TSENAT:::.calculate_tsallis_entropy(x, q = q, norm = FALSE,
            pseudocount = pc, effective_length = leff)
        boot <- TSENAT:::.bootstrap_resample_optimized(x = x, q = q, norm = FALSE,
            nboot = nboot, log_base = exp(1), pseudocount = pc, what = "S",
            paired = FALSE, effective_length = leff)
        expect_true(all(is.finite(boot)))
        expect_lt(abs(mean(boot) - point), 0.02)
    }
})

test_that("bootstrap CI point estimate equals the estimator on raw input (pc > 0, leff)", {
    x <- c(150, 60, 12, 1, 250, 0)
    leff <- c(1300, 950, 700, 450, 1400, 800)
    pc <- 0.5

    ci_data <- TSENAT:::.bootstrap_compute_ci(x = x, q = 1, norm = FALSE,
        nboot = 200, ci = 0.95, method = "percentile", log_base = exp(1),
        pseudocount = pc, what = "S", paired = FALSE, effective_length = leff,
        min_valid_frac = 0.75, resample_by = "read", counts_matrix = NULL)

    point <- TSENAT:::.calculate_tsallis_entropy(x, q = 1, norm = FALSE,
        pseudocount = pc, effective_length = leff)
    expect_equal(as.numeric(ci_data$point_est), as.numeric(point), tolerance = 1e-12)
})

test_that("q=0 bootstrap carries the support distribution of the draws (pc = 0)", {
    set.seed(11)
    x <- c(150, 60, 12, 1, 250, 0)
    leff <- c(1300, 950, 700, 450, 1400, 800)

    boot <- TSENAT:::.bootstrap_resample_optimized(x = x, q = 0, norm = FALSE,
        nboot = 3000, log_base = exp(1), pseudocount = 0, what = "S",
        paired = FALSE, effective_length = leff)

    # The zero isoform never enters the resampling probabilities (pc = 0), so
    # each replicate's S0 is the DRAWN support minus one. Its expectation is
    # sum_i [1 - (1 - p_i)^N] - 1 over the 5 observed isoforms (the rare
    # 1-count isoform is legitimately missed in some draws).
    p_hat <- (x/leff)/sum(x/leff)
    N <- sum(x)
    expected <- sum(1 - (1 - p_hat)^N) - 1
    expect_lt(abs(mean(boot) - expected), 0.05)
    expect_true(length(unique(boot)) > 1)  # not degenerate
})

test_that("paired block bootstrap is centered on the point estimate (pc > 0, leff)", {
    set.seed(23)
    x <- rep(c(52, 48, 45, 50, 47), each = 2) + rep(c(-1, 1, 2, -2, 0), each = 2)
    x <- as.numeric(x)
    leff <- rep(c(1100, 950, 900, 1200, 1000), each = 2)
    pc <- 0.5

    point <- TSENAT:::.calculate_tsallis_entropy(x, q = 1, norm = FALSE,
        pseudocount = pc, effective_length = leff)
    boot <- TSENAT:::.bootstrap_resample_optimized(x = x, q = 1, norm = FALSE,
        nboot = 5000, log_base = exp(1), pseudocount = pc, what = "S",
        paired = TRUE, effective_length = leff)

    expect_true(all(is.finite(boot)))
    expect_lt(abs(mean(boot) - point), 0.03)
})
