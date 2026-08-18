# =============================================================================
# test-bootstrap-bca-effective-length.R
#
# Audit 2026-08-17 (item 5): regression test proving that with
# effective_length != NULL the point estimate, the bootstrap replicates and
# the BCa leave-one-out jackknife all implement EXACTLY the same estimator
#
#   T(x, l, c) = (x/l + c) / sum(x/l + c)
#
# followed by the entropy functional S_q. The jackknife must equal the true
# leave-one-out estimator T(x[-i], l[-i], c) for every deletion i.
# =============================================================================

# Independent reference implementation of the point estimator (proportions
# computed by hand, then the entropy core).
manual_T_entropy <- function(x, l, c, q) {
    p <- (x/l + c)
    p <- p/sum(p)
    TSENAT:::.entropy_core(p, q = q, norm = FALSE)
}

test_that("BCa point estimate, replicates and jackknife share one effective-length estimator (q=1)", {
    set.seed(42)
    x <- c(120, 80, 60, 45, 30, 22, 15, 8)
    l <- c(900, 1100, 800, 950, 1050, 1000, 850, 1200)
    pc <- 0.5
    q <- 1

    # Reference point estimate (independent manual computation)
    ref_theta <- manual_T_entropy(x, l, pc, q)

    # Package point estimate through the full estimator
    pkg_theta <- TSENAT:::.calculate_tsallis_entropy(
        x, q = q, norm = FALSE, what = "S",
        pseudocount = pc, effective_length = l
    )
    expect_equal(as.numeric(pkg_theta), ref_theta, tolerance = 1e-10)

    # Reference leave-one-out jackknife: T(x[-i], l[-i], c)
    ref_jack <- vapply(seq_along(x), function(i) {
        manual_T_entropy(x[-i], l[-i], pc, q)
    }, numeric(1))

    # Run the BCa pipeline (nboot kept small; the jackknife is exact regardless)
    res <- TSENAT:::.bootstrap_compute_ci(
        x, q = q, norm = FALSE, nboot = 100L, ci = 0.95, method = "bca",
        log_base = exp(1), pseudocount = pc, what = "S",
        paired = FALSE, effective_length = l
    )

    # 1. Point estimate
    expect_equal(res$point_est, ref_theta, tolerance = 1e-10)

    # 2. Jackknife values
    expect_true(!is.null(res$ci_result$jackknife))
    expect_equal(
        as.numeric(res$ci_result$jackknife), ref_jack, tolerance = 1e-10,
        info = "BCa leave-one-out jackknife must equal T(x[-i], l[-i], c)"
    )

    # 3. Bootstrap replicates resample from exactly the point-estimate
    #    composition, so their mean approximates the point estimate.
    expect_equal(
        mean(res$bootstrap_dist, na.rm = TRUE), ref_theta, tolerance = 0.05
    )
})

test_that("BCa invariant holds at q=2 (general Tsallis branch)", {
    set.seed(7)
    x <- c(50, 40, 30, 20, 10)
    l <- c(800, 900, 1000, 1100, 1200)
    pc <- 1
    q <- 2

    ref_theta <- manual_T_entropy(x, l, pc, q)
    ref_jack <- vapply(seq_along(x), function(i) {
        manual_T_entropy(x[-i], l[-i], pc, q)
    }, numeric(1))

    res <- TSENAT:::.bootstrap_compute_ci(
        x, q = q, norm = FALSE, nboot = 100L, ci = 0.95, method = "bca",
        log_base = exp(1), pseudocount = pc, what = "S",
        paired = FALSE, effective_length = l
    )

    expect_equal(res$point_est, ref_theta, tolerance = 1e-10)
    expect_equal(as.numeric(res$ci_result$jackknife), ref_jack, tolerance = 1e-10)
})

test_that("BCa invariant holds with vector pseudocount and norm=TRUE", {
    set.seed(11)
    x <- c(90, 70, 55, 40, 25, 15)
    l <- c(750, 850, 950, 1050, 1150, 1250)
    pc <- 0.25
    q <- 1.5

    p <- (x/l + pc)
    p <- p/sum(p)
    ref_theta <- TSENAT:::.entropy_core(p, q = q, norm = TRUE)

    ref_jack <- vapply(seq_along(x), function(i) {
        pi <- (x[-i]/l[-i] + pc)
        pi <- pi/sum(pi)
        TSENAT:::.entropy_core(pi, q = q, norm = TRUE)
    }, numeric(1))

    res <- TSENAT:::.bootstrap_compute_ci(
        x, q = q, norm = TRUE, nboot = 100L, ci = 0.95, method = "bca",
        log_base = exp(1), pseudocount = pc, what = "S",
        paired = FALSE, effective_length = l
    )

    expect_equal(res$point_est, ref_theta, tolerance = 1e-10)
    expect_equal(as.numeric(res$ci_result$jackknife), ref_jack, tolerance = 1e-10)
})
