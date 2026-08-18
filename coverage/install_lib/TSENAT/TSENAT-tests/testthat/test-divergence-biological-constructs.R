# ============================================================================
# test-divergence-biological-constructs.R
#
# Additional property tests requested in review_divergence.md (sections 10-12):
#   A. Coarse-graining behavior
#   B. (q=2 analytic formula lives in test-divergence-kernel-validation.R)
#   C. Biological-replicate replication invariance
#   D. Bootstrap consistency (mean of bootstrap ~ point estimate)
#   E. Support behavior in BOTH directions (divergence is directed)
#   F. Biological construct validity: abundance-only change, pure isoform
#      switch, rare-isoform remodeling, dominant-isoform remodeling
# ============================================================================

test_that("Coarse-graining: merging categories changes the divergence", {
    P <- c(0.4, 0.3, 0.2, 0.1)
    Q <- c(0.2, 0.3, 0.4, 0.1)
    # Merge the first two categories
    Pc <- c(0.7, 0.2, 0.1)
    Qc <- c(0.5, 0.4, 0.1)

    for (q in c(0.5, 1, 2)) {
        d_fine <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = q,
            pseudocount = 0)
        d_coarse <- TSENAT:::.tsallis_divergence_scalar(Pc * 100, Qc * 100,
            q_val = q, pseudocount = 0)
        expect_true(is.finite(d_fine) && d_fine >= 0)
        expect_true(is.finite(d_coarse) && d_coarse >= 0)
        # Coarse-graining genuinely changes the divergence for q > 0
        expect_true(abs(d_fine - d_coarse) > 1e-10, info = paste("q =", q))
    }
})

test_that("Replication invariance: duplicated replicates preserve pooled divergence", {
    # Pooling the same composition across duplicated biological replicates
    # must not change the divergence (pseudocount = 0: no per-bin
    # regularization that would depend on the number of bins).
    x <- c(40, 10, 5)
    y <- c(15, 25, 8)
    x2 <- rep(x, times = 3)
    y2 <- rep(y, times = 3)

    for (q in c(0, 0.5, 1, 2)) {
        d1 <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        d2 <- TSENAT:::.tsallis_divergence_scalar(x2, y2, q_val = q, pseudocount = 0)
        expect_equal(d1, d2, tolerance = 1e-10, info = paste("q =", q))
    }
})

test_that("Bootstrap consistency: bootstrap mean centers on the point estimate", {
    set.seed(42)
    x_mat <- matrix(rpois(6 * 4, lambda = 60), nrow = 6)  # 6 isoforms x 4 ctrl
    y_mat <- matrix(rpois(6 * 4, lambda = 55), nrow = 6)  # 6 isoforms x 4 trt
    colnames(x_mat) <- paste0("C", 1:4)
    colnames(y_mat) <- paste0("T", 1:4)

    x <- rowSums(x_mat)
    y <- rowSums(y_mat)

    point <- TSENAT:::.tsallis_divergence_vector(x, y, 1, pseudocount = 0.5)

    plan <- TSENAT:::.bootstrap_sample_indices(x_mat, y_mat, 500, NULL)
    dist <- vapply(seq_len(500), function(b) {
        xb <- rowSums(x_mat[, plan$ctrl[[b]], drop = FALSE])
        yb <- rowSums(y_mat[, plan$trt[[b]], drop = FALSE])
        TSENAT:::.tsallis_divergence_vector(xb, yb, 1, pseudocount = 0.5)
    }, numeric(1))

    expect_true(is.finite(point))
    expect_gt(mean(dist, na.rm = TRUE), point - 0.05)
    expect_lt(mean(dist, na.rm = TRUE), point + 0.05)
})

test_that("Support is directional: P>0/Q=0 differs from P=0/Q>0", {
    # Directedness: support violations depend on the direction.
    P <- c(0.7, 0.3, 0.0)
    Q <- c(0.4, 0.0, 0.6)

    d_pq <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 2,
        pseudocount = 0)
    d_qp <- TSENAT:::.tsallis_divergence_scalar(Q * 100, P * 100, q_val = 2,
        pseudocount = 0)

    expect_true(is.finite(d_pq) && is.finite(d_qp))
    expect_true(abs(d_pq - d_qp) > 1e-8)  # directed divergence
    # q=0 raw support: D_0(P||Q) = R-mass on P's zero support = Q_3 = 0.6;
    # D_0(Q||P) = P-mass on Q's zero support = P_2 = 0.3
    expect_equal(TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 0,
        pseudocount = 0.5), 0.6)
    expect_equal(TSENAT:::.tsallis_divergence_scalar(Q * 100, P * 100, q_val = 0,
        pseudocount = 0.5), 0.3)
})

test_that("Scenario 1: abundance change with identical isoform usage gives D=0", {
    x <- c(80, 20)
    y <- c(800, 200)  # same composition, 10x expression
    for (q in c(0, 0.5, 1, 2)) {
        # Exact equality requires pseudocount = 0 (with a pseudocount the
        # regularized compositions differ by an O(pc/N) amount).
        d <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        expect_equal(d, 0, tolerance = 1e-10, info = paste("q =", q))
    }
    # With the default pseudocount the divergence remains negligible
    d_pc <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 1, pseudocount = 0.5)
    expect_lt(d_pc, 0.01)
})

test_that("Scenario 2: pure isoform switch gives D>0 at all q", {
    x <- c(80, 20)
    y <- c(20, 80)
    for (q in c(0.5, 1, 2)) {
        d <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        expect_true(d > 0, info = paste("q =", q))
    }
})

test_that("Scenarios 3-4: rare vs dominant remodeling show the expected q-sensitivity", {
    # Scenario 3: rare-isoform remodeling (low-abundance tails change)
    S3_P <- c(0.90, 0.09, 0.01)
    S3_Q <- c(0.90, 0.05, 0.05)
    # Scenario 4: dominant-isoform remodeling (abundant components change)
    S4_P <- c(0.8, 0.1, 0.1)
    S4_Q <- c(0.5, 0.25, 0.25)

    d_low <- function(P, Q) TSENAT:::.tsallis_divergence_scalar(P * 1000, Q * 1000,
        q_val = 0.5, pseudocount = 0)
    d_high <- function(P, Q) TSENAT:::.tsallis_divergence_scalar(P * 1000, Q * 1000,
        q_val = 2, pseudocount = 0)

    # Rare remodeling is RELATIVELY more visible at low q than dominant
    # remodeling: ratio D_0.5 / D_2 must be larger for scenario 3.
    ratio_S3 <- d_low(S3_P, S3_Q)/d_high(S3_P, S3_Q)
    ratio_S4 <- d_low(S4_P, S4_Q)/d_high(S4_P, S4_Q)

    expect_true(d_low(S3_P, S3_Q) > 0)
    expect_true(d_high(S3_P, S3_Q) > 0)
    expect_true(d_low(S4_P, S4_Q) > 0)
    expect_true(d_high(S4_P, S4_Q) > 0)
    expect_true(ratio_S3 > ratio_S4)
})
