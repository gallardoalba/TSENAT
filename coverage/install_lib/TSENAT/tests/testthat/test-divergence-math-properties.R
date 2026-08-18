# =============================================================================
# test-divergence-math-properties.R
#
# Mathematical property tests for the Tsallis divergence, independent of the
# implementation history. Reference distributions:
#   P = (0.7, 0.2, 0.1), Q = (0.4, 0.4, 0.2)
#
# Properties checked (following the divergence review):
#   1. Identity:            D_q(P||P) = 0 for all q
#   2. Non-negativity:      D_q(P||Q) >= 0
#   3. KL limit:            lim_{q->1} D_q(P||Q) = D_KL(P||Q)
#   4. Asymmetry:           D_q(P||Q) != D_q(Q||P) in general
#   5. Permutation invariance: joint permutation of categories
#   6. Count-scale invariance: D(x, y) = D(c*x, c*y) for c > 0
#   7. q=0 convention:      D_0 = 1 - sum_{i: P_i>0} Q_i
#   8. Support handling:    P_i>0, Q_i=0 finite under regularization
#   9. Construct validity:  identical isoform usage with different total
#                           expression -> D_q = 0
#  10. Unequal group sizes: valid as long as the isoform state space matches
# =============================================================================

test_that("Identity: D_q(P||P) = 0 for all q", {
    x <- c(70, 20, 10)
    for (q in c(0, 0.5, 1, 1.5, 2, 3)) {
        expect_equal(
            TSENAT:::.tsallis_divergence_scalar(x, x, q_val = q, pseudocount = 0),
            0, tolerance = 1e-10,
            info = paste("q =", q)
        )
    }
})

test_that("Non-negativity: D_q(P||Q) >= 0", {
    P <- c(0.7, 0.2, 0.1) * 100
    Q <- c(0.4, 0.4, 0.2) * 100
    for (q in c(0.5, 1, 1.5, 2, 3)) {
        div <- TSENAT:::.tsallis_divergence_scalar(P, Q, q_val = q, pseudocount = 0)
        expect_true(is.finite(div) && div >= 0, info = paste("q =", q))
    }
})

test_that("KL limit: D_q -> D_KL as q -> 1", {
    P <- c(0.7, 0.2, 0.1)
    Q <- c(0.4, 0.4, 0.2)
    kl_ref <- sum(P * log(P/Q))
    d_kl <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 1,
        pseudocount = 0)
    expect_equal(d_kl, kl_ref, tolerance = 1e-10)

    # Approaching q = 1 from both sides must converge to KL
    d_below <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 0.999,
        pseudocount = 0)
    d_above <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 1.001,
        pseudocount = 0)
    expect_equal(d_below, kl_ref, tolerance = 1e-3)
    expect_equal(d_above, kl_ref, tolerance = 1e-3)

    # q = 0.99 must NOT be treated as the KL special case (no ±0.01 band)
    d_099 <- TSENAT:::.tsallis_divergence_scalar(P * 100, Q * 100, q_val = 0.99,
        pseudocount = 0)
    expect_false(abs(d_099 - kl_ref) < 1e-6)
})

test_that("Asymmetry: D_q(P||Q) != D_q(Q||P) for q > 1", {
    P <- c(0.7, 0.2, 0.1) * 100
    Q <- c(0.4, 0.4, 0.2) * 100
    # Note: q = 0.5 is symmetric by construction (Bhattacharyya-type form),
    # so asymmetry is tested for q > 1 only.
    for (q in c(1.5, 2, 3)) {
        d_pq <- TSENAT:::.tsallis_divergence_scalar(P, Q, q_val = q, pseudocount = 0)
        d_qp <- TSENAT:::.tsallis_divergence_scalar(Q, P, q_val = q, pseudocount = 0)
        expect_true(abs(d_pq - d_qp) > 1e-10, info = paste("q =", q))
    }
})

test_that("Permutation invariance: joint reordering of categories", {
    P <- c(0.7, 0.2, 0.1) * 100
    Q <- c(0.4, 0.4, 0.2) * 100
    perm <- c(2, 3, 1)
    for (q in c(0.5, 1, 2)) {
        d1 <- TSENAT:::.tsallis_divergence_scalar(P, Q, q_val = q, pseudocount = 0)
        d2 <- TSENAT:::.tsallis_divergence_scalar(P[perm], Q[perm], q_val = q,
            pseudocount = 0)
        expect_equal(d1, d2, tolerance = 1e-10, info = paste("q =", q))
    }
})

test_that("Count-scale invariance: D(x, y) = D(c*x, c*y)", {
    x <- c(90, 10)
    y <- c(50, 50)
    for (q in c(0.5, 1, 2)) {
        d1 <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        d2 <- TSENAT:::.tsallis_divergence_scalar(10 * x, 1000 * y, q_val = q,
            pseudocount = 0)
        expect_equal(d1, d2, tolerance = 1e-10, info = paste("q =", q))
    }
})

test_that("q=0 convention: Q-mass on P's zero support", {
    # pseudocount = 0: supports differ
    x <- c(5, 5, 0)   # P = (0.5, 0.5, 0)
    y <- c(3, 0, 7)   # Q = (0.3, 0.0, 0.7)
    d0 <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 0, pseudocount = 0)
    expect_equal(d0, 0.7, tolerance = 1e-10)  # 1 - (0.3 + 0.0)

    # Identical supports -> 0
    d0_same <- TSENAT:::.tsallis_divergence_scalar(x, c(3, 7, 0), q_val = 0,
        pseudocount = 0)
    expect_equal(d0_same, 0, tolerance = 1e-10)

    # Under pseudocount regularization the q=0 limit is evaluated on the RAW
    # (pre-pseudocount) support (review Option A): D_0 keeps its
    # support-difference value instead of collapsing to 0.
    d0_pc <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 0, pseudocount = 0.5)
    expect_equal(d0_pc, 0.7, tolerance = 1e-12)
})

test_that("Support violation (P_i>0, Q_i=0) stays finite under regularization", {
    x <- c(90, 10, 0)
    y <- c(0, 100, 0)
    for (q in c(0.5, 1, 2)) {
        div_pc <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0.5)
        expect_true(is.finite(div_pc), info = paste("q =", q, "pseudocount=0.5"))
    }
})

test_that("Construct validity: identical isoform usage, different depth -> D=0", {
    # Control (90, 10) and treatment (900, 100): same isoform proportions
    x <- c(90, 10)
    y <- c(900, 100)
    for (q in c(0, 0.5, 1, 1.5, 2)) {
        d <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = q, pseudocount = 0)
        expect_equal(d, 0, tolerance = 1e-10, info = paste("q =", q))
    }
    # With the default pseudocount the divergence stays small
    d_pc <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 1, pseudocount = 0.5)
    expect_lt(d_pc, 0.01)
})

test_that("Pipeline: unequal group sizes produce valid divergence", {
    skip_if_not_installed("SummarizedExperiment")
    # 2 isoforms, 2 control + 3 treatment samples, same isoform usage
    counts <- rbind(
        iso1 = c(90, 90, 900, 900, 900),
        iso2 = c(10, 10, 100, 100, 100)
    )
    colnames(counts) <- paste0("S", 1:5)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(gene_name = c("GENE1", "GENE1")),
        colData = data.frame(group = c("A", "A", "B", "B", "B"),
            row.names = paste0("S", 1:5))
    )

    result <- TSENAT:::.calculate_divergence(
        se = se, group_col = "group", control_group = "A",
        q = 1, pseudocount = 0, norm = "none", progress = FALSE
    )

    est <- SummarizedExperiment::assay(result)[1, 1]
    expect_true(is.finite(est))
    expect_lt(est, 1e-6)  # identical usage -> ~0 divergence
})
