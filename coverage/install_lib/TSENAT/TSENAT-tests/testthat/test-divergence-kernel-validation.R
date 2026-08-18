# ============================================================================
# test-divergence-kernel-validation.R
#
# "Definitive test" requested in div_review.md (sections 8-11): validate the
# numerical kernel independently of the implementation.
#
#  1. Manual D_q recomputation (pure R, from the EXACT P/R vectors that reach
#     the kernel) vs .tsallis_divergence_scalar (pipeline) and
#     tsallis_divergence_cpp (C++ kernel) for q = 0, 0.25, ..., 2
#  2. Closed-form identities: q = 1 -> KL, q = 2 -> sum P^2/R - 1
#  3. Forward/reverse spectra must differ (Tsallis is not symmetric)
#  4. Finiteness sanity checks before the kernel (review section 11)
#  5. Documented support-violation sentinel (1e10, bootstrap safeguard)
#  6. End-to-end: D_q stored by the pipeline == independent recomputation
#     from the raw per-condition isoform counts
# ============================================================================

# Independent pure-R reference implementation (no shared helpers)
.dq_manual <- function(p, r, q, log_base = exp(1)) {
    if (abs(q) < 1e-10) {
        # Limit convention: D_0 = 1 - sum_{i: p_i > 0} r_i
        return(1 - sum(r[p > 1e-15]))
    }
    if (abs(q - 1) < 1e-10) {
        # KL divergence (support violation -> Inf, mirrored by 1e10 sentinel)
        if (any(p > 1e-15 & r <= 1e-15)) return(1e10)
        return(sum(p[p > 1e-15] * log(p[p > 1e-15] / r[p > 1e-15]) / log(log_base)))
    }
    if (q > 0) {
        if (any(p > 1e-15 & r <= 1e-15)) return(1e10)
        return((sum(p^q * r^(1 - q)) - 1) / (q - 1))
    }
    NA_real_
}

# Pipeline-equivalent normalization (pseudocount formula)
.normalize_counts <- function(x, pseudocount) {
    (x + pseudocount) / (sum(x) + length(x) * pseudocount)
}

test_that("Kernel: manual D_q matches pipeline scalar and C++ for a q grid", {
    # Asymmetric, multi-category counts with heterogeneous log-ratios
    # (mimics the transcriptomic heterogeneity discussed in div_review.md)
    x <- c(412, 89, 15, 3, 0)
    y <- c(201, 180, 44, 12, 2)
    q_grid <- seq(0, 2, by = 0.25)

    # Raw (pre-pseudocount) normalized vectors: q = 0 is evaluated on the RAW
    # support (review_divergence.md, Option A), so the reference computation
    # for q = 0 must use pseudocount-free probabilities.
    p0 <- .normalize_counts(x, 0)
    r0 <- .normalize_counts(y, 0)

    for (pc in c(0.5, 0)) {
        p <- .normalize_counts(x, pc)
        r <- .normalize_counts(y, pc)

        # Sanity checks BEFORE the kernel (review section 11)
        expect_true(all(is.finite(p)) && all(is.finite(r)))
        expect_true(min(p) >= 0 && min(r) >= 0)

        for (qv in q_grid) {
            p_use <- if (qv == 0) p0 else p
            r_use <- if (qv == 0) r0 else r

            d_manual <- .dq_manual(p_use, r_use, qv)
            d_scalar <- TSENAT:::.tsallis_divergence_scalar(x, y, q_val = qv,
                pseudocount = pc)
            d_cpp <- TSENAT:::tsallis_divergence_cpp(p_use, r_use, q = qv)

            expect_true(is.finite(d_manual), info = paste("pc", pc, "q", qv))
            # Review section 9: all.equal(tolerance = 1e-10). The C++ kernel
            # must match the independent recomputation exactly.
            expect_equal(d_cpp, d_manual, tolerance = 1e-10,
                info = paste("cpp pc", pc, "q", qv))
            # The pipeline scalar applies a documented min_prob = 1e-10 clamp
            # (with renormalization) ONLY when pseudocount == 0, as a safety
            # net against log(0); with pc > 0 it matches exactly. At q = 0
            # both kernels evaluate the raw support regardless of pc.
            tol_scalar <- if (pc < 1e-10) 1e-04 else 1e-10
            expect_equal(d_scalar, d_manual, tolerance = tol_scalar,
                info = paste("scalar pc", pc, "q", qv))
        }
    }
})

test_that("Kernel: q = 1 equals KL and q = 2 equals sum P^2/R - 1", {
    x <- c(70, 20, 10)
    y <- c(40, 40, 20)
    p <- .normalize_counts(x, 0)
    r <- .normalize_counts(y, 0)

    kl_ref <- sum(p * log(p / r))
    chi2_ref <- sum(p^2 / r) - 1

    expect_equal(TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 1,
        pseudocount = 0), kl_ref, tolerance = 1e-10)
    expect_equal(TSENAT:::tsallis_divergence_cpp(p, r, q = 1), kl_ref,
        tolerance = 1e-10)
    expect_equal(TSENAT:::.tsallis_divergence_scalar(x, y, q_val = 2,
        pseudocount = 0), chi2_ref, tolerance = 1e-10)
    expect_equal(TSENAT:::tsallis_divergence_cpp(p, r, q = 2), chi2_ref,
        tolerance = 1e-10)
})

test_that("Kernel: forward and reverse spectra differ (asymmetry)", {
    x <- c(412, 89, 15, 3, 0)
    y <- c(201, 180, 44, 12, 2)
    p <- .normalize_counts(x, 0.5)
    r <- .normalize_counts(y, 0.5)

    q_grid <- seq(0.25, 2, by = 0.25)
    d_fwd <- vapply(q_grid, function(qv) TSENAT:::tsallis_divergence_cpp(p, r,
        q = qv), numeric(1))
    d_rev <- vapply(q_grid, function(qv) TSENAT:::tsallis_divergence_cpp(r, p,
        q = qv), numeric(1))

    # In general the curves must differ (at least at some q; KL included).
    # q < 1 may be close by construction, but KL (q = 1) and q = 2 must differ.
    expect_true(abs(d_fwd[q_grid == 1] - d_rev[q_grid == 1]) > 1e-8)
    expect_true(abs(d_fwd[q_grid == 2] - d_rev[q_grid == 2]) > 1e-8)
})

test_that("Kernel: support violation returns the documented 1e10 sentinel", {
    p <- c(0.7, 0.2, 0.1)
    r <- c(0.5, 0.5, 0.0)  # category 3 missing in r

    expect_equal(TSENAT:::tsallis_divergence_cpp(p, r, q = 1), 1e10)
    expect_equal(TSENAT:::tsallis_divergence_cpp(p, r, q = 1.5), 1e10)
    # q = 0 limit convention is finite: p has full support here, so the
    # r-mass on p's support is 1 and D_0 = 1 - 1 = 0
    expect_equal(TSENAT:::tsallis_divergence_cpp(p, r, q = 0), 0)
})

test_that("End-to-end: stored D_q equals independent recomputation from raw counts", {
    # Two genes (3 isoforms each) and 2 conditions: one sample per condition.
    # The pipeline aggregates per-gene per-condition isoform counts, so the
    # per-sample distribution is irrelevant for the point estimate.
    set.seed(42)
    counts <- matrix(c(
        400, 350,   # gene A isoform 1: control, treatment
        120, 100,   # gene A isoform 2
         30,  25,   # gene A isoform 3
        180, 150,   # gene B isoform 1
        200, 220,   # gene B isoform 2
         60,  70    # gene B isoform 3
    ), nrow = 6, byrow = TRUE)
    colnames(counts) <- c("ctrl", "trt")
    rownames(counts) <- paste0("tx", 1:6)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(
            gene_id = c(rep("geneA", 3), rep("geneB", 3)),
            row.names = rownames(counts)
        ),
        colData = S4Vectors::DataFrame(
            condition = c("control", "treatment"),
            row.names = colnames(counts)
        )
    )

    q_vals <- c(0.5, 1, 2)
    res <- TSENAT:::.calculate_divergence_impl(
        se, group_col = "condition", control_group = "control",
        q = q_vals, bootstrap = FALSE, pseudocount = 0.5,
        nthreads = 1, progress = FALSE
    )

    expect_s4_class(res, "SummarizedExperiment")
    est <- SummarizedExperiment::assay(res, "divergence")
    expect_equal(rownames(est), c("geneA", "geneB"))

    # Independent recomputation: P = per-condition isoform fractions
    for (g in c("geneA", "geneB")) {
        idx <- which(SummarizedExperiment::rowData(se)$gene_id == g)
        x_counts <- counts[idx, "ctrl"]
        y_counts <- counts[idx, "trt"]
        p <- .normalize_counts(x_counts, 0.5)
        r <- .normalize_counts(y_counts, 0.5)

        d_manual <- vapply(q_vals, function(qv) .dq_manual(p, r, qv), numeric(1))
        expect_equal(as.numeric(est[g, ]), d_manual, tolerance = 1e-8,
            info = paste("gene", g))
    }
})
