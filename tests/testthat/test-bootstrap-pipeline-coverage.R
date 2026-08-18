context("Bootstrap pipeline: additional branch coverage")

# ============================================================================
# .bootstrap_validate_inputs — warning branches for q ranges
# ============================================================================

test_that(".bootstrap_validate_inputs warns on very small q", {
    expect_warning(
        TSENAT:::.bootstrap_validate_inputs(
            x = c(10, 5, 2),
            q = 0.0005,
            nboot = 100,
            ci = 0.95,
            paired = FALSE
        ),
        "very small q"
    )
})

test_that(".bootstrap_validate_inputs warns on extreme q > 100", {
    expect_warning(
        TSENAT:::.bootstrap_validate_inputs(
            x = c(10, 5, 2),
            q = 101,
            nboot = 100,
            ci = 0.95,
            paired = FALSE
        ),
        "extreme q values > 100"
    )
})

# ============================================================================
# .bootstrap_process_matrix — 0-row error branch
# ============================================================================

test_that(".bootstrap_process_matrix errors on 0-row matrix", {
    mat <- matrix(numeric(0), nrow = 0, ncol = 3)
    expect_error(
        TSENAT:::.bootstrap_process_matrix(
            x = mat,
            q = 1,
            norm = FALSE,
            nboot = 10,
            ci = 0.95,
            method = "percentile",
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            gene_name = NULL,
            verbose = FALSE,
            include_diagnostics = FALSE,
            use_job = FALSE,
            nthreads = 1,
            paired = FALSE
        ),
        "0 rows"
    )
})

# ============================================================================
# .bootstrap_process_se — input validation branches
# ============================================================================

test_that(".bootstrap_process_se rejects non-SummarizedExperiment", {
    expect_error(
        TSENAT:::.bootstrap_process_se(
            se = list(),
            res = data.frame(gene_id = "g1"),
            top_n = 1,
            q = 1,
            norm = FALSE,
            nboot = 10,
            ci = 0.95,
            method = "percentile",
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            gene_name = NULL,
            verbose = FALSE,
            include_diagnostics = FALSE,
            use_job = FALSE,
            paired = FALSE
        ),
        "'se' must be SummarizedExperiment"
    )
})

test_that(".bootstrap_process_se rejects non-data.frame res", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2))
    )
    expect_error(
        TSENAT:::.bootstrap_process_se(
            se = se,
            res = "not_a_data_frame",
            top_n = 1,
            q = 1,
            norm = FALSE,
            nboot = 10,
            ci = 0.95,
            method = "percentile",
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            gene_name = NULL,
            verbose = FALSE,
            include_diagnostics = FALSE,
            use_job = FALSE,
            paired = FALSE
        ),
        "'res' must be data.frame"
    )
})

test_that(".bootstrap_process_se rejects res without gene_id", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2))
    )
    expect_error(
        TSENAT:::.bootstrap_process_se(
            se = se,
            res = data.frame(other = "g1"),
            top_n = 1,
            q = 1,
            norm = FALSE,
            nboot = 10,
            ci = 0.95,
            method = "percentile",
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            gene_name = NULL,
            verbose = FALSE,
            include_diagnostics = FALSE,
            use_job = FALSE,
            paired = FALSE
        ),
        "gene_id"
    )
})

# ============================================================================
# .bootstrap_process_multiple_q — fast Hill-number (what = "D") path
# ============================================================================

test_that(".bootstrap_process_multiple_q fast path handles Hill numbers", {
    x <- c(100, 50, 25, 10, 5, 2)
    res <- TSENAT:::.bootstrap_process_multiple_q(
        x = x,
        q = c(1, 2, 3),
        norm = FALSE,
        nboot = 10,
        ci = 0.95,
        method = "percentile",
        log_base = exp(1),
        pseudocount = 0,
        what = "D",
        gene_name = NULL,
        verbose = FALSE,
        include_diagnostics = FALSE,
        use_job = FALSE,
        paired = FALSE
    )
    expect_s3_class(res, "tsenat_bootstrap_ci_list")
    expect_equal(names(res), paste0("q=", c(1, 2, 3)))
    # q = 1 uses the exp(H) conversion; q > 1 uses the Hill-number power form
    expect_true(all(vapply(res, function(x) is.finite(x$estimate), logical(1))))
})

# ============================================================================
# .bootstrap_resample_with_quality_control — CRITICAL stop + late return
# ============================================================================

test_that(".bootstrap_resample_with_quality_control stops when valid_frac < 0.5", {
    # Tiny values make replicate resampling produce all-NA replicates, which
    # regeneration cannot repair -> valid_frac stays 0 -> CRITICAL stop.
    expect_error(
        TSENAT:::.bootstrap_resample_with_quality_control(
            x = c(1e-11, 1e-11),
            q = 1,
            norm = FALSE,
            nboot = 10L,
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            paired = FALSE,
            effective_length = NULL,
            min_valid_frac = 0.75,
            resample_by = "replicate"
        ),
        "CRITICAL"
    )
})

test_that(".bootstrap_resample_with_quality_control returns with partial NA when threshold already met", {
    # Some resamples are all-tiny (NA) but the valid fraction already meets the
    # threshold, so the function returns the raw distribution (final return).
    set.seed(1)
    res <- TSENAT:::.bootstrap_resample_with_quality_control(
        x = c(1e-11, 1e-11, 1),
        q = 1,
        norm = FALSE,
        nboot = 30L,
        log_base = exp(1),
        pseudocount = 0,
        what = "S",
        paired = FALSE,
        effective_length = NULL,
        min_valid_frac = 0.6,
        resample_by = "replicate"
    )
    expect_is(res, "numeric")
    expect_length(res, 30)
    # Some NAs may remain since no regeneration was needed
    expect_true(sum(!is.na(res)) / length(res) >= 0.6)
})

# ============================================================================
# .bootstrap_compute_ci — replicate counts_matrix branch + all-NA branch
# ============================================================================

test_that(".bootstrap_compute_ci applies effective_length to counts_matrix", {
    counts_matrix <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, ncol = 3)
    effective_length <- c(100, 200)
    res <- TSENAT:::.bootstrap_compute_ci(
        x = c(30, 70),
        q = 1,
        norm = FALSE,
        nboot = 8,
        ci = 0.95,
        method = "percentile",
        log_base = exp(1),
        pseudocount = 1,
        what = "S",
        paired = FALSE,
        effective_length = effective_length,
        min_valid_frac = 0.5,
        resample_by = "replicate",
        counts_matrix = counts_matrix
    )
    expect_true(is.finite(res$point_est))
    expect_true(is.finite(res$ci_result$lower))
    expect_true(is.finite(res$ci_result$upper))
})

test_that(".bootstrap_compute_ci warns and returns NA CI when all replicates invalid", {
    # min_valid_frac = 0 lets the QC step return the all-NA distribution, which
    # triggers the "all replicates produced NA/NaN" branch.
    expect_warning(
        res <- TSENAT:::.bootstrap_compute_ci(
            x = c(1e-11, 1e-11),
            q = 1,
            norm = FALSE,
            nboot = 10,
            ci = 0.95,
            method = "percentile",
            log_base = exp(1),
            pseudocount = 0,
            what = "S",
            paired = FALSE,
            effective_length = NULL,
            min_valid_frac = 0,
            resample_by = "replicate"
        ),
        "All bootstrap replicates produced NA/NaN"
    )
    expect_true(is.na(res$ci_result$lower))
    expect_true(is.na(res$ci_result$upper))
})
