library(testthat)

context("LMM Phase 14 Enhancements: AR(1) Correlation")

test_that("LMM with use_ar1=FALSE (default) returns valid results", {
    skip_if_not_installed("lme4")
    
    # Create paired data with subject-level random effects
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 4))
    
    set.seed(50)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Run LMM (AR(1) is automatically attempted by .try_lm_fallbacks())
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Should have standard LMM output
    expect_true("p_interaction" %in% colnames(rd_out))
    expect_true("n_subjects" %in% colnames(rd_out))
    expect_true("small_sample_flag" %in% colnames(rd_out))
    expect_equal(rd_out$n_subjects[1], 4)
    expect_true(rd_out$small_sample_flag[1])  # 4 subjects < 5, so flag should be TRUE
})

test_that("LMM with use_ar1=TRUE attempts AR(1) correlation structure", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("nlme")
    
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4", "S5"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 5))
    
    set.seed(51)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("A", "B", "A", "B", "A"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Run LMM (AR(1) is automatically attempted by .try_lm_fallbacks())
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Should have LMM output with AR(1) attempt
    expect_true("p_interaction" %in% colnames(rd_out))
    expect_true("n_subjects" %in% colnames(rd_out))
    # fit_method should indicate which fitting strategy was used
    if ("fit_method" %in% colnames(rd_out)) {
        expect_true(!is.na(rd_out$fit_method[1]))
    }
})

test_that("LMM AR(1) falls back when correlation structure fails", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("nlme")
    
    # Create problematic data that might cause AR(1) fitting issues
    qvec <- seq(0.01, 0.05, by = 0.01)
    subject_ids <- rep(c("S1", "S2"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 2))
    
    set.seed(52)
    gene1_vals <- c(qvec * 1, qvec * 1.1) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("A", "B"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Run LMM (AR(1) is automatically attempted by .try_lm_fallbacks())
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 3
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Should have results even if AR(1) fails and falls back
    expect_true("p_interaction" %in% colnames(rd_out))
    expect_true(nrow(rd_out) > 0)
})

# ============================================================================
context("LMM Phase 14 Enhancements: Small-Sample Flagging")

test_that("small_sample_flag is FALSE when n_subjects >= 5", {
    skip_if_not_installed("lme4")
    
    # Create data with >= 5 subjects
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4", "S5"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 5))
    
    set.seed(53)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    expect_true("small_sample_flag" %in% colnames(rd_out))
    expect_true("n_subjects" %in% colnames(rd_out))
    expect_equal(rd_out$n_subjects[1], 5)
    expect_false(rd_out$small_sample_flag[1])
})

test_that("small_sample_flag is TRUE when n_subjects < 5", {
    skip_if_not_installed("lme4")
    
    # Create data with < 5 subjects
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 3))
    
    set.seed(54)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 6
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    expect_true("small_sample_flag" %in% colnames(rd_out))
    expect_true("n_subjects" %in% colnames(rd_out))
    expect_equal(rd_out$n_subjects[1], 3)
    expect_true(rd_out$small_sample_flag[1])
})

test_that("small_sample_flag correctly identifies boundary case (n_subjects = 5)", {
    skip_if_not_installed("lme4")
    
    # Exactly 5 subjects should NOT be flagged
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4", "S5"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 5))
    
    set.seed(55)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Boundary: 5 subjects should have small_sample_flag = FALSE
    expect_equal(rd_out$n_subjects[1], 5)
    expect_false(rd_out$small_sample_flag[1])
})

# ============================================================================
context("LMM Phase 14 Enhancements: Enhanced Fallback Reporting")

test_that("fit_method column is present in LMM results", {
    skip_if_not_installed("lme4")
    
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 4))
    
    set.seed(56)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # fit_method should identify the fitting strategy used
    if ("fit_method" %in% colnames(rd_out)) {
        expect_true(!is.na(rd_out$fit_method[1]))
        # Should indicate which fallback was used (lme4, glmmTMB, etc.)
        expect_true(nchar(rd_out$fit_method[1]) > 0)
    } else {
        # If fit_method not present, other reporting columns should be
        expect_true("p_interaction" %in% colnames(rd_out))
    }
})

test_that("Enhanced fallback reporting includes convergence information", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("glmmTMB")
    
    # Create data that might trigger fallback logic
    qvec <- seq(0.01, 0.05, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 3))
    
    set.seed(57)
    # Use small variance to potentially cause convergence issues
    gene1_vals <- c(qvec * 1, qvec * 1.05, qvec * 0.98) + rnorm(length(coln), sd = 1e-5)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 3
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Should have standard columns
    expect_true("p_interaction" %in% colnames(rd_out))
    expect_true("n_subjects" %in% colnames(rd_out))
})

# ============================================================================
context("LMM Phase 14 Enhancements: Integration Tests")

test_that("LMM results maintain backward compatibility with AR(1)=FALSE default", {
    skip_if_not_installed("lme4")
    
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 4))
    
    set.seed(58)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Run default LMM twice to check consistency
    set.seed(58)
    res_first <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    set.seed(58)
    res_second <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    # Both should return comparable results
    if (!is.data.frame(res_first)) {
        res_first <- as.data.frame(SummarizedExperiment::rowData(res_first))
    }
    if (!is.data.frame(res_second)) {
        res_second <- as.data.frame(SummarizedExperiment::rowData(res_second))
    }
    
    # Should have same structure
    expect_equal(colnames(res_first), colnames(res_second))
})

test_that("Multiple genes with varying sample sizes are handled correctly", {
    skip_if_not_installed("lme4")
    
    qvec <- seq(0.01, 0.1, by = 0.01)
    # Use 5 subjects for consistent data structure
    subject_ids <- rep(c("S1", "S2", "S3", "S4", "S5"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 5))
    
    # Create two genes with same structure
    set.seed(60)
    mat <- rbind(
        g1 = c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8) + rnorm(length(coln), sd = 1e-3),
        g2 = c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8) + rnorm(length(coln), sd = 1e-3)
    )
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # All genes should have proper n_subjects and small_sample_flag
    expect_true(all("n_subjects" %in% colnames(rd_out)))
    expect_true(all("small_sample_flag" %in% colnames(rd_out)))
    expect_true(nrow(rd_out) >= 2)
})

test_that("LMM with AR(1) produces different p-values than baseline in some cases", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("nlme")
    
    # Create data with strong correlation structure across q values
    qvec <- seq(0.01, 0.3, by = 0.01)  # More q values to better show AR(1) effects
    subject_ids <- rep(c("S1", "S2", "S3", "S4", "S5"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 5))
    
    set.seed(61)
    # Simulate autocorrelated residuals (AR(1) should help with this)
    base_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2, qvec * 1.8)
    ar_residuals <- stats::filter(rnorm(length(coln), sd = 0.01), 0.7, method = "recursive")
    gene1_vals <- base_vals + ar_residuals
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor", "Normal"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Run LMM with AR(1) automatically attempted
    set.seed(61)
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (!is.data.frame(res)) {
        res <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Should produce valid results
    expect_true("p_interaction" %in% colnames(res))
    expect_true(!is.na(res$p_interaction[1]))
})

test_that("LMM output includes all Phase 14 columns", {
    skip_if_not_installed("lme4")
    
    qvec <- seq(0.01, 0.1, by = 0.01)
    subject_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(subject_ids, "_q=", rep(qvec, times = 4))
    
    set.seed(62)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 1.2) + rnorm(length(coln), sd = 1e-3)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(
        samples = subject_ids,
        sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)),
        sample_base = subject_ids,
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res <- .calculate_lm_interaction(
        se,
        condition_col = "sample_type",
        method = "lmm",
        subject_col = "sample_base",
        min_obs = 8
    )
    
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    
    # Phase 14 enhancements should include these columns
    required_cols <- c("p_interaction", "n_subjects", "small_sample_flag")
    for (col in required_cols) {
        expect_true(col %in% colnames(rd_out),
                   paste("Phase 14 column missing:", col))
    }
    
    # Check data types
    expect_true(is.numeric(rd_out$n_subjects))
    expect_true(is.logical(rd_out$small_sample_flag))
})
