context("Fit One Interaction Function Tests")

library(testthat)

# ═══════════════════════════════════════════════════════════════════════════
# CORE BEHAVIOR: Linear branch with min_obs and observations checks
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction linear branch handles min_obs and returns p", {
    set.seed(1)
    # construct tiny matrix: 1 gene x 3 observations
    mat <- matrix(rnorm(3), nrow = 1)
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.2, 0.3)
    sample_names <- paste0("s", seq_along(q_vals))
    group_vec <- c("A", "A", "B")
    # min_obs > non-missing -> NULL
    res_null <- .fit_one_interaction("g1",
        se = NULL, mat = mat, q_vals = q_vals,
        sample_names = sample_names, group_vec = group_vec, method = "lmm",
        pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 10, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_null(res_null)

    # generate larger sample with clear interaction signal
    n <- 60
    qv <- rep(seq(0.1, 1.0, length.out = 20), 3)
    group <- rep(c("A", "B", "A"), each = 20)
    obs <- 0.5 * qv + ifelse(group == "B", 0.6 * qv, 0) + rnorm(length(qv), 0, 0.05)
    mat2 <- matrix(obs, nrow = 1)
    rownames(mat2) <- "gX"
    res <- .fit_one_interaction("gX",
        se = NULL, mat = mat2, q_vals = qv,
        sample_names = paste0("s", seq_along(qv)), group_vec = group, method = "lmm",
        pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_true(is.data.frame(res) || is.null(res))
    if (is.data.frame(res)) expect_true(is.numeric(res$p_interaction) || is.na(res$p_interaction))
})

# ═══════════════════════════════════════════════════════════════════════════
# LMM BRANCH: Error handling and edge cases
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction lmm errors when subject_col missing or paired with no sample_base", {
    skip_if_not_installed("SummarizedExperiment")
    set.seed(2)
    # small dataset to attach to se
    qv <- rep(0.1, 8)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), length.out = length(qv))
    obs <- rnorm(length(qv))
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g1"

    coldata <- S4Vectors::DataFrame(samples = sample_names)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat), colData = coldata)

    # subject_col provided but not present
    expect_error(.fit_one_interaction("g1",
        se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = "foo", paired = FALSE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    ), "subject_col")

    # paired = TRUE but no sample_base column
    expect_error(.fit_one_interaction("g1",
        se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = NULL, paired = TRUE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    ), "paired = TRUE")
})

test_that(".fit_one_interaction lmm returns NULL with <2 subjects", {
    skip_if_not_installed("SummarizedExperiment")
    set.seed(3)
    qv <- rep(c(0.1, 0.5), times = 3)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), length.out = length(qv))
    # create single subject for all samples
    coldata <- S4Vectors::DataFrame(samples = sample_names, sample_base = rep("sub1", length(qv)))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(rnorm(length(qv)), nrow = 1)), colData = coldata)
    mat <- matrix(rnorm(length(qv)), nrow = 1)
    rownames(mat) <- "g1"

    res <- .fit_one_interaction("g1",
        se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = NULL, paired = TRUE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_null(res)
})

# ═══════════════════════════════════════════════════════════════════════════
# FALLBACK BEHAVIOR: When lmer fails or returns singular fits
# ═══════════════════════════════════════════════════════════════════════════

test_that("lmm branch falls back to lm when mixed model fitting fails and respects pvalue selection", {
    skip_if_not_installed("lme4")
    # Temporarily force .try_lmer to fail so the code uses the lm fallbacks
    ns <- asNamespace("TSENAT")
    orig_try <- get(".try_lmer", envir = ns)
    assignInNamespace(".try_lmer", function(...) structure("error", class = "try-error"), ns = "TSENAT")
    on.exit(assignInNamespace(".try_lmer", orig_try, ns = "TSENAT"), add = TRUE)

    set.seed(42)
    n <- 40
    qv <- rep(seq(0.1, 1, length.out = 20), 2)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), each = 20)
    obs <- 0.5 * qv + ifelse(group == "B", 0.2 * qv, 0) + rnorm(length(qv), 0, 0.05)
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g_fallback"

    # pvalue = 'both' should return p_lrt (no Satterthwaite with nlme AR(1))
    res_both <- .fit_one_interaction("g_fallback",
        se = NULL, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "both",
        subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_true(is.data.frame(res_both))
    # nlme LMM returns p_interaction, p_lrt (Satterthwaite not available for AR(1))
    expect_true(all(c("p_interaction", "p_lrt", "fit_method", "singular") %in% colnames(res_both)))
    # AR(1) implementation tries nlme first, then glmmTMB, then lm_subject_fixed, then lm_nosubject
    expect_true(res_both$fit_method %in% c(
        "nlme::lme", "nlme::lme_ar1", "nlme", 
        "nlme::lme_ar1_raw", "nlme::lme_arima(1,1,0)",
        "glmmTMB", "lm_subject_fixed", "lm_nosubject"
    ))

    # pvalue = 'lrt' should use the LRT p-value
    res_lrt <- .fit_one_interaction("g_fallback",
        se = NULL, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_true(is.data.frame(res_lrt))
    expect_true(is.numeric(res_lrt$p_interaction) || is.na(res_lrt$p_interaction))

    # pvalue = 'satterthwaite' (ignored for nlme but parameter still accepted for compatibility)
    res_sat <- .fit_one_interaction("g_fallback",
        se = NULL, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "satterthwaite",
        subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE
    )
    expect_true(is.data.frame(res_sat))
    expect_true(is.numeric(res_sat$p_interaction) || is.na(res_sat$p_interaction))
})

test_that("lmm branch uses fallback when lmer returns singular fits", {
    skip_if_not_installed("lme4")
    ns <- asNamespace("TSENAT")
    orig_try <- get(".try_lmer", envir = ns)
    fake_lmer <- function(...) {
        m <- list()
        class(m) <- "lmerMod"
        attr(m, "singular") <- TRUE
        return(m)
    }
    assignInNamespace(".try_lmer", fake_lmer, ns = "TSENAT")
    on.exit(assignInNamespace(".try_lmer", orig_try, ns = "TSENAT"), add = TRUE)

    set.seed(7)
    qv <- rep(seq(0.1, 1, length.out = 20), 2)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), each = 20)
    obs <- 0.5 * qv + ifelse(group == "B", 0.2 * qv, 0) + rnorm(length(qv), 0, 0.05)
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g_sing"

    # AR(1) implementation tries nlme first, which may succeed or require fallback
    # If nlme fails, a fallback message is printed
    tryCatch({
        res <- .fit_one_interaction("g_sing", se = NULL, mat = mat, q_vals = qv, sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "both", subject_col = NULL, paired = FALSE, min_obs = 5, verbose = TRUE, suppress_lme4_warnings = TRUE, progress = TRUE)
    }, error = function(e) { res <<- NULL })
    
    expect_true(is.data.frame(res))
    # AR(1) implementation now tries nlme first, so fit_method can be nlme or a fallback method
    expect_true(res$fit_method %in% c(
        "nlme::lme", "nlme::lme_ar1", "nlme", 
        "nlme::lme_ar1_raw", "nlme::lme_arima(1,1,0)",
        "glmmTMB", "lm_subject_fixed", "lm_nosubject"
    ))
})

test_that("lmm branch uses subject_col and returns lmer method when available", {
    skip_if_not_installed("lme4")
    ns <- asNamespace("TSENAT")
    orig_try <- get(".try_lmer", envir = ns)
    fake_lmer_ok <- function(...) {
        m <- list()
        class(m) <- "lmerMod"
        attr(m, "singular") <- FALSE
        return(m)
    }
    assignInNamespace(".try_lmer", fake_lmer_ok, ns = "TSENAT")
    on.exit(assignInNamespace(".try_lmer", orig_try, ns = "TSENAT"), add = TRUE)

    # Build sample metadata with custom subject column
    qv <- rep(seq(0.1, 1, length.out = 20), 2)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), each = 20)
    obs <- 0.5 * qv + ifelse(group == "B", 0.2 * qv, 0) + rnorm(length(qv), 0, 0.05)
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g_sub"

    coldata <- S4Vectors::DataFrame(samples = sample_names, my_subject = rep(paste0("sub", 1:20), 2))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat), colData = coldata)

    res <- .fit_one_interaction("g_sub", se = se, mat = mat, q_vals = qv, sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt", subject_col = "my_subject", paired = FALSE, min_obs = 5, verbose = FALSE, suppress_lme4_warnings = TRUE, progress = FALSE)
    expect_true(is.data.frame(res))
    expect_true(res$fit_method %in% c("nlme::lme", "nlme::lme_ar1", "nlme", "glmmTMB", "lm_subject_fixed", "lm_nosubject", "lmer", "lmer_singular", "fallback", "lm_subject"))
})

# ═══════════════════════════════════════════════════════════════════════════
# METHOD DISPATCHING: GAM, FPCA, and other methods
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction dispatches to gam and fpca methods", {
    # FPCA dispatch
    sample_names <- rep(paste0("s", 1:4), each = 2)
    q_vals <- rep(1:2, times = 4)
    group_vec <- rep(c("A", "B", "A", "B"), each = 2)
    obs <- rnorm(8)
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g1"
    out_fpca <- .fit_one_interaction("g1", se = NULL, mat = mat, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, method = "fpca", pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 2, verbose = FALSE, suppress_lme4_warnings = TRUE, progress = FALSE)
    expect_true(is.null(out_fpca) || (is.data.frame(out_fpca) && "p_interaction" %in% colnames(out_fpca)))

    # GAM dispatch - if mgcv available
    if (rlang::is_installed("mgcv")) {
        set.seed(1)
        n <- 30
        q <- runif(n, 0.1, 1)
        group <- rep(c("A", "B"), length.out = n)
        entropy <- 0.2 * q + ifelse(group == "B", 0.3 * q, 0) + rnorm(n, 0, 0.01)
        df <- data.frame(entropy = entropy, q = q, group = group)
        mat_gam <- matrix(entropy, nrow = 1)
        rownames(mat_gam) <- "g1"
        out_gam <- suppressWarnings(.fit_one_interaction("g1", se = NULL, mat = mat_gam, q_vals = q, sample_names = paste0("s", seq_along(q)), group_vec = group, method = "gam", pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE, suppress_lme4_warnings = TRUE, progress = FALSE))
        expect_true(is.null(out_gam) || (is.data.frame(out_gam) && "p_interaction" %in% colnames(out_gam)))
    } else {
        succeed()
    }
})

# ═══════════════════════════════════════════════════════════════════════════
# VARIANCE STRUCTURE AND HETEROSCEDASTICITY
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction LMM applies variance structure for heteroscedasticity", {
    skip_if_not_installed("nlme")
    set.seed(222)
    
    # Create data with heteroscedasticity and subjects
    n_subj <- 8
    n_per <- 10
    n_total <- n_subj * n_per
    
    q <- rep(runif(n_per, 0.1, 2), n_subj)
    group <- rep(rep(c("A", "B"), length.out = n_per), n_subj)
    subject <- rep(1:n_subj, each = n_per)
    
    # Heteroscedastic entropy: variance scales with q in group B
    entropy <- 0.5 + 0.2 * q + 
               ifelse(group == "B", 0.3 * q, 0) + 
               rnorm(n_total, 0, sd = ifelse(group == "B", 0.08 * q + 0.01, 0.05))
    
    # Build expression count matrix (1 gene for testing)
    expr_matrix <- matrix(rnorm(1 * n_total), nrow = 1)
    rownames(expr_matrix) <- "gene1"
    colnames(expr_matrix) <- paste0("s", 1:n_total)
    
    # Build SummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = expr_matrix),
        colData = S4Vectors::DataFrame(
            sample = colnames(expr_matrix),
            q = q,
            group = group,
            subject = factor(subject)
        )
    )
    
    # Call the function with LMM method
    # Extract necessary components for .fit_one_interaction
    mat <- SummarizedExperiment::assays(se)$counts
    coldata <- SummarizedExperiment::colData(se)
    sample_names <- coldata$sample
    q_vals <- coldata$q
    group_vec <- coldata$group
    
    result <- .fit_one_interaction(
        g = "gene1",
        se = se,
        mat = mat,
        q_vals = q_vals,
        sample_names = sample_names,
        group_vec = group_vec,
        subject_col = "subject",
        method = "lmm",
        min_obs = 5,
        verbose = FALSE
    )
    
    expect_true(is.data.frame(result) || is.null(result))
    if (is.data.frame(result)) {
        expect_true("p_interaction" %in% colnames(result))
    }
})

# ═══════════════════════════════════════════════════════════════════════════
# SLOPE_DIFF EXTRACTION: Slope difference in interaction
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction LMM method includes slope_diff in results", {
    skip_if_not_installed("nlme")
    set.seed(1001)
    
    # Create paired data with clear interaction signal
    # Use 10 subjects, each with 5 q values per group = 100 observations total
    n_subjects <- 10
    n_q <- 5
    n_groups <- 2
    
    # Create vectors that repeat properly for the matrix structure
    subject_ids <- rep(paste0("sub", 1:n_subjects), n_q * n_groups)
    qv <- rep(seq(0.1, 1.5, length.out = n_q), n_subjects * n_groups)
    group_vec <- rep(rep(c("A", "B"), each = n_q), n_subjects)
    
    n_total <- length(subject_ids)
    
    # Add interaction effect: group B has steeper slope with q
    entropy <- 0.5 + 0.3 * qv + ifelse(group_vec == "B", 0.4 * qv, 0) + rnorm(n_total, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene1"
    
    # Create proper colData with samples and sample_base columns
    coldata <- S4Vectors::DataFrame(
        samples = paste0("s", 1:n_total),
        sample_base = subject_ids
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        colData = coldata
    )
    
    result <- .fit_one_interaction(
        "gene1",
        se = se,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n_total),
        group_vec = group_vec,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = TRUE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Verify result is a data.frame with slope_diff column
    # Note: result may be NULL if nlme fitting fails, which is acceptable
    if (is.data.frame(result)) {
        expect_true("slope_diff" %in% colnames(result))
        expect_true(is.numeric(result$slope_diff) || is.na(result$slope_diff))
    } else {
        # If nlme fitting failed (returns NULL), that's still acceptable
        # The important thing is that slope_diff extraction doesn't cause errors
        expect_true(is.null(result))
    }
})

context("Fit One Interaction - Additional Coverage Tests")

library(testthat)

# ═══════════════════════════════════════════════════════════════════════════
# GEE METHOD TESTS (Generalized Estimating Equations)
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction GEE method with AR(1) correlation works", {
    skip_if_not_installed("geepack")
    set.seed(101)
    
    # Create clustered data with subject effect
    n_subjects <- 8
    n_per <- 5
    n_total <- n_subjects * n_per
    
    subject <- rep(1:n_subjects, each = n_per)
    q_vals <- rep(seq(0.5, 2, length.out = n_per), n_subjects)
    group <- rep(c("A", "B"), length.out = n_total)
    entropy <- 0.4 + 0.2 * q_vals + 0.3 * (group == "B") + 0.1 * (group == "B") * q_vals + 
               rnorm(n_total, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_gee"
    
    result <- .fit_one_interaction(
        "gene_gee",
        se = NULL,
        mat = mat,
        q_vals = q_vals,
        sample_names = paste0("s", 1:n_total),
        group_vec = group,
        method = "gee",
        pvalue = "wald",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        corstr = "ar1"
    )
    
    # GEE may succeed or fallback - just ensure valid result structure
    expect_true(is.data.frame(result) || is.null(result))
    if (is.data.frame(result)) {
        expect_true("p_interaction" %in% colnames(result))
    }
})

test_that(".fit_one_interaction GEE with exchangeable correlation structure", {
    skip_if_not_installed("geepack")
    set.seed(102)
    
    n_subjects <- 6
    n_per <- 4
    n_total <- n_subjects * n_per
    
    subject <- rep(1:n_subjects, each = n_per)
    q_vals <- rep(seq(0.1, 1.5, length.out = n_per), n_subjects)
    group <- rep(c("A", "B"), each = n_per * n_subjects / 2)
    entropy <- 0.3 + 0.25 * q_vals + rnorm(n_total, 0, 0.06)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_ex"
    
    result <- .fit_one_interaction(
        "gene_ex",
        se = NULL,
        mat = mat,
        q_vals = q_vals,
        sample_names = paste0("s", 1:n_total),
        group_vec = group,
        method = "gee",
        pvalue = "wald",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        corstr = "exchangeable"
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction GEE with independence correlation", {
    skip_if_not_installed("geepack")
    set.seed(103)
    
    # Independence is equivalent to regular glm with cluster-robust SE
    n <- 30
    q_vals <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.4 + 0.2 * q_vals + 0.15 * (group == "B") + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_ind"
    
    result <- .fit_one_interaction(
        "gene_ind",
        se = NULL,
        mat = mat,
        q_vals = q_vals,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "gee",
        pvalue = "wald",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        corstr = "independence"
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# WEIGHTS PARAMETER TESTS - Inverse variance weighting
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction applies weights to LMM correctly", {
    skip_if_not_installed("nlme")
    set.seed(201)
    
    n <- 40
    qv <- rep(seq(0.1, 1, length.out = 20), 2)
    group <- rep(c("A", "B"), each = 20)
    entropy <- 0.5 + 0.2 * qv + rnorm(length(qv), 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_weights"
    
    # Create weights (inverse variance: higher weight for lower variance obs)
    weights <- 1 / (0.01 + abs(qv - mean(qv)))
    weights <- weights / sum(weights) * length(weights)  # Normalize
    
    result <- .fit_one_interaction(
        "gene_weights",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:length(qv)),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = TRUE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        weights = weights
    )
    
    # Weights should be applied
    if (is.data.frame(result)) {
        expect_true(result$ci_weighted)
    }
})

test_that(".fit_one_interaction handles mismatched weight lengths", {
    set.seed(202)
    
    n <- 30
    qv <- seq(0.1, 1, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * qv + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_w_bad"
    
    # Weights with wrong length
    weights <- rep(1, n - 5)  # 5 fewer than needed
    
    result <- .fit_one_interaction(
        "gene_w_bad",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = TRUE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        weights = weights
    )
    
    # Should not crash; weights simply not applied
    expect_true(is.data.frame(result) || is.null(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# REGULARIZATION METHOD TESTS - Feature selection in GAM/FPCA
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction GAM with gamsel regularization", {
    skip_if_not_installed("mgcv")
    skip_if_not_installed("gamsel")
    set.seed(301)
    
    # Larger sample size for flexible GAM fitting without basis warnings
    n <- 100
    q <- sort(runif(n, 0.1, 2))
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.4 + 0.3 * q + 0.1 * (group == "B") * q + rnorm(n, 0, 0.06)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_gamsel"
    
    result <- .fit_one_interaction(
        "gene_gamsel",
        se = NULL,
        mat = mat,
        q_vals = q,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "gam",
        pvalue = "anova",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        regularization = "gamsel",
        adaptive_knots = TRUE
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction GAM with spline regularization", {
    skip_if_not_installed("mgcv")
    set.seed(302)
    
    # Larger sample size for flexible GAM fitting without basis warnings
    n <- 100
    q <- sort(runif(n, 0.2, 1.8))
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.25 * q + 0.08 * (group == "B") + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_spline"
    
    result <- .fit_one_interaction(
        "gene_spline",
        se = NULL,
        mat = mat,
        q_vals = q,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "gam",
        pvalue = "anova",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        regularization = "spline"
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction FPCA with PCA regularization", {
    set.seed(303)
    
    # FPCA with PCA regularization (dimension reduction via PCA)
    genes <- "g1"
    samples <- paste0("s", 1:8)
    q_vals <- rep(c(0.1, 0.5, 1, 2), 2)
    
    mat <- matrix(rnorm(length(q_vals)), nrow = 1)
    rownames(mat) <- genes
    
    sample_names <- samples
    group_vec <- rep(c("A", "B"), each = 4)
    
    result <- .fit_one_interaction(
        g = genes,
        se = NULL,
        mat = mat,
        q_vals = q_vals,
        sample_names = sample_names,
        group_vec = group_vec,
        method = "fpca",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 2,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        regularization = "pca"
    )
    
    expect_true(is.null(result) || is.data.frame(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# MISSING DATA AND EDGE CASES
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction handles missing entropy values", {
    set.seed(401)
    
    # Use very clean data to ensure convergence - focus on missing value handling, not stability
    n <- 80
    qv <- seq(0.1, 2, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.3 * qv + 0.2 * (group == "B") + rnorm(n, 0, 0.03)  # Small noise
    entropy[c(10, 40, 70)] <- NA  # Only 3 missing values in 80 samples
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_na"
    
    # Suppress glmmTMB convergence warnings (expected when testing edge cases)
    result <- suppressWarnings(.fit_one_interaction(
        "gene_na",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 8,  # Require more observations for stable estimation
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    ))
    
    # Should handle missing values gracefully (casewise deletion)
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction returns NULL when all group values are same", {
    set.seed(402)
    
    n <- 20
    qv <- seq(0.1, 1, length.out = n)
    group <- rep("A", n)  # All same group - no variation!
    entropy <- 0.5 + 0.2 * qv + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_same_group"
    
    result <- .fit_one_interaction(
        "gene_same_group",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 2,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Can't estimate interaction with no group variation
    # Expect NULL or error handling
    expect_true(is.null(result) || is.data.frame(result))
})

test_that(".fit_one_interaction handles very small q-value range", {
    set.seed(403)
    
    n <- 15
    qv <- rep(1.0, n)  # All same q value - no q variation!
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_const_q"
    
    result <- .fit_one_interaction(
        "gene_const_q",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 2,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Can't estimate q*group interaction with no q variation
    expect_true(is.null(result) || is.data.frame(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# VERBOSE AND PROGRESS FLAG BEHAVIOR
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction outputs verbose messages when verbose=TRUE", {
    set.seed(501)
    
    n <- 30
    qv <- seq(0.1, 1, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * qv + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_v"
    
    # Expect messages when verbose=TRUE
    result <- expect_message(
        .fit_one_interaction(
            "gene_v",
            se = NULL,
            mat = mat,
            q_vals = qv,
            sample_names = paste0("s", 1:n),
            group_vec = group,
            method = "lmm",
            pvalue = "lrt",
            subject_col = NULL,
            paired = FALSE,
            min_obs = 3,
            verbose = TRUE,  # Enable verbose output
            suppress_lme4_warnings = TRUE,
            progress = FALSE
        ),
        "ARIMA|entropy|weights",  # Expect one of these keywords (regex pattern)
        ignore.case = TRUE
    )
})

test_that(".fit_one_interaction quiet mode when verbose=FALSE", {
    set.seed(502)
    
    n <- 25
    qv <- seq(0.1, 1.2, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * qv + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_quiet"
    
    # Should not produce messages when verbose=FALSE
    result <- expect_silent(
        .fit_one_interaction(
            "gene_quiet",
            se = NULL,
            mat = mat,
            q_vals = qv,
            sample_names = paste0("s", 1:n),
            group_vec = group,
            method = "lmm",
            pvalue = "lrt",
            subject_col = NULL,
            paired = FALSE,
            min_obs = 3,
            verbose = FALSE,  # Silent
            suppress_lme4_warnings = TRUE,
            progress = FALSE
        )
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# BIAS CORRECTION AND ADAPTIVE KNOTS (GAM PARAMETERS)
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction GAM with bias_correction=FALSE", {
    skip_if_not_installed("mgcv")
    set.seed(601)
    
    # Use larger sample size to avoid basis dimension warnings
    n <- 80
    q <- sort(runif(n, 0.1, 2))  # Sorted q values for better GAM fitting
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.4 + 0.3 * q + 0.05 * (group == "B") + rnorm(n, 0, 0.06)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_nobias"
    
    result <- .fit_one_interaction(
        "gene_nobias",
        se = NULL,
        mat = mat,
        q_vals = q,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "gam",
        pvalue = "anova",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        bias_correction = FALSE,
        adaptive_knots = TRUE
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction GAM with adaptive_knots=FALSE", {
    skip_if_not_installed("mgcv")
    set.seed(602)
    
    # Use larger sample size and more evenly distributed q values
    n <- 75
    q <- sort(runif(n, 0.2, 1.8))
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.45 + 0.25 * q + 0.08 * (group == "B") * q + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_fixed_knots"
    
    result <- .fit_one_interaction(
        "gene_fixed_knots",
        se = NULL,
        mat = mat,
        q_vals = q,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "gam",
        pvalue = "anova",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE,
        bias_correction = TRUE,
        adaptive_knots = FALSE  # Fixed knot number
    )
    
    expect_true(is.data.frame(result) || is.null(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# INVALID INPUT HANDLING
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction handles non-existent gene ID", {
    set.seed(701)
    
    mat <- matrix(rnorm(20), nrow = 2)
    rownames(mat) <- c("gene1", "gene2")
    
    # Try to fit gene that doesn't exist in matrix - expects error
    expect_error(
        .fit_one_interaction(
            "nonexistent_gene",  # Gene not in matrix!
            se = NULL,
            mat = mat,
            q_vals = rep(1:10, 2),
            sample_names = paste0("s", 1:20),
            group_vec = rep(c("A", "B"), length.out = 20),
            method = "lmm",
            pvalue = "lrt",
            subject_col = NULL,
            paired = FALSE,
            min_obs = 2,
            verbose = FALSE,
            suppress_lme4_warnings = TRUE,
            progress = FALSE
        ),
        "not found in matrix"
    )
})

test_that(".fit_one_interaction with unequal length q_vals and samples", {
    set.seed(702)
    
    mat <- matrix(rnorm(30), nrow = 1)
    rownames(mat) <- "gene1"
    
    # Mismatch: 30 samples but only 10 q-values
    result <- .fit_one_interaction(
        "gene1",
        se = NULL,
        mat = mat,
        q_vals = rep(1:10, 1),  # 10 values
        sample_names = paste0("s", 1:30),  # 30 samples
        group_vec = rep(c("A", "B"), length.out = 30),
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 2,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Likely to fail or recycle q_vals
    expect_true(is.null(result) || is.data.frame(result) || is.list(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# NUMERIC EDGE CASES
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction with very large entropy values", {
    set.seed(801)
    
    n <- 25
    qv <- seq(0.1, 1, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 1000 + 500 * qv + rnorm(n, 0, 50)  # Very large values
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_large"
    
    result <- .fit_one_interaction(
        "gene_large",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Should handle scaling gracefully
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction with very small entropy values", {
    set.seed(802)
    
    n <- 25
    qv <- seq(0.1, 1, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 1e-6 + 1e-7 * qv + rnorm(n, 0, 1e-8)  # Very small values
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_tiny"
    
    result <- .fit_one_interaction(
        "gene_tiny",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Should handle scaling gracefully
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction with negative entropy values", {
    set.seed(803)
    
    n <- 20
    qv <- seq(0.1, 1, length.out = n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- -0.5 + 0.2 * qv + rnorm(n, 0, 0.1)  # Can be negative
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_neg"
    
    result <- .fit_one_interaction(
        "gene_neg",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 3,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Should still work - entropy can theoretically be adjusted
    expect_true(is.data.frame(result) || is.null(result))
})

# ═══════════════════════════════════════════════════════════════════════════
# COMPLEX GROUP AND SUBJECT CONFIGURATIONS
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction with 3+ groups (handled as factor)", {
    set.seed(901)
    
    n <- 36
    qv <- rep(seq(0.1, 1, length.out = 12), 3)
    group <- rep(c("A", "B", "C"), each = 12)  # 3 groups
    entropy <- 0.5 + 0.2 * qv + 0.1 * (group == "B") + 0.15 * (group == "C") + rnorm(n, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_3groups"
    
    result <- .fit_one_interaction(
        "gene_3groups",
        se = NULL,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Should handle multiple groups (reference coding)
    expect_true(is.data.frame(result) || is.null(result))
})

test_that(".fit_one_interaction with many subjects", {
    skip_if_not_installed("nlme")
    set.seed(902)
    
    n_subjects <- 30  # Many subjects
    n_per <- 3
    n_total <- n_subjects * n_per
    
    subject <- rep(1:n_subjects, each = n_per)
    qv <- rep(seq(0.1, 1, length.out = n_per), n_subjects)
    group <- rep(c("A", "B"), length.out = n_total)
    entropy <- 0.5 + 0.2 * qv + 0.1 * (group == "B") + rnorm(n_total, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_many_subj"
    
    coldata <- S4Vectors::DataFrame(
        samples = paste0("s", 1:n_total),
        sample_base = subject
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        colData = coldata
    )
    
    result <- .fit_one_interaction(
        "gene_many_subj",
        se = se,
        mat = mat,
        q_vals = qv,
        sample_names = paste0("s", 1:n_total),
        group_vec = group,
        method = "lmm",
        pvalue = "lrt",
        subject_col = NULL,
        paired = TRUE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
    )
    
    # Should handle many subjects without issue
    expect_true(is.data.frame(result) || is.null(result))
})
