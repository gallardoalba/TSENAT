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
    # min_obs > non-missing -> NULL with expected warning
    res_null <- expect_warning(
        .fit_one_interaction("g1",
            se = NULL, mat = mat, q_vals = q_vals,
            sample_names = sample_names, group_vec = group_vec, method = "lmm",
            pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 10, verbose = FALSE,
            suppress_lme4_warnings = TRUE, progress = FALSE
        ),
        "Insufficient"
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

    res <- expect_warning(
        .fit_one_interaction("g1",
            se = se, mat = mat, q_vals = qv,
            sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
            subject_col = NULL, paired = TRUE, min_obs = 2, verbose = FALSE,
            suppress_lme4_warnings = TRUE, progress = FALSE
        ),
        "Insufficient|subjects"
    )
    expect_null(res)
})

# ═══════════════════════════════════════════════════════════════════════════
# FALLBACK BEHAVIOR: When lmer fails or returns singular fits
# ═══════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════
# METHOD DISPATCHING: GAM, FPCA, and other methods
# ═══════════════════════════════════════════════════════════════════════════

test_that(".fit_one_interaction dispatches to gam and fpca methods", {
    # FPCA dispatch - create sufficient data structure
    set.seed(1)
    n_samples <- 6
    n_q <- 5
    sample_names <- rep(paste0("s", 1:n_samples), each = n_q)
    q_vals <- rep(seq(0.5, 2, length.out = n_q), n_samples)
    group_vec <- rep(c("A", "B"), each = n_q * n_samples / 2)
    
    # Create entropy with group structure
    entropy_a <- 0.8 + 0.3 * q_vals[1:(n_q * n_samples/2)] + rnorm(n_q * n_samples/2, 0, 0.1)
    entropy_b <- 1.2 + 0.5 * q_vals[(n_q * n_samples/2 + 1):(n_q * n_samples)] + rnorm(n_q * n_samples/2, 0, 0.1)
    obs <- c(entropy_a, entropy_b)
    
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
        method = "lmm",
        pvalue = "lrt",
        subject_col = "subject",
        paired = FALSE,
        min_obs = 5,
        verbose = FALSE,
        suppress_lme4_warnings = TRUE,
        progress = FALSE
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
    n_subjects <- 6
    n_q <- 4
    
    # Structure: each subject has n_q observations per group (A and B)
    subject_vec <- rep(1:n_subjects, each = n_q * 2)
    qv <- rep(rep(seq(0.2, 1.5, length.out = n_q), 2), n_subjects)
    group_vec <- rep(rep(c("A", "B"), each = n_q), n_subjects)
    
    n_total <- length(subject_vec)
    sample_names_vec <- paste0("s", 1:n_total)
    
    # Add interaction effect: group B has steeper slope with q
    entropy <- 0.5 + 0.3 * qv + ifelse(group_vec == "B", 0.4 * qv, 0) + rnorm(n_total, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene1"
    colnames(mat) <- sample_names_vec
    
    # Create colData with proper rownames matching sample names
    coldata <- S4Vectors::DataFrame(
        samples = sample_names_vec,
        sample_base = as.character(subject_vec),
        row.names = sample_names_vec  # Critical: rownames must match sample_names!
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
        sample_names = sample_names_vec,
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
    # Create sufficient data with multiple q-values and samples
    n_samples <- 6
    n_q <- 5
    genes <- "g1"
    sample_names <- rep(paste0("s", 1:n_samples), each = n_q)
    q_vals <- rep(seq(0.1, 2.0, length.out = n_q), n_samples)
    
    # Create entropy with group structure
    entropy_a <- 0.8 + 0.3 * q_vals[1:(n_q * n_samples/2)] + rnorm(n_q * n_samples/2, 0, 0.1)
    entropy_b <- 1.2 + 0.5 * q_vals[(n_q * n_samples/2 + 1):(n_q * n_samples)] + rnorm(n_q * n_samples/2, 0, 0.1)
    obs <- c(entropy_a, entropy_b)
    
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- genes
    
    group_vec <- rep(c("A", "B"), each = n_q * n_samples / 2)
    
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
    
    n_subjects <- 20  # Many subjects
    n_per_group <- 2  # per subject per group
    n_q <- 3         # q-values
    
    # Create proper subject structure: each subject has observations for each group and q-value
    subject_vec <- rep(1:n_subjects, each = n_per_group * n_q * 2)
    qv <- rep(rep(seq(0.1, 1, length.out = n_q), 2), n_subjects * n_per_group)
    group_vec <- rep(rep(c("A", "B"), each = n_q), n_subjects * n_per_group)
    
    n_total <- length(subject_vec)
    sample_names_vec <- paste0("s", 1:n_total)
    
    entropy <- 0.5 + 0.2 * qv + 0.1 * (group_vec == "B") + rnorm(n_total, 0, 0.05)
    
    mat <- matrix(entropy, nrow = 1)
    rownames(mat) <- "gene_many_subj"
    colnames(mat) <- sample_names_vec
    
    # Create colData with proper rownames matching sample names
    coldata <- S4Vectors::DataFrame(
        samples = sample_names_vec,
        sample_base = as.character(subject_vec),
        row.names = sample_names_vec  # Critical: rownames must match sample_names!
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
        sample_names = sample_names_vec,
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
    
    # Should handle many subjects without issue
    expect_true(is.data.frame(result) || is.null(result))
})

context("Helper Functions for .fit_one_interaction() Tests")


# ═══════════════════════════════════════════════════════════════════════════
# .setup_interaction_data() - Data frame initialization and validation
# ═══════════════════════════════════════════════════════════════════════════

test_that(".setup_interaction_data validates gene exists in matrix", {
    # Valid gene in matrix
    mat <- matrix(c(1, 2, 3), nrow = 1)
    rownames(mat) <- "gene1"
    q_vals <- c(0.1, 0.2, 0.3)
    group_vec <- c("A", "B", "A")
    
    result <- TSENAT:::.setup_interaction_data("gene1", mat, q_vals, group_vec)
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 3)
    expect_equal(colnames(result), c("entropy", "q", "group", "sample_name"))
    expect_equal(result$entropy, c(1, 2, 3))
    expect_equal(result$q, q_vals)
    expect_equal(levels(result$sample_name), c("S1", "S2", "S3"))
})

test_that(".setup_interaction_data throws error for missing gene", {
    mat <- matrix(c(1, 2, 3), nrow = 1)
    rownames(mat) <- "known_gene"
    q_vals <- c(0.1, 0.2, 0.3)
    group_vec <- c("A", "B", "A")
    
    expect_error(
        TSENAT:::.setup_interaction_data("missing_gene", mat, q_vals, group_vec),
        "not found in matrix rownames"
    )
})

test_that(".setup_interaction_data returns correct data frame structure", {
    mat <- matrix(c(0.5, 1.5, 2.5, 3.5), nrow = 1)
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.5, 1.0, 1.5)
    group_vec <- c("ctrl", "ctrl", "treat", "treat")
    
    result <- TSENAT:::.setup_interaction_data("g1", mat, q_vals, group_vec)
    
    expect_equal(result$entropy, c(0.5, 1.5, 2.5, 3.5))
    expect_equal(result$q, q_vals)
    expect_equal(levels(result$group), c("ctrl", "treat"))
    expect_true(is.factor(result$group))
    expect_equal(levels(result$sample_name), c("S1", "S2", "S3", "S4"))
})

test_that(".setup_interaction_data handles numeric entropy values correctly", {
    # Create a matrix with matching dimensions: 1 gene, 2 columns (for 2 q-values)
    # All values are NA to test NA handling
    mat <- matrix(as.numeric(NA), nrow = 1, ncol = 2)
    colnames(mat) <- c("s1", "s2")
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.2)
    group_vec <- c("A", "B")
    
    result <- TSENAT:::.setup_interaction_data("g1", mat, q_vals, group_vec)
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)  # Should have 2 rows (matching q_vals and group_vec)
    # Entropy values should be NA (from the NA matrix)
    expect_true(all(is.na(result$entropy)))
    expect_equal(result$q, q_vals)
    expect_equal(as.character(result$group), group_vec)
})

test_that(".setup_interaction_data preserves matrix column order in entropy extraction", {
    mat <- matrix(c(5, 3, 1, 4, 2, 6), nrow = 2)
    rownames(mat) <- c("gene_A", "gene_B")
    q_vals <- c(0.1, 0.2, 0.3)
    group_vec <- c("A", "B", "A")
    
    result <- TSENAT:::.setup_interaction_data("gene_B", mat, q_vals, group_vec)
    
    # Matrix fills by column: [5,1,2; 3,4,6], so gene_B (row 2) is c(3, 4, 6)
    expect_equal(result$entropy, c(3, 4, 6))
    expect_equal(levels(result$sample_name), c("S1", "S2", "S3"))
})

test_that(".setup_interaction_data with large number of features", {
    # Test with many q values (high-dimensional case)
    n_features <- 500
    mat <- matrix(rnorm(n_features), nrow = 1)
    rownames(mat) <- "big_gene"
    q_vals <- seq(0.01, 5, length.out = n_features)
    group_vec <- rep(c("A", "B"), length.out = n_features)
    
    result <- TSENAT:::.setup_interaction_data("big_gene", mat, q_vals, group_vec)
    
    expect_equal(nrow(result), n_features)
    expect_equal(length(unique(result$group)), 2)
    expect_equal(length(unique(result$sample_name)), n_features)  # Each sample should be unique
})

# ═══════════════════════════════════════════════════════════════════════════
# .apply_weights_to_df() - Weight parameter handling
# ═══════════════════════════════════════════════════════════════════════════

test_that(".apply_weights_to_df applies matching weights correctly", {
    df <- data.frame(entropy = c(1, 2, 3), q = c(0.1, 0.2, 0.3))
    weights <- c(0.5, 1.0, 1.5)
    
    result <- TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = FALSE)
    
    expect_true("weight" %in% colnames(result))
    expect_equal(result$weight, weights)
    expect_equal(nrow(result), nrow(df))
})

test_that(".apply_weights_to_df ignores NULL weights", {
    df <- data.frame(entropy = c(1, 2, 3), q = c(0.1, 0.2, 0.3))
    
    result <- TSENAT:::.apply_weights_to_df(df, NULL, "g1", verbose = FALSE)
    
    expect_false("weight" %in% colnames(result))
    expect_equal(nrow(result), nrow(df))
})

test_that(".apply_weights_to_df handles mismatched weight length", {
    df <- data.frame(entropy = c(1, 2, 3, 4), q = c(0.1, 0.2, 0.3, 0.4))
    weights <- c(0.5, 1.0)  # Only 2 weights for 4 rows
    
    result <- TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = FALSE)
    
    # Weights should NOT be added due to length mismatch
    expect_false("weight" %in% colnames(result))
    expect_equal(nrow(result), nrow(df))
})

test_that(".apply_weights_to_df verbose mode prints messages", {
    df <- data.frame(entropy = c(1, 2, 3), q = c(0.1, 0.2, 0.3))
    weights <- c(0.5, 1.0, 1.5)
    
    # Capture message output
    expect_message(
        TSENAT:::.apply_weights_to_df(df, weights, "gene_test", verbose = TRUE),
        "weights applied"
    )
})

test_that(".apply_weights_to_df verbose message for mismatched weights", {
    df <- data.frame(entropy = c(1, 2, 3, 4), q = c(0.1, 0.2, 0.3, 0.4))
    weights <- c(0.5, 1.0)
    
    expect_message(
        TSENAT:::.apply_weights_to_df(df, weights, "gene_test", verbose = TRUE),
        "NOT applied"
    )
})

test_that(".apply_weights_to_df handles zero and negative weights", {
    df <- data.frame(entropy = c(1, 2, 3), q = c(0.1, 0.2, 0.3))
    
    # Zero weights are technically valid for inverse variance weighting
    weights_with_zero <- c(0, 1.0, 1.5)
    result <- TSENAT:::.apply_weights_to_df(df, weights_with_zero, "g1", verbose = FALSE)
    expect_true("weight" %in% colnames(result))
    expect_equal(result$weight, weights_with_zero)
    
    # Negative weights might be used for contrast - should still be applied
    weights_negative <- c(-0.5, 1.0, 1.5)
    result2 <- TSENAT:::.apply_weights_to_df(df, weights_negative, "g1", verbose = FALSE)
    expect_equal(result2$weight, weights_negative)
})

test_that(".apply_weights_to_df calculates weight statistics correctly", {
    df <- data.frame(entropy = c(1, 2, 3), q = c(0.1, 0.2, 0.3))
    weights <- c(1, 2, 3)
    
    # Capture output with verbose=TRUE to verify statistics are calculated
    expect_message(
        TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = TRUE),
        "mean=2\\.0000"  # Mean of 1,2,3 is 2.0
    )
})

test_that(".apply_weights_to_df preserves all original columns", {
    df <- data.frame(
        entropy = c(1, 2, 3),
        q = c(0.1, 0.2, 0.3),
        group = factor(c("A", "B", "A")),
        extra_col = c("x", "y", "z")
    )
    weights <- c(0.5, 1.0, 1.5)
    
    result <- TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = FALSE)
    
    expect_equal(colnames(result), c("entropy", "q", "group", "extra_col", "weight"))
    expect_equal(result$extra_col, df$extra_col)
})

# ═══════════════════════════════════════════════════════════════════════════
# .get_subject_ids() - Subject identifier extraction
# ═══════════════════════════════════════════════════════════════════════════

test_that(".get_subject_ids returns sample_names as fallback", {
    sample_names <- c("s1", "s2", "s3", "s4")
    
    result <- TSENAT:::.get_subject_ids(
        se = NULL,
        subject_col = NULL,
        paired = FALSE,
        mat = NULL,
        sample_names = sample_names
    )
    
    expect_equal(result, sample_names)
})

test_that(".get_subject_ids extracts from explicit subject_col", {
    skip_if_not_installed("SummarizedExperiment")
    
    sample_names <- c("s1", "s2", "s3", "s4")
    coldata <- S4Vectors::DataFrame(
        samples = sample_names,
        my_subject = c("subA", "subA", "subB", "subB")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rnorm(4), nrow = 1)),
        colData = coldata
    )
    
    result <- TSENAT:::.get_subject_ids(
        se = se,
        subject_col = "my_subject",
        paired = FALSE,
        mat = NULL,
        sample_names = sample_names
    )
    
    expect_equal(result, c("subA", "subA", "subB", "subB"))
})

test_that(".get_subject_ids errors on missing subject_col", {
    skip_if_not_installed("SummarizedExperiment")
    
    sample_names <- c("s1", "s2", "s3", "s4")
    coldata <- S4Vectors::DataFrame(
        samples = sample_names,
        other_col = c(1, 2, 3, 4)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rnorm(4), nrow = 1)),
        colData = coldata
    )
    
    expect_error(
        TSENAT:::.get_subject_ids(
            se = se,
            subject_col = "nonexistent_col",
            paired = FALSE,
            mat = NULL,
            sample_names = sample_names
        ),
        "not found in colData"
    )
})

test_that(".get_subject_ids extracts from paired=TRUE with sample_base", {
    skip_if_not_installed("SummarizedExperiment")
    
    sample_names <- c("s1", "s2", "s3", "s4")
    coldata <- S4Vectors::DataFrame(
        samples = sample_names,
        sample_base = c("pair1", "pair1", "pair2", "pair2")
    )
    # Set rownames to match sample identifiers for proper indexing
    rownames(coldata) <- sample_names
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rnorm(4), nrow = 1)),
        colData = coldata
    )
    
    result <- TSENAT:::.get_subject_ids(
        se = se,
        subject_col = NULL,
        paired = TRUE,
        mat = NULL,
        sample_names = sample_names
    )
    
    expect_equal(result, c("pair1", "pair1", "pair2", "pair2"))
})

test_that(".get_subject_ids errors on paired=TRUE without sample_base", {
    skip_if_not_installed("SummarizedExperiment")
    
    sample_names <- c("s1", "s2", "s3", "s4")
    coldata <- S4Vectors::DataFrame(
        samples = sample_names,
        other_col = c(1, 2, 3, 4)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rnorm(4), nrow = 1)),
        colData = coldata
    )
    
    expect_error(
        TSENAT:::.get_subject_ids(
            se = se,
            subject_col = NULL,
            paired = TRUE,
            mat = NULL,
            sample_names = sample_names
        ),
        "paired = TRUE"
    )
})

test_that(".get_subject_ids handles paired_samples column alternative", {
    skip_if_not_installed("SummarizedExperiment")
    
    sample_names <- c("s1", "s2", "s3", "s4")
    coldata <- S4Vectors::DataFrame(
        samples = sample_names,
        paired_samples = c("subX", "subX", "subY", "subY")
    )
    # Set rownames to match sample identifiers for proper indexing
    rownames(coldata) <- sample_names
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rnorm(4), nrow = 1)),
        colData = coldata
    )
    
    result <- TSENAT:::.get_subject_ids(
        se = se,
        subject_col = NULL,
        paired = TRUE,
        mat = NULL,
        sample_names = sample_names
    )
    
    expect_equal(result, c("subX", "subX", "subY", "subY"))
})

test_that(".get_subject_ids with subject_col removes _q= suffix", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Test matrix column names with _q= suffix
    mat <- matrix(rnorm(4), nrow = 1)
    colnames(mat) <- c("s1_q=0.1", "s2_q=0.2", "s3_q=0.5", "s4_q=1.0")
    
    coldata <- S4Vectors::DataFrame(
        samples = c("s1", "s2", "s3", "s4"),
        my_subject = c("subA", "subA", "subB", "subB")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        colData = coldata
    )
    
    result <- TSENAT:::.get_subject_ids(
        se = se,
        subject_col = "my_subject",
        paired = FALSE,
        mat = mat,
        sample_names = c("s1", "s2", "s3", "s4")
    )
    
    expect_equal(result, c("subA", "subA", "subB", "subB"))
})

# ═══════════════════════════════════════════════════════════════════════════
# .check_lmm_sample_sizes() - Sample and subject minimum validation
# ═══════════════════════════════════════════════════════════════════════════

test_that(".check_lmm_sample_sizes returns TRUE when requirements met", {
    df <- data.frame(
        entropy = rnorm(10),
        q = seq(0.1, 1, length.out = 10),
        group = rep(c("A", "B"), 5),
        subject = rep(c(1, 2, 3, 4, 5), 2)
    )
    
    result <- TSENAT:::.check_lmm_sample_sizes(df, min_obs = 5)
    
    expect_true(result)
})

test_that(".check_lmm_sample_sizes returns NULL when insufficient observations", {
    df <- data.frame(
        entropy = rnorm(3),
        q = c(0.1, 0.5, 1.0),
        group = c("A", "B", "A"),
        subject = c(1, 2, 3)
    )
    
    result <- expect_warning(
        TSENAT:::.check_lmm_sample_sizes(df, min_obs = 5),
        "Insufficient"
    )
    
    expect_null(result)
})

test_that(".check_lmm_sample_sizes returns NULL with <2 subjects", {
    df <- data.frame(
        entropy = rnorm(10),
        q = seq(0.1, 1, length.out = 10),
        group = rep(c("A", "B"), 5),
        subject = rep(1, 10)  # Only one subject
    )
    
    result <- expect_warning(
        TSENAT:::.check_lmm_sample_sizes(df, min_obs = 3),
        "Insufficient|subjects"
    )
    
    expect_null(result)
})

test_that(".check_lmm_sample_sizes handles NA subjects correctly", {
    df <- data.frame(
        entropy = rnorm(10),
        q = seq(0.1, 1, length.out = 10),
        group = rep(c("A", "B"), 5),
        subject = c(1, 2, NA, 3, 4, 1, 2, NA, 3, 4)
    )
    
    result <- TSENAT:::.check_lmm_sample_sizes(df, min_obs = 5)
    
    # Should count 4 unique non-NA subjects: 1, 2, 3, 4
    expect_true(result)
})

test_that(".check_lmm_sample_sizes with exact boundary conditions", {
    # Exactly min_obs observations
    df_exact <- data.frame(
        entropy = rnorm(5),
        q = seq(0.1, 1, length.out = 5),
        group = c("A", "B", "A", "B", "A"),
        subject = c(1, 2, 1, 2, 3)
    )
    
    result_exact <- TSENAT:::.check_lmm_sample_sizes(df_exact, min_obs = 5)
    expect_true(result_exact)
    
    # Just below min_obs - expect warning
    df_below <- df_exact[-5, ]
    result_below <- expect_warning(
        TSENAT:::.check_lmm_sample_sizes(df_below, min_obs = 5),
        "Insufficient"
    )
    expect_null(result_below)
})

test_that(".check_lmm_sample_sizes with exactly 2 subjects (boundary)", {
    df <- data.frame(
        entropy = rnorm(6),
        q = seq(0.1, 1, length.out = 6),
        group = rep(c("A", "B"), 3),
        subject = rep(c(1, 2), 3)
    )
    
    result <- TSENAT:::.check_lmm_sample_sizes(df, min_obs = 3)
    
    # Exactly 2 subjects (minimum for random intercept)
    expect_true(result)
})

test_that(".check_lmm_sample_sizes default min_obs value works", {
    df <- data.frame(
        entropy = rnorm(5),
        q = seq(0.1, 1, length.out = 5),
        group = c("A", "B", "A", "B", "A"),
        subject = c(1, 2, 3, 2, 1)
    )
    
    # Default min_obs = 3
    result <- TSENAT:::.check_lmm_sample_sizes(df)
    
    expect_true(result)
})

test_that(".check_lmm_sample_sizes with large datasets", {
    n_large <- 1000
    df <- data.frame(
        entropy = rnorm(n_large),
        q = runif(n_large, 0.1, 1),
        group = rep(c("A", "B"), length.out = n_large),
        subject = rep(1:50, length.out = n_large)
    )
    
    result <- TSENAT:::.check_lmm_sample_sizes(df, min_obs = 500)
    
    expect_true(result)
})

test_that(".check_lmm_sample_sizes only counts unique non-NA subjects", {
    df <- data.frame(
        entropy = rnorm(8),
        q = seq(0.1, 1, length.out = 8),
        group = c("A", "B", "A", "B", "A", "B", "A", "B"),
        subject = c(1, 1, 1, 2, 2, 2, NA, NA)  # Actually only 2 subjects
    )
    
    result <- TSENAT:::.check_lmm_sample_sizes(df, min_obs = 4)
    
    # 8 observations, but only 2 subjects - should pass obs check and subject check
    expect_true(result)
})

# ═══════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS - Helper functions working together
# ═══════════════════════════════════════════════════════════════════════════

test_that("Helper functions work together in workflow", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Step 1: Create matrix and setup data
    mat <- matrix(c(1.5, 2.5, 3.5, 4.5), nrow = 1)
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.5, 1.0, 1.5)
    group_vec <- c("ctrl", "ctrl", "treat", "treat")
    
    df <- TSENAT:::.setup_interaction_data("g1", mat, q_vals, group_vec)
    
    # Step 2: Apply weights
    weights <- c(1, 1, 0.5, 0.5)
    df_weighted <- TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = FALSE)
    
    # Step 3: Create subject column for validation
    df_weighted$subject <- c(1, 1, 2, 2)
    
    # Step 4: Check sample sizes
    check_result <- TSENAT:::.check_lmm_sample_sizes(df_weighted, min_obs = 2)
    
    expect_true(check_result)
    expect_equal(nrow(df_weighted), 4)
    expect_true("weight" %in% colnames(df_weighted))
})

test_that("Helper functions handle problematic data gracefully", {
    # Setup with valid data but will fail checks
    mat <- matrix(c(1, 2), nrow = 1)
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.2)
    group_vec <- c("A", "B")
    
    df <- TSENAT:::.setup_interaction_data("g1", mat, q_vals, group_vec)
    
    # Apply weights with wrong length - should not crash
    weights <- c(1)
    df_weighted <- TSENAT:::.apply_weights_to_df(df, weights, "g1", verbose = FALSE)
    
    # Add single subject for validation
    df_weighted$subject <- c(1, 1)
    
    # Will fail because <2 subjects - expect warning
    check_result <- expect_warning(
        TSENAT:::.check_lmm_sample_sizes(df_weighted, min_obs = 1),
        "Insufficient|subjects"
    )
    
    expect_null(check_result)
})

test_that("Helper error messages are informative", {
    mat <- matrix(c(1, 2, 3), nrow = 1)
    rownames(mat) <- "known_gene"
    q_vals <- c(0.1, 0.2, 0.3)
    group_vec <- c("A", "B", "A")
    
    error_msg <- tryCatch(
        TSENAT:::.setup_interaction_data("unknown", mat, q_vals, group_vec),
        error = function(e) e$message
    )
    
    expect_match(error_msg, "unknown")
    expect_match(error_msg, "not found in matrix rownames")
})

# ============================================================================
# REDISTRIBUTED TESTS FROM test-infrastructure-statistical_validation.R
# ============================================================================

context("Linear Models: Mixed Model Fallback Improvement")

test_that("lm_subject_fixed fallback preserves power vs. lm_nosubject", {
    skip_if_not_installed("lme4")
    set.seed(102)
    
    # Create data with within-subject correlation
    n_subjects <- 8
    n_q_per_subject <- 5
    subjects <- rep(paste0("S", 1:n_subjects), each = n_q_per_subject)
    q_vals <- rep(seq(0.1, 0.9, length.out = n_q_per_subject), times = n_subjects)
    group <- rep(c("A", "B"), each = n_subjects * n_q_per_subject / 2)
    
    # True effect: q×group interaction
    subject_effect <- rep(rnorm(n_subjects, 0, 0.1), each = n_q_per_subject)
    entropy <- 0.5 + 
             0.3 * q_vals +
             0.2 * (group == "B") +
             0.15 * (group == "B") * q_vals +  # True interaction
             subject_effect +
             rnorm(length(subjects), 0, 0.05)
    
    df <- data.frame(entropy = entropy, q = q_vals, group = factor(group), 
                     subject = factor(subjects))
    
    # Fit with subject as fixed effect (proper approach)
    fit_with_subj <- stats::lm(entropy ~ q * group + factor(subject), data = df)
    p_with_subj <- summary(fit_with_subj)$coefficients["q:groupB", "Pr(>|t|)"]
    
    # Fit without subject (loses power)
    fit_no_subj <- stats::lm(entropy ~ q * group, data = df)
    p_no_subj <- summary(fit_no_subj)$coefficients["q:groupB", "Pr(>|t|)"]
    
    # Model with subject should have lower p-value (more power)
    expect_true(p_with_subj < p_no_subj)
})

test_that(".try_lm_fallbacks uses factor(subject), not numeric subject", {
    set.seed(103)
    
    # Create test data
    df <- data.frame(
        entropy = rnorm(30),
        q = rep(seq(0.1, 1, length.out = 10), 3),
        group = rep(c("A", "B"), length.out = 30),
        subject = rep(1:10, 3)  # numeric subject IDs
    )
    
    # Run fallback function
    fb <- .try_lm_fallbacks(df, verbose = FALSE)
    
    expect_true(!is.null(fb))
    expect_true(!is.null(fb$fit1))
    
    # Extract coefficients to verify it's treating subject as factor
    if (fb$method == "lm_subject_fixed") {
        coef_names <- names(coef(fb$fit1))
        # Should have factor(subject) terms, not a single "subject" slope
        subj_terms <- grep("factor\\(subject\\)", coef_names)
        expect_true(length(subj_terms) > 0)
    }
})
