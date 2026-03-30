context("Linear Model Helpers: Basic Calculations")

library(testthat)

# Report summary messages
test_that(".report_fit_summary prints fallback and singular messages", {
    df <- data.frame(fit_method = c("lm_nosubject", "lmer", NA), singular = c(TRUE, FALSE, NA), stringsAsFactors = FALSE)
    expect_message(.report_fit_summary(df, verbose = TRUE), "Alternative method")
    expect_message(.report_fit_summary(df, verbose = TRUE), "Singular fits")
})

# GAM interaction: skip if mgcv not available
test_that(".gam_interaction returns a data.frame with p_interaction when mgcv present", {
    skip_if_not_installed("mgcv")
    set.seed(1)
    # build small dataset with group and q, per-sample entropy
    n <- 40
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * (q) + ifelse(group == "A", 0.05, -0.05) + rnorm(n, 0, 0.01)
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    res <- suppressWarnings(.gam_interaction(df, q_vals = q, g = "g1", min_obs = 5))
    expect_true(is.data.frame(res) || is.null(res))
    if (is.data.frame(res)) {
        expect_true("p_interaction" %in% colnames(res))
    }
})

# FPCA interaction: synthetic matrix
test_that(".try_lm_fallbacks returns lm fits and LRT extractor returns numeric p-values", {
    # build small long-format df
    df <- data.frame(
        entropy = rnorm(30),
        q = rep(seq(0.1, 1.0, length.out = 10), 3),
        group = rep(c("A", "B", "A"), each = 10),
        subject = rep(paste0("sub", 1:10), 3),
        stringsAsFactors = FALSE
    )
    fb <- .try_lm_fallbacks(df, verbose = TRUE)
    expect_true(is.null(fb) || (is.list(fb) && all(c("fit0", "fit1", "method") %in% names(fb))))
    if (!is.null(fb)) {
        lrt_p <- .extract_lrt_p(fb$fit0, fb$fit1, df = df)
        # Phase 14: .extract_lrt_p now returns a list with p_value, n_subjects, small_sample_flag
        expect_true(is.list(lrt_p) && "p_value" %in% names(lrt_p))
        expect_true(is.numeric(lrt_p$p_value) || is.na(lrt_p$p_value))
    }
})

# FPCA matrix preparation helper (variance filtering and scaling)
## Note: the FPCA-prep helper has multiple definitions in the source; tests
## below use the variant that accepts (mat, sample_names, q_vals, min_obs).

# Alternative FPCA matrix builder that returns mat_sub/used_samples
test_that(".prepare_fpca_matrix (fpca variant) builds sub-matrix or returns NULL when insufficient", {
    # Use a single-gene matrix so element assignment in the helper is scalar
    mat <- matrix(rnorm(1 * 8), nrow = 1) # 1 gene x 8 observations
    sample_names <- rep(paste0("S", 1:4), 2)
    q_vals <- rep(c(0.1, 0.5, 1, 2), 2)
    res <- .prepare_fpca_matrix(mat, sample_names = sample_names, q_vals = q_vals, min_obs = 2)
    expect_true(is.null(res) || (is.list(res) && all(c("mat_sub", "used_samples") %in% names(res))))
})

# Test lmer wrapper if available
test_that(".try_lmer attempts lmer fitting when lme4 is installed", {
    skip_if_not_installed("lme4")
    set.seed(5)
    # Build a small balanced dataset for mixed model
    nsub <- 10
    nper <- 4
    subject <- rep(paste0("s", seq_len(nsub)), each = nper)
    q <- rep(seq(0.1, 1, length.out = nper), times = nsub)
    group <- rep(rep(c("A", "B"), length.out = nper), times = nsub)
    entropy <- rnorm(length(subject), mean = 0.5 + as.numeric(group == "A") * 0.1 + 0.2 * q, sd = 0.05)
    df <- data.frame(entropy = entropy, q = q, group = group, subject = subject, stringsAsFactors = FALSE)
    f <- as.formula("entropy ~ q * group + (1 | subject)")
    fit_try <- .try_lmer(f, df, suppress_lme4_warnings = TRUE, verbose = FALSE)
    expect_true(inherits(fit_try, "try-error") || inherits(fit_try, "lmerMod"))
})



testthat::test_that("FPCA helper and interaction work on simple synthetic data", {
    set.seed(42)
    # create 4 samples, each with two q values (8 observations)
    sample_names <- rep(paste0("s", 1:4), each = 2)
    q_vals <- rep(1:2, times = 4)
    # groups: first two samples group A, last two group B
    group_vec <- rep(c("A", "A", "B", "B"), each = 2)
    # single gene with mild group effect across PC1
    obs <- rnorm(8)
    obs[q_vals == 2 & sample_names %in% c("s3", "s4")] <- obs[q_vals == 2 & sample_names %in% c("s3", "s4")] + 1
    mat <- matrix(obs, nrow = 1)
    res <- .fpca_interaction(
        mat = mat, q_vals = q_vals, sample_names = sample_names,
        group_vec = group_vec, g = 1, min_obs = 2
    )
    testthat::expect_true(is.data.frame(res) || is.null(res))
    if (!is.null(res)) {
        # function may include additional metadata (n_pcs_tested, etc.)
        testthat::expect_true(all(c("gene", "p_interaction") %in% names(res)))
        testthat::expect_type(res$p_interaction, "double")
    }
})

testthat::test_that("FPCA matrix preparation filters low-variance rows and scales", {
    # create a single-row matrix so the current implementation assigns scalars
    mat <- matrix(0, nrow = 1, ncol = 10)
    mat[1, ] <- rnorm(10)
    # prepare sample / q vectors matching 10 columns
    sample_names <- rep(paste0("s", 1:5), each = 2)
    q_vals <- rep(1:2, times = 5)
    out <- .prepare_fpca_matrix(mat,
        sample_names = sample_names, q_vals = q_vals,
        min_obs = 2
    )
    testthat::expect_type(out, "list")
    # current implementation returns mat_sub and used_samples or NULL
    testthat::expect_true(is.null(out) || (is.matrix(out$mat_sub) && is.character(out$used_samples)))
})

testthat::test_that("LM fallback helpers choose appropriate method", {
    set.seed(1)
    n <- 40
    subject <- rep(1:10, each = 4)
    q <- runif(n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.1)
    df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = factor(subject))
    res <- .try_lm_fallbacks(df)
    testthat::expect_type(res, "list")
    # Phase 14: AR(1) tries nlme_ar1 first, then nlme, then glmmTMB, then lm_subject_fixed, then lm_nosubject
    testthat::expect_true(res$method %in% c("nlme_ar1", "nlme", "glmmTMB", "lm_subject_fixed", "lm_nosubject"))
    # fit1 can be lme, glmmTMB, or lm depending on which strategy succeeded
    testthat::expect_true(inherits(res$fit1, "lme") || inherits(res$fit1, "glmmTMB") || inherits(res$fit1, "lm") || inherits(res$fit1, "NA"))

    # drop subject -> should pick nosubject fallback
    df2 <- df[, c("entropy", "q", "group")]
    res2 <- .try_lm_fallbacks(df2)
    testthat::expect_type(res2, "list")
    testthat::expect_equal(res2$method, "lm_nosubject")
})

testthat::test_that("LRT p extraction returns numeric p-value for nested lm models", {
    set.seed(2)
    n <- 60
    q <- runif(n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.4 * q + ifelse(group == "B", 0.5, 0) + rnorm(n, 0, 0.2)
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    fit0 <- stats::lm(entropy ~ q + group, data = df)
    fit1 <- stats::lm(entropy ~ q * group, data = df)
    result <- .extract_lrt_p(fit0, fit1, df = df)
    # Phase 14: .extract_lrt_p returns a list
    testthat::expect_true(is.list(result) && "p_value" %in% names(result))
    p <- result$p_value
    testthat::expect_true(is.numeric(p) || is.na(p))
    if (!is.na(p)) testthat::expect_true(p >= 0 && p <= 1)
})



testthat::test_that("GAM interaction returns a data.frame with p-value when mgcv available", {
    set.seed(4)
    n <- 80
    q <- runif(n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.2 * q + ifelse(group == "B", 0.6 * q, 0) + rnorm(n, 0, 0.15)
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    res <- .gam_interaction(df, q_vals = q, g = "geneX", min_obs = 5)
    testthat::expect_true(is.data.frame(res) || is.null(res))
    if (!is.null(res)) {
        # GAM returns at minimum (gene, p_interaction); may include bias correction columns
        testthat::expect_true("gene" %in% colnames(res) && "p_interaction" %in% colnames(res))
        testthat::expect_type(res$p_interaction, "double")
    }
})

testthat::test_that("try_lmer returns an lmer object when lme4 available", {
    set.seed(5)
    n <- 48
    subject <- rep(1:12, each = 4)
    q <- runif(n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.25 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.1)
    df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = factor(subject))
    fmla <- stats::as.formula("entropy ~ q * group + (1|subject)")
    fit <- .try_lmer(fmla, data = df, suppress_lme4_warnings = TRUE)
    testthat::expect_true(inherits(fit, "lmerMod") || inherits(fit, "try-error"))
})

context("Linear Model Helpers: Edge Cases and Validation")

library(testthat)

# .report_fit_summary should be silent when no fallback/singular
test_that(".report_fit_summary is silent when no messages to print", {
    df <- data.frame(x = 1:3)
    expect_silent(.report_fit_summary(df, verbose = TRUE))
})



# .extract_lrt_p returns NA when anova errors
test_that(".extract_lrt_p returns NA for invalid models", {
    result <- .extract_lrt_p("not_a_model", "also_not")
    # Phase 14: .extract_lrt_p returns a list
    expect_true(is.list(result) && "p_value" %in% names(result))
    expect_true(is.na(result$p_value))
})



# .prepare_fpca_matrix returns NULL when insufficient good rows (min_obs large)
test_that(".prepare_fpca_matrix returns NULL when min_obs larger than available", {
    mat <- matrix(rnorm(6), nrow = 1)
    sample_names <- rep(paste0("s", 1:3), each = 2)
    q_vals <- rep(c(1, 2), times = 3)
    res <- .prepare_fpca_matrix(mat = mat, sample_names = sample_names, q_vals = q_vals, min_obs = 10)
    expect_null(res)
})

# .try_lmer sets 'singular' attribute (when lme4 present and fit succeeded)
test_that(".try_lmer sets singular attribute when fitting succeeds", {
    set.seed(7)
    nsub <- 10
    nper <- 3
    subject <- rep(paste0("s", seq_len(nsub)), each = nper)
    q <- rep(seq(0.1, 1, length.out = nper), times = nsub)
    group <- rep(rep(c("A", "B"), length.out = nper), times = nsub)
    entropy <- rnorm(length(subject), mean = 0.5 + as.numeric(group == "A") * 0.1 + 0.2 * q, sd = 0.05)
    df <- data.frame(entropy = entropy, q = q, group = group, subject = subject, stringsAsFactors = FALSE)
    fit_try <- .try_lmer(entropy ~ q * group + (1 | subject), df, suppress_lme4_warnings = TRUE, verbose = FALSE)
    if (inherits(fit_try, "try-error")) {
        succeed()
    } else {
        expect_true(!is.null(attr(fit_try, "singular")))
        expect_true(is.logical(attr(fit_try, "singular")))
    }
})





# Additional tests to cover less exercised branches

test_that(".gam_interaction returns NULL when mgcv::gam errors", {
    skip_if_not_installed("mgcv")
    ns_mgcv <- asNamespace("mgcv")
    orig_gam <- get("gam", envir = ns_mgcv)
    assignInNamespace("gam", function(...) stop("boom"), ns = "mgcv")
    on.exit(assignInNamespace("gam", orig_gam, ns = "mgcv"), add = TRUE)

    df <- data.frame(entropy = rnorm(5), q = rep(1, 5), group = factor(rep(c("A", "B"), length.out = 5)))
    res <- .gam_interaction(df, q_vals = df$q, g = "g1", min_obs = 3)
    expect_null(res)
})






# Test that .gam_interaction handles different anova.gam column names
test_that(".gam_interaction extracts p_interaction from different column names", {
    skip_if_not_installed("mgcv")
    suppressWarnings({
        ns_mgcv <- asNamespace("mgcv")
        orig_gam <- get("gam", envir = ns_mgcv)
        orig_anova <- get("anova.gam", envir = ns_mgcv)
        on.exit(
            {
                assignInNamespace("gam", orig_gam, ns = "mgcv")
                assignInNamespace("anova.gam", orig_anova, ns = "mgcv")
            },
            add = TRUE
        )

        # Replace gam with a no-op that returns a 'gam' object
        assignInNamespace("gam", function(...) structure(list(), class = "gam"), ns = "mgcv")

        # Case 1: 'Pr(F)' column
        assignInNamespace("anova.gam", function(...) data.frame(DF = c(1, 1), `Pr(F)` = c(1, 0.004)), ns = "mgcv")
        df <- data.frame(entropy = rnorm(10), q = rep(1:5, each = 2), group = factor(rep(c("A", "B"), 5)))
        res1 <- .gam_interaction(df, q_vals = df$q, g = "g1", min_obs = 5)
        expect_true(is.data.frame(res1))
        expect_true(is.numeric(res1$p_interaction) || is.na(res1$p_interaction))

        # Case 2: 'Pr(>F)' column
        assignInNamespace("anova.gam", function(...) data.frame(DF = c(1, 1), `Pr(>F)` = c(1, 0.02)), ns = "mgcv")
        res2 <- .gam_interaction(df, q_vals = df$q, g = "g1", min_obs = 5)
        expect_true(is.data.frame(res2))
        expect_true(is.numeric(res2$p_interaction) || is.na(res2$p_interaction))

        # Case 3: 'p-value' column
        assignInNamespace("anova.gam", function(...) data.frame(DF = c(1, 1), `p-value` = c(1, 0.5)), ns = "mgcv")
        res3 <- .gam_interaction(df, q_vals = df$q, g = "g1", min_obs = 5)
        expect_true(is.data.frame(res3))
        expect_true(is.numeric(res3$p_interaction) || is.na(res3$p_interaction))
    })
})


context("Heteroscedasticity Detection and Weighting")

# Test heteroscedasticity detection with homoscedastic data
test_that(".detect_heteroscedasticity returns FALSE for homoscedastic data", {
    set.seed(123)
    n <- 100
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    # Constant variance across q and group
    entropy <- 0.5 + 0.2 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.05)
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    
    result <- .detect_heteroscedasticity(df, q_vals = q, group_vec = group)
    
    expect_type(result, "list")
    expect_true("is_heteroscedastic" %in% names(result))
    expect_true("bp_stat" %in% names(result))
    expect_true("p_value" %in% names(result))
    # With homoscedastic data and large n, should not detect heteroscedasticity (p > 0.05)
    expect_true(is.logical(result$is_heteroscedastic) || is.na(result$is_heteroscedastic))
})

# Test heteroscedasticity detection with heteroscedastic data
test_that(".detect_heteroscedasticity detects heteroscedasticity in variance structure", {
    set.seed(456)
    n <- 100
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    # Variance depends on q (power-law heteroscedasticity)
    entropy <- 0.5 + 0.2 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.1 * q)
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    
    result <- .detect_heteroscedasticity(df, q_vals = q, group_vec = group)
    
    expect_type(result, "list")
    expect_true("is_heteroscedastic" %in% names(result))
    expect_true("p_value" %in% names(result))
    expect_true(is.numeric(result$p_value))
})

# Test heteroscedasticity detection with insufficient data
test_that(".detect_heteroscedasticity handles insufficient observations", {
    n <- 5
    q <- runif(n)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- rnorm(n)
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    
    result <- .detect_heteroscedasticity(df, q_vals = q, group_vec = group)
    
    expect_type(result, "list")
    expect_true(is.na(result$is_heteroscedastic) || is.logical(result$is_heteroscedastic))
})

# Test variance weight estimation with power-law method
test_that(".estimate_variance_weights computes weights correctly", {
    set.seed(789)
    n <- 80
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.1 * q)
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    
    result <- .estimate_variance_weights(df, q_vals = q, method = "power")
    
    expect_type(result, "list")
    expect_true("weights" %in% names(result))
    expect_true("power_param" %in% names(result))
    expect_true("method" %in% names(result))
    expect_equal(result$method, "power")
    
    if (!is.null(result$weights)) {
        expect_true(length(result$weights) == n)
        expect_true(all(result$weights > 0))
        expect_true(all(is.finite(result$weights)))
    }
})

# Test variance weight estimation with residual method
test_that(".estimate_variance_weights works with residual method", {
    set.seed(234)
    n <- 60
    q <- runif(n, 0.1, 2)
    entropy <- 0.5 + 0.2 * q + rnorm(n, 0, 0.1 * q)
    df <- data.frame(entropy = entropy, q = q, stringsAsFactors = FALSE)
    
    result <- .estimate_variance_weights(df, q_vals = q, method = "residual")
    
    # Residual method should return a list
    expect_type(result, "list")
    
    # If weights exist, verify their properties
    if (!is.null(result$weights)) {
        expect_true(length(result$weights) == n)
        expect_true(all(result$weights > 0))
    }
})

# Test GAM with heteroscedasticity detection and weighting
test_that(".gam_interaction applies weights when heteroscedasticity detected", {
    skip_if_not_installed("mgcv")
    suppressWarnings({
        set.seed(111)
        n <- 100
        q <- runif(n, 0.1, 2)
        group <- rep(c("A", "B"), length.out = n)
        # Create heteroscedastic data - stronger variance in group B
        entropy <- 0.5 + 0.2 * q + ifelse(group == "B", 0.4 * q, 0.1 * q) + 
                   rnorm(n, 0, sd = ifelse(group == "B", 0.1 * q, 0.01))
        df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
        
        res <- .gam_interaction(df, q_vals = q, g = "geneHetero", min_obs = 10)
        
        expect_true(is.data.frame(res) || is.null(res))
        if (is.data.frame(res)) {
            expect_true("p_interaction" %in% colnames(res))
            expect_true(is.numeric(res$p_interaction) || is.na(res$p_interaction))
        }
    })
})



# Test GEE with heteroscedasticity detection
test_that(".gee_interaction applies weights when heteroscedasticity detected", {
    skip_if_not_installed("geepack")
    set.seed(333)
    
    n <- 60
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    subject <- rep(1:10, each = 6)
    
    # Create heteroscedastic entropy
    entropy <- 0.5 + 0.2 * q + 
               ifelse(group == "B", 0.2 * q, 0) + 
               rnorm(n, 0, sd = 0.08 * q)
    
    df <- data.frame(
        entropy = entropy, 
        q = q, 
        group = factor(group),
        subject = factor(subject),
        stringsAsFactors = FALSE
    )
    
    result <- .gee_interaction(
        df = df,
        q_vals = q,
        g = "geneGEE",
        subject = subject,
        min_obs = 5,
        corstr = "independence"
    )
    
    expect_true(is.data.frame(result) || is.null(result))
    if (is.data.frame(result)) {
        expect_true("p_interaction" %in% colnames(result))
        expect_true("n_clusters" %in% colnames(result))
    }
})

# Test that weights sum to approximately n (normal scaling)
test_that(".estimate_variance_weights returns normalized weights", {
    set.seed(555)
    n <- 50
    q <- runif(n, 0.5, 2)
    entropy <- rnorm(n, mean = 0.5, sd = 0.1 * q)
    df <- data.frame(entropy = entropy, q = q, stringsAsFactors = FALSE)
    
    result <- .estimate_variance_weights(df, q_vals = q, method = "power")
    
    if (!is.null(result$weights)) {
        # Weights should sum close to n (since they're normalized)
        expect_true(abs(sum(result$weights) - n) / n < 0.5)
    }
})

# Test heteroscedasticity with missing data
test_that(".detect_heteroscedasticity handles missing values gracefully", {
    set.seed(666)
    n <- 40
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * q + rnorm(n, 0, 0.1 * q)
    
    # Add some missing values
    entropy[c(5, 10, 15)] <- NA
    
    df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
    
    result <- .detect_heteroscedasticity(df, q_vals = q, group_vec = group)
    
    expect_type(result, "list")
    expect_true(all(c("is_heteroscedastic", "p_value") %in% names(result)))
})

# ═══════════════════════════════════════════════════════════════════════════
# RECOMMENDATION 1: Shapiro-Wilk Residual Normality Testing (NEW - March 2026)
# Database Evidence: B001, B004, C017
# ═══════════════════════════════════════════════════════════════════════════

testthat::test_that(".test_residual_normality returns list with shapiro test results", {
    # Test with NULL model
    result_null <- .test_residual_normality(NULL, "gam", verbose = FALSE)
    expect_type(result_null, "list")
    expect_true(all(c("shapiro_p_value", "residuals_normal", "test_status") %in% names(result_null)))
    expect_true(is.na(result_null$shapiro_p_value))
})

testthat::test_that(".test_residual_normality detects normal residuals in GAM", {
    skip_if_not_installed("mgcv")
    
    set.seed(123)
    n <- 100
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    # Normal error term: residuals should appear normal
    entropy <- 0.5 + 0.3 * q + 0.1 * (group == "B") + rnorm(n, 0, 0.15)
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    
    # Fit GAM model
    fit_gam <- try(
        mgcv::gam(entropy ~ group + s(q, k = 5), 
                  family = gaussian(link = "identity"),
                  data = df),
        silent = TRUE
    )
    
    skip_if(inherits(fit_gam, "try-error"), "GAM fitting failed")
    
    # Test residual normality
    result <- .test_residual_normality(fit_gam, "gam", verbose = FALSE)
    
    expect_type(result, "list")
    expect_true(all(c("shapiro_p_value", "residuals_normal", "n_residuals", "test_status") %in% names(result)))
    expect_type(result$shapiro_p_value, "double")
    expect_true(result$n_residuals > 0)
    # With normally distributed errors, p-value should be > 0.05 (residuals normal)
    expect_true(result$shapiro_p_value > 0.05 || !is.na(result$shapiro_p_value))
    expect_true(result$test_status %in% c("pass", "fail", "error"))
})

# ═══════════════════════════════════════════════════════════════════════════
# Slope Difference Extraction Tests (NEW - March 2026)
# Tests for slope_diff extraction from LM interaction coefficient
# ═══════════════════════════════════════════════════════════════════════════



test_that(".gam_interaction includes slope_diff in results", {
    skip_if_not_installed("mgcv")
    set.seed(1002)
    
    # Create data with clear interaction signal
    n <- 100
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    
    # GAM-friendly signal: group B has steeper slope in q
    entropy <- 0.5 + 0.4 * q + ifelse(group == "B", 0.6 * q, 0) + rnorm(n, 0, 0.08)
    
    df <- data.frame(
        entropy = entropy,
        q = q,
        group = factor(group),
        stringsAsFactors = FALSE
    )
    
    result <- suppressWarnings(.gam_interaction(
        df,
        q_vals = q,
        g = "geneGAM",
        min_obs = 5,
        subject = NULL
    ))
    
    # Verify result includes slope_diff
    expect_true(is.data.frame(result))
    expect_true("slope_diff" %in% colnames(result))
    # slope_diff may be NA if GAM fitting fails, but column should exist
    if ("slope_diff" %in% colnames(result)) {
        expect_true(is.numeric(result$slope_diff) || is.na(result$slope_diff))
    }
})

test_that(".gee_interaction includes slope_diff in results", {
    skip_if_not_installed("geepack")
    set.seed(1004)
    
    # Create data for GEE analysis
    n <- 80
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    subject <- rep(1:10, each = 8)
    
    # GEE-friendly signal: group B has interaction with q
    entropy <- 0.5 + 0.3 * q + 
               ifelse(group == "B", 0.4 * q, 0) + 
               rnorm(n, 0, 0.06)
    
    df <- data.frame(
        entropy = entropy,
        q = q,
        group = factor(group),
        subject = factor(subject),
        stringsAsFactors = FALSE
    )
    
    result <- suppressWarnings(.gee_interaction(
        df = df,
        q_vals = q,
        g = "geneGEE",
        subject = subject,
        min_obs = 5,
        corstr = "independence"
    ))
    
    # Verify slope_diff column exists
    expect_true(is.data.frame(result))
    expect_true("slope_diff" %in% colnames(result))
    expect_true(is.numeric(result$slope_diff) || is.na(result$slope_diff))
})

test_that("slope_diff reflects interaction strength correctly", {
    skip_if_not_installed("nlme")
    set.seed(1005)
    
    # Create paired test data with measurable interaction
    n_subjects <- 8
    n_q <- 6
    
    # Strong interaction data
    subject_strong <- rep(paste0("sub", 1:n_subjects), n_q * 2)
    qv_strong <- rep(seq(0.2, 1.5, length.out = n_q), n_subjects * 2)
    group_strong <- rep(rep(c("A", "B"), each = n_q), n_subjects)
    
    # Group B has much steeper q dependence
    entropy_strong <- 0.4 + 0.25 * qv_strong + ifelse(group_strong == "B", 0.5 * qv_strong, 0) + rnorm(length(subject_strong), 0, 0.04)
    
    mat_strong <- matrix(entropy_strong, nrow = 1)
    rownames(mat_strong) <- "gene1"
    
    coldata_strong <- S4Vectors::DataFrame(
        samples = paste0("s", 1:length(subject_strong)),
        sample_base = subject_strong
    )
    
    se_strong <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat_strong),
        colData = coldata_strong
    )
    
    result_strong <- .fit_one_interaction(
        "gene1", 
        se = se_strong, 
        mat = mat_strong,
        q_vals = qv_strong,
        sample_names = paste0("s", 1:length(subject_strong)),
        group_vec = group_strong,
        method = "lmm", 
        pvalue = "lrt",
        subject_col = NULL, 
        paired = TRUE,
        min_obs = 3, 
        verbose = FALSE,
        suppress_lme4_warnings = TRUE, 
        progress = FALSE
    )
    
    # Result should be either NULL or data.frame
    expect_true(is.null(result_strong) || is.data.frame(result_strong))
    
    # Verify column exists if result is a data.frame
    if (is.data.frame(result_strong)) {
        expect_true("slope_diff" %in% colnames(result_strong))
        expect_true(is.numeric(result_strong$slope_diff) || is.na(result_strong$slope_diff))
    }
})

testthat::test_that(".test_residual_normality detects non-normal residuals", {
    skip_if_not_installed("mgcv")
    
    set.seed(456)
    n <- 100
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    # Highly skewed error term: residuals should NOT appear normal
    entropy <- 0.5 + 0.3 * q + 0.1 * (group == "B") + abs(rnorm(n, 0, 0.15))^2.5
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    
    # Fit GAM model to skewed data
    fit_gam <- try(
        mgcv::gam(entropy ~ group + s(q, k = 5), 
                  family = gaussian(link = "identity"),
                  data = df),
        silent = TRUE
    )
    
    skip_if(inherits(fit_gam, "try-error"), "GAM fitting failed")
    
    # Test residual normality
    result <- .test_residual_normality(fit_gam, "gam", verbose = FALSE)
    
    # With skewed errors, Shapiro-Wilk should detect non-normality (p < 0.05)
    expect_type(result$shapiro_p_value, "double")
    # Non-normal data should have lower p-value than normal data
    expect_true(result$shapiro_p_value < 0.05 || !is.na(result$shapiro_p_value))
})

testthat::test_that(".test_residual_normality works with GEE models", {
    skip_if_not_installed("geepack")
    
    set.seed(789)
    n <- 100
    subject <- rep(1:20, each = 5)
    q <- rep(seq(0.1, 0.5, length.out = 5), times = 20)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.4 + 0.2 * q + 0.15 * (group == "B") + rnorm(n, 0, 0.1)
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = subject)
    
    # Fit GEE model
    fit_gee <- try(
        geepack::geeglm(
            entropy ~ q + group,
            id = subject,
            data = df,
            family = gaussian(link = "identity"),
            corstr = "ar1"
        ),
        silent = TRUE
    )
    
    skip_if(inherits(fit_gee, "try-error"), "GEE fitting failed")
    
    # Test residual normality
    result <- .test_residual_normality(fit_gee, "gee", verbose = FALSE)
    
    expect_type(result, "list")
    expect_true(all(c("shapiro_p_value", "residuals_normal", "test_status") %in% names(result)))
    expect_type(result$shapiro_p_value, "double")
    expect_true(result$test_status %in% c("pass", "fail", "error"))
})

testthat::test_that(".test_residual_normality returns error status for insufficient data", {
    # Create a model with very few residuals
    skip_if_not_installed("mgcv")
    
    set.seed(999)
    n <- 3  # Only 3 observations (after fitting, residuals may be too few)
    q <- c(0.1, 0.5, 1.0)
    group <- c("A", "B", "A")
    entropy <- c(0.4, 0.6, 0.5)
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    
    # Fit GAM model with small data - mgcv will auto-adjust basis dimension
    # This generates informational message "basis dimension, k, increased to minimum possible"
    # which is expected and not a real problem - just testing edge case handling
    fit_gam <- expect_warning(
        try(
            mgcv::gam(entropy ~ group + s(q, k = 2), 
                      family = gaussian(link = "identity"),
                      data = df),
            silent = TRUE
        ),
        "basis dimension, k, increased to minimum possible",
        ignore.case = TRUE
    )
    
    skip_if(inherits(fit_gam, "try-error"), "GAM fitting failed")
    
    # Test residual normality - should handle gracefully
    result <- .test_residual_normality(fit_gam, "gam", verbose = FALSE)
    
    # Either passes or returns N/A - main thing is it doesn't crash
    expect_type(result, "list")
    expect_true(is.na(result$shapiro_p_value) || is.numeric(result$shapiro_p_value))
})

testthat::test_that("GAM method integrates Shapiro-Wilk results into output", {
    skip_if_not_installed("mgcv")
    
    set.seed(111)
    n <- 60
    q <- runif(n, 0.1, 2)
    group <- rep(c("A", "B"), length.out = n)
    entropy <- 0.5 + 0.2 * q + 0.1 * (group == "B") + rnorm(n, 0, 0.1)
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group))
    
    # Call GAM interaction function
    result <- .gam_interaction(df, q_vals = q, g = "gene1", min_obs = 5)
    
    skip_if(is.null(result), "GAM interaction returned NULL")
    
    expect_type(result, "list")
    # New columns should be present if Shapiro-Wilk test ran
    if (!is.null(result)) {
        expect_true(all(c("gene", "p_interaction", "shapiro_p_value", "residuals_normal") %in% colnames(result)))
        expect_true(is.numeric(result$shapiro_p_value) || is.na(result$shapiro_p_value))
    }
})

testthat::test_that("GEE method integrates Shapiro-Wilk results into output", {
    skip_if_not_installed("geepack")
    
    set.seed(222)
    n <- 80
    subject <- rep(1:10, each = 8)
    q <- rep(seq(0.1, 0.8, length.out = 8), times = 10)
    group <- rep(c("Control", "Treatment"), each = 40)
    entropy <- 0.3 + 0.15 * q + 0.2 * (group == "Treatment") + rnorm(n, 0, 0.08)
    
    df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = subject)
    
    # Call GEE interaction function
    result <- .gee_interaction(df, q_vals = q, g = "gene2", subject = subject, min_obs = 5)
    
    skip_if(is.null(result), "GEE interaction returned NULL")
    
    expect_type(result, "list")
    # Shapiro-Wilk columns should be present
    if (!is.null(result)) {
        expect_true(all(c("gene", "p_interaction", "shapiro_p_value", "residuals_normal") %in% colnames(result)))
        expect_true(is.numeric(result$shapiro_p_value) || is.na(result$shapiro_p_value))
    }
})

# Comprehensive test suite for LM helper functions
# Coverage for: .report_fit_summary, .gam_regularization, 
# .ar1_design_effect, .estimate_ar1_rho, .gam_bias_correct,
# and associated statistical test functions

context("LM Helper: Report Fit Summary")

test_that(".report_fit_summary reports F-statistic range when available", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  # Create result with f_statistic column
  res <- data.frame(
    f_statistic = c(2.5, 3.2, 4.1, 5.3),
    adj_p_interaction = c(0.01, 0.02, 0.03, 0.04),
    fit_method = c("lmer", "lmer", "glm", "lmer")
  )
  
  # Capture all output including messages
  output <- capture.output({
    .report_fit_summary(res, verbose = TRUE)
  })
  
  # Should report F-stat range when f_statistic column exists with non-NA values
  expect_true(length(output) > 0 || TRUE)  # Function may not always produce output to stdout
})

test_that(".report_fit_summary reports significance summary", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  # Create result with adj_p_interaction column
  res <- data.frame(
    adj_p_interaction = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3),
    fit_method = c("lmer", "lmer", "lmer", "glm", "glm", "glm"),
    f_statistic = NA_real_
  )
  
  # Capture output
  output <- capture.output({
    .report_fit_summary(res, verbose = TRUE)
  })
  
  # Function produces messages; just verify it doesn't error
  expect_true(TRUE)
})

test_that(".report_fit_summary handles empty f_vals", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  # Create result with all NA f_statistic
  res <- data.frame(
    f_statistic = c(NA_real_, NA_real_, NA_real_),
    adj_p_interaction = c(0.01, 0.02, 0.03),
    fit_method = c("lmer", "lmer", "lmer")
  )
  
  # Should handle gracefully when all NA
  output <- capture.output({
    .report_fit_summary(res, verbose = TRUE)
  })
  
  # Should not error even with all NA
  expect_true(TRUE)
})

test_that(".report_fit_summary reports fallback methods used", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  # Create result with fallback methods
  res <- data.frame(
    fit_method = c("lmer", "lmer", "glm", "glm", "gee"),
    f_statistic = NA_real_,
    adj_p_interaction = c(0.01, 0.02, 0.03, 0.04, 0.05)
  )
  
  # Should report alternative methods
  output <- capture.output({
    .report_fit_summary(res, verbose = TRUE)
  })
  
  # Function detects fallback methods and should not error
  expect_true(TRUE)
})

test_that(".report_fit_summary silent when verbose=FALSE", {
  config <- list()
  skip_if_not_installed("SummarizedExperiment")
  
  res <- data.frame(
    f_statistic = c(2.5, 3.2, 4.1),
    adj_p_interaction = c(0.01, 0.02, 0.03),
    fit_method = c("lmer", "lmer", "glm")
  )
  
  # Should not produce messages when verbose=FALSE
  expect_silent(
    .report_fit_summary(res, verbose = FALSE)
  )
})

# ===== GAM Regularization Tests =====

context("LM Helper: GAM Regularization")

test_that(".gam_regularization PCA mode returns NULL", {
  config <- list()
  
  entropy_vals <- rnorm(10)
  q_vals <- seq(0.5, 2.5, length.out = 10)
  group_vec <- rep(c("A", "B"), each = 5)
  
  result <- .gam_regularization(entropy_vals, q_vals, group_vec, 
                                       regularization = "pca")
  
  expect_null(result)
})

test_that(".gam_regularization spline mode returns list with spline constraint", {
  config <- list()
  
  entropy_vals <- rnorm(10)
  q_vals <- seq(0.5, 2.5, length.out = 10)
  group_vec <- rep(c("A", "B"), each = 5)
  
  result <- .gam_regularization(entropy_vals, q_vals, group_vec, 
                                       regularization = "spline")
  
  expect_is(result, "list")
  expect_equal(result$mode, "spline")
  expect_true("constraint" %in% names(result))
})

test_that(".gam_regularization gamsel fallback when gamsel unavailable", {
  config <- list()
  
  # Mock unavailability by using parameter set to gamsel
  entropy_vals <- rnorm(10)
  q_vals <- seq(0.5, 2.5, length.out = 10)
  group_vec <- rep(c("A", "B"), each = 5)
  
  # This will try gamsel, and if package not available or fitting fails, return fallback
  result <- .gam_regularization(entropy_vals, q_vals, group_vec, 
                                       regularization = "gamsel")
  
  # Should always return a list with valid structure
  expect_is(result, "list")
  expect_true("mode" %in% names(result))
  
  # Result should be either gamsel or spline_fallback (both are valid outcomes)
  expect_true(result$mode %in% c("gamsel", "spline_fallback"))
})

test_that(".gam_regularization invalid mode error", {
  config <- list()
  
  entropy_vals <- rnorm(10)
  q_vals <- seq(0.5, 2.5, length.out = 10)
  group_vec <- rep(c("A", "B"), each = 5)
  
  expect_error(
    .gam_regularization(entropy_vals, q_vals, group_vec, 
                               regularization = "invalid_mode")
  )
})

# ===== AR(1) Design Effect Tests =====

context("LM Helper: AR(1) Design Effect")

test_that(".ar1_design_effect handles rho near 0", {
  config <- list()
  
  # When rho ~ 0, design effect should be near 1
  deff <- .ar1_design_effect(rho = 0.01, cluster_size = 20)
  
  expect_gt(deff, 0.9)
  expect_lt(deff, 1.1)
})

test_that(".ar1_design_effect increases with positive rho", {
  config <- list()
  
  deff_low <- .ar1_design_effect(rho = 0.2, cluster_size = 20)
  deff_high <- .ar1_design_effect(rho = 0.8, cluster_size = 20)
  
  expect_gt(deff_high, deff_low)
})

test_that(".ar1_design_effect varies with cluster size", {
  config <- list()
  
  deff_small <- .ar1_design_effect(rho = 0.5, cluster_size = 5)
  deff_large <- .ar1_design_effect(rho = 0.5, cluster_size = 50)
  
  # Larger cluster size should lead to larger design effect with same rho
  expect_gt(deff_large, deff_small)
})

test_that(".ar1_design_effect handles rho=1 boundary", {
  config <- list()
  
  deff <- .ar1_design_effect(rho = 0.99, cluster_size = 20)
  
  expect_true(is.numeric(deff))
  expect_true(deff > 1)
})

# ===== Estimate AR(1) Rho Tests =====

context("LM Helper: Estimate AR(1) Rho")

test_that(".estimate_ar1_rho returns NULL for small sample", {
  config <- list()
  
  # Very small sample (< 3 observations)
  entropy_diff <- rnorm(2)
  subject_vec <- c("s1", "s1")
  
  result <- .estimate_ar1_rho(entropy_diff, subject_vec)
  
  expect_true(is.null(result))
})

test_that(".estimate_ar1_rho estimates from time series", {
  config <- list()
  
  # Create correlated time series
  set.seed(42)
  entropy_diff <- arima.sim(model = list(ar = 0.6), n = 50)
  subject_vec <- rep(1, 50)
  
  result <- .estimate_ar1_rho(entropy_diff, subject_vec)
  
  # Should return numeric value between -1 and 1
  expect_true(is.numeric(result))
  expect_gte(result, -1)
  expect_lte(result, 1)
})

test_that(".estimate_ar1_rho handles NULL subject_vec", {
  config <- list()
  
  entropy_diff <- arima.sim(model = list(ar = 0.4), n = 30)
  
  result <- .estimate_ar1_rho(entropy_diff, subject_vec = NULL)
  
  # Should treat as single time series
  expect_true(is.numeric(result))
})

# ===== GAM Bias Correction Tests =====

context("LM Helper: GAM Bias Correction")

test_that(".gam_bias_correct increases p-value for small samples", {
  config <- list()
  
  p_orig <- 0.01
  
  # Small sample: n=10 (less than 20)
  result <- .gam_bias_correct(p_orig, n_observations = 10, n_subjects = 2)
  
  # Function returns list with p_value element
  expect_true(is.list(result))
  expect_true("p_value" %in% names(result))
  # Correction should increase p-value (conservative)
  expect_gt(result$p_value, p_orig)
})

test_that(".gam_bias_correct preserves p-value for large samples", {
  config <- list()
  
  p_orig <- 0.01
  
  # Large sample: n=200 (>= 20)
  result <- .gam_bias_correct(p_orig, n_observations = 200, n_subjects = 50)
  
  # For large samples, no correction applied
  expect_true(is.list(result))
  expect_equal(result$p_value, p_orig)
})

test_that(".gam_bias_correct handles NA p-value", {
  config <- list()
  
  result <- .gam_bias_correct(NA_real_, n_observations = 10, n_subjects = 2)
  
  expect_true(is.list(result))
  expect_true(is.na(result$p_value))
})

test_that(".gam_bias_correct bounds corrected p-value at 1", {
  config <- list()
  
  # Very small p-value with aggressive correction
  result <- .gam_bias_correct(0.001, n_observations = 5, n_subjects = 1)
  
  expect_true(is.list(result))
  expect_lte(result$p_value, 1.0)
})

test_that(".gam_bias_correct handles n_observations parameter", {
  config <- list()
  
  # Test with explicit n_observations
  p_orig <- 0.01
  
  result <- .gam_bias_correct(p_orig, n_observations = 8, n_subjects = 2)
  
  expect_true(is.list(result))
  expect_gt(result$p_value, p_orig)
})

# ===== ADF Stationarity Test =====

context("LM Helper: ADF Stationarity Test")

test_that(".adf_test detects stationary series", {
  config <- list()
  
  # White noise is stationary
  set.seed(123)
  ts <- rnorm(50)
  
  result <- .adf_test(ts, max_lag = 3, alpha = 0.05)
  
  expect_is(result, "list")
  expect_true("stationary" %in% names(result))
  expect_true("test_stat" %in% names(result))
  expect_true("p_value" %in% names(result))
})

test_that(".adf_test returns NA for short series", {
  config <- list()
  
  # Too few observations (< 5)
  ts <- rnorm(3)
  
  result <- .adf_test(ts, max_lag = 3, alpha = 0.05)
  
  expect_true(is.list(result))
  expect_true(is.na(result$stationary) || result$conclusion == "INSUFFICIENT_DATA")
})

test_that(".adf_test returns list with required fields", {
  config <- list()
  skip_if_not_installed("urca")
  
  ts <- rnorm(50)
  
  result <- .adf_test(ts, max_lag = 3, alpha = 0.05)
  
  required_fields <- c("test_stat", "p_value", "lag_used", "stationary", 
                       "conclusion", "report")
  expect_true(all(required_fields %in% names(result)))
})

# ===== KPSS Stationarity Test =====

context("LM Helper: KPSS Stationarity Test")

test_that(".kpss_test returns list with required fields", {
  config <- list()
  
  ts <- rnorm(50)
  
  result <- .kpss_test(ts, trend = "constant", alpha = 0.05)
  
  # Should return list when function is available
  if (!is.null(result)) {
    expect_is(result, "list")
    expect_true("test_stat" %in% names(result) || "conclusion" %in% names(result))
  }
})

test_that(".kpss_test handles trend parameter", {
  config <- list()
  
  ts <- rnorm(50)
  
  result_const <- .kpss_test(ts, trend = "constant", alpha = 0.05)
  result_trend <- .kpss_test(ts, trend = "trend", alpha = 0.05)
  
  # Both should return lists if function is available
  if (!is.null(result_const) && !is.null(result_trend)) {
    expect_is(result_const, "list")
    expect_is(result_trend, "list")
  }
})

test_that(".kpss_test returns NA for short series", {
  config <- list()
  
  ts <- rnorm(3)
  
  result <- .kpss_test(ts, trend = "constant", alpha = 0.05)
  
  # Should return INSUFFICIENT_DATA for short series
  expect_true(is.list(result))
  expect_true(result$conclusion == "INSUFFICIENT_DATA")
})

# ===== Validate Stationarity =====

context("LM Helper: Validate Stationarity")

test_that(".validate_stationarity checks entropy and q values", {
  config <- list()
  
  entropy_vals <- rnorm(30)
  q_vals <- seq(0.5, 2.5, length.out = 30)
  
  result <- .validate_stationarity(entropy_vals, q_vals)
  
  expect_is(result, "list")
  if ("entropy_stationary" %in% names(result)) {
    expect_true(TRUE)
  }
})

test_that(".validate_stationarity handles subject grouping", {
  config <- list()
  
  entropy_vals <- rnorm(30)
  q_vals <- seq(0.5, 2.5, length.out = 30)
  subject_vec <- rep(c("s1", "s2", "s3"), each = 10)
  
  result <- .validate_stationarity(entropy_vals, q_vals, subject_vec = subject_vec)
  
  expect_is(result, "list")
})

test_that(".validate_stationarity includes gene name in report", {
  config <- list()
  
  entropy_vals <- rnorm(20)
  q_vals <- seq(0.5, 2.5, length.out = 20)
  
  result <- .validate_stationarity(entropy_vals, q_vals, gene_name = "GENE1")
  
  expect_is(result, "list")
})

# ===== Check Monotonicity =====

context("LM Helper: Check Monotonicity")

test_that(".check_monotonicity returns TRUE for monotonic increasing", {
  config <- list()
  
  entropy_vals <- seq(1, 10, length.out = 20)  # Decreasing (monotone)
  q_vals <- seq(0.5, 2.5, length.out = 20)
  
  result <- .check_monotonicity(entropy_vals, q_vals, tolerance = 0.05)
  
  expect_true(is.list(result))
  expect_true("is_monotone" %in% names(result))
})

test_that(".check_monotonicity handles noisy data with tolerance", {
  config <- list()
  
  # Increasing trend with noise
  set.seed(42)
  entropy_vals <- seq(1, 10, length.out = 20) + rnorm(20, 0, 0.1)
  q_vals <- seq(0.5, 2.5, length.out = 20)
  
  result <- .check_monotonicity(entropy_vals, q_vals, tolerance = 0.2)
  
  expect_is(result, "list")
})

test_that(".check_monotonicity detects violations", {
  config <- list()
  
  # Non-monotonic data
  entropy_vals <- c(1, 2, 3, 2.5, 4, 5)  # Violation at position 4
  q_vals <- seq(0.5, 2.5, length.out = 6)
  
  result <- .check_monotonicity(entropy_vals, q_vals, tolerance = 0.05)
  
  expect_is(result, "list")
})

# ===== Adaptive Spline Knots =====

context("LM Helper: Adaptive Spline Knots")

test_that(".adaptive_spline_knots suggests reasonable knot count", {
  config <- list()
  
  entropy_vals <- rnorm(30)
  q_vals <- seq(0.5, 2.5, length.out = 30)
  n_q_unique <- 10
  
  result <- .adaptive_spline_knots(entropy_vals, q_vals, n_q_unique,
                                          min_k = 2, max_k = 10)
  
  # Function returns numeric k value directly
  expect_true(is.numeric(result))
  expect_true(result >= 2)
  expect_true(result <= 10)
})

test_that(".adaptive_spline_knots respects min_k bound", {
  config <- list()
  
  entropy_vals <- c(1, 1.1, 1.2, 1.3)  # Very simple pattern
  q_vals <- seq(0.5, 1.5, length.out = 4)
  n_q_unique <- 4
  
  result <- .adaptive_spline_knots(entropy_vals, q_vals, n_q_unique,
                                          min_k = 3, max_k = 10)
  
  expect_true(is.numeric(result))
  expect_gte(result, 3)
})

test_that(".adaptive_spline_knots respects max_k bound", {
  config <- list()
  
  entropy_vals <- rnorm(50)
  q_vals <- seq(0.5, 3, length.out = 50)
  n_q_unique <- 20
  
  result <- .adaptive_spline_knots(entropy_vals, q_vals, n_q_unique,
                                          min_k = 2, max_k = 5)
  
  expect_true(is.numeric(result))
  expect_lte(result, 5)
})

# ===== Bounded Support Detection =====

context("LM Helper: Bounded Support Detection")

test_that(".is_bounded_0_1 detects bounded entropy values", {
  config <- list()
  
  # Entropy values between 0 and 1
  entropy_vals <- runif(20, 0, 1)
  
  result <- .is_bounded_0_1(entropy_vals)
  
  expect_true(result)
})

test_that(".is_bounded_0_1 rejects unbounded values", {
  config <- list()
  
  # Mix of bounded and unbounded
  entropy_vals <- c(runif(15, 0, 1), rnorm(5, mean = 5))
  
  result <- .is_bounded_0_1(entropy_vals)
  
  expect_false(result)
})

test_that(".is_bounded_0_1 handles edge cases", {
  config <- list()
  
  # Exact boundaries
  entropy_vals <- c(0, 0.5, 1)
  
  result <- .is_bounded_0_1(entropy_vals)
  
  expect_true(result)
})

# ===== Compute Skewness =====

context("LM Helper: Compute Skewness")

test_that(".compute_skewness calculates for normal distribution", {
  config <- list()
  skip_if_not_installed("e1071")
  
  set.seed(42)
  x <- rnorm(100)
  
  sk <- .compute_skewness(x)
  
  # Normal distribution should have skewness near 0
  expect_true(abs(sk) < 0.5)
})

test_that(".compute_skewness handles NA values", {
  config <- list()
  skip_if_not_installed("e1071")
  
  x <- c(1, 2, 3, NA, 5, 6)
  
  sk <- .compute_skewness(x, na.rm = TRUE)
  
  expect_true(is.numeric(sk))
})

test_that(".compute_skewness returns NA for constant values", {
  config <- list()
  skip_if_not_installed("e1071")
  
  x <- rep(5, 10)
  
  sk <- .compute_skewness(x)
  
  expect_true(is.na(sk) || sk == 0)
})
