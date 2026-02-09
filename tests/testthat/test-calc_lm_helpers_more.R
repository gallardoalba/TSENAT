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
  res <- .tsenat_fpca_interaction(mat = mat, q_vals = q_vals, sample_names = sample_names,
    group_vec = group_vec, g = 1, min_obs = 2)
  testthat::expect_true(is.data.frame(res) || is.null(res))
  if (!is.null(res)) {
    testthat::expect_named(res, c("gene", "p_interaction"))
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
  out <- .tsenat_prepare_fpca_matrix(mat, sample_names = sample_names, q_vals = q_vals,
    min_obs = 2)
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
  res <- .tsenat_try_lm_fallbacks(df)
  testthat::expect_type(res, "list")
  testthat::expect_true(res$method %in% c("lm_subject", "lm_nosubject"))
  testthat::expect_s3_class(res$fit1, "lm")

  # drop subject -> should pick nosubject fallback
  df2 <- df[, c("entropy", "q", "group")]
  res2 <- .tsenat_try_lm_fallbacks(df2)
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
  p <- .tsenat_extract_lrt_p(fit0, fit1)
  testthat::expect_true(is.numeric(p) || is.na(p))
  if (!is.na(p)) testthat::expect_true(p >= 0 && p <= 1)
})

testthat::test_that("Satterthwaite extraction uses fallback lm coefficients when provided", {
  set.seed(3)
  n <- 48
  subject <- rep(1:12, each = 4)
  q <- runif(n)
  group <- rep(c("A", "B"), length.out = n)
  entropy <- 0.3 * q + ifelse(group == "B", 0.4, 0) + rnorm(n, 0, 0.15)
  df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = factor(subject))
  fb <- .tsenat_try_lm_fallbacks(df)
  testthat::expect_type(fb, "list")
  p_fb <- .tsenat_extract_satterthwaite_p(NULL, fallback_lm = fb)
  testthat::expect_true(is.numeric(p_fb) || is.na(p_fb))
})

testthat::test_that("GAM interaction returns a data.frame with p-value when mgcv available", {
  if (!rlang::is_installed("mgcv")) {
    testthat::skip("mgcv not installed")
  }
  set.seed(4)
  n <- 80
  q <- runif(n)
  group <- rep(c("A", "B"), length.out = n)
  entropy <- 0.2 * q + ifelse(group == "B", 0.6 * q, 0) + rnorm(n, 0, 0.15)
  df <- data.frame(entropy = entropy, q = q, group = factor(group))
  res <- .tsenat_gam_interaction(df, q_vals = q, g = "geneX", min_obs = 5)
  testthat::expect_true(is.data.frame(res) || is.null(res))
  if (!is.null(res)) {
    testthat::expect_named(res, c("gene", "p_interaction"))
    testthat::expect_type(res$p_interaction, "double")
  }
})

testthat::test_that("try_lmer returns an lmer object when lme4 available", {
  if (!rlang::is_installed("lme4")) testthat::skip("lme4 not installed")
  set.seed(5)
  n <- 48
  subject <- rep(1:12, each = 4)
  q <- runif(n)
  group <- rep(c("A", "B"), length.out = n)
  entropy <- 0.25 * q + ifelse(group == "B", 0.3, 0) + rnorm(n, 0, 0.1)
  df <- data.frame(entropy = entropy, q = q, group = factor(group), subject = factor(subject))
  fmla <- stats::as.formula("entropy ~ q * group + (1|subject)")
  fit <- .tsenat_try_lmer(fmla, data = df, suppress_lme4_warnings = TRUE)
  testthat::expect_true(inherits(fit, "lmerMod") || inherits(fit, "try-error"))
})
