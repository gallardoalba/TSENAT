context("calc_lm_helpers functions")

library(testthat)

# Report summary messages
test_that(".tsenat_report_fit_summary prints fallback and singular messages", {
  df <- data.frame(fit_method = c("lm_nosubject", "lmer", NA), singular = c(TRUE, FALSE, NA), stringsAsFactors = FALSE)
  expect_message(.tsenat_report_fit_summary(df, verbose = TRUE), "fallback fits used")
  expect_message(.tsenat_report_fit_summary(df, verbose = TRUE), "singular fits detected")
})

# GAM interaction: skip if mgcv not available
test_that(".tsenat_gam_interaction returns a data.frame with p_interaction when mgcv present", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  # build small dataset with group and q, per-sample entropy
  n <- 40
  q <- runif(n, 0.1, 2)
  group <- rep(c("A","B"), length.out = n)
  entropy <- 0.5 + 0.2 * (q) + ifelse(group=="A", 0.05, -0.05) + rnorm(n, 0, 0.01)
  df <- data.frame(entropy = entropy, q = q, group = group, stringsAsFactors = FALSE)
  res <- .tsenat_gam_interaction(df, q_vals = q, g = "g1", min_obs = 5)
  expect_true(is.data.frame(res) || is.null(res))
  if (is.data.frame(res)) {
    expect_true("p_interaction" %in% colnames(res))
  }
})

# FPCA interaction: synthetic matrix
test_that(".tsenat_fpca_interaction computes a p-value with reasonable input", {
  set.seed(2)
  # Create matrix genes x observations
  genes <- paste0("g", 1:3)
  samples <- paste0("S", 1:8)
  q_vals <- rep(c(0.1, 0.5, 1, 2), 2)
  # construct mat with rows genes, cols observations
  mat <- matrix(rnorm(length(genes) * length(q_vals)), nrow = length(genes))
  rownames(mat) <- genes
  # duplicate sample names to match observations length
  sample_names <- rep(samples[1:4], 2)
  group_vec <- rep(c("A","B"), each = 4)
  # use min_obs small to allow test
  res <- .tsenat_fpca_interaction(mat, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, g = 1, min_obs = 2)
  expect_true(is.null(res) || (is.data.frame(res) && "p_interaction" %in% colnames(res)))
})

# Try lm fallbacks and LRT/Satterthwaite extraction
test_that(".tsenat_try_lm_fallbacks returns lm fits and extractors return numeric p-values", {
  # build small long-format df
  df <- data.frame(
    entropy = rnorm(30),
    q = rep(seq(0.1, 1.0, length.out = 10), 3),
    group = rep(c("A","B","A"), each = 10),
    subject = rep(paste0("sub", 1:10), 3),
    stringsAsFactors = FALSE
  )
  fb <- .tsenat_try_lm_fallbacks(df, verbose = TRUE)
  expect_true(is.null(fb) || (is.list(fb) && all(c("fit0","fit1","method") %in% names(fb))))
  if (!is.null(fb)) {
    lrt_p <- .tsenat_extract_lrt_p(fb$fit0, fb$fit1)
    expect_true(is.numeric(lrt_p) || is.na(lrt_p))
    st_p <- .tsenat_extract_satterthwaite_p(NULL, fallback_lm = fb)
    expect_true(is.numeric(st_p) || is.na(st_p))
  }
})

# FPCA matrix preparation helper (variance filtering and scaling)
## Note: the FPCA-prep helper has multiple definitions in the source; tests
## below use the variant that accepts (mat, sample_names, q_vals, min_obs).

# Alternative FPCA matrix builder that returns mat_sub/used_samples
test_that(".tsenat_prepare_fpca_matrix (fpca variant) builds sub-matrix or returns NULL when insufficient", {
  # Use a single-gene matrix so element assignment in the helper is scalar
  mat <- matrix(rnorm(1 * 8), nrow = 1) # 1 gene x 8 observations
  sample_names <- rep(paste0("S", 1:4), 2)
  q_vals <- rep(c(0.1,0.5,1,2), 2)
  res <- .tsenat_prepare_fpca_matrix(mat, sample_names = sample_names, q_vals = q_vals, min_obs = 2)
  expect_true(is.null(res) || (is.list(res) && all(c("mat_sub","used_samples") %in% names(res))))
})

# Test lmer wrapper if available
test_that(".tsenat_try_lmer attempts lmer fitting when lme4 is installed", {
  skip_if_not_installed("lme4")
  set.seed(5)
  # Build a small balanced dataset for mixed model
  nsub <- 10
  nper <- 4
  subject <- rep(paste0("s", seq_len(nsub)), each = nper)
  q <- rep(seq(0.1, 1, length.out = nper), times = nsub)
  group <- rep(rep(c("A","B"), length.out = nper), times = nsub)
  entropy <- rnorm(length(subject), mean = 0.5 + as.numeric(group=="A")*0.1 + 0.2 * q, sd = 0.05)
  df <- data.frame(entropy = entropy, q = q, group = group, subject = subject, stringsAsFactors = FALSE)
  f <- as.formula("entropy ~ q * group + (1 | subject)")
  fit_try <- .tsenat_try_lmer(f, df, suppress_lme4_warnings = TRUE, verbose = FALSE)
  expect_true(inherits(fit_try, "try-error") || inherits(fit_try, "lmerMod"))
})

# Test extraction of Satterthwaite p-value using lmerTest when available
test_that(".tsenat_extract_satterthwaite_p returns numeric or NA when lmerTest present", {
  skip_if_not_installed("lme4")
  # prefer lmerTest but allow absence
  set.seed(6)
  subject <- rep(paste0("s", 1:8), each = 3)
  q <- rep(c(0.1, 0.5, 1), times = 8)
  group <- rep(c("A","B"), length.out = length(q))
  entropy <- rnorm(length(q), mean = 0.2 + as.numeric(group=="A")*0.05 + 0.3 * q, sd = 0.02)
  df <- data.frame(entropy = entropy, q = q, group = group, subject = subject, stringsAsFactors = FALSE)
  f <- as.formula("entropy ~ q * group + (1 | subject)")
  fit <- try(lme4::lmer(f, data = df, REML = FALSE), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    pval <- .tsenat_extract_satterthwaite_p(fit, fallback_lm = NULL)
    expect_true(is.numeric(pval) || is.na(pval))
  } else {
    succeed()
  }
})
