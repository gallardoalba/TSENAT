context("Bootstrap Divergence: Refactored Helper Functions Tests")

# Suppress nboot warnings for this test file
options(TSENAT.suppress_nboot_warning = TRUE)

# Safety check: ensure test factory is loaded
if (!exists("create_test_se_simple", mode = "function")) {
  factory_file <- file.path(dirname(getwd()), "testthat", "tests-factory.R")
  if (!file.exists(factory_file)) {
    factory_file <- "tests/testthat/tests-factory.R"
  }
  if (file.exists(factory_file)) {
    source(factory_file, local = FALSE)
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 1: .bootstrap_divergence_validate_inputs()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_validate_inputs accepts valid inputs", {
  expect_no_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80, 60, 40),
      y = c(90, 70, 50, 40),
      ci = 0.95,
      q = 1
    )
  )
})

test_that(".bootstrap_divergence_validate_inputs rejects NULL x or y", {
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = NULL,
      y = c(100, 80),
      ci = 0.95,
      q = 1
    ),
    "Must provide"
  )
  
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80),
      y = NULL,
      ci = 0.95,
      q = 1
    ),
    "Must provide"
  )
})

test_that(".bootstrap_divergence_validate_inputs rejects non-numeric inputs", {
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c("a", "b"),
      y = c(100, 80),
      ci = 0.95,
      q = 1
    ),
    "numeric"
  )
})

test_that(".bootstrap_divergence_validate_inputs rejects negative counts", {
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, -80, 60),
      y = c(90, 70, 50),
      ci = 0.95,
      q = 1
    ),
    "non-negative"
  )
})

test_that(".bootstrap_divergence_validate_inputs rejects invalid ci", {
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80),
      y = c(90, 70),
      ci = 0,
      q = 1
    ),
    "ci must be between"
  )
  
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80),
      y = c(90, 70),
      ci = 1,
      q = 1
    ),
    "ci must be between"
  )
})

test_that(".bootstrap_divergence_validate_inputs rejects non-positive q", {
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80),
      y = c(90, 70),
      ci = 0.95,
      q = -0.5
    ),
    "q must be"
  )
  
  expect_error(
    TSENAT:::.bootstrap_divergence_validate_inputs(
      x = c(100, 80),
      y = c(90, 70),
      ci = 0.95,
      q = 0
    ),
    "q must be"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 2: .bootstrap_divergence_extract_data()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_extract_data returns vectors when x and y provided", {
  result <- TSENAT:::.bootstrap_divergence_extract_data(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    se = NULL,
    res = NULL,
    top_n = 1,
    group_col = "group",
    control_group = "Normal"
  )
  
  expect_equal(result$x, c(100, 80, 60))
  expect_equal(result$y, c(90, 70, 50))
  expect_null(result$gene_name)
})

test_that(".bootstrap_divergence_extract_data extracts from SummarizedExperiment", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(100, 80, 60, 40, 90, 70, 50, 40),
      nrow = 2,
      ncol = 4,
      dimnames = list(c("GENE1", "GENE2"), NULL)
    )),
    colData = data.frame(group = c("Normal", "Normal", "Treatment", "Treatment"))
  )
  
  res <- data.frame(row.names = c("GENE1", "GENE2"), pval = c(0.01, 0.5))
  
  result <- TSENAT:::.bootstrap_divergence_extract_data(
    x = NULL,
    y = NULL,
    se = se,
    res = res,
    top_n = 1,
    group_col = "group",
    control_group = "Normal"
  )
  
  # Matrix fills column-first: [100,80,60,40,90,70,50,40] in 2x4 = [[100,60,90,50],[80,40,70,40]]
  # GENE1 is row 1, samples 1-2 are Normal, 3-4 are Treatment
  expect_equal(result$x, c(100, 60))  # GENE1 Normal group (row 1, cols 1-2)
  expect_equal(result$y, c(90, 50))   # GENE1 Treatment group (row 1, cols 3-4)
})

test_that(".bootstrap_divergence_extract_data rejects out-of-range top_n", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(100, 80, 60, 40, 90, 70, 50, 40),
      nrow = 2,
      ncol = 4,
      dimnames = list(c("GENE1", "GENE2"), NULL)
    )),
    colData = data.frame(group = c("Normal", "Normal", "Treatment", "Treatment"))
  )
  
  res <- data.frame(row.names = c("GENE1", "GENE2"), pval = c(0.01, 0.5))
  
  expect_error(
    TSENAT:::.bootstrap_divergence_extract_data(
      x = NULL,
      y = NULL,
      se = se,
      res = res,
      top_n = 100,
      group_col = "group",
      control_group = "Normal"
    ),
    "out of range"
  )
})

test_that(".bootstrap_divergence_extract_data rejects non-SE object", {
  expect_error(
    TSENAT:::.bootstrap_divergence_extract_data(
      x = NULL,
      y = NULL,
      se = list(data = "not_se"),
      res = data.frame(),
      top_n = 1,
      group_col = "group",
      control_group = "Normal"
    ),
    "SummarizedExperiment"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 3: .bootstrap_divergence_compute_ci()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_compute_ci returns percentile CI correctly", {
  set.seed(42)
  bootstrap_dist <- rnorm(1000, mean = 0.5, sd = 0.1)
  
  result <- TSENAT:::.bootstrap_divergence_compute_ci(
    valid_divs = bootstrap_dist,
    estimate = 0.5,
    method = "percentile",
    nboot = 1000,
    ci = 0.95
  )
  
  expect_true(is.numeric(result$lower))
  expect_true(is.numeric(result$upper))
  expect_true(result$lower < result$upper)
  expect_true(result$lower > 0.2)  # Expected approximate range
  expect_true(result$upper < 0.8)
})

test_that(".bootstrap_divergence_compute_ci handles nboot=0", {
  result <- TSENAT:::.bootstrap_divergence_compute_ci(
    valid_divs = numeric(0),
    estimate = 0.5,
    method = "percentile",
    nboot = 0,
    ci = 0.95
  )
  
  expect_true(is.na(result$lower))
  expect_true(is.na(result$upper))
})

test_that(".bootstrap_divergence_compute_ci returns BCA CI correctly", {
  set.seed(42)
  bootstrap_dist <- rbeta(1000, shape1 = 2, shape2 = 5)
  
  result <- TSENAT:::.bootstrap_divergence_compute_ci(
    valid_divs = bootstrap_dist,
    estimate = mean(bootstrap_dist),
    method = "bca",
    nboot = 1000,
    ci = 0.95
  )
  
  expect_true(is.numeric(result$lower))
  expect_true(is.numeric(result$upper))
  expect_true(result$lower < result$upper)
})

test_that(".bootstrap_divergence_compute_ci respects different CI levels", {
  set.seed(42)
  bootstrap_dist <- rnorm(1000, mean = 0.5, sd = 0.1)
  
  ci_90 <- TSENAT:::.bootstrap_divergence_compute_ci(
    valid_divs = bootstrap_dist,
    estimate = 0.5,
    method = "percentile",
    nboot = 1000,
    ci = 0.90
  )
  
  ci_99 <- TSENAT:::.bootstrap_divergence_compute_ci(
    valid_divs = bootstrap_dist,
    estimate = 0.5,
    method = "percentile",
    nboot = 1000,
    ci = 0.99
  )
  
  # 99% CI should be wider than 90% CI
  expect_true((ci_99$upper - ci_99$lower) > (ci_90$upper - ci_90$lower))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 4: .bootstrap_divergence_print_single()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_print_single prints without error", {
  result <- list(
    estimate = 0.5,
    lower_ci = 0.4,
    upper_ci = 0.6,
    ci_level = 0.95,
    method = "percentile",
    nboot = 1000,
    bootstrap_dist = rnorm(1000),
    q = 1,
    gene_name = "TEST_GENE"
  )
  
  expect_no_error(
    TSENAT:::.bootstrap_divergence_print_single(
      result = result,
      gene_name = "TEST_GENE",
      q = 1,
      nboot = 1000,
      method = "percentile",
      ci = 0.95
    )
  )
})

test_that(".bootstrap_divergence_print_single handles various q values", {
  result <- list(
    estimate = 0.25,
    lower_ci = 0.15,
    upper_ci = 0.35
  )
  
  for (q_val in c(0.5, 1, 1.5, 2)) {
    expect_no_error(
      TSENAT:::.bootstrap_divergence_print_single(
        result = result,
        gene_name = "GENE_X",
        q = q_val,
        nboot = 500,
        method = "percentile",
        ci = 0.95
      )
    )
  }
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 5: .bootstrap_divergence_compute_paired_samples()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_compute_paired_samples extracts pair_id correctly", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      rpois(100, 10),
      nrow = 5,
      ncol = 20
    )),
    colData = data.frame(
      group = rep(c("Control", "Treatment"), 10),
      pair_id = rep(1:10, each = 2)
    )
  )
  
  result <- TSENAT:::.bootstrap_divergence_compute_paired_samples(
    se = se,
    pair_id_col = "pair_id",
    group_col = "group",
    control_group = "Control",
    x = numeric(10),
    y = numeric(10)
  )
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 20)
  expect_equal(length(unique(result)), 10)  # 10 unique pair IDs
})

test_that(".bootstrap_divergence_compute_paired_samples auto-detects pair_id", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      rpois(100, 10),
      nrow = 5,
      ncol = 20
    )),
    colData = data.frame(
      group = rep(c("Control", "Treatment"), 10),
      subject_id = rep(1:10, each = 2)
    )
  )
  
  result <- TSENAT:::.bootstrap_divergence_compute_paired_samples(
    se = se,
    pair_id_col = NULL,
    group_col = "group",
    control_group = "Control",
    x = numeric(10),
    y = numeric(10)
  )
  
  expect_true(is.numeric(result))
  expect_equal(length(unique(result)), 10)
})

test_that(".bootstrap_divergence_compute_paired_samples rejects insufficient pairs", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      rpois(20, 10),
      nrow = 5,
      ncol = 4
    )),
    colData = data.frame(
      group = c("Control", "Treatment", "Control", "Treatment"),
      pair_id = c(1, 1, 1, 1)
    )
  )
  
  expect_error(
    TSENAT:::.bootstrap_divergence_compute_paired_samples(
      se = se,
      pair_id_col = "pair_id",
      group_col = "group",
      control_group = "Control",
      x = numeric(2),
      y = numeric(2)
    ),
    "requires at least 2 pairs"
  )
})

test_that(".bootstrap_divergence_compute_paired_samples missing pair_id column", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      rpois(50, 10),
      nrow = 5,
      ncol = 10
    )),
    colData = data.frame(
      group = rep(c("Control", "Treatment"), 5)
    )
  )
  
  expect_error(
    TSENAT:::.bootstrap_divergence_compute_paired_samples(
      se = se,
      pair_id_col = NULL,
      group_col = "group",
      control_group = "Control",
      x = numeric(5),
      y = numeric(5)
    ),
    "pair_id column found in colData"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 6: .bootstrap_divergence_compute_samples()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_compute_samples returns numeric vector", {
  result <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    nboot = 500,
    q = 1,
    pseudocount = 0.5,
    log_base = exp(1),
    paired = FALSE,
    se = NULL,
    pair_id_col = NULL,
    group_col = "group",
    control_group = "Normal"
  )
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 500)
})

test_that(".bootstrap_divergence_compute_samples produces different values each call", {
  set.seed(123)
  result1 <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    nboot = 500,
    q = 1,
    pseudocount = 0.5,
    log_base = exp(1),
    paired = FALSE,
    se = NULL,
    pair_id_col = NULL,
    group_col = "group",
    control_group = "Normal"
  )
  
  set.seed(456)
  result2 <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    nboot = 500,
    q = 1,
    pseudocount = 0.5,
    log_base = exp(1),
    paired = FALSE,
    se = NULL,
    pair_id_col = NULL,
    group_col = "group",
    control_group = "Normal"
  )
  
  expect_false(identical(result1, result2))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 7: .bootstrap_divergence_handle_multiple_q()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_handle_multiple_q returns list with multiple q", {
  result <- TSENAT:::.bootstrap_divergence_handle_multiple_q(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = c(0.5, 1, 1.5),
    norm = FALSE,
    nboot = 500,
    ci = 0.95,
    method = "percentile",
    log_base = exp(1),
    pseudocount = 0.5,
    gene_name = "TEST_GENE",
    verbose = FALSE,
    paired = FALSE,
    pair_id_col = NULL,
    se = NULL
  )
  
  expect_equal(class(result), "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 3)
  expect_equal(names(result), c("q_0.5", "q_1", "q_1.5"))
})

test_that(".bootstrap_divergence_handle_multiple_q result contains valid CIs", {
  result <- TSENAT:::.bootstrap_divergence_handle_multiple_q(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = c(0.5, 1),
    norm = FALSE,
    nboot = 500,
    ci = 0.95,
    method = "percentile",
    log_base = exp(1),
    pseudocount = 0.5,
    gene_name = "TEST_GENE",
    verbose = FALSE,
    paired = FALSE,
    pair_id_col = NULL,
    se = NULL
  )
  
  for (i in 1:2) {
    expect_true(is.numeric(result[[i]]$estimate))
    expect_true(is.numeric(result[[i]]$lower_ci))
    expect_true(is.numeric(result[[i]]$upper_ci))
    expect_true(result[[i]]$lower_ci < result[[i]]$upper_ci)
  }
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 8: .bootstrap_divergence() - Main Function Integration Tests
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence returns valid single-q result", {
  result <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    ci = 0.95,
    method = "percentile",
    verbose = FALSE
  )
  
  expect_equal(class(result), "tsenat_divergence_bootstrap_ci")
  expect_true(is.numeric(result$estimate))
  expect_true(is.numeric(result$lower_ci))
  expect_true(is.numeric(result$upper_ci))
  expect_equal(result$ci_level, 0.95)
  expect_equal(result$q, 1)
  expect_equal(result$nboot, 500)
  expect_equal(result$method, "percentile")
})

test_that(".bootstrap_divergence returns list for multiple q values", {
  result <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = c(0.5, 1, 1.5),
    nboot = 500,
    ci = 0.95,
    method = "percentile",
    verbose = FALSE
  )
  
  expect_equal(class(result), "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 3)
})

test_that(".bootstrap_divergence uses BCA method correctly", {
  result_percentile <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    ci = 0.95,
    method = "percentile",
    verbose = FALSE
  )
  
  result_bca <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    ci = 0.95,
    method = "bca",
    verbose = FALSE
  )
  
  expect_equal(result_percentile$method, "percentile")
  expect_equal(result_bca$method, "bca")
  # Both should return valid CIs
  expect_true(result_percentile$lower_ci < result_percentile$upper_ci)
  expect_true(result_bca$lower_ci < result_bca$upper_ci)
})

test_that(".bootstrap_divergence accepts SummarizedExperiment input", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(100, 80, 60, 40, 90, 70, 50, 40),
      nrow = 2,
      ncol = 4,
      dimnames = list(c("GENE1", "GENE2"), NULL)
    )),
    colData = data.frame(group = c("Normal", "Normal", "Treatment", "Treatment"))
  )
  
  res <- data.frame(row.names = c("GENE1", "GENE2"), stat = c(5, 2))
  
  result <- TSENAT:::.bootstrap_divergence(
    se = se,
    res = res,
    top_n = 1,
    q = 1,
    nboot = 500,
    ci = 0.95,
    group_col = "group",
    control_group = "Normal",
    verbose = FALSE
  )
  
  expect_equal(class(result), "tsenat_divergence_bootstrap_ci")
  expect_equal(result$gene_name, "GENE1")
})

test_that(".bootstrap_divergence rejects invalid SE input", {
  expect_error(
    TSENAT:::.bootstrap_divergence(
      se = list(data = "not_se"),
      res = data.frame(),
      q = 1,
      nboot = 100
    ),
    "SummarizedExperiment"
  )
})

test_that(".bootstrap_divergence rejects missing x and y with no se", {
  expect_error(
    TSENAT:::.bootstrap_divergence(
      x = NULL,
      y = NULL,
      se = NULL,
      res = NULL,
      q = 1,
      nboot = 100
    ),
    "Must provide"
  )
})

test_that(".bootstrap_divergence handles different pseudocount values", {
  result_pc0 <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    pseudocount = 0,
    verbose = FALSE
  )
  
  result_pc05 <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    pseudocount = 0.5,
    verbose = FALSE
  )
  
  # Both should produce valid results
  expect_true(is.numeric(result_pc0$estimate))
  expect_true(is.numeric(result_pc05$estimate))
})

test_that(".bootstrap_divergence normalizes divergence when requested", {
  result_norm_false <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    norm = FALSE,
    verbose = FALSE
  )
  
  result_norm_true <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    norm = TRUE,
    verbose = FALSE
  )
  
  # Normalized result should be in [0, 1]
  expect_true(result_norm_true$estimate >= 0 && result_norm_true$estimate <= 1)
})

test_that(".bootstrap_divergence prints when verbose=TRUE and gene_name provided", {
  expect_no_error(
    TSENAT:::.bootstrap_divergence(
      x = c(100, 80, 60),
      y = c(90, 70, 50),
      q = 1,
      nboot = 500,
      gene_name = "TEST_GENE",
      verbose = TRUE
    )
  )
})

test_that(".bootstrap_divergence respects CI level parameter", {
  result_90 <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    ci = 0.90,
    verbose = FALSE
  )
  
  result_99 <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    ci = 0.99,
    verbose = FALSE
  )
  
  expect_equal(result_90$ci_level, 0.90)
  expect_equal(result_99$ci_level, 0.99)
  # 99% CI should be wider than 90% CI
  width_90 <- result_90$upper_ci - result_90$lower_ci
  width_99 <- result_99$upper_ci - result_99$lower_ci
  expect_true(width_99 >= width_90)
})

test_that(".bootstrap_divergence contains bootstrap distribution", {
  result <- TSENAT:::.bootstrap_divergence(
    x = c(100, 80, 60),
    y = c(90, 70, 50),
    q = 1,
    nboot = 500,
    verbose = FALSE
  )
  
  expect_true(is.numeric(result$bootstrap_dist))
  expect_true(length(result$bootstrap_dist) > 0)
  expect_true(length(result$bootstrap_dist) <= 500)
})

# ==============================================================================
# .bootstrap_divergence_handle_multiple_q(): Tests for multi-Q divergence (47.4%)
# ==============================================================================

test_that(".bootstrap_divergence_handle_multiple_q processes multiple q-values", {
  x <- c(10, 15, 20, 25, 30, 18, 22, 19, 21, 23)
  y <- c(12, 18, 22, 28, 32, 20, 24, 21, 23, 25)
  q_vals <- c(0.5, 1.0, 2.0)
  
  result <- .bootstrap_divergence_handle_multiple_q(
    x = x, y = y, q = q_vals, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, gene_name = "test_gene", verbose = FALSE,
    paired = FALSE, pair_id_col = NULL, se = NULL
  )
  
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), length(q_vals))
})

test_that(".bootstrap_divergence_handle_multiple_q returns list with results", {
  x <- c(10, 20, 30, 40, 50)
  y <- c(15, 25, 35, 45, 55)
  q_vals <- c(0.5, 1)
  
  result <- .bootstrap_divergence_handle_multiple_q(
    x = x, y = y, q = q_vals, norm = FALSE, nboot = 25,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 1e-6, gene_name = "test_gene2", verbose = FALSE,
    paired = FALSE, pair_id_col = NULL, se = NULL
  )
  
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 2)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST SUITE FOR .bootstrap_divergence_handle_multiple_q (47.3% coverage)
# ════════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_handle_multiple_q processes single q value", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result <- TSENAT:::.bootstrap_divergence_handle_multiple_q(
    x = x, y = y, q = 1.0, norm = FALSE, nboot = 5, ci = 0.95,
    method = "percentile", log_base = 10, pseudocount = 1, 
    gene_name = NULL, verbose = FALSE, paired = FALSE, pair_id_col = NULL, se = NULL
  )
  
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 1)
  expect_true("q_1" %in% names(result))
})

test_that(".bootstrap_divergence_handle_multiple_q processes multiple q values", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result <- TSENAT:::.bootstrap_divergence_handle_multiple_q(
    x = x, y = y, q = c(0.5, 1.0, 2.0), norm = FALSE, nboot = 5, ci = 0.95,
    method = "percentile", log_base = 10, pseudocount = 1, 
    gene_name = NULL, verbose = FALSE, paired = FALSE, pair_id_col = NULL, se = NULL
  )
  
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 3)
  expect_true(all(c("q_0.5", "q_1", "q_2") %in% names(result)))
})

test_that(".bootstrap_divergence_handle_multiple_q returns with gene_name and verbose", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  # When verbose=TRUE with gene_name, function should produce diagnostic output (messages)
  # We suppress warnings/messages to keep test output clean while verifying execution
  result <- suppressMessages({
    TSENAT:::.bootstrap_divergence_handle_multiple_q(
      x = x, y = y, q = c(1.0, 2.0), norm = FALSE, nboot = 5, ci = 0.95,
      method = "percentile", log_base = 10, pseudocount = 1, 
      gene_name = "TEST_GENE", verbose = TRUE, paired = FALSE, pair_id_col = NULL, se = NULL
    )
  })
  
  # Verify result structure is correct
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 2)
  expect_true(all(c("q_1", "q_2") %in% names(result)))
})

test_that(".bootstrap_divergence_handle_multiple_q with BCA method", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result <- TSENAT:::.bootstrap_divergence_handle_multiple_q(
    x = x, y = y, q = 1.0, norm = FALSE, nboot = 10, ci = 0.95,
    method = "bca", log_base = 10, pseudocount = 1, 
    gene_name = NULL, verbose = FALSE, paired = FALSE, pair_id_col = NULL, se = NULL
  )
  
  expect_is(result, "tsenat_divergence_bootstrap_list")
  expect_equal(length(result), 1)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST SUITE FOR .bootstrap_divergence_compute_samples (50% coverage)
# ════════════════════════════════════════════════════════════════════════════════

test_that(".bootstrap_divergence_compute_samples processes unpaired data", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 10, q = 1.0, pseudocount = 1,
    log_base = 10, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  # Should return numeric vector of bootstrap samples
  expect_is(result, "numeric")
  expect_equal(length(result), 10)
  expect_true(all(!is.na(result)))
})

test_that(".bootstrap_divergence_compute_samples handles different q values", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  # Test with q = 2.0 (different from q = 1.0 which is KL divergence)
  result <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 5, q = 2.0, pseudocount = 1,
    log_base = 10, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  # Should return numeric vector
  expect_is(result, "numeric")
  # 5 bootstrap replicates
  expect_equal(length(result), 5)
  expect_true(all(is.finite(result)))
})

test_that(".bootstrap_divergence_compute_samples with different log bases", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result_log10 <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 5, q = 1.0, pseudocount = 1,
    log_base = 10, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  result_log2 <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 5, q = 1.0, pseudocount = 1,
    log_base = 2, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  expect_is(result_log10, "numeric")
  expect_is(result_log2, "numeric")
  # Results should differ due to different log bases
  expect_false(isTRUE(all.equal(result_log10, result_log2)))
})

test_that(".bootstrap_divergence_compute_samples with pseudocount effect", {
  x <- c(100, 80, 60, 40, 20)
  y <- c(90, 70, 50, 40, 30)
  
  result_low_pc <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 5, q = 1.0, pseudocount = 0.1,
    log_base = 10, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  result_high_pc <- TSENAT:::.bootstrap_divergence_compute_samples(
    x = x, y = y, nboot = 5, q = 1.0, pseudocount = 10,
    log_base = 10, paired = FALSE, se = NULL, pair_id_col = NULL,
    group_col = NULL, control_group = NULL
  )
  
  expect_is(result_low_pc, "numeric")
  expect_is(result_high_pc, "numeric")
  expect_equal(length(result_low_pc), length(result_high_pc))
})
# M7: Paired bootstrap — warns when paired=TRUE but no pair column
# ============================================================================ 
# ============================================================================
# M7: Paired bootstrap — warns when paired=TRUE but no pair column# ============================================================================

test_that("M7: Paired bootstrap warns when no pair column detected", {
    skip_if_not_installed("SummarizedExperiment")
    set.seed(42)

    se <- create_count_se(n_genes = 4, n_samples = 6, n_control = 3, lambda = 100, seed = 42)
    colnames(se) <- paste0("Sample_", seq_len(ncol(se)))

    # paired=TRUE but no pairing column → should warn
    expect_warning(
        TSENAT:::.prepare_divergence_execution(
            se = se, bootstrap = TRUE, paired = TRUE,
            nboot = 10, method = "percentile", nthreads = 1, progress = FALSE
        ),
        "paired.*TRUE.*no pair"
    )
})
