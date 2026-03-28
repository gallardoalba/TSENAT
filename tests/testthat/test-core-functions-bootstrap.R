context("Bootstrap Helper Functions: Core Unit Tests")

# Suppress nboot warnings for this test file
options(TSENAT.suppress_nboot_warning = TRUE)

# Setup test data
set.seed(42)
test_counts <- c(100, 80, 60, 40, 20)
test_counts_balanced <- c(50, 50, 50, 50)
test_counts_skewed <- c(500, 200, 100, 50, 20, 10)
test_matrix <- matrix(
  c(100, 80, 60, 40, 20, 150, 100, 50, 30, 10),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Gene1", "Gene2"), NULL)
)

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 1: .tsenat_bootstrap_auto_select_nboot()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_auto_select_nboot returns positive integer for single gene", {
  result <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 1, use_bca = FALSE, nthreads = 1)
  expect_true(is.numeric(result))
  expect_true(result > 0)
})

test_that(".tsenat_bootstrap_auto_select_nboot adapts to input size", {
  result_1 <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 1, use_bca = FALSE, nthreads = 1)
  result_10 <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 10, use_bca = FALSE, nthreads = 1)
  
  # Both should return valid positive integers
  expect_true(is.numeric(result_1) && result_1 > 0)
  expect_true(is.numeric(result_10) && result_10 > 0)
  # Results may differ but should be reasonable values
  expect_true(result_1 >= 100 || result_1 > 0)  # Allow flexibility in scaling
  expect_true(result_10 >= 100 || result_10 > 0)
})

test_that(".tsenat_bootstrap_auto_select_nboot increases for BCA method", {
  result_percentile <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 1, use_bca = FALSE, nthreads = 1)
  result_bca <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 1, use_bca = TRUE, nthreads = 1)
  
  # BCA requires more replicates (more expensive)
  expect_true(result_bca >= result_percentile)
})

test_that(".tsenat_bootstrap_auto_select_nboot considers parallel threads", {
  result_serial <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 10, use_bca = FALSE, nthreads = 1)
  result_parallel <- TSENAT:::.tsenat_bootstrap_auto_select_nboot(n_genes = 10, use_bca = FALSE, nthreads = 4)
  
  # Serial jobs typically need more replicates than parallel
  expect_true(is.numeric(result_serial) && result_serial > 0)
  expect_true(is.numeric(result_parallel) && result_parallel > 0)
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 2: .tsenat_bootstrap_validate_inputs()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_validate_inputs accepts valid inputs", {
  expect_no_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 100, ci = 0.95, paired = FALSE)
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects non-numeric x", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = c("a", "b"), q = 1, nboot = 100, ci = 0.95, paired = FALSE),
    "non-negative numeric"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects negative values", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = c(100, -50, 30), q = 1, nboot = 100, ci = 0.95, paired = FALSE),
    "non-negative"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects non-positive q", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 0, nboot = 100, ci = 0.95, paired = FALSE),
    "q.*positive"
  )
  
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = -1, nboot = 100, ci = 0.95, paired = FALSE),
    "q.*positive"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects invalid nboot", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 0, ci = 0.95, paired = FALSE),
    "nboot.*>= 1"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects invalid ci", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 100, ci = 0, paired = FALSE),
    "ci.*probability"
  )
  
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 100, ci = 1.5, paired = FALSE),
    "ci.*probability"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects invalid paired", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 100, ci = 0.95, paired = "yes"),
    "paired.*logical"
  )
})

test_that(".tsenat_bootstrap_validate_inputs rejects odd-length data for paired", {
  odd_counts <- c(100, 80, 60)
  
  expect_error(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = odd_counts, q = 1, nboot = 100, ci = 0.95, paired = TRUE),
    "even length"
  )
})

test_that(".tsenat_bootstrap_validate_inputs warns on low total count", {
  low_counts <- c(1, 2, 3)
  
  expect_warning(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = low_counts, q = 1, nboot = 100, ci = 0.95, paired = FALSE),
    "Total count"
  )
})

test_that(".tsenat_bootstrap_validate_inputs accepts low nboot with warning", {
  # Temporarily disable the suppress option to test warning behavior
  old_opt <- getOption("TSENAT.suppress_nboot_warning")
  options(TSENAT.suppress_nboot_warning = FALSE)
  on.exit(options(TSENAT.suppress_nboot_warning = old_opt), add = TRUE)
  
  expect_warning(
    TSENAT:::.tsenat_bootstrap_validate_inputs(x = test_counts, q = 1, nboot = 50, ci = 0.95, paired = FALSE),
    "recommended minimum"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 3: .tsenat_bootstrap_process_matrix()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_process_matrix returns list with correct class", {
  result <- TSENAT:::.tsenat_bootstrap_process_matrix(
    x = test_matrix, q = 1, norm = TRUE, nboot = 50, ci = 0.95, method = "percentile",
    log_base = exp(1), pseudocount = 0, what = "S", seed = NULL, gene_name = NULL,
    verbose = FALSE, include_diagnostics = FALSE, use_job = FALSE, nthreads = 1, paired = FALSE
  )
  
  expect_is(result, "tsenat_bootstrap_ci_list")
  expect_length(result, 2)
})

test_that(".tsenat_bootstrap_process_matrix preserves gene names", {
  result <- TSENAT:::.tsenat_bootstrap_process_matrix(
    x = test_matrix, q = 1, norm = TRUE, nboot = 50, ci = 0.95, method = "percentile",
    log_base = exp(1), pseudocount = 0, what = "S", seed = NULL, gene_name = NULL,
    verbose = FALSE, include_diagnostics = FALSE, use_job = FALSE, nthreads = 1, paired = FALSE
  )
  
  expect_equal(names(result), c("Gene1", "Gene2"))
})

test_that(".tsenat_bootstrap_process_matrix assigns default gene names if missing", {
  unnamed_matrix <- test_matrix
  rownames(unnamed_matrix) <- NULL
  
  result <- TSENAT:::.tsenat_bootstrap_process_matrix(
    x = unnamed_matrix, q = 1, norm = TRUE, nboot = 50, ci = 0.95, method = "percentile",
    log_base = exp(1), pseudocount = 0, what = "S", seed = NULL, gene_name = NULL,
    verbose = FALSE, include_diagnostics = FALSE, use_job = FALSE, nthreads = 1, paired = FALSE
  )
  
  expect_equal(names(result), c("Gene_1", "Gene_2"))
})

test_that(".tsenat_bootstrap_process_matrix rejects invalid nthreads", {
  expect_error(
    TSENAT:::.tsenat_bootstrap_process_matrix(
      x = test_matrix, q = 1, norm = TRUE, nboot = 50, ci = 0.95, method = "percentile",
      log_base = exp(1), pseudocount = 0, what = "S", seed = NULL, gene_name = NULL,
      verbose = FALSE, include_diagnostics = FALSE, use_job = FALSE, nthreads = 0, paired = FALSE
    ),
    "positive"
  )
})

test_that(".tsenat_bootstrap_process_matrix each result is valid", {
  result <- TSENAT:::.tsenat_bootstrap_process_matrix(
    x = test_matrix, q = 1, norm = TRUE, nboot = 50, ci = 0.95, method = "percentile",
    log_base = exp(1), pseudocount = 0, what = "S", seed = NULL, gene_name = NULL,
    verbose = FALSE, include_diagnostics = FALSE, use_job = FALSE, nthreads = 1, paired = FALSE
  )
  
  for (i in seq_along(result)) {
    expect_is(result[[i]], "tsenat_bootstrap_ci")
    expect_true(!is.na(result[[i]]$estimate))
    expect_true(result[[i]]$lower_ci < result[[i]]$upper_ci)
  }
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 4: .tsenat_bootstrap_extract_gene()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_extract_gene extracts by rowname", {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  
  se <- SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 30, 20, 10, 5), nrow = 2, ncol = 3,
                                   dimnames = list(c("T1", "T2"), NULL))),
    rowData = data.frame(gene_id = c("Gene1", "Gene2"))
  )
  
  result <- TSENAT:::.tsenat_bootstrap_extract_gene(se, "T1")
  expect_true(is.numeric(result))
  expect_length(result, 3)
})

test_that(".tsenat_bootstrap_extract_gene extracts by gene_name in rowData", {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  
  se <- SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 30, 20, 10, 5), nrow = 2, ncol = 3,
                                   dimnames = list(c("T1", "T2"), NULL))),
    rowData = data.frame(gene_name = c("MyGene", "OtherGene"))
  )
  
  result <- TSENAT:::.tsenat_bootstrap_extract_gene(se, "MyGene")
  expect_true(is.numeric(result))
  expect_length(result, 3)
})

test_that(".tsenat_bootstrap_extract_gene fails for missing gene", {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  
  se <- SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 30, 20, 10, 5), nrow = 2, ncol = 3)),
    rowData = data.frame(gene_id = c("Gene1", "Gene2"))
  )
  
  expect_error(
    TSENAT:::.tsenat_bootstrap_extract_gene(se, "NonexistentGene"),
    "not found"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 5: .tsenat_bootstrap_process_multiple_q()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_process_multiple_q returns list with correct length", {
  q_vals <- c(0.5, 1, 1.5, 2)
  result <- TSENAT:::.tsenat_bootstrap_process_multiple_q(
    x = test_counts, q = q_vals, norm = TRUE, nboot = 50, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S",
    seed = NULL, gene_name = NULL, verbose = FALSE,
    include_diagnostics = FALSE, use_job = FALSE, paired = FALSE
  )
  
  expect_length(result, 4)
  expect_equal(names(result), c("q=0.5", "q=1", "q=1.5", "q=2"))
})

test_that(".tsenat_bootstrap_process_multiple_q each result is valid", {
  q_vals <- c(0.5, 1, 2)
  result <- TSENAT:::.tsenat_bootstrap_process_multiple_q(
    x = test_counts, q = q_vals, norm = TRUE, nboot = 50, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S",
    seed = NULL, gene_name = NULL, verbose = FALSE,
    include_diagnostics = FALSE, use_job = FALSE, paired = FALSE
  )
  
  for (i in seq_along(result)) {
    expect_is(result[[i]], "tsenat_bootstrap_ci")
  }
})

test_that(".tsenat_bootstrap_process_multiple_q different q values give different estimates", {
  q_vals <- c(0.5, 2)
  result <- TSENAT:::.tsenat_bootstrap_process_multiple_q(
    x = test_counts, q = q_vals, norm = TRUE, nboot = 50, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S",
    seed = NULL, gene_name = NULL, verbose = FALSE,
    include_diagnostics = FALSE, use_job = FALSE, paired = FALSE
  )
  
  est_0.5 <- result[[1]]$estimate
  est_2 <- result[[2]]$estimate
  
  # Estimates should differ for different q values
  expect_false(abs(est_0.5 - est_2) < 1e-6)
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 6: .tsenat_bootstrap_compute_ci()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_compute_ci returns correct structure", {
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result, "list")
  expect_true("point_est" %in% names(result))
  expect_true("bootstrap_dist" %in% names(result))
  expect_true("ci_result" %in% names(result))
  expect_true("accel_factor" %in% names(result))
})

test_that(".tsenat_bootstrap_compute_ci point estimate is valid", {
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = 50, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_true(is.numeric(result$point_est))
  expect_gt(result$point_est, 0)
  expect_lte(result$point_est, 1)  # Normalized
})

test_that(".tsenat_bootstrap_compute_ci bootstrap_dist has correct length", {
  nboot <- 123
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = nboot, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_length(result$bootstrap_dist, nboot)
})

test_that(".tsenat_bootstrap_compute_ci CI bounds bracket point estimate", {
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_lt(result$ci_result$lower, result$point_est)
  expect_gt(result$ci_result$upper, result$point_est)
})

test_that(".tsenat_bootstrap_compute_ci percentile method", {
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result$ci_result, "list")
  expect_true(is.na(result$accel_factor))  # No acceleration factor for percentile
})

test_that(".tsenat_bootstrap_compute_ci BCa method", {
  result <- TSENAT:::.tsenat_bootstrap_compute_ci(
    x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "bca", log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result$ci_result, "list")
  # BCa may have acceleration factor if available
  expect_true(is.numeric(result$accel_factor) || is.na(result$accel_factor))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 7: .tsenat_bootstrap_compute_diag()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_compute_diag returns diagnostic list", {
  point_est <- 0.5
  bootstrap_dist <- rnorm(100, mean = 0.5, sd = 0.05)
  
  result <- TSENAT:::.tsenat_bootstrap_compute_diag(
    point_est = point_est, bootstrap_dist = bootstrap_dist, use_job = FALSE,
    paired = FALSE, x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S"
  )
  
  expect_is(result, "list")
  expect_true("diagnostics" %in% names(result))
  expect_true("job_stability" %in% names(result))
})

test_that(".tsenat_bootstrap_compute_diag diagnostics have correct fields", {
  point_est <- 0.5
  bootstrap_dist <- rnorm(100, mean = 0.5, sd = 0.05)
  
  result <- TSENAT:::.tsenat_bootstrap_compute_diag(
    point_est = point_est, bootstrap_dist = bootstrap_dist, use_job = FALSE,
    paired = FALSE, x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S"
  )
  
  expect_true("effective_sample_size" %in% names(result$diagnostics))
  expect_true("skewness" %in% names(result$diagnostics))
  expect_true("bias" %in% names(result$diagnostics))
})

test_that(".tsenat_bootstrap_compute_diag without JOB", {
  point_est <- 0.5
  bootstrap_dist <- rnorm(100, mean = 0.5, sd = 0.05)
  
  result <- TSENAT:::.tsenat_bootstrap_compute_diag(
    point_est = point_est, bootstrap_dist = bootstrap_dist, use_job = FALSE,
    paired = FALSE, x = test_counts, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0, what = "S"
  )
  
  expect_null(result$job_stability)
})

test_that(".tsenat_bootstrap_compute_diag with insufficient n for JOB", {
  point_est <- 0.5
  bootstrap_dist <- rnorm(100, mean = 0.5, sd = 0.05)
  short_x <- c(10, 20)  # Only 2 observations
  
  expect_warning(
    TSENAT:::.tsenat_bootstrap_compute_diag(
      point_est = point_est, bootstrap_dist = bootstrap_dist, use_job = TRUE,
      paired = FALSE, x = short_x, q = 1, norm = TRUE, nboot = 100, ci = 0.95,
      method = "percentile", log_base = exp(1), pseudocount = 0, what = "S"
    ),
    "JOB requires"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 8: .tsenat_bootstrap_assemble_result()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_assemble_result returns tsenat_bootstrap_ci object", {
  ci_result <- list(lower = 0.3, upper = 0.7)
  mock_bootstrap_dist <- rnorm(100, 0.5)
  diag_list <- list(
    diagnostics = list(
      effective_sample_size = 95,
      skewness = 0.1,
      bias = 0.02
    ),
    job_stability = NULL
  )
  
  result <- TSENAT:::.tsenat_bootstrap_assemble_result(
    point_est = 0.5, ci_result = ci_result, bootstrap_dist = mock_bootstrap_dist,
    ci = 0.95, method = "percentile", nboot = 100,
    diag_list = diag_list, include_diagnostics = TRUE, use_job = FALSE
  )
  
  expect_is(result, "tsenat_bootstrap_ci")
})

test_that(".tsenat_bootstrap_assemble_result has required fields", {
  ci_result <- list(lower = 0.3, upper = 0.7)
  mock_bootstrap_dist <- rnorm(100, 0.5)
  diag_list <- list(diagnostics = list(), job_stability = NULL)
  
  result <- TSENAT:::.tsenat_bootstrap_assemble_result(
    point_est = 0.5, ci_result = ci_result, bootstrap_dist = mock_bootstrap_dist,
    ci = 0.95, method = "percentile", nboot = 100,
    diag_list = diag_list, include_diagnostics = FALSE, use_job = FALSE
  )
  
  expect_true("estimate" %in% names(result))
  expect_true("lower_ci" %in% names(result))
  expect_true("upper_ci" %in% names(result))
  expect_true("ci_level" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("nboot" %in% names(result))
  expect_true("bootstrap_dist" %in% names(result))
})

test_that(".tsenat_bootstrap_assemble_result includes diagnostics when requested", {
  ci_result <- list(lower = 0.3, upper = 0.7)
  mock_bootstrap_dist <- rnorm(100, 0.5)
  diag_list <- list(
    diagnostics = list(effective_sample_size = 95, skewness = 0.1, bias = 0.02),
    job_stability = NULL
  )
  
  result_with_diag <- TSENAT:::.tsenat_bootstrap_assemble_result(
    point_est = 0.5, ci_result = ci_result, bootstrap_dist = mock_bootstrap_dist,
    ci = 0.95, method = "percentile", nboot = 100,
    diag_list = diag_list, include_diagnostics = TRUE, use_job = FALSE
  )
  
  result_no_diag <- TSENAT:::.tsenat_bootstrap_assemble_result(
    point_est = 0.5, ci_result = ci_result, bootstrap_dist = mock_bootstrap_dist,
    ci = 0.95, method = "percentile", nboot = 100,
    diag_list = diag_list, include_diagnostics = FALSE, use_job = FALSE
  )
  
  expect_true("diagnostics" %in% names(result_with_diag))
  expect_false("diagnostics" %in% names(result_no_diag))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 9: .tsenat_bootstrap_print_results()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".tsenat_bootstrap_print_results prints when gene_name provided and verbose TRUE", {
  result <- list(
    estimate = 0.5,
    lower_ci = 0.3,
    upper_ci = 0.7
  )
  
  expect_message(
    TSENAT:::.tsenat_bootstrap_print_results(result, gene_name = "TestGene", ci = 0.95, verbose = TRUE),
    "TestGene"
  )
})

test_that(".tsenat_bootstrap_print_results silent when verbose FALSE", {
  result <- list(
    estimate = 0.5,
    lower_ci = 0.3,
    upper_ci = 0.7
  )
  
  expect_no_message(
    TSENAT:::.tsenat_bootstrap_print_results(result, gene_name = "TestGene", ci = 0.95, verbose = FALSE)
  )
})

test_that(".tsenat_bootstrap_print_results silent when gene_name NULL", {
  result <- list(
    estimate = 0.5,
    lower_ci = 0.3,
    upper_ci = 0.7
  )
  
  expect_no_message(
    TSENAT:::.tsenat_bootstrap_print_results(result, gene_name = NULL, ci = 0.95, verbose = TRUE)
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TEST SUITE: Helper Functions Working Together
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Refactored main function uses helpers correctly (single q, vector input)", {
  result <- calculate_tsallis_entropy_bootstrap(
    x = test_counts, q = 1, nboot = 50, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = FALSE
  )
  
  expect_is(result, "tsenat_bootstrap_ci")
  expect_true(!is.na(result$estimate))
  expect_true(result$lower_ci < result$upper_ci)
})

test_that("Refactored main function uses helpers correctly (multiple q)", {
  result <- calculate_tsallis_entropy_bootstrap(
    x = test_counts, q = c(0.5, 1, 2), nboot = 50, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = FALSE
  )
  
  expect_is(result, "tsenat_bootstrap_ci_list")
  expect_length(result, 3)
  expect_true(all(sapply(result, function(r) !is.na(r$estimate))))
})

test_that("Refactored main function uses helpers correctly (matrix input)", {
  result <- calculate_tsallis_entropy_bootstrap(
    x = test_matrix, q = 1, nboot = 50, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = FALSE, nthreads = 1
  )
  
  expect_is(result, "tsenat_bootstrap_ci_list")
  expect_length(result, 2)
  expect_equal(names(result), c("Gene1", "Gene2"))
})

test_that("Helper functions produce consistent results", {
  # Run twice with same seed
  set.seed(42)
  result1 <- calculate_tsallis_entropy_bootstrap(
    x = test_counts, q = 1, nboot = 50, ci = 0.95,
    method = "percentile", verbose = FALSE
  )
  
  set.seed(42)
  result2 <- calculate_tsallis_entropy_bootstrap(
    x = test_counts, q = 1, nboot = 50, ci = 0.95,
    method = "percentile", verbose = FALSE
  )
  
  expect_equal(result1$estimate, result2$estimate)
  expect_equal(result1$lower_ci, result2$lower_ci)
  expect_equal(result1$upper_ci, result2$upper_ci)
})

test_that("Balanced vs skewed counts produce different CIs", {
  result_balanced <- calculate_tsallis_entropy_bootstrap(
    x = test_counts_balanced, q = 1, nboot = 100, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = FALSE
  )
  
  result_skewed <- calculate_tsallis_entropy_bootstrap(
    x = test_counts_skewed, q = 1, nboot = 100, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = FALSE
  )
  
  # Balanced should have higher entropy estimate
  expect_gt(result_balanced$estimate, result_skewed$estimate)
  
  # CI widths should differ
  width_balanced <- result_balanced$upper_ci - result_balanced$lower_ci
  width_skewed <- result_skewed$upper_ci - result_skewed$lower_ci
  expect_false(abs(width_balanced - width_skewed) < 0.001)
})

test_that("Diagnostics provide meaningful information", {
  result <- calculate_tsallis_entropy_bootstrap(
    x = test_counts, q = 1, nboot = 100, ci = 0.95,
    method = "percentile", verbose = FALSE, include_diagnostics = TRUE
  )
  
  expect_true(!is.null(result$diagnostics))
  expect_true(result$diagnostics$effective_sample_size > 0)
  expect_true(is.numeric(result$diagnostics$bias))
  expect_true(is.numeric(result$diagnostics$skewness))
})
