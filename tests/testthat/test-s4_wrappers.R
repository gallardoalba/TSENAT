context("S4 Wrappers: Input Validation and Structure")

# Helper to create test analysis
make_test_se_for_wrappers <- function(n_genes = 20, n_samples = 8) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10), nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Sample", 1:n_samples)
  SummarizedExperiment(assays = list(counts = counts))
}

# ============================================================================
# DIVERSITY WRAPPER TESTS
# ============================================================================

test_that("calculate_diversity_s4 validates input is TSENATAnalysis", {
  expect_error(calculate_diversity_s4("not_analysis", q = 1.0), "TSENATAnalysis")
})

test_that("calculate_diversity_s4 validates q is numeric", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(calculate_diversity_s4(analysis, q = "not_numeric"), "must be numeric")
})

test_that("calculate_diversity_s4 rejects empty SummarizedExperiment", {
  empty_se <- SummarizedExperiment(assays = list(counts = matrix(0, 0, 0)))
  analysis <- TSENATAnalysis(empty_se)
  expect_error(calculate_diversity_s4(analysis), "empty")
})

# ============================================================================
# LM INTERACTION WRAPPER TESTS
# ============================================================================

test_that("calculate_lm_interaction_s4 validates input", {
  expect_error(calculate_lm_interaction_s4("not_analysis"), "TSENATAnalysis")
})

test_that("calculate_lm_interaction_s4 accepts formula from config", {
  se <- make_test_se_for_wrappers()
  cfg <- list(formula = ~ treatment)
  analysis <- TSENATAnalysis(se, config = cfg)
  expect_equal(analysis@config$formula, ~ treatment)
})

# ============================================================================
# JACKKNIFE WRAPPER TESTS
# ============================================================================

test_that("jackknife_tsallis_entropy_s4 requires diversity results", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    jackknife_tsallis_entropy_s4(analysis, q = 1.0),
    "Diversity results required"
  )
})

test_that("jackknife_tsallis_entropy_s4 errors on unavailable q-value", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)

  # Manually add diversity for q=1.0
  analysis@diversity_results$q_1.0 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )

  # Jackknife for q=2.0 should error
  expect_error(
    jackknife_tsallis_entropy_s4(analysis, q = 2.0),
    "not calculated"
  )
})

# ============================================================================
# DIVERGENCE WRAPPER TESTS
# ============================================================================

test_that("calculate_divergence_s4 requires diversity results", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    calculate_divergence_s4(analysis),
    "Diversity results required"
  )
})

# ============================================================================
# Q-GENE INTERACTIONS WRAPPER TESTS
# ============================================================================

test_that("detect_q_gene_interactions_s4 requires diversity", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    detect_q_gene_interactions_s4(analysis),
    "Diversity results required"
  )
})

# ============================================================================
# DIFFERENCE WRAPPER TESTS
# ============================================================================

test_that("calculate_difference_s4 requires control specification", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    calculate_difference_s4(analysis),
    "must be specified"
  )
})

# ============================================================================
# METADATA TRACKING TESTS
# ============================================================================

test_that("Wrappers initialize metadata tracking structure", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)

  expect_true("function_calls" %in% names(analysis@metadata))
  expect_true(is.character(analysis@metadata$function_calls))
  expect_equal(length(analysis@metadata$function_calls), 0)
})
