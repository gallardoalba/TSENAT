context("Pseudocount Transparency: calculate_diversity_s4 and get_pseudocount")

# ===========================================================================
# Setup: Helper to create test analysis
# ===========================================================================

make_test_analysis <- function(n_genes = 10, n_samples = 4, seed = 123) {
  set.seed(seed)
  readcounts <- matrix(rpois(n_genes * n_samples, lambda = 10), 
                       nrow = n_genes, ncol = n_samples)
  colnames(readcounts) <- paste0("S", 1:n_samples)
  genes <- rep(paste0("G", 1:(n_genes/2)), each = 2)
  
  TSENAT:::build_analysis_s4(readcounts = readcounts, tx2gene = data.frame(
    txid = rownames(readcounts) <- paste0("T", 1:n_genes),
    geneid = genes
  ))
}

# ===========================================================================
# TEST GROUP 1: Basic Pseudocount Retrieval with Numeric Values
# ===========================================================================

test_that("get_pseudocount returns NULL when diversity not yet calculated", {
  analysis <- make_test_analysis()
  
  # Should return NULL before any diversity calculation
  pc <- TSENAT::get_pseudocount(analysis)
  expect_null(pc,
              info = "get_pseudocount should return NULL when diversity not calculated")
})

test_that("get_pseudocount returns numeric value when pseudocount is a number", {
  analysis <- make_test_analysis()
  
  # Calculate diversity with explicit numeric pseudocount
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0.5,
    verbose = FALSE
  )
  
  # Should retrieve the same value
  pc <- TSENAT::get_pseudocount(analysis)
  expect_is(pc, "numeric")
  expect_equal(pc, 0.5,
               info = "Should retrieve the exact pseudocount value passed in")
})

test_that("get_pseudocount(original=TRUE) returns numeric when input was numeric", {
  analysis <- make_test_analysis()
  
  # Calculate diversity with explicit numeric pseudocount
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0.5,
    verbose = FALSE
  )
  
  # Original should also be numeric
  pc_orig <- TSENAT::get_pseudocount(analysis, original = TRUE)
  expect_is(pc_orig, "numeric")
  expect_equal(pc_orig, 0.5)
})

# ===========================================================================
# TEST GROUP 2: Pseudocount "auto" Mode Resolution
# ===========================================================================

test_that('get_pseudocount returns numeric value when pseudocount="auto"', {
  analysis <- make_test_analysis()
  
  # Calculate diversity with auto pseudocount
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = "auto",
    verbose = FALSE
  )
  
  # Should retrieve a numeric value (resolved from "auto")
  pc <- TSENAT::get_pseudocount(analysis)
  expect_is(pc, "numeric",
            info = 'get_pseudocount should return numeric even when pseudocount="auto"')
  expect_gt(pc, 0,
            info = "Resolved pseudocount should be positive")
})

test_that('get_pseudocount(original=TRUE) returns "auto" when input was "auto"', {
  analysis <- make_test_analysis()
  
  # Calculate diversity with auto pseudocount
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = "auto",
    verbose = FALSE
  )
  
  # Original should be the "auto" string
  pc_orig <- TSENAT::get_pseudocount(analysis, original = TRUE)
  expect_is(pc_orig, "character")
  expect_equal(tolower(pc_orig), "auto",
               info = "original=TRUE should return the original 'auto' parameter")
})

test_that("auto pseudocount is different from zero pseudocount", {
  analysis1 <- make_test_analysis(seed = 123)
  analysis2 <- make_test_analysis(seed = 123)
  
  # One with auto, one with explicit 0
  analysis1 <- TSENAT::calculate_diversity_s4(
    analysis1,
    q = 1.0,
    pseudocount = "auto",
    verbose = FALSE
  )
  
  analysis2 <- TSENAT::calculate_diversity_s4(
    analysis2,
    q = 1.0,
    pseudocount = 0,
    verbose = FALSE
  )
  
  pc_auto <- TSENAT::get_pseudocount(analysis1)
  pc_zero <- TSENAT::get_pseudocount(analysis2)
  
  # Auto should be greater than zero
  expect_gt(pc_auto, pc_zero,
            info = "Auto-computed pseudocount should be > 0 for non-empty data")
})

# ===========================================================================
# TEST GROUP 3: Metadata Storage and Audit Trail
# ===========================================================================

test_that("pseudocount is stored in diversity_combined metadata", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0.75,
    verbose = FALSE
  )
  
  # Check metadata structure
  meta <- TSENAT::getMeta(analysis)
  expect_true("diversity_combined" %in% names(meta),
              info = "diversity_combined should be in metadata")
  
  div_combined <- meta$diversity_combined
  expect_true("computation_params" %in% names(div_combined),
              info = "computation_params should be in diversity_combined")
  
  params <- div_combined$computation_params
  expect_true("pseudocount" %in% names(params),
              info = "pseudocount should be stored in computation_params")
  expect_equal(params$pseudocount, 0.75)
})

test_that("both resolved and original pseudocount stored with auto", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = "auto",
    verbose = FALSE
  )
  
  # Check metadata
  meta <- TSENAT::getMeta(analysis)
  params <- meta$diversity_combined$computation_params
  
  # Should have both fields
  expect_true("pseudocount" %in% names(params),
              info = "pseudocount (resolved) should be stored")
  expect_true("pseudocount_original" %in% names(params),
              info = "pseudocount_original should be stored")
  
  # Resolved should be numeric
  expect_is(params$pseudocount, "numeric")
  
  # Original should be "auto" string
  expect_equal(tolower(params$pseudocount_original), "auto")
})

test_that("pseudocount stored in config last_diversity_run", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0.5,
    verbose = FALSE
  )
  
  # Check config
  config <- TSENAT::getConfig(analysis)
  expect_true("last_diversity_run" %in% names(config),
              info = "last_diversity_run should be in config")
  
  last_run <- config$last_diversity_run
  expect_true("parameters_used" %in% names(last_run),
              info = "parameters_used should be in last_diversity_run")
  
  params_used <- last_run$parameters_used
  expect_equal(params_used$pseudocount, 0.5,
               info = "pseudocount should be in parameters_used")
})

# ===========================================================================
# TEST GROUP 4: Multiple Q-values with Pseudocount
# ===========================================================================

test_that("pseudocount consistent across multiple q-values", {
  analysis <- make_test_analysis()
  
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_values,
    pseudocount = 0.25,
    verbose = FALSE
  )
  
  # Pseudocount should be the same regardless of q-values
  pc <- TSENAT::get_pseudocount(analysis)
  expect_equal(pc, 0.25,
               info = "Pseudocount should be consistent across q-values")
  
  # Verify all q-value results were computed
  div_results <- TSENAT::divResults(analysis)
  expect_equal(length(div_results), length(q_values),
               info = "Should have results for all q-values")
})

test_that("auto pseudocount consistent across multiple q-values", {
  analysis <- make_test_analysis()
  
  q_values <- c(0.5, 1.0, 1.5)
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_values,
    pseudocount = "auto",
    verbose = FALSE
  )
  
  # Auto-resolved pseudocount should be single value
  pc <- TSENAT::get_pseudocount(analysis)
  expect_is(pc, "numeric")
  expect_length(pc, 1,
                info = "Pseudocount should be single value even with multiple q")
  
  # Original should be "auto"
  pc_orig <- TSENAT::get_pseudocount(analysis, original = TRUE)
  expect_equal(tolower(pc_orig), "auto")
})

# ===========================================================================
# TEST GROUP 5: Verbose Output Messages
# ===========================================================================

test_that("verbose=TRUE shows pseudocount resolution message", {
  analysis <- make_test_analysis()
  
  # Capture messages when verbose=TRUE
  expect_message(
    {
      analysis <- TSENAT::calculate_diversity_s4(
        analysis,
        q = 1.0,
        pseudocount = "auto",
        verbose = TRUE
      )
    },
    "Resolved pseudocount",
    info = "Should show message about resolved pseudocount when verbose=TRUE"
  )
})

test_that("verbose=FALSE suppresses pseudocount resolution message", {
  analysis <- make_test_analysis()
  
  # Should not produce message when verbose=FALSE
  expect_silent(
    {
      analysis <- TSENAT::calculate_diversity_s4(
        analysis,
        q = 1.0,
        pseudocount = "auto",
        verbose = FALSE
      )
    },
    info = "Should not show message when verbose=FALSE"
  )
})

# ===========================================================================
# TEST GROUP 6: Integration with Accessor Functions
# ===========================================================================

test_that("get_pseudocount works with output from calculate_diversity_s4", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0, 2.0),
    pseudocount = "auto",
    verbose = FALSE
  )
  
  # Should work smoothly with diversity calculation output
  pc <- TSENAT::get_pseudocount(analysis)
  expect_is(pc, "numeric")
  
  # Verify diversity results also exist
  div_res <- TSENAT::divResults(analysis)
  expect_gt(length(div_res), 0,
            info = "Diversity results should be present")
})

test_that("pseudocount survives metadata access patterns", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 1.5,
    verbose = FALSE
  )
  
  # Get through different access patterns
  pc1 <- TSENAT::get_pseudocount(analysis)
  
  meta <- TSENAT::getMeta(analysis)
  pc2 <- meta$diversity_combined$computation_params$pseudocount
  
  # Should be identical
  expect_equal(pc1, pc2,
               info = "Pseudocount should be accessible through different metadata paths")
})

# ===========================================================================
# TEST GROUP 7: Edge Cases
# ===========================================================================

test_that("get_pseudocount handles empty/NULL analysis gracefully", {
  analysis <- make_test_analysis()
  
  # Without calling calculate_diversity_s4
  pc <- TSENAT::get_pseudocount(analysis)
  expect_null(pc,
              info = "Should return NULL for analysis without diversity calculation")
})

test_that("get_pseudocount parameter validation", {
  analysis <- make_test_analysis()
  
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0.5,
    verbose = FALSE
  )
  
  # Should accept boolean for original parameter
  pc1 <- TSENAT::get_pseudocount(analysis, original = FALSE)
  pc2 <- TSENAT::get_pseudocount(analysis, original = TRUE)
  
  expect_is(pc1, "numeric")
  expect_is(pc2, "numeric")
  expect_equal(pc1, pc2,
               info = "original parameter should not affect numeric input")
})

test_that("pseudocount zero is preserved correctly", {
  analysis <- make_test_analysis()
  
  # Explicitly set pseudocount to 0
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    pseudocount = 0,
    verbose = FALSE
  )
  
  pc <- TSENAT::get_pseudocount(analysis)
  expect_equal(pc, 0,
               info = "Pseudocount of 0 should be preserved exactly")
})
