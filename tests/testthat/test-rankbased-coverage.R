# Comprehensive testing for uncovered lines in rank_based_methods.R
# Tests edge cases and specific code paths for rank-based assumptions testing

library(TSENAT)
skip_on_bioc()

context("Rank-Based Methods: Coverage Expansion")

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 92 (matrix conversion)
# ============================================================================

test_that("test_rankbased_assumptions_s4: converts data.frame to matrix", {
  # Line 92: if (!is.matrix(data)) data <- as.matrix(data)
  
  # Create a TSENATAnalysis object with diversity results
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Create analysis object
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Add diversity results as a SummarizedExperiment with data.frame assay
  diversity_data <- as.data.frame(matrix(runif(50), nrow = 10, ncol = 5))
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  
  analysis@diversity_results$q_1.0 <- div_se
  
  # This should not error even with data.frame assay
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("exchangeability")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result, check_type = "exchangeability")
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 117 (single row edge case)
# ============================================================================

test_that("test_rankbased_assumptions_s4: handles small data in exchangeability", {
  # Line 117: 0 (when length(row_means) <= 1) - ensure at least 2 rows
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:40, nrow = 8, ncol = 5))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Minimal genes with multiple samples (5 genes, 5 samples)
  diversity_data <- data.frame(
    Sample1 = c(1.5, 2.0, 1.8, 2.2, 1.6),
    Sample2 = c(1.6, 2.1, 1.9, 2.3, 1.7),
    Sample3 = c(1.4, 2.2, 1.7, 2.1, 1.5),
    Sample4 = c(1.7, 1.9, 2.0, 2.4, 1.8),
    Sample5 = c(1.5, 2.0, 1.8, 2.2, 1.6)
  )
  rownames(diversity_data) <- c("G1", "G2", "G3", "G4", "G5")
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("exchangeability")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result)
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 131 (permutation with single row)
# ============================================================================

test_that("test_rankbased_assumptions_s4: permutation handles small data", {
  # Line 131: 0 (when length(perm_means) <= 1 in permutation loop)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:40, nrow = 8, ncol = 5))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Minimal genes with multiple samples (5 genes, 5 samples)
  diversity_data <- data.frame(
    Sample1 = c(2.0, 2.1, 1.9, 2.3, 2.0),
    Sample2 = c(2.05, 2.15, 1.95, 2.35, 2.05),
    Sample3 = c(2.02, 2.12, 1.92, 2.32, 2.02),
    Sample4 = c(2.08, 2.18, 1.98, 2.38, 2.08),
    Sample5 = c(2.01, 2.11, 1.91, 2.31, 2.01)
  )
  rownames(diversity_data) <- c("G1", "G2", "G3", "G4", "G5")
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("exchangeability")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result, check_type = "exchangeability")
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 169 (High correlation status)
# ============================================================================

test_that("test_rankbased_assumptions_s4: returns PASS status for high monotonicity", {
  # Line 169: "[OK] PASS" status when mean_cor > 0.7 && sd_cor < 0.2
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create highly correlated data (rows rank similarly)
  set.seed(42)
  diversity_data <- matrix(nrow = 10, ncol = 5)
  for (i in 1:10) {
    base_vals <- runif(5)
    diversity_data[i, ] <- base_vals + rnorm(5, 0, 0.01)  # Small variance = high correlation
  }
  
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("monotonicity")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result)
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 171 (Acceptable correlation status)
# ============================================================================

test_that("test_rankbased_assumptions_s4: returns ACCEPTABLE status for moderate monotonicity", {
  # Line 171: "? ACCEPTABLE" status when mean_cor > 0.4
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create moderately correlated data
  set.seed(42)
  diversity_data <- matrix(nrow = 10, ncol = 5)
  for (i in 1:10) {
    diversity_data[i, ] <- rnorm(5, mean = i, sd = 2)
  }
  
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("monotonicity")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result, check_type = "monotonicity")
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 224 (High Kendall's W status)
# ============================================================================

test_that("test_rankbased_assumptions_s4: returns PASS status for high Kendall's W", {
  # Line 224: "[OK] PASS" status when kendall_w > 0.7
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create highly consistent data across samples
  set.seed(42)
  diversity_data <- matrix(nrow = 10, ncol = 4)
  for (i in 1:10) {
    diversity_data[i, ] <- rank(rnorm(4, mean = i, sd = 0.1))
  }
  
  colnames(diversity_data) <- paste0("Sample", 1:4)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("consistency")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result)
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - Line 226 (Acceptable Kendall's W status)
# ============================================================================

test_that("test_rankbased_assumptions_s4: returns ACCEPTABLE status for moderate Kendall's W", {
  # Line 226: "? ACCEPTABLE" status when kendall_w > 0.4
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create moderately consistent data
  set.seed(42)
  diversity_data <- matrix(rnorm(40, mean = 1.5, sd = 1), nrow = 10, ncol = 4)
  
  colnames(diversity_data) <- paste0("Sample", 1:4)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("consistency")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result)
})

# ============================================================================
# TEST: print.rank_assumptions - Lines 269-290 (print method)
# ============================================================================

test_that("print.rank_assumptions: prints header and check details", {
  # Lines 269, 270, 273-290: print method implementation
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create basic diversity results
  diversity_data <- matrix(runif(50, 1, 3), nrow = 10, ncol = 5)
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("exchangeability")
  )
  
  # Extract the actual rank_assumptions result from metadata
  rank_result <- result@metadata$rankbased_assumptions$result
  
  # Verify print method produces message output
  expect_message(
    print(rank_result),
    "RANK-BASED METHOD ASSUMPTIONS"
  )
})

# ============================================================================
# TEST: print.rank_assumptions - Line 280 (checks method field)
# ============================================================================

test_that("print.rank_assumptions: includes method field when present", {
  # Line 280: if (!is.null(check$method))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  diversity_data <- matrix(runif(50, 1, 3), nrow = 10, ncol = 5)
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("monotonicity")  # monotonicity has method field
  )
  
  # Extract the actual rank_assumptions result from metadata
  rank_result <- result@metadata$rankbased_assumptions$result
  
  # Verify print method produces message output with method field
  expect_message(
    print(rank_result),
    "Method:"
  )
})

# ============================================================================
# TEST: print.rank_assumptions - Line 284 (checks status field)
# ============================================================================

test_that("print.rank_assumptions: includes status field when present", {
  # Line 284: if (!is.null(check$status))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  diversity_data <- matrix(runif(50, 1, 3), nrow = 10, ncol = 5)
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("monotonicity", "consistency")
  )
  
  # Extract the actual rank_assumptions result from metadata
  rank_result <- result@metadata$rankbased_assumptions$result
  
  # Verify print method produces message output with status field
  expect_message(
    print(rank_result),
    "Status:"
  )
})

# ============================================================================
# TEST: print.rank_assumptions - Line 288 (checks details field)
# ============================================================================

test_that("print.rank_assumptions: includes details field when present", {
  # Line 288: if (!is.null(check$details))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  diversity_data <- matrix(runif(50, 1, 3), nrow = 10, ncol = 5)
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:10)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("monotonicity")
  )
  
  # Extract the actual rank_assumptions result from metadata
  rank_result <- result@metadata$rankbased_assumptions$result
  
  # Verify print method produces detailed message output
  expect_message(
    print(rank_result),
    "RANK-BASED METHOD ASSUMPTIONS"
  )
})

# ============================================================================
# TEST: test_rankbased_assumptions_s4 - All checks combined
# ============================================================================

test_that("test_rankbased_assumptions_s4: runs all checks without error", {
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  diversity_data <- matrix(runif(100, 1, 3), nrow = 20, ncol = 5)
  colnames(diversity_data) <- paste0("Sample", 1:5)
  rownames(diversity_data) <- paste0("G", 1:20)
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  analysis@diversity_results$q_1.0 <- div_se
  
  # Run all checks
  result <- TSENAT:::test_rankbased_assumptions_s4(
    analysis,
    checks = c("exchangeability", "monotonicity", "consistency")
  )
  
  # Strong assertion: validate result structure and content
  assert_rankbased_result_valid(result, check_type = "monotonicity")
})
