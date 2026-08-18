context("TSENATAnalysis S4 Class: Subsetting and Validation")

# ===========================================================================
# Helper function: Create standardized test data
# ===========================================================================
.create_test_analysis <- function(n_genes = 20, n_samples = 10, seed = 42) {
  set.seed(seed)
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 50, size = 2),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("GENE_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("A", "B"), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = colData
  )
  
  # Create TSENATAnalysis object
  analysis <- new("TSENATAnalysis", se = se)
  
  # Add some results to test subsetting of different result types
  analysis@config <- list(q_values = c(0.5, 1.0), group_col = "condition")
  
  # Add diversity results (list of matrices)
  analysis@diversity_results$q_0.5 <- matrix(rnorm(n_genes * n_samples), 
                                              nrow = n_genes, ncol = n_samples,
                                              dimnames = list(rownames(counts), colnames(counts)))
  analysis@diversity_results$q_1.0 <- matrix(rnorm(n_genes * n_samples), 
                                              nrow = n_genes, ncol = n_samples,
                                              dimnames = list(rownames(counts), colnames(counts)))
  
  # Add jackknife results with resamples and influence scores
  analysis@jackknife_results$q_0.5 <- list(
    resamples = matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples),
    influence_scores = matrix(rnorm(n_samples * 5), nrow = n_samples, ncol = 5),
    ci_matrix = array(rnorm(n_genes * n_samples * 2), dim = c(n_genes, n_samples, 2))
  )
  
  # Add SAIT results with residuals
  analysis@sait_results$sait_interaction <- list(
    results = data.frame(gene_id = rownames(counts), p_value = runif(n_genes)),
    residuals = matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples),
    fitted = matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  )
  
  # Add divergence results (SummarizedExperiment)
  div_counts <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(div_counts) <- rownames(counts)
  colnames(div_counts) <- colnames(counts)
  analysis@divergence_results$tsallis <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = div_counts)
  )
  
  analysis@metadata$function_calls <- c("TSENAT_config", "calculate_diversity")
  
  return(analysis)
}

# ===========================================================================
# Test 1: Basic subsetting with numeric indices
# ===========================================================================
test_that("Numeric subsetting works correctly for rows and columns", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Subset first 10 genes and first 5 samples
  subset_obj <- analysis[1:10, 1:5]
  
  expect_equal(nrow(se(subset_obj)), 10, info = "Row count should be 10")
  expect_equal(ncol(se(subset_obj)), 5, info = "Column count should be 5")
  expect_equal(nrow(subset_obj@diversity_results$q_0.5), 10)
  expect_equal(ncol(subset_obj@diversity_results$q_0.5), 5)
})

# ===========================================================================
# Test 2: Logical subsetting for indices
# ===========================================================================
test_that("Logical subsetting works for rows and columns", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Logical indexing
  row_logic <- rep(c(TRUE, FALSE), length.out = 20)  # Select every other gene
  col_logic <- rep(c(TRUE, FALSE), length.out = 10)   # Select every other sample
  
  subset_obj <- analysis[row_logic, col_logic]
  
  expect_equal(nrow(se(subset_obj)), 10, info = "Logical row subset should give 10 genes")
  expect_equal(ncol(se(subset_obj)), 5, info = "Logical column subset should give 5 samples")
})

# ===========================================================================
# Test 3: Character subsetting by gene names
# ===========================================================================
test_that("Character subsetting works by gene names", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  gene_names <- c("GENE_1", "GENE_5", "GENE_10")
  subset_obj <- analysis[gene_names, ]
  
  expect_equal(nrow(se(subset_obj)), 3, info = "Should have 3 selected genes")
  expect_equal(nrow(subset_obj@diversity_results$q_0.5), 3)
})

# ===========================================================================
# Test 4: Character subsetting by sample names
# ===========================================================================
test_that("Character subsetting works by sample names", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  sample_names <- c("Sample_1", "Sample_3", "Sample_7")
  subset_obj <- analysis[, sample_names]
  
  expect_equal(ncol(se(subset_obj)), 3, info = "Should have 3 selected samples")
  expect_equal(ncol(subset_obj@diversity_results$q_0.5), 3)
})

# ===========================================================================
# Test 5: ERROR - Invalid gene names not in object
# ===========================================================================
test_that("Error when gene names not found in object", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  invalid_genes <- c("GENE_1", "INVALID_GENE", "GENE_5")
  
  expect_error(
    analysis[invalid_genes, ],
    "Some gene names not found in object",
    info = "Should throw error for invalid gene names"
  )
})

# ===========================================================================
# Test 6: ERROR - Invalid sample names not in object
# ===========================================================================
test_that("Error when sample names not found in object", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  invalid_samples <- c("Sample_1", "INVALID_SAMPLE", "Sample_5")
  
  expect_error(
    analysis[, invalid_samples],
    "Some sample names not found in object",
    info = "Should throw error for invalid sample names"
  )
})

# ===========================================================================
# Test 7: ERROR - Row indices out of bounds (exceed nrow)
# ===========================================================================
test_that("Error when row indices exceed bounds", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  expect_error(
    analysis[1:25, ],  # Only 20 genes, asking for 25
    "Row indices out of bounds",
    info = "Should error when row index > nrow"
  )
})

# ===========================================================================
# Test 8: ERROR - Row indices out of bounds (zero or negative)
# ===========================================================================
test_that("Error when row indices are <= 0 or > nrow", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  expect_error(
    analysis[c(0, 1, 2), ],  # 0 is invalid
    "Row indices out of bounds",
    info = "Should error when row index is 0 or negative"
  )
})

# ===========================================================================
# Test 9: ERROR - Column indices out of bounds
# ===========================================================================
test_that("Error when column indices exceed bounds", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  expect_error(
    analysis[, 1:15],  # Only 10 samples, asking for 15
    "Column indices out of bounds",
    info = "Should error when column index > ncol"
  )
})

# ===========================================================================
# Test 10: ERROR - Column indices <= 0 or out of bounds
# ===========================================================================
test_that("Error when column indices are <= 0 or > ncol", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  expect_error(
    analysis[, c(-1, 1, 2)],  # Negative index
    "Column indices out of bounds",
    info = "Should error when column index is negative"
  )
})

# ===========================================================================
# Test 11: Default subsetting (no indices) returns full object
# ===========================================================================
test_that("Subsetting with missing indices returns full object", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Subset with no row or column indices
  full_subset <- analysis[, ]
  
  expect_equal(nrow(se(full_subset)), 20)
  expect_equal(ncol(se(full_subset)), 10)
})

# ===========================================================================
# Test 12: Subsetting with only row indices (no column index)
# ===========================================================================
test_that("Subsetting with only row indices works", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:5, ]
  
  expect_equal(nrow(se(subset_obj)), 5)
  expect_equal(ncol(se(subset_obj)), 10, info = "Should keep all samples")
})

# ===========================================================================
# Test 13: Subsetting with only column indices (no row index)
# ===========================================================================
test_that("Subsetting with only column indices works", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[, 1:3]
  
  expect_equal(nrow(se(subset_obj)), 20, info = "Should keep all genes")
  expect_equal(ncol(se(subset_obj)), 3)
})

# ===========================================================================
# Test 14: Subsetting preserves diversity results correctly
# ===========================================================================
test_that("Diversity results are subsetted correctly by both genes and samples", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:5, 1:3]
  
  # Check diversity_results subsetting
  expect_equal(nrow(subset_obj@diversity_results$q_0.5), 5)
  expect_equal(ncol(subset_obj@diversity_results$q_0.5), 3)
  expect_equal(nrow(subset_obj@diversity_results$q_1.0), 5)
  expect_equal(ncol(subset_obj@diversity_results$q_1.0), 3)
})

# ===========================================================================
# Test 15: Subsetting preserves jackknife results correctly
# ===========================================================================
test_that("Jackknife results are subsetted correctly (sample-level diagnostics)", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:5, 1:3]
  
  # Check jackknife resamples (should subset columns by samples)
  expect_equal(ncol(subset_obj@jackknife_results$q_0.5$resamples), 3)
  
  # Check influence scores (should subset rows to match samples)
  expect_equal(nrow(subset_obj@jackknife_results$q_0.5$influence_scores), 3)
  
  # Check CI matrix
  expect_equal(dim(subset_obj@jackknife_results$q_0.5$ci_matrix)[1], 5)  # genes
  expect_equal(dim(subset_obj@jackknife_results$q_0.5$ci_matrix)[2], 3)  # samples
})

# ===========================================================================
# Test 16: Subsetting preserves SAIT results correctly
# ===========================================================================
test_that("SAIT results are subsetted correctly (residuals and fitted values)", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:5, 1:3]
  
  # Check residuals (sample-level diagnostics)
  expect_equal(ncol(subset_obj@sait_results$sait_interaction$residuals), 3)
  
  # Check fitted values
  expect_equal(ncol(subset_obj@sait_results$sait_interaction$fitted), 3)
})

# ===========================================================================
# Test 17: Subsetting preserves divergence results (SummarizedExperiment)
# ===========================================================================
test_that("Divergence results (SummarizedExperiment) are subsetted correctly", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  ncol_div_before <- ncol(analysis@divergence_results$tsallis)
  
  subset_obj <- analysis[1:5, 1:3]
  
  # Divergence rows = genes (subset by i), columns = q-VALUES (preserved:
  # the second index refers to SAMPLES and must not shrink the q axis)
  expect_equal(nrow(subset_obj@divergence_results$tsallis), 5)
  expect_equal(ncol(subset_obj@divergence_results$tsallis), ncol_div_before)
})

# ===========================================================================
# Test 18: Subsetting preserves metadata and config
# ===========================================================================
test_that("Config and metadata are preserved after subsetting", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:5, 1:3]
  
  expect_equal(subset_obj@config$q_values, c(0.5, 1.0))
  expect_equal(subset_obj@metadata$function_calls, c("TSENAT_config", "calculate_diversity"))
})

# ===========================================================================
# Test 19: Single gene and single sample subsetting
# ===========================================================================
test_that("Subsetting to single gene or single sample works", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Single gene
  single_gene <- analysis[5, ]
  expect_equal(nrow(se(single_gene)), 1)
  expect_equal(ncol(se(single_gene)), 10)
  
  # Single sample
  single_sample <- analysis[, 3]
  expect_equal(nrow(se(single_sample)), 20)
  expect_equal(ncol(se(single_sample)), 1)
  
  # Single gene and sample
  single_both <- analysis[5, 3]
  expect_equal(nrow(se(single_both)), 1)
  expect_equal(ncol(se(single_both)), 1)
})

# ===========================================================================
# Test 20: Validation of subsetted object
# ===========================================================================
test_that("Subsetted object passes validity checks", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  subset_obj <- analysis[1:10, 1:5]
  
  # Should not throw error on validObject
  expect_true(validObject(subset_obj))
})

# ===========================================================================
# Test 21: Handle edge case: other result structures are preserved
# ===========================================================================
test_that("Other result structures (non-matrix, non-SE) are preserved", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Add a numeric vector result
  analysis@diversity_results$numeric_result <- rnorm(100)
  
  subset_obj <- analysis[1:5, 1:3]
  
  # Should preserve vector as-is
  expect_equal(length(subset_obj@diversity_results$numeric_result), 100)
})

# ===========================================================================
# Test 22: Handle edge case: unknown result structure in lapply
# ===========================================================================
test_that("Unknown result structures are preserved as-is", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Add a vector result
  analysis@diversity_results$vector_result <- 1:10
  
  subset_obj <- analysis[1:5, ]
  
  # Vector should be preserved unchanged
  expect_equal(subset_obj@diversity_results$vector_result, 1:10)
})

# ===========================================================================
# Test 23: Character subsetting preserves order
# ===========================================================================
test_that("Character subsetting preserves order of selected items", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Select genes in specific order (not sorted)
  selected_genes <- c("GENE_10", "GENE_2", "GENE_15")
  subset_obj <- analysis[selected_genes, ]
  
  result_names <- rownames(se(subset_obj))
  expect_equal(result_names, selected_genes)
})

# ===========================================================================
# Test 24: Mixed logical/character subsetting paths
# ===========================================================================
test_that("Mixed index types are converted correctly", {
  analysis <- .create_test_analysis(n_genes = 20, n_samples = 10)
  
  # Rows by character, columns by logical
  row_names <- paste0("GENE_", 1:5)
  col_logic <- rep(c(TRUE, FALSE), length.out = 10)
  
  subset_obj <- analysis[row_names, col_logic]
  
  expect_equal(nrow(se(subset_obj)), 5)
  expect_equal(ncol(se(subset_obj)), 5)
})

# ============================================================================
# BUG FIX 5: Verify pairwise_results slot removed from TSENATAnalysis
# ============================================================================

test_that("BUG 5: TSENATAnalysis does not have pairwise_results slot", {
  # The pairwise_results slot was removed as it was never read or written
  slot_names <- methods::slotNames("TSENATAnalysis")
  expect_false("pairwise_results" %in% slot_names,
               info = "pairwise_results slot should have been removed from S4 class")
})
