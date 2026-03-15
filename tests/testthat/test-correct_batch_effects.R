# Tests for unified correct_batch_effects function
library(SummarizedExperiment)
library(TSENAT)

context("correct_batch_effects: Unified Batch Correction Interface")

# Helper to create test SummarizedExperiment with batch variable
.make_batch_test_se <- function(n_genes = 50, n_samples = 20, seed = 42) {
  set.seed(seed)
  counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 2), 
                  nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("gene_", 1:n_genes)
  colnames(counts) <- paste0("sample_", 1:n_samples)
  
  batch <- factor(c(rep("batch1", n_samples/2), rep("batch2", n_samples/2)))
  
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = DataFrame(batch = batch)
  )
  
  return(se)
}

# ============================================================================
# Test: correct_batch_effects with ComBat method
# ============================================================================

test_that("correct_batch_effects with combat method returns valid SE", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  
  # Check output type
  expect_s4_class(result, "SummarizedExperiment")
  
  # Check dimensions preserved
  expect_equal(nrow(result), nrow(se))
  expect_equal(ncol(result), ncol(se))
  
  # Check batch_corrected assay exists
  expect_true("batch_corrected" %in% names(assays(result)))
  
  # Check original counts assay preserved
  expect_true("counts" %in% names(assays(result)))
  
  # Check colData preserved
  expect_equal(colData(result)$batch, colData(se)$batch)
})

test_that("correct_batch_effects with combat_seq method returns valid SE", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat_seq"
  )
  
  # Check output type
  expect_s4_class(result, "SummarizedExperiment")
  
  # Check batch_corrected assay exists
  expect_true("batch_corrected" %in% names(assays(result)))
  
  # Verify corrected matrix is numeric
  corrected <- assay(result, "batch_corrected")
  expect_true(is.numeric(corrected))
  expect_equal(dim(corrected), c(50, 20))
})

test_that("correct_batch_effects with combat_ref method returns valid SE", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat_ref",
    reference_batch = "batch1"
  )
  
  # Check output type
  expect_s4_class(result, "SummarizedExperiment")
  
  # Check batch_corrected assay exists
  expect_true("batch_corrected" %in% names(assays(result)))
  
  # Verify corrected matrix dimensions
  corrected <- assay(result, "batch_corrected")
  expect_equal(dim(corrected), dim(assay(se)))
})

test_that("correct_batch_effects with ranking method returns valid SE", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "ranking"
  )
  
  # Check output type
  expect_s4_class(result, "SummarizedExperiment")
  
  # Check batch_corrected assay exists
  expect_true("batch_corrected" %in% names(assays(result)))
  
  # Check original counts preserved
  expect_true("counts" %in% names(assays(result)))
  
  # Verify corrected matrix is numeric
  corrected <- assay(result, "batch_corrected")
  expect_true(is.numeric(corrected))
})

# ============================================================================
# Test: Metadata tracking
# ============================================================================

test_that("correct_batch_effects stores metadata correctly", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  
  # Check metadata exists and contains required fields
  meta <- S4Vectors::metadata(result)
  expect_true("batch_correction_method" %in% names(meta))
  expect_true("batch_correction_date" %in% names(meta))
  expect_true("batch_column" %in% names(meta))
  
  # Verify method is recorded
  expect_equal(meta$batch_correction_method, "combat")
  expect_equal(meta$batch_column, "batch")
  
  # Verify date is today
  expect_equal(meta$batch_correction_date, Sys.Date())
})

# ============================================================================
# Test: Error handling and validation
# ============================================================================

test_that("correct_batch_effects throws error for non-SE input", {
  counts <- matrix(rnbinom(1000, mu = 100, size = 2), nrow = 50, ncol = 20)
  batch <- factor(c(rep("B1", 10), rep("B2", 10)))
  
  # Should error on matrix input (requires SE)
  expect_error(
    correct_batch_effects(x = counts, batch_column = "batch", method = "combat"),
    "must be a SummarizedExperiment"
  )
})

test_that("correct_batch_effects throws error for missing batch column", {
  se <- .make_batch_test_se()
  
  # Try to use non-existent column
  expect_error(
    correct_batch_effects(x = se, batch_column = "nonexistent", method = "combat"),
    "not found in colData"
  )
})

test_that("correct_batch_effects throws error for invalid method", {
  se <- .make_batch_test_se()
  
  # Try invalid method (this should trigger match.arg error)
  expect_error(
    correct_batch_effects(x = se, batch_column = "batch", method = "invalid_method"),
    "should be one of"  # match.arg error message
  )
})

# ============================================================================
# Test: Assay structure and data preservation
# ============================================================================

test_that("correct_batch_effects creates batch_corrected assay alongside original", {
  se <- .make_batch_test_se()
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  
  # Both assays should exist
  expect_true("counts" %in% names(assays(result)))
  expect_true("batch_corrected" %in% names(assays(result)))
  
  # batch_corrected should be numeric
  corrected <- assay(result, "batch_corrected")
  expect_true(is.numeric(corrected))
})

test_that("correct_batch_effects produces different corrected matrix", {
  se <- .make_batch_test_se()
  original_counts <- assay(se, "counts")
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  
  # Batch-corrected matrix should be different from original
  corrected <- assay(result, "batch_corrected")
  expect_false(identical(original_counts, corrected))
  
  # But should have same dimensions
  expect_equal(dim(original_counts), dim(corrected))
})

# ============================================================================
# Test: Multiple batches handling
# ============================================================================

test_that("correct_batch_effects handles multiple batches", {
  set.seed(42)
  counts <- matrix(rnbinom(500, mu = 100, size = 2), nrow = 50, ncol = 30)
  rownames(counts) <- paste0("gene_", 1:50)
  colnames(counts) <- paste0("sample_", 1:30)
  
  # Create 3 batches
  batch <- factor(c(rep("B1", 10), rep("B2", 10), rep("B3", 10)))
  
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = DataFrame(batch = batch)
  )
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("batch_corrected" %in% names(assays(result)))
  expect_equal(ncol(result), 30)
})

# ============================================================================
# Test: Unbalanced batch design
# ============================================================================

test_that("correct_batch_effects handles unbalanced batches", {
  set.seed(42)
  counts <- matrix(rnbinom(400, mu = 100, size = 2), nrow = 50, ncol = 8)
  rownames(counts) <- paste0("gene_", 1:50)
  colnames(counts) <- paste0("sample_", 1:8)
  
  # Unbalanced batches (5 and 3)
  batch <- factor(c(rep("B1", 5), rep("B2", 3)))
  
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = DataFrame(batch = batch)
  )
  
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "ranking"
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("batch_corrected" %in% names(assays(result)))
})

# ============================================================================
# Test: Reference batch specification
# ============================================================================

test_that("correct_batch_effects handles reference_batch parameter for combat_ref", {
  se <- .make_batch_test_se()
  
  # Should work with explicit reference
  result <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat_ref",
    reference_batch = "batch1"
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("batch_corrected" %in% names(assays(result)))
})

test_that("correct_batch_effects requires reference_batch for combat_ref method", {
  se <- .make_batch_test_se()
  
  # combat_ref requires an explicit reference_batch
  expect_error(
    correct_batch_effects(
      x = se,
      batch_column = "batch",
      method = "combat_ref",
      reference_batch = NULL
    ),
    "ref_batch"
  )
})

# ============================================================================
# Test: Method consistency (same data with same method should be consistent)
# ============================================================================

test_that("correct_batch_effects is deterministic for same input", {
  set.seed(42)
  se <- .make_batch_test_se(seed = 42)
  
  set.seed(123)  # Different seed for function execution
  result1 <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "ranking"
  )
  
  set.seed(456)  # Different seed again
  result2 <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "ranking"
  )
  
  # Results should be identical (ranking is deterministic)
  corrected1 <- assay(result1, "batch_corrected")
  corrected2 <- assay(result2, "batch_corrected")
  expect_identical(corrected1, corrected2)
})

# ============================================================================
# Test: Batch column case sensitivity
# ============================================================================

test_that("correct_batch_effects uses specified batch column correctly", {
  se <- .make_batch_test_se()
  
  # Add another batch column
  colData(se)$batch_id <- colData(se)$batch
  
  # Should work with correct column name
  result1 <- correct_batch_effects(
    x = se,
    batch_column = "batch",
    method = "combat"
  )
  expect_s4_class(result1, "SummarizedExperiment")
  
  # Should also work with the other batch column
  result2 <- correct_batch_effects(
    x = se,
    batch_column = "batch_id",
    method = "combat"
  )
  expect_s4_class(result2, "SummarizedExperiment")
  
  # Verify different columns were used (metadata should reflect this)
  meta1 <- S4Vectors::metadata(result1)
  meta2 <- S4Vectors::metadata(result2)
  expect_equal(meta1$batch_column, "batch")
  expect_equal(meta2$batch_column, "batch_id")
})
