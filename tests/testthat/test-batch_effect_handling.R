# Add required libraries for tests
library(SummarizedExperiment)
library(TSENAT)

context("batch_effect_handling: Detection and Correction Methods")

# Helper to convert readcounts data.frame to SummarizedExperiment
.make_test_se <- function() {
  # Load salmon_dataset from package data
  load(system.file("data", "readcounts.RData", package = "TSENAT"))
  
  # Create metadata
  sample_ids <- colnames(salmon_dataset)
  col_data <- data.frame(
    sample = sample_ids,
    sample_type = factor(c(rep("normal", 8), rep("tumor", 8))),
    row.names = sample_ids
  )
  
  # Create SE
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset)),
    colData = col_data
  )
  
  return(se)
}

test_that("detect_batch_effects runs without error", {
  se <- .make_test_se()
  SummarizedExperiment::colData(se)$batch_id <- factor(c(rep("A", 8), rep("B", 8)))
  
  result <- detect_batch_effects(
    se = se,
    batch = "batch_id",
    method = "pca",
    n_components = 3,
    n_permutations = 10
  )
  
  expect_s3_class(result, "batch_detection")
  expect_true(all(!is.na(result$batch_variance_pct)))
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
  expect_true(length(result$batch_variance_pct) == 3)
})

test_that("adjust_batch_effects runs standard method without error", {
  se <- .make_test_se()
  SummarizedExperiment::colData(se)$batch_id <- factor(c(rep("A", 8), rep("B", 8)))
  
  result <- adjust_batch_effects(
    se = se,
    batch = "batch_id",
    method = "standard",
    shrinkage = TRUE,
    par.prior = TRUE
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_equal(dim(result), dim(se))
  expect_true("batch_correction" %in% names(metadata(result)))
  expect_true(metadata(result)$batch_correction$method == "ComBat-seq")
})

test_that("adjust_batch_effects runs reference method without error", {
  se <- .make_test_se()
  SummarizedExperiment::colData(se)$batch_id <- factor(c(rep("A", 8), rep("B", 8)))
  
  result <- adjust_batch_effects(
    se = se,
    batch = "batch_id",
    method = "ref",
    ref_batch = "A",
    mean.only = FALSE
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_equal(dim(result), dim(se))
  expect_true(metadata(result)$batch_correction$method == "ComBat-ref")
})

test_that("invalid inputs raise errors", {
  se <- .make_test_se()
  
  expect_error(
    detect_batch_effects(se, batch = "nonexistent")
  )
  
  SummarizedExperiment::colData(se)$batch_id <- factor(c(rep("A", 8), rep("B", 8)))
  expect_error(
    adjust_batch_effects(se, batch = "batch_id", method = "ref", ref_batch = "Z")
  )
})

# ============================================================================
# Rank-Based Batch Effect Detection and Correction Tests
# ============================================================================

test_that("detect_batch_structure handles matrix input", {
  # Create test entropy matrix
  entropy_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:40)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  result <- detect_batch_structure(
    entropy_lists = entropy_matrix,
    sample_metadata = NULL,
    n_pcs = 3
  )
  
  expect_s3_class(result, "batch_pca")
  # Note: n_pcs is adjusted to min(requested, n_samples-1), so with 5 samples, min(3, 4) = 3
  # But the function may return up to n_pcs components
  expect_true(length(result$variance_explained) >= 1)
  expect_true(length(result$variance_explained) <= 5)
  expect_equal(nrow(result$batch_pca_scores), 5)
  expect_true(all(!is.na(result$variance_explained)))
  expect_true(all(result$variance_explained >= 0))
  expect_true(sum(result$variance_explained) <= 1.01)  # Allow small numerical error
})

test_that("detect_batch_structure with sample metadata", {
  entropy_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:40)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  metadata <- data.frame(
    sample_id = paste0("sample", 1:5),
    batch = factor(c("A", "A", "B", "B", "A")),
    condition = factor(c("normal", "tumor", "normal", "tumor", "normal"))
  )
  
  result <- detect_batch_structure(
    entropy_lists = entropy_matrix,
    sample_metadata = metadata,
    color_by = "batch"
  )
  
  expect_s3_class(result, "batch_pca")
  expect_true("batch" %in% colnames(result$batch_pca_scores))
  expect_true("condition" %in% colnames(result$batch_pca_scores))
  expect_true(is.logical(result$is_batch_confounded))
  expect_true(is.numeric(result$batch_effect_strength))
})

test_that("detect_batch_structure handles list of matrices", {
  # Create multiple entropy matrices for different q-values
  q01_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  q05_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  q10_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  
  rownames(q01_matrix) <- paste0("gene", 1:40)
  colnames(q01_matrix) <- paste0("sample", 1:5)
  rownames(q05_matrix) <- rownames(q01_matrix)
  colnames(q05_matrix) <- colnames(q01_matrix)
  rownames(q10_matrix) <- rownames(q01_matrix)
  colnames(q10_matrix) <- colnames(q01_matrix)
  
  entropy_list <- list(q01 = q01_matrix, q05 = q05_matrix, q10 = q10_matrix)
  
  result <- detect_batch_structure(
    entropy_lists = entropy_list,
    n_pcs = 2
  )
  
  expect_s3_class(result, "batch_pca")
  # Check that variance explained is reasonable
  expect_true(length(result$variance_explained) >= 1)
  expect_true(length(result$variance_explained) <= 5)
  expect_equal(nrow(result$entropy_data), 40)
  expect_equal(ncol(result$entropy_data), 5)
})

test_that("detect_batch_structure detects strong batch effects", {
  # Create data with strong batch signal
  set.seed(42)
  entropy_matrix <- matrix(rnorm(200, sd = 0.5), nrow = 40, ncol = 5)
  
  # Add strong batch effect to first 3 samples
  entropy_matrix[, 1:3] <- entropy_matrix[, 1:3] + 3
  
  rownames(entropy_matrix) <- paste0("gene", 1:40)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  metadata <- data.frame(
    sample_id = paste0("sample", 1:5),
    batch = factor(c("A", "A", "A", "B", "B"))
  )
  
  result <- detect_batch_structure(
    entropy_lists = entropy_matrix,
    sample_metadata = metadata,
    color_by = "batch",
    n_pcs = 2
  )
  
  expect_s3_class(result, "batch_pca")
  expect_true(is.logical(result$is_batch_confounded))
  expect_true(result$batch_effect_strength > 0)
})

test_that("apply_batch_correction_ranking processes matrix correctly", {
  # Create test entropy matrix
  entropy_matrix <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:20)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  batch_factor <- factor(c("A", "A", "A", "B", "B"))
  
  result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor
  )
  
  expect_is(result, "list")
  expect_true("entropy_corrected" %in% names(result))
  expect_true("batch_effects" %in% names(result))
  expect_true("model_fits" %in% names(result))
  expect_true("r_squared_by_gene" %in% names(result))
  expect_equal(dim(result$entropy_corrected), dim(entropy_matrix))
})

test_that("apply_batch_correction_ranking with condition factor", {
  entropy_matrix <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:20)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  batch_factor <- factor(c("A", "A", "A", "B", "B"))
  condition_factor <- factor(c("control", "case", "control", "case", "control"))
  
  result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor,
    condition_factor = condition_factor
  )
  
  expect_is(result, "list")
  expect_equal(nrow(result$entropy_corrected), 20)
  expect_equal(ncol(result$entropy_corrected), 5)
  expect_true(length(result$batch_effects) <= 20)
})

test_that("apply_batch_correction_ranking r_squared values are valid", {
  entropy_matrix <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:20)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  batch_factor <- factor(c("A", "A", "B", "B", "A"))
  
  result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor
  )
  
  expect_true(length(result$r_squared_by_gene) == 20)
  expect_true(all(result$r_squared_by_gene >= 0, na.rm = TRUE))
  expect_true(all(result$r_squared_by_gene <= 1, na.rm = TRUE))
  expect_true(is.numeric(result$mean_r_squared))
})

test_that("batch correction preserves matrix structure", {
  entropy_matrix <- matrix(1:30, nrow = 6, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:6)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  batch_factor <- factor(c("A", "A", "B", "B", "B"))
  
  result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor
  )
  
  # Same dimensions
  expect_equal(dim(result$entropy_corrected), dim(entropy_matrix))
  
  # Same rownames and colnames
  expect_equal(rownames(result$entropy_corrected), rownames(entropy_matrix))
  expect_equal(colnames(result$entropy_corrected), colnames(entropy_matrix))
})

test_that("batch correction handles different batch level counts", {
  entropy_matrix <- matrix(rnorm(150), nrow = 30, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:30)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  # Unbalanced batches
  batch_factor <- factor(c("A", "A", "B", "B", "B"))
  
  result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor
  )
  
  expect_equal(length(result$batch_levels), 2)
  expect_true(all(c("A", "B") %in% result$batch_levels))
  expect_equal(length(result$r_squared_by_gene), 30)
})

test_that("batch correction requires multiple batch levels", {
  entropy_matrix <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:20)
  
  # All same batch - should raise error because lm needs contrast levels
  batch_factor <- factor(rep("A", 5))
  
  # Single batch level should cause an error in lm
  expect_error(
    apply_batch_correction_ranking(
      entropy_matrix = entropy_matrix,
      batch_factor = batch_factor
    )
  )
})

test_that("batch detection works with SummarizedExperiment input", {
  # Create SE with entropy in assay
  entropy_matrix <- matrix(rnorm(200), nrow = 40, ncol = 5)
  rownames(entropy_matrix) <- paste0("gene", 1:40)
  colnames(entropy_matrix) <- paste0("sample", 1:5)
  
  col_data <- data.frame(
    sample_id = paste0("sample", 1:5),
    batch = factor(c("A", "A", "B", "B", "A")),
    row.names = paste0("sample", 1:5)
  )
  
  se <- SummarizedExperiment(
    assays = list(entropy = entropy_matrix),
    colData = col_data
  )
  
  result <- detect_batch_structure(
    entropy_lists = se,
    sample_metadata = as.data.frame(colData(se)),
    color_by = "batch"
  )
  
  expect_s3_class(result, "batch_pca")
  expect_equal(nrow(result$batch_pca_scores), 5)
  expect_true("batch" %in% colnames(result$batch_pca_scores))
})

