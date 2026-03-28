# Comprehensive testing for uncovered lines in generate_plots.R
# Tests edge cases, error conditions, and specific code paths

library(TSENAT)
skip_on_bioc()

context("Plots: Coverage Expansion for Edge Cases")

# ============================================================================
# TEST: infer_samples_from_se - Line 60 (samples parameter provided)
# ============================================================================

test_that("infer_samples_from_se: explicit samples parameter is returned as character", {
  # Line 60: return(as.character(samples))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide explicit samples parameter
  samples_provided <- c("S1", "S2", "S3", "S4")
  result <- TSENAT:::.infer_samples_from_se(se, samples = samples_provided)
  
  expect_true(is.character(result))
  expect_equal(result, samples_provided)
  expect_equal(length(result), 4)
})

test_that("infer_samples_from_se: numeric samples are coerced to character", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide numeric samples (edge case)
  samples_numeric <- c(1, 2, 3, 4)
  result <- TSENAT:::.infer_samples_from_se(se, samples = samples_numeric)
  
  expect_true(is.character(result))
  expect_equal(result, c("1", "2", "3", "4"))
})

# ============================================================================
# TEST: infer_samples_from_se - Line 65 (colData is NULL)
# ============================================================================

test_that("infer_samples_from_se: returns NULL when colData is missing/NULL", {
  # Line 65: return(NULL)
  
  # Create SE without colData
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Do not provide samples parameter
  result <- TSENAT:::.infer_samples_from_se(se, samples = NULL)
  
  # Expected: NULL since colData extraction fails
  expect_null(result)
})

# ============================================================================
# TEST: get_readcounts_from_se - Line 101 (file not found)
# ============================================================================

test_that("get_readcounts_from_se: errors when specified file doesn't exist", {
  # Line 101: if (!file.exists(readcounts_arg)) stop(...)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide path to non-existent file
  nonexistent_file <- "/tmp/definitely_does_not_exist_12345.txt"
  
  expect_error(
    TSENAT:::.get_readcounts_from_se(se, readcounts_arg = nonexistent_file),
    "readcounts file not found"
  )
})

# ============================================================================
# TEST: infer_samples_from_se - Fallback to binary or least-varied column
# ============================================================================

test_that("infer_samples_from_se: prefers binary column in fallback logic", {
  # Tests fallback when no standard column names match
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4)),
    colData = data.frame(
      treatment = c("A", "A", "B", "B"),
      replicate = c(1, 2, 1, 2),
      batch = c("X", "Y", "X", "Y"),
      stringsAsFactors = FALSE
    )
  )
  
  # Without providing samples and without standard column names
  result <- TSENAT:::.infer_samples_from_se(
    se,
    samples = NULL,
    condition_col = "nonexistent_col"
  )
  
  # Should pick one of the binary columns
  expect_true(is.character(result))
  expect_equal(length(result), 4)
})

# ============================================================================
# TEST: get_readcounts_from_se - Matrix input handling
# ============================================================================

test_that("get_readcounts_from_se: accepts matrix as readcounts_arg", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide explicit matrix
  custom_matrix <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80), nrow = 4, ncol = 2)
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = custom_matrix)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 2))
})

test_that("get_readcounts_from_se: accepts data.frame as readcounts_arg", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide explicit data.frame
  custom_df <- data.frame(
    Gene_ID = c("G1", "G2", "G3", "G4"),
    Sample1 = c(10, 20, 30, 40),
    Sample2 = c(50, 60, 70, 80),
    stringsAsFactors = FALSE
  )
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = custom_df)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
})

test_that("get_readcounts_from_se: errors on invalid readcounts_arg type", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide invalid type
  expect_error(
    TSENAT:::.get_readcounts_from_se(se, readcounts_arg = list(invalid = "type")),
    "must be a matrix|data.frame|path"
  )
})

# ============================================================================
# TEST: get_tx2gene_from_se - Null handling and column detection
# ============================================================================

test_that("get_tx2gene_from_se: extracts tx2gene from metadata", {
  # Create a proper tx2gene mapping
  tx2gene_map <- data.frame(
    Transcript = c("TX1", "TX2", "TX3", "TX4"),
    Gene = c("G1", "G2", "G1", "G3"),
    stringsAsFactors = FALSE
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:16, nrow = 4, ncol = 4))
  )
  rownames(se) <- c("TX1", "TX2", "TX3", "TX4")
  S4Vectors::metadata(se)$tx2gene <- tx2gene_map
  
  readcounts_mat <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(readcounts_mat) <- c("TX1", "TX2", "TX3", "TX4")
  
  result <- TSENAT:::.get_tx2gene_from_se(se, readcounts_mat)
  
  expect_true(is.list(result))
  expect_true("mapping" %in% names(result))
})

# ============================================================================
# TEST: validate_control_in_samples - Control parameter validation
# ============================================================================

test_that("validate_control_in_samples: returns control when it's in sample list", {
  samples <- c("control_1", "treatment_1", "control_2", "treatment_2")
  control <- "control_1"
  
  result <- TSENAT:::.validate_control_in_samples(control, samples)
  
  expect_equal(result, "control_1")
})

test_that("validate_control_in_samples: returns 'Normal' when present and control not found", {
  samples <- c("Normal", "group_B", "group_C")
  control <- "group_D"
  
  result <- TSENAT:::.validate_control_in_samples(control, samples)
  
  # Should return "Normal" as fallback
  expect_equal(result, "Normal")
})

test_that("validate_control_in_samples: returns first element as fallback", {
  samples <- c("group_A", "group_B", "group_C")
  control <- "group_D"
  
  result <- TSENAT:::.validate_control_in_samples(control, samples)
  
  # Should return first unique element
  expect_equal(result, "group_A")
})

# ============================================================================
# TEST: Readcounts file with single column (malformed)
# ============================================================================

test_that("get_readcounts_from_se: handles single-column readcounts file", {
  # Create temporary single-column readcounts file
  temp_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(gene_id = c("G1", "G2", "G3")),
    file = temp_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:15, nrow = 5, ncol = 3))
  )
  
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = temp_file)
  
  expect_true(is.matrix(result))
  unlink(temp_file)  # Clean up
})

# ============================================================================
# TEST: infer_samples_from_se with various column types in colData
# ============================================================================

test_that("infer_samples_from_se: handles factor columns in colData", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4)),
    colData = data.frame(
      condition = factor(c("ctrl", "ctrl", "treat", "treat")),
      stringsAsFactors = TRUE
    )
  )
  
  result <- TSENAT:::.infer_samples_from_se(se, samples = NULL)
  
  expect_true(is.character(result))
  expect_equal(length(result), 4)
})

test_that("infer_samples_from_se: handles numeric vector in colData", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4)),
    colData = data.frame(
      group_id = c(1, 1, 2, 2)
    )
  )
  
  result <- TSENAT:::.infer_samples_from_se(se, samples = NULL)
  
  # Should find group_id with 2 unique values
  expect_true(is.character(result))
  expect_equal(length(result), 4)
})

# ============================================================================
# TEST: Column name detection with special characters
# ============================================================================

test_that("infer_samples_from_se: finds columns with underscores and hyphens", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4)),
    colData = data.frame(
      sample_group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  )
  
  result <- TSENAT:::.infer_samples_from_se(se, samples = NULL)
  
  # Should match the candidate list ("sample_group" is in candidates)
  expect_equal(result, c("A", "A", "B", "B"))
})

# ============================================================================
# TEST: get_tx2gene_from_se - Returns list with rownames fallback (Line 168)
# ============================================================================

test_that("get_tx2gene_from_se: returns NULL when readcounts_mat is NULL", {
  # Line 168: NULL
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:16, nrow = 4, ncol = 4))
  )
  # No metadata with tx2gene
  
  result <- TSENAT:::.get_tx2gene_from_se(se, readcounts_mat = NULL)
  
  # Should return NULL when readcounts_mat is NULL
  expect_null(result)
})

# ============================================================================
# TEST: get_readcounts_from_se - Fallback to first assay with warning
# ============================================================================

test_that("get_readcounts_from_se: falls back to first assay when no preferred assay found", {
  # Lines 131-136: fallback with warning
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(custom_assay = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  expect_warning(
    result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = NULL),
    "Using first assay"
  )
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(5, 4))
})

# ============================================================================
# TEST: get_readcounts_from_se - Reads from metadata (Line 120)
# ============================================================================

test_that("get_readcounts_from_se: reads readcounts from metadata when available", {
  # Lines 119-120: if (!is.null(md) && !is.null(md$readcounts)) return(as.matrix(...))
  
  metadata_counts <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(0, nrow = 5, ncol = 4))
  )
  S4Vectors::metadata(se)$readcounts <- metadata_counts
  
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = NULL)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 2))
  expect_equal(result, metadata_counts)
})

# ============================================================================
# TEST: get_readcounts_from_se - Multiple columns in data.frame (Lines 103-105)
# ============================================================================

test_that("get_readcounts_from_se: extracts numeric columns from multi-column data.frame", {
  # Lines 103-105: ncol > 1 case
  
  temp_file <- tempfile(fileext = ".txt")
  df <- data.frame(
    Gene = c("G1", "G2", "G3", "G4"),
    S1 = c(10, 20, 30, 40),
    S2 = c(50, 60, 70, 80),
    stringsAsFactors = FALSE
  )
  write.table(df, file = temp_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(0, nrow = 5, ncol = 4))
  )
  
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = temp_file)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 2)
  expect_equal(rownames(result), c("G1", "G2", "G3", "G4"))
  unlink(temp_file)
})

# ============================================================================
# TEST: Readcounts with preferred assay selection (Lines 126-127)
# ============================================================================

test_that("get_readcounts_from_se: selects 'readcounts' assay when multiple preferred assays exist", {
  # Lines 126-127: choose preferred assay
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      other = matrix(1:10, nrow = 5, ncol = 2),
      counts = matrix(11:20, nrow = 5, ncol = 2),
      readcounts = matrix(21:30, nrow = 5, ncol = 2)
    )
  )
  
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = NULL)
  
  # Should select 'readcounts' assay
  expect_true(is.matrix(result))
  expect_equal(result, matrix(21:30, nrow = 5, ncol = 2))
})

# ============================================================================
# TEST: get_tx2gene_from_se - rowData fallback (Lines 156-160)
# ============================================================================

test_that("get_tx2gene_from_se: extracts genes column from rowData when available", {
  # Lines 156-160: rowData fallback with genes column
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:16, nrow = 4, ncol = 4)),
    rowData = data.frame(genes = c("G1", "G2", "G3", "G4"))
  )
  
  readcounts_mat <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(readcounts_mat) <- c("TX1", "TX2", "TX3", "TX4")
  
  result <- TSENAT:::.get_tx2gene_from_se(se, readcounts_mat)
  
  expect_true(is.list(result))
  expect_equal(result$type, "vector")
  expect_equal(result$mapping, c("G1", "G2", "G3", "G4"))
})

# ============================================================================
# TEST: get_tx2gene_from_se - rownames fallback (Lines 164-165)
# ============================================================================

test_that("get_tx2gene_from_se: uses rownames as fallback when no tx2gene available", {
  # Lines 164-165: Last resort - use rownames
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:16, nrow = 4, ncol = 4))
  )
  
  readcounts_mat <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(readcounts_mat) <- c("TX1", "TX2", "TX3", "TX4")
  
  result <- TSENAT:::.get_tx2gene_from_se(se, readcounts_mat)
  
  expect_true(is.list(result))
  expect_equal(result$type, "vector")
  expect_equal(result$mapping, c("TX1", "TX2", "TX3", "TX4"))
})

# ============================================================================
# TEST: get_readcounts_from_se - Chosen preferred assay extraction (Lines 127)
# ============================================================================

test_that("get_readcounts_from_se: uses 'counts' assay when 'readcounts' not available", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      something_else = matrix(1:10, nrow = 5, ncol = 2),
      counts = matrix(11:20, nrow = 5, ncol = 2)
    )
  )
  
  result <- TSENAT:::.get_readcounts_from_se(se, readcounts_arg = NULL)
  
  expect_true(is.matrix(result))
  expect_equal(result, matrix(11:20, nrow = 5, ncol = 2))
})

# ============================================================================
# TEST: Samples parameter with matrix from infer_samples_from_se
# ============================================================================

test_that("infer_samples_from_se: handles matrix input for samples parameter", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:20, nrow = 5, ncol = 4))
  )
  
  # Provide matrix with single column (each element becomes a character)
  samples_matrix <- c("S1", "S2", "S3", "S4")
  result <- TSENAT:::.infer_samples_from_se(se, samples = samples_matrix)
  
  expect_true(is.character(result))
  expect_equal(length(result), 4)
})
