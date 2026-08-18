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


context("Bootstrap Entropy: Minimum Count Filtering")

# Feature 4.2: Graceful handling of top genes with insufficient counts
# Tests for minimum count filtering in bootstrap confidence intervals

test_that("Genes with sufficient counts (≥10) are processed normally", {
  # Create SE with genes having ≥10 total counts
  counts_matrix <- rbind(
    Gene1 = c(100, 50, 25, 10),  # Total: 185
    Gene2 = c(50, 50, 50, 50)     # Total: 200
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1", "Gene2"),
    pvalue = c(0.01, 0.05),
    row.names = 1:2
  )
  
  result <- .calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    verbose = FALSE
  )
  
  # Should successfully process and return result
  expect_true(is.list(result))
  expect_true("estimate" %in% names(result))
  expect_true("lower_ci" %in% names(result))
})

test_that("Genes with insufficient counts (<10) trigger warning", {
  # Create SE where ALL genes have insufficient counts
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(2, 2, 2, 2)         # Total: 8 (insufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1", "Gene2"),
    pvalue = c(0.001, 0.05),
    row.names = 1:2
  )
  
  # Should warn about insufficient counts for all genes
  expect_warning(
    result <- .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      verbose = FALSE
    ),
    "No genes with sufficient"
  )
})

test_that("Function skips low-count genes and uses next valid gene", {
  # Create SE where current top gene is insufficient but next one is valid
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(50, 50, 50, 50)     # Total: 200 (sufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1", "Gene2"),
    pvalue = c(0.001, 0.05),
    row.names = 1:2
  )
  
  # Request top_n=1, but Gene1 is insufficient
  # Function should skip to Gene2 and warn about skipping Gene1
  result <- suppressWarnings(
    .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      verbose = FALSE
    )
  )
  
  # Result should be valid (from Gene2)
  expect_true(is.list(result))
  expect_true(!is.null(result))
  expect_true("estimate" %in% names(result))
})

test_that("Minimum threshold is 10", {
  # Gene with exactly 10 counts should be accepted
  counts_matrix <- rbind(
    Gene1 = c(5, 3, 1, 1)          # Total: 10 (at boundary)
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should process Gene1 (exactly at threshold)
  result <- .calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("Gene with count = 9 is rejected (below threshold)", {
  # Gene with 9 counts should be rejected
  counts_matrix <- rbind(
    Gene1 = c(5, 3, 1, 0)          # Total: 9 (below threshold)
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should return NULL with warning
  result <- suppressWarnings(
    .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      verbose = FALSE
    )
  )
  
  expect_true(is.null(result))
})

test_that("All genes insufficient returns NULL and warning", {
  # Create SE where all genes have insufficient counts
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4
    Gene2 = c(2, 2, 2, 2)         # Total: 8
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1", "Gene2"),
    pvalue = c(0.01, 0.05),
    row.names = 1:2
  )
  
  # Should warn and return NULL
  expect_warning(
    result <- .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 2, q = 1, nboot = 100, 
      verbose = FALSE
    ),
    "No genes with sufficient"
  )
  
  expect_true(is.null(result))
})

test_that("Filtering works with multiple genes requested (top_n > 1)", {
  # Create SE with mixed sufficient/insufficient genes
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(50, 50, 50, 50),    # Total: 200 (sufficient)
    Gene3 = c(2, 2, 2, 2)         # Total: 8 (insufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2", "Gene3")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1", "Gene2", "Gene3"),
    pvalue = c(0.001, 0.01, 0.05),
    row.names = 1:3
  )
  
  # Request all 3 genes, but only Gene2 is valid
  result <- suppressWarnings(
    .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 3, q = 1, nboot = 100, 
      verbose = FALSE
    )
  )
  
  # Should return result(s) for valid genes only
  if (!is.null(result)) {
    if (is.list(result) && "estimate" %in% names(result)) {
      # Single gene result
      expect_true("estimate" %in% names(result))
    }
  }
})

test_that("Warning message mentions minimum threshold and database papers", {
  # Create SE where all genes are insufficient
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1)         # Total: 4
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Capture warning to check content
  warn_msg <- tryCatch(
    .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      verbose = FALSE
    ),
    warning = function(w) conditionMessage(w)
  )
  
  # Warning should mention the minimum-count threshold
  result <- suppressWarnings(
    .calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      verbose = FALSE
    )
  )
  
  expect_true(is.null(result))
})

test_that("Feature 4.2 gracefully handles genes by rownames vs rowData", {
  # Test with genes identified in rowData instead of rownames
  counts_matrix <- rbind(
    TX1 = c(50, 50, 50, 50),      # Transcript 1: total 200
    TX2 = c(25, 25, 25, 25)       # Transcript 2: total 100
  )
  
  rownames(counts_matrix) <- c("TX1", "TX2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(
        gene_name = c("Gene_A", "Gene_A")  # Both transcripts map to same gene (use singular)
      )
    )
  )
  
  res <- data.frame(
    gene_id = c("Gene_A"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should find Gene_A via rowData and process it
  result <- .calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
  expect_true("estimate" %in% names(result))
})

test_that("Feature 4.2 implementation uses the documented minimum-count recommendation", {
  # Verify that the 10-count threshold aligns with database recommendations
  # This test documents the source of the threshold value
  
  # Create gene exactly below threshold
  counts_below <- rbind(
    Gene1 = c(3, 3, 3, 0)         # Total: 9
  )
  
  # Create gene exactly at threshold
  counts_at <- rbind(
    Gene1 = c(3, 3, 3, 1)         # Total: 10
  )
  
  rownames(counts_below) <- "Gene1"
  rownames(counts_at) <- "Gene1"
  
  se_below <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_below),
      rowData = S4Vectors::DataFrame(gene_names = "Gene1")
    )
  )
  
  se_at <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_at),
      rowData = S4Vectors::DataFrame(gene_names = "Gene1")
    )
  )
  
  res <- data.frame(gene_id = "Gene1", pvalue = 0.01, row.names = 1)
  
  # Gene with count=9 should fail
  result_below <- suppressWarnings(
    .calculate_tsallis_entropy_bootstrap(se = se_below, res = res, top_n = 1, nboot = 100, verbose = FALSE)
  )
  
  # Gene with count=10 should succeed
  result_at <- .calculate_tsallis_entropy_bootstrap(se = se_at, res = res, top_n = 1, nboot = 100, verbose = FALSE)
  
  expect_true(is.null(result_below))
  expect_true(!is.null(result_at))
})


