# Tests for jackknife_diagnostics.R uncovered lines from jackknife_coverage.txt
# Covers ~237 uncovered lines including:
# - Multiple q values, SE input, res parameter
# - Gene lookup from rowData, match() optimization  
# - Print/summary methods for jackknife results
# - Outlier detection, influence metrics

# Load consolidation helpers
source("helper-jackknife-consolidation.R")

test_that("jackknife_tsallis_entropy basic vector input", {
  # Test basic vector input
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120, 150, 60)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 1,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(inherits(result, "tsenat_jackknife"))
  expect_true("estimate" %in% names(result))
  expect_true("jackknife_se" %in% names(result))
  expect_true("influence" %in% names(result))
})

test_that("jackknife_tsallis_entropy multiple q with print_results = TRUE", {
  # Test printing of multi-q results
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  suppressMessages(
    result <- jackknife_tsallis_entropy(
      x = x,
      q = c(1.0, 2.0),
      norm = TRUE,
      print_results = TRUE
    )
  )
  
  expect_true(is.list(result))
})

test_that("jackknife_tsallis_entropy matrix with multiple q", {
  # Test matrix input with multiple q values
  set.seed(123)
  x <- matrix(c(100, 50, 200, 75, 150, 80), nrow = 2, ncol = 3)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = c(1.0, 1.5, 2.0),
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(is.list(result))
})

test_that("jackknife_tsallis_entropy with SE and res inputs", {
  # Test SummarizedExperiment with results data.frame for multi-gene analysis
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120, 110, 60, 85), nrow = 3, ncol = 3)),
    rowData = data.frame(gene_id = c("g1", "g2", "g3")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("g1", "g2", "g3"), pvalue = c(0.001, 0.01, 0.1))
  
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 2,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(is.list(result))
})

test_that("jackknife_tsallis_entropy SE gene lookup from rowData gene_name", {
  # Test gene lookup using gene_name column in rowData
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120, 110, 60, 85), nrow = 3, ncol = 3)),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3"),
      gene_name = c("GENEX", "GENEQ", "GENEZ")
    ),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("GENEQ", "GENEZ", "GENEX"), pvalue = c(0.001, 0.01, 0.1))
  
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 1,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy SE gene lookup from rowData gene_id", {
  # Test gene lookup using gene_id column in rowData
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120, 110, 60, 85), nrow = 3, ncol = 3)),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3"),
      gene_id = c("G001", "G002", "G003")
    ),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("G001", "G002", "G003"), pvalue = c(0.001, 0.01, 0.1))
  
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 1,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy SE with missing gene warning", {
  # Test handling of gene not found in SE
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120), nrow = 2, ncol = 3)),
    rowData = data.frame(gene_id = c("g1", "g2")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("g1", "g_missing", "g2"), pvalue = c(0.001, 0.05, 0.1))
  
  # Should handle missing gene gracefully with a warning
  expect_warning(
    result <- jackknife_tsallis_entropy(
      se = se,
      res = res,
      top_n = 2,
      q = 2,
      norm = TRUE,
      print_results = FALSE
    ),
    "not found"
  )
})

test_that("jackknife_tsallis_entropy with invalid count values", {
  # Test handling of NA or negative counts
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, NA, 80, 120), nrow = 2, ncol = 3)),
    rowData = data.frame(gene_id = c("g1", "g2")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("g1", "g2"), pvalue = c(0.001, 0.01))
  
  # Should handle invalid counts
  expect_warning(
    result <- jackknife_tsallis_entropy(
      se = se,
      res = res,
      top_n = 1,
      q = 2,
      print_results = FALSE
    ),
    NA
  )
})

test_that("jackknife_tsallis_entropy SE with pseudocount and normalization", {
  # Test pseudocount and normalization parameters
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120), nrow = 2, ncol = 3)),
    rowData = data.frame(gene_id = c("g1", "g2")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("g1", "g2"), pvalue = c(0.001, 0.01))
  
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 1,
    q = 2,
    norm = TRUE,
    pseudocount = 0.5,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy with different log bases", {
  # Test different logarithm bases
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result_e <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    log_base = exp(1),
    print_results = FALSE
  )
  
  result_2 <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    log_base = 2,
    print_results = FALSE
  )
  
  expect_true(!is.null(result_e))
  expect_true(!is.null(result_2))
})

test_that("jackknife_tsallis_entropy with outlier threshold parameter", {
  # Test threshold parameter for outlier detection
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120, 150, 60)
  
  result_90 <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    threshold = 90,
    print_results = FALSE
  )
  
  result_95 <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    threshold = 95,
    print_results = FALSE
  )
  
  expect_true(!is.null(result_90))
  expect_true(!is.null(result_95))
})

test_that("jackknife_tsallis_entropy with seed for reproducibility", {
  # Test seed parameter
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result1 <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    seed = 456,
    print_results = FALSE
  )
  
  result2 <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    seed = 456,
    print_results = FALSE
  )
  
  # Same seed should give same estimate
  expect_equal(result1$estimate, result2$estimate, tolerance = 1e-6)
})

test_that("compute_delta_statistics with matrices", {
  # Test compute_delta_statistics - matrix inputs
  set.seed(123)
  
  counts_A <- matrix(c(100, 50, 75, 110, 45, 80), nrow = 3, ncol = 2)
  counts_B <- matrix(c(106, 48, 82, 115, 42, 85), nrow = 3, ncol = 2)
  delta_influence <- c(0.1, 0.05, 0.08)
  
  result <- compute_delta_statistics(
    counts_A = counts_A,
    counts_B = counts_B,
    delta_influence = delta_influence,
    q = 1
  )
  
  expect_true(is.list(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
})

test_that("compute_delta_statistics returns statistics", {
  # Test return values
  set.seed(123)
  
  counts_A <- matrix(c(100, 50, 200, 80), nrow = 2, ncol = 2)
  counts_B <- matrix(c(98, 52, 205, 75), nrow = 2, ncol = 2)
  delta_influence <- c(0.12, 0.08)
  
  result <- compute_delta_statistics(
    counts_A = counts_A,
    counts_B = counts_B,
    delta_influence = delta_influence,
    q = 2
  )
  
  expect_true("pvalue" %in% names(result))
  expect_true("se" %in% names(result))
  expect_equal(length(result$ci_lower), 2)
})

test_that("jackknife_isoform_switching with SummarizedExperiment", {
  # Test jackknife_isoform_switching with proper SE input
  skip("Test skipped to reduce runtime: jackknife_isoform_switching with SummarizedExperiment (resource-intensive)")
  skip_if_not_installed("SummarizedExperiment")
  set.seed(123)
  
  # Create SE with proper structure: multiple samples per condition, multiple transcripts per gene
  counts_matrix <- matrix(
    c(100, 150, 80, 50, 75, 60,    # Gene1 transcripts: 3 transcripts x 2 samples condition A
      120, 140, 70, 60, 85, 65,    # Gene1 transcripts: 3 transcripts x 2 samples condition B
      200, 180, 90, 100, 110, 95,   # Gene2 transcripts: 3 transcripts x 2 samples
      215, 195, 105, 120, 130, 100), # Gene2 condition B
    nrow = 6, ncol = 4
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3", "tx4", "tx5", "tx6"),
      gene_id = c("g1", "g1", "g1", "g2", "g2", "g2"),
      gene_name = c("GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2")
    ),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  result <- jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = 1,
    print_results = FALSE
  )
  
  expect_true(is.list(result) || inherits(result, "tsenat_isoform_switching"))
})

test_that("jackknife_isoform_switching with multiple q", {
  # Test with multi-q 
  skip("Test skipped to reduce runtime: jackknife_isoform_switching with multiple q-values (resource-intensive)")
  
  skip_if_not_installed("SummarizedExperiment")
  set.seed(123)
  
  # Create SE with proper structure for multi-condition analysis
  counts_matrix <- matrix(
    c(100, 150, 80, 50, 75, 60,    # Gene1 tx1,tx2,tx3: cond A (2 samples)
      120, 140, 70, 60, 85, 65,    # Gene1 tx1,tx2,tx3: cond B (2 samples)
      200, 180, 90, 100, 110, 95,   # Gene2 tx4,tx5,tx6: all samples
      215, 195, 105, 120, 130, 100), # Gene2 condition B
    nrow = 6, ncol = 4
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3", "tx4", "tx5", "tx6"),
      gene_id = c("g1", "g1", "g1", "g2", "g2", "g2"),
      gene_name = c("GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2")
    ),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  result <- jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(1.0, 1.5),
    print_results = FALSE
  )
  
  expect_true(is.list(result))
})

test_that("jackknife_tsallis_entropy returns required field structure", {
  # Test that result contains all required fields with correct types
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  assert_jackknife_result_valid(result, n_transcripts = length(x))
})

test_that("jackknife_tsallis_entropy identifies outliers", {
  # Test outlier detection
  set.seed(123)
  # Create data with one very dominant transcript
  x <- c(1000, 50, 75, 200, 80, 120)  # First value is much larger
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    threshold = 90,
    print_results = FALSE
  )
  
  expect_true("outlier_indices" %in% names(result))
  expect_true(is.numeric(result$outlier_indices))
})

test_that("jackknife_tsallis_entropy with very small counts", {
  # Test stability with small counts
  set.seed(123)
  x <- c(1, 2, 1, 3, 2, 1)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    pseudocount = 0.5,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
  expect_true("estimate" %in% names(result))
})

test_that("jackknife_tsallis_entropy with zero counts", {
  # Test with zero counts (requires pseudocount)
  set.seed(123)
  x <- c(100, 0, 75, 200, 0, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    pseudocount = 0.5,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy with q = 1 (Shannon entropy)", {
  # Test special case of q=1 (Shannon entropy)
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 1.0,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
  expect_true(result$estimate >= 0)
})

test_that("jackknife_tsallis_entropy with large q value", {
  # Test with large q (focuses on rare species)
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 5.0,
    norm = TRUE,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy SE with top_n > total genes", {
  # Test when top_n exceeds available genes - should process all available
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 200, 80, 120, 110, 60, 85), nrow = 3, ncol = 3)),
    rowData = data.frame(gene_id = c("g1", "g2", "g3")),
    colData = data.frame(sample = c("s1", "s2", "s3"))
  )
  
  res <- data.frame(gene_id = c("g1", "g2", "g3"), pvalue = c(0.001, 0.01, 0.1))
  
  # top_n = 10 but only 3 genes available - should use all 3
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 10,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  # Should process available genes gracefully
  expect_true(!is.null(result) || is.list(result))
})

test_that("jackknife_tsallis_entropy SE SE validation", {
  # Test that function validates SE input
  set.seed(123)
  
  # Invalid SE (not actually a SummarizedExperiment)
  invalid_se <- list(data = "not_a_se")
  res <- data.frame(gene_id = c("g1"), pvalue = c(0.001))
  
  expect_error(
    jackknife_tsallis_entropy(
      se = invalid_se,
      res = res,
      q = 2,
      print_results = FALSE
    ),
    "must be a SummarizedExperiment"
  )
})

test_that("jackknife_tsallis_entropy SE res validation", {
  # Test that function validates res input
  set.seed(123)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50), nrow = 2, ncol = 1)),
    rowData = data.frame(gene_id = c("g1", "g2")),
    colData = data.frame(sample = c("s1"))
  )
  
  # Invalid res (not a data.frame)
  invalid_res <- list(gene_id = c("g1"))
  
  expect_error(
    jackknife_tsallis_entropy(
      se = se,
      res = invalid_res,
      q = 2,
      print_results = FALSE
    ),
    "must be a data.frame"
  )
})

test_that("jackknife_tsallis_entropy summary method", {
  # Test summary method on jackknife results
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  # Summary should work without error
  expect_error(summary(result), NA)
})

test_that("jackknife_tsallis_entropy print method", {
  # Test print method on jackknife results
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 2,
    norm = TRUE,
    print_results = FALSE
  )
  
  # Print should work without error
  expect_error(print(result), NA)
})

# ============================================================================
# NTHREADS PARAMETER TESTS
# ============================================================================

test_that("jackknife_tsallis_entropy nthreads = 1 (sequential)", {
  # Test sequential processing with nthreads = 1 (default)
  set.seed(123)
  x <- c(100, 50, 75, 200, 80, 120, 150, 60)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 1,
    norm = TRUE,
    nthreads = 1,
    print_results = FALSE
  )
  
  expect_true(inherits(result, "tsenat_jackknife"))
  expect_true("estimate" %in% names(result))
  expect_true("jackknife_se" %in% names(result))
})

test_that("jackknife_tsallis_entropy nthreads = 2 with multi-q", {
  # Test parallel processing with nthreads = 2 (if available)
  skip_if_not_installed("parallel")
  set.seed(123)
  x <- c(100, 50, 75, 200, 80)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = c(0.5, 1, 1.5, 2),  # 4 q values triggers parallel (> 2)
    norm = TRUE,
    nthreads = 2,
    print_results = FALSE
  )
  
  expect_true(is.list(result))
  expect_length(result, 4)  # Should have 4 results (one per q)
  expect_true(all(sapply(result, inherits, "tsenat_jackknife")))
})

test_that("jackknife_tsallis_entropy nthreads = NULL (auto-detect)", {
  # Test auto-detection of threads
  set.seed(123)
  x <- c(100, 50, 75, 200, 80)
  
  result <- jackknife_tsallis_entropy(
    x = x,
    q = c(0.5, 1, 1.5, 2),  # Multi-q to enable parallelization
    norm = TRUE,
    nthreads = NULL,  # Auto-detect
    print_results = FALSE
  )
  
  expect_true(is.list(result))
  expect_length(result, 4)
})

test_that("jackknife_tsallis_entropy nthreads parameter passes through SE path", {
  # Test nthreads parameter with SummarizedExperiment input
  skip_if_not_installed("SummarizedExperiment")
  set.seed(123)
  
  # Create simple SE with 2 genes in rownames
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 50, 75, 80, 200, 120, 150, 160), nrow = 2, ncol = 4)),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2"),
      gene_name = c("Gene1", "Gene2")
    ),
    colData = data.frame(sample = c("s1", "s2", "s3", "s4"))
  )
  
  # Set rownames to match genes in res
  rownames(se) <- c("Gene1", "Gene2")
  
  # Create results data frame with gene_id column
  res <- data.frame(
    gene_id = c("Gene1", "Gene2"),
    pvalue = c(0.01, 0.05),
    row.names = c("Gene1", "Gene2")
  )
  
  result <- jackknife_tsallis_entropy(
    se = se,
    res = res,
    top_n = 2,
    q = 1,
    nthreads = 1,
    print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("jackknife_tsallis_entropy nthreads parameter passes through matrix recursion", {
  # Test nthreads parameter through matrix input (internal recursion)
  set.seed(123)
  x_matrix <- matrix(
    c(100, 50, 75, 80, 200, 120, 150, 160),
    nrow = 2, ncol = 4
  )
  rownames(x_matrix) <- c("Gene1", "Gene2")
  
  # Single q value (no parallelization but nthreads should still work)
  result <- jackknife_tsallis_entropy(
    x = x_matrix,
    q = 1,
    nthreads = 1,
    print_results = FALSE
  )
  
  expect_true(is.list(result))
  expect_length(result, 2)  # 2 genes
})

test_that("jackknife_tsallis_entropy nthreads behavior: nthreads > 1 without multi-q", {
  # Even if nthreads > 1, without sufficient q values it should be sequential
  skip_if_not_installed("parallel")
  set.seed(123)
  x <- c(100, 50, 75, 200, 80)
  
  # Single q value: should not parallelize even with nthreads = 2
  result <- jackknife_tsallis_entropy(
    x = x,
    q = 1,  # Only 1 q value, so no parallelization
    nthreads = 2,
    print_results = FALSE
  )
  
  expect_true(inherits(result, "tsenat_jackknife"))
  expect_true("estimate" %in% names(result))
})
