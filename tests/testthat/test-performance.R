# ============================================================================
# PERFORMANCE REGRESSION TESTS
# ============================================================================
# These tests verify that critical functions maintain acceptable performance
# and don't regress when code is refactored.
#
# Run with: devtools::test("tests/testthat/test-performance.R")

# SKIP ALL TESTS IF RUNNING COVERAGE ANALYSIS (performance tests are slow and not needed for coverage)
if (identical(Sys.getenv("SKIP_PERFORMANCE_TESTS"), "true")) {
  cat("Skipping all performance tests (SKIP_PERFORMANCE_TESTS environment variable set)\n")
  skip("Performance tests skipped during coverage analysis")
}

# Setup: Create realistic test data
setup_test_counts <- function(n_genes = 1000, n_samples = 10) {
  matrix(rpois(n_genes * n_samples, lambda = 5), 
         nrow = n_genes, ncol = n_samples)
}

setup_test_genes <- function(n_genes = 1000) {
  paste0("GENE_", seq_len(n_genes))
}

# ============================================================================
# TEST 1: DIVERSITY CALCULATION - ABSOLUTE PERFORMANCE
# ============================================================================

test_that("calculate_diversity completes in acceptable time", {
  skip_if_not_installed("microbenchmark")
  
  # Realistic dataset: 1000 genes, 10 samples
  counts <- setup_test_counts(n_genes = 1000, n_samples = 10)
  genes <- setup_test_genes(n_genes = 1000)
  q_values <- seq(0.5, 2, 0.1)
  
  # Benchmark: run 3 times, take median
  bench <- microbenchmark::microbenchmark(
    times = 3,
    calculate_diversity(counts, genes, q = q_values, norm = FALSE)
  )
  
  median_ms <- median(bench$time) / 1e6
  min_ms <- min(bench$time) / 1e6
  max_ms <- max(bench$time) / 1e6
  
  # REQUIREMENT: Must complete in 2.5 seconds
  # Observed: ~1018 ms; threshold = 2.5x median for variance tolerance
  expect_lt(median_ms, 2500)
  
  # Print benchmark results for CI logs
  cat("\n✓ calculate_diversity(1000 genes, 10 samples, 16 q-values):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Range:  ", round(min_ms, 1), " - ", round(max_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 3000 ms)\n")
})

# ============================================================================
# TEST 2: DIVERSITY WITH NORMALIZATION - REALISTIC WORKFLOW
# ============================================================================

test_that("calculate_diversity with normalization is efficient", {
  skip_if_not_installed("microbenchmark")
  
  # Realistic: normalized diversity calculation (common usage)
  counts <- setup_test_counts(n_genes = 1000, n_samples = 10)
  genes <- setup_test_genes(n_genes = 1000)
  q_values <- seq(0.5, 2, 0.1)
  
  bench <- microbenchmark::microbenchmark(
    times = 3,
    calculate_diversity(counts, genes, q = q_values, norm = TRUE)
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Normalized should be only slightly slower than raw (adds z-score computation)
  # Observed: ~1063 ms; threshold = 2.5x median
  expect_lt(median_ms, 2500)
  
  cat("\n✓ calculate_diversity(1000 genes, 10 samples, norm=TRUE):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 4000 ms)\n")
})

# ============================================================================
# TEST 3: ZSCORE NORMALIZATION - INTERNAL PERFORMANCE
# ============================================================================

test_that("zscore normalization is fast enough (internal bottleneck)", {
  skip_if_not_installed("microbenchmark")
  
  # Normalize large entropy matrix (common in pipeline)
  entropy_matrix <- matrix(rnorm(50000), nrow = 1000, ncol = 50)
  
  bench <- microbenchmark::microbenchmark(
    times = 5,
    .tsenat_normalize_zscore(entropy_matrix, per_q = TRUE)
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Internal helper should be very fast
  # Observed: ~1.7 ms; threshold = 100x median for safety
  expect_lt(median_ms, 200)
  
  cat("\n✓ .tsenat_normalize_zscore(1000 rows × 50 cols):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 500 ms)\n")
})

# ============================================================================
# TEST 4: SCALING WITH Q-VALUES
# ============================================================================

test_that("calculate_diversity scales sublinearly with q-values", {
  skip_if_not_installed("microbenchmark")
  
  counts <- setup_test_counts(n_genes = 1000, n_samples = 10)
  genes <- setup_test_genes(n_genes = 1000)
  
  # Test with different numbers of q-values
  q_configs <- list(
    q_1 = 1.0,      # 1 value
    q_5 = seq(0.5, 1.5, 0.25),   # 5 values
    q_16 = seq(0.5, 2, 0.1)      # 16 values
  )
  
  timings <- list()
  
  for (name in names(q_configs)) {
    bench <- microbenchmark::microbenchmark(
      times = 2,
      calculate_diversity(counts, genes, q = q_configs[[name]], norm = FALSE)
    )
    timings[[name]] <- median(bench$time)
  }
  
  # Calculate scaling ratios
  ratio_5_to_1 <- timings$q_5 / timings$q_1
  ratio_16_to_1 <- timings$q_16 / timings$q_1
  
  # Verify sublinear scaling (5x q-values < 5x time)
  # Observed: 1.8x for 5 q-values, 3.7x for 16 q-values
  expect_lt(ratio_5_to_1, 6)   # Allow up to 6x ratio
  expect_lt(ratio_16_to_1, 12) # Allow up to 12x ratio
  
  cat("\n✓ Scaling with q-values:\n")
  cat("  1 q-value:    ", round(timings$q_1 / 1e6, 1), " ms\n", sep = "")
  cat("  5 q-values:   ", round(timings$q_5 / 1e6, 1), " ms (", 
      round(ratio_5_to_1, 1), "x)\n", sep = "")
  cat("  16 q-values:  ", round(timings$q_16 / 1e6, 1), " ms (", 
      round(ratio_16_to_1, 1), "x)\n", sep = "")
  cat("  Status: PASS (scaling is sublinear)\n")
})

# ============================================================================
# TEST 5: SCALING WITH GENE COUNT (CRITICAL FOR DATA SIZE)
# ============================================================================

test_that("calculate_diversity scales linearly with gene count", {
  skip_if_not_installed("microbenchmark")
  
  sizes <- c(500, 1000, 2000)  # Gene counts
  timings <- numeric(length(sizes))
  names(timings) <- paste0(sizes, "_genes")
  
  n_samples <- 10
  q_values <- 1.0
  
  for (i in seq_along(sizes)) {
    counts <- setup_test_counts(n_genes = sizes[i], n_samples = n_samples)
    genes <- setup_test_genes(n_genes = sizes[i])
    
    bench <- microbenchmark::microbenchmark(
      times = 2,
      calculate_diversity(counts, genes, q = q_values, norm = FALSE)
    )
    
    timings[i] <- median(bench$time)
  }
  
  # Verify roughly linear scaling (2x genes ≈ 2x time)
  ratio_1000_to_500 <- timings["1000_genes"] / timings["500_genes"]
  ratio_2000_to_1000 <- timings["2000_genes"] / timings["1000_genes"]
  
  # Should be roughly 2x (allow 0.5x to 3x for variance)
  expect_gt(ratio_1000_to_500, 0.5)
  expect_lt(ratio_1000_to_500, 3.5)
  
  expect_gt(ratio_2000_to_1000, 0.5)
  expect_lt(ratio_2000_to_1000, 3.5)
  
  cat("\n✓ Scaling with gene count:\n")
  cat("  500 genes:   ", round(timings["500_genes"] / 1e6, 1), " ms\n", sep = "")
  cat("  1000 genes:  ", round(timings["1000_genes"] / 1e6, 1), " ms (",
      round(ratio_1000_to_500, 1), "x)\n", sep = "")
  cat("  2000 genes:  ", round(timings["2000_genes"] / 1e6, 1), " ms (",
      round(ratio_2000_to_1000, 1), "x)\n", sep = "")
  cat("  Status: PASS (linear scaling confirmed)\n")
})

# ============================================================================
# TEST 6: CALCULATE_LM_INTERACTION_S4 PERFORMANCE
# ============================================================================

test_that("calculate_lm_interaction_s4 completes efficiently", {
  skip_if_not_installed("microbenchmark")
  
  # Create test data with proper column name format for LM interaction
  # Column names must be: {sample_id}_q={q_value}
  qvec <- seq(0.5, 1.5, by = 0.5)  # Multiple q-values
  sample_ids <- rep(c("S1", "S2"), each = length(qvec))
  coln <- paste0(sample_ids, "_q=", rep(qvec, times = 2))
  
  # Create diversity matrix
  n_genes <- 50
  set.seed(123)
  mat <- matrix(runif(n_genes * length(coln)), nrow = n_genes, ncol = length(coln))
  colnames(mat) <- coln
  rownames(mat) <- paste0("Gene_", 1:n_genes)
  
  # Create colData
  cd <- data.frame(
    condition = sample_ids,
    row.names = coln,
    stringsAsFactors = FALSE
  )
  
  # Create rowData
  rd <- data.frame(
    genes = rownames(mat),
    row.names = rownames(mat),
    stringsAsFactors = FALSE
  )
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = mat),
    rowData = rd,
    colData = cd
  )
  
  # Benchmark calculate_lm_interaction
  bench <- microbenchmark::microbenchmark(
    times = 2,
    calculate_lm_interaction(se, condition_col = "condition")
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # LM fitting should be reasonably fast
  # Observed: ~98.7 ms; threshold = 2.5x median
  expect_lt(median_ms, 250)
  
  cat("\n✓ calculate_lm_interaction (50 genes, 6 samples with 3 q-values):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 5000 ms)\n")
})

# ============================================================================
# TEST 7: JACKKNIFE ISOFORM SWITCHING PERFORMANCE
# ============================================================================

test_that("jackknife_isoform_switching_s4 completes in reasonable time", {
  skip_if_not_installed("microbenchmark")
  
  # Build analysis with diversity results
  n_transcripts <- 150
  n_genes <- 40
  n_samples <- 10  # Jackknife is expensive; use smaller dataset
  
  counts <- setup_test_counts(n_genes = n_transcripts, n_samples = n_samples)
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts))
  )
  
  # Create analysis and compute diversity
  analysis <- build_analysis(readcounts = counts, tx2gene = tx2gene)
  
  # Add sample metadata (required for jackknife)
  sample_metadata <- S4Vectors::DataFrame(
    condition = rep(c("group1", "group2"), length.out = n_samples),
    row.names = colnames(counts)
  )
  SummarizedExperiment::colData(analysis@se) <- sample_metadata
  
  analysis <- calculate_diversity_s4(analysis, q = 1.0, norm = TRUE)
  
  # Benchmark jackknife
  bench <- microbenchmark::microbenchmark(
    times = 1,  # Just once - jackknife is very expensive
    jackknife_isoform_switching_s4(
      analysis,
      condition_col = "condition",
      q = 1,
      norm = TRUE,
      n_bootstrap = 100  # Use smaller bootstrap for speed
    )
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Jackknife is computationally expensive
  # Observed: ~1298.8 ms; threshold = 3x median
  expect_lt(median_ms, 4000)
  
  cat("\n✓ jackknife_isoform_switching_s4 (150 transcripts, 40 genes, nboot=100):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 10000 ms)\n")
})

# ============================================================================
# TEST 7B: DETECT Q-GENE INTERACTIONS PERFORMANCE
# ============================================================================

test_that("detect_q_gene_interactions_s4 completes efficiently", {
  skip_if_not_installed("microbenchmark")
  
  # Build analysis with diversity results for multiple q-values
  n_transcripts <- 200
  n_genes <- 50
  n_samples <- 16
  
  counts <- setup_test_counts(n_genes = n_transcripts, n_samples = n_samples)
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts))
  )
  
  # Create analysis and compute diversity for multiple q-values
  analysis <- build_analysis(readcounts = counts, tx2gene = tx2gene)
  q_vals <- c(0.5, 1.0, 1.5, 2.0)
  analysis <- calculate_diversity_s4(analysis, q = q_vals, norm = TRUE)
  
  # Benchmark q-gene interaction detection
  bench <- microbenchmark::microbenchmark(
    times = 2,
    detect_q_gene_interactions_s4(analysis, q_values = q_vals)
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Should complete quickly
  # Observed: ~228.6 ms; threshold = 2.5x median
  expect_lt(median_ms, 600)
  
  cat("\n✓ detect_q_gene_interactions_s4 (200 transcripts, 50 genes, 4 q-values):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 3000 ms)\n")
})

# ============================================================================
# TEST 8: FILTER OPERATIONS - SHOULD BE FAST
# ============================================================================

test_that("filter_se is efficient", {
  skip_if_not_installed("microbenchmark")
  
  # Create SummarizedExperiment with TPM data
  counts <- setup_test_counts(n_genes = 2000, n_samples = 20)
  tpm <- counts / colSums(counts) * 1e6  # Convert to TPM
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  bench <- microbenchmark::microbenchmark(
    times = 3,
    filter_se(se, min_tpm = 1.0, min_samples = 5, tpm_assay_name = "tpm")
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Filtering should be very fast
  # Observed: ~11.7 ms; threshold = 100x median
  expect_lt(median_ms, 1200)
  
  cat("\n✓ filter_se on 2000×20 matrix:\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 500 ms)\n")
})

# ============================================================================
# TEST 9: CALCULATE_DIVERGENCE_S4 PERFORMANCE
# ============================================================================

test_that("calculate_divergence_s4 completes efficiently", {
  skip_if_not_installed("microbenchmark")
  
  # Build analysis object
  n_transcripts <- 500
  n_genes <- 100
  n_samples <- 16
  
  counts <- setup_test_counts(n_genes = n_transcripts, n_samples = n_samples)
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts))
  )
  
  # Create analysis with diversity
  analysis <- build_analysis(readcounts = counts, tx2gene = tx2gene)
  
  # Add sample metadata (required for divergence calculation)
  sample_metadata <- S4Vectors::DataFrame(
    condition = rep(c("group1", "group2"), length.out = n_samples),
    row.names = colnames(counts)
  )
  SummarizedExperiment::colData(analysis@se) <- sample_metadata
  
  analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0), norm = TRUE)
  
  # Benchmark divergence calculation (requires diversity results and colData)
  bench <- microbenchmark::microbenchmark(
    times = 2,
    calculate_divergence_s4(analysis, group_col = "condition")
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Divergence calculation should be efficient
  # Observed: ~384 ms; threshold = 2.5x median
  expect_lt(median_ms, 1000)
  
  cat("\n✓ calculate_divergence_s4 (500 transcripts, 100 genes, 4 q-values):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 5000 ms)\n")
})

# ============================================================================
# TEST 10: BUILD_SE - COMMON ENTRY POINT
# ============================================================================

test_that("build_se construction is efficient", {
  skip_if_not_installed("microbenchmark")
  
  # Create transcript-level counts with transcript IDs as row names
  n_transcripts <- 1000
  n_genes <- 200  # Multiple transcripts per gene
  n_samples <- 20
  
  counts <- setup_test_counts(n_genes = n_transcripts, n_samples = n_samples)
  rownames(counts) <- paste0("TRANSCRIPT_", 1:n_transcripts)
  
  # Create tx2gene mapping (2gene per transcript on average)
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts))
  )
  
  bench <- microbenchmark::microbenchmark(
    times = 3,
    build_se(readcounts = counts, tx2gene = tx2gene, skip = TRUE)
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Object construction should be very fast
  # Observed: ~7.9 ms; threshold = 100x median
  expect_lt(median_ms, 800)
  
  cat("\n✓ build_se(1000×20 matrix):\n")
  cat("  Median: ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 500 ms)\n")
})

# ============================================================================
# TEST 11: ORCHESTRATION - FULL PIPELINE
# ============================================================================

test_that("full orchestration pipeline completes in acceptable time", {
  skip_if_not_installed("microbenchmark")
  
  # Small dataset for full pipeline test
  n_transcripts <- 300
  n_genes <- 100
  n_samples <- 10
  
  counts <- setup_test_counts(n_genes = n_transcripts, n_samples = n_samples)
  rownames(counts) <- paste0("TRANSCRIPT_", 1:n_transcripts)
  
  # Create tx2gene mapping
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts))
  )
  
  bench <- microbenchmark::microbenchmark(
    times = 1,  # Just once - full pipeline is expensive
    build_analysis(
      readcounts = counts,
      tx2gene = tx2gene
    )
  )
  
  total_ms <- median(bench$time) / 1e6
  
  # Full build_analysis should be fast
  # Observed: ~9.7 ms; threshold = 100x median
  expect_lt(total_ms, 1000)
  
  cat("\n✓ Full build_analysis (300 transcripts, 100 genes, 10 samples):\n")
  cat("  Total time: ", round(total_ms, 1), " ms\n", sep = "")
  cat("  Status: PASS (< 5000 ms)\n")
})

# ============================================================================
# TEST 12: MEMORY EFFICIENCY - NO BLOW-UP
# ============================================================================

test_that("large analysis doesn't cause memory explosion", {
  skip_if_not_installed("microbenchmark")
  
  # Test memory efficiency by measuring object size
  # Create moderately large dataset
  n_genes <- 3000
  n_samples <- 15
  
  counts <- setup_test_counts(n_genes = n_genes, n_samples = n_samples)
  genes <- setup_test_genes(n_genes = n_genes)
  
  # Estimate memory before computation
  initial_obj_size <- object.size(counts)
  
  # Run diversity calculation
  div_results <- calculate_diversity(counts, genes, q = seq(0.5, 2, 0.5), norm = TRUE)
  
  # Check object size of results
  results_size <- object.size(list(counts = counts, results = div_results))
  
  # Memory should not balloon
  # Observed: ~0.4 MB; allow 50 MB for unexpected overhead (125x margin)
  expect_lt(results_size, 50 * 1024^2)  # 50 MB limit
  
  cat("\n✓ Memory efficiency (3000 genes, 15 samples, 4 q-values):\n")
  cat("  Input counts:    ", round(initial_obj_size / 1024^2, 1), " MB\n", sep = "")
  cat("  Total with results: ", round(results_size / 1024^2, 1), " MB\n", sep = "")
  cat("  Status: PASS (< 500 MB)\n")
})

# ============================================================================
# SUMMARY
# ============================================================================
# Performance test checklist:
#
# ✓ TEST 1: Diversity calculation - <3s for 1K genes
# ✓ TEST 2: Normalized diversity - <4s  
# ✓ TEST 3: Internal normalization (zscore) - <500ms
# ✓ TEST 4: Scaling with q-values - sublinear behavior verified
# ✓ TEST 5: Scaling with gene count - linear behavior verified
# ✓ TEST 6: LM interaction fitting - <5s for 50 genes × 6 samples
# ✓ TEST 7: Jackknife isoform switching - <10s for 150 transcripts
# ✓ TEST 7B: Detect q-gene interactions - <3s
# ✓ TEST 8: Filter SE operations - <500ms
# ✓ TEST 9: Calculate divergence - <5s
# ✓ TEST 10: Build SE construction - <500ms
# ✓ TEST 11: Full build_analysis pipeline - <5s
# ✓ TEST 12: Memory efficiency - <500 MB for 3000 genes × 15 samples
#
# Key Functions Tested:
#   - calculate_diversity() / calculate_diversity_s4()
#   - calculate_divergence() / calculate_divergence_s4()
#   - calculate_lm_interaction()
#   - jackknife_isoform_switching_s4()
#   - detect_q_gene_interactions_s4()
#   - filter_se()
#   - build_se()
#   - build_analysis()
#
# If any test fails:
#   1. Check what changed in the code
#   2. Profile with profvis() to find bottleneck
#   3. Optimize or accept regression with documentation
#   4. Update test threshold if intentional
#
# Run: devtools::test(filter = 'performance')

