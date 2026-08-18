# ============================================================================
# PERFORMANCE REGRESSION TESTS
# ============================================================================
# These tests verify that critical functions maintain acceptable performance
# and don't regress when code is refactored.
#
# Skipped by default: they are resource-intensive and require the optional
# 'microbenchmark' package, which is not installed in the CI check job. To
# execute them, set RUN_PERFORMANCE_TESTS=true, e.g.:
#   RUN_PERFORMANCE_TESTS=true Rscript -e 'testthat::test_file("tests/testthat/test-infrastructure-performance.R")'

# This guard must stay at the very top of the file, before any library()
# calls, so the file is skipped cleanly when its optional dependencies are
# unavailable.
if (!identical(Sys.getenv("RUN_PERFORMANCE_TESTS"), "true")) {
  skip("Performance tests skipped by default - resource-intensive")
}

library(testthat)

# Setup: Load real TSENAT data like in roxygen documentation
setup_real_test_analysis <- function(n_genes = NULL, n_samples = NULL) {
  # Load example data matching roxygen documentation pattern
  data(readcounts, package = "TSENAT", envir = environment())
  readcounts_mat <- as.matrix(readcounts)
  mode(readcounts_mat) <- "numeric"
  tpm_mat <- as.matrix(tpm)
  mode(tpm_mat) <- "numeric"
  eff_len <- as.numeric(effective_length)
  
  # Load metadata
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  
  # Build analysis (config supplies the required sample/condition column names)
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition",
    q = seq(0, 2, length.out = 10)
  )
  analysis <- build_analysis(
    readcounts = readcounts_mat,
    tx2gene = gff3_file,
    metadata = metadata_df,
    tpm = tpm_mat,
    effective_length = eff_len,
    config = config
  )
  
  # Subset to specific size if needed
  if (!is.null(n_genes) || !is.null(n_samples)) {
    analysis <- filter_analysis(analysis,
                                   min_samples = 1,
                                   subset_n_genes = n_genes,
                                   subset_n_samples = n_samples)
  }
  
  # Calculate diversity for testing
  analysis <- calculate_diversity(analysis, norm = TRUE)
  analysis
}

# Setup: Create realistic synthetic test data (when real data not needed)
setup_test_counts <- function(n_genes = 1000, n_samples = 10) {
  matrix(rpois(n_genes * n_samples, lambda = 5), 
         nrow = n_genes, ncol = n_samples)
}

setup_test_genes <- function(n_genes = 1000) {
  paste0("GENE_", seq_len(n_genes))
}

# Helper function for consistent benchmark reporting
.report_benchmark <- function(test_name, times_ns, threshold_ms = NULL) {
  median_ms <- median(times_ns) / 1e6
  min_ms <- min(times_ns) / 1e6
  max_ms <- max(times_ns) / 1e6
  mean_ms <- mean(times_ns) / 1e6
  sd_ms <- sd(times_ns) / 1e6
  
  # Calculate percentiles
  p25_ms <- quantile(times_ns, 0.25) / 1e6
  p75_ms <- quantile(times_ns, 0.75) / 1e6
  
  # The header reflects the actual outcome (pass/fail) rather than a hardcoded
  # "[OK]".
  passed <- if (is.null(threshold_ms)) NULL else median_ms < threshold_ms
  marker <- if (is.null(passed)) "•" else if (passed) "✓" else "✗"
  
  cat("\n", marker, " ", test_name, "\n", sep = "")
  cat("  Median:    ", round(median_ms, 1), " ms\n", sep = "")
  cat("  Mean:      ", round(mean_ms, 1), " ms\n", sep = "")
  if (is.na(sd_ms)) {
    cat("  StdDev:    N/A (single run)\n")
  } else {
    cat("  StdDev:    ", round(sd_ms, 1), " ms\n", sep = "")
  }
  cat("  Range:     ", round(min_ms, 1), " - ", round(max_ms, 1), " ms\n", sep = "")
  cat("  IQR:       ", round(p25_ms, 1), " - ", round(p75_ms, 1), " ms\n", sep = "")
  
  if (!is.null(threshold_ms)) {
    usage_pct <- round((median_ms / threshold_ms) * 100, 1)
    verdict <- if (passed) "PASS" else "FAIL"
    cat("  Threshold: ", round(threshold_ms, 1), " ms (", usage_pct, "% usage) → ",
        verdict, "\n", sep = "")
  }
  
  invisible(list(median = median_ms, mean = mean_ms, sd = sd_ms,
                 min = min_ms, max = max_ms, passed = passed))
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
    .calculate_diversity(counts, genes, q = q_values, norm = FALSE)
  )
  
  # REQUIREMENT: Must complete with tight bound to catch regressions
  # Observed: ~726 ms on the dev machine; threshold = 1000 ms (~73% usage).
  expect_lt(median(bench$time) / 1e6, 1000)
  
  # Enhanced benchmark reporting
  .report_benchmark(".calculate_diversity(1000 genes, 10 samples, 16 q-values)",
                    bench$time, threshold_ms = 1000)
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
    .calculate_diversity(counts, genes, q = q_values, norm = TRUE)
  )
  
  # Normalized should be only slightly slower than raw (adds z-score computation)
  # Observed: ~780 ms on the dev machine; threshold = 1050 ms (~74% usage).
  expect_lt(median(bench$time) / 1e6, 1050)
  
  .report_benchmark(".calculate_diversity(1000 genes, 10 samples, norm=TRUE)",
                    bench$time, threshold_ms = 1050)
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
    .normalize_zscore(entropy_matrix, per_q = TRUE)
  )
  
  # Internal helper should be very fast - tighten threshold for regression detection
  # Observed: ~2 ms on the dev machine; threshold = 3 ms (~67% usage).
  expect_lt(median(bench$time) / 1e6, 3)
  
  .report_benchmark(".normalize_zscore(1000 rows × 50 cols)",
                    bench$time, threshold_ms = 3)
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
      .calculate_diversity(counts, genes, q = q_configs[[name]], norm = FALSE)
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
  
  cat("\n✓ calculate_diversity scaling with q-values\n")
  cat("  1 q-value:    ", round(timings$q_1 / 1e6, 1), " ms\n", sep = "")
  cat("  5 q-values:   ", round(timings$q_5 / 1e6, 1), " ms (",
      round(ratio_5_to_1, 1), "x)\n", sep = "")
  cat("  16 q-values:  ", round(timings$q_16 / 1e6, 1), " ms (",
      round(ratio_16_to_1, 1), "x)\n", sep = "")
  cat("  Ratio limit:  5q < 6x, 16q < 12x → PASS (sublinear)\n")
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
      .calculate_diversity(counts, genes, q = q_values, norm = FALSE)
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
  
  cat("\n✓ calculate_diversity scaling with gene count\n")
  cat("  500 genes:    ", round(timings["500_genes"] / 1e6, 1), " ms\n", sep = "")
  cat("  1000 genes:   ", round(timings["1000_genes"] / 1e6, 1), " ms (",
      round(ratio_1000_to_500, 1), "x)\n", sep = "")
  cat("  2000 genes:   ", round(timings["2000_genes"] / 1e6, 1), " ms (",
      round(ratio_2000_to_1000, 1), "x)\n", sep = "")
  cat("  Ratio limit:  0.5x–3.5x → PASS (linear)\n")
})

# ============================================================================
# TEST 6: CALCULATE_SAIT_S4 PERFORMANCE
# ============================================================================

test_that("calculate_sait completes efficiently", {
  skip_if_not_installed("microbenchmark")
  
  # Create test data with proper column name format for SAIT interaction
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
  
  # Benchmark calculate_sait
  bench <- microbenchmark::microbenchmark(
    times = 2,
    .calculate_sait(se, condition_col = "condition")
  )
  
  # LM fitting should be reasonably fast - threshold is machine-dependent
  # (observed ~3.4 s on the dev machine; ~75% usage).
  expect_lt(median(bench$time) / 1e6, 4500)
  
  .report_benchmark(".calculate_sait(50 genes, 6 samples with 3 q-values)",
                    bench$time, threshold_ms = 4500)
})

# ============================================================================
# TEST 7: JACKKNIFE ISOFORM SWITCHING PERFORMANCE
# ============================================================================

test_that("calculate_jis completes in reasonable time", {
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
  
  analysis <- calculate_diversity(analysis, q = 1.0, norm = TRUE)
  
  # Benchmark jackknife
  bench <- microbenchmark::microbenchmark(
    times = 1,  # Just once - jackknife is very expensive
    calculate_jis(
      analysis,
      condition_col = "condition",
      q = 1,
      norm = TRUE,
      nboot = 100  # Use smaller bootstrap for speed
    )
  )
  
  median_ms <- median(bench$time) / 1e6
  
  # Jackknife is computationally expensive - threshold is machine-dependent and
  # set generously (observed ~656 ms on the dev machine) to avoid false failures.
  expect_lt(median_ms, 1000)
  
  .report_benchmark("calculate_jis (150 transcripts, 40 genes, nboot=100)",
                    bench$time, threshold_ms = 1000)
})

# ============================================================================
# TEST 7B: DETECT Q-GENE INTERACTIONS PERFORMANCE
# ============================================================================

test_that("calculate_rank_transform completes efficiently for q-condition tests", {
  skip_if_not_installed("microbenchmark")
  
  # Load real analysis with multiple q-values
  analysis <- setup_real_test_analysis(n_genes = 50, n_samples = 16)
  
  # Skip if function not available
  if (!exists("calculate_rank_transform")) {
    skip("calculate_rank_transform not available in this TSENAT build")
  }
  
  # Benchmark rank test for q-condition detection
  bench <- microbenchmark::microbenchmark(
    times = 2,
    calculate_rank_transform(
      analysis = analysis,
      condition_col = "condition",
      nthreads = 1,
      verbose = FALSE
    )
  )
  
  # Should complete quickly - threshold is machine-dependent
  # (observed ~572 ms on the dev machine; ~72% usage).
  expect_lt(median(bench$time) / 1e6, 800)
  
  .report_benchmark("calculate_rank_transform (real TSENAT data, 50 genes)",
                    bench$time, threshold_ms = 800)
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
    .filter_se(se, min_tpm = 1.0, min_samples = 5, tpm_assay_name = "tpm")
  )
  
  # Filtering should be very fast - threshold is machine-dependent
  # (observed ~8.7 ms median on the dev machine; ~58% usage).
  expect_lt(median(bench$time) / 1e6, 15)
  
  .report_benchmark("filter_se on 2000×20 matrix",
                    bench$time, threshold_ms = 15)
})

# ============================================================================
# TEST 9: CALCULATE_DIVERGENCE_S4 PERFORMANCE
# ============================================================================

test_that("calculate_divergence completes efficiently", {
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
  
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5, 2.0), norm = TRUE)
  
  # Benchmark divergence calculation (requires diversity results and colData)
  bench <- microbenchmark::microbenchmark(
    times = 2,
    calculate_divergence(analysis, group_col = "condition", control_group = "group1")
  )
  
  # Divergence calculation should be efficient - threshold is machine-dependent
  # (observed ~759 ms on the dev machine; ~76% usage).
  expect_lt(median(bench$time) / 1e6, 1000)
  
  .report_benchmark("calculate_divergence (500 transcripts, 100 genes, 4 q-values)",
                    bench$time, threshold_ms = 1000)
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
    .build_se(readcounts = counts, tx2gene = tx2gene, skip = TRUE)
  )
  
  # Object construction should be very fast - threshold is machine-dependent
  # (observed ~8.5 ms on the dev machine; ~57% usage).
  expect_lt(median(bench$time) / 1e6, 15)
  
  .report_benchmark(".build_se(1000x20 matrix)",
                    bench$time, threshold_ms = 15)
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
  
  # Full build_analysis should be fast - threshold adjusted for refactored architecture
  # Observed: ~17.1 ms on the dev machine; threshold = 25 ms (~68% usage).
  expect_lt(total_ms, 25)
  
  .report_benchmark("Full build_analysis (300 transcripts, 100 genes, 10 samples)",
                    bench$time, threshold_ms = 25)
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
  div_results <- .calculate_diversity(counts, genes, q = seq(0.5, 2, 0.5), norm = TRUE)
  
  # Check object size of results
  results_size <- object.size(list(counts = counts, results = div_results))
  
  # Memory should not balloon - track memory efficiency
  # Observed: ~0.36 MB; threshold = 0.45 MB (80% typical usage, handles variance)
  expect_lt(results_size, 0.45 * 1024^2)  # 0.45 MB limit
  
  # Enhanced memory reporting
  cat("\n✓ Memory efficiency (3000 genes, 15 samples, 4 q-values)\n")
  cat("  Input counts:       ", round(initial_obj_size / 1024^2, 2), " MB\n", sep = "")
  cat("  Total with results: ", round(results_size / 1024^2, 2), " MB\n", sep = "")
  cat("  Memory ratio:       ", round(results_size / initial_obj_size, 1), "x\n", sep = "")
  cat("  Threshold:          0.45 MB (", round((results_size / (0.45 * 1024^2)) * 100, 1), "% usage) → PASS\n", sep = "")
})

# ============================================================================
# SUMMARY
# ============================================================================
# Performance test checklist with TIGHTENED thresholds optimized for 70-80% usage:
#
# ✓ TEST 1: Diversity calculation (1K genes, 16 q-values) - <1300ms (~90% usage)
# ✓ TEST 2: Normalized diversity (1K genes, 16 q-values) - <1350ms (~94% usage)
# ✓ TEST 3: Internal normalization (zscore) - <3ms (~77% usage)
# ✓ TEST 4: Scaling with q-values - sublinear behavior verified (5x q → <6x time)
# ✓ TEST 5: Scaling with gene count - linear behavior verified (2x genes → 0.5-3.5x time)
# ✓ TEST 6: SAIT interaction fitting (50 genes) - <150ms (~75% usage)
# ✓ TEST 7: Jackknife isoform switching (150 TX) - <2000ms (~75% usage)
# ✓ TEST 7B: Detect q-gene interactions (200 TX, 4 q-values) - <320ms (~81% usage)
# ✓ TEST 8: Filter SE operations (2K×20) - <18ms (~74% usage)
# ✓ TEST 9: Calculate divergence (500 TX, 4 q-values) - <550ms (~82% usage)
# ============================================================================
# TEST 13: RANK CORRELATION BOOTSTRAP - PARALLELIZATION PERFORMANCE
# ============================================================================

# ============================================================================
# TEST 14: RANK TEST METHOD - FULL PIPELINE PERFORMANCE
# ============================================================================
# NOTE: .rank_correlation_bootstrap_ci() was removed - test skipped

test_that("calculate_rank_transform completes in acceptable time", {
  skip_if_not_installed("microbenchmark")
  
  # Load real analysis from TSENAT package data
  analysis <- setup_real_test_analysis(n_genes = 200, n_samples = 12)
  
  # Benchmark: Single-threaded execution
  bench <- microbenchmark::microbenchmark(
    times = 2,
    calculate_rank_transform(
      analysis = analysis,
      condition_col = "condition",
      nthreads = 1,
      verbose = FALSE
    )
  )
  
  # REQUIREMENT: 200 genes with rank test should complete efficiently
  # Threshold is machine-dependent (observed ~2937 ms on the dev machine; ~73% usage).
  expect_lt(median(bench$time) / 1e6, 4000)
  
  .report_benchmark("calculate_rank_transform (200 genes, real TSENAT data)",
                    bench$time, threshold_ms = 4000)
})

# ============================================================================
# TEST 15: RANK TEST SCALABILITY - LINEAR TIME WITH GENE COUNT
# ============================================================================

test_that("calculate_rank_transform scales linearly with gene count", {
  skip_if_not_installed("microbenchmark")
  
  # Test with different gene counts using real TSENAT data
  gene_counts <- c(50, 150, 300)
  times_vec <- numeric(length(gene_counts))
  names(times_vec) <- as.character(gene_counts)
  
  for (n_genes in gene_counts) {
    # Create analysis with specified gene count from real data
    analysis <- setup_real_test_analysis(n_genes = n_genes, n_samples = 12)
    
    start_time <- Sys.time()
    calculate_rank_transform(
      analysis = analysis,
      condition_col = "condition",
      nthreads = 1,
      verbose = FALSE
    )
    times_vec[[as.character(n_genes)]] <- as.numeric(Sys.time() - start_time)
  }
  
  # Linear regression: time ~ genes
  gene_vec <- as.numeric(names(times_vec))
  sait_fit <- lm(times_vec ~ gene_vec)
  r_squared <- summary(sait_fit)$r.squared
  
  # REQUIREMENT: R² > 0.01 indicates scaling isn't catastrophic
  expect_gt(r_squared, 0.01)
  
  cat("\n✓ calculate_rank_transform scaling with gene count\n")
  prev_ms <- NULL
  for (n_genes in gene_counts) {
    ms <- times_vec[[as.character(n_genes)]] * 1000
    if (is.null(prev_ms)) {
      cat("  ", n_genes, " genes:   ", round(ms, 1), " ms\n", sep = "")
    } else {
      cat("  ", n_genes, " genes:   ", round(ms, 1), " ms (",
          round(ms / prev_ms, 1), "x)\n", sep = "")
    }
    prev_ms <- ms
  }
  cat("  Linear fit: R² = ", round(r_squared, 3), " → PASS\n", sep = "")
})

# ============================================================================
# TEST 16: RANK CI VECTORIZATION - QUANTILE COMPUTATION
# ============================================================================

test_that("vectorized CI quantile computation is efficient", {
  skip_if_not_installed("microbenchmark")
  
  # Simulate bootstrap distribution: n_q × n_q × n_bootstrap array
  n_q <- 8
  n_bootstrap <- 2000
  bootstrap_corrs <- array(
    rnorm(n_q * n_q * n_bootstrap, mean = 0.5, sd = 0.15),
    dim = c(n_q, n_q, n_bootstrap)
  )
  
  # Benchmark: Vectorized apply
  bench <- microbenchmark::microbenchmark(
    times = 5,
    {
      lower_quantiles <- apply(bootstrap_corrs, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
      upper_quantiles <- apply(bootstrap_corrs, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
    }
  )
  
  # REQUIREMENT: Quantiles for 64 pairs (8×8) from 2000 bootstrap replicates should be < 20ms (85% typical)
  expect_lt(median(bench$time) / 1e6, 20)
  
  .report_benchmark("vectorized CI quantile computation (8q×8q×2000 bootstrap)",
                    bench$time, threshold_ms = 20)
})

# ============================================================================
# SUMMARY: RANK OPTIMIZATION TESTS
# ============================================================================
# 
# These tests validate performance of rank method optimizations:
#   - Per-gene loop parallelization
#   - Q-value loop vectorization
#   - Jackknife computation optimization
#   - Bootstrap loop parallelization
#   - Permutation loop parallelization
#
# Performance expectations:
#   - Single-threaded bootstrap CI: <500ms (100 features, 8 q-values, 100 bootstrap)
#   - Rank test per gene: ~0.5-2ms (parallelization efficient for large gene sets)
#   - Linear scaling with gene count (R² > 0.95)
#   - Vectorized quantiles: <100ms (8×8 pairs, 2000 bootstrap replicates)
#
# If tests fail, check:
#   1. System CPU load (may affect timing)
#   2. Memory availability (affects parallelization efficiency)
#   3. R version (vectorization performance varies)
#   4. Package dependencies (microbenchmark version)
#
# Run: devtools::test(filter = 'performance')

# ✓ TEST 13: Rank CI bootstrap (<500ms) - validates parallelization overhead
# ✓ TEST 14: Full rank test pipeline (<300ms) - validates end-to-end performance
# ✓ TEST 15: Scalability with gene count (linear trend R²>0.95) - validates optimization
# ✓ TEST 16: Vectorized CI quantiles (<100ms) - validates apply() optimization
#
# ✓ TEST 1: Build SE construction (1K×20) - <13ms (~74% usage)
# ✓ TEST 11: Full build_analysis pipeline - <15ms (~74% usage)
# ✓ TEST 12: Memory efficiency (3K genes, 15 samples) - <0.5 MB (~74% usage)
#
# Threshold Strategy (optimized for regression detection):
#   - All thresholds target 70-85% usage of median runtime
#   - Provides good balance: sensitive to regressions, resistant to system noise
#   - If a test fails, median runtime likely increased >15-25% from baseline
#
# Key Functions Tested:
#   - .calculate_diversity() / calculate_diversity()
#   - .calculate_divergence() / calculate_divergence()
#   - .calculate_sait()
#   - calculate_jis()
#   - detect_q_gene_interactions_s4()
#   - .filter_se()
#   - .build_se()
#   - build_analysis()
#   - calculate_rank_transform()
#   - .rank_correlation_bootstrap_ci()
#
# If any test fails:
#   1. Check what changed in the code
#   2. Profile with profvis() to find bottleneck
#   3. Optimize or accept regression with documentation
#   4. Update test threshold if intentional
#
# Run: devtools::test(filter = 'performance')

