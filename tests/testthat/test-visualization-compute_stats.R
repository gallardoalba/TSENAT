context("compute_stats: Data Transformation and Statistical Computation Helpers")

# ============================================================================
# TEST: Gene Selection Functions
# ============================================================================

testthat::test_that("select_top_genes returns correct number of genes", {
  # Create sample results data frame
  results <- data.frame(
    gene_id = c("G001", "G002", "G003", "G004", "G005"),
    padj = c(0.001, 0.01, 0.05, 0.1, 0.2),
    stringsAsFactors = FALSE
  )

  # Select top 3
  top_genes <- .select_top_genes(results, n_genes = 3)
  testthat::expect_equal(length(top_genes), 3)
  testthat::expect_equal(top_genes, c("G001", "G002", "G003"))
})

testthat::test_that("select_top_genes auto-detects gene column", {
  results <- data.frame(
    gene = c("A", "B", "C"),
    p_col = c(0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  top_genes <- .select_top_genes(results, p_col = "p_col", gene_col = "gene", n_genes = 2)
  testthat::expect_equal(length(top_genes), 2)
  testthat::expect_equal(top_genes, c("A", "B"))
})

testthat::test_that("select_top_genes handles empty results", {
  results <- data.frame(
    gene_id = character(0),
    padj = numeric(0),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .select_top_genes(results),
    "must be a non-empty data frame"
  )
})

testthat::test_that("select_top_genes requests more genes than available", {
  results <- data.frame(
    gene_id = c("G001", "G002"),
    padj = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )

  top_genes <- .select_top_genes(results, n_genes = 5)
  testthat::expect_equal(length(top_genes), 2)
})

# ============================================================================
# TEST: Gene Filtering by P-Value
# ============================================================================

testthat::test_that("filter_genes_by_pvalue returns significant genes", {
  results <- data.frame(
    gene_id = c("G001", "G002", "G003", "G004"),
    padj = c(0.001, 0.01, 0.1, 0.5),
    stringsAsFactors = FALSE
  )

  sig_genes <- .filter_genes_by_pvalue(results, p_threshold = 0.05)
  testthat::expect_equal(length(sig_genes), 2)
  testthat::expect_equal(sig_genes, c("G001", "G002"))
})

testthat::test_that("filter_genes_by_pvalue returns empty when no significant genes", {
  results <- data.frame(
    gene_id = c("G001", "G002"),
    padj = c(0.1, 0.5),
    stringsAsFactors = FALSE
  )

  sig_genes <- .filter_genes_by_pvalue(results, p_threshold = 0.05)
  testthat::expect_equal(length(sig_genes), 0)
})

testthat::test_that("filter_genes_by_pvalue auto-detects columns", {
  results <- data.frame(
    gene = c("A", "B", "C"),
    pvalue = c(0.001, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  sig_genes <- .filter_genes_by_pvalue(results, p_threshold = 0.06)
  testthat::expect_equal(length(sig_genes), 2)
})

# ============================================================================
# TEST: Validation Functions
# ============================================================================

testthat::test_that("validate_diversity_se checks for SummarizedExperiment class", {
  not_se <- data.frame(x = 1:10)

  testthat::expect_error(
    .validate_diversity_se(not_se),
    "must be a SummarizedExperiment"
  )
})

testthat::test_that("validate_diversity_se checks for diversity assay", {
  create_test_se_without_diversity <- function() {
    require_pkgs("SummarizedExperiment")
    mat <- matrix(1:10, nrow = 5)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
    se
  }

  se <- create_test_se_without_diversity()
  testthat::expect_error(
    .validate_diversity_se(se),
    "diversity.*assay not found"
  )
})

testthat::test_that("validate_diversity_se passes with valid diversity SE", {
  require_pkgs("SummarizedExperiment")

  mat <- matrix(rnorm(50), nrow = 5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = mat)
  )

  result <- .validate_diversity_se(se, check_metadata = FALSE)
  testthat::expect_true(result)
})

testthat::test_that("validate_results_df checks for gene column", {
  results <- data.frame(
    x = 1:5,
    padj = 0.05,
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .validate_results_df(results),
    "No gene identifier column found"
  )
})

testthat::test_that("validate_results_df checks for p-value column", {
  results <- data.frame(
    gene_id = paste0("G", 1:5),
    x = 1:5,
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .validate_results_df(results),
    "No p-value column found"
  )
})

testthat::test_that("validate_results_df passes with valid data frame", {
  results <- data.frame(
    gene_id = paste0("G", 1:5),
    padj = seq(0.001, 0.05, length.out = 5),
    stringsAsFactors = FALSE
  )

  result <- .validate_results_df(results)
  testthat::expect_true(result)
})

# ============================================================================
# TEST: Formatting Functions
# ============================================================================

testthat::test_that("format_pvalue handles different thresholds", {
  testthat::expect_equal(.format_pvalue(0.0001, threshold = 0.001), "< 0.001")
  testthat::expect_match(.format_pvalue(0.01, threshold = 0.001), "[0-9]")
})

testthat::test_that("format_pvalue handles NA values", {
  testthat::expect_equal(.format_pvalue(NA), "NA")
})

testthat::test_that("format_q_label formats q values correctly", {
  label <- .format_q_label(1.5)
  testthat::expect_match(label, "q = 1\\.50")
})

testthat::test_that("format_q_label handles NA", {
  testthat::expect_equal(.format_q_label(NA), "NA")
})

testthat::test_that("format_label removes underscores and capitalizes", {
  testthat::expect_equal(.format_label("fold_change"), "Fold change")
  testthat::expect_equal(.format_label("adjusted_p_values"), "Adjusted p values")
  testthat::expect_equal(.format_label("X"), "X")
})

testthat::test_that("format_label handles empty strings", {
  testthat::expect_equal(.format_label(""), "")
})

# ============================================================================
# TEST: Read tx2gene
# ============================================================================

testthat::test_that("read_tx2gene validates required columns", {
  bad_mapping <- data.frame(
    tx_id = c("tx1", "tx2"),
    gene_id = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .read_tx2gene(bad_mapping),
    "must have columns 'Transcript' and 'Gen'"
  )
})

testthat::test_that("read_tx2gene accepts valid data frame", {
  mapping <- data.frame(
    Transcript = c("tx1", "tx2", "tx3"),
    Gen = c("G1", "G1", "G2"),
    stringsAsFactors = FALSE
  )

  result <- .read_tx2gene(mapping)
  testthat::expect_equal(nrow(result), 3)
  testthat::expect_equal(colnames(result), c("Transcript", "Gen"))
})

# ============================================================================
# TEST: Infer Samples from Coldata
# ============================================================================

testthat::test_that("infer_samples_from_coldata handles row-indexed coldata", {
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  coldata <- data.frame(
    sample_type = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(coldata) <- colnames(counts)

  samples <- .infer_samples_from_coldata(coldata, counts, condition_col = "sample_type")
  testthat::expect_equal(samples, c("A", "B", "A", "B"))
})

testthat::test_that("infer_samples_from_coldata handles sample ID column", {
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  coldata <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4"),
    sample_type = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  samples <- .infer_samples_from_coldata(coldata, counts, condition_col = "sample_type")
  testthat::expect_equal(samples, c("A", "B", "A", "B"))
})

testthat::test_that("infer_samples_from_coldata errors on mismatched samples", {
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  coldata <- data.frame(
    sample_id = c("X1", "X2", "X3"),
    sample_type = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .infer_samples_from_coldata(coldata, counts, condition_col = "sample_type"),
    "doesn't match"
  )
})

# ============================================================================
# TEST: Create Aggregation Function
# ============================================================================

testthat::test_that("create_aggregation_function creates median function by default", {
  result <- .create_aggregation_function(metric = "median")

  testthat::expect_is(result$agg_fun, "function")
  testthat::expect_equal(result$metric_choice, "median")
  test_data <- c(1, 2, 3, 4, 5)
  expected_median <- median(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_median)
})

testthat::test_that("create_aggregation_function creates mean function", {
  result <- .create_aggregation_function(metric = "mean")

  testthat::expect_equal(result$metric_choice, "mean")
  test_data <- c(1, 2, 3, 4, 5)
  expected_mean <- mean(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_mean)
})

testthat::test_that("create_aggregation_function creates iqr function", {
  result <- .create_aggregation_function(metric = "iqr")

  testthat::expect_equal(result$metric_choice, "iqr")
  test_data <- c(1, 2, 3, 4, 5)
  expected_iqr <- IQR(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_iqr)
})

testthat::test_that("create_aggregation_function generates appropriate label", {
  result <- .create_aggregation_function(metric = "median")

  testthat::expect_match(result$agg_label_unique, "median")
})

# ============================================================================
# TEST: Build Transcript Long Format
# ============================================================================

testthat::test_that("build_transcript_long creates long-format data", {
  require_pkgs("tidyr")

  # Create test data
  mapping <- data.frame(
    Transcript = c("tx1", "tx2", "tx3"),
    Gen = c("G1", "G1", "G2"),
    stringsAsFactors = FALSE
  )

  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- mapping$Transcript
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  samples <- c("A", "B", "A", "B")

  result <- .build_transcript_long(
    gene_single = "G1",
    mapping = mapping,
    counts = counts,
    samples = samples
  )

  testthat::expect_is(result$df_long, "data.frame")
  testthat::expect_equal(ncol(result$df_long), 4)  # tx, sample, expr, group
  testthat::expect_equal(length(result$txs), 2)  # 2 transcripts for G1
})

testthat::test_that("build_transcript_long respects top_n parameter", {
  require_pkgs("tidyr")

  mapping <- data.frame(
    Transcript = c("tx1", "tx2", "tx3", "tx4"),
    Gen = c("G1", "G1", "G1", "G1"),
    stringsAsFactors = FALSE
  )

  counts <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(counts) <- mapping$Transcript
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  samples <- c("A", "B", "A", "B")

  result <- .build_transcript_long(
    gene_single = "G1",
    mapping = mapping,
    counts = counts,
    samples = samples,
    top_n = 2
  )

  testthat::expect_equal(length(result$txs), 2)
})

testthat::test_that("build_transcript_long errors on missing gene", {
  require_pkgs("tidyr")

  mapping <- data.frame(
    Transcript = c("tx1", "tx2"),
    Gen = c("G1", "G1"),
    stringsAsFactors = FALSE
  )

  counts <- matrix(1:8, nrow = 2, ncol = 4)
  rownames(counts) <- mapping$Transcript
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  samples <- c("A", "B", "A", "B")

  testthat::expect_error(
    .build_transcript_long(
      gene_single = "G999",
      mapping = mapping,
      counts = counts,
      samples = samples
    ),
    "No transcripts found"
  )
})

# ============================================================================
# TEST: Aggregate Transcript Data
# ============================================================================

testthat::test_that("aggregate_transcript_data computes log2 expression", {
  df_long <- data.frame(
    tx = c("tx1", "tx1", "tx2", "tx2"),
    group = c("A", "B", "A", "B"),
    expr = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )

  agg_fun <- function(x) median(x, na.rm = TRUE)
  result <- .aggregate_transcript_data(df_long, agg_fun, pseudocount = 0)

  testthat::expect_is(result, "data.frame")
  testthat::expect_true("log2expr" %in% colnames(result))
  # Check log2 calculation
  expected_log2 <- log2(c(10, 30, 20, 40) + 0)
  testthat::expect_equal(result$log2expr, expected_log2, tolerance = 1e-6)
})

testthat::test_that("aggregate_transcript_data applies pseudocount", {
  df_long <- data.frame(
    tx = c("tx1", "tx1"),
    group = c("A", "B"),
    expr = c(0, 0),
    stringsAsFactors = FALSE
  )

  agg_fun <- function(x) median(x, na.rm = TRUE)
  result <- .aggregate_transcript_data(df_long, agg_fun, pseudocount = 1)

  # log2(0 + 1) = 0
  expected_log2 <- log2(0 + 1)
  testthat::expect_equal(result$log2expr[1], expected_log2)
})

# ============================================================================
# TEST: Select Genes from Results
# ============================================================================

testthat::test_that("select_genes_from_results orders by p-value", {
  res <- data.frame(
    genes = c("G1", "G2", "G3"),
    padj = c(0.1, 0.001, 0.05),
    stringsAsFactors = FALSE
  )

  top <- .select_genes_from_results(res, top_n = 2)
  testthat::expect_equal(top, c("G2", "G3"))
})

testthat::test_that("select_genes_from_results removes duplicates", {
  res <- data.frame(
    genes = c("G1", "G1", "G2"),
    padj = c(0.001, 0.01, 0.05),
    stringsAsFactors = FALSE
  )

  top <- .select_genes_from_results(res, top_n = 2)
  testthat::expect_equal(length(top), 2)
  testthat::expect_equal(top[1], "G1")
})

testthat::test_that("select_genes_from_results errors on missing genes column", {
  res <- data.frame(
    gene_id = c("G1", "G2"),
    padj = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    .select_genes_from_results(res, top_n = 1),
    "must contain a 'genes' column"
  )
})

# ============================================================================
# TEST: Plot Tsallis Q-Curve Helpers
# ============================================================================

testthat::test_that(".prepare_combined_se converts TSENATAnalysis to SE", {
  # Create test TSENATAnalysis object
  analysis <- TSENAT:::.create_test_analysis(
    n_genes = 5, n_samples_per_group = 3,
    q_values = c(1, 2, 3), seed = 42
  )
  
  # Call helper directly
  se_combined <- TSENAT:::.prepare_combined_se(analysis)
  
  # Verify it's a SummarizedExperiment
  testthat::expect_true(methods::is(se_combined, "SummarizedExperiment"))
  
  # Verify it has diversity assay
  testthat::expect_true("diversity" %in% SummarizedExperiment::assayNames(se_combined))
  
  # Verify dimensions: rows = genes, cols = samples_per_group * 2_groups * q_values
  testthat::expect_equal(nrow(se_combined), 5)  # 5 genes
  # 3 samples per group × 2 groups × 3 q_values = 18 columns
  testthat::expect_equal(ncol(se_combined), 18)
  
  # Verify column names have q suffixes
  colnames_combined <- colnames(se_combined)
  testthat::expect_true(any(grepl("_q=", colnames_combined)))
})

testthat::test_that(".prepare_combined_se creates colData with q column", {
  analysis <- TSENAT:::.create_test_analysis(
    n_genes = 3, n_samples_per_group = 2,
    q_values = c(1, 2), seed = 42
  )
  
  se_combined <- TSENAT:::.prepare_combined_se(analysis)
  
  # Check colData has q column
  coldata <- SummarizedExperiment::colData(se_combined)
  testthat::expect_true("q" %in% colnames(coldata))
  
  # Verify q values are present (2 q_values, each repeated for n_samples)
  q_vals <- unique(coldata$q)
  testthat::expect_equal(length(q_vals), 2)
  testthat::expect_true(all(q_vals %in% c(1, 2)))
})

testthat::test_that(".prepare_combined_se creates rowData with gene_id", {
  analysis <- TSENAT:::.create_test_analysis(
    n_genes = 4, n_samples_per_group = 2,
    q_values = c(1, 2), seed = 42
  )
  
  se_combined <- TSENAT:::.prepare_combined_se(analysis)
  
  # Check rowData has gene_id column
  rowdata <- SummarizedExperiment::rowData(se_combined)
  testthat::expect_true("gene_id" %in% colnames(rowdata))
  testthat::expect_equal(nrow(rowdata), 4)
})

testthat::test_that(".prepare_combined_se handles dimension mismatches", {
  # Create small test analysis with 2 q values
  analysis <- TSENAT:::.create_test_analysis(
    n_genes = 3, n_samples_per_group = 2,
    q_values = c(1, 2), seed = 42
  )
  
  # Should handle conversion without error
  se_combined <- TSENAT:::.prepare_combined_se(analysis)
  testthat::expect_true(methods::is(se_combined, "SummarizedExperiment"))
  
  # Verify structure is correct
  testthat::expect_equal(nrow(se_combined), 3)
  # 2 samples per group × 2 groups × 2 q_values = 8 columns
  testthat::expect_equal(ncol(se_combined), 8)
})

testthat::test_that(".compute_gene_group_stats computes median and sd", {
  # Create sample long-format data
  long_data <- data.frame(
    Gene = c("G1", "G1", "G1", "G1", "G1", "G1"),
    q = c(1, 1, 1, 2, 2, 2),
    group = c("A", "A", "A", "B", "B", "B"),
    tsallis = c(0.5, 0.6, 0.7, 1.0, 1.1, 1.2),
    stringsAsFactors = FALSE
  )
  
  stats <- TSENAT:::.compute_gene_group_stats(long_data)
  
  # Verify output structure
  testthat::expect_true(is.data.frame(stats))
  testthat::expect_equal(nrow(stats), 2)  # 2 groups
  testthat::expect_true("central" %in% colnames(stats))
  testthat::expect_true("spread" %in% colnames(stats))
  testthat::expect_true("qnum" %in% colnames(stats))
})

testthat::test_that(".compute_gene_group_stats calculates correct values", {
  long_data <- data.frame(
    Gene = c("G1", "G1", "G1"),
    q = c(1, 1, 1),
    group = c("A", "A", "A"),
    tsallis = c(1.0, 2.0, 3.0),  # median = 2.0, sd = 1.0
    stringsAsFactors = FALSE
  )
  
  stats <- TSENAT:::.compute_gene_group_stats(long_data)
  
  # Median should be 2.0
  testthat::expect_equal(stats$central[1], 2.0)
  
  # SD should be 1.0 (var([1,2,3]) = 1, sd = sqrt(1) = 1)
  testthat::expect_equal(stats$spread[1], 1.0)
})

testthat::test_that(".compute_gene_group_stats handles multiple groups and q values", {
  long_data <- data.frame(
    Gene = rep("G1", 12),
    q = rep(c(1, 2), each = 6),
    group = rep(c("A", "A", "A", "B", "B", "B"), 2),
    tsallis = rnorm(12, mean = 1, sd = 0.1),
    stringsAsFactors = FALSE
  )
  
  stats <- TSENAT:::.compute_gene_group_stats(long_data)
  
  # Should have 4 rows: 2 groups × 2 q_values (A-1, A-2, B-1, B-2)
  testthat::expect_equal(nrow(stats), 4)
  testthat::expect_equal(nrow(stats[stats$group == "A", ]), 2)
  testthat::expect_equal(nrow(stats[stats$group == "B", ]), 2)
})

testthat::test_that(".bootstrap_aggregate_ci aggregates CI bounds", {
  # Create SE with CI assays
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = matrix(rnorm(20), nrow = 5, ncol = 4),
      ci_lower = matrix(rnorm(20, mean = 0.5), nrow = 5, ncol = 4),
      ci_upper = matrix(rnorm(20, mean = 1.5), nrow = 5, ncol = 4)
    ),
    colData = data.frame(
      sample = c("S1", "S2", "S1", "S2"),
      group = c("A", "A", "B", "B"),
      q = c(1, 1, 1, 1)
    )
  )
  
  # Create long format data
  long <- data.frame(
    sample = c("S1", "S1", "S2", "S2"),
    group = c("A", "A", "B", "B"),
    q = c(1, 1, 1, 1),
    tsallis = c(0.5, 0.6, 1.0, 1.1),
    stringsAsFactors = FALSE
  )
  
  plot_df <- TSENAT:::.bootstrap_aggregate_ci(se, long)
  
  # Verify output structure
  testthat::expect_true(is.data.frame(plot_df))
  testthat::expect_true("ci_lower" %in% colnames(plot_df))
  testthat::expect_true("ci_upper" %in% colnames(plot_df))
  testthat::expect_true("median" %in% colnames(plot_df))
  testthat::expect_true("group" %in% colnames(plot_df))
})

testthat::test_that(".bootstrap_aggregate_ci handles multiple q values", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = matrix(rnorm(24), nrow = 4, ncol = 6),
      ci_lower = matrix(rnorm(24, mean = 0.5), nrow = 4, ncol = 6),
      ci_upper = matrix(rnorm(24, mean = 1.5), nrow = 4, ncol = 6)
    ),
    colData = data.frame(
      sample = c("S1", "S2", "S1", "S2", "S1", "S2"),
      group = c("A", "A", "B", "B", "A", "B"),
      q = c(1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
    )
  )
  
  long <- data.frame(
    sample = c("S1", "S2", "S1", "S2", "S1", "S2"),
    group = c("A", "A", "B", "B", "A", "B"),
    q = c(1.0, 1.0, 1.0, 1.0, 2.0, 2.0),
    tsallis = rnorm(6),
    stringsAsFactors = FALSE
  )
  
  plot_df <- TSENAT:::.bootstrap_aggregate_ci(se, long)
  
  # Should have rows for each group-q combination
  testthat::expect_true(nrow(plot_df) >= 2)
  testthat::expect_true(all(plot_df$ci_lower <= plot_df$ci_upper))
})

testthat::test_that(".bootstrap_aggregate_ci CI bounds are valid", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = matrix(seq(1, 20), nrow = 5, ncol = 4),
      ci_lower = matrix(seq(0.1, 2.0, length.out = 20), nrow = 5, ncol = 4),
      ci_upper = matrix(seq(0.2, 2.1, length.out = 20), nrow = 5, ncol = 4)
    ),
    colData = data.frame(
      sample = c("S1", "S2", "S1", "S2"),
      group = c("A", "A", "B", "B"),
      q = c(1.0, 1.0, 1.0, 1.0),
      stringsAsFactors = FALSE
    )
  )
  
  long <- data.frame(
    sample = c("S1", "S2", "S1", "S2"),
    group = c("A", "A", "B", "B"),
    q = c(1.0, 1.0, 1.0, 1.0),
    tsallis = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  
  plot_df <- TSENAT:::.bootstrap_aggregate_ci(se, long)
  
  # All ci_lower should be less than or equal to ci_upper
  testthat::expect_true(all(plot_df$ci_lower <= plot_df$ci_upper))
  
  # All values should be finite
  testthat::expect_true(all(is.finite(plot_df$ci_lower)))
  testthat::expect_true(all(is.finite(plot_df$ci_upper)))
})

# TEST: GAM Interaction Helpers
# ============================================================================

testthat::test_that(".prepare_sample_group_mapping builds correct mapping", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = matrix(1:12, nrow = 3, ncol = 4)),
    colData = S4Vectors::DataFrame(
      sample_type = c("Normal", "Normal", "Tumor", "Tumor"),
      row.names = c("S1_q=1", "S2_q=1", "S3_q=1", "S4_q=1")
    )
  )
  
  cdata <- SummarizedExperiment::colData(se)
  mapping <- TSENAT:::.prepare_sample_group_mapping(cdata, "sample_type")
  
  # Should have 4 unique samples
  testthat::expect_equal(length(mapping), 4)
  
  # Check mapping values (unname to compare just values)
  testthat::expect_equal(unname(mapping["S1"]), "Normal")
  testthat::expect_equal(unname(mapping["S2"]), "Normal")
  testthat::expect_equal(unname(mapping["S3"]), "Tumor")
  testthat::expect_equal(unname(mapping["S4"]), "Tumor")
})

testthat::test_that(".plot_gam_prepare_gene_data creates valid plot data", {
  # Create matrix with samples x q-values columns
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3", "gene4", "gene5")
  colnames(mat) <- c("S1_q=1", "S2_q=1", "S3_q=2", "S4_q=2")
  
  sample_to_group <- c(S1 = "Normal", S2 = "Normal", S3 = "Tumor", S4 = "Tumor")
  
  plot_df <- TSENAT:::.plot_gam_prepare_gene_data("gene1", mat, sample_to_group)
  
  # Should have 4 rows (one per column)
  testthat::expect_equal(nrow(plot_df), 4)
  
  # Should have correct columns
  testthat::expect_true(all(c("sample", "group", "q", "entropy") %in% colnames(plot_df)))
  
  # Check values
  testthat::expect_equal(plot_df$group, c("Normal", "Normal", "Tumor", "Tumor"))
  testthat::expect_equal(plot_df$q, c(1, 1, 2, 2))
})

testthat::test_that(".plot_gam_prepare_gene_data returns NULL for invalid gene", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3", "gene4", "gene5")
  colnames(mat) <- c("S1_q=1", "S2_q=1", "S3_q=2", "S4_q=2")
  
  sample_to_group <- c(S1 = "Normal", S2 = "Normal", S3 = "Tumor", S4 = "Tumor")
  
  # Query non-existent gene
  result <- TSENAT:::.plot_gam_prepare_gene_data("invalid_gene", mat, sample_to_group)
  
  testthat::expect_null(result)
})

testthat::test_that(".plot_gam_fit_group generates predictions", {
  # Create data with more q values and points per group for GAM fitting
  q_seq <- seq(0.5, 2, by = 0.25)  # More q values for GAM smoothing
  
  plot_df <- data.frame(
    sample = rep(c("S1", "S2", "S3", "S4"), times = length(q_seq)),
    group = rep(c("Normal", "Normal", "Tumor", "Tumor"), times = length(q_seq)),
    q = rep(q_seq, each = 4),
    entropy = rnorm(4 * length(q_seq), mean = 2, sd = 0.3),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.plot_gam_fit_group(plot_df)
  
  # Should not be NULL
  testthat::expect_false(is.null(result))
  
  # Should have three components
  testthat::expect_true(all(c("plot_data", "pred_data", "group_levels") %in% names(result)))
  
  # pred_data should have predictions (100 per group)
  testthat::expect_true(nrow(result$pred_data) >= 100)
  
  # pred_data should have required columns
  testthat::expect_true(all(c("group", "q", "entropy_fit", "se") %in% colnames(result$pred_data)))
  
  # group_levels should contain both groups
  testthat::expect_equal(sort(result$group_levels), c("Normal", "Tumor"))
})

testthat::test_that(".plot_gam_fit_group returns NULL for insufficient data", {
  # Only 1 group - insufficient for fitting
  plot_df <- data.frame(
    sample = c("S1", "S2"),
    group = c("Normal", "Normal"),
    q = c(1, 2),
    entropy = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.plot_gam_fit_group(plot_df)
  
  testthat::expect_null(result)
})

testthat::test_that(".plot_select_genes returns user-specified genes", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3", "gene4"),
    adj_p_interaction = c(0.001, 0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )
  
  genes_specified <- c("gene2", "gene4")
  selected <- TSENAT:::.plot_select_genes(lm_res, genes = genes_specified, n_top = 2)
  
  testthat::expect_equal(selected, genes_specified)
})

testthat::test_that(".plot_select_genes filters by significance", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3", "gene4"),
    adj_p_interaction = c(0.001, 0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )
  
  # Select top 2 with sig_alpha = 0.05 should get gene1, gene2, gene3
  selected <- TSENAT:::.plot_select_genes(lm_res, genes = NULL, n_top = 2, sig_alpha = 0.05)
  
  # Should return exactly 2 genes (top 2 by p-value)
  testthat::expect_equal(length(selected), 2)
  
  # Should be the most significant
  testthat::expect_true("gene1" %in% selected)
  testthat::expect_true("gene2" %in% selected)
})

testthat::test_that(".plot_select_genes returns NULL when no significant genes", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  
  selected <- TSENAT:::.plot_select_genes(lm_res, genes = NULL, n_top = 2, sig_alpha = 0.05)
  
  testthat::expect_null(selected)
})
