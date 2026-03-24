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
  top_genes <- select_top_genes(results, n_genes = 3)
  testthat::expect_equal(length(top_genes), 3)
  testthat::expect_equal(top_genes, c("G001", "G002", "G003"))
})

testthat::test_that("select_top_genes auto-detects gene column", {
  results <- data.frame(
    gene = c("A", "B", "C"),
    p_col = c(0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  top_genes <- select_top_genes(results, p_col = "p_col", gene_col = "gene", n_genes = 2)
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
    select_top_genes(results),
    "must be a non-empty data frame"
  )
})

testthat::test_that("select_top_genes requests more genes than available", {
  results <- data.frame(
    gene_id = c("G001", "G002"),
    padj = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )

  top_genes <- select_top_genes(results, n_genes = 5)
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

  sig_genes <- filter_genes_by_pvalue(results, p_threshold = 0.05)
  testthat::expect_equal(length(sig_genes), 2)
  testthat::expect_equal(sig_genes, c("G001", "G002"))
})

testthat::test_that("filter_genes_by_pvalue returns empty when no significant genes", {
  results <- data.frame(
    gene_id = c("G001", "G002"),
    padj = c(0.1, 0.5),
    stringsAsFactors = FALSE
  )

  sig_genes <- filter_genes_by_pvalue(results, p_threshold = 0.05)
  testthat::expect_equal(length(sig_genes), 0)
})

testthat::test_that("filter_genes_by_pvalue auto-detects columns", {
  results <- data.frame(
    gene = c("A", "B", "C"),
    pvalue = c(0.001, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  sig_genes <- filter_genes_by_pvalue(results, p_threshold = 0.06)
  testthat::expect_equal(length(sig_genes), 2)
})

# ============================================================================
# TEST: Validation Functions
# ============================================================================

testthat::test_that("validate_diversity_se checks for SummarizedExperiment class", {
  not_se <- data.frame(x = 1:10)

  testthat::expect_error(
    validate_diversity_se(not_se),
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
    validate_diversity_se(se),
    "diversity.*assay not found"
  )
})

testthat::test_that("validate_diversity_se passes with valid diversity SE", {
  require_pkgs("SummarizedExperiment")

  mat <- matrix(rnorm(50), nrow = 5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = mat)
  )

  result <- validate_diversity_se(se, check_metadata = FALSE)
  testthat::expect_true(result)
})

testthat::test_that("validate_results_df checks for gene column", {
  results <- data.frame(
    x = 1:5,
    padj = 0.05,
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    validate_results_df(results),
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
    validate_results_df(results),
    "No p-value column found"
  )
})

testthat::test_that("validate_results_df passes with valid data frame", {
  results <- data.frame(
    gene_id = paste0("G", 1:5),
    padj = seq(0.001, 0.05, length.out = 5),
    stringsAsFactors = FALSE
  )

  result <- validate_results_df(results)
  testthat::expect_true(result)
})

# ============================================================================
# TEST: Formatting Functions
# ============================================================================

testthat::test_that("format_pvalue handles different thresholds", {
  testthat::expect_equal(format_pvalue(0.0001, threshold = 0.001), "< 0.001")
  testthat::expect_match(format_pvalue(0.01, threshold = 0.001), "[0-9]")
})

testthat::test_that("format_pvalue handles NA values", {
  testthat::expect_equal(format_pvalue(NA), "NA")
})

testthat::test_that("format_q_label formats q values correctly", {
  label <- format_q_label(1.5)
  testthat::expect_match(label, "q = 1\\.50")
})

testthat::test_that("format_q_label handles NA", {
  testthat::expect_equal(format_q_label(NA), "NA")
})

testthat::test_that("format_label removes underscores and capitalizes", {
  testthat::expect_equal(format_label("fold_change"), "Fold change")
  testthat::expect_equal(format_label("adjusted_p_values"), "Adjusted p values")
  testthat::expect_equal(format_label("X"), "X")
})

testthat::test_that("format_label handles empty strings", {
  testthat::expect_equal(format_label(""), "")
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
    read_tx2gene(bad_mapping),
    "must have columns 'Transcript' and 'Gen'"
  )
})

testthat::test_that("read_tx2gene accepts valid data frame", {
  mapping <- data.frame(
    Transcript = c("tx1", "tx2", "tx3"),
    Gen = c("G1", "G1", "G2"),
    stringsAsFactors = FALSE
  )

  result <- read_tx2gene(mapping)
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

  samples <- infer_samples_from_coldata(coldata, counts, condition_col = "sample_type")
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

  samples <- infer_samples_from_coldata(coldata, counts, condition_col = "sample_type")
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
    infer_samples_from_coldata(coldata, counts, condition_col = "sample_type"),
    "doesn't match"
  )
})

# ============================================================================
# TEST: Create Aggregation Function
# ============================================================================

testthat::test_that("create_aggregation_function creates median function by default", {
  result <- create_aggregation_function(metric = "median")

  testthat::expect_is(result$agg_fun, "function")
  testthat::expect_equal(result$metric_choice, "median")
  test_data <- c(1, 2, 3, 4, 5)
  expected_median <- median(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_median)
})

testthat::test_that("create_aggregation_function creates mean function", {
  result <- create_aggregation_function(metric = "mean")

  testthat::expect_equal(result$metric_choice, "mean")
  test_data <- c(1, 2, 3, 4, 5)
  expected_mean <- mean(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_mean)
})

testthat::test_that("create_aggregation_function creates iqr function", {
  result <- create_aggregation_function(metric = "iqr")

  testthat::expect_equal(result$metric_choice, "iqr")
  test_data <- c(1, 2, 3, 4, 5)
  expected_iqr <- IQR(test_data)
  testthat::expect_equal(result$agg_fun(test_data), expected_iqr)
})

testthat::test_that("create_aggregation_function generates appropriate label", {
  result <- create_aggregation_function(metric = "median")

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

  result <- build_transcript_long(
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

  result <- build_transcript_long(
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
    build_transcript_long(
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
  result <- aggregate_transcript_data(df_long, agg_fun, pseudocount = 0)

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
  result <- aggregate_transcript_data(df_long, agg_fun, pseudocount = 1)

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

  top <- select_genes_from_results(res, top_n = 2)
  testthat::expect_equal(top, c("G2", "G3"))
})

testthat::test_that("select_genes_from_results removes duplicates", {
  res <- data.frame(
    genes = c("G1", "G1", "G2"),
    padj = c(0.001, 0.01, 0.05),
    stringsAsFactors = FALSE
  )

  top <- select_genes_from_results(res, top_n = 2)
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
    select_genes_from_results(res, top_n = 1),
    "must contain a 'genes' column"
  )
})
