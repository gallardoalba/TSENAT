library(testthat)

context("High-Level Analysis Workflow: Method Concordance")

# =============================================================================
# TEST: compute_method_concordance - Basic functionality
# =============================================================================

test_that("compute_method_concordance computes correlation correctly", {
  # Create sample data
  gam_results <- data.frame(
    gene = paste0("GENE_", 1:20),
    p_interaction = runif(20),
    adj_p_interaction = runif(20),
    effect_size = rnorm(20)
  )
  
  kw_results <- data.frame(
    gene = paste0("GENE_", 1:20),
    p_value = runif(20),
    adj_p_value = runif(20),
    effect_size_eta2 = runif(20)
  )
  
  result <- compute_method_concordance(gam_results, kw_results)
  
  expect_true(is.list(result))
  expect_true("comparison_df" %in% names(result))
  expect_true("spearman_rho" %in% names(result))
  expect_true("high_conf" %in% names(result))
  expect_true("agreement_table" %in% names(result))
  
  # Check correlation is valid
  expect_true(is.numeric(result$spearman_rho))
  expect_true(result$spearman_rho >= -1 && result$spearman_rho <= 1)
})

test_that("compute_method_concordance identifies agreement categories", {
  gam_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_interaction = c(0.001, 0.01, 0.5, 0.5, 0.001, 0.01, 0.5, 0.5, 0.001, 0.01),
    adj_p_interaction = c(0.01, 0.05, 0.5, 0.5, 0.01, 0.05, 0.5, 0.5, 0.01, 0.05)
  )
  
  kw_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = c(0.001, 0.5, 0.01, 0.5, 0.5, 0.5, 0.001, 0.5, 0.5, 0.05),
    adj_p_value = c(0.01, 0.5, 0.05, 0.5, 0.5, 0.5, 0.01, 0.5, 0.5, 0.05)
  )
  
  result <- compute_method_concordance(gam_results, kw_results)
  
  expect_true(nrow(result$comparison_df) > 0)
  expect_true("agreement" %in% colnames(result$comparison_df))
  
  # Check agreement categories exist
  categories <- unique(result$comparison_df$agreement)
  expect_true(any(c("Both significant", "GAM only", "Friedman only", "Neither significant") %in% categories))
})

test_that("compute_method_concordance extracts high-confidence genes", {
  gam_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_interaction = c(0.001, 0.5, rep(0.5, 8)),
    adj_p_interaction = c(0.01, 0.5, rep(0.5, 8))
  )
  
  kw_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = c(0.001, 0.5, rep(0.5, 8)),
    adj_p_value = c(0.01, 0.5, rep(0.5, 8))
  )
  
  result <- compute_method_concordance(gam_results, kw_results)
  
  # First gene should be in high_conf
  expect_true(nrow(result$high_conf) >= 1)
  expect_true("GENE_1" %in% result$high_conf$gene)
})

test_that("compute_method_concordance handles missing p_interaction column", {
  gam_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = runif(10)
  )
  
  kw_results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = runif(10)
  )
  
  expect_error(
    compute_method_concordance(gam_results, kw_results),
    "p_interaction"
  )
})

test_that("compute_method_concordance handles non-data.frame input", {
  gam_results <- list(a = 1, b = 2)
  kw_results <- data.frame(gene = 1:10, p_value = runif(10))
  
  expect_error(
    compute_method_concordance(gam_results, kw_results),
    "data.frame"
  )
})

test_that("compute_method_concordance handles mismatched genes", {
  gam_results <- data.frame(
    gene = paste0("GENE_A_", 1:10),
    p_interaction = runif(10),
    adj_p_interaction = runif(10)
  )
  
  kw_results <- data.frame(
    gene = paste0("GENE_B_", 1:10),
    p_value = runif(10),
    adj_p_value = runif(10)
  )
  
  result <- compute_method_concordance(gam_results, kw_results)
  
  # No common genes - expect NULL or empty results
  expect_true(is.null(result$comparison_df) || nrow(result$comparison_df) == 0)
})

# =============================================================================
# TEST: plot_method_concordance - Visualization
# =============================================================================

context("High-Level Analysis Workflow: Concordance Plotting")

test_that("plot_method_concordance creates valid plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("cowplot")
  
  comparison_df <- data.frame(
    gene = paste0("GENE_", 1:20),
    p_gam = runif(20),
    p_friedman = runif(20),
    agreement = sample(c("Both significant", "GAM only", "Friedman only", "Neither significant"),
                      size = 20, replace = TRUE)
  )
  
  plot <- plot_method_concordance(comparison_df)
  
  # Check that result is a grob object
  expect_true(methods::is(plot, "grob") || methods::is(plot, "gtable"))
})

test_that("plot_method_concordance rejects empty data.frame", {
  skip_if_not_installed("ggplot2")
  
  comparison_df <- data.frame()
  
  expect_error(
    plot_method_concordance(comparison_df),
    "non-empty"
  )
})

test_that("plot_method_concordance rejects NULL input", {
  skip_if_not_installed("ggplot2")
  
  expect_error(
    plot_method_concordance(NULL),
    "non-empty"
  )
})

test_that("plot_method_concordance requires required columns", {
  skip_if_not_installed("ggplot2")
  
  comparison_df <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_gam = runif(10)
    # Missing p_friedman and agreement columns
  )
  
  expect_error(
    plot_method_concordance(comparison_df),
    "required columns"
  )
})

# =============================================================================
# CONTEXT: Filter Analysis Workflow
# =============================================================================

context("Analysis Workflow: Filter Analysis")

# Helper functions
make_test_se <- function(n_genes = 100, n_samples = 20) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50), nrow = n_genes)
  rownames(counts) <- paste0("TX_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Add TPM data to avoid warnings during filtering
  # Simple TPM calculation: scale counts to sum to 1 million per sample
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
}

test_that("filter_analysis modifies SE in analysis object", {
  se <- make_test_se(n_genes = 100, n_samples = 20)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:10, each = 2),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  analysis <- TSENATAnalysis(se)
  
  # Apply filtering
  filtered_analysis <- filter_analysis_s4(analysis, stringency = "severe", verbose = FALSE)
  
  # Check that analysis is returned
  expect_true(inherits(filtered_analysis, "TSENATAnalysis"))
  
  # Check that SE was modified
  expect_true(nrow(filtered_analysis@se) <= nrow(analysis@se))
})

test_that("filter_analysis preserves colData", {
  se <- make_test_se(n_genes = 100, n_samples = 20)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:10, each = 2),
    condition = rep(c("control", "treatment"), 10),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  
  analysis <- TSENATAnalysis(se)
  filtered_analysis <- filter_analysis_s4(analysis, stringency = "soft", verbose = FALSE)
  
  # colData should preserve original columns (sample_id is added by constructor)
  filtered_coldata <- SummarizedExperiment::colData(filtered_analysis@se)
  expect_true("pair_id" %in% colnames(filtered_coldata))
  expect_true("condition" %in% colnames(filtered_coldata))
  expect_true("sample_id" %in% colnames(filtered_coldata))
  # Check expected column count: pair_id + condition + sample_id (added by constructor)
  expect_equal(ncol(filtered_coldata), 3)
})

test_that("filter_analysis validates input type", {
  bad_input <- "not_an_analysis"
  
  expect_error(
    filter_analysis_s4(bad_input),
    "TSENATAnalysis"
  )
})

test_that("filter_analysis accepts stringency parameter", {
  se <- make_test_se(n_genes = 200, n_samples = 30)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:15, each = 2),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  analysis <- TSENATAnalysis(se)
  
  # Test different stringency levels
  for (stringency in c("soft", "medium", "severe")) {
    result <- filter_analysis_s4(analysis, stringency = stringency, verbose = FALSE)
    expect_true(inherits(result, "TSENATAnalysis"))
  }
})

# =============================================================================
# CONTEXT: Build Analysis Workflow
# =============================================================================

context("Analysis Workflow: Build Analysis")

test_that("build_analysis creates valid TSENATAnalysis object", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  analysis <- build_analysis_s4(readcounts = counts, tx2gene = tx2gene)
  
  expect_true(inherits(analysis, "TSENATAnalysis"))
  expect_true(inherits(analysis@se, "SummarizedExperiment"))
})

test_that("build_analysis includes metadata in config when provided", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  metadata <- data.frame(
    condition = rep(c("control", "treatment"), 5),
    row.names = colnames(counts)
  )
  
  analysis <- build_analysis_s4(readcounts = counts, tx2gene = tx2gene, metadata = metadata)
  
  expect_true("metadata" %in% names(analysis@config))
})

test_that("build_analysis accepts custom config parameters", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  config <- list(analysis_id = "TEST001", version = "1.0")
  analysis <- build_analysis_s4(readcounts = counts, tx2gene = tx2gene, config = config)
  
  expect_equal(analysis@config$analysis_id, "TEST001")
  expect_equal(analysis@config$version, "1.0")
})

test_that("build_analysis initializes empty result slots", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  analysis <- build_analysis_s4(readcounts = counts, tx2gene = tx2gene)
  
  expect_true(is.list(analysis@diversity_results) && length(analysis@diversity_results) == 0)
  expect_true(is.list(analysis@divergence_results) && length(analysis@divergence_results) == 0)
  expect_true(is.list(analysis@lm_results) && length(analysis@lm_results) == 0)
  expect_true(is.list(analysis@jackknife_results) && length(analysis@jackknife_results) == 0)
})

test_that("build_analysis accepts TPM and effective_length", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  eff_len <- rep(1000, nrow(counts))
  
  analysis <- build_analysis_s4(
    readcounts = counts,
    tx2gene = tx2gene,
    tpm = tpm,
    effective_length = eff_len
  )
  
  expect_true(inherits(analysis, "TSENATAnalysis"))
})

test_that("build_analysis with custom assay name", {
  counts <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts) <- paste0("TX_", 1:10)
  colnames(counts) <- paste0("Sample_", 1:10)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gene = paste0("GENE_", rep(1:5, 2))
  )
  
  analysis <- build_analysis_s4(
    readcounts = counts,
    tx2gene = tx2gene,
    assay_name = "normalized_counts"
  )
  
  expect_true("normalized_counts" %in% names(SummarizedExperiment::assays(analysis@se)))
})
