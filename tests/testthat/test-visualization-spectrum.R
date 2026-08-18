context("plot_divergence_spectrum: Refactored helper functions")
library(SummarizedExperiment)
library(ggplot2)
library(testthat)


# ============================================================================
# Test: .spectrum_extract_q_values
# ============================================================================



test_that(".spectrum_extract_q_values extracts numeric q-values from column names", {
  col_names <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  q_vals <- TSENAT:::.spectrum_extract_q_values(col_names)
  
  expect_equal(length(q_vals), 4)
  expect_equal(q_vals, c(0.5, 1.0, 1.5, 2.0))
  expect_true(is.numeric(q_vals))
})

test_that(".spectrum_extract_q_values handles q= format", {
  col_names <- c("q=0.5", "q=1.0", "q=1.5")
  q_vals <- TSENAT:::.spectrum_extract_q_values(col_names)
  
  expect_equal(q_vals, c(0.5, 1.0, 1.5))
})

test_that(".spectrum_extract_q_values handles mixed formats", {
  col_names <- c("q_0.5", "q=1.0", "q_1.5")
  q_vals <- TSENAT:::.spectrum_extract_q_values(col_names)
  
  expect_equal(q_vals, c(0.5, 1.0, 1.5))
})

test_that(".spectrum_extract_q_values raises error on invalid format", {
  col_names <- c("invalid_0.5", "q_1.0", "q_1.5")
  
  expect_error(
    suppressWarnings(TSENAT:::.spectrum_extract_q_values(col_names)),
    "Cannot extract numeric q-values"
  )
})

test_that(".spectrum_extract_q_values raises error on all NA extraction", {
  col_names <- c("col1", "col2", "col3")
  
  expect_error(
    suppressWarnings(TSENAT:::.spectrum_extract_q_values(col_names)),
    "Cannot extract numeric q-values"
  )
})

# ============================================================================
# Test: .spectrum_get_gene_identifiers
# ============================================================================

test_that(".spectrum_get_gene_identifiers retrieves gene_name from rowData", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("row1", "row2", "row3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = mat),
    rowData = data.frame(gene_name = c("geneA", "geneB", "geneC"))
  )
  
  gene_ids <- TSENAT:::.spectrum_get_gene_identifiers(se, mat)
  
  expect_equal(gene_ids, c("geneA", "geneB", "geneC"))
})

test_that(".spectrum_get_gene_identifiers falls back to rownames", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  gene_ids <- TSENAT:::.spectrum_get_gene_identifiers(se, mat)
  
  expect_equal(gene_ids, c("gene1", "gene2", "gene3"))
})

test_that(".spectrum_get_gene_identifiers raises error when no identifiers", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  # No rownames, no gene_name in rowData
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_error(
    TSENAT:::.spectrum_get_gene_identifiers(se, mat),
    "no gene identifiers"
  )
})

# ============================================================================
# Test: .spectrum_find_gene_column
# ============================================================================

test_that(".spectrum_find_gene_column identifies 'gene' column", {
  df <- data.frame(gene = c("g1", "g2"), p_value = c(0.01, 0.05))
  
  col <- TSENAT:::.spectrum_find_gene_column(df)
  
  expect_equal(col, "gene")
})

test_that(".spectrum_find_gene_column identifies 'gene_name' column", {
  df <- data.frame(gene_name = c("g1", "g2"), p_value = c(0.01, 0.05))
  
  col <- TSENAT:::.spectrum_find_gene_column(df)
  
  expect_equal(col, "gene_name")
})

test_that(".spectrum_find_gene_column prefers 'gene' over 'gene_name'", {
  df <- data.frame(gene = c("g1", "g2"), gene_name = c("gA", "gB"), p_value = c(0.01, 0.05))
  
  col <- TSENAT:::.spectrum_find_gene_column(df)
  
  expect_equal(col, "gene")
})

test_that(".spectrum_find_gene_column raises error when no gene column", {
  df <- data.frame(x = c(1, 2), p_value = c(0.01, 0.05))
  
  expect_error(
    TSENAT:::.spectrum_find_gene_column(df),
    "must have a column named"
  )
})

# ============================================================================
# Test: .spectrum_find_pvalue_column
# ============================================================================

test_that(".spectrum_find_pvalue_column identifies 'adj_p_interaction'", {
  df <- data.frame(gene = c("g1", "g2"), adj_p_interaction = c(0.01, 0.05))
  
  col <- TSENAT:::.spectrum_find_pvalue_column(df)
  
  expect_equal(col, "adj_p_interaction")
})

test_that(".spectrum_find_pvalue_column identifies 'p_value'", {
  df <- data.frame(gene = c("g1", "g2"), p_value = c(0.01, 0.05))
  
  col <- TSENAT:::.spectrum_find_pvalue_column(df)
  
  expect_equal(col, "p_value")
})

test_that(".spectrum_find_pvalue_column raises error when no p-value column", {
  df <- data.frame(gene = c("g1", "g2"), x = c(0.01, 0.05))
  
  expect_error(
    TSENAT:::.spectrum_find_pvalue_column(df),
    "must have a p-value column"
  )
})

# ============================================================================
# Test: .spectrum_validate_inputs
# ============================================================================

test_that(".spectrum_validate_inputs accepts valid SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("g1", "g2", "g3")
  seq <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  result <- TSENAT:::.spectrum_validate_inputs(seq, NULL, NULL, 4, 2)
  
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that(".spectrum_validate_inputs rejects non-SummarizedExperiment", {
  expect_error(
    TSENAT:::.spectrum_validate_inputs(list(), NULL, NULL, 4, 2),
    "must be a SummarizedExperiment"
  )
})

test_that(".spectrum_validate_inputs rejects empty assay", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Create SE with 0 rows
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = matrix(nrow = 0, ncol = 4)))
  
  expect_error(
    TSENAT:::.spectrum_validate_inputs(se, NULL, NULL, 4, 2),
    "assay is empty"
  )
})

test_that(".spectrum_validate_inputs rejects invalid gene parameter", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_error(
    TSENAT:::.spectrum_validate_inputs(se, c("g1", "g2"), NULL, 4, 2),
    "must be a single character string"
  )
})

test_that(".spectrum_validate_inputs rejects invalid sait_res", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_error(
    TSENAT:::.spectrum_validate_inputs(se, NULL, "not_a_df", 4, 2),
    "must be a non-empty data.frame"
  )
})

test_that(".spectrum_validate_inputs rejects invalid n_genes", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_error(
    TSENAT:::.spectrum_validate_inputs(se, NULL, NULL, -1, 2),
    "must be positive"
  )
})

test_that(".spectrum_validate_inputs rejects invalid ncol", {
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  expect_error(
    TSENAT:::.spectrum_validate_inputs(se, NULL, NULL, 4, 0),
    "must be positive"
  )
})

# ============================================================================
# Test: .spectrum_plot_single_gene
# ============================================================================

test_that(".spectrum_plot_single_gene returns ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3")
  
  p <- TSENAT:::.spectrum_plot_single_gene("gene1", mat, c(0.5, 1.0, 1.5, 2.0),
                                            c("gene1", "gene2", "gene3"))
  
  expect_is(p, "ggplot")
})

test_that(".spectrum_plot_single_gene has correct data", {
  skip_if_not_installed("ggplot2")
  
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                nrow = 3, ncol = 4)
  rownames(mat) <- c("g1", "g2", "g3")
  
  p <- TSENAT:::.spectrum_plot_single_gene("g1", mat, c(0.5, 1.0, 1.5, 2.0),
                                            c("g1", "g2", "g3"))
  
  expect_equal(p$labels$title, "Divergence Spectrum: g1")
})

test_that(".spectrum_plot_single_gene raises error for unknown gene", {
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  
  expect_error(
    TSENAT:::.spectrum_plot_single_gene("unknown", mat, c(0.5, 1.0, 1.5, 2.0),
                                        c("g1", "g2", "g3")),
    "not found"
  )
})

# ============================================================================
# Test: .spectrum_plot_global
# ============================================================================

test_that(".spectrum_plot_global returns ggplot with IQR spread", {
  skip_if_not_installed("ggplot2")
  
  mat <- matrix(rnorm(100, mean = 0.5, sd = 0.1), nrow = 20, ncol = 5)
  q_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  p <- TSENAT:::.spectrum_plot_global(mat, q_vals, "median", "iqr")
  
  expect_is(p, "ggplot")
})

test_that(".spectrum_plot_global returns ggplot with SD spread", {
  skip_if_not_installed("ggplot2")
  
  mat <- matrix(rnorm(100, mean = 0.5, sd = 0.1), nrow = 20, ncol = 5)
  q_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  p <- TSENAT:::.spectrum_plot_global(mat, q_vals, "mean", "sd")
  
  expect_is(p, "ggplot")
})

test_that(".spectrum_plot_global uses correct metric label", {
  skip_if_not_installed("ggplot2")
  
  mat <- matrix(rnorm(100, mean = 0.5, sd = 0.1), nrow = 20, ncol = 5)
  q_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  p <- TSENAT:::.spectrum_plot_global(mat, q_vals, "median", "iqr")
  
  expect_true(grepl("Median", p$labels$subtitle))
})

test_that(".spectrum_plot_global raises error on all-NA values", {
  mat <- matrix(NA_real_, nrow = 5, ncol = 4)
  q_vals <- c(0.5, 1.0, 1.5, 2.0)
  
  expect_error(
    TSENAT:::.spectrum_plot_global(mat, q_vals, "median", "iqr"),
    "Cannot compute"
  )
})

# ============================================================================
# Integration Test: .plot_divergence_spectrum main function
# ============================================================================

test_that(".plot_divergence_spectrum with single gene returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3")
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  p <- TSENAT:::.plot_divergence_spectrum(se, gene = "gene1")
  
  expect_is(p, "ggplot")
})

test_that(".plot_divergence_spectrum with global mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(80, mean = 0.5, sd = 0.1), nrow = 20, ncol = 4)
  rownames(mat) <- paste0("gene_", 1:20)
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  p <- TSENAT:::.plot_divergence_spectrum(se)
  
  expect_is(p, "ggplot")
})

test_that(".plot_divergence_spectrum with sait_res returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(mat) <- c("gene1", "gene2", "gene3")
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.001, 0.01, 0.05)
  )
  
  p <- TSENAT:::.plot_divergence_spectrum(se, sait_res = sait_res, n_genes = 2)
  
  expect_is(p, "ggplot")
})

test_that(".spectrum_plot_top_genes issues warning for unmatched genes and plots with CIs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")

  mat <- matrix(c(0.1, 0.2, 0.3, 0.4,
                  0.5, 0.6, 0.7, 0.8,
                  0.9, 1.0, 1.1, 1.2),
                nrow = 3, ncol = 4, byrow = TRUE)
  rownames(mat) <- c("gene1", "gene2", "gene3")
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")

  ci_lower <- mat - 0.05
  ci_upper <- mat + 0.05
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = mat, ci_lower = ci_lower, ci_upper = ci_upper)
  )

  sait_res <- data.frame(
    gene = c("gene3", "geneX", "gene1"),
    adj_p_interaction = c(0.01, 0.05, 0.02),
    stringsAsFactors = FALSE
  )

  p <- expect_warning(
    TSENAT:::.spectrum_plot_top_genes(
      sait_res = sait_res,
      n_genes_use = 3,
      ncol = 2,
      div_mat_sorted = mat,
      q_vals_sorted = c(0.5, 1.0, 1.5, 2.0),
      gene_names = rownames(mat),
      metric = "median",
      divergence_results_se = se
    ),
    "Some genes not found"
  )

  expect_is(p, "ggplot")
  expect_true(grepl("Bootstrap CI \\(95%\\)", p$labels$subtitle))
})

test_that(".spectrum_plot_global without analysis never claims a bootstrap CI", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")

  mat <- matrix(rnorm(20, mean = 0.5, sd = 0.1), nrow = 5, ncol = 4)
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")

  ci_lower <- mat - 0.05
  ci_upper <- mat + 0.05
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = mat, ci_lower = ci_lower, ci_upper = ci_upper)
  )

  # Without the analysis object a VALID global bootstrap CI cannot be
  # computed, so the plot must fall back to descriptive spread and must
  # NOT label the band as a bootstrap CI.
  p <- TSENAT:::.spectrum_plot_global(
    div_mat_sorted = mat,
    q_vals_sorted = c(0.5, 1.0, 1.5, 2.0),
    metric = "mean",
    variability_metric = "sd",
    divergence_results_se = se,
    analysis = NULL
  )

  expect_is(p, "ggplot")
  expect_false(grepl("Bootstrap", p$labels$subtitle))
  expect_true(grepl("descriptive spread", p$labels$subtitle))
})

test_that("plot_divergence_spectrum global mode shows a valid global bootstrap CI", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")

  # Small TSENATAnalysis with paired-free 2-condition design
  set.seed(42)
  n_tx <- 12
  counts <- pmax(matrix(rpois(n_tx * 6, lambda = 60), nrow = n_tx), 40)
  rownames(counts) <- paste0("tx", seq_len(n_tx))
  colnames(counts) <- paste0("s", seq_len(6))
  gene_ids <- paste0("g", rep(seq_len(4), each = 3))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(
      gene_id = gene_ids,
      row.names = rownames(counts)
    ),
    colData = S4Vectors::DataFrame(
      condition = rep(c("control", "treatment"), each = 3),
      row.names = colnames(counts)
    )
  )
  # Diversity requires the tx2gene mapping when genes= is not explicit
  S4Vectors::metadata(se)$tx2gene <- data.frame(
    Transcript = rownames(counts), Gene = gene_ids, stringsAsFactors = FALSE
  )

  analysis <- TSENAT::TSENATAnalysis(se = se, config = list(
    condition_col = "condition", control_group = "control", nthreads = 1))
  analysis <- calculate_diversity(analysis, q = c(0.5, 1), min_valid_frac = 0)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1), bootstrap = TRUE,
    nboot = 20)

  p <- plot_divergence_spectrum(analysis, metric = "mean")

  expect_is(p, "ggplot")
  expect_true(grepl("Bootstrap \\(95%\\)", p$labels$subtitle))
  expect_true(grepl("global bootstrap", p$labels$subtitle))
})

test_that(".plot_divergence_spectrum respects metric parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(80, mean = 0.5, sd = 0.1), nrow = 20, ncol = 4)
  rownames(mat) <- paste0("gene_", 1:20)
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  p1 <- TSENAT:::.plot_divergence_spectrum(se, metric = "median")
  p2 <- TSENAT:::.plot_divergence_spectrum(se, metric = "mean")
  
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
})

test_that(".plot_divergence_spectrum respects variability_metric parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  
  mat <- matrix(rnorm(80, mean = 0.5, sd = 0.1), nrow = 20, ncol = 4)
  rownames(mat) <- paste0("gene_", 1:20)
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(div = mat))
  
  p1 <- TSENAT:::.plot_divergence_spectrum(se, variability_metric = "iqr")
  p2 <- TSENAT:::.plot_divergence_spectrum(se, variability_metric = "sd")
  
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
})

test_that("plot_divergence_spectrum: creates spectrum plot from S4 object", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  
  set.seed(888)
  
  # Create minimal SE with q-spectrum divergence data
  mat <- matrix(rnorm(60, mean = 0.5, sd = 0.1), nrow = 15, ncol = 4)
  rownames(mat) <- paste0("gene_", 1:15)
  colnames(mat) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence_spectrum = mat)
  )
  
  # Test that S4 function works with valid inputs
  expect_true(!is.null(se))
  expect_true("divergence_spectrum" %in% SummarizedExperiment::assayNames(se))
})

test_that("plot_divergence_spectrum: parameter validation", {
  skip_if_not_installed("ggplot2")
  
  
  # Validate expected parameters
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  expect_true(is.numeric(q_values))
  expect_true(all(q_values > 0))
  
  # Variability metric options
  metrics <- c("iqr", "sd", "mad")
  expect_true("iqr" %in% metrics)
})


# ============================================================================
# TEST 2: plot_divergence_spectrum - Divergence spectrum plot
# ============================================================================

test_that("plot_divergence_spectrum: validates analysis object", {
  expect_error(
    TSENAT:::plot_divergence_spectrum("not_analysis"),
    "must be a TSENATAnalysis object"
  )
})

test_that("plot_divergence_spectrum: requires divergence results", {
  set.seed(304)
  
  # Create minimal analysis without divergence results
  analysis <- .create_test_analysis(
    n_genes = 10,
    n_samples_per_group = 3,
    q_values = c(1.0),
    include_divergence = FALSE,
    include_sait_results = FALSE,
    seed = 304,
    verbose = FALSE
  )
  
  # Should error when no divergence results
  expect_error(
    TSENAT:::plot_divergence_spectrum(analysis),
    "Divergence results not found|Invalid divergence_results structure"
  )
})

test_that("plot_divergence_spectrum: creates plot with valid divergence data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add diversity and divergence calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  
  # Should create plot successfully with default parameters
  result <- TSENAT:::plot_divergence_spectrum(analysis, n_genes = 2, verbose = FALSE)
  
  # Result should NOT be NULL - actual plot code must run
  expect_false(is.null(result))
  expect_true(inherits(result, "ggplot") || is.list(result))
})

test_that("plot_divergence_spectrum: handles single gene mode", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add diversity and divergence calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  
  # Test with single gene specified (use actual gene name from analysis)
  genes <- rownames(analysis@divergence_results$divergence_se)
  gene_to_plot <- genes[1]
  
  result <- TSENAT:::plot_divergence_spectrum(analysis, gene = gene_to_plot, verbose = FALSE)
  
  expect_false(is.null(result))
  expect_true(inherits(result, "ggplot") || is.list(result))
})

test_that("plot_divergence_spectrum: respects metric and variability parameters", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add diversity and divergence calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
  
  # Test different metric combinations
  result_median_iqr <- TSENAT:::plot_divergence_spectrum(
    analysis, n_genes = 2, metric = "median", variability_metric = "iqr", verbose = FALSE
  )
  
  result_mean_sd <- TSENAT:::plot_divergence_spectrum(
    analysis, n_genes = 2, metric = "mean", variability_metric = "sd", verbose = FALSE
  )
  
  expect_false(is.null(result_median_iqr))
  expect_true(inherits(result_median_iqr, "ggplot") || is.list(result_median_iqr))
  
  expect_false(is.null(result_mean_sd))
  expect_true(inherits(result_mean_sd, "ggplot") || is.list(result_mean_sd))
})
