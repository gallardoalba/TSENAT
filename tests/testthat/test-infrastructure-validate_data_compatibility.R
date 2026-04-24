library(testthat)

context("Data Validation: Gene Name Alignment")

# Helper functions
make_test_se <- function(n_genes = 10, n_samples = 5) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10), nrow = n_genes)
  rownames(counts) <- paste0("GENE_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
}

make_test_results <- function(n_genes = 10, gene_col = "gene") {
  df <- data.frame(
    p_value = runif(n_genes),
    effect_size = rnorm(n_genes),
    row.names = paste0("GENE_", 1:n_genes)
  )
  if (gene_col == "gene") {
    df$gene <- rownames(df)
  } else if (gene_col == "gene_id") {
    df$gene_id <- rownames(df)
  }
  df
}

# =============================================================================
# TEST: validate_gene_names - Aligned case
# =============================================================================

test_that("validate_gene_names returns TRUE for perfectly aligned genes", {
  se <- make_test_se(n_genes = 10)
  results <- make_test_results(n_genes = 10, gene_col = "gene")
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    se_name = "SE",
    results_name = "Results",
    verbose = FALSE
  )
  
  expect_true(validation$is_aligned)
  expect_equal(validation$n_se_genes, 10)
  expect_equal(validation$n_results_genes, 10)
  expect_equal(validation$n_matched, 10)
  expect_equal(validation$n_se_only, 0)
  expect_equal(validation$n_results_only, 0)
})

# =============================================================================
# TEST: validate_gene_names - Misaligned genes
# =============================================================================

test_that("validate_gene_names detects subset mismatch", {
  se <- make_test_se(n_genes = 10)
  results <- make_test_results(n_genes = 7, gene_col = "gene")
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    se_name = "SE",
    results_name = "Results",
    verbose = FALSE
  )
  
  expect_false(validation$is_aligned)
  expect_equal(validation$n_se_genes, 10)
  expect_equal(validation$n_results_genes, 7)
  expect_equal(validation$n_matched, 7)
  expect_equal(validation$n_se_only, 3)
  expect_equal(validation$n_results_only, 0)
})

test_that("validate_gene_names detects superset mismatch", {
  se <- make_test_se(n_genes = 5)
  results <- make_test_results(n_genes = 10, gene_col = "gene")
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    se_name = "SE",
    results_name = "Results",
    verbose = FALSE
  )
  
  expect_false(validation$is_aligned)
  expect_equal(validation$n_se_genes, 5)
  expect_equal(validation$n_results_genes, 10)
  expect_equal(validation$n_se_only, 0)
  expect_equal(validation$n_results_only, 5)
})

test_that("validate_gene_names handles completely disjoint genes", {
  se <- make_test_se(n_genes = 5)
  rownames(se) <- paste0("SETGENE_", 1:5)
  
  results <- make_test_results(n_genes = 5, gene_col = "gene")
  results$gene <- paste0("RESGENE_", 1:5)
  rownames(results) <- paste0("RESGENE_", 1:5)
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    se_name = "SE",
    results_name = "Results",
    verbose = FALSE
  )
  
  expect_false(validation$is_aligned)
  expect_equal(validation$n_matched, 0)
  expect_equal(validation$n_se_only, 5)
  expect_equal(validation$n_results_only, 5)
})

# =============================================================================
# TEST: validate_gene_names - NULL/empty inputs
# =============================================================================

test_that("validate_gene_names handles NULL SE", {
  results <- make_test_results(n_genes = 10, gene_col = "gene")
  
  validation <- .validate_gene_names(
    se = NULL,
    results = results,
    verbose = FALSE
  )
  
  expect_false(validation$is_aligned)
  expect_true("error" %in% names(validation))
})

test_that("validate_gene_names handles NULL results", {
  se <- make_test_se(n_genes = 10)
  
  validation <- .validate_gene_names(
    se = se,
    results = NULL,
    verbose = FALSE
  )
  
  expect_false(validation$is_aligned)
  expect_true("error" %in% names(validation))
})

# =============================================================================
# TEST: validate_gene_names - Flexible gene column detection
# =============================================================================

test_that("validate_gene_names detects 'gene_id' column", {
  se <- make_test_se(n_genes = 10)
  results <- make_test_results(n_genes = 10, gene_col = "gene_id")
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    verbose = FALSE
  )
  
  expect_true(validation$is_aligned)
})

test_that("validate_gene_names detects 'gene_name' column", {
  se <- make_test_se(n_genes = 10)
  results <- make_test_results(n_genes = 10, gene_col = "gene_id")
  colnames(results)[which(colnames(results) == "gene_id")] <- "gene_name"
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    verbose = FALSE
  )
  
  expect_true(validation$is_aligned)
})

test_that("validate_gene_names uses rownames when gene column absent", {
  se <- make_test_se(n_genes = 10)
  results <- data.frame(
    p_value = runif(10),
    effect_size = rnorm(10),
    row.names = paste0("GENE_", 1:10)
  )
  
  validation <- .validate_gene_names(
    se = se,
    results = results,
    verbose = FALSE
  )
  
  expect_true(validation$is_aligned)
})

# =============================================================================
# CONTEXT: Data Validation - SE Dimensions
# =============================================================================

context("Data Validation: SummarizedExperiment Dimensions")

test_that("validate_se_dimensions passes for valid SE", {
  se <- make_test_se(n_genes = 20, n_samples = 10)
  
  validation <- .validate_se_dimensions(se, verbose = FALSE)
  
  expect_true(validation$is_valid)
  expect_equal(validation$n_genes, 20)
  expect_equal(validation$n_samples, 10)
  expect_length(validation$issues, 0)
})

test_that("validate_se_dimensions checks expected sample count", {
  se <- make_test_se(n_genes = 20, n_samples = 10)
  
  validation <- .validate_se_dimensions(
    se = se,
    expected_n_samples = 5,
    verbose = FALSE
  )
  
  expect_false(validation$is_valid)
  expect_length(validation$issues, 1)
  expect_match(validation$issues[1], "Sample count mismatch")
})

test_that("validate_se_dimensions checks expected assays", {
  se <- make_test_se(n_genes = 20, n_samples = 10)
  
  validation <- .validate_se_dimensions(
    se = se,
    expected_assays = c("counts", "missing_assay"),
    verbose = FALSE
  )
  
  expect_false(validation$is_valid)
  expect_length(validation$issues, 1)
  expect_match(validation$issues[1], "Missing assays")
})

test_that("validate_se_dimensions passes with properly aligned colData", {
  se <- make_test_se(n_genes = 20, n_samples = 10)
  # Set valid colData with correct sample count
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
    sample_id = 1:10,
    condition = rep(c("A", "B"), 5),
    row.names = colnames(se)
  )
  
  validation <- .validate_se_dimensions(se, verbose = FALSE)
  
  expect_true(validation$is_valid)
  expect_length(validation$issues, 0)
})

test_that("validate_se_dimensions passes with properly aligned rowData", {
  se <- make_test_se(n_genes = 20, n_samples = 10)
  # Set valid rowData with correct gene count
  SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(
    gene_id = paste0("GENEID_", 1:20),
    biotype = rep(c("protein_coding", "lncRNA"), 10),
    row.names = rownames(se)
  )
  
  validation <- .validate_se_dimensions(se, verbose = FALSE)
  
  expect_true(validation$is_valid)
  expect_length(validation$issues, 0)
})

# =============================================================================
# CONTEXT: Data Validation - SAIT Results
# =============================================================================

context("Data Validation: SAIT Results Structure")

test_that("validate_sait_results passes for valid results", {
  results <- make_test_results(n_genes = 10, gene_col = "gene")
  
  validation <- .validate_sait_results(results, verbose = FALSE)
  
  expect_true(validation$is_valid)
  expect_equal(validation$n_results, 10)
  expect_length(validation$issues, 0)
})

test_that("validate_sait_results rejects non-data.frame", {
  results <- list(a = 1, b = 2)
  
  validation <- .validate_sait_results(results, verbose = FALSE)
  
  expect_false(validation$is_valid)
  expect_length(validation$issues, 1)
})

test_that("validate_sait_results accepts list with 'results' element", {
  results_df <- make_test_results(n_genes = 10, gene_col = "gene")
  results_list <- list(results = results_df)
  
  validation <- .validate_sait_results(results_list, verbose = FALSE)
  
  expect_true(validation$is_valid)
})

test_that("validate_sait_results detects missing p-value column", {
  results <- data.frame(
    gene = paste0("GENE_", 1:10),
    effect_size = rnorm(10)
  )
  
  validation <- .validate_sait_results(results, verbose = FALSE)
  
  expect_false(validation$is_valid)
  expect_match(validation$issues[1], "Missing p-value column")
})

test_that("validate_sait_results accepts various p-value column names", {
  for (pval_col in c("p_value", "p_interaction", "p_raw", "adj_p_interaction")) {
    results <- data.frame(
      gene = paste0("GENE_", 1:10)
    )
    results[[pval_col]] <- runif(10)
    
    validation <- .validate_sait_results(results, verbose = FALSE)
    
    expect_true(validation$is_valid, info = paste("Failed for column:", pval_col))
  }
})

test_that("validate_sait_results detects invalid p-values", {
  results <- data.frame(
    gene = paste0("GENE_", 1:10),
    p_value = c(runif(5), -0.1, 1.5, 2.0, NA, 0.05)
  )
  
  validation <- .validate_sait_results(results, verbose = FALSE)
  
  expect_false(validation$is_valid)
  expect_match(validation$issues[1], "Invalid p-values")
})

test_that("validate_sait_results checks gene alignment when expected_genes provided", {
  results <- data.frame(
    gene = paste0("GENE_", 1:5),
    p_value = runif(5)
  )
  expected_genes <- paste0("GENE_", 1:10)
  
  validation <- .validate_sait_results(
    results,
    expected_genes = expected_genes,
    verbose = FALSE
  )
  
  expect_false(validation$is_valid)
  expect_match(validation$issues[1], "Gene mismatch")
})

# =============================================================================
# CONTEXT: Data Validation - Plot Data
# =============================================================================

context("Data Validation: Comprehensive Plot Data Validation")

test_that("validate_plot_data passes for valid SE and results", {
  se <- make_test_se(n_genes = 10, n_samples = 5)
  results <- make_test_results(n_genes = 10, gene_col = "gene")
  
  validation <- .validate_plot_data(se, results, verbose = FALSE, stop_on_error = FALSE)
  
  expect_length(validation, 0)
})

test_that("validate_plot_data works without results", {
  se <- make_test_se(n_genes = 10, n_samples = 5)
  
  validation <- .validate_plot_data(se, sait_results = NULL, verbose = FALSE, stop_on_error = FALSE)
  
  expect_length(validation, 0)
})

test_that("validate_plot_data stops on error when requested", {
  se <- make_test_se(n_genes = 10, n_samples = 5)
  results <- make_test_results(n_genes = 5, gene_col = "gene")  # Mismatched
  
  expect_error(
    .validate_plot_data(se, results, verbose = FALSE, stop_on_error = TRUE)
  )
})

test_that("validate_plot_data returns issues without stopping if requested", {
  se <- make_test_se(n_genes = 10, n_samples = 5)
  results <- make_test_results(n_genes = 5, gene_col = "gene")  # Mismatched
  
  expect_warning(
    validation <- .validate_plot_data(se, results, verbose = FALSE, stop_on_error = FALSE)
  )
  
  expect_length(validation, 1)
  expect_true("gene_alignment" %in% names(validation))
})

test_that("validate_plot_data consolidates multiple issues", {
  # Create SE with proper colData and results with gene mismatch
  se <- make_test_se(n_genes = 10, n_samples = 5)
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
    sample_id = 1:5,
    row.names = colnames(se)
  )
  
  results <- make_test_results(n_genes = 5, gene_col = "gene")  # Gene mismatch
  
  expect_warning(
    validation <- .validate_plot_data(se, results, verbose = FALSE, stop_on_error = FALSE)
  )
  
  # Should report both SE dimension issue and gene alignment issue
  expect_true(length(validation) > 0)
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE: .validate_se_dimensions() - NEWLY ADDED FOR COVERAGE
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".validate_se_dimensions accepts valid SummarizedExperiment", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    result <- TSENAT:::.validate_se_dimensions(
        se = se,
        verbose = FALSE
    )
    
    expect_true(result$is_valid)
})

test_that(".validate_se_dimensions checks sample count", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    expect_message(
        result <- TSENAT:::.validate_se_dimensions(
            se = se,
            expected_n_samples = 3,
            verbose = TRUE
        ),
        "dimension checks passed"
    )
    
    expect_true(result$is_valid)
})

test_that(".validate_se_dimensions detects sample count mismatch", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    expect_message(
        result <- TSENAT:::.validate_se_dimensions(
            se = se,
            expected_n_samples = 5,
            verbose = TRUE
        ),
        "mismatch"
    )
    
    expect_false(result$is_valid)
})

test_that(".validate_se_dimensions checks for expected assays", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    assay_names <- names(SummarizedExperiment::assays(se))
    
    result <- TSENAT:::.validate_se_dimensions(
        se = se,
        expected_assays = assay_names,
        verbose = FALSE
    )
    
    expect_true(result$is_valid)
})

test_that(".validate_se_dimensions detects missing assays", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    expect_message(
        result <- TSENAT:::.validate_se_dimensions(
            se = se,
            expected_assays = c("missing_assay"),
            verbose = TRUE
        ),
        "Missing assays"
    )
    
    expect_false(result$is_valid)
})

test_that(".validate_se_dimensions checks colData consistency", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    # colData should be consistent with ncol(se)
    result <- TSENAT:::.validate_se_dimensions(
        se = se,
        verbose = FALSE
    )
    
    expect_true(result$is_valid)
})

test_that(".validate_se_dimensions checks rowData consistency", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    # rowData should be consistent with nrow(se)
    result <- TSENAT:::.validate_se_dimensions(
        se = se,
        verbose = FALSE
    )
    
    expect_true(result$is_valid)
})

test_that(".validate_se_dimensions rejects non-SummarizedExperiment", {
    not_se <- data.frame(x = 1:3)
    
    expect_error(
        TSENAT:::.validate_se_dimensions(se = not_se),
        "must be a SummarizedExperiment"
    )
})

test_that(".validate_se_dimensions provides verbose output on success", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    expect_message(
        TSENAT:::.validate_se_dimensions(
            se = se,
            verbose = TRUE
        ),
        "All dimension checks passed"
    )
})

test_that(".validate_se_dimensions provides verbose output on failure", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    
    expect_message(
        TSENAT:::.validate_se_dimensions(
            se = se,
            expected_n_samples = 999,
            verbose = TRUE
        ),
        "Dimension mismatches found"
    )
})

test_that(".validate_se_dimensions handles expected_assays and expected_n_samples together", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- make_test_se(n_genes = 5, n_samples = 3)
    assay_names <- names(SummarizedExperiment::assays(se))
    
    result <- TSENAT:::.validate_se_dimensions(
        se = se,
        expected_n_samples = 3,
        expected_assays = assay_names,
        verbose = FALSE
    )
    
    expect_true(result$is_valid)
})
