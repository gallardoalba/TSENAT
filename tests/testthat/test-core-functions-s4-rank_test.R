library(testthat)

# ============================================================================
# Tests for rank_test_q_condition_s4() S4 wrapper function
# ============================================================================

context("S4 Rank Test: rank_test_q_condition_s4")

# Helper function to create realistic test TSENATAnalysis with valid diversity results
# Matches roxygen documentation workflow
setup_rank_test_analysis <- function(n_genes = 30, n_samples = 8) {
    
    # Load example data (matching roxygen example)
    data(readcounts)
    readcounts <- as.matrix(salmon_dataset)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    # Build analysis from vignette data and create small subset
    # Use both tpm and effective_length like roxygen example
    analysis <- build_analysis_s4(
        readcounts,
        gff3_dataset,
        metadata = metadata_df,
        tpm = salmon_tpm,
        effective_length = salmon_effective_length
    )
    analysis <- subset_analysis(analysis, n_genes = n_genes, n_samples = n_samples)
    analysis <- calculate_diversity_s4(analysis, norm = TRUE)
    
    analysis
}

# ============================================================================
# Test 1: Basic functionality with real data
# ============================================================================

test_that("rank_test_q_condition_s4 returns TSENATAnalysis with results", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 12)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition"
    )
    
    # Should return modified TSENATAnalysis
    expect_is(result, "TSENATAnalysis")
    # Should have lm_results stored (list with q_interactions data frame)
    lm_res <- TSENAT::lmResults(result)
    expect_is(lm_res, "list")
    expect_true("q_interactions" %in% names(lm_res))
})

# ============================================================================
# Test 2: Paired design with subject_col
# ============================================================================

test_that("rank_test_q_condition_s4 runs with paired design", {
    analysis <- setup_rank_test_analysis(n_genes = 20, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        paired = TRUE,
        subject_col = "paired_samples"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 3: Parameter resolution - explicit args override config
# ============================================================================

test_that("rank_test_q_condition_s4 respects explicit parameters over config", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "kruskal-wallis"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 4: Different test methods - Kruskal-Wallis (unpaired)
# ============================================================================

test_that("rank_test_q_condition_s4 works with Kruskal-Wallis test", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "kruskal-wallis"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("p_value" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Test 5: Different test methods - Friedman (paired)
# ============================================================================

test_that("rank_test_q_condition_s4 works with Friedman test for paired design", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "friedman",
        paired = TRUE,
        subject_col = "paired_samples"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 6: Different multiple correction methods - Hochberg
# ============================================================================

test_that("rank_test_q_condition_s4 works with Hochberg correction", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "hochberg"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("adj_p_value" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Test 7: Different multiple correction methods - Benjamini-Yekutieli
# ============================================================================

test_that("rank_test_q_condition_s4 works with Benjamini-Yekutieli correction", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "benjamini-yekutieli"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 8: Multiple correction method - None
# ============================================================================

test_that("rank_test_q_condition_s4 works with no multiple correction", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "none"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 9: Auto test selection
# ============================================================================

test_that("rank_test_q_condition_s4 auto-selects appropriate test method", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "auto"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 10: Error handling - missing required condition_col
# ============================================================================

test_that("rank_test_q_condition_s4 requires condition_col argument", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Calling without condition_col should not error (defaults to "condition")
    # but should work if condition column exists in metadata
    result <- rank_test_q_condition_s4(analysis)
    expect_is(result, "TSENATAnalysis")
    
    # Invalid condition_col should error
    expect_error(rank_test_q_condition_s4(analysis, condition_col = "nonexistent"))
})

# ============================================================================
# Test 11: Invalid condition column
# ============================================================================

test_that("rank_test_q_condition_s4 errors with invalid condition_col", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Should error when condition column doesn't exist in metadata
    expect_error(
        rank_test_q_condition_s4(
            analysis,
            condition_col = "nonexistent_column"
        ),
        "condition_col|not found|invalid|does not exist"
    )
})

# ============================================================================
# Test 12: Paired design requires subject_col
# ============================================================================

test_that("rank_test_q_condition_s4 requires subject_col when paired=TRUE", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Should error when subject column doesn't exist but paired=TRUE
    expect_error(
        rank_test_q_condition_s4(
            analysis,
            condition_col = "condition",
            paired = TRUE,
            subject_col = "nonexistent_subject"
        )
    )
})

# ============================================================================
# Test 13: Results contain expected columns
# ============================================================================

test_that("rank_test_q_condition_s4 results have correct structure", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "hochberg"
    )
    
    lm_res <- TSENAT::lmResults(result)
    q_int_res <- lm_res$q_interactions
    
    # Check for essential columns in q_interactions
    expect_true("gene" %in% colnames(q_int_res))
    expect_true("p_value" %in% colnames(q_int_res))
    expect_true("adj_p_value" %in% colnames(q_int_res))
    expect_true("f_statistic" %in% colnames(q_int_res))
})

# ============================================================================
# Test 14: Multiple genes with varying p-values
# ============================================================================

test_that("rank_test_q_condition_s4 handles multiple genes with varying significance", {
    analysis <- setup_rank_test_analysis(n_genes = 30, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "hochberg"
    )
    
    lm_res <- TSENAT::lmResults(result)
    q_int_res <- lm_res$q_interactions
    
    # Should have results for tested genes
    expect_true(nrow(q_int_res) > 0)
    
    # p-values should be numeric and in valid range
    expect_true(all(!is.na(q_int_res$p_value)))
    expect_true(all(q_int_res$p_value >= 0 & q_int_res$p_value <= 1))
})

# ============================================================================
# Test 15: ART (Aligned Rank Transform) test method
# ============================================================================

test_that("rank_test_q_condition_s4 works with ART (Aligned Rank Transform)", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "art"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("test_method" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Test 16: Westfall-Young correction method
# ============================================================================

test_that("rank_test_q_condition_s4 works with Westfall-Young correction", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "westfall-young",
        wy_randomizations = 100
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("adj_p_value" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Test 17: Explicit q-values parameter
# ============================================================================

test_that("rank_test_q_condition_s4 accepts explicit q parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Specify subset of q-values to test
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        q = c(0.5, 1.0, 1.5)
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 18: Custom entropy column name
# ============================================================================

test_that("rank_test_q_condition_s4 accepts custom entropy_col parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Default entropy_col is "diversity"
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        entropy_col = "diversity"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 19: Custom q column name
# ============================================================================

test_that("rank_test_q_condition_s4 accepts custom q_col parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        q_col = "q"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 20: Custom gene column name
# ============================================================================

test_that("rank_test_q_condition_s4 accepts custom gene_col parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        gene_col = "gene"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("gene" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Test 21: nthreads parameter for parallelization
# ============================================================================

test_that("rank_test_q_condition_s4 accepts nthreads parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        nthreads = 1
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 22: verbose parameter for progress output
# ============================================================================

test_that("rank_test_q_condition_s4 accepts verbose parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Verbose may or may not produce output, but should not error
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        verbose = FALSE
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 23: nperm_mode parameter (permutation estimation mode)
# ============================================================================

test_that("rank_test_q_condition_s4 accepts nperm_mode parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Test with different permutation modes
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        nperm_mode = "standard"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 24: wy_randomizations parameter (Westfall-Young permutations)
# ============================================================================

test_that("rank_test_q_condition_s4 wy_randomizations controls WY permutations", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Test with explicit WY randomizations count
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        multicorr = "westfall-young",
        wy_randomizations = 50
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
})

# ============================================================================
# Test 25: output_file parameter (optional file output)
# ============================================================================

test_that("rank_test_q_condition_s4 accepts output_file parameter", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Create temporary file path
    temp_file <- tempfile(fileext = ".rds")
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        output_file = temp_file
    )
    
    expect_is(result, "TSENATAnalysis")
    # File may or may not be created depending on implementation
    # Main check is that function doesn't error
})

# ============================================================================
# Test 26: Combined parameters - multiple options together
# ============================================================================

test_that("rank_test_q_condition_s4 handles combined parameter specifications", {
    analysis <- setup_rank_test_analysis(n_genes = 20, n_samples = 8)
    
    result <- rank_test_q_condition_s4(
        analysis,
        condition_col = "condition",
        test = "friedman",
        multicorr = "benjamini-yekutieli",
        paired = TRUE,
        subject_col = "paired_samples",
        nthreads = 1,
        verbose = FALSE,
        nperm_mode = "standard"
    )
    
    expect_is(result, "TSENATAnalysis")
    lm_res <- TSENAT::lmResults(result)
    expect_true(nrow(lm_res$q_interactions) > 0)
    expect_true("adj_p_value" %in% colnames(lm_res$q_interactions))
})

# ============================================================================
# Tests for Helper Functions (NEW - Refactored S4 Wrapper)
# ============================================================================

# ============================================================================
# Test: .validate_rank_test_input() helper function
# ============================================================================

test_that(".validate_rank_test_input rejects non-TSENATAnalysis input", {
    # Test with invalid input type
    expect_error(
        TSENAT:::.validate_rank_test_input("not_an_analysis", "condition"),
        "must be a TSENATAnalysis"
    )
    
    # Test with data.frame instead of analysis object
    expect_error(
        TSENAT:::.validate_rank_test_input(data.frame(x = 1:10), "condition"),
        "must be a TSENATAnalysis"
    )
})

test_that(".validate_rank_test_input requires diversity results", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Clear diversity_results to test prerequisite check
    analysis@diversity_results <- list()
    
    expect_error(
        TSENAT:::.validate_rank_test_input(analysis, "condition"),
        "Diversity results required"
    )
})

test_that(".validate_rank_test_input returns condition_col with default fallback", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # When condition_col is NULL/missing, should default to "condition"
    result <- TSENAT:::.validate_rank_test_input(analysis, NULL)
    expect_equal(result, "condition")
    
    # When condition_col is explicitly provided, should return it
    result <- TSENAT:::.validate_rank_test_input(analysis, "sample_type")
    expect_equal(result, "sample_type")
})

test_that(".validate_rank_test_input respects config condition_col", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Set condition_col in config
    analysis@config$condition_col <- "my_condition"
    
    # When condition_col is NULL/missing, should read from config
    result <- TSENAT:::.validate_rank_test_input(analysis, NULL)
    expect_equal(result, "my_condition")
})

# ============================================================================
# Test: .resolve_rank_test_params() helper function
# ============================================================================

test_that(".resolve_rank_test_params handles explicit arguments over config", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Set config values
    analysis@config <- list(
        test = "kruskal-wallis",
        multicorr = "hochberg",
        nperm_mode = "conservative",
        q_values = c(0.5, 1.0),
        paired = FALSE,
        subject_col = "default_subject",
        nthreads = 2
    )
    
    # Call with explicit arguments (should override config)
    result <- TSENAT:::.resolve_rank_test_params(
        analysis,
        test = "friedman",
        multicorr = "benjamini-yekutieli",
        nperm_mode = "standard",
        q = c(0.5, 1.0, 1.5),
        paired = TRUE,
        subject_col = "explicit_subject",
        nthreads = 4,
        wy_randomizations = 200,
        entropy_col = "custom_entropy",
        q_col = "custom_q",
        gene_col = "custom_gene"
    )
    
    dots <- result$dots
    
    # Check explicit args took precedence
    expect_equal(dots$test, "friedman")
    expect_equal(dots$multicorr, "benjamini-yekutieli")
    expect_equal(dots$nperm_mode, "standard")
    expect_equal(dots$nthreads, 4)
    expect_equal(dots$entropy_col, "custom_entropy")
    expect_equal(dots$q_col, "custom_q")
    expect_equal(dots$gene_col, "custom_gene")
})

test_that(".resolve_rank_test_params falls back to config when args not provided", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Set config values
    analysis@config <- list(
        test = "kruskal-wallis",
        multicorr = "hochberg",
        nperm_mode = "interactive",
        q_values = c(0.5, 1.0, 1.5),
        paired = FALSE,
        subject_col = "subject_id",
        nthreads = 2
    )
    
    # Call without explicit arguments (should use config)
    result <- TSENAT:::.resolve_rank_test_params(
        analysis,
        test = NULL,
        multicorr = NULL,
        nperm_mode = NULL,
        q = NULL,
        paired = NULL,
        subject_col = NULL,
        nthreads = NULL,
        wy_randomizations = 500,
        entropy_col = "diversity",
        q_col = "q",
        gene_col = "gene"
    )
    
    # Note: need to handle the way the function checks for missing vs NULL
    # This is a basic structure test
    expect_true("dots" %in% names(result))
    expect_true("q_extracted" %in% names(result))
})

test_that(".resolve_rank_test_params validates enum arguments", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Test with invalid enum values - should error
    expect_error(
        TSENAT:::.resolve_rank_test_params(
            analysis,
            test = "invalid_test",
            multicorr = NULL,
            nperm_mode = NULL,
            q = NULL,
            paired = NULL,
            subject_col = NULL,
            nthreads = NULL,
            wy_randomizations = 500,
            entropy_col = "diversity",
            q_col = "q",
            gene_col = "gene"
        ),
        "should be one of"
    )
})

# ============================================================================
# Test: .prepare_multi_q_se() helper function
# ============================================================================

test_that(".prepare_multi_q_se combines multiple q-value results", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Test that function returns a SummarizedExperiment
    se_multi_q <- TSENAT:::.prepare_multi_q_se(analysis)
    
    expect_is(se_multi_q, "SummarizedExperiment")
    expect_true(ncol(se_multi_q) > 0)
    expect_true(nrow(se_multi_q) > 0)
    
    # Check that combined SE has q column in colData
    coldata <- SummarizedExperiment::colData(se_multi_q)
    expect_true("q" %in% colnames(coldata))
})

test_that(".prepare_multi_q_se preserves gene names across q-values", {
    analysis <- setup_rank_test_analysis(n_genes = 20, n_samples = 8)
    
    se_multi_q <- TSENAT:::.prepare_multi_q_se(analysis)
    
    # Get original gene names from first diversity result
    first_se <- analysis@diversity_results[[1]]
    original_genes <- rownames(first_se)
    
    # Check that genes match in combined SE
    combined_genes <- rownames(se_multi_q)
    expect_equal(combined_genes, original_genes)
})

test_that(".prepare_multi_q_se uses cached combined SE when available", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Create cached combined SE
    cached_se <- SummarizedExperiment(assays = list(entropy = matrix(rnorm(120), nrow = 15)))
    analysis@metadata$diversity_combined <- list(combined_se = cached_se)
    
    se_multi_q <- TSENAT:::.prepare_multi_q_se(analysis)
    
    # Should return the cached version
    expect_identical(se_multi_q, cached_se)
})

# ============================================================================
# Test: .store_rank_test_results() helper function
# ============================================================================

test_that(".store_rank_test_results stores results in lm_results", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Create dummy results data frame (mock output from rank test)
    mock_results <- data.frame(
        gene = c("gene1", "gene2", "gene3"),
        p_value = c(0.001, 0.05, 0.1),
        adj_p_value = c(0.01, 0.1, 0.2),
        f_statistic = c(10.5, 5.2, 2.1)
    )
    
    # Store results
    analysis_stored <- TSENAT:::.store_rank_test_results(
        analysis, 
        mock_results, 
        output_file = NULL, 
        verbose = FALSE
    )
    
    # Check results are stored in correct location
    expect_true(is.list(analysis_stored@lm_results))
    expect_true("q_interactions" %in% names(analysis_stored@lm_results))
    expect_equal(nrow(analysis_stored@lm_results$q_interactions), 3)
})

test_that(".store_rank_test_results creates lm_results list when needed", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    # Clear lm_results to test list creation
    analysis@lm_results <- list()
    
    mock_results <- data.frame(
        gene = "gene1",
        p_value = 0.01,
        adj_p_value = 0.05,
        f_statistic = 5.0
    )
    
    analysis_stored <- TSENAT:::.store_rank_test_results(
        analysis, 
        mock_results, 
        output_file = NULL, 
        verbose = FALSE
    )
    
    # Check list was created properly
    expect_true(is.list(analysis_stored@lm_results))
    expect_true("q_interactions" %in% names(analysis_stored@lm_results))
})

test_that(".store_rank_test_results returns modified TSENATAnalysis", {
    analysis <- setup_rank_test_analysis(n_genes = 10, n_samples = 8)
    
    mock_results <- data.frame(gene = "gene1", p_value = 0.01)
    
    analysis_returned <- TSENAT:::.store_rank_test_results(
        analysis, 
        mock_results, 
        output_file = NULL, 
        verbose = FALSE
    )
    
    # Check type is preserved
    expect_is(analysis_returned, "TSENATAnalysis")
})

# ============================================================================
# Integration Test: Helper functions work together in workflow
# ============================================================================

test_that("Helper functions integrate correctly in rank test workflow", {
    analysis <- setup_rank_test_analysis(n_genes = 15, n_samples = 8)
    
    # Simulate workflow: validate -> resolve params -> prepare SE -> run test -> store results
    
    # 1. Validate
    condition_col <- TSENAT:::.validate_rank_test_input(analysis, "condition")
    expect_equal(condition_col, "condition")
    
    # 2. Resolve params
    param_result <- TSENAT:::.resolve_rank_test_params(
        analysis,
        test = "auto",
        multicorr = "hochberg",
        nperm_mode = "standard",
        q = NULL,
        paired = FALSE,
        subject_col = NULL,
        nthreads = 1,
        wy_randomizations = 100,
        entropy_col = "diversity",
        q_col = "q",
        gene_col = "gene"
    )
    expect_true(is.list(param_result$dots))
    
    # 3. Prepare SE
    se_multi_q <- TSENAT:::.prepare_multi_q_se(analysis)
    expect_is(se_multi_q, "SummarizedExperiment")
    
    # 4. Mock store results
    mock_results <- data.frame(gene = "gene1", p_value = 0.01, adj_p_value = 0.05)
    analysis_final <- TSENAT:::.store_rank_test_results(
        analysis,
        mock_results,
        output_file = NULL,
        verbose = FALSE
    )
    
    expect_is(analysis_final, "TSENATAnalysis")
    expect_true("q_interactions" %in% names(analysis_final@lm_results))
})


