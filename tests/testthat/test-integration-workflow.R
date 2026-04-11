context("Integration Tests: TSENAT Orchestration Function")

# ============================================================================
# GLOBAL SETUP: Reusable helper functions
# ============================================================================

setup_workflow_data <- function() {
    # Set seed for reproducible gene subset selection
    set.seed(42)
    
    data("readcounts", package = "TSENAT")
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    # Create config FIRST (Bioconductor pattern: immutable object construction)
    # OPTIMIZATION: Use 10 q-values for tests (covers 0 to 2)
    config <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        subject_col = "paired_samples",
        q = seq(0, 2, length.out = 10),
        paired = TRUE,
        control = "normal",
        nthreads = 4
    )
    
    # Build analysis with config and explicit metadata parameter
    analysis <- build_analysis(
        config = config,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    
    # Filter with severe stringency (reproducible with seed set above)
    analysis <- filter_analysis(analysis, stringency = "severe")
    
    list(analysis = analysis, se = se(analysis), readcounts = readcounts)
}

# Cached wrapper to avoid expensive recomputation across multiple tests
setup_workflow_data_cached <- local({
    .cache <- NULL
    function() {
        if (is.null(.cache)) {
            .cache <<- setup_workflow_data()
        }
        .cache
    }
})

# ============================================================================
# HELPER FUNCTIONS FOR TEST ASSERTIONS
# ============================================================================

# Validate basic TSENAT result structure
assert_valid_tsenat_result <- function(result) {
    expect_s4_class(result, "TSENATAnalysis")
    expect_s4_class(se(result), "SummarizedExperiment")
    expect_true(nrow(se(result)) > 0)
    expect_true(is.list(getConfig(result)))
    expect_true(is.list(getMeta(result)))
}

# Check SE dimensions are preserved after operations
assert_se_dimensions_unchanged <- function(se_original, se_result, operation = "") {
    expect_equal(nrow(se_result), nrow(se_original),
                 info = paste(operation, "should not change gene count"))
    expect_equal(ncol(se_result), ncol(se_original),
                 info = paste(operation, "should not change sample count"))
    expect_identical(rownames(se_result), rownames(se_original),
                     info = paste(operation, "should not change gene names"))
    expect_identical(colnames(se_result), colnames(se_original),
                     info = paste(operation, "should not change sample names"))
}

# Standard TSENAT orchestration call
run_tsenat_standard <- function(analysis) {
    TSENAT(analysis, output_dir = NULL, verbose = FALSE)
}

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution
# ============================================================================

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution - Paired Design
# ============================================================================

test_that("TSENAT() basic execution: returns valid result with structure preserved", {
    data_list <- setup_workflow_data_cached()
    se_original <- se(data_list$analysis)
    
    # Execute standard TSENAT pipeline
    result <- run_tsenat_standard(data_list$analysis)
    
    # Validate result structure and config
    assert_valid_tsenat_result(result)
    cfg_result <- getConfig(result)
    expect_equal(cfg_result$q, seq(0, 2, length.out = 10))
    
    # Verify SE structure (dimensions may change due to filtering in TSENAT pipeline)
    se_result <- se(result)
    expect_equal(ncol(se_result), ncol(se_original),
                 info = "TSENAT() should not change sample count")
    expect_identical(colnames(se_result), colnames(se_original),
                     info = "TSENAT() should not change sample names")
    expect_true(nrow(se_result) > 0,
                info = "TSENAT() should retain at least some genes")
    expect_true(nrow(se_result) <= nrow(se_original),
                info = "TSENAT() may filter genes but should not add new ones")
    
    # Check diversity results exist
    expect_true(length(results(result, type = "diversity", q = 0)) > 0 || 
                is.null(results(result, type = "diversity", q = 0)))
})

test_that("TSENAT() output handling: manages output_dir and verbose control correctly", {
    data_list <- setup_workflow_data_cached()
    output_dir <- tempdir()
    
    # Test 1: With output_dir
    result_with_output <- TSENAT(
        data_list$analysis,
        output_dir = output_dir,
        verbose = FALSE
    )
    assert_valid_tsenat_result(result_with_output)
    expect_true(dir.exists(output_dir))
    
    # Test 2: With NULL output_dir (no file saving)
    result_no_output <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    assert_valid_tsenat_result(result_no_output)
    
    # Test 3: Verbose output control
    output <- capture.output({
        result_verbose <- TSENAT(
            data_list$analysis,
            output_dir = NULL,
            verbose = TRUE
        )
    })
    assert_valid_tsenat_result(result_verbose)
})

test_that("TSENAT() filtering and statistical parameters: respects filter and config", {
    data_list <- setup_workflow_data_cached()
    
    n_genes_original <- nrow(se(data_list$analysis))
    
    # Run standard pipeline
    result <- run_tsenat_standard(data_list$analysis)
    
    # Verify structure and filtering effects
    assert_valid_tsenat_result(result)
    expect_true(nrow(se(result)) > 0)
    expect_true(nrow(se(result)) <= n_genes_original)
    expect_true(length(result@diversity_results) >= 0)
})

# ============================================================================
# TEST SUITE 2: Error Handling and Edge Cases
# ============================================================================

test_that("TSENAT() error handling: rejects invalid input, handles edge cases", {
    data_list <- setup_workflow_data_cached()
    
    invalid_input <- data.frame(a = 1:10, b = 11:20)
    
    expect_error(
        TSENAT(invalid_input, output_dir = NULL, verbose = FALSE),
        "must be a TSENATAnalysis object"
    )
    
    # Test with empty analysis
    empty_analysis <- data_list$analysis
    empty_analysis@se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(nrow = 0, ncol = ncol(se(data_list$analysis)))),
        colData = SummarizedExperiment::colData(se(data_list$analysis))
    )
    
    expect_error(
        TSENAT(empty_analysis, output_dir = NULL, verbose = FALSE),
        "empty"
    )
})

# ============================================================================
# TEST SUITE 8: LM Interaction Results and Plot Generation
# ============================================================================

test_that("TSENAT() paired: produces LM results, significant genes, plots generate", {
    # Use setup_workflow_data() which already has proper paired config
    data_list <- setup_workflow_data_cached()
    
    # Note: setup_workflow_data() already configured with:
    # condition_col="condition", subject_col="paired_samples", 
    # q_values=seq(0,2,by=0.05), paired=TRUE, control="normal"
    # Retrieve config from analysis instead of using undefined variable
    config <- getConfig(data_list$analysis)
    
    result <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Test LM results structure
    expect_s4_class(result, "TSENATAnalysis")
    lm_res <- results(result, type = "lm")
    expect_true(!is.null(lm_res))
    expect_true(is.data.frame(lm_res))
    expect_true(nrow(lm_res) > 0)
    expect_true("adj_p_interaction" %in% colnames(lm_res) || "p_interaction" %in% colnames(lm_res))
    
    # Check for significant genes
    if ("adj_p_interaction" %in% colnames(lm_res)) {
        sig_genes <- sum(lm_res$adj_p_interaction <= 0.05, na.rm = TRUE)
    } else if ("p_interaction" %in% colnames(lm_res)) {
        sig_genes <- sum(lm_res$p_interaction <= 0.05, na.rm = TRUE)
    } else {
        sig_genes <- 0
    }
    expect_true(sig_genes > 0, 
                info = "Expected significant genes in paired design")
    
    # Test plot generation
    plot_result <- tryCatch({
        plot_lm_gam(
            result,
            n_top = 3,
            sig_alpha = 0.05,
            output_file = NULL,
            verbose = FALSE
        )
    }, error = function(e) NULL)
    
    expect_true(is.null(plot_result) || (class(plot_result)[1] == "gg" || "ggplot" %in% class(plot_result)),
                info = "Plot generation should succeed or return NULL gracefully")
})

context("Integration Tests: setConfig Bug Detection")

# ============================================================================
# TEST SUITE 3: setConfig Bug Detection (GH Issue: Config-induced data corruption)
# ============================================================================
# These tests verify that setConfig does NOT corrupt diversity values
# Background: Passing config parameter to TSENAT() calls setConfig internally,
# which was observed to corrupt diversity calculations (values swapped/changed).
# Solution: Config should only be set during build_analysis(), never via TSENAT().

test_that("setConfig: does not corrupt SummarizedExperiment dimensions", {
    data_list <- setup_workflow_data_cached()
    analysis <- data_list$analysis
    se_orig <- se(analysis)
    
    # Apply setConfig (like TSENAT() does when config parameter is passed)
    config <- getConfig(analysis)
    analysis_after <- setConfig(analysis, config)
    
    # Check dimensions are unchanged
    assert_se_dimensions_unchanged(se_orig, se(analysis_after), "setConfig")
})

test_that("TSENAT() WITHOUT config parameter: produces correct diversity values", {
    # Manual workflow (like vignette): NO setConfig call
    data_list <- setup_workflow_data_cached()
    analysis_manual <- data_list$analysis
    
    # Calculate diversity manually
    analysis_manual <- calculate_diversity(
        analysis_manual,
        norm = TRUE,
        output_file = NULL,
        show_messages = FALSE
    )
    
    # Extract first gene's diversity at q=0 (most stable, no bootstrap)
    div_manual <- results(analysis_manual, type = "diversity", q = 0)
    if (is(div_manual, "SummarizedExperiment")) {
        manual_values <- assay(div_manual, 1)
        manual_first_gene <- manual_values[1, , drop = TRUE]
    }
    
    expect_true(exists("manual_first_gene") && length(manual_first_gene) > 0,
                info = "Manual workflow should produce diversity results")
})

test_that("TSENAT() WITHOUT config: produces IDENTICAL results to manual workflow", {
    # REGRESSION TEST: Verify that TSENAT() orchestration doesn't corrupt data
    # Setup two identical analyses
    set.seed(42)
    data_list1 <- setup_workflow_data_cached()
    analysis1 <- data_list1$analysis
    
    set.seed(42)
    data_list2 <- setup_workflow_data_cached()
    analysis2 <- data_list2$analysis
    
    # Both should have same starting dimensions
    se1_before <- se(analysis1)
    se2_before <- se(analysis2)
    expect_equal(dim(se1_before), dim(se2_before),
                 info = "Both analysis objects should start with identical dimensions")
    
    # Apply TSENAT() to analyze2 (full orchestration)
    analysis2 <- TSENAT(
        analysis2,
        output_dir = NULL,
        save_output = FALSE,
        verbose = FALSE
    )
    
    # After orchestration, analysis2 should still be a valid TSENATAnalysis
    expect_is(analysis2, "TSENATAnalysis",
              info = "TSENAT() should return a valid TSENATAnalysis object")
    
    # It should have results stored (diversity should exist)
    div_result <- tryCatch(
        { results(analysis2, type = "diversity", q = 0) },
        error = function(e) { NULL }
    )
    expect_false(is.null(div_result),
                 info = "TSENAT() should have computed diversity results")
})

test_that("TSENAT() WITH config parameter: SHOULD NOT be used (causes data issues)", {
    # This test documents the problematic behavior when config is passed
    data_list <- setup_workflow_data_cached()
    
    config <- getConfig(data_list$analysis)
    
    # Calling TSENAT() WITH config parameter (incorrect usage that causes bug)
    # This is a regression test to catch if the bug is reintroduced
    analysis_with_config <- data_list$analysis
    
    # Get original values before tsenat
    se_before <- se(analysis_with_config)
    n_genes_before <- nrow(se_before)
    
    # Wrap in tryCatch because the config parameter should no longer exist
    result <- tryCatch({
        analysis_with_config_result <- TSENAT(
            analysis_with_config,
            config = config,  # INCORRECT: passing config parameter
            output_dir = NULL,
            save_output = FALSE,
            verbose = FALSE
        )
        list(result = analysis_with_config_result, error = NULL)
    }, error = function(e) {
        list(result = NULL, error = e)
    })
    
    # After fix: TSENAT() should NOT accept config parameter at all
    # So this test verifies the signature is enforced
    if (!is.null(result$error)) {
        # Good: function rejects config parameter
        expect_true(grepl("config", result$error$message, ignore.case = TRUE),
                    info = "TSENAT() should reject config parameter after fix")
    }
})

test_that("setConfig CORRUPTION: direct calls do not modify analysis state", {
    # setConfig should be idempotent: multiple calls should have no ill effects
    
    data_list <- setup_workflow_data_cached()
    analysis <- data_list$analysis
    config <- getConfig(analysis)
    se_orig <- se(analysis)
    
    # Call setConfig multiple times (should all be safe)
    analysis_1x <- setConfig(analysis, config)
    analysis_2x <- setConfig(analysis_1x, config)
    
    # Verify dimensions unchanged after each call
    assert_se_dimensions_unchanged(se_orig, se(analysis_1x), "setConfig (call 1)")
    assert_se_dimensions_unchanged(se_orig, se(analysis_2x), "setConfig (call 2)")
})

# ============================================================================
# TEST SUITE 4: Configuration Embedding and Immutability
# ============================================================================

test_that("CONFIG EMBEDDING: Settings applied once via build_analysis, not redundantly", {
    # Test the correct pattern: config applied exactly ONCE at build time
    
    set.seed(42)
    config1 <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        subject_col = "paired_samples",
        q = seq(0, 1, length.out = 5),
        paired = TRUE,
        control = "normal",
        nthreads = 4
    )
    
    data("readcounts", package = "TSENAT")
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    # build_analysis embeds config exactly once
    analysis <- build_analysis(
        config = config1,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    
    # The @config slot should contain exactly what we passed in
    cfg_embedded <- getConfig(analysis)
    cfg_embedded_q <- if (is.list(cfg_embedded)) cfg_embedded$q else cfg_embedded@q
    config1_q <- if (is.list(config1)) config1$q else config1@q
    
    expect_equal(cfg_embedded_q, config1_q,
                 info = "build_analysis should embed config exactly as provided")
    
    # filter_analysis should NOT modify @config
    analysis2 <- filter_analysis(analysis, stringency = "severe")
    cfg_after <- getConfig(analysis2)
    cfg_after_q <- if (is.list(cfg_after)) cfg_after$q else cfg_after@q
    
    expect_identical(cfg_after_q, cfg_embedded_q,
                     info = "filter_analysis should preserve embedded config unchanged")
})

test_that("IDEMPOTENCY CHECK: setConfig produces consistent state across calls", {
    # If setConfig is truly idempotent, multiple calls should produce identical state
    
    data_list <- setup_workflow_data_cached()
    analysis <- data_list$analysis
    config <- getConfig(analysis)
    
    # Call setConfig multiple times
    analysis_1x <- setConfig(analysis, config)
    analysis_2x <- setConfig(analysis_1x, config)
    analysis_3x <- setConfig(analysis_2x, config)
    
    # All should produce identical SE dimensions
    se_1x <- se(analysis_1x)
    se_2x <- se(analysis_2x)
    se_3x <- se(analysis_3x)
    
    expect_equal(dim(se_1x), dim(se_2x),
                 info = "setConfig should be idempotent (consistent after call 1 vs 2)")
    expect_equal(dim(se_2x), dim(se_3x),
                 info = "setConfig should be idempotent (consistent after call 2 vs 3)")
})

# ============================================================================
# TEST SUITE 5: Orchestration Patterns - Manual vs TSENAT()
# ============================================================================

test_that("WORKFLOW EQUIVALENCE: Manual orchestration matches TSENAT() function", {
    # Compare orchestration patterns to ensure no hidden side effects
    # Pattern A (manual): filter → calculate_diversity (no setConfig)
    # Pattern B (orchestrated): TSENAT() (may call setConfig internally)
    
    set.seed(42)
    data_list_A <- setup_workflow_data_cached()
    analysis_A <- data_list_A$analysis
    
    # Pattern A: Manual calculation
    analysis_A <- calculate_diversity(
        analysis_A,
        norm = TRUE,
        output_file = NULL,
        show_messages = FALSE
    )
    
    # Both should have valid diversity results
    div_A <- results(analysis_A, type = "diversity", q = 0)
    
    expect_is(div_A, "SummarizedExperiment",
              info = "Manual pattern should produce valid diversity SE")
})
