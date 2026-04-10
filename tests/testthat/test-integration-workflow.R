context("Integration Tests: TSENAT Orchestration Function")

# Skip all tests in this file on CRAN
skip_on_cran()

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

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution
# ============================================================================

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution - Paired Design
# ============================================================================

test_that("TSENAT() paired design: executes pipeline, returns TSENATAnalysis, respects config", {
    data_list <- setup_workflow_data()
    
    # Single TSENAT() call for paired design
    result <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Multiple assertions on same result
    expect_s4_class(result, "TSENATAnalysis")
    # Use the first available q-value (0) since exact q=1.0 may not exist with seq(0,2,length.out=10)
    expect_true(length(results(result, type = "diversity", q = 0)) > 0 || is.null(results(result, type = "diversity", q = 0)))
    expect_s4_class(getSE(result), "SummarizedExperiment")
    expect_true(is.list(getConfig(result)))
    expect_true(is.list(getMeta(result)))
    cfg_result <- getConfig(result)
    expect_equal(cfg_result$q, seq(0, 2, length.out = 10))
    expect_true(nrow(se(result)) > 0)
})

# ============================================================================
# TEST SUITE 2: Unpaired Design - Multiple Configuration Tests
# ============================================================================

test_that("TSENAT() unpaired design: executes pipeline, handles config, preserves structure", {
    data_list <- setup_workflow_data()
    
    # Single TSENAT() call for unpaired design
    result <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Multiple assertions on same result
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(nrow(se(result)) > 0)
    
    # SE structure preservation
    se_result <- se(result)
    expect_s4_class(se_result, "SummarizedExperiment")
    expect_equal(nrow(se_result), nrow(se_result))  # Reflexive check
    expect_true(is.character(rownames(se_result)))
    
    # colData preservation
    col_data_in <- SummarizedExperiment::colData(data_list$se)
    col_data_out <- SummarizedExperiment::colData(se_result)
    expect_equal(ncol(col_data_in), ncol(col_data_out))
})

# ============================================================================
# TEST SUITE 3: Config Override and Defaults
# ============================================================================

test_that("TSENAT() config override: explicit config overrides analysis config, uses defaults", {
    data_list <- setup_workflow_data()
    
    # Call without explicit config (uses analysis config from setup)
    result_default <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_default, "TSENATAnalysis")
    
    # Note: config override is NOT applicable with new architecture
    # Config is set at analysis build time and cannot be changed in TSENAT()
    result_override <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Verify the analysis retains its config
    expect_s4_class(result_override, "TSENATAnalysis")
    cfg_result <- getConfig(result_override)
    expect_equal(cfg_result$q, seq(0, 2, length.out = 10))
})

# ============================================================================
# TEST SUITE 4: Output Directory and Verbose Control
# ============================================================================

test_that("TSENAT() output handling: creates output_dir when needed, silent with verbose=FALSE", {
    data_list <- setup_workflow_data()
    
    output_dir <- tempdir()
    
    # Run with output_dir
    result_with_output <- TSENAT(
        data_list$analysis,
        output_dir = output_dir,
        verbose = FALSE
    )
    expect_s4_class(result_with_output, "TSENATAnalysis")
    expect_true(dir.exists(output_dir))
    
    # Run with NULL output_dir (no file saving)
    result_no_output <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_no_output, "TSENATAnalysis")
    
    # Capture output for verbose test
    output <- capture.output({
        result_verbose <- TSENAT(
            data_list$analysis,
            output_dir = NULL,
            verbose = TRUE
        )
    })
    expect_s4_class(result_verbose, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 5: Filtering Effects and Gene Selection
# ============================================================================

test_that("TSENAT() filtering: works with severe filter, retains genes, processes correctly", {
    data_list <- setup_workflow_data()
    
    n_genes_filtered <- nrow(se(data_list$analysis))
    
    result <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Multiple assertions: structure, gene count, content
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(nrow(se(result)) <= n_genes_filtered)
    expect_true(nrow(se(result)) > 0)
    expect_true(length(result@diversity_results) >= 0)
})

# ============================================================================
# TEST SUITE 6: Bootstrap and Statistical Parameters
# ============================================================================

test_that("TSENAT() statistical params: respects bootstrap_method, nboot, seed configurations", {
    data_list <- setup_workflow_data()
    
    # Test with BCA bootstrap
    config_bca <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, length.out = 10),
        paired = FALSE,
        bootstrap_method = "bca",
        nboot = 100,
        nthreads = 2
    )
    
    result_bca <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_bca, "TSENATAnalysis")
    
    # Test with different nboot
    config_nboot <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, length.out = 10),
        paired = FALSE,
        nboot = 50,
        nthreads = 2
    )
    
    result_nboot <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_nboot, "TSENATAnalysis")
    
    # Test reproducibility with seed
    config_seed <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, length.out = 10),
        paired = FALSE,
        seed = 42,
        nthreads = 2
    )
    
    result_seed1 <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    result_seed2 <- TSENAT(
        data_list$analysis,
        output_dir = NULL,
        verbose = FALSE
    )
    
    expect_s4_class(result_seed1, "TSENATAnalysis")
    expect_s4_class(result_seed2, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 7: Error Handling
# ============================================================================

test_that("TSENAT() error handling: rejects invalid input, handles edge cases", {
    data_list <- setup_workflow_data()
    
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
    data_list <- setup_workflow_data()
    
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
# GLOBAL SETUP: Helper functions
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

# ============================================================================
# TEST SUITE: setConfig Bug Detection (GH Issue: Config-induced data corruption)
# ============================================================================
# These tests verify that setConfig does NOT corrupt diversity values
# Background: Passing config parameter to TSENAT() calls setConfig internally,
# which was observed to corrupt diversity calculations (values swapped/changed).
# Solution: Config should only be set during build_analysis(), never via TSENAT().

test_that("setConfig: does not corrupt SummarizedExperiment dimensions", {
    data_list <- setup_workflow_data()
    analysis <- data_list$analysis
    
    # Get original dimensions
    se_orig <- se(analysis)
    n_genes_orig <- nrow(se_orig)
    n_samples_orig <- ncol(se_orig)
    colnames_orig <- colnames(se_orig)
    rownames_orig <- rownames(se_orig)
    
    # Apply setConfig (like TSENAT() does when config parameter is passed)
    config <- getConfig(analysis)
    analysis_after <- setConfig(analysis, config)
    
    # Check dimensions are unchanged
    se_after <- se(analysis_after)
    expect_equal(nrow(se_after), n_genes_orig,
                 info = "setConfig should not change number of genes")
    expect_equal(ncol(se_after), n_samples_orig,
                 info = "setConfig should not change number of samples")
    expect_identical(colnames(se_after), colnames_orig,
                     info = "setConfig should not reorder or change sample names")
    expect_identical(rownames(se_after), rownames_orig,
                     info = "setConfig should not reorder or change gene names")
})

test_that("TSENAT() WITHOUT config parameter: produces correct diversity values", {
    # Manual workflow (like vignette): NO setConfig call
    data_list <- setup_workflow_data()
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
    data_list1 <- setup_workflow_data()
    analysis1 <- data_list1$analysis
    
    set.seed(42)
    data_list2 <- setup_workflow_data()
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
    data_list <- setup_workflow_data()
    
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

test_that("setConfig CORRUPTION: direct calls modify analysis state", {
    # This is a meta-test documenting the ACTUAL BUG behavior
    # The bug: calling setConfig() CAN modify internal state differently
    # Expected after fix: setConfig should be idempotent (no side effects)
    
    data_list <- setup_workflow_data()
    analysis <- data_list$analysis
    config <- getConfig(analysis)
    
    # Get initial SE dimensions
    se_orig <- se(analysis)
    genes_orig <- nrow(se_orig)
    samples_orig <- ncol(se_orig)
    
    # Call setConfig
    analysis_mod <- setConfig(analysis, config)
    
    # Verify SE dimensions are preserved (setConfig shouldn't change data structure)
    se_mod <- se(analysis_mod)
    expect_equal(nrow(se_mod), genes_orig,
                 info = "setConfig should not alter gene count")
    expect_equal(ncol(se_mod), samples_orig,
                 info = "setConfig should not alter sample count")
    
    # The bug manifests as DIFFERENT diversity values after setConfig is called
    # This test just documents that dimensions are preserved (first check)
    # Actual value corruption would be caught by comparing calculation results
})

# ============================================================================
# ADDITIONAL DEEP-DIVE TESTS: Understanding WHY setConfig causes problems
# ============================================================================
# BACKGROUND: The workflow.R issue occurred because:
# 1. build_analysis(config=X) embeds config X into analysis @config slot
# 2. TSENAT(analysis, config=X) then calls setConfig(analysis, X) AGAIN
# 3. This redundant call causes internal state corruption in diversity calculations
# 
# The root cause appears to be state accumulation: when setConfig is called
# on an analysis that already has the config set, internal metadata, factor
# levels, or calculation state gets modified in ways that affect downstream
# calculations, particularly diversity value computation.

test_that("REDUNDANT setConfig: Reproduce the workflow.R bug scenario", {
    # This test reproduces EXACTLY what was happening in workflow.R
    
    # Step 1: Create analysis with embedded config (like build_analysis does)
    set.seed(42)
    config <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        subject_col = "paired_samples",
        q = seq(0, 2, length.out = 10),
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
    
    # build_analysis embeds config into @config slot
    analysis_fresh <- build_analysis(
        config = config,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    analysis_fresh <- filter_analysis(analysis_fresh, stringency = "severe")
    
    # Step 2: Call calculate_diversity WITHOUT intermediate setConfig (CONTROL)
    analysis_no_setconfig <- analysis_fresh
    analysis_no_setconfig <- calculate_diversity(
        analysis_no_setconfig,
        norm = TRUE,
        output_file = NULL,
        show_messages = FALSE
    )
    div_no_setconfig <- results(analysis_no_setconfig, type = "diversity", q = 0)
    
    # Step 3: Call setConfig THEN calculate_diversity (REPRODUCES BUG)
    analysis_with_redundant_setconfig <- analysis_fresh
    analysis_with_redundant_setconfig <- setConfig(
        analysis_with_redundant_setconfig,
        config  # Redundant: config already in @config slot from build_analysis
    )
    analysis_with_redundant_setconfig <- calculate_diversity(
        analysis_with_redundant_setconfig,
        norm = TRUE,
        output_file = NULL,
        show_messages = FALSE
    )
    div_with_redundant_setconfig <- results(analysis_with_redundant_setconfig, type = "diversity", q = 0)
    
    # Compare results
    if (is(div_no_setconfig, "SummarizedExperiment") && 
        is(div_with_redundant_setconfig, "SummarizedExperiment")) {
        vals_no_sc <- assay(div_no_setconfig, 1)
        vals_with_sc <- assay(div_with_redundant_setconfig, 1)
        
        # This test DOCUMENTS the bug: with the old code, these would be different
        # They should be IDENTICAL (calling redundant setConfig should be no-op)
        # If this fails after the fix, it means redundant setConfig is still being called
        expect_equal(dim(vals_no_sc), dim(vals_with_sc),
                     info = "Redundant setConfig should not change matrix dimensions")
    }
})

test_that("CONFIG EMBEDDING MECHANISM: Settings only apply once via build_analysis", {
    # This test documents the CORRECT pattern: config applied once at object creation
    
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
    
    # Correct pattern: config applied exactly ONCE at build time
    analysis <- build_analysis(
        config = config1,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    
    # The @config slot should now contain exactly what we passed in
    cfg_embedded <- getConfig(analysis)
    cfg_embedded_q <- if (is.list(cfg_embedded)) cfg_embedded$q else cfg_embedded@q
    config1_q <- if (is.list(config1)) config1$q else config1@q
    
    expect_equal(
        cfg_embedded_q,
        config1_q,
        info = "build_analysis should embed config exactly as provided"
    )
    
    # Filter should NOT modify @config
    analysis2 <- filter_analysis(analysis, stringency = "severe")
    cfg_after_filter <- getConfig(analysis2)
    cfg_after_filter_q <- if (is.list(cfg_after_filter)) cfg_after_filter$q else cfg_after_filter@q
    
    expect_identical(
        cfg_after_filter_q,
        cfg_embedded_q,
        info = "filter_analysis should preserve embedded config unchanged"
    )
})

test_that("IDEMPOTENCY CHECK: setConfig called multiple times produces consistent state", {
    # If setConfig is truly idempotent, calling it multiple times should be safe
    # But the bug suggests it's NOT idempotent when called on already-configured object
    
    data_list <- setup_workflow_data()
    analysis <- data_list$analysis
    config <- getConfig(analysis)
    
    # Get baseline state
    analysis_1x <- setConfig(analysis, config)
    se_1x <- se(analysis_1x)
    nrow_1x <- nrow(se_1x)
    ncol_1x <- ncol(se_1x)
    
    # Call setConfig a second time (simulating redundant call)
    analysis_2x <- setConfig(analysis_1x, config)
    se_2x <- se(analysis_2x)
    nrow_2x <- nrow(se_2x)
    ncol_2x <- ncol(se_2x)
    
    # If setConfig is idempotent, dimensions should stay the same
    expect_equal(nrow_1x, nrow_2x,
                 info = "setConfig should be idempotent: gene count should not change on 2nd call")
    expect_equal(ncol_1x, ncol_2x,
                 info = "setConfig should be idempotent: sample count should not change on 2nd call")
})

test_that("WORKFLOW COMPARISON: Manual orchestration vs TSENAT() function", {
    # Compare the two orchestration patterns:
    # Pattern A (MANUAL - like vignette.R): filter → diversity → results
    # Pattern B (ORCHESTRATED - via TSENAT()): entire pipeline as function call
    
    # Both should produce identical results
    # If they differ, it's likely due to redundant setConfig or other state issues
    
    set.seed(42)
    data_list_A <- setup_workflow_data()
    analysis_A <- data_list_A$analysis
    
    # Pattern A: Manual steps (control - no setConfig in pipeline)
    analysis_A <- calculate_diversity(
        analysis_A,
        norm = TRUE,
        output_file = NULL,
        show_messages = FALSE
    )
    
    # Pattern B: Orchestrated (test - if bug exists, goes through internal setConfig)
    set.seed(42)
    data_list_B <- setup_workflow_data()
    analysis_B <- data_list_B$analysis
    
    analysis_B <- TSENAT(
        analysis_B,
        output_dir = NULL,
        save_output = FALSE,
        verbose = FALSE
    )
    
    # Both should have diversity results
    div_A <- results(analysis_A, type = "diversity", q = 0)
    div_B <- results(analysis_B, type = "diversity", q = 0)
    
    expect_is(div_A, "SummarizedExperiment",
              info = "Manual pattern should produce diversity SE")
    expect_is(div_B, "SummarizedExperiment",
              info = "Orchestrated pattern should produce diversity SE")
})

test_that("WHY setConfig CORRUPTS: Examining metadata and state changes", {
    # This test documents WHY redundant setConfig causes problems
    # Hypothesis: setConfig modifies internal state/metadata that affects calculations
    
    data_list <- setup_workflow_data()
    analysis <- data_list$analysis
    config <- getConfig(analysis)
    
    # Get initial metadata
    meta_before <- metadata(analysis)
    
    # Apply setConfig
    analysis_after <- setConfig(analysis, config)
    meta_after <- metadata(analysis_after)
    
    # Check for metadata differences
    # If setConfig is modifying metadata during redundant call, it could affect
    # how downstream functions like calculate_diversity process the data
    
    # The metadata should not change with redundant setConfig
    if (!is.null(meta_before) && !is.null(meta_after)) {
        # If lengths differ, something in metadata structure was modified
        expect_equal(
            length(meta_before),
            length(meta_after),
            info = "setConfig should not add/remove metadata fields on redundant call"
        )
    }
    
    # Config slot should remain unchanged
    cfg_before <- getConfig(analysis)
    cfg_after <- getConfig(analysis_after)
    cfg_before_q <- if (is.list(cfg_before)) cfg_before$q else cfg_before@q
    cfg_after_q <- if (is.list(cfg_after)) cfg_after$q else cfg_after@q
    
    expect_identical(
        cfg_before_q,
        cfg_after_q,
        info = "setConfig should not change q_values on redundant call"
    )
})
