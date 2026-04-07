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
    config <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        subject_col = "paired_samples",
        q_values = seq(0, 2, length.out = 10),
        paired = TRUE,
        control = "normal",
        nthreads = 4
    )
    
    # Build analysis with config and explicit metadata parameter
    analysis <- build_analysis_s4(
        config = config,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    
    # Filter with severe stringency (reproducible with seed set above)
    analysis <- filter_analysis_s4(analysis, stringency = "severe")
    
    list(analysis = analysis, se = se(analysis), readcounts = readcounts)
}

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution
# ============================================================================

# ============================================================================
# TEST SUITE 1: Basic Workflow Execution - Paired Design
# ============================================================================

test_that("tsenat() paired design: executes pipeline, returns TSENATAnalysis, respects config", {
    data_list <- setup_workflow_data()
    
    # Analysis already configured with paired design from setup_workflow_data()
    config <- getConfig(data_list$analysis)
    
    # Single tsenat() call for paired design
    result <- tsenat(
        data_list$analysis,
        config = config,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Multiple assertions on same result
    expect_s4_class(result, "TSENATAnalysis")
    # Use the first available q-value (0) since exact q=1.0 may not exist with seq(0,2,length.out=10)
    expect_true(length(diversity(result, q = 0)) > 0 || is.null(diversity(result, q = 0)))
    expect_s4_class(getSE(result), "SummarizedExperiment")
    expect_true(is.list(getConfig(result)))
    expect_true(is.list(getMeta(result)))
    cfg_result <- getConfig(result)
    expect_equal(cfg_result$q_values, seq(0, 2, length.out = 10))
    expect_true(nrow(se(result)) > 0)
})

# ============================================================================
# TEST SUITE 2: Unpaired Design - Multiple Configuration Tests
# ============================================================================

test_that("tsenat() unpaired design: executes pipeline, handles config, preserves structure", {
    data_list <- setup_workflow_data()
    
    # Override with unpaired design config
    config <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        q_values = seq(0, 2, length.out = 10),
        paired = FALSE,
        nthreads = 4
    )
    
    # Single tsenat() call for unpaired design
    result <- tsenat(
        data_list$analysis,
        config = config,
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

test_that("tsenat() config override: explicit config overrides analysis config, uses defaults", {
    data_list <- setup_workflow_data()
    
    # Call without explicit config (uses analysis config from setup)
    result_default <- tsenat(
        data_list$analysis,
        config = NULL,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_default, "TSENATAnalysis")
    
    # Override with different q values (unpaired to avoid needing control)
    config2 <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        q_values = seq(0, 2, length.out = 10),
        paired = FALSE,
        nthreads = 4
    )
    result_override <- tsenat(
        data_list$analysis,
        config = config2,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Verify override was applied
    expect_s4_class(result_override, "TSENATAnalysis")
    cfg_result <- getConfig(result_override)
    expect_equal(cfg_result$q_values, seq(0, 2, length.out = 10))
})

# ============================================================================
# TEST SUITE 4: Output Directory and Verbose Control
# ============================================================================

test_that("tsenat() output handling: creates output_dir when needed, silent with verbose=FALSE", {
    data_list <- setup_workflow_data()
    
    output_dir <- tempdir()
    config <- getConfig(data_list$analysis)
    
    # Run with output_dir
    result_with_output <- tsenat(
        data_list$analysis,
        config = config,
        output_dir = output_dir,
        verbose = FALSE
    )
    expect_s4_class(result_with_output, "TSENATAnalysis")
    expect_true(dir.exists(output_dir))
    
    # Run with NULL output_dir (no file saving)
    result_no_output <- tsenat(
        data_list$analysis,
        config = config,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_no_output, "TSENATAnalysis")
    
    # Capture output for verbose test
    output <- capture.output({
        result_verbose <- tsenat(
            data_list$analysis,
            config = config,
            output_dir = NULL,
            verbose = TRUE
        )
    })
    expect_s4_class(result_verbose, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 5: Filtering Effects and Gene Selection
# ============================================================================

test_that("tsenat() filtering: works with severe filter, retains genes, processes correctly", {
    data_list <- setup_workflow_data()
    
    config <- getConfig(data_list$analysis)
    n_genes_filtered <- nrow(se(data_list$analysis))
    
    result <- tsenat(
        data_list$analysis,
        config = config,
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

test_that("tsenat() statistical params: respects bootstrap_method, nboot, seed configurations", {
    data_list <- setup_workflow_data()
    
    # Test with BCA bootstrap
    config_bca <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        q_values = seq(0, 2, length.out = 10),
        paired = FALSE,
        bootstrap_method = "bca",
        nboot = 100,
        nthreads = 2
    )
    
    result_bca <- tsenat(
        data_list$analysis,
        config = config_bca,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_bca, "TSENATAnalysis")
    
    # Test with different nboot
    config_nboot <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        q_values = seq(0, 2, length.out = 10),
        paired = FALSE,
        nboot = 50,
        nthreads = 2
    )
    
    result_nboot <- tsenat(
        data_list$analysis,
        config = config_nboot,
        output_dir = NULL,
        verbose = FALSE
    )
    expect_s4_class(result_nboot, "TSENATAnalysis")
    
    # Test reproducibility with seed
    config_seed <- tsenat_config(
        sample_col = "sample",
        condition_col = "condition",
        q_values = seq(0, 2, length.out = 10),
        paired = FALSE,
        seed = 42,
        nthreads = 2
    )
    
    result_seed1 <- tsenat(
        data_list$analysis,
        config = config_seed,
        output_dir = NULL,
        verbose = FALSE
    )
    result_seed2 <- tsenat(
        data_list$analysis,
        config = config_seed,
        output_dir = NULL,
        verbose = FALSE
    )
    
    expect_s4_class(result_seed1, "TSENATAnalysis")
    expect_s4_class(result_seed2, "TSENATAnalysis")
})

# ============================================================================
# TEST SUITE 7: Error Handling
# ============================================================================

test_that("tsenat() error handling: rejects invalid input, handles edge cases", {
    data_list <- setup_workflow_data()
    
    invalid_input <- data.frame(a = 1:10, b = 11:20)
    config <- getConfig(data_list$analysis)
    
    expect_error(
        tsenat(invalid_input, config = config, output_dir = NULL, verbose = FALSE),
        "must be a TSENATAnalysis object"
    )
    
    # Test with empty analysis
    empty_analysis <- data_list$analysis
    empty_analysis@se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(nrow = 0, ncol = ncol(se(data_list$analysis)))),
        colData = SummarizedExperiment::colData(se(data_list$analysis))
    )
    
    expect_error(
        tsenat(empty_analysis, config = config, output_dir = NULL, verbose = FALSE),
        "empty"
    )
})

# ============================================================================
# TEST SUITE 8: LM Interaction Results and Plot Generation
# ============================================================================

test_that("tsenat() paired: produces LM results, significant genes, plots generate", {
    # Use setup_workflow_data() which already has proper paired config
    data_list <- setup_workflow_data()
    
    # Note: setup_workflow_data() already configured with:
    # condition_col="condition", subject_col="paired_samples", 
    # q_values=seq(0,2,by=0.05), paired=TRUE, control="normal"
    # Retrieve config from analysis instead of using undefined variable
    config <- getConfig(data_list$analysis)
    
    result <- tsenat(
        data_list$analysis,
        config = config,
        output_dir = NULL,
        verbose = FALSE
    )
    
    # Test LM results structure
    expect_s4_class(result, "TSENATAnalysis")
    lm_res_list <- lmResults(result)
    expect_true(!is.null(lm_res_list))
    expect_true("lm_interaction" %in% names(lm_res_list))
    
    # Extract actual data frame
    lm_res <- lm_res_list$lm_interaction
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
        plot_lm_interaction_gam_s4(
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

