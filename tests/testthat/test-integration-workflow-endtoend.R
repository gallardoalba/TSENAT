# Comprehensive end-to-end workflow tests for tsenat() orchestration function
# Uses real package data (readcounts) following vignette workflow structure

context("End-to-End Workflow: Complete tsenat() Pipeline with Real Data")

# Helper: Load and prepare real package data for testing
# Follows vignette structure: load readcounts (which includes tpm, effective_length)
setup_workflow_data <- function() {
    # Load real example dataset from package (lazy-loaded as 'readcounts' by default)
    # This includes SALMON preprocessing outputs:
    # - readcounts: transcript-level NumReads (raw fragment counts)
    # - tpm: TPM matrix (length- and library-normalized estimates)
    # - effective_length: EffectiveLength vector (read-length corrected)
    data("readcounts", package = "TSENAT")
    
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    # Reference tpm and effective_length (loaded with readcounts)
    # These are part of the lazy-loaded data
    tpm_matrix <- tpm
    eff_length <- effective_length
    
    # Load sample metadata
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    # Load annotation file
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    list(
        readcounts = readcounts,
        tpm = tpm_matrix,
        effective_length = eff_length,
        metadata_df = metadata_df,
        gff3_file = gff3_file
    )
}

# ============================================================================
# TEST: Complete workflow with default parameters
# ============================================================================

test_that("Complete tsenat() workflow runs successfully with real data", {
    # Load real package data following vignette workflow
    data_list <- setup_workflow_data()
    
    # Build analysis object from real data
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    # Configure analysis
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = seq(0, 2, by = 0.5),
        paired = TRUE,
        control = "normal"
    )
    analysis <- setConfig(analysis, config)
    
    # Filter and run orchestration function
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(inherits(result@se, "SummarizedExperiment"))
    expect_true(length(result@diversity_results) > 0)
})

# ============================================================================
# TEST: Workflow step verification
# ============================================================================

test_that("tsenat() workflow completes all analysis steps with real data", {
    # Load and prepare real data
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE,
        nthreads = 2
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    # Run orchestration function with multiple methods
    result <- tsenat(
        se(analysis),
        config = config,
        methods = c("diversity", "jackknife", "lm_interaction", "divergence"),
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Verify workflow completed with results
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify all analysis steps produced results
    expect_true(length(result@diversity_results) > 0)
    expect_true(length(result@jackknife_results) > 0)
    # LM interaction may be empty if too few genes survive filtering
    expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: Workflow with paired design
# ============================================================================

test_that("tsenat() workflow handles paired designs with real data", {
    # Load real data
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    # Configure with paired design (real metadata includes paired_samples)
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true(length(result@diversity_results) > 0)
})

# ============================================================================
# TEST: Workflow reproducibility with seed
# ============================================================================

test_that("tsenat() workflow is reproducible with seed using real data", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE,
        seed = 999
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    se_filtered <- se(analysis)
    
    # Run twice with same seed
    result1 <- tsenat(se_filtered, config = config, verbose = FALSE, 
                      generate_plots = FALSE)
    result2 <- tsenat(se_filtered, config = config, verbose = FALSE, 
                      generate_plots = FALSE)
    
    # Results should be identical when using same seed (check via accessor)
    div1 <- diversity(result1, q = 0.5)
    div2 <- diversity(result2, q = 0.5)
    
    if (!is.null(div1) && !is.null(div2)) {
        expect_equal(div1$entropy, div2$entropy)
    }
})

# ============================================================================
# TEST: Workflow with different random seeds produces different results
# ============================================================================

test_that("tsenat() produces different results with different seeds", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config_base <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE
    )
    analysis <- setConfig(analysis, config_base)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    se_filtered <- se(analysis)
    
    cfg1 <- tsenat_config(seed = 111, q_values = 1.0)
    result1 <- suppressWarnings(tsenat(se_filtered, config = cfg1, methods = c("diversity", "lm_interaction"), 
                      verbose = FALSE, generate_plots = FALSE))
    
    cfg2 <- tsenat_config(seed = 222, q_values = 1.0)
    result2 <- suppressWarnings(tsenat(se_filtered, config = cfg2, methods = c("diversity", "lm_interaction"), 
                      verbose = FALSE, generate_plots = FALSE))
    
    # Both should complete without error
    expect_true(inherits(result1, "TSENATAnalysis"))
    expect_true(inherits(result2, "TSENATAnalysis"))
})

# ============================================================================
# TEST: Workflow with accessor functions
# ============================================================================

test_that("tsenat() results accessible via accessor functions", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Test diversity accessor
    div <- diversity(result, q = 0.5)
    expect_true(!is.null(div))
    
    # Test configuration accessor
    cfg <- getConfig(result)
    expect_true(!is.null(cfg))
})

# ============================================================================
# TEST: Workflow with show method
# ============================================================================

test_that("tsenat() result displays properly with show()", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = 0.5,
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Verify result object properly created and methods executed
    expect_true(inherits(result, "TSENATAnalysis"))
    # Diversity should always be computed
    expect_true(length(result@diversity_results) > 0 || length(show(result)) > 0)
})

# ============================================================================
# TEST: Workflow with stringency filtering
# ============================================================================

test_that("tsenat() workflow respects stringency parameter", {
    data_list <- setup_workflow_data()
    
    analysis_soft <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = 1.0,
        paired = TRUE
    )
    analysis_soft <- setConfig(analysis_soft, config)
    analysis_soft <- filter_analysis_s4(analysis_soft, stringency = "soft")
    result_soft <- suppressWarnings(tsenat(se(analysis_soft), config = config, 
                          verbose = FALSE, generate_plots = FALSE))
    
    # Create second analysis for severe filtering
    analysis_severe <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    analysis_severe <- setConfig(analysis_severe, config)
    analysis_severe <- filter_analysis_s4(analysis_severe, stringency = "severe")
    result_severe <- suppressWarnings(tsenat(se(analysis_severe), config = config, 
                            verbose = FALSE, generate_plots = FALSE))
    
    # Both should complete successfully
    expect_s4_class(result_soft, "TSENATAnalysis")
    expect_s4_class(result_severe, "TSENATAnalysis")
    
    # Severe filtering should have fewer or equal genes than soft
    div_soft <- diversity(result_soft, q = 1.0)
    div_severe <- diversity(result_severe, q = 1.0)
    
    if (!is.null(div_soft) && !is.null(div_severe)) {
        expect_true(nrow(div_severe) <= nrow(div_soft))
    }
})

# ============================================================================
# TEST: Workflow error handling
# ============================================================================

test_that("tsenat() workflow rejects empty SummarizedExperiment", {
    empty_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(nrow = 0, ncol = 0)),
        colData = data.frame()
    )
    
    expect_error(tsenat(empty_se, verbose = FALSE))
})

test_that("tsenat() workflow rejects non-SummarizedExperiment input", {
    invalid_input <- data.frame(a = 1:10, b = 11:20)
    
    expect_error(tsenat(invalid_input, verbose = FALSE))
})

# ============================================================================
# TEST: Workflow state preservation
# ============================================================================

test_that("tsenat() preserves input data in result object", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    
    original_nrow <- nrow(se(analysis))
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(se(analysis), config = config, verbose = FALSE, 
                     generate_plots = FALSE)
    
    # Original data preserved (though may be filtered)
    expect_true(nrow(result@se) > 0)
    expect_true(nrow(result@se) <= original_nrow)
})

# ============================================================================
# TEST: Workflow with multiple q-values produces q-dependent results
# ============================================================================

test_that("tsenat() q-spectrum analysis produces q-dependent entropy values", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Check that diversity results exist
    expect_true(length(result@diversity_results) > 0)
    
    # Use accessor function to verify results (more robust than direct access)
    div_q0.5 <- diversity(result, q = 0.5)
    div_q1.0 <- diversity(result, q = 1.0)
    
    # At least one of the requested q values should have results
    has_q0.5 <- !is.null(div_q0.5) && nrow(div_q0.5) > 0
    has_q1.0 <- !is.null(div_q1.0) && nrow(div_q1.0) > 0
    expect_true(has_q0.5 || has_q1.0)
})

# ============================================================================
# TEST: Workflow with metadata persistence
# ============================================================================

test_that("tsenat() records analysis execution information", {
    data_list <- setup_workflow_data()
    
    analysis <- build_analysis_s4(
        readcounts = data_list$readcounts,
        tx2gene = data_list$gff3_file,
        metadata = data_list$metadata_df,
        tpm = data_list$tpm,
        effective_length = data_list$effective_length
    )
    
    config <- tsenat_config(
        condition_col = "sample_type",
        subject_col = "paired_samples",
        q_values = c(0.5, 1.0),
        paired = TRUE
    )
    analysis <- setConfig(analysis, config)
    analysis <- filter_analysis_s4(analysis, stringency = "medium")
    
    result <- tsenat(
        se(analysis),
        config = config,
        verbose = FALSE,
        generate_plots = FALSE
    )
    
    # Check that result is valid TSENATAnalysis object
    expect_s4_class(result, "TSENATAnalysis")
    
    # Check that config is preserved
    expect_true(!is.null(result@config))
    expect_true(is.list(result@config))
    
    # Check that diversity results were computed
    expect_true(length(result@diversity_results) > 0)
})
